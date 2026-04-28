%SCANMFCPIPPERFMAP Scan CP/IP performance over SNR and window length.
% Dev scan script. Edit the configuration section directly before running.
% It compares CP-K / CP-U from the current simple-flow bundle against
% IP-K / IP-U reruns from the same static seed, and stores only lightweight
% scanData tables and plot data.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfCpIpPerfMap";
saveSnapshot = true;
optVerbose = false;
baseSeed = 253;
numRepeat = 3;
snrDbList = [0; 5; 10; 15];
frameCountList = [8; 10; 20];
toothResidualTolHz = 50;

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);
snrDbList = reshape(double(snrDbList), [], 1);
frameCountList = reshape(double(frameCountList), [], 1);

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.optVerbose = logical(optVerbose);
scanConfig.baseSeed = baseSeed;
scanConfig.numRepeat = numRepeat;
scanConfig.seedList = seedList;
scanConfig.snrDbList = snrDbList;
scanConfig.frameCountList = frameCountList;
scanConfig.toothResidualTolHz = toothResidualTolHz;

runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
localPrintScanHeader(scanName, scanConfig);

try
  %% Run scan batch
  taskList = localBuildTaskList(frameCountList, snrDbList, seedList);
  taskOutCell = cell(numel(taskList), 1);
  progressTracker = localCreateScanProgressTracker(numel(taskList));
  progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>

  useParfor = localCanUseParfor(numel(taskList));
  if useParfor && progressTracker.enabled && ~isempty(progressTracker.queue)
    progressQueue = progressTracker.queue;
    parfor iTask = 1:numel(taskList)
      taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanConfig);
      send(progressQueue, 1);
    end
  elseif useParfor
    parfor iTask = 1:numel(taskList)
      taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanConfig);
    end
  else
    for iTask = 1:numel(taskList)
      taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanConfig);
      localAdvanceScanProgress(progressTracker, false);
    end
  end

  localCloseScanProgressTracker(progressTracker);
  clear progressCleanup;

  [perfTable, repeatOutCell] = localCollectTaskOutput(taskOutCell);
  aggregateTable = localBuildAggregateTable(perfTable, toothResidualTolHz);

  %% Data storage
  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.perfTable = perfTable;
  scanData.aggregateTable = aggregateTable;
  scanData.repeatOutCell = repeatOutCell;
  scanData.plotData = localBuildPlotData(aggregateTable);

  if saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end
catch ME
  if exist('progressTracker', 'var')
    localCloseScanProgressTracker(progressTracker);
  end
  if exist('progressCleanup', 'var')
    clear progressCleanup;
  end
  fprintf('Scan failed while building CP/IP performance map data.\n');
  rethrow(ME);
end

%% Summary output and plotting
if ~exist('scanData', 'var') || ~isstruct(scanData)
  error('scanMfCpIpPerfMap:MissingScanData', ...
    'Scan data is missing. Run the scan batch sections or load a snapshot containing scanData.');
end

fprintf('\n========== CP/IP performance-map aggregate ==========\n');
fprintf('Frames                          : %s\n', localFormatRow(scanData.config.frameCountList));
fprintf('SNR list (dB)                   : %s\n', localFormatRow(scanData.config.snrDbList));
fprintf('Repeats per config              : %d\n', scanData.config.numRepeat);
fprintf('Truth-tooth rule                : toothIdx == 0 and |residual| <= %.2f Hz\n', ...
  scanData.config.toothResidualTolHz);
disp(scanData.aggregateTable);
scanData.plotData = localPlotScan(scanData);

%% Local helpers

function taskList = localBuildTaskList(frameCountList, snrDbList, seedList)
%LOCALBUILDTASKLIST Build independent scan tasks for the outer batch loop.

taskTemplate = struct('numFrame', NaN, 'snrDb', NaN, 'taskSeed', NaN);
taskList = repmat(taskTemplate, numel(frameCountList) * numel(snrDbList) * numel(seedList), 1);
idxTask = 0;
for iP = 1:numel(frameCountList)
  for iSnr = 1:numel(snrDbList)
    for iSeed = 1:numel(seedList)
      idxTask = idxTask + 1;
      taskList(idxTask).numFrame = frameCountList(iP);
      taskList(idxTask).snrDb = snrDbList(iSnr);
      taskList(idxTask).taskSeed = seedList(iSeed);
    end
  end
end
end

function taskOut = localRunOneTask(task, scanConfig)
%LOCALRUNONETASK Run one frame-count, SNR, and seed configuration.

numFrame = task.numFrame;
offsetIdx = localCenteredOffsets(numFrame);
[subsetOffsetCell, subsetLabelList] = localBuildFrameCountSubsetBank(numFrame);
masterOffsetIdx = localBuildMasterOffsets(offsetIdx, subsetOffsetCell);
contextOpt = struct();
contextOpt.periodicOffsetIdx = offsetIdx;
contextOpt.masterOffsetIdx = masterOffsetIdx;
contextOpt.subsetOffsetCell = subsetOffsetCell;
contextOpt.subsetLabelList = subsetLabelList;
context = buildDynamicDualSatEciContext(contextOpt);
flowOpt = localBuildFlowOpt(subsetLabelList);
repeatOut = localRunOneRepeat(context, flowOpt, task.snrDb, task.taskSeed, scanConfig.optVerbose);

taskOut = struct();
taskOut.caseRowList = repeatOut.caseRowList;
taskOut.repeatSlim = localStripRepeatOut(repeatOut);
end

function flowOpt = localBuildFlowOpt(subsetLabelList)
%LOCALBUILDFLOWOPT Build the simple dynamic flow options used by this scan.

subsetLabelList = reshape(string(subsetLabelList), 1, []);
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'subsetDefaultBankLabelList', subsetLabelList, ...
  'subsetMarginFallbackBankLabelList', strings(1, 0), ...
  'subsetEnableMarginFallback', false, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti", ...
  'periodicPolishEnableWhenMulti', true, ...
  'periodicPolishDoaHalfWidthDeg', [0.002; 0.002]));
flowOpt = applyDynamicTransitionFlowDefaults(flowOpt);
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function [subsetOffsetCell, subsetLabelList] = localBuildFrameCountSubsetBank(numFrame)
%LOCALBUILDFRAMECOUNTSUBSETBANK Build valid subset schedules for one frame count.

numFrame = round(numFrame);
if numFrame < 2
  error('scanMfCpIpPerfMap:InvalidFrameCount', ...
    'Frame count must be at least two for multi-frame CP/IP scanning.');
end
subsetSize = min(numFrame, 10);
baseOffsetCell = { ...
  [-8, -7, -5, -4, -2, -1, 0, 4, 7, 9], ...
  [-7, -4, -1, 0, 3, 5, 7, 8, 9, 10]};
subsetOffsetCell = cell(1, numel(baseOffsetCell));
for iSubset = 1:numel(baseOffsetCell)
  subsetOffsetCell{iSubset} = localResizeSubsetOffsets(baseOffsetCell{iSubset}, subsetSize);
end
subsetLabelList = ["curated1"; "curated2"];
end

function offsetIdx = localResizeSubsetOffsets(baseOffsetIdx, subsetSize)
%LOCALRESIZESUBSETOFFSETS Select a fixed-size schedule while keeping the reference frame.

baseOffsetIdx = unique(reshape(double(baseOffsetIdx), 1, []), 'stable');
if subsetSize >= numel(baseOffsetIdx)
  offsetIdx = sort(baseOffsetIdx);
  return;
end
if ~any(baseOffsetIdx == 0)
  error('scanMfCpIpPerfMap:MissingReferenceFrame', ...
    'The scan subset template must include one zero-offset frame.');
end
nonzeroOffset = baseOffsetIdx(baseOffsetIdx ~= 0);
numNonzero = subsetSize - 1;
if numNonzero < 1 || numNonzero > numel(nonzeroOffset)
  error('scanMfCpIpPerfMap:InvalidSubsetSize', ...
    'Invalid subset size requested for CP/IP scan subset construction.');
end
pickIdx = unique(round(linspace(1, numel(nonzeroOffset), numNonzero)), 'stable');
if numel(pickIdx) < numNonzero
  for iCand = 1:numel(nonzeroOffset)
    if numel(pickIdx) >= numNonzero
      break;
    end
    if ~any(pickIdx == iCand)
      pickIdx(end + 1) = iCand; %#ok<AGROW>
    end
  end
end
pickIdx = sort(pickIdx(1:numNonzero));
offsetIdx = sort([0, nonzeroOffset(pickIdx)]);
end

function masterOffsetIdx = localBuildMasterOffsets(periodicOffsetIdx, subsetOffsetCell)
%LOCALBUILDMASTEROFFSETS Build one contiguous master window covering all scan fixtures.

allOffset = reshape(periodicOffsetIdx, 1, []);
for iSubset = 1:numel(subsetOffsetCell)
  allOffset = [allOffset, reshape(subsetOffsetCell{iSubset}, 1, [])]; %#ok<AGROW>
end
minOffset = min(allOffset);
maxOffset = max(allOffset);
masterOffsetIdx = minOffset:maxOffset;
if ~any(masterOffsetIdx == 0)
  masterOffsetIdx = sort([masterOffsetIdx, 0]);
end
end

function repeatOut = localRunOneRepeat(context, flowOpt, snrDb, taskSeed, optVerbose)
%LOCALRUNONEREPEAT Run one repeat and collect CP/IP dynamic cases.

repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);
bundle = buildDoaDopplerDynamicTransitionBundle( ...
  repeatData.periodicFixture, repeatData.subsetFixtureCell, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  optVerbose, flowOpt);
caseIpKnown = localRunIpKnownCase(repeatData.periodicFixture, bundle.bestStaticMsCase, ...
  context, flowOpt, optVerbose);
caseIpUnknown = localRunIpUnknownCase(repeatData.periodicFixture, bundle.bestStaticMsCase, ...
  bundle.caseBundle.staticMsOpt, caseIpKnown, context, flowOpt, optVerbose);
caseRowList = repmat(localEmptyCaseRow(), 4, 1);
caseRowList(1) = localBuildCaseRow("CP-K", "continuous", "known", bundle.caseDynMsKnown, ...
  repeatData.truth, snrDb, taskSeed, context);
caseRowList(2) = localBuildCaseRow("CP-U", "continuous", "unknown", bundle.caseDynMsUnknownFinal, ...
  repeatData.truth, snrDb, taskSeed, context);
caseRowList(3) = localBuildCaseRow("IP-K", "independent", "known", caseIpKnown, ...
  repeatData.truth, snrDb, taskSeed, context);
caseRowList(4) = localBuildCaseRow("IP-U", "independent", "unknown", caseIpUnknown, ...
  repeatData.truth, snrDb, taskSeed, context);
repeatOut = struct();
repeatOut.taskSeed = taskSeed;
repeatOut.snrDb = snrDb;
repeatOut.truth = repeatData.truth;
repeatOut.bundle = bundle;
repeatOut.caseIpKnown = caseIpKnown;
repeatOut.caseIpUnknown = caseIpUnknown;
repeatOut.caseRowList = caseRowList;
end

function caseUse = localRunIpKnownCase(periodicFixture, bestStaticMsCase, context, flowOpt, optVerbose)
%LOCALRUNIPKNOWNCASE Rerun the known-rate case with independent phase tying.

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = 'independent';
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth(:);
dynOpt.enableFdAliasUnwrap = true;
initParam = buildDynamicInitParamFromCase(bestStaticMsCase, true, periodicFixture.truth.fdRateFit);
caseUse = runDynamicDoaDopplerCase("MS-MF-IP-K", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, periodicFixture.fdRange, periodicFixture.fdRateRange, ...
  optVerbose, dynOpt, true, struct(), initParam);
end

function caseUse = localRunIpUnknownCase(periodicFixture, bestStaticMsCase, staticMsOpt, caseIpKnown, context, flowOpt, optVerbose)
%LOCALRUNIPUNKNOWNCASE Rerun the unknown-rate case with independent phase tying.

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = 'independent';
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynOpt.initDoaHalfWidth = flowOpt.msUnknownDoaHalfWidth(:);
dynOpt.enableFdAliasUnwrap = true;
initParamKnown = buildDynamicInitParamFromCase(caseIpKnown, false, caseIpKnown.estResult.fdRateEst);
initParamStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, periodicFixture.truth.fdRateFit);
initCandidate = buildUnknownInitCandidateSet(caseIpKnown, bestStaticMsCase, initParamKnown, initParamStatic, ...
  dynOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseUse = runDynamicDoaDopplerCase("MS-MF-IP-U", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, periodicFixture.fdRange, periodicFixture.fdRateRange, ...
  optVerbose, dynOpt, false, struct(), initCandidate);
end

function row = localBuildCaseRow(displayName, phaseMode, fdRateMode, caseUse, truth, snrDb, taskSeed, context)
%LOCALBUILDCASEMETRICROW Convert one estimator case into a compact scan row.

row = localEmptyCaseRow();
row.displayName = string(displayName);
row.phaseMode = string(phaseMode);
row.fdRateMode = string(fdRateMode);
row.snrDb = snrDb;
row.taskSeed = taskSeed;
row.numFrame = numel(context.periodicOffsetIdx);
row.frameIntvlSec = context.frameIntvlSec;
row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
[angleErrDeg, fdErrHz, fdRateErrHzPerSec, isResolved] = extractDynamicCaseMetric(caseUse, truth);
row.angleErrDeg = angleErrDeg;
row.fdRefAbsErrHz = fdErrHz;
row.fdRateAbsErrHzPerSec = fdRateErrHzPerSec;
row.isResolved = isResolved;
if isfield(caseUse, 'estResult')
  row.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
  row.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
  row.fdRefSignedErrHz = localGetFieldOrDefault(caseUse.estResult, 'fdRefEst', NaN) - truth.fdRefFit;
  row.fdRateSignedErrHzPerSec = localGetFieldOrDefault(caseUse.estResult, 'fdRateEst', NaN) - truth.fdRateFit;
  row.fdRefAbsErrHz = abs(row.fdRefSignedErrHz);
  row.fdRateAbsErrHzPerSec = abs(row.fdRateSignedErrHzPerSec);
end
row.toothIdx = round(row.fdRefSignedErrHz / (1 / context.frameIntvlSec));
row.toothResidualHz = row.fdRefSignedErrHz - row.toothIdx * (1 / context.frameIntvlSec);
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for per-repeat case metrics.

row = struct('displayName', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'taskSeed', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, ...
  'windowSec', NaN, 'angleErrDeg', NaN, 'fdRefAbsErrHz', NaN, ...
  'fdRateAbsErrHzPerSec', NaN, 'fdRefSignedErrHz', NaN, ...
  'fdRateSignedErrHzPerSec', NaN, 'toothIdx', NaN, 'toothResidualHz', NaN, ...
  'isResolved', false, 'solveVariant', "", 'runTimeMs', NaN);
end

function [perfTable, repeatOutCell] = localCollectTaskOutput(taskOutCell)
%LOCALCOLLECTTASKOUTPUT Merge per-task rows and slim repeat diagnostics.

rowList = repmat(localEmptyCaseRow(), 0, 1);
repeatOutCell = cell(numel(taskOutCell), 1);
for iTask = 1:numel(taskOutCell)
  taskOut = taskOutCell{iTask};
  if isempty(taskOut)
    continue;
  end
  rowList = [rowList; taskOut.caseRowList(:)]; %#ok<AGROW>
  repeatOutCell{iTask} = taskOut.repeatSlim;
end
perfTable = struct2table(rowList(:));
end

function aggregateTable = localBuildAggregateTable(perfTable, toothResidualTolHz)
%LOCALBUILDAGGREGATETABLE Aggregate RMSE and tooth-hit metrics by case and grid point.

groupKey = unique(perfTable(:, {'displayName', 'phaseMode', 'fdRateMode', 'snrDb', 'numFrame', 'frameIntvlSec'}), 'rows', 'stable');
rowList = repmat(localEmptyAggregateRow(), height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iRow) & ...
    perfTable.snrDb == groupKey.snrDb(iRow) & ...
    perfTable.numFrame == groupKey.numFrame(iRow) & ...
    abs(perfTable.frameIntvlSec - groupKey.frameIntvlSec(iRow)) < eps;
  row = localEmptyAggregateRow();
  row.displayName = groupKey.displayName(iRow);
  row.phaseMode = groupKey.phaseMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.numFrame = groupKey.numFrame(iRow);
  row.frameIntvlSec = groupKey.frameIntvlSec(iRow);
  row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
  row.numRepeat = nnz(mask);
  row.resolvedRate = mean(double(perfTable.isResolved(mask)), 'omitnan');
  row.truthToothHitRate = mean(double(abs(perfTable.toothIdx(mask)) == 0 & ...
    abs(perfTable.toothResidualHz(mask)) <= toothResidualTolHz), 'omitnan');
  row.angleRmseDeg = sqrt(mean(perfTable.angleErrDeg(mask).^2, 'omitnan'));
  row.angleP95Deg = prctile(perfTable.angleErrDeg(mask), 95);
  row.fdRefRmseHz = sqrt(mean(perfTable.fdRefAbsErrHz(mask).^2, 'omitnan'));
  row.fdRateRmseHzPerSec = sqrt(mean(perfTable.fdRateAbsErrHzPerSec(mask).^2, 'omitnan'));
  row.meanRunTimeMs = mean(perfTable.runTimeMs(mask), 'omitnan');
  rowList(iRow) = row;
end
aggregateTable = struct2table(rowList);
aggregateTable = localAddCpIpRatios(aggregateTable);
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return a typed empty row for aggregate metrics.

row = struct('displayName', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, 'windowSec', NaN, ...
  'numRepeat', NaN, 'resolvedRate', NaN, 'truthToothHitRate', NaN, ...
  'angleRmseDeg', NaN, 'angleP95Deg', NaN, 'fdRefRmseHz', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'meanRunTimeMs', NaN, 'angleRatioToCp', NaN, ...
  'fdRatioToCp', NaN);
end

function aggregateTable = localAddCpIpRatios(aggregateTable)
%LOCALADDCPIPRATIOS Add IP-to-CP ratios for matching SNR and frame-count rows.

aggregateTable.angleRatioToCp = nan(height(aggregateTable), 1);
aggregateTable.fdRatioToCp = nan(height(aggregateTable), 1);
for iRow = 1:height(aggregateTable)
  if aggregateTable.phaseMode(iRow) ~= "independent"
    continue;
  end
  cpName = "CP-K";
  if aggregateTable.fdRateMode(iRow) == "unknown"
    cpName = "CP-U";
  end
  maskCp = aggregateTable.displayName == cpName & aggregateTable.snrDb == aggregateTable.snrDb(iRow) & ...
    aggregateTable.numFrame == aggregateTable.numFrame(iRow) & ...
    abs(aggregateTable.frameIntvlSec - aggregateTable.frameIntvlSec(iRow)) < eps;
  idxCp = find(maskCp, 1, 'first');
  if ~isempty(idxCp)
    aggregateTable.angleRatioToCp(iRow) = aggregateTable.angleRmseDeg(iRow) / aggregateTable.angleRmseDeg(idxCp);
    aggregateTable.fdRatioToCp(iRow) = aggregateTable.fdRefRmseHz(iRow) / aggregateTable.fdRefRmseHz(idxCp);
  end
end
end

function plotData = localBuildPlotData(aggregateTable)
%LOCALBUILDPLOTDATA Build lightweight tables used to redraw scan figures.

plotData = struct();
plotData.figureFiles = strings(0, 1);
plotData.angleVsSnrTable = localBuildAngleVsSnrPlotTable(aggregateTable);
plotData.angleRatioTable = aggregateTable(aggregateTable.phaseMode == "independent", :);
plotData.toothHitVsFrameTable = localBuildToothHitVsFramePlotTable(aggregateTable);
end

function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw the default CP/IP performance-map figures.

plotData = localBuildPlotData(scanData.aggregateTable);
angleTable = plotData.angleVsSnrTable;
caseList = unique(angleTable.displayName, 'stable');
fig = figure('Name', 'CP/IP angle RMSE vs SNR'); hold on;
for iCase = 1:numel(caseList)
  mask = angleTable.displayName == caseList(iCase);
  plot(angleTable.snrDb(mask), angleTable.angleRmseDeg(mask), '-o');
end
grid on; xlabel('SNR (dB)'); ylabel('angle RMSE (deg)'); legend(cellstr(caseList), 'Location', 'best');
drawnow;

ratioTable = plotData.angleRatioTable;
figure('Name', 'IP/CP angle ratio');
scatter(ratioTable.windowSec * 1e3, ratioTable.angleRatioToCp, 60, ratioTable.snrDb, 'filled');
grid on; xlabel('window span (ms)'); ylabel('IP / CP angle RMSE');
drawnow;

toothTable = plotData.toothHitVsFrameTable;
caseList = unique(toothTable.displayName, 'stable');
figure('Name', 'CP/IP truth-tooth hit rate'); hold on;
for iCase = 1:numel(caseList)
  mask = toothTable.displayName == caseList(iCase);
  plot(toothTable.numFrame(mask), toothTable.truthToothHitRate(mask), '-o');
end
grid on; xlabel('number of frames'); ylabel('truth-tooth hit rate'); legend(cellstr(caseList), 'Location', 'best');
drawnow;
end

function plotTable = localBuildAngleVsSnrPlotTable(aggregateTable)
%LOCALBUILDANGLEVSSNRPLOTTABLE Select rows for the angle-vs-SNR figure.

if isempty(aggregateTable)
  plotTable = aggregateTable;
  return;
end
maxFrame = max(aggregateTable.numFrame);
plotTable = aggregateTable(aggregateTable.numFrame == maxFrame, :);
end

function plotTable = localBuildToothHitVsFramePlotTable(aggregateTable)
%LOCALBUILDTOOTHHITVSFRAMEPLOTTABLE Select rows for the tooth-hit figure.

if isempty(aggregateTable)
  plotTable = aggregateTable;
  return;
end
maxSnr = max(aggregateTable.snrDb);
plotTable = aggregateTable(aggregateTable.snrDb == maxSnr, :);
end

function repeatSlim = localStripRepeatOut(repeatOut)
%LOCALSTRIPREPEATOUT Keep only lightweight repeat diagnostics for scanData.

repeatSlim = struct();
repeatSlim.taskSeed = repeatOut.taskSeed;
repeatSlim.snrDb = repeatOut.snrDb;
repeatSlim.truth = repeatOut.truth;
repeatSlim.caseRowList = repeatOut.caseRowList;
repeatSlim.selectedSubsetLabel = localGetFieldOrDefault(repeatOut.bundle, 'selectedSubsetLabel', "");
repeatSlim.selectedFinalTag = localGetFieldOrDefault(repeatOut.bundle, 'selectedFinalTag', "");
repeatSlim.subsetCandidateTable = localGetFieldOrDefault(repeatOut.bundle, 'subsetCandidateTable', table());
repeatSlim.finalSelectTable = localGetFieldOrDefault(repeatOut.bundle, 'finalSelectTable', table());
end

function txt = localFormatRow(x)
%LOCALFORMATROW Format a numeric vector as one compact comma-separated string.

txt = strjoin(compose('%.6g', x(:).'), ', ');
end

function localPrintScanHeader(scanName, config)
%LOCALPRINTSCANHEADER Print compact scan configuration before the batch run.

fprintf('Running %s ...\n', char(scanName));
fprintf('  frame count list                : %s\n', localFormatRow(config.frameCountList));
fprintf('  SNR list (dB)                   : %s\n', localFormatRow(config.snrDbList));
fprintf('  repeats per config             : %d\n', config.numRepeat);
fprintf('  base seed                       : %d\n', config.baseSeed);
fprintf('  save snapshot                   : %d\n', config.saveSnapshot);
fprintf('  estimator verbose               : %d\n', config.optVerbose);
end

function progressTracker = localCreateScanProgressTracker(totalCount)
%LOCALCREATESCANPROGRESSTRACKER Create a best-effort scan progress tracker.

progressTracker = struct('enabled', false, 'queue', []);
if totalCount <= 0
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  progressTracker.enabled = true;
  try
    progressTracker.queue = parallel.pool.DataQueue;
    afterEach(progressTracker.queue, @(~) progressbar('advance'));
  catch
    progressTracker.queue = [];
  end
catch
  progressTracker.enabled = false;
  progressTracker.queue = [];
end
end

function localAdvanceScanProgress(progressTracker, isParforWorker)
%LOCALADVANCESCANPROGRESS Advance progress from serial or parfor execution.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
if isParforWorker
  if isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
    send(progressTracker.queue, 1);
  end
else
  try
    progressbar('advance');
  catch
  end
end
end

function localCloseScanProgressTracker(progressTracker)
%LOCALCLOSESCANPROGRESSTRACKER Close a best-effort scan progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Return true when an outer scan loop can use parfor.

useParfor = false;
if numTask <= 1
  return;
end
try
  useParfor = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && ...
    license('test', 'Distrib_Computing_Toolbox');
catch
  useParfor = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct field with a fallback value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
