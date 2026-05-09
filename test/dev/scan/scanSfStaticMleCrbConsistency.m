%SCANSFSTATICMLECRBCONSISTENCY Scan single-frame static MLE consistency with CRB.
% Dev scan script. Edit the configuration section directly before running.
% It mirrors doaDopplerStatDualSatUraEciPerf, keeps the existing SF static
% estimator path, and reports full/resolved/CRB-local RMSE against the
% pilot/static DoA-Doppler CRB for the static paper-facing anchor.

clear; close all; clc;

%% Scan configuration

scanName = "scanSfStaticMleCrbConsistency";
saveSnapshot = true;
notifyTelegramEnable = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
optVerbose = false;
baseSeed = 253;
numRepeat = 2000;
snrDbList = (-20:3:15).';
trimEnable = true;
trimNormCap = 5;
trimMinKeepRate = 0.8;
topTailCountPerGroup = 5;
topTailPrintMaxRows = 40;
numUsr = 1;
numSym = 512;
sampleRate = 512e6;
symbolRate = 128e6;
carrierFreq = 11.7e9;
gridSize = [50, 50];
searchMarginDeg = 5;
fdRangeBase = [-2e5, 2e5];
staticMsHalfWidthDeg = [0.002; 0.002];
enableWeightSweep = false;
weightSweepAlpha = [0; 0.25; 0.5; 1];
methodNameList = [
  "SS-SF-DoA"
  "MS-SF-DoA"
  "SS-SF-Static"
  "MS-SF-Static"
];
% Keep the four canonical static-anchor cases above. Enable the optional
% weight sweep only for static diagnostic scans, not for the default paper
% anchor run.

if ~enableWeightSweep
  weightSweepAlpha = zeros(0, 1);
end
seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
staticMsHalfWidthDeg = reshape(double(staticMsHalfWidthDeg), [], 1);
methodNameList = reshape(string(methodNameList), [], 1);
methodNameList = methodNameList(strlength(methodNameList) > 0);
numRepeat = numel(seedList);

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.baseSeed = baseSeed;
scanConfig.numRepeat = numRepeat;
scanConfig.seedList = seedList;
scanConfig.snrDbList = snrDbList;
scanConfig.trimEnable = logical(trimEnable);
scanConfig.trimNormCap = trimNormCap;
scanConfig.trimMinKeepRate = trimMinKeepRate;
scanConfig.topTailCountPerGroup = topTailCountPerGroup;
scanConfig.topTailPrintMaxRows = topTailPrintMaxRows;
scanConfig.numUsr = numUsr;
scanConfig.numSym = numSym;
scanConfig.sampleRate = sampleRate;
scanConfig.symbolRate = symbolRate;
scanConfig.carrierFreq = carrierFreq;
scanConfig.gridSize = reshape(double(gridSize), 1, []);
scanConfig.searchMarginDeg = searchMarginDeg;
scanConfig.fdRangeBase = reshape(double(fdRangeBase), 1, []);
scanConfig.staticMsHalfWidthDeg = staticMsHalfWidthDeg;
scanConfig.enableWeightSweep = logical(enableWeightSweep);
scanConfig.weightSweepAlpha = reshape(double(weightSweepAlpha), [], 1);
scanConfig.methodNameList = methodNameList;
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.checkpointResume = logical(checkpointResume);
scanConfig.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);
scanConfig.optVerbose = logical(optVerbose);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
checkpointDir = "";
checkpointRunState = struct();
checkpointCleanupReport = struct([]);
scanData = struct();

try
  %% Build context and scan tasks

  taskList = localBuildTaskList(scanConfig.snrDbList, scanConfig.seedList);
  scanConfig.numTask = numel(taskList);
  scanConfig.runKey = runKey;
  useParfor = localCanUseParfor(numel(taskList));
  fprintf('Preparing static scan runtime cache...\n');
  scanRuntime = localBuildScanRuntime(scanConfig);
  if scanConfig.checkpointEnable
    checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, scanRuntime, taskList, useParfor);
    checkpointDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
  else
    checkpointOpt = struct();
  end
  printMfScanHeader(char(scanName), scanConfig, checkpointDir);

  %% Run scan batch

  if scanConfig.checkpointEnable
    pendingCount = localCountPendingCheckpointTasks(checkpointDir, numel(taskList));
    progressTracker = localCreateScanProgressTracker(pendingCount);
    progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>
    if progressTracker.enabled
      checkpointOpt.progressFcn = @(~) localAdvanceScanProgress(progressTracker, false);
    end
    checkpointRunState = runPerfTaskGridWithCheckpoint(taskList, scanRuntime, ...
      @localRunCheckpointTask, checkpointOpt);
    taskOutCell = checkpointRunState.resultCell;
    localCloseScanProgressTracker(progressTracker);
    clear progressCleanup;
  else
    taskOutCell = cell(numel(taskList), 1);
    progressTracker = localCreateScanProgressTracker(numel(taskList));
    progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>
    if useParfor && progressTracker.enabled && ~isempty(progressTracker.queue)
      progressQueue = progressTracker.queue;
      parfor iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
        send(progressQueue, 1);
      end
    elseif useParfor
      parfor iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
      end
    else
      for iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
        localAdvanceScanProgress(progressTracker, false);
      end
    end
    localCloseScanProgressTracker(progressTracker);
    clear progressCleanup;
  end

  repeatTable = localCollectTaskOutput(taskOutCell);
  crbTable = localBuildCrbTable(scanRuntime, scanConfig);
  perfTable = localAttachCrbMetrics(repeatTable, crbTable);
  aggregateTable = localBuildAggregateTable(perfTable, scanConfig);
  crbLocalSummaryTable = localBuildCrbLocalSummaryTable(aggregateTable);
  topTailTable = localBuildTopTailTable(perfTable, scanConfig.topTailCountPerGroup);
  topTailExportTable = localBuildTopTailExportTable(topTailTable);
  if scanConfig.checkpointEnable && localGetFieldOrDefault(checkpointRunState, 'isComplete', false) && ...
      scanConfig.checkpointCleanupOnSuccess
    checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(checkpointRunState, ...
      struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
  end
  checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunState, ...
    checkpointDir, checkpointCleanupReport);

  %% Data storage

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.truth = scanRuntime.truth;
  scanData.repeatTable = repeatTable;
  scanData.crbTable = crbTable;
  scanData.perfTable = perfTable;
  scanData.aggregateTable = aggregateTable;
  scanData.crbLocalSummaryTable = crbLocalSummaryTable;
  scanData.topTailTable = topTailTable;
  scanData.topTailExportTable = topTailExportTable;
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupReport = checkpointCleanupReport;
  scanData.plotData = localBuildPlotData(aggregateTable, crbTable);
  scanData.elapsedSec = toc(runTic);

  if ~scanConfig.checkpointEnable
    scanData = finalizeMfScanResult(scanData, "");
  end

  if scanConfig.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output and plotting

  printMfScanSection('SF static MLE/CRB consistency aggregate', scanData.aggregateTable);
  printMfScanSection('SF static CRB-local compact summary', scanData.crbLocalSummaryTable);
  if height(scanData.topTailTable) > 0
    printMfScanSection('Top CRB-normalized tail preview', ...
      localTablePreview(scanData.topTailTable, scanData.config.topTailPrintMaxRows));
    printMfScanSection('Top-tail compact export preview', ...
      localTablePreview(scanData.topTailExportTable, scanData.config.topTailPrintMaxRows));
  end
  if height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  fprintf('SNR list (dB)                   : %s\n', localFormatRow(scanData.config.snrDbList));
  fprintf('Active method list              : %s\n', localFormatStringRow(scanData.config.methodNameList));
  fprintf('Repeats per config              : %d\n', scanData.config.numRepeat);
  fprintf('Trim normalized-error cap       : %.2f sigma\n', scanData.config.trimNormCap);
  fprintf('Static MS DoA half-width        : %s deg\n', localFormatRow(scanData.config.staticMsHalfWidthDeg));
  fprintf('Selected satellites (global)    : %s\n', localFormatRow(scanData.truth.selectedSatIdxGlobal));
  scanData.plotData = localPlotScan(scanData);

  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "DONE", ...
    'subtitleText', "Single-frame static DoA-Doppler MLE/CRB consistency scan", ...
    'config', scanData.config, ...
    'snapshotFile', scanData.snapshotFile, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Static scan completed with the existing SF estimator path."; ...
      "Use full/resolved/CRB-local gaps as the dynamic comparison anchor."]));

catch ME
  if exist('progressTracker', 'var')
    localCloseScanProgressTracker(progressTracker);
  end
  if exist('progressCleanup', 'var')
    clear progressCleanup;
  end
  localPrintCheckpointFailureHint(checkpointDir);
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'subtitleText', "Single-frame static DoA-Doppler MLE/CRB consistency scan", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function taskList = localBuildTaskList(snrDbList, seedList)
%LOCALBUILDTASKLIST Build one flat SNR-repeat task list.

taskTemplate = struct('taskId', 0, 'snrDb', NaN, 'snrIndex', 0, ...
  'repeatIndex', 0, 'taskSeed', 0);
taskList = repmat(taskTemplate, numel(snrDbList) * numel(seedList), 1);
taskId = 0;
for iRepeat = 1:numel(seedList)
  for iSnr = 1:numel(snrDbList)
    taskId = taskId + 1;
    taskList(taskId).taskId = taskId;
    taskList(taskId).snrDb = snrDbList(iSnr);
    taskList(taskId).snrIndex = iSnr;
    taskList(taskId).repeatIndex = iRepeat;
    taskList(taskId).taskSeed = seedList(iRepeat);
  end
end
end

function taskOut = localRunCheckpointTask(task, scanRuntime)
%LOCALRUNCHECKPOINTTASK Run one checkpointed task.

taskOut = localRunOneTask(task, scanRuntime);
end

function taskOut = localRunOneTask(task, scanRuntime)
%LOCALRUNONETASK Run one independent static SNR-repeat task.

scanConfig = scanRuntime.config;
rng(task.taskSeed, 'twister');
pwrNoise = scanRuntime.pwrSource / (10^(task.snrDb / 10));
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  scanRuntime.steeringInfo, scanRuntime.pilotWave, scanRuntime.linkParam, ...
  scanRuntime.carrierFreq, scanRuntime.sampleRate, pwrNoise, ...
  scanRuntime.pathGain, scanRuntime.snapOpt);

viewRef = scanRuntime.viewRefTemplate;
viewRef.rxSigSf = selectRxSigBySat(rxSig, scanRuntime.refSatIdxLocal, 'singleFrame');
viewOther = scanRuntime.viewOtherTemplate;
viewOther.rxSigSf = selectRxSigBySat(rxSig, scanRuntime.otherSatIdxLocal, 'singleFrame');
viewMs = scanRuntime.viewMsTemplate;
viewMs.rxSigSf = rxSig;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRef, viewOther, viewMs, scanRuntime.wavelen, scanRuntime.pilotWave, ...
  scanRuntime.carrierFreq, scanRuntime.sampleRate, scanRuntime.fdRange, ...
  scanRuntime.truth, scanRuntime.otherSatIdxGlobal, scanConfig.optVerbose, ...
  scanRuntime.doaOnlyOpt, scanRuntime.staticOptBase, ...
  scanConfig.weightSweepAlpha, scanConfig.staticMsHalfWidthDeg);

caseList = [ ...
  caseBundle.caseRefDoa, ...
  caseBundle.caseMsDoa, ...
  caseBundle.caseStaticRefOnly, ...
  caseBundle.caseStaticMs];
truthFdList = [NaN, NaN, ...
  scanRuntime.truth.fdSatTrueHz(scanRuntime.refSatIdxLocal), ...
  scanRuntime.truth.fdRefTrueHz];
if ~isempty(caseBundle.weightCase)
  caseList = [caseList, caseBundle.weightCase]; %#ok<AGROW>
  truthFdList = [truthFdList, repmat(scanRuntime.truth.fdRefTrueHz, 1, numel(caseBundle.weightCase))]; %#ok<AGROW>
end
caseNameVec = string({caseList.displayName}).';
activeMask = ismember(caseNameVec, scanConfig.methodNameList) | startsWith(caseNameVec, "MS-SF-Static-W");
caseList = caseList(activeMask);
truthFdList = truthFdList(activeMask);

rowList = repmat(localEmptyRepeatRow(), numel(caseList), 1);
for iCase = 1:numel(caseList)
  rowList(iCase) = localBuildRepeatRow(caseList(iCase), task, truthFdList(iCase), scanRuntime.truth);
end

taskOut = struct();
taskOut.repeatTable = struct2table(rowList);
end

function scanRuntime = localBuildScanRuntime(scanConfig)
%LOCALBUILDSCANRUNTIME Build the static context shared by all MC tasks.

rng(scanConfig.baseSeed);
lightSpeed = 299792458;
osf = scanConfig.sampleRate / scanConfig.symbolRate;
wavelen = lightSpeed / scanConfig.carrierFreq;
pwrSource = 1;
E = referenceEllipsoid('sphere');
usrLla = [37.78; 36.59; 0];
truthLatlon = usrLla(1:2, 1);
searchRange = [truthLatlon(1) - scanConfig.searchMarginDeg, truthLatlon(1) + scanConfig.searchMarginDeg; ...
               truthLatlon(2) - scanConfig.searchMarginDeg, truthLatlon(2) + scanConfig.searchMarginDeg];
utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);

[~, satAccess] = findVisibleSatFromTle(utc, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccess, 2, 1, "available");
refSatIdxGlobal = satIdx(1);
scene = genMultiSatScene(utc, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal);
linkParam = getLinkParam(scene, wavelen);
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
[refState, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);
if scene.numSat ~= 2
  error('scanSfStaticMleCrbConsistency:InvalidNumSat', ...
    'This scan expects exactly two selected satellites.');
end
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

truth = struct();
truth.utc = scene.utc;
truth.latlonTrueDeg = truthLatlon;
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.usrElevationDeg = reshape(scene.accessInfo.usrElevationDeg(:, 1), 1, []);
truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
truth.deltaFdTrueHz = reshape(truthDopplerState.deltaFdRefFrame, [], 1);
truth.localDoaRef = reshape(scene.localDoa(:, refSatIdxLocal), 2, 1);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;
truth.fdRateTrueHzPerSec = NaN;
fdRange = expandRangeToTruth(scanConfig.fdRangeBase, [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], scanConfig.gridSize, searchRange, E);
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], scanConfig.gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(scene, [], scanConfig.gridSize, searchRange, E);

[pilotSym, ~] = genPilotSymbol(scanConfig.numUsr, scanConfig.numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, scanConfig.symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;
staticOptBase = struct();
staticOptBase.useLogObjective = true;

scanRuntime = struct();
scanRuntime.config = scanConfig;
scanRuntime.scene = scene;
scanRuntime.sceneRefOnly = sceneRefOnly;
scanRuntime.sceneOtherOnly = sceneOtherOnly;
scanRuntime.truth = truth;
scanRuntime.truthDopplerState = truthDopplerState;
scanRuntime.linkParam = linkParam;
scanRuntime.steeringInfo = steeringInfo;
scanRuntime.pilotWave = pilotWave;
scanRuntime.sampleRate = waveInfo.sampleRate;
scanRuntime.carrierFreq = scanConfig.carrierFreq;
scanRuntime.pwrSource = pwrSource;
scanRuntime.pathGain = pathGain;
scanRuntime.snapOpt = snapOpt;
scanRuntime.fdRange = fdRange;
scanRuntime.wavelen = wavelen;
scanRuntime.refSatIdxLocal = refSatIdxLocal;
scanRuntime.otherSatIdxLocal = otherSatIdxLocal;
scanRuntime.otherSatIdxGlobal = otherSatIdxGlobal;
scanRuntime.viewRefTemplate = viewRefTemplate;
scanRuntime.viewOtherTemplate = viewOtherTemplate;
scanRuntime.viewMsTemplate = viewMsTemplate;
scanRuntime.doaOnlyOpt = doaOnlyOpt;
scanRuntime.staticOptBase = staticOptBase;
end

function row = localBuildRepeatRow(caseInfo, task, truthFdHz, truth)
%LOCALBUILDREPEATROW Build one repeat-level metric row.

row = localEmptyRepeatRow();
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.repeatIndex = task.repeatIndex;
row.displayName = string(caseInfo.displayName);
row.satMode = string(caseInfo.satMode);
row.frameMode = string(caseInfo.frameMode);
row.paramMode = string(caseInfo.paramMode);
row.dynamicMode = string(caseInfo.dynamicMode);
row.truthFdRefHz = truthFdHz;
if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  row.tailSubtype = "missing-estimate";
  return;
end

estResult = caseInfo.estResult;
latlonEst = getDoaDopplerLatlonEst(estResult);
row.latEstDeg = latlonEst(1);
row.lonEstDeg = latlonEst(2);
if isfinite(latlonEst(1)) && isfinite(latlonEst(2))
  row.angleErrDeg = calcLatlonAngleError(latlonEst, truth.latlonTrueDeg);
end

fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
row.fdRefEstHz = fdRefEst;
if isscalar(fdRefEst) && isfinite(fdRefEst) && isfinite(truthFdHz)
  row.fdRefErrHz = fdRefEst - truthFdHz;
  row.fdRefAbsErrHz = abs(row.fdRefErrHz);
end

if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  row.isResolved = logical(estResult.isResolved);
else
  row.isResolved = isfinite(row.angleErrDeg);
end
if ~row.isResolved
  row.tailSubtype = "unresolved";
elseif ~isfinite(row.angleErrDeg)
  row.tailSubtype = "angle-nonfinite";
else
  row.tailSubtype = "resolved";
end
end

function row = localEmptyRepeatRow()
%LOCALEMPTYREPEATROW Return one typed repeat-row template.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'repeatIndex', NaN, ...
  'displayName', "", 'satMode', "", 'frameMode', "", ...
  'paramMode', "", 'dynamicMode', "", 'latEstDeg', NaN, ...
  'lonEstDeg', NaN, 'truthFdRefHz', NaN, 'fdRefEstHz', NaN, ...
  'angleErrDeg', NaN, 'fdRefErrHz', NaN, 'fdRefAbsErrHz', NaN, ...
  'isResolved', false, 'tailSubtype', "");
end

function crbTable = localBuildCrbTable(scanRuntime, scanConfig)
%LOCALBUILDCRBTABLE Evaluate the static pilot DoA-Doppler CRB over SNR.

snrDbList = reshape(scanConfig.snrDbList, [], 1);
rowTemplate = struct('snrDb', NaN, 'angleCrbSingleDeg', NaN, ...
  'angleCrbJointDeg', NaN, 'fdRefCrbSingleHz', NaN, 'fdRefCrbJointHz', NaN);
rowList = repmat(rowTemplate, numel(snrDbList), 1);
crbDdOpt = struct();
crbDdOpt.doaType = 'latlon';

fprintf('Running static CRB SNR sweep: %d points.\n', numel(snrDbList));
useParfor = localCanUseParfor(numel(snrDbList));
progressTracker = localCreateScanProgressTracker(numel(snrDbList));
progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>
if useParfor && progressTracker.enabled && ~isempty(progressTracker.queue)
  progressQueue = progressTracker.queue;
  parfor iSnr = 1:numel(snrDbList)
    rowList(iSnr) = localRunOneCrbPoint(scanRuntime, snrDbList(iSnr), crbDdOpt);
    send(progressQueue, 1);
  end
elseif useParfor
  parfor iSnr = 1:numel(snrDbList)
    rowList(iSnr) = localRunOneCrbPoint(scanRuntime, snrDbList(iSnr), crbDdOpt);
  end
else
  for iSnr = 1:numel(snrDbList)
    rowList(iSnr) = localRunOneCrbPoint(scanRuntime, snrDbList(iSnr), crbDdOpt);
    localAdvanceScanProgress(progressTracker, false);
  end
end
localCloseScanProgressTracker(progressTracker);
clear progressCleanup;
crbTable = struct2table(rowList);
end

function row = localRunOneCrbPoint(scanRuntime, snrDb, crbDdOpt)
%LOCALRUNONECRBPOINT Run one static CRB SNR point.

pwrNoise = scanRuntime.pwrSource / (10^(snrDb / 10));
[crbDdSingle, ~] = crbPilotSfDoaDoppler(scanRuntime.sceneRefOnly, ...
  scanRuntime.pilotWave, scanRuntime.carrierFreq, scanRuntime.sampleRate, ...
  scanRuntime.truth.latlonTrueDeg, scanRuntime.truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);
[crbDdJoint, ~] = crbPilotSfDoaDoppler(scanRuntime.scene, ...
  scanRuntime.pilotWave, scanRuntime.carrierFreq, scanRuntime.sampleRate, ...
  scanRuntime.truth.latlonTrueDeg, scanRuntime.truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);

row = struct();
row.snrDb = snrDb;
row.angleCrbSingleDeg = projectCrbToAngleMetric(crbDdSingle(1:2, 1:2), ...
  scanRuntime.truth.latlonTrueDeg, 'latlon');
row.angleCrbJointDeg = projectCrbToAngleMetric(crbDdJoint(1:2, 1:2), ...
  scanRuntime.truth.latlonTrueDeg, 'latlon');
row.fdRefCrbSingleHz = sqrt(crbDdSingle(3, 3));
row.fdRefCrbJointHz = sqrt(crbDdJoint(3, 3));
end

function perfTable = localAttachCrbMetrics(repeatTable, crbTable)
%LOCALATTACHCRBMETRICS Attach CRB std and normalized errors to repeat rows.

if isempty(repeatTable) || height(repeatTable) == 0
  perfTable = repeatTable;
  return;
end
perfTable = repeatTable;
perfTable.angleCrbStdDeg = nan(height(perfTable), 1);
perfTable.fdRefCrbStdHz = nan(height(perfTable), 1);
for iRow = 1:height(perfTable)
  crbIdx = find(crbTable.snrDb == perfTable.snrDb(iRow), 1, 'first');
  if isempty(crbIdx)
    continue;
  end
  if perfTable.satMode(iRow) == "multi"
    perfTable.angleCrbStdDeg(iRow) = crbTable.angleCrbJointDeg(crbIdx);
    perfTable.fdRefCrbStdHz(iRow) = crbTable.fdRefCrbJointHz(crbIdx);
  else
    perfTable.angleCrbStdDeg(iRow) = crbTable.angleCrbSingleDeg(crbIdx);
    perfTable.fdRefCrbStdHz(iRow) = crbTable.fdRefCrbSingleHz(crbIdx);
  end
end
perfTable.angleNormErr = abs(perfTable.angleErrDeg) ./ perfTable.angleCrbStdDeg;
perfTable.fdRefNormErr = abs(perfTable.fdRefAbsErrHz) ./ perfTable.fdRefCrbStdHz;
perfTable.crbLocalMask = perfTable.isResolved & isfinite(perfTable.angleNormErr);
fdFiniteMask = isfinite(perfTable.fdRefNormErr);
perfTable.crbLocalMask(fdFiniteMask) = perfTable.crbLocalMask(fdFiniteMask) & isfinite(perfTable.fdRefNormErr(fdFiniteMask));
for iRow = 1:height(perfTable)
  if perfTable.tailSubtype(iRow) ~= "resolved"
    continue;
  end
  angleTail = isfinite(perfTable.angleNormErr(iRow)) && perfTable.angleNormErr(iRow) > 5;
  fdTail = isfinite(perfTable.fdRefNormErr(iRow)) && perfTable.fdRefNormErr(iRow) > 5;
  if angleTail && fdTail
    perfTable.tailSubtype(iRow) = "angle-fd-crb-tail";
  elseif angleTail
    perfTable.tailSubtype(iRow) = "angle-crb-tail";
  elseif fdTail
    perfTable.tailSubtype(iRow) = "fdref-crb-tail";
  else
    perfTable.tailSubtype(iRow) = "crb-local";
  end
end
end

function aggregateTable = localBuildAggregateTable(perfTable, scanConfig)
%LOCALBUILDAGGREGATETABLE Build full/resolved/CRB-local aggregate metrics.

rowTemplate = localEmptyAggregateRow();
if isempty(perfTable) || height(perfTable) == 0
  aggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'satMode', 'frameMode', 'paramMode', 'dynamicMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iGroup = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iGroup) & ...
    perfTable.snrDb == groupKey.snrDb(iGroup);
  row = localBuildOneAggregateRow(perfTable, mask, scanConfig, rowTemplate);
  row.displayName = groupKey.displayName(iGroup);
  row.satMode = groupKey.satMode(iGroup);
  row.frameMode = groupKey.frameMode(iGroup);
  row.paramMode = groupKey.paramMode(iGroup);
  row.dynamicMode = groupKey.dynamicMode(iGroup);
  row.snrDb = groupKey.snrDb(iGroup);
  rowList(iGroup) = row;
end
aggregateTable = struct2table(rowList);
end

function row = localBuildOneAggregateRow(perfTable, mask, scanConfig, row)
%LOCALBUILDONEAGGREGATEROW Fill one aggregate row.

numRepeat = nnz(mask);
angleErr = perfTable.angleErrDeg(mask);
fdErr = perfTable.fdRefAbsErrHz(mask);
angleNorm = perfTable.angleNormErr(mask);
fdNorm = perfTable.fdRefNormErr(mask);
isResolved = perfTable.isResolved(mask);
row.numRepeat = numRepeat;
row.angleCrbStdDeg = median(perfTable.angleCrbStdDeg(mask), 'omitnan');
row.fdRefCrbStdHz = median(perfTable.fdRefCrbStdHz(mask), 'omitnan');
row.fullFiniteRate = mean(isfinite(angleErr));
row.resolvedRate = mean(isResolved);
row.fullAngleRmseDeg = localRmse(angleErr(isfinite(angleErr)));
row.fullAngleP95Deg = localPercentile(abs(angleErr(isfinite(angleErr))), 95);
row.fullAngleRmseOverCrb = localSafeRatio(row.fullAngleRmseDeg, row.angleCrbStdDeg);
row.resolvedAngleRmseDeg = localRmse(angleErr(isResolved & isfinite(angleErr)));
row.resolvedAngleP95Deg = localPercentile(abs(angleErr(isResolved & isfinite(angleErr))), 95);
row.resolvedAngleRmseOverCrb = localSafeRatio(row.resolvedAngleRmseDeg, row.angleCrbStdDeg);
if any(isfinite(fdErr))
  row.fullFdRefRmseHz = localRmse(fdErr(isfinite(fdErr)));
  row.fullFdRefP95Hz = localPercentile(abs(fdErr(isfinite(fdErr))), 95);
  row.fullFdRefRmseOverCrb = localSafeRatio(row.fullFdRefRmseHz, row.fdRefCrbStdHz);
  row.resolvedFdRefRmseHz = localRmse(fdErr(isResolved & isfinite(fdErr)));
  row.resolvedFdRefP95Hz = localPercentile(abs(fdErr(isResolved & isfinite(fdErr))), 95);
  row.resolvedFdRefRmseOverCrb = localSafeRatio(row.resolvedFdRefRmseHz, row.fdRefCrbStdHz);
end

crbLocalMask = isResolved & isfinite(angleNorm);
if scanConfig.trimEnable
  crbLocalMask = crbLocalMask & angleNorm <= scanConfig.trimNormCap;
  fdFinite = isfinite(fdNorm);
  crbLocalMask(fdFinite) = crbLocalMask(fdFinite) & fdNorm(fdFinite) <= scanConfig.trimNormCap;
end
row.crbLocalKeepRate = mean(crbLocalMask);
row.outlierRate = 1 - row.crbLocalKeepRate;
row.crbLocalAngleRmseDeg = localRmse(angleErr(crbLocalMask & isfinite(angleErr)));
row.crbLocalAngleP95Deg = localPercentile(abs(angleErr(crbLocalMask & isfinite(angleErr))), 95);
row.crbLocalAngleRmseOverCrb = localSafeRatio(row.crbLocalAngleRmseDeg, row.angleCrbStdDeg);
if any(isfinite(fdErr))
  row.crbLocalFdRefRmseHz = localRmse(fdErr(crbLocalMask & isfinite(fdErr)));
  row.crbLocalFdRefP95Hz = localPercentile(abs(fdErr(crbLocalMask & isfinite(fdErr))), 95);
  row.crbLocalFdRefRmseOverCrb = localSafeRatio(row.crbLocalFdRefRmseHz, row.fdRefCrbStdHz);
end
if row.crbLocalKeepRate < scanConfig.trimMinKeepRate
  row.crbLocalStatus = "low-keep-rate";
else
  row.crbLocalStatus = "ok";
end
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return one typed aggregate-row template.

row = struct('displayName', "", 'satMode', "", 'frameMode', "", ...
  'paramMode', "", 'dynamicMode', "", 'snrDb', NaN, 'numRepeat', NaN, ...
  'angleCrbStdDeg', NaN, 'fdRefCrbStdHz', NaN, 'fullFiniteRate', NaN, ...
  'resolvedRate', NaN, 'fullAngleRmseDeg', NaN, 'fullAngleP95Deg', NaN, ...
  'fullAngleRmseOverCrb', NaN, 'resolvedAngleRmseDeg', NaN, ...
  'resolvedAngleP95Deg', NaN, 'resolvedAngleRmseOverCrb', NaN, ...
  'crbLocalKeepRate', NaN, 'crbLocalAngleRmseDeg', NaN, ...
  'crbLocalAngleP95Deg', NaN, 'crbLocalAngleRmseOverCrb', NaN, ...
  'fullFdRefRmseHz', NaN, 'fullFdRefP95Hz', NaN, ...
  'fullFdRefRmseOverCrb', NaN, 'resolvedFdRefRmseHz', NaN, ...
  'resolvedFdRefP95Hz', NaN, 'resolvedFdRefRmseOverCrb', NaN, ...
  'crbLocalFdRefRmseHz', NaN, 'crbLocalFdRefP95Hz', NaN, ...
  'crbLocalFdRefRmseOverCrb', NaN, 'outlierRate', NaN, ...
  'crbLocalStatus', "");
end

function summaryTable = localBuildCrbLocalSummaryTable(aggregateTable)
%LOCALBUILDCRBLOCALSUMMARYTABLE Build compact paper-facing CRB-local rows.

if isempty(aggregateTable) || height(aggregateTable) == 0
  summaryTable = table();
  return;
end
keepVar = {'displayName', 'snrDb', 'numRepeat', 'crbLocalKeepRate', ...
  'crbLocalAngleRmseDeg', 'angleCrbStdDeg', 'crbLocalAngleRmseOverCrb', ...
  'crbLocalFdRefRmseHz', 'fdRefCrbStdHz', 'crbLocalFdRefRmseOverCrb', ...
  'resolvedRate', 'outlierRate', 'crbLocalStatus'};
summaryTable = aggregateTable(:, keepVar);
end

function topTailTable = localBuildTopTailTable(perfTable, topCountPerGroup)
%LOCALBUILDTOPTAILTABLE Select top normalized-error tail rows per group.

rowTemplate = localEmptyTopTailRow();
if isempty(perfTable) || height(perfTable) == 0 || topCountPerGroup <= 0
  topTailTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, 0, 1);
for iGroup = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iGroup) & perfTable.snrDb == groupKey.snrDb(iGroup);
  subTable = perfTable(mask, :);
  tailScore = max([subTable.angleNormErr, subTable.fdRefNormErr], [], 2, 'omitnan');
  tailScore(~isfinite(tailScore)) = -inf;
  [~, orderIdx] = sort(tailScore, 'descend');
  orderIdx = orderIdx(1:min(topCountPerGroup, numel(orderIdx)));
  for iRank = 1:numel(orderIdx)
    src = subTable(orderIdx(iRank), :);
    row = rowTemplate;
    row.displayName = src.displayName;
    row.snrDb = src.snrDb;
    row.rankInGroup = iRank;
    row.taskSeed = src.taskSeed;
    row.repeatIndex = src.repeatIndex;
    row.angleErrDeg = src.angleErrDeg;
    row.angleNormErr = src.angleNormErr;
    row.fdRefAbsErrHz = src.fdRefAbsErrHz;
    row.fdRefNormErr = src.fdRefNormErr;
    row.isResolved = src.isResolved;
    row.tailSubtype = src.tailSubtype;
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
topTailTable = struct2table(rowList);
end

function row = localEmptyTopTailRow()
%LOCALEMPTYTOPTAILROW Return one typed top-tail row.

row = struct('displayName', "", 'snrDb', NaN, 'rankInGroup', NaN, ...
  'taskSeed', NaN, 'repeatIndex', NaN, 'angleErrDeg', NaN, ...
  'angleNormErr', NaN, 'fdRefAbsErrHz', NaN, 'fdRefNormErr', NaN, ...
  'isResolved', false, 'tailSubtype', "");
end

function exportTable = localBuildTopTailExportTable(topTailTable)
%LOCALBUILDTOPTAILEXPORTTABLE Build compact top-tail export columns.

if isempty(topTailTable) || height(topTailTable) == 0
  exportTable = table();
  return;
end
keepVar = {'displayName', 'snrDb', 'rankInGroup', 'taskSeed', ...
  'angleNormErr', 'fdRefNormErr', 'tailSubtype'};
exportTable = topTailTable(:, keepVar);
end

function plotData = localBuildPlotData(aggregateTable, crbTable)
%LOCALBUILDPLOTDATA Store lightweight data needed to redraw figures.

plotData = struct();
plotData.aggregateTable = aggregateTable;
plotData.crbTable = crbTable;
end

function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw static RMSE/CRB and normalized-error figures.

plotData = scanData.plotData;
aggregateTable = scanData.aggregateTable;
if isempty(aggregateTable) || height(aggregateTable) == 0
  return;
end
methodList = unique(aggregateTable.displayName, 'stable');
crbTable = scanData.crbTable;

figure();
subplot(1, 2, 1);
hold on;
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName;
  semilogy(aggregateTable.snrDb(mask), aggregateTable.fullAngleRmseDeg(mask), '-o');
end
semilogy(crbTable.snrDb, crbTable.angleCrbSingleDeg, '-.');
semilogy(crbTable.snrDb, crbTable.angleCrbJointDeg, ':');
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Angle RMSE / CRB (deg)');
title('SF Static Angle Error');
legend([cellstr(methodList); {'SS static CRB'; 'MS static CRB'}], 'Location', 'southwest');

subplot(1, 2, 2);
hold on;
fdLegend = strings(0, 1);
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName & isfinite(aggregateTable.fullFdRefRmseHz);
  if any(mask)
    semilogy(aggregateTable.snrDb(mask), aggregateTable.fullFdRefRmseHz(mask), '-o');
    fdLegend(end + 1, 1) = methodName; %#ok<AGROW>
  end
end
semilogy(crbTable.snrDb, crbTable.fdRefCrbSingleHz, '-.');
semilogy(crbTable.snrDb, crbTable.fdRefCrbJointHz, ':');
fdLegend = [fdLegend; "SS static CRB"; "MS static CRB"];
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Reference Doppler RMSE / CRB (Hz)');
title('SF Static Reference Doppler Error');
legend(cellstr(fdLegend), 'Location', 'southwest');

figure();
subplot(1, 2, 1);
hold on;
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName;
  plot(aggregateTable.snrDb(mask), aggregateTable.crbLocalAngleRmseOverCrb(mask), '-o');
end
hold off;
grid on;
yline(1, '--');
xlabel('SNR (dB)');
ylabel('CRB-local angle RMSE / CRB');
title('SF Static CRB-local Angle Ratio');
legend(cellstr(methodList), 'Location', 'northeast');

subplot(1, 2, 2);
hold on;
fdRatioLegend = strings(0, 1);
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName & isfinite(aggregateTable.crbLocalFdRefRmseOverCrb);
  if any(mask)
    plot(aggregateTable.snrDb(mask), aggregateTable.crbLocalFdRefRmseOverCrb(mask), '-o');
    fdRatioLegend(end + 1, 1) = methodName; %#ok<AGROW>
  end
end
hold off;
grid on;
yline(1, '--');
xlabel('SNR (dB)');
ylabel('CRB-local fdRef RMSE / CRB');
title('SF Static CRB-local Doppler Ratio');
legend(cellstr(fdRatioLegend), 'Location', 'northeast');

figure();
hold on;
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName;
  plot(aggregateTable.snrDb(mask), aggregateTable.crbLocalKeepRate(mask), '-o');
end
hold off;
grid on;
ylim([0, 1.05]);
xlabel('SNR (dB)');
ylabel('CRB-local keep rate');
title('SF Static CRB-local Keep Rate');
legend(cellstr(methodList), 'Location', 'southeast');
end

function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build compact notification metrics.

metricLineList = strings(0, 1);
summaryTable = scanData.crbLocalSummaryTable;
if isempty(summaryTable) || height(summaryTable) == 0
  return;
end
snrMax = max(summaryTable.snrDb);
methodList = ["SS-SF-Static"; "MS-SF-Static"];
for iMethod = 1:numel(methodList)
  mask = summaryTable.displayName == methodList(iMethod) & summaryTable.snrDb == snrMax;
  if ~any(mask)
    continue;
  end
  row = summaryTable(find(mask, 1, 'first'), :);
  metricLineList(end + 1, 1) = sprintf('<code>%s @ %.1f dB</code>: angle %.3gx, fdRef %.3gx, keep %.3f', ...
    char(methodList(iMethod)), row.snrDb, row.crbLocalAngleRmseOverCrb, ...
    row.crbLocalFdRefRmseOverCrb, row.crbLocalKeepRate); %#ok<AGROW>
end
end

function checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, scanRuntime, taskList, useParfor)
%LOCALBUILDCHECKPOINTOPT Build checkpoint runner options under repo tmp/.

repoRoot = localGetRepoRoot();
checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, scanRuntime);
checkpointOpt.outputRoot = fullfile(repoRoot, 'tmp');
checkpointOpt.useParfor = logical(useParfor);
checkpointOpt.resume = scanConfig.checkpointResume;
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.meta = localBuildCheckpointMeta(scanConfig, scanRuntime, taskList);
end

function runKey = localBuildCheckpointRunKey(scanConfig, scanRuntime)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable, compact checkpoint key.

satToken = "sat" + strjoin(string(scanRuntime.truth.selectedSatIdxGlobal), "_");
snrToken = sprintf('snr%gto%g_n%d', scanConfig.snrDbList(1), ...
  scanConfig.snrDbList(end), numel(scanConfig.snrDbList));
methodToken = sprintf('case%d', numel(scanConfig.methodNameList) + numel(scanConfig.weightSweepAlpha));
runKey = sprintf('seed%d_%s_%s_rep%d_%s_%s', scanConfig.baseSeed, ...
  char(satToken), snrToken, scanConfig.numRepeat, methodToken, ...
  localShortHash(localBuildConfigSignature(scanConfig)));
runKey = regexprep(runKey, '[^A-Za-z0-9_\-]+', '');
end

function meta = localBuildCheckpointMeta(scanConfig, scanRuntime, taskList)
%LOCALBUILDCHECKPOINTMETA Store the task-defining scan signature.

meta = struct();
meta.scanName = string(scanConfig.scanName);
meta.snrDbList = reshape(scanConfig.snrDbList, 1, []);
meta.seedList = reshape(scanConfig.seedList, 1, []);
meta.numTask = numel(taskList);
meta.methodNameList = reshape(scanConfig.methodNameList, 1, []);
meta.selectedSatIdxGlobal = reshape(scanRuntime.truth.selectedSatIdxGlobal, 1, []);
meta.fdRange = scanRuntime.fdRange;
meta.staticMsHalfWidthDeg = reshape(scanConfig.staticMsHalfWidthDeg, 1, []);
meta.trimEnable = scanConfig.trimEnable;
meta.trimNormCap = scanConfig.trimNormCap;
meta.enableWeightSweep = scanConfig.enableWeightSweep;
meta.weightSweepAlpha = reshape(scanConfig.weightSweepAlpha, 1, []);
end

function signatureText = localBuildConfigSignature(scanConfig)
%LOCALBUILDCONFIGSIGNATURE Build text signature for stable checkpoint hashing.

signatureText = strjoin([ ...
  string(scanConfig.scanName), ...
  "seed=" + string(scanConfig.baseSeed), ...
  "snr=" + strjoin(string(scanConfig.snrDbList(:).'), ','), ...
  "repeat=" + string(scanConfig.numRepeat), ...
  "method=" + strjoin(scanConfig.methodNameList(:).', ','), ...
  "weight=" + strjoin(string(scanConfig.weightSweepAlpha(:).'), ','), ...
  "trim=" + string(scanConfig.trimNormCap)], '|');
end

function hashText = localShortHash(textValue)
%LOCALSHORTHASH Build a small deterministic non-cryptographic hash token.

textValue = char(string(textValue));
hashValue = uint32(2166136261);
for iChar = 1:numel(textValue)
  hashValue = bitxor(hashValue, uint32(textValue(iChar)));
  hashValue = uint32(mod(uint64(hashValue) * uint64(16777619), uint64(2^32)));
end
hashText = lower(dec2hex(hashValue, 8));
end

function pendingCount = localCountPendingCheckpointTasks(checkpointDir, numTask)
%LOCALCOUNTPENDINGCHECKPOINTTASKS Estimate remaining task count for progress display.

pendingCount = numTask;
if strlength(string(checkpointDir)) == 0 || ~isfolder(checkpointDir)
  return;
end
doneFile = dir(fullfile(checkpointDir, 'task', 'task_*.mat'));
pendingCount = max(0, numTask - numel(doneFile));
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runState, checkpointDir, cleanupReport)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Build a compact checkpoint status table.

if isempty(runState) || ~isstruct(runState) || ~isfield(runState, 'numTask')
  checkpointSummaryTable = table();
  return;
end
numTask = localGetFieldOrDefault(runState, 'numTask', NaN);
numDone = localGetFieldOrDefault(runState, 'numDone', NaN);
isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
runDir = string(checkpointDir);
cleanupStatus = "not-requested";
if ~isempty(cleanupReport)
  cleanupStatus = "requested";
end
checkpointSummaryTable = table(numTask, numDone, isComplete, runDir, cleanupStatus);
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print retained checkpoint path after failure.

if strlength(string(checkpointDir)) > 0
  fprintf('Checkpoint directory retained for resume/debug: %s\n', char(checkpointDir));
end
end

function repeatTable = localCollectTaskOutput(taskOutCell)
%LOCALCOLLECTTASKOUTPUT Concatenate repeat rows from task output cells.

repeatTable = table();
for iTask = 1:numel(taskOutCell)
  taskOut = taskOutCell{iTask};
  if isempty(taskOut) || ~isstruct(taskOut) || ~isfield(taskOut, 'repeatTable')
    continue;
  end
  if isempty(repeatTable)
    repeatTable = taskOut.repeatTable;
  else
    repeatTable = [repeatTable; taskOut.repeatTable]; %#ok<AGROW>
  end
end
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Return true when parallel execution is available.

if nargin < 1
  numTask = 2;
end
useParfor = false;
if numTask <= 1 || isempty(ver('parallel'))
  return;
end
try
  poolObj = gcp('nocreate');
  if isempty(poolObj)
    poolObj = parpool;
  end
  useParfor = ~isempty(poolObj);
catch
  useParfor = false;
end
end

function progressTracker = localCreateScanProgressTracker(numStep)
%LOCALCREATESCANPROGRESSTRACKER Create optional progressbar/DataQueue tracker.

progressTracker = struct('enabled', false, 'queue', [], 'numStep', max(0, numStep), ...
  'numDone', 0, 'listener', []);
if numStep <= 0 || exist('progressbar', 'file') ~= 2
  return;
end
try
  progressbar('reset', numStep);
  progressTracker.enabled = true;
  if ~isempty(ver('parallel'))
    progressTracker.queue = parallel.pool.DataQueue;
    progressTracker.listener = afterEach(progressTracker.queue, ...
      @(~) localAdvanceScanProgress(progressTracker, true)); %#ok<NASGU>
  end
catch
  progressTracker.enabled = false;
  progressTracker.queue = [];
end
end

function localAdvanceScanProgress(progressTracker, fromQueue)
%LOCALADVANCESCANPROGRESS Advance optional progressbar.

if nargin < 2
  fromQueue = false; %#ok<NASGU>
end
if ~isstruct(progressTracker) || ~isfield(progressTracker, 'enabled') || ~progressTracker.enabled
  return;
end
try
  progressbar('advance', 1);
catch
end
end

function localCloseScanProgressTracker(progressTracker)
%LOCALCLOSESCANPROGRESSTRACKER Close optional progressbar.

if ~isstruct(progressTracker) || ~isfield(progressTracker, 'enabled') || ~progressTracker.enabled
  return;
end
try
  progressbar('end');
catch
end
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve repository root from this scan location.

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(fileparts(thisFile))));
end

function tableOut = localTablePreview(tableIn, maxRows)
%LOCALTABLEPREVIEW Return a compact table preview.

if isempty(tableIn) || height(tableIn) <= maxRows
  tableOut = tableIn;
else
  tableOut = tableIn(1:maxRows, :);
end
end

function value = localRmse(x)
%LOCALRMSE Root-mean-square of finite values.

x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = sqrt(mean(x(:).^2));
end
end

function value = localPercentile(x, pct)
%LOCALPERCENTILE Compute a percentile without toolbox-specific calls.

x = sort(x(isfinite(x)));
if isempty(x)
  value = NaN;
  return;
end
if numel(x) == 1
  value = x;
  return;
end
pos = 1 + (numel(x) - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
  value = x(lo);
else
  value = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide while protecting empty or zero denominators.

if isfinite(numerator) && isfinite(denominator) && denominator > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function textValue = localFormatRow(value)
%LOCALFORMATROW Format numeric values for compact logs.

value = reshape(double(value), 1, []);
if isempty(value)
  textValue = '';
else
  textValue = strjoin(compose('%.6g', value), ', ');
end
end

function textValue = localFormatStringRow(value)
%LOCALFORMATSTRINGROW Format string values for compact logs.

value = reshape(string(value), 1, []);
if isempty(value)
  textValue = '';
else
  textValue = char(strjoin(value, ', '));
end
end
