%SCANMFMSMLECRBINTOOTHCONSISTENCY Scan MS in-tooth MLE consistency with CRB.
% Dev scan script. Edit the configuration section directly before running.
% It focuses on MS-MF-CP-K / MS-MF-CP-U over the paper SNR range and reports
% full, resolved, CRB-local, MS error-type, and representative-tail seed
% tables. Heavy bank/rescue probes stay in replay, not in this scan.

clear; close all; clc;

%% Scan configuration

scanName = "scanMfMsMleCrbInToothConsistency";
saveSnapshot = true;
notifyTelegramEnable = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
optVerbose = false;
baseSeed = 253;
numRepeat = 1000;
snrDbList = (-15:5:10).';
frameCountList = 10;
oracleFdHalfToothFraction = 0.4;
oracleFdRateHalfWidthHzPerSec = 1000;
staticLocalDoaHalfWidthDeg = [0.002; 0.002];
initMode = "auto";
% Use "auto" for each MF estimator to build its own MUSIC/Bartlett DoA and
% reference-satellite fd-line initializer inside the in-tooth frequency box.
% Use "staticTruthFreq" only as a diagnostic oracle baseline that seeds DoA
% from the static transition bundle and seeds fdRef/fdRate at truth.
resolvedToothHalfWidthFraction = 0.25;
resolvedFdRateAbsTolHzPerSec = 250;
boundaryTolFraction = 0.02;
coherenceCollapseFloor = 0.25;
coreCoherenceFloor = 0.8;
trimEnable = true;
trimNormCap = 5;
trimMinKeepRate = 0.8;
topTailCountPerGroup = 5;
representativeSeedCountPerType = 3;
topTailPrintMaxRows = 60;
diagnosticVersion = "ms-crb-normalized-error-type-v2";
methodNameList = [
  "MS-MF-CP-K"
  "MS-MF-CP-U"
];
% This MS-focused scan intentionally omits SS methods. Use
% scanMfMleCrbInToothConsistency when a mixed SS/MS comparison is needed.

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
frameCountList = reshape(double(frameCountList), [], 1);
staticLocalDoaHalfWidthDeg = reshape(double(staticLocalDoaHalfWidthDeg), [], 1);
initMode = string(initMode);
methodNameList = reshape(string(methodNameList), [], 1);
methodNameList = methodNameList(strlength(methodNameList) > 0);
numRepeat = numel(seedList);

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.baseSeed = baseSeed;
scanConfig.numRepeat = numRepeat;
scanConfig.seedList = seedList;
scanConfig.snrDbList = snrDbList;
scanConfig.frameCountList = frameCountList;
scanConfig.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
scanConfig.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
scanConfig.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg;
scanConfig.initMode = initMode;
scanConfig.resolvedToothHalfWidthFraction = resolvedToothHalfWidthFraction;
scanConfig.resolvedFdRateAbsTolHzPerSec = resolvedFdRateAbsTolHzPerSec;
scanConfig.boundaryTolFraction = boundaryTolFraction;
scanConfig.coherenceCollapseFloor = coherenceCollapseFloor;
scanConfig.coreCoherenceFloor = coreCoherenceFloor;
scanConfig.trimEnable = logical(trimEnable);
scanConfig.trimNormCap = trimNormCap;
scanConfig.trimMinKeepRate = trimMinKeepRate;
scanConfig.topTailCountPerGroup = topTailCountPerGroup;
scanConfig.representativeSeedCountPerType = representativeSeedCountPerType;
scanConfig.topTailPrintMaxRows = topTailPrintMaxRows;
scanConfig.diagnosticVersion = string(diagnosticVersion);
scanConfig.methodNameList = methodNameList;
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.checkpointResume = logical(checkpointResume);
scanConfig.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);
scanConfig.optVerbose = logical(optVerbose);
localValidateInitMode(scanConfig.initMode);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
checkpointDir = "";
checkpointRunState = struct();
checkpointCleanupReport = struct([]);
scanData = struct();

try
  %% Build context and scan tasks

  taskList = localBuildTaskList(frameCountList, snrDbList, seedList);
  scanConfig.numTask = numel(taskList);
  useParfor = localCanUseParfor(numel(taskList));
  if scanConfig.checkpointEnable
    checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, taskList, useParfor);
    checkpointDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
  else
    checkpointOpt = struct();
  end
  printMfScanHeader(char(scanName), scanConfig, checkpointDir);
  fprintf('Preparing runtime cache...\n');
  scanRuntime = localBuildScanRuntime(scanConfig);

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

  [perfTable, repeatOutCell] = localCollectTaskOutput(taskOutCell);
  aggregateTable = localBuildAggregateTable(perfTable, scanConfig);
  crbLocalSummaryTable = localBuildCrbLocalSummaryTable(aggregateTable);
  failureSummaryTable = localBuildFailureSummaryTable(perfTable);
  msErrorTypeSummaryTable = localBuildMsErrorTypeSummaryTable(perfTable, scanConfig);
  representativeSeedTable = localBuildRepresentativeSeedTable(perfTable, scanConfig);
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
  scanData.perfTable = perfTable;
  scanData.aggregateTable = aggregateTable;
  scanData.crbLocalSummaryTable = crbLocalSummaryTable;
  scanData.failureSummaryTable = failureSummaryTable;
  scanData.msErrorTypeSummaryTable = msErrorTypeSummaryTable;
  scanData.representativeSeedTable = representativeSeedTable;
  scanData.topTailTable = topTailTable;
  scanData.topTailExportTable = topTailExportTable;
  scanData.repeatOutCell = repeatOutCell;
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupReport = checkpointCleanupReport;
  scanData.checkpointDir = checkpointDir;
  scanData.plotData = localBuildPlotData(aggregateTable, failureSummaryTable, topTailTable);
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

  if ~exist('scanData', 'var') || ~isstruct(scanData)
    error('scanMfMsMleCrbInToothConsistency:MissingScanData', ...
      'scanData is missing. Load a snapshot containing scanData before running this section.');
  end
  scanData = localValidateScanDataForSummary(scanData);
  scanNameForReport = string(localGetFieldOrDefault(scanData, 'scanName', "scanMfMsMleCrbInToothConsistency"));
  scanConfigForReport = localGetFieldOrDefault(scanData, 'config', struct());
  scanSnapshotFile = string(localGetFieldOrDefault(scanData, 'snapshotFile', ""));
  scanCheckpointDir = string(localGetFieldOrDefault(scanData, 'checkpointDir', ""));
  scanElapsedSec = localGetFieldOrDefault(scanData, 'elapsedSec', NaN);
  topTailPrintMaxRowsForReport = localGetFieldOrDefault(scanConfigForReport, 'topTailPrintMaxRows', 60);

  printMfScanSection('MS in-tooth MLE/CRB consistency aggregate', scanData.aggregateTable);
  printMfScanSection('MS CRB-normalized compact summary', scanData.crbLocalSummaryTable);
  printMfScanSection('MS error type summary', scanData.msErrorTypeSummaryTable);
  printMfScanSection('MS representative tail seeds', scanData.representativeSeedTable);
  printMfScanSection('Resolved failure reason summary', scanData.failureSummaryTable);
  if height(scanData.topTailTable) > 0
    printMfScanSection('Top CRB-normalized tail preview', ...
      localTablePreview(scanData.topTailTable, topTailPrintMaxRowsForReport));
    printMfScanSection('Top-tail compact export preview', ...
      localTablePreview(scanData.topTailExportTable, topTailPrintMaxRowsForReport));
  end
  fprintf('Frames                           : %s\n', localFormatRow(localGetFieldOrDefault(scanConfigForReport, 'frameCountList', [])));
  fprintf('SNR list (dB)                    : %s\n', localFormatRow(localGetFieldOrDefault(scanConfigForReport, 'snrDbList', [])));
  fprintf('Active method list               : %s\n', localFormatStringRow(localGetFieldOrDefault(scanConfigForReport, 'methodNameList', strings(0, 1))));
  fprintf('Runtime cache                    : context per frame, CRB per frame/SNR.\n');
  fprintf('Dynamic init mode                : %s\n', char(string(localGetFieldOrDefault(scanConfigForReport, 'initMode', ""))));
  fprintf('Diagnostic version               : %s\n', char(string(localGetFieldOrDefault(scanConfigForReport, 'diagnosticVersion', ""))));
  fprintf('Representative seeds per type    : %d\n', localGetFieldOrDefault(scanConfigForReport, 'representativeSeedCountPerType', 0));
  fprintf('Repeats per config               : %d\n', localGetFieldOrDefault(scanConfigForReport, 'numRepeat', NaN));
  fprintf('fdRef oracle half-tooth fraction : %.3f\n', localGetFieldOrDefault(scanConfigForReport, 'oracleFdHalfToothFraction', NaN));
  fprintf('fdRate oracle half-width         : %.2f Hz/s\n', localGetFieldOrDefault(scanConfigForReport, 'oracleFdRateHalfWidthHzPerSec', NaN));
  fprintf('Resolved tooth half-width        : %.3f tooth\n', localGetFieldOrDefault(scanConfigForReport, 'resolvedToothHalfWidthFraction', NaN));
  fprintf('Resolved fdRate tolerance        : %.2f Hz/s\n', localGetFieldOrDefault(scanConfigForReport, 'resolvedFdRateAbsTolHzPerSec', NaN));
  fprintf('Core coherence floor             : %.3f\n', localGetFieldOrDefault(scanConfigForReport, 'coreCoherenceFloor', NaN));
  fprintf('Trim normalized-error cap        : %.2f sigma\n', localGetFieldOrDefault(scanConfigForReport, 'trimNormCap', NaN));
  fprintf('Trim minimum keep rate           : %.3f\n', localGetFieldOrDefault(scanConfigForReport, 'trimMinKeepRate', NaN));
  fprintf('Resolved rule                    : solver-valid, in-tooth, no frequency-boundary hit, fdRate healthy for CP-U.\n');
  fprintf('Core rule                        : resolved plus non-ref coherence floor for multi-sat cases.\n');
  fprintf('Trim rule                        : core plus fixed CRB-normalized angle/fdRef cap; report keep rate separately.\n');
  fprintf('CRB-local rule                   : alias of trimmed core; use as paper-facing local-regime subset.\n');
  fprintf('CRB-normalized summary           : reports RMSE/CRB and MSE/CRB first; keep rates remain diagnostic.\n');
  fprintf('Angle error is report-only and is not used to define loose resolved samples.\n');
  if height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  scanData.plotData = localPlotScan(scanData);

  notifyMfScanStatus(struct( ...
    'scanName', scanNameForReport, ...
    'statusText', "DONE", ...
    'config', scanConfigForReport, ...
    'snapshotFile', scanSnapshotFile, ...
    'checkpointDir', scanCheckpointDir, ...
    'elapsedSec', scanElapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "MS in-tooth MLE-vs-CRB error-type scan completed."; ...
      "Use representative seeds in replayMfMsMleCrbFlowDiagnose for mechanism probes."]));

catch ME
  if exist('progressTracker', 'var')
    localCloseScanProgressTracker(progressTracker);
  end
  if exist('progressCleanup', 'var')
    clear progressCleanup;
  end
  localPrintCheckpointFailureHint(checkpointDir);
  fprintf('Scan failed while building in-tooth MLE/CRB consistency data.\n');
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function taskOut = localRunCheckpointTask(task, scanRuntime)
%LOCALRUNCHECKPOINTTASK Run one checkpointed scan task.

taskOut = localRunOneTask(task, scanRuntime);
end

function checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, taskList, useParfor)
%LOCALBUILDCHECKPOINTOPT Build task-level checkpoint runner options.

checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, taskList);
checkpointOpt.outputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.useParfor = logical(useParfor);
checkpointOpt.resume = logical(scanConfig.checkpointResume);
checkpointOpt.meta = localBuildCheckpointMeta(scanConfig, taskList);
checkpointOpt.progressFcn = [];
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end

function runKey = localBuildCheckpointRunKey(scanConfig, taskList)
%LOCALBUILDCHECKPOINTRUNKEY Build a compact stable key for interrupted scans.

frameList = unique([taskList.numFrame], 'stable');
snrList = unique([taskList.snrDb], 'stable');
seedList = unique([taskList.taskSeed], 'stable');
signature = localBuildCheckpointSignature(scanConfig, taskList);
signatureHash = localCompactStringHash(signature);
runKey = sprintf('f%s_snr%.6gto%.6g_seed%dto%d_rep%d_%s', ...
  localJoinNumericToken(frameList), snrList(1), snrList(end), ...
  seedList(1), seedList(end), numel(seedList), signatureHash);
runKey = localPathSafeToken(runKey);
end

function signature = localBuildCheckpointSignature(scanConfig, taskList)
%LOCALBUILDCHECKPOINTSIGNATURE Build the full semantic signature behind runKey.

frameList = reshape(unique([taskList.numFrame], 'stable'), 1, []);
snrList = reshape(unique([taskList.snrDb], 'stable'), 1, []);
seedList = reshape(unique([taskList.taskSeed], 'stable'), 1, []);
signature = strjoin(string({ ...
  sprintf('frames%s', mat2str(frameList, 8)), ...
  sprintf('snr%s', mat2str(snrList, 8)), ...
  sprintf('seed%dto%d_rep%d', seedList(1), seedList(end), numel(seedList)), ...
  sprintf('oracleFdFrac%.8g', scanConfig.oracleFdHalfToothFraction), ...
  sprintf('oracleFdRate%.8g', scanConfig.oracleFdRateHalfWidthHzPerSec), ...
  sprintf('resolvedToothFrac%.8g', scanConfig.resolvedToothHalfWidthFraction), ...
  sprintf('resolvedFdRate%.8g', scanConfig.resolvedFdRateAbsTolHzPerSec), ...
  sprintf('boundaryTol%.8g', scanConfig.boundaryTolFraction), ...
  sprintf('cohCollapse%.8g', scanConfig.coherenceCollapseFloor), ...
  sprintf('coreCoh%.8g', scanConfig.coreCoherenceFloor), ...
  sprintf('trimEnable%d', scanConfig.trimEnable), ...
  sprintf('trimCap%.8g', scanConfig.trimNormCap), ...
  sprintf('trimMin%.8g', scanConfig.trimMinKeepRate), ...
  sprintf('topTail%d', scanConfig.topTailCountPerGroup), ...
  sprintf('reprSeed%d', scanConfig.representativeSeedCountPerType), ...
  sprintf('diag%s', char(string(scanConfig.diagnosticVersion))), ...
  sprintf('methods%s', char(strjoin(reshape(string(scanConfig.methodNameList), 1, []), ','))), ...
  sprintf('init%s', char(string(scanConfig.initMode))), ...
  sprintf('doaHalf%s', mat2str(reshape(scanConfig.staticLocalDoaHalfWidthDeg(:), 1, []), 8)) ...
  }), '|');
end

function hashText = localCompactStringHash(textValue)
%LOCALCOMPACTSTRINGHASH Return a short deterministic hash for path-safe names.

byteValues = double(char(textValue));
hashValue = 5381;
modValue = 2147483647;
for iByte = 1:numel(byteValues)
  hashValue = mod(hashValue * 33 + byteValues(iByte), modValue);
end
hashText = sprintf('%08x', uint32(hashValue));
end

function token = localPathSafeToken(value)
%LOCALPATHSAFETOKEN Convert a short semantic token into a path-safe string.

token = string(value);
token = replace(token, '.', 'p');
token = replace(token, '-', 'm');
token = replace(token, ' ', '');
token = replace(token, '+', '');
token = regexprep(char(token), '[^A-Za-z0-9_]', '');
token = string(token);
end

function localValidateInitMode(initMode)
%LOCALVALIDATEINITMODE Validate the scan-level dynamic initializer mode.

validModeList = ["auto"; "staticTruthFreq"];
if ~any(string(initMode) == validModeList)
  error('scanMfMsMleCrbInToothConsistency:InvalidInitMode', ...
    'initMode must be one of: %s.', char(strjoin(validModeList, ', ')));
end
end

function meta = localBuildCheckpointMeta(scanConfig, taskList)
%LOCALBUILDCHECKPOINTMETA Store the task-defining scan signature.

meta = struct();
meta.scanName = string(scanConfig.scanName);
meta.frameCountList = reshape(unique([taskList.numFrame], 'stable'), 1, []);
meta.snrDbList = reshape(unique([taskList.snrDb], 'stable'), 1, []);
meta.seedList = reshape(unique([taskList.taskSeed], 'stable'), 1, []);
meta.oracleFdHalfToothFraction = scanConfig.oracleFdHalfToothFraction;
meta.oracleFdRateHalfWidthHzPerSec = scanConfig.oracleFdRateHalfWidthHzPerSec;
meta.resolvedToothHalfWidthFraction = scanConfig.resolvedToothHalfWidthFraction;
meta.resolvedFdRateAbsTolHzPerSec = scanConfig.resolvedFdRateAbsTolHzPerSec;
meta.staticLocalDoaHalfWidthDeg = reshape(scanConfig.staticLocalDoaHalfWidthDeg(:), 1, []);
meta.initMode = string(scanConfig.initMode);
meta.coreCoherenceFloor = scanConfig.coreCoherenceFloor;
meta.trimEnable = scanConfig.trimEnable;
meta.trimNormCap = scanConfig.trimNormCap;
meta.trimMinKeepRate = scanConfig.trimMinKeepRate;
meta.representativeSeedCountPerType = scanConfig.representativeSeedCountPerType;
meta.methodNameList = reshape(string(scanConfig.methodNameList), 1, []);
meta.diagnosticVersion = string(scanConfig.diagnosticVersion);
meta.numTask = numel(taskList);
end

function pendingCount = localCountPendingCheckpointTasks(checkpointDir, numTask)
%LOCALCOUNTPENDINGCHECKPOINTTASKS Estimate remaining task count for progress display.

pendingCount = numTask;
if strlength(string(checkpointDir)) == 0 || ~isfolder(checkpointDir)
  return;
end
doneFile = dir(fullfile(checkpointDir, 'task', 'task_*.mat'));
if isempty(doneFile)
  doneFile = dir(fullfile(checkpointDir, 'task_*.mat'));
end
pendingCount = max(0, numTask - numel(doneFile));
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runState, checkpointDir, cleanupReport)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Build a compact checkpoint status table.

checkpointSummaryTable = table();
if isempty(runState) || ~isstruct(runState) || ~isfield(runState, 'manifest')
  return;
end
numTask = localGetFieldOrDefault(runState.manifest, 'numTask', NaN);
numDone = numel(localGetFieldOrDefault(runState, 'resultCell', {}));
isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
runDir = string(checkpointDir);
cleanupStatus = "not-requested";
if ~isempty(cleanupReport)
  cleanupStatus = "requested";
end
checkpointSummaryTable = table(numTask, numDone, isComplete, runDir, cleanupStatus);
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print the retained checkpoint path after failures.

if strlength(string(checkpointDir)) > 0
  fprintf('Checkpoint directory retained for resume/debug: %s\n', char(checkpointDir));
end
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve repository root from this scan location.

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(fileparts(thisFile))));
end

function token = localJoinNumericToken(value)
%LOCALJOINNUMERICTOKEN Format a numeric vector for stable run keys.

value = reshape(double(value), 1, []);
token = strjoin(compose('%.6g', value), '_');
end

function token = localJoinStringToken(value)
%LOCALJOINSTRINGTOKEN Format string labels for stable run keys.

value = reshape(string(value), 1, []);
if isempty(value)
  token = 'none';
else
  token = char(strjoin(value, '_'));
  token = regexprep(token, '[^A-Za-z0-9_]+', '');
end
end

function taskList = localBuildTaskList(frameCountList, snrDbList, seedList)
%LOCALBUILDTASKLIST Build one flat task list from frame, SNR, and seed axes.

taskList = repmat(struct('taskId', 0, 'numFrame', 0, 'snrDb', NaN, 'taskSeed', 0), 0, 1);
taskId = 0;
for iFrame = 1:numel(frameCountList)
  for iSnr = 1:numel(snrDbList)
    for iSeed = 1:numel(seedList)
      taskId = taskId + 1;
      taskList(taskId, 1).taskId = taskId; %#ok<AGROW>
      taskList(taskId, 1).numFrame = frameCountList(iFrame); %#ok<AGROW>
      taskList(taskId, 1).snrDb = snrDbList(iSnr); %#ok<AGROW>
      taskList(taskId, 1).taskSeed = seedList(iSeed); %#ok<AGROW>
    end
  end
end
end

function taskOut = localRunOneTask(task, scanRuntime)
%LOCALRUNONETASK Run one frame-count, SNR, and seed configuration.

scanConfig = scanRuntime.config;
context = localGetRuntimeContext(scanRuntime, task.numFrame);
flowOpt = localGetRuntimeFlowOpt(scanRuntime, task.numFrame);
methodList = scanRuntime.methodList;
crbTable = localGetRuntimeCrbTable(scanRuntime, task.numFrame, task.snrDb);
repeatOut = localRunOneRepeat(context, flowOpt, methodList, scanConfig, ...
  task.snrDb, task.taskSeed, crbTable);

taskOut = struct();
taskOut.caseRowList = repeatOut.caseRowList;
taskOut.repeatSlim = localStripRepeatOut(repeatOut);
end

function scanRuntime = localBuildScanRuntime(scanConfig)
%LOCALBUILDSCANRUNTIME Cache frame-level context and CRB tables outside task workers.

localValidateInitMode(scanConfig.initMode);
frameCountList = reshape(unique(scanConfig.frameCountList, 'stable'), [], 1);
contextBank = repmat(struct('numFrame', NaN, 'context', struct(), 'flowOpt', struct()), numel(frameCountList), 1);
crbTableBank = repmat(struct('numFrame', NaN, 'snrDb', NaN, 'crbTable', table()), 0, 1);
parallelOpt = localBuildSerialInnerParallelOpt();
for iFrame = 1:numel(frameCountList)
  numFrame = frameCountList(iFrame);
  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = parallelOpt;
  context = buildDynamicDualSatEciContext(contextOpt);
  context = localDisableSubsetBankForInTooth(context);
  flowOpt = localBuildFlowOpt(scanConfig.staticLocalDoaHalfWidthDeg, context.parallelOpt);
  contextBank(iFrame).numFrame = numFrame;
  contextBank(iFrame).context = context;
  contextBank(iFrame).flowOpt = flowOpt;

  refRepeat = buildDynamicRepeatData(context, 0, scanConfig.baseSeed);
  refFixture = refRepeat.periodicFixture;
  for iSnr = 1:numel(scanConfig.snrDbList)
    snrDb = scanConfig.snrDbList(iSnr);
    noiseVar = 1 / (10^(snrDb / 10));
    crbBundle = localBuildCrbBundleQuiet(refFixture, context.pilotWave, ...
      context.carrierFreq, context.waveInfo.sampleRate, noiseVar);
    crbTable = buildDynamicCrbSummaryTable(crbBundle.truth, ...
      crbBundle.crbSfRef, crbBundle.auxCrbSfRef, ...
      crbBundle.crbSfMs, crbBundle.auxCrbSfMs, ...
      crbBundle.crbMfRefKnown, crbBundle.auxCrbMfRefKnown, ...
      crbBundle.crbMfMsKnown, crbBundle.auxCrbMfMsKnown, ...
      crbBundle.crbMfRefUnknown, crbBundle.auxCrbMfRefUnknown, ...
      crbBundle.crbMfMsUnknown, crbBundle.auxCrbMfMsUnknown, snrDb);
    crbTableBank(end + 1, 1).numFrame = numFrame; %#ok<AGROW>
    crbTableBank(end, 1).snrDb = snrDb;
    crbTableBank(end, 1).crbTable = crbTable;
  end
end

scanRuntime = struct();
scanRuntime.config = scanConfig;
scanRuntime.contextBank = contextBank;
scanRuntime.crbTableBank = crbTableBank;
scanRuntime.methodList = localBuildMethodList(scanConfig.staticLocalDoaHalfWidthDeg, scanConfig.methodNameList);
end

function context = localGetRuntimeContext(scanRuntime, numFrame)
%LOCALGETRUNTIMECONTEXT Return the cached context for one frame count.

idx = find([scanRuntime.contextBank.numFrame] == numFrame, 1, 'first');
if isempty(idx)
  error('scanMfMsMleCrbInToothConsistency:MissingRuntimeContext', ...
    'No cached runtime context found for numFrame=%g.', numFrame);
end
context = scanRuntime.contextBank(idx).context;
end

function flowOpt = localGetRuntimeFlowOpt(scanRuntime, numFrame)
%LOCALGETRUNTIMEFLOWOPT Return the cached flow options for one frame count.

idx = find([scanRuntime.contextBank.numFrame] == numFrame, 1, 'first');
if isempty(idx)
  error('scanMfMsMleCrbInToothConsistency:MissingRuntimeFlowOpt', ...
    'No cached runtime flow options found for numFrame=%g.', numFrame);
end
flowOpt = scanRuntime.contextBank(idx).flowOpt;
end

function crbTable = localGetRuntimeCrbTable(scanRuntime, numFrame, snrDb)
%LOCALGETRUNTIMECRBTABLE Return the cached CRB summary for one frame/SNR pair.

frameVec = [scanRuntime.crbTableBank.numFrame];
snrVec = [scanRuntime.crbTableBank.snrDb];
idx = find(frameVec == numFrame & abs(snrVec - snrDb) < 1e-12, 1, 'first');
if isempty(idx)
  error('scanMfMsMleCrbInToothConsistency:MissingRuntimeCrbTable', ...
    'No cached CRB table found for numFrame=%g, snrDb=%g.', numFrame, snrDb);
end
crbTable = scanRuntime.crbTableBank(idx).crbTable;
end

function parallelOpt = localBuildSerialInnerParallelOpt()
%LOCALBUILDSERIALINNERPARALLELOPT Disable inner parallelism for outer parfor scans.

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
end

function context = localDisableSubsetBankForInTooth(context)
%LOCALDISABLESUBSETBANKFORINTOOTH Remove subset fixtures from this local scan.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function flowOpt = localBuildFlowOpt(staticLocalDoaHalfWidthDeg, parallelOpt)
%LOCALBUILDFLOWOPT Build static and dynamic options shared by local methods.

flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false, ...
  'staticMsHalfWidth', staticLocalDoaHalfWidthDeg(:)));
flowOpt = applyDynamicTransitionFlowDefaults(flowOpt);
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if numFrame < 2
  error('scanMfMsMleCrbInToothConsistency:InvalidFrameCount', ...
    'Frame count must be at least two for multi-frame consistency scanning.');
end
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function methodList = localBuildMethodList(staticLocalDoaHalfWidthDeg, methodNameList)
%LOCALBUILDMETHODLIST Define and filter the tooth-resolved CP-K / CP-U method bank.

fullBank = repmat(localMakeMethod("", "", "", "", false, "", "", []), 0, 1);
fullBank(end + 1, 1) = localMakeMethod("SS-MF-CP-K", "ref", "continuous", "known", true, ...
  "truthHalfTooth", "known", staticLocalDoaHalfWidthDeg);
fullBank(end + 1, 1) = localMakeMethod("SS-MF-CP-U", "ref", "continuous", "unknown", false, ...
  "truthHalfTooth", "truthLocal", staticLocalDoaHalfWidthDeg);
fullBank(end + 1, 1) = localMakeMethod("MS-MF-CP-K", "ms", "continuous", "known", true, ...
  "truthHalfTooth", "known", staticLocalDoaHalfWidthDeg);
fullBank(end + 1, 1) = localMakeMethod("MS-MF-CP-U", "ms", "continuous", "unknown", false, ...
  "truthHalfTooth", "truthLocal", staticLocalDoaHalfWidthDeg);

requestedNameList = reshape(string(methodNameList), [], 1);
requestedNameList = requestedNameList(strlength(requestedNameList) > 0);
availableNameList = reshape([fullBank.displayName], [], 1);
if isempty(requestedNameList)
  error('scanMfMsMleCrbInToothConsistency:EmptyMethodNameList', ...
    'methodNameList must contain at least one method. Available methods: %s.', ...
    char(strjoin(availableNameList, ', ')));
end
if numel(unique(requestedNameList, 'stable')) ~= numel(requestedNameList)
  error('scanMfMsMleCrbInToothConsistency:DuplicateMethodName', ...
    'methodNameList contains duplicate method names: %s.', char(strjoin(requestedNameList, ', ')));
end
[isKnownMethod, methodIdx] = ismember(requestedNameList, availableNameList);
if any(~isKnownMethod)
  error('scanMfMsMleCrbInToothConsistency:UnknownMethodName', ...
    'Unknown methodNameList entry %s. Available methods: %s.', ...
    char(strjoin(requestedNameList(~isKnownMethod), ', ')), char(strjoin(availableNameList, ', ')));
end
methodList = fullBank(methodIdx);
end

function method = localMakeMethod(displayName, viewMode, phaseMode, fdRateMode, isKnownRate, fdRangeMode, fdRateRangeMode, doaHalfWidthDeg)
%LOCALMAKEMETHOD Construct one local MLE method descriptor.

method = struct();
method.displayName = string(displayName);
method.viewMode = string(viewMode);
method.satMode = "multi";
if method.viewMode == "ref"
  method.satMode = "single";
end
method.phaseMode = string(phaseMode);
method.fdRateMode = string(fdRateMode);
method.isKnownRate = logical(isKnownRate);
method.fdRangeMode = string(fdRangeMode);
method.fdRateRangeMode = string(fdRateRangeMode);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
end

function repeatOut = localRunOneRepeat(context, flowOpt, methodList, scanConfig, snrDb, taskSeed, crbTable)
%LOCALRUNONEREPEAT Run one noisy repeat and collect MLE/CRB rows.

lastwarn('', '');
tRepeat = tic;
repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('scanMfMsMleCrbInToothConsistency:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for the controlled in-tooth scan.');
end
fdHalfWidthHz = scanConfig.oracleFdHalfToothFraction * toothStepHz;
fdRangeOracle = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeOracle = truthFdRateHzPerSec + scanConfig.oracleFdRateHalfWidthHzPerSec * [-1, 1];

if scanConfig.initMode == "staticTruthFreq"
  staticBundle = buildDoaDopplerStaticTransitionBundle( ...
    fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
    context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
    fixture.fdRange, truth, context.otherSatIdxGlobal, scanConfig.optVerbose, ...
    flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);
else
  staticBundle = struct();
end

caseRowList = repmat(localEmptyCaseRow(), numel(methodList), 1);
for iMethod = 1:numel(methodList)
  method = methodList(iMethod);
  tMethod = tic;
  caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, context, flowOpt, ...
    fdRangeOracle, fdRateRangeOracle, toothStepHz, scanConfig.optVerbose, scanConfig.initMode);
  wallTimeMs = 1000 * toc(tMethod);
  caseRowList(iMethod) = localBuildCaseRow(method, caseUse, truth, snrDb, taskSeed, ...
    fixture, toothStepHz, wallTimeMs, crbTable, scanConfig);
end
repeatTotalMs = 1000 * toc(tRepeat);
[warningMessage, warningId] = lastwarn;

repeatOut = struct();
repeatOut.taskSeed = taskSeed;
repeatOut.snrDb = snrDb;
repeatOut.truth = truth;
repeatOut.toothStepHz = toothStepHz;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.fdRangeOracle = fdRangeOracle;
repeatOut.fdRateRangeOracle = fdRateRangeOracle;
repeatOut.crbTable = crbTable;
repeatOut.caseRowList = caseRowList;
repeatOut.repeatTotalMs = repeatTotalMs;
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function crbBundle = localBuildCrbBundleQuiet(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%LOCALBUILDCRBBUNDLEQUIET Build CRB while suppressing expected full-FIM warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
end

function caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, context, ...
  flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz, optVerbose, initMode)
%LOCALRUNDYNAMICMETHOD Run one tooth-resolved dynamic method.

[viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle, initMode);
fdRangeUse = fixture.fdRange;
fdRateRangeUse = fixture.fdRateRange;
if method.fdRangeMode == "truthHalfTooth"
  fdRangeUse = fdRangeOracle;
end
if method.fdRateRangeMode == "truthLocal"
  fdRateRangeUse = fdRateRangeOracle;
end

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.enableFdAliasUnwrap = true;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
initOverride = [];
if initMode == "staticTruthFreq"
  truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
  truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
  seedDoaParam = staticCaseUse.estResult.doaParamEst(:);
  dynOpt.initDoaParam = seedDoaParam(:);
  dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
  initParam = localBuildTruthFreqInitParam(seedDoaParam, truthFdRefHz, truthFdRateHzPerSec, method.isKnownRate);
  initOverride = struct('startTag', lower(method.displayName) + "-static-truth-freq", ...
    'initParam', initParam, 'initDoaParam', seedDoaParam(:), ...
    'initDoaHalfWidth', method.doaHalfWidthDeg(:), 'freezeDoa', false);
end

caseUse = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
  viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
  debugTruthUse, initOverride);
caseUse.inToothMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, ...
  'fdRateRangeUse', fdRateRangeUse);
end

function [viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle, initMode)
%LOCALRESOLVEMETHODVIEW Resolve the view and optional static DoA seed for one method.

staticCaseUse = struct();
switch method.viewMode
  case "ref"
    viewUse = fixture.viewRefOnly;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());
    if initMode == "staticTruthFreq"
      staticCaseUse = staticBundle.caseStaticRefOnly;
    end
  case "ms"
    viewUse = fixture.viewMs;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());
    if initMode == "staticTruthFreq"
      staticCaseUse = staticBundle.caseStaticMs;
    end
  otherwise
    error('scanMfMsMleCrbInToothConsistency:UnknownViewMode', ...
      'Unknown method viewMode "%s".', method.viewMode);
end
end

function initParam = localBuildTruthFreqInitParam(doaParam, fdRefHz, fdRateHzPerSec, isKnownRate)
%LOCALBUILDTRUTHFREQINITPARAM Build an oracle frequency initializer.

if isKnownRate
  initParam = [doaParam(:); fdRefHz];
else
  initParam = [doaParam(:); fdRefHz; fdRateHzPerSec];
end
end

function row = localBuildCaseRow(method, caseUse, truth, snrDb, taskSeed, fixture, toothStepHz, wallTimeMs, crbTable, scanConfig)
%LOCALBUILDCASEMETRICROW Convert one estimator case into a compact scan row.

row = localEmptyCaseRow();
row.displayName = method.displayName;
row.satMode = method.satMode;
row.phaseMode = method.phaseMode;
row.fdRateMode = method.fdRateMode;
row.snrDb = snrDb;
row.taskSeed = taskSeed;
row.numFrame = fixture.sceneSeq.numFrame;
row.frameIntvlSec = median(diff(fixture.sceneSeq.timeOffsetSec));
row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
row.toothStepHz = toothStepHz;
row.inToothMode = method.fdRangeMode;
row.fdRateRangeMode = method.fdRateRangeMode;
if scanConfig.initMode == "staticTruthFreq"
  row.doaSeedMode = "static";
  row.initMode = "static-truth-freq";
else
  row.doaSeedMode = "auto";
  row.initMode = "auto-mf";
end
info = summarizeDoaDopplerCase(caseUse, truth);
dynSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
initEval = localGetFieldOrDefault(estDebug, 'initEval', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
row.solverResolved = logical(localGetFieldOrDefault(info, 'isResolved', false));
row.angleErrDeg = localGetFieldOrDefault(info, 'angleErrDeg', NaN);
[latlonEstDeg, truthLatlonDeg] = localResolveCaseLatlon(caseUse, truth);
row.latEstDeg = latlonEstDeg(1);
row.lonEstDeg = latlonEstDeg(2);
row.latTrueDeg = truthLatlonDeg(1);
row.lonTrueDeg = truthLatlonDeg(2);
if all(isfinite(latlonEstDeg(1:2))) && all(isfinite(truthLatlonDeg(1:2)))
  row.angleErrDeg = calcLatlonAngleError(latlonEstDeg(1:2), truthLatlonDeg(1:2));
end
row.initAngleErrDeg = localResolveInitAngleErrDeg(estResult, truth);
row.finalMinusInitAngleDeg = localResolveFinalMinusInitAngleDeg(estResult, row.angleErrDeg);
row.exitflag = localGetFieldOrDefault(info, 'exitflag', NaN);
row.iterations = localGetFieldOrDefault(info, 'iterations', NaN);
row.funcCount = localGetFieldOrDefault(info, 'funcCount', NaN);
row.initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objectiveImprove = row.initObj - row.finalObj;
row.firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
row.stepSize = localGetFieldOrDefault(optimInfo, 'stepsize', NaN);
row.fdRefEstHz = localGetFieldOrDefault(info, 'fdRefEstHz', NaN);
row.fdRefErrHz = localGetFieldOrDefault(info, 'fdRefErrHz', NaN);
row.fdRefAbsErrHz = abs(row.fdRefErrHz);
row.fdRateEstHzPerSec = localGetFieldOrDefault(info, 'fdRateEstHzPerSec', NaN);
row.fdRateErrHzPerSec = localGetFieldOrDefault(info, 'fdRateErrHzPerSec', NaN);
row.fdRateAbsErrHzPerSec = abs(row.fdRateErrHzPerSec);
row.initFdRefHz = localResolveInitFdRefHz(estResult);
row.initFdRateHzPerSec = localResolveInitFdRateHzPerSec(estResult);
row.fdRefInitMoveHz = row.fdRefEstHz - row.initFdRefHz;
row.fdRateInitMoveHzPerSec = row.fdRateEstHzPerSec - row.initFdRateHzPerSec;
row.fdLineRmseHz = localGetFieldOrDefault(info, 'fdLineRmseHz', NaN);
row.fdSatRmseHz = localGetFieldOrDefault(info, 'fdSatRmseHz', NaN);
row.deltaFdRmseHz = localGetFieldOrDefault(info, 'deltaFdRmseHz', NaN);
row.toothIdx = localGetFieldOrDefault(dynSummary, 'toothIdx', NaN);
row.toothResidualHz = localGetFieldOrDefault(dynSummary, 'toothResidualHz', NaN);
row.nonRefCoherenceFloor = localGetFieldOrDefault(dynSummary, 'nonRefCoherenceFloor', NaN);
row.nonRefRmsPhaseResidRad = localGetFieldOrDefault(dynSummary, 'nonRefRmsPhaseResidRad', NaN);
row.solveVariant = string(localGetFieldOrDefault(dynSummary, 'solveVariant', ""));
row.runTimeMs = localGetFieldOrDefault(info, 'runTimeMs', NaN);
row.wallTimeMs = wallTimeMs;
crbRow = localFindCrbRow(crbTable, method.displayName);
row.angleCrbStdDeg = crbRow.angleCrbStdDeg;
row.fdRefCrbStdHz = crbRow.fdRefCrbStdHz;
row.angleNormErr = localSafeRatio(row.angleErrDeg, row.angleCrbStdDeg);
row.fdRefNormErr = localSafeRatio(row.fdRefAbsErrHz, row.fdRefCrbStdHz);
inToothMeta = localGetFieldOrDefault(caseUse, 'inToothMeta', struct());
fdRangeUse = localGetFieldOrDefault(inToothMeta, 'fdRangeUse', [NaN, NaN]);
fdRateRangeUse = localGetFieldOrDefault(inToothMeta, 'fdRateRangeUse', [NaN, NaN]);
row.freqBoundaryHit = localIsNearRangeBoundary(row.fdRefEstHz, fdRangeUse, scanConfig.boundaryTolFraction);
if method.isKnownRate
  row.fdRateResolved = true;
  row.fdRateBoundaryHit = false;
else
  row.fdRateResolved = isfinite(row.fdRateAbsErrHzPerSec) && ...
    row.fdRateAbsErrHzPerSec <= scanConfig.resolvedFdRateAbsTolHzPerSec;
  row.fdRateBoundaryHit = localIsNearRangeBoundary(row.fdRateEstHzPerSec, fdRateRangeUse, scanConfig.boundaryTolFraction);
end
row.toothResolved = localIsToothResolved(row, scanConfig.resolvedToothHalfWidthFraction);
row.coherenceCollapsed = isfinite(row.nonRefCoherenceFloor) && ...
  row.nonRefCoherenceFloor < scanConfig.coherenceCollapseFloor;
row.isResolved = row.solverResolved && row.toothResolved && row.fdRateResolved && ...
  ~(row.freqBoundaryHit || row.fdRateBoundaryHit);
row.coreResolved = row.isResolved && localIsCoreCoherenceResolved(row, scanConfig.coreCoherenceFloor);
row.trimmedCore = localIsTrimmedCore(row, scanConfig.trimEnable, scanConfig.trimNormCap);
row.trimmedOutlier = row.coreResolved && ~row.trimmedCore;
row.failureReason = localBuildFailureReason(row);
row.tailSubtype = localBuildTailSubtype(row, scanConfig.trimNormCap);
end

function [latlonEstDeg, truthLatlonDeg] = localResolveCaseLatlon(caseUse, truth)
%LOCALRESOLVECASELATLON Resolve estimated and true lat/lon for angle metrics.

latlonEstDeg = [NaN; NaN];
truthLatlonDeg = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', [NaN; NaN]), [], 1);
if numel(truthLatlonDeg) < 2
  truthLatlonDeg = [NaN; NaN];
else
  truthLatlonDeg = truthLatlonDeg(1:2);
end
if ~(isstruct(caseUse) && isfield(caseUse, 'estResult') && ~isempty(caseUse.estResult))
  return;
end
try
  latlonUse = getDoaDopplerLatlonEst(caseUse.estResult);
  latlonUse = reshape(double(latlonUse), [], 1);
  if numel(latlonUse) >= 2
    latlonEstDeg = latlonUse(1:2);
  end
catch
  latlonEstDeg = [NaN; NaN];
end
end

function tf = localIsCoreCoherenceResolved(row, coreCoherenceFloor)
%LOCALISCORECOHERENCERESOLVED Apply the offline core-resolved coherence condition.

tf = true;
if row.satMode ~= "multi"
  return;
end
tf = isfinite(row.nonRefCoherenceFloor) && row.nonRefCoherenceFloor >= coreCoherenceFloor;
end

function tf = localIsTrimmedCore(row, trimEnable, trimNormCap)
%LOCALISTRIMMEDCORE Apply a fixed CRB-normalized tail trimming rule.

if ~row.coreResolved
  tf = false;
  return;
end
if ~trimEnable
  tf = true;
  return;
end
tf = isfinite(row.angleNormErr) && isfinite(row.fdRefNormErr) && ...
  abs(row.angleNormErr) <= trimNormCap && abs(row.fdRefNormErr) <= trimNormCap;
end

function crbRow = localFindCrbRow(crbTable, displayName)
%LOCALFINDCRBROW Read the CRB standard deviations for one display name.

crbRow = struct('angleCrbStdDeg', NaN, 'fdRefCrbStdHz', NaN);
if isempty(crbTable) || ~istable(crbTable)
  return;
end
mask = string(crbTable.displayName) == string(displayName);
idx = find(mask, 1, 'first');
if isempty(idx)
  return;
end
crbRow.angleCrbStdDeg = crbTable.angleCrbStdDeg(idx);
crbRow.fdRefCrbStdHz = crbTable.fdRefCrbStdHz(idx);
end

function tf = localIsToothResolved(row, resolvedToothHalfWidthFraction)
%LOCALISTOOTHRESOLVED Check the offline tooth-resolved condition.

toothTolHz = resolvedToothHalfWidthFraction * row.toothStepHz;
if isfinite(row.toothIdx) && isfinite(row.toothResidualHz)
  tf = abs(row.toothIdx) == 0 && abs(row.toothResidualHz) <= toothTolHz;
else
  tf = isfinite(row.fdRefAbsErrHz) && row.fdRefAbsErrHz <= toothTolHz;
end
end

function tf = localIsNearRangeBoundary(value, rangeValue, boundaryTolFraction)
%LOCALISNEARRANGEBOUNDARY Check if one value is too close to a frequency bound.

tf = false;
rangeValue = reshape(double(rangeValue), 1, []);
if numel(rangeValue) ~= 2 || any(~isfinite(rangeValue)) || ~isfinite(value)
  return;
end
rangeWidth = abs(rangeValue(2) - rangeValue(1));
if rangeWidth <= 0
  return;
end
boundaryTol = max(boundaryTolFraction * rangeWidth, eps(rangeWidth));
tf = abs(value - min(rangeValue)) <= boundaryTol || abs(value - max(rangeValue)) <= boundaryTol;
end

function reason = localBuildFailureReason(row)
%LOCALBUILDFAILUREREASON Assign an offline resolved/outlier reason label.

if row.isResolved
  if row.coherenceCollapsed
    reason = "resolved-coherence-low";
  else
    reason = "resolved";
  end
elseif ~row.solverResolved
  reason = "solver-invalid";
elseif ~row.toothResolved
  reason = "outside-resolved-tooth";
elseif row.freqBoundaryHit || row.fdRateBoundaryHit
  reason = "frequency-boundary-hit";
elseif ~row.fdRateResolved
  reason = "fdRate-unresolved";
else
  reason = "unclassified-outlier";
end
end

function subtype = localBuildTailSubtype(row, normCap)
%LOCALBUILDTAILSUBTYPE Classify the dominant outlier mechanism for diagnostics.

subtype = "crb-local";
if row.trimmedCore
  return;
end
if ~row.solverResolved
  subtype = "solver-tail";
elseif ~row.toothResolved
  subtype = "outside-tooth-tail";
elseif row.freqBoundaryHit || row.fdRateBoundaryHit
  subtype = "fd-boundary-tail";
elseif ~row.fdRateResolved
  subtype = "fdrate-unresolved-tail";
elseif row.coreResolved && isfinite(row.fdRefNormErr) && abs(row.fdRefNormErr) > normCap
  subtype = "same-tooth-fdref-tail";
elseif row.coreResolved && isfinite(row.angleNormErr) && abs(row.angleNormErr) > normCap
  subtype = "angle-local-tail";
elseif row.coherenceCollapsed
  subtype = "coherence-collapse-tail";
else
  subtype = "other-tail";
end
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for per-repeat case metrics.

row = struct('displayName', "", 'satMode', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'taskSeed', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, ...
  'windowSec', NaN, 'toothStepHz', NaN, 'inToothMode', "", ...
  'fdRateRangeMode', "", 'doaSeedMode', "", 'initMode', "", ...
  'latEstDeg', NaN, 'lonEstDeg', NaN, 'latTrueDeg', NaN, 'lonTrueDeg', NaN, ...
  'angleErrDeg', NaN, 'initAngleErrDeg', NaN, 'finalMinusInitAngleDeg', NaN, ...
  'exitflag', NaN, 'iterations', NaN, 'funcCount', NaN, 'initObj', NaN, ...
  'finalObj', NaN, 'objectiveImprove', NaN, 'firstOrderOpt', NaN, 'stepSize', NaN, ...
  'fdRefEstHz', NaN, 'fdRefErrHz', NaN, ...
  'fdRefAbsErrHz', NaN, 'fdRateEstHzPerSec', NaN, 'fdRateErrHzPerSec', NaN, ...
  'fdRateAbsErrHzPerSec', NaN, 'initFdRefHz', NaN, 'initFdRateHzPerSec', NaN, ...
  'fdRefInitMoveHz', NaN, 'fdRateInitMoveHzPerSec', NaN, ...
  'fdLineRmseHz', NaN, 'fdSatRmseHz', NaN, ...
  'deltaFdRmseHz', NaN, 'toothIdx', NaN, 'toothResidualHz', NaN, ...
  'nonRefCoherenceFloor', NaN, 'nonRefRmsPhaseResidRad', NaN, ...
  'angleCrbStdDeg', NaN, 'fdRefCrbStdHz', NaN, 'angleNormErr', NaN, ...
  'fdRefNormErr', NaN, 'solverResolved', false, 'toothResolved', false, ...
  'fdRateResolved', false, 'freqBoundaryHit', false, 'fdRateBoundaryHit', false, ...
  'coherenceCollapsed', false, 'isResolved', false, 'coreResolved', false, ...
  'trimmedCore', false, 'trimmedOutlier', false, 'failureReason', "", ...
  'tailSubtype', "", 'solveVariant', "", 'runTimeMs', NaN, 'wallTimeMs', NaN);
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

function aggregateTable = localBuildAggregateTable(perfTable, scanConfig)
%LOCALBUILDAGGREGATETABLE Aggregate full, resolved, core, and trimmed metrics.

if isempty(perfTable) || height(perfTable) == 0
  aggregateTable = struct2table(repmat(localEmptyAggregateRow(), 0, 1));
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'satMode', 'phaseMode', 'fdRateMode', ...
  'snrDb', 'numFrame', 'frameIntvlSec'}), 'rows', 'stable');
rowList = repmat(localEmptyAggregateRow(), height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iRow) & ...
    perfTable.snrDb == groupKey.snrDb(iRow) & ...
    perfTable.numFrame == groupKey.numFrame(iRow) & ...
    abs(perfTable.frameIntvlSec - groupKey.frameIntvlSec(iRow)) < eps;
  resolvedMask = mask & perfTable.isResolved;
  coreMask = mask & perfTable.coreResolved;
  trimMask = mask & perfTable.trimmedCore;
  crbMask = mask;
  row = localEmptyAggregateRow();
  row.displayName = groupKey.displayName(iRow);
  row.satMode = groupKey.satMode(iRow);
  row.phaseMode = groupKey.phaseMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.numFrame = groupKey.numFrame(iRow);
  row.frameIntvlSec = groupKey.frameIntvlSec(iRow);
  row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
  row.numRepeat = nnz(mask);
  row.resolvedCount = nnz(resolvedMask);
  row.resolvedRate = localSafeRatio(row.resolvedCount, row.numRepeat);
  row.outlierRate = 1 - row.resolvedRate;
  row.coreResolvedCount = nnz(coreMask);
  row.coreResolvedRate = localSafeRatio(row.coreResolvedCount, row.numRepeat);
  row.coreRejectedRate = row.resolvedRate - row.coreResolvedRate;
  row.trimmedCoreCount = nnz(trimMask);
  row.trimmedCoreRate = localSafeRatio(row.trimmedCoreCount, row.numRepeat);
  row.trimKeepRate = localSafeRatio(row.trimmedCoreCount, row.coreResolvedCount);
  row.trimmedOutlierRate = localSafeRatio(nnz(mask & perfTable.trimmedOutlier), row.numRepeat);
  row.trimKeepRateBelowMin = isfinite(row.trimKeepRate) && ...
    row.trimKeepRate < scanConfig.trimMinKeepRate;
  row.toothResolvedRate = mean(double(perfTable.toothResolved(mask)), 'omitnan');
  row.freqBoundaryHitRate = mean(double(perfTable.freqBoundaryHit(mask) | perfTable.fdRateBoundaryHit(mask)), 'omitnan');
  row.fdRateResolvedRate = mean(double(perfTable.fdRateResolved(mask)), 'omitnan');
  row.coherenceCollapseRate = mean(double(perfTable.coherenceCollapsed(mask)), 'omitnan');
  row.fullAngleMseDeg2 = localMse(perfTable.angleErrDeg(mask));
  row.fullAngleRmseDeg = sqrt(row.fullAngleMseDeg2);
  row.fullAngleP95Deg = localPrctileFinite(perfTable.angleErrDeg(mask), 95);
  row.fullAngleP99Deg = localPrctileFinite(perfTable.angleErrDeg(mask), 99);
  row.resolvedAngleMseDeg2 = localMse(perfTable.angleErrDeg(resolvedMask));
  row.resolvedAngleRmseDeg = sqrt(row.resolvedAngleMseDeg2);
  row.resolvedAngleP95Deg = localPrctileFinite(perfTable.angleErrDeg(resolvedMask), 95);
  row.resolvedAngleP99Deg = localPrctileFinite(perfTable.angleErrDeg(resolvedMask), 99);
  row.coreAngleMseDeg2 = localMse(perfTable.angleErrDeg(coreMask));
  row.coreAngleRmseDeg = sqrt(row.coreAngleMseDeg2);
  row.coreAngleP95Deg = localPrctileFinite(perfTable.angleErrDeg(coreMask), 95);
  row.coreAngleP99Deg = localPrctileFinite(perfTable.angleErrDeg(coreMask), 99);
  row.trimAngleMseDeg2 = localMse(perfTable.angleErrDeg(trimMask));
  row.trimAngleRmseDeg = sqrt(row.trimAngleMseDeg2);
  row.trimAngleP95Deg = localPrctileFinite(perfTable.angleErrDeg(trimMask), 95);
  row.trimAngleP99Deg = localPrctileFinite(perfTable.angleErrDeg(trimMask), 99);
  row.crbLocalCount = row.trimmedCoreCount;
  row.crbLocalRate = row.trimmedCoreRate;
  row.crbLocalKeepRate = row.trimKeepRate;
  row.crbLocalAngleMseDeg2 = row.trimAngleMseDeg2;
  row.crbLocalAngleRmseDeg = row.trimAngleRmseDeg;
  row.crbLocalAngleP95Deg = row.trimAngleP95Deg;
  row.crbLocalAngleP99Deg = row.trimAngleP99Deg;
  row.angleCrbStdDeg = median(perfTable.angleCrbStdDeg(crbMask), 'omitnan');
  row.fullAngleRmseOverCrb = localSafeRatio(row.fullAngleRmseDeg, row.angleCrbStdDeg);
  row.resolvedAngleRmseOverCrb = localSafeRatio(row.resolvedAngleRmseDeg, row.angleCrbStdDeg);
  row.coreAngleRmseOverCrb = localSafeRatio(row.coreAngleRmseDeg, row.angleCrbStdDeg);
  row.trimAngleRmseOverCrb = localSafeRatio(row.trimAngleRmseDeg, row.angleCrbStdDeg);
  row.crbLocalAngleRmseOverCrb = row.trimAngleRmseOverCrb;
  row.fullAngleMseOverCrb = localSafeRatio(row.fullAngleMseDeg2, row.angleCrbStdDeg.^2);
  row.resolvedAngleMseOverCrb = localSafeRatio(row.resolvedAngleMseDeg2, row.angleCrbStdDeg.^2);
  row.coreAngleMseOverCrb = localSafeRatio(row.coreAngleMseDeg2, row.angleCrbStdDeg.^2);
  row.trimAngleMseOverCrb = localSafeRatio(row.trimAngleMseDeg2, row.angleCrbStdDeg.^2);
  row.crbLocalAngleMseOverCrb = row.trimAngleMseOverCrb;
  row.resolvedAngleNormP95 = localPrctileFinite(perfTable.angleNormErr(resolvedMask), 95);
  row.coreAngleNormP95 = localPrctileFinite(perfTable.angleNormErr(coreMask), 95);
  row.trimAngleNormP95 = localPrctileFinite(perfTable.angleNormErr(trimMask), 95);
  row.crbLocalAngleNormP95 = row.trimAngleNormP95;
  row.fullFdRefMseHz2 = localMse(perfTable.fdRefAbsErrHz(mask));
  row.fullFdRefRmseHz = sqrt(row.fullFdRefMseHz2);
  row.fullFdRefP95Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(mask), 95);
  row.fullFdRefP99Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(mask), 99);
  row.resolvedFdRefMseHz2 = localMse(perfTable.fdRefAbsErrHz(resolvedMask));
  row.resolvedFdRefRmseHz = sqrt(row.resolvedFdRefMseHz2);
  row.resolvedFdRefP95Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(resolvedMask), 95);
  row.resolvedFdRefP99Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(resolvedMask), 99);
  row.coreFdRefMseHz2 = localMse(perfTable.fdRefAbsErrHz(coreMask));
  row.coreFdRefRmseHz = sqrt(row.coreFdRefMseHz2);
  row.coreFdRefP95Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(coreMask), 95);
  row.coreFdRefP99Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(coreMask), 99);
  row.trimFdRefMseHz2 = localMse(perfTable.fdRefAbsErrHz(trimMask));
  row.trimFdRefRmseHz = sqrt(row.trimFdRefMseHz2);
  row.trimFdRefP95Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(trimMask), 95);
  row.trimFdRefP99Hz = localPrctileFinite(perfTable.fdRefAbsErrHz(trimMask), 99);
  row.crbLocalFdRefMseHz2 = row.trimFdRefMseHz2;
  row.crbLocalFdRefRmseHz = row.trimFdRefRmseHz;
  row.crbLocalFdRefP95Hz = row.trimFdRefP95Hz;
  row.crbLocalFdRefP99Hz = row.trimFdRefP99Hz;
  row.fdRefCrbStdHz = median(perfTable.fdRefCrbStdHz(crbMask), 'omitnan');
  row.fullFdRefRmseOverCrb = localSafeRatio(row.fullFdRefRmseHz, row.fdRefCrbStdHz);
  row.resolvedFdRefRmseOverCrb = localSafeRatio(row.resolvedFdRefRmseHz, row.fdRefCrbStdHz);
  row.coreFdRefRmseOverCrb = localSafeRatio(row.coreFdRefRmseHz, row.fdRefCrbStdHz);
  row.trimFdRefRmseOverCrb = localSafeRatio(row.trimFdRefRmseHz, row.fdRefCrbStdHz);
  row.crbLocalFdRefRmseOverCrb = row.trimFdRefRmseOverCrb;
  row.fullFdRefMseOverCrb = localSafeRatio(row.fullFdRefMseHz2, row.fdRefCrbStdHz.^2);
  row.resolvedFdRefMseOverCrb = localSafeRatio(row.resolvedFdRefMseHz2, row.fdRefCrbStdHz.^2);
  row.coreFdRefMseOverCrb = localSafeRatio(row.coreFdRefMseHz2, row.fdRefCrbStdHz.^2);
  row.trimFdRefMseOverCrb = localSafeRatio(row.trimFdRefMseHz2, row.fdRefCrbStdHz.^2);
  row.crbLocalFdRefMseOverCrb = row.trimFdRefMseOverCrb;
  row.resolvedFdRefNormP95 = localPrctileFinite(perfTable.fdRefNormErr(resolvedMask), 95);
  row.coreFdRefNormP95 = localPrctileFinite(perfTable.fdRefNormErr(coreMask), 95);
  row.trimFdRefNormP95 = localPrctileFinite(perfTable.fdRefNormErr(trimMask), 95);
  row.crbLocalFdRefNormP95 = row.trimFdRefNormP95;
  row.fdRateRmseHzPerSec = localRmse(perfTable.fdRateAbsErrHzPerSec(mask));
  row.nonRefCoherenceFloorMedian = median(perfTable.nonRefCoherenceFloor(mask), 'omitnan');
  row.meanWallTimeMs = mean(perfTable.wallTimeMs(mask), 'omitnan');
  rowList(iRow) = row;
end
aggregateTable = struct2table(rowList);
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return a typed empty row for aggregate metrics.

row = struct('displayName', "", 'satMode', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, 'windowSec', NaN, ...
  'numRepeat', NaN, 'resolvedCount', NaN, 'resolvedRate', NaN, 'outlierRate', NaN, ...
  'coreResolvedCount', NaN, 'coreResolvedRate', NaN, 'coreRejectedRate', NaN, ...
  'trimmedCoreCount', NaN, 'trimmedCoreRate', NaN, 'trimKeepRate', NaN, ...
  'trimmedOutlierRate', NaN, 'trimKeepRateBelowMin', false, ...
  'crbLocalCount', NaN, 'crbLocalRate', NaN, 'crbLocalKeepRate', NaN, ...
  'toothResolvedRate', NaN, 'freqBoundaryHitRate', NaN, 'fdRateResolvedRate', NaN, ...
  'coherenceCollapseRate', NaN, ...
  'fullAngleMseDeg2', NaN, 'fullAngleRmseDeg', NaN, 'fullAngleP95Deg', NaN, ...
  'fullAngleP99Deg', NaN, 'resolvedAngleMseDeg2', NaN, 'resolvedAngleRmseDeg', NaN, ...
  'resolvedAngleP95Deg', NaN, 'resolvedAngleP99Deg', NaN, ...
  'coreAngleMseDeg2', NaN, 'coreAngleRmseDeg', NaN, 'coreAngleP95Deg', NaN, ...
  'coreAngleP99Deg', NaN, 'trimAngleMseDeg2', NaN, 'trimAngleRmseDeg', NaN, ...
  'trimAngleP95Deg', NaN, 'trimAngleP99Deg', NaN, ...
  'crbLocalAngleMseDeg2', NaN, 'crbLocalAngleRmseDeg', NaN, ...
  'crbLocalAngleP95Deg', NaN, 'crbLocalAngleP99Deg', NaN, 'angleCrbStdDeg', NaN, ...
  'fullAngleRmseOverCrb', NaN, 'resolvedAngleRmseOverCrb', NaN, ...
  'coreAngleRmseOverCrb', NaN, 'trimAngleRmseOverCrb', NaN, ...
  'crbLocalAngleRmseOverCrb', NaN, ...
  'fullAngleMseOverCrb', NaN, 'resolvedAngleMseOverCrb', NaN, ...
  'coreAngleMseOverCrb', NaN, 'trimAngleMseOverCrb', NaN, ...
  'crbLocalAngleMseOverCrb', NaN, ...
  'resolvedAngleNormP95', NaN, 'coreAngleNormP95', NaN, 'trimAngleNormP95', NaN, ...
  'crbLocalAngleNormP95', NaN, ...
  'fullFdRefMseHz2', NaN, 'fullFdRefRmseHz', NaN, 'fullFdRefP95Hz', NaN, ...
  'fullFdRefP99Hz', NaN, 'resolvedFdRefMseHz2', NaN, 'resolvedFdRefRmseHz', NaN, ...
  'resolvedFdRefP95Hz', NaN, 'resolvedFdRefP99Hz', NaN, ...
  'coreFdRefMseHz2', NaN, 'coreFdRefRmseHz', NaN, 'coreFdRefP95Hz', NaN, ...
  'coreFdRefP99Hz', NaN, 'trimFdRefMseHz2', NaN, 'trimFdRefRmseHz', NaN, ...
  'trimFdRefP95Hz', NaN, 'trimFdRefP99Hz', NaN, ...
  'crbLocalFdRefMseHz2', NaN, 'crbLocalFdRefRmseHz', NaN, ...
  'crbLocalFdRefP95Hz', NaN, 'crbLocalFdRefP99Hz', NaN, 'fdRefCrbStdHz', NaN, ...
  'fullFdRefRmseOverCrb', NaN, 'resolvedFdRefRmseOverCrb', NaN, ...
  'coreFdRefRmseOverCrb', NaN, 'trimFdRefRmseOverCrb', NaN, ...
  'crbLocalFdRefRmseOverCrb', NaN, ...
  'fullFdRefMseOverCrb', NaN, 'resolvedFdRefMseOverCrb', NaN, ...
  'coreFdRefMseOverCrb', NaN, 'trimFdRefMseOverCrb', NaN, ...
  'crbLocalFdRefMseOverCrb', NaN, ...
  'resolvedFdRefNormP95', NaN, 'coreFdRefNormP95', NaN, 'trimFdRefNormP95', NaN, ...
  'crbLocalFdRefNormP95', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'nonRefCoherenceFloorMedian', NaN, 'meanWallTimeMs', NaN);
end

function topTailTable = localBuildTopTailTable(perfTable, topCountPerGroup)
%LOCALBUILDTOPTAILTABLE Keep the worst CRB-normalized tails for diagnostics.

if isempty(perfTable) || height(perfTable) == 0 || topCountPerGroup <= 0
  topTailTable = table();
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'snrDb', 'numFrame'}), 'rows', 'stable');
tailTableCell = {};
for iGroup = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iGroup) & ...
    perfTable.snrDb == groupKey.snrDb(iGroup) & perfTable.numFrame == groupKey.numFrame(iGroup);
  groupTable = perfTable(mask, :);
  tailScoreMat = [abs(groupTable.angleNormErr), abs(groupTable.fdRefNormErr)];
  tailScoreMat(~isfinite(tailScoreMat)) = -Inf;
  tailScore = max(tailScoreMat, [], 2);
  [~, orderIdx] = sort(tailScore, 'descend');
  keepCount = min(topCountPerGroup, numel(orderIdx));
  if keepCount <= 0
    continue;
  end
  keepIdx = orderIdx(1:keepCount);
  topTable = groupTable(keepIdx, {'displayName', 'snrDb', 'numFrame', 'taskSeed', ...
    'angleErrDeg', 'angleNormErr', 'fdRefAbsErrHz', 'fdRefNormErr', ...
    'fdRateAbsErrHzPerSec', 'nonRefCoherenceFloor', 'initAngleErrDeg', ...
    'finalMinusInitAngleDeg', 'fdRefInitMoveHz', 'fdRateInitMoveHzPerSec', ...
    'objectiveImprove', 'iterations', 'exitflag', 'firstOrderOpt', ...
    'failureReason', 'tailSubtype', 'coreResolved', 'trimmedCore', 'trimmedOutlier'});
  topTable.tailScore = tailScore(keepIdx);
  tailTableCell{end + 1, 1} = topTable; %#ok<AGROW>
end
if isempty(tailTableCell)
  topTailTable = table();
  return;
end
topTailTable = vertcat(tailTableCell{:});
topTailTable = movevars(topTailTable, 'tailScore', 'After', 'numFrame');
end

function previewTable = localTablePreview(dataTable, maxRows)
%LOCALTABLEPREVIEW Return a compact top slice of a table for command output.

previewTable = dataTable;
if isempty(dataTable) || height(dataTable) <= maxRows
  return;
end
previewTable = dataTable(1:maxRows, :);
end

function topTailExportTable = localBuildTopTailExportTable(topTailTable)
%LOCALBUILDTOPTAILEXPORTTABLE Keep compact tail rows for replay seed selection.

topTailExportTable = table();
if isempty(topTailTable) || height(topTailTable) == 0
  return;
end
fieldList = {'displayName', 'snrDb', 'numFrame', 'taskSeed', 'tailScore', ...
  'angleErrDeg', 'angleNormErr', 'fdRefAbsErrHz', 'fdRefNormErr', ...
  'fdRateAbsErrHzPerSec', 'nonRefCoherenceFloor', 'initAngleErrDeg', ...
  'finalMinusInitAngleDeg', 'fdRefInitMoveHz', 'fdRateInitMoveHzPerSec', ...
  'objectiveImprove', 'iterations', 'exitflag', 'firstOrderOpt', ...
  'failureReason', 'tailSubtype', 'coreResolved', 'trimmedCore', 'trimmedOutlier'};
keepField = fieldList(ismember(fieldList, topTailTable.Properties.VariableNames));
topTailExportTable = topTailTable(:, keepField);
end

function initAngleErrDeg = localResolveInitAngleErrDeg(estResult, truth)
%LOCALRESOLVEINITANGLEERRDEG Read estimator init DoA and compare to truth.

initAngleErrDeg = NaN;
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
truthDoa = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
if numel(initParam) >= 2 && numel(truthDoa) >= 2 && all(isfinite(initParam(1:2)))
  initAngleErrDeg = calcLatlonAngleError(initParam(1:2), truthDoa(1:2));
end
end

function finalMinusInitDeg = localResolveFinalMinusInitAngleDeg(estResult, fallbackAngleErrDeg)
%LOCALRESOLVEFINALMINUSINITANGLEDEG Measure final DoA movement from estimator init.

finalMinusInitDeg = NaN;
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
finalDoa = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(initParam) >= 2 && numel(finalDoa) >= 2 && all(isfinite(initParam(1:2))) && all(isfinite(finalDoa(1:2)))
  finalMinusInitDeg = calcLatlonAngleError(finalDoa(1:2), initParam(1:2));
elseif isfinite(fallbackAngleErrDeg)
  finalMinusInitDeg = NaN;
end
end

function initFdRefHz = localResolveInitFdRefHz(estResult)
%LOCALRESOLVEINITFDREFHZ Read fdRef from estimator initParam when available.

initFdRefHz = NaN;
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
if numel(initParam) >= 3
  initFdRefHz = initParam(3);
end
end

function initFdRateHzPerSec = localResolveInitFdRateHzPerSec(estResult)
%LOCALRESOLVEINITFDRATEHZPERSEC Read fdRate from estimator initParam when available.

initFdRateHzPerSec = NaN;
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
if numel(initParam) >= 4
  initFdRateHzPerSec = initParam(4);
end
end

function scanData = localValidateScanDataForSummary(scanData)
%LOCALVALIDATESCANDATAFORSUMMARY Rebuild summary-only fields from saved scanData.

if ~isfield(scanData, 'aggregateTable') || ~istable(scanData.aggregateTable)
  error('scanMfMsMleCrbInToothConsistency:InvalidScanData', ...
    'scanData.aggregateTable is missing. Load a scan snapshot that stores summary inputs.');
end
if ~isfield(scanData, 'scanName') || strlength(string(scanData.scanName)) == 0
  scanData.scanName = "scanMfMsMleCrbInToothConsistency";
end
if ~isfield(scanData, 'config') || ~isstruct(scanData.config)
  scanData.config = struct();
end
scanData.config = localCompleteSummaryConfig(scanData.config, scanData.aggregateTable);
if ~isfield(scanData, 'snapshotFile')
  scanData.snapshotFile = "";
end
if ~isfield(scanData, 'checkpointDir')
  scanData.checkpointDir = "";
end
if ~isfield(scanData, 'elapsedSec')
  scanData.elapsedSec = NaN;
end

scanData.crbLocalSummaryTable = localBuildCrbLocalSummaryTable(scanData.aggregateTable);
if ~isfield(scanData, 'perfTable') || ~istable(scanData.perfTable)
  scanData.perfTable = table();
end
if ~isfield(scanData, 'failureSummaryTable') || ~istable(scanData.failureSummaryTable)
  scanData.failureSummaryTable = localBuildFailureSummaryTable(scanData.perfTable);
end
if ~isfield(scanData, 'msErrorTypeSummaryTable') || ~istable(scanData.msErrorTypeSummaryTable)
  scanData.msErrorTypeSummaryTable = localBuildMsErrorTypeSummaryTable(scanData.perfTable, scanData.config);
end
if ~isfield(scanData, 'representativeSeedTable') || ~istable(scanData.representativeSeedTable)
  scanData.representativeSeedTable = localBuildRepresentativeSeedTable(scanData.perfTable, scanData.config);
end
if ~isfield(scanData, 'topTailTable') || ~istable(scanData.topTailTable)
  scanData.topTailTable = localBuildTopTailTable(scanData.perfTable, scanData.config.topTailCountPerGroup);
end
if ~isfield(scanData, 'topTailExportTable') || ~istable(scanData.topTailExportTable)
  scanData.topTailExportTable = localBuildTopTailExportTable(scanData.topTailTable);
end
if ~isfield(scanData, 'checkpointSummaryTable') || ~istable(scanData.checkpointSummaryTable)
  scanData.checkpointSummaryTable = table();
end
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  scanData.plotData = localBuildPlotData(scanData.aggregateTable, ...
    scanData.failureSummaryTable, scanData.topTailTable);
end
end

function scanConfig = localCompleteSummaryConfig(scanConfig, aggregateTable)
%LOCALCOMPLETESUMMARYCONFIG Fill reporting defaults needed by snapshot summaries.

if ~isstruct(scanConfig)
  scanConfig = struct();
end
scanConfig = localSetMissingField(scanConfig, 'frameCountList', localUniqueNumericColumn(aggregateTable, 'numFrame'));
scanConfig = localSetMissingField(scanConfig, 'snrDbList', localUniqueNumericColumn(aggregateTable, 'snrDb'));
scanConfig = localSetMissingField(scanConfig, 'methodNameList', localUniqueStringColumn(aggregateTable, 'displayName'));
scanConfig = localSetMissingField(scanConfig, 'initMode', "");
scanConfig = localSetMissingField(scanConfig, 'diagnosticVersion', "");
scanConfig = localSetMissingField(scanConfig, 'representativeSeedCountPerType', 3);
scanConfig = localSetMissingField(scanConfig, 'numRepeat', localResolveRepeatCountForSummary(aggregateTable));
scanConfig = localSetMissingField(scanConfig, 'oracleFdHalfToothFraction', NaN);
scanConfig = localSetMissingField(scanConfig, 'oracleFdRateHalfWidthHzPerSec', NaN);
scanConfig = localSetMissingField(scanConfig, 'resolvedToothHalfWidthFraction', NaN);
scanConfig = localSetMissingField(scanConfig, 'resolvedFdRateAbsTolHzPerSec', NaN);
scanConfig = localSetMissingField(scanConfig, 'coreCoherenceFloor', 0.8);
scanConfig = localSetMissingField(scanConfig, 'trimNormCap', 5);
scanConfig = localSetMissingField(scanConfig, 'trimMinKeepRate', NaN);
scanConfig = localSetMissingField(scanConfig, 'topTailCountPerGroup', 5);
scanConfig = localSetMissingField(scanConfig, 'topTailPrintMaxRows', 60);
end

function dataStruct = localSetMissingField(dataStruct, fieldName, defaultValue)
%LOCALSETMISSINGFIELD Set one default only when a field is absent or empty.

if ~isfield(dataStruct, fieldName) || isempty(dataStruct.(fieldName))
  dataStruct.(fieldName) = defaultValue;
end
end

function value = localUniqueNumericColumn(dataTable, fieldName)
%LOCALUNIQUENUMERICCOLUMN Return a stable numeric column list when present.

value = [];
if istable(dataTable) && ismember(fieldName, dataTable.Properties.VariableNames)
  value = unique(double(dataTable.(fieldName)), 'stable');
  value = reshape(value(isfinite(value)), [], 1);
end
end

function value = localUniqueStringColumn(dataTable, fieldName)
%LOCALUNIQUESTRINGCOLUMN Return a stable string column list when present.

value = strings(0, 1);
if istable(dataTable) && ismember(fieldName, dataTable.Properties.VariableNames)
  value = unique(string(dataTable.(fieldName)), 'stable');
  value = reshape(value(strlength(value) > 0), [], 1);
end
end

function value = localResolveRepeatCountForSummary(aggregateTable)
%LOCALRESOLVEREPEATCOUNTFORSUMMARY Resolve repeat count for summary-only runs.

value = NaN;
if istable(aggregateTable) && ismember('numRepeat', aggregateTable.Properties.VariableNames) && ...
    ~isempty(aggregateTable.numRepeat)
  value = max(double(aggregateTable.numRepeat), [], 'omitnan');
end
end

function crbLocalSummaryTable = localBuildCrbLocalSummaryTable(aggregateTable)
%LOCALBUILDCRBLOCALSUMMARYTABLE Build a compact paper-facing local-regime table.

crbLocalSummaryTable = table();
if isempty(aggregateTable) || height(aggregateTable) == 0
  return;
end
fieldList = {'displayName', 'snrDb', 'numFrame', 'numRepeat', ...
  'crbLocalAngleRmseOverCrb', 'crbLocalFdRefRmseOverCrb', ...
  'crbLocalAngleMseOverCrb', 'crbLocalFdRefMseOverCrb', ...
  'coreAngleRmseOverCrb', 'coreFdRefRmseOverCrb', ...
  'coreAngleMseOverCrb', 'coreFdRefMseOverCrb', ...
  'resolvedAngleRmseOverCrb', 'resolvedFdRefRmseOverCrb', ...
  'fullAngleRmseOverCrb', 'fullFdRefRmseOverCrb', ...
  'crbLocalAngleNormP95', 'crbLocalFdRefNormP95', ...
  'resolvedRate', 'coreResolvedRate', 'crbLocalRate', 'crbLocalKeepRate', ...
  'outlierRate', 'freqBoundaryHitRate', 'fdRateResolvedRate'};
keepField = fieldList(ismember(fieldList, aggregateTable.Properties.VariableNames));
crbLocalSummaryTable = aggregateTable(:, keepField);
end

function failureSummaryTable = localBuildFailureSummaryTable(perfTable)
%LOCALBUILDFAILURESUMMARYTABLE Count offline resolved/outlier reasons.

if isempty(perfTable) || height(perfTable) == 0
  failureSummaryTable = table();
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'snrDb', 'numFrame', 'failureReason'}), 'rows', 'stable');
count = zeros(height(groupKey), 1);
rate = zeros(height(groupKey), 1);
for iRow = 1:height(groupKey)
  baseMask = perfTable.displayName == groupKey.displayName(iRow) & ...
    perfTable.snrDb == groupKey.snrDb(iRow) & perfTable.numFrame == groupKey.numFrame(iRow);
  reasonMask = baseMask & perfTable.failureReason == groupKey.failureReason(iRow);
  count(iRow) = nnz(reasonMask);
  rate(iRow) = localSafeRatio(count(iRow), nnz(baseMask));
end
failureSummaryTable = groupKey;
failureSummaryTable.count = count;
failureSummaryTable.rate = rate;
end


function msErrorTypeSummaryTable = localBuildMsErrorTypeSummaryTable(perfTable, scanConfig)
%LOCALBUILDMSERRORTYPESUMMARYTABLE Count MS offline error types by SNR.

msErrorTypeSummaryTable = table();
if isempty(perfTable) || height(perfTable) == 0
  return;
end
msErrorType = strings(height(perfTable), 1);
for iRow = 1:height(perfTable)
  msErrorType(iRow) = localClassifyMsErrorType(perfTable(iRow, :), scanConfig);
end
typeKey = table(perfTable.displayName, perfTable.snrDb, perfTable.numFrame, msErrorType, ...
  'VariableNames', {'displayName', 'snrDb', 'numFrame', 'msErrorType'});
groupKey = unique(typeKey, 'rows', 'stable');
count = zeros(height(groupKey), 1);
rate = NaN(height(groupKey), 1);
badOnlyRate = NaN(height(groupKey), 1);
medianAngleNormErr = NaN(height(groupKey), 1);
medianFdRefNormErr = NaN(height(groupKey), 1);
medianNonRefCoherenceFloor = NaN(height(groupKey), 1);
medianIterations = NaN(height(groupKey), 1);
medianFirstOrderOpt = NaN(height(groupKey), 1);
for iGroup = 1:height(groupKey)
  baseMask = perfTable.displayName == groupKey.displayName(iGroup) & ...
    perfTable.snrDb == groupKey.snrDb(iGroup) & perfTable.numFrame == groupKey.numFrame(iGroup);
  typeMask = baseMask & msErrorType == groupKey.msErrorType(iGroup);
  badMask = baseMask & msErrorType ~= "crb-local";
  count(iGroup) = nnz(typeMask);
  rate(iGroup) = localSafeRatio(count(iGroup), nnz(baseMask));
  if groupKey.msErrorType(iGroup) ~= "crb-local"
    badOnlyRate(iGroup) = localSafeRatio(count(iGroup), nnz(badMask));
  end
  medianAngleNormErr(iGroup) = median(abs(perfTable.angleNormErr(typeMask)), 'omitnan');
  medianFdRefNormErr(iGroup) = median(abs(perfTable.fdRefNormErr(typeMask)), 'omitnan');
  medianNonRefCoherenceFloor(iGroup) = median(perfTable.nonRefCoherenceFloor(typeMask), 'omitnan');
  medianIterations(iGroup) = median(perfTable.iterations(typeMask), 'omitnan');
  medianFirstOrderOpt(iGroup) = median(perfTable.firstOrderOpt(typeMask), 'omitnan');
end
msErrorTypeSummaryTable = groupKey;
msErrorTypeSummaryTable.count = count;
msErrorTypeSummaryTable.rate = rate;
msErrorTypeSummaryTable.badOnlyRate = badOnlyRate;
msErrorTypeSummaryTable.medianAngleNormErr = medianAngleNormErr;
msErrorTypeSummaryTable.medianFdRefNormErr = medianFdRefNormErr;
msErrorTypeSummaryTable.medianNonRefCoherenceFloor = medianNonRefCoherenceFloor;
msErrorTypeSummaryTable.medianIterations = medianIterations;
msErrorTypeSummaryTable.medianFirstOrderOpt = medianFirstOrderOpt;
end

function representativeSeedTable = localBuildRepresentativeSeedTable(perfTable, scanConfig)
%LOCALBUILDREPRESENTATIVESEEDTABLE Pick replay seeds for each MS error type.

representativeSeedTable = table();
if isempty(perfTable) || height(perfTable) == 0 || scanConfig.representativeSeedCountPerType <= 0
  return;
end
msErrorType = strings(height(perfTable), 1);
tailScore = NaN(height(perfTable), 1);
for iRow = 1:height(perfTable)
  msErrorType(iRow) = localClassifyMsErrorType(perfTable(iRow, :), scanConfig);
  tailScore(iRow) = localBuildMsTailScore(perfTable(iRow, :), msErrorType(iRow), scanConfig);
end
tailMask = msErrorType ~= "crb-local";
if ~any(tailMask)
  return;
end
typeKey = table(perfTable.displayName(tailMask), perfTable.snrDb(tailMask), ...
  perfTable.numFrame(tailMask), msErrorType(tailMask), ...
  'VariableNames', {'displayName', 'snrDb', 'numFrame', 'msErrorType'});
groupKey = unique(typeKey, 'rows', 'stable');
rowCell = cell(0, 1);
for iGroup = 1:height(groupKey)
  groupMask = perfTable.displayName == groupKey.displayName(iGroup) & ...
    perfTable.snrDb == groupKey.snrDb(iGroup) & perfTable.numFrame == groupKey.numFrame(iGroup) & ...
    msErrorType == groupKey.msErrorType(iGroup);
  idx = find(groupMask);
  [~, orderLocal] = sort(tailScore(idx), 'descend');
  keepCount = min(scanConfig.representativeSeedCountPerType, numel(orderLocal));
  if keepCount <= 0
    continue;
  end
  keepIdx = idx(orderLocal(1:keepCount));
  outTable = perfTable(keepIdx, {'displayName', 'snrDb', 'numFrame', 'taskSeed', ...
    'angleErrDeg', 'angleNormErr', 'fdRefAbsErrHz', 'fdRefNormErr', ...
    'fdRateAbsErrHzPerSec', 'nonRefCoherenceFloor', 'iterations', ...
    'firstOrderOpt', 'failureReason', 'tailSubtype', 'coreResolved', 'trimmedCore'});
  outTable.msErrorType = msErrorType(keepIdx);
  outTable.tailScore = tailScore(keepIdx);
  outTable.rankInType = reshape((1:keepCount), [], 1);
  outTable = movevars(outTable, {'msErrorType', 'tailScore', 'rankInType'}, 'After', 'numFrame');
  rowCell{end + 1, 1} = outTable; %#ok<AGROW>
end
if isempty(rowCell)
  return;
end
representativeSeedTable = vertcat(rowCell{:});
end

function msErrorType = localClassifyMsErrorType(row, scanConfig)
%LOCALCLASSIFYMSERRORTYPE Assign one offline MS error label.

msErrorType = "other-tail";
if logical(row.trimmedCore)
  msErrorType = "crb-local";
elseif ~logical(row.solverResolved)
  msErrorType = "solver-conditioning-tail";
elseif ~logical(row.toothResolved) || ...
    (isfinite(row.fdRefNormErr) && abs(row.fdRefNormErr) > scanConfig.trimNormCap)
  msErrorType = "fdref-branch-tail";
elseif ~logical(row.fdRateResolved) || logical(row.fdRateBoundaryHit) || logical(row.freqBoundaryHit)
  msErrorType = "fd-rate-tail";
elseif row.satMode == "multi" && logical(row.coherenceCollapsed)
  msErrorType = "coherence-collapse-tail";
elseif isfinite(row.angleNormErr) && abs(row.angleNormErr) > scanConfig.trimNormCap && ...
    isfinite(row.nonRefCoherenceFloor) && row.nonRefCoherenceFloor >= scanConfig.coreCoherenceFloor
  msErrorType = "doa-basin-limited";
elseif isfinite(row.angleNormErr) && abs(row.angleNormErr) > scanConfig.trimNormCap
  msErrorType = "angle-local-tail";
elseif (isfinite(row.firstOrderOpt) && row.firstOrderOpt > 10) || ...
    (isfinite(row.iterations) && row.iterations >= 100)
  msErrorType = "solver-conditioning-tail";
end
end

function tailScore = localBuildMsTailScore(row, msErrorType, scanConfig)
%LOCALBUILDMSTAILSCORE Build a sortable offline tail severity score.

scoreVec = [abs(row.angleNormErr), abs(row.fdRefNormErr)];
scoreVec = scoreVec(isfinite(scoreVec));
if isempty(scoreVec)
  tailScore = 0;
else
  tailScore = max(scoreVec);
end
if msErrorType == "fd-rate-tail"
  tailScore = max(tailScore, scanConfig.trimNormCap + 0.5);
elseif msErrorType == "coherence-collapse-tail"
  if isfinite(row.nonRefCoherenceFloor)
    tailScore = max(tailScore, scanConfig.trimNormCap + 1 - row.nonRefCoherenceFloor);
  else
    tailScore = max(tailScore, scanConfig.trimNormCap + 1);
  end
elseif msErrorType == "solver-conditioning-tail"
  tailScore = max(tailScore, scanConfig.trimNormCap + 0.25);
end
end

function plotData = localBuildPlotData(aggregateTable, failureSummaryTable, topTailTable)
%LOCALBUILDPLOTDATA Store lightweight data needed to redraw scan figures.

plotData = struct();
plotData.aggregateTable = aggregateTable;
plotData.failureSummaryTable = failureSummaryTable;
plotData.topTailTable = topTailTable;
end

function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw paired angle/fdRef consistency, MSE/CRB, and keep-rate figures.

plotData = scanData.plotData;
aggregateTable = scanData.aggregateTable;
if isempty(aggregateTable) || height(aggregateTable) == 0
  return;
end
frameList = unique(aggregateTable.numFrame, 'stable');
for iFrame = 1:numel(frameList)
  frameMask = aggregateTable.numFrame == frameList(iFrame);
  plotTable = aggregateTable(frameMask, :);
  methodList = unique(plotTable.displayName, 'stable');

  figure('Name', sprintf('%s RMSE and CRB P=%d', char(scanData.scanName), frameList(iFrame)));
  subplot(1, 2, 1);
  hold on;
  for iMethod = 1:numel(methodList)
    rowMask = plotTable.displayName == methodList(iMethod);
    semilogy(plotTable.snrDb(rowMask), plotTable.coreAngleRmseDeg(rowMask), '-o', ...
      'DisplayName', char(methodList(iMethod) + " core RMSE"));
    semilogy(plotTable.snrDb(rowMask), plotTable.trimAngleRmseDeg(rowMask), '--s', ...
      'DisplayName', char(methodList(iMethod) + " CRB-local RMSE"));
    semilogy(plotTable.snrDb(rowMask), plotTable.angleCrbStdDeg(rowMask), ':x', ...
      'DisplayName', char(methodList(iMethod) + " CRB"));
  end
  grid on;
  xlabel('SNR (dB)');
  ylabel('Global angle RMSE and CRB (deg)');
  title(sprintf('Global angle error, P=%d', frameList(iFrame)));
  legend('Location', 'best');

  subplot(1, 2, 2);
  hold on;
  for iMethod = 1:numel(methodList)
    rowMask = plotTable.displayName == methodList(iMethod);
    semilogy(plotTable.snrDb(rowMask), plotTable.coreFdRefRmseHz(rowMask), '-o', ...
      'DisplayName', char(methodList(iMethod) + " core RMSE"));
    semilogy(plotTable.snrDb(rowMask), plotTable.trimFdRefRmseHz(rowMask), '--s', ...
      'DisplayName', char(methodList(iMethod) + " CRB-local RMSE"));
    semilogy(plotTable.snrDb(rowMask), plotTable.fdRefCrbStdHz(rowMask), ':x', ...
      'DisplayName', char(methodList(iMethod) + " CRB"));
  end
  grid on;
  xlabel('SNR (dB)');
  ylabel('Reference Doppler RMSE and CRB (Hz)');
  title(sprintf('Reference Doppler error, P=%d', frameList(iFrame)));
  legend('Location', 'best');

  figure('Name', sprintf('%s MSE over CRB P=%d', char(scanData.scanName), frameList(iFrame)));
  subplot(1, 2, 1);
  hold on;
  for iMethod = 1:numel(methodList)
    rowMask = plotTable.displayName == methodList(iMethod);
    plot(plotTable.snrDb(rowMask), plotTable.coreAngleMseOverCrb(rowMask), '-o', ...
      'DisplayName', char(methodList(iMethod) + " core"));
    plot(plotTable.snrDb(rowMask), plotTable.trimAngleMseOverCrb(rowMask), '--s', ...
      'DisplayName', char(methodList(iMethod) + " CRB-local"));
  end
  grid on;
  xlabel('SNR (dB)');
  ylabel('Angle MSE / CRB');
  title(sprintf('Global angle MSE/CRB, P=%d', frameList(iFrame)));
  legend('Location', 'best');

  subplot(1, 2, 2);
  hold on;
  for iMethod = 1:numel(methodList)
    rowMask = plotTable.displayName == methodList(iMethod);
    plot(plotTable.snrDb(rowMask), plotTable.coreFdRefMseOverCrb(rowMask), '-o', ...
      'DisplayName', char(methodList(iMethod) + " core"));
    plot(plotTable.snrDb(rowMask), plotTable.trimFdRefMseOverCrb(rowMask), '--s', ...
      'DisplayName', char(methodList(iMethod) + " CRB-local"));
  end
  grid on;
  xlabel('SNR (dB)');
  ylabel('fdRef MSE / CRB');
  title(sprintf('Reference Doppler MSE/CRB, P=%d', frameList(iFrame)));
  legend('Location', 'best');

  figure('Name', sprintf('%s keep rates P=%d', char(scanData.scanName), frameList(iFrame)));
  hold on;
  for iMethod = 1:numel(methodList)
    rowMask = plotTable.displayName == methodList(iMethod);
    plot(plotTable.snrDb(rowMask), plotTable.resolvedRate(rowMask), '-o', ...
      'DisplayName', char(methodList(iMethod) + " loose"));
    plot(plotTable.snrDb(rowMask), plotTable.coreResolvedRate(rowMask), '--s', ...
      'DisplayName', char(methodList(iMethod) + " core"));
    plot(plotTable.snrDb(rowMask), plotTable.trimmedCoreRate(rowMask), ':x', ...
      'DisplayName', char(methodList(iMethod) + " CRB-local"));
  end
  grid on;
  ylim([0, 1.05]);
  xlabel('SNR (dB)');
  ylabel('Sample keep rate');
  title(sprintf('Resolved/core/CRB-local keep rates, P=%d', frameList(iFrame)));
  legend('Location', 'best');
end
plotData.aggregateTable = aggregateTable;
end

function repeatSlim = localStripRepeatOut(repeatOut)
%LOCALSTRIPREPEATOUT Keep lightweight per-repeat diagnostics in scanData.

repeatSlim = struct();
repeatSlim.taskSeed = repeatOut.taskSeed;
repeatSlim.snrDb = repeatOut.snrDb;
repeatSlim.toothStepHz = repeatOut.toothStepHz;
repeatSlim.truthFdRefHz = repeatOut.truthFdRefHz;
repeatSlim.truthFdRateHzPerSec = repeatOut.truthFdRateHzPerSec;
repeatSlim.fdRangeOracle = repeatOut.fdRangeOracle;
repeatSlim.fdRateRangeOracle = repeatOut.fdRateRangeOracle;
repeatSlim.crbTable = repeatOut.crbTable;
repeatSlim.caseRowList = repeatOut.caseRowList;
repeatSlim.repeatTotalMs = repeatOut.repeatTotalMs;
repeatSlim.warningSeen = repeatOut.warningSeen;
repeatSlim.warningId = repeatOut.warningId;
repeatSlim.warningMessage = repeatOut.warningMessage;
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve the Doppler tooth spacing from time offsets.

timeOffsetSec = [];
if isfield(fixture, 'sceneSeq') && isfield(fixture.sceneSeq, 'timeOffsetSec')
  timeOffsetSec = fixture.sceneSeq.timeOffsetSec(:);
end
if numel(timeOffsetSec) < 2
  toothStepHz = 750;
  return;
end
dt = median(diff(sort(timeOffsetSec)));
toothStepHz = 1 / dt;
end

function value = localResolveTruthScalar(truth, fieldNameList)
%LOCALRESOLVETRUTHSCALAR Resolve one scalar truth value from fallback fields.

value = NaN;
for iField = 1:numel(fieldNameList)
  fieldName = fieldNameList{iField};
  if isfield(truth, fieldName) && isscalar(truth.(fieldName)) && isfinite(truth.(fieldName))
    value = double(truth.(fieldName));
    return;
  end
end
end

function txt = localFormatRow(x)
%LOCALFORMATROW Format a numeric row for compact logs.

x = reshape(double(x), 1, []);
if isempty(x)
  txt = "";
else
  txt = strjoin(compose('%.6g', x), ', ');
end
end

function txt = localFormatStringRow(x)
%LOCALFORMATSTRINGROW Format a string row for compact logs.

x = reshape(string(x), 1, []);
if isempty(x)
  txt = "";
else
  txt = strjoin(x, ', ');
end
end

function lineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build HTML-ready key metric lines for notifications.

aggregateTable = scanData.aggregateTable;
if isempty(aggregateTable) || height(aggregateTable) == 0
  lineList = strings(0, 1);
  return;
end
framePick = max(aggregateTable.numFrame);
snrPick = max(aggregateTable.snrDb);
methodPick = "MS-MF-CP-U";
row = localFindAggregateRow(aggregateTable, methodPick, snrPick, framePick);
if isempty(row)
  row = aggregateTable(1, :);
end
lineList = strings(0, 1);
lineList(end + 1, 1) = sprintf('<code>%s</code> P=%d SNR=%.1f dB', ...
  char(row.displayName(1)), row.numFrame(1), row.snrDb(1));
lineList(end + 1, 1) = sprintf('loose/core/CRB-local rates=<code>%.3f</code>/<code>%.3f</code>/<code>%.3f</code>', ...
  row.resolvedRate(1), row.coreResolvedRate(1), row.crbLocalRate(1));
lineList(end + 1, 1) = sprintf('core angle/fdRef MSE-CRB=<code>%.3f</code>/<code>%.3f</code>', ...
  row.coreAngleMseOverCrb(1), row.coreFdRefMseOverCrb(1));
lineList(end + 1, 1) = sprintf('CRB-local angle/fdRef MSE-CRB=<code>%.3f</code>/<code>%.3f</code>', ...
  row.crbLocalAngleMseOverCrb(1), row.crbLocalFdRefMseOverCrb(1));
if isfield(scanData, 'msErrorTypeSummaryTable') && istable(scanData.msErrorTypeSummaryTable) && ...
    height(scanData.msErrorTypeSummaryTable) > 0
  typeTable = scanData.msErrorTypeSummaryTable;
  typeMask = typeTable.displayName == row.displayName(1) & typeTable.snrDb == row.snrDb(1) & ...
    typeTable.msErrorType ~= "crb-local";
  if any(typeMask)
    typeSub = typeTable(typeMask, :);
    [~, typeIdx] = max(typeSub.rate);
    lineList(end + 1, 1) = sprintf('top MS error=<code>%s</code> rate=<code>%.3f</code>', ...
      char(typeSub.msErrorType(typeIdx)), typeSub.rate(typeIdx));
  end
end
end

function row = localFindAggregateRow(aggregateTable, displayName, snrDb, numFrame)
%LOCALFINDAGGREGATEROW Return one aggregate row for compact notifications.

mask = aggregateTable.displayName == string(displayName) & aggregateTable.snrDb == snrDb & ...
  aggregateTable.numFrame == numFrame;
idx = find(mask, 1, 'first');
if isempty(idx)
  row = [];
else
  row = aggregateTable(idx, :);
end
end

function value = localMse(x)
%LOCALMSE Compute MSE over finite values only.

x = double(x(:));
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = mean(x.^2);
end
end

function value = localRmse(x)
%LOCALRMSE Compute RMSE over finite values only.

x = double(x(:));
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = sqrt(mean(x.^2));
end
end

function value = localPrctileFinite(x, p)
%LOCALPRCTILEFINITE Compute percentile over finite values only.

x = double(x(:));
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = prctile(x, p);
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide only when the denominator is finite and nonzero.

if isfinite(numerator) && isfinite(denominator) && abs(denominator) > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
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
%LOCALCANUSEPARFOR Check whether the outer task loop can use parfor.

useParfor = false;
if numTask < 2 || ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  useParfor = ~isempty(ver('parallel'));
catch
  useParfor = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field from a struct with a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end
