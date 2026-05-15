% scanMfMsMleCrbCleanBoundScale
% Purpose: scan truth-centered DoA and fdRef CRB-scale search boxes for
% clean SS/MS MF MLE-vs-CRB local-bound consistency.
%
% Usage: edit the Scan configuration block, then run this script directly.
% The script leaves scanData in the workspace; checkpointEnable=true resumes
% interrupted task runs from per-task files under the repo-root tmp.
% saveSnapshot=true saves only scanData via saveExpSnapshot. Telegram notice
% is best-effort only and never affects numerical results.

clear; close all; clc;

%% Scan configuration

scanName = "scanMfMsMleCrbCleanBoundScale";
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = true;
optVerbose = false;

% Dynamic MF estimator dialect. The default is the clean MLE core: no CP
% engineering penalties, no fd alias unwrap, and no estimator-internal
% basin-entry / warm-anchor route. Keep this scan on core-only unless
% temporarily running a manual route diagnostic.
dynamicObjectiveMode = "pure-mle";
dynamicRouteMode = "core-only";
% Paper-model override for this scan only: fixed steering at the reference
% frame and affine reference-tied Doppler over absolute time, matching
% replayMfMsMleCrbCleanTrim.
cleanMfSignalTimeModel = "linearRef";
cleanMfSteeringMode = "frozenRef";

numFrame = 10;
contextUsrLla = [55; 36.59; 0];
contextUtcRef = datetime([2026, 03, 18, 17, 08, 00], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
contextTleFileName = "statlink_20260318.tle";
contextSelectedSatIdxGlobal = [5259, 1243, 348, 5652, 14, 4437];
contextRefSatIdxGlobal = 5259;

% This scan is a range-surface diagnostic, not an SNR curve by default.
% Increase numRepeat or expand snrDbList after the surface shape is stable.
snrDbList = 0;
baseSeed = 253;
numRepeat = 30;
seedList = baseSeed + (0:(numRepeat - 1));
methodNameList = [
  "SS-MF-CP-K"
  "SS-MF-CP-U"
  "MS-MF-CP-K"
  "MS-MF-CP-U"
];

% CRB-scaled truth-centered hard boxes. The scan axis is intentionally small
% and sparse first; refine locally after finding the knee / transition region.
truthLocalDoaCrbScaleList = [1.75; 2; 2.5; 3];
truthLocalFdRefCrbScaleList = [2.25; 2.5; 2.75; 3; 3.5];
truthLocalFdRateHalfWidthHzPerSec = 1000;
fdRefRangeMode = "crb-scale";
fdRefMaxHalfToothFraction = 0.45;

boundaryTolFraction = 0.02;
healthFdRefToothFraction = 0.25;
healthFdRateAbsMaxHzPerSec = 250;
trimAngleNormMax = 5;
trimFdRefNormMax = 5;
crbTargetMseRatio = 1.1;

seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
methodNameList = reshape(string(methodNameList), [], 1);
truthLocalDoaCrbScaleList = reshape(double(truthLocalDoaCrbScaleList), [], 1);
truthLocalFdRefCrbScaleList = reshape(double(truthLocalFdRefCrbScaleList), [], 1);
if isempty(truthLocalDoaCrbScaleList) || any(~isfinite(truthLocalDoaCrbScaleList) | truthLocalDoaCrbScaleList <= 0)
  error('scanMfMsMleCrbCleanBoundScale:InvalidDoaScaleList', ...
    'truthLocalDoaCrbScaleList must contain positive finite values.');
end
if isempty(truthLocalFdRefCrbScaleList) || any(~isfinite(truthLocalFdRefCrbScaleList) | truthLocalFdRefCrbScaleList <= 0)
  error('scanMfMsMleCrbCleanBoundScale:InvalidFdRefScaleList', ...
    'truthLocalFdRefCrbScaleList must contain positive finite values.');
end

runTic = tic;
scanData = struct();
config = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and scan options

  config.scanName = string(scanName);
  config.numFrame = numFrame;
  config.snrDbList = snrDbList;
  config.seedList = seedList;
  config.baseSeed = seedList(1);
  config.numRepeat = numel(seedList);
  config.methodNameList = methodNameList;
  config.truthLocalDoaCrbScaleList = truthLocalDoaCrbScaleList;
  config.truthLocalFdRefCrbScaleList = truthLocalFdRefCrbScaleList;
  config.truthLocalFdRateHalfWidthHzPerSec = double(truthLocalFdRateHalfWidthHzPerSec);
  config.fdRefRangeMode = string(fdRefRangeMode);
  config.fdRefMaxHalfToothFraction = double(fdRefMaxHalfToothFraction);
  config.scaleProfileMode = "doa-fdref-2d";
  config.boundaryTolFraction = boundaryTolFraction;
  config.healthFdRefToothFraction = healthFdRefToothFraction;
  config.healthFdRateAbsMaxHzPerSec = healthFdRateAbsMaxHzPerSec;
  config.trimAngleNormMax = trimAngleNormMax;
  config.trimFdRefNormMax = trimFdRefNormMax;
  config.crbTargetMseRatio = crbTargetMseRatio;
  config.saveSnapshot = logical(saveSnapshot);
  config.checkpointEnable = logical(checkpointEnable);
  config.checkpointResume = logical(checkpointResume);
  config.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
  config.notifyTelegramEnable = logical(notifyTelegramEnable);
  config.optVerbose = logical(optVerbose);
  config.dynamicObjectiveMode = localNormalizeDynamicObjectiveMode(dynamicObjectiveMode);
  config.dynamicRouteMode = localNormalizeDynamicRouteMode(dynamicRouteMode);
  config.initMode = "truth-crb-scale-sweep-static-seed-fixed-bound";
  config.sweepMode = "snr-seed-doa-fdref-scale";
  config.routeTag = "clean-bound-scale-l" + string(numel(contextSelectedSatIdxGlobal)) + ...
    "-ref" + string(contextRefSatIdxGlobal) + "-" + ...
    config.dynamicObjectiveMode + "-" + config.dynamicRouteMode + "-" + ...
    string(cleanMfSignalTimeModel) + "-" + string(cleanMfSteeringMode);
  config.contextUsrLla = reshape(double(contextUsrLla), [], 1);
  config.contextUtcRef = contextUtcRef;
  config.contextTleFileName = string(contextTleFileName);
  config.contextSelectedSatIdxGlobal = reshape(double(contextSelectedSatIdxGlobal), 1, []);
  config.contextRefSatIdxGlobal = double(contextRefSatIdxGlobal);
  config.cleanMfSignalTimeModel = localNormalizeCleanMfSignalTimeModel(cleanMfSignalTimeModel);
  config.cleanMfSteeringMode = localNormalizeCleanMfSteeringMode(cleanMfSteeringMode);

  taskList = localBuildTaskList(config.snrDbList, config.seedList, ...
    config.truthLocalDoaCrbScaleList, config.truthLocalFdRefCrbScaleList);
  config.numTask = numel(taskList);
  config.useParfor = localCanUseParfor(config.numTask);
  config.runKey = localBuildRunKey(config);
  checkpointOpt = localBuildCheckpointOpt(scanName, config, config.useParfor);
  if config.checkpointEnable
    checkpointRunDir = string(checkpointOpt.runDir);
    config.checkpointRunDir = checkpointRunDir;
  end

  printMfScanHeader(char(scanName), config, char(checkpointRunDir));
  localPrintReplayConfig(config);

  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(config.numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.usrLla = config.contextUsrLla;
  contextOpt.utcRef = config.contextUtcRef;
  contextOpt.tleFileName = config.contextTleFileName;
  contextOpt.selectedSatIdxGlobal = config.contextSelectedSatIdxGlobal;
  contextOpt.refSatIdxGlobal = config.contextRefSatIdxGlobal;
  contextOpt.satPickCount = numel(config.contextSelectedSatIdxGlobal);
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = localBuildSerialInnerParallelOpt();
  context = buildDynamicMultiSatEciContext(contextOpt);
  context = localApplyCleanSignalModelOverride(context, config);
  context.subsetOffsetCell = {};
  context.subsetLabelList = strings(0, 1);
  context.numSubsetRandomTrial = 0;
  flowOpt = localBuildFlowOpt(context.parallelOpt);

  %% Run clean scan batch

  try
    [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, scanName);
  catch ME
    if strlength(string(checkpointRunDir)) > 0
      fprintf('Scan failed. Checkpoint artifacts kept at: %s\n', char(checkpointRunDir));
    end
    rethrow(ME);
  end

  caseTable = localBuildCaseTable(repeatCell);
  caseTable = localAnnotateHealthFlags(caseTable, config);
  scaleRawAggregateTable = localBuildScaleAggregateTable(caseTable, "raw", true(height(caseTable), 1));
  scaleHealthAggregateTable = localBuildScaleAggregateTable(caseTable, "health", caseTable.healthResolved);
  scaleJointTrimAggregateTable = localBuildScaleAggregateTable(caseTable, "joint-trim", caseTable.jointTrimKeep);
  scaleReadableRmseCrbTable = localBuildScaleReadableRmseCrbTable( ...
    scaleRawAggregateTable, scaleJointTrimAggregateTable);
  scaleSurfaceSummaryTable = localBuildScaleSurfaceSummaryTable(scaleJointTrimAggregateTable);
  scaleSearchRangeAuditTable = localBuildScaleSearchRangeAuditTable(caseTable);
  scaleOutlierTable = localBuildScaleOutlierTable(caseTable);
  runtimeTable = buildMfRuntimeTable(repeatCell);
  runtimeAggregateTable = buildMfRuntimeAggregateTable(runtimeTable);
  topSlowRuntimeTable = buildMfTopSlowRuntimeTable(runtimeTable, 12);

  %% Data storage

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = config.runKey;
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = config;
  scanData.contextSummary = localBuildContextSummary(context);
  scanData.caseTable = caseTable;
  scanData.scaleRawAggregateTable = scaleRawAggregateTable;
  scanData.scaleHealthAggregateTable = scaleHealthAggregateTable;
  scanData.scaleJointTrimAggregateTable = scaleJointTrimAggregateTable;
  scanData.scaleReadableRmseCrbTable = scaleReadableRmseCrbTable;
  scanData.scaleSurfaceSummaryTable = scaleSurfaceSummaryTable;
  scanData.scaleSearchRangeAuditTable = scaleSearchRangeAuditTable;
  scanData.scaleOutlierTable = scaleOutlierTable;
  scanData.runtimeTable = runtimeTable;
  scanData.runtimeAggregateTable = runtimeAggregateTable;
  scanData.topSlowRuntimeTable = topSlowRuntimeTable;
  scanData.repeatCell = localStripRepeatCell(repeatCell);
  scanData.checkpointDir = string(checkpointRunDir);
  scanData.checkpointCleanupReport = struct([]);
  scanData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  if config.checkpointEnable && localGetFieldOrDefault(runState, 'isComplete', false) && ...
      config.checkpointCleanupOnSuccess
    scanData.checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(runState, struct('verbose', false));
    scanData.checkpointSummary.cleanupReport = scanData.checkpointCleanupReport;
    scanData.checkpointSummary.cleanedOnSuccess = true;
  elseif isfield(scanData.checkpointSummary, 'cleanedOnSuccess')
    scanData.checkpointSummary.cleanedOnSuccess = false;
  end
  scanData.checkpointSummaryTable = localBuildCheckpointSummaryTable(scanData.checkpointSummary, ...
    scanData.checkpointDir);
  scanData.plotData = localBuildScalePlotData(scaleJointTrimAggregateTable);
  scanData.elapsedSec = toc(runTic);
  scanData = finalizeMfScanResult(scanData, "");

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output

  if ~exist('scanData', 'var') || ~isstruct(scanData)
    error('scanMfMsMleCrbCleanBoundScale:MissingScanData', ...
      'Load a snapshot containing scanData before running the summary section.');
  end
  scanData = localValidateScanDataForSummary(scanData);
  scanConfigForReport = localGetFieldOrDefault(scanData, 'config', struct());
  scanNameForReport = string(localGetFieldOrDefault(scanData, 'scanName', "scanMfMsMleCrbCleanBoundScale"));
  scanSnapshotFile = string(localGetFieldOrDefault(scanData, 'snapshotFile', ""));
  scanElapsedSec = localGetFieldOrDefault(scanData, 'elapsedSec', NaN);
  checkpointDirForReport = string(localGetFieldOrDefault(scanData, 'checkpointDir', ""));

  printMfScanSection('Scale raw aggregate', scanData.scaleRawAggregateTable);
  printMfScanSection('Scale health aggregate', scanData.scaleHealthAggregateTable);
  printMfScanSection('Scale joint-trim aggregate', scanData.scaleJointTrimAggregateTable);
  printMfScanSection('Scale surface summary', scanData.scaleSurfaceSummaryTable);
  printMfScanSection('Scale search-range audit', scanData.scaleSearchRangeAuditTable);
  localPrintCompactTable('Scale outlier table', scanData.scaleOutlierTable, 12);
  printMfScanSection('Scale final readable RMSE/CRB summary', scanData.scaleReadableRmseCrbTable);
  printMfScanSection('Clean runtime aggregate', scanData.runtimeAggregateTable);
  localPrintCompactTable('Top slow clean tasks', scanData.topSlowRuntimeTable, 12);
  if isfield(scanData, 'checkpointSummaryTable') && istable(scanData.checkpointSummaryTable) && ...
      height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  fprintf('Scan purpose                    : clean truth-centered DoA/fdRef CRB-scale surface.\n');
  trimAngleNormMaxForReport = localGetFieldOrDefault(scanConfigForReport, 'trimAngleNormMax', NaN);
  trimFdRefNormMaxForReport = localGetFieldOrDefault(scanConfigForReport, 'trimFdRefNormMax', NaN);
  fprintf('Trim definition                   : health + angle <= %.2f CRB + fdRef <= %.2f CRB.\n', ...
    trimAngleNormMaxForReport, trimFdRefNormMaxForReport);
  fprintf('Main decision table               : scaleReadableRmseCrbTable.\n');
  fprintf('Surface table                     : scaleSurfaceSummaryTable.\n');
  scanData.plotData = localPlotScaleScan(scanData);

  notifyMfScanStatus(struct( ...
    'scanName', scanNameForReport, ...
    'statusText', "DONE", ...
    'config', scanConfigForReport, ...
    'snapshotFile', scanSnapshotFile, ...
    'checkpointDir', checkpointDirForReport, ...
    'elapsedSec', scanElapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Clean-bound scale scan completed; inspect scaleReadableRmseCrbTable before changing estimator logic."]));

catch ME
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'snapshotFile', "", ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'metricLineList', string(ME.identifier), ...
    'commentLineList', string(ME.message), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function scanData = localValidateScanDataForSummary(scanData)
%LOCALVALIDATESCANDATAFORSUMMARY Validate stored scan data for rerunnable summary.

if ~isfield(scanData, 'caseTable') || ~istable(scanData.caseTable)
  error('scanMfMsMleCrbCleanBoundScale:MissingCaseTable', ...
    'scanData.caseTable is missing. Store caseTable inside scanData before saving.');
end
if ~isfield(scanData, 'config') || ~isstruct(scanData.config)
  scanData.config = struct();
end
scanData.caseTable = localAnnotateHealthFlags(scanData.caseTable, scanData.config);
if ~isfield(scanData, 'scaleRawAggregateTable') || ~istable(scanData.scaleRawAggregateTable)
  scanData.scaleRawAggregateTable = localBuildScaleAggregateTable(scanData.caseTable, ...
    "raw", true(height(scanData.caseTable), 1));
end
if ~isfield(scanData, 'scaleHealthAggregateTable') || ~istable(scanData.scaleHealthAggregateTable)
  scanData.scaleHealthAggregateTable = localBuildScaleAggregateTable(scanData.caseTable, ...
    "health", scanData.caseTable.healthResolved);
end
if ~isfield(scanData, 'scaleJointTrimAggregateTable') || ~istable(scanData.scaleJointTrimAggregateTable)
  scanData.scaleJointTrimAggregateTable = localBuildScaleAggregateTable(scanData.caseTable, ...
    "joint-trim", scanData.caseTable.jointTrimKeep);
end
if ~isfield(scanData, 'scaleReadableRmseCrbTable') || ~istable(scanData.scaleReadableRmseCrbTable)
  scanData.scaleReadableRmseCrbTable = localBuildScaleReadableRmseCrbTable( ...
    scanData.scaleRawAggregateTable, scanData.scaleJointTrimAggregateTable);
end
if ~isfield(scanData, 'scaleSurfaceSummaryTable') || ~istable(scanData.scaleSurfaceSummaryTable)
  scanData.scaleSurfaceSummaryTable = localBuildScaleSurfaceSummaryTable(scanData.scaleJointTrimAggregateTable);
end
if ~isfield(scanData, 'scaleSearchRangeAuditTable') || ~istable(scanData.scaleSearchRangeAuditTable)
  scanData.scaleSearchRangeAuditTable = localBuildScaleSearchRangeAuditTable(scanData.caseTable);
end
if ~isfield(scanData, 'scaleOutlierTable') || ~istable(scanData.scaleOutlierTable)
  scanData.scaleOutlierTable = localBuildScaleOutlierTable(scanData.caseTable);
end
if ~isfield(scanData, 'runtimeAggregateTable') || ~istable(scanData.runtimeAggregateTable)
  scanData.runtimeAggregateTable = table();
end
if ~isfield(scanData, 'topSlowRuntimeTable') || ~istable(scanData.topSlowRuntimeTable)
  scanData.topSlowRuntimeTable = table();
end
if ~isfield(scanData, 'checkpointDir')
  scanData.checkpointDir = "";
end
if ~isfield(scanData, 'checkpointSummaryTable')
  checkpointSummary = localGetFieldOrDefault(scanData, 'checkpointSummary', struct());
  scanData.checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointSummary, scanData.checkpointDir);
end
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  scanData.plotData = localBuildScalePlotData(scanData.scaleJointTrimAggregateTable);
end
end

function taskList = localBuildTaskList(snrDbList, seedList, doaScaleList, fdRefScaleList)
%LOCALBUILDTASKLIST Build one task per SNR, seed, DoA scale, and fdRef scale.

taskList = repmat(struct('taskId', 0, 'snrDb', NaN, 'taskSeed', NaN, ...
  'doaCrbScale', NaN, 'fdRefCrbScale', NaN), 0, 1);
taskId = 0;
for iSnr = 1:numel(snrDbList)
  for iDoa = 1:numel(doaScaleList)
    for iFd = 1:numel(fdRefScaleList)
      for iSeed = 1:numel(seedList)
        taskId = taskId + 1;
        taskList(taskId, 1).taskId = taskId; %#ok<AGROW>
        taskList(taskId, 1).snrDb = snrDbList(iSnr); %#ok<AGROW>
        taskList(taskId, 1).taskSeed = seedList(iSeed); %#ok<AGROW>
        taskList(taskId, 1).doaCrbScale = doaScaleList(iDoa); %#ok<AGROW>
        taskList(taskId, 1).fdRefCrbScale = fdRefScaleList(iFd); %#ok<AGROW>
      end
    end
  end
end
end

function rangeCfg = localResolveMethodRangeConfig(config, methodName, crbBundle, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz, task)
%LOCALRESOLVEMETHODRANGECONFIG Resolve one scale-sweep local search range.

methodName = string(methodName);
doaCrbScale = double(task.doaCrbScale);
fdRefCrbScale = double(task.fdRefCrbScale);
fdRateHalfWidthHzPerSec = double(localGetFieldOrDefault(config, ...
  'truthLocalFdRateHalfWidthHzPerSec', 1000));
if ~(isfinite(doaCrbScale) && doaCrbScale > 0 && ...
    isfinite(fdRefCrbScale) && fdRefCrbScale > 0 && ...
    isfinite(fdRateHalfWidthHzPerSec) && fdRateHalfWidthHzPerSec > 0)
  error('scanMfMsMleCrbCleanBoundScale:InvalidScaleTask', ...
    'Each scale-sweep task must have positive finite DoA, fdRef, and fdRate bounds.');
end
fdRangeUse = localBuildFdRefCrbScaledRange(truthFdRefHz, crbBundle, methodName, ...
  fdRefCrbScale, toothStepHz, localGetFieldOrDefault(config, 'fdRefMaxHalfToothFraction', 0.45));
rangeCfg = struct();
rangeCfg.methodName = methodName;
rangeCfg.doaCrbScale = doaCrbScale;
rangeCfg.fdRefHalfToothFraction = localSafeRatio(0.5 * abs(fdRangeUse(2) - fdRangeUse(1)), toothStepHz);
rangeCfg.fdRefCrbScale = fdRefCrbScale;
rangeCfg.fdRateHalfWidthHzPerSec = fdRateHalfWidthHzPerSec;
rangeCfg.fdRangeRequested = fdRangeUse;
rangeCfg.fdRangeUse = fdRangeUse;
rangeCfg.fdRateRangeUse = truthFdRateHzPerSec + fdRateHalfWidthHzPerSec * [-1, 1];
end

function [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, scanName)
%LOCALRUNTASKBATCH Run clean tasks with optional checkpoint/resume.

numTask = numel(taskList);
runState = struct();
progressTracker = localCreateProgressTracker(numTask);
progressCleanup = onCleanup(@() localCloseProgressTracker(progressTracker)); %#ok<NASGU>
try
  if logical(localGetFieldOrDefault(config, 'checkpointEnable', false))
    checkpointOpt = localBuildCheckpointOpt(scanName, config, localGetFieldOrDefault(config, 'useParfor', false));
    taskDir = fullfile(checkpointOpt.runDir, 'task');
    numDoneTask = localCountCheckpointTaskFile(taskDir, numTask);
    localAdvanceProgressByStep(progressTracker, numDoneTask);
    sharedData = struct('context', context, 'flowOpt', flowOpt, 'config', config);
    runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker);
    runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, @localCheckpointTaskRunner, runnerOpt);
    repeatCell = runState.resultCell;
    return;
  end

  repeatCell = cell(numTask, 1);
  useParfor = logical(localGetFieldOrDefault(config, 'useParfor', false));
  if useParfor && isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
    progressQueue = progressTracker.queue;
    parfor iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
      send(progressQueue, 1);
    end
  elseif useParfor
    parfor iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
    end
  else
    for iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
      localAdvanceProgress(progressTracker);
    end
  end
catch ME
  localCloseProgressTracker(progressTracker);
  rethrow(ME);
end
end

function runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only fields accepted by checkpoint runner.

runnerOpt = struct();
runnerOpt.runName = checkpointOpt.runName;
runnerOpt.runKey = checkpointOpt.runKey;
runnerOpt.outputRoot = checkpointOpt.outputRoot;
runnerOpt.useParfor = logical(localGetFieldOrDefault(checkpointOpt, 'useParfor', false));
runnerOpt.resume = checkpointOpt.resume;
runnerOpt.meta = checkpointOpt.meta;
runnerOpt.progressFcn = @(step) localAdvanceProgressByStep(progressTracker, step);
runnerOpt.cleanupOnSuccess = false;
runnerOpt.cleanupOpt = struct();
end

function numDone = localCountCheckpointTaskFile(taskDir, numTask)
%LOCALCOUNTCHECKPOINTTASKFILE Count completed checkpoint task files before progress reset.

numDone = 0;
if ~isfolder(taskDir)
  return;
end
for iTask = 1:numTask
  if isfile(fullfile(taskDir, sprintf('task_%06d.mat', iTask)))
    numDone = numDone + 1;
  end
end
end

function localAdvanceProgressByStep(progressTracker, step)
%LOCALADVANCEPROGRESSBYSTEP Advance progressbar by a task count.

for iStep = 1:max(0, round(step))
  localAdvanceProgress(progressTracker);
end
end

function repeatOut = localCheckpointTaskRunner(taskInfo, sharedData)
%LOCALCHECKPOINTTASKRUNNER Run one checkpointed clean task.

repeatOut = localRunOneTask(sharedData.context, sharedData.flowOpt, sharedData.config, taskInfo);
end

function repeatOut = localRunOneTask(context, flowOpt, config, task)
%LOCALRUNONETASK Run one SNR/seed clean truth-centered task.

lastwarn('', '');
tRepeat = tic;
runtimeRows = repmat(emptyMfRuntimeRow(), 0, 1);

tStage = tic;
repeatData = buildDynamicRepeatData(context, task.snrDb, task.taskSeed);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "repeat-data", "setup", 1000 * toc(tStage)); %#ok<AGROW>
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('scanMfMsMleCrbCleanBoundScale:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for this truth-centered clean scan.');
end

noiseVar = 1 / (10^(task.snrDb / 10));
tStage = tic;
crbBundle = localBuildCrbBundleQuiet(fixture, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, noiseVar);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "crb-bundle", "crb", 1000 * toc(tStage)); %#ok<AGROW>

tStage = tic;
staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz, config, task);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "static-bundle", "static", 1000 * toc(tStage)); %#ok<AGROW>

caseRowList = repmat(localEmptyCaseRow(), numel(config.methodNameList), 1);
caseMap = struct();
for iMethod = 1:numel(config.methodNameList)
  methodName = config.methodNameList(iMethod);
  rangeCfg = localResolveMethodRangeConfig(config, methodName, crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz, task);
  tMethod = tic;
  lastwarn('', '');
  if methodName == "SS-SF-Static"
    caseUse = staticBundle.caseStaticRefOnly;
    initMeta = staticBundle.initMetaStaticRefOnly;
    methodRuntimeRows = repmat(emptyMfRuntimeRow(), 0, 1);
  elseif methodName == "MS-SF-Static"
    caseUse = staticBundle.caseStaticMs;
    initMeta = staticBundle.initMetaStaticMs;
    methodRuntimeRows = repmat(emptyMfRuntimeRow(), 0, 1);
  else
    method = localBuildDynamicMethod(methodName, rangeCfg.doaCrbScale, config.dynamicRouteMode);
    method = localApplyDoaCrbScaledWidth(method, crbBundle);
    staticSeedCase = localSelectStaticSeedCase(method, staticBundle);
    [caseUse, initMeta, methodRuntimeRows] = localRunDynamicMethod(method, fixture, truth, context, flowOpt, ...
      rangeCfg.fdRangeUse, rangeCfg.fdRateRangeUse, toothStepHz, config.optVerbose, staticSeedCase, config, task);
  end
  runtimeRows = [runtimeRows; methodRuntimeRows(:)]; %#ok<AGROW>
  row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
    rangeCfg.fdRangeUse, rangeCfg.fdRateRangeUse, crbBundle, config, initMeta, rangeCfg);
  [warningMessage, warningId] = lastwarn;
  row.warningFlag = ~(isempty(warningMessage) && isempty(warningId));
  caseRowList(iMethod) = row;
  caseMap.(localMethodKey(methodName)) = caseUse; %#ok<NASGU>
  runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, methodName, "method-total", "method", 1000 * toc(tMethod)); %#ok<AGROW>
end

[warningMessage, warningId] = lastwarn;
repeatOut = struct();
repeatOut.taskSeed = task.taskSeed;
repeatOut.snrDb = task.snrDb;
repeatOut.doaCrbScale = task.doaCrbScale;
repeatOut.fdRefCrbScale = task.fdRefCrbScale;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.toothStepHz = toothStepHz;
repeatOut.caseRowList = caseRowList;
repeatOut.runtimeRowList = runtimeRows;
repeatOut.repeatTotalMs = 1000 * toc(tRepeat);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz, config, task)
%LOCALBUILDCLEANSTATICBUNDLE Run only the static seeds required by selected methods.

truthDoaParam = reshape(truth.latlonTrueDeg, [], 1);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundScale:MissingTruthDoa', ...
    'Truth lat/lon is required to build the clean truth-centered static search box.');
end
truthDoaParam = truthDoaParam(1:2);

methodNameList = reshape(string(config.methodNameList), [], 1);
needStaticRefOnly = any(methodNameList == "SS-SF-Static" | startsWith(methodNameList, "SS-MF-"));
needStaticMs = any(methodNameList == "MS-SF-Static" | startsWith(methodNameList, "MS-MF-"));

staticBundle = struct();
if needStaticRefOnly
  rangeCfgRef = localResolveMethodRangeConfig(config, "SS-SF-Static", crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz, task);
  [staticBundle.caseStaticRefOnly, staticBundle.initMetaStaticRefOnly] = localRunCleanStaticCase( ...
    "SS-SF-Static", "single", fixture.viewRefOnly, context, flowOpt.staticBaseOpt, ...
    crbBundle, truthDoaParam, rangeCfgRef.fdRangeUse, rangeCfgRef.doaCrbScale, config);
end
if needStaticMs
  rangeCfgMs = localResolveMethodRangeConfig(config, "MS-SF-Static", crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz, task);
  [staticBundle.caseStaticMs, staticBundle.initMetaStaticMs] = localRunCleanStaticCase( ...
    "MS-SF-Static", "multi", fixture.viewMs, context, flowOpt.staticBaseOpt, ...
    crbBundle, truthDoaParam, rangeCfgMs.fdRangeUse, rangeCfgMs.doaCrbScale, config);
  staticBundle.bestStaticMsCase = staticBundle.caseStaticMs;
end
end

function [caseInfo, initMeta] = localRunCleanStaticCase(displayName, satMode, viewUse, context, staticBaseOpt, ...
  crbBundle, truthDoaParam, fdRangeUse, doaCrbScale, config)
%LOCALRUNCLEANSTATICCASE Run one SF static estimator inside the CRB-scaled truth-centered box.

modelOpt = staticBaseOpt;
modelOpt.initDoaParam = truthDoaParam(:);
modelOpt.initDoaHalfWidth = localResolveDoaCrbScaledHalfWidth(displayName, doaCrbScale, crbBundle);
[estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  viewUse.sceneRef, viewUse.rxSigSf, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, viewUse.doaGrid, fdRangeUse, 1, [], ...
  config.optVerbose, modelOpt);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa-doppler", "static", estResult);
initMeta = localBuildStaticInitMeta(estResult, modelOpt);
end

function initMeta = localBuildStaticInitMeta(estResult, modelOpt)
%LOCALBUILDSTATICINITMETA Record static estimator initial state and local box.

initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
initMeta = struct();
initMeta.requestedInitMode = "truth-crb-scaled-static";
initMeta.requestedInitParam = initParam(:);
initMeta.requestedInitDoaParam = reshape(modelOpt.initDoaParam, [], 1);
initMeta.initParam = initParam(:);
initMeta.initDoaParam = reshape(localGetFieldOrDefault(modelOpt, 'initDoaParam', NaN(2, 1)), [], 1);
if numel(initParam) >= 2 && all(isfinite(initParam(1:2)))
  initMeta.initDoaParam = initParam(1:2);
end
initMeta.initFdRefHz = NaN;
if numel(initParam) >= 3
  initMeta.initFdRefHz = initParam(3);
end
initMeta.initFdRateHzPerSec = NaN;
initMeta.searchHalfWidthDeg = reshape(localGetFieldOrDefault(modelOpt, 'initDoaHalfWidth', NaN(2, 1)), [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
end

function [caseUse, initMeta, runtimeRows] = localRunDynamicMethod(method, fixture, truth, context, flowOpt, ...
  fdRangeUse, fdRateRangeUse, toothStepHz, optVerbose, staticSeedCase, config, task)
%LOCALRUNDYNAMICMETHOD Run one clean MF estimator path.

runtimeRows = repmat(emptyMfRuntimeRow(), 0, 1);
if string(method.satMode) == "multi"
  viewUse = fixture.viewMs;
  debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());
else
  viewUse = fixture.viewRefOnly;
  debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());
end

fdRateRangeForMethod = fdRateRangeUse;
if method.isKnownRate
  fdRateRangeForMethod = [];
end

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.fdRateMode = char(method.fdRateMode);
if method.isKnownRate
  dynOpt.fdRateKnown = truth.fdRateFit;
end
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt = localApplyDynamicObjectiveMode(dynOpt, config);
dynOpt = localApplyDynamicRouteMode(dynOpt, method.dynamicRouteMode);
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
[initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, staticSeedCase, truth, config);
if isfield(initRequestMeta, 'initDoaParam') && ~isempty(initRequestMeta.initDoaParam)
  dynOpt.initDoaParam = initRequestMeta.initDoaParam(:);
  dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
  initRequestMeta.initDoaHalfWidth = method.doaHalfWidthDeg(:);
end

tStage = tic;
caseUse = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
  viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fdRangeUse, fdRateRangeForMethod, optVerbose, dynOpt, method.isKnownRate, ...
  debugTruthUse, initParamOverride);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, method.displayName, "baseline-solve", "estimator", 1000 * toc(tStage)); %#ok<AGROW>
caseUse.inToothMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, ...
  'fdRateRangeUse', fdRateRangeForMethod);
initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta);
end


function mode = localNormalizeCleanMfSignalTimeModel(mode)
%LOCALNORMALIZECLEANMFSIGNALTIMEMODE Normalize the scan signal phase model.

mode = lower(strtrim(string(mode)));
switch mode
  case {"linearref", "linear-ref", "linear_ref"}
    mode = "linearRef";
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidSignalTimeModel', ...
      'Unsupported cleanMfSignalTimeModel: %s.', char(mode));
end
end

function mode = localNormalizeCleanMfSteeringMode(mode)
%LOCALNORMALIZECLEANMFSTEERINGMODE Normalize the scan steering model.

mode = lower(strtrim(string(mode)));
switch mode
  case {"frozenref", "frozen-ref", "frozen_ref"}
    mode = "frozenRef";
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidSteeringMode', ...
      'Unsupported cleanMfSteeringMode: %s.', char(mode));
end
end

function context = localApplyCleanSignalModelOverride(context, config)
%LOCALAPPLYCLEANSIGNALMODELOVERRIDE Force the scan-only paper signal model.

if ~isfield(context, 'simOpt') || isempty(context.simOpt)
  context.simOpt = struct();
end
if isfield(context.simOpt, 'spatial')
  context.simOpt = rmfield(context.simOpt, 'spatial');
end
if ~isfield(context.simOpt, 'phase') || isempty(context.simOpt.phase)
  context.simOpt.phase = struct();
end
context.simOpt.phase.timeModel = char(config.cleanMfSignalTimeModel);
context.simOpt.phase.frameModel = 'shared';
context.simOpt.steeringMode = char(config.cleanMfSteeringMode);
context.simOpt.accessMode = 'framewise';
context.simOpt.steeringRefFrameIdx = context.sceneSeqMaster.refFrameIdx;
if ~isfield(context.simOpt, 'wave') || isempty(context.simOpt.wave)
  context.simOpt.wave = struct();
end
context.simOpt.wave.delayModel = 'phaseOnly';
context.simOpt.wave.timeRef = 'zero';
context.simOpt.wave.carrierPhaseModel = 'none';
end

function mode = localNormalizeDynamicObjectiveMode(mode)
%LOCALNORMALIZEDYNAMICOBJECTIVEMODE Normalize the dynamic objective dialect tag.

mode = lower(strtrim(string(mode)));
switch mode
  case {"pure-mle", "puremle", "pure"}
    mode = "pure-mle";
  case {"flow-default", "flowdefault", "robust", "engineering"}
    mode = "flow-default";
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidDynamicObjectiveMode', ...
      'Unsupported dynamicObjectiveMode: %s.', char(mode));
end
end

function dynOpt = localApplyDynamicObjectiveMode(dynOpt, config)
%LOCALAPPLYDYNAMICOBJECTIVEMODE Apply pure-MLE or robust objective options.

switch string(config.dynamicObjectiveMode)
  case "pure-mle"
    dynOpt.continuousPhaseConsistencyWeight = 0;
    dynOpt.continuousPhaseCollapsePenaltyWeight = 0;
    dynOpt.continuousPhaseNegativeProjectionPenaltyWeight = 0;
    dynOpt.continuousPhaseNonRefFitFloorWeight = 0;
    dynOpt.continuousPhaseNonRefSupportFloorWeight = 0;
    dynOpt.enableFdAliasUnwrap = false;
  case "flow-default"
    dynOpt.enableFdAliasUnwrap = true;
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidDynamicObjectiveMode', ...
      'Unsupported dynamic objective mode: %s.', char(config.dynamicObjectiveMode));
end
end

function mode = localNormalizeDynamicRouteMode(mode)
%LOCALNORMALIZEDYNAMICROUTEMODE Normalize the estimator route tag.

mode = lower(strtrim(string(mode)));
switch mode
  case {"core-only", "coreonly", "core"}
    mode = "core-only";
  case {"robust-estimator", "robustestimator", "robust", "default"}
    mode = "robust-estimator";
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidDynamicRouteMode', ...
      'Unsupported dynamicRouteMode: %s.', char(mode));
end
end

function dynOpt = localApplyDynamicRouteMode(dynOpt, routeMode)
%LOCALAPPLYDYNAMICROUTEMODE Apply one-shot core or robust branch route.

routeMode = localNormalizeDynamicRouteMode(routeMode);
switch string(routeMode)
  case "core-only"
    dynOpt.disableDoaBasinEntry = true;
    dynOpt.doaBasinEntryScope = 'off';
    dynOpt.disableKnownEmbedded = true;
    dynOpt.disableUnknownWarmAnchor = true;
  case "robust-estimator"
    dynOpt.disableDoaBasinEntry = false;
    dynOpt.disableKnownEmbedded = false;
    dynOpt.disableUnknownWarmAnchor = false;
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:InvalidDynamicRouteMode', ...
      'Unsupported dynamic route mode: %s.', char(routeMode));
end
end


function method = localBuildDynamicMethod(methodName, doaCrbScale, defaultRouteMode)
%LOCALBUILDDYNAMICMETHOD Build one SS/MS MF method descriptor.

if nargin < 3 || isempty(defaultRouteMode)
  defaultRouteMode = "core-only";
end
methodName = string(methodName);
baseMethodName = localBaseDynamicMethodName(methodName);
routeMode = localNormalizeDynamicRouteMode(defaultRouteMode);
if endsWith(methodName, "-Robust")
  routeMode = "robust-estimator";
elseif endsWith(methodName, "-Core")
  routeMode = "core-only";
end

if contains(baseMethodName, "-K")
  fdRateMode = "known";
  isKnownRate = true;
elseif contains(baseMethodName, "-U")
  fdRateMode = "unknown";
  isKnownRate = false;
else
  error('scanMfMsMleCrbCleanBoundScale:UnknownMethod', 'Unsupported method name "%s".', char(methodName));
end
if contains(baseMethodName, "-IP-")
  phaseMode = "independent";
elseif contains(baseMethodName, "-CP-")
  phaseMode = "continuous";
else
  error('scanMfMsMleCrbCleanBoundScale:UnknownPhaseMethod', 'Unsupported method phase in "%s".', char(methodName));
end
method = struct();
method.displayName = methodName;
method.baseDisplayName = baseMethodName;
method.dynamicRouteMode = routeMode;
if startsWith(baseMethodName, "MS-")
  method.satMode = "multi";
else
  method.satMode = "single";
end
method.phaseMode = phaseMode;
method.fdRateMode = fdRateMode;
method.isKnownRate = logical(isKnownRate);
method.doaHalfWidthDeg = reshape(double(doaCrbScale), [], 1);
end

function baseMethodName = localBaseDynamicMethodName(methodName)
%LOCALBASEDYNAMICMETHODNAME Strip replay-only route suffix from a method name.

baseMethodName = string(methodName);
baseMethodName = erase(baseMethodName, "-Robust");
baseMethodName = erase(baseMethodName, "-Core");
end

function staticSeedCase = localSelectStaticSeedCase(method, staticBundle)
%LOCALSELECTSTATICSEEDCASE Select the matching static seed for SS/MS MF.

if string(method.satMode) == "multi"
  if ~isfield(staticBundle, 'caseStaticMs') || isempty(staticBundle.caseStaticMs)
    error('scanMfMsMleCrbCleanBoundScale:MissingStaticMsSeed', ...
      'A multi-satellite static seed is required by method "%s".', char(method.displayName));
  end
  staticSeedCase = staticBundle.caseStaticMs;
else
  if ~isfield(staticBundle, 'caseStaticRefOnly') || isempty(staticBundle.caseStaticRefOnly)
    error('scanMfMsMleCrbCleanBoundScale:MissingStaticRefSeed', ...
      'A reference-satellite static seed is required by method "%s".', char(method.displayName));
  end
  staticSeedCase = staticBundle.caseStaticRefOnly;
end
end

function [initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, staticSeedCase, truth, config)
%LOCALBUILDDYNAMICINITOVERRIDE Build a static seed with a fixed CRB-scaled MF box.

initMode = "truth-crb-scaled-static-seed-fixed-bound";
initParamOverride = [];
initRequestMeta = struct('initMode', initMode, 'initParam', [], ...
  'initDoaParam', [], 'initDoaHalfWidth', [], 'seedDoaParam', [], ...
  'boundCenterDoaParam', []);

truthDoaParam = reshape(truth.latlonTrueDeg, [], 1);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundScale:MissingTruthDoa', ...
    'Truth DoA is required to build the fixed clean MF search box.');
end
truthDoaParam = truthDoaParam(1:2);

estResult = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
seedDoaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(seedDoaParam) < 2 || any(~isfinite(seedDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundScale:MissingStaticDoaSeed', ...
    'A finite truth-centered static DoA seed is required for the MF refinement.');
end
seedDoaParam = seedDoaParam(1:2);
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});

% The override is the solver seed. The DoA bound center below intentionally
% stays at truth so that static seed errors cannot expand the clean envelope.
if method.isKnownRate
  initParamOverride = [seedDoaParam(:); truthFdRefHz];
else
  initParamOverride = [seedDoaParam(:); truthFdRefHz; truthFdRateHzPerSec];
end
initRequestMeta.initParam = initParamOverride(:);
initRequestMeta.initDoaParam = truthDoaParam(:);
initRequestMeta.initDoaHalfWidth = method.doaHalfWidthDeg(:);
initRequestMeta.seedDoaParam = seedDoaParam(:);
initRequestMeta.boundCenterDoaParam = truthDoaParam(:);
end

function initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta)
%LOCALBUILDDYNAMICINITMETA Extract requested and internal initializer.

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
initMeta = struct();
initMeta.requestedInitMode = string(localGetFieldOrDefault(initRequestMeta, 'initMode', ""));
initMeta.requestedInitParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initParam', []), [], 1);
initMeta.requestedInitDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaParam', []), [], 1);
initMeta.seedDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'seedDoaParam', []), [], 1);
initMeta.boundCenterDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'boundCenterDoaParam', []), [], 1);
initMeta.initParam = initParam(:);
initMeta.initDoaParam = NaN(2, 1);
initMeta.initFdRefHz = NaN;
initMeta.initFdRateHzPerSec = NaN;
if numel(initParam) >= 2
  initMeta.initDoaParam = initParam(1:2);
end
if numel(initParam) >= 3
  initMeta.initFdRefHz = initParam(3);
end
if ~method.isKnownRate && numel(initParam) >= 4
  initMeta.initFdRateHzPerSec = initParam(4);
end
initMeta.dynamicRouteMode = string(localGetFieldOrDefault(method, 'dynamicRouteMode', ""));
initMeta.baseDisplayName = string(localGetFieldOrDefault(method, 'baseDisplayName', method.displayName));
initMeta.searchHalfWidthDeg = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaHalfWidth', NaN(2, 1)), [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
end

function row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
  fdRangeUse, fdRateRangeUse, crbBundle, config, initMeta, rangeCfg)
%LOCALBUILDCASEMETRICROW Convert one clean case into a compact metric row.

row = localEmptyCaseRow();
methodName = string(methodName);
info = summarizeDoaDopplerCase(caseUse, truth);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
initEval = localGetFieldOrDefault(estDebug, 'initEval', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());

row.displayName = methodName;
row.dynamicObjectiveMode = string(config.dynamicObjectiveMode);
row.dynamicRouteMode = string(localGetFieldOrDefault(initMeta, 'dynamicRouteMode', config.dynamicRouteMode));
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.doaCrbScale = localGetFieldOrDefault(rangeCfg, 'doaCrbScale', NaN);
row.fdRefHalfToothFraction = localGetFieldOrDefault(rangeCfg, 'fdRefHalfToothFraction', NaN);
row.fdRefCrbScale = localGetFieldOrDefault(rangeCfg, 'fdRefCrbScale', NaN);
row.fdRateHalfWidthHzPerSec = localGetFieldOrDefault(rangeCfg, 'fdRateHalfWidthHzPerSec', NaN);
row.numFrame = fixture.sceneSeq.numFrame;
row.frameIntvlSec = median(diff(fixture.sceneSeq.timeOffsetSec));
row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
row.toothStepHz = toothStepHz;
row.truthLatDeg = truth.latlonTrueDeg(1);
row.truthLonDeg = truth.latlonTrueDeg(2);
row.finalLatDeg = info.latEstDeg;
row.finalLonDeg = info.lonEstDeg;
row.latErrDeg = info.latErrDeg;
row.lonErrDeg = info.lonErrDeg;
row.sphericalAngleErrDeg = info.angleErrDeg;
row.finalAngleErrDeg = info.angleErrDeg;
row.fdRefTrueHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
row.fdRateTrueHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
row.finalFdRefHz = info.fdRefEstHz;
row.fdRefErrHz = info.fdRefErrHz;
row.fdRefAbsErrHz = abs(info.fdRefErrHz);
row.finalFdRateHzPerSec = info.fdRateEstHzPerSec;
row.fdRateErrHzPerSec = info.fdRateErrHzPerSec;
row.fdRateAbsErrHzPerSec = abs(info.fdRateErrHzPerSec);
row.fdLineRmseHz = info.fdLineRmseHz;
row.solverResolved = logical(info.isResolved);
row.exitflag = info.exitflag;
row.iterations = info.iterations;
row.funcCount = info.funcCount;
row.runTimeMs = info.runTimeMs;

if endsWith(methodName, "SF-Static")
  row.satMode = string(ternary(startsWith(methodName, "MS-"), "multi", "single"));
  row.frameMode = "single";
  row.phaseMode = "single-frame";
  row.fdRateMode = "zero";
  row.initLatDeg = initMeta.initDoaParam(1);
  row.initLonDeg = initMeta.initDoaParam(2);
  row.initAngleErrDeg = calcLatlonAngleError(initMeta.initDoaParam(:), truth.latlonTrueDeg(:));
  row.finalMinusInitAngleDeg = calcLatlonAngleError([row.finalLatDeg; row.finalLonDeg], initMeta.initDoaParam(:));
  row.initFdRefHz = initMeta.initFdRefHz;
  row.initFdRateHzPerSec = NaN;
  row.searchHalfWidthLatDeg = initMeta.searchHalfWidthDeg(1);
  row.searchHalfWidthLonDeg = initMeta.searchHalfWidthDeg(min(2, numel(initMeta.searchHalfWidthDeg)));
  row.fdRateRangeLowHzPerSec = NaN;
  row.fdRateRangeHighHzPerSec = NaN;
  row.mfInitMode = initMeta.requestedInitMode;
else
  baseMethodName = localBaseDynamicMethodName(methodName);
  row.satMode = string(ternary(startsWith(baseMethodName, "MS-"), "multi", "single"));
  row.frameMode = "multi";
  if contains(baseMethodName, "-IP-")
    row.phaseMode = "independent";
  else
    row.phaseMode = "continuous";
  end
  row.fdRateMode = string(ternary(endsWith(baseMethodName, "-K"), "known", "unknown"));
  row.initLatDeg = initMeta.initDoaParam(1);
  row.initLonDeg = initMeta.initDoaParam(2);
  row.initAngleErrDeg = calcLatlonAngleError(initMeta.initDoaParam(:), truth.latlonTrueDeg(:));
  row.finalMinusInitAngleDeg = calcLatlonAngleError([row.finalLatDeg; row.finalLonDeg], initMeta.initDoaParam(:));
  row.initFdRefHz = initMeta.initFdRefHz;
  row.initFdRateHzPerSec = initMeta.initFdRateHzPerSec;
  row.searchHalfWidthLatDeg = initMeta.searchHalfWidthDeg(1);
  row.searchHalfWidthLonDeg = initMeta.searchHalfWidthDeg(min(2, numel(initMeta.searchHalfWidthDeg)));
  row.mfInitMode = initMeta.requestedInitMode;
end

row.seedLatDeg = NaN;
row.seedLonDeg = NaN;
seedDoaParam = reshape(localGetFieldOrDefault(initMeta, 'seedDoaParam', []), [], 1);
if numel(seedDoaParam) >= 2 && all(isfinite(seedDoaParam(1:2)))
  row.seedLatDeg = seedDoaParam(1);
  row.seedLonDeg = seedDoaParam(2);
end
row.boundCenterLatDeg = row.truthLatDeg;
row.boundCenterLonDeg = row.truthLonDeg;
boundCenterDoaParam = reshape(localGetFieldOrDefault(initMeta, 'boundCenterDoaParam', []), [], 1);
if numel(boundCenterDoaParam) >= 2 && all(isfinite(boundCenterDoaParam(1:2)))
  row.boundCenterLatDeg = boundCenterDoaParam(1);
  row.boundCenterLonDeg = boundCenterDoaParam(2);
end
boxTolDeg = 1e-10;
row.doaLatBoxViolation = abs(row.finalLatDeg - row.boundCenterLatDeg) > row.searchHalfWidthLatDeg + boxTolDeg;
row.doaLonBoxViolation = abs(row.finalLonDeg - row.boundCenterLonDeg) > row.searchHalfWidthLonDeg + boxTolDeg;
row.doaTruthBoxViolation = row.doaLatBoxViolation || row.doaLonBoxViolation;

[crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
[angleCrbStdDeg, angleAux] = projectCrbToAngleMetric(crbDoa, truth.latlonTrueDeg, 'latlon');
row.crbLatStdDeg = sqrt(max(real(crbDoa(1, 1)), 0));
row.crbLonStdDeg = sqrt(max(real(crbDoa(2, 2)), 0));
row.crbTraceStdDeg = sqrt(max(real(trace(crbDoa)), 0));
row.crbSphericalApproxStdDeg = angleCrbStdDeg;
row.crbAngleMetricVarRad2 = angleAux.angleVarRad2;
row.crbMetricTraceRatio = localSafeRatio(row.crbTraceStdDeg, row.crbSphericalApproxStdDeg);
row.angleErrOverSphericalCrb = localSafeRatio(row.sphericalAngleErrDeg, row.crbSphericalApproxStdDeg);
row.angleErrOverTraceCrb = localSafeRatio(row.sphericalAngleErrDeg, row.crbTraceStdDeg);
row.fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
row.fdRefErrOverCrb = localSafeRatio(row.fdRefAbsErrHz, row.fdRefCrbStdHz);

row.fdRangeLowHz = fdRangeUse(1);
row.fdRangeHighHz = fdRangeUse(2);
row.fdRangeHalfWidthHz = 0.5 * abs(fdRangeUse(2) - fdRangeUse(1));
if row.frameMode == "multi" && row.fdRateMode == "unknown"
  row.fdRateRangeLowHzPerSec = fdRateRangeUse(1);
  row.fdRateRangeHighHzPerSec = fdRateRangeUse(2);
end
row.fdRefBoundaryHit = localIsNearRangeBoundary(row.finalFdRefHz, fdRangeUse, config.boundaryTolFraction);
row.fdRateBoundaryHit = row.frameMode == "multi" && row.fdRateMode == "unknown" && ...
  localIsNearRangeBoundary(row.finalFdRateHzPerSec, fdRateRangeUse, config.boundaryTolFraction);

row.initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objectiveImprove = row.initObj - row.finalObj;
row.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
[row.refCoherence, row.nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval);
row.stepSize = localGetFieldOrDefault(optimInfo, 'stepsize', NaN);
row.firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
row.constrViolation = localGetFieldOrDefault(optimInfo, 'constrviolation', NaN);
row.algorithm = string(localGetFieldOrDefault(optimInfo, 'algorithm', ""));
row.usedScaledSolve = logical(localGetFieldOrDefault(optimInfo, 'usedScaledSolve', false));
row.fallbackTriggered = logical(localGetFieldOrDefault(optimInfo, 'fallbackTriggered', false));
row.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));
row.candidateCount = numel(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.candidateVariantText = strjoin(cellstr(reshape(string(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1))), 1, [])), ',');
row.fdRefInitMoveHz = row.finalFdRefHz - row.initFdRefHz;
row.fdRateInitMoveHzPerSec = row.finalFdRateHzPerSec - row.initFdRateHzPerSec;
end

function caseTable = localBuildCaseTable(repeatCell)
%LOCALBUILDCASETABLE Collect all per-method rows into one table.

rowList = repmat(localEmptyCaseRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isempty(repeatOut)
    continue;
  end
  rowList = [rowList; repeatOut.caseRowList(:)]; %#ok<AGROW>
end
caseTable = struct2table(rowList(:));
end

function caseTable = localAnnotateHealthFlags(caseTable, config)
%LOCALANNOTATEHEALTHFLAGS Add clean health and joint-trim flags.

if isempty(caseTable) || height(caseTable) == 0
  return;
end
numRow = height(caseTable);
isDynamic = contains(caseTable.displayName, "-MF-");
isUnknown = endsWith(caseTable.displayName, "-U") & isDynamic;

solverHealthy = true(numRow, 1);
solverHealthy(isDynamic) = isfinite(caseTable.finalObj(isDynamic)) & ...
  isfinite(caseTable.iterations(isDynamic)) & caseTable.iterations(isDynamic) > 0 & ...
  isfinite(caseTable.objectiveImprove(isDynamic));

freqHealthy = true(numRow, 1);
fdRefThresholdHz = config.healthFdRefToothFraction .* caseTable.toothStepHz;
freqHealthy(isDynamic) = caseTable.fdRefAbsErrHz(isDynamic) <= fdRefThresholdHz(isDynamic);

rateHealthy = true(numRow, 1);
rateHealthy(isUnknown) = caseTable.fdRateAbsErrHzPerSec(isUnknown) <= config.healthFdRateAbsMaxHzPerSec;

boundaryHealthy = true(numRow, 1);
boundaryHealthy(isDynamic) = ~(caseTable.fdRefBoundaryHit(isDynamic) | caseTable.fdRateBoundaryHit(isDynamic));

notNoSolve = true(numRow, 1);
noSolve = isDynamic & (caseTable.iterations <= 0 | ~isfinite(caseTable.iterations)) & ...
  (~isfinite(caseTable.objectiveImprove) | abs(caseTable.objectiveImprove) <= eps);
notNoSolve(noSolve) = false;

healthResolved = solverHealthy & freqHealthy & rateHealthy & boundaryHealthy & notNoSolve;
angleNormOk = ~isDynamic | caseTable.angleErrOverSphericalCrb <= config.trimAngleNormMax;
fdRefNormOk = ~isDynamic | caseTable.fdRefErrOverCrb <= config.trimFdRefNormMax;
trimKeep = healthResolved & angleNormOk;
jointTrimKeep = healthResolved & angleNormOk & fdRefNormOk;
fdRefNormTail = isDynamic & healthResolved & angleNormOk & ~fdRefNormOk;

caseTable.isDynamicMethod = isDynamic;
caseTable.solverHealthy = solverHealthy;
caseTable.freqHealthy = freqHealthy;
caseTable.rateHealthy = rateHealthy;
caseTable.boundaryHealthy = boundaryHealthy;
caseTable.notNoSolve = notNoSolve;
caseTable.healthResolved = healthResolved;
caseTable.trimKeep = trimKeep;
caseTable.jointTrimKeep = jointTrimKeep;
caseTable.fdRefNormTail = fdRefNormTail;
caseTable.rejectReason = localBuildRejectReason(caseTable, angleNormOk, fdRefNormOk);
end


function aggregateTable = localBuildScaleAggregateTable(caseTable, filterName, keepMask)
%LOCALBUILDSCALEAGGREGATETABLE Aggregate clean metrics by method and scale.

if isempty(caseTable) || height(caseTable) == 0
  aggregateTable = table();
  return;
end
keepMask = reshape(logical(keepMask), [], 1);
if numel(keepMask) ~= height(caseTable)
  error('scanMfMsMleCrbCleanBoundScale:InvalidKeepMask', ...
    'Scale aggregate keep mask must match caseTable height.');
end
[groupId, displayName, snrDb, doaCrbScale, fdRefCrbScale] = findgroups( ...
  caseTable.displayName, caseTable.snrDb, caseTable.doaCrbScale, caseTable.fdRefCrbScale);
numGroup = max(groupId);
rowList = repmat(localEmptyScaleAggregateRow(), numGroup, 1);
for iGroup = 1:numGroup
  rawMask = groupId == iGroup;
  mask = rawMask & keepMask;
  row = localEmptyScaleAggregateRow();
  row.filterName = string(filterName);
  row.displayName = displayName(iGroup);
  row.snrDb = snrDb(iGroup);
  row.doaCrbScale = doaCrbScale(iGroup);
  row.fdRefCrbScale = fdRefCrbScale(iGroup);
  row.rawCount = nnz(rawMask);
  row.keepCount = nnz(mask);
  row.keepRate = localSafeRatio(row.keepCount, row.rawCount);
  row.angleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(mask));
  row.angleMseOverCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(mask), ...
    caseTable.crbSphericalApproxStdDeg(mask));
  row.angleRmseOverCrb = sqrt(row.angleMseOverCrb);
  row.fdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(mask));
  row.fdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(mask), ...
    caseTable.fdRefCrbStdHz(mask));
  row.fdRefRmseOverCrb = sqrt(row.fdRefMseOverCrb);
  row.fdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(mask));
  row.boundaryHitRate = mean(double(caseTable.fdRefBoundaryHit(mask) | caseTable.fdRateBoundaryHit(mask)), 'omitnan');
  row.doaTruthBoxViolationRate = mean(double(caseTable.doaTruthBoxViolation(mask)), 'omitnan');
  row.iterationsMedian = median(caseTable.iterations(mask), 'omitnan');
  row.firstOrderOptMedian = median(caseTable.firstOrderOpt(mask), 'omitnan');
  row.topRejectReason = localTopRejectReason(caseTable.rejectReason(rawMask & ~keepMask));
  rowList(iGroup) = row;
end
aggregateTable = struct2table(rowList(:));
aggregateTable = sortrows(aggregateTable, {'displayName', 'snrDb', 'doaCrbScale', 'fdRefCrbScale'});
end

function row = localEmptyScaleAggregateRow()
%LOCALEMTYSCALEAGGREGATEROW Return one typed scale aggregate row.

row = struct('filterName', "", 'displayName', "", 'snrDb', NaN, ...
  'doaCrbScale', NaN, 'fdRefCrbScale', NaN, 'rawCount', 0, ...
  'keepCount', 0, 'keepRate', NaN, 'angleRmseDeg', NaN, ...
  'angleRmseOverCrb', NaN, 'angleMseOverCrb', NaN, ...
  'fdRefRmseHz', NaN, 'fdRefRmseOverCrb', NaN, 'fdRefMseOverCrb', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'boundaryHitRate', NaN, ...
  'doaTruthBoxViolationRate', NaN, 'iterationsMedian', NaN, ...
  'firstOrderOptMedian', NaN, 'topRejectReason', "");
end

function summaryTable = localBuildScaleReadableRmseCrbTable(rawTable, trimTable)
%LOCALBUILDSCALEREADABLERMSCRBTABLE Build final human-readable scale table.

if isempty(rawTable) || height(rawTable) == 0
  summaryTable = table();
  return;
end
rowTemplate = struct('displayName', "", 'snrDb', NaN, 'doaCrbScale', NaN, ...
  'fdRefCrbScale', NaN, 'rawCount', 0, 'trimCount', 0, ...
  'trimKeepRate', NaN, 'trimStatus', "", ...
  'rawDoaRmseOverCrb', NaN, 'trimDoaRmseOverCrb', NaN, ...
  'rawFdRefRmseOverCrb', NaN, 'trimFdRefRmseOverCrb', NaN, ...
  'rawDoaRmseDeg', NaN, 'trimDoaRmseDeg', NaN, ...
  'rawFdRefRmseHz', NaN, 'trimFdRefRmseHz', NaN, 'topTrimRejectReason', "");
rowList = repmat(rowTemplate, height(rawTable), 1);
for iRow = 1:height(rawTable)
  rawRow = rawTable(iRow, :);
  row = rowTemplate;
  row.displayName = rawRow.displayName;
  row.snrDb = rawRow.snrDb;
  row.doaCrbScale = rawRow.doaCrbScale;
  row.fdRefCrbScale = rawRow.fdRefCrbScale;
  row.rawCount = rawRow.rawCount;
  row.rawDoaRmseOverCrb = rawRow.angleRmseOverCrb;
  row.rawFdRefRmseOverCrb = rawRow.fdRefRmseOverCrb;
  row.rawDoaRmseDeg = rawRow.angleRmseDeg;
  row.rawFdRefRmseHz = rawRow.fdRefRmseHz;
  trimRow = localFindScaleAggregateRow(trimTable, rawRow);
  if ~isempty(trimRow)
    row.trimCount = trimRow.keepCount;
    row.trimKeepRate = trimRow.keepRate;
    row.trimDoaRmseOverCrb = trimRow.angleRmseOverCrb;
    row.trimFdRefRmseOverCrb = trimRow.fdRefRmseOverCrb;
    row.trimDoaRmseDeg = trimRow.angleRmseDeg;
    row.trimFdRefRmseHz = trimRow.fdRefRmseHz;
    row.topTrimRejectReason = trimRow.topRejectReason;
  end
  row.trimStatus = localClassifyScaleTrimStatus(row);
  rowList(iRow) = row;
end
summaryTable = struct2table(rowList(:));
summaryTable = sortrows(summaryTable, {'displayName', 'snrDb', 'doaCrbScale', 'fdRefCrbScale'});
end

function trimRow = localFindScaleAggregateRow(tableIn, rawRow)
%LOCALFINDSCALEAGGREGATEROW Find matching scale aggregate row.

trimRow = table();
if isempty(tableIn) || height(tableIn) == 0
  return;
end
mask = tableIn.displayName == rawRow.displayName & tableIn.snrDb == rawRow.snrDb & ...
  tableIn.doaCrbScale == rawRow.doaCrbScale & tableIn.fdRefCrbScale == rawRow.fdRefCrbScale;
if any(mask)
  trimRow = tableIn(find(mask, 1, 'first'), :);
end
end

function status = localClassifyScaleTrimStatus(row)
%LOCALCLASSIFYSCALETRIMSTATUS Classify one raw-vs-trim row.

if row.rawCount <= 0
  status = "empty";
elseif row.trimCount <= 0
  status = "all-rejected";
elseif row.trimKeepRate < 0.8
  status = "partial-trim";
elseif isfinite(row.trimDoaRmseOverCrb) && isfinite(row.trimFdRefRmseOverCrb) && ...
    row.trimDoaRmseOverCrb <= 1.1 && row.trimFdRefRmseOverCrb <= 1.1
  status = "near-crb";
else
  status = "kept";
end
end

function surfaceTable = localBuildScaleSurfaceSummaryTable(trimTable)
%LOCALBUILDSCALESURFACESUMMARYTABLE Pick the best trim point per method/SNR.

if isempty(trimTable) || height(trimTable) == 0
  surfaceTable = table();
  return;
end
[groupId, displayName, snrDb] = findgroups(trimTable.displayName, trimTable.snrDb);
numGroup = max(groupId);
rowTemplate = struct('displayName', "", 'snrDb', NaN, 'numGridPoint', 0, ...
  'bestDoaCrbScale', NaN, 'bestFdRefCrbScale', NaN, ...
  'bestAngleRmseOverCrb', NaN, 'bestFdRefRmseOverCrb', NaN, ...
  'bestKeepRate', NaN, 'targetPass', false);
rowList = repmat(rowTemplate, numGroup, 1);
for iGroup = 1:numGroup
  mask = groupId == iGroup;
  sub = trimTable(mask, :);
  score = sub.angleRmseOverCrb + 0.25 * sub.fdRefRmseOverCrb;
  score(~isfinite(score)) = Inf;
  [~, idx] = min(score);
  row = rowTemplate;
  row.displayName = displayName(iGroup);
  row.snrDb = snrDb(iGroup);
  row.numGridPoint = height(sub);
  if isfinite(score(idx))
    row.bestDoaCrbScale = sub.doaCrbScale(idx);
    row.bestFdRefCrbScale = sub.fdRefCrbScale(idx);
    row.bestAngleRmseOverCrb = sub.angleRmseOverCrb(idx);
    row.bestFdRefRmseOverCrb = sub.fdRefRmseOverCrb(idx);
    row.bestKeepRate = sub.keepRate(idx);
    row.targetPass = row.bestKeepRate >= 0.8 && row.bestAngleRmseOverCrb <= 1.1 && ...
      (isnan(row.bestFdRefRmseOverCrb) || row.bestFdRefRmseOverCrb <= 1.1);
  end
  rowList(iGroup) = row;
end
surfaceTable = struct2table(rowList(:));
end

function auditTable = localBuildScaleSearchRangeAuditTable(caseTable)
%LOCALBUILDSCALESEARCHRANGEAUDITTABLE Summarize actual CRB-scaled range coverage.

if isempty(caseTable) || height(caseTable) == 0
  auditTable = table();
  return;
end
[groupId, displayName, snrDb, doaCrbScale, fdRefCrbScale] = findgroups( ...
  caseTable.displayName, caseTable.snrDb, caseTable.doaCrbScale, caseTable.fdRefCrbScale);
numGroup = max(groupId);
rowTemplate = struct('displayName', "", 'snrDb', NaN, 'doaCrbScale', NaN, ...
  'fdRefCrbScale', NaN, 'numRepeat', 0, 'minDoaLatHalfOverCrb', NaN, ...
  'minDoaLonHalfOverCrb', NaN, 'minFdRefHalfOverCrb', NaN, ...
  'doaBelowCrb', false, 'fdRefBelowCrb', false, ...
  'doaTruthBoxViolationRate', NaN, 'fdRefBoundaryHitRate', NaN);
rowList = repmat(rowTemplate, numGroup, 1);
for iGroup = 1:numGroup
  mask = groupId == iGroup;
  row = rowTemplate;
  row.displayName = displayName(iGroup);
  row.snrDb = snrDb(iGroup);
  row.doaCrbScale = doaCrbScale(iGroup);
  row.fdRefCrbScale = fdRefCrbScale(iGroup);
  row.numRepeat = nnz(mask);
  row.minDoaLatHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLatDeg(mask), caseTable.crbLatStdDeg(mask)), [], 'omitnan');
  row.minDoaLonHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLonDeg(mask), caseTable.crbLonStdDeg(mask)), [], 'omitnan');
  row.minFdRefHalfOverCrb = min(localSafeVectorRatio(caseTable.fdRangeHalfWidthHz(mask), caseTable.fdRefCrbStdHz(mask)), [], 'omitnan');
  row.doaBelowCrb = row.minDoaLatHalfOverCrb < 1 || row.minDoaLonHalfOverCrb < 1;
  row.fdRefBelowCrb = row.minFdRefHalfOverCrb < 1;
  row.doaTruthBoxViolationRate = mean(double(caseTable.doaTruthBoxViolation(mask)), 'omitnan');
  row.fdRefBoundaryHitRate = mean(double(caseTable.fdRefBoundaryHit(mask)), 'omitnan');
  rowList(iGroup) = row;
end
auditTable = struct2table(rowList(:));
auditTable = sortrows(auditTable, {'displayName', 'snrDb', 'doaCrbScale', 'fdRefCrbScale'});
end

function outlierTable = localBuildScaleOutlierTable(caseTable)
%LOCALBUILDSCALEOUTLIERTABLE Keep rejected dynamic rows with scale coordinates.

if isempty(caseTable) || height(caseTable) == 0
  outlierTable = table();
  return;
end
mask = contains(caseTable.displayName, "-MF-") & ~caseTable.jointTrimKeep;
columnList = {'displayName', 'snrDb', 'taskSeed', 'doaCrbScale', 'fdRefCrbScale', ...
  'sphericalAngleErrDeg', 'angleErrOverSphericalCrb', 'fdRefAbsErrHz', ...
  'fdRefErrOverCrb', 'fdRateAbsErrHzPerSec', 'iterations', 'firstOrderOpt', ...
  'objectiveImprove', 'fdRefBoundaryHit', 'fdRateBoundaryHit', 'healthResolved', ...
  'jointTrimKeep', 'rejectReason'};
columnList = columnList(ismember(columnList, caseTable.Properties.VariableNames));
outlierTable = caseTable(mask, columnList);
if height(outlierTable) > 0
  outlierTable = sortrows(outlierTable, {'angleErrOverSphericalCrb'}, {'descend'});
end
end

function reason = localTopRejectReason(reasonVec)
%LOCALTOPREJECTREASON Return the most frequent nonempty reject reason.

reasonVec = reshape(string(reasonVec), [], 1);
reasonVec = reasonVec(strlength(reasonVec) > 0 & reasonVec ~= "kept" & reasonVec ~= "baseline");
if isempty(reasonVec)
  reason = "";
  return;
end
[grp, reasonName] = findgroups(reasonVec);
count = splitapply(@numel, reasonVec, grp);
[~, idx] = max(count);
reason = reasonName(idx);
end

function rejectReason = localBuildRejectReason(caseTable, angleNormOk, fdRefNormOk)
%LOCALBUILDREJECTREASON Build compact health / trim reason tags.

numRow = height(caseTable);
rejectReason = strings(numRow, 1);
for iRow = 1:numRow
  if ~caseTable.isDynamicMethod(iRow)
    rejectReason(iRow) = "baseline";
    continue;
  end
  reasonList = strings(0, 1);
  if ~caseTable.solverHealthy(iRow)
    reasonList(end + 1, 1) = "solver"; %#ok<AGROW>
  end
  if ~caseTable.freqHealthy(iRow)
    reasonList(end + 1, 1) = "fdRef"; %#ok<AGROW>
  end
  if ~caseTable.rateHealthy(iRow)
    reasonList(end + 1, 1) = "fdRate"; %#ok<AGROW>
  end
  if ~caseTable.boundaryHealthy(iRow)
    reasonList(end + 1, 1) = "boundary"; %#ok<AGROW>
  end
  if ~caseTable.notNoSolve(iRow)
    reasonList(end + 1, 1) = "no-solve"; %#ok<AGROW>
  end
  if caseTable.healthResolved(iRow) && ~angleNormOk(iRow)
    reasonList(end + 1, 1) = "trim-angle"; %#ok<AGROW>
  end
  if caseTable.healthResolved(iRow) && ~fdRefNormOk(iRow)
    reasonList(end + 1, 1) = "trim-fdRef"; %#ok<AGROW>
  end
  if isempty(reasonList)
    rejectReason(iRow) = "kept";
  else
    rejectReason(iRow) = strjoin(cellstr(reasonList), "+");
  end
end
end

function compareTable = localBuildCleanSsMsCompareTable(caseTable, config)
%LOCALBUILDCLEANSSMSCOMPARETABLE Compare SS-MF and MS-MF under each filter.

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'phaseMode', "", ...
  'fdRateMode', "", 'ssName', "", 'msName', "", ...
  'ssKeepRate', NaN, 'msKeepRate', NaN, ...
  'ssAngleRmseDeg', NaN, 'msAngleRmseDeg', NaN, 'msAbsGainOverSs', NaN, ...
  'ssAngleMseOverOwnCrb', NaN, 'msAngleMseOverOwnCrb', NaN, 'msAngleMseOverSsCrb', NaN, ...
  'ssFdRefMseOverOwnCrb', NaN, 'msFdRefMseOverOwnCrb', NaN, ...
  'ssFdRateRmseHzPerSec', NaN, 'msFdRateRmseHzPerSec', NaN, ...
  'msBeatsSsAbs', false, 'msOwnCrbTargetPass', false, 'targetClass', "");
filterNameList = ["raw"; "health"; "joint-trim"];
snrList = unique(caseTable.snrDb, 'stable');
phaseTagList = ["CP"; "IP"];
fdRateTagList = ["K"; "U"];
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList) * ...
  numel(phaseTagList) * numel(fdRateTagList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iPhase = 1:numel(phaseTagList)
    phaseTag = phaseTagList(iPhase);
    phaseMode = string(ternary(phaseTag == "IP", "independent", "continuous"));
    for iFdRate = 1:numel(fdRateTagList)
      fdRateTag = fdRateTagList(iFdRate);
      fdRateMode = string(ternary(fdRateTag == "K", "known", "unknown"));
      ssName = "SS-MF-" + phaseTag + "-" + fdRateTag;
      msName = "MS-MF-" + phaseTag + "-" + fdRateTag;
      ssRaw = caseTable.displayName == ssName & caseTable.snrDb == snrDb;
      msRaw = caseTable.displayName == msName & caseTable.snrDb == snrDb;
      for iFilter = 1:numel(filterNameList)
        filterName = filterNameList(iFilter);
        idx = idx + 1;
        row = rowTemplate;
        row.filterName = filterName;
        row.snrDb = snrDb;
        row.phaseMode = phaseMode;
        row.fdRateMode = fdRateMode;
        row.ssName = ssName;
        row.msName = msName;
        ssMask = ssRaw & localFilterMask(caseTable, filterName);
        msMask = msRaw & localFilterMask(caseTable, filterName);
        row.ssKeepRate = localSafeRatio(nnz(ssMask), nnz(ssRaw));
        row.msKeepRate = localSafeRatio(nnz(msMask), nnz(msRaw));
        row.ssAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(ssMask));
        row.msAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(msMask));
        row.msAbsGainOverSs = localSafeRatio(row.ssAngleRmseDeg - row.msAngleRmseDeg, row.ssAngleRmseDeg);
        row.ssAngleMseOverOwnCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(ssMask), caseTable.crbSphericalApproxStdDeg(ssMask));
        row.msAngleMseOverOwnCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(msMask), caseTable.crbSphericalApproxStdDeg(msMask));
        row.msAngleMseOverSsCrb = localMsErrorOverMatchedSsCrb(caseTable, msMask, ssRaw);
        row.ssFdRefMseOverOwnCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(ssMask), caseTable.fdRefCrbStdHz(ssMask));
        row.msFdRefMseOverOwnCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(msMask), caseTable.fdRefCrbStdHz(msMask));
        row.ssFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(ssMask));
        row.msFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(msMask));
        row.msBeatsSsAbs = isfinite(row.msAbsGainOverSs) && row.msAbsGainOverSs > 0;
        row.msOwnCrbTargetPass = row.msAngleMseOverOwnCrb <= config.crbTargetMseRatio && ...
          row.msFdRefMseOverOwnCrb <= config.crbTargetMseRatio;
        if row.msOwnCrbTargetPass
          row.targetClass = "crb-level";
        elseif row.msBeatsSsAbs
          row.targetClass = "beats-ss-only";
        else
          row.targetClass = "not-yet";
        end
        rowList(idx) = row;
      end
    end
  end
end
compareTable = struct2table(rowList(:));
end

function compareTable = localBuildCleanKnownUnknownCompareTable(caseTable)
%LOCALBUILDCLEANKNOWNUNKNOWNCOMPARETABLE Compare K/U clean MF rows for CP and IP.

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'satMode', "", ...
  'phaseMode', "", 'knownName', "", 'unknownName', "", ...
  'knownKeepRate', NaN, 'unknownKeepRate', NaN, ...
  'knownAngleRmseDeg', NaN, 'unknownAngleRmseDeg', NaN, ...
  'unknownRmseGainOverKnown', NaN, ...
  'knownAngleMseOverCrb', NaN, 'unknownAngleMseOverCrb', NaN, ...
  'unknownMinusKnownAngleMseOverCrb', NaN, ...
  'knownFdRefMseOverCrb', NaN, 'unknownFdRefMseOverCrb', NaN, ...
  'unknownMinusKnownFdRefMseOverCrb', NaN, ...
  'knownFdRefRmseHz', NaN, 'unknownFdRefRmseHz', NaN, ...
  'knownFdRateRmseHzPerSec', NaN, 'unknownFdRateRmseHzPerSec', NaN, ...
  'comparisonClass', "");
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  compareTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
filterNameList = ["raw"; "health"; "joint-trim"];
snrList = unique(caseTable.snrDb, 'stable');
satPrefixList = ["SS"; "MS"];
phaseTagList = ["CP"; "IP"];
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList) * numel(satPrefixList) * numel(phaseTagList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iPhase = 1:numel(phaseTagList)
    phaseTag = phaseTagList(iPhase);
    for iSat = 1:numel(satPrefixList)
      satPrefix = satPrefixList(iSat);
      knownName = satPrefix + "-MF-" + phaseTag + "-K";
      unknownName = satPrefix + "-MF-" + phaseTag + "-U";
      knownRaw = caseTable.displayName == knownName & caseTable.snrDb == snrDb;
      unknownRaw = caseTable.displayName == unknownName & caseTable.snrDb == snrDb;
      for iFilter = 1:numel(filterNameList)
        filterName = filterNameList(iFilter);
        idx = idx + 1;
        row = rowTemplate;
        row.filterName = filterName;
        row.snrDb = snrDb;
        row.satMode = string(ternary(satPrefix == "MS", "multi", "single"));
        row.phaseMode = string(ternary(phaseTag == "IP", "independent", "continuous"));
        row.knownName = knownName;
        row.unknownName = unknownName;
        knownMask = knownRaw & localFilterMask(caseTable, filterName);
        unknownMask = unknownRaw & localFilterMask(caseTable, filterName);
        row.knownKeepRate = localSafeRatio(nnz(knownMask), nnz(knownRaw));
        row.unknownKeepRate = localSafeRatio(nnz(unknownMask), nnz(unknownRaw));
        row.knownAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(knownMask));
        row.unknownAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(unknownMask));
        row.unknownRmseGainOverKnown = localSafeRatio(row.knownAngleRmseDeg - row.unknownAngleRmseDeg, row.knownAngleRmseDeg);
        row.knownAngleMseOverCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(knownMask), caseTable.crbSphericalApproxStdDeg(knownMask));
        row.unknownAngleMseOverCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(unknownMask), caseTable.crbSphericalApproxStdDeg(unknownMask));
        row.unknownMinusKnownAngleMseOverCrb = row.unknownAngleMseOverCrb - row.knownAngleMseOverCrb;
        row.knownFdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(knownMask), caseTable.fdRefCrbStdHz(knownMask));
        row.unknownFdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(unknownMask), caseTable.fdRefCrbStdHz(unknownMask));
        row.unknownMinusKnownFdRefMseOverCrb = row.unknownFdRefMseOverCrb - row.knownFdRefMseOverCrb;
        row.knownFdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(knownMask));
        row.unknownFdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(unknownMask));
        row.knownFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(knownMask));
        row.unknownFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(unknownMask));
        row.comparisonClass = localClassifyKnownUnknownComparison(row);
        rowList(idx) = row;
      end
    end
  end
end
compareTable = struct2table(rowList(:));
end

function classText = localClassifyKnownUnknownComparison(row)
%LOCALCLASSIFYKNOWNUNKNOWNCOMPARISON Classify finite-MC CP-K/CP-U comparison.

if ~(isfinite(row.knownAngleRmseDeg) && isfinite(row.unknownAngleRmseDeg))
  classText = "missing";
  return;
end
relTol = 0.02;
if row.unknownAngleRmseDeg < (1 - relTol) * row.knownAngleRmseDeg
  classText = "unknown-smaller-rmse";
elseif row.knownAngleRmseDeg < (1 - relTol) * row.unknownAngleRmseDeg
  classText = "known-smaller-rmse";
else
  classText = "similar";
end
end

function curveTable = localBuildCleanSnrCurveTable(compareTable)
%LOCALBUILDCLEANSNRCURVETABLE Keep compact joint-trim SNR metrics for summary.

rowTemplate = struct('snrDb', NaN, 'phaseMode', "", 'fdRateMode', "", ...
  'ssName', "", 'msName', "", ...
  'ssAngleMseOverCrb', NaN, 'msAngleMseOverCrb', NaN, ...
  'ssFdRefMseOverCrb', NaN, 'msFdRefMseOverCrb', NaN, 'ssAngleRmseDeg', NaN, ...
  'msAngleRmseDeg', NaN, 'ssKeepRate', NaN, 'msKeepRate', NaN, ...
  'msAbsGainOverSs', NaN, 'targetClass', "");
if isempty(compareTable) || ~istable(compareTable) || height(compareTable) == 0
  curveTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
mask = compareTable.filterName == "joint-trim";
if ~any(mask)
  mask = compareTable.filterName == "raw";
end
sortNameList = {'snrDb'};
if ismember('phaseMode', compareTable.Properties.VariableNames)
  sortNameList{end + 1} = 'phaseMode';
end
if ismember('fdRateMode', compareTable.Properties.VariableNames)
  sortNameList{end + 1} = 'fdRateMode';
end
curveSource = sortrows(compareTable(mask, :), sortNameList);
rowList = repmat(rowTemplate, height(curveSource), 1);
for iRow = 1:height(curveSource)
  row = rowTemplate;
  row.snrDb = curveSource.snrDb(iRow);
  if ismember('phaseMode', curveSource.Properties.VariableNames)
    row.phaseMode = curveSource.phaseMode(iRow);
  end
  if ismember('fdRateMode', curveSource.Properties.VariableNames)
    row.fdRateMode = curveSource.fdRateMode(iRow);
  end
  if ismember('ssName', curveSource.Properties.VariableNames)
    row.ssName = curveSource.ssName(iRow);
  end
  if ismember('msName', curveSource.Properties.VariableNames)
    row.msName = curveSource.msName(iRow);
  end
  row.ssAngleMseOverCrb = curveSource.ssAngleMseOverOwnCrb(iRow);
  row.msAngleMseOverCrb = curveSource.msAngleMseOverOwnCrb(iRow);
  row.ssFdRefMseOverCrb = curveSource.ssFdRefMseOverOwnCrb(iRow);
  row.msFdRefMseOverCrb = curveSource.msFdRefMseOverOwnCrb(iRow);
  row.ssAngleRmseDeg = curveSource.ssAngleRmseDeg(iRow);
  row.msAngleRmseDeg = curveSource.msAngleRmseDeg(iRow);
  row.ssKeepRate = curveSource.ssKeepRate(iRow);
  row.msKeepRate = curveSource.msKeepRate(iRow);
  row.msAbsGainOverSs = curveSource.msAbsGainOverSs(iRow);
  row.targetClass = curveSource.targetClass(iRow);
  rowList(iRow) = row;
end
curveTable = struct2table(rowList(:));
end

function curveTable = localBuildCleanEstimateCrbCurveTable(trimAggregateTable, healthAggregateTable, rawAggregateTable)
%LOCALBUILDCLEANESTIMATECRBCURVETABLE Build absolute RMSE/CRB curves for plotting.

rowTemplate = struct('filterName', "", 'displayName', "", 'snrDb', NaN, ...
  'angleRmseDeg', NaN, 'angleCrbSphericalMedianDeg', NaN, ...
  'fdRefRmseHz', NaN, 'fdRefCrbMedianHz', NaN);
sourceTable = trimAggregateTable;
filterName = "joint-trim";
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  sourceTable = healthAggregateTable;
  filterName = "health";
end
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  sourceTable = rawAggregateTable;
  filterName = "raw";
end
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  curveTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
methodNameList = ["SS-SF-Static"; "MS-SF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "MS-MF-CP-K"; "MS-MF-CP-U"; "SS-MF-IP-K"; "SS-MF-IP-U"; "MS-MF-IP-K"; "MS-MF-IP-U"];
mask = ismember(sourceTable.displayName, methodNameList);
sourceTable = sortrows(sourceTable(mask, :), {'snrDb', 'displayName'});
rowList = repmat(rowTemplate, height(sourceTable), 1);
for iRow = 1:height(sourceTable)
  row = rowTemplate;
  row.filterName = filterName;
  row.displayName = string(sourceTable.displayName(iRow));
  row.snrDb = sourceTable.snrDb(iRow);
  row.angleRmseDeg = sourceTable.angleRmseDeg(iRow);
  row.angleCrbSphericalMedianDeg = sourceTable.angleCrbSphericalMedianDeg(iRow);
  row.fdRefRmseHz = sourceTable.fdRefRmseHz(iRow);
  row.fdRefCrbMedianHz = sourceTable.fdRefCrbMedianHz(iRow);
  rowList(iRow) = row;
end
curveTable = struct2table(rowList(:));
end

function plotData = localBuildCleanPlotData(cleanEstimateCrbCurveTable)
%LOCALBUILDCLEANPLOTDATA Store all summary plot inputs in scanData.

plotData = struct();
plotData.cleanEstimateCrbCurveTable = cleanEstimateCrbCurveTable;
end

function plotData = localBuildScalePlotData(scaleJointTrimAggregateTable)
%LOCALBUILDSCALEPLOTDATA Store lightweight scale surface plot inputs.

plotData = struct();
plotData.scaleJointTrimAggregateTable = scaleJointTrimAggregateTable;
end

function plotData = localPlotScaleScan(scanData)
%LOCALPLOTSCALESCAN Plot compact scale-surface heatmaps for the first SNR point.

plotData = localGetFieldOrDefault(scanData, 'plotData', struct());
if ~isfield(plotData, 'scaleJointTrimAggregateTable')
  plotData = localBuildScalePlotData(scanData.scaleJointTrimAggregateTable);
end
aggTable = plotData.scaleJointTrimAggregateTable;
if isempty(aggTable) || ~istable(aggTable) || height(aggTable) == 0
  return;
end
try
  snrList = unique(aggTable.snrDb, 'stable');
  snrUse = snrList(1);
  methodList = unique(aggTable.displayName, 'stable');
  for iMethod = 1:numel(methodList)
    methodMask = aggTable.displayName == methodList(iMethod) & aggTable.snrDb == snrUse;
    methodTable = aggTable(methodMask, :);
    if height(methodTable) == 0
      continue;
    end
    doaList = unique(methodTable.doaCrbScale, 'stable');
    fdList = unique(methodTable.fdRefCrbScale, 'stable');
    angleMat = NaN(numel(fdList), numel(doaList));
    fdMat = NaN(numel(fdList), numel(doaList));
    for iRow = 1:height(methodTable)
      iDoa = find(doaList == methodTable.doaCrbScale(iRow), 1, 'first');
      iFd = find(fdList == methodTable.fdRefCrbScale(iRow), 1, 'first');
      angleMat(iFd, iDoa) = methodTable.angleRmseOverCrb(iRow);
      fdMat(iFd, iDoa) = methodTable.fdRefRmseOverCrb(iRow);
    end
    figure('Name', char("Scale surface " + methodList(iMethod)));
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    imagesc(doaList, fdList, angleMat);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('DoA CRB scale');
    ylabel('fdRef CRB scale');
    title(sprintf('%s angle RMSE/CRB, SNR %.1f dB', char(methodList(iMethod)), snrUse), 'Interpreter', 'none');
    nexttile;
    imagesc(doaList, fdList, fdMat);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('DoA CRB scale');
    ylabel('fdRef CRB scale');
    title(sprintf('%s fdRef RMSE/CRB, SNR %.1f dB', char(methodList(iMethod)), snrUse), 'Interpreter', 'none');
  end
catch ME
  warning('scanMfMsMleCrbCleanBoundScale:PlotFailed', ...
    'Scale-surface plotting failed: %s', ME.message);
end
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointSummary, checkpointDir)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Convert stored checkpoint summary to a table.

checkpointSummaryTable = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
  zeros(0, 1), zeros(0, 1), false(0, 1), false(0, 1), ...
  'VariableNames', {'runName', 'runKey', 'runDir', 'numTask', 'numDone', ...
  'isComplete', 'cleanedOnSuccess'});
if nargin < 1 || ~isstruct(checkpointSummary) || isempty(checkpointSummary)
  return;
end
runNameList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runName', strings(0, 1))), [], 1);
runKeyList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runKey', strings(0, 1))), [], 1);
runDirList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runDir', strings(0, 1))), [], 1);
numTaskList = reshape(double(localGetFieldOrDefault(checkpointSummary, 'numTask', zeros(0, 1))), [], 1);
numDoneList = reshape(double(localGetFieldOrDefault(checkpointSummary, 'numDone', zeros(0, 1))), [], 1);
isCompleteList = reshape(logical(localGetFieldOrDefault(checkpointSummary, 'isComplete', false(0, 1))), [], 1);
if isempty(runNameList) && nargin >= 2 && strlength(string(checkpointDir)) > 0
  runNameList = "";
  runKeyList = "";
  runDirList = string(checkpointDir);
  numTaskList = 0;
  numDoneList = 0;
  isCompleteList = false;
end
numRow = max([numel(runNameList), numel(runKeyList), numel(runDirList), ...
  numel(numTaskList), numel(numDoneList), numel(isCompleteList)]);
if numRow == 0
  return;
end
runNameList = localPadColumn(runNameList, numRow, "");
runKeyList = localPadColumn(runKeyList, numRow, "");
runDirList = localPadColumn(runDirList, numRow, "");
numTaskList = localPadColumn(numTaskList, numRow, 0);
numDoneList = localPadColumn(numDoneList, numRow, 0);
isCompleteList = localPadColumn(isCompleteList, numRow, false);
cleanedOnSuccess = false(numRow, 1);
if isfield(checkpointSummary, 'cleanedOnSuccess') && ~isempty(checkpointSummary.cleanedOnSuccess)
  cleanedOnSuccess(:) = logical(checkpointSummary.cleanedOnSuccess);
end
checkpointSummaryTable = table(runNameList, runKeyList, runDirList, numTaskList, ...
  numDoneList, isCompleteList, cleanedOnSuccess, 'VariableNames', ...
  {'runName', 'runKey', 'runDir', 'numTask', 'numDone', 'isComplete', 'cleanedOnSuccess'});
end

function value = localPadColumn(value, numRow, fillValue)
%LOCALPADCOLUMN Pad or trim a column vector for lightweight summary tables.

value = reshape(value, [], 1);
if numel(value) >= numRow
  value = value(1:numRow);
  return;
end
value(end + 1:numRow, 1) = fillValue;
end

function localPlotCleanEstimateVsCrbCurves(curveTable)
%LOCALPLOTCLEANESTIMATEVSCRBCURVES Plot static/MF estimates and CRB on log axes.

if isempty(curveTable) || ~istable(curveTable) || height(curveTable) == 0
  return;
end
methodNameList = ["SS-SF-Static"; "MS-SF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "MS-MF-CP-K"; "MS-MF-CP-U"; "SS-MF-IP-K"; "SS-MF-IP-U"; "MS-MF-IP-K"; "MS-MF-IP-U"];
methodLabelList = ["SS Static"; "MS Static"; "SS-CP-K"; "SS-CP-U"; "MS-CP-K"; "MS-CP-U"; "SS-IP-K"; "SS-IP-U"; "MS-IP-K"; "MS-IP-U"];
markerList = {'o', 's', '^', 'v', 'd', 'p', '>', '<', 'x', '+'};
filterLabel = string(curveTable.filterName(1));
try
  figure('Name', 'Clean static/MF estimate vs CRB');
  tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

  nexttile;
  hold on;
  for iMethod = 1:numel(methodNameList)
    methodMask = curveTable.displayName == methodNameList(iMethod);
    methodTable = sortrows(curveTable(methodMask, :), 'snrDb');
    if height(methodTable) == 0
      continue;
    end
    marker = markerList{1 + mod(iMethod - 1, numel(markerList))};
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.angleRmseDeg), ...
      ['-' marker], 'DisplayName', methodLabelList(iMethod) + " est");
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.angleCrbSphericalMedianDeg), ...
      ['--' marker], 'DisplayName', methodLabelList(iMethod) + " CRB");
  end
  xlabel('SNR (dB)');
  ylabel('Spherical angle RMSE and CRB std (deg)');
  title(sprintf('Angle, %s', char(filterLabel)));
  set(gca, 'YScale', 'log');
  grid on;
  legend('Location', 'best');

  nexttile;
  hold on;
  for iMethod = 1:numel(methodNameList)
    methodMask = curveTable.displayName == methodNameList(iMethod);
    methodTable = sortrows(curveTable(methodMask, :), 'snrDb');
    if height(methodTable) == 0
      continue;
    end
    marker = markerList{1 + mod(iMethod - 1, numel(markerList))};
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.fdRefRmseHz), ...
      ['-' marker], 'DisplayName', methodLabelList(iMethod) + " est");
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.fdRefCrbMedianHz), ...
      ['--' marker], 'DisplayName', methodLabelList(iMethod) + " CRB");
  end
  xlabel('SNR (dB)');
  ylabel('fdRef RMSE and CRB std (Hz)');
  title(sprintf('Reference Doppler, %s', char(filterLabel)));
  set(gca, 'YScale', 'log');
  grid on;
  legend('Location', 'best');
catch ME
  warning('scanMfMsMleCrbCleanBoundScale:PlotFailed', ...
    'Clean estimate-vs-CRB plot failed: %s', ME.message);
end
end

function value = localPositiveForLog(value)
%LOCALPOSITIVEFORLOG Replace nonpositive values with NaN for log-scale plots.

value(value <= 0 | ~isfinite(value)) = NaN;
end

function mask = localFilterMask(caseTable, filterName)
%LOCALFILTERMASK Return row filter for clean comparison.

filterName = string(filterName);
switch filterName
  case "raw"
    mask = true(height(caseTable), 1);
  case "health"
    mask = caseTable.healthResolved;
  case "joint-trim"
    mask = caseTable.jointTrimKeep;
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:UnknownFilter', 'Unsupported filter "%s".', char(filterName));
end
end

function value = localMsErrorOverMatchedSsCrb(caseTable, msMask, ssRaw)
%LOCALMSERROROVERMATCHEDSSCRB Normalize MS angle error with same-seed SS-MF CRB.

msIdx = find(msMask(:));
ratioList = NaN(numel(msIdx), 1);
for i = 1:numel(msIdx)
  idx = msIdx(i);
  ssIdx = find(ssRaw & caseTable.taskSeed == caseTable.taskSeed(idx), 1, 'first');
  if isempty(ssIdx)
    continue;
  end
  ssCrb = caseTable.crbSphericalApproxStdDeg(ssIdx);
  if isfinite(ssCrb) && ssCrb > 0
    ratioList(i) = (caseTable.sphericalAngleErrDeg(idx) ./ ssCrb).^2;
  end
end
value = mean(ratioList(isfinite(ratioList)), 'omitnan');
end

function auditTable = localBuildCleanSearchRangeAuditTable(caseTable)
%LOCALBUILDCLEANSEARCHRANGEAUDITTABLE Check clean local boxes against CRB scale.

if isempty(caseTable) || height(caseTable) == 0
  auditTable = table();
  return;
end
methodList = unique(caseTable.displayName, 'stable');
rowTemplate = struct('displayName', "", 'snrDb', NaN, 'numRepeat', NaN, ...
  'minDoaLatHalfOverCrb', NaN, 'minDoaLonHalfOverCrb', NaN, ...
  'minFdRefHalfOverCrb', NaN, 'doaBelowCrb', false, 'fdRefBelowCrb', false, ...
  'doaTruthBoxViolationRate', NaN);
snrList = unique(caseTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(methodList) * numel(snrList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iMethod = 1:numel(methodList)
    idx = idx + 1;
    mask = caseTable.displayName == methodList(iMethod) & caseTable.snrDb == snrDb;
    row = rowTemplate;
    row.displayName = methodList(iMethod);
    row.snrDb = snrDb;
    row.numRepeat = nnz(mask);
    row.minDoaLatHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLatDeg(mask), caseTable.crbLatStdDeg(mask)), [], 'omitnan');
    row.minDoaLonHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLonDeg(mask), caseTable.crbLonStdDeg(mask)), [], 'omitnan');
    row.minFdRefHalfOverCrb = min(localSafeVectorRatio(caseTable.fdRangeHalfWidthHz(mask), caseTable.fdRefCrbStdHz(mask)), [], 'omitnan');
    row.doaBelowCrb = row.minDoaLatHalfOverCrb < 1 || row.minDoaLonHalfOverCrb < 1;
    row.fdRefBelowCrb = row.minFdRefHalfOverCrb < 1;
    if ismember('doaTruthBoxViolation', caseTable.Properties.VariableNames)
      row.doaTruthBoxViolationRate = mean(double(caseTable.doaTruthBoxViolation(mask)), 'omitnan');
    end
    rowList(idx) = row;
  end
end
auditTable = struct2table(rowList(:));
end

function outlierTable = localBuildCleanOutlierTable(caseTable)
%LOCALBUILDCLEANOUTLIERTABLE Keep MS-MF unknown-rate rows rejected by clean trim.

if isempty(caseTable) || height(caseTable) == 0
  outlierTable = table();
  return;
end
mask = startsWith(caseTable.displayName, "MS-MF-") & endsWith(caseTable.displayName, "-U") & ~caseTable.jointTrimKeep;
columnList = {'displayName', 'snrDb', 'taskSeed', 'sphericalAngleErrDeg', ...
  'angleErrOverSphericalCrb', 'fdRefAbsErrHz', 'fdRefErrOverCrb', ...
  'fdRateAbsErrHzPerSec', 'iterations', 'firstOrderOpt', 'objectiveImprove', ...
  'fdRefBoundaryHit', 'fdRateBoundaryHit', 'healthResolved', 'jointTrimKeep', 'rejectReason'};
columnList = columnList(ismember(columnList, caseTable.Properties.VariableNames));
outlierTable = caseTable(mask, columnList);
if height(outlierTable) > 0
  outlierTable = sortrows(outlierTable, {'angleErrOverSphericalCrb'}, {'descend'});
end
end

function localPrintReplayConfig(config)
%LOCALPRINTREPLAYCONFIG Print clean scale-scan axes.

fprintf('  %-32s : %s\n', 'SNR list (dB)', localFormatRow(config.snrDbList));
fprintf('  %-32s : %d\n', 'frame count', config.numFrame);
fprintf('  %-32s : %s\n', 'user LLA', localFormatRow(config.contextUsrLla));
fprintf('  %-32s : %s\n', 'selected satellites', localFormatRow(config.contextSelectedSatIdxGlobal));
fprintf('  %-32s : %.0f\n', 'reference satellite', config.contextRefSatIdxGlobal);
fprintf('  %-32s : %s\n', 'method list', strjoin(cellstr(config.methodNameList), ', '));
fprintf('  %-32s : %s\n', 'DoA CRB scale list', localFormatRow(config.truthLocalDoaCrbScaleList));
fprintf('  %-32s : %s\n', 'fdRef CRB scale list', localFormatRow(config.truthLocalFdRefCrbScaleList));
fprintf('  %-32s : %s / %.2f tooth cap\n', 'fdRef range mode', ...
  char(config.fdRefRangeMode), config.fdRefMaxHalfToothFraction);
fprintf('  %-32s : %.2f Hz/s\n', 'fdRate truth-local half-width', ...
  config.truthLocalFdRateHalfWidthHzPerSec);
fprintf('  %-32s : %s\n', 'estimation chain', char(config.initMode));
fprintf('  %-32s : %s / %s\n', 'MF objective / route', ...
  char(config.dynamicObjectiveMode), char(config.dynamicRouteMode));
fprintf('  %-32s : %s / %s\n', 'MF signal / steering', ...
  char(config.cleanMfSignalTimeModel), char(config.cleanMfSteeringMode));
fprintf('  %-32s : %.2f / %.2f CRB\n', 'trim angle/fdRef cap', config.trimAngleNormMax, config.trimFdRefNormMax);
fprintf('  %-32s : %d / %d\n', 'task count / outer parfor', config.numTask, logical(config.useParfor));
end

function localPrintCompactTable(sectionTitle, dataTable, edgeCount)
%LOCALPRINTCOMPACTTABLE Print a compact scan table through common report helper.

if nargin < 3
  edgeCount = 12;
end
printMfReportTableSection(char(sectionTitle), dataTable, ...
  'MaxPreviewRow', edgeCount, 'SectionPrinter', @printMfScanSection);
end

function lineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build short notification metrics.

lineList = strings(0, 1);
summaryTable = localGetFieldOrDefault(scanData, 'scaleSurfaceSummaryTable', table());
if istable(summaryTable) && height(summaryTable) > 0
  passRate = mean(double(summaryTable.targetPass), 'omitnan');
  lineList(end + 1, 1) = sprintf('scale surface: methods=%d, targetPassRate=%.3g', ...
    height(summaryTable), passRate); %#ok<AGROW>
end
readableTable = localGetFieldOrDefault(scanData, 'scaleReadableRmseCrbTable', table());
if istable(readableTable) && height(readableTable) > 0
  lineList(end + 1, 1) = sprintf('grid rows=%d, median trim keep=%.3g', ...
    height(readableTable), median(readableTable.trimKeepRate, 'omitnan')); %#ok<AGROW>
end
lineList(end + 1, 1) = sprintf('elapsed %.1f min', scanData.elapsedSec / 60); %#ok<AGROW>
end

function checkpointOpt = localBuildCheckpointOpt(scanName, config, useParfor)
%LOCALBUILDCHECKPOINTOPT Build common checkpoint options with a stable run key.

meta = struct('snrDbList', reshape(config.snrDbList, 1, []), ...
  'seedList', reshape(config.seedList, 1, []), ...
  'truthLocalDoaCrbScaleList', reshape(config.truthLocalDoaCrbScaleList, 1, []), ...
  'truthLocalFdRefCrbScaleList', reshape(config.truthLocalFdRefCrbScaleList, 1, []), ...
  'methodNameList', reshape(config.methodNameList, 1, []), ...
  'dynamicObjectiveMode', char(config.dynamicObjectiveMode), ...
  'dynamicRouteMode', char(config.dynamicRouteMode), ...
  'routeTag', char(config.routeTag));
checkpointOpt = buildMfReplayCheckpointOpt(scanName, config, struct( ...
  'runKey', config.runKey, 'meta', meta, 'useParfor', logical(useParfor)));
end

function runKey = localBuildRunKey(config)
%LOCALBUILDRUNKEY Build path-safe run key.

runKey = sprintf('frames%d_snr%s_seed%dto%d_rep%d_%s', config.numFrame, ...
  localCompactNumberList(config.snrDbList), config.seedList(1), config.seedList(end), ...
  numel(config.seedList), localCompactStringHash(localBuildRunSignature(config)));
runKey = string(runKey);
end

function signature = localBuildRunSignature(config)
%LOCALBUILDRUNSIGNATURE Build deterministic signature for checkpoint compatibility.

signature = strjoin(string({ ...
  sprintf('frames%d', config.numFrame), ...
  sprintf('snr%s', mat2str(reshape(double(config.snrDbList), 1, []), 8)), ...
  sprintf('seed%s', mat2str(reshape(double(config.seedList), 1, []), 8)), ...
  sprintf('methods%s', char(strjoin(reshape(string(config.methodNameList), 1, []), ','))), ...
  sprintf('doaScale%s', mat2str(reshape(double(config.truthLocalDoaCrbScaleList), 1, []), 8)), ...
  sprintf('fdRefScale%s', mat2str(reshape(double(config.truthLocalFdRefCrbScaleList), 1, []), 8)), ...
  sprintf('fdRefRangeMode%s', char(string(config.fdRefRangeMode))), ...
  sprintf('init%s', char(string(config.initMode))), ...
  sprintf('objective%s', char(string(config.dynamicObjectiveMode))), ...
  sprintf('route%s', char(string(config.dynamicRouteMode))), ...
  sprintf('trimAngle%.8g', config.trimAngleNormMax), ...
  sprintf('trimFdRef%.8g', config.trimFdRefNormMax), ...
  char(config.routeTag) ...
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

function textValue = localCompactNumberList(valueList)
%LOCALCOMPACTNUMBERLIST Format numeric axes for a checkpoint run key.

valueList = reshape(double(valueList), 1, []);
textValue = strjoin(compose('%.6g', valueList), '_');
textValue = replace(textValue, '-', 'm');
textValue = replace(textValue, '.', 'p');
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Check whether the outer scan task loop can use parfor.

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

function flowOpt = localBuildFlowOpt(parallelOpt)
%LOCALBUILDFLOWOPT Build static and dynamic options shared by clean methods.

flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false, ...
  'staticMsHalfWidth', [1e-8; 1e-8]));
flowOpt = applyDynamicTransitionFlowDefaults(flowOpt);
end

function parallelOpt = localBuildSerialInnerParallelOpt()
%LOCALBUILDSERIALINNERPARALLELOPT Disable estimator-inner parallelism.

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if numFrame < 2
  error('scanMfMsMleCrbCleanBoundScale:InvalidFrameCount', ...
    'Frame count must be at least two for multi-frame diagnostics.');
end
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function crbBundle = localBuildCrbBundleQuiet(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%LOCALBUILDCRBBUNDLEQUIET Build CP/IP CRB while suppressing expected warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
crbBundle = localAddIndependentPhaseCrb(crbBundle, periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
end

function crbBundle = localAddIndependentPhaseCrb(crbBundle, periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%LOCALADDINDEPENDENTPHASECRB Add local IP-K/IP-U CRB entries for this scan.

truth = periodicFixture.truth;
pathGainRef = ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
noiseVarRef = noiseVar * ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
pathGainMs = ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);
noiseVarMs = noiseVar * ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);

ipKnownOpt = struct();
ipKnownOpt.doaType = 'latlon';
ipKnownOpt.phaseMode = 'independent';
ipKnownOpt.fdRateMode = 'known';
ipKnownOpt.steeringMode = 'framewise';
ipKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
[crbBundle.crbMfRefIpKnown, crbBundle.auxCrbMfRefIpKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, ipKnownOpt);
[crbBundle.crbMfMsIpKnown, crbBundle.auxCrbMfMsIpKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, ipKnownOpt);

ipUnknownOpt = ipKnownOpt;
ipUnknownOpt.fdRateMode = 'unknown';
[crbBundle.crbMfRefIpUnknown, crbBundle.auxCrbMfRefIpUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, ipUnknownOpt);
[crbBundle.crbMfMsIpUnknown, crbBundle.auxCrbMfMsIpUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, ipUnknownOpt);
end

function [crb, aux] = localTryBuildDynamicCrb(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt)
%LOCALTRYBUILDDYNAMICCRB Build one dynamic CRB and return NaNs on failure.

try
  [crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = getDoaDopplerFieldOrDefault(crbOpt, 'fdRateMode', "unknown");
  aux.phaseMode = getDoaDopplerFieldOrDefault(crbOpt, 'phaseMode', "independent");
  aux.errorMessage = string(ME.message);
end
end

function [crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName)
%LOCALRESOLVECRBFORMETHOD Return the CRB block matching one method.

methodName = string(methodName);
methodName = localBaseDynamicMethodName(methodName);
switch methodName
  case "SS-SF-Static"
    crbFull = crbBundle.crbSfRef;
  case "MS-SF-Static"
    crbFull = crbBundle.crbSfMs;
  case "SS-MF-CP-K"
    crbFull = crbBundle.crbMfRefKnown;
  case "SS-MF-CP-U"
    crbFull = crbBundle.crbMfRefUnknown;
  case "MS-MF-CP-K"
    crbFull = crbBundle.crbMfMsKnown;
  case "MS-MF-CP-U"
    crbFull = crbBundle.crbMfMsUnknown;
  case "SS-MF-IP-K"
    crbFull = crbBundle.crbMfRefIpKnown;
  case "SS-MF-IP-U"
    crbFull = crbBundle.crbMfRefIpUnknown;
  case "MS-MF-IP-K"
    crbFull = crbBundle.crbMfMsIpKnown;
  case "MS-MF-IP-U"
    crbFull = crbBundle.crbMfMsIpUnknown;
  otherwise
    error('scanMfMsMleCrbCleanBoundScale:UnknownCrbMethod', 'Unsupported CRB method "%s".', char(methodName));
end
crbDoa = real(crbFull(1:2, 1:2));
end

function method = localApplyDoaCrbScaledWidth(method, crbBundle)
%LOCALAPPLYDOACRBSCALEDWIDTH Set the local DoA box from the method CRB.

method.doaHalfWidthDeg = localResolveDoaCrbScaledHalfWidth( ...
  method.baseDisplayName, method.doaHalfWidthDeg, crbBundle);
end

function doaHalfWidthDeg = localResolveDoaCrbScaledHalfWidth(methodName, doaCrbScale, crbBundle)
%LOCALRESOLVEDOACRBSCALEDHALFWIDTH Build a method-specific CRB-scaled DoA box.

doaCrbScale = double(doaCrbScale);
if isempty(doaCrbScale) || ~isfinite(doaCrbScale(1)) || doaCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundScale:InvalidDoaCrbScale', ...
    'The clean DoA CRB scale must be a positive finite scalar.');
end
try
  [crbDoa, ~] = localResolveCrbForMethod(crbBundle, methodName);
catch
  doaHalfWidthDeg = NaN(2, 1);
  return;
end
doaStdDeg = sqrt(max(real(diag(crbDoa)), 0));
if numel(doaStdDeg) < 2 || any(~isfinite(doaStdDeg(1:2)))
  doaHalfWidthDeg = NaN(2, 1);
  return;
end
doaHalfWidthDeg = doaCrbScale(1) * reshape(doaStdDeg(1:2), [], 1);
end

function fdRangeUse = localBuildFdRefCrbScaledRange(truthFdRefHz, crbBundle, ...
  methodName, fdRefCrbScale, toothStepHz, maxHalfToothFraction)
%LOCALBUILDFDREFCRBSCALEDRANGE Build a pure CRB-scaled fdRef box.

fdRefCrbScale = double(fdRefCrbScale);
if isempty(fdRefCrbScale) || ~isfinite(fdRefCrbScale(1)) || fdRefCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundScale:InvalidFdRefCrbScale', ...
    'The fdRef CRB scale must be a positive finite scalar.');
end
[~, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
if ~(isfinite(fdRefCrbStdHz) && fdRefCrbStdHz > 0)
  error('scanMfMsMleCrbCleanBoundScale:InvalidFdRefCrbStd', ...
    'Cannot build a CRB-scaled fdRef box for method "%s".', char(methodName));
end
fdHalfWidthHz = fdRefCrbScale(1) * fdRefCrbStdHz;
maxHalfToothFraction = double(maxHalfToothFraction);
if isfinite(maxHalfToothFraction) && maxHalfToothFraction > 0 && isfinite(toothStepHz) && toothStepHz > 0
  fdHalfWidthHz = min(fdHalfWidthHz, maxHalfToothFraction * toothStepHz);
end
fdRangeUse = truthFdRefHz + fdHalfWidthHz * [-1, 1];
end

function fdRangeUse = localApplyFdRefCrbScaledFloor(fdRangeUse, truthFdRefHz, crbBundle, methodNameList, fdRefCrbScale)
%LOCALAPPLYFDREFCRBSCALEDFLOOR Ensure the fdRef box is not a one-sigma CRB truncation.

fdRefCrbScale = double(fdRefCrbScale);
if isempty(fdRefCrbScale) || ~isfinite(fdRefCrbScale(1)) || fdRefCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundScale:InvalidFdRefCrbScale', ...
    'The clean fdRef CRB floor scale must be a positive finite scalar.');
end
fdHalfWidthHz = 0.5 * abs(fdRangeUse(2) - fdRangeUse(1));
fdRefCrbFloorHz = NaN;
for iMethod = 1:numel(methodNameList)
  methodName = string(methodNameList(iMethod));
  try
    [~, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
    fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
    if isfinite(fdRefCrbStdHz)
      fdRefCrbFloorHz = max([fdRefCrbFloorHz; fdRefCrbScale(1) * fdRefCrbStdHz], [], 'omitnan');
    end
  catch
  end
end
if isfinite(fdRefCrbFloorHz) && fdHalfWidthHz < fdRefCrbFloorHz
  fdRangeUse = truthFdRefHz + fdRefCrbFloorHz * [-1, 1];
end
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve Doppler tooth spacing from frame offsets.

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

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for clean case metrics.

row = struct('displayName', "", 'satMode', "", 'frameMode', "", 'phaseMode', "", ...
  'fdRateMode', "", 'mfInitMode', "", 'dynamicObjectiveMode', "", ...
  'dynamicRouteMode', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'doaCrbScale', NaN, 'fdRefHalfToothFraction', NaN, ...
  'fdRefCrbScale', NaN, 'fdRateHalfWidthHzPerSec', NaN, 'numFrame', NaN, ...
  'frameIntvlSec', NaN, 'windowSec', NaN, 'toothStepHz', NaN, ...
  'truthLatDeg', NaN, 'truthLonDeg', NaN, 'initLatDeg', NaN, 'initLonDeg', NaN, ...
  'finalLatDeg', NaN, 'finalLonDeg', NaN, 'latErrDeg', NaN, 'lonErrDeg', NaN, ...
  'initAngleErrDeg', NaN, 'finalAngleErrDeg', NaN, 'sphericalAngleErrDeg', NaN, ...
  'finalMinusInitAngleDeg', NaN, 'searchHalfWidthLatDeg', NaN, ...
  'searchHalfWidthLonDeg', NaN, 'seedLatDeg', NaN, 'seedLonDeg', NaN, ...
  'boundCenterLatDeg', NaN, 'boundCenterLonDeg', NaN, ...
  'doaLatBoxViolation', false, 'doaLonBoxViolation', false, ...
  'doaTruthBoxViolation', false, 'crbLatStdDeg', NaN, 'crbLonStdDeg', NaN, ...
  'crbTraceStdDeg', NaN, 'crbSphericalApproxStdDeg', NaN, ...
  'crbAngleMetricVarRad2', NaN, 'crbMetricTraceRatio', NaN, ...
  'angleErrOverSphericalCrb', NaN, 'angleErrOverTraceCrb', NaN, ...
  'fdRefTrueHz', NaN, 'initFdRefHz', NaN, 'finalFdRefHz', NaN, ...
  'fdRefErrHz', NaN, 'fdRefAbsErrHz', NaN, 'fdRefCrbStdHz', NaN, ...
  'fdRefErrOverCrb', NaN, 'fdRefInitMoveHz', NaN, 'fdRangeLowHz', NaN, ...
  'fdRangeHighHz', NaN, 'fdRangeHalfWidthHz', NaN, 'fdRefBoundaryHit', false, ...
  'fdRateTrueHzPerSec', NaN, 'initFdRateHzPerSec', NaN, ...
  'finalFdRateHzPerSec', NaN, 'fdRateErrHzPerSec', NaN, ...
  'fdRateAbsErrHzPerSec', NaN, 'fdRateInitMoveHzPerSec', NaN, ...
  'fdRateRangeLowHzPerSec', NaN, 'fdRateRangeHighHzPerSec', NaN, ...
  'fdRateBoundaryHit', false, 'fdLineRmseHz', NaN, 'solverResolved', false, ...
  'exitflag', NaN, 'iterations', NaN, 'funcCount', NaN, 'stepSize', NaN, ...
  'firstOrderOpt', NaN, 'constrViolation', NaN, 'algorithm', "", ...
  'usedScaledSolve', false, 'fallbackTriggered', false, 'solveVariant', "", ...
  'candidateCount', NaN, 'candidateVariantText', "", 'initObj', NaN, ...
  'finalObj', NaN, 'objectiveImprove', NaN, 'finalResidualNorm', NaN, ...
  'refCoherence', NaN, 'nonRefCoherenceFloor', NaN, 'runTimeMs', NaN, ...
  'isDynamicMethod', false, 'solverHealthy', false, 'freqHealthy', false, ...
  'rateHealthy', false, 'boundaryHealthy', false, 'notNoSolve', false, ...
  'healthResolved', false, 'trimKeep', false, 'jointTrimKeep', false, ...
  'fdRefNormTail', false, 'rejectReason', "", 'warningFlag', false);
end

function progressTracker = localCreateProgressTracker(totalCount)
%LOCALCREATEPROGRESSTRACKER Create a best-effort scan progress tracker.

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

function localAdvanceProgress(progressTracker)
%LOCALADVANCEPROGRESS Advance progress when progressbar is available.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('advance');
catch
end
end

function localCloseProgressTracker(progressTracker)
%LOCALCLOSEPROGRESSTRACKER Close a best-effort scan progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a struct field or return a fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide only when denominator is finite and positive.

if isfinite(numerator) && isfinite(denominator) && denominator > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function ratio = localSafeVectorRatio(numerator, denominator)
%LOCALSAFEVECTORRATIO Element-wise safe ratio for table columns.

numerator = reshape(double(numerator), [], 1);
denominator = reshape(double(denominator), [], 1);
ratio = NaN(size(numerator));
mask = isfinite(numerator) & isfinite(denominator) & denominator > 0;
ratio(mask) = numerator(mask) ./ denominator(mask);
end

function value = localMse(x)
%LOCALMSE Mean-square of finite entries.

x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = mean(x(:).^2);
end
end

function value = localRmse(x)
%LOCALRMSE Root-mean-square of finite entries.

value = sqrt(localMse(x));
end

function value = localNormalizedMse(err, crbStd)
%LOCALNORMALIZEDMSE Mean squared CRB-normalized error over finite entries.

err = reshape(double(err), [], 1);
crbStd = reshape(double(crbStd), [], 1);
mask = isfinite(err) & isfinite(crbStd) & crbStd > 0;
if ~any(mask)
  value = NaN;
else
  value = mean((err(mask) ./ crbStd(mask)).^2);
end
end

function key = localMethodKey(methodName)
%LOCALMETHODKEY Convert a display name into a struct field key.

key = matlab.lang.makeValidName(char(regexprep(string(methodName), '[^A-Za-z0-9]', '_')));
end

function out = ternary(conditionValue, trueValue, falseValue)
%TERNARY Return one of two values based on a scalar condition.

if conditionValue
  out = trueValue;
else
  out = falseValue;
end
end

function textValue = localFormatRow(valueList)
%LOCALFORMATROW Format numeric/string row values.

if isstring(valueList) || ischar(valueList)
  textValue = strjoin(cellstr(reshape(string(valueList), 1, [])), ', ');
  return;
end
valueList = reshape(double(valueList), 1, []);
if isempty(valueList)
  textValue = '';
else
  textValue = strjoin(compose('%.6g', valueList), ', ');
end
end

function contextSummary = localBuildContextSummary(context)
%LOCALBUILDCONTEXTSUMMARY Keep a lightweight context summary.

contextSummary = struct();
contextSummary.carrierFreq = localGetFieldOrDefault(context, 'carrierFreq', NaN);
contextSummary.wavelen = localGetFieldOrDefault(context, 'wavelen', NaN);
contextSummary.otherSatIdxGlobal = localGetFieldOrDefault(context, 'otherSatIdxGlobal', NaN);
if isfield(context, 'waveInfo')
  contextSummary.sampleRate = localGetFieldOrDefault(context.waveInfo, 'sampleRate', NaN);
end
end

function repeatCell = localStripRepeatCell(repeatCell)
%LOCALSTRIPREPEATCELL Drop heavyweight task data before snapshot saving.

for iRepeat = 1:numel(repeatCell)
  if isempty(repeatCell{iRepeat})
    continue;
  end
  if isfield(repeatCell{iRepeat}, 'truth')
    repeatCell{iRepeat} = rmfield(repeatCell{iRepeat}, 'truth');
  end
end
end
