% replayMfMsMleCrbCleanTrim
% Purpose: run a clean MS-vs-SS truth-centered MF MLE/CRB replay and report
% raw, health-trimmed, and CRB-local trimmed metrics from the formal
% estimator path only.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted task runs from per-task files under the repo-root tmp.
% saveSnapshot=true saves only replayData via saveExpSnapshot. Telegram notice
% is best-effort only and never affects numerical results.

clear; close all; clc;

%% Replay configuration

replayName = "replayMfMsMleCrbCleanTrim";
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = true;
optVerbose = false;

% Dynamic MF estimator dialect. The default is the clean MLE core: no CP
% engineering penalties, no fd alias unwrap, and no estimator-internal
% basin-entry / warm-anchor route. Keep this replay on core-only unless
% temporarily running a manual route diagnostic.
dynamicObjectiveMode = "pure-mle";
dynamicRouteMode = "core-only";

numFrame = 10;
% Use this replay as the lightweight clean SNR sweep. Keep numRepeat at 100
% for stable K/U baseline comparison; lower it temporarily for smoke tests.
snrDbList = (-15:5:10).';
baseSeed = 253;
numRepeat = 200;
seedList = baseSeed + (0:(numRepeat - 1));
methodNameList = [
  "SS-SF-Static"
  "MS-SF-Static"
  "SS-MF-CP-K"
  "SS-MF-CP-U"
  "MS-MF-CP-K"
  "MS-MF-CP-U"
];

% Truth-local frequency range for estimator modification diagnostics.
% MF keeps the requested in-tooth box. Static uses the same request plus a
% CRB-scaled floor, because a 1-sigma truth-centered fdRef box truncates the
% static fdRef error and can make MSE/CRB look artificially below one.
truthLocalFdHalfToothFraction = 0.3;
truthLocalFdRefCrbScale = 3.0;
truthLocalFdRateHalfWidthHzPerSec = 1000;

% DoA hard envelope used by the clean static-to-MF chain. Each method and
% SNR uses a truth-centered box whose half-width is this common multiplier
% times the corresponding CRB standard deviation.
truthLocalDoaCrbScale = 2;
boundaryTolFraction = 0.02;
healthFdRefToothFraction = 0.25;
healthFdRateAbsMaxHzPerSec = 250;
trimAngleNormMax = 5;
trimFdRefNormMax = 5;
crbTargetMseRatio = 1.1;

seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
methodNameList = reshape(string(methodNameList), [], 1);
truthLocalDoaCrbScale = double(truthLocalDoaCrbScale);
truthLocalFdRefCrbScale = double(truthLocalFdRefCrbScale);

runTic = tic;
replayData = struct();
config = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and replay options

  config.replayName = string(replayName);
  config.numFrame = numFrame;
  config.snrDbList = snrDbList;
  config.seedList = seedList;
  config.baseSeed = seedList(1);
  config.numRepeat = numel(seedList);
  config.methodNameList = methodNameList;
  config.truthLocalFdHalfToothFraction = truthLocalFdHalfToothFraction;
  config.truthLocalFdRefCrbScale = truthLocalFdRefCrbScale;
  config.truthLocalFdRateHalfWidthHzPerSec = truthLocalFdRateHalfWidthHzPerSec;
  config.truthLocalDoaCrbScale = truthLocalDoaCrbScale;
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
  config.initMode = "truth-crb-scaled-static-seed-fixed-bound";
  config.sweepMode = "snr-list";
  config.routeTag = "clean-truth-crb-scaled-static-mf-ku-snr-v5-fdref-floor-crb-audit" + ...
    string(config.truthLocalFdRefCrbScale) + "-" + ...
    config.dynamicObjectiveMode + "-" + config.dynamicRouteMode;

  taskList = localBuildTaskList(config.snrDbList, config.seedList);
  config.numTask = numel(taskList);
  config.useParfor = localCanUseParfor(config.numTask);
  config.runKey = localBuildRunKey(config);
  checkpointOpt = localBuildCheckpointOpt(replayName, config, config.useParfor);
  if config.checkpointEnable
    checkpointRunDir = string(checkpointOpt.runDir);
    config.checkpointRunDir = checkpointRunDir;
  end

  printMfReplayHeader(char(replayName), config, char(checkpointRunDir));
  localPrintReplayConfig(config);

  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(config.numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = localBuildSerialInnerParallelOpt();
  context = buildDynamicDualSatEciContext(contextOpt);
  context.subsetOffsetCell = {};
  context.subsetLabelList = strings(0, 1);
  context.numSubsetRandomTrial = 0;
  flowOpt = localBuildFlowOpt(context.parallelOpt);

  %% Run clean replay batch

  try
    [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, replayName);
  catch ME
    if strlength(string(checkpointRunDir)) > 0
      fprintf('Replay failed. Checkpoint artifacts kept at: %s\n', char(checkpointRunDir));
    end
    rethrow(ME);
  end

  caseTable = localBuildCaseTable(repeatCell);
  caseTable = localAnnotateHealthFlags(caseTable, config);
  cleanAggregateTable = buildMfMetricAggregateTable(caseTable);
  cleanHealthAggregateTable = buildMfFilteredMetricAggregateTable(caseTable, "health", caseTable.healthResolved);
  cleanTrimAggregateTable = buildMfFilteredMetricAggregateTable(caseTable, "joint-trim", caseTable.jointTrimKeep);
  cleanSsMsCompareTable = localBuildCleanSsMsCompareTable(caseTable, config);
  cleanKnownUnknownCompareTable = localBuildCleanKnownUnknownCompareTable(caseTable);
  cleanSnrCurveTable = localBuildCleanSnrCurveTable(cleanSsMsCompareTable);
  cleanEstimateCrbCurveTable = localBuildCleanEstimateCrbCurveTable( ...
    cleanTrimAggregateTable, cleanHealthAggregateTable, cleanAggregateTable);
  cleanSearchRangeAuditTable = localBuildCleanSearchRangeAuditTable(caseTable);
  cleanCrbPathAuditTable = localBuildCleanCrbPathAuditTable(repeatCell);
  cleanCrbPathAuditAggregateTable = localBuildCleanCrbPathAuditAggregateTable(cleanCrbPathAuditTable);
  cleanOutlierTable = localBuildCleanOutlierTable(caseTable);
  runtimeTable = buildMfRuntimeTable(repeatCell);
  runtimeAggregateTable = buildMfRuntimeAggregateTable(runtimeTable);
  topSlowRuntimeTable = buildMfTopSlowRuntimeTable(runtimeTable, 12);

  %% Data storage

  replayData = struct();
  replayData.replayName = string(replayName);
  replayData.runKey = config.runKey;
  replayData.utcRun = datetime('now', 'TimeZone', 'local');
  replayData.config = config;
  replayData.contextSummary = localBuildContextSummary(context);
  replayData.caseTable = caseTable;
  replayData.cleanAggregateTable = cleanAggregateTable;
  replayData.cleanHealthAggregateTable = cleanHealthAggregateTable;
  replayData.cleanTrimAggregateTable = cleanTrimAggregateTable;
  replayData.cleanSsMsCompareTable = cleanSsMsCompareTable;
  replayData.cleanKnownUnknownCompareTable = cleanKnownUnknownCompareTable;
  replayData.cleanSnrCurveTable = cleanSnrCurveTable;
  replayData.cleanEstimateCrbCurveTable = cleanEstimateCrbCurveTable;
  replayData.cleanSearchRangeAuditTable = cleanSearchRangeAuditTable;
  replayData.cleanCrbPathAuditTable = cleanCrbPathAuditTable;
  replayData.cleanCrbPathAuditAggregateTable = cleanCrbPathAuditAggregateTable;
  replayData.cleanOutlierTable = cleanOutlierTable;
  replayData.runtimeTable = runtimeTable;
  replayData.runtimeAggregateTable = runtimeAggregateTable;
  replayData.topSlowRuntimeTable = topSlowRuntimeTable;
  replayData.repeatCell = localStripRepeatCell(repeatCell);
  if config.checkpointEnable
    replayData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  end
  replayData.elapsedSec = toc(runTic);
  replayData = finalizeMfReplayResult(replayData, "");

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'replayData'}}, ...
      'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
    replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
  else
    replayData.snapshotFile = "";
  end

  if config.checkpointEnable
    if config.checkpointCleanupOnSuccess
      replayData.checkpointSummary.cleanupReport = cleanupPerfTaskGridCheckpoint(runState, struct('verbose', false));
      replayData.checkpointSummary.cleanedOnSuccess = true;
    else
      replayData.checkpointSummary.cleanedOnSuccess = false;
    end
  end

  %% Summary output

  if ~exist('replayData', 'var') || ~isstruct(replayData)
    error('replayMfMsMleCrbCleanTrim:MissingReplayData', ...
      'Load a snapshot containing replayData before running the summary section.');
  end
  replayData = localValidateReplayDataForSummary(replayData);
  replayConfigForReport = localGetFieldOrDefault(replayData, 'config', struct());
  replayNameForReport = string(localGetFieldOrDefault(replayData, 'replayName', "replayMfMsMleCrbCleanTrim"));
  replaySnapshotFile = string(localGetFieldOrDefault(replayData, 'snapshotFile', ""));
  replayElapsedSec = localGetFieldOrDefault(replayData, 'elapsedSec', NaN);
  checkpointDirForReport = string(localGetFieldOrDefault(replayConfigForReport, 'checkpointRunDir', ""));

  printMfReplaySection('Clean raw aggregate', replayData.cleanAggregateTable);
  printMfReplaySection('Clean health aggregate', replayData.cleanHealthAggregateTable);
  printMfReplaySection('Clean joint-trim aggregate', replayData.cleanTrimAggregateTable);
  printMfReplaySection('Clean MS-vs-SS comparison', replayData.cleanSsMsCompareTable);
  printMfReplaySection('Clean known-vs-unknown comparison', replayData.cleanKnownUnknownCompareTable);
  printMfReplaySection('Clean SNR curve summary', replayData.cleanSnrCurveTable);
  printMfReplaySection('Clean estimate-vs-CRB curve summary', replayData.cleanEstimateCrbCurveTable);
  printMfReplaySection('Clean search-range audit', replayData.cleanSearchRangeAuditTable);
  printMfReplaySection('Clean CRB path audit aggregate', replayData.cleanCrbPathAuditAggregateTable);
  localPrintCompactTable('Clean CRB path audit detail', replayData.cleanCrbPathAuditTable, 12);
  localPrintCompactTable('Clean outlier table', replayData.cleanOutlierTable, 12);
  printMfReplaySection('Clean runtime aggregate', replayData.runtimeAggregateTable);
  localPrintCompactTable('Top slow clean tasks', replayData.topSlowRuntimeTable, 12);
  fprintf('Replay purpose                    : clean truth-centered CRB-scaled estimator + SNR sweep / trim statistics.\n');
  trimAngleNormMaxForReport = localGetFieldOrDefault(replayConfigForReport, 'trimAngleNormMax', NaN);
  trimFdRefNormMaxForReport = localGetFieldOrDefault(replayConfigForReport, 'trimFdRefNormMax', NaN);
  fprintf('Trim definition                   : health + angle <= %.2f CRB + fdRef <= %.2f CRB.\n', ...
    trimAngleNormMaxForReport, trimFdRefNormMaxForReport);
  fprintf('Main decision table               : cleanSsMsCompareTable.\n');
  fprintf('Known-rate baseline table         : cleanKnownUnknownCompareTable.\n');
  localPlotCleanEstimateVsCrbCurves(replayData.cleanEstimateCrbCurveTable);

  notifyMfReplayStatus(struct( ...
    'replayName', replayNameForReport, ...
    'statusText', "DONE", ...
    'config', replayConfigForReport, ...
    'snapshotFile', replaySnapshotFile, ...
    'checkpointDir', checkpointDirForReport, ...
    'elapsedSec', replayElapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(replayData), ...
    'commentLineList', [ ...
      "Clean MS/SS trim replay completed; inspect cleanSsMsCompareTable and cleanKnownUnknownCompareTable before changing estimator logic."]));

catch ME
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'snapshotFile', "", ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'metricLineList', string(ME.identifier), ...
    'commentLineList', string(ME.message)));
  rethrow(ME);
end

%% Local helpers

function replayData = localValidateReplayDataForSummary(replayData)
%LOCALVALIDATEREPLAYDATAFORSUMMARY Validate stored replay data for rerunnable summary.

requiredFieldList = [
  "config"
  "caseTable"
  "cleanAggregateTable"
  "cleanHealthAggregateTable"
  "cleanTrimAggregateTable"
  "cleanSsMsCompareTable"
  "cleanSearchRangeAuditTable"
  "cleanOutlierTable"
  "runtimeAggregateTable"
  "topSlowRuntimeTable"
];
for iField = 1:numel(requiredFieldList)
  fieldName = char(requiredFieldList(iField));
  if ~isfield(replayData, fieldName)
    error('replayMfMsMleCrbCleanTrim:MissingReplayDataField', ...
      'replayData.%s is missing. Re-run the replay before rerunning this summary section.', fieldName);
  end
end
if ~isfield(replayData, 'cleanKnownUnknownCompareTable')
  replayData.cleanKnownUnknownCompareTable = localBuildCleanKnownUnknownCompareTable(replayData.caseTable);
end
if ~isfield(replayData, 'cleanSnrCurveTable')
  replayData.cleanSnrCurveTable = localBuildCleanSnrCurveTable(replayData.cleanSsMsCompareTable);
end
if ~isfield(replayData, 'cleanEstimateCrbCurveTable')
  replayData.cleanEstimateCrbCurveTable = localBuildCleanEstimateCrbCurveTable( ...
    replayData.cleanTrimAggregateTable, replayData.cleanHealthAggregateTable, replayData.cleanAggregateTable);
end
if ~isfield(replayData, 'cleanCrbPathAuditTable')
  replayData.cleanCrbPathAuditTable = table();
end
if ~isfield(replayData, 'cleanCrbPathAuditAggregateTable')
  replayData.cleanCrbPathAuditAggregateTable = localBuildCleanCrbPathAuditAggregateTable(replayData.cleanCrbPathAuditTable);
end
end

function taskList = localBuildTaskList(snrDbList, seedList)
%LOCALBUILDTASKLIST Build one task per SNR and seed pair.

taskList = repmat(struct('taskId', 0, 'snrDb', NaN, 'taskSeed', NaN), 0, 1);
taskId = 0;
for iSnr = 1:numel(snrDbList)
  for iSeed = 1:numel(seedList)
    taskId = taskId + 1;
    taskList(taskId, 1).taskId = taskId; %#ok<AGROW>
    taskList(taskId, 1).snrDb = snrDbList(iSnr); %#ok<AGROW>
    taskList(taskId, 1).taskSeed = seedList(iSeed); %#ok<AGROW>
  end
end
end

function [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, replayName)
%LOCALRUNTASKBATCH Run clean tasks with optional checkpoint/resume.

numTask = numel(taskList);
runState = struct();
progressTracker = localCreateProgressTracker(numTask);
progressCleanup = onCleanup(@() localCloseProgressTracker(progressTracker)); %#ok<NASGU>
try
  if logical(localGetFieldOrDefault(config, 'checkpointEnable', false))
    checkpointOpt = localBuildCheckpointOpt(replayName, config, localGetFieldOrDefault(config, 'useParfor', false));
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
  error('replayMfMsMleCrbCleanTrim:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for this truth-centered clean replay.');
end

fdHalfWidthHz = config.truthLocalFdHalfToothFraction * toothStepHz;
fdRangeRequested = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeUse = truthFdRateHzPerSec + config.truthLocalFdRateHalfWidthHzPerSec * [-1, 1];

noiseVar = 1 / (10^(task.snrDb / 10));
tStage = tic;
crbBundle = localBuildCrbBundleQuiet(fixture, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, noiseVar);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "crb-bundle", "crb", 1000 * toc(tStage)); %#ok<AGROW>
staticMethodNameList = config.methodNameList(endsWith(config.methodNameList, "SF-Static"));
dynamicMethodNameList = config.methodNameList(contains(config.methodNameList, "-MF-"));
fdRangeStaticUse = localApplyFdRefCrbScaledFloor(fdRangeRequested, truthFdRefHz, ...
  crbBundle, staticMethodNameList, config.truthLocalFdRefCrbScale);
fdRangeMfUse = localApplyFdRefCrbScaledFloor(fdRangeRequested, truthFdRefHz, ...
  crbBundle, dynamicMethodNameList, config.truthLocalFdRefCrbScale);

tStage = tic;
staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, fdRangeStaticUse, config);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "static-bundle", "static", 1000 * toc(tStage)); %#ok<AGROW>

caseRowList = repmat(localEmptyCaseRow(), numel(config.methodNameList), 1);
caseMap = struct();
for iMethod = 1:numel(config.methodNameList)
  methodName = config.methodNameList(iMethod);
  tMethod = tic;
  lastwarn('', '');
  if methodName == "SS-SF-Static"
    caseUse = staticBundle.caseStaticRefOnly;
    initMeta = staticBundle.initMetaStaticRefOnly;
  elseif methodName == "MS-SF-Static"
    caseUse = staticBundle.caseStaticMs;
    initMeta = staticBundle.initMetaStaticMs;
  else
    method = localBuildDynamicMethod(methodName, config.truthLocalDoaCrbScale, config.dynamicRouteMode);
    method = localApplyDoaCrbScaledWidth(method, crbBundle);
    staticSeedCase = localSelectStaticSeedCase(method, staticBundle);
    [caseUse, initMeta, methodRuntimeRows] = localRunDynamicMethod(method, fixture, truth, context, flowOpt, ...
      fdRangeMfUse, fdRateRangeUse, toothStepHz, config.optVerbose, staticSeedCase, config, task);
    runtimeRows = [runtimeRows; methodRuntimeRows(:)]; %#ok<AGROW>
  end
  if endsWith(methodName, "SF-Static")
    fdRangeForRow = fdRangeStaticUse;
  else
    fdRangeForRow = fdRangeMfUse;
  end
  row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
    fdRangeForRow, fdRateRangeUse, crbBundle, config, initMeta);
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
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.toothStepHz = toothStepHz;
repeatOut.fdRangeUse = fdRangeMfUse;
repeatOut.fdRangeRequested = fdRangeRequested;
repeatOut.fdRangeStaticUse = fdRangeStaticUse;
repeatOut.fdRangeMfUse = fdRangeMfUse;
repeatOut.fdRateRangeUse = fdRateRangeUse;
repeatOut.caseRowList = caseRowList;
repeatOut.crbPathAuditRowList = localBuildCleanCrbPathAuditRows(task, fixture, truth, crbBundle);
repeatOut.runtimeRowList = runtimeRows;
repeatOut.repeatTotalMs = 1000 * toc(tRepeat);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, fdRangeUse, config)
%LOCALBUILDCLEANSTATICBUNDLE Run truth-centered SS/MS static seeds for MF init.

truthDoaParam = reshape(truth.latlonTrueDeg, [], 1);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2)))
  error('replayMfMsMleCrbCleanTrim:MissingTruthDoa', ...
    'Truth lat/lon is required to build the clean truth-centered static search box.');
end
truthDoaParam = truthDoaParam(1:2);

staticBundle = struct();
[staticBundle.caseStaticRefOnly, staticBundle.initMetaStaticRefOnly] = localRunCleanStaticCase( ...
  "SS-SF-Static", "single", fixture.viewRefOnly, context, flowOpt.staticBaseOpt, ...
  crbBundle, truthDoaParam, fdRangeUse, config);
[staticBundle.caseStaticMs, staticBundle.initMetaStaticMs] = localRunCleanStaticCase( ...
  "MS-SF-Static", "multi", fixture.viewMs, context, flowOpt.staticBaseOpt, ...
  crbBundle, truthDoaParam, fdRangeUse, config);
staticBundle.bestStaticMsCase = staticBundle.caseStaticMs;
end

function [caseInfo, initMeta] = localRunCleanStaticCase(displayName, satMode, viewUse, context, staticBaseOpt, ...
  crbBundle, truthDoaParam, fdRangeUse, config)
%LOCALRUNCLEANSTATICCASE Run one SF static estimator inside the CRB-scaled truth-centered box.

modelOpt = staticBaseOpt;
modelOpt.initDoaParam = truthDoaParam(:);
modelOpt.initDoaHalfWidth = localResolveDoaCrbScaledHalfWidth(displayName, config.truthLocalDoaCrbScale, crbBundle);
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

function mode = localNormalizeDynamicObjectiveMode(mode)
%LOCALNORMALIZEDYNAMICOBJECTIVEMODE Normalize the dynamic objective dialect tag.

mode = lower(strtrim(string(mode)));
switch mode
  case {"pure-mle", "puremle", "pure"}
    mode = "pure-mle";
  case {"flow-default", "flowdefault", "robust", "engineering"}
    mode = "flow-default";
  otherwise
    error('replayMfMsMleCrbCleanTrim:InvalidDynamicObjectiveMode', ...
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
    error('replayMfMsMleCrbCleanTrim:InvalidDynamicObjectiveMode', ...
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
    error('replayMfMsMleCrbCleanTrim:InvalidDynamicRouteMode', ...
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
    error('replayMfMsMleCrbCleanTrim:InvalidDynamicRouteMode', ...
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

switch baseMethodName
  case {"SS-MF-CP-K", "MS-MF-CP-K"}
    fdRateMode = "known";
    isKnownRate = true;
  case {"SS-MF-CP-U", "MS-MF-CP-U"}
    fdRateMode = "unknown";
    isKnownRate = false;
  otherwise
    error('replayMfMsMleCrbCleanTrim:UnknownMethod', 'Unsupported method name "%s".', char(methodName));
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
method.phaseMode = "continuous";
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
  staticSeedCase = staticBundle.caseStaticMs;
else
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
  error('replayMfMsMleCrbCleanTrim:MissingTruthDoa', ...
    'Truth DoA is required to build the fixed clean MF search box.');
end
truthDoaParam = truthDoaParam(1:2);

estResult = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
seedDoaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(seedDoaParam) < 2 || any(~isfinite(seedDoaParam(1:2)))
  error('replayMfMsMleCrbCleanTrim:MissingStaticDoaSeed', ...
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
  fdRangeUse, fdRateRangeUse, crbBundle, config, initMeta)
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
  row.phaseMode = "continuous";
  row.fdRateMode = string(ternary(endsWith(baseMethodName, "CP-K"), "known", "unknown"));
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
isUnknown = contains(caseTable.displayName, "CP-U") & isDynamic;

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

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'ssKeepRate', NaN, 'msKeepRate', NaN, ...
  'ssAngleRmseDeg', NaN, 'msAngleRmseDeg', NaN, 'msAbsGainOverSs', NaN, ...
  'ssAngleMseOverOwnCrb', NaN, 'msAngleMseOverOwnCrb', NaN, 'msAngleMseOverSsCrb', NaN, ...
  'ssFdRefMseOverOwnCrb', NaN, 'msFdRefMseOverOwnCrb', NaN, ...
  'ssFdRateRmseHzPerSec', NaN, 'msFdRateRmseHzPerSec', NaN, ...
  'msBeatsSsAbs', false, 'msOwnCrbTargetPass', false, 'targetClass', "");
filterNameList = ["raw"; "health"; "joint-trim"];
snrList = unique(caseTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iFilter = 1:numel(filterNameList)
    filterName = filterNameList(iFilter);
    idx = idx + 1;
    row = rowTemplate;
    row.filterName = filterName;
    row.snrDb = snrDb;
    ssRaw = caseTable.displayName == "SS-MF-CP-U" & caseTable.snrDb == snrDb;
    msRaw = caseTable.displayName == "MS-MF-CP-U" & caseTable.snrDb == snrDb;
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
compareTable = struct2table(rowList(:));
end

function compareTable = localBuildCleanKnownUnknownCompareTable(caseTable)
%LOCALBUILDCLEANKNOWNUNKNOWNCOMPARETABLE Compare CP-K and CP-U clean MF rows.

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'satMode', "", ...
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
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList) * numel(satPrefixList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iSat = 1:numel(satPrefixList)
    satPrefix = satPrefixList(iSat);
    knownName = satPrefix + "-MF-CP-K";
    unknownName = satPrefix + "-MF-CP-U";
    knownRaw = caseTable.displayName == knownName & caseTable.snrDb == snrDb;
    unknownRaw = caseTable.displayName == unknownName & caseTable.snrDb == snrDb;
    for iFilter = 1:numel(filterNameList)
      filterName = filterNameList(iFilter);
      idx = idx + 1;
      row = rowTemplate;
      row.filterName = filterName;
      row.snrDb = snrDb;
      row.satMode = string(ternary(satPrefix == "MS", "multi", "single"));
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

rowTemplate = struct('snrDb', NaN, 'ssAngleMseOverCrb', NaN, 'msAngleMseOverCrb', NaN, ...
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
curveSource = sortrows(compareTable(mask, :), 'snrDb');
rowList = repmat(rowTemplate, height(curveSource), 1);
for iRow = 1:height(curveSource)
  row = rowTemplate;
  row.snrDb = curveSource.snrDb(iRow);
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
methodNameList = ["SS-SF-Static"; "MS-SF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "MS-MF-CP-K"; "MS-MF-CP-U"];
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

function localPlotCleanEstimateVsCrbCurves(curveTable)
%LOCALPLOTCLEANESTIMATEVSCRBCURVES Plot static/MF estimates and CRB on log axes.

if isempty(curveTable) || ~istable(curveTable) || height(curveTable) == 0
  return;
end
methodNameList = ["SS-SF-Static"; "MS-SF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "MS-MF-CP-K"; "MS-MF-CP-U"];
methodLabelList = ["SS Static"; "MS Static"; "SS-MF-K"; "SS-MF-U"; "MS-MF-K"; "MS-MF-U"];
markerList = {'o', 's', '^', 'v', 'd', 'p'};
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
  warning('replayMfMsMleCrbCleanTrim:PlotFailed', ...
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
    error('replayMfMsMleCrbCleanTrim:UnknownFilter', 'Unsupported filter "%s".', char(filterName));
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


function auditTable = localBuildCleanCrbPathAuditTable(repeatCell)
%LOCALBUILDCLEANCRBPATHAUDITTABLE Collect K/U CRB path audit rows from tasks.

rowList = repmat(localEmptyCrbPathAuditRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isempty(repeatOut) || ~isfield(repeatOut, 'crbPathAuditRowList')
    continue;
  end
  rowList = [rowList; repeatOut.crbPathAuditRowList(:)]; %#ok<AGROW>
end
auditTable = struct2table(rowList(:));
end

function aggregateTable = localBuildCleanCrbPathAuditAggregateTable(auditTable)
%LOCALBUILDCLEANCRBPATHAUDITAGGREGATETABLE Summarize K/U CRB path health.

rowTemplate = struct('satMode', "", 'numRow', NaN, ...
  'crbMonotonicPassRate', NaN, 'efimMonotonicPassRate', NaN, ...
  'lossPsdPassRate', NaN, 'sharedBlockPassRate', NaN, ...
  'fdRefInvariantPassRate', NaN, 'fdRateInvariantPassRate', NaN, ...
  'gammaCouplingFiniteRate', NaN, 'angleCrbUnknownOverKnownMedian', NaN, ...
  'fdRefCrbUnknownOverKnownMedian', NaN, 'lossTraceRatioMedian', NaN, ...
  'gammaCouplingMedian', NaN, 'sharedBlockRelErrMax', NaN, ...
  'minEigEfimKnownMinusUnknownMin', NaN, 'minEigCrbUnknownMinusKnownMin', NaN, ...
  'primaryAuditClass', "");
if isempty(auditTable) || ~istable(auditTable) || height(auditTable) == 0
  aggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
satModeList = unique(auditTable.satMode, 'stable');
rowList = repmat(rowTemplate, numel(satModeList), 1);
for iMode = 1:numel(satModeList)
  mask = auditTable.satMode == satModeList(iMode);
  row = rowTemplate;
  row.satMode = satModeList(iMode);
  row.numRow = nnz(mask);
  row.crbMonotonicPassRate = mean(double(auditTable.crbMonotonicOk(mask)), 'omitnan');
  row.efimMonotonicPassRate = mean(double(auditTable.efimMonotonicOk(mask)), 'omitnan');
  row.lossPsdPassRate = mean(double(auditTable.lossPsdOk(mask)), 'omitnan');
  row.sharedBlockPassRate = mean(double(auditTable.sharedBlockOk(mask)), 'omitnan');
  row.fdRefInvariantPassRate = mean(double(auditTable.fdRefInvariantOk(mask)), 'omitnan');
  row.fdRateInvariantPassRate = mean(double(auditTable.fdRateInvariantOk(mask)), 'omitnan');
  row.gammaCouplingFiniteRate = mean(double(isfinite(auditTable.gammaCoupling(mask)) & auditTable.gammaCoupling(mask) > 0), 'omitnan');
  row.angleCrbUnknownOverKnownMedian = median(auditTable.angleCrbUnknownOverKnown(mask), 'omitnan');
  row.fdRefCrbUnknownOverKnownMedian = median(auditTable.fdRefCrbUnknownOverKnown(mask), 'omitnan');
  row.lossTraceRatioMedian = median(auditTable.lossTraceRatio(mask), 'omitnan');
  row.gammaCouplingMedian = median(auditTable.gammaCoupling(mask), 'omitnan');
  row.sharedBlockRelErrMax = max(auditTable.sharedBlockRelErr(mask), [], 'omitnan');
  row.minEigEfimKnownMinusUnknownMin = min(auditTable.minEigEfimKnownMinusUnknown(mask), [], 'omitnan');
  row.minEigCrbUnknownMinusKnownMin = min(auditTable.minEigCrbUnknownMinusKnown(mask), [], 'omitnan');
  row.primaryAuditClass = localClassifyCrbAuditAggregate(row);
  rowList(iMode) = row;
end
aggregateTable = struct2table(rowList(:));
end

function classText = localClassifyCrbAuditAggregate(row)
%LOCALCLASSIFYCRBAUDITAGGREGATE Classify compact CRB path audit outcome.

if ~(isfinite(row.crbMonotonicPassRate) && isfinite(row.efimMonotonicPassRate) && ...
    isfinite(row.lossPsdPassRate) && isfinite(row.sharedBlockPassRate) && ...
    isfinite(row.fdRefInvariantPassRate) && isfinite(row.fdRateInvariantPassRate) && ...
    isfinite(row.gammaCouplingFiniteRate))
  classText = "audit-missing";
elseif row.crbMonotonicPassRate < 1
  classText = "crb-monotonic-fail";
elseif row.efimMonotonicPassRate < 1
  classText = "efim-monotonic-fail";
elseif row.lossPsdPassRate < 1
  classText = "schur-loss-fail";
elseif row.sharedBlockPassRate < 1
  classText = "shared-block-mismatch";
elseif row.fdRefInvariantPassRate < 1 || row.fdRateInvariantPassRate < 1
  classText = "reference-state-invariant-fail";
elseif row.gammaCouplingFiniteRate < 1
  classText = "gamma-coupling-weak-or-missing";
else
  classText = "path-ok";
end
end

function rowList = localBuildCleanCrbPathAuditRows(task, fixture, truth, crbBundle)
%LOCALBUILDCLEANCRBPATHAUDITROWS Build SS/MS K/U CRB path audit rows.

rowList = repmat(localEmptyCrbPathAuditRow(), 2, 1);
rowList(1) = localBuildOneCrbPathAuditRow(task, fixture.sceneSeqRefOnly, truth, ...
  "single", "SS-MF-CP-K", "SS-MF-CP-U", ...
  crbBundle.crbMfRefKnown, crbBundle.crbMfRefUnknown, ...
  crbBundle.auxCrbMfRefKnown, crbBundle.auxCrbMfRefUnknown);
rowList(2) = localBuildOneCrbPathAuditRow(task, fixture.sceneSeq, truth, ...
  "multi", "MS-MF-CP-K", "MS-MF-CP-U", ...
  crbBundle.crbMfMsKnown, crbBundle.crbMfMsUnknown, ...
  crbBundle.auxCrbMfMsKnown, crbBundle.auxCrbMfMsUnknown);
end

function row = localBuildOneCrbPathAuditRow(task, sceneSeq, truth, satMode, knownName, unknownName, crbKnown, crbUnknown, auxKnown, auxUnknown)
%LOCALBUILDONECRBPATHAUDITROW Audit one known-rate / unknown-rate CRB pair.

row = localEmptyCrbPathAuditRow();
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.satMode = string(satMode);
row.knownName = string(knownName);
row.unknownName = string(unknownName);
row.refSatIdxLocal = localResolveRefSatIdxLocal(sceneSeq);
row.gammaParamIdx = localFindParamNameIndex(auxUnknown, "fdRate");

[angleCrbKnownDeg, fdRefCrbKnownHz] = localExtractCrbStdMetrics(crbKnown, truth);
[angleCrbUnknownDeg, fdRefCrbUnknownHz] = localExtractCrbStdMetrics(crbUnknown, truth);
row.angleCrbKnownDeg = angleCrbKnownDeg;
row.angleCrbUnknownDeg = angleCrbUnknownDeg;
row.angleCrbUnknownOverKnown = localSafeRatio(angleCrbUnknownDeg, angleCrbKnownDeg);
row.fdRefCrbKnownHz = fdRefCrbKnownHz;
row.fdRefCrbUnknownHz = fdRefCrbUnknownHz;
row.fdRefCrbUnknownOverKnown = localSafeRatio(fdRefCrbUnknownHz, fdRefCrbKnownHz);

fimKnown = localGetFiniteMatrixField(auxKnown, 'fimInterest', 3, 3);
fimUnknown = localGetFiniteMatrixField(auxUnknown, 'fimInterest', 3, 3);
crbKnown = localEnsureMatrixSize(crbKnown, 3, 3);
crbUnknown = localEnsureMatrixSize(crbUnknown, 3, 3);
row.minEigEfimKnownMinusUnknown = localMinSymEig(fimKnown - fimUnknown);
row.minEigCrbUnknownMinusKnown = localMinSymEig(crbUnknown - crbKnown);

schurAudit = localBuildGammaSchurAudit(auxUnknown, fimKnown);
row.lossMinEig = schurAudit.lossMinEig;
row.lossTraceRatio = schurAudit.lossTraceRatio;
row.gammaCoupling = schurAudit.gammaCoupling;
row.gammaInfo = schurAudit.gammaInfo;
row.sharedBlockRelErr = schurAudit.sharedBlockRelErr;

fdInvKnown = localCheckCrbReferenceStateInvariant(auxKnown, row.refSatIdxLocal, truth);
fdInvUnknown = localCheckCrbReferenceStateInvariant(auxUnknown, row.refSatIdxLocal, truth);
row.fdSatRefErrKnownHz = fdInvKnown.fdSatRefErrHz;
row.fdSatRefErrUnknownHz = fdInvUnknown.fdSatRefErrHz;
row.fdRateRefErrKnownHzPerSec = fdInvKnown.fdRateRefErrHzPerSec;
row.fdRateRefErrUnknownHzPerSec = fdInvUnknown.fdRateRefErrHzPerSec;
row.fdRefInvariantOk = fdInvKnown.fdRefInvariantOk && fdInvUnknown.fdRefInvariantOk;
row.fdRateInvariantOk = fdInvKnown.fdRateInvariantOk && fdInvUnknown.fdRateInvariantOk;

relEigTol = 1e-8;
row.crbMonotonicOk = localIsNonnegativeWithTol(row.minEigCrbUnknownMinusKnown, crbKnown, relEigTol);
row.efimMonotonicOk = localIsNonnegativeWithTol(row.minEigEfimKnownMinusUnknown, fimKnown, relEigTol);
row.lossPsdOk = localIsNonnegativeWithTol(row.lossMinEig, fimKnown, relEigTol);
row.sharedBlockOk = isfinite(row.sharedBlockRelErr) && row.sharedBlockRelErr <= 1e-6;
row.auditClass = localClassifyCrbAuditRow(row);
end

function row = localEmptyCrbPathAuditRow()
%LOCALEMPTYCRBPATHAUDITROW Return a typed empty CRB path audit row.

row = struct('satMode', "", 'knownName', "", 'unknownName', "", ...
  'snrDb', NaN, 'taskSeed', NaN, 'refSatIdxLocal', NaN, 'gammaParamIdx', NaN, ...
  'angleCrbKnownDeg', NaN, 'angleCrbUnknownDeg', NaN, ...
  'angleCrbUnknownOverKnown', NaN, 'fdRefCrbKnownHz', NaN, ...
  'fdRefCrbUnknownHz', NaN, 'fdRefCrbUnknownOverKnown', NaN, ...
  'minEigEfimKnownMinusUnknown', NaN, 'minEigCrbUnknownMinusKnown', NaN, ...
  'lossMinEig', NaN, 'lossTraceRatio', NaN, 'gammaCoupling', NaN, ...
  'gammaInfo', NaN, 'sharedBlockRelErr', NaN, ...
  'fdSatRefErrKnownHz', NaN, 'fdSatRefErrUnknownHz', NaN, ...
  'fdRateRefErrKnownHzPerSec', NaN, 'fdRateRefErrUnknownHzPerSec', NaN, ...
  'crbMonotonicOk', false, 'efimMonotonicOk', false, 'lossPsdOk', false, ...
  'sharedBlockOk', false, 'fdRefInvariantOk', false, 'fdRateInvariantOk', false, ...
  'auditClass', "");
end

function [angleCrbStdDeg, fdRefCrbStdHz] = localExtractCrbStdMetrics(crbFull, truth)
%LOCALEXTRACTCRBSTDMETRICS Extract spherical angle and fdRef CRB std.

crbFull = localEnsureMatrixSize(crbFull, 3, 3);
angleCrbStdDeg = NaN;
fdRefCrbStdHz = NaN;
if all(isfinite(crbFull(:)))
  [angleCrbStdDeg, ~] = projectCrbToAngleMetric(real(crbFull(1:2, 1:2)), truth.latlonTrueDeg, 'latlon');
  fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
end
end

function schurAudit = localBuildGammaSchurAudit(auxUnknown, fimKnown)
%LOCALBUILDGAMMASCHURAUDIT Rebuild gamma Schur loss after common nuisance removal.

schurAudit = struct('lossMinEig', NaN, 'lossTraceRatio', NaN, ...
  'gammaCoupling', NaN, 'gammaInfo', NaN, 'sharedBlockRelErr', NaN);
fimFull = localGetFieldOrDefault(auxUnknown, 'fimFull', []);
gammaIdx = localFindParamNameIndex(auxUnknown, "fdRate");
if isempty(fimFull) || ~isnumeric(fimFull) || ~isfinite(gammaIdx) || gammaIdx < 4
  return;
end
fimFull = real(fimFull);
numParam = size(fimFull, 1);
if size(fimFull, 2) ~= numParam || numParam < gammaIdx || numParam < 4
  return;
end
interestGammaIdx = [1, 2, 3, gammaIdx];
commonNuisanceIdx = setdiff(1:numParam, interestGammaIdx);
reducedFim = localSchurReduce(fimFull, interestGammaIdx, commonNuisanceIdx);
if size(reducedFim, 1) ~= 4 || any(~isfinite(reducedFim(:)))
  return;
end
etaEta = reducedFim(1:3, 1:3);
etaGamma = reducedFim(1:3, 4);
gammaGamma = reducedFim(4, 4);
if isfinite(gammaGamma) && abs(gammaGamma) > eps(max(1, norm(reducedFim, 'fro')))
  lossMat = (etaGamma * etaGamma') ./ gammaGamma;
  schurAudit.lossMinEig = localMinSymEig(lossMat);
  schurAudit.lossTraceRatio = localSafeRatio(real(trace(lossMat)), real(trace(etaEta)));
  schurAudit.gammaCoupling = norm(etaGamma) ./ sqrt(max(real(trace(etaEta)), 0) * abs(gammaGamma));
  schurAudit.gammaInfo = gammaGamma;
end
if all(isfinite(fimKnown(:))) && any(fimKnown(:) ~= 0)
  schurAudit.sharedBlockRelErr = norm(etaEta - fimKnown, 'fro') ./ max(norm(fimKnown, 'fro'), eps);
end
end

function reducedFim = localSchurReduce(fimFull, keepIdx, nuisanceIdx)
%LOCALSCHURREDUCE Eliminate nuisance indices from a symmetric FIM.

fimFull = localSymPart(real(fimFull));
keepIdx = reshape(double(keepIdx), 1, []);
nuisanceIdx = reshape(double(nuisanceIdx), 1, []);
if isempty(nuisanceIdx)
  reducedFim = fimFull(keepIdx, keepIdx);
  return;
end
fimKeepKeep = fimFull(keepIdx, keepIdx);
fimKeepNuisance = fimFull(keepIdx, nuisanceIdx);
fimNuisanceNuisance = fimFull(nuisanceIdx, nuisanceIdx);
reducedFim = fimKeepKeep - fimKeepNuisance * (pinv(fimNuisanceNuisance) * fimKeepNuisance');
reducedFim = localSymPart(reducedFim);
end

function invInfo = localCheckCrbReferenceStateInvariant(aux, refSatIdxLocal, truth)
%LOCALCHECKCRBREFERENCESTATEINVARIANT Check CRB-side reference Doppler state.

invInfo = struct('fdSatRefErrHz', NaN, 'fdRateRefErrHzPerSec', NaN, ...
  'fdRefInvariantOk', false, 'fdRateInvariantOk', false);
if ~isfinite(refSatIdxLocal) || refSatIdxLocal < 1
  return;
end
refIdx = round(refSatIdxLocal);
fdSat = reshape(localGetFieldOrDefault(aux, 'fdSat', []), [], 1);
fdRateSat = reshape(localGetFieldOrDefault(aux, 'fdRateSat', []), [], 1);
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if numel(fdSat) >= refIdx && isfinite(fdSat(refIdx)) && isfinite(truthFdRefHz)
  invInfo.fdSatRefErrHz = fdSat(refIdx) - truthFdRefHz;
  fdTol = max(1e-5, 1e-10 * max(1, abs(truthFdRefHz)));
  invInfo.fdRefInvariantOk = abs(invInfo.fdSatRefErrHz) <= fdTol;
end
if numel(fdRateSat) >= refIdx && isfinite(fdRateSat(refIdx)) && isfinite(truthFdRateHzPerSec)
  invInfo.fdRateRefErrHzPerSec = fdRateSat(refIdx) - truthFdRateHzPerSec;
  fdRateTol = max(1e-5, 1e-10 * max(1, abs(truthFdRateHzPerSec)));
  invInfo.fdRateInvariantOk = abs(invInfo.fdRateRefErrHzPerSec) <= fdRateTol;
end
end

function refSatIdxLocal = localResolveRefSatIdxLocal(sceneSeq)
%LOCALRESOLVEREFSATIDXLOCAL Resolve local reference-satellite index for an audit row.

refSatIdxLocal = NaN;
if isstruct(sceneSeq) && isfield(sceneSeq, 'refSatIdxLocal') && ~isempty(sceneSeq.refSatIdxLocal)
  refSatIdxLocal = double(sceneSeq.refSatIdxLocal(1));
elseif isstruct(sceneSeq) && isfield(sceneSeq, 'ref') && isfield(sceneSeq.ref, 'satIdxLocal') && ~isempty(sceneSeq.ref.satIdxLocal)
  refSatIdxLocal = double(sceneSeq.ref.satIdxLocal(1));
elseif isstruct(sceneSeq) && isfield(sceneSeq, 'sceneCell') && ~isempty(sceneSeq.sceneCell)
  refFrameIdx = round(localGetFieldOrDefault(sceneSeq, 'refFrameIdx', 1));
  refFrameIdx = min(max(refFrameIdx, 1), numel(sceneSeq.sceneCell));
  sceneRef = sceneSeq.sceneCell{refFrameIdx};
  if isstruct(sceneRef) && isfield(sceneRef, 'ref') && isfield(sceneRef.ref, 'satIdxLocal') && ~isempty(sceneRef.ref.satIdxLocal)
    refSatIdxLocal = double(sceneRef.ref.satIdxLocal(1));
  end
end
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && refSatIdxLocal >= 1)
  refSatIdxLocal = 1;
end
end

function idx = localFindParamNameIndex(aux, paramName)
%LOCALFINDPARAMNAMEINDEX Find a full-parameter index by exact or contains match.

idx = NaN;
nameList = string(localGetFieldOrDefault(aux, 'paramNameFull', strings(0, 1)));
if isempty(nameList)
  return;
end
matchIdx = find(nameList == string(paramName), 1, 'first');
if isempty(matchIdx)
  matchIdx = find(contains(lower(nameList), lower(string(paramName))), 1, 'first');
end
if ~isempty(matchIdx)
  idx = matchIdx;
end
end

function matrixValue = localGetFiniteMatrixField(dataStruct, fieldName, numRow, numCol)
%LOCALGETFINITEMATRIXFIELD Read one finite matrix field or return NaNs.

matrixValue = NaN(numRow, numCol);
value = localGetFieldOrDefault(dataStruct, fieldName, []);
if isnumeric(value) && isequal(size(value), [numRow, numCol])
  matrixValue = real(value);
end
end

function matrixValue = localEnsureMatrixSize(matrixValue, numRow, numCol)
%LOCALENSUREMATRIXSIZE Return a matrix with the requested size or NaNs.

if ~(isnumeric(matrixValue) && isequal(size(matrixValue), [numRow, numCol]))
  matrixValue = NaN(numRow, numCol);
else
  matrixValue = real(matrixValue);
end
end

function minEig = localMinSymEig(matrixValue)
%LOCALMINSYMEIG Return the smallest eigenvalue of a symmetric matrix.

if isempty(matrixValue) || any(~isfinite(matrixValue(:))) || size(matrixValue, 1) ~= size(matrixValue, 2)
  minEig = NaN;
  return;
end
minEig = min(real(eig(localSymPart(matrixValue))));
end

function matrixValue = localSymPart(matrixValue)
%LOCALSYMPART Force a matrix to its real symmetric part.

matrixValue = 0.5 * (real(matrixValue) + real(matrixValue)');
end

function tf = localIsNonnegativeWithTol(minEig, scaleMatrix, relTol)
%LOCALISNONNEGATIVEWITHTOL Check PSD with a relative numerical tolerance.

if ~(isfinite(minEig) && isnumeric(scaleMatrix))
  tf = false;
  return;
end
scaleValue = max(1, norm(real(scaleMatrix), 'fro'));
tf = minEig >= -relTol * scaleValue;
end

function classText = localClassifyCrbAuditRow(row)
%LOCALCLASSIFYCRBAUDITROW Classify one K/U CRB path audit row.

if ~(isfinite(row.minEigCrbUnknownMinusKnown) && isfinite(row.minEigEfimKnownMinusUnknown) && ...
    isfinite(row.lossMinEig) && isfinite(row.sharedBlockRelErr))
  classText = "audit-missing";
elseif ~row.crbMonotonicOk
  classText = "crb-monotonic-fail";
elseif ~row.efimMonotonicOk
  classText = "efim-monotonic-fail";
elseif ~row.lossPsdOk
  classText = "schur-loss-fail";
elseif ~row.sharedBlockOk
  classText = "shared-block-mismatch";
elseif ~(row.fdRefInvariantOk && row.fdRateInvariantOk)
  classText = "reference-state-invariant-fail";
elseif ~(isfinite(row.gammaCoupling) && row.gammaCoupling > 0)
  classText = "gamma-coupling-weak-or-missing";
else
  classText = "path-ok";
end
end

function outlierTable = localBuildCleanOutlierTable(caseTable)
%LOCALBUILDCLEANOUTLIERTABLE Keep the MS-MF-CP-U rows rejected by clean trim.

if isempty(caseTable) || height(caseTable) == 0
  outlierTable = table();
  return;
end
mask = caseTable.displayName == "MS-MF-CP-U" & ~caseTable.jointTrimKeep;
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
%LOCALPRINTREPLAYCONFIG Print clean replay axes.

fprintf('  %-32s : %s\n', 'SNR list (dB)', localFormatRow(config.snrDbList));
fprintf('  %-32s : %d\n', 'frame count', config.numFrame);
fprintf('  %-32s : %s\n', 'method list', strjoin(cellstr(config.methodNameList), ', '));
fprintf('  %-32s : %.3f tooth\n', 'fdRef requested half-width', config.truthLocalFdHalfToothFraction);
fprintf('  %-32s : %.2f x CRB\n', 'static fdRef CRB floor scale', config.truthLocalFdRefCrbScale);
fprintf('  %-32s : %.2f Hz/s\n', 'fdRate truth-local half-width', config.truthLocalFdRateHalfWidthHzPerSec);
fprintf('  %-32s : %.2f x CRB\n', 'DoA CRB half-width scale', config.truthLocalDoaCrbScale);
fprintf('  %-32s : %s\n', 'estimation chain', char(config.initMode));
fprintf('  %-32s : %s / %s\n', 'MF objective / route', ...
  char(config.dynamicObjectiveMode), char(config.dynamicRouteMode));
fprintf('  %-32s : %.2f / %.2f CRB\n', 'trim angle/fdRef cap', config.trimAngleNormMax, config.trimFdRefNormMax);
fprintf('  %-32s : %d / %d\n', 'task count / outer parfor', config.numTask, logical(config.useParfor));
end

function localPrintCompactTable(sectionTitle, dataTable, edgeCount)
%LOCALPRINTCOMPACTTABLE Print compact table preview.

if nargin < 3 || isempty(edgeCount)
  edgeCount = 12;
end
if isempty(dataTable) || height(dataTable) <= edgeCount
  printMfReplaySection(sectionTitle, dataTable);
  return;
end
printMfReplaySection(sectionTitle, table(height(dataTable), 'VariableNames', {'numRows'}));
fprintf('%s preview:\n', sectionTitle);
dispMfReplayTablePreview(dataTable, edgeCount);
end

function lineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build short notification metrics.

lineList = strings(0, 1);
compareTable = replayData.cleanSsMsCompareTable;
knownUnknownTable = localGetFieldOrDefault(replayData, 'cleanKnownUnknownCompareTable', table());
if istable(knownUnknownTable) && height(knownUnknownTable) > 0
  trimKuRow = knownUnknownTable(knownUnknownTable.filterName == "joint-trim" & knownUnknownTable.satMode == "multi", :);
  if ~isempty(trimKuRow)
    row = trimKuRow(1, :);
    lineList(end + 1, 1) = sprintf('MS K/U: dMSE/CRB=%.3g, U gain=%.3g, class=%s', ...
      row.unknownMinusKnownAngleMseOverCrb, row.unknownRmseGainOverKnown, char(row.comparisonClass)); %#ok<AGROW>
  end
end
crbAuditAggregate = localGetFieldOrDefault(replayData, 'cleanCrbPathAuditAggregateTable', table());
if istable(crbAuditAggregate) && height(crbAuditAggregate) > 0
  msAuditRow = crbAuditAggregate(crbAuditAggregate.satMode == "multi", :);
  if ~isempty(msAuditRow)
    row = msAuditRow(1, :);
    lineList(end + 1, 1) = sprintf('MS CRB audit: mono=%.2f, gamma=%.2f, shared=%.2f, class=%s', ...
      row.crbMonotonicPassRate, row.gammaCouplingFiniteRate, row.sharedBlockPassRate, char(row.primaryAuditClass)); %#ok<AGROW>
  end
end
compareTable = replayData.cleanSsMsCompareTable;
if istable(compareTable) && height(compareTable) > 0
  trimRow = compareTable(compareTable.filterName == "joint-trim", :);
  if ~isempty(trimRow)
    row = trimRow(1, :);
    lineList(end + 1, 1) = sprintf('joint-trim: ms/ownCRB=%.3g, ms/ssCRB=%.3g, msAbsGain=%.3g, class=%s', ...
      row.msAngleMseOverOwnCrb, row.msAngleMseOverSsCrb, row.msAbsGainOverSs, char(row.targetClass)); %#ok<AGROW>
  end
end
lineList(end + 1, 1) = sprintf('elapsed %.1f min', replayData.elapsedSec / 60); %#ok<AGROW>
end

function checkpointOpt = localBuildCheckpointOpt(replayName, config, useParfor)
%LOCALBUILDCHECKPOINTOPT Build common checkpoint options with a stable run key.

meta = struct('snrDbList', reshape(config.snrDbList, 1, []), ...
  'seedList', reshape(config.seedList, 1, []), ...
  'dynamicObjectiveMode', char(config.dynamicObjectiveMode), ...
  'dynamicRouteMode', char(config.dynamicRouteMode), ...
  'routeTag', char(config.routeTag));
checkpointOpt = buildMfReplayCheckpointOpt(replayName, config, struct( ...
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
  sprintf('fdFrac%.8g', config.truthLocalFdHalfToothFraction), ...
  sprintf('fdRefCrbScale%.8g', config.truthLocalFdRefCrbScale), ...
  sprintf('fdRate%.8g', config.truthLocalFdRateHalfWidthHzPerSec), ...
  sprintf('doaCrbScale%.6g', config.truthLocalDoaCrbScale), ...
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
%LOCALCANUSEPARFOR Check whether the outer replay task loop can use parfor.

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
  error('replayMfMsMleCrbCleanTrim:InvalidFrameCount', ...
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
%LOCALBUILDCRBBUNDLEQUIET Build CRB while suppressing expected full-FIM warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
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
  otherwise
    error('replayMfMsMleCrbCleanTrim:UnknownCrbMethod', 'Unsupported CRB method "%s".', char(methodName));
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
  error('replayMfMsMleCrbCleanTrim:InvalidDoaCrbScale', ...
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

function fdRangeUse = localApplyFdRefCrbScaledFloor(fdRangeUse, truthFdRefHz, crbBundle, methodNameList, fdRefCrbScale)
%LOCALAPPLYFDREFCRBSCALEDFLOOR Ensure the fdRef box is not a one-sigma CRB truncation.

fdRefCrbScale = double(fdRefCrbScale);
if isempty(fdRefCrbScale) || ~isfinite(fdRefCrbScale(1)) || fdRefCrbScale(1) <= 0
  error('replayMfMsMleCrbCleanTrim:InvalidFdRefCrbScale', ...
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
  'dynamicRouteMode', "", 'snrDb', NaN, 'taskSeed', NaN, 'numFrame', NaN, ...
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
%LOCALCREATEPROGRESSTRACKER Create a best-effort replay progress tracker.

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
%LOCALCLOSEPROGRESSTRACKER Close a best-effort replay progress tracker.

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
