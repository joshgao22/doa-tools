% replayMfSsMleCrbMetricDiagnose
% Purpose: diagnose single-satellite SF/MF in-tooth MLE-vs-CRB metric gaps.
% This replay mirrors scanMfMleCrbInToothConsistency on a small seed/SNR set
% and prints estimator metric, CRB projection, requested/init/final, frequency
% lock, health/trim filtering, CP-U release, objective probes, solve probes,
% and DoA-path probes without changing the paper-facing scan path.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted task runs from per-task files under the repo-root tmp. The outer
% task loop uses parfor automatically when available; estimator-inner parallelism
% remains disabled for numerical-path consistency.
% saveSnapshot=true saves only replayData via saveExpSnapshot. Telegram notice is best-effort only.

clear; close all; clc;

%% Replay configuration

replayName = "replayMfSsMleCrbMetricDiagnose";
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = true;
optVerbose = false;

numFrame = 10;
snrDbList = -15:3:10;
% Append 15/20 dB only for numerical stress; the paper-facing satellite regime
% is low-to-medium SNR.
seedList = [266; 271; 347; 253; 260];
methodNameList = [
  "SS-SF-Static"
  "SS-MF-CP-K"
  "SS-MF-CP-U"
];

oracleFdHalfToothFraction = 0.4;
oracleFdRateHalfWidthHzPerSec = 1000;
staticLocalDoaHalfWidthDeg = [0.002; 0.002];
boundaryTolFraction = 0.02;
mfInitMode = "static-doa-truth-fd";
healthFdRefToothFraction = 0.25;
healthFdRateAbsMaxHzPerSec = 250;
trimAngleNormMax = 5;
trimFdRefNormMax = 5;
objectiveProbeEnable = true;
solveProbeEnable = true;
pathProbeEnable = true;
probeWideDoaHalfWidthDegList = 0.012;
probeTruthDoaHalfWidthDeg = [0.002; 0.002];
objectivePathAlphaList = [0; 0.25; 0.5; 0.75; 1];

seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
methodNameList = reshape(string(methodNameList), [], 1);
staticLocalDoaHalfWidthDeg = reshape(double(staticLocalDoaHalfWidthDeg), [], 1);
probeWideDoaHalfWidthDegList = reshape(double(probeWideDoaHalfWidthDegList), [], 1);
probeTruthDoaHalfWidthDeg = reshape(double(probeTruthDoaHalfWidthDeg), [], 1);
objectivePathAlphaList = reshape(double(objectivePathAlphaList), [], 1);

runTic = tic;
replayData = struct();
config = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and flow options

  config.replayName = string(replayName);
  config.numFrame = numFrame;
  config.snrDbList = snrDbList;
  config.seedList = seedList;
  config.baseSeed = seedList(1);
  config.numRepeat = numel(seedList);
  config.methodNameList = methodNameList;
  config.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
  config.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
  config.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg;
  config.boundaryTolFraction = boundaryTolFraction;
  config.mfInitMode = string(mfInitMode);
  config.healthFdRefToothFraction = healthFdRefToothFraction;
  config.healthFdRateAbsMaxHzPerSec = healthFdRateAbsMaxHzPerSec;
  config.trimAngleNormMax = trimAngleNormMax;
  config.trimFdRefNormMax = trimFdRefNormMax;
  config.objectiveProbeEnable = logical(objectiveProbeEnable);
  config.solveProbeEnable = logical(solveProbeEnable);
  config.pathProbeEnable = logical(pathProbeEnable);
  config.probeWideDoaHalfWidthDegList = probeWideDoaHalfWidthDegList;
  config.probeTruthDoaHalfWidthDeg = probeTruthDoaHalfWidthDeg;
  config.objectivePathAlphaList = objectivePathAlphaList;
  config.saveSnapshot = logical(saveSnapshot);
  config.checkpointEnable = logical(checkpointEnable);
  config.checkpointResume = logical(checkpointResume);
  config.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
  config.notifyTelegramEnable = logical(notifyTelegramEnable);
  config.optVerbose = logical(optVerbose);
  config.initMode = config.mfInitMode;
  taskList = localBuildTaskList(config.snrDbList, config.seedList);
  config.numTask = numel(taskList);
  config.useParfor = localCanUseParfor(config.numTask);
  config.runKey = localBuildCheckpointRunKey(config);
  checkpointRunDir = "";
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
  context = localDisableSubsetBankForDiagnose(context);
  flowOpt = localBuildFlowOpt(config.staticLocalDoaHalfWidthDeg, context.parallelOpt);

  %% Run replay batch

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
  aggregateTable = localBuildAggregateTable(caseTable);
  healthAggregateTable = localBuildFilteredAggregateTable(caseTable, "health", caseTable.healthResolved);
  trimAggregateTable = localBuildFilteredAggregateTable(caseTable, "trim", caseTable.trimKeep);
  tailTable = localBuildTailTable(caseTable);
  crbMetricTable = localBuildCrbMetricTable(caseTable);
  releaseCompareTable = localBuildReleaseCompareTable(caseTable);
  rangeTable = localBuildRangeTable(repeatCell);
  objectiveProbeTable = localBuildObjectiveProbeTable(repeatCell);
  solveProbeTable = localBuildSolveProbeTable(repeatCell);
  pathProbeTable = localBuildPathProbeTable(repeatCell);
  diagnosisTable = localBuildDiagnosisTable(caseTable, objectiveProbeTable, solveProbeTable);
  plotData = localBuildPlotData(caseTable, aggregateTable, releaseCompareTable);

  %% Data storage

  replayData = struct();
  replayData.replayName = string(replayName);
  replayData.runKey = config.runKey;
  replayData.utcRun = datetime('now', 'TimeZone', 'local');
  replayData.config = config;
  replayData.contextSummary = localBuildContextSummary(context);
  replayData.caseTable = caseTable;
  replayData.aggregateTable = aggregateTable;
  replayData.healthAggregateTable = healthAggregateTable;
  replayData.trimAggregateTable = trimAggregateTable;
  replayData.tailTable = tailTable;
  replayData.crbMetricTable = crbMetricTable;
  replayData.releaseCompareTable = releaseCompareTable;
  replayData.rangeTable = rangeTable;
  replayData.objectiveProbeTable = objectiveProbeTable;
  replayData.solveProbeTable = solveProbeTable;
  replayData.pathProbeTable = pathProbeTable;
  replayData.diagnosisTable = diagnosisTable;
  replayData.repeatCell = localStripRepeatCell(repeatCell);
  if config.checkpointEnable
    replayData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  end
  replayData.plotData = plotData;
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

  %% Summary output and plotting

  printMfReplaySection('SS MLE/CRB diagnostic raw aggregate', replayData.aggregateTable);
  printMfReplaySection('Health-filtered aggregate', replayData.healthAggregateTable);
  printMfReplaySection('CRB-normalized trimmed aggregate', replayData.trimAggregateTable);
  if height(replayData.tailTable) <= 12
    printMfReplaySection('Rejected / trimmed tail cases', replayData.tailTable);
  else
    printMfReplaySection('Rejected / trimmed tail cases', table(height(replayData.tailTable), ...
      'VariableNames', {'numTailRows'}));
    fprintf('Rejected / trimmed tail cases preview:\n');
    dispMfReplayTablePreview(replayData.tailTable, 6);
  end
  localPrintProbeSection('MF objective probe', replayData.objectiveProbeTable, 18);
  localPrintProbeSection('MF solve probe', replayData.solveProbeTable, 18);
  localPrintProbeSection('MF DoA path probe', replayData.pathProbeTable, 18);
  localPrintProbeSection('MF diagnosis decision table', replayData.diagnosisTable, 18);
  printMfReplaySection('CRB metric projection check', replayData.crbMetricTable);
  printMfReplaySection('CP-K / CP-U release compare', replayData.releaseCompareTable);
  printMfReplaySection('Oracle in-tooth range summary', replayData.rangeTable);
  fprintf('Replay purpose                    : diagnose metric/init/solver causes, not paper-facing statistics.\n');
  fprintf('Angle error metric                : spherical great-circle distance in degrees.\n');
  fprintf('CRB spherical metric              : projectCrbToAngleMetric(latlon) trace(J*C*J^T).\n');
fprintf('fdRef CRB ratio                   : fdRefRmseOverCrb keeps legacy median denominator; fdRefNormRmse is seed-normalized.\n');
  fprintf('Raw aggregate filtering           : not applied; all rows are diagnostic rows.\n');
  fprintf('Health aggregate filtering        : solver/frequency/boundary/no-solve health only.\n');
  fprintf('Trim aggregate filtering          : health + CRB-normalized angle/fdRef cap.\n');
  fprintf('Objective probe                   : compare static seed, final, truth and mixed points in the MF objective.\n');
  fprintf('Solve probe                       : rerun truth-DoA and wide-static-start probes to locate basin / solver misses.\n');
  replayData.plotData = localPlotReplay(replayData);

  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "DONE", ...
    'config', config, ...
    'snapshotFile', replayData.snapshotFile, ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', replayData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(replayData), ...
    'commentLineList', [ ...
      "SS metric diagnose completed; inspect objectiveProbeTable, solveProbeTable, pathProbeTable, and diagnosisTable before changing estimator logic."; ...
      "No estimator, scan, or regression path is modified by this replay."]));

catch ME
  if exist('progressTracker', 'var')
    localCloseProgressTracker(progressTracker);
  end
  if exist('progressCleanup', 'var')
    clear progressCleanup;
  end
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME, ...
    'commentLineList', "SS MLE/CRB metric diagnose failed."));
  rethrow(ME);
end

%% Local helpers

function checkpointOpt = localBuildCheckpointOpt(replayName, config, useParfor)
%LOCALBUILDCHECKPOINTOPT Build stable per-task checkpoint options for this replay.

checkpointMeta = struct( ...
  'snrDbList', reshape(config.snrDbList, 1, []), ...
  'seedList', reshape(config.seedList, 1, []), ...
  'methodNameList', reshape(config.methodNameList, 1, []), ...
  'numFrame', config.numFrame, ...
  'oracleFdHalfToothFraction', config.oracleFdHalfToothFraction, ...
  'oracleFdRateHalfWidthHzPerSec', config.oracleFdRateHalfWidthHzPerSec, ...
  'initMode', config.initMode, ...
  'healthFdRefToothFraction', config.healthFdRefToothFraction, ...
  'healthFdRateAbsMaxHzPerSec', config.healthFdRateAbsMaxHzPerSec, ...
  'trimAngleNormMax', config.trimAngleNormMax, ...
  'trimFdRefNormMax', config.trimFdRefNormMax, ...
  'objectiveProbeEnable', config.objectiveProbeEnable, ...
  'solveProbeEnable', config.solveProbeEnable, ...
  'pathProbeEnable', config.pathProbeEnable, ...
  'probeWideDoaHalfWidthDegList', reshape(config.probeWideDoaHalfWidthDegList, 1, []), ...
  'probeTruthDoaHalfWidthDeg', reshape(config.probeTruthDoaHalfWidthDeg, 1, []), ...
  'objectivePathAlphaList', reshape(config.objectivePathAlphaList, 1, []));
checkpointOpt = buildMfReplayCheckpointOpt(replayName, config, struct( ...
  'runKey', localBuildCheckpointRunKey(config), ...
  'meta', checkpointMeta, ...
  'useParfor', logical(useParfor)));
end

function runKey = localBuildCheckpointRunKey(config)
%LOCALBUILDCHECKPOINTRUNKEY Build a compact stable key for interrupted replays.

seedList = reshape(double(config.seedList), 1, []);
snrList = reshape(double(config.snrDbList), 1, []);
signature = localBuildCheckpointSignature(config);
signatureHash = localCompactStringHash(signature);
runKey = sprintf('frames%d_snr%.6gto%.6g_seed%dto%d_rep%d_%s', ...
  config.numFrame, snrList(1), snrList(end), seedList(1), seedList(end), ...
  numel(seedList), signatureHash);
runKey = localPathSafeToken(runKey);
end

function signature = localBuildCheckpointSignature(config)
%LOCALBUILDCHECKPOINTSIGNATURE Build the full semantic signature behind runKey.

signature = strjoin(string({ ...
  sprintf('frames%d', config.numFrame), ...
  sprintf('snr%s', mat2str(reshape(double(config.snrDbList), 1, []), 8)), ...
  sprintf('seed%s', mat2str(reshape(double(config.seedList), 1, []), 8)), ...
  sprintf('methods%s', char(strjoin(reshape(string(config.methodNameList), 1, []), ','))), ...
  sprintf('oracleFdFrac%.8g', config.oracleFdHalfToothFraction), ...
  sprintf('oracleFdRate%.8g', config.oracleFdRateHalfWidthHzPerSec), ...
  sprintf('doaHalf%s', mat2str(reshape(double(config.staticLocalDoaHalfWidthDeg), 1, []), 8)), ...
  sprintf('init%s', char(string(config.initMode))), ...
  sprintf('healthFdRef%.8g', config.healthFdRefToothFraction), ...
  sprintf('healthFdRate%.8g', config.healthFdRateAbsMaxHzPerSec), ...
  sprintf('trimAngle%.8g', config.trimAngleNormMax), ...
  sprintf('trimFdRef%.8g', config.trimFdRefNormMax), ...
  sprintf('objectiveProbe%d', logical(config.objectiveProbeEnable)), ...
  sprintf('solveProbe%d', logical(config.solveProbeEnable)), ...
  sprintf('pathProbe%d', logical(config.pathProbeEnable)), ...
  sprintf('wideDoa%s', mat2str(reshape(double(config.probeWideDoaHalfWidthDegList), 1, []), 8)), ...
  sprintf('truthDoaHalf%s', mat2str(reshape(double(config.probeTruthDoaHalfWidthDeg), 1, []), 8)), ...
  sprintf('pathAlpha%s', mat2str(reshape(double(config.objectivePathAlphaList), 1, []), 8)) ...
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
token = regexprep(char(token), '[^A-Za-z0-9_]', '');
token = string(token);
end

function [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, replayName)
%LOCALRUNTASKBATCH Run SNR/seed tasks with optional checkpoint/resume.

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
    checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker);
    runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, @localCheckpointTaskRunner, checkpointRunnerOpt);
    repeatCell = runState.resultCell;
  else
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
  end
  localCloseProgressTracker(progressTracker);
catch ME
  localCloseProgressTracker(progressTracker);
  rethrow(ME);
end
end

function repeatOut = localCheckpointTaskRunner(taskInfo, sharedData)
%LOCALCHECKPOINTTASKRUNNER Run one checkpointed SNR/seed task.

repeatOut = localRunOneTask(sharedData.context, sharedData.flowOpt, sharedData.config, taskInfo);
end

function checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only fields accepted by checkpoint runner.

checkpointRunnerOpt = struct();
checkpointRunnerOpt.runName = checkpointOpt.runName;
checkpointRunnerOpt.runKey = checkpointOpt.runKey;
checkpointRunnerOpt.outputRoot = checkpointOpt.outputRoot;
checkpointRunnerOpt.useParfor = logical(localGetFieldOrDefault(checkpointOpt, 'useParfor', false));
checkpointRunnerOpt.resume = checkpointOpt.resume;
checkpointRunnerOpt.meta = checkpointOpt.meta;
checkpointRunnerOpt.progressFcn = @(step) localAdvanceProgressByStep(progressTracker, step);
checkpointRunnerOpt.cleanupOnSuccess = false;
checkpointRunnerOpt.cleanupOpt = struct();
end

function numDone = localCountCheckpointTaskFile(taskDir, numTask)
%LOCALCOUNTCHECKPOINTTASKFILE Count existing completed checkpoint task files.

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

for iStep = 1:step
  localAdvanceProgress(progressTracker);
end
end

function textValue = localCompactNumberList(valueList)
%LOCALCOMPACTNUMBERLIST Format numeric axes for a checkpoint run key.

valueList = reshape(double(valueList), 1, []);
textValue = strjoin(compose('%.6g', valueList), '_');
textValue = replace(textValue, '-', 'm');
textValue = replace(textValue, '.', 'p');
end

function localPrintReplayConfig(config)
%LOCALPRINTREPLAYCONFIG Print replay-specific axes and oracle settings.

fprintf('  %-32s : %s\n', 'seed list', localFormatRow(config.seedList));
fprintf('  %-32s : %s\n', 'SNR list (dB)', localFormatRow(config.snrDbList));
fprintf('  %-32s : %d\n', 'frame count', config.numFrame);
fprintf('  %-32s : %s\n', 'method list', strjoin(cellstr(config.methodNameList), ', '));
fprintf('  %-32s : %.3f tooth\n', 'fdRef oracle half-width', config.oracleFdHalfToothFraction);
fprintf('  %-32s : %.2f Hz/s\n', 'fdRate oracle half-width', config.oracleFdRateHalfWidthHzPerSec);
fprintf('  %-32s : [%s] deg\n', 'DoA release half-width', localFormatRow(config.staticLocalDoaHalfWidthDeg));
fprintf('  %-32s : %s\n', 'MF init mode', char(config.initMode));
fprintf('  %-32s : %.3f tooth\n', 'health fdRef threshold', config.healthFdRefToothFraction);
fprintf('  %-32s : %.2f Hz/s\n', 'health fdRate threshold', config.healthFdRateAbsMaxHzPerSec);
fprintf('  %-32s : %.2f / %.2f\n', 'trim norm angle/fdRef', config.trimAngleNormMax, config.trimFdRefNormMax);
fprintf('  %-32s : %d / %d / %d\n', 'objective/solve/path probes', ...
  config.objectiveProbeEnable, config.solveProbeEnable, config.pathProbeEnable);
fprintf('  %-32s : %d / %d\n', 'task count / outer parfor', ...
  config.numTask, logical(localGetFieldOrDefault(config, 'useParfor', false)));
fprintf('  %-32s : %s deg\n', 'wide DoA probe half-width', localFormatRow(config.probeWideDoaHalfWidthDegList));
end


function localPrintProbeSection(sectionTitle, probeTable, maxPreviewRow)
%LOCALPRINTPROBESECTION Print probe tables compactly.

if nargin < 3 || isempty(maxPreviewRow)
  maxPreviewRow = 18;
end
if isempty(probeTable) || height(probeTable) <= maxPreviewRow
  printMfReplaySection(sectionTitle, probeTable);
  return;
end
printMfReplaySection(sectionTitle, table(height(probeTable), 'VariableNames', {'numRows'}));
fprintf('%s preview:\n', sectionTitle);
dispMfReplayTablePreview(probeTable, maxPreviewRow);
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

function repeatOut = localRunOneTask(context, flowOpt, config, task)
%LOCALRUNONETASK Run one SNR/seed diagnostic task.

lastwarn('', '');
tRepeat = tic;
repeatData = buildDynamicRepeatData(context, task.snrDb, task.taskSeed);
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('replayMfSsMleCrbMetricDiagnose:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for this in-tooth diagnostic replay.');
end

fdHalfWidthHz = config.oracleFdHalfToothFraction * toothStepHz;
fdRangeOracle = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeOracle = truthFdRateHzPerSec + config.oracleFdRateHalfWidthHzPerSec * [-1, 1];

staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fixture.fdRange, truth, context.otherSatIdxGlobal, config.optVerbose, ...
  flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);

noiseVar = 1 / (10^(task.snrDb / 10));
crbBundle = localBuildCrbBundleQuiet(fixture, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, noiseVar);

caseRowList = repmat(localEmptyCaseRow(), numel(config.methodNameList), 1);
objectiveProbeRowList = repmat(localEmptyObjectiveProbeRow(), 0, 1);
solveProbeRowList = repmat(localEmptySolveProbeRow(), 0, 1);
pathProbeRowList = repmat(localEmptyPathProbeRow(), 0, 1);
for iMethod = 1:numel(config.methodNameList)
  methodName = config.methodNameList(iMethod);
  if methodName == "SS-SF-Static"
    caseUse = staticBundle.caseStaticRefOnly;
    caseRowList(iMethod) = localBuildCaseRow(methodName, caseUse, truth, task, ...
      fixture, toothStepHz, fdRangeOracle, fdRateRangeOracle, crbBundle, config, []);
  else
    method = localBuildDynamicMethod(methodName, config.staticLocalDoaHalfWidthDeg);
    [caseUse, initMeta, objectiveRows, solveRows, pathRows] = localRunDynamicMethod(method, fixture, truth, context, ...
      flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz, config.optVerbose, ...
      staticBundle.caseStaticRefOnly, config, task);
    objectiveProbeRowList = [objectiveProbeRowList; objectiveRows(:)]; %#ok<AGROW>
    solveProbeRowList = [solveProbeRowList; solveRows(:)]; %#ok<AGROW>
    pathProbeRowList = [pathProbeRowList; pathRows(:)]; %#ok<AGROW>
    caseRowList(iMethod) = localBuildCaseRow(methodName, caseUse, truth, task, ...
      fixture, toothStepHz, fdRangeOracle, fdRateRangeOracle, crbBundle, config, initMeta);
  end
end
[warningMessage, warningId] = lastwarn;

repeatOut = struct();
repeatOut.taskSeed = task.taskSeed;
repeatOut.snrDb = task.snrDb;
repeatOut.truth = truth;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.toothStepHz = toothStepHz;
repeatOut.fdRangeOracle = fdRangeOracle;
repeatOut.fdRateRangeOracle = fdRateRangeOracle;
repeatOut.caseRowList = caseRowList;
repeatOut.objectiveProbeRowList = objectiveProbeRowList;
repeatOut.solveProbeRowList = solveProbeRowList;
repeatOut.pathProbeRowList = pathProbeRowList;
repeatOut.repeatTotalMs = 1000 * toc(tRepeat);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function method = localBuildDynamicMethod(methodName, doaHalfWidthDeg)
%LOCALBUILDDYNAMICMETHOD Build a single-satellite MF method descriptor.

methodName = string(methodName);
switch methodName
  case "SS-MF-CP-K"
    fdRateMode = "known";
    isKnownRate = true;
  case "SS-MF-CP-U"
    fdRateMode = "unknown";
    isKnownRate = false;
  otherwise
    error('replayMfSsMleCrbMetricDiagnose:UnknownMethod', ...
      'Unsupported method name "%s".', char(methodName));
end
method = struct();
method.displayName = methodName;
method.satMode = "single";
method.phaseMode = "continuous";
method.fdRateMode = fdRateMode;
method.isKnownRate = logical(isKnownRate);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
end

function [caseUse, initMeta, objectiveProbeRows, solveProbeRows, pathProbeRows] = localRunDynamicMethod(method, fixture, truth, context, ...
  flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz, optVerbose, staticSeedCase, config, task)
%LOCALRUNDYNAMICMETHOD Run one single-satellite in-tooth dynamic method.
% The default diagnostic mode uses the static DoA seed and truth-centered
% frequency seed to inspect local in-tooth MLE behavior. Internal estimator
% initialization remains available as an explicit stress mode.

viewUse = fixture.viewRefOnly;
debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());

fdRateRangeUse = fdRateRangeOracle;
if method.isKnownRate
  fdRateRangeUse = [];
end

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.fdRateMode = char(method.fdRateMode);
if method.isKnownRate
  dynOpt.fdRateKnown = truth.fdRateFit;
end
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.enableFdAliasUnwrap = true;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);

[initParamOverride, initRequestMeta] = localBuildDynamicInitOverride( ...
  method, staticSeedCase, truth, config);
if isfield(initRequestMeta, 'initDoaParam') && ~isempty(initRequestMeta.initDoaParam)
  dynOpt.initDoaParam = initRequestMeta.initDoaParam(:);
  dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
end

caseUse = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
  viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fdRangeOracle, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
  debugTruthUse, initParamOverride);
caseUse.inToothMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeOracle, ...
  'fdRateRangeUse', fdRateRangeUse);

initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta);
[objectiveProbeRows, pathProbeRows] = localBuildObjectiveAndPathProbeRows(method, caseUse, staticSeedCase, truth, task, ...
  viewUse, context, fdRangeOracle, fdRateRangeUse, dynOpt, initParamOverride, config);
solveProbeRows = localBuildSolveProbeRows(method, caseUse, staticSeedCase, truth, task, ...
  viewUse, context, fdRangeOracle, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, config);
end

function [initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, staticSeedCase, truth, config)
%LOCALBUILDDYNAMICINITOVERRIDE Build optional external local-MLE seed.

initMode = string(localGetFieldOrDefault(config, 'mfInitMode', "static-doa-truth-fd"));
initParamOverride = [];
initRequestMeta = struct('initMode', initMode, 'initParam', [], ...
  'initDoaParam', [], 'initDoaHalfWidth', []);

switch initMode
  case "internal-estimator"
    return;

  case "static-doa-truth-fd"
    estResult = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
    seedDoaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
    if numel(seedDoaParam) < 2 || any(~isfinite(seedDoaParam(1:2)))
      error('replayMfSsMleCrbMetricDiagnose:MissingStaticDoaSeed', ...
        'Static DoA seed is required for mfInitMode="static-doa-truth-fd".');
    end
    seedDoaParam = seedDoaParam(1:2);
    truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
    truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
    if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
      error('replayMfSsMleCrbMetricDiagnose:MissingTruthFrequencySeed', ...
        'Truth fdRef/fdRate seeds are required for mfInitMode="static-doa-truth-fd".');
    end
    if method.isKnownRate
      initParamOverride = [seedDoaParam(:); truthFdRefHz];
    else
      initParamOverride = [seedDoaParam(:); truthFdRefHz; truthFdRateHzPerSec];
    end
    initRequestMeta.initParam = initParamOverride(:);
    initRequestMeta.initDoaParam = seedDoaParam(:);
    initRequestMeta.initDoaHalfWidth = reshape(config.staticLocalDoaHalfWidthDeg, [], 1);

  otherwise
    error('replayMfSsMleCrbMetricDiagnose:UnknownInitMode', ...
      'Unsupported mfInitMode "%s".', char(initMode));
end
end

function initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta)
%LOCALBUILDDYNAMICINITMETA Extract requested and estimator-internal initializer for summary.

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
initMeta = struct();
initMeta.requestedInitMode = string(localGetFieldOrDefault(initRequestMeta, 'initMode', ""));
initMeta.requestedInitParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initParam', []), [], 1);
initMeta.requestedInitDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaParam', []), [], 1);
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
initMeta.searchHalfWidthDeg = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaHalfWidth', NaN(2, 1)), [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
end


function [objectiveProbeRows, pathProbeRows] = localBuildObjectiveAndPathProbeRows(method, caseUse, staticSeedCase, truth, task, ...
  viewUse, context, fdRangeOracle, fdRateRangeUse, dynOpt, initParamOverride, config)
%LOCALBUILDOBJECTIVEANDPATHPROBEROWS Evaluate fixed points in the exact MF objective.

objectiveProbeRows = repmat(localEmptyObjectiveProbeRow(), 0, 1);
pathProbeRows = repmat(localEmptyPathProbeRow(), 0, 1);
if ~logical(localGetFieldOrDefault(config, 'objectiveProbeEnable', true)) && ...
    ~logical(localGetFieldOrDefault(config, 'pathProbeEnable', true))
  return;
end

try
  [model, ~, ~] = buildDoaDopplerMfModel(viewUse.sceneSeq, viewUse.rxSigMf, ...
    context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
    viewUse.doaGrid, fdRangeOracle, fdRateRangeUse, dynOpt);
  [initParamUse, ~, model] = buildDoaDopplerMfInit(model, initParamOverride);
catch ME
  row = localEmptyObjectiveProbeRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = "model-build-failed";
  row.evalOk = false;
  row.evalMessage = string(ME.message);
  objectiveProbeRows = row;
  return;
end

pointMap = localBuildObjectiveProbePointMap(method, caseUse, staticSeedCase, truth, initParamUse);
if logical(localGetFieldOrDefault(config, 'objectiveProbeEnable', true))
  tagList = string(fieldnames(pointMap));
  tagList = tagList(:);
  rowList = repmat(localEmptyObjectiveProbeRow(), numel(tagList), 1);
  finalObj = localEvaluateProbeObjective(model, localGetFieldOrDefault(pointMap, 'finalPoint', []));
  for iTag = 1:numel(tagList)
    optVar = pointMap.(char(tagList(iTag)));
    rowList(iTag) = localBuildObjectiveProbeRow(model, method, truth, task, tagList(iTag), optVar, finalObj);
  end
  objectiveProbeRows = rowList(:);
end

if logical(localGetFieldOrDefault(config, 'pathProbeEnable', true))
  pathProbeRows = localBuildPathProbeRows(model, method, truth, task, pointMap, config);
end
end

function pointMap = localBuildObjectiveProbePointMap(method, caseUse, staticSeedCase, truth, initParamUse)
%LOCALBUILDOBJECTIVEPROBEPOINTMAP Build static/final/truth/mixed objective points.

truthDoa = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
truthDoa = truthDoa(1:min(2, numel(truthDoa)));
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});

staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
staticDoa = reshape(localGetFieldOrDefault(staticEst, 'doaParamEst', []), [], 1);
staticDoa = staticDoa(1:min(2, numel(staticDoa)));

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
finalDoa = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
finalDoa = finalDoa(1:min(2, numel(finalDoa)));
finalFdRefHz = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
finalFdRateHzPerSec = localGetFieldOrDefault(estResult, 'fdRateEst', truthFdRateHzPerSec);

if method.isKnownRate
  pointMap = struct( ...
    'initPoint', localMakeMfOptVar(method, localTakeTwo(initParamUse), initParamUse(3), truthFdRateHzPerSec), ...
    'staticSeedPoint', localMakeMfOptVar(method, staticDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'finalPoint', localMakeMfOptVar(method, finalDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'truthPoint', localMakeMfOptVar(method, truthDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'truthDoaFinalFdPoint', localMakeMfOptVar(method, truthDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'finalDoaTruthFdPoint', localMakeMfOptVar(method, finalDoa, truthFdRefHz, truthFdRateHzPerSec));
else
  pointMap = struct( ...
    'initPoint', localMakeMfOptVar(method, localTakeTwo(initParamUse), initParamUse(3), initParamUse(4)), ...
    'staticSeedPoint', localMakeMfOptVar(method, staticDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'finalPoint', localMakeMfOptVar(method, finalDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'truthPoint', localMakeMfOptVar(method, truthDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'truthDoaFinalFdPoint', localMakeMfOptVar(method, truthDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'finalDoaTruthFdPoint', localMakeMfOptVar(method, finalDoa, truthFdRefHz, truthFdRateHzPerSec));
end
end

function optVar = localMakeMfOptVar(method, doaParam, fdRefHz, fdRateHzPerSec)
%LOCALMAKEMFOPTVAR Build one MF optimization vector.

doaParam = localTakeTwo(doaParam);
if method.isKnownRate
  optVar = [doaParam(:); fdRefHz];
else
  optVar = [doaParam(:); fdRefHz; fdRateHzPerSec];
end
end

function twoValue = localTakeTwo(value)
%LOCALTAKETWO Return a finite two-vector or NaNs.

twoValue = NaN(2, 1);
value = reshape(double(value), [], 1);
if numel(value) >= 2
  twoValue = value(1:2);
end
end

function row = localBuildObjectiveProbeRow(model, method, truth, task, probeTag, optVar, finalObj)
%LOCALBUILDOBJECTIVEPROBEROW Build one point-evaluation row.

row = localEmptyObjectiveProbeRow();
row.displayName = method.displayName;
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.probeTag = string(probeTag);
row.objective = localEvaluateProbeObjective(model, optVar);
row.finalObjective = finalObj;
row.objMinusFinal = row.objective - finalObj;
row.evalOk = isfinite(row.objective);
if ~row.evalOk
  row.evalMessage = "objective-eval-failed";
end
row.angleErrDeg = localProbeAngleErr(optVar, truth);
[row.fdRefErrHz, row.fdRateErrHzPerSec] = localProbeFdErr(method, optVar, truth);
end

function obj = localEvaluateProbeObjective(model, optVar)
%LOCALEVALUATEPROBEOBJECTIVE Safely evaluate one MF objective point.

obj = NaN;
try
  optVar = reshape(double(optVar), [], 1);
  if numel(optVar) ~= localGetProbeNumVar(model) || any(~isfinite(optVar))
    return;
  end
  obj = evaluateDoaDopplerMfObjective(model, optVar);
catch
  obj = NaN;
end
end

function numVar = localGetProbeNumVar(model)
%LOCALGETPROBENUMVAR Number of optimizer variables for an MF model.

numVar = 3;
if isfield(model, 'fdRateMode') && strcmp(model.fdRateMode, 'unknown')
  numVar = 4;
end
end

function angleErrDeg = localProbeAngleErr(optVar, truth)
%LOCALPROBEANGLEERR Compute angular error for one probe point.

angleErrDeg = NaN;
truthDoa = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
optVar = reshape(double(optVar), [], 1);
if numel(optVar) >= 2 && numel(truthDoa) >= 2 && all(isfinite(optVar(1:2)))
  angleErrDeg = calcLatlonAngleError(optVar(1:2), truthDoa(1:2));
end
end

function [fdRefErrHz, fdRateErrHzPerSec] = localProbeFdErr(method, optVar, truth)
%LOCALPROBEFDERR Compute fdRef/fdRate errors for one probe point.

fdRefErrHz = NaN;
fdRateErrHzPerSec = NaN;
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
optVar = reshape(double(optVar), [], 1);
if numel(optVar) >= 3
  fdRefErrHz = optVar(3) - truthFdRefHz;
end
if method.isKnownRate
  fdRateErrHzPerSec = 0;
elseif numel(optVar) >= 4
  fdRateErrHzPerSec = optVar(4) - truthFdRateHzPerSec;
end
end

function pathRows = localBuildPathProbeRows(model, method, truth, task, pointMap, config)
%LOCALBUILDPATHPROBEROWS Evaluate compact DoA interpolation paths.

alphaList = reshape(double(config.objectivePathAlphaList), [], 1);
pathRows = repmat(localEmptyPathProbeRow(), 0, 1);
pathRows = [pathRows; localBuildOnePathProbeRows(model, method, truth, task, ...
  "final-to-truth-doa-final-fd", pointMap.finalPoint, pointMap.truthDoaFinalFdPoint, alphaList)];
pathRows = [pathRows; localBuildOnePathProbeRows(model, method, truth, task, ...
  "static-to-truth-doa-truth-fd", pointMap.staticSeedPoint, pointMap.truthPoint, alphaList)];
end

function pathRows = localBuildOnePathProbeRows(model, method, truth, task, pathTag, startPoint, endPoint, alphaList)
%LOCALBUILDONEPATHPROBEROWS Evaluate one path with linear interpolation in opt-var space.

numAlpha = numel(alphaList);
pathRows = repmat(localEmptyPathProbeRow(), numAlpha, 1);
finalObj = localEvaluateProbeObjective(model, startPoint);
for iAlpha = 1:numAlpha
  alpha = alphaList(iAlpha);
  optVar = (1 - alpha) .* startPoint(:) + alpha .* endPoint(:);
  row = localEmptyPathProbeRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.pathTag = string(pathTag);
  row.alpha = alpha;
  row.objective = localEvaluateProbeObjective(model, optVar);
  row.objMinusStart = row.objective - finalObj;
  row.angleErrDeg = localProbeAngleErr(optVar, truth);
  [row.fdRefErrHz, row.fdRateErrHzPerSec] = localProbeFdErr(method, optVar, truth);
  row.evalOk = isfinite(row.objective);
  pathRows(iAlpha) = row;
end
end

function solveRows = localBuildSolveProbeRows(method, baselineCase, staticSeedCase, truth, task, ...
  viewUse, context, fdRangeOracle, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, config)
%LOCALBUILDSOLVEPROBEROWS Rerun controlled init/bounds probes for MF diagnosis.

solveRows = repmat(localEmptySolveProbeRow(), 0, 1);
solveRows = [solveRows; localBuildSolveProbeRowFromCase(method, baselineCase, truth, task, ...
  "baseline", config.staticLocalDoaHalfWidthDeg, NaN, "baseline")];
if ~logical(localGetFieldOrDefault(config, 'solveProbeEnable', true))
  return;
end

truthDoa = localTakeTwo(localGetFieldOrDefault(truth, 'latlonTrueDeg', []));
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
staticDoa = localTakeTwo(localGetFieldOrDefault(staticEst, 'doaParamEst', []));
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});

probeList = repmat(struct('tag', "", 'centerDoa', NaN(2, 1), 'halfWidth', NaN(2, 1)), 0, 1);
probeItem = struct('tag', "truth-doa-truth-fd", 'centerDoa', truthDoa, ...
  'halfWidth', localExpandDoaHalfWidth(config.probeTruthDoaHalfWidthDeg));
probeList = [probeList; probeItem]; %#ok<AGROW>
wideList = reshape(double(config.probeWideDoaHalfWidthDegList), [], 1);
for iWide = 1:numel(wideList)
  probeItem = struct('tag', string(sprintf('static-doa-wide%.4g', wideList(iWide))), ...
    'centerDoa', staticDoa, 'halfWidth', [wideList(iWide); wideList(iWide)]);
  probeList = [probeList; probeItem]; %#ok<AGROW>
end

for iProbe = 1:numel(probeList)
  initParam = localMakeMfOptVar(method, probeList(iProbe).centerDoa, truthFdRefHz, truthFdRateHzPerSec);
  initCandidate = struct('initParam', initParam, 'initDoaParam', probeList(iProbe).centerDoa, ...
    'initDoaHalfWidth', probeList(iProbe).halfWidth, 'startTag', probeList(iProbe).tag);
  try
    probeCase = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
      viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
      fdRangeOracle, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
      debugTruthUse, initCandidate);
    solveRows = [solveRows; localBuildSolveProbeRowFromCase(method, probeCase, truth, task, ...
      probeList(iProbe).tag, probeList(iProbe).halfWidth, solveRows(1).finalObj, "probe")]; %#ok<AGROW>
  catch ME
    row = localEmptySolveProbeRow();
    row.displayName = method.displayName;
    row.snrDb = task.snrDb;
    row.taskSeed = task.taskSeed;
    row.probeTag = probeList(iProbe).tag;
    row.probeGroup = "probe";
    row.doaHalfWidthLatDeg = probeList(iProbe).halfWidth(1);
    row.doaHalfWidthLonDeg = probeList(iProbe).halfWidth(2);
    row.evalOk = false;
    row.message = string(ME.message);
    solveRows = [solveRows; row]; %#ok<AGROW>
  end
end
end

function halfWidth = localExpandDoaHalfWidth(value)
%LOCALEXPANDDOAHALFWIDTH Convert scalar/two-vector half-width to 2x1.

value = reshape(double(value), [], 1);
if isempty(value)
  halfWidth = [NaN; NaN];
elseif numel(value) == 1
  halfWidth = repmat(value, 2, 1);
else
  halfWidth = value(1:2);
end
end

function row = localBuildSolveProbeRowFromCase(method, caseUse, truth, task, probeTag, halfWidthDeg, baselineObj, probeGroup)
%LOCALBUILDSOLVEPROBEROWFROMCASE Convert one probe solve to a compact row.

row = localEmptySolveProbeRow();
row.displayName = method.displayName;
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.probeTag = string(probeTag);
row.probeGroup = string(probeGroup);
row.doaHalfWidthLatDeg = halfWidthDeg(1);
row.doaHalfWidthLonDeg = halfWidthDeg(min(2, numel(halfWidthDeg)));
info = summarizeDoaDopplerCase(caseUse, truth);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
row.finalAngleErrDeg = info.angleErrDeg;
row.finalFdRefErrHz = info.fdRefErrHz;
row.finalFdRateErrHzPerSec = info.fdRateErrHzPerSec;
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objMinusBaseline = row.finalObj - baselineObj;
row.iterations = info.iterations;
row.exitflag = info.exitflag;
row.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));
row.candidateCount = numel(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.boundaryHit = false;
row.evalOk = isfinite(row.finalObj);
row.message = "";
end

function row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
  fdRangeOracle, fdRateRangeOracle, crbBundle, config, initMeta)
%LOCALBUILDCASEMETRICROW Convert one static or dynamic case into a diagnostic row.

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
row.fdRefTrueHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
row.fdRateTrueHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
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

if methodName == "SS-SF-Static"
  row.satMode = "single";
  row.frameMode = "single";
  row.phaseMode = "single-frame";
  row.fdRateMode = "zero";
  row.initLatDeg = row.finalLatDeg;
  row.initLonDeg = row.finalLonDeg;
  row.initAngleErrDeg = row.finalAngleErrDeg;
  row.finalMinusInitAngleDeg = 0;
  row.initFdRefHz = row.finalFdRefHz;
  row.initFdRateHzPerSec = NaN;
  row.searchHalfWidthLatDeg = NaN;
  row.searchHalfWidthLonDeg = NaN;
  row.fdRateRangeLowHzPerSec = NaN;
  row.fdRateRangeHighHzPerSec = NaN;
  row.mfInitMode = "static-baseline";
else
  row.satMode = "single";
  row.frameMode = "multi";
  row.phaseMode = "continuous";
  if methodName == "SS-MF-CP-K"
    row.fdRateMode = "known";
  else
    row.fdRateMode = "unknown";
  end
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
row.fdRefSignedErrOverCrb = localSafeRatio(row.fdRefErrHz, row.fdRefCrbStdHz);

row.fdRangeLowHz = fdRangeOracle(1);
row.fdRangeHighHz = fdRangeOracle(2);
row.fdRangeHalfWidthHz = 0.5 * abs(fdRangeOracle(2) - fdRangeOracle(1));
if methodName ~= "SS-SF-Static"
  row.fdRateRangeLowHzPerSec = fdRateRangeOracle(1);
  row.fdRateRangeHighHzPerSec = fdRateRangeOracle(2);
end
row.fdRefBoundaryHit = localIsNearRangeBoundary(row.finalFdRefHz, fdRangeOracle, config.boundaryTolFraction);
row.fdRateBoundaryHit = methodName == "SS-MF-CP-U" && ...
  localIsNearRangeBoundary(row.finalFdRateHzPerSec, fdRateRangeOracle, config.boundaryTolFraction);

row.initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objectiveImprove = row.initObj - row.finalObj;
row.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
row.refCoherence = NaN;
row.nonRefCoherenceFloor = NaN;
if isfield(finalEval, 'coherenceSat') && ~isempty(finalEval.coherenceSat)
  row.refCoherence = finalEval.coherenceSat(1);
end
row.stepSize = localGetFieldOrDefault(optimInfo, 'stepsize', NaN);
row.firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
row.constrViolation = localGetFieldOrDefault(optimInfo, 'constrviolation', NaN);
row.algorithm = string(localGetFieldOrDefault(optimInfo, 'algorithm', ""));
row.usedScaledSolve = logical(localGetFieldOrDefault(optimInfo, 'usedScaledSolve', false));
row.fallbackTriggered = logical(localGetFieldOrDefault(optimInfo, 'fallbackTriggered', false));
row.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));
row.candidateCount = numel(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.candidateVariantText = localJoinStringList(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.fdRefInitMoveHz = row.finalFdRefHz - row.initFdRefHz;
row.fdRateInitMoveHzPerSec = row.finalFdRateHzPerSec - row.initFdRateHzPerSec;
row.warningFlag = false;
end

function [crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName)
%LOCALRESOLVECRBFORMETHOD Return the CRB block matching one method.

methodName = string(methodName);
switch methodName
  case "SS-SF-Static"
    crbFull = crbBundle.crbSfRef;
  case "SS-MF-CP-K"
    crbFull = crbBundle.crbMfRefKnown;
  case "SS-MF-CP-U"
    crbFull = crbBundle.crbMfRefUnknown;
  otherwise
    error('replayMfSsMleCrbMetricDiagnose:UnknownCrbMethod', ...
      'Unsupported CRB method "%s".', char(methodName));
end
crbDoa = real(crbFull(1:2, 1:2));
end

function crbBundle = localBuildCrbBundleQuiet(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%LOCALBUILDCRBBUNDLEQUIET Build CRB while suppressing expected full-FIM warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
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

function aggregateTable = localBuildAggregateTable(caseTable)
%LOCALBUILDAGGREGATETABLE Aggregate diagnostic rows by method and SNR.

if isempty(caseTable) || height(caseTable) == 0
  aggregateTable = struct2table(repmat(localEmptyAggregateRow(), 0, 1));
  return;
end
groupKey = unique(caseTable(:, {'displayName', 'frameMode', 'fdRateMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(localEmptyAggregateRow(), height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = caseTable.displayName == groupKey.displayName(iRow) & caseTable.snrDb == groupKey.snrDb(iRow);
  row = localEmptyAggregateRow();
  row.displayName = groupKey.displayName(iRow);
  row.frameMode = groupKey.frameMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.numRepeat = nnz(mask);
  row.angleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(mask));
  row.angleMedianDeg = median(caseTable.sphericalAngleErrDeg(mask), 'omitnan');
  row.angleMaxDeg = max(caseTable.sphericalAngleErrDeg(mask), [], 'omitnan');
  row.angleCrbSphericalMedianDeg = median(caseTable.crbSphericalApproxStdDeg(mask), 'omitnan');
  row.angleCrbTraceMedianDeg = median(caseTable.crbTraceStdDeg(mask), 'omitnan');
  row.angleRmseOverSphericalCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbSphericalMedianDeg);
  row.angleRmseOverTraceCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbTraceMedianDeg);
  row.fdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(mask));
  row.fdRefCrbMedianHz = median(caseTable.fdRefCrbStdHz(mask), 'omitnan');
  row.fdRefRmseOverCrb = localSafeRatio(row.fdRefRmseHz, row.fdRefCrbMedianHz);
  row = localFillFdRefCrbAggregate(row, caseTable.fdRefErrHz(mask), caseTable.fdRefCrbStdHz(mask));
  row.fdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(mask));
  row.initAngleMedianDeg = median(caseTable.initAngleErrDeg(mask), 'omitnan');
  row.finalMinusInitMedianDeg = median(caseTable.finalMinusInitAngleDeg(mask), 'omitnan');
  row.fdRefInitMoveMedianHz = median(abs(caseTable.fdRefInitMoveHz(mask)), 'omitnan');
  row.fdRateInitMoveMedianHzPerSec = median(abs(caseTable.fdRateInitMoveHzPerSec(mask)), 'omitnan');
  row.objectiveImproveMedian = median(caseTable.objectiveImprove(mask), 'omitnan');
  row.boundaryHitRate = mean(double(caseTable.fdRefBoundaryHit(mask) | caseTable.fdRateBoundaryHit(mask)), 'omitnan');
  row.candidateCountMedian = median(caseTable.candidateCount(mask), 'omitnan');
  row.iterationsMedian = median(caseTable.iterations(mask), 'omitnan');
  row.firstOrderOptMedian = median(caseTable.firstOrderOpt(mask), 'omitnan');
  rowList(iRow) = row;
end
aggregateTable = struct2table(rowList);
end

function row = localFillFdRefCrbAggregate(row, fdRefErrHz, fdRefCrbStdHz)
%LOCALFILLFDREFCRBAGGREGATE Add seed-normalized fdRef CRB diagnostics.

fdRefErrHz = double(fdRefErrHz(:));
fdRefCrbStdHz = double(fdRefCrbStdHz(:));
validMask = isfinite(fdRefErrHz) & isfinite(fdRefCrbStdHz) & fdRefCrbStdHz > 0;
if ~any(validMask)
  return;
end
errUse = fdRefErrHz(validMask);
crbUse = fdRefCrbStdHz(validMask);
crbMeanSquare = mean(crbUse.^2);
row.fdRefCrbMeanHz = mean(crbUse);
row.fdRefCrbRmsHz = sqrt(crbMeanSquare);
row.fdRefMseOverMeanCrb = mean(errUse.^2) ./ crbMeanSquare;
row.fdRefRmseOverMeanCrbStd = sqrt(row.fdRefMseOverMeanCrb);
row.fdRefNormRmse = sqrt(mean((errUse ./ crbUse).^2));
row.fdRefBiasHz = mean(errUse);
row.fdRefBiasOverMeanCrbStd = localSafeRatio(row.fdRefBiasHz, row.fdRefCrbRmsHz);
end

function caseTable = localAnnotateHealthFlags(caseTable, config)
%LOCALANNOTATEHEALTHFLAGS Add scan-like health and trim flags without hiding raw rows.

if isempty(caseTable) || height(caseTable) == 0
  return;
end
numRow = height(caseTable);
isDynamic = startsWith(caseTable.displayName, "SS-MF");
isUnknown = caseTable.displayName == "SS-MF-CP-U";

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
trimKeep = healthResolved & angleNormOk & fdRefNormOk;

caseTable.isDynamicMethod = isDynamic;
caseTable.solverHealthy = solverHealthy;
caseTable.freqHealthy = freqHealthy;
caseTable.rateHealthy = rateHealthy;
caseTable.boundaryHealthy = boundaryHealthy;
caseTable.notNoSolve = notNoSolve;
caseTable.healthResolved = healthResolved;
caseTable.trimKeep = trimKeep;
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

function filteredTable = localBuildFilteredAggregateTable(caseTable, filterName, keepMask)
%LOCALBUILDFILTEREDAGGREGATETABLE Aggregate rows after a health / trim filter.

if isempty(caseTable) || height(caseTable) == 0
  filteredTable = struct2table(repmat(localEmptyFilteredAggregateRow(), 0, 1));
  return;
end
keepMask = logical(keepMask(:));
groupKey = unique(caseTable(:, {'displayName', 'frameMode', 'fdRateMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(localEmptyFilteredAggregateRow(), height(groupKey), 1);
for iRow = 1:height(groupKey)
  groupMask = caseTable.displayName == groupKey.displayName(iRow) & caseTable.snrDb == groupKey.snrDb(iRow);
  mask = groupMask & keepMask;
  row = localEmptyFilteredAggregateRow();
  row.filterName = string(filterName);
  row.displayName = groupKey.displayName(iRow);
  row.frameMode = groupKey.frameMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.rawCount = nnz(groupMask);
  row.keepCount = nnz(mask);
  row.keepRate = localSafeRatio(row.keepCount, row.rawCount);
  row.rejectRate = 1 - row.keepRate;
  row.numRepeat = row.keepCount;
  if row.keepCount > 0
    row.angleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(mask));
    row.angleMedianDeg = median(caseTable.sphericalAngleErrDeg(mask), 'omitnan');
    row.angleMaxDeg = max(caseTable.sphericalAngleErrDeg(mask), [], 'omitnan');
    row.angleCrbSphericalMedianDeg = median(caseTable.crbSphericalApproxStdDeg(mask), 'omitnan');
    row.angleCrbTraceMedianDeg = median(caseTable.crbTraceStdDeg(mask), 'omitnan');
    row.angleRmseOverSphericalCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbSphericalMedianDeg);
    row.angleRmseOverTraceCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbTraceMedianDeg);
    row.fdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(mask));
    row.fdRefCrbMedianHz = median(caseTable.fdRefCrbStdHz(mask), 'omitnan');
    row.fdRefRmseOverCrb = localSafeRatio(row.fdRefRmseHz, row.fdRefCrbMedianHz);
    row = localFillFdRefCrbAggregate(row, caseTable.fdRefErrHz(mask), caseTable.fdRefCrbStdHz(mask));
    row.fdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(mask));
    row.initAngleMedianDeg = median(caseTable.initAngleErrDeg(mask), 'omitnan');
    row.finalMinusInitMedianDeg = median(caseTable.finalMinusInitAngleDeg(mask), 'omitnan');
    row.fdRefInitMoveMedianHz = median(abs(caseTable.fdRefInitMoveHz(mask)), 'omitnan');
    row.fdRateInitMoveMedianHzPerSec = median(abs(caseTable.fdRateInitMoveHzPerSec(mask)), 'omitnan');
    row.objectiveImproveMedian = median(caseTable.objectiveImprove(mask), 'omitnan');
    row.boundaryHitRate = mean(double(caseTable.fdRefBoundaryHit(mask) | caseTable.fdRateBoundaryHit(mask)), 'omitnan');
    row.candidateCountMedian = median(caseTable.candidateCount(mask), 'omitnan');
    row.iterationsMedian = median(caseTable.iterations(mask), 'omitnan');
    row.firstOrderOptMedian = median(caseTable.firstOrderOpt(mask), 'omitnan');
  end
  rowList(iRow) = row;
end
filteredTable = struct2table(rowList);
end

function tailTable = localBuildTailTable(caseTable)
%LOCALBUILDTAILTABLE Keep rows rejected by health or CRB-normalized trim filters.

if isempty(caseTable) || height(caseTable) == 0
  tailTable = table();
  return;
end
tailMask = caseTable.isDynamicMethod & (~caseTable.healthResolved | ~caseTable.trimKeep);
columnList = {'displayName', 'snrDb', 'taskSeed', 'mfInitMode', 'sphericalAngleErrDeg', ...
  'angleErrOverSphericalCrb', 'fdRefAbsErrHz', 'fdRefErrOverCrb', 'fdRefSignedErrOverCrb', ...
  'fdRateAbsErrHzPerSec', 'iterations', 'objectiveImprove', 'fdRefBoundaryHit', ...
  'fdRateBoundaryHit', 'healthResolved', 'trimKeep', 'rejectReason'};
columnList = columnList(ismember(columnList, caseTable.Properties.VariableNames));
tailTable = caseTable(tailMask, columnList);
if height(tailTable) > 0
  tailTable = sortrows(tailTable, {'snrDb', 'displayName', 'taskSeed'});
end
end

function crbMetricTable = localBuildCrbMetricTable(caseTable)
%LOCALBUILDCRBMETRICTABLE Summarize whether angle error and CRB use comparable metrics.

if isempty(caseTable) || height(caseTable) == 0
  crbMetricTable = table();
  return;
end
groupKey = unique(caseTable(:, {'displayName', 'snrDb'}), 'rows', 'stable');
crbMetricTable = groupKey;
crbMetricTable.latStdMedianDeg = NaN(height(groupKey), 1);
crbMetricTable.lonStdMedianDeg = NaN(height(groupKey), 1);
crbMetricTable.traceStdMedianDeg = NaN(height(groupKey), 1);
crbMetricTable.sphericalStdMedianDeg = NaN(height(groupKey), 1);
crbMetricTable.traceOverSphericalMedian = NaN(height(groupKey), 1);
crbMetricTable.errOverSphericalMedian = NaN(height(groupKey), 1);
crbMetricTable.errOverTraceMedian = NaN(height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = caseTable.displayName == groupKey.displayName(iRow) & caseTable.snrDb == groupKey.snrDb(iRow);
  crbMetricTable.latStdMedianDeg(iRow) = median(caseTable.crbLatStdDeg(mask), 'omitnan');
  crbMetricTable.lonStdMedianDeg(iRow) = median(caseTable.crbLonStdDeg(mask), 'omitnan');
  crbMetricTable.traceStdMedianDeg(iRow) = median(caseTable.crbTraceStdDeg(mask), 'omitnan');
  crbMetricTable.sphericalStdMedianDeg(iRow) = median(caseTable.crbSphericalApproxStdDeg(mask), 'omitnan');
  crbMetricTable.traceOverSphericalMedian(iRow) = median(caseTable.crbMetricTraceRatio(mask), 'omitnan');
  crbMetricTable.errOverSphericalMedian(iRow) = median(caseTable.angleErrOverSphericalCrb(mask), 'omitnan');
  crbMetricTable.errOverTraceMedian(iRow) = median(caseTable.angleErrOverTraceCrb(mask), 'omitnan');
end
end

function releaseCompareTable = localBuildReleaseCompareTable(caseTable)
%LOCALBUILDRELEASECOMPARETABLE Compare CP-U against CP-K for the same seed and SNR.

releaseCompareTable = table();
if isempty(caseTable) || height(caseTable) == 0
  return;
end
kMaskAll = caseTable.displayName == "SS-MF-CP-K";
uMaskAll = caseTable.displayName == "SS-MF-CP-U";
keyTable = unique(caseTable(uMaskAll, {'taskSeed', 'snrDb'}), 'rows', 'stable');
if isempty(keyTable) || height(keyTable) == 0
  return;
end
seed = keyTable.taskSeed;
snrDb = keyTable.snrDb;
angleDeltaDeg = NaN(height(keyTable), 1);
fdRefAbsErrDeltaHz = NaN(height(keyTable), 1);
fdRateAbsErrHzPerSec = NaN(height(keyTable), 1);
finalObjDelta = NaN(height(keyTable), 1);
fdRateMoveAbsHzPerSec = NaN(height(keyTable), 1);
cpUCandidateCount = NaN(height(keyTable), 1);
cpUSolveVariant = strings(height(keyTable), 1);
cpUBoundaryHit = false(height(keyTable), 1);
for iRow = 1:height(keyTable)
  kIdx = find(kMaskAll & caseTable.taskSeed == seed(iRow) & caseTable.snrDb == snrDb(iRow), 1, 'first');
  uIdx = find(uMaskAll & caseTable.taskSeed == seed(iRow) & caseTable.snrDb == snrDb(iRow), 1, 'first');
  if isempty(kIdx) || isempty(uIdx)
    continue;
  end
  angleDeltaDeg(iRow) = caseTable.sphericalAngleErrDeg(uIdx) - caseTable.sphericalAngleErrDeg(kIdx);
  fdRefAbsErrDeltaHz(iRow) = caseTable.fdRefAbsErrHz(uIdx) - caseTable.fdRefAbsErrHz(kIdx);
  fdRateAbsErrHzPerSec(iRow) = caseTable.fdRateAbsErrHzPerSec(uIdx);
  finalObjDelta(iRow) = caseTable.finalObj(uIdx) - caseTable.finalObj(kIdx);
  fdRateMoveAbsHzPerSec(iRow) = abs(caseTable.fdRateInitMoveHzPerSec(uIdx));
  cpUCandidateCount(iRow) = caseTable.candidateCount(uIdx);
  cpUSolveVariant(iRow) = caseTable.solveVariant(uIdx);
  cpUBoundaryHit(iRow) = caseTable.fdRateBoundaryHit(uIdx);
end
releaseCompareTable = table(seed, snrDb, angleDeltaDeg, fdRefAbsErrDeltaHz, ...
  fdRateAbsErrHzPerSec, finalObjDelta, fdRateMoveAbsHzPerSec, cpUCandidateCount, ...
  cpUSolveVariant, cpUBoundaryHit);
end

function rangeTable = localBuildRangeTable(repeatCell)
%LOCALBUILDRANGETABLE Record the truth-centered oracle boxes used by the replay.

rowList = repmat(struct('taskSeed', NaN, 'snrDb', NaN, 'toothStepHz', NaN, ...
  'truthFdRefHz', NaN, 'fdRangeLowHz', NaN, 'fdRangeHighHz', NaN, ...
  'truthFdRateHzPerSec', NaN, 'fdRateRangeLowHzPerSec', NaN, ...
  'fdRateRangeHighHzPerSec', NaN), numel(repeatCell), 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  rowList(iRepeat).taskSeed = repeatOut.taskSeed;
  rowList(iRepeat).snrDb = repeatOut.snrDb;
  rowList(iRepeat).toothStepHz = repeatOut.toothStepHz;
  rowList(iRepeat).truthFdRefHz = repeatOut.truthFdRefHz;
  rowList(iRepeat).fdRangeLowHz = repeatOut.fdRangeOracle(1);
  rowList(iRepeat).fdRangeHighHz = repeatOut.fdRangeOracle(2);
  rowList(iRepeat).truthFdRateHzPerSec = repeatOut.truthFdRateHzPerSec;
  rowList(iRepeat).fdRateRangeLowHzPerSec = repeatOut.fdRateRangeOracle(1);
  rowList(iRepeat).fdRateRangeHighHzPerSec = repeatOut.fdRateRangeOracle(2);
end
rangeTable = struct2table(rowList);
end


function objectiveProbeTable = localBuildObjectiveProbeTable(repeatCell)
%LOCALBUILDOBJECTIVEPROBETABLE Collect objective point probes.

rowList = repmat(localEmptyObjectiveProbeRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isfield(repeatOut, 'objectiveProbeRowList') && ~isempty(repeatOut.objectiveProbeRowList)
    rowList = [rowList; repeatOut.objectiveProbeRowList(:)]; %#ok<AGROW>
  end
end
objectiveProbeTable = struct2table(rowList(:));
end

function solveProbeTable = localBuildSolveProbeTable(repeatCell)
%LOCALBUILDSOLVEPROBETABLE Collect controlled solve probes.

rowList = repmat(localEmptySolveProbeRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isfield(repeatOut, 'solveProbeRowList') && ~isempty(repeatOut.solveProbeRowList)
    rowList = [rowList; repeatOut.solveProbeRowList(:)]; %#ok<AGROW>
  end
end
solveProbeTable = struct2table(rowList(:));
end

function pathProbeTable = localBuildPathProbeTable(repeatCell)
%LOCALBUILDPATHPROBETABLE Collect DoA-path objective probes.

rowList = repmat(localEmptyPathProbeRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isfield(repeatOut, 'pathProbeRowList') && ~isempty(repeatOut.pathProbeRowList)
    rowList = [rowList; repeatOut.pathProbeRowList(:)]; %#ok<AGROW>
  end
end
pathProbeTable = struct2table(rowList(:));
end

function diagnosisTable = localBuildDiagnosisTable(caseTable, objectiveProbeTable, solveProbeTable)
%LOCALBUILDDIAGNOSISTABLE Classify why MF does or does not align to its CRB.

rowList = repmat(localEmptyDiagnosisRow(), 0, 1);
if isempty(caseTable) || height(caseTable) == 0 || isempty(objectiveProbeTable) || height(objectiveProbeTable) == 0
  diagnosisTable = struct2table(rowList(:));
  return;
end
mfMask = startsWith(caseTable.displayName, "SS-MF");
mfRows = caseTable(mfMask, :);
for iCase = 1:height(mfRows)
  methodName = mfRows.displayName(iCase);
  snrDb = mfRows.snrDb(iCase);
  taskSeed = mfRows.taskSeed(iCase);
  objMask = objectiveProbeTable.displayName == methodName & ...
    objectiveProbeTable.snrDb == snrDb & objectiveProbeTable.taskSeed == taskSeed;
  objRows = objectiveProbeTable(objMask, :);
  solveMask = solveProbeTable.displayName == methodName & ...
    solveProbeTable.snrDb == snrDb & solveProbeTable.taskSeed == taskSeed;
  solveRows = solveProbeTable(solveMask, :);
  row = localEmptyDiagnosisRow();
  row.displayName = methodName;
  row.snrDb = snrDb;
  row.taskSeed = taskSeed;
  row.baselineAngleErrDeg = mfRows.finalAngleErrDeg(iCase);
  row.angleErrOverCrb = mfRows.angleErrOverSphericalCrb(iCase);
  row.fdRefErrOverCrb = mfRows.fdRefErrOverCrb(iCase);
  row.truthMinusFinalObj = localProbeObjDelta(objRows, "truthPoint", "finalPoint");
  row.staticSeedMinusFinalObj = localProbeObjDelta(objRows, "staticSeedPoint", "finalPoint");
  row.truthDoaFinalFdMinusFinalObj = localProbeObjDelta(objRows, "truthDoaFinalFdPoint", "finalPoint");
  row.finalDoaTruthFdMinusFinalObj = localProbeObjDelta(objRows, "finalDoaTruthFdPoint", "finalPoint");
  row.bestSolveTag = localFindBestSolveTag(solveRows);
  [row.bestSolveObjDelta, row.bestSolveAngleErrDeg] = localFindBestSolveMetric(solveRows);
  finalObjForTol = localProbeObj(objRows, "finalPoint");
  if ~isfinite(finalObjForTol)
    finalObjForTol = 1;
  end
  objTol = max(1e-6, 1e-9 * max(abs(finalObjForTol), 1));
  truthBetter = isfinite(row.truthMinusFinalObj) && row.truthMinusFinalObj < -objTol;
  truthDoaBetter = isfinite(row.truthDoaFinalFdMinusFinalObj) && row.truthDoaFinalFdMinusFinalObj < -objTol;
  bestSolveBetter = isfinite(row.bestSolveObjDelta) && row.bestSolveObjDelta < -objTol;
  if truthBetter && bestSolveBetter
    row.diagnosis = "basin-or-local-box-limited";
  elseif truthBetter && ~bestSolveBetter
    row.diagnosis = "solver-miss-even-with-good-objective-point";
  elseif ~truthBetter && truthDoaBetter
    row.diagnosis = "frequency-coupled-objective-mismatch";
  elseif ~truthBetter
    row.diagnosis = "objective-or-crb-statistical-check";
  else
    row.diagnosis = "undetermined";
  end
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
diagnosisTable = struct2table(rowList(:));
end

function obj = localProbeObj(objRows, probeTag)
%LOCALPROBEOBJ Return objective for a tag.

obj = NaN;
if isempty(objRows) || height(objRows) == 0
  return;
end
idx = find(objRows.probeTag == string(probeTag), 1, 'first');
if ~isempty(idx)
  obj = objRows.objective(idx);
end
end

function delta = localProbeObjDelta(objRows, tagA, tagB)
%LOCALPROBEOBJDELTA Return objective(tagA) - objective(tagB).

objA = localProbeObj(objRows, tagA);
objB = localProbeObj(objRows, tagB);
if isfinite(objA) && isfinite(objB)
  delta = objA - objB;
else
  delta = NaN;
end
end

function bestTag = localFindBestSolveTag(solveRows)
%LOCALFINDBESTSOLVETAG Return the best final-objective solve probe tag.

bestTag = "";
if isempty(solveRows) || height(solveRows) == 0
  return;
end
validMask = isfinite(solveRows.finalObj);
if ~any(validMask)
  return;
end
validIdx = find(validMask);
[~, relIdx] = min(solveRows.finalObj(validMask));
bestTag = solveRows.probeTag(validIdx(relIdx));
end

function [bestObjDelta, bestAngleErrDeg] = localFindBestSolveMetric(solveRows)
%LOCALFINDBESTSOLVEMETRIC Return best objective delta and associated angle.

bestObjDelta = NaN;
bestAngleErrDeg = NaN;
if isempty(solveRows) || height(solveRows) == 0
  return;
end
validMask = isfinite(solveRows.finalObj);
if ~any(validMask)
  return;
end
validIdx = find(validMask);
[~, relIdx] = min(solveRows.finalObj(validMask));
bestIdx = validIdx(relIdx);
bestObjDelta = solveRows.objMinusBaseline(bestIdx);
bestAngleErrDeg = solveRows.finalAngleErrDeg(bestIdx);
end

function plotData = localBuildPlotData(caseTable, aggregateTable, releaseCompareTable)
%LOCALBUILDPLOTDATA Store lightweight data needed to redraw diagnostic figures.

plotData = struct();
plotData.caseTable = caseTable;
plotData.aggregateTable = aggregateTable;
plotData.releaseCompareTable = releaseCompareTable;
end

function plotData = localPlotReplay(replayData)
%LOCALPLOTREPLAY Draw compact diagnostic figures from replayData.

plotData = replayData.plotData;
caseTable = replayData.caseTable;
aggregateTable = replayData.aggregateTable;
if isempty(caseTable) || height(caseTable) == 0
  return;
end

methodList = unique(caseTable.displayName, 'stable');
figure('Name', char(replayData.replayName + " angle metric diagnose"));
hold on;
for iMethod = 1:numel(methodList)
  mask = aggregateTable.displayName == methodList(iMethod);
  plot(aggregateTable.snrDb(mask), aggregateTable.angleRmseDeg(mask), '-o', ...
    'DisplayName', char(methodList(iMethod) + " RMSE"));
  plot(aggregateTable.snrDb(mask), aggregateTable.angleCrbSphericalMedianDeg(mask), '--', ...
    'DisplayName', char(methodList(iMethod) + " spherical CRB"));
end
grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Angle metric (deg)');
title('SS angle RMSE and projected angular CRB');
legend('Location', 'best');

figure('Name', char(replayData.replayName + " init final angle"));
hold on;
for iMethod = 1:numel(methodList)
  mask = caseTable.displayName == methodList(iMethod);
  scatter(caseTable.snrDb(mask), caseTable.initAngleErrDeg(mask), 'DisplayName', char(methodList(iMethod) + " init"));
  scatter(caseTable.snrDb(mask), caseTable.finalAngleErrDeg(mask), 'DisplayName', char(methodList(iMethod) + " final"));
end
grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Angle error (deg)');
title('Init vs final angle error floor check');
legend('Location', 'best');

figure('Name', char(replayData.replayName + " fdRef metric diagnose"));
hold on;
for iMethod = 1:numel(methodList)
  mask = aggregateTable.displayName == methodList(iMethod);
  plot(aggregateTable.snrDb(mask), aggregateTable.fdRefRmseHz(mask), '-o', ...
    'DisplayName', char(methodList(iMethod) + " fdRef RMSE"));
  crbPlotValue = aggregateTable.fdRefCrbMedianHz(mask);
  crbLabel = " fdRef CRB median";
  if ismember('fdRefCrbRmsHz', aggregateTable.Properties.VariableNames)
    crbPlotValue = aggregateTable.fdRefCrbRmsHz(mask);
    crbLabel = " fdRef CRB rms";
  end
  plot(aggregateTable.snrDb(mask), crbPlotValue, '--', ...
    'DisplayName', char(methodList(iMethod) + crbLabel));
end
grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('fdRef metric (Hz)');
title('SS fdRef RMSE and CRB floor check');
legend('Location', 'best');

plotData.caseTable = caseTable;
plotData.aggregateTable = aggregateTable;
end

function contextSummary = localBuildContextSummary(context)
%LOCALBUILDCONTEXTSUMMARY Build a lightweight context description.

contextSummary = struct();
contextSummary.carrierFreq = context.carrierFreq;
contextSummary.sampleRate = context.waveInfo.sampleRate;
contextSummary.refSatIdxGlobal = localGetFieldOrDefault(context, 'refSatIdxGlobal', NaN);
contextSummary.otherSatIdxGlobal = localGetFieldOrDefault(context, 'otherSatIdxGlobal', NaN);
contextSummary.periodicOffsetIdx = localGetFieldOrDefault(context, 'periodicOffsetIdx', []);
end

function repeatSlimCell = localStripRepeatCell(repeatCell)
%LOCALSTRIPREPEATCELL Keep lightweight repeat diagnostics in replayData.

repeatSlimCell = cell(size(repeatCell));
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  slim = struct();
  slim.taskSeed = repeatOut.taskSeed;
  slim.snrDb = repeatOut.snrDb;
  slim.truthFdRefHz = repeatOut.truthFdRefHz;
  slim.truthFdRateHzPerSec = repeatOut.truthFdRateHzPerSec;
  slim.toothStepHz = repeatOut.toothStepHz;
  slim.fdRangeOracle = repeatOut.fdRangeOracle;
  slim.fdRateRangeOracle = repeatOut.fdRateRangeOracle;
  slim.caseRowList = repeatOut.caseRowList;
  slim.objectiveProbeRowList = localGetFieldOrDefault(repeatOut, 'objectiveProbeRowList', repmat(localEmptyObjectiveProbeRow(), 0, 1));
  slim.solveProbeRowList = localGetFieldOrDefault(repeatOut, 'solveProbeRowList', repmat(localEmptySolveProbeRow(), 0, 1));
  slim.pathProbeRowList = localGetFieldOrDefault(repeatOut, 'pathProbeRowList', repmat(localEmptyPathProbeRow(), 0, 1));
  slim.repeatTotalMs = repeatOut.repeatTotalMs;
  slim.warningSeen = repeatOut.warningSeen;
  slim.warningId = repeatOut.warningId;
  slim.warningMessage = repeatOut.warningMessage;
  repeatSlimCell{iRepeat} = slim;
end
end

function lineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build HTML-ready key metric lines for notifications.

aggregateTable = replayData.aggregateTable;
if isempty(aggregateTable) || height(aggregateTable) == 0
  lineList = strings(0, 1);
  return;
end
snrPick = max(aggregateTable.snrDb);
rowK = localFindAggregateRow(aggregateTable, "SS-MF-CP-K", snrPick);
rowU = localFindAggregateRow(aggregateTable, "SS-MF-CP-U", snrPick);
lineList = strings(0, 1);
if ~isempty(rowK)
  lineList(end + 1, 1) = sprintf('CP-K @ %.1f dB: angle/CRB=<code>%.3f</code>, fdRef/CRB=<code>%.3f</code>', ...
    rowK.snrDb(1), rowK.angleRmseOverSphericalCrb(1), localGetTableScalar(rowK, 'fdRefNormRmse', rowK.fdRefRmseOverCrb(1)));
end
if ~isempty(rowU)
  lineList(end + 1, 1) = sprintf('CP-U @ %.1f dB: angle/CRB=<code>%.3f</code>, fdRef/CRB=<code>%.3f</code>', ...
    rowU.snrDb(1), rowU.angleRmseOverSphericalCrb(1), localGetTableScalar(rowU, 'fdRefNormRmse', rowU.fdRefRmseOverCrb(1)));
end
healthTable = localGetFieldOrDefault(replayData, 'healthAggregateTable', table());
if ~isempty(healthTable) && height(healthTable) > 0
  healthRowK = localFindAggregateRow(healthTable, "SS-MF-CP-K", snrPick);
  if ~isempty(healthRowK) && ismember('keepRate', healthTable.Properties.VariableNames)
    lineList(end + 1, 1) = sprintf('health keep @ %.1f dB: CP-K=<code>%.2f</code>', ...
      healthRowK.snrDb(1), healthRowK.keepRate(1));
  end
end
diagnosisTable = localGetFieldOrDefault(replayData, 'diagnosisTable', table());
if ~isempty(diagnosisTable) && height(diagnosisTable) > 0 && ismember('diagnosis', diagnosisTable.Properties.VariableNames)
  diagPick = diagnosisTable(diagnosisTable.snrDb == min(diagnosisTable.snrDb), :);
  if height(diagPick) > 0
    lineList(end + 1, 1) = "diagnosis sample: <code>" + localHtmlEscape(string(diagPick.diagnosis(1))) + "</code>";
  end
end
lineList(end + 1, 1) = sprintf('seeds=<code>%d</code>, snr points=<code>%d</code>', ...
  replayData.config.numRepeat, numel(replayData.config.snrDbList));
end

function value = localGetTableScalar(rowTable, variableName, defaultValue)
%LOCALGETTABLESCALAR Read one scalar table value with fallback.

value = defaultValue;
if istable(rowTable) && height(rowTable) > 0 && ismember(variableName, rowTable.Properties.VariableNames)
  candidate = rowTable.(variableName)(1);
  if isnumeric(candidate) || islogical(candidate)
    value = candidate;
  end
end
end

function row = localFindAggregateRow(aggregateTable, displayName, snrDb)
%LOCALFINDAGGREGATEROW Return one aggregate row for compact notification.

idx = find(aggregateTable.displayName == string(displayName) & aggregateTable.snrDb == snrDb, 1, 'first');
if isempty(idx)
  row = [];
else
  row = aggregateTable(idx, :);
end
end


function row = localEmptyObjectiveProbeRow()
%LOCALEMPTYOBJECTIVEPROBEROW Return one typed objective probe row.

row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'probeTag', "", 'objective', NaN, 'finalObjective', NaN, ...
  'objMinusFinal', NaN, 'angleErrDeg', NaN, 'fdRefErrHz', NaN, ...
  'fdRateErrHzPerSec', NaN, 'evalOk', false, 'evalMessage', "");
end

function row = localEmptySolveProbeRow()
%LOCALEMPTYSOLVEPROBEROW Return one typed controlled solve probe row.

row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'probeTag', "", 'probeGroup', "", 'doaHalfWidthLatDeg', NaN, ...
  'doaHalfWidthLonDeg', NaN, 'finalAngleErrDeg', NaN, ...
  'finalFdRefErrHz', NaN, 'finalFdRateErrHzPerSec', NaN, ...
  'finalObj', NaN, 'objMinusBaseline', NaN, 'iterations', NaN, ...
  'exitflag', NaN, 'solveVariant', "", 'candidateCount', NaN, ...
  'boundaryHit', false, 'evalOk', false, 'message', "");
end

function row = localEmptyPathProbeRow()
%LOCALEMPTYPATHPROBEROW Return one typed DoA path probe row.

row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'pathTag', "", 'alpha', NaN, 'objective', NaN, 'objMinusStart', NaN, ...
  'angleErrDeg', NaN, 'fdRefErrHz', NaN, 'fdRateErrHzPerSec', NaN, ...
  'evalOk', false);
end

function row = localEmptyDiagnosisRow()
%LOCALEMPTYDIAGNOSISROW Return one typed diagnosis decision row.

row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'baselineAngleErrDeg', NaN, 'angleErrOverCrb', NaN, ...
  'fdRefErrOverCrb', NaN, 'truthMinusFinalObj', NaN, ...
  'staticSeedMinusFinalObj', NaN, 'truthDoaFinalFdMinusFinalObj', NaN, ...
  'finalDoaTruthFdMinusFinalObj', NaN, 'bestSolveTag', "", ...
  'bestSolveObjDelta', NaN, 'bestSolveAngleErrDeg', NaN, 'diagnosis', "");
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for diagnostic case metrics.

row = struct('displayName', "", 'satMode', "", 'frameMode', "", 'phaseMode', "", ...
  'fdRateMode', "", 'mfInitMode', "", 'snrDb', NaN, 'taskSeed', NaN, 'numFrame', NaN, ...
  'frameIntvlSec', NaN, 'windowSec', NaN, 'toothStepHz', NaN, ...
  'truthLatDeg', NaN, 'truthLonDeg', NaN, 'initLatDeg', NaN, 'initLonDeg', NaN, ...
  'finalLatDeg', NaN, 'finalLonDeg', NaN, 'latErrDeg', NaN, 'lonErrDeg', NaN, ...
  'initAngleErrDeg', NaN, 'finalAngleErrDeg', NaN, 'sphericalAngleErrDeg', NaN, ...
  'finalMinusInitAngleDeg', NaN, 'searchHalfWidthLatDeg', NaN, ...
  'searchHalfWidthLonDeg', NaN, 'crbLatStdDeg', NaN, 'crbLonStdDeg', NaN, ...
  'crbTraceStdDeg', NaN, 'crbSphericalApproxStdDeg', NaN, ...
  'crbAngleMetricVarRad2', NaN, 'crbMetricTraceRatio', NaN, ...
  'angleErrOverSphericalCrb', NaN, 'angleErrOverTraceCrb', NaN, ...
  'fdRefTrueHz', NaN, 'initFdRefHz', NaN, 'finalFdRefHz', NaN, ...
  'fdRefErrHz', NaN, 'fdRefAbsErrHz', NaN, 'fdRefCrbStdHz', NaN, ...
  'fdRefErrOverCrb', NaN, 'fdRefSignedErrOverCrb', NaN, ...
  'fdRefInitMoveHz', NaN, 'fdRangeLowHz', NaN, ...
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
  'healthResolved', false, 'trimKeep', false, 'rejectReason', "", ...
  'warningFlag', false);
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return a typed empty row for aggregate metrics.

row = struct('displayName', "", 'frameMode', "", 'fdRateMode', "", 'snrDb', NaN, ...
  'numRepeat', NaN, 'angleRmseDeg', NaN, 'angleMedianDeg', NaN, ...
  'angleMaxDeg', NaN, 'angleCrbSphericalMedianDeg', NaN, ...
  'angleCrbTraceMedianDeg', NaN, 'angleRmseOverSphericalCrb', NaN, ...
  'angleRmseOverTraceCrb', NaN, 'fdRefRmseHz', NaN, ...
  'fdRefCrbMedianHz', NaN, 'fdRefRmseOverCrb', NaN, ...
  'fdRefCrbMeanHz', NaN, 'fdRefCrbRmsHz', NaN, ...
  'fdRefRmseOverMeanCrbStd', NaN, 'fdRefNormRmse', NaN, ...
  'fdRefMseOverMeanCrb', NaN, 'fdRefBiasHz', NaN, ...
  'fdRefBiasOverMeanCrbStd', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'initAngleMedianDeg', NaN, ...
  'finalMinusInitMedianDeg', NaN, 'fdRefInitMoveMedianHz', NaN, ...
  'fdRateInitMoveMedianHzPerSec', NaN, 'objectiveImproveMedian', NaN, ...
  'boundaryHitRate', NaN, 'candidateCountMedian', NaN, ...
  'iterationsMedian', NaN, 'firstOrderOptMedian', NaN);
end

function row = localEmptyFilteredAggregateRow()
%LOCALEMPTYFILTEREDAGGREGATEROW Return a typed empty row for filtered aggregates.

row = localEmptyAggregateRow();
row.filterName = "";
row.rawCount = NaN;
row.keepCount = NaN;
row.keepRate = NaN;
row.rejectRate = NaN;
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

function parallelOpt = localBuildSerialInnerParallelOpt()
%LOCALBUILDSERIALINNERPARALLELOPT Disable inner parallelism for diagnostic replay consistency.

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
end

function context = localDisableSubsetBankForDiagnose(context)
%LOCALDISABLESUBSETBANKFORDIAGNOSE Remove subset fixtures from this local metric replay.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if numFrame < 2
  error('replayMfSsMleCrbMetricDiagnose:InvalidFrameCount', ...
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

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide only when the denominator is finite and nonzero.

if isfinite(numerator) && isfinite(denominator) && abs(denominator) > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
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

function textValue = localJoinStringList(value)
%LOCALJOINSTRINGLIST Join a short list of variant tags for display.

value = string(value(:));
value = value(strlength(value) > 0);
if isempty(value)
  textValue = "";
else
  textValue = strjoin(cellstr(value), ',');
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


function textValue = localHtmlEscape(textValue)
%LOCALHTMLESCAPE Escape dynamic text for Telegram HTML.

textValue = string(textValue);
textValue = replace(textValue, "&", "&amp;");
textValue = replace(textValue, "<", "&lt;");
textValue = replace(textValue, ">", "&gt;");
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

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field from a struct with a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
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
