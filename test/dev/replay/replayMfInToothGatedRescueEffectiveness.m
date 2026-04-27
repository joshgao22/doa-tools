% replayMfInToothGatedRescueEffectiveness
% Purpose: run a small Monte Carlo validation of no-truth gated in-tooth
% rescue candidates. The replay assumes the fdRef tooth is already correct by
% using a truth-centered half-tooth fd range, then compares disabled and gated
% wide/single-MF basin-entry banks. Offline case labels may use truth, but the
% rescue trigger and selected candidate use only default-estimate diagnostics.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted repeat runs from per-repeat files under the repo-root tmp.
% saveSnapshot=true saves only replayData via saveExpSnapshot.

clear; close all; clc;

%% Replay configuration
replayName = "replayMfInToothGatedRescueEffectiveness";

snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % First seed of the small replay batch.
numRepeat = 100;                    % Number of consecutive seeds to replay.
saveSnapshot = true;                % true saves the lightweight replayData only.
optVerbose = false;                 % true enables compact estimator / branch trace.
checkpointEnable = true;             % true enables per-repeat checkpoint/resume under repo-root tmp/.
oracleFdHalfToothFraction = 0.49;   % Half-width fraction of the 1/T_f tooth step.
oracleFdRateHalfWidthHzPerSec = 1000; % Local truth-centered fdRate half-width.
staticLocalDoaHalfWidthDeg = [0.002; 0.002]; % Local DoA box around the static seed.
staticWideDoaHalfWidthDeg = [0.010; 0.010];  % Wider same-tooth DoA box for rescue center.
truthLocalDoaHalfWidthDeg = [0.002; 0.002];  % Truth DoA oracle box used only for offline labeling.
coarseDoaStepDegList = [-0.006; -0.004; -0.002; 0; 0.002; 0.004; 0.006]; % Implementable rescue DoA offsets.
fdRefHealthyAbsHz = 1;              % Offline fdRef health threshold.
fdRateHealthyAbsHzPerSec = 50;      % Offline fdRate health threshold.
coherenceHealthyThreshold = 0.95;   % Offline selected non-ref coherence recovery threshold.
truthDoaGapThresholdDeg = 1e-3;     % Offline default gap to truth-DoA oracle threshold.
angleGainThresholdDeg = 5e-4;       % Offline rescue gain threshold.
damageThresholdDeg = 5e-4;          % Offline damage threshold.
gatedRescueCoherenceThreshold = 0.20; % No-truth trigger: default non-ref coherence has clearly collapsed.
gatedRescuePhaseResidThresholdRad = 1.0; % No-truth trigger: default non-ref phase residual is large.
gatedRescueGateVersion = "coherence-or-phase-v2"; % Checkpoint signature for the no-truth trigger logic.

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

%% Build context and flow options
config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.saveSnapshot = saveSnapshot;
config.optVerbose = optVerbose;
config.checkpointEnable = checkpointEnable;
config.checkpointResume = true;
config.checkpointCleanupOnSuccess = true;
config.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
config.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
config.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg(:);
config.staticWideDoaHalfWidthDeg = staticWideDoaHalfWidthDeg(:);
config.truthLocalDoaHalfWidthDeg = truthLocalDoaHalfWidthDeg(:);
config.coarseDoaStepDegList = coarseDoaStepDegList(:);
config.fdRefHealthyAbsHz = fdRefHealthyAbsHz;
config.fdRateHealthyAbsHzPerSec = fdRateHealthyAbsHzPerSec;
config.coherenceHealthyThreshold = coherenceHealthyThreshold;
config.truthDoaGapThresholdDeg = truthDoaGapThresholdDeg;
config.angleGainThresholdDeg = angleGainThresholdDeg;
config.damageThresholdDeg = damageThresholdDeg;
config.gatedRescueCoherenceThreshold = gatedRescueCoherenceThreshold;
config.gatedRescuePhaseResidThresholdRad = gatedRescuePhaseResidThresholdRad;
config.gatedRescueGateVersion = string(gatedRescueGateVersion);

fprintf('Running %s ...\n', char(replayName));
fprintf('  repeats                         : %d\n', config.numRepeat);
fprintf('  snr (dB)                        : %.2f\n', config.snrDb);
fprintf('  base seed                       : %d\n', config.baseSeed);
fprintf('  repeat mode                     : %s\n', 'parfor-auto');
fprintf('  save snapshot                   : %d\n', config.saveSnapshot);
fprintf('  checkpoint                      : %d\n', config.checkpointEnable);
fprintf('  fd oracle half-tooth fraction   : %.3f\n', config.oracleFdHalfToothFraction);
fprintf('  fdRate oracle half-width        : %.2f Hz/s\n', config.oracleFdRateHalfWidthHzPerSec);
fprintf('  coarse rescue DoA steps         : %s deg\n', mat2str(config.coarseDoaStepDegList(:).'));
fprintf('  gated rescue coherence threshold: %.3f\n', config.gatedRescueCoherenceThreshold);
fprintf('  gated rescue phase threshold    : %.3f rad\n', config.gatedRescuePhaseResidThresholdRad);
fprintf('  gated rescue gate version       : %s\n', char(config.gatedRescueGateVersion));

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
contextOpt = struct( ...
  'baseSeed', config.baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt);
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false));
methodList = localBuildMethodList(config);

checkpointOpt = localBuildCheckpointOpt(replayName, config);
checkpointRunDir = "";
if config.checkpointEnable
  checkpointRunDir = checkpointOpt.runDir;
  config.checkpointRunDir = string(checkpointRunDir);
  fprintf('  checkpoint resume               : %d\n', config.checkpointResume);
  fprintf('  checkpoint dir                  : %s\n', checkpointRunDir);
end

fprintf('[%s] Build shared dynamic context.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
context = buildDynamicDualSatEciContext(contextOpt);
context = localDisableSubsetBankForInTooth(context);
fprintf('[%s] Shared dynamic context built.\n', char(datetime('now', 'Format', 'HH:mm:ss')));

%% Run replay batch
repeatCell = cell(numRepeat, 1);
useParfor = localCanUseParfor() && numRepeat > 1;
runState = struct();
try
  if config.checkpointEnable
    taskGrid = localBuildRepeatTaskGrid(config.seedList);
    sharedData = struct('context', context, 'flowOpt', flowOpt, 'methodList', methodList, 'config', config);
    numDoneTask = localCountCheckpointTaskFile(fullfile(checkpointRunDir, 'task'), numRepeat);
    numTodoTask = numRepeat - numDoneTask;
    tracker = localCreateProgressTracker(sprintf('%s repeat batch (%s, resume %d/%d)', ...
      char(replayName), localModeText(useParfor), numDoneTask, numRepeat), numTodoTask, false);
    checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, tracker, useParfor);
    runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, @localCheckpointTaskRunner, checkpointRunnerOpt);
    localCloseProgressTracker(tracker);
    if ~runState.isComplete
      error('replayMfInToothGatedRescueEffectiveness:CheckpointIncomplete', ...
        'Checkpoint replay returned an incomplete repeat batch.');
    end
    repeatCell = runState.resultCell;
  else
    tracker = localCreateProgressTracker(sprintf('%s repeat batch (%s)', char(replayName), localModeText(useParfor)), numRepeat, useParfor);
    if useParfor
      progressQueue = tracker.queue;
      parfor iRepeat = 1:numRepeat
        repeatCell{iRepeat} = localRunOneRepeat(iRepeat, config.seedList(iRepeat), context, flowOpt, methodList, config);
        if ~isempty(progressQueue)
          send(progressQueue, iRepeat);
        end
      end
    else
      for iRepeat = 1:numRepeat
        repeatCell{iRepeat} = localRunOneRepeat(iRepeat, config.seedList(iRepeat), context, flowOpt, methodList, config);
        localAdvanceProgressTracker(tracker);
      end
    end
    localCloseProgressTracker(tracker);
  end
catch ME
  if exist('tracker', 'var')
    localCloseProgressTracker(tracker);
  end
  if config.checkpointEnable
    fprintf('Replay failed. Checkpoint artifacts kept at: %s\n', checkpointRunDir);
  end
  rethrow(ME);
end

methodTable = localBuildMethodTable(repeatCell);
methodAggregateTable = localBuildMethodAggregateTable(methodTable);
candidateProbeTable = localBuildCandidateProbeTable(repeatCell);
rescueBankDecisionTable = localBuildRescueEffectCaseTable(repeatCell);
rescueEffectCaseTable = localSelectGatedEffectCaseTable(rescueBankDecisionTable);
triggerReasonTable = localBuildTriggerReasonTable(rescueEffectCaseTable);
rescueEffectAggregateTable = localBuildRescueEffectAggregateTable(rescueBankDecisionTable);
rescueEffectVerdictTable = localBuildRescueEffectVerdictTable(rescueEffectAggregateTable);
rangeTable = localBuildRangeTable(repeatCell);
timingTable = localBuildTimingTable(repeatCell);
timingAggregateTable = localBuildTimingAggregateTable(timingTable);
plotData = localBuildPlotData(rescueBankDecisionTable, rescueEffectCaseTable, triggerReasonTable);

%% Data storage
replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, methodList);
replayData.methodTable = methodTable;
replayData.methodAggregateTable = methodAggregateTable;
replayData.candidateProbeTable = candidateProbeTable;
replayData.rescueBankDecisionTable = rescueBankDecisionTable;
replayData.rescueEffectCaseTable = rescueEffectCaseTable;
replayData.triggerReasonTable = triggerReasonTable;
replayData.rescueEffectAggregateTable = rescueEffectAggregateTable;
replayData.rescueEffectVerdictTable = rescueEffectVerdictTable;
replayData.rangeTable = rangeTable;
replayData.timingTable = timingTable;
replayData.timingAggregateTable = timingAggregateTable;
replayData.plotData = plotData;
replayData.methodList = localStripMethodList(methodList);
if config.checkpointEnable
  replayData.checkpointSummary = localBuildCheckpointSummary(runState);
end

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
if ~exist('replayData', 'var') || ~isstruct(replayData)
  error('Replay data is missing. Run the replay batch sections or load a snapshot containing replayData.');
end
fprintf('Running replayMfInToothGatedRescueEffectiveness ...\n');
fprintf('  repeats                         : %d\n', numel(unique(replayData.rescueEffectCaseTable.taskSeed)));
fprintf('  snr (dB)                        : %.2f\n', replayData.config.snrDb);
fprintf('  base seed                       : %d\n', replayData.config.baseSeed);
if isfield(replayData.config, 'checkpointEnable')
  fprintf('  checkpoint                      : %d\n', replayData.config.checkpointEnable);
end
if isfield(replayData.config, 'checkpointRunDir')
  fprintf('  checkpoint dir                  : %s\n', char(replayData.config.checkpointRunDir));
end
fprintf('\n========== In-tooth method aggregate ==========\n');
disp(replayData.methodAggregateTable);
fprintf('\n========== Gated rescue aggregate ==========\n');
disp(replayData.rescueEffectAggregateTable);
if isfield(replayData, 'rescueEffectVerdictTable') && ~isempty(replayData.rescueEffectVerdictTable)
  fprintf('\n========== Gated rescue verdict ==========\n');
  disp(replayData.rescueEffectVerdictTable);
end
fprintf('\n========== Trigger reason table ==========\n');
localDispTablePreview(replayData.triggerReasonTable, 8);
if isfield(replayData, 'rescueBankDecisionTable')
  fprintf('\n========== Rescue bank decision table ==========\n');
  localDispTablePreview(replayData.rescueBankDecisionTable, 8);
end
fprintf('\n========== Gated rescue case table ==========\n');
localDispTablePreview(replayData.rescueEffectCaseTable, 8);
if ~isempty(replayData.timingAggregateTable)
  fprintf('\n========== Runtime timing summary ==========\n');
  disp(replayData.timingAggregateTable);
end
localPlotReplay(replayData.plotData);

%% Local helpers

function checkpointOpt = localBuildCheckpointOpt(replayName, config)
%LOCALBUILDCHECKPOINTOPT Build stable per-repeat checkpoint options for this replay.

checkpointOpt = struct('enable', config.checkpointEnable);
if ~config.checkpointEnable
  checkpointOpt.runDir = "";
  return;
end
checkpointRunKey = localBuildCheckpointRunKey(config);
checkpointOutputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.resume = config.checkpointResume;
checkpointOpt.runName = string(replayName);
checkpointOpt.runKey = checkpointRunKey;
checkpointOpt.outputRoot = checkpointOutputRoot;
checkpointOpt.runDir = fullfile(checkpointOutputRoot, char(replayName), char(checkpointRunKey));
checkpointOpt.meta = struct( ...
  'snrDb', config.snrDb, ...
  'seedList', reshape(config.seedList, 1, []), ...
  'gateVersion', char(config.gatedRescueGateVersion), ...
  'coherenceThreshold', config.gatedRescueCoherenceThreshold, ...
  'phaseResidThresholdRad', config.gatedRescuePhaseResidThresholdRad, ...
  'oracleFdHalfToothFraction', config.oracleFdHalfToothFraction, ...
  'oracleFdRateHalfWidthHzPerSec', config.oracleFdRateHalfWidthHzPerSec);
end

function runKey = localBuildCheckpointRunKey(config)
%LOCALBUILDCHECKPOINTRUNKEY Build a checkpoint key that changes with gate semantics.

seedList = reshape(double(config.seedList), 1, []);
runKey = string(sprintf('seed%dto%d_rep%d_snr%.2f_%s_coh%.3f_phase%.3f', ...
  seedList(1), seedList(end), numel(seedList), config.snrDb, ...
  char(config.gatedRescueGateVersion), config.gatedRescueCoherenceThreshold, ...
  config.gatedRescuePhaseResidThresholdRad));
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve the repository root from the replay script location.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
end

function taskGrid = localBuildRepeatTaskGrid(seedList)
%LOCALBUILDREPEATTASKGRID Build one checkpoint task per repeat seed.

numTask = numel(seedList);
taskGrid = repmat(struct('taskIndex', 0, 'repeatIndex', 0, 'taskSeed', 0), numTask, 1);
for iTask = 1:numTask
  taskGrid(iTask).taskIndex = iTask;
  taskGrid(iTask).repeatIndex = iTask;
  taskGrid(iTask).taskSeed = seedList(iTask);
end
end

function repeatOut = localCheckpointTaskRunner(taskInfo, sharedData)
%LOCALCHECKPOINTTASKRUNNER Run one checkpointed repeat task.

repeatOut = localRunOneRepeat(taskInfo.repeatIndex, taskInfo.taskSeed, ...
  sharedData.context, sharedData.flowOpt, sharedData.methodList, sharedData.config);
end

function checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, tracker, useParfor)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only options accepted by the checkpoint runner.

checkpointRunnerOpt = struct();
checkpointRunnerOpt.runName = checkpointOpt.runName;
checkpointRunnerOpt.runKey = checkpointOpt.runKey;
checkpointRunnerOpt.outputRoot = checkpointOpt.outputRoot;
checkpointRunnerOpt.useParfor = useParfor;
checkpointRunnerOpt.resume = checkpointOpt.resume;
checkpointRunnerOpt.meta = checkpointOpt.meta;
checkpointRunnerOpt.progressFcn = @(step) localAdvanceProgressByStep(tracker, step);
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

function checkpointSummary = localBuildCheckpointSummary(runState)
%LOCALBUILDCHECKPOINTSUMMARY Store lightweight checkpoint status only.

checkpointSummary = struct();
checkpointSummary.runDir = string(localGetFieldOrDefault(runState, 'runDir', ""));
checkpointSummary.numTask = localGetFieldOrDefault(runState, 'numTask', 0);
checkpointSummary.numDone = localGetFieldOrDefault(runState, 'numDone', 0);
checkpointSummary.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
checkpointSummary.cleanedOnSuccess = false;
end


function context = localDisableSubsetBankForInTooth(context)
%LOCALDISABLESUBSETBANKFORINTOOTH Skip curated and random subset fixtures.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function methodList = localBuildMethodList(config)
%LOCALBUILDMETHODLIST Define the compact in-tooth methods needed by the rescue replay.

methodList = repmat(localMakeMethod("", "", "", false, "", "", [], false, "", "", false), 0, 1);
methodList(end + 1, 1) = localMakeMethod("ss-mf-cp-u-in-tooth", "single-center", "ref", ...
  false, "static", "truthFreq", config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal", true);
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-in-tooth", "default-target", "ms", ...
  false, "static", "truthFreq", config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal", true);
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-wide-doa-in-tooth", "wide-center", "ms", ...
  false, "static", "truthFreq", config.staticWideDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal", true);
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-truth-doa-oracle", "truth-doa-label", "ms", ...
  false, "truth", "truthFreq", config.truthLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal", true);
end

function method = localMakeMethod(label, family, viewMode, isKnownRate, doaSeedMode, initMode, doaHalfWidthDeg, freezeDoa, fdRangeMode, fdRateRangeMode, isOracle)
%LOCALMAKEMETHOD Construct one method descriptor with a stable field set.

method = struct();
method.label = string(label);
method.family = string(family);
method.viewMode = string(viewMode);
method.isKnownRate = logical(isKnownRate);
method.doaSeedMode = string(doaSeedMode);
method.initMode = string(initMode);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
method.freezeDoa = logical(freezeDoa);
method.fdRangeMode = string(fdRangeMode);
method.fdRateRangeMode = string(fdRateRangeMode);
method.isOracle = logical(isOracle);
end

function repeatOut = localRunOneRepeat(iRepeat, taskSeed, context, flowOpt, methodList, config)
%LOCALRUNONEREPEAT Build one repeat, evaluate rescue candidates, and pack compact outputs.

lastwarn('', '');
tRepeat = tic;
tBuildData = tic;
repeatData = buildDynamicRepeatData(context, config.snrDb, taskSeed);
buildDataMs = 1000 * toc(tBuildData);
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('replayMfInToothGatedRescueEffectiveness:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for the in-tooth rescue replay.');
end
fdHalfWidthHz = config.oracleFdHalfToothFraction * toothStepHz;
fdRangeOracle = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeOracle = truthFdRateHzPerSec + config.oracleFdRateHalfWidthHzPerSec * [-1, 1];

tStaticBundle = tic;
staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fixture.fdRange, truth, context.otherSatIdxGlobal, config.optVerbose, ...
  flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);
staticBundleMs = 1000 * toc(tStaticBundle);

rowCell = cell(numel(methodList) + 2, 1);
rowCell{1} = localBuildSummaryRow(iRepeat, taskSeed, "ss-sf-static", "static-baseline", ...
  "single", "static", false, "static", "static", "default", "none", NaN, false, ...
  staticBundle.caseStaticRefOnly, truth, toothStepHz, NaN);
rowCell{2} = localBuildSummaryRow(iRepeat, taskSeed, "ms-sf-static", "static-baseline", ...
  "multi", "static", false, "static", "static", "default", "none", NaN, false, ...
  staticBundle.caseStaticMs, truth, toothStepHz, NaN);

methodCaseCell = cell(numel(methodList), 1);
dynamicMethodWallTimeMs = nan(numel(methodList), 1);
for iMethod = 1:numel(methodList)
  method = methodList(iMethod);
  tMethod = tic;
  caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, ...
    context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, config.optVerbose, ...
    flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz);
  dynamicMethodWallTimeMs(iMethod) = 1000 * toc(tMethod);
  methodCaseCell{iMethod} = caseUse;
  rowCell{iMethod + 2} = localBuildSummaryRow(iRepeat, taskSeed, method.label, method.family, ...
    localInferSatMode(method.viewMode), "dynamic-cp", method.isKnownRate, method.doaSeedMode, ...
    method.initMode, method.fdRangeMode, method.fdRateRangeMode, localMaxAbs(method.doaHalfWidthDeg), ...
    method.freezeDoa, caseUse, truth, toothStepHz, dynamicMethodWallTimeMs(iMethod));
end
candidateProbeTable = localBuildCandidateProbeForRepeat(iRepeat, taskSeed, fixture, staticBundle, truth, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, flowOpt, methodList, ...
  methodCaseCell, fdRangeOracle, fdRateRangeOracle, toothStepHz, config);
repeatTotalMs = 1000 * toc(tRepeat);
rescueCaseTable = localBuildRescueCaseTableForRepeat(taskSeed, candidateProbeTable, config, repeatTotalMs);

[warningMessage, warningId] = lastwarn;
repeatOut = struct();
repeatOut.iRepeat = iRepeat;
repeatOut.taskSeed = taskSeed;
repeatOut.snrDb = config.snrDb;
repeatOut.toothStepHz = toothStepHz;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.fdRangeOracle = fdRangeOracle;
repeatOut.fdRateRangeOracle = fdRateRangeOracle;
repeatOut.methodSummary = struct2table([rowCell{:}].');
repeatOut.candidateProbeTable = candidateProbeTable;
repeatOut.rescueCaseTable = rescueCaseTable;
repeatOut.timing = struct( ...
  'buildDataMs', buildDataMs, ...
  'staticBundleMs', staticBundleMs, ...
  'dynamicMethodTotalMs', sum(dynamicMethodWallTimeMs, 'omitnan'), ...
  'dynamicMethodMaxMs', max(dynamicMethodWallTimeMs, [], 'omitnan'), ...
  'repeatTotalMs', repeatTotalMs);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz)
%LOCALRUNDYNAMICMETHOD Run one in-tooth dynamic method descriptor.

[viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle);
fdRangeUse = fixture.fdRange;
fdRateRangeUse = fixture.fdRateRange;
if method.fdRangeMode == "truthHalfTooth"
  fdRangeUse = fdRangeOracle;
end
if method.fdRateRangeMode == "truthLocal"
  fdRateRangeUse = fdRateRangeOracle;
end
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
seedDoaParam = localResolveDoaSeed(method, staticCaseUse, truth);

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = method.freezeDoa;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);

switch method.initMode
  case "truthFreq"
    initParam = localBuildTruthFreqInitParam(seedDoaParam, truthFdRefHz, truthFdRateHzPerSec, method.isKnownRate);
    initOverride = struct('startTag', method.label, 'initParam', initParam, ...
      'initDoaParam', seedDoaParam(:), 'initDoaHalfWidth', method.doaHalfWidthDeg(:), ...
      'freezeDoa', method.freezeDoa);
  otherwise
    error('replayMfInToothGatedRescueEffectiveness:UnknownInitMode', ...
      'Unknown method initMode "%s".', method.initMode);
end

displayPrefix = "MS";
if method.viewMode == "ref"
  displayPrefix = "SS";
end
displayName = displayPrefix + "-MF-" + method.label;
caseUse = runDynamicDoaDopplerCase(displayName, localInferSatMode(method.viewMode), viewUse, truth, ...
  pilotWave, carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, optVerbose, ...
  dynOpt, method.isKnownRate, debugTruthUse, initOverride);
caseUse.oracleMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, 'fdRateRangeUse', fdRateRangeUse);
end

function [viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle)
%LOCALRESOLVEMETHODVIEW Resolve view, debug truth, and static seed for one method.

switch method.viewMode
  case "ref"
    viewUse = fixture.viewRefOnly;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());
    staticCaseUse = staticBundle.caseStaticRefOnly;
  case "ms"
    viewUse = fixture.viewMs;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());
    staticCaseUse = staticBundle.caseStaticMs;
  otherwise
    error('replayMfInToothGatedRescueEffectiveness:UnknownViewMode', ...
      'Unknown method viewMode "%s".', method.viewMode);
end
end

function satMode = localInferSatMode(viewMode)
%LOCALINFERSATMODE Convert a view-mode label to runDynamicDoaDopplerCase sat mode.

if string(viewMode) == "ref"
  satMode = "single";
else
  satMode = "multi";
end
end

function seedDoaParam = localResolveDoaSeed(method, staticCase, truth)
%LOCALRESOLVEDOASEED Choose the DoA seed for one replay-only method.

switch method.doaSeedMode
  case "static"
    seedDoaParam = staticCase.estResult.doaParamEst(:);
  case "truth"
    seedDoaParam = reshape(truth.latlonTrueDeg, [], 1);
  otherwise
    error('replayMfInToothGatedRescueEffectiveness:UnknownDoaSeedMode', ...
      'Unknown DoA seed mode "%s".', method.doaSeedMode);
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

function candidateProbeTable = localBuildCandidateProbeForRepeat(iRepeat, taskSeed, fixture, staticBundle, truth, pilotWave, carrierFreq, sampleRate, flowOpt, methodList, methodCaseCell, fdRangeOracle, fdRateRangeOracle, toothStepHz, config)
%LOCALBUILDCANDIDATEPROBEFORREPEAT Evaluate implementable rescue grids without solver adoption.

candidateProbeTable = table();
defaultIdx = find([methodList.label] == "ms-mf-cp-u-in-tooth", 1, 'first');
if isempty(defaultIdx) || defaultIdx > numel(methodCaseCell)
  return;
end
defaultMethod = methodList(defaultIdx);
defaultCase = methodCaseCell{defaultIdx};
defaultOptVar = localResolveOptVarFromCase(defaultCase);
if numel(defaultOptVar) < 4
  return;
end
[viewUse, ~, staticCaseUse] = localResolveMethodView(defaultMethod, fixture, staticBundle);
seedDoaParam = localResolveDoaSeed(defaultMethod, staticCaseUse, truth);
dynOpt = localBuildDynOptForProbe(defaultMethod, fixture, flowOpt, seedDoaParam);
[model, ~, ~] = buildDoaDopplerMfModel(viewUse.sceneSeq, viewUse.rxSigMf, pilotWave, ...
  carrierFreq, sampleRate, viewUse.doaGrid, fdRangeOracle, fdRateRangeOracle, dynOpt);
[~, ~, model] = buildDoaDopplerMfInit(model, defaultOptVar);

rowCell = cell(0, 1);
rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
  "default-final", "default-final", "ms-mf-cp-u-in-tooth", defaultOptVar, [0; 0], ...
  model, truth, toothStepHz);

truthIdx = find([methodList.label] == "ms-mf-cp-u-truth-doa-oracle", 1, 'first');
if ~isempty(truthIdx) && truthIdx <= numel(methodCaseCell)
  truthOptVar = localResolveOptVarFromCase(methodCaseCell{truthIdx});
  if numel(truthOptVar) == numel(defaultOptVar)
    rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
      "truth-doa-oracle-final", "truth-doa-oracle-final", "ms-mf-cp-u-truth-doa-oracle", ...
      truthOptVar, truthOptVar(1:2) - defaultOptVar(1:2), model, truth, toothStepHz);
  end
end

wideIdx = find([methodList.label] == "ms-mf-cp-u-wide-doa-in-tooth", 1, 'first');
if ~isempty(wideIdx) && wideIdx <= numel(methodCaseCell)
  wideOptVar = localResolveOptVarFromCase(methodCaseCell{wideIdx});
  rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
    "wide-coarse-doa-grid", "ms-mf-cp-u-wide-doa-in-tooth", wideOptVar, defaultOptVar, ...
    config.coarseDoaStepDegList, model, truth, toothStepHz);
end

singleIdx = find([methodList.label] == "ss-mf-cp-u-in-tooth", 1, 'first');
if ~isempty(singleIdx) && singleIdx <= numel(methodCaseCell)
  singleOptVar = localResolveOptVarFromCase(methodCaseCell{singleIdx});
  rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
    "single-mf-coarse-doa-grid", "ss-mf-cp-u-in-tooth", singleOptVar, defaultOptVar, ...
    config.coarseDoaStepDegList, model, truth, toothStepHz);
end

candidateProbeTable = struct2table([rowCell{:}].');
defaultMask = candidateProbeTable.candidateFamily == "default-final";
if any(defaultMask)
  defaultObj = candidateProbeTable.obj(find(defaultMask, 1, 'first'));
  candidateProbeTable.objGainFromDefault = defaultObj - candidateProbeTable.obj;
end
end

function dynOpt = localBuildDynOptForProbe(method, fixture, flowOpt, seedDoaParam)
%LOCALBUILDDYNOPTFORPROBE Build the same evaluator options used by the default in-tooth method.

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = method.freezeDoa;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.verbose = false;
end

function rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, candidateFamily, sourceMethodLabel, centerOptVar, defaultOptVar, probeDoaStepDegList, model, truth, toothStepHz)
%LOCALAPPENDPURECENTEREDPROBEGRID Add a DoA-only grid around one implementable center.

centerOptVar = reshape(double(centerOptVar), [], 1);
defaultOptVar = reshape(double(defaultOptVar), [], 1);
probeDoaStepDegList = reshape(double(probeDoaStepDegList), [], 1);
if numel(centerOptVar) ~= numel(defaultOptVar) || numel(centerOptVar) < 4
  return;
end
if any(~isfinite(centerOptVar(1:4))) || any(~isfinite(defaultOptVar(1:4)))
  return;
end
centerLabel = string(sourceMethodLabel) + "-coarse-center";
for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    doaStepLocal = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
    optVar = centerOptVar;
    optVar(1:2) = optVar(1:2) + doaStepLocal;
    doaStepFromDefault = optVar(1:2) - defaultOptVar(1:2);
    tag = sprintf('coarse-center-dlat%+.4g-dlon%+.4g', doaStepLocal(1), doaStepLocal(2));
    row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, candidateFamily, string(tag), ...
      sourceMethodLabel, optVar, doaStepFromDefault, model, truth, toothStepHz);
    row.probeCenterFamily = centerLabel;
    rowCell{end + 1, 1} = row; %#ok<AGROW>
  end
end
end

function row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, candidateFamily, candidateTag, sourceMethodLabel, optVar, doaStepDeg, model, truth, toothStepHz)
%LOCALEVALUATECANDIDATEPROBEPOINT Evaluate one candidate optVar and pack compact metrics.

optVar = reshape(double(optVar), [], 1);
[obj, ~, noiseVar, evalAux] = evaluateDoaDopplerMfObjective(model, optVar);
evalDiag = buildDoaDopplerMfEvalDiag(model, optVar, obj, noiseVar, evalAux);
truthLatlon = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', [NaN; NaN]), [], 1);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
[nonRefCoherenceFloor, nonRefRmsPhaseResidRad] = localExtractProbeNonRefMetrics(evalDiag);
fdRefErrHz = optVar(3) - truthFdRefHz;
fdRateErrHzPerSec = optVar(4) - truthFdRateHzPerSec;
toothIdx = NaN;
toothResidualHz = NaN;
if isfinite(toothStepHz) && toothStepHz > 0 && isfinite(fdRefErrHz)
  toothIdx = round(fdRefErrHz / toothStepHz);
  toothResidualHz = fdRefErrHz - toothIdx * toothStepHz;
end
row = struct();
row.iRepeat = iRepeat;
row.taskSeed = taskSeed;
row.candidateFamily = string(candidateFamily);
row.candidateTag = string(candidateTag);
row.sourceMethodLabel = string(sourceMethodLabel);
row.probeCenterFamily = "";
row.doaStep1Deg = localVectorElem(doaStepDeg, 1, NaN);
row.doaStep2Deg = localVectorElem(doaStepDeg, 2, NaN);
row.doa1Deg = optVar(1);
row.doa2Deg = optVar(2);
row.fdRefHz = optVar(3);
row.fdRateHzPerSec = optVar(4);
row.obj = obj;
row.objGainFromDefault = NaN;
row.residualNorm = localGetFieldOrDefault(evalDiag, 'residualNorm', NaN);
row.angleErrDeg = calcLatlonAngleError(optVar(1:2), truthLatlon(1:2));
row.fdRefErrHz = fdRefErrHz;
row.fdRateErrHzPerSec = fdRateErrHzPerSec;
row.toothIdx = toothIdx;
row.toothResidualHz = toothResidualHz;
row.nonRefCoherenceFloor = nonRefCoherenceFloor;
row.nonRefRmsPhaseResidRad = nonRefRmsPhaseResidRad;
row.nonRefFitRatioFloor = localGetFieldOrDefault(evalDiag, 'nonRefFitRatioFloor', NaN);
row.nonRefSupportRatioFloor = localGetFieldOrDefault(evalDiag, 'nonRefSupportRatioFloor', NaN);
end

function [nonRefCoherenceFloor, nonRefRmsPhaseResidRad] = localExtractProbeNonRefMetrics(evalDiag)
%LOCALEXTRACTPROBENONREFMETRICS Extract non-reference coherence and phase residual metrics.

nonRefCoherenceFloor = NaN;
nonRefRmsPhaseResidRad = NaN;
refSatIdxLocal = localGetFieldOrDefault(evalDiag, 'refSatIdxLocal', NaN);
coherenceSat = reshape(localGetFieldOrDefault(evalDiag, 'coherenceSat', []), [], 1);
if isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= numel(coherenceSat)
  nonRefMask = true(size(coherenceSat));
  nonRefMask(refSatIdxLocal) = false;
  nonRefCoherenceFloor = min(coherenceSat(nonRefMask), [], 'omitnan');
else
  nonRefCoherenceFloor = min(coherenceSat, [], 'omitnan');
end
blockSummary = localGetFieldOrDefault(evalDiag, 'blockSummary', []);
if istable(blockSummary) && ismember('rmsPhaseResidRad', blockSummary.Properties.VariableNames)
  phaseResid = blockSummary.rmsPhaseResidRad;
  if ismember('isRefSat', blockSummary.Properties.VariableNames)
    phaseResid = phaseResid(~blockSummary.isRefSat);
  end
  nonRefRmsPhaseResidRad = max(phaseResid, [], 'omitnan');
end
end

function optVar = localResolveOptVarFromCase(caseUse)
%LOCALRESOLVEOPTVARFROMCASE Read a final optimizer vector from one case result.

optVar = [];
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if ~isstruct(estResult)
  return;
end
optVar = reshape(localGetFieldOrDefault(estResult, 'optVarEst', []), [], 1);
if ~isempty(optVar)
  return;
end
doaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
fdRef = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRate = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
if numel(doaParam) >= 2 && isfinite(fdRef) && isfinite(fdRate)
  optVar = [doaParam(1:2); fdRef; fdRate];
end
end

function row = localBuildSummaryRow(iRepeat, taskSeed, methodLabel, methodFamily, satMode, modelClass, isKnownRate, doaSeedMode, initMode, fdRangeMode, fdRateRangeMode, doaHalfWidthMaxDeg, freezeDoa, caseUse, truth, toothStepHz, wallTimeMs)
%LOCALBUILDSUMMARYROW Convert one method result into a stable table row.

info = summarizeDoaDopplerCase(caseUse, truth);
dynSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth);
row = struct();
row.iRepeat = iRepeat;
row.taskSeed = taskSeed;
row.methodLabel = string(methodLabel);
row.methodFamily = string(methodFamily);
row.satMode = string(satMode);
row.modelClass = string(modelClass);
row.rateMode = localInferRateMode(modelClass, isKnownRate);
row.oracleLevel = localInferOracleLevel(fdRangeMode, fdRateRangeMode, doaSeedMode);
row.isKnownRate = logical(isKnownRate);
row.doaSeedMode = string(doaSeedMode);
row.initMode = string(initMode);
row.fdRangeMode = string(fdRangeMode);
row.fdRateRangeMode = string(fdRateRangeMode);
row.doaHalfWidthMaxDeg = doaHalfWidthMaxDeg;
row.freezeDoa = logical(freezeDoa);
row.solveVariant = string(localGetFieldOrDefault(dynSummary, 'solveVariant', ""));
row.isResolved = logical(localGetFieldOrDefault(info, 'isResolved', false));
row.angleErrDeg = localGetFieldOrDefault(info, 'angleErrDeg', NaN);
row.fdRefErrHz = localGetFieldOrDefault(info, 'fdRefErrHz', NaN);
row.fdRateErrHzPerSec = localGetFieldOrDefault(info, 'fdRateErrHzPerSec', NaN);
row.fdLineRmseHz = localGetFieldOrDefault(info, 'fdLineRmseHz', NaN);
row.fdSatRmseHz = localGetFieldOrDefault(info, 'fdSatRmseHz', NaN);
row.deltaFdRmseHz = localGetFieldOrDefault(info, 'deltaFdRmseHz', NaN);
row.toothIdx = localGetFieldOrDefault(dynSummary, 'toothIdx', NaN);
row.toothResidualHz = localGetFieldOrDefault(dynSummary, 'toothResidualHz', NaN);
row.finalObj = localGetFieldOrDefault(dynSummary, 'finalObj', NaN);
row.finalResidualNorm = localGetFieldOrDefault(dynSummary, 'finalResidualNorm', NaN);
row.nonRefCoherenceFloor = localGetFieldOrDefault(dynSummary, 'nonRefCoherenceFloor', NaN);
row.nonRefRmsPhaseResidRad = localGetFieldOrDefault(dynSummary, 'nonRefRmsPhaseResidRad', NaN);
row.runTimeMs = localGetFieldOrDefault(info, 'runTimeMs', NaN);
row.wallTimeMs = wallTimeMs;
end

function rateMode = localInferRateMode(modelClass, isKnownRate)
%LOCALINFERRATEMODE Return a compact rate-mode label for summary tables.

if string(modelClass) == "static"
  rateMode = "static";
elseif isKnownRate
  rateMode = "known-rate";
else
  rateMode = "unknown-rate";
end
end

function oracleLevel = localInferOracleLevel(fdRangeMode, fdRateRangeMode, doaSeedMode)
%LOCALINFERORACLELEVEL Return the oracle scope used by one replay row.

labelPart = strings(0, 1);
if string(fdRangeMode) == "truthHalfTooth"
  labelPart(end + 1, 1) = "fd";
end
if string(fdRateRangeMode) == "truthLocal"
  labelPart(end + 1, 1) = "fdRate";
end
if string(doaSeedMode) == "truth"
  labelPart(end + 1, 1) = "truthDoA";
end
if isempty(labelPart)
  oracleLevel = "none";
else
  oracleLevel = strjoin(labelPart, "+");
end
end

function methodTable = localBuildMethodTable(repeatCell)
%LOCALBUILDMETHODTABLE Stack per-repeat method summaries into one table.

summaryCell = cell(numel(repeatCell), 1);
for iRepeat = 1:numel(repeatCell)
  summaryCell{iRepeat} = repeatCell{iRepeat}.methodSummary;
end
methodTable = vertcat(summaryCell{:});
end

function aggregateTable = localBuildMethodAggregateTable(methodTable)
%LOCALBUILDMETHODAGGREGATETABLE Build compact method-level distribution metrics.

methodList = unique(methodTable.methodLabel, 'stable');
rowCell = cell(numel(methodList), 1);
for iMethod = 1:numel(methodList)
  methodNow = methodList(iMethod);
  mask = methodTable.methodLabel == methodNow;
  idxFirst = find(mask, 1, 'first');
  row = struct();
  row.methodLabel = methodNow;
  row.methodFamily = methodTable.methodFamily(idxFirst);
  row.satMode = methodTable.satMode(idxFirst);
  row.modelClass = methodTable.modelClass(idxFirst);
  row.numRepeat = nnz(mask);
  row.resolvedRate = mean(double(methodTable.isResolved(mask)), 'omitnan');
  row.toothHitRate = mean(double(methodTable.toothIdx(mask) == 0), 'omitnan');
  row.angleRmseDeg = sqrt(mean(methodTable.angleErrDeg(mask) .^ 2, 'omitnan'));
  row.angleMedianDeg = median(methodTable.angleErrDeg(mask), 'omitnan');
  row.angleP95Deg = localPercentile(methodTable.angleErrDeg(mask), 95);
  row.fdRefAbsMedianHz = median(abs(methodTable.fdRefErrHz(mask)), 'omitnan');
  row.fdRateAbsMedianHzPerSec = median(abs(methodTable.fdRateErrHzPerSec(mask)), 'omitnan');
  row.nonRefCoherenceFloorMedian = median(methodTable.nonRefCoherenceFloor(mask), 'omitnan');
  row.wallTimeMedianMs = median(methodTable.wallTimeMs(mask), 'omitnan');
  rowCell{iMethod} = row;
end
aggregateTable = struct2table([rowCell{:}].');
end

function candidateProbeTable = localBuildCandidateProbeTable(repeatCell)
%LOCALBUILDCANDIDATEPROBETABLE Stack per-repeat candidate objective probes.

probeCell = cell(0, 1);
for iRepeat = 1:numel(repeatCell)
  if isfield(repeatCell{iRepeat}, 'candidateProbeTable') && ~isempty(repeatCell{iRepeat}.candidateProbeTable)
    probeCell{end + 1, 1} = repeatCell{iRepeat}.candidateProbeTable; %#ok<AGROW>
  end
end
if isempty(probeCell)
  candidateProbeTable = table();
else
  candidateProbeTable = vertcat(probeCell{:});
end
end

function rescueEffectCaseTable = localBuildRescueEffectCaseTable(repeatCell)
%LOCALBUILDRESCUEEFFECTCASETABLE Stack per-repeat rescue bank decisions.

caseCell = cell(0, 1);
for iRepeat = 1:numel(repeatCell)
  if isfield(repeatCell{iRepeat}, 'rescueCaseTable') && ~isempty(repeatCell{iRepeat}.rescueCaseTable)
    caseCell{end + 1, 1} = repeatCell{iRepeat}.rescueCaseTable; %#ok<AGROW>
  end
end
if isempty(caseCell)
  rescueEffectCaseTable = table();
else
  rescueEffectCaseTable = vertcat(caseCell{:});
end
end

function rescueCaseTable = localBuildRescueCaseTableForRepeat(taskSeed, candidateProbeTable, config, repeatWallTimeMs)
%LOCALBUILDRESCUECASETABLEFORREPEAT Select each rescue bank without using truth in the trigger.

bankLabelList = ["disabled"; "gated-wide-only"; "gated-single-mf-only"; ...
  "gated-wide-single-bank"; "blanket-wide-single-bank"];
rowCell = cell(numel(bankLabelList), 1);
defaultRow = localBestProbeRow(candidateProbeTable, "default-final");
truthRow = localBestProbeRow(candidateProbeTable, "truth-doa-oracle-final");
wideRow = localBestProbeRow(candidateProbeTable, "wide-coarse-doa-grid");
singleRow = localBestProbeRow(candidateProbeTable, "single-mf-coarse-doa-grid");
defaultAngle = localProbeValue(defaultRow, 'angleErrDeg', NaN);
defaultCoherence = localProbeValue(defaultRow, 'nonRefCoherenceFloor', NaN);
defaultPhaseResid = localProbeValue(defaultRow, 'nonRefRmsPhaseResidRad', NaN);
defaultFdRefAbs = abs(localProbeValue(defaultRow, 'fdRefErrHz', NaN));
defaultFdRateAbs = abs(localProbeValue(defaultRow, 'fdRateErrHzPerSec', NaN));
defaultToothIdx = localProbeValue(defaultRow, 'toothIdx', NaN);
truthAngle = localProbeValue(truthRow, 'angleErrDeg', NaN);
caseRole = localClassifyCaseRole(defaultAngle, truthAngle, defaultCoherence, ...
  defaultFdRefAbs, defaultFdRateAbs, defaultToothIdx, config);
[triggered, gateReason, triggeredByCoherence, triggeredByPhaseResid] = ...
  localShouldTriggerGatedRescue(defaultCoherence, defaultPhaseResid, config);
for iBank = 1:numel(bankLabelList)
  bankLabel = bankLabelList(iBank);
  switch bankLabel
    case "disabled"
      selectedRow = defaultRow;
      rescueTriggered = false;
      selectedGateReason = "disabled";
    case "gated-wide-only"
      rescueTriggered = triggered;
      selectedGateReason = gateReason;
      selectedRow = defaultRow;
      if rescueTriggered
        selectedRow = localBestProbeRowAny(candidateProbeTable, ["default-final"; "wide-coarse-doa-grid"]);
      end
    case "gated-single-mf-only"
      rescueTriggered = triggered;
      selectedGateReason = gateReason;
      selectedRow = defaultRow;
      if rescueTriggered
        selectedRow = localBestProbeRowAny(candidateProbeTable, ["default-final"; "single-mf-coarse-doa-grid"]);
      end
    case "gated-wide-single-bank"
      rescueTriggered = triggered;
      selectedGateReason = gateReason;
      selectedRow = defaultRow;
      if rescueTriggered
        selectedRow = localBestProbeRowAny(candidateProbeTable, ["default-final"; "wide-coarse-doa-grid"; "single-mf-coarse-doa-grid"]);
      end
    case "blanket-wide-single-bank"
      rescueTriggered = true;
      selectedGateReason = "blanket-reference";
      selectedRow = localBestProbeRowAny(candidateProbeTable, ["default-final"; "wide-coarse-doa-grid"; "single-mf-coarse-doa-grid"]);
    otherwise
      selectedRow = defaultRow;
      rescueTriggered = false;
      selectedGateReason = "unknown-bank";
  end
  rowCell{iBank} = localMakeRescueCaseRow(taskSeed, caseRole, bankLabel, selectedRow, defaultRow, truthRow, ...
    rescueTriggered, selectedGateReason, triggeredByCoherence, triggeredByPhaseResid, ...
    defaultCoherence, defaultPhaseResid, repeatWallTimeMs, config);
end
rescueCaseTable = struct2table([rowCell{:}].');
end

function row = localMakeRescueCaseRow(taskSeed, caseRole, bankLabel, selectedRow, defaultRow, truthRow, rescueTriggered, gateReason, triggeredByCoherence, triggeredByPhaseResid, defaultCoherence, defaultPhaseResid, repeatWallTimeMs, config)
%LOCALMAKERESCUECASEROW Convert one selected bank candidate into the user-facing case table.

defaultAngle = localProbeValue(defaultRow, 'angleErrDeg', NaN);
selectedAngle = localProbeValue(selectedRow, 'angleErrDeg', NaN);
truthAngle = localProbeValue(truthRow, 'angleErrDeg', NaN);
angleGain = defaultAngle - selectedAngle;
row = struct();
row.taskSeed = taskSeed;
row.bankLabel = string(bankLabel);
row.caseRole = string(caseRole);
row.defaultAngleErrDeg = defaultAngle;
row.defaultNonRefCoherenceFloor = defaultCoherence;
row.defaultNonRefRmsPhaseResidRad = defaultPhaseResid;
row.defaultFdRefErrHz = localProbeValue(defaultRow, 'fdRefErrHz', NaN);
row.defaultFdRateErrHzPerSec = localProbeValue(defaultRow, 'fdRateErrHzPerSec', NaN);
row.defaultToothIdx = localProbeValue(defaultRow, 'toothIdx', NaN);
row.truthDoaOracleAngleErrDeg = truthAngle;
row.rescueTriggered = logical(rescueTriggered);
row.triggeredByCoherence = logical(triggeredByCoherence);
row.triggeredByPhaseResid = logical(triggeredByPhaseResid);
row.phaseResidAvailable = isfinite(defaultPhaseResid);
row.triggerReason = string(gateReason);
row.selectedBank = string(bankLabel);
row.selectedCandidateFamily = string(localProbeValue(selectedRow, 'candidateFamily', ""));
row.selectedCandidateTag = string(localProbeValue(selectedRow, 'candidateTag', ""));
row.selectedSourceMethodLabel = string(localProbeValue(selectedRow, 'sourceMethodLabel', ""));
row.selectedAngleErrDeg = selectedAngle;
row.selectedNonRefCoherenceFloor = localProbeValue(selectedRow, 'nonRefCoherenceFloor', NaN);
row.selectedNonRefRmsPhaseResidRad = localProbeValue(selectedRow, 'nonRefRmsPhaseResidRad', NaN);
row.selectedObj = localProbeValue(selectedRow, 'obj', NaN);
row.selectedObjGainFromDefault = localProbeValue(selectedRow, 'objGainFromDefault', NaN);
row.selectedFdRefErrHz = localProbeValue(selectedRow, 'fdRefErrHz', NaN);
row.selectedFdRateErrHzPerSec = localProbeValue(selectedRow, 'fdRateErrHzPerSec', NaN);
row.selectedToothIdx = localProbeValue(selectedRow, 'toothIdx', NaN);
row.selectedDoaStep1Deg = localProbeValue(selectedRow, 'doaStep1Deg', NaN);
row.selectedDoaStep2Deg = localProbeValue(selectedRow, 'doaStep2Deg', NaN);
row.angleGainDeg = angleGain;
row.wallTimeMs = repeatWallTimeMs;
row.gapToTruthDoaOracleDeg = selectedAngle - truthAngle;
row.isHardRescued = row.caseRole == "hard-collapse" && ...
  isfinite(angleGain) && angleGain >= config.angleGainThresholdDeg && ...
  isfinite(row.selectedNonRefCoherenceFloor) && row.selectedNonRefCoherenceFloor >= config.coherenceHealthyThreshold;
row.isDamaged = isfinite(angleGain) && angleGain <= -config.damageThresholdDeg;
end

function caseRole = localClassifyCaseRole(defaultAngle, truthAngle, defaultCoherence, defaultFdRefAbs, defaultFdRateAbs, defaultToothIdx, config)
%LOCALCLASSIFYCASEROLE Assign an offline evaluation label for one repeat.

sameTooth = isfinite(defaultToothIdx) && defaultToothIdx == 0;
fdHealthy = sameTooth && isfinite(defaultFdRefAbs) && defaultFdRefAbs <= config.fdRefHealthyAbsHz && ...
  isfinite(defaultFdRateAbs) && defaultFdRateAbs <= config.fdRateHealthyAbsHzPerSec;
truthGap = defaultAngle - truthAngle;
coherenceCollapsed = isfinite(defaultCoherence) && defaultCoherence < config.coherenceHealthyThreshold;
if ~sameTooth
  caseRole = "wrong-tooth-negative";
elseif ~fdHealthy
  caseRole = "fd-not-healthy-negative";
elseif coherenceCollapsed && isfinite(truthGap) && truthGap >= config.truthDoaGapThresholdDeg
  caseRole = "hard-collapse";
else
  caseRole = "easy-negative";
end
end

function [triggerRescue, gateReason, triggeredByCoherence, triggeredByPhaseResid] = localShouldTriggerGatedRescue(defaultCoherenceFloor, defaultRmsPhaseResidRad, config)
%LOCALSHOULDTRIGGERGATEDRESCUE Gate the implementable bank using no-truth diagnostics.

cohThreshold = config.gatedRescueCoherenceThreshold;
phaseThreshold = config.gatedRescuePhaseResidThresholdRad;
phaseGateEnabled = isfinite(phaseThreshold) && phaseThreshold > 0;
phaseResidAvailable = phaseGateEnabled && isfinite(defaultRmsPhaseResidRad);
triggeredByCoherence = isfinite(defaultCoherenceFloor) && defaultCoherenceFloor < cohThreshold;
triggeredByPhaseResid = phaseResidAvailable && defaultRmsPhaseResidRad >= phaseThreshold;
triggerRescue = triggeredByCoherence || triggeredByPhaseResid;

if triggeredByCoherence && triggeredByPhaseResid
  gateReason = "coherence-collapse-and-phase-residual-large";
elseif triggeredByCoherence && ~phaseResidAvailable
  gateReason = "coherence-collapse-phase-unavailable";
elseif triggeredByCoherence
  gateReason = "coherence-collapse";
elseif triggeredByPhaseResid
  gateReason = "phase-residual-large";
elseif ~isfinite(defaultCoherenceFloor)
  gateReason = "coherence-unavailable";
elseif ~phaseResidAvailable
  gateReason = "coherence-not-collapsed-phase-unavailable";
else
  gateReason = "coherence-not-collapsed";
end
end

function rescueEffectCaseTable = localSelectGatedEffectCaseTable(rescueBankDecisionTable)
%LOCALSELECTGATEDEFFECTCASETABLE Keep one gated wide+single-MF result row per seed.

if isempty(rescueBankDecisionTable)
  rescueEffectCaseTable = table();
  return;
end
mask = rescueBankDecisionTable.bankLabel == "gated-wide-single-bank";
rescueEffectCaseTable = rescueBankDecisionTable(mask, :);
end

function triggerReasonTable = localBuildTriggerReasonTable(rescueEffectCaseTable)
%LOCALBUILDTRIGGERREASONTABLE Extract one no-truth trigger row per seed for the gated bank.

if isempty(rescueEffectCaseTable)
  triggerReasonTable = table();
  return;
end
mask = rescueEffectCaseTable.bankLabel == "gated-wide-single-bank";
base = rescueEffectCaseTable(mask, :);
keepFields = {'taskSeed', 'defaultNonRefCoherenceFloor', 'defaultNonRefRmsPhaseResidRad', ...
  'triggeredByCoherence', 'triggeredByPhaseResid', 'phaseResidAvailable', ...
  'rescueTriggered', 'triggerReason'};
triggerReasonTable = base(:, keepFields);
end

function rescueEffectAggregateTable = localBuildRescueEffectAggregateTable(rescueEffectCaseTable)
%LOCALBUILDRESCUEEFFECTAGGREGATETABLE Build rescue and damage rates for each bank.

if isempty(rescueEffectCaseTable)
  rescueEffectAggregateTable = table();
  return;
end
bankList = unique(rescueEffectCaseTable.bankLabel, 'stable');
rowCell = cell(numel(bankList), 1);
for iBank = 1:numel(bankList)
  bankNow = bankList(iBank);
  mask = rescueEffectCaseTable.bankLabel == bankNow;
  hardMask = mask & rescueEffectCaseTable.caseRole == "hard-collapse";
  easyMask = mask & rescueEffectCaseTable.caseRole == "easy-negative";
  fdMask = mask & rescueEffectCaseTable.caseRole == "fd-not-healthy-negative";
  row = struct();
  row.bankLabel = bankNow;
  row.numRepeat = nnz(mask);
  row.numHardCase = nnz(hardMask);
  row.numEasyNegativeCase = nnz(easyMask);
  row.numFdNegativeCase = nnz(fdMask);
  row.triggerRate = localMeanLogical(rescueEffectCaseTable.rescueTriggered(mask));
  row.hardTriggerRate = localMeanLogical(rescueEffectCaseTable.rescueTriggered(hardMask));
  row.easyTriggerRate = localMeanLogical(rescueEffectCaseTable.rescueTriggered(easyMask));
  row.fdNegativeTriggerRate = localMeanLogical(rescueEffectCaseTable.rescueTriggered(fdMask));
  row.hardRescueRate = localMeanLogical(rescueEffectCaseTable.isHardRescued(hardMask));
  row.hardMedianAngleGainDeg = median(rescueEffectCaseTable.angleGainDeg(hardMask), 'omitnan');
  row.hardP95SelectedAngleDeg = localPercentile(rescueEffectCaseTable.selectedAngleErrDeg(hardMask), 95);
  row.easyDamageRate = localMeanLogical(rescueEffectCaseTable.isDamaged(easyMask));
  row.fdNegativeDamageRate = localMeanLogical(rescueEffectCaseTable.isDamaged(fdMask));
  row.overallMedianAngleDeg = median(rescueEffectCaseTable.selectedAngleErrDeg(mask), 'omitnan');
  row.overallP95AngleDeg = localPercentile(rescueEffectCaseTable.selectedAngleErrDeg(mask), 95);
  row.overallMaxAngleDeg = max(rescueEffectCaseTable.selectedAngleErrDeg(mask), [], 'omitnan');
  row.medianWallTimeMs = median(rescueEffectCaseTable.wallTimeMs(mask), 'omitnan');
  rowCell{iBank} = row;
end
rescueEffectAggregateTable = struct2table([rowCell{:}].');
end

function verdictTable = localBuildRescueEffectVerdictTable(aggregateTable)
%LOCALBUILDRESCUEEFFECTVERDICTTABLE Summarize whether each bank satisfies replay promotion checks.

if isempty(aggregateTable)
  verdictTable = table();
  return;
end
disabledRow = localSelectAggregateRow(aggregateTable, "disabled");
blanketRow = localSelectAggregateRow(aggregateTable, "blanket-wide-single-bank");
bankMask = aggregateTable.bankLabel ~= "disabled";
bankList = aggregateTable.bankLabel(bankMask);
rowCell = cell(numel(bankList), 1);
for iBank = 1:numel(bankList)
  bankNow = bankList(iBank);
  bankRow = localSelectAggregateRow(aggregateTable, bankNow);
  numHardCase = localTableScalar(bankRow, 'numHardCase', NaN);
  numFdNegativeCase = localTableScalar(bankRow, 'numFdNegativeCase', NaN);
  hardRescueRate = localTableScalar(bankRow, 'hardRescueRate', NaN);
  easyDamageRate = localTableScalar(bankRow, 'easyDamageRate', NaN);
  fdDamageRate = localTableScalar(bankRow, 'fdNegativeDamageRate', NaN);
  disabledP95 = localTableScalar(disabledRow, 'overallP95AngleDeg', NaN);
  disabledMax = localTableScalar(disabledRow, 'overallMaxAngleDeg', NaN);
  bankP95 = localTableScalar(bankRow, 'overallP95AngleDeg', NaN);
  bankMax = localTableScalar(bankRow, 'overallMaxAngleDeg', NaN);
  blanketFdDamage = localTableScalar(blanketRow, 'fdNegativeDamageRate', NaN);
  p95Improvement = disabledP95 - bankP95;
  maxImprovement = disabledMax - bankMax;
  hardRescuePass = isfinite(hardRescueRate) && hardRescueRate >= 0.8;
  easyDamagePass = ~isfinite(easyDamageRate) || easyDamageRate <= 0;
  fdDamagePass = ~isfinite(fdDamageRate) || fdDamageRate <= 0.05 || ...
    (isfinite(blanketFdDamage) && fdDamageRate < blanketFdDamage);
  p95Pass = ~isfinite(p95Improvement) || p95Improvement >= 0;
  maxPass = ~isfinite(maxImprovement) || maxImprovement >= 0;
  enoughHardCase = isfinite(numHardCase) && numHardCase >= 3;
  recommendedPass = enoughHardCase && hardRescuePass && easyDamagePass && fdDamagePass && p95Pass && maxPass;
  row = struct();
  row.bankLabel = bankNow;
  row.numHardCase = numHardCase;
  row.numFdNegativeCase = numFdNegativeCase;
  row.hardRescueRate = hardRescueRate;
  row.easyDamageRate = easyDamageRate;
  row.fdNegativeDamageRate = fdDamageRate;
  row.p95ImprovementVsDisabledDeg = p95Improvement;
  row.maxImprovementVsDisabledDeg = maxImprovement;
  row.hardRescuePass = hardRescuePass;
  row.easyDamagePass = easyDamagePass;
  row.fdDamagePass = fdDamagePass;
  row.p95NotWorseThanDisabled = p95Pass;
  row.maxNotWorseThanDisabled = maxPass;
  row.enoughHardCase = enoughHardCase;
  row.recommendedPass = recommendedPass;
  row.recommendation = localBuildVerdictRecommendation(bankNow, recommendedPass, enoughHardCase);
  rowCell{iBank} = row;
end
verdictTable = struct2table([rowCell{:}].');
end

function row = localSelectAggregateRow(aggregateTable, bankLabel)
%LOCALSELECTAGGREGATEROW Return one aggregate row for a bank label.

row = aggregateTable([], :);
if isempty(aggregateTable) || ~ismember('bankLabel', aggregateTable.Properties.VariableNames)
  return;
end
mask = aggregateTable.bankLabel == string(bankLabel);
if any(mask)
  idx = find(mask, 1, 'first');
  row = aggregateTable(idx, :);
end
end

function value = localTableScalar(dataTable, fieldName, defaultValue)
%LOCALTABLESCALAR Read one scalar value from a table row.

value = defaultValue;
if isempty(dataTable) || height(dataTable) == 0 || ~ismember(fieldName, dataTable.Properties.VariableNames)
  return;
end
rawValue = dataTable.(fieldName);
if ~isempty(rawValue)
  value = rawValue(1);
end
end

function recommendation = localBuildVerdictRecommendation(bankLabel, recommendedPass, enoughHardCase)
%LOCALBUILDVERDICTRECOMMENDATION Return a compact replay-only next-step label.

if ~enoughHardCase
  recommendation = "needs-more-hard-cases";
elseif recommendedPass && string(bankLabel) == "gated-wide-single-bank"
  recommendation = "ready-for-flow-like-replay";
elseif recommendedPass
  recommendation = "candidate-but-not-target";
else
  recommendation = "hold";
end
end

function bestRow = localBestProbeRow(seedProbe, candidateFamily)
%LOCALBESTPROBEROW Return the lowest-objective row for one candidate family.

bestRow = seedProbe([], :);
if isempty(seedProbe) || ~ismember('candidateFamily', seedProbe.Properties.VariableNames)
  return;
end
mask = seedProbe.candidateFamily == string(candidateFamily) & isfinite(seedProbe.obj);
if ~any(mask)
  return;
end
idx = find(mask);
[~, idxRel] = min(seedProbe.obj(idx));
bestRow = seedProbe(idx(idxRel), :);
end

function bestRow = localBestProbeRowAny(seedProbe, candidateFamilyList)
%LOCALBESTPROBEROWANY Return the lowest-objective row across several candidate families.

bestRow = seedProbe([], :);
if isempty(seedProbe) || ~ismember('candidateFamily', seedProbe.Properties.VariableNames)
  return;
end
candidateFamilyList = string(candidateFamilyList(:));
mask = ismember(seedProbe.candidateFamily, candidateFamilyList) & isfinite(seedProbe.obj);
if ~any(mask)
  return;
end
idx = find(mask);
[~, idxRel] = min(seedProbe.obj(idx));
bestRow = seedProbe(idx(idxRel), :);
end

function value = localProbeValue(row, fieldName, defaultValue)
%LOCALPROBEVALUE Read a scalar field from one probe table row.

value = defaultValue;
if isempty(row) || height(row) == 0 || ~ismember(fieldName, row.Properties.VariableNames)
  return;
end
rawValue = row.(fieldName);
if ~isempty(rawValue)
  value = rawValue(1);
end
end

function value = localMeanLogical(maskValue)
%LOCALMEANLOGICAL Return the mean of a logical vector with empty fallback.

if isempty(maskValue)
  value = NaN;
else
  value = mean(double(maskValue), 'omitnan');
end
end

function rangeTable = localBuildRangeTable(repeatCell)
%LOCALBUILDRANGETABLE Summarize the oracle fd / fdRate boxes used per repeat.

numRepeat = numel(repeatCell);
taskSeed = nan(numRepeat, 1);
toothStepHz = nan(numRepeat, 1);
truthFdRefHz = nan(numRepeat, 1);
truthFdRateHzPerSec = nan(numRepeat, 1);
fdRangeLowerHz = nan(numRepeat, 1);
fdRangeUpperHz = nan(numRepeat, 1);
fdRateRangeLowerHzPerSec = nan(numRepeat, 1);
fdRateRangeUpperHzPerSec = nan(numRepeat, 1);
for iRepeat = 1:numRepeat
  item = repeatCell{iRepeat};
  taskSeed(iRepeat) = item.taskSeed;
  toothStepHz(iRepeat) = item.toothStepHz;
  truthFdRefHz(iRepeat) = item.truthFdRefHz;
  truthFdRateHzPerSec(iRepeat) = item.truthFdRateHzPerSec;
  fdRangeLowerHz(iRepeat) = item.fdRangeOracle(1);
  fdRangeUpperHz(iRepeat) = item.fdRangeOracle(2);
  fdRateRangeLowerHzPerSec(iRepeat) = item.fdRateRangeOracle(1);
  fdRateRangeUpperHzPerSec(iRepeat) = item.fdRateRangeOracle(2);
end
rangeTable = table(taskSeed, toothStepHz, truthFdRefHz, truthFdRateHzPerSec, ...
  fdRangeLowerHz, fdRangeUpperHz, fdRateRangeLowerHzPerSec, fdRateRangeUpperHzPerSec);
end

function timingTable = localBuildTimingTable(repeatCell)
%LOCALBUILDTIMINGTABLE Stack per-repeat wall-clock timing summaries.

numRepeat = numel(repeatCell);
taskSeed = nan(numRepeat, 1);
buildDataMs = nan(numRepeat, 1);
staticBundleMs = nan(numRepeat, 1);
dynamicMethodTotalMs = nan(numRepeat, 1);
dynamicMethodMaxMs = nan(numRepeat, 1);
repeatTotalMs = nan(numRepeat, 1);
for iRepeat = 1:numRepeat
  item = repeatCell{iRepeat};
  taskSeed(iRepeat) = item.taskSeed;
  timing = localGetFieldOrDefault(item, 'timing', struct());
  buildDataMs(iRepeat) = localGetFieldOrDefault(timing, 'buildDataMs', NaN);
  staticBundleMs(iRepeat) = localGetFieldOrDefault(timing, 'staticBundleMs', NaN);
  dynamicMethodTotalMs(iRepeat) = localGetFieldOrDefault(timing, 'dynamicMethodTotalMs', NaN);
  dynamicMethodMaxMs(iRepeat) = localGetFieldOrDefault(timing, 'dynamicMethodMaxMs', NaN);
  repeatTotalMs(iRepeat) = localGetFieldOrDefault(timing, 'repeatTotalMs', NaN);
end
timingTable = table(taskSeed, buildDataMs, staticBundleMs, dynamicMethodTotalMs, dynamicMethodMaxMs, repeatTotalMs);
end

function timingAggregateTable = localBuildTimingAggregateTable(timingTable)
%LOCALBUILDTIMINGAGGREGATETABLE Build compact timing summary across repeats.

if isempty(timingTable)
  timingAggregateTable = table();
  return;
end
stageLabel = ["build-repeat-data"; "static-transition-bundle"; "dynamic-methods-total"; ...
  "slowest-dynamic-method"; "repeat-total"];
valueCell = {timingTable.buildDataMs; timingTable.staticBundleMs; timingTable.dynamicMethodTotalMs; ...
  timingTable.dynamicMethodMaxMs; timingTable.repeatTotalMs};
medianMs = nan(numel(stageLabel), 1);
p95Ms = nan(numel(stageLabel), 1);
totalSec = nan(numel(stageLabel), 1);
for iStage = 1:numel(stageLabel)
  value = valueCell{iStage};
  medianMs(iStage) = median(value, 'omitnan');
  p95Ms(iStage) = localPercentile(value, 95);
  totalSec(iStage) = sum(value, 'omitnan') / 1000;
end
timingAggregateTable = table(stageLabel, medianMs, p95Ms, totalSec);
end

function contextSummary = localBuildContextSummary(context, methodList)
%LOCALBUILDCONTEXTSUMMARY Keep only lightweight context metadata.

contextSummary = struct();
contextSummary.frameIntvlSec = context.frameIntvlSec;
contextSummary.periodicOffsetIdx = reshape(context.periodicOffsetIdx, 1, []);
contextSummary.masterOffsetIdx = reshape(context.masterOffsetIdx, 1, []);
contextSummary.selectedSatIdxGlobal = reshape(context.selectedSatIdxGlobal, 1, []);
contextSummary.refSatIdxGlobal = context.refSatIdxGlobal;
contextSummary.otherSatIdxGlobal = context.otherSatIdxGlobal;
contextSummary.contextBaseSeed = context.baseSeed;
contextSummary.tleFileName = string(context.tleFileName);
contextSummary.methodLabelList = reshape([methodList.label], [], 1);
contextSummary.subsetInitialization = "disabled: in-tooth rescue replay uses static/single-MF centers only";
end

function methodListOut = localStripMethodList(methodList)
%LOCALSTRIPMETHODLIST Store method descriptors without heavy fields.

methodListOut = methodList;
end

function plotData = localBuildPlotData(rescueBankDecisionTable, rescueEffectCaseTable, triggerReasonTable)
%LOCALBUILDPLOTDATA Pack lightweight plotting inputs.

plotData = struct();
plotData.rescueBankDecisionTable = rescueBankDecisionTable;
plotData.rescueEffectCaseTable = rescueEffectCaseTable;
plotData.triggerReasonTable = triggerReasonTable;
end

function localPlotReplay(plotData)
%LOCALPLOTREPLAY Draw per-seed rescue changes with line plots.

if isempty(plotData) || ~isstruct(plotData)
  return;
end
bankTable = localGetFieldOrDefault(plotData, 'rescueBankDecisionTable', table());
if isempty(bankTable)
  bankTable = localGetFieldOrDefault(plotData, 'rescueEffectCaseTable', table());
end
if isempty(bankTable)
  return;
end
gatedTable = bankTable(bankTable.bankLabel == "gated-wide-single-bank", :);
if isempty(gatedTable)
  return;
end
blanketTable = bankTable(bankTable.bankLabel == "blanket-wide-single-bank", :);
x = double(gatedTable.taskSeed);

figure('Name', 'In-tooth gated rescue seed trace');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
plot(x, max(gatedTable.defaultAngleErrDeg, eps), '-o');
hold on;
plot(x, max(gatedTable.selectedAngleErrDeg, eps), '-x');
if ~isempty(blanketTable)
  plot(double(blanketTable.taskSeed), max(blanketTable.selectedAngleErrDeg, eps), '-s');
end
plot(x, max(gatedTable.truthDoaOracleAngleErrDeg, eps), ':d');
set(gca, 'YScale', 'log');
ylabel('angle error (deg)');
title('Per-seed angle error before and after rescue');
grid on;
legend(localPlotLegend({'default MS-MF', 'gated selected', 'blanket reference', 'truth-DoA oracle'}, ~isempty(blanketTable)), ...
  'Location', 'best', 'Interpreter', 'none');

nexttile;
plot(x, gatedTable.angleGainDeg, '-o');
hold on;
if ~isempty(blanketTable)
  plot(double(blanketTable.taskSeed), blanketTable.angleGainDeg, '-s');
end
yline(0, '--');
ylabel('angle gain (deg)');
title('Per-seed rescue gain');
grid on;
legend(localPlotLegend({'gated gain', 'blanket gain'}, ~isempty(blanketTable)), ...
  'Location', 'best', 'Interpreter', 'none');

nexttile;
plot(x, gatedTable.defaultNonRefCoherenceFloor, '-o');
hold on;
plot(x, gatedTable.selectedNonRefCoherenceFloor, '-x');
if ~isempty(blanketTable)
  plot(double(blanketTable.taskSeed), blanketTable.selectedNonRefCoherenceFloor, '-s');
end
yline(0.20, '--', 'gate');
yline(0.95, ':', 'healthy');
xlabel('task seed');
ylabel('non-ref coherence floor');
title('Per-seed non-ref coherence before and after rescue');
grid on;
legend(localPlotLegend({'default MS-MF', 'gated selected', 'blanket reference'}, ~isempty(blanketTable)), ...
  'Location', 'best', 'Interpreter', 'none');

figure('Name', 'In-tooth gated rescue trigger trace');
plot(x, double(gatedTable.rescueTriggered), '-o');
hold on;
plot(x, double(gatedTable.caseRole == "hard-collapse"), '-x');
plot(x, double(gatedTable.caseRole == "fd-not-healthy-negative"), '-s');
plot(x, double(gatedTable.isDamaged), '-d');
ylim([-0.1, 1.1]);
xlabel('task seed');
ylabel('indicator');
title('No-truth trigger and offline case-role trace');
grid on;
legend({'rescue triggered', 'hard collapse', 'fd-negative', 'damaged'}, ...
  'Location', 'best', 'Interpreter', 'none');
end

function legendText = localPlotLegend(baseText, includeOptional)
%LOCALPLOTLEGEND Drop optional blanket labels when the table is unavailable.

legendText = string(baseText(:));
if ~includeOptional
  legendText = legendText(~contains(legendText, "blanket"));
end
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve the fd tooth spacing from fixture timing.

toothStepHz = NaN;
if isfield(fixture, 'frameIntvlSec') && isfinite(fixture.frameIntvlSec) && fixture.frameIntvlSec > 0
  toothStepHz = 1 / fixture.frameIntvlSec;
  return;
end
if isfield(fixture, 'sceneSeq') && isfield(fixture.sceneSeq, 'timeOffsetSec')
  dt = diff(reshape(fixture.sceneSeq.timeOffsetSec, 1, []));
  dt = dt(isfinite(dt) & dt > 0);
  if ~isempty(dt)
    toothStepHz = 1 / median(dt);
  end
end
if ~(isfinite(toothStepHz) && toothStepHz > 0)
  error('replayMfInToothGatedRescueEffectiveness:MissingToothStep', ...
    'Cannot resolve tooth step from fixture timing.');
end
end

function value = localResolveTruthScalar(truth, fieldNameList)
%LOCALRESOLVETRUTHSCALAR Resolve a scalar truth field using ordered fallbacks.

value = NaN;
for iField = 1:numel(fieldNameList)
  fieldName = fieldNameList{iField};
  rawValue = localGetFieldOrDefault(truth, fieldName, NaN);
  if isscalar(rawValue) && isfinite(rawValue)
    value = rawValue;
    return;
  end
end
end

function maxVal = localMaxAbs(value)
%LOCALMAXABS Return max abs for a possibly empty numeric vector.

if isempty(value)
  maxVal = NaN;
else
  maxVal = max(abs(value(:)));
end
end

function value = localVectorElem(vec, idx, defaultValue)
%LOCALVECTORELEM Read one vector element with fallback.

value = defaultValue;
if isempty(vec)
  return;
end
vec = reshape(vec, [], 1);
if idx >= 1 && idx <= numel(vec) && isfinite(vec(idx))
  value = vec(idx);
end
end

function percentileValue = localPercentile(valueVec, percentile)
%LOCALPERCENTILE Compute one percentile without requiring Statistics Toolbox.

valueVec = reshape(valueVec, [], 1);
valueVec = valueVec(isfinite(valueVec));
if isempty(valueVec)
  percentileValue = NaN;
  return;
end
valueVec = sort(valueVec);
if numel(valueVec) == 1
  percentileValue = valueVec;
  return;
end
rankPos = 1 + (numel(valueVec) - 1) * percentile / 100;
lo = floor(rankPos);
hi = ceil(rankPos);
if lo == hi
  percentileValue = valueVec(lo);
else
  weightHi = rankPos - lo;
  percentileValue = (1 - weightHi) * valueVec(lo) + weightHi * valueVec(hi);
end
end

function localDispTablePreview(dataTable, edgeCount)
%LOCALDISPTABLEPREVIEW Display short first/last preview for long tables.

if nargin < 2 || isempty(edgeCount)
  edgeCount = 4;
end
if isempty(dataTable) || height(dataTable) <= 2 * edgeCount
  disp(dataTable);
  return;
end
fprintf('  showing first %d and last %d of %d rows\n', edgeCount, edgeCount, height(dataTable));
disp(dataTable(1:edgeCount, :));
fprintf('  ...\n');
disp(dataTable((height(dataTable) - edgeCount + 1):height(dataTable), :));
end

function tracker = localCreateProgressTracker(titleText, totalCount, useParfor)
%LOCALCREATEPROGRESSTRACKER Create a client-side progressbar tracker.

tracker = struct('active', false, 'queue', [], 'totalCount', totalCount);
if totalCount <= 0
  return;
end
fprintf('%s\n', titleText);
if exist('progressbar', 'file') ~= 2
  fprintf('[%s] progressbar.m not found; continuing without progress bar.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  tracker.active = true;
  if useParfor
    tracker.queue = parallel.pool.DataQueue;
    afterEach(tracker.queue, @(~) progressbar('advance'));
  end
catch ME
  tracker.active = false;
  tracker.queue = [];
  fprintf('[%s] Progress bar disabled: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), ME.message);
end
end

function localAdvanceProgressTracker(tracker)
%LOCALADVANCEPROGRESSTRACKER Advance the progressbar in a serial loop.

if tracker.active && isempty(tracker.queue)
  progressbar('advance');
end
end

function localAdvanceProgressByStep(tracker, step)
%LOCALADVANCEPROGRESSBYSTEP Advance the progressbar by a checkpoint callback step.

if nargin < 2 || isempty(step)
  step = 1;
end
if tracker.active
  for iStep = 1:step
    progressbar('advance');
  end
end
end

function localCloseProgressTracker(tracker)
%LOCALCLOSEPROGRESSTRACKER Close the progressbar if it is active.

if tracker.active
  pause(0.05);
  progressbar('end');
end
end

function textOut = localModeText(useParfor)
%LOCALMODETEXT Return a short run-mode tag.

if useParfor
  textOut = 'parfor';
else
  textOut = 'serial';
end
end

function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether an outer parfor repeat loop can be used.

tf = false;
try
  tf = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or the supplied default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
