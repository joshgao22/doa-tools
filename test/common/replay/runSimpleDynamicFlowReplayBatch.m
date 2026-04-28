function [repeatTable, repeatCell, context, flowOpt, runState] = runSimpleDynamicFlowReplayBatch(batchOpt)
%RUNSIMPLEDYNAMICFLOWREPLAYBATCH Run a small fixed-seed simple-flow batch.
% This helper is shared by flow-like replay and scan entries. It keeps repeat
% construction, static seed construction, and simple-flow execution identical
% across mechanism scripts without turning them into regression.

arguments
  batchOpt (1,1) struct
end

runState = localEmptyRunState();
batchTimer = tic;

resolveTimer = tic;
localLog('Resolve replay batch options.');
batchOpt = localResolveBatchOpt(batchOpt);
localLog(sprintf('Replay batch options resolved in %.2f s.', toc(resolveTimer)));

contextTimer = tic;
localLog('Build shared dynamic context: scene, waveform, and fixture defaults.');
context = buildDynamicDualSatEciContext(batchOpt.contextOpt);
localLog(sprintf('Shared dynamic context built in %.2f s.', toc(contextTimer)));
flowOpt = batchOpt.flowOpt;

numRepeat = numel(batchOpt.seedList);
repeatCell = cell(numRepeat, 1);
parallelAvailable = localCanUseParfor();
useParfor = parallelAvailable && numRepeat > 1;
if numRepeat > 1 && ~parallelAvailable
  localLog('Parallel toolbox is unavailable; falling back to serial repeat loop.');
end
if numRepeat <= 1
  useParfor = false;
end

modeText = 'serial';
if useParfor
  modeText = 'parfor';
end
localLog(sprintf('Prepare repeat batch: numRepeat=%d, mode=%s.', numRepeat, modeText));
repeatOpt = localBuildRepeatOpt(batchOpt);
seedList = batchOpt.seedList;
repeatTimer = tic;

if localCheckpointEnabled(batchOpt)
  checkpointOpt = localResolveCheckpointOpt(batchOpt.checkpointOpt, batchOpt, useParfor);
  numDoneTask = localCountCheckpointTaskFile(fullfile(checkpointOpt.runDir, 'task'), numRepeat);
  numTodoTask = numRepeat - numDoneTask;
  trackerTitle = sprintf('%s (%s, resume %d/%d)', char(batchOpt.progressTitle), modeText, numDoneTask, numRepeat);
  tracker = localCreateProgressTracker(trackerTitle, numTodoTask, false);
  checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, tracker);
  taskGrid = localBuildRepeatTaskGrid(seedList);
  sharedData = struct('context', context, 'repeatOpt', repeatOpt, 'flowOpt', flowOpt);
  if checkpointRunnerOpt.useParfor && numTodoTask > 0
    localLog('Enter checkpointed parfor repeat loop. MATLAB may spend time starting or attaching the parallel pool before worker iterations begin.');
  elseif numTodoTask > 0
    localLog('Enter checkpointed serial repeat loop.');
  else
    localLog('All checkpoint repeat tasks are already complete; loading cached repeat results.');
  end
  try
    runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, @localCheckpointTaskRunner, checkpointRunnerOpt);
  catch ME
    localCloseProgressTracker(tracker);
    rethrow(ME);
  end
  localCloseProgressTracker(tracker);
  if ~runState.isComplete
    error('runSimpleDynamicFlowReplayBatch:CheckpointIncomplete', ...
      'Checkpoint replay returned an incomplete repeat batch.');
  end
  repeatCell = runState.resultCell;
else
  trackerTitle = sprintf('%s (%s)', char(batchOpt.progressTitle), modeText);
  tracker = localCreateProgressTracker(trackerTitle, numRepeat, useParfor);
  if useParfor
    localLog('Enter parfor repeat loop. MATLAB may spend time starting or attaching the parallel pool before worker iterations begin.');
    progressQueue = tracker.queue;
    parfor iRepeat = 1:numRepeat
      taskSeed = seedList(iRepeat);
      repeatCell{iRepeat} = localRunOneRepeat(iRepeat, taskSeed, context, repeatOpt, flowOpt);
      if ~isempty(progressQueue)
        send(progressQueue, iRepeat);
      end
    end
  else
    localLog('Enter serial repeat loop.');
    for iRepeat = 1:numRepeat
      taskSeed = seedList(iRepeat);
      repeatCell{iRepeat} = localRunOneRepeat(iRepeat, taskSeed, context, repeatOpt, flowOpt);
      localAdvanceProgressTracker(tracker);
    end
  end
  localCloseProgressTracker(tracker);
end
repeatLoopSec = toc(repeatTimer);
localLog(sprintf('Repeat batch finished in %.2f s.', repeatLoopSec));

repeatTable = localBuildRepeatTable(repeatCell);
localLog(sprintf('Replay batch complete in %.2f s.', toc(batchTimer)));
end


function batchOpt = localResolveBatchOpt(batchOpt)
% Fill missing batch options with replay-safe defaults.
batchOpt = localFillField(batchOpt, 'replayName', "mfReplay");
batchOpt = localFillField(batchOpt, 'snrDb', 10);
batchOpt = localFillField(batchOpt, 'baseSeed', 253);
batchOpt = localFillField(batchOpt, 'numRepeat', 8);
batchOpt = localFillField(batchOpt, 'seedList', batchOpt.baseSeed + (0:(batchOpt.numRepeat - 1)));
batchOpt = localFillField(batchOpt, 'optVerbose', false);
batchOpt = localFillField(batchOpt, 'runPeriodicWide', false);
batchOpt = localFillField(batchOpt, 'progressTitle', sprintf('%s repeat batch', char(batchOpt.replayName)));
batchOpt = localFillField(batchOpt, 'contextOpt', struct());
batchOpt = localFillField(batchOpt, 'checkpointOpt', struct('enable', false));

parallelOpt = localGetFieldOrDefault(batchOpt, 'parallelOpt', struct());
parallelOpt = localMergeStruct(struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false), parallelOpt);

contextOpt = batchOpt.contextOpt;
contextOpt = localMergeStruct(struct( ...
  'baseSeed', batchOpt.baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt), contextOpt);
batchOpt.contextOpt = contextOpt;
batchOpt.parallelOpt = parallelOpt;

if ~isfield(batchOpt, 'flowOpt') || isempty(batchOpt.flowOpt)
  batchOpt.flowOpt = buildSimpleDynamicFlowOpt(struct( ...
    'parallelOpt', parallelOpt, ...
    'periodicRefineFdHalfWidthHz', 50, ...
    'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
    'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
    'periodicRefineFreezeDoa', true, ...
    'periodicRefineDoaSeedMode', "dualWhenMulti", ...
    'periodicPolishEnableWhenMulti', true, ...
    'periodicPolishDoaHalfWidthDeg', [0.002; 0.002]));
end

batchOpt.seedList = reshape(double(batchOpt.seedList), [], 1);
end
function repeatOpt = localBuildRepeatOpt(batchOpt)
% Keep only the repeat-level fields that workers need.
repeatOpt = struct();
repeatOpt.replayName = batchOpt.replayName;
repeatOpt.snrDb = batchOpt.snrDb;
repeatOpt.optVerbose = batchOpt.optVerbose;
repeatOpt.runPeriodicWide = logical(localGetFieldOrDefault(batchOpt, 'runPeriodicWide', false));
end


function tf = localCheckpointEnabled(batchOpt)
% Return true when the caller requests per-repeat checkpoint files.
tf = false;
if isfield(batchOpt, 'checkpointOpt') && isstruct(batchOpt.checkpointOpt) && ...
    isfield(batchOpt.checkpointOpt, 'enable') && ~isempty(batchOpt.checkpointOpt.enable)
  tf = logical(batchOpt.checkpointOpt.enable);
end
end


function checkpointOpt = localResolveCheckpointOpt(rawOpt, batchOpt, useParfor)
% Resolve replay checkpoint options before calling the generic checkpoint runner.
checkpointOpt = rawOpt;
checkpointOpt = localFillField(checkpointOpt, 'runName', string(batchOpt.replayName));
checkpointOpt = localFillField(checkpointOpt, 'runKey', localBuildCheckpointRunKey(batchOpt));
checkpointOpt = localFillField(checkpointOpt, 'outputRoot', fullfile(localGetRepoRoot(), 'tmp'));
checkpointOpt = localFillField(checkpointOpt, 'resume', true);
checkpointOpt = localFillField(checkpointOpt, 'useParfor', useParfor);
checkpointOpt = localFillField(checkpointOpt, 'meta', localBuildCheckpointMeta(batchOpt));
checkpointOpt = localFillField(checkpointOpt, 'cleanupOnSuccess', false);
checkpointOpt = localFillField(checkpointOpt, 'cleanupOpt', struct());
checkpointOpt.runDir = fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey);
end


function checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, tracker)
% Keep only options accepted by runPerfTaskGridWithCheckpoint.
checkpointRunnerOpt = struct();
checkpointRunnerOpt.runName = checkpointOpt.runName;
checkpointRunnerOpt.runKey = checkpointOpt.runKey;
checkpointRunnerOpt.outputRoot = checkpointOpt.outputRoot;
checkpointRunnerOpt.useParfor = checkpointOpt.useParfor;
checkpointRunnerOpt.resume = checkpointOpt.resume;
checkpointRunnerOpt.meta = checkpointOpt.meta;
checkpointRunnerOpt.progressFcn = @(step) localAdvanceProgressByStep(tracker, step);
checkpointRunnerOpt.cleanupOnSuccess = checkpointOpt.cleanupOnSuccess;
checkpointRunnerOpt.cleanupOpt = checkpointOpt.cleanupOpt;
end


function taskGrid = localBuildRepeatTaskGrid(seedList)
% Build one checkpoint task per replay repeat seed.
numRepeat = numel(seedList);
taskGrid = repmat(struct('taskIndex', 0, 'repeatIndex', 0, 'taskSeed', 0), numRepeat, 1);
for iRepeat = 1:numRepeat
  taskGrid(iRepeat).taskIndex = iRepeat;
  taskGrid(iRepeat).repeatIndex = iRepeat;
  taskGrid(iRepeat).taskSeed = seedList(iRepeat);
end
end


function job = localCheckpointTaskRunner(taskInfo, sharedData)
% Run one checkpointed replay repeat task.
job = localRunOneRepeat(taskInfo.repeatIndex, taskInfo.taskSeed, ...
  sharedData.context, sharedData.repeatOpt, sharedData.flowOpt);
end


function meta = localBuildCheckpointMeta(batchOpt)
% Store compact configuration fields used to validate checkpoint compatibility.
meta = struct();
meta.replayName = string(batchOpt.replayName);
meta.snrDb = batchOpt.snrDb;
meta.baseSeed = batchOpt.baseSeed;
meta.seedList = reshape(double(batchOpt.seedList), 1, []);
meta.numRepeat = numel(batchOpt.seedList);
meta.runPeriodicWide = logical(localGetFieldOrDefault(batchOpt, 'runPeriodicWide', false));
meta.contextOpt = batchOpt.contextOpt;
meta.flowSignature = localBuildFlowSignature(batchOpt.flowOpt);
end


function flowSignature = localBuildFlowSignature(flowOpt)
% Keep scalar/string flow settings that affect replay branch behavior.
fieldList = { ...
  'periodicRefineFdHalfWidthHz', ...
  'periodicRefineFdRateHalfWidthHzPerSec', ...
  'periodicRefineDoaHalfWidthDeg', ...
  'periodicRefineFreezeDoa', ...
  'periodicRefineDoaSeedMode', ...
  'periodicPolishEnableWhenMulti', ...
  'periodicPolishDoaHalfWidthDeg', ...
  'subsetDefaultBankLabelList', ...
  'subsetMarginFallbackBankLabelList'};
flowSignature = struct();
for iField = 1:numel(fieldList)
  fieldName = fieldList{iField};
  if isfield(flowOpt, fieldName)
    flowSignature.(fieldName) = flowOpt.(fieldName);
  end
end
end


function runKey = localBuildCheckpointRunKey(batchOpt)
% Build a stable run key from the replay seed and SNR configuration.
seedList = reshape(double(batchOpt.seedList), 1, []);
runKey = string(sprintf('seed%dto%d_rep%d_snr%.2f', ...
  seedList(1), seedList(end), numel(seedList), batchOpt.snrDb));
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
end


function repoRoot = localGetRepoRoot()
% Resolve the repository root from this common replay helper location.
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(fileparts(thisFile))));
end


function runState = localEmptyRunState()
% Build an empty checkpoint run-state placeholder for non-checkpoint batches.
runState = struct();
runState.runName = "";
runState.runKey = "";
runState.outputRoot = "";
runState.runDir = "";
runState.taskDir = "";
runState.manifestFile = "";
runState.resultCell = cell(0, 1);
runState.doneMask = false(0, 1);
runState.numTask = 0;
runState.numDone = 0;
runState.isComplete = true;
runState.cleanupOnSuccess = false;
runState.cleanedOnSuccess = false;
runState.cleanupReport = struct([]);
end


function job = localRunOneRepeat(iRepeat, taskSeed, context, repeatOpt, flowOpt)
% Generate one noisy repeat and run the periodic-subset flow.
lastwarn('', '');
repeatData = buildDynamicRepeatData(context, repeatOpt.snrDb, taskSeed);
staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  repeatData.periodicFixture.viewRefOnly, repeatData.periodicFixture.viewOtherOnly, repeatData.periodicFixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  repeatData.periodicFixture.fdRange, repeatData.truth, context.otherSatIdxGlobal, ...
  repeatOpt.optVerbose, flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);

msFlow = runSimpleDynamicSubsetPeriodicFlow(repeatOpt.replayName, "multi", ...
  repeatData.periodicFixture, repeatData.subsetFixtureCell, staticBundle.caseStaticMs, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  repeatOpt.optVerbose, flowOpt);

periodicWideCase = struct();
periodicWideSummary = struct();
if repeatOpt.runPeriodicWide
  [periodicWideCase, periodicWideSummary] = localRunPeriodicWideCase( ...
    repeatData.periodicFixture, staticBundle.caseStaticMs, context, repeatOpt, flowOpt);
end

job = struct();
job.iRepeat = iRepeat;
job.taskSeed = taskSeed;
job.snrDb = repeatOpt.snrDb;
job.truth = repeatData.truth;
job.staticBundle = localStripStaticBundle(staticBundle);
job.msFlow = msFlow;
job.periodicWideCase = periodicWideCase;
job.periodicWideSummary = periodicWideSummary;
[warningMessage, warningId] = lastwarn;
job.warning = localBuildWarningSummary(warningMessage, warningId);
job.summary = localBuildJobSummary(msFlow, repeatData.truth, taskSeed, repeatOpt.snrDb, job.warning);
if ~isempty(fieldnames(periodicWideSummary))
  job.summary.periodicWideToothIdx = localGetFieldOrDefault(periodicWideSummary, 'toothIdx', NaN);
  job.summary.periodicWideAngleErrDeg = localGetFieldOrDefault(periodicWideSummary, 'angleErrDeg', NaN);
  job.summary.periodicWideFdRefErrHz = localGetFieldOrDefault(periodicWideSummary, 'fdRefErrHz', NaN);
  job.summary.periodicWideFdRateErrHzPerSec = localGetFieldOrDefault(periodicWideSummary, 'fdRateErrHzPerSec', NaN);
else
  job.summary.periodicWideToothIdx = NaN;
  job.summary.periodicWideAngleErrDeg = NaN;
  job.summary.periodicWideFdRefErrHz = NaN;
  job.summary.periodicWideFdRateErrHzPerSec = NaN;
end
end


function [caseUse, summaryUse] = localRunPeriodicWideCase(periodicFixture, staticSeedCase, context, repeatOpt, flowOpt)
% Run the full periodic fixture before subset selection for comparison.
dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = staticSeedCase.estResult.doaParamEst(:);
dynOpt.initDoaHalfWidth = flowOpt.subsetDoaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;

initParamSeed = buildDynamicInitParamFromCase(staticSeedCase, false, 0);
initCandidate = struct( ...
  'startTag', "periodic-wide-fromStatic", ...
  'initParam', initParamSeed, ...
  'initDoaParam', staticSeedCase.estResult.doaParamEst(:), ...
  'initDoaHalfWidth', flowOpt.subsetDoaHalfWidthDeg(:));

caseUse = runDynamicDoaDopplerCase("MS-MF-Dyn-periodic-wide", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, periodicFixture.fdRange, periodicFixture.fdRateRange, ...
  repeatOpt.optVerbose, dynOpt, false, struct(), initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, localResolveToothStepHz(periodicFixture), periodicFixture.truth);
summaryUse.seedSource = "periodic-wide";
summaryUse.startTag = "periodic-wide-fromStatic";
end


function toothStepHz = localResolveToothStepHz(periodicFixture)
% Resolve the fdRef tooth spacing from frame timing metadata.
toothStepHz = NaN;
if isfield(periodicFixture, 'frameIntvlSec') && isfinite(periodicFixture.frameIntvlSec) && periodicFixture.frameIntvlSec > 0
  toothStepHz = 1 / periodicFixture.frameIntvlSec;
  return;
end
if isfield(periodicFixture, 'sceneSeq') && isfield(periodicFixture.sceneSeq, 'timeOffsetSec')
  dt = diff(reshape(periodicFixture.sceneSeq.timeOffsetSec, 1, []));
  dt = dt(isfinite(dt) & dt > 0);
  if ~isempty(dt)
    toothStepHz = 1 / median(dt);
  end
end
end


function staticBundle = localStripStaticBundle(staticBundle)
% Drop heavy static fields that are not needed by replay summaries.
if isfield(staticBundle, 'weightCase')
  staticBundle.weightCase = [];
end
end


function warningSummary = localBuildWarningSummary(warningMessage, warningId)
% Store the final MATLAB warning seen during one repeat without changing results.
warningSummary = struct();
warningSummary.hasWarning = ~(isempty(warningMessage) && isempty(warningId));
warningSummary.warningId = string(warningId);
warningSummary.warningMessage = string(warningMessage);
end


function summary = localBuildJobSummary(msFlow, truth, taskSeed, snrDb, warningSummary)
% Build the per-repeat compact summary table row.
finalSummary = msFlow.finalSummary;
periodicDoaSeed = localGetFieldOrDefault(msFlow, 'periodicDoaSeed', struct());
summary = struct();
summary.taskSeed = taskSeed;
summary.snrDb = snrDb;
summary.selectedSubsetLabel = string(localGetFieldOrDefault(msFlow, 'selectedSubsetLabel', "unknown"));
summary.selectedSubsetTrusted = logical(localGetFieldOrDefault(msFlow, 'selectedSubsetTrusted', false));
summary.subsetSelectionReason = string(localGetFieldOrDefault(msFlow, 'subsetSelectionReason', "unknown"));
summary.selectedFinalTag = string(localGetFieldOrDefault(periodicDoaSeed, 'selectedReplayTag', ...
  localGetFieldOrDefault(finalSummary, 'startTag', "unknown")));
summary.selectedFinalReason = string(localGetFieldOrDefault(periodicDoaSeed, 'selectedReplaySelectionReason', "unknown"));
summary.usedPolish = logical(localGetFieldOrDefault(periodicDoaSeed, 'usedDoaPolish', false));
summary.polishSelected = logical(localGetFieldOrDefault(periodicDoaSeed, 'polishSelected', false));
summary.polishGateReason = string(localGetFieldOrDefault(periodicDoaSeed, 'polishGateReason', "unknown"));
summary.solveVariant = string(localGetFieldOrDefault(finalSummary, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(finalSummary, 'isResolved', false));
summary.angleErrDeg = localGetFieldOrDefault(finalSummary, 'angleErrDeg', NaN);
summary.fdRefErrHz = localGetFieldOrDefault(finalSummary, 'fdRefErrHz', NaN);
summary.fdRateErrHzPerSec = localGetFieldOrDefault(finalSummary, 'fdRateErrHzPerSec', NaN);
summary.toothIdx = localGetFieldOrDefault(finalSummary, 'toothIdx', NaN);
summary.toothResidualHz = localGetFieldOrDefault(finalSummary, 'toothResidualHz', NaN);
summary.finalObj = localGetFieldOrDefault(finalSummary, 'finalObj', NaN);
summary.finalResidualNorm = localGetFieldOrDefault(finalSummary, 'finalResidualNorm', NaN);
summary.nonRefCoherenceFloor = localGetFieldOrDefault(finalSummary, 'nonRefCoherenceFloor', NaN);
summary.nonRefRmsPhaseResidRad = localGetFieldOrDefault(finalSummary, 'nonRefRmsPhaseResidRad', NaN);
summary.runTimeMs = localGetFieldOrDefault(finalSummary, 'runTimeMs', NaN);
summary.warningSeen = logical(localGetFieldOrDefault(warningSummary, 'hasWarning', false));
summary.warningId = string(localGetFieldOrDefault(warningSummary, 'warningId', ""));
summary.warningMessage = string(localGetFieldOrDefault(warningSummary, 'warningMessage', ""));
summary.truthFdRefHz = localGetFieldOrDefault(truth, 'fdRefFit', localGetFieldOrDefault(truth, 'fdRefTrueHz', NaN));
summary.truthFdRateHzPerSec = localGetFieldOrDefault(truth, 'fdRateFit', localGetFieldOrDefault(truth, 'fdRateTrueHzPerSec', NaN));
end


function repeatTable = localBuildRepeatTable(repeatCell)
% Convert repeat summaries to a compact table for replay scripts.
numRepeat = numel(repeatCell);
rowCell = cell(numRepeat, 1);
for iRepeat = 1:numRepeat
  rowCell{iRepeat} = repeatCell{iRepeat}.summary;
  rowCell{iRepeat} = localFillField(rowCell{iRepeat}, 'warningSeen', false);
  rowCell{iRepeat} = localFillField(rowCell{iRepeat}, 'warningId', "");
  rowCell{iRepeat} = localFillField(rowCell{iRepeat}, 'warningMessage', "");
end
repeatTable = struct2table([rowCell{:}].');
end


function tracker = localCreateProgressTracker(titleText, totalCount, useParfor)
% Create a client-side progressbar tracker for serial or parfor repeats.
tracker = struct('active', false, 'queue', [], 'totalCount', totalCount);
if totalCount <= 0
  return;
end
fprintf('%s\n', titleText);
if exist('progressbar', 'file') ~= 2
  localLog('progressbar.m not found; continuing without progress bar.');
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
  localLog(sprintf('Progress bar disabled: %s', ME.message));
end
end


function localAdvanceProgressTracker(tracker)
% Advance the progressbar only in the serial client loop.
if tracker.active && isempty(tracker.queue)
  progressbar('advance');
end
end


function localCloseProgressTracker(tracker)
% Close an active progressbar without touching replay results.
if tracker.active
  pause(0.05);
  progressbar('end');
end
end


function localAdvanceProgressByStep(tracker, step)
% Advance the progressbar by a checkpoint runner callback step count.
if ~tracker.active
  return;
end
for iStep = 1:max(1, step)
  progressbar('advance');
end
end


function numDoneTask = localCountCheckpointTaskFile(taskDir, numTask)
% Count completed checkpoint task files before resetting the progressbar.
numDoneTask = 0;
if ~isfolder(taskDir)
  return;
end
for iTask = 1:numTask
  taskFile = fullfile(taskDir, sprintf('task_%06d.mat', iTask));
  if isfile(taskFile)
    numDoneTask = numDoneTask + 1;
  end
end
end


function outStruct = localFillField(inStruct, fieldName, defaultValue)
% Assign a default field value only when it is missing or empty.
outStruct = inStruct;
if ~isfield(outStruct, fieldName) || isempty(outStruct.(fieldName))
  outStruct.(fieldName) = defaultValue;
end
end


function outStruct = localMergeStruct(baseStruct, overrideStruct)
% Apply a shallow override struct while preserving unspecified defaults.
outStruct = baseStruct;
if nargin < 2 || isempty(overrideStruct)
  return;
end
fieldList = fieldnames(overrideStruct);
for iField = 1:numel(fieldList)
  outStruct.(fieldList{iField}) = overrideStruct.(fieldList{iField});
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
% Read a nonempty struct field or return the supplied default value.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end


function localLog(messageText)
% Print a compact timestamped replay orchestration log line.
logTime = char(datetime('now', 'Format', 'HH:mm:ss'));
fprintf('[%s] %s\n', logTime, char(messageText));
end


function tf = localCanUseParfor()
% Check whether an outer parfor repeat loop can be used safely.
tf = false;
try
  tf = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end
