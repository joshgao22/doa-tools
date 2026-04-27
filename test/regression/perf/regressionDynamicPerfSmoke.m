function regressionDynamicPerfSmoke(varargin)
%REGRESSIONDYNAMICPERFSMOKE Smoke-test dynamic perf orchestration.
% This regression keeps the task grid tiny and focuses on the orchestration
% contracts:
%   1) checkpoint manifest/task files are written;
%   2) resume continues from partial task output instead of restarting;
%   3) compact summary and CRB tables can be built from the resumed result;
%   4) a small saved result file is emitted successfully.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

baseSeed = 253;
snrDbList = 10;
numRepeat = 2;
optVerbose = verbose;
parallelOpt = struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', 4);
context = buildDynamicDualSatEciContext(struct( ...
  'baseSeed', baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt));
flowOpt = localBuildPerfFlowOpt();
caseName = localBuildCaseName(flowOpt);
numCase = numel(caseName);
idxMsDynUnknown = 8;

[repeatGrid, snrGrid] = ndgrid(1:numRepeat, 1:numel(snrDbList));
repeatTaskIdx = repeatGrid(:);
snrTaskIdx = snrGrid(:);
numTask = numel(repeatTaskIdx);
taskGrid = repmat(struct('taskIndex', 0, 'snrIndex', 0, 'repeatIndex', 0, 'taskSeed', 0), numTask, 1);
for iTask = 1:numTask
  taskGrid(iTask).taskIndex = iTask;
  taskGrid(iTask).snrIndex = snrTaskIdx(iTask);
  taskGrid(iTask).repeatIndex = repeatTaskIdx(iTask);
  taskGrid(iTask).taskSeed = baseSeed + iTask - 1;
end

sharedData = struct();
sharedData.context = context;
sharedData.flowOpt = flowOpt;
sharedData.snrDbList = snrDbList(:).';
sharedData.caseName = caseName(:).';
sharedData.numCase = numCase;
sharedData.idxMsDynUnknown = idxMsDynUnknown;
sharedData.optVerbose = optVerbose;
sharedData.weightSweepAlpha = reshape(flowOpt.weightSweepAlpha, [], 1);

repoRoot = localGetRepoRoot();
outputRoot = fullfile(repoRoot, 'tmp');
runKey = localBuildCheckpointRunKey(baseSeed, snrDbList, numRepeat, numCase, context.selectedSatIdxGlobal);
checkpointOpt = struct();
checkpointOpt.runName = "regressionDynamicPerfSmoke";
checkpointOpt.runKey = runKey;
checkpointOpt.outputRoot = outputRoot;
checkpointOpt.useParfor = false;
checkpointOpt.resume = true;
checkpointOpt.meta = struct('snrDbList', snrDbList(:).', 'numRepeat', numRepeat, ...
  'numCase', numCase, 'selectedSatIdxGlobal', context.selectedSatIdxGlobal(:).');

runDir = fullfile(outputRoot, checkpointOpt.runName, checkpointOpt.runKey);
if isfolder(runDir)
  cleanupRunArtifacts(runDir, struct('requiredPrefix', 'seed', 'verbose', false));
end

  fprintf('Running regressionDynamicPerfSmoke\n');

firstPassFailed = false;
try
  runPerfTaskGridWithCheckpoint(taskGrid, sharedData, @localPartialTaskRunner, checkpointOpt);
catch ME
  firstPassFailed = true;
    fprintf('  expected partial-run failure   : %s\n', char(string(ME.identifier)));
end
if ~firstPassFailed
  error('regressionDynamicPerfSmoke:ExpectedPartialFailure', ...
    'The first-pass smoke run should stop early to test resume.');
end

taskDir = fullfile(runDir, 'task');
numDoneTask = localCountCheckpointTaskFile(taskDir, numTask);
if numDoneTask ~= 1
  error('regressionDynamicPerfSmoke:UnexpectedPartialTaskCount', ...
    'Expected exactly one completed task after the partial run, got %d.', numDoneTask);
end

runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, @runDynamicPerfTask, checkpointOpt);
if ~runState.isComplete
  error('regressionDynamicPerfSmoke:ResumeIncomplete', ...
    'Checkpoint resume did not finish the tiny task grid.');
end
if localCountCheckpointTaskFile(taskDir, numTask) ~= numTask
  error('regressionDynamicPerfSmoke:MissingTaskFilesAfterResume', ...
    'Resume finished but not all task files exist.');
end

angleErrDeg = nan(numRepeat, numCase);
fdErrHz = nan(numRepeat, numCase);
fdRateErrHzPerSec = nan(numRepeat, numCase);
isResolved = false(numRepeat, numCase);
for iTask = 1:numTask
  taskOut = runState.resultCell{iTask};
  if isempty(taskOut)
    continue;
  end
  iRepeat = taskGrid(iTask).repeatIndex;
  angleErrDeg(iRepeat, :) = taskOut.angleVec;
  fdErrHz(iRepeat, :) = taskOut.fdVec;
  fdRateErrHzPerSec(iRepeat, :) = taskOut.fdRateVec;
  isResolved(iRepeat, :) = taskOut.resolvedVec;
end

summaryTable = buildDynamicCaseSummaryTable(caseName(:), snrDbList(1), ...
  angleErrDeg, fdErrHz, fdRateErrHzPerSec, isResolved);
if height(summaryTable) ~= numCase
  error('regressionDynamicPerfSmoke:UnexpectedSummaryHeight', ...
    'Summary table height mismatch: expected %d rows, got %d.', numCase, height(summaryTable));
end
if ~any(summaryTable.displayName == "MS-MF-CP-U")
  error('regressionDynamicPerfSmoke:MissingMsUnknownSummaryRow', ...
    'Summary table is missing the MS-MF-CP-U row.');
end
if ~any(summaryTable.resolveRate > 0)
  error('regressionDynamicPerfSmoke:NoResolvedSummaryRow', ...
    'Smoke summary unexpectedly contains no resolved rows.');
end

anchorRepeat = buildDynamicRepeatData(context, snrDbList(1), baseSeed);
crbBundle = buildDynamicCrbBundle(anchorRepeat.periodicFixture, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, 1 / (10^(snrDbList(1) / 10)));
crbTable = buildDynamicCrbSummaryTable(crbBundle.truth, ...
  crbBundle.crbSfRef, crbBundle.auxCrbSfRef, ...
  crbBundle.crbSfMs, crbBundle.auxCrbSfMs, ...
  crbBundle.crbMfRefKnown, crbBundle.auxCrbMfRefKnown, ...
  crbBundle.crbMfMsKnown, crbBundle.auxCrbMfMsKnown, ...
  crbBundle.crbMfRefUnknown, crbBundle.auxCrbMfRefUnknown, ...
  crbBundle.crbMfMsUnknown, crbBundle.auxCrbMfMsUnknown, snrDbList(1));
if height(crbTable) ~= 6
  error('regressionDynamicPerfSmoke:UnexpectedCrbHeight', ...
    'CRB table height mismatch: expected 6 rows, got %d.', height(crbTable));
end

smokeResultFile = fullfile(runDir, 'dynamic_perf_smoke_result.mat');
save(smokeResultFile, 'summaryTable', 'crbTable', 'runState');
if exist(smokeResultFile, 'file') ~= 2
  error('regressionDynamicPerfSmoke:SmokeSaveMissing', ...
    'Smoke result file was not written: %s', smokeResultFile);
end

  fprintf('  resumed tasks                  : %d / %d\n', runState.numDone, runState.numTask);
  fprintf('  summary rows                   : %d\n', height(summaryTable));
  fprintf('  CRB rows                       : %d\n', height(crbTable));
  fprintf('  saved result                   : %s\n', smokeResultFile);
  fprintf('PASS: regressionDynamicPerfSmoke\n');

if isfolder(runDir)
  cleanupRunArtifacts(runDir, struct('requiredPrefix', 'seed', 'verbose', false));
end
end

function flowOpt = localBuildPerfFlowOpt()
flowOpt = struct();
flowOpt.enableWeightSweep = false;
flowOpt.weightSweepAlpha = zeros(0, 1);
flowOpt.staticMsHalfWidth = [0.002; 0.002];
flowOpt.doaOnlyOpt = struct('useLogObjective', true);
flowOpt.staticBaseOpt = struct('useLogObjective', true);
flowOpt.dynBaseOpt = struct( ...
  'useLogObjective', true, ...
  'initFdCount', 81, ...
  'useAccessMask', false, ...
  'phaseMode', 'continuous', ...
  'steeringMode', 'framewise', ...
  'continuousPhaseConsistencyWeight', 0.05, ...
  'continuousPhaseCollapsePenaltyWeight', 0.10, ...
  'continuousPhaseNegativeProjectionPenaltyWeight', 0.10, ...
  'unknownWarmAnchorUseScaledSolve', true, ...
  'unknownWarmAnchorFallbackSqp', true, ...
  'debugEnable', true, ...
  'debugStoreEvalTrace', false, ...
  'debugMaxEvalTrace', 120);
flowOpt.msContinuousPhaseConsistencyWeight = flowOpt.dynBaseOpt.continuousPhaseConsistencyWeight;
flowOpt.msContinuousPhaseCollapsePenaltyWeight = flowOpt.dynBaseOpt.continuousPhaseCollapsePenaltyWeight;
flowOpt.msContinuousPhaseNegativeProjectionPenaltyWeight = flowOpt.dynBaseOpt.continuousPhaseNegativeProjectionPenaltyWeight;
flowOpt.unknownWarmAnchorUseScaledSolve = flowOpt.dynBaseOpt.unknownWarmAnchorUseScaledSolve;
flowOpt.unknownWarmAnchorFallbackSqp = flowOpt.dynBaseOpt.unknownWarmAnchorFallbackSqp;
flowOpt.refKnownDoaHalfWidth = [0.005; 0.005];
flowOpt.refUnknownDoaHalfWidth = [0.003; 0.003];
flowOpt.msKnownDoaHalfWidth = [0.003; 0.003];
flowOpt.msUnknownDoaHalfWidth = [0.002; 0.002];
flowOpt.subsetSelectDoaHalfWidthDeg = [0.01; 0.01];
flowOpt.inToothFdHalfWidthHz = 20;
flowOpt.inToothFdRateHalfWidthHzPerSec = 20;
flowOpt.inToothDoaHalfWidthDeg = [2e-4; 2e-4];
flowOpt.enableAnchorDoaPolish = false;
flowOpt.anchorDoaPolishDoaHalfWidthDeg = [0.005; 0.005];
flowOpt.anchorDoaPolishFdHalfWidthHz = 20;
flowOpt.anchorDoaPolishFdRateHalfWidthHzPerSec = 20;
flowOpt.anchorDoaPolishCentralResidualTolHz = 50;
flowOpt.enableUnknownWideBranches = true;
flowOpt.enableUnknownWideRefine = true;
flowOpt.deferUnknownWideBranches = true;
flowOpt.enableFastWideFallback = false;
flowOpt.enableFastSubsetEscalation = true;
[~, ~, rescueSubsetOffsetCell, rescueSubsetLabelList] = getDynamicCuratedSubsetBank();
flowOpt.fastRescueSubsetOffsetCell = rescueSubsetOffsetCell;
flowOpt.fastRescueSubsetLabelList = rescueSubsetLabelList;
flowOpt.fastNumRandomSubsetTrialFallback = 0;
flowOpt.fastEscalateOnToothDisagreement = false;
flowOpt.fastToothIdxThreshold = 0;
flowOpt.fastToothResidualHzThreshold = 50;
flowOpt.fastUnknownFdDriftHzThreshold = 1000;
flowOpt.fastUnknownDoaDriftDegThreshold = 0.005;
flowOpt.fastWideToothIdxThreshold = 0;
flowOpt.fastWideToothResidualHzThreshold = 50;
flowOpt.fastWideFdDriftHzThreshold = 1000;
flowOpt.fastWideDoaDriftDegThreshold = 0.005;
flowOpt.enableSubsetInToothRefine = true;
flowOpt.enableInToothCentralSkip = true;
flowOpt.inToothCentralResidualTolHz = 5;
flowOpt.enableSubsetCaseSeedReuse = true;
flowOpt.enableConditionalRandomSubsetRescue = false;
flowOpt.conditionalRandomSubsetTrialCount = 0;
flowOpt.conditionalRandomToothResidualHzThreshold = 20;
flowOpt.conditionalRandomMaxHealthBucket = 1;
flowOpt.conditionalRandomMaxFdDriftHz = 250;
flowOpt.conditionalRandomMaxDoaDriftDeg = 0.003;
flowOpt.conditionalRandomRescueOnToothDisagreement = true;
flowOpt.enableAdaptiveInToothDoaEscalation = false;
flowOpt.inToothEscalateWhenRandomSelected = true;
flowOpt.inToothEscalateMinAnchorDoaDriftFromKnownDeg = 0.006;
flowOpt.inToothEscalatedDoaHalfWidthDeg = [0.01; 0.01];
flowOpt.inToothEscalatedOptimOpt = struct('MaxIterations', 180, 'MaxFunctionEvaluations', 1800);
flowOpt.parallelOpt = struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', 4);
end

function caseName = localBuildCaseName(flowOpt)
caseName = ["SS-SF-DoA", "MS-SF-DoA", "SS-SF-Static", "MS-SF-Static", ...
  "SS-MF-CP-K", "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"];
if logical(localGetFieldOrDefault(flowOpt, 'enableWeightSweep', false))
  caseName = [caseName, compose('MS-SF-Static-W%.2f', flowOpt.weightSweepAlpha.')];
end
end

function taskOut = localPartialTaskRunner(taskInfo, sharedData)
if taskInfo.taskIndex > 1
  error('regressionDynamicPerfSmoke:IntentionalPartialStop', ...
    'Intentional stop after one task to test checkpoint resume.');
end
taskOut = runDynamicPerfTask(taskInfo, sharedData);
end

function repoRoot = localGetRepoRoot()
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
end

function runKey = localBuildCheckpointRunKey(baseSeed, snrDb, numRepeat, numCase, selectedSatIdxGlobal)
satToken = sprintf('sat%d', selectedSatIdxGlobal(1));
for iSat = 2:numel(selectedSatIdxGlobal)
  satToken = sprintf('%s_%d', satToken, selectedSatIdxGlobal(iSat));
end
runKey = string(sprintf('seed%d_%s_snr%dto%d_n%d_rep%d_case%d', ...
  baseSeed, satToken, snrDb(1), snrDb(end), numel(snrDb), numRepeat, numCase));
end

function numDoneTask = localCountCheckpointTaskFile(taskDir, numTask)
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

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end
