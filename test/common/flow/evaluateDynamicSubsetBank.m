function [subsetCaseCell, subsetSummaryCell, subsetRunTimeSec] = evaluateDynamicSubsetBank( ...
  subsetFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  fdRangeUse, fdRateRangeUse, subsetSeedInfo, truth, toothStepHz, ...
  runSubsetCaseFn, buildSummaryFn)
%EVALUATEDYNAMICSUBSETBANK Evaluate one subset-fixture bank with optional parfor.
% This helper keeps the dynamic transition bundle entry focused on flow
% orchestration while centralizing the subset-bank loop, early-stop guard,
% and safe parfor gating. The actual CP-U subset solve and summary build
% stay injectable through function handles so the shared flow can keep its
% current behavior.

arguments
  subsetFixtureCell cell
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical
  flowOpt (1, 1) struct
  fdRangeUse (1, 2) double
  fdRateRangeUse (1, 2) double
  subsetSeedInfo (1, 1) struct
  truth (1, 1) struct
  toothStepHz (1, 1) double
  runSubsetCaseFn (1, 1) function_handle
  buildSummaryFn (1, 1) function_handle
end

numSubset = numel(subsetFixtureCell);
subsetCaseCell = cell(numSubset, 1);
subsetSummaryCell = cell(numSubset, 1);
subsetRunTimeSec = nan(numSubset, 1);
useParfor = localShouldUseSubsetEvalParfor(flowOpt, numSubset);

if useParfor
  parfor iSubset = 1:numSubset
    fixtureUse = subsetFixtureCell{iSubset};
    fixtureUse.debugTruthMs = buildTruthDebugForView(fixtureUse.viewMs, fixtureUse.truth);
    subsetTimer = tic;
    subsetCaseCell{iSubset} = runSubsetCaseFn(fixtureUse, pilotWave, ...
      carrierFreq, sampleRate, optVerbose, flowOpt, fdRangeUse, fdRateRangeUse, subsetSeedInfo);
    subsetRunTimeSec(iSubset) = toc(subsetTimer);
    subsetSummaryCell{iSubset} = buildSummaryFn(subsetCaseCell{iSubset}, fixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(fixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(fixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
  return;
end

numEval = 0;
for iSubset = 1:numSubset
  if localShouldEarlyStopSubsetBank(subsetFixtureCell, subsetSummaryCell, numEval, flowOpt)
    break;
  end
  fixtureUse = subsetFixtureCell{iSubset};
  fixtureUse.debugTruthMs = buildTruthDebugForView(fixtureUse.viewMs, fixtureUse.truth);
  subsetTimer = tic;
  subsetCaseCell{iSubset} = runSubsetCaseFn(fixtureUse, pilotWave, ...
    carrierFreq, sampleRate, optVerbose, flowOpt, fdRangeUse, fdRateRangeUse, subsetSeedInfo);
  subsetRunTimeSec(iSubset) = toc(subsetTimer);
  subsetSummaryCell{iSubset} = buildSummaryFn(subsetCaseCell{iSubset}, fixtureUse.truth, toothStepHz);
  subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(fixtureUse, ...
    'subsetLabel', "subset" + string(iSubset)));
  subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(fixtureUse, ...
    'subsetOffsetIdx', []), 1, []);
  numEval = iSubset;
end
subsetCaseCell = subsetCaseCell(1:numEval);
subsetSummaryCell = subsetSummaryCell(1:numEval);
subsetRunTimeSec = subsetRunTimeSec(1:numEval);
end


function tf = localShouldEarlyStopSubsetBank(subsetFixtureCell, subsetSummaryCell, numEval, flowOpt)
%LOCALSHOULDEARLYSTOPSUBSETBANK Skip late random fallback when curated bank is already healthy.
% Fast statistics often append one or more random-labeled subset fixtures as
% a last-resort tooth selector. When the deterministic curated bank has
% already produced one trusted central-tooth candidate with healthy
% non-reference phase metrics, evaluating the trailing random fallback
% rarely changes the winner but still costs another full CP-U solve.

tf = false;
if numEval <= 0
  return;
end
if ~logical(localGetFieldOrDefault(flowOpt, 'enableSubsetRandomEarlyStop', true))
  return;
end
if numEval >= numel(subsetFixtureCell)
  return;
end
minEval = localGetFieldOrDefault(flowOpt, 'subsetRandomEarlyStopMinEvaluated', 4);
if numEval < minEval
  return;
end
nextFixture = subsetFixtureCell{numEval + 1};
nextLabel = lower(strtrim(char(string(localGetFieldOrDefault(nextFixture, 'subsetLabel', "")))));
if ~(startsWith(nextLabel, 'random') || startsWith(nextLabel, 'rescue-random'))
  return;
end

summaryEvalCell = subsetSummaryCell(1:numEval);
[bestIdx, isTrusted] = selectBestDynamicSubsetSummary(summaryEvalCell);
if ~isTrusted || ~isfinite(bestIdx)
  return;
end
bestSummary = summaryEvalCell{bestIdx};
toothIdx = localGetFieldOrDefault(bestSummary, 'toothIdx', NaN);
toothResidualHz = abs(localGetFieldOrDefault(bestSummary, 'toothResidualHz', inf));
if ~(isfinite(toothIdx) && abs(toothIdx) == 0)
  return;
end
if ~(isfinite(toothResidualHz) && toothResidualHz <= ...
    localGetFieldOrDefault(flowOpt, 'subsetRandomEarlyStopMaxToothResidualHz', 20))
  return;
end
if localBuildSubsetHealthBucket(bestSummary) > ...
    localGetFieldOrDefault(flowOpt, 'subsetRandomEarlyStopMaxHealthBucket', 0)
  return;
end

tf = true;
end


function bucket = localBuildSubsetHealthBucket(summary)
%LOCALBUILDSUBSETHEALTHBUCKET Count coarse non-reference quality failures.

bucket = 0;
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefSupportRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefFitRatioFloor', 1) < 0.95);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefConsistencyRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefCoherenceFloor', 1) < 0.995);
rmsPhaseResid = localGetFieldOrDefault(summary, 'nonRefRmsPhaseResidRad', NaN);
if isfinite(rmsPhaseResid)
  bucket = bucket + double(rmsPhaseResid > 0.003);
end
maxPhaseResid = localGetFieldOrDefault(summary, 'nonRefMaxAbsPhaseResidRad', NaN);
if isfinite(maxPhaseResid)
  bucket = bucket + double(maxPhaseResid > 0.005);
end
negativeRatio = localGetFieldOrDefault(summary, 'maxNonRefNegativeProjectionRatio', NaN);
if isfinite(negativeRatio)
  bucket = bucket + double(negativeRatio > 0.05);
end
end


function useParfor = localShouldUseSubsetEvalParfor(flowOpt, numSubset)
%LOCALSHOULDUSESUBSETEVALPARFOR Decide whether subset evaluation should use parfor.
% Only use parfor outside an existing task. This keeps the Monte Carlo outer
% loop free of nested parfor while still accelerating representative
% full-diag subset-bank evaluation.

parallelOpt = localGetFieldOrDefault(flowOpt, 'parallelOpt', struct());
useParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableSubsetEvalParfor', false));
if ~useParfor
  return;
end
minSubset = localGetFieldOrDefault(parallelOpt, 'minSubsetEvalParfor', 4);
useParfor = (numSubset >= minSubset) && localCanUseParfor();
end


function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether parfor can be used safely here.

tf = false;
if ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  tf = ~isempty(ver('parallel'));
catch
  tf = false;
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or a default value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
