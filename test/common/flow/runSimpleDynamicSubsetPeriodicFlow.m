function flow = runSimpleDynamicSubsetPeriodicFlow(displayName, satMode, periodicFixture, subsetFixtureCell, ...
  staticSeedCase, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt)
%RUNSIMPLEDYNAMICSUBSETPERIODICFLOW Run the simplified unknown-rate dynamic flow.
% The simplified flow keeps only three stages:
%   1) take one static seed;
%   2) run unknown-rate free search on anti-periodic subsets;
%   3) return to the periodic full-data window and compare narrow frozen
%      same-tooth replays before optionally running one very-small DoA
%      polish on clearly trusted same-tooth hard cases.
% No known-rate dynamic anchor, wide periodic fallback, or blanket
% same-tooth rescue is used here.

arguments
  displayName (1,1) string
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  periodicFixture (1,1) struct
  subsetFixtureCell cell
  staticSeedCase (1,1) struct
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
end

bundleTimer = tic;
toothStepHz = localResolveToothStepHz(periodicFixture);
[truthUse, debugTruth] = localResolveFlowRole(satMode, periodicFixture);

stageTimer = tic;
[subsetCaseCell, subsetSummaryCell, subsetDecisionSummaryCell, subsetRunTimeSec, useParfor] = localEvaluateSubsetBank( ...
  displayName, satMode, subsetFixtureCell, staticSeedCase, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse, toothStepHz);
stageTiming.subsetSelectSec = toc(stageTimer);

[bestSubsetIdx, subsetTrusted, subsetSelectionReason] = selectBestDynamicSubsetSummary(subsetDecisionSummaryCell);
subsetOrderIdx = localRankSubsetSummaryCell(subsetDecisionSummaryCell);
selectedSubsetFixture = subsetFixtureCell{bestSubsetIdx};
selectedSubsetCase = subsetCaseCell{bestSubsetIdx};
selectedSubsetSummary = subsetSummaryCell{bestSubsetIdx};
selectedSubsetDecisionSummary = subsetDecisionSummaryCell{bestSubsetIdx};
subsetCandidateTable = localBuildSubsetCandidateTable(subsetSummaryCell, subsetOrderIdx);
subsetTrustDiag = localBuildSubsetTrustDiag(subsetDecisionSummaryCell, subsetOrderIdx, subsetTrusted, subsetSelectionReason);

stageTimer = tic;
fdRangePeriodicRefine = [selectedSubsetCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangePeriodicRefine = [selectedSubsetCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];
[seedCandidateList, periodicCaseCell, periodicSummaryCell, periodicDecisionSummaryCell] = runSimplePeriodicRefineCandidates( ...
  displayName, satMode, periodicFixture, staticSeedCase, selectedSubsetCase, ...
  toothStepHz, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, truthUse);
[finalCase, periodicDoaSeed, periodicCaseCell, periodicSummaryCell, periodicDecisionSummaryCell] = selectSimpleDynamicFinalCase( ...
  displayName, satMode, periodicFixture, selectedSubsetCase, subsetTrustDiag, toothStepHz, ...
  seedCandidateList, periodicCaseCell, periodicSummaryCell, periodicDecisionSummaryCell, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse);
periodicCandidateTable = localBuildPeriodicCandidateTable(periodicSummaryCell);
finalSummary = buildDynamicUnknownCaseSummary(finalCase, toothStepHz, truthUse);
stageTiming.periodicRefineSec = toc(stageTimer);

stageTiming.totalFlowSec = toc(bundleTimer);

flow = struct();
flow.displayName = displayName;
flow.satMode = satMode;
flow.caseSubsetCell = subsetCaseCell;
flow.subsetSummaryCell = subsetSummaryCell;
flow.subsetDecisionSummaryCell = subsetDecisionSummaryCell;
flow.subsetRunTimeSec = subsetRunTimeSec;
flow.bestSubsetIdx = bestSubsetIdx;
flow.selectedSubsetLabel = string(localGetFieldOrDefault(selectedSubsetFixture, 'subsetLabel', "subset" + string(bestSubsetIdx)));
flow.selectedSubsetOffsets = reshape(localGetFieldOrDefault(selectedSubsetFixture, 'subsetOffsetIdx', []), 1, []);
flow.selectedSubsetCase = selectedSubsetCase;
flow.selectedSubsetSummary = selectedSubsetSummary;
flow.selectedSubsetDecisionSummary = selectedSubsetDecisionSummary;
flow.subsetCandidateTable = subsetCandidateTable;
flow.subsetTrustDiag = subsetTrustDiag;
flow.subsetSelectionReason = string(localGetFieldOrDefault(subsetTrustDiag, 'selectionReason', "unknown"));
flow.selectedSubsetTrusted = logical(localGetFieldOrDefault(subsetTrustDiag, 'isTrusted', false));
flow.fdRangePeriodicRefine = fdRangePeriodicRefine;
flow.fdRateRangePeriodicRefine = fdRateRangePeriodicRefine;
flow.periodicDoaSeed = periodicDoaSeed;
flow.periodicCaseCell = periodicCaseCell;
flow.periodicSummaryCell = periodicSummaryCell;
flow.periodicDecisionSummaryCell = periodicDecisionSummaryCell;
flow.periodicCandidateTable = periodicCandidateTable;
flow.caseFinal = finalCase;
flow.finalSummary = finalSummary;
flow.stageTiming = stageTiming;
flow.usedSubsetParfor = useParfor;
end

function caseUse = localRunSubsetUnknownCase(displayName, satMode, fixtureUse, staticSeedCase, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth)
[viewUse, truthUse] = localResolveViewAndTruth(satMode, fixtureUse);

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = fixtureUse.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = staticSeedCase.estResult.doaParamEst(:);
dynOpt.initDoaHalfWidth = flowOpt.subsetDoaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
if satMode == "multi"
  dynOpt.continuousPhaseConsistencyWeight = localGetFieldOrDefault(flowOpt.dynBaseOpt, 'continuousPhaseConsistencyWeight', 0.05);
  dynOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt.dynBaseOpt, 'continuousPhaseCollapsePenaltyWeight', 0.10);
  dynOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt.dynBaseOpt, 'continuousPhaseNegativeProjectionPenaltyWeight', 0.10);
end

initParamSeed = buildDynamicInitParamFromCase(staticSeedCase, false, 0);
initCandidate = struct( ...
  'startTag', "fromStatic", ...
  'initParam', initParamSeed, ...
  'initDoaParam', staticSeedCase.estResult.doaParamEst(:), ...
  'initDoaHalfWidth', flowOpt.subsetDoaHalfWidthDeg(:));

caseUse = runDynamicDoaDopplerCase(displayName, satMode, viewUse, truthUse, pilotWave, ...
  carrierFreq, sampleRate, fixtureUse.fdRange, fixtureUse.fdRateRange, ...
  optVerbose, dynOpt, false, debugTruth, initCandidate);
end

function satMode = localInferSatMode(viewUse)
satMode = "multi";
if isstruct(viewUse) && isfield(viewUse, 'numSat')
  if isequal(viewUse.numSat, 1)
    satMode = "single";
  end
end
end

function subsetTrustDiag = localBuildSubsetTrustDiag(subsetSummaryCell, orderIdx, isTrusted, selectionReason)
subsetTrustDiag = struct();
subsetTrustDiag.isTrusted = logical(isTrusted);
subsetTrustDiag.selectionReason = string(selectionReason);
subsetTrustDiag.selectedSubsetFinalObj = NaN;
subsetTrustDiag.runnerUpFinalObj = NaN;
subsetTrustDiag.relativeMargin = NaN;
subsetTrustDiag.runnerUpLabel = "none";
subsetTrustDiag.selectedSubsetHealthBucket = NaN;
subsetTrustDiag.runnerUpHealthBucket = NaN;
if isempty(orderIdx)
  return;
end
bestIdx = orderIdx(1);
selectedSummary = subsetSummaryCell{bestIdx};
subsetTrustDiag.selectedSubsetFinalObj = localGetFieldOrDefault(selectedSummary, 'finalObj', NaN);
subsetTrustDiag.selectedSubsetHealthBucket = localBuildPeriodicHealthBucket(selectedSummary);
if numel(orderIdx) < 2
  return;
end
runnerUpIdx = orderIdx(2);
runnerUpSummary = subsetSummaryCell{runnerUpIdx};
subsetTrustDiag.runnerUpFinalObj = localGetFieldOrDefault(runnerUpSummary, 'finalObj', NaN);
subsetTrustDiag.relativeMargin = localCalcRelativeObjGap(subsetTrustDiag.selectedSubsetFinalObj, subsetTrustDiag.runnerUpFinalObj);
subsetTrustDiag.runnerUpLabel = string(localGetFieldOrDefault(runnerUpSummary, 'subsetLabel', "subset" + string(runnerUpIdx)));
subsetTrustDiag.runnerUpHealthBucket = localBuildPeriodicHealthBucket(runnerUpSummary);
end

function relGap = localCalcRelativeObjGap(objA, objB)
relGap = NaN;
if ~isfinite(objA) || ~isfinite(objB)
  return;
end
scale = max([1, abs(objA), abs(objB)]);
relGap = abs(objA - objB) / scale;
end

function toothStepHz = localResolveToothStepHz(periodicFixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve 1/T_f from the fixture without assuming one field path.

toothStepHz = NaN;
if isstruct(periodicFixture) && isfield(periodicFixture, 'frameIntvlSec')
  rawValue = periodicFixture.frameIntvlSec;
  if ~isempty(rawValue) && isfinite(rawValue) && (rawValue > 0)
    toothStepHz = 1 / rawValue;
    return;
  end
end

if isstruct(periodicFixture) && isfield(periodicFixture, 'sceneSeq')
  sceneSeq = periodicFixture.sceneSeq;
  if isstruct(sceneSeq) && isfield(sceneSeq, 'frameIntvlSec')
    rawValue = sceneSeq.frameIntvlSec;
    if ~isempty(rawValue) && isfinite(rawValue) && (rawValue > 0)
      toothStepHz = 1 / rawValue;
      return;
    end
  end

  if isstruct(sceneSeq) && isfield(sceneSeq, 'timeOffsetSec')
    timeOffsetSec = reshape(sceneSeq.timeOffsetSec, 1, []);
    dt = diff(timeOffsetSec);
    dt = dt(isfinite(dt) & (dt > 0));
    if ~isempty(dt)
      frameIntvlSec = median(dt);
      if isfinite(frameIntvlSec) && (frameIntvlSec > 0)
        toothStepHz = 1 / frameIntvlSec;
        return;
      end
    end
  end
end

error('runSimpleDynamicSubsetPeriodicFlow:MissingFrameInterval', ...
  'Cannot resolve frame interval from periodicFixture. Expected fixture.frameIntvlSec, sceneSeq.frameIntvlSec, or sceneSeq.timeOffsetSec.');
end

function [truthUse, debugTruth] = localResolveFlowRole(satMode, periodicFixture)
[~, truthUse] = localResolveViewAndTruth(satMode, periodicFixture);
debugTruth = struct();
end

function [viewUse, truthUse] = localResolveViewAndTruth(satMode, fixtureUse)
if satMode == "single"
  viewUse = fixtureUse.viewRefOnly;
  truthUse = fixtureUse.truth;
else
  viewUse = fixtureUse.viewMs;
  truthUse = fixtureUse.truth;
end
end

function tf = localShouldUseSubsetParfor(parallelOpt, numSubset)
enableParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableSubsetEvalParfor', false));
minSubset = localGetFieldOrDefault(parallelOpt, 'minSubsetEvalParfor', 4);
tf = enableParfor && (numSubset >= minSubset) && localCanUseParfor();
end

function tableOut = localBuildSubsetCandidateTable(summaryCell, orderIdx)
numCase = numel(summaryCell);
subsetLabel = strings(numCase, 1);
fdRefEst = nan(numCase, 1);
fdRateEst = nan(numCase, 1);
finalObj = nan(numCase, 1);
finalResidualNorm = nan(numCase, 1);
runTimeMs = nan(numCase, 1);
healthBucket = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
for iCase = 1:numCase
  summaryUse = summaryCell{iCase};
  subsetLabel(iCase) = string(localGetFieldOrDefault(summaryUse, 'subsetLabel', "subset" + string(iCase)));
  fdRefEst(iCase) = localGetFieldOrDefault(summaryUse, 'fdRefEst', NaN);
  fdRateEst(iCase) = localGetFieldOrDefault(summaryUse, 'fdRateEst', NaN);
  finalObj(iCase) = localGetFieldOrDefault(summaryUse, 'finalObj', NaN);
  finalResidualNorm(iCase) = localGetFieldOrDefault(summaryUse, 'finalResidualNorm', NaN);
  runTimeMs(iCase) = localGetFieldOrDefault(summaryUse, 'runTimeMs', NaN);
  healthBucket(iCase) = localBuildPeriodicHealthBucket(summaryUse);
  toothIdx(iCase) = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
  toothResidualHz(iCase) = localGetFieldOrDefault(summaryUse, 'toothResidualHz', NaN);
end
rankVec = nan(numCase, 1);
rankVec(orderIdx) = (1:numCase).';
tableOut = table(rankVec, subsetLabel, toothIdx, toothResidualHz, healthBucket, fdRefEst, fdRateEst, finalObj, finalResidualNorm, runTimeMs, ...
  'VariableNames', {'rank', 'subsetLabel', 'toothIdx', 'toothResidualHz', 'healthBucket', 'fdRefEstHz', 'fdRateEstHzPerSec', 'finalObj', 'finalResidualNorm', 'runTimeMs'});
end

function tableOut = localBuildPeriodicCandidateTable(summaryCell)
numCase = numel(summaryCell);
rankVec = nan(numCase, 1);
seedSource = strings(numCase, 1);
startTag = strings(numCase, 1);
fdRefEst = nan(numCase, 1);
fdRateEst = nan(numCase, 1);
finalObj = nan(numCase, 1);
finalResidualNorm = nan(numCase, 1);
runTimeMs = nan(numCase, 1);
scoreMat = nan(numCase, 18);
for iCase = 1:numCase
  scoreMat(iCase, :) = buildDynamicSelectionScoreVector(summaryCell{iCase});
  seedSource(iCase) = string(localGetFieldOrDefault(summaryCell{iCase}, 'seedSource', "unknown"));
  startTag(iCase) = string(localGetFieldOrDefault(summaryCell{iCase}, 'startTag', "unknown"));
  fdRefEst(iCase) = localGetFieldOrDefault(summaryCell{iCase}, 'fdRefEst', NaN);
  fdRateEst(iCase) = localGetFieldOrDefault(summaryCell{iCase}, 'fdRateEst', NaN);
  finalObj(iCase) = localGetFieldOrDefault(summaryCell{iCase}, 'finalObj', NaN);
  finalResidualNorm(iCase) = localGetFieldOrDefault(summaryCell{iCase}, 'finalResidualNorm', NaN);
  runTimeMs(iCase) = localGetFieldOrDefault(summaryCell{iCase}, 'runTimeMs', NaN);
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
rankVec(orderIdx) = (1:numCase).';
tableOut = table(rankVec, seedSource, startTag, fdRefEst, fdRateEst, finalObj, finalResidualNorm, runTimeMs, ...
  'VariableNames', {'rank', 'seedSource', 'startTag', 'fdRefEstHz', 'fdRateEstHzPerSec', 'finalObj', 'finalResidualNorm', 'runTimeMs'});
end

function [subsetCaseCell, subsetSummaryCell, subsetDecisionSummaryCell, subsetRunTimeSec, useParfor] = localEvaluateSubsetBank( ...
  displayName, satMode, subsetFixtureCell, staticSeedCase, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse, toothStepHz)

runSubsetCaseFn = @(fixtureUse, pilotWaveUse, carrierFreqUse, sampleRateUse, optVerboseUse, flowOptUse, fdRangeIgnored, fdRateRangeIgnored, subsetSeedInfoIgnored) ...
  localRunSubsetUnknownCase(displayName + "-subset", satMode, fixtureUse, staticSeedCase, ...
    pilotWaveUse, carrierFreqUse, sampleRateUse, optVerboseUse, flowOptUse, debugTruth);
buildSummaryFn = @(caseUse, truthLocal, toothStepLocal) ...
  buildDynamicUnknownCaseSummary(caseUse, toothStepLocal, struct());
subsetSeedInfo = struct();
numSubsetRequested = numel(subsetFixtureCell);
useParfor = localShouldUseSubsetParfor(localGetFieldOrDefault(flowOpt, 'parallelOpt', struct()), numSubsetRequested);
[subsetCaseCell, subsetDecisionSummaryCell, subsetRunTimeSec] = evaluateDynamicSubsetBank( ...
  subsetFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  [NaN, NaN], [NaN, NaN], subsetSeedInfo, truthUse, toothStepHz, ...
  runSubsetCaseFn, buildSummaryFn);
subsetSummaryCell = localBuildSubsetEvalSummaryCell(subsetCaseCell, subsetFixtureCell, toothStepHz);
end



function summaryCell = localBuildSubsetEvalSummaryCell(caseCell, subsetFixtureCell, toothStepHz)
%LOCALBUILDSUBSETEVALSUMMARYCELL Build truth-aware summaries for reporting only.

numCase = numel(caseCell);
summaryCell = cell(numCase, 1);
for iCase = 1:numCase
  fixtureUse = subsetFixtureCell{iCase};
  summaryCell{iCase} = buildDynamicUnknownCaseSummary(caseCell{iCase}, toothStepHz, fixtureUse.truth);
  summaryCell{iCase}.subsetLabel = string(localGetFieldOrDefault(fixtureUse, ...
    'subsetLabel', "subset" + string(iCase)));
  summaryCell{iCase}.subsetOffsetIdx = reshape(localGetFieldOrDefault(fixtureUse, ...
    'subsetOffsetIdx', []), 1, []);
end
end

function orderIdx = localRankSubsetSummaryCell(summaryCell)
numCase = numel(summaryCell);
orderIdx = nan(numCase, 1);
remainingIdx = 1:numCase;
for iRank = 1:numCase
  [bestLocalIdx, ~] = selectBestDynamicSubsetSummary(summaryCell(remainingIdx));
  orderIdx(iRank) = remainingIdx(bestLocalIdx);
  remainingIdx(bestLocalIdx) = [];
end
end


function bucket = localBuildPeriodicHealthBucket(summary)
%LOCALBUILDPERIODICHEALTHBUCKET Count coarse non-reference health failures.

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


function tf = localCanUseParfor()
tf = false;
try
  tf = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
