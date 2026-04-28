function [caseUse, periodicDoaSeed, periodicCaseCell, periodicSummaryCell, periodicDecisionSummaryCell] = selectSimpleDynamicFinalCase( ...
  displayName, satMode, periodicFixture, selectedSubsetCase, subsetTrustDiag, toothStepHz, ...
  seedCandidateList, periodicCaseCell, periodicSummaryCell, periodicDecisionSummaryCell, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse)
%SELECTSIMPLEDYNAMICFINALCASE Select the final same-tooth replay winner.
% This helper keeps the simplified flow entry focused on orchestration. It
% chooses between frozen replay candidates, optionally runs one very-small
% DoA polish, and returns the final flow case together with compact replay
% diagnostics. Selection and polish gates use decision summaries built
% without truth; truth-aware summaries are only kept for reporting.

arguments
  displayName (1,1) string
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  periodicFixture (1,1) struct
  selectedSubsetCase (1,1) struct
  subsetTrustDiag (1,1) struct
  toothStepHz (1,1) double
  seedCandidateList struct
  periodicCaseCell cell
  periodicSummaryCell cell
  periodicDecisionSummaryCell cell
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
  debugTruth (1,1) struct
  truthUse (1,1) struct
end

scoreMat = nan(numel(periodicDecisionSummaryCell), 18);
for iCand = 1:numel(periodicDecisionSummaryCell)
  scoreMat(iCand, :) = buildDynamicSelectionScoreVector(periodicDecisionSummaryCell{iCand});
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
bestFrozenIdx = orderIdx(1);
runnerUpFrozenIdx = localGetRunnerUpIndex(orderIdx);
periodicDoaSeed = localBuildPeriodicSeedDiag(seedCandidateList, periodicCaseCell, periodicDecisionSummaryCell, ...
  bestFrozenIdx, runnerUpFrozenIdx, subsetTrustDiag);

[selectedReplayIdx, replaySelectionReason] = localResolvePeriodicReplayWinner(seedCandidateList, periodicDoaSeed, flowOpt);
selectedReplaySummary = periodicSummaryCell{selectedReplayIdx};
selectedReplayDecisionSummary = periodicDecisionSummaryCell{selectedReplayIdx};
periodicDoaSeed.selectedReplayIdx = selectedReplayIdx;
periodicDoaSeed.selectedReplaySeedSource = string(seedCandidateList(selectedReplayIdx).seedSource);
periodicDoaSeed.selectedReplayTag = string(seedCandidateList(selectedReplayIdx).startTag);
periodicDoaSeed.selectedReplaySelectionReason = string(replaySelectionReason);
periodicDoaSeed.selectedReplaySummary = selectedReplaySummary;
periodicDoaSeed.selectedReplayDecisionSummary = selectedReplayDecisionSummary;
periodicDoaSeed.selectedReplayTruthToothIdx = localGetFieldOrDefault(selectedReplaySummary, 'truthToothIdx', ...
  localGetFieldOrDefault(selectedReplaySummary, 'toothIdx', NaN));
periodicDoaSeed.selectedReplayTruthToothResidualHz = localGetFieldOrDefault(selectedReplaySummary, 'truthToothResidualHz', ...
  localGetFieldOrDefault(selectedReplaySummary, 'toothResidualHz', NaN));
periodicDoaSeed.selectedReplayHealthBucket = localBuildPeriodicHealthBucket(selectedReplayDecisionSummary);

caseUse = periodicCaseCell{selectedReplayIdx};
periodicDoaSeed.usedDoaPolish = false;
[runPolish, polishGateReason, polishGateDiag] = shouldRunSimpleVerySmallDoaPolish( ...
  satMode, periodicDoaSeed, selectedReplayDecisionSummary, flowOpt);
periodicDoaSeed.polishGateReason = string(polishGateReason);
periodicDoaSeed.polishGateDiag = polishGateDiag;
if runPolish
  [polishCase, polishSummary, polishDecisionSummary] = localRunPeriodicPolishCandidate( ...
    displayName, periodicFixture, caseUse, periodicDoaSeed.selectedReplaySeedSource, ...
    periodicDoaSeed.selectedReplayTag, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, debugTruth, toothStepHz, truthUse);
  periodicCaseCell{end + 1, 1} = polishCase; %#ok<AGROW>
  periodicSummaryCell{end + 1, 1} = polishSummary; %#ok<AGROW>
  periodicDecisionSummaryCell{end + 1, 1} = polishDecisionSummary; %#ok<AGROW>
  candidateScore = buildDynamicSelectionScoreVector(polishDecisionSummary);
  selectedScore = buildDynamicSelectionScoreVector(selectedReplayDecisionSummary);
  choosePolish = localIsScoreBetter(candidateScore, selectedScore);
  periodicDoaSeed.usedDoaPolish = true;
  periodicDoaSeed.polishSelected = logical(choosePolish);
  periodicDoaSeed.polishCandidateSummary = polishSummary;
  periodicDoaSeed.polishCandidateDecisionSummary = polishDecisionSummary;
  if choosePolish
    caseUse = polishCase;
    periodicDoaSeed.selectedReplaySeedSource = string(localGetFieldOrDefault(polishSummary, 'seedSource', periodicDoaSeed.selectedReplaySeedSource));
    periodicDoaSeed.selectedReplayTag = string(localGetFieldOrDefault(polishSummary, 'startTag', "periodic-polish"));
    periodicDoaSeed.selectedReplaySelectionReason = "periodic-polish-better";
  else
    periodicDoaSeed.selectedReplaySelectionReason = "frozen-better-than-polish";
  end
else
  periodicDoaSeed.usedDoaPolish = false;
  periodicDoaSeed.polishSelected = false;
end
end

function [caseUse, summaryUse, decisionSummary] = localRunPeriodicPolishCandidate(displayName, periodicFixture, selectedReplayCase, ...
  replaySeedSource, replayTag, pilotWave, carrierFreq, sampleRate, optVerbose, ...
  flowOpt, debugTruth, toothStepHz, truthUse)
fdRangeUse = [selectedReplayCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedReplayCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangeUse = [selectedReplayCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedReplayCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = selectedReplayCase.estResult.doaParamEst(:);
dynOpt.initDoaHalfWidth = flowOpt.periodicRefinePolishDoaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = false;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = flowOpt.periodicRefinePolishDoaHalfWidthDeg(:);

initParamSeed = buildDynamicInitParamFromCase(selectedReplayCase, false, selectedReplayCase.estResult.fdRateEst);
initCandidate = struct( ...
  'startTag', "periodic-polish-" + string(replaySeedSource), ...
  'initParam', initParamSeed, ...
  'initDoaParam', selectedReplayCase.estResult.doaParamEst(:), ...
  'initDoaHalfWidth', flowOpt.periodicRefinePolishDoaHalfWidthDeg(:), ...
  'freezeDoa', false);

viewUse = localGetViewUse(periodicFixture);
caseUse = runDynamicDoaDopplerCase(displayName, localInferSatMode(viewUse), ...
  viewUse, truthUse, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, false, debugTruth, initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truthUse);
decisionSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, struct());
summaryUse.seedSource = string(replaySeedSource) + "-polish";
summaryUse.startTag = "periodic-polish-" + string(replaySeedSource);
summaryUse.replayBaseTag = string(replayTag);
summaryUse.subsetDriftFromStaticDeg = NaN;
decisionSummary.seedSource = summaryUse.seedSource;
decisionSummary.startTag = summaryUse.startTag;
decisionSummary.replayBaseTag = summaryUse.replayBaseTag;
decisionSummary.subsetDriftFromStaticDeg = NaN;
end

function viewUse = localGetViewUse(periodicFixture)
viewUse = periodicFixture.viewMs;
if isfield(periodicFixture, 'viewRefOnly') && isfield(periodicFixture.viewMs, 'numSat') && isequal(periodicFixture.viewMs.numSat, 1)
  viewUse = periodicFixture.viewRefOnly;
end
end

function periodicDoaSeed = localBuildPeriodicSeedDiag(seedCandidateList, periodicCaseCell, periodicSummaryCell, ...
  bestFrozenIdx, runnerUpFrozenIdx, subsetTrustDiag)
periodicDoaSeed = struct();
bestCandidate = seedCandidateList(bestFrozenIdx);
bestSummary = periodicSummaryCell{bestFrozenIdx};
periodicDoaSeed.doaInitParam = reshape(bestCandidate.doaInitParam(:), 1, []);
periodicDoaSeed.seedSource = string(bestCandidate.seedSource);
periodicDoaSeed.subsetDriftFromStaticDeg = bestCandidate.subsetDriftFromStaticDeg;
periodicDoaSeed.bestFrozenIdx = bestFrozenIdx;
periodicDoaSeed.bestFrozenSeedSource = string(bestCandidate.seedSource);
periodicDoaSeed.bestFrozenTag = string(bestCandidate.startTag);
periodicDoaSeed.bestFrozenFinalObj = localGetFieldOrDefault(bestSummary, 'finalObj', NaN);
periodicDoaSeed.selectedReplaySelectionReason = "best-frozen";
periodicDoaSeed.runnerUpFrozenIdx = runnerUpFrozenIdx;
periodicDoaSeed.subsetRankMarginRelative = localGetFieldOrDefault(subsetTrustDiag, 'relativeMargin', NaN);
periodicDoaSeed.subsetRunnerUpLabel = string(localGetFieldOrDefault(subsetTrustDiag, 'runnerUpLabel', "none"));
periodicDoaSeed.frozenDoaDisagreementDeg = NaN;
periodicDoaSeed.frozenRelativeObjGap = NaN;
periodicDoaSeed.frozenRunnerUpSeedSource = "none";
periodicDoaSeed.frozenRunnerUpTag = "none";
if ~isnan(runnerUpFrozenIdx)
  runnerUpSummary = periodicSummaryCell{runnerUpFrozenIdx};
  runnerUpCase = periodicCaseCell{runnerUpFrozenIdx};
  bestCase = periodicCaseCell{bestFrozenIdx};
  periodicDoaSeed.frozenRunnerUpSeedSource = string(seedCandidateList(runnerUpFrozenIdx).seedSource);
  periodicDoaSeed.frozenRunnerUpTag = string(seedCandidateList(runnerUpFrozenIdx).startTag);
  periodicDoaSeed.frozenRelativeObjGap = localCalcRelativeObjGap( ...
    localGetFieldOrDefault(bestSummary, 'finalObj', NaN), ...
    localGetFieldOrDefault(runnerUpSummary, 'finalObj', NaN));
  periodicDoaSeed.frozenDoaDisagreementDeg = localCalcCaseDoaDriftDeg(bestCase, runnerUpCase);
end
end

function [selectedReplayIdx, selectionReason] = localResolvePeriodicReplayWinner(seedCandidateList, periodicDoaSeed, flowOpt)
selectedReplayIdx = localGetFieldOrDefault(periodicDoaSeed, 'bestFrozenIdx', 1);
selectionReason = "best-frozen";
staticIdx = localFindSeedSourceIndex(seedCandidateList, "static-seed");
if isnan(staticIdx)
  selectionReason = "no-static-frozen-candidate";
  return;
end
if selectedReplayIdx == staticIdx
  selectionReason = "best-frozen-static";
  return;
end
bestSeedSource = string(seedCandidateList(selectedReplayIdx).seedSource);
if bestSeedSource ~= "selected-subset"
  selectionReason = "best-frozen-nonstatic";
  return;
end
subsetDriftDeg = localGetFieldOrDefault(periodicDoaSeed, 'subsetDriftFromStaticDeg', NaN);
subsetMargin = localGetFieldOrDefault(periodicDoaSeed, 'subsetRankMarginRelative', NaN);
frozenDoaGap = localGetFieldOrDefault(periodicDoaSeed, 'frozenDoaDisagreementDeg', NaN);
frozenObjGap = localGetFieldOrDefault(periodicDoaSeed, 'frozenRelativeObjGap', NaN);
if localGetFieldOrDefault(flowOpt, 'periodicRefinePreferStaticWhenSubsetDriftLarge', true) && ...
    isfinite(subsetDriftDeg) && subsetDriftDeg > flowOpt.periodicRefineMaxSubsetDoaDriftDeg
  selectedReplayIdx = staticIdx;
  selectionReason = "subset-static-drift-too-large";
  return;
end
if isfinite(subsetMargin) && subsetMargin < flowOpt.periodicRefineSubsetTrustMinRelativeMargin
  selectedReplayIdx = staticIdx;
  selectionReason = "subset-rank-margin-too-small";
  return;
end
if isfinite(frozenDoaGap) && (frozenDoaGap > flowOpt.periodicRefineMaxFrozenDoaDisagreementDeg) && ...
    isfinite(frozenObjGap) && (frozenObjGap < flowOpt.periodicRefineMaxFrozenRelativeObjGapForTie)
  selectedReplayIdx = staticIdx;
  selectionReason = "frozen-subset-static-tie-with-large-doa-gap";
  return;
end
selectionReason = "best-frozen-subset-trusted";
end

function idx = localFindSeedSourceIndex(seedCandidateList, seedSource)
idx = NaN;
for iCand = 1:numel(seedCandidateList)
  if string(seedCandidateList(iCand).seedSource) == string(seedSource)
    idx = iCand;
    return;
  end
end
end

function tf = localIsScoreBetter(scoreA, scoreB)
tf = false;
for iElem = 1:min(numel(scoreA), numel(scoreB))
  a = scoreA(iElem);
  b = scoreB(iElem);
  if ~isfinite(a) && ~isfinite(b)
    continue;
  end
  if ~isfinite(a)
    return;
  end
  if ~isfinite(b)
    tf = true;
    return;
  end
  if a < b
    tf = true;
    return;
  end
  if a > b
    return;
  end
end
end

function idx = localGetRunnerUpIndex(orderIdx)
idx = NaN;
if numel(orderIdx) >= 2
  idx = orderIdx(2);
end
end

function driftDeg = localCalcCaseDoaDriftDeg(caseA, caseB)
driftDeg = NaN;
doaA = localGetCaseDoa(caseA);
doaB = localGetCaseDoa(caseB);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  driftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  driftDeg = NaN;
end
end

function doaParam = localGetCaseDoa(caseUse)
doaParam = [];
if ~isstruct(caseUse) || ~isfield(caseUse, 'estResult') || isempty(caseUse.estResult)
  return;
end
if ~isfield(caseUse.estResult, 'doaParamEst') || isempty(caseUse.estResult.doaParamEst)
  return;
end
rawValue = reshape(caseUse.estResult.doaParamEst, [], 1);
if any(~isfinite(rawValue))
  return;
end
doaParam = rawValue;
end

function relGap = localCalcRelativeObjGap(objA, objB)
relGap = NaN;
if ~isfinite(objA) || ~isfinite(objB)
  return;
end
scale = max([1, abs(objA), abs(objB)]);
relGap = abs(objA - objB) / scale;
end

function bucket = localBuildPeriodicHealthBucket(summary)
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

function satMode = localInferSatMode(viewUse)
satMode = "multi";
if isstruct(viewUse) && isfield(viewUse, 'numSat')
  if isequal(viewUse.numSat, 1)
    satMode = "single";
  end
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
