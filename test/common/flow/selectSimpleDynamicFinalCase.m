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
periodicDoaSeed.usedBasinEntryRescue = false;
periodicDoaSeed.basinEntryRescueTriggered = false;
periodicDoaSeed.basinEntryRescueSelected = false;
[runBasinEntryRescue, basinEntryGateReason, basinEntryGateDiag] = shouldRunSimpleSameToothBasinEntryRescue( ...
  satMode, selectedReplayDecisionSummary, subsetTrustDiag, flowOpt);
periodicDoaSeed.basinEntryGateReason = string(basinEntryGateReason);
periodicDoaSeed.basinEntryGateDiag = basinEntryGateDiag;
if runBasinEntryRescue
  [rescueCaseCell, rescueSummaryCell, rescueDecisionSummaryCell] = localRunBasinEntryRescueCandidates( ...
    displayName, satMode, periodicFixture, caseUse, seedCandidateList, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, debugTruth, toothStepHz, truthUse);
  periodicDoaSeed.usedBasinEntryRescue = ~isempty(rescueCaseCell);
  periodicDoaSeed.basinEntryRescueTriggered = ~isempty(rescueCaseCell);
  periodicDoaSeed.basinEntryCandidateCount = numel(rescueCaseCell);
  periodicDoaSeed.basinEntryCandidateSummaryCell = rescueSummaryCell;
  for iRescue = 1:numel(rescueCaseCell)
    periodicCaseCell{end + 1, 1} = rescueCaseCell{iRescue}; %#ok<AGROW>
    periodicSummaryCell{end + 1, 1} = rescueSummaryCell{iRescue}; %#ok<AGROW>
    periodicDecisionSummaryCell{end + 1, 1} = rescueDecisionSummaryCell{iRescue}; %#ok<AGROW>
    candidateScore = buildDynamicSelectionScoreVector(rescueDecisionSummaryCell{iRescue});
    selectedScore = buildDynamicSelectionScoreVector(selectedReplayDecisionSummary);
    if localIsScoreBetter(candidateScore, selectedScore)
      caseUse = rescueCaseCell{iRescue};
      selectedReplaySummary = rescueSummaryCell{iRescue};
      selectedReplayDecisionSummary = rescueDecisionSummaryCell{iRescue};
      selectedReplayIdx = numel(periodicCaseCell);
      periodicDoaSeed.basinEntryRescueSelected = true;
      periodicDoaSeed.selectedReplaySeedSource = string(localGetFieldOrDefault(selectedReplaySummary, 'seedSource', "same-tooth-rescue"));
      periodicDoaSeed.selectedReplayTag = string(localGetFieldOrDefault(selectedReplaySummary, 'startTag', "same-tooth-rescue"));
      periodicDoaSeed.selectedReplaySelectionReason = "same-tooth-basin-entry-better";
    end
  end
  if ~periodicDoaSeed.basinEntryRescueSelected
    periodicDoaSeed.selectedReplaySelectionReason = "frozen-better-than-basin-entry";
  end
else
  periodicDoaSeed.basinEntryCandidateCount = 0;
  periodicDoaSeed.basinEntryCandidateSummaryCell = cell(0, 1);
end
periodicDoaSeed.selectedReplayIdx = selectedReplayIdx;
periodicDoaSeed.selectedReplaySummary = selectedReplaySummary;
periodicDoaSeed.selectedReplayDecisionSummary = selectedReplayDecisionSummary;
periodicDoaSeed.selectedReplayHealthBucket = localBuildPeriodicHealthBucket(selectedReplayDecisionSummary);
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

function [rescueCaseCell, rescueSummaryCell, rescueDecisionSummaryCell] = localRunBasinEntryRescueCandidates(displayName, satMode, periodicFixture, selectedReplayCase, ...
  seedCandidateList, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz, truthUse)
%LOCALRUNBASINENTRYRESCUECANDIDATES Evaluate gated wide and single-MF basin-entry candidates.

rescueCaseCell = cell(0, 1);
rescueSummaryCell = cell(0, 1);
rescueDecisionSummaryCell = cell(0, 1);
if satMode ~= "multi"
  return;
end

bankMode = lower(strtrim(char(string(localGetFieldOrDefault(flowOpt, 'sameToothRescueBankMode', "wide-single")))));
runWide = any(strcmp(bankMode, {'wide', 'wide-only', 'wide-single', 'wide-single-bank'}));
runSingle = any(strcmp(bankMode, {'single', 'single-mf', 'single-mf-only', 'single-only', 'wide-single', 'wide-single-bank'}));
staticDoa = localResolveSeedDoa(seedCandidateList, "static-seed", selectedReplayCase.estResult.doaParamEst(:));

if runWide
  [caseWide, summaryWide, decisionWide] = localRunBasinEntryDynamicCandidate( ...
    displayName, periodicFixture, selectedReplayCase, staticDoa, flowOpt.sameToothRescueWideDoaHalfWidthDeg(:), ...
    "same-tooth-rescue-wide", "wide-center", periodicFixture.viewMs, "multi", pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, debugTruth, toothStepHz, truthUse);
  rescueCaseCell{end + 1, 1} = caseWide; %#ok<AGROW>
  rescueSummaryCell{end + 1, 1} = summaryWide; %#ok<AGROW>
  rescueDecisionSummaryCell{end + 1, 1} = decisionWide; %#ok<AGROW>
end

if runSingle
  [caseSingleCenter, ~, ~] = localRunBasinEntryDynamicCandidate( ...
    displayName, periodicFixture, selectedReplayCase, staticDoa, flowOpt.sameToothRescueSingleDoaHalfWidthDeg(:), ...
    "same-tooth-rescue-single-center", "single-mf-center", periodicFixture.viewRefOnly, "single", pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, debugTruth, toothStepHz, truthUse);
  singleCenterDoa = selectedReplayCase.estResult.doaParamEst(:);
  if isfield(caseSingleCenter, 'estResult') && isfield(caseSingleCenter.estResult, 'doaParamEst') && ...
      ~isempty(caseSingleCenter.estResult.doaParamEst) && all(isfinite(caseSingleCenter.estResult.doaParamEst(:)))
    singleCenterDoa = caseSingleCenter.estResult.doaParamEst(:);
  end
  [caseSingle, summarySingle, decisionSingle] = localRunBasinEntryDynamicCandidate( ...
    displayName, periodicFixture, selectedReplayCase, singleCenterDoa, flowOpt.sameToothRescueMultiDoaHalfWidthDeg(:), ...
    "same-tooth-rescue-single-mf", "single-mf-center", periodicFixture.viewMs, "multi", pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, debugTruth, toothStepHz, truthUse);
  rescueCaseCell{end + 1, 1} = caseSingle; %#ok<AGROW>
  rescueSummaryCell{end + 1, 1} = summarySingle; %#ok<AGROW>
  rescueDecisionSummaryCell{end + 1, 1} = decisionSingle; %#ok<AGROW>
end
end

function [caseUse, summaryUse, decisionSummary] = localRunBasinEntryDynamicCandidate(displayName, periodicFixture, selectedReplayCase, ...
  seedDoaParam, doaHalfWidthDeg, startTag, seedSource, viewUse, satModeUse, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, debugTruth, toothStepHz, truthUse)
%LOCALRUNBASINENTRYDYNAMICCANDIDATE Run one same-tooth rescue branch.

fdRangeUse = [selectedReplayCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedReplayCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangeUse = [selectedReplayCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedReplayCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = false;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = doaHalfWidthDeg(:);

initParamSeed = buildDynamicInitParamFromCase(selectedReplayCase, false, selectedReplayCase.estResult.fdRateEst);
if ~isempty(initParamSeed) && numel(initParamSeed) >= numel(seedDoaParam)
  initParamSeed(1:numel(seedDoaParam)) = seedDoaParam(:);
end
initCandidate = struct( ...
  'startTag', string(startTag), ...
  'initParam', initParamSeed, ...
  'initDoaParam', seedDoaParam(:), ...
  'initDoaHalfWidth', doaHalfWidthDeg(:), ...
  'freezeDoa', false);

caseUse = runDynamicDoaDopplerCase(displayName, satModeUse, viewUse, truthUse, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, false, debugTruth, initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truthUse);
decisionSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, struct());
summaryUse.seedSource = string(seedSource);
summaryUse.startTag = string(startTag);
summaryUse.subsetDriftFromStaticDeg = NaN;
decisionSummary.seedSource = summaryUse.seedSource;
decisionSummary.startTag = summaryUse.startTag;
decisionSummary.subsetDriftFromStaticDeg = NaN;
end

function doaParam = localResolveSeedDoa(seedCandidateList, seedSource, fallbackDoa)
%LOCALRESOLVESEEDDOA Resolve a named periodic seed or fall back to the selected replay DoA.

doaParam = fallbackDoa(:);
for iCand = 1:numel(seedCandidateList)
  if string(seedCandidateList(iCand).seedSource) ~= string(seedSource)
    continue;
  end
  rawDoa = seedCandidateList(iCand).doaInitParam;
  if ~isempty(rawDoa) && all(isfinite(rawDoa(:)))
    doaParam = rawDoa(:);
    return;
  end
end
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
