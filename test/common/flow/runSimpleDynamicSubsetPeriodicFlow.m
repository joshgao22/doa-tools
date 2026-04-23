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
[subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, useParfor] = localEvaluateSubsetBank( ...
  displayName, satMode, subsetFixtureCell, staticSeedCase, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse, toothStepHz);
stageTiming.subsetSelectSec = toc(stageTimer);

[bestSubsetIdx, subsetTrusted, subsetSelectionReason] = selectBestDynamicSubsetSummary(subsetSummaryCell);
subsetOrderIdx = localRankSubsetSummaryCell(subsetSummaryCell);
selectedSubsetFixture = subsetFixtureCell{bestSubsetIdx};
selectedSubsetCase = subsetCaseCell{bestSubsetIdx};
selectedSubsetSummary = subsetSummaryCell{bestSubsetIdx};
subsetCandidateTable = localBuildSubsetCandidateTable(subsetSummaryCell, subsetOrderIdx);
subsetTrustDiag = localBuildSubsetTrustDiag(subsetSummaryCell, subsetOrderIdx, subsetTrusted, subsetSelectionReason);

stageTimer = tic;
fdRangePeriodicRefine = [selectedSubsetCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangePeriodicRefine = [selectedSubsetCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];
[finalCase, periodicDoaSeed, periodicCaseCell, periodicSummaryCell, periodicCandidateTable] = ...
  localRunPeriodicRefineCase(displayName, satMode, periodicFixture, staticSeedCase, selectedSubsetCase, ...
  subsetTrustDiag, toothStepHz, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, truthUse);
finalSummary = buildDynamicUnknownCaseSummary(finalCase, toothStepHz, truthUse);
stageTiming.periodicRefineSec = toc(stageTimer);

stageTiming.totalFlowSec = toc(bundleTimer);

flow = struct();
flow.displayName = displayName;
flow.satMode = satMode;
flow.caseSubsetCell = subsetCaseCell;
flow.subsetSummaryCell = subsetSummaryCell;
flow.subsetRunTimeSec = subsetRunTimeSec;
flow.bestSubsetIdx = bestSubsetIdx;
flow.selectedSubsetLabel = string(localGetFieldOrDefault(selectedSubsetFixture, 'subsetLabel', "subset" + string(bestSubsetIdx)));
flow.selectedSubsetOffsets = reshape(localGetFieldOrDefault(selectedSubsetFixture, 'subsetOffsetIdx', []), 1, []);
flow.selectedSubsetCase = selectedSubsetCase;
flow.selectedSubsetSummary = selectedSubsetSummary;
flow.subsetCandidateTable = subsetCandidateTable;
flow.subsetTrustDiag = subsetTrustDiag;
flow.subsetSelectionReason = string(localGetFieldOrDefault(subsetTrustDiag, 'selectionReason', "unknown"));
flow.selectedSubsetTrusted = logical(localGetFieldOrDefault(subsetTrustDiag, 'isTrusted', false));
flow.fdRangePeriodicRefine = fdRangePeriodicRefine;
flow.fdRateRangePeriodicRefine = fdRateRangePeriodicRefine;
flow.periodicDoaSeed = periodicDoaSeed;
flow.periodicCaseCell = periodicCaseCell;
flow.periodicSummaryCell = periodicSummaryCell;
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

function [caseUse, periodicDoaSeed, periodicCaseCell, periodicSummaryCell, periodicCandidateTable] = localRunPeriodicRefineCase( ...
  displayName, satMode, periodicFixture, staticSeedCase, selectedSubsetCase, subsetTrustDiag, toothStepHz, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, truthUse)
[viewUse, ~] = localResolveViewAndTruth(satMode, periodicFixture);

seedCandidateList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOpt, satMode, staticSeedCase, selectedSubsetCase);
numFrozenCandidate = numel(seedCandidateList);
periodicCaseCell = cell(numFrozenCandidate, 1);
periodicSummaryCell = cell(numFrozenCandidate, 1);
scoreMat = nan(numFrozenCandidate, 18);
for iCand = 1:numFrozenCandidate
  [periodicCaseCell{iCand}, periodicSummaryCell{iCand}] = localRunPeriodicReplayCandidate( ...
    displayName, viewUse, truthUse, periodicFixture, selectedSubsetCase, seedCandidateList(iCand), ...
    pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz);
  scoreMat(iCand, :) = buildDynamicSelectionScoreVector(periodicSummaryCell{iCand});
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
bestFrozenIdx = orderIdx(1);
runnerUpFrozenIdx = localGetRunnerUpIndex(orderIdx);
periodicDoaSeed = localBuildPeriodicSeedDiag(seedCandidateList, periodicCaseCell, periodicSummaryCell, ...
  bestFrozenIdx, runnerUpFrozenIdx, subsetTrustDiag);

[selectedReplayIdx, replaySelectionReason] = localResolvePeriodicReplayWinner(seedCandidateList, periodicDoaSeed, flowOpt);
selectedReplaySummary = periodicSummaryCell{selectedReplayIdx};
periodicDoaSeed.selectedReplayIdx = selectedReplayIdx;
periodicDoaSeed.selectedReplaySeedSource = string(seedCandidateList(selectedReplayIdx).seedSource);
periodicDoaSeed.selectedReplayTag = string(seedCandidateList(selectedReplayIdx).startTag);
periodicDoaSeed.selectedReplaySelectionReason = string(replaySelectionReason);
periodicDoaSeed.selectedReplaySummary = selectedReplaySummary;
periodicDoaSeed.selectedReplayToothIdx = localGetFieldOrDefault(selectedReplaySummary, 'toothIdx', NaN);
periodicDoaSeed.selectedReplayToothResidualHz = localGetFieldOrDefault(selectedReplaySummary, 'toothResidualHz', NaN);
periodicDoaSeed.selectedReplayHealthBucket = localBuildPeriodicHealthBucket(selectedReplaySummary);

caseUse = periodicCaseCell{selectedReplayIdx};
periodicDoaSeed.usedDoaPolish = false;
[runPolish, polishGateReason, polishGateDiag] = localShouldRunPeriodicDoaPolish( ...
  satMode, periodicDoaSeed, selectedReplaySummary, flowOpt);
periodicDoaSeed.polishGateReason = string(polishGateReason);
periodicDoaSeed.polishGateDiag = polishGateDiag;
if runPolish
  [polishCase, polishSummary] = localRunPeriodicPolishCandidate( ...
    displayName, viewUse, truthUse, periodicFixture, caseUse, ...
    periodicDoaSeed.selectedReplaySeedSource, periodicDoaSeed.selectedReplayTag, ...
    pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz);
  periodicCaseCell{end + 1, 1} = polishCase; %#ok<AGROW>
  periodicSummaryCell{end + 1, 1} = polishSummary; %#ok<AGROW>
  candidateScore = buildDynamicSelectionScoreVector(polishSummary);
  selectedScore = buildDynamicSelectionScoreVector(selectedReplaySummary);
  choosePolish = localIsScoreBetter(candidateScore, selectedScore);
  periodicDoaSeed.usedDoaPolish = true;
  periodicDoaSeed.polishSelected = logical(choosePolish);
  periodicDoaSeed.polishCandidateSummary = polishSummary;
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
periodicCandidateTable = localBuildPeriodicCandidateTable(periodicSummaryCell);
end

function [caseUse, summaryUse] = localRunPeriodicReplayCandidate(displayName, viewUse, truthUse, periodicFixture, selectedSubsetCase, ...
  seedCandidate, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz)
fdRangeUse = [selectedSubsetCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangeUse = [selectedSubsetCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedCandidate.doaInitParam(:);
dynOpt.initDoaHalfWidth = seedCandidate.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = logical(seedCandidate.freezeDoa);
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = seedCandidate.doaHalfWidthDeg(:);

initParamSeed = buildDynamicInitParamFromCase(selectedSubsetCase, false, selectedSubsetCase.estResult.fdRateEst);
if ~isempty(initParamSeed) && numel(initParamSeed) >= numel(seedCandidate.doaInitParam)
  initParamSeed(1:numel(seedCandidate.doaInitParam)) = seedCandidate.doaInitParam(:);
end
initCandidate = struct( ...
  'startTag', string(seedCandidate.startTag), ...
  'initParam', initParamSeed, ...
  'initDoaParam', seedCandidate.doaInitParam(:), ...
  'initDoaHalfWidth', seedCandidate.doaHalfWidthDeg(:), ...
  'freezeDoa', logical(seedCandidate.freezeDoa));

caseUse = runDynamicDoaDopplerCase(displayName, localInferSatMode(viewUse), viewUse, truthUse, pilotWave, ...
  carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, ...
  optVerbose, dynOpt, false, debugTruth, initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truthUse);
summaryUse.seedSource = string(seedCandidate.seedSource);
summaryUse.startTag = string(seedCandidate.startTag);
summaryUse.subsetDriftFromStaticDeg = seedCandidate.subsetDriftFromStaticDeg;
end

function [caseUse, summaryUse] = localRunPeriodicPolishCandidate(displayName, viewUse, truthUse, periodicFixture, selectedReplayCase, ...
  replaySeedSource, replayTag, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz)
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

caseUse = runDynamicDoaDopplerCase(displayName, localInferSatMode(viewUse), viewUse, truthUse, pilotWave, ...
  carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, ...
  optVerbose, dynOpt, false, debugTruth, initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truthUse);
summaryUse.seedSource = string(replaySeedSource) + "-polish";
summaryUse.startTag = "periodic-polish-" + string(replaySeedSource);
summaryUse.replayBaseTag = string(replayTag);
summaryUse.subsetDriftFromStaticDeg = NaN;
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

function [runPolish, gateReason, gateDiag] = localShouldRunPeriodicDoaPolish(satMode, periodicDoaSeed, selectedReplaySummary, flowOpt)
runPolish = false;
gateReason = "disabled";
gateDiag = struct();
gateDiag.selectedReplaySeedSource = string(localGetFieldOrDefault(periodicDoaSeed, 'selectedReplaySeedSource', "unknown"));
gateDiag.subsetDriftFromStaticDeg = localGetFieldOrDefault(periodicDoaSeed, 'subsetDriftFromStaticDeg', NaN);
gateDiag.subsetRankMarginRelative = localGetFieldOrDefault(periodicDoaSeed, 'subsetRankMarginRelative', NaN);
gateDiag.frozenDoaDisagreementDeg = localGetFieldOrDefault(periodicDoaSeed, 'frozenDoaDisagreementDeg', NaN);
gateDiag.frozenRelativeObjGap = localGetFieldOrDefault(periodicDoaSeed, 'frozenRelativeObjGap', NaN);
gateDiag.selectedReplayHealthBucket = localBuildPeriodicHealthBucket(selectedReplaySummary);
gateDiag.selectedReplayToothIdx = localGetFieldOrDefault(selectedReplaySummary, 'toothIdx', NaN);
gateDiag.selectedReplayToothResidualHz = abs(localGetFieldOrDefault(selectedReplaySummary, 'toothResidualHz', inf));
if satMode ~= "multi"
  gateReason = "single-mode";
  return;
end
if ~localGetFieldOrDefault(flowOpt, 'periodicRefineEnableVerySmallDoaPolish', false)
  gateReason = "polish-disabled";
  return;
end
if string(gateDiag.selectedReplaySeedSource) ~= "selected-subset"
  gateReason = "replay-not-from-subset";
  return;
end
if ~logical(localGetFieldOrDefault(selectedReplaySummary, 'isResolved', false))
  gateReason = "selected-replay-unresolved";
  return;
end
if ~(isfinite(gateDiag.selectedReplayToothIdx) && abs(gateDiag.selectedReplayToothIdx) == 0)
  gateReason = "selected-replay-not-central-tooth";
  return;
end
maxToothResidualHz = localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMaxSelectedToothResidualHz', 50);
if ~(isfinite(gateDiag.selectedReplayToothResidualHz) && gateDiag.selectedReplayToothResidualHz <= maxToothResidualHz)
  gateReason = "selected-replay-tooth-residual-too-large";
  return;
end
if localBuildPeriodicHealthBucket(selectedReplaySummary) > localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMaxHealthBucket', 0)
  gateReason = "fd-health-not-ready";
  return;
end
subsetDriftDeg = gateDiag.subsetDriftFromStaticDeg;
if isfinite(subsetDriftDeg) && subsetDriftDeg > flowOpt.periodicRefinePolishTriggerDoaDriftDeg
  gateReason = "subset-static-drift-too-large";
  return;
end
subsetMargin = gateDiag.subsetRankMarginRelative;
if ~(isfinite(subsetMargin) && subsetMargin >= flowOpt.periodicRefineSubsetTrustMinRelativeMargin)
  gateReason = "subset-rank-margin-too-small";
  return;
end
frozenDoaGap = gateDiag.frozenDoaDisagreementDeg;
minFrozenDoaGap = localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMinFrozenDoaDisagreementDeg', 5e-4);
if ~isfinite(frozenDoaGap)
  gateReason = "missing-frozen-doa-gap";
  return;
end
if frozenDoaGap < minFrozenDoaGap
  gateReason = "frozen-doa-gap-too-small";
  return;
end
if frozenDoaGap > flowOpt.periodicRefineMaxFrozenDoaDisagreementDeg
  gateReason = "frozen-doa-disagreement-too-large";
  return;
end
frozenObjGap = gateDiag.frozenRelativeObjGap;
if ~(isfinite(frozenObjGap) && frozenObjGap <= flowOpt.periodicRefinePolishMaxFrozenRelativeObjGap)
  gateReason = "frozen-objective-gap-too-large";
  return;
end
runPolish = true;
gateReason = "triggered";
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
end
rankVec = nan(numCase, 1);
rankVec(orderIdx) = (1:numCase).';
tableOut = table(rankVec, subsetLabel, toothIdx, healthBucket, fdRefEst, fdRateEst, finalObj, finalResidualNorm, runTimeMs, ...
  'VariableNames', {'rank', 'subsetLabel', 'toothIdx', 'healthBucket', 'fdRefEstHz', 'fdRateEstHzPerSec', 'finalObj', 'finalResidualNorm', 'runTimeMs'});
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

function [subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, useParfor] = localEvaluateSubsetBank( ...
  displayName, satMode, subsetFixtureCell, staticSeedCase, pilotWave, carrierFreq, ...
  sampleRate, optVerbose, flowOpt, debugTruth, truthUse, toothStepHz)

runSubsetCaseFn = @(fixtureUse, pilotWaveUse, carrierFreqUse, sampleRateUse, optVerboseUse, flowOptUse, fdRangeIgnored, fdRateRangeIgnored, subsetSeedInfoIgnored) ...
  localRunSubsetUnknownCase(displayName + "-subset", satMode, fixtureUse, staticSeedCase, ...
    pilotWaveUse, carrierFreqUse, sampleRateUse, optVerboseUse, flowOptUse, debugTruth);
buildSummaryFn = @(caseUse, truthLocal, toothStepLocal) ...
  buildDynamicUnknownCaseSummary(caseUse, toothStepLocal, truthLocal);
subsetSeedInfo = struct();
[subsetCaseCell, subsetSummaryCell, subsetRunTimeSec] = evaluateDynamicSubsetBank( ...
  subsetFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  [NaN, NaN], [NaN, NaN], subsetSeedInfo, truthUse, toothStepHz, ...
  runSubsetCaseFn, buildSummaryFn);
useParfor = localShouldUseSubsetParfor(flowOpt.parallelOpt, numel(subsetCaseCell));
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
