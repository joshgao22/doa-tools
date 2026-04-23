function flow = runMsToothSelectFlow(modeTag, periodicFixture, subsetFixtureInput, caseWide, truth, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, weightSweepAlpha, ...
  doaOnlyOpt, staticBaseOpt, dynBaseOpt, staticMsOpt, subsetDoaHalfWidthDeg, ...
  inToothFdHalfWidthHz, inToothFdRateHalfWidthHzPerSec, inToothDoaHalfWidthDeg, ...
  toothStepHz, parallelOpt)
%RUNMSTOOTHSELECTFLOW Run one reusable MS-MF tooth-selection flow.
%
% The dev flow is:
%   1) Keep one caller-provided periodic wide result as the baseline.
%   2) Run one bank of anti-periodic subset solves to pick one reliable
%      fdRef tooth / branch.
%   3) Return to the periodic full-data window and refine inside one narrow
%      in-tooth fdRef box while keeping DoA nearly frozen.
%
% modeTag:
%   "known"   -> MS-MF-CP-K on every subset and on the final in-tooth replay.
%   "unknown" -> MS-MF-CP-U on every subset and on the final in-tooth replay.

arguments
  modeTag (1, 1) string
  periodicFixture (1, 1) struct
  subsetFixtureInput
  caseWide (1, 1) struct
  truth (1, 1) struct
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical
  weightSweepAlpha (:, 1) double
  doaOnlyOpt (1, 1) struct
  staticBaseOpt (1, 1) struct
  dynBaseOpt (1, 1) struct
  staticMsOpt (1, 1) struct
  subsetDoaHalfWidthDeg (:, 1) double
  inToothFdHalfWidthHz (1, 1) double
  inToothFdRateHalfWidthHzPerSec (1, 1) double
  inToothDoaHalfWidthDeg (:, 1) double
  toothStepHz (1, 1) double
  parallelOpt (1, 1) struct
end

modeTag = lower(string(modeTag));
if modeTag ~= "known" && modeTag ~= "unknown"
  error('runMsToothSelectFlow:InvalidMode', ...
    'modeTag must be "known" or "unknown".');
end

subsetFixtureCell = localNormalizeSubsetFixtureBank(subsetFixtureInput);
numSubset = numel(subsetFixtureCell);
subsetSummaryCell = cell(numSubset, 1);
subsetCaseCell = cell(numSubset, 1);
useParfor = localShouldUseParfor(parallelOpt, numSubset, 'minSubsetBankForParfor');
if useParfor
  parfor iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    subsetCaseCell{iSubset} = localRunMsCaseFromFixture(modeTag, ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
      subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, subsetDoaHalfWidthDeg(:));
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
else
  for iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    subsetCaseCell{iSubset} = localRunMsCaseFromFixture(modeTag, ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
      subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, subsetDoaHalfWidthDeg(:));
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
end

[bestSubsetIdx, ~] = selectBestDynamicSubsetSummary(subsetSummaryCell);
bestSubsetIdx = localMaybePromoteCentralRandomSubset(subsetFixtureCell, subsetSummaryCell, bestSubsetIdx);
bestSubsetCase = subsetCaseCell{bestSubsetIdx};
bestSubsetSummary = subsetSummaryCell{bestSubsetIdx};
fdWideSummary = localBuildCaseSummary(caseWide, truth, toothStepHz);

wideRefineCase = localRunMsCaseFromFixture(modeTag, ...
  periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
  weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
  [caseWide.estResult.fdRefEst - inToothFdHalfWidthHz, caseWide.estResult.fdRefEst + inToothFdHalfWidthHz], ...
  localBuildFdRateRange(caseWide, modeTag, inToothFdRateHalfWidthHzPerSec), ...
  inToothDoaHalfWidthDeg(:), caseWide.estResult.doaParamEst(:), caseWide, ...
  struct('disableUnknownDoaReleaseFloor', true, ...
    'unknownDoaReleaseHalfWidth', inToothDoaHalfWidthDeg(:)));
wideSummary = localBuildCaseSummary(wideRefineCase, truth, toothStepHz);

fdRangeSubsetAnchor = [bestSubsetCase.estResult.fdRefEst - inToothFdHalfWidthHz, ...
  bestSubsetCase.estResult.fdRefEst + inToothFdHalfWidthHz];
fdRangeInTooth = fdRangeSubsetAnchor;
if modeTag == "unknown"
  fdRateRangeSubsetAnchor = [bestSubsetCase.estResult.fdRateEst - inToothFdRateHalfWidthHzPerSec, ...
    bestSubsetCase.estResult.fdRateEst + inToothFdRateHalfWidthHzPerSec];
  fdRateRangeInTooth = fdRateRangeSubsetAnchor;
else
  fdRateRangeSubsetAnchor = [];
  fdRateRangeInTooth = [];
end

anchorDoaHalfWidth = repmat(1e-8, size(inToothDoaHalfWidthDeg(:)));
anchorCase = localRunMsCaseFromFixture(modeTag, ...
  periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
  weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
  fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, anchorDoaHalfWidth, ...
  bestSubsetCase.estResult.doaParamEst(:), bestSubsetCase, struct('freezeDoa', true));
anchorSummary = localBuildCaseSummary(anchorCase, truth, toothStepHz);

if localShouldRunInToothRefine(bestSubsetSummary)
  refineCase = localRunMsCaseFromFixture(modeTag, ...
    periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
    weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
    localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
    fdRangeInTooth, fdRateRangeInTooth, inToothDoaHalfWidthDeg(:), ...
    bestSubsetCase.estResult.doaParamEst(:), bestSubsetCase, ...
    struct('disableUnknownDoaReleaseFloor', true, ...
      'unknownDoaReleaseHalfWidth', inToothDoaHalfWidthDeg(:)));
  refineSummary = localBuildCaseSummary(refineCase, truth, toothStepHz);
else
  refineCase = anchorCase;
  refineSummary = anchorSummary;
end

selectedSummaryCell = {wideSummary; fdWideSummary; anchorSummary; refineSummary};
selectedTagList = ["periodic-wide"; "periodic-fd-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
finalSelectOpt = struct( ...
  'enableSubsetAnchorGuard', true, ...
  'subsetAnchorIdx', 3, ...
  'subsetAnchorTag', "periodic-fd-anchor", ...
  'maxFdRefDriftHz', inToothFdHalfWidthHz, ...
  'maxFdRateDriftHzPerSec', max(inToothFdRateHalfWidthHzPerSec, 0), ...
  'maxDoaDriftDeg', 0.01, ...
  'minObjGainToLeaveAnchor', 1e4, ...
  'minResidualGainToLeaveAnchor', 1e3);
[selectedFinalIdx, selectedFinalTag, finalSelectReason, finalScoreMat] = selectFinalDynamicUnknownSummary( ...
  selectedSummaryCell, selectedTagList, finalSelectOpt);
caseFinalCell = {wideRefineCase; caseWide; anchorCase; refineCase};
finalCase = caseFinalCell{selectedFinalIdx};

toothFlowTable = localBuildToothFlowSummaryTable(fdWideSummary, wideSummary, bestSubsetSummary, anchorSummary, refineSummary);
subsetCandidateTable = localBuildSubsetCandidateSummaryTable(subsetSummaryCell);
finalSelectTable = localBuildFinalSelectTable(selectedSummaryCell, selectedTagList, selectedFinalTag, finalScoreMat);

flow = struct();
flow.modeTag = modeTag;
flow.caseFdWide = caseWide;
flow.caseWide = wideRefineCase;
flow.caseSubset = bestSubsetCase;
flow.caseAnchor = anchorCase;
flow.caseFinal = finalCase;
flow.fdWideSummary = fdWideSummary;
flow.wideSummary = wideSummary;
flow.subsetSummary = bestSubsetSummary;
flow.anchorSummary = anchorSummary;
flow.refineSummary = refineSummary;
flow.toothFlowTable = toothFlowTable;
flow.subsetCandidateTable = subsetCandidateTable;
flow.finalSelectTable = finalSelectTable;
flow.bestSubsetIdx = bestSubsetIdx;
flow.selectedSubsetLabel = bestSubsetSummary.subsetLabel;
flow.selectedSubsetOffsets = bestSubsetSummary.subsetOffsetIdx;
flow.selectedFinalTag = selectedFinalTag;
flow.finalSelectReason = finalSelectReason;
flow.fdRangeSubsetAnchor = fdRangeSubsetAnchor;
flow.fdRateRangeSubsetAnchor = fdRateRangeSubsetAnchor;
flow.fdRangeInTooth = fdRangeInTooth;
flow.fdRateRangeInTooth = fdRateRangeInTooth;
flow.usedParfor = useParfor;
end

function bestIdx = localMaybePromoteCentralRandomSubset(subsetFixtureCell, subsetSummaryCell, selectedIdx)
%LOCALMAYBEPROMOTECENTRALRANDOMSUBSET Prefer one central random rescue when the current winner is weak.

bestIdx = selectedIdx;
if isempty(subsetSummaryCell) || selectedIdx < 1 || selectedIdx > numel(subsetSummaryCell)
  return;
end
selectedLabel = lower(strtrim(char(string(localGetFieldOrDefault( ...
  subsetFixtureCell{selectedIdx}, 'subsetLabel', "")))));
if startsWith(selectedLabel, 'random') || startsWith(selectedLabel, 'rescue-random')
  return;
end
selectedSummary = subsetSummaryCell{selectedIdx};
selectedResidualHz = abs(localGetFieldOrDefault(selectedSummary, 'toothResidualHz', inf));
selectedHealthBucket = localBuildSubsetHealthBucket(selectedSummary);
selectedNeedsRescue = localSummaryNeedsTrustedCentralRescue(selectedSummary, 20, inf);
randomIdxList = [];
randomSummaryCell = {};
for iCase = 1:numel(subsetFixtureCell)
  labelUse = lower(strtrim(char(string(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetLabel', "")))));
  if ~(startsWith(labelUse, 'random') || startsWith(labelUse, 'rescue-random'))
    continue;
  end
  summaryUse = subsetSummaryCell{iCase};
  if ~logical(localGetFieldOrDefault(summaryUse, 'isResolved', false))
    continue;
  end
  toothIdxUse = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
  toothResidualHzUse = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
  if ~(isfinite(toothIdxUse) && abs(toothIdxUse) == 0 && isfinite(toothResidualHzUse) && toothResidualHzUse <= 20)
    continue;
  end
  randomIdxList(end + 1, 1) = iCase; %#ok<AGROW>
  randomSummaryCell{end + 1, 1} = summaryUse; %#ok<AGROW>
end
if isempty(randomIdxList)
  return;
end
[bestRandomLocalIdx, ~] = selectBestDynamicSubsetSummary(randomSummaryCell);
bestRandomSummary = randomSummaryCell{bestRandomLocalIdx};
bestRandomResidualHz = abs(localGetFieldOrDefault(bestRandomSummary, 'toothResidualHz', inf));
bestRandomHealthBucket = localBuildSubsetHealthBucket(bestRandomSummary);
canUseHealthOverride = isfinite(selectedResidualHz) && isfinite(bestRandomResidualHz) && ...
  bestRandomResidualHz <= selectedResidualHz + 20 && bestRandomHealthBucket + 1 < selectedHealthBucket;
if ~(selectedNeedsRescue || canUseHealthOverride)
  return;
end
bestIdx = randomIdxList(bestRandomLocalIdx);
end


function tf = localSummaryNeedsTrustedCentralRescue(summaryUse, residualTolHz, maxHealthBucket)
%LOCALSUMMARYNEEDSTRUSTEDCENTRALRESCUE Return true when one summary still looks unsafe.

tf = true;
if ~(isstruct(summaryUse) && ~isempty(summaryUse) && ...
    logical(localGetFieldOrDefault(summaryUse, 'isResolved', false)))
  return;
end
toothIdx = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
if ~(isfinite(toothIdx) && abs(toothIdx) == 0)
  return;
end
toothResidualHz = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
if ~(isfinite(toothResidualHz) && toothResidualHz <= residualTolHz)
  return;
end
if nargin >= 3 && isfinite(maxHealthBucket)
  if localBuildSubsetHealthBucket(summaryUse) > maxHealthBucket
    return;
  end
end
tf = false;
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


function tf = localShouldRunInToothRefine(bestSubsetSummary)
%LOCALSHOULDRUNINTOOTHREFINE Skip redundant in-tooth replay on the central tooth.
% When the best subset already lands on tooth 0 with only a tiny residual,
% the nearly frozen in-tooth replay rarely changes the final choice but still
% costs a full extra solve.

tf = true;
subsetLabel = lower(strtrim(char(string(localGetFieldOrDefault(bestSubsetSummary, 'subsetLabel', "")))));
if startsWith(subsetLabel, 'random') || startsWith(subsetLabel, 'rescue-random')
  return;
end
toothIdx = localGetFieldOrDefault(bestSubsetSummary, 'toothIdx', NaN);
if ~isfinite(toothIdx) || abs(toothIdx) > 0
  return;
end
residualHz = abs(localGetFieldOrDefault(bestSubsetSummary, 'toothResidualHz', inf));
if residualHz <= 5
  tf = false;
end
end


function caseResult = localRunMsCaseFromFixture(modeTag, fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, weightSweepAlpha, doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  staticInitDoaHalfWidth, fdRangeUse, fdRateRangeUse, initDoaHalfWidth, initDoaParam, knownSeedCase, dynOptOverride)
%LOCALRUNMSCASEFROMFIXTURE Build one fixture-local MS-MF CP-K or CP-U case.

if nargin < 15 || isempty(initDoaHalfWidth)
  initDoaHalfWidth = [0.01; 0.01];
end
if nargin < 16
  initDoaParam = [];
end
if nargin < 17
  knownSeedCase = struct();
end
if nargin < 18
  dynOptOverride = struct();
end

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  fixture.wavelen, pilotWave, carrierFreq, sampleRate, fdRangeUse, ...
  fixture.truth, fixture.otherSatIdxGlobal, optVerbose, doaOnlyOpt, ...
  staticBaseOpt, weightSweepAlpha, staticInitDoaHalfWidth);
bestStaticMsCase = caseBundle.bestStaticMsCase;
truthDebug = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());

if isempty(initDoaParam)
  initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
end

if modeTag == "known"
  dynMsKnownOpt = dynBaseOpt;
  dynMsKnownOpt.initDoaParam = initDoaParam(:);
  dynMsKnownOpt.initDoaHalfWidth = initDoaHalfWidth(:);
  dynMsKnownOpt.enableFdAliasUnwrap = true;
  dynMsKnownOpt.continuousPhaseConsistencyWeight = ...
    localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseConsistencyWeight', 0.05);
  dynMsKnownOpt.continuousPhaseCollapsePenaltyWeight = ...
    localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseCollapsePenaltyWeight', 0);
  dynMsKnownOpt.continuousPhaseNegativeProjectionPenaltyWeight = ...
    localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNegativeProjectionPenaltyWeight', 0);
  dynMsKnownOpt.continuousPhaseNonRefFitFloorWeight = ...
    localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNonRefFitFloorWeight', 0);
  dynMsKnownOpt.continuousPhaseNonRefSupportFloorWeight = ...
    localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNonRefSupportFloorWeight', 0);
  dynMsKnownOpt.unknownWarmAnchorUseScaledSolve = ...
    localGetFieldOrDefault(dynBaseOpt, 'unknownWarmAnchorUseScaledSolve', true);
  dynMsKnownOpt.unknownWarmAnchorFallbackSqp = ...
    localGetFieldOrDefault(dynBaseOpt, 'unknownWarmAnchorFallbackSqp', true);
  dynMsKnownOpt = localApplyDynOptOverride(dynMsKnownOpt, dynOptOverride);
  caseResult = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
    fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
    fdRangeUse, [], optVerbose, dynMsKnownOpt, true, ...
    truthDebug, buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));
  return;
end

% unknown-rate: always seed from one known-rate solve built on the same fixture.
if isempty(knownSeedCase)
  knownSeedCase = localRunMsCaseFromFixture("known", fixture, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, weightSweepAlpha, doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
    staticInitDoaHalfWidth, fdRangeUse, [], initDoaHalfWidth, initDoaParam);
end

fdRateSeedKnown = localGetCaseFdRateEst(knownSeedCase, fixture.truth.fdRateFit);
initParamMsUnknownCpK = buildDynamicInitParamFromCase(knownSeedCase, false, fdRateSeedKnown);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = initDoaParam(:);
dynMsUnknownOpt.initDoaHalfWidth = initDoaHalfWidth(:);
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.continuousPhaseConsistencyWeight = ...
  localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseConsistencyWeight', 0.05);
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = ...
  localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = ...
  localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = ...
  localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = ...
  localGetFieldOrDefault(dynBaseOpt, 'continuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = ...
  localGetFieldOrDefault(dynBaseOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = ...
  localGetFieldOrDefault(dynBaseOpt, 'unknownWarmAnchorFallbackSqp', true);
dynMsUnknownOpt = localApplyDynOptOverride(dynMsUnknownOpt, dynOptOverride);
msUnknownCand = buildUnknownInitCandidateSet( ...
  knownSeedCase, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, caseBundle.staticMsOpt.initDoaHalfWidth);
caseResult = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, ...
  truthDebug, msUnknownCand);
end

function summary = localBuildCaseSummary(caseUse, truth, toothStepHz)
%LOCALBUILDCASESUMMARY Build one compact tooth-selection summary.

summary = struct();
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.doaParamEst = reshape(caseUse.estResult.doaParamEst(:), 1, []);
summary.fdRefEst = caseUse.estResult.fdRefEst;
summary.fdRateEst = caseUse.estResult.fdRateEst;
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseUse.estResult, 'aux', struct()), 'debug', struct());
finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
summary.finalObj = localGetFieldOrDefault(finalEval, 'obj', NaN);
summary.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
summary.refFitRatio = localGetFieldOrDefault(finalEval, 'refFitRatio', NaN);
summary.nonRefFitRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefFitRatioFloor', NaN);
summary.refSupportRatio = localGetFieldOrDefault(finalEval, 'refSupportRatio', NaN);
summary.nonRefSupportRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefSupportRatioFloor', NaN);
summary.refConsistencyNorm = localGetFieldOrDefault(finalEval, 'refConsistencyNorm', NaN);
summary.nonRefConsistencyRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefConsistencyRatioFloor', NaN);
summary.maxNonRefNegativeProjectionRatio = localGetFieldOrDefault(finalEval, 'maxNonRefNegativeProjectionRatio', NaN);
summary.isDoaFrozenLike = localIsFrozenLikeSummary(summary.solveVariant);
[summary.refCoherence, summary.nonRefCoherenceFloor, summary.nonRefMaxAbsPhaseResidRad, ...
  summary.nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval);
summary.subsetLabel = string.empty(1, 0);
summary.subsetOffsetIdx = [];
end

function toothTable = localBuildToothFlowSummaryTable(fdWideSummary, wideSummary, subsetSummary, anchorSummary, refineSummary)
%LOCALBUILDTOOTHFLOWSUMMARYTABLE Build one compact anti-periodic flow table.

stageName = ["periodic-fd-wide"; "periodic-wide"; "subset-select"; "periodic-fd-anchor"; "periodic-in-tooth"];
summaryCell = {fdWideSummary; wideSummary; subsetSummary; anchorSummary; refineSummary};
numStage = numel(summaryCell);
solveVariant = strings(numStage, 1);
isResolved = false(numStage, 1);
angleErrDeg = nan(numStage, 1);
fdRefErrHz = nan(numStage, 1);
fdRateErrHzPerSec = nan(numStage, 1);
toothIdx = nan(numStage, 1);
toothResidualHz = nan(numStage, 1);
runTimeMs = nan(numStage, 1);
for iStage = 1:numStage
  s = summaryCell{iStage};
  solveVariant(iStage) = s.solveVariant;
  isResolved(iStage) = s.isResolved;
  angleErrDeg(iStage) = s.angleErrDeg;
  fdRefErrHz(iStage) = s.fdRefErrHz;
  fdRateErrHzPerSec(iStage) = s.fdRateErrHzPerSec;
  toothIdx(iStage) = s.toothIdx;
  toothResidualHz(iStage) = s.toothResidualHz;
  runTimeMs(iStage) = s.runTimeMs;
end

toothTable = table(stageName, solveVariant, isResolved, angleErrDeg, ...
  fdRefErrHz, fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs);
end

function subsetTable = localBuildSubsetCandidateSummaryTable(summaryCell)
%LOCALBUILDSUBSETCANDIDATESUMMARYTABLE Build one compact subset-bank table.

numCase = numel(summaryCell);
subsetLabel = strings(numCase, 1);
subsetOffsetStr = strings(numCase, 1);
isResolved = false(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
for iCase = 1:numCase
  s = summaryCell{iCase};
  subsetLabel(iCase) = s.subsetLabel;
  subsetOffsetStr(iCase) = localFormatIntegerRow(s.subsetOffsetIdx);
  isResolved(iCase) = logical(s.isResolved);
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
end
subsetTable = table(subsetLabel, subsetOffsetStr, isResolved, angleErrDeg, ...
  fdRefErrHz, fdRateErrHzPerSec, toothIdx, toothResidualHz);
end


function finalSelectTable = localBuildFinalSelectTable(summaryCell, candidateTag, selectedTag, scoreMat)
%LOCALBUILDFINALSELECTTABLE Build one compact final-select table.

numCase = numel(summaryCell);
isSelected = (candidateTag == string(selectedTag));
isResolved = false(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
finalObj = nan(numCase, 1);
finalResidualNorm = nan(numCase, 1);
for iCase = 1:numCase
  s = summaryCell{iCase};
  isResolved(iCase) = s.isResolved;
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
  finalObj(iCase) = localGetFieldOrDefault(s, 'finalObj', NaN);
  finalResidualNorm(iCase) = localGetFieldOrDefault(s, 'finalResidualNorm', NaN);
end
baseVarName = {'candidateTag', 'isSelected', 'isResolved', 'angleErrDeg', 'fdRefErrHz', ...
  'fdRateErrHzPerSec', 'toothIdx', 'toothResidualHz', 'finalObj', 'finalResidualNorm'};
baseValue = {candidateTag, isSelected, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, finalObj, finalResidualNorm};
scoreCell = cell(1, size(scoreMat, 2));
scoreVarName = localBuildScoreVarNames(size(scoreMat, 2));
for iCol = 1:size(scoreMat, 2)
  scoreCell{iCol} = scoreMat(:, iCol);
end
finalSelectTable = table(baseValue{:}, scoreCell{:}, 'VariableNames', [baseVarName, scoreVarName]);
end


function scoreVarName = localBuildScoreVarNames(numScore)
%LOCALBUILDSCOREVARNAMES Build score-column names for the final-select table.

scoreVarName = cellstr("score" + string(1:numScore));
templateName = {'scoreResolvedPenalty', 'scoreNonRefSupportPenalty', 'scoreNonRefFitPenalty', ...
  'scoreNonRefConsistencyPenalty', 'scoreNonRefCoherencePenalty', ...
  'scoreNonRefRmsPhaseResidRad', 'scoreNonRefMaxAbsPhaseResidRad', ...
  'scoreNonRefNegativePenalty', 'scoreFinalObj', 'scoreFinalResidualNorm', 'scoreRunTimeMs'};
numTemplate = min(numScore, numel(templateName));
scoreVarName(1:numTemplate) = templateName(1:numTemplate);
end


function [refCoherence, nonRefCoherenceFloor, nonRefMaxAbsPhaseResidRad, nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval)
%LOCALEXTRACTNONREFEVALMETRICS Build compact truth-free non-reference metrics.

refCoherence = NaN;
nonRefCoherenceFloor = NaN;
nonRefMaxAbsPhaseResidRad = NaN;
nonRefRmsPhaseResidRad = NaN;
if ~isstruct(finalEval)
  return;
end
refSatIdxLocal = localGetFieldOrDefault(finalEval, 'refSatIdxLocal', NaN);
coherenceSat = reshape(localGetFieldOrDefault(finalEval, 'coherenceSat', []), [], 1);
if isfinite(refSatIdxLocal) && ~isempty(coherenceSat) && refSatIdxLocal >= 1 && refSatIdxLocal <= numel(coherenceSat)
  refCoherence = coherenceSat(refSatIdxLocal);
  nonRefMask = true(size(coherenceSat));
  nonRefMask(refSatIdxLocal) = false;
  nonRefVals = coherenceSat(nonRefMask & isfinite(coherenceSat));
  if ~isempty(nonRefVals)
    nonRefCoherenceFloor = min(nonRefVals);
  end
end
blockPhaseResidMat = localGetFieldOrDefault(finalEval, 'blockPhaseResidMat', []);
if isempty(blockPhaseResidMat)
  return;
end
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= size(blockPhaseResidMat, 1))
  return;
end
nonRefMask = true(size(blockPhaseResidMat, 1), 1);
nonRefMask(refSatIdxLocal) = false;
nonRefPhaseResid = abs(blockPhaseResidMat(nonRefMask, :));
nonRefPhaseResid = nonRefPhaseResid(isfinite(nonRefPhaseResid));
if isempty(nonRefPhaseResid)
  return;
end
nonRefMaxAbsPhaseResidRad = max(nonRefPhaseResid);
nonRefRmsPhaseResidRad = sqrt(mean(nonRefPhaseResid .^ 2));
end


function tf = localIsFrozenLikeSummary(solveVariant)
%LOCALISFROZENLIKESUMMARY Return true for fixed-DoA-style summary variants.

tf = contains(lower(strtrim(string(solveVariant))), "fixeddoa");
end


function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row for compact console printing.

if isempty(valueRow)
  textOut = "[]";
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '));
textOut = string(textOut);
end

function useParfor = localShouldUseParfor(parallelOpt, numCase, minFieldName)
%LOCALSHOULDUSEPARFOR Decide whether the current loop should use parfor.

useParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableSubsetSolveParfor', ...
  localGetFieldOrDefault(parallelOpt, 'enableParfor', true)));
if ~useParfor
  return;
end
minCase = localGetFieldOrDefault(parallelOpt, minFieldName, 12);
useParfor = useParfor && (numCase >= minCase) && localCanUseParfor();
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

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end

function fdRateEst = localGetCaseFdRateEst(caseUse, defaultValue)
%LOCALGETCASEFDRATEEST Safely read fdRateEst from one case struct.

if nargin < 2 || isempty(defaultValue)
  defaultValue = 0;
end
fdRateEst = defaultValue;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult)
  return;
end
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', defaultValue);
if isempty(fdRateEst) || ~isscalar(fdRateEst) || ~isfinite(fdRateEst)
  fdRateEst = defaultValue;
end
end
