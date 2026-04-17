function bundle = buildDoaDopplerDynamicTransitionBundle(periodicFixture, subsetFixtureCell, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt)
%BUILDDOADOPPLERDYNAMICTRANSITIONBUNDLE Build one shared static->dynamic MF bundle.
% This helper keeps the dynamic dev/perf scripts aligned with the static
% dual-satellite scripts: the entry script only prepares the scene/snapshot
% fixture and this helper runs the shared SF transition plus the current MF
% CP-K / CP-U flow.

arguments
  periodicFixture (1, 1) struct
  subsetFixtureCell cell
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical = false
  flowOpt (1, 1) struct = struct()
end

flowOpt = localApplyDynamicFlowDefaults(flowOpt);
truth = periodicFixture.truth;
viewRefOnly = periodicFixture.viewRefOnly;
viewOtherOnly = periodicFixture.viewOtherOnly;
viewMs = periodicFixture.viewMs;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, periodicFixture.wavelen, pilotWave, ...
  carrierFreq, sampleRate, periodicFixture.fdRange, truth, ...
  periodicFixture.otherSatIdxGlobal, optVerbose, flowOpt.doaOnlyOpt, ...
  flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));

caseRefDoa = caseBundle.caseRefDoa;
caseMsDoa = caseBundle.caseMsDoa;
caseStaticRefOnly = caseBundle.caseStaticRefOnly;
caseStaticRefAbl = caseBundle.caseStaticRefAbl;
caseStaticOtherOnly = caseBundle.caseStaticOtherOnly;
caseStaticMs = caseBundle.caseStaticMs;
weightCase = caseBundle.weightCase;
bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

initParamStaticRef = buildDynamicInitParamFromCase(caseStaticRefOnly, true, truth.fdRateFit);

periodicFixture.debugTruthRef = localBuildTruthDebugForView(viewRefOnly, truth);
periodicFixture.debugTruthMs = localBuildTruthDebugForView(viewMs, truth);
caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef);
caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly);

caseDynMsKnown = localRunMsKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, bestStaticMsCase);
caseDynMsUnknownWide = localRunMsUnknownWideCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynMsKnown, bestStaticMsCase, staticMsOpt);

subsetCaseCell = cell(numel(subsetFixtureCell), 1);
subsetSummaryCell = cell(numel(subsetFixtureCell), 1);
toothStepHz = 1 / median(diff(periodicFixture.sceneSeq.timeOffsetSec));
wideUnknownSummary = localBuildUnknownCaseSummary(caseDynMsUnknownWide, truth, toothStepHz);
for iSubset = 1:numel(subsetFixtureCell)
  fixtureUse = subsetFixtureCell{iSubset};
  fixtureUse.debugTruthMs = localBuildTruthDebugForView(fixtureUse.viewMs, fixtureUse.truth);
  subsetCaseCell{iSubset} = localRunMsUnknownCaseFromFixture(fixtureUse, pilotWave, ...
    carrierFreq, sampleRate, optVerbose, flowOpt, periodicFixture.fdRange, periodicFixture.fdRateRange);
  subsetSummaryCell{iSubset} = localBuildUnknownCaseSummary( ...
    subsetCaseCell{iSubset}, fixtureUse.truth, toothStepHz);
end

if isempty(subsetFixtureCell)
  bestSubsetIdx = NaN;
  subsetIsTrusted = false;
  selectedSubsetFixture = struct('subsetLabel', "none", 'subsetOffsetIdx', []);
  selectedSubsetSummary = wideUnknownSummary;
  subsetCandidateTable = table();
else
  [bestSubsetIdx, subsetIsTrusted] = localSelectBestSubsetForDev(subsetSummaryCell, toothStepHz);
  selectedSubsetFixture = subsetFixtureCell{bestSubsetIdx};
  selectedSubsetSummary = subsetSummaryCell{bestSubsetIdx};
  subsetCandidateTable = localBuildSubsetCandidateTable(subsetFixtureCell, subsetSummaryCell);
end

fdRangeInTooth = [caseDynMsUnknownWide.estResult.fdRefEst - flowOpt.inToothFdHalfWidthHz, ...
  caseDynMsUnknownWide.estResult.fdRefEst + flowOpt.inToothFdHalfWidthHz];
fdRateRangeInTooth = [caseDynMsUnknownWide.estResult.fdRateEst - flowOpt.inToothFdRateHalfWidthHzPerSec, ...
  caseDynMsUnknownWide.estResult.fdRateEst + flowOpt.inToothFdRateHalfWidthHzPerSec];

if subsetIsTrusted
  refinedCaseDynMsUnknown = localRunMsUnknownInTooth(periodicFixture, pilotWave, ...
    carrierFreq, sampleRate, optVerbose, flowOpt, ...
    fdRangeInTooth, fdRateRangeInTooth, subsetCaseCell{bestSubsetIdx});
  refinedUnknownSummary = localBuildUnknownCaseSummary(refinedCaseDynMsUnknown, truth, toothStepHz);
else
  refinedCaseDynMsUnknown = caseDynMsUnknownWide;
  refinedUnknownSummary = wideUnknownSummary;
end

[caseDynMsUnknown, selectedFinalTag] = localSelectFinalUnknownCaseForDev( ...
  caseDynMsUnknownWide, wideUnknownSummary, refinedCaseDynMsUnknown, refinedUnknownSummary);

toothFlowTable = localBuildToothFlowSummaryTable( ...
  wideUnknownSummary, selectedSubsetSummary, refinedUnknownSummary);

baseCase = [caseRefDoa, caseMsDoa, caseStaticRefOnly, caseStaticMs, ...
  caseDynRefKnown, caseDynRefUnknown, caseDynMsKnown, caseDynMsUnknown];
caseResult = [baseCase, weightCase];

dynSummary = summarizeDynamicEstimatorDebug(baseCase, periodicFixture.sceneSeq, truth);
dynMultiStartTable = summarizeDynamicMultiStart(baseCase, truth);

bundle = struct();
bundle.caseBundle = caseBundle;
bundle.caseRefDoa = caseRefDoa;
bundle.caseMsDoa = caseMsDoa;
bundle.caseStaticRefOnly = caseStaticRefOnly;
bundle.caseStaticRefAbl = caseStaticRefAbl;
bundle.caseStaticOtherOnly = caseStaticOtherOnly;
bundle.caseStaticMs = caseStaticMs;
bundle.weightCase = weightCase;
bundle.bestStaticMsCase = bestStaticMsCase;
bundle.caseDynRefKnown = caseDynRefKnown;
bundle.caseDynRefUnknown = caseDynRefUnknown;
bundle.caseDynMsKnown = caseDynMsKnown;
bundle.caseDynMsUnknownWide = caseDynMsUnknownWide;
bundle.caseDynMsUnknown = caseDynMsUnknown;
bundle.baseCase = baseCase;
bundle.caseResult = caseResult;
bundle.truth = truth;
bundle.toothStepHz = toothStepHz;
bundle.wideUnknownSummary = wideUnknownSummary;
bundle.selectedSubsetSummary = selectedSubsetSummary;
bundle.refinedUnknownSummary = refinedUnknownSummary;
bundle.selectedSubsetLabel = string(localGetFieldOrDefault(selectedSubsetFixture, 'subsetLabel', "none"));
bundle.selectedSubsetOffsetIdx = reshape(localGetFieldOrDefault(selectedSubsetFixture, 'subsetOffsetIdx', []), 1, []);
bundle.bestSubsetIdx = bestSubsetIdx;
bundle.subsetIsTrusted = subsetIsTrusted;
bundle.selectedFinalTag = string(selectedFinalTag);
bundle.fdRangeInTooth = fdRangeInTooth;
bundle.fdRateRangeInTooth = fdRateRangeInTooth;
bundle.toothFlowTable = toothFlowTable;
bundle.subsetCandidateTable = subsetCandidateTable;
bundle.dynSummary = dynSummary;
bundle.dynMultiStartTable = dynMultiStartTable;
end


function flowOpt = localApplyDynamicFlowDefaults(flowOpt)
%LOCALAPPLYDYNAMICFLOWDEFAULTS Fill one dynamic-flow option struct.

flowOpt = localMergeStruct(struct( ...
  'weightSweepAlpha', [0; 0.25; 0.5; 1], ...
  'staticMsHalfWidth', [0.002; 0.002], ...
  'doaOnlyOpt', struct('useLogObjective', true), ...
  'staticBaseOpt', struct('useLogObjective', true), ...
  'dynBaseOpt', struct( ...
    'useLogObjective', true, ...
    'initFdCount', 81, ...
    'useAccessMask', false, ...
    'phaseMode', 'continuous', ...
    'steeringMode', 'framewise', ...
    'debugEnable', true, ...
    'debugStoreEvalTrace', false, ...
    'debugMaxEvalTrace', 120), ...
  'refKnownDoaHalfWidth', [0.005; 0.005], ...
  'refUnknownDoaHalfWidth', [0.003; 0.003], ...
  'msKnownDoaHalfWidth', [0.003; 0.003], ...
  'msUnknownDoaHalfWidth', [0.002; 0.002], ...
  'subsetSelectDoaHalfWidthDeg', [0.01; 0.01], ...
  'inToothFdHalfWidthHz', 50, ...
  'inToothFdRateHalfWidthHzPerSec', 100, ...
  'inToothDoaHalfWidthDeg', [2e-4; 2e-4]), flowOpt);

flowOpt.weightSweepAlpha = reshape(flowOpt.weightSweepAlpha, [], 1);
flowOpt.staticMsHalfWidth = reshape(flowOpt.staticMsHalfWidth, [], 1);
flowOpt.refKnownDoaHalfWidth = reshape(flowOpt.refKnownDoaHalfWidth, [], 1);
flowOpt.refUnknownDoaHalfWidth = reshape(flowOpt.refUnknownDoaHalfWidth, [], 1);
flowOpt.msKnownDoaHalfWidth = reshape(flowOpt.msKnownDoaHalfWidth, [], 1);
flowOpt.msUnknownDoaHalfWidth = reshape(flowOpt.msUnknownDoaHalfWidth, [], 1);
flowOpt.subsetSelectDoaHalfWidthDeg = reshape(flowOpt.subsetSelectDoaHalfWidthDeg, [], 1);
flowOpt.inToothDoaHalfWidthDeg = reshape(flowOpt.inToothDoaHalfWidthDeg, [], 1);
end


function caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef)
%LOCALRUNREFKNOWNCASE Run one SS-MF-CP-K case.

dynRefKnownOpt = flowOpt.dynBaseOpt;
dynRefKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefKnownOpt.initDoaHalfWidth = flowOpt.refKnownDoaHalfWidth;
dynRefKnownOpt.initDoaParam = reshape(initParamStaticRef(1:2), [], 1);

caseDynRefKnown = runDynamicDoaDopplerCase("SS-MF-CP-K", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefKnownOpt, true, ...
  periodicFixture.debugTruthRef, initParamStaticRef);
end


function caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly)
%LOCALRUNREFUNKNOWNCASE Run one SS-MF-CP-U case.

dynRefUnknownOpt = flowOpt.dynBaseOpt;
dynRefUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefUnknownOpt.initDoaParam = caseStaticRefOnly.estResult.doaParamEst(:);
dynRefUnknownOpt.initDoaHalfWidth = flowOpt.refUnknownDoaHalfWidth;

initParamRefUnknownCpK = buildDynamicInitParamFromCase(caseDynRefKnown, false, caseDynRefKnown.estResult.fdRateEst);
initParamRefUnknownStatic = buildDynamicInitParamFromCase(caseStaticRefOnly, false, periodicFixture.truth.fdRateFit);
refUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynRefKnown, caseStaticRefOnly, initParamRefUnknownCpK, initParamRefUnknownStatic, ...
  dynRefUnknownOpt.initDoaHalfWidth, flowOpt.refKnownDoaHalfWidth);
caseDynRefUnknown = runDynamicDoaDopplerCase("SS-MF-CP-U", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefUnknownOpt, false, ...
  periodicFixture.debugTruthRef, refUnknownCand);
end


function caseDynMsKnown = localRunMsKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, bestStaticMsCase)
%LOCALRUNMSKNOWNCASE Run one MS-MF-CP-K case.

dynMsKnownOpt = flowOpt.dynBaseOpt;
dynMsKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth;
dynMsKnownOpt.enableFdAliasUnwrap = true;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynMsKnownOpt, true, ...
  periodicFixture.debugTruthMs, buildDynamicInitParamFromCase(bestStaticMsCase, true, periodicFixture.truth.fdRateFit));
end


function caseDynMsUnknownWide = localRunMsUnknownWideCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynMsKnown, bestStaticMsCase, staticMsOpt)
%LOCALRUNMSUNKNOWNWIDECASE Run the periodic wide-box MS-MF-CP-U case.

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.msUnknownDoaHalfWidth;
dynMsUnknownOpt.enableFdAliasUnwrap = true;

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, periodicFixture.truth.fdRateFit);
msUnknownCandWide = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknownWide = runDynamicDoaDopplerCase("MS-MF-CP-U-wide", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, msUnknownCandWide);
end


function caseDynMsUnknown = localRunMsUnknownCaseFromFixture(fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeUse, fdRateRangeUse)
%LOCALRUNMSUNKNOWNCASEFROMFIXTURE Run one MS-MF-CP-U case on one frame subset.

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, fixture.wavelen, ...
  pilotWave, carrierFreq, sampleRate, fdRangeUse, fixture.truth, ...
  fixture.otherSatIdxGlobal, optVerbose, flowOpt.doaOnlyOpt, ...
  flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

dynMsKnownOpt = flowOpt.dynBaseOpt;
dynMsKnownOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth;
dynMsKnownOpt.enableFdAliasUnwrap = true;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsKnownOpt, true, fixture.debugTruthMs, ...
  buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.subsetSelectDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, fixture.debugTruthMs, msUnknownCand);
end


function caseDynMsUnknown = localRunMsUnknownInTooth(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeInTooth, fdRateRangeInTooth, selectedSubsetCase)
%LOCALRUNMSUNKNOWNINTOOTH Run the periodic in-tooth MS-MF-CP-U refine.

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  periodicFixture.viewRefOnly, periodicFixture.viewOtherOnly, periodicFixture.viewMs, ...
  periodicFixture.wavelen, pilotWave, carrierFreq, sampleRate, fdRangeInTooth, ...
  periodicFixture.truth, periodicFixture.otherSatIdxGlobal, optVerbose, ...
  flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

dynMsKnownOpt = flowOpt.dynBaseOpt;
dynMsKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth;
dynMsKnownOpt.enableFdAliasUnwrap = true;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeInTooth, fdRateRangeInTooth, optVerbose, dynMsKnownOpt, true, periodicFixture.debugTruthMs, ...
  buildDynamicInitParamFromCase(bestStaticMsCase, true, periodicFixture.truth.fdRateFit));

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, periodicFixture.truth.fdRateFit);

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = selectedSubsetCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.inToothDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeInTooth, fdRateRangeInTooth, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, msUnknownCand);
end


function summary = localBuildUnknownCaseSummary(caseUse, truth, toothStepHz)
%LOCALBUILDUNKNOWNCASESUMMARY Build one compact tooth-selection summary.

summary = struct();
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
end


function toothTable = localBuildToothFlowSummaryTable(wideSummary, subsetSummary, refineSummary)
%LOCALBUILDTOOTHFLOWSUMMARYTABLE Build one compact anti-periodic flow table.

stageName = ["periodic-wide"; "subset-select"; "periodic-in-tooth"];
solveVariant = [wideSummary.solveVariant; subsetSummary.solveVariant; refineSummary.solveVariant];
isResolved = [wideSummary.isResolved; subsetSummary.isResolved; refineSummary.isResolved];
angleErrDeg = [wideSummary.angleErrDeg; subsetSummary.angleErrDeg; refineSummary.angleErrDeg];
fdRefErrHz = [wideSummary.fdRefErrHz; subsetSummary.fdRefErrHz; refineSummary.fdRefErrHz];
fdRateErrHzPerSec = [wideSummary.fdRateErrHzPerSec; subsetSummary.fdRateErrHzPerSec; refineSummary.fdRateErrHzPerSec];
toothIdx = [wideSummary.toothIdx; subsetSummary.toothIdx; refineSummary.toothIdx];
toothResidualHz = [wideSummary.toothResidualHz; subsetSummary.toothResidualHz; refineSummary.toothResidualHz];
runTimeMs = [wideSummary.runTimeMs; subsetSummary.runTimeMs; refineSummary.runTimeMs];

toothTable = table(stageName, solveVariant, isResolved, angleErrDeg, ...
  fdRefErrHz, fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs);
end


function [bestIdx, isTrusted] = localSelectBestSubsetForDev(summaryCell, toothStepHz)
%LOCALSELECTBESTSUBSETFORDEV Pick one subset using truth-aware dev ranking.

numCase = numel(summaryCell);
scoreMat = inf(numCase, 5);
for iCase = 1:numCase
  s = summaryCell{iCase};
  if ~s.isResolved
    continue;
  end
  scoreMat(iCase, 1) = abs(s.fdRefErrHz);
  scoreMat(iCase, 2) = abs(s.fdRateErrHzPerSec);
  scoreMat(iCase, 3) = s.angleErrDeg;
  scoreMat(iCase, 4) = abs(s.toothIdx);
  scoreMat(iCase, 5) = abs(s.toothResidualHz);
end
[~, orderIdx] = sortrows(scoreMat, [1, 2, 3, 4, 5]);
bestIdx = orderIdx(1);
bestSummary = summaryCell{bestIdx};
isTrusted = bestSummary.isResolved && ...
  abs(bestSummary.fdRefErrHz) <= max(100, 0.2 * toothStepHz) && ...
  abs(bestSummary.fdRateErrHzPerSec) <= 250;
end


function subsetTable = localBuildSubsetCandidateTable(subsetFixtureCell, summaryCell)
%LOCALBUILDSUBSETCANDIDATETABLE Build one compact subset ranking table.

numCase = numel(summaryCell);
labelOut = strings(numCase, 1);
offsetText = strings(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
isResolved = false(numCase, 1);
for iCase = 1:numCase
  fixtureUse = subsetFixtureCell{iCase};
  s = summaryCell{iCase};
  labelOut(iCase) = string(localGetFieldOrDefault(fixtureUse, 'subsetLabel', sprintf('cand%d', iCase)));
  offsetText(iCase) = string(localFormatIntegerRow(localGetFieldOrDefault(fixtureUse, 'subsetOffsetIdx', [])));
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
  isResolved(iCase) = s.isResolved;
end
subsetTable = table(labelOut, offsetText, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, ...
  'VariableNames', {'label', 'offsets', 'isResolved', 'angleErrDeg', ...
  'fdRefErrHz', 'fdRateErrHzPerSec', 'toothIdx', 'toothResidualHz'});
end


function [caseUse, selectedTag] = localSelectFinalUnknownCaseForDev(wideCase, wideSummary, refineCase, refineSummary)
%LOCALSELECTFINALUNKNOWNCASEFORDEV Keep the safer unknown result in the dev flow.

scoreWide = localBuildDevScoreVector(wideSummary);
scoreRefine = localBuildDevScoreVector(refineSummary);
isRefineBetter = false;
for iElem = 1:numel(scoreWide)
  if scoreRefine(iElem) < scoreWide(iElem)
    isRefineBetter = true;
    break;
  elseif scoreRefine(iElem) > scoreWide(iElem)
    break;
  end
end

if isRefineBetter
  caseUse = refineCase;
  selectedTag = "periodic-in-tooth";
else
  caseUse = wideCase;
  selectedTag = "periodic-wide";
end
end


function scoreVec = localBuildDevScoreVector(summary)
%LOCALBUILDDEVSCOREVECTOR Build one lexicographic score vector.

scoreVec = [double(~summary.isResolved), abs(summary.fdRefErrHz), ...
  abs(summary.fdRateErrHzPerSec), summary.angleErrDeg, abs(summary.toothResidualHz)];
end


function debugTruth = localBuildTruthDebugForView(view, truth)
%LOCALBUILDTRUTHDEBUGFORVIEW Build one debug-truth struct aligned with one view.

sceneSeqUse = localGetFieldOrDefault(view, 'sceneSeq', []);
if isempty(sceneSeqUse)
  sceneRefUse = view.sceneRef;
  debugTruth = struct();
  debugTruth.localDoa = reshape(sceneRefUse.localDoa, 2, sceneRefUse.numSat, 1);
  debugTruth.fdLocal = reshape(localBuildTruthFdLocalForView(sceneRefUse, truth), sceneRefUse.numSat, 1);
  debugTruth.deltaFdRef = reshape(truth.deltaFdTrueHz(:), sceneRefUse.numSat, 1);
  debugTruth.deltaFdRate = zeros(sceneRefUse.numSat, 1);
  debugTruth.refFrameIdx = 1;
  debugTruth.refSatIdxLocal = localResolveRefSatIdxLocalForView(sceneRefUse, truth);
  return;
end

sceneRefUse = sceneSeqUse.sceneCell{sceneSeqUse.refFrameIdx};
debugTruth = struct();
debugTruth.localDoa = localExtractSceneLocalDoa(sceneSeqUse);
debugTruth.fdLocal = localBuildTruthSeriesForView(sceneRefUse, truth, 'fdSatSeries', sceneSeqUse.numFrame, NaN);
debugTruth.deltaFdRef = localBuildTruthSeriesForView(sceneRefUse, truth, 'deltaFdSeries', sceneSeqUse.numFrame, NaN);
debugTruth.deltaFdRate = localBuildTruthSeriesForView(sceneRefUse, truth, 'deltaFdRate', sceneSeqUse.numFrame, 0);
debugTruth.refFrameIdx = sceneSeqUse.refFrameIdx;
debugTruth.refSatIdxLocal = localResolveRefSatIdxLocalForView(sceneRefUse, truth);
end


function seriesView = localBuildTruthSeriesForView(sceneRef, truth, fieldName, numFrame, defaultValue)
%LOCALBUILDTRUTHSERIESFORVIEW Align one truth series field to one view.

numSat = sceneRef.numSat;
seriesView = defaultValue * ones(numSat, numFrame);
seriesRaw = localNormalizeTruthSeriesField(localGetFieldOrDefault(truth, fieldName, []), numFrame);
selectedSatIdxGlobal = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
if isempty(seriesRaw) || isempty(selectedSatIdxGlobal) || isempty(sceneSatIdx)
  return;
end

for iSat = 1:min(numSat, numel(sceneSatIdx))
  matchIdx = find(selectedSatIdxGlobal == sceneSatIdx(iSat), 1, 'first');
  if ~isempty(matchIdx) && matchIdx <= size(seriesRaw, 1)
    seriesView(iSat, :) = seriesRaw(matchIdx, 1:numFrame);
  end
end
end


function seriesNorm = localNormalizeTruthSeriesField(fieldValue, numFrame)
%LOCALNORMALIZETRUTHSERIESFIELD Normalize truth field to Nsat-by-Nframe.

if isempty(fieldValue)
  seriesNorm = [];
  return;
end
fieldValue = reshape(fieldValue, size(fieldValue, 1), []);
if size(fieldValue, 2) == 1
  seriesNorm = repmat(fieldValue, 1, numFrame);
else
  seriesNorm = fieldValue(:, 1:min(numFrame, size(fieldValue, 2)));
  if size(seriesNorm, 2) < numFrame
    seriesNorm(:, end + 1:numFrame) = repmat(seriesNorm(:, end), 1, numFrame - size(seriesNorm, 2));
  end
end
end


function refSatIdxLocal = localResolveRefSatIdxLocalForView(sceneRef, truth)
%LOCALRESOLVEREFSATIDXLOCALFORVIEW Resolve one view-local reference satellite index.

refSatIdxLocal = localGetFieldOrDefault(truth, 'refSatIdxLocal', 1);
refSatIdxGlobal = localGetFieldOrDefault(truth, 'refSatIdxGlobal', NaN);
sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
if isempty(sceneSatIdx) || isnan(refSatIdxGlobal)
  return;
end
matchIdx = find(sceneSatIdx == refSatIdxGlobal, 1, 'first');
if ~isempty(matchIdx)
  refSatIdxLocal = matchIdx;
end
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
%LOCALBUILDTRUTHFDLOCALFORVIEW Build one single-frame truth local-frequency vector.

numSat = sceneRef.numSat;
fdLocalTrue = nan(numSat, 1);
selectedSatIdxGlobal = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
fdSatTrueHz = reshape(localGetFieldOrDefault(truth, 'fdSatTrueHz', []), [], 1);
if isempty(selectedSatIdxGlobal) || isempty(fdSatTrueHz)
  return;
end

sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
for iSat = 1:min(numSat, numel(sceneSatIdx))
  matchIdx = find(selectedSatIdxGlobal == sceneSatIdx(iSat), 1, 'first');
  if ~isempty(matchIdx) && matchIdx <= numel(fdSatTrueHz)
    fdLocalTrue(iSat) = fdSatTrueHz(matchIdx);
  end
end
end


function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
%LOCALEXTRACTSCENELOCALDOA Stack one scene sequence local-DoA tensor.

numSat = sceneSeq.numSat;
numFrame = sceneSeq.numFrame;
truthLocalDoa = nan(2, numSat, numFrame);
for iFrame = 1:numFrame
  sceneUse = sceneSeq.sceneCell{iFrame};
  if isfield(sceneUse, 'localDoa') && ~isempty(sceneUse.localDoa)
    localDoaUse = sceneUse.localDoa;
    truthLocalDoa(:, 1:size(localDoaUse, 2), iFrame) = localDoaUse;
  end
end
end


function out = localMergeStruct(base, override)
%LOCALMERGESTRUCT Merge two scalar structs with override precedence.

out = base;
if nargin < 2 || isempty(override)
  return;
end
if ~isstruct(override)
  error('buildDoaDopplerDynamicTransitionBundle:InvalidOverride', ...
    'Override data must be a struct.');
end
fieldList = fieldnames(override);
for iField = 1:numel(fieldList)
  out.(fieldList{iField}) = override.(fieldList{iField});
end
end


function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row compactly.

if isempty(valueRow)
  textOut = '[]';
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '), ']'];
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
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
