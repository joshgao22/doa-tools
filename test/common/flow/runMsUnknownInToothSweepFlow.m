function flow = runMsUnknownInToothSweepFlow(periodicFixture, subsetFixtureInput, caseWide, truth, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, weightSweepAlpha, ...
  doaOnlyOpt, staticBaseOpt, dynBaseOpt, staticInitDoaHalfWidth, ...
  subsetSelectDoaHalfWidthDeg, inToothFdHalfWidthHz, inToothFdRateHalfWidthHzPerSec, ...
  doaHalfWidthListDeg, toothStepHz, parallelOpt)
%RUNMSUNKNOWNINTOOTHSWEEPFLOW Run the regression-style CP-U in-tooth flow.
% This helper follows the validated regression path:
%   1) keep the caller-provided periodic-wide CP-U result as baseline;
%   2) run one deterministic-plus-random subset bank to select the tooth;
%   3) return to the periodic full-data window and refine inside one narrow
%      in-tooth [fdRef, fdRate] box;
%   4) sweep a small list of DoA half-widths and keep the best stable case.
%
% The helper is intentionally regression-oriented. It does not add any
% extra release bank or top-level source arbitration.

arguments
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
  staticInitDoaHalfWidth (:, 1) double
  subsetSelectDoaHalfWidthDeg (:, 1) double
  inToothFdHalfWidthHz (1, 1) double
  inToothFdRateHalfWidthHzPerSec (1, 1) double
  doaHalfWidthListDeg (:, 1) double
  toothStepHz (1, 1) double
  parallelOpt (1, 1) struct
end

subsetFixtureCell = localNormalizeSubsetFixtureBank(subsetFixtureInput);
numSubset = numel(subsetFixtureCell);
subsetSummaryCell = cell(numSubset, 1);
subsetCaseCell = cell(numSubset, 1);
useParfor = localShouldUseParfor(parallelOpt, numSubset, 'minSubsetBankForParfor');
if useParfor
  parfor iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    unknownOverride = struct();
    unknownOverride.initDoaHalfWidth = subsetSelectDoaHalfWidthDeg(:);
    unknownOverride.useStaticOnly = true;
    subsetCaseCell{iSubset} = localRunMsUnknownCase( ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      staticInitDoaHalfWidth(:), subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, ...
      struct(), unknownOverride);
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
else
  for iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    unknownOverride = struct();
    unknownOverride.initDoaHalfWidth = subsetSelectDoaHalfWidthDeg(:);
    unknownOverride.useStaticOnly = true;
    subsetCaseCell{iSubset} = localRunMsUnknownCase( ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      staticInitDoaHalfWidth(:), subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, ...
      struct(), unknownOverride);
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
end

wideSummary = localBuildCaseSummary(caseWide, truth, toothStepHz);
bestSubsetIdx = localSelectBestSubsetSummary(subsetSummaryCell, parallelOpt);
bestSubsetCase = subsetCaseCell{bestSubsetIdx};
bestSubsetSummary = subsetSummaryCell{bestSubsetIdx};
fdRefCenterInTooth = localResolveInToothFdRefCenter(caseWide, wideSummary, bestSubsetSummary, toothStepHz);
fdRateCenterInTooth = bestSubsetCase.estResult.fdRateEst;
fdRangeInTooth = [fdRefCenterInTooth - inToothFdHalfWidthHz, ...
  fdRefCenterInTooth + inToothFdHalfWidthHz];
fdRateRangeInTooth = [fdRateCenterInTooth - inToothFdRateHalfWidthHzPerSec, ...
  fdRateCenterInTooth + inToothFdRateHalfWidthHzPerSec];
directInitParamInTooth = buildDynamicInitParamFromCase(bestSubsetCase, false, fdRateCenterInTooth);
directInitParamInTooth = localOverwriteInToothCenter(directInitParamInTooth, fdRefCenterInTooth, fdRateCenterInTooth);

numWidth = numel(doaHalfWidthListDeg);
sweepSummaryCell = cell(numWidth, 1);
sweepCaseCell = cell(numWidth, 1);
if useParfor
  parfor iWidth = 1:numWidth
    halfWidthUse = doaHalfWidthListDeg(iWidth) * ones(2, 1);
    unknownOverride = struct();
    unknownOverride.initDoaParam = bestSubsetCase.estResult.doaParamEst(:);
    unknownOverride.initDoaHalfWidth = halfWidthUse;
    unknownOverride.directInitParam = directInitParamInTooth;
    unknownOverride.candidateDoaParam = bestSubsetCase.estResult.doaParamEst(:);
    sweepCaseCell{iWidth} = localRunMsUnknownCase( ...
      periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      staticInitDoaHalfWidth(:), fdRangeInTooth, fdRateRangeInTooth, struct(), unknownOverride);
    sweepSummaryCell{iWidth} = localBuildCaseSummary(sweepCaseCell{iWidth}, truth, toothStepHz);
    sweepSummaryCell{iWidth}.doaHalfWidthDeg = doaHalfWidthListDeg(iWidth);
  end
else
  for iWidth = 1:numWidth
    halfWidthUse = doaHalfWidthListDeg(iWidth) * ones(2, 1);
    unknownOverride = struct();
    unknownOverride.initDoaParam = bestSubsetCase.estResult.doaParamEst(:);
    unknownOverride.initDoaHalfWidth = halfWidthUse;
    unknownOverride.directInitParam = directInitParamInTooth;
    unknownOverride.candidateDoaParam = bestSubsetCase.estResult.doaParamEst(:);
    sweepCaseCell{iWidth} = localRunMsUnknownCase( ...
      periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      staticInitDoaHalfWidth(:), fdRangeInTooth, fdRateRangeInTooth, struct(), unknownOverride);
    sweepSummaryCell{iWidth} = localBuildCaseSummary(sweepCaseCell{iWidth}, truth, toothStepHz);
    sweepSummaryCell{iWidth}.doaHalfWidthDeg = doaHalfWidthListDeg(iWidth);
  end
end

[bestAngleIdx, bestStableIdx] = localSelectBestSweep(sweepSummaryCell);
bestAngleCase = sweepCaseCell{bestAngleIdx};
bestStableCase = sweepCaseCell{bestStableIdx};
bestAngleSummary = sweepSummaryCell{bestAngleIdx};
bestStableSummary = sweepSummaryCell{bestStableIdx};

toothFlowTable = localBuildToothFlowSummaryTable(wideSummary, bestSubsetSummary, bestStableSummary);
subsetCandidateTable = localBuildSubsetCandidateSummaryTable(subsetSummaryCell);
sweepSummaryTable = localBuildSweepSummaryTable(sweepSummaryCell);

flow = struct();
flow.modeTag = "unknown";
flow.caseWide = caseWide;
flow.caseSubset = bestSubsetCase;
flow.caseBestAngle = bestAngleCase;
flow.caseBestStable = bestStableCase;
flow.caseFinal = bestStableCase;
flow.wideSummary = wideSummary;
flow.subsetSummary = bestSubsetSummary;
flow.bestAngleSummary = bestAngleSummary;
flow.bestStableSummary = bestStableSummary;
flow.refineSummary = bestStableSummary;
flow.toothFlowTable = toothFlowTable;
flow.subsetCandidateTable = subsetCandidateTable;
flow.sweepSummaryTable = sweepSummaryTable;
flow.bestSubsetIdx = bestSubsetIdx;
flow.bestAngleIdx = bestAngleIdx;
flow.bestStableIdx = bestStableIdx;
flow.bestAngleWidthDeg = bestAngleSummary.doaHalfWidthDeg;
flow.bestStableWidthDeg = bestStableSummary.doaHalfWidthDeg;
flow.selectedSubsetLabel = bestSubsetSummary.subsetLabel;
flow.selectedSubsetOffsets = bestSubsetSummary.subsetOffsetIdx;
flow.fdRangeInTooth = fdRangeInTooth;
flow.fdRateRangeInTooth = fdRateRangeInTooth;
flow.usedParfor = useParfor;
end

function caseDynMsUnknown = localRunMsUnknownCase(fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, weightSweepAlpha, doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  staticInitDoaHalfWidth, fdRangeUse, fdRateRangeUse, dynKnownOverride, dynUnknownOverride)
%LOCALRUNMSUNKNOWNCASE Run the regression-style MS-MF-CP-U path on one fixture.

fixture = localEnsureFixtureViews(fixture);

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, fixture.wavelen, ...
  pilotWave, carrierFreq, sampleRate, fdRangeUse, fixture.truth, ...
  fixture.otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha(:).', staticInitDoaHalfWidth(:));
bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;
truthDebug = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = [0.003; 0.003];
dynMsKnownOpt = localMergeStruct(dynMsKnownOpt, dynKnownOverride);

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsKnownOpt, true, ...
  truthDebug, buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

useStaticOnly = logical(localGetFieldOrDefault(dynUnknownOverride, 'useStaticOnly', false));
directInitParam = localGetFieldOrDefault(dynUnknownOverride, 'directInitParam', []);
candidateDoaParam = localGetFieldOrDefault(dynUnknownOverride, 'candidateDoaParam', []);
dynUnknownOptUse = localStripControlFields(dynUnknownOverride, {'useStaticOnly', 'directInitParam', 'candidateDoaParam'});

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];
dynMsUnknownOpt = localMergeStruct(dynMsUnknownOpt, dynUnknownOptUse);

if ~isempty(candidateDoaParam)
  dynMsUnknownOpt.initDoaParam = reshape(candidateDoaParam, [], 1);
end

if ~isempty(directInitParam)
  initParamUse = directInitParam;
elseif useStaticOnly
  initParamUse = initParamMsUnknownStatic;
else
  initParamUse = buildUnknownInitCandidateSet( ...
    caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
    dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
end

caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, ...
  truthDebug, initParamUse);
end

function summary = localBuildCaseSummary(caseUse, truth, toothStepHz)
summary = struct();
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
summary.finalObj = localResolveCaseFinalObj(caseUse);
summary.subsetLabel = string.empty(1, 0);
summary.subsetOffsetIdx = [];
summary.doaHalfWidthDeg = NaN;
end

function fdRefCenterHz = localResolveInToothFdRefCenter(caseWide, wideSummary, subsetSummary, toothStepHz)
wideFdRefHz = caseWide.estResult.fdRefEst;
targetToothIdx = subsetSummary.toothIdx;
wideToothIdx = wideSummary.toothIdx;
fdRefCenterHz = wideFdRefHz + (targetToothIdx - wideToothIdx) * toothStepHz;
end

function finalObj = localResolveCaseFinalObj(caseUse)
finalObj = NaN;
aux = localGetFieldOrDefault(caseUse.estResult, 'aux', struct());
debugInfo = localGetFieldOrDefault(aux, 'debug', struct());
finalEval = localGetFieldOrDefault(debugInfo, 'finalEval', struct());
finalObj = localGetFieldOrDefault(finalEval, 'obj', NaN);
if ~isfinite(finalObj)
  finalObj = localGetFieldOrDefault(caseUse.estResult, 'fval', NaN);
end
end

function initParamOut = localOverwriteInToothCenter(initParamIn, fdRefCenterHz, fdRateCenterHz)
initParamOut = initParamIn;
if isempty(initParamOut)
  return;
end
initParamOut = reshape(initParamOut, [], 1);
if numel(initParamOut) >= 3 && isfinite(fdRefCenterHz)
  initParamOut(3) = fdRefCenterHz;
end
if numel(initParamOut) >= 4 && isfinite(fdRateCenterHz)
  initParamOut(4) = fdRateCenterHz;
end
end

function bestIdx = localSelectBestSubsetSummary(summaryCell, parallelOpt)
numCase = numel(summaryCell);
forceSubsetLabel = string(localGetFieldOrDefault(parallelOpt, 'forceSubsetLabel', ""));
selectMask = true(numCase, 1);
if strlength(forceSubsetLabel) > 0
  labelList = strings(numCase, 1);
  for iCase = 1:numCase
    labelList(iCase) = string(localGetFieldOrDefault(summaryCell{iCase}, 'subsetLabel', ""));
  end
  matchedMask = strcmp(labelList, forceSubsetLabel);
  if any(matchedMask)
    selectMask = matchedMask;
  end
end

scoreMat = inf(numCase, 5);
for iCase = 1:numCase
  if ~selectMask(iCase)
    continue;
  end
  s = summaryCell{iCase};
  if ~logical(s.isResolved)
    continue;
  end
  scoreMat(iCase, :) = [abs(s.toothIdx), abs(s.toothResidualHz), s.angleErrDeg, ...
    abs(s.fdRateErrHzPerSec), abs(s.fdRefErrHz)];
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
bestIdx = orderIdx(1);
end

function [bestAngleIdx, bestStableIdx] = localSelectBestSweep(summaryCell)
numCase = numel(summaryCell);
scoreAngle = inf(numCase, 5);
scoreStable = inf(numCase, 5);
for iCase = 1:numCase
  s = summaryCell{iCase};
  if ~logical(s.isResolved)
    continue;
  end
  scoreAngle(iCase, :) = [abs(s.toothIdx), s.angleErrDeg, abs(s.toothResidualHz), ...
    abs(s.fdRefErrHz), abs(s.fdRateErrHzPerSec)];
  finalObjUse = localGetFieldOrDefault(s, 'finalObj', NaN);
  if ~isfinite(finalObjUse)
    finalObjUse = abs(s.fdRefErrHz);
  end
  scoreStable(iCase, :) = [abs(s.toothIdx), finalObjUse, abs(s.toothResidualHz), ...
    abs(s.fdRateErrHzPerSec), s.angleErrDeg];
end
[~, orderAngle] = sortrows(scoreAngle, 1:size(scoreAngle, 2));
[~, orderStable] = sortrows(scoreStable, 1:size(scoreStable, 2));
bestAngleIdx = orderAngle(1);
bestStableIdx = orderStable(1);
end

function toothTable = localBuildToothFlowSummaryTable(wideSummary, subsetSummary, refineSummary)
stageName = ["periodic-wide"; "subset-select"; "periodic-in-tooth"];
solveVariant = [wideSummary.solveVariant; subsetSummary.solveVariant; refineSummary.solveVariant];
isResolved = [wideSummary.isResolved; subsetSummary.isResolved; refineSummary.isResolved];
angleErrDeg = [wideSummary.angleErrDeg; subsetSummary.angleErrDeg; refineSummary.angleErrDeg];
fdRefErrHz = [wideSummary.fdRefErrHz; subsetSummary.fdRefErrHz; refineSummary.fdRefErrHz];
fdRateErrHzPerSec = [wideSummary.fdRateErrHzPerSec; subsetSummary.fdRateErrHzPerSec; refineSummary.fdRateErrHzPerSec];
toothIdx = [wideSummary.toothIdx; subsetSummary.toothIdx; refineSummary.toothIdx];
toothResidualHz = [wideSummary.toothResidualHz; subsetSummary.toothResidualHz; refineSummary.toothResidualHz];
runTimeMs = [wideSummary.runTimeMs; subsetSummary.runTimeMs; refineSummary.runTimeMs];
toothTable = table(stageName, solveVariant, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs);
end

function subsetTable = localBuildSubsetCandidateSummaryTable(summaryCell)
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

function sweepTable = localBuildSweepSummaryTable(summaryCell)
numCase = numel(summaryCell);
doaHalfWidthDeg = nan(numCase, 1);
solveVariant = strings(numCase, 1);
isResolved = false(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
runTimeMs = nan(numCase, 1);
for iCase = 1:numCase
  s = summaryCell{iCase};
  doaHalfWidthDeg(iCase) = s.doaHalfWidthDeg;
  solveVariant(iCase) = s.solveVariant;
  isResolved(iCase) = logical(s.isResolved);
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
  runTimeMs(iCase) = s.runTimeMs;
end
sweepTable = table(doaHalfWidthDeg, solveVariant, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs);
end

function subsetFixtureCell = localNormalizeSubsetFixtureBank(subsetFixtureInput)
if iscell(subsetFixtureInput)
  subsetFixtureCell = reshape(subsetFixtureInput, [], 1);
else
  subsetFixtureCell = {subsetFixtureInput};
end
if isempty(subsetFixtureCell)
  error('runMsUnknownInToothSweepFlow:EmptySubsetBank', ...
    'The subset fixture bank must contain at least one candidate.');
end
end

function textOut = localFormatIntegerRow(valueRow)
if isempty(valueRow)
  textOut = "[]";
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = string(sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', ')));
end

function useParfor = localShouldUseParfor(parallelOpt, numCase, minFieldName)
useParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableSubsetSolveParfor', ...
  localGetFieldOrDefault(parallelOpt, 'enableParfor', true)));
if ~useParfor
  return;
end
minCase = localGetFieldOrDefault(parallelOpt, minFieldName, 12);
useParfor = useParfor && (numCase >= minCase) && localCanUseParfor();
end

function tf = localCanUseParfor()
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

function fixture = localEnsureFixtureViews(fixture)
if isfield(fixture, 'viewRefOnly') && isfield(fixture, 'viewOtherOnly') && isfield(fixture, 'viewMs') && ...
    ~isempty(fixture.viewRefOnly) && ~isempty(fixture.viewOtherOnly) && ~isempty(fixture.viewMs)
  return;
end
requiredFieldList = {'sceneSeq','sceneRef','rxSigCell','rxSigRef','refSatIdxLocal','otherSatIdxLocal', ...
  'gridSize','searchRange','E'};
for iField = 1:numel(requiredFieldList)
  fieldName = requiredFieldList{iField};
  if ~isfield(fixture, fieldName) || isempty(fixture.(fieldName))
    error('runMsUnknownInToothSweepFlow:MissingFixtureField', ...
      'The fixture is missing required field "%s" needed to rebuild viewRefOnly/viewOtherOnly/viewMs.', fieldName);
  end
end
sceneSeqRefOnly = selectSatSceneSeq(fixture.sceneSeq, fixture.refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(fixture.rxSigCell, fixture.refSatIdxLocal, 'multiFrame');
fixture.viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  fixture.gridSize, fixture.searchRange, fixture.E, struct('sceneSeq', sceneSeqRefOnly));
sceneOtherOnly = selectSatScene(fixture.sceneRef, fixture.otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(fixture.rxSigRef, fixture.otherSatIdxLocal, 'singleFrame');
fixture.viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  fixture.gridSize, fixture.searchRange, fixture.E);
fixture.viewMs = buildDoaDopplerEstView(fixture.sceneRef, fixture.rxSigCell, ...
  fixture.gridSize, fixture.searchRange, fixture.E, struct('sceneSeq', fixture.sceneSeq));
fixture.sceneSeqRefOnly = sceneSeqRefOnly;
fixture.sceneRefOnly = sceneRefOnly;
fixture.sceneOtherOnly = sceneOtherOnly;
end

function out = localMergeStruct(base, override)
out = base;
if isempty(override)
  return;
end
fieldList = fieldnames(override);
for iField = 1:numel(fieldList)
  out.(fieldList{iField}) = override.(fieldList{iField});
end
end

function out = localStripControlFields(inStruct, fieldNameList)
out = inStruct;
if isempty(out)
  return;
end
for iField = 1:numel(fieldNameList)
  fieldName = fieldNameList{iField};
  if isfield(out, fieldName)
    out = rmfield(out, fieldName);
  end
end
end

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
fieldValue = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
