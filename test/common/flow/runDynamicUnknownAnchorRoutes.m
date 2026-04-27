function [caseDynMsUnknownAnchor, anchorUnknownSummary, caseDynMsUnknownAnchorDoaPolish, ...
  anchorDoaPolishSummary, fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, runMeta] = ...
  runDynamicUnknownAnchorRoutes(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, selectedSubsetCase, selectedSubsetSummary, knownSummary, ...
  truth, toothStepHz, buildUnknownSummaryFn)
%RUNDYNAMICUNKNOWNANCHORROUTES Run the subset-anchor route family.
% Group the fixed-DoA subset anchor and the optional anchor-DoA polish into
% one functional module so the transition bundle only decides whether to run
% the anchor family.

arguments
  periodicFixture (1,1) struct
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
  selectedSubsetCase (1,1) struct
  selectedSubsetSummary (1,1) struct
  knownSummary (1,1) struct
  truth (1,1) struct
  toothStepHz (1,1) double
  buildUnknownSummaryFn (1,1) function_handle
end

fdRangeSubsetAnchor = [selectedSubsetCase.estResult.fdRefEst - flowOpt.subsetAnchorFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.subsetAnchorFdHalfWidthHz];
fdRateRangeSubsetAnchor = [selectedSubsetCase.estResult.fdRateEst - flowOpt.subsetAnchorFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.subsetAnchorFdRateHalfWidthHzPerSec];

runMeta = struct();
runMeta.hasAnchorDoaPolish = false;

stageTimer = tic;
caseDynMsUnknownAnchor = localRunMsUnknownSubsetAnchor(periodicFixture, pilotWave, ...
  carrierFreq, sampleRate, optVerbose, flowOpt, ...
  fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, selectedSubsetCase, ...
  selectedSubsetSummary, knownSummary);
anchorUnknownSummary = localBuildRouteSummary(caseDynMsUnknownAnchor, ...
  "periodic-fd-anchor", "anchor", buildUnknownSummaryFn, truth, toothStepHz);
runMeta.subsetAnchorSec = toc(stageTimer);

caseDynMsUnknownAnchorDoaPolish = struct();
anchorDoaPolishSummary = localBuildSkippedUnknownSummary("periodic-fd-anchor-doa-polish", "anchor");
stageTimer = tic;
if localShouldRunAnchorDoaPolish(anchorUnknownSummary, flowOpt)
  caseDynMsUnknownAnchorDoaPolish = localRunMsUnknownAnchorDoaPolish(periodicFixture, pilotWave, ...
    carrierFreq, sampleRate, optVerbose, flowOpt, caseDynMsUnknownAnchor);
  anchorDoaPolishSummary = localBuildRouteSummary(caseDynMsUnknownAnchorDoaPolish, ...
    "periodic-fd-anchor-doa-polish", "anchor", buildUnknownSummaryFn, truth, toothStepHz);
  runMeta.hasAnchorDoaPolish = localCaseHasUsableEstimate(caseDynMsUnknownAnchorDoaPolish);
end
runMeta.anchorDoaPolishSec = toc(stageTimer);
end


function caseDynMsUnknown = localRunMsUnknownSubsetAnchor(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, selectedSubsetCase, ...
  selectedSubsetSummary, knownSummary)
%LOCALRUNMSUNKNOWNSUBSETANCHOR Run one same-tooth subset anchor solve.

runLabel = string(localGetFieldOrDefault(flowOpt, 'subsetAnchorRunLabel', "MS-MF-CP-U-anchor"));
dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
[subsetAnchorDoaInitParam, ~] = resolveRecoveredToothDoaSeed( ...
  flowOpt, selectedSubsetCase, selectedSubsetSummary, knownSummary);
dynMsUnknownOpt.initDoaParam = subsetAnchorDoaInitParam(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.subsetAnchorDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.freezeDoa = logical(flowOpt.subsetAnchorFreezeDoa);
dynMsUnknownOpt.disableUnknownDoaReleaseFloor = true;
dynMsUnknownOpt.unknownDoaReleaseHalfWidth = flowOpt.subsetAnchorDoaHalfWidthDeg(:);
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

initParamMsUnknown = buildDynamicInitParamFromCase(selectedSubsetCase, false, selectedSubsetCase.estResult.fdRateEst);
caseDynMsUnknown = runDynamicDoaDopplerCase(runLabel, "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, initParamMsUnknown);
end


function caseDynMsUnknown = localRunMsUnknownAnchorDoaPolish(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, anchorCase)
%LOCALRUNMSUNKNOWNANCHORDOAPOLISH Release DoA around one trusted anchor.

runLabel = string(localGetFieldOrDefault(flowOpt, 'anchorDoaPolishRunLabel', "MS-MF-CP-U-anchorDoaPolish"));
dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = anchorCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.anchorDoaPolishDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.freezeDoa = false;
dynMsUnknownOpt.disableUnknownDoaReleaseFloor = logical(localGetFieldOrDefault(flowOpt, 'anchorDoaPolishDisableUnknownDoaReleaseFloor', true));
dynMsUnknownOpt.unknownDoaReleaseHalfWidth = flowOpt.anchorDoaPolishDoaHalfWidthDeg(:);
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);
dynMsUnknownOpt.optimOpt = localMergeOptimOpt(localGetFieldOrDefault(dynMsUnknownOpt, 'optimOpt', struct()), ...
  localGetFieldOrDefault(flowOpt, 'anchorDoaPolishOptimOpt', struct()));

fdRefCenter = anchorCase.estResult.fdRefEst;
fdRateCenter = localGetCaseFdRateEst(anchorCase, periodicFixture.truth.fdRateFit);
fdRangeAnchor = [fdRefCenter - flowOpt.anchorDoaPolishFdHalfWidthHz, ...
  fdRefCenter + flowOpt.anchorDoaPolishFdHalfWidthHz];
fdRateRangeAnchor = [fdRateCenter - flowOpt.anchorDoaPolishFdRateHalfWidthHzPerSec, ...
  fdRateCenter + flowOpt.anchorDoaPolishFdRateHalfWidthHzPerSec];
initParamMsUnknown = buildDynamicInitParamFromCase(anchorCase, false, fdRateCenter);
caseDynMsUnknown = runDynamicDoaDopplerCase(runLabel, "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeAnchor, fdRateRangeAnchor, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, initParamMsUnknown);
end

function tf = localShouldRunAnchorDoaPolish(anchorUnknownSummary, flowOpt)
%LOCALSHOULDRUNANCHORDOAPOLISH Decide whether one anchor-centered DoA polish is needed.

tf = false;
if ~logical(localGetFieldOrDefault(flowOpt, 'enableAnchorDoaPolish', false))
  return;
end
if ~logical(localGetFieldOrDefault(anchorUnknownSummary, 'isResolved', false))
  return;
end
toothIdx = localGetFieldOrDefault(anchorUnknownSummary, 'toothIdx', NaN);
if ~(isfinite(toothIdx) && abs(toothIdx) == 0)
  return;
end
toothResidualHz = abs(localGetFieldOrDefault(anchorUnknownSummary, 'toothResidualHz', inf));
if ~(isfinite(toothResidualHz) && toothResidualHz <= localGetFieldOrDefault(flowOpt, 'anchorDoaPolishCentralResidualTolHz', 50))
  return;
end
fitFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefFitRatioFloor', NaN);
supportFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefSupportRatioFloor', NaN);
fitReady = ~(isfinite(fitFloor) && fitFloor >= localGetFieldOrDefault(flowOpt, 'anchorDoaPolishMinNonRefFitFloorToSkip', 0.9));
supportReady = ~(isfinite(supportFloor) && supportFloor >= localGetFieldOrDefault(flowOpt, 'anchorDoaPolishMinNonRefSupportFloorToSkip', 0.8));
tf = fitReady || supportReady;
end


function summary = localBuildRouteSummary(caseUse, stageTag, routeFamily, buildUnknownSummaryFn, truth, toothStepHz)
%LOCALBUILDROUTESUMMARY Normalize one route summary to the shared field set.

baseSummary = buildDynamicUnknownCaseSummary(caseUse, NaN, struct());
summary = baseSummary;
summaryUse = buildUnknownSummaryFn(caseUse, truth, toothStepHz);
summary = localOverlaySummary(summary, summaryUse);
summary.stageTag = string(stageTag);
summary.routeFamily = string(routeFamily);
end


function summary = localBuildSkippedUnknownSummary(stageTag, routeFamily)
%LOCALBUILDSKIPPEDUNKNOWNSUMMARY Build one placeholder summary for skipped routes.

summary = buildDynamicUnknownCaseSummary(struct(), NaN, struct());
summary.solveVariant = "skipped";
summary.isResolved = false;
summary.runTimeMs = 0;
summary.stageTag = string(stageTag);
summary.routeFamily = string(routeFamily);
end


function summary = localOverlaySummary(baseSummary, overrideSummary)
%LOCALOVERLAYSUMMARY Overlay one summary struct on top of a canonical base.

summary = baseSummary;
if ~isstruct(overrideSummary)
  return;
end
fieldList = fieldnames(overrideSummary);
for iField = 1:numel(fieldList)
  summary.(fieldList{iField}) = overrideSummary.(fieldList{iField});
end
end


function tf = localCaseHasUsableEstimate(caseUse)
%LOCALCASEHASUSABLEESTIMATE Return true when one case contains finite key estimates.

tf = false;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult) || ~isstruct(estResult)
  return;
end
fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
doaParamEst = localGetFieldOrDefault(estResult, 'doaParamEst', nan(2, 1));
tf = isfinite(fdRefEst) && isfinite(fdRateEst) && all(isfinite(doaParamEst(:)));
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

function optimOpt = localMergeOptimOpt(baseOptimOpt, overrideOptimOpt)
%LOCALMERGEOPTIMOPT Merge one local optimizer override struct.

optimOpt = baseOptimOpt;
if nargin < 1 || isempty(optimOpt)
  optimOpt = struct();
end
if nargin < 2 || isempty(overrideOptimOpt)
  return;
end
if ~isstruct(overrideOptimOpt)
  error('runDynamicUnknownAnchorRoutes:InvalidLocalOptimOpt', ...
    'Optimizer overrides must be provided as a struct.');
end
fieldList = fieldnames(overrideOptimOpt);
for iField = 1:numel(fieldList)
  optimOpt.(fieldList{iField}) = overrideOptimOpt.(fieldList{iField});
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
