function [caseDynMsUnknownFdWide, caseDynMsUnknownWide, fdWideUnknownSummary, wideUnknownSummary] = ...
  runDynamicUnknownWideRoutes(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, wideSeedCase, truth, toothStepHz, runWideRefine, runFdWideFirst, ...
  buildUnknownSummaryFn)
%RUNDYNAMICUNKNOWNWIDEROUTES Run the optional wide unknown-route family.
% Keep the fd-wide tooth escape and the following periodic-wide refine in
% one helper so the transition bundle only decides whether this route family
% should run.

arguments
  periodicFixture (1,1) struct
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
  wideSeedCase (1,1) struct
  truth (1,1) struct
  toothStepHz (1,1) double
  runWideRefine (1,1) logical
  runFdWideFirst (1,1) logical
  buildUnknownSummaryFn (1,1) function_handle
end

caseDynMsUnknownFdWide = struct();
fdWideUnknownSummary = localBuildSkippedUnknownSummary("periodic-fd-wide");

if runFdWideFirst
  caseDynMsUnknownFdWide = localRunMsUnknownWideCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, wideSeedCase);
  fdWideUnknownSummary = buildUnknownSummaryFn(caseDynMsUnknownFdWide, truth, toothStepHz);
  wideRefineSeedCase = caseDynMsUnknownFdWide;
else
  wideRefineSeedCase = wideSeedCase;
end

if runWideRefine
  caseDynMsUnknownWide = localRunMsUnknownWideRefineCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, wideRefineSeedCase);
  wideUnknownSummary = buildUnknownSummaryFn(caseDynMsUnknownWide, truth, toothStepHz);
elseif runFdWideFirst
  caseDynMsUnknownWide = caseDynMsUnknownFdWide;
  wideUnknownSummary = fdWideUnknownSummary;
else
  caseDynMsUnknownWide = wideRefineSeedCase;
  wideUnknownSummary = buildUnknownSummaryFn(caseDynMsUnknownWide, truth, toothStepHz);
end
end


function caseDynMsUnknownWide = localRunMsUnknownWideCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, wideSeedCase)
%LOCALRUNMSUNKNOWNWIDECASE Run one fixed-DoA fd-wide tooth escape.

runLabel = string(localGetFieldOrDefault(flowOpt, 'fdWideRunLabel', "MS-MF-CP-U-fdWide"));
dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = wideSeedCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.inToothDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.freezeDoa = logical(flowOpt.fdWideFreezeDoa);
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

fdRateSeed = localGetCaseFdRateEst(wideSeedCase, periodicFixture.truth.fdRateFit);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(wideSeedCase, false, fdRateSeed);
caseDynMsUnknownWide = runDynamicDoaDopplerCase(runLabel, "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, initParamMsUnknownStatic);
end


function caseDynMsUnknownWide = localRunMsUnknownWideRefineCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, wideSeedCase)
%LOCALRUNMSUNKNOWNWIDEREFINECASE Run one periodic-wide refine from a wide seed.

runLabel = string(localGetFieldOrDefault(flowOpt, 'wideRefineRunLabel', "MS-MF-CP-U-wide"));
dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = wideSeedCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.wideRefineDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.freezeDoa = false;
dynMsUnknownOpt.disableUnknownDoaReleaseFloor = logical(localGetFieldOrDefault(flowOpt, 'wideRefineDisableUnknownDoaReleaseFloor', true));
dynMsUnknownOpt.unknownDoaReleaseHalfWidth = flowOpt.wideRefineDoaHalfWidthDeg(:);
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

fdRefCenter = wideSeedCase.estResult.fdRefEst;
fdRateCenter = localGetCaseFdRateEst(wideSeedCase, periodicFixture.truth.fdRateFit);
fdRangeWide = [fdRefCenter - flowOpt.wideRefineFdHalfWidthHz, fdRefCenter + flowOpt.wideRefineFdHalfWidthHz];
fdRateRangeWide = [fdRateCenter - flowOpt.wideRefineFdRateHalfWidthHzPerSec, ...
  fdRateCenter + flowOpt.wideRefineFdRateHalfWidthHzPerSec];
initParamMsUnknown = buildDynamicInitParamFromCase(wideSeedCase, false, fdRateCenter);
caseDynMsUnknownWide = runDynamicDoaDopplerCase(runLabel, "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeWide, fdRateRangeWide, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, initParamMsUnknown);
end


function summary = localBuildSkippedUnknownSummary(stageTag)
%LOCALBUILDSKIPPEDUNKNOWNSUMMARY Build one placeholder summary for skipped routes.

summary = struct();
summary.solveVariant = "skipped";
summary.isResolved = false;
summary.doaParamEst = nan(1, 2);
summary.fdRefEst = NaN;
summary.fdRateEst = NaN;
summary.angleErrDeg = NaN;
summary.fdRefErrHz = NaN;
summary.fdRateErrHzPerSec = NaN;
summary.toothIdx = NaN;
summary.toothResidualHz = NaN;
summary.runTimeMs = 0;
summary.finalObj = NaN;
summary.finalResidualNorm = NaN;
summary.stageTag = string(stageTag);
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
