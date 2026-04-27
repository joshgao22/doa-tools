function caseDynMsUnknown = runDynamicMsUnknownSubsetCase(fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeUse, fdRateRangeUse, subsetSeedInfo)
%RUNDYNAMICMSUNKNOWNSUBSETCASE Run one MS-MF-CP-U case on one frame subset.

if nargin < 9
  subsetSeedInfo = struct();
end
reusePeriodicSeed = logical(localGetFieldOrDefault(subsetSeedInfo, 'reusePeriodicSeeds', false));
if reusePeriodicSeed
  bestStaticMsCase = localGetFieldOrDefault(subsetSeedInfo, 'bestStaticMsCase', struct());
  knownSeedCase = localGetFieldOrDefault(subsetSeedInfo, 'knownSeedCase', struct());
  staticInitDoaHalfWidth = localGetFieldOrDefault(subsetSeedInfo, 'staticInitDoaHalfWidth', [0.01; 0.01]);
else
  caseBundle = buildDoaDopplerStaticTransitionBundle( ...
    fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, fixture.wavelen, ...
    pilotWave, carrierFreq, sampleRate, fdRangeUse, fixture.truth, ...
    fixture.otherSatIdxGlobal, optVerbose, flowOpt.doaOnlyOpt, ...
    flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));
  bestStaticMsCase = caseBundle.bestStaticMsCase;
  knownSeedCase = struct();
  staticInitDoaHalfWidth = localGetFieldOrDefault(caseBundle.staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]);
end

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.subsetSelectDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);
initSeed = initParamMsUnknownStatic;
if reusePeriodicSeed && localCaseHasUsableEstimate(knownSeedCase)
  fdRateSeedKnown = localGetCaseFdRateEst(knownSeedCase, fixture.truth.fdRateFit);
  initParamMsUnknownCpK = buildDynamicInitParamFromCase(knownSeedCase, false, fdRateSeedKnown);
  initSeed = buildUnknownInitCandidateSet( ...
    knownSeedCase, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
    dynMsUnknownOpt.initDoaHalfWidth, staticInitDoaHalfWidth);
end
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, fixture.debugTruthMs, initSeed);
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
