function caseDynMsKnown = runDynamicMsKnownTransitionCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, bestStaticMsCase)
%RUNDYNAMICMSKNOWNTRANSITIONCASE Run one MS-MF-CP-K case.

dynMsKnownOpt = flowOpt.dynBaseOpt;
dynMsKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth;
dynMsKnownOpt.enableFdAliasUnwrap = true;
dynMsKnownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsKnownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsKnownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsKnownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsKnownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsKnownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsKnownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynMsKnownOpt, true, ...
  periodicFixture.debugTruthMs, buildDynamicInitParamFromCase(bestStaticMsCase, true, periodicFixture.truth.fdRateFit));
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
