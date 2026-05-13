function [obj, pathGain, noiseVar, aux] = evaluateDoaDopplerMfObjective(model, optVar)
%EVALUATEDOADOPPLERMFOBJECTIVE Concentrated MF dynamic pilot objective.
% Keep hypothesis construction and profile-likelihood evaluation in one
% functional module so the branch solver only handles branch orchestration.

arguments
  model (1,1) struct
  optVar (:,1) double
end

hyp = buildDoaDopplerMfHypothesis(model, optVar);
[obj, prof, profAux] = evalDoaDopplerDynProfileLike(model, hyp);

if nargout <= 1
  return;
end

pathGain = prof.pathGainEst;
noiseVar = prof.noiseVarEst;

aux = struct();
aux.doaParam = hyp.doaParam;
aux.eciAngle = hyp.eciAngle;
aux.eciDirection = hyp.eciDirection;
aux.latlon = hyp.latlon;
aux.userPosEci = hyp.userPosEci;
aux.userVelEci = hyp.userVelEci;
aux.refSatIdxLocal = localGetStructField(profAux, 'refSatIdxLocal', NaN);
aux.deltaFdSeq = hyp.deltaFdSeq;
aux.deltaFdRef = hyp.deltaFdRef;
aux.deltaFdRefRaw = localGetStructField(profAux, 'deltaFdRefRaw', hyp.deltaFdRef);
aux.deltaFdRefEval = localGetStructField(profAux, 'deltaFdRefEval', hyp.deltaFdRef);
aux.deltaFdRate = hyp.deltaFdRate;
aux.deltaFdRateRaw = localGetStructField(profAux, 'deltaFdRateRaw', hyp.deltaFdRate);
aux.deltaFdRateEval = localGetStructField(profAux, 'deltaFdRateEval', hyp.deltaFdRate);
aux.fdRef = hyp.fdRef;
aux.fdRate = hyp.fdRate;
aux.fdSat = hyp.fdSat;
aux.fdSatRaw = localGetStructField(profAux, 'fdSatRaw', hyp.fdSat);
aux.fdRateSat = hyp.fdRateSat;
aux.fdLocal = profAux.fdLocal;
aux.fdLocalRaw = localGetStructField(profAux, 'fdLocalRaw', profAux.fdLocal);
aux.fdSatEval = profAux.fdSatEval;
aux.fdLocalEval = profAux.fdLocalEval;
aux.fdAliasStepHz = profAux.fdAliasStepHz;
aux.fdAliasIndex = profAux.fdAliasIndex;
aux.fdAliasShiftHz = profAux.fdAliasShiftHz;
aux.localDoaArrUsed = profAux.localDoaArrUsed;
aux.localDoaArr = profAux.localDoaArrUsed;
aux.phaseSat = prof.phaseSatEst;
aux.framePhase = prof.framePhaseEst;
aux.ampEst = prof.ampEst;
aux.noiseVarGlobal = prof.noiseVarGlobal;
aux.residualNorm = prof.residualNorm;
aux.countPerSat = prof.countPerSat;
aux.residualSat = localGetStructField(prof, 'residualSat', []);
aux.fitValueSat = localGetStructField(prof, 'fitValueSat', []);
aux.objectiveSat = localGetStructField(prof, 'objectiveSat', []);
aux.effectiveFrameSupportSat = localGetStructField(prof, 'effectiveFrameSupportSat', []);
aux.effectiveFrameSupportRatioSat = localGetStructField(prof, 'effectiveFrameSupportRatioSat', []);
aux.negativeProjectionRatioSat = localGetStructField(prof, 'negativeProjectionRatioSat', []);
aux.collapsedFrameCountSat = localGetStructField(prof, 'collapsedFrameCountSat', []);
aux.satFitRatio = localGetStructField(prof, 'satFitRatio', []);
aux.refFitRatio = localGetStructField(prof, 'refFitRatio', NaN);
aux.nonRefFitRatioFloor = localGetStructField(prof, 'nonRefFitRatioFloor', NaN);
aux.refSupportRatio = localGetStructField(prof, 'refSupportRatio', NaN);
aux.nonRefSupportRatioFloor = localGetStructField(prof, 'nonRefSupportRatioFloor', NaN);
aux.refConsistencyNorm = localGetStructField(prof, 'refConsistencyNorm', NaN);
aux.nonRefConsistencyRatioFloor = localGetStructField(prof, 'nonRefConsistencyRatioFloor', NaN);
aux.maxNonRefNegativeProjectionRatio = localGetStructField(prof, 'maxNonRefNegativeProjectionRatio', NaN);
aux.nonRefFitFloorPenalty = localGetStructField(prof, 'nonRefFitFloorPenalty', 0);
aux.nonRefSupportFloorPenalty = localGetStructField(prof, 'nonRefSupportFloorPenalty', 0);
aux.additionalObjectivePenalty = localGetStructField(prof, 'additionalObjectivePenalty', 0);
aux.blockValue = prof.blockValue;
aux.blockNorm2 = prof.blockNorm2;
aux.countPerBlock = localGetStructField(prof, 'countPerBlock', []);
aux.zMat = prof.zMat;
aux.etaMat = prof.etaMat;
end



function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct/object field with a default value.

fieldValue = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end
if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
