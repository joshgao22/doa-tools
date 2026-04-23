function evalDiag = buildDoaDopplerMfEvalDiag(model, optVar, obj, noiseVar, evalAux)
%BUILDDOADOPPLERMFEVALDIAG Build one compact MF evaluation summary.

arguments
  model (1,1) struct
  optVar (:,1) double
  obj (1,1) double
  noiseVar (:,1) double
  evalAux (1,1) struct
end

evalDiag = struct();
evalDiag.optVar = optVar(:);
evalDiag.obj = obj;
evalDiag.doaType = model.doaType;
evalDiag.phaseMode = model.phaseMode;
evalDiag.fdRateMode = model.fdRateMode;
evalDiag.refSatIdxLocal = localGetStructField(evalAux, 'refSatIdxLocal', NaN);
evalDiag.doaParam = evalAux.doaParam;
evalDiag.eciAngle = evalAux.eciAngle;
evalDiag.latlon = evalAux.latlon;
evalDiag.fdRef = evalAux.fdRef;
evalDiag.fdRate = evalAux.fdRate;
evalDiag.deltaFdRef = evalAux.deltaFdRef(:);
evalDiag.deltaFdRefRaw = localGetStructField(evalAux, 'deltaFdRefRaw', evalAux.deltaFdRef(:));
evalDiag.deltaFdRefEval = localGetStructField(evalAux, 'deltaFdRefEval', evalAux.deltaFdRef(:));
evalDiag.deltaFdRate = evalAux.deltaFdRate(:);
evalDiag.deltaFdRateRaw = localGetStructField(evalAux, 'deltaFdRateRaw', evalAux.deltaFdRate(:));
evalDiag.deltaFdRateEval = localGetStructField(evalAux, 'deltaFdRateEval', evalAux.deltaFdRate(:));
evalDiag.fdSat = evalAux.fdSat;
evalDiag.fdSatRaw = localGetStructField(evalAux, 'fdSatRaw', evalAux.fdSat);
evalDiag.fdSatEval = evalAux.fdSatEval;
evalDiag.fdRateSat = evalAux.fdRateSat;
evalDiag.fdLocal = evalAux.fdLocal;
evalDiag.fdLocalRaw = localGetStructField(evalAux, 'fdLocalRaw', evalAux.fdLocal);
evalDiag.fdLocalEval = evalAux.fdLocalEval;
evalDiag.localDoaArrUsed = localGetStructField(evalAux, 'localDoaArrUsed', []);
evalDiag.fdAliasStepHz = evalAux.fdAliasStepHz;
evalDiag.fdAliasIndex = evalAux.fdAliasIndex;
evalDiag.fdAliasShiftHz = evalAux.fdAliasShiftHz;
evalDiag.noiseVar = noiseVar(:);
evalDiag.countPerSat = evalAux.countPerSat(:);
evalDiag.residualSat = localGetStructField(evalAux, 'residualSat', noiseVar(:) .* evalAux.countPerSat(:));
evalDiag.fitValueSat = localGetStructField(evalAux, 'fitValueSat', nan(size(evalDiag.residualSat)));
evalDiag.objectiveSat = localGetStructField(evalAux, 'objectiveSat', nan(size(evalDiag.residualSat)));
evalDiag.effectiveFrameSupportSat = localGetStructField(evalAux, 'effectiveFrameSupportSat', []);
evalDiag.effectiveFrameSupportRatioSat = localGetStructField(evalAux, 'effectiveFrameSupportRatioSat', []);
evalDiag.negativeProjectionRatioSat = localGetStructField(evalAux, 'negativeProjectionRatioSat', []);
evalDiag.collapsedFrameCountSat = localGetStructField(evalAux, 'collapsedFrameCountSat', []);
evalDiag.satFitRatio = localGetStructField(evalAux, 'satFitRatio', []);
evalDiag.refFitRatio = localGetStructField(evalAux, 'refFitRatio', NaN);
evalDiag.nonRefFitRatioFloor = localGetStructField(evalAux, 'nonRefFitRatioFloor', NaN);
evalDiag.refSupportRatio = localGetStructField(evalAux, 'refSupportRatio', NaN);
evalDiag.nonRefSupportRatioFloor = localGetStructField(evalAux, 'nonRefSupportRatioFloor', NaN);
evalDiag.refConsistencyNorm = localGetStructField(evalAux, 'refConsistencyNorm', NaN);
evalDiag.nonRefConsistencyRatioFloor = localGetStructField(evalAux, 'nonRefConsistencyRatioFloor', NaN);
evalDiag.maxNonRefNegativeProjectionRatio = localGetStructField(evalAux, 'maxNonRefNegativeProjectionRatio', NaN);
evalDiag.nonRefFitFloorPenalty = localGetStructField(evalAux, 'nonRefFitFloorPenalty', 0);
evalDiag.nonRefSupportFloorPenalty = localGetStructField(evalAux, 'nonRefSupportFloorPenalty', 0);
evalDiag.additionalObjectivePenalty = localGetStructField(evalAux, 'additionalObjectivePenalty', 0);
evalDiag.residualNorm = evalAux.residualNorm;
evalDiag.noiseVarGlobal = evalAux.noiseVarGlobal;
evalDiag.coherenceSat = localCalcSatCoherence(evalAux.zMat);
evalDiag.zMat = localGetStructField(evalAux, 'zMat', []);
evalDiag.etaMat = localGetStructField(evalAux, 'etaMat', []);
evalDiag.blockValueMat = localGetStructField(evalAux, 'blockValue', []);
evalDiag.blockNorm2Mat = localGetStructField(evalAux, 'blockNorm2', []);
evalDiag.blockResidualMat = evalDiag.blockNorm2Mat - evalDiag.blockValueMat;
evalDiag.countPerBlock = localGetStructField(evalAux, 'countPerBlock', []);
evalDiag.blockAbsZMat = abs(evalDiag.zMat);
evalDiag.blockPhaseResidMat = localCalcBlockPhaseResid(evalDiag.zMat, evalAux.framePhase);
end


function coherenceSat = localCalcSatCoherence(zMat)
%LOCALCALCSATCOHERENCE Compute |sum(z)| / sum(|z|) per satellite.

if isempty(zMat)
  coherenceSat = [];
  return;
end
magSum = sum(abs(zMat), 2);
coherenceSat = nan(size(magSum));
validMask = magSum > 0;
coherenceSat(validMask) = abs(sum(zMat(validMask, :), 2)) ./ magSum(validMask);
coherenceSat = coherenceSat(:);
end


function blockPhaseResidMat = localCalcBlockPhaseResid(zMat, framePhaseMat)
%LOCALCALCBLOCKPHASERESID Compute angle(z)-phase for each valid block.

blockPhaseResidMat = nan(size(zMat));
if isempty(zMat) || isempty(framePhaseMat)
  return;
end
validMask = isfinite(zMat) & isfinite(framePhaseMat);
blockPhaseResidMat(validMask) = angle(zMat(validMask)) - framePhaseMat(validMask);
blockPhaseResidMat(validMask) = mod(blockPhaseResidMat(validMask) + pi, 2 * pi) - pi;
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
