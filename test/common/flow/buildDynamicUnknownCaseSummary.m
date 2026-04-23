function summary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth)
%BUILDDYNAMICUNKNOWNCASESUMMARY Build one compact dynamic-case summary.
% This helper keeps summary construction reusable outside the main bundle.
% The truth struct is optional: when omitted, the summary remains selection-
% safe and only contains internal objective / health metrics.

arguments
  caseUse (1,1) struct
  toothStepHz (1,1) double = NaN
  truth (1,1) struct = struct()
end

summary = struct();
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
summary.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(estResult, 'isResolved', false));
summary.doaParamEst = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', nan(2, 1)), 1, []);
summary.fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
summary.fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
summary.runTimeMs = localGetFieldOrDefault(estResult, 'runTimeMs', NaN);

debugAux = localGetFieldOrDefault(localGetFieldOrDefault(estResult, 'aux', struct()), 'debug', struct());
finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
summary.finalObj = localGetFieldOrDefault(finalEval, 'obj', NaN);
summary.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
summary.refFitRatio = localGetFieldOrDefault(finalEval, 'refFitRatio', NaN);
summary.nonRefFitRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefFitRatioFloor', NaN);
summary.refSupportRatio = localGetFieldOrDefault(finalEval, 'refSupportRatio', NaN);
summary.nonRefSupportRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefSupportRatioFloor', NaN);
summary.maxNonRefNegativeProjectionRatio = localGetFieldOrDefault(finalEval, 'maxNonRefNegativeProjectionRatio', NaN);
summary.nonRefFitFloorPenalty = localGetFieldOrDefault(finalEval, 'nonRefFitFloorPenalty', 0);
summary.nonRefSupportFloorPenalty = localGetFieldOrDefault(finalEval, 'nonRefSupportFloorPenalty', 0);
summary.additionalObjectivePenalty = localGetFieldOrDefault(finalEval, 'additionalObjectivePenalty', 0);
summary.refConsistencyNorm = localGetFieldOrDefault(finalEval, 'refConsistencyNorm', NaN);
summary.nonRefConsistencyRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefConsistencyRatioFloor', NaN);
summary.isDoaFrozenLike = contains(lower(strtrim(summary.solveVariant)), "fixeddoa");
[summary.refCoherence, summary.nonRefCoherenceFloor, summary.nonRefMaxAbsPhaseResidRad, ...
  summary.nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval);

summary.angleErrDeg = NaN;
summary.fdRefErrHz = NaN;
summary.fdRateErrHzPerSec = NaN;
summary.toothIdx = NaN;
summary.toothResidualHz = NaN;
if ~isempty(fieldnames(truth))
  truthLatlon = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
  if numel(truthLatlon) >= 2 && numel(summary.doaParamEst) >= 2 && all(isfinite(summary.doaParamEst(1:2)))
    summary.angleErrDeg = calcLatlonAngleError(summary.doaParamEst(1:2).', truthLatlon(1:2));
  end
  truthFdRef = localGetFieldOrDefault(truth, 'fdRefTrueHz', localGetFieldOrDefault(truth, 'fdRefFit', NaN));
  truthFdRate = localGetFieldOrDefault(truth, 'fdRateTrueHzPerSec', localGetFieldOrDefault(truth, 'fdRateFit', NaN));
  summary.fdRefErrHz = summary.fdRefEst - truthFdRef;
  summary.fdRateErrHzPerSec = summary.fdRateEst - truthFdRate;
  if isfinite(toothStepHz) && toothStepHz > 0 && isfinite(summary.fdRefErrHz)
    summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
    summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
  end
end
end

function [refCoherence, nonRefCoherenceFloor, nonRefMaxAbsPhaseResidRad, nonRefRmsPhaseResidRad] = ...
  localExtractNonRefEvalMetrics(finalEval)
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

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
