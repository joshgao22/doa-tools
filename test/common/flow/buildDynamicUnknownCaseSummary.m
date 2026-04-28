function summary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth)
%BUILDDYNAMICUNKNOWNCASESUMMARY Build one compact dynamic-case summary.
% This helper keeps summary construction reusable outside the main bundle.
% The truth struct is optional: when omitted, the summary remains decision-
% safe and only contains internal objective / health metrics. Truth-relative
% fields are evaluation-only aliases and must not be used by flow selectors.

arguments
  caseUse (1,1) struct
  toothStepHz (1,1) double = NaN
  truth (1,1) struct = struct()
end

summary = localBuildSummarySkeleton();
summary.stageTag = localInferStageTag(caseUse);
summary.startTag = string(localGetFieldOrDefault(caseUse, 'startTag', ""));
summary.routeFamily = string(localGetFieldOrDefault(caseUse, 'routeFamily', ""));

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult) || ~isstruct(estResult)
  return;
end

summary.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', summary.solveVariant));
summary.isResolved = logical(localGetFieldOrDefault(estResult, 'isResolved', false));
summary.doaParamEst = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', nan(2, 1)), 1, []);
summary.fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
summary.fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
summary.runTimeMs = localGetFieldOrDefault(estResult, 'runTimeMs', NaN);

if strlength(summary.stageTag) == 0
  summary.stageTag = string(localGetFieldOrDefault(estResult, 'stageTag', summary.solveVariant));
end

if strlength(summary.startTag) == 0
  summary.startTag = string(localGetFieldOrDefault(estResult, 'startTag', ""));
end
if strlength(summary.routeFamily) == 0
  summary.routeFamily = string(localGetFieldOrDefault(estResult, 'routeFamily', ""));
end

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
summary.isDoaFrozenLike = localIsFrozenLikeSummary(summary.solveVariant, summary.stageTag);
[summary.refCoherence, summary.nonRefCoherenceFloor, summary.nonRefMaxAbsPhaseResidRad, ...
  summary.nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval);

if isempty(fieldnames(truth))
  return;
end

truthLatlon = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
if numel(truthLatlon) >= 2 && numel(summary.doaParamEst) >= 2 && all(isfinite(summary.doaParamEst(1:2)))
  summary.truthAngleErrDeg = calcLatlonAngleError(summary.doaParamEst(1:2).', truthLatlon(1:2));
  summary.angleErrDeg = summary.truthAngleErrDeg;
end
truthFdRef = localGetFieldOrDefault(truth, 'fdRefTrueHz', localGetFieldOrDefault(truth, 'fdRefFit', NaN));
truthFdRate = localGetFieldOrDefault(truth, 'fdRateTrueHzPerSec', localGetFieldOrDefault(truth, 'fdRateFit', NaN));
summary.truthFdRefErrHz = summary.fdRefEst - truthFdRef;
summary.truthFdRateErrHzPerSec = summary.fdRateEst - truthFdRate;
summary.fdRefErrHz = summary.truthFdRefErrHz;
summary.fdRateErrHzPerSec = summary.truthFdRateErrHzPerSec;
if isfinite(toothStepHz) && toothStepHz > 0 && isfinite(summary.truthFdRefErrHz)
  summary.truthToothIdx = round(summary.truthFdRefErrHz / toothStepHz);
  summary.truthToothResidualHz = summary.truthFdRefErrHz - summary.truthToothIdx * toothStepHz;
  summary.toothIdx = summary.truthToothIdx;
  summary.toothResidualHz = summary.truthToothResidualHz;
end
end


function summary = localBuildSummarySkeleton()
%LOCALBUILDSUMMARYSKELETON Build one stable default summary shape.

summary = struct();
summary.solveVariant = "skipped";
summary.isResolved = false;
summary.doaParamEst = nan(1, 2);
summary.fdRefEst = NaN;
summary.fdRateEst = NaN;
summary.runTimeMs = NaN;
summary.finalObj = NaN;
summary.finalResidualNorm = NaN;
summary.refFitRatio = NaN;
summary.nonRefFitRatioFloor = NaN;
summary.refSupportRatio = NaN;
summary.nonRefSupportRatioFloor = NaN;
summary.maxNonRefNegativeProjectionRatio = NaN;
summary.nonRefFitFloorPenalty = 0;
summary.nonRefSupportFloorPenalty = 0;
summary.additionalObjectivePenalty = 0;
summary.refConsistencyNorm = NaN;
summary.nonRefConsistencyRatioFloor = NaN;
summary.isDoaFrozenLike = false;
summary.refCoherence = NaN;
summary.nonRefCoherenceFloor = NaN;
summary.nonRefMaxAbsPhaseResidRad = NaN;
summary.nonRefRmsPhaseResidRad = NaN;
summary.angleErrDeg = NaN;
summary.truthAngleErrDeg = NaN;
summary.fdRefErrHz = NaN;
summary.fdRateErrHzPerSec = NaN;
summary.toothIdx = NaN;
summary.toothResidualHz = NaN;
summary.truthFdRefErrHz = NaN;
summary.truthFdRateErrHzPerSec = NaN;
summary.truthToothIdx = NaN;
summary.truthToothResidualHz = NaN;
summary.stageTag = "";
summary.startTag = "";
summary.routeFamily = "";
end


function stageTag = localInferStageTag(caseUse)
%LOCALINFERSTAGETAG Infer one stable stage tag from a case or summary struct.

stageTag = "";
if ~isstruct(caseUse)
  return;
end
stageTag = string(localGetFieldOrDefault(caseUse, 'stageTag', stageTag));
if strlength(stageTag) > 0
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
stageTag = string(localGetFieldOrDefault(estResult, 'stageTag', stageTag));
if strlength(stageTag) > 0
  return;
end
stageTag = string(localGetFieldOrDefault(estResult, 'solveVariant', stageTag));
end


function tf = localIsFrozenLikeSummary(solveVariant, stageTag)
%LOCALISFROZENLIKESUMMARY Return true when one summary looks fixed-DoA.

solveVariant = lower(string(solveVariant));
stageTag = lower(string(stageTag));
tf = contains(solveVariant, "fixeddoa") || ...
  (contains(stageTag, "anchor") && ~contains(stageTag, "polish") && ~contains(stageTag, "wide"));
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
