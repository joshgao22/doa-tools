function [angleErrDeg, fdErrHz, fdRateErrHzPerSec, isResolved] = extractDynamicCaseMetric(caseInfo, truth)
%EXTRACTDYNAMICCASEMETRIC Extract one compact metric tuple from one dynamic case.
% This helper keeps dev, perf, and smoke regression summaries on the same
% error convention without touching estimator result fields.

angleErrDeg = NaN;
fdErrHz = NaN;
fdRateErrHzPerSec = NaN;
isResolved = false;

if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
latlonEst = getDoaDopplerLatlonEst(estResult);
if all(isfinite(latlonEst(1:2)))
  angleErrDeg = calcLatlonAngleError(latlonEst(1:2), truth.latlonTrueDeg(:));
end

fdRefTruth = truth.fdRefTrueHz;
if isfield(caseInfo, 'dynamicMode') && startsWith(string(caseInfo.dynamicMode), "cp-")
  fdRefTruth = truth.fdRefFit;
end
fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
if isfinite(fdRefEst) && isfinite(fdRefTruth)
  fdErrHz = abs(fdRefEst - fdRefTruth);
end

fdRateTruth = NaN;
if isfield(caseInfo, 'dynamicMode') && startsWith(string(caseInfo.dynamicMode), "cp-")
  fdRateTruth = truth.fdRateFit;
end
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
if isfinite(fdRateEst) && isfinite(fdRateTruth)
  fdRateErrHzPerSec = abs(fdRateEst - fdRateTruth);
end

if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  isResolved = logical(estResult.isResolved);
else
  isResolved = isfinite(angleErrDeg);
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one nonempty field with default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end
