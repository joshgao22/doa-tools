function initParam = buildDynamicInitParamFromCase(caseInfo, isKnownRate, fdRateSeed)
%BUILDDYNAMICINITPARAMFROMCASE Build one MF initializer from one case.

arguments
  caseInfo (1, 1) struct
  isKnownRate (1, 1) logical
  fdRateSeed (1, 1) double = NaN
end

initParam = [];
if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
if ~localCaseResolved(caseInfo) || ~isfield(estResult, 'doaParamEst') || ...
    ~isfield(estResult, 'fdRefEst')
  return;
end

doaParam = reshape(estResult.doaParamEst, [], 1);
fdRef = estResult.fdRefEst;
if ~isnumeric(doaParam) || numel(doaParam) ~= 2 || any(~isfinite(doaParam)) || ...
    ~isscalar(fdRef) || ~isfinite(fdRef)
  return;
end

if isKnownRate
  initParam = [doaParam; fdRef];
  return;
end

if ~isfinite(fdRateSeed)
  fdRateEst = getDoaDopplerFieldOrDefault(estResult, 'fdRateEst', NaN);
  if isscalar(fdRateEst) && isfinite(fdRateEst)
    fdRateSeed = fdRateEst;
  else
    return;
  end
end
initParam = [doaParam; fdRef; fdRateSeed];
end


function isResolved = localCaseResolved(caseInfo)
%LOCALCASERESOLVED Return true when one case contains a usable estimate.

isResolved = false;
if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  isResolved = logical(estResult.isResolved);
else
  isResolved = true;
end
if ~isResolved
  return;
end

if ~isfield(estResult, 'doaParamEst') || numel(estResult.doaParamEst) ~= 2 || ...
    any(~isfinite(estResult.doaParamEst(:)))
  isResolved = false;
  return;
end
if ~isfield(estResult, 'fdRefEst') || ~isscalar(estResult.fdRefEst) || ~isfinite(estResult.fdRefEst)
  isResolved = false;
end
end
