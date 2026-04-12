function bestCase = selectBestStaticSeedCase(baseCase, weightCase, truth)
%SELECTBESTSTATICSEEDCASE Choose the most reliable static multi-sat seed.
% Dynamic MF refinement should start from the best resolved static MS case
% rather than always from the equal-weight branch.

arguments
  baseCase (1, 1) struct
  weightCase (1, :) struct
  truth
end

candidateCase = [baseCase, weightCase];
bestCase = baseCase;
bestScore = [inf, inf];

truthLatlon = reshape(getDoaDopplerFieldOrDefault(truth, 'latlonTrueDeg', nan(2, 1)), [], 1);
truthFdRef = localResolveTruthFdRef(truth);

for iCase = 1:numel(candidateCase)
  currentCase = candidateCase(iCase);
  if ~localCaseResolved(currentCase)
    continue;
  end

  estResult = currentCase.estResult;
  angleErrDeg = calcLatlonAngleError(estResult.doaParamEst(:), truthLatlon);
  fdErrHz = abs(estResult.fdRefEst - truthFdRef);
  scoreNow = [angleErrDeg, fdErrHz];
  if localLexicoLess(scoreNow, bestScore)
    bestScore = scoreNow;
    bestCase = currentCase;
  end
end
end


function truthFdRef = localResolveTruthFdRef(truth)
%LOCALRESOLVETRUTHFDREF Resolve reference Doppler truth with fallback.

truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefTrueHz', NaN);
if ~isfinite(truthFdRef)
  truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefFit', NaN);
end
end


function isResolved = localCaseResolved(caseInfo)
%LOCALCASERESOLVED Return true when one case contains a usable estimate.

isResolved = false;
if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
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


function isLess = localLexicoLess(scoreNow, scoreBest)
%LOCALLEXICOLESS Lexicographic compare for [angleErr, fdErr].

if isempty(scoreBest) || any(~isfinite(scoreBest))
  isLess = all(isfinite(scoreNow));
  return;
end
if any(~isfinite(scoreNow))
  isLess = false;
  return;
end
if scoreNow(1) < scoreBest(1) - 1e-12
  isLess = true;
  return;
end
if scoreNow(1) > scoreBest(1) + 1e-12
  isLess = false;
  return;
end
isLess = scoreNow(2) < scoreBest(2);
end
