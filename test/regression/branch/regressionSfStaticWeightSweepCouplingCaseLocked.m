function regressionSfStaticWeightSweepCouplingCaseLocked(varargin)
% Regression check for the locked SF static sat2-weight sweep after the
% DoA-anchor fallback update.
% The contract is now:
%   1) alpha=0 must remain a resolved near-ref-only anchor;
%   2) at least one positive sat2 weight must improve |fdRef error| relative
%      to the zero-weight anchor;
%   3) there must exist a positive-weight tradeoff branch that improves
%      |fdRef error| without paying a material angle penalty. The old
%      monotone pattern (fd improves only by sacrificing angle) should stay
%      suppressed, but the single most fd-improving branch is allowed to pay
%      a small extra angle cost.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
localPrint = @(varargin) fprintf(varargin{:});
fixture = buildSfStaticJointCouplingFixture();

localPrint('Running regressionSfStaticWeightSweepCouplingCaseLocked ...\n');

caseBundle = fixture.caseBundle;
refOnlyCase = caseBundle.caseStaticRefOnly;
weightCase = reshape(caseBundle.weightCase, 1, []);
alphaList = reshape(fixture.weightSweepAlpha, [], 1);
zeroIdx = localResolveZeroWeightIndex(alphaList);
numCase = numel(weightCase);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
absFdRefErrHz = nan(numCase, 1);
fvalList = nan(numCase, 1);
resolvedMask = false(numCase, 1);

for iCase = 1:numCase
  currentCase = weightCase(iCase);
  resolvedMask(iCase) = localCaseResolved(currentCase);
  if ~resolvedMask(iCase)
    continue;
  end

  estResult = currentCase.estResult;
  angleErrDeg(iCase) = calcLatlonAngleError(estResult.doaParamEst(:), fixture.truth.latlonTrueDeg(:));
  fdRefErrHz(iCase) = estResult.fdRefEst - fixture.truth.fdRefTrueHz;
  absFdRefErrHz(iCase) = abs(fdRefErrHz(iCase));
  fvalList(iCase) = getDoaDopplerFieldOrDefault(estResult, 'fval', NaN);
end

if ~all(resolvedMask)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:UnresolvedWeightCase', ...
    'All locked static weight-sweep cases must resolve before coupling checks.');
end

refAngleErrDeg = calcLatlonAngleError(refOnlyCase.estResult.doaParamEst(:), fixture.truth.latlonTrueDeg(:));
refFdRefErrHz = refOnlyCase.estResult.fdRefEst - fixture.truth.fdRefTrueHz;
refFval = getDoaDopplerFieldOrDefault(refOnlyCase.estResult, 'fval', NaN);

zeroAngleDiffDeg = abs(angleErrDeg(zeroIdx) - refAngleErrDeg);
zeroFdDiffHz = abs(fdRefErrHz(zeroIdx) - refFdRefErrHz);
zeroFvalDiff = abs(fvalList(zeroIdx) - refFval);

angleAnchorTolDeg = 1e-3;
fdAnchorTolHz = 5;
fvalAnchorTol = max(10, 1e-4 * max(1, abs(refFval)));
if zeroAngleDiffDeg > angleAnchorTolDeg
  error('regressionSfStaticWeightSweepCouplingCaseLocked:ZeroWeightAngleDrift', ...
    ['The zero-weight MS static branch drifted too far from the ref-only ', ...
     'static angle anchor.']);
end
if zeroFdDiffHz > fdAnchorTolHz
  error('regressionSfStaticWeightSweepCouplingCaseLocked:ZeroWeightFdDrift', ...
    ['The zero-weight MS static branch drifted too far from the ref-only ', ...
     'static fdRef anchor.']);
end
if zeroFvalDiff > fvalAnchorTol
  error('regressionSfStaticWeightSweepCouplingCaseLocked:ZeroWeightObjectiveDrift', ...
    ['The zero-weight MS static branch objective drifted too far from the ', ...
     'ref-only static anchor.']);
end

positiveIdx = find(alphaList > 0);
if isempty(positiveIdx)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:MissingPositiveWeight', ...
    'The locked weight sweep must include at least one positive sat2 weight.');
end

[bestAbsFdErrHz, bestFdPosLocalIdx] = min(absFdRefErrHz(positiveIdx));
bestFdIdx = positiveIdx(bestFdPosLocalIdx);
baseAbsFdErrHz = absFdRefErrHz(zeroIdx);
fdImproveVecHz = baseAbsFdErrHz - absFdRefErrHz(positiveIdx);
angleDeltaVecDeg = angleErrDeg(positiveIdx) - angleErrDeg(zeroIdx);
fdImproveHz = fdImproveVecHz(bestFdPosLocalIdx);
angleDeltaDeg = angleDeltaVecDeg(bestFdPosLocalIdx);

fdImproveTolHz = 1e-3;
angleSafetyTolDeg = max(1e-6, 0.1 * angleErrDeg(zeroIdx));
tradeoffMask = fdImproveVecHz > fdImproveTolHz;
safeTradeoffMask = tradeoffMask & (angleDeltaVecDeg <= angleSafetyTolDeg);

tradeoffScore = [angleErrDeg(positiveIdx), absFdRefErrHz(positiveIdx)];
tradeoffScore(~safeTradeoffMask, :) = inf;
bestSafeIdx = NaN;
bestSafeFdImproveHz = NaN;
bestSafeAngleDeltaDeg = NaN;
if any(safeTradeoffMask)
  [bestSafeLocalIdx, ~] = localArgMinLexicographic(tradeoffScore);
  bestSafeIdx = positiveIdx(bestSafeLocalIdx);
  bestSafeFdImproveHz = baseAbsFdErrHz - absFdRefErrHz(bestSafeIdx);
  bestSafeAngleDeltaDeg = angleErrDeg(bestSafeIdx) - angleErrDeg(zeroIdx);
end

bestAngleIdx = localArgMin(angleErrDeg);
[bestOverallFdIdx, ~] = localArgMin(absFdRefErrHz);

localPrint('  alpha list                    : %s\n', localFormatNumericRow(alphaList));
localPrint('  angle err sweep (deg)         : %s\n', localFormatNumericRow(angleErrDeg));
localPrint('  fdRef err sweep (Hz)          : %s\n', localFormatNumericRow(fdRefErrHz));
localPrint('  objective sweep               : %s\n', localFormatNumericRow(fvalList));
localPrint('  zero-weight angle err (deg)   : %.6g\n', angleErrDeg(zeroIdx));
localPrint('  zero-weight fdRef err (Hz)    : %.6f\n', fdRefErrHz(zeroIdx));
localPrint('  ref-vs-W0 angle diff (deg)    : %.6g\n', zeroAngleDiffDeg);
localPrint('  ref-vs-W0 fdRef diff (Hz)     : %.6f\n', zeroFdDiffHz);
localPrint('  ref-vs-W0 fval diff           : %.6g\n', zeroFvalDiff);
localPrint('  best-angle alpha              : %.2f\n', alphaList(bestAngleIdx));
localPrint('  best-fd alpha                 : %.2f\n', alphaList(bestOverallFdIdx));
localPrint('  strongest-fd alpha            : %.2f\n', alphaList(bestFdIdx));
localPrint('  strongest-fd improve (Hz)     : %.6f\n', fdImproveHz);
localPrint('  strongest-fd angle delta (deg): %.6g\n', angleDeltaDeg);
if isnan(bestSafeIdx)
  bestSafeAlpha = NaN;
else
  bestSafeAlpha = alphaList(bestSafeIdx);
end
localPrint('  best-safe alpha               : %.2f\n', bestSafeAlpha);
localPrint('  best-safe fd improve (Hz)     : %.6f\n', bestSafeFdImproveHz);
localPrint('  best-safe angle delta (deg)   : %.6g\n', bestSafeAngleDeltaDeg);
localPrint('  angle safety tol (deg)        : %.6g\n', angleSafetyTolDeg);

if ~any(tradeoffMask)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:MissingFdTradeoff', ...
    ['The locked static sat2-weight sweep must expose a positive-weight branch ', ...
     'that improves |fdRef error| relative to the zero-weight anchor.']);
end
if ~any(safeTradeoffMask)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:ResidualAnglePenalty', ...
    ['After enabling the SF static DoA-anchor fallback, the locked static ', ...
     'sat2-weight sweep must contain at least one fd-improving positive-weight ', ...
     'branch whose angle penalty stays within the material-safety tolerance.']);
end

localPrint('PASS: regressionSfStaticWeightSweepCouplingCaseLocked\n');

end

function idx = localResolveZeroWeightIndex(alphaList)
idx = find(abs(alphaList) <= eps(max(1, max(abs(alphaList)))), 1, 'first');
if isempty(idx)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:MissingZeroWeight', ...
    'The locked static weight sweep must include alpha = 0.');
end
end

function tf = localCaseResolved(caseInfo)
tf = false;
if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end
estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  tf = logical(estResult.isResolved);
else
  tf = true;
end
if ~tf
  return;
end

tf = isfield(estResult, 'doaParamEst') && numel(estResult.doaParamEst) == 2 && ...
  all(isfinite(estResult.doaParamEst(:))) && isfield(estResult, 'fdRefEst') && ...
  isscalar(estResult.fdRefEst) && isfinite(estResult.fdRefEst);
end

function [idx, value] = localArgMin(metric)
[value, idx] = min(metric);
if isempty(idx) || ~isfinite(value)
  error('regressionSfStaticWeightSweepCouplingCaseLocked:NonFiniteMetric', ...
    'The locked static weight-sweep metric must stay finite.');
end
end

function [idx, scoreRow] = localArgMinLexicographic(scoreMat)
scoreMat = reshape(scoreMat, size(scoreMat, 1), []);
[~, order] = sortrows(scoreMat);
idx = order(1);
scoreRow = scoreMat(idx, :);
if any(~isfinite(scoreRow))
  error('regressionSfStaticWeightSweepCouplingCaseLocked:NoFiniteTradeoff', ...
    'The locked static weight sweep must contain one finite safe tradeoff branch.');
end
end

function text = localFormatNumericRow(value)
value = reshape(value, 1, []);
text = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%.6g', x), value, ...
  'UniformOutput', false), ', '));
end
