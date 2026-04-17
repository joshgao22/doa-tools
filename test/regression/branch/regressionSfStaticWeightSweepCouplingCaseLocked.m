% Regression check for the locked SF static sat2-weight sweep after the
% DoA-anchor fallback update.
% The contract is now:
%   1) alpha=0 must remain a resolved near-ref-only anchor;
%   2) at least one positive sat2 weight must improve |fdRef error| relative
%      to the zero-weight anchor;
%   3) the best fd-improving positive-weight branch must not pay a material
%      angle penalty anymore. A small numerical drift is allowed, but the
%      old monotone angle degradation pattern should be suppressed.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticJointCouplingFixture();

fprintf('Running regressionSfStaticWeightSweepCouplingCaseLocked ...\n');

caseBundle = fixture.caseBundle;
refOnlyCase = caseBundle.caseStaticRefOnly;
weightCase = reshape(caseBundle.weightCase, 1, []);
alphaList = reshape(fixture.weightSweepAlpha, [], 1);
zeroIdx = localResolveZeroWeightIndex(alphaList);
zeroCase = weightCase(zeroIdx);

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
fdImproveHz = baseAbsFdErrHz - bestAbsFdErrHz;
angleDeltaDeg = angleErrDeg(bestFdIdx) - angleErrDeg(zeroIdx);

fdImproveTolHz = 1e-3;
angleSafetyTolDeg = max(1e-6, 0.1 * angleErrDeg(zeroIdx));
if fdImproveHz <= fdImproveTolHz
  error('regressionSfStaticWeightSweepCouplingCaseLocked:MissingFdTradeoff', ...
    ['The locked static sat2-weight sweep must expose a positive-weight branch ', ...
     'that improves |fdRef error| relative to the zero-weight anchor.']);
end
if angleDeltaDeg > angleSafetyTolDeg
  error('regressionSfStaticWeightSweepCouplingCaseLocked:ResidualAnglePenalty', ...
    ['After enabling the SF static DoA-anchor fallback, the best fd-improving ', ...
     'sat2-weight branch should not incur a material angle penalty relative ', ...
     'to the zero-weight anchor.']);
end

bestAngleIdx = localArgMin(angleErrDeg);
[bestOverallFdIdx, ~] = localArgMin(absFdRefErrHz);

fprintf('  alpha list                    : %s\n', localFormatNumericRow(alphaList));
fprintf('  angle err sweep (deg)         : %s\n', localFormatNumericRow(angleErrDeg));
fprintf('  fdRef err sweep (Hz)          : %s\n', localFormatNumericRow(fdRefErrHz));
fprintf('  objective sweep               : %s\n', localFormatNumericRow(fvalList));
fprintf('  zero-weight angle err (deg)   : %.6g\n', angleErrDeg(zeroIdx));
fprintf('  zero-weight fdRef err (Hz)    : %.6f\n', fdRefErrHz(zeroIdx));
fprintf('  ref-vs-W0 angle diff (deg)    : %.6g\n', zeroAngleDiffDeg);
fprintf('  ref-vs-W0 fdRef diff (Hz)     : %.6f\n', zeroFdDiffHz);
fprintf('  ref-vs-W0 fval diff           : %.6g\n', zeroFvalDiff);
fprintf('  best-angle alpha              : %.2f\n', alphaList(bestAngleIdx));
fprintf('  best-fd alpha                 : %.2f\n', alphaList(bestOverallFdIdx));
fprintf('  fd-improving alpha            : %.2f\n', alphaList(bestFdIdx));
fprintf('  fd improvement from W0 (Hz)   : %.6f\n', fdImproveHz);
fprintf('  angle delta from W0 (deg)     : %.6g\n', angleDeltaDeg);
fprintf('  angle safety tol (deg)        : %.6g\n', angleSafetyTolDeg);
fprintf('PASS: regressionSfStaticWeightSweepCouplingCaseLocked\n');

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

function text = localFormatNumericRow(value)
value = reshape(value, 1, []);
text = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%.6g', x), value, ...
  'UniformOutput', false), ', '));
end

function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
