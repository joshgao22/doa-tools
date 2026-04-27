function regressionSfStaticReferenceDopplerInvariant(varargin)
% Regression check for SF static reference-Doppler invariants.
% This script focuses on one narrow contract only:
%   1) every resolved SF static estimator result must satisfy
%         deltaFdRef(refSat) = 0
%         fdSat(refSat) = fdRef
%         fdSat = fdRef + deltaFdRef
%      in the estimator output aux fields;
%   2) the ref-only subset must remap refSatIdxLocal to 1, while the full
%      multi-satellite cases must preserve the full-scene reference index.

opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = buildSfStaticRegressionFixture();

  fprintf('Running regressionSfStaticReferenceDopplerInvariant ...\n');

caseList = [ ...
  fixture.caseBundle.caseStaticRefOnly, ...
  fixture.caseStaticZeroWeightRefInit, ...
  fixture.caseBundle.caseStaticMs, ...
  fixture.caseBundle.weightCase ...
  ];
expectedRefIdxList = [ ...
  1, ...
  fixture.refSatIdxLocal, ...
  fixture.refSatIdxLocal, ...
  repmat(fixture.refSatIdxLocal, 1, numel(fixture.caseBundle.weightCase)) ...
  ];
expectedNumSatList = [ ...
  1, ...
  fixture.scene.numSat, ...
  fixture.scene.numSat, ...
  repmat(fixture.scene.numSat, 1, numel(fixture.caseBundle.weightCase)) ...
  ];

maxComposeErrHz = 0;
maxRefMatchErrHz = 0;
maxDeltaRefErrHz = 0;
for iCase = 1:numel(caseList)
  [deltaRefErrHz, refMatchErrHz, composeErrHz] = localCheckCaseInvariant( ...
    caseList(iCase), expectedRefIdxList(iCase), expectedNumSatList(iCase), verbose);
  maxDeltaRefErrHz = max(maxDeltaRefErrHz, deltaRefErrHz);
  maxRefMatchErrHz = max(maxRefMatchErrHz, refMatchErrHz);
  maxComposeErrHz = max(maxComposeErrHz, composeErrHz);
end

  fprintf('  checked case count           : %d\n', numel(caseList));
  fprintf('  max |deltaFdRef(refSat)| Hz  : %.6g\n', maxDeltaRefErrHz);
  fprintf('  max |fdSat(refSat)-fdRef| Hz : %.6g\n', maxRefMatchErrHz);
  fprintf('  max composition error Hz     : %.6g\n', maxComposeErrHz);
  fprintf('PASS: regressionSfStaticReferenceDopplerInvariant\n');


end

function [deltaRefErrHz, refMatchErrHz, composeErrHz] = localCheckCaseInvariant(caseInfo, expectedRefIdxLocal, expectedNumSat, verbose)
%LOCALCHECKCASEINVARIANT Check one resolved SF static case.

if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  error('regressionSfStaticReferenceDopplerInvariant:MissingEstimate', ...
    'Each checked case must contain one estimator result.');
end

estResult = caseInfo.estResult;
if ~logical(getDoaDopplerFieldOrDefault(estResult, 'isResolved', false))
  caseTag = char(string(caseInfo.displayName));
  error('regressionSfStaticReferenceDopplerInvariant:UnresolvedCase', ...
    'Case %s must resolve before checking reference-Doppler invariants.', ...
    caseTag);
end

caseTag = char(string(caseInfo.displayName));
aux = getDoaDopplerFieldOrDefault(estResult, 'aux', struct());
fdSat = reshape(getDoaDopplerFieldOrDefault(aux, 'fdSatEst', []), [], 1);
deltaFd = reshape(getDoaDopplerFieldOrDefault(aux, 'deltaFdRefEst', []), [], 1);
refSatIdxLocal = getDoaDopplerFieldOrDefault(aux, 'refSatIdxLocal', NaN);
fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);

if numel(fdSat) ~= expectedNumSat || numel(deltaFd) ~= expectedNumSat
  error('regressionSfStaticReferenceDopplerInvariant:InvalidSatLength', ...
    ['Case %s must return fdSatEst and deltaFdRefEst with %d entries. ', ...
     'Got %d and %d.'], ...
    caseTag, expectedNumSat, numel(fdSat), numel(deltaFd));
end
if ~isscalar(refSatIdxLocal) || refSatIdxLocal ~= expectedRefIdxLocal
  error('regressionSfStaticReferenceDopplerInvariant:InvalidReferenceIndex', ...
    'Case %s must keep refSatIdxLocal = %d, got %g.', ...
    caseTag, expectedRefIdxLocal, refSatIdxLocal);
end
if ~isscalar(fdRefEst) || ~isfinite(fdRefEst)
  error('regressionSfStaticReferenceDopplerInvariant:InvalidFdRef', ...
    'Case %s must return one finite fdRefEst.', caseTag);
end

composeVec = fdRefEst + deltaFd;
deltaRefErrHz = abs(deltaFd(refSatIdxLocal));
refMatchErrHz = abs(fdSat(refSatIdxLocal) - fdRefEst);
composeErrHz = max(abs(fdSat - composeVec));

tolHz = 1e-9 * max([1; abs(fdSat(:)); abs(fdRefEst)]);
if deltaRefErrHz > tolHz
  error('regressionSfStaticReferenceDopplerInvariant:ReferenceDeltaMismatch', ...
    'Case %s violates deltaFdRef(refSat)=0 by %.3e Hz (tol %.3e Hz).', ...
    caseTag, deltaRefErrHz, tolHz);
end
if refMatchErrHz > tolHz
  error('regressionSfStaticReferenceDopplerInvariant:ReferenceFdMismatch', ...
    'Case %s violates fdSat(refSat)=fdRef by %.3e Hz (tol %.3e Hz).', ...
    caseTag, refMatchErrHz, tolHz);
end
if composeErrHz > tolHz
  error('regressionSfStaticReferenceDopplerInvariant:FdCompositionMismatch', ...
    'Case %s violates fdSat=fdRef+deltaFdRef by %.3e Hz (tol %.3e Hz).', ...
    caseTag, composeErrHz, tolHz);
end

if verbose
  fprintf('  %-28s refIdx=%d numSat=%d composeErr=%.3e Hz\n', ...
    caseTag, refSatIdxLocal, numel(fdSat), composeErrHz);
end
end
