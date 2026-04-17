% Regression check for SF static ref-only equivalence.
% This script focuses on one narrow contract only:
%   1) the dedicated SS-SF-Static case and the multi-sat static case with
%      sat2 weight = 0 AND the same DoA initializer must produce the same
%      resolved estimate;
%   2) adding the second satellite with zero weight must not perturb the
%      reference-only static objective value when the initialization is held
%      fixed.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionSfStaticRefOnlyEquivalence ...\n');

caseRefOnly = fixture.caseBundle.caseStaticRefOnly;
caseWeightZero = fixture.caseStaticZeroWeightRefInit;

estRefOnly = caseRefOnly.estResult;
estWeightZero = caseWeightZero.estResult;

if ~logical(estRefOnly.isResolved) || ~logical(estWeightZero.isResolved)
  error('regressionSfStaticRefOnlyEquivalence:UnresolvedCase', ...
    'Both SS-SF-Static and zero-weight MS-SF-Static must resolve successfully.');
end

angleDiffDeg = calcLatlonAngleError(estRefOnly.doaParamEst(:), estWeightZero.doaParamEst(:));
fdRefDiffHz = abs(estRefOnly.fdRefEst - estWeightZero.fdRefEst);
fvalDiff = abs(estRefOnly.fval - estWeightZero.fval);

if angleDiffDeg > 1e-9
  error('regressionSfStaticRefOnlyEquivalence:DoaMismatch', ...
    'Zero-weight multi-sat static must match ref-only static DoA when the same DoA initializer is used.');
end
if fdRefDiffHz > 1e-6
  error('regressionSfStaticRefOnlyEquivalence:FdRefMismatch', ...
    'Zero-weight multi-sat static must match ref-only static fdRef when the same DoA initializer is used.');
end
if fvalDiff > 1e-6
  error('regressionSfStaticRefOnlyEquivalence:FvalMismatch', ...
    'Zero-weight multi-sat static must match ref-only static objective value when the same DoA initializer is used.');
end
if estRefOnly.exitflag ~= estWeightZero.exitflag
  error('regressionSfStaticRefOnlyEquivalence:ExitflagMismatch', ...
    'Zero-weight multi-sat static must keep the same exitflag as ref-only static.');
end

fprintf('  ref-only angle err (deg) : %.6g\n', ...
  calcLatlonAngleError(estRefOnly.doaParamEst(:), fixture.truth.latlonTrueDeg(:)));
fprintf('  zero-weight angle err    : %.6g\n', ...
  calcLatlonAngleError(estWeightZero.doaParamEst(:), fixture.truth.latlonTrueDeg(:)));
fprintf('  fdRef diff (Hz)          : %.6g\n', fdRefDiffHz);
fprintf('  fval diff                : %.6g\n', fvalDiff);
fprintf('PASS: regressionSfStaticRefOnlyEquivalence\n');


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
