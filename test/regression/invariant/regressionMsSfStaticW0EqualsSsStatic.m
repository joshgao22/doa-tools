% Regression check for the MS-SF-Static zero-weight anchor.
% Contract:
%   1) MS-SF-Static-W0.00 must match SS-SF-Static when both use the same
%      ref-sat DoA initializer;
%   2) the same zero-weight case must also match the ref-only ablation
%      branch kept in the shared transition bundle.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionMsSfStaticW0EqualsSsStatic ...\n');

caseBundle = fixture.caseBundle;
caseRefOnly = caseBundle.caseStaticRefOnly;
caseRefAbl = caseBundle.caseStaticRefAbl;
caseWeightZero = fixture.caseStaticZeroWeightRefInit;

localAssertResolved(caseRefOnly, 'SS-SF-Static');
localAssertResolved(caseRefAbl, 'MS-SF-Static-refSat');
localAssertResolved(caseWeightZero, 'MS-SF-Static-W0.00-RefInit');

truthLatlon = fixture.truth.latlonTrueDeg(:);
angleErrRefDeg = calcLatlonAngleError(caseRefOnly.estResult.doaParamEst(:), truthLatlon);
angleErrAblDeg = calcLatlonAngleError(caseRefAbl.estResult.doaParamEst(:), truthLatlon);
angleErrW0Deg = calcLatlonAngleError(caseWeightZero.estResult.doaParamEst(:), truthLatlon);

angleDiffRefVsW0Deg = calcLatlonAngleError( ...
  caseRefOnly.estResult.doaParamEst(:), caseWeightZero.estResult.doaParamEst(:));
angleDiffAblVsW0Deg = calcLatlonAngleError( ...
  caseRefAbl.estResult.doaParamEst(:), caseWeightZero.estResult.doaParamEst(:));
fdDiffRefVsW0Hz = abs(caseRefOnly.estResult.fdRefEst - caseWeightZero.estResult.fdRefEst);
fdDiffAblVsW0Hz = abs(caseRefAbl.estResult.fdRefEst - caseWeightZero.estResult.fdRefEst);
fvalDiffRefVsW0 = abs(localGetFval(caseRefOnly) - localGetFval(caseWeightZero));
fvalDiffAblVsW0 = abs(localGetFval(caseRefAbl) - localGetFval(caseWeightZero));
tolAngleDeg = 1e-6;
tolFdRefHz = 1e-3;
tolFval = 1e-3;

if angleDiffRefVsW0Deg > tolAngleDeg
  error('regressionMsSfStaticW0EqualsSsStatic:RefOnlyMismatch', ...
    'MS-SF-Static-W0.00 must match SS-SF-Static DoA under the same ref-only initializer.');
end
if angleDiffAblVsW0Deg > tolAngleDeg
  error('regressionMsSfStaticW0EqualsSsStatic:RefAblMismatch', ...
    'MS-SF-Static-W0.00 must match the ref-sat ablation DoA exactly.');
end
if fdDiffRefVsW0Hz > tolFdRefHz || fdDiffAblVsW0Hz > tolFdRefHz
  error('regressionMsSfStaticW0EqualsSsStatic:FdRefMismatch', ...
    'MS-SF-Static-W0.00 must preserve the ref-only fdRef estimate exactly.');
end
if fvalDiffRefVsW0 > tolFval || fvalDiffAblVsW0 > tolFval
  error('regressionMsSfStaticW0EqualsSsStatic:FvalMismatch', ...
    'MS-SF-Static-W0.00 must preserve the ref-only objective value exactly.');
end

fprintf('  SS-SF-Static angle err (deg)      : %.6g\n', angleErrRefDeg);
fprintf('  ref-sat ablation angle err (deg)  : %.6g\n', angleErrAblDeg);
fprintf('  W0 angle err (deg)                : %.6g\n', angleErrW0Deg);
fprintf('  ref-vs-W0 angle diff (deg)        : %.6g\n', angleDiffRefVsW0Deg);
fprintf('  ref-vs-W0 fdRef diff (Hz)         : %.6g\n', fdDiffRefVsW0Hz);
fprintf('  ref-vs-W0 objective diff          : %.6g\n', fvalDiffRefVsW0);
fprintf('PASS: regressionMsSfStaticW0EqualsSsStatic\n');


function localAssertResolved(caseInfo, displayName)
%LOCALASSERTRESOLVED Ensure one regression case resolved successfully.

if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  error('regressionMsSfStaticW0EqualsSsStatic:MissingCase', ...
    'Missing regression case: %s.', displayName);
end
estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  if ~logical(estResult.isResolved)
    error('regressionMsSfStaticW0EqualsSsStatic:UnresolvedCase', ...
      'Regression case did not resolve: %s.', displayName);
  end
end
end

function fval = localGetFval(caseInfo)
%LOCALGETFVAL Resolve one objective value with default fallback.

fval = NaN;
if isstruct(caseInfo) && isfield(caseInfo, 'estResult') && isstruct(caseInfo.estResult)
  fval = getDoaDopplerFieldOrDefault(caseInfo.estResult, 'fval', NaN);
end
end

function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
