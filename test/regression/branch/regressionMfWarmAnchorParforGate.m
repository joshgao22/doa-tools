function regressionMfWarmAnchorParforGate(varargin)
% Regression check for the warm-anchor parfor gate.
% This script verifies one branch-level contract only: the estimator default
% path, verbose trace path, fixed-DoA tooth-guard path, and too-small release
% family must all stay serial.  The opt-in non-frozen parfor path remains
% environment-dependent and is not required to be available in quick tests.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fprintf('Running regressionMfWarmAnchorParforGate ...\n');

baseModel = struct();
baseModel.unknownWarmAnchorUseParfor = false;
baseModel.unknownWarmAnchorMinParforSeed = 3;
baseModel.freezeDoa = false;
localCheckFalse(baseModel, 8, false, 'default opt-out');

optInModel = baseModel;
optInModel.unknownWarmAnchorUseParfor = true;
localCheckFalse(optInModel, 8, true, 'verbose trace');

fixedDoaModel = optInModel;
fixedDoaModel.freezeDoa = true;
localCheckFalse(fixedDoaModel, 8, false, 'fixed-DoA tooth guard');

smallFamilyModel = optInModel;
smallFamilyModel.unknownWarmAnchorMinParforSeed = 5;
localCheckFalse(smallFamilyModel, 4, false, 'release family too small');

maybeParfor = shouldUseDoaDopplerMfWarmAnchorParfor(optInModel, 8, false);
if ~(islogical(maybeParfor) && isscalar(maybeParfor))
  error('regressionMfWarmAnchorParforGate:InvalidOptInReturn', ...
    'The opt-in warm-anchor parfor gate must return one logical scalar.');
end

fprintf('  default opt-out gate       : serial\n');
fprintf('  verbose trace gate         : serial\n');
fprintf('  fixed-DoA tooth guard gate : serial\n');
fprintf('  small release family gate  : serial\n');
fprintf('  opt-in environment result  : %d\n', maybeParfor);
fprintf('PASS: regressionMfWarmAnchorParforGate\n');
end

function localCheckFalse(model, numReleaseSeed, verboseFlag, caseName)
%LOCALCHECKFALSE Assert that one gate case stays serial.

useParfor = shouldUseDoaDopplerMfWarmAnchorParfor(model, numReleaseSeed, verboseFlag);
if useParfor
  error('regressionMfWarmAnchorParforGate:UnexpectedParfor', ...
    'Warm-anchor parfor must stay disabled for case: %s.', caseName);
end
end
