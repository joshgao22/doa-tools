% Regression check for SF static zero-weight fdRef initializer invariance.
% This script guards one narrow contract only:
%   1) with the same fixed DoA anchor, the multi-sat static fdRef initializer
%      under satWeight = [1; 0] must reduce to the ref-only static fdRef
%      initializer;
%   2) the zero-weight reduction must hold at the initializer stage as well,
%      not only after the final free solve.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionSfStaticFdInitZeroWeightInvariant ...\n');

msDoaAnchor = fixture.caseBundle.caseMsDoa.estResult.doaParamEst(:);

refOpt = fixture.staticBaseOpt;
refOpt.initDoaParam = msDoaAnchor;
msZeroOpt = fixture.staticBaseOpt;
msZeroOpt.initDoaParam = msDoaAnchor;
msZeroOpt.satWeight = [1; 0];

[modelRef, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.viewRefOnly.sceneRef, fixture.viewRefOnly.rxSigSf, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewRefOnly.doaGrid, ...
  fixture.fdRange, 1, refOpt);
[modelMsZero, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.viewMs.sceneRef, fixture.viewMs.rxSigSf, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, 1, msZeroOpt);

initRef = buildDoaDopplerSfInit(modelRef, []);
initMsZero = buildDoaDopplerSfInit(modelMsZero, []);

fdRefInitRefHz = initRef(end);
fdRefInitMsZeroHz = initMsZero(end);
fdRefInitDiffHz = fdRefInitMsZeroHz - fdRefInitRefHz;

tolHz = 1e-9 * max(1, abs(fdRefInitRefHz));
if abs(fdRefInitDiffHz) > tolHz
  error('regressionSfStaticFdInitZeroWeightInvariant:FdInitMismatch', ...
    ['With identical DoA anchors, the zero-weight multi-sat static fdRef ', ...
     'initializer must reduce to the ref-only initializer.']);
end

fprintf('  fixed DoA anchor (deg)         : [%0.6f, %0.6f]\n', msDoaAnchor(1), msDoaAnchor(2));
fprintf('  ref-only fdRef init (Hz)       : %.12f\n', fdRefInitRefHz);
fprintf('  zero-weight fdRef init (Hz)    : %.12f\n', fdRefInitMsZeroHz);
fprintf('  zero-weight init diff (Hz)     : %.6g\n', fdRefInitDiffHz);
fprintf('PASS: regressionSfStaticFdInitZeroWeightInvariant\n');


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
