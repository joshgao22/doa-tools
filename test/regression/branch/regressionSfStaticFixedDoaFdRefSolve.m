% Regression check for SF static fixed-DoA fdRef solves.
% This script focuses on one diagnostic contract only:
%   1) when DoA is fixed to a trusted anchor, the multi-sat static 1-D fdRef
%      solve must recover a physically reasonable reference Doppler;
%   2) if this fixed-DoA solve stays healthy while the free multi-sat static
%      solve is worse, the remaining issue is in the joint DoA-fdRef coupling
%      rather than in the static reference-Doppler chain itself.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionSfStaticFixedDoaFdRefSolve ...\n');

caseBundle = fixture.caseBundle;
freeSsCase = caseBundle.caseStaticRefOnly;
freeMsCase = caseBundle.caseStaticMs;
msDoaAnchor = caseBundle.caseMsDoa.estResult.doaParamEst(:);
truthDoaAnchor = fixture.truth.latlonTrueDeg(:);

solveSsTruth = localSolveFixedDoaFdRef(fixture.viewRefOnly, fixture, truthDoaAnchor, 1);
solveMsTruth = localSolveFixedDoaFdRef(fixture.viewMs, fixture, truthDoaAnchor, [1; 1]);
solveSsMsDoa = localSolveFixedDoaFdRef(fixture.viewRefOnly, fixture, msDoaAnchor, 1);
solveMsMsDoa = localSolveFixedDoaFdRef(fixture.viewMs, fixture, msDoaAnchor, [1; 1]);

if abs(solveMsTruth.fdRefErrHz) > 5
  error('regressionSfStaticFixedDoaFdRefSolve:TruthFixedMsFdRefBad', ...
    'With truth DoA fixed, the multi-sat static fdRef solve should stay close to truth.');
end
if abs(solveMsMsDoa.fdRefErrHz) > 25
  error('regressionSfStaticFixedDoaFdRefSolve:DoaFixedMsFdRefBad', ...
    'With MS-SF-DoA fixed, the multi-sat static fdRef solve should stay near truth.');
end
if abs(solveMsTruth.fdRefErrHz) > abs(solveSsTruth.fdRefErrHz) + 5
  error('regressionSfStaticFixedDoaFdRefSolve:TruthFixedMultiWorseThanSingle', ...
    'With truth DoA fixed, the multi-sat fdRef solve should not be materially worse than the ref-only solve.');
end
if ~all(isfinite([solveMsTruth.obj, solveMsMsDoa.obj]))
  error('regressionSfStaticFixedDoaFdRefSolve:NonFiniteObjective', ...
    'The fixed-DoA static fdRef solve must keep a finite objective.');
end

freeMsAngleErrDeg = calcLatlonAngleError(freeMsCase.estResult.doaParamEst(:), fixture.truth.latlonTrueDeg(:));
freeSsAngleErrDeg = calcLatlonAngleError(freeSsCase.estResult.doaParamEst(:), fixture.truth.latlonTrueDeg(:));

fprintf('  free SS angle err (deg)        : %.6g\n', freeSsAngleErrDeg);
fprintf('  free MS angle err (deg)        : %.6g\n', freeMsAngleErrDeg);
fprintf('  free SS fdRef err (Hz)         : %.6f\n', freeSsCase.estResult.fdRefEst - fixture.truth.fdRefTrueHz);
fprintf('  free MS fdRef err (Hz)         : %.6f\n', freeMsCase.estResult.fdRefEst - fixture.truth.fdRefTrueHz);
fprintf('  fixed truth DoA fdRef err SS   : %.6f\n', solveSsTruth.fdRefErrHz);
fprintf('  fixed truth DoA fdRef err MS   : %.6f\n', solveMsTruth.fdRefErrHz);
fprintf('  fixed MS-SF-DoA fdRef err SS   : %.6f\n', solveSsMsDoa.fdRefErrHz);
fprintf('  fixed MS-SF-DoA fdRef err MS   : %.6f\n', solveMsMsDoa.fdRefErrHz);
fprintf('PASS: regressionSfStaticFixedDoaFdRefSolve\n');


function solveInfo = localSolveFixedDoaFdRef(view, fixture, doaAnchor, satWeight)
%LOCALSOLVEFIXEDDOAFDREF Solve one 1-D static fdRef fit with DoA fixed.

modelOpt = fixture.staticBaseOpt;
modelOpt.useLogObjective = false;
modelOpt.satWeight = satWeight;
[model, ~, ~] = buildDoaDopplerSfModel( ...
  view.sceneRef, view.rxSigSf, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, view.doaGrid, fixture.fdRange, 1, modelOpt);

objFun = @(fdRef) localEvalFixedDoaObjective(model, doaAnchor, fdRef);
[fdRefBest, objBest] = fminbnd(objFun, model.fdRange(1), model.fdRange(2));
[~, ~, ~, auxBest] = evalDoaDopplerSfProfileLike(model, [doaAnchor(:); fdRefBest]);

solveInfo = struct();
solveInfo.fdRefEstHz = fdRefBest;
solveInfo.fdRefErrHz = fdRefBest - fixture.truth.fdRefTrueHz;
solveInfo.obj = objBest;
solveInfo.objectiveSat = auxBest.objectiveSat(:);
solveInfo.fdSatEstHz = auxBest.fd(:);
solveInfo.deltaFdRefEstHz = auxBest.deltaFd(:);
end


function obj = localEvalFixedDoaObjective(model, doaAnchor, fdRef)
%LOCALEVALFIXEDDOAOBJECTIVE Evaluate one fixed-DoA static objective slice.

obj = evalDoaDopplerSfProfileLike(model, [doaAnchor(:); fdRef]);
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
