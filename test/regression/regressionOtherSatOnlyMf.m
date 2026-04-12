% Regression check for the other-satellite-only MF dynamic path.
% This script focuses on one narrow contract only: after selecting the
% non-reference satellite as a single-satellite subset, the MF generator,
% model builder, initializer, and known-rate solver should still remain
% self-consistent and recover that subset truth.
clear(); close all;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionOtherSatOnlyMf ...\n');

truthHyp = struct();
truthHyp.localDoaArr = localExtractTruthLocalDoa(fixture.sceneSeqOtherOnly);
truthHyp.fdSat = fixture.truthOther.fdSatFitHz(:);
truthHyp.fdRateSat = fixture.truthOther.fdRateSatTrueHzPerSec(:);

modelTruthOpt = struct();
modelTruthOpt.useLogObjective = false;
modelTruthOpt.phaseMode = 'continuous';
modelTruthOpt.fdRateMode = 'known';
modelTruthOpt.fdRateKnown = fixture.truthOther.fdRateFit;
modelTruthOpt.steeringMode = 'framewise';
modelTruthOpt.steeringRefFrameIdx = fixture.sceneSeqOtherOnly.refFrameIdx;
[modelTruth, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeqOtherOnly, fixture.rxSigMfOtherOnly, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewOtherOnly.doaGrid, ...
  fixture.fdRangeOther, fixture.fdRateRangeOther, modelTruthOpt);
[objTruth, profTruth] = evalDoaDopplerDynProfileLike(modelTruth, truthHyp);
truthNorm2 = max(sum(profTruth.blockNorm2(:)), eps);
residualRatioTruth = profTruth.residualNorm / truthNorm2;

solveOpt = modelTruthOpt;
solveOpt.initDoaParam = fixture.truthOther.latlonTrueDeg(:) + [0.05; -0.04];
solveOpt.initDoaHalfWidth = [0.10; 0.10];
solveOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-10);
[modelSolve, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeqOtherOnly, fixture.rxSigMfOtherOnly, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewOtherOnly.doaGrid, ...
  fixture.fdRangeOther, fixture.fdRateRangeOther, solveOpt);
[initParamSolve, ~] = buildDoaDopplerMfInit(modelSolve, []);
[solveUse, optimInfo, initEvalDiag] = solveDoaDopplerMfBranches(modelSolve, initParamSolve, false);

latlonEst = reshape(solveUse.finalEvalDiag.latlon, 2, []);
angleErrDeg = calcLatlonAngleError(latlonEst, fixture.truthOther.latlonTrueDeg(:));
fdRefErrHz = abs(solveUse.finalEvalDiag.fdRef - fixture.truthOther.fdRefFit);
fdRateErrHzPerSec = abs(solveUse.finalEvalDiag.fdRate - fixture.truthOther.fdRateFit);

objTol = max(1e-8, 1e-8 * max(abs([solveUse.fval, initEvalDiag.obj])));
if residualRatioTruth > 1e-8
  error('regressionOtherSatOnlyMf:LargeTruthResidual', ...
    'The other-sat-only truth residual ratio %.3e exceeds tolerance.', residualRatioTruth);
end
if solveUse.fval > initEvalDiag.obj + objTol
  error('regressionOtherSatOnlyMf:OptimizationDidNotImprove', ...
    'The other-sat-only known-rate solve did not improve the initialized objective.');
end
if angleErrDeg > 0.01 || fdRefErrHz > 1 || fdRateErrHzPerSec > 50
  error('regressionOtherSatOnlyMf:LargeRecoveryError', ...
    ['The other-sat-only MF path did not recover the subset truth ', ...
     '(angle %.6f deg, fdRef %.6f Hz, fdRate %.6f Hz/s).'], ...
    angleErrDeg, fdRefErrHz, fdRateErrHzPerSec);
end

fprintf('  otherSat global idx     : %d\n', fixture.otherSatIdxGlobal);
fprintf('  truth residual ratio    : %.3e\n', residualRatioTruth);
fprintf('  solveVariant            : %s\n', string(optimInfo.solveVariant));
fprintf('  angle err (deg)         : %.6f\n', angleErrDeg);
fprintf('  fdRef err (Hz)          : %.6f\n', fdRefErrHz);
fprintf('  fdRate err (Hz/s)       : %.6f\n', fdRateErrHzPerSec);
fprintf('PASS: regressionOtherSatOnlyMf\n');


function localDoa = localExtractTruthLocalDoa(sceneSeq)
%LOCALEXTRACTTRUTHLOCALDOA Extract 2 x Ns x Nf local DoA truth.

localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
  return;
end
if ndims(localDoaRaw) == 4
  localDoa = localDoaRaw(:, :, 1, :);
  return;
end

error('regressionOtherSatOnlyMf:InvalidLocalDoaSize', ...
  'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one other-satellite-only MF case.

rng(253);

numUsr = 1;
numFrame = 7;
refFrameIdx = ceil(numFrame / 2);
frameIntvlSec = 1 / 750;

sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
numSym = 32;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);
usrLla = [37.78; 36.59; 0];

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = sceneRef.satIdx(otherSatIdxLocal);

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.spatial.model = 'dynamic';
snapOpt.spatial.refFrameIdx = sceneSeq.refFrameIdx;
snapOpt.phase.timeModel = 'global';
snapOpt.phase.frameModel = 'shared';
snapOpt.phase.sharedPhase = 2 * pi * rand(sceneSeq.numSat, sceneSeq.numUser);
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
snapOpt.precomp.linkParamCell = linkParamCell;

[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, 0, repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame), snapOpt);

sceneSeqOtherOnly = selectSatSceneSeq(sceneSeq, otherSatIdxLocal);
rxSigMfOtherOnly = selectRxSigBySat(rxSigCell, otherSatIdxLocal, 'multiFrame');
sceneOtherRef = sceneSeqOtherOnly.sceneCell{sceneSeqOtherOnly.refFrameIdx};

linkParamCellOther = cell(1, sceneSeqOtherOnly.numFrame);
for iFrame = 1:sceneSeqOtherOnly.numFrame
  linkParamCellOther{iFrame} = getLinkParam(sceneSeqOtherOnly.sceneCell{iFrame}, wavelen);
end
truthOther = buildDynTruthFromLinkParam(linkParamCellOther, sceneSeqOtherOnly.timeOffsetSec, sampleRate, 1);
truthOther.latlonTrueDeg = usrLla(1:2, 1);

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewOtherOnly = buildDoaDopplerEstView(sceneOtherRef, rxSigMfOtherOnly, [41, 41], searchRange, E, ...
  struct('sceneSeq', sceneSeqOtherOnly));

fdRangeOther = expandRangeToTruth([-2e5, 2e5], [truthOther.fdRefFit; truthOther.fdSatFitHz(:)], 0.1, 2e4);
fdRateRangeOther = expandRangeToTruth([-1e4, 1e4], ...
  [truthOther.fdRateFit; truthOther.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

fixture = struct();
fixture.sceneSeqOtherOnly = sceneSeqOtherOnly;
fixture.rxSigMfOtherOnly = rxSigMfOtherOnly;
fixture.viewOtherOnly = viewOtherOnly;
fixture.truthOther = truthOther;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRangeOther = fdRangeOther;
fixture.fdRateRangeOther = fdRateRangeOther;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
