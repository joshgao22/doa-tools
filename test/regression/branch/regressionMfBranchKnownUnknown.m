% Regression check for MF known-rate / unknown-rate branch solving.
% This script focuses on one narrow contract only:
%   1) the unknown-rate branch must build continuation candidates;
%   2) both CP-K and CP-U must improve upon the same initialized state;
%   3) CP-U must not fit worse than CP-K on the same noiseless data.
clear(); close all;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionMfBranchKnownUnknown ...\n');

commonOpt = struct();
commonOpt.useLogObjective = false;
commonOpt.phaseMode = 'continuous';
commonOpt.steeringMode = 'framewise';
commonOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
commonOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.05; -0.04];
commonOpt.initDoaHalfWidth = [0.10; 0.10];
commonOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-10);

optKnown = commonOpt;
optKnown.fdRateMode = 'known';
optKnown.fdRateKnown = fixture.truth.fdRateFit;
[modelKnown, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optKnown);
[initParamKnown, ~] = buildDoaDopplerMfInit(modelKnown, []);
[solveKnown, optimInfoKnown, initEvalDiagKnown] = solveDoaDopplerMfBranches(modelKnown, initParamKnown, false);

optUnknown = commonOpt;
optUnknown.fdRateMode = 'unknown';
[modelUnknown, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optUnknown);
[initParamUnknown, ~] = buildDoaDopplerMfInit(modelUnknown, []);
[solveUnknown, optimInfoUnknown, initEvalDiagUnknown] = solveDoaDopplerMfBranches(modelUnknown, initParamUnknown, false);

latlonKnown = reshape(solveKnown.finalEvalDiag.latlon, 2, []);
latlonUnknown = reshape(solveUnknown.finalEvalDiag.latlon, 2, []);
angleErrKnownDeg = calcLatlonAngleError(latlonKnown, fixture.truth.latlonTrueDeg(:));
angleErrUnknownDeg = calcLatlonAngleError(latlonUnknown, fixture.truth.latlonTrueDeg(:));
fdRefErrKnownHz = abs(solveKnown.finalEvalDiag.fdRef - fixture.truth.fdRefFit);
fdRefErrUnknownHz = abs(solveUnknown.finalEvalDiag.fdRef - fixture.truth.fdRefFit);
fdRateErrUnknownHzPerSec = abs(solveUnknown.finalEvalDiag.fdRate - fixture.truth.fdRateFit);

objImproveKnown = initEvalDiagKnown.obj - solveKnown.fval;
objImproveUnknown = initEvalDiagUnknown.obj - solveUnknown.fval;
objTol = max(1e-8, 1e-8 * max(abs([solveKnown.fval, solveUnknown.fval, initEvalDiagKnown.obj, initEvalDiagUnknown.obj])));

if objImproveKnown < -objTol
  error('regressionMfBranchKnownUnknown:KnownDidNotImprove', ...
    'CP-K final objective is worse than its initialized objective.');
end
if objImproveUnknown < -objTol
  error('regressionMfBranchKnownUnknown:UnknownDidNotImprove', ...
    'CP-U final objective is worse than its initialized objective.');
end
if solveUnknown.fval > solveKnown.fval + objTol
  error('regressionMfBranchKnownUnknown:UnknownWorseThanKnown', ...
    'CP-U should not fit worse than CP-K on the same data.');
end
if ~isfield(optimInfoUnknown, 'candidateVariant') || isempty(optimInfoUnknown.candidateVariant)
  error('regressionMfBranchKnownUnknown:MissingUnknownCandidates', ...
    'CP-U must expose at least one continuation candidate in the multi-satellite CP path.');
end
candidateVariant = string(optimInfoUnknown.candidateVariant);
allowedTagMask = strcmp(candidateVariant, "main") | strcmp(candidateVariant, "mainWarmAnchor");
if ~all(allowedTagMask)
  error('regressionMfBranchKnownUnknown:InvalidUnknownCandidateTag', ...
    'CP-U continuation candidates expose unexpected solveVariant tags.');
end
if ~any(strcmp(candidateVariant, "mainWarmAnchor"))
  error('regressionMfBranchKnownUnknown:MissingContinuationAnchor', ...
    'CP-U continuation candidates must include mainWarmAnchor.');
end
angleImproveTolDeg = 1e-4;
if angleErrUnknownDeg > angleErrKnownDeg + angleImproveTolDeg
  error('regressionMfBranchKnownUnknown:UnknownAngleWorseThanKnown', ...
    'CP-U angle error should not be worse than CP-K after the branch simplification.');
end
if abs(angleErrKnownDeg - angleErrUnknownDeg) <= angleImproveTolDeg
  error('regressionMfBranchKnownUnknown:MissingUnknownAngleGain', ...
    'CP-U should provide a measurable angle improvement over CP-K in this noiseless regression case.');
end
if fdRefErrKnownHz > 1 || fdRefErrUnknownHz > 1 || fdRateErrUnknownHzPerSec > 50
  error('regressionMfBranchKnownUnknown:LargeFrequencyError', ...
    'CP-K/CP-U frequency errors are too large for the noiseless regression case.');
end

fprintf('  CP-K solveVariant      : %s\n', string(optimInfoKnown.solveVariant));
fprintf('  CP-U solveVariant      : %s\n', string(optimInfoUnknown.solveVariant));
fprintf('  CP-U candidate count   : %d\n', numel(candidateVariant));
fprintf('  CP-K angle err (deg)   : %.6f\n', angleErrKnownDeg);
fprintf('  CP-U angle err (deg)   : %.6f\n', angleErrUnknownDeg);
fprintf('  CP-K fdRef err (Hz)    : %.6f\n', fdRefErrKnownHz);
fprintf('  CP-U fdRef err (Hz)    : %.6f\n', fdRefErrUnknownHz);
fprintf('  CP-U fdRate err (Hz/s) : %.6f\n', fdRateErrUnknownHzPerSec);
fprintf('PASS: regressionMfBranchKnownUnknown\n');


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one noiseless dynamic multi-sat case.

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

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.latlonTrueDeg = usrLla(1:2, 1);

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

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [41, 41], searchRange, E, ...
  struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = expandRangeToTruth([-1e4, 0], ...
  [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.truth = truth;
fixture.viewMs = viewMs;
fixture.rxSigCell = rxSigCell;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
