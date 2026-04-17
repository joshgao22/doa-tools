% Strategy comparison for the CP-U DoA-release effect.
% This script compares two practical release modes around the same CP-K seed:
%   1) with an almost zero DoA search box, CP-U should only release fdRate
%      while leaving the CP-K DoA seed essentially unchanged;
%   2) with a wider DoA search box around the same CP-K seed, CP-U should be
%      allowed to move in the DoA variables as well;
%   3) the printed trace should expose the actual DoA bound widths and the
%      margins to those bounds, so the next debugging step can tell whether
%      DoA is pinned by constraints or by the objective itself.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running compareMfUnknownDoaReleaseEffect ...\n');

commonOpt = struct();
commonOpt.useLogObjective = false;
commonOpt.phaseMode = 'continuous';
commonOpt.steeringMode = 'framewise';
commonOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
commonOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.05; -0.04];
commonOpt.initDoaHalfWidth = [0.10; 0.10];
commonOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-10);

fdRateBiasHzPerSec = 150;
fdRateSeed = min(max(fixture.truth.fdRateFit + fdRateBiasHzPerSec, ...
  fixture.fdRateRange(1)), fixture.fdRateRange(2));

optKnown = commonOpt;
optKnown.fdRateMode = 'known';
optKnown.fdRateKnown = fdRateSeed;
[modelKnown, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optKnown);
[initParamKnown, ~] = buildDoaDopplerMfInit(modelKnown, []);
[solveKnown, optimInfoKnown] = solveDoaDopplerMfBranches(modelKnown, initParamKnown, false);

seedLatlon = reshape(solveKnown.finalEvalDiag.latlon, 2, 1);
initParamRelease = [solveKnown.optVar(:); fdRateSeed];
seedAngleErrDeg = calcLatlonAngleError(seedLatlon, fixture.truth.latlonTrueDeg(:));

optUnknownFdOnly = commonOpt;
optUnknownFdOnly.fdRateMode = 'unknown';
optUnknownFdOnly.initDoaParam = seedLatlon;
optUnknownFdOnly.initDoaHalfWidth = [1e-8; 1e-8];
[modelUnknownFdOnly, lbFdOnly, ubFdOnly] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optUnknownFdOnly);
[solveUnknownFdOnly, optimInfoUnknownFdOnly] = solveDoaDopplerMfBranches( ...
  modelUnknownFdOnly, initParamRelease, false);

optUnknownDoa = commonOpt;
optUnknownDoa.fdRateMode = 'unknown';
optUnknownDoa.initDoaParam = seedLatlon;
optUnknownDoa.initDoaHalfWidth = [0.03; 0.03];
[modelUnknownDoa, lbDoa, ubDoa] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optUnknownDoa);
[solveUnknownDoa, optimInfoUnknownDoa] = solveDoaDopplerMfBranches( ...
  modelUnknownDoa, initParamRelease, false);

latlonFdOnly = reshape(solveUnknownFdOnly.finalEvalDiag.latlon, 2, 1);
latlonDoa = reshape(solveUnknownDoa.finalEvalDiag.latlon, 2, 1);

doaMoveFdOnlyDeg = calcLatlonAngleError(latlonFdOnly, seedLatlon);
doaMoveDoaDeg = calcLatlonAngleError(latlonDoa, seedLatlon);
angleErrFdOnlyDeg = calcLatlonAngleError(latlonFdOnly, fixture.truth.latlonTrueDeg(:));
angleErrDoaDeg = calcLatlonAngleError(latlonDoa, fixture.truth.latlonTrueDeg(:));
fdRateGapFdOnly = abs(solveUnknownFdOnly.finalEvalDiag.fdRate - fixture.truth.fdRateFit);
fdRateGapDoa = abs(solveUnknownDoa.finalEvalDiag.fdRate - fixture.truth.fdRateFit);
moveNormFdOnly = norm(solveUnknownFdOnly.optVar(:) - initParamRelease(:));
moveNormDoa = norm(solveUnknownDoa.optVar(:) - initParamRelease(:));
objTol = max(1e-8, 1e-8 * max(abs([solveUnknownFdOnly.fval, solveUnknownDoa.fval])));
angleImproveDeg = angleErrFdOnlyDeg - angleErrDoaDeg;

doaTraceFdOnly = localBuildDoaConstraintTrace(lbFdOnly, ubFdOnly, initParamRelease, solveUnknownFdOnly.optVar(:));
doaTraceDoa = localBuildDoaConstraintTrace(lbDoa, ubDoa, initParamRelease, solveUnknownDoa.optVar(:));

if doaMoveFdOnlyDeg > 1e-5
  error('compareMfUnknownDoaReleaseEffect:FdOnlyDoaMoved', ...
    ['The fd-only CP-U branch moved DoA by %.6e deg, which is too large ', ...
     'for the near-zero DoA search box.'], doaMoveFdOnlyDeg);
end
if ~isfield(optimInfoUnknownFdOnly, 'candidateVariant') || isempty(optimInfoUnknownFdOnly.candidateVariant) || ...
    ~isfield(optimInfoUnknownDoa, 'candidateVariant') || isempty(optimInfoUnknownDoa.candidateVariant)
  error('compareMfUnknownDoaReleaseEffect:MissingUnknownCandidates', ...
    'Both CP-U branches must expose continuation candidates.');
end
if doaMoveDoaDeg <= doaMoveFdOnlyDeg + 1e-6 && angleImproveDeg <= 1e-6 && ...
    solveUnknownDoa.fval >= solveUnknownFdOnly.fval - objTol
  error('compareMfUnknownDoaReleaseEffect:NoDoaReleaseEffect', ...
    ['The wider DoA-release branch produced no additional DoA motion and ', ...
     'no angle improvement relative to the fd-only release branch.']);
end

fprintf('  CP-K seed solveVariant     : %s\n', string(optimInfoKnown.solveVariant));
fprintf('  CP-U fd-only solveVariant  : %s\n', string(optimInfoUnknownFdOnly.solveVariant));
fprintf('  CP-U DoA solveVariant      : %s\n', string(optimInfoUnknownDoa.solveVariant));
fprintf('  seed angle err (deg)       : %.6f\n', seedAngleErrDeg);
fprintf('  fd-only DoA move (deg)     : %.6e\n', doaMoveFdOnlyDeg);
fprintf('  DoA-release move (deg)     : %.6e\n', doaMoveDoaDeg);
fprintf('  fd-only angle err (deg)    : %.6f\n', angleErrFdOnlyDeg);
fprintf('  DoA-release angle err (deg): %.6f\n', angleErrDoaDeg);
fprintf('  angle improvement (deg)    : %.6f\n', angleImproveDeg);
fprintf('  fd-only fdRate gap (Hz/s)  : %.6f\n', fdRateGapFdOnly);
fprintf('  DoA-release fdRate gap     : %.6f\n', fdRateGapDoa);
fprintf('  fd-only fval               : %.6e\n', solveUnknownFdOnly.fval);
fprintf('  DoA-release fval           : %.6e\n', solveUnknownDoa.fval);
fprintf('  fd-only move norm          : %.6f\n', moveNormFdOnly);
fprintf('  DoA-release move norm      : %.6f\n', moveNormDoa);
fprintf('  fd-only DoA lb             : %s\n', localFormatNumericRow(doaTraceFdOnly.lb));
fprintf('  fd-only DoA ub             : %s\n', localFormatNumericRow(doaTraceFdOnly.ub));
fprintf('  fd-only DoA width          : %s\n', localFormatNumericRow(doaTraceFdOnly.width));
fprintf('  fd-only init margin        : %s\n', localFormatNumericRow(doaTraceFdOnly.initMargin));
fprintf('  fd-only final margin       : %s\n', localFormatNumericRow(doaTraceFdOnly.finalMargin));
fprintf('  DoA-release lb             : %s\n', localFormatNumericRow(doaTraceDoa.lb));
fprintf('  DoA-release ub             : %s\n', localFormatNumericRow(doaTraceDoa.ub));
fprintf('  DoA-release width          : %s\n', localFormatNumericRow(doaTraceDoa.width));
fprintf('  DoA-release init margin    : %s\n', localFormatNumericRow(doaTraceDoa.initMargin));
fprintf('  DoA-release final margin   : %s\n', localFormatNumericRow(doaTraceDoa.finalMargin));
fprintf('PASS: compareMfUnknownDoaReleaseEffect\n');


function doaTrace = localBuildDoaConstraintTrace(lb, ub, initParam, optVar)
%LOCALBUILDDOACONSTRAINTTRACE Summarize the first two DoA-variable bounds.

idxDoa = 1:2;
doaTrace = struct();
doaTrace.lb = reshape(lb(idxDoa), 1, []);
doaTrace.ub = reshape(ub(idxDoa), 1, []);
doaTrace.width = doaTrace.ub - doaTrace.lb;
doaTrace.initVal = reshape(initParam(idxDoa), 1, []);
doaTrace.finalVal = reshape(optVar(idxDoa), 1, []);
doaTrace.initMargin = min(doaTrace.initVal - doaTrace.lb, doaTrace.ub - doaTrace.initVal);
doaTrace.finalMargin = min(doaTrace.finalVal - doaTrace.lb, doaTrace.ub - doaTrace.finalVal);
end


function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one numeric row for compact logging.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


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

fdRange = localExpandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = localExpandRangeToTruth([-1e4, 0], ...
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


function value = localExpandRangeToTruth(defaultRange, truthVals, fracPad, absPad)
%LOCALEXPANDRANGETOTRUTH Expand a nominal range so that truth is included.

truthVals = truthVals(isfinite(truthVals));
if isempty(truthVals)
  value = defaultRange;
  return;
end

minTruth = min(truthVals);
maxTruth = max(truthVals);
spanTruth = max(maxTruth - minTruth, eps);
pad = max(absPad, fracPad * spanTruth);
value = [min(defaultRange(1), minTruth - pad), ...
         max(defaultRange(2), maxTruth + pad)];
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
