function regressionMfUnknownFixedDoaWarmAnchor(varargin)
% Regression check for the fixed-DoA CP-U warm-anchor branch.
% This regression checks one narrow contract only:
%   1) when the caller marks one CP-U solve as freezeDoa, the solver must
%      represent that intent with an exact DoA constraint rather than a
%      nearly collapsed DoA box;
%   2) the unknown-rate branch must still release fdRate away from the
%      biased seed;
%   3) the final DoA must stay exactly on the trusted warm anchor;
%   4) the frozen-DoA solve must keep a valid resolved fit. This regression
%      guards solver-path semantics rather than truth-monotonic improvement.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

fprintf('Running regressionMfUnknownFixedDoaWarmAnchor ...\n');
  fprintf('  verbose trace enabled for fixed-DoA warm-anchor regression.\n');

commonOpt = struct();
commonOpt.useLogObjective = false;
commonOpt.phaseMode = 'continuous';
commonOpt.steeringMode = 'framewise';
commonOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
commonOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.05; -0.04];
commonOpt.initDoaHalfWidth = [0.10; 0.10];
commonOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-10);
commonOpt.verbose = verbose;

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
[solveKnown, ~] = solveDoaDopplerMfBranches(modelKnown, initParamKnown, verbose);

doAAnchor = solveKnown.optVar(1:2);
optUnknown = commonOpt;
optUnknown.fdRateMode = 'unknown';
optUnknown.initDoaParam = doAAnchor(:);
optUnknown.initDoaHalfWidth = [1e-8; 1e-8];
optUnknown.freezeDoa = true;
[modelUnknown, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optUnknown);
initParamRelease = [doAAnchor(:); solveKnown.optVar(3); fdRateSeed];
[solveUnknown, optimInfoUnknown] = solveDoaDopplerMfBranches(modelUnknown, initParamRelease, verbose);

objTol = max(1e-8, 1e-8 * max(abs([solveKnown.fval, solveUnknown.fval])));
doaMove = norm(solveUnknown.optVar(1:2) - doAAnchor(:));
fdRateMoveHzPerSec = abs(solveUnknown.finalEvalDiag.fdRate - fdRateSeed);
seedTruthGap = abs(fdRateSeed - fixture.truth.fdRateFit);
finalTruthGap = abs(solveUnknown.finalEvalDiag.fdRate - fixture.truth.fdRateFit);

if string(optimInfoUnknown.solveVariant) ~= "mainWarmAnchorFixedDoa"
  error('regressionMfUnknownFixedDoaWarmAnchor:WrongSolveVariant', ...
    'freezeDoa CP-U must report the fixed-DoA warm-anchor branch.');
end
if doaMove > 1e-12
  error('regressionMfUnknownFixedDoaWarmAnchor:DoaWasReleased', ...
    'freezeDoa CP-U changed DoA instead of staying on the trusted anchor.');
end
if fdRateMoveHzPerSec <= 1e-3
  error('regressionMfUnknownFixedDoaWarmAnchor:FdRateDidNotRelease', ...
    'freezeDoa CP-U did not release fdRate away from the biased seed.');
end
if solveUnknown.fval > solveKnown.fval + objTol
  error('regressionMfUnknownFixedDoaWarmAnchor:UnknownWorseThanSeed', ...
    'freezeDoa CP-U fit is worse than the biased CP-K seed.');
end
finalResidualNorm = localGetFieldOrDefault(localGetFieldOrDefault(solveUnknown, 'finalEvalDiag', struct()), 'residualNorm', NaN);
if ~(isfinite(finalResidualNorm) && finalResidualNorm > 0)
  error('regressionMfUnknownFixedDoaWarmAnchor:InvalidFinalResidual', ...
    'freezeDoa CP-U must return one finite final residual norm.');
end

fprintf('  CP-U solveVariant           : %s\n', string(optimInfoUnknown.solveVariant));
fprintf('  DoA anchor (deg)            : [%.6f, %.6f]\n', doAAnchor(1), doAAnchor(2));
fprintf('  final DoA move norm         : %.3e\n', doaMove);
fprintf('  biased fdRate seed (Hz/s)   : %.6f\n', fdRateSeed);
fprintf('  final fdRate (Hz/s)         : %.6f\n', solveUnknown.finalEvalDiag.fdRate);
fprintf('  fdRate move from seed       : %.6f\n', fdRateMoveHzPerSec);
fprintf('  seed gap to truth (Hz/s)    : %.6f\n', seedTruthGap);
fprintf('  final gap to truth (Hz/s)   : %.6f\n', finalTruthGap);
fprintf('PASS: regressionMfUnknownFixedDoaWarmAnchor\n');

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


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field with default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
