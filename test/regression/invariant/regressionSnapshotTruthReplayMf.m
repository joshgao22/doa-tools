function regressionSnapshotTruthReplayMf(varargin)
% Regression check for MF truth replay consistency.
% This script verifies that the dynamic generator and the MF evaluator stay
% self-consistent under the project continuous-phase model. It uses one
% noiseless dynamic scene and checks:
%   1) the truth hypothesis yields a near-zero residual;
%   2) small DoA or Doppler perturbations degrade the objective.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

  fprintf('Running regressionSnapshotTruthReplayMf ...\n');

modelOpt = struct();
modelOpt.useLogObjective = false;
modelOpt.phaseMode = 'continuous';
modelOpt.fdRateMode = 'unknown';
modelOpt.steeringMode = 'framewise';
modelOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
modelOpt.debugEnable = false;

[model, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, modelOpt);

truthHyp = localBuildTruthHypothesis(fixture.sceneSeq, fixture.truth);
[objTruth, profTruth] = evalDoaDopplerDynProfileLike(model, truthHyp);

hypDoaPert = truthHyp;
hypDoaPert.localDoaArr = truthHyp.localDoaArr + repmat(reshape(deg2rad([0.05; 0.03]), 2, 1, 1), ...
  1, size(truthHyp.localDoaArr, 2), size(truthHyp.localDoaArr, 3));
[objDoaPert, profDoaPert] = evalDoaDopplerDynProfileLike(model, hypDoaPert);

hypFdPert = truthHyp;
hypFdPert.fdSat = hypFdPert.fdSat + [25; -35];
hypFdPert.fdRateSat = hypFdPert.fdRateSat + [120; -160];
[objFdPert, profFdPert] = evalDoaDopplerDynProfileLike(model, hypFdPert);

truthNorm2 = max(sum(profTruth.blockNorm2(:)), eps);
residualRatioTruth = profTruth.residualNorm / truthNorm2;

if ~isfinite(objTruth)
  error('regressionSnapshotTruthReplayMf:InvalidTruthObjective', ...
    'Truth objective must be finite.');
end
if residualRatioTruth > 1e-8
  error('regressionSnapshotTruthReplayMf:LargeTruthResidual', ...
    'Truth residual ratio %.3e exceeds the regression tolerance.', residualRatioTruth);
end
if ~(objTruth < objDoaPert)
  error('regressionSnapshotTruthReplayMf:DoaPerturbationNotWorse', ...
    'Truth objective must be lower than the DoA-perturbed objective.');
end
if ~(objTruth < objFdPert)
  error('regressionSnapshotTruthReplayMf:DopplerPerturbationNotWorse', ...
    'Truth objective must be lower than the Doppler-perturbed objective.');
end

  fprintf('  truth residual ratio : %.3e\n', residualRatioTruth);
  fprintf('  truth objective      : %.6f\n', objTruth);
  fprintf('  DoA-perturbed obj    : %.6f\n', objDoaPert);
  fprintf('  fd-perturbed obj     : %.6f\n', objFdPert);
  fprintf('  DoA-perturbed resid  : %.6f\n', profDoaPert.residualNorm);
  fprintf('  fd-perturbed resid   : %.6f\n', profFdPert.residualNorm);
  fprintf('PASS: regressionSnapshotTruthReplayMf\n');



end

function hyp = localBuildTruthHypothesis(sceneSeq, truth)
%LOCALBUILDTRUTHHYPOTHESIS Build one evaluator hypothesis from stored truth.

hyp = struct();
hyp.localDoaArr = localExtractTruthLocalDoa(sceneSeq);
hyp.fdSat = truth.fdSatFitHz(:);
hyp.fdRateSat = truth.fdRateSatTrueHzPerSec(:);
end


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

error('regressionSnapshotTruthReplayMf:InvalidLocalDoaSize', ...
  'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one noiseless dynamic MF regression case.

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

pathGainCell = repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame);
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, 0, pathGainCell, snapOpt);

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
