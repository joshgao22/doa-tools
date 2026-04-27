function regressionSfObjectiveShape(varargin)
% Regression check for the SF profile likelihood shape.
% This script focuses on one narrow contract only: in a noiseless static
% multi-satellite case, the SF objective must attain its best value at the
% truth hypothesis and must worsen under small DoA or fdRef perturbations.

opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

  fprintf('Running regressionSfObjectiveShape ...\n');

modelOpt = struct();
modelOpt.useLogObjective = false;
[model, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.scene, fixture.rxSig, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, 1, modelOpt);

truthOptVar = [fixture.truth.latlonTrueDeg(:); fixture.truth.fdRefTrueHz];
[objTruth, ~, ~, auxTruth] = evalDoaDopplerSfProfileLike(model, truthOptVar);
truthNorm2 = max(sum(auxTruth.countPerSat(:)) * max(auxTruth.noiseVarGlobal, eps), eps);
residualRatioTruth = auxTruth.residualNorm / max(sum(auxTruth.countPerSat(:)) * max(abs(fixture.pathGain).^2), eps);

latOffsetsDeg = [-0.05, 0, 0.05];
latObj = zeros(size(latOffsetsDeg));
for iVal = 1:numel(latOffsetsDeg)
  optVarCur = truthOptVar;
  optVarCur(1) = optVarCur(1) + latOffsetsDeg(iVal);
  latObj(iVal) = evalDoaDopplerSfProfileLike(model, optVarCur);
end

fdOffsetsHz = [-40, 0, 40];
fdObj = zeros(size(fdOffsetsHz));
for iVal = 1:numel(fdOffsetsHz)
  optVarCur = truthOptVar;
  optVarCur(3) = optVarCur(3) + fdOffsetsHz(iVal);
  fdObj(iVal) = evalDoaDopplerSfProfileLike(model, optVarCur);
end

[~, latBestIdx] = min(latObj);
[~, fdBestIdx] = min(fdObj);
if ~isfinite(objTruth)
  error('regressionSfObjectiveShape:InvalidTruthObjective', ...
    'Truth objective must be finite.');
end
if residualRatioTruth > 1e-8
  error('regressionSfObjectiveShape:LargeTruthResidual', ...
    'Truth residual ratio %.3e exceeds the regression tolerance.', residualRatioTruth);
end
if latBestIdx ~= 2
  error('regressionSfObjectiveShape:LatSliceMismatch', ...
    'The truth point is not the best point on the local latitude slice.');
end
if fdBestIdx ~= 2
  error('regressionSfObjectiveShape:FdSliceMismatch', ...
    'The truth point is not the best point on the local fdRef slice.');
end
if ~(latObj(2) < latObj(1) && latObj(2) < latObj(3))
  error('regressionSfObjectiveShape:LatPerturbationNotWorse', ...
    'Small latitude perturbations must worsen the SF objective.');
end
if ~(fdObj(2) < fdObj(1) && fdObj(2) < fdObj(3))
  error('regressionSfObjectiveShape:FdPerturbationNotWorse', ...
    'Small fdRef perturbations must worsen the SF objective.');
end

  fprintf('  truth objective       : %.6f\n', objTruth);
  fprintf('  truth residual ratio  : %.3e\n', residualRatioTruth);
  fprintf('  lat slice objective   : [%.6f, %.6f, %.6f]\n', latObj(1), latObj(2), latObj(3));
  fprintf('  fd slice objective    : [%.6f, %.6f, %.6f]\n', fdObj(1), fdObj(2), fdObj(3));
  fprintf('PASS: regressionSfObjectiveShape\n');


end

function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one noiseless static multi-sat case.

rng(253);

numUsr = 1;
numSym = 32;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

gridSize = [41, 41];
searchMarginDeg = 5;
usrLla = [37.78; 36.59; 0];

utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);

scene = genMultiSatScene(utc, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1);
linkParam = getLinkParam(scene, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, 0, ...
  pathGain, snapOpt);

searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];
viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);
fdRange = expandRangeToTruth([-1e5, 0], [truthDopplerState.fdRefRefFrame; truthDopplerState.fdSatRefFrame(:)], 0.1, 2e4);

fixture = struct();
fixture.scene = scene;
fixture.rxSig = rxSig;
fixture.viewMs = viewMs;
fixture.truth = struct();
fixture.truth.latlonTrueDeg = usrLla(1:2, 1);
fixture.truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.pathGain = pathGain;
end
