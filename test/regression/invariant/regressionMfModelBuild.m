function regressionMfModelBuild(varargin)
% Regression check for MF model construction.
% This script verifies one narrow contract only: the MF model builder must
% preserve the reference-satellite identity and must enforce the requested
% local DoA box through its bounds for both CP-K and CP-U paths.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

  fprintf('Running regressionMfModelBuild ...\n');

initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.08; -0.06];
initDoaHalfWidth = [0.12; 0.10];

optKnown = struct();
optKnown.useLogObjective = true;
optKnown.phaseMode = 'continuous';
optKnown.fdRateMode = 'known';
optKnown.fdRateKnown = fixture.truth.fdRateFit;
optKnown.steeringMode = 'framewise';
optKnown.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
optKnown.initDoaParam = initDoaParam;
optKnown.initDoaHalfWidth = initDoaHalfWidth;

optUnknown = optKnown;
optUnknown.fdRateMode = 'unknown';
optUnknown = rmfield(optUnknown, 'fdRateKnown');

[modelKnown, lbKnown, ubKnown] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optKnown);
[modelUnknown, lbUnknown, ubUnknown] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, optUnknown);

[~, refSatIdxLocal] = resolveReferenceSatState(fixture.sceneRef, ...
  fixture.sceneRef.satPosEci, fixture.sceneRef.satVelEci);
if modelKnown.refSatIdxLocal ~= refSatIdxLocal || modelUnknown.refSatIdxLocal ~= refSatIdxLocal
  error('regressionMfModelBuild:ReferenceMismatch', ...
    'MF model builder changed the resolved reference-satellite index.');
end

if numel(modelKnown.pilotPadCell) ~= fixture.sceneSeq.numFrame || ...
    numel(modelKnown.timeSecCell) ~= fixture.sceneSeq.numFrame
  error('regressionMfModelBuild:FrameCellSizeMismatch', ...
    'pilotPadCell and timeSecCell must contain one entry per frame.');
end
if any(abs(modelKnown.timeOffsetSec(:) - fixture.sceneSeq.timeOffsetSec(:)) > 1e-12)
  error('regressionMfModelBuild:TimeOffsetMismatch', ...
    'MF model timeOffsetSec must match sceneSeq.timeOffsetSec exactly.');
end

baseRange = fixture.viewMs.doaGrid{1}.range;
expDoaLb = max(baseRange(:, 1), initDoaParam - initDoaHalfWidth);
expDoaUb = min(baseRange(:, 2), initDoaParam + initDoaHalfWidth);

localCheckBounds(lbKnown(1:2), ubKnown(1:2), expDoaLb, expDoaUb, 'known');
localCheckBounds(lbUnknown(1:2), ubUnknown(1:2), expDoaLb, expDoaUb, 'unknown');

if numel(lbKnown) ~= 3 || numel(ubKnown) ~= 3
  error('regressionMfModelBuild:KnownBoundSize', ...
    'Known-rate MF bounds must contain exactly 3 entries.');
end
if numel(lbUnknown) ~= 4 || numel(ubUnknown) ~= 4
  error('regressionMfModelBuild:UnknownBoundSize', ...
    'Unknown-rate MF bounds must contain exactly 4 entries.');
end
if abs(lbKnown(3) - fixture.fdRange(1)) > 1e-12 || abs(ubKnown(3) - fixture.fdRange(2)) > 1e-12
  error('regressionMfModelBuild:KnownFdRangeMismatch', ...
    'Known-rate fdRef bounds must match the configured fdRange.');
end
if abs(lbUnknown(3) - fixture.fdRange(1)) > 1e-12 || abs(ubUnknown(3) - fixture.fdRange(2)) > 1e-12
  error('regressionMfModelBuild:UnknownFdRangeMismatch', ...
    'Unknown-rate fdRef bounds must match the configured fdRange.');
end
if abs(lbUnknown(4) - fixture.fdRateRange(1)) > 1e-12 || ...
    abs(ubUnknown(4) - fixture.fdRateRange(2)) > 1e-12
  error('regressionMfModelBuild:FdRateRangeMismatch', ...
    'Unknown-rate fdRate bounds must match the configured fdRateRange.');
end

  fprintf('  refSatIdxLocal  : %d\n', modelKnown.refSatIdxLocal);
  fprintf('  known DoA bounds: [%.6f, %.6f] x [%.6f, %.6f]\n', ...
    lbKnown(1), ubKnown(1), lbKnown(2), ubKnown(2));
  fprintf('  unknown vars    : %d\n', numel(lbUnknown));
  fprintf('PASS: regressionMfModelBuild\n');



end

function localCheckBounds(doaLb, doaUb, expDoaLb, expDoaUb, caseTag)
%LOCALCHECKBOUNDS Check the enforced local DoA box.

tol = 1e-12;
if any(abs(doaLb(:) - expDoaLb(:)) > tol) || any(abs(doaUb(:) - expDoaUb(:)) > tol)
  error('regressionMfModelBuild:DoaBoundMismatch', ...
    '[%s] local DoA bounds do not match the requested anchored box.', caseTag);
end
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact MF model-build case.

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
noiseVar = 1e-6;
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, noiseVar, pathGainCell, snapOpt);

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [41, 41], searchRange, E, ...
  struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = expandRangeToTruth([-1e4, 0], ...
  [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.sceneRef = sceneRef;
fixture.truth = truth;
fixture.viewMs = viewMs;
fixture.rxSigCell = rxSigCell;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
end
