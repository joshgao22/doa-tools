function fixture = buildSfStaticJointCouplingFixture()
%BUILDSFSTATICJOINTCOUPLINGFIXTURE Build the locked SF static coupling case.
% This fixture mirrors the main setup of doaDopplerStatDualSatUraEci so the
% branch regression can probe the exact static SS/MS coupling behavior
% currently seen in the dev script. The helper is intentionally narrow and
% keeps the regression file focused on contract checks and diagnostics.

rng(253);

numUsr = 1;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

usrLla = [37.78; 36.59; 0];
gridSize = [50, 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];
fdRange = [-2e5, 2e5];
weightSweepAlpha = [0; 0.25; 0.5; 1];
staticMsHalfWidth = [0.002; 0.002];

utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);

[~, satAccess] = findVisibleSatFromTle(utc, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccess, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

scene = genMultiSatScene(utc, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal);
linkParam = getLinkParam(scene, wavelen);
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
[refState, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);

if scene.numSat ~= 2
  error('buildSfStaticJointCouplingFixture:InvalidNumSat', ...
    'The locked SF static coupling fixture expects exactly two satellites.');
end
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

truth = struct();
truth.utc = scene.utc;
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.usrElevationDeg = reshape(scene.accessInfo.usrElevationDeg(:, 1), 1, []);
truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
truth.deltaFdTrueHz = reshape(truthDopplerState.deltaFdRefFrame, [], 1);
truth.localDoaRef = reshape(scene.localDoa(:, refSatIdxLocal), 2, 1);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;
truth.fdRateTrueHzPerSec = NaN;
fdRange = expandRangeToTruth(fdRange, [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, ...
  1e-6, ones(scene.numSat, scene.numUser), snapOpt);

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
rxSigRefOnly = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigRefOnly, gridSize, searchRange, E);

sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, gridSize, searchRange, E);

viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, false, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha, [0.01; 0.01]);

fixture = struct();
fixture.scene = scene;
fixture.viewRefOnly = viewRefOnly;
fixture.viewOtherOnly = viewOtherOnly;
fixture.viewMs = viewMs;
fixture.rxSig = rxSig;
fixture.truth = truth;
fixture.truthDopplerState = truthDopplerState;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.wavelen = wavelen;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.refSatIdxGlobal = refSatIdxGlobal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
fixture.staticBaseOpt = staticBaseOpt;
fixture.doaOnlyOpt = doaOnlyOpt;
fixture.weightSweepAlpha = weightSweepAlpha;
fixture.caseBundle = caseBundle;
end
