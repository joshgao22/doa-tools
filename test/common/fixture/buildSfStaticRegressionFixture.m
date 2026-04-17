function fixture = buildSfStaticRegressionFixture()
%BUILDSFSTATICREGRESSIONFIXTURE Build one reusable SF static regression case.
% This helper centralizes the single-frame two-satellite fixture shared by
% the SF static regression scripts. It keeps the dev-script setup out of the
% regression files so the contracts stay focused on one behavior each.

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
[refState, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = scene.satIdx(otherSatIdxLocal);

truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
truth = struct();
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
truth.fdRateTrueHzPerSec = 0;
truth.refSatIdxLocal = refSatIdxLocal;
truth.refSatIdxGlobal = scene.satIdx(refSatIdxLocal);
truth.otherSatIdxLocal = otherSatIdxLocal;
truth.otherSatIdxGlobal = otherSatIdxGlobal;
truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
truth.deltaFdTrueHz = reshape(truthDopplerState.deltaFdRefFrame, [], 1);
truth.refStateSource = string(refState.source);

steeringInfo = getSceneSteering(scene, wavelen);
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
  steeringInfo, pilotWave, getLinkParam(scene, wavelen), carrierFreq, ...
  waveInfo.sampleRate, 1e-6, ones(scene.numSat, scene.numUser), snapOpt);

searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];
viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
rxSigRefOnly = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigRefOnly, gridSize, searchRange, E);

sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, gridSize, searchRange, E);

fdRange = expandRangeToTruth([-1e5, 0], ...
  [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;
weightSweepAlpha = [0; 0.25; 0.5; 1];
staticMsHalfWidth = [0.002; 0.002];

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, fdRange, truth, otherSatIdxGlobal, false, ...
  doaOnlyOpt, staticBaseOpt, weightSweepAlpha, staticMsHalfWidth);

staticZeroWeightRefOpt = staticBaseOpt;
staticZeroWeightRefOpt.satWeight = [1; 0];
staticZeroWeightRefOpt.initDoaParam = caseBundle.caseRefDoa.estResult.doaParamEst(:);
caseStaticZeroWeightRefInit = localRunStaticCase( ...
  "MS-SF-Static-W0.00-RefInit", "multi", viewMs, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, fdRange, false, staticZeroWeightRefOpt);

fixture = struct();
fixture.scene = scene;
fixture.sceneRefOnly = sceneRefOnly;
fixture.sceneOtherOnly = sceneOtherOnly;
fixture.rxSig = rxSig;
fixture.rxSigRefOnly = rxSigRefOnly;
fixture.rxSigOtherOnly = rxSigOtherOnly;
fixture.viewMs = viewMs;
fixture.viewRefOnly = viewRefOnly;
fixture.viewOtherOnly = viewOtherOnly;
fixture.truth = truth;
fixture.truthDopplerState = truthDopplerState;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.E = E;
fixture.wavelen = wavelen;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.caseBundle = caseBundle;
fixture.caseStaticZeroWeightRefInit = caseStaticZeroWeightRefInit;
fixture.staticBaseOpt = staticBaseOpt;
fixture.doaOnlyOpt = doaOnlyOpt;
fixture.weightSweepAlpha = weightSweepAlpha;
end


function caseInfo = localRunStaticCase(displayName, satMode, view, pilotWave, carrierFreq, ...
  sampleRate, fdRange, verbose, modelOpt)
%LOCALRUNSTATICCASE Run one SF static estimator case.

[estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  view.sceneRef, view.rxSigSf, pilotWave, carrierFreq, sampleRate, ...
  view.doaGrid, fdRange, 1, [], verbose, modelOpt);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa-doppler", "static", estResult);
end
