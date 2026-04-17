function fixture = buildSfStaticSingleSnrRepeatFixture(repeatIdx, snrDb)
%BUILDSFSTATICSINGLESNRREPEATFIXTURE Build one replayable SF static repeat.
% This helper mirrors the fixed-SNR setup used by doaDopplerStatDualSatUraEci
% and materializes one specific Monte Carlo repeat so regression scripts can
% probe the exact SS/MS static solve path on a known failing seed.

arguments
  repeatIdx (1, 1) double
  snrDb (1, 1) double = 10
end

if repeatIdx < 1 || repeatIdx ~= floor(repeatIdx)
  error('buildSfStaticSingleSnrRepeatFixture:InvalidRepeatIdx', ...
    'repeatIdx must be a positive integer scalar.');
end

numUsr = 1;
numSym = 512;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
baseSeed = 253;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));
E = referenceEllipsoid('sphere');

weightSweepAlpha = [0; 0.25; 0.5; 1];
staticMsHalfWidth = [0.002; 0.002];
fdRange = [-2e5, 2e5];

taskSeed = baseSeed + repeatIdx - 1;

usrLla = [37.78; 36.59; 0];
gridSize = [50, 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

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
  error('buildSfStaticSingleSnrRepeatFixture:InvalidNumSat', ...
    'The replay fixture expects exactly two satellites.');
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

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);

rng(taskSeed, 'twister');
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGain, snapOpt);

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
  weightSweepAlpha, staticMsHalfWidth);

mainCase = [caseBundle.caseRefDoa, caseBundle.caseMsDoa, ...
  caseBundle.caseStaticRefOnly, caseBundle.caseStaticMs, caseBundle.weightCase];
caseTable = buildDoaDopplerSummaryTable(mainCase, truth, struct('mode', 'static'));
ablationTruth = [ ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneRefOnly), ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneOtherOnly)];
ablationTable = buildDoaDopplerCaseSummaryTable( ...
  [caseBundle.caseStaticRefAbl, caseBundle.caseStaticOtherOnly], ablationTruth);

fixture = struct();
fixture.repeatIdx = repeatIdx;
fixture.taskSeed = taskSeed;
fixture.baseSeed = baseSeed;
fixture.snrDb = snrDb;
fixture.pwrNoise = pwrNoise;
fixture.scene = scene;
fixture.sceneRefOnly = sceneRefOnly;
fixture.sceneOtherOnly = sceneOtherOnly;
fixture.truth = truth;
fixture.truthDopplerState = truthDopplerState;
fixture.linkParam = linkParam;
fixture.steeringInfo = steeringInfo;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.pathGain = pathGain;
fixture.snapOpt = snapOpt;
fixture.fdRange = fdRange;
fixture.wavelen = wavelen;
fixture.searchRange = searchRange;
fixture.gridSize = gridSize;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.refSatIdxGlobal = refSatIdxGlobal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
fixture.viewRefOnly = viewRefOnly;
fixture.viewOtherOnly = viewOtherOnly;
fixture.viewMs = viewMs;
fixture.doaOnlyOpt = doaOnlyOpt;
fixture.staticBaseOpt = staticBaseOpt;
fixture.weightSweepAlpha = weightSweepAlpha;
fixture.staticMsHalfWidth = staticMsHalfWidth;
fixture.caseBundle = caseBundle;
fixture.caseTable = caseTable;
fixture.ablationTable = ablationTable;
end
