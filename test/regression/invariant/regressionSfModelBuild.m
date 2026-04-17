% Regression check for SF model construction.
% This script focuses on one narrow contract only: the SF model builder must
% preserve the resolved reference-satellite identity, pass satWeight into
% the model, and enforce the requested local DoA box through its bounds.
clear(); close all;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionSfModelBuild ...\n');

satWeight = [1; 0.35];
initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.10; -0.07];
initDoaHalfWidth = [0.12; 0.09];

modelOpt = struct();
modelOpt.useLogObjective = true;
modelOpt.initFdCount = 65;
modelOpt.initDoaParam = initDoaParam;
modelOpt.initDoaHalfWidth = initDoaHalfWidth;
modelOpt.satWeight = satWeight;

[modelMulti, lbMulti, ubMulti] = buildDoaDopplerSfModel( ...
  fixture.scene, fixture.rxSig, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, 1, modelOpt);

sceneRefOnly = selectSatScene(fixture.scene, fixture.refSatIdxLocal);
rxSigRefOnly = selectRxSigBySat(fixture.rxSig, fixture.refSatIdxLocal, 'singleFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigRefOnly, fixture.gridSize, ...
  fixture.searchRange, fixture.E);
modelOptRefOnly = modelOpt;
modelOptRefOnly.satWeight = 1;
[modelRefOnly, lbRefOnly, ubRefOnly] = buildDoaDopplerSfModel( ...
  sceneRefOnly, rxSigRefOnly, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, viewRefOnly.doaGrid, fixture.fdRange, 1, modelOptRefOnly);

[~, refSatIdxLocalFull] = resolveReferenceSatState(fixture.scene, ...
  fixture.scene.satPosEci, fixture.scene.satVelEci);
if modelMulti.refSatIdxLocal ~= refSatIdxLocalFull
  error('regressionSfModelBuild:ReferenceMismatch', ...
    'The SF model builder changed the resolved reference-satellite index.');
end
if modelRefOnly.refSatIdxLocal ~= 1
  error('regressionSfModelBuild:SubsetReferenceMismatch', ...
    'The ref-only SF subset must keep refSatIdxLocal = 1.');
end
if any(abs(modelMulti.satWeight(:) - satWeight(:)) > 1e-12)
  error('regressionSfModelBuild:SatWeightMismatch', ...
    'model.satWeight does not match the configured satellite weights.');
end

baseRange = fixture.viewMs.doaGrid{1}.range;
expDoaLb = max(baseRange(:, 1), initDoaParam - initDoaHalfWidth);
expDoaUb = min(baseRange(:, 2), initDoaParam + initDoaHalfWidth);
localCheckBounds(lbMulti(1:2), ubMulti(1:2), expDoaLb, expDoaUb, 'multi');
localCheckBounds(lbRefOnly(1:2), ubRefOnly(1:2), expDoaLb, expDoaUb, 'refOnly');

if numel(lbMulti) ~= 3 || numel(ubMulti) ~= 3 || numel(lbRefOnly) ~= 3 || numel(ubRefOnly) ~= 3
  error('regressionSfModelBuild:BoundSizeMismatch', ...
    'Single-frame DoA-Doppler bounds must contain exactly 3 entries.');
end
if abs(lbMulti(3) - fixture.fdRange(1)) > 1e-12 || abs(ubMulti(3) - fixture.fdRange(2)) > 1e-12
  error('regressionSfModelBuild:FdRangeMismatch', ...
    'The SF fdRef bounds must match the configured fdRange.');
end

fprintf('  refSatIdxLocal (full) : %d\n', modelMulti.refSatIdxLocal);
fprintf('  satWeight             : [%.3f, %.3f]\n', modelMulti.satWeight(1), modelMulti.satWeight(2));
fprintf('  DoA bounds            : [%.6f, %.6f] x [%.6f, %.6f]\n', ...
  lbMulti(1), ubMulti(1), lbMulti(2), ubMulti(2));
fprintf('PASS: regressionSfModelBuild\n');


function localCheckBounds(doaLb, doaUb, expDoaLb, expDoaUb, caseTag)
%LOCALCHECKBOUNDS Check the enforced local DoA box.

tol = 1e-12;
if any(abs(doaLb(:) - expDoaLb(:)) > tol) || any(abs(doaUb(:) - expDoaUb(:)) > tol)
  error('regressionSfModelBuild:DoaBoundMismatch', ...
    '[%s] local DoA bounds do not match the requested anchored box.', caseTag);
end
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact single-frame static case.

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
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
[~, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);

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
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, 1e-6, ...
  ones(scene.numSat, scene.numUser), snapOpt);

searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];
viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);
fdRange = expandRangeToTruth([-1e5, 0], [truthDopplerState.fdRefRefFrame; truthDopplerState.fdSatRefFrame(:)], 0.1, 2e4);

fixture = struct();
fixture.scene = scene;
fixture.rxSig = rxSig;
fixture.viewMs = viewMs;
fixture.truth = struct('latlonTrueDeg', usrLla(1:2, 1));
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.E = E;
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
