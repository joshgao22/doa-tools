function regressionSfDoaOnlyPilotModelCrb(varargin)
% Regression check for SF DoA-only pilot-model CRB scaling.
% This case only protects the CRB mapping contract: the pilot-model
% DoA-only wrapper must reduce each per-satellite FIM contribution by the
% no-Doppler coherent projection gain squared, while the unit-gain mode
% remains unchanged.

opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
fixture = localBuildRegressionFixture();

  fprintf('Running regressionSfDoaOnlyPilotModelCrb ...\n');

pwrNoise = 1;
unitOpt = struct('doaType', 'latlon', 'effectiveGainMode', 'unit');
pilotOpt = unitOpt;
pilotOpt.effectiveGainMode = 'pilotModel';

sceneRefOnly = selectSatScene(fixture.scene, fixture.refSatIdxLocal);
[crbUnitRef, auxUnitRef] = crbPilotSfDoaOnlyEffective(sceneRefOnly, ...
  fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdSatTrueHz(fixture.refSatIdxLocal), ...
  1, pwrNoise, unitOpt);
[crbPilotRef, auxPilotRef] = crbPilotSfDoaOnlyEffective(sceneRefOnly, ...
  fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdSatTrueHz(fixture.refSatIdxLocal), ...
  1, pwrNoise, pilotOpt);
[crbUnitJoint, auxUnitJoint] = crbPilotSfDoaOnlyEffective(fixture.scene, ...
  fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdRefTrueHz, 1, pwrNoise, unitOpt);
[crbPilotJoint, auxPilotJoint] = crbPilotSfDoaOnlyEffective(fixture.scene, ...
  fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdRefTrueHz, 1, pwrNoise, pilotOpt);

localAssertFiniteCrb(crbUnitRef, 'unit ref');
localAssertFiniteCrb(crbPilotRef, 'pilot ref');
localAssertFiniteCrb(crbUnitJoint, 'unit joint');
localAssertFiniteCrb(crbPilotJoint, 'pilot joint');

localAssertUnitGain(auxUnitRef, 'ref');
localAssertUnitGain(auxUnitJoint, 'joint');
localAssertFimScaling(auxUnitRef, auxPilotRef, 'ref');
localAssertFimScaling(auxUnitJoint, auxPilotJoint, 'joint');

unitAngleRef = projectCrbToAngleMetric(crbUnitRef, fixture.truth.latlonTrueDeg, 'latlon');
pilotAngleRef = projectCrbToAngleMetric(crbPilotRef, fixture.truth.latlonTrueDeg, 'latlon');
unitAngleJoint = projectCrbToAngleMetric(crbUnitJoint, fixture.truth.latlonTrueDeg, 'latlon');
pilotAngleJoint = projectCrbToAngleMetric(crbPilotJoint, fixture.truth.latlonTrueDeg, 'latlon');
if pilotAngleRef + 1e-12 < unitAngleRef
  error('regressionSfDoaOnlyPilotModelCrb:PilotRefImprovedUnit', ...
    'Pilot-model ref CRB must not improve over the unit-gain CRB.');
end
if pilotAngleJoint + 1e-12 < unitAngleJoint
  error('regressionSfDoaOnlyPilotModelCrb:PilotJointImprovedUnit', ...
    'Pilot-model joint CRB must not improve over the unit-gain CRB.');
end

if verbose
  fprintf('  ref effective gain            : %.6g\n', auxPilotRef.effectivePilotGainAbs(1));
  fprintf('  joint effective gain          : %s\n', mat2str(auxPilotJoint.effectivePilotGainAbs(:).', 6));
  fprintf('  ref angle unit/pilot (deg)    : %.6g / %.6g\n', unitAngleRef, pilotAngleRef);
  fprintf('  joint angle unit/pilot (deg)  : %.6g / %.6g\n', unitAngleJoint, pilotAngleJoint);
end
fprintf('PASS: regressionSfDoaOnlyPilotModelCrb\n');

end

function localAssertFiniteCrb(crb, caseTag)
if ~isequal(size(crb), [2, 2]) || any(~isfinite(crb(:)))
  error('regressionSfDoaOnlyPilotModelCrb:InvalidCrb', ...
    '%s CRB must be a finite 2x2 matrix.', caseTag);
end
end

function localAssertUnitGain(auxUnit, caseTag)
gainAbs = reshape(auxUnit.effectivePilotGainAbs, [], 1);
if any(abs(gainAbs - 1) > 1e-12)
  error('regressionSfDoaOnlyPilotModelCrb:UnitGainMismatch', ...
    '%s unit-gain mode must report effective gain equal to one.', caseTag);
end
end

function localAssertFimScaling(auxUnit, auxPilot, caseTag)
unitTrace = cellfun(@(x) trace(x), auxUnit.fimCell(:));
pilotTrace = cellfun(@(x) trace(x), auxPilot.fimCell(:));
gainSq = reshape(auxPilot.effectivePilotGainAbs(:).^2, [], 1);
ratio = pilotTrace(:) ./ unitTrace(:);
tol = max(1e-10, 1e-8 * max(abs(gainSq(:))));
if any(abs(ratio(:) - gainSq(:)) > tol)
  error('regressionSfDoaOnlyPilotModelCrb:FimScaleMismatch', ...
    '%s per-sat FIM scaling does not match effectiveGainAbs^2.', caseTag);
end
if any(gainSq(:) > 1 + 1e-12) || any(gainSq(:) < -1e-12)
  error('regressionSfDoaOnlyPilotModelCrb:InvalidGainRange', ...
    '%s pilot-model effective gains must stay within [0, 1].', caseTag);
end
end

function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact two-satellite static scene.

rng(253);
numUsr = 1;
numSym = 128;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
usrLla = [37.78; 36.59; 0];
utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
scene = genMultiSatScene(utc, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1);
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct('rolloff', 0.25, 'span', 8);
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
[~, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);

fixture = struct();
fixture.scene = scene;
fixture.truth = struct();
fixture.truth.latlonTrueDeg = usrLla(1:2, 1);
fixture.truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
fixture.truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.refSatIdxLocal = refSatIdxLocal;
end
