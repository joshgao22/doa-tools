function regressionMfInitFdLine(varargin)
% Regression check for the MF fd-line initializer.
% This script isolates buildDoaDopplerMfInit and checks one narrow contract:
% with a trusted DoA anchor, the reference-satellite fd-line initializer
% should stay near the dynamic truth, should avoid range clipping, and
% should give nearly the same answer for the ref-only and multi-satellite
% views because the current initializer is reference-satellite-only.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

  fprintf('Running regressionMfInitFdLine ...\n');

initFdCountList = [65, 81];
resultCell = cell(numel(initFdCountList), 1);

for iCase = 1:numel(initFdCountList)
  initFdCount = initFdCountList(iCase);
  resultCell{iCase} = localRunOneInitCase(fixture, initFdCount);
  localCheckInitCase(resultCell{iCase}, fixture.truth, fixture.fdRange, fixture.fdRateRange);
end

  for iCase = 1:numel(resultCell)
    resultCur = resultCell{iCase};
    fprintf('  initFdCount = %d\n', resultCur.initFdCount);
    fprintf('    ref-only fdRefInit (Hz) : %.6f\n', resultCur.fdRefInitRefOnly);
    fprintf('    multi-sat fdRefInit (Hz): %.6f\n', resultCur.fdRefInitMulti);
    fprintf('    ref-only fdRateInit     : %.6f\n', resultCur.fdRateInitRefOnly);
    fprintf('    multi-sat fdRateInit    : %.6f\n', resultCur.fdRateInitMulti);
  end
  fprintf('PASS: regressionMfInitFdLine\n');



end

function result = localRunOneInitCase(fixture, initFdCount)
%LOCALRUNONEINITCASE Run one pair of ref-only / multi-sat fd-line inits.

commonOpt = struct();
commonOpt.useLogObjective = true;
commonOpt.phaseMode = 'continuous';
commonOpt.fdRateMode = 'unknown';
commonOpt.steeringMode = 'framewise';
commonOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
commonOpt.initFdCount = initFdCount;
commonOpt.initDoaParam = fixture.truth.latlonTrueDeg(:);
commonOpt.debugTruth = struct('fdRef', fixture.truth.fdRefFit, 'fdRate', fixture.truth.fdRateFit);

[modelRefOnly, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeqRefOnly, fixture.rxSigMfRefOnly, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewRefOnly.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, commonOpt);
[initParamRefOnly, initDiagRefOnly] = buildDoaDopplerMfInit(modelRefOnly, []);

[modelMulti, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, commonOpt);
[initParamMulti, initDiagMulti] = buildDoaDopplerMfInit(modelMulti, []);

result = struct();
result.initFdCount = initFdCount;
result.initParamRefOnly = initParamRefOnly;
result.initParamMulti = initParamMulti;
result.fdRefInitRefOnly = initDiagRefOnly.fdLine.fdRefInit;
result.fdRateInitRefOnly = initDiagRefOnly.fdLine.fdRateInit;
result.fdRefInitMulti = initDiagMulti.fdLine.fdRefInit;
result.fdRateInitMulti = initDiagMulti.fdLine.fdRateInit;
result.fdGridStepHzRefOnly = initDiagRefOnly.fdLine.fdGridStepHz;
result.fdGridStepHzMulti = initDiagMulti.fdLine.fdGridStepHz;
end


function localCheckInitCase(result, truth, fdRange, fdRateRange)
%LOCALCHECKINITCASE Check one fd-line initialization result.

fdRefTolHz = max([2 * result.fdGridStepHzRefOnly, 2 * result.fdGridStepHzMulti, 50]);
fdRateTolHzPerSec = max(0.05 * diff(fdRateRange), 500);
fdRefCrossTolHz = max([0.5 * result.fdGridStepHzRefOnly, 0.5 * result.fdGridStepHzMulti, 5]);
fdRateCrossTolHzPerSec = 250;
rangeTol = 1e-9;

if abs(result.fdRefInitRefOnly - truth.fdRefFit) > fdRefTolHz
  error('regressionMfInitFdLine:RefOnlyFdRefInitMismatch', ...
    ['ref-only fdRefInit %.3f Hz is too far from truth %.3f Hz ', ...
     '(tol %.3f Hz).'], result.fdRefInitRefOnly, truth.fdRefFit, fdRefTolHz);
end
if abs(result.fdRefInitMulti - truth.fdRefFit) > fdRefTolHz
  error('regressionMfInitFdLine:MultiFdRefInitMismatch', ...
    ['multi-sat fdRefInit %.3f Hz is too far from truth %.3f Hz ', ...
     '(tol %.3f Hz).'], result.fdRefInitMulti, truth.fdRefFit, fdRefTolHz);
end
if abs(result.fdRateInitRefOnly - truth.fdRateFit) > fdRateTolHzPerSec
  error('regressionMfInitFdLine:RefOnlyFdRateInitMismatch', ...
    ['ref-only fdRateInit %.3f Hz/s is too far from truth %.3f Hz/s ', ...
     '(tol %.3f Hz/s).'], result.fdRateInitRefOnly, truth.fdRateFit, fdRateTolHzPerSec);
end
if abs(result.fdRateInitMulti - truth.fdRateFit) > fdRateTolHzPerSec
  error('regressionMfInitFdLine:MultiFdRateInitMismatch', ...
    ['multi-sat fdRateInit %.3f Hz/s is too far from truth %.3f Hz/s ', ...
     '(tol %.3f Hz/s).'], result.fdRateInitMulti, truth.fdRateFit, fdRateTolHzPerSec);
end

if abs(result.fdRefInitRefOnly - result.fdRefInitMulti) > fdRefCrossTolHz
  error('regressionMfInitFdLine:CrossViewFdRefMismatch', ...
    'ref-only and multi-sat fdRefInit disagree more than %.3f Hz.', fdRefCrossTolHz);
end
if abs(result.fdRateInitRefOnly - result.fdRateInitMulti) > fdRateCrossTolHzPerSec
  error('regressionMfInitFdLine:CrossViewFdRateMismatch', ...
    'ref-only and multi-sat fdRateInit disagree more than %.3f Hz/s.', fdRateCrossTolHzPerSec);
end

if abs(result.fdRefInitRefOnly - fdRange(1)) <= rangeTol || abs(result.fdRefInitRefOnly - fdRange(2)) <= rangeTol
  error('regressionMfInitFdLine:RefOnlyFdRefClipped', ...
    'ref-only fdRefInit is clipped to the fdRange boundary.');
end
if abs(result.fdRefInitMulti - fdRange(1)) <= rangeTol || abs(result.fdRefInitMulti - fdRange(2)) <= rangeTol
  error('regressionMfInitFdLine:MultiFdRefClipped', ...
    'multi-sat fdRefInit is clipped to the fdRange boundary.');
end
if abs(result.fdRateInitRefOnly - fdRateRange(1)) <= rangeTol || ...
    abs(result.fdRateInitRefOnly - fdRateRange(2)) <= rangeTol
  error('regressionMfInitFdLine:RefOnlyFdRateClipped', ...
    'ref-only fdRateInit is clipped to the fdRateRange boundary.');
end
if abs(result.fdRateInitMulti - fdRateRange(1)) <= rangeTol || ...
    abs(result.fdRateInitMulti - fdRateRange(2)) <= rangeTol
  error('regressionMfInitFdLine:MultiFdRateClipped', ...
    'multi-sat fdRateInit is clipped to the fdRateRange boundary.');
end
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one high-SNR dynamic MF init case.

rng(253);

numUsr = 1;
numFrame = 9;
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
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);

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

sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal);
rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame');
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, [41, 41], searchRange, E, ...
  struct('sceneSeq', sceneSeqRefOnly));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = expandRangeToTruth([-1e4, 0], ...
  [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.sceneSeqRefOnly = sceneSeqRefOnly;
fixture.truth = truth;
fixture.viewMs = viewMs;
fixture.viewRefOnly = viewRefOnly;
fixture.rxSigCell = rxSigCell;
fixture.rxSigMfRefOnly = rxSigMfRefOnly;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
end
