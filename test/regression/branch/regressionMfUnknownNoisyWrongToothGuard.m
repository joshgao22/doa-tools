function regressionMfUnknownNoisyWrongToothGuard(varargin)
% Regression check for one noisy MF CP-U guard against wrong-tooth collapse.
% This guard should follow the same fast subset-recovery flow used by the
% dynamic dev/perf scripts: curated primary subsets first, then one fixed
% curated rescue bank, and finally a conditional periodic-wide rescue on the
% periodic full-data window. The goal is to keep the noisy wrong-tooth case
% pinned to the correct tooth neighborhood without relying on true-random
% fallback schedules.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

fprintf('Running regressionMfUnknownNoisyWrongToothGuard ...%s', newline);
  fprintf('  verbose trace enabled for noisy wrong-tooth guard regression.%s', newline);

flowOpt = fixture.flowOpt;
flowOpt.verbose = verbose;
bundle = buildDoaDopplerDynamicTransitionBundle( ...
  fixture.periodicFixture, fixture.subsetFixtureCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, false, flowOpt);

caseUnknown = bundle.caseResult(8);
estResultUnknown = caseUnknown.estResult;
optimInfoUnknown = localGetFieldOrDefault(estResultUnknown, 'optimInfo', struct());
solveVariantUnknown = string(localGetFieldOrDefault(estResultUnknown, 'solveVariant', ""));
latlonUnknown = reshape(estResultUnknown.doaParamEst, 2, []);
angleErrUnknownDeg = calcLatlonAngleError(latlonUnknown, fixture.truth.latlonTrueDeg(:));
fdRefErrUnknownHz = abs(estResultUnknown.fdRefEst - fixture.truth.fdRefFit);
fdRateErrUnknownHzPerSec = abs(estResultUnknown.fdRateEst - fixture.truth.fdRateFit);
selectedToothIdx = localGetFieldOrDefault(bundle.selectedSubsetSummary, 'toothIdx', NaN);
selectedToothResidualHz = abs(localGetFieldOrDefault(bundle.selectedSubsetSummary, 'toothResidualHz', NaN));
selectedSubsetLabel = string(localGetFieldOrDefault(bundle, 'selectedSubsetLabel', ""));
selectedFinalTag = string(localGetFieldOrDefault(bundle, 'selectedFinalTag', ""));

if ~estResultUnknown.isResolved
  error('regressionMfUnknownNoisyWrongToothGuard:UnknownDidNotResolve', ...
    'CP-U must remain resolved on the noisy regression fixture.');
end
if strlength(solveVariantUnknown) <= 0
  error('regressionMfUnknownNoisyWrongToothGuard:MissingSolveVariant', ...
    'The final CP-U result must expose a valid solveVariant.');
end
if angleErrUnknownDeg > 0.05
  error('regressionMfUnknownNoisyWrongToothGuard:LargeAngleError', ...
    'CP-U angle error is too large on the noisy regression fixture.');
end
if fdRefErrUnknownHz >= 0.5 * fixture.toothStepHz
  error('regressionMfUnknownNoisyWrongToothGuard:WrongToothSelected', ...
    'CP-U fdRef landed outside the correct tooth neighborhood.');
end
if fdRateErrUnknownHzPerSec > 500
  error('regressionMfUnknownNoisyWrongToothGuard:LargeFdRateError', ...
    'CP-U fdRate error is too large on the noisy regression fixture.');
end
if ~(selectedFinalTag == "periodic-wide" || selectedFinalTag == "periodic-fd-anchor" || ...
    selectedFinalTag == "periodic-in-tooth")
  error('regressionMfUnknownNoisyWrongToothGuard:UnexpectedFinalTag', ...
    'The noisy flow selected an unexpected final branch tag.');
end
if startsWith(selectedSubsetLabel, "random")
  error('regressionMfUnknownNoisyWrongToothGuard:RandomFallbackUsed', ...
    'The fast noisy guard should use fixed curated rescue subsets, not true-random fallback.');
end
if isfinite(selectedToothIdx) && abs(selectedToothIdx) > 0 && selectedToothResidualHz > 50
  error('regressionMfUnknownNoisyWrongToothGuard:BadSubsetTooth', ...
    'The selected subset winner should not stay on a clearly bad tooth.');
end

fprintf('  selected subset label    : %s%s', selectedSubsetLabel, newline);
fprintf('  selected final tag       : %s%s', selectedFinalTag, newline);
fprintf('  CP-U angle err (deg)     : %.6f%s', angleErrUnknownDeg, newline);
fprintf('  CP-U fdRef err (Hz)      : %.6f%s', fdRefErrUnknownHz, newline);
fprintf('  CP-U fdRate err (Hz/s)   : %.6f%s', fdRateErrUnknownHzPerSec, newline);
fprintf('  CP-U solveVariant        : %s%s', solveVariantUnknown, newline);
fprintf('PASS: regressionMfUnknownNoisyWrongToothGuard%s', newline);

end

function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one noisy fast-flow dynamic case.

rng(280);

numUsr = 1;
frameIntvlSec = 1 / 750;
masterOffsetIdx = -9:10;
periodicOffsetIdx = -4:5;
[primarySubsetOffsetCellFull, primarySubsetLabelListFull, rescueSubsetOffsetCellFull, rescueSubsetLabelListFull] = ...
  getDynamicCuratedSubsetBank();
% This regression only guards the noisy wrong-tooth sentinel. Keep one strong
% primary schedule plus one deterministic rescue schedule so the check stays
% representative without paying for the full curated bank.
primarySubsetOffsetCell = {primarySubsetOffsetCellFull{2}};
primarySubsetLabelList = primarySubsetLabelListFull(2);
rescueSubsetOffsetCell = {rescueSubsetOffsetCellFull{2}};
rescueSubsetLabelList = rescueSubsetLabelListFull(2);

snrDb = 10;
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
utcVecMaster = utcRef + seconds(masterOffsetIdx * frameIntvlSec);
usrLla = [37.78; 36.59; 0];

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
sceneSeqMaster = genMultiFrameScene(utcVecMaster, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1, find(masterOffsetIdx == 0, 1, 'first'));

linkParamCellMaster = cell(1, sceneSeqMaster.numFrame);
for iFrame = 1:sceneSeqMaster.numFrame
  linkParamCellMaster{iFrame} = getLinkParam(sceneSeqMaster.sceneCell{iFrame}, wavelen);
end

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.spatial.model = 'dynamic';
snapOpt.spatial.refFrameIdx = sceneSeqMaster.refFrameIdx;
snapOpt.phase.timeModel = 'global';
snapOpt.phase.frameModel = 'shared';
snapOpt.phase.sharedPhase = 2 * pi * rand(sceneSeqMaster.numSat, sceneSeqMaster.numUser);
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
snapOpt.precomp.linkParamCell = linkParamCellMaster;

pwrNoise = 1 / (10^(snrDb / 10));
[rxSigCellMaster, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeqMaster, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, pwrNoise, repmat({ones(sceneSeqMaster.numSat, sceneSeqMaster.numUser)}, 1, sceneSeqMaster.numFrame), snapOpt);

gridSize = [41, 41];
searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
fdRangeDefault = [-2e5, 2e5];
fdRateRangeDefault = [-1e4, 0];

[periodicFixture, subsetFixtureCell] = buildDynamicRepeatFixtures( ...
  sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, masterOffsetIdx, periodicOffsetIdx, ...
  primarySubsetOffsetCell, primarySubsetLabelList, 0, 280, gridSize, searchRange, E, ...
  wavelen, waveInfo.sampleRate, fdRangeDefault, fdRateRangeDefault, struct());
truth = periodicFixture.truth;
truth.latlonTrueDeg = usrLla(1:2, 1);

flowOpt = struct();
flowOpt.enableWeightSweep = false;
flowOpt.weightSweepAlpha = zeros(0, 1);
flowOpt.staticMsHalfWidth = [0.002; 0.002];
flowOpt.doaOnlyOpt = struct('useLogObjective', true);
flowOpt.staticBaseOpt = struct('useLogObjective', true);
flowOpt.dynBaseOpt = struct( ...
  'useLogObjective', true, ...
  'initFdCount', 81, ...
  'useAccessMask', false, ...
  'phaseMode', 'continuous', ...
  'steeringMode', 'framewise', ...
  'continuousPhaseConsistencyWeight', 0.05, ...
  'continuousPhaseCollapsePenaltyWeight', 0.10, ...
  'continuousPhaseNegativeProjectionPenaltyWeight', 0.10, ...
  'unknownWarmAnchorUseScaledSolve', true, ...
  'unknownWarmAnchorFallbackSqp', true, ...
  'debugEnable', true, ...
  'debugStoreEvalTrace', false, ...
  'debugMaxEvalTrace', 120);
flowOpt.msContinuousPhaseConsistencyWeight = flowOpt.dynBaseOpt.continuousPhaseConsistencyWeight;
flowOpt.msContinuousPhaseCollapsePenaltyWeight = flowOpt.dynBaseOpt.continuousPhaseCollapsePenaltyWeight;
flowOpt.msContinuousPhaseNegativeProjectionPenaltyWeight = flowOpt.dynBaseOpt.continuousPhaseNegativeProjectionPenaltyWeight;
flowOpt.unknownWarmAnchorUseScaledSolve = true;
flowOpt.unknownWarmAnchorFallbackSqp = true;
flowOpt.refKnownDoaHalfWidth = [0.005; 0.005];
flowOpt.refUnknownDoaHalfWidth = [0.003; 0.003];
flowOpt.msKnownDoaHalfWidth = [0.003; 0.003];
flowOpt.msUnknownDoaHalfWidth = [0.002; 0.002];
flowOpt.subsetSelectDoaHalfWidthDeg = [0.01; 0.01];
flowOpt.subsetAnchorFdHalfWidthHz = 50;
flowOpt.subsetAnchorFdRateHalfWidthHzPerSec = 100;
flowOpt.subsetAnchorDoaHalfWidthDeg = [1e-8; 1e-8];
flowOpt.subsetAnchorFreezeDoa = true;
flowOpt.fdWideFreezeDoa = true;
flowOpt.wideRefineFdHalfWidthHz = 100;
flowOpt.wideRefineFdRateHalfWidthHzPerSec = 200;
flowOpt.wideRefineDoaHalfWidthDeg = [0.003; 0.003];
flowOpt.wideRefineDisableUnknownDoaReleaseFloor = true;
flowOpt.inToothFdHalfWidthHz = 50;
flowOpt.inToothFdRateHalfWidthHzPerSec = 100;
flowOpt.inToothDoaHalfWidthDeg = [2e-4; 2e-4];
flowOpt.inToothDisableUnknownDoaReleaseFloor = true;
flowOpt.enableUnknownWideBranches = false;
flowOpt.enableUnknownWideRefine = false;
flowOpt.deferUnknownWideBranches = false;
flowOpt.enableFastWideFallback = false;
flowOpt.enableFastFdWideFallback = false;
flowOpt.enableFastSubsetEscalation = false;
flowOpt.fastRescueSubsetOffsetCell = rescueSubsetOffsetCell;
flowOpt.fastRescueSubsetLabelList = rescueSubsetLabelList;
flowOpt.fastNumRandomSubsetTrialFallback = 0;
flowOpt.fastToothIdxThreshold = 0;
flowOpt.fastToothResidualHzThreshold = 50;
flowOpt.fastUnknownFdDriftHzThreshold = 1000;
flowOpt.fastUnknownDoaDriftDegThreshold = 0.005;
flowOpt.fastWideToothIdxThreshold = 0;
flowOpt.fastWideToothResidualHzThreshold = 250;
flowOpt.fastWideFdDriftHzThreshold = 2500;
flowOpt.fastWideDoaDriftDegThreshold = 0.01;
flowOpt.enableSubsetInToothRefine = false;
flowOpt.enableInToothCentralSkip = true;
flowOpt.inToothCentralResidualTolHz = 5;
flowOpt.enableSubsetCaseSeedReuse = true;
flowOpt.fastRandomSeedOffset = 0;
flowOpt.enableConditionalRandomSubsetRescue = false;
flowOpt.parallelOpt = struct('enableSubsetEvalParfor', false, 'enableFixtureBankParfor', false);
flowOpt.finalSelectMaxFdRefDriftHz = 100;
flowOpt.finalSelectMaxFdRateDriftHzPerSec = 100;
flowOpt.finalSelectMaxDoaDriftDeg = 0.01;
flowOpt.finalSelectMinObjGainToLeaveAnchor = 1e4;
flowOpt.finalSelectMinResidualGainToLeaveAnchor = 1e3;

toothStepHz = 1 / median(diff(periodicFixture.sceneSeq.timeOffsetSec));
fixture = struct();
fixture.periodicFixture = periodicFixture;
fixture.subsetFixtureCell = subsetFixtureCell;
fixture.truth = truth;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.toothStepHz = toothStepHz;
fixture.flowOpt = flowOpt;
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field with a default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
