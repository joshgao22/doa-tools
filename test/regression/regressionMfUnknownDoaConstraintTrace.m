% Trace script for CP-U DoA constraint behavior.
% This script prints the default unknown-branch DoA bounds, the near-zero
% fd-only bounds, and the wider DoA-release bounds, then compares where the
% seed and final solutions sit relative to those bounds.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionMfUnknownDoaConstraintTrace ...\n');

commonOpt = struct();
commonOpt.useLogObjective = false;
commonOpt.phaseMode = 'continuous';
commonOpt.steeringMode = 'framewise';
commonOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
commonOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.05; -0.04];
commonOpt.initDoaHalfWidth = [0.10; 0.10];
commonOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

fdRateBiasHzPerSec = 150;
fdRateSeed = min(max(fixture.truth.fdRateFit + fdRateBiasHzPerSec, fixture.fdRateRange(1)), fixture.fdRateRange(2));

optKnown = commonOpt; optKnown.fdRateMode = 'known'; optKnown.fdRateKnown = fdRateSeed;
[modelKnown, ~, ~] = buildDoaDopplerMfModel(fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, optKnown);
[initParamKnown, ~] = buildDoaDopplerMfInit(modelKnown, []);
[solveKnown, ~] = solveDoaDopplerMfBranches(modelKnown, initParamKnown, false);
seedLatlon = reshape(solveKnown.finalEvalDiag.latlon, 2, 1);
initParamRelease = [solveKnown.optVar(:); fdRateSeed];

optUnknownDefault = commonOpt; optUnknownDefault.fdRateMode = 'unknown'; optUnknownDefault.initDoaParam = seedLatlon; optUnknownDefault.initDoaHalfWidth = [0.002; 0.002];
[modelDefault, lbDefault, ubDefault] = buildDoaDopplerMfModel(fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, optUnknownDefault);
[solveDefault, optimDefault] = solveDoaDopplerMfBranches(modelDefault, initParamRelease, false);

optUnknownFdOnly = commonOpt; optUnknownFdOnly.fdRateMode = 'unknown'; optUnknownFdOnly.initDoaParam = seedLatlon; optUnknownFdOnly.initDoaHalfWidth = [1e-8; 1e-8];
[modelFdOnly, lbFdOnly, ubFdOnly] = buildDoaDopplerMfModel(fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, optUnknownFdOnly);
[solveFdOnly, optimFdOnly] = solveDoaDopplerMfBranches(modelFdOnly, initParamRelease, false);

optUnknownDoa = commonOpt; optUnknownDoa.fdRateMode = 'unknown'; optUnknownDoa.initDoaParam = seedLatlon; optUnknownDoa.initDoaHalfWidth = [0.03; 0.03];
[modelDoa, lbDoa, ubDoa] = buildDoaDopplerMfModel(fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, optUnknownDoa);
[solveDoa, optimDoa] = solveDoaDopplerMfBranches(modelDoa, initParamRelease, false);

traceDefault = localBuildDoaConstraintTrace(lbDefault, ubDefault, initParamRelease, solveDefault.optVar(:));
traceFdOnly = localBuildDoaConstraintTrace(lbFdOnly, ubFdOnly, initParamRelease, solveFdOnly.optVar(:));
traceDoa = localBuildDoaConstraintTrace(lbDoa, ubDoa, initParamRelease, solveDoa.optVar(:));

fprintf('  seed DoA (deg)                : %s\n', localFormatNumericRow(seedLatlon));
localPrintTrace('default unknown', traceDefault, solveDefault, fixture.truth.latlonTrueDeg(:), optimDefault);
localPrintTrace('fd-only unknown', traceFdOnly, solveFdOnly, fixture.truth.latlonTrueDeg(:), optimFdOnly);
localPrintTrace('DoA-release unknown', traceDoa, solveDoa, fixture.truth.latlonTrueDeg(:), optimDoa);
fprintf('DONE: regressionMfUnknownDoaConstraintTrace\n');

function localPrintTrace(labelText, traceInfo, solveInfo, doaTrue, optimInfo)
fprintf('  -- %s --\n', labelText);
fprintf('     DoA lb                     : %s\n', localFormatNumericRow(traceInfo.lb));
fprintf('     DoA ub                     : %s\n', localFormatNumericRow(traceInfo.ub));
fprintf('     DoA width                  : %s\n', localFormatNumericRow(traceInfo.width));
fprintf('     init DoA                   : %s\n', localFormatNumericRow(traceInfo.initVal));
fprintf('     final DoA                  : %s\n', localFormatNumericRow(traceInfo.finalVal));
fprintf('     init margin to bounds      : %s\n', localFormatNumericRow(traceInfo.initMargin));
fprintf('     final margin to bounds     : %s\n', localFormatNumericRow(traceInfo.finalMargin));
fprintf('     final angle err (deg)      : %.6f\n', calcLatlonAngleError(traceInfo.finalVal(:), doaTrue));
fprintf('     solveVariant               : %s\n', string(optimInfo.solveVariant));
fprintf('     fval                       : %.6e\n', solveInfo.fval);
end

function traceInfo = localBuildDoaConstraintTrace(lb, ub, initParam, optVar)
idxDoa = 1:2;
traceInfo = struct();
traceInfo.lb = reshape(lb(idxDoa), 1, []);
traceInfo.ub = reshape(ub(idxDoa), 1, []);
traceInfo.width = traceInfo.ub - traceInfo.lb;
traceInfo.initVal = reshape(initParam(idxDoa), 1, []);
traceInfo.finalVal = reshape(optVar(idxDoa), 1, []);
traceInfo.initMargin = min(traceInfo.initVal - traceInfo.lb, traceInfo.ub - traceInfo.initVal);
traceInfo.finalMargin = min(traceInfo.finalVal - traceInfo.lb, traceInfo.ub - traceInfo.finalVal);
end

function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false), ', '), ']'];
end

function fixture = localBuildRegressionFixture()
rng(253);
numUsr = 1; numFrame = 7; refFrameIdx = ceil(numFrame / 2); frameIntvlSec = 1 / 750;
sampleRate = 122.88e6; symbolRate = 3.84e6; osf = sampleRate / symbolRate; numSym = 32; carrierFreq = 11.7e9; lightSpeed = 299792458; wavelen = lightSpeed / carrierFreq; E = referenceEllipsoid('sphere');
utcRef = datetime([2026, 03, 18, 17, 08, 0], 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec); usrLla = [37.78; 36.59; 0];
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle")); arrUpa = createUpa([4, 4], wavelen / 2);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, 15, 55, "satellite", 1, refFrameIdx); sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
linkParamCell = cell(1, sceneSeq.numFrame); for iFrame = 1:sceneSeq.numFrame, linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen); end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1); truth.latlonTrueDeg = usrLla(1:2, 1);
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1); pulseOpt = struct(); pulseOpt.rolloff = 0.25; pulseOpt.span = 8; [pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
snapOpt = struct(); snapOpt.spatial.model = 'dynamic'; snapOpt.spatial.refFrameIdx = sceneSeq.refFrameIdx; snapOpt.phase.timeModel = 'global'; snapOpt.phase.frameModel = 'shared'; snapOpt.phase.sharedPhase = 2 * pi * rand(sceneSeq.numSat, sceneSeq.numUser); snapOpt.wave.delayModel = 'phaseOnly'; snapOpt.wave.timeRef = 'zero'; snapOpt.wave.carrierPhaseModel = 'none'; snapOpt.precomp.linkParamCell = linkParamCell;
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, 0, repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame), snapOpt);
searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [41, 41], searchRange, E, struct('sceneSeq', sceneSeq));
fdRange = localExpandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4); fdRateRange = localExpandRangeToTruth([-1e4, 0], [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);
fixture = struct(); fixture.sceneSeq = sceneSeq; fixture.truth = truth; fixture.viewMs = viewMs; fixture.rxSigCell = rxSigCell; fixture.pilotWave = pilotWave; fixture.sampleRate = waveInfo.sampleRate; fixture.carrierFreq = carrierFreq; fixture.fdRange = fdRange; fixture.fdRateRange = fdRateRange;
end

function value = localExpandRangeToTruth(defaultRange, truthVals, fracPad, absPad)
truthVals = truthVals(isfinite(truthVals)); if isempty(truthVals), value = defaultRange; return; end
minTruth = min(truthVals); maxTruth = max(truthVals); spanTruth = max(maxTruth - minTruth, eps); pad = max(absPad, fracPad * spanTruth); value = [min(defaultRange(1), minTruth - pad), max(defaultRange(2), maxTruth + pad)];
end

function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath')); projectRoot = fileparts(fileparts(scriptDir)); addpath(genpath(projectRoot));
end
