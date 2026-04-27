% Strategy sweep for CP-U DoA box width.
% This script keeps the same CP-K seed and compares several CP-U DoA local
% boxes so we can choose a practical default width. It reports the final
% angle error, fval, boundary margin, and fdRate gap for each box.
clear(); close all; clc;
fixture = localBuildRegressionFixture();

fprintf('Running sweepMfUnknownDoaBox ...\n');

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

halfWidthList = [1e-8, 0.002, 0.005, 0.015, 0.03];
resultCell = cell(numel(halfWidthList), 1);
for iWidth = 1:numel(halfWidthList)
  halfWidthUse = halfWidthList(iWidth) * ones(2, 1);
  optUnknown = commonOpt;
  optUnknown.fdRateMode = 'unknown';
  optUnknown.initDoaParam = seedLatlon;
  optUnknown.initDoaHalfWidth = halfWidthUse;
  [modelUse, lbUse, ubUse] = buildDoaDopplerMfModel( ...
    fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, ...
    fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, optUnknown);
  [solveUse, optimUse] = solveDoaDopplerMfBranches(modelUse, initParamRelease, false);
  traceUse = localBuildDoaConstraintTrace(lbUse, ubUse, initParamRelease, solveUse.optVar(:));
  resultCell{iWidth} = localBuildSweepRow(halfWidthUse, solveUse, optimUse, traceUse, fixture.truth.latlonTrueDeg(:), fixture.truth.fdRateFit);
end

fprintf('  seed DoA (deg)                : %s\n', localFormatNumericRow(seedLatlon));
for iWidth = 1:numel(resultCell)
  row = resultCell{iWidth};
  fprintf('  -- half width %.6f deg --\n', row.halfWidthDeg);
  fprintf('     DoA width total            : %s\n', localFormatNumericRow(row.totalWidthDeg));
  fprintf('     final DoA                  : %s\n', localFormatNumericRow(row.finalDoaDeg));
  fprintf('     final margin to bounds     : %s\n', localFormatNumericRow(row.finalMarginDeg));
  fprintf('     final angle err (deg)      : %.6f\n', row.angleErrDeg);
  fprintf('     final fdRate gap (Hz/s)    : %.6f\n', row.fdRateGapHzPerSec);
  fprintf('     fval                       : %.6e\n', row.fval);
  fprintf('     solveVariant               : %s\n', row.solveVariant);
end

angleErrList = cellfun(@(s) s.angleErrDeg, resultCell);
fvalList = cellfun(@(s) s.fval, resultCell);
[~, bestAngleIdx] = min(angleErrList);
[~, bestFvalIdx] = min(fvalList);

fprintf('  best angle width (deg)        : %.6f\n', resultCell{bestAngleIdx}.halfWidthDeg);
fprintf('  best fval width (deg)         : %.6f\n', resultCell{bestFvalIdx}.halfWidthDeg);
fprintf('DONE: sweepMfUnknownDoaBox\n');

function row = localBuildSweepRow(halfWidthUse, solveUse, optimUse, traceUse, doaTrue, fdRateTrue)
row = struct();
row.halfWidthDeg = halfWidthUse(1);
row.totalWidthDeg = traceUse.width;
row.finalDoaDeg = traceUse.finalVal;
row.finalMarginDeg = traceUse.finalMargin;
row.angleErrDeg = calcLatlonAngleError(traceUse.finalVal(:), doaTrue);
row.fdRateGapHzPerSec = abs(solveUse.optVar(end) - fdRateTrue);
row.fval = solveUse.fval;
row.solveVariant = string(optimUse.solveVariant);
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
