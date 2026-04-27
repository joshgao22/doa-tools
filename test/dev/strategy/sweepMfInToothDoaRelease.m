% Strategy sweep for in-tooth DoA release after subset tooth selection.
% The goal is to isolate the second-stage local basin inside the correct
% tooth after non-periodic subset selection has already removed the obvious
% 1/T_f comb ambiguity.
% Steps:
%   1) Run the periodic full-data baseline.
%   2) Run one deterministic and several random non-periodic subsets.
%   3) Select the best subset that lands on the truth tooth.
%   4) Return to the periodic full-data window, keep fdRef/fdRate inside a
%      narrow in-tooth box around the subset result, and sweep the DoA
%      half-width.
% This answers: once the tooth is fixed, how much DoA release can we allow
% before the periodic full-data refinement falls into an in-tooth bad basin?
clear(); close all; clc;

periodicOffsetIdx = -4:5;
masterOffsetIdx = -9:10;
deterministicOffsetIdx = [-9, -7, -6, -5, -4, 0, 1, 3, 7, 8];
numRandomTrial = 24;
subsetSelectDoaHalfWidthDeg = [0.01; 0.01];
fdOnlyHalfWidthHz = 50;
fdOnlyHalfWidthHzPerSec = 100;
doaHalfWidthListDeg = [1e-8, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2];
rng(253);

fixturePeriodic = localBuildPipelineFixture(periodicOffsetIdx);
fixtureMaster = localBuildPipelineFixture(masterOffsetIdx);
toothStepHz = 1 / fixturePeriodic.frameIntvlSec;

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;
staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;
dynBaseOpt = struct();
dynBaseOpt.useLogObjective = true;
dynBaseOpt.initFdCount = 81;
dynBaseOpt.useAccessMask = false;
dynBaseOpt.phaseMode = 'continuous';
dynBaseOpt.steeringMode = 'framewise';
dynBaseOpt.steeringRefFrameIdx = fixturePeriodic.sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = false;
dynBaseOpt.enableFdAliasUnwrap = true;
weightSweepAlpha = [0, 0.25, 0.5, 1];


baselineCase = localRunMsUnknownCase( ...
  fixturePeriodic, fixturePeriodic.fdRange, fixturePeriodic.fdRateRange, struct(), struct());

subsetBankParallelOpt = struct();
subsetBankParallelOpt.enableParfor = false;
subsetBankParallelOpt.enableFixtureBankParfor = false;
subsetBankParallelOpt.minFixtureBankForParfor = inf;
subsetFixtureCell = buildDynamicSubsetFixtureBank( ...
  fixtureMaster.sceneSeq, fixtureMaster.linkParamCell, fixtureMaster.rxSigCell, ...
  fixtureMaster.masterOffsetIdx, deterministicOffsetIdx, numRandomTrial, 253, ...
  fixtureMaster.gridSize, fixtureMaster.searchRange, fixtureMaster.E, ...
  fixtureMaster.wavelen, fixtureMaster.sampleRate, fixtureMaster.fdRange, ...
  fixtureMaster.fdRateRange, subsetBankParallelOpt);

helperParallelOpt = struct();
helperParallelOpt.enableParfor = false;
helperParallelOpt.enableSubsetSolveParfor = false;
helperParallelOpt.minSubsetBankForParfor = inf;
helperParallelOpt.forceSubsetLabel = "deterministic";

flow = runMsUnknownInToothSweepFlow( ...
  fixturePeriodic, subsetFixtureCell, baselineCase, fixturePeriodic.truth, ...
  fixturePeriodic.pilotWave, fixturePeriodic.carrierFreq, fixturePeriodic.sampleRate, false, weightSweepAlpha(:), ...
  doaOnlyOpt, staticBaseOpt, dynBaseOpt, [0.01; 0.01], ...
  subsetSelectDoaHalfWidthDeg(:), fdOnlyHalfWidthHz, fdOnlyHalfWidthHzPerSec, ...
  doaHalfWidthListDeg(:), toothStepHz, helperParallelOpt);

baselineSummary = flow.wideSummary;
subsetSummary = flow.subsetSummary;
fdRangeInTooth = flow.fdRangeInTooth;
fdRateRangeInTooth = flow.fdRateRangeInTooth;
bestAngleSummary = flow.bestAngleSummary;
bestStableSummary = flow.bestStableSummary;

fprintf('Running sweepMfInToothDoaRelease\n');
fprintf('  tooth step (Hz)                    : %.6f\n', toothStepHz);
localPrintCaseSummary('periodic wide baseline', baselineSummary);
fprintf('  selected subset label              : %s\n', flow.selectedSubsetLabel);
fprintf('  selected subset offsets            : %s\n', localFormatIntegerRow(flow.selectedSubsetOffsets));
localPrintCaseSummary('selected subset', subsetSummary);
fprintf('  in-tooth periodic fdRange          : %s\n', localFormatNumericRow(fdRangeInTooth));
fprintf('  in-tooth periodic fdRateRange      : %s\n', localFormatNumericRow(fdRateRangeInTooth));
for iWidth = 1:height(flow.sweepSummaryTable)
  s = localSummaryFromSweepTable(flow.sweepSummaryTable(iWidth, :));
  localPrintCaseSummary(sprintf('DoA half width %.6g deg', s.doaHalfWidthDeg), s);
end
fprintf('  best angle width (deg)             : %.6g\n', flow.bestAngleWidthDeg);
localPrintCaseSummary('best-angle summary', bestAngleSummary);
fprintf('  best stable width (deg)            : %.6g\n', flow.bestStableWidthDeg);
localPrintCaseSummary('best-stable summary', bestStableSummary);
localPrintImprovement('best-angle vs baseline', baselineSummary, bestAngleSummary);
localPrintImprovement('best-stable vs baseline', baselineSummary, bestStableSummary);
fprintf('PASS: sweepMfInToothDoaRelease\n');

function [summary, caseUse] = localRunSubsetSchedule(offsetIdxVec, toothStepHz, doaHalfWidthDeg)
fixtureUse = localBuildPipelineFixture(offsetIdxVec);
unknownOverride = struct();
unknownOverride.initDoaHalfWidth = doaHalfWidthDeg(:);
caseUse = localRunMsUnknownCase(fixtureUse, fixtureUse.fdRange, fixtureUse.fdRateRange, struct(), unknownOverride);
summary = localBuildCaseSummary(caseUse, fixtureUse.truth, toothStepHz);
end


function [bestAngleIdx, bestStableIdx] = localSelectBestSweep(summaryCell)
numCase = numel(summaryCell);
scoreAngle = nan(numCase, 5);
scoreStable = nan(numCase, 5);
for iCase = 1:numCase
  s = summaryCell{iCase};
  scoreAngle(iCase, :) = [abs(s.toothIdx), abs(s.toothResidualHz), s.angleErrDeg, abs(s.fdRateErrHzPerSec), abs(s.fdRefErrHz)];
  scoreStable(iCase, :) = [abs(s.toothIdx), abs(s.toothResidualHz), abs(s.fdRefErrHz), abs(s.fdRateErrHzPerSec), s.angleErrDeg];
end
  [~, orderAngle] = sortrows(scoreAngle, [1, 2, 3, 4, 5]);
  [~, orderStable] = sortrows(scoreStable, [1, 2, 3, 4, 5]);
  bestAngleIdx = orderAngle(1);
  bestStableIdx = orderStable(1);
end


function bestIdx = localSelectBestSubset(summaryCell)
numCase = numel(summaryCell);
scoreMat = nan(numCase, 4);
for iCase = 1:numCase
  s = summaryCell{iCase};
  scoreMat(iCase, 1) = abs(s.toothIdx);
  scoreMat(iCase, 2) = abs(s.toothResidualHz);
  scoreMat(iCase, 3) = s.angleErrDeg;
  scoreMat(iCase, 4) = abs(s.fdRateErrHzPerSec);
end
[~, orderIdx] = sortrows(scoreMat, [1, 2, 3, 4]);
bestIdx = orderIdx(1);
end


function localPrintImprovement(tag, baselineSummary, refineSummary)
fprintf('  %s angle improvement    : %.6f deg\n', tag, baselineSummary.angleErrDeg - refineSummary.angleErrDeg);
fprintf('  %s fdRef improvement    : %.6f Hz\n', tag, abs(baselineSummary.fdRefErrHz) - abs(refineSummary.fdRefErrHz));
fprintf('  %s fdRate improvement   : %.6f Hz/s\n', tag, abs(baselineSummary.fdRateErrHzPerSec) - abs(refineSummary.fdRateErrHzPerSec));
end


function localPrintCaseSummary(caseName, summary)
fprintf('  -- %s --\n', caseName);
fprintf('     solveVariant                    : %s\n', summary.solveVariant);
fprintf('     isResolved                      : %d\n', summary.isResolved);
fprintf('     angle err (deg)                 : %.6f\n', summary.angleErrDeg);
fprintf('     fdRef err (Hz)                  : %.6f\n', summary.fdRefErrHz);
fprintf('     fdRate err (Hz/s)               : %.6f\n', summary.fdRateErrHzPerSec);
fprintf('     nearest tooth index             : %d\n', summary.toothIdx);
fprintf('     residual to nearest tooth (Hz)  : %.6f\n', summary.toothResidualHz);
fprintf('     runTimeMs                       : %.6f\n', summary.runTimeMs);
end


function summary = localSummaryFromSweepTable(rowUse)
summary = struct();
summary.solveVariant = string(rowUse.solveVariant);
summary.isResolved = logical(rowUse.isResolved);
summary.angleErrDeg = rowUse.angleErrDeg;
summary.fdRefErrHz = rowUse.fdRefErrHz;
summary.fdRateErrHzPerSec = rowUse.fdRateErrHzPerSec;
summary.toothIdx = rowUse.toothIdx;
summary.toothResidualHz = rowUse.toothResidualHz;
summary.runTimeMs = rowUse.runTimeMs;
summary.doaHalfWidthDeg = rowUse.doaHalfWidthDeg;
end

function summary = localBuildCaseSummary(caseUse, truth, toothStepHz)
summary = struct();
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
end


function caseDynMsUnknown = localRunMsUnknownCase(fixture, fdRangeUse, fdRateRangeUse, dynKnownOverride, dynUnknownOverride)
[viewRefOnly, viewOtherOnly] = localBuildSingleSatViews(fixture);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;
staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;
weightSweepAlpha = [0, 0.25, 0.5, 1];

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, fixture.viewMs, fixture.wavelen, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fdRangeUse, fixture.truth, ...
  fixture.otherSatIdxGlobal, false, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha, [0.01; 0.01]);

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

dynBaseOpt = struct();
dynBaseOpt.useLogObjective = true;
dynBaseOpt.initFdCount = 81;
dynBaseOpt.useAccessMask = false;
dynBaseOpt.phaseMode = 'continuous';
dynBaseOpt.steeringMode = 'framewise';
dynBaseOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = false;
dynBaseOpt.enableFdAliasUnwrap = true;

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = [0.003; 0.003];
dynMsKnownOpt = localMergeStruct(dynMsKnownOpt, dynKnownOverride);

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fdRangeUse, fdRateRangeUse, false, dynMsKnownOpt, true, ...
  fixture.debugTruthMs, buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];
dynMsUnknownOpt = localMergeStruct(dynMsUnknownOpt, dynUnknownOverride);

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);

caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fdRangeUse, fdRateRangeUse, false, dynMsUnknownOpt, false, ...
  fixture.debugTruthMs, msUnknownCand);
end


function [viewRefOnly, viewOtherOnly] = localBuildSingleSatViews(fixture)
sceneSeqRefOnly = selectSatSceneSeq(fixture.sceneSeq, fixture.refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(fixture.rxSigCell, fixture.refSatIdxLocal, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  fixture.gridSize, fixture.searchRange, fixture.E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = selectSatScene(fixture.sceneRef, fixture.otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(fixture.rxSigRef, fixture.otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  fixture.gridSize, fixture.searchRange, fixture.E);
end


function fixture = localBuildPipelineFixture(offsetIdxVec)
rng(253);

numUsr = 1;
refFrameIdx = find(offsetIdxVec == 0, 1);
if isempty(refFrameIdx)
  error('sweepMfInToothDoaRelease:MissingReferenceFrame', ...
    'offsetIdxVec must contain one zero-offset reference frame.');
end

frameIntvlSec = 1 / 750;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

snrDb = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));

usrLla = [37.78; 36.59; 0];
utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(offsetIdxVec(:).' * frameIntvlSec);

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, ~] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.utcRef = sceneSeq.utcRef;
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.usrElevationDeg = reshape(sceneRef.accessInfo.usrElevationDeg(:, 1), 1, []);
truth.refWeight = sceneRef.ref.weight(:);
truth.fdRefTrueHz = truth.fdRefSeries(refFrameIdx);
truth.fdRateTrueHzPerSec = truth.fdRateFit;
truth.fdSatTrueHz = reshape(truth.fdSatSeries(:, refFrameIdx), [], 1);
truth.deltaFdTrueHz = reshape(truth.deltaFdSeries(:, refFrameIdx), [], 1);
truth.localDoaRef = reshape(sceneRef.localDoa(:, refSatIdxLocal), 2, 1);
fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatTrueHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth([-1e4, 0], fdRateTruthCand, 0.1, 5e2);

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
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
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots( ...
  sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGainCell, snapOpt);
rxSigRef = rxSigCell{sceneSeq.refFrameIdx};

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; usrLla(2, 1) - 5, usrLla(2, 1) + 5];
gridSize = [50 50];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

debugTruthMs = localBuildTruthDebugForView(viewMs, truth);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.sceneRef = sceneRef;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
fixture.truth = truth;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
fixture.linkParamCell = linkParamCell;
fixture.masterOffsetIdx = reshape(offsetIdxVec, 1, []);
fixture.rxSigCell = rxSigCell;
fixture.rxSigRef = rxSigRef;
fixture.viewMs = viewMs;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.wavelen = wavelen;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.E = E;
fixture.frameIntvlSec = frameIntvlSec;
fixture.debugTruthMs = debugTruthMs;
end


function debugTruth = localBuildTruthDebugForView(view, truth)
sceneRef = view.sceneRef;
sceneSeq = view.sceneSeq;

debugTruth = struct();
debugTruth.doaParam = truth.latlonTrueDeg(:);
debugTruth.fdRef = truth.fdRefFit;
debugTruth.fdRate = truth.fdRateFit;
debugTruth.localDoaTrue = localExtractSceneLocalDoa(sceneSeq);
debugTruth.fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth);
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
debugTruth.refSatIdxLocal = refSatIdxLocal;
debugTruth.satIdxGlobal = reshape(sceneRef.satIdx, [], 1);
debugTruth.refWeight = reshape(sceneRef.ref.weight, [], 1);
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
numSat = sceneRef.numSat;
numFrame = size(truth.fdSatSeries, 2);
fdLocalTrue = nan(numSat, numFrame);
sceneSatIdx = reshape(sceneRef.satIdx, 1, []);
truthSatIdx = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
[isFound, loc] = ismember(sceneSatIdx, truthSatIdx);
if any(isFound)
  fdLocalTrue(isFound, :) = truth.fdSatSeries(loc(isFound), :);
end
end


function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
localDoa = sceneSeq.localDoa;
if ndims(localDoa) == 3
  truthLocalDoa = localDoa;
elseif ndims(localDoa) == 4
  truthLocalDoa = localDoa(:, :, 1, :);
else
  error('sweepMfInToothDoaRelease:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end
end


function out = localMergeStruct(base, override)
out = base;
if isempty(override)
  return;
end
fieldList = fieldnames(override);
for iField = 1:numel(fieldList)
  out.(fieldList{iField}) = override.(fieldList{iField});
end
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
fieldValue = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end


function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false), ', '), ']'];
end


function textOut = localFormatIntegerRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueVec, 'UniformOutput', false), ', '), ']'];
end