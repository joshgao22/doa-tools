% Diagnostic statistics for periodic vs non-periodic frame schedules.
% This script extends compareMfPeriodicVsRandomSubset by repeating the
% non-periodic subset experiment over multiple deterministic/random
% schedules while preserving the original frame timestamps.
% The goal is to measure how often strict 1/T_f comb attraction survives
% once the periodic sampling pattern is broken.
clear(); close all; clc;

localAddProjectPath();

periodicOffsetIdx = -4:5;
masterOffsetIdx = -9:10;
deterministicOffsetIdx = [-9, -7, -6, -5, -4, 0, 1, 3, 7, 8];
numRandomTrial = 24;
rng(253);

fixturePeriodic = localBuildPipelineFixture(periodicOffsetIdx);
periodicCase = localRunMsUnknownCase(fixturePeriodic, fixturePeriodic.fdRange, fixturePeriodic.fdRateRange);
toothStepHz = 1 / fixturePeriodic.frameIntvlSec;
periodicSummary = localBuildCaseSummary(periodicCase, fixturePeriodic.truth, toothStepHz);

summaryCell = cell(numRandomTrial + 1, 1);
scheduleCell = cell(numRandomTrial + 1, 1);
labelCell = strings(numRandomTrial + 1, 1);

scheduleCell{1} = deterministicOffsetIdx;
labelCell(1) = "deterministic";

masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
for iTrial = 1:numRandomTrial
  pickUse = sort([0, masterNoZero(randperm(numel(masterNoZero), 9))]);
  scheduleCell{iTrial + 1} = pickUse;
  labelCell(iTrial + 1) = "random" + string(iTrial);
end

parfor iTrial = 1:numel(scheduleCell)
  summaryCell{iTrial} = localRunOneSchedule(scheduleCell{iTrial}, toothStepHz);
end

isResolvedList = cellfun(@(s) s.isResolved, summaryCell);
if ~all(isResolvedList)
  badIdx = find(~isResolvedList, 1, 'first');
  error('scanMfPeriodicVsRandomSubsetStats:UnresolvedCase', ...
    'A non-periodic schedule failed to resolve (index %d, label %s).', ...
    badIdx, labelCell(badIdx));
end

angleErrList = cellfun(@(s) s.angleErrDeg, summaryCell);
fdRefErrList = cellfun(@(s) s.fdRefErrHz, summaryCell);
fdRateErrList = cellfun(@(s) s.fdRateErrHzPerSec, summaryCell);
toothIdxList = cellfun(@(s) s.toothIdx, summaryCell);
toothResidualList = cellfun(@(s) s.toothResidualHz, summaryCell);
absResidualList = abs(toothResidualList);
truthToothMask = toothIdxList == 0;

uniqueTooth = unique([periodicSummary.toothIdx; toothIdxList(:)]);
countText = strings(0, 1);
for iVal = 1:numel(uniqueTooth)
  idxUse = uniqueTooth(iVal);
  countUse = sum(toothIdxList == idxUse);
  countText(end + 1) = sprintf('%d:%d', idxUse, countUse); %#ok<AGROW>
end

bestIdx = find(absResidualList == min(absResidualList), 1, 'first');
worstIdx = find(absResidualList == max(absResidualList), 1, 'first');

fprintf('Running scanMfPeriodicVsRandomSubsetStats ...\n');
localPrintCaseSummary('periodic baseline', periodicSummary);
fprintf('  deterministic offsets              : %s\n', localFormatIntegerRow(deterministicOffsetIdx));
localPrintCaseSummary('deterministic subset', summaryCell{1});
fprintf('  number of random trials            : %d\n', numRandomTrial);
fprintf('  truth-tooth hit count              : %d / %d\n', sum(truthToothMask), numel(summaryCell));
fprintf('  tooth count map                    : %s\n', strjoin(cellstr(countText), ', '));
fprintf('  mean angle err (deg)               : %.6f\n', mean(angleErrList));
fprintf('  median angle err (deg)             : %.6f\n', median(angleErrList));
fprintf('  mean fdRef err (Hz)                : %.6f\n', mean(fdRefErrList));
fprintf('  median fdRef err (Hz)              : %.6f\n', median(fdRefErrList));
fprintf('  mean fdRate err (Hz/s)             : %.6f\n', mean(fdRateErrList));
fprintf('  median fdRate err (Hz/s)           : %.6f\n', median(fdRateErrList));
fprintf('  mean |tooth residual| (Hz)         : %.6f\n', mean(absResidualList));
fprintf('  median |tooth residual| (Hz)       : %.6f\n', median(absResidualList));
fprintf('  best subset label                  : %s\n', labelCell(bestIdx));
fprintf('  best subset offsets                : %s\n', localFormatIntegerRow(scheduleCell{bestIdx}));
localPrintCaseSummary('best subset', summaryCell{bestIdx});
fprintf('  worst subset label                 : %s\n', labelCell(worstIdx));
fprintf('  worst subset offsets               : %s\n', localFormatIntegerRow(scheduleCell{worstIdx}));
localPrintCaseSummary('worst subset', summaryCell{worstIdx});
if sum(truthToothMask) <= floor(numel(summaryCell) / 2)
  fprintf('  note                               : truth tooth is not dominant across random non-periodic schedules.\n');
end
fprintf('PASS: scanMfPeriodicVsRandomSubsetStats\n');


function summary = localRunOneSchedule(offsetIdxVec, toothStepHz)
fixtureUse = localBuildPipelineFixture(offsetIdxVec);
caseUse = localRunMsUnknownCase(fixtureUse, fixtureUse.fdRange, fixtureUse.fdRateRange);
summary = localBuildCaseSummary(caseUse, fixtureUse.truth, toothStepHz);
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


function caseDynMsUnknown = localRunMsUnknownCase(fixture, fdRangeUse, fdRateRangeUse)
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
caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fdRangeUse, fdRateRangeUse, false, dynMsKnownOpt, true, ...
  fixture.debugTruthMs, buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];

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
  error('scanMfPeriodicVsRandomSubsetStats:MissingReferenceFrame', ...
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
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
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
truth.pickAux = satPickAux;
truth.fdRefTrueHz = truth.fdRefSeries(refFrameIdx);
truth.fdRateTrueHzPerSec = truth.fdRateFit;
truth.fdSatTrueHz = reshape(truth.fdSatSeries(:, refFrameIdx), [], 1);
truth.deltaFdTrueHz = reshape(truth.deltaFdSeries(:, refFrameIdx), [], 1);

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

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
gridSize = [50, 50];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth([-1e4, 0], fdRateTruthCand, 0.1, 5e2);

debugTruthMs = localBuildTruthDebugForView(viewMs, truth);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.sceneRef = sceneRef;
fixture.rxSigCell = rxSigCell;
fixture.rxSigRef = rxSigRef;
fixture.viewMs = viewMs;
fixture.truth = truth;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
fixture.frameIntvlSec = frameIntvlSec;
fixture.wavelen = wavelen;
fixture.E = E;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
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
if isfield(sceneRef, 'satIdx') && ~isempty(sceneRef.satIdx)
  debugTruth.satIdxGlobal = reshape(sceneRef.satIdx, [], 1);
else
  debugTruth.satIdxGlobal = reshape(truth.selectedSatIdxGlobal(1:sceneRef.numSat), [], 1);
end
if isfield(sceneRef, 'ref') && isstruct(sceneRef.ref) && isfield(sceneRef.ref, 'weight') && ~isempty(sceneRef.ref.weight)
  debugTruth.refWeight = reshape(sceneRef.ref.weight, [], 1);
else
  debugTruth.refWeight = nan(sceneRef.numSat, 1);
end
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
numSat = sceneRef.numSat;
numFrame = size(truth.fdSatSeries, 2);
fdLocalTrue = nan(numSat, numFrame);
sceneSatIdx = [];
if isfield(sceneRef, 'satIdx') && ~isempty(sceneRef.satIdx)
  sceneSatIdx = reshape(sceneRef.satIdx, 1, []);
end
truthSatIdx = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
if ~isempty(sceneSatIdx) && ~isempty(truthSatIdx)
  [isFound, loc] = ismember(sceneSatIdx, truthSatIdx);
  if any(isFound)
    fdLocalTrue(isFound, :) = truth.fdSatSeries(loc(isFound), :);
  end
  return;
end
numCopy = min(numSat, size(truth.fdSatSeries, 1));
fdLocalTrue(1:numCopy, :) = truth.fdSatSeries(1:numCopy, :);
end


function localDoa = localExtractSceneLocalDoa(sceneSeq)
localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
elseif ndims(localDoaRaw) == 4
  localDoa = localDoaRaw(:, :, 1, :);
else
  error('scanMfPeriodicVsRandomSubsetStats:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    value = dataStruct.(fieldName);
  end
  return;
end
if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    value = dataStruct.(fieldName);
  end
  return;
end
try
  value = dataStruct.(fieldName);
catch
end
end


function textOut = localFormatIntegerRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueVec, 'UniformOutput', false), ', '), ']'];
end


function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
