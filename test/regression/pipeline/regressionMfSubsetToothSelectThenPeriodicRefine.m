% Pipeline regression test for a two-stage anti-periodic strategy:
%   1) use a non-periodic 10-frame subset to pick the Doppler tooth;
%   2) refine on the original periodic 10-frame window with a narrow
%      [fdRef, fdRate] box centered at the selected tooth.
% This script does not change the main algorithm. It is only meant to test
% whether "non-periodic tooth selection + periodic full-data refinement"
% is a useful way to decouple the 1/T_f comb from the final coherent run.
clear(); close all;

localAddProjectPath();

periodicOffsetIdx = -4:5;
masterOffsetIdx = -9:10;
deterministicOffsetIdx = [-9, -7, -6, -5, -4, 0, 1, 3, 7, 8];
numRandomTrial = 24;
truthFdHalfWidthHz = 200;
truthFdRateHalfWidthHzPerSec = 500;
refineDoaHalfWidthDeg = [0.01; 0.01];
rng(253);

fixturePeriodic = localBuildPipelineFixture(periodicOffsetIdx);
toothStepHz = 1 / fixturePeriodic.frameIntvlSec;

baselineCase = localRunMsUnknownCase(fixturePeriodic, fixturePeriodic.fdRange, fixturePeriodic.fdRateRange, struct(), struct());
baselineSummary = localBuildCaseSummary(baselineCase, fixturePeriodic.truth, toothStepHz);

scheduleCell = cell(numRandomTrial + 1, 1);
labelCell = strings(numRandomTrial + 1, 1);
scheduleCell{1} = deterministicOffsetIdx;
labelCell(1) = "deterministic";
masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
parfor iTrial = 1:numRandomTrial
  pickUse = sort([0, masterNoZero(randperm(numel(masterNoZero), 9))]);
  scheduleCell{iTrial + 1} = pickUse;
  labelCell(iTrial + 1) = "random" + string(iTrial);
end

subsetSummaryCell = cell(numel(scheduleCell), 1);
subsetCaseCell = cell(numel(scheduleCell), 1);
parfor iTrial = 1:numel(scheduleCell)
  [subsetSummaryCell{iTrial}, subsetCaseCell{iTrial}] = localRunSubsetSchedule(scheduleCell{iTrial}, toothStepHz);
end

isResolvedList = cellfun(@(s) s.isResolved, subsetSummaryCell);
if ~all(isResolvedList)
  badIdx = find(~isResolvedList, 1, 'first');
  error('regressionMfSubsetToothSelectThenPeriodicRefine:UnresolvedSubset', ...
    'A subset run failed to resolve (index %d, label %s).', badIdx, labelCell(badIdx));
end

bestSubsetIdx = localSelectBestSubset(subsetSummaryCell);
bestSubsetSummary = subsetSummaryCell{bestSubsetIdx};
bestSubsetCase = subsetCaseCell{bestSubsetIdx};

fdRangeRefine = [bestSubsetCase.estResult.fdRefEst - truthFdHalfWidthHz, ...
  bestSubsetCase.estResult.fdRefEst + truthFdHalfWidthHz];
fdRateRangeRefine = [bestSubsetCase.estResult.fdRateEst - truthFdRateHalfWidthHzPerSec, ...
  bestSubsetCase.estResult.fdRateEst + truthFdRateHalfWidthHzPerSec];

dynOverride = struct();
dynOverride.initDoaParam = bestSubsetCase.estResult.doaParamEst(:);
dynOverride.initDoaHalfWidth = refineDoaHalfWidthDeg(:);

refineCase = localRunMsUnknownCase(fixturePeriodic, fdRangeRefine, fdRateRangeRefine, dynOverride, dynOverride);
refineSummary = localBuildCaseSummary(refineCase, fixturePeriodic.truth, toothStepHz);

if ~refineSummary.isResolved
  error('regressionMfSubsetToothSelectThenPeriodicRefine:PeriodicRefineUnresolved', ...
    'The periodic refinement stage must resolve.');
end
if abs(refineSummary.toothIdx) > abs(baselineSummary.toothIdx)
  error('regressionMfSubsetToothSelectThenPeriodicRefine:ToothSelectionWorse', ...
    'The selected-tooth periodic refinement moved farther away from the truth tooth.');
end

fprintf('Running regressionMfSubsetToothSelectThenPeriodicRefine ...\n');
fprintf('  tooth step (Hz)                    : %.6f\n', toothStepHz);
localPrintCaseSummary('periodic wide baseline', baselineSummary);
fprintf('  number of subset candidates        : %d\n', numel(scheduleCell));
fprintf('  selected subset label              : %s\n', labelCell(bestSubsetIdx));
fprintf('  selected subset offsets            : %s\n', localFormatIntegerRow(scheduleCell{bestSubsetIdx}));
localPrintCaseSummary('selected subset', bestSubsetSummary);
fprintf('  periodic refine fdRange            : %s\n', localFormatNumericRow(fdRangeRefine));
fprintf('  periodic refine fdRateRange        : %s\n', localFormatNumericRow(fdRateRangeRefine));
fprintf('  periodic refine init DoA (deg)     : %s\n', localFormatNumericRow(bestSubsetCase.estResult.doaParamEst(:).')); 
localPrintCaseSummary('periodic full-data refine', refineSummary);
fprintf('  angle improvement over baseline    : %.6f deg\n', baselineSummary.angleErrDeg - refineSummary.angleErrDeg);
fprintf('  fdRef improvement over baseline    : %.6f Hz\n', abs(baselineSummary.fdRefErrHz) - abs(refineSummary.fdRefErrHz));
fprintf('  fdRate improvement over baseline   : %.6f Hz/s\n', abs(baselineSummary.fdRateErrHzPerSec) - abs(refineSummary.fdRateErrHzPerSec));
fprintf('PASS: regressionMfSubsetToothSelectThenPeriodicRefine\n');


function [summary, caseUse] = localRunSubsetSchedule(offsetIdxVec, toothStepHz)
fixtureUse = localBuildPipelineFixture(offsetIdxVec);
caseUse = localRunMsUnknownCase(fixtureUse, fixtureUse.fdRange, fixtureUse.fdRateRange, struct(), struct());
summary = localBuildCaseSummary(caseUse, fixtureUse.truth, toothStepHz);
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
  error('regressionMfSubsetToothSelectThenPeriodicRefine:MissingReferenceFrame', ...
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

gridSize = [50, 50];
searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatTrueHz(:)], 0.1, 2e4);
fdRateRange = expandRangeToTruth([-1e4, 0], [truth.fdRateFit; truth.fdRateFit], 0.1, 5e2);
debugTruthMs = localBuildTruthDebugForView(viewMs, truth);

fixture = struct();
fixture.truth = truth;
fixture.sceneSeq = sceneSeq;
fixture.sceneRef = sceneRef;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
fixture.rxSigCell = rxSigCell;
fixture.rxSigRef = rxSigRef;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.wavelen = wavelen;
fixture.E = E;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.viewMs = viewMs;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
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
debugTruth.satIdxGlobal = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), [], 1);
debugTruth.refWeight = reshape(localGetFieldOrDefault(truth, 'refWeight', nan(sceneRef.numSat, 1)), [], 1);
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
numSat = sceneRef.numSat;
numFrame = size(truth.fdSatSeries, 2);
fdLocalTrue = nan(numSat, numFrame);
sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
truthSatIdx = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
if ~isempty(sceneSatIdx) && ~isempty(truthSatIdx)
  [isFound, loc] = ismember(sceneSatIdx, truthSatIdx);
  if any(isFound)
    fdLocalTrue(isFound, :) = truth.fdSatSeries(loc(isFound), :);
    return;
  end
end
numCopy = min(numSat, size(truth.fdSatSeries, 1));
fdLocalTrue(1:numCopy, :) = truth.fdSatSeries(1:numCopy, :);
end


function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
truthLocalDoa = [];
if isempty(sceneSeq) || ~isfield(sceneSeq, 'localDoa') || isempty(sceneSeq.localDoa)
  return;
end
localDoa = sceneSeq.localDoa;
if ndims(localDoa) == 3
  truthLocalDoa = localDoa;
elseif ndims(localDoa) == 4
  truthLocalDoa = localDoa(:, :, 1, :);
else
  error('regressionMfSubsetToothSelectThenPeriodicRefine:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end
end


function dataOut = localMergeStruct(baseIn, overrideIn)
dataOut = baseIn;
if isempty(overrideIn)
  return;
end
fieldList = fieldnames(overrideIn);
for iField = 1:numel(fieldList)
  dataOut.(fieldList{iField}) = overrideIn.(fieldList{iField});
end
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
try
  fieldValue = dataStruct.(fieldName);
catch
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


function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
