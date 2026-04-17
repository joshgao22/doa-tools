% Diagnostic compare for periodic vs non-periodic frame schedules.
% This script keeps the current MS-MF-CP pipeline but changes only the
% multi-frame time schedule:
%   1) one periodic 10-frame window with uniform 1/750 spacing;
%   2) one deterministic random 10-frame subset taken from a 20-frame grid,
%      while preserving the original frame timestamps.
% The goal is to check whether breaking strict periodic sampling weakens the
% 1/T_f branch attraction.
clear(); close all; clc;

localAddProjectPath();

periodicOffsetIdx = -4:5;
masterOffsetIdx = -9:10;
rng(253);
masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
randPick = sort([0, masterNoZero(randperm(numel(masterNoZero), 9))]);
nonPeriodicOffsetIdx = sort(randPick);

fixturePeriodic = localBuildPipelineFixture(periodicOffsetIdx);
fixtureNonPeriodic = localBuildPipelineFixture(nonPeriodicOffsetIdx);

toothStepHz = 1 / fixturePeriodic.frameIntvlSec;

fprintf('Running compareMfPeriodicVsRandomSubset ...\n');

periodicCase = localRunMsUnknownCase(fixturePeriodic, fixturePeriodic.fdRange, fixturePeriodic.fdRateRange);
nonPeriodicCase = localRunMsUnknownCase(fixtureNonPeriodic, fixtureNonPeriodic.fdRange, fixtureNonPeriodic.fdRateRange);

periodicSummary = localBuildCaseSummary(periodicCase, fixturePeriodic.truth, toothStepHz);
nonPeriodicSummary = localBuildCaseSummary(nonPeriodicCase, fixtureNonPeriodic.truth, toothStepHz);

periodicGap = diff(periodicOffsetIdx);
nonPeriodicGap = diff(nonPeriodicOffsetIdx);

if ~all(abs(periodicGap - periodicGap(1)) < 1e-12)
  error('compareMfPeriodicVsRandomSubset:PeriodicScheduleInvalid', ...
    'The periodic baseline schedule must stay uniformly spaced.');
end
if all(abs(nonPeriodicGap - nonPeriodicGap(1)) < 1e-12)
  error('compareMfPeriodicVsRandomSubset:NonPeriodicScheduleInvalid', ...
    'The random subset schedule must not collapse back to uniform spacing.');
end
if ~periodicSummary.isResolved || ~nonPeriodicSummary.isResolved
  error('compareMfPeriodicVsRandomSubset:UnresolvedCase', ...
    'Both periodic and non-periodic runs must resolve in this regression script.');
end

fprintf('  frame interval (s)                 : %.9f\n', fixturePeriodic.frameIntvlSec);
fprintf('  tooth step (Hz)                    : %.6f\n', toothStepHz);
fprintf('  periodic offsets                   : %s\n', localFormatIntegerRow(periodicOffsetIdx));
fprintf('  periodic gaps                      : %s\n', localFormatIntegerRow(periodicGap));
fprintf('  non-periodic offsets               : %s\n', localFormatIntegerRow(nonPeriodicOffsetIdx));
fprintf('  non-periodic gaps                  : %s\n', localFormatIntegerRow(nonPeriodicGap));
localPrintCaseSummary('periodic', periodicSummary);
localPrintCaseSummary('non-periodic subset', nonPeriodicSummary);
fprintf('  angle improvement (deg)            : %.6f\n', periodicSummary.angleErrDeg - nonPeriodicSummary.angleErrDeg);
fprintf('  |tooth residual| improvement (Hz)  : %.6f\n', abs(periodicSummary.toothResidualHz) - abs(nonPeriodicSummary.toothResidualHz));
fprintf('PASS: compareMfPeriodicVsRandomSubset\n');


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
  error('compareMfPeriodicVsRandomSubset:MissingReferenceFrame', ...
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
fixture.wavelen = wavelen;
fixture.E = E;
fixture.gridSize = gridSize;
fixture.searchRange = searchRange;
fixture.refSatIdxLocal = refSatIdxLocal;
fixture.otherSatIdxLocal = otherSatIdxLocal;
fixture.otherSatIdxGlobal = otherSatIdxGlobal;
fixture.debugTruthMs = debugTruthMs;
fixture.frameIntvlSec = frameIntvlSec;
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
sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
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
  error('compareMfPeriodicVsRandomSubset:InvalidLocalDoaSize', ...
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
