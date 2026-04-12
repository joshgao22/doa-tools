% Development check for the formal MF route under the selected static pair
% This script follows the same case-building pattern as
% `doaDopplerStatDualSatUraEci` and then extends it to the dynamic
% multi-frame CP model:
%   1) Build SS-SF and MS-SF baselines on the reference frame.
%   2) Keep the static ablation and sat2-weight sweep as a transition
%      diagnostic before entering the dynamic model.
%   3) Use MS-SF-DoA as the coarse global DoA initializer.
%   4) Run SS-MF / MS-MF continuous-phase refinement with local DoA
%      boxes and compare them against the static baselines.
clear(); close all; clc;

%% Parameters
numUsr = 1;
numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = ceil(numFrame / 2);

sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 11.7e9;
wavelen = 299792458 / carrierFreq;
rng(253);

elemSpace = wavelen / 2;
numElem = [4 4];
snrDb = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);

utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa(numElem, elemSpace);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

fdRange = [-1e5, 0];
fdRateRange = [-1e4, 0];
optVerbose = false;
weightSweepAlpha = [0, 0.25, 0.5, 1];

if numUsr ~= 1
  error('doaDopplerDynDualSatUraEci:OnlySingleUserSupported', ...
    'This script currently supports one user only.');
end

%% Select Two Visible Satellites And The Reference Satellite
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

%% Multi-frame scene with satellite reference
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
refFrameIdx = sceneSeq.refFrameIdx;
sceneRef = sceneSeq.sceneCell{refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);

if sceneSeq.numSat ~= 2
  error('doaDopplerDynDualSatUraEci:InvalidNumSat', ...
    'This script expects exactly two selected satellites.');
end
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
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;
truth.fdRefTrueHz = truth.fdRefSeries(refFrameIdx);
truth.fdRateTrueHzPerSec = truth.fdRateFit;
truth.fdSatTrueHz = reshape(truth.fdSatSeries(:, refFrameIdx), [], 1);
truth.deltaFdTrueHz = reshape(truth.deltaFdSeries(:, refFrameIdx), [], 1);
truth.localDoaRef = reshape(sceneRef.localDoa(:, refSatIdxLocal), 2, 1);
fdRange = expandRangeToTruth(fdRange, [truth.fdRefFit; truth.fdSatTrueHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth(fdRateRange, fdRateTruthCand, 0.1, 5e2);

if truth.fdRefFit < fdRange(1) || truth.fdRefFit > fdRange(2)
  warning('doaDopplerDynDualSatUraEci:FdRefRangeMiss', ...
    'truth.fdRefFit = %.3f Hz is outside fdRange = [%.3f, %.3f] Hz.', ...
    truth.fdRefFit, fdRange(1), fdRange(2));
end
if truth.fdRateFit < fdRateRange(1) || truth.fdRateFit > fdRateRange(2)
  warning('doaDopplerDynDualSatUraEci:FdRateRangeMiss', ...
    'truth.fdRateFit = %.3f Hz/s is outside fdRateRange = [%.3f, %.3f] Hz/s.', ...
    truth.fdRateFit, fdRateRange(1), fdRateRange(2));
end

steerDiag = buildSteeringDriftDiag(sceneSeq, refFrameIdx, truth.selectedSatIdxGlobal);
fdDiag = buildDeltaFdFitDiag(truth);

%% Pilot waveform
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

%% Multi-frame snapshots
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

%% Single-satellite and multi-satellite views
sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = selectSatScene(sceneRef, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSigRef, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  gridSize, searchRange, E);

viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

%% Common estimator options
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
dynBaseOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = true;
dynBaseOpt.debugStoreEvalTrace = false;
dynBaseOpt.debugMaxEvalTrace = 120;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha, [0.01; 0.01]);

caseRefDoa = caseBundle.caseRefDoa;
caseMsDoa = caseBundle.caseMsDoa;
caseOtherDoa = caseBundle.caseOtherDoa;
staticRefOpt = caseBundle.staticRefOpt;
staticOtherOpt = caseBundle.staticOtherOpt;
staticMsOpt = caseBundle.staticMsOpt;
caseStaticRefOnly = caseBundle.caseStaticRefOnly;
caseStaticRefAbl = caseBundle.caseStaticRefAbl;
caseStaticOtherOnly = caseBundle.caseStaticOtherOnly;
caseStaticMs = caseBundle.caseStaticMs;
weightCase = caseBundle.weightCase;
bestStaticMsCase = caseBundle.bestStaticMsCase;

initParamStaticRef = buildDynamicInitParamFromCase(caseStaticRefOnly, true, truth.fdRateFit);
initParamStaticMs = buildDynamicInitParamFromCase(bestStaticMsCase, true, truth.fdRateFit);
%% Dynamic CP cases

% Use the single-frame static estimates as the dynamic coarse state. This
% mirrors the intended pipeline more closely than restarting MF directly
% from the DoA-only coarse point. For CP-U, continue from the CP-K result
% so that fdRate is released from a physically meaningful basin instead of
% restarting near zero.

dynRefKnownOpt = dynBaseOpt;
dynRefKnownOpt.initDoaParam = caseStaticRefOnly.estResult.doaParamEst(:);
dynRefKnownOpt.initDoaHalfWidth = [0.005; 0.005];

dynRefUnknownOpt = dynBaseOpt;
dynRefUnknownOpt.initDoaParam = caseStaticRefOnly.estResult.doaParamEst(:);
dynRefUnknownOpt.initDoaHalfWidth = [0.003; 0.003];

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = [0.003; 0.003];
dynMsKnownOpt.enableFdAliasUnwrap = true;

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];
dynMsUnknownOpt.enableFdAliasUnwrap = true;

debugTruthRef = localBuildTruthDebugForView(viewRefOnly, truth);
debugTruthMs = localBuildTruthDebugForView(viewMs, truth);

caseDynRefKnown = runDynamicDoaDopplerCase("SS-MF-CP-K", "single", ...
  viewRefOnly, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynRefKnownOpt, true, debugTruthRef, initParamStaticRef);

initParamRefUnknownCpK = buildDynamicInitParamFromCase(caseDynRefKnown, false, caseDynRefKnown.estResult.fdRateEst);
initParamRefUnknownStatic = buildDynamicInitParamFromCase(caseStaticRefOnly, false, truth.fdRateFit);
refUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynRefKnown, caseStaticRefOnly, initParamRefUnknownCpK, initParamRefUnknownStatic, ...
  dynRefUnknownOpt.initDoaHalfWidth, dynRefKnownOpt.initDoaHalfWidth);
caseDynRefUnknown = runDynamicDoaDopplerCase("SS-MF-CP-U", "single", ...
  viewRefOnly, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynRefUnknownOpt, false, debugTruthRef, refUnknownCand);

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsKnownOpt, true, debugTruthMs, initParamStaticMs);

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, truth.fdRateFit);
msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsUnknownOpt, false, debugTruthMs, msUnknownCand);

baseCase = [ ...
  caseRefDoa, ...
  caseMsDoa, ...
  caseStaticRefOnly, ...
  caseStaticMs, ...
  caseDynRefKnown, ...
  caseDynRefUnknown, ...
  caseDynMsKnown, ...
  caseDynMsUnknown];
caseResult = [baseCase, weightCase];

dynSummary = summarizeDynamicEstimatorDebug(baseCase, sceneSeq, truth);
dynDiagTable = dynSummary.diagTable;
dynObjTable = dynSummary.objectiveTable;
dynSatTable = dynSummary.perSatTable;
dynBlockTable = dynSummary.blockTable;
dynLocalStateTable = dynSummary.localStateTable;
dynTruthStateTable = dynSummary.truthStateTable;
dynSatOrderTable = dynSummary.satOrderTable;
dynMultiStartTable = summarizeDynamicMultiStart(baseCase, truth);

%% CRB calculation
crbSfOpt = struct();
crbSfOpt.doaType = 'latlon';
[crbSfRef, auxCrbSfRef] = localTryBuildStaticCrb( ...
  sceneRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbSfOpt);
[crbSfMs, auxCrbSfMs] = localTryBuildStaticCrb( ...
  sceneRef, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbSfOpt);

pathGainRef = ones(sceneSeqRefOnly.numSat, sceneSeqRefOnly.numFrame);
noiseVarRef = pwrNoise * ones(sceneSeqRefOnly.numSat, sceneSeqRefOnly.numFrame);
pathGainMs = ones(sceneSeq.numSat, sceneSeq.numFrame);
noiseVarMs = pwrNoise * ones(sceneSeq.numSat, sceneSeq.numFrame);

crbMfKnownOpt = struct();
crbMfKnownOpt.doaType = 'latlon';
crbMfKnownOpt.phaseMode = 'continuous';
crbMfKnownOpt.fdRateMode = 'known';
crbMfKnownOpt.steeringMode = dynBaseOpt.steeringMode;
crbMfKnownOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
[crbMfRefKnown, auxCrbMfRefKnown] = localTryBuildDynamicCrb( ...
  sceneSeqRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, ...
  pathGainRef, noiseVarRef, crbMfKnownOpt);
[crbMfMsKnown, auxCrbMfMsKnown] = localTryBuildDynamicCrb( ...
  sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, ...
  pathGainMs, noiseVarMs, crbMfKnownOpt);

crbMfUnknownOpt = crbMfKnownOpt;
crbMfUnknownOpt.fdRateMode = 'unknown';
[crbMfRefUnknown, auxCrbMfRefUnknown] = localTryBuildDynamicCrb( ...
  sceneSeqRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, ...
  pathGainRef, noiseVarRef, crbMfUnknownOpt);
[crbMfMsUnknown, auxCrbMfMsUnknown] = localTryBuildDynamicCrb( ...
  sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, ...
  pathGainMs, noiseVarMs, crbMfUnknownOpt);

crbSummary = localBuildDynamicCrbSummary(truth, ...
  crbSfRef, auxCrbSfRef, crbSfMs, auxCrbSfMs, ...
  crbMfRefKnown, auxCrbMfRefKnown, crbMfMsKnown, auxCrbMfMsKnown, ...
  crbMfRefUnknown, auxCrbMfRefUnknown, crbMfMsUnknown, auxCrbMfMsUnknown);

%% Summaries
estTable = buildDoaDopplerSummaryTable(caseResult, truth, struct('mode', 'dynamic'));

ablationTruth = [ ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneRefOnly), ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneOtherOnly)];
ablationTable = buildDoaDopplerCaseSummaryTable([caseStaticRefAbl, caseStaticOtherOnly], ablationTruth);
weightTable = buildDoaDopplerWeightSweepTable(weightSweepAlpha, weightCase, truth);

fprintf('\n========== Selected satellites ==========%s', newline);
disp(table((1:sceneSeq.numSat).', truth.selectedSatIdxGlobal(:), truth.usrElevationDeg(:), ...
  truth.fdSatTrueHz(:), truth.deltaFdTrueHz(:), ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg', 'fdSatRefFrameHz', 'deltaFdRefFrameHz'}));

fprintf('\n========== Truth ==========%s', newline);
disp(table(truth.latlonTrueDeg(1), truth.latlonTrueDeg(2), truth.fdRefFit, truth.fdRateFit, ...
  truth.refSatIdxLocal, truth.refSatIdxGlobal, otherSatIdxLocal, otherSatIdxGlobal, ...
  'VariableNames', {'latTrueDeg', 'lonTrueDeg', 'fdRefTrueHz', 'fdRateTrueHzPerSec', ...
  'refSatIdxLocal', 'refSatIdxGlobal', 'otherSatIdxLocal', 'otherSatIdxGlobal'}));

fprintf('\n========== Steering drift summary ==========%s', newline);
disp(steerDiag);

fprintf('\n========== Reference-Doppler linear-fit summary ==========%s', newline);
disp(fdDiag.refTable);

fprintf('\n========== Differential-Doppler linear-fit summary ==========%s', newline);
disp(fdDiag.deltaTable);

fprintf('\n========== Estimator summary ==========%s', newline);
disp(estTable);

fprintf('\n========== Static ablation summary ==========%s', newline);
disp(ablationTable);

fprintf('\n========== Static weight-sweep summary ==========%s', newline);
disp(weightTable);

fprintf('\n========== Dynamic estimator diagnostics ==========%s', newline);
disp(dynDiagTable);

fprintf('\n========== Dynamic objective summary ==========%s', newline);
disp(dynObjTable);

fprintf('\n========== Dynamic per-satellite summary ==========%s', newline);
disp(dynSatTable);

fprintf('\n========== Dynamic block summary ==========%s', newline);
disp(dynBlockTable);

fprintf('\n========== Dynamic local-state compare ==========%s', newline);
disp(dynLocalStateTable);

fprintf('\n========== Dynamic truth-state compare ==========%s', newline);
disp(dynTruthStateTable);

fprintf('\n========== Dynamic satellite-order summary ==========%s', newline);
disp(dynSatOrderTable);

if ~isempty(dynMultiStartTable)
  fprintf('\n========== Dynamic multi-start summary ==========%s', newline);
  disp(dynMultiStartTable);
end

fprintf('\n========== Dynamic CRB summary ==========%s', newline);
disp(crbSummary);

badMask = ~estTable.isResolved;
if any(badMask)
  warning('doaDopplerDynDualSatUraEci:EstimatorUnresolved', ...
    'Some estimators are not resolved:\n%s', evalc('disp(estTable(badMask, :))'));
end

%% Plots
plotDoaDopplerComparison(sceneSeq.timeOffsetSec, truth, baseCase);
plotDoaDopplerGeometryComparison(truth.latlonTrueDeg, baseCase);
localPlotWeightSweep(weightTable);

%% Local functions
function [crb, aux] = localTryBuildStaticCrb(scene, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRefHz, numSource, noiseVar, crbOpt)
%LOCALTRYBUILDSTATICCRB Build one static CRB with graceful fallback.

try
  [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, numSource, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = "zero";
  aux.phaseMode = "single-frame";
  aux.errorMessage = string(ME.message);
end
end

function [crb, aux] = localTryBuildDynamicCrb(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt)
%LOCALTRYBUILDDYNAMICCRB Build one dynamic CRB with graceful fallback.

try
  [crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = localGetFieldOrDefault(crbOpt, 'fdRateMode', "unknown");
  aux.phaseMode = localGetFieldOrDefault(crbOpt, 'phaseMode', "continuous");
  aux.errorMessage = string(ME.message);
end
end

function crbTable = localBuildDynamicCrbSummary(truth, ...
  crbSfRef, auxCrbSfRef, crbSfMs, auxCrbSfMs, ...
  crbMfRefKnown, auxCrbMfRefKnown, crbMfMsKnown, auxCrbMfMsKnown, ...
  crbMfRefUnknown, auxCrbMfRefUnknown, crbMfMsUnknown, auxCrbMfMsUnknown)
%LOCALBUILDDYNAMICCRBSUMMARY Build one compact CRB summary table.

caseName = ["SS-SF-Static"; "MS-SF-Static"; "SS-MF-CP-K"; ...
  "MS-MF-CP-K"; "SS-MF-CP-U"; "MS-MF-CP-U"];
satMode = ["single"; "multi"; "single"; "multi"; "single"; "multi"];
frameMode = ["single"; "single"; "multi"; "multi"; "multi"; "multi"];
angleStdDeg = [ ...
  projectCrbToAngleMetric(crbSfRef(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbSfMs(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMfRefKnown(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMfMsKnown(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMfRefUnknown(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMfMsUnknown(1:2, 1:2), truth.latlonTrueDeg, 'latlon')];
fdStdHz = [sqrt(max(real(crbSfRef(3, 3)), 0)); ...
  sqrt(max(real(crbSfMs(3, 3)), 0)); ...
  sqrt(max(real(crbMfRefKnown(3, 3)), 0)); ...
  sqrt(max(real(crbMfMsKnown(3, 3)), 0)); ...
  sqrt(max(real(crbMfRefUnknown(3, 3)), 0)); ...
  sqrt(max(real(crbMfMsUnknown(3, 3)), 0))];
phaseMode = string({auxCrbSfRef.phaseMode; auxCrbSfMs.phaseMode; ...
  auxCrbMfRefKnown.phaseMode; auxCrbMfMsKnown.phaseMode; ...
  auxCrbMfRefUnknown.phaseMode; auxCrbMfMsUnknown.phaseMode});
fdRateMode = string({auxCrbSfRef.fdRateMode; auxCrbSfMs.fdRateMode; ...
  auxCrbMfRefKnown.fdRateMode; auxCrbMfMsKnown.fdRateMode; ...
  auxCrbMfRefUnknown.fdRateMode; auxCrbMfMsUnknown.fdRateMode});

crbTable = table(caseName, satMode, frameMode, phaseMode, fdRateMode, ...
  angleStdDeg, fdStdHz, ...
  'VariableNames', {'displayName', 'satMode', 'frameMode', 'phaseMode', ...
  'fdRateMode', 'angleCrbStdDeg', 'fdRefCrbStdHz'});
end

function localPlotWeightSweep(weightTable)
%LOCALPLOTWEIGHTSWEEP Plot one compact sat2-weight sweep summary.

figure();
plot(weightTable.alphaSat2, weightTable.angleErrDeg, 'o-', 'LineWidth', 1.3, ...
  'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Angle error (deg)');
title('MS-SF-Static angle error under sat2 weighting');
end

function debugTruth = localBuildTruthDebugForView(view, truth)
%LOCALBUILDTRUTHDEBUGFORVIEW Build compact truth probes for one scene view.

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
elseif isfield(truth, 'selectedSatIdxGlobal') && ~isempty(truth.selectedSatIdxGlobal)
  debugTruth.satIdxGlobal = reshape(truth.selectedSatIdxGlobal(1:min(sceneRef.numSat, numel(truth.selectedSatIdxGlobal))), [], 1);
else
  debugTruth.satIdxGlobal = nan(sceneRef.numSat, 1);
end
if isfield(sceneRef, 'ref') && isstruct(sceneRef.ref) && isfield(sceneRef.ref, 'weight') && ~isempty(sceneRef.ref.weight)
  debugTruth.refWeight = reshape(sceneRef.ref.weight, [], 1);
else
  debugTruth.refWeight = nan(sceneRef.numSat, 1);
end
end

function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
%LOCALBUILDTRUTHFDLOCALFORVIEW Build view-aligned true local Doppler series.

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

function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
%LOCALEXTRACTSCENELOCALDOA Extract 2xNsxNf truth local DoA array.

truthLocalDoa = [];
if isempty(sceneSeq) || ~isfield(sceneSeq, 'localDoa') || isempty(sceneSeq.localDoa)
  return;
end

localDoa = sceneSeq.localDoa;
if ndims(localDoa) == 3
  truthLocalDoa = localDoa;
  return;
end

if ndims(localDoa) == 4
  truthLocalDoa = localDoa(:, :, 1, :);
  return;
end

error('doaDopplerDynDualSatUraEci:InvalidLocalDoaSize', ...
  'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

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

