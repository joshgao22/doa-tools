% Development script for diagnosing frame-interval induced fdRef branches.
% The script scans the MF dynamic objective along one dimension:
%   fdRef = fdRefCenter + deltaFd
% while keeping DoA and fdRate fixed at selected centers. The goal is to
% check whether the objective exhibits a comb-like near-periodicity with
% spacing close to 1 / frameIntvlSec.
%
% In addition to the baseline uniform-time evaluator, the script can build a
% probe-only jittered evaluator by perturbing model.timeOffsetSec. This is a
% deliberate diagnostic mismatch used to test whether the comb structure is
% tied to a perfectly uniform frame schedule.
clear(); close all; clc;

%% Parameters
numUsr = 1;
numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = ceil(numFrame / 2);

sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
numSym = 32;
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
weightSweepAlpha = [0; 0.25; 0.5; 1];

scanCfg = struct();
scanCfg.runKnown = true;
scanCfg.runUnknown = true;
scanCfg.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanCfg.centerListUnknown = ["truth"; "cpKnownSeed"; "finalEstimate"];
scanCfg.numAliasSide = 2;
scanCfg.numGrid = 801;
scanCfg.scanHalfWidthHz = [];
scanCfg.jitterEnable = true;
scanCfg.jitterStdSec = 2e-6;
scanCfg.jitterMode = "zeroMeanRandom";
scanCfg.showDeltaFigure = true;
scanCfg.showPeakFigure = true;
scanCfg.showFoldedFigure = true;
scanCfg.saveTag = "doaDopplerDynFdRefCombScan";
scanCfg.useParfor = true;
scanCfg.autoStartParpool = true;


%% Select satellites and build scene
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
refFrameIdx = sceneSeq.refFrameIdx;
sceneRef = sceneSeq.sceneCell{refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);

if sceneSeq.numSat ~= 2
  error('doaDopplerDynFdRefCombScan:InvalidNumSat', ...
    'This script expects exactly two selected satellites.');
end

otherSatIdxLocal = setdiff(1:sceneSeq.numSat, refSatIdxLocal);
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

fdRange = expandRangeToTruth(fdRange, [truth.fdRefFit; truth.fdSatTrueHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; ...
  truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth(fdRateRange, fdRateTruthCand, 0.1, 5e2);

%% Pilot waveform and snapshots
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

%% Views
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

%% Static transition and dynamic centers
staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha, [0.01; 0.01]);

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;
initParamStaticMs = buildDynamicInitParamFromCase(bestStaticMsCase, true, truth.fdRateFit);

dynBaseOpt = struct();
dynBaseOpt.useLogObjective = true;
dynBaseOpt.initFdCount = 81;
dynBaseOpt.useAccessMask = false;
dynBaseOpt.phaseMode = 'continuous';
dynBaseOpt.steeringMode = 'framewise';
dynBaseOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = false;
dynBaseOpt.debugStoreEvalTrace = false;
dynBaseOpt.debugMaxEvalTrace = 1;

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = [0.003; 0.003];
dynMsKnownOpt.enableFdAliasUnwrap = true;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsKnownOpt, true, [], initParamStaticMs);

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];
dynMsUnknownOpt.enableFdAliasUnwrap = true;

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsUnknownOpt, false, [], msUnknownCand);

centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown);

localMaybeStartParpool(scanCfg);

%% Build MF debug models
baseDynOpt = struct();
baseDynOpt.useLogObjective = true;
baseDynOpt.initFdCount = 81;
baseDynOpt.useAccessMask = false;
baseDynOpt.phaseMode = 'continuous';
baseDynOpt.steeringMode = 'framewise';
baseDynOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
baseDynOpt.enableFdAliasUnwrap = true;
baseDynOpt.debugEnable = false;
baseDynOpt.debugStoreEvalTrace = false;
baseDynOpt.debugMaxEvalTrace = 1;
baseDynOpt.initDoaParam = [];
baseDynOpt.initDoaHalfWidth = [];

scanResult = struct();
scanResult.truth = truth;
scanResult.scanCfg = scanCfg;
scanResult.centerRegistry = centerRegistry;
scanResult.frameInfo = struct( ...
  'frameIntvlSec', frameIntvlSec, ...
  'timeOffsetSec', sceneSeq.timeOffsetSec(:), ...
  'numFrame', sceneSeq.numFrame, ...
  'numPilotSample', numel(pilotWave), ...
  'refFrameIdx', sceneSeq.refFrameIdx, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'refStateSource', string(refState.source));

if scanCfg.runKnown
  modelOptKnown = baseDynOpt;
  modelOptKnown.fdRateMode = 'known';
  modelOptKnown.fdRateKnown = truth.fdRateFit;
  [modelKnown, ~, ~, modelOptKnown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, [], modelOptKnown);
  scanResult.cpKnown = localRunFdRefCombSet(modelKnown, centerRegistry, scanCfg, 'known');
  scanResult.cpKnown.modelOpt = modelOptKnown;
end

if scanCfg.runUnknown
  modelOptUnknown = baseDynOpt;
  modelOptUnknown.fdRateMode = 'unknown';
  [modelUnknown, ~, ~, modelOptUnknown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, fdRateRange, modelOptUnknown);
  scanResult.cpUnknown = localRunFdRefCombSet(modelUnknown, centerRegistry, scanCfg, 'unknown');
  scanResult.cpUnknown.modelOpt = modelOptUnknown;
end

%% Report and plot
if isfield(scanResult, 'cpKnown')
  disp('========== CP-K fdRef comb summary ==========');
  localDisplayModeSummary('CP-K', scanResult.cpKnown);
end
if isfield(scanResult, 'cpUnknown')
  disp('========== CP-U fdRef comb summary ==========');
  localDisplayModeSummary('CP-U', scanResult.cpUnknown);
end

if localGetFieldOrDefault(scanCfg, 'showDeltaFigure', true)
  if isfield(scanResult, 'cpKnown')
    localPlotModeDeltaFigure(scanResult.cpKnown, 'CP-K');
  end
  if isfield(scanResult, 'cpUnknown')
    localPlotModeDeltaFigure(scanResult.cpUnknown, 'CP-U');
  end
end
if localGetFieldOrDefault(scanCfg, 'showPeakFigure', true)
  if isfield(scanResult, 'cpKnown')
    localPlotModePeakFigure(scanResult.cpKnown, 'CP-K');
  end
  if isfield(scanResult, 'cpUnknown')
    localPlotModePeakFigure(scanResult.cpUnknown, 'CP-U');
  end
end
if localGetFieldOrDefault(scanCfg, 'showFoldedFigure', true)
  if isfield(scanResult, 'cpKnown')
    localPlotModeFoldedFigure(scanResult.cpKnown, 'CP-K');
  end
  if isfield(scanResult, 'cpUnknown')
    localPlotModeFoldedFigure(scanResult.cpUnknown, 'CP-U');
  end
end

scanResult.utcRun = datetime('now', 'TimeZone', 'local');
scanResult.scriptTag = scanCfg.saveTag;
saveFile = fullfile(fileparts(mfilename('fullpath')), ...
  sprintf('%s_%s.mat', scanCfg.saveTag, datestr(now, 'yyyymmdd-HHMMSS')));
save(saveFile, 'scanResult', '-v7.3');
fprintf('Saved scanResult to %s\n', saveFile);


function modeResult = localRunFdRefCombSet(modelBase, centerRegistry, scanCfg, modeTag)
%LOCALRUNFDREFCOMBSET Run one-dimensional fdRef scans for one mode.

centerNameList = localResolveCenterNameList(scanCfg, modeTag);
modeResult = struct();
modeResult.modeTag = string(modeTag);
modeResult.centerNameList = centerNameList;
modeResult.aliasStepHz = localResolveAliasStep(modelBase);
modeResult.centerScan = repmat(struct( ...
  'centerName', "", ...
  'basePoint', struct(), ...
  'uniform', struct(), ...
  'jittered', struct()), numel(centerNameList), 1);

for iCenter = 1:numel(centerNameList)
  centerName = centerNameList(iCenter);
  basePoint = localResolveCenterPoint(centerRegistry, modeTag, centerName);
  halfWidthHz = localResolveHalfWidth(scanCfg, modeResult.aliasStepHz);
  fdRefGrid = linspace(basePoint.fdRef - halfWidthHz, basePoint.fdRef + halfWidthHz, scanCfg.numGrid);

  centerScan = struct();
  centerScan.centerName = centerName;
  centerScan.basePoint = basePoint;
  centerScan.uniform = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanCfg, false);

  if localGetFieldOrDefault(scanCfg, 'jitterEnable', false)
    centerScan.jittered = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanCfg, true);
  end

  modeResult.centerScan(iCenter) = centerScan;
end
end

function lineInfo = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanCfg, useJitter)
%LOCALSCANFDREFLINE Evaluate one fdRef comb scan.

modelScan = modelBase;
lineTag = "uniform";
if useJitter
  lineTag = "jittered";
  modelScan = localBuildJitteredModel(modelBase, scanCfg);
end

numGrid = numel(fdRefGrid);
objVec = nan(numGrid, 1);
residualVec = nan(numGrid, 1);
fdRefEvalVec = nan(numGrid, 1);
fdRateEvalVec = nan(numGrid, 1);
cohSatMat = nan(numGrid, modelScan.numSat);
objSatMat = nan(numGrid, modelScan.numSat);

probeTemplate = struct('doaParam', basePoint.latlon(:), ...
  'fdRef', basePoint.fdRef, 'fdRate', basePoint.fdRate, ...
  'obj', NaN, 'residualNorm', NaN);

useParfor = localResolveUseParfor(scanCfg, numGrid);
if useParfor
  parfor iGrid = 1:numGrid
    [objCur, residualCur, fdRefCur, fdRateCur, cohCur, objSatCur] = ...
      localEvaluateScanPoint(modelScan, basePoint, fdRefGrid(iGrid), probeTemplate);
    objVec(iGrid) = objCur;
    residualVec(iGrid) = residualCur;
    fdRefEvalVec(iGrid) = fdRefCur;
    fdRateEvalVec(iGrid) = fdRateCur;
    cohSatMat(iGrid, :) = cohCur;
    objSatMat(iGrid, :) = objSatCur;
  end
else
  for iGrid = 1:numGrid
    [objCur, residualCur, fdRefCur, fdRateCur, cohCur, objSatCur] = ...
      localEvaluateScanPoint(modelScan, basePoint, fdRefGrid(iGrid), probeTemplate);
    objVec(iGrid) = objCur;
    residualVec(iGrid) = residualCur;
    fdRefEvalVec(iGrid) = fdRefCur;
    fdRateEvalVec(iGrid) = fdRateCur;
    cohSatMat(iGrid, :) = cohCur;
    objSatMat(iGrid, :) = objSatCur;
  end
end

[minObj, minIdx] = min(objVec);
aliasStepHz = localResolveAliasStep(modelScan);
deltaFdRef = fdRefGrid(:) - basePoint.fdRef;
foldedOffset = localFoldOffset(deltaFdRef, aliasStepHz);
reciprocalPeak = localBuildReciprocalPeak(objVec, minObj);

lineInfo = struct();
lineInfo.lineTag = lineTag;
lineInfo.useJitter = useJitter;
lineInfo.aliasStepHz = aliasStepHz;
lineInfo.basePoint = basePoint;
lineInfo.fdRefGrid = fdRefGrid(:);
lineInfo.deltaFdRef = deltaFdRef;
lineInfo.foldedOffset = foldedOffset;
lineInfo.objVec = objVec;
lineInfo.deltaObjVec = objVec - minObj;
lineInfo.reciprocalPeak = reciprocalPeak;
lineInfo.residualVec = residualVec;
lineInfo.fdRefEvalVec = fdRefEvalVec;
lineInfo.fdRateEvalVec = fdRateEvalVec;
lineInfo.cohSatMat = cohSatMat;
lineInfo.objectiveSatMat = objSatMat;
lineInfo.minObj = minObj;
lineInfo.minIdx = minIdx;
lineInfo.minFdRef = fdRefGrid(minIdx);
lineInfo.minDeltaFdRef = fdRefGrid(minIdx) - basePoint.fdRef;
lineInfo.minAliasIndex = round(lineInfo.minDeltaFdRef / aliasStepHz);
lineInfo.centerIdx = localFindNearestGridIndex(fdRefGrid, basePoint.fdRef);
lineInfo.centerObj = objVec(lineInfo.centerIdx);
lineInfo.centerDeltaObj = lineInfo.centerObj - minObj;
lineInfo.jitterInfo = localBuildJitterInfo(modelScan, modelBase, scanCfg, useJitter);
end

function [objCur, residualCur, fdRefCur, fdRateCur, cohCur, objSatCur] = ...
    localEvaluateScanPoint(modelScan, basePoint, fdRefCurIn, probeTemplate)
%LOCALEVALUATESCANPOINT Evaluate one fdRef grid point.

pointCur = basePoint;
pointCur.fdRef = fdRefCurIn;
probeEval = localEvalProbeAtPoint(modelScan, pointCur, probeTemplate);
objCur = probeEval.obj;
residualCur = localGetFieldOrDefault(probeEval, 'residualNorm', NaN);
fdRefCur = localGetFieldOrDefault(probeEval, 'fdRef', NaN);
fdRateCur = localGetFieldOrDefault(probeEval, 'fdRate', NaN);

cohCur = nan(1, modelScan.numSat);
objSatCur = nan(1, modelScan.numSat);
cohSat = reshape(localGetFieldOrDefault(probeEval, 'coherenceSat', []), 1, []);
objSat = reshape(localGetFieldOrDefault(probeEval, 'objectiveSat', []), 1, []);
numSatCopy = min(modelScan.numSat, numel(cohSat));
if numSatCopy > 0
  cohCur(1, 1:numSatCopy) = cohSat(1, 1:numSatCopy);
end
numSatCopy = min(modelScan.numSat, numel(objSat));
if numSatCopy > 0
  objSatCur(1, 1:numSatCopy) = objSat(1, 1:numSatCopy);
end
end

function jitterInfo = localBuildJitterInfo(modelJit, modelBase, scanCfg, useJitter)
%LOCALBUILDJITTERINFO Build compact time-offset diagnostics.

jitterInfo = struct();
jitterInfo.isEnabled = useJitter;
jitterInfo.timeOffsetSecBase = modelBase.timeOffsetSec(:);
jitterInfo.timeOffsetSecEval = modelJit.timeOffsetSec(:);
jitterInfo.fdAliasStepHzBase = localResolveAliasStep(modelBase);
jitterInfo.fdAliasStepHzEval = localResolveAliasStep(modelJit);
jitterInfo.jitterStdSec = localGetFieldOrDefault(scanCfg, 'jitterStdSec', NaN);
end

function modelJit = localBuildJitteredModel(modelBase, scanCfg)
%LOCALBUILDJITTEREDMODEL Build a probe-only jittered evaluator model.

modelJit = modelBase;
timeOffsetSec = modelBase.timeOffsetSec(:);
refFrameIdx = modelBase.refFrameIdx;
numFrame = numel(timeOffsetSec);
if numFrame <= 1
  return;
end

jitterStdSec = localGetFieldOrDefault(scanCfg, 'jitterStdSec', 0);
if ~(isscalar(jitterStdSec) && isfinite(jitterStdSec) && jitterStdSec > 0)
  return;
end

switch lower(char(localGetFieldOrDefault(scanCfg, 'jitterMode', 'zeroMeanRandom')))
  case 'zeromeanrandom'
    jitterVec = jitterStdSec * randn(numFrame, 1);
    jitterVec(refFrameIdx) = 0;
    jitterVec = jitterVec - mean(jitterVec);
    jitterVec(refFrameIdx) = 0;
  case 'alternating'
    jitterVec = jitterStdSec * (-1) .^ ((1:numFrame).' - refFrameIdx);
    jitterVec(refFrameIdx) = 0;
  otherwise
    error('doaDopplerDynFdRefCombScan:InvalidJitterMode', ...
      'Unsupported jitterMode: %s', string(scanCfg.jitterMode));
end

modelJit.timeOffsetSec = timeOffsetSec + jitterVec;
diffTimeSec = diff(modelJit.timeOffsetSec);
validDiff = diffTimeSec(isfinite(diffTimeSec) & diffTimeSec > 0);
if ~isempty(validDiff)
  modelJit.fdAliasStepHz = 1 / median(validDiff);
end
end

function probeEval = localEvalProbeAtPoint(model, pointCur, probeTemplate)
%LOCALEVALPROBEATPOINT Evaluate one probe point through the shared helper.

probeEval = evalDoaDopplerMfProbePoint(model, pointCur, struct('probeTemplate', probeTemplate, 'debugEnable', false));
end

function reciprocalPeak = localBuildReciprocalPeak(objVec, minObj)
%LOCALBUILDRECIPROCALPEAK Build one reciprocal peak curve from delta objective.

deltaObj = objVec - minObj;
positiveDelta = deltaObj(isfinite(deltaObj) & deltaObj > 0);
if isempty(positiveDelta)
  floorVal = 1;
else
  floorVal = max(min(positiveDelta), 1e-12);
end
reciprocalPeak = 1 ./ max(deltaObj + floorVal, floorVal);
reciprocalPeak = reciprocalPeak / max(reciprocalPeak);
end

function foldedOffset = localFoldOffset(deltaFdRef, aliasStepHz)
%LOCALFOLDOFFSET Fold offsets into one alias period.

if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  foldedOffset = deltaFdRef(:);
  return;
end
foldedOffset = mod(deltaFdRef(:) + aliasStepHz / 2, aliasStepHz) - aliasStepHz / 2;
end

function aliasStepHz = localResolveAliasStep(model)
%LOCALRESOLVEALIASSTEP Resolve one alias step in Hz.

aliasStepHz = localGetFieldOrDefault(model, 'fdAliasStepHz', NaN);
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  diffTimeSec = diff(model.timeOffsetSec(:));
  diffTimeSec = diffTimeSec(isfinite(diffTimeSec) & diffTimeSec > 0);
  if isempty(diffTimeSec)
    aliasStepHz = NaN;
  else
    aliasStepHz = 1 / median(diffTimeSec);
  end
end
end

function centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown)
%LOCALBUILDCENTERREGISTRY Build named center points.

centerRegistry = struct();
centerRegistry.truth = localBuildPointStruct( ...
  truth.latlonTrueDeg(:), truth.fdRefFit, truth.fdRateFit);

staticEst = localGetFieldOrDefault(bestStaticMsCase, 'estResult', struct());
centerRegistry.staticSeed = localBuildPointStruct( ...
  localResolveDoaParam(staticEst, truth.latlonTrueDeg(:)), ...
  localResolveScalarField(staticEst, 'fdRefEst', truth.fdRefFit), ...
  truth.fdRateFit);

knownEst = localGetFieldOrDefault(caseDynMsKnown, 'estResult', struct());
centerRegistry.cpKnownSeed = localBuildPointStruct( ...
  localResolveDoaParam(knownEst, centerRegistry.staticSeed.latlon(:)), ...
  localResolveScalarField(knownEst, 'fdRefEst', centerRegistry.staticSeed.fdRef), ...
  localResolveScalarField(knownEst, 'fdRateEst', truth.fdRateFit));

unknownEst = localGetFieldOrDefault(caseDynMsUnknown, 'estResult', struct());
centerRegistry.finalEstimateKnown = centerRegistry.cpKnownSeed;
centerRegistry.finalEstimateUnknown = localBuildPointStruct( ...
  localResolveDoaParam(unknownEst, centerRegistry.staticSeed.latlon(:)), ...
  localResolveScalarField(unknownEst, 'fdRefEst', centerRegistry.staticSeed.fdRef), ...
  localResolveScalarField(unknownEst, 'fdRateEst', truth.fdRateFit));
end

function centerNameList = localResolveCenterNameList(scanCfg, modeTag)
%LOCALRESOLVECENTERNAMELIST Resolve the center list for one mode.

if strcmpi(modeTag, 'known')
  centerNameList = reshape(string(scanCfg.centerListKnown), [], 1);
else
  centerNameList = reshape(string(scanCfg.centerListUnknown), [], 1);
end
end

function centerPoint = localResolveCenterPoint(centerRegistry, modeTag, centerName)
%LOCALRESOLVECENTERPOINT Resolve one named scan center.

switch lower(char(centerName))
  case 'truth'
    centerPoint = centerRegistry.truth;
  case 'staticseed'
    centerPoint = centerRegistry.staticSeed;
  case 'cpknownseed'
    centerPoint = centerRegistry.cpKnownSeed;
  case 'finalestimate'
    if strcmpi(modeTag, 'known')
      centerPoint = centerRegistry.finalEstimateKnown;
    else
      centerPoint = centerRegistry.finalEstimateUnknown;
    end
  otherwise
    error('doaDopplerDynFdRefCombScan:UnknownCenter', ...
      'Unsupported center name: %s', centerName);
end
end

function pointInfo = localBuildPointStruct(doaParam, fdRef, fdRate)
%LOCALBUILDPOINTSTRUCT Build one compact point representation.

doaParam = reshape(doaParam, [], 1);
pointInfo = struct();
pointInfo.lat = doaParam(1);
pointInfo.lon = doaParam(2);
pointInfo.latlon = doaParam(1:2);
pointInfo.fdRef = fdRef;
pointInfo.fdRate = fdRate;
end

function localDisplayModeSummary(modeLabel, modeResult)
%LOCALDISPLAYMODESUMMARY Print a compact one-line summary per center.

for iCenter = 1:numel(modeResult.centerScan)
  scanInfo = modeResult.centerScan(iCenter).uniform;
  fprintf('%s | %s | uniform: minDeltaFd = %.6f Hz, minAliasIdx = %d, centerDeltaObj = %.6g\n', ...
    modeLabel, modeResult.centerScan(iCenter).centerName, ...
    scanInfo.minDeltaFdRef, scanInfo.minAliasIndex, scanInfo.centerDeltaObj);
  if isfield(modeResult.centerScan(iCenter), 'jittered') && ...
      isstruct(modeResult.centerScan(iCenter).jittered) && ...
      ~isempty(fieldnames(modeResult.centerScan(iCenter).jittered))
    scanInfoJit = modeResult.centerScan(iCenter).jittered;
    fprintf('%s | %s | jittered: minDeltaFd = %.6f Hz, minAliasIdx = %d, centerDeltaObj = %.6g\n', ...
      modeLabel, modeResult.centerScan(iCenter).centerName, ...
      scanInfoJit.minDeltaFdRef, scanInfoJit.minAliasIndex, scanInfoJit.centerDeltaObj);
  end
end
end

function localPlotModeDeltaFigure(modeResult, modeLabel)
%LOCALPLOTMODEDELTAFIGURE Plot delta objective versus fdRef offset.

numCenter = numel(modeResult.centerScan);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s fdRef delta scan', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  scanU = modeResult.centerScan(iCenter).uniform;
  plot(scanU.deltaFdRef, scanU.deltaObjVec, 'LineWidth', 1.4, 'DisplayName', 'uniform');
  if isfield(modeResult.centerScan(iCenter), 'jittered') && ...
      isstruct(modeResult.centerScan(iCenter).jittered) && ...
      ~isempty(fieldnames(modeResult.centerScan(iCenter).jittered))
    scanJ = modeResult.centerScan(iCenter).jittered;
    plot(scanJ.deltaFdRef, scanJ.deltaObjVec, '--', 'LineWidth', 1.2, 'DisplayName', 'jittered');
  end
  xline(0, ':', 'truth/center', 'LabelVerticalAlignment', 'bottom');
  localOverlayAliasLines(modeResult.aliasStepHz, scanU.deltaFdRef);
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('\Delta objective');
  title(sprintf('%s | %s | aliasStep %.6f Hz', ...
    modeLabel, modeResult.centerScan(iCenter).centerName, modeResult.aliasStepHz), ...
    'Interpreter', 'none');
  legend('Location', 'best');
end
end

function localPlotModePeakFigure(modeResult, modeLabel)
%LOCALPLOTMODEPEAKFIGURE Plot reciprocal peak metric versus fdRef offset.

numCenter = numel(modeResult.centerScan);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s fdRef peak scan', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  scanU = modeResult.centerScan(iCenter).uniform;
  plot(scanU.deltaFdRef, scanU.reciprocalPeak, 'LineWidth', 1.4, 'DisplayName', 'uniform');
  if isfield(modeResult.centerScan(iCenter), 'jittered') && ...
      isstruct(modeResult.centerScan(iCenter).jittered) && ...
      ~isempty(fieldnames(modeResult.centerScan(iCenter).jittered))
    scanJ = modeResult.centerScan(iCenter).jittered;
    plot(scanJ.deltaFdRef, scanJ.reciprocalPeak, '--', 'LineWidth', 1.2, 'DisplayName', 'jittered');
  end
  xline(0, ':', 'truth/center', 'LabelVerticalAlignment', 'bottom');
  localOverlayAliasLines(modeResult.aliasStepHz, scanU.deltaFdRef);
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('normalized reciprocal peak');
  title(sprintf('%s | %s | reciprocal peak', ...
    modeLabel, modeResult.centerScan(iCenter).centerName), ...
    'Interpreter', 'none');
  legend('Location', 'best');
end
end

function localPlotModeFoldedFigure(modeResult, modeLabel)
%LOCALPLOTMODEFOLDEDFIGURE Plot folded offsets modulo the alias step.

numCenter = numel(modeResult.centerScan);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s fdRef folded scan', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  scanU = modeResult.centerScan(iCenter).uniform;
  scatter(scanU.foldedOffset, scanU.deltaObjVec, 18, 'filled', 'DisplayName', 'uniform');
  if isfield(modeResult.centerScan(iCenter), 'jittered') && ...
      isstruct(modeResult.centerScan(iCenter).jittered) && ...
      ~isempty(fieldnames(modeResult.centerScan(iCenter).jittered))
    scanJ = modeResult.centerScan(iCenter).jittered;
    scatter(scanJ.foldedOffset, scanJ.deltaObjVec, 18, 'DisplayName', 'jittered');
  end
  grid on;
  xlabel(sprintf('wrapped \\Delta fdRef mod %.6f Hz', modeResult.aliasStepHz));
  ylabel('\Delta objective');
  title(sprintf('%s | %s | folded modulo scan', ...
    modeLabel, modeResult.centerScan(iCenter).centerName), ...
    'Interpreter', 'none');
  legend('Location', 'best');
end
end

function localOverlayAliasLines(aliasStepHz, xGrid)
%LOCALOVERLAYALIASLINES Overlay integer alias-step guides.

if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  return;
end

xMin = min(xGrid);
xMax = max(xGrid);
numSide = ceil(max(abs([xMin, xMax])) / aliasStepHz) + 1;
for k = -numSide:numSide
  xVal = k * aliasStepHz;
  if xVal < xMin || xVal > xMax
    continue;
  end
  xline(xVal, ':');
end
end

function localMaybeStartParpool(scanCfg)
%LOCALMAYBESTARTPARPOOL Start a parallel pool when requested.

if ~localGetFieldOrDefault(scanCfg, 'useParfor', false)
  return;
end
if ~localGetFieldOrDefault(scanCfg, 'autoStartParpool', false)
  return;
end
poolObj = gcp('nocreate');
if isempty(poolObj)
  parpool();
end
end

function useParfor = localResolveUseParfor(scanCfg, numGrid)
%LOCALRESOLVEUSEPARFOR Decide whether to use parfor for the scan.

useParfor = localGetFieldOrDefault(scanCfg, 'useParfor', false);
if ~useParfor
  return;
end
if numGrid < 64
  useParfor = false;
  return;
end
poolObj = gcp('nocreate');
useParfor = ~isempty(poolObj);
end

function idx = localFindNearestGridIndex(valGrid, valTarget)
%LOCALFINDNEARESTGRIDINDEX Return the nearest grid index for one scalar.

[~, idx] = min(abs(valGrid - valTarget));
end

function halfWidthHz = localResolveHalfWidth(scanCfg, aliasStepHz)
%LOCALRESOLVEHALFWIDTH Resolve the fdRef scan half width.

halfWidthHz = localGetFieldOrDefault(scanCfg, 'scanHalfWidthHz', []);
if isempty(halfWidthHz)
  numAliasSide = localGetFieldOrDefault(scanCfg, 'numAliasSide', 2);
  halfWidthHz = numAliasSide * aliasStepHz;
end
end

function doaParam = localResolveDoaParam(estResult, fallbackDoa)
%LOCALRESOLVEDOAPARAM Resolve one 2x1 DoA vector with fallback.

doaParam = reshape(fallbackDoa, [], 1);
cand = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(cand) == 2 && all(isfinite(cand))
  doaParam = cand;
end
end

function scalarValue = localResolveScalarField(estResult, fieldName, fallbackValue)
%LOCALRESOLVESCALARFIELD Resolve one finite scalar field with fallback.

scalarValue = fallbackValue;
cand = localGetFieldOrDefault(estResult, fieldName, []);
if isscalar(cand) && isnumeric(cand) && isfinite(cand)
  scalarValue = cand;
end
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

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
