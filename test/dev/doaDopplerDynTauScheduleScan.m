% Development script for diagnosing how the tau schedule changes the fdRef
% comb structure in MF dynamic DoA-Doppler estimation.
%
% The script keeps the observation block length fixed at the paper-aligned
% standard block length (Np = 2112) and only changes model.timeOffsetSec.
% For each tau schedule, it scans fdRef around selected centers and compares
% how the 1 / Tf comb changes.
clear(); close all; clc;

%% Parameters
numUsr = 1;
numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = ceil(numFrame / 2);

sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
standardBlockLen = 2112;
numSym = round(standardBlockLen / osf);
carrierFreq = 18e9;
wavelength = 299792458 / carrierFreq;
rng(253);

elemSpace = wavelength / 2;
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

doaLocalHalfWidthKnown = [0.003; 0.003];
doaLocalHalfWidthUnknown = [0.002; 0.002];

scanCfg = struct();
scanCfg.runKnown = true;
scanCfg.runUnknown = true;
scanCfg.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanCfg.centerListUnknown = ["truth"; "cpKnownSeed"; "finalEstimate"];
scanCfg.tauModeList = ["uniform"; "zeroMeanJitter"; "alternatingJitter"; "oneGapLate"];
scanCfg.jitterStdSec = 2e-6;
scanCfg.alternatingJitterSec = 5e-6;
scanCfg.oneGapLateSec = 25e-6;
scanCfg.numAliasSide = 2;
scanCfg.numGrid = 801;
scanCfg.scanHalfWidthHz = [];
scanCfg.useParfor = true;
scanCfg.autoStartParpool = true;
scanCfg.minGridForParfor = 161;
scanCfg.showDeltaFigure = true;
scanCfg.showPeakFigure = true;
scanCfg.showFoldedFigure = true;
scanCfg.showToothTable = true;
scanCfg.saveTag = "doaDopplerDynTauScheduleScan";

fprintf('[TauScheduleScan] Building scene and selecting satellites ...\n');
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
if sceneSeq.numSat ~= 2
  error('doaDopplerDynTauScheduleScan:InvalidNumSat', ...
    'This script expects exactly two selected satellites.');
end
otherSatIdxLocal = setdiff(1:sceneSeq.numSat, refSatIdxLocal);
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelength);
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
truth.fdRefTrueHz = truth.fdRefSeries(sceneSeq.refFrameIdx);
truth.fdRateTrueHzPerSec = truth.fdRateFit;
truth.fdSatTrueHz = reshape(truth.fdSatSeries(:, sceneSeq.refFrameIdx), [], 1);
truth.deltaFdTrueHz = reshape(truth.deltaFdSeries(:, sceneSeq.refFrameIdx), [], 1);

fdRange = expandRangeToTruth(fdRange, [truth.fdRefFit; truth.fdSatTrueHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth(fdRateRange, fdRateTruthCand, 0.1, 5e2);

fprintf('[TauScheduleScan] Generating standard pilot waveform and snapshots ...\n');
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
pilotWave = localTrimPilotWave(pilotWave, standardBlockLen);
if localGetPilotLength(pilotWave) ~= standardBlockLen
  error('doaDopplerDynTauScheduleScan:InvalidPilotLength', ...
    'Generated pilot length (%d) does not match the standard block length (%d).', ...
    localGetPilotLength(pilotWave), standardBlockLen);
end

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

fprintf('[TauScheduleScan] Building static seed and dynamic centers ...\n');
staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelength, pilotWave, ...
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
dynMsKnownOpt.initDoaHalfWidth = doaLocalHalfWidthKnown;
dynMsKnownOpt.enableFdAliasUnwrap = true;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsKnownOpt, true, [], initParamStaticMs);

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = doaLocalHalfWidthUnknown;
dynMsUnknownOpt.enableFdAliasUnwrap = true;

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  fdRange, fdRateRange, optVerbose, dynMsUnknownOpt, false, [], msUnknownCand);

centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown);

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
  'numPilotSample', localGetPilotLength(pilotWave), ...
  'refFrameIdx', sceneSeq.refFrameIdx, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'refStateSource', string(refState.source));

if scanCfg.autoStartParpool && scanCfg.useParfor
  localMaybeStartParpool();
end

if scanCfg.runKnown
  fprintf('[TauScheduleScan] Building cp-known model ...\n');
  modelOptKnown = baseDynOpt;
  modelOptKnown.fdRateMode = 'known';
  modelOptKnown.fdRateKnown = truth.fdRateFit;
  [modelKnown, ~, ~, modelOptKnown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, [], modelOptKnown);
  scanResult.cpKnown = localRunTauScheduleSet(modelKnown, centerRegistry, scanCfg, 'known');
  scanResult.cpKnown.modelOpt = modelOptKnown;
end

if scanCfg.runUnknown
  fprintf('[TauScheduleScan] Building cp-unknown model ...\n');
  modelOptUnknown = baseDynOpt;
  modelOptUnknown.fdRateMode = 'unknown';
  [modelUnknown, ~, ~, modelOptUnknown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, fdRateRange, modelOptUnknown);
  scanResult.cpUnknown = localRunTauScheduleSet(modelUnknown, centerRegistry, scanCfg, 'unknown');
  scanResult.cpUnknown.modelOpt = modelOptUnknown;
end

if isfield(scanResult, 'cpKnown')
  disp('========== CP-K tau schedule summary ==========' );
  localDisplayModeSummary('CP-K', scanResult.cpKnown);
end
if isfield(scanResult, 'cpUnknown')
  disp('========== CP-U tau schedule summary ==========' );
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
if localGetFieldOrDefault(scanCfg, 'showToothTable', true)
  if isfield(scanResult, 'cpKnown')
    disp('========== Tau-schedule tooth table (CP-K) ==========' );
    disp(localBuildTauToothTable(scanResult.cpKnown));
  end
  if isfield(scanResult, 'cpUnknown')
    disp('========== Tau-schedule tooth table (CP-U) ==========' );
    disp(localBuildTauToothTable(scanResult.cpUnknown));
  end
end

scanResult.utcRun = datetime('now', 'TimeZone', 'local');
scanResult.scriptTag = scanCfg.saveTag;
saveFile = fullfile(fileparts(mfilename('fullpath')), ...
  sprintf('%s_%s.mat', scanCfg.saveTag, datestr(now, 'yyyymmdd-HHMMSS')));
save(saveFile, 'scanResult');
fprintf('[TauScheduleScan] Saved scanResult to %s\n', saveFile);

%% Local functions
function modeResult = localRunTauScheduleSet(modelBase, centerRegistry, scanCfg, modeTag)
centerNameList = localResolveCenterNameList(scanCfg, modeTag);
tauModeList = reshape(string(scanCfg.tauModeList), [], 1);
modeResult = struct();
modeResult.modeTag = string(modeTag);
modeResult.centerNameList = centerNameList;
modeResult.tauModeList = tauModeList;
modeResult.baseAliasStepHz = localResolveAliasStep(modelBase);
modeResult.tauScan = repmat(struct( ...
  'tauMode', "", ...
  'modelInfo', struct(), ...
  'centerScan', []), numel(tauModeList), 1);

for iTau = 1:numel(tauModeList)
  tauMode = tauModeList(iTau);
  modelTau = localBuildTauModifiedModel(modelBase, tauMode, scanCfg);
  tauScan = struct();
  tauScan.tauMode = tauMode;
  tauScan.modelInfo = localBuildTauModelInfo(modelBase, modelTau, tauMode);
  tauScan.centerScan = repmat(struct( ...
    'centerName', "", ...
    'basePoint', struct(), ...
    'lineInfo', struct()), numel(centerNameList), 1);

  for iCenter = 1:numel(centerNameList)
    centerName = centerNameList(iCenter);
    basePoint = localResolveCenterPoint(centerRegistry, modeTag, centerName);
    halfWidthHz = localResolveHalfWidth(scanCfg, localResolveAliasStep(modelTau));
    fdRefGrid = linspace(basePoint.fdRef - halfWidthHz, basePoint.fdRef + halfWidthHz, scanCfg.numGrid);

    fprintf('[TauScheduleScan]   mode=%s | center=%s | grid=%d\n', tauMode, centerName, numel(fdRefGrid));
    centerScan = struct();
    centerScan.centerName = centerName;
    centerScan.basePoint = basePoint;
    centerScan.lineInfo = localScanFdRefLine(modelTau, basePoint, fdRefGrid, scanCfg);
    tauScan.centerScan(iCenter) = centerScan;
  end

  modeResult.tauScan(iTau) = tauScan;
end
end

function lineInfo = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanCfg)
modelScan = modelBase;
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
end

function [objCur, residualCur, fdRefCur, fdRateCur, cohCur, objSatCur] = ...
    localEvaluateScanPoint(modelScan, basePoint, fdRefCurIn, probeTemplate)
pointCur = basePoint;
pointCur.fdRef = fdRefCurIn;
probeEval = evalDoaDopplerMfProbePoint(modelScan, pointCur, ...
  struct('probeTemplate', probeTemplate, 'debugEnable', false));
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

function modelTau = localBuildTauModifiedModel(modelBase, tauMode, scanCfg)
modelTau = modelBase;
baseTau = modelBase.timeOffsetSec(:);
refFrameIdx = modelBase.refFrameIdx;
modelTau.timeOffsetSec = localBuildTauOffsets(baseTau, refFrameIdx, tauMode, scanCfg);
modelTau.fdAliasStepHz = localResolveAliasStep(modelTau);
end

function tauOut = localBuildTauOffsets(baseTau, refFrameIdx, tauMode, scanCfg)
baseTau = baseTau(:);
numFrame = numel(baseTau);
if numFrame <= 1
  tauOut = baseTau;
  return;
end
nominalStep = median(diff(baseTau));
stepFwd = nominalStep * ones(numFrame - refFrameIdx, 1);
stepBwd = nominalStep * ones(refFrameIdx - 1, 1);

switch lower(char(tauMode))
  case 'uniform'
    % keep nominal steps
  case 'zeromeanjitter'
    stream = RandStream('mt19937ar', 'Seed', 3101);
    stepFwd = stepFwd + scanCfg.jitterStdSec * randn(stream, size(stepFwd));
    stepBwd = stepBwd + scanCfg.jitterStdSec * randn(stream, size(stepBwd));
  case 'alternatingjitter'
    altFwd = (-1) .^ (1:numel(stepFwd)).';
    altBwd = (-1) .^ (1:numel(stepBwd)).';
    stepFwd = stepFwd + scanCfg.alternatingJitterSec * altFwd;
    stepBwd = stepBwd + scanCfg.alternatingJitterSec * altBwd;
  case 'onegaplate'
    if ~isempty(stepFwd)
      idx = max(1, ceil(numel(stepFwd) / 2));
      stepFwd(idx) = stepFwd(idx) + scanCfg.oneGapLateSec;
    end
  otherwise
    error('doaDopplerDynTauScheduleScan:UnknownTauMode', ...
      'Unsupported tauMode: %s', tauMode);
end

minStep = max(abs(nominalStep) / 20, 1e-7);
stepFwd = max(stepFwd, minStep);
stepBwd = max(stepBwd, minStep);

tauOut = zeros(numFrame, 1);
for i = refFrameIdx + 1:numFrame
  tauOut(i) = tauOut(i - 1) + stepFwd(i - refFrameIdx);
end
for i = refFrameIdx - 1:-1:1
  tauOut(i) = tauOut(i + 1) - stepBwd(refFrameIdx - i);
end
end

function info = localBuildTauModelInfo(modelBase, modelTau, tauMode)
info = struct();
info.tauMode = string(tauMode);
info.timeOffsetSecBase = modelBase.timeOffsetSec(:);
info.timeOffsetSecEval = modelTau.timeOffsetSec(:);
info.fdAliasStepHzBase = localResolveAliasStep(modelBase);
info.fdAliasStepHzEval = localResolveAliasStep(modelTau);
end

function reciprocalPeak = localBuildReciprocalPeak(objVec, minObj)
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
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  foldedOffset = deltaFdRef(:);
  return;
end
foldedOffset = mod(deltaFdRef(:) + aliasStepHz / 2, aliasStepHz) - aliasStepHz / 2;
end

function aliasStepHz = localResolveAliasStep(model)
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
if strcmpi(modeTag, 'known')
  centerNameList = reshape(string(scanCfg.centerListKnown), [], 1);
else
  centerNameList = reshape(string(scanCfg.centerListUnknown), [], 1);
end
end

function centerPoint = localResolveCenterPoint(centerRegistry, modeTag, centerName)
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
    error('doaDopplerDynTauScheduleScan:UnknownCenter', ...
      'Unsupported center name: %s', centerName);
end
end

function pointInfo = localBuildPointStruct(doaParam, fdRef, fdRate)
doaParam = reshape(doaParam, [], 1);
pointInfo = struct();
pointInfo.lat = doaParam(1);
pointInfo.lon = doaParam(min(2, numel(doaParam)));
pointInfo.latlon = doaParam;
pointInfo.doaParam = doaParam;
pointInfo.fdRef = fdRef;
pointInfo.fdRate = fdRate;
end

function doaParam = localResolveDoaParam(estResult, fallbackDoa)
doaParam = localGetFieldOrDefault(estResult, 'doaParamEst', []);
if isempty(doaParam)
  doaParam = fallbackDoa(:);
else
  doaParam = reshape(doaParam, [], 1);
end
end

function val = localResolveScalarField(s, fieldName, fallbackVal)
val = localGetFieldOrDefault(s, fieldName, fallbackVal);
if isempty(val)
  val = fallbackVal;
end
end

function halfWidthHz = localResolveHalfWidth(scanCfg, aliasStepHz)
halfWidthHz = localGetFieldOrDefault(scanCfg, 'scanHalfWidthHz', []);
if isempty(halfWidthHz)
  halfWidthHz = localGetFieldOrDefault(scanCfg, 'numAliasSide', 2) * aliasStepHz;
end
if ~(isscalar(halfWidthHz) && isfinite(halfWidthHz) && halfWidthHz > 0)
  halfWidthHz = 2 * max(abs(aliasStepHz), 1);
end
end

function useParfor = localResolveUseParfor(scanCfg, numGrid)
useParfor = logical(localGetFieldOrDefault(scanCfg, 'useParfor', false));
if ~useParfor
  return;
end
minGridForParfor = localGetFieldOrDefault(scanCfg, 'minGridForParfor', inf);
useParfor = numGrid >= minGridForParfor;
end

function localMaybeStartParpool()
poolObj = gcp('nocreate');
if isempty(poolObj)
  parpool('Processes');
end
end

function idx = localFindNearestGridIndex(gridVec, targetVal)
[~, idx] = min(abs(gridVec(:) - targetVal));
end

function tbl = localBuildTauToothTable(modeResult)
rows = struct([]);
for iTau = 1:numel(modeResult.tauScan)
  tauScan = modeResult.tauScan(iTau);
  for iCenter = 1:numel(tauScan.centerScan)
    centerScan = tauScan.centerScan(iCenter);
    lineInfo = centerScan.lineInfo;
    row.tauMode = tauScan.tauMode; %#ok<STRNU>
    row.centerName = centerScan.centerName;
    row.aliasStepHz = lineInfo.aliasStepHz;
    row.centerObj = lineInfo.centerObj;
    row.minObj = lineInfo.minObj;
    row.centerDeltaObj = lineInfo.centerDeltaObj;
    row.minDeltaFdRef = lineInfo.minDeltaFdRef;
    row.minAliasIndex = lineInfo.minAliasIndex;
    if isempty(rows)
      rows = row;
    else
      rows(end + 1, 1) = row; %#ok<AGROW>
    end
  end
end
if isempty(rows)
  tbl = table();
else
  tbl = struct2table(rows);
end
end

function localDisplayModeSummary(displayName, modeResult)
for iTau = 1:numel(modeResult.tauScan)
  tauScan = modeResult.tauScan(iTau);
  fprintf('%s | tauMode=%s | aliasStep=%.6f Hz\n', ...
    displayName, tauScan.tauMode, tauScan.modelInfo.fdAliasStepHzEval);
  for iCenter = 1:numel(tauScan.centerScan)
    lineInfo = tauScan.centerScan(iCenter).lineInfo;
    fprintf('  center=%s | centerDeltaObj=%.6g | minDeltaFdRef=%.6f Hz | aliasIndex=%d\n', ...
      tauScan.centerScan(iCenter).centerName, ...
      lineInfo.centerDeltaObj, lineInfo.minDeltaFdRef, lineInfo.minAliasIndex);
  end
end
end

function localPlotModeDeltaFigure(modeResult, modeName)
figure('Name', sprintf('%s tau-schedule delta objective', modeName));
tiledlayout(numel(modeResult.tauScan), numel(modeResult.centerNameList), 'TileSpacing', 'compact', 'Padding', 'compact');
for iTau = 1:numel(modeResult.tauScan)
  tauScan = modeResult.tauScan(iTau);
  for iCenter = 1:numel(tauScan.centerScan)
    nexttile;
    lineInfo = tauScan.centerScan(iCenter).lineInfo;
    plot(lineInfo.deltaFdRef, lineInfo.deltaObjVec, 'LineWidth', 1.1);
    hold on;
    xline(0, '--');
    if isfinite(lineInfo.aliasStepHz)
      for k = -2:2
        xline(k * lineInfo.aliasStepHz, ':');
      end
    end
    hold off;
    grid on;
    xlabel('\Delta f_{ref} (Hz)');
    ylabel('\Delta objective');
    title(sprintf('%s | %s', tauScan.tauMode, tauScan.centerScan(iCenter).centerName), 'Interpreter', 'none');
  end
end
end

function localPlotModePeakFigure(modeResult, modeName)
figure('Name', sprintf('%s tau-schedule reciprocal peak', modeName));
tiledlayout(numel(modeResult.tauScan), numel(modeResult.centerNameList), 'TileSpacing', 'compact', 'Padding', 'compact');
for iTau = 1:numel(modeResult.tauScan)
  tauScan = modeResult.tauScan(iTau);
  for iCenter = 1:numel(tauScan.centerScan)
    nexttile;
    lineInfo = tauScan.centerScan(iCenter).lineInfo;
    plot(lineInfo.deltaFdRef, lineInfo.reciprocalPeak, 'LineWidth', 1.1);
    hold on;
    xline(0, '--');
    if isfinite(lineInfo.aliasStepHz)
      for k = -2:2
        xline(k * lineInfo.aliasStepHz, ':');
      end
    end
    hold off;
    grid on;
    xlabel('\Delta f_{ref} (Hz)');
    ylabel('normalized reciprocal peak');
    title(sprintf('%s | %s', tauScan.tauMode, tauScan.centerScan(iCenter).centerName), 'Interpreter', 'none');
  end
end
end

function localPlotModeFoldedFigure(modeResult, modeName)
figure('Name', sprintf('%s tau-schedule folded modulo', modeName));
tiledlayout(numel(modeResult.tauScan), numel(modeResult.centerNameList), 'TileSpacing', 'compact', 'Padding', 'compact');
for iTau = 1:numel(modeResult.tauScan)
  tauScan = modeResult.tauScan(iTau);
  for iCenter = 1:numel(tauScan.centerScan)
    nexttile;
    lineInfo = tauScan.centerScan(iCenter).lineInfo;
    scatter(lineInfo.foldedOffset, lineInfo.deltaObjVec, 10, lineInfo.deltaFdRef, 'filled');
    grid on;
    xlabel('folded \Delta f_{ref} (Hz)');
    ylabel('\Delta objective');
    title(sprintf('%s | %s', tauScan.tauMode, tauScan.centerScan(iCenter).centerName), 'Interpreter', 'none');
    colorbar;
  end
end
end

function pilotWave = localTrimPilotWave(pilotWaveFull, blockLen)
if size(pilotWaveFull, 1) >= blockLen
  pilotWave = pilotWaveFull(1:blockLen, :);
elseif size(pilotWaveFull, 2) >= blockLen
  pilotWave = pilotWaveFull(:, 1:blockLen);
else
  error('doaDopplerDynTauScheduleScan:InvalidPilotLength', ...
    'Requested block length %d exceeds both dimensions of the generated pilot.', blockLen);
end
end

function len = localGetPilotLength(pilotWave)
if size(pilotWave, 1) >= size(pilotWave, 2)
  len = size(pilotWave, 1);
else
  len = size(pilotWave, 2);
end
end

function val = localGetFieldOrDefault(s, fieldName, defaultVal)
if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
