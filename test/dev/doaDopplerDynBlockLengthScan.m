% Development script for diagnosing how block length affects the fdRef comb
% structure in MF continuous-phase DoA-Doppler estimation.
%
% The script keeps the frame interval fixed and only changes the number of
% pilot samples retained from each frame. For each block length it rebuilds
% the single-frame / multi-frame views, reruns the compact dynamic cases to
% obtain representative centers, and then performs one-dimensional fdRef
% comb scans through the shared probe helper evalDoaDopplerMfProbePoint.
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

doBaseHalfWidth = [0.01; 0.01];

doaLocalHalfWidthKnown = [0.003; 0.003];
doaLocalHalfWidthUnknown = [0.002; 0.002];

scanCfg = struct();
scanCfg.runKnown = true;
scanCfg.runUnknown = true;
scanCfg.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanCfg.centerListUnknown = ["truth"; "cpKnownSeed"; "finalEstimate"];
scanCfg.blockLenList = [528; 1056; 2112];
scanCfg.numAliasSide = 2;
scanCfg.numGrid = 801;
scanCfg.scanHalfWidthHz = [];
scanCfg.useParfor = true;
scanCfg.autoStartParpool = true;
scanCfg.minGridForParfor = 161;
scanCfg.showDeltaFigure = true;
scanCfg.showFoldedFigure = true;
scanCfg.showPeakFigure = true;
scanCfg.showToothTable = true;
scanCfg.saveTag = "doaDopplerDynBlockLengthScan";

%% Select satellites and build scene
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
if sceneSeq.numSat ~= 2
  error('doaDopplerDynBlockLengthScan:InvalidNumSat', ...
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

%% Pilot waveform and snapshots
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWaveFull, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

numPilotSampleFull = localGetSignalLength(pilotWaveFull);
if numPilotSampleFull ~= standardBlockLen
  error('doaDopplerDynBlockLengthScan:UnexpectedPilotLength', ...
    ['Generated pilot length is %d samples, but this script is aligned to the paper ' ...
     'definition N_p = %d samples per frame. Please keep genPilotSymbol / genPilotWaveform ' ...
     'consistent with the standard block definition before running the scan.'], ...
    numPilotSampleFull, standardBlockLen);
end

blockLenListReq = reshape(scanCfg.blockLenList, [], 1);
blockLenList = blockLenListReq(blockLenListReq > 0 & blockLenListReq <= standardBlockLen);
blockLenList = unique(round(blockLenList), 'stable');
if isempty(blockLenList)
  error('doaDopplerDynBlockLengthScan:EmptyBlockLenList', ...
    'Requested blockLenList is empty after clipping to the standard block length.');
end
scanCfg.blockLenList = blockLenList;

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
[rxSigCellFull, ~, ~, ~, ~] = genMultiFrameSnapshots( ...
  sceneSeq, pilotWaveFull, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGainCell, snapOpt);

%% Auto-start parpool if requested
localEnsureParpool(scanCfg);

%% Run block-length scans
scanResult = struct();
scanResult.truth = truth;
scanResult.scanCfg = scanCfg;
scanResult.frameInfo = struct( ...
  'frameIntvlSec', frameIntvlSec, ...
  'timeOffsetSec', sceneSeq.timeOffsetSec(:), ...
  'numFrame', sceneSeq.numFrame, ...
  'numPilotSampleFull', numPilotSampleFull, ...
  'refFrameIdx', sceneSeq.refFrameIdx, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'refStateSource', string(refState.source));
scanResult.blockLenList = blockLenList;

blockResult = repmat(struct( ...
  'blockLen', NaN, ...
  'pilotWave', [], ...
  'viewMs', struct(), ...
  'centerRegistry', struct(), ...
  'cpKnown', struct(), ...
  'cpUnknown', struct()), numel(blockLenList), 1);

for iBlock = 1:numel(blockLenList)
  blockLen = blockLenList(iBlock);
  fprintf('Running block-length scan: %d samples\n', blockLen);

  pilotWave = localTrimPilotWave(pilotWaveFull, blockLen);
  rxSigCell = localTrimRxSigCell(rxSigCellFull, blockLen);
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

  staticBaseOpt = struct();
  staticBaseOpt.useLogObjective = true;
  doaOnlyOpt = struct();
  doaOnlyOpt.useLogObjective = true;

  caseBundle = buildDoaDopplerStaticTransitionBundle( ...
    viewRefOnly, viewOtherOnly, viewMs, wavelength, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
    otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, ...
    weightSweepAlpha, doBaseHalfWidth);

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
  dynMsKnownOpt.enableFdAliasUnwrap = false;

  caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
    viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    fdRange, fdRateRange, optVerbose, dynMsKnownOpt, true, [], initParamStaticMs);

  initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
  initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, truth.fdRateFit);

  dynMsUnknownOpt = dynBaseOpt;
  dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
  dynMsUnknownOpt.initDoaHalfWidth = doaLocalHalfWidthUnknown;
  dynMsUnknownOpt.enableFdAliasUnwrap = false;

  msUnknownCand = buildUnknownInitCandidateSet( ...
    caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
    dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
  caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
    viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    fdRange, fdRateRange, optVerbose, dynMsUnknownOpt, false, [], msUnknownCand);

  centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown);

  blockResult(iBlock).blockLen = blockLen;
  blockResult(iBlock).pilotWave = pilotWave;
  blockResult(iBlock).viewMs = viewMs;
  blockResult(iBlock).centerRegistry = centerRegistry;

  if scanCfg.runKnown
    blockResult(iBlock).cpKnown = localRunBlockSet(viewMs, pilotWave, waveInfo.sampleRate, ...
      carrierFreq, fdRange, fdRateRange, centerRegistry, scanCfg, truth, blockLen, 'known');
  end
  if scanCfg.runUnknown
    blockResult(iBlock).cpUnknown = localRunBlockSet(viewMs, pilotWave, waveInfo.sampleRate, ...
      carrierFreq, fdRange, fdRateRange, centerRegistry, scanCfg, truth, blockLen, 'unknown');
  end
end

scanResult.blockResult = blockResult;

%% Report and plot
if scanCfg.runKnown
  disp('========== CP-K block-length summary ==========');
  disp(localBuildBlockSummaryTable(blockResult, 'known'));
end
if scanCfg.runUnknown
  disp('========== CP-U block-length summary ==========');
  disp(localBuildBlockSummaryTable(blockResult, 'unknown'));
end

if localGetFieldOrDefault(scanCfg, 'showToothTable', true)
  if scanCfg.runKnown
    disp('========== CP-K block-length alias-tooth table ==========');
    disp(localBuildMergedToothTable(blockResult, 'known'));
  end
  if scanCfg.runUnknown
    disp('========== CP-U block-length alias-tooth table ==========');
    disp(localBuildMergedToothTable(blockResult, 'unknown'));
  end
end

if localGetFieldOrDefault(scanCfg, 'showDeltaFigure', true)
  if scanCfg.runKnown
    localPlotBlockDeltaFigure(blockResult, 'known', 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotBlockDeltaFigure(blockResult, 'unknown', 'CP-U');
  end
end

if localGetFieldOrDefault(scanCfg, 'showFoldedFigure', true)
  if scanCfg.runKnown
    localPlotBlockFoldedFigure(blockResult, 'known', 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotBlockFoldedFigure(blockResult, 'unknown', 'CP-U');
  end
end

if localGetFieldOrDefault(scanCfg, 'showPeakFigure', true)
  if scanCfg.runKnown
    localPlotBlockPeakFigure(blockResult, 'known', 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotBlockPeakFigure(blockResult, 'unknown', 'CP-U');
  end
end

%% Save
saveName = sprintf('%s_%s.mat', scanCfg.saveTag, datestr(now, 'yyyymmdd-HHMMSS'));
save(saveName, 'scanResult', '-v7.3');
fprintf('Saved block-length scan result to %s\n', saveName);

function modeResult = localRunBlockSet(viewMs, pilotWave, sampleRate, carrierFreq, ...
    fdRange, fdRateRange, centerRegistry, scanCfg, truth, blockLen, fdRateMode)
%LOCALRUNBLOCKSET Run one block-length scan set for known/unknown fdRate.

centerNameList = localResolveCenterNameList(scanCfg, fdRateMode);
modelOpt = localBuildModelOpt(fdRateMode, truth);
[modelPhase, ~, ~, modelOpt] = buildDoaDopplerMfModel( ...
  viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
  viewMs.doaGrid, fdRange, fdRateRange, modelOpt);

modeResult = struct();
modeResult.blockLen = blockLen;
modeResult.fdRateMode = string(fdRateMode);
modeResult.modelOpt = modelOpt;
modeResult.aliasStepHz = localResolveAliasStep(modelPhase);
modeResult.centerNameList = centerNameList;
modeResult.centerScan = repmat(struct( ...
  'centerName', "", ...
  'basePoint', struct(), ...
  'lineInfo', struct()), numel(centerNameList), 1);

iRow = 0;
toothCell = cell(numel(centerNameList), 1);
for iCenter = 1:numel(centerNameList)
  centerName = centerNameList(iCenter);
  basePoint = localResolveCenterPoint(centerRegistry, fdRateMode, centerName);
  halfWidthHz = localResolveHalfWidth(scanCfg, modeResult.aliasStepHz);
  fdRefGrid = linspace(basePoint.fdRef - halfWidthHz, basePoint.fdRef + halfWidthHz, scanCfg.numGrid);
  lineInfo = localScanFdRefLine(modelPhase, basePoint, fdRefGrid, scanCfg);

  modeResult.centerScan(iCenter).centerName = centerName;
  modeResult.centerScan(iCenter).basePoint = basePoint;
  modeResult.centerScan(iCenter).lineInfo = lineInfo;

  iRow = iRow + 1;
  toothCell{iRow} = localBuildAliasToothRow(blockLen, fdRateMode, centerName, lineInfo, scanCfg.numAliasSide);
end
modeResult.toothTable = vertcat(toothCell{1:iRow});
end

function modelOpt = localBuildModelOpt(fdRateMode, truth)
%LOCALBUILDMODELOPT Build one compact MF modelOpt for block-length scans.

modelOpt = struct();
modelOpt.useLogObjective = true;
modelOpt.initFdCount = 81;
modelOpt.useAccessMask = false;
modelOpt.phaseMode = 'continuous';
modelOpt.steeringMode = 'framewise';
modelOpt.steeringRefFrameIdx = [];
modelOpt.enableFdAliasUnwrap = false;
modelOpt.debugEnable = false;
modelOpt.debugStoreEvalTrace = false;
modelOpt.debugMaxEvalTrace = 1;
modelOpt.initDoaParam = [];
modelOpt.initDoaHalfWidth = [];

switch lower(char(fdRateMode))
  case 'known'
    modelOpt.fdRateMode = 'known';
    modelOpt.fdRateKnown = truth.fdRateFit;
  case 'unknown'
    modelOpt.fdRateMode = 'unknown';
  otherwise
    error('doaDopplerDynBlockLengthScan:InvalidFdRateMode', ...
      'Unsupported fdRateMode: %s', string(fdRateMode));
end
end

function lineInfo = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanCfg)
%LOCALSCANFDREFLINE Evaluate one fdRef comb line.

numGrid = numel(fdRefGrid);
objVec = nan(numGrid, 1);
residualVec = nan(numGrid, 1);
fdRefEvalVec = nan(numGrid, 1);
fdRateEvalVec = nan(numGrid, 1);
cohSatMat = nan(numGrid, modelBase.numSat);
objSatMat = nan(numGrid, modelBase.numSat);

probeTemplate = struct('doaParam', basePoint.latlon(:), ...
  'fdRef', basePoint.fdRef, 'fdRate', basePoint.fdRate, ...
  'obj', NaN, 'residualNorm', NaN);

useParfor = localResolveParforFlag(scanCfg, numGrid);
if useParfor
  parfor iGrid = 1:numGrid
    [objVec(iGrid), residualVec(iGrid), fdRefEvalVec(iGrid), fdRateEvalVec(iGrid), ...
      cohSatMat(iGrid, :), objSatMat(iGrid, :)] = ...
      localEvaluateScanPoint(modelBase, basePoint, fdRefGrid(iGrid), probeTemplate);
  end
else
  for iGrid = 1:numGrid
    [objVec(iGrid), residualVec(iGrid), fdRefEvalVec(iGrid), fdRateEvalVec(iGrid), ...
      cohSatMat(iGrid, :), objSatMat(iGrid, :)] = ...
      localEvaluateScanPoint(modelBase, basePoint, fdRefGrid(iGrid), probeTemplate);
  end
end

[minObj, minIdx] = min(objVec);
aliasStepHz = localResolveAliasStep(modelBase);
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

function [objVal, residualVal, fdRefEvalVal, fdRateEvalVal, cohSatRow, objSatRow] = ...
    localEvaluateScanPoint(modelBase, basePoint, fdRefVal, probeTemplate)
%LOCALEVALUATESCANPOINT Evaluate one scan point through the shared helper.

pointCur = basePoint;
pointCur.fdRef = fdRefVal;
probeEval = evalDoaDopplerMfProbePoint(modelBase, pointCur, ...
  struct('probeTemplate', probeTemplate, 'debugEnable', false));

objVal = localGetFieldOrDefault(probeEval, 'obj', NaN);
residualVal = localGetFieldOrDefault(probeEval, 'residualNorm', NaN);
fdRefEvalVal = localGetFieldOrDefault(probeEval, 'fdRef', NaN);
fdRateEvalVal = localGetFieldOrDefault(probeEval, 'fdRate', NaN);

cohSatRow = nan(1, modelBase.numSat);
objSatRow = nan(1, modelBase.numSat);
cohSat = reshape(localGetFieldOrDefault(probeEval, 'coherenceSat', []), 1, []);
objSat = reshape(localGetFieldOrDefault(probeEval, 'objectiveSat', []), 1, []);
numSatCopy = min(modelBase.numSat, numel(cohSat));
cohSatRow(1, 1:numSatCopy) = cohSat(1:numSatCopy);
numSatCopy = min(modelBase.numSat, numel(objSat));
objSatRow(1, 1:numSatCopy) = objSat(1:numSatCopy);
end

function toothTable = localBuildAliasToothRow(blockLen, fdRateMode, centerName, lineInfo, numAliasSide)
%LOCALBUILDALIASTOOTHROW Build one compact tooth summary row set.

aliasIdxVec = (-numAliasSide:numAliasSide).';
numTooth = numel(aliasIdxVec);
fdToothHz = nan(numTooth, 1);
deltaObjTooth = nan(numTooth, 1);
deltaFdTooth = nan(numTooth, 1);
for iTooth = 1:numTooth
  targetFd = lineInfo.basePoint.fdRef + aliasIdxVec(iTooth) * lineInfo.aliasStepHz;
  idx = localFindNearestGridIndex(lineInfo.fdRefGrid, targetFd);
  fdToothHz(iTooth) = lineInfo.fdRefGrid(idx);
  deltaObjTooth(iTooth) = lineInfo.objVec(idx) - lineInfo.minObj;
  deltaFdTooth(iTooth) = lineInfo.fdRefGrid(idx) - lineInfo.basePoint.fdRef;
end

truthGap = localGetAliasGap(deltaObjTooth, aliasIdxVec, 0);
aliasGap1 = localGetAliasGap(deltaObjTooth, aliasIdxVec, 1);
aliasGap2 = localGetAliasGap(deltaObjTooth, aliasIdxVec, 2);

toothTable = table(repmat(blockLen, numTooth, 1), repmat(string(fdRateMode), numTooth, 1), ...
  repmat(string(centerName), numTooth, 1), aliasIdxVec, fdToothHz, deltaFdTooth, deltaObjTooth, ...
  repmat(truthGap, numTooth, 1), repmat(aliasGap1, numTooth, 1), repmat(aliasGap2, numTooth, 1), ...
  'VariableNames', {'blockLen', 'fdRateMode', 'centerName', 'aliasIndex', ...
  'fdRefSampleHz', 'deltaFdRefHz', 'deltaObj', 'truthGap', 'aliasGap1', 'aliasGap2'});
end

function gapVal = localGetAliasGap(deltaObjTooth, aliasIdxVec, absAliasIdx)
%LOCALGETALIASGAP Extract one alias-tooth gap by absolute index.

mask = abs(aliasIdxVec) == absAliasIdx;
if ~any(mask)
  gapVal = NaN;
else
  gapVal = min(deltaObjTooth(mask));
end
end

function centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown)
%LOCALBUILDCENTERREGISTRY Build named centers shared by all scans.

centerRegistry = struct();
centerRegistry.truth = localBuildPointStruct(truth.latlonTrueDeg(:), truth.fdRefFit, truth.fdRateFit);

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

function centerNameList = localResolveCenterNameList(scanCfg, fdRateMode)
%LOCALRESOLVECENTERNAMELIST Resolve center names for known/unknown scans.

switch lower(char(fdRateMode))
  case 'known'
    centerNameList = reshape(string(scanCfg.centerListKnown), [], 1);
  case 'unknown'
    centerNameList = reshape(string(scanCfg.centerListUnknown), [], 1);
  otherwise
    error('doaDopplerDynBlockLengthScan:InvalidFdRateMode', ...
      'Unsupported fdRateMode: %s', string(fdRateMode));
end
end

function centerPoint = localResolveCenterPoint(centerRegistry, fdRateMode, centerName)
%LOCALRESOLVECENTERPOINT Resolve one named center point.

switch lower(char(centerName))
  case 'truth'
    centerPoint = centerRegistry.truth;
  case 'staticseed'
    centerPoint = centerRegistry.staticSeed;
  case 'cpknownseed'
    centerPoint = centerRegistry.cpKnownSeed;
  case 'finalestimate'
    if strcmpi(fdRateMode, 'known')
      centerPoint = centerRegistry.finalEstimateKnown;
    else
      centerPoint = centerRegistry.finalEstimateUnknown;
    end
  otherwise
    error('doaDopplerDynBlockLengthScan:UnknownCenter', ...
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

function summaryTable = localBuildBlockSummaryTable(blockResult, fdRateMode)
%LOCALBUILDBLOCKSUMMARYTABLE Build one compact summary across block lengths.

numBlock = numel(blockResult);
blockLen = nan(numBlock, 1);
truthGap = nan(numBlock, 1);
aliasGap1 = nan(numBlock, 1);
aliasGap2 = nan(numBlock, 1);
centerDeltaTruth = nan(numBlock, 1);
centerDeltaFinal = nan(numBlock, 1);
minDeltaTruth = nan(numBlock, 1);
minDeltaFinal = nan(numBlock, 1);

for iBlock = 1:numBlock
  blockLen(iBlock) = blockResult(iBlock).blockLen;
  modeResult = localSelectModeResult(blockResult(iBlock), fdRateMode);
  truthLine = localSelectCenterLine(modeResult, 'truth');
  finalLine = localSelectCenterLine(modeResult, 'finalEstimate');
  truthGap(iBlock) = localGetAliasGapFromLine(truthLine, 0);
  aliasGap1(iBlock) = localGetAliasGapFromLine(truthLine, 1);
  aliasGap2(iBlock) = localGetAliasGapFromLine(truthLine, 2);
  centerDeltaTruth(iBlock) = truthLine.centerDeltaObj;
  centerDeltaFinal(iBlock) = finalLine.centerDeltaObj;
  minDeltaTruth(iBlock) = truthLine.minDeltaFdRef;
  minDeltaFinal(iBlock) = finalLine.minDeltaFdRef;
end

summaryTable = table(blockLen, truthGap, aliasGap1, aliasGap2, centerDeltaTruth, ...
  centerDeltaFinal, minDeltaTruth, minDeltaFinal, ...
  'VariableNames', {'blockLen', 'truthGap', 'aliasGap1', 'aliasGap2', ...
  'centerDeltaTruth', 'centerDeltaFinal', 'minDeltaTruthHz', 'minDeltaFinalHz'});
end

function toothTable = localBuildMergedToothTable(blockResult, fdRateMode)
%LOCALBUILDMERGEDTOOTHTABLE Merge all tooth tables across block lengths.

cellTable = cell(numel(blockResult), 1);
for iBlock = 1:numel(blockResult)
  modeResult = localSelectModeResult(blockResult(iBlock), fdRateMode);
  cellTable{iBlock} = modeResult.toothTable;
end
toothTable = vertcat(cellTable{:});
end

function localPlotBlockDeltaFigure(blockResult, fdRateMode, modeLabel)
%LOCALPLOTBLOCKDELTAFIGURE Plot delta-objective curves vs block length.

centerNameList = localSelectModeResult(blockResult(1), fdRateMode).centerNameList;
numCenter = numel(centerNameList);
figure('Name', sprintf('%s block scan: delta objective', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  for iBlock = 1:numel(blockResult)
    modeResult = localSelectModeResult(blockResult(iBlock), fdRateMode);
    lineInfo = modeResult.centerScan(iCenter).lineInfo;
    plot(lineInfo.deltaFdRef, lineInfo.deltaObjVec, 'LineWidth', 1.2, ...
      'DisplayName', sprintf('N_p=%d', blockResult(iBlock).blockLen));
  end
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom');
  localOverlayAliasLines(localSelectModeResult(blockResult(1), fdRateMode).aliasStepHz, ...
    localSelectModeResult(blockResult(1), fdRateMode).centerScan(iCenter).lineInfo.deltaFdRef);
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('\Delta objective');
  title(sprintf('%s | %s', modeLabel, centerNameList(iCenter)), 'Interpreter', 'none');
  legend('Location', 'best');
  hold off;
end
end

function localPlotBlockFoldedFigure(blockResult, fdRateMode, modeLabel)
%LOCALPLOTBLOCKFOLDEDFIGURE Plot folded modulo-one-period views.

centerNameList = localSelectModeResult(blockResult(1), fdRateMode).centerNameList;
numCenter = numel(centerNameList);
figure('Name', sprintf('%s block scan: folded', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  for iBlock = 1:numel(blockResult)
    modeResult = localSelectModeResult(blockResult(iBlock), fdRateMode);
    lineInfo = modeResult.centerScan(iCenter).lineInfo;
    [xFold, sortIdx] = sort(lineInfo.foldedOffset);
    yFold = lineInfo.deltaObjVec(sortIdx);
    plot(xFold, yFold, '.', 'DisplayName', sprintf('N_p=%d', blockResult(iBlock).blockLen));
  end
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom');
  grid on;
  xlabel(sprintf('folded \\Delta fdRef (Hz), period = %.6f', ...
    localSelectModeResult(blockResult(1), fdRateMode).aliasStepHz));
  ylabel('\\Delta objective');
  title(sprintf('%s | %s | folded modulo one alias period', modeLabel, centerNameList(iCenter)), ...
    'Interpreter', 'none');
  legend('Location', 'best');
  hold off;
end
end

function localPlotBlockPeakFigure(blockResult, fdRateMode, modeLabel)
%LOCALPLOTBLOCKPEAKFIGURE Plot reciprocal-peak curves vs block length.

centerNameList = localSelectModeResult(blockResult(1), fdRateMode).centerNameList;
numCenter = numel(centerNameList);
figure('Name', sprintf('%s block scan: peak', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  for iBlock = 1:numel(blockResult)
    modeResult = localSelectModeResult(blockResult(iBlock), fdRateMode);
    lineInfo = modeResult.centerScan(iCenter).lineInfo;
    plot(lineInfo.deltaFdRef, lineInfo.reciprocalPeak, 'LineWidth', 1.2, ...
      'DisplayName', sprintf('N_p=%d', blockResult(iBlock).blockLen));
  end
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom');
  localOverlayAliasLines(localSelectModeResult(blockResult(1), fdRateMode).aliasStepHz, ...
    localSelectModeResult(blockResult(1), fdRateMode).centerScan(iCenter).lineInfo.deltaFdRef);
  grid on;
  xlabel('\\Delta fdRef (Hz)');
  ylabel('normalized reciprocal peak');
  title(sprintf('%s | %s | reciprocal peak', modeLabel, centerNameList(iCenter)), ...
    'Interpreter', 'none');
  legend('Location', 'best');
  hold off;
end
end

function modeResult = localSelectModeResult(blockEntry, fdRateMode)
%LOCALSELECTMODERESULT Select known or unknown mode result from one block entry.

switch lower(char(fdRateMode))
  case 'known'
    modeResult = blockEntry.cpKnown;
  case 'unknown'
    modeResult = blockEntry.cpUnknown;
  otherwise
    error('doaDopplerDynBlockLengthScan:InvalidFdRateMode', ...
      'Unsupported fdRateMode: %s', string(fdRateMode));
end
end

function lineInfo = localSelectCenterLine(modeResult, centerName)
%LOCALSELECTCENTERLINE Select one center line by name.

centerName = string(centerName);
for iCenter = 1:numel(modeResult.centerScan)
  if strcmpi(modeResult.centerScan(iCenter).centerName, centerName)
    lineInfo = modeResult.centerScan(iCenter).lineInfo;
    return;
  end
end
error('doaDopplerDynBlockLengthScan:CenterNotFound', ...
  'Center %s not found.', centerName);
end

function gapVal = localGetAliasGapFromLine(lineInfo, absAliasIdx)
%LOCALGETALIASGAPFROMLINE Extract one alias gap from lineInfo.

if ~(isscalar(lineInfo.aliasStepHz) && isfinite(lineInfo.aliasStepHz) && lineInfo.aliasStepHz > 0)
  gapVal = NaN;
  return;
end
targetFd = lineInfo.basePoint.fdRef + absAliasIdx * lineInfo.aliasStepHz;
idxPos = localFindNearestGridIndex(lineInfo.fdRefGrid, targetFd);
idxNeg = localFindNearestGridIndex(lineInfo.fdRefGrid, lineInfo.basePoint.fdRef - absAliasIdx * lineInfo.aliasStepHz);
if absAliasIdx == 0
  gapVal = lineInfo.objVec(idxPos) - lineInfo.minObj;
else
  gapVal = min(lineInfo.objVec([idxPos; idxNeg])) - lineInfo.minObj;
end
end



function numSample = localGetSignalLength(sig)
%LOCALGETSIGNALLENGTH Return waveform length along the sample dimension.
if isvector(sig)
  numSample = numel(sig);
  return;
end
[numRow, numCol] = size(sig);
numSample = max(numRow, numCol);
end

function blockLenList = localBuildDefaultBlockLenList(numPilotSampleFull)
%LOCALBUILDDEFAULTBLOCKLENLIST Build a fallback scan list from waveform length.
ratioList = [1/8; 1/4; 1/2; 3/4; 1];
blockLenList = round(ratioList * numPilotSampleFull);
blockLenList = unique(max(1, min(numPilotSampleFull, blockLenList)), 'stable');
end

function pilotWaveTrim = localTrimPilotWave(pilotWave, blockLen)
%LOCALTRIMPILOTWAVE Trim one pilot waveform along its sample dimension.

if isempty(pilotWave)
  pilotWaveTrim = pilotWave;
  return;
end

numRow = size(pilotWave, 1);
numCol = size(pilotWave, 2);
if numRow >= numCol
  lenUse = min(blockLen, numRow);
  pilotWaveTrim = pilotWave(1:lenUse, :, :, :, :, :);
else
  lenUse = min(blockLen, numCol);
  pilotWaveTrim = pilotWave(:, 1:lenUse, :, :, :, :);
end
end

function rxSigCellTrim = localTrimRxSigCell(rxSigCell, blockLen)
%LOCALTRIMRXSIGCELL Trim each frame signal to the first blockLen samples.

rxSigCellTrim = rxSigCell;
for iFrame = 1:numel(rxSigCell)
  frameSig = rxSigCell{iFrame};
  if iscell(frameSig)
    for iSat = 1:numel(frameSig)
      frameSig{iSat} = localTrimOneSig(frameSig{iSat}, blockLen);
    end
    rxSigCellTrim{iFrame} = frameSig;
  else
    rxSigCellTrim{iFrame} = localTrimOneSig(frameSig, blockLen);
  end
end
end

function sigTrim = localTrimOneSig(sig, blockLen)
%LOCALTRIMONESIG Trim one signal matrix along the sample dimension.

if isempty(sig)
  sigTrim = sig;
  return;
end
numRow = size(sig, 1);
lenUse = min(blockLen, numRow);
sigTrim = sig(1:lenUse, :, :, :, :, :);
end

function localEnsureParpool(scanCfg)
%LOCALENSUREPARPOOL Start parpool on demand for scan loops.

useParfor = logical(localGetFieldOrDefault(scanCfg, 'useParfor', false));
autoStartParpool = logical(localGetFieldOrDefault(scanCfg, 'autoStartParpool', false));
if ~(useParfor && autoStartParpool)
  return;
end
if isempty(gcp('nocreate'))
  parpool;
end
end

function useParfor = localResolveParforFlag(scanCfg, numGrid)
%LOCALRESOLVEPARFORFLAG Resolve whether one line scan should use parfor.

useParfor = logical(localGetFieldOrDefault(scanCfg, 'useParfor', false));
minGridForParfor = localGetFieldOrDefault(scanCfg, 'minGridForParfor', inf);
if ~(isscalar(minGridForParfor) && isnumeric(minGridForParfor))
  minGridForParfor = inf;
end
useParfor = useParfor && (numGrid >= minGridForParfor);
if useParfor
  useParfor = ~isempty(gcp('nocreate'));
end
end

function halfWidthHz = localResolveHalfWidth(scanCfg, aliasStepHz)
%LOCALRESOLVEHALFWIDTH Resolve scan half-width around one center.

halfWidthHz = localGetFieldOrDefault(scanCfg, 'scanHalfWidthHz', []);
if isempty(halfWidthHz)
  numAliasSide = localGetFieldOrDefault(scanCfg, 'numAliasSide', 2);
  halfWidthHz = max(1, numAliasSide) * aliasStepHz;
end
end

function reciprocalPeak = localBuildReciprocalPeak(objVec, minObj)
%LOCALBUILDRECIPROCALPEAK Build one reciprocal-peak curve.

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

function idx = localFindNearestGridIndex(gridVec, targetVal)
%LOCALFINDNEARESTGRIDINDEX Find the nearest index in one grid.

[~, idx] = min(abs(gridVec(:) - targetVal));
end

function val = localResolveDoaParam(estResult, defaultVal)
%LOCALRESOLVEDOAPARAM Resolve one DoA parameter vector from estResult.

val = localGetFieldOrDefault(estResult, 'doaParamEst', defaultVal);
val = reshape(val, [], 1);
end

function val = localResolveScalarField(s, fieldName, defaultVal)
%LOCALRESOLVESCALARFIELD Resolve one scalar field with fallback.

val = localGetFieldOrDefault(s, fieldName, defaultVal);
if isempty(val)
  val = defaultVal;
end
end

function localOverlayAliasLines(aliasStepHz, deltaFdRef)
%LOCALOVERLAYALIASLINES Overlay ±k/Tf reference lines.

if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  return;
end
spanMax = max(abs(deltaFdRef(:)));
numSide = max(1, ceil(spanMax / aliasStepHz));
for k = -numSide:numSide
  xline(k * aliasStepHz, ':');
end
end

function val = localGetFieldOrDefault(s, fieldName, defaultVal)
%LOCALGETFIELDORDEFAULT Get one struct field with fallback.

if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
