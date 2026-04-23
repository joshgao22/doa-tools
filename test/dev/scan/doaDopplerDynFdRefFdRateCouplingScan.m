% Development script for diagnosing fdRef-fdRate coupling in MF
% DoA-Doppler estimation. The script keeps the paper-aligned observation
% block definition and evaluates two-dimensional objective maps around a few
% representative centers. It is intended to answer whether unknown fdRate
% mainly creates a slanted valley that connects adjacent 1/T_f comb teeth.
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
rng(271);

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
scanCfg.numGridFdRef = 241;
scanCfg.numGridFdRate = 121;
scanCfg.numAliasSide = 2;
scanCfg.fdRefHalfWidthHz = [];
scanCfg.fdRateHalfWidthHzPerSec = 600;
scanCfg.useParfor = true;
scanCfg.autoStartParpool = true;
scanCfg.minGridForParfor = 1000;
scanCfg.showDeltaFigure = true;
scanCfg.showPeakFigure = true;
scanCfg.showMeshFigure = true;
scanCfg.showSummaryTable = true;
scanCfg.saveTag = "doaDopplerDynFdRefFdRateCouplingScan";
scanCfg.logEnable = true;
scanCfg.logPrefix = "[FdRefFdRateCoupling]";

%% Select satellites and build scene
localLog(scanCfg, 'Building scene and selecting satellites ...');
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
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

%% Pilot waveform and snapshots aligned to N_p = 2112
localLog(scanCfg, 'Generating pilot waveform and snapshots ...');
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
numPilotSample = localGetSignalLength(pilotWave);
if numPilotSample ~= standardBlockLen
  error('doaDopplerDynFdRefFdRateCouplingScan:UnexpectedPilotLength', ...
    ['Generated pilot length is %d samples, but this script is aligned to the paper ' ...
     'definition N_p = %d samples per frame.'], numPilotSample, standardBlockLen);
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

%% Auto-start parpool if requested
localEnsureParpool(scanCfg);
localLog(scanCfg, 'Preparing estimation views and dynamic centers ...');

%% Build views and dynamic centers
sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = selectSatScene(sceneRef, otherSatIdxLocal);
rxSigRef = rxSigCell{sceneSeq.refFrameIdx};
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

%% Run 2-D scans
localLog(scanCfg, 'Starting fdRef-fdRate coupling scans ...');
scanResult = struct();
scanResult.truth = truth;
scanResult.scanCfg = scanCfg;
scanResult.frameInfo = struct( ...
  'frameIntvlSec', frameIntvlSec, ...
  'timeOffsetSec', sceneSeq.timeOffsetSec(:), ...
  'numFrame', sceneSeq.numFrame, ...
  'numPilotSample', numPilotSample, ...
  'refFrameIdx', sceneSeq.refFrameIdx, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'refStateSource', string(refState.source));
scanResult.centerRegistry = centerRegistry;

if scanCfg.runKnown
  scanResult.cpKnown = localRunCouplingSet(viewMs, pilotWave, waveInfo.sampleRate, ...
    carrierFreq, fdRange, fdRateRange, centerRegistry, scanCfg, truth, 'known');
end
if scanCfg.runUnknown
  scanResult.cpUnknown = localRunCouplingSet(viewMs, pilotWave, waveInfo.sampleRate, ...
    carrierFreq, fdRange, fdRateRange, centerRegistry, scanCfg, truth, 'unknown');
end

%% Report and plot
if localGetFieldOrDefault(scanCfg, 'showSummaryTable', true)
  if scanCfg.runKnown
    disp('========== CP-K fdRef-fdRate coupling summary ==========');
    disp(localBuildCouplingSummaryTable(scanResult.cpKnown));
  end
  if scanCfg.runUnknown
    disp('========== CP-U fdRef-fdRate coupling summary ==========');
    disp(localBuildCouplingSummaryTable(scanResult.cpUnknown));
  end
end

if localGetFieldOrDefault(scanCfg, 'showDeltaFigure', true)
  if scanCfg.runKnown
    localPlotCouplingDeltaFigure(scanResult.cpKnown, 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotCouplingDeltaFigure(scanResult.cpUnknown, 'CP-U');
  end
end

if localGetFieldOrDefault(scanCfg, 'showPeakFigure', true)
  if scanCfg.runKnown
    localPlotCouplingPeakFigure(scanResult.cpKnown, 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotCouplingPeakFigure(scanResult.cpUnknown, 'CP-U');
  end
end

if localGetFieldOrDefault(scanCfg, 'showMeshFigure', true)
  if scanCfg.runKnown
    localPlotCouplingMeshFigure(scanResult.cpKnown, 'CP-K');
  end
  if scanCfg.runUnknown
    localPlotCouplingMeshFigure(scanResult.cpUnknown, 'CP-U');
  end
end

%% Save
saveName = sprintf('%s_%s.mat', scanCfg.saveTag, datestr(now, 'yyyymmdd-HHMMSS'));
save(saveName, 'scanResult', '-v7.3');
fprintf('Saved fdRef-fdRate coupling scan result to %s\n', saveName);

function modeResult = localRunCouplingSet(viewMs, pilotWave, sampleRate, carrierFreq, ...
    fdRange, fdRateRange, centerRegistry, scanCfg, truth, fdRateMode)
%LOCALRUNCOUPLINGSET Run one 2-D coupling scan set for known/unknown fdRate.

centerNameList = localResolveCenterNameList(scanCfg, fdRateMode);
localLog(scanCfg, sprintf('Building MF model for fdRateMode=%s ...', string(fdRateMode)));
modelOpt = localBuildModelOpt(fdRateMode, truth);
[modelPhase, ~, ~, modelOpt] = buildDoaDopplerMfModel( ...
  viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
  viewMs.doaGrid, fdRange, fdRateRange, modelOpt);
localLog(scanCfg, sprintf('Built MF model for fdRateMode=%s. aliasStepHz = %.6g', ...
  string(fdRateMode), localResolveAliasStep(modelPhase)));

modeResult = struct();
modeResult.fdRateMode = string(fdRateMode);
modeResult.modelOpt = modelOpt;
modeResult.aliasStepHz = localResolveAliasStep(modelPhase);
modeResult.centerNameList = centerNameList;
modeResult.centerScan = repmat(struct( ...
  'centerName', "", ...
  'basePoint', struct(), ...
  'mapInfo', struct()), numel(centerNameList), 1);

for iCenter = 1:numel(centerNameList)
  centerName = centerNameList(iCenter);
  basePoint = localResolveCenterPoint(centerRegistry, fdRateMode, centerName);
  fdRefHalfWidthHz = localResolveFdRefHalfWidth(scanCfg, modeResult.aliasStepHz);
  fdRateHalfWidth = localGetFieldOrDefault(scanCfg, 'fdRateHalfWidthHzPerSec', 600);
  fdRefGrid = linspace(basePoint.fdRef - fdRefHalfWidthHz, basePoint.fdRef + fdRefHalfWidthHz, ...
    scanCfg.numGridFdRef);
  fdRateGrid = linspace(basePoint.fdRate - fdRateHalfWidth, basePoint.fdRate + fdRateHalfWidth, ...
    scanCfg.numGridFdRate);
  localLog(scanCfg, sprintf(['Running %s center=%s | grid=%d x %d | fdRef span=[%.3f, %.3f] Hz | ' ...
    'fdRate span=[%.3f, %.3f] Hz/s'], string(fdRateMode), centerName, ...
    numel(fdRefGrid), numel(fdRateGrid), fdRefGrid(1), fdRefGrid(end), fdRateGrid(1), fdRateGrid(end)));
  tCenter = tic;
  mapInfo = localScanFdRefFdRateMap(modelPhase, basePoint, fdRefGrid, fdRateGrid, scanCfg);
  localLog(scanCfg, sprintf(['Finished %s center=%s in %.2f s | minΔfdRef=%.3f Hz | ' ...
    'minΔfdRate=%.3f Hz/s | centerΔobj=%.6g'], string(fdRateMode), centerName, toc(tCenter), ...
    mapInfo.minDeltaFdRef, mapInfo.minDeltaFdRate, mapInfo.centerDeltaObj));

  modeResult.centerScan(iCenter).centerName = centerName;
  modeResult.centerScan(iCenter).basePoint = basePoint;
  modeResult.centerScan(iCenter).mapInfo = mapInfo;
end
end

function mapInfo = localScanFdRefFdRateMap(modelBase, basePoint, fdRefGrid, fdRateGrid, scanCfg)
%LOCALSCANFDREFFDRATEMAP Evaluate one two-dimensional fdRef-fdRate map.

numFdRef = numel(fdRefGrid);
numFdRate = numel(fdRateGrid);
numPoint = numFdRef * numFdRate;

objVec = nan(numPoint, 1);
residualVec = nan(numPoint, 1);
fdRefEvalVec = nan(numPoint, 1);
fdRateEvalVec = nan(numPoint, 1);
cohSatMat = nan(numPoint, modelBase.numSat);

useParfor = localResolveParforFlag(scanCfg, numPoint);
localLog(scanCfg, sprintf('  Evaluating %d map points (%s) ...', numPoint, ternary(useParfor, 'parfor', 'for')));
if useParfor
  parfor iPoint = 1:numPoint
    [iRate, iRef] = ind2sub([numFdRate, numFdRef], iPoint);
    [objVec(iPoint), residualVec(iPoint), fdRefEvalVec(iPoint), fdRateEvalVec(iPoint), ...
      cohSatMat(iPoint, :)] = localEvaluateMapPoint(modelBase, basePoint, ...
      fdRefGrid(iRef), fdRateGrid(iRate));
  end
else
  for iPoint = 1:numPoint
    [iRate, iRef] = ind2sub([numFdRate, numFdRef], iPoint);
    [objVec(iPoint), residualVec(iPoint), fdRefEvalVec(iPoint), fdRateEvalVec(iPoint), ...
      cohSatMat(iPoint, :)] = localEvaluateMapPoint(modelBase, basePoint, ...
      fdRefGrid(iRef), fdRateGrid(iRate));
  end
end

objMat = reshape(objVec, [numFdRate, numFdRef]);
residualMat = reshape(residualVec, [numFdRate, numFdRef]);
fdRefEvalMat = reshape(fdRefEvalVec, [numFdRate, numFdRef]);
fdRateEvalMat = reshape(fdRateEvalVec, [numFdRate, numFdRef]);
cohSatCube = reshape(cohSatMat, [numFdRate, numFdRef, modelBase.numSat]);
[minObj, minIdx] = min(objVec);
[iRateMin, iRefMin] = ind2sub([numFdRate, numFdRef], minIdx);
centerIdxRef = localFindNearestGridIndex(fdRefGrid, basePoint.fdRef);
centerIdxRate = localFindNearestGridIndex(fdRateGrid, basePoint.fdRate);

mapInfo = struct();
mapInfo.aliasStepHz = localResolveAliasStep(modelBase);
mapInfo.basePoint = basePoint;
mapInfo.fdRefGrid = fdRefGrid(:);
mapInfo.fdRateGrid = fdRateGrid(:);
mapInfo.objMat = objMat;
mapInfo.deltaObjMat = objMat - minObj;
mapInfo.reciprocalPeakMat = localBuildReciprocalPeak(objMat, minObj);
mapInfo.residualMat = residualMat;
mapInfo.fdRefEvalMat = fdRefEvalMat;
mapInfo.fdRateEvalMat = fdRateEvalMat;
mapInfo.cohSatCube = cohSatCube;
mapInfo.minObj = minObj;
mapInfo.minFdRef = fdRefGrid(iRefMin);
mapInfo.minFdRate = fdRateGrid(iRateMin);
mapInfo.minDeltaFdRef = fdRefGrid(iRefMin) - basePoint.fdRef;
mapInfo.minDeltaFdRate = fdRateGrid(iRateMin) - basePoint.fdRate;
mapInfo.centerObj = objMat(centerIdxRate, centerIdxRef);
mapInfo.centerDeltaObj = mapInfo.centerObj - minObj;
mapInfo.centerIdxRef = centerIdxRef;
mapInfo.centerIdxRate = centerIdxRate;
mapInfo.bestAliasIndex = round(mapInfo.minDeltaFdRef / mapInfo.aliasStepHz);
mapInfo.rowAtCenterRate = objMat(centerIdxRate, :);
mapInfo.colAtCenterRef = objMat(:, centerIdxRef);
localLog(scanCfg, sprintf('  Map ready | bestAliasIndex=%d | minObj=%.6g', mapInfo.bestAliasIndex, mapInfo.minObj));
end

function [objVal, residualVal, fdRefEvalVal, fdRateEvalVal, cohSatRow] = ...
    localEvaluateMapPoint(modelBase, basePoint, fdRefVal, fdRateVal)
%LOCALEVALUATEMAPPOINT Evaluate one map point through the shared helper.

pointCur = basePoint;
pointCur.fdRef = fdRefVal;
pointCur.fdRate = fdRateVal;
probeEval = evalDoaDopplerMfProbePoint(modelBase, pointCur);

objVal = localGetFieldOrDefault(probeEval, 'obj', NaN);
residualVal = localGetFieldOrDefault(probeEval, 'residualNorm', NaN);
fdRefEvalVal = localGetFieldOrDefault(probeEval, 'fdRef', NaN);
fdRateEvalVal = localGetFieldOrDefault(probeEval, 'fdRate', NaN);
cohSatRow = nan(1, modelBase.numSat);
cohSat = reshape(localGetFieldOrDefault(probeEval, 'coherenceSat', []), 1, []);
numSatCopy = min(modelBase.numSat, numel(cohSat));
cohSatRow(1, 1:numSatCopy) = cohSat(1:numSatCopy);
end

function summaryTable = localBuildCouplingSummaryTable(modeResult)
%LOCALBUILDCOUPLINGSUMMARYTABLE Build one compact summary table.

numCenter = numel(modeResult.centerScan);
centerName = strings(numCenter, 1);
centerDeltaObj = nan(numCenter, 1);
minDeltaFdRef = nan(numCenter, 1);
minDeltaFdRate = nan(numCenter, 1);
aliasGap1 = nan(numCenter, 1);
aliasGap2 = nan(numCenter, 1);
fdRateValleyGap = nan(numCenter, 1);

for iCenter = 1:numCenter
  mapInfo = modeResult.centerScan(iCenter).mapInfo;
  centerName(iCenter) = modeResult.centerScan(iCenter).centerName;
  centerDeltaObj(iCenter) = mapInfo.centerDeltaObj;
  minDeltaFdRef(iCenter) = mapInfo.minDeltaFdRef;
  minDeltaFdRate(iCenter) = mapInfo.minDeltaFdRate;
  aliasGap1(iCenter) = localGetAliasGapFromRow(mapInfo, 1);
  aliasGap2(iCenter) = localGetAliasGapFromRow(mapInfo, 2);
  fdRateValleyGap(iCenter) = min(mapInfo.colAtCenterRef) - mapInfo.minObj;
end

summaryTable = table(centerName, centerDeltaObj, minDeltaFdRef, minDeltaFdRate, ...
  aliasGap1, aliasGap2, fdRateValleyGap, ...
  'VariableNames', {'centerName', 'centerDeltaObj', 'minDeltaFdRefHz', ...
  'minDeltaFdRateHzPerSec', 'aliasGap1', 'aliasGap2', 'fdRateValleyGap'});
end

function gapVal = localGetAliasGapFromRow(mapInfo, absAliasIdx)
%LOCALGETALIASGAPFROMROW Get one alias gap from the row through center fdRate.

if ~(isscalar(mapInfo.aliasStepHz) && isfinite(mapInfo.aliasStepHz) && mapInfo.aliasStepHz > 0)
  gapVal = NaN;
  return;
end
targetPos = mapInfo.basePoint.fdRef + absAliasIdx * mapInfo.aliasStepHz;
targetNeg = mapInfo.basePoint.fdRef - absAliasIdx * mapInfo.aliasStepHz;
idxPos = localFindNearestGridIndex(mapInfo.fdRefGrid, targetPos);
idxNeg = localFindNearestGridIndex(mapInfo.fdRefGrid, targetNeg);
if absAliasIdx == 0
  gapVal = mapInfo.rowAtCenterRate(idxPos) - mapInfo.minObj;
else
  gapVal = min(mapInfo.rowAtCenterRate([idxPos; idxNeg])) - mapInfo.minObj;
end
end

function localPlotCouplingDeltaFigure(modeResult, modeLabel)
%LOCALPLOTCOUPLINGDELTAFIGURE Plot 2-D delta-objective maps.

numCenter = numel(modeResult.centerScan);
figure('Name', sprintf('%s fdRef-fdRate coupling: delta objective', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  mapInfo = modeResult.centerScan(iCenter).mapInfo;
  imagesc(mapInfo.fdRefGrid, mapInfo.fdRateGrid, mapInfo.deltaObjMat);
  axis xy;
  hold on;
  plot(mapInfo.basePoint.fdRef, mapInfo.basePoint.fdRate, 'wo', 'MarkerSize', 8, 'LineWidth', 1.4);
  plot(mapInfo.minFdRef, mapInfo.minFdRate, 'wx', 'MarkerSize', 9, 'LineWidth', 1.8);
  grid on;
  xlabel('fdRef (Hz)');
  ylabel('fdRate (Hz/s)');
  title(sprintf('%s | %s | \Delta objective', modeLabel, modeResult.centerScan(iCenter).centerName), ...
    'Interpreter', 'none');
  colorbar;
  hold off;
end
end

function localPlotCouplingPeakFigure(modeResult, modeLabel)
%LOCALPLOTCOUPLINGPEAKFIGURE Plot reciprocal-peak maps.

numCenter = numel(modeResult.centerScan);
figure('Name', sprintf('%s fdRef-fdRate coupling: reciprocal peak', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  mapInfo = modeResult.centerScan(iCenter).mapInfo;
  imagesc(mapInfo.fdRefGrid, mapInfo.fdRateGrid, mapInfo.reciprocalPeakMat);
  axis xy;
  hold on;
  plot(mapInfo.basePoint.fdRef, mapInfo.basePoint.fdRate, 'wo', 'MarkerSize', 8, 'LineWidth', 1.4);
  plot(mapInfo.minFdRef, mapInfo.minFdRate, 'wx', 'MarkerSize', 9, 'LineWidth', 1.8);
  grid on;
  xlabel('fdRef (Hz)');
  ylabel('fdRate (Hz/s)');
  title(sprintf('%s | %s | reciprocal peak', modeLabel, modeResult.centerScan(iCenter).centerName), ...
    'Interpreter', 'none');
  colorbar;
  hold off;
end
end

function localPlotCouplingMeshFigure(modeResult, modeLabel)
%LOCALPLOTCOUPLINGMESHFIGURE Plot 3-D mesh views.

numCenter = numel(modeResult.centerScan);
figure('Name', sprintf('%s fdRef-fdRate coupling: mesh', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  mapInfo = modeResult.centerScan(iCenter).mapInfo;
  mesh(mapInfo.fdRefGrid, mapInfo.fdRateGrid, mapInfo.deltaObjMat);
  hold on;
  plot3(mapInfo.basePoint.fdRef, mapInfo.basePoint.fdRate, mapInfo.centerDeltaObj, ...
    'ko', 'MarkerSize', 7, 'LineWidth', 1.2);
  plot3(mapInfo.minFdRef, mapInfo.minFdRate, 0, ...
    'kx', 'MarkerSize', 8, 'LineWidth', 1.6);
  grid on;
  xlabel('fdRef (Hz)');
  ylabel('fdRate (Hz/s)');
  zlabel('\Delta objective');
  title(sprintf('%s | %s | mesh', modeLabel, modeResult.centerScan(iCenter).centerName), ...
    'Interpreter', 'none');
  hold off;
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
    error('doaDopplerDynFdRefFdRateCouplingScan:InvalidFdRateMode', ...
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
    error('doaDopplerDynFdRefFdRateCouplingScan:UnknownCenter', ...
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
pointInfo.doaParam = doaParam(1:2);
pointInfo.fdRef = fdRef;
pointInfo.fdRate = fdRate;
end

function modelOpt = localBuildModelOpt(fdRateMode, truth)
%LOCALBUILDMODELOPT Build one compact MF modelOpt for coupling scans.

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
    error('doaDopplerDynFdRefFdRateCouplingScan:InvalidFdRateMode', ...
      'Unsupported fdRateMode: %s', string(fdRateMode));
end
end

function halfWidthHz = localResolveFdRefHalfWidth(scanCfg, aliasStepHz)
%LOCALRESOLVEFDREFHALFWIDTH Resolve fdRef half-width around one center.

halfWidthHz = localGetFieldOrDefault(scanCfg, 'fdRefHalfWidthHz', []);
if isempty(halfWidthHz)
  numAliasSide = localGetFieldOrDefault(scanCfg, 'numAliasSide', 2);
  halfWidthHz = max(1, numAliasSide) * aliasStepHz;
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


function localLog(scanCfg, msg)
%LOCALLOG Lightweight progress logger for long scans.

if ~localGetFieldOrDefault(scanCfg, 'logEnable', true)
  return;
end
prefix = string(localGetFieldOrDefault(scanCfg, 'logPrefix', '[doaDopplerDynFdRefFdRateCouplingScan]'));
fprintf('%s %s\n', prefix, msg);
end

function out = ternary(cond, a, b)
%TERNARY Small helper for readable log strings.
if cond
  out = a;
else
  out = b;
end
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

function useParfor = localResolveParforFlag(scanCfg, numPoint)
%LOCALRESOLVEPARFORFLAG Resolve whether one map scan should use parfor.

useParfor = logical(localGetFieldOrDefault(scanCfg, 'useParfor', false));
minGridForParfor = localGetFieldOrDefault(scanCfg, 'minGridForParfor', inf);
if ~(isscalar(minGridForParfor) && isnumeric(minGridForParfor))
  minGridForParfor = inf;
end
useParfor = useParfor && (numPoint >= minGridForParfor);
if useParfor
  useParfor = ~isempty(gcp('nocreate'));
end
end

function reciprocalPeak = localBuildReciprocalPeak(objArr, minObj)
%LOCALBUILDRECIPROCALPEAK Build one reciprocal-peak surface.

deltaObj = objArr - minObj;
positiveDelta = deltaObj(isfinite(deltaObj) & deltaObj > 0);
if isempty(positiveDelta)
  floorVal = 1;
else
  floorVal = max(min(positiveDelta), 1e-12);
end
reciprocalPeak = 1 ./ max(deltaObj + floorVal, floorVal);
reciprocalPeak = reciprocalPeak / max(reciprocalPeak(:));
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

function idx = localFindNearestGridIndex(gridVec, value)
%LOCALFINDNEARESTGRIDINDEX Find the nearest grid index to one scalar value.

[~, idx] = min(abs(gridVec(:) - value));
end

function val = localResolveDoaParam(estResult, defaultDoa)
%LOCALRESOLVEDOAPARAM Resolve DoA parameter from an estimate struct.

val = localGetFieldOrDefault(estResult, 'doaParamEst', []);
if isempty(val)
  val = defaultDoa;
end
val = reshape(val, [], 1);
end

function val = localResolveScalarField(s, fieldName, defaultVal)
%LOCALRESOLVESCALARFIELD Resolve one scalar field with a default fallback.

val = localGetFieldOrDefault(s, fieldName, defaultVal);
if isempty(val)
  val = defaultVal;
end
end

function val = localGetFieldOrDefault(s, fieldName, defaultVal)
%LOCALGETFIELDORDEFAULT Return one struct field with default fallback.

if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
