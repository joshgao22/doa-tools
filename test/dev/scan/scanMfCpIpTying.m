%SCANMFCPIPTYING Scan how CP/IP phase tying changes fdRef comb structure.
% Dev scan script. Edit the configuration section directly before running.
% Data storage saves only scanData through saveExpSnapshot; plots can be
% regenerated from scanData after loading.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfCpIpTying";
saveSnapshot = true;
notifyTelegramEnable = true;
checkpointEnable = false;

numUsr = 1;
numFrame = 20;
frameIntvlSec = 1 / 750;
refFrameIdx = ceil(numFrame / 2);

sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 128;
carrierFreq = 11.7e9;
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

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.runKnown = true;
scanConfig.runUnknown = true;
scanConfig.phaseList = ["continuous"; "relaxed"; "independent"];
scanConfig.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanConfig.centerListUnknown = ["truth"; "cpKnownSeed"; "finalEstimate"];
scanConfig.numAliasSide = 2;
scanConfig.numGrid = 801;
scanConfig.scanHalfWidthHz = [];
scanConfig.numFrame = numFrame;
scanConfig.frameCountList = numFrame;
scanConfig.frameIntvlSec = frameIntvlSec;
scanConfig.snrDb = snrDb;
scanConfig.numTask = scanConfig.numGrid * numel(scanConfig.phaseList) * ...
  (scanConfig.runKnown * numel(scanConfig.centerListKnown) + ...
   scanConfig.runUnknown * numel(scanConfig.centerListUnknown));

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
checkpointDir = "";
scanJustRan = false;
scanConfig.runKey = runKey;
printMfScanHeader(char(scanName), scanConfig, checkpointDir);

try
%% Build context and scan tasks
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);

if sceneSeq.numSat ~= 2
  error('doaDopplerDynCpIpTyingScan:InvalidNumSat', ...
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

%% Run scan batch
scanResult = struct();
scanResult.truth = truth;
scanResult.scanConfig = scanConfig;
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

if scanConfig.runKnown
  scanResult.cpKnown = localRunPhaseCompareSet(viewMs, pilotWave, waveInfo.sampleRate, ...
    carrierFreq, fdRange, fdRateRange, centerRegistry, scanConfig, truth, 'known');
end
if scanConfig.runUnknown
  scanResult.cpUnknown = localRunPhaseCompareSet(viewMs, pilotWave, waveInfo.sampleRate, ...
    carrierFreq, fdRange, fdRateRange, centerRegistry, scanConfig, truth, 'unknown');
end

%% Data storage
scanData = scanResult;
scanData.scanName = string(scanName);
scanData.runKey = string(runKey);
scanData.config = scanConfig;
scanData.utcRun = datetime('now', 'TimeZone', 'local');
scanData.elapsedSec = toc(runTic);
scanData = finalizeMfScanResult(scanData, checkpointDir);
scanData.snapshotFile = "";
if scanConfig.saveSnapshot
  saveOpt = struct('includeVars', {{'scanData'}}, ...
    'extraMeta', struct('scanName', char(scanName)), ...
    'verbose', true);
  scanData.snapshotFile = string(saveExpSnapshot(char(scanName), saveOpt));
end
scanJustRan = true;
catch ME
  fprintf('Scan failed while building CP/IP tying data.\n');
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Summary output and plotting
scanResult = scanData;
scanConfig = scanData.config;

if isfield(scanResult, 'cpKnown')
  printMfScanSection('CP-K phase-tying summary');
  localDisplayPhaseSummary('CP-K', scanResult.cpKnown);
end
if isfield(scanResult, 'cpUnknown')
  printMfScanSection('CP-U phase-tying summary');
  localDisplayPhaseSummary('CP-U', scanResult.cpUnknown);
end

if isfield(scanResult, 'cpKnown')
  printMfScanSection('CP-K alias-tooth table', scanResult.cpKnown.toothTable);
end
if isfield(scanResult, 'cpUnknown')
  printMfScanSection('CP-U alias-tooth table', scanResult.cpUnknown.toothTable);
end

if isfield(scanResult, 'cpKnown')
  localPlotPhaseDeltaFigure(scanResult.cpKnown, 'CP-K');
end
if isfield(scanResult, 'cpUnknown')
  localPlotPhaseDeltaFigure(scanResult.cpUnknown, 'CP-U');
end
if isfield(scanResult, 'cpKnown')
  localPlotPhaseFoldedFigure(scanResult.cpKnown, 'CP-K');
end
if isfield(scanResult, 'cpUnknown')
  localPlotPhaseFoldedFigure(scanResult.cpUnknown, 'CP-U');
end

if exist('scanJustRan', 'var') && scanJustRan
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "DONE", ...
    'config', scanConfig, ...
    'snapshotFile', scanData.snapshotFile, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "CP/IP phase-tying scan completed."; ...
      "Detailed tables remain in scanData and scan results docs."]));
end

%% Local helpers

function modeResult = localRunPhaseCompareSet(viewMs, pilotWave, sampleRate, carrierFreq, ...
    fdRange, fdRateRange, centerRegistry, scanConfig, truth, fdRateMode)
%LOCALRUNPHASECOMPARESET Run one fdRate-mode scan across phase modes.

phaseList = reshape(string(scanConfig.phaseList), [], 1);
centerNameList = localResolveCenterNameList(scanConfig, fdRateMode);
modeResult = struct();
modeResult.fdRateMode = string(fdRateMode);
modeResult.phaseList = phaseList;
modeResult.centerNameList = centerNameList;
modeResult.phaseScan = repmat(struct( ...
  'phaseMode', "", ...
  'modelOpt', struct(), ...
  'aliasStepHz', NaN, ...
  'centerScan', []), numel(phaseList), 1);

toothCell = cell(numel(phaseList) * numel(centerNameList), 1);
iRow = 0;
totalEvalCount = numel(phaseList) * numel(centerNameList) * scanConfig.numGrid;
fprintf('[CpIpTying] %s: %d phase modes x %d centers x %d fdRef points\n', ...
  upper(char(fdRateMode)), numel(phaseList), numel(centerNameList), scanConfig.numGrid);
progressTracker = localCreateScanProgressTracker(totalEvalCount);

try
for iPhase = 1:numel(phaseList)
  phaseMode = phaseList(iPhase);
  modelOpt = localBuildPhaseModelOpt(phaseMode, fdRateMode, truth);
  [modelPhase, ~, ~, modelOpt] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
    viewMs.doaGrid, fdRange, fdRateRange, modelOpt);

  phaseScan = struct();
  phaseScan.phaseMode = phaseMode;
  phaseScan.modelOpt = modelOpt;
  phaseScan.aliasStepHz = localResolveAliasStep(modelPhase);
  phaseScan.centerScan = repmat(struct( ...
    'centerName', "", ...
    'basePoint', struct(), ...
    'lineInfo', struct()), numel(centerNameList), 1);

  for iCenter = 1:numel(centerNameList)
    centerName = centerNameList(iCenter);
    basePoint = localResolveCenterPoint(centerRegistry, fdRateMode, centerName);
    halfWidthHz = localResolveHalfWidth(scanConfig, phaseScan.aliasStepHz);
    fdRefGrid = linspace(basePoint.fdRef - halfWidthHz, basePoint.fdRef + halfWidthHz, scanConfig.numGrid);
    lineInfo = localScanFdRefLine(modelPhase, basePoint, fdRefGrid, scanConfig, progressTracker);

    phaseScan.centerScan(iCenter).centerName = centerName;
    phaseScan.centerScan(iCenter).basePoint = basePoint;
    phaseScan.centerScan(iCenter).lineInfo = lineInfo;

    iRow = iRow + 1;
    toothCell{iRow} = localBuildAliasToothRow(fdRateMode, phaseMode, centerName, lineInfo, scanConfig.numAliasSide);
  end

  modeResult.phaseScan(iPhase) = phaseScan;
end
catch ME
  localCloseScanProgressTracker(progressTracker);
  rethrow(ME);
end

localCloseScanProgressTracker(progressTracker);
modeResult.summaryTable = localBuildPhaseSummaryTable(modeResult);
modeResult.toothTable = vertcat(toothCell{1:iRow});
end

function summaryTable = localBuildPhaseSummaryTable(modeResult)
%LOCALBUILDPHASESUMMARYTABLE Build one compact phase / center summary table.

numRow = 0;
for iPhase = 1:numel(modeResult.phaseScan)
  numRow = numRow + numel(modeResult.phaseScan(iPhase).centerScan);
end

fdRateMode = strings(numRow, 1);
phaseMode = strings(numRow, 1);
centerName = strings(numRow, 1);
aliasStepHz = nan(numRow, 1);
minDeltaFdRefHz = nan(numRow, 1);
minAliasIndex = nan(numRow, 1);
centerDeltaObj = nan(numRow, 1);
minFdRefHz = nan(numRow, 1);
centerFdRefHz = nan(numRow, 1);

rowIdx = 0;
for iPhase = 1:numel(modeResult.phaseScan)
  phaseScan = modeResult.phaseScan(iPhase);
  for iCenter = 1:numel(phaseScan.centerScan)
    rowIdx = rowIdx + 1;
    lineInfo = phaseScan.centerScan(iCenter).lineInfo;
    fdRateMode(rowIdx) = string(modeResult.fdRateMode);
    phaseMode(rowIdx) = string(phaseScan.phaseMode);
    centerName(rowIdx) = string(phaseScan.centerScan(iCenter).centerName);
    aliasStepHz(rowIdx) = lineInfo.aliasStepHz;
    minDeltaFdRefHz(rowIdx) = lineInfo.minDeltaFdRef;
    minAliasIndex(rowIdx) = lineInfo.minAliasIndex;
    centerDeltaObj(rowIdx) = lineInfo.centerDeltaObj;
    minFdRefHz(rowIdx) = lineInfo.minFdRef;
    centerFdRefHz(rowIdx) = lineInfo.basePoint.fdRef;
  end
end

summaryTable = table(fdRateMode, phaseMode, centerName, aliasStepHz, ...
  minDeltaFdRefHz, minAliasIndex, centerDeltaObj, minFdRefHz, centerFdRefHz);
end

function progressTracker = localCreateScanProgressTracker(totalCount)
%LOCALCREATESCANPROGRESSTRACKER Create a best-effort progress tracker.

progressTracker = struct('enabled', false, 'queue', []);
if totalCount <= 0
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  progressTracker.enabled = true;
  try
    progressTracker.queue = parallel.pool.DataQueue;
    afterEach(progressTracker.queue, @(~) progressbar('advance'));
  catch
    progressTracker.queue = [];
  end
catch
  progressTracker.enabled = false;
  progressTracker.queue = [];
end
end

function localAdvanceScanProgress(progressTracker, isParforWorker)
%LOCALADVANCESCANPROGRESS Advance progress from serial or parfor execution.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
if isParforWorker
  if isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
    send(progressTracker.queue, 1);
  end
else
  try
    progressbar('advance');
  catch
  end
end
end

function localCloseScanProgressTracker(progressTracker)
%LOCALCLOSESCANPROGRESSTRACKER Close a best-effort progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function modelOpt = localBuildPhaseModelOpt(phaseMode, fdRateMode, truth)
%LOCALBUILDPHASEMODELOPT Build one compact MF modelOpt for probe scans.

modelOpt = struct();
modelOpt.useLogObjective = true;
modelOpt.initFdCount = 81;
modelOpt.useAccessMask = false;
modelOpt.phaseMode = char(phaseMode);
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
    error('doaDopplerDynCpIpTyingScan:InvalidFdRateMode', ...
      'Unsupported fdRateMode: %s', string(fdRateMode));
end
end

function lineInfo = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanConfig, progressTracker)
%LOCALSCANFDREFLINE Evaluate one fdRef line for one model / center pair.

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

useParfor = localCanUseParfor(numGrid);
if useParfor
  parfor iGrid = 1:numGrid
    [objVec(iGrid), residualVec(iGrid), fdRefEvalVec(iGrid), fdRateEvalVec(iGrid), ...
      cohSatMat(iGrid, :), objSatMat(iGrid, :)] = ...
      localEvaluateScanPoint(modelBase, basePoint, fdRefGrid(iGrid), probeTemplate);
    localAdvanceScanProgress(progressTracker, true);
  end
else
  for iGrid = 1:numGrid
    [objVec(iGrid), residualVec(iGrid), fdRefEvalVec(iGrid), fdRateEvalVec(iGrid), ...
      cohSatMat(iGrid, :), objSatMat(iGrid, :)] = ...
      localEvaluateScanPoint(modelBase, basePoint, fdRefGrid(iGrid), probeTemplate);
    localAdvanceScanProgress(progressTracker, false);
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

function toothTable = localBuildAliasToothRow(fdRateMode, phaseMode, centerName, lineInfo, numAliasSide)
%LOCALBUILDALIASTOOTHROW Build one compact alias-tooth summary table.

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

phaseMode = string(phaseMode);
centerName = string(centerName);
fdRateMode = string(fdRateMode);
toothTable = table(repmat(fdRateMode, numTooth, 1), repmat(phaseMode, numTooth, 1), ...
  repmat(centerName, numTooth, 1), aliasIdxVec, fdToothHz, deltaFdTooth, deltaObjTooth, ...
  repmat(truthGap, numTooth, 1), repmat(aliasGap1, numTooth, 1), repmat(aliasGap2, numTooth, 1), ...
  'VariableNames', {'fdRateMode', 'phaseMode', 'centerName', 'aliasIndex', ...
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
%LOCALBUILDCENTERREGISTRY Build named scan centers.

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

function centerNameList = localResolveCenterNameList(scanConfig, fdRateMode)
%LOCALRESOLVECENTERNAMELIST Resolve the center list for one fdRate mode.

switch lower(char(fdRateMode))
  case 'known'
    centerNameList = reshape(string(scanConfig.centerListKnown), [], 1);
  case 'unknown'
    centerNameList = reshape(string(scanConfig.centerListUnknown), [], 1);
  otherwise
    error('doaDopplerDynCpIpTyingScan:InvalidFdRateMode', ...
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
    error('doaDopplerDynCpIpTyingScan:UnknownCenter', ...
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

function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build compact scan-specific notification lines.

metricLineList = strings(0, 1);
config = scanData.config;
metricLineList(end + 1, 1) = sprintf('• Phase modes: <code>%s</code>', ...
  strjoin(string(config.phaseList(:)).', ', '));
metricLineList(end + 1, 1) = sprintf('• fdRef grid: <code>%d</code>', config.numGrid);
if isfield(scanData, 'cpKnown') && isfield(scanData.cpKnown, 'summaryTable')
  row = localFindMetricRow(scanData.cpKnown.summaryTable, "continuous", "truth");
  if ~isempty(row)
    metricLineList(end + 1, 1) = sprintf('• CP-K truth min tooth: <code>%+.0f</code>', row.minAliasIndex);
  end
end
if isfield(scanData, 'cpUnknown') && isfield(scanData.cpUnknown, 'summaryTable')
  row = localFindMetricRow(scanData.cpUnknown.summaryTable, "continuous", "cpKnownSeed");
  if ~isempty(row)
    metricLineList(end + 1, 1) = sprintf('• CP-U cpKnownSeed min tooth: <code>%+.0f</code>', row.minAliasIndex);
  end
end
end

function row = localFindMetricRow(summaryTable, phaseMode, centerName)
%LOCALFINDMETRICROW Return one summary row for notification metrics.

idx = summaryTable.phaseMode == string(phaseMode) & summaryTable.centerName == string(centerName);
if any(idx)
  row = summaryTable(find(idx, 1, 'first'), :);
else
  row = [];
end
end

function localDisplayPhaseSummary(modeLabel, modeResult)
%LOCALDISPLAYPHASESUMMARY Print a compact one-line summary per phase / center.

for iPhase = 1:numel(modeResult.phaseScan)
  phaseScan = modeResult.phaseScan(iPhase);
  for iCenter = 1:numel(phaseScan.centerScan)
    lineInfo = phaseScan.centerScan(iCenter).lineInfo;
    fprintf('%s | %s | %s | minDeltaFd = %.6f Hz | aliasIdx = %d | centerDeltaObj = %.6g\n', ...
      modeLabel, phaseScan.phaseMode, phaseScan.centerScan(iCenter).centerName, ...
      lineInfo.minDeltaFdRef, lineInfo.minAliasIndex, lineInfo.centerDeltaObj);
  end
end
end

function localPlotPhaseDeltaFigure(modeResult, modeLabel)
%LOCALPLOTPHASEDELTAFIGURE Plot delta-objective scans for all phase modes.

numCenter = numel(modeResult.centerNameList);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s tying scan: delta objective', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  for iPhase = 1:numel(modeResult.phaseScan)
    phaseScan = modeResult.phaseScan(iPhase);
    lineInfo = phaseScan.centerScan(iCenter).lineInfo;
    plot(lineInfo.deltaFdRef, lineInfo.deltaObjVec, 'LineWidth', 1.35, ...
      'DisplayName', char(phaseScan.phaseMode));
  end
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
  localOverlayAliasLines(modeResult.phaseScan(1).aliasStepHz, ...
    modeResult.phaseScan(1).centerScan(iCenter).lineInfo.deltaFdRef);
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('\Delta objective');
  title(sprintf('%s | %s', modeLabel, modeResult.centerNameList(iCenter)), 'Interpreter', 'none');
  legend('Location', 'best');
  hold off;
end
end

function localPlotPhaseFoldedFigure(modeResult, modeLabel)
%LOCALPLOTPHASEFOLDEDFIGURE Plot folded modulo-one-period views.

numCenter = numel(modeResult.centerNameList);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s tying scan: folded', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  subplot(numCenter, 1, iCenter);
  hold on;
  for iPhase = 1:numel(modeResult.phaseScan)
    phaseScan = modeResult.phaseScan(iPhase);
    lineInfo = phaseScan.centerScan(iCenter).lineInfo;
    [xFold, sortIdx] = sort(lineInfo.foldedOffset);
    yFold = lineInfo.deltaObjVec(sortIdx);
    plot(xFold, yFold, '.', 'DisplayName', char(phaseScan.phaseMode));
  end
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
  grid on;
  xlabel(sprintf('folded \\Delta fdRef (Hz), period = %.6f', modeResult.phaseScan(1).aliasStepHz));
  ylabel('\Delta objective');
  title(sprintf('%s | %s | folded modulo one alias period', ...
    modeLabel, modeResult.centerNameList(iCenter)), 'Interpreter', 'none');
  legend('Location', 'best');
  hold off;
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
  xVal = k * aliasStepHz;
  xline(xVal, ':', 'HandleVisibility', 'off');
end
end

function halfWidthHz = localResolveHalfWidth(scanConfig, aliasStepHz)
%LOCALRESOLVEHALFWIDTH Resolve one scan half-width around the center.

halfWidthHz = localGetFieldOrDefault(scanConfig, 'scanHalfWidthHz', []);
if isempty(halfWidthHz)
  numAliasSide = localGetFieldOrDefault(scanConfig, 'numAliasSide', 2);
  halfWidthHz = max(1, numAliasSide) * aliasStepHz;
end
end

function reciprocalPeak = localBuildReciprocalPeak(objVec, minObj)
%LOCALBUILDRECIPROCALPEAK Build one reciprocal peak curve.

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
%LOCALRESOLVEDOAPARAM Resolve one DoA parameter vector from an estResult.

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

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Return true when an outer scan loop can use parfor.

useParfor = false;
if numTask <= 1
  return;
end
try
  useParfor = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && ...
    license('test', 'Distrib_Computing_Toolbox');
catch
  useParfor = false;
end
end

