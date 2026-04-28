%SCANMFBLOCKLENGTH Scan pilot block length effects on MF comb structure.
% This scan trims a common generated pilot/snapshot to several block lengths
% and reuses the CP evaluator to compare fdRef tooth separation under
% known-rate and unknown-rate centers. The default list keeps the 2112-sample
% baseline and adds 2x/4x longer blocks for long-block sensitivity.
% Data storage saves only lightweight scanData through saveExpSnapshot; plots
% can be regenerated from scanData after loading.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfBlockLength";

numUsr = 1;
numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = ceil(numFrame / 2);

sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
baseBlockLen = 2112;
standardBlockLen = 4 * baseBlockLen;
numSym = round(standardBlockLen / osf);
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

doBaseHalfWidth = [0.01; 0.01];

doaLocalHalfWidthKnown = [0.003; 0.003];
doaLocalHalfWidthUnknown = [0.002; 0.002];

scanConfig = struct();
scanConfig.saveSnapshot = false;
scanConfig.checkpointEnable = true;
scanConfig.checkpointResume = true;
scanConfig.checkpointCleanupOnSuccess = true;
scanConfig.runKnown = true;
scanConfig.runUnknown = true;
scanConfig.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanConfig.centerListUnknown = ["truth"; "cpKnownSeed"; "finalEstimate"];
scanConfig.blockLenList = unique(round([baseBlockLen / 4; baseBlockLen / 2; ...
  baseBlockLen; 2 * baseBlockLen; standardBlockLen]), 'stable');
scanConfig.numAliasSide = 2;
scanConfig.numGrid = 801;
scanConfig.scanHalfWidthHz = [];

runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
%% Build context and flow options
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

%% Generate pilot waveform and snapshots
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWaveFull, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

numPilotSampleFull = localGetSignalLength(pilotWaveFull);
if numPilotSampleFull ~= standardBlockLen
  error('doaDopplerDynBlockLengthScan:UnexpectedPilotLength', ...
    ['Generated pilot length is %d samples, but this script expected the generated max ' ...
     'pilot block length to be %d samples. Keep genPilotSymbol / genPilotWaveform ' ...
     'consistent with the long-block scan configuration before running the scan.'], ...
    numPilotSampleFull, standardBlockLen);
end

blockLenListReq = reshape(scanConfig.blockLenList, [], 1);
blockLenList = blockLenListReq(blockLenListReq > 0 & blockLenListReq <= standardBlockLen);
blockLenList = unique(round(blockLenList), 'stable');
if isempty(blockLenList)
  error('doaDopplerDynBlockLengthScan:EmptyBlockLenList', ...
    'Requested blockLenList is empty after clipping to the standard block length.');
end
scanConfig.blockLenList = blockLenList;
localPrintScanHeader(scanName, scanConfig, numFrame, frameIntvlSec, snrDb, ...
  waveInfo.sampleRate, symbolRate, baseBlockLen, standardBlockLen, numSym);


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

%% Run scan batch
scanResult = struct();
scanResult.truth = truth;
scanResult.scanConfig = scanConfig;
scanResult.frameInfo = struct( ...
  'frameIntvlSec', frameIntvlSec, ...
  'timeOffsetSec', sceneSeq.timeOffsetSec(:), ...
  'numFrame', sceneSeq.numFrame, ...
  'baseBlockLen', baseBlockLen, ...
  'standardBlockLen', standardBlockLen, ...
  'numPilotSampleFull', numPilotSampleFull, ...
  'refFrameIdx', sceneSeq.refFrameIdx, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'refStateSource', string(refState.source));
scanResult.blockLenList = blockLenList;

blockTaskGrid = localBuildBlockTaskGrid(blockLenList);
blockSharedData = localBuildBlockSharedData(sceneSeq, sceneRef, refSatIdxLocal, ...
  otherSatIdxLocal, otherSatIdxGlobal, pilotWaveFull, rxSigCellFull, gridSize, ...
  searchRange, E, wavelength, carrierFreq, waveInfo.sampleRate, fdRange, ...
  fdRateRange, truth, optVerbose, weightSweepAlpha, doBaseHalfWidth, ...
  doaLocalHalfWidthKnown, doaLocalHalfWidthUnknown, scanConfig, baseBlockLen, ...
  numPilotSampleFull);

checkpointDir = "";
checkpointRunState = struct();
try
  if scanConfig.checkpointEnable
    checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, blockLenList, snrDb, ...
      waveInfo.sampleRate, symbolRate, carrierFreq);
    checkpointDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
    checkpointRunState = runPerfTaskGridWithCheckpoint(blockTaskGrid, blockSharedData, ...
      @localRunBlockCheckpointTask, checkpointOpt);
    blockResult = localCollectBlockTaskResults(checkpointRunState.resultCell, blockLenList);
  else
    blockResult = localRunBlockTaskList(blockTaskGrid, blockSharedData);
  end
catch ME
  localPrintCheckpointFailureHint(checkpointDir);
  rethrow(ME);
end

scanResult.blockResult = blockResult;
scanResult.summaryKnown = table();
scanResult.summaryUnknown = table();
scanResult.toothKnown = table();
scanResult.toothUnknown = table();
if scanConfig.runKnown
  scanResult.summaryKnown = localBuildBlockSummaryTable(blockResult, 'known');
  scanResult.toothKnown = localBuildMergedToothTable(blockResult, 'known');
end
if scanConfig.runUnknown
  scanResult.summaryUnknown = localBuildBlockSummaryTable(blockResult, 'unknown');
  scanResult.toothUnknown = localBuildMergedToothTable(blockResult, 'unknown');
end
scanResult.aggregateTable = localBuildAggregateTable(scanResult.summaryKnown, scanResult.summaryUnknown);
checkpointCleanupReport = struct([]);
if scanConfig.checkpointEnable
  if localGetFieldOrDefault(checkpointRunState, 'isComplete', false) && scanConfig.checkpointCleanupOnSuccess
    checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(checkpointRunState, ...
      struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
  end
end
scanResult.checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunState, ...
  checkpointDir, checkpointCleanupReport);
scanResult.checkpointCleanupReport = checkpointCleanupReport;

%% Data storage
scanData = scanResult;
scanData.scanName = scanName;
scanData.runKey = runKey;
scanData.utcRun = datetime('now', 'TimeZone', 'local');
scanData.snapshotFile = "";
if scanConfig.saveSnapshot
  saveOpt = struct('includeVars', {{'scanData'}}, ...
    'extraMeta', struct('scanName', char(scanName)), ...
    'verbose', true);
  scanData.snapshotFile = string(saveExpSnapshot(char(scanName), saveOpt));
end

%% Summary output and plotting
if ~exist('scanData', 'var') || ~isstruct(scanData)
  error('scanMfBlockLength:MissingScanData', ...
    'Run the scan batch sections or load a snapshot containing scanData.');
end

scanConfig = scanData.scanConfig;
blockResult = scanData.blockResult;

fprintf('\n========== Block-length aggregate summary ==========%s', newline);
disp(scanData.aggregateTable);
if isfield(scanData, 'checkpointSummaryTable') && istable(scanData.checkpointSummaryTable) && ...
    height(scanData.checkpointSummaryTable) > 0
  fprintf('\n========== Checkpoint summary ==========%s', newline);
  disp(scanData.checkpointSummaryTable);
end
if scanConfig.runKnown
  fprintf('\n========== CP-K block-length alias-tooth table ==========%s', newline);
  disp(scanData.toothKnown);
end
if scanConfig.runUnknown
  fprintf('\n========== CP-U block-length alias-tooth table ==========%s', newline);
  disp(scanData.toothUnknown);
end

if scanConfig.runKnown
  localPlotBlockDeltaFigure(blockResult, 'known', 'CP-K');
  localPlotBlockFoldedFigure(blockResult, 'known', 'CP-K');
  localPlotBlockPeakFigure(blockResult, 'known', 'CP-K');
end
if scanConfig.runUnknown
  localPlotBlockDeltaFigure(blockResult, 'unknown', 'CP-U');
  localPlotBlockFoldedFigure(blockResult, 'unknown', 'CP-U');
  localPlotBlockPeakFigure(blockResult, 'unknown', 'CP-U');
end


function taskGrid = localBuildBlockTaskGrid(blockLenList)
%LOCALBUILDBLOCKTASKGRID Build one independent checkpoint task per block length.

blockLenList = reshape(blockLenList, [], 1);
taskGrid = repmat(struct('taskIndex', NaN, 'blockLen', NaN), numel(blockLenList), 1);
for iTask = 1:numel(blockLenList)
  taskGrid(iTask).taskIndex = iTask;
  taskGrid(iTask).blockLen = blockLenList(iTask);
end
end

function sharedData = localBuildBlockSharedData(sceneSeq, sceneRef, refSatIdxLocal, ...
    otherSatIdxLocal, otherSatIdxGlobal, pilotWaveFull, rxSigCellFull, gridSize, ...
    searchRange, E, wavelength, carrierFreq, sampleRate, fdRange, fdRateRange, ...
    truth, optVerbose, weightSweepAlpha, doBaseHalfWidth, doaLocalHalfWidthKnown, ...
    doaLocalHalfWidthUnknown, scanConfig, baseBlockLen, numPilotSampleFull)
%LOCALBUILDBLOCKSHAREDDATA Pack read-only inputs for block-level tasks.

sharedData = struct();
sharedData.sceneSeq = sceneSeq;
sharedData.sceneRef = sceneRef;
sharedData.refSatIdxLocal = refSatIdxLocal;
sharedData.otherSatIdxLocal = otherSatIdxLocal;
sharedData.otherSatIdxGlobal = otherSatIdxGlobal;
sharedData.pilotWaveFull = pilotWaveFull;
sharedData.rxSigCellFull = rxSigCellFull;
sharedData.gridSize = gridSize;
sharedData.searchRange = searchRange;
sharedData.E = E;
sharedData.wavelength = wavelength;
sharedData.carrierFreq = carrierFreq;
sharedData.sampleRate = sampleRate;
sharedData.fdRange = fdRange;
sharedData.fdRateRange = fdRateRange;
sharedData.truth = truth;
sharedData.optVerbose = optVerbose;
sharedData.weightSweepAlpha = weightSweepAlpha;
sharedData.doBaseHalfWidth = doBaseHalfWidth;
sharedData.doaLocalHalfWidthKnown = doaLocalHalfWidthKnown;
sharedData.doaLocalHalfWidthUnknown = doaLocalHalfWidthUnknown;
sharedData.scanConfig = scanConfig;
sharedData.baseBlockLen = baseBlockLen;
sharedData.numPilotSampleFull = numPilotSampleFull;
end

function blockResult = localRunBlockTaskList(taskGrid, sharedData)
%LOCALRUNBLOCKTASKLIST Run block tasks without checkpoint files.

blockResult = repmat(localEmptyBlockResult(), numel(taskGrid), 1);
for iTask = 1:numel(taskGrid)
  blockResult(iTask) = localRunBlockCheckpointTask(taskGrid(iTask), sharedData);
end
end

function blockEntry = localRunBlockCheckpointTask(taskInfo, sharedData)
%LOCALRUNBLOCKCHECKPOINTTASK Run one checkpointable block-length task.

blockEntry = localRunOneBlockScan(taskInfo.blockLen, sharedData);
end

function blockResult = localCollectBlockTaskResults(resultCell, blockLenList)
%LOCALCOLLECTBLOCKTASKRESULTS Collect checkpoint task results in blockLen order.

blockResult = repmat(localEmptyBlockResult(), numel(blockLenList), 1);
for iTask = 1:numel(blockLenList)
  if iTask > numel(resultCell) || isempty(resultCell{iTask})
    error('scanMfBlockLength:MissingCheckpointTask', ...
      'Checkpoint task %d for blockLen=%d did not produce a result.', iTask, blockLenList(iTask));
  end
  blockEntry = resultCell{iTask};
  if isstruct(blockEntry) && isfield(blockEntry, 'blockResult')
    blockEntry = blockEntry.blockResult;
  end
  blockResult(iTask) = blockEntry;
end
end

function blockEntry = localEmptyBlockResult()
%LOCALEMPTYBLOCKRESULT Return the scalar block result template.

blockEntry = struct( ...
  'blockLen', NaN, ...
  'blockInfo', struct(), ...
  'centerRegistry', struct(), ...
  'cpKnown', struct(), ...
  'cpUnknown', struct());
end

function blockEntry = localRunOneBlockScan(blockLen, d)
%LOCALRUNONEBLOCKSCAN Run static seeds, dynamic centers, and fdRef lines for one block.

fprintf('Running block-length scan: %d samples\n', blockLen);

pilotWave = localTrimPilotWave(d.pilotWaveFull, blockLen);
rxSigCell = localTrimRxSigCell(d.rxSigCellFull, blockLen);
rxSigRef = rxSigCell{d.sceneSeq.refFrameIdx};

sceneSeqRefOnly = selectSatSceneSeq(d.sceneSeq, d.refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(rxSigCell, d.refSatIdxLocal, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  d.gridSize, d.searchRange, d.E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = selectSatScene(d.sceneRef, d.otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSigRef, d.otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  d.gridSize, d.searchRange, d.E);

viewMs = buildDoaDopplerEstView(d.sceneRef, rxSigCell, ...
  d.gridSize, d.searchRange, d.E, struct('sceneSeq', d.sceneSeq));

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;
doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, d.wavelength, pilotWave, ...
  d.carrierFreq, d.sampleRate, d.fdRange, d.truth, ...
  d.otherSatIdxGlobal, d.optVerbose, doaOnlyOpt, staticBaseOpt, ...
  d.weightSweepAlpha, d.doBaseHalfWidth);

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;
initParamStaticMs = buildDynamicInitParamFromCase(bestStaticMsCase, true, d.truth.fdRateFit);

dynBaseOpt = struct();
dynBaseOpt.useLogObjective = true;
dynBaseOpt.initFdCount = 81;
dynBaseOpt.useAccessMask = false;
dynBaseOpt.phaseMode = 'continuous';
dynBaseOpt.steeringMode = 'framewise';
dynBaseOpt.steeringRefFrameIdx = d.sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = false;
dynBaseOpt.debugStoreEvalTrace = false;
dynBaseOpt.debugMaxEvalTrace = 1;

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = d.doaLocalHalfWidthKnown;
dynMsKnownOpt.enableFdAliasUnwrap = false;

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  viewMs, d.truth, pilotWave, d.carrierFreq, d.sampleRate, ...
  d.fdRange, d.fdRateRange, d.optVerbose, dynMsKnownOpt, true, [], initParamStaticMs);

initParamMsUnknownCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, d.truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = d.doaLocalHalfWidthUnknown;
dynMsUnknownOpt.enableFdAliasUnwrap = false;

msUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, staticMsOpt.initDoaHalfWidth);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  viewMs, d.truth, pilotWave, d.carrierFreq, d.sampleRate, ...
  d.fdRange, d.fdRateRange, d.optVerbose, dynMsUnknownOpt, false, [], msUnknownCand);

centerRegistry = localBuildCenterRegistry(d.truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown);

blockEntry = localEmptyBlockResult();
blockEntry.blockLen = blockLen;
blockEntry.blockInfo = localBuildBlockInfo(blockLen, d.baseBlockLen, ...
  d.numPilotSampleFull, d.sampleRate, d.truth);
blockEntry.centerRegistry = centerRegistry;

if d.scanConfig.runKnown
  blockEntry.cpKnown = localRunBlockSet(viewMs, pilotWave, d.sampleRate, ...
    d.carrierFreq, d.fdRange, d.fdRateRange, centerRegistry, d.scanConfig, ...
    d.truth, blockLen, 'known');
end
if d.scanConfig.runUnknown
  blockEntry.cpUnknown = localRunBlockSet(viewMs, pilotWave, d.sampleRate, ...
    d.carrierFreq, d.fdRange, d.fdRateRange, centerRegistry, d.scanConfig, ...
    d.truth, blockLen, 'unknown');
end
end

function checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, blockLenList, snrDb, ...
    sampleRate, symbolRate, carrierFreq)
%LOCALBUILDCHECKPOINTOPT Build block-level checkpoint runner options.

checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, blockLenList, snrDb, ...
  sampleRate, symbolRate, carrierFreq);
checkpointOpt.outputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.useParfor = false;
checkpointOpt.resume = logical(scanConfig.checkpointResume);
checkpointOpt.meta = localBuildCheckpointMeta(scanConfig, blockLenList, snrDb, ...
  sampleRate, symbolRate, carrierFreq);
checkpointOpt.progressFcn = [];
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end

function runKey = localBuildCheckpointRunKey(scanConfig, blockLenList, snrDb, ...
    sampleRate, symbolRate, carrierFreq)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable key for interrupted block scans.

blockText = char(strjoin(compose('%d', reshape(blockLenList, 1, [])), 'x'));
runKey = sprintf('blocks%s_grid%d_snr%.2f_fs%.0f_sym%.0f_fc%.0f_known%d_unknown%d', ...
  blockText, scanConfig.numGrid, snrDb, sampleRate, symbolRate, carrierFreq, ...
  scanConfig.runKnown, scanConfig.runUnknown);
runKey = string(runKey);
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
runKey = replace(runKey, ' ', '');
end

function meta = localBuildCheckpointMeta(scanConfig, blockLenList, snrDb, ...
    sampleRate, symbolRate, carrierFreq)
%LOCALBUILDCHECKPOINTMETA Store the task-defining scan signature.

meta = struct();
meta.blockLenList = reshape(blockLenList, 1, []);
meta.snrDb = snrDb;
meta.sampleRate = sampleRate;
meta.symbolRate = symbolRate;
meta.carrierFreq = carrierFreq;
meta.runKnown = logical(scanConfig.runKnown);
meta.runUnknown = logical(scanConfig.runUnknown);
meta.centerListKnown = reshape(string(scanConfig.centerListKnown), 1, []);
meta.centerListUnknown = reshape(string(scanConfig.centerListUnknown), 1, []);
meta.numAliasSide = scanConfig.numAliasSide;
meta.numGrid = scanConfig.numGrid;
meta.scanHalfWidthHz = scanConfig.scanHalfWidthHz;
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runState, checkpointDir, cleanupReport)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Summarize block checkpoint resume and cleanup state.

if ~isstruct(runState) || ~isfield(runState, 'numTask')
  checkpointSummaryTable = table();
  return;
end
cleanupApplied = ~isempty(cleanupReport);
checkpointSummaryTable = table(string(checkpointDir), runState.numTask, ...
  runState.numDone, logical(runState.isComplete), cleanupApplied, ...
  'VariableNames', {'expectedRunDir', 'numTask', 'numDone', 'isComplete', 'cleanupApplied'});
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print preserved checkpoint path after failures.

checkpointDir = string(checkpointDir);
if strlength(checkpointDir) == 0
  fprintf('Scan failed. No checkpoint directory was created.\n');
  return;
end
fprintf('Scan failed. Completed checkpoint tasks were kept under tmp for resume.\n');
fprintf('  checkpoint dir: %s\n', char(checkpointDir));
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve the repository root from this scan script location.

scanDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scanDir)));
end

function localPrintScanHeader(scanName, scanConfig, numFrame, frameIntvlSec, snrDb, ...
    sampleRate, symbolRate, baseBlockLen, standardBlockLen, numSym)
%LOCALPRINTSCANHEADER Print compact scan configuration.

fprintf('Running %s ...\n', char(scanName));
fprintf('  frame count                     : %d\n', numFrame);
fprintf('  frame interval (ms)             : %.6f\n', frameIntvlSec * 1e3);
fprintf('  sample rate (MHz)               : %.3f\n', sampleRate / 1e6);
fprintf('  symbol rate (MHz)               : %.3f\n', symbolRate / 1e6);
fprintf('  baseline block length             : %d samples\n', baseBlockLen);
fprintf('  generated max block length         : %d samples (%.1fx baseline)\n', ...
  standardBlockLen, standardBlockLen / baseBlockLen);
fprintf('  generated pilot symbols            : %d\n', numSym);
fprintf('  snr (dB)                        : %.2f\n', snrDb);
fprintf('  block lengths (samples)         : %s\n', char(localFormatIntegerRow(scanConfig.blockLenList)));
fprintf('  run known / unknown             : %d / %d\n', scanConfig.runKnown, scanConfig.runUnknown);
fprintf('  fdRef scan grid count           : %d\n', scanConfig.numGrid);
fprintf('  alias side count                : %d\n', scanConfig.numAliasSide);
fprintf('  checkpoint enable               : %d (resume=%d, cleanupOnSuccess=%d)\n', ...
  scanConfig.checkpointEnable, scanConfig.checkpointResume, scanConfig.checkpointCleanupOnSuccess);
end

function blockInfo = localBuildBlockInfo(blockLen, baseBlockLen, fullBlockLen, sampleRate, truth)
%LOCALBUILDBLOCKINFO Build physical block-duration metrics for summary tables.

blockDurationSec = blockLen / sampleRate;
fdRateAbsHzPerSec = abs(localResolveScalarField(truth, 'fdRateFit', NaN));
blockInfo = struct();
blockInfo.blockLen = blockLen;
blockInfo.baseBlockLen = baseBlockLen;
blockInfo.fullBlockLen = fullBlockLen;
blockInfo.blockLenRatio = blockLen / fullBlockLen;
blockInfo.blockLenBaseRatio = blockLen / baseBlockLen;
blockInfo.blockDurationSec = blockDurationSec;
blockInfo.blockDurationUs = blockDurationSec * 1e6;
blockInfo.inBlockFdDriftHz = fdRateAbsHzPerSec * blockDurationSec;
blockInfo.inBlockQuadPhaseRad = pi * fdRateAbsHzPerSec * blockDurationSec^2;
end

function modeResult = localRunBlockSet(viewMs, pilotWave, sampleRate, carrierFreq, ...
    fdRange, fdRateRange, centerRegistry, scanConfig, truth, blockLen, fdRateMode)
%LOCALRUNBLOCKSET Run one block-length scan set for known/unknown fdRate.

centerNameList = localResolveCenterNameList(scanConfig, fdRateMode);
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
  halfWidthHz = localResolveHalfWidth(scanConfig, modeResult.aliasStepHz);
  fdRefGrid = linspace(basePoint.fdRef - halfWidthHz, basePoint.fdRef + halfWidthHz, scanConfig.numGrid);
  lineInfo = localScanFdRefLine(modelPhase, basePoint, fdRefGrid, scanConfig, ...
    blockLen, fdRateMode, centerName);

  modeResult.centerScan(iCenter).centerName = centerName;
  modeResult.centerScan(iCenter).basePoint = basePoint;
  modeResult.centerScan(iCenter).lineInfo = lineInfo;

  iRow = iRow + 1;
  toothCell{iRow} = localBuildAliasToothRow(blockLen, fdRateMode, centerName, lineInfo, scanConfig.numAliasSide);
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

function lineInfo = localScanFdRefLine(modelBase, basePoint, fdRefGrid, scanConfig, blockLen, fdRateMode, centerName)
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

progressLabel = sprintf('block=%d, mode=%s, center=%s', blockLen, char(fdRateMode), char(centerName));
progressTracker = localStartProgress(numGrid, progressLabel);
cleanupObj = onCleanup(@() localFinishProgress(progressTracker)); %#ok<NASGU>

useParfor = localCanUseParfor(numGrid);
if useParfor && progressTracker.isActive
  progressQueue = parallel.pool.DataQueue;
  afterEach(progressQueue, @(~) progressbar('advance'));
  parfor iGrid = 1:numGrid
    [objVec(iGrid), residualVec(iGrid), fdRefEvalVec(iGrid), fdRateEvalVec(iGrid), ...
      cohSatMat(iGrid, :), objSatMat(iGrid, :)] = ...
      localEvaluateScanPoint(modelBase, basePoint, fdRefGrid(iGrid), probeTemplate);
    send(progressQueue, 0);
  end
elseif useParfor
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
    localAdvanceProgress(progressTracker);
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

function centerNameList = localResolveCenterNameList(scanConfig, fdRateMode)
%LOCALRESOLVECENTERNAMELIST Resolve center names for known/unknown scans.

switch lower(char(fdRateMode))
  case 'known'
    centerNameList = reshape(string(scanConfig.centerListKnown), [], 1);
  case 'unknown'
    centerNameList = reshape(string(scanConfig.centerListUnknown), [], 1);
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
fdRateModeCol = repmat(string(fdRateMode), numBlock, 1);
blockLen = nan(numBlock, 1);
blockLenRatio = nan(numBlock, 1);
blockLenBaseRatio = nan(numBlock, 1);
blockDurationUs = nan(numBlock, 1);
inBlockFdDriftHz = nan(numBlock, 1);
inBlockQuadPhaseRad = nan(numBlock, 1);
truthGap = nan(numBlock, 1);
aliasGap1 = nan(numBlock, 1);
aliasGap2 = nan(numBlock, 1);
centerDeltaTruth = nan(numBlock, 1);
centerDeltaFinal = nan(numBlock, 1);
minDeltaTruth = nan(numBlock, 1);
minDeltaFinal = nan(numBlock, 1);

for iBlock = 1:numBlock
  blockInfo = localGetFieldOrDefault(blockResult(iBlock), 'blockInfo', struct());
  blockLen(iBlock) = blockResult(iBlock).blockLen;
  blockLenRatio(iBlock) = localGetFieldOrDefault(blockInfo, 'blockLenRatio', NaN);
  blockLenBaseRatio(iBlock) = localGetFieldOrDefault(blockInfo, 'blockLenBaseRatio', NaN);
  blockDurationUs(iBlock) = localGetFieldOrDefault(blockInfo, 'blockDurationUs', NaN);
  inBlockFdDriftHz(iBlock) = localGetFieldOrDefault(blockInfo, 'inBlockFdDriftHz', NaN);
  inBlockQuadPhaseRad(iBlock) = localGetFieldOrDefault(blockInfo, 'inBlockQuadPhaseRad', NaN);
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

summaryTable = table(fdRateModeCol, blockLen, blockLenRatio, blockLenBaseRatio, blockDurationUs, ...
  inBlockFdDriftHz, inBlockQuadPhaseRad, truthGap, aliasGap1, aliasGap2, ...
  centerDeltaTruth, centerDeltaFinal, minDeltaTruth, minDeltaFinal, ...
  'VariableNames', {'fdRateMode', 'blockLen', 'blockLenRatio', 'blockLenBaseRatio', ...
  'blockDurationUs', 'inBlockFdDriftHz', 'inBlockQuadPhaseRad', ...
  'truthGap', 'aliasGap1', 'aliasGap2', 'centerDeltaTruth', ...
  'centerDeltaFinal', 'minDeltaTruthHz', 'minDeltaFinalHz'});
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

function aggregateTable = localBuildAggregateTable(summaryKnown, summaryUnknown)
%LOCALBUILDAGGREGATETABLE Merge available known/unknown summary tables.

tableCell = {};
if istable(summaryKnown) && height(summaryKnown) > 0
  tableCell{end + 1, 1} = summaryKnown;
end
if istable(summaryUnknown) && height(summaryUnknown) > 0
  tableCell{end + 1, 1} = summaryUnknown;
end
if isempty(tableCell)
  aggregateTable = table();
else
  aggregateTable = vertcat(tableCell{:});
end
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
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
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
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
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
  xline(0, ':', 'center', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
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

function halfWidthHz = localResolveHalfWidth(scanConfig, aliasStepHz)
%LOCALRESOLVEHALFWIDTH Resolve scan half-width around one center.

halfWidthHz = localGetFieldOrDefault(scanConfig, 'scanHalfWidthHz', []);
if isempty(halfWidthHz)
  numAliasSide = localGetFieldOrDefault(scanConfig, 'numAliasSide', 2);
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
  xline(k * aliasStepHz, ':', 'HandleVisibility', 'off');
end
end

function tracker = localStartProgress(totalCount, progressLabel)
%LOCALSTARTPROGRESS Start a client-side progressbar when available.

tracker = struct('isActive', false, 'label', string(progressLabel));
fprintf('  fdRef line scan                  : %s (%d points)\n', progressLabel, totalCount);
if totalCount <= 1 || exist('progressbar', 'file') ~= 2
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  tracker.isActive = true;
catch
  tracker.isActive = false;
end
end

function localAdvanceProgress(tracker)
%LOCALADVANCEPROGRESS Advance an active serial progressbar.

if isstruct(tracker) && localGetFieldOrDefault(tracker, 'isActive', false)
  progressbar('advance');
end
end

function localFinishProgress(tracker)
%LOCALFINISHPROGRESS Close an active progressbar without touching results.

if isstruct(tracker) && localGetFieldOrDefault(tracker, 'isActive', false)
  try
    progressbar('end');
  catch
  end
end
end

function textVal = localFormatIntegerRow(value)
%LOCALFORMATINTEGERROW Format an integer vector for compact logging.

value = reshape(value, 1, []);
if isempty(value)
  textVal = "(empty)";
else
  textVal = string(strtrim(sprintf('%d ', value)));
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

