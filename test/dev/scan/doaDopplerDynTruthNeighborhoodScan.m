% Development scan for the MF dynamic objective around multiple centers.
% This script extends the truth-neighborhood scan so the same 2-D slices can
% be evaluated around the truth point, the best static seed, and the final
% MF estimate. The goal is to compare whether the optimizer remains inside
% the truth basin or gets trapped in a different local peak.
%
% The plotted metric uses a reciprocal transform of the relative objective,
% so local minima of the original cost appear as peaks:
%   peakMetric = 1 / (deltaObj + displayFloor)
% Optional normalization scales the slice maximum to one.
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
scanCfg.numGrid = 41;
scanCfg.doaSpanDeg = [0.05; 0.05];
scanCfg.fdRefSpanHz = 80;
scanCfg.fdRateSpanHzPerSec = 1200;
scanCfg.runKnown = true;
scanCfg.runUnknown = true;
scanCfg.useParfor = true;
scanCfg.autoStartParpool = true;
scanCfg.showContourFigure = true;
scanCfg.showMeshFigure = true;
scanCfg.displayMode = "reciprocalDelta";
scanCfg.displayFloor = [];
scanCfg.displayNormalize = true;
scanCfg.displayContourLevels = 24;
scanCfg.centerListKnown = ["truth"; "staticSeed"; "finalEstimate"];
scanCfg.centerListUnknown = ["truth"; "staticSeed"; "cpKnownSeed"; "finalEstimate"];

%% Select two visible satellites and the reference satellite
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
refFrameIdx = sceneSeq.refFrameIdx;
sceneRef = sceneSeq.sceneCell{refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);

if sceneSeq.numSat ~= 2
  error('doaDopplerDynTruthNeighborhoodScan:InvalidNumSat', ...
    'This script expects exactly two selected satellites.');
end

otherSatIdxLocal = setdiff(1:sceneSeq.numSat, refSatIdxLocal);
if ~isscalar(otherSatIdxLocal)
  error('doaDopplerDynTruthNeighborhoodScan:InvalidOtherSat', ...
    'This script expects exactly one non-reference satellite.');
end
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

%% Pilot waveform and MF snapshots
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

%% Build single-satellite and multi-satellite views
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

%% Static transition cases and dynamic final estimates
% Reuse the same static-to-dynamic pipeline as the main dev script so the
% scan centers correspond to the actual MF initialization path.
doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

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

%% Dynamic model options for objective probing
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
scanResult.sceneInfo = struct( ...
  'selectedSatIdxGlobal', satIdx(:).', ...
  'refSatIdxGlobal', refSatIdxGlobal, ...
  'refSatIdxLocal', refSatIdxLocal, ...
  'otherSatIdxGlobal', otherSatIdxGlobal, ...
  'otherSatIdxLocal', otherSatIdxLocal, ...
  'refStateSource', string(refState.source));
scanResult.centerRegistry = centerRegistry;
scanResult.caseSummary = struct( ...
  'bestStaticDisplayName', string(bestStaticMsCase.displayName), ...
  'cpKnownDisplayName', string(caseDynMsKnown.displayName), ...
  'cpUnknownDisplayName', string(caseDynMsUnknown.displayName));

useParfor = localResolveParforFlag(scanCfg);
scanResult.parallel = struct('useParfor', useParfor, ...
  'poolSize', localGetParpoolSize());

%% CP-K scans around truth / static seed / final estimate
if scanCfg.runKnown
  modelOptKnown = baseDynOpt;
  modelOptKnown.fdRateMode = 'known';
  modelOptKnown.fdRateKnown = truth.fdRateFit;
  [modelKnown, ~, ~, modelOptKnown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, [], modelOptKnown);

  scanResult.cpKnown = localRunCenterScanSet(modelKnown, ...
    scanCfg.centerListKnown, centerRegistry, useParfor, scanCfg, 'known');
  scanResult.cpKnown.modelOpt = modelOptKnown;
end

%% CP-U scans around truth / static seed / final estimate
if scanCfg.runUnknown
  modelOptUnknown = baseDynOpt;
  modelOptUnknown.fdRateMode = 'unknown';
  [modelUnknown, ~, ~, modelOptUnknown] = buildDoaDopplerMfModel( ...
    viewMs.sceneSeq, viewMs.rxSigMf, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    viewMs.doaGrid, fdRange, fdRateRange, modelOptUnknown);

  scanResult.cpUnknown = localRunCenterScanSet(modelUnknown, ...
    scanCfg.centerListUnknown, centerRegistry, useParfor, scanCfg, 'unknown');
  scanResult.cpUnknown.modelOpt = modelOptUnknown;
end

%% Report and plot
if isfield(scanResult, 'cpKnown')
  disp('========== CP-K center-compare minima ==========');
  localDisplayModeSummary('CP-K', scanResult.cpKnown);
end
if isfield(scanResult, 'cpUnknown')
  disp('========== CP-U center-compare minima ==========');
  localDisplayModeSummary('CP-U', scanResult.cpUnknown);
end

if localGetFieldOrDefault(scanCfg, 'showContourFigure', true)
  if isfield(scanResult, 'cpKnown')
    localPlotModeContourFigure(scanResult.cpKnown, 'CP-K', scanCfg);
  end
  if isfield(scanResult, 'cpUnknown')
    localPlotModeContourFigure(scanResult.cpUnknown, 'CP-U', scanCfg);
  end
end
if localGetFieldOrDefault(scanCfg, 'showMeshFigure', true)
  if isfield(scanResult, 'cpKnown')
    localPlotModeMeshFigure(scanResult.cpKnown, 'CP-K');
  end
  if isfield(scanResult, 'cpUnknown')
    localPlotModeMeshFigure(scanResult.cpUnknown, 'CP-U');
  end
end

scanResult.utcRun = datetime('now', 'TimeZone', 'local');
scanResult.scriptTag = "doaDopplerDynTruthNeighborhoodScan_centerCompare";
saveFile = fullfile(fileparts(mfilename('fullpath')), ...
  sprintf('doaDopplerDynTruthNeighborhoodScan_%s.mat', ...
  datestr(now, 'yyyymmdd-HHMMSS')));
save(saveFile, 'scanResult', '-v7.3');
fprintf('Saved scanResult to %s\n', saveFile);


function modeResult = localRunCenterScanSet(model, centerNameList, centerRegistry, useParfor, scanCfg, modeTag)
%LOCALRUNCENTERSCANSET Run all 2-D slices for one mode across multiple centers.

centerNameList = reshape(string(centerNameList), [], 1);
modeResult = struct();
modeResult.modeTag = string(modeTag);
modeResult.centerNameList = centerNameList;
modeResult.comparePointList = localBuildComparePointList(centerRegistry, modeTag);
if strcmpi(modeTag, 'known')
  centerTemplate = struct( ...
    'centerName', "", ...
    'centerPoint', struct(), ...
    'comparePointList', struct([]), ...
    'doaSlice', struct(), ...
    'fdSlice', struct());
else
  centerTemplate = struct( ...
    'centerName', "", ...
    'centerPoint', struct(), ...
    'comparePointList', struct([]), ...
    'doaSlice', struct(), ...
    'fdRateSlice', struct());
end
modeResult.centerScan = repmat(centerTemplate, numel(centerNameList), 1);

for iCenter = 1:numel(centerNameList)
  centerName = centerNameList(iCenter);
  centerPoint = localResolveCenterPoint(centerRegistry, modeTag, centerName);
  centerScan = struct();
  centerScan.centerName = centerName;
  centerScan.centerPoint = centerPoint;
  centerScan.comparePointList = modeResult.comparePointList;

  latGrid = linspace(centerPoint.lat - scanCfg.doaSpanDeg(1), ...
    centerPoint.lat + scanCfg.doaSpanDeg(1), scanCfg.numGrid);
  lonGrid = linspace(centerPoint.lon - scanCfg.doaSpanDeg(2), ...
    centerPoint.lon + scanCfg.doaSpanDeg(2), scanCfg.numGrid);
  fdRefGrid = linspace(centerPoint.fdRef - scanCfg.fdRefSpanHz, ...
    centerPoint.fdRef + scanCfg.fdRefSpanHz, scanCfg.numGrid);
  fdRateGrid = linspace(centerPoint.fdRate - scanCfg.fdRateSpanHzPerSec, ...
    centerPoint.fdRate + scanCfg.fdRateSpanHzPerSec, scanCfg.numGrid);

  if strcmpi(modeTag, 'known')
    centerScan.doaSlice = localScanNeighborhood2d( ...
      model, centerPoint, 'lat', latGrid, 'lon', lonGrid, ...
      useParfor, scanCfg, centerName, modeResult.comparePointList);
    centerScan.fdSlice = localScanNeighborhood2d( ...
      model, centerPoint, 'lat', latGrid, 'fdRef', fdRefGrid, ...
      useParfor, scanCfg, centerName, modeResult.comparePointList);
  else
    centerScan.doaSlice = localScanNeighborhood2d( ...
      model, centerPoint, 'lat', latGrid, 'lon', lonGrid, ...
      useParfor, scanCfg, centerName, modeResult.comparePointList);
    centerScan.fdRateSlice = localScanNeighborhood2d( ...
      model, centerPoint, 'fdRef', fdRefGrid, 'fdRate', fdRateGrid, ...
      useParfor, scanCfg, centerName, modeResult.comparePointList);
  end

  modeResult.centerScan(iCenter) = centerScan;
end
end

function centerRegistry = localBuildCenterRegistry(truth, bestStaticMsCase, caseDynMsKnown, caseDynMsUnknown)
%LOCALBUILDCENTERREGISTRY Build the points used as scan centers and markers.

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

function pointList = localBuildComparePointList(centerRegistry, modeTag)
%LOCALBUILDCOMPAREPOINTLIST Build the landmark set shown on each plot.

pointList = repmat(struct('tag', "", 'point', struct()), 0, 1);
pointList(end + 1, 1) = struct('tag', "truth", 'point', centerRegistry.truth); %#ok<AGROW>
pointList(end + 1, 1) = struct('tag', "staticSeed", 'point', centerRegistry.staticSeed); %#ok<AGROW>
if strcmpi(modeTag, 'known')
  pointList(end + 1, 1) = struct('tag', "finalEstimate", 'point', centerRegistry.finalEstimateKnown); %#ok<AGROW>
else
  pointList(end + 1, 1) = struct('tag', "cpKnownSeed", 'point', centerRegistry.cpKnownSeed); %#ok<AGROW>
  pointList(end + 1, 1) = struct('tag', "finalEstimate", 'point', centerRegistry.finalEstimateUnknown); %#ok<AGROW>
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
    error('doaDopplerDynTruthNeighborhoodScan:UnknownCenter', ...
      'Unsupported center name: %s', centerName);
end
end

function scanInfo = localScanNeighborhood2d(model, centerPoint, xName, xGrid, yName, yGrid, ...
  useParfor, scanCfg, centerName, comparePointList)
%LOCALSCANNEIGHBORHOOD2D Evaluate one 2-D objective slice around one center.

nx = numel(xGrid);
ny = numel(yGrid);
objMat = nan(ny, nx);
residualMat = nan(ny, nx);
fdRefEvalMat = nan(ny, nx);
fdRateEvalMat = nan(ny, nx);
cohSatMat = nan(ny, nx, model.numSat);
objSatMat = nan(ny, nx, model.numSat);

probeTemplate = struct('doaParam', centerPoint.latlon(:), ...
  'fdRef', centerPoint.fdRef, 'fdRate', centerPoint.fdRate, ...
  'obj', NaN, 'residualNorm', NaN);

if useParfor
  parfor iy = 1:ny
    [objRow, residualRow, fdRefRow, fdRateRow, cohRow, objSatRow] = ...
      localScanNeighborhoodRow(model, centerPoint, probeTemplate, xName, xGrid, yName, yGrid(iy));
    objMat(iy, :) = objRow;
    residualMat(iy, :) = residualRow;
    fdRefEvalMat(iy, :) = fdRefRow;
    fdRateEvalMat(iy, :) = fdRateRow;
    cohSatMat(iy, :, :) = reshape(cohRow, [1, nx, model.numSat]);
    objSatMat(iy, :, :) = reshape(objSatRow, [1, nx, model.numSat]);
  end
else
  for iy = 1:ny
    [objRow, residualRow, fdRefRow, fdRateRow, cohRow, objSatRow] = ...
      localScanNeighborhoodRow(model, centerPoint, probeTemplate, xName, xGrid, yName, yGrid(iy));
    objMat(iy, :) = objRow;
    residualMat(iy, :) = residualRow;
    fdRefEvalMat(iy, :) = fdRefRow;
    fdRateEvalMat(iy, :) = fdRateRow;
    cohSatMat(iy, :, :) = reshape(cohRow, [1, nx, model.numSat]);
    objSatMat(iy, :, :) = reshape(objSatRow, [1, nx, model.numSat]);
  end
end

[minObj, minIdx] = min(objMat(:));
[minRow, minCol] = ind2sub(size(objMat), minIdx);
centerEval = localEvalProbeAtPoint(model, centerPoint, probeTemplate);

scanInfo = struct();
scanInfo.centerName = string(centerName);
scanInfo.centerPoint = centerPoint;
scanInfo.comparePointList = comparePointList;
scanInfo.xName = string(xName);
scanInfo.yName = string(yName);
scanInfo.xGrid = xGrid(:).';
scanInfo.yGrid = yGrid(:).';
scanInfo.objMat = objMat;
scanInfo.deltaObjMat = objMat - minObj;
scanInfo.residualMat = residualMat;
scanInfo.fdRefEvalMat = fdRefEvalMat;
scanInfo.fdRateEvalMat = fdRateEvalMat;
scanInfo.cohSatMat = cohSatMat;
scanInfo.objectiveSatMat = objSatMat;
scanInfo.minObj = minObj;
scanInfo.minPoint = centerPoint;
scanInfo.minPoint.(xName) = xGrid(minCol);
scanInfo.minPoint.(yName) = yGrid(minRow);
scanInfo.centerEval = centerEval;
scanInfo.centerDeltaObj = centerEval.obj - minObj;
scanInfo.centerGridIdx = localFindNearestGridIndex(xGrid, centerPoint.(xName));
scanInfo.centerGridIdy = localFindNearestGridIndex(yGrid, centerPoint.(yName));
scanInfo.display = localBuildDisplayMetric(scanInfo, scanCfg);
end

function [objRow, residualRow, fdRefRow, fdRateRow, cohRow, objSatRow] = ...
  localScanNeighborhoodRow(model, centerPoint, probeTemplate, xName, xGrid, yName, yVal)
%LOCALSCANNEIGHBORHOODROW Evaluate one y-slice row for one 2-D scan.

nx = numel(xGrid);
objRow = nan(1, nx);
residualRow = nan(1, nx);
fdRefRow = nan(1, nx);
fdRateRow = nan(1, nx);
cohRow = nan(nx, model.numSat);
objSatRow = nan(nx, model.numSat);

for ix = 1:nx
  pointCur = centerPoint;
  pointCur.(xName) = xGrid(ix);
  pointCur.(yName) = yVal;
  pointCur.latlon = [pointCur.lat; pointCur.lon];
  probeEval = localEvalProbeAtPoint(model, pointCur, probeTemplate);

  objRow(ix) = probeEval.obj;
  residualRow(ix) = probeEval.residualNorm;
  fdRefRow(ix) = localGetFieldOrDefault(probeEval, 'fdRef', NaN);
  fdRateRow(ix) = localGetFieldOrDefault(probeEval, 'fdRate', NaN);

  cohSat = reshape(localGetFieldOrDefault(probeEval, 'coherenceSat', []), 1, []);
  objSat = reshape(localGetFieldOrDefault(probeEval, 'objectiveSat', []), 1, []);
  numSatCopy = min(model.numSat, numel(cohSat));
  cohRow(ix, 1:numSatCopy) = cohSat(1:numSatCopy);
  numSatCopy = min(model.numSat, numel(objSat));
  objSatRow(ix, 1:numSatCopy) = objSat(1:numSatCopy);
end
end

function useParfor = localResolveParforFlag(scanCfg)
%LOCALRESOLVEPARFORFLAG Resolve the parfor execution mode for this scan.

useParfor = isfield(scanCfg, 'useParfor') && logical(scanCfg.useParfor);
if ~useParfor
  return;
end

useParfor = license('test', 'Distrib_Computing_Toolbox');
if ~useParfor
  warning('doaDopplerDynTruthNeighborhoodScan:NoPctLicense', ...
    'Parallel Computing Toolbox is unavailable. Fall back to serial scan.');
  return;
end

if isfield(scanCfg, 'autoStartParpool') && logical(scanCfg.autoStartParpool)
  poolObj = gcp('nocreate');
  if isempty(poolObj)
    parpool('Processes');
  end
end
end

function poolSize = localGetParpoolSize()
%LOCALGETPARPOOLSIZE Return the current parallel pool size.

poolSize = 0;
poolObj = gcp('nocreate');
if ~isempty(poolObj)
  poolSize = poolObj.NumWorkers;
end
end

function probeEval = localEvalProbeAtPoint(model, pointCur, probeTemplate)
%LOCALEVALPROBEATPOINT Evaluate one probe point through the shared helper.

probeEval = evalDoaDopplerMfProbePoint(model, pointCur, struct('probeTemplate', probeTemplate, 'debugEnable', false));
end

function idx = localFindNearestGridIndex(valGrid, valTarget)
%LOCALFINDNEARESTGRIDINDEX Return the nearest grid index for one scalar.

[~, idx] = min(abs(valGrid - valTarget));
end

function localDisplayModeSummary(modeLabel, modeResult)
%LOCALDISPLAYMODESUMMARY Print compact minima summaries for one mode.

for iCenter = 1:numel(modeResult.centerScan)
  centerScan = modeResult.centerScan(iCenter);
  fprintf('%s | %s | doaSlice: minObj = %.6g, centerDeltaObj = %.6g, min(%s,%s) = (%.6f, %.6f)\n', ...
    modeLabel, centerScan.centerName, centerScan.doaSlice.minObj, centerScan.doaSlice.centerDeltaObj, ...
    centerScan.doaSlice.xName, centerScan.doaSlice.yName, ...
    centerScan.doaSlice.minPoint.(char(centerScan.doaSlice.xName)), ...
    centerScan.doaSlice.minPoint.(char(centerScan.doaSlice.yName)));

  if isfield(centerScan, 'fdSlice')
    scanInfo = centerScan.fdSlice;
    sliceTag = 'fdSlice';
  else
    scanInfo = centerScan.fdRateSlice;
    sliceTag = 'fdRateSlice';
  end
  fprintf('%s | %s | %s: minObj = %.6g, centerDeltaObj = %.6g, min(%s,%s) = (%.6f, %.6f)\n', ...
    modeLabel, centerScan.centerName, sliceTag, scanInfo.minObj, scanInfo.centerDeltaObj, ...
    scanInfo.xName, scanInfo.yName, ...
    scanInfo.minPoint.(char(scanInfo.xName)), ...
    scanInfo.minPoint.(char(scanInfo.yName)));
end
end

function localPlotModeContourFigure(modeResult, modeLabel, scanCfg)
%LOCALPLOTMODECONTOURFIGURE Plot reciprocal-peak contour maps for one mode.

numCenter = numel(modeResult.centerScan);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s center-compare contour', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  centerScan = modeResult.centerScan(iCenter);
  sliceCell = localBuildModeSliceCell(centerScan);
  for iSlice = 1:numel(sliceCell)
    scanInfo = sliceCell{iSlice};
    subplot(numCenter, numel(sliceCell), (iCenter - 1) * numel(sliceCell) + iSlice);
    contourf(scanInfo.xGrid, scanInfo.yGrid, scanInfo.display.metricMat, ...
      scanCfg.displayContourLevels, 'LineColor', 'none');
    colorbar;
    hold on;
    localOverlayComparePoints2d(scanInfo);
    xlabel(char(scanInfo.xName), 'Interpreter', 'none');
    ylabel(char(scanInfo.yName), 'Interpreter', 'none');
    title(sprintf('%s | %s | %s-%s', ...
      modeLabel, scanInfo.centerName, scanInfo.xName, scanInfo.yName), ...
      'Interpreter', 'none');
    grid on;
    axis tight;
  end
end
end

function localPlotModeMeshFigure(modeResult, modeLabel)
%LOCALPLOTMODEMESHFIGURE Plot 3-D reciprocal peak meshes for one mode.

numCenter = numel(modeResult.centerScan);
if numCenter == 0
  return;
end

figure('Name', sprintf('%s center-compare mesh', modeLabel), 'Color', 'w');
for iCenter = 1:numCenter
  centerScan = modeResult.centerScan(iCenter);
  sliceCell = localBuildModeSliceCell(centerScan);
  for iSlice = 1:numel(sliceCell)
    scanInfo = sliceCell{iSlice};
    subplot(numCenter, numel(sliceCell), (iCenter - 1) * numel(sliceCell) + iSlice);
    [xMesh, yMesh] = meshgrid(scanInfo.xGrid, scanInfo.yGrid);
    zMesh = scanInfo.display.metricMat;
    mesh(xMesh, yMesh, zMesh);
    hold on;
    localOverlayComparePoints3d(scanInfo);
    xlabel(char(scanInfo.xName), 'Interpreter', 'none');
    ylabel(char(scanInfo.yName), 'Interpreter', 'none');
    zlabel(scanInfo.display.zLabel, 'Interpreter', 'none');
    title(sprintf('%s | %s | mesh %s-%s', ...
      modeLabel, scanInfo.centerName, scanInfo.xName, scanInfo.yName), ...
      'Interpreter', 'none');
    grid on;
    axis tight;
    view(45, 28);
  end
end
end

function sliceCell = localBuildModeSliceCell(centerScan)
%LOCALBUILDMODESLICECELL Collect the two slices for one center.

sliceCell = {centerScan.doaSlice};
if isfield(centerScan, 'fdSlice')
  sliceCell{end + 1} = centerScan.fdSlice;
else
  sliceCell{end + 1} = centerScan.fdRateSlice;
end
end

function localOverlayComparePoints2d(scanInfo)
%LOCALOVERLAYCOMPAREPOINTS2D Overlay truth/seed/final markers on one 2-D plot.

for iPoint = 1:numel(scanInfo.comparePointList)
  pointTag = scanInfo.comparePointList(iPoint).tag;
  pointVal = scanInfo.comparePointList(iPoint).point;
  [markerSpec, markerLabel] = localResolveMarkerStyle(pointTag);
  if ~localPointInsideGrid(scanInfo, pointVal)
    continue;
  end
  xVal = pointVal.(char(scanInfo.xName));
  yVal = pointVal.(char(scanInfo.yName));
  plot(xVal, yVal, markerSpec, 'LineWidth', 1.5, 'MarkerSize', 9);
  text(xVal, yVal, ['  ' char(markerLabel)], ...
    'Color', 'w', 'FontWeight', 'bold', 'Interpreter', 'none');
end

zMin = localGetGridZVal(scanInfo, scanInfo.minPoint.(char(scanInfo.xName)), ...
  scanInfo.minPoint.(char(scanInfo.yName)));
plot(scanInfo.minPoint.(char(scanInfo.xName)), ...
  scanInfo.minPoint.(char(scanInfo.yName)), 'wo', 'LineWidth', 1.5, 'MarkerSize', 8);
text(scanInfo.minPoint.(char(scanInfo.xName)), ...
  scanInfo.minPoint.(char(scanInfo.yName)), ...
  sprintf('  min %.3g', zMin), ...
  'Color', 'w', 'FontWeight', 'bold', 'Interpreter', 'none');
end

function localOverlayComparePoints3d(scanInfo)
%LOCALOVERLAYCOMPAREPOINTS3D Overlay truth/seed/final markers on one 3-D plot.

for iPoint = 1:numel(scanInfo.comparePointList)
  pointTag = scanInfo.comparePointList(iPoint).tag;
  pointVal = scanInfo.comparePointList(iPoint).point;
  [markerSpec, markerLabel] = localResolveMarkerStyle(pointTag);
  if ~localPointInsideGrid(scanInfo, pointVal)
    continue;
  end
  xVal = pointVal.(char(scanInfo.xName));
  yVal = pointVal.(char(scanInfo.yName));
  zVal = localGetGridZVal(scanInfo, xVal, yVal);
  plot3(xVal, yVal, zVal, markerSpec, 'LineWidth', 1.5, 'MarkerSize', 9);
  text(xVal, yVal, zVal, ['  ' char(markerLabel)], ...
    'Color', 'k', 'FontWeight', 'bold', 'Interpreter', 'none');
end

xMin = scanInfo.minPoint.(char(scanInfo.xName));
yMin = scanInfo.minPoint.(char(scanInfo.yName));
zMin = localGetGridZVal(scanInfo, xMin, yMin);
plot3(xMin, yMin, zMin, 'wo', 'LineWidth', 1.5, 'MarkerSize', 8);
text(xMin, yMin, zMin, sprintf('  min %.3g', zMin), ...
  'Color', 'k', 'FontWeight', 'bold', 'Interpreter', 'none');
end

function [markerSpec, markerLabel] = localResolveMarkerStyle(pointTag)
%LOCALRESOLVEMARKERSTYLE Resolve marker shape for one landmark.

switch lower(char(pointTag))
  case 'truth'
    markerSpec = 'wx';
    markerLabel = "truth";
  case 'staticseed'
    markerSpec = 'ws';
    markerLabel = "static";
  case 'finalestimate'
    markerSpec = 'wd';
    markerLabel = "final";
  otherwise
    markerSpec = 'w+';
    markerLabel = string(pointTag);
end
end

function isInside = localPointInsideGrid(scanInfo, pointVal)
%LOCALPOINTINSIDEGRID Return true when one point lies inside the current slice.

xVal = pointVal.(char(scanInfo.xName));
yVal = pointVal.(char(scanInfo.yName));
isInside = xVal >= min(scanInfo.xGrid) && xVal <= max(scanInfo.xGrid) && ...
  yVal >= min(scanInfo.yGrid) && yVal <= max(scanInfo.yGrid);
end

function displayInfo = localBuildDisplayMetric(scanInfo, scanCfg)
%LOCALBUILDDISPLAYMETRIC Build a reciprocal display surface from delta objective.

deltaObj = scanInfo.deltaObjMat;
modeName = string(localGetFieldOrDefault(scanCfg, 'displayMode', 'reciprocalDelta'));

displayFloor = localGetFieldOrDefault(scanCfg, 'displayFloor', []);
if isempty(displayFloor)
  positiveDelta = deltaObj(isfinite(deltaObj) & deltaObj > 0);
  if isempty(positiveDelta)
    displayFloor = 1;
  else
    displayFloor = min(positiveDelta);
    displayFloor = max(displayFloor, 1e-12);
  end
end

switch lower(char(modeName))
  case 'reciprocaldelta'
    metricRaw = 1 ./ max(deltaObj + displayFloor, displayFloor);
    label = sprintf('1/(\\DeltaJ + %.3g)', displayFloor);
    zLabel = 'reciprocal peak';
  otherwise
    error('doaDopplerDynTruthNeighborhoodScan:UnknownDisplayMode', ...
      'Unsupported displayMode: %s', modeName);
end

if localGetFieldOrDefault(scanCfg, 'displayNormalize', true)
  metricScale = max(metricRaw(:), [], 'omitnan');
  if isempty(metricScale) || ~isfinite(metricScale) || metricScale <= 0
    metricScale = 1;
  end
  metricMat = metricRaw / metricScale;
  label = sprintf('%s (normalized)', label);
else
  metricMat = metricRaw;
end

displayInfo = struct();
displayInfo.mode = modeName;
displayInfo.floor = displayFloor;
displayInfo.metricMat = metricMat;
displayInfo.metricRaw = metricRaw;
displayInfo.label = label;
displayInfo.zLabel = zLabel;
end

function zVal = localGetGridZVal(scanInfo, xVal, yVal)
%LOCALGETGRIDZVAL Read the plotted z value at the nearest grid point.

ix = localFindNearestGridIndex(scanInfo.xGrid, xVal);
iy = localFindNearestGridIndex(scanInfo.yGrid, yVal);
zVal = scanInfo.display.metricMat(iy, ix);
end

function doaParam = localResolveDoaParam(estResult, fallbackDoa)
%LOCALRESOLVEDOAPARAM Resolve one 2x1 DoA parameter vector with fallback.

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
