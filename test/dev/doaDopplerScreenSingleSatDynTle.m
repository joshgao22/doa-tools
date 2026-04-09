% doaDopplerScreenSingleSatDynTle
% Single-satellite dynamic screening for reference-satellite selection.
%
% Purpose:
% 1) Screen all satellites visible at the reference epoch and keep candidates with
%    reliable single-satellite estimation behavior over the selected dynamic window;
% 2) Provide a reference-satellite shortlist before any pairwise multi-satellite
%    fusion study, instead of directly sweeping all satellite pairs.
%
% Basic idea:
% 1) Filter satellites by visibility at the reference epoch and by continuous access
%    over the selected multi-frame window;
% 2) For each candidate satellite, run single-satellite Monte Carlo evaluation on
%    the same dynamic window;
% 3) Evaluate both SS-SF-DoA and SS-SF-Static metrics and record the dynamic
%    strength of the link;
% 4) Export two summary tables:
%    - Best single-satellite accuracy;
%    - Recommended reference-satellite candidates.
%
% Main criteria:
% - elevRefDeg / elevMinDeg: higher and more stable elevation usually indicates a
%   cleaner line of sight and more favorable geometry;
% - doaAngleRmseDeg: stability of the single-frame DoA-only coarse estimate;
% - staticAngleRmseDeg: single-frame static DoA-Doppler joint estimation accuracy;
% - staticFdRmseHz: stability of the reference-link fdRef estimate;
% - fdRateFitHzPerSec and quadPhaseWinRad: dynamic strength indicators over the
%   selected observation window.
%
% Practical notes:
% - A good reference satellite is not determined by the smallest angle error alone.
%   High elevation, stable fdRef estimation, and moderate dynamic stress should be
%   considered together;
% - If a non-trivial dual-satellite pair is desired later, it is better to first
%   shortlist several strong reference candidates here and then rank the second
%   satellite with a separate fixed-reference script, rather than directly choosing
%   two satellites that are both close to zenith.
clear(); close all; clc;
%%
%[text] ## Parameters
numUsr = 1;
numSym = 32;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 18e9;
wavelen = 299792458 / carrierFreq;
rng(253);

elemSpace = wavelen / 2;
numElem = [4 4];
snrDb = 10;
numRepeat = 30;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);

utc0 = datetime([2026, 03, 18, 17, 08, 00], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(localResolveTlePath("statlink_20260318.tle"));
arrUpa = createUpa(numElem, elemSpace);

numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = round((numFrame + 1) / 2);
utcVec = utc0 + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

minUsrElevationDeg = 15;
maxSatOffAxisDeg = 55;
fdRange = [-4e5, 4e5];
numTopReport = 8;
minElevRefCandidateDeg = 45;
maxQuadPhaseRefCandidateRad = 0.9;
maxDoaRmseRefCandidateDeg = 0.05;
maxStaticFdRmseRefCandidateHz = 150;
optVerbose = false;

if numUsr ~= 1
  error('doaDopplerScreenSingleSatDynTle:OnlySingleUserSupported', ...
    'This script currently supports one user only.');
end
%%
%[text] ## Visible satellites at the reference epoch
[satIdxRef, satAccessRef] = findVisibleSatFromTle(utc0, tle, usrLla, ...
  minUsrElevationDeg, maxSatOffAxisDeg);
candidateSatIdxGlobal = reshape(satIdxRef.available, 1, []);
numCandidate = numel(candidateSatIdxGlobal);

if numCandidate == 0
  error('doaDopplerScreenSingleSatDynTle:NoVisibleSatellite', ...
    'No available satellite is found at the requested UTC epoch.');
end

fprintf('\n========== Reference-epoch available satellites ==========%s', newline);
disp(table(candidateSatIdxGlobal(:), satAccessRef.usrElevationDeg(candidateSatIdxGlobal, 1), ...
  'VariableNames', {'globalSatIdx', 'usrElevationDeg'}));
%%
%[text] ## Pilot waveform
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
%%
%[text] ## Shared snapshot / estimator options
snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticOptBase = struct();
staticOptBase.useLogObjective = true;
%%
%[text] ## Screen satellites one by one under a dynamic window
screenResult = repmat(localBuildScreenResultTemplate(), 1, numCandidate);
progressbar('reset', numCandidate);

parfor iCand = 1:numCandidate
  satIdxGlobal = candidateSatIdxGlobal(iCand);
  % progressbar('advance');
  fprintf('Evaluating candidate %d/%d (globalSatIdx = %d) ...\n', ...
    iCand, numCandidate, satIdxGlobal);

  screenResult(iCand) = localEvaluateSingleSatCandidate( ...
    satIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, ...
    pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
    gridSize, searchRange, E, fdRange, ...
    doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef);
end
%%
%[text] ## Build report tables
statusTable = struct2table(screenResult);
statusTable = movevars(statusTable, {'globalSatIdx', 'satName', 'isAllAccess', 'skipReason'}, ...
  'Before', 1);

rankMask = statusTable.isAllAccess & isfinite(statusTable.staticAngleRmseDeg) & ...
  statusTable.staticResolveRate > 0;
staticRankTable = statusTable(rankMask, :);

if ~isempty(staticRankTable)
  staticRankTable = sortrows(staticRankTable, ...
    {'staticAngleRmseDeg', 'doaAngleRmseDeg', 'staticFdRmseHz', 'elevRefDeg'}, ...
    {'ascend', 'ascend', 'ascend', 'descend'});
  staticRankTable.rank = (1:height(staticRankTable)).';
  staticRankTable = movevars(staticRankTable, 'rank', 'Before', 'globalSatIdx');
else
  staticRankTable = table();
end

refCandidateMask = rankMask & ...
  statusTable.elevRefDeg >= minElevRefCandidateDeg & ...
  statusTable.quadPhaseWinRad <= maxQuadPhaseRefCandidateRad & ...
  statusTable.doaAngleRmseDeg <= maxDoaRmseRefCandidateDeg & ...
  statusTable.staticFdRmseHz <= maxStaticFdRmseRefCandidateHz;
statusTable.isRefCandidate = refCandidateMask;
refRankTable = statusTable(refCandidateMask, :);

if ~isempty(refRankTable)
  refRankTable = sortrows(refRankTable, ...
    {'staticAngleRmseDeg', 'staticFdRmseHz', 'quadPhaseWinRad', 'elevRefDeg'}, ...
    {'ascend', 'ascend', 'ascend', 'descend'});
  refRankTable.refRank = (1:height(refRankTable)).';
  refRankTable = movevars(refRankTable, 'refRank', 'Before', 'globalSatIdx');
else
  refRankTable = table();
end

staticTopCount = min(numTopReport, height(staticRankTable));
if staticTopCount > 0
  topStaticTable = staticRankTable(1:staticTopCount, :);
else
  topStaticTable = table();
end

refTopCount = min(numTopReport, height(refRankTable));
if refTopCount > 0
  topRefTable = refRankTable(1:refTopCount, :);
else
  topRefTable = table();
end

fprintf('\n========== Dynamic single-satellite screening summary ==========%s', newline);
disp(statusTable(:, {'globalSatIdx', 'satName', 'isAllAccess', 'isRefCandidate', 'elevRefDeg', ...
  'elevMinDeg', 'fdRefTrueHz', 'fdRateFitHzPerSec', 'deltaFdWinHz', ...
  'quadPhaseWinRad', 'localDirDriftMaxDeg', 'doaAngleRmseDeg', ...
  'staticAngleRmseDeg', 'staticFdRmseHz', 'staticAngleCrbStdDeg', ...
  'staticResolveRate', 'skipReason'}));

fprintf('\n========== Best single-satellite accuracy (sorted) ==========%s', newline);
if isempty(topStaticTable)
  disp('No continuously-accessible satellite passes the screening conditions.');
else
  disp(topStaticTable(:, {'rank', 'globalSatIdx', 'satName', 'elevRefDeg', ...
    'elevMinDeg', 'fdRateFitHzPerSec', 'quadPhaseWinRad', ...
    'doaAngleRmseDeg', 'staticAngleRmseDeg', 'staticFdRmseHz', ...
    'staticAngleCrbStdDeg', 'staticResolveRate'}));
end

fprintf('\n========== Recommended reference-satellite candidates ==========%s', newline);
if isempty(topRefTable)
  disp('No satellite satisfies the stricter reference-candidate conditions.');
else
  disp(topRefTable(:, {'refRank', 'globalSatIdx', 'satName', 'elevRefDeg', ...
    'elevMinDeg', 'fdRateFitHzPerSec', 'quadPhaseWinRad', ...
    'doaAngleRmseDeg', 'staticAngleRmseDeg', 'staticFdRmseHz', ...
    'staticAngleCrbStdDeg', 'staticResolveRate'}));
end
%%
%[text] ## Quick plots
if ~isempty(staticRankTable)
  localPlotRankSummary(staticRankTable, refRankTable, numTopReport);
end

function satResult = localEvaluateSingleSatCandidate( ...
  satIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
  minUsrElevationDeg, maxSatOffAxisDeg, ...
  pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
  gridSize, searchRange, E, fdRange, ...
  doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef)
%LOCALEVALUATESINGLESATCANDIDATE Evaluate one satellite under a dynamic window.

satResult = localBuildScreenResultTemplate();
satResult.globalSatIdx = satIdxGlobal;
satResult.satName = localExtractSatName(satAccessRef, satIdxGlobal);
satResult.numRepeat = numRepeat;

try
  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdxGlobal, [], arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, "satellite", satIdxGlobal);
  sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

  elevSeqDeg = localExtractElevationSeq(sceneSeq);
  satResult.elevRefDeg = elevSeqDeg(sceneSeq.refFrameIdx);
  satResult.elevMeanDeg = mean(elevSeqDeg);
  satResult.elevMinDeg = min(elevSeqDeg);
  satResult.elevMaxDeg = max(elevSeqDeg);
  satResult.isAllAccess = all(sceneSeq.access(:));

  [satResult.localDirDriftMaxDeg, satResult.localDirDriftRmsDeg] = ...
    localCalcLocalDirDrift(sceneSeq);

  if ~satResult.isAllAccess
    satResult.skipReason = "notAllAccessOverWindow";
    return;
  end

  linkParamCell = cell(1, sceneSeq.numFrame);
  for iFrame = 1:sceneSeq.numFrame
    linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
  end

  truthDyn = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, waveInfo.sampleRate, 1);
  linkParamRef = linkParamCell{sceneSeq.refFrameIdx};

  satResult.fdRefTrueHz = truthDyn.fdRefSeries(sceneSeq.refFrameIdx);
  satResult.fdRateFitHzPerSec = truthDyn.fdRateFit;
  satResult.deltaFdWinHz = truthDyn.deltaFdWin;
  satResult.quadPhaseWinRad = truthDyn.quadPhaseWin;

  steeringRef = getSceneSteering(sceneRef, wavelen);
  pathGain = ones(sceneRef.numSat, sceneRef.numUser);
  view = buildDoaDopplerEstView(sceneRef, [], gridSize, searchRange, E);
  doaGrid = view.doaGrid;
  latlonTrueDeg = usrLla(1:2, 1);

  crbOpt = struct();
  crbOpt.doaType = 'latlon';
  [crbStatic, ~] = crbPilotSfDoaDoppler( ...
    sceneRef, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    latlonTrueDeg, satResult.fdRefTrueHz, 1, pwrNoise, crbOpt);
  satResult.staticAngleCrbStdDeg = projectCrbToAngleMetric(crbStatic(1:2, 1:2), ...
    latlonTrueDeg, 'latlon');
  satResult.fdRefCrbStdHz = sqrt(max(real(crbStatic(3, 3)), 0));

  doaErrSq = nan(numRepeat, 1);
  staticErrSq = nan(numRepeat, 1);
  fdErrSq = nan(numRepeat, 1);
  doaResolved = false(numRepeat, 1);
  staticResolved = false(numRepeat, 1);

  for iRepeat = 1:numRepeat
    [rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
      steeringRef, pilotWave, linkParamRef, carrierFreq, waveInfo.sampleRate, ...
      pwrNoise, pathGain, snapOpt);

    estDoa = estimatorDoaMlePilotOpt( ...
      sceneRef.array, wavelen, rxSig, pilotWave, doaGrid, [], optVerbose, doaOnlyOpt);
    latlonDoaDeg = getDoaDopplerLatlonEst(estDoa);
    doaResolved(iRepeat) = localReadResolvedFlag(estDoa);

    if all(isfinite(latlonDoaDeg))
      doaAngleErrDeg = calcLatlonAngleError(latlonDoaDeg, latlonTrueDeg);
      doaErrSq(iRepeat) = doaAngleErrDeg .^ 2;
    end

    staticOpt = staticOptBase;
    if all(isfinite(latlonDoaDeg))
      staticOpt.initDoaParam = latlonDoaDeg(:);
    end

    estStatic = estimatorDoaDopplerMlePilotSfOpt( ...
      sceneRef, rxSig, pilotWave, carrierFreq, waveInfo.sampleRate, ...
      doaGrid, fdRange, 1, [], optVerbose, staticOpt);
    latlonStaticDeg = getDoaDopplerLatlonEst(estStatic);
    staticResolved(iRepeat) = localReadResolvedFlag(estStatic);

    if all(isfinite(latlonStaticDeg))
      staticAngleErrDeg = calcLatlonAngleError(latlonStaticDeg, latlonTrueDeg);
      staticErrSq(iRepeat) = staticAngleErrDeg .^ 2;
    end

    fdRefEstHz = getDoaDopplerFieldOrDefault(estStatic, 'fdRefEst', NaN);
    if isfinite(fdRefEstHz)
      fdErrSq(iRepeat) = (fdRefEstHz - satResult.fdRefTrueHz) .^ 2;
    end
  end

  satResult.doaAngleRmseDeg = localCalcRmse(doaErrSq);
  satResult.staticAngleRmseDeg = localCalcRmse(staticErrSq);
  satResult.staticFdRmseHz = localCalcRmse(fdErrSq);
  satResult.doaResolveRate = mean(doaResolved);
  satResult.staticResolveRate = mean(staticResolved);

catch ME
  satResult.skipReason = "runtimeError";
  satResult.errorMessage = string(ME.message);
end
end


function result = localBuildScreenResultTemplate()
%LOCALBUILDSCREENRESULTTEMPLATE Default result container for one satellite.

result = struct();
result.globalSatIdx = NaN;
result.satName = "";
result.isAllAccess = false;
result.skipReason = "";
result.errorMessage = "";

result.elevRefDeg = NaN;
result.elevMeanDeg = NaN;
result.elevMinDeg = NaN;
result.elevMaxDeg = NaN;

result.fdRefTrueHz = NaN;
result.fdRateFitHzPerSec = NaN;
result.deltaFdWinHz = NaN;
result.quadPhaseWinRad = NaN;
result.localDirDriftMaxDeg = NaN;
result.localDirDriftRmsDeg = NaN;

result.doaAngleRmseDeg = NaN;
result.staticAngleRmseDeg = NaN;
result.staticFdRmseHz = NaN;
result.doaResolveRate = NaN;
result.staticResolveRate = NaN;

result.staticAngleCrbStdDeg = NaN;
result.fdRefCrbStdHz = NaN;
result.numRepeat = NaN;
end


function value = localReadResolvedFlag(estResult)
%LOCALREADRESOLVEDFLAG Read one estimator resolved flag with fallback.

value = false;
if isempty(estResult)
  return;
end

value = logical(getDoaDopplerFieldOrDefault(estResult, 'isResolved', false));
end


function rmseValue = localCalcRmse(errSq)
%LOCALCALCRMSE RMSE with NaN protection.

validMask = isfinite(errSq);
if ~any(validMask)
  rmseValue = NaN;
  return;
end

rmseValue = sqrt(mean(errSq(validMask)));
end


function elevSeqDeg = localExtractElevationSeq(sceneSeq)
%LOCALEXTRACTELEVATIONSEQ Extract one satellite elevation sequence in degrees.

numFrame = sceneSeq.numFrame;
elevSeqDeg = zeros(1, numFrame);
for iFrame = 1:numFrame
  accessInfo = sceneSeq.sceneCell{iFrame}.accessInfo;
  elevSeqDeg(iFrame) = accessInfo.usrElevationDeg(1, 1);
end
end


function [maxDriftDeg, rmsDriftDeg] = localCalcLocalDirDrift(sceneSeq)
%LOCALCALCLOCALDIRDRIFT Build local-direction drift metrics for one satellite.

localDoaSeries = sceneSeq.localDoa;
if ndims(localDoaSeries) == 2
  localDoaSeries = reshape(localDoaSeries, 2, 1, 1);
end

numFrame = sceneSeq.numFrame;
refLocalDoa = reshape(localDoaSeries(:, 1, sceneSeq.refFrameIdx), 2, 1);
refDir = localAngleToUnitVec(refLocalDoa);

driftDeg = zeros(1, numFrame);
for iFrame = 1:numFrame
  currentLocalDoa = reshape(localDoaSeries(:, 1, iFrame), 2, 1);
  currentDir = localAngleToUnitVec(currentLocalDoa);
  cosineVal = max(-1, min(1, real(refDir' * currentDir)));
  driftDeg(iFrame) = acosd(cosineVal);
end

maxDriftDeg = max(driftDeg);
rmsDriftDeg = sqrt(mean(driftDeg .^ 2));
end


function dirVec = localAngleToUnitVec(localDoa)
%LOCALANGLETOTUNITVEC Convert [az; el] in radians to a 3x1 unit vector.

azRad = localDoa(1);
elRad = localDoa(2);

dirVec = [cos(elRad) * cos(azRad); ...
  cos(elRad) * sin(azRad); ...
  sin(elRad)];
end


function satName = localExtractSatName(satAccessRef, satIdxGlobal)
%LOCALEXTRACTSATNAME Read one TLE satellite name when available.

satName = "";
if ~isfield(satAccessRef, 'satName') || isempty(satAccessRef.satName)
  return;
end

satNameList = satAccessRef.satName;
if iscell(satNameList)
  if satIdxGlobal >= 1 && satIdxGlobal <= numel(satNameList)
    satName = string(satNameList{satIdxGlobal});
  end
  return;
end

if isstring(satNameList) || ischar(satNameList)
  satNameList = string(satNameList);
  if satIdxGlobal >= 1 && satIdxGlobal <= numel(satNameList)
    satName = satNameList(satIdxGlobal);
  end
end
end


function localPlotRankSummary(staticRankTable, refRankTable, numTopReport)
%LOCALPLOTRANKSUMMARY Quick-look plots for satellite screening results.

staticTopCount = min(numTopReport, height(staticRankTable));
if staticTopCount <= 0
  return;
end

staticTopTable = staticRankTable(1:staticTopCount, :);
labelText = compose('%d', staticTopTable.globalSatIdx);

figure('Name', 'Single-satellite screening: angle RMSE');
bar(staticTopTable.staticAngleRmseDeg);
grid on;
xlabel('Ranked satellite');
ylabel('SS-SF-Static angle RMSE (deg)');
title('Top satellites by single-satellite single-frame static accuracy');
set(gca, 'XTick', 1:staticTopCount, 'XTickLabel', labelText);

figure('Name', 'Single-satellite screening: RMSE vs elevation');
scatter(staticRankTable.elevRefDeg, staticRankTable.staticAngleRmseDeg, 48, ...
  staticRankTable.quadPhaseWinRad, 'filled');
grid on;
xlabel('Reference-frame elevation (deg)');
ylabel('SS-SF-Static angle RMSE (deg)');
title('Screening summary under the dynamic window');
colorbar;

if isempty(refRankTable)
  return;
end

refTopCount = min(numTopReport, height(refRankTable));
refTopTable = refRankTable(1:refTopCount, :);
refLabelText = compose('%d', refTopTable.globalSatIdx);

figure('Name', 'Reference-satellite candidate screening');
bar(refTopTable.staticAngleRmseDeg);
grid on;
xlabel('Reference-candidate satellite');
ylabel('SS-SF-Static angle RMSE (deg)');
title('Top stable reference-satellite candidates');
set(gca, 'XTick', 1:refTopCount, 'XTickLabel', refLabelText);
end

function tlePath = localResolveTlePath(fileName)
%LOCALRESOLVETLEPATH Resolve one TLE file from common test locations.

candidateList = { ...
  fullfile('/tle', fileName), ...
  fullfile('tle', fileName), ...
  fullfile('test', 'data', 'tle', fileName), ...
  fullfile('..', 'data', 'tle', fileName)};

for iPath = 1:numel(candidateList)
  if exist(candidateList{iPath}, 'file')
    tlePath = candidateList{iPath};
    return;
  end
end

error('doaDopplerScreenSingleSatDynTle:TleFileNotFound', ...
  'Unable to locate the TLE file "%s".', fileName);
end