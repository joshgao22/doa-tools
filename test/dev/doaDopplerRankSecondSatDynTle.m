% doaDopplerRankSecondSatDynTle
% Fixed-reference screening for the second satellite in a dual-satellite pair.
%
% Purpose:
% 1) Assume the reference satellite has already been selected, and evaluate the gain
%    brought by each candidate second satellite in static two-satellite fusion;
% 2) Rank candidate second satellites while avoiding pairs that are too close to the
%    reference geometry and therefore too trivial.
%
% Basic idea:
% 1) Fix one reference satellite refSatIdxGlobal;
% 2) Enumerate all other satellites visible at the reference epoch as candidate
%    second satellites;
% 3) For each pair, sweep the sat2 weight alpha and compare the pair result against
%    the reference-only baseline;
% 4) Record the best sat2 weight, the best angle RMSE, the fd RMSE, and the
%    reference-frame differential Doppler between the two satellites.
%
% Main criteria:
% - bestAngleRmseDeg: best static angle error achieved by the pair after sat2 weight
%   selection;
% - angleGainVsRefDeg: angle-error improvement relative to the reference-only case;
% - bestFdRmseHz: stability of fdRef after fusion;
% - deltaFdRefHz: differential Doppler of the second satellite with respect to the
%   reference satellite at the reference frame.
%
% Pair-selection principles:
% - If only the most stable result is needed, bestAngleRmseDeg and angleGainVsRefDeg
%   can be prioritized;
% - If the pair should not be overly trivial, combinations with very small
%   deltaFdRefHz should not be favored blindly;
% - A more meaningful second satellite usually satisfies three conditions:
%   1) it provides positive angle gain;
%   2) it is not an obviously weak single-satellite link;
%   3) it has a reasonable separation from the reference satellite in elevation,
%      Doppler, or both.
%
% Practical notes:
% - This script is intended for fixed-reference ranking of the second satellite;
% - For a paper-quality pair selection, it is usually helpful to keep two outcomes:
%   one pair for the most stable overall performance, and another pair that is less
%   trivial from the geometry / Doppler-separation viewpoint.
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
numRepeat = 20;
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
weightSweepAlpha = [0.25, 0.5, 1];
numTopReport = 8;
optVerbose = false;

refSatIdxGlobal = 4154;

if numUsr ~= 1
  error('doaDopplerRankSecondSatDynTle:OnlySingleUserSupported', ...
    'This script currently supports one user only.');
end
%%
%[text] ## Visible satellites at the reference epoch
[satIdxRef, satAccessRef] = findVisibleSatFromTle(utc0, tle, usrLla, ...
  minUsrElevationDeg, maxSatOffAxisDeg);
availableSatIdxGlobal = reshape(satIdxRef.available, 1, []);

if ~ismember(refSatIdxGlobal, availableSatIdxGlobal)
  error('doaDopplerRankSecondSatDynTle:ReferenceNotVisible', ...
    'The requested reference satellite %d is not visible at the reference epoch.', refSatIdxGlobal);
end

candidateSatIdxGlobal = setdiff(availableSatIdxGlobal, refSatIdxGlobal, 'stable');
numCandidate = numel(candidateSatIdxGlobal);

if numCandidate == 0
  error('doaDopplerRankSecondSatDynTle:NoSecondSatellite', ...
    'No second-satellite candidate is available at the requested UTC epoch.');
end

fprintf('\n========== Fixed reference satellite ==========%s', newline);
disp(table(refSatIdxGlobal, satAccessRef.usrElevationDeg(refSatIdxGlobal, 1), ...
  'VariableNames', {'globalSatIdx', 'usrElevationDeg'}));

fprintf('\n========== Candidate second satellites ==========%s', newline);
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
staticOptBase.initDoaHalfWidth = [0.01; 0.01];
%%
%[text] ## Reference-only baseline
refOnlyResult = localEvaluateReferenceOnly( ...
  refSatIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
  minUsrElevationDeg, maxSatOffAxisDeg, ...
  pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
  gridSize, searchRange, E, fdRange, ...
  doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef);

if ~refOnlyResult.isAllAccess
  error('doaDopplerRankSecondSatDynTle:ReferenceAccessLost', ...
    'The reference satellite is not continuously accessible over the dynamic window.');
end
%%
%[text] ## Rank second satellites for the fixed reference
pairResult = repmat(localBuildPairResultTemplate(weightSweepAlpha), 1, numCandidate);

parfor iCand = 1:numCandidate
  satIdxGlobal = candidateSatIdxGlobal(iCand);
  fprintf('Evaluating pair %d/%d (ref = %d, pair = %d) ...\n', ...
    iCand, numCandidate, refSatIdxGlobal, satIdxGlobal);

  pairResult(iCand) = localEvaluatePairCandidate( ...
    refSatIdxGlobal, satIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, ...
    pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
    gridSize, searchRange, E, fdRange, weightSweepAlpha, ...
    doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef);
end
%%
%[text] ## Report tables
refTable = struct2table(refOnlyResult);
pairTable = struct2table(pairResult);
pairTable = movevars(pairTable, {'pairSatIdxGlobal', 'pairSatName', 'isAllAccess', 'skipReason'}, ...
  'Before', 1);

rankMask = pairTable.isAllAccess & isfinite(pairTable.bestAngleRmseDeg) & pairTable.bestResolveRate > 0;
rankTable = pairTable(rankMask, :);

if ~isempty(rankTable)
  rankTable.angleGainVsRefDeg = refOnlyResult.staticAngleRmseDeg - rankTable.bestAngleRmseDeg;
  rankTable = sortrows(rankTable, ...
    {'bestAngleRmseDeg', 'bestFdRmseHz', 'pairElevRefDeg'}, ...
    {'ascend', 'ascend', 'descend'});
  rankTable.rank = (1:height(rankTable)).';
  rankTable = movevars(rankTable, 'rank', 'Before', 'pairSatIdxGlobal');
else
  rankTable = table();
end

topCount = min(numTopReport, height(rankTable));
if topCount > 0
  topTable = rankTable(1:topCount, :);
else
  topTable = table();
end

fprintf('\n========== Reference-only baseline ==========%s', newline);
disp(refTable(:, {'refSatIdxGlobal', 'satName', 'elevRefDeg', 'elevMinDeg', ...
  'fdRateFitHzPerSec', 'quadPhaseWinRad', 'doaAngleRmseDeg', ...
  'staticAngleRmseDeg', 'staticFdRmseHz', 'staticResolveRate'}));

fprintf('\n========== Fixed-reference pair screening summary ==========%s', newline);
disp(pairTable(:, {'pairSatIdxGlobal', 'pairSatName', 'isAllAccess', ...
  'pairElevRefDeg', 'pairElevMinDeg', 'deltaFdRefHz', 'bestAlphaSat2', ...
  'bestAngleRmseDeg', 'bestFdRmseHz', 'bestResolveRate', 'skipReason'}));

fprintf('\n========== Recommended second satellites (sorted) ==========%s', newline);
if isempty(topTable)
  disp('No second-satellite candidate passes the fixed-reference screening conditions.');
else
  disp(topTable(:, {'rank', 'pairSatIdxGlobal', 'pairSatName', 'pairElevRefDeg', ...
    'pairElevMinDeg', 'deltaFdRefHz', 'bestAlphaSat2', 'angleGainVsRefDeg', ...
    'bestAngleRmseDeg', 'bestFdRmseHz', 'bestResolveRate'}));
end

if ~isempty(rankTable)
  localPlotPairRankSummary(rankTable, numTopReport);
end

function result = localEvaluateReferenceOnly( ...
  refSatIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
  minUsrElevationDeg, maxSatOffAxisDeg, ...
  pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
  gridSize, searchRange, E, fdRange, ...
  doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef)
%LOCALEVALUATEREFERENCEONLY Evaluate the fixed reference satellite alone.

result = localBuildRefOnlyResultTemplate();
result.refSatIdxGlobal = refSatIdxGlobal;
result.satName = localExtractSatName(satAccessRef, refSatIdxGlobal);
result.numRepeat = numRepeat;

try
  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, refSatIdxGlobal, [], arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, "satellite", refSatIdxGlobal);
  sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

  elevSeqDeg = localExtractElevationSeq(sceneSeq);
  result.elevRefDeg = elevSeqDeg(sceneSeq.refFrameIdx);
  result.elevMinDeg = min(elevSeqDeg);
  result.isAllAccess = all(sceneSeq.access(:));
  if ~result.isAllAccess
    result.skipReason = "notAllAccessOverWindow";
    return;
  end

  linkParamCell = cell(1, sceneSeq.numFrame);
  for iFrame = 1:sceneSeq.numFrame
    linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
  end
  truthDyn = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, waveInfo.sampleRate, 1);
  linkParamRef = linkParamCell{sceneSeq.refFrameIdx};
  steeringRef = getSceneSteering(sceneRef, wavelen);

  result.fdRefTrueHz = truthDyn.fdRefSeries(sceneSeq.refFrameIdx);
  result.fdRateFitHzPerSec = truthDyn.fdRateFit;
  result.quadPhaseWinRad = truthDyn.quadPhaseWin;

  pathGain = ones(sceneRef.numSat, sceneRef.numUser);
  view = buildDoaDopplerEstView(sceneRef, [], gridSize, searchRange, E);
  latlonTrueDeg = usrLla(1:2, 1);

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
      sceneRef.array, wavelen, rxSig, pilotWave, view.doaGrid, [], optVerbose, doaOnlyOpt);
    latlonDoaDeg = getDoaDopplerLatlonEst(estDoa);
    doaResolved(iRepeat) = localReadResolvedFlag(estDoa);
    if all(isfinite(latlonDoaDeg))
      doaErrSq(iRepeat) = calcLatlonAngleError(latlonDoaDeg, latlonTrueDeg) .^ 2;
    end

    staticOpt = staticOptBase;
    if all(isfinite(latlonDoaDeg))
      staticOpt.initDoaParam = latlonDoaDeg(:);
    end

    estStatic = estimatorDoaDopplerMlePilotSfOpt( ...
      sceneRef, rxSig, pilotWave, carrierFreq, waveInfo.sampleRate, ...
      view.doaGrid, fdRange, 1, [], optVerbose, staticOpt);
    latlonStaticDeg = getDoaDopplerLatlonEst(estStatic);
    staticResolved(iRepeat) = localReadResolvedFlag(estStatic);
    if all(isfinite(latlonStaticDeg))
      staticErrSq(iRepeat) = calcLatlonAngleError(latlonStaticDeg, latlonTrueDeg) .^ 2;
    end

    fdRefEstHz = getDoaDopplerFieldOrDefault(estStatic, 'fdRefEst', NaN);
    if isfinite(fdRefEstHz)
      fdErrSq(iRepeat) = (fdRefEstHz - result.fdRefTrueHz) .^ 2;
    end
  end

  result.doaAngleRmseDeg = localCalcRmse(doaErrSq);
  result.staticAngleRmseDeg = localCalcRmse(staticErrSq);
  result.staticFdRmseHz = localCalcRmse(fdErrSq);
  result.doaResolveRate = mean(doaResolved);
  result.staticResolveRate = mean(staticResolved);

catch ME
  result.skipReason = "runtimeError";
  result.errorMessage = string(ME.message);
end
end


function result = localEvaluatePairCandidate( ...
  refSatIdxGlobal, pairSatIdxGlobal, tle, utcVec, usrLla, arrUpa, ...
  minUsrElevationDeg, maxSatOffAxisDeg, ...
  pilotWave, waveInfo, wavelen, carrierFreq, pwrNoise, ...
  gridSize, searchRange, E, fdRange, weightSweepAlpha, ...
  doaOnlyOpt, staticOptBase, snapOpt, optVerbose, numRepeat, satAccessRef)
%LOCALEVALUATEPAIRCANDIDATE Evaluate one fixed-reference satellite pair.

result = localBuildPairResultTemplate(weightSweepAlpha);
result.refSatIdxGlobal = refSatIdxGlobal;
result.refSatName = localExtractSatName(satAccessRef, refSatIdxGlobal);
result.pairSatIdxGlobal = pairSatIdxGlobal;
result.pairSatName = localExtractSatName(satAccessRef, pairSatIdxGlobal);
result.numRepeat = numRepeat;

try
  satIdxGlobal = [refSatIdxGlobal, pairSatIdxGlobal];
  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdxGlobal, [], arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, "satellite", refSatIdxGlobal);
  sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

  if sceneRef.numSat ~= 2
    result.skipReason = "invalidNumSat";
    return;
  end

  result.isAllAccess = all(sceneSeq.access(:));
  if ~result.isAllAccess
    result.skipReason = "notAllAccessOverWindow";
    return;
  end

  [refState, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
  pairSatIdxLocal = setdiff(1:sceneRef.numSat, refSatIdxLocal, 'stable');

  accessInfoRef = sceneRef.accessInfo;
  result.refElevRefDeg = accessInfoRef.usrElevationDeg(refSatIdxLocal, 1);
  result.pairElevRefDeg = accessInfoRef.usrElevationDeg(pairSatIdxLocal, 1);
  result.refElevMinDeg = min(localExtractElevationSeq(sceneSeq, refSatIdxLocal));
  result.pairElevMinDeg = min(localExtractElevationSeq(sceneSeq, pairSatIdxLocal));

  linkParamCell = cell(1, sceneSeq.numFrame);
  for iFrame = 1:sceneSeq.numFrame
    linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
  end
  truthDyn = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, waveInfo.sampleRate, 1);
  result.fdRefTrueHz = truthDyn.fdRefSeries(sceneSeq.refFrameIdx);
  result.fdRateFitHzPerSec = truthDyn.fdRateFit;
  result.deltaFdRefHz = truthDyn.deltaFdSeries(pairSatIdxLocal, sceneSeq.refFrameIdx);

  steeringRef = getSceneSteering(sceneRef, wavelen);
  linkParamRef = linkParamCell{sceneSeq.refFrameIdx};
  pathGain = ones(sceneRef.numSat, sceneRef.numUser);

  viewPair = buildDoaDopplerEstView(sceneRef, [], gridSize, searchRange, E);
  latlonTrueDeg = usrLla(1:2, 1);

  numWeight = numel(weightSweepAlpha);
  angleErrSq = nan(numRepeat, numWeight);
  fdErrSq = nan(numRepeat, numWeight);
  resolveMat = false(numRepeat, numWeight);

  for iRepeat = 1:numRepeat
    [rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
      steeringRef, pilotWave, linkParamRef, carrierFreq, waveInfo.sampleRate, ...
      pwrNoise, pathGain, snapOpt);

    estDoa = estimatorDoaMlePilotOpt( ...
      sceneRef.array, wavelen, rxSig, pilotWave, viewPair.doaGrid, [], optVerbose, doaOnlyOpt);
    latlonDoaDeg = getDoaDopplerLatlonEst(estDoa);

    staticOpt = staticOptBase;
    if all(isfinite(latlonDoaDeg))
      staticOpt.initDoaParam = latlonDoaDeg(:);
    end

    for iWeight = 1:numWeight
      currentOpt = staticOpt;
      satWeight = ones(sceneRef.numSat, 1);
      satWeight(pairSatIdxLocal) = weightSweepAlpha(iWeight);
      currentOpt.satWeight = satWeight;

      estStatic = estimatorDoaDopplerMlePilotSfOpt( ...
        sceneRef, rxSig, pilotWave, carrierFreq, waveInfo.sampleRate, ...
        viewPair.doaGrid, fdRange, 1, [], optVerbose, currentOpt);

      latlonStaticDeg = getDoaDopplerLatlonEst(estStatic);
      resolveMat(iRepeat, iWeight) = localReadResolvedFlag(estStatic);
      if all(isfinite(latlonStaticDeg))
        angleErrSq(iRepeat, iWeight) = calcLatlonAngleError(latlonStaticDeg, latlonTrueDeg) .^ 2;
      end

      fdRefEstHz = getDoaDopplerFieldOrDefault(estStatic, 'fdRefEst', NaN);
      if isfinite(fdRefEstHz)
        fdErrSq(iRepeat, iWeight) = (fdRefEstHz - result.fdRefTrueHz) .^ 2;
      end
    end
  end

  result.weightSweepAlpha = reshape(weightSweepAlpha, 1, []);
  result.angleRmsePerWeightDeg = arrayfun(@(idx) localCalcRmse(angleErrSq(:, idx)), 1:numWeight);
  result.fdRmsePerWeightHz = arrayfun(@(idx) localCalcRmse(fdErrSq(:, idx)), 1:numWeight);
  result.resolveRatePerWeight = mean(resolveMat, 1);

  validMask = isfinite(result.angleRmsePerWeightDeg) & result.resolveRatePerWeight > 0;
  if ~any(validMask)
    result.skipReason = "noResolvedWeight";
    return;
  end

  validIdx = find(validMask);
  bestIdx = validIdx(1);
  for iValid = 2:numel(validIdx)
    idx = validIdx(iValid);
    if result.angleRmsePerWeightDeg(idx) < result.angleRmsePerWeightDeg(bestIdx) || ...
        (result.angleRmsePerWeightDeg(idx) == result.angleRmsePerWeightDeg(bestIdx) && ...
         result.fdRmsePerWeightHz(idx) < result.fdRmsePerWeightHz(bestIdx))
      bestIdx = idx;
    end
  end

  result.bestAlphaSat2 = weightSweepAlpha(bestIdx);
  result.bestAngleRmseDeg = result.angleRmsePerWeightDeg(bestIdx);
  result.bestFdRmseHz = result.fdRmsePerWeightHz(bestIdx);
  result.bestResolveRate = result.resolveRatePerWeight(bestIdx);
  result.refStateSource = string(refState.source);

catch ME
  result.skipReason = "runtimeError";
  result.errorMessage = string(ME.message);
end
end


function result = localBuildRefOnlyResultTemplate()
%LOCALBUILDREFONLYRESULTTEMPLATE Default result container for the reference baseline.

result = struct();
result.refSatIdxGlobal = NaN;
result.satName = "";
result.isAllAccess = false;
result.skipReason = "";
result.errorMessage = "";
result.elevRefDeg = NaN;
result.elevMinDeg = NaN;
result.fdRefTrueHz = NaN;
result.fdRateFitHzPerSec = NaN;
result.quadPhaseWinRad = NaN;
result.doaAngleRmseDeg = NaN;
result.staticAngleRmseDeg = NaN;
result.staticFdRmseHz = NaN;
result.doaResolveRate = NaN;
result.staticResolveRate = NaN;
result.numRepeat = NaN;
end


function result = localBuildPairResultTemplate(weightSweepAlpha)
%LOCALBUILDPAIRRESULTTEMPLATE Default result container for one second-satellite candidate.

numWeight = numel(weightSweepAlpha);
result = struct();
result.refSatIdxGlobal = NaN;
result.refSatName = "";
result.pairSatIdxGlobal = NaN;
result.pairSatName = "";
result.isAllAccess = false;
result.skipReason = "";
result.errorMessage = "";
result.refElevRefDeg = NaN;
result.refElevMinDeg = NaN;
result.pairElevRefDeg = NaN;
result.pairElevMinDeg = NaN;
result.fdRefTrueHz = NaN;
result.fdRateFitHzPerSec = NaN;
result.deltaFdRefHz = NaN;
result.weightSweepAlpha = reshape(weightSweepAlpha, 1, []);
result.angleRmsePerWeightDeg = nan(1, numWeight);
result.fdRmsePerWeightHz = nan(1, numWeight);
result.resolveRatePerWeight = nan(1, numWeight);
result.bestAlphaSat2 = NaN;
result.bestAngleRmseDeg = NaN;
result.bestFdRmseHz = NaN;
result.bestResolveRate = NaN;
result.refStateSource = "";
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


function elevSeqDeg = localExtractElevationSeq(sceneSeq, localSatIdx)
%LOCALEXTRACTELEVATIONSEQ Extract one satellite elevation sequence in degrees.

if nargin < 2 || isempty(localSatIdx)
  localSatIdx = 1;
end

numFrame = sceneSeq.numFrame;
elevSeqDeg = zeros(1, numFrame);
for iFrame = 1:numFrame
  accessInfo = sceneSeq.sceneCell{iFrame}.accessInfo;
  elevSeqDeg(iFrame) = accessInfo.usrElevationDeg(localSatIdx, 1);
end
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


function localPlotPairRankSummary(rankTable, numTopReport)
%LOCALPLOTPAIRRANKSUMMARY Quick-look plots for fixed-reference pair screening.

topCount = min(numTopReport, height(rankTable));
if topCount <= 0
  return;
end

topTable = rankTable(1:topCount, :);
labelText = compose('%d', topTable.pairSatIdxGlobal);

figure('Name', 'Fixed-reference pair screening: angle RMSE');
bar(topTable.bestAngleRmseDeg);
grid on;
xlabel('Second satellite');
ylabel('Best MS-SF-Static angle RMSE (deg)');
title('Top second-satellite candidates for the fixed reference');
set(gca, 'XTick', 1:topCount, 'XTickLabel', labelText);

figure('Name', 'Fixed-reference pair screening: gain vs reference');
bar(topTable.angleGainVsRefDeg);
grid on;
xlabel('Second satellite');
ylabel('Angle RMSE gain vs ref-only (deg)');
title('Pair gain relative to the reference-only baseline');
set(gca, 'XTick', 1:topCount, 'XTickLabel', labelText);
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

error('doaDopplerRankSecondSatDynTle:TleFileNotFound', ...
  'Unable to locate the TLE file "%s".', fileName);
end