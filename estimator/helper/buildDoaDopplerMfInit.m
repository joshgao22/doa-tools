function [initParam, initDiag, model] = buildDoaDopplerMfInit(model, initParam)
%BUILDDOADOPPLERMFINIT Build and validate a multi-frame initializer.
% This helper reuses the reference-satellite fd-line initialization,
% external DoA seeding, and alias-prior finalization used by
% estimatorDoaDopplerMlePilotMfOpt.

arguments
  model (1,1) struct
  initParam = []
end

[lb, ub] = localBuildBounds(model);

if isempty(initParam)
  [initParam, initDiag] = localBuildInit(model);
else
  initParam = localCheckInitParam(initParam, lb, ub, model);
  initDiag = localBuildUserInitDiag(model, initParam);
end

[model, initDiag] = localFinalizeFdAliasPrior(model, initParam, initDiag);
end

function uniqueVec = localUniqueStableTol(valVec, tol)
%LOCALUNIQUESTABLETOL Stable unique with tolerance for warm-start seeds.

if nargin < 2 || isempty(tol)
  tol = 0;
end
valVec = reshape(valVec, 1, []);
if isempty(valVec)
  uniqueVec = valVec;
  return;
end

keepMask = false(size(valVec));
uniqueList = zeros(1, 0);
for iVal = 1:numel(valVec)
  if isempty(uniqueList) || all(abs(valVec(iVal) - uniqueList) > tol)
    uniqueList(end + 1) = valVec(iVal); %#ok<AGROW>
    keepMask(iVal) = true;
  end
end
uniqueVec = valVec(keepMask);
end

function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build the DoA search box used by all MF refinements.
% The local DoA refinement must be enforced by the optimizer bounds rather
% than only stored as a nominal debug setting. Precomputing the DoA box here
% makes every path reuse the same local basin definition.

if isfield(model, 'cachedBoundsReady') && logical(model.cachedBoundsReady) && ...
    isfield(model, 'doaLb') && isfield(model, 'doaUb') && ...
    ~isempty(model.doaLb) && ~isempty(model.doaUb)
  doaLb = model.doaLb;
  doaUb = model.doaUb;
  return;
end

baseRange = model.doaGrid{1}.range;
doaLb = baseRange(:, 1);
doaUb = baseRange(:, 2);

if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  return;
end

if isempty(model.initDoaParam) || isempty(model.initDoaHalfWidth)
  return;
end

doaCenter = model.initDoaParam(:);
doaHalfWidth = model.initDoaHalfWidth(:);

doaLbLocal = doaCenter - doaHalfWidth;
doaUbLocal = doaCenter + doaHalfWidth;

if strcmp(model.doaType, 'angle')
  doaCenter(1) = mod(doaCenter(1), 2 * pi);
  doaLbLocal(1) = doaCenter(1) - doaHalfWidth(1);
  doaUbLocal(1) = doaCenter(1) + doaHalfWidth(1);
  if doaLbLocal(1) < baseRange(1, 1) || doaUbLocal(1) > baseRange(1, 2)
    doaLbLocal(1) = baseRange(1, 1);
    doaUbLocal(1) = baseRange(1, 2);
  end
end

doaLb = max(doaLb, doaLbLocal);
doaUb = min(doaUb, doaUbLocal);
invalidMask = doaLb > doaUb;
doaLb(invalidMask) = baseRange(invalidMask, 1);
doaUb(invalidMask) = baseRange(invalidMask, 2);
end

function [lb, ub] = localBuildBounds(model)
%LOCALBUILDBOUNDS Build box constraints for fmincon.

if isfield(model, 'lb') && isfield(model, 'ub') && ...
    ~isempty(model.lb) && ~isempty(model.ub)
  lb = model.lb;
  ub = model.ub;
  return;
end

[doaLb, doaUb] = localBuildDoaBounds(model);

lb = [doaLb; model.fdRange(1)];
ub = [doaUb; model.fdRange(2)];

if strcmp(model.fdRateMode, 'unknown')
  lb = [lb; model.fdRateRange(1)];
  ub = [ub; model.fdRateRange(2)];
end
if isfield(model, 'cachedBoundsReady') && ~logical(model.cachedBoundsReady)
  model.lb = lb;
  model.ub = ub;
end
end

function initParam = localCheckInitParam(initParam, lb, ub, model)
%LOCALCHECKINITPARAM Validate one user-provided initial point.

initParam = initParam(:);
numVar = localGetNumOptVar(model);

if numel(initParam) ~= numVar
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidInitParamSize', ...
    'initParam must contain %d entries for fdRateMode = ''%s''.', ...
    numVar, model.fdRateMode);
end
if any(~isfinite(initParam))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidInitParamValue', ...
    'initParam must contain finite values.');
end

initParam = min(max(initParam, lb), ub);
if strcmp(model.doaType, 'angle')
  initParam(1) = mod(initParam(1), 2*pi);
end
end

function [model, initDiag] = localFinalizeFdAliasPrior(model, initParam, initDiag)
%LOCALFINALIZEFDALIASPRIOR Finalize the per-satellite Doppler prior.

if ~model.enableFdAliasUnwrap
  return;
end

if ~isempty(model.fdSatPriorHz)
  if numel(model.fdSatPriorHz) ~= model.numSat
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdSatPriorHzSize', ...
      'modelOpt.fdSatPriorHz must contain one Doppler prior per satellite.');
  end
  model.fdSatPriorHz = reshape(model.fdSatPriorHz, [], 1);
  priorSource = 'user';
else
  hypInit = buildDoaDopplerMfHypothesis(model, initParam);
  model.fdSatPriorHz = hypInit.fdSat(:);
  priorSource = 'initHyp';
end

if isstruct(initDiag)
  if ~isfield(initDiag, 'fdLine') || ~isstruct(initDiag.fdLine)
    initDiag.fdLine = struct();
  end
  initDiag.fdLine.fdAliasPriorSource = priorSource;
  initDiag.fdLine.fdAliasStepHz = model.fdAliasStepHz;
  initDiag.fdLine.fdSatPriorHz = model.fdSatPriorHz(:);
end
end

function [initParam, initDiag] = localBuildInit(model)
%LOCALBUILDINIT Build a two-stage initializer.
% The mainline first forms a DoA coarse point, then estimates the
% reference-satellite fd-line, and finally lets the multi-frame likelihood
% refine the result inside the configured search box.

[doaParam0, doaDiag] = localBuildDoaInit(model);
[fdRef0, fdRate0, fdDiag] = localBuildFdLineInit(model, doaParam0);

switch model.fdRateMode
  case 'unknown'
    initParam = [doaParam0(:); fdRef0; fdRate0];
  case {'known', 'zero'}
    initParam = [doaParam0(:); fdRef0];
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end

[lb, ub] = localBuildBounds(model);
initParam = min(max(initParam, lb), ub);

[doaLb, doaUb] = localBuildDoaBounds(model);

initDiag = struct();
initDiag.source = 'auto';
initDiag.initParam = initParam(:);
initDiag.doa = doaDiag;
initDiag.fdLine = fdDiag;
initDiag.doaLb = doaLb;
initDiag.doaUb = doaUb;
end

function initDiag = localBuildUserInitDiag(model, initParam)
%LOCALBUILDUSERINITDIAG Build lightweight debug info for a user init point.

fdRateVal = nan;
if strcmp(model.fdRateMode, 'unknown') && numel(initParam) >= 4
  fdRateVal = initParam(4);
elseif strcmp(model.fdRateMode, 'known')
  fdRateVal = model.fdRateKnown;
elseif strcmp(model.fdRateMode, 'zero')
  fdRateVal = 0;
end

[doaLb, doaUb] = localBuildDoaBounds(model);

initDiag = struct();
initDiag.source = 'user';
initDiag.initParam = initParam(:);
initDiag.doa = struct('selectedDoaParam', initParam(1:2), 'selectedGridIdx', nan);
initDiag.fdLine = struct('fdRefInit', initParam(3), 'fdRateInit', fdRateVal);
initDiag.doaLb = doaLb;
initDiag.doaUb = doaUb;
end

function [doaParam0, doaDiag] = localBuildDoaInit(model)
%LOCALBUILDDOAINIT Build DoA initializer from frame-aggregated MUSIC.
% If model.initDoaParam is provided, it takes priority over the automatic
% MUSIC/Bartlett route so that the multi-frame solver can reuse a reliable
% coarse DoA estimate from an upstream single-frame stage.

[musicSpec, bartlettMetric] = localBuildDoaInitMetric(model);

try
  musicResult = estimatorMusic(model.array, model.wavelength, model.sampleCovAgg, 1, model.doaGrid);
catch
  musicResult = struct();
end

if ~isempty(model.initDoaParam)
  doaParam0 = model.initDoaParam(:);
  peakIdx = nan;
  initSource = 'userDoa';
elseif isfield(musicResult, 'isResolved') && musicResult.isResolved && ...
    isfield(musicResult, 'peakIndices') && ~isempty(musicResult.peakIndices)
  peakIdx = musicResult.peakIndices(1);
  initSource = 'music';
  doaParam0 = localGridIndexToDoaParam(model, peakIdx);
else
  peakIdx = localFallbackDoaSearch(model, bartlettMetric);
  initSource = 'bartlettFallback';
  doaParam0 = localGridIndexToDoaParam(model, peakIdx);
end

if strcmp(model.doaType, 'angle')
  doaParam0(1) = mod(doaParam0(1), 2 * pi);
end

doaDiag = struct();
doaDiag.initSource = initSource;
doaDiag.selectedGridIdx = peakIdx;
doaDiag.selectedDoaParam = doaParam0(:);
doaDiag.musicResolved = isfield(musicResult, 'isResolved') && logical(musicResult.isResolved);
doaDiag.musicTop = localBuildTopGridDiag(model, musicSpec, 3);
doaDiag.bartlettTop = localBuildTopGridDiag(model, bartlettMetric, 3);
end

function peakIdx = localFallbackDoaSearch(model, metric)
%LOCALFALLBACKDOASEARCH Fallback DoA initializer by Bartlett fusion.

if nargin < 2 || isempty(metric)
  [~, metric] = localBuildDoaInitMetric(model);
end

peakResult = findDoaFromSpectrum(model.doaGrid{1}, metric, 1);
if ~isfield(peakResult, 'isResolved') || ~peakResult.isResolved || isempty(peakResult.peakIndices)
  [~, peakIdx] = max(metric);
else
  peakIdx = peakResult.peakIndices(1);
end
end

function [musicSpec, bartlettMetric] = localBuildDoaInitMetric(model)
%LOCALBUILDDOAINITMETRIC Build MUSIC and Bartlett spectra on the search grid.

numGrid = size(model.doaGrid{1}.angleGrid, 2);
musicInv = zeros(1, numGrid);
bartlettMetric = zeros(1, numGrid);

for iSat = 1:model.numSat
  currentGrid = model.doaGrid{iSat}.angleGrid;
  A = steeringMatrix(model.array{iSat}, model.wavelength, currentGrid);
  R = 0.5 * (model.sampleCovAgg{iSat} + model.sampleCovAgg{iSat}');

  [eigVec, eigVal] = eig(R, 'vector');
  [~, sortIdx] = sort(real(eigVal), 'ascend');
  noiseSub = eigVec(:, sortIdx(1:end-1));
  proj = noiseSub' * A;
  musicInv = musicInv + real(sum(conj(proj) .* proj, 1));

  RA = R * A;
  numerator = real(sum(conj(A) .* RA, 1));
  denominator = max(real(sum(conj(A) .* A, 1)), eps);
  bartlettMetric = bartlettMetric + numerator ./ denominator;
end

musicSpec = 1 ./ max(musicInv, eps);
end

function topDiag = localBuildTopGridDiag(model, metric, topCount)
%LOCALBUILDTOPGRIDDIAG Pack the strongest grid candidates for debugging.

if nargin < 3 || isempty(topCount)
  topCount = 5;
end

metric = reshape(real(metric), 1, []);
numKeep = min(topCount, numel(metric));
[metricSort, idxSort] = sort(metric, 'descend');
idxKeep = idxSort(1:numKeep);

topDiag = struct();
topDiag.gridIdx = idxKeep(:);
topDiag.metric = metricSort(1:numKeep).';

topDoa = zeros(2, numKeep);
for iGrid = 1:numKeep
  topDoa(:, iGrid) = localGridIndexToDoaParam(model, idxKeep(iGrid));
end
topDiag.doaParam = topDoa;
end

function doaParam = localGridIndexToDoaParam(model, gridIdx)
%LOCALGRIDINDEXTODOAPARAM Convert one grid index to continuous DoA parameters.

switch model.doaType
  case 'angle'
    doaParam = model.eciAngleGrid(:, gridIdx);

  case 'latlon'
    doaParam = model.latlonGrid(:, gridIdx);

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end
end

function [fdRef0, fdRate0, fdDiag] = localBuildFdLineInit(model, doaParam0)
%LOCALBUILDFDLINEINIT Build [fdRef, fdRate] initializer.
%
% This initializer intentionally uses only the reference-satellite branch.
% The DoA coarse estimate is still shared, but the Doppler line is seeded
% from a single-satellite per-frame search to avoid multi-satellite coupling
% errors during initialization.

doaState = buildDoaDopplerMfHypothesis(model, doaParam0, "doa");
refSatIdx = localResolveRefSatIdx(model);

if model.initFdCount == 1
  fdGrid = mean(model.fdRange);
else
  fdGrid = linspace(model.fdRange(1), model.fdRange(2), model.initFdCount);
end

fdFrame = nan(1, model.numFrame);
frameWeight = zeros(1, model.numFrame);
frameBestScore = nan(1, model.numFrame);
framePeakToMedian = nan(1, model.numFrame);
deltaFdRefInit = nan(1, model.numFrame);

for iFrame = 1:model.numFrame
  if ~model.frameMask(refSatIdx, iFrame)
    continue;
  end

  currentDoa = doaState.localDoaArr(:, refSatIdx, iFrame);
  aVec = steeringMatrix(model.array{refSatIdx}, model.wavelength, currentDoa);
  if isempty(aVec)
    continue;
  end

  beamSig = (aVec' * model.rxSig{iFrame}{refSatIdx}) / max(real(aVec' * aVec), eps);
  deltaFd = doaState.deltaFdSeq(refSatIdx, iFrame);
  deltaFdRefInit(iFrame) = deltaFd;

  currentScore = localEvalFdGrid(model.pilotPadCell{iFrame}, beamSig, ...
    model.timeSecCell{iFrame}, fdGrid, deltaFd);

  validMask = isfinite(currentScore);
  if ~any(validMask)
    continue;
  end

  [bestScore, bestIdx] = max(currentScore(validMask));
  validIdx = find(validMask);
  bestGridIdx = validIdx(bestIdx);
  fdFrameRefined = localRefineFdPeakFromGrid(currentScore, fdGrid, bestGridIdx);

  fdFrame(iFrame) = fdFrameRefined;
  frameWeight(iFrame) = 1;
  frameBestScore(iFrame) = bestScore;

  medianScore = median(currentScore(validMask));
  if isfinite(medianScore) && medianScore > 0
    framePeakToMedian(iFrame) = bestScore / medianScore;
  else
    framePeakToMedian(iFrame) = NaN;
  end
end

[fdRef0, fdRate0] = localFitFdLineByMode(model, fdFrame, frameWeight);

fdRef0 = min(max(fdRef0, model.fdRange(1)), model.fdRange(2));
switch model.fdRateMode
  case 'known'
    fdRate0 = model.fdRateKnown;
  case 'zero'
    fdRate0 = 0;
  otherwise
    fdRate0 = min(max(fdRate0, model.fdRateRange(1)), model.fdRateRange(2));
    [fdRef0, fdRate0, fdRefineDiag] = localRefineRefFdLineInit( ...
      model, doaState, refSatIdx, fdRef0, fdRate0, fdGrid);
end

fdDiag = struct();
fdDiag.method = 'singleSatFrameLine';
fdDiag.refSatIdxLocal = refSatIdx;
fdDiag.doaParamInit = doaParam0(:);
fdDiag.fdGridRange = model.fdRange(:);
fdDiag.fdGridCount = numel(fdGrid);
fdDiag.fdGridStepHz = localGetGridStep(fdGrid);
fdDiag.fdFrame = fdFrame(:);
fdDiag.frameWeight = frameWeight(:);
fdDiag.frameBestScore = frameBestScore(:);
fdDiag.framePeakToMedian = framePeakToMedian(:);
fdDiag.deltaFdRefInit = deltaFdRefInit(:);
fdDiag.fdRefInit = fdRef0;
fdDiag.fdRateInit = fdRate0;
fdDiag.fdRefTruth = localGetStructField(model.debugTruth, 'fdRef', NaN);
fdDiag.fdRateTruth = localGetStructField(model.debugTruth, 'fdRate', NaN);
if exist('fdRefineDiag', 'var')
  fdDiag.refine = fdRefineDiag;
else
  fdDiag.refine = struct();
end
end

function [fdRefBest, fdRateBest, refineDiag] = localRefineRefFdLineInit(model, doaState, refSatIdx, fdRefInit, fdRateInit, fdGrid)
%LOCALREFINEREFFDLINEINIT Refine the reference-satellite fd line.
% The coarse frame-wise fd grid is good at seeding fdRef but can collapse
% fdRate to zero when every frame picks the same coarse fd bin. This helper
% therefore performs a reference-satellite-only continuous-phase search on
% [fdRef, fdRate] in two stages: a broad rate scan to escape the zero-rate
% seed, followed by a compact local refinement around the best coarse pair.

fdRefBest = fdRefInit;
fdRateBest = fdRateInit;
refineDiag = struct('isEnabled', false);

if ~strcmp(model.phaseMode, 'continuous') || model.numFrame <= 1
  return;
end
if model.numSat < 1 || refSatIdx < 1 || refSatIdx > model.numSat
  return;
end
if isempty(fdGrid) || numel(fdGrid) < 2
  fdGridStep = 0;
else
  fdGridStep = abs(fdGrid(min(end, 2)) - fdGrid(1));
end
if ~(isfinite(fdGridStep) && fdGridStep > 0)
  fdGridStep = 0.01 * max(model.fdRange(2) - model.fdRange(1), 1);
end

beamSigCell = cell(1, model.numFrame);
pilotCell = cell(1, model.numFrame);
timeGlobalCell = cell(1, model.numFrame);
frameMaskUse = false(1, model.numFrame);
energyVec = zeros(1, model.numFrame);
for iFrame = 1:model.numFrame
  if ~model.frameMask(refSatIdx, iFrame)
    continue;
  end
  currentDoa = doaState.localDoaArr(:, refSatIdx, iFrame);
  aVec = steeringMatrix(model.array{refSatIdx}, model.wavelength, currentDoa);
  if isempty(aVec)
    continue;
  end
  beamSig = (aVec' * model.rxSig{iFrame}{refSatIdx}) / max(real(aVec' * aVec), eps);
  pilotNow = reshape(model.pilotPadCell{iFrame}, 1, []);
  timeGlobal = reshape(model.timeOffsetSec(iFrame) + model.timeSecCell{iFrame}, 1, []);
  beamSigCell{iFrame} = reshape(beamSig, 1, []);
  pilotCell{iFrame} = pilotNow;
  timeGlobalCell{iFrame} = timeGlobal;
  energyVec(iFrame) = max(real(pilotNow * pilotNow'), eps);
  frameMaskUse(iFrame) = true;
end

if nnz(frameMaskUse) < 2
  return;
end

fdRefHalfWidthCoarse = max(1.5 * fdGridStep, 1);
fdRefGridCoarse = linspace(max(model.fdRange(1), fdRefInit - fdRefHalfWidthCoarse), ...
  min(model.fdRange(2), fdRefInit + fdRefHalfWidthCoarse), 17);
fdRefGridCoarse = localUniqueStableTol(fdRefGridCoarse, 1e-12);
fdRateGridCoarse = localBuildRefFdRateInitGrid(model, fdRateInit);
if isempty(fdRateGridCoarse)
  return;
end

[fdRefBest, fdRateBest, scoreBest] = localSearchRefFdLineGrid( ...
  beamSigCell, pilotCell, timeGlobalCell, energyVec, frameMaskUse, ...
  fdRefGridCoarse, fdRateGridCoarse, fdRefBest, fdRateBest);

fdRefHalfWidthFine = max(0.5 * fdGridStep, 0.25);
fdRefGridFine = linspace(max(model.fdRange(1), fdRefBest - fdRefHalfWidthFine), ...
  min(model.fdRange(2), fdRefBest + fdRefHalfWidthFine), 17);
fdRefGridFine = localUniqueStableTol(fdRefGridFine, 1e-12);
fdRateHalfWidthFine = localGetRefFdRateFineHalfWidth(model, fdRateGridCoarse, fdRateBest);
fdRateGridFine = linspace(max(model.fdRateRange(1), fdRateBest - fdRateHalfWidthFine), ...
  min(model.fdRateRange(2), fdRateBest + fdRateHalfWidthFine), 21);
fdRateGridFine = localUniqueStableTol(fdRateGridFine, 1e-9);

[fdRefBest, fdRateBest, scoreBest] = localSearchRefFdLineGrid( ...
  beamSigCell, pilotCell, timeGlobalCell, energyVec, frameMaskUse, ...
  fdRefGridFine, fdRateGridFine, fdRefBest, fdRateBest, scoreBest);

fdRefBest = min(max(fdRefBest, model.fdRange(1)), model.fdRange(2));
fdRateBest = min(max(fdRateBest, model.fdRateRange(1)), model.fdRateRange(2));

refineDiag.isEnabled = true;
refineDiag.fdRefInit = fdRefInit;
refineDiag.fdRateInit = fdRateInit;
refineDiag.fdRefBest = fdRefBest;
refineDiag.fdRateBest = fdRateBest;
refineDiag.fdRefGridCoarse = fdRefGridCoarse(:);
refineDiag.fdRateGridCoarse = fdRateGridCoarse(:);
refineDiag.fdRefGridFine = fdRefGridFine(:);
refineDiag.fdRateGridFine = fdRateGridFine(:);
refineDiag.refSatIdxLocal = refSatIdx;
refineDiag.numValidFrame = nnz(frameMaskUse);
refineDiag.bestScore = scoreBest;
end

function fdRateGrid = localBuildRefFdRateInitGrid(model, fdRateSeed)
%LOCALBUILDREFFDRATEINITGRID Build a broad-plus-local fdRate grid.
% The reference-satellite refiner must be able to escape a collapsed zero
% seed. We therefore combine a broad coarse scan over the full fdRate range
% with a denser local grid around the current seed. This stays reference-only
% and does not reintroduce the failed multi-satellite joint initializer.

fdRateGrid = [];
if isempty(model.fdRateRange) || numel(model.fdRateRange) ~= 2
  return;
end
fdMin = model.fdRateRange(1);
fdMax = model.fdRateRange(2);
if ~(isfinite(fdMin) && isfinite(fdMax) && fdMin <= fdMax)
  return;
end

rangeSpan = fdMax - fdMin;
if rangeSpan <= 0
  fdRateGrid = fdMin;
  return;
end

if ~isfinite(fdRateSeed)
  fdRateSeed = 0.5 * (fdMin + fdMax);
end

localHalfWidth = min(0.2 * rangeSpan, max(50, 0.08 * max(abs(fdRateSeed), rangeSpan)));
localGrid = linspace(max(fdMin, fdRateSeed - localHalfWidth), ...
  min(fdMax, fdRateSeed + localHalfWidth), 11);
coarseCount = 17;
coarseGrid = linspace(fdMin, fdMax, coarseCount);
fdRateGrid = [fdRateSeed, localGrid, coarseGrid];
fdRateGrid = localUniqueStableTol(min(max(fdRateGrid, fdMin), fdMax), 1e-9);
end

function [fdRefBest, fdRateBest, scoreBest] = localSearchRefFdLineGrid( ...
  beamSigCell, pilotCell, timeGlobalCell, energyVec, frameMaskUse, ...
  fdRefGrid, fdRateGrid, fdRefBest, fdRateBest, scoreBest)
%LOCALSEARCHREFFDLINEGRID Search one [fdRef, fdRate] grid on ref-sat data.

if nargin < 10 || isempty(scoreBest)
  scoreBest = -Inf;
end

for iRate = 1:numel(fdRateGrid)
  rateVal = fdRateGrid(iRate);
  for iFd = 1:numel(fdRefGrid)
    fdVal = fdRefGrid(iFd);
    scoreNow = 0;
    for iFrame = 1:numel(frameMaskUse)
      if ~frameMaskUse(iFrame)
        continue;
      end
      phaseVec = exp(1j * 2 * pi * (fdVal * timeGlobalCell{iFrame} + ...
        0.5 * rateVal * (timeGlobalCell{iFrame} .^ 2)));
      refSig = pilotCell{iFrame} .* phaseVec;
      corrVal = beamSigCell{iFrame} * refSig';
      scoreNow = scoreNow + abs(corrVal)^2 / energyVec(iFrame);
    end
    if scoreNow > scoreBest
      scoreBest = scoreNow;
      fdRefBest = fdVal;
      fdRateBest = rateVal;
    end
  end
end
end

function fdRateHalfWidth = localGetRefFdRateFineHalfWidth(model, fdRateGridCoarse, fdRateBest)
%LOCALGETREFFDRATEFINEHALFWIDTH Build the local fdRate fine-search radius.

if nargin < 2 || isempty(fdRateGridCoarse)
  coarseStep = 0.1 * max(model.fdRateRange(2) - model.fdRateRange(1), 1);
else
  coarseStep = median(abs(diff(sort(fdRateGridCoarse(:)))));
end
if ~(isfinite(coarseStep) && coarseStep > 0)
  coarseStep = 0.1 * max(model.fdRateRange(2) - model.fdRateRange(1), 1);
end
fdRateHalfWidth = max(0.5 * coarseStep, max(20, 0.05 * max(abs(fdRateBest), 1)));
fdRateHalfWidth = min(fdRateHalfWidth, 0.15 * max(model.fdRateRange(2) - model.fdRateRange(1), 1));
end

function refSatIdx = localResolveRefSatIdx(model)
%LOCALRESOLVEREFSATIDX Resolve the local reference-satellite index.

refSatIdx = [];

if isfield(model, 'refSatIdxLocal') && ~isempty(model.refSatIdxLocal) && ...
    isscalar(model.refSatIdxLocal) && isfinite(model.refSatIdxLocal)
  refSatIdx = model.refSatIdxLocal;
end

if isempty(refSatIdx) && isfield(model, 'sceneCell') && ~isempty(model.sceneCell) && ...
    isfield(model, 'refFrameIdx') && isscalar(model.refFrameIdx) && ...
    model.refFrameIdx >= 1 && model.refFrameIdx <= numel(model.sceneCell)
  refScene = model.sceneCell{model.refFrameIdx};
  if isstruct(refScene) && isfield(refScene, 'ref') && isstruct(refScene.ref) && ...
      isfield(refScene.ref, 'weight') && ~isempty(refScene.ref.weight)
    refWeight = refScene.ref.weight(:);
    refSatIdx = find(refWeight > 0.5, 1, 'first');
  end
end

if isempty(refSatIdx)
  refSatIdx = localGetStructField(model.debugTruth, 'refSatIdxLocal', []);
end

if isempty(refSatIdx) || ~isscalar(refSatIdx) || ~isfinite(refSatIdx)
  refSatIdx = 1;
end

refSatIdx = max(1, min(model.numSat, round(refSatIdx)));
end

function score = localEvalFdGrid(pilotPad, beamSig, timeSec, fdGrid, deltaFd)
%LOCALEVALFDGRID Evaluate a concentrated 1D fdRef grid score.

pilotRow = reshape(pilotPad, 1, []);
beamRow = reshape(beamSig, 1, []);
pilotEnergy = real(pilotRow * pilotRow');

score = nan(1, numel(fdGrid));
for iFd = 1:numel(fdGrid)
  currentFd = fdGrid(iFd) + deltaFd;
  refSig = pilotRow .* exp(1j * 2*pi * currentFd * timeSec);
  corrVal = beamRow * refSig';
  score(iFd) = abs(corrVal)^2 / max(pilotEnergy, eps);
end
end

function [fdRef, fdRate] = localFitFdLineByMode(model, fdFrame, frameWeight)
%LOCALFITFDLINEBYMODE Fit fdFrame ~ fdRef + fdRate * timeOffsetSec.

validMask = isfinite(fdFrame) & frameWeight > 0;
validCount = nnz(validMask);

if validCount == 0
  fdRef = mean(model.fdRange);
  switch model.fdRateMode
    case 'known'
      fdRate = model.fdRateKnown;
    case 'zero'
      fdRate = 0;
    otherwise
      fdRate = min(max(0, model.fdRateRange(1)), model.fdRateRange(2));
  end
  return;
end

timeVec = reshape(model.timeOffsetSec(validMask), [], 1);
fdVec = reshape(fdFrame(validMask), [], 1);
weightVec = reshape(frameWeight(validMask), [], 1);

switch model.fdRateMode
  case 'known'
    fdRate = model.fdRateKnown;
    fdRef = sum(weightVec .* (fdVec - fdRate * timeVec)) / max(sum(weightVec), eps);
    return;

  case 'zero'
    fdRate = 0;
    fdRef = sum(weightVec .* fdVec) / max(sum(weightVec), eps);
    return;
end

if validCount == 1 || max(timeVec) - min(timeVec) <= 0
  fdRef = fdVec(1);
  fdRate = min(max(0, model.fdRateRange(1)), model.fdRateRange(2));
  return;
end

designMat = [ones(validCount, 1), timeVec];
weightMat = diag(weightVec);
coef = (designMat' * weightMat * designMat) \ (designMat' * weightMat * fdVec);
fdRef = coef(1);
fdRate = coef(2);
end

function gridStep = localGetGridStep(gridVec)
%LOCALGETGRIDSTEP Get the uniform grid spacing when available.

if numel(gridVec) < 2
  gridStep = NaN;
else
  stepVec = diff(gridVec(:));
  if any(abs(stepVec - stepVec(1)) > 1e-12 * max(abs(stepVec(1)), 1))
    gridStep = NaN;
  else
    gridStep = stepVec(1);
  end
end
end

function fdPeak = localRefineFdPeakFromGrid(scoreVec, fdGrid, bestGridIdx)
%LOCALREFINEFDPEAKFROMGRID Refine one coarse fd-grid peak by parabola fit.

fdPeak = fdGrid(bestGridIdx);
if bestGridIdx <= 1 || bestGridIdx >= numel(fdGrid)
  return;
end

scoreTriplet = real(scoreVec(bestGridIdx + (-1:1)));
if any(~isfinite(scoreTriplet)) || any(scoreTriplet <= 0)
  return;
end

stepLeft = fdGrid(bestGridIdx) - fdGrid(bestGridIdx - 1);
stepRight = fdGrid(bestGridIdx + 1) - fdGrid(bestGridIdx);
if ~(isfinite(stepLeft) && isfinite(stepRight) && abs(stepLeft - stepRight) <= 1e-9 * max(abs(stepLeft), 1))
  return;
end

gridStep = 0.5 * (stepLeft + stepRight);
logScore = log(max(scoreTriplet, realmin));
denom = logScore(1) - 2 * logScore(2) + logScore(3);
if ~isfinite(denom) || abs(denom) < eps
  return;
end

delta = 0.5 * (logScore(1) - logScore(3)) / denom;
delta = min(max(delta, -1), 1);
fdPeak = fdGrid(bestGridIdx) + delta * gridStep;
end

function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct/object field with a default value.

fieldValue = defaultValue;

if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end

function numVar = localGetNumOptVar(model)
%LOCALGETNUMOPTVAR Get the current optimization vector size.

switch model.fdRateMode
  case 'unknown'
    numVar = 4;
  case {'known', 'zero'}
    numVar = 3;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end
end
