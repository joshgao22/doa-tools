function initParam = buildDoaDopplerSfInit(model, initParam)
%BUILDDOADOPPLERSFINIT Build and validate a single-frame initializer.
% Reuses the coarse DoA search and fdRef profile refinement used by
% estimatorDoaDopplerMlePilotSfOpt.

arguments
  model (1,1) struct
  initParam = []
end

[lb, ub] = localBuildBounds(model);
if isempty(initParam)
  initParam = localBuildInit(model);
else
  initParam = localCheckInitParam(initParam, lb, ub, model);
end
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
lb = [doaLb(:); repmat(model.fdRange(1), model.numSource, 1)];
ub = [doaUb(:); repmat(model.fdRange(2), model.numSource, 1)];
if isfield(model, 'cachedBoundsReady') && ~logical(model.cachedBoundsReady)
  model.lb = lb;
  model.ub = ub;
end
end

function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build a per-source DoA search box.

if isfield(model, 'cachedBoundsReady') && logical(model.cachedBoundsReady) && ...
    isfield(model, 'doaLb') && isfield(model, 'doaUb') && ...
    ~isempty(model.doaLb) && ~isempty(model.doaUb)
  doaLb = model.doaLb;
  doaUb = model.doaUb;
  return;
end

baseRange = model.doaGrid{1}.range;
doaLb = repmat(baseRange(:, 1), 1, model.numSource);
doaUb = repmat(baseRange(:, 2), 1, model.numSource);

if isempty(model.initDoaParam) || isempty(model.initDoaHalfWidth)
  return;
end

doaCenter = reshape(model.initDoaParam, 2, []);
if size(doaCenter, 2) == 1 && model.numSource > 1
  doaCenter = repmat(doaCenter, 1, model.numSource);
end
if size(doaCenter, 2) ~= model.numSource
  error('estimatorDoaDopplerMlePilotSfOpt:InitDoaSourceMismatch', ...
    'modelOpt.initDoaParam must provide one 2x1 initializer per source.');
end

doaHalfWidth = reshape(model.initDoaHalfWidth, 2, 1);
for iSrc = 1:model.numSource
  currentCenter = doaCenter(:, iSrc);
  currentLb = currentCenter - doaHalfWidth;
  currentUb = currentCenter + doaHalfWidth;

  if strcmp(model.doaType, 'angle')
    currentCenter(1) = mod(currentCenter(1), 2 * pi);
    currentLb(1) = currentCenter(1) - doaHalfWidth(1);
    currentUb(1) = currentCenter(1) + doaHalfWidth(1);
    if currentLb(1) < baseRange(1, 1) || currentUb(1) > baseRange(1, 2)
      currentLb(1) = baseRange(1, 1);
      currentUb(1) = baseRange(1, 2);
    end
  end

  doaLb(:, iSrc) = max(doaLb(:, iSrc), currentLb);
  doaUb(:, iSrc) = min(doaUb(:, iSrc), currentUb);
end

invalidMask = doaLb > doaUb;
for iSrc = 1:model.numSource
  badIdx = invalidMask(:, iSrc);
  doaLb(badIdx, iSrc) = baseRange(badIdx, 1);
  doaUb(badIdx, iSrc) = baseRange(badIdx, 2);
end
end

function initParam = localCheckInitParam(initParam, lb, ub, model)
%LOCALCHECKINITPARAM Validate user-provided initial point.

numDoaVar = 2 * model.numSource;
numOptVar = 3 * model.numSource;

initParam = initParam(:);
if numel(initParam) ~= numOptVar
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidInitParamSize', ...
    'initParam must contain %d entries.', numOptVar);
end
if any(~isfinite(initParam))
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidInitParamValue', ...
    'initParam must contain finite values.');
end

initParam = min(max(initParam, lb), ub);
initParam = localSortOptVar(model, initParam);

doas = reshape(initParam(1:numDoaVar), 2, model.numSource);
if strcmp(model.doaType, 'angle')
  doas(1, :) = mod(doas(1, :), 2*pi);
  initParam(1:numDoaVar) = doas(:);
end
end

function initParam = localBuildInit(model)
%LOCALBUILDINIT Build a two-stage initializer.

[doaParam0, ~] = localBuildDoaInit(model);
fdRef0 = localBuildFdInit(model, doaParam0);
initParam = [doaParam0(:); fdRef0(:)];
initParam = localSortOptVar(model, initParam);
end

function [doaParam0, peakIdx] = localBuildDoaInit(model)
%LOCALBUILDDOAINIT Build the DoA initializer from MUSIC or an external hint.

if ~isempty(model.initDoaParam)
  doaParam0 = reshape(model.initDoaParam, 2, []);
  if size(doaParam0, 2) == 1 && model.numSource > 1
    doaParam0 = repmat(doaParam0, 1, model.numSource);
  end
  if size(doaParam0, 2) ~= model.numSource
    error('estimatorDoaDopplerMlePilotSfOpt:InitDoaSourceMismatch', ...
      'modelOpt.initDoaParam must provide one 2x1 initializer per source.');
  end
  if strcmp(model.doaType, 'angle')
    doaParam0(1, :) = mod(doaParam0(1, :), 2 * pi);
  end
  peakIdx = [];
  return;
end

try
  musicResult = estimatorMusic(model.array, model.wavelength, model.sampleCov, ...
    model.numSource, model.doaGrid);
catch
  musicResult = struct();
end

if isfield(musicResult, 'isResolved') && musicResult.isResolved && ...
    isfield(musicResult, 'peakIndices') && numel(musicResult.peakIndices) >= model.numSource
  peakIdx = reshape(musicResult.peakIndices(1:model.numSource), 1, []);
else
  peakIdx = localFallbackDoaSearch(model);
end

doaParam0 = localGridIndexToDoaParam(model, peakIdx);
end

function peakIdx = localFallbackDoaSearch(model)
%LOCALFALLBACKDOASEARCH Fallback DoA initializer by Bartlett fusion.

numGrid = size(model.doaGrid{1}.angleGrid, 2);
metric = zeros(1, numGrid);

for iSat = 1:model.numSat
  currentGrid = model.doaGrid{iSat}.angleGrid;
  A = steeringMatrix(model.array{iSat}, model.wavelength, currentGrid);
  R = 0.5 * (model.sampleCov{iSat} + model.sampleCov{iSat}');
  RA = R * A;

  numerator = real(sum(conj(A) .* RA, 1));
  denominator = max(real(sum(conj(A) .* A, 1)), eps);
  metric = metric + numerator ./ denominator;
end

peakResult = findDoaFromSpectrum(model.doaGrid{1}, metric, model.numSource);
if ~isfield(peakResult, 'isResolved') || ~peakResult.isResolved || ...
    numel(peakResult.peakIndices) < model.numSource
  [~, sortIdx] = sort(metric, 'descend');
  peakIdx = sort(sortIdx(1:model.numSource), 'ascend');
else
  peakIdx = reshape(peakResult.peakIndices(1:model.numSource), 1, []);
end
end

function doaParam = localGridIndexToDoaParam(model, gridIdx)
%LOCALGRIDINDEXTODOAPARAM Convert grid indices to continuous DoA parameters.

gridIdx = reshape(gridIdx, 1, []);
switch model.doaType
  case 'angle'
    doaParam = model.eciAngleGrid(:, gridIdx);

  case 'latlon'
    doaParam = model.latlonGrid(:, gridIdx);

  otherwise
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end

doaParam = localSortDoaParam(doaParam);
end

function fdRef0 = localBuildFdInit(model, doaParam0)
%LOCALBUILDFDINIT Build per-source reference-satellite Doppler initializers.
% The main parameterization is defined on the reference satellite. The
% initializer therefore searches fdRef only on that branch, then performs a
% small all-satellite profile refinement so that the multi-satellite static
% solver starts from a reference-Doppler state that is already consistent
% with the shared fdRef / deltaFd / fdSat semantics.

doaState = localBuildDoaState(model, doaParam0);
refSatIdx = model.refSatIdxLocal;

if model.initFdCount == 1
  fdGrid = mean(model.fdRange);
else
  fdGrid = linspace(model.fdRange(1), model.fdRange(2), model.initFdCount);
end

score = zeros(model.numSource, numel(fdGrid));
currentDoa = doaState.localDoaCell{refSatIdx};
A = steeringMatrix(model.array{refSatIdx}, model.wavelength, currentDoa);
if isempty(A)
  fdRef0 = repmat(mean(model.fdRange), model.numSource, 1);
  return;
end

sepSig = pinv(A) * model.rxSig{refSatIdx};
for iSrc = 1:model.numSource
  beamSig = sepSig(iSrc, :) .* exp(-1j * 2*pi * doaState.deltaFd(refSatIdx, iSrc) * model.timeSec);
  currentScore = localEvalFdGrid(model.pilotPad(iSrc, :), beamSig, model.timeSec, fdGrid);
  score(iSrc, :) = currentScore;
  [~, bestGridIdx] = max(currentScore);
  fdRef0(iSrc, 1) = localRefineFdPeakFromGrid(currentScore, fdGrid, bestGridIdx);
end

fdRef0 = localRefineFdInitProfile(model, doaState, fdRef0);
end

function fdPeak = localRefineFdPeakFromGrid(score, fdGrid, bestIdx)
%LOCALREFINEFDPEAKFROMGRID Refine one 1-D Doppler peak by parabola fitting.

fdPeak = fdGrid(bestIdx);
if numel(fdGrid) < 3 || bestIdx <= 1 || bestIdx >= numel(fdGrid)
  return;
end

yPrev = score(bestIdx - 1);
yMid = score(bestIdx);
yNext = score(bestIdx + 1);
if ~all(isfinite([yPrev, yMid, yNext]))
  return;
end

denom = yPrev - 2 * yMid + yNext;
if abs(denom) <= eps(max(abs([yPrev, yMid, yNext])))
  return;
end

delta = 0.5 * (yPrev - yNext) / denom;
delta = min(max(delta, -1), 1);
fdStep = 0.5 * (fdGrid(bestIdx + 1) - fdGrid(bestIdx - 1));
fdPeak = fdGrid(bestIdx) + delta * fdStep;
end

function fdRef0 = localRefineFdInitProfile(model, doaState, fdRefInit)
%LOCALREFINEFDINITPROFILE Refine fdRef by the shared multi-sat profile cost.
% The coarse ref-satellite search is robust, but multi-satellite static can
% still start from a slightly biased fdRef that pulls the joint optimizer to
% a worse local point. A small all-satellite profile refinement keeps the
% reference semantics fixed while nudging fdRef toward the best joint start.

fdRef0 = fdRefInit(:);
if model.numSat <= 1
  return;
end

fdRangeSpan = model.fdRange(2) - model.fdRange(1);
fdGridStep = localEstimateFdGridStep(model.fdRange, model.initFdCount);
coarseHalfWidth = max(2 * fdGridStep, max(1, 0.01 * max(fdRangeSpan, 1)));
fineHalfWidth = max(fdGridStep, max(0.25, 0.003 * max(fdRangeSpan, 1)));

for iSrc = 1:model.numSource
  coarseGrid = linspace(max(model.fdRange(1), fdRef0(iSrc) - coarseHalfWidth), ...
    min(model.fdRange(2), fdRef0(iSrc) + coarseHalfWidth), 17);
  fdRef0(iSrc) = localSearchFdRefByProfile(model, doaState, fdRef0, iSrc, coarseGrid);

  fineGrid = linspace(max(model.fdRange(1), fdRef0(iSrc) - fineHalfWidth), ...
    min(model.fdRange(2), fdRef0(iSrc) + fineHalfWidth), 17);
  fdRef0(iSrc) = localSearchFdRefByProfile(model, doaState, fdRef0, iSrc, fineGrid);
end
end

function fdBest = localSearchFdRefByProfile(model, doaState, fdRefBase, srcIdx, fdGrid)
%LOCALSEARCHFDREFBYPROFILE Search one local fdRef grid using the main cost.

fdBest = fdRefBase(srcIdx);
objBest = inf;
fdGrid = unique(fdGrid(:).');

for iFd = 1:numel(fdGrid)
  fdCand = fdRefBase;
  fdCand(srcIdx) = fdGrid(iFd);
  hyp = doaState;
  hyp.fdRef = fdCand(:);
  hyp.fd = doaState.deltaFd + reshape(fdCand, 1, []);
  [~, noiseVar, residualNorm] = localEstimateGainAndNoise(model, hyp);
  if model.useLogObjective
    objNow = sum((model.numElement * model.numSample) .* log(max(noiseVar, eps)));
  else
    objNow = residualNorm;
  end
  if isfinite(objNow) && objNow < objBest
    objBest = objNow;
    fdBest = fdGrid(iFd);
  end
end
end

function fdStep = localEstimateFdGridStep(fdRange, initFdCount)
%LOCALESTIMATEFDGRIDSTEP Estimate the coarse fd grid spacing.

if isempty(fdRange) || numel(fdRange) ~= 2 || initFdCount <= 1
  fdStep = 0;
  return;
end
fdStep = abs(fdRange(2) - fdRange(1)) / max(initFdCount - 1, 1);
if ~isfinite(fdStep)
  fdStep = 0;
end
end

function score = localEvalFdGrid(pilotPad, beamSig, timeSec, fdGrid)
%LOCALEVALFDGRID Evaluate concentrated 1D Doppler scores.

numFd = numel(fdGrid);
score = zeros(1, numFd);
beamSig = reshape(beamSig, 1, []);
pilotRow = reshape(pilotPad, 1, []);
pilotEnergy = real(pilotRow * pilotRow');

for iFd = 1:numFd
  refSig = pilotRow .* exp(1j * 2*pi * fdGrid(iFd) * timeSec);
  corrVal = beamSig * refSig';
  score(iFd) = abs(corrVal)^2 / max(pilotEnergy, eps);
end
end

function doaState = localBuildDoaState(model, doaParam)
%LOCALBUILDDOASTATE Map continuous DoA parameters to geometry quantities.

doaParam = reshape(doaParam, 2, model.numSource);

doas = doaParam;
if strcmp(model.doaType, 'angle')
  doas(1, :) = mod(doas(1, :), 2*pi);
end

switch model.doaType
  case 'angle'
    eciAngle = doas;
    eciDirection = doa2dir(eciAngle);
    localDoaCell = eciToAngleGrid(eciDirection, model.rotMat);
    if ~iscell(localDoaCell)
      localDoaCell = {localDoaCell};
    end

    deltaFd = ((eciDirection.' * model.relSatVelEci) / model.wavelength).';
    latlon = [];
    userPosEci = [];
    userVelEci = [];

  case 'latlon'
    latlon = doas;
    [userPosEci, userVelEci] = localLatlonToUserState(latlon, model.userStateRef);
    eciDirection = userPosEci ./ vecnorm(userPosEci, 2, 1);
    eciAngle = localDirectionToAngle(eciDirection);

    localDoaCell = cell(1, model.numSat);
    for iSat = 1:model.numSat
      localDoaCell{iSat} = globalToLocalDoa(userPosEci, model.satPosEci(:, iSat), model.rotMat{iSat});
    end

    refDopplerState = localComputeExactReferenceDopplerState(model, userPosEci, userVelEci);
    deltaFd = refDopplerState.deltaFd;

  otherwise
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end

localDoaArr = zeros(2, model.numSource, model.numSat);
for iSat = 1:model.numSat
  localDoaArr(:, :, iSat) = localDoaCell{iSat};
end

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = latlon;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFd = deltaFd;
if exist('refDopplerState', 'var') && isstruct(refDopplerState)
  doaState.fdRefGeom = refDopplerState.fdRefRefFrame;
  doaState.fdSatGeom = refDopplerState.fdSatRefFrame;
  doaState.deltaFdGeom = refDopplerState.deltaFdRefFrame;
end
doaState.localDoaCell = localDoaCell;
doaState.localDoaArr = localDoaArr;
end

function refDopplerState = localComputeExactReferenceDopplerState(model, userPosEci, userVelEci)
%LOCALCOMPUTEEXACTREFERENCEDOPPLERSTATE Compute exact reference Doppler state.

refDopplerState = buildReferenceDopplerState( ...
  model.sceneRef, model.satPosEci, model.satVelEci, ...
  userPosEci, userVelEci, model.wavelength, [], model.refSatIdxLocal);
end

function [pathGain, noiseVar, residualNorm, gainDiag] = localEstimateGainAndNoise(model, hyp)
%LOCALESTIMATEGAINANDNOISE Concentrate per-satellite gains and noise powers.

pathGain = zeros(model.numSat, model.numSource);
noiseVar = zeros(model.numSat, 1);
residualNorm = 0;
residualEnergySat = zeros(model.numSat, 1);

for iSat = 1:model.numSat
  dictMat = localBuildPilotDict(model, hyp, iSat);
  yVec = model.rxSig{iSat}(:);

  colEnergy = real(sum(conj(dictMat) .* dictMat, 1));
  if any(colEnergy <= eps)
    pathGain = nan(model.numSat, model.numSource);
    noiseVar = inf(model.numSat, 1);
    residualNorm = inf;
    gainDiag = struct('residualEnergySat', inf(model.numSat, 1));
    return;
  end

  gainVec = dictMat \ yVec;
  residualVec = yVec - dictMat * gainVec;
  currentVar = real(residualVec' * residualVec) / numel(yVec);
  currentVar = max(currentVar, eps);

  currentResidualEnergy = real(residualVec' * residualVec);
  pathGain(iSat, :) = gainVec.';
  noiseVar(iSat) = currentVar;
  residualEnergySat(iSat) = currentResidualEnergy;
  residualNorm = residualNorm + model.satWeight(iSat) * currentResidualEnergy;
end

gainDiag = struct();
gainDiag.residualEnergySat = residualEnergySat;
end

function dictMat = localBuildPilotDict(model, hyp, satIdx)
%LOCALBUILDPILOTDICT Build the per-satellite pilot dictionary.

numEntry = model.numElement(satIdx) * model.numSample;
dictMat = zeros(numEntry, model.numSource);
currentDoa = hyp.localDoaCell{satIdx};

for iSrc = 1:model.numSource
  [atomMat, ~] = buildPilotAtomSfKernel(model.array{satIdx}, ...
    model.wavelength, currentDoa(:, iSrc), model.pilotPad(iSrc, :), ...
    model.timeSec, hyp.fd(satIdx, iSrc));
  dictMat(:, iSrc) = atomMat(:);
end
end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2*pi); ...
         asin(max(min(direction(3, :), 1), -1))];
end

function [userPosEci, userVelEci] = localLatlonToUserState(latlon, userStateRef)
%LOCALLATLONTOUSERSTATE Convert ground [lat; lon] to motion-aware ECI states.

[userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, userStateRef);
end

function doaParam = localSortDoaParam(doaParam)
%LOCALSORTDOAPARAM Sort sources by the first and second DOA components.

if isempty(doaParam)
  return;
end

[~, sortIdx] = sortrows(doaParam.', [1, 2]);
doaParam = doaParam(:, sortIdx);
end

function optVar = localSortOptVar(model, optVar)
%LOCALSORTOPTVAR Sort source-related variables in one optimization vector.

optVar = optVar(:);
if isempty(optVar)
  return;
end

numDoaVar = 2 * model.numSource;
doaParam = reshape(optVar(1:numDoaVar), 2, model.numSource);
fdRef = optVar(numDoaVar + 1:end);
[~, sortIdx] = sortrows(doaParam.', [1, 2]);
doaParam = doaParam(:, sortIdx);
fdRef = fdRef(sortIdx);
optVar = [doaParam(:); fdRef(:)];
end
