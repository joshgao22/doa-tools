function hyp = buildDoaDopplerMfHypothesis(model, param, buildMode)
%BUILDDOADOPPLERMFHYPOTHESIS Map MF parameters to geometry and Doppler.
% This helper is the single estimator-side construction path for the
% dynamic DoA state and the reference-satellite Doppler hypothesis. Use
% buildMode="doa" when only the DoA-dependent state is needed during
% initialization.

arguments
  model (1,1) struct
  param (:,1) double
  buildMode (1,1) string = "optvar"
end

buildMode = lower(strtrim(buildMode));
switch buildMode
  case "optvar"
    hyp = localBuildOptVarHypothesis(model, param);
  case "doa"
    hyp = localBuildDoaState(model, param);
    hyp.doaParam = reshape(param(1:2), [], 1);
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidHypothesisBuildMode', ...
      'Unsupported MF hypothesis build mode: %s.', char(buildMode));
end
end

function hyp = localBuildOptVarHypothesis(model, optVar)
%LOCALBUILDOPTVARHYPOTHESIS Map one optimization point to full hypothesis.

optVar = optVar(:);
numVar = localGetNumOptVar(model);
if numel(optVar) ~= numVar
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidOptVarSize', ...
    'The optimization vector must contain %d entries for fdRateMode = ''%s''.', ...
    numVar, model.fdRateMode);
end

doaParam = optVar(1:2);
fdRef = optVar(3);

switch model.fdRateMode
  case 'unknown'
    fdRate = optVar(4);
  case 'known'
    fdRate = model.fdRateKnown;
  case 'zero'
    fdRate = 0;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end

doaState = localBuildDoaState(model, doaParam);
fdSat = doaState.deltaFdRef + fdRef;
fdRateSat = doaState.deltaFdRate + fdRate;

hyp = doaState;
hyp.doaParam = doaParam(:);
hyp.fdRef = fdRef;
hyp.fdRate = fdRate;
hyp.fdSat = fdSat;
hyp.fdRateSat = fdRateSat;
end

function doaState = localBuildDoaState(model, doaParam)
%LOCALBUILDDOASTATE Map continuous DoA parameters to dynamic geometry.

doaParam = reshape(doaParam, 2, 1);

switch model.doaType
  case 'angle'
    doaState = localBuildAngleState(model, doaParam);
  case 'latlon'
    doaState = localBuildLatlonState(model, doaParam);
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end
end

function doaState = localBuildAngleState(model, doaParam)
%LOCALBUILDANGLESTATE Build frame-wise states for angle parameterization.

eciAngle = doaParam;
eciAngle(1) = mod(eciAngle(1), 2*pi);
eciDirection = doa2dir(eciAngle);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refRotMat = model.rotMatCell{model.steeringRefFrameIdx};
    localDoaCell = eciToAngleGrid(eciDirection, refRotMat);
    if ~iscell(localDoaCell)
      localDoaCell = {localDoaCell};
    end
    for iFrame = 1:model.numFrame
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end
  case 'framewise'
    for iFrame = 1:model.numFrame
      localDoaCell = eciToAngleGrid(eciDirection, model.rotMatCell{iFrame});
      if ~iscell(localDoaCell)
        localDoaCell = {localDoaCell};
      end
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

deltaFdSeq = zeros(model.numSat, model.numFrame);
for iFrame = 1:model.numFrame
  relSatVel = model.satVelEci(:, :, iFrame) - model.refVelEci(:, iFrame);
  deltaFdSeq(:, iFrame) = ((eciDirection.' * relSatVel) / model.wavelength).';
end

[deltaFdRef, deltaFdRate] = localFitSatFdLine(model.timeOffsetSec, deltaFdSeq);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = [];
doaState.userPosEci = [];
doaState.userVelEci = [];
doaState.deltaFdSeq = deltaFdSeq;
doaState.deltaFdRef = deltaFdRef;
doaState.deltaFdRate = deltaFdRate;
doaState.localDoaArr = localDoaArr;
end

function doaState = localBuildLatlonState(model, doaParam)
%LOCALBUILDLATLONSTATE Build frame-wise states for latlon parameterization.

latlon = doaParam;
[userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, model.userStateRef);

eciDirection = userPosEci ./ vecnorm(userPosEci, 2, 1);
eciAngle = localDirectionToAngle(eciDirection);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refFrameIdx = model.steeringRefFrameIdx;
    refUserPos = userPosEci(:, refFrameIdx);
    for iSat = 1:model.numSat
      localDoaArr(:, iSat, :) = repmat( ...
        globalToLocalDoa(refUserPos, model.satPosEci(:, iSat, refFrameIdx), ...
        model.rotMatCell{refFrameIdx}{iSat}), ...
        1, 1, model.numFrame);
    end
  case 'framewise'
    for iFrame = 1:model.numFrame
      currentUserPos = userPosEci(:, iFrame);
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = globalToLocalDoa( ...
          currentUserPos, model.satPosEci(:, iSat, iFrame), ...
          model.rotMatCell{iFrame}{iSat});
      end
    end
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

refDopplerState = buildReferenceDopplerState( ...
  model.refScene, model.satPosEci, model.satVelEci, ...
  userPosEci, userVelEci, model.wavelength, model.timeOffsetSec, model.refSatIdxLocal);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = latlon;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFdSeq = refDopplerState.deltaFd;
doaState.deltaFdRef = refDopplerState.deltaFdRef;
doaState.deltaFdRate = refDopplerState.deltaFdRate;
doaState.fdSatSeq = refDopplerState.fdSat;
doaState.fdRefSeq = refDopplerState.fdRef;
doaState.fdRefGeom = refDopplerState.fdRefRefFrame;
doaState.fdSatGeom = refDopplerState.fdSatRefFrame;
doaState.deltaFdGeom = refDopplerState.deltaFdRefFrame;
doaState.localDoaArr = localDoaArr;
end

function [fdRefPerSat, fdRatePerSat] = localFitSatFdLine(timeOffsetSec, fdSeq)
%LOCALFITSATFDLINE Fit one affine Doppler trend per satellite.

numSat = size(fdSeq, 1);
fdRefPerSat = zeros(numSat, 1);
fdRatePerSat = zeros(numSat, 1);

timeVec = reshape(timeOffsetSec, [], 1);
if numel(timeVec) < 2 || max(timeVec) - min(timeVec) <= 0
  fdRefPerSat = fdSeq(:, 1);
  return;
end

designMat = [ones(numel(timeVec), 1), timeVec];
for iSat = 1:numSat
  fdVec = reshape(fdSeq(iSat, :), [], 1);
  validMask = isfinite(fdVec) & isfinite(timeVec);
  if nnz(validMask) == 0
    continue;
  elseif nnz(validMask) == 1
    fdRefPerSat(iSat) = fdVec(validMask);
    continue;
  end
  coef = designMat(validMask, :) \ fdVec(validMask);
  fdRefPerSat(iSat) = coef(1);
  fdRatePerSat(iSat) = coef(2);
end
end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2*pi); ...
         asin(max(min(direction(3, :), 1), -1))];
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
