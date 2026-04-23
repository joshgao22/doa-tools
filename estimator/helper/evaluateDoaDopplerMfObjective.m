function [obj, pathGain, noiseVar, aux] = evaluateDoaDopplerMfObjective(model, optVar)
%EVALUATEDOADOPPLERMFOBJECTIVE Concentrated MF dynamic pilot objective.
% Keep hypothesis construction and profile-likelihood evaluation in one
% functional module so the branch solver only handles branch orchestration.

arguments
  model (1,1) struct
  optVar (:,1) double
end

hyp = localBuildHypothesis(model, optVar);
[obj, prof, profAux] = evalDoaDopplerDynProfileLike(model, hyp);

if nargout <= 1
  return;
end

pathGain = prof.pathGainEst;
noiseVar = prof.noiseVarEst;

aux = struct();
aux.doaParam = hyp.doaParam;
aux.eciAngle = hyp.eciAngle;
aux.eciDirection = hyp.eciDirection;
aux.latlon = hyp.latlon;
aux.userPosEci = hyp.userPosEci;
aux.userVelEci = hyp.userVelEci;
aux.refSatIdxLocal = localGetStructField(profAux, 'refSatIdxLocal', NaN);
aux.deltaFdSeq = hyp.deltaFdSeq;
aux.deltaFdRef = hyp.deltaFdRef;
aux.deltaFdRefRaw = localGetStructField(profAux, 'deltaFdRefRaw', hyp.deltaFdRef);
aux.deltaFdRefEval = localGetStructField(profAux, 'deltaFdRefEval', hyp.deltaFdRef);
aux.deltaFdRate = hyp.deltaFdRate;
aux.deltaFdRateRaw = localGetStructField(profAux, 'deltaFdRateRaw', hyp.deltaFdRate);
aux.deltaFdRateEval = localGetStructField(profAux, 'deltaFdRateEval', hyp.deltaFdRate);
aux.fdRef = hyp.fdRef;
aux.fdRate = hyp.fdRate;
aux.fdSat = hyp.fdSat;
aux.fdSatRaw = localGetStructField(profAux, 'fdSatRaw', hyp.fdSat);
aux.fdRateSat = hyp.fdRateSat;
aux.fdLocal = profAux.fdLocal;
aux.fdLocalRaw = localGetStructField(profAux, 'fdLocalRaw', profAux.fdLocal);
aux.fdSatEval = profAux.fdSatEval;
aux.fdLocalEval = profAux.fdLocalEval;
aux.fdAliasStepHz = profAux.fdAliasStepHz;
aux.fdAliasIndex = profAux.fdAliasIndex;
aux.fdAliasShiftHz = profAux.fdAliasShiftHz;
aux.localDoaArrUsed = profAux.localDoaArrUsed;
aux.localDoaArr = profAux.localDoaArrUsed;
aux.phaseSat = prof.phaseSatEst;
aux.framePhase = prof.framePhaseEst;
aux.ampEst = prof.ampEst;
aux.noiseVarGlobal = prof.noiseVarGlobal;
aux.residualNorm = prof.residualNorm;
aux.countPerSat = prof.countPerSat;
aux.residualSat = localGetStructField(prof, 'residualSat', []);
aux.fitValueSat = localGetStructField(prof, 'fitValueSat', []);
aux.objectiveSat = localGetStructField(prof, 'objectiveSat', []);
aux.effectiveFrameSupportSat = localGetStructField(prof, 'effectiveFrameSupportSat', []);
aux.effectiveFrameSupportRatioSat = localGetStructField(prof, 'effectiveFrameSupportRatioSat', []);
aux.negativeProjectionRatioSat = localGetStructField(prof, 'negativeProjectionRatioSat', []);
aux.collapsedFrameCountSat = localGetStructField(prof, 'collapsedFrameCountSat', []);
aux.satFitRatio = localGetStructField(prof, 'satFitRatio', []);
aux.refFitRatio = localGetStructField(prof, 'refFitRatio', NaN);
aux.nonRefFitRatioFloor = localGetStructField(prof, 'nonRefFitRatioFloor', NaN);
aux.refSupportRatio = localGetStructField(prof, 'refSupportRatio', NaN);
aux.nonRefSupportRatioFloor = localGetStructField(prof, 'nonRefSupportRatioFloor', NaN);
aux.refConsistencyNorm = localGetStructField(prof, 'refConsistencyNorm', NaN);
aux.nonRefConsistencyRatioFloor = localGetStructField(prof, 'nonRefConsistencyRatioFloor', NaN);
aux.maxNonRefNegativeProjectionRatio = localGetStructField(prof, 'maxNonRefNegativeProjectionRatio', NaN);
aux.nonRefFitFloorPenalty = localGetStructField(prof, 'nonRefFitFloorPenalty', 0);
aux.nonRefSupportFloorPenalty = localGetStructField(prof, 'nonRefSupportFloorPenalty', 0);
aux.additionalObjectivePenalty = localGetStructField(prof, 'additionalObjectivePenalty', 0);
aux.blockValue = prof.blockValue;
aux.blockNorm2 = prof.blockNorm2;
aux.countPerBlock = localGetStructField(prof, 'countPerBlock', []);
aux.zMat = prof.zMat;
aux.etaMat = prof.etaMat;
end


function hyp = localBuildHypothesis(model, optVar)
%LOCALBUILDHYPOTHESIS Map one optimization point to geometry and Doppler.

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
%LOCALBUILDANGLESTATE Build frame-wise states for angle-mode parameterization.

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
[userPosEci, userVelEci] = localLatlonToUserState(latlon, model.userStateRef);

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


function [userPosEci, userVelEci] = localLatlonToUserState(latlon, userStateRef)
%LOCALLATLONTOUSERSTATE Convert one ground point into motion-aware ECI states.

[userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, userStateRef);
end


function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2*pi); ...
         asin(max(min(direction(3, :), 1), -1))];
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
