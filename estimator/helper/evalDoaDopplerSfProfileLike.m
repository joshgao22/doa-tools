function [obj, pathGain, noiseVar, aux] = evalDoaDopplerSfProfileLike(model, optVar)
%EVALDOADOPPLERSFPROFILELIKE Evaluate the single-frame profiled likelihood.
% This is the extracted single-frame profile-likelihood core used by
% estimatorDoaDopplerMlePilotSfOpt.

arguments
  model (1,1) struct
  optVar (:,1) double
end

[obj, pathGain, noiseVar, aux] = localNllPilotStoCon(model, optVar);
end

function [obj, pathGain, noiseVar, aux] = localNllPilotStoCon(model, optVar)
%LOCALNLLPILOTSTOCON Concentrated negative log-likelihood.

hyp = localBuildHypothesis(model, optVar);
[pathGain, noiseVar, residualNorm, gainDiag] = localEstimateGainAndNoise(model, hyp);

countPerSat = model.numElement(:) * model.numSample;
if model.useLogObjective
  objectiveSat = model.satWeight(:) .* countPerSat .* log(max(noiseVar, eps));
  obj = sum(objectiveSat);
else
  objectiveSat = model.satWeight(:) .* gainDiag.residualEnergySat(:);
  obj = residualNorm;
end

if nargout <= 1
  return;
end

aux = struct();
aux.doaParam = hyp.doaParam;
aux.eciAngle = hyp.eciAngle;
aux.eciDirection = hyp.eciDirection;
aux.latlon = hyp.latlon;
aux.userPosEci = hyp.userPosEci;
aux.userVelEci = hyp.userVelEci;
aux.fdRef = hyp.fdRef;
aux.deltaFd = hyp.deltaFd;
aux.fd = hyp.fd;
aux.localDoaArr = hyp.localDoaArr;
aux.residualNorm = residualNorm;
aux.residualNormSat = gainDiag.residualEnergySat(:);
aux.countPerSat = countPerSat;
aux.satWeight = model.satWeight(:);
aux.objectiveSat = objectiveSat(:);
aux.noiseVarSat = noiseVar(:);
aux.noiseVarGlobal = residualNorm / max(sum(aux.countPerSat), eps);
end

function hyp = localBuildHypothesis(model, optVar)
%LOCALBUILDHYPOTHESIS Map optimization variables to per-satellite parameters.

optVar = optVar(:);
numDoaVar = 2 * model.numSource;
numOptVar = 3 * model.numSource;

if numel(optVar) ~= numOptVar
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidOptVarSize', ...
    'The optimization vector must contain %d entries.', numOptVar);
end

doaParam = reshape(optVar(1:numDoaVar), 2, model.numSource);
fdRef = optVar(numDoaVar + 1:end);
doaState = localBuildDoaState(model, doaParam);
fd = doaState.deltaFd + reshape(fdRef, 1, []);

hyp = doaState;
hyp.doaParam = doaParam;
hyp.fdRef = fdRef(:);
hyp.fd = fd;
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
