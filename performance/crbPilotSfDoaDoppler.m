function [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRef, pathGain, noiseVar, modelOpt)
%CRBPILOTSFDOADOPPLER Deterministic-model CRB for single-frame pilot DoA-Doppler estimation.
%
%Syntax:
%  crb = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
%    doaParam, fdRef, pathGain, noiseVar)
%
%  [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
%    doaParam, fdRef, pathGain, noiseVar, modelOpt)
%
%Inputs:
%  scene       - Multi-satellite scene structure.
%                Required fields:
%                  .array, .rotMat, .satPosEci, .satVelEci
%                Optional fields:
%                  .ref, .utc
%
%  pilotWave   - Known pilot waveform, size 1xN or Nx1.
%
%  carrierFreq - Carrier frequency in Hz.
%
%  sampleRate  - Sampling rate in Hz.
%                Used when modelOpt.timeInput is not provided.
%
%  doaParam    - 2x1 global DoA parameter vector.
%                angle mode : [eciAz; eciEl] in radians
%                latlon mode: [lat; lon] in degrees
%
%  fdRef       - Scalar reference Doppler in Hz.
%
%  pathGain    - Deterministic per-satellite complex path gains.
%                - scalar : same gain for all satellites
%                - vector : Ns x 1 or 1 x Ns
%
%  noiseVar    - Noise power.
%                - scalar : same noise power for all satellites
%                - vector : Ns x 1 or 1 x Ns
%
%  modelOpt    - Optional settings structure.
%                .doaType    : 'angle' or 'latlon', default 'angle'
%                .lightSpeed : propagation speed in m/s,
%                              default 299792458
%                .timeInput  : scalar sample rate or 1xN time axis,
%                              default sampleRate
%                .jacOpt     : options passed to buildDoaDopplerJacobian,
%                              default struct()
%
%Outputs:
%  crb         - CRB matrix of the interest parameters
%                [doaParam(1); doaParam(2); fdRef].
%
%  aux         - Auxiliary structure with fields:
%                .modelType
%                .doaType
%                .phaseMode
%                .fdRateMode
%                .paramNameModel
%                .paramNameFull
%                .paramNameInterest
%                .paramName          (alias of .paramNameInterest)
%                .activeParamIdx
%                .nuisanceParamIdx
%                .doaState
%                .jacobianModel
%                .jacobianFull
%                .jacobianInterest
%                .jacobian          (alias of .jacobianFull)
%                .fimFull
%                .fimInterest
%                .fim               (alias of .fimInterest)
%                .crbFull
%                .frameMask
%                .pathGain
%                .pathAmp
%                .pathPhase
%                .noiseVar
%                .wavelength
%
%Description:
%  Builds the deterministic-model CRB for the single-source single-frame
%  multi-satellite pilot model
%
%    y_l = g_l * vec(a_l(doa) * s_l(fd_l)) + w_l,
%
%  where g_l is the complex path gain of satellite l and fd_l is induced by
%  the global DoA and the reference Doppler fdRef. The returned CRB is for
%  the structural parameters [doaParam; fdRef]. Satellite-wise phase and
%  amplitude terms are treated as nuisance parameters and eliminated using
%  the Schur complement of the full FIM.
%
%Notes:
%  - This function is single-source only.
%  - The first two CRB dimensions inherit the units of doaParam.
%  - The third CRB dimension is always in Hz.
%
%See also:
%  crbPilotMfDoaDoppler, buildDoaDopplerJacobian,
%  estimatorDoaDopplerMlePilotSfOpt, steeringMatrix

arguments
  scene (1,1) struct
  pilotWave {mustBeNumeric, mustBeFinite}
  carrierFreq (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  sampleRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  doaParam {mustBeNumeric, mustBeFinite}
  fdRef (1,1) {mustBeNumeric, mustBeFinite}
  pathGain {mustBeNumeric, mustBeFinite}
  noiseVar {mustBeNumeric, mustBeFinite}
  modelOpt (1,1) struct = struct()
end

modelOpt = localParseModelOpt(modelOpt, sampleRate);
model = localBuildModel(scene, carrierFreq, modelOpt);
numSat = model.numSat;

pilotRow = localParsePilotWave(pilotWave);
numSample = numel(pilotRow);
if ~isscalar(modelOpt.timeInput) && numel(modelOpt.timeInput) ~= numSample
  error('crbPilotSfDoaDoppler:TimeAxisLengthMismatch', ...
    'modelOpt.timeInput length must match the pilotWave length.');
end

doaParam = reshape(doaParam, [], 1);
if numel(doaParam) ~= 2
  error('crbPilotSfDoaDoppler:InvalidDoaParamSize', ...
    'doaParam must contain exactly two entries.');
end

pathGain = localExpandPerSat(pathGain, numSat, 'pathGain');
noiseVar = localExpandPerSat(noiseVar, numSat, 'noiseVar');
if any(noiseVar <= 0)
  error('crbPilotSfDoaDoppler:InvalidNoiseVar', ...
    'noiseVar must be positive.');
end

[jacobianModel, doaState] = buildDoaDopplerJacobian( ...
  model.geom, doaParam, fdRef, modelOpt.jacOpt);
numInterest = numel(jacobianModel.paramName);

numPhase = numSat;
numAmp = numSat;
paramNamePhase = arrayfun(@(k) sprintf('phiSat%d', k), 1:numSat, 'UniformOutput', false);
paramNameAmp = arrayfun(@(k) sprintf('ampSat%d', k), 1:numSat, 'UniformOutput', false);
paramNameFull = [jacobianModel.paramName, paramNamePhase, paramNameAmp];

blockLen = model.numElement(:) * numSample;
totalDim = sum(blockLen);
numParamFull = numel(paramNameFull);
fullJac = zeros(totalDim, numParamFull);
atomCell = cell(1, numSat);
dAtomCell = cell(numSat, numInterest);

pathAmp = abs(pathGain);
pathPhase = angle(pathGain);

startIdx = 1;
for iSat = 1:numSat
  currentIdx = startIdx:(startIdx + blockLen(iSat) - 1);
  startIdx = currentIdx(end) + 1;

  currentDoa = doaState.localDoaCell{iSat};
  [steeringVec, dSteering] = steeringMatrix(model.array{iSat}, model.wavelength, currentDoa);
  [pilotFd, dPilotDfd] = buildDopplerPilotJacobian(pilotRow, doaState.fd(iSat), modelOpt.timeInput);

  atomMat = steeringVec * pilotFd;
  atomVec = atomMat(:);
  atomCell{iSat} = atomMat;

  noiseScale = 1 / sqrt(noiseVar(iSat));
  gainSat = pathGain(iSat);
  meanVec = noiseScale * gainSat * atomVec;

  for iParam = 1:numInterest
    localDoaJac = zeros(2, 1);
    if iParam <= jacobianModel.numDoaParam
      localDoaJac = jacobianModel.localDoa(:, iParam, iSat);
    end
    dFd = jacobianModel.fd(iSat, iParam);

    dSteeringMat = dSteering.az * localDoaJac(1) + dSteering.el * localDoaJac(2);
    dPilotMat = dPilotDfd * dFd;
    dAtomMat = dSteeringMat * pilotFd + steeringVec * dPilotMat;
    dAtomCell{iSat, iParam} = dAtomMat;

    fullJac(currentIdx, iParam) = noiseScale * gainSat * dAtomMat(:);
  end

  fullJac(currentIdx, numInterest + iSat) = 1j * meanVec;
  fullJac(currentIdx, numInterest + numPhase + iSat) = ...
    noiseScale * exp(1j * pathPhase(iSat)) * atomVec;
end

fimFull = 2 * real(fullJac' * fullJac);
fimFull = localSymmetrize(fimFull);
interestIdx = 1:numInterest;
nuisanceIdx = (numInterest + 1):numParamFull;
[fimInterest, crb, crbFull] = localExtractInterestCrb(fimFull, interestIdx, nuisanceIdx);

if nargout >= 2
  aux = struct();
  aux.modelType = 'sfStatic';
  aux.doaType = modelOpt.doaType;
  aux.phaseMode = 'singleFrame';
  aux.fdRateMode = 'zero';
  aux.paramNameModel = jacobianModel.paramName;
  aux.paramNameFull = paramNameFull;
  aux.paramNameInterest = jacobianModel.paramName;
  aux.paramName = aux.paramNameInterest;
  aux.activeParamIdx = interestIdx;
  aux.nuisanceParamIdx = nuisanceIdx;
  aux.doaState = doaState;
  aux.jacobianModel = jacobianModel;
  aux.jacobianFull = fullJac;
  aux.jacobianInterest = fullJac(:, interestIdx);
  aux.jacobian = aux.jacobianFull;
  aux.fimFull = fimFull;
  aux.fimInterest = fimInterest;
  aux.fim = fimInterest;
  aux.crbFull = crbFull;
  aux.frameMask = true(numSat, 1);
  aux.atomCell = atomCell;
  aux.dAtomCell = dAtomCell;
  aux.pathGain = pathGain;
  aux.pathAmp = pathAmp;
  aux.pathPhase = pathPhase;
  aux.noiseVar = noiseVar;
  aux.wavelength = model.wavelength;
else
  aux = [];
end

end


function modelOpt = localParseModelOpt(modelOpt, sampleRate)
%LOCALPARSEMODELOPT Parse optional CRB settings.

if ~isfield(modelOpt, 'doaType') || isempty(modelOpt.doaType)
  modelOpt.doaType = 'angle';
end
if ~isfield(modelOpt, 'lightSpeed') || isempty(modelOpt.lightSpeed)
  modelOpt.lightSpeed = 299792458;
end
if ~isfield(modelOpt, 'timeInput') || isempty(modelOpt.timeInput)
  modelOpt.timeInput = sampleRate;
end
if ~isfield(modelOpt, 'jacOpt') || isempty(modelOpt.jacOpt)
  modelOpt.jacOpt = struct();
end

modelOpt.doaType = validatestring(modelOpt.doaType, {'angle', 'latlon'}, ...
  mfilename, 'modelOpt.doaType');
if ~isscalar(modelOpt.lightSpeed) || ~isfinite(modelOpt.lightSpeed) || modelOpt.lightSpeed <= 0
  error('crbPilotSfDoaDoppler:InvalidLightSpeed', ...
    'modelOpt.lightSpeed must be a positive finite scalar.');
end
if ~(isscalar(modelOpt.timeInput) || isvector(modelOpt.timeInput))
  error('crbPilotSfDoaDoppler:InvalidTimeInput', ...
    'modelOpt.timeInput must be a scalar sample rate or a sample-time vector.');
end
if ~isstruct(modelOpt.jacOpt)
  error('crbPilotSfDoaDoppler:InvalidJacOpt', ...
    'modelOpt.jacOpt must be a struct.');
end

end


function model = localBuildModel(scene, carrierFreq, modelOpt)
%LOCALBUILDMODEL Build fixed scene quantities used by the CRB.

arrayCell = localParseSceneArray(scene);
rotMatCell = localParseRotMat(scene);
[satPosEci, satVelEci] = localParseSatState(scene);
refState = localParseReferenceState(scene, satPosEci, satVelEci);

numSat = numel(arrayCell);
if numel(rotMatCell) ~= numSat
  error('crbPilotSfDoaDoppler:RotMatCountMismatch', ...
    'The number of rotation matrices must match the number of arrays.');
end
if size(satPosEci, 2) ~= numSat || size(satVelEci, 2) ~= numSat
  error('crbPilotSfDoaDoppler:SatCountMismatch', ...
    'scene.satPosEci and scene.satVelEci must match the number of arrays.');
end

numElement = zeros(numSat, 1);
for iSat = 1:numSat
  currentArray = arrayCell{iSat};
  if ~isstruct(currentArray) || ~isfield(currentArray, 'positions') || isempty(currentArray.positions)
    error('crbPilotSfDoaDoppler:InvalidArrayEntry', ...
      'Each scene.array entry must contain a non-empty positions field.');
  end
  numElement(iSat) = size(currentArray.positions, 2);
end

wavelength = modelOpt.lightSpeed / carrierFreq;
geom = struct();
geom.doaType = modelOpt.doaType;
geom.rotMat = rotMatCell;
geom.wavelength = wavelength;
geom.ref = refState;
if isfield(scene, 'utc') && ~isempty(scene.utc)
  geom.sceneUtc = scene.utc;
end

switch modelOpt.doaType
  case 'angle'
    geom.relSatVelEci = satVelEci - refState.velEci;

  case 'latlon'
    geom.satPosEci = satPosEci;
    geom.satVelEci = satVelEci;

  otherwise
    error('crbPilotSfDoaDoppler:InvalidDoaType', ...
      'Unsupported doaType: %s.', modelOpt.doaType);
end

model = struct();
model.array = arrayCell;
model.numSat = numSat;
model.numElement = numElement;
model.wavelength = wavelength;
model.geom = geom;

end


function pilotRow = localParsePilotWave(pilotWave)
%LOCALPARSEPILOTWAVE Normalize pilot waveform to a row vector.

if ~isvector(pilotWave) || isempty(pilotWave)
  error('crbPilotSfDoaDoppler:InvalidPilotWave', ...
    'pilotWave must be a non-empty vector.');
end
pilotRow = reshape(pilotWave, 1, []);
if real(pilotRow * pilotRow') <= eps
  error('crbPilotSfDoaDoppler:ZeroPilotEnergy', ...
    'pilotWave must have non-zero energy.');
end

end


function valueVec = localExpandPerSat(valueIn, numSat, nameText)
%LOCALEXPANDPERSAT Expand scalar or vector input to an Ns x 1 vector.

if isscalar(valueIn)
  valueVec = repmat(valueIn, numSat, 1);
  return;
end

if ~isvector(valueIn) || numel(valueIn) ~= numSat
  error('crbPilotSfDoaDoppler:PerSatSizeMismatch', ...
    '%s must be a scalar or a vector with one entry per satellite.', nameText);
end

valueVec = reshape(valueIn, [], 1);

end


function arrayCell = localParseSceneArray(scene)
%LOCALPARSESCENEARRAY Convert scene.array to a 1xNs cell array.

if ~isfield(scene, 'array') || isempty(scene.array)
  error('crbPilotSfDoaDoppler:MissingArrayField', ...
    'scene.array is required.');
end

if iscell(scene.array)
  arrayCell = reshape(scene.array, 1, []);
elseif isstruct(scene.array)
  if isscalar(scene.array)
    arrayCell = {scene.array};
  else
    arrayCell = num2cell(reshape(scene.array, 1, []));
  end
else
  error('crbPilotSfDoaDoppler:InvalidArrayField', ...
    'scene.array must be a struct, struct array, or cell array.');
end

end


function rotMatCell = localParseRotMat(scene)
%LOCALPARSEROTMAT Convert scene.rotMat to a 1xNs cell array.

if ~isfield(scene, 'rotMat') || isempty(scene.rotMat)
  error('crbPilotSfDoaDoppler:MissingRotMatField', ...
    'scene.rotMat is required.');
end

if isnumeric(scene.rotMat)
  if ~isequal(size(scene.rotMat), [3, 3])
    error('crbPilotSfDoaDoppler:InvalidRotMatSize', ...
      'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  rotMatCell = {scene.rotMat};
  return;
end

if ~iscell(scene.rotMat) || isempty(scene.rotMat)
  error('crbPilotSfDoaDoppler:InvalidRotMatType', ...
    'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

rotMatCell = reshape(scene.rotMat, 1, []);
for iSat = 1:numel(rotMatCell)
  if ~isnumeric(rotMatCell{iSat}) || ~isequal(size(rotMatCell{iSat}), [3, 3])
    error('crbPilotSfDoaDoppler:InvalidRotMatCell', ...
      'Each scene.rotMat entry must be a 3x3 numeric matrix.');
  end
end

end


function [satPosEci, satVelEci] = localParseSatState(scene)
%LOCALPARSESATSTATE Validate satellite ECI state arrays.

if ~isfield(scene, 'satPosEci') || ~isfield(scene, 'satVelEci') || ...
    isempty(scene.satPosEci) || isempty(scene.satVelEci)
  error('crbPilotSfDoaDoppler:MissingSatState', ...
    'scene.satPosEci and scene.satVelEci are required.');
end

satPosEci = scene.satPosEci;
satVelEci = scene.satVelEci;
if ~isnumeric(satPosEci) || ~isnumeric(satVelEci) || ...
    ~isequal(size(satPosEci), size(satVelEci)) || size(satPosEci, 1) ~= 3
  error('crbPilotSfDoaDoppler:InvalidSatStateSize', ...
    'scene.satPosEci and scene.satVelEci must both have size 3xNs.');
end

end


function ref = localParseReferenceState(scene, satPosEci, satVelEci)
%LOCALPARSEREFERENCESTATE Parse the reference state used in Doppler mapping.

numSat = size(satPosEci, 2);

if ~isfield(scene, 'ref') || isempty(scene.ref)
  refWeight = ones(numSat, 1) / numSat;
  ref = struct();
  ref.type = 'centroid';
  ref.weight = refWeight;
  ref.posEci = satPosEci * refWeight;
  ref.velEci = satVelEci * refWeight;
  return;
end

refIn = scene.ref;
if ~isstruct(refIn) || ~isfield(refIn, 'posEci') || ~isfield(refIn, 'velEci')
  error('crbPilotSfDoaDoppler:MissingReferenceState', ...
    'scene.ref.posEci and scene.ref.velEci are required when scene.ref is provided.');
end
if ~isnumeric(refIn.posEci) || ~isequal(size(refIn.posEci), [3, 1])
  error('crbPilotSfDoaDoppler:InvalidReferencePosSize', ...
    'scene.ref.posEci must have size 3x1.');
end
if ~isnumeric(refIn.velEci) || ~isequal(size(refIn.velEci), [3, 1])
  error('crbPilotSfDoaDoppler:InvalidReferenceVelSize', ...
    'scene.ref.velEci must have size 3x1.');
end

ref = struct();
ref.type = 'custom';
if isfield(refIn, 'type') && ~isempty(refIn.type)
  ref.type = char(refIn.type);
end
ref.weight = [];
if isfield(refIn, 'weight') && ~isempty(refIn.weight)
  ref.weight = refIn.weight;
end
ref.posEci = refIn.posEci;
ref.velEci = refIn.velEci;

end


function [fimInterest, crbInterest, crbFull] = localExtractInterestCrb(fimFull, interestIdx, nuisanceIdx)
%LOCALEXTRACTINTERESTCRB Extract interest-parameter CRB from the full FIM.

fimFull = localSymmetrize(fimFull);
if isempty(nuisanceIdx)
  fimInterest = fimFull(interestIdx, interestIdx);
else
  fimTheta = fimFull(interestIdx, interestIdx);
  fimCross = fimFull(interestIdx, nuisanceIdx);
  fimNuis = fimFull(nuisanceIdx, nuisanceIdx);
  fimInterest = fimTheta - fimCross * localSafeSolve(fimNuis, fimCross.');
  fimInterest = localSymmetrize(fimInterest);
end

crbInterest = localSafeInverse(fimInterest, ...
  'crbPilotSfDoaDoppler:IllConditionedInterestFim', ...
  'The interest-parameter FIM is singular or ill-conditioned. Using pinv.');
crbFull = localSafeInverse(fimFull, ...
  'crbPilotSfDoaDoppler:IllConditionedFullFim', ...
  'The full FIM is singular or ill-conditioned. Using pinv.');

end


function xMat = localSafeSolve(aMat, bMat)
%LOCALSAFESOLVE Solve a linear system with a pinv fallback.

aMat = localSymmetrize(aMat);
if rcond(aMat + eye(size(aMat)) * eps) < 1e-12
  warning('crbPilotSfDoaDoppler:IllConditionedSolve', ...
    'Nuisance FIM is singular or ill-conditioned. Using pinv in Schur complement.');
  xMat = pinv(aMat) * bMat;
else
  xMat = aMat \ bMat;
end

end


function invMat = localSafeInverse(aMat, warnId, warnMsg)
%LOCALSAFEINVERSE Invert a symmetric matrix with a pinv fallback.

aMat = localSymmetrize(aMat);
if rcond(aMat + eye(size(aMat)) * eps) < 1e-12
  warning(warnId, '%s', warnMsg);
  invMat = pinv(aMat);
else
  invMat = aMat \ eye(size(aMat));
end
invMat = localSymmetrize(invMat);

end


function aMat = localSymmetrize(aMat)
%LOCALSYMMETRIZE Force exact Hermitian symmetry on a real matrix.

aMat = 0.5 * (aMat + aMat.');

end
