function [crb, aux] = crbPilotDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRef, pathGain, noiseVar, modelOpt)
%CRBPILOTDOADOPPLER Deterministic-model CRB for joint DoA-Doppler parameters.
%
%Syntax:
%  crb = crbPilotDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
%    doaParam, fdRef, pathGain, noiseVar)
%
%  [crb, aux] = crbPilotDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
%    doaParam, fdRef, pathGain, noiseVar, modelOpt)
%
%Inputs:
%  scene       - Multi-satellite scene structure
%                required fields:
%                  .array, .rotMat, .satPosEci, .satVelEci
%                optional fields:
%                  .ref, .utc
%
%  pilotWave   - Known pilot waveform, size 1xN or Nx1
%                This waveform is used directly as the deterministic source
%                template. If zero padding is needed, it should already be
%                included in pilotWave.
%
%  carrierFreq - Carrier frequency in Hz
%
%  sampleRate  - Sampling rate in Hz
%                Ignored when modelOpt.timeInput is provided.
%
%  doaParam    - 2x1 global DoA parameter vector
%                angle mode : [eciAz; eciEl] in radians
%                latlon mode: [lat; lon] in degrees
%
%  fdRef       - Scalar reference Doppler in Hz
%
%  pathGain    - Deterministic per-satellite complex path gains
%                - scalar : same gain for all satellites
%                - vector : Ns x 1 or 1 x Ns
%
%  noiseVar    - Noise power
%                - scalar : same noise power for all satellites
%                - vector : Ns x 1 or 1 x Ns
%
%  modelOpt    - Optional settings structure
%                .doaType    : 'angle' or 'latlon', default 'angle'
%                .lightSpeed : propagation speed in m/s,
%                              default 299792458
%                .timeInput  : scalar sample rate or 1xN time axis,
%                              default sampleRate
%                .jacOpt     : options passed to buildDoaDopplerJacobian,
%                              default struct()
%
%Outputs:
%  crb         - CRB matrix of the structural parameters
%                angle mode : 3x3 for [eciAz; eciEl; fdRef]
%                latlon mode: 3x3 for [lat; lon; fdRef]
%
%  aux         - Auxiliary structure with fields:
%                .paramNameModel
%                .paramNameFull
%                .paramNameInterest
%                .paramName          (alias of .paramNameInterest)
%                .fdRateMode
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
%                .fimCore
%                .atomCell
%                .dAtomCell
%                .pathGain
%                .noiseVar
%                .wavelength
%
%Description:
%  Builds the deterministic-model CRB for the single-source multi-satellite
%  DoA-Doppler model
%
%    y_k = alpha_k * vec(a_k(doa) * s_k(fd_k)) + w_k,
%
%  where each satellite has its own complex path gain alpha_k, local DoA,
%  and Doppler fd_k = fdRef + deltaFd_k. The nuisance path gains are
%  eliminated using the orthogonal projector onto the complement of the
%  whitened atom subspace.
%
%Notes:
%  - This function is single-source only.
%  - The first two CRB dimensions inherit the parameter units of doaParam:
%      * angle mode  -> radians
%      * latlon mode -> degrees
%    The third dimension is always in Hz.
%  - If noiseVar is unknown but satellite-wise white, treating it as an
%    additional nuisance parameter does not change this structural CRB.
%
%See also:
%  crbDetDoa, buildDoaDopplerJacobian, buildDopplerPilotJacobian,
%  estimatorDoaDopplerMlePilotOpt, steeringMatrix

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
  error('crbDetDoaDoppler:TimeAxisLengthMismatch', ...
    'modelOpt.timeInput length must match the pilotWave length.');
end

doaParam = reshape(doaParam, [], 1);
if numel(doaParam) ~= 2
  error('crbDetDoaDoppler:InvalidDoaParamSize', ...
    'doaParam must contain exactly two entries.');
end

pathGain = localExpandPerSat(pathGain, numSat, 'pathGain');
noiseVar = localExpandPerSat(noiseVar, numSat, 'noiseVar');
if any(noiseVar <= 0)
  error('crbDetDoaDoppler:InvalidNoiseVar', ...
    'noiseVar must be positive.');
end

[jacobian, doaState] = buildDoaDopplerJacobian(model.geom, doaParam, fdRef, modelOpt.jacOpt);
numParam = numel(jacobian.paramName);

blockLen = model.numElement(:) * numSample;
totalDim = sum(blockLen);

atomMat = zeros(totalDim, numSat);
derivMat = zeros(totalDim, numParam);
atomCell = cell(1, numSat);
dAtomCell = cell(numSat, numParam);

startIdx = 1;
for iSat = 1:numSat
  currentBlock = startIdx:(startIdx + blockLen(iSat) - 1);
  startIdx = currentBlock(end) + 1;

  currentDoa = doaState.localDoaCell{iSat};
  [steeringVec, dSteering] = steeringMatrix(model.array{iSat}, model.wavelength, currentDoa);
  [pilotFd, dPilotDfd] = buildDopplerPilotJacobian(pilotRow, doaState.fd(iSat), modelOpt.timeInput);

  atomCurrent = steeringVec * pilotFd;
  atomVec = atomCurrent(:);
  noiseScale = 1 / sqrt(noiseVar(iSat));

  atomCell{iSat} = atomCurrent;
  atomMat(currentBlock, iSat) = noiseScale * atomVec;

  for iParam = 1:numParam
    localDoaJac = jacobian.localDoa(:, min(iParam, 2), iSat);
    if iParam > 2
      localDoaJac = [0; 0];
    end
    dFd = jacobian.fd(iSat, iParam);

    dSteeringCurrent = dSteering.az * localDoaJac(1) + dSteering.el * localDoaJac(2);
    dPilotCurrent = dPilotDfd * dFd;
    dAtomCurrent = dSteeringCurrent * pilotFd + steeringVec * dPilotCurrent;

    dAtomCell{iSat, iParam} = dAtomCurrent;
    derivMat(currentBlock, iParam) = derivMat(currentBlock, iParam) + ...
      noiseScale * pathGain(iSat) * dAtomCurrent(:);
  end
end

gramAtom = atomMat' * atomMat;
gramAtom = 0.5 * (gramAtom + gramAtom');
if rcond(gramAtom) < 1e-12
  warning('crbDetDoaDoppler:IllConditionedAtom', ...
    'The whitened atom Gram matrix is singular or ill-conditioned. CRB may be unreliable.');
end

crossTerm = atomMat' * derivMat;
projGram = derivMat' * derivMat - crossTerm' * (gramAtom \ crossTerm);
projGram = 0.5 * (projGram + projGram');

fim = 2 * real(projGram);
fim = 0.5 * (fim + fim.');
if rcond(fim) < 1e-12
  warning('crbDetDoaDoppler:IllConditionedFim', ...
    'FIM is singular or ill-conditioned. CRB may be unreliable.');
end

crb = fim \ eye(numParam);
crb = 0.5 * (crb + crb.');

if nargout >= 2
  aux = struct();
  aux.paramNameModel = jacobian.paramName;
  aux.paramNameFull = jacobian.paramName;
  aux.paramNameInterest = jacobian.paramName;
  aux.paramName = jacobian.paramName;
  aux.fdRateMode = 'fixedzero';
  aux.activeParamIdx = 1:numParam;
  aux.nuisanceParamIdx = [];
  aux.doaState = doaState;
  aux.jacobianModel = jacobian;
  aux.jacobianFull = jacobian;
  aux.jacobianInterest = jacobian;
  aux.jacobian = jacobian;
  aux.fimFull = fim;
  aux.fimInterest = fim;
  aux.fim = fim;
  aux.crbFull = crb;
  aux.fimCore = projGram;
  aux.atomCell = atomCell;
  aux.dAtomCell = dAtomCell;
  aux.pathGain = pathGain;
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
  error('crbDetDoaDoppler:InvalidLightSpeed', ...
    'modelOpt.lightSpeed must be a positive finite scalar.');
end
if ~(isscalar(modelOpt.timeInput) || isvector(modelOpt.timeInput))
  error('crbDetDoaDoppler:InvalidTimeInput', ...
    'modelOpt.timeInput must be a scalar sample rate or a sample-time vector.');
end
if ~isstruct(modelOpt.jacOpt)
  error('crbDetDoaDoppler:InvalidJacOpt', ...
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
  error('crbDetDoaDoppler:RotMatCountMismatch', ...
    'The number of rotation matrices must match the number of arrays.');
end
if size(satPosEci, 2) ~= numSat || size(satVelEci, 2) ~= numSat
  error('crbDetDoaDoppler:SatCountMismatch', ...
    'scene.satPosEci and scene.satVelEci must match the number of arrays.');
end

numElement = zeros(numSat, 1);
for iSat = 1:numSat
  currentArray = arrayCell{iSat};
  if ~isstruct(currentArray) || ~isfield(currentArray, 'positions') || isempty(currentArray.positions)
    error('crbDetDoaDoppler:InvalidArrayEntry', ...
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
geom.scene = scene;
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
    error('crbDetDoaDoppler:InvalidDoaType', ...
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
  error('crbDetDoaDoppler:InvalidPilotWave', ...
    'pilotWave must be a non-empty vector.');
end
pilotRow = reshape(pilotWave, 1, []);

pilotEnergy = real(pilotRow * pilotRow');
if pilotEnergy <= eps
  error('crbDetDoaDoppler:ZeroPilotEnergy', ...
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
  error('crbDetDoaDoppler:PerSatSizeMismatch', ...
    '%s must be a scalar or a vector with one entry per satellite.', nameText);
end

valueVec = reshape(valueIn, [], 1);

end

function arrayCell = localParseSceneArray(scene)
%LOCALPARSESCENEARRAY Convert scene.array to a 1xNs cell array.

if ~isfield(scene, 'array') || isempty(scene.array)
  error('crbDetDoaDoppler:MissingArrayField', ...
    'scene.array is required.');
end

arrayInput = scene.array;
if iscell(arrayInput)
  arrayCell = reshape(arrayInput, 1, []);
elseif isstruct(arrayInput)
  if isscalar(arrayInput)
    arrayCell = {arrayInput};
  else
    arrayCell = num2cell(reshape(arrayInput, 1, []));
  end
else
  error('crbDetDoaDoppler:InvalidArrayField', ...
    'scene.array must be a struct, struct array, or cell array.');
end

end

function rotMatCell = localParseRotMat(scene)
%LOCALPARSEROTMAT Convert scene.rotMat to a 1xNs cell array.

if ~isfield(scene, 'rotMat') || isempty(scene.rotMat)
  error('crbDetDoaDoppler:MissingRotMatField', ...
    'scene.rotMat is required.');
end

rotMat = scene.rotMat;
if isnumeric(rotMat)
  if ~isequal(size(rotMat), [3, 3])
    error('crbDetDoaDoppler:InvalidRotMatSize', ...
      'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  rotMatCell = {rotMat};
  return;
end

if ~iscell(rotMat) || isempty(rotMat)
  error('crbDetDoaDoppler:InvalidRotMatType', ...
    'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

rotMatCell = reshape(rotMat, 1, []);
for iSat = 1:numel(rotMatCell)
  if ~isnumeric(rotMatCell{iSat}) || ~isequal(size(rotMatCell{iSat}), [3, 3])
    error('crbDetDoaDoppler:InvalidRotMatCell', ...
      'Each scene.rotMat entry must be a 3x3 numeric matrix.');
  end
end

end

function [satPosEci, satVelEci] = localParseSatState(scene)
%LOCALPARSESATSTATE Validate satellite ECI state arrays.

if ~isfield(scene, 'satPosEci') || ~isfield(scene, 'satVelEci') || ...
    isempty(scene.satPosEci) || isempty(scene.satVelEci)
  error('crbDetDoaDoppler:MissingSatState', ...
    'scene.satPosEci and scene.satVelEci are required.');
end

satPosEci = scene.satPosEci;
satVelEci = scene.satVelEci;
if ~isnumeric(satPosEci) || ~isnumeric(satVelEci) || ...
    ~isequal(size(satPosEci), size(satVelEci)) || size(satPosEci, 1) ~= 3
  error('crbDetDoaDoppler:InvalidSatStateSize', ...
    'scene.satPosEci and scene.satVelEci must both have size 3xNs.');
end

end

function ref = localParseReferenceState(scene, satPosEci, satVelEci)
%LOCALPARSEREFERENCESTATE Parse the Doppler reference state.

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
  error('crbDetDoaDoppler:MissingReferenceState', ...
    'scene.ref.posEci and scene.ref.velEci are required when scene.ref is provided.');
end
if ~isnumeric(refIn.posEci) || ~isequal(size(refIn.posEci), [3, 1])
  error('crbDetDoaDoppler:InvalidReferencePosSize', ...
    'scene.ref.posEci must have size 3x1.');
end
if ~isnumeric(refIn.velEci) || ~isequal(size(refIn.velEci), [3, 1])
  error('crbDetDoaDoppler:InvalidReferenceVelSize', ...
    'scene.ref.velEci must have size 3x1.');
end

ref = struct();
if isfield(refIn, 'type') && ~isempty(refIn.type)
  ref.type = char(refIn.type);
else
  ref.type = 'custom';
end
if isfield(refIn, 'weight') && ~isempty(refIn.weight)
  ref.weight = refIn.weight;
else
  ref.weight = [];
end
ref.posEci = refIn.posEci;
ref.velEci = refIn.velEci;

end
