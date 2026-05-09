function [crb, aux] = crbPilotSfDoaOnlyEffective(scene, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRef, pathGain, noiseVar, modelOpt)
%CRBPILOTSFDOAONLYEFFECTIVE DoA-only deterministic CRB for pilot SF data.
%
%Syntax:
%  crb = crbPilotSfDoaOnlyEffective(scene, pilotWave, carrierFreq, sampleRate, ...
%    doaParam, fdRef, pathGain, noiseVar)
%  [crb, aux] = crbPilotSfDoaOnlyEffective(..., modelOpt)
%
%Description:
%  Builds a DoA-only deterministic CRB using the array-only CRB kernel with
%  a per-satellite source power matched to the pilot model. This helper is
%  intended for DoA-only estimators that ignore Doppler in the atom while
%  the received known pilot is Doppler shifted. In the default
%  "pilotModel" mode, each satellite contribution is scaled by the
%  coherent projection of the Doppler-shifted pilot onto the no-Doppler
%  pilot atom. In "unit" mode this projection gain is fixed to one.
%
%Inputs:
%  scene       - Single-frame satellite scene structure.
%  pilotWave   - Known pilot waveform vector.
%  carrierFreq - Carrier frequency in Hz.
%  sampleRate  - Sampling rate in Hz.
%  doaParam    - 2x1 DoA parameter vector.
%                angle mode : [eciAz; eciEl] in radians
%                latlon mode: [lat; lon] in degrees
%  fdRef       - Reference Doppler in Hz used to induce per-sat Doppler.
%  pathGain    - Scalar or per-sat complex deterministic path gain.
%  noiseVar    - Scalar or per-sat noise variance.
%  modelOpt    - Optional settings structure.
%                .doaType           : 'angle' or 'latlon', default 'angle'
%                .effectiveGainMode : 'pilotModel' or 'unit', default 'pilotModel'
%                .lightSpeed        : propagation speed, default 299792458
%                .timeInput         : sample rate or time vector, default sampleRate
%                .jacOpt            : options passed to buildDoaDopplerJacobian
%
%Outputs:
%  crb         - 2x2 CRB of the DoA parameters.
%  aux         - Auxiliary structure with FIM, effective gains and source
%                powers passed to crbDetDoa.
%
%Notes:
%  - This is not a replacement for crbPilotSfDoaDoppler. Matched
%    DoA-Doppler estimators should still use crbPilotSfDoaDoppler.
%  - crbDetDoa remains the generic array-only CRB. This wrapper only maps
%    the pilot / Doppler / path-gain model to an effective deterministic
%    source power for that CRB.
%
%See also:
%  crbDetDoa, crbPilotSfDoaDoppler, buildDopplerPilot

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
pathGain = localExpandPerSat(pathGain, numSat, 'pathGain');
noiseVar = localExpandPerSat(noiseVar, numSat, 'noiseVar');
if any(noiseVar <= 0)
  error('crbPilotSfDoaOnlyEffective:InvalidNoiseVar', ...
    'noiseVar must be positive.');
end

doaParam = reshape(doaParam, [], 1);
if numel(doaParam) ~= 2
  error('crbPilotSfDoaOnlyEffective:InvalidDoaParamSize', ...
    'doaParam must contain exactly two entries.');
end

[jacobianModel, doaState] = buildDoaDopplerJacobian( ...
  model.geom, doaParam, fdRef, modelOpt.jacOpt);

localDoaCell = reshape(doaState.localDoaCell, 1, []);
localDoaJacCell = cell(1, numSat);
for iSat = 1:numSat
  localDoaJacCell{iSat} = jacobianModel.localDoa(:, :, iSat);
end

effectivePilotGainAbs = localBuildEffectivePilotGainAbs(pilotRow, doaState.fd(:), ...
  modelOpt.timeInput, modelOpt.effectiveGainMode);
pilotMeanPower = mean(abs(pilotRow(:)).^2);
pwrSourceCell = cell(1, numSat);
for iSat = 1:numSat
  pwrSourceCell{iSat} = pilotMeanPower * abs(pathGain(iSat)).^2 * ...
    effectivePilotGainAbs(iSat)^2;
end

crbDetOpt = struct();
crbDetOpt.localDoaJac = localDoaJacCell;
crbDetOpt.paramName = jacobianModel.paramNameDoa;
[crb, detAux] = crbDetDoa(model.array, model.wavelength, localDoaCell, ...
  pwrSourceCell, noiseVar, numel(pilotRow), crbDetOpt);

if nargout >= 2
  aux = detAux;
  aux.modelType = 'sfDoaOnlyEffective';
  aux.doaType = modelOpt.doaType;
  aux.effectiveGainMode = modelOpt.effectiveGainMode;
  aux.doaState = doaState;
  aux.jacobianModel = jacobianModel;
  aux.pathGain = pathGain;
  aux.pathAmp = abs(pathGain);
  aux.noiseVar = noiseVar;
  aux.pilotMeanPower = pilotMeanPower;
  aux.effectivePilotGainAbs = effectivePilotGainAbs;
  aux.pwrSourceCell = pwrSourceCell;
  aux.wavelength = model.wavelength;
else
  aux = [];
end

end

function modelOpt = localParseModelOpt(modelOpt, sampleRate)
%LOCALPARSEMODELOPT Parse optional CRB settings.

if ~isfield(modelOpt, 'doaType') || isempty(modelOpt.doaType)
  modelOpt.doaType = 'angle';
else
  modelOpt.doaType = char(lower(string(modelOpt.doaType)));
end

if ~ismember(modelOpt.doaType, {'angle', 'latlon'})
  error('crbPilotSfDoaOnlyEffective:InvalidDoaType', ...
    'modelOpt.doaType must be ''angle'' or ''latlon''.');
end

if ~isfield(modelOpt, 'effectiveGainMode') || isempty(modelOpt.effectiveGainMode)
  modelOpt.effectiveGainMode = 'pilotModel';
else
  modelOpt.effectiveGainMode = char(lower(string(modelOpt.effectiveGainMode)));
end

if ~ismember(modelOpt.effectiveGainMode, {'pilotmodel', 'unit'})
  error('crbPilotSfDoaOnlyEffective:InvalidEffectiveGainMode', ...
    'modelOpt.effectiveGainMode must be ''pilotModel'' or ''unit''.');
end
if strcmp(modelOpt.effectiveGainMode, 'pilotmodel')
  modelOpt.effectiveGainMode = 'pilotModel';
end

if ~isfield(modelOpt, 'lightSpeed') || isempty(modelOpt.lightSpeed)
  modelOpt.lightSpeed = 299792458;
end
if ~isscalar(modelOpt.lightSpeed) || ~isfinite(modelOpt.lightSpeed) || modelOpt.lightSpeed <= 0
  error('crbPilotSfDoaOnlyEffective:InvalidLightSpeed', ...
    'modelOpt.lightSpeed must be a positive finite scalar.');
end

if ~isfield(modelOpt, 'timeInput') || isempty(modelOpt.timeInput)
  modelOpt.timeInput = sampleRate;
end
if ~isnumeric(modelOpt.timeInput) || any(~isfinite(modelOpt.timeInput(:))) || ...
    ~(isscalar(modelOpt.timeInput) || isvector(modelOpt.timeInput))
  error('crbPilotSfDoaOnlyEffective:InvalidTimeInput', ...
    'modelOpt.timeInput must be a scalar sample rate or a sample-time vector.');
end

if ~isfield(modelOpt, 'jacOpt') || isempty(modelOpt.jacOpt)
  modelOpt.jacOpt = struct();
end
if ~isstruct(modelOpt.jacOpt)
  error('crbPilotSfDoaOnlyEffective:InvalidJacOpt', ...
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
  error('crbPilotSfDoaOnlyEffective:RotMatCountMismatch', ...
    'The number of rotation matrices must match the number of arrays.');
end
if size(satPosEci, 2) ~= numSat || size(satVelEci, 2) ~= numSat
  error('crbPilotSfDoaOnlyEffective:SatCountMismatch', ...
    'scene.satPosEci and scene.satVelEci must match the number of arrays.');
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
end

model = struct();
model.array = arrayCell;
model.numSat = numSat;
model.wavelength = wavelength;
model.geom = geom;

end

function gainAbs = localBuildEffectivePilotGainAbs(pilotRow, fdVec, timeInput, gainMode)
%LOCALBUILDEFFECTIVEPILOTGAINABS Build no-Doppler projection gain.

numSat = numel(fdVec);
if strcmp(gainMode, 'unit')
  gainAbs = ones(numSat, 1);
  return;
end

pilot0 = repmat(reshape(pilotRow, 1, []), numSat, 1);
pilotFd = buildDopplerPilot(pilot0, fdVec(:), timeInput);
atomEnergy = sum(abs(pilot0).^2, 2);
gainAbs = abs(sum(conj(pilot0) .* pilotFd, 2) ./ atomEnergy);
gainAbs(~isfinite(gainAbs)) = NaN;

end

function pilotRow = localParsePilotWave(pilotWave)
%LOCALPARSEPILOTWAVE Normalize pilot waveform to a row vector.

if ~isvector(pilotWave) || isempty(pilotWave)
  error('crbPilotSfDoaOnlyEffective:InvalidPilotWave', ...
    'pilotWave must be a non-empty vector.');
end
pilotRow = reshape(pilotWave, 1, []);
if real(pilotRow * pilotRow') <= eps
  error('crbPilotSfDoaOnlyEffective:ZeroPilotEnergy', ...
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
  error('crbPilotSfDoaOnlyEffective:PerSatSizeMismatch', ...
    '%s must be a scalar or a vector with one entry per satellite.', nameText);
end
valueVec = reshape(valueIn, [], 1);

end

function arrayCell = localParseSceneArray(scene)
%LOCALPARSESCENEARRAY Convert scene.array to a 1xNs cell array.

if ~isfield(scene, 'array') || isempty(scene.array)
  error('crbPilotSfDoaOnlyEffective:MissingArrayField', ...
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
  error('crbPilotSfDoaOnlyEffective:InvalidArrayField', ...
    'scene.array must be a struct, struct array, or cell array.');
end

end

function rotMatCell = localParseRotMat(scene)
%LOCALPARSEROTMAT Convert scene.rotMat to a 1xNs cell array.

if ~isfield(scene, 'rotMat') || isempty(scene.rotMat)
  error('crbPilotSfDoaOnlyEffective:MissingRotMatField', ...
    'scene.rotMat is required.');
end

if isnumeric(scene.rotMat)
  if ~isequal(size(scene.rotMat), [3, 3])
    error('crbPilotSfDoaOnlyEffective:InvalidRotMatSize', ...
      'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  rotMatCell = {scene.rotMat};
  return;
end

if ~iscell(scene.rotMat) || isempty(scene.rotMat)
  error('crbPilotSfDoaOnlyEffective:InvalidRotMatType', ...
    'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

rotMatCell = reshape(scene.rotMat, 1, []);
for iSat = 1:numel(rotMatCell)
  if ~isnumeric(rotMatCell{iSat}) || ~isequal(size(rotMatCell{iSat}), [3, 3])
    error('crbPilotSfDoaOnlyEffective:InvalidRotMatCell', ...
      'Each scene.rotMat entry must be a 3x3 numeric matrix.');
  end
end

end

function [satPosEci, satVelEci] = localParseSatState(scene)
%LOCALPARSESATSTATE Validate satellite ECI state arrays.

if ~isfield(scene, 'satPosEci') || ~isfield(scene, 'satVelEci') || ...
    isempty(scene.satPosEci) || isempty(scene.satVelEci)
  error('crbPilotSfDoaOnlyEffective:MissingSatState', ...
    'scene.satPosEci and scene.satVelEci are required.');
end

satPosEci = scene.satPosEci;
satVelEci = scene.satVelEci;
if ~isnumeric(satPosEci) || ~isnumeric(satVelEci) || ...
    ~isequal(size(satPosEci), size(satVelEci)) || size(satPosEci, 1) ~= 3
  error('crbPilotSfDoaOnlyEffective:InvalidSatStateSize', ...
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
  error('crbPilotSfDoaOnlyEffective:MissingReferenceState', ...
    'scene.ref.posEci and scene.ref.velEci are required when scene.ref is provided.');
end
if ~isnumeric(refIn.posEci) || ~isequal(size(refIn.posEci), [3, 1])
  error('crbPilotSfDoaOnlyEffective:InvalidReferencePosSize', ...
    'scene.ref.posEci must have size 3x1.');
end
if ~isnumeric(refIn.velEci) || ~isequal(size(refIn.velEci), [3, 1])
  error('crbPilotSfDoaOnlyEffective:InvalidReferenceVelSize', ...
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
