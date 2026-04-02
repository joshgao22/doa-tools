function [jacobian, doaState] = buildDynDoaDopplerJacobian(model, doaParam, varargin)
%BUILDDYNDOADOPPLERJACOBIAN Build geometry Jacobians for dynamic DoA-Doppler models.
%
%Syntax:
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam)
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam, fdRef)
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam, fdRef, fdRate)
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam, jacOpt)
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam, fdRef, jacOpt)
%  jacobian = buildDynDoaDopplerJacobian(model, doaParam, fdRef, fdRate, jacOpt)
%  [jacobian, doaState] = buildDynDoaDopplerJacobian(...)
%
%Inputs:
%  model    - Dynamic DoA-Doppler geometry model structure
%             Common fields:
%               .doaType            : 'angle' or 'latlon'
%               .rotMat             : rotation matrices
%                                     - Nf x Ns cell, each 3x3
%                                     - 1 x Ns cell replicated to all frames
%                                     - 3x3 numeric for single-satellite use
%               .wavelength         : carrier wavelength in meters
%               .timeOffsetSec      : 1xNf or Nfx1 frame offsets in seconds
%
%             Optional common fields:
%               .steeringMode       : 'frozenRef' or 'framewise'
%                                     default: 'frozenRef'
%               .steeringRefFrameIdx: reference frame for frozen steering
%                                     default: round((Nf + 1) / 2)
%
%             angle mode fields:
%               .relSatVelEci       : 3xNsxNf relative satellite velocity
%                                     or
%               .satVelEci          : 3xNsxNf satellite velocity
%               .refVelEci          : 3xNf reference velocity
%
%             latlon mode fields:
%               .satPosEci          : 3xNsxNf satellite positions in ECI
%               .satVelEci          : 3xNsxNf satellite velocities in ECI
%               .refPosEci          : 3xNf reference positions in ECI
%               .refVelEci          : 3xNf reference velocities in ECI
%               .utcVec / .sceneUtcVec / .utc : frame UTC epochs used by
%                                               lla2eci
%
%  doaParam - 2x1 continuous geometry parameter vector
%             angle mode : [eciAz; eciEl]
%             latlon mode: [lat; lon]
%
%  fdRef    - (optional) scalar Doppler intercept in Hz
%
%  fdRate   - (optional) scalar Doppler rate in Hz/s
%             If fdRef is provided but fdRate is omitted, fdRate defaults
%             to zero.
%
%  jacOpt   - (optional) Jacobian options structure
%             .diffStep : finite-difference step for latlon mode
%                         scalar or 2x1 vector, default 1e-6
%
%Outputs:
%  jacobian - Jacobian structure
%             .doaType
%             .numSat
%             .numFrame
%             .numDoaParam
%             .paramNameDoa
%             .paramName
%             .steeringMode
%             .steeringRefFrameIdx
%             .timeOffsetSec
%             .localDoa      : 2x2xNsxNf local DoA Jacobians
%             .localDoaCell  : Nf x Ns cell, each 2x2
%             .deltaFd       : NsxNfx2 differential Doppler Jacobians
%             .fdRef         : NsxNf ones array when fdRef is provided
%             .fdRate        : NsxNf time-offset array when fdRef is provided
%             .fd            : NsxNfx4 total Doppler Jacobians for
%                              [doaParam; fdRef; fdRate]
%
%  doaState - Dynamic geometry state structure
%             .eciAngle
%             .eciDirection
%             .latlon
%             .userPosEci
%             .userVelEci
%             .deltaFd       : NsxNf
%             .localDoaArr   : 2xNsxNf
%             .localDoaCell  : Nf x Ns cell, each 2x1
%             .fdRef
%             .fdRate
%             .fdTrend       : 1xNf
%             .fd            : NsxNf total Doppler
%
%Description:
%  Builds the Jacobians required by dynamic pilot-based DoA-Doppler CRB and
%  MLE formulations under the model
%
%    fd(k,n) = fdRef + fdRate * timeOffsetSec(n) + deltaFd(k,n, doaParam).
%
%  The DoA part can either be frozen at one steering reference frame or be
%  updated frame by frame, controlled by model.steeringMode.
%
%Notes:
%  - In angle mode, the local DoA and differential Doppler Jacobians are
%    computed analytically.
%  - In latlon mode, the Jacobians are obtained using central finite
%    differences around the exact dynamic geometry map.
%  - The dynamic model is single-source only.
%
%See also:
%  buildDoaDopplerJacobian, doa2dir, doa2dirJacobian,
%  estimatorDoaDopplerMlePilotDynOpt, eciToAngleGrid, globalToLocalDoa

narginchk(2, 5);

validateattributes(model, {'struct'}, {'scalar'}, mfilename, 'model', 1);
validateattributes(doaParam, {'numeric'}, {'finite', 'real', 'vector', 'numel', 2}, ...
  mfilename, 'doaParam', 2);
doaParam = reshape(doaParam, [], 1);

[fdRef, fdRate, jacOpt] = localParseOptionalInput(varargin{:});
[doaType, paramNameDoa] = localParseDoaType(model);
timeOffsetSec = localParseTimeOffset(model);
numFrame = numel(timeOffsetSec);
[steeringMode, steeringRefFrameIdx] = localParseSteeringOpt(model, numFrame);
rotMatCell = localParseRotMat(model, numFrame);
numSat = size(rotMatCell, 2);

switch doaType
  case 'angle'
    doaState = localBuildAngleState(model, doaParam, rotMatCell, steeringMode, steeringRefFrameIdx);
    jacobian = localBuildAngleJacobian(model, doaParam, doaState, rotMatCell, ...
      steeringMode, steeringRefFrameIdx, paramNameDoa);

  case 'latlon'
    doaState = localBuildLatlonState(model, doaParam, rotMatCell, steeringMode, steeringRefFrameIdx);
    jacobian = localBuildLatlonJacobian(model, doaParam, doaState, rotMatCell, ...
      steeringMode, steeringRefFrameIdx, jacOpt, paramNameDoa);

  otherwise
    error('buildDynDoaDopplerJacobian:InvalidDoaType', ...
      'Unsupported doaType: %s.', doaType);
end

jacobian.doaType = doaType;
jacobian.numSat = numSat;
jacobian.numFrame = numFrame;
jacobian.numDoaParam = numel(doaParam);
jacobian.paramNameDoa = paramNameDoa;
jacobian.steeringMode = steeringMode;
jacobian.steeringRefFrameIdx = steeringRefFrameIdx;
jacobian.timeOffsetSec = timeOffsetSec;

if isempty(fdRef)
  jacobian.fdRef = [];
  jacobian.fdRate = [];
  jacobian.fd = [];
  jacobian.paramName = paramNameDoa;

  doaState.fdRef = [];
  doaState.fdRate = [];
  doaState.fdTrend = [];
  doaState.fd = [];
else
  validateattributes(fdRef, {'numeric'}, {'finite', 'real', 'scalar'}, ...
    mfilename, 'fdRef');
  validateattributes(fdRate, {'numeric'}, {'finite', 'real', 'scalar'}, ...
    mfilename, 'fdRate');

  jacobian.fdRef = ones(numSat, numFrame);
  jacobian.fdRate = repmat(reshape(timeOffsetSec, 1, []), numSat, 1);
  jacobian.fd = cat(3, jacobian.deltaFd, jacobian.fdRef, jacobian.fdRate);
  jacobian.paramName = [paramNameDoa, {'fdRef', 'fdRate'}];

  doaState.fdRef = fdRef;
  doaState.fdRate = fdRate;
  doaState.fdTrend = fdRef + fdRate * reshape(timeOffsetSec, 1, []);
  doaState.fd = doaState.deltaFd + repmat(doaState.fdTrend, numSat, 1);
end

end

function [fdRef, fdRate, jacOpt] = localParseOptionalInput(varargin)
%LOCALPARSEOPTIONALINPUT Parse optional fdRef, fdRate, and jacOpt inputs.

fdRef = [];
fdRate = [];
jacOpt = struct();

switch numel(varargin)
  case 0
    % keep defaults

  case 1
    if isstruct(varargin{1})
      jacOpt = varargin{1};
    else
      fdRef = varargin{1};
      fdRate = 0;
    end

  case 2
    if isstruct(varargin{2})
      fdRef = varargin{1};
      fdRate = 0;
      jacOpt = varargin{2};
    else
      fdRef = varargin{1};
      fdRate = varargin{2};
    end

  case 3
    fdRef = varargin{1};
    fdRate = varargin{2};
    jacOpt = varargin{3};

  otherwise
    error('buildDynDoaDopplerJacobian:TooManyInputs', ...
      'At most five input arguments are supported.');
end

if ~isempty(jacOpt) && ~isstruct(jacOpt)
  error('buildDynDoaDopplerJacobian:InvalidJacOpt', ...
    'jacOpt must be a structure.');
end

if ~isfield(jacOpt, 'diffStep') || isempty(jacOpt.diffStep)
  jacOpt.diffStep = 1e-6;
end

validateattributes(jacOpt.diffStep, {'numeric'}, {'finite', 'real', 'vector'}, ...
  mfilename, 'jacOpt.diffStep');

if isscalar(jacOpt.diffStep)
  jacOpt.diffStep = repmat(jacOpt.diffStep, 2, 1);
else
  if numel(jacOpt.diffStep) ~= 2
    error('buildDynDoaDopplerJacobian:InvalidDiffStep', ...
      'jacOpt.diffStep must be a scalar or a 2-element vector.');
  end
  jacOpt.diffStep = reshape(jacOpt.diffStep, 2, 1);
end

if any(jacOpt.diffStep <= 0)
  error('buildDynDoaDopplerJacobian:InvalidDiffStep', ...
    'jacOpt.diffStep must contain positive values.');
end

end

function [doaType, paramNameDoa] = localParseDoaType(model)
%LOCALPARSEDOATYPE Parse doaType and parameter names.

if ~isfield(model, 'doaType') || isempty(model.doaType)
  doaType = 'angle';
else
  doaType = lower(string(model.doaType));
  doaType = char(doaType);
end

switch doaType
  case 'angle'
    paramNameDoa = {'eciAz', 'eciEl'};
  case 'latlon'
    paramNameDoa = {'lat', 'lon'};
  otherwise
    error('buildDynDoaDopplerJacobian:InvalidDoaType', ...
      'model.doaType must be ''angle'' or ''latlon''.');
end

end

function timeOffsetSec = localParseTimeOffset(model)
%LOCALPARSETIMEOFFSET Normalize frame time offsets to a row vector.

if ~isfield(model, 'timeOffsetSec') || isempty(model.timeOffsetSec)
  error('buildDynDoaDopplerJacobian:MissingTimeOffset', ...
    'model.timeOffsetSec is required for dynamic Jacobian construction.');
end

validateattributes(model.timeOffsetSec, {'numeric'}, {'finite', 'real', 'vector', 'nonempty'}, ...
  mfilename, 'model.timeOffsetSec');
timeOffsetSec = reshape(model.timeOffsetSec, 1, []);

end

function [steeringMode, steeringRefFrameIdx] = localParseSteeringOpt(model, numFrame)
%LOCALPARSESTEERINGOPT Parse steering mode and reference frame index.

if ~isfield(model, 'steeringMode') || isempty(model.steeringMode)
  steeringMode = 'frozenref';
else
  steeringMode = lower(string(model.steeringMode));
  steeringMode = char(steeringMode);
end

if ~ismember(steeringMode, {'frozenref', 'framewise'})
  error('buildDynDoaDopplerJacobian:InvalidSteeringMode', ...
    'model.steeringMode must be ''frozenRef'' or ''framewise''.');
end

if ~isfield(model, 'steeringRefFrameIdx') || isempty(model.steeringRefFrameIdx)
  steeringRefFrameIdx = round((numFrame + 1) / 2);
else
  steeringRefFrameIdx = model.steeringRefFrameIdx;
end

if ~isscalar(steeringRefFrameIdx) || ~isnumeric(steeringRefFrameIdx) || ...
    ~isfinite(steeringRefFrameIdx) || mod(steeringRefFrameIdx, 1) ~= 0 || ...
    steeringRefFrameIdx < 1 || steeringRefFrameIdx > numFrame
  error('buildDynDoaDopplerJacobian:InvalidSteeringRefFrameIdx', ...
    'model.steeringRefFrameIdx must be an integer within [1, numFrame].');
end

end

function rotMatCell = localParseRotMat(model, numFrame)
%LOCALPARSEROTMAT Normalize dynamic rotation input to an Nf x Ns cell array.

if ~isfield(model, 'rotMat') || isempty(model.rotMat)
  error('buildDynDoaDopplerJacobian:MissingRotMat', ...
    'model.rotMat is required.');
end

rotMat = model.rotMat;

if isnumeric(rotMat)
  if ~isequal(size(rotMat), [3, 3])
    error('buildDynDoaDopplerJacobian:InvalidRotMatSize', ...
      'Numeric model.rotMat must be a 3x3 matrix.');
  end
  rotMatCell = repmat({rotMat}, numFrame, 1);
  return;
end

if ~iscell(rotMat) || isempty(rotMat)
  error('buildDynDoaDopplerJacobian:InvalidRotMatType', ...
    'model.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

cellSize = size(rotMat);
for iCell = 1:numel(rotMat)
  currentRot = rotMat{iCell};
  if ~isnumeric(currentRot) || ~isequal(size(currentRot), [3, 3])
    error('buildDynDoaDopplerJacobian:InvalidRotMatCell', ...
      'Each model.rotMat cell must contain one 3x3 numeric matrix.');
  end
end

if numel(cellSize) == 2 && cellSize(1) == numFrame
  rotMatCell = rotMat;
elseif numel(cellSize) == 2 && cellSize(1) == 1
  rotMatCell = repmat(reshape(rotMat, 1, []), numFrame, 1);
elseif numel(cellSize) == 2 && cellSize(2) == 1 && numFrame == 1
  rotMatCell = reshape(rotMat, 1, []);
else
  error('buildDynDoaDopplerJacobian:InvalidRotMatGrid', ...
    'model.rotMat must be Nf x Ns or 1 x Ns when numFrame = Nf.');
end

end

function doaState = localBuildAngleState(model, doaParam, rotMatCell, steeringMode, steeringRefFrameIdx)
%LOCALBUILDANGLESTATE Build frame-wise states for angle-mode parameterization.

localCheckWavelength(model);
relSatVelEci = localParseRelSatVel(model, size(rotMatCell, 2), size(rotMatCell, 1));

eciAngle = doaParam;
eciAngle(1) = mod(eciAngle(1), 2 * pi);
eciDirection = doa2dir(eciAngle);

localDoaArr = localBuildAngleLocalDoa(eciDirection, rotMatCell, steeringMode, steeringRefFrameIdx);
deltaFd = zeros(size(rotMatCell, 2), size(rotMatCell, 1));
for iFrame = 1:size(rotMatCell, 1)
  deltaFd(:, iFrame) = ((eciDirection.' * relSatVelEci(:, :, iFrame)) / model.wavelength).';
end

localDoaCell = localStateCellFromArray(localDoaArr);

doaState = struct();
doaState.doaParam = doaParam;
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = [];
doaState.userPosEci = [];
doaState.userVelEci = [];
doaState.deltaFd = deltaFd;
doaState.localDoaArr = localDoaArr;
doaState.localDoaCell = localDoaCell;

end

function doaState = localBuildLatlonState(model, doaParam, rotMatCell, steeringMode, steeringRefFrameIdx)
%LOCALBUILDLATLONSTATE Build frame-wise states for latlon parameterization.

localCheckWavelength(model);
[numSat, numFrame, satPosEci, satVelEci, refPosEci, refVelEci, utcVec] = localParseLatlonStateInput(model, size(rotMatCell, 2), size(rotMatCell, 1));

[userPosEci, userVelEci] = localLatlonToUserState(doaParam, utcVec);
eciDirection = userPosEci ./ vecnorm(userPosEci, 2, 1);
eciAngle = localDirectionToAngle(eciDirection);

localDoaArr = zeros(2, numSat, numFrame);
switch steeringMode
  case 'frozenref'
    refUserPos = userPosEci(:, steeringRefFrameIdx);
    for iSat = 1:numSat
      frozenDoa = globalToLocalDoa(refUserPos, satPosEci(:, iSat, steeringRefFrameIdx), ...
        rotMatCell{steeringRefFrameIdx, iSat});
      localDoaArr(:, iSat, :) = repmat(frozenDoa, 1, 1, numFrame);
    end

  case 'framewise'
    for iFrame = 1:numFrame
      currentUserPos = userPosEci(:, iFrame);
      for iSat = 1:numSat
        localDoaArr(:, iSat, iFrame) = globalToLocalDoa(currentUserPos, ...
          satPosEci(:, iSat, iFrame), rotMatCell{iFrame, iSat});
      end
    end

  otherwise
    error('buildDynDoaDopplerJacobian:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', steeringMode);
end

deltaFd = localComputeExactDeltaFd(satPosEci, satVelEci, refPosEci, refVelEci, ...
  userPosEci, userVelEci, model.wavelength);
localDoaCell = localStateCellFromArray(localDoaArr);

doaState = struct();
doaState.doaParam = doaParam;
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = doaParam;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFd = deltaFd;
doaState.localDoaArr = localDoaArr;
doaState.localDoaCell = localDoaCell;

end

function jacobian = localBuildAngleJacobian(model, doaParam, doaState, rotMatCell, steeringMode, steeringRefFrameIdx, paramNameDoa)
%LOCALBUILDANGLEJACOBIAN Build analytical Jacobians for angle mode.

numFrame = size(rotMatCell, 1);
numSat = size(rotMatCell, 2);
numParam = numel(doaParam);
relSatVelEci = localParseRelSatVel(model, numSat, numFrame);

[~, duDir] = doa2dirJacobian(doaParam);
duCell = {duDir.az, duDir.el};

localDoa = zeros(2, numParam, numSat, numFrame);
localDoaCell = cell(numFrame, numSat);
deltaFd = zeros(numSat, numFrame, numParam);

switch steeringMode
  case 'frozenref'
    refRotRow = rotMatCell(steeringRefFrameIdx, :);
    for iSat = 1:numSat
      currentJac = localBuildOneAngleJacobian(doaState.eciDirection, refRotRow{iSat}, duCell, iSat);
      for iFrame = 1:numFrame
        localDoa(:, :, iSat, iFrame) = currentJac;
        localDoaCell{iFrame, iSat} = currentJac;
      end
    end

  case 'framewise'
    for iFrame = 1:numFrame
      for iSat = 1:numSat
        currentJac = localBuildOneAngleJacobian(doaState.eciDirection, rotMatCell{iFrame, iSat}, duCell, iSat);
        localDoa(:, :, iSat, iFrame) = currentJac;
        localDoaCell{iFrame, iSat} = currentJac;
      end
    end

  otherwise
    error('buildDynDoaDopplerJacobian:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', steeringMode);
end

for iFrame = 1:numFrame
  relVel = relSatVelEci(:, :, iFrame);
  deltaFd(:, iFrame, 1) = (duCell{1}.' * relVel / model.wavelength).';
  deltaFd(:, iFrame, 2) = (duCell{2}.' * relVel / model.wavelength).';
end

jacobian = struct();
jacobian.paramNameDoa = paramNameDoa;
jacobian.localDoa = localDoa;
jacobian.localDoaCell = localDoaCell;
jacobian.deltaFd = deltaFd;

end

function jacobian = localBuildLatlonJacobian(model, doaParam, doaState, rotMatCell, steeringMode, steeringRefFrameIdx, jacOpt, paramNameDoa)
%LOCALBUILDLATLONJACOBIAN Build finite-difference Jacobians for latlon mode.

numFrame = size(rotMatCell, 1);
numSat = size(rotMatCell, 2);
numParam = numel(doaParam);

localDoa = zeros(2, numParam, numSat, numFrame);
localDoaCell = cell(numFrame, numSat);
deltaFd = zeros(numSat, numFrame, numParam);

for iParam = 1:numParam
  stepVal = jacOpt.diffStep(iParam);
  paramPlus = doaParam;
  paramMinus = doaParam;
  paramPlus(iParam) = paramPlus(iParam) + stepVal;
  paramMinus(iParam) = paramMinus(iParam) - stepVal;

  statePlus = localBuildLatlonState(model, paramPlus, rotMatCell, steeringMode, steeringRefFrameIdx);
  stateMinus = localBuildLatlonState(model, paramMinus, rotMatCell, steeringMode, steeringRefFrameIdx);

  deltaFd(:, :, iParam) = (statePlus.deltaFd - stateMinus.deltaFd) / (2 * stepVal);

  for iFrame = 1:numFrame
    for iSat = 1:numSat
      azPlus = statePlus.localDoaArr(1, iSat, iFrame);
      azMinus = stateMinus.localDoaArr(1, iSat, iFrame);
      elPlus = statePlus.localDoaArr(2, iSat, iFrame);
      elMinus = stateMinus.localDoaArr(2, iSat, iFrame);

      dAz = localAngleDiff(azPlus, azMinus) / (2 * stepVal);
      dEl = (elPlus - elMinus) / (2 * stepVal);
      localDoa(:, iParam, iSat, iFrame) = [dAz; dEl];
    end
  end
end

for iFrame = 1:numFrame
  for iSat = 1:numSat
    localDoaCell{iFrame, iSat} = localDoa(:, :, iSat, iFrame);
  end
end

jacobian = struct();
jacobian.paramNameDoa = paramNameDoa;
jacobian.localDoa = localDoa;
jacobian.localDoaCell = localDoaCell;
jacobian.deltaFd = deltaFd;

if isempty(doaState)
  error('buildDynDoaDopplerJacobian:UnexpectedState', ...
    'Internal state must not be empty.');
end

end

function localDoaArr = localBuildAngleLocalDoa(eciDirection, rotMatCell, steeringMode, steeringRefFrameIdx)
%LOCALBUILDANGLELOCALDOA Build frame-wise local DoA array in angle mode.

numFrame = size(rotMatCell, 1);
numSat = size(rotMatCell, 2);
localDoaArr = zeros(2, numSat, numFrame);

switch steeringMode
  case 'frozenref'
    refRotRow = rotMatCell(steeringRefFrameIdx, :);
    localDoaRef = eciToAngleGrid(eciDirection, refRotRow);
    if ~iscell(localDoaRef)
      localDoaRef = {localDoaRef};
    end

    for iFrame = 1:numFrame
      for iSat = 1:numSat
        localDoaArr(:, iSat, iFrame) = localDoaRef{iSat};
      end
    end

  case 'framewise'
    for iFrame = 1:numFrame
      localDoaCell = eciToAngleGrid(eciDirection, rotMatCell(iFrame, :));
      if ~iscell(localDoaCell)
        localDoaCell = {localDoaCell};
      end
      for iSat = 1:numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end

  otherwise
    error('buildDynDoaDopplerJacobian:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', steeringMode);
end

end

function currentJac = localBuildOneAngleJacobian(eciDirection, rotMat, duCell, satIdx)
%LOCALBUILDONEANGLEJACOBIAN Build local DoA Jacobian for one satellite.

uLocal = rotMat' * eciDirection;
xCoord = uLocal(1);
yCoord = uLocal(2);
zCoord = uLocal(3);
xyPow = xCoord^2 + yCoord^2;
xyNorm = sqrt(max(xyPow, eps));

if xyPow <= 1e-14
  warning('buildDynDoaDopplerJacobian:AzimuthPole', ...
    'Local azimuth Jacobian is poorly conditioned near elevation pole for satellite %d.', satIdx);
end

currentJac = zeros(2, numel(duCell));
for iParam = 1:numel(duCell)
  duLocal = rotMat' * duCell{iParam};
  dx = duLocal(1);
  dy = duLocal(2);
  dz = duLocal(3);

  currentJac(1, iParam) = (xCoord * dy - yCoord * dx) / max(xyPow, eps);
  currentJac(2, iParam) = dz / xyNorm;
end

end

function relSatVelEci = localParseRelSatVel(model, numSat, numFrame)
%LOCALPARSERELSATVEL Parse dynamic relative satellite velocities.

if isfield(model, 'relSatVelEci') && ~isempty(model.relSatVelEci)
  relSatVelEci = model.relSatVelEci;
  if ~isnumeric(relSatVelEci) || ~isequal(size(relSatVelEci), [3, numSat, numFrame])
    error('buildDynDoaDopplerJacobian:InvalidRelSatVelSize', ...
      'model.relSatVelEci must have size 3x%d x %d.', numSat, numFrame);
  end
  return;
end

if ~isfield(model, 'satVelEci') || isempty(model.satVelEci) || ...
    ~isfield(model, 'refVelEci') || isempty(model.refVelEci)
  error('buildDynDoaDopplerJacobian:MissingRelSatVel', ...
    'Angle mode requires model.relSatVelEci, or model.satVelEci and model.refVelEci.');
end

satVelEci = model.satVelEci;
refVelEci = model.refVelEci;
if ~isnumeric(satVelEci) || ~isequal(size(satVelEci), [3, numSat, numFrame])
  error('buildDynDoaDopplerJacobian:InvalidSatVelSize', ...
    'model.satVelEci must have size 3x%d x %d.', numSat, numFrame);
end
if ~isnumeric(refVelEci) || ~isequal(size(refVelEci), [3, numFrame])
  error('buildDynDoaDopplerJacobian:InvalidRefVelSize', ...
    'model.refVelEci must have size 3x%d.', numFrame);
end

relSatVelEci = satVelEci - reshape(refVelEci, 3, 1, numFrame);

end

function [numSat, numFrame, satPosEci, satVelEci, refPosEci, refVelEci, utcVec] = ...
  localParseLatlonStateInput(model, numSatExp, numFrameExp)
%LOCALPARSELATLONSTATEINPUT Parse dynamic latlon geometry inputs.

requiredField = {'satPosEci', 'satVelEci', 'refPosEci', 'refVelEci'};
for iField = 1:numel(requiredField)
  fieldName = requiredField{iField};
  if ~isfield(model, fieldName) || isempty(model.(fieldName))
    error('buildDynDoaDopplerJacobian:MissingLatlonField', ...
      'model.%s is required for latlon mode.', fieldName);
  end
end

satPosEci = model.satPosEci;
satVelEci = model.satVelEci;
refPosEci = model.refPosEci;
refVelEci = model.refVelEci;

if ~isnumeric(satPosEci) || ~isequal(size(satPosEci), [3, numSatExp, numFrameExp])
  error('buildDynDoaDopplerJacobian:InvalidSatPosSize', ...
    'model.satPosEci must have size 3x%d x %d.', numSatExp, numFrameExp);
end
if ~isnumeric(satVelEci) || ~isequal(size(satVelEci), [3, numSatExp, numFrameExp])
  error('buildDynDoaDopplerJacobian:InvalidSatVelSize', ...
    'model.satVelEci must have size 3x%d x %d.', numSatExp, numFrameExp);
end
if ~isnumeric(refPosEci) || ~isequal(size(refPosEci), [3, numFrameExp])
  error('buildDynDoaDopplerJacobian:InvalidRefPosSize', ...
    'model.refPosEci must have size 3x%d.', numFrameExp);
end
if ~isnumeric(refVelEci) || ~isequal(size(refVelEci), [3, numFrameExp])
  error('buildDynDoaDopplerJacobian:InvalidRefVelSize', ...
    'model.refVelEci must have size 3x%d.', numFrameExp);
end

utcVec = localGetUtcVec(model, numFrameExp);
numSat = numSatExp;
numFrame = numFrameExp;

end

function utcVec = localGetUtcVec(model, numFrame)
%LOCALGETUTCVEC Extract frame UTC epochs from common field names.

utcVec = [];

candidateField = {'utcVec', 'sceneUtcVec', 'sceneUtc', 'utc'};
for iField = 1:numel(candidateField)
  fieldName = candidateField{iField};
  if isfield(model, fieldName) && ~isempty(model.(fieldName))
    utcVec = model.(fieldName);
    break;
  end
end

if isempty(utcVec) && isfield(model, 'scene') && isstruct(model.scene)
  for iField = 1:numel(candidateField)
    fieldName = candidateField{iField};
    if isfield(model.scene, fieldName) && ~isempty(model.scene.(fieldName))
      utcVec = model.scene.(fieldName);
      break;
    end
  end
end

if isempty(utcVec)
  error('buildDynDoaDopplerJacobian:MissingUtcVec', ...
    'utcVec is required for latlon mode.');
end

if isdatetime(utcVec)
  if ~isvector(utcVec) || numel(utcVec) ~= numFrame
    error('buildDynDoaDopplerJacobian:InvalidUtcVecLength', ...
      'utcVec must be a datetime vector with numFrame entries.');
  end
  utcVec = reshape(utcVec, 1, []);

elseif isnumeric(utcVec)
  if ~(ismatrix(utcVec) && size(utcVec, 2) == 6 && size(utcVec, 1) == numFrame)
    error('buildDynDoaDopplerJacobian:InvalidUtcMatSize', ...
      'Numeric utcVec must have size numFrame x 6.');
  end

else
  error('buildDynDoaDopplerJacobian:InvalidUtcVecType', ...
    'utcVec must be a datetime vector or an numFrame x 6 datevec matrix.');
end

end

function localDoaCell = localStateCellFromArray(localDoaArr)
%LOCALSTATECELLFROMARRAY Convert 2xNsxNf local DoA array to Nf x Ns cells.

numSat = size(localDoaArr, 2);
numFrame = size(localDoaArr, 3);
localDoaCell = cell(numFrame, numSat);
for iFrame = 1:numFrame
  for iSat = 1:numSat
    localDoaCell{iFrame, iSat} = localDoaArr(:, iSat, iFrame);
  end
end

end

function [userPosEci, userVelEci] = localLatlonToUserState(latlon, utcVec)
%LOCALLATLONTOUSERSTATE Convert one ground point [lat; lon] into ECI states.

numFrame = localGetNumFrameFromUtc(utcVec);
llaMat = [repmat(latlon(1).', numFrame, 1), ...
          repmat(latlon(2).', numFrame, 1), ...
          zeros(numFrame, 1)];

if isdatetime(utcVec)
  utcMat = datevec(utcVec);
else
  utcMat = utcVec;
end

userPosEci = lla2eci(llaMat, utcMat).';
userVelEci = localComputeEarthFixedVel(userPosEci);

end

function numFrame = localGetNumFrameFromUtc(utcVec)
%LOCALGETNUMFRAMEFROMUTC Get the number of frames from UTC input.

if isdatetime(utcVec)
  numFrame = numel(utcVec);
else
  numFrame = size(utcVec, 1);
end

end

function velEci = localComputeEarthFixedVel(posEci)
%LOCALCOMPUTEEARTHFIXEDVEL Compute ECI velocity of Earth-fixed points.

omegaEarth = 7.2921150e-5;
omegaVec = repmat([0; 0; omegaEarth], 1, size(posEci, 2));
velEci = cross(omegaVec, posEci, 1);

end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2 * pi); ...
         asin(max(min(direction(3, :), 1), -1))];

end

function deltaFd = localComputeExactDeltaFd(satPosEci, satVelEci, refPosEci, refVelEci, ...
  userPosEci, userVelEci, wavelength)
%LOCALCOMPUTEEXACTDELTAFD Compute exact dynamic differential Doppler.

numSat = size(satPosEci, 2);
numFrame = size(satPosEci, 3);
deltaFd = zeros(numSat, numFrame);

for iFrame = 1:numFrame
  currentUserPos = userPosEci(:, iFrame);
  currentUserVel = userVelEci(:, iFrame);

  refLos = refPosEci(:, iFrame) - currentUserPos;
  refRange = norm(refLos);
  if refRange <= 0
    error('buildDynDoaDopplerJacobian:InvalidReferenceRange', ...
      'Reference receiver coincides with a candidate user position.');
  end
  refLosUnit = refLos / refRange;
  refRelVel = refVelEci(:, iFrame) - currentUserVel;
  refRangeRate = refLosUnit' * refRelVel;
  refFd = -refRangeRate / wavelength;

  for iSat = 1:numSat
    satLos = satPosEci(:, iSat, iFrame) - currentUserPos;
    satRange = norm(satLos);
    if satRange <= 0
      error('buildDynDoaDopplerJacobian:InvalidSatRange', ...
        'Satellite %d coincides with a candidate user position at frame %d.', ...
        iSat, iFrame);
    end

    satLosUnit = satLos / satRange;
    satRelVel = satVelEci(:, iSat, iFrame) - currentUserVel;
    satRangeRate = satLosUnit' * satRelVel;
    satFd = -satRangeRate / wavelength;

    deltaFd(iSat, iFrame) = satFd - refFd;
  end
end

end

function localCheckWavelength(model)
%LOCALCHECKWAVELENGTH Validate wavelength field.

if ~isfield(model, 'wavelength') || isempty(model.wavelength)
  error('buildDynDoaDopplerJacobian:MissingWavelength', ...
    'model.wavelength is required.');
end
validateattributes(model.wavelength, {'numeric'}, {'finite', 'real', 'scalar', 'positive'}, ...
  mfilename, 'model.wavelength');

end

function diffVal = localAngleDiff(angleLeft, angleRight)
%LOCALANGLEDIFF Wrapped difference angleLeft - angleRight.

diffVal = angle(exp(1j * (angleLeft - angleRight)));

end
