function [jacobian, doaState] = buildDoaDopplerJacobian(model, doaParam, varargin)
%BUILDDOADOPPLERJACOBIAN Build geometry state Jacobians for DoA-Doppler models.
%
%Syntax:
%  jacobian = buildDoaDopplerJacobian(model, doaParam)
%  jacobian = buildDoaDopplerJacobian(model, doaParam, fdRef)
%  jacobian = buildDoaDopplerJacobian(model, doaParam, jacOpt)
%  jacobian = buildDoaDopplerJacobian(model, doaParam, fdRef, jacOpt)
%  [jacobian, doaState] = buildDoaDopplerJacobian(...)
%
%Inputs:
%  model    - DoA-Doppler geometry model structure
%             Required fields depend on model.doaType:
%
%             Common fields:
%               .doaType        : 'angle' or 'latlon'
%               .rotMat         : 3x3 matrix or 1xNs cell of 3x3 matrices
%               .wavelength     : carrier wavelength in meters
%
%             angle mode fields:
%               .relSatVelEci   : 3xNs relative satellite velocity matrix
%
%             latlon mode fields:
%               .satPosEci      : 3xNs satellite positions in ECI
%               .satVelEci      : 3xNs satellite velocities in ECI
%               .ref.posEci     : 3x1 reference position in ECI
%               .ref.velEci     : 3x1 reference velocity in ECI
%               .sceneUtc / .utc / .scene.utc : UTC epoch for lla2eci
%
%  doaParam - 2x1 continuous geometry parameter vector
%             angle mode : [eciAz; eciEl]
%             latlon mode: [lat; lon]
%
%  fdRef    - (optional) scalar reference Doppler in Hz
%             When provided, jacobian.fd and doaState.fd are also returned.
%
%  jacOpt   - (optional) Jacobian options structure
%             .diffStep : finite-difference step for latlon mode
%                         scalar or 2x1 vector, default 1e-6
%
%Outputs:
%  jacobian - Jacobian structure
%             .doaType        : copied doa type
%             .numSat         : number of satellites
%             .numDoaParam    : number of geometry parameters
%             .paramNameDoa   : {'eciAz','eciEl'} or {'lat','lon'}
%             .paramName      : geometry names, plus 'fdRef' when available
%             .localDoa       : 2x2xNs local DoA Jacobians
%             .localDoaCell   : 1xNs cells, each 2x2
%             .deltaFd        : Nsx2 Jacobian of differential Doppler
%             .fdRef          : Nsx1 ones vector when fdRef is provided
%             .fd             : Nsx3 Jacobian of total Doppler [doaParam; fdRef]
%
%  doaState - Geometry state structure
%             .eciAngle
%             .eciDirection
%             .latlon
%             .userPosEci
%             .userVelEci
%             .deltaFd
%             .localDoaCell
%             .localDoaMat
%             .fdRef          : copied fdRef when provided
%             .fd             : fdRef + deltaFd when provided
%
%Description:
%  Builds the Jacobians needed by joint DoA-Doppler estimators and CRB
%  derivations. For angle mode, the Jacobians are computed analytically.
%  For latlon mode, the local DoA and differential Doppler Jacobians are
%  obtained using central finite differences around the exact geometry map.
%
%Notes:
%  - The local DoA Jacobian of satellite k is stored as
%        jacobian.localDoa(:, :, k)
%    with rows [az; el] and columns corresponding to doaParam.
%
%  - The differential Doppler follows the same convention as the existing
%    estimators:
%        fd(k) = fdRef + deltaFd(k).
%
%  - In angle mode, deltaFd(k) is modeled as
%        deltaFd(k) = u.' * relSatVelEci(:, k) / wavelength.
%
%See also:
%  doa2dir, doa2dirJacobian, eciToAngleGrid, globalToLocalDoa

narginchk(2, 4);

validateattributes(model, {'struct'}, {'scalar'}, mfilename, 'model', 1);
validateattributes(doaParam, {'numeric'}, {'finite', 'real', 'vector', 'numel', 2}, ...
  mfilename, 'doaParam', 2);
doaParam = reshape(doaParam, [], 1);

[fdRef, jacOpt] = localParseOptionalInput(varargin{:});
rotMatCell = localParseRotMat(model);
numSat = numel(rotMatCell);

[doaType, paramNameDoa] = localParseDoaType(model);
numDoaParam = numel(doaParam);

switch doaType
  case 'angle'
    doaState = localBuildAngleState(model, doaParam, rotMatCell);
    jacobian = localBuildAngleJacobian(model, doaParam, doaState, rotMatCell, paramNameDoa);

  case 'latlon'
    doaState = localBuildLatlonState(model, doaParam, rotMatCell);
    jacobian = localBuildLatlonJacobian(model, doaParam, doaState, rotMatCell, jacOpt, paramNameDoa);

  otherwise
    error('buildDoaDopplerJacobian:InvalidDoaType', ...
      'Unsupported doaType: %s.', doaType);
end

jacobian.doaType = doaType;
jacobian.numSat = numSat;
jacobian.numDoaParam = numDoaParam;
jacobian.paramNameDoa = paramNameDoa;

if isempty(fdRef)
  jacobian.fdRef = [];
  jacobian.fd = [];
  jacobian.paramName = paramNameDoa;
  doaState.fdRef = [];
  doaState.fd = [];
else
  validateattributes(fdRef, {'numeric'}, {'finite', 'real', 'scalar'}, ...
    mfilename, 'fdRef');
  jacobian.fdRef = ones(numSat, 1);
  jacobian.fd = [jacobian.deltaFd, jacobian.fdRef];
  jacobian.paramName = [paramNameDoa, {'fdRef'}];

  doaState.fdRef = fdRef;
  doaState.fd = fdRef + doaState.deltaFd(:);
end

end

function [fdRef, jacOpt] = localParseOptionalInput(varargin)
%LOCALPARSEOPTIONALINPUT Parse optional fdRef and jacOpt inputs.

fdRef = [];
jacOpt = struct();

if isempty(varargin)
  return;
elseif numel(varargin) == 1
  if isstruct(varargin{1})
    jacOpt = varargin{1};
  else
    fdRef = varargin{1};
  end
elseif numel(varargin) == 2
  fdRef = varargin{1};
  jacOpt = varargin{2};
else
  error('buildDoaDopplerJacobian:TooManyInputs', ...
    'At most four input arguments are supported.');
end

if ~isempty(jacOpt) && ~isstruct(jacOpt)
  error('buildDoaDopplerJacobian:InvalidJacOpt', ...
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
    error('buildDoaDopplerJacobian:InvalidDiffStep', ...
      'jacOpt.diffStep must be a scalar or a 2-element vector.');
  end
  jacOpt.diffStep = reshape(jacOpt.diffStep, 2, 1);
end

if any(jacOpt.diffStep <= 0)
  error('buildDoaDopplerJacobian:InvalidDiffStep', ...
    'jacOpt.diffStep must contain positive values.');
end

end

function rotMatCell = localParseRotMat(model)
%LOCALPARSEROTMAT Normalize rotation input to a 1xNs cell array.

if ~isfield(model, 'rotMat') || isempty(model.rotMat)
  error('buildDoaDopplerJacobian:MissingRotMat', ...
    'model.rotMat is required.');
end

rotMat = model.rotMat;
if isnumeric(rotMat)
  if ~isequal(size(rotMat), [3, 3])
    error('buildDoaDopplerJacobian:InvalidRotMatSize', ...
      'model.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  rotMatCell = {rotMat};

elseif iscell(rotMat)
  if isempty(rotMat)
    error('buildDoaDopplerJacobian:EmptyRotMatCell', ...
      'model.rotMat must not be empty.');
  end

  rotMatCell = reshape(rotMat, 1, []);
  for iSat = 1:numel(rotMatCell)
    currentRot = rotMatCell{iSat};
    if ~isnumeric(currentRot) || ~isequal(size(currentRot), [3, 3])
      error('buildDoaDopplerJacobian:InvalidRotMatCell', ...
        'Each model.rotMat cell must contain one 3x3 numeric matrix.');
    end
  end

else
  error('buildDoaDopplerJacobian:InvalidRotMatType', ...
    'model.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
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
    error('buildDoaDopplerJacobian:InvalidDoaType', ...
      'model.doaType must be ''angle'' or ''latlon''.');
end

end

function doaState = localBuildAngleState(model, doaParam, rotMatCell)
%LOCALBUILDANGLESTATE Build geometry state for ECI angle parameters.

localCheckWavelength(model);
numSat = numel(rotMatCell);

if ~isfield(model, 'relSatVelEci') || isempty(model.relSatVelEci)
  error('buildDoaDopplerJacobian:MissingRelSatVel', ...
    'model.relSatVelEci is required for angle mode.');
end
if ~isequal(size(model.relSatVelEci), [3, numSat])
  error('buildDoaDopplerJacobian:InvalidRelSatVelSize', ...
    'model.relSatVelEci must have size 3x%d.', numSat);
end

eciAngle = doaParam;
eciDirection = doa2dir(eciAngle);
localDoaCell = eciToAngleGrid(eciDirection, rotMatCell);
localDoaMat = localCellToMatrix(localDoaCell);
deltaFd = (eciDirection.' * model.relSatVelEci) / model.wavelength;

doaState = struct();
doaState.doaParam = doaParam;
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = [];
doaState.userPosEci = [];
doaState.userVelEci = [];
doaState.deltaFd = deltaFd(:);
doaState.localDoaCell = localDoaCell;
doaState.localDoaMat = localDoaMat;

end

function doaState = localBuildLatlonState(model, doaParam, rotMatCell)
%LOCALBUILDLATLONSTATE Build geometry state for [lat; lon] parameters.

localCheckWavelength(model);
numSat = numel(rotMatCell);
[sceneUtc, hasUtc] = localGetSceneUtc(model);
if ~hasUtc
  error('buildDoaDopplerJacobian:MissingUtc', ...
    'scene.utc, model.sceneUtc, or model.utc is required for latlon mode.');
end

if ~isfield(model, 'satPosEci') || ~isequal(size(model.satPosEci), [3, numSat])
  error('buildDoaDopplerJacobian:InvalidSatPosSize', ...
    'model.satPosEci must have size 3x%d.', numSat);
end
if ~isfield(model, 'satVelEci') || ~isequal(size(model.satVelEci), [3, numSat])
  error('buildDoaDopplerJacobian:InvalidSatVelSize', ...
    'model.satVelEci must have size 3x%d.', numSat);
end
if ~isfield(model, 'ref') || ~isstruct(model.ref) || ...
    ~isfield(model.ref, 'posEci') || ~isfield(model.ref, 'velEci')
  error('buildDoaDopplerJacobian:MissingReferenceState', ...
    'model.ref.posEci and model.ref.velEci are required for latlon mode.');
end

utcMat = localExpandUtc(sceneUtc, 1);
userPosEci = lla2eci([doaParam(1), doaParam(2), 0], utcMat).';
userVelEci = localComputeEarthFixedVel(userPosEci);

eciDirection = userPosEci - model.ref.posEci;
eciDirection = eciDirection / norm(eciDirection);
eciAngle = localDirectionToAngle(eciDirection);

localDoaCell = cell(1, numSat);
for iSat = 1:numSat
  localDoaCell{iSat} = globalToLocalDoa(userPosEci, model.satPosEci(:, iSat), rotMatCell{iSat});
end
localDoaMat = localCellToMatrix(localDoaCell);
deltaFd = localComputeExactDeltaFd(model, userPosEci, userVelEci);

doaState = struct();
doaState.doaParam = doaParam;
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = doaParam;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFd = deltaFd(:);
doaState.localDoaCell = localDoaCell;
doaState.localDoaMat = localDoaMat;

end

function jacobian = localBuildAngleJacobian(model, doaParam, doaState, rotMatCell, paramNameDoa)
%LOCALBUILDANGLEJACOBIAN Build analytical Jacobians for angle mode.

numSat = numel(rotMatCell);
numParam = numel(doaParam);

[~, duDir] = doa2dirJacobian(doaParam);
duCell = {duDir.az, duDir.el};

localDoa = zeros(2, numParam, numSat);
localDoaCell = cell(1, numSat);
deltaFd = zeros(numSat, numParam);

for iSat = 1:numSat
  currentRot = rotMatCell{iSat};
  uLocal = currentRot' * doaState.eciDirection;

  xCoord = uLocal(1);
  yCoord = uLocal(2);
  zCoord = uLocal(3);
  xyPow = xCoord^2 + yCoord^2;
  xyNorm = sqrt(max(xyPow, eps));

  if xyPow <= 1e-14
    warning('buildDoaDopplerJacobian:AzimuthPole', ...
      'Local azimuth Jacobian is poorly conditioned near elevation pole for satellite %d.', iSat);
  end

  currentJac = zeros(2, numParam);
  for iParam = 1:numParam
    duLocal = currentRot' * duCell{iParam};
    dx = duLocal(1);
    dy = duLocal(2);
    dz = duLocal(3);

    currentJac(1, iParam) = (xCoord * dy - yCoord * dx) / max(xyPow, eps);
    currentJac(2, iParam) = dz / xyNorm;
  end

  localDoa(:, :, iSat) = currentJac;
  localDoaCell{iSat} = currentJac;
end

relSatVelEci = model.relSatVelEci;
deltaFd(:, 1) = (duCell{1}.' * relSatVelEci / model.wavelength).';
deltaFd(:, 2) = (duCell{2}.' * relSatVelEci / model.wavelength).';

jacobian = struct();
jacobian.paramNameDoa = paramNameDoa;
jacobian.localDoa = localDoa;
jacobian.localDoaCell = localDoaCell;
jacobian.deltaFd = deltaFd;

end

function jacobian = localBuildLatlonJacobian(model, doaParam, doaState, rotMatCell, jacOpt, paramNameDoa)
%LOCALBUILDLATLONJACOBIAN Build finite-difference Jacobians for latlon mode.

numSat = numel(rotMatCell);
numParam = numel(doaParam);
localDoa = zeros(2, numParam, numSat);
localDoaCell = cell(1, numSat);
deltaFd = zeros(numSat, numParam);

for iParam = 1:numParam
  stepVal = jacOpt.diffStep(iParam);
  paramPlus = doaParam;
  paramMinus = doaParam;
  paramPlus(iParam) = paramPlus(iParam) + stepVal;
  paramMinus(iParam) = paramMinus(iParam) - stepVal;

  statePlus = localBuildLatlonState(model, paramPlus, rotMatCell);
  stateMinus = localBuildLatlonState(model, paramMinus, rotMatCell);

  deltaFd(:, iParam) = (statePlus.deltaFd - stateMinus.deltaFd) / (2 * stepVal);

  for iSat = 1:numSat
    azPlus = statePlus.localDoaMat(1, iSat);
    azMinus = stateMinus.localDoaMat(1, iSat);
    elPlus = statePlus.localDoaMat(2, iSat);
    elMinus = stateMinus.localDoaMat(2, iSat);

    dAz = localAngleDiff(azPlus, azMinus) / (2 * stepVal);
    dEl = (elPlus - elMinus) / (2 * stepVal);

    localDoa(:, iParam, iSat) = [dAz; dEl];
  end
end

for iSat = 1:numSat
  localDoaCell{iSat} = localDoa(:, :, iSat);
end

jacobian = struct();
jacobian.paramNameDoa = paramNameDoa;
jacobian.localDoa = localDoa;
jacobian.localDoaCell = localDoaCell;
jacobian.deltaFd = deltaFd;

% Keep the current state as an explicit dependency to avoid unused-input
% warnings in code analyzers.
if isempty(doaState)
  error('buildDoaDopplerJacobian:UnexpectedState', ...
    'Internal state must not be empty.');
end

end

function localDoaMat = localCellToMatrix(localDoaCell)
%LOCALCELLTOMATRIX Convert a 1xNs local DoA cell array to 2xNs matrix.

numSat = numel(localDoaCell);
localDoaMat = zeros(2, numSat);
for iSat = 1:numSat
  localDoaMat(:, iSat) = localDoaCell{iSat};
end

end

function localCheckWavelength(model)
%LOCALCHECKWAVELENGTH Validate wavelength field.

if ~isfield(model, 'wavelength') || isempty(model.wavelength)
  error('buildDoaDopplerJacobian:MissingWavelength', ...
    'model.wavelength is required.');
end
validateattributes(model.wavelength, {'numeric'}, {'finite', 'real', 'scalar', 'positive'}, ...
  mfilename, 'model.wavelength');

end

function [sceneUtc, hasUtc] = localGetSceneUtc(model)
%LOCALGETSCENEUTC Extract UTC epoch from common field names.

sceneUtc = [];
hasUtc = false;

if isfield(model, 'sceneUtc') && ~isempty(model.sceneUtc)
  sceneUtc = model.sceneUtc;
  hasUtc = true;
elseif isfield(model, 'utc') && ~isempty(model.utc)
  sceneUtc = model.utc;
  hasUtc = true;
elseif isfield(model, 'scene') && isstruct(model.scene) && ...
    isfield(model.scene, 'utc') && ~isempty(model.scene.utc)
  sceneUtc = model.scene.utc;
  hasUtc = true;
end

end


function utcMat = localExpandUtc(utc, numPoint)
%LOCALEXPUTC Expand UTC input to an numPoint-by-6 numeric date matrix.

if isdatetime(utc)
  utc = datevec(utc);
end

if isnumeric(utc) && isequal(size(utc), [1, 6])
  utcMat = repmat(utc, numPoint, 1);
elseif isnumeric(utc) && size(utc, 2) == 6 && size(utc, 1) == numPoint
  utcMat = utc;
else
  error('buildDoaDopplerJacobian:InvalidUtc', ...
    'utc must be a datetime scalar, 1x6 date vector, or %d-by-6 date matrix.', numPoint);
end

end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit direction to [az; el].

direction = direction / norm(direction);
zCoord = max(min(direction(3), 1), -1);
angle = [mod(atan2(direction(2), direction(1)), 2*pi); asin(zCoord)];

end

function velEci = localComputeEarthFixedVel(posEci)
%LOCALCOMPUTEEARTHFIXEDVEL Compute ECI velocity of an Earth-fixed point.

omegaEarth = 7.2921150e-5;
omegaVec = [0; 0; omegaEarth];
velEci = cross(omegaVec, posEci);

end

function deltaFd = localComputeExactDeltaFd(model, userPosEci, userVelEci)
%LOCALCOMPUTEEXACTDELTAFD Compute exact geometric differential Doppler.

numSat = size(model.satPosEci, 2);
deltaFd = zeros(numSat, 1);

refLos = model.ref.posEci - userPosEci;
refRange = norm(refLos);
if refRange <= 0
  error('buildDoaDopplerJacobian:InvalidReferenceRange', ...
    'Reference receiver coincides with the candidate user position.');
end

refLosUnit = refLos / refRange;
refRelVel = model.ref.velEci - userVelEci;
refRangeRate = refLosUnit' * refRelVel;
refFd = -refRangeRate / model.wavelength;

for iSat = 1:numSat
  satLos = model.satPosEci(:, iSat) - userPosEci;
  satRange = norm(satLos);
  if satRange <= 0
    error('buildDoaDopplerJacobian:InvalidSatRange', ...
      'Satellite %d coincides with the candidate user position.', iSat);
  end

  satLosUnit = satLos / satRange;
  satRelVel = model.satVelEci(:, iSat) - userVelEci;
  satRangeRate = satLosUnit' * satRelVel;
  satFd = -satRangeRate / model.wavelength;

  deltaFd(iSat) = satFd - refFd;
end

end

function diffVal = localAngleDiff(angleLeft, angleRight)
%LOCALANGLEDIFF Wrapped difference angleLeft - angleRight.

diffVal = angle(exp(1j * (angleLeft - angleRight)));

end
