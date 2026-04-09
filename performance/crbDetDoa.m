function [crb, aux] = crbDetDoa(array, wavelength, localDoa, pwrSource, ...
  noiseVar, numSnaps, crbOpt)
%CRBDETDOA Deterministic-model CRB for single-array or multi-array DOA parameters.
%
%Syntax:
%  crb = crbDetDoa(array, wavelength, localDoa, pwrSource, noiseVar)
%  crb = crbDetDoa(array, wavelength, localDoa, pwrSource, noiseVar, numSnaps)
%  [crb, aux] = crbDetDoa(array, wavelength, localDoa, pwrSource, ...
%    noiseVar, numSnaps, crbOpt)
%
%Inputs:
%  array      - Array structure, or 1xL cell array of array structures
%
%  wavelength - Positive scalar wavelength
%
%  localDoa   - Local DOA used to build steering matrices
%               single-array:
%                 * 1xK / Kx1 vector
%                 * 2xK matrix [az; el]
%               multi-array:
%                 * numeric : same local DOA applied to all arrays
%                 * cell    : 1xL cell, one DOA entry per array
%
%  pwrSource  - Deterministic source power matrix / matrices
%               single-array:
%                 * KxK matrix
%               multi-array:
%                 * KxK matrix applied to all arrays
%                 * 1xL cell, one KxK matrix per array
%
%  noiseVar   - Noise power
%               * scalar : same for all arrays
%               * vector : one value per array
%
%  numSnaps   - Snapshot count
%               * scalar : same for all arrays
%               * vector : one value per array
%               default is 1
%
%  crbOpt     - Optional settings structure
%               .localDoaJac : Jacobian(s) from shared structural parameter
%                              vector theta to vec(localDoa_i)
%                              empty : estimate localDoa itself
%                              numeric / cell supported
%
%                              1D localDoa ordering:
%                                [theta1, ..., thetaK].'
%
%                              2D localDoa ordering:
%                                [az1, ..., azK, el1, ..., elK].'
%
%                              Size of each Jacobian:
%                                numLocalParam x numStructParam
%
%               .paramName    : optional parameter names stored in aux
%
%Outputs:
%  crb        - CRB matrix
%               * if crbOpt.localDoaJac is empty:
%                   CRB of the local DOA parameter vector
%               * otherwise:
%                   CRB of the shared structural parameter vector
%
%  aux        - Auxiliary structure with fields:
%               .fim
%               .fimCell
%               .localFimCell
%               .localDoaCell
%               .localDoaJacCell
%               .pwrSourceCell
%               .noiseVar
%               .numSnaps
%               .numArray
%               .numSource
%               .paramName
%               .paramMode
%
%Description:
%  For array i, the deterministic model is
%
%    X_i = A_i(localDoa_i) * S_i + N_i ,
%
%  where S_i is an unknown deterministic source signal matrix and
%
%    P_i = (S_i * S_i') / numSnaps_i
%
%  is supplied through pwrSource. In multi-array mode, the total FIM is the
%  sum of array-wise FIMs. If localDoa_i is generated from a shared
%  structural parameter vector theta, the chain rule is applied through
%  crbOpt.localDoaJac.
%
%Notes:
%  - Path gain is not a separate input here. For deterministic CRB, it can
%    be absorbed into pwrSource:
%       single source : P_i = abs(alpha_i)^2 * P0
%       multi source  : P_i = diag(alpha_i) * P0 * diag(alpha_i)'
%  - When array is a cell and localDoa is also a cell, but localDoaJac is
%    empty, all localDoa entries must be identical. Otherwise there is no
%    common parameter vector to jointly accumulate the FIM.
%
%See also:
%  steeringMatrix, crbPilotDoaDoppler, buildDoaDopplerJacobian

arguments
  array (1,:) {mustBeA(array, ["struct", "cell"])}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  localDoa {mustBeA(localDoa, ["double", "cell"])}
  pwrSource {mustBeA(pwrSource, ["double", "cell"])}
  noiseVar {mustBeNumeric}
  numSnaps {mustBeNumeric} = 1
  crbOpt (1,1) struct = struct()
end

crbOpt = localParseCrbOpt(crbOpt);

arrayCell = localParseArrayCell(array);
numArray = numel(arrayCell);

[localDoaCell, numSource, paramMode, numLocalParam] = ...
  localParseDoaCell(localDoa, numArray);

pwrSourceCell = localParseSourcePower(pwrSource, numArray, numSource);
noiseVar = localExpandRealPositive(noiseVar, numArray, 'noiseVar');
numSnaps = localExpandRealPositive(numSnaps, numArray, 'numSnaps');

localDoaJacCell = localParseDoaJac(crbOpt.localDoaJac, numArray, numLocalParam);

if isempty(localDoaJacCell)
  if iscell(localDoa)
    refDoa = localDoaCell{1};
    for iArray = 2:numArray
      if ~isequal(size(localDoaCell{iArray}), size(refDoa)) || ...
          any(abs(localDoaCell{iArray}(:) - refDoa(:)) > 0)
        error('crbDetDoa:MissingDoaJacobian', ...
          ['When localDoa is cell-valued and differs across arrays, ' ...
           'crbOpt.localDoaJac must be provided.']);
      end
    end
  end
  numParam = numLocalParam;
  paramName = localDefaultParamName(paramMode, numSource);
else
  numParam = size(localDoaJacCell{1}, 2);
  if isempty(crbOpt.paramName)
    paramName = "param" + string(1:numParam);
  else
    paramName = reshape(string(crbOpt.paramName), 1, []);
    if numel(paramName) ~= numParam
      error('crbDetDoa:ParamNameSizeMismatch', ...
        'crbOpt.paramName length must match the number of structural parameters.');
    end
  end
end

fim = zeros(numParam, numParam);
fimCell = cell(1, numArray);
localFimCell = cell(1, numArray);

for iArray = 1:numArray
  [~, dA] = steeringMatrix(arrayCell{iArray}, wavelength, localDoaCell{iArray});
  A = steeringMatrix(arrayCell{iArray}, wavelength, localDoaCell{iArray});

  gramA = A' * A;
  gramA = 0.5 * (gramA + gramA');
  if rcond(gramA) < 1e-12
    warning('crbDetDoa:IllConditionedSteering', ...
      'A'' * A is singular or ill-conditioned at array %d. CRB may be unreliable.', iArray);
  end

  numSensor = size(A, 1);
  projOrth = eye(numSensor) - A * (gramA \ A');
  projOrth = 0.5 * (projOrth + projOrth');

  localFimCore = localBuildLocalFimCore(dA, projOrth, pwrSourceCell{iArray});
  localFimCore = 0.5 * (localFimCore + localFimCore.');

  localFimCell{iArray} = localFimCore;

  if isempty(localDoaJacCell)
    currentFim = (2 * numSnaps(iArray) / noiseVar(iArray)) * localFimCore;
  else
    currentJac = localDoaJacCell{iArray};
    currentFim = (2 * numSnaps(iArray) / noiseVar(iArray)) * ...
      (currentJac' * localFimCore * currentJac);
  end

  currentFim = 0.5 * (currentFim + currentFim.');
  fimCell{iArray} = currentFim;
  fim = fim + currentFim;
end

fim = 0.5 * (fim + fim.');
if rcond(fim) < 1e-12
  warning('crbDetDoa:IllConditionedFim', ...
    'FIM is singular or ill-conditioned. CRB may be unreliable.');
end

crb = fim \ eye(size(fim, 1));
crb = 0.5 * (crb + crb.');

if nargout >= 2
  aux = struct();
  aux.fim = fim;
  aux.fimCell = fimCell;
  aux.localFimCell = localFimCell;
  aux.localDoaCell = localDoaCell;
  aux.localDoaJacCell = localDoaJacCell;
  aux.pwrSourceCell = pwrSourceCell;
  aux.noiseVar = noiseVar;
  aux.numSnaps = numSnaps;
  aux.numArray = numArray;
  aux.numSource = numSource;
  aux.paramName = cellstr(paramName);
  aux.paramMode = paramMode;
else
  aux = [];
end

end

function crbOpt = localParseCrbOpt(crbOpt)
%LOCALPARSECRBOPT Parse optional CRB settings.

if ~isfield(crbOpt, 'localDoaJac') || isempty(crbOpt.localDoaJac)
  crbOpt.localDoaJac = [];
end
if ~isfield(crbOpt, 'paramName') || isempty(crbOpt.paramName)
  crbOpt.paramName = [];
end

end

function arrayCell = localParseArrayCell(array)
%LOCALPARSEARRAYCELL Normalize array input to a 1xL cell array.

if ~iscell(array)
  arrayCell = {array};
else
  arrayCell = reshape(array, 1, []);
end

if isempty(arrayCell)
  error('crbDetDoa:EmptyArrayInput', ...
    'array must not be empty.');
end

for iArray = 1:numel(arrayCell)
  currentArray = arrayCell{iArray};
  if ~isstruct(currentArray) || ~isfield(currentArray, 'positions') || isempty(currentArray.positions)
    error('crbDetDoa:InvalidArrayEntry', ...
      'Each array entry must be a struct with a non-empty positions field.');
  end
end

end

function [localDoaCell, numSource, paramMode, numLocalParam] = localParseDoaCell(localDoa, numArray)
%LOCALPARSEDOACELL Normalize localDoa input to a 1xL cell array.

if ~iscell(localDoa)
  localDoaCell = repmat({localNormalizeDoa(localDoa)}, 1, numArray);
else
  localDoaCell = reshape(localDoa, 1, []);
  if numel(localDoaCell) ~= numArray
    error('crbDetDoa:DoaCountMismatch', ...
      'The number of localDoa entries must match the number of arrays.');
  end
  localDoaCell = cellfun(@localNormalizeDoa, localDoaCell, 'UniformOutput', false);
end

refDoa = localDoaCell{1};
[paramMode, numSource, numLocalParam] = localGetDoaMeta(refDoa);

for iArray = 2:numArray
  [paramModeCur, numSourceCur, numLocalParamCur] = localGetDoaMeta(localDoaCell{iArray});
  if ~strcmp(paramModeCur, paramMode) || numSourceCur ~= numSource || numLocalParamCur ~= numLocalParam
    error('crbDetDoa:DoaShapeMismatch', ...
      'All localDoa entries must have consistent DOA mode and source count.');
  end
end

end

function doa = localNormalizeDoa(doa)
%LOCALNORMALIZEDOA Normalize one DOA entry.

if ~isnumeric(doa) || isempty(doa)
  error('crbDetDoa:InvalidDoaEntry', ...
    'Each localDoa entry must be a non-empty numeric array.');
end

if size(doa, 1) == 2
  doa = reshape(doa, 2, []);
elseif isvector(doa)
  doa = reshape(doa, 1, []);
else
  error('crbDetDoa:InvalidDoaSize', ...
    'localDoa must be a vector or a 2xK matrix.');
end

end

function [paramMode, numSource, numLocalParam] = localGetDoaMeta(doa)
%LOCALGETDOAMETA Get local DOA mode and parameter count.

if size(doa, 1) == 2
  paramMode = 'azel';
  numSource = size(doa, 2);
  numLocalParam = 2 * numSource;
elseif size(doa, 1) == 1
  paramMode = 'theta';
  numSource = numel(doa);
  numLocalParam = numSource;
else
  error('%s:InvalidDoaSizeInternal', ...
    'Normalized localDoa must be 1xK or 2xK.', mfilename);
end

if numSource < 1
  error('crbDetDoa:InvalidSourceCount', ...
    'At least one source is required.');
end

end

function pwrSourceCell = localParseSourcePower(pwrSource, numArray, numSource)
%LOCALPARSESOURCEPOWER Normalize source power input to a 1xL cell array.

if ~iscell(pwrSource)
  pwrSourceCell = repmat({pwrSource}, 1, numArray);
else
  pwrSourceCell = reshape(pwrSource, 1, []);
  if numel(pwrSourceCell) ~= numArray
    error('crbDetDoa:SourcePowerCountMismatch', ...
      'The number of pwrSource entries must match the number of arrays.');
  end
end

for iArray = 1:numArray
  currentPwr = pwrSourceCell{iArray};
  if ~isnumeric(currentPwr) || ~isequal(size(currentPwr), [numSource, numSource])
    error('crbDetDoa:SourcePowerSizeMismatch', ...
      'Each pwrSource entry must be a %d x %d matrix.', numSource, numSource);
  end
  pwrSourceCell{iArray} = 0.5 * (currentPwr + currentPwr');
end

end

function value = localExpandRealPositive(value, numArray, valueName)
%LOCALEXPANDREALPOSITIVE Expand scalar or vector input to array-wise vector.

if ~isnumeric(value) || isempty(value)
  error('crbDetDoa:InvalidNumericInput', ...
    '%s must be numeric and non-empty.', valueName);
end

if isscalar(value)
  value = repmat(double(value), numArray, 1);
else
  value = double(value(:));
  if numel(value) ~= numArray
    error('crbDetDoa:InputLengthMismatch', ...
      '%s must be scalar or have one entry per array.', valueName);
  end
end

if any(~isfinite(value)) || any(value <= 0)
  error('crbDetDoa:NonPositiveInput', ...
    '%s must contain positive finite values.', valueName);
end

end

function localDoaJacCell = localParseDoaJac(localDoaJac, numArray, numLocalParam)
%LOCALPARSEDOAJAC Normalize local DOA Jacobians.

if isempty(localDoaJac)
  localDoaJacCell = {};
  return;
end

if iscell(localDoaJac)
  localDoaJacCell = reshape(localDoaJac, 1, []);
  if numel(localDoaJacCell) ~= numArray
    error('crbDetDoa:DoaJacCountMismatch', ...
      'crbOpt.localDoaJac must have one entry per array.');
  end

elseif isnumeric(localDoaJac)
  if ndims(localDoaJac) == 2
    localDoaJacCell = repmat({localDoaJac}, 1, numArray);
  elseif ndims(localDoaJac) == 3
    if size(localDoaJac, 3) ~= numArray
      error('crbDetDoa:DoaJacCountMismatch', ...
        'The third dimension of crbOpt.localDoaJac must match the number of arrays.');
    end
    localDoaJacCell = cell(1, numArray);
    for iArray = 1:numArray
      localDoaJacCell{iArray} = localDoaJac(:, :, iArray);
    end
  else
    error('crbDetDoa:InvalidDoaJacType', ...
      'crbOpt.localDoaJac must be a cell, a 2-D matrix, or a 3-D array.');
  end

else
  error('crbDetDoa:InvalidDoaJacType', ...
    'crbOpt.localDoaJac must be numeric or cell.');
end

numParam = size(localDoaJacCell{1}, 2);
for iArray = 1:numArray
  currentJac = localDoaJacCell{iArray};
  if ~isnumeric(currentJac) || size(currentJac, 1) ~= numLocalParam || size(currentJac, 2) ~= numParam
    error('crbDetDoa:DoaJacSizeMismatch', ...
      'Each localDoaJac entry must be of size %d x P.', numLocalParam);
  end
end

end

function fimCore = localBuildLocalFimCore(dA, projOrth, pwrSource)
%LOCALBUILDLOCALFIMCORE Build local deterministic-model FIM core.

if isfield(dA, 'theta')
  gradTheta = dA.theta;
  fimCore = real((gradTheta' * projOrth * gradTheta) .* pwrSource.');

elseif isfield(dA, 'az') && isfield(dA, 'el')
  gradAz = dA.az;
  gradEl = dA.el;

  coreAzAz = real((gradAz' * projOrth * gradAz) .* pwrSource.');
  coreAzEl = real((gradAz' * projOrth * gradEl) .* pwrSource.');
  coreElAz = real((gradEl' * projOrth * gradAz) .* pwrSource.');
  coreElEl = real((gradEl' * projOrth * gradEl) .* pwrSource.');

  fimCore = [coreAzAz, coreAzEl;
    coreElAz, coreElEl];

else
  error('crbDetDoa:UnsupportedDerivative', ...
    'Unsupported derivative fields returned by steeringMatrix.');
end

fimCore = 0.5 * (fimCore + fimCore.');

end

function paramName = localDefaultParamName(paramMode, numSource)
%LOCALDEFAULTPARAMNAME Build default parameter names.

switch paramMode
  case 'theta'
    paramName = "theta" + string(1:numSource);

  case 'azel'
    paramName = [ ...
      "az" + string(1:numSource), ...
      "el" + string(1:numSource)];

  otherwise
    error('crbDetDoa:InvalidParamMode', ...
      'Unsupported parameter mode: %s.', paramMode);
end

end
