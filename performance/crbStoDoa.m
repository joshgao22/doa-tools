function [crb, aux] = crbStoDoa(array, wavelength, localDoa, pwrSource, ...
  noiseVar, numSnaps, crbOpt)
%CRBSTODOA Stochastic-model CRB for single-array or multi-array DOA parameters.
%
%Syntax:
%  crb = crbStoDoa(array, wavelength, localDoa, pwrSource, noiseVar)
%  crb = crbStoDoa(array, wavelength, localDoa, pwrSource, noiseVar, numSnaps)
%  [crb, aux] = crbStoDoa(array, wavelength, localDoa, pwrSource, ...
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
%  pwrSource  - Stochastic source covariance matrix / matrices
%               single-array:
%                 * scalar : pwrSource * eye(K)
%                 * vector : diag(pwrSource)
%                 * matrix : KxK covariance matrix
%               multi-array:
%                 * scalar / vector / matrix applied to all arrays
%                 * 1xL cell, one covariance entry per array
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
%  For array i, the stochastic model is
%
%    x_i(t) = A_i(localDoa_i) * s_i(t) + n_i(t),
%
%  with source covariance
%
%    P_i = E[s_i(t) s_i(t)'] .
%
%  In multi-array mode, the total FIM is the sum of array-wise FIMs. If
%  localDoa_i is generated from a shared structural parameter vector theta,
%  the chain rule is applied through crbOpt.localDoaJac.
%
%Reference:
%  P. Stoica and A. Nehorai, "Performance study of conditional and
%  unconditional direction-of-arrival estimation," IEEE Transactions on
%  Acoustics, Speech and Signal Processing, vol. 38, no. 10,
%  pp. 1783-1795, Oct. 1990.
%
%See also:
%  steeringMatrix, crbDetDoa, crbPilotDoaDoppler

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
        error('crbStoDoa:MissingDoaJacobian', ...
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
      error('crbStoDoa:ParamNameSizeMismatch', ...
        'crbOpt.paramName length must match the number of structural parameters.');
    end
  end
end

fim = zeros(numParam, numParam);
fimCell = cell(1, numArray);
localFimCell = cell(1, numArray);

for iArray = 1:numArray
  [A, dA] = steeringMatrix(arrayCell{iArray}, wavelength, localDoaCell{iArray});
  numSensor = size(A, 1);

  gramA = A' * A;
  gramA = 0.5 * (gramA + gramA');
  if rcond(gramA) < 1e-12
    warning('crbStoDoa:IllConditionedSteering', ...
      'A'' * A is singular or ill-conditioned at array %d. CRB may be unreliable.', iArray);
  end

  projOrth = eye(numSensor) - A * (gramA \ A');
  projOrth = 0.5 * (projOrth + projOrth');

  covSig = A * pwrSourceCell{iArray} * A';
  covMat = covSig + noiseVar(iArray) * eye(numSensor);
  covMat = 0.5 * (covMat + covMat');

  if rcond(covMat) < 1e-12
    warning('crbStoDoa:IllConditionedCovariance', ...
      'Array covariance matrix is singular or ill-conditioned at array %d. CRB may be unreliable.', iArray);
  end

  qMat = pwrSourceCell{iArray} * (A' * (covMat \ A)) * pwrSourceCell{iArray};
  qMat = 0.5 * (qMat + qMat');

  localFimCore = localBuildLocalFimCore(dA, projOrth, qMat);
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
  warning('crbStoDoa:IllConditionedFim', ...
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
  error('crbStoDoa:EmptyArrayInput', ...
    'array must not be empty.');
end

for iArray = 1:numel(arrayCell)
  currentArray = arrayCell{iArray};
  if ~isstruct(currentArray) || ~isfield(currentArray, 'positions') || isempty(currentArray.positions)
    error('crbStoDoa:InvalidArrayEntry', ...
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
    error('crbStoDoa:DoaCountMismatch', ...
      'The number of localDoa entries must match the number of arrays.');
  end
  localDoaCell = cellfun(@localNormalizeDoa, localDoaCell, 'UniformOutput', false);
end

refDoa = localDoaCell{1};
[paramMode, numSource, numLocalParam] = localGetDoaMeta(refDoa);

for iArray = 2:numArray
  [paramModeCur, numSourceCur, numLocalParamCur] = localGetDoaMeta(localDoaCell{iArray});
  if ~strcmp(paramModeCur, paramMode) || numSourceCur ~= numSource || numLocalParamCur ~= numLocalParam
    error('crbStoDoa:DoaShapeMismatch', ...
      'All localDoa entries must have consistent DOA mode and source count.');
  end
end

end

function doa = localNormalizeDoa(doa)
%LOCALNORMALIZEDOA Normalize one DOA entry.

if ~isnumeric(doa) || isempty(doa)
  error('crbStoDoa:InvalidDoaEntry', ...
    'Each localDoa entry must be a non-empty numeric array.');
end

if isvector(doa)
  doa = reshape(doa, 1, []);
elseif size(doa, 1) == 2
  % keep 2xK
else
  error('crbStoDoa:InvalidDoaSize', ...
    'localDoa must be a vector or a 2xK matrix.');
end

end

function [paramMode, numSource, numLocalParam] = localGetDoaMeta(doa)
%LOCALGETDOAMETA Get local DOA mode and parameter count.

if isvector(doa) || size(doa, 1) == 1
  paramMode = 'theta';
  numSource = numel(doa);
  numLocalParam = numSource;
else
  paramMode = 'azel';
  numSource = size(doa, 2);
  numLocalParam = 2 * numSource;
end

if numSource < 1
  error('crbStoDoa:InvalidSourceCount', ...
    'At least one source is required.');
end

end

function pwrSourceCell = localParseSourcePower(pwrSource, numArray, numSource)
%LOCALPARSESOURCEPOWER Normalize source covariance input to a 1xL cell array.

if ~iscell(pwrSource)
  pwrSourceCell = repmat({localUnifySourcePower(pwrSource, numSource)}, 1, numArray);
else
  pwrSourceCell = reshape(pwrSource, 1, []);
  if numel(pwrSourceCell) ~= numArray
    error('crbStoDoa:SourcePowerCountMismatch', ...
      'The number of pwrSource entries must match the number of arrays.');
  end
  for iArray = 1:numArray
    pwrSourceCell{iArray} = localUnifySourcePower(pwrSourceCell{iArray}, numSource);
  end
end

for iArray = 1:numArray
  pwrSourceCell{iArray} = 0.5 * (pwrSourceCell{iArray} + pwrSourceCell{iArray}');
end

end

function pwrSource = localUnifySourcePower(pwrSource, numSource)
%LOCALUNIFYSOURCEPOWER Convert scalar / vector / matrix source power input.

if ~isnumeric(pwrSource) || isempty(pwrSource)
  error('crbStoDoa:InvalidSourcePower', ...
    'pwrSource must be numeric and non-empty.');
end

if isscalar(pwrSource)
  pwrSource = pwrSource * eye(numSource);

elseif isvector(pwrSource)
  if numel(pwrSource) ~= numSource
    error('crbStoDoa:SourcePowerSizeMismatch', ...
      'Source power vector length must equal the number of sources.');
  end
  pwrSource = diag(pwrSource(:));

else
  if ~isequal(size(pwrSource), [numSource, numSource])
    error('crbStoDoa:SourcePowerSizeMismatch', ...
      'Source power matrix must be a %d x %d matrix.', numSource, numSource);
  end
end

end

function value = localExpandRealPositive(value, numArray, valueName)
%LOCALEXPANDREALPOSITIVE Expand scalar or vector input to array-wise vector.

if ~isnumeric(value) || isempty(value)
  error('crbStoDoa:InvalidNumericInput', ...
    '%s must be numeric and non-empty.', valueName);
end

if isscalar(value)
  value = repmat(double(value), numArray, 1);
else
  value = double(value(:));
  if numel(value) ~= numArray
    error('crbStoDoa:InputLengthMismatch', ...
      '%s must be scalar or have one entry per array.', valueName);
  end
end

if any(~isfinite(value)) || any(value <= 0)
  error('crbStoDoa:NonPositiveInput', ...
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
    error('crbStoDoa:DoaJacCountMismatch', ...
      'crbOpt.localDoaJac must have one entry per array.');
  end

elseif isnumeric(localDoaJac)
  if ndims(localDoaJac) == 2
    localDoaJacCell = repmat({localDoaJac}, 1, numArray);
  elseif ndims(localDoaJac) == 3
    if size(localDoaJac, 3) ~= numArray
      error('crbStoDoa:DoaJacCountMismatch', ...
        'The third dimension of crbOpt.localDoaJac must match the number of arrays.');
    end
    localDoaJacCell = cell(1, numArray);
    for iArray = 1:numArray
      localDoaJacCell{iArray} = localDoaJac(:, :, iArray);
    end
  else
    error('crbStoDoa:InvalidDoaJacType', ...
      'crbOpt.localDoaJac must be a cell, a 2-D matrix, or a 3-D array.');
  end

else
  error('crbStoDoa:InvalidDoaJacType', ...
    'crbOpt.localDoaJac must be numeric or cell.');
end

numParam = size(localDoaJacCell{1}, 2);
for iArray = 1:numArray
  currentJac = localDoaJacCell{iArray};
  if ~isnumeric(currentJac) || size(currentJac, 1) ~= numLocalParam || size(currentJac, 2) ~= numParam
    error('crbStoDoa:DoaJacSizeMismatch', ...
      'Each localDoaJac entry must be of size %d x P.', numLocalParam);
  end
end

end

function fimCore = localBuildLocalFimCore(dA, projOrth, qMat)
%LOCALBUILDLOCALFIMCORE Build local stochastic-model FIM core.

if isfield(dA, 'theta')
  gradTheta = dA.theta;
  fimCore = real((gradTheta' * projOrth * gradTheta) .* qMat.');

elseif isfield(dA, 'az') && isfield(dA, 'el')
  gradAz = dA.az;
  gradEl = dA.el;

  coreAzAz = real((gradAz' * projOrth * gradAz) .* qMat.');
  coreAzEl = real((gradAz' * projOrth * gradEl) .* qMat.');
  coreElAz = real((gradEl' * projOrth * gradAz) .* qMat.');
  coreElEl = real((gradEl' * projOrth * gradEl) .* qMat.');

  fimCore = [coreAzAz, coreAzEl;
    coreElAz, coreElEl];

else
  error('crbStoDoa:UnsupportedDerivative', ...
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
    error('crbStoDoa:InvalidParamMode', ...
      'Unsupported parameter mode: %s.', paramMode);
end

end
