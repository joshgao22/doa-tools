function [snaps, sampleCov, sourceSigs] = snapGenStochastic(array, doas, ...
  wavelength, numSnaps, noiseCov, sourceCov, useNormalizedAngles, modType)
%SNAPGENSTOCHASTIC Generate array snapshots with configurable source model.
% Supports single array struct and cell array of array structs.
%
%Syntax:
%   [snaps, sampleCov, sourceSigs] = snapGenStochastic(array, doas, wavelength, type)
%   ... = snapGenStochastic(array, doas, wavelength, type, numSnaps, noiseCov, sourceCov)
%   ... = snapGenStochastic(..., useNormalizedAngles, modType)
%
%Inputs:
%   array  - Array struct or cell array of structs. Required fields:
%            - positions : [dim, M]
%            - count     : M
%            - dim       : 1 / 2 / 3
%            Optional error fields supported by steeringMatrix():
%            - positionErrors, gainErrors, phaseErrors, spacing
%   doas   - DOAs (radians):
%            - 1D: 1×D or D×1 broadside angles
%            - 2D: 2×D az/el
%   wavelength - Wavelength (scalar, meters)
%   numSnaps   - Number of snapshots (default: 1)
%   noiseCov   - Noise covariance spec: scalar σ², vector (M×1), or matrix (M×M)
%   sourceCov  - Source covariance spec: scalar σ², vector (D×1), or matrix (D×D)
%   useNormalizedAngles - (Optional) pass-through to steeringMatrix (default: false)
%   modType    - (Optional, type="symbol") 'bpsk','qpsk','8psk','16qam','64qam',...
%                (default: "qpsk")
%
%Outputs:
%   snaps     - Snapshot matrix (M×T) or cell array for multiple arrays
%   sampleCov - Sample covariance (M×M) or cell array
%   sourceSigs- Source signals (D×T)

arguments
  array 
  doas {mustBeNumeric}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  numSnaps (1,1) {mustBePositive, mustBeInteger} = 1
  noiseCov {mustBeNumeric} = 1
  sourceCov {mustBeNumeric} = 1
  useNormalizedAngles (1,1) logical = false
  modType (1,:) char = 'none'
end

modType = lower(modType);

% -------------------------------------------------------------------------
% Determine number of sources
% -------------------------------------------------------------------------
if isrow(doas)
  numSources = numel(doas);
else
  if size(doas, 1) ~= 2
    error('snapGenStochastic:InvalidDoas', ...
      'For 2D DOA, doas must be 2×D.');
  end
  numSources = size(doas, 2);
end

% -------------------------------------------------------------------------
% Generate source signals (D×T)
% -------------------------------------------------------------------------
switch modType
  case 'none'
    sourceSigs = genComplexGaussianSamples(numSources, numSnaps, sourceCov);

  otherwise
    % genModulatedSymbols() should map sourceCov to symbol power/covariance.
    sourceSigs = genModulatedSymbols(numSources, numSnaps, modType, sourceCov);
end

% -------------------------------------------------------------------------
% Single vs multiple arrays
% -------------------------------------------------------------------------
if iscell(array)
  numArrays = numel(array);
  snaps = cell(1, numArrays);
  sampleCov = cell(1, numArrays);

  for k = 1:numArrays
    arr = array{k};

    A = steeringMatrix(arr, wavelength, doas, useNormalizedAngles); % M×D
    m = arr.count;

    signalPart = A * sourceSigs;                                  % M×T
    noisePart  = genComplexGaussianSamples(m, numSnaps, noiseCov); % M×T

    snaps{k} = signalPart + noisePart;
    sampleCov{k} = (snaps{k} * snaps{k}') / numSnaps;
  end

else
  A = steeringMatrix(array, wavelength, doas, useNormalizedAngles);
  m = array.count;

  signalPart = A * sourceSigs;
  noisePart  = genComplexGaussianSamples(m, numSnaps, noiseCov);

  snaps = signalPart + noisePart;
  sampleCov = (snaps * snaps') / numSnaps;
end

end


function samples = genComplexGaussianSamples(numRows, numCols, covSpec)
%GENCOMPLEXGAUSSIANSAMPLES Generate complex circularly symmetric Gaussian samples.
%
%Syntax:
%   samples = genComplexGaussianSamples(numRows, numCols, covSpec)
%
%Inputs:
%   numRows  - D
%   numCols  - N
%   covSpec  - Covariance specification across rows:
%             - scalar σ²      : IID rows, var σ²
%             - vector (D×1)   : diagonal covariance
%             - matrix (D×D)   : full covariance
%
%Output:
%   samples  - D×N CCSG samples

arguments
  numRows (1,1) {mustBePositive, mustBeInteger}
  numCols (1,1) {mustBePositive, mustBeInteger}
  covSpec {mustBeNumeric} = 1
end

base = randn(numRows, numCols) + 1j * randn(numRows, numCols);

if isscalar(covSpec)
  samples = sqrt(covSpec/2) * base;

elseif isvector(covSpec)
  v = covSpec(:);
  if numel(v) ~= numRows
    error('genComplexGaussianSamples:SizeMismatch', ...
      'Vector covariance must have length numRows.');
  end
  samples = bsxfun(@times, base, sqrt(v));

else
  if ~isequal(size(covSpec), [numRows, numRows])
    error('genComplexGaussianSamples:SizeMismatch', ...
      'Full covariance must be numRows×numRows.');
  end
  % covSpec is for complex vector; scale by 1/2 for real+imag parts
  T = sqrtm(covSpec/2);
  samples = T * base;
end

end


function symbols = genModulatedSymbols(numSym, numSnap, modType, sourceCov)
% sourceCov: scalar / (numSym×1) / (numSym×numSym)

arguments
  numSym (1,1) {mustBePositive, mustBeInteger}
  numSnap (1,1) {mustBePositive, mustBeInteger}
  modType (1,:) char {mustBeMember(modType, ['bpsk', 'qpsk', '8psk', '16qam', '64qam'])}
  sourceCov {mustBeNumeric} = 1
end

% --- generate IID unit-power symbols U ---
modType = lower(modType);
switch modType
  case 'bpsk'
    M = 2; u = pskmod(randi([0 M-1], numSym, numSnap), M);
  case 'qpsk'
    M = 4; u = pskmod(randi([0 M-1], numSym, numSnap), M);
  case '8psk'
    M = 8; u = pskmod(randi([0 M-1], numSym, numSnap), M);
  case '16qam'
    M = 16; u = qammod(randi([0 M-1], numSym, numSnap), M, 'UnitAveragePower', true);
  case '64qam'
    M = 64; u = qammod(randi([0 M-1], numSym, numSnap), M, 'UnitAveragePower', true);
  otherwise
    error('genModulatedSymbols:InvalidInput', 'Unsupported modulation type: %s', modType);
end

% --- map covariance: S = L U so that E[SS^H] = sourceCov ---
if isscalar(sourceCov)
  symbols = sqrt(sourceCov) * u;

elseif isvector(sourceCov)
  v = sourceCov(:);
  if numel(v) ~= numSym
    error('genModulatedSymbols:SizeMismatch', ...
      'Vector sourceCov must have length numSym.');
  end
  symbols = u .* sqrt(v);

else
  if ~isequal(size(sourceCov), [numSym, numSym])
    error('genModulatedSymbols:SizeMismatch', ...
      'Full sourceCov must be numSym×numSym.');
  end
  % Ensure Hermitian PSD assumed; chol may fail if not PSD
  L = chol((sourceCov + sourceCov')/2, 'lower');
  symbols = L * u;
end
end
