function [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, ...
  wavelen, sampleCovMat, numSource, doaGrid, initParam, verbose)
%ESTIMATORMLESTOCONOPT Stochastic ML (uncorrelated) estimator via joint optimization.
% Solves the stochastic uncorrelated-source ML problem by directly minimizing
%   min_{Theta, p, sigma} logdet(S) + tr(S^{-1} R)
%   s.t. S = A(Theta) diag(p) A(Theta)^H + sigma I,
%        p >= 0, sigma >= 0,
% where Theta contains the DOA parameters per source (1D or 2D, read from doaGrid.dimension).
% The result is a continuous-domain estimate (non-grid/discrete).
%
%Syntax:
%  [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, wavelength, sampleCovMat, ...
%    numSource, doaGrid)
%  ... = estimatorMLEStoConOpt(..., verbose)
%
%Inputs:
%  array        - Array struct for steeringMatrix().
%  wavelen      - Wavelength (scalar).
%  sampleCovMat - Sample covariance (MxM), Hermitian.
%  numSource    - Number of sources (K).
%  doaGrid      - DOA grid struct (used for bounds and initial guess).
%                - dimension: 1 or 2
%                - type: 'local' recommended
%                - angleGrid: dim x N (radians), optional for MVDR initialization
%  verbose      - (Optional) true/false, show fmincon iterations if true (default: false).
%
%Outputs:
%  doaResult    - DOA result struct:
%                - doaGrid, spectrum (empty), peakIndices (empty)
%                - estimatedDoas: dim x K DOAs (radians)
%                - isResolved: true if optimization succeeds
%                - isDiscrete: true (discrete solution)
%  signalPower  - Estimated source powers (Kx1), nonnegative.
%  noiseVar     - Estimated noise variance (scalar), nonnegative.

arguments (Input)
  array (1,:) {mustBeA(array, ["struct", "cell"])}
  wavelen (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, ["double", "cell"])}
  numSource (1,1) {mustBePositive, mustBeInteger}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  initParam (1,:) = []
  verbose (1,1) logical = false
end

isArrayCell = iscell(array);
numArray = numel(array);

if ~isArrayCell
  array = {array};
  sampleCovMat = {sampleCovMat};
  doaGrid = {doaGrid};
end

if numel(sampleCovMat) ~= numArray || numel(doaGrid) ~= numArray
  error('estimatorMusic:InvalidInput', ...
    'Cell inputs "array", "sampleCovMat", and "doaGrid" must have the same number of elements.');
end

baseGrid = doaGrid{1};
dim = baseGrid.dimension;
numDoaVars = dim * numSource;
numOptVars = numDoaVars + numSource*numArray + numArray;
for ii = 1:numArray
  sampleCovMat{ii} = 0.5 * (sampleCovMat{ii} + sampleCovMat{ii}');  % Hermitianize
end

if isempty(initParam)
  initParam = zeros(numOptVars, 1);

  % -------------------------------------------------------------------------
  % Initial guess from MVDR (dim x K, radians)
  % -------------------------------------------------------------------------
  mvdrResult = estimatorMvdr(array, wavelen, sampleCovMat, numSource, doaGrid);

  % choose init DOA parameterization based on gridType
  switch baseGrid.type
    case 'angle'
      if mvdrResult.isResolved
        initParam(1:numDoaVars) = mvdrResult.angleEst(:); % radians, dim x K
      else
        warning('estimatorMleStoConOpt: Failed to obtain a good initial guess.');
        % we just assume the doas is uniformly placed
        switch dim
          case 1
            initParam(1:numDoaVars) = ...
              linspace(baseGrid.range(1,1), baseGrid.range(1,2), numSource);
          case 2
            initParam(1:2:numDoaVars-1) = ...
              linspace(baseGrid.range(1,1), baseGrid.range(1,2), numSource);
            initParam(2:2:numDoaVars) = ...
              linspace(baseGrid.range(2,1), baseGrid.range(2,2), numSource);
        end
      end
      Theta0 = reshape(initParam(1:numDoaVars), dim, numSource);
      A = steeringMatrix(array{1}, wavelen, Theta0);
      B = [khatri_rao(conj(A), A), reshape(eye(size(A, 1)), [], 1)];
      z = real(B \ sampleCovMat{1}(:));
      z(z < 0) = 0;
      initParam(numDoaVars + 1:end) = z;

    case 'latlon'
      if mvdrResult.isResolved
        initParam(1:dim*numSource) = mvdrResult.latlonEst(:); % degrees, 2 x K  (lat; lon)
      else
        warning('estimatorMleStoConOpt: Failed to obtain a good initial guess.');
        initParam(1:2:numDoaVars-1) = ...
          linspace(baseGrid.range(1,1), baseGrid.range(1,2), numSource);
        initParam(2:2:numDoaVars) = ...
          linspace(baseGrid.range(2,1), baseGrid.range(2,2), numSource);
      end

      for ii = 1:numArray
        % altitude fixed at 0
        [x,y,z] = geodetic2ecef(doaGrid{ii}.spheroid, initParam(1:2:numDoaVars-1), ...
          initParam(2:2:numDoaVars), zeros(numSource, 1));
        ecef = [x'; y'; z'];  % 3xK

        % local angles in radians, columns are [az; el]
        Theta0 = ecefToAngleGrid(ecef, doaGrid{ii}.arrayCenter); % 2xK
        A = steeringMatrix(array{ii}, wavelen, Theta0);
        B = [khatri_rao(conj(A), A), reshape(eye(size(A, 1)), [], 1)];
        z = real(B \ sampleCovMat{ii}(:));
        z(z < 0) = 0;

        idxP = numDoaVars + (ii-1)*numSource + (1:numSource);
        idxN = numDoaVars + numSource*numArray + ii;
        initParam(idxP) = z(1:end-1);
        initParam(idxN) = z(end);
      end
  end
end

% -------------------------------------------------------------------------
% Bounds
% -------------------------------------------------------------------------
lb = zeros(numOptVars, 1);
ub = inf(numOptVars, 1);

% DOA bounds
[lb, ub] = setDoaBounds(lb, ub, baseGrid, dim, numSource);

% p and noiseVar already [0, inf)

% -------------------------------------------------------------------------
% fmincon options
% -------------------------------------------------------------------------
if verbose
  options = optimoptions('fmincon', 'Display', 'iter');
else
  options = optimoptions('fmincon', 'Display', 'off');
end

% -------------------------------------------------------------------------
% Solve
% -------------------------------------------------------------------------
objFun = @(optVar) nllStoUc(array, wavelen, sampleCovMat, doaGrid, ...
  dim, numSource, optVar);

[optVar, ~, exitflag] = fmincon(objFun, initParam, [], [], [], [], lb, ub, [], options);

% -------------------------------------------------------------------------
% Unpack and post-process
% -------------------------------------------------------------------------
estDoaParam = reshape(optVar(1:numDoaVars), dim, numSource);  % angle(rad) or latlon(deg)
signalPower = optVar(numDoaVars + (1:numSource*numArray));
noiseVar    = optVar(end-numArray+1:end);

% Convert to angles for output consistency
switch baseGrid.type
  case 'angle'
    estAngle = estDoaParam; % radians
    estLatlon = [];
  case 'latlon'
    estLatlon = estDoaParam; % degrees
    if numArray == 1
      [x,y,z] = geodetic2ecef(baseGrid.spheroid, estLatlon(1,:), estLatlon(2,:), zeros(1,numSource));
      estEcef = [x;y;z];
      estAngle = ecefToAngleGrid(estEcef, baseGrid.arrayCenter); % radians (az;el)
      [~, sortIdx] = sort(estAngle(1,:));
      estAngle     = estAngle(:, sortIdx);
      signalPower  = signalPower(sortIdx);
    else
      estAngle = [];
    end
end

% Ensure nonnegativity
signalPower = max(signalPower, 0);
noiseVar = max(noiseVar, 0);

% -------------------------------------------------------------------------
% Build doaResult (MUSIC-like output struct)
% -------------------------------------------------------------------------
doaResult = struct();
doaResult.doaGrid = doaGrid;
doaResult.spectrum = [];
doaResult.peakIndices = [];
doaResult.angleEst = estAngle;
doaResult.latlonEst = estLatlon;
doaResult.isResolved = exitflag > 0 && all(isfinite(estAngle(:)));
doaResult.isDiscrete = true;

end


% ===== helpers =====

function [lb, ub] = setDoaBounds(lb, ub, doaGrid, dim, numSource)
%SETDOABOUNDS Set box constraints for continuous DOA optimization variables.
%
%This function fills the lower/upper bounds of the DOA-related entries in the
%optimization vector, according to the grid domain stored in doaGrid.
%
%DOA parameterization:
%  - 1D DOA: theta (broadside angle)
%  - 2D DOA: azimuth–elevation pairs [az; el]
%
%Bounds:
%  - 1D: theta ∈ [doaGrid.range(1), doaGrid.range(2)]
%  - 2D: az ∈ [doaGrid.range(1,1), doaGrid.range(1,2)]
%         el ∈ [doaGrid.range(2,1), doaGrid.range(2,2)]
%
%Optimization vector layout (MATLAB column-major stacking):
%  Let K = numSource, dim ∈ {1,2}.
%
%  - dim = 1:
%      Theta(:) = [theta_1; theta_2; ...; theta_K]
%
%  - dim = 2:
%      Theta(:) = [az_1; el_1; az_2; el_2; ...; az_K; el_K]
%
%Notes:
%  - This function assumes doaGrid.range has already been validated.
%  - Only the DOA-related entries are modified; bounds for power/noise
%    variables should be handled separately.

if dim == 1
  lb(1:numSource) = doaGrid.range(1);
  ub(1:numSource) = doaGrid.range(2);
  return;
end

azIdx = 1:2:(2*numSource-1);
elIdx = 2:2:(2*numSource);

switch doaGrid.type
  case 'angle'
    % range = [azMin azMax; elMin elMax] (rad)
    lb(azIdx) = doaGrid.range(1,1);  ub(azIdx) = doaGrid.range(1,2);
    lb(elIdx) = doaGrid.range(2,1);  ub(elIdx) = doaGrid.range(2,2);

  case 'latlon'
    % range = [latMin latMax; lonMin lonMax] (deg)
    lb(azIdx) = doaGrid.range(1,1);  ub(azIdx) = doaGrid.range(1,2); % lat
    lb(elIdx) = doaGrid.range(2,1);  ub(elIdx) = doaGrid.range(2,2); % lon

  otherwise
    error('setDoaBounds:InvalidType', 'doaGrid.type must be ''angle'' or ''latlon''.');
end
end


function obj = nllStoUc(array, wavelength, sampleCovMat, doaGrid, dim, numSource, optVar)
%NLLSTOUC Stochastic uncorrelated-source negative log-likelihood.
%
%Computes the objective
%    f(Theta, p, sigma) = log det(S) + tr(S^{-1} Rhat),
%where
%    S = A(Theta) diag(p) A(Theta)^H + sigma I,
%    Rhat = sampleCovMat,
%    p >= 0 is the source power vector, and sigma >= 0 is the noise variance.
%
%DOA parameterization and optVar layout:
%  optVar = [doaParam(:); p; sigma]
%  - doaParam depends on doaGrid.type:
%      * 'angle'  : doaParam(:) is Theta(:) with Theta = [az; el] (rad) for dim=2
%                  (MATLAB column-major stacking: [az1; el1; az2; el2; ...]).
%      * 'latlon' : doaParam(:) stores [lat1; lon1; lat2; lon2; ...] (deg).
%                  Each (lat,lon) is mapped to local (az,el) (rad) using
%                  geodetic2ecef() + ecefToAngleGrid() w.r.t. doaGrid.arrayCenter.
%  - p is a Kx1 vector (K = numSource).
%  - sigma is a scalar.
%
%Inputs:
%  array        : array struct for steeringMatrix()
%  wavelength   : wavelength (scalar)
%  sampleCovMat : sample covariance Rhat (MxM, Hermitian)
%  doaGrid      : grid struct with fields:
%                 - type in {'angle','latlon'}
%                 - spheroid, arrayCenter (required for 'latlon')
%  dim          : DOA dimension (1 or 2)
%  numSource    : number of sources K
%  optVar       : optimization variable vector as above
%
%Output:
%  obj          : scalar objective value; returns +inf if S is not PD

arguments
  array {mustBeA(array, "cell")}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, "cell")}
  doaGrid {mustBeA(doaGrid, "cell")}
  dim (1,1) {mustBePositive, mustBeNumeric}
  numSource (1,1) {mustBePositive, mustBeInteger}
  optVar (1,:) {mustBeNumeric}
end

numArray = numel(array);
numDoaVars = dim * numSource;
optVar = optVar(:);
baseGrid = doaGrid{1};
gridType = baseGrid.type;

% ---- DOA decoding (always produce doaVec = Theta(:) in radians) ----
switch gridType
  case 'angle'
    doaVec = {optVar(1:numDoaVars)};

  case 'latlon'
    % doaParam(:) = [lat1; lon1; lat2; lon2; ...] (deg)
    latlonVec = optVar(1:numDoaVars);
    lat = latlonVec(1:2:end).';
    lon = latlonVec(2:2:end).';

    % altitude fixed at 0
    [x,y,z] = geodetic2ecef(baseGrid.spheroid, lat, lon, zeros(size(lat)));
    ecef = [x; y; z];  % 3xK

    doaVec = cell(size(array));
    for ii = 1:numArray
      doaVec{ii} = zeros(numDoaVars,1);
      doa = ecefToAngleGrid(ecef, doaGrid{ii}.arrayCenter);
      doaVec{ii}(1:2:end) = doa(1,:).';
      doaVec{ii}(2:2:end) = doa(2,:).';
    end

  otherwise
    error('nllStoUc:InvalidType', 'doaGrid.type must be ''angle'' or ''latlon''.');
end

% ---- power/noise decoding ----
signalPower = optVar(numDoaVars + (1:numSource*numArray));
noiseVar    = optVar(end-numArray+1:end);

% ---- covariance model ----
obj = 0;
for ii = 1:numArray
  Theta = reshape(doaVec{ii}, dim, numSource);                 % dim x K, radians
  A     = steeringMatrix(array{ii}, wavelength, Theta); % M x K
  M     = size(A, 1);

  idxP = (ii-1)*numSource + (1:numSource);
  S = A * bsxfun(@times, signalPower(idxP), A') + noiseVar(ii) * eye(M);
  S = 0.5 * (S + S'); % Hermitianize

  % ---- stable logdet/tr via Cholesky ----
  [L, flag] = chol(S, 'lower');
  if flag > 0
    obj = inf;
    return;
  end

  logdetS = 2 * sum(log(diag(L)));
  SinvR   = L' \ (L \ sampleCovMat{ii});
  obj     = obj + logdetS + real(trace(SinvR));
end
end
