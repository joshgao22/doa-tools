function [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, ...
  wavelen, sampleCovMat, numSource, doaGrid, initParam, verbose)
%ESTIMATORMLESTOCONOPT Stochastic ML (concentrated) DOA estimation via continuous optimization.
%
% This function estimates directions of arrival (DOAs) by directly minimizing
% the stochastic uncorrelated-source negative log-likelihood:
%   f(Theta, p, sigma^2) = logdet(S) + tr(S^{-1} Rhat)
% subject to nonnegativity constraints:
%   p >= 0,  sigma^2 >= 0,
% with the covariance model
%   S = A(Theta) diag(p) A(Theta)^H + sigma^2 I.
%
% The DOA parameters Theta are optimized in the continuous domain (not limited
% to discrete grid points). The doaGrid input is used to define parameter
% bounds and (optionally) construct an initial guess.
%
% Supported modes:
%   1) Single-array optimization: array/sampleCovMat/doaGrid are structs.
%   2) Multi-array joint optimization: array/sampleCovMat/doaGrid are cells.
%      The objective is the sum over arrays (joint likelihood up to constants).
%
% Initialization:
%   - If initParam is empty, an MVDR-based initializer is attempted.
%   - For grid type 'angle', Theta is parameterized directly in radians.
%   - For grid type 'latlon', optimization variables are (lat, lon) in degrees,
%     which are mapped to local (az, el) angles (radians) per array when forming
%     steering matrices.
%
%% Syntax
%   [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, wavelen, ...
%       sampleCovMat, numSource, doaGrid)
%   [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, wavelen, ...
%       sampleCovMat, numSource, doaGrid, initParam)
%   [doaResult, signalPower, noiseVar] = estimatorMleStoConOpt(array, wavelen, ...
%       sampleCovMat, numSource, doaGrid, initParam, verbose)
%
%% Inputs
%   array        - Array geometry object(s). Either:
%                  * struct/object representing one array, or
%                  * 1-by-N cell array for multiple arrays.
%                  Each array must be compatible with steeringMatrix().
%   wavelen      - Signal wavelength (scalar, same unit as array geometry).
%   sampleCovMat - Sample covariance matrix/matrices. Either:
%                  * M-by-M complex matrix (single array), or
%                  * 1-by-N cell array, each is Mi-by-Mi (multi-array).
%                  Matrices are Hermitianized internally as 0.5*(R + R').
%   numSource    - Number of sources (positive integer, K).
%   doaGrid      - Grid definition struct(s) used for:
%                  * parameterization selection via .type in {'angle','latlon'}
%                  * DOA dimension via .dimension in {1,2}
%                  * box bounds via .range
%                  Additional fields required for 'latlon':
%                    .spheroid, .arrayCenter
%   initParam    - (Optional) Initial optimization vector. If empty, an
%                  initializer is constructed internally.
%                  Layout: [doaParam(:); p; noiseVar], where doaParam depends
%                  on doaGrid.type (see Notes).
%   verbose      - (Optional) Logical flag to show fmincon iterations.
%                  Default: false.
%
%% Outputs
%   doaResult    - MUSIC-like result struct with fields:
%                  .doaGrid     : input grid cell/struct (for metadata)
%                  .spectrum    : [] (continuous optimization, no scan spectrum)
%                  .peakIndices : [] (not applicable)
%                  .angleEst    : dim-by-K DOA estimates in radians (for 'angle';
%                                for 'latlon' only set when numArray==1)
%                  .latlonEst   : 2-by-K (lat;lon) in degrees when grid type is 'latlon'
%                  .isResolved  : true if optimization exitflag indicates success
%                  .isDiscrete  : flag indicating discrete vs continuous output
%                                (see Notes regarding current implementation)
%   signalPower  - Estimated source powers stacked across arrays.
%                  Size: (K*numArray)-by-1, with K entries per array in array order.
%   noiseVar     - Estimated noise variance per array.
%                  Size: numArray-by-1.
%
%% Notes
%   - Objective: this implements the standard stochastic ML criterion for
%     uncorrelated sources. The objective is summed across arrays in multi-array mode.
%   - DOA parameterization:
%       * doaGrid.type == 'angle'  : doaParam(:) stores angles in radians.
%         For dim=2, stacking is [az1; el1; az2; el2; ...] (column-major).
%       * doaGrid.type == 'latlon' : doaParam(:) stores [lat1; lon1; lat2; lon2; ...]
%         in degrees; each array maps (lat,lon) to local (az,el) radians using
%         geodetic2ecef + ecefToAngleGrid relative to doaGrid{ii}.arrayCenter.
%   - Bounds: only box constraints are used (see setDoaBounds). Power and noise
%     variables are constrained by [0, +inf).
%   - Output consistency: for grid type 'latlon' and numArray>1, angleEst is left
%     empty because local (az,el) depends on arrayCenter; latlonEst is the shared
%     global estimate.
%   - Implementation note: doaResult.isDiscrete is currently set to true in code,
%     but this solver produces continuous-domain estimates. Consider setting it
%     to false for semantic consistency (not changed here).
%
%% See also
%   fmincon, steeringMatrix, estimatorMvdr, setDoaBounds, nllStoUc, ecefToAngleGrid

arguments (Input)
  array (1,:) {mustBeA(array, ["struct", "cell"])}
  wavelen (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, ["double", "cell"])}
  numSource (1,1) {mustBePositive, mustBeInteger}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  initParam (1,:) = []
  verbose (1,1) logical = false
end

% Normalize inputs to cell containers
% Wrap single-array inputs into 1-by-1 cells to unify downstream processing.
if ~iscell(array)
  array        = {array};
  sampleCovMat = {sampleCovMat};
  doaGrid      = {doaGrid};
end

% Validate multi-array container lengths
numArray = numel(array);
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
if numArray == 1
  doaResult.doaGrid = doaGrid{1};
else
  doaResult.doaGrid = doaGrid;
end
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
% This helper fills the lower/upper bounds of the DOA-related entries in the
% optimization vector, according to the domain specified by doaGrid.range.
% Only DOA variables are modified; bounds for source power and noise variance
% are handled outside this function.
%
% DOA variable stacking (MATLAB column-major):
%   Let K = numSource, dim ∈ {1,2}.
%   - dim = 1 (1D): doaParam(:) = [theta_1; theta_2; ...; theta_K]
%   - dim = 2 (2D): doaParam(:) = [v1_1; v2_1; v1_2; v2_2; ...; v1_K; v2_K]
%     where (v1,v2) is:
%       * (az, el) in radians when doaGrid.type == 'angle'
%       * (lat, lon) in degrees when doaGrid.type == 'latlon'
%
% Bounds:
%   - dim = 1:
%       theta ∈ [doaGrid.range(1), doaGrid.range(2)]
%   - dim = 2:
%       v1 ∈ [doaGrid.range(1,1), doaGrid.range(1,2)]
%       v2 ∈ [doaGrid.range(2,1), doaGrid.range(2,2)]
%
%% See also
%   estimatorMleStoConOpt, nllStoUc

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
%NLLSTOUC Joint stochastic uncorrelated-source negative log-likelihood.
%
% This helper evaluates the per-array stochastic ML objective and sums across
% arrays (joint criterion up to constants):
%   obj = sum_i [ logdet(S_i) + tr(S_i^{-1} Rhat_i) ].
%
% Per array i, the covariance model is
%   S_i = A_i(Theta_i) diag(p_i) A_i(Theta_i)^H + sigma_i^2 I,
% where p_i >= 0 (source powers) and sigma_i^2 >= 0 (noise variance).
%
% The optimization vector layout is:
%   optVar = [doaParam(:); signalPower; noiseVar]
% where:
%   - doaParam(:) contains K DOA hypotheses (dim-by-K) in either radians ('angle')
%     or degrees ('latlon'), stacked column-wise.
%   - signalPower stacks K powers per array: (K*numArray)-by-1.
%   - noiseVar stacks one variance per array: numArray-by-1.
%
% DOA decoding:
%   - doaGrid{1}.type == 'angle':
%       doaParam(:) is interpreted directly as local DOA parameters (radians).
%   - doaGrid{1}.type == 'latlon':
%       doaParam(:) stores [lat1; lon1; lat2; lon2; ...] (degrees). The shared
%       global (lat,lon) is mapped to local angles (az,el) radians for each array
%       using geodetic2ecef() + ecefToAngleGrid() w.r.t. doaGrid{ii}.arrayCenter.
%
% Numerical evaluation:
%   - Cholesky factorization is used for stable logdet and linear solves.
%   - If any S_i is not positive definite, obj returns +Inf.
%
%% See also
%   steeringMatrix, geodetic2ecef, ecefToAngleGrid, setDoaBounds

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
