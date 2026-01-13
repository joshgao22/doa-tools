function [doaResult, signalPower, noiseVar] = estimatorMleStoConGrid(array, ...
  wavelength, sampleCovMat, numSource, doaGrid, refineEstimates, maxIterations, subgridSize)
%ESTIMATORMLESTOCONGRID Stochastic ML (concentrated) DOA estimation on a predefined grid.
%
% This function estimates directions of arrival (DOAs) using a stochastic
% maximum-likelihood (SML) objective with closed-form (profile) elimination
% of nuisance parameters (source power and noise variance). A scanning
% spectrum is evaluated on a DOA grid, coarse peaks are detected, and
% optionally refined via iterative local subgrid search.
%
% Supported modes:
%   1) Single-array estimation: array/sampleCovMat/doaGrid are structs.
%   2) Multi-array joint estimation: array/sampleCovMat/doaGrid are cells.
%      The scanning spectrum is formed by summing per-array concentrated
%      scores (joint log-likelihood up to additive constants).
%
% After peak selection (and optional refinement), per-array closed-form
% estimates of:
%   - noise variance (sigma^2), and
%   - per-source power (diagonal of the estimated signal covariance)
% are computed at the final K=numSource DOA estimates.
%
%% Syntax
%   [doaResult, signalPower, noiseVar] = estimatorMleStoConGrid(array, wavelength, ...
%       sampleCovMat, numSource, doaGrid)
%   [doaResult, signalPower, noiseVar] = estimatorMleStoConGrid(array, wavelength, ...
%       sampleCovMat, numSource, doaGrid, refineEstimates, maxIterations, subgridSize)
%
%% Inputs
%   array           - Array geometry object(s). Either:
%                     * struct/object representing one array, or
%                     * 1-by-N cell array for multiple arrays.
%                     Each array must be compatible with steeringMatrix().
%   wavelength      - Signal wavelength (same unit as element positions).
%   sampleCovMat    - Sample covariance matrix/matrices. Either:
%                     * M-by-M complex matrix (single array), or
%                     * 1-by-N cell array, each is Mi-by-Mi (multi-array).
%                     Matrices are Hermitianized internally as 0.5*(R + R').
%   numSource       - Number of sources to estimate (positive integer, K).
%   doaGrid         - Scanning grid definition. Either:
%                     * grid struct (single array), or
%                     * 1-by-N cell array of grid structs (multi-array).
%                     Each grid struct must contain an angleGrid field for
%                     steering evaluation.
%   refineEstimates - Logical flag enabling local refinement of coarse peaks.
%                     Default: false.
%   maxIterations   - Refinement iterations (positive integer).
%                     Default: 10.
%   subgridSize     - Local grid resolution per refinement step.
%                     Default: 10 (scalar or 2-element vector depending on grid).
%
%% Outputs
%   doaResult       - Structure returned by findDoaFromSpectrum, optionally
%                     augmented with:
%                       .angleEst  : refined DOA estimates (if refined)
%                       .latlonEst : refined ground estimates (latlon grid only)
%                     It may also contain:
%                       .peakIndices, .isResolved, and other status/meta fields.
%   signalPower     - Stacked per-source power estimates across arrays.
%                     Size: (K*numArray)-by-1, with K entries per array in
%                     array order: [p_array1; p_array2; ...].
%                     Empty if doaResult.isResolved is false.
%   noiseVar        - Per-array noise variance estimates.
%                     Size: numArray-by-1. Empty if doaResult.isResolved is false.
%
%% Notes
%   - Rank condition: for each array i, Mi > numSource is required to estimate
%     noise variance using the closed-form projector expression.
%   - Multi-array fusion: the joint scanning spectrum is the sum of per-array
%     concentrated scores, assuming a common DOA hypothesis grid across arrays.
%   - Concentrated score: this implementation uses a monotone surrogate of the
%     (negative) log-likelihood, namely -logdet(S(theta)), where S is the
%     model covariance built from profile ML estimates of {p, sigma^2}.
%   - Grid evaluation: computeSpec must evaluate the steering vector at each
%     grid point gg (i.e., a(theta_gg)). Ensure the steeringMatrix call is
%     indexed per grid point if angleGrid stores multiple hypotheses.
%
%% See also
%   steeringMatrix, findDoaFromSpectrum, refineGridEstimate

arguments (Input)
  array (1,:) {mustBeA(array, ["struct","cell"])}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, ["double","cell"])}
  numSource (1,1) {mustBePositive, mustBeInteger}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct","cell"])}
  refineEstimates (1,1) logical = false
  maxIterations (1,1) {mustBePositive, mustBeInteger} = 10
  subgridSize (1,1) {mustBePositive, mustBeInteger} = 10
end

%% Normalize inputs to cell containers
% Wrap single-array inputs into 1-by-1 cells to unify downstream processing.
if ~iscell(array)
  array        = {array};
  sampleCovMat = {sampleCovMat};
  doaGrid      = {doaGrid};
end

%% Validate multi-array container lengths
numArray = numel(array);
if numel(sampleCovMat) ~= numArray || numel(doaGrid) ~= numArray
  error('estimatorMleStoConGrid:InvalidInput', ...
    'Cell inputs "array", "sampleCovMat", and "doaGrid" must have the same number of elements.');
end

%% Hermitianize covariance for each array
% Reduce asymmetry caused by finite-snapshot estimation.
for ii = 1:numArray
  sampleCovMat{ii} = 0.5 * (sampleCovMat{ii} + sampleCovMat{ii}');
end

%% Evaluate concentrated SML scanning spectrum on the DOA grid
spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid);

%% Detect coarse peaks on the grid
% In multi-array mode, doaGrid is a common hypothesis set; use doaGrid{1} for
% indexing and peak reporting.
doaResult = findDoaFromSpectrum(doaGrid{1}, spectrum, numSource);

%% Optional: local refinement around coarse peaks
% Refinement is applied only when the coarse solution is marked as resolved.
if isfield(doaResult, 'isResolved') && doaResult.isResolved && refineEstimates
  objFunc = @(grid) computeSpec(array, sampleCovMat, wavelength, grid);

  [angleEst, latlonEst] = refineGridEstimate(objFunc, doaGrid, doaResult, ...
    maxIterations, subgridSize);

  %% Post-process refined estimates
  % For single-array mode, sort by the first DOA component for stable ordering.
  if numArray == 1
    [~, sortIdx] = sort(angleEst(1,:));
    angleEst = angleEst(:, sortIdx);

    if ~isempty(latlonEst)
      latlonEst = latlonEst(:, sortIdx);
    end
  end

  doaResult.angleEst  = angleEst;
  doaResult.latlonEst = latlonEst;
end

%% Closed-form estimates of signal powers and noise variance at final DOAs
% Per-array estimates are computed given the final K=numSource angle estimates.
signalPower = [];
noiseVar    = [];

if isfield(doaResult,'isResolved') && doaResult.isResolved
  estAngles = doaResult.angleEst;

  noiseVar     = zeros(numArray, 1);
  signalPower  = zeros(numSource*numArray, 1);

  for ii = 1:numArray
    R = sampleCovMat{ii};
    steeringMat = steeringMatrix(array{ii}, wavelength, estAngles); % MxK
    M = size(R, 1);

    if M <= numSource
      error('estimatorMleStoConGrid:InvalidInput', ...
        'Need numSensors > numSource to estimate noise variance.');
    end

    G  = steeringMat' * steeringMat;  % KxK
    iG = pinv(G);

    % Noise variance estimator (projector-based)
    sigma2 = real(trace((eye(M) - steeringMat*iG*steeringMat') * R)) / (M - numSource);

    % Signal covariance and per-source power estimates
    Rs = real(iG * steeringMat' * R * steeringMat * iG - sigma2 * iG);
    p  = diag(Rs);
    p  = max(p, 0);

    noiseVar(ii) = max(sigma2, 0);

    idxP = (ii-1)*numSource + (1:numSource);
    signalPower(idxP) = p(:);
  end
end

end


function spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid)
%COMPUTESPEC Joint concentrated stochastic ML scanning spectrum over a DOA grid.
%
% This helper evaluates a concentrated (profile) stochastic ML objective at
% each grid point and forms a joint spectrum by summing per-array scores:
%   spectrum(theta) = sum_i score_i(theta)
% where each score_i(theta) is a monotone function of the SML log-likelihood
% after closed-form elimination of source power and noise variance (k=1 profile).
%
% Per array i and grid hypothesis theta:
%   Rhat_i           : sample covariance (Hermitian)
%   a_i(theta)       : steering vector
%   p_i(theta) >= 0  : profile ML source power (k=1)
%   sigma_i^2(theta) : profile ML noise variance
%   S_i(theta)       : model covariance = p_i a_i a_i^H + sigma_i^2 I
% Score (monotone in likelihood, up to constants):
%   score_i(theta) = -logdet(S_i(theta))
%
%% Syntax
%   spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid)
%
%% Inputs
%   array        - 1-by-N cell array of array structs/objects (steeringMatrix compatible)
%   sampleCovMat - 1-by-N cell array of sample covariance matrices (Mi-by-Mi, Hermitian)
%   wavelength   - Signal wavelength (scalar)
%   doaGrid      - 1-by-N cell array of grid structs containing .angleGrid
%
%% Output
%   spectrum     - 1-by-Ngrid joint scanning spectrum aligned with
%                  doaGrid{1}.angleGrid linearization.
%
%% Notes
%   - Common hypothesis grid: this function assumes all grids share identical
%     grid-point ordering across arrays.
%   - Numerical robustness: Cholesky factorization is used to compute logdet;
%     hypotheses yielding non-PD S are treated as -Inf.
%   - Grid indexing: if doaGrid{ii}.angleGrid stores multiple hypotheses, the
%     steering vector must be evaluated at the current grid point (gg).
%
%% See also
%   steeringMatrix

arguments
  array (1,:) cell
  sampleCovMat (1,:) cell
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  doaGrid (1,:) cell
end

%% Initialize spectrum accumulation
% All arrays are assumed to share a common hypothesis grid with identical
% point ordering.
numArray = numel(array);
numGrid  = size(doaGrid{1}.angleGrid, 2);
spectrum = zeros(1, numGrid);

for gg = 1:numGrid

  jointScore = 0;
  for ii = 1:numArray
    R    = sampleCovMat{ii};
    M    = size(R, 1);
    eyeM = eye(M);

    % Steering vector at the current grid point (theta_gg).
    % NOTE: ensure steeringMatrix is evaluated for the gg-th hypothesis if
    % angleGrid contains multiple grid points.
    theta = doaGrid{ii}.angleGrid(:, gg);
    a = steeringMatrix(array{ii}, wavelength, theta); % Mx1
    aHa = real(a' * a);

    if aHa <= 0
      jointScore = -Inf;
      break;
    end

    invaHa = 1 / aHa;

    %% k=1 profile ML estimates (p, sigma^2) under S = p a a^H + sigma^2 I
    noiseVar = real(trace((eyeM - a*(invaHa*a')) * R)) / (M - 1);
    p        = real(invaHa * (a' * R * a) * invaHa - noiseVar * invaHa);

    p        = max(p, 0);
    noiseVar = max(noiseVar, 0);

    %% Build model covariance and evaluate concentrated score
    S = a * p * a' + noiseVar * eyeM;
    S = 0.5 * (S + S');  % enforce Hermitian symmetry

    % Concentrated score (monotone in likelihood): -logdet(S)
    % Use Cholesky for stability; if not PD, treat as -Inf.
    [L, flag] = chol(S, 'lower');
    if flag > 0
      jointScore = -Inf;
      break;
    end
    logdetS   = 2 * sum(log(diag(L)));
    jointScore = jointScore - real(logdetS);
  end

  spectrum(gg) = jointScore;
end

end
