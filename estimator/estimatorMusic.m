function doaResult = estimatorMusic(array, wavelength, sampleCovMat, ...
  numSource, doaGrid, refineEstimates, maxIterations, subgridSize)
%ESTIMATORMUSIC MUSIC-based DOA estimation on a predefined grid.
%
% This function estimates directions of arrival (DOAs) using the classical
% MUSIC spatial spectrum:
%   P(theta) = 1 / || Un^H a(theta) ||^2
% where Un is the noise-subspace eigenvector matrix and a(theta) is the
% steering vector/matrix evaluated on a scanning grid.
%
% The function supports:
%   1) Single-array MUSIC: array/sampleCovMat/doaGrid are structs.
%   2) Multi-array incoherent MUSIC fusion: array/sampleCovMat/doaGrid are cells,
%      and the spectrum is formed by summing denominator terms.
%
% Optionally, coarse grid peaks can be refined by iterative local subgrid
% search (zoom-in refinement) using refineGridEstimate.
%
%% Syntax
%   doaResult = estimatorMusic(array, wavelength, sampleCovMat, numSource, doaGrid)
%   doaResult = estimatorMusic(array, wavelength, sampleCovMat, numSource, doaGrid, ...
%                              refineEstimates, maxIterations, subgridSize)
%
%% Inputs
%   array           - Array geometry object(s). Either:
%                     * struct/object representing one array, or
%                     * 1-by-N cell array for multiple arrays.
%                     Each array must provide:
%                       .count (number of elements) and geometry for steeringMatrix().
%   wavelength      - Signal wavelength (same unit as element positions).
%   sampleCovMat    - Sample covariance matrix/matrices. Either:
%                     * M-by-M complex matrix (single array), or
%                     * 1-by-N cell array, each is Mi-by-Mi (multi-array).
%                     Matrices are Hermitianized internally as 0.5*(R + R').
%   numSource       - Number of sources to estimate (positive integer).
%   doaGrid         - Scanning grid definition. Either:
%                     * grid struct (single array), or
%                     * 1-by-N cell array of grid structs (multi-array).
%                     Each grid struct must be compatible with genDoaGrid()
%                     and must contain an angleGrid field for steering evaluation.
%   refineEstimates - Logical flag enabling local refinement of coarse peaks.
%                     Default: false.
%   maxIterations   - Refinement iterations (positive integer).
%                     Default: 10.
%   subgridSize     - Local grid resolution per refinement step.
%                     Default: 10 (scalar or 2-element vector depending on grid).
%
%% Outputs
%   doaResult       - Structure returned by findDoaFromSpectrum, augmented with:
%                       .angleEst  : refined local DOA estimates (if refined)
%                       .latlonEst : refined ground estimates (latlon grid only)
%                     It may also contain:
%                       .peakIndices, .isResolved, and other status/meta fields
%                     depending on the findDoaFromSpectrum implementation.
%
%% Notes
%   - Requirement: for each array, the number of elements must exceed numSource.
%   - Eigenvalue sorting: eigenvalues are sorted in ascending order to extract
%     the noise subspace, rather than relying on eig() output ordering.
%   - Multi-array fusion: an incoherent fusion is performed on the denominator:
%       denom_total = sum_i || Un_i^H a_i(theta) ||^2
%     This assumes a common DOA hypothesis (same theta/latlon) and combines
%     evidence without phase alignment across arrays.
%   - Refinement: when enabled and doaResult.isResolved is true, refinement is
%     performed around coarse peaks using refineGridEstimate. After refinement,
%     estimates are sorted by the first DOA component for stable ordering.
%
%% See also
%   computeSpec, steeringMatrix, findDoaFromSpectrum, refineGridEstimate

arguments (Input)
  array (1,:) {mustBeA(array, ["struct", "cell"])}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, ["double", "cell"])}
  numSource (1,1) {mustBePositive, mustBeInteger}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  refineEstimates (1,1) logical = false
  maxIterations (1,1) {mustBePositive, mustBeInteger} = 10
  subgridSize {mustBePositive, mustBeInteger} = 10
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
  error('estimatorMusic:InvalidInput', ...
    'Cell inputs "array", "sampleCovMat", and "doaGrid" must have the same number of elements.');
end

%% Extract noise subspace by eigendecomposition
% Hermitianize sample covariance to reduce asymmetry from finite snapshots.
noiseSubspace = cell(size(array));

for ii = 1:numArray
  R = 0.5 * (sampleCovMat{ii} + sampleCovMat{ii}');

  % Eigen-decomposition with vector eigenvalue output for explicit sorting.
  [eigVecs, eigVals] = eig(R, 'vector');

  % Sort eigenvalues in ascending order and take the smallest (M - numSource).
  [~, sortIdx] = sort(real(eigVals), 'ascend');
  noiseSubspace{ii} = eigVecs(:, sortIdx(1:end - numSource));
end

%% Evaluate MUSIC spectrum on the scanning grid
spectrum = computeSpec(array, noiseSubspace, wavelength, doaGrid);

%% Detect coarse peaks on the grid
% In multi-array mode, doaGrid is a common hypothesis set; use doaGrid{1} for
% indexing and peak reporting.
doaResult = findDoaFromSpectrum(doaGrid{1}, spectrum, numSource);

%% Optional: local refinement around coarse peaks
% Refinement is applied only when the coarse solution is marked as resolved.
if isfield(doaResult, 'isResolved') && doaResult.isResolved && refineEstimates
  objFunc = @(grid) computeSpec(array, noiseSubspace, wavelength, grid);

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

end


function spec = computeSpec(array, noiseSubspace, wavelength, doaGrid)
%COMPUTESPEC Compute MUSIC spatial spectrum over a DOA grid.
%
% This helper computes the MUSIC pseudo-spectrum:
%   P = 1 ./ max(denom, eps)
% where denom = || Un^H A ||^2 evaluated per grid point.
%
% Supported modes:
%   1) Single array: inputs are structs/matrices (wrapped as cells upstream).
%   2) Multi-array incoherent fusion: inputs are cells of equal length, and the
%      combined denominator is formed as a sum over arrays.
%
%% Syntax
%   spec = computeSpec(array, noiseSubspace, wavelength, doaGrid)
%
%% Inputs
%   array          - 1-by-N cell array of array structs/objects
%   noiseSubspace  - 1-by-N cell array of Un matrices
%   wavelength     - Signal wavelength (scalar)
%   doaGrid        - 1-by-N cell array of grid structs containing .angleGrid
%
%% Output
%   spec           - 1-by-Ngrid spectrum values aligned with the linearization
%                    used by doaGrid{1}.angleGrid (and genDoaGrid).
%
%% Notes
%   - Common hypothesis grid: this function assumes all grids share the same
%     number and ordering of grid points across arrays.
%   - Numerical robustness: real(sum(conj(proj).*proj)) removes small imaginary
%     parts caused by numerical roundoff.
%
%% See also
%   steeringMatrix

arguments
  array (1,:) cell
  noiseSubspace (1,:) cell
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  doaGrid (1,:) cell
end

%% Initialize spectrum accumulation
% All arrays are assumed to share a common hypothesis grid with identical
% point ordering.
numArray = numel(array);
numGrid  = size(doaGrid{1}.angleGrid, 2);
invSpec  = zeros(1, numGrid);

%% Accumulate MUSIC denominators across arrays
for ii = 1:numArray
  % Steering matrix evaluated on the local grid realization.
  A = steeringMatrix(array{ii}, wavelength, doaGrid{ii}.angleGrid);

  % Projection onto the noise subspace:
  %   Un^H A -> (M - numSource)-by-Ngrid
  % Denominator contribution is the squared 2-norm per grid point.
  proj    = noiseSubspace{ii}' * A;
  invSpec = invSpec + real(sum(conj(proj) .* proj, 1));
end

%% Form MUSIC pseudo-spectrum
% Use eps to avoid division by zero at exact nulls.
spec = 1 ./ max(invSpec, eps);

end
