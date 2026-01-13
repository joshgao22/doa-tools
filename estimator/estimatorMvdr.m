function doaResult = estimatorMvdr(array, wavelength, sampleCovMat, ...
  numSource, doaGrid, refineEstimates, maxIterations, subgridSize)
%ESTIMATORMVDR Perform MVDR-based DOA estimation (top-level interface).
% Estimates the directions of arrival (DOAs) of multiple sources using the
% classical MVDR (Capon) beamformer spectrum. It supports grid search and
% optional local refinement of the estimated DOAs.
%
%Syntax:
%   doaResult = estimatorMvdr(array, wavelength, sampleCovMat, ...
%                numSource, doaGrid, refineEstimates, maxIterations, subgridSize)
%
%Inputs:
%   array           - Array geometry struct, providing element positions and count.
%   wavelength      - Signal wavelength (same unit as array positions).
%   sampleCovMat    - Sample covariance matrix (M×M), Hermitian.
%   numSource       - Number of sources to estimate (positive integer).
%   doaGrid         - DOA grid structure defining the scanning grid.
%   refineEstimates - (Optional) Logical flag to enable local refinement
%                     around coarse estimates (default: false).
%   maxIterations   - (Optional) Maximum number of iterations for refinement
%                     (default: 10).
%   subgridSize     - (Optional) Subgrid resolution per dimension for refinement
%                     (default: 10).
%
%Output:
%   doaResult       - DOA result structure containing the estimated DOAs and status flag.
%
%Notes:
%   - MVDR requires a well-conditioned covariance estimate; diagonal loading is used.
%   - 1D "ground" grid is not supported by design (use 2D lat/lon ground grid).
%
%See also: computeSpec, findDoaFromSpectrum, refineGridEstimate

arguments (Input)
  array {mustBeA(array, ["struct", "cell"])}
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  sampleCovMat {mustBeA(sampleCovMat, ["double", "cell"])}
  numSource (1,1) {mustBePositive, mustBeInteger}
  doaGrid {mustBeA(doaGrid, ["struct", "cell"])}
  refineEstimates (1,1) logical = false
  maxIterations (1,1) {mustBePositive, mustBeInteger} = 10
  subgridSize {mustBePositive, mustBeInteger} = 10
end

%% Input consistency checks
% Enforce consistent container types (struct vs cell) and matching lengths.
if ~iscell(array)
  array = {array};
  sampleCovMat = {sampleCovMat};
  doaGrid = {doaGrid};
end

numArray = numel(array);

if numel(sampleCovMat) ~= numArray || numel(doaGrid) ~= numArray
  error('estimatorMusic:InvalidInput', ...
    'Cell inputs "array", "sampleCovMat", and "doaGrid" must have the same number of elements.');
end

% Hermitianize first to reduce numerical asymmetry in sample covariance.
for ii = 1:numArray
  sampleCovMat{ii} = 0.5 * (sampleCovMat{ii} + sampleCovMat{ii}'); 
end

%% Compute MVDR spectrum on the scanning grid
spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid);

%% Peak detection on the coarse grid
% In multi-array mode, the spectrum corresponds to the common hypothesis grid.
% Here we use doaGrid{1} as the representative grid for peak indexing.
doaResult = findDoaFromSpectrum(doaGrid{1}, spectrum, numSource);

%% Optional: refine coarse peaks by local zoom-in search
if isfield(doaResult, 'isResolved') && doaResult.isResolved && refineEstimates
  objFunc = @(grid) computeSpec(array, sampleCovMat, wavelength, grid);

  [angleEst, latlonEst] = refineGridEstimate(objFunc, doaGrid, doaResult, ...
    maxIterations, subgridSize);

  if numArray == 1
    % Sort by the first DOA component for consistency
    [~, sortIdx] = sort(angleEst(1,:));
    angleEst     = angleEst(:, sortIdx);
  
    if ~isempty(latlonEst)
      latlonEst  = latlonEst(:, sortIdx);
    end
  end

  doaResult.angleEst  = angleEst;
  doaResult.latlonEst = latlonEst;
end

end


function spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid)
%COMPUTESPEC Compute MVDR spatial spectrum over a DOA grid.
% spectrum = computeSpec(array, sampleCovMat, wavelength, doaGrid)
%
%Inputs:
%   array        - Array struct
%   sampleCovMat - Sample covariance matrix (M×M), Hermitian
%   wavelength   - Signal wavelength (scalar)
%   doaGrid      - DOA grid structure (provides .angleGrid)
%
%Output:
%   spectrum     - MVDR spectrum (1×N), N matches number of grid points
%
%Notes:
%   - MVDR spectrum: P = 1 / real(a^H * R^{-1} * a)
%   - Linear solves are used instead of explicit matrix inversion.

arguments
  array (1,:) cell
  sampleCovMat (1,:) cell
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  doaGrid (1,:) cell
end

numArray = numel(array);

%% MVDR spectrum evaluation
% Assume a common hypothesis grid across arrays (same number/order of points).
numGrid = size(doaGrid{1}.angleGrid, 2);
invSpec = zeros(1, numGrid);

for ii = 1:numArray
  A = steeringMatrix(array{ii}, wavelength, doaGrid{ii}.angleGrid);
  
  % a^H R^{-1} a via linear solve
  X = sampleCovMat{ii} \ A;                         
  invSpec = invSpec + real(sum(conj(A) .* X, 1));  
end

invSpec(invSpec <= 0) = eps;
spectrum = 1 ./ invSpec;

end
