function [refinedAngleEst, refinedLatlonEst] = refineGridEstimate( ...
  objFunction, doaGrid, doaResult, maxIterations, subgridSize)
%REFINEGRIDESTIMATE Refine DOA peak estimates via iterative local grid search.
%
% This function refines coarse peak locations (given by doaResult.peakIndices)
% by repeatedly constructing a small local subgrid around the current peak,
% evaluating an objective spectrum on that subgrid, and shrinking the search
% bounds toward the dominant peak.
%
%% Syntax
%   refinedAngleEst = refineGridEstimate(objFunction, doaGrid, doaResult)
%   [refinedAngleEst, refinedLatlonEst] = refineGridEstimate( ...
%       objFunction, doaGrid, doaResult, maxIterations, subgridSize)
%% Description
%   refinedAngleEst = refineGridEstimate(...) refines the DOA estimates in
%   local array coordinates (angle grid). For 1D grids, the output is a
%   1-by-K vector of angles. For 2D angle grids, the output is a 2-by-K
%   matrix [az; el].
%
%   [refinedAngleEst, refinedLatlonEst] = refineGridEstimate(...) also
%   returns refined estimates in ground coordinates (lat/lon) when the grid
%   type is 'latlon'. For a single latlon grid, refinedAngleEst is 2-by-K
%   and refinedLatlonEst is 2-by-K. For a cell array of latlon grids
%   (multiple arrays), refinedAngleEst is a 1-by-N cell, each cell contains
%   a 2-by-K local-angle estimate, and refinedLatlonEst is 2-by-K.
%% Inputs
%   objFunction    - Objective function handle that evaluates the spectrum
%                    over a candidate grid.
%                    Signature:
%                      spectrum = objFunction(candidateGrid)
%                    where candidateGrid is:
%                      - a grid struct (angle/latlon), or
%                      - a cell array of grid structs (latlon only)
%   doaGrid        - Coarse DOA grid definition. Either:
%                    (1) Struct for a single grid, or
%                    (2) Cell array of structs for multi-array latlon grids.
%                    Required fields:
%                      .dimension : 1 or 2
%                      .type      : 'angle' or 'latlon'
%                    For latlon grids, genDoaGrid may additionally require:
%                      .arrayCenter, .spheroid
%   doaResult      - Coarse search result containing detected peak indices.
%                    Required field:
%                      .peakIndices : linear indices of detected peaks
%   maxIterations  - Number of refinement iterations (positive integer).
%                    Default: 10.
%   subgridSize    - Number of grid points used at each refinement step.
%                    For 1D: scalar N.
%                    For 2D: scalar N (interpreted as [N, N]) or [n1, n2].
%                    Interpretation depends on grid type:
%                      - angle : [nAz, nEl]
%                      - latlon: [nLat, nLon]
%% Outputs
%   refinedAngleEst  - Refined estimates in local array coordinates.
%                      If doaGrid.type = 'angle':
%                        - dim=1: 1-by-K
%                        - dim=2: 2-by-K  ([az; el])
%                      If doaGrid.type = 'latlon':
%                        - single grid: 2-by-K local [az; el]
%                        - cell grid  : 1-by-N cell, each is 2-by-K
%   refinedLatlonEst - Refined estimates in ground coordinates [lat; lon].
%                      Only populated when doaGrid.type = 'latlon'.
%                      Otherwise returned as [].
%% Notes
%   - 1D refinement: starts from the coarse peak index, builds a local 1D
%     subgrid within neighbor bounds, selects the max, then shrinks bounds
%     to the new peak's neighbors; repeats for maxIterations.
%   - 2D refinement: evaluates the objective on a local 2D subgrid and selects
%     the dominant peak using regional maxima (imregionalmax). If no regional
%     maximum exists, falls back to the global maximum.
%   - Bound handling: neighbor bounds are clipped to grid limits; this
%     implementation does not perform angle wrapping (e.g., azimuth periodicity).
%   - Grid ordering: this function reshapes the spectrum to match the
%     grid generation convention used in genDoaGrid. In this implementation:
%       * angle grid uses reshape(objSpectrum, numEl, numAz) -> [el, az]
%       * latlon grid uses reshape(objSpectrum, numLon, numLat) -> [lon, lat]
%     Ensure objFunction returns spectrum in the same linearization order as
%     genDoaGrid outputs its grid points.
%   - Multi-array latlon (cell doaGrid): all latlon candidate grids are built
%     using the same [lat, lon] bounds, while local angle grids differ per array.
%% See also
%   genDoaGrid, getNeighbourGrid, imregionalmax, reshape, sub2ind

%% Argument validation (MATLAB function arguments block)
% - Use arguments block to enforce basic type/shape constraints early.
% - objFunction must be a scalar function handle.
% - doaGrid can be either:
%     (1) a struct  : single grid (angle or latlon)
%     (2) a cell    : multiple grids (only supported for 'latlon' mode)
% - doaResult must contain doaResult.peakIndices (linear indices on doaGrid).
% - maxIterations controls how many zoom-in steps to run.
% - subgridSize controls the resolution of each zoom-in grid:
%     * 1D: scalar
%     * 2D: scalar or 2-element vector
arguments
  objFunction (1,1) function_handle
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  doaResult (1,1) struct
  maxIterations (1,1) {mustBePositive, mustBeInteger} = 10
  subgridSize {mustBePositive, mustBeInteger} = 10
end

%% Grid meta parsing: single grid vs multi-grid (cell)
% isGridCell indicates whether we are refining over:
%   - a single grid struct (angle / latlon), or
%   - a cell array of grids (multi-array latlon refinement)
if ~iscell(doaGrid)
  doaGrid = {doaGrid};
end

% Multi-array mode: doaGrid is a cell array, each element corresponds to
% one array (or one sensor platform). We assume all grids share the same
% (lat, lon) search region and the same discretization in lat/lon.
%
% NOTE: In this implementation, cell-grid refinement is restricted to
%       'latlon' only. Local angle grids (angleGrid) are array-dependent.
numGridCell = numel(doaGrid);
baseGrid = doaGrid{1};
dim = baseGrid.dimension;
gridType = lower(string(baseGrid.type));
if gridType == "angle" && numGridCell ~= 1
  error('only ''latlon'' supported with multiple cell grid')
end

%% Initialize outputs (filled depending on gridType)
% Convention:
%   - If gridType == "angle": refinedAngleEst is populated; refinedLatlonEst=[]
%   - If gridType == "latlon": refinedLatlonEst is populated; refinedAngleEst
%     is also populated as local DOA (for single grid) or a cell array (for
%     multi-array grids).
refinedAngleEst  = [];
refinedLatlonEst = [];

%% Parse coarse peaks to be refined
% doaResult.peakIndices are linear indices w.r.t. the coarse grid ordering.
% Convert to a column vector for consistent loop handling.
peakIndices = doaResult.peakIndices(:);
numEst = numel(peakIndices);

%% Refinement by grid dimension
% dim = 1:
%   - only supports 'angle' grid; refinement is interval shrinking.
% dim = 2:
%   - supports 'angle' and 'latlon'; refinement is 2D subgrid evaluation
%     + dominant-peak selection (regional max preferred).
switch dim
  case 1
    %% SD refinement: angle only
    % 1D search is interpreted as a single angular parameter (e.g., azimuth
    % or a scan angle). latlon is not meaningful in 1D here.
    if gridType ~= "angle"
      error('refineGridEstimate:InvalidGrid', ...
        '1D refinement only supports doaGrid.type = ''angle''.');
    end
    if ~isscalar(subgridSize)
      error('refineGridEstimate:InputDimensionNotMatch', ...
        'subgridSize must be scalar for 1D grid.');
    end

    % Output: 1-by-K refined angles
    refinedAngleEst = zeros(1, numEst);

    for estIdx = 1:numEst
      % Initialize neighbor bounds around the coarse peak.
      % Expected format from getNeighbourGrid (1D): [lower, upper]
      neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

      for iter = 1:maxIterations
        % Construct a local grid within current bounds and evaluate spectrum.
        candidateGrid = genDoaGrid('angle', 1, subgridSize, neighbourGrid);

        objSpectrum = objFunction({candidateGrid});

        % Pick the best point on this local grid (global max in 1D).
        [~, peakIdx] = max(objSpectrum);

        % Shrink bounds to the neighbors of the new peak (zoom-in).
        neighbourGrid = getNeighbourGrid(candidateGrid, peakIdx);
      end

      % Store final refined estimate (scalar angle at the last peakIdx).
      refinedAngleEst(estIdx) = candidateGrid.angleGrid(peakIdx);
    end

  case 2
    %% 2D refinement: angle or latlon
    switch gridType
      case "angle"
        %% 2D angle refinement: [az, el]
        % subgridSize handling:
        %   - scalar N -> use NxN grid (numAz=numEl=N)
        %   - [nAz, nEl] -> rectangular local grid
        if isscalar(subgridSize)
          numAz = subgridSize; numEl = subgridSize;
        elseif isvector(subgridSize) && numel(subgridSize) == 2
          numAz = subgridSize(1); numEl = subgridSize(2);
        else
          error('refineGridEstimate:InputDimensionNotMatch', ...
            'subgridSize must be scalar or 2-element vector for 2D angle grid.');
        end

        % Output: 2-by-K, column is [az; el]
        refinedAngleEst = zeros(2, numEst);

        for estIdx = 1:numEst
          % Neighbor bounds around coarse peak.
          % Expected format from getNeighbourGrid (2D): 2x2 bounds,
          % e.g., [azMin azMax; elMin elMax] (depending on your implementation).
          neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

          for iter = 1:maxIterations
            % Local 2D angle grid within neighbor bounds
            candidateGrid = genDoaGrid('angle', 2, subgridSize, neighbourGrid);

            % Evaluate objective spectrum on this local grid
            objSpectrum = objFunction({candidateGrid});

            % Reshape to 2D matrix for peak picking.
            % Convention here: reshape to [el, az] = (numEl x numAz).
            % This must match the linear indexing order used by genDoaGrid.
            objSpectrum2D = reshape(objSpectrum, numEl, numAz);

            % Pick dominant peak: prefer regional maxima; fallback to global max
            linearIdx = pickDominantPeak(objSpectrum2D, [numEl, numAz]);

            % Recenter (shrink) bounds around the selected peak
            neighbourGrid = getNeighbourGrid(candidateGrid, linearIdx);
          end

          % Save refined DOA at final iteration
          refinedAngleEst(:, estIdx) = candidateGrid.angleGrid(:, linearIdx);
        end

      case "latlon"
        %% 2D latlon refinement: [lat, lon] with optional multi-array support
        % subgridSize handling:
        %   - scalar N -> use NxN grid (numLat=numLon=N)
        %   - [nLat, nLon] -> rectangular local grid
        if isscalar(subgridSize)
          numLat = subgridSize; numLon = subgridSize;
        elseif isvector(subgridSize) && numel(subgridSize) == 2
          numLat = subgridSize(1); numLon = subgridSize(2);
        else
          error('refineGridEstimate:InputDimensionNotMatch', ...
            'subgridSize must be scalar or 2-element vector for 2D latlon grid.');
        end

        % Preallocate outputs depending on single-grid vs cell-grid
        % Ground estimate is shared across arrays -> 2-by-K [lat; lon]
        refinedLatlonEst = zeros(2, numEst);

        % Local angle estimate is array-dependent -> 1-by-N cell, each 2-by-K
        refinedAngleEst = cell(1, numGridCell);
        for ii = 1:numGridCell
          refinedAngleEst{ii} = zeros(2, numEst);
        end

        for estIdx = 1:numEst
          % Use the first grid to define (lat, lon) neighbor bounds.
          % Assumption: all doaGrid{ii} share identical lat/lon lattice and bounds.
          neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

          for iter = 1:maxIterations
            % Build candidate grids for each array using the SAME lat/lon bounds,
            % but different arrayCenter/spheroid as needed.
            candidateGrid = cell(size(doaGrid));
            for ii = 1:numGridCell
              candidateGrid{ii} = genDoaGrid('latlon', 2, subgridSize, neighbourGrid, ...
                doaGrid{ii}.arrayCenter, doaGrid{ii}.spheroid);
            end

            % Evaluate objective using the collection of grids
            objSpectrum = objFunction(candidateGrid);

            % Reshape to 2D for peak picking.
            % Convention here: reshape to [lon, lat] = (numLon x numLat),
            % consistent with your original code and downstream indexing.
            objSpectrum2D = reshape(objSpectrum, numLon, numLat);

            % Pick dominant peak on the ground grid
            linearIdx = pickDominantPeak(objSpectrum2D, [numLon, numLat]);

            % Update bounds using the first candidate grid (lat/lon identical across arrays)
            neighbourGrid = getNeighbourGrid(candidateGrid{1}, linearIdx);
          end

          % Save refined ground coordinate (shared peak)
          refinedLatlonEst(:, estIdx) = candidateGrid{1}.latlonGrid(:, linearIdx);

          % Save refined local DOA for each array (array-dependent mapping)
          for ii = 1:numGridCell
            refinedAngleEst{ii}(:, estIdx) = candidateGrid{ii}.angleGrid(:, linearIdx);
          end
        end

        if numGridCell == 1
          refinedAngleEst = refinedAngleEst{1};
        end
    end

  otherwise
    %% Unsupported dimension
    error('refineGridEstimate:InvalidGridList', ...
      'Only 1D or 2D DOA grids are supported.');
end

end


function linearIdx = pickDominantPeak(Z, shape)
%PICKDOMINANTPEAK Select the dominant peak from a 2D spectrum.
%
% This function identifies local maxima in a 2D objective spectrum using
% regional-maximum detection. Among all detected local peaks, it selects
% the one with the largest spectral value. If no regional maximum exists,
% it falls back to the global maximum of the spectrum.
%
%% Syntax
%   linearIdx = pickDominantPeak(Z, shape)
%% Inputs
%   Z      - 2D real-valued spectrum matrix.
%            Typical examples include:
%              * spatial spectrum in azimuth–elevation grid
%              * objective surface in latitude–longitude grid
%   shape  - Size of the 2D grid used for linear indexing, i.e.,
%            shape = [numRows, numCols].
%            This must be consistent with the reshape convention used
%            before calling this function.
%% Outputs
%   linearIdx - Linear index (column-major, MATLAB convention) of the
%               selected dominant peak in the original vectorized grid.
%% Description
%   The function proceeds in two stages:
%     1) Detect all regional maxima of Z using imregionalmax.
%     2) Among those maxima, select the one with the largest value.
%   If step (1) yields no regional maxima (e.g., flat or monotonic surfaces),
%   the function selects the global maximum of Z instead.
%% Notes
%   - Regional maxima are defined with respect to 8-connected neighborhoods
%     (default behavior of imregionalmax for 2D arrays).
%   - The returned linearIdx is computed using sub2ind(shape, row, col)
%     and therefore follows MATLAB's column-major linearization.
%   - This function assumes that Z does not contain NaN values. If NaNs are
%     present, the behavior of imregionalmax and max may be undefined.
%   - No tie-breaking rule is explicitly defined when multiple peaks share
%     the same maximum value; the first encountered peak is returned.
%% See also
%   imregionalmax, max, sub2ind

%% Detect regional maxima
peakMask = imregionalmax(Z);

% Extract values at regional maxima locations
peakVals = Z(peakMask);

% Row/column indices of all regional maxima
[r, c] = find(peakMask);

%% Fallback: no regional maximum detected
if isempty(peakVals)
  % Use global maximum of the entire spectrum
  [~, linearIdx] = max(Z(:));
  return;
end

%% Select dominant peak among all regional maxima
[~, k] = max(peakVals);

% Convert (row, col) index to linear index
linearIdx = sub2ind(shape, r(k), c(k));

end
