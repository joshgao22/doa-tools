function [refinedAngleEst, refinedLatlonEst] = refineGridEstimate( ...
  objFunction, doaGrid, doaResult, maxIterations, subgridSize)
%REFINEGRIDESTIMATE Refine coarse DOA peaks by iterative local grid search.
% Starting from coarse peak indices, this function repeatedly builds a
% smaller local subgrid around each peak, evaluates the objective spectrum,
% and updates the search window toward the dominant local peak.
%
%Syntax:
%  refinedAngleEst = refineGridEstimate(objFunction, doaGrid, doaResult)
%  [refinedAngleEst, refinedLatlonEst] = refineGridEstimate( ...
%    objFunction, doaGrid, doaResult, maxIterations, subgridSize)
%
%Inputs:
%  objFunction    - objective function handle
%                   Signature:
%                     objSpectrum = objFunction(candidateGrid)
%                   where candidateGrid is either:
%                     - a single grid struct, or
%                     - a cell array of grid structs
%
%  doaGrid        - coarse grid definition
%                   - struct      : single-grid refinement
%                   - 1-by-N cell : multi-grid refinement (latlon only)
%
%  doaResult      - coarse estimation result, must contain:
%                   .peakIndices
%
%  maxIterations  - number of refinement iterations
%                   default: 10
%
%  subgridSize    - local refinement grid resolution
%                   - 1D : scalar
%                   - 2D : scalar or 2-element vector
%                   default: 10
%
%Outputs:
%  refinedAngleEst  - refined local DOA estimates
%                     - angle grid:
%                       * 1D      : 1-by-K
%                       * 2D      : 2-by-K
%                     - latlon grid:
%                       * single  : 2-by-K
%                       * multi   : 1-by-N cell, each 2-by-K
%
%  refinedLatlonEst - refined [lat; lon] estimates for latlon grid
%                     empty for angle grid
%
%Notes:
%  - For latlon refinement, the candidate subgrid is reconstructed using
%    the same global-frame metadata stored in each doaGrid cell:
%      globalFrame, utc, arrayCenter, rotMat, spheroid
%  - This allows both ECEF and ECI latlon grids to be refined through the
%    same interface.
%
%See also:
%  genDoaGrid, getNeighbourGrid, findDoaFromSpectrum, imregionalmax

arguments
  objFunction (1,1) function_handle
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  doaResult (1,1) struct
  maxIterations (1,1) {mustBePositive, mustBeInteger} = 10
  subgridSize {mustBePositive, mustBeInteger} = 10
end

% -------------------------------------------------------------------------
% Normalize grid input to cell container
% -------------------------------------------------------------------------
if ~iscell(doaGrid)
  doaGrid = {doaGrid};
end

numGridCell = numel(doaGrid);
baseGrid = doaGrid{1};

if ~isfield(baseGrid, 'dimension') || isempty(baseGrid.dimension)
  error('refineGridEstimate:MissingProperty', ...
    'doaGrid.dimension is required.');
end
if ~isfield(baseGrid, 'type') || isempty(baseGrid.type)
  error('refineGridEstimate:MissingProperty', ...
    'doaGrid.type is required.');
end
if ~isfield(doaResult, 'peakIndices') || isempty(doaResult.peakIndices)
  error('refineGridEstimate:MissingProperty', ...
    'doaResult.peakIndices is required.');
end

dim = baseGrid.dimension;
gridType = lower(string(baseGrid.type));

if gridType == "angle" && numGridCell ~= 1
  error('refineGridEstimate:InvalidGridList', ...
    'Multiple grid cells are only supported for ''latlon'' mode.');
end

% -------------------------------------------------------------------------
% Initialize outputs
% -------------------------------------------------------------------------
refinedAngleEst = [];
refinedLatlonEst = [];

peakIndices = doaResult.peakIndices(:);
numEst = numel(peakIndices);

% -------------------------------------------------------------------------
% Refinement by grid dimension
% -------------------------------------------------------------------------
switch dim
  case 1
    % =====================================================================
    % 1D refinement: angle grid only
    % =====================================================================
    if gridType ~= "angle"
      error('refineGridEstimate:InvalidGrid', ...
        '1D refinement only supports doaGrid.type = ''angle''.');
    end

    if ~isscalar(subgridSize)
      error('refineGridEstimate:InputDimensionNotMatch', ...
        'subgridSize must be scalar for 1D grid.');
    end

    refinedAngleEst = zeros(1, numEst);

    for estIdx = 1:numEst
      neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

      for iter = 1:maxIterations
        candidateGrid = genDoaGrid('angle', 1, subgridSize, neighbourGrid);
        objSpectrum = objFunction({candidateGrid});
        [~, peakIdx] = max(objSpectrum);
        neighbourGrid = getNeighbourGrid(candidateGrid, peakIdx);
      end

      refinedAngleEst(estIdx) = candidateGrid.angleGrid(peakIdx);
    end

  case 2
    switch gridType
      case "angle"
        % =================================================================
        % 2D refinement: local angle grid
        % =================================================================
        if isscalar(subgridSize)
          numAz = subgridSize;
          numEl = subgridSize;
        elseif isvector(subgridSize) && numel(subgridSize) == 2
          numAz = subgridSize(1);
          numEl = subgridSize(2);
        else
          error('refineGridEstimate:InputDimensionNotMatch', ...
            'subgridSize must be scalar or 2-element vector for 2D angle grid.');
        end

        refinedAngleEst = zeros(2, numEst);

        for estIdx = 1:numEst
          neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

          for iter = 1:maxIterations
            candidateGrid = genDoaGrid('angle', 2, subgridSize, neighbourGrid);
            objSpectrum = objFunction({candidateGrid});

            % angle grid ordering: reshape to [numEl, numAz]
            objSpectrum2D = reshape(objSpectrum, numEl, numAz);
            linearIdx = pickDominantPeak(objSpectrum2D, [numEl, numAz]);

            neighbourGrid = getNeighbourGrid(candidateGrid, linearIdx);
          end

          refinedAngleEst(:, estIdx) = candidateGrid.angleGrid(:, linearIdx);
        end

      case "latlon"
        % =================================================================
        % 2D refinement: latlon grid (single or multiple arrays)
        % =================================================================
        if isscalar(subgridSize)
          numLat = subgridSize;
          numLon = subgridSize;
        elseif isvector(subgridSize) && numel(subgridSize) == 2
          numLat = subgridSize(1);
          numLon = subgridSize(2);
        else
          error('refineGridEstimate:InputDimensionNotMatch', ...
            'subgridSize must be scalar or 2-element vector for 2D latlon grid.');
        end

        refinedLatlonEst = zeros(2, numEst);
        refinedAngleEst = cell(1, numGridCell);
        for ii = 1:numGridCell
          refinedAngleEst{ii} = zeros(2, numEst);
        end

        for estIdx = 1:numEst
          % Shared neighbour bounds in [lat; lon]
          neighbourGrid = getNeighbourGrid(baseGrid, peakIndices(estIdx));

          for iter = 1:maxIterations
            candidateGrid = cell(size(doaGrid));

            for ii = 1:numGridCell
              candidateGrid{ii} = genDoaGrid( ...
                'latlon', 2, subgridSize, neighbourGrid, ...
                doaGrid{ii}.globalFrame, doaGrid{ii}.utc, ...
                doaGrid{ii}.arrayCenter, doaGrid{ii}.rotMat, ...
                doaGrid{ii}.spheroid);
            end

            objSpectrum = objFunction(candidateGrid);

            % latlon grid ordering: reshape to [numLon, numLat]
            objSpectrum2D = reshape(objSpectrum, numLon, numLat);
            linearIdx = pickDominantPeak(objSpectrum2D, [numLon, numLat]);

            % All candidate grids share the same latlon lattice
            neighbourGrid = getNeighbourGrid(candidateGrid{1}, linearIdx);
          end

          refinedLatlonEst(:, estIdx) = candidateGrid{1}.latlonGrid(:, linearIdx);

          for ii = 1:numGridCell
            refinedAngleEst{ii}(:, estIdx) = candidateGrid{ii}.angleGrid(:, linearIdx);
          end
        end

        if numGridCell == 1
          refinedAngleEst = refinedAngleEst{1};
        end

      otherwise
        error('refineGridEstimate:InvalidGrid', ...
          'Unsupported doaGrid.type for 2D refinement.');
    end

  otherwise
    error('refineGridEstimate:InvalidGrid', ...
      'Only 1D or 2D DOA grids are supported.');
end

end


function linearIdx = pickDominantPeak(objSpectrum2D, gridShape)
%PICKDOMINANTPEAK Pick the dominant peak from a 2D objective surface.
% Regional maxima are preferred. If none exists, the global maximum is used.

peakMask = imregionalmax(objSpectrum2D);
peakVals = objSpectrum2D(peakMask);
[rowIdx, colIdx] = find(peakMask);

if isempty(peakVals)
  [~, linearIdx] = max(objSpectrum2D(:));
  return;
end

[~, bestIdx] = max(peakVals);
linearIdx = sub2ind(gridShape, rowIdx(bestIdx), colIdx(bestIdx));

end
