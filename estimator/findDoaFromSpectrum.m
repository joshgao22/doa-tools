function doaResult = findDoaFromSpectrum(doaGrid, spectrum, numTargets)
%FINDDOAFROMSPECTRUM Find DOA peaks from a spatial spectrum on a grid.
% Finds the N highest peaks in a 1D or 2D spatial spectrum, returning
% estimated directions, peak indices, and resolution status as a struct.
%
%Syntax:
%   doaResult = findDoaFromSpectrum(doaGrid, spectrum, numTargets)
%
%Inputs:
%   doaGrid    - DOA grid structure with fields:
%               - type       : 'angle' or 'latlon'
%               - dimension  : 1 or 2
%               - resolution : scalar or 2-element vector (2D)
%               - angleGrid  : 1×N (1D) or 2×N (2D) angles (radians)
%               - latlonGrid : 2×N (ground only, degrees)
%
%   spectrum   - 1×N vector of spectrum values (must match number of grid points)
%   numTargets - Number of DOAs to estimate (positive integer)
%
%Output:
%   doaResult  - Result structure with fields:
%               - doaGrid        : input doaGrid
%               - spectrum       : input spectrum
%               - peakIndices    : 1×K linear indices (K=numTargets if resolved)
%               - estimatedDoas  : estimated directions
%               - isResolved     : logical
%               - isDiscrete     : logical (kept for compatibility; set false here)
%
%Notes:
%   - For 2D, spectrum is reshaped as [numEl, numAz] (meshgrid output shape).
%   - If detected peak count < numTargets, peakIndices/estimatedDoas are empty
%     and isResolved=false.
%   - Output indices follow MATLAB column-major order.
%
%Example:
%   res = findDoaFromSpectrum(g, s, 2);
%
%See also: findpeaks, imregionalmax, sub2ind

arguments (Input)
  doaGrid (1,1) struct
  spectrum (1,:) {mustBeNumeric, mustBeVector}
  numTargets (1,1) {mustBeInteger, mustBePositive}
end

% -------------------------------------------------------------------------
% Initialize output structure (DoaResult as struct)
% -------------------------------------------------------------------------
doaResult = struct();
doaResult.doaGrid     = doaGrid;
doaResult.spectrum    = spectrum;
doaResult.peakIndices = [];
doaResult.angleEst    = [];
doaResult.latlonEst   = [];
doaResult.isResolved  = false;
doaResult.isDiscrete  = false;
% -------------------------------------------------------------------------
% Plot payload (directly plottable spectrum)
% -------------------------------------------------------------------------
doaResult.plot = struct();
doaResult.plot.dimension = doaGrid.dimension;
doaResult.plot.x = []; doaResult.plot.y = []; doaResult.plot.z = [];
doaResult.plot.xLabel = ""; doaResult.plot.yLabel = "";
doaResult.plot.gridShape = [];


% -------------------------------------------------------------------------
% Size check
% -------------------------------------------------------------------------
if numel(spectrum) ~= size(doaGrid.angleGrid, 2)
  error('findDoaFromSpectrum:SizeMismatch', ...
    'Length of spectrum must match number of grid points.');
end

switch doaGrid.dimension
  case 1
    [peakValues, peakIndices] = findpeaks(spectrum);

    if numel(peakIndices) < numTargets
      return;
    end

    [~, idxSorted] = sort(peakValues, 'descend');
    topIndices = peakIndices(idxSorted(1:numTargets));
    topIdxSorted = sort(topIndices, 'ascend');

    angleEst = doaGrid.angleGrid(topIdxSorted);

    doaResult.peakIndices = topIdxSorted;
    doaResult.angleEst    = angleEst;
    doaResult.isResolved  = true;

    % 1D plot data
    doaResult.plot.gridShape = size(spectrum);

    switch string(doaGrid.type)
      case 'angle'
        doaResult.plot.x = doaGrid.angleGrid(:).';       % 1xN
        doaResult.plot.xLabel = "Angle (rad)";
        doaResult.plot.y = spectrum(:).';
        doaResult.plot.yLabel = "Spectrum";
      otherwise
        error('findDoaFromSpectrum:UnsupportedType', ...
          '1D DOA grid does not support doaGrid.type = ''ground''.');
    end

  case 2
    % Determine grid shape
    if ~isfield(doaGrid, 'resolution') || isempty(doaGrid.resolution)
      error('findDoaFromSpectrum:MissingProperty', ...
        'doaGrid.resolution is required for 2D grid.');
    end

    if isscalar(doaGrid.resolution)
      numAz = doaGrid.resolution;
      numEl = doaGrid.resolution;
    elseif isvector(doaGrid.resolution) && numel(doaGrid.resolution) == 2
      % Keep the original assumption in your code:
      %   [numAz, numEl] = doaGrid.resolution
      numAz = doaGrid.resolution(1);
      numEl = doaGrid.resolution(2);
    else
      error('findDoaFromSpectrum:InvalidResolution', ...
        'doaGrid.resolution must be scalar or [numAz, numEl].');
    end

    % slow-varying azimuth is assumed
    spectrum2D = reshape(spectrum, numEl, numAz); % [el, az]
    peakMask   = imregionalmax(spectrum2D);

    peakValues = spectrum2D(peakMask);
    [peakElIdx, peakAzIdx] = find(peakMask);

    if numel(peakValues) < numTargets
      return;
    end

    [~, idxSorted] = sort(peakValues, 'descend');
    topElIdx = peakElIdx(idxSorted(1:numTargets));
    topAzIdx = peakAzIdx(idxSorted(1:numTargets));

    linearIndices = sub2ind([numEl, numAz], topElIdx, topAzIdx);
    linearIdxSorted = sort(linearIndices, 'ascend');

    switch string(doaGrid.type)
      case 'angle'
        angleEst = doaGrid.angleGrid(:, linearIdxSorted);
        latlonEst = [];

      case 'latlon'
        angleEst = doaGrid.angleGrid(:, linearIdxSorted);
        latlonEst = doaGrid.latlonGrid(:, linearIdxSorted);

      otherwise
        error('findDoaFromSpectrum:InvalidType', ...
          'doaGrid.type must be ''angle'' or ''latlon''.');
    end

    doaResult.peakIndices = linearIdxSorted;
    doaResult.angleEst    = angleEst;
    doaResult.latlonEst   = latlonEst;
    doaResult.isResolved  = true;

    % 2D plot data (meshgrid + Z)
    doaResult.plot.gridShape = [numEl, numAz];
    doaResult.plot.z = spectrum2D;
    
    switch doaGrid.type
      case 'angle'
        % doaGrid.angleGrid assumed to be 2xN:
        %   row 1: azimuth samples
        %   row 2: elevation samples (or vice versa, but consistent)
        azVec = reshape(doaGrid.angleGrid(1,:), numEl, numAz);
        elVec = reshape(doaGrid.angleGrid(2,:), numEl, numAz);
    
        % Extract representative grid vectors for mesh generation
        % Assumes a regular 2D grid
        az = azVec(1,:);      % 1 x numAz
        el = elVec(:,1);      % numEl x 1
    
        % Generate plotting mesh
        [doaResult.plot.x, doaResult.plot.y] = meshgrid(az, el);
        doaResult.plot.xLabel = "Azimuth (rad)";
        doaResult.plot.yLabel = "Elevation (rad)";
    
      case 'latlon'
        % latlonGrid: 2xN vectorized ground grid
        % Reshaped to [numEl, numAz] for plotting
        latMat = reshape(doaGrid.latlonGrid(1,:), numEl, numAz);
        lonMat = reshape(doaGrid.latlonGrid(2,:), numEl, numAz);
    
        % Directly use latitude/longitude matrices for mesh plotting
        % Note: non-orthogonal grids may introduce geometric distortion
        doaResult.plot.x = lonMat;
        doaResult.plot.y = latMat;
        doaResult.plot.xLabel = "Longitude (deg)";
        doaResult.plot.yLabel = "Latitude (deg)";
    
      otherwise
        error('findDoaFromSpectrum:InvalidType', ...
          'doaGrid.type must be ''local'' or ''ground''.');
    end

  otherwise
    error('findDoaFromSpectrum:UnsupportedDimension', ...
      'Only 1D or 2D DOA grids are supported.');
end
end
