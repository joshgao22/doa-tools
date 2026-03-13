function [angleGrid, range] = genAngleGrid(dimension, resolution, range)
%GENANGLEGRID Uniform angular grid for array-local DOA scanning (az, el).
%
%Syntax:
%  angleGrid = genAngleGrid(dimension, resolution)
%  angleGrid = genAngleGrid(dimension, resolution, range)
%  [angleGrid, range] = genAngleGrid(...)
%
%Inputs:
%  dimension  : 1 (theta) or 2 (azimuth-elevation)
%  resolution : 1D -> N (integer)
%               2D -> N or [Naz Nel] (integers)
%  range      : (optional, radians)
%               1D -> [thetaMin thetaMax], default [-pi/2 pi/2]
%               2D -> [azMin azMax; elMin elMax], default [0 2*pi; 0 pi/2]
%
%Outputs:
%  angleGrid  : 1D -> 1xN
%               2D -> 2x(Naz*Nel), each column is [az; el] (radians)
%  range      : validated scan range (radians)
%
%Notes:
%  - 2D grid uses meshgrid(az, el): AZ/EL size is (Nel-by-Naz).
%  - Linear indexing (:) corresponds to subscript size [Nel, Naz] (row=el, col=az).
%  - If azimuth range is exactly [0, 2*pi], endpoint 2*pi is excluded (covers [0, 2*pi)).

arguments
  dimension (1,1) {mustBeMember(dimension,[1 2])}
  resolution {mustBeNumeric, mustBePositive}
  range double = []
end

% ---- enforce integer resolution (avoid silent rounding) ----
if any(resolution(:) ~= round(resolution(:)))
  error('genAngleGrid:InvalidResolution', ...
    'resolution must be integer-valued (1D: N, 2D: N or [Naz Nel]).');
end

switch dimension
  case 1
    if ~isscalar(resolution)
      error('genAngleGrid:InvalidResolution', 'For 1D grid, resolution must be a scalar N.');
    end
    N = double(resolution);

    thetaMin = -pi/2; thetaMax = pi/2;

    if isempty(range)
      range = [thetaMin thetaMax];
    else
      valid = isvector(range) && numel(range) == 2 && ...
              range(1) < range(2) && ...
              range(1) >= thetaMin && range(2) <= thetaMax;
      if ~valid
        error('genAngleGrid:InvalidRange', ...
          '1D range must be [thetaMin thetaMax] within [%g, %g] rad.', ...
          thetaMin, thetaMax);
      end
      range = double(range(:)).'; % row
    end

    angleGrid = linspace(range(1), range(2), N);

  case 2
    if isscalar(resolution)
      Naz = double(resolution);
      Nel = double(resolution);
    elseif isvector(resolution) && numel(resolution) == 2
      Naz = double(resolution(1));  % az samples
      Nel = double(resolution(2));  % el samples
    else
      error('genAngleGrid:InvalidResolution', ...
        'For 2D grid, resolution must be N or [Naz Nel].');
    end

    azMin = 0; azMax = 2*pi;
    elMin = 0; elMax = pi/2;

    if isempty(range)
      range = [azMin azMax; elMin elMax];
    else
      valid = isequal(size(range),[2 2]) && ...
              range(1,1) < range(1,2) && range(2,1) < range(2,2) && ...
              range(1,1) >= azMin && range(1,2) <= azMax && ...
              range(2,1) >= elMin && range(2,2) <= elMax;
      if ~valid
        error('genAngleGrid:InvalidRange', ...
          ['2D range must be [azMin azMax; elMin elMax] (rad), ', ...
           'az∈[%g,%g], el∈[%g,%g], with min < max.'], ...
          azMin, azMax, elMin, elMax);
      end
      range = reshape(double(range), 2, 2);
    end

    azR = range(1,:); 
    elR = range(2,:);

    if isequal(azR, [0 2*pi])
      az = linspace(0, 2*pi, Naz + 1);
      az = az(1:end-1);  % avoid duplicate 2*pi
    else
      az = linspace(azR(1), azR(2), Naz);
    end
    el = linspace(elR(1), elR(2), Nel);

    [AZ, EL] = meshgrid(az, el);      % size: Nel-by-Naz
    angleGrid = [AZ(:)'; EL(:)'];     % 2x(Naz*Nel), columns [az; el]
end
end
