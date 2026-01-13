function [dg_rad, dg_display, dg_range, dg_list] = default_doa_grid(n, unit, dim, dg_range)
%DEFAULT_DOA_GRID Create a uniform DOA search grid.
%
% [dg_rad, dg_display, dg_range, dg_list] = DEFAULT_DOA_GRID(n, unit, dim)
% [dg_rad, dg_display, dg_range, dg_list] = DEFAULT_DOA_GRID(n, unit, dim, dg_range)
%
% This function generates a uniform DOA search grid using half-open
% intervals (the right endpoint is excluded). For 2D DOA, the grid is 
% constructed by forming per-axis coordinate vectors, generating a 
% meshgrid, and flattening the 2D grid for algorithmic use.
%
% Inputs:
%   n        Number of grid points.
%            - If dim == 1: scalar, number of points for 1D DOA.
%            - If dim == 2:
%              · scalar: uses [n n] for azimuth and elevation
%              · [n_az n_el]: number of points per axis.
%   unit     Unit of the grid values:
%            - 'radian': input ranges are given in radians; output grid is
%                        constructed in radians.
%            - 'degree': input ranges are given in degrees; internal
%                        computations use radians; dg_display and dg_list
%                        are in degrees.
%            - 'sin'   : grid is generated in the sine domain:
%                        · 1D: sin(theta) ∈ [-1, 1)
%                        · 2D: [sin(az - pi/2); sin(el)] where az maps to [0, pi].
%   dim      DOA dimension:
%            - 1 : single-angle DOA
%            - 2 : azimuth + elevation DOA
%   dg_range (optional) DOA range.
%            If empty, a default range is used depending on dim and unit.
%            For dim == 1:
%              - 'radian' : [az_min az_max] in radians (default [-pi/2, pi/2))
%              - 'degree' : [az_min az_max] in degrees (default [-90, 90))
%              - 'sin'    : [s_min s_max]    (default [-1, 1))
%            For dim == 2:
%              - 'radian' : [az_min az_max; el_min el_max] in radians
%                           (default [0 2*pi; 0 pi/2))
%              - 'degree' : same format in degrees
%                           (default [0 360; 0 90])
%              - 'sin'    : same format in sine domain
%                           (default [-1 1; 0 1])
%            All ranges are interpreted as half-open intervals.
% Outputs:
%   dg_rad      Grid values in radians.
%               - dim == 1: 1×n vector
%               - dim == 2: 2×(n_az*n_el) matrix
%                   · first row: azimuth (flattened)
%                   · second row: elevation (flattened)
%   dg_display  Grid values in the selected display unit ('radian',
%               'degree', or 'sin'). Same layout as dg_rad.
%   dg_range    Effective range used for grid construction.
%               Returned in the same unit as the user input:
%               radians for 'radian', degrees for 'degree', and sine domain
%               values for 'sin'.
%   dg_list     Per-axis coordinate vectors for meshgrid reconstruction.
%               - dim == 1: 1D display-unit vector
%               - dim == 2: {az_list, el_list}, each a 1D vector in the
%                           display unit ('radian', 'degree', or 'sin').
% Notes:
%   - In 2D mode, dg_rad is flattened for array-processing algorithms
%     (e.g., MUSIC), while dg_list allows easy reshaping for visualization.
%   - In 'degree' mode, user-provided ranges are interpreted as degrees,
%     but all internal computations for dg_rad use radians.
%   - The right endpoint of every interval is excluded to avoid duplicated
%     boundary points for periodic DOA models.

arguments
  n {mustBeInteger, mustBePositive, mustBeVector}
  unit {mustBeMember(unit, ["radian", "degree", "sin"])}
  dim (1,1) {mustBeMember(dim, [1, 2])}
  dg_range {mustBeNumeric} = []
end

% Validate dg_range shape if provided
if ~isempty(dg_range)
  switch dim
    case 1  % 1D grid: two-element vector
      if ~(isvector(dg_range) && numel(dg_range) == 2)
        error('dg_range must be a 2-element vector when dim == 1.');
      end
    case 2  % 2D grid: 2x2 matrix
      if ~isequal(size(dg_range), [2 2])
        error('dg_range must be a 2x2 matrix when dim == 2.');
      end
    otherwise
      error('Unsupported grid dimension dim = %d.', dim);
  end
end

% Normalize n for 2D case
if dim == 2 && isscalar(n)
  n = [n n];
end

switch lower(unit)
  case 'radian'
    if dim == 1
      % 1D, radians
      if isempty(dg_range); dg_range = [-pi/2, pi/2]; end
      dg_step = (dg_range(2) - dg_range(1)) / n;
      dg_display = dg_range(1):dg_step:(dg_range(2) - dg_step);
      dg_rad = dg_display;
      dg_list = dg_display;
    else
      % 2D, radians
      if isempty(dg_range); dg_range = [0, 2*pi; 0, pi/2]; end

      % Azimuth axis (half-open)
      az_min = dg_range(1,1);
      az_max = dg_range(1,2);
      az_step = (az_max - az_min) / n(1);
      az_list = az_min:az_step:(az_max - az_step);
      az_list_rad = az_list;

      % Elevation axis (half-open)
      el_min = dg_range(2,1);
      el_max = dg_range(2,2);
      el_step = (el_max - el_min) / n(2);
      el_list = el_min:el_step:(el_max - el_step);
      el_list_rad = el_list;

      % 2D meshgrid in radians
      [Az_rad, El_rad] = meshgrid(az_list_rad, el_list_rad); % n_el × n_az

      % Flatten for algorithmic use
      dg_rad = [Az_rad(:).'; El_rad(:).'];  % 2 × (n_az*n_el)

      % Display in rad
      dg_display = dg_rad;

      % Axis lists in display unit (radian)
      dg_list = {az_list, el_list};
    end
  case 'degree'
    if dim == 1
      % 1D, degrees (internal radians)
      if isempty(dg_range); dg_range = [-90, 90]; end
      dg_step = (dg_range(2) - dg_range(1)) / n;
      dg_display = dg_range(1):dg_step:(dg_range(2) - dg_step);
      dg_rad = deg2rad(dg_display);     % map to radians
      dg_list = dg_display;
    else
      % 2D, degrees (internal radians)
      if isempty(dg_range); dg_range = [0, 360; 0, 90]; end

      % Axes in radians
      az_min = dg_range(1,1);
      az_max = dg_range(1,2);
      az_step = (az_max - az_min) / n(1);
      az_list = az_min:az_step:(az_max - az_step);
      az_list_rad = deg2rad(az_list);

      el_min = dg_range(2,1);
      el_max = dg_range(2,2);
      el_step = (el_max - el_min) / n(2);
      el_list = el_min:el_step:(el_max - el_step);
      el_list_rad = deg2rad(el_list);

      % 2D meshgrid in radians
      [Az_rad, El_rad] = meshgrid(az_list_rad, el_list_rad);

      % Flatten in radians
      dg_rad = [Az_rad(:).'; El_rad(:).'];      % 2 × (n_az*n_el)

      % Display in degrees
      dg_display = rad2deg(dg_rad);

      % Axis lists in degrees
      dg_list = {az_list, el_list};
    end
  case 'sin'
    if dim == 1
      % 1D, sin(theta)
      if isempty(dg_range); dg_range = [-1, 1]; end
      dg_step = (dg_range(2) - dg_range(1)) / n;
      dg_display = dg_range(1):dg_step:(dg_range(2) - dg_step); % sin-domain
      dg_rad = asin(dg_display);      % map to radians
      dg_list = dg_display;
    else
      % 2D, sin-domain
      if isempty(dg_range); dg_range = [-1, 1; 0, 1]; end

      % Axes in sin-domain
      s_az_min = dg_range(1,1);
      s_az_max = dg_range(1,2);
      s_az_step = (s_az_max - s_az_min) / n(1);
      s_az_list = s_az_min:s_az_step:(s_az_max - s_az_step);

      s_el_min = dg_range(2,1);
      s_el_max = dg_range(2,2);
      s_el_step = (s_el_max - s_el_min) / n(2);
      s_el_list = s_el_min:s_el_step:(s_el_max - s_el_step);

      % 2D meshgrid in sin-domain
      [Saz, Sel] = meshgrid(s_az_list, s_el_list); % n_el × n_az

      % Flatten sin-domain representation
      dg_display = [Saz(:).'; Sel(:).'];          % 2 × (n_az*n_el)

      % Map to radians:
      % dg_display(1,:) = sin(az - pi/2)
      % dg_display(2,:) = sin(el)
      dg_rad = asin(dg_display);
      dg_rad(1,:) = dg_rad(1,:) + pi/2;           % adjust az to [0, pi]

      % Axis lists in display unit (sin-domain)
      az_list = s_az_list;
      el_list = s_el_list;
      dg_list = {az_list, el_list};
    end

  otherwise
    error('Unknown unit ''%s''.', unit);
end
end
