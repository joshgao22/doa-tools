function [globalAngleGrid, eciDirectionGrid] = genEciDirectionGrid(resolution, range)
%GENECIDIRECTIONGRID Generate a global direction grid and its ECI unit vectors.
% Builds a uniform angular grid over a given (azimuth, elevation) domain and
% converts each grid point to a unit direction vector in the ECI frame.
%
%Syntax:
%  [globalAngleGrid, eciDirectionGrid] = genEciDirectionGrid(resolution, range)
%
%Inputs:
%  resolution       : grid size
%                     - scalar N      : N_az = N_el = N
%                     - [N_az N_el]   : number of azimuth / elevation samples
%
%  range            : 2x2 angular domain in radians
%                     [az_min az_max; el_min el_max]
%
%Outputs:
%  globalAngleGrid  : 2x(N_az*N_el), each column is [az; el] (radians)
%
%  eciDirectionGrid : 3x(N_az*N_el), each column is a unit direction vector
%                     [x; y; z] in the ECI frame
%
%Notes:
%  - The angular ranges include both endpoints.
%  - The unit direction vector is parameterized as:
%        u = [cos(el)*cos(az);
%             cos(el)*sin(az);
%             sin(el)]
%  - Grid construction uses:
%        azValue = linspace(az_min, az_max, N_az)
%        elValue = linspace(el_min, el_max, N_el)
%        [azGrid, elGrid] = meshgrid(azValue, elValue)
%    so size(azGrid) = size(elGrid) = (N_el-by-N_az).
%    Columns index azimuth samples; rows index elevation samples.
%  - Column-wise vectorization azGrid(:) / elGrid(:) follows MATLAB
%    column-major order, consistent with subscript size [N_el, N_az].
%
%Example:
%  [globalAngleGrid, eciDirectionGrid] = genEciDirectionGrid( ...
%    [181 91], [-pi pi; -pi/2 pi/2]);
%
%See also: meshgrid, linspace, cos, sin

arguments
  resolution {mustBePositive, mustBeNumeric}
  range (2,2) {mustBeNumeric}
end

% -------------------------------------------------------------------------
% Azimuth / elevation range validation
% -------------------------------------------------------------------------
azRange = range(1,:);
elRange = range(2,:);

% Validate azimuth range
if numel(azRange) ~= 2 || azRange(1) >= azRange(2)
  error('genEciDirectionGrid:InvalidAzRange', ...
    'az_range must be a 2-element vector with az_min < az_max.');
end
if azRange(2) - azRange(1) > 2*pi
  error('genEciDirectionGrid:ExcessiveAzSpan', ...
    'Azimuth span must not exceed 2*pi radians.');
end

% Validate elevation range
if numel(elRange) ~= 2 || elRange(1) >= elRange(2)
  error('genEciDirectionGrid:InvalidElRange', ...
    'el_range must be a 2-element vector with el_min < el_max.');
end
if any(elRange < -pi/2) || any(elRange > pi/2)
  error('genEciDirectionGrid:OutOfBoundsEl', ...
    'Elevation values must lie within [-pi/2, pi/2] radians.');
end

% -------------------------------------------------------------------------
% Parse grid resolution
% -------------------------------------------------------------------------
if isscalar(resolution)
  numGridAz = resolution;
  numGridEl = resolution;
elseif isvector(resolution) && numel(resolution) == 2
  numGridAz = resolution(1);
  numGridEl = resolution(2);
else
  error('genEciDirectionGrid:InvalidGridSize', ...
    'resolution must be a scalar or [N_az, N_el].');
end

% -------------------------------------------------------------------------
% Generate azimuth/elevation mesh and map to ECI unit vectors
% -------------------------------------------------------------------------
azValue = linspace(azRange(1), azRange(2), numGridAz);
elValue = linspace(elRange(1), elRange(2), numGridEl);

% meshgrid convention:
%   size(azGrid) = [N_el, N_az]
[azGrid, elGrid] = meshgrid(azValue, elValue);

xCoordinate = cos(elGrid) .* cos(azGrid);
yCoordinate = cos(elGrid) .* sin(azGrid);
zCoordinate = sin(elGrid);

% Column-oriented outputs
globalAngleGrid  = [azGrid(:)'; elGrid(:)'];
eciDirectionGrid = [xCoordinate(:)'; yCoordinate(:)'; zCoordinate(:)'];
end