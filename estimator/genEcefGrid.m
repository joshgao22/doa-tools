function [latlonGrid, ecefGrid] = genEcefGrid(resolution, range, spheroid)
%GENECEFGRID Generate a latitude-longitude grid and its ECEF coordinates.
% Builds a uniform geodetic grid over a given (lat, lon) domain and converts
% each grid point to ECEF coordinates (meters) at altitude 0 using the given
% reference spheroid.
%
%Syntax:
%  [latlonGrid, ecefGrid] = genEcefGrid(resolution, range)
%  [latlonGrid, ecefGrid] = genEcefGrid(resolution, range, spheroid)
%
%Inputs:
%  resolution : grid size
%               - scalar N      : N_lat = N_lon = N
%               - [N_lat N_lon] : number of latitude / longitude samples
%  range      : 2x2 domain in degrees
%               [lat_min lat_max; lon_min lon_max]
%  spheroid   : (optional) reference spheroid for geodetic2ecef
%               default: wgs84Ellipsoid("meter")
%
%Outputs:
%  latlonGrid : 2x(N_lat*N_lon), each column is [lat; lon] (degrees)
%  ecefGrid   : 3x(N_lat*N_lon), each column is [x; y; z]   (meters, ECEF)
%
%Notes:
%  - Altitude is fixed to 0 (sea level).
%  - Latitude/longitude ranges include both endpoints.
%  - Grid construction uses:
%        latValue = linspace(lat_min, lat_max, N_lat)
%        lonValue = linspace(lon_min, lon_max, N_lon)
%        [latGrid, lonGrid] = meshgrid(latValue, lonValue)
%    so size(latGrid) = size(lonGrid) = (N_lon-by-N_lat).
%    Columns index latitude samples; rows index longitude samples.
%  - Column-wise vectorization latGrid(:) / lonGrid(:) follows MATLAB
%    column-major order, consistent with subscript size [N_lon, N_lat].
%
%Example:
%  sph = wgs84Ellipsoid("meter");
%  [latlonGrid, ecefGrid] = genEcefGrid([100 200], [-60 60; 0 180], sph);
%
%See also: geodetic2ecef, wgs84Ellipsoid, meshgrid, linspace

arguments
  resolution {mustBePositive, mustBeNumeric}
  range (2,2) {mustBeNumeric}
  spheroid {mustBeA(spheroid, ["referenceEllipsoid", "oblateSpheroid", ...
    "referenceSphere"])} = wgs84Ellipsoid("meter")
end

% -------------------------------------------------------------------------
% Latitude / longitude range validation
% -------------------------------------------------------------------------
latRange = range(1,:);
lonRange = range(2,:);

% Validate latitude range
if numel(latRange) ~= 2 || latRange(1) >= latRange(2)
  error('genGroundGrid:InvalidLatRange', ...
    'lat_range must be a 2-element vector with lat_min < lat_max.');
end
if any(latRange < -90) || any(latRange > 90)
  error('genGroundGrid:OutOfBoundsLat', ...
    'Latitude values must lie within [-90, 90] degrees.');
end

% Validate longitude range
if numel(lonRange) ~= 2 || lonRange(1) >= lonRange(2)
  error('genGroundGrid:InvalidLonRange', ...
    'lon_range must be a 2-element vector with lon_min < lon_max.');
end
if lonRange(2) - lonRange(1) > 360
  error('genGroundGrid:ExcessiveLonSpan', ...
    'Longitude span must not exceed 360 degrees.');
end

% -------------------------------------------------------------------------
% Parse grid resolution
% -------------------------------------------------------------------------
if isscalar(resolution)
  numGridLat = resolution;
  numGridLon = resolution;
elseif isvector(resolution) && numel(resolution) == 2
  numGridLat = resolution(1);
  numGridLon = resolution(2);
else
  error('genGroundGrid:InvalidGridSize', ...
    'resolution must be a scalar or [N_lat, N_lon].');
end

% -------------------------------------------------------------------------
% Generate latitude/longitude mesh and map to ECEF
% -------------------------------------------------------------------------
latValue = linspace(latRange(1), latRange(2), numGridLat);
lonValue = linspace(lonRange(1), lonRange(2), numGridLon);

% meshgrid convention:
%   size(latGrid) = [N_lon, N_lat]
[latGrid, lonGrid] = meshgrid(latValue, lonValue);

% Altitude fixed at sea level
[x, y, z] = geodetic2ecef(spheroid, latGrid, lonGrid, zeros(size(latGrid)));

% Column-oriented outputs
latlonGrid = [latGrid(:)'; lonGrid(:)'];
ecefGrid   = [x(:)'; y(:)'; z(:)'];
end
