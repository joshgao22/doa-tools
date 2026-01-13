function [latlonGrid, ecefGrid] = default_ground_grid(range, resolution, spheroid)
%GENGROUNDGRID Generate latitude-longitude grid and corresponding ECEF coordinates.
%   [latlonGrid, ecefGrid] = default_ground_grid(range, resolution, spheroid)
%Inputs:
%   range - 2×2 matrix:
%     range(1,:) = [lat_min lat_max] in degrees
%     range(2,:) = [lon_min lon_max] in degrees
%   resolution - scalar or 2-element vector:
%     - scalar: N_lat = N_lon = resolution
%     - [N_lat N_lon]: number of points along latitude / longitude
%   spheroid - reference spheroid object, e.g. wgs84Ellipsoid('meters')
%Outputs:
%   latlonGrid - 2×N matrix, each column is [lat; lon] (degrees)
%   ecefGrid - 3×N matrix, each column is [x; y; z] (meters, ECEF)
%Notes:
%   - All grid points are generated at altitude 0 (sea level).
%   - Both endpoints of [lat_min, lat_max] and [lon_min, lon_max] are included.
%   - Mesh is generated as [latGrid, lonGrid] = meshgrid(lat, lon);
%     rows: latitude, columns: longitude.

arguments
  range (2,2) {mustBeNumeric}
  resolution {mustBeInteger, mustBePositive}
  spheroid {mustBeA(spheroid, ["referenceEllipsoid", "oblateSpheroid", ...
    "referenceSphere"])} = wgs84Ellipsoid("meter")
end

%% parameter check
latRange = range(1,:);
lonRange = range(2,:);

% validate lat_range
if latRange(1) >= latRange(2)
  error('lat_range must be a 2-element vector with lat_min < lat_max.');
end
if any(latRange < -90) || any(latRange > 90)
  error('Latitude values must lie within [-90, 90] degrees.');
end

% validate lon_range
if lonRange(1) >= lonRange(2)
  error('lon_range must be a 2-element vector with lon_min < lon_max.');
end
if lonRange(2) - lonRange(1) > 360
  error('Longitude span must not exceed 360 degrees.');
end

% parse grid size
if isscalar(resolution)
  numGridLat = resolution;
  numGridLon = resolution;
elseif isvector(resolution) && numel(resolution) == 2
  numGridLat = resolution(1);
  numGridLon = resolution(2);
else
  error('resolution must be a scalar or 2-element vector [N_lat, N_lon].');
end

%% generate mesh grid (degrees)
latValue = linspace(latRange(1), latRange(2), numGridLat);
lonValue = linspace(lonRange(1), lonRange(2), numGridLon);
[latGrid, lonGrid] = meshgrid(latValue, lonValue);   % rows: lat, cols: lon

% convert to ECEF (meters, altitude 0)
[x, y, z] = geodetic2ecef(spheroid, latGrid, lonGrid, zeros(size(latGrid)));

latlonGrid = [latGrid(:)'; lonGrid(:)'];
ecefGrid   = [x(:)'; y(:)'; z(:)'];
end
