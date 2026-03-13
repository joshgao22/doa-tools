function doaGrid = genDoaGrid(type, dimension, resolution, range, arrayCenter, spheroid)
%GENDOAGRID Generate a DOA grid structure for local or ground-based directions.
% Creates a direction-of-arrival (DOA) grid represented as a structure,
% supporting both array-centric (local angular) grids and geodetic
% (ground/ECEF) grids for satellite and remote-sensing applications.
%
%Syntax:
%   doaGrid = genDoaGrid('angle', dimension, resolution)
%   doaGrid = genDoaGrid('angle', dimension, resolution, range)
%   doaGrid = genDoaGrid('latlon', 2, resolution, range, arrayCenter)
%   doaGrid = genDoaGrid('latlon', 2, resolution, range, arrayCenter, spheroid)
%
%Inputs:
%   type        - Grid type: 'angle' or 'latlon'
%                 'angle'  : array-centric angular grid
%                 'latlon' : latitude/longitude grid mapped to ECEF
%
%   dimension   - Grid dimensionality
%                 1: 1D grid (theta, local mode only)
%                 2: 2D grid (azimuth-elevation for local,
%                              latitude-longitude for ground)
%
%   resolution  - Grid resolution
%                 Scalar: number of points per dimension
%                 2-element vector: [N1, N2] for 2D grid
%
%   range       - Grid domain
%                 local, 1D: [theta_min, theta_max] (radians)
%                 local, 2D: [az_min az_max; el_min el_max] (radians)
%                 ground   : [lat_min lat_max; lon_min lon_max] (degrees)
%
%   arrayCenter - (Required for 'latlon')
%                 3×1 ECEF coordinates of array/satellite (meters)
%
%   spheroid    - (Optional, ground only)
%                 Reference spheroid for geodetic2ecef conversion
%                 Default: wgs84Ellipsoid("meter")
%
%Output:
%   doaGrid - Structure with fields:
%     - type        : 'angle' or 'latlon'
%     - dimension   : 1 or 2
%     - resolution  : grid resolution
%     - range       : angular or lat/lon domain
%     - angleGrid   : angle grid
%                     local  : angles in radians
%                     ground : [az; el] (radians) from arrayCenter
%     - latlonGrid  : 2×N [latitude; longitude] (degrees), ground only
%     - ecefGrid    : 3×N [x; y; z] (meters, ECEF), ground only
%     - arrayCenter : ECEF reference position (ground only)
%     - spheroid    : reference spheroid object (ground only)
%
%Notes:
%   - For "angle" mode, only angular grids are generated.
%   - For "latlon" mode, a latitude-longitude grid is constructed and
%     mapped to ECEF; angleGrid gives the local azimuth/elevation from the
%     array center to each grid point.
%   - All grid points are stored column-wise.
%
%Example:
%   g1 = genDoaGrid("angle", 1, 180);
%   g2 = genDoaGrid("angle", 2, [90 45], [0 2*pi; 0 pi/4]);
%   g3 = genDoaGrid("latlon", 2, [100 200], [-60 60; 0 180], [7e6;0;0]);
%
%See also:
%   genAngleGrid, genGroundGrid, ground2angleGrid,
%   geodetic2ecef, wgs84Ellipsoid

arguments
  type char {mustBeMember(type, ["angle", "latlon"])}
  dimension (1,1) {mustBeMember(dimension, [1, 2])}
  resolution (1,:) {mustBePositive, mustBeNumeric}
  range {mustBeNumeric} = []
  arrayCenter (3,:) {mustBeNumeric} = []
  spheroid (1,1) {mustBeA(spheroid, ["referenceEllipsoid", "oblateSpheroid", ...
    "referenceSphere"])} = wgs84Ellipsoid("meter")
end

% -------------------------------------------------------------------------
% Initialize output structure
% -------------------------------------------------------------------------
doaGrid = struct();
doaGrid.type       = type;
doaGrid.dimension  = dimension;
doaGrid.resolution = resolution;

switch type
  case 'angle'
    % Local (array-centric) angular grid

    [doaGrid.angleGrid, doaGrid.range] = genAngleGrid(dimension, ...
        resolution, range);
    doaGrid.latlonGrid  = [];
    doaGrid.ecefGrid    = [];
    doaGrid.arrayCenter = [];
    doaGrid.spheroid    = [];

  case 'latlon'
    % Ground-based grid (latitude/longitude mapped to ECEF)

    % ---- Input validation ----
    if isempty(arrayCenter)
      error("genDoaGrid:MissingECEF", ...
        "Ground grid requires arrayCenter (ECEF position).");
    end

    if ~isequal(size(range), [2 2])
      error("genDoaGrid:InvalidRange", ...
        "Ground grid requires range: [lat_min lat_max; lon_min lon_max].");
    end

    doaGrid.range       = range;
    doaGrid.arrayCenter = arrayCenter;
    doaGrid.spheroid    = spheroid;

    % ---- Grid generation ----
    % Each column corresponds to one ground point
    [doaGrid.latlonGrid, doaGrid.ecefGrid] = ...
      genEcefGrid(resolution, range, spheroid);

    % ---- Angle mapping ----
    % Convert ECEF grid to local azimuth/elevation w.r.t. arrayCenter
    doaGrid.angleGrid = ...
      ecefToAngleGrid(doaGrid.ecefGrid, arrayCenter);
end
end
