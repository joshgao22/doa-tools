function doaGrid = genDoaGrid( ...
  type, dimension, resolution, range, globalFrame, utc, arrayCenter, rotMat, spheroid)
%GENDOAGRID Generate a DOA grid structure for local or ground-based search.
% Creates a direction-of-arrival (DOA) grid represented as a structure,
% supporting both array-centric angular grids and latitude-longitude grids
% under either the ECEF or ECI global frame.
%
%Syntax:
%  doaGrid = genDoaGrid('angle', dimension, resolution)
%  doaGrid = genDoaGrid('angle', dimension, resolution, range)
%  doaGrid = genDoaGrid('latlon', 2, resolution, range, ...
%    'ecef', [], arrayCenter, [], spheroid)
%  doaGrid = genDoaGrid('latlon', 2, resolution, range, ...
%    'eci', utc, arrayCenter, rotMat, spheroid)
%
%Inputs:
%  type         - Grid type
%                 'angle'  : array-centric angular grid
%                 'latlon' : latitude/longitude grid
%
%  dimension    - Grid dimensionality
%                 1: 1D angular grid
%                 2: 2D angular or lat/lon grid
%
%  resolution   - Grid resolution
%                 scalar          : same number of points per dimension
%                 2-element vector: [N1, N2] for 2D grid
%
%  range        - Grid domain
%                 angle, 1D : [thetaMin, thetaMax] (rad)
%                 angle, 2D : [azMin azMax; elMin elMax] (rad)
%                 latlon    : [latMin latMax; lonMin lonMax] (deg)
%
%  globalFrame  - Global frame for 'latlon' mode
%                 'ecef' : map lat/lon grid to ECEF, then to local DOA
%                 'eci'  : map lat/lon grid to ECI at utc, then to local DOA
%                 For 'angle' mode, this input is ignored and can be [].
%
%  utc          - UTC time for 'latlon' + 'eci' mode
%                 Supported forms:
%                 - datetime scalar
%                 - 1x6 numeric date vector [Y M D h m s]
%                 Otherwise use [].
%
%  arrayCenter  - Reference array position
%                 'latlon' + 'ecef' : 3x1 ECEF array center
%                 'latlon' + 'eci'  : 3x1 ECI  array center
%
%  rotMat       - Local-to-global rotation matrix for 'latlon' + 'eci'
%                 3x3 matrix, mapping local coordinates to ECI
%
%  spheroid     - Reference spheroid for geodetic conversion
%                 Default: wgs84Ellipsoid("meter")
%
%Output:
%  doaGrid      - Structure with fields:
%                 .type
%                 .dimension
%                 .resolution
%                 .range
%                 .globalFrame
%                 .angleGrid
%                 .latlonGrid
%                 .ecefGrid
%                 .eciGrid
%                 .utc
%                 .arrayCenter
%                 .rotMat
%                 .spheroid
%
%Notes:
%  - For 'angle' mode, only angleGrid is generated.
%  - For 'latlon' + 'ecef', angleGrid is computed by ecefToAngleGrid.
%  - For 'latlon' + 'eci', angleGrid is computed by globalToLocalDoa after
%    converting each lat/lon grid point to ECI at the given utc.
%  - All grid points are stored column-wise.
%
%See also:
%  genAngleGrid, genEcefGrid, ecefToAngleGrid, globalToLocalDoa, lla2eci

arguments
  type {mustBeMember(type, ["angle", "latlon"])}
  dimension (1,1) {mustBeMember(dimension, [1, 2])}
  resolution (1,:) {mustBePositive, mustBeNumeric}
  range {mustBeNumeric} = []
  globalFrame {mustBeMember(globalFrame, ["ecef", "eci"])} = []
  utc (:,6) {mustBeNumeric} = []
  arrayCenter (3,1) {mustBeNumeric} = zeros(3,1)
  rotMat {mustBeNumeric} = []
  spheroid (1,1) {mustBeA(spheroid, ["referenceEllipsoid", ...
    "oblateSpheroid", "referenceSphere"])} = wgs84Ellipsoid("meter")
end

% -------------------------------------------------------------------------
% Initialize output structure
% -------------------------------------------------------------------------
doaGrid = struct();
doaGrid.type        = char(string(type));
doaGrid.dimension   = dimension;
doaGrid.resolution  = resolution;
doaGrid.range       = range;
doaGrid.globalFrame = [];
doaGrid.angleGrid   = [];
doaGrid.latlonGrid  = [];
doaGrid.ecefGrid    = [];
doaGrid.eciGrid     = [];
doaGrid.utc         = [];
doaGrid.arrayCenter = [];
doaGrid.rotMat      = [];
doaGrid.spheroid    = [];

switch lower(string(type))
  case "angle"
    % ---------------------------------------------------------------------
    % Local angular grid
    % ---------------------------------------------------------------------
    [doaGrid.angleGrid, doaGrid.range] = genAngleGrid(dimension, ...
      resolution, range);

  case "latlon"
    % ---------------------------------------------------------------------
    % Latitude / longitude grid
    % ---------------------------------------------------------------------
    if dimension ~= 2
      error('genDoaGrid:InvalidDimension', ...
        '''latlon'' grid only supports dimension = 2.');
    end

    if ~isequal(size(range), [2, 2])
      error('genDoaGrid:InvalidRange', ...
        ['''latlon'' grid requires range = ', ...
         '[latMin latMax; lonMin lonMax].']);
    end

    if isempty(globalFrame)
      error('genDoaGrid:MissingGlobalFrame', ...
        ['''latlon'' grid requires globalFrame = ', ...
         '''ecef'' or ''eci''.']);
    end

    globalFrame = lower(string(globalFrame));

    doaGrid.globalFrame = char(globalFrame);
    doaGrid.arrayCenter = arrayCenter;
    doaGrid.spheroid    = spheroid;

    % ---------------------------------------------------------------------
    % Generate latitude / longitude grid and ECEF positions
    % ---------------------------------------------------------------------
    [doaGrid.latlonGrid, doaGrid.ecefGrid] = genEcefGrid( ...
      resolution, range, spheroid);

    switch globalFrame
      case "ecef"
        % ---------------------------------------------------------------
        % ECEF mode
        % ---------------------------------------------------------------
        if all(arrayCenter == 0)
          error('genDoaGrid:MissingArrayCenter', ...
            ['''latlon'' + ''ecef'' grid requires non-empty ', ...
             'arrayCenter in ECEF.']);
        end

        doaGrid.angleGrid = ecefToAngleGrid(doaGrid.ecefGrid, arrayCenter);

      case "eci"
        % ---------------------------------------------------------------
        % ECI mode
        % ---------------------------------------------------------------
        if isempty(utc)
          error('genDoaGrid:MissingUtc', ...
            '''latlon'' + ''eci'' grid requires utc.');
        end

        if all(arrayCenter == 0)
          error('genDoaGrid:MissingArrayCenter', ...
            ['''latlon'' + ''eci'' grid requires non-empty ', ...
             'arrayCenter in ECI.']);
        end

        if isempty(rotMat) || ~isequal(size(rotMat), [3, 3])
          error('genDoaGrid:InvalidRotMat', ...
            ['''latlon'' + ''eci'' grid requires rotMat to be a ', ...
             'non-empty 3x3 matrix.']);
        end

        doaGrid.utc    = utc;
        doaGrid.rotMat = rotMat;
        doaGrid.eciGrid = localLatlonToEciGrid(doaGrid.latlonGrid, utc);

        doaGrid.angleGrid = globalToLocalDoa( ...
          doaGrid.eciGrid, arrayCenter, rotMat);
    end
end

end

% =========================================================================
% Local functions
% =========================================================================
function eciGrid = localLatlonToEciGrid(latlonGrid, utc)
%LOCALLATLONTOECIGRID Convert [lat; lon] grid to ECI positions at utc.

numGrid = size(latlonGrid, 2);
eciGrid = zeros(3, numGrid);

if isa(utc, 'datetime')
  utc = datevec(utc);
end
if isnumeric(utc) && isequal(size(utc), [1, 6])
  utcMat = repmat(utc, numGrid, 1);
elseif isnumeric(utc) && size(utc, 2) == 6 && size(utc, 1) == numGrid
  utcMat = utc;
else
  error('genDoaGrid:InvalidUtcFormat', ...
    'utc must be a datetime scalar, 1x6 date vector, or N-by-6 date matrix.');
end

llaTmp = [latlonGrid.', zeros(numGrid, 1)];
eciGrid = lla2eci(llaTmp, utcMat).';

end
