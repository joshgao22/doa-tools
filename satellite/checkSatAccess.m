function satAccess = checkSatAccess(usrPos, satPos, minUsrElevationDeg, maxSatOffAxisDeg)
%CHECKSATACCESS Evaluate satellite visibility and beam availability.
% Computes user-side elevation, satellite-side nadir off-axis angle, and 
% availability indicators for a set of satellites with respect to one user position.
%
%Syntax:
%  satAccess = checkSatAccess(usrPos, satPos)
%  satAccess = checkSatAccess(usrPos, satPos, minUsrElevationDeg, maxSatOffAxisDeg)
%
%Inputs:
%  usrPos             - 3x1 user position vector
%
%  satPos             - satellite positions, supported formats:
%                       - 3xN numeric matrix
%                       - 3x1xN numeric array
%                       - 1xK / Kx1 cell, each cell is:
%                           * 3xN matrix, or
%                           * 3x1xN array
%
%  minUsrElevationDeg - (optional) minimum user elevation angle in degrees
%                       default: 15
%
%  maxSatOffAxisDeg   - (optional) maximum satellite nadir off-axis angle
%                       in degrees
%                       default: 55
%
%Output:
%  satAccess          - structure with fields:
%                       .isAvailable       : Nsx1 logical
%                       .isVisible         : Nsx1 logical
%                       .isInBeam          : Nsx1 logical
%                       .usrElevationDeg   : Nsx1 user elevation angle
%                       .satOffAxisDeg     : Nsx1 satellite nadir off-axis angle
%                       .slantRange        : Nsx1 slant range
%                       .numSat            : scalar
%
%Notes:
%  - User elevation is defined by the LOS direction relative to the local
%    zenith at the user position.
%  - Satellite off-axis angle is measured between the satellite-to-user
%    direction and the satellite nadir direction.
%  - For cell input, all satellite positions are concatenated column-wise.
%
%See also:
%  vecnorm, asind, acosd

arguments
  usrPos (3,1) {mustBeNumeric, mustBeReal, mustBeFinite}
  satPos
  minUsrElevationDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 15
  maxSatOffAxisDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 55
end

% -------------------------------------------------------------------------
% Basic input check
% -------------------------------------------------------------------------
if minUsrElevationDeg < 0 || minUsrElevationDeg > 90
  error('checkSatAccess:InvalidMinUsrElevation', ...
    'minUsrElevationDeg must lie within [0, 90].');
end

if maxSatOffAxisDeg < 0 || maxSatOffAxisDeg > 90
  error('checkSatAccess:InvalidMaxSatOffAxis', ...
    'maxSatOffAxisDeg must lie within [0, 90].');
end

% -------------------------------------------------------------------------
% Parse satellite positions to 3xNs matrix
% -------------------------------------------------------------------------
satPosMat = localParseSatPos(satPos);
numSat = size(satPosMat, 2);

% -------------------------------------------------------------------------
% User local zenith direction
% -------------------------------------------------------------------------
usrNorm = norm(usrPos);
if usrNorm == 0
  error('checkSatAccess:InvalidUsrPos', ...
    'usrPos must be a non-zero 3x1 vector.');
end
usrZenithHat = usrPos / usrNorm;

% -------------------------------------------------------------------------
% Line-of-sight vectors from user to satellites
% -------------------------------------------------------------------------
losMat = satPosMat - usrPos;                   % 3xNs
slantRange = vecnorm(losMat, 2, 1);            % 1xNs

if any(slantRange == 0)
  error('checkSatAccess:CoincidentPosition', ...
    'At least one satellite position coincides with usrPos.');
end

losHatMat = losMat ./ slantRange;              % 3xNs

% -------------------------------------------------------------------------
% User-side elevation angles
%   sin(elev) = dot(losHat, usrZenithHat)
% -------------------------------------------------------------------------
sinUsrElevation = usrZenithHat.' * losHatMat;  % 1xNs
sinUsrElevation = max(-1, min(1, sinUsrElevation));
usrElevationDeg = asind(sinUsrElevation).';    % Nsx1

% -------------------------------------------------------------------------
% Satellite-side nadir off-axis angles
%   satToUsrHat   = -losHat
%   satNadirHat   = -satPosHat
%   cos(offAxis)  = dot(satToUsrHat, satNadirHat)
%                 = dot(losHat, satPosHat)
% -------------------------------------------------------------------------
satNorm = vecnorm(satPosMat, 2, 1);            % 1xNs
if any(satNorm == 0)
  error('checkSatAccess:InvalidSatPos', ...
    'Satellite positions must be non-zero vectors.');
end

satRadialHatMat = satPosMat ./ satNorm;        % 3xNs
cosSatOffAxis = sum(losHatMat .* satRadialHatMat, 1);
cosSatOffAxis = max(-1, min(1, cosSatOffAxis));
satOffAxisDeg = acosd(cosSatOffAxis).';        % Nsx1

% -------------------------------------------------------------------------
% Visibility and beam constraints
% -------------------------------------------------------------------------
isVisible   = usrElevationDeg >= minUsrElevationDeg;
isInBeam    = satOffAxisDeg <= maxSatOffAxisDeg;
isAvailable = isVisible & isInBeam;

% -------------------------------------------------------------------------
% Pack outputs
% -------------------------------------------------------------------------
satAccess = struct( ...
  'isAvailable',      isAvailable, ...
  'isVisible',        isVisible, ...
  'isInBeam',         isInBeam, ...
  'usrElevationDeg',  usrElevationDeg, ...
  'satOffAxisDeg',    satOffAxisDeg, ...
  'slantRange',       slantRange.', ...
  'numSat',           numSat);
end

function satPosMat = localParseSatPos(satPos)
%LOCALPARSESATPOS Convert supported satellite position inputs to 3xNs matrix.

if isnumeric(satPos)
  satPosMat = localNumericSatPosToMat(satPos);
  return;
end

if iscell(satPos)
  if isempty(satPos)
    satPosMat = zeros(3, 0);
    return;
  end

  satPosCell = cell(size(satPos));
  for ii = 1:numel(satPos)
    if ~isnumeric(satPos{ii}) || ~isreal(satPos{ii})
      error('checkSatAccess:InvalidSatPosCell', ...
        'Each cell element in satPos must be a real numeric array.');
    end
    satPosCell{ii} = localNumericSatPosToMat(satPos{ii});
  end

  satPosMat = cat(2, satPosCell{:});
  return;
end

error('checkSatAccess:UnsupportedSatPosType', ...
  'satPos must be a numeric array or a cell array.');
end

function satPosMat = localNumericSatPosToMat(satPos)
%LOCALNUMERICSATPOSTOMAT Convert numeric satellite positions to 3xNs matrix.

if ~isreal(satPos)
  error('checkSatAccess:InvalidSatPosType', ...
    'satPos must be real-valued.');
end

if ismatrix(satPos)
  if size(satPos, 1) ~= 3
    error('checkSatAccess:InvalidSatPosSize', ...
      'Numeric satPos must have size 3xN or 3x1xN.');
  end
  satPosMat = satPos;
  return;
end

if ndims(satPos) == 3
  if size(satPos, 1) ~= 3 || size(satPos, 2) ~= 1
    error('checkSatAccess:InvalidSatPosPageSize', ...
      'Paged satPos must have size 3x1xN.');
  end
  satPosMat = reshape(satPos, 3, []);
  return;
end

error('checkSatAccess:InvalidSatPosDim', ...
  'Numeric satPos must be 2-D or 3-D.');
end