function satAccess = checkSatAccess(usrPos, satPos, minUsrElevationDeg, maxSatOffAxisDeg)
%CHECKSATACCESS Evaluate satellite visibility and beam availability.
% Computes user-side elevation, satellite-side nadir off-axis angle, and
% availability indicators for a set of satellites with respect to one or
% more user positions.
%Syntax:
%  satAccess = checkSatAccess(usrPos, satPos)
%  satAccess = checkSatAccess(usrPos, satPos, minUsrElevationDeg, maxSatOffAxisDeg)
%Inputs:
%  usrPos             - 3xNu user position matrix
%  satPos             - satellite positions, supported formats:
%                       - 3xNs numeric matrix
%                       - 3x1xNs numeric array
%                       - 1xK / Kx1 cell, each cell is:
%                           * 3xNs matrix, or
%                           * 3x1xNs array
%  minUsrElevationDeg - (optional) minimum user elevation angle in degrees
%                       default: 15
%  maxSatOffAxisDeg   - (optional) maximum satellite nadir off-axis angle
%                       in degrees
%                       default: 55
%Output:
%  satAccess          - structure with fields:
%                       .isAvailable       : NsxNu logical
%                       .isVisible         : NsxNu logical
%                       .isInBeam          : NsxNu logical
%                       .usrElevationDeg   : NsxNu user elevation angle
%                       .satOffAxisDeg     : NsxNu satellite nadir off-axis angle
%                       .slantRange        : NsxNu slant range
%                       .numSat            : scalar
%                       .numUser           : scalar
%Definitions:
%  los(:,k,u)         = satPos(:,k) - usrPos(:,u)
%  slantRange(k,u)    = || los(:,k,u) ||
%  usrElevationDeg    = asind(dot(losHat, usrZenithHat))
%  satOffAxisDeg      = acosd(dot(losHat, satRadialHat))
%Notes:
%  - User elevation is defined by the LOS direction relative to the local
%    zenith at the user position.
%  - Satellite off-axis angle is measured between the satellite-to-user
%    direction and the satellite nadir direction.
%  - For cell input, all satellite positions are concatenated column-wise.
%  - If usrPos contains a single user, the output still keeps the user
%    dimension in the form Nsx1.
%See also:
%  vecnorm, asind, acosd

arguments
  usrPos (3,:) {mustBeNumeric, mustBeReal, mustBeFinite}
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
numUser = size(usrPos, 2);

% -------------------------------------------------------------------------
% User local zenith directions
% -------------------------------------------------------------------------
usrNorm = vecnorm(usrPos, 2, 1);
if any(usrNorm == 0)
  error('checkSatAccess:InvalidUsrPos', ...
    'Each user position must be a non-zero 3x1 vector.');
end
usrZenithHat = usrPos ./ usrNorm;

% -------------------------------------------------------------------------
% Line-of-sight vectors from users to satellites
% -------------------------------------------------------------------------
satPosExp = reshape(satPosMat, 3, numSat, 1);           % 3xNsx1
usrPosExp = reshape(usrPos, 3, 1, numUser);             % 3x1xNu
losMat = satPosExp - usrPosExp;                         % 3xNsxNu
slantRange = reshape(vecnorm(losMat, 2, 1), numSat, numUser);

if any(slantRange(:) == 0)
  error('checkSatAccess:CoincidentPosition', ...
    'At least one satellite position coincides with a user position.');
end

losHatMat = losMat ./ reshape(slantRange, 1, numSat, numUser);

% -------------------------------------------------------------------------
% User-side elevation angles
%   sin(elev) = dot(losHat, usrZenithHat)
% -------------------------------------------------------------------------
usrZenithExp = reshape(usrZenithHat, 3, 1, numUser);   % 3x1xNu
sinUsrElevation = reshape(sum(losHatMat .* usrZenithExp, 1), numSat, numUser);
sinUsrElevation = max(-1, min(1, sinUsrElevation));
usrElevationDeg = asind(sinUsrElevation);

% -------------------------------------------------------------------------
% Satellite-side nadir off-axis angles
%   satToUsrHat   = -losHat
%   satNadirHat   = -satPosHat
%   cos(offAxis)  = dot(satToUsrHat, satNadirHat)
%                 = dot(losHat, satPosHat)
% -------------------------------------------------------------------------
satNorm = vecnorm(satPosMat, 2, 1);
if any(satNorm == 0)
  error('checkSatAccess:InvalidSatPos', ...
    'Satellite positions must be non-zero vectors.');
end

satRadialHatMat = satPosMat ./ satNorm;
satRadialExp = reshape(satRadialHatMat, 3, numSat, 1); % 3xNsx1
cosSatOffAxis = reshape(sum(losHatMat .* satRadialExp, 1), numSat, numUser);
cosSatOffAxis = max(-1, min(1, cosSatOffAxis));
satOffAxisDeg = acosd(cosSatOffAxis);

% -------------------------------------------------------------------------
% Visibility and beam constraints
% -------------------------------------------------------------------------
isVisible = usrElevationDeg >= minUsrElevationDeg;
isInBeam = satOffAxisDeg <= maxSatOffAxisDeg;
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
  'slantRange',       slantRange, ...
  'numSat',           numSat, ...
  'numUser',          numUser);
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