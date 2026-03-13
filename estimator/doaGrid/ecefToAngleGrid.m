function angleGrid = ecefToAngleGrid(ecefGrid, arrayCenter)
%ECEFTOANGLEGRID Convert ECEF ground points to local azimuth/elevation angles.
% Computes local azimuth and elevation angles (in radians) from ECEF ground
% grid points to a given array/satellite center, using a local coordinate
% frame defined at the array center.
%
%Syntax:
%   angleGrid = ecefToAngleGrid(ecefGrid, arrayCenter)
%
%Inputs:
%   ecefGrid    - 3×N matrix of ground points in ECEF coordinates (meters),
%                 each column corresponds to one ground point
%
%   arrayCenter - 3×1 or 3×N ECEF coordinates of the array/satellite center(s) (meters)
%
%Output:
%   angleGrid   - If one array center: 2×N matrix (azimuth, elevation)
%                 If multiple array centers: 1xN cell array of 2×N matrices
%
%Notes:
%   - Direction vectors are computed from the array center to each ground point.
%   - The local coordinate frame at the array center is defined as:
%       z-axis: pointing toward the Earth's center (−arrayCenter direction)
%       x-axis: orthogonal to z-axis, chosen deterministically
%       y-axis: completes a right-handed coordinate system
%   - Azimuth is measured counterclockwise from the local x-axis in the
%     x–y plane; elevation is measured from the local horizontal plane.
%   - Returned angles are suitable for DOA estimation, beam steering,
%     and satellite coverage analysis.
%
%Example:
%   angleGrid = ecefToAngleGrid(doaGrid.ecefGrid, doaGrid.arrayCenter);
%
%See also:
%   atan2, asin, cross, vecnorm

arguments
  ecefGrid (3,:) {mustBeNumeric}
  arrayCenter (3,:) {mustBeNumeric}
end

numGrid = size(ecefGrid, 2);
numArray = size(arrayCenter, 2);  % Number of array centers (columns)

% If only one array center, return a 2xN matrix; otherwise, return a cell array
if numArray == 1
    angleGrid = zeros(2, numGrid);  % 2×N matrix for a single array center
else
    angleGrid = cell(1, numArray);  % 1×N cell array for multiple array centers
end

% Loop over each array center
for aa = 1:numArray
  % -------------------------------------------------------------------------
  % Direction vectors from array center to ground points (ECEF)
  % -------------------------------------------------------------------------
  dirGlobal = ecefGrid - arrayCenter(:, aa);             % 3×N
  dirGlobal = dirGlobal ./ vecnorm(dirGlobal, 2, 1);      % normalize to unit vectors

  % -------------------------------------------------------------------------
  % Construct local coordinate frame at array center
  %   z-axis points toward Earth center
  % -------------------------------------------------------------------------
  zAxis = -arrayCenter(:, aa) / norm(arrayCenter(:, aa));

  % Choose a temporary axis not collinear with zAxis
  if abs(zAxis(1)) < 0.9
    temp = [1; 0; 0];
  else
    temp = [0; 1; 0];
  end

  xAxis = cross(temp, zAxis);
  xAxis = xAxis / norm(xAxis);
  yAxis = cross(zAxis, xAxis); % right-handed system

  rotationMatrix = [xAxis, yAxis, zAxis]; % 3×3

  % -------------------------------------------------------------------------
  % Project global directions into local frame
  % -------------------------------------------------------------------------
  dirLocal = rotationMatrix' * dirGlobal; % 3×N

  xCoordinate = dirLocal(1, :);
  yCoordinate = dirLocal(2, :);
  zCoordinate = dirLocal(3, :);

  % -------------------------------------------------------------------------
  % Compute azimuth and elevation
  % -------------------------------------------------------------------------
  azimuth   = atan2(yCoordinate, xCoordinate); % [-pi, pi]
  azimuth   = mod(azimuth, 2*pi);              % [0, 2pi)
  elevation = asin(zCoordinate);               % [-pi/2, pi/2]

  % Store results in angleGrid for current array center
  if numArray == 1
      angleGrid(1, :) = azimuth;  % Store azimuth for single array center
      angleGrid(2, :) = elevation;  % Store elevation for single array center
  else
      angleGrid{aa} = [azimuth; elevation];  % Store azimuth and elevation in cell array
  end
end

end
