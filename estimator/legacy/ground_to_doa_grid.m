function angleGrid = ground_to_doa_grid(ecefGrid, arrayCenter)
%GROUNDTOANGLEGRID Convert ECEF ground points to local azimuth/elevation.
%   angleGrid = ground_to_doa_grid(ecefGrid, arrayCenter, dimension)
%Inputs:
%   ecefGrid - 3×N ECEF coordinates of ground points (meters)
%   arrayCenter - 3×1 ECEF coordinate of the satellite/array center (meters)
%Output:
%   angleGrid:
%       if dimension = 1 : 1×N azimuth angles (radian)
%       if dimension = 2 : 2×N [azimuth; elevation]
%Notes:
%   - Local frame:
%     z-axis: pointing toward Earth's center (i.e., –arrayCenter direction)
%     x-axis: orthogonal to z (constructed automatically)
%     y-axis: completes right-handed frame
%   - azimuth ∈ [0, 2π)
%   - elevation ∈ [–π/2, π/2]

arguments
  ecefGrid (3,:) {mustBeNumeric}
  arrayCenter {mustBeNumeric, mustBeVector}
end

%% parameter check
if isempty(arrayCenter) || numel(arrayCenter) ~= 3
  error("ECEF grid requires a 3-element ecefCenter vector (satellite position).");
else
  arrayCenter = arrayCenter(:);
end

%% global direction vectors
dirGlobal = ecefGrid - arrayCenter;           % 3×N
dirGlobal = dirGlobal ./ vecnorm(dirGlobal,2,1);  % Normalize

%% build local coordinate system
zAxis = -arrayCenter / norm(arrayCenter);     % toward Earth center

% choose a stable vector for cross product
if abs(zAxis(1)) < 0.9
  temp = [1;0;0];
else
  temp = [0;1;0];
end

xAxis = cross(temp, zAxis);  xAxis = xAxis / norm(xAxis);
yAxis = cross(zAxis, xAxis);

R = [xAxis, yAxis, zAxis];    % 3×3 (local-to-global)

% project to local frame
dirLocal = R' * dirGlobal;   % 3×N

x = dirLocal(1,:);
y = dirLocal(2,:);
z = dirLocal(3,:);

% compute angles
az = atan2(y, x);       % [-pi, pi]
az = mod(az, 2*pi);     % [0, 2pi)
el = asin(z);

angleGrid = [az; el];
end
