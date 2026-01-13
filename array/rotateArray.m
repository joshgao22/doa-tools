function [arrayRot, rotMat] = rotateArray(array, tgtN, rotCtr)
%ROTATEARRAY Rotate an array to align its normal vector to a target direction.
%
%   [arrayRot, rotMat] = rotateArray(array, tgtN)
%   [arrayRot, rotMat] = rotateArray(array, tgtN, rotCtr)
%
% Rotates the array described by struct `array` to align its normal vector
% (assumed along +z for planar/embedded arrays) to the target direction `tgtN`.
% The rotation is performed by translating the array centroid to the origin,
% applying rotation, then translating to `rotCtr`.
%
%Inputs:
%   array  - Array struct with fields:
%            positions : [dim, M] element positions
%            dim       : 1 / 2 / 3
%            count     : M
%   tgtN   - Target normal vector:
%            2D: [nx; ny] or [nx ny]
%            3D: [nx; ny; nz] or [nx ny nz]
%   rotCtr - (Optional) rotation center (same dimension as tgtN), default zeros
%
%Outputs:
%   arrayRot - Rotated array struct with updated positions/dim
%   rotMat   - Rotation matrix (2×2 for 2D, 3×3 for 3D)
%
%Notes:
%   - For dim==1 and tgtN is 2D: performs 2D rotation and embeds the 1D array
%     on the y-axis before rotation (consistent with the original code).
%   - For dim==1 and tgtN is 3D: embeds the 1D array on the y-axis in 3D then rotates.
%   - For dim==2: embedded into 3D (z=0) then rotated.
%   - For dim==3: rotated directly in 3D.
%
%See also: cross, dot, acos

arguments
  array (1,1) struct
  tgtN {mustBeNumeric, mustBeVector}
  rotCtr {mustBeNumeric, mustBeVector} = []
end

% -------------------------------------------------------------------------
% Required field checks
% -------------------------------------------------------------------------
if ~isfield(array, 'positions') || isempty(array.positions)
  error('rotateArray:MissingProperty', 'array.positions is required.');
end
if ~isfield(array, 'dim') || isempty(array.dim)
  error('rotateArray:MissingProperty', 'array.dim is required.');
end
if ~isfield(array, 'count') || isempty(array.count)
  error('rotateArray:MissingProperty', 'array.count is required.');
end

% -------------------------------------------------------------------------
% Normalize tgtN and handle rotCtr
% -------------------------------------------------------------------------
tgtN = tgtN(:);
if ~(numel(tgtN) == 2 || numel(tgtN) == 3)
  error('rotateArray:InvalidTargetNormal', ...
    'tgtN must be a 2D or 3D vector.');
end
if norm(tgtN) < eps
  error('rotateArray:InvalidTargetNormal', 'tgtN must be non-zero.');
end
tgtN = tgtN / norm(tgtN);

if isempty(rotCtr)
  rotCtr = zeros(size(tgtN));
else
  rotCtr = rotCtr(:);
  if numel(rotCtr) ~= numel(tgtN)
    error('rotateArray:CenterDimMismatch', ...
      'rotCtr and tgtN must have the same dimension.');
  end
end

% -------------------------------------------------------------------------
% Dispatch by array.dim
% -------------------------------------------------------------------------
switch array.dim
  case 1
    if numel(tgtN) == 2
      [arrayRot, rotMat] = rotate1dTo2d(array, tgtN, rotCtr);
      arrayRot.dim = 2;
    else
      [arrayRot, rotMat] = rotateTo3d(array, tgtN, rotCtr);
      arrayRot.dim = 3;
    end

  case 2
    if numel(tgtN) ~= 3
      error('rotateArray:InvalidTargetNormal', ...
        'Only 3D rotation is allowed for 2D array.');
    end
    [arrayRot, rotMat] = rotateTo3d(array, tgtN, rotCtr);
    arrayRot.dim = 3;

  case 3
    if numel(tgtN) ~= 3
      error('rotateArray:InvalidTargetNormal', ...
        'Only 3D rotation is allowed for 3D array.');
    end
    [arrayRot, rotMat] = rotateTo3d(array, tgtN, rotCtr);
    arrayRot.dim = 3;

  otherwise
    error('rotateArray:UnsupportedDim', ...
      'array.dim must be 1, 2, or 3.');
end

end


% =====================================================================
% Helpers
% =====================================================================

function [arrayRot, rotMat] = rotate1dTo2d(array, tgtN, rotCtr)
%ROTATE1DTO2D Embed 1D array on y-axis and rotate in 2D to align to tgtN.

theta = atan2(tgtN(2), tgtN(1));
rotMat = [cos(theta), -sin(theta);
          sin(theta),  cos(theta)];

% Embed 1D positions on y-axis, centroid at origin
y = array.positions(:).'; % 1×M
y = y - mean(y);
pos2 = [zeros(1, array.count);
        y];

posRot = rotMat * pos2 + rotCtr; % 2×M

arrayRot = array;
arrayRot.positions = posRot;
end


function [arrayRot, rotMat] = rotateTo3d(array, tgtN, rotCtr)
%ROTATETO3D Embed array into 3D (if needed) and rotate so +z aligns to tgtN.

% Rotation: +z -> tgtN (axis-angle)
z0 = [0; 0; 1];
axisRot = cross(z0, tgtN);
axisNorm = norm(axisRot);

if axisNorm < 1e-8
  % Parallel or anti-parallel
  if dot(z0, tgtN) > 0
    rotMat = eye(3);
  else
    % 180-degree rotation (one valid choice)
    rotMat = [-1 0  0;
               0 -1 0;
               0 0  1];
  end
else
  axisRot = axisRot / axisNorm;
  ang = acos(max(min(dot(z0, tgtN), 1), -1));
  rotMat = axang2rotm([axisRot.' ang]); % Aerospace Toolbox / Robotics Sys Toolbox
end

% Build 3D positions with centroid removal
switch array.dim
  case 1
    y = array.positions(:).';
    y = y - mean(y);
    pos3 = [zeros(1, array.count);
            y;
            zeros(1, array.count)];

  case 2
    pos2 = array.positions;
    pos2 = pos2 - mean(pos2, 2);
    pos3 = [pos2;
            zeros(1, array.count)];

  case 3
    pos3 = array.positions;
    pos3 = pos3 - mean(pos3, 2);

  otherwise
    error('rotateTo3d:UnsupportedDim', 'array.dim must be 1, 2, or 3.');
end

posRot = rotMat * pos3 + rotCtr; % 3×M

arrayRot = array;
arrayRot.positions = posRot;
end
