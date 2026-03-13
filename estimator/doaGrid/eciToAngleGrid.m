function angleGrid = eciToAngleGrid(eciDirectionGrid, rotMat)
%ECITOANGLEGRID Convert ECI direction vectors to local azimuth/elevation angles.
% Computes local azimuth and elevation angles (in radians) from ECI direction
% vectors using a local coordinate frame defined by the given rotation matrix.
%
%Syntax:
%   angleGrid = eciToAngleGrid(eciDirectionGrid, rotMat)
%
%Inputs:
%   eciDirectionGrid - 3×N matrix of unit direction vectors in the ECI frame,
%                      each column corresponds to one candidate direction
%
%   rotMat           - local-to-ECI rotation matrix / matrices
%                      - single array  : 3×3 numeric matrix
%                      - multiple arrays: cell array, each cell contains one
%                        3×3 numeric rotation matrix
%
%Output:
%   angleGrid        - If one rotation matrix: 2×N matrix [azimuth; elevation]
%                      If multiple rotation matrices: cell array, each cell
%                      contains one 2×N angle grid
%
%Notes:
%   - Each input direction vector is normalized internally for robustness.
%   - The local direction vector is computed as:
%         dirLocal = rotMat' * dirGlobal
%     since rotMat maps local coordinates to ECI.
%   - Azimuth is measured counterclockwise from the local x-axis in the
%     x-y plane; elevation is measured from the local horizontal plane.
%   - Returned angles are suitable for DOA estimation and beam steering
%     under the direction-only ECI measurement model.
%
%Example:
%   angleGrid = eciToAngleGrid(doaGrid.eciDirectionGrid, qMatrix);
%
%   rotMat = {qMatrix1, qMatrix2, qMatrix3};
%   angleGrid = eciToAngleGrid(doaGrid.eciDirectionGrid, rotMat);
%
%See also:
%   atan2, asin, vecnorm

arguments
  eciDirectionGrid (3,:) {mustBeNumeric}
  rotMat
end

% -------------------------------------------------------------------------
% Validate rotation matrix input
% -------------------------------------------------------------------------
if isnumeric(rotMat)
  if ~isequal(size(rotMat), [3, 3])
    error('eciToAngleGrid:InvalidRotationSize', ...
      'rotationMatrix must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  numArray = 1;
  isCellInput = false;

elseif iscell(rotMat)
  numArray = numel(rotMat);
  isCellInput = true;

  if isempty(rotMat)
    error('eciToAngleGrid:EmptyRotationCell', ...
      'rotationMatrix cell array must not be empty.');
  end

  for aa = 1:numArray
    currentRotation = rotMat{aa};
    if ~isnumeric(currentRotation) || ~isequal(size(currentRotation), [3, 3])
      error('eciToAngleGrid:InvalidRotationCell', ...
        'Each cell in rotationMatrix must contain one 3x3 numeric matrix.');
    end
  end

else
  error('eciToAngleGrid:InvalidRotationType', ...
    'rotationMatrix must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

numGrid = size(eciDirectionGrid, 2);

% -------------------------------------------------------------------------
% Normalize ECI direction vectors for robustness
% -------------------------------------------------------------------------
dirGlobal = eciDirectionGrid ./ vecnorm(eciDirectionGrid, 2, 1);

% -------------------------------------------------------------------------
% Initialize output
% -------------------------------------------------------------------------
if isCellInput
  angleGrid = cell(size(rotMat));
else
  angleGrid = zeros(2, numGrid);
end

% -------------------------------------------------------------------------
% Loop over each array
% -------------------------------------------------------------------------
for aa = 1:numArray
  if isCellInput
    currentRotation = rotMat{aa};
  else
    currentRotation = rotMat;
  end

  % Project ECI directions into local frame
  dirLocal = currentRotation' * dirGlobal;  % 3×N

  xCoordinate = dirLocal(1, :);
  yCoordinate = dirLocal(2, :);
  zCoordinate = dirLocal(3, :);

  % Clamp z-coordinate for numerical robustness before asin
  zCoordinate = max(min(zCoordinate, 1), -1);

  % Compute azimuth and elevation
  azimuth = atan2(yCoordinate, xCoordinate); % [-pi, pi]
  azimuth = mod(azimuth, 2*pi);              % [0, 2pi)
  elevation = asin(zCoordinate);             % [-pi/2, pi/2]

  if isCellInput
    angleGrid{aa} = [azimuth; elevation];
  else
    angleGrid = [azimuth; elevation];
  end
end

end