function rotationMat = get_rotation_mat(zAxisGlobal)
%GET_ROTATION_MAT Compute rotation matrix from a given global z-axis.
% Constructs a right-handed orthonormal coordinate frame whose local
% z-axis aligns with the specified global direction vector.
%
%Syntax:
%  rotationMat = get_rotation_mat(zAxisGlobal)
%
%Input:
%  zAxisGlobal : 3x1 vector specifying the desired local z-axis direction
%                expressed in the global coordinate frame
%
%Output:
%  rotationMat : 3x3 rotation matrix mapping local coordinates to the
%                global frame
%
%Construction:
%  zAxis = zAxisGlobal / ||zAxisGlobal||
%
%  Choose an auxiliary vector v that is not parallel to zAxis.
%  xAxis = normalize(cross(v, zAxis))
%  yAxis = cross(zAxis, xAxis)
%
%The rotation matrix is assembled as
%
%  rotationMat = [xAxis, yAxis, zAxis]
%
%which maps vectors from the local frame to the global frame.
%
%Example:
%  rotationMat = get_rotation_mat([0;0;1]);
%
%See also: cross, norm

arguments
  zAxisGlobal (3,1) {mustBeNumeric}
end

% -------------------------------------------------------------------------
% Normalize input direction to obtain local z-axis
% -------------------------------------------------------------------------
zAxis = zAxisGlobal / norm(zAxisGlobal);

% -------------------------------------------------------------------------
% Select an auxiliary vector not parallel to zAxis
% -------------------------------------------------------------------------
if abs(zAxis(1)) < 0.9
  auxVec = [1; 0; 0];
else
  auxVec = [0; 1; 0];
end

% -------------------------------------------------------------------------
% Construct orthonormal x- and y-axes
% -------------------------------------------------------------------------
xAxis = cross(auxVec, zAxis);
xAxis = xAxis / norm(xAxis);

yAxis = cross(zAxis, xAxis);

% -------------------------------------------------------------------------
% Assemble rotation matrix (local → global)
% -------------------------------------------------------------------------
rotationMat = [xAxis, yAxis, zAxis];

end