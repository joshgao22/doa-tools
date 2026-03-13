function [u_local, rotMat_l2g] = dir2local(u_global, z_axis_global)
%DIR2LOCAL Convert a direction vector from the global frame to a local frame.
%   [u_local, rotMat_l2g] = dir2local(u_global, z_axis_global)
%
%Inputs:
%   u_global      : 2×1 or 3×1 direction vector expressed in the global frame
%   z_axis_global :
%       Direction of the local frame's z-axis expressed in the global frame.
%       - In 2D, this argument specifies the local x-axis (e.g., [cosθ; sinθ]).
%       - In 3D, this argument specifies the desired local z-axis direction.
%Outputs:
%   u_local       : u_global expressed in the local coordinate frame
%   rotMat_l2g    : Rotation matrix whose columns are the local axes
%                   expressed in the global frame (local → global).
%Usage:
%   u_local = rotMat_l2g' * u_global
%Notes:
%   get_rotation_mat(z_axis_global) constructs a complete orthonormal
%   basis of the local frame consistent with the given axis direction.

arguments
  u_global {mustBeVector, mustBeNumeric}
  z_axis_global {mustBeVector, mustBeNumeric}
end

% Convert to column vectors
u_global = u_global(:);
z_axis_global = z_axis_global(:);

% Dimension check
d = length(z_axis_global);
if length(u_global) ~= d
  error('u_global and z_axis_global must have the same dimension.');
end

rotMat_l2g = get_rotation_mat(z_axis_global);
u_local = rotMat_l2g' * u_global;

end
