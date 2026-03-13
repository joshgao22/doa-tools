function [u_global, rotMat_l2g] = dir2global(u_local, z_axis_global)
%DIR2GLOBAL Convert a direction vector from a local frame to the global frame.
%   [u_global, rotMat_l2g] = dir2global(u_local, z_axis_global)
%
%Inputs:
%   u_local      : 2×1 or 3×1 direction vector expressed in the local frame
%   z_axis_global:
%       Direction of the local frame's z-axis expressed in the global frame.
%       - In 2D, this argument specifies the local x-axis
%         (e.g., [cos(theta); sin(theta)]).
%       - In 3D, this argument specifies the desired local z-axis direction.
%Outputs:
%   u_global     : u_local expressed in the global coordinate frame
%   rotMat_l2g   : Rotation matrix whose columns are the local axes
%                  expressed in the global frame (local → global).
%Usage:
%   u_global = rotMat_l2g * u_local
%Note:
%   This function uses get_rotation_mat(z_axis_global) to construct
%   the local-to-global rotation matrix consistent with dir2local.

arguments
  u_local {mustBeVector, mustBeNumeric}
  z_axis_global {mustBeVector, mustBeNumeric}
end

% Ensure column vectors
u_local = u_local(:);
z_axis_global = z_axis_global(:);

% Dimension check
d = length(z_axis_global);
if length(u_local) ~= d
  error('u_local and z_axis_global must have the same dimension.');
end

rotMat_l2g = get_rotation_mat(z_axis_global);
u_global = rotMat_l2g * u_local;

end
