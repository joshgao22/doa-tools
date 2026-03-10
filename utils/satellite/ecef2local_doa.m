function doas_l = ecef2local_doa(pos_g, pos_l, Q_l2g)
%ECEF2LOCAL_DOA Convert user positions in ECEF coordinates to local DOAs
%as observed by an array or satellite.
%
%   doas_l = ecef2local_doa(pos_g, pos_l, Q_l2g)
%
%Inputs:
%   pos_g : 3 × K
%       User positions in the global ECEF frame (meters). Each column
%       corresponds to one user.
%
%   pos_l : 3 × 1
%       Array/satellite position in the global ECEF frame (meters).
%
%   Q_l2g : 3 × 3
%       Local-to-global rotation matrix of the array.
%       Each column is a unit local axis expressed in the global frame.
%       Thus, for a vector u_local:    u_global = Q_l2g * u_local.
%       This function uses its transpose to obtain local coordinates:
%           u_local = Q_l2g' * u_global.
%
%Output:
%   doas_l : 2 × K
%       Local azimuth–elevation DOAs in radians:
%           az = atan2(u_y, u_x)
%           el = asin(u_z)
%       where [u_x; u_y; u_z] is the line-of-sight unit vector expressed
%       in the array's local coordinate frame.
%
%Description:
%   For each user position pos_g(:,k), the function computes the line-of-sight
%   direction from the array at pos_l to the user, normalizes it, rotates the
%   resulting vector into the local array frame using Q_l2g', and finally
%   converts it to azimuth and elevation angles. The implementation is fully
%   vectorized over K users.
%
%Notes:
%   - The azimuth convention follows atan2(y, x) measured counterclockwise
%     from the local +x axis toward the local +y axis.
%   - The elevation is measured from the local x–y plane toward +z.
%

% -------- Basic dimension checks --------
if size(pos_g,1) ~= 3
  error('pos_g must be 3 × K.');
end
if ~isequal(size(pos_l), [3 1])
  error('pos_l must be 3 × 1.');
end
if ~isequal(size(Q_l2g), [3 3])
  error('Q_l2g must be 3 × 3.');
end

% Global line-of-sight vectors: array -> user (un-normalized)
V_glob = pos_g - pos_l;      % 3 × K

% Normalize column-wise
V_norm = V_glob ./ vecnorm(V_glob, 2, 1);   % 3 × K

% Transform to local coordinates: u_local = Q_l2g' * v_global
U_local = Q_l2g' * V_norm;                 % 3 × K

Ux = U_local(1, :);
Uy = U_local(2, :);
Uz = U_local(3, :);

% Azimuth and elevation
az_l = mod(atan2(Uy, Ux), 2*pi);  % 1 × K
el_l = asin(Uz);                  % 1 × K

doas_l = [az_l; el_l];            % 2 × K
end
