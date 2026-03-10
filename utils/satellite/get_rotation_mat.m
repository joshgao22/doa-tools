function Q = get_rotation_mat(z_axis_global)
%GET_ROTATION_MAT Construct a local-to-global rotation matrix.
%   Q = GET_ROTATION_MAT(z_axis_global)
%
%   This function constructs an orthonormal rotation matrix whose columns
%   represent the local coordinate axes expressed in the global frame.
%Input:
%   z_axis_global : 2×1 or 3×1 vector (column)
%     - For 2D: interpreted as the desired local x-axis direction
%       (e.g., [cos(theta); sin(theta)]).
%     - For 3D: interpreted as the desired local z-axis direction
%       (e.g., pointing toward Earth center).
%Output:
%   Q: d×d orthonormal rotation matrix, where d = length(z_axis_global)
%     - Columns of R are the local coordinate axes expressed in global
%       coordinates:
%           2D: R = [x_hat, y_hat]
%           3D: R = [x_hat, y_hat, z_hat]
%     - Transformation usage:
%           u_global = R * u_local
%           u_local  = R.' * u_global
%Notes
%   - In 2D, the local frame is constructed by rotating the x-axis by +90°
%     to obtain the y-axis.
%   - In 3D, the function automatically selects a helper vector that is
%     not aligned with the input z-axis in order to construct x_hat via
%     cross products.
%

arguments
  z_axis_global {mustBeVector, mustBeNumeric}
end

% Ensure column vector
if ismember(numel(z_axis_global), [2, 3])
  z_axis_global = z_axis_global(:);
  d = length(z_axis_global);
else
  error('Only 2D and 3D vectors are supported.');
end

if norm(z_axis_global) < 1e-8
  error('Input vector must be non-zero.');
end

% 2D case
if d == 2
  % Normalize the input (treated as local x-axis)
  x_hat = z_axis_global / norm(z_axis_global);

  % Local y-axis: +90 degree perpendicular rotation
  y_hat = [-x_hat(2); x_hat(1)];

  Q = [x_hat, y_hat];
  return
end

% 3D case
% Normalize desired local z-axis
z_hat = z_axis_global / norm(z_axis_global);

% Choose helper vector not parallel to z_hat
if abs(z_hat(1)) < 0.9
  v = [1; 0; 0];
else
  v = [0; 1; 0];
end

% Local x-axis via cross product
x_hat = cross(v, z_hat);
if norm(x_hat) < 1e-8
  error('Failed to construct orthogonal x-axis.');
end
x_hat = x_hat / norm(x_hat);

% Local y-axis completes right-handed frame
y_hat = cross(z_hat, x_hat);

Q = [x_hat, y_hat, z_hat];
end
