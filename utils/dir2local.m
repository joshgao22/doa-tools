function u_local = dir2local(u_global, z_axis_local)
%DIR2LOCAL Convert a global direction vector to a local frame.
%
%   u_local = dir2local(u_global, z_axis_local)
%
%   Inputs:
%     u_global       : 2x1 or 3x1 direction vector in global frame
%     z_axis_local   : 2x1 or 3x1 target "z-axis" (direction of local frame)
%
%   Output:
%     u_local        : direction vector expressed in the local frame
%
%   For 2D: z_axis_local = [cos(theta); sin(theta)]
%   For 3D: z_axis_local is normal vector; full orthonormal frame is constructed.

    % Convert to column vectors
    u_global = u_global(:);
    z_axis_local = z_axis_local(:);

    % Dimension check
    d = length(z_axis_local);
    if length(u_global) ~= d
        error('u_global and z_axis_local must have the same dimension.');
    end

    % 2D Case
    if d == 2
        % Normalize direction
        z_hat = z_axis_local / norm(z_axis_local);
        % Construct orthogonal x/y axes: x = z, y = rotate 90°
        x_hat = z_hat;
        y_hat = [-z_hat(2); z_hat(1)];  % Rotate +90°
        R = [x_hat, y_hat];  % Columns are local axes
        u_local = R' * u_global;

    % 3D Case
    elseif d == 3
        if norm(z_axis_local) < 1e-8
            error('z_axis_local must be non-zero.');
        end
        z_hat = z_axis_local / norm(z_axis_local);
        % Choose helper vector not parallel
        if abs(z_hat(1)) < 0.9
            v = [1; 0; 0];
        else
            v = [0; 1; 0];
        end
        x_hat = cross(v, z_hat);
        if norm(x_hat) < 1e-8
            error('Failed to compute orthogonal x-axis (z_hat may be aligned with helper).');
        end
        x_hat = x_hat / norm(x_hat);
        y_hat = cross(z_hat, x_hat);
        R = [x_hat, y_hat, z_hat];
        u_local = R' * u_global;

    else
        error('Only 2D or 3D vectors are supported.');
    end
end
