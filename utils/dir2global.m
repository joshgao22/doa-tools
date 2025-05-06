function u_global = dir2global(u_local, z_axis_local)
%DIR2GLOBAL Convert a direction vector from local frame to global frame.
%
%   u_global = dir2global(u_local, z_axis_local)
%
%   Inputs:
%     u_local       : 2x1 or 3x1 direction vector in local frame
%     z_axis_local  : 2x1 or 3x1 vector specifying local z-axis (same as in dir2local)
%
%   Output:
%     u_global      : direction vector expressed in global frame

    % Ensure column vectors
    u_local = u_local(:);
    z_axis_local = z_axis_local(:);

    % Dimension check
    d = length(z_axis_local);
    if length(u_local) ~= d
        error('u_local and z_axis_local must have the same dimension.');
    end

    % 2D case
    if d == 2
        z_hat = z_axis_local / norm(z_axis_local);
        x_hat = z_hat;
        y_hat = [-z_hat(2); z_hat(1)];
        R = [x_hat, y_hat];  % Local basis in global coordinates
        u_global = R * u_local;

    % 3D case
    elseif d == 3
        if norm(z_axis_local) < 1e-8
            error('z_axis_local must be non-zero.');
        end
        z_hat = z_axis_local / norm(z_axis_local);
        if abs(z_hat(1)) < 0.9
            v = [1; 0; 0];
        else
            v = [0; 1; 0];
        end
        x_hat = cross(v, z_hat);
        if norm(x_hat) < 1e-8
            error('Failed to compute orthogonal x-axis.');
        end
        x_hat = x_hat / norm(x_hat);
        y_hat = cross(z_hat, x_hat);
        R = [x_hat, y_hat, z_hat];
        u_global = R * u_local;

    else
        error('Only 2D or 3D vectors are supported.');
    end
end
