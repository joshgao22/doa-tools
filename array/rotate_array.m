function design_rot = rotate_array(design, n)
%ROTATE_ARRAY Rotate an array to align its normal vector to a target direction.
%Inputs:
%   design - Array structure. Fields:
%       design.element_positions - Positions of array elements.
%       design.element_count - Number of elements.
%       design.dim - Dimension of the array (1 for 1D, 2 for 2D, 3 for 3D).
%   n - Target normal vector.
%       2D case: 2-element vector [nx, ny].
%       3D case: 3-element vector [nx, ny, nz].
%Outputs:
%   design_rot - Rotated array structure.

% Validate inputs
if ~(isvector(n) && (length(n) == 2 || length(n) == 3))
    error('Target normal vector n must be 2D or 3D.');
end

% Determine array dimension from design.dim
switch design.dim
    case 1
        % 1D array: check if 2D or 3D rotation is needed
        if length(n) == 2
            % 2D rotation for 1D array
            design_rot = rotate_1d_array(design, n, '2D');
            design_rot.dim = 2;
        elseif length(n) == 3
            % 3D rotation for 1D array
            design_rot = rotate_1d_array(design, n, '3D');
            design_rot.dim = 3;
        else
            error('Invalid normal vector length for 1D array.');
        end
    case 2
        % 2D array: only 3D rotation is allowed
        if length(n) == 3
            design_rot = rotate_2d_array(design, n);
            design_rot.dim = 3;
        else
            error('Only 3D rotation is allowed for 2D array.');
        end
    case 3
        % 3D array: only 3D rotation is allowed
        if length(n) == 3
            design_rot = rotate_3d_array(design, n);
            design_rot.dim = 3; % No change in dimension
        else
            error('Only 3D rotation is allowed for 3D array.');
        end
    otherwise
        error('Unsupported design.dim. Only 1D, 2D, or 3D arrays are supported.');
end

% Update element positions after rotation
design_rot.element_positions = design_rot.element_positions;

end

% General helper function for 1D and 2D array rotation
function design_rot = rotate_1d_array(design, n, type)
% Rotate a 1D or 2D array to align its normal vector to target direction (2D or 3D).

% If 2D rotation, calculate the angle in 2D
if strcmp(type, '2D')
    theta = atan2(n(2), n(1));  % compute rotation angle
    rot_mat = [cos(theta), -sin(theta); 
               sin(theta),  cos(theta)];  % 2D rotation matrix
    % Arrange elements along Y-axis and center at origin
    elem_pos_centroid = [zeros(1, design.element_count);
                         design.element_positions - mean(design.element_positions)];
    % Apply rotation
    elem_pos_rot = rot_mat * elem_pos_centroid;
    
elseif strcmp(type, '3D')
    % 3D rotation, we need the full rotation logic
    design_rot = rotate_general_array(design, n);
    return;
else
    error('Unsupported rotation type.');
end

design_rot.element_positions = elem_pos_rot;

end

% General rotation function for both 1D/2D arrays when rotating to 3D
function design_rot = rotate_general_array(design, n)
    % Rotate to a 3D direction using axis-angle method
    original_n = [0; 0; 1];  % Original normal is along +Z axis
    axis_rot = cross(original_n, n);
    axis_norm = norm(axis_rot);
    
    if axis_norm < 1e-8
        % If vectors are parallel or anti-parallel
        if dot(original_n, n) > 0
            rot_mat = eye(3);  % No rotation
        else
            rot_mat = [-1 0 0; 0 -1 0; 0 0 1];  % 180-degree rotation
        end
    else
        axis_rot = axis_rot / axis_norm;
        angle_rot = acos(dot(original_n, n));
        rot_mat = axis_angle_to_matrix(axis_rot, angle_rot);
    end
    
    % Apply rotation for array
    elem_pos_centroid = [zeros(3, design.element_count);
                         design.element_positions - mean(design.element_positions)];
    elem_pos_rot = rot_mat * elem_pos_centroid;
    
    design_rot.element_positions = elem_pos_rot;
end

% Helper function for 2D array rotation to 3D direction
function design_rot = rotate_2d_array(design, n)
    design_rot = rotate_general_array(design, n);
end

% Helper function for 3D rotation of 3D array
function design_rot = rotate_3d_array(design, n)
    design_rot = rotate_general_array(design, n);
end

% Helper function to generate 3D rotation matrix from axis and angle
function R = axis_angle_to_matrix(axis, angle)
    % Generate rotation matrix from axis and angle
    x = axis(1);
    y = axis(2);
    z = axis(3);
    c = cos(angle);
    s = sin(angle);
    C = 1 - c;

    R = [x*x*C + c,   x*y*C - z*s, x*z*C + y*s;
         y*x*C + z*s, y*y*C + c,   y*z*C - x*s;
         z*x*C - y*s, z*y*C + x*s, z*z*C + c];
end