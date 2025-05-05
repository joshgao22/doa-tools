function [design_rot, rot_mat] = rotate_array(design, n, center)
%ROTATE_ARRAY Rotate an array to align its normal vector to a target direction.
%
%   design_rot = rotate_array(design, n, center) rotates the array described 
%   by the structure `design` to align its normal vector to the target direction 
%   described by the vector `n`. The rotation is performed by first translating 
%   the array's centroid to the origin, performing the rotation, and then translating 
%   the array back to the specified center. The rotation can be performed either in 
%   2D or 3D based on the input `n` and the dimension of the array (`design.dim`). 
%
%Inputs:
%   design - Array structure with the following fields:
%       design.element_positions - Positions of array elements (column vector).
%       design.element_count - Number of elements in the array.
%       design.dim - Dimension of the array (1 for 1D, 2 for 2D, 3 for 3D).
%   n - Target normal vector to align the array's normal to. 
%       2D case: 2-element vector [nx, ny].
%       3D case: 3-element vector [nx, ny, nz].
%   center - (Optional) A 2D or 3D vector specifying the center of rotation.
%            Default is the origin if not provided, based on the dimension of `n`.
%
%Outputs:
%   design_rot - Rotated array structure with updated element positions.

% Ensure n is a 2D or 3D vector
if ~(isvector(n) && (length(n) == 2 || length(n) == 3))
    error('Target normal vector n must be a 2D or 3D vector.');
end

% Normalize the target normal vector
n = n(:);  % Ensure column vector
n = n / norm(n);

% Handle center input
if nargin < 3 || isempty(center)
    center = zeros(size(n));  % Default center is origin of same dimension as n
else
    center = center(:);  % Ensure column vector
    if length(center) ~= length(n)
        error('Center and target normal vector n must have the same dimension (2D or 3D).');
    end
end

% Determine array dimension from design.dim
switch design.dim
    case 1
        % 1D array: check if 2D or 3D rotation is needed
        if length(n) == 2
            % 2D rotation for 1D array
            [design_rot, rot_mat] = rotate_1d_array(design, n, center, '2D');
            design_rot.dim = 2;
        elseif length(n) == 3
            % 3D rotation for 1D array
            [design_rot, rot_mat] = rotate_1d_array(design, n, center, '3D');
            design_rot.dim = 3;
        else
            error('Invalid normal vector length for 1D array.');
        end
    case 2
        % 2D array: only 3D rotation is allowed
        if length(n) == 3
            [design_rot, rot_mat] = rotate_2d_array(design, n, center);
            design_rot.dim = 3;
        else
            error('Only 3D rotation is allowed for 2D array.');
        end
    case 3
        % 3D array: only 3D rotation is allowed
        if length(n) == 3
            [design_rot, rot_mat] = rotate_3d_array(design, n, center);
            design_rot.dim = 3;  % No change in dimension
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
function [design_rot, rot_mat] = rotate_1d_array(design, n, center, type)
% Rotate a 1D or 2D array to align its normal vector to target direction (2D or 3D).

% If 2D rotation, calculate the angle in 2D
if strcmp(type, '2D')
    theta = atan2(n(2), n(1));  % compute rotation angle
    rot_mat = [cos(theta), -sin(theta); 
               sin(theta),  cos(theta)];  % 2D rotation matrix
    % Arrange elements along Y-axis and center at origin
    elem_pos_centroid = [zeros(1, design.element_count);
                         design.element_positions - mean(design.element_positions)];
    % Apply rotation and translate to the specified center
    elem_pos_rot = rot_mat * elem_pos_centroid + center;
    
elseif strcmp(type, '3D')
    % 3D rotation, we need the full rotation logic
    [design_rot, rot_mat] = rotate_general_array(design, n, center);
    return;
else
    error('Unsupported rotation type.');
end

design_rot.element_positions = elem_pos_rot;

end

% Helper function for 2D array rotation to 3D direction
function [design_rot, rot_mat] = rotate_2d_array(design, n, center)
    [design_rot, rot_mat] = rotate_general_array(design, n, center);
end

% Helper function for 3D rotation of 3D array
function [design_rot, rot_mat] = rotate_3d_array(design, n, center)
    [design_rot, rot_mat] = rotate_general_array(design, n, center);
end

% General rotation function for both 1D/2D arrays when rotating to 3D
function [design_rot, rot_mat] = rotate_general_array(design, n, center)
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
        rot_mat = axang2rot(axis_rot, angle_rot);
    end
    
    % Apply rotation for array
    switch design.dim
      case 1
        elem_pos_centroid = [zeros(1, design.element_count);
                     design.element_positions - mean(design.element_positions);
                     zeros(1, design.element_count)];
      case 2
        elem_pos_centroid = [design.element_positions - mean(design.element_positions, 2);
                     zeros(1, design.element_count)];
      case 3
        elem_pos_centroid = design.element_positions - mean(design.element_positions); 
    end
    elem_pos_rot = rot_mat * elem_pos_centroid + center;
    
    design_rot = design;
    design_rot.element_positions = elem_pos_rot;
end
