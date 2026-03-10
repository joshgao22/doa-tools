function u = doa2dir(ang)
%DOA2DIR Convert 1-angle or 2-angle DOA representation into a unit
%direction vector.
%   u = doa2dir(ang)
%
%Description:
%   This function converts angular DOA parameters into a Cartesian
%   unit-norm direction vector. 
%Input:
%   ang :
%     - 1 × n  : Single angle θ (radians). Interpreted as a 2-D DOA, 
%                producing direction vectors [cosθ; sinθ].
%     - 2 × n  : [az; el] (radians). Interpreted as a 3-D DOA with
%                azimuth and elevation.
%Angle conventions:
%   - Azimuth  az : angle measured counterclockwise from +x toward +y.
%   - Elevation el: angle measured upward from the x–y plane.
%   - Resulting direction vector follows the standard right-handed frame.
%Output:
%   u :
%     - If input is 1 × n : output is 2 × n unit vectors.
%     - If input is 2 × n : output is 3 × n unit vectors.
%Examples:
%   u = doa2dir(pi/4);          % 2-D direction vector
%   u = doa2dir([az; el]);      % 3-D direction vector
%

arguments
  ang {mustBeNumeric}
end

[d, ~] = size(ang);

switch d
  case 1
    theta = ang;
    u = [cos(theta);
         sin(theta)];

  case 2
    az = ang(1, :);
    el = ang(2, :);

    % Standard 3-D spherical coordinate mapping (right-handed)
    u = [cos(el).*cos(az);
         cos(el).*sin(az);
         sin(el)];

  otherwise
    error('doa2dir: Input must have 1 or 2 rows (single angle or [az; el]).');
end

end
