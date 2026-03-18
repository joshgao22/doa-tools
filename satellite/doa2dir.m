function u = doa2dir(ang)
%DOA2DIR Convert DOA parameters into unit direction vectors.
%   u = doa2dir(ang)
%
%Description:
%   Converts either broadside angles or azimuth/elevation pairs into
%   Cartesian unit direction vectors, using angle conventions consistent
%   with steeringMatrix.
%
%Input:
%   ang :
%     - 1 × n  : Broadside angles theta (radians).
%                For consistency with steeringMatrix, theta is measured
%                from +y toward +x in the x-y plane, so that
%                u = [sin(theta); cos(theta)].
%
%     - 2 × n  : [az; el] (radians).
%                az is measured counterclockwise from +x toward +y in the
%                x-y plane, and el is measured upward from the x-y plane.
%
%Output:
%   u :
%     - If input is 1 × n : output is 2 × n unit vectors.
%     - If input is 2 × n : output is 3 × n unit vectors.

arguments
  ang {mustBeNumeric}
end

[d, ~] = size(ang);

switch d
  case 1
    theta = ang;
    u = [sin(theta);
         cos(theta)];

  case 2
    az = ang(1, :);
    el = ang(2, :);

    u = [cos(el).*cos(az);
         cos(el).*sin(az);
         sin(el)];

  otherwise
    error('doa2dir:InputSize', ...
      'Input must have 1 or 2 rows (broadside angle or [az; el]).');
end

end