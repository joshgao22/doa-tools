function [u, du] = doa2dirJacobian(ang)
%DOA2DIRJACOBIAN Convert DOA parameters into unit direction vectors and
%their Jacobians.
%
%Syntax:
%   [u, du] = doa2dirJacobian(ang)
%
%Inputs:
%   ang - DOA parameter matrix
%         - 1×N : broadside angles theta (radians)
%         - 2×N : [az; el] (radians)
%
%Outputs:
%   u  - Unit direction vectors
%        - 2×N for broadside input
%        - 3×N for [az; el] input
%
%   du - Jacobian structure
%        - broadside input : du.theta (same size as u)
%        - [az; el] input  : du.az, du.el (same size as u)
%
%Description:
%   Computes Cartesian unit direction vectors and their first-order
%   derivatives with respect to the DOA parameters, using angle conventions
%   consistent with doa2dir and steeringMatrix.
%
%Notes:
%   - For broadside input theta, the direction vector is
%         u = [sin(theta); cos(theta)]
%     so theta is measured from +y toward +x in the x-y plane.
%
%   - For [az; el] input, the direction vector is
%         u = [cos(el)cos(az); cos(el)sin(az); sin(el)]
%     where az is measured counterclockwise from +x toward +y in the x-y
%     plane, and el is measured upward from the x-y plane.
%
%See also:
%   doa2dir, steeringMatrix

arguments
  ang {mustBeNumeric}
end

[numRow, ~] = size(ang);

du = struct();

switch numRow
  case 1
    theta = ang;

    sinTheta = sin(theta);
    cosTheta = cos(theta);

    u = [sinTheta;
         cosTheta];

    du.theta = [cosTheta;
                -sinTheta];

  case 2
    az = ang(1, :);
    el = ang(2, :);

    cosAz = cos(az);
    sinAz = sin(az);
    cosEl = cos(el);
    sinEl = sin(el);

    u = [cosEl .* cosAz;
         cosEl .* sinAz;
         sinEl];

    du.az = [-cosEl .* sinAz;
              cosEl .* cosAz;
              zeros(1, size(ang, 2))];

    du.el = [-sinEl .* cosAz;
             -sinEl .* sinAz;
              cosEl];

  otherwise
    error('doa2dirJacobian:InputSize', ...
      'Input must have 1 or 2 rows (broadside angle or [az; el]).');
end

end
