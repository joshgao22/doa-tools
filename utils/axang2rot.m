function R = axang2rot(axis, angle)
%AXANG2ROT Generates a 3x3 rotation matrix from an axis-angle pair.
%
%   R = axang2rot(axis, angle)
%
%   Inputs:
%     axis  : 3x1 rotation axis (not necessarily unit length)
%     angle : rotation angle in radians
%
%   Output:
%     R     : 3x3 rotation matrix (right-hand rule)

    axis = axis(:);
    if norm(axis) == 0
        error('Rotation axis must be non-zero.');
    end
    u = axis / norm(axis);

    ux = u(1);
    uy = u(2);
    uz = u(3);

    c = cos(angle);
    s = sin(angle);
    C = 1 - c;

    R = [c + ux^2*C,     ux*uy*C - uz*s, ux*uz*C + uy*s;
         uy*ux*C + uz*s, c + uy^2*C,     uy*uz*C - ux*s;
         uz*ux*C - uy*s, uz*uy*C + ux*s, c + uz^2*C];
end
