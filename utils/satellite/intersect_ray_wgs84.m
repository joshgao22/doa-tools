function [p_g, t_star, is_resolved] = intersect_ray_wgs84(p0, d, E)
%INTERSECT_RAY_WGS84  Compute intersections between a ray (or rays) and the
%WGS-84 ellipsoid (zero-height surface).
%   [p_g, t_star, is_resolved] = intersect_ray_wgs84(p0, d, E)
%
%Description:
%   Given a ray (or multiple rays) parameterized as
%       r(t) = p0 + t * d ,   t ∈ ℝ,
%   this function computes the forward intersection point(s) with the WGS-84
%   reference ellipsoid defined by:
%       x^2 / a^2 + y^2 / a^2 + z^2 / b^2 = 1,
%   where a and b denote the semi-major (equatorial) and semi-minor (polar)
%   axes, respectively.
%Inputs:
%   p0 : 3×1 or 3×N
%        Ray origin(s) in ECEF coordinates (meters).
%        - If 3×1, it is replicated for each ray direction.
%        - If 3×N, each column is the origin of the corresponding ray.
%   d  : 3×N
%        Ray direction(s) in ECEF coordinates. Normalization is not required;
%        the quadratic formulation remains valid for arbitrary d.
%   E  : referenceEllipsoid / oblateSpheroid / referenceSphere
%        Provides the semi-major and semi-minor axes of the ellipsoid.
%Outputs:
%   p_g        : 3×N
%                Intersection points on the ellipsoid surface (ECEF, meters).
%                Columns with no forward intersection are filled with NaN.
%   t_star     : 1×N
%                Smallest positive root of the ray–ellipsoid intersection
%                equation for each ray. NaN if no forward intersection.
%   is_resolved: 1×N logical
%                True if a valid forward intersection exists for the ray.
%Notes:
%   - If the discriminant is negative, the ray does not intersect the ellipsoid.
%   - If both solutions t₁ and t₂ are non-positive, the intersection lies
%     entirely behind the ray origin.
%   - For each ray, the closest forward intersection is returned.

arguments
  p0 {mustBeNumeric}
  d  {mustBeNumeric}
  E  {mustBeA(E, ["referenceEllipsoid", "oblateSpheroid", ...
                  "referenceSphere"])} = wgs84Ellipsoid("meter")
end

% ----- size checks -----
if size(d,1) ~= 3
  error('d must be a 3xN array of direction vectors.');
end

[~, N] = size(d);  % number of rays

if numel(p0) == 3
  % Single origin replicated for all rays
  p0 = p0(:);
  p0 = repmat(p0, 1, N);   % 3 x N
elseif isequal(size(p0), [3, N])
  % One origin per ray
  % keep as is
else
  error('p0 must be either 3x1 or 3xN, consistent with d.');
end

% ----- ellipsoid parameters -----
a = E.SemimajorAxis;
b = E.SemiminorAxis;

% Components (each row is 1×N)
dx = d(1, :);
dy = d(2, :);
dz = d(3, :);

x0 = p0(1, :);
y0 = p0(2, :);
z0 = p0(3, :);

% Quadratic coefficients A t^2 + B t + C = 0 for each ray
A = (dx.^2 + dy.^2)/a^2 + (dz.^2)/b^2;
B = 2*((x0.*dx + y0.*dy)/a^2 + (z0.*dz)/b^2);
C = (x0.^2 + y0.^2)/a^2 + (z0.^2)/b^2 - 1;

% Discriminant
disc = B.^2 - 4.*A.*C;

% Initialize outputs
p_g        = nan(3, N);
t_star     = nan(1, N);
is_resolved = false(1, N);

% Rays with real intersection(s)
mask_disc = disc >= 0;
if ~any(mask_disc)
  return;  % no intersections at all
end

% Restrict to valid discriminants
A_v   = A(mask_disc);
B_v   = B(mask_disc);
disc_v = disc(mask_disc);

sqrt_disc_v = sqrt(disc_v);

t1 = (-B_v - sqrt_disc_v) ./ (2.*A_v);
t2 = (-B_v + sqrt_disc_v) ./ (2.*A_v);

% For each ray: pick smallest positive root
ts = [t1; t2];           % 2 x Nv
ts(ts <= 0) = NaN;       % discard non-forward intersections

[t_star_v, ~] = min(ts, [], 1, 'omitnan');  % 1 x Nv

valid_forward = ~isnan(t_star_v);
if ~any(valid_forward)
  return;
end

% Indices in the full set
idx_all = find(mask_disc);
idx_ok  = idx_all(valid_forward);

% Fill outputs
t_star(idx_ok)      = t_star_v(valid_forward);
is_resolved(idx_ok) = true;

% Compute intersection points for valid rays
t_row = t_star(idx_ok);            % 1 x N_ok
p0_ok = p0(:, idx_ok);             % 3 x N_ok
d_ok  = d(:, idx_ok);              % 3 x N_ok

p_g(:, idx_ok) = p0_ok + d_ok .* t_row;   % implicit expansion: 3xN_ok

end
