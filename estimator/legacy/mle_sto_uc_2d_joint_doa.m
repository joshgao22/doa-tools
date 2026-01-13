function [sp, p, noise_var] = mle_sto_uc_2d_joint_doa(R_all, n, design_all, ...
  wavelength, geom, x0, varargin)
%MLE_STO_UC_2D_JOINT_DOA Multi-array ML estimator (stochastic model,
%uncorrelated sources, 2D DOA).
%   [sp, p, noise_var] = mle_sto_uc_2d_joint(R_all, n, design_all, ...
%       wavelength, geom, x0, ...)
%Maximum-likelihood problem:
%   min_{Theta_ref, p, sigma_1,...,sigma_L}  sum_{ell=1}^L
%       [ logdet(S_ell) + tr(S_ell^{-1} R_ell) ]
%   s.t.
%       S_ell = A_ell(Theta_ell) diag(p) A_ell(Theta_ell)^H + sigma_ell I,
%       Theta_ell = P_ell(Theta_ref)  (geometry-based mapping),
%       p >= 0,  sigma_ell >= 0.
%   Here, Theta_ref collects the DOAs in the reference array's local frame,
%   and all other arrays share the same source powers p but have their own
%   noise variances sigma_ell and local DOAs obtained via geometric mapping.
%
%Syntax:
%   sp = mle_sto_uc_2d_multi(R_all, n, design_all, wavelength, geom, x0, ...);
%   [sp, p, noise_var] = mle_sto_uc_2d_multi(...);
%Inputs:
%   R_all      - 1×L cell array, {R1, ..., RL}.
%                R_ell is the sample covariance matrix of the ell-th array,
%                of size M_ell × M_ell.
%   n          - Number of sources.
%   design_all - 1×L cell array, {design1, ..., designL}.
%                Each entry is either:
%                  - an array design structure compatible with steering_matrix, or
%                  - a steering function handle:
%                        A = design_l(wavelength, doas_l),
%                    where doas_l is 2×n, [az; el] in radians.
%   wavelength - Wavelength of the impinging signals.
%   geom       - Geometry information structure:
%                  geom.L         : number of arrays/satellites (L)
%                  geom.ref_index : index of the reference array (e.g., 1)
%                  geom.pos{ell}  : 3×1 ECEF position of array/satellite ell
%                  geom.Q_l2g{ell}: 3×3 local-to-global rotation matrix for ell
%                  geom.E         : wgs84Ellipsoid object (Mapping Toolbox)
%   x0         - (Optional) initial parameter vector:
%                  x = [az_ref(1:n); el_ref(1:n); p(1:n); sigma(1:L)]
%                If empty or omitted, an internal initialization is used.
%Name–value pairs:
%   'Unit'    - 'radian' (default), 'degree', or 'sin'.
%               Only affects the DOA representation at the input (if x0 is
%               provided) and output; the optimization is always performed
%               in radians internally.
%
%   'Verbose' - true/false, whether to display fmincon output.
%Outputs:
%   sp - Structure summarizing the DOA estimates:
%          sp.x_est   : 2×n estimated DOAs in the reference array frame
%                       (row 1: azimuth, row 2: elevation)
%          sp.x       : alias of x_est
%          sp.x_unit  : unit of x_est ('radian' / 'degree' / 'sin')
%          sp.resolved, sp.discrete : kept for compatibility (dummy fields)
%   p  - n×1 vector of estimated source powers (shared by all arrays).
%   noise_var - L×1 vector of estimated noise variances for each array.
%

% ---------- Basic consistency checks ----------
L = geom.L;

if numel(R_all) ~= L || numel(design_all) ~= L
  error('R_all and design_all must be 1xL cell arrays consistent with geom.L.');
end

if nargin < 6
  x0 = [];
end

% ---------- Parse name–value options ----------
unit    = 'radian';
verbose = false;
for ii = 1:2:(nargin - 6)
  option_name  = varargin{ii};
  option_value = varargin{ii + 1};
  switch lower(option_name)
    case 'unit'
      unit = option_value;
    case 'verbose'
      verbose = option_value;
    otherwise
      error('Unknown option ''%s''.', option_name);
  end
end

% ---------- Dimension bookkeeping ----------
num_theta = 2*n;   % Reference array 2D DOAs: [az1; el1; ...; azn; eln]
num_p     = n;     % Source powers
num_sigma = L;     % One noise variance per array

len_x = num_theta + num_p + num_sigma;

% ---------- Initialization of x ----------
if isempty(x0)
  x0 = zeros(len_x, 1);

  % 1) Reference array DOA initialization via MVDR_2D
  ell0 = geom.ref_index;
  R0   = R_all{ell0};
  des0 = design_all{ell0};

  sp_mvdr = mvdr_2d(R0, n, des0, wavelength, 50);
  if sp_mvdr.resolved
    doas_ref0 = sp_mvdr.x_est;   % 2×n, radians
  else
    warning('MVDR_2D failed. Using a uniform DOA grid as initial guess.');
    az0 = linspace(0, 2*pi, n);
    el0 = pi/6 * ones(1, n);
    doas_ref0 = [az0; el0];
  end

  x0(1:num_theta) = doas_ref0(:);

  % 2) Initialize p and sigma from the reference array via a simple LS fit
  if isa(des0, 'function_handle')
    A0 = des0(wavelength, doas_ref0);
  else
    A0 = steering_matrix(des0, wavelength, doas_ref0);
  end
  B  = [khatri_rao(conj(A0), A0), reshape(eye(size(A0, 1)), [], 1)];
  z  = real(B \ R0(:));
  z(z < 0) = 0;
  p0     = z(1:n);
  sigma0 = z(end);

  x0(num_theta+1 : num_theta+num_p)         = p0;
  x0(num_theta+num_p+1 : num_theta+num_p+L) = sigma0 * ones(L, 1);
else
  x0 = x0(:);
  if length(x0) ~= len_x
    error('Initial x0 must have length 2*n + n + L = %d.', len_x);
  end
  % Convert initial DOAs to radians if necessary
  switch (unit)
    case 'degree'
      x0(1:num_theta) = deg2rad(x0(1:num_theta));
    case 'sin'
      x0(1:num_theta) = asin(x0(1:num_theta));
    case 'radian'
      % already in radians
    otherwise
      error('Unexpected unit ''%s''.', unit);
  end
end

% ---------- fmincon options ----------
if verbose
  options = optimoptions('fmincon', 'Display', 'iter');
else
  options = optimoptions('fmincon', 'Display', 'off');
end

% ---------- Bounds on optimization variables ----------
lb = zeros(len_x, 1);
ub = inf(len_x, 1);

% Azimuth: [0, 2*pi]
lb(1:2:2*n-1) = 0;
ub(1:2:2*n-1) = 2*pi;

% Elevation: [0, pi/2]
lb(2:2:2*n) = 0;
ub(2:2:2*n) = pi/2;

% p, sigma: [0, +inf] already covered by default lb/ub

% ---------- Solve the ML problem ----------
x_opt = fmincon(@(x) nll_2d_multi(R_all, design_all, wavelength, geom, n, x), ...
  x0, [], [], [], [], lb, ub, [], options);

% ---------- Unpack solution ----------
doas_ref_vec = x_opt(1:num_theta);
p            = x_opt(num_theta+1 : num_theta+num_p);
sigma_vec    = x_opt(num_theta+num_p+1 : end);   % L×1

doas_ref_mat = reshape(doas_ref_vec, 2, []);

% Sort sources by azimuth in the reference frame (for a consistent ordering)
[~, idx]     = sort(doas_ref_mat(1, :));
doas_ref_mat = doas_ref_mat(:, idx);
p            = p(idx);

% ---------- Output unit conversion ----------
switch (unit)
  case 'degree'
    doas_out = rad2deg(doas_ref_mat);
  case 'sin'
    doas_out = sin(doas_ref_mat);
  case 'radian'
    doas_out = doas_ref_mat;
  otherwise
    error('Unexpected unit ''%s''.', unit);
end

% ---------- Pack output structure ----------
sp = struct();
sp.x_est    = doas_out;   % 2×n, from the reference array viewpoint
sp.x        = sp.x_est;
sp.x_unit   = unit;
sp.y        = ones(1, n);  % dummy (kept for compatibility)
sp.resolved = true;
sp.discrete = true;

noise_var = sigma_vec;    % L×1

end

function obj = nll_2d_multi(R_all, design_all, wavelength, geom, n, x)
%NLL_2D_MULTI Negative log-likelihood for the multi-array stochastic
%2D-DOA model.
%
%   x = [az_ref(1:n); el_ref(1:n); p(1:n); sigma(1:L)]
%
%   The reference DOAs are first mapped to ground intersection points, and
%   then to local DOAs for each array via the known geometry. For each
%   array, the conditional covariance S_ell is constructed and the
%   contribution logdet(S_ell) + tr(S_ell^{-1} R_ell) is accumulated.

L = geom.L;
num_theta = 2*n;
num_p     = n;
num_sigma = L;

doas_ref_vec = x(1:num_theta);
p            = x(num_theta+1 : num_theta+num_p);
sigma_vec    = x(num_theta+num_p+1 : num_theta+num_p+num_sigma);

doas_ref = reshape(doas_ref_vec, 2, []);   % 2×n

% ---------- 1) Reference-array DOAs -> ground intersection points ----------
ref = geom.ref_index;
pos_ref = geom.pos{ref};       % 3×1 ECEF position of reference array
Q_ref   = geom.Q_l2g{ref};     % 3×3 local->global rotation
E       = geom.E;              % wgs84Ellipsoid

% Local unit direction vectors (reference array frame)
u_ref_local  = doa2dir(doas_ref);    % 3×n
% Corresponding global unit direction vectors
u_ref_global = Q_ref * u_ref_local;  % 3×n

% Ray–WGS84 ellipsoid intersections (ground points)
[pos_g, ~, is_resolved] = intersect_ray_wgs84(pos_ref, u_ref_global, E);

if ~all(is_resolved)
    obj = 1e20;
    return;
end

% ---------- 2) Ground points -> local DOAs at each array + NLL accumulation ----------
obj = 0;

for ell = 1:L
  Rl      = R_all{ell};
  designl = design_all{ell};

  pos_l = geom.pos{ell};
  Q_l   = geom.Q_l2g{ell};

  % 2.1 Ground points -> local DOAs at array ell
  doas_l = ecef2local_doa(pos_g, pos_l, Q_l);

  % 2.2 Steering matrix for array ell
  if isa(designl, 'function_handle')
    A_l = designl(wavelength, doas_l);
  else
    A_l = steering_matrix(designl, wavelength, doas_l);
  end

  % 2.3 Model covariance matrix S_l = A_l diag(p) A_l^H + sigma_l I
  sigma_l = sigma_vec(ell);
  S_l = A_l * (bsxfun(@times, p, A_l')) + sigma_l * eye(size(A_l, 1));
  S_l = 0.5 * (S_l + S_l');   % enforce Hermitian numerically

  % 2.4 Use Cholesky factorization to compute logdet and tr(S^{-1} R)
  [L_chol, flag] = chol(S_l, 'lower');
  if flag > 0
    obj = inf;
    return;
  end

  logdetS = 2 * sum(log(diag(L_chol)));
  SinvR   = L_chol' \ (L_chol \ Rl);
  obj     = obj + logdetS + real(trace(SinvR));
end

end
