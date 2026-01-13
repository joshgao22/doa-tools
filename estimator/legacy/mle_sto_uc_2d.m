function [sp, p, noise_var] = mle_sto_uc_2d(R, n, design, wavelength, x0, varargin)
%MLE_STO_UC_2D ML estimator for the stochastic model with uncorrelated sources (2D DOA).
%   min_{Theta, p, sigma} logdet(S) + tr(S^{-1} R)
%   s.t. S = A(Theta) diag(p) A(Theta)^H + sigma I,
%        Theta = [az; el], p >= 0, sigma >= 0.
%Syntax:
%   sp = MLE_STO_UC_2D(R, n, design, wavelength, x0, ...);
%   sp = MLE_STO_UC_2D(R, n, f_steering, [], x0, ...);
%   [sp, p, noise_var] = MLE_STO_UC_2D(R, n, design, wavelength, x0, ...);
%   [sp, p, noise_var] = MLE_STO_UC_2D(R, n, f_steering, [], x0, ...);
%Inputs:
%   R          - Sample covariance matrix.
%   n          - Number of sources.
%   design     - Array design, or a steering-matrix function handle.
%                If function handle, it must be: A = design(wavelength, doas),
%                where doas is a 2 x n matrix [az; el] in radians.
%   wavelength - Wavelength.
%   x0         - (Optional) (3n+1)x1 vector:
%                x0 = [az1; el1; ...; azn; eln; p(1:n); noise_var].
%                If given in 'degree' or 'sin', it will be converted to radians.
%   Name–value options:
%       'Unit'    - 'radian' (default), 'degree', or 'sin'.
%       'Verbose' - true/false, display fmincon output if true.
%Outputs:
%   sp - Structure:
%           x        - 2 x n estimated DOAs (same as x_est).
%           y        - 1 x n (dummy, all ones).
%           x_est    - 2 x n estimated DOAs (row 1 az, row 2 el).
%           x_unit   - Unit of x_est.
%           resolved - true.
%           discrete - true.
%   p  - n x 1 vector of estimated source powers.
%   noise_var - Estimated noise power.

if nargin < 5
  x0 = [];
else
  if ~isempty(x0) && length(x0) ~= 3*n + 1
    error('The vector of the initial guess must have length 3*n+1 for 2D DOA.');
  end
end

unit    = 'radian';
verbose = false;
for ii = 1:2:nargin - 5
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

% prepare initial guess
if isempty(x0)
  x0 = zeros(3*n + 1, 1);

  % 使用 MVDR_2D 作为 DOA 初值（始终在弧度下工作）
  sp_mvdr = mvdr_2d(R, n, design, wavelength, 50);
  if sp_mvdr.resolved
    % x_est 是 2 x n，行 1 为 az，行 2 为 el
    doas0 = sp_mvdr.x_est;    % 2 x n in radians
    x0(1:2*n) = doas0(:);
  else
    warning('Failed to obtain a good initial guess from MVDR_2D. Using uniform DOAs.');
    az0 = linspace(0, 2*pi, n);
    el0 = pi/6 * ones(1, n);
    x0(1:2:2*n-1) = az0;
    x0(2:2:2*n)   = el0;
    doas0 = [az0; el0];
  end

  % 用当前 DOA 做一次最小二乘，初始化 p 和 noise_var
  if isa(design, 'function_handle')
    A = design(wavelength, doas0);
  else
    A = steering_matrix(design, wavelength, doas0);
  end
  B = [khatri_rao(conj(A), A), reshape(eye(size(A, 1)), [], 1)];
  z = real(B \ R(:));
  z(z < 0) = 0;
  x0(2*n + 1:end) = z;
else
  x0 = x0(:);
  % Unify to radians.
  switch (unit)
    case 'degree'
      x0(1:2*n) = deg2rad(x0(1:2*n));
    case 'sin'
      x0(1:2*n) = asin(x0(1:2*n));
    case 'radian'
      % do nothing
    otherwise
      error('Unexpected unit ''%s''.', unit);
  end
end

% prepare the optimization problem
if verbose
  options = optimoptions('fmincon', 'Display', 'iter');
else
  options = optimoptions('fmincon', 'Display', 'off');
end

% x = [az1; el1; ...; azn; eln; p(1:n); noise_var]
lb = zeros(3*n + 1, 1);
ub = inf(3*n + 1, 1);

% azimuth: [0, 2*pi]
lb(1:2:2*n-1) = 0;
ub(1:2:2*n-1) = 2*pi;

% elevation: [0, pi/2]
lb(2:2:2*n) = 0;
ub(2:2:2*n) = pi/2;

% p, noise_var: already [0, inf)

% -------- solve --------
x_opt = fmincon(@(x) nll_2d(R, design, wavelength, ...
  x(1:2*n), x(2*n + 1:2*n + n), x(end)), ...
  x0, [], [], [], [], lb, ub, [], options);

% -------- unpack solution --------
doas_vec  = x_opt(1:2*n);
p         = x_opt(2*n + 1:2*n + n);
noise_var = x_opt(end);

% reshape to 2 x n: row1=az, row2=el
doas_mat = reshape(doas_vec, 2, []);

% 按 azimuth 排序
[~, idx] = sort(doas_mat(1, :));
doas_mat = doas_mat(:, idx);
p        = p(idx);

% 输出单位转换
switch (unit)
  case 'degree'
    doas_out = rad2deg(doas_mat);
  case 'sin'
    doas_out = sin(doas_mat);
  case 'radian'
    doas_out = doas_mat;
  otherwise
    error('Unexpected unit ''%s''.', unit);
end

% -------- store results --------
sp = struct();
sp.x_est   = doas_out;     % 2 x n
sp.x       = sp.x_est;
sp.x_unit  = unit;
sp.y       = ones(1, n);   % dummy
sp.resolved = true;
sp.discrete = true;
end

function obj = nll_2d(R, design, wavelength, doas_vec, p, noise_var)
%NLL_2D Negative log-likelihood for 2D DOA.
%   doas_vec: 2n x 1, packed as [az1; el1; az2; el2; ...].

doas = reshape(doas_vec, 2, []);   % 2 x n

if isa(design, 'function_handle')
  A = design(wavelength, doas);
else
  A = steering_matrix(design, wavelength, doas);
end

S = A * bsxfun(@times, p, A') + noise_var * eye(size(A, 1));
S = 0.5 * (S + S');   % enforce Hermitian

[L, flag] = chol(S, 'lower');
if flag > 0
  obj = inf;
  return;
end

% logdet(S) = 2 * sum(log(diag(L)))
logdetS = 2 * sum(log(diag(L)));

% trace(S^{-1} R) via S^{-1}R = L'\(L\R)
SinvR = L' \ (L \ R);
obj   = logdetS + real(trace(SinvR));
end
