%[text] # 随机信号模型下的单源 DoA 集中似然估计（网格搜索）
%[text] 随机信号（stochastic）信号模型下，DoA 估计问题可以表示为如下优化问题
%[text]{"align":"center"} $\\underset{\\mathbf{\\theta}, \\mathbf{p}, \\sigma}{\\mathrm{minimize}}\\quad \\log\\det(\\mathbf{R}) + \\mathrm{tr}\\left(\\mathbf{R}^{-1} \\hat{\\mathbf{R}}\\right), \\\\\n\\text{subject to}\\quad \\mathbf{R} = \\mathbf{A}(\\mathbf{\\theta}) \\mathrm{diag}(\\mathbf{p})  \\mathbf{A}^\\mathrm{H}(\\mathbf{\\theta}) + \\sigma \\mathbf{I}\\\\\n\\qquad\\qquad\\,\\,\\,\\,\\theta\_k \\in\[-\\pi/2,\\pi/2\] \\\\\n\\qquad\\qquad\\,\\,\\,\\,\\mathbf{p}\_k, \\sigma \\geq 0,$
%[text] 其中$\\mathbf{R}$和$\\hat{\\mathbf{R}}$分别是接收信号的协方差矩阵及其估计（也即样本协方差矩阵），$\\mathbf{p}$和$\\mathbf{\\theta}$分别是$K$个信号源的功率和来向，$\\sigma\n$是噪声功率。为了降低问题求解的维度，将$\\mathbf{p}$和$\\sigma\n$的估计分别用$\\mathbf{\\theta}$表示为
%[text]{"align":"center"} $\\hat{\\sigma}(\\mathbf{\\theta}) = \\frac{1}{M-K} \\mathrm{tr}\\left\[ \\left(\\mathbf{I}\_M - \\mathbf{A}(\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{-1} \\mathbf{A}^\\mathrm{H} \\right) \\hat{\\mathbf{R}} \\right\],$
%[text]{"align":"center"} $\\hat{\\mathbf{P}}(\\mathbf{\\theta}) = (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{-1} \\mathbf{A}^\\mathrm{H} \\hat{\\mathbf{R}} \\mathbf{A} (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{-1} - \\hat{\\sigma} (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{-1},\n$
%[text] 则可以得到$\\mathbf{\\theta}$的集中似然估计问题
%[text]{"align":"center"} $\\mathbf{\\theta}\_\\mathrm{CSML} = \\arg\\min\_{\\mathbf{\\theta}}\\log \\det \\left\[ \\mathbf{A} \\hat{\\mathbf{P}} \\mathbf{A}^\\mathrm{H} + \\hat{\\sigma} \\mathbf{I}\_M \\right\].$
%[text] 本函数使用网格搜索求解该问题，由于搜索维度随信号源$K$指数增长，因此本函数仅对单个信号源进行求解。
%[text] ## Syntax
%[text] `sp = mle_sto_con_grid_1d(R, k, design, wavelength, grid_size, ...)`
%[text] `sp = mle_sto_con_grid_1d(R, k, f_steering, [], grid_size, ...)`
%[text] `[sp, p, noise_var] = mle_sto_con_grid_1d(R, k, design, wavelength, grid_size, ...)`
%[text] `[sp, p, noise_var] = mle_sto_con_grid_1d(R, k, f_steering, [], grid_size, ...)`
%[text] ## Inputs
%[text] `R` — Sample covariance matrix. 
%[text] `k` — Number of sources. Can only be 1.
%[text] `design` — Array design. Can also be a function handle that generates a steering matrix. This function must take two arguments, wavelength and the doa vector. 
%[text] `wavelength` — Wavelength. 
%[text] `grid_size` —  Number of grid points used.
%[text] `...` — Options: 
%[text]  `'Unit'` — Can be `'radian'`, `'degree'`, or `'sin'`. Default value is `'radian'`. 
%[text]  `'RefineEstimates'` - If set to true, will refine the estimated direction of arrivals around the grid.
%[text] ## Output
%[text] `sp` — Spectrum structure with the following fields: 
%[text]     x - An 1 x grid\_size vector.
%[text]     y - An 1 x grid\_size vector. Calling `plot(x, y)` will plot the spectrum.
%[text]  `x_est` — An $1\\times k$ vector storing the estimated DOAs. 
%[text]  `x_unit` — The same as the unit specified by `'Unit'`.
%[text]  `resolved` — True if the number of peaks in the spectrum is greater or equal to the number of sources.
%[text]  `discrete` — Constant value `false`.
%[text] `p` — An $k\\times 1$ vector of estimated source powers.
%[text] `noise_var` — Estimated noise power.
%[text] ## Function Definition
function [sp, p, noise_var] = mle_sto_con_grid_1d(R, k, design, wavelength, grid_size, varargin)
unit = 'radian';
refine_estimates = false;
for ii = 1:2:nargin-5
  option_name = varargin{ii};
  option_value = varargin{ii+1};
  switch lower(option_name)
    case 'unit'
      unit = lower(option_value);
    case 'refineestimates'
      refine_estimates = true;
    otherwise
      error('Unknown option ''%s''.', option_name);
  end
end
%[text] ### compute spectrum
% discretize and create the corresponding steering matrix
[doa_grid_rad, doa_grid_display, ~] = default_doa_grid(grid_size, unit, 1);
% compute spectrum
sp_intl = compute_spectrum(R, design, wavelength, doa_grid_rad);
[x_est, x_est_idx, resolved] = find_doa_from_spectrum(doa_grid_display, sp_intl, k);
% refine
if resolved && refine_estimates
  switch unit
    case 'radian'
      f_obj = @(x) compute_spectrum(R, design, wavelength, x);
    case 'degree'
      f_obj = @(x) compute_spectrum(R, design, wavelength, deg2rad(x));
    case 'sin'
      f_obj = @(x) compute_spectrum(R, design, wavelength, asin(x));
    otherwise
      error('Invalid unit ''%s''.', unit);
  end
  x_est = refine_grid_estimates(f_obj, doa_grid_display, x_est_idx);
end
% Compute final estimates of A, P, and noise variance
if ishandle(design)
  A = design(wavelength, x_est);
else
  A = steering_matrix(design, wavelength, x_est);
end
G = A' * A;                                     % Gram matrix
Ginv = inv(G);                                  % Inverse Gram matrix

m = design.element_count;

% Estimate signal covariance matrix P and noise variance
noise_var = real(trace((eye(m) - A * Ginv * A') * R)) / (m - k);
p = diag(real(Ginv * A' * R * A * Ginv - noise_var * Ginv));
%[text] ### store results
sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est;
sp.x_unit = unit;
sp.y = sp_intl;
sp.resolved = resolved;
sp.discrete = false;
end
%%
%[text] # Concentrated log-likelihood for fixed DOAs
%[text] This function evaluates the concentrated log-likelihood cost given DOAs theta, using closed-form estimates of signal covariance and noise variance.
function v = compute_spectrum(R, design, wavelength, theta)
n_theta = length(theta);
v = zeros(1,n_theta);
for ii = 1:n_theta
  if ishandle(design)
    A = design(wavelength, theta(ii));
  else
    A = steering_matrix(design, wavelength, theta(ii));
  end

  m = size(R, 1);
  G = A' * A;
  Ginv = inv(G);

  % Closed-form noise variance estimator
  noise_var_est = real(trace((eye(m) - A * Ginv * A') * R)) / (m - 1);

  % Estimate signal powers (diagonal of P)
  P_est = real(Ginv * A' * R * A * Ginv - noise_var_est * Ginv);
  P_est = max(P_est, 0);

  % Construct model covariance matrix
  S = A * P_est * A' + noise_var_est * eye(m);

  % Negative log-likelihood (only log-det term remains in profile form)
  v(ii) = -real(log(det(S)));
end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:styles]
%   data: {"code":{"fontFamily":"consolaslxgw"}}
%---
