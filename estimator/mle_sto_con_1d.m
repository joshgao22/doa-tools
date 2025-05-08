%[text] # 随机信号模型下的 DoA 集中似然估计
%[text] 随机信号（stochastic）信号模型下，DoA 估计问题可以表示为如下优化问题
%[text]{"align":"center"} $\\underset{\\mathbf{\\theta}, \\mathbf{p}, \\sigma}{\\mathrm{minimize}}\\quad \\log\\det(\\mathbf{R}) \+ \\mathrm{tr}\\left(\\mathbf{R}^{\-1} \\hat{\\mathbf{R}}\\right), \\\\\n\\text{subject to}\\quad \\mathbf{R} = \\mathbf{A}(\\mathbf{\\theta}) \\mathrm{diag}(\\mathbf{p})  \\mathbf{A}^\\mathrm{H}(\\mathbf{\\theta}) \+ \\sigma \\mathbf{I}\\\\\n\\qquad\\qquad\\,\\,\\,\\,\\theta\_k \\in[\-\\pi/2,\\pi/2] \\\\\n\\qquad\\qquad\\,\\,\\,\\,\\mathbf{p}\_k, \\sigma \\geq 0,$
%[text] 其中$\\mathbf{R}$和$\\hat{\\mathbf{R}}$分别是接收信号的协方差矩阵及其估计（也即样本协方差矩阵），$\\mathbf{p}$和$\\mathbf{\\theta}$分别是$K$个信号源的功率和来向，$\\sigma\n$是噪声功率。为了降低问题求解的维度，将$\\mathbf{p}$和$\\sigma\n$的估计分别用$\\mathbf{\\theta}$表示为
%[text]{"align":"center"} $\\hat{\\sigma}(\\mathbf{\\theta}) = \\frac{1}{M\-K} \\mathrm{tr}\\left[ \\left(\\mathbf{I}\_M \- \\mathbf{A}(\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{\-1} \\mathbf{A}^\\mathrm{H} \\right) \\hat{\\mathbf{R}} \\right],$
%[text]{"align":"center"} $\\hat{\\mathbf{P}}(\\mathbf{\\theta}) = (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{\-1} \\mathbf{A}^\\mathrm{H} \\hat{\\mathbf{R}} \\mathbf{A} (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{\-1} \- \\hat{\\sigma} (\\mathbf{A}^\\mathrm{H} \\mathbf{A})^{\-1},\n$
%[text] 则可以得到$\\mathbf{\\theta}$的集中似然估计问题
%[text]{"align":"center"} $\\mathbf{\\theta}\_\\mathrm{CSML} = \\arg\\min\_{\\mathbf{\\theta}}\\log \\det \\left[ \\mathbf{A} \\hat{\\mathbf{P}} \\mathbf{A}^\\mathrm{H} \+ \\hat{\\sigma} \\mathbf{I}\_M \\right].$
%[text] This optimization problem is solved by MATLAB's built\-in function `fmincon`. It is adviced to provide a good initial guess to obtain good results. 
%[text] ## Syntax
%[text] `sp = mle_sto_con_1d(R, n, design, wavelength, [x0, ...])`
%[text] `sp = mle_sto_con_1d(R, n, f_steering, [], [x0, ...])`
%[text] `[sp, P, noise_var] = mle_sto_con_1d(R, n, design, wavelength, [x0, ...])`
%[text] `[sp, P, noise_var] = mle_sto_con_1d(R, n, f_steering, [], [x0, ...])`
%[text] ## Inputs
%[text] `R` — Sample covariance matrix. 
%[text] `k` — Number of sources. 
%[text] `design` — Array design. Can also be a function handle that generates a steering matrix. This function must take two arguments, wavelength and the doa vector. 
%[text] `wavelength` — Wavelength. 
%[text] `x0` — (Optional) an $k \\times1$ vector storing the initial guess in the following order: `[doas; p; noise_var]`. If absent, an initial guess will be obtained using the MVDR beamformer. 
%[text] `...` — Options: 
%[text]     `'Unit'` — Can be `'radian'`, `'degree'`, or `'sin'`. Default value is `'radian'`. 
%[text]     `'Verbose'` — If set to true, will display detailed solver outputs. 
%[text] ## Output
%[text] `sp` — Spectrum structure with the following fields: 
%[text]     `x` — Designed to be a 1 × grid\_size vector (spectrum axis). In this function, it equals `x_est`.
%[text]     `y` — Designed to be an 1 × grid\_size spectrum value vector. In this function, it is set to all ones.
%[text]     `x_est` — An $1\\times k$ vector storing the estimated DOAs. 
%[text]     `x_unit` — The same as the unit specified by `'Unit'`.
%[text]     `resolved` — Constant value `true`. 
%[text]     `discrete` — Constant value `true`. 
%[text] `P` — An $k\\times k$ matrix of estimated source power covariance. 
%[text] `noise_var` — Estimated noise power.
%[text] ## Function Definition
function [sp, P, noise_var] = mle_sto_con_1d(R, k, design, wavelength, x0, varargin)
if nargin < 5
    x0 = [];
else
    if length(x0) ~= k
        error('The vector of the initial guess has incorrect length.');
    end
end
unit = 'radian';
verbose = false;
for ii = 1:2:nargin - 5
    option_name = varargin{ii};
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
%[text] ### prepare initial guess
if isempty(x0)
    % generate initial guess
    x0 = zeros(k, 1);
    sp_mvdr = mvdr_1d(R, k, design, wavelength, 180);
    if sp_mvdr.resolved
        x0(1:k) = sp_mvdr.x_est;
    else
        warning('Failed to obtain a good initial guess.');
        % we just assume the doas is uniformly placed
        x0(1:k) = linspace(-pi/3, pi/3, k);
    end
else
    x0 = x0(:);
    % Unify to radians.
    switch (unit)
        case 'degree'
            x0(1:k) = deg2rad(x0(1:k));
        case 'sin'
            x0(1:k) = asin(x0(1:k));
        case 'radian'
        otherwise
            error('Unexpected unit ''%s''.', unit);
    end
end
%[text] ### prepare the optimization problem
if verbose
    options = optimoptions('fmincon', 'Display', 'iter');
else
    options = optimoptions('fmincon', 'Display', 'off');
end

% Set bounds for DOA search: [-pi/2, pi/2] for each DOA
lb = -pi/2 * ones(k, 1);
ub =  pi/2 * ones(k, 1);
%[text] ### solve the optimization problem
doas = fmincon(@(doas_est) cnll(R, design, wavelength, doas_est), ...
                   x0(:), [], [], [], [], lb, ub, [], options);

% Compute final estimates of A, P, and noise variance
if ishandle(design)
    A = design(wavelength, doas);
else
    A = steering_matrix(design, wavelength, doas);
end
G = A' * A;                                     % Gram matrix
Ginv = inv(G);                                  % Inverse Gram matrix

m = size(R, 1);

% Estimate signal covariance matrix P and noise variance
noise_var = real(trace((eye(m) - A * Ginv * A') * R)) / (m - k);
P = Ginv * A' * R * A * Ginv - noise_var * Ginv;

switch (unit)
    case 'degree'
        doas = rad2deg(doas);
    case 'sin'
        doas = sin(doas);
    case 'radian'
    otherwise
        error('Unexpected unit ''%s''.', unit);
end
%[text] ### store results
sp = struct();
sp.x_est = doas';
sp.x = sp.x_est;
sp.x_unit = unit;
sp.y = ones(1, k);
sp.resolved = true;
sp.discrete = true;
end
%%
%[text] # Concentrated negative log\-likelihood for fixed DOAs
%[text] This function evaluates the concentrated negative log\-likelihood cost given DOAs theta, using closed\-form estimates of signal covariance and noise variance.
function obj = cnll(R, design, wavelength, doas)
    if ishandle(design)
        A = design(wavelength, doas);
    else
        A = steering_matrix(design, wavelength, doas);
    end

    m = design.element_count;
    k = length(doas);
    
    G = A' * A;
    if rcond(G) < 1e-10  % 检查 Gram 矩阵奇异性
        obj = Inf; return;
    end
    Ginv = inv(G);

    noise_var_est = real(trace((eye(m) - A * Ginv * A') * R)) / (m - k);
    if noise_var_est <= 0 || ~isfinite(noise_var_est)
        obj = Inf; return;
    end

    P_est = Ginv * A' * R * A * Ginv - noise_var_est * Ginv;

    S = A * P_est * A' + noise_var_est * eye(m);
    S = 0.5 * (S + S');  % 强制 Hermitian

    % 使用 Cholesky 分解判断是否正定
    [L, p] = chol(S, 'lower');
    if p > 0
        obj = Inf; return;
    end

    logdetS = sum(log(diag(L)));  % 更稳健地计算 logdet
    if ~isfinite(logdetS)
        obj = Inf; return;
    end

    obj = logdetS;
end


%[appendix]
%---
%[metadata:styles]
%   data: {"code":{"fontFamily":"consolaslxgw"}}
%---
