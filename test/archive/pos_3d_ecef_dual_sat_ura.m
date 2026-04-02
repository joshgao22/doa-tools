%[text] %[text:anchor:3998] # Perform position estimation in 2D\-ECEF coordinate
%[text] This script demonstrates 2D user positioning using dual\-satellite direction\-of\-arrival (DOA) estimation. Each satellite is equipped with a uniform linear array (ULA) and estimates the signal's direction based on spatial snapshots. Unlike isolated estimation, the two satellites share the same transmitted source signal, allowing more accurate localization through intersecting DOA lines in a global 2D coordinate system centered at the origin. The setup provides a practical basis for extending to ECEF\-based or 3D satellite positioning systems.
%[text] **🔧 Workflow Overview**
%[text] 1. **Parameter Initialization**: Define signal, array, and geometry parameters including wavelength, element count, and satellite/user positions.
%[text] 2. **Source Signal Generation**: Simulate a single stochastic user signal shared by both satellite receivers.
%[text] 3. **Dual DOA Estimation**: For each satellite: 1) Construct and rotate the ULA to face the origin. 2) Generate synthetic signal snapshots with the shared signal and independent noise. 3) Estimate the signal's direction using the MUSIC algorithm. 4) Convert the estimated angle into a global 2D unit direction vector.
%[text] 4. **Position Estimation**: 1) Model each DOA estimate as a line in space. 2) Solve a least\-squares system to find the intersection point of the two lines—this serves as the estimated user position.
%[text] 5. **Error Analysis**: 1) Compare the estimated position with the ground truth. 2) Report the localization error in kilometers. \
%[text] **✅ Final Output**
%[text] - **Estimated User Coordinates** in the global xOy plane.
%[text] - **Localization Error** (in kilometers), computed as the Euclidean distance between estimated and true user positions.
%[text] - This output reflects the accuracy of dual\-satellite DOA fusion and validates the shared\-signal modeling approach under noisy conditions. \
clear(); close all;
%[text] %[text:anchor:8871] ## Parameters
wavelength = 1;            % normalized wavelength
d = wavelength / 2;        % half-wavelength spacing
m = 12;                    % number of array elements
snapshot_count = 2000;     % number of signal snapshots
snr_db = 0;               % signal-to-noise ratio in dB
power_source = 1;          % source signal power
power_noise = power_source / (10^(snr_db / 10)); % noise power

% User position (ground truth)
r_usr = [0; 300e3];

% Satellite positions
r_sat1 = [-400e3; -800e3];
r_sat2 = [500e3; -700e3];
%[text] %[text:anchor:0ece] ## DOA Estimation
% Create shared source signal
S_shared = gen_ccsg(1, snapshot_count, power_source);  % k=1 source

% DOA estimation for both satellites using the shared signal
u_est1 = estimate_doa(r_sat1, r_usr, S_shared, m, d, wavelength, snapshot_count, power_noise);
u_est2 = estimate_doa(r_sat2, r_usr, S_shared, m, d, wavelength, snapshot_count, power_noise);
%[text] ## Position Estimation
%[text] 在二维或三维空间中，一条直线可以由起点$\\mathbf{r}\_0$和方向矢量$\\mathbf{u}$唯一确定，也即直线可以表示为
%[text]{"align":"center"} $\\mathbf{r}(t) = \\mathbf{r}\_0 \+ \\mathbf{u}$
%[text] 也即从起点出发，沿着方向$\\mathbf{u}$移动$t$个单位，得到直线上任意一点。那么回到双星定位问题中，用户的位置在两条方向线的交点附近，也即
%[text]{"align":"center"} $\\mathbf{r}\_\\mathrm{user} \\approx \\mathbf{r}\_1 \+ t\_1 \\mathbf{u}\_1 \\approx \\mathbf{r}\_2 \+ t\_2 \\mathbf{u}\_2$
%[text] 通过求解上述等式，可以获得用户的位置。在二维中，如果两条方向线不是平行的，那么这个方程是可以唯一解的。实际中我们会写成 $\\mathbf{Ax = b}$ 的形式
%[text]{"align":"center"} $\\left[ \\matrix{\\mathbf{u}\_1 & \-\\mathbf{u}\_2} \\right] \\left[ \\matrix{t\_1 \\cr t\_2} \\right] = \\mathbf{r}\_2 \- \\mathbf{r}\_1$
%[text] 在 MATLAB 中可以使用左除运算很方便地求解。
% Least-squares intersection point estimation
A = [u_est1, -u_est2];
b = r_sat2 - r_sat1;
t = A \ b;
r_usr_est = r_sat1 + u_est1 * t(1);

% Compute localization error
loc_error = norm(r_usr_est - r_usr);
fprintf('Estimated user position: (%.2f km, %.2f km)\n', r_usr_est(1)/1e3, r_usr_est(2)/1e3); %[output:31f098f0]
fprintf('True user position:      (%.2f km, %.2f km)\n', r_usr(1)/1e3, r_usr(2)/1e3); %[output:84d2aa76]
fprintf('Localization error: %.2f km\n', loc_error/1e3); %[output:35ff2a19]
%[text] ## Functions
function u_est = estimate_doa(r_sat, r_usr, S_shared, m, d, wavelength, snapshot_count, power_noise)
  % Create ULA and rotate toward origin
  design_ula = design_array_1d('ula', m, d);
  design_ula = array_rot_1d(design_ula, -r_sat);

  % Generate user direction and steering vector
  u = (r_sat - r_usr) / norm(r_sat - r_usr);
  wavenum = 2 * pi / wavelength;
  A = exp(1j * wavenum * design_ula.elememt_positions' * u);

  % Generate snapshots using shared source and unique noise
  X = A * S_shared + gen_ccsg(m, snapshot_count, power_noise);
  R = (X * X') / snapshot_count;

  % Estimate DOA
  sp_music = music_1d(R, 1, design_ula, wavelength, 180, 'RefineEstimates', true);

  % Convert to unit vector in global coordinate
  theta_sat = atan2(r_sat(2), r_sat(1));
  theta_user = theta_sat - sp_music.x_est;
  u_est = [cos(theta_user); sin(theta_user)];
end

function X = gen_ccsg(m, n, cov)
  X0 = randn(m, n) + 1j * randn(m, n);
  if isscalar(cov)
    X = sqrt(cov / 2) * X0;
  elseif isvector(cov)
    X = bsxfun(@times, X0, cov(:));
  else
    C = sqrtm(cov / 2);
    X = C * X0;
  end
end

function design_rot = array_rot_1d(design, n)
% Rotate a 1D array to align its normal vector to a target direction
  theta = atan2(n(2), n(1));                % compute rotation angle
  rot_mat = [cos(theta), -sin(theta);
             sin(theta),  cos(theta)];     % construct 2D rotation matrix

  % Arrange element positions along Y-axis and center at origin
  elem_pos_centroid = [zeros(1, design.element_count);
                       design.element_positions - mean(design.element_positions)];

  % Apply rotation and translate to satellite direction
  elem_pos_rot = rot_mat * elem_pos_centroid + n;

  design_rot = design;
  design_rot.elememt_positions = elem_pos_rot;  % update element positions
end

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":14.8}
%---
%[output:31f098f0]
%   data: {"dataType":"text","outputData":{"text":"Estimated user position: (0.20 km, 299.45 km)\n","truncated":false}}
%---
%[output:84d2aa76]
%   data: {"dataType":"text","outputData":{"text":"True user position:      (0.00 km, 300.00 km)\n","truncated":false}}
%---
%[output:35ff2a19]
%   data: {"dataType":"text","outputData":{"text":"Localization error: 0.58 km\n","truncated":false}}
%---
