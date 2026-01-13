%[text] # Perform DOA estimation in 2D coordinate
%[text] This model provides a 2D basis for DOA analysis. Instead of using a local array-centered coordinate system, it operates in a global coordinate system centered at the origin. The array is rotated to face the origin, and the direction of arrival is represented as a 2D unit vector in the $xOy$ plane.
clear(); close all;
%[text] ## Parameters
wavelength = 1;     % normalized wavelength
d = wavelength / 2; % half-wavelength spacing
m = 12;             % number of array elements
snapshot_count = 2000;  % number of signal snapshots
snr_db = 10;            % signal-to-noise ratio in dB
power_source = 1;       % source signal power
power_noise = power_source / (10^(snr_db / 10)); % noise power from SNR

% Set positions (in meters)
doas_deg = [30 40];
doas_rad = doas_deg/180*pi;
source_count = length(doas_deg); % number of signal sources
%[text] ## Array formulation
% design a ULA with m elements
design_ula = design_array_1d('ula', m, d);
%[text] ## Signal Generation & DOA estimation
% Stochastic signal model
[~, R, S] = snapshot_gen_sto(design_ula, doas_rad, wavelength, snapshot_count, power_noise, power_source);

% Perform 1D music with 2D array
sp_music = music_1d(R, source_count, design_ula, wavelength, 180, 'RefineEstimates', true);
sp_mle = mle_sto_con_grid_1d(R, source_count, design_ula, wavelength, 180, 'RefineEstimates', true);

% Compute estimated angles
[doas_est, ~] = sort(sp_music.x_est, 'ascend');

% Compute the angular error in degrees
angle_error_deg = (doas_est - doas_rad)/pi*180;

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":14.8}
%---
