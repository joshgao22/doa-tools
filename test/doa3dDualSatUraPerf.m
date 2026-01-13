%[text] # Perform DOA estimation in 3D-ECEF coordinate
%[text] This script performs 3D direction-of-arrival (DOA) estimation in the Earth-Centered Earth-Fixed (ECEF) coordinate system using a 2D uniform rectangular array (URA) mounted on a satellite.
%[text] **Key Features:**
%[text] - Operates fully in the global (ECEF) coordinate system without transforming to a local frame.
%[text] - The array is rotated to point toward the origin to align with the ECEF geometry.
%[text] - The direction of arrival (DOA) is represented as a unit vector in 3D space.
%[text] - MUSIC algorithm is used for 2D DOA estimation (azimuth & elevation) in the local frame.
%[text] - The estimated DOA is converted back to the global frame to estimate user position. \
%[text] **Inputs:**
%[text] - Satellite position (r\_sat) and user position (r\_usr) in meters (ECEF frame)
%[text] - Array geometry (URA), wavelength, SNR, and snapshot count \
%[text] **Outputs:**
%[text] - Estimated direction vector u\_est
%[text] - Estimated user position r\_usr\_est (based on known z\_usr)
%[text] - Estimation error norm(r\_usr\_est - r\_usr) \
%[text] *Note:* The user's z-coordinate is assumed to be known. The full 3D position is then computed by solving a line-plane intersection along the estimated direction.
clear(); close all;
%[text] ## Parameters
wavelen = 1;              % normalized wavelength
elemSpace = wavelen / 2;  % half-wavelength spacing
size = [4 4];             % number of array elements
numSnap = 2000;           % number of signal snapshots
snr_db = 10;              % signal-to-noise ratio in dB
pwrSource = 1;            % source signal power
pwrNoise = pwrSource / (10^(snr_db / 10)); % noise power from SNR
E = wgs84Ellipsoid('meters');
numSource = 1;            % number of signal sources

% Positions
lat = [30   31   28];
lon = [120  123  134];
h   = [0    550e3 550e3];
[x, y, z] = geodetic2ecef(E, lat, lon, h);

posUsr1 = [x(1); y(1); z(1)];
posSat1 = [x(2); y(2); z(2)];
posSat2 = [x(3); y(3); z(3)];

% direction vector
rotMat1_l2g = get_rotation_mat(-posSat1);
rotMat2_l2g = get_rotation_mat(-posSat2);

doa1_local = ecef2local_doa(posUsr1, posSat1, rotMat1_l2g);
doa2_local = ecef2local_doa(posUsr1, posSat2, rotMat2_l2g);

% design a URA with m elements
arrUpa = createUpa(size, elemSpace);
%%
%[text] ## DOA estimation
% Signal strength
snr_db = -20:2:10; % signal-to-noise ratio in dB

numParam = length(snr_db);
numRepeat = 2000; % 200 Monte Carlo runs for each parameter

% outcome preparation
mseMusSingl = zeros(numParam, 1);
mseMusMulti = zeros(numParam, 1);
mseMusJoint = zeros(numParam, 1);
mseMleSingl = zeros(numParam, 1);
mseMleMulti = zeros(numParam, 1);
mseMleJoint = zeros(numParam, 1);

progressbar('reset', numParam);

doaGrid1 = genDoaGrid("latlon", 2, [50 50], [min(lat)-5, max(lat)+5; min(lon)-5, max(lon)+5], posSat1, E);
doaGrid2 = genDoaGrid("latlon", 2, [50 50], [min(lat)-5, max(lat)+5; min(lon)-5, max(lon)+5], posSat2, E);

for ii = 1:numParam
  pwrNoise = pwrSource / (10^(snr_db(ii) / 10));

  currSeMusSingl = 0;
  currSeMusMulti = 0;
  currSeMusJoint = 0;
  currSeMleSingl = 0;
  currSeMleMulti = 0;
  currSeMleJoint = 0;

  for rr = 1:numRepeat
    % Generate signal snapshots with noise
    [~, R1, ~] = snapGenStochastic(arrUpa, doa1_local, wavelen, numSnap, pwrNoise, pwrSource);
    [~, R2, ~] = snapGenStochastic(arrUpa, doa2_local, wavelen, numSnap, pwrNoise, pwrSource);

    % Perform doa estimation
    spMusArr1 = estimatorMusic(arrUpa, wavelen, R1, numSource, doaGrid1, true);
    spMusArr2 = estimatorMusic(arrUpa, wavelen, R2, numSource, doaGrid2, true);
    spMusJoin = estimatorMusic({arrUpa, arrUpa}, wavelen, {R1, R2}, numSource, {doaGrid1, doaGrid2}, true);

    spMleArr1 = estimatorMleStoConOpt(arrUpa, wavelen, R1, numSource, doaGrid1);
    spMleArr2 = estimatorMleStoConOpt(arrUpa, wavelen, R2, numSource, doaGrid2);
    spMleJoin = estimatorMleStoConOpt({arrUpa, arrUpa}, wavelen, {R1, R2}, numSource, {doaGrid1, doaGrid2});

    latlonEstMusSingl = spMusArr1.latlonEst;
    latlonEstMusMulti = (spMusArr1.latlonEst + spMusArr2.latlonEst) / 2;
    latlonEstMusJoint = spMusJoin.latlonEst;
    latlonEstMleSingl = spMleArr1.latlonEst;
    latlonEstMleMulti = (spMleArr1.latlonEst + spMleArr2.latlonEst) / 2;
    latlonEstMleJoint = spMleJoin.latlonEst;

    currSeMusSingl = currSeMusSingl + sum((posUsr1(1:2) - latlonEstMusSingl).^2);
    currSeMusMulti = currSeMusMulti + sum((posUsr1(1:2) - latlonEstMusMulti).^2);
    currSeMusJoint = currSeMusJoint + sum((posUsr1(1:2) - latlonEstMusJoint).^2);
    currSeMleSingl = currSeMleSingl + sum((posUsr1(1:2) - latlonEstMleSingl).^2);
    currSeMleMulti = currSeMleMulti + sum((posUsr1(1:2) - latlonEstMleMulti).^2);
    currSeMleJoint = currSeMleJoint + sum((posUsr1(1:2) - latlonEstMleJoint).^2);
  end

  mseMusSingl = currSeMusSingl / (numSource * numRepeat);
  mseMusMulti = currSeMusMulti / (numSource * numRepeat);
  mseMusJoint = currSeMusJoint / (numSource * numRepeat);
  mseMleSingl = currSeMleSingl / (numSource * numRepeat);
  mseMleMulti = currSeMleMulti / (numSource * numRepeat);
  mseMleJoint = currSeMleJoint / (numSource * numRepeat);

  progressbar('advance');
end
%%
% semilogy( ...
%   snr_db, mleMusSingl, '-o', snr_db, mleMusMulti, '--o', ...
%   snr_db, mleMleSingl, '-x', snr_db, mleMleMulti, '--x', ...
%   snr_db, mleMleJoint, '-.x');
% xlabel('SNR (dB)'); ylabel('MSE / rad^2'); grid on;
% legend({'MUSIC-Single', 'MUSIC-Multi', 'CSML-Single', 'CSML-Multi', 'CSML-Joint'}, Location="southwest");
% title(['DoA Est with UPA (Sinlge Source, ', num2str(numSnap), ' snaps)']);

%[appendix]{"version":"1.0"}
%---
%[metadata:styles]
%   data: {"code":{"fontFamily":"consolaslxgw"}}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":14.8}
%---
