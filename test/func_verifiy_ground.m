% compares csml and music performance of ULA with the CRB.
clear(); close all;

wavelen = 1; % normalized
d = wavelen / 2;

snr_db = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^snr_db / 10);
numSnap = 20000;
E = wgs84Ellipsoid('meters');

% Set position
h = 550e3;
lat_sat = 35; 
lon_sat = 120;  
[x, y, z] = geodetic2ecef(E, lat_sat, lon_sat, h);
posSat = [x; y; z];

lat_usr = lat_sat + 10;   % degrees
lon_usr = lon_sat + 10;  % degrees
h_usr   = 0;    % m
[xu, yu, zu] = geodetic2ecef(E, lat_usr, lon_usr, h_usr);
posUsr1 = [xu; yu; zu];

lat_usr = lat_sat - 15;   % degrees
lon_usr = lon_sat + 15;  % degrees
h_usr   = 0;    % m
[xu, yu, zu] = geodetic2ecef(E, lat_usr, lon_usr, h_usr);
posUsr2 = [xu; yu; zu];

lat_usr = lat_sat - 10;   % degrees
lon_usr = lon_sat - 10;  % degrees
h_usr   = 0;    % m
[xu, yu, zu] = geodetic2ecef(E, lat_usr, lon_usr, h_usr);
posUsr3 = [xu; yu; zu];

% Compute the local doa (relative to the satellite)
rotMat_l2g = get_rotation_mat(-posSat);

doa1_local = ecef2local_doa(posUsr1, posSat, rotMat_l2g);
doa2_local = ecef2local_doa(posUsr2, posSat, rotMat_l2g);
doa3_local = ecef2local_doa(posUsr3, posSat, rotMat_l2g);

doa_local = [doa1_local, doa2_local, doa3_local];
doa_local = doa_local + (rand(size(doa_local))-0.5)*pi/3600;
numSource = size(doa_local, 2);

%% Angle based DoA Search

% design a URA with m elements
design_upa = ura_2d([4, 4], d);

[~, R, ~] = snapshot_gen_sto(design_upa, doa_local, wavelen, numSnap, pwrNoise, pwrSource);
sp_music_2d_old = music_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mvdr_2d_old = mvdr_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mle_sto_grid_2d_old = mle_sto_con_grid_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mle_sto_opt_2d_old = mle_sto_uc_2d(R, numSource, design_upa, wavelen);

%% Ground-Based DoA Search

arrayUpa = createUpa([4, 4], d);
[~, R, ~] = snapGenStochastic(arrayUpa, doa_local, wavelen, numSnap, pwrNoise, pwrSource);
doagrid = genDoaGrid("latlon", 2, [100 100], [lat_sat-30, lat_sat+30; lon_sat-30, lon_sat+30], posSat, E);
sp_music_2d_new = estimatorMusic(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mvdr_2d_new = estimatorMvdr(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mle_sto_grid_2d_new = estimatorMleStoConGrid(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mle_sto_opt_2d_new = estimatorMleStoConOpt(arrayUpa, wavelen, R, numSource, doagrid);

