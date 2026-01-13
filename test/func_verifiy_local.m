% compares csml and music performance of ULA with the CRB.
clear(); close all;

wavelen = 1; % normalized
d = wavelen / 2;

snr_db = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^snr_db / 10);
numSnap = 20000;

%% 1D
doas = [-pi/3, pi/6, pi/4];
doas = doas + (rand(size(doas))-0.5)*pi/3600;
numSource = length(doas);

design_ula = ula_1d(12, d);
[~, R, ~] = snapshot_gen_sto(design_ula, doas, wavelen, numSnap, pwrNoise, pwrSource);
sp_music_1d_old = music_1d(R, numSource, design_ula, wavelen, 180, 'refineestimates', true);
sp_mvdr_1d_old = mvdr_1d(R, numSource, design_ula, wavelen, 180, 'refineestimates', true);
sp_mle_sto_grid_1d_old = mle_sto_con_grid_1d(R, numSource, design_ula, wavelen, 180, 'refineestimates', true);
sp_mle_sto_opt_1d_old = mle_sto_uc_1d(R, numSource, design_ula, wavelen);

arrayUla = createUla(12, d);
[~, R, ~] = snapGenStochastic(arrayUla, doas, wavelen, numSnap, pwrNoise, pwrSource);
doagrid = genDoaGrid("angle", 1, 180);
sp_music_1d_new = estimatorMusic(arrayUla, wavelen, R, numSource, doagrid, true);
sp_mvdr_1d_new = estimatorMvdr(arrayUla, wavelen, R, numSource, doagrid, true);
sp_mle_sto_grid_1d_new = estimatorMleStoConGrid(arrayUla, wavelen, R, numSource, doagrid, true);
sp_mle_sto_opt_1d_new = estimatorMleStoConOpt(arrayUla, wavelen, R, numSource, doagrid);

%% 2D
doas = [pi/3, pi/2, 7*pi/6; pi/6, pi/4, pi/3];
doas = doas + (rand(size(doas))-0.5)*pi/3600;
numSource = size(doas, 2);

design_upa = ura_2d([4, 4], d);
[~, R, ~] = snapshot_gen_sto(design_upa, doas, wavelen, numSnap, pwrNoise, pwrSource);
sp_music_2d_old = music_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mvdr_2d_old = mvdr_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mle_sto_grid_2d_old = mle_sto_con_grid_2d(R, numSource, design_upa, wavelen, [180 45], 'refineestimates', true);
sp_mle_sto_opt_2d_old = mle_sto_uc_2d(R, numSource, design_upa, wavelen);

arrayUpa = createUpa([4, 4], d);
[~, R, ~] = snapGenStochastic(arrayUpa, doas, wavelen, numSnap, pwrNoise, pwrSource);
doagrid = genDoaGrid("angle", 2, [180 45]);
sp_music_2d_new = estimatorMusic(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mvdr_2d_new = estimatorMvdr(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mle_sto_grid_2d_new = estimatorMleStoConGrid(arrayUpa, wavelen, R, numSource, doagrid, true);
sp_mle_sto_opt_2d_new = estimatorMleStoConOpt(arrayUpa, wavelen, R, numSource, doagrid);

