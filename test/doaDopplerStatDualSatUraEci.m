clear(); close all;
%[text] ## Parameters
% parameters
numUsr = 1;
numSym = 32;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 18e9;
wavelen = 299792458 / carrierFreq;
rng(253)

elemSpace = wavelen / 2;  % half-wavelength spacing
numElem = [4 4];          % number of array elements
snr_db = 10;              % signal-to-noise ratio in dB
pwrSource = 1;            % source signal power
pwrNoise = pwrSource / (10^(snr_db / 10)); % noise power from SNR
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]',[37.58, 37.51, 0]'];
usrLla = usrLla(:,1:numUsr);
data = [2026, 03, 12, 08, 53, 0];
utc = datetime(data);
tle = tleread("./tle/starlink_20260312.tle");
arrUpa = createUpa(numElem, elemSpace);
%%
%[text] ## Multi-Satellite Scene
% satellite scene
scene = genMultiSatScene(utc, tle, usrLla, [], [], arrUpa);
linkParam = getLinkParam(scene, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
%%
%[text] ## Signal Generation
[pilotSym, pilotInfo] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);

% pulse shaping
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;

[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

% rx signal generation
modelOpt = struct();
modelOpt.delayMode = 'phaseonly';
modelOpt.carrierPhaseMode = 'none';

[rxSig, sampleCov, cleanSig, noiseSig, meta] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, [], modelOpt);
%%
%[text] ## DOA Estimation
% s-t music
doaGrid1 = genDoaGrid("latlon", 2, [50 50], [usrLla(1)-5, usrLla(1)+5; usrLla(2)-5, usrLla(2)+5], 'eci', datevec(utc), scene.satPosEci(:,1), scene.rotMat{1}, E);
doaGrid2 = genDoaGrid("latlon", 2, [50 50], [usrLla(1)-5, usrLla(1)+5; usrLla(2)-5, usrLla(2)+5], 'eci', datevec(utc), scene.satPosEci(:,2), scene.rotMat{2}, E);
fdRange = [0 2e5];

% mle
[ddEstResult, ddPathGain, ddNoiseVar] = estimatorDoaDopplerMlePilotOpt(scene, ...
  rxSig, pilotWave, carrierFreq, sampleRate, {doaGrid1, doaGrid2}, fdRange);

[estResult, pathGain, noiseVar] = estimatorDoaMlePilotOpt(scene.array, ...
  wavelen, rxSig, pilotWave, {doaGrid1, doaGrid2});
doa1_local = scene.localDoa(:,1);
doa2_local = scene.localDoa(:,2);

[crb, aux] = crbPilotDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
  usrLla(1:2), linkParam.ref.fdGeom, 1, pwrNoise);

ddEstResult.doaParamEst %[output:36d99c43]
ddEstResult.fdRefEst %[output:9a2c33b1]
estResult.doaParamEst %[output:7973bdae]

%[appendix]{"version":"1.0"}
%---
%[metadata:styles]
%   data: {"code":{"fontFamily":"consolaslxgw"}}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":14.8}
%---
%[output:36d99c43]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":2,"type":"double","value":[["37.7794"],["36.5929"]]}}
%---
%[output:9a2c33b1]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"1.0487e+05"}}
%---
%[output:7973bdae]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":2,"type":"double","value":[["37.7795"],["36.5942"]]}}
%---
