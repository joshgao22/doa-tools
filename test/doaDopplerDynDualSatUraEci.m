clear(); close all;
%[text] ## Parameters
numUsr = 1;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 1e9;
wavelen = 299792458 / carrierFreq;
rng(253);

elemSpace = wavelen / 2;
numElem = [4 4];
snrDb = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);

utc0 = datetime([2026, 03, 12, 08, 53, 0], ...
    'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread("./tle/starlink_20260312.tle");
arrUpa = createUpa(numElem, elemSpace);

% dynamic frame settings
numFrame = 100;
frameIntvlSec = 1 / 750;
% refFrameIdx = round((numFrame + 1) / 2);
refFrameIdx = 1;
utcVec = utc0 + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);

% search settings
gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

% estimator settings
useFrozenSteering = true;
optVerbose = true;
%%
%[text] ## Select Two Visible Satellites
[satIdxInfo, satAccess] = findVisibleSatFromTle(utc0, tle, usrLla);
satIdx = localPickDualSat(satIdxInfo, satAccess, 1);
%%
%[text] ## Multi-Frame Scene
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "centroid", [], refFrameIdx);
refFrameIdx = sceneSeq.refFrameIdx;
sceneRef = sceneSeq.sceneCell{refFrameIdx};

if sceneSeq.numSat ~= 2
  error('doaDopplerDynDualSatUraEci:InvalidNumSat', ...
    'This script expects exactly two satellites.');
end

% frame-wise link truth
linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end

truth = localGetDynTruth(linkParamCell, sceneSeq.timeOffsetSec);
%%
%[text] ## Pilot Waveform
[pilotSym, pilotInfo] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);

pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;

[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
%%
%[text] ## Dynamic Snapshot Generation
snapOpt = struct();
snapOpt.delayMode = 'phaseonly';
snapOpt.carrierPhaseMode = 'none';

simOpt = struct();
simOpt.steeringMode = 'framewise';      % use true frame-wise geometry
simOpt.accessMode = 'framewise';
simOpt.linkParamCell = linkParamCell;
simOpt.snapshotModelOpt = snapOpt;

% optional: independent nuisance phase per frame / per satellite
pathGainCell = localBuildFramePathGain(sceneSeq.numFrame, sceneSeq.numSat, sceneSeq.numUser);

[rxSigCell, sampleCovCell, cleanSigCell, noiseSigCell, meta] = ...
  genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGainCell, simOpt);
%%
%[text] ### Reference-Frame DoA Grid
doaGridCell = cell(1, sceneSeq.numSat);
for iSat = 1:sceneSeq.numSat
  doaGridCell{iSat} = genDoaGrid("latlon", 2, gridSize, searchRange, ...
    'eci', datevec(sceneSeq.utcRef), sceneRef.satPosEci(:, iSat), ...
    sceneRef.rotMat{iSat}, E);
end
%%
%[text] ### Dynamic Joint DoA-Doppler Estimation
fdRange = [1e3 1e4];
fdRateRange = [-200 0];

estOpt = struct();
estOpt.steeringMode = 'framewise';
if useFrozenSteering
  estOpt.steeringMode = 'frozenRef';
end
estOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
estOpt.initFdCount = 81;
estOpt.useLogObjective = true;

[ddEstResult, ddPathGain, ddNoiseVar] = estimatorDoaDopplerMlePilotDynOpt( ... %[output:group:5c591c6f] %[output:719dcf89]
  sceneSeq, rxSigCell, pilotWave, carrierFreq, waveInfo.sampleRate, ... %[output:719dcf89]
  doaGridCell, fdRange, fdRateRange, [], optVerbose, estOpt); %[output:group:5c591c6f] %[output:719dcf89]
%%
%[text] ### Reference-Frame Static Baseline
statFdRange = fdRange;
[statEstResult, statPathGain, statNoiseVar] = estimatorDoaDopplerMlePilotOpt( ...
  sceneRef, rxSigCell{refFrameIdx}, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  doaGridCell, statFdRange);
%%
%[text] ### Display Truth and Estimates
latlonTrue = usrLla(1:2, 1);
latlonDynEst = ddEstResult.latlonEst(:);
latlonStatEst = statEstResult.latlonEst(:);

fdRefTrue = truth.fdRefFit;
fdRateTrue = truth.fdRateFit;

fdRefDynEst = ddEstResult.fdRefEst;
fdRateDynEst = ddEstResult.fdRateEst;

disp('========== Dynamic truth / estimate =========='); %[output:6924ce25]
disp(table(latlonTrue(1), latlonDynEst(1), latlonStatEst(1), ... %[output:group:1e453f80] %[output:80981154]
  'VariableNames', {'latTrueDeg', 'latDynEstDeg', 'latStatEstDeg'})); %[output:group:1e453f80] %[output:80981154]
disp(table(latlonTrue(2), latlonDynEst(2), latlonStatEst(2), ... %[output:group:92a9ccf3] %[output:75cb5ef4]
  'VariableNames', {'lonTrueDeg', 'lonDynEstDeg', 'lonStatEstDeg'})); %[output:group:92a9ccf3] %[output:75cb5ef4]
disp(table(fdRefTrue, fdRefDynEst, fdRateTrue, fdRateDynEst, ... %[output:group:891a465d] %[output:46f10a04]
  'VariableNames', {'fdRefTrueHz', 'fdRefDynEstHz', 'fdRateTrueHzPerSec', 'fdRateDynEstHzPerSec'})); %[output:group:891a465d] %[output:46f10a04]
%%
%[text] ### Plot Doppler Trajectories
timeMs = 1e3 * sceneSeq.timeOffsetSec(:).';
fdRefFitTrue = truth.fdRefFit + truth.fdRateFit * sceneSeq.timeOffsetSec(:).';
fdRefFitEst = ddEstResult.fdRefEst + ddEstResult.fdRateEst * sceneSeq.timeOffsetSec(:).';

figure(); %[output:2b9b03cf]
tiledlayout(1, 3); %[output:2b9b03cf]

nexttile(); %[output:2b9b03cf]
plot(timeMs, truth.refFd, '-o', 'LineWidth', 1.2); hold on; %[output:2b9b03cf]
plot(timeMs, fdRefFitTrue, '--', 'LineWidth', 1.2); %[output:2b9b03cf]
plot(timeMs, fdRefFitEst, '-x', 'LineWidth', 1.2); %[output:2b9b03cf]
grid on; %[output:2b9b03cf]
xlabel('Time offset (ms)'); %[output:2b9b03cf]
ylabel('Doppler (Hz)'); %[output:2b9b03cf]
legend({'ref fd true', 'ref fd fitted true', 'ref fd fitted est'}, 'Location', 'best'); %[output:2b9b03cf]
title('Reference Doppler'); %[output:2b9b03cf]

nexttile(); %[output:2b9b03cf]
plot(timeMs, truth.deltaFd(1, :), '-o', 'LineWidth', 1.2); hold on; %[output:2b9b03cf]
plot(timeMs, ddEstResult.deltaFdEst(1, :), '--x', 'LineWidth', 1.2); %[output:2b9b03cf]
plot(timeMs, truth.deltaFd(2, :), '-s', 'LineWidth', 1.2); %[output:2b9b03cf]
plot(timeMs, ddEstResult.deltaFdEst(2, :), '--d', 'LineWidth', 1.2); %[output:2b9b03cf]
grid on; %[output:2b9b03cf]
xlabel('Time offset (ms)'); %[output:2b9b03cf]
ylabel('\Delta Doppler (Hz)'); %[output:2b9b03cf]
legend({'\Delta fd sat1 true', '\Delta fd sat1 est', ... %[output:2b9b03cf]
        '\Delta fd sat2 true', '\Delta fd sat2 est'}, ... %[output:2b9b03cf]
        'Location', 'best'); %[output:2b9b03cf]
title('Differential Doppler'); %[output:2b9b03cf]

nexttile(); %[output:2b9b03cf]
plot(timeMs, truth.satFd(1, :), '-o', 'LineWidth', 1.2); hold on; %[output:2b9b03cf]
plot(timeMs, ddEstResult.fdEst(1, :), '--x', 'LineWidth', 1.2); %[output:2b9b03cf]
plot(timeMs, truth.satFd(2, :), '-s', 'LineWidth', 1.2); %[output:2b9b03cf]
plot(timeMs, ddEstResult.fdEst(2, :), '--d', 'LineWidth', 1.2); %[output:2b9b03cf]
grid on; %[output:2b9b03cf]
xlabel('Time offset (ms)'); %[output:2b9b03cf]
ylabel('Doppler (Hz)'); %[output:2b9b03cf]
legend({'fd sat1 true', 'fd sat1 est', 'fd sat2 true', 'fd sat2 est'}, ... %[output:2b9b03cf]
  'Location', 'best'); %[output:2b9b03cf]
title('Per-Satellite Doppler'); %[output:2b9b03cf]
%%
%[text] ### Optional snapshot save
% saveExpSnapshot("doaDopplerDynDualSatUraEci");
%%
function satIdx = localPickDualSat(satIdxInfo, satAccess, userIdx)
%LOCALPICKDUALSAT Pick two available satellites for one user.

if nargin < 3 || isempty(userIdx)
  userIdx = 1;
end

if iscell(satIdxInfo.available)
  satIdx = satIdxInfo.available{userIdx};
else
  satIdx = satIdxInfo.available;
end

if numel(satIdx) < 2
  error('doaDopplerDynDualSatUraEci:NotEnoughVisibleSat', ...
    'Need at least two available satellites at the reference epoch.');
end

usrElev = satAccess.usrElevationDeg(satIdx, userIdx);
[~, sortIdx] = sort(usrElev, 'descend');
satIdx = satIdx(sortIdx(1:2));
end

function truth = localGetDynTruth(linkParamCell, timeOffsetSec)
%LOCALGETDYNTRUTH Collect dynamic Doppler truth and fit fdRef / fdRate.

numFrame = numel(linkParamCell);
numSat = linkParamCell{1}.numSat;

refFd = zeros(1, numFrame);
satFd = zeros(numSat, numFrame);
deltaFd = zeros(numSat, numFrame);

for iFrame = 1:numFrame
  linkParam = linkParamCell{iFrame};
  refFd(iFrame) = linkParam.ref.fdGeom(1);
  satFd(:, iFrame) = linkParam.fdGeom(:, 1);
  deltaFd(:, iFrame) = linkParam.ref.deltaFdGeom(:, 1);
end

[fdRefFit, fdRateFit] = localFitFdLine(timeOffsetSec, refFd);

truth = struct();
truth.refFd = refFd;
truth.satFd = satFd;
truth.deltaFd = deltaFd;
truth.fdRefFit = fdRefFit;
truth.fdRateFit = fdRateFit;
end

function [fdRefFit, fdRateFit] = localFitFdLine(timeOffsetSec, fdRefSeq)
%LOCALFITFDLINE Fit fdRefSeq ~= fdRefFit + fdRateFit * timeOffsetSec.

timeVec = timeOffsetSec(:);
fdVec = fdRefSeq(:);

if numel(timeVec) <= 1 || max(timeVec) - min(timeVec) <= 0
  fdRefFit = fdVec(1);
  fdRateFit = 0;
  return;
end

designMat = [ones(numel(timeVec), 1), timeVec];
coef = designMat \ fdVec;

fdRefFit = coef(1);
fdRateFit = coef(2);
end

function pathGainCell = localBuildFramePathGain(numFrame, numSat, numUser)
%LOCALBUILDFRAMEPATHGAIN Generate block-wise complex gains.

pathGainCell = cell(1, numFrame);
for iFrame = 1:numFrame
  phaseMat = 2 * pi * rand(numSat, numUser);
  pathGainCell{iFrame} = exp(1j * phaseMat);
end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:719dcf89]
%   data: {"dataType":"text","outputData":{"text":"初始点 X0 不在边界 LB 和 UB 之间；\nFMINCON 对 X0 进行了偏移以严格满足边界条件。\n\n                                            First-order      Norm of\n Iter F-count            f(x)  Feasibility   optimality         step\n    0       5   -1.474160e+07    0.000e+00    4.596e+06\n    1      18   -1.504099e+07    0.000e+00    1.477e+06    1.203e-01\n    2      26   -1.508746e+07    0.000e+00    3.368e+05    6.014e-02\n    3      31   -1.508771e+07    0.000e+00    3.521e+05    1.241e+00\n    4      36   -1.508888e+07    0.000e+00    4.743e+02    4.142e-01\n    5      41   -1.508888e+07    0.000e+00    6.168e+01    1.000e+00\n    6      46   -1.508888e+07    0.000e+00    5.821e+02    6.274e+00\n    7      51   -1.508890e+07    0.000e+00    3.612e+02    1.197e+01\n    8      56   -1.508895e+07    0.000e+00    1.590e+03    9.230e+01\n    9      61   -1.508895e+07    0.000e+00    7.219e+01    3.302e+01\n   10      66   -1.508895e+07    0.000e+00    2.092e+02    5.416e-01\n   11      71   -1.508895e+07    0.000e+00    2.301e+01    1.709e-01\n   12      77   -1.508895e+07    0.000e+00    1.979e+03    3.534e+00\n   13      82   -1.508895e+07    0.000e+00    1.243e+03    3.463e-01\n   14      87   -1.508895e+07    0.000e+00    7.167e+02    7.238e-01\n   15      92   -1.508895e+07    0.000e+00    3.182e+01    2.101e-01\n   16      97   -1.508895e+07    0.000e+00    3.256e+02    2.767e-01\n   17     106   -1.508895e+07    0.000e+00    1.547e+02    6.929e-06\n   18     115   -1.508895e+07    0.000e+00    5.470e-02    6.749e-06\n   19     135   -1.508895e+07    0.000e+00    1.036e+01    1.352e-03\n   20     140   -1.508895e+07    0.000e+00    3.164e+00    9.452e-02\n   21     154   -1.508895e+07    0.000e+00    2.201e+00    2.119e-05\n   22     160   -1.508895e+07    0.000e+00    2.164e+00    2.169e-05\n   23     169   -1.508895e+07    0.000e+00    1.592e+00    1.769e-07\n   24     174   -1.508895e+07    0.000e+00    6.533e-02    5.609e-04\n   25     183   -1.508895e+07    0.000e+00    4.777e+00    1.119e-03\n   26     188   -1.508895e+07    0.000e+00    1.467e+00    9.573e-04\n   27     193   -1.508895e+07    0.000e+00    8.927e-01    5.401e-04\n   28     198   -1.508895e+07    0.000e+00    1.004e+00    3.029e-03\n   29     211   -1.508895e+07    0.000e+00    5.395e-01    4.618e-06\n   30     219   -1.508895e+07    0.000e+00    5.361e-01    1.182e-06\n\n                                            First-order      Norm of\n Iter F-count            f(x)  Feasibility   optimality         step\n   31     229   -1.508895e+07    0.000e+00    1.204e+00    5.523e-08\n   32     237   -1.508895e+07    0.000e+00    1.204e+00    1.004e-07\n\n<a href = \"matlab: helpview('optim','feasible_better_objective','CSHelpWindow');\">找到具有较低目标函数值的可行点，但不满足最优性条件。请参阅 output.bestfeasible。<\/a>\n\n\n<a href = \"matlab: helpview('optim','local_min_poss_with_constr','CSHelpWindow');\">可能存在局部最小值。满足约束<\/a>。\n\nfmincon 已停止，因为<a href = \"matlab: helpview('optim','norm_curr_step_simple_fminconip','CSHelpWindow');\">当前步长<\/a>小于\n<a href = \"matlab: helpview('optim','step_size_tol','CSHelpWindow');\">步长容差<\/a>值并且在<a href = \"matlab: helpview('optim','constraint_tolerance','CSHelpWindow');\">约束容差<\/a>值范围内满足约束。\n\n<<a href = \"matlab: createExitMsg({'optimlib:sqpLineSearch:Exit2basic','fmincon'},{'optimlib:sqpLineSearch:Exit2detailed','1.000000e-10','0.000000e+00','1.000000e-06'},true,true);;\">停止条件详细信息<\/a>>\n","truncated":false}}
%---
%[output:6924ce25]
%   data: {"dataType":"text","outputData":{"text":"========== Dynamic truth \/ estimate ==========\n","truncated":false}}
%---
%[output:80981154]
%   data: {"dataType":"text","outputData":{"text":"    latTrueDeg    latDynEstDeg    latStatEstDeg\n    __________    ____________    _____________\n\n      37.78          37.779          37.779    \n\n","truncated":false}}
%---
%[output:75cb5ef4]
%   data: {"dataType":"text","outputData":{"text":"    lonTrueDeg    lonDynEstDeg    lonStatEstDeg\n    __________    ____________    _____________\n\n      36.59          36.586           36.59    \n\n","truncated":false}}
%---
%[output:46f10a04]
%   data: {"dataType":"text","outputData":{"text":"    fdRefTrueHz    fdRefDynEstHz    fdRateTrueHzPerSec    fdRateDynEstHzPerSec\n    ___________    _____________    __________________    ____________________\n\n      5826.9          5808.8             -134.43                -4.5718       \n\n","truncated":false}}
%---
%[output:2b9b03cf]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tvX+8VsV1N7oSaRVRo4fEAkGBJNDa3ltUNMVj8mIkhtTKySc0efHQeilBPVUjVkk4gDZAbvipx1SNNmgp8cZytDb6CSdNykUUXhWNDRFivVaMyK+ANYLmBwFaEu+7Np2H9czZ+9mzZ8\/az9r7rOcf8Twzs9d85zuzv8+aNbPe9c4777wD+lEEFAFFQBFQBBQBRaBECLxLBUyJRktNVQQUAUVAEVAEFIEIARUwSgRFQBFQBBQBRUARKB0CKmBKN2RqsCKgCCgCioAioAiogFEOKAKKgCKgCCgCikDpEFABU7ohU4MVAUVAEVAEFAFFQAWMckARUAQUAUVAEVAESoeACpjSDZkarAgoAoqAIqAIKAIqYJQDioAioAgoAoqAIlA6BFTAlG7I1GBFQBFQBBQBRUARUAGjHFAEFAFFQBFQBBSB0iGgAqZ0Q6YGKwKKgCKgCCgCioAKGOWAIqAIKAKKgCKgCJQOARUwpRsyNVgRUAQUAUVAEVAEVMAkcOCVV16BadOmwZ49e+pKdHR0QGdnZypzDh48CLNnz4aenh5wrZPaKHOBuD4PGTIEVq5cCSNHjmR7+tKlS2H48OEwefJktmdUqWGXcXr22Wehu7sblixZArt3765xGcfykUceKZSXdHypXf37948dFjN32tvbYezYsXVlHnroIZgzZ07d3yZOnBj1M6m9vGPfyJ68bVelfty4hFw79u\/fD9OnT4ctW7ZEkI0ePRpWrFgBLS0tqRC6ri+Um3feeWe0JrW1tUXrOHJx1KhRMHPmTJg7d26m9dBlvqZ2wqOAy1zzaFZUFRUwDQTMrFmzYNmyZTWymkmEAsZeWO1msOyCBQtg3rx5TpNMAitwotl9xknwhS98gVXEuC4wEjCSYEPWccIx3LBhQyS8m8HLrOObJmA2btxYJ1iwffyhwSViVMCksx4FjD0uoV6gZt29\/PLLaz9y8HkPPvigk4hx5V+cgKE\/qtAOXwHTjHU1FP7po9+8EipgMggYLGpPBiTJlClTolbML0Fc8MyvBfMrBL83Hh3668G8jIwZ+KsCf0ka7w3+fdWqVZFgwrKLFi2CU045JfoFTb\/Df9NfKfQZ1Btk16Hdj3sxmj7jf43nKa7PaDNiM2DAAFi3bl30S4l6ntK+Mx6YJFvN399++2148skna5g0b+o078ku42QWr0mTJkW8w88nPvEJ2Lt3L7zwwgtgeDl06NBYrsXhjZyK46Upi20uX748epYZe\/rLfPHixTBs2LCaZwjL0fbM\/DF\/T\/LA2C\/KuB8WyDdjCz7XvIjSeBjHX1vAJHE0bi67eAiax6RwT44TMGa96urqin7EJa0b+Hccl5\/+9KcwZsyYXkI0ju\/2mNsemiT+IQ8a2WG8lrYHBueR8VxSz1ISz0Kuq2Yun3zyydG6R9d2l++MsI+ztezrqgqYDALGnjR0YpkXAZLb\/NI1an3gwIGRoDGeG\/qL0bj3b7vttppXB783goF6QIwIuv7666MFmbZjL\/q4oGzfvj16ZlJ79rZQ0osxbjsC7TUvNNNnfA4KK9ymiOtzo++MgEmy1cY33NJbvpZcxgkFpFmM6b9xwaK\/IrPgnVYWkcTFEp9HvXZJW0j4kjA8p3PL8MpVwBiRbThEf53v27cvEnBmfvlw1LYnbX7SuVw+dvlZHCdg6N9wnTNeCHsup3l5qWA0P+aoleb71tbWaF205wflX6M1O+sWkt0\/yjMXAeO6rhrBZYS4zT\/8AZ30nVkDVq9eXfOQ0XeOvYb7jX7zaqmAaSBg4mJg6K85e9JSQtIXBS6iSDqzZ0t\/meB31L1oiyT66w9FAS1rTwD0zphfO6Zbjdqzt8GyvhjR62J+PWHf7rvvvtoLCf+RJKLivsOXzyWXXFIn9GjfzUQzi1Tzpkzzn5x1nJIEjC0yG+HdiEf22NiudpcYmLhn+woY20tKF3z6b1eOzpgxoy4Ogv4YaTQ\/m8+U4iyIi4GhnoJGayXyk66PSVZTzwmWiRMz+PdG\/EuzI8kDY8fAGBFG16Okraqs89VeV7du3Zr4\/kj7DvvzpS99Cb785S9DnK0mxqes66oKGEcPDH1RG7dwo0mLzZpfuihSzDaTeRzdWqLCw3aFmvIonM4999xoC4m6ZOmv7LhFoFF7dtBs0kSjk54qeZxoVIyhgKHBuLaASfvOCBgTqEf7ToPp0uKPilu2m\/Mkl3Fy8cAYAeOCdyMe2WPjKmBocLFB0rjITeBkXBCvvYVEf4HHLchZeWjmhakXJ2DiMLPnZ3PY0ZynxnlgqCWN1kp8CZt1zGxF21uRdq\/sOWCLG3ubh3rn7CBwI7SoHWlBvPbWq7Ev7sCGy3xttK7a+NA11+U7I2BM2AG1lXK7jOuqChhHAYPF7GBBujDazdBFHAUMnaC0rL1P3ChQLG5PmZ40SfLAuAaeZYmtMPuqjTwwjX752t8ZD0ySrRpIeYw1WcepkQfGFe9GvLTHxkXA2L8K83hgbO9QFg+MC0dtAZOEmT0\/myMlmvNUFwFjtrRtC9OCTWkQuqlL+YKng6hXLM0D42KHq4CJ8xKmiS3zvb0VRN8RdF21vSxZvqMemDhby76uqoDJIGAaxcBgPAkVOHQLyXbV23v0tvCgxDZH8HBfHdtJ8sCYwGETZ0MXFBprQNtz2UKy96dp\/bgYGCyP20n4seN+Gn0XFwMT9yyXBaM5S3hxT3U5hURfCq4xMGl4J\/HSjhHxETDmF7SPBybuh4U5oRIXA5OVo41iYBrNz+IY0fwnpQmYuLgUc3KM8jPuKHxckDZdl+LW17vuuqt2cjIpBsZes6kdaQIG69I+m\/WXnpQyo+IyXxutq3Ru4JodFwNjttOSRBH18FBby+7ZVgGTQcBgUfv4HnVd0j1fexE3BMVJS92bcb\/a7FMOJu6mkQfGbOfEnXRKai\/ul4Id9xN3l0NSFD9Onl\/84hewfv366Firffqj0XdJp5BMG2X\/pRDyFUO5ZNq1x8lVwCRxIw5v17I29832gX0KCV8YZmsVXe\/4se\/e8L0HBrmYdAopK0fTTiElzc+QYy69rTQBg\/YnrZVpHhisa29h2nynW1R\/\/dd\/DZs3b47ubkH+UP7Zp5DiTvSgdzlOwBghu2nTpjpx5LLdlWddRXzMM\/A6BHrvUdp3aVtzZV9XSylg7IXUDhYze5x0P5JOHpv8SYud9EVDmn1JQWxm+62ql9XFLcB0QaXBhpSHVOBRTtuXdCW1JW38y2CPcvTohX+N1lDlqCwmNxJ4LuJPVm\/CWlNKAZO0H09dddStaN+gSL0odH\/RPhEUFurqt9YXXw5GXNBfRUlBdnRbEdlAtwOTAk0bBexVn1Hhe6gcPSpgktZQ+nflaHj++bSoAiYZtVIKmKRgOTr5UMC4BNuhaDEvEvrvvnIBlc+E0jpHEcCFZceOHdG\/6ckYe2+cnkYzJ8XMZYXo5raPIdMFy967dg3I1jFSBBpxNGkNpQGiylHlkHQESilg7CN51BVv9koR+KRcGfbpIRNTEHcLpPQBVPuaj4C9\/0\/5RQMQ0VKzJ43\/xqPCeP+CHUhHXyJr166tXUgYF8zY\/N6rBWVAII6j9DixWUOpeA7J0R0bf1kGmCpt47DWkyrXv1IKGDoK9gvC\/MLFMvQUjKlju+PoxKa3RXImL6wci\/p4h1TA9HEClKD7jYJsuUX2z3b9J9x93r+XAKVqm\/hnj34AqiZiSi9g6CVWSD96xr\/RrZuGqnRP3L6S2pT5wAc+UG1mC+\/dtm3bRFsYJ2DMllLSjcwh3fPKz+bSQzo\/EZ1GAoaue5irKuQ2J4qXBz69DfC\/+mkuAu8547fhuh\/8XnONCPz0UgoYnIz4sfNe0Cv7zQvCXJFM61AM487y29mm8QUhaZGSZg\/iyWkTZ9sh5pP9cig6iFcaPtLs6ev8jBMwSWsojR3Eer6B5ovnXQ8nnHACrPqzw7Bz44EQ00zbCIDAma0DYMo\/HA\/vPnFogNaa30QpBUyje01ofIw5Rh13bwY9qpp2jFragizNnr7+goj7dWt4GHc\/i7n7hB6vppymJ5rMywfjFeLu5OHG3meJUn76oMZbx+ZoozWUHqP25eiI0\/bB47O\/DT\/oaeftmLaeGYHzJnbDx++8rhIippQCJvOI5awgbUF+7bXXYMSIETl7FbY6p03S8A+LXP7WpOHDyQVftDhtkoa\/L0Yh6+1\/eRt8\/X9o4G5ITEO2VZV4GBUwDqyQtkBxLsYOcMQW4bRJGv6+GHHVk4YPJxd8MeS0SRr+vhiFqqdxL6GQ5GunKvEwKmAcOCJtgeJcjB3gUAHjCxJTPeVnOrCcc0Ya\/ulo8Jb4zoxd8KOH3uJ9iLaeG4E\/nHwaXHbnGbnbaWYDKmAc0Je2QHEuxg5wqIDxBYmpnvIzHVjOOSMN\/3Q0+EvgsWk9ecSPs+8T0APz549+APC\/Zf6ogHEYPWkLFOdi7ACHChhfkJjqKT\/TgeWcM9LwT0eDv8Rrj70M3X92mP9B+oTMCJw88A1o+7uxlbgTpk8KGL03IzPnS1VB0pH3IoCT9gLlFAu+eHLaJA1\/X4xC1jvy5rN6CikkoAHb0lNIAcFsRlO64DQD9WKe2RfHVlqfOcWCL4s4bZKGvy9GoeuhiHnwqvfpPTChgc3RHt4Dc\/l9P4V+7x2boxU5VfusB4bjV\/rO\/Ydg1\/5DcOGHTpUzwn3Mkr74MpHWZ06x4EtnTpuk4e+LEVc9jYfhQjZbu1WJe6G9VgGTjQOxpVG4fL77JXjq1bdr35\/ZcgLcfflZKmYC4Julib74MrlwWBus6l6VBSbWsq\/vfR0GDR7E+oysjXPaNKV9Cjy9Y3VWk\/pMeUzk+PVJz8LQoXJufz1y5Aj069dP1Bhw2oR5\/v7ykWrEvaiACZgaAMXL2V95BlCwtJ8\/OPov\/m3pmtcinHuuPadwEWNu2dy0aROsXLkSTGLKpEzdeEsn5kAZO\/aYWxHbeOKJJ+DSSy8VNcnTjOlrAkYT5aUxopjvq3IxGBdaKrLTkVWRnY6RXUI9MNkxq6vRdvfzsPOtQ7D62nMi8WI+Rth85IOnwurrzsn5lGzVUajcd999MGPGDMCcUOaDKRXWrFkDn\/\/85+sajBMwSWWzWVJ86b4kYPTCsOL5lfTEqlwMxoWotHnJuaXoiyGnTdLw98VIBUxM4sGnf3xs6ycrsBPveR7azx8EU84f3KsqbimhJ+bu9rPgzNOOiZusz6AxNZin5G\/\/9m\/hySefBMxTgh+TWwdzP6FomT17NvT09ADNqWO8L1u2bAEsN27cuKge5te56KKL4LLLLqt5YGielMWLF0feGfPMr3\/96\/D888\/XxNHXvvY1mDBhQuQeNs9NytmTtd8+5as6USkWv\/nV7uh\/NVGeD0P46mCA5J8\/+kG+B5S4ZWnzklMs+A4Tp03S8PfFSAVMjIBpuemJUHiytLP\/9o\/V2kUBs2HDBsCM2ShKFixYAPPmzYOWlhYwGWYvueSShh4YFC6mnsna3d7eXreFRD0w9jOpd8cImB\/+8IeRjZgh3LaLBZSERqs6UWl39YhqkYzK9qwLr\/svGPelMdkq9YHS0uYlp1jwHU5Om6Th74uRCpgCPTCr\/nUvdP\/r68E9MDt27IiEQlyWbfSuXHXVVQ0FDHpM7r\/\/frj55pujLaa0LSQUMOaZ9vaUETCPPPIILF++vMavZnlhqjpR6cTVRHmhlj+edjQepjeu0uYlp1jwZRWnTdLw98VIBUyMgMkDJsbA4FaRHayLMTBt9zwfNb35lgvyPKKuri0mqAfGFEyLgaEeGPTcLF26NNpSokG8tgeGChjqvVm4cCFMnToVqAcmWGc9GqrqRDVQaNyLBykKrqLxMCpgfCinAiY7ahrEmx2zuhoYP4NxMHGnkDiOUlMBg4bg\/5sYGPx\/jFlJ20LCIF5a76Mf\/Shcc801vQTMtGnTojgaFDdGwOAzUPCgtwW9LB\/84AfhlltuqYuBwTKjR4+GFStWRFtbRX6qLmA0UV6RbPJ\/VhUS5fn3XgWMD3YqYLKjpgImO2a9ahhvC\/7XfFC82CeTAjxKm0hBoOoCBruvF4PJngZVvDAsL+LS5iWnWPDFitMmafj7YmTXUwETCkmA6P4XvIkXP3obb0BgMzRV1YlKIdBEeRkIUXDRKiXKCwmdtHnJKRZ8ceO0SRr+vhipgAkcAxNqILSdMAhUdaJSdPQUUhiucLRSpUR5IfGRNi85xYIvbpw2ScPfF6NKCBh6Twl2iMZb4ImaOXPmRP3EEzl43Bg\/dqwIRy6kUIOi7fgjUNWJaiOiifL8OcJVc8jv\/hu0feFmOPXjT8K7T5RzbT5Xf7O0K21ecoqFLLjUeVZfew1GjBjhW71hPWn4h+pkKbeQ8JTNzJkzYe7cubVr8hEQPDkza9YsWLZsGQwcOBCmT58eCZhRo0bVla\/qYBpShEolgBfTYTqCG2+8ETBPBx7dfvzxx+Hss8+Ojl9nSTVQVLqCqo+tPfE1HibUUpivHYx7af+H34ZTBr5RmUy\/+RCpry1tXqqACTm6zWurlAIGhcqiRYugq6ur7pQLFTYoYOJEDkItbTKFHv60Y9QuqQTi2kBhdOedd0b3zOzbty82LUFSX9LumgmFQdXH1sZJE+W5MUcT5bnhxFVK2rxUAcM10sW2W0oBQ7eJEC48OozeAfwkJSyksEqbTFmGvOhUAniMetKkSfD6669HZuL23J\/8yZ9E\/\/7nf\/7n2rFt9HZhmgK6nWe27VzSFWB73\/3ud+HVV1+FO+64Ax577DGvdAVlHtssPKBlNVFeOnKaKC8dI84S0ualChjO0S6u7VIKGAqPESwm1gXvKMH7R\/BjtpDoBW3mpUpjYPbOP3ZVfxr0g+f3TjvgWv+EP7gITvvsvF6PSKtPn1lUKgHqgUFhgvfAtLW1xXpg6EV46B3DW37xXhncysM0By7pCkwaBJOKwDddgbSFMo1PIb6X1md9OYQY1Wq1oRxNH0\/OeSMN\/3Q03EqUXsCYeI\/W1taox9u3b68F7uKLFT9G3OBLEr\/HS9h8BMxJF02Fky\/6izpk33p4ARx6cb0T2nHiBytmFTBFpBJwFTAoEk0SRwMCemHwcrtHH33UOV0B3WLKk66gqhO1EcGk9ZlzIXaaaDGFOG2Shr8vRpz1pGHEyQdfHDltkoa\/L0Z2vVIKGPprnQbuYlyG8cCYX\/0obPBXPa1T5sEsKpWAq4DBeJq4VAR2Qse0dAW2gPFNV1DmsfWd1NL6zLkQ+2LEaZM0\/H0x4qwnDSNOPvjiyGmTNPx9MaqEgLGPUdMYmLhj1HFJD8t6jLqoVAJJAga9LTt37ow8LDfccEOUagCDek0MDBIM\/7ZkyZIoJsakOUhLVzB8+HAYNmxYLZ2Bb7oCKRPV2I94UH7S4\/z075TTdhoGyulVq1bVpXzA9qX02SwunAux78LHaZM0\/H0x4qwnDSNOPvjiyGmTNPx9MaqEgMnb+aoOZl5cqlBfwtiiSDGeQPQK0qP95mQcYk1P0pntTdzuxLoo6EzGcVNu69at0N3dHYlD9DCaj4Q+U+5wLsS+HOW0SRr+vhhx1pOGEScffHHktEka\/r4YqYAR+Is11GBqOzK8EfSYPwoYKkDsLc729vbo5BZ6tsx2JwogI1RWr14NGzdujEQLemnirgaQtjhxLsS+HOe0SRr+vhhx1pOGEScffHHktEka\/r4YqYBRAROKOyLbkTJRzbblmDFjah4TKkwQPCNa8HQX\/hvFDJ6Yox6ctWvX1gLT6Yk7erJOSp91C0nklBBhlHI0fRhUwKRjpAJGBUx2lpSohoSFEreDjNdk9+7dtS0k9MYYz0poAYPtrVu3TsRIYZ+HDpV1lT6HTePHj6\/hXdaYuqIII2Fe6jbntqKGu7DnlPIUUl50pE2mvP2x65cllUCWfqP3YfPmzXDxxRc3rCZhbGkMCz3mj0HKuoWUZdTDldVft+Gw9GlJwrxUAaMCxoe74uqEnky\/+dVuOPD8F+HkC7ujvuL\/H971T\/BbA8c2JS9KWVIJZCGGffoqqW7osc1ioylLPTAoYOJycmFZDeL1QdevjgoYP9xC1ZIwL1XAqIAJxeemthN6Mv3i6Xb4r33PRoJlwDm3RuLl4Mt3wPFnfCb6\/5CfqqQSQEzo0WFMN7By5cooOSc9aowBrrfeeit88YtfjI5lxx0jpviGHlvfsXM5Rk37QrEwx9DNSSNzjJpiJLHPxiZOseA7Hpw2SeGcLzZF1JOGEScffPHktEka\/r4Y2fV0CwkAjrz5bG48f7n5i5Hn5d0nDo3+i2Km\/+\/ekLtdbKDfe8fW2qlSKgE7fYC5vA5f7BMmTKjLNF4mD0yQQc\/QiLTFiXMhzgBLXVFOm6Th74sRZz1pGHHywRdHTpuk4e+LkQqYmCDe\/atHhMKTpZ2WttfqBExVUglQLwV20HgY8N\/Tpk2DPXv21JJD4h0opt+NQK7qRC1TnzkXYt8JxmlTX+Rc1nGQhhEnH7JiU4TnUhr+vhipgIkRMHk9ML8+uBv+c+e3om0k80FPzElnh9k+sj0w5kVuX9dvnp0WA4O34xpvR0tLS2IqAJNQMUsyR2pD2jOoByaJ0KYMBsCqgIlHSdripC+HUMtzddpRjqaPJee8kYZ\/OhpuJXQLyQ2nhqUOvvw3dTEvJiaGKwaGvshpvAgaidfTX3LJJUCzORvj8W6SNWvWAOYvovXirvkvIpWAnRLCXKGPtmPCTfyYv6EHBoUXvX4\/blCqOlHVA5NvourLIR9+eWtLm5ecfPDFitMmafj7YqQeGIZ7YNCDc3jXt+oCdvFU0vFn\/GlTTiGFIkcZ26nqRFUBk4+N+nLIh1\/e2tLmJScffLHitEka\/r4YqYBhEDChBkPbyY9AVSeqCph83NCXQz788taWNi85+eCLFadN0vD3xUgFjAqYUNwR2U5VJ6oKmHx005dDPvzy1pY2Lzn54IsVp03S8PfFSAWMCphQ3BHZTlUnqgqYfHTTl0M+\/PLWljYvOfngixWnTdLw98VIBcx\/C5hQAGo78hDoa3lppC1OnAuxL9s4bZKGvy9GnPWkYcTJB18cOW2Shr8vRipgPJBrNPhtdz8PT736NrSfPwi6\/\/V1mHhgDczbvxS6P\/UEdP3wndrT7m4\/KyoT4sNJdF\/7OG2q6uTzxdquJw0fTi74YsZpkzT8fTHirCcNI04++OLIaZM0\/H0xUgHjgVyjwX\/6x2\/DxHuer7XaOWEEfOQbF8DgI6\/DeWc8Xvv7Rz54Kqy+7hyPp\/euwkl0XwM5barq5PPFWgVMduSUn9kxC1lD2hzm5IMvbpw2ScPfFyMVMB7IpQ2+8cKYpo0X5psXrIQ7dg+rPbHn2nPgwg+d6mFBfRVOovsax2lTGv6+NlelnjR8OLngO2acNknD3xcjznrSMOLkgy+OnDZJw98Xo0oImKQL0PBWWZP4Djva0dEBnZ2dUZ\/pxW1JSfHwanv8mDoGrLTBx62j67pfqmGLW0UzH\/l9+M6ACTC\/pbO2vRTKC8NJdF9icdqUhr+vzVWpJw0fTi74jhmnTdLw98WIs540jDj54Isjp03S8PfFqBICBm+JnTlzJsydO7cu4R\/eNDtr1ixYtmwZDBw4EKZPnx6JkVGjRtWVR5Hz4IMPwooVKwBFDxU4VPS4Chgs13LTE3XYzt+\/FC47sAbaBq+CPf2Oxb5svuUCOLPlhFzjx0l0X8M4barq5PPF2q4nDR9OLvhixmmTNPx9MeKsJw0jTj744shpkzT8fTGqhIBBobJo0SLo6uqqCRDsGBU2KGDiRA6Ws+ubnEIjR46EAwcOZPbAYJt2MO+nB7wMN\/\/7NfAvwzvgll9PruEeIpiXk+i+xOK0qaqTzxdrFTDZkVN+ZscsZA1pc5iTD764cdokDX9fjCohYOg2EXaI5sdBMYKeF\/xQDwvtONbfvn17Tajg1tG4ceOiZIH071k8MHHBvB\/45qdhzOHNwYN5OYnuSyxOm6o6+XyxVgGTHTnlZ3bMQtaQNoc5+eCLG6dN0vD3xagSAoZ2wggWE7eCYgSFC37MFtLYsWNrVTAWpru7G5YsWQL9+\/ePYmM2bNgQiRlb2FABg\/9et25dQ9yvfuR12PSTQ7UyJpj3O+Pug\/nbPlj7+72TBsGY9\/tvI+3evRuGDh0aigNB2uGwafz48TXb+trdLlkGRdrixLkQZ8GFluW0SRr+vhhx1pOGEScffHHktEka\/r4YVU7AmIDe1tbWqG+2ZwX\/ZsRNnEBBwWMyHxtw7DgY18FPC+Y17ecN5uUkui+xOG1yxd\/X9rLXk4YPJxd8x4rTJmn4+2LEWU8aRpx88MWR0yZp+PtiVAkBg0IEP5MnT47iWUzg7r59+8B4YNC7Mnv2bEBhg+VonSTwGnlgXD0AdjCveRbeD\/P0j9+KLr3DT54j1ZxE9yUWp01VnXy+WNv1pOHDyQVfzDhtkoa\/L0ac9aRhxMkHXxw5bZKGvy9GlRAw9jFqGgMTd4waRc60adNgz549tf6PHj26V4xMCAFD74TB00YoXPCINR6tfvrVt2Hn\/qNbTHm8MJxE9yUWp01VnXy+WKuAyY6c8jM7ZiFrSJvDnHzwxY3TJmn4+2JUCQETqvOu7WQZfBrMiyKl\/cOD6+6Ioc\/0PVLNSXRXTOxynDZlwd\/X\/jLXk4YPJxd8x4nTJmn4+2LEWU8aRpx88MWR0yZp+PtipALGA7msg2\/fzGseafIlmf\/3PVLNSXQPeKIqnDZlxd+3D2WtJw0fTi74jhGnTdLw98WIs540jDj54Isjp03S8PfFSAWMB3JZB98+Um2Ey\/I3boKO02+vWeC7jcRJdA94VMD4ghaoXlZ+BnpsYjPKT26Ey9e+cjR9zDjnjTT809FwK\/Gud95551jKZLc6fa6Uz+DbXhhzM6+dpdonmJeT6L6Dy2mTD\/6+\/ShjPWn4cHLBd3w4bZKGvy9GnPWkYcTJB18cOW2Shr8vRurIG+7MAAAgAElEQVSB8UDOZ\/DtI9Uzz30XtH\/7Y0HyI3ES3QMe9cD4ghaong8\/Az06thnlJye65WxbOZo+bpzzRhr+6Wi4lVAPjANOvoNvH6nGLST7Zl58fNZgXk6iO8BR+EvLF3\/fvoSoR5OHJiUVpafn6Mk6+4QcPVm3atUqoBczoq3S8FF+hmBQtdpQjqaPJ+e8kYZ\/OhpuJVTAOODkO\/hJ+ZG+ecFKuGP3sNqTswbzchLdAQ4VMCkg0VxbeB\/RwoULYerUqVGCUZOfC5ug+bzoEX68y2j48OG1e45Mua1bt9bdIm3M8OWn71in1VN+piHU975XjqaPOee8kYZ\/OhpuJVTAOODkO\/hx+ZE+8o0LYM9xg3IF83IS3QEOFTApICVdmoheGfuixfb2dkCPC710kaa7WL16NWzcuDFKfYFemrgEpb789B3rtHrKzzSE+t73ytH0MeecN9LwT0fDrYQKGAec8gw+RzAvJ9Ed4FAB4yBgXnjhBVi\/fn10eaLZQqLCBJswoqWtrS36N4oZ3B6iQmft2rW19Bg07xfdRsrDT9\/xblRP+cmBarnbVI6mjx\/nvJGGfzoabiVUwDjglGfwOfIjcRLdAQ4VMCkgoZcFhYvxmpikoljNJBJVAePLNL96nHMmz\/rg15vy1ZKGEScffEeH0yZp+PtiZNdTAeOAZN7BTwrm\/f6svdD93N7M+ZE4ie4AhwoYBw+MSSpKk40OGzaMbQsJTUrLlu47tlnrcWQmz2qDXZ7DJs2W7j4qeddQ9ye5ldQ11A0n6aVUwDiMUN7JZ+dH+seLfgb97vhjOHLD9+B\/rn9P5vxIOvkcBq2JRWgQL5phPDCjRo3SIN4mjQvnnMm7PjQJkkIfKw0jTj74AstpkzT8fTFSD4wHcnkHPy4\/0pm3nwMTB6\/qZY3LkWpOonvAE1XhtCkv\/r59ylOPHn2mx6Xp8Wp6JJoeo544cWK0\/YQnmPBj2hoyZAisXLkSRo4cWWeaNHw4ueA7Jpw2ScPfFyPOetIw4uSDL46cNknD3xcjFTAeyIUY\/JD5kTiJ7gGPChhf0ALVC8HPQKawc8HXTs45Iw1\/X4w460nDiJMPvjhy2iQNf1+MVMB4IBdi8JPyI9nmuORH4iS6BzzsL60Q+Pv2qwz1pOGj\/CwDa4q1UTmajjfnvJGGfzoabiU0BsYBp1CD7+qFScuPxEl0Bzhii3DaFAp\/375JrycNH04u+I4Fp03S8PfFiLOeNIw4+eCLI6dN0vD3xUg9MB7IhRp8+0g1elueevXtmkUma3WaF4aT6B7wqAfGF7RA9ULxM5A5rPFQvjZyzhlp+PtixFlPGkacfPDFkdMmafj7YqQCxgO5kINvH6lOMqdRMC8n0T3gUQHjC1qGeshB\/ZQXgW3btpXX+ACWh1xDA5ijIjsEiALa0C0kh0EIOflC5EdSAeMwaBUrEpKDFYNGfHd07DThqAtJOdf1qnKwlAKGHjlFYtDsvfT4alIWYHoc1VzPvmXLlohj9MirIV3IwY\/Lj\/TZvxsB3xkwAea3dNZ43mgbiZPoLhMtrgynTSHx9+1fs+spBs0eAf\/n98WxU4+hP1+KqFkVj2ApBQyKjrikdniB2KxZs2DZsmVR5t+4C8TwDg0UOQ8++CCsWLEC7rvvvrrMv9OmTYPbbrstyknDIWCwzbz5kTjFgu\/k4bSpL74A7HHgxGDn\/kOwa\/8huPBDp\/oOv3M9OvdaWlrq6pkfJps2baq77wbn9Zo1a+Dzn\/+883NoQfxxcuKJJ9bdn4PpHsaNG1c3z7EOri2bN2+Giy++2OtZcZU4xy6YkYEb6ot9DgwhW3NVGptSChh60yldBKmwQQETJ3KQFUn1zQJqkupxCRg7mHfmue+C9m9\/DO59z1S495SpqV4YTrHgO2s4barShPPFlwMDFC6f736pLpD8zJYT4O7Lz2IRMzi\/7rzzThgwYACcd955seIBf1DMmDGjdomfma95BMzXvvY1mDBhQiRgqPeWXiRoxgUvGtyxYwdMnjzZd6h61eMYu2DGMTXE1ecixTYTNE1vlmtsmtGxUgoYuk2EoNFtH7MlhH9HD4v9Kw\/\/jvVNrhoKepKw4RhwO5i3Z+8UGHzkdXDJj8QpFnxJyGkTB\/6+\/WxWvdAY4Ivg7K88AyhY2s8fHP0X\/7Z0zWtRF9OO8vvggOJgw4YNMGnSpF4eFfsmYvSkoscFvSdTpkyBQYMG1XlgcK6itxSTZpotYeN1pdvBmH8K65tt5l27dkXemEceeaSXB4ZuJ+ONx5gJ\/MUXX4QzzzwzElwogNAzSz1C9GZlumVN8Qk9dj7YF10ndJ+LFttF41Xk80KPTZG2288qpYChnTCLTmfn0fgRdA2jcMGP2UKi20G44JiMwOaqdiybtC2F35n93JDJ8q5+5HXY9JNDkZ1DTukHy0fvgMGrroC9U74JHVuGwZ6fH4m+G\/P+E+DeSYPqxo0jMV1eEnLYpMnyjo1K6EUHtzF3vnUIVl97TiRezMcIm7Sj\/D58Mds2KCYWLlwIU6dOrdvWwTloPDDoqTFbPHE\/OPBvKE7o3H711VfhtNNOi360oMi4\/\/774eabb47WA+OBMXYnbSEZD0xbW1udjfR5RsCgMFqwYAHMmzcveiaWwY\/tvQk9dj7YF10nZJ+LFttcW5nY7hNPPAGXXnppbTgoT+n7CAvEbX2GGMeQYxPCnjxtlF7A0Gy\/CAT1rOAihR8jbrJ6XgywHAMelx\/pj5YN7hXMizbYR6o5vR2+ZOK0iQN\/3342q56NAQqQPB+8fwiFy5mnHRMvNRHz1qHIG4Mixvez+rpz6qpScYILddzL3pS58sor4atf\/WpN4MTFwFBvCQ3Kxzm\/fPny6Nkmp5SvgEERddVVV9XEiRFMxh4URcYLZDob54Xpi\/wN2eeixbbNVTO2eWOx7PrGe2fnPjPPo1ufvvMwrl7IsQlpl09bpRQwdPGjgbv79u2reWBwkZw9eza0trZGv4iSfh01UsCcAgbbdg3mvbv9LMBL7syHUyz4kAjrcNpUpQnni6+NAW7\/+H5QnOCHel5oW2nfuzwXRTf92PMv7iWR5IExW0\/mh4j9fPP98OHDa16ZEB4YW8Dgc3EtMc9DcUM9MEm49EX+2n3GH2y+n4n3PB+tf1POH9yrCRTiuO2Ja2ScGHd9Jg1gv+GGG6CnpycSwGlbmfaJWBPOQMMc0OOIYgQPl2C7WAY9fOiNwW1J4ymkHhgjbrDurbfeCnfddRdgcDvaZPPcxIclJZC1MagSH0spYJJIgwMVd4ya7pebwbRJRQfZDu7jGnD7SPVnhx+Gzif\/OPVINadYcJ3wdjlOm7jw9+1rM+qFxgBjsDonjIDOCcN7dWfpmu3RSyFUHIw9X5PmGhUwWAe3gNGNjnMVtxPpKSQ6z40HBn\/A4LYOfnD76f3vfz\/cdNNN0bYUvjhoJu9GW0hmawjja4wHhq4h1B4aA4PP5b6GoRnc83mmzVfXCzx9nhWizv7bP1ZrJstWZpxXBrm7devWiLf4MVzDGC07GL3RD2jjgcF69lZlnDeQCqEkfqM9odeSEPj7tlFKAePbWd96nANue2GWv3ETdJx+e\/RyMQGVaDd9mXCKBV+MOG3ixN+3v0XXC42B4Z0tUtD70nbP0e0p24tSdJ+r8rzQY1cGXOw+48lL38913S9FHpiPfOi0Xk10P7c3OkWH62WSR9HludTDnWUr0wgUs21pBKz9oxl\/FOcRMPR0Xlw8FgaZG\/Fu+hsnplXAuLChYmU4FyCf\/EicYsF36Dht4sTft79F1wuNgfH+xZ1C4jxKXTRuEp4Xeux8+hR3eIF6sajXmXqV6EuQetLo5aFoj90WvkxDXZZWtNj23crEeugpuf766+Fb3\/pW7TqAEB4YW8Ag5nQ7E0\/2xW1FxXFFAh99OBxXRz0wDkhyD7ire9UE83KKBQc4Yotw2sSNv2+fi6zHgYHxtpiYF+wPihf7ZFKR\/aziszjGLgtORlzQYFF6ZQRud5iTmShSzP1Z+IxFixZBV1dXLZDZHJLAlzLGYuBLNK4t3LILJWCKFttZtjJtT4uJb8H4S8QAPxjYjVide+65UdA3joOJ52q0hYQYYxt33HEHPPbYYzVBlLSdaV8vEnfPkXpgssycipTlXoDs\/EgGNpOd2vy\/CeblFAu+Q8ZpEzf+vn0ush4nBuZyMOxPEbfxFombhGdxjl1a\/8yxcCy3ceNGWLJkSXRJIL7szP9T0RJ3EAIv9kSPCz0UQT06q1ev7tUW3vcTSsCg7Sq200ba\/ftm8tHdSreS6oFxwIl7wO1gXjy+ivu69sfczcEpFhzgUA+ML0g56nFzMIdpWjUFAQljRwWLETDGm2LfpWW8MdgtI1rw1Az+29xSjgLG3LmFF\/7ZbWHwdUgBYyBWsZ1\/ukngY\/5eHG1BBYwDkkUMeFww75jDm+HhK1\/rFcw75Li3YMSIEQ6WF1eEU1QVgX9xSPk9STHww01CLQljVxUBI2E8y26DBD6GwlAFjAOSRQx4Un6kBS2d0PO\/M1WbD3ph7rj0VBUwDuNWpSJFcLAIvKQnc8yCweOPPw5nn312bLoS2o6EsYsTMGXaQsoyLlq2MQIS+BhqjFTAOCBZ1IDH5Ufac9wgOHLD98AcF0RzMbXAZz5yloPlxRVRDwwv1hwcPPLT7fDTu6fB4PlPRMbj\/\/9i\/f3Q\/\/fHwQl\/cFHwDpUhmaNrp01fzD0xjepxjJ2rnaacLWDKFMSbta95y0tJJZClH65iGtuUwMcsfWtUVgWMA5JFDTjdRsLTIPe\/79tw8mNLYedNz8NfrTsYBbLhB\/MjrZ1Zf9OpQzdYi6iAYYWXZdHZO\/9jcPDF9dD\/Dy6C9123MhIvb\/3jfDj5or+I\/j\/0R3oyR+wvTVFgjgqbW73NqRI8aYKfOXPm1NIV2HlsKHZFrR+NxssWMFjWnFqhqRjw7\/QYNT3JYifcNAHBcW1hmoWQMTBFim0pqQRc518WMa0CxhXVCpUragGy8yPhzbwXfuMCp\/xIzYZbBQzvCNgcxAXd9dPvfb1v2zX10QODIgbL4N+MmKFtN6qfZENcnTIkc6Q3mJojrnhUGE\/V0JuAs7w0ilo\/XPlQRLnQfS5SbEtIJWALSZNjK+4W+ixiWgVMEewX9ozQk7FR91yDee38SM2GTAUM7wjYHNz22Xc5P\/ADD7\/Tqyx3ffuZZUjmGJfywHhh8CIxetsqnsqhuZIaDUaR64czKZgL2n1Goez6ifP+Yf1frP9G1IQR2\/hv9BbWie3Th8Npn53X61Fpz6fPlJBKwFyKZ2c6x7tk7HQEWcS0ChhXFlaoXJELkH2k+oahO+CKZ6ZBXDCvnfG3mZCrgOFF3+bgWw8vcH5g3IKO9Y+8sT3yvtjenNP+5\/y6tpPqNzLArlOWZI6Ncshgf82LBZP8PfDAA7VcSSpg6hHwFdxx25dUvKSRPk6sY500wU7rSUglEJe\/z3hhaMZ1c3Geq5hWAZPGoAp+X6SAQfhsL0zP3imAwbzbrng0MT9Ss2Evs4CxE\/IZLJNussyCNS42+DE3b+a5pj1kTAHahCKGxrwYN33oGJgyJXOkMTCIEd6aeuWVV8I111wDmNwRP\/S21Z07d8KKFSsankQqev3Iwk+usqH7bALMka\/mg9udJgA9ZD8kpBKwPTBx\/fMR0ypgQjKlJG2Fnoxp3baPVF97\/Ab43I8XQMfpX4VNx4+GiWedBD0v\/RLMxXZp7RXxfRkFjBET5pcNxYm+dH2FjN0+Lji+17SHvtkU+3roxfVR4C51n+Ov3ZMvmspyCqkIHkp8RtHrhwQMQve5KLGN2ElJJRCX6dykI\/AV0ypgJMyOgm0IPRldzLePVP9g18Wig3nLJmCSThrEjQ0mScNf4i0tLS5DF5Uxv45GjhwJBw4ciDww9PZSc7LF9Zr2kLllnDuhBYMg0Iz1I4jhORoJ3WcV2zkGw6oaemzCWZa9JT1G7YBZMwbczo+0\/I2boOP026O08UvXvFazWkowb9kEjMOw5ypiYil27NhRu2ad5o\/BxrNc044BpKG3kHJ1UCs7I9CM9cPZOKaCfbHPTFAGb7ZKY6MCxoEezRjwrPmRHLrBWqSsAsbEPJx55pm1RHfUg2JOAWQBz9x3gl4XDF41eWJUwGRBsTplQ60fnLFaodEO1efQdml7epFdUA6UYVI2azLawbwG+I4PnwrLnzuW7LHn2nOankW4zAJmwYKjJ3poQKZLEF3SRKCnBGpj1tEB48aNqyXA0y2koMuI6Mbyrh\/csVoc4OXtM4dN2uZRBKo0Nk3zwJRpUjZrwO1g3vbzBwH+zXzM\/0sI5i27gEFPy9atW+ELX\/gCrFy5EgYOHAgobHw8MHShpB4YaUG8zVjQm50Lif5gwrgmepus8bxt3rwZLr744mDw5Fk\/uGO1gnXSaihPn22bfvOr3Q3NfPeJQ7m6Ucl2Q45NswFqioDJOyntY5nmsikMsjTXYyOw9HQJXbjsq7NpnbgTJ80ccDuYN4kwm2+5ADD9QLM+VRAwyB9z\/8LcuXOjC6NCChgcmzzXtIeMgWnGS6HZuZBsUWrfTWPGB+OW8PbdUJ9mrh+h+pC1nZB9PvLms\/Dzje2JJpz68SchpIgpIheS\/Q6Le+9s2bIFTjzxRMCDACE\/IccmpF0+bTVFwPgYSuvQX7J0cPHlgxdMLVu2LPoFPX369Oj0x6hRo2rHV7E8\/RW4b98+WLRoEXR1dUW\/wLu7u3v9KmvmgNvBvO\/tdxj+5bU\/hu5PPQFdPzx2w2qzg3mrImDMr3DkDn7S7vjIy2XX+qE5WPRLwYgDPA4+adKk2NtEMagZT1uhZwTnMV7dj4v4lClTYNCgQXVX+dOLvswPEjPnsQ5+8L6WYcOGRfXpjxyDOYpJKlboHTDohVu7di28+OKLgPFR5513XvQiGTt2bCRyzW2oVJDGHcfHZ4UYO45YLVfu+ZQL0WfzXMPVlrZjhxfMd\/tXj4DQAqaIXEhUPCdtV3\/ta18DzCmlAiaZgU0VML6TkmZSpUdbqbDBxczcuWETgNbHRSourTytE3IyZl0M7GDeJYO\/Dx9\/bo64m3nLJmDsC8vixiXupZd1\/EKVD83Bol8KiIOkXEg0oJomYjSiBlMFLFy4EKZOnRq9QPCFg2KIChgURnSbMc6jE1LAhI7VCsXNuHZC8rVorhaRC4l6+OJSARhhjGvQrbfeCnfddRds2rQpEvfDhw\/vxUUU+3QnAcV7khcx5Nhwcsil7aYLGJ9JSQfK\/NIyg2VeTI1+PdO4BDtGwXhtcKEyn2YPeNzNvIOPvA4PX\/mamJt5yyZg7MmR9EJzmURFlLE5ePDlv8n12N\/86idweNc\/Qf\/fvaFXOwdfviP271ke2P93\/6quuKRcSHTOx\/EAvTJ2rqM4AYO\/jqdNm1a7oRfbivPChFg\/6K90jlitLGPrUtbuM3pKJH+od6fIXEhmKwnvg6LvHMTKeGDitj5tMY1cxPuqbr75ZkBB3iglRgg+ShlLEQImTwClESzmqnYcOHT74ydOjNgvKlcBg+2tW7euKeOGt+7Of+zN2rOnnrwFrv\/\/bqzdzGu+GPP+E+DeSYOaYuPu3bth6NCwwXTjx4+v9SVk\/EccQGUTMGV6ISDeUnIh4fqAv2CTfp1SDwzNL0PtN8fkr7rqKqdA7xAvDHubIXSsVuhFw+7zgee\/6P0IjNf6r33PwvFnfKZXGyjC4\/6e9WEDzrm1VqWoXEhJoRDGECpgMJnojBkzInESJ6ZxixM9gvST5IUJwces+HKVFyNgfAMojYJtbW2NMDJ3buC\/7Tw0cb+88G+St5DMwMfdzLvljE\/BofavQ\/dze+GpV48eq27WkWr1wHBN0aPt2otOWhBumjX4QsCXCsYP2J+3H\/sonNLanSswkgZVSsmF9OUvfxm+9KUv1XlM7EXeuO7RM4zXtaNIoWsT\/g3d+iiu0W0fd927LY5CvDDi4iRcvM1pPOD6PkSfjW1FbyEVkQspKQyCjkcjAYPlkGdGTGNcGfXANBrXkGPDxR\/XdkUJGDTaZVLSX0M0cBcDco0HxtyzgcIGBzppf5oSSWIQrxlIuo005JR+cP9x98BpP3oQ4CsvwaSH34Kd+w9FRZt1pFoFjOuU8ysXetEp+qXg1+tq1Ao9dmVAJWSfiw44LyIXEo4h3q5tPnHxdvguw6D2O+64Ax577LGaB4YGsFMxbYdWJOVwCzk2zeZiUwWMb+ftX3T0V1TcMeq41ORxR6\/t49XGPgkDToN5cauofeAOuOif2+He90yFe0+ZWgdlM45Uq4DxZbNbvdAcLPql4NbLapbKM3ZlCzbnWjORr0mffu89Fq9YTQaF7VUePoa1JH9rTREwZZuUUgbcNZi3GUeqyyZglIMA+lLIv4C6tBBy\/ZAeq8UlYCjOOzb+Eoa1nuQCvZaJQSAkH5sNcFMEjN1p6ZNSyoDbR6pvGLoDrnhmWq9g3mZsI5VNwDR74mV9PjcH9aWQdUTcy4ccO+lrJbeA+dmu\/4S7z\/t3+LNHP6Aixp2CdSVD8tHThGDVVMA4QClpwG0vzA92XQwYzLt1wt809Uh12QQMemBwfxnv+Uj7YHAc3r9A7xxKqxP6e04O6ksh9GjVtxdy7PqygEGePvDpbYD\/fc8Zvw3X\/eD3eAeuoq2H5GOzIVIB4zACkgbczo+E6QNMAC92pVn5kcomYBArzccF0ctAXwoOi0COIiHXj74sYB749Kvws13\/FQkX\/Dd+\/vzRD+YYmfiqUlIJZOnY448\/DmeffbbTj6yQfMxiI0dZFTAOqEobcIn5kcooYMzQ9+WM6EW9FAzW0pM5OiwHtSKuL42Q60dfFTBP3vof8KOH3oLL7hwabR0Zr+FHv\/A78NEv\/k6WYUstKyWVQKqh\/10g7ibfRnVD8tHVRq5yKmAckJU24Jd0PQObfnKo5m0xXeicMKJuG6nIYN4yCxgHCjS9CAcHi3wpIIBlSOboOtBZXhp5xq5sweYGvzx9tscgSaxg3NY\/fHpb8HgYCakEqIcY\/21ueI47gYvfz5kzJ9rmtrOrx\/E55Ni4zheuck0RMGWblNIG\/J+eegmufuT1GicwaNdcZEeJUmQwrwoYril6tN3QHCz6pWAWZMnJHNFGujaZqxbMnVIYM4UfvLYhy0sj9NjxMi1M63af0Xvi+vnDyafVipotTvxDnKcFRbj9Ha1vGkp7Pq0jIZWAfXGhucfs3HPPjU2ESm+NTsO5SnxsioBJA1ja99IGHMXCDd99u5dowXiY9vMHNyWYVwUML2tDctC8FIa1DoDL7jyjl+F0Wylkr8qQzJHmkMH7ozCAGy\/CROGFN++aT5EemDIFmxt8bL4u+p0fOVHpzNYBdXEt35mxK9o6cv1gfAwG+NqftOfP\/Y8\/rFWRkEog7u4y44VBjppL8FBM23m70rAKuZakPYv7+6YIGD0Bkm9YUSw8+2Z\/uK77pVpDXzvu72Hs9gegbfAq2NPvWD6korwwKmDyjWlabXvRwaOkrh\/7tIZ5KcQt9KZNc9ID\/z\/utEfa8+06ZUjmGJfywHhhMBeN70sj7wujTMHmSQIGt3tcPnH3uxiuYfyL\/UEPDAb2mu+S7odJez6tJyGVQFzqCLvvpsysWbPggQceqKW9SMM5Lx\/T2i\/y+6YIGOxgmSaltAE3YoEG8445vAWWv3FjdKS6GfmRVMDwTlubg+YUhstT405qmF+k+IvX\/uzceCD6FfueM34r+iquftrz7TplSebYKIsvYuHz0gi1fpQh2DxJwLjwNKlMUqwLemZQwJjA3jzPoHWlpBKIy7OFW0g0A7rxwMyePRt27twZJTJOu+4hFB9D4Z2nnaYJGGN0GSaltAE3YoHeCYPbR987MhcOvrge+i0\/CG33PF9ofqSyChiXXzp5JliouqE5WORLoUzJHO34PAyMvPLKK+Gaa66pJYHM+tIIPXahOMXZTug+o1B58rb\/qAXscp5C4sRFQtuhx6aZfWq6gGlm512fLW3AjVigN\/PiVtG1x2+AUWtuhAUtndAzYEJd97jzI5VVwCBI+Kt7+PDhUayD1A8HB\/WlUMxoc4xdMZb7P4WjzyY2688f\/QD0zNiV6B30t7pv1OQYm2YhpwLGAXlpA07FQlx+pD3HDYJtVzxa6JHqsgqYRifi4jLEOtCFpQgXB\/WlwDJcdY1yjR2\/5f5P4Oozjb1CIdMojsvf+mrX5BqbZqCmAsYBdWkDTsWCnR\/pK8c9BJ\/cvrzw\/EhlFTAOwy+iCCcH9aXAO8ScY8druX\/rXH3Grc\/vzNgdPO7Fv6flq8k1Ns1AQgWMA+rSBtwWCxLyI5VZwNAYDYx5wKj+W265BebOnQsjR450YAh\/EU4O6kuBd\/w4x47Xcv\/Wsc\/6kYvAtm3b5BqXwTIVMA5gSVuAbLFg50fCYN737X0Wzjvj8ah3ReRHKquAMeJlyJAhMGnSpOjej5tvvhlWr14NGzdudLrZ0oFCuYtwc1CzUeceosQGQoxdWYLNk0AIgUHIEeJcr3zt5LRJGv6+GNn1mipgyjIppQ1+HNHpkeqJB9bAvP1LCw3mLevkoxzct29fTcCgsFmwYAHMmzcv9VhiqMnYqB39RVsEynzPCPGLtwzB5ipg\/DlU1jXUv8f5azZVwKD5ZZiUZRAwZhvJeFt69k6Bfu8bDuv\/pLuQYN4yTz7k4J49e+Bzn\/scPPzww9GRWbx1dezYsdDZ2Zl\/ljG3UAZ+MkOQ2nyZ+YmdK0uwuQqYVComFig7R\/177l+zqQKmLJOyDC8IO5gX74XZuf9QL2Zw3cxb9skXd2mU77Fqm9d4b4hpiz6H\/p3G4dinn\/ASOEzWhp9Vq1ZFwop+ysBP\/yUqTM2y8zMMCs1rRTmajr1yNB0jUVtI2c09WsO+GIsu+HSxN7kj6HPi0tHbuSXsF1dZJp8dzGv6bbwy5v97rj0HLvzQqb7wx9bTyXcMFupVNDlNbrvtNhg1ahTMnGDcQcUAACAASURBVDkzCg7Gz6JFi6CrqyvaokLebt++PfL42PVNua1bt0J3d3evuJyy8DMo4TI2pvzMCFjg4srRdECVo+kYVULA4C9c8yKgp0TwZYEnSJYtWwYDBw6E6dOnRy8E84vViBuadhwFDb4w8ApmjIEw9Wm7ZZl8djBvEh04vDA6+eLRNmK7vb09KmC4ZjIc499RgONV4K2trZGnhopsGkyMbcXxviz8zL48hatRBX6W4bRc0ogpR9O5XAWOpvcybImmbiHZ3hTXI6woVOivVwMJFTYoYOhijy+FHTt2REXp6RLaFgqYuHbLNPloMK\/BJS5Ldeibecs2+RptXxrcQlxkR\/lFPSj4DCNaMJss\/hvFDIptKqrXrl1b88wYm6kox3bKxM+wy5d7a2Xjp92zspyWUwHjzkm7ZNk56t9z\/5pNFTC+k5JuE2HXaSyBWeTx73GJrbCufTzWuPnHjBkTe2y2TC8Is42EouXuy8+C6x48mrHa\/NvExYTeRtLJ13sS2p5C6llRAeO\/aPnULDs\/y3JaTgWMDzuP1ik7R\/177l+zqQImxKSkv0oRBuOix3\/bW0j4N1vA0P\/fvXt34hYS1l23bp0\/0gFrop1Dh\/ZOLY+P6HnplzD\/sTdhyCn9YOLvnQTvbPoWXL19Ya+bee+dNAjGvP+EYFY1ssn3IePHj69VDXEM1dcOn3pxXkLqWQm9hVQWfvpgGaJOFfhZ5tNyZfoRGIJvPm2ogMmOWlMFjBEceY6wGi8Oxg\/gxwRCmrbxv\/QorC1gaMAkbYsG8pZp8tmnkbD\/P9h1MXxnwAR47PzF8NSrb0c43d1+VnTBXahPmSdfXLZkGieVFSMUL+ZCPBQq5kM9Mvg3DeLNiqx\/+TLzk\/Y65Gk5fzSz1yzTGpq9d2FqVIWjYdBwa6XpAgbNzDopUYTgB0UGDdzFGBY7SNIERho4Gnlg8EUW57Up2+SjcTAYsPvxf50Dlx1YA22DV8GefkdFS+hA3rJOPrqNaYSu+RvitGTJEqAiJG1axYkhrGOOP1Ou0yPRdoAmfa7ZMsXbgleuXNkrvUHZ+JmGIcf3ZeUnBxbNaFM5mo66cjQdI7uECAGT1Wz7JUFjYNKOUcfFwFTlGLXB0T5OnXQzb8g4mLJOvqTboMtySzSOub4c0leQsvIzvWflKKEcTR8n5Wg6RpUQMNm7ma9G2SaffZwat4o6n7wUjvx0e5QfiSM3UpknX9LdQMOHD69dQJePQby1y8ZPXjTiWy8rP4s6Lcc9JsrRdITLytH0nvGVaIoHpmyTsoyTzz5ObbwwHad\/FTYdP7rGqFDHqcs6+crGxbiloIz85FvSqiVgisaJ63nK0XRky7qGpveMr0RTBAxfd3haLuPks3Mj4bHqR7a0wpYzPgXT4YYaUKGCeXXy8XDPpdUy8tOlXyHLKD9Dopm9LeVoOmbK0XSMdAspO0aljDGwTyN1ThgBo9b8FYze9W2WYN4yT77Qp5A8KJarir4c0uErMz9N78rMU+Vo3+Boei\/Dlmi6B6YMk7Ksk88O5h1zeAssf+NG+JfhHXDLryfXmBQimLesL4jQp5DCTk+31srKT7fehSlVVn7a4gVPooU4LRcGVfdWlKPpWJWdo+k9DF+iqQKmLC+Psk6+IoN5yzr59BRS+EWFkwu+1nLaVMT6UHaeFoFRFm5w8iGLHbQsp03S8PfFSNQWUlkmpbTBz0L0uNxIceTJG8ybxaas5OXGX08hZR2RxuU5ueBrKadN3Pw0fS4zT4vCyJUfnHxwtcEux2mTNPx9MRIlYNCYMkxKaYOfhehxwbyYDwmDek1eJByHvMG8WWzKSl5O\/PUUUtbRSC\/PyYX0p8eX4LSJk5+mN2XnaREYZeEGJx+y2KEeGF+0jtZr6hZSWSZlmSdfXDBv97\/urRMvSIS8N\/NyLgjS8M835cLXloYPJxd80eO0SRr+vhhx1pOGEScffHHktEka\/r4YifPAhOoIZzvSBj8r0e1gXoOV7YXJE8yb1aYs4yUN\/yy2F1FWGj6cXPDFk9Mmafj7YsRZTxpGnHzwxZHTJmn4+2KkAsYDOWmDn5XoccG8+DdM8og385pPHi9MVpuyDAM3\/mU4CdcIL258sowVluXkQlZbTHlOm4rCv8w8LQojV35w8sHVBrscp03S8PfFSJyAKcOklDb4PkRPupl364Svwj2Hx9WyVPt6YXxsciUxJ\/5xJ+HQLsyPhVnSsyZzdO1TyHKc+PjYyckFH3u4RVUR+JflxGbS+BSBURZuKEezoCW3bFNjYMoyKasw+eg2Em4d3X35WTB4zmnw1h9eDtPfNaMWE+PrheFcEDjxL8tJOPXA5FtEy8pP0+uy85RzDvswg5MPPvZUQWT79jtPvaYKmLJMyipMPhrMiyKl\/cOD4YTuv+x1My+SyedINeeCwI0\/elt6enpg5cqVMHLkSHjllVdg2rRpMHHixNqlYXkmGXddbnyy2s\/Jhay2VGkLqQwnNtUD48tQ3q1XaWuEP0r1NZsqYNCUMkxKaYPv+4Kwg3mHHHkdVu+dAs+dORmufaejxgyfI9W+NrkQuQj8H3roIZgzZ07NnMWLF5ciEzUaXAQ+LuNUhFjIYgctW3Z+luXEpgoYX4aqgPFBrqkCpiyTsiovCPtIdfv5g+Cif54CYw5vzh3MW\/YXhM\/kkVKnKvzkxFP5yYluetvK0XSMlKPpGNklmipgspvbnBpVmny2F2bigTUwb\/9SeOzDi2H23j+qAZw1mFcnX3O4qR4YN9yrwM8yHHhQD4wbH+NKVYGj\/r33q9l0AVOGSVklAZOUH+n7Pz8VOk6\/HdArg2WyBvOWdfKZeBc8cWQ+o0ePhhUrVkBLS4vfrCq4VpX4yQVdWflp8Cj7aTnlaDqzy87R9B6GL9FUAeM7KW3RQ184NJaho6OjVxBmXMwN\/m3KlCkRunF1qjb57CPV8\/cvhcsOrIGO078Km44fXWNZlmDeMk4+M+52vIvh0KpVq2Ds2LHhZ13gFqvGz8DwRM2VkZ8Uh7IceFAPjD97y85R\/57712yqgPGdlFhv5syZMHfu3OjUiPngr+lZs2bBsmXLYODAgTB9+vRIwJiXkHkx4ekSc78H1lm0aBF0dXVB\/\/79YeHChTB16tS6dqv2gojLj\/TIllbYcsanYDrcUMMzSzBv2SafEcHt7e2xIiVO6PpPM96aVeMnB1pl42ccBmU+LaccTWd1FTia3suwJZoqYLArPpOSig7q5qfCBgUMFTn4QtqxY0eE3saNG2sCBkUNfiZPnpyIbNUmX1x+pA9889PRFhL9ZNlGKtvkSxLPpv9p34edhvlaqxo\/86ERX7ts\/EzCoKyn5ZSj6ayuCkfTexquRNMFDHYl66RsVN6cbMJ24+IYsK4tYF544QVYv359dPNqX9hCQmyS8iOhaHnq1bdrDHMN5i3b5EsTKGnfh5uC+VvSl0M6hmXjZ3qPylVCOZo+XsrRdIzsEiIETHazj9UwggW3ioxHB4ULfuwtJCOWqIChV8bjtkJcnSpOvqT8SPZYuHphyjb50gRK2vd5OBu6bhX5GRqjsvEzdP+b3Z5yNH0ElKPpGFVOwJhYhtbW1qhv27dvrwXuojjBjxE3cQIGPTKmDm2Lbinh5MPPunXrsiPMUGP37t0wdOjQ3C2PuWt7bBsTzzoJ9vz8CGz6yaHo+3snDYIx7z+h4fNC2UQfMn78+Nr\/btu2LXd\/aQNluYPIpdP6ckhHqewvh7KfllOOVp+j6T0MX6KpHhjfSUnjVmjg7r59+6KYGvTAYEDu7NmzAYUNFSP2FhKNp0F4+4oHBvsalx9p4j3Pw\/8YeRps33cwU36ksr8gwk+t4lrUl0M61mXmZxVOyylHq83R9N7xlGiagMkzKe1j1PQYbNoxalvAGK+MuUY+7gr5qk6+uPxIS9e8FgkXTDOwp9+gGuvSjlSX+QXBM7XCtUo5HXe0u6r8DIdgeY9RV+W0nHI0nc26hqZjZJdoioAp26Ss8uRLupl3\/Z90wxd+9Ds1vqQdqdbJl33yudSgHsKtW7dCd3d37QSdqV9lfrpg5FKmrPxMi8VK+94FmyLKKEfTUS4rR9N7xleiKQImbdKlfc8HR3zLVZ58cfmRpqy+OAJi4uBVNUDSgnl18vGwknoMUfjH3X9UZX6GQrWs\/ExbC9O+D4Vf3naUo+kIlpWj6T3jK6ECxgHbqk++JC\/MNy9YCXfsHlZDqNGR6ipNPrNlIyGlAA0ypyfu6A3BVeenwxRNLVJWfqYJlLTvU4EpqIByNB3osnI0vWd8JVTAOGBb9ckXd6R65iO\/D98ZMAHmt3Q65Ucq++SjAeUoXG677Ta48847Yd68eU3NiaQCxmGCOhQpKz+rclqu6muoAwVTi5SVo6kdYyzQNAGDp322bNmS2DUJv36NcX1h8iXlR2obvMopmLeskw9PrS1fvhxsvkn5Zeu6hYRcrdox\/5DrXtmO+Yfsu4S2+sIamhfnsq6hefudp35TBEweg5tRty9MPjs\/0sQDa2De\/qXQ\/aknoOuH79RgTwrmLdvkS9qOMR2VImA0iDfMjC8bP8P0Wk4rfWENzYu2cjQ7gipgHDDrC5MvLj\/SR75xAQw+8jqcd8bjNZSSgnnLOPnotpF9fF6KgEHgTUzOkCFDYOXKlXWJRvH7vsBPh2nasEgZ+Zm3z5LqK0fTR0M5mo6RXUIFjANmfWXy5QnmLfvks+9aGTVqFCxYsKDpMTAO9FQB4wBS2fnp0EXRRfrKGppnEJSj2dFTAeOAWV+ZfGnBvAaqOC9MlSZfUlyMA1WaUqSv8DMPuFXiZx4cmlVXOZqOvHI0HSP1wGTHqE\/9wo0L5sWTSJ0TRsDTP36rlqnaPlKtk8+DWIGq6MshHUjlZzpGnCWUo+noKkfTMVIBkx2jPiVg7PxIKFyu634pOkr99KtvJ+ZH0snnQaxAVfTlkA6k8rMeIzsdCz2FZ9K8YA0aG0br2Kf2NN1FOgfTSihH0xDq\/b1uITlg1pdeEHH5kVDAxH1ofiSdfA5EYirSl\/jpC6Hysx45DFKPu9WZ\/h1rLFq0CLq6uqK7kOidRLjNOnz48ChRrp6U82VlfT3laHYcVcA4YNbXXhB2MK+BCL0wGCdjPvRItU4+ByIxFelr\/PSBUflZjxoVHShOzAe9LyhOVqxYAf3794fZs2dDe3t7dE8S\/ru1tTUSLVjO5OVavXo1bNy4McrRpekufNh5tI5yNDt2KmAcMOtrL4i4\/EhUuBjIaDCvTj4HIjEV6Wv89IFR+VmPGt3ywW\/MVhEVJvh3I1ra2tpqYgbTWFChs3btWti+fTt0dnaCprvwYacKGF\/UVMA4INcXXxC2F2b+\/qVw2YE1vS62M8G8+oJwIBJTkb7Iz6xQKj+TEaOiA0sZz4oKmKwsy1deOZodPxUwDpj1xReEfaT6hqE74IpnpiXmR9LJ50AkpiJ9kZ9ZoVR+JiNmgnNxe2jYsGFsW0hogaa7SB4HTXeRdVYDqIBxwKyvviDsI9XL37gJxhzeXHczL8KHwby\/\/tleGDFihAOa2YtIwz97D3hrSMOHUyz4IslpkzT8XTDCLST8mCDcWbNmwbJly2DgwIG14F78XoN4XdAMU0Y5mh1HFTAOmElboDiJTuHIkh9p7HsPqoBx4BJHkb7KzyxYcs4Zafi74GIfo6bHpekx6lWrVgHGvOCH1pk4cWIUtIuBvvjRdBcuqDcuoxzNjqEKGAfMpC1QnESncMTlR\/rs342ATcefDR2n314risG8d1x6qgoYBy5xFOmr\/MyCJeeckYZ\/FlyKKisNI04++GLKaZM0\/H0xsuuVUsA0uoSJRtd3dHREkfH0Q6Psza8H8z0eH8SPXUfa4HMS3SaIazDvvZMGwWc+clYoXta1Iw1\/lk7maFQaPkXy0xU2Tpuk4e+KSZHlpGHEyQdfXDltkoa\/L0aVEDBJlzDh3QZ0L3f69OmRGDEuUCNubPcngmLcpnGiR9rgcxLdJogdzDvz3HdB+7c\/1iuYd8z7T4C1My8IxUsVMBmQ7Mv8dIWJc85Iw98VkyLLScOIkw++uHLaJA1\/X4wqIWCSLmGiwoYGo40cOTISKDt27Ij6by5dMh4YrIeZh7HcgQMH1ANjsSRLMO+ZLSeE4matnapOvlBAScOHcyH2xYzTJmn4+2LEWU8aRpx88MWR0yZp+PtiVAkBk3QJE3bO3GmA\/8bbJOktk\/g3rGsLGNw6GjduXCRwzIVMFChpg89J9Dhi2cG8nx7wMtz879fANy9YCXfsHlarQm\/mDUVQbEca\/iH7FqItafgUzU8XDDltkoa\/Cx5Fl5GGEScffLHltEka\/r4YVULA0E7YlzCZa7CxjL2FFCdg0DOzYcOGyOtCc32ogDmGQFww70e+cQHsOW5Qr2De1dedE4qb6oFxRFLa4sS5EDtC0qsYp03S8PfFiLOeNIw4+eCLI6dN0vD3xahyAoZewoSdox6UuKBc2wODZZYvX16Hix0Hg4OPn758CdPVj7wOm35yqIYT3sw7v6UTJp51EvS89Mva3zGYF+NhQnzGjx9fa2bbtm0hmqxkG9IWJ86F2HcAOW2Shr8vRpz1pGHEyQdfHDltkoa\/L0aVEDBJlzDt27ev1y2SJvmY6XjcFhL9TreQ4qllB\/PaiR1NLZofKRRJqzr5qooP50LsixmnTcrP9FGRhhEnH9LRiC\/BaZM0\/H0xqoSAaXQJU9oxahUw\/tSxg3lNSxj78vf\/67Wah8bkR\/J\/Un3Nqk6+quLDuRD7YsZpk\/IzfVSkYcTJh3Q0VMD4YlQJAROq867t6OQ7ihS9EwZPG3VOGAHXdb8U\/febz+yCPT8\/EpUL7YWRhr8rb4oqJw0ffTkUNfLleY5yNH2sOOeNNPzT0XArUcqL7Ny6Fq6UtMHnJHoj1GgwL4qU9g8PjgRM3AfzI4U6Ui0N\/3DMCtOSNHyaxc9GaHLaJA3\/MKwK24o0jDj54Iscp03S8PfFSD0wHshJG3xOoqfBY9\/Mi+V\/sOtiuPfi1XDvKyfVqoc8Ui0N\/zSMiv5eGj7N5GcS9pw2ScO\/aP65PE8aRpx8cMEjrgynTdLw98VIBYwHctIGn5PoafDYR6oxmHfmI79fu5nX1A+5jSQN\/zSMiv5eGj7N5KcKmKLZ5\/Y85Wg6TpzzRhr+6Wi4ldAtJAecpA0+J9Ed4KiLhcHyeKT6sgNroPtTT0DXD9+pNREqmFca\/i4YFVlGGj7N5qf+ui2SfW7PUo6m48Q5b6Thn46GWwkVMA44SRt8TqI7wAH2keobhu6AK56ZBgtaOqFnwAQwR6xDeWGk4e+CUZFlpOHTbH6qgCmSfW7PUo6m48Q5b6Thn46GWwkVMA44SRt8TqI7wBEVsY9U9+ydEv194uBVdU2ECOaVhn8aRuZ26C1btkRFFy9eDJMnT47+bZKG2n+nVwOMHj26Lg0GvRpg1apVteSkxg5p+Ejgpz1GnDZJwz+Nn834XhpGnHzwxZfTJmn4+2Jk11MB44CktMHnJLoDHFEROz\/SxN9+Eea9ej1LfiRp+KdhhLc7Dx8+PBItmHh02rRpcNttt8GoUaNg5syZMHfu3KiJRYsWQVdXV5Svi6axsOubclu3boXu7m5YsmQJmESk2I40fCTwUwVMGkuL\/V45mo4357yRhn86Gm4lVMA44CRt8DmJ7gBHVMQO5u348Klw1bfOZQnmlYa\/K0ZYznhW2tvbo2omVxcKkNmzZwP+HT0u+G9zazR6aYxQWb16dS35KLZlBBBmTlcPjPtIcM6ZMvPTHcF8JaVhxMkHX6Q4bZKGvy9G6oHxQE7a4HMSPQs89pFqrmBeafhnwQg9MHEeFGzDiJa2traamBk7dmy0zWSEztq1a2v5vWjiUiynAsZ9JDjnTJn56Y5gvpLSMOLkgy9SnDZJw98XIxUwHshJG3xOomeBxw7m\/ezww9D55B\/XgnlNW3mDeaXh74oRCg7qMaGeldACBtvry8lG08Zk9+7dMHTo0LRimb7XZKPucEmbw1LWUIogp03S8HdnTuOSuoXkgKS0weckugMcdUXignkHH3kdvj9rL3Q\/txeeevXtqHyeI9XS8LcxohnNTZAt9bxgjAt+qGdFt5CyMi1fec45I52f+ZALU1saRpx88EWM0yZp+PtipB4YD+SkDT4n0bPCk5QfCW\/iXbrmNdi5\/1DUZB4vjDT80zBC8XL\/\/ffDzTffXBdsSz0y2IYG8aYhGe57zjlTNn6GQ9W9JWkYcfLBHZX6kpw2ScPfFyMVMB7ISRt8TqJnhaeI\/EjS8G+EkZ0p3ZQ1nhl6jJoeiab1Jk6cWHfSyByjHjJkCKxcuRJoAC+2Lw0fSfw0+HPaJA3\/rHO4iPLSMOLkgy+enDZJw98XIxUwHshJG3xOonvA0+tmXtMGZqlGL4z5+OZHkoa\/D0acdaThI42fiD2nTdLw5+Sab9vSMOLkgy9GnDZJw98XIxUwHshJG3xOonvA0+tINW4XmdgX2p7vNpI0\/H0w4qwjDR9p\/FQBw8k+t7aVo+k4cc4bafino+FWQoN4HXCSNvicRHeAI7bIJV3PwKafHI13oR\/bC+MTzCsNf1+MuOpJw0ciPzltkoY\/F8\/ytCsNI04++OLEaZM0\/H0xUg+MB3LSBp+T6B7wRFXu\/O4LMP+xN2vVQ+ZHkoa\/L0Zc9aThI5GfnDZJw5+LZ3nalYYRJx98ceK0SRr+vhipgPFATtrgcxLdA56oCto05q7tddV\/sOviXjfzYoGs+ZGk4e+LEVc9afhI5eeIESNYhkAa\/iydzNmoNIyUozkHVEj1Um4h2Sc9aPI7mviuo6MDOjs766Cml4nhXRyNEu+Zijr50tmKC8IN3307in0x2ai\/ctxD8Mnty6H7U09A1w\/fqTWSNZhXGv7paBRbQho++nIodvzL8DTlaPoocc4bafino+FWopQCxr7h1HQV79+YNWsWLFu2DAYOHAjTp0+PBIy5dt2IG3pMNSnxnl7V7kYgUwon355fnwYT73m+VtFsI31nwASY33JMSGYN5q3q5MuGcHJpafhwLsS+mHHaJA1\/X4w460nDiJMPvjhy2iQNf1+M7HqlFDBxt5xix6iwQQFjX+O+Y8eOqP8bN27sldEX\/04T76mAyUYxM\/ns\/EjL37gJxhzeDA9f+VrdkeoswbxVnXzZEFYBkwcvfTnkQS9\/XWlzmJMPvmhx2iQNf1+MKiFg6DYRdmjx4sUwefLkqG9mSwj\/vWLFCjDXuJuOY90kAZMkjKQNPifRfYllbLLzI9nBvKb9LF4Yafj7YsRVTxo+kvnJMQbS8OfoY942pWGkHM07ojLql9IDQ6GjGXrx7yaLL\/7b3kLCvyUJmKRtKayjky+drHRBiMuPtOe4QXDkhu955UeShn86GsWWkIaPvhyKHf8yPE05mj5KnPNGGv7paLiVKL2AMds+ra2tUY+3b99eC9xFMYMfGsgbJ2CSPC8GQhx8\/Gi232RS0Wy\/Vz\/yeu1OmCGn9IPlJz0Eg79\/N+y97nHoWN8P9vz8SNTQmPefAPdOGpTYqGb7dZvE0hYnzoXYDZHepThtkoa\/L0ac9aRhxMkHXxw5bZKGvy9Gdr1SChgUIfjBbSMauLtv376aB8Zk+0VhY7aX4jwwSYn3KFDSBp+T6L7EojbZ+ZE+O\/wwXPiNC7yPVEvD3xcjrnrS8JHOz9DjIA3\/0P0L0Z40jJSjIUa1+W2UUsDYx6hpDEzaMWrqgUH4Z8+eDT09PXUjQZPs4Rc6+dKJai8IccG8Q379Ojz1F89kzo8kDf90NIotIQ0ffTkUO\/5leJpyNH2UOOeNNPzT0XArUUoB49a1cKWkDT4n0X1Rs22iXhhsM09+JGn4+2LEVU8aPmXgZ8ixkIZ\/yL6FaksaRsrRUCPb3HZUwDjgr5MvHaS4BcH2wphWsuZHkoZ\/OhrFlpCGj74cih3\/MjxNOZo+SpzzRhr+6Wi4lVAB44CTtMHnJLoDHLFF4myyj1Sf2XIC7Nx\/LOGjubE37Ui1NPx9MeKqJw2fsvAz1HhIwz9Uv0K2Iw0j5WjI0W1eWypgHLDXyZcOUtKCYB+pTmqpUX4kafino1FsCWn46Muh2PEvw9OUo+mjxDlvpOGfjoZbCRUwDjhJG3xOojvA4eyBwYJmG8l4W4wXxt5GapQfSRr+vhhx1ZOGT5n4GWJMpOEfok+h25CGkXI09Ag3pz0VMA646+RLBylpQbCDeZcM\/j58\/Lk5sKClE3r+d44k82m0jSQN\/3Q0ii0hDR99ORQ7\/mV4mnI0fZQ45400\/NPRcCuhAsYBJ2mDz0l0BzgyeWCoF8ZU7Nk7BfBm3m1XPOqUH0ka\/r4YcdWThk\/Z+Jl3XKThn7c\/HPWlYaQc5Rjl4ttUAeOAuU6+dJAaLQhJ+ZE6Tv8qbDp+NKQF80rDPx2NYktIw0dfDsWOfxmephxNHyXOeSMN\/3Q03EqogHHASdrgcxLdAY7MHhisYAfz\/mDXxc4380rD3xcjrnrS8CkjP\/OMjTT88\/SFq640jJSjXCNdbLsqYBzw1smXDlLagmAH887bvwwmHvgX6P7UE9D1w3dqD4gL5pWGfzoaxZaQhk8aF4pF5+jTOG2Shn8z8E17pjSMOPmQhkXS95w2ScPfFyO7ngoYBySlDT4n0R3g8PLA5AnmlYa\/L0Zc9aThU0Z+5hkbafjn6QtXXWkYKUe5RrrYdlXAOOCtky8dJJcFwb6ZF4N5Bx95HR6+8rWGwbzS8E9Ho9gS0vBx4UKxCKkHpmi87ecpR9NHgHPeSMM\/HQ23EipgHHCSNvicRHeAw8sDg5XSgnlNw\/aRamn4Z8Fo6dKlUfHOzs7ov88++yxMmTIl+jdNQkoTlI4ePRpWrFgBLS0tUTmaoNRONIrfS8OnrPzMMq60rDT8ffvBWU8aRspRztEurm0VMA5Y6+RLB8l1QYgL5t1yxqfgUPvXofu5vfDUq29HD+u59hy48EOnRv+Whn86i2sNVwAAFl9JREFUGkdLGLHS0dERCZj9+\/fDzJkzYe7cudH3ixYtgq6urkiooEjZvn17VA5Fz\/Dhw2Hy5Mnwyiuv1Mpt3boVuru7YcmSJdC\/f\/+aGdLwceWCK44hynHaJA3\/EHiFbkMaRpx88MWO0yZp+PtiZNdTAeOApLTB5yS6AxyxRVxtottIeCvvinfuhNN+9CDsXfwWXPfgS7VcSdQLIw1\/F4xQrCxYsABGjhwJBw4ciIQJChoUJ+hdQQEye\/ZsaG9vB\/S44L9bW1sj0YLljFBZvXo1bNy4MRIt6KUxAgjbNR9p+LhywQXHUGU4bZKGfyjMQrYjDSNOPvjixmmTNPx9MVIB44GctMHnJLoHPFEVV5toMC+KlDNaToi2luI+Jj+SNPxdMEKhMm7cONixY0fNs0KFCbZhREtbW1tNzIwdO7ZO6Kxdu7ZWH0XR9OnTIzGE5VTAuIzE0TKu\/HRv8VjJMvLTp5956kjDiJMPvjhx2iQNf1+MVMB4ICdt8DmJ7gFP5heEHcxrnpmUH0ka\/mkYoVDZsGFDJDTo1pAKmDTk+L7nnDNl4ycfysktS8OIkw+++HLaJA1\/X4xUwHggJ23wOYnuAU9mAWMfqUZPjIl9oc8320jS8LcxQm\/L8uXLoz9jkC2KF\/P\/pizGwaBHhmsLCZ+zbt063+ELWm\/37t0wdOjQoG3mbYzDpvHjx9fM2rZtW14TK11f2hwu+xqalSzS8M9qf1J5jYFxQFLa4Fdh8sV5YTAmpv38wb2OVF\/xiXOhrC8I6oHRIF6HycZUhHPOSFsfmCDM1aw0jDj54AsUp03S8PfFqBIeGHrkFDtEj53SI6fm9AftNHXjm5Mcekw1P52yTj77SDWKl7tfnAQd77sd9vQbVJcf6d+6\/rQSAgZRpseo6ZFoyumJEyfWnTQy\/BwyZAisXLkyCgymH2mLU1Yu5GdfegucNknDPx2N4ktIw4iTD77octokDX9fjCohYOgvWbqY45HTWbNmwbJly2DgwIG9Ah7Ni4C+IPSYahgq+Uw+eqR6zOEtsPyNG2PzI53y\/3bC9n97LoyhFWxF2uLkwwXuYeG0SRr+3Fj6tC8NI04++OCDdThtkoa\/L0aVEDBUdJjLvrBjVNiggKFHTvGXL54IwY85looeGBQ1ekw1P518Jp+dH2n5GzfBmMObe93Me+IP\/x52r\/9mfiMr2oK0xcmHC9xDw2mTNPy5sfRpXxpGnHzwwUcFjB9qpYyBoVs+2G16o6k5aop\/p7eZGnioYDECxlwgpsdU\/UjkO\/lc8yP1e\/NleOP\/+Ut\/4ypeU18O6QPM+cKShn86GsWXkIYRJx980eW0SRr+vhhVwgNDO0FFB\/7dnPLAf8fdmaECJhR16tvxnXx2MO8Pdl0MeDPv1gl\/0zA\/Ek8vytmqtMXJlwuc6HPaJA1\/Thx925aGEScffDHitEka\/r4YVU7AmOBHvMUUP8abYsQM\/tfkocF\/xwkYly0krKvHVJNp53tMteelX8L8x96sNbz0nb+H8bsfgLbBq6JgXvOx8yOFmgBVaEfa4sS5EPuOF6dN0vD3xYiznjSMOPngiyOnTdLw98WoEgIGRQh+TK4YE7i7b9++XvdsmOvZTcdtAaNBvGGolGfyxQXzpuVHCmN1NVqRtjjl4QLXiHDaJA1\/LgzztCsNI04++OLEaZM0\/H0xqoSAsY9R0xiYtGPUtoAxXpk5c+aAHlP1p1WeyWfnR\/rekblw8MX10G\/5QWi75\/nY\/Ej+llavprTFKQ8XuEaH0yZp+HNhmKddaRhx8sEXJ06bpOHvi1ElBEyozru2I23wOYnuioldLo9Ndn6ka4\/fAKPW3AgLWjqhZ8CEukeZ\/Ei+dlaxnvIzfVTz8DOtdWn4p9nbjO+lYcTJB198OW2Shr8vRipgPJCTNvicRPeAJ6qS1yY7mLdn7xTYc9wg2HbFo3XBvHe3nxVdcqefYwgoP9PZkJefjZ4gDf90NIovIQ0jTj74ostpkzT8fTFSAeOBnLTB5yS6BzxBBIx9pHoF3AHT4YZe5mgwb+8RUn6ms5ZzzkjDPx2N4ktIw4iTD77octokDX9fjFTAeCAnbfA5ie4BTxABg40k5Ud6Y+M\/wqHfa6uZ1nPtOXDhh071NbVy9ZSf6UPKOWek4Z+ORvElpGHEyQdfdDltkoa\/L0YqYDyQkzb4nET3gCeYgInLj7Rz\/6GofcyVdOEHTwUso16Y+lFSfqazlnPOSMM\/HY3iS0jDiJMPvuhy2iQNf1+MVMB4ICdt8DmJ7gFPMAGDDdEj1Y1s0WDeY+goP9NZyzlnpOGfjkbxJaRhxMkHX3Q5bZKGvy9GKmA8kJM2+JxE94AnqICxj1SjB+bdv3oTJo\/7PyLvi\/noNpIKmCxc5Zwz0taHLLgUVVYaRpx88MWU0yZp+PtipALGAzlpg89JdA94ggoYGsyLp4127T8ET736dvQM3EYyW0oqYFTAZOEq55yRtj5kwaWostIw4uSDL6acNknD3xcjFTAeyEkbfE6ie8ATVMCgQDn7K8\/UmYH5kSaN3gh3X34WXPfgS5GI0S0kFTBZuMo5Z6StD1lwKaqsNIw4+eCLKadN0vD3xUgFjAdy0gafk+ge8AQVMPZxagzYfd\/jX4ZRV39VkzsmDI7yM521nHNGGv7paBRfQhpGnHzwRZfTJmn4+2KkAsYRud\/8ajcceP6LcPKF3YCD\/+N\/+19weNc\/wW8NHAv93jvWsRWeYpxE97U4pE0YyIvbR0+\/+nZtywjtoieRdAtJPTBZuBqSn\/Zzq\/pyyIJvWllpGHHyIQ2LpO85bZKGvy9GKmAckfvF0+3wX\/uejQTLs99\/Fs7\/w6GAogY\/7z5xaGIrWJ5+jj\/jT4MLHk6iO8LTq1hIm0wgL4qUM1pOgI988tPQvao7+jfmRsIPbiHp5ygC0hankFwINcacNknDPxRmIduRhhEnH3xx47RJGv6+GKmAyYCcETEZqsQWbSR4sEJW0cNJdN++hrTJbCOhx6X9\/MFw1+K\/huvn\/N\/RFhL+DWNh9CK7YyMlbXEKyQVfPtr1OG2Shn8ozEK2Iw0jTj744sZpkzT8fTFSAZMBOfS4\/Hxje83zkqFq7qKNRM+h4\/8ATj7p5NozOLw8WTsQevJhoC7NRI32oHhZfe050X\/1cxQBxMl4qDgw8RGKobkQol+cNlX15RACd9OGNIw4+eCLG6dN0vD3xUgFjCNyKF4w5uXgy3c41mhusSTBk9W749sLrslHX9A+L1Pf\/kivh7h8vvul2hFzTntRMHZOGAFnnuYmHPe+vhcGDxrMaVLmtjltap\/SDjuf+15mm\/pSBWkvUK71Ks+YctokDf88ONG673rnnXfeCdVYldo5+PLflEa8ZMWdI4ZHJ1\/WUfAvj+Ll0vmPwrz9y+oa2XTCaBhyhFz2N+CTsOn40f4P0pqZENDtzWS4pL1AOderTKQhhTltkoa\/L0bqgXFE7sibz0bbR\/pxC1r+xS9\/EW1rcWxnVXXy+XILg5yveOZzMObw5tQm9vYblFpGC\/gjYAvEngGfhPk3\/oXGaBFIubc5fUYvySPXTC+vCpjsI6kemAaYoYjBLSQ8hTT2j46eLvr1waMnkeI+5pRS9mGoVo3QHh4VMMf4QS\/6W\/7GTU4iplrskt2b7wyYAI+dvxhWX3eObEMLsK7Ibc6Q3cm6ZRrq2brNmR1JFTAOmLm+QI3goU02EjxYTkVPuofnW9\/6J\/jTP\/1MBCuHh8eBAmKK0Iv+cLto+U9vgsFk20iMoX3QEBQv81s6o5739Zui427U7oOUENflqm1zllLAHDx4EGbPng09PT0RQUaPHg0rVqyAlpYWeOihh2DOnDnR3zs6OqCz8+iCQuvQ8vjd0qVLYfny5VG5xYsXw+TJk+uI5ypgfNiaVfT0dcFz\/BmfgQHn3OoDdaF1nn32WZgyZUovHtK\/U6414ifl9KpVq+DSfzwYxbpc9qs1cPXP7i+0X\/qweARwq27i4FW1L\/vaRYs2RxdtOh62vfzvveK0hvz6WIyWcqkYBKq8zVlKAbN\/\/36YOXMmzJ07F0aOHFljwSuvvAKzZs2CZcuWwcCBA2H69OmRgBk7dmwkbLZv3x79PwqW4cOHR0IFXyj4\/yiA9u3bV6tP2+UUMFkpjIJn\/5al0P+EYydCkrw8VRM7ZREvyMNFixZBV1cX9O\/fHxYuXAhTp06NOGl4i+NuyhjhHcdP2tbWrVuhu7sbfvyh\/wv+zy13qXjJOnmYy286\/mzoOP326Cl9ScDYHL2v+9vw6HGfAN3iZCacZ\/NV2uYspYChEwYXf\/Ohwoa+LIYOHRp5bFpbW2uiBV8ES5Ysgd27d9deJChg6EvFtCtJwKBNrsFeVfLuYFzNqR9\/0nPKFlsNxTJ+bE8eFcsobJCT7e3tkQcxiZ+rV6+GjRs3RlxFLw0KoE\/+xUz4u398rE7AuAT0FotC33wavhz+\/kPz+tRN0ch3ytGpnbfC0yd\/IiKAihhZ86Bq25ylFDDUXYn0oK54FDHoecGP2VYy7nl8WaA3hr5IUAChIJo2bRqMGTMmelHgy4V+pAkYLnukCx680wZzU0n\/ID9feOEFWL9+PezZs6e2lYm8M8IZ+2BES1tbW03M2Pxcu3ZtzXNouI1exCGjzq676G\/M4S29PDLqrudnyp7j6k95\/cuIDvjc1Mv71Ckk6t02HN0wYkYEvsZp8XPQ9QlV3OYspYChA0YXdfy72Q7Cf5stJPMLN07A4AvC\/HpAb4zZgrK3kFxJ0pfKjfkQwNWfrL9GaMgxh1gvKAY3+M4VtzJsIyEHUbgYr4nhIfYxlIBBoYMfc0TVFT\/Xcr8693PwmxPf61pcy\/Xhm6LjBMz+D\/8V7Nu7R+O0hM2Mqm1zll7AGO8Kbg\/hx8QRGDGD\/50xY0aii\/7OO++sxcPQtmz3vzAeltKcrB4e7ORx\/esTZ\/b\/3RuCJ8fMCyYNAscg2x07dtR4SDk1bNiwmsDOs4Vkx37ltT+pPoqjXfsPcTVfuXabeYdIM8G0t5DMNufz37xN47SaOTAJz67SNmcpBQyNMaCBuxjDYjww5gVh4l6SgnjtyUcDfwVyT00qAQI0Rot6AkeNGhUkiDdum7MEsKiJFUUgLtAcOfqubd+HnsU31vVatzWLJ0GVtzlLKWDsY9Q0BsblGPXEiRPrYl3SjlEXTzl9YtkRoDyk\/KTHqNFbY7aCKKdtfpq2hgwZAitXrqw7eVd2nNT+aiDQiKNl8ORd9+BL0XZs1T9VS4hbSgFTdZJp\/xQBRUARUASKRaAMQisvIlXb5lQB04AR9uVM5tdyXhJlqd\/o0r6kS9GytJ+lLD1FY05qJWFkgqu3bNkCtkchyzO1bDICys96bJSf8maLclQ5yslKFTAJ6Cbt69pHrDkHB9tOurSP\/h3Lxd1fE9I2sxBRMdIII3NZoH1EOKRNfbkt5Wf96Cs\/5c0G5ahylJuVKmASEI6LrC\/q9Ac1KenSvqRL0Ti8RPgsPFmDH3PkHIVcEkaNbkHmJnRfaV\/5eWyklZ8yWa8cVY5yM1MFTAMBY45k07tmOARCo0FOurQv6VI0zuPfdEEyAiYOI3raBu\/Tsetxk7ovtB9394ZJm1Fk\/5WfRaJdrmcpR3uPl66hYTmsAka4gKHm2Zf2xV2KpgIm7ASR2pqUl4PyUypDmm+XclQFDDcLVcAI30Ki5rlcisbpIYr79WDn6cFtNt1C4p62kLh9R2+Q5rei\/gnKz6IRl\/08KVtIuobK5kke61TAJKAnJQAt6dK+RpmN8xAibTuLxsBoEC8X0untKj\/Tf90qP9N5xFlCOaoc5eQXtq0CpgHCEi4Qa3RpX9KlaFykiYtlScKIHqPu6OgAjM\/QT1gElJ\/1eCo\/w\/IrRGvKUeVoCB4ltaEChhNdbVsRUAQUAUVAEVAEWBBQAcMCqzaqCCgCioAioAgoApwIqIDhRFfbVgQUAUVAEVAEFAEWBFTAsMCqjSoCioAioAgoAooAJwIqYDjR1bYVAUVAEVAEFAFFgAUBFTAesNLTP3b1VatWwYYNG2D48OHAealcktmYg2j58uVRAsVJkybBtGnTYMiQIbBy5UrIekdIUhoD+mx83rhx4yDE\/TMh2\/IY1spUUX4eG8qQnArZVmXI5tkR5ahy1JM6ddVUwOREERc1\/Eg4JoxHrhcuXAhTp06NxAq9Q8anm2kCBhchFGuh+o5HrxcsWADz5s2DlpYWH5O1joWA8lP5KX1SKEeVo74cVQHji9x\/14ubfHYm5rPPPhtWrFgBe\/bsgcWLFwPmD0IvyejRo6O\/m5c1\/VVCsz7bJqKwQM8KtmfawNxEs2fPhp6enqj41VdfDffee2\/0b7yHZcaMGXXfox3GQxTXHtabPn06bNmypZed+J0tlozYGTNmDNx+++2R1+e2224DxALboHfB0H7aGOQVXTmHs3LVlZ9HxbzyUy61laPKUV92qoDxRS6DgMGiS5YsiV7kU6ZMiURMW1tbJChaW1sjIYEL7KxZs2DZsmUwdOjQ6DsUAbZ3w1wQd\/nll0f1cPKjkMH28YP12tvboy0dmovEzktiPB1GqJhEgLS93bt3w6JFi6Crq6uXRwTtvf\/+++Hmm28GFE9GBF1\/\/fU1u1BM4dYVfkzf6A3CcV4iu92cw9Pnq7u8HJSfys9mThTlqK6hvvxTAeOLXAYBQ0UKFQTGU4NCxL5FlGabRoFgPva2TpzwSRMwtMv4HLTDeIJo+\/v27UsUMPb2EbXDCBOaqXrmzJlg8iSZf8fF5Og2Uk5CWtVdXg7Kz\/2g\/AzLuyytKUeP3lKua2gW1hwtqwImO2Z1NVwmnxEUtviwBcycOXPq2ra3V\/BLW3DgC98svsZzEydgsK4J8MV\/Y7AxemnigulM0C+WS\/LAUI+OmXy0rO3xoaKFblnZW2UqYHIS0kPAKD+PzSGz3WS2aJWfYfkY15quoccEjK6h2fimAiYbXr1Kh5x8xmPRyCRfDwxtk4oe9LJ0d3dHW1DU0xMnSmzPDQ3gte1qJGBoO7YQUgGTk5CMAkb5CaD8DMtP88MK\/0u3y+04QleRrRztWxxVAZNzPoYSMLb7kMaiUGHhGwNDvT10EcbuY7CuiYFBQfHggw9GW0qNtpDiYmBcfj3YXh07aFdjYHISkknAKD8nR8gqP8PyM6SAUY72PY6qgMk5H0MJGLM9hEG++InbPjKmxp0awpNMJnN13BYSzQ6N7ZgtJLMoG5c5vTPG1MEy9LQU\/n\/SKSQT8NvIA2My1GI7toteTyHlJCSTgFF+Ho1DU36G5WdIAaMc7XscVQETfj72mRb1Hpg+M9Sl7Kjys5TD1qeMVo7mG24VMPnw6\/O1Q95OGrKtPj8wCkCEQEhOhWxLh0cRMAiE5FXItsowQipgyjBKaqMioAgoAoqAIqAI1CGgAkYJoQgoAoqAIqAIKAKlQ0AFTOmGTA1WBBQBRUARUAQUARUwygFFQBFQBBQBRUARKB0CKmBKN2RqsCKgCCgCioAioAiogFEOKAKKgCKgCCgCikDpEFABU7ohU4MVAUVAEVAEFAFFQAWMckARUAQUAUVAEVAESoeACpjSDZkarAgoAoqAIqAIKAIqYJQDioAioAgoAoqAIlA6BFTAlG7I1GBFQBFQBBQBRUARUAGjHFAEFAFFQBFQBBSB0iHw\/wMnsDFSy6WexwAAAABJRU5ErkJggg==","height":337,"width":560}}
%---
