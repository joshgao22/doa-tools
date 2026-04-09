%[text] DOASTATDUALSATURAECIPERF
%[text] Static single-frame dual-satellite DoA-only performance check.
%[text] This script compares single-satellite and dual-satellite pilot-based
%[text] DoA MLE under a static single-frame model in ECI geometry. The signal
%[text] generation explicitly removes Doppler phase so that the Monte Carlo
%[text] snapshots are consistent with estimatorDoaMlePilotOpt and crbDetDoa.
%[text] The reference satellite is selected by user elevation, and the single-
%[text] satellite baseline always uses that same reference satellite.
clear(); close all;
%[text] ## Parameters
numUsr = 1;
numSym = 64;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 18e9;
wavelen = 299792458 / carrierFreq;
rngSeed = 253;
rng(rngSeed);

elemSpace = wavelen / 2;
numElem = [4 4];
pwrSource = 1;
pathGain = 1;
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);
truthLatlon = usrLla(1:2, 1);

utc = datetime([2026, 03, 12, 08, 53, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(localResolveTlePath("starlink_20260312.tle"));
arrUpa = createUpa(numElem, elemSpace);
minUsrElevationDeg = 15;
maxSatOffAxisDeg = 55;
numPickSat = 2;
searchMarginDeg = 5;
gridSize = [50 50];
snrDb = (-20:3:10).';
numRepeat = 200;
showTable = true;
%%
%[text] ## Multi-Satellite Scene With Elevation-Based Satellite Reference
[~, satAccess] = findVisibleSatFromTle( ...
  utc, tle, usrLla, minUsrElevationDeg, maxSatOffAxisDeg);
[satIdxSel, satPickAux] = pickVisibleSatByElevation( ...
  satAccess, numPickSat, 1, "available");
refSatIdxGlobal = satIdxSel(1);

scene = genMultiSatScene( ...
  utc, tle, usrLla, satIdxSel, [], arrUpa, ...
  minUsrElevationDeg, maxSatOffAxisDeg, "satellite", refSatIdxGlobal);
refSatIdxLocal = localGetReferenceSatIdx(scene);
sceneSingle = selectSatScene(scene, refSatIdxLocal);

steeringInfo = getSceneSteering(scene, wavelen);
linkParam = getLinkParam(scene, wavelen);
truthLocalDoaCell = localBuildLocalDoaCell(scene);
localDoaJacCell = localBuildLatlonToLocalDoaJac( ...
  scene, utc, tle, usrLla, arrUpa, satIdxSel, refSatIdxGlobal, ...
  minUsrElevationDeg, maxSatOffAxisDeg);

if showTable %[output:group:9d236d8b]
  fprintf('\n========== Selected satellites ==========%s', newline); %[output:64882147]
  disp(table((1:scene.numSat).', satIdxSel(:), ... %[output:5f44fb45]
    reshape(scene.accessInfo.usrElevationDeg(:, 1), [], 1), ... %[output:5f44fb45]
    'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg'})); %[output:5f44fb45]

  fprintf('\n========== Reference satellite ==========%s', newline); %[output:054f418e]
  disp(table(refSatIdxLocal, refSatIdxGlobal, ... %[output:080025cc]
    'VariableNames', {'refSatIdxLocal', 'refSatIdxGlobal'})); %[output:080025cc]
end %[output:group:9d236d8b]
%%
%[text] ## Pilot And Zero-Doppler DoA-Only Signal Model
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, ~] = genPilotWaveform( ...
  pilotSym, symbolRate, osf, 'rrc', pulseOpt);
numSnap = numel(pilotWave);

pathGainMat = localBuildPathGainMat(scene.numSat, numUsr, pathGain);
modelOpt = struct();
modelOpt.delayMode = 'phaseonly';
modelOpt.carrierPhaseMode = 'none';
modelOpt.precomp.phaseSeqCell = ...
  localBuildUnitPhaseSeqCell(scene.numSat, numUsr, numSnap);

searchRange = [ ...
  truthLatlon(1) + [-searchMarginDeg, searchMarginDeg]; ...
  truthLatlon(2) + [-searchMarginDeg, searchMarginDeg]];
doaGridSingleCell = localBuildDoaGridCell(sceneSingle, gridSize, searchRange, E);
doaGridSingle = doaGridSingleCell{1};
doaGridJoint = localBuildDoaGridCell(scene, gridSize, searchRange, E);
%%
%[text] ## Monte Carlo Loop
numSnr = numel(snrDb);
mseSingleAng = zeros(numSnr, 1);
mseJointAng = zeros(numSnr, 1);
crbSingleDeg = zeros(numSnr, 1);
crbJointDeg = zeros(numSnr, 1);
idealJointCrbDeg = zeros(numSnr, 1);

crbSingleOpt = struct();
crbSingleOpt.localDoaJac = localDoaJacCell{refSatIdxLocal};
crbSingleOpt.paramName = ["lat", "lon"];

crbJointOpt = struct();
crbJointOpt.localDoaJac = localDoaJacCell;
crbJointOpt.paramName = ["lat", "lon"];

singlePwrSource = abs(pathGainMat(refSatIdxLocal, 1))^2 * pwrSource;
jointPwrSource = cell(1, scene.numSat);
for iSat = 1:scene.numSat
  jointPwrSource{iSat} = abs(pathGainMat(iSat, 1))^2 * pwrSource;
end

progressbar('reset', numSnr);
for iSnr = 1:numSnr %[output:group:0a2d994f]
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));
  seSingleAng = zeros(numRepeat, 1);
  seJointAng = zeros(numRepeat, 1);

  parfor iRepeat = 1:numRepeat
    [rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
      steeringInfo, pilotWave, linkParam, ...
      carrierFreq, sampleRate, pwrNoise, pathGainMat, modelOpt);

    rxSigSingle = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');

    estSingle = estimatorDoaMlePilotOpt( ...
      sceneSingle.array{1}, wavelen, rxSigSingle, pilotWave, doaGridSingle);
    estJoint = estimatorDoaMlePilotOpt( ...
      scene.array, wavelen, rxSig, pilotWave, doaGridJoint);

    singleLatlon = localGetLatlonEstimate(estSingle);
    jointLatlon = localGetLatlonEstimate(estJoint);

    seSingleAng(iRepeat) = deg2rad(calcLatlonAngleError(singleLatlon, truthLatlon))^2;
    seJointAng(iRepeat) = deg2rad(calcLatlonAngleError(jointLatlon, truthLatlon))^2;
  end

  mseSingleAng(iSnr) = mean(seSingleAng);
  mseJointAng(iSnr) = mean(seJointAng);

  [singleCrb, ~] = crbDetDoa(sceneSingle.array{1}, wavelen, ...
    truthLocalDoaCell{refSatIdxLocal}, singlePwrSource, pwrNoise, numSnap, crbSingleOpt);
  [jointCrb, ~] = crbDetDoa(scene.array, wavelen, ...
    truthLocalDoaCell, jointPwrSource, pwrNoise, numSnap, crbJointOpt);

  crbSingleDeg(iSnr) = projectCrbToAngleMetric(singleCrb, truthLatlon, 'latlon');
  crbJointDeg(iSnr) = projectCrbToAngleMetric(jointCrb, truthLatlon, 'latlon');
  idealJointCrbDeg(iSnr) = crbSingleDeg(iSnr) / sqrt(2);

  progressbar('advance'); %[output:8fe1bd02]
end %[output:group:0a2d994f]
%%
%[text] ## Metrics And Plot
rmseSingleDeg = rad2deg(sqrt(mseSingleAng));
rmseJointDeg = rad2deg(sqrt(mseJointAng));
idealJointRmseDeg = rmseSingleDeg / sqrt(2);

rmseGainDb = 20 * log10(rmseSingleDeg ./ rmseJointDeg);
crbGainDb = 20 * log10(crbSingleDeg ./ crbJointDeg);
crbGapToIdealDb = 20 * log10(idealJointCrbDeg ./ crbJointDeg);

result = struct();
result.cfg = struct();
result.cfg.numUsr = numUsr;
result.cfg.numSym = numSym;
result.cfg.sampleRate = sampleRate;
result.cfg.symbolRate = symbolRate;
result.cfg.carrierFreq = carrierFreq;
result.cfg.rngSeed = rngSeed;
result.cfg.numElem = numElem;
result.cfg.pwrSource = pwrSource;
result.cfg.pathGain = pathGain;
result.cfg.usrLla = usrLla;
result.cfg.utc = utc;
result.cfg.minUsrElevationDeg = minUsrElevationDeg;
result.cfg.maxSatOffAxisDeg = maxSatOffAxisDeg;
result.cfg.numPickSat = numPickSat;
result.cfg.searchMarginDeg = searchMarginDeg;
result.cfg.gridSize = gridSize;
result.cfg.snrDb = snrDb;
result.cfg.numRepeat = numRepeat;

result.truth = struct();
result.truth.latlon = truthLatlon;
result.truth.selectedSatIdxGlobal = reshape(satIdxSel, 1, []);
result.truth.refSatIdxLocal = refSatIdxLocal;
result.truth.refSatIdxGlobal = refSatIdxGlobal;
result.truth.usrElevationDeg = reshape(scene.accessInfo.usrElevationDeg(:, 1), 1, []);
result.truth.pickAux = satPickAux;

result.metric = struct();
result.metric.snrDb = snrDb;
result.metric.numSnap = numSnap;
result.metric.rmseSingleDeg = rmseSingleDeg;
result.metric.rmseJointDeg = rmseJointDeg;
result.metric.idealJointRmseDeg = idealJointRmseDeg;
result.metric.crbSingleDeg = crbSingleDeg;
result.metric.crbJointDeg = crbJointDeg;
result.metric.idealJointCrbDeg = idealJointCrbDeg;
result.metric.rmseGainDb = rmseGainDb;
result.metric.crbGainDb = crbGainDb;
result.metric.crbGapToIdealDb = crbGapToIdealDb;

result.aux = struct();
result.aux.scene = scene;
result.aux.sceneSingle = sceneSingle;
result.aux.localDoaJacCell = localDoaJacCell;
result.aux.pathGainMat = pathGainMat;
result.aux.searchRange = searchRange;

figure(); %[output:0e74613e]
subplot(1, 2, 1); %[output:0e74613e]
semilogy(snrDb, rmseSingleDeg, '-o', ... %[output:0e74613e]
  snrDb, rmseJointDeg, '-s', ... %[output:0e74613e]
  snrDb, crbSingleDeg, '--o', ... %[output:0e74613e]
  snrDb, crbJointDeg, '--s', ... %[output:0e74613e]
  snrDb, idealJointCrbDeg, ':', ... %[output:0e74613e]
  snrDb, idealJointRmseDeg, '-.'); %[output:0e74613e]
grid on; %[output:0e74613e]
xlabel('SNR (dB)'); %[output:0e74613e]
ylabel('Global angle RMSE / CRB (deg)'); %[output:0e74613e]
title(sprintf('Static DoA-only (%d snaps)', numSnap)); %[output:0e74613e]
legend({ ... %[output:0e74613e]
  'SS-SF DoA-only RMSE', ... %[output:0e74613e]
  'MS-SF DoA-only RMSE', ... %[output:0e74613e]
  'SS-SF DoA-only CRB', ... %[output:0e74613e]
  'MS-SF DoA-only CRB', ... %[output:0e74613e]
  'Ideal MS CRB (-3 dB)', ... %[output:0e74613e]
  'Ideal MS RMSE (-3 dB)'}, ... %[output:0e74613e]
  'Location', 'southwest'); %[output:0e74613e]

subplot(1, 2, 2); %[output:0e74613e]
plot(snrDb, rmseGainDb, '-s', ... %[output:0e74613e]
  snrDb, crbGainDb, '--o'); %[output:0e74613e]
grid on; %[output:0e74613e]
hold on; %[output:0e74613e]
yline(20 * log10(sqrt(2)), ':'); %[output:0e74613e]
hold off; %[output:0e74613e]
xlabel('SNR (dB)'); %[output:0e74613e]
ylabel('Joint gain over single (dB)'); %[output:0e74613e]
title('Dual-satellite gain'); %[output:0e74613e]
legend({ ... %[output:group:006d945c] %[output:0e74613e]
  'RMSE gain', ... %[output:0e74613e]
  'CRB gain', ... %[output:0e74613e]
  'Ideal gain (3 dB)'}, ...
  'Location', 'southeast'); %[output:group:006d945c] %[output:0e74613e]
%%
% saveExpSnapshot("doaStatDualSatUraEciPerf");
%%
function localDoaCell = localBuildLocalDoaCell(scene)
%LOCALBUILDLOCALDOACELL Convert scene.localDoa into a 1xNs cell array.

numSat = scene.numSat;
localDoaCell = cell(1, numSat);
for iSat = 1:numSat
  localDoaCell{iSat} = reshape(scene.localDoa(:, iSat), 2, 1);
end
end


function jacCell = localBuildLatlonToLocalDoaJac(scene, utc, tle, usrLla, arrUpa, ...
  satIdxSel, refSatIdxGlobal, minUsrElevationDeg, maxSatOffAxisDeg)
%LOCALBUILDLATLONTOLOCALDOAJAC Build numerical Jacobian for each satellite.

numSat = scene.numSat;
baseLocalDoa = scene.localDoa;
deltaDeg = 1e-6;
jacCell = cell(1, numSat);
for iSat = 1:numSat
  jacCell{iSat} = zeros(2, 2);
end

for iParam = 1:2
  usrPert = usrLla;
  usrPert(iParam, 1) = usrPert(iParam, 1) + deltaDeg;
  scenePert = genMultiSatScene( ...
    utc, tle, usrPert, satIdxSel, [], arrUpa, ...
    minUsrElevationDeg, maxSatOffAxisDeg, 'satellite', refSatIdxGlobal);
  deltaLocalDoa = (scenePert.localDoa - baseLocalDoa) / deltaDeg;

  for iSat = 1:numSat
    jacCell{iSat}(:, iParam) = reshape(deltaLocalDoa(:, iSat), 2, 1);
  end
end
end


function refSatIdxLocal = localGetReferenceSatIdx(scene)
%LOCALGETREFERENCESATIDX Resolve the local reference satellite index.

if ~isfield(scene, 'ref') || ~isfield(scene.ref, 'weight') || isempty(scene.ref.weight)
  error('doaStatDualSatUraEciPerf:MissingReferenceWeight', ...
    'scene.ref.weight is required to resolve the reference satellite.');
end

refWeight = scene.ref.weight(:);
refSatIdxLocal = find(refWeight > 0.5, 1, 'first');
if isempty(refSatIdxLocal)
  error('doaStatDualSatUraEciPerf:InvalidReferenceWeight', ...
    'The scene does not contain a unique reference satellite.');
end
end


function doaGridCell = localBuildDoaGridCell(scene, gridSize, searchRange, E)
%LOCALBUILDDOAGRIDCELL Build one global lat-lon grid per satellite.

numSat = scene.numSat;
doaGridCell = cell(1, numSat);
for iSat = 1:numSat
  doaGridCell{iSat} = genDoaGrid('latlon', 2, gridSize, searchRange, ...
    'eci', datevec(scene.utc), scene.satPosEci(:, iSat), scene.rotMat{iSat}, E);
end
end


function phaseSeqCell = localBuildUnitPhaseSeqCell(numSat, numUsr, numSnap)
%LOCALBUILDUNITPHASESEQCELL Build per-satellite unit phase sequences.

phaseSeqCell = cell(1, numSat);
unitSeq = ones(numUsr, numSnap);
for iSat = 1:numSat
  phaseSeqCell{iSat} = unitSeq;
end
end


function pathGainMat = localBuildPathGainMat(numSat, numUsr, pathGain)
%LOCALBUILDPATHGAINMAT Expand pathGain into a numSat x numUsr matrix.

if isscalar(pathGain)
  pathGainMat = repmat(pathGain, numSat, numUsr);
  return;
end

if isvector(pathGain)
  if numel(pathGain) == numSat
    pathGainMat = repmat(reshape(pathGain, numSat, 1), 1, numUsr);
  elseif numel(pathGain) == numUsr
    pathGainMat = repmat(reshape(pathGain, 1, numUsr), numSat, 1);
  else
    error('doaStatDualSatUraEciPerf:InvalidPathGainVector', ...
      'pathGain vector length must match numSat or numUsr.');
  end
elseif isequal(size(pathGain), [numSat, numUsr])
  pathGainMat = pathGain;
else
  error('doaStatDualSatUraEciPerf:InvalidPathGainSize', ...
    'pathGain must be scalar, numSat vector, numUsr vector, or numSat x numUsr matrix.');
end
end


function latlon = localGetLatlonEstimate(estResult)
%LOCALGETLATLONESTIMATE Read global DoA estimate from estimator output.

if isfield(estResult, 'doaParamEst') && ~isempty(estResult.doaParamEst)
  latlon = estResult.doaParamEst;
elseif isfield(estResult, 'latlonEst') && ~isempty(estResult.latlonEst)
  latlon = estResult.latlonEst;
else
  error('doaStatDualSatUraEciPerf:MissingDoaField', ...
    'Estimator output does not contain doaParamEst or latlonEst.');
end

latlon = reshape(latlon(:, 1), 2, 1);
end


function tlePath = localResolveTlePath(fileName)
%LOCALRESOLVETLEPATH Resolve one TLE file from common test locations.

candidateList = { ...
  fullfile('/tle', fileName), ...
  fullfile('tle', fileName), ...
  fullfile('test', 'data', 'tle', fileName), ...
  fullfile('..', 'data', 'tle', fileName)};

for iPath = 1:numel(candidateList)
  if exist(candidateList{iPath}, 'file')
    tlePath = candidateList{iPath};
    return;
  end
end

error('doaStatDualSatUraEciPerf:TleFileNotFound', ...
  'Unable to locate the TLE file "%s".', fileName);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:64882147]
%   data: {"dataType":"text","outputData":{"text":"\n========== Selected satellites ==========\n","truncated":false}}
%---
%[output:5f44fb45]
%   data: {"dataType":"text","outputData":{"text":"    localSatIdx    globalSatIdx    usrElevationDeg\n    ___________    ____________    _______________\n\n         1              2              71.717     \n         2              1              42.036     \n\n","truncated":false}}
%---
%[output:054f418e]
%   data: {"dataType":"text","outputData":{"text":"\n========== Reference satellite ==========\n","truncated":false}}
%---
%[output:080025cc]
%   data: {"dataType":"text","outputData":{"text":"    refSatIdxLocal    refSatIdxGlobal\n    ______________    _______________\n\n          1                  2       \n\n","truncated":false}}
%---
%[output:8fe1bd02]
%   data: {"dataType":"text","outputData":{"text":"[====================] (11\/11) 100.0% [2 m 0 s]","truncated":false}}
%---
%[output:0e74613e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7svQmYFsW1Pn5UVIZFAUURWQYVFLkqqEQc5S9eYohGwDUwZEFAg4piruxoZEhkFbiIokFFNEYGLkYjYzDEoPhDR4wLEhdcYVgcUGQRkEFB\/T+nSH3WV9PdVdV96uv6hurn8ZH5+nTVqfecqn67lnMO+v77778Hf3kEPAIeAY+AR8Aj4BHIIwQO8gQmj6zlVfUIeAQ8Ah4Bj4BHgCHgCYx3BI+AR8Aj4BHwCHgE8g4BT2DyzmReYY+AR8Aj4BHwCHgEPIHxPuAR8Ah4BDwCHgGPQN4h4AlM3pnMK+wR8Ah4BDwCHgGPgCcw3gc8Ah4Bj4BHwCPgEcg7BDyByTuTeYU9Ah4Bj4BHwCPgEfAExvuAR8Aj4BHwCHgEPAJ5h4AnMHlnMq+wR8Aj4BHwCHgEPAKewBD6wPz582HUqFGZEgcOHAgjRoxIXMPy5cth0qRJMHv2bGjUqJFWeR999BH069cPKisrs+Tj6oT1ox4mOmgpKghhHXjFwQx1e\/HFF9mzctvlNot2mjBhAvTq1auaqlVVVTBy5EgoKirKuo86zpo1i8nHxdIUl6TyqPMFF1wAnTp1SlqUfz4FBORxhaswd+5cEpvGGV8oYcD+Onz4cJg8eTIrNujfrVu3pqxSq6ywMUDrYS+UEwQ8gSGCGQeZefPmZb3gxRfy1q1bYcCAAewFq3qRyLJxBhhxUOCdn5eL9ZuQBHxu7NixDKni4mKl\/nEhjUtguH5jxoxhVYs48zb37t2bEZGwwVIeIPlLQyQ4+Ft5eTlMnDgRcHDDeni5cduci+dEfHQJcC708nXoISD6XUFBAXsIx4Q+ffoABYmJM76oNDcZ73QIjEhsckVmPIFRWTn9+57AENkg6OWLHXP8+PEwderUai\/WqGpNOn9YOUEEBmXDfo\/SBwfQiooKKCwszLzA+UBKBB8rJi6BQf3wCppJkcuVXwZYJ7ZLfBbxHzJkCGzfvj2LoMj6xdWXEjPdslQY6Zbj5XKPQBCBQS3CfjfV0FUCIxKVOOOWKQ5ePv8Q8ASGyGZBMzC8aM7ky8rK2E\/8q4l\/RXE5\/Nrv0aMHW7oQZflLmC\/fiFPK3bt3ZzMCMqEI6\/BB5ChqWYTrjjMvRx11VGZ6N+oriNexcuXKrPbysurXrw9Lly5ly1tnnHFGZtaKE4IrrrgiQ\/xwxkDUQZ69wnvjxo2Dvn37QphOItHQISEo07FjR3j66aezlpCCZmDCZtRE2zZt2hTmzJnD9OMvixYtWmRsLM7yyL4i47Nz507A\/9A\/xHvyc\/LyFtZbWloa6CtEXcAXYwmBMKIiEo8tW7Zk9U1VPw\/yybDlYXF8kMcb3TEM+624tCvWr5qBuf322+HOO+8EHE\/E53TGQW4SeYwTl8OD2oAfNOIMzJlnnsnwbdeuHRu\/8cqXJWRLbulEsZ7AEJlBfoHIHV0eUPhX\/ujRo7NebDiI4CUug4gD1Ycffpi1HyZsFiCMwMjToqplEXEWCclE0IyFCKFcPuo+dOhQ9gJv1qwZI2fr1q1jpAVJF\/6NgxISAd6WwYMHs9\/5cpWsg1hf1D2U44PmlClT2NKXrD+fXeJLaij\/6KOPwq233sqWzeQ9MHywEwdS2YWCbMvJAw7COPXPSYuIDxIcUR8ZSz4IcwKMfyMJRAK7cOFCNkuG7QgifLJORG7vi8kBAmEEJmo5NGgZWiSw4rgRNQMjE1+x\/5iMYVHL4iL5Qjh19sPIOkfNhoqyfMzhYxDWh7OtQeMwl8UxAAkM7inEcV3cZ8fHlRy4ga8iAAFPYCy4hfilwVm6allInMGJIjAPPvgg01i1h0WHwPDZHvElLQ+W8gte9SUv1yu+hIPqC5sdEeuVdRBNFqVP0J6fKAIjzuZwsiViIxIG1CFoky\/+HrXXSB54o6bG5f07Yv04uIrPvvnmm3DPPfdkZnqCSJXuHiwLXcIXmQABCgIT9JHBPxxUBEZ3r41qDBMJVNjMri6Bkftx1IeMTG6i2iu2QSYwnFjhh4ZqPE9gbv+oAQKewBiAZSoqfv23adMma1ZFnrHp3Lkz7Nixg81MhBGYe++9l+3Ul2cFgvTSWUKSdcJygmZk+FKQWA\/OAvDn+X2cVWjZsmXWDFFcAiMOSEjawk7RhBGYMBIRtYQknmSSZz+CBqyogVBeRuNEVkVg5KWnxo0bZ\/bhyLrLX8DilLo8QxS1DGfq114+twhQLCHJyzfoV3wzv+iT2Nf4KTs+i8xnDXmr+QygyRgWdpIKxwy+PBN1CkkkNvzjgi+zc72CZkWDNuIGzcjwssRx2BOY3Pp5nNo8gYmDmvRM2PS8+NKTyYL8IhP\/tjkDE7SkEzYDgwNX0L6JqOlaqhkYPvB069YNFi9eDHjCKOgETRCBkWcuRHPJsznil5y4Ti4+g+Tjuuuuq3aKLIrAiM+LRBZ\/x3r4fgMRL3nWJ2gJiS8ZyTMw8v4febbGfzESdPSUiggjMKKNN2zYELkHRkXcRZ+MaqY4QyEvZ6vGsLA9WKo9MDKxQV9XLWWLbYiagYlqgycwKTm8QbWewBiAFSWKnUTcGMZnM\/jR6iBSErY3JIrAyB1OZ3o56hh11B6YsEEi6sWtswdGXpbB9op7YPjyGP9qC9uojM\/JU8e8fj49LttM9xh10Jdb0BJSUD0yiRMJLq73qwiMOLUv7pfB5\/BLkW8IFvWZMWMGa6qIHd8Tg7\/7PTBEHT2FYoL6OJ+p47MhMmmX78tkR9zPoVpS4aED8IUufjDgB47uGCYTaJHUi4cDxJmWsH+Lm+HFgw1yGAtuqqg9MDieqvbo8T0wfgkpBedXVOkJDKFN5GlS+cXLv\/Dl00b4wsJNZDh9i18b\/AsDp3KDlmXEesSTKEFf\/XIgu6DAbeLMg7hnR9zcJpYdNcPBX5a43yLsFJIugeGD3M033xx6RFo+hSROlYs6i7YQ8TMJZKc66SPWJy4F4e\/ipt0wAsMHZiQteKEt+MUJHtoTNyAitlGnkGS\/UO1dIuwGvihiBIKWX4KWS0Q59CE8rSZuhOeBLbEvtG\/fHt566y22ARx9KWoGRhwfxHrF\/qAaw3gMJjG4Ju8TOjMw4rIRJ21RS6ayCXgbUE8cm3DTu3iQAD8M5DaIM6LiMpffA0Ps4AmK8wQmAXj+UbsI6M4a4ECGV1gcGLta5q70qKU7lRYmU+6qsvx9j0A+I+DJfD5bL1t3T2Bqji1rXEt0BxokOnjkOWyfTE0BJi6BOVDwqSl29u2gRSAoNEHYEjNtzb402wh4AmMbYV9+LATk\/R6qQsQTRCrZfL0fl8Dgcz4XUr5a3eudFAH5RGDUnrqkdfnnc4uAJzC5xdvX5hHwCHgEPAIeAY8AAQKewBCA6IvwCHgEPAIeAY+ARyC3CHgCo4H3CSecoCHlRTwC+Y3A6tWr87sBhtr7fm0ImBfPWwRqat\/2BEbDJXGgo3AAqnI0VNYWcVEnVN7rpW1CMqxcxVwfCTNJ2+21Xb5Za8OlvZ5USPpxiw5JvZI8gdHAiaqDr1mzBlq1aqVRY+5EXNQJW+\/10vcBKqyo\/Fxf83QlbbfXdvlU6Hk9qZD0BIYOSb2SPIHRwImqg1O9aDRU1hZxUSdPYLTNxwSpbEjl52bapydtu71UdrGNkNeTDmFXsbTt63QImpXkCYwGXlTGd9G5XdSJ8qWsYV4jERfxotKJys+NAE1R2HZ7qexiGyKvJx3CrmJp29fpEDQryRMYDbyojO+ic7uokycwGk4piFDZkMrPzbRPT9p2e6nsYhshrycdwq5iadvX6RA0K8kTGA28qIzvonO7qJMnMBpO6QmMGUgB0lT9OkwRV\/uWrK\/XM7ErZQpwFUvbvk6HoFlJnsBo4MWNLzsn\/o0X35irur9s2TLo3LlzlrObPK8qP879DRs2ZHSK87wt\/RErTKYmYivXlcbfvM4g26ehD+LDdUlaf9euXUlO22l0KSdEbA\/qrr7MPIGx536u2ty2r9tDNLpkT2A0kKcyvovO7aJOfgZGwyn9DIwZSH4GJhQvV8eAfCRarmJJ9Q5L3OmIC\/AERgNQKuO76Nwu6uQJjIZTegJjBpInMJ7AJPYYdQGujqdU7zA1ArmV8ARGA28q47vo3C7q5AmMhlN6AmMGkicwnsAk9hh1Aa6Op1TvMDUCuZXwBEYDbyrju+jcLurkCYyGU3oCYwaSJzCewCT2GHUBro6nVO8wNQK5lfAERgNvKuO76Nwu6uQJjIZTegJjBpInMJ7AJPYYdQGujqdU7zA1ArmV8ARGA28q47vo3C7q5AmMhlN6AmMGkicwnsAk9hh1Aa6Op1TvMDUCuZXwBEYDbyrju+jcLurkCYyGU3oCYwaSJzCewCT2GHUBro6nVO8wNQK5lfAERgNvKuO76Nwu6uQJjIZTegJjBpInMJ7AJPYYdQGujqdU7zA1ArmV8ARGA28fyG5\/Bu1cBrrzgeyqB0lEGwQF9jtQAtktX74c+vTpk+mxc+fOhU6dOlXrwVVVVTBy5EgoKyuDpk2bwpw5c6B169bV5GwP6q6+zGQgvJ4aLwFNEVextO3rmvCQi3kCowEplfFddG4XdfIzMBpOeYDNwGzduhXGjh0LY8aMgUaNGgGSmUmTJsHs2bPZ3+KFv+M1YsSISDmqfh1mLVf7licwZv3LRNpVm9v2dROMKGU9gdFAk8r4Ljq3izp5AqPhlAcYgZERQUIzYMAARlLEWRg++1JUVAS9evWCjz76CIYPHw6TJ0+uNgtD1a89gTHz17jSro5VYntc1dG2r8e1adLnPIHRQJDK+C46t4s6eQKj4ZQHOIGJIibz58+H8vJymDhxIqxcuTJ0poaqX3sCY+avcaVdHas8gYlr0eTPeQKjgSHVQOdiB3RRJ09gNJzyACYw8ixLEFpIYkaNGgXdu3dnRKagoKCaGPZrvJYsWWIGuKY0JkrFhKSuX15POgu5hiUmaOXX6tWr6RrqSEmewGgYwhMYDZCIRTyx0geUCisqP9fX3FySkxfcnIvLR\/Lll5DMMaXyH\/OazZ7IBz1d1TEf+raZN+yX9gRGAzUq47vo3C7q5GdgNJzyAJyB0Zl5kZeWop6h6tdh1nK1b8n6ej3N+luUtKtY2vZ1OgTNSvIERgMvKuO76Nwu6uQJjIZTHmAERoe8ICRBMzD9+vWDKVOmVDtyTdWvPYEx89e40q6OVWJ7XNXRtq\/HtWnS5zyB0UCQyvguOreLOnkCo+GUBxiBwZkVJCKVlZVZ4GAsmDZt2sCQIUNg9OjR7KQRP6GEG3jxmjBhAjuRJF9U\/doTGDN\/jSvt6ljlCUxciyZ\/zhMYDQypBjoXO6CLOnkCo+GUBxiBMUNET5qqX3sCo4d3UilXxypPYJJaNv7znsBoYEc10LnYAV3UyRMYDaf0BMYMpABpqn7tCUxiU2gV4OpY5QmMlvmsCHkCowGrTyXgUwng4MkvDOfP\/w4K7c8JGP4\/F\/dl3eLWj0cua+JRy7Au7gnMfmTygRjki56uYmnb1zVeo1ZEPIHRgJXK+C46t4s6uTxYuYgXlU5Ufq7RpZwQsd1eKrvYBsvrSYewq1ja9nU6BM1K8gRGAy8q47vo3C7q5AmMhlMKIlQ2pPJzM+3Tk7bdXiq72EbI60mHsKtY2vZ1OgTNSvIERgMvKuO76Nwu6uQJjIZTegJjBlKANFW\/DlPE1b4l6+v1TOxKmQJcxdK2r9MhaFaSJzAaeDXr8ivoOXB0RrJPx+PgvJMaaDyZLeKic7uokycwZq5FZcOaOsiFoWm7vVR2MfMGc2mvpzlm+UZabfs6HYJmJXkCo8Br0uIKmLR4DYzo1gpGdCtk0vjbuq1VMLO4rRHaLg4ULurkCYyRW5Ftwqypg5wnMNH+5OoYkI8zRa5iWVP7ticwEX375Y+3w0ufbId1ZTdAh9OPguKOxzHpw5tfCVPfaALnn9jAaCbGRed2USdPYDyBMUMgnrTtQd3VvpWPxMDlMUHE01Wb2\/b1eD0w+VOewERgOKh0FUw7\/1343wnD4KCeL7CZmLIbO8CZ3z4C3+3+FIa\/399oFsZF53ZRJ5cHKxfxotKppg5yfgbGz8Akf1XqlUDVF\/Vq05eqqX3bE5gIHxj3+Dy49exN0O7iGSw+BhKa0tc2ZUjMb59rALNuukbbi1x0bhd18gRG26WYIJUNa+og5wmMJzBmPSq+NFVfjK9B8JM1tW97AhPhKauW3Axtu94DaPx3nx0MVR\/cDb94awRsPPgMmNm7LTRa+zt2X\/dy0bld1InypaxrG105F\/Gi0qmmDnKewHgCo9u\/k8pR9cWkesjP19S+7QlMhKd8tWIY3Pv5IJh1w3+zGRj8++v1T2RIzOPtJ0O7S5\/S9jUXndtFnTyB0XYpPwNjBlWWtO1B3dW+JUPm9UzgRNKjrmJp29fpEDQryRMYBYHBfS5PPfYAVC76Xya58+ViqPj0Y5hR0RN+e+obnsCY+Zu2tKsDgYt6UelUUwc5PwPjZ2C0B56EglR9MaEa1R6vqX3bE5gIT9n3xXLYu2U5nPyrxTDoD\/ezo9N4TWzxe0Zihr8\/AA496hxYOKiDlr+56Nwu6uRnYLTcKSNEZcOaOsh5AuMJjFmPii9N1RfjaxD8ZE3t257AKDwFl43uefAJGHnf\/mR+3+3ewJaRvl7\/F1i\/dQ9csHxyVoyYqOJcdG4XdfIExmz4orKh64Pc8uXLoU+fPhlw5s6dC506dQoEa9KkSTBr1ix2b8KECdCrV6+cf5VS2cXMG8ylvZ7mmIU94SqWrvftuBbwBEaB3LobWwEcvAka9r4evvq+CurXqw\/7PtsDe95+Hw45bRe89nkD+MVbw9nJpKaHbMtkH+YvYfw\/z0i8bNky6Ny5c9aXs3hfdn782\/b9DRs2ZHRKo\/6w9iFWzZo1y0k25yBbBWHPDeezUccdbuI\/t3XrVhg7diyMGTMGGjVqBEhmkKTMnj2b\/S1e8+fPh\/Lycpg4cSJUVVXBkCFDYPTo0dC6dessOduDuqsvM9kKXs\/4fpkvWNr2dToEzUryBEaB1+aZ\/eDjf\/wZCjucD8eVvJCR3rZgLBz5swGw\/Z+d4dXtp2RITFSKARcHChd14oSCEz8zl7Yr7SJeVDrl0yCHhGbAgAEwYsSIrFkYJCzjxo2Dvn37ViMssmfYbi+VXex6NN0xfK+nu1ja9nXbtg8r3xMYBfJIVC753Wx4\/LT1UNCuS4bE4O8Nrx4DuE9mR3kxPLnpPLYnBmdiwkiMiwOaizp5AmM2HFDZMJ8GuY8++giGDx8OkydPziIqfKYGZ1umTZvGgPRLSNH+ROU\/Zl5rLp0PerqqYz71bRPPSJ3A4BfTyJEjoaysrJre3bt3Z9PABQUFJm0ilX1m8Hp45Mn7YVCfnrBz6SNQq3Eh1O9yDXzx5odweekPy0FYaY+ZK2Ddtj3w1u3nBurgonO7qJMnMGYuTGXDfBnk+JhRVFRUbW8Ln5nBvTE4OxNGdBBh2+2lsouZN5hLez3NMQt7wlUsbfs6HYJmJaVGYHBg6devH9N2zpw5gdO9fNNe06ZNQ2XMmmsuveyuz6Dv\/eeyODB73l0KlSUXMgLz\/u6J0HnYsdUKjCIxLjq3izp5AmPmp1Q2zIdBjpMXHBOQoMgXEhhxz0sU2cH24rVkyRIzwDWlcX8Z7uNy\/fJ60lnINSy7du2aaRy+w2ralQqBwUHmwQcfhMGDB2vNruAgNGPGjMABy7ZBOIH564Q3AGdjrrxzPRz6j0vhnbXXQ49\/3V+t+nVb90D7O19hiR7l49VULxrKNruokycwZhamsqHrBCaKjHDEuExxcTHbG6MiMDYHdSq7mHmDubTX0xwzPwNDh1mSklIhMEkUzvWzLwx\/FV4s\/ynccstv4e9\/OATWfNYDLil+DNaV74ILJ5\/DZmPES4zWK8eIcXGgcFEnT2DMvJzKhi4TGB3ywlHDU0gVFRXsgwdncYcOHRo4g2u7vVR2MfMGc2mvpzlmnsDQYZakpNQJDF+zXrlyZWg7XFlCQgX\/fPknsK78K2hRVBeKap0EjQfNqUZiMFovBsDDvEldzvkpjOhWyNrm4kDhok6uYuWqXlQ2tP1CTzJQ8SXnysrKrGIwFkybNm2qHZUW48CExYux3V4quyTBTedZr6cOSnoyrmJp29f10KGXSp3AYJNwsCksLMzakCd+RfG4DnfffTc9AooSxT0wXJSTmF9PeRW+WfibUBKD0XoxRsz9\/S9iJ5NcdG4XdXKVKLiqF5UNa+ogF9bFbbeXyi62Bz2vJx3CrmJp29fpEDQrKXUCI2+64+rjF9f48eNh6tSpsGXLFvZv3Oyb64ufQrrllluyql425TM4svlh0GfEfNj2fyXQtOQFqN2uS5YMz5uE0XqDAt3lui1B9bna4bxe+t5BhVVNHeQ8gYn2JSr\/0ffYeJL5oKerOtbUvp06geEzMBj2m0\/18tNHAwcOZOvYac7AoH5hxseZGLy6\/fhOdsQ6jMS89Ml2NhPzwBVN4Krz28brvZaecrXDeb30DU6FVU0d5DyB8QRGvzclk6Tqi8m0qP50Te3bThAYhFtc4xb3vIgzMXLIcGojmw50X67\/Bmae\/T7bDxNGYjB3Eo\/WO2z1aHinJDt2TK7aEFaPqx3O66XvGVRY1dRBzrRf6yPviQEVVjrlUPm5Tl1xZVzVsab2bWcITFyHycVzUcY3JTH37xyvnb06F21ztcN5vfStT4VVTR3kPIHxREu\/NyWTpOqLybTwMzDU+OV1eSYD+8aSC6Hq3aXVlpPElAML9w5xhsS42uG8XvpdhgorEz\/X185dSdvtpbKLbQS9nnQIu4qlbV+nQ9CsJD8Do4GXqfHDSMz6d\/4KdVf\/D4x4fwAc1vxKmFmc\/n4YVzuc10vDMf8jQoWVqZ\/ra+impO32UtnFNnpeTzqEXcXStq\/TIWhWkhMERsyHhPmP+vfvD9OnT2cnkNLa9yLCGMf4QSSGO\/fLH2+H7vetgBHdWmVixJiZjU7a1Q7n9dK3MRVWcfxcX0v3JG23l8outpHzetIh7CqWtn2dDkGzklInMGKEzZYtW0JpaSlL4Lhw4UIoLy9PPZkjwsmNLzsn\/o1Xq1at2P\/x7w1P1AE8Yv2Lp06AbxdeBLB1AzS95TF2xHrZsmXQufP+Tbylr22CQaWrYOCPGsCO72uz33bt2gW\/ubB1Jpt1UPm8Ll6fXL\/pfczdwXXSaZ9p+XH1Q6wwj4yIrVxWGn\/zOoNsk4Y+iA\/XJWn9mDfFZmh9s6HJvrTtQd3Vl5mMrNeTztdcxdK2r9MhaFZS6gRGjAOD8V44gcEXK48Dk\/YsTOuLr4a6wx+EhzrUg7Mb1FIijLFj\/j1\/GyMxh\/2jO+z9vAKOGTQHNtZpmXkhYyGYM6njYUtg+kVfZsqc\/t5ZUHnwGTlbXnK1w3m9lG6WEaDCqqYOcmFI2m4vlV30PSGepNczHm5BT7mKpW1fp0PQrKTUCQyqi5F4MUR4z5494emnn2ZLSIMGDQJcTgrKOGvWxOTSaPzhz78Hf1yzR5vE8Gi9Ion59soJ0Kprb6YQLiMtWfoQNKv9BQx9sxsLdIfReqs+mA6Pv\/g6tP3vezMzMclbEF6Cqx3O66VvdSqsauog5wlMtC9R+Y++x8aTzAc9XdWxpvZtJwgMujMPXsdde8KECVmpBeK5PM1T3PjXrtgFr2\/fZ0Rivly\/Fwa9fgqwPTGVH8MJD6xnSo17fB7cevYm2PfFq4ApB3714Qx46\/Zz2T0kMQtWF8KvL76MpgERpbja4bxe+qanwqqmDnKewHgCo9+bkklS9cVkWlR\/uqb2bWcIDLXBKMtD4187\/9ewqKIMqo6cBsfW7wCLzj1CqwqciUESM\/DZ2rBuajHU2rEJWty3BjBrdd0Od7EyeMqB2zY\/kDlevfiJa6HbVQ9p1ZFEyNUO5\/XStyoVVjV1kPMExhMY\/d6UTJKqLybTwhMYavzyujx5YMeZmMo932mRGAx0VzZ4\/6zLeb+rhNoL72B7Yr668lJo2\/WeDC5IYjDlAAa6u7e4LewsL4Z2lz5lHTdXO5zXS9\/0VFjlgsCIJw7lFuKSMW7gLygo0G98Aknb7aWyS4Imaj3q9dSCSUvIVSxt+7oWOBaEUpmBwY27AwYMgJUrV0Y26YwzzoDZs2enfpRaNj6Sl0te2cE29OLGXtXFo\/Uec+YhcM2f6sK6G1tBnaKj4dFTX8scoxZTDmCgu8mnzIb655Wqik5839UO5\/XSNy0VVjYHOZ4qBFuFSVlbt25drYF8GVlMJaKPgrmkzfaiNlR2MW+Z2RNeTzO8oqRdxdK2r9MhaFZSKgRGVhE38RYWFmbteeEJHHP5RRYGnWz8m5f+hu2FweUkUxKDeZN6P3AIVI5tCytbdoL3Tp+TITEVGz6GI968CD7dczRsafkH6HLOT82sGUPa1Q7n9dI3JhVWtgY5\/GB58MEHYfDgwVqzKzhLM2PGDOsb+G21l1uOyi76nhBP0usZD7egp1zF0rav0yFoVlLqBEY8Ri1+lbmQxJFDSWX8V\/\/yMSy5cTd0HnosnHtNFXzx546MrNRq1JNV9d3uT+E4+DscdmI9mFHREy66pMT6SSRXO5zXS78jU2FF5ef6mptJyhv9efb6qFLw4wivoNOMtttLZRczlMylvZ7mmIU94SqWtn2dDkGzklInMKguDjKzZs0CPiDxgWrgwIHWv8J04IoyPs7E4J4YnRgxsnNvntkPqlbNg9qnnQJ1f\/Rjpkqto8+Bb9Z8DHu3z2IpB6654iarJMbVDuf10vHM\/TJUWNkc5OS9L9jX8erTpw\/7v2r\/C37ojB07FsaMGcOWlHGMwHEjaolZNY7YbC+lXfQ9IZ4klf\/Eq13\/qXzQ01Udbfu6vhVpJZ0gMNgkvkaO8WBytQauC2WY8Td9VQlXLeoO39TpC9\/U+bWSxMjOvW3BWCg49QKoLLkQCtp1geNKXmAq4e8Nrx4DPWauYBt7eYwYXX1N5FztcF7yXoo7AAAgAElEQVQvfStSYWVzkBNnQmRiwckN9nvduE98Hx3Kd+rUqRpYnPDgjbBybbbXExh9\/9WVpPJz3friyLmqo21fj4MVxTPOEBiKxtgqQ8f4OjFigggMEhW8Vl99EDT8eQkjLpzA4O9IYtZt25OJEUPdRlc7nNdL39JUWOn4ub5WP0jKy8Ri+pBevXoxQdMlY5QfPnw4TJ48OXBDMN9XV1FRwcr3S0jhlqPynzi+YfJMPujpqo62+raJ\/WzIpkJgXN3UFwawyvgrNr8BuLH3to4lcHFh91A7yc6NKQeObH4Yk9+3uQJ2Ln0E6ne5Brav\/wYuL92fM4mTGPz\/wkEdyH3A1Q7n9dI3NRVWKj\/X1yhbMozAFBcXZ2ZPTAhMEAESa8SyHn30UbjtttvYZuAoAoP3lixZErdpkc9hOhTM5+X65fWks5BrWGJ+M37VxDxnqRAYBNTFY5VxCQw+x0nMPV0egA6NzwosSn7RLLvrM1hbvosFurt0RjOWN6nq3aWwttVyuHDyOZky1m3dw\/ImnX9iA3ISQ\/XyoxsC9pfk9dJHlAqrfCAwquUmvD9u3Djo27cvm5nxm3jVfkTlP+qakknkg56u6mirbyezaPKnUyMwXHWXAlslITD82agYMbJzvzD8VWjf8e\/wl9ubw9f1i6DPiPmw7f9K4J2110OPf92fpQ6P1vvHb\/9EmujR1Q7n9dLv3FRY2RrkqGZgVDMv4ocR7qUTr6ADAbbay+ulsou+J8ST9HrGwy3oKVextO3rdAialZQ6gTFTNx1pU+Pzk0k9mhwGv29bJ6N00AxM52HHsvtiyoFFP5sEnYc1YfthxIuTmEVHLsjEjkmKiKsdzuulb1kqrEz9XFdDisCVOuQlSB8\/A6O2EpX\/qGtKJpEPerqqo62+ncyiyZ\/2BEYDQ1Pji4HukMAgkcErisCIKQdaFtWD4\/\/VBJqWvAC123XJaIjRejGHEp5MOqKolOR4tasdzuul4Zj\/EaHCytTP9TVMLimeUhRLw+PYbdq0gSFDhsDo0aOrbej1BEaNPZX\/qGtKJpEPerqqo8t9O4lXeAKjgV5c4y\/c9A3csWp35nh1FIFBNXjKAfz3Dfc+xZaTMPFjrcaFWSRm+z87w6vbTyEhMa52OK+XhmMeQARGHw19ybj9WrcGV31Y1t\/rqWtRtZyrWNr2dTUydiQ8gdHANa7xcSZm+Z6usPfwbozEHLVtPbRq1SpTo3gKif+IJObf87fB6b0awlmNfs1+5vFhuMy+L5bDjvJieHLTedD2v+9NNBPjaofzemk4Zp4QGIolJH009CXj9mvdGlz1YU9gdC1oLueqzW37ujlSNE94AqOBY1zjY6C7ca+VQPk3v4BvDz0D\/tBkJ3Rv21xZI55Mevzy1XBu3ypouaZTJj6M+CAnMX\/d0Qv6\/3KisswwAVc7nNdL36RUWMX1c31N90u6kvvMdnup7GKKr6m819MUsXB5V7G07et0CJqV5AmMBl5Uxo\/j3HveXcoi9cr7YVBtTmIwe\/U1V96k0ZLqInF0ilWR4UNeL33AqLCi8vMozV3KfWa7vVR20feEeJJez3i4BT3lKpa2fZ0OQbOSUiMwSXOjmDUzmTQ3vuyc+DdefFko6j4er75l+UY48sgj2XISXrrPN3j9T2w\/zEGj\/h+0OvOHAHf4\/GHbFkPBhklw\/84J0KfonKwlKp3yMfBS5877y0zSPurnly1bxoKAidjKWKfxN68zCNs09EF8uC5J68egV7kIduVK7jPbg7qrLzN5NPR6Jns\/iE+7iqVtX6dD0Kyk1AgMdW4Us2abSVMYHwPd3bTsDviq4eNwdoNaGRKjq8nGkguZqLwfBn9LEujO1Q7n9dL1DLqgfxR+rqu1C7nPbLfXVR\/2BEbXS83lXLW5bV83R4rmiVQIjI3cKDRwBJdCZXx07i0Nm7Ps1de3qg3XF9bWVhtTDay7sVXgfhgs5OWPt0P3+1bAiG6tjGLEuNrhvF7arkEWtZjKz\/U1T1fSdntd9WFPYOz5nas2t+3r9hCNLtkpAhM3N4pt8KiMj879\/kHvsI29USkHgtqDge62vPkRnHPy76D9\/96VFR+Gy3MSY5K92tUO5\/XS92oqrKj8XF\/zdCVtt5fKLrZR8nrSIewqlrZ9nQ5Bs5I8gdHAi8r43Lkffu8BePjdWbFJTPcfXRy4qRebMmlxBUxavAaQxGAW65c+3pZpYZ+Ox1U7cu1qh\/N6aTjmf0SosKLyc1lzV5O32movbz+VXfQ9IZ6k1zMebkFPuYqlbV+nQ9CsJE9gNPCiMr7s3BjkDoPd4aZe3Bejc+FMzJ73XoSfFT8WuB8Gy8CUAxit97VvusKQsz7LFDv9vbOg8uAzsnIpudrhvF463rBfhgorKj8P0tzF5K0220tpF31PiCdJ5T\/xatd\/Kh\/0dFVH276ub0VaydQIzIABA2DlypWRrTnjjDNg9uzZ0KhRI9pWG5ZGZfwg58b9MHhCCUlM09oHKzXj0XqPafA69By8DBoPmlPtmf\/39jvQcFU\/9vtpP38tc7\/qg+nw+IuvZwW\/c7XDeb2UrpARoMKKys+jNHcpeavt9lLZRd8T4kl6PePhFvSUq1ja9nU6BM1KSoXAmKmYvjSV8YOcG6P1YqC7Y+t3gEXnHqHVWJHE\/PLJE6rthxn3+Dy49exNUPXB3bDx4PbQ7tKnskjMgtWF8OuLLyP9etdS3EDI1YHARb2odKLycwMzpypqu71UdrENkteTDmFXsbTt63QImpXkNIF5\/vnnoX379jV6Boab65JXdrAZGB4jRmVGHq231bEL4fLS87NIDCZ8rNvhrkygu8\/qXgptu96TKXLxE9dCt6se8gRGBXLAfRcHKCqdauogF2Zm2+2lsksMNzV6xOtpBFeksKtY2vZ1OgTNSkqNwMj5UTCrbKdOnZj2fL28cePGNX4JiZvryr9fC58cPsUoRgwnMWdfsgp+MqdXxvKrltycISw8Wm\/BybdAwcm\/ZTFjdpYXZ2ZlXO1wXi\/9jkyFVU0d5DyBifYlKv\/R99h4kvmgp6s61tS+nQqB4evgTZs2hREjRjDCMnz4cCgpKYGHH34YysrKYODAgeyeCxeV8aOcG\/MmXf7c7VB15DSjGDGcxPS8ZRm0Gz2IwYUzMPd+PigTD4aTmCOKSmH+u7WgR62pUP+8Uj8DE8O5XBygqHSi8vMYsKbyiO32UtnFNjheTzqEXcXStq\/TIWhWUioEJigfCg8v3r17d5g4cSIUFBSYtcSiNJXxbTm3nC8JCcudc+exGZcR3QoZMh+v\/DM0Wvs7lsH6mituglpH75\/tsqVTUnN4vfQRpMKKys\/1NTeTXL58OfTp0yfzkDhrK5akK2e7vVR2MUPJXNrraY5Z2BOuYmnb1+kQNCvJGQIzf\/58KC8vd468IJxUxtd17j9W7IE\/rtljdLx624KxLF8ST\/qIszCvVNaFN2pdA+u2VsF3uz+FCS3Hwne7NwDOxHgCY9ZRuLSuDeOVHu8pKp2o\/DxeK6Kfwo+esWPHwpgxY9ieOCQp+NEjn1LUlaPs1\/n2MpP1pfIfG3YXy8wHPV3V0eW+ncRvnCIw2JBevX7Yy5GkYZTPUhlf17kxb9Jv3vgY9h7ezYjEYL6kvZ9XQIv79ieZxJmYvVuWM\/KCV62jz4GrF57I\/r1wUAf2f12dKPHUKcvrpYPSfhkqrKj8XKW5eJQaZ1z79+8P06dPh6lTp2pv2Od76HCZme+dC6o3Ss52e6nsosIz6X2vZ1IEf3jeVSxt+zodgmYleQKjgReV8U2dm8eI+X3bOlqB7ni+pIJ2XUKD3PHEjzxnkqlOGnCRiHi99GGkworKz6M05+SlqKgIWrZsCaWlpWzWdeHChUYzsHzf3OTJk6F169ahVUbJ2W4vlV30PSGepNczHm5BT7mKpW1fp0PQrKTUCIwPZKc21PjXSmBRRZlRygHcD\/PnK1bD1\/WLYNDrpwRWIuZManrINmjVqpVamRxLuDoQuKgXlU65GOTE\/W9btmzJEJgNGzbA+PHjtWZhRBIUNWOrkrPdXiq72O56Xk86hF3F0rav0yFoVlIqBMZMxfSlqYwfx7kx0B0uKZkkf9zw2IPwl9ubQ63GhXB6r4ZZAGIQvEtnNM\/kTHrgiiZw1flt0wdZ0iAOVrlohIt6UelE5ecqO+DelcrKSujZsyc8\/fTTbAlp0KBBgMtJqpOH8gnGsLp05LC9eC1ZskSlcqz7SMqaNWsW69lcPuT1pEPbNSy7du2aadzq1avpGupISZ7AaBiCamCP+6LBVAMY6A7zJekGulvw48fgo7fPgKatvwA8Ys2v1xadkokZc8Psf8K4xtfBoUd1yhyr1oAjJyJxsbKtnIt6UelE5ec6NpBPCk2YMEG5\/001o8Lr1ZWz3V4qu+jgmUTG65kEvexnXcXStq\/TIWhWkicwGnhRGT+Jc7++fR\/gnpgeTQ4D3BOjul4Y\/iocveoKKPvXs9CiqC788qn9m3eX3fUZdB52bObxfjPLYGrzwc6RmCRYqbBJct9Fvah0ovLzJPiGPatLSnTlsB7b7aWyiw08xTK9nnQIu4qlbV+nQ9CspFQJDB6dnjdvXtZRSPwyGzp0KMyZMydyc55ZM5NJUxk\/qXPf+PJ9UP5NH0ZgkMhEXUhUOl6yCt76n2Hw\/MrZbCkJl45kAvPyyo\/gzqeXwOPtJ8Hhza9iKQhcuJJiZasNLupFpROVn9vAnkfnxqUn8cJYMG3atIEhQ4bA6NGj2a1+\/fqxJSpZTj6tZLu9VHaxgacnMHZQddXmtn3dDprqUlMjMEHkhavLB6ubb75ZOa2sbmJyCSrjUzj3wk3fwB2rdiuPV3OigvFhPnjwGUZikMB8ue6brBkY1Kny24Ys8B2SGCQwSGTSviiwstEGF\/Wi0onKz2Xc5bQhYXbJdfZ5W+3l7aOyiw0\/9gTGDqqu2ty2r9tBU11qKgQmKBKvrGpYoCp1k+glqIxP4dyYcqD7smXKGDHiTAvGh3nv+SPg1Q\/+kLWchEhxnSYtroClr\/6dkRgx0B09mnolUmClV5OZlIt6UelE5edmiKYnbbu9VHaxjZDXkw5hV7G07et0CJqV5CyBwVkY3WOVZk02l6YyPpVz46kk3A\/z7aFnhM7EPDN4PRzZ\/Idlpm0LSmDNZz3hqz1N4RdPnQAti+oxIESdesxcAWd++ygMLnw6dRJDhZW5taOfcFEvKp2o\/Jwac1vl2W4vlV1std\/PFNEj66rNbfs6PZJ6JaZCYFA1PEqJV9ixSZdSC1AZP5fOvXlmP6h1zP48SHjt+7wCdi59hP37hAXfZ36XdUIS0+PQqXBFk5dTJTG5xEqvq+yXclEvKp2o\/DwKT53lJEzymos9cLbbS2UXE\/+MI+v1jINa8DOuYmnb1+kQNCspNQIjDmTiEUokLqNGjYJcDWI6cFEZn9q5cSYG48ToxojhSR8b\/rwEGl49JvCFzCP1Pt5+MpzT4P3USAw1Vjp21pFxUS8qnaj8XIUjfrwUFhZm7W\/Dfl9RUcE+aPjHy913360qKtF92+2lskuiRmo87PXUAElTxFUsbfu6JjzkYqkRGN4SOR6ES8SF60hlfBvOPfv9p+GejRdox4iRM1cH6SRG6j3vpAbkTqdToA2sdOpVybioF5VOVH6umoHhJ4bEFADikjFG6MXlY5yFsXnZbi+VXWxi4OqsYlCb8wFPV3W07eu2fTSs\/NQJTFoNN6mXyvi2nNs0RgwuL+FyEmau3linZWAqAdzUO2nxGii7sQOkQWJsYWVi93wZRKmwovJzFcY4AzNr1izAI9B4tJl\/xAwcONDPwKjAs3Cfyn8sqJZVZD7o6aqOuerbtn1ALt8TGA3EqYxv07lNUw7gyaQv3vgQXt33D6hdu3Ym0J0IB+6HWbdtD7x1+7kaKNGK2MQqiaYu6kWlE5Wf6+ArxnURZ11zuXnfdnup7KKDZxIZr2cS9LKfdRVL275Oh6BZSZ7AaOBFZXzbzo0kZtPujXDPBbOgSd2mypYxErNmN5T9\/cFqx6v5w0hi8Fo4qIOyPEoB21jF1dVFvah0ovLzuNjm+jnb7aWyi21cvJ50CLuKpW1fp0PQrCRPYDTw4saXnRP\/xotnc1bdX7ZsGXTu3DlTo+nzqvLx\/uKDjoM\/rtnDjlcftW19tH5vLoPvJ\/x\/8PmXZ7NAd8eceQhc+2y7LP0qd+yD7o9ugCFnfQY31B+VlXLAhv4cS8QKE+GJ2MpYp\/E3rzOo7Wnog\/hwXZLWj4nfamLCt7AubntQd\/VlJuPh9dR4CWiKuIqlbV\/XhIdcLBUCw08g4fq3KvsseYtjFEhl\/Fw4Nwa6u\/3DIwD3xSCJwQSQUdeaJfPg+z8Ww47jh8Ki0l9B56HHZkXqxWf5pt4ne30Np392Q7VIvYc3vxJqHd0pBrLhj+QCqzgKu6gXlU5Ufq7CVd64z+V9JF4VcnbuU\/mPHe1+KDUf9HRVx1z1bds+IJefCoHhSogDWa4HLxOgqYyfK+dGEvOzlz+CY+t3gEXnHhFNYNasgQav\/wm2\/V8J7P3JM\/CX25uzlAOYO0m8+KbeV65+FRpvnpWVcqDqg+nw3e5PSfMo5QorEz9AWRf1otKJys+jMOUfL7179049TYjt9lLZxdRHTeW9nqaI5d+Hl21fp0PQrKRUCYyoqhzgip9QMGuOHWkq4+d6oMBovXjhTEzYxXXC\/TBV7y6FL9o+Cf+cdWJWtF7+7LjH58HeLcth1OXnw1crhmXFiEESc+hRnchmYnKNla7nuKgXlU5Ufq4iMEHHqHXxp5Sz3V4qu1C2Oagsrycdwq5iadvX6RA0K8kZAiOrjUctcYZm9uzZ0KhRI7NWEUtTGT\/Xzv3mlg3Q\/9\/7l5HCSIyoE5KYvZ9XwNpWy2HZlM+qkRgkLcXlfRi6lx05Hy47Yj48c+j\/womtu7Cj1nifKpt1rrHSdRkX9aLSicrPVViKQetUsjbv224vlV1sYoBlez3pEHYVS9u+ToegWUnOEhizZtiVpjJ+Gs5duec7uOSVHXB9q9pwfWHtakCJOu3bXAGVYy6EQ48phCUrH4J15V9lkRgkKPd+PojFhxnRrRXcdMxM+Hr9E3D\/zglQefAZMLHF76H+eaUkxkgDKx3FXdSLSicqP1fNwAwYMABWrlxZTSzXy8i220tlFx2\/TCLj9UyCXvazrmJp29fpEDQryRMYDbyojO+ic8s68Ui9Be26wHElL2ShwwnM+Sc2gO73rWBB7nBT77dVG+CvX\/aCK5u87AmMhj9Ri1D5FZWfU7fPtDx5k3DYcrTt9lLZxbT9pvJeT1PEwuVdxdK2r9MhaFaSJzAaeFEZP03nfvi9B+Dhd2dVy5sUpBMnMY0HzYH6Xa7JIPSnZ\/8KV59QAQUn\/xbESL2cxCzcOwR+ffFlGoiqRdLEKko7F\/Wi0onKz9XWtSeBe+nGjh0LY8aMYUvPSGZwOTpoKdp2e6nsYg+t\/SV7PekQdhVL275Oh6BZSZ7AaOBFZfy0nRs39eKSkngyKUynbQvGspNJmG6gdrsuDKVBpatg8ikPw8F1jmckBoPc4ezLtPPehRaNakOvf15AFvAubazC3MJFvah0ovJzGTu+QX\/z5s0wZcoURihytYTE68ZwDRi2QbxstZfXQWUXjSEqkYjXMxF8WQ+7iqVtX6dD0KwkT2A08KIyvgvOLZ9MitKJn0ziJAYJzMzitrDvi+XsNBIenX7pk23wxOpC+N21Q+Gm0lWewGj4E7UIlV9R+Tl1+5KUh6kJhg8fDpMnTwYxeSSWabu9VHZJ0n6dZ72eOijpybiKpW1f10OHXioVAoNfRfJRSvwqKywszMSGyGVOFBWsVMZ3xbl5yoEFl5Qpp4\/5yaQjR30AM89+PzDlAM+ZNLN3W7LEj65gJfuGi3pR6UTl56r+lKv7VVVVMHLkSCgqKgqMOWO7vVR2sY2X15MOYVextO3rdAialeQJjAZeVMZ30blVOuHJpHU3tgLc1PvNT8rg8ctXw5HND8sEusNkj6X\/2ghHfHUI7PplLT8Do+FP1CIqG+rWR+XnqvrEfSnPPfccjBo1CsSEjqrnde5z8oLlhkX7xvbitWTJEp0ijWU2bNjA0mG4fnk96SzkGpaYHoRfNTFNiCcwGr5LNbBTvWg0VFaK\/LFiD8uZ9IcmO6F72+aR8nxTb8Ofl8Df5v6SHa8+6PLDYcMF37Hnzj+pITRbehD03Pge4AklTDmwo7wYCk6+he2ViXO5hJWov4t6UelE5edR9ubEori4GNq0aQN4pJoTjNLSUpg4cSIUFBTEcZnMM6qZFy5ou71UdkkEhsbDXk8NkDRFXMXStq9rwkMu5gmMBqRUxnfNuXE\/zJtbN8B5h0xjp5OiLk5iMMjdCZe1YzMxYsqBZXd9Bgdffjg7Xl3csQncfdGXjMRgYLvDm1+lgXK2iGtYce1c1ItKJyo\/jzK2uHy8ZcuWzAkh\/Pf48eNh6tSpiQJX6pIX1NF2e6nsYtx5DB\/wehoCFiHuKpa2fZ0OQbOSPIHRwIvK+C469y9f+QJq164dmW6AQ4Qnk5bdtQl+MqcXfPjOGfDM4PWZQHdIYDoPOzaT+BED3f321NerpRzQgJuJuIiVq3pRYUXl57ozMGvXroXy8nI267Jw4cLMv5PMwODeuX79+kFlZWWWGkGxYGy3l8ouun0mrpzXMy5y1Z9zFUvbvk6HoFlJnsBo4EVlfBed+9WP1sLADUeGRuqV4Vn0s0lwSp2RLD7MSy\/1h1XPHwk9b1kGn350NCM2ePHs1XhiqUetqSxa7xFFpUZ5klzEyhMYjc6iIcIDzfF9L\/hI2EkhjeJii1D16zAFXPVhWV+vZ2wXqvagq1ja9nU6BM1KSo3AhIUTF9XPdWjxMOiojO+ic6NOWxo2B1xOKvjyVrij\/eVwcWH3UC\/CmZaTvurNcia1uG8N\/PnyT9iemBZFdeGXT52Yea70tU0sbgyP1ovHrk1IjItYeQJjNri4Lk3Vrz2ByY2lXR0TxNa7qqNtX8+NB1SvJRUCk1Zj49ZLZXwXnZvrxDf1YtJHTP4YduGyUd3albBz6SNQq3Ehm4nBxI94DXr9FHZCiV9B0Xrrtb9LaybGRaw8gYnbg9x8jqpfewKTG\/u6OiZ4ApMb+wfVkiqBwdgvuFbNTx5gllo8UonXwIEDQ48\/5houqoHOxQ4o6hQUqVfGGvfBNLx6DIgnk\/Bv\/rssj7MwOBsj5k1q8ONlShO6iJUnMEqz5ZUAVb\/2BCY3Znd1TPAEJjf2d4rAIFnhG\/hw054YuA7\/xgBUUTEccgkZ1UDnYgcM0gkD3a3Y\/Ea1vEmIuUhUxJxJ+zavZcQm6BID3WHeJLxUWatdxMoTmFz2Ovt1UfVrT2Ds28rVvie33NVxy7av58YDqteSygwMHnUcN24c9O3bNxPeGwlNRUVFZtbFR+LNjUvIHQ5zJV3yyo7QTb3yTAsuJW2e2Y8FupOzV4st4CTmrdvP1WqYqwOBi3pR6VRTB7kwh7PdXiq7aHWYBEJezwTgSY+6iqVtX6dD0KykVAiMnEogKHaDJzBmhowrHdThXt++j23qDdoPg2Sl1jGFWdVh0ke8MGfS10cUKVMO6JAYVwcCF\/Wi0ikXg1xQGpG4vpv0OdvtpbJL0naqnvd6qhDSv+8qlrZ9XR8hWkknCExQwjU8akkVmTMpZFTGd9G5VTrhchJeqkB3YuLHz748mwW6O71XQxbsjl\/rtu5hCR\/xWjioQ6RZVHoltWnc513Ui0onKj9XYSvnPVPJ27pvu71UdrHVfl6u15MOYVextO3rdAialZQKgUEVcRDDC8OIy5t58SsNj1n37t07MAmbWROTS1MZ30XnptSJJ35sOvYFFhdGjtaLlkAS0\/7OV1jKgSgSQ6lXcg\/4oQQX9aLSicrPo\/DmfXvlypXVxHIdNsF2e6nsQum\/QWV5PekQdhVL275Oh6BZSakRGL5sVFZWlpXIjQe5mjBhghPkBeGkMr6Lzq3SCZeS8MLlJJ2LkxiMEbO2fBcjMb946gRoWfTD8yKJ4XmT5BgxKr10dLEh46JeVDpR+bkN3G2Uabu9VHax0XaxTK8nHcKuYmnb1+kQNCspNQJjpma60lTGd9G5VTrxTb3n1H0X3l07GC4p7A6jO+7f8xJ0YfZq3CeDF27qxcB3GCdGJjE8Wi\/mTZpy5mKo+uDurEB3Kr3S8ggX9aLSicrP07KNab2220tlF9N2mcp7PU0RC5d3FUvbvk6HoFlJnsBo4EVlfBedW0enqE29YSRm3Y2tMieTMPjdv+dvCyUxYsoBzGD93e5PYeeunVC\/Xn04vPmVWoHvNMxIIqKDF0lFBoVQ6UTl5yrVxdnX7t27Q\/\/+\/WH69OmJEzmq6pXv224vlV1M22Uq7\/U0RcwTGDrEkpWUCoGJWgcXm5PrNfEwKKkGOhcHCl2ddCP1cgx5jBh+vJqnHAibicFAd+0qusN3uzewmZj1O4+FVq1aQdUH0xmhwazWLly6eOVSVyqdqPw8qu3iicOWLVtmNupTJXM0wd12e6nsYtKmOLJezzioBT\/jKpa2fZ0OQbOSUicwrpCUKNiojO+ic5voxCP13nzcizDutZLAQHcijnK0XiQxX67fy1IOiBemHFj66t9h+kXbYd+WV+H42l\/ALWumQb169aBPx+PgzG8fgUOP6uTETIwJXmZdMb40lU5Ufh7VEvEY9ZYtWzIEZsOGDTB+\/PiczsLYbi+VXeJ7ht6TXk89nHSkXMXStq\/rYGNDJhUCIzaEb9rF31wlM1TGd9G5TXUy3dQrRuvFvElhJGbxE9dCcfkvAPfETGzxe\/h6ZwUc\/dNXAMnNuq1VMPmUh52YhTHFy0anlcuk0onKz1Vt5qcOe\/bsCU8\/\/TRbQho0aBDgchKeSszVZbu9VHaxjYfXkw5hV7G07et0CJqVlDqByfoSnzQJZs2axX7CwYznSDJrEr00lfFddCh\/QU4AACAASURBVG5TnVSReoPQ59F6MdBd7XZdAg20asnNUPxKH3bMGpeTWq\/tBwW1a7OUA0hiLt37P9Du0qfojWtYoilehsXHEqfSicrPdRohfrigfBqnDm23l8ouOngmkfF6JkEv+1lXsbTt63QImpXkFIHhqvM9Mvj37NmzoVGjRmatIpamMr6Lzp1UJwx0t2n3RrjnglnQpG7TUOTxZBISmTASgwSmbdd7gKccuP+n++Dkyl+zpSPc\/\/LVimHK\/EnEZg8sLileNnSk0onKz220USxTFaU7LESDrJft9lLZxTaeXk86hF3F0rav0yFoVpIzBEYcdLAJPhu1mSHjSueyw4nReuWZmD89+1e4+oQKKDj5t6wp\/1WyDN4p6cz+jRt5f\/tcA5h10zVxm0n2XC7x0lWaSqd8GOSQvPTr1w8aN24c+nEjBsnEmR78O+hDyHZ7qeyi6wdx5byecZGr\/pyrWNr2dToEzUpKncDg4MKXjVwiLSKMVMZ30bmT6MQ39S469whtrxMD3YkPDSpdxfa5HFzneEZiLpr6CpzU4CuYdv47bFlp2JvdlOkHtJVIIJgErwTVRj5KpROVn9tqJ5KRoUOHwujRo2HBggWBG37lvGpBaUq4frbbS2UXW3jycr2edAi7iqVtX6dD0KykVAiMeIzaVdLiCYyeI3ES8\/Dpu+CqRd2hQ+OztPIm7f28AsRovQ0H14MdZx8MQ87aBHu3LIf3P34H3vh0D1Qe3B4Ob34VSz1w3kkN9JSyKOXiAEWlUy4GOd73W7RoEXuPm2oJCTPbl5eXs\/IxZYGfgYnuEFT+Y7HbsaLzQU9XdcxF37Zt\/6DyUycwUY125VQSlfFddO40dJKj9WKQOwx298rPv4HzLz0aRnQrZIPV8i8KAGdmkLws6PEJ2wcjpxzIdadJAy9VG6l0ovJzlb58GaiysjIjarKJV0VgsFAkMaNGjYo8DIDtxWvJkiUqlWPdx6PhzZo1i\/VsLh\/yetKh7RqWXbt2zTRu9erVdA11pKRUCIwjbddWg2pgp3rRaCuuIZiWTkhixGi9POVAg58UwJHND4NdO3dCvfr1WQv+9sLnMOix1iwejJxyQKOJpCJp4RXVCCqdqPw8DuBIOObNm6e1aT+KwPglJHP0qfzHvGazJ\/JBT1d1TLNvm1nZTNoTGA28qIzvonNT6cSPV5\/doBYMavYB4Omk2zqWwMWF3UMR5jFiMD5M40FzYObZ78OX679hKQe+O24zi8SL14TrP4C76lSy49Xtd46Dr9c\/kdpMDBVeGm6nLUKlE5WfqxQPmoFp2rQpzJkzB1q3bq16HKIIjLznRSY0YuG220tlFyUgCQW8ngkBFB53FUvbvk6HoFlJThOY559\/Htq3b++PUZvZ1EiassPxnEk9mhwGv29bR0sPMVrvOxXXswzWGK337FG14JwrT2Jl4OwMEph12\/bAwhs7QMNV\/eDbqg1Qr\/1dOY\/OS4mXFkAaQlQ65WKQowiRYDoDg6eWpkyZAp06dcpC03Z7qeyi4QKJRLyeieDLethVLG37Oh2CZiWlRmDkfEhz587NDDA6RyXNmplMmsr4Ljo3tU6cxCCBQSKjc3ESs7bVcrhw8jnskbvbvwO3vPVfGQLTedixmRgxb91+Lux8uZiRmAY\/XqZTBZkMNV4UilHpROXnFG2KKkMmMGJ6ApzBkceWsP01tttLZRfbeHo96RB2FUvbvk6HoFlJqRAYPq2L08YYOpxP+5aUlMDDDz8MZWVlPg6MmR1jS9vocGL26n98dCcsqihT5k3atmAsLLtrE5zbdw\/UOqYQ\/jK1E3y9ryX0GTEfXlt0CvxkTi92nLr9na+wTb0LB3VIhcTYwCu28f7zIJVOuRzk+CZb3naTTbxJ8eLP224vlV2o2htWjteTDmFXsbTt63QImpWUCoGRv5hQZR4PxqUUAtQDnYvObUsn0+zViPWin02CU+qMZNF6N9ZpmdkDg0tIOAOD18sfb4fu961gOZPu6VmfnUzCC1MO5OKyhVcS3al0ytUgh30dY7rw4HJ8xgSXeHwupCSeEO9ZKv+JV7v+U\/mgp6s65qpv61uTRtIZAiPGbigoKKBpHVEpVMZ30blt6nTHqt2wcNM38FCHeoCbe1UXHqX+bvn4\/WIX3QINGzRk\/8R9Mb986sTM45zEzCxuC73a7YPt\/+zMUg7kgsTYxEuFj+0vaCo\/j2pH0McLyuscjY6LT9hzttvroq8EYeH1pPMsV7G07et0CJqV5BSBQdV79epl1oIcSFMZ30XnzqVOPG\/SgkvKAq2Gy0gNrx7D7q3+TXM44YH17N\/i7\/xBTPA4afEadjIplwHucomXrmtT6UTl5yq9\/QyMCqHc3qfyH9ta54OeruqYq75t2wfk8g9oArNv3z549tln4cUXX2QbiC+77DKoVav6TAGV8V107lzoxKP14kxM09oHh\/p4FoEZeS7U2rGJResNIjBYCCZ+fOmT7TklMbnAy3QQoNKJys919Pd7YHRQyo0Mlf\/Y1jYf9HRVx1z2bdt+IJafGoEZMGAAC\/MdddmOxLtq1Sp48803oU+fPvDMM89AkyZNoGPHjtVUojK+i86dK510Ug6IRGXNm8vgkId+DYceUwiY+HHH8UPg8ctXw6UzmsPpvfYvLXESg8er8WRSLq5c4WXSFiqdqPzcRPc0ZW23l8outjHyetIh7CqWtn2dDkGzklIhMGYq6knjF11FRUXWJkAxYFZUzqXvv\/8eHnnkETYL07ZtW09g9CA3lkISgxcPdCfnTdo8sx87gYTX9m3boN7eL2Hn0kfY3ycs+B54ygEMdNeyqF4WicE\/8GSS7cvFAYpKp5o6yIX5hO32UtnlQPTpoDbnA56u6mjb1237aFj5NYLA8OlokaTwkw14qgFnckaOHAlFRUXsv7Vr10L9+vWhXbt2DBd8vk6dOmwJ6aCDDvIEJi1vFOrlA4EcrZenHBBJjHy8uuqD6dZSDrg4QFHpVFMHOU9gojs0lf\/YHjbyQU9XdaypfTvvCQxuCMSkcBix9\/PPP8\/MwMinGvgppzvuuIORFNzrgqQF486ceuqpcN555wWSF\/b1f8IJQJEIy0XndlEnxFzUi5MYTDeAaQf+fPknsK78K5ZygM\/E8JNJI7q1YskgMdAdZrWmTv7oIl5UOlH5ue0XGVX5tttLZReq9oaV4\/WkQ9hVLG37Oh2CZiXlPYHhzZWXkDDORGlpKUycOBHwWDb+jWSHx57A59544w2YNWtWZtmoZ8+ejKzIF5XxXXTuNHQSUw6cW7AExr1WUi3QnawXJzEYIwb3xCCJwZQDl85oVo3E8JNJPFovZcqBNPBSdWkqnaj8XKWvK\/dtt5fKLrbx8nrSIewqlrZ9nQ5Bs5JqLIGR48oEERhdqDipWbJkie4jgXKupVpHJdPS6Z09teB3m+pD7wZ7oFeDqmp4Ber13N3w\/XN3w0HXlwKccA4suXE3fLXxe+jxVN3M8yX\/\/ALKVu2CB65oAmcdXxvqrr4VDt67CXaePBcO27YYDvnqrYzs3oY\/hX11zzCyaVp4RSlJoVPXrl1ZFRQzjUaApihse1B39WUmQ+71pHNCV7G07et0CJqVVGMJjM4MjC5UVMZ30bnT1EnMm\/TO+vFZKQfC9MKNvrixV5yJQTuKge7weDVP\/NiiUe1MyoHDm18JBSf\/NmN23Cvz3e5PoW6Hu3RdIWtpS\/shy4JUNqTyc1VzsW\/iyT\/5sn3qUK7Pdnup7KLCM+l9r2dSBH943lUsbfs6HYJmJTlBYHhuJMyBhKkE+vfvD9OnT4epU6dqZ6KWl5DC9sDwJSUTmKiM76Jzp60TRurFiL0YI+bRt26ETbs3wuiOJdBgV6NMKgHZVhtLLoS9n1fAMYPmwNdHFMGfL18NRzY\/NJDE4PHqr9c\/wVIOBEXrRRKDv9c6OjtTcZh\/pI1XkF5UOlH5eVTf4pvre\/funXrQStvtpbKLyVgVR9brGQe14GdcxdK2r9MhaFZS6gSGkxc8HdSyZcvMvpWFCxdCeXl5Zg+LqlkygQk7hRQn0i+V8V10bhd0CsqbpNKLkxgMdBd0iSeTSovmQsHJt4SmHEByozsLo9JL5ac27lPpROXnKgIzZMgQGD16NGDm6DQv2+2lsottjLye0QgPKl0FLRpVT2+zbmsVYDoT8XIVS9u+bttHw8pPncCIuVG2bNmSITC4rj9+\/HjtWZgkcWBU4FMZ30XndkUnOW+Sjl5IYvA6ruSFQBPyk0mTT5kN1\/7mz7Dvi+Wwo7wYDm9+VRZhwc2+unmUdPRS+RP1fSqdqPxc1b6gvqp6xsZ92+2lsouNtufDS1dud1p4YtoSPNkoX0G\/p6Wjykds+7qqflv3Uycw2DB+FBpPAT399NNsCWnQoEFsOSmX2WnDQKYyvovO7aJOaIfLFl4MzY9ozk4nhV37NlfAuhtbQUG7LpEkZtXzN0HdDlNYMS99vC1TXJ+Ox8G5TXex5SVPYOjCBUQNVnxmNCgKt+4eGJ3Ejzy7PeoyYcKEwOUqqn4d1l5X+5YrxMD0pZYWnp7AmFoqd\/JOEBhsrryxL2zQyR00P9RENdCl1QGjMHNNJ55yYGaTLaF7YMT26JCYpa\/+HfC\/zY0HZk354sB01r5HoMs5P\/V7YAjjHdnsozy6duPGjbNCIoh1iicQcYk6bMmKql97AmPT4j+UndZY5QlMbuwbpxZnCEwc5XP1DNVAl1YHzCcCg7oiiVm36xuYdMoawAzWlxR2Zxt7wy4eI6bhz0tYNms55QAuJe16ayi8UlkPLrqkhGWv\/m73Bra5F+\/hySTdjNY12YZUfm6rX+JHztChQ9n+mQULFgQuLyNhGTduHPTt21e5x8Z2e130lSDbeD2jPZYTGNwLU\/raJvYRVNyxCfglJFs9Xb9cT2A0sOIDndzR8W+8WrVqxf6vur9s2TLo3LlzpkbT51Xlx7mPe424TnGeN2m\/Sfn\/vfRzaHr4wTDkpLWMxFzR9OdwedOfZ2GdVfeSefD9H4uBR+st7b8K1vxtL4vWO3ntehjaqTbMfv5FOOHw96FXu29h566d8G3d9nD8mdcDDkx4X25L0N\/8tyDb6Twv+gqVPNclaXkYC8ZGHBi+bLR582aYMmUKWzK2tYSEdY0dO5aRl2nTpjGb+iWk6EHOExg1gTn\/xAbQ\/b4VgGEZzjuxASMxnsBovDwti6RCYKLWwcX26q6JW8bIpxKwDXBA+XEGVYwPg3FixBgxmHJgw\/UHw31j\/4vVgjFi8OKJH01TDsTRyzZ8VDrZnpGgwiFqDwwfWzAxK+6fQ9nhw4fD5MmTq83I2G4vlV2ocAsrx+sZjTCPK4VSLRrWhpc+2Q6YsuTlj7dVSyDrKpa2fd22j4aVnwqBSauxceulMr6Lzu2iTmgn1GtLw+ZsOen6VrXh+sL9MySqKyjQ3dsrd8LoF0+FI5sfliExGOhuZu+20LxRbXj3mcvZZl7Mm3RwnWaRVbiIF5VOVH6uslHS+yoCI+55EcM0yCEUqCJsh7WHIkJyUqx0nvd6RqOE0b3f+HQPzLq8CTQ9ohacdU8FlPz4aOjetl61B13DkkfYRkVtzK7q+JdNGU9gNNClGtipXjQaKmuLuKgTJzC43CLGiPnHR3dmResNayQer656d2nWTAzmTTq9V0P2CJIX3PuC1xFfHQzd724Op392A3xbtQEa\/HiZJzDa3pOOYBSB4YSluLgYcBZGRWBsDuqu9i3Zal7PcD\/moRj4vheUxBkZ\/PCRY8CI41Y6PSO8Vqp3mGvtSp3A6CwnNW3aFObMmaPclGcLXCrjuzhQuKiTPBDIMWJ07CwGuvty\/TfVovVioLse962AZksPglF\/PCVznBrLjjpS7SJeVDpR+bnKPklTCaiOUYtxZvjG36Dxw3Z7qeyiwjPpfa9nMIJ8jOB7XrgUbuTFfXNbp+2PQyVermJp29eT+mDc51MnMKg4buorLCzMitUgDkL8WOTdd98dt52JnqMyvovO7aJOQV8y\/Hj1onOPYJt6ecqBDo3PCrWtSGKeGbyenU466JRD2J4YvM4\/qSF899TXcPOe1WxNe1jnWqHRenklLuJFpROVn0d1Nv7BcvPNN8Pzzz\/PTgs1a9YMRo4cCRiNWydStkxgxGCYPLqvGAdm7ty5bDZGvmy3l8ouiQYvjYe9nsEgIUl5+ZPtsPDGDmzzrng1uvUF4FnvPYHRcDJLIqkTmKDBB9sqDlIYoRej8uJXVBoX1UDn4kDhok5BBAZ\/QxKDF+ZN4iRmwSVloS6BMWIqx1wIhx5TCB\/XnQen924IM89+H1oU1c3kTVp212ewocv37IsKj0be07N+JIlxES8qnaj8XEVg+B6VJ598MvPhoppVsdHvbbeXyi422p4PL1253bnEM2jpSNSn\/Z2vZE4j5QOWtn3dto+GlZ86gUHF+NcS\/1LiU8wDBw5kJwn8DIw998jloGDSiiC9Kvd8B5e8sgPOblCLkRidiwe6+7jePPjJnF7VHkEC03nYsWxPDB6TxOOST\/b6OjDlQBix0tHDpgyVDXMxyPE9KbgsfMEFF7Cs1A899BCLwL1u3brQ4HQ28LPdXiq72Gh7Prx00yIwYUtHoj48Joy8jOSqzW37um0fdZrAoHI8wmZlZSWIe17S+DKTwfJxYPTi3MidF\/\/GSzdOjvg8xszBpQU5bgo\/mdS7wR7o1aAqPC6MUDcGuvtHv\/nwX+2fgYYXXct02r5tf0qBj145lxEbrLtyxz7o\/ugGNl1cdtlnUHf1\/7CcSZg7ibeFtyeobXJbc\/W3rBsnWqb124oDEzT43HvvvdCtWzfA2VUkMZg2JE6m+CQDp+1B3dWXWVrEIImtcvnxELV0xNvA98FgtntxeclVm9v29aS2jfu8EzMwcZXP1XNUxnfRuV3UyWSw2vRVJVy1qDvgXpiovElPFS+DQz6ZCw2vzo7ou7Z8V2Y5CevFr6+bSlexk0p4zFqO0OsiXlQ6Ufl5rvpl0npst5fKLknbqXre6\/kDQqqlIxFLXEbCvXO49MwvV7G07esqH7N13xMYDWSpjO+ic7uokw6BEY9X45KS6tq2YCwUnHoBVJZcmDlejc\/g75h+QL7wqCQGrJI36rmIF5VOVH6usoUr9223l8outvHyeu5HWGfpSLQFztTgJR6ndhVL275u20fDyvcERgN5KuO76Nwu6qRDYFAGN\/W+vn0f2w9zyN6VbGNv\/3YDof+pv6lmVZGorL76oEzKAZyZWfX8kSzlQMui7H01cu4TXb00XIpUhMqGVH4e1TiXwibYbi+VXUidJaAwr+d+UHSWjkT4go5Tu4qlbV+37aNOE5ikcSFsg0dlfBed20WdTIiCeLx6xeY3GInBpST5eLVIYHjKAcybtG\/zWvjb3F8CBrq7dEazaiQG851MWryGTRWP6FZYLd+Vbd\/TKZ\/KhlR+rtLZlbAJtttLZRcVnknvez0hs4lfDFingysepxafcRVL276ug5UNmdRnYPgXWe\/evbViQNgAQVUmlfFddG4XdTIhMCrb8fvyUhEnMbUaF0KL+9bAny\/\/hJGYQa+fUq1Ivi7+ePvJ0PGY7XD0T1\/RrTYnclQ2pPJz1QyMGOqfy6YRNsF2e6nsYtuJDnQ9+dIR5jniOdJ0MZej8rqKpW1f18WLWs4JAhM0oFE3NEl5VMZ30bld1MkGgcEcSbWOKcxyg32fVwASmfpdrmFLSrokBvMmqVIOJPE302epbEjl5yr9XQmbYLu9VHZR4Zn0\/oGuJ86ylr62MTBgnQpbeRnJVSxt+7oKJ1v3Uycw2DAx6q6thiYpl8r4Ljq3izrFITC4F0ZM\/IhLSbikFLScJPoCHrHGjb0F7bpAnQGLoWzwenb7l0+dWM1l8Evtugeehd82fYDFi4lKOZDE30yfpbIhlZ\/r6O9C2ATb7aWyiw6eSWQOZD1NTh2FYYzLSPw4tatY2vb1JP6X5NnUCUzUpr4zzjgjp4GtwoCkMr6Lzu2iTnEIDD5jejKJ21smMXK0XtEvXl75Efz59bUwrvF1cOhRnZwgMVQ2pPLzJANSLp+13V4qu9jG5EDVM8nSkWgTMSqvq1ja9nXbPhpWfuoEJq2Gm9RLZXwXndtFneISGHxOTvyok3IAn+MRe3Em5puflMHjl6\/OSjnA\/YXjNe7xeXBD\/VGMxBzW4krY98WrGZc6vPmVUOvo6rl3THzORJbKhlR+bqJ7mrK220tlF9sYHah6Jlk6Em0iHqd2FUvbvm7bRz2BSYAwN77snPg3XrqRZjG6bOfOnTOamD6ftP6g5zds2JDRyUb5JviI9YdF4pUj8wb9jUtJ63Z9A39oshNaNj0Uxr1WAnv27IHRJ4+Njty7bQMc8tCvWe6kdac9CEtu3A2dhx7LUg1wW\/H24N\/Pvv4y9D78d7Dx4PZQp920jC9UfTAdvvz8fahqNkIrUjAnbDJWun\/LusUtL5eReBN0R7JHbQ\/qrr7MZAAPRD0plo44jmJU3m+\/3Jjp82SOSlCQbV8nUDFWEU7MwPhj1LFsR\/KQq4NXUr3ExI8mQOFMDG743ft5BTudFDXY7\/tiOXy9\/i\/w9fon4K36t8F\/X7g\/TQFeSGJwdiYXMzFJseI619RBLsz+tttLZRcT\/40je6DpSbV0JGLNj1N3OvqH9CZxbGHrGdu+bktvVbmpExi+B+bmm2+G559\/Hvr27cty4IwcORKKioqcOFpNZXwXBwoXdeIzCHx2ReXEQfflxI+6KQewLJHEHDNoDtRu1yVThYjXVyuGsVxJSGDw30cUlWYRFn4\/jv4mz1DZkMrPTXRPU9Z2e6nsYhujA01PqqUj0S58GWlop9p+Bsa2wwrlO0Fg+DHqJ598EgoLCxlpcSGJo+6XKQ6E\/vIIuIrA6tWrtVSz\/UJHJfCDRRUHplGjRlr6JhWy3d4DjRgktYfqeQo8KZeORH35MtIbNxd6AqMyJOH91AlMVVUVm23BDNQXXHABy0z70EMPwdNPPw3r1q3Li1NItgdCQnv7og4wBEx800TWFEbx2HTYs7oZqU0+bjDmDF4jRoyoVq3N9mJlFC9cU5zjyB8oetpYOhLxxmWkkh8fDYMvOS2OGaw+Y9vXrSofUXjqBIbrdu+990K3bt1gy5YtjMToDma5AE5lfNX9XOjo6\/AIBCFg4psmsnHRDpuB0S2PE6HGjRsrP2743rqBAwd6AhMB8IFCYGwsHYmw4nHqM46tBY9e11HXnXMml4u+nbPGCBU5Q2DSaLxunSrjq+7r1iPK4ZTkSx9vy\/zUp+NxcN5JDeIU5Z85gBEw8U0T2TQgRUIydOhQGD16NCxYsACmTp0KYctNSJTGjh3L1MTZXT8DE26xA4HA2Fo6ksds3AuzddqFaXSPyDpd79txAfMERgM5lfFV9zWqyBJB8oLTnZg8kF\/49bBua1VW6nbTcsPk5VNg4herHGgwKrhgWDl8mbCsrCxLBfnLOOg0WtjXc1BbcLmgsrISJk6cCAUFBYnh4XoXFxdDp07RcV1UbQybDcAo1KNGjYK5c+eyOkwxVM1UmvimiWxicBMUoLOExBNGVlRUsJo8gTlwCYztpSOOrHicukWj2gk8nP7RfOnbpi33BEYDMZXxVfc1qsiI4JfCS59szyIvIonBEPaUMzGcoOAAjy9QObmmvIcgjCREldOjRw+tU2X48sbyZ8+ezb6sxf1RQS8gEVf+xV2\/fn12kq1169YmsAfKxiEw4sk5fNEOHz4cJk+ezJZGsW3t2rWD2267jREskfQggWnTpg0MGDCAvWxlW+hiKDfExDdNZJOAmzRsgorA4P1HH32U4TxjxoxIAoM3lyxZot2ckn9+AU3r16omX7lzH9v\/IF4YYwlPVLp+1XQ9H3h1O5S9vwtmXd4Emh5R3XaU9rl49lq4sego6N62HmWxscvC2E780t3MH7uyFB70BEYDdNXArrqvUUVGBKcgMT172KW6b1IXyoovWfmlz1+wOsfZKcqRCQzXb\/z48Wy5AF\/6uOEbZ3JwWWDOnDkZooLPvvjii2wjOP4\/jPCIpIGXwY\/t49+zZs1iECKhwNkmrO+KK64APCHHZ2KwrtLS0qyZniCsRAKEZeIzWAeWh1gjZrj3a8eOHYAzTUcddVSG8CSxhegDJr5pImvqZ1yeIvt8FIFBzMeNG5chsdSbeHEmVJwZ5e0K+v1AWJqJ6wdxnouDZy6WjsS29H3wNahXr56VmfI4mPFnctG3k+gX91lPYDSQUxlfvI870V2+gtZn+TIG6j1hwoSs2Dvi6RHVckVYOUHLK0FLUUEERtz0iSQCLyQnsizfBI4kAPc+jBkzJnB\/hPhC42XgszhLwvdKYDvKy8vhjjvugN\/\/\/veMuKxdu5bVjUf8+fIE\/ptfQQRGnoFBAtO+fXuoU6cOKwfr2b17N7z11lsZcmSCIdatWmJT+W5cshPXx5Nu4pVJrbwHJuy0UxBOJtgEEZUeM1fAvcVtAZcLPIGJ6xH6z5kSmFwtHYktmLHobcBZOp7cUb91diXj+LpdjWhK9wRGA0du\/LBQ+5Qh2KNmWLBD3lS6ChYO6qChtbmIuN+F78kQS8EXN85QqF6acjl8JkM1kxNFYG699VZ2vJ7PgogvQiQtDz74IAwePJjN0qCeOBMj71uRl4R4Gbxsrh+fYREJDNaByxIoO23atGrLVGF7YDgh5GX279+fbUDl5Vx99dXw8MMPZ9rF8Y6LoWx10XfxXlQqBko\/jvK+pNnnVUtIss9y0huGjUlPEYkKnjrZv1et1X5iLexZw79NX7gmelDK1lQ9bZ86CrIBYnnWPRVsBqa4YxNKMyUqyxOYRPBlPxyVgVqUzJds1JTOkes9MEFmjXrB8JfHTTfdBLfccgvbNBs2M8PLQWKhE1k5agnpzjvvZLMkQQTmzTffZBthxQtJFq+XLzndf\/\/9gSRIh8CgL+LSBJIcnJ3h+1h4nfIMjLykxgkMzgwhAeLlYN04jXp5dwAAIABJREFUYxS2UdgUwyQvaUo\/DhsuKLLPywQmalbH1hISzr40b1Qbzj+pIUxavCZDZEQSU1OJAeGrwKgoEzxzvXTEG4I6Tlm+h\/0ZtRXAqOEEwrno2wRqGhfhZ2A0IFMZX3Vfo4osEZyFadGoIPNFh195uMM96CvPtGxZnh9NFfeT4EsTr4suuohtKu3du3dmWUncICme9IkqR3cDqmoTb9Tyj0gAdF9o8hJS1AwMPyGEJOiSSy6pluIiaAlJbM+HH36Y2TezcOFCWLRoEZx66qkZkoX644XHhINsoYuh6wQmqb9SPh+n3+JXPX5Z4+xL2Y0dMhvqsc9iH8XlpOKOx7G+a\/LCpWyXaVk1Tc80lo5EArP8iwJw7Th1HF839aM05D2B0UBdZXzVfY0qqonwmRg8Oo0XfunZmpIU911gXeISkbynIGpWLKwc3c3AqmPUQRtwUV++yVfcDxG0TwVlozbxqghMFDEKayPqge264YYb4G9\/+xvb+IunPvr16wdTpkzJbBTmBEyFoXwUXd7M7CKB4bMumzdvZm1GTFauXFnN53M94xqn3+KLaf3WPbBu2x5GVPj18sfb2H4YJDE4I4NEpv+Z9ZyMyioDX9MITBpLRyKBwWVa3AspEtw47wDKZ+L4OmX9tspygsAkPVZpCxxersr4qvu29fPl5waBoNNHuak5fi0mvmkiG18jd56M0178usfZl6g9DnzGlBOZmb3bkoY+oECQz\/JiWdu3bYMGDRuyYm3FmqLQOYhoie3g+iOJxI+9NJZwuI7oI+ed2CAVHYKwjuPrFDazXUbqBKYmZKOuqc5h2\/nyqXw54Fy+6G7imyay+dL+KD3jtJfPsOApE9X18sqPYNJLu1hcJ4zfxE8sqZ7LxX1xM7JIDMKOiedCJ1UdQQRG1FdcOjrvpIaBx91VdSS9z3XkS4quROWN4+tJscjF804QmJqQjbomBgnKhQP6OuwiYDJwmcja1To3pcdpL9+8q\/N1z19m\/PTg\/gCVrVJ5scqI1kQCIy4dIdEMitdj27O4zV2LyhvH121jRVF+6gTGZ6OmMKMvwyMQjIDJwGUimxRvea+PHH8oafk6z5u2ly8f6e5tkGcMcF\/boHmrAk8s6ehLKYMzBJhfDffyPPz\/1kBJz\/3LXPk6A4N645IdX9pLqx2izXEZCQmrrb2LJv5g6usmZacpmzqB4Y332ajTcYOkOXqCtA57OYUFGZNjzvB4M2LZQXFpwhCjzolkEnxN1UaOjfzCRp1xgy4\/gWSKYRgBMBm4TGSTeCvf2MxTRvBlZDzppUoZkaRe+VnT9vINvLpxmMI2x\/KXrXhiibJdQWUh+Xr5k+0sQSxuQsbZILxQB8ygvPKzfYB7dcLSmNjWT6f8qD0wpa9tZMHj+OUCgUF\/wUtntk6n\/UlkTH09SV25fNYZApPLRpvWpTK+6r5pfZtn9oNax\/yQyJE\/v+\/zCmg8aI5pcZHy\/Khv3Bw9YjRarEiOf8Jf6Hj6JCpUvqikHLtDLEOVVNFGTiRTAsPzH\/F0ADyyL55AwiPU8+bNAwwah7F08BJP6SCBwUssIw6GHE8T3zSRjeuEYViaBKeLW3cSAqOzeVcuP+p0j7zRlx+9pmgblo0XJyz4f\/yNJxjEzaV48Qz3qOcti\/bnYMO9OroEjUJXkzJkPPnSHD8RJie\/TXMJCdvl0jJSLvq2iS2pZD2B0UBSZXzVfY0qskS2LRgLDa8eU+2xsN9Nyxfl+cmauDl65LqDgtFxmah8SVEEBu+JwfWCjkJzsqCbE0kMqMaP72I9uB\/r+OOPZ7mQ+BFlJF74+69+9St47LHHWF4mPLIddFQ7qI3iy\/m5556Dt99+GzDp5HXXXcfK4TbAnEijR4\/OJH7kMxQiNroYukpgUK98nIEx2bzLsdc5nhx0YmnuaxtZHCj5CjshhGXwWRU8zi3OrrRoWBtwQ6ucBDboFBLOYvDIwmm8\/FXjmIgn37CLzyAhk\/FCHNIgYrLN8Ti1C1F5qd9RKlvl6r4TBCYoDLsq706uAMJ6VMZX3TfVNQ0CEzdHT1DbxCUgceknaHklKC1BUPRU8QizmGE4Tk4kOWYLnyHBWQ+cFeGB+7geSDSQwGDUXJTFbNec1CDhEBMvBhEMeQamoqKCwcbTHWA9p512GksxwMszwRDLCltiM\/FNE1lTn5bl820PjMnmXRMCw2XxhYx7OJAohZ1Y4ssiUctB+DLHmFH7icv+mRbVJb50eQRbF2dixE3RPe5bwZq18MYOmZklVTtzcV8mMHH8xoaeuezbNvQPKzN1AiNu4uXr3\/w3VBqn3cWIr7kEh9elMr54f\/XVB6WhonadJyz4Pks2aY6eqCUdTliwwqClkTClowiMmJ8I646TE0lexuAzJJgeYNiwYWwfBpbNZ304gUFygWkL8DrzzDNZbiQ5pUAQSRODzfEyedZsLBtTCWDaA0yXEESIMOidKYa6vivaQOXn2k6WJ4K67eUvdd3Nu3EIjEhk8OXMZ0JwAygPnMdTFqAsLgchSeHpDJJsFA3abNz9vhWsDp3j4rkyN+pZ+W1DQN1cJFiIg4wlX0ZK+zi1rq\/nypZU9aROYFxaEw8DlRtfdk78Gy\/qJHhpzMDEzdGDL+I+ffowHMKSPPKlFnzpy\/tDgjAPIjByTiCTnEhXXHEFi3zL8zahDrfffnuGLJgQGNQXicuJJ56YySottkGegRFnX5CIi6QIiUu3bt3YkpJIksQZHV62KYYygeG+mkYyR53cZ2nMuOoO6qabd5MQGHwWZ1rwBc1PLImEpfhHxxnNrui8KIKWusT9Ja4E4hs1bwXM+pc7R9GDsA3C0oWovLq+ruMvLsmkTmAQDBzUcWMjX\/OXswanDZjK+Kr7pvqnQWD4BlPTHD3yJl7ZlmxA\/k92aKpNvElzIqmWkKJmYJBc4Im51157jZEgmWwELSGJ+op7efD39957j+VWwrxTPB4SzvKI\/SEOhjKB0fFBaj\/WqVMkaPhv104hmR6dFtusswcmkMAvrgiMYWLrZE2UnrgEgntqTGefTGyvkhX3CbkSRydM5yAsXYjKm2bfVtk3yf1UCIzOFxk2Kte5UcKAVBlfdd\/UQGGnkPa8uxSOK3nBtLhIeXFvSZwcPYED8KRJMGvWrMwtfsRXdwOq6hg1RU6ksE28mLxSRWCiUgoEtVE8JlxYWAi4BwbrEBNgyntqZAxEDPlskoh92AyGiW+ayJI64X9OrwXltKKuRyxPp71xNu\/yOmoCgcG28KiyaZAHTl5wg\/HFJ9WGCb072HSJxGWHHfXGk2BpLsfp+HrixqdQQCoEJoV2JqpSZXzV\/USV+4edQyAsUaRzimpsQDd9odtqo6vHqJNswoxLYOT8PhxzWydrdPTkezlyufdEPGmEy1hND9kGfAnUlh8mLTcISxeOU9fUd5QnMBoeqzK+6r5GFV4kDxDgMz\/r1q3LLHe6rraJb5rIUrdb3isUVb6K7MjJYeOe0Iq7eTfpDAw1tqrydAgMlpHLE0oieeEnjXT1VLXX5v0wHdOOyptm37aJtycwGuiqjK+6r1GFF\/EIWEHAxDdNZKmUNSWF\/JRX48aNA0kkD2SIm9J5jB2cMQuKqaNqb9jm3Z1LH4E9776YgaB+l75Qu12XapAkeeHq1pHEDryOnbt2Qv169SGsHWId8qyI6qh2nHbwfUfyKagwPOPUYYqbbh1hOupE5dWtw1R3lFf5epwyXXjGGQLjQlyIMIOojK+674KhvQ4HJgImvmkimwaafM8QHjXHmDk8oGCULnz\/Ed\/XJMpGtTds8y6+ZPZtXpsVaBI33QdFyY5LYEzqiGsHsQ6uZ1g75DrE5JRRm3vjtIOnWQjabxOEZ5w6TDEzqSPM5qrj1CZ1mOrvCUwcxAyecSUyp0sE5uv1T8C+L17NqHR48yuh1tGdDFD1oh4Bsy8v1wkMt6dqCUm0e9TG8aj2Bm3exU30Ve+9GBolu+DUC7JmYuIQGNM64vi4XIeoJ5IYuR1hdfATSkGRZuO0I4q8oA4ynnHqMMXLtI4om4dF5TWtw7QNnsDEQUzzmXyKAxOX4GhCkRFD8vLd7g1QcPJvM79VfTAdvtv9KdTtcJdpcUp5ec+AGM9FPjEWdTIsrJygSMuolBw3Rn4+SCaqMdRJHOMc5w\/bf6FK8pi07RS+WdMIjHxcXsYoqr1Bm3fxdGBULjL5fhwCY1qHsnMHCKj0VOkgFslJBwbSE5MWqsoQ74snjaJyQsl4mtQRByd8xrQO1ZF0DDwoJ3c0rSNOW\/Klb5u2zYklJD8D84PZ9n2xHPZuWZ5FXvhdJDGHHtWJdCZGnmLnf8vh9Hl8jjCSEFVOjx49YOTIkVBUVARy3BjRYeW0AEFRmsMc3EYSR1MCI8fAETEJioEjblxduXIli5cjZmjmcWGCAtvpdnSTgctEVrd+G3I6MzA6voPtxWvJkiUZNb8fvv83neugyauzxL6f1Qfgk+U6jzIZ+Xn8LUn9cZ\/H8AnNmjVjOsWp\/41P98BvntwEZx1fG2YtLzJqf+WOfVDyzy9g48598JsfNYBL55yu\/TycfSUc9PPsD7o4+osVmj6PNj9o4NxE\/vP9\/w3LtMO0fhmsqOflKOz6QLsr6QSBQXj8Hpj9TvLVimGRsyyq+6auFjXFrvqCFeuiKCcoEaT4ssJItkiEysrKMokWTZM4BsWQwcEby8WQ\/zx+DZ5cwdkm\/B0j+T755JPAo\/8GxYEJIztIYrZs2cKgkqMQi2378MMPPYHRdF4VgdH12zDChrMveMnJAKO+lPdtrmBf62KcJuoZmKA6NCHLEouagYlbh7jxFk8NFZTeEDpbxevYO\/hZMMlpZDIDE7cdMp6UNucYYTwYnhkc6zOtI47N8+XjxLRtzhAYU8VzKa8yvnh\/68JW2qo16rE\/FYF4JX0ey4oqI6hOkTzygGlcJ3HpQxXuPaycoCWkoKWoIAIjLjEiicCLB4ETZywwOi6G5ceZDgzRz0+hyPhGRfFFAoNl85kRMe\/S2rVrWVE4g6SbhVpF8OQZGJ6SgT8XdvxX28EMTx+o\/NykXpuyUQRGl7ygfkHtjYq8a7pXIQ6BMa0jDs5Ue2DkusX0A1NO\/wzO+npl6H6hzcedAz9feuR+oqiZkDHf98BgW\/E4tbxMlgub50vfNvXnVAiMj8QbbqaoGRbcF4P3659XampnLXnRLkEvTx4dNiznEa9ELofPZJguITEytnVrJhP0Qw89lJkFiZPEUZ4lkcvm+vEZFpHAIDHCHEiYkXratGksI7W4tBO2lyuICPLfxCSPMnmLOj2jZcz\/CJkMXCayJjpQy8oERsQe6wqKVBzkz0HtVUXe5VGyG149hjULv\/R3Ln2U\/Zv\/xtsbh8DgsyZ1xMVWrAP1bF7voNB2mNbBN\/eWbJ0Ev\/rZuRlcOFal\/9oIt3\/bC+R9M6p6gvDMNVZJbR52nNp2O\/Klb6t8QL6fCoExVTJteZXxVfdN9M\/1Hpgg3cR8PfJ9\/vK46aab4JZbbskkSAzKGi4nYIxDYHh9mKl58uTJgQQGcweNGjUqS1UkWZjhWVxyuv\/++yGIBCEpwd+jCAySsHHjxjGZ8vLyalmow5aQOBm69tprAQkRtgGJjxy4LWj2iSLir4lvmsia+LSrskHt1Ym8y7+Y8eg0XrXbXQD1u1xTrZlxCQwWpFtHEmx5HdtWv83iwIS1I04dfHMvzsJMOX0TNN33GStm6L+PhbK63SBOWoIwPHOJVVKbR0XltdmOmtq3PYHR6J0q46vua1SRJYKzLAfXOT6zkRdnXvBkEl7iySTTcoPkxXw8fEYBX654YYJBzA3EN\/Tib0gocCbitttuA9yTwq+ocqg28SZN4oi6RpURRWA6derEchchCcLki0GbkU038QbpIm7iFfMyxbW1iW+ayMbVx6Xn5PYmjbwrty0JgcklTrb05DMxmH7g3uK2gC9vzGnUomHtavuLdNprS0+dunVldHQMO06tW0ccuZratz2B0fAGlfFV9zWqqCbCZ2Lw6DRetY4+Bw5vflWcopTPyBuoxSUi+fhv1DHqsHJ09yWojhJTJHGM2sSrIjCqZSIEWsaAL12okjxecMEFIO+BUS3VKQ1bQ\/fA6LRbR0but2Gbd3XKCpLReZnFLZvyOVt64iwMkpfu961gy0WY0BD3f+A1oluhcRNs6WmsSMQDOjrqROWl1AnLsvGOotYxTnmpEhhxLRtPaohr16oNo3EaG\/cZbnzZOfFvvLp27QqrV2cfqYxbl3\/OXQSislC7qrXou6gjT4bHfVf8+0DzY3FQj9q8G9e2Oi+zuGXz5\/hR8KTl+OfzH4God5AnMMT25V\/b4hfq+PHjWXhwflxWtWeCWKXQ4lTGV93PlZ6+HnsI8JkVipNB9rSsXrKJb5rI5rINtuoS26vavBtHh1wRGP\/xFMc6NesZVd9V3c9XNFKZgQna8CifLHDpa1dlfNX9fHUOr3f+I2Dimyay+Y9M9rS6zuZd0zZ7AmOKmJePi4Cq76rux6037edSITBBewlkAqMKVpVL4FTGV93Ppa6+Lo+AiICJb5rI1gSUeXtx8+6geatgZu+2oMqubNJuT2BM0PKySRBQ9V3V\/SR1p\/lsKgQGZ2DwSKocS0MEws\/A5MYt+FKevGFUXjKJypckaxoWVVmVD4iXw+PNiOWaLN1Q50TS2bwrYyC2Qdz4HLRRWYwHk7Ttsh44cK16cbLWBvCaOsiF9STeXurNu7w+T2CixzDVeMDHETm4JvYRjMY9Z86cTEgCMYwClw8rXy7PdKSNMx6Y1mEqr+q7qvum9bkinwqBwcbLcTBEQHRPreQKRJXxVfdN9Xxm8Ho4svlh1R77cv03cOmM5qbFRcrz+CPt2rXLHI0WT+ogcWjTpg07To1RavE4sZwvSaxAPm3DB5EpU6awKLlyOP0g5cTjxXhfLAPrj7ps5EQyGbCCcvCIMV7klAHYFrG9ctuTEnn0zdenfw9HFJUqc2hR+zGpo1ooDNu79PX3WHTUshs7kM6+oLr5SGDwhEyLRj+ER+Cwr9taVS0JYVKTBJ3ME98LCxcuhHnz5rFDEhh3Ci8+9mzevJkRGLzEMSXOeJO0HS48r+q7qvsutCGODqkRGPElKX79c9bt4imkMICpnWPZXZ9B52HHVqsu7Pc4hufP8BckzgJgzh+MBYODAIbm37FjB8sYrUs8sMyggGy8rqh8SWIb5Jc43hOD6wUdhTbNiSRGC+YzJFgPJlA8\/vjjAYkbnxnB9uPvv\/rVr+Cxxx5jG80bNWoUmlKAb0ZHGbHtWE4QgRFJyowZM9gjPHkmBYF54J5fw5nfPqokMdR+nMQvc\/Estnfc\/HKYtHgNYH4a6isfCQwefQ464hz2exLMgsYDcevAc889B2+\/\/TbUr18frrvuOtbneH\/AsWn06NEsz5iYUkTUR3e84eMWhjHAPt+lSxdWJ09ZIoY34B90OB7wAJhyDjXVR1YSzOK+g2pq306NwHBDyGkFxOl0G4aOU6bK+Kr7pnWmQWDat28PderUYQHakCzs3r0b3nrrrUzk26h8SXL7xGUQceknaEo3KNZJEIEJe8nLhEknJ5I8w8e\/+vBLDr\/05EzcOHjyAQtlcemTkxocRMWUAlEzi2EET56B4QklUT4q7o6OX6FvHtH\/cVh24QyW5TxqJobaj3X0S1MG2\/tfQ\/4C553UMFZcEpXunsBEI6QzA1NRsT\/aMcZJQmKAfeW0006DBQsWMAKDfc9kvMGy5OVoMW0HT3vC86Lh+MIJCf+IEscDjOAt51ALikqu8pWk91V9V3U\/af1pPZ86gUmr4Sb1qowv3h9\/7L+1ix79WfXU8Umfx8qjypDr5MSgf\/\/+bFDguX6uvvpqePjhhzMEJohwqvalcMKCzwZN94YBFUVgxPxEfDkLyQUOZkgqHnzwQZZCAI\/iYzl84BPrkpeE+FcfRhceNmxYZqlMHrCwDkxbgNeZZ54ZGpG4tLQUwgaxoD0wIkkJWkIK+8LUcTT+kkbZx9tPhm+rNkC99ncFLiep\/FynvnySafGji6FRjzvIN+9yDHJNYDDCq8vX1mkXZqkX9EEjfsDy\/od9+MUXX2SzMJioFfs3phaRPx7ijjfyR5CcSkXss\/jBJRMYOQCmJzC580JPYDSwVg3sqvsaVWSJpDEDg9mbMUkhz\/WDRAYHi+Li4swXiKikOLjwKdawyLE8pw++9OPsgcF65bxKXC+RjITlRMKlMR4kEZcmUYfbb789MwCaEBjUBVMpnHjiiZkZKxGXsNNzHIOWLVtmTXnLhEX+22T\/TZCfoW8+9v+3d\/3BXhXX\/Tg1DURplQxOQwk8tFJrNGqj1VBnQO2vMaN22liEZtRGDEHzcFoZ5Gn9WaGImEZpKnR8GkwjdGj8RYLVCRH88TRBjLTVmGIMvDAvVlSIQGCUse3nMuebffv23t177+69e+879x\/l7d7ds59zdr+fe3b3nCdeTKKhIv\/MFYf3JCTmiD94OklPceCt73Veu\/Sab9LqJ4ZPQMajLl5GZ5x+eqGw9i5zumoC4yKTrU6dW0i691L9gMBahGzz2FJiAqETGB5b3vUmjcBwLjW0C1KCMznwCAmBsVlRdeVCYBywthEUW7lDF7UTGJ6ga9eupeOPP76TCBFEAc\/cuXM7p\/6ZUOC\/ek4gPR8Q6rAnxPUsje0Qb9mcSLYtJD6sbPLAwGWNbaqNGzcmJEjdPsJY8x7i1WUJ4YFBoDM9zw\/IC3Jsqbm1Fl0xkbov\/ywddsrteU22cfVDRN7VQRACk20Wpi0k1f5VTwj+\/sorryR5yJCjjb2u+GjBQV\/OIVZkvUnbQmICgzVQ3VoSAhPPdBcC46ALG0GxlTt0MahK2i2kbX176HMPHZO3ucz66tmS7du3J54K3BjiCcuejqx8SXoH+lVg9VqjqwdGPQeC9tXtKh85kdIO8aq3rdIIjO1grSofZNevUetbQuz67u7uTr7wssaeV\/mqbeK6cP\/O\/fTCFw9JzsPoiUFR9+XH5tCHPnqG9cZSXjliq4\/Iu929T9Nbyy4MJloTCUzaLaRnX9vp3VNlyw\/W1dWVzAc+TMsfUfr5s6z1Rk1Pw4o2XRDhbSL9EC+ve\/g71gacC2QPbloWe9lCCjalhjQsBMYBaxtBsZU7dCFVGoQAu6hNGaljG4ZumyAx53\/oDpr5hX8ZIirXRTb0NnhhbMEwQ8\/bJhKY2Oy3anliC+HhOn6bLdvKXfuJrZ4QGAeN2JRvK3foQqo0AAFe3Pr7+we5rGMWXbdNbJ088W8z6e0Jt9L40SPomdd2dsRf8w9\/Q\/3ff4x2PzudRv3+ypiHZZWNvVpjxoxJ1VXoeSsExqqmKCroN2FjCuHhCpDNlm3lrv3EVk8IjINGbMq3lTt0IVUEgSAIwDZXzH5uUFwheFg+\/vVzafppv0F3\/uHP6d2+6TTyt6+iY\/56G3Vf\/uf0paO+2mgCg+0AbDfgkCdu1nHcHh3g0PNWCEwQk5ZGDQjYbNlW3lRQhcA4aM6mfFu5QxdSRRAIggBsc+beh+kvHzqaJkw+POnj9S3r6f7HHqa7tl6QRKA9\/YhXExIz6x8PoVlXXUxHT5pKRx87NYg8VTYqW0hVoi191YmA7TfIVl6n7GX6FgLjgJ5N+bZyhy6kiiAQBAHY5rIpf08vrP2dDonBQc3Fx91Ltz\/9Pj387rQkCu2+H32F9v3oziRn0rxXP+89bHyQwVkaFQJTDPXd679G+1\/e0Hl51NRLaMQnwhHaWPKG8YD1ODDFUKz2LdtvkK28Wmn99SYExgFLm\/Jt5Q5d1FYlLdw2n\/dIiwOTJbDplg63t2nTpkHXsXn\/efz48UmsBTzz589PkrXhyYrMrAfCUpO0mZIi8t623odLP48\/\/ngnHwuPXe3fJWKuisvmzZtJDVGONjmODupt27ZtyBX1IkYC2\/zOp35Cm3beT1v+86SExCze9tOEoBx463m6\/98fpt8c8RadecyRdO1t36Q7V71OOOj76JWnFOkuqndcCAwEXrduXRC5catv3LhxQdrmRpEnCNfkfT0gLwd2bKMjL7yx0+TO1TfTgTe30pgrD+Ye8vXkDTmAfkPmDfM1rjrawTw32THsgx+fdlLHGE19CoFx0ISNoNjKHbqorUrVBGbXrl00e\/bsTnA8vr7I5ALBovr6+jqRbE1xZQAWn3PgjLS2eCqmWA8cQRPtpX11oV3kJ+JcLKwoneCZogfrStUJjHqdWg9YxykR9DgzeQ0FtvnsxRNo38vr6bv\/0Utv7jqVnvuL92jN0lOTpjg+DILcLZ99dpLc8EsrfzhsCEzIRb1pZ2D2v7ye9r2yYRB5YXsDiRl5\/BSvnpg0gom\/+8wbpoY1wDrDeZRALtWPJXUNUgPW6bnRys7JvHPYpb7tN8hW7tJHjHWEwDhoxaZ8W7lDF7VV0QkMey4QvhsPvAII2a\/GgVEj7prcv0hYqIfT50UEOZf27t3b8WbwD\/\/AwEBCWvRkhmnAmK4ygwQguRvnR8G7nBSRv94QW+L8889PFi4XAuPqDUm7Wq16aYAbjxMemCwCg\/dMXp+8hsK2+bObzuqQmC3v\/y51LRlNs\/9s\/CAS85EX76Xunr+jM485wntm5rxy+6jv4oERAvNLpHd89a8yvSy28rw6qypvmNoP5h3Hk0HuNzwgKHpEbyYwiP2i50ZT15S8Yw5V3\/YbZCsPJVfodoXAOCBsU75afubqTyUtPnPhpuS\/+r\/5b0XLfbSvDlklMGpmV\/w\/B7VDfSYk+H\/8+GNradKkSbRz584krL7qkVDrIycRHi5HOHD8MCN1AR4OEY6\/gcCgHhYNLDRp2zKqNyUt86spJD8HqeMgfbxNBTnSrk6m5VNiDLOu65q8QiqB0beQ1GB9GCOwAU5qVmsHcx1URbVNJjFjb3qSbpr3UXr3sA+Sm0iJvf54FyFY2QkfPowWffvEvN1EWX+4EZjXLzzEWQ9Hr\/7fIXWBjwgIAAAT\/klEQVTLvo8Gs9rQ+7QFhPSRN0yfg7qn0xQUE9F9VQKjR+YWAuNsZsErCoFxgJh\/BHSXMP6Nx\/c+tINI3qqoBObBBx9M2sUEVQkJzmP09PQM6jPNC4Mf4SwCA+KDxGzs4eH\/NyVAZK+PTi5czueYzsCwzKZgVaavQdRbsGBBkn3a5jbW86kAA927ZdtCUqMA5+k7yxhU20W9ESs+T++\/uZVeonvovRlj6cd7R1L\/O\/toz549tO4bd9HKzywZdOXam6FF2JDtw6SsyE3bQsrysBzYsZVQ\/rH\/J7++niryhmURGHyk4SMCaxY+xtT0BEJgfGk5bDtCYBzwtS10tnKHLmqr4kpgOKS3KqgaAp+3ZTh3UtoWEpeDuOBhIpOWwTktmaFpy0ZdEJGVmsmY7rExERjTYqqTCJUU6Zm4TWeJ8hAYyKqOyTeBUfUGT8wPNv4JnXbuq3ToUV2dojvv\/ApNPnYZ\/dF902qzxyo7Dj1vm0Zgqj4Dk\/cQr+2cm+kjAvaUtoUEAsPrjrq1JB6YKmdhub6EwDjgZ1vobOUOXdRWxXULifeN+eAbzo8gu\/TChQuTQGFpW076FhJvPeFr54MPPkgSIvJCcsMNN9Att9yS3DxiN23aVkreQ7zqONUxcDqAtP34tC0knWCkeXDUszZoK+sMjOqBCbGFpBrZ2s\/cRidMWEZjb36SDh1zkMTAjvWgd7UZZgUdh563TSMwgBxeFpBavoUEz8vu9SsSbag3k3ypp4q8YWmHeHFQmLer4YnZvXt3sjUOj7N4YHxpOGw7QmAc8LUtdLZyhy5qq1LmEK++MIwaNYpOPPFEmjBhQuohXr6Wrf+Y85eQegYGoGRdo9b3yPVr1OyBYXBBMpYuXUp333033XPPPZ2r2ln9oA94i0z73i7XqLMO8aZdo4Y8Wf3mMZY023xy3vfouI\/MJyQI\/e7mXvrktCNpzremCIHJA66lbhMJDIbEnhhcncYz4hNTaNTUSz0iU19TaV6a+iTy07PtN8hW7keK6lsRAuOAuU35tnKHLqRKpAikXaMOLa7Pa9SmmzZP3\/4\/9OlL9yVf3ExiNnz4Lvrb7gVyBsaTcptKYDwNP5pm1BuUEErf\/o1G0BKC2H6DbOUluq71VSEwDvDblG8rd+hCqkSMgC9viOsQffaXZpvfmvNT+vWP\/yod3CL4Gr3589OSGDHjJx9Gn3voGFdRG10v9LwVAtNo82iU8DZbtpU3arCKsEJgHDRnU76t3KELqSIIBEEgzTYRmEw959B\/xUTavHsETb1sfpCzDkEGV7LR0PNWCExJBcnrzgjYbNlW7txRZBWFwDgoxKZ8W7lDF1JFEAiCQJpt8mFN7hTnHeCJwYHe8f90MDxA25\/Q81YITNstKJ7x2WzZVh7PSPJJIgTGAS+b8m3lDl3UVqXqVAK+ciG55Doq05ep\/bTAelnnZNR21EPGJoWr8W1wQwKBBHFriR81l5MpvYGpzTy2madubQbrsePQ4xUCk62sqtaevCbjEmdKb1NdA\/hGJuZu3hxpmNfLly8f1Dyf2cmKCm6zZVt5XoxiqS8ExkETNuXbyh26qK1KVYsILwq+ciG5RNot05caR8amHBwSxM0rPSqwGrgOif3mzZtHixcvTg2KpxMYtb6+qLqmGoBtjr3tSFo69Z\/plDEHo0SnPU22Y5uOypK7Iu0LgWkmgSmia14DOMq3etsS7WVF71XXCT2Vih6tOO1wP+burj\/tTSJrn\/lbR3YibPNY2jq3hcA4WKtN+bZyhy5qq1LmGjWErjMXkr4wcCA4DqpXJu+SS3JG9I94LSA7c+bMIY55Y1JmVlh7U\/4peGCyCAz6cLmpBNs8efkMenPXY1YS02Q7LjKBQo9XCEw+AhMqDxukSAtngA8L1dMJL6kalBMxYZBctr+\/P0lvokYg59FlrQF5c6TZCEzahwtsecG\/9tEzr+2klRvfoPGjR9D00z5G1\/zxL2M8hcz7VWT++XhHCIwDiraFzlbu0EVtVVwD2cWUC4mJk0pgTNmmy+RdytqiUomKy40hbsu0haTGpVCDAZq2kPT38dWHh4PxpXkZjvvGS7T2079mtbEm27F1cIYKoccrBMadwFSRh42Tt6oxqJB4FmkEMKfZ28EBNTmo3apVq6i3tzcJuGnyoprWgKI50vQtJH0LKi3ApWrL\/e\/sT0jMbY8fPMuG5Kwvrv4ybV\/\/9SLTJOp3hMA4qMe20KnlJz+5y9riS2cdkdRR6\/Lf9L+bGivzvtoPf5nwpGxKLiQmMPpesZ7rqEzeJdes2C4kAvKa0hfoRMxlC0nNoO1Cnmy2q9pXnrpWI29AhdDjrZrAzJw5M0EdQRrx6P\/mv9VZrpqFaxqTsnnYslJ6gLio3hmcM9MJDDwwnGiW8yWpudGy1oC8OdJMHhg1a31aipE0WwaZAZEBoXnny2c1YFbmE1EIjANetoXOVu7QRW1VXBeRGHMhATTsLaflOiqTd0klc1nKURcvfSHEoqd6a0zbUurfsggMZEBfqh6KEJju9V+gN37xM1p97pohw2qyHReZQKHHWzWBKYJBne9UtfZkERgQBswjeFjgjYGn2SeByZsjTScwpuzZpgSzNlvuOuH3aOt\/fb9OdQfpWwiMA6w247CVO3RRWxXXLaTYciHpZMCU66hM3iXXQ7xZJEIt00kWKzxrC8l0Bkb1wLh4f9g2H33jPbrhh7+ge045nE494lCjvTXZjotMoNDjFQKTrZWq1x7TFhIIQ1dXV7INy1tLeQmMOs\/L5kizeWBctpBMqIe29SLzz8c7QmAcULQp31bu0EVtVcoc4q0zF5LJm6HnOiqTd8l0BgZK0sOQ25IuulyjTjvEm3aNGl4d12zVqm3O\/MEeGtj\/Ad37yT302bXnJbeScDuJnybbcZEJFHq8QmDcCQy2ZPIc4s2z9kAKl5xk119\/PT311FN06623JrcF+QyMbQtJXwPK5EjLukaNcaR9MNls2VZeZP7E8I4QGAct2JRvK3foQqo0GIG0a9Qhh5TnGrXp9sEbeweGkJgm2DH\/yGUl+VR\/3LLqhR6vEJiQM6B42643DPP0UNUakHWNOuuWUWhbz4OVz7pCYBzQtCnfVu7QhVRpMAI2L4zvobl6X9BvHtvMU9f3mFzaw48Efw3jSqt6uFF9X\/2Byso+HHq8QmBctBq+jkpo0ZtLcLm8UlWxBmRtV9ts2Vaed7yx1BcC46AJm\/Jt5Q5dSBVBIAgCJtt8YdcBwnYSzsP8yvubCQd7rzvtJrry7Kso5lgRakyNtDNF+m2vtECNecldEeUIgSmCmrxTBAHbb5CtvEifMbwjBMZBCzblo1weQSBWBEykZNnW\/bTsJ\/uT+DCb3\/g2Ldh4E721fDe9+p3\/jnIYeiTitGvpEN7VU2Ob12WBqIrAlJVT3m8HArKF1A49eh+Fr4WuigUt7+BjlAljELncNVkUK3hh8MATU4VHwn1EQ2umXZVXb2Wpb4HEIH4I548yRUnmD49169aVES31XUR5HTduXJC2fTYqcvpDMzYszznnnM7gYvauFtWAeGAckBMC4wCS5ypFf5Q9izGkuRjlKioTbiSd+9y7ybVqkBhfdh5CB64emOG2heQD66L246PvPG00Qc5YZYx5buexAb2uEBgH9HwpP0bjjlEm8cA4GKVSpYwO+TzMFyeOoMVnH9\/4MzD6mZesrSZf8zpNW2X0ks8CytUWOcvhp74dK5ahbd0fgvlaEgLjgJcv5cdo3DHKJATGwSg9ERg0wyRm7+LLactjq\/N1XmFtl7MtJg8M4uksWbJkSLZwX\/NaCEw1RhDrWiUEphr9m3oRAuOAva+FLsYJGKNMQmAcjNIjgUFTfKg3K1JvPqnC1DbFgdHDrfN5GVy1xmNKoom\/+5rXQmDC6FpvNda1SghMNfoXAlMQZ18LXYwTMEaZhMDkM1RfOvRl5\/mkr6926PH60ktohEROfwjHimVoW\/eHYL6WxAPjgJcv5cdo3DHKJATGwSg9e2Cq8EjkG1X42r7mtXhgwusq5jVBPDDV6F88MAVx9rXQxUgWYpQp5sUqRrx8yeTLzgtOs8pfCz1eX3oJDYzI6Q\/hWLEMbev+EMzXknhgHPDypfwYjTtGmYTAOBileGDygWSo7WteiwemtCqcGoh1rRIPjJP6glQSAuMAKy90+gTCv\/FMnDgx+a+tXF8w875va79IOQIdcYCjIu\/nGX+e9oEVAoyp2Op91fFvxsukuzrkAT6MVdn+VVtwmBaNrxKawIRu35cCRE5fSIY\/GF5U0qboOO\/4hMA4IMYROx2qShVBoNEItDFaZ5pCZF432lRF+JwItHFuC4HJaQRSXRAQBAQBQUAQEATqR0AITP06EAkEAUFAEBAEBAFBICcCQmByAibVBQFBQBAQBAQBQaB+BITA1K8DkUAQEAQEAUFAEBAEciIgBCYnYFJdEBAEBAFBQBAQBOpHQAhM\/ToQCQQBQUAQEAQEAUEgJwJCYHICVqT6li1baOHChXTHHXfQ6NGjkyaQWbenp6fT3KxZs+iaa64p0nzhd0xyqYnwzjvvPFq0aBGNHDmycB9lXuTEfdxGWlK+Mn24vguskNV4YGCA6tBVmpwxYeSK5XCoF\/vc0nUQsx3FOveahGFb55wQmMCa5ck3ZswY6u3t7RAYLBhdXV00bdq0wBKYmzfJtW\/fPpo\/fz5Nnjw5kQsy4qmaWKFPlmX69Ol0xhln1IIRd8qkDjicdNJJgzCqU7CYMKoTh9j6jn1u6XjFbEexzr0mYRjb\/PApjxAYn2hqbT3\/\/PM0d+5cuvbaa2n16tUdD0zdC0aaXFgsrr766kTeY489llAPJEYlXgHhGtS0LktV\/Zr60b+m4T3r6+ur1TsFOWPCqE79xNR3E+aWjlfMdhTr3GsShjHND9+yCIHxjaihPX0Sqts0qI6v+jpIgi6X6d\/z5s2jxYsXJ4Smykd1G6PfOrez8KO0cuXKDmGpk9ipOogJoyptowl9xTy3dPxitqNY516TMGzCfCkqoxCYosjleM9GDODlwNmKqs+b6HLpP8wor4vAqLLgDA62tsaOHVvLdpbucYmFwMSEUY7pMCyqxjy3dAXEbEexzr0mYdjmCScExpN2sQjMmDEjaU33qJgO9Olf0qGIQh65bETLE1RDmuEttTVr1iRlpsO6dZKGpnwF1olRKNuIud0mzC0dv9jnmokYxOj9tNmlzEUbQn7KhcD4wTGzFRcCo99SqkAsMm1txXIGxraQVYEP99GUfXidaFWJkfQ1GIEmza2Y5pouS1PmXswYtnluCoGpQLsmd7L6VVHXbR9drphuIcF1vHXr1mTLSJerApUN6iLWmxAxYVS1TmLvL+a5pWMXsx3FOveahGHsc6WMfEJgyqDn+K4tDkxdB1Rjj1WhxqaoO\/ZKrLEoYsLIcToMi2qxzy1dCTHbUaxzr0kYtnXSCYFpq2ZlXIKAICAICAKCQIsREALTYuXK0AQBQUAQEAQEgbYiIASmrZqVcQkCgoAgIAgIAi1GQAhMi5UrQxMEBAFBQBAQBNqKgBCYtmpWxiUICAKCgCAgCLQYASEwLVZuiKHpaRAQHfe+++7rpBrAlcylS5cO+hvkUGOU4Fr0ZZddRps3bx4i4gMPPJCZvBHvLliwgC655JIh6Q3UK5eTJk0y9qHeZoJMGzZsqCW6bwjdSJuCQFEEZF4XRU7eqxMBITB1ot+wvnmRu+iiizpZtPV0AyAwPT09pF97NhEYxHhRM03j3VWrVmXmhUIdPKYs3iYCo\/Zhkh\/XR6dMmVJ7xuuGmYKI2yIEZF63SJnDbChCYIaZwssMNy03EkhAV1dXQiqYhOzYsYOWLFnSIQYuBEYlICqxYZlRfvPNN9ONN95Io0ePHjIUG4HBC3rQQIxpxYoVdN111xFyLskjCAw3BGReDzeNt2e8QmDao8vgIzF9qemdcvK1k08+mR599NGON8UHgTFt+aj5aBAQsL+\/P9kS4i0kkwdG\/5uaPiE4iNKBIBAZAjKvI1OIiOOMgBAYZ6ikIhBQCQP+rZ9ZYQIDLwm8JZMnT048My4ExraFpHp60LfuseHtK8iUdgbGlChSb1c0LQgMNwRkXg83jbdjvEJg2qHHWkbB4cfVg7xMYBYtWkTbt28nzrL99ttvE+d\/SjvEqx8IVgfF+ZCmT58+aFsKMvT29iZbSrYtJG4D\/cALw09duahqUZp0KghYEJB5LSbSFASEwDRFUxHLiQVvYGCAQFqwbdTX15f8P86UcNkFF1xAjzzySPJ3JjC8lZNGLGwERiVL6EslOaYtJPYgqaQHf1OT2UUMs4gmCFSKgMzrSuGWzgogIASmAGjD9ZW0H3qVSOgEhr0iOJTLJEcnMOp2kHrDyUZg4PbO44FJIzDigRmuFi3jziLwMq\/FPmJHQAhM7BqKSD7OCtvd3d25xqwfANS9IkwaZsyYQZx120RguN7cuXOHxJBhCNLOwDDp4X189QyMemA3awuJb1FFBLeIIghUgoDM60pglk4CICAEJgCobW5SD3iFsaoHeU0EhokD6pq2kBgvroebRHyuRcXSdAuJF194dxDPBQ9i0KQd4mUSxVemMR65hdRmi5WxuSAg89oFJakTGwJCYGLTiMiTioAtDkwR6NTbURIHpgiC8o4gUA4Bmdfl8BvObwuBGc7ab+DYsyLxFhmOROItgpq8Iwj4RUDmtV88h0trQmCGi6ZbMs6sXEh5hyi5kPIiJvUFgTAIyLwOg2vbWxUC03YNy\/gEAUFAEBAEBIEWIiAEpoVKlSEJAoKAICAICAJtR0AITNs1LOMTBAQBQUAQEARaiIAQmBYqVYYkCAgCgoAgIAi0HYH\/A87jxejyoU30AAAAAElFTkSuQmCC","height":337,"width":560}}
%---
