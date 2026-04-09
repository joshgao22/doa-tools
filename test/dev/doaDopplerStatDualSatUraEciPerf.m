% DOADOPPLERSTATDUALSATURAECIPERF
% Static single-frame performance sweep for one fixed dual-satellite pair.
% The script follows the same reference-satellite parameterization used by
% doaDopplerStatDualSatUraEci and compares
%   1) SS-SF-DoA
%   2) MS-SF-DoA
%   3) SS-SF-Static
%   4) MS-SF-Static
% together with the static multi-satellite weight sweep
%   alphaSat2 in {0, 0.25, 0.5, 1}.
% The performance curves are reported by angle RMSE, reference-Doppler
% RMSE, and estimator resolve rate.
clear(); close all; clc;

%% Parameters
numUsr = 1;
numSym = 32;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 18e9;
wavelen = 299792458 / carrierFreq;
rng(253);

elemSpace = wavelen / 2;
numElem = [4 4];
pwrSource = 1;
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);
truthLatlon = usrLla(1:2, 1);

utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(localResolveTlePath("starlink_pair_4154_1165_20260318_170800.tle"));
arrUpa = createUpa(numElem, elemSpace);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [truthLatlon(1) - searchMarginDeg, truthLatlon(1) + searchMarginDeg; ...
               truthLatlon(2) - searchMarginDeg, truthLatlon(2) + searchMarginDeg];

fdRange = [-2e5, 0];
weightSweepAlpha = [0, 0.25, 0.5, 1];
optVerbose = false;

snrDb = -20:3:10;
numParam = numel(snrDb);
numRepeat = 200;

%% Select satellites and build the reference scene
[~, satAccess] = findVisibleSatFromTle(utc, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccess, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

scene = genMultiSatScene(utc, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal);
linkParam = getLinkParam(scene, wavelen);
truthDopplerState = buildReferenceDopplerState( ...
  scene, scene.satPosEci, scene.satVelEci, scene.usrPosEci, scene.usrVelEci, wavelen);
steeringInfo = getSceneSteering(scene, wavelen);
[refState, refSatIdxLocal] = resolveReferenceSatState(scene, scene.satPosEci, scene.satVelEci);

if scene.numSat ~= 2
  error('doaDopplerStatDualSatUraEciPerf:InvalidNumSat', ...
    'This script expects exactly two selected satellites.');
end
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

truth = struct();
truth.utc = scene.utc;
truth.latlonTrueDeg = truthLatlon;
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.usrElevationDeg = reshape(scene.accessInfo.usrElevationDeg(:, 1), 1, []);
truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
truth.deltaFdTrueHz = reshape(truthDopplerState.deltaFdRefFrame, [], 1);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);

viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], gridSize, searchRange, E);
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(scene, [], gridSize, searchRange, E);

truthLocalDoaSingle = reshape(sceneRefOnly.localDoa(:, 1), 2, 1);
truthLocalDoaCell = localBuildLocalDoaCell(scene);
localDoaJacCell = localBuildLatlonToLocalDoaJac( ...
  scene, utc, tle, usrLla, arrUpa, satIdx, refSatIdxGlobal, 15, 55);

%% Pilot waveform and snapshot model
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
numSnap = length(pilotWave);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);

%% Common estimator options
caseName = ["SS-SF-DoA", "MS-SF-DoA", "SS-SF-Static", "MS-SF-Static", ...
  "MS-SF-Static-W0.00", "MS-SF-Static-W0.25", "MS-SF-Static-W0.50", "MS-SF-Static-W1.00"];
numCase = numel(caseName);
idxSsDoa = 1;
idxMsDoa = 2;
idxSsStat = 3;
idxMsStat = 4;
idxWeightStart = 5;

doAOnlyOpt = struct();
doAOnlyOpt.useLogObjective = true;

staticOptBase = struct();
staticOptBase.useLogObjective = true;

%% Selected satellites
fprintf('\n========== Selected satellites ==========%s', newline);
disp(table((1:scene.numSat).', truth.selectedSatIdxGlobal(:), truth.usrElevationDeg(:), ...
  truth.fdSatTrueHz(:), truth.deltaFdTrueHz(:), ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg', 'fdSatTrueHz', 'deltaFdTrueHz'}));

fprintf('\n========== Truth ==========%s', newline);
disp(table(truth.latlonTrueDeg(1), truth.latlonTrueDeg(2), truth.fdRefTrueHz, ...
  truth.refSatIdxLocal, truth.refSatIdxGlobal, otherSatIdxLocal, otherSatIdxGlobal, ...
  'VariableNames', {'latTrueDeg', 'lonTrueDeg', 'fdRefTrueHz', ...
  'refSatIdxLocal', 'refSatIdxGlobal', 'otherSatIdxLocal', 'otherSatIdxGlobal'}));

%% Monte Carlo loop
angleErrDeg = nan(numRepeat, numParam, numCase);
fdErrHz = nan(numRepeat, numParam, numCase);
isResolved = false(numRepeat, numParam, numCase);
progressbar('reset', numParam);

for iSnr = 1:numParam
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));

  angleErrCur = nan(numRepeat, numCase);
  fdErrCur = nan(numRepeat, numCase);
  isResolvedCur = false(numRepeat, numCase);

  parfor iRepeat = 1:numRepeat
    [rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
      steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, ...
      pwrNoise, pathGain, snapOpt);

    rxSigRefOnly = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');
    rxSigOtherOnly = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');

    viewRef = viewRefTemplate;
    viewRef.rxSigSf = rxSigRefOnly;
    viewOther = viewOtherTemplate;
    viewOther.rxSigSf = rxSigOtherOnly;
    viewMs = viewMsTemplate;
    viewMs.rxSigSf = rxSig;

    caseRefDoa = runDoaOnlyCase(caseName(idxSsDoa), "single", ...
      viewRef, wavelen, pilotWave, optVerbose, doAOnlyOpt);
    caseMsDoa = runDoaOnlyCase(caseName(idxMsDoa), "multi", ...
      viewMs, wavelen, pilotWave, optVerbose, doAOnlyOpt);

    staticRefOpt = staticOptBase;
    staticRefOpt.initDoaParam = caseRefDoa.estResult.doaParamEst(:);

    staticMsOpt = staticOptBase;
    staticMsOpt.initDoaParam = caseMsDoa.estResult.doaParamEst(:);
    staticMsOpt.initDoaHalfWidth = [0.01; 0.01];

    caseSsStat = runStaticDoaDopplerCase(caseName(idxSsStat), "single", ...
      viewRef, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticRefOpt);
    caseMsStat = runStaticDoaDopplerCase(caseName(idxMsStat), "multi", ...
      viewMs, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticMsOpt);

    weightCase = repmat(caseMsStat, 1, numel(weightSweepAlpha));
    for iAlpha = 1:numel(weightSweepAlpha)
      alpha = weightSweepAlpha(iAlpha);
      if abs(alpha - 1) < eps
        weightCase(iAlpha).displayName = caseName(idxWeightStart + iAlpha - 1);
        continue;
      end

      currentOpt = staticMsOpt;
      currentOpt.satWeight = [1; alpha];
      weightCase(iAlpha) = runStaticDoaDopplerCase( ...
        caseName(idxWeightStart + iAlpha - 1), "multi", ...
        viewMs, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, currentOpt);
    end

    caseList = [caseRefDoa, caseMsDoa, caseSsStat, caseMsStat, weightCase];
    truthFdList = [NaN, NaN, truth.fdSatTrueHz(refSatIdxLocal), truth.fdRefTrueHz, ...
      repmat(truth.fdRefTrueHz, 1, numel(weightSweepAlpha))];

    angleVec = nan(1, numCase);
    fdVec = nan(1, numCase);
    resolvedVec = false(1, numCase);
    for iCase = 1:numCase
      [angleVec(iCase), fdVec(iCase), resolvedVec(iCase)] = localExtractCaseMetric( ...
        caseList(iCase), truth.latlonTrueDeg, truthFdList(iCase));
    end

    angleErrCur(iRepeat, :) = angleVec;
    fdErrCur(iRepeat, :) = fdVec;
    isResolvedCur(iRepeat, :) = resolvedVec;
  end

  angleErrDeg(:, iSnr, :) = angleErrCur;
  fdErrHz(:, iSnr, :) = fdErrCur;
  isResolved(:, iSnr, :) = isResolvedCur;
  progressbar('advance');
end
progressbar('end');

%% Monte Carlo summaries
rmseAngleDeg = nan(numParam, numCase);
rmseFdHz = nan(numParam, numCase);
resolveRate = nan(numParam, numCase);
p95AngleDeg = nan(numParam, numCase);
p95FdHz = nan(numParam, numCase);

for iCase = 1:numCase
  for iSnr = 1:numParam
    failMaskAngle = ~isResolved(:, iSnr, iCase) | ~isfinite(angleErrDeg(:, iSnr, iCase));
    angleStat = summarizeMonteCarloStat(angleErrDeg(:, iSnr, iCase), failMaskAngle);
    rmseAngleDeg(iSnr, iCase) = angleStat.rmse;
    p95AngleDeg(iSnr, iCase) = angleStat.p95;
    resolveRate(iSnr, iCase) = 1 - angleStat.failRate;

    if iCase >= idxSsStat
      failMaskFd = ~isResolved(:, iSnr, iCase) | ~isfinite(fdErrHz(:, iSnr, iCase));
      fdStat = summarizeMonteCarloStat(fdErrHz(:, iSnr, iCase), failMaskFd);
      rmseFdHz(iSnr, iCase) = fdStat.rmse;
      p95FdHz(iSnr, iCase) = fdStat.p95;
    end
  end
end

weightSummary = localBuildWeightSummaryTable(snrDb, weightSweepAlpha, ...
  rmseAngleDeg(:, idxWeightStart:end), rmseFdHz(:, idxWeightStart:end), ...
  resolveRate(:, idxWeightStart:end));

snrSelectIdx = numParam;
caseSummaryTable = localBuildCaseSummaryTable(caseName, snrDb(snrSelectIdx), ...
  rmseAngleDeg(snrSelectIdx, :), rmseFdHz(snrSelectIdx, :), ...
  resolveRate(snrSelectIdx, :), p95AngleDeg(snrSelectIdx, :), p95FdHz(snrSelectIdx, :));

%% CRB curves
crbDoaSingleOpt = struct();
crbDoaSingleOpt.localDoaJac = localDoaJacCell{refSatIdxLocal};
crbDoaSingleOpt.paramName = ["lat", "lon"];

crbDoaJointOpt = struct();
crbDoaJointOpt.localDoaJac = localDoaJacCell;
crbDoaJointOpt.paramName = ["lat", "lon"];

crbDdOpt = struct();
crbDdOpt.doaType = 'latlon';

crbDoaOnlySingleDeg = zeros(numParam, 1);
crbDoaOnlyJointDeg = zeros(numParam, 1);
crbDdSingleDeg = zeros(numParam, 1);
crbDdJointDeg = zeros(numParam, 1);
crbDdSingleFd = zeros(numParam, 1);
crbDdJointFd = zeros(numParam, 1);

for iSnr = 1:numParam
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));

  [crbDoaSingle, ~] = crbDetDoa(sceneRefOnly.array{1}, wavelen, ...
    truthLocalDoaSingle, pwrSource, pwrNoise, numSnap, crbDoaSingleOpt);
  [crbDoaJoint, ~] = crbDetDoa(scene.array, wavelen, truthLocalDoaCell, ...
    pwrSource, pwrNoise, numSnap, crbDoaJointOpt);

  [crbDdSingle, ~] = crbPilotSfDoaDoppler(sceneRefOnly, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);
  [crbDdJoint, ~] = crbPilotSfDoaDoppler(scene, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);

  crbDoaOnlySingleDeg(iSnr) = projectCrbToAngleMetric(crbDoaSingle, truthLatlon, 'latlon');
  crbDoaOnlyJointDeg(iSnr) = projectCrbToAngleMetric(crbDoaJoint, truthLatlon, 'latlon');
  crbDdSingleDeg(iSnr) = projectCrbToAngleMetric(crbDdSingle(1:2, 1:2), truthLatlon, 'latlon');
  crbDdJointDeg(iSnr) = projectCrbToAngleMetric(crbDdJoint(1:2, 1:2), truthLatlon, 'latlon');
  crbDdSingleFd(iSnr) = sqrt(crbDdSingle(3, 3));
  crbDdJointFd(iSnr) = sqrt(crbDdJoint(3, 3));
end

%% Summary tables
fprintf('\n========== Weight-sweep best summary ==========%s', newline);
disp(weightSummary);

fprintf('\n========== Case summary at highest SNR ==========%s', newline);
disp(caseSummaryTable);

%% Plots
figure();
subplot(1, 2, 1);
semilogy(snrDb, rmseAngleDeg(:, idxSsDoa), '-o', ...
  snrDb, rmseAngleDeg(:, idxMsDoa), '-s', ...
  snrDb, rmseAngleDeg(:, idxSsStat), '--o', ...
  snrDb, rmseAngleDeg(:, idxMsStat), '--s', ...
  snrDb, rmseAngleDeg(:, idxWeightStart + 1), '--d', ...
  snrDb, crbDoaOnlySingleDeg, '-.', ...
  snrDb, crbDoaOnlyJointDeg, ':', ...
  snrDb, crbDdSingleDeg, '-.x', ...
  snrDb, crbDdJointDeg, ':x');
grid on;
xlabel('SNR (dB)');
ylabel('Angle RMSE / CRB (deg)');
title(sprintf('Static SF Global Angle Error (%d snaps)', numSnap));
legend({ ...
  'SS-SF-DoA', ...
  'MS-SF-DoA', ...
  'SS-SF-Static', ...
  'MS-SF-Static', ...
  'MS-SF-Static-W0.25', ...
  'SS-SF-DoA CRB', ...
  'MS-SF-DoA CRB', ...
  'SS-SF-Static CRB', ...
  'MS-SF-Static CRB'}, ...
  'Location', 'southwest');

subplot(1, 2, 2);
semilogy(snrDb, rmseFdHz(:, idxSsStat), '--o', ...
  snrDb, rmseFdHz(:, idxMsStat), '--s', ...
  snrDb, rmseFdHz(:, idxWeightStart + 1), '--d', ...
  snrDb, crbDdSingleFd, '-.x', ...
  snrDb, crbDdJointFd, ':x');
grid on;
xlabel('SNR (dB)');
ylabel('Reference Doppler RMSE / CRB (Hz)');
title('Static SF Reference Doppler Error');
legend({ ...
  'SS-SF-Static', ...
  'MS-SF-Static', ...
  'MS-SF-Static-W0.25', ...
  'SS-SF-Static CRB', ...
  'MS-SF-Static CRB'}, ...
  'Location', 'southwest');

figure();
subplot(1, 2, 1);
plot(snrDb, resolveRate(:, idxSsStat), '-o', ...
  snrDb, resolveRate(:, idxMsStat), '-s', ...
  snrDb, resolveRate(:, idxWeightStart + 1), '-d');
grid on;
ylim([0, 1.05]);
xlabel('SNR (dB)');
ylabel('Resolve rate');
title('Static estimator resolve rate');
legend({'SS-SF-Static', 'MS-SF-Static', 'MS-SF-Static-W0.25'}, ...
  'Location', 'southeast');

subplot(1, 2, 2);
plot(snrDb, rmseAngleDeg(:, idxWeightStart:end), '-o');
grid on;
xlabel('SNR (dB)');
ylabel('Angle RMSE (deg)');
title('Static multi-satellite weight sweep');
legend(compose('alpha = %.2f', weightSweepAlpha), 'Location', 'southwest');

%% Optional snapshot save
saveExpSnapshot("doaDopplerStatDualSatUraEciPerf");

%% Local function
function [angleErrDeg, fdErrHz, isResolved] = localExtractCaseMetric(caseInfo, truthLatlon, truthFdHz)
%LOCALEXTRACTCASEMETRIC Extract one compact metric tuple from a case result.

estResult = getDoaDopplerFieldOrDefault(caseInfo, 'estResult', struct());
latlonEst = getDoaDopplerLatlonEst(estResult);
[angleErrDeg, ~, ~] = calcLatlonAngleError(latlonEst, truthLatlon);

fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
if isfinite(truthFdHz) && isfinite(fdRefEst)
  fdErrHz = fdRefEst - truthFdHz;
else
  fdErrHz = NaN;
end

isResolved = logical(getDoaDopplerFieldOrDefault(estResult, 'isResolved', false));
end

function summaryTable = localBuildWeightSummaryTable(snrDb, weightAlpha, rmseAngleDeg, rmseFdHz, resolveRate)
%LOCALBUILDWEIGHTSUMMARYTABLE Build one compact best-weight table per SNR.

numSnr = numel(snrDb);
bestAlphaAngle = nan(numSnr, 1);
bestAngleRmseDeg = nan(numSnr, 1);
bestAlphaFd = nan(numSnr, 1);
bestFdRmseHz = nan(numSnr, 1);
bestResolveRate = nan(numSnr, 1);

for iSnr = 1:numSnr
  [bestAngleRmseDeg(iSnr), idxAngle] = min(rmseAngleDeg(iSnr, :), [], 'omitnan');
  if ~isempty(idxAngle) && isfinite(bestAngleRmseDeg(iSnr))
    bestAlphaAngle(iSnr) = weightAlpha(idxAngle);
    bestResolveRate(iSnr) = resolveRate(iSnr, idxAngle);
  end

  [bestFdRmseHz(iSnr), idxFd] = min(rmseFdHz(iSnr, :), [], 'omitnan');
  if ~isempty(idxFd) && isfinite(bestFdRmseHz(iSnr))
    bestAlphaFd(iSnr) = weightAlpha(idxFd);
  end
end

summaryTable = table(snrDb(:), bestAlphaAngle, bestAngleRmseDeg, ...
  bestAlphaFd, bestFdRmseHz, bestResolveRate, ...
  'VariableNames', {'snrDb', 'bestAlphaAngle', 'bestAngleRmseDeg', ...
  'bestAlphaFd', 'bestFdRmseHz', 'bestResolveRate'});
end

function summaryTable = localBuildCaseSummaryTable(caseName, snrDb, rmseAngleDeg, rmseFdHz, resolveRate, p95AngleDeg, p95FdHz)
%LOCALBUILDCASESUMMARYTABLE Build one case summary row per estimator.

summaryTable = table(caseName(:), repmat(snrDb, numel(caseName), 1), ...
  reshape(rmseAngleDeg, [], 1), reshape(rmseFdHz, [], 1), ...
  reshape(resolveRate, [], 1), reshape(p95AngleDeg, [], 1), ...
  reshape(p95FdHz, [], 1), ...
  'VariableNames', {'displayName', 'snrDb', 'angleRmseDeg', ...
  'fdRefRmseHz', 'resolveRate', 'angleP95Deg', 'fdRefP95Hz'});
end

function doaGridCell = localBuildDoaGridCell(sceneRef, gridSize, searchRange, E)
%LOCALBUILDDOAGRIDCELL Build one lat-lon search grid per selected satellite.

numSat = sceneRef.numSat;
doaGridCell = cell(1, numSat);
for iSat = 1:numSat
  doaGridCell{iSat} = genDoaGrid("latlon", 2, gridSize, searchRange, ...
    'eci', datevec(sceneRef.utc), sceneRef.satPosEci(:, iSat), ...
    sceneRef.rotMat{iSat}, E);
end
end

function localDoaCell = localBuildLocalDoaCell(scene)
%LOCALBUILDLOCALDOACELL Convert scene.localDoa into a 1xL cell array.

numSat = numel(scene.array);
localDoaCell = cell(1, numSat);
for iSat = 1:numSat
  localDoaCell{iSat} = reshape(scene.localDoa(:, iSat), 2, 1);
end
end

function jacCell = localBuildLatlonToLocalDoaJac(scene, utc, tle, usrLla, arrUpa, ...
  satIdxSel, refSatIdxGlobal, minUsrElevationDeg, maxSatOffAxisDeg)
%LOCALBUILDLATLONTOLOCALDOAJAC Build numerical Jacobians for multi-array CRB.

numSat = numel(scene.array);
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
    minUsrElevationDeg, maxSatOffAxisDeg, "satellite", refSatIdxGlobal);
  deltaLocalDoa = (scenePert.localDoa - baseLocalDoa) / deltaDeg;

  for iSat = 1:numSat
    jacCell{iSat}(:, iParam) = reshape(deltaLocalDoa(:, iSat), 2, 1);
  end
end
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

error('doaDopplerStatDualSatUraEciPerf:TleFileNotFound', ...
  'Unable to locate the TLE file "%s".', fileName);
end
