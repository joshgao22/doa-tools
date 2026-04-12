% DOADOPPLERSTATDUALSATURAECIPERF
% Static single-frame performance sweep for one fixed dual-satellite pair.
% The script now reuses the shared dev/common helpers so that the static
% performance path stays aligned with the main dev scripts while keeping
% Monte Carlo bookkeeping local to this perf entry.
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
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa(numElem, elemSpace);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [truthLatlon(1) - searchMarginDeg, truthLatlon(1) + searchMarginDeg; ...
               truthLatlon(2) - searchMarginDeg, truthLatlon(2) + searchMarginDeg];

fdRange = [-2e5, 0];
weightSweepAlpha = [0; 0.25; 0.5; 1];
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
fdRange = expandRangeToTruth(fdRange, [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);

viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], gridSize, searchRange, E);
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(scene, [], gridSize, searchRange, E);

truthLocalDoaSingle = reshape(sceneRefOnly.localDoa(:, 1), 2, 1);
truthLocalDoaCell = localBuildLocalDoaCell(scene);
% Note: do not build an array-only DoA CRB here. For this script, the
% meaningful theoretical baseline is the pilot/static DoA-Doppler CRB,
% because the estimated global DoA is coupled with reference-link Doppler
% in the actual signal model.

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

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

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

    viewRef = viewRefTemplate;
    viewRef.rxSigSf = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');

    viewOther = viewOtherTemplate;
    viewOther.rxSigSf = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');

    viewMs = viewMsTemplate;
    viewMs.rxSigSf = rxSig;

    caseBundle = buildDoaDopplerStaticTransitionBundle( ...
      viewRef, viewOther, viewMs, wavelen, pilotWave, carrierFreq, ...
      waveInfo.sampleRate, fdRange, optVerbose, truth, otherSatIdxGlobal, ...
      doaOnlyOpt, staticOptBase, weightSweepAlpha, [0.01; 0.01]);

    caseList = [ ...
      caseBundle.caseRefDoa, ...
      caseBundle.caseMsDoa, ...
      caseBundle.caseStaticRefOnly, ...
      caseBundle.caseStaticMs, ...
      caseBundle.weightCase];
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
crbDdOpt = struct();
crbDdOpt.doaType = 'latlon';

% Use the pilot/static DoA-Doppler CRB as the meaningful baseline for both
% the DoA-only and static estimators. An array-only DoA CRB from crbDetDoa
% is too optimistic here because it omits the reference-link Doppler
% nuisance that is present in the actual single-frame signal model.
crbDdSingleDeg = zeros(numParam, 1);
crbDdJointDeg = zeros(numParam, 1);
crbDdSingleFd = zeros(numParam, 1);
crbDdJointFd = zeros(numParam, 1);

for iSnr = 1:numParam
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));

  [crbDdSingle, ~] = crbPilotSfDoaDoppler(sceneRefOnly, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);
  [crbDdJoint, ~] = crbPilotSfDoaDoppler(scene, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);

  crbDdSingleDeg(iSnr) = projectCrbToAngleMetric(crbDdSingle(1:2, 1:2), truthLatlon, 'latlon');
  crbDdJointDeg(iSnr) = projectCrbToAngleMetric(crbDdJoint(1:2, 1:2), truthLatlon, 'latlon');
  crbDdSingleFd(iSnr) = sqrt(crbDdSingle(3, 3));
  crbDdJointFd(iSnr) = sqrt(crbDdJoint(3, 3));
end

% Reuse the pilot/static CRB curves when plotting the DoA-only cases as a
% conservative and model-consistent lower bound in the same signal model.
crbDoaOnlySingleDeg = crbDdSingleDeg;
crbDoaOnlyJointDeg = crbDdJointDeg;

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
  'SS/MF-consistent DoA lower bound (single)', ...
  'SS/MF-consistent DoA lower bound (joint)', ...
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

%% Local functions
function [angleErrDeg, fdErrHz, isResolved] = localExtractCaseMetric(caseInfo, truthLatlon, truthFdHz)
%LOCALEXTRACTCASEMETRIC Extract one compact metric triple from one case.

angleErrDeg = NaN;
fdErrHz = NaN;
isResolved = false;

if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
latlonEst = getDoaDopplerLatlonEst(estResult);
if isfinite(latlonEst(1)) && isfinite(latlonEst(2))
  angleErrDeg = calcLatlonAngleError(latlonEst, truthLatlon);
end

fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
if isscalar(fdRefEst) && isfinite(fdRefEst) && isfinite(truthFdHz)
  fdErrHz = abs(fdRefEst - truthFdHz);
end

if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  isResolved = logical(estResult.isResolved);
else
  isResolved = isfinite(angleErrDeg);
end
end


function summaryTable = localBuildWeightSummaryTable(snrDb, weightAlpha, rmseAngleDeg, rmseFdHz, resolveRate)
%LOCALBUILDWEIGHTSUMMARYTABLE Build one compact summary for the sat-weight sweep.

numWeight = numel(weightAlpha);
[bestAngleDeg, bestIdx] = min(rmseAngleDeg, [], 2, 'omitnan');
bestFdHz = nan(size(bestAngleDeg));
bestResolve = nan(size(bestAngleDeg));
for iRow = 1:numel(bestIdx)
  idx = bestIdx(iRow);
  if isfinite(idx) && idx >= 1 && idx <= numWeight
    bestFdHz(iRow) = rmseFdHz(iRow, idx);
    bestResolve(iRow) = resolveRate(iRow, idx);
  end
end
summaryTable = table(snrDb(:), weightAlpha(bestIdx(:)), bestAngleDeg(:), bestFdHz(:), ...
  bestResolve(:), 'VariableNames', {'snrDb', 'bestAlphaSat2', 'bestAngleRmseDeg', ...
  'bestFdRmseHz', 'bestResolveRate'});
end


function summaryTable = localBuildCaseSummaryTable(caseName, snrDb, rmseAngleDeg, rmseFdHz, resolveRate, p95AngleDeg, p95FdHz)
%LOCALBUILDCASESUMMARYTABLE Build one compact case summary row per estimator.

summaryTable = table(caseName(:), repmat(snrDb, numel(caseName), 1), rmseAngleDeg(:), ...
  rmseFdHz(:), resolveRate(:), p95AngleDeg(:), p95FdHz(:), ...
  'VariableNames', {'displayName', 'snrDb', 'angleRmseDeg', 'fdRmseHz', ...
  'resolveRate', 'angleP95Deg', 'fdP95Hz'});
end


function localDoaCell = localBuildLocalDoaCell(scene)
%LOCALBUILDLOCALDOACELL Collect truth local DoA for CRB evaluation.

numSat = scene.numSat;
localDoaCell = cell(1, numSat);
for iSat = 1:numSat
  localDoaCell{iSat} = reshape(scene.localDoa(:, iSat), 2, []);
end
end


function jacCell = localBuildLatlonToLocalDoaJac(scene, utc, tle, usrLla, arrUpa, ...
  satIdx, refSatIdxGlobal, minElevDeg, maxElevDeg)
%LOCALBUILDLATLONTOLOCALDOAJAC Build local Jacobians for static DoA CRB.
% The Jacobian is with respect to global lat/lon in degrees so that it is
% consistent with paramName=["lat","lon"] and projectCrbToAngleMetric(...,'latlon').

numSat = scene.numSat;
jacCell = cell(1, numSat);
truthLatlon = usrLla(1:2, 1);
stepDeg = 1e-4;

for iSat = 1:numSat
  jacMat = zeros(2, 2);
  for iDim = 1:2
    latlonPlus = truthLatlon;
    latlonMinus = truthLatlon;
    latlonPlus(iDim) = latlonPlus(iDim) + stepDeg;
    latlonMinus(iDim) = latlonMinus(iDim) - stepDeg;

    usrLlaPlus = [latlonPlus; usrLla(3, 1)];
    usrLlaMinus = [latlonMinus; usrLla(3, 1)];

    scenePlus = genMultiSatScene(utc, tle, usrLlaPlus, satIdx, [], arrUpa, ...
      minElevDeg, maxElevDeg, "satellite", refSatIdxGlobal);
    sceneMinus = genMultiSatScene(utc, tle, usrLlaMinus, satIdx, [], arrUpa, ...
      minElevDeg, maxElevDeg, "satellite", refSatIdxGlobal);

    doaPlus = scenePlus.localDoa(:, iSat);
    doaMinus = sceneMinus.localDoa(:, iSat);
    jacMat(:, iDim) = (doaPlus - doaMinus) / (2 * stepDeg);
  end
  jacCell{iSat} = jacMat;
end
end
