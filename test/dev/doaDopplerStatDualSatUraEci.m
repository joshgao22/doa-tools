% DOADOPPLERSTATDUALSATURAECI
% Fixed-SNR repeated static SS/MS DoA-Doppler diagnosis for one visible
% dual-satellite pair. This script is a lightweight single-SNR version of
% doaDopplerStatDualSatUraEciPerf: it keeps the same shared SF static flow,
% but replaces the full SNR sweep with repeated trials at one SNR so the
% remaining multi-satellite DoA issue can be inspected with compact
% repeat-level and representative-trial diagnostics.
clear(); close all; clc;

%% Parameters
numUsr = 1;
numSym = 512;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
carrierFreq = 11.7e9;
wavelen = 299792458 / carrierFreq;
baseSeed = 253;
rng(baseSeed);

elemSpace = wavelen / 2;
numElem = [4 4];
pwrSource = 1;
E = referenceEllipsoid('sphere');

snrDb = 10;
numRepeat = 40;
staticMsHalfWidth = [0.002; 0.002];
weightSweepAlpha = [0; 0.25; 0.5; 1];
optVerbose = false;

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

fdRange = [-2e5, 2e5];
pwrNoise = pwrSource / (10^(snrDb / 10));

if numUsr ~= 1
  error('doaDopplerStatDualSatUraEci:OnlySingleUserSupported', ...
    'This script currently supports one user.');
end

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
  error('doaDopplerStatDualSatUraEci:InvalidNumSat', ...
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
truth.localDoaRef = reshape(scene.localDoa(:, refSatIdxLocal), 2, 1);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;
truth.fdRateTrueHzPerSec = NaN;
fdRange = expandRangeToTruth(fdRange, [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);

%% Pilot waveform and snapshot model
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);

viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], gridSize, searchRange, E);
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(scene, [], gridSize, searchRange, E);

%% Common estimator options
caseName = ["SS-SF-DoA", "MS-SF-DoA", "SS-SF-Static", "MS-SF-Static", ...
  compose('MS-SF-Static-W%.2f', weightSweepAlpha.')];
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

%% Selected satellites and truth
fprintf('\n========== Selected satellites ==========%s', newline);
disp(table((1:scene.numSat).', truth.selectedSatIdxGlobal(:), truth.usrElevationDeg(:), ...
  truth.fdSatTrueHz(:), truth.deltaFdTrueHz(:), ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg', 'fdSatTrueHz', 'deltaFdTrueHz'}));

fprintf('\n========== Truth ==========%s', newline);
disp(table(truth.latlonTrueDeg(1), truth.latlonTrueDeg(2), truth.fdRefTrueHz, ...
  truth.refSatIdxLocal, truth.refSatIdxGlobal, otherSatIdxLocal, otherSatIdxGlobal, ...
  'VariableNames', {'latTrueDeg', 'lonTrueDeg', 'fdRefTrueHz', ...
  'refSatIdxLocal', 'refSatIdxGlobal', 'otherSatIdxLocal', 'otherSatIdxGlobal'}));

%% Fixed-SNR repeated scan
angleErrDeg = nan(numRepeat, numCase);
fdErrHz = nan(numRepeat, numCase);
isResolved = false(numRepeat, numCase);
bestWeightIdx = nan(numRepeat, 1);
bestWeightAlpha = nan(numRepeat, 1);
bestWeightAngleErrDeg = nan(numRepeat, 1);
bestWeightFdErrHz = nan(numRepeat, 1);
taskSeedList = nan(numRepeat, 1);

fprintf('Running fixed-SNR static diagnosis: %d repeats at %.1f dB.%s', ...
  numRepeat, snrDb, newline);
progressbar('reset', numRepeat);
progressQueue = parallel.pool.DataQueue;
afterEach(progressQueue, @(~) progressbar('advance'));
parfor iRepeat = 1:numRepeat
  taskSeed = baseSeed + iRepeat - 1;
  repeatOut = localRunStaticRepeatTask(taskSeed, pwrNoise, steeringInfo, pilotWave, ...
    linkParam, carrierFreq, waveInfo.sampleRate, pathGain, snapOpt, ...
    viewRefTemplate, viewOtherTemplate, viewMsTemplate, refSatIdxLocal, ...
    otherSatIdxLocal, wavelen, fdRange, optVerbose, truth, otherSatIdxGlobal, ...
    doaOnlyOpt, staticOptBase, weightSweepAlpha, staticMsHalfWidth, numCase, idxWeightStart);

  angleErrDeg(iRepeat, :) = repeatOut.angleErrDeg;
  fdErrHz(iRepeat, :) = repeatOut.fdErrHz;
  isResolved(iRepeat, :) = repeatOut.isResolved;
  bestWeightIdx(iRepeat) = repeatOut.bestWeightIdx;
  bestWeightAlpha(iRepeat) = repeatOut.bestWeightAlpha;
  bestWeightAngleErrDeg(iRepeat) = repeatOut.bestWeightAngleErrDeg;
  bestWeightFdErrHz(iRepeat) = repeatOut.bestWeightFdErrHz;
  taskSeedList(iRepeat) = taskSeed;
  send(progressQueue, iRepeat);
end
progressbar('end');

%% Aggregate summaries at one SNR
caseSummaryTable = localBuildCaseSummaryTable(caseName, snrDb, angleErrDeg, fdErrHz, isResolved);
gainSummaryTable = localBuildStaticGainSummaryTable(angleErrDeg, fdErrHz, isResolved, ...
  bestWeightAngleErrDeg, bestWeightFdErrHz, idxSsStat, idxMsStat);
weightPrefTable = localBuildWeightPreferenceTable(weightSweepAlpha, bestWeightIdx, ...
  angleErrDeg(:, idxWeightStart:end), fdErrHz(:, idxWeightStart:end), ...
  isResolved(:, idxWeightStart:end));
repeatCompareTable = localBuildRepeatCompareTable((1:numRepeat).', taskSeedList, snrDb, ...
  angleErrDeg, fdErrHz, isResolved, bestWeightAlpha, bestWeightAngleErrDeg, ...
  bestWeightFdErrHz, idxSsStat, idxMsStat);

worstLossMask = isfinite(repeatCompareTable.msMinusSsAngleDeg);
worstLossTable = repeatCompareTable(worstLossMask, :);
worstLossTable = sortrows(worstLossTable, 'msMinusSsAngleDeg', 'descend');
worstLossTable = worstLossTable(1:min(10, height(worstLossTable)), :);

fprintf('\n========== Single-SNR case summary ==========%s', newline);
disp(caseSummaryTable);

fprintf('\n========== MS vs SS static gain summary ==========%s', newline);
disp(gainSummaryTable);

fprintf('\n========== Weight preference summary ==========%s', newline);
disp(weightPrefTable);

fprintf('\n========== Worst MS angle-loss repeats ==========%s', newline);
disp(worstLossTable);

%% Representative repeat diagnostics
repTable = localSelectRepresentativeRepeats(repeatCompareTable);
repCount = height(repTable);
repCaseSummary = repmat(struct('repeatIdx', NaN, 'taskSeed', NaN, 'label', "", 'caseTable', table(), ...
  'ablationTable', table(), 'weightPathDiagTable', table(), 'msSatDiagTable', table(), ...
  'bestWeightSatDiagTable', table()), repCount, 1);

for iRep = 1:repCount
  repIdx = repTable.repeatIdx(iRep);
  repSeed = repTable.taskSeed(iRep);
  repBundle = localRunStaticRepresentativeBundle(repSeed, pwrNoise, steeringInfo, ...
    pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, pathGain, snapOpt, ...
    viewRefTemplate, viewOtherTemplate, viewMsTemplate, refSatIdxLocal, ...
    otherSatIdxLocal, wavelen, fdRange, optVerbose, truth, otherSatIdxGlobal, ...
    doaOnlyOpt, staticOptBase, weightSweepAlpha, staticMsHalfWidth);

  mainCase = [repBundle.caseBundle.caseRefDoa, repBundle.caseBundle.caseMsDoa, ...
    repBundle.caseBundle.caseStaticRefOnly, repBundle.caseBundle.caseStaticMs, ...
    repBundle.caseBundle.weightCase];
  caseTable = buildDoaDopplerSummaryTable(mainCase, truth, struct('mode', 'static'));
  ablationTruth = [ ...
    buildDoaDopplerCaseTruthFromScene(truth, sceneRefOnly), ...
    buildDoaDopplerCaseTruthFromScene(truth, sceneOtherOnly)];
  ablationTable = buildDoaDopplerCaseSummaryTable( ...
    [repBundle.caseBundle.caseStaticRefAbl, repBundle.caseBundle.caseStaticOtherOnly], ...
    ablationTruth);

  msSatDiagTable = localBuildStaticSatDiagTable( ...
    repBundle.caseBundle.caseStaticMs.estResult, truth, scene);

  currentBestWeightIdx = bestWeightIdx(repIdx);
  if isfinite(currentBestWeightIdx) && currentBestWeightIdx >= 1 && ...
      currentBestWeightIdx <= numel(weightSweepAlpha)
    bestWeightSatDiagTable = localBuildStaticSatDiagTable( ...
      repBundle.caseBundle.weightCase(currentBestWeightIdx).estResult, truth, scene);
  else
    bestWeightSatDiagTable = table();
  end

  repCaseSummary(iRep).repeatIdx = repIdx;
  repCaseSummary(iRep).taskSeed = repTable.taskSeed(iRep);
  repCaseSummary(iRep).label = repTable.label(iRep);
  repCaseSummary(iRep).caseTable = caseTable;
  repCaseSummary(iRep).ablationTable = ablationTable;
  repCaseSummary(iRep).weightPathDiagTable = localBuildWeightSolvePathDiagTable( ...
    repBundle.caseBundle.caseMsDoa, repBundle.caseBundle.caseStaticRefOnly, ...
    repBundle.caseBundle.weightCase, weightSweepAlpha, truth);
  repCaseSummary(iRep).msSatDiagTable = msSatDiagTable;
  repCaseSummary(iRep).bestWeightSatDiagTable = bestWeightSatDiagTable;
end

fprintf('\n========== Representative repeats ==========%s', newline);
disp(repTable);

for iRep = 1:repCount
  fprintf('\n========== Representative repeat: %s (repeat %d, taskSeed %d) ==========%s', ...
    repCaseSummary(iRep).label, repCaseSummary(iRep).repeatIdx, repCaseSummary(iRep).taskSeed, newline);
  disp(repCaseSummary(iRep).caseTable);

  fprintf('---------- Static ablation ----------%s', newline);
  disp(repCaseSummary(iRep).ablationTable);

  fprintf('---------- Weight-sweep solve-path diagnostics ----------%s', newline);
  disp(repCaseSummary(iRep).weightPathDiagTable);

  fprintf('---------- MS-SF-Static per-sat diagnostics ----------%s', newline);
  disp(repCaseSummary(iRep).msSatDiagTable);

  if ~isempty(repCaseSummary(iRep).bestWeightSatDiagTable) && ...
      height(repCaseSummary(iRep).bestWeightSatDiagTable) > 0
    fprintf('---------- Best-weight MS-SF-Static per-sat diagnostics ----------%s', newline);
    disp(repCaseSummary(iRep).bestWeightSatDiagTable);
  end
end

%% One-shot CRB anchor
crbOpt = struct();
crbOpt.doaType = 'latlon';
[crbSs, auxCrbSs] = crbPilotSfDoaDoppler( ...
  sceneRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);
[crbMs, auxCrbMs] = crbPilotSfDoaDoppler( ...
  scene, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);
crbSummary = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs, snrDb);

fprintf('\n========== Static CRB summary ==========%s', newline);
disp(crbSummary);

%% Plots
localPlotStaticRepeatScatter(repeatCompareTable);
localPlotWeightSweep(weightSweepAlpha, angleErrDeg(:, idxWeightStart:end), ...
  fdErrHz(:, idxWeightStart:end), isResolved(:, idxWeightStart:end));
localPlotRepresentativeGain(repeatCompareTable);

%% Optional snapshot save
% saveExpSnapshot("doaDopplerStatDualSatUraEci");

%% Local functions
function repeatOut = localRunStaticRepeatTask(taskSeed, pwrNoise, steeringInfo, pilotWave, ...
  linkParam, carrierFreq, sampleRate, pathGain, snapOpt, viewRefTemplate, ...
  viewOtherTemplate, viewMsTemplate, refSatIdxLocal, otherSatIdxLocal, wavelen, ...
  fdRange, optVerbose, truth, otherSatIdxGlobal, doaOnlyOpt, staticOptBase, ...
  weightSweepAlpha, staticMsHalfWidth, numCase, idxWeightStart)
%LOCALRUNSTATICREPEATTASK Run one fixed-SNR static repeat.

rng(taskSeed, 'twister');
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, sampleRate, ...
  pwrNoise, pathGain, snapOpt);

viewRef = viewRefTemplate;
viewRef.rxSigSf = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');

viewOther = viewOtherTemplate;
viewOther.rxSigSf = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');

viewMs = viewMsTemplate;
viewMs.rxSigSf = rxSig;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRef, viewOther, viewMs, wavelen, pilotWave, carrierFreq, ...
  sampleRate, fdRange, truth, otherSatIdxGlobal, optVerbose, ...
  doaOnlyOpt, staticOptBase, weightSweepAlpha, staticMsHalfWidth);
caseList = [caseBundle.caseRefDoa, caseBundle.caseMsDoa, ...
  caseBundle.caseStaticRefOnly, caseBundle.caseStaticMs, caseBundle.weightCase];
truthFdList = [NaN, NaN, truth.fdSatTrueHz(truth.refSatIdxLocal), truth.fdRefTrueHz, ...
  repmat(truth.fdRefTrueHz, 1, numel(weightSweepAlpha))];

angleVec = nan(1, numCase);
fdVec = nan(1, numCase);
resolvedVec = false(1, numCase);
for iCase = 1:numCase
  [angleVec(iCase), fdVec(iCase), resolvedVec(iCase)] = localExtractCaseMetric( ...
    caseList(iCase), truth.latlonTrueDeg, truthFdList(iCase));
end

weightAngle = angleVec(idxWeightStart:end);
weightFd = fdVec(idxWeightStart:end);
weightResolved = resolvedVec(idxWeightStart:end) & isfinite(weightAngle);
bestWeightIdx = NaN;
bestWeightAlpha = NaN;
bestWeightAngleErrDeg = NaN;
bestWeightFdErrHz = NaN;
if any(weightResolved)
  validIdx = find(weightResolved);
  [~, bestLocal] = min(weightAngle(weightResolved));
  bestWeightIdx = validIdx(bestLocal);
  bestWeightAlpha = weightSweepAlpha(bestWeightIdx);
  bestWeightAngleErrDeg = weightAngle(bestWeightIdx);
  bestWeightFdErrHz = weightFd(bestWeightIdx);
end

repeatOut = struct();
repeatOut.angleErrDeg = angleVec;
repeatOut.fdErrHz = fdVec;
repeatOut.isResolved = resolvedVec;
repeatOut.bestWeightIdx = bestWeightIdx;
repeatOut.bestWeightAlpha = bestWeightAlpha;
repeatOut.bestWeightAngleErrDeg = bestWeightAngleErrDeg;
repeatOut.bestWeightFdErrHz = bestWeightFdErrHz;
end


function repOut = localRunStaticRepresentativeBundle(taskSeed, pwrNoise, steeringInfo, ...
  pilotWave, linkParam, carrierFreq, sampleRate, pathGain, snapOpt, ...
  viewRefTemplate, viewOtherTemplate, viewMsTemplate, refSatIdxLocal, ...
  otherSatIdxLocal, wavelen, fdRange, optVerbose, truth, otherSatIdxGlobal, ...
  doaOnlyOpt, staticOptBase, weightSweepAlpha, staticMsHalfWidth)
%LOCALRUNSTATICREPRESENTATIVEBUNDLE Rerun one repeat and keep the full case bundle.

rng(taskSeed, 'twister');
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, sampleRate, ...
  pwrNoise, pathGain, snapOpt);

viewRef = viewRefTemplate;
viewRef.rxSigSf = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');

viewOther = viewOtherTemplate;
viewOther.rxSigSf = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');

viewMs = viewMsTemplate;
viewMs.rxSigSf = rxSig;

repOut = struct();
repOut.caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRef, viewOther, viewMs, wavelen, pilotWave, carrierFreq, ...
  sampleRate, fdRange, truth, otherSatIdxGlobal, optVerbose, ...
  doaOnlyOpt, staticOptBase, weightSweepAlpha, staticMsHalfWidth);
end


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


function summaryTable = localBuildCaseSummaryTable(caseName, snrDb, angleErrDeg, fdErrHz, isResolved)
%LOCALBUILDCASESUMMARYTABLE Build one compact single-SNR case summary.

numCase = numel(caseName);
angleRmseDeg = nan(numCase, 1);
angleP95Deg = nan(numCase, 1);
fdRmseHz = nan(numCase, 1);
fdP95Hz = nan(numCase, 1);
resolveRate = nan(numCase, 1);

for iCase = 1:numCase
  angleFail = ~isResolved(:, iCase) | ~isfinite(angleErrDeg(:, iCase));
  angleStat = summarizeMonteCarloStat(angleErrDeg(:, iCase), angleFail);
  angleRmseDeg(iCase) = angleStat.rmse;
  angleP95Deg(iCase) = angleStat.p95;
  resolveRate(iCase) = 1 - angleStat.failRate;

  if any(isfinite(fdErrHz(:, iCase)))
    fdFail = ~isResolved(:, iCase) | ~isfinite(fdErrHz(:, iCase));
    fdStat = summarizeMonteCarloStat(fdErrHz(:, iCase), fdFail);
    fdRmseHz(iCase) = fdStat.rmse;
    fdP95Hz(iCase) = fdStat.p95;
  end
end

summaryTable = table(caseName(:), repmat(snrDb, numCase, 1), angleRmseDeg, angleP95Deg, ...
  fdRmseHz, fdP95Hz, resolveRate, ...
  'VariableNames', {'displayName', 'snrDb', 'angleRmseDeg', 'angleP95Deg', ...
  'fdRmseHz', 'fdP95Hz', 'resolveRate'});
end


function gainTable = localBuildStaticGainSummaryTable(angleErrDeg, fdErrHz, isResolved, ...
  bestWeightAngleErrDeg, bestWeightFdErrHz, idxSsStat, idxMsStat)
%LOCALBUILDSTATICGAINSUMMARYTABLE Build one compact SS/MS static gain table.

ssAngle = angleErrDeg(:, idxSsStat);
msAngle = angleErrDeg(:, idxMsStat);
ssFd = fdErrHz(:, idxSsStat);
msFd = fdErrHz(:, idxMsStat);
maskStatic = isResolved(:, idxSsStat) & isResolved(:, idxMsStat) & ...
  isfinite(ssAngle) & isfinite(msAngle);
maskFd = maskStatic & isfinite(ssFd) & isfinite(msFd);
maskBestWeight = isfinite(bestWeightAngleErrDeg) & maskStatic;
maskBestWeightFd = isfinite(bestWeightFdErrHz) & maskFd;

angleDiff = msAngle - ssAngle;
fdDiff = msFd - ssFd;
bestWeightVsSs = bestWeightAngleErrDeg - ssAngle;
bestWeightVsMs = bestWeightAngleErrDeg - msAngle;

metricName = ["numRepeat"; "msBetterAngleRate"; "msBetterFdRate"; ...
  "medianMsMinusSsAngleDeg"; "medianMsMinusSsFdHz"; ...
  "bestWeightBetterThanSsRate"; "bestWeightBetterThanMsRate"; ...
  "medianBestWeightMinusSsAngleDeg"; "medianBestWeightMinusMsAngleDeg"];
metricValue = [ ...
  numel(ssAngle); ...
  localMeanLogical(angleDiff(maskStatic) < 0); ...
  localMeanLogical(fdDiff(maskFd) < 0); ...
  median(angleDiff(maskStatic), 'omitnan'); ...
  median(fdDiff(maskFd), 'omitnan'); ...
  localMeanLogical(bestWeightVsSs(maskBestWeight) < 0); ...
  localMeanLogical(bestWeightVsMs(maskBestWeight) < 0); ...
  median(bestWeightVsSs(maskBestWeight), 'omitnan'); ...
  median(bestWeightVsMs(maskBestWeight), 'omitnan')];

gainTable = table(metricName, metricValue, ...
  'VariableNames', {'metricName', 'metricValue'});
end


function weightPrefTable = localBuildWeightPreferenceTable(weightAlpha, bestWeightIdx, ...
  weightAngleErrDeg, weightFdErrHz, weightResolved)
%LOCALBUILDWEIGHTPREFERENCETABLE Build one compact alpha-preference summary.

numWeight = numel(weightAlpha);
numBest = zeros(numWeight, 1);
angleRmseDeg = nan(numWeight, 1);
fdRmseHz = nan(numWeight, 1);
resolveRate = nan(numWeight, 1);

bestWeightIdx = reshape(bestWeightIdx, [], 1);
for iWeight = 1:numWeight
  numBest(iWeight) = sum(bestWeightIdx == iWeight);

  failAngle = ~weightResolved(:, iWeight) | ~isfinite(weightAngleErrDeg(:, iWeight));
  angleStat = summarizeMonteCarloStat(weightAngleErrDeg(:, iWeight), failAngle);
  angleRmseDeg(iWeight) = angleStat.rmse;
  resolveRate(iWeight) = 1 - angleStat.failRate;

  failFd = ~weightResolved(:, iWeight) | ~isfinite(weightFdErrHz(:, iWeight));
  fdStat = summarizeMonteCarloStat(weightFdErrHz(:, iWeight), failFd);
  fdRmseHz(iWeight) = fdStat.rmse;
end

weightPrefTable = table(weightAlpha(:), numBest, angleRmseDeg, fdRmseHz, resolveRate, ...
  'VariableNames', {'alphaSat2', 'numBestRepeat', 'angleRmseDeg', 'fdRmseHz', 'resolveRate'});
end


function repeatTable = localBuildRepeatCompareTable(repeatIdx, taskSeed, snrDb, angleErrDeg, fdErrHz, isResolved, ...
  bestWeightAlpha, bestWeightAngleErrDeg, bestWeightFdErrHz, idxSsStat, idxMsStat)
%LOCALBUILDREPEATCOMPARETABLE Build one compact repeat-level compare table.

repeatTable = table();
repeatTable.repeatIdx = repeatIdx;
repeatTable.taskSeed = taskSeed;
repeatTable.snrDb = repmat(snrDb, numel(repeatIdx), 1);
repeatTable.ssAngleErrDeg = angleErrDeg(:, idxSsStat);
repeatTable.msAngleErrDeg = angleErrDeg(:, idxMsStat);
repeatTable.msMinusSsAngleDeg = angleErrDeg(:, idxMsStat) - angleErrDeg(:, idxSsStat);
repeatTable.ssFdErrHz = fdErrHz(:, idxSsStat);
repeatTable.msFdErrHz = fdErrHz(:, idxMsStat);
repeatTable.msMinusSsFdHz = fdErrHz(:, idxMsStat) - fdErrHz(:, idxSsStat);
repeatTable.ssResolved = isResolved(:, idxSsStat);
repeatTable.msResolved = isResolved(:, idxMsStat);
repeatTable.bestAlphaSat2 = bestWeightAlpha;
repeatTable.bestWeightAngleErrDeg = bestWeightAngleErrDeg;
repeatTable.bestWeightFdErrHz = bestWeightFdErrHz;
repeatTable.bestWeightMinusSsAngleDeg = bestWeightAngleErrDeg - angleErrDeg(:, idxSsStat);
repeatTable.bestWeightMinusMsAngleDeg = bestWeightAngleErrDeg - angleErrDeg(:, idxMsStat);
end


function repTable = localSelectRepresentativeRepeats(repeatTable)
%LOCALSELECTREPRESENTATIVEREPEATS Select a few representative repeats.

labelList = strings(0, 1);
idxList = zeros(0, 1);

[bestGainIdx, hasBestGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'min');
if hasBestGain
  labelList(end+1, 1) = "bestMsAngleGain";
  idxList(end+1, 1) = repeatTable.repeatIdx(bestGainIdx);
end

[worstGainIdx, hasWorstGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'max');
if hasWorstGain
  labelList(end+1, 1) = "worstMsAngleGain";
  idxList(end+1, 1) = repeatTable.repeatIdx(worstGainIdx);
end

[medianMsIdx, hasMedianMs] = localSelectMedianResolved(repeatTable.msAngleErrDeg, repeatTable.msResolved);
if hasMedianMs
  labelList(end+1, 1) = "medianMsStatic";
  idxList(end+1, 1) = repeatTable.repeatIdx(medianMsIdx);
end

[bestWeightLiftIdx, hasWeightLift] = localSelectByMetric(repeatTable.bestWeightMinusMsAngleDeg, 'min');
if hasWeightLift
  labelList(end+1, 1) = "bestWeightLift";
  idxList(end+1, 1) = repeatTable.repeatIdx(bestWeightLiftIdx);
end

[idxUnique, keepIdx] = unique(idxList, 'stable');
taskSeedList = nan(numel(idxUnique), 1);
for iRow = 1:numel(idxUnique)
  matchIdx = find(repeatTable.repeatIdx == idxUnique(iRow), 1, 'first');
  if ~isempty(matchIdx)
    taskSeedList(iRow) = repeatTable.taskSeed(matchIdx);
  end
end
repTable = table(labelList(keepIdx), idxUnique, taskSeedList, ...
  'VariableNames', {'label', 'repeatIdx', 'taskSeed'});
end


function [selIdx, hasValue] = localSelectByMetric(metricVec, modeName)
%LOCALSELECTBYMETRIC Select one index by min or max finite metric.

metricVec = reshape(metricVec, [], 1);
finiteMask = isfinite(metricVec);
selIdx = NaN;
hasValue = any(finiteMask);
if ~hasValue
  return;
end
switch modeName
  case 'min'
    [~, idxLocal] = min(metricVec(finiteMask));
  case 'max'
    [~, idxLocal] = max(metricVec(finiteMask));
  otherwise
    error('doaDopplerStatDualSatUraEci:InvalidSelectMode', ...
      'Unsupported select mode: %s.', modeName);
end
finiteIdx = find(finiteMask);
selIdx = finiteIdx(idxLocal);
end


function [selIdx, hasValue] = localSelectMedianResolved(metricVec, resolvedMask)
%LOCALSELECTMEDIANRESOLVED Select one repeat closest to the resolved median.

metricVec = reshape(metricVec, [], 1);
resolvedMask = reshape(resolvedMask, [], 1) & isfinite(metricVec);
selIdx = NaN;
hasValue = any(resolvedMask);
if ~hasValue
  return;
end
validVal = metricVec(resolvedMask);
medianVal = median(validVal, 'omitnan');
validIdx = find(resolvedMask);
[~, idxLocal] = min(abs(validVal - medianVal));
selIdx = validIdx(idxLocal);
end


function pathDiagTable = localBuildWeightSolvePathDiagTable(caseMsDoa, caseStaticRefOnly, weightCase, weightAlpha, truth)
%LOCALBUILDWEIGHTSOLVEPATHDIAGTABLE Build one weight-sweep solve-path table.

numWeight = numel(weightAlpha);
msDoaAnchor = caseMsDoa.estResult.doaParamEst(:);
w0Doa = caseStaticRefOnly.estResult.doaParamEst(:);
w0FdRef = caseStaticRefOnly.estResult.fdRefEst;
truthLatlon = truth.latlonTrueDeg(:);

alphaSat2 = reshape(weightAlpha, [], 1);
angleErrDeg = nan(numWeight, 1);
fdRefEstHz = nan(numWeight, 1);
fdRefShiftVsW0Hz = nan(numWeight, 1);
angleToMsDoaAnchorDeg = nan(numWeight, 1);
angleToW0Deg = nan(numWeight, 1);
funcCount = nan(numWeight, 1);
iterations = nan(numWeight, 1);

for iWeight = 1:numWeight
  caseInfo = weightCase(iWeight);
  if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
    continue;
  end

  estResult = caseInfo.estResult;
  doaNow = getDoaDopplerLatlonEst(estResult);
  if all(isfinite(doaNow))
    angleErrDeg(iWeight) = calcLatlonAngleError(doaNow(:), truthLatlon);
    angleToMsDoaAnchorDeg(iWeight) = calcLatlonAngleError(doaNow(:), msDoaAnchor);
    angleToW0Deg(iWeight) = calcLatlonAngleError(doaNow(:), w0Doa);
  end

  fdRefNow = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
  fdRefEstHz(iWeight) = fdRefNow;
  if isfinite(fdRefNow)
    fdRefShiftVsW0Hz(iWeight) = fdRefNow - w0FdRef;
  end

  optimInfo = getDoaDopplerFieldOrDefault(estResult, 'optimInfo', struct());
  funcCount(iWeight) = getDoaDopplerFieldOrDefault(optimInfo, 'funcCount', NaN);
  iterations(iWeight) = getDoaDopplerFieldOrDefault(optimInfo, 'iterations', NaN);
end

pathDiagTable = table(alphaSat2, angleErrDeg, fdRefEstHz, fdRefShiftVsW0Hz, ...
  angleToMsDoaAnchorDeg, angleToW0Deg, funcCount, iterations, ...
  'VariableNames', {'alphaSat2', 'angleErrDeg', 'fdRefEstHz', ...
  'fdRefShiftVsW0Hz', 'angleToMsDoaAnchorDeg', 'angleToW0Deg', ...
  'funcCount', 'iterations'});
end


function satDiagTable = localBuildStaticSatDiagTable(estResult, truth, scene)
%LOCALBUILDSTATICSATDIAGTABLE Build one per-satellite static diagnostic table.

if isempty(estResult) || ~isstruct(estResult)
  satDiagTable = table();
  return;
end

aux = getDoaDopplerFieldOrDefault(estResult, 'aux', struct());
numSat = scene.numSat;
truthLocalDoa = reshape(scene.localDoa, 2, []);
estLocalDoa = localAlignLocalDoa(aux, numSat);
fdSatEst = reshape(getDoaDopplerFieldOrDefault(aux, 'fdSatEst', nan(numSat, 1)), [], 1);
deltaFdEst = reshape(getDoaDopplerFieldOrDefault(aux, 'deltaFdRefEst', nan(numSat, 1)), [], 1);
objectiveSat = reshape(getDoaDopplerFieldOrDefault(aux, 'objectiveSat', nan(numSat, 1)), [], 1);
residualNormSat = reshape(getDoaDopplerFieldOrDefault(aux, 'residualNormSat', nan(numSat, 1)), [], 1);
noiseVarSat = reshape(getDoaDopplerFieldOrDefault(aux, 'noiseVarSat', nan(numSat, 1)), [], 1);
satWeight = reshape(getDoaDopplerFieldOrDefault(aux, 'satWeight', nan(numSat, 1)), [], 1);
fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);

localDirErrDeg = nan(numSat, 1);
azErrDeg = nan(numSat, 1);
elErrDeg = nan(numSat, 1);
for iSat = 1:numSat
  if all(isfinite(estLocalDoa(:, iSat))) && all(isfinite(truthLocalDoa(:, iSat)))
    localDirErrDeg(iSat) = calcDoaAngleError(estLocalDoa(:, iSat), truthLocalDoa(:, iSat), 'angle');
    azErrDeg(iSat) = rad2deg(estLocalDoa(1, iSat) - truthLocalDoa(1, iSat));
    elErrDeg(iSat) = rad2deg(estLocalDoa(2, iSat) - truthLocalDoa(2, iSat));
  end
end

fdConsistencyErrHz = fdSatEst - (fdRefEst + deltaFdEst);
refSatMask = false(numSat, 1);
refSatMask(truth.refSatIdxLocal) = true;

satDiagTable = table((1:numSat).', truth.selectedSatIdxGlobal(:), refSatMask, ...
  truth.usrElevationDeg(:), rad2deg(truthLocalDoa(1, :).'), rad2deg(truthLocalDoa(2, :).'), ...
  rad2deg(estLocalDoa(1, :).'), rad2deg(estLocalDoa(2, :).'), azErrDeg, elErrDeg, ...
  localDirErrDeg, truth.fdSatTrueHz(:), fdSatEst, fdSatEst - truth.fdSatTrueHz(:), ...
  truth.deltaFdTrueHz(:), deltaFdEst, deltaFdEst - truth.deltaFdTrueHz(:), ...
  fdConsistencyErrHz, objectiveSat, residualNormSat, noiseVarSat, satWeight, ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'isRefSat', 'usrElevationDeg', ...
  'truthAzDeg', 'truthElDeg', 'estAzDeg', 'estElDeg', 'azErrDeg', 'elErrDeg', ...
  'localDirErrDeg', 'fdSatTrueHz', 'fdSatEstHz', 'fdSatErrHz', 'deltaFdTrueHz', ...
  'deltaFdEstHz', 'deltaFdErrHz', 'fdConsistencyErrHz', 'objectiveSat', ...
  'residualNormSat', 'noiseVarSat', 'satWeight'});
end


function estLocalDoa = localAlignLocalDoa(aux, numSat)
%LOCALALIGNLOCALDOA Normalize static local-DoA storage to 2xNs.

estLocalDoa = nan(2, numSat);
localDoaEst = getDoaDopplerFieldOrDefault(aux, 'localDoaEst', []);
if isempty(localDoaEst)
  return;
end

if isequal(size(localDoaEst), [2, numSat])
  estLocalDoa = localDoaEst;
  return;
end
if ndims(localDoaEst) == 3 && size(localDoaEst, 1) == 2 && size(localDoaEst, 3) == numSat
  estLocalDoa = reshape(localDoaEst(:, 1, :), 2, numSat);
  return;
end
if ndims(localDoaEst) == 3 && size(localDoaEst, 1) == 2 && size(localDoaEst, 2) == numSat
  estLocalDoa = reshape(localDoaEst(:, :, 1), 2, numSat);
end
end


function crbTable = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs, snrDb)
%LOCALBUILDCRBSUMMARY Build one compact static CRB summary table.

caseName = ["SS-SF-Static"; "MS-SF-Static"];
satMode = ["single"; "multi"];
angleStdDeg = [ ...
  projectCrbToAngleMetric(crbSs(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMs(1:2, 1:2), truth.latlonTrueDeg, 'latlon')];
fdStdHz = [sqrt(max(real(crbSs(3, 3)), 0)); ...
  sqrt(max(real(crbMs(3, 3)), 0))];

crbTable = table(repmat(snrDb, 2, 1), caseName, satMode, angleStdDeg, fdStdHz, ...
  {auxCrbSs.fdRateMode; auxCrbMs.fdRateMode}, ...
  'VariableNames', {'snrDb', 'displayName', 'satMode', 'angleCrbStdDeg', ...
  'fdRefCrbStdHz', 'fdRateMode'});
end


function localPlotStaticRepeatScatter(repeatTable)
%LOCALPLOTSTATICREPEATSCATTER Plot one compact SS/MS repeat scatter view.

maskAngle = repeatTable.ssResolved & repeatTable.msResolved & ...
  isfinite(repeatTable.ssAngleErrDeg) & isfinite(repeatTable.msAngleErrDeg);
maskFd = maskAngle & isfinite(repeatTable.ssFdErrHz) & isfinite(repeatTable.msFdErrHz);

figure();
subplot(1, 2, 1);
scatter(repeatTable.ssAngleErrDeg(maskAngle), repeatTable.msAngleErrDeg(maskAngle), 28, 'filled');
hold on;
xyMax = max([repeatTable.ssAngleErrDeg(maskAngle); repeatTable.msAngleErrDeg(maskAngle)]);
if isempty(xyMax) || ~isfinite(xyMax)
  xyMax = 1;
end
plot([0, xyMax], [0, xyMax], '--k', 'LineWidth', 1.1);
grid on;
xlabel('SS-SF-Static angle error (deg)');
ylabel('MS-SF-Static angle error (deg)');
title('Repeat-level static angle compare');

subplot(1, 2, 2);
scatter(repeatTable.ssFdErrHz(maskFd), repeatTable.msFdErrHz(maskFd), 28, 'filled');
hold on;
xyMax = max([repeatTable.ssFdErrHz(maskFd); repeatTable.msFdErrHz(maskFd)]);
if isempty(xyMax) || ~isfinite(xyMax)
  xyMax = 1;
end
plot([0, xyMax], [0, xyMax], '--k', 'LineWidth', 1.1);
grid on;
xlabel('SS-SF-Static fd error (Hz)');
ylabel('MS-SF-Static fd error (Hz)');
title('Repeat-level static fd compare');
end


function localPlotWeightSweep(weightAlpha, weightAngleErrDeg, weightFdErrHz, weightResolved)
%LOCALPLOTWEIGHTSWEEP Plot one compact fixed-SNR weight sweep summary.

numWeight = numel(weightAlpha);
angleRmseDeg = nan(numWeight, 1);
fdRmseHz = nan(numWeight, 1);
for iWeight = 1:numWeight
  angleFail = ~weightResolved(:, iWeight) | ~isfinite(weightAngleErrDeg(:, iWeight));
  angleStat = summarizeMonteCarloStat(weightAngleErrDeg(:, iWeight), angleFail);
  angleRmseDeg(iWeight) = angleStat.rmse;

  fdFail = ~weightResolved(:, iWeight) | ~isfinite(weightFdErrHz(:, iWeight));
  fdStat = summarizeMonteCarloStat(weightFdErrHz(:, iWeight), fdFail);
  fdRmseHz(iWeight) = fdStat.rmse;
end

figure();
subplot(1, 2, 1);
plot(weightAlpha, angleRmseDeg, 'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Angle RMSE (deg)');
title('Fixed-SNR static angle vs sat2 weight');

subplot(1, 2, 2);
plot(weightAlpha, fdRmseHz, 'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Reference Doppler RMSE (Hz)');
title('Fixed-SNR static fd vs sat2 weight');
end


function localPlotRepresentativeGain(repeatTable)
%LOCALPLOTREPRESENTATIVEGAIN Plot repeat-wise MS-SS gain traces.

figure();
subplot(2, 1, 1);
plot(repeatTable.repeatIdx, repeatTable.msMinusSsAngleDeg, '-o', 'LineWidth', 1.1, 'MarkerSize', 4);
yline(0, '--k', 'LineWidth', 1.0);
grid on;
xlabel('Repeat index');
ylabel('MS - SS angle error (deg)');
title('Repeat-wise static angle gain/loss');

subplot(2, 1, 2);
plot(repeatTable.repeatIdx, repeatTable.msMinusSsFdHz, '-o', 'LineWidth', 1.1, 'MarkerSize', 4);
yline(0, '--k', 'LineWidth', 1.0);
grid on;
xlabel('Repeat index');
ylabel('MS - SS fd error (Hz)');
title('Repeat-wise static fd gain/loss');
end


function value = localMeanLogical(mask)
%LOCALMEANLOGICAL Return mean(mask) with empty handling.

if isempty(mask)
  value = NaN;
  return;
end
value = mean(double(mask));
end
