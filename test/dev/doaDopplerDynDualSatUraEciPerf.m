% DOADOPPLERDYNDUALSATURAECIPERF Performance comparison under dynamic geometry.
% Keep the Monte Carlo outputs local, but reuse the shared dev/common case
% builders so the performance script follows the same static->dynamic path
% as the current dual-satellite debugging scripts.
clear(); close all;

%% Parameters
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
pwrSource = 1;
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);

utc0 = datetime([2026, 03, 12, 08, 53, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_20260312.tle", "tle"));
arrUpa = createUpa(numElem, elemSpace);

numFrame = 10;
frameIntvlSec = 1 / 750;
refFrameIdx = round((numFrame + 1) / 2);
utcVec = utc0 + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

optVerbose = false;

snrDb = -20:2:10;
numParam = numel(snrDb);
numRepeat = 1000;

%% Robust evaluation options
robustOpt = struct();
robustOpt.numSigma = 4;
robustOpt.corrMode = 'winsor';
robustOpt.plotDiag = true;
robustOpt.reportStat = true;

%% Multi-frame scene
[~, satAccessRef] = findVisibleSatFromTle(utc0, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
refFrameIdx = sceneSeq.refFrameIdx;
sceneRef = sceneSeq.sceneCell{refFrameIdx};
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = satIdx(otherSatIdxLocal);

if sceneSeq.numSat ~= 2
  error('doaDopplerDynDualSatUraEciPerf:InvalidNumSat', ...
    'This script expects exactly two satellites.');
end

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end

truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.pickAux = satPickAux;
latlonTrue = usrLla(1:2, 1);
llaTrue = usrLla(:, 1);
fdRefTrueDyn = truth.fdRefFit;
fdRateTrueDyn = truth.fdRateFit;
fdRefTrueStat = truth.fdRefSeries(refFrameIdx);
fdRange = expandRangeToTruth([0 1e3], [truth.fdRefFit; reshape(truth.fdSatSeries(:, refFrameIdx), [], 1)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(getDoaDopplerFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth([-500 50], fdRateTruthCand, 0.1, 5e2);

%% Pilot waveform
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

%% Snapshot and estimator options
snapOpt = struct();
snapOpt.delayMode = 'phaseonly';
snapOpt.carrierPhaseMode = 'none';

simOpt = struct();
simOpt.steeringMode = 'framewise';
simOpt.accessMode = 'framewise';
simOpt.linkParamCell = linkParamCell;
simOpt.snapshotModelOpt = snapOpt;

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

dynKnownOpt = struct();
dynKnownOpt.useLogObjective = true;
dynKnownOpt.initFdCount = 81;
dynKnownOpt.useAccessMask = false;
dynKnownOpt.phaseMode = 'continuous';
dynKnownOpt.steeringMode = 'framewise';
dynKnownOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
dynKnownOpt.enableFdAliasUnwrap = true;
dynKnownOpt.initDoaHalfWidth = [0.003; 0.003];

dynUnknownOpt = dynKnownOpt;
dynUnknownOpt.initDoaHalfWidth = [0.002; 0.002];

sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
sceneOtherOnly = selectSatScene(sceneRef, otherSatIdxLocal);

viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(sceneRef, [], gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

%% Outcome buffers
msePosDoaSingle = zeros(numParam, 1);
msePosDoaJoint = zeros(numParam, 1);
msePosDdStat = zeros(numParam, 1);
msePosDdDyn = zeros(numParam, 1);

mseFdDdStat = zeros(numParam, 1);
mseFdDdDyn = zeros(numParam, 1);

crbPosDdStat = zeros(numParam, 1);
crbPosDdDyn = zeros(numParam, 1);
crbFdDdStat = zeros(numParam, 1);
crbFdDdDyn = zeros(numParam, 1);

medPosDdDyn = nan(numParam, 1);
p95PosDdDyn = nan(numParam, 1);
trimRmsePosDdDyn = nan(numParam, 1);
winRmsePosDdDyn = nan(numParam, 1);
corrRmsePosDdDyn = nan(numParam, 1);
failCntPosDdDyn = zeros(numParam, 1);
outlierCntPosDdDyn = zeros(numParam, 1);
failRatePosDdDyn = zeros(numParam, 1);
outlierRatePosDdDyn = zeros(numParam, 1);

medFdDdDyn = nan(numParam, 1);
p95FdDdDyn = nan(numParam, 1);
trimRmseFdDdDyn = nan(numParam, 1);
winRmseFdDdDyn = nan(numParam, 1);
corrRmseFdDdDyn = nan(numParam, 1);
failCntFdDdDyn = zeros(numParam, 1);
outlierCntFdDdDyn = zeros(numParam, 1);
failRateFdDdDyn = zeros(numParam, 1);
outlierRateFdDdDyn = zeros(numParam, 1);

%% Monte Carlo loop
progressbar('reset', numParam);

for iSnr = 1:numParam
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));

  sePosDoaSingle = nan(numRepeat, 1);
  sePosDoaJoint = nan(numRepeat, 1);
  sePosDdStat = nan(numRepeat, 1);
  sePosDdDyn = nan(numRepeat, 1);

  seFdDdStat = nan(numRepeat, 1);
  seFdDdDyn = nan(numRepeat, 1);

  isResolvedDdDyn = false(numRepeat, 1);
  exitflagDdDyn = nan(numRepeat, 1);
  fvalDdDyn = nan(numRepeat, 1);
  latlonValidDdDyn = false(numRepeat, 1);
  fdValidDdDyn = false(numRepeat, 1);

  parfor iRepeat = 1:numRepeat
    pathGainCell = buildFramePathGain(sceneSeq.numFrame, sceneSeq.numSat, sceneSeq.numUser);

    [rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots( ...
      sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
      pwrNoise, pathGainCell, simOpt);

    rxSigRef = rxSigCell{refFrameIdx};
    rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame');
    rxSigOtherOnly = selectRxSigBySat(rxSigRef, otherSatIdxLocal, 'singleFrame');

    viewRef = localUpdateMultiFrameView(viewRefTemplate, rxSigMfRefOnly);
    viewOther = viewOtherTemplate;
    viewOther.rxSigSf = rxSigOtherOnly;
    viewMs = localUpdateMultiFrameView(viewMsTemplate, rxSigCell);

    caseBundle = buildDoaDopplerStaticTransitionBundle( ...
      viewRef, viewOther, viewMs, wavelen, pilotWave, carrierFreq, ...
      waveInfo.sampleRate, fdRange, optVerbose, truth, otherSatIdxGlobal, ...
      doaOnlyOpt, staticBaseOpt, [0; 0.25; 0.5; 1], [0.01; 0.01]);

    caseDoaSingle = caseBundle.caseRefDoa;
    caseDoaJoint = caseBundle.caseMsDoa;
    caseDdStat = caseBundle.caseStaticMs;
    bestStaticMsCase = caseBundle.bestStaticMsCase;

    initParamStaticMs = buildDynamicInitParamFromCase(bestStaticMsCase, true, truth.fdRateFit);
    dynKnownOptCur = dynKnownOpt;
    dynKnownOptCur.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
    caseDynKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
      viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
      fdRange, fdRateRange, optVerbose, dynKnownOptCur, true, [], initParamStaticMs);

    initParamUnknownCpK = buildDynamicInitParamFromCase(caseDynKnown, false, caseDynKnown.estResult.fdRateEst);
    initParamUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, truth.fdRateFit);
    initCandUnknown = buildUnknownInitCandidateSet( ...
      caseDynKnown, bestStaticMsCase, initParamUnknownCpK, initParamUnknownStatic, ...
      dynUnknownOpt.initDoaHalfWidth, caseBundle.staticMsOpt.initDoaHalfWidth);

    dynUnknownOptCur = dynUnknownOpt;
    dynUnknownOptCur.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
    estDdDynCase = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
      viewMs, truth, pilotWave, carrierFreq, waveInfo.sampleRate, ...
      fdRange, fdRateRange, optVerbose, dynUnknownOptCur, false, [], initCandUnknown);
    estDdDyn = estDdDynCase.estResult;

    latlonDoaSingle = getDoaDopplerLatlonEst(caseDoaSingle.estResult);
    if all(isfinite(latlonDoaSingle))
      sePosDoaSingle(iRepeat) = localCalcPosSe(latlonDoaSingle, llaTrue, E);
    end

    latlonDoaJoint = getDoaDopplerLatlonEst(caseDoaJoint.estResult);
    if all(isfinite(latlonDoaJoint))
      sePosDoaJoint(iRepeat) = localCalcPosSe(latlonDoaJoint, llaTrue, E);
    end

    latlonDdStat = getDoaDopplerLatlonEst(caseDdStat.estResult);
    if all(isfinite(latlonDdStat))
      sePosDdStat(iRepeat) = localCalcPosSe(latlonDdStat, llaTrue, E);
    end

    fdRefDdStat = getDoaDopplerFieldOrDefault(caseDdStat.estResult, 'fdRefEst', NaN);
    if isfinite(fdRefDdStat)
      seFdDdStat(iRepeat) = (fdRefDdStat - fdRefTrueStat)^2;
    end

    latlonDdDyn = getDoaDopplerLatlonEst(estDdDyn);
    fdRefDdDyn = getDoaDopplerFieldOrDefault(estDdDyn, 'fdRefEst', NaN);

    if all(isfinite(latlonDdDyn))
      sePosDdDyn(iRepeat) = localCalcPosSe(latlonDdDyn, llaTrue, E);
      latlonValidDdDyn(iRepeat) = true;
    end

    if isfinite(fdRefDdDyn)
      seFdDdDyn(iRepeat) = (fdRefDdDyn - fdRefTrueDyn)^2;
      fdValidDdDyn(iRepeat) = true;
    end

    isResolvedDdDyn(iRepeat) = logical(getDoaDopplerFieldOrDefault(estDdDyn, 'isResolved', false));
    exitflagDdDyn(iRepeat) = getDoaDopplerFieldOrDefault(estDdDyn, 'exitflag', NaN);
    fvalDdDyn(iRepeat) = getDoaDopplerFieldOrDefault(estDdDyn, 'fval', NaN);
  end

  msePosDoaSingle(iSnr) = localMeanFinite(sePosDoaSingle);
  msePosDoaJoint(iSnr) = localMeanFinite(sePosDoaJoint);
  msePosDdStat(iSnr) = localMeanFinite(sePosDdStat);
  msePosDdDyn(iSnr) = localMeanFinite(sePosDdDyn);

  mseFdDdStat(iSnr) = localMeanFinite(seFdDdStat);
  mseFdDdDyn(iSnr) = localMeanFinite(seFdDdDyn);

  failMaskDdDyn = ~isResolvedDdDyn | ...
                  ~isfinite(exitflagDdDyn) | ...
                  ~isfinite(fvalDdDyn) | ...
                  ~latlonValidDdDyn | ...
                  ~fdValidDdDyn | ...
                  ~isfinite(sePosDdDyn) | ...
                  ~isfinite(seFdDdDyn);

  posStatDdDyn = localSummarizeMcSe(sePosDdDyn, failMaskDdDyn, robustOpt.numSigma);
  fdStatDdDyn = localSummarizeMcSe(seFdDdDyn, failMaskDdDyn, robustOpt.numSigma);

  medPosDdDyn(iSnr) = posStatDdDyn.median;
  p95PosDdDyn(iSnr) = posStatDdDyn.p95;
  trimRmsePosDdDyn(iSnr) = posStatDdDyn.trimRmse;
  winRmsePosDdDyn(iSnr) = posStatDdDyn.winRmse;
  corrRmsePosDdDyn(iSnr) = localPickCorrRmse(posStatDdDyn, robustOpt.corrMode);
  failCntPosDdDyn(iSnr) = posStatDdDyn.numFail;
  outlierCntPosDdDyn(iSnr) = posStatDdDyn.numOutlier;
  failRatePosDdDyn(iSnr) = posStatDdDyn.numFail / numRepeat;
  outlierRatePosDdDyn(iSnr) = posStatDdDyn.numOutlier / numRepeat;

  medFdDdDyn(iSnr) = fdStatDdDyn.median;
  p95FdDdDyn(iSnr) = fdStatDdDyn.p95;
  trimRmseFdDdDyn(iSnr) = fdStatDdDyn.trimRmse;
  winRmseFdDdDyn(iSnr) = fdStatDdDyn.winRmse;
  corrRmseFdDdDyn(iSnr) = localPickCorrRmse(fdStatDdDyn, robustOpt.corrMode);
  failCntFdDdDyn(iSnr) = fdStatDdDyn.numFail;
  outlierCntFdDdDyn(iSnr) = fdStatDdDyn.numOutlier;
  failRateFdDdDyn(iSnr) = fdStatDdDyn.numFail / numRepeat;
  outlierRateFdDdDyn(iSnr) = fdStatDdDyn.numOutlier / numRepeat;

  progressbar('advance');
end

%% CRB evaluation
for iSnr = 1:numParam
  pwrNoise = pwrSource / (10^(snrDb(iSnr) / 10));

  crbStatOpt = struct();
  crbStatOpt.doaType = 'latlon';

  crbStat = crbPilotSfDoaDoppler( ...
    sceneRef, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    latlonTrue, fdRefTrueStat, ones(sceneRef.numSat, 1), pwrNoise, crbStatOpt);

  crbDynOpt = struct();
  crbDynOpt.doaType = 'latlon';
  crbDynOpt.phaseMode = 'continuous';
  crbDynOpt.fdRateMode = 'unknown';
  crbDynOpt.steeringMode = dynUnknownOpt.steeringMode;
  crbDynOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;

  crbDyn = crbPilotMfDoaDoppler( ...
    sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
    latlonTrue, fdRefTrueDyn, fdRateTrueDyn, ...
    ones(sceneSeq.numSat, sceneSeq.numFrame), pwrNoise, crbDynOpt);

  crbPosDdStat(iSnr) = localCalcPosCrb(crbStat(1:2, 1:2), llaTrue, E);
  crbPosDdDyn(iSnr) = localCalcPosCrb(crbDyn(1:2, 1:2), llaTrue, E);

  crbFdDdStat(iSnr) = real(crbStat(3, 3));
  crbFdDdDyn(iSnr) = real(crbDyn(3, 3));
end

%% Plot performance curves
rmsePosDoaSingle = sqrt(msePosDoaSingle);
rmsePosDoaJoint = sqrt(msePosDoaJoint);
rmsePosDdStat = sqrt(msePosDdStat);
rmsePosDdDyn = sqrt(msePosDdDyn);

rmseFdDdStat = sqrt(mseFdDdStat);
rmseFdDdDyn = sqrt(mseFdDdDyn);

rmseCrbPosDdStat = sqrt(crbPosDdStat);
rmseCrbPosDdDyn = sqrt(crbPosDdDyn);
rmseCrbFdDdStat = sqrt(crbFdDdStat);
rmseCrbFdDdDyn = sqrt(crbFdDdDyn);

figure();
tiledlayout(1, 2);

nexttile();
semilogy(snrDb, rmsePosDoaSingle, '-o', ...
         snrDb, rmsePosDoaJoint, '-s', ...
         snrDb, rmsePosDdStat, '-d', ...
         snrDb, rmsePosDdDyn, '-^', ...
         snrDb, rmseCrbPosDdStat, '--d', ...
         snrDb, rmseCrbPosDdDyn, '--^');
grid on;
xlabel('SNR (dB)');
ylabel('Position RMSE (m)');
legend({'Single-sat DoA', 'Dual-sat DoA', ...
        'Dual-sat Static DoA-Doppler', 'Dual-sat Dynamic DoA-Doppler', ...
        'CRB Static DoA-Doppler', 'CRB Dynamic DoA-Doppler'}, ...
        'Location', 'best');
title(sprintf('Position RMSE (%d frames)', sceneSeq.numFrame));

nexttile();
semilogy(snrDb, rmseFdDdStat, '-d', ...
         snrDb, rmseFdDdDyn, '-^', ...
         snrDb, rmseCrbFdDdStat, '--d', ...
         snrDb, rmseCrbFdDdDyn, '--^');
grid on;
xlabel('SNR (dB)');
ylabel('fdRef RMSE (Hz)');
legend({'Static DoA-Doppler', 'Dynamic DoA-Doppler', ...
        'CRB Static DoA-Doppler', 'CRB Dynamic DoA-Doppler'}, ...
        'Location', 'best');
title(sprintf('Reference Doppler RMSE (%d frames)', sceneSeq.numFrame));

if robustOpt.plotDiag
  corrLabel = localCorrLabel(robustOpt.corrMode);

  figure();
  tiledlayout(1, 2);

  nexttile();
  semilogy(snrDb, rmsePosDdDyn, '-^', ...
           snrDb, corrRmsePosDdDyn, '-v', ...
           snrDb, medPosDdDyn, '-o', ...
           snrDb, p95PosDdDyn, '-s', ...
           snrDb, rmseCrbPosDdDyn, '--^');
  grid on;
  xlabel('SNR (dB)');
  ylabel('Position Error (m)');
  legend({'Dynamic raw RMSE', corrLabel, 'Dynamic median', ...
          'Dynamic p95', 'CRB Dynamic DoA-Doppler'}, ...
          'Location', 'best');
  title('Dynamic Position Robust Diagnostics');

  nexttile();
  semilogy(snrDb, rmseFdDdDyn, '-^', ...
           snrDb, corrRmseFdDdDyn, '-v', ...
           snrDb, medFdDdDyn, '-o', ...
           snrDb, p95FdDdDyn, '-s', ...
           snrDb, rmseCrbFdDdDyn, '--^');
  grid on;
  xlabel('SNR (dB)');
  ylabel('fdRef Error (Hz)');
  legend({'Dynamic raw RMSE', corrLabel, 'Dynamic median', ...
          'Dynamic p95', 'CRB Dynamic DoA-Doppler'}, ...
          'Location', 'best');
  title('Dynamic fdRef Robust Diagnostics');

  figure();
  tiledlayout(1, 2);

  nexttile();
  plot(snrDb, 100 * failRatePosDdDyn, '-o', ...
       snrDb, 100 * outlierRatePosDdDyn, '-s');
  grid on;
  xlabel('SNR (dB)');
  ylabel('Rate (%)');
  legend({'Fail rate', 'Outlier rate'}, 'Location', 'best');
  title('Dynamic Position Fail / Outlier Rate');

  nexttile();
  plot(snrDb, 100 * failRateFdDdDyn, '-o', ...
       snrDb, 100 * outlierRateFdDdDyn, '-s');
  grid on;
  xlabel('SNR (dB)');
  ylabel('Rate (%)');
  legend({'Fail rate', 'Outlier rate'}, 'Location', 'best');
  title('Dynamic fdRef Fail / Outlier Rate');
end

if robustOpt.reportStat
  localPrintRobustSummary(snrDb, rmsePosDdDyn, corrRmsePosDdDyn, ...
    medPosDdDyn, p95PosDdDyn, failCntPosDdDyn, outlierCntPosDdDyn, ...
    'Dynamic position');

  localPrintRobustSummary(snrDb, rmseFdDdDyn, corrRmseFdDdDyn, ...
    medFdDdDyn, p95FdDdDyn, failCntFdDdDyn, outlierCntFdDdDyn, ...
    'Dynamic fdRef');
end

%% Optional snapshot save
saveExpSnapshot("doaDopplerDynDualSatUraEciPerf");

%% Local functions
function view = localUpdateMultiFrameView(viewTemplate, rxSigMf)
%LOCALUPDATEMULTIFRAMEVIEW Refresh one multi-frame view with current data.

view = viewTemplate;
view.rxSigMf = rxSigMf;
view.rxSigSf = rxSigMf{view.sceneSeq.refFrameIdx};
end


function meanVal = localMeanFinite(valVec)
%LOCALMEANFINITE Mean over finite samples only.

validMask = isfinite(valVec);
if ~any(validMask)
  meanVal = NaN;
  return;
end
meanVal = mean(valVec(validMask));
end


function posSe = localCalcPosSe(latlonEst, llaTrue, spheroid)
%LOCALCALCPOSSE Position squared error in meters.

posEstEcef = localLla2Ecef([latlonEst(:); llaTrue(3)], spheroid);
posTrueEcef = localLla2Ecef(llaTrue, spheroid);
errVec = posEstEcef - posTrueEcef;
posSe = errVec.' * errVec;
end


function posCrb = localCalcPosCrb(crbLatLon, llaTrue, spheroid)
%LOCALCALCPOSCRB Project one lat-lon CRB to position variance.

lat0 = llaTrue(1);
lon0 = llaTrue(2);
h0 = llaTrue(3);
stepDeg = 1e-4;

ecefLatPlus = localLla2Ecef([lat0 + stepDeg; lon0; h0], spheroid);
ecefLatMinus = localLla2Ecef([lat0 - stepDeg; lon0; h0], spheroid);
ecefLonPlus = localLla2Ecef([lat0; lon0 + stepDeg; h0], spheroid);
ecefLonMinus = localLla2Ecef([lat0; lon0 - stepDeg; h0], spheroid);

jacPos = [ ...
  (ecefLatPlus - ecefLatMinus) / (2 * stepDeg), ...
  (ecefLonPlus - ecefLonMinus) / (2 * stepDeg)];

posCrbMat = jacPos * crbLatLon * jacPos.';
posCrb = max(real(trace(posCrbMat)), 0);
end


function posEcef = localLla2Ecef(llaDeg, spheroid)
%LOCALLLA2ECEF Convert one LLA vector in degrees to ECEF coordinates.

[x, y, z] = geodetic2ecef(spheroid, llaDeg(1), llaDeg(2), llaDeg(3));
posEcef = [x; y; z];
end


function stat = localSummarizeMcSe(seVec, failMask, numSigma)
%LOCALSUMMARIZEMCSE Build robust diagnostics on one squared-error sequence.

errAbs = sqrt(seVec(:));
validMask = isfinite(errAbs) & ~failMask(:);
validErr = errAbs(validMask);

stat = struct();
stat.rmse = NaN;
stat.median = NaN;
stat.p95 = NaN;
stat.trimRmse = NaN;
stat.winRmse = NaN;
stat.numFail = sum(failMask(:) | ~isfinite(errAbs));
stat.numOutlier = 0;
stat.validMask = validMask;
stat.outlierMask = false(size(errAbs));

if isempty(validErr)
  return;
end

stat.rmse = sqrt(mean(validErr .^ 2));
stat.median = median(validErr);
stat.p95 = prctile(validErr, 95);

medErr = stat.median;
madErr = median(abs(validErr - medErr));
robustScale = 1.4826 * madErr;
if ~isfinite(robustScale) || robustScale <= 0
  robustScale = std(validErr, 'omitnan');
end
if ~isfinite(robustScale) || robustScale <= 0
  stat.trimRmse = stat.rmse;
  stat.winRmse = stat.rmse;
  return;
end

thresh = medErr + numSigma * robustScale;
outlierMaskValid = validErr > thresh;
stat.numOutlier = sum(outlierMaskValid);
outlierMask = false(size(errAbs));
validIdx = find(validMask);
outlierMask(validIdx(outlierMaskValid)) = true;
stat.outlierMask = outlierMask;

trimErr = validErr(~outlierMaskValid);
if isempty(trimErr)
  stat.trimRmse = NaN;
else
  stat.trimRmse = sqrt(mean(trimErr .^ 2));
end

winErr = validErr;
winErr(outlierMaskValid) = thresh;
stat.winRmse = sqrt(mean(winErr .^ 2));
end


function corrRmse = localPickCorrRmse(stat, corrMode)
%LOCALPICKCORRRMSE Pick one corrected RMSE metric.

switch lower(corrMode)
  case 'trim'
    corrRmse = stat.trimRmse;
  case 'winsor'
    corrRmse = stat.winRmse;
  otherwise
    error('doaDopplerDynDualSatUraEciPerf:InvalidCorrMode', ...
      'Unknown corrMode: %s', corrMode);
end
end


function corrLabel = localCorrLabel(corrMode)
%LOCALCORRLABEL Build one legend label for the corrected RMSE.

switch lower(corrMode)
  case 'trim'
    corrLabel = 'Dynamic trimmed RMSE';
  case 'winsor'
    corrLabel = 'Dynamic winsor RMSE';
  otherwise
    corrLabel = 'Dynamic corrected RMSE';
end
end


function localPrintRobustSummary(snrDb, rawRmse, corrRmse, medErr, p95Err, ...
  failCnt, outlierCnt, metricName)
%LOCALPRINTROBUSTSUMMARY Print one compact robust summary table.

summaryTable = table(snrDb(:), rawRmse(:), corrRmse(:), medErr(:), p95Err(:), ...
  failCnt(:), outlierCnt(:), ...
  'VariableNames', {'snrDb', 'rawRmse', 'corrRmse', 'median', 'p95', ...
  'failCount', 'outlierCount'});
fprintf('\n========== %s robust summary ==========%s', metricName, newline);
disp(summaryTable);
end
