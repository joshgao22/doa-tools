clear(); close all; clc;

%% Parameters
numUsr = 1;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 11.7e9;
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

utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
% tle = tleread(resolveTestDataPath("starlink_pair_4154_742_20260318_170800.tle", "tle"));
arrUpa = createUpa(numElem, elemSpace);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

fdRange = [-2e5, 2e5];
optVerbose = false;
weightSweepAlpha = [0, 0.25, 0.5, 1];

if numUsr ~= 1
  error('doaDopplerStatDualSatUraEci:OnlySingleUserSupported', ...
    'This script currently supports one user only.');
end

%% Select Two Visible Satellites And The Reference Satellite
[~, satAccess] = findVisibleSatFromTle(utc, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccess, 2, 1, "available");
refSatIdxGlobal = satIdx(1);

%% Multi-Satellite Scene With Satellite Reference
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
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.refSatIdxGlobal = refSatIdxGlobal;
truth.refSatIdxLocal = refSatIdxLocal;
truth.selectedSatIdxGlobal = satIdx(:).';
truth.usrElevationDeg = reshape(scene.accessInfo.usrElevationDeg(:, 1), 1, []);
truth.fdRefTrueHz = truthDopplerState.fdRefRefFrame;
truth.fdSatTrueHz = reshape(truthDopplerState.fdSatRefFrame, [], 1);
truth.deltaFdTrueHz = reshape(truthDopplerState.deltaFdRefFrame, [], 1);
truth.localDoaRef = reshape(scene.localDoa(:, refSatIdxLocal), 2, 1);
fdRange = expandRangeToTruth(fdRange, [truth.fdRefTrueHz; truth.fdSatTrueHz(:)], 0.1, 2e4);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;
truth.fdRateTrueHzPerSec = NaN;

%% Pilot Waveform
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

%% Single-Frame Snapshots
snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';

pathGain = ones(scene.numSat, scene.numUser);
[rxSig, ~, ~, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGain, snapOpt);

%% Single-Satellite And Multi-Satellite Views
sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
rxSigRefOnly = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigRefOnly, gridSize, searchRange, E);

sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, gridSize, searchRange, E);

viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);

%% Common estimator options
doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticOptBase = struct();
staticOptBase.useLogObjective = true;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticOptBase, ...
  weightSweepAlpha, [0.01; 0.01]);

caseRefDoa = caseBundle.caseRefDoa;
caseMsDoa = caseBundle.caseMsDoa;
caseOtherDoa = caseBundle.caseOtherDoa;
staticRefOpt = caseBundle.staticRefOpt;
staticOtherOpt = caseBundle.staticOtherOpt;
staticMsOpt = caseBundle.staticMsOpt;
caseStaticRefOnly = caseBundle.caseStaticRefOnly;
caseStaticRefAbl = caseBundle.caseStaticRefAbl;
caseStaticOtherOnly = caseBundle.caseStaticOtherOnly;
caseStaticMs = caseBundle.caseStaticMs;
weightCase = caseBundle.weightCase;

mainCase = [caseRefDoa, caseMsDoa, caseStaticRefOnly, caseStaticMs, weightCase];

%% ## CRB Calculation
crbOpt = struct();
crbOpt.doaType = 'latlon';
[crbSs, auxCrbSs] = crbPilotSfDoaDoppler( ...
  sceneRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);

[crbMs, auxCrbMs] = crbPilotSfDoaDoppler( ...
  scene, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);

crbSummary = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs);

%% Summaries
estTable = buildDoaDopplerSummaryTable(mainCase, truth, struct('mode', 'static'));

ablationTruth = [ ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneRefOnly), ...
  buildDoaDopplerCaseTruthFromScene(truth, sceneOtherOnly)];
ablationTable = buildDoaDopplerCaseSummaryTable([caseStaticRefAbl, caseStaticOtherOnly], ablationTruth);

weightTable = buildDoaDopplerWeightSweepTable(weightSweepAlpha, weightCase, truth);

fprintf('\n========== Selected satellites ==========%s', newline);
disp(table((1:scene.numSat).', truth.selectedSatIdxGlobal(:), truth.usrElevationDeg(:), ...
  truth.fdSatTrueHz(:), truth.deltaFdTrueHz(:), ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg', 'fdSatTrueHz', 'deltaFdTrueHz'}));

fprintf('\n========== Truth ==========%s', newline);
disp(table(truth.latlonTrueDeg(1), truth.latlonTrueDeg(2), truth.fdRefTrueHz, ...
  truth.refSatIdxLocal, truth.refSatIdxGlobal, otherSatIdxLocal, otherSatIdxGlobal, ...
  'VariableNames', {'latTrueDeg', 'lonTrueDeg', 'fdRefTrueHz', ...
  'refSatIdxLocal', 'refSatIdxGlobal', 'otherSatIdxLocal', 'otherSatIdxGlobal'}));

fprintf('\n========== Estimator summary ==========%s', newline);
disp(estTable);

fprintf('\n========== Static ablation summary ==========%s', newline);
disp(ablationTable);

fprintf('\n========== Static weight-sweep summary ==========%s', newline);
disp(weightTable);

fprintf('\n========== Static CRB summary ==========%s', newline);
disp(crbSummary);

badMask = ~estTable.isResolved;
if any(badMask)
  warning('doaDopplerStatDualSatUraEci:EstimatorUnresolved', ...
    'Some estimators are not resolved:\n%s', evalc('disp(estTable(badMask, :))'));
end

%% Plots
plotDoaDopplerGeometryComparison(truth.latlonTrueDeg, mainCase);
localPlotFdComparison(truth.fdRefTrueHz, mainCase, crbSummary);
localPlotWeightSweep(weightTable);

%% Local functions
function crbTable = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs)
%LOCALBUILDCRBSUMMARY Build a compact static CRB summary table.

caseName = ["SS-SF-Static"; "MS-SF-Static"];
satMode = ["single"; "multi"];
angleStdDeg = [ ...
  projectCrbToAngleMetric(crbSs(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMs(1:2, 1:2), truth.latlonTrueDeg, 'latlon')];
fdStdHz = [sqrt(max(real(crbSs(3, 3)), 0)); ...
  sqrt(max(real(crbMs(3, 3)), 0))];

crbTable = table(caseName, satMode, angleStdDeg, fdStdHz, ...
  {auxCrbSs.fdRateMode; auxCrbMs.fdRateMode}, ...
  'VariableNames', {'displayName', 'satMode', 'angleCrbStdDeg', ...
  'fdRefCrbStdHz', 'fdRateMode'});
end

function localPlotFdComparison(fdRefTrueHz, caseResult, crbSummary)
%LOCALPLOTFDCOMPARISON Plot one-shot reference-Doppler estimates.

caseName = strings(1, numel(caseResult));
fdEst = nan(1, numel(caseResult));
for iCase = 1:numel(caseResult)
  caseName(iCase) = caseResult(iCase).displayName;
  fdEst(iCase) = getDoaDopplerFieldOrDefault(caseResult(iCase).estResult, 'fdRefEst', NaN);
end

xPos = 1:numel(caseName);
figure();
plot(xPos, fdEst, 'o', 'LineWidth', 1.3, 'MarkerSize', 7); hold on;
yline(fdRefTrueHz, '--k', 'LineWidth', 1.2);

grid on;
xticks(xPos);
xticklabels(cellstr(caseName));
xtickangle(20);
xlabel('Case');
ylabel('Reference Doppler (Hz)');
title('Reference-Doppler estimates');

staticMask = ismember(caseName, crbSummary.displayName.');
if any(staticMask)
  staticIdx = find(staticMask);
  staticName = caseName(staticMask);
  staticFdEst = fdEst(staticMask);
  fdStdHz = zeros(size(staticFdEst));
  for iCase = 1:numel(staticName)
    idxCrb = find(crbSummary.displayName == staticName(iCase), 1, 'first');
    fdStdHz(iCase) = crbSummary.fdRefCrbStdHz(idxCrb);
  end
  errorbar(staticIdx, staticFdEst, fdStdHz, 'LineStyle', 'none', ...
    'LineWidth', 1.1, 'CapSize', 10);
  legend({'estimate', 'truth', 'static CRB std'}, 'Location', 'best');
else
  legend({'estimate', 'truth'}, 'Location', 'best');
end
end

function localPlotWeightSweep(weightTable)
%LOCALPLOTWEIGHTSWEEP Plot one compact sat2-weight sweep summary.

figure();
plot(weightTable.alphaSat2, weightTable.angleErrDeg, 'o-', 'LineWidth', 1.3, ...
  'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Angle error (deg)');
title('MS-SF-Static angle error under sat2 weighting');
end

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end

