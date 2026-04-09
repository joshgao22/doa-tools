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
snrDb = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));
E = referenceEllipsoid('sphere');

usrLla = [[37.78, 36.59, 0]', [37.58, 37.51, 0]'];
usrLla = usrLla(:, 1:numUsr);

utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(localResolveTlePath("starlink_pair_4154_1165_20260318_170800.tle"));
% tle = tleread(localResolveTlePath("starlink_pair_4154_742_20260318_170800.tle"));
arrUpa = createUpa(numElem, elemSpace);

gridSize = [50 50];
searchMarginDeg = 5;
searchRange = [usrLla(1, 1) - searchMarginDeg, usrLla(1, 1) + searchMarginDeg; ...
               usrLla(2, 1) - searchMarginDeg, usrLla(2, 1) + searchMarginDeg];

fdRange = [-2e5, 2e5];
optVerbose = false;
weightSweepAlpha = [0, 0.25, 0.5, 1];
sliceHalfWidthDeg = [0.02; 0.02];
sliceGridSize = [21 21];

if numUsr ~= 1
  error('doaDopplerStatDualSatUraEci:OnlySingleUserSupported', ...
    'This script currently supports one user only.');
end

%% Select Two Visible Satellites And The Reference Satellite
[~, satAccess] = findVisibleSatFromTle(utc, tle, usrLla);
[satIdx, satPickAux] = pickVisibleSatByElevation(satAccess, 2, 1, "available");
refSatIdxGlobal = satIdx(1);
%[text] ## Multi-Satellite Scene With Satellite Reference
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
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;

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
%[text] ## Single-Satellite And Multi-Satellite Views
sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
rxSigRefOnly = selectRxSigBySat(rxSig, refSatIdxLocal, 'singleFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigRefOnly, gridSize, searchRange, E);

sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSig, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, gridSize, searchRange, E);

sceneDupRef = selectSatScene(scene, [refSatIdxLocal, refSatIdxLocal]);
rxSigDupRef = selectRxSigBySat(rxSig, [refSatIdxLocal, refSatIdxLocal], 'singleFrame');
viewDupRef = buildDoaDopplerEstView(sceneDupRef, rxSigDupRef, gridSize, searchRange, E);

viewMs = buildDoaDopplerEstView(scene, rxSig, gridSize, searchRange, E);

%% Common Estimator Options
doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticOptBase = struct();
staticOptBase.useLogObjective = true;
%[text] ## DoA Initializers
caseRefDoa = runDoaOnlyCase("SS-SF-DoA", "single", viewRefOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);
caseMsDoa = runDoaOnlyCase("MS-SF-DoA", "multi",  viewMs, wavelen, pilotWave, optVerbose, doaOnlyOpt);
caseOtherDoa = runDoaOnlyCase("SS-SF-DoA-sat2Only", "single", viewOtherOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);

staticRefOpt = staticOptBase;
staticRefOpt.initDoaParam = caseRefDoa.estResult.doaParamEst(:);

staticOtherOpt = staticOptBase;
staticOtherOpt.initDoaParam = caseOtherDoa.estResult.doaParamEst(:);

staticMsOpt = staticOptBase;
staticMsOpt.initDoaParam = caseMsDoa.estResult.doaParamEst(:);
staticMsOpt.initDoaHalfWidth = [0.01; 0.01];

staticDupOpt = staticMsOpt;

%% Static baseline and ablation cases
caseStaticRefOnly = runStaticDoaDopplerCase("SS-SF-Static", "single", ...
  viewRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticRefOpt);
caseStaticRefAbl = caseStaticRefOnly;
caseStaticRefAbl.displayName = "MS-SF-Static-sat1Only";
caseStaticRefAbl.satMode = "multi";
caseStaticOtherOnly = runStaticDoaDopplerCase("MS-SF-Static-sat2Only", "multi", ...
  viewOtherOnly, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticOtherOpt);
caseStaticMs = runStaticDoaDopplerCase("MS-SF-Static", "multi", ...
  viewMs, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticMsOpt);
[caseStaticDupRef, dupRefDiag] = localTryRunStaticDoaDopplerCase("MS-SF-Static-DupRef", "multi", ...
  viewDupRef, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, staticDupOpt);

weightCase = repmat(buildDoaDopplerCaseResult("", "multi", "single", ...
  "doa-doppler", "static", struct()), 1, numel(weightSweepAlpha));
for iCase = 1:numel(weightSweepAlpha)
  alpha = weightSweepAlpha(iCase);
  currentOpt = staticMsOpt;
  currentOpt.satWeight = [1; alpha];
  weightCase(iCase) = runStaticDoaDopplerCase( ...
    sprintf('MS-SF-Static-W%.2f', alpha), "multi", viewMs, ...
    pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, optVerbose, currentOpt);
end

mainCase = [caseRefDoa, caseMsDoa, caseStaticRefOnly, caseStaticMs, caseStaticDupRef, weightCase];

%% ## CRB Calculation
crbOpt = struct();
crbOpt.doaType = 'latlon';
[crbSs, auxCrbSs] = crbPilotSfDoaDoppler( ...
  sceneRefOnly, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);

[crbMs, auxCrbMs] = crbPilotSfDoaDoppler( ...
  scene, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);

[crbDup, auxCrbDup] = localTryBuildStaticCrb( ...
  sceneDupRef, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefTrueHz, 1, pwrNoise, crbOpt);

crbSummary = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs, crbDup, auxCrbDup);

%% Summaries
estTable = buildDoaDopplerSummaryTable(mainCase, truth, struct('mode', 'static'));

ablationTruth = [ ...
  localBuildCaseTruth(truth.latlonTrueDeg, truth.fdSatTrueHz(refSatIdxLocal)), ...
  localBuildCaseTruth(truth.latlonTrueDeg, truth.fdSatTrueHz(otherSatIdxLocal))];
ablationTable = localBuildCaseSummaryTable([caseStaticRefAbl, caseStaticOtherOnly], ablationTruth);

weightTable = localBuildWeightSweepTable(weightSweepAlpha, weightCase, truth);

dupCaseTruth = localBuildCaseTruth(truth.latlonTrueDeg, truth.fdRefTrueHz);
dupTable = localBuildCaseSummaryTable(caseStaticDupRef, dupCaseTruth);
if isfield(dupRefDiag, 'errorMessage') && strlength(dupRefDiag.errorMessage) > 0
  dupTable.errorMessage = repmat(dupRefDiag.errorMessage, height(dupTable), 1);
end

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

fprintf('\n========== Duplicate-reference summary ==========%s', newline);
disp(dupTable);

fprintf('\n========== Static CRB summary ==========%s', newline);
disp(crbSummary);

badMask = ~estTable.isResolved;
if any(badMask)
  warning('doaDopplerStatDualSatUraEci:EstimatorUnresolved', ...
    'Some estimators are not resolved:\n%s', evalc('disp(estTable(badMask, :))'));
end
%[text] ## Objective slices
sliceDiag = cell(1, 4);
sliceDiag{1} = localBuildStaticObjectiveSlice("Ref-only", viewRefOnly, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth.fdSatTrueHz(refSatIdxLocal), ...
  truth.latlonTrueDeg, sliceHalfWidthDeg, sliceGridSize, staticRefOpt);
sliceDiag{2} = localBuildStaticObjectiveSlice("Other-only", viewOtherOnly, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth.fdSatTrueHz(otherSatIdxLocal), ...
  truth.latlonTrueDeg, sliceHalfWidthDeg, sliceGridSize, staticOtherOpt);
sliceDiag{3} = localBuildStaticObjectiveSlice("Multi-sat", viewMs, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth.fdRefTrueHz, ...
  truth.latlonTrueDeg, sliceHalfWidthDeg, sliceGridSize, staticMsOpt);
sliceDiag{4} = localTryBuildStaticObjectiveSlice("Dup-ref", viewDupRef, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth.fdRefTrueHz, ...
  truth.latlonTrueDeg, sliceHalfWidthDeg, sliceGridSize, staticDupOpt, dupRefDiag);
%[text] ## Plots
plotDoaDopplerGeometryComparison(truth.latlonTrueDeg, mainCase);
localPlotFdComparison(truth.fdRefTrueHz, mainCase, crbSummary);
localPlotWeightSweep(weightTable);
localPlotObjectiveSlice(sliceDiag, truth.latlonTrueDeg);

%% Local functions
function [caseInfo, diagInfo] = localTryRunStaticDoaDopplerCase(displayName, satMode, view, ...
  pilotWave, carrierFreq, sampleRate, fdRange, verbose, modelOpt)
%LOCALTRYRUNSTATICDOADOPPLERCASE Run one static case with graceful fallback.

diagInfo = struct();
diagInfo.errorMessage = "";

try
  caseInfo = runStaticDoaDopplerCase(displayName, satMode, view, ...
    pilotWave, carrierFreq, sampleRate, fdRange, verbose, modelOpt);
catch ME
  diagInfo.errorMessage = string(ME.message);
  estResult = struct();
  estResult.doaParamEst = [NaN; NaN];
  estResult.fdRefEst = NaN;
  estResult.fdRateEst = NaN;
  estResult.exitflag = NaN;
  estResult.isResolved = false;
  estResult.optimInfo = struct('runTimeSec', NaN, 'funcCount', NaN, 'iterations', NaN);
  estResult.aux = struct('errorMessage', diagInfo.errorMessage);
  caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
    "doa-doppler", "static", estResult);
end
end


function [crb, aux] = localTryBuildStaticCrb(scene, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRefHz, numSource, noiseVar, crbOpt)
%LOCALTRYBUILDSTATICCRB Build one static CRB with graceful fallback.

try
  [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, numSource, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = "static";
  aux.errorMessage = string(ME.message);
end
end


function sliceDiag = localTryBuildStaticObjectiveSlice(displayName, view, pilotWave, ...
  carrierFreq, sampleRate, fdRange, fdRefHz, latlonCenter, halfWidthDeg, gridSize, modelOpt, diagInfo)
%LOCALTRYBUILDSTATICOBJECTIVESLICE Build one objective slice or a NaN placeholder.

if nargin >= 11 && isfield(diagInfo, 'errorMessage') && strlength(diagInfo.errorMessage) > 0
  latGrid = linspace(latlonCenter(1) - halfWidthDeg(1), ...
    latlonCenter(1) + halfWidthDeg(1), gridSize(1));
  lonGrid = linspace(latlonCenter(2) - halfWidthDeg(2), ...
    latlonCenter(2) + halfWidthDeg(2), gridSize(2));
  sliceDiag = struct();
  sliceDiag.displayName = string(displayName) + " (skipped)";
  sliceDiag.latGrid = latGrid;
  sliceDiag.lonGrid = lonGrid;
  sliceDiag.objGrid = nan(numel(latGrid), numel(lonGrid));
  sliceDiag.objSatGrid = nan(numel(latGrid), numel(lonGrid), max(view.sceneRef.numSat, 1));
  sliceDiag.fdRefHz = fdRefHz;
  return;
end

sliceDiag = localBuildStaticObjectiveSlice(displayName, view, pilotWave, ...
  carrierFreq, sampleRate, fdRange, fdRefHz, latlonCenter, halfWidthDeg, gridSize, modelOpt);
end


function caseTruth = localBuildCaseTruth(latlonTrueDeg, fdRefTrueHz)
%LOCALBUILDCASETRUTH Build a compact truth struct for one static case.

caseTruth = struct();
caseTruth.latlonTrueDeg = latlonTrueDeg(:);
caseTruth.fdRefTrueHz = fdRefTrueHz;
end


function summaryTable = localBuildCaseSummaryTable(caseList, truthList)
%LOCALBUILDCASESUMMARYTABLE Build one summary table with per-case truth.

caseCell = num2cell(caseList);
truthCell = num2cell(truthList);
infoCell = cell(1, numel(caseCell));
for iCase = 1:numel(caseCell)
  infoCell{iCase} = summarizeDoaDopplerCase(caseCell{iCase}, truthCell{iCase});
end
summaryTable = struct2table([infoCell{:}], 'AsArray', true);
end


function weightTable = localBuildWeightSweepTable(alphaList, caseList, truth)
%LOCALBUILDWEIGHTSWEEPTABLE Build one compact weight-sweep summary table.

summaryTable = localBuildCaseSummaryTable(caseList, repmat(localBuildCaseTruth( ...
  truth.latlonTrueDeg, truth.fdRefTrueHz), 1, numel(caseList)));
weightTable = table(alphaList(:), summaryTable.displayName, summaryTable.angleErrDeg, ...
  summaryTable.fdRefErrHz, summaryTable.runTimeMs, summaryTable.isResolved, ...
  'VariableNames', {'alphaSat2', 'displayName', 'angleErrDeg', ...
  'fdRefErrHz', 'runTimeMs', 'isResolved'});
end


function sliceDiag = localBuildStaticObjectiveSlice(displayName, view, pilotWave, ...
  carrierFreq, sampleRate, fdRange, fdRefHz, latlonCenter, halfWidthDeg, gridSize, modelOpt)
%LOCALBUILDSTATICOBJECTIVESLICE Build one fixed-fd DoA objective slice.

latGrid = linspace(latlonCenter(1) - halfWidthDeg(1), ...
  latlonCenter(1) + halfWidthDeg(1), gridSize(1));
lonGrid = linspace(latlonCenter(2) - halfWidthDeg(2), ...
  latlonCenter(2) + halfWidthDeg(2), gridSize(2));
objGrid = zeros(numel(latGrid), numel(lonGrid));
objSatGrid = nan(numel(latGrid), numel(lonGrid), view.sceneRef.numSat);

evalOpt = modelOpt;
evalOpt.evalOnly = true;
evalOpt.initDoaHalfWidth = [];

for iLat = 1:numel(latGrid)
  for iLon = 1:numel(lonGrid)
    initParam = [latGrid(iLat); lonGrid(iLon); fdRefHz];
    [estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
      view.sceneRef, view.rxSigSf, pilotWave, carrierFreq, sampleRate, ...
      view.doaGrid, fdRange, 1, initParam, false, evalOpt);
    objGrid(iLat, iLon) = estResult.fval;

    objectiveSat = getDoaDopplerFieldOrDefault(estResult.aux, 'objectiveSat', []);
    if numel(objectiveSat) == view.sceneRef.numSat
      objSatGrid(iLat, iLon, :) = reshape(objectiveSat, 1, 1, []);
    end
  end
end

sliceDiag = struct();
sliceDiag.displayName = string(displayName);
sliceDiag.latGrid = latGrid;
sliceDiag.lonGrid = lonGrid;
sliceDiag.objGrid = objGrid;
sliceDiag.objSatGrid = objSatGrid;
sliceDiag.fdRefHz = fdRefHz;
end


function crbTable = localBuildCrbSummary(truth, crbSs, auxCrbSs, crbMs, auxCrbMs, crbDup, auxCrbDup)
%LOCALBUILDCRBSUMMARY Build a compact static CRB summary table.

caseName = ["SS-SF-Static"; "MS-SF-Static"; "MS-SF-Static-DupRef"];
satMode = ["single"; "multi"; "multi"];
angleStdDeg = [ ...
  projectCrbToAngleMetric(crbSs(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMs(1:2, 1:2), truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbDup(1:2, 1:2), truth.latlonTrueDeg, 'latlon')];
fdStdHz = [sqrt(max(real(crbSs(3, 3)), 0)); ...
  sqrt(max(real(crbMs(3, 3)), 0)); ...
  sqrt(max(real(crbDup(3, 3)), 0))];

crbTable = table(caseName, satMode, angleStdDeg, fdStdHz, ...
  {auxCrbSs.fdRateMode; auxCrbMs.fdRateMode; auxCrbDup.fdRateMode}, ...
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


function localPlotObjectiveSlice(sliceDiag, truthLatlon)
%LOCALPLOTOBJECTIVESLICE Plot fixed-fd objective slices for static cases.

if iscell(sliceDiag)
  sliceList = sliceDiag;
else
  sliceList = num2cell(sliceDiag);
end

figure();
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for iDiag = 1:numel(sliceList)
  diagItem = sliceList{iDiag};
  nexttile();
  imagesc(diagItem.lonGrid, diagItem.latGrid, diagItem.objGrid);
  axis xy;
  hold on;
  plot(truthLatlon(2), truthLatlon(1), 'wx', 'LineWidth', 1.8, 'MarkerSize', 9);
  hold off;
  xlabel('Longitude (deg)');
  ylabel('Latitude (deg)');
  title(sprintf('%s, fixed fd=%.1f Hz', diagItem.displayName, diagItem.fdRefHz));
  colorbar();
end
end


function tlePath = localResolveTlePath(fileName)
%LOCALRESOLVETLEPATH Resolve one TLE file from common test locations.

candidateList = [ ...
  fullfile('/tle', fileName); ...
  fullfile('tle', fileName); ...
  fullfile('test', 'data', 'tle', fileName); ...
  fullfile('..', 'data', 'tle', fileName)];

for iPath = 1:numel(candidateList)
  if exist(candidateList{iPath}, 'file')
    tlePath = candidateList{iPath};
    return;
  end
end

error('doaDopplerStatDualSatUraEci:TleFileNotFound', ...
  'Unable to locate the TLE file "%s".', fileName);
end
