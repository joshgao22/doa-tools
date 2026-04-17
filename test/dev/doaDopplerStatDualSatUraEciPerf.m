% DOADOPPLERSTATDUALSATURAECIPERF
% Static single-frame Monte Carlo performance sweep for one fixed dual-satellite
% pair. The estimation flow mirrors doaDopplerStatDualSatUraEci: same shared
% static transition bundle, same sat2-weight ablation, and the same SF static
% estimator behavior including the current DoA-anchor fallback path.
%
% Parallelization policy:
% - The Monte Carlo outer task grid (SNR x repeat) is evaluated with parfor.
% - The per-SNR CRB sweep is also evaluated with parfor.
% - Lightweight summary/plot loops remain serial to avoid parfor overhead.
% - Each Monte Carlo task uses a deterministic task-wise RNG seed so the perf
%   run remains reproducible under parfor scheduling.
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
enableWeightSweep = false;
weightSweepAlpha = [0; 0.25; 0.5; 1];
if ~enableWeightSweep
  weightSweepAlpha = zeros(0, 1);
end
staticMsHalfWidth = [0.002; 0.002];
optVerbose = false;

snrDb = -20:3:10;
numParam = numel(snrDb);
numRepeat = 20;
pwrNoiseList = pwrSource ./ (10.^(snrDb(:) / 10));

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

% Note: do not build an array-only DoA CRB here. For this script, the
% meaningful theoretical baseline is the pilot/static DoA-Doppler CRB,
% because the estimated global DoA is coupled with reference-link Doppler
% in the actual single-frame signal model.

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
caseName = ["SS-SF-DoA", "MS-SF-DoA", "SS-SF-Static", "MS-SF-Static"];
idxSsDoa = 1;
idxMsDoa = 2;
idxSsStat = 3;
idxMsStat = 4;
idxWeightStart = [];
if enableWeightSweep
  idxWeightStart = numel(caseName) + 1;
  caseName = [caseName, compose('MS-SF-Static-W%.2f', weightSweepAlpha.')];
end
numCase = numel(caseName);

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

%% Monte Carlo loop (parallel over the full SNR x repeat task grid)
angleErrDeg = nan(numRepeat, numParam, numCase);
fdErrHz = nan(numRepeat, numParam, numCase);
isResolved = false(numRepeat, numParam, numCase);

[repeatGrid, snrGrid] = ndgrid(1:numRepeat, 1:numParam);
repeatTaskIdx = repeatGrid(:);
snrTaskIdx = snrGrid(:);
numTask = numel(repeatTaskIdx);

taskGrid = repmat(struct('taskIndex', 0, 'snrIndex', 0, 'repeatIndex', 0, 'taskSeed', 0), numTask, 1);
for iTask = 1:numTask
  taskGrid(iTask).taskIndex = iTask;
  taskGrid(iTask).snrIndex = snrTaskIdx(iTask);
  taskGrid(iTask).repeatIndex = repeatTaskIdx(iTask);
  taskGrid(iTask).taskSeed = baseSeed + iTask - 1;
end

sharedData = struct();
sharedData.pwrNoiseList = pwrNoiseList;
sharedData.steeringInfo = steeringInfo;
sharedData.pilotWave = pilotWave;
sharedData.linkParam = linkParam;
sharedData.carrierFreq = carrierFreq;
sharedData.sampleRate = waveInfo.sampleRate;
sharedData.pathGain = pathGain;
sharedData.snapOpt = snapOpt;
sharedData.viewRefTemplate = viewRefTemplate;
sharedData.viewOtherTemplate = viewOtherTemplate;
sharedData.viewMsTemplate = viewMsTemplate;
sharedData.refSatIdxLocal = refSatIdxLocal;
sharedData.otherSatIdxLocal = otherSatIdxLocal;
sharedData.wavelen = wavelen;
sharedData.fdRange = fdRange;
sharedData.optVerbose = optVerbose;
sharedData.truth = truth;
sharedData.otherSatIdxGlobal = otherSatIdxGlobal;
sharedData.doaOnlyOpt = doaOnlyOpt;
sharedData.staticOptBase = staticOptBase;
sharedData.weightSweepAlpha = weightSweepAlpha;
sharedData.staticMsHalfWidth = staticMsHalfWidth;
sharedData.numCase = numCase;

repoRoot = localGetRepoRoot();
checkpointMeta = struct();
checkpointMeta.snrDb = snrDb(:).';
checkpointMeta.numRepeat = numRepeat;
checkpointMeta.numCase = numCase;
checkpointMeta.baseSeed = baseSeed;
checkpointMeta.selectedSatIdxGlobal = truth.selectedSatIdxGlobal(:).';
checkpointMeta.fdRange = fdRange;
checkpointMeta.caseName = caseName(:).';
checkpointMeta.enableWeightSweep = enableWeightSweep;

checkpointOpt = struct();
checkpointOpt.runName = "doaDopplerStatDualSatUraEciPerf";
checkpointOpt.runKey = localBuildCheckpointRunKey(baseSeed, snrDb, numRepeat, numCase, truth.selectedSatIdxGlobal);
checkpointOpt.outputRoot = fullfile(repoRoot, "out", "checkpoint_tmp");
checkpointOpt.useParfor = true;
checkpointOpt.resume = true;
checkpointOpt.meta = checkpointMeta;

numProgressStep = numTask + numParam;
progressbar('reset', numProgressStep);
checkpointOpt.progressFcn = @(step) progressbar('advance', step);

fprintf('\nRunning Monte Carlo task grid with checkpoint runner: %d SNR x %d repeats = %d tasks.%s', ...
  numParam, numRepeat, numTask, newline);
runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, @localPerfTaskRunner, checkpointOpt);
if ~runState.isComplete
  error('doaDopplerStatDualSatUraEciPerf:CheckpointIncomplete', ...
    'Checkpoint runner returned an incomplete task grid.');
end

angleErrTask = nan(numTask, numCase);
fdErrTask = nan(numTask, numCase);
isResolvedTask = false(numTask, numCase);
for iTask = 1:numTask
  taskOut = runState.resultCell{iTask};
  if isempty(taskOut)
    continue;
  end
  angleErrTask(iTask, :) = taskOut.angleVec;
  fdErrTask(iTask, :) = taskOut.fdVec;
  isResolvedTask(iTask, :) = taskOut.resolvedVec;
end

for iCase = 1:numCase
  angleErrDeg(:, :, iCase) = reshape(angleErrTask(:, iCase), numRepeat, numParam);
  fdErrHz(:, :, iCase) = reshape(fdErrTask(:, iCase), numRepeat, numParam);
  isResolved(:, :, iCase) = reshape(isResolvedTask(:, iCase), numRepeat, numParam);
end

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

weightSummary = table();
if enableWeightSweep
  weightSummary = localBuildWeightSummaryTable(snrDb, weightSweepAlpha, ...
    rmseAngleDeg(:, idxWeightStart:end), rmseFdHz(:, idxWeightStart:end), ...
    resolveRate(:, idxWeightStart:end));
end

snrSelectIdx = numParam;
caseSummaryTable = localBuildCaseSummaryTable(caseName, snrDb(snrSelectIdx), ...
  rmseAngleDeg(snrSelectIdx, :), rmseFdHz(snrSelectIdx, :), ...
  resolveRate(snrSelectIdx, :), p95AngleDeg(snrSelectIdx, :), p95FdHz(snrSelectIdx, :));

%% CRB curves (parallel over SNR)
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

fprintf('Running CRB SNR sweep with parfor: %d points.%s', numParam, newline);
crbProgressQueue = parallel.pool.DataQueue;
afterEach(crbProgressQueue, @(~) progressbar('advance'));
parfor iSnr = 1:numParam
  pwrNoise = pwrNoiseList(iSnr);

  [crbDdSingle, ~] = crbPilotSfDoaDoppler(sceneRefOnly, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);
  [crbDdJoint, ~] = crbPilotSfDoaDoppler(scene, pilotWave, ...
    carrierFreq, waveInfo.sampleRate, truthLatlon, truth.fdRefTrueHz, 1, pwrNoise, crbDdOpt);

  crbDdSingleDeg(iSnr) = projectCrbToAngleMetric(crbDdSingle(1:2, 1:2), truthLatlon, 'latlon');
  crbDdJointDeg(iSnr) = projectCrbToAngleMetric(crbDdJoint(1:2, 1:2), truthLatlon, 'latlon');
  crbDdSingleFd(iSnr) = sqrt(crbDdSingle(3, 3));
  crbDdJointFd(iSnr) = sqrt(crbDdJoint(3, 3));
  send(crbProgressQueue, iSnr);
end
progressbar('end');

% Reuse the pilot/static CRB curves when plotting the DoA-only cases as a
% conservative and model-consistent lower bound in the same signal model.
crbDoaOnlySingleDeg = crbDdSingleDeg;
crbDoaOnlyJointDeg = crbDdJointDeg;

%% Summary tables
if enableWeightSweep
  fprintf('\n========== Weight-sweep best summary ==========%s', newline);
  disp(weightSummary);
end

fprintf('\n========== Case summary at highest SNR ==========%s', newline);
disp(caseSummaryTable);

%% Plots
figure();
subplot(1, 2, 1);
semilogy(snrDb, rmseAngleDeg(:, idxSsDoa), '-o', ...
  snrDb, rmseAngleDeg(:, idxMsDoa), '-s', ...
  snrDb, rmseAngleDeg(:, idxSsStat), '--o', ...
  snrDb, rmseAngleDeg(:, idxMsStat), '--s', ...
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
  'SS/MF-consistent DoA lower bound (single)', ...
  'SS/MF-consistent DoA lower bound (joint)', ...
  'SS-SF-Static CRB', ...
  'MS-SF-Static CRB'}, ...
  'Location', 'southwest');

subplot(1, 2, 2);
semilogy(snrDb, rmseFdHz(:, idxSsStat), '--o', ...
  snrDb, rmseFdHz(:, idxMsStat), '--s', ...
  snrDb, crbDdSingleFd, '-.x', ...
  snrDb, crbDdJointFd, ':x');
grid on;
xlabel('SNR (dB)');
ylabel('Reference Doppler RMSE / CRB (Hz)');
title('Static SF Reference Doppler Error');
legend({ ...
  'SS-SF-Static', ...
  'MS-SF-Static', ...
  'SS-SF-Static CRB', ...
  'MS-SF-Static CRB'}, ...
  'Location', 'southwest');

figure();
plot(snrDb, resolveRate(:, idxSsStat), '-o', ...
  snrDb, resolveRate(:, idxMsStat), '-s');
grid on;
ylim([0, 1.05]);
xlabel('SNR (dB)');
ylabel('Resolve rate');
title('Static estimator resolve rate');
legend({'SS-SF-Static', 'MS-SF-Static'}, ...
  'Location', 'southeast');

if enableWeightSweep
  figure();
  plot(snrDb, rmseAngleDeg(:, idxWeightStart:end), '-o');
  grid on;
  xlabel('SNR (dB)');
  ylabel('Angle RMSE (deg)');
  title('Static multi-satellite weight sweep');
  legend(compose('alpha = %.2f', weightSweepAlpha), 'Location', 'southwest');
end

%% Optional snapshot save
checkpointRunDir = runState.runDir;
clear sharedData taskGrid checkpointMeta checkpointOpt runState;
saveExpSnapshot("doaDopplerStatDualSatUraEciPerf");
if isfolder(checkpointRunDir)
  rmdir(checkpointRunDir, 's');
end

%% Local functions
function taskOut = localPerfTaskRunner(taskInfo, sharedData)
%LOCALPERFTASKRUNNER Thin adapter from one outer task to one MC result struct.

[angleVec, fdVec, resolvedVec] = localRunMonteCarloTask( ...
  taskInfo.taskSeed, sharedData.pwrNoiseList(taskInfo.snrIndex), ...
  sharedData.steeringInfo, sharedData.pilotWave, sharedData.linkParam, ...
  sharedData.carrierFreq, sharedData.sampleRate, sharedData.pathGain, ...
  sharedData.snapOpt, sharedData.viewRefTemplate, sharedData.viewOtherTemplate, ...
  sharedData.viewMsTemplate, sharedData.refSatIdxLocal, sharedData.otherSatIdxLocal, ...
  sharedData.wavelen, sharedData.fdRange, sharedData.optVerbose, ...
  sharedData.truth, sharedData.otherSatIdxGlobal, sharedData.doaOnlyOpt, ...
  sharedData.staticOptBase, sharedData.weightSweepAlpha, ...
  sharedData.staticMsHalfWidth, sharedData.numCase);

taskOut = struct();
taskOut.angleVec = angleVec;
taskOut.fdVec = fdVec;
taskOut.resolvedVec = resolvedVec;
end


function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Return the repository root for temporary checkpoint output.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
end


function runKey = localBuildCheckpointRunKey(baseSeed, snrDb, numRepeat, numCase, selectedSatIdxGlobal)
%LOCALBUILDCHECKPOINTRUNKEY Build one stable checkpoint directory key.

satToken = sprintf('sat%d', selectedSatIdxGlobal(1));
for iSat = 2:numel(selectedSatIdxGlobal)
  satToken = sprintf('%s_%d', satToken, selectedSatIdxGlobal(iSat));
end
runKey = string(sprintf('seed%d_%s_snr%dto%d_n%d_rep%d_case%d', ...
  baseSeed, satToken, snrDb(1), snrDb(end), numel(snrDb), numRepeat, numCase));
end


function [angleVec, fdVec, resolvedVec] = localRunMonteCarloTask( ...
  taskSeed, pwrNoise, steeringInfo, pilotWave, linkParam, carrierFreq, sampleRate, ...
  pathGain, snapOpt, viewRefTemplate, viewOtherTemplate, viewMsTemplate, ...
  refSatIdxLocal, otherSatIdxLocal, wavelen, fdRange, optVerbose, truth, ...
  otherSatIdxGlobal, doaOnlyOpt, staticOptBase, weightSweepAlpha, ...
  staticMsHalfWidth, numCase)
%LOCALRUNMONTECARLOTASK Run one independent SNR-repeat Monte Carlo task.

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

caseList = [ ...
  caseBundle.caseRefDoa, ...
  caseBundle.caseMsDoa, ...
  caseBundle.caseStaticRefOnly, ...
  caseBundle.caseStaticMs];
truthFdList = [NaN, NaN, truth.fdSatTrueHz(refSatIdxLocal), truth.fdRefTrueHz];
if ~isempty(caseBundle.weightCase)
  caseList = [caseList, caseBundle.weightCase];
  truthFdList = [truthFdList, repmat(truth.fdRefTrueHz, 1, numel(caseBundle.weightCase))];
end

angleVec = nan(1, numCase);
fdVec = nan(1, numCase);
resolvedVec = false(1, numCase);
for iCase = 1:numCase
  [angleVec(iCase), fdVec(iCase), resolvedVec(iCase)] = localExtractCaseMetric( ...
    caseList(iCase), truth.latlonTrueDeg, truthFdList(iCase));
end
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
