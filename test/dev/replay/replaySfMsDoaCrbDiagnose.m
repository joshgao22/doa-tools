% replaySfMsDoaCrbDiagnose
% Purpose: diagnose the MS-SF-DoA MLE-vs-CRB gap in the static anchor.
% This replay keeps scanSfStaticMleCrbConsistency as the paper-facing scan
% and runs a small MC plus local probes to separate CRB / metric mismatch
% from DoA-only MLE solver or local-basin issues.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace. saveSnapshot=true saves only
% replayData via saveExpSnapshot. checkpointEnable=true writes resumable
% per-task artifacts under tmp/ and cleans them after a successful run.
% Telegram notice is best-effort only.

clear; close all; clc;

%% Replay configuration

replayName = "replaySfMsDoaCrbDiagnose";
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
notifyTelegramEnable = true;
optVerbose = false;

baseSeed = 253;
numRepeat = 500;
snrDbList = -20:5:10;
trimEnable = true;
trimNormCap = 5;
topTailCountPerGroup = 5;
probeSeedCountPerSnr = 4;
hessianStepDeg = [1e-4; 1e-4];
hessianMeanRepeat = 8;
alphaSat2List = [0; 0.1; 0.25; 0.5; 1; 2; 5];
numUsr = 1;
numSym = 512;
sampleRate = 512e6;
symbolRate = 128e6;
carrierFreq = 11.7e9;
gridSize = [50, 50];
searchMarginDeg = 5;

seedList = reshape(double(baseSeed + (0:(numRepeat - 1))), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
hessianStepDeg = reshape(double(hessianStepDeg), [], 1);
numRepeat = numel(seedList);

config = struct();
config.replayName = string(replayName);
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.snrDbList = snrDbList;
config.trimEnable = logical(trimEnable);
config.trimNormCap = trimNormCap;
config.topTailCountPerGroup = topTailCountPerGroup;
config.probeSeedCountPerSnr = probeSeedCountPerSnr;
config.hessianStepDeg = hessianStepDeg;
config.hessianMeanRepeat = hessianMeanRepeat;
config.alphaSat2List = reshape(double(alphaSat2List), [], 1);
config.numUsr = numUsr;
config.numSym = numSym;
config.sampleRate = sampleRate;
config.symbolRate = symbolRate;
config.carrierFreq = carrierFreq;
config.gridSize = reshape(double(gridSize), 1, []);
config.searchMarginDeg = searchMarginDeg;
config.saveSnapshot = logical(saveSnapshot);
config.checkpointEnable = logical(checkpointEnable);
config.checkpointResume = logical(checkpointResume);
config.notifyTelegramEnable = logical(notifyTelegramEnable);
config.optVerbose = logical(optVerbose);

runTic = tic;
replayData = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and flow options

  config.runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
  taskList = localBuildTaskList(config.snrDbList, config.seedList, config.probeSeedCountPerSnr);
  config.numTask = numel(taskList);
  config.useParfor = localCanUseParfor(config.numTask);
  checkpointOpt = localBuildCheckpointOpt(replayName, config, config.useParfor);
  checkpointRunDir = string(localGetFieldOrDefault(checkpointOpt, 'runDir', ""));
  runtime = localBuildReplayRuntime(config);
  printMfReplayHeader(char(replayName), config, checkpointRunDir);
  localPrintReplayConfig(config, runtime.truth);

  %% Run replay batch

  fprintf('Running MS-SF-DoA CRB diagnostic tasks: %d tasks.\n', numel(taskList));
  [taskOutCell, runState] = localRunTaskBatch(config, taskList, runtime, replayName, checkpointOpt);

  caseTable = localCollectTaskTable(taskOutCell, 'caseTable');
  truthInitPairTable = localCollectTaskTable(taskOutCell, 'truthInitPairTable');
  objectiveProbeTable = localCollectTaskTable(taskOutCell, 'objectiveProbeTable');
  hessianProbeTable = localCollectTaskTable(taskOutCell, 'hessianProbeTable');
  noiseFreeCaseTable = localCollectTaskTable(taskOutCell, 'noiseFreeCaseTable');
  alphaSweepCaseTable = localCollectTaskTable(taskOutCell, 'alphaSweepCaseTable');
  crbTable = localBuildCrbTable(runtime, config);
  perfTable = localAttachCrbMetrics(caseTable, crbTable);
  aggregateTable = localBuildAggregateTable(perfTable, config);
  topTailTable = localBuildTopTailTable(perfTable, config.topTailCountPerGroup);
  crbCompareTable = localBuildCrbCompareTable(crbTable);
  truthInitDiffTable = localBuildTruthInitDiffTable(truthInitPairTable, crbTable);
  doaGainTable = localBuildDoaGainTable(aggregateTable, crbTable);
  otherDoaCrbAuditTable = localBuildOtherDoaCrbAuditTable(perfTable, crbTable, runtime.truth, config);
  normalizationAuditTable = localBuildNormalizationAuditTable(perfTable, crbTable, hessianProbeTable, runtime, config);
  alphaWeightSweepTable = localBuildAlphaWeightSweepTable(alphaSweepCaseTable, crbTable, aggregateTable);
  empiricalGainTable = localBuildEmpiricalGainTable(aggregateTable, crbTable);
  ampAwareCrbAuditTable = localBuildAmpAwareCrbAuditTable(aggregateTable, crbTable);
  hessianAggregateTable = localBuildHessianAggregateTable(hessianProbeTable);
  plotData = localBuildPlotData(aggregateTable, crbTable, crbCompareTable, truthInitDiffTable, ...
    doaGainTable, otherDoaCrbAuditTable, normalizationAuditTable, alphaWeightSweepTable, ...
    empiricalGainTable, ampAwareCrbAuditTable, noiseFreeCaseTable);

  %% Data storage

  replayData = struct();
  replayData.replayName = string(replayName);
  replayData.runKey = config.runKey;
  replayData.utcRun = datetime('now', 'TimeZone', 'local');
  replayData.config = config;
  replayData.truth = runtime.truth;
  replayData.caseTable = caseTable;
  replayData.truthInitPairTable = truthInitPairTable;
  replayData.objectiveProbeTable = objectiveProbeTable;
  replayData.hessianProbeTable = hessianProbeTable;
  replayData.noiseFreeCaseTable = noiseFreeCaseTable;
  replayData.crbTable = crbTable;
  replayData.crbCompareTable = crbCompareTable;
  replayData.perfTable = perfTable;
  replayData.aggregateTable = aggregateTable;
  replayData.truthInitDiffTable = truthInitDiffTable;
  replayData.doaGainTable = doaGainTable;
  replayData.otherDoaCrbAuditTable = otherDoaCrbAuditTable;
  replayData.normalizationAuditTable = normalizationAuditTable;
  replayData.alphaSweepCaseTable = alphaSweepCaseTable;
  replayData.alphaWeightSweepTable = alphaWeightSweepTable;
  replayData.empiricalGainTable = empiricalGainTable;
  replayData.ampAwareCrbAuditTable = ampAwareCrbAuditTable;
  replayData.hessianAggregateTable = hessianAggregateTable;
  replayData.topTailTable = topTailTable;
  replayData.plotData = plotData;
  replayData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  replayData.elapsedSec = toc(runTic);
  replayData = finalizeMfReplayResult(replayData, checkpointRunDir);

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'replayData'}}, ...
      'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
    replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
  else
    replayData.snapshotFile = "";
  end

  %% Summary output and plotting

  replayData = localValidateReplayDataForSummary(replayData);
  replayNameForReport = string(localGetFieldOrDefault(replayData, 'replayName', "replaySfMsDoaCrbDiagnose"));
  replayConfigForReport = localGetFieldOrDefault(replayData, 'config', struct());
  replaySnapshotFile = localGetFieldOrDefault(replayData, 'snapshotFile', "");
  replayElapsedSec = localGetFieldOrDefault(replayData, 'elapsedSec', NaN);

  printMfReplaySection('MS-SF-DoA CRB diagnostic aggregate', replayData.aggregateTable);
  printMfReplaySection('DoA gain realization summary', replayData.doaGainTable);
  printMfReplaySection('Empirical joint-gain prediction', replayData.empiricalGainTable);
  printMfReplaySection('Pilot-model amp-aware CRB audit', replayData.ampAwareCrbAuditTable);
  printMfReplaySection('Other-sat DoA CRB audit', replayData.otherDoaCrbAuditTable);
  printMfReplaySection('Per-sat normalization audit', replayData.normalizationAuditTable);
  printMfReplaySection('CRB definition comparison', replayData.crbCompareTable);
  printMfReplaySection('Baseline vs truth-init compact diff', replayData.truthInitDiffTable);
  if height(replayData.alphaWeightSweepTable) > 0
    printMfReplaySection('MS-SF-DoA alpha weight sweep', replayData.alphaWeightSweepTable);
  end
  if isfield(replayData, 'hessianAggregateTable') && height(replayData.hessianAggregateTable) > 0
    printMfReplaySection('Per-sat numerical Hessian aggregate (mean-noisy only)', ...
      localFilterHessianAggregateForPrint(replayData.hessianAggregateTable));
  end
  if height(replayData.noiseFreeCaseTable) > 0
    printMfReplaySection('Noise-free DoA sanity preview', ...
      localProbePreviewBySnrSeed(replayData.noiseFreeCaseTable, 1));
  end
  if height(replayData.objectiveProbeTable) > 0
    printMfReplaySection('Truth/final objective probe preview', ...
      localProbePreviewBySnrSeed(replayData.objectiveProbeTable, 1));
  end
  if height(replayData.topTailTable) > 0
    printMfReplaySection('Top normalized-error tail preview', ...
      localTablePreview(replayData.topTailTable, 24));
  end
  localPrintReplaySummaryConfig(replayData);
  replayData.plotData = localPlotReplay(replayData);

  notifyMfReplayStatus(struct( ...
    'replayName', replayNameForReport, ...
    'statusText', "DONE", ...
    'subtitleText', "MS-SF-DoA MLE/CRB diagnostic replay", ...
    'config', replayConfigForReport, ...
    'snapshotFile', replaySnapshotFile, ...
    'elapsedSec', replayElapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(replayData), ...
    'commentLineList', [ ...
      "MS-SF-DoA diagnostic replay is closed as a CRB-scale audit."; ...
      "Use the saved replayData tables to formalize the amp-aware deterministic CRB before changing scan or estimator defaults."]));

catch ME
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'subtitleText', "MS-SF-DoA MLE/CRB diagnostic replay", ...
    'config', config, ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function checkpointOpt = localBuildCheckpointOpt(replayName, config, useParfor)
%LOCALBUILDCHECKPOINTOPT Build stable per-task checkpoint options for this replay.

checkpointMeta = struct( ...
  'snrDbList', reshape(config.snrDbList, 1, []), ...
  'seedList', reshape(config.seedList, 1, []), ...
  'probeSeedCountPerSnr', config.probeSeedCountPerSnr, ...
  'hessianStepDeg', reshape(config.hessianStepDeg, 1, []), ...
  'hessianMeanRepeat', config.hessianMeanRepeat, ...
  'alphaSat2List', reshape(config.alphaSat2List, 1, []), ...
  'numSym', config.numSym, ...
  'sampleRate', config.sampleRate, ...
  'symbolRate', config.symbolRate, ...
  'carrierFreq', config.carrierFreq, ...
  'gridSize', reshape(config.gridSize, 1, []), ...
  'searchMarginDeg', config.searchMarginDeg);
checkpointOpt = buildMfReplayCheckpointOpt(replayName, config, struct( ...
  'runKey', localBuildCheckpointRunKey(config), ...
  'meta', checkpointMeta, ...
  'useParfor', logical(useParfor)));
end

function runKey = localBuildCheckpointRunKey(config)
%LOCALBUILDCHECKPOINTRUNKEY Build a compact stable key for interrupted replays.

seedList = reshape(double(config.seedList), 1, []);
snrList = reshape(double(config.snrDbList), 1, []);
signatureHash = localCompactStringHash(localBuildCheckpointSignature(config));
runKey = sprintf('snr%.6gto%.6g_seed%dto%d_rep%d_%s', ...
  snrList(1), snrList(end), seedList(1), seedList(end), numel(seedList), signatureHash);
runKey = localPathSafeToken(runKey);
end

function signature = localBuildCheckpointSignature(config)
%LOCALBUILDCHECKPOINTSIGNATURE Build the semantic signature behind runKey.

signature = strjoin(string({ ...
  sprintf('snr%s', mat2str(reshape(double(config.snrDbList), 1, []), 8)), ...
  sprintf('seed%s', mat2str(reshape(double(config.seedList), 1, []), 8)), ...
  sprintf('probe%d', config.probeSeedCountPerSnr), ...
  sprintf('hStep%s', mat2str(reshape(double(config.hessianStepDeg), 1, []), 8)), ...
  sprintf('hRepeat%d', config.hessianMeanRepeat), ...
  sprintf('alpha%s', mat2str(reshape(double(config.alphaSat2List), 1, []), 8)), ...
  sprintf('numSym%d', config.numSym), ...
  sprintf('sample%.8g', config.sampleRate), ...
  sprintf('symbol%.8g', config.symbolRate), ...
  sprintf('carrier%.8g', config.carrierFreq), ...
  sprintf('grid%s', mat2str(reshape(double(config.gridSize), 1, []), 8)), ...
  sprintf('margin%.8g', config.searchMarginDeg) ...
  }), '|');
end

function hashText = localCompactStringHash(textValue)
%LOCALCOMPACTSTRINGHASH Return a short deterministic hash for path-safe names.

byteValues = double(char(textValue));
hashValue = 5381;
modValue = 2147483647;
for iByte = 1:numel(byteValues)
  hashValue = mod(hashValue * 33 + byteValues(iByte), modValue);
end
hashText = sprintf('%08x', uint32(hashValue));
end

function token = localPathSafeToken(value)
%LOCALPATHSAFETOKEN Convert a short semantic token into a path-safe string.

token = string(value);
token = replace(token, '.', 'p');
token = replace(token, '-', 'm');
token = replace(token, ' ', '');
token = regexprep(char(token), '[^A-Za-z0-9_]', '');
token = string(token);
end

function [taskOutCell, runState] = localRunTaskBatch(config, taskList, runtime, replayName, checkpointOpt)
%LOCALRUNTASKBATCH Run SNR/seed tasks with checkpoint and progress support.

numTask = numel(taskList);
runState = struct();
progressTracker = localCreateProgressTracker(numTask);
progressCleanup = onCleanup(@() localCloseProgressTracker(progressTracker)); %#ok<NASGU>
try
  if logical(localGetFieldOrDefault(config, 'checkpointEnable', false))
    taskDir = fullfile(checkpointOpt.runDir, 'task');
    numDoneTask = localCountCheckpointTaskFile(taskDir, numTask);
    localAdvanceProgressByStep(progressTracker, numDoneTask);
    sharedData = struct('runtime', runtime);
    checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker);
    runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, @localCheckpointTaskRunner, checkpointRunnerOpt);
    taskOutCell = runState.resultCell;
  else
    taskOutCell = cell(numTask, 1);
    useParfor = logical(localGetFieldOrDefault(config, 'useParfor', false));
    if useParfor && isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
      progressQueue = progressTracker.queue;
      parfor iTask = 1:numTask
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), runtime);
        send(progressQueue, 1);
      end
    elseif useParfor
      parfor iTask = 1:numTask
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), runtime);
      end
    else
      for iTask = 1:numTask
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), runtime);
        localAdvanceProgress(progressTracker);
      end
    end
  end
  localCloseProgressTracker(progressTracker);
catch ME
  localCloseProgressTracker(progressTracker);
  rethrow(ME);
end
end

function taskOut = localCheckpointTaskRunner(taskInfo, sharedData)
%LOCALCHECKPOINTTASKRUNNER Run one checkpointed task.

taskOut = localRunOneTask(taskInfo, sharedData.runtime);
end

function checkpointRunnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only fields accepted by checkpoint runner.

checkpointRunnerOpt = struct();
checkpointRunnerOpt.runName = checkpointOpt.runName;
checkpointRunnerOpt.runKey = checkpointOpt.runKey;
checkpointRunnerOpt.outputRoot = checkpointOpt.outputRoot;
checkpointRunnerOpt.useParfor = logical(localGetFieldOrDefault(checkpointOpt, 'useParfor', false));
checkpointRunnerOpt.resume = checkpointOpt.resume;
checkpointRunnerOpt.meta = checkpointOpt.meta;
checkpointRunnerOpt.progressFcn = @(step) localAdvanceProgressByStep(progressTracker, step);
checkpointRunnerOpt.cleanupOnSuccess = false;
checkpointRunnerOpt.cleanupOpt = struct();
end

function numDone = localCountCheckpointTaskFile(taskDir, numTask)
%LOCALCOUNTCHECKPOINTTASKFILE Count existing completed checkpoint task files.

numDone = 0;
if ~isfolder(taskDir)
  return;
end
for iTask = 1:numTask
  if isfile(fullfile(taskDir, sprintf('task_%06d.mat', iTask)))
    numDone = numDone + 1;
  end
end
end

function progressTracker = localCreateProgressTracker(totalCount)
%LOCALCREATEPROGRESSTRACKER Create a best-effort replay progress tracker.

progressTracker = struct('enabled', false, 'queue', []);
if totalCount <= 0
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  progressTracker.enabled = true;
  try
    progressTracker.queue = parallel.pool.DataQueue;
    afterEach(progressTracker.queue, @(~) progressbar('advance'));
  catch
    progressTracker.queue = [];
  end
catch
  progressTracker.enabled = false;
  progressTracker.queue = [];
end
end

function localAdvanceProgressByStep(progressTracker, step)
%LOCALADVANCEPROGRESSBYSTEP Advance progressbar by a task count.

for iStep = 1:step
  localAdvanceProgress(progressTracker);
end
end

function localAdvanceProgress(progressTracker)
%LOCALADVANCEPROGRESS Advance progress when progressbar is available.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('advance');
catch
end
end

function localCloseProgressTracker(progressTracker)
%LOCALCLOSEPROGRESSTRACKER Close a best-effort replay progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function replayData = localValidateReplayDataForSummary(replayData)
%LOCALVALIDATEREPLAYDATAFORSUMMARY Validate replayData contents for summary reruns.

requiredTable = {'aggregateTable', 'crbCompareTable', 'truthInitDiffTable', ...
  'doaGainTable', 'otherDoaCrbAuditTable', 'normalizationAuditTable', ...
  'empiricalGainTable', 'ampAwareCrbAuditTable', 'hessianAggregateTable'};
for iField = 1:numel(requiredTable)
  fieldName = requiredTable{iField};
  if ~isfield(replayData, fieldName) || ~istable(replayData.(fieldName))
    error('replaySfMsDoaCrbDiagnose:MissingReplayTable', ...
      'replayData.%s is missing. Store summary inputs before saving.', fieldName);
  end
end
end

function taskList = localBuildTaskList(snrDbList, seedList, probeSeedCountPerSnr)
%LOCALBUILDTASKLIST Build one flat SNR-repeat task list.

taskTemplate = struct('taskId', 0, 'snrDb', NaN, 'snrIndex', 0, ...
  'repeatIndex', 0, 'taskSeed', 0, 'probeEnable', false);
taskList = repmat(taskTemplate, numel(snrDbList) * numel(seedList), 1);
taskId = 0;
for iRepeat = 1:numel(seedList)
  for iSnr = 1:numel(snrDbList)
    taskId = taskId + 1;
    taskList(taskId).taskId = taskId;
    taskList(taskId).snrDb = snrDbList(iSnr);
    taskList(taskId).snrIndex = iSnr;
    taskList(taskId).repeatIndex = iRepeat;
    taskList(taskId).taskSeed = seedList(iRepeat);
    taskList(taskId).probeEnable = iRepeat <= probeSeedCountPerSnr;
  end
end
end

function runtime = localBuildReplayRuntime(config)
%LOCALBUILDREPLAYRUNTIME Build the shared static scene and waveform context.

rng(config.baseSeed);
lightSpeed = 299792458;
osf = config.sampleRate / config.symbolRate;
wavelen = lightSpeed / config.carrierFreq;
pwrSource = 1;
E = referenceEllipsoid('sphere');
usrLla = [37.78; 36.59; 0];
truthLatlon = usrLla(1:2, 1);
searchRange = [truthLatlon(1) - config.searchMarginDeg, truthLatlon(1) + config.searchMarginDeg; ...
               truthLatlon(2) - config.searchMarginDeg, truthLatlon(2) + config.searchMarginDeg];
utc = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);

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
  error('replaySfMsDoaCrbDiagnose:InvalidNumSat', ...
    'This replay expects exactly two selected satellites.');
end
otherSatIdxLocal = 3 - refSatIdxLocal;

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
truth.fdOtherTrueHz = truth.fdSatTrueHz(otherSatIdxLocal);
truth.otherSatIdxGlobal = satIdx(otherSatIdxLocal);
truth.otherSatIdxLocal = otherSatIdxLocal;
truth.localDoaRef = reshape(scene.localDoa(:, refSatIdxLocal), 2, 1);
truth.refWeight = scene.ref.weight(:);
truth.refStateSource = string(refState.source);
truth.pickAux = satPickAux;

sceneRefOnly = selectSatScene(scene, refSatIdxLocal);
sceneOtherOnly = selectSatScene(scene, otherSatIdxLocal);
viewRefTemplate = buildDoaDopplerEstView(sceneRefOnly, [], config.gridSize, searchRange, E);
viewOtherTemplate = buildDoaDopplerEstView(sceneOtherOnly, [], config.gridSize, searchRange, E);
viewMsTemplate = buildDoaDopplerEstView(scene, [], config.gridSize, searchRange, E);

[pilotSym, ~] = genPilotSymbol(config.numUsr, config.numSym, 'pn', pwrSource);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, config.symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
pathGain = ones(scene.numSat, scene.numUser);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

[~, ~, cleanSigPilotModel, ~, ~] = genMultiSatSnapshots( ...
  steeringInfo, pilotWave, linkParam, config.carrierFreq, waveInfo.sampleRate, [], pathGain, snapOpt);
effectivePilotGainAbs = localEstimateTruthPilotGainAbs(viewMsTemplate, cleanSigPilotModel, ...
  pilotWave, wavelen, truthLatlon);
truth.effectivePilotGainAbs = effectivePilotGainAbs(:);
truth.effectivePilotGainAbsRef = effectivePilotGainAbs(refSatIdxLocal);
truth.effectivePilotGainAbsOther = effectivePilotGainAbs(otherSatIdxLocal);

runtime = struct();
runtime.config = config;
runtime.scene = scene;
runtime.sceneRefOnly = sceneRefOnly;
runtime.sceneOtherOnly = sceneOtherOnly;
runtime.truth = truth;
runtime.truthDopplerState = truthDopplerState;
runtime.linkParam = linkParam;
runtime.steeringInfo = steeringInfo;
runtime.pilotWave = pilotWave;
runtime.sampleRate = waveInfo.sampleRate;
runtime.carrierFreq = config.carrierFreq;
runtime.pwrSource = pwrSource;
runtime.pathGain = pathGain;
runtime.snapOpt = snapOpt;
runtime.wavelen = wavelen;
runtime.refSatIdxLocal = refSatIdxLocal;
runtime.otherSatIdxLocal = otherSatIdxLocal;
runtime.viewRefTemplate = viewRefTemplate;
runtime.viewOtherTemplate = viewOtherTemplate;
runtime.viewMsTemplate = viewMsTemplate;
runtime.doaOnlyOpt = doaOnlyOpt;
end

function taskOut = localRunOneTask(task, runtime)
%LOCALRUNONETASK Run one independent static DoA-only diagnostic task.

config = runtime.config;
rng(task.taskSeed, 'twister');
pwrNoise = runtime.pwrSource / (10^(task.snrDb / 10));
[rxSig, ~, cleanSig, ~, ~] = genMultiSatSnapshots( ...
  runtime.steeringInfo, runtime.pilotWave, runtime.linkParam, ...
  runtime.carrierFreq, runtime.sampleRate, pwrNoise, runtime.pathGain, ...
  runtime.snapOpt);

viewRef = runtime.viewRefTemplate;
viewRef.rxSigSf = selectRxSigBySat(rxSig, runtime.refSatIdxLocal, 'singleFrame');
viewOther = runtime.viewOtherTemplate;
viewOther.rxSigSf = selectRxSigBySat(rxSig, runtime.otherSatIdxLocal, 'singleFrame');
viewMs = runtime.viewMsTemplate;
viewMs.rxSigSf = rxSig;

caseList = [ ...
  localRunDoaOnlyCase("SS-SF-DoA", "single", viewRef, runtime.wavelen, ...
    runtime.pilotWave, config.optVerbose, runtime.doaOnlyOpt, []), ...
  localRunDoaOnlyCase("Other-SF-DoA", "other", viewOther, runtime.wavelen, ...
    runtime.pilotWave, config.optVerbose, runtime.doaOnlyOpt, []), ...
  localRunDoaOnlyCase("MS-SF-DoA", "multi", viewMs, runtime.wavelen, ...
    runtime.pilotWave, config.optVerbose, runtime.doaOnlyOpt, [])];
truthInitCase = localRunDoaOnlyCase("MS-SF-DoA-TruthInit", "multi", viewMs, ...
  runtime.wavelen, runtime.pilotWave, config.optVerbose, runtime.doaOnlyOpt, ...
  runtime.truth.latlonTrueDeg);

objectiveViewList = {viewRef, viewOther, viewMs};
rowList = repmat(localEmptyCaseRow(), numel(caseList), 1);
for iCase = 1:numel(caseList)
  rowList(iCase) = localBuildCaseRow(caseList(iCase), task, runtime.truth, ...
    objectiveViewList{iCase}, runtime.wavelen, runtime.pilotWave);
end
truthInitPairRows = localBuildTruthInitPairRows(caseList(3), truthInitCase, task, runtime.truth);

objectiveProbeRows = localBuildObjectiveProbeRows(caseList, objectiveViewList, task, ...
  runtime.truth, runtime.wavelen, runtime.pilotWave);
hessianProbeRows = repmat(localEmptyHessianProbeRow(), 0, 1);
noiseFreeRows = repmat(localEmptyCaseRow(), 0, 1);
alphaSweepRows = repmat(localEmptyAlphaSweepRow(), 0, 1);
if task.probeEnable
  cleanViewRef = viewRef;
  cleanViewRef.rxSigSf = selectRxSigBySat(cleanSig, runtime.refSatIdxLocal, 'singleFrame');
  cleanViewOther = viewOther;
  cleanViewOther.rxSigSf = selectRxSigBySat(cleanSig, runtime.otherSatIdxLocal, 'singleFrame');
  cleanViewMs = viewMs;
  cleanViewMs.rxSigSf = cleanSig;
  hessianProbeRows = localBuildHessianProbeRows(viewRef, viewOther, viewMs, ...
    cleanViewRef, cleanViewOther, cleanViewMs, task, runtime, pwrNoise);
  noiseFreeRows = localBuildNoiseFreeCaseRows(cleanViewRef, cleanViewOther, cleanViewMs, ...
    task, runtime);
  alphaSweepRows = localBuildAlphaSweepRows(viewMs, caseList(3), task, runtime);
end

taskOut = struct();
taskOut.caseTable = struct2table(rowList);
taskOut.truthInitPairTable = struct2table(truthInitPairRows);
taskOut.objectiveProbeTable = struct2table(objectiveProbeRows);
taskOut.hessianProbeTable = struct2table(hessianProbeRows);
taskOut.noiseFreeCaseTable = struct2table(noiseFreeRows);
taskOut.alphaSweepCaseTable = struct2table(alphaSweepRows);
end

function caseInfo = localRunDoaOnlyCase(displayName, satMode, view, wavelen, ...
  pilotWave, verbose, modelOpt, initParam)
%LOCALRUNDOAONLYCASE Run one DoA-only estimator case.

[estRaw, ~, ~] = estimatorDoaMlePilotOpt( ...
  view.sceneRef.array, wavelen, view.rxSigSf, pilotWave, ...
  view.doaGrid, initParam, verbose, modelOpt);

estResult = localWrapDoaOnlyResult(estRaw, view);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa", "none", estResult);
end

function estOut = localWrapDoaOnlyResult(estIn, view)
%LOCALWRAPDOAONLYRESULT Convert DoA-only output to the common summary format.

estOut = estIn;
estOut.modelType = 'doa-only';
estOut.fdRefEst = NaN;
estOut.fdRateEst = NaN;
if ~isfield(estOut, 'aux') || ~isstruct(estOut.aux)
  estOut.aux = struct();
end
if isfield(estIn, 'latlonEst') && ~isfield(estOut.aux, 'latlonEst')
  estOut.aux.latlonEst = estIn.latlonEst;
end
estOut.timeOffsetSec = 0;
if isfield(view, 'timeOffsetSec') && ~isempty(view.timeOffsetSec)
  estOut.timeOffsetSec = view.timeOffsetSec(1);
end
end

function row = localBuildCaseRow(caseInfo, task, truth, objectiveView, wavelen, pilotWave)
%LOCALBUILDCASE Build one repeat-level case row.

row = localEmptyCaseRow();
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.repeatIndex = task.repeatIndex;
row.displayName = string(caseInfo.displayName);
row.satMode = string(caseInfo.satMode);
row.frameMode = string(caseInfo.frameMode);
row.paramMode = string(caseInfo.paramMode);
row.dynamicMode = string(caseInfo.dynamicMode);
if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  row.tailSubtype = "missing-estimate";
  return;
end

estResult = caseInfo.estResult;
latlonEst = getDoaDopplerLatlonEst(estResult);
row.latEstDeg = latlonEst(1);
row.lonEstDeg = latlonEst(2);
if isfinite(latlonEst(1)) && isfinite(latlonEst(2))
  row.latErrDeg = latlonEst(1) - truth.latlonTrueDeg(1);
  row.lonErrDeg = latlonEst(2) - truth.latlonTrueDeg(2);
  row.angleErrDeg = calcLatlonAngleError(latlonEst, truth.latlonTrueDeg);
end
if isfield(estResult, 'initParam') && numel(estResult.initParam) == 2
  row.initAngleErrDeg = calcLatlonAngleError(estResult.initParam(:), truth.latlonTrueDeg);
end
row.exitflag = localGetFieldOrDefault(estResult, 'exitflag', NaN);
row.fval = localGetFieldOrDefault(estResult, 'fval', NaN);
if isfield(estResult, 'noiseVarEst')
  row.noiseVarEstMean = mean(estResult.noiseVarEst(:), 'omitnan');
end
if isfield(estResult, 'pathGainEst')
  row.pathGainAbsMean = mean(abs(estResult.pathGainEst(:)), 'omitnan');
end
row.rxSignalPowerMean = localMeanRxSignalPower(objectiveView.rxSigSf);
row.pilotEnergy = sum(abs(pilotWave(:)).^2);
row.numObsMean = localMeanNumObs(objectiveView.rxSigSf);
row.isResolved = localGetFieldOrDefault(estResult, 'isResolved', isfinite(row.angleErrDeg));
row.objectiveAtFinal = localEvaluateDoaOnlyObjective(objectiveView, objectiveView.rxSigSf, pilotWave, ...
  wavelen, latlonEst, true);
row.objectiveAtTruth = localEvaluateDoaOnlyObjective(objectiveView, objectiveView.rxSigSf, pilotWave, ...
  wavelen, truth.latlonTrueDeg, true);
row.objectiveTruthMinusFinal = row.objectiveAtTruth - row.objectiveAtFinal;
if ~row.isResolved
  row.tailSubtype = "unresolved";
elseif ~isfinite(row.angleErrDeg)
  row.tailSubtype = "angle-nonfinite";
else
  row.tailSubtype = "resolved";
end
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return one typed case-row template.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'repeatIndex', NaN, ...
  'displayName', "", 'satMode', "", 'frameMode', "", ...
  'paramMode', "", 'dynamicMode', "", 'latEstDeg', NaN, ...
  'lonEstDeg', NaN, 'latErrDeg', NaN, 'lonErrDeg', NaN, ...
  'angleErrDeg', NaN, 'initAngleErrDeg', NaN, 'exitflag', NaN, ...
  'fval', NaN, 'noiseVarEstMean', NaN, 'pathGainAbsMean', NaN, ...
  'rxSignalPowerMean', NaN, 'pilotEnergy', NaN, 'numObsMean', NaN, ...
  'objectiveAtFinal', NaN, ...
  'objectiveAtTruth', NaN, 'objectiveTruthMinusFinal', NaN, ...
  'isResolved', false, 'tailSubtype', "");
end

function rowList = localBuildObjectiveProbeRows(caseList, objectiveViewList, task, truth, wavelen, pilotWave)
%LOCALBUILDOBJECTIVEPROBEROWS Compare final and truth objective values.

rowTemplate = localEmptyObjectiveProbeRow();
rowList = repmat(rowTemplate, numel(caseList), 1);
for iCase = 1:numel(caseList)
  estResult = caseList(iCase).estResult;
  latlonEst = getDoaDopplerLatlonEst(estResult);
  row = rowTemplate;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.displayName = string(caseList(iCase).displayName);
  row.angleErrDeg = calcLatlonAngleError(latlonEst, truth.latlonTrueDeg);
  objectiveView = objectiveViewList{iCase};
  row.objectiveAtFinal = localEvaluateDoaOnlyObjective(objectiveView, objectiveView.rxSigSf, ...
    pilotWave, wavelen, latlonEst, true);
  row.objectiveAtTruth = localEvaluateDoaOnlyObjective(objectiveView, objectiveView.rxSigSf, ...
    pilotWave, wavelen, truth.latlonTrueDeg, true);
  row.truthMinusFinal = row.objectiveAtTruth - row.objectiveAtFinal;
  row.truthBetter = row.truthMinusFinal < 0;
  rowList(iCase) = row;
end
end

function row = localEmptyObjectiveProbeRow()
%LOCALEMPTYOBJECTIVEPROBEROW Return one typed objective-probe row.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'displayName', "", ...
  'angleErrDeg', NaN, 'objectiveAtFinal', NaN, 'objectiveAtTruth', NaN, ...
  'truthMinusFinal', NaN, 'truthBetter', false);
end

function rowList = localBuildTruthInitPairRows(baseCase, truthInitCase, task, truth)
%LOCALBUILDTRUTHINITPAIRROWS Store a compact per-task truth-init diff row.

row = localEmptyTruthInitPairRow();
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.repeatIndex = task.repeatIndex;
baseLatlon = getDoaDopplerLatlonEst(baseCase.estResult);
truthInitLatlon = getDoaDopplerLatlonEst(truthInitCase.estResult);
row.baselineAngleErrDeg = calcLatlonAngleError(baseLatlon, truth.latlonTrueDeg);
row.truthInitAngleErrDeg = calcLatlonAngleError(truthInitLatlon, truth.latlonTrueDeg);
row.estimateDiffDeg = calcLatlonAngleError(baseLatlon, truthInitLatlon);
row.baselineLatEstDeg = baseLatlon(1);
row.baselineLonEstDeg = baseLatlon(2);
row.truthInitLatEstDeg = truthInitLatlon(1);
row.truthInitLonEstDeg = truthInitLatlon(2);
row.truthInitBetter = abs(row.truthInitAngleErrDeg) < abs(row.baselineAngleErrDeg);
rowList = row;
end

function row = localEmptyTruthInitPairRow()
%LOCALEMPTYTRUTHINITPAIRROW Return one typed truth-init diff row.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'repeatIndex', NaN, ...
  'baselineAngleErrDeg', NaN, 'truthInitAngleErrDeg', NaN, ...
  'estimateDiffDeg', NaN, 'baselineLatEstDeg', NaN, ...
  'baselineLonEstDeg', NaN, 'truthInitLatEstDeg', NaN, ...
  'truthInitLonEstDeg', NaN, 'truthInitBetter', false);
end

function truthInitDiffTable = localBuildTruthInitDiffTable(truthInitPairTable, crbTable)
%LOCALBUILDTRUTHINITDIFFTABLE Aggregate the truth-init diagnostic.

rowTemplate = struct('snrDb', NaN, 'numRepeat', NaN, ...
  'baselineRmseDeg', NaN, 'truthInitRmseDeg', NaN, ...
  'baselineRmseOverDetCrb', NaN, 'truthInitRmseOverDetCrb', NaN, ...
  'truthInitBetterRate', NaN, 'medianEstimateDiffDeg', NaN, ...
  'maxEstimateDiffDeg', NaN);
if isempty(truthInitPairTable) || height(truthInitPairTable) == 0
  truthInitDiffTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
snrList = unique(truthInitPairTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(snrList), 1);
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  mask = truthInitPairTable.snrDb == snrDb;
  subTable = truthInitPairTable(mask, :);
  row = rowTemplate;
  row.snrDb = snrDb;
  row.numRepeat = height(subTable);
  row.baselineRmseDeg = localRmse(subTable.baselineAngleErrDeg);
  row.truthInitRmseDeg = localRmse(subTable.truthInitAngleErrDeg);
  crbIdx = find(crbTable.snrDb == snrDb, 1, 'first');
  if ~isempty(crbIdx)
    row.baselineRmseOverDetCrb = localSafeRatio(row.baselineRmseDeg, crbTable.angleCrbDetJointDeg(crbIdx));
    row.truthInitRmseOverDetCrb = localSafeRatio(row.truthInitRmseDeg, crbTable.angleCrbDetJointDeg(crbIdx));
  end
  row.truthInitBetterRate = mean(subTable.truthInitBetter);
  row.medianEstimateDiffDeg = median(abs(subTable.estimateDiffDeg), 'omitnan');
  row.maxEstimateDiffDeg = max(abs(subTable.estimateDiffDeg), [], 'omitnan');
  rowList(iSnr) = row;
end
truthInitDiffTable = struct2table(rowList);
end

function rowList = localBuildAlphaSweepRows(viewMs, baseCase, task, runtime)
%LOCALBUILDALPHASWEEPROWS Run probe-only multi-sat DoA weight sweep.

alphaList = reshape(runtime.config.alphaSat2List, [], 1);
rowTemplate = localEmptyAlphaSweepRow();
rowList = repmat(rowTemplate, numel(alphaList), 1);
baseLatlon = getDoaDopplerLatlonEst(baseCase.estResult);
for iAlpha = 1:numel(alphaList)
  alpha = alphaList(iAlpha);
  [~, numArray] = localWrapArray(viewMs.sceneRef.array);
  weightVec = ones(numArray, 1);
  weightVec(runtime.refSatIdxLocal) = 1;
  weightVec(runtime.otherSatIdxLocal) = alpha;
  estResult = localRunWeightedDoaOnly(viewMs, runtime.wavelen, runtime.pilotWave, ...
    baseLatlon, weightVec, runtime.config.optVerbose);
  latlonEst = estResult.latlonEst;
  row = rowTemplate;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.repeatIndex = task.repeatIndex;
  row.alphaSat2 = alpha;
  row.latEstDeg = latlonEst(1);
  row.lonEstDeg = latlonEst(2);
  row.angleErrDeg = calcLatlonAngleError(latlonEst, runtime.truth.latlonTrueDeg);
  row.objectiveAtFinal = estResult.fval;
  row.exitflag = estResult.exitflag;
  rowList(iAlpha) = row;
end
end

function row = localEmptyAlphaSweepRow()
%LOCALEMPTYALPHASWEEPROW Return one typed alpha-sweep row.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'repeatIndex', NaN, ...
  'alphaSat2', NaN, 'latEstDeg', NaN, 'lonEstDeg', NaN, ...
  'angleErrDeg', NaN, 'objectiveAtFinal', NaN, 'exitflag', NaN);
end

function estResult = localRunWeightedDoaOnly(view, wavelen, pilotWave, initParam, weightVec, verbose)
%LOCALRUNWEIGHTEDDOAONLY Run a replay-local weighted DoA objective.

if iscell(view.doaGrid)
  baseGrid = view.doaGrid{1};
else
  baseGrid = view.doaGrid;
end
lb = baseGrid.range(:, 1);
ub = baseGrid.range(:, 2);
initParam = reshape(initParam, [], 1);
if numel(initParam) ~= 2 || any(~isfinite(initParam))
  initParam = mean(baseGrid.range, 2);
end
initParam = min(max(initParam, lb), ub);
if verbose
  displayMode = 'iter';
else
  displayMode = 'off';
end
optimOpt = optimoptions('fmincon', 'Display', displayMode, 'Algorithm', 'interior-point');
objFun = @(x) localEvaluateDoaOnlyObjective(view, view.rxSigSf, pilotWave, ...
  wavelen, x, true, weightVec);
[optVar, fval, exitflag] = fmincon(objFun, initParam, [], [], [], [], lb, ub, [], optimOpt);
estResult = struct();
estResult.latlonEst = optVar(:);
estResult.doaParamEst = optVar(:);
estResult.initParam = initParam;
estResult.fval = fval;
estResult.exitflag = exitflag;
end

function rowList = localBuildNoiseFreeCaseRows(cleanViewRef, cleanViewOther, cleanViewMs, task, runtime)
%LOCALBUILDNOISEFREECASEROWS Run a compact noise-free DoA sanity probe.

caseList = [ ...
  localRunDoaOnlyCase("SS-SF-DoA-Clean", "single", cleanViewRef, runtime.wavelen, ...
    runtime.pilotWave, runtime.config.optVerbose, runtime.doaOnlyOpt, []), ...
  localRunDoaOnlyCase("Other-SF-DoA-Clean", "other", cleanViewOther, runtime.wavelen, ...
    runtime.pilotWave, runtime.config.optVerbose, runtime.doaOnlyOpt, []), ...
  localRunDoaOnlyCase("MS-SF-DoA-Clean", "multi", cleanViewMs, runtime.wavelen, ...
    runtime.pilotWave, runtime.config.optVerbose, runtime.doaOnlyOpt, [])];
objectiveViewList = {cleanViewRef, cleanViewOther, cleanViewMs};
rowList = repmat(localEmptyCaseRow(), numel(caseList), 1);
for iCase = 1:numel(caseList)
  rowList(iCase) = localBuildCaseRow(caseList(iCase), task, runtime.truth, ...
    objectiveViewList{iCase}, runtime.wavelen, runtime.pilotWave);
end
end

function rowList = localBuildHessianProbeRows(viewRefNoisy, viewOtherNoisy, viewMsNoisy, ...
  viewRefClean, viewOtherClean, viewMsClean, task, runtime, pwrNoise)
%LOCALBUILDHESSIANPROBEROWS Build per-sat numerical Hessian rows.

rowTemplate = localEmptyHessianProbeRow();
rowList = repmat(rowTemplate, 0, 1);
sourceTagList = ["noisy"; "clean"];
viewSet = { ...
  {viewRefNoisy, viewOtherNoisy, viewMsNoisy}; ...
  {viewRefClean, viewOtherClean, viewMsClean}};
scopeNameList = ["ref"; "other"; "joint"];
satModeList = ["single"; "other"; "multi"];
sceneList = {runtime.sceneRefOnly; runtime.sceneOtherOnly; runtime.scene};
satIdxGlobalList = [runtime.truth.refSatIdxGlobal; runtime.truth.otherSatIdxGlobal; NaN];
crbDetOpt = struct('doaType', 'latlon');
for iSource = 1:numel(sourceTagList)
  currentViewList = viewSet{iSource};
  for iScope = 1:numel(scopeNameList)
    currentView = currentViewList{iScope};
    row = rowTemplate;
    row.snrDb = task.snrDb;
    row.taskSeed = task.taskSeed;
    row.sourceTag = sourceTagList(iSource);
    row.scopeName = scopeNameList(iScope);
    row.satMode = satModeList(iScope);
    row.satIdxGlobal = satIdxGlobalList(iScope);
    H = localNumericalDoaHessian(currentView, currentView.rxSigSf, runtime.pilotWave, ...
      runtime.wavelen, runtime.truth.latlonTrueDeg, runtime.config.hessianStepDeg);
    H = 0.5 * (H + H.');
    row.hessianCond = cond(H);
    row.hessianMinEig = min(eig(H));
    row.hessianMaxEig = max(eig(H));
    row.hessianTrace = trace(H);
    row.hessianDet = det(H);
    row.objectiveAtTruth = localEvaluateDoaOnlyObjective(currentView, currentView.rxSigSf, ...
      runtime.pilotWave, runtime.wavelen, runtime.truth.latlonTrueDeg, true);
    row.hessianRepeat = 1;

    [~, detAux] = localRunOneDetDoaCrb(sceneList{iScope}, runtime, pwrNoise, crbDetOpt);
    fimDet = detAux.fim;
    row.detFimTrace = trace(fimDet);
    row.detFimCond = cond(fimDet);
    row.detFimAdditivityRelErr = localFimAdditivityRelErr(detAux);
    row.hessianToDetFimTraceRatio = trace(H) / trace(fimDet);
    row.hessianToDetFimShapeErr = localShapeError(H, fimDet);
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
if runtime.config.hessianMeanRepeat > 1
  meanRows = localBuildMeanNoisyHessianProbeRows(task, runtime, pwrNoise, crbDetOpt, rowTemplate);
  rowList = [rowList; meanRows];
end
end

function rowList = localBuildMeanNoisyHessianProbeRows(task, runtime, pwrNoise, crbDetOpt, rowTemplate)
%LOCALBUILDMEANNOISYHESSIANPROBEROWS Average noisy Hessian over extra noise draws.

scopeNameList = ["ref"; "other"; "joint"];
satModeList = ["single"; "other"; "multi"];
sceneList = {runtime.sceneRefOnly; runtime.sceneOtherOnly; runtime.scene};
satIdxGlobalList = [runtime.truth.refSatIdxGlobal; runtime.truth.otherSatIdxGlobal; NaN];
HsumCell = repmat({zeros(2, 2)}, numel(scopeNameList), 1);
objSum = zeros(numel(scopeNameList), 1);
numMean = max(1, round(runtime.config.hessianMeanRepeat));
for iMean = 1:numMean
  rng(task.taskSeed + 1000003 + 7919 * iMean + 101 * task.snrIndex, 'twister');
  [rxSigMean, ~, ~, ~, ~] = genMultiSatSnapshots( ...
    runtime.steeringInfo, runtime.pilotWave, runtime.linkParam, ...
    runtime.carrierFreq, runtime.sampleRate, pwrNoise, runtime.pathGain, ...
    runtime.snapOpt);
  viewRef = runtime.viewRefTemplate;
  viewRef.rxSigSf = selectRxSigBySat(rxSigMean, runtime.refSatIdxLocal, 'singleFrame');
  viewOther = runtime.viewOtherTemplate;
  viewOther.rxSigSf = selectRxSigBySat(rxSigMean, runtime.otherSatIdxLocal, 'singleFrame');
  viewMs = runtime.viewMsTemplate;
  viewMs.rxSigSf = rxSigMean;
  viewList = {viewRef; viewOther; viewMs};
  for iScope = 1:numel(scopeNameList)
    currentView = viewList{iScope};
    H = localNumericalDoaHessian(currentView, currentView.rxSigSf, runtime.pilotWave, ...
      runtime.wavelen, runtime.truth.latlonTrueDeg, runtime.config.hessianStepDeg);
    HsumCell{iScope} = HsumCell{iScope} + 0.5 * (H + H.');
    objSum(iScope) = objSum(iScope) + localEvaluateDoaOnlyObjective(currentView, ...
      currentView.rxSigSf, runtime.pilotWave, runtime.wavelen, ...
      runtime.truth.latlonTrueDeg, true);
  end
end
rowList = repmat(rowTemplate, numel(scopeNameList), 1);
for iScope = 1:numel(scopeNameList)
  H = HsumCell{iScope} / numMean;
  row = rowTemplate;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.sourceTag = "mean-noisy";
  row.scopeName = scopeNameList(iScope);
  row.satMode = satModeList(iScope);
  row.satIdxGlobal = satIdxGlobalList(iScope);
  row.objectiveAtTruth = objSum(iScope) / numMean;
  row.hessianRepeat = numMean;
  row.hessianCond = cond(H);
  row.hessianMinEig = min(eig(H));
  row.hessianMaxEig = max(eig(H));
  row.hessianTrace = trace(H);
  row.hessianDet = det(H);
  [~, detAux] = localRunOneDetDoaCrb(sceneList{iScope}, runtime, pwrNoise, crbDetOpt);
  fimDet = detAux.fim;
  row.detFimTrace = trace(fimDet);
  row.detFimCond = cond(fimDet);
  row.detFimAdditivityRelErr = localFimAdditivityRelErr(detAux);
  row.hessianToDetFimTraceRatio = trace(H) / trace(fimDet);
  row.hessianToDetFimShapeErr = localShapeError(H, fimDet);
  rowList(iScope) = row;
end
end

function row = localEmptyHessianProbeRow()
%LOCALEMPTYHESSIANPROBEROW Return one typed Hessian-probe row.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'sourceTag', "", ...
  'scopeName', "", 'satMode', "", 'satIdxGlobal', NaN, ...
  'objectiveAtTruth', NaN, 'hessianRepeat', NaN, 'hessianTrace', NaN, 'hessianDet', NaN, ...
  'hessianCond', NaN, 'hessianMinEig', NaN, 'hessianMaxEig', NaN, ...
  'detFimTrace', NaN, 'detFimCond', NaN, 'detFimAdditivityRelErr', NaN, ...
  'hessianToDetFimTraceRatio', NaN, 'hessianToDetFimShapeErr', NaN);
end

function crbTable = localBuildCrbTable(runtime, config)
%LOCALBUILDCRBTABLE Evaluate static DoA-Doppler and DoA-only CRB over SNR.

snrDbList = reshape(config.snrDbList, [], 1);
rowTemplate = struct('snrDb', NaN, 'angleCrbDdSingleDeg', NaN, ...
  'angleCrbDdOtherDeg', NaN, 'angleCrbDdJointDeg', NaN, ...
  'angleCrbDetSingleDeg', NaN, 'angleCrbDetOtherDeg', NaN, ...
  'angleCrbDetJointDeg', NaN, 'angleCrbAmpSingleDeg', NaN, ...
  'angleCrbAmpOtherDeg', NaN, 'angleCrbAmpJointDeg', NaN, ...
  'ampToDetSingleRatio', NaN, 'ampToDetOtherRatio', NaN, ...
  'ampToDetJointRatio', NaN, 'ddToDetSingleRatio', NaN, ...
  'ddToDetOtherRatio', NaN, 'ddToDetJointRatio', NaN, ...
  'jointDetGainOverSingle', NaN, 'jointDetGainOverOther', NaN, ...
  'jointDdGainOverSingle', NaN, 'jointDdGainOverOther', NaN, ...
  'detFimAdditivityRelErr', NaN, 'ampFimAdditivityRelErr', NaN, ...
  'detCrbSingleMat', {[]}, 'detCrbOtherMat', {[]}, 'detCrbJointMat', {[]}, ...
  'ampCrbSingleMat', {[]}, 'ampCrbOtherMat', {[]}, 'ampCrbJointMat', {[]}, ...
  'detFimTraceSingle', NaN, 'detFimTraceOther', NaN, ...
  'detFimTraceJoint', NaN, 'ampFimTraceSingle', NaN, ...
  'ampFimTraceOther', NaN, 'ampFimTraceJoint', NaN);
rowList = repmat(rowTemplate, numel(snrDbList), 1);
crbOpt = struct('doaType', 'latlon');
for iSnr = 1:numel(snrDbList)
  pwrNoise = runtime.pwrSource / (10^(snrDbList(iSnr) / 10));
  [crbDdSingle, ~] = crbPilotSfDoaDoppler(runtime.sceneRefOnly, ...
    runtime.pilotWave, runtime.carrierFreq, runtime.sampleRate, ...
    runtime.truth.latlonTrueDeg, runtime.truth.fdRefTrueHz, 1, pwrNoise, crbOpt);
  [crbDdOther, ~] = crbPilotSfDoaDoppler(runtime.sceneOtherOnly, ...
    runtime.pilotWave, runtime.carrierFreq, runtime.sampleRate, ...
    runtime.truth.latlonTrueDeg, runtime.truth.fdOtherTrueHz, 1, pwrNoise, crbOpt);
  [crbDdJoint, ~] = crbPilotSfDoaDoppler(runtime.scene, ...
    runtime.pilotWave, runtime.carrierFreq, runtime.sampleRate, ...
    runtime.truth.latlonTrueDeg, runtime.truth.fdRefTrueHz, 1, pwrNoise, crbOpt);
  [crbDetSingle, detAuxSingle] = localRunOneDetDoaCrb(runtime.sceneRefOnly, runtime, pwrNoise, crbOpt);
  [crbDetOther, detAuxOther] = localRunOneDetDoaCrb(runtime.sceneOtherOnly, runtime, pwrNoise, crbOpt);
  [crbDetJoint, detAuxJoint] = localRunOneDetDoaCrb(runtime.scene, runtime, pwrNoise, crbOpt);
  [crbAmpSingle, ampAuxSingle] = localRunOneDetDoaCrb(runtime.sceneRefOnly, runtime, pwrNoise, crbOpt, ...
    runtime.truth.effectivePilotGainAbsRef);
  [crbAmpOther, ampAuxOther] = localRunOneDetDoaCrb(runtime.sceneOtherOnly, runtime, pwrNoise, crbOpt, ...
    runtime.truth.effectivePilotGainAbsOther);
  [crbAmpJoint, ampAuxJoint] = localRunOneDetDoaCrb(runtime.scene, runtime, pwrNoise, crbOpt, ...
    runtime.truth.effectivePilotGainAbs);

  row = rowTemplate;
  row.snrDb = snrDbList(iSnr);
  row.angleCrbDdSingleDeg = projectCrbToAngleMetric(crbDdSingle(1:2, 1:2), ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbDdOtherDeg = projectCrbToAngleMetric(crbDdOther(1:2, 1:2), ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbDdJointDeg = projectCrbToAngleMetric(crbDdJoint(1:2, 1:2), ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbDetSingleDeg = projectCrbToAngleMetric(crbDetSingle, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbDetOtherDeg = projectCrbToAngleMetric(crbDetOther, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbDetJointDeg = projectCrbToAngleMetric(crbDetJoint, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbAmpSingleDeg = projectCrbToAngleMetric(crbAmpSingle, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbAmpOtherDeg = projectCrbToAngleMetric(crbAmpOther, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.angleCrbAmpJointDeg = projectCrbToAngleMetric(crbAmpJoint, ...
    runtime.truth.latlonTrueDeg, 'latlon');
  row.ampToDetSingleRatio = localSafeRatio(row.angleCrbAmpSingleDeg, row.angleCrbDetSingleDeg);
  row.ampToDetOtherRatio = localSafeRatio(row.angleCrbAmpOtherDeg, row.angleCrbDetOtherDeg);
  row.ampToDetJointRatio = localSafeRatio(row.angleCrbAmpJointDeg, row.angleCrbDetJointDeg);
  row.ddToDetSingleRatio = localSafeRatio(row.angleCrbDdSingleDeg, row.angleCrbDetSingleDeg);
  row.ddToDetOtherRatio = localSafeRatio(row.angleCrbDdOtherDeg, row.angleCrbDetOtherDeg);
  row.ddToDetJointRatio = localSafeRatio(row.angleCrbDdJointDeg, row.angleCrbDetJointDeg);
  row.jointDetGainOverSingle = localSafeRatio(row.angleCrbDetSingleDeg, row.angleCrbDetJointDeg);
  row.jointDetGainOverOther = localSafeRatio(row.angleCrbDetOtherDeg, row.angleCrbDetJointDeg);
  row.jointDdGainOverSingle = localSafeRatio(row.angleCrbDdSingleDeg, row.angleCrbDdJointDeg);
  row.jointDdGainOverOther = localSafeRatio(row.angleCrbDdOtherDeg, row.angleCrbDdJointDeg);
  row.detFimAdditivityRelErr = localFimAdditivityRelErr(detAuxJoint);
  row.ampFimAdditivityRelErr = localFimAdditivityRelErr(ampAuxJoint);
  row.detCrbSingleMat = {crbDetSingle};
  row.detCrbOtherMat = {crbDetOther};
  row.detCrbJointMat = {crbDetJoint};
  row.ampCrbSingleMat = {crbAmpSingle};
  row.ampCrbOtherMat = {crbAmpOther};
  row.ampCrbJointMat = {crbAmpJoint};
  row.detFimTraceSingle = trace(detAuxSingle.fim);
  row.detFimTraceOther = trace(detAuxOther.fim);
  row.detFimTraceJoint = trace(detAuxJoint.fim);
  row.ampFimTraceSingle = trace(ampAuxSingle.fim);
  row.ampFimTraceOther = trace(ampAuxOther.fim);
  row.ampFimTraceJoint = trace(ampAuxJoint.fim);
  rowList(iSnr) = row;
end
crbTable = struct2table(rowList);
end

function [crb, aux] = localRunOneDetDoaCrb(scene, runtime, pwrNoise, crbOpt, gainAbs)
%LOCALRUNONEDETDOACRB Run deterministic DoA-only CRB for one scene.

[arrayCell, numArray] = localWrapArray(scene.array);
if nargin < 5 || isempty(gainAbs)
  gainAbs = ones(numArray, 1);
else
  gainAbs = reshape(double(abs(gainAbs)), [], 1);
  if isscalar(gainAbs)
    gainAbs = repmat(gainAbs, numArray, 1);
  elseif numel(gainAbs) ~= numArray
    error('replaySfMsDoaCrbDiagnose:GainCountMismatch', ...
      'gainAbs must be scalar or contain one gain per array.');
  end
end
localDoaCell = cell(1, numArray);
localDoaJacCell = cell(1, numArray);
for iArray = 1:numArray
  doaGrid = genDoaGrid("latlon", 2, runtime.config.gridSize, ...
    [runtime.truth.latlonTrueDeg(1) - runtime.config.searchMarginDeg, ...
     runtime.truth.latlonTrueDeg(1) + runtime.config.searchMarginDeg; ...
     runtime.truth.latlonTrueDeg(2) - runtime.config.searchMarginDeg, ...
     runtime.truth.latlonTrueDeg(2) + runtime.config.searchMarginDeg], ...
    'eci', datevec(scene.utc), scene.satPosEci(:, iArray), scene.rotMat{iArray}, ...
    referenceEllipsoid('sphere'));
  localDoaCell{iArray} = localLatlonToAngle(doaGrid, runtime.truth.latlonTrueDeg(1), ...
    runtime.truth.latlonTrueDeg(2));
  localDoaJacCell{iArray} = localFiniteDiffLocalDoaJac(doaGrid, runtime.truth.latlonTrueDeg);
end
pilotMeanPower = mean(abs(runtime.pilotWave(:)).^2);
numSnap = numel(runtime.pilotWave);
pwrSourceCell = arrayfun(@(g) pilotMeanPower * (abs(g).^2), gainAbs(:).', ...
  'UniformOutput', false);
crbDetOpt = struct();
crbDetOpt.localDoaJac = localDoaJacCell;
crbDetOpt.paramName = {'latDeg', 'lonDeg'};
[crb, aux] = crbDetDoa(arrayCell, runtime.wavelen, localDoaCell, ...
  pwrSourceCell, pwrNoise, numSnap, crbDetOpt);
aux.requestedDopplerCrbOpt = crbOpt;
aux.effectivePilotGainAbs = gainAbs;
end

function perfTable = localAttachCrbMetrics(caseTable, crbTable)
%LOCALATTACHCRBMETRICS Attach CRB std and normalized errors to repeat rows.

perfTable = caseTable;
if isempty(perfTable) || height(perfTable) == 0
  return;
end
perfTable.angleCrbDdStdDeg = nan(height(perfTable), 1);
perfTable.angleCrbDetStdDeg = nan(height(perfTable), 1);
for iRow = 1:height(perfTable)
  crbIdx = find(crbTable.snrDb == perfTable.snrDb(iRow), 1, 'first');
  if isempty(crbIdx)
    continue;
  end
  if perfTable.satMode(iRow) == "multi"
    perfTable.angleCrbDdStdDeg(iRow) = crbTable.angleCrbDdJointDeg(crbIdx);
    perfTable.angleCrbDetStdDeg(iRow) = crbTable.angleCrbDetJointDeg(crbIdx);
  elseif perfTable.satMode(iRow) == "other"
    perfTable.angleCrbDdStdDeg(iRow) = crbTable.angleCrbDdOtherDeg(crbIdx);
    perfTable.angleCrbDetStdDeg(iRow) = crbTable.angleCrbDetOtherDeg(crbIdx);
  else
    perfTable.angleCrbDdStdDeg(iRow) = crbTable.angleCrbDdSingleDeg(crbIdx);
    perfTable.angleCrbDetStdDeg(iRow) = crbTable.angleCrbDetSingleDeg(crbIdx);
  end
end
perfTable.angleNormErrDd = abs(perfTable.angleErrDeg) ./ perfTable.angleCrbDdStdDeg;
perfTable.angleNormErrDet = abs(perfTable.angleErrDeg) ./ perfTable.angleCrbDetStdDeg;
perfTable.crbLocalMaskDd = perfTable.isResolved & isfinite(perfTable.angleNormErrDd);
perfTable.crbLocalMaskDet = perfTable.isResolved & isfinite(perfTable.angleNormErrDet);
for iRow = 1:height(perfTable)
  if perfTable.tailSubtype(iRow) ~= "resolved"
    continue;
  end
  ddTail = isfinite(perfTable.angleNormErrDd(iRow)) && perfTable.angleNormErrDd(iRow) > 5;
  detTail = isfinite(perfTable.angleNormErrDet(iRow)) && perfTable.angleNormErrDet(iRow) > 5;
  if ddTail && detTail
    perfTable.tailSubtype(iRow) = "angle-dd-det-crb-tail";
  elseif ddTail
    perfTable.tailSubtype(iRow) = "angle-dd-crb-tail";
  elseif detTail
    perfTable.tailSubtype(iRow) = "angle-det-crb-tail";
  else
    perfTable.tailSubtype(iRow) = "crb-local";
  end
end
end

function aggregateTable = localBuildAggregateTable(perfTable, config)
%LOCALBUILDAGGREGATETABLE Build aggregate metrics for DD and DoA-only CRB.

rowTemplate = localEmptyAggregateRow();
if isempty(perfTable) || height(perfTable) == 0
  aggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'satMode', 'frameMode', ...
  'paramMode', 'dynamicMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iGroup = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iGroup) & ...
    perfTable.snrDb == groupKey.snrDb(iGroup);
  row = localBuildOneAggregateRow(perfTable, mask, config, rowTemplate);
  row.displayName = groupKey.displayName(iGroup);
  row.satMode = groupKey.satMode(iGroup);
  row.frameMode = groupKey.frameMode(iGroup);
  row.paramMode = groupKey.paramMode(iGroup);
  row.dynamicMode = groupKey.dynamicMode(iGroup);
  row.snrDb = groupKey.snrDb(iGroup);
  rowList(iGroup) = row;
end
aggregateTable = struct2table(rowList);
end

function row = localBuildOneAggregateRow(perfTable, mask, config, row)
%LOCALBUILDONEAGGREGATEROW Fill one aggregate row.

angleErr = perfTable.angleErrDeg(mask);
angleNormDd = perfTable.angleNormErrDd(mask);
angleNormDet = perfTable.angleNormErrDet(mask);
isResolved = perfTable.isResolved(mask);
row.numRepeat = nnz(mask);
row.angleCrbDdStdDeg = median(perfTable.angleCrbDdStdDeg(mask), 'omitnan');
row.angleCrbDetStdDeg = median(perfTable.angleCrbDetStdDeg(mask), 'omitnan');
row.fullFiniteRate = mean(isfinite(angleErr));
row.resolvedRate = mean(isResolved);
row.fullAngleRmseDeg = localRmse(angleErr(isfinite(angleErr)));
row.fullAngleP95Deg = localPercentile(abs(angleErr(isfinite(angleErr))), 95);
row.fullAngleRmseOverDdCrb = localSafeRatio(row.fullAngleRmseDeg, row.angleCrbDdStdDeg);
row.fullAngleRmseOverDetCrb = localSafeRatio(row.fullAngleRmseDeg, row.angleCrbDetStdDeg);
row.resolvedAngleRmseDeg = localRmse(angleErr(isResolved & isfinite(angleErr)));
row.resolvedAngleP95Deg = localPercentile(abs(angleErr(isResolved & isfinite(angleErr))), 95);
row.resolvedAngleRmseOverDdCrb = localSafeRatio(row.resolvedAngleRmseDeg, row.angleCrbDdStdDeg);
row.resolvedAngleRmseOverDetCrb = localSafeRatio(row.resolvedAngleRmseDeg, row.angleCrbDetStdDeg);

crbLocalMaskDd = isResolved & isfinite(angleNormDd);
crbLocalMaskDet = isResolved & isfinite(angleNormDet);
if config.trimEnable
  crbLocalMaskDd = crbLocalMaskDd & angleNormDd <= config.trimNormCap;
  crbLocalMaskDet = crbLocalMaskDet & angleNormDet <= config.trimNormCap;
end
row.crbLocalKeepRateDd = mean(crbLocalMaskDd);
row.crbLocalKeepRateDet = mean(crbLocalMaskDet);
row.crbLocalAngleRmseDegDd = localRmse(angleErr(crbLocalMaskDd & isfinite(angleErr)));
row.crbLocalAngleRmseDegDet = localRmse(angleErr(crbLocalMaskDet & isfinite(angleErr)));
row.crbLocalAngleRmseOverDdCrb = localSafeRatio(row.crbLocalAngleRmseDegDd, row.angleCrbDdStdDeg);
row.crbLocalAngleRmseOverDetCrb = localSafeRatio(row.crbLocalAngleRmseDegDet, row.angleCrbDetStdDeg);
row.outlierRateDd = 1 - row.crbLocalKeepRateDd;
row.outlierRateDet = 1 - row.crbLocalKeepRateDet;
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return one typed aggregate-row template.

row = struct('displayName', "", 'satMode', "", 'frameMode', "", ...
  'paramMode', "", 'dynamicMode', "", 'snrDb', NaN, 'numRepeat', NaN, ...
  'angleCrbDdStdDeg', NaN, 'angleCrbDetStdDeg', NaN, ...
  'fullFiniteRate', NaN, 'resolvedRate', NaN, 'fullAngleRmseDeg', NaN, ...
  'fullAngleP95Deg', NaN, 'fullAngleRmseOverDdCrb', NaN, ...
  'fullAngleRmseOverDetCrb', NaN, 'resolvedAngleRmseDeg', NaN, ...
  'resolvedAngleP95Deg', NaN, 'resolvedAngleRmseOverDdCrb', NaN, ...
  'resolvedAngleRmseOverDetCrb', NaN, 'crbLocalKeepRateDd', NaN, ...
  'crbLocalKeepRateDet', NaN, 'crbLocalAngleRmseDegDd', NaN, ...
  'crbLocalAngleRmseDegDet', NaN, 'crbLocalAngleRmseOverDdCrb', NaN, ...
  'crbLocalAngleRmseOverDetCrb', NaN, 'outlierRateDd', NaN, ...
  'outlierRateDet', NaN);
end

function topTailTable = localBuildTopTailTable(perfTable, topCountPerGroup)
%LOCALBUILDTOPTAILTABLE Select top normalized-error tail rows per group.

rowTemplate = localEmptyTopTailRow();
if isempty(perfTable) || height(perfTable) == 0 || topCountPerGroup <= 0
  topTailTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(perfTable(:, {'displayName', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, 0, 1);
for iGroup = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iGroup) & perfTable.snrDb == groupKey.snrDb(iGroup);
  subTable = perfTable(mask, :);
  tailScore = max([subTable.angleNormErrDd, subTable.angleNormErrDet], [], 2, 'omitnan');
  tailScore(~isfinite(tailScore)) = -inf;
  [~, orderIdx] = sort(tailScore, 'descend');
  orderIdx = orderIdx(1:min(topCountPerGroup, numel(orderIdx)));
  for iRank = 1:numel(orderIdx)
    src = subTable(orderIdx(iRank), :);
    row = rowTemplate;
    row.displayName = src.displayName;
    row.snrDb = src.snrDb;
    row.rankInGroup = iRank;
    row.taskSeed = src.taskSeed;
    row.repeatIndex = src.repeatIndex;
    row.angleErrDeg = src.angleErrDeg;
    row.angleNormErrDd = src.angleNormErrDd;
    row.angleNormErrDet = src.angleNormErrDet;
    row.tailSubtype = src.tailSubtype;
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
topTailTable = struct2table(rowList);
end

function row = localEmptyTopTailRow()
%LOCALEMPTYTOPTAILROW Return one typed top-tail row.

row = struct('displayName', "", 'snrDb', NaN, 'rankInGroup', NaN, ...
  'taskSeed', NaN, 'repeatIndex', NaN, 'angleErrDeg', NaN, ...
  'angleNormErrDd', NaN, 'angleNormErrDet', NaN, 'tailSubtype', "");
end

function crbCompareTable = localBuildCrbCompareTable(crbTable)
%LOCALBUILDCRBCOMPARETABLE Build compact CRB comparison table.

if isempty(crbTable) || height(crbTable) == 0
  crbCompareTable = table();
  return;
end
keepVar = {'snrDb', 'angleCrbDdSingleDeg', 'angleCrbDetSingleDeg', ...
  'ddToDetSingleRatio', 'angleCrbAmpSingleDeg', 'ampToDetSingleRatio', ...
  'angleCrbDdOtherDeg', 'angleCrbDetOtherDeg', 'ddToDetOtherRatio', ...
  'angleCrbAmpOtherDeg', 'ampToDetOtherRatio', 'angleCrbDdJointDeg', ...
  'angleCrbDetJointDeg', 'ddToDetJointRatio', 'angleCrbAmpJointDeg', ...
  'ampToDetJointRatio', 'jointDdGainOverSingle', 'jointDetGainOverSingle', ...
  'jointDdGainOverOther', 'jointDetGainOverOther', 'detFimAdditivityRelErr', ...
  'ampFimAdditivityRelErr'};
crbCompareTable = crbTable(:, keepVar);
end

function doaGainTable = localBuildDoaGainTable(aggregateTable, crbTable)
%LOCALBUILDDOAGAINTABLE Compare expected CRB gain with realized MLE gain.

rowTemplate = struct('snrDb', NaN, 'ssRmseDeg', NaN, 'otherRmseDeg', NaN, ...
  'msRmseDeg', NaN, 'expectedJointGainDet', NaN, ...
  'expectedJointGainDd', NaN, 'realizedMleGainVsRef', NaN, ...
  'realizedMleGainVsBestSingle', NaN, 'gainRealizationRateDet', NaN, ...
  'gainRealizationRateBestDet', NaN);
if isempty(crbTable) || height(crbTable) == 0
  doaGainTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
rowList = repmat(rowTemplate, height(crbTable), 1);
for iSnr = 1:height(crbTable)
  snrDb = crbTable.snrDb(iSnr);
  ssRmse = localFindAggregateMetric(aggregateTable, "SS-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  otherRmse = localFindAggregateMetric(aggregateTable, "Other-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  msRmse = localFindAggregateMetric(aggregateTable, "MS-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  row = rowTemplate;
  row.snrDb = snrDb;
  row.ssRmseDeg = ssRmse;
  row.otherRmseDeg = otherRmse;
  row.msRmseDeg = msRmse;
  row.expectedJointGainDet = crbTable.jointDetGainOverSingle(iSnr);
  row.expectedJointGainDd = crbTable.jointDdGainOverSingle(iSnr);
  row.realizedMleGainVsRef = localSafeRatio(ssRmse, msRmse);
  row.realizedMleGainVsBestSingle = localSafeRatio(min([ssRmse, otherRmse], [], 'omitnan'), msRmse);
  row.gainRealizationRateDet = localSafeRatio(row.realizedMleGainVsRef, row.expectedJointGainDet);
  expectedBestGain = localSafeRatio(min([crbTable.angleCrbDetSingleDeg(iSnr), ...
    crbTable.angleCrbDetOtherDeg(iSnr)], [], 'omitnan'), crbTable.angleCrbDetJointDeg(iSnr));
  row.gainRealizationRateBestDet = localSafeRatio(row.realizedMleGainVsBestSingle, expectedBestGain);
  rowList(iSnr) = row;
end
doaGainTable = struct2table(rowList);
end

function value = localFindAggregateMetric(aggregateTable, displayName, snrDb, metricName)
%LOCALFINDAGGREGATEMETRIC Read one aggregate metric by method and SNR.

value = NaN;
if isempty(aggregateTable) || height(aggregateTable) == 0 || ~ismember(metricName, aggregateTable.Properties.VariableNames)
  return;
end
mask = aggregateTable.displayName == displayName & aggregateTable.snrDb == snrDb;
if any(mask)
  metricVec = aggregateTable.(metricName);
  value = metricVec(find(mask, 1, 'first'));
end
end

function otherAuditTable = localBuildOtherDoaCrbAuditTable(perfTable, crbTable, truth, config)
%LOCALBUILDOTHERDOACRBAUDITTABLE Audit other-sat DoA-only errors vs CRB shape.

rowTemplate = struct('snrDb', NaN, 'numRepeat', NaN, ...
  'otherRmseDeg', NaN, 'otherCrbDetDeg', NaN, 'otherRmseOverDetCrb', NaN, ...
  'crbLocalKeepRateDet', NaN, 'crbLocalRmseOverDetCrb', NaN, ...
  'biasLatDeg', NaN, 'biasLonDeg', NaN, 'errCovEig1', NaN, ...
  'errCovEig2', NaN, 'detCrbEig1', NaN, 'detCrbEig2', NaN, ...
  'ampCrbEig1', NaN, 'ampCrbEig2', NaN, ...
  'eigRatio1', NaN, 'eigRatio2', NaN, ...
  'ampEigRatio1', NaN, 'ampEigRatio2', NaN);
if isempty(perfTable) || height(perfTable) == 0
  otherAuditTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
snrList = unique(perfTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(snrList), 1);
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  mask = perfTable.displayName == "Other-SF-DoA" & perfTable.snrDb == snrDb;
  subTable = perfTable(mask, :);
  row = rowTemplate;
  row.snrDb = snrDb;
  row.numRepeat = height(subTable);
  if height(subTable) > 0
    latErr = subTable.latEstDeg - truth.latlonTrueDeg(1);
    lonErr = subTable.lonEstDeg - truth.latlonTrueDeg(2);
    errMat = [latErr(:), lonErr(:)];
    finiteMask = all(isfinite(errMat), 2);
    row.biasLatDeg = mean(latErr(finiteMask), 'omitnan');
    row.biasLonDeg = mean(lonErr(finiteMask), 'omitnan');
    row.otherRmseDeg = localRmse(subTable.angleErrDeg);
    row.otherCrbDetDeg = median(subTable.angleCrbDetStdDeg, 'omitnan');
    row.otherRmseOverDetCrb = localSafeRatio(row.otherRmseDeg, row.otherCrbDetDeg);
    crbMask = subTable.isResolved & isfinite(subTable.angleNormErrDet);
    if config.trimEnable
      crbMask = crbMask & subTable.angleNormErrDet <= config.trimNormCap;
    end
    row.crbLocalKeepRateDet = mean(crbMask);
    row.crbLocalRmseOverDetCrb = localSafeRatio(localRmse(subTable.angleErrDeg(crbMask)), row.otherCrbDetDeg);
    if nnz(finiteMask) >= 2
      errCov = (errMat(finiteMask, :)' * errMat(finiteMask, :)) / nnz(finiteMask);
      errEig = sort(real(eig(0.5 * (errCov + errCov.'))), 'ascend');
      row.errCovEig1 = errEig(1);
      row.errCovEig2 = errEig(2);
    end
  end
  crbIdx = find(crbTable.snrDb == snrDb, 1, 'first');
  if ~isempty(crbIdx) && ismember('detCrbOtherMat', crbTable.Properties.VariableNames)
    crbMat = localTableCellValue(crbTable.detCrbOtherMat, crbIdx);
    crbEig = sort(real(eig(0.5 * (crbMat + crbMat.'))), 'ascend');
    row.detCrbEig1 = crbEig(1);
    row.detCrbEig2 = crbEig(2);
    row.eigRatio1 = localSafeRatio(row.errCovEig1, row.detCrbEig1);
    row.eigRatio2 = localSafeRatio(row.errCovEig2, row.detCrbEig2);
  end
  if ~isempty(crbIdx) && ismember('ampCrbOtherMat', crbTable.Properties.VariableNames)
    ampCrbMat = localTableCellValue(crbTable.ampCrbOtherMat, crbIdx);
    ampCrbEig = sort(real(eig(0.5 * (ampCrbMat + ampCrbMat.'))), 'ascend');
    row.ampCrbEig1 = ampCrbEig(1);
    row.ampCrbEig2 = ampCrbEig(2);
    row.ampEigRatio1 = localSafeRatio(row.errCovEig1, row.ampCrbEig1);
    row.ampEigRatio2 = localSafeRatio(row.errCovEig2, row.ampCrbEig2);
  end
  rowList(iSnr) = row;
end
otherAuditTable = struct2table(rowList);
end

function normalizationAuditTable = localBuildNormalizationAuditTable(perfTable, crbTable, hessianProbeTable, runtime, config)
%LOCALBUILDNORMALIZATIONAUDITTABLE Compare per-sat data scale, CRB scale, and Hessian scale.

rowTemplate = struct('snrDb', NaN, 'scopeName', "", 'satMode', "", ...
  'satIdxGlobal', NaN, 'pilotEnergy', NaN, 'numObsMean', NaN, ...
  'noiseSigma2Generator', NaN, 'noiseSigma2EstimatorMedian', NaN, ...
  'rxSignalPowerMedian', NaN, 'pathGainAbsMedian', NaN, ...
  'weightInEstimator', NaN, 'weightInCrb', NaN, 'effectivePilotGainAbs', NaN, ...
  'detFimTrace', NaN, 'ampFimTrace', NaN, 'ampToDetFimTraceRatio', NaN, ...
  'noisyHessianTraceMedian', NaN, 'meanNoisyHessianTraceMedian', NaN, ...
  'hessianToFimTraceRatioMedian', NaN, 'meanHessianToFimTraceRatioMedian', NaN);
scopeNameList = ["ref"; "other"; "joint"];
displayNameList = ["SS-SF-DoA"; "Other-SF-DoA"; "MS-SF-DoA"];
satModeList = ["single"; "other"; "multi"];
satIdxList = [runtime.truth.refSatIdxGlobal; runtime.truth.otherSatIdxGlobal; NaN];
effectiveGainList = [runtime.truth.effectivePilotGainAbsRef; ...
  runtime.truth.effectivePilotGainAbsOther; median(runtime.truth.effectivePilotGainAbs, 'omitnan')];
rowList = repmat(rowTemplate, 0, 1);
for iSnr = 1:height(crbTable)
  snrDb = crbTable.snrDb(iSnr);
  for iScope = 1:numel(scopeNameList)
    row = rowTemplate;
    row.snrDb = snrDb;
    row.scopeName = scopeNameList(iScope);
    row.satMode = satModeList(iScope);
    row.satIdxGlobal = satIdxList(iScope);
    row.pilotEnergy = sum(abs(runtime.pilotWave(:)).^2);
    row.noiseSigma2Generator = runtime.pwrSource / (10^(snrDb / 10));
    row.effectivePilotGainAbs = effectiveGainList(iScope);
    switch row.scopeName
      case "ref"
        row.detFimTrace = crbTable.detFimTraceSingle(iSnr);
        row.ampFimTrace = crbTable.ampFimTraceSingle(iSnr);
      case "other"
        row.detFimTrace = crbTable.detFimTraceOther(iSnr);
        row.ampFimTrace = crbTable.ampFimTraceOther(iSnr);
      otherwise
        row.detFimTrace = crbTable.detFimTraceJoint(iSnr);
        row.ampFimTrace = crbTable.ampFimTraceJoint(iSnr);
    end
    row.ampToDetFimTraceRatio = localSafeRatio(row.ampFimTrace, row.detFimTrace);
    caseMask = perfTable.displayName == displayNameList(iScope) & perfTable.snrDb == snrDb;
    if any(caseMask)
      row.numObsMean = median(perfTable.numObsMean(caseMask), 'omitnan');
      row.noiseSigma2EstimatorMedian = median(perfTable.noiseVarEstMean(caseMask), 'omitnan');
      row.rxSignalPowerMedian = median(perfTable.rxSignalPowerMean(caseMask), 'omitnan');
      row.pathGainAbsMedian = median(perfTable.pathGainAbsMean(caseMask), 'omitnan');
    end
    row.weightInEstimator = 1;
    row.weightInCrb = 1;
    if ~isempty(hessianProbeTable) && istable(hessianProbeTable) && ...
        ismember('snrDb', hessianProbeTable.Properties.VariableNames)
      noisyMask = hessianProbeTable.snrDb == snrDb & ...
        hessianProbeTable.scopeName == row.scopeName & hessianProbeTable.sourceTag == "noisy";
      meanMask = hessianProbeTable.snrDb == snrDb & ...
        hessianProbeTable.scopeName == row.scopeName & hessianProbeTable.sourceTag == "mean-noisy";
      if any(noisyMask)
        row.noisyHessianTraceMedian = median(hessianProbeTable.hessianTrace(noisyMask), 'omitnan');
        row.hessianToFimTraceRatioMedian = median(hessianProbeTable.hessianToDetFimTraceRatio(noisyMask), 'omitnan');
      end
      if any(meanMask)
        row.meanNoisyHessianTraceMedian = median(hessianProbeTable.hessianTrace(meanMask), 'omitnan');
        row.meanHessianToFimTraceRatioMedian = median(hessianProbeTable.hessianToDetFimTraceRatio(meanMask), 'omitnan');
      end
    end
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
normalizationAuditTable = struct2table(rowList);
end

function alphaWeightSweepTable = localBuildAlphaWeightSweepTable(alphaSweepCaseTable, crbTable, aggregateTable)
%LOCALBUILDALPHAWEIGHTSWEEPTABLE Aggregate probe-only alpha-sweep diagnostics.

rowTemplate = struct('snrDb', NaN, 'alphaSat2', NaN, 'numProbe', NaN, ...
  'rmseDeg', NaN, 'rmseOverJointDetCrb', NaN, 'realizedGainVsRef', NaN, ...
  'damageRateVsRef', NaN, 'medianObjective', NaN);
if isempty(alphaSweepCaseTable) || height(alphaSweepCaseTable) == 0
  alphaWeightSweepTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(alphaSweepCaseTable(:, {'snrDb', 'alphaSat2'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iGroup = 1:height(groupKey)
  snrDb = groupKey.snrDb(iGroup);
  alpha = groupKey.alphaSat2(iGroup);
  mask = alphaSweepCaseTable.snrDb == snrDb & alphaSweepCaseTable.alphaSat2 == alpha;
  subTable = alphaSweepCaseTable(mask, :);
  crbIdx = find(crbTable.snrDb == snrDb, 1, 'first');
  jointCrb = NaN;
  if ~isempty(crbIdx)
    jointCrb = crbTable.angleCrbDetJointDeg(crbIdx);
  end
  ssRmse = localFindAggregateMetric(aggregateTable, "SS-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  row = rowTemplate;
  row.snrDb = snrDb;
  row.alphaSat2 = alpha;
  row.numProbe = height(subTable);
  row.rmseDeg = localRmse(subTable.angleErrDeg);
  row.rmseOverJointDetCrb = localSafeRatio(row.rmseDeg, jointCrb);
  row.realizedGainVsRef = localSafeRatio(ssRmse, row.rmseDeg);
  row.damageRateVsRef = mean(abs(subTable.angleErrDeg) > ssRmse);
  row.medianObjective = median(subTable.objectiveAtFinal, 'omitnan');
  rowList(iGroup) = row;
end
alphaWeightSweepTable = struct2table(rowList);
end

function empiricalGainTable = localBuildEmpiricalGainTable(aggregateTable, crbTable)
%LOCALBUILDEMPIRICALGAINTABLE Predict joint gain from measured single-sat variances.

rowTemplate = struct('snrDb', NaN, 'ssVar', NaN, 'otherVar', NaN, ...
  'analyticSingleCrbVar', NaN, 'analyticOtherCrbVar', NaN, ...
  'analyticJointCrbVar', NaN, 'empiricalJointStdDeg', NaN, ...
  'expectedGainAnalytic', NaN, 'expectedGainEmpirical', NaN, ...
  'realizedGain', NaN);
if isempty(crbTable) || height(crbTable) == 0
  empiricalGainTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
rowList = repmat(rowTemplate, height(crbTable), 1);
for iSnr = 1:height(crbTable)
  snrDb = crbTable.snrDb(iSnr);
  ssRmse = localFindAggregateMetric(aggregateTable, "SS-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  otherRmse = localFindAggregateMetric(aggregateTable, "Other-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  msRmse = localFindAggregateMetric(aggregateTable, "MS-SF-DoA", snrDb, 'crbLocalAngleRmseDegDet');
  row = rowTemplate;
  row.snrDb = snrDb;
  row.ssVar = ssRmse.^2;
  row.otherVar = otherRmse.^2;
  row.analyticSingleCrbVar = crbTable.angleCrbDetSingleDeg(iSnr).^2;
  row.analyticOtherCrbVar = crbTable.angleCrbDetOtherDeg(iSnr).^2;
  row.analyticJointCrbVar = crbTable.angleCrbDetJointDeg(iSnr).^2;
  row.empiricalJointStdDeg = localCombineIndependentStd(ssRmse, otherRmse);
  row.expectedGainAnalytic = crbTable.jointDetGainOverSingle(iSnr);
  row.expectedGainEmpirical = localSafeRatio(ssRmse, row.empiricalJointStdDeg);
  row.realizedGain = localSafeRatio(ssRmse, msRmse);
  rowList(iSnr) = row;
end
empiricalGainTable = struct2table(rowList);
end

function jointStd = localCombineIndependentStd(stdA, stdB)
%LOCALCOMBINEINDEPENDENTSTD Combine two independent scalar standard deviations.

if ~isfinite(stdA) || ~isfinite(stdB) || stdA <= 0 || stdB <= 0
  jointStd = NaN;
  return;
end
jointStd = 1 / sqrt(1 / stdA.^2 + 1 / stdB.^2);
end


function ampAuditTable = localBuildAmpAwareCrbAuditTable(aggregateTable, crbTable)
%LOCALBUILDAMPAWARECRBAUDITTABLE Compare unit-gain and effective-gain CRB diagnostics.

rowTemplate = struct('snrDb', NaN, ...
  'unitOtherCrbDeg', NaN, 'ampOtherCrbDeg', NaN, ...
  'otherFullRmseDeg', NaN, 'otherFullRmseOverUnitCrb', NaN, ...
  'otherFullRmseOverAmpCrb', NaN, 'unitJointCrbDeg', NaN, ...
  'ampJointCrbDeg', NaN, 'msFullRmseDeg', NaN, ...
  'msFullRmseOverUnitCrb', NaN, 'msFullRmseOverAmpCrb', NaN, ...
  'expectedGainUnit', NaN, 'expectedGainAmp', NaN, ...
  'realizedGainVsRef', NaN, 'ampOtherScaleVsUnit', NaN, ...
  'ampJointScaleVsUnit', NaN);
if isempty(crbTable) || height(crbTable) == 0
  ampAuditTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
rowList = repmat(rowTemplate, height(crbTable), 1);
for iSnr = 1:height(crbTable)
  snrDb = crbTable.snrDb(iSnr);
  ssRmse = localFindAggregateMetric(aggregateTable, "SS-SF-DoA", snrDb, 'fullAngleRmseDeg');
  otherRmse = localFindAggregateMetric(aggregateTable, "Other-SF-DoA", snrDb, 'fullAngleRmseDeg');
  msRmse = localFindAggregateMetric(aggregateTable, "MS-SF-DoA", snrDb, 'fullAngleRmseDeg');
  row = rowTemplate;
  row.snrDb = snrDb;
  row.unitOtherCrbDeg = crbTable.angleCrbDetOtherDeg(iSnr);
  row.ampOtherCrbDeg = crbTable.angleCrbAmpOtherDeg(iSnr);
  row.otherFullRmseDeg = otherRmse;
  row.otherFullRmseOverUnitCrb = localSafeRatio(otherRmse, row.unitOtherCrbDeg);
  row.otherFullRmseOverAmpCrb = localSafeRatio(otherRmse, row.ampOtherCrbDeg);
  row.unitJointCrbDeg = crbTable.angleCrbDetJointDeg(iSnr);
  row.ampJointCrbDeg = crbTable.angleCrbAmpJointDeg(iSnr);
  row.msFullRmseDeg = msRmse;
  row.msFullRmseOverUnitCrb = localSafeRatio(msRmse, row.unitJointCrbDeg);
  row.msFullRmseOverAmpCrb = localSafeRatio(msRmse, row.ampJointCrbDeg);
  row.expectedGainUnit = crbTable.jointDetGainOverSingle(iSnr);
  row.expectedGainAmp = localSafeRatio(crbTable.angleCrbAmpSingleDeg(iSnr), row.ampJointCrbDeg);
  row.realizedGainVsRef = localSafeRatio(ssRmse, msRmse);
  row.ampOtherScaleVsUnit = localSafeRatio(row.ampOtherCrbDeg, row.unitOtherCrbDeg);
  row.ampJointScaleVsUnit = localSafeRatio(row.ampJointCrbDeg, row.unitJointCrbDeg);
  rowList(iSnr) = row;
end
ampAuditTable = struct2table(rowList);
end

function hessianAggregateTable = localBuildHessianAggregateTable(hessianProbeTable)
%LOCALBUILDHESSIANAGGREGATETABLE Aggregate per-seed Hessian diagnostics for printing.

rowTemplate = struct('snrDb', NaN, 'sourceTag', "", 'scopeName', "", ...
  'numProbe', NaN, 'medianHessianToDetFimTraceRatio', NaN, ...
  'medianHessianToDetFimShapeErr', NaN, 'medianHessianTrace', NaN, ...
  'medianDetFimTrace', NaN, 'medianHessianCond', NaN);
if isempty(hessianProbeTable) || height(hessianProbeTable) == 0
  hessianAggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
groupKey = unique(hessianProbeTable(:, {'snrDb', 'sourceTag', 'scopeName'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iGroup = 1:height(groupKey)
  mask = hessianProbeTable.snrDb == groupKey.snrDb(iGroup) & ...
    hessianProbeTable.sourceTag == groupKey.sourceTag(iGroup) & ...
    hessianProbeTable.scopeName == groupKey.scopeName(iGroup);
  row = rowTemplate;
  row.snrDb = groupKey.snrDb(iGroup);
  row.sourceTag = groupKey.sourceTag(iGroup);
  row.scopeName = groupKey.scopeName(iGroup);
  row.numProbe = nnz(mask);
  row.medianHessianToDetFimTraceRatio = median(hessianProbeTable.hessianToDetFimTraceRatio(mask), 'omitnan');
  row.medianHessianToDetFimShapeErr = median(hessianProbeTable.hessianToDetFimShapeErr(mask), 'omitnan');
  row.medianHessianTrace = median(hessianProbeTable.hessianTrace(mask), 'omitnan');
  row.medianDetFimTrace = median(hessianProbeTable.detFimTrace(mask), 'omitnan');
  row.medianHessianCond = median(hessianProbeTable.hessianCond(mask), 'omitnan');
  rowList(iGroup) = row;
end
hessianAggregateTable = struct2table(rowList);
end

function plotData = localBuildPlotData(aggregateTable, crbTable, crbCompareTable, ...
  truthInitDiffTable, doaGainTable, otherDoaCrbAuditTable, normalizationAuditTable, ...
  alphaWeightSweepTable, empiricalGainTable, ampAwareCrbAuditTable, noiseFreeCaseTable)
%LOCALBUILDPLOTDATA Store lightweight plot data.

plotData = struct();
plotData.aggregateTable = aggregateTable;
plotData.crbTable = crbTable;
plotData.crbCompareTable = crbCompareTable;
plotData.truthInitDiffTable = truthInitDiffTable;
plotData.doaGainTable = doaGainTable;
plotData.otherDoaCrbAuditTable = otherDoaCrbAuditTable;
plotData.normalizationAuditTable = normalizationAuditTable;
plotData.alphaWeightSweepTable = alphaWeightSweepTable;
plotData.empiricalGainTable = empiricalGainTable;
plotData.ampAwareCrbAuditTable = ampAwareCrbAuditTable;
plotData.noiseFreeCaseTable = noiseFreeCaseTable;
end

function plotData = localPlotReplay(replayData)
%LOCALPLOTREPLAY Draw compact diagnostic figures.

plotData = replayData.plotData;
aggregateTable = replayData.aggregateTable;
if isempty(aggregateTable) || height(aggregateTable) == 0
  return;
end
methodList = unique(aggregateTable.displayName, 'stable');

figure();
hold on;
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = aggregateTable.displayName == methodName;
  plot(aggregateTable.snrDb(mask), aggregateTable.crbLocalAngleRmseOverDetCrb(mask), '-o');
end
hold off;
grid on;
yline(1, '--');
xlabel('SNR (dB)');
ylabel('CRB-local angle RMSE / DoA-only CRB');
title('SF DoA-only MLE vs deterministic DoA CRB');
legend(cellstr(methodList), 'Location', 'northeast');

figure();
subplot(1, 2, 1);
plot(replayData.crbCompareTable.snrDb, replayData.crbCompareTable.ddToDetJointRatio, '-o');
grid on;
yline(1, '--');
xlabel('SNR (dB)');
ylabel('DD CRB / DoA-only CRB');
title('MS CRB definition gap');

subplot(1, 2, 2);
plot(replayData.truthInitDiffTable.snrDb, replayData.truthInitDiffTable.baselineRmseOverDetCrb, '-o');
hold on;
plot(replayData.truthInitDiffTable.snrDb, replayData.truthInitDiffTable.truthInitRmseOverDetCrb, '-s');
hold off;
grid on;
yline(1, '--');
xlabel('SNR (dB)');
ylabel('RMSE / DoA-only CRB');
title('MS baseline vs truth-init');
legend({'baseline', 'truth-init'}, 'Location', 'northeast');
end


function gainAbs = localEstimateTruthPilotGainAbs(view, rxSig, pilotWave, wavelen, latlon)
%LOCALESTIMATETRUTHPILOTGAINABS Estimate no-Doppler pilot-model gain at truth.

[arrayCell, numArray] = localWrapArray(view.sceneRef.array);
rxCell = localWrapRxSig(rxSig, numArray);
doaGridCell = localWrapGrid(view.doaGrid, numArray);
pilotPad = localBuildPilotPad(pilotWave, size(rxCell{1}, 2));
gainAbs = nan(numArray, 1);
for iArray = 1:numArray
  localDoa = localLatlonToAngle(doaGridCell{iArray}, latlon(1), latlon(2));
  steeringVec = steeringMatrix(arrayCell{iArray}, wavelen, localDoa(:, 1));
  if size(steeringVec, 2) ~= 1
    steeringVec = steeringVec(:, 1);
  end
  atomMat = steeringVec * pilotPad;
  hVec = atomMat(:);
  yVec = rxCell{iArray}(:);
  atomEnergy = real(hVec' * hVec);
  if atomEnergy > eps && isfinite(atomEnergy)
    gainAbs(iArray) = abs((hVec' * yVec) / atomEnergy);
  end
end
end

function value = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlon, useLogObjective, weightVec)
%LOCALEVALUATEDOAONLYOBJECTIVE Evaluate the same concentrated DoA-only objective.

value = NaN;
if nargin < 7 || isempty(weightVec)
  weightVec = [];
end
if numel(latlon) ~= 2 || any(~isfinite(latlon(:)))
  return;
end
[arrayCell, numArray] = localWrapArray(view.sceneRef.array);
rxCell = localWrapRxSig(rxSig, numArray);
doaGridCell = localWrapGrid(view.doaGrid, numArray);
if isempty(weightVec)
  weightVec = ones(numArray, 1);
else
  weightVec = reshape(double(weightVec), [], 1);
  if numel(weightVec) ~= numArray
    error('replaySfMsDoaCrbDiagnose:WeightCountMismatch', ...
      'weightVec must have one entry per array.');
  end
end
pilotPad = localBuildPilotPad(pilotWave, size(rxCell{1}, 2));
numObs = zeros(numArray, 1);
noiseVar = zeros(numArray, 1);
residualNorm = 0;
for iArray = 1:numArray
  localDoa = localLatlonToAngle(doaGridCell{iArray}, latlon(1), latlon(2));
  steeringVec = steeringMatrix(arrayCell{iArray}, wavelen, localDoa(:, 1));
  if size(steeringVec, 2) ~= 1
    steeringVec = steeringVec(:, 1);
  end
  atomMat = steeringVec * pilotPad;
  hVec = atomMat(:);
  yVec = rxCell{iArray}(:);
  atomEnergy = real(hVec' * hVec);
  if atomEnergy <= eps || ~isfinite(atomEnergy)
    value = inf;
    return;
  end
  pathGain = (hVec' * yVec) / atomEnergy;
  residualVec = yVec - pathGain * hVec;
  currentResidual = max(real(residualVec' * residualVec), eps);
  numObs(iArray) = numel(yVec);
  noiseVar(iArray) = currentResidual / numObs(iArray);
  residualNorm = residualNorm + weightVec(iArray) * currentResidual;
end
if useLogObjective
  value = sum(weightVec(:) .* numObs .* log(max(noiseVar, eps)));
else
  value = residualNorm;
end
end

function H = localNumericalDoaHessian(view, rxSig, pilotWave, wavelen, latlonCenter, stepDeg)
%LOCALNUMERICALDOAHESSIAN Central-difference Hessian of the DoA-only objective.

latlonCenter = reshape(latlonCenter, [], 1);
stepDeg = reshape(stepDeg, [], 1);
if numel(stepDeg) == 1
  stepDeg = repmat(stepDeg, 2, 1);
end
H = zeros(2, 2);
f0 = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter, true);
for iDim = 1:2
  ei = zeros(2, 1);
  ei(iDim) = stepDeg(iDim);
  fPlus = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter + ei, true);
  fMinus = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter - ei, true);
  H(iDim, iDim) = (fPlus - 2 * f0 + fMinus) / (stepDeg(iDim)^2);
end
for iDim = 1:2
  for jDim = (iDim + 1):2
    ei = zeros(2, 1); ej = zeros(2, 1);
    ei(iDim) = stepDeg(iDim);
    ej(jDim) = stepDeg(jDim);
    fpp = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter + ei + ej, true);
    fpm = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter + ei - ej, true);
    fmp = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter - ei + ej, true);
    fmm = localEvaluateDoaOnlyObjective(view, rxSig, pilotWave, wavelen, latlonCenter - ei - ej, true);
    H(iDim, jDim) = (fpp - fpm - fmp + fmm) / (4 * stepDeg(iDim) * stepDeg(jDim));
    H(jDim, iDim) = H(iDim, jDim);
  end
end
end

function localDoaJac = localFiniteDiffLocalDoaJac(doaGrid, latlonCenter)
%LOCALFINITEDIFFLOCALDOAJAC Finite-difference local-DoA Jacobian wrt lat/lon degrees.

latlonCenter = reshape(latlonCenter, [], 1);
stepDeg = [1e-5; 1e-5];
localDoaJac = zeros(2, 2);
for iDim = 1:2
  stepVec = zeros(2, 1);
  stepVec(iDim) = stepDeg(iDim);
  doaPlus = localLatlonToAngle(doaGrid, latlonCenter(1) + stepVec(1), ...
    latlonCenter(2) + stepVec(2));
  doaMinus = localLatlonToAngle(doaGrid, latlonCenter(1) - stepVec(1), ...
    latlonCenter(2) - stepVec(2));
  localDoaJac(:, iDim) = (doaPlus(:) - doaMinus(:)) / (2 * stepDeg(iDim));
end
end

function doa = localLatlonToAngle(doaGrid, lat, lon)
%LOCALLATLONTOANGLE Map a shared lat/lon hypothesis to local DOA.

globalFrame = lower(string(doaGrid.globalFrame));
switch globalFrame
  case "ecef"
    [x, y, z] = geodetic2ecef(doaGrid.spheroid, lat, lon, 0);
    posGlob = [x; y; z];
    doa = ecefToAngleGrid(posGlob, doaGrid.arrayCenter);
  case "eci"
    utcMat = localExpandUtc(doaGrid.utc, 1);
    posGlob = lla2eci([lat, lon, 0], utcMat).';
    doa = globalToLocalDoa(posGlob, doaGrid.arrayCenter, doaGrid.rotMat);
  otherwise
    error('replaySfMsDoaCrbDiagnose:InvalidGlobalFrame', ...
      'Unsupported global frame: %s.', globalFrame);
end
end

function utcMat = localExpandUtc(utc, numPoint)
%LOCALEXPUTC Expand utc to an N-by-6 numeric date matrix.

if isa(utc, 'datetime')
  utc = datevec(utc);
end
if isnumeric(utc) && isequal(size(utc), [1, 6])
  utcMat = repmat(utc, numPoint, 1);
elseif isnumeric(utc) && size(utc, 2) == 6 && size(utc, 1) == numPoint
  utcMat = utc;
else
  error('replaySfMsDoaCrbDiagnose:InvalidUtc', ...
    'utc must be a datetime scalar, 1x6 date vector, or N-by-6 date matrix.');
end
end

function relErr = localFimAdditivityRelErr(detAuxJoint)
%LOCALFIMADDITIVITYRELERR Check MS FIM equals the sum of per-sat FIMs.

relErr = NaN;
if ~isfield(detAuxJoint, 'fimCell') || isempty(detAuxJoint.fimCell)
  return;
end
fimSum = zeros(size(detAuxJoint.fim));
for iCell = 1:numel(detAuxJoint.fimCell)
  fimSum = fimSum + detAuxJoint.fimCell{iCell};
end
relErr = norm(fimSum - detAuxJoint.fim, 'fro') / max(norm(detAuxJoint.fim, 'fro'), eps);
end

function shapeErr = localShapeError(A, B)
%LOCALSHAPEERROR Compare two positive-semidefinite matrices after trace scaling.

shapeErr = NaN;
if any(~isfinite(A(:))) || any(~isfinite(B(:))) || trace(A) <= 0 || trace(B) <= 0
  return;
end
As = A / trace(A);
Bs = B / trace(B);
shapeErr = norm(As - Bs, 'fro') / max(norm(Bs, 'fro'), eps);
end

function value = localTableCellValue(tableVar, rowIdx)
%LOCALTABLECELLVALUE Read a table variable that may store matrices in cells.

if iscell(tableVar)
  value = tableVar{rowIdx};
else
  value = tableVar(rowIdx, :);
end
end

function tableOut = localCollectTaskTable(taskOutCell, fieldName)
%LOCALCOLLECTTASKTABLE Concatenate a table field from task outputs.

tableCell = cell(numel(taskOutCell), 1);
for iTask = 1:numel(taskOutCell)
  if isstruct(taskOutCell{iTask}) && isfield(taskOutCell{iTask}, fieldName)
    tableCell{iTask} = taskOutCell{iTask}.(fieldName);
  else
    tableCell{iTask} = table();
  end
end
tableOut = vertcat(tableCell{:});
end

function [arrayCell, numArray] = localWrapArray(array)
%LOCALWRAPARRAY Normalize array input to a cell array.

if iscell(array)
  arrayCell = reshape(array, 1, []);
elseif isstruct(array)
  if isscalar(array)
    arrayCell = {array};
  else
    arrayCell = num2cell(reshape(array, 1, []));
  end
else
  error('replaySfMsDoaCrbDiagnose:InvalidArrayType', ...
    'array must be a struct, struct array, or cell array.');
end
numArray = numel(arrayCell);
end

function rxCell = localWrapRxSig(rxSig, numArray)
%LOCALWRAPRXSIG Normalize received snapshots to a cell array.

if iscell(rxSig)
  rxCell = reshape(rxSig, 1, []);
elseif isnumeric(rxSig)
  rxCell = {rxSig};
else
  error('replaySfMsDoaCrbDiagnose:InvalidRxSigType', ...
    'rxSig must be numeric or a cell array.');
end
if numel(rxCell) ~= numArray
  error('replaySfMsDoaCrbDiagnose:RxSigCountMismatch', ...
    'rxSig count must match array count.');
end
end

function gridCell = localWrapGrid(doaGrid, numArray)
%LOCALWRAPGRID Normalize DoA grid input to a cell array.

if iscell(doaGrid)
  gridCell = reshape(doaGrid, 1, []);
else
  gridCell = {doaGrid};
end
if numel(gridCell) ~= numArray
  error('replaySfMsDoaCrbDiagnose:GridCountMismatch', ...
    'DoA grid count must match array count.');
end
end

function pilotPad = localBuildPilotPad(pilotWave, numSample)
%LOCALBUILDPILOTPAD Build row pilot padded to the received length.

pilotRow = reshape(pilotWave, 1, []);
if numel(pilotRow) > numSample
  error('replaySfMsDoaCrbDiagnose:PilotTooLong', ...
    'pilotWave length must not exceed the received snapshot length.');
end
pilotPad = zeros(1, numSample, 'like', pilotRow);
pilotPad(1:numel(pilotRow)) = pilotRow;
end

function value = localMeanRxSignalPower(rxSig)
%LOCALMEANRXSIGNALPOWER Mean received signal power over one or more arrays.

rxCell = localWrapRxAny(rxSig);
powerVec = zeros(numel(rxCell), 1);
for iCell = 1:numel(rxCell)
  currentRx = rxCell{iCell};
  powerVec(iCell) = mean(abs(currentRx(:)).^2, 'omitnan');
end
value = mean(powerVec, 'omitnan');
end

function value = localMeanNumObs(rxSig)
%LOCALMEANNUMOBS Mean observation count over one or more arrays.

rxCell = localWrapRxAny(rxSig);
numObs = zeros(numel(rxCell), 1);
for iCell = 1:numel(rxCell)
  numObs(iCell) = numel(rxCell{iCell});
end
value = mean(numObs, 'omitnan');
end

function rxCell = localWrapRxAny(rxSig)
%LOCALWRAPRXANY Normalize received snapshots without requiring array count.

if iscell(rxSig)
  rxCell = reshape(rxSig, 1, []);
elseif isnumeric(rxSig)
  rxCell = {rxSig};
else
  error('replaySfMsDoaCrbDiagnose:InvalidRxSigType', ...
    'rxSig must be numeric or a cell array.');
end
end

function value = localRmse(x)
%LOCALRMSE Root mean square with finite filtering.

x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = sqrt(mean(x(:).^2));
end
end

function value = localPercentile(x, pct)
%LOCALPERCENTILE Percentile with finite filtering.

x = sort(x(isfinite(x)));
if isempty(x)
  value = NaN;
  return;
end
idx = 1 + (numel(x) - 1) * pct / 100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
  value = x(lo);
else
  value = x(lo) + (idx - lo) * (x(hi) - x(lo));
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide only when denominator is finite and positive.

if isfinite(numerator) && isfinite(denominator) && denominator > 0
  ratio = numerator / denominator;
else
  ratio = NaN;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Check whether task-level parfor can be used.

useParfor = false;
if numTask < 2 || ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  useParfor = ~isempty(ver('parallel'));
catch
  useParfor = false;
end
end

function tableOut = localFilterHessianAggregateForPrint(tableIn)
%LOCALFILTERHESSIANAGGREGATEFORPRINT Keep the command-line Hessian summary compact.

if isempty(tableIn) || height(tableIn) == 0
  tableOut = tableIn;
  return;
end
if any(strcmp(tableIn.Properties.VariableNames, 'sourceTag'))
  tableOut = tableIn(tableIn.sourceTag == "mean-noisy", :);
  if height(tableOut) > 0
    return;
  end
end
tableOut = localTablePreview(tableIn, 12);
end

function tableOut = localProbePreviewBySnrSeed(tableIn, maxSeedPerSnr)
%LOCALPROBEPREVIEWBYSNRSEED Return one compact probe seed per SNR by default.

if isempty(tableIn) || height(tableIn) == 0
  tableOut = tableIn;
  return;
end
if ~(any(strcmp(tableIn.Properties.VariableNames, 'snrDb')) && ...
    any(strcmp(tableIn.Properties.VariableNames, 'taskSeed')))
  tableOut = localTablePreview(tableIn, 12);
  return;
end
snrList = unique(tableIn.snrDb, 'stable');
keepMask = false(height(tableIn), 1);
for iSnr = 1:numel(snrList)
  snrMask = tableIn.snrDb == snrList(iSnr);
  seedList = unique(tableIn.taskSeed(snrMask), 'stable');
  seedList = seedList(1:min(maxSeedPerSnr, numel(seedList)));
  keepMask = keepMask | (snrMask & ismember(tableIn.taskSeed, seedList));
end
tableOut = tableIn(keepMask, :);
end

function tableOut = localTablePreview(tableIn, maxRows)
%LOCALTABLEPREVIEW Return a head/tail preview of a table.

if isempty(tableIn) || height(tableIn) <= maxRows
  tableOut = tableIn;
  return;
end
headCount = ceil(maxRows / 2);
tailCount = floor(maxRows / 2);
tableOut = [tableIn(1:headCount, :); tableIn((end - tailCount + 1):end, :)];
end

function textValue = localFormatRow(value)
%LOCALFORMATROW Format a numeric vector for compact printing.

value = reshape(value, 1, []);
if isempty(value)
  textValue = '[]';
  return;
end
textValue = strjoin(arrayfun(@(x) sprintf('%.6g', x), value, 'UniformOutput', false), ', ');
end

function metricLineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build compact HTML-ready metric lines.

metricLineList = strings(0, 1);
agg = replayData.aggregateTable;
msMask = agg.displayName == "MS-SF-DoA";
otherMask = agg.displayName == "Other-SF-DoA";
if any(msMask)
  metricLineList(end + 1, 1) = "• MS baseline median RMSE/DetCRB=<code>" + ...
    string(sprintf('%.3f', median(agg.crbLocalAngleRmseOverDetCrb(msMask), 'omitnan'))) + "</code>";
end
if any(otherMask)
  metricLineList(end + 1, 1) = "• Other-sat median RMSE/DetCRB=<code>" + ...
    string(sprintf('%.3f', median(agg.crbLocalAngleRmseOverDetCrb(otherMask), 'omitnan'))) + "</code>";
end
if isfield(replayData, 'truthInitDiffTable') && height(replayData.truthInitDiffTable) > 0
  metricLineList(end + 1, 1) = "• max truth-init diff deg=<code>" + ...
    string(sprintf('%.3g', max(replayData.truthInitDiffTable.maxEstimateDiffDeg, [], 'omitnan'))) + "</code>";
end
if isfield(replayData, 'doaGainTable') && height(replayData.doaGainTable) > 0
  metricLineList(end + 1, 1) = "• median realized/expected gain=<code>" + ...
    string(sprintf('%.3f', median(replayData.doaGainTable.gainRealizationRateDet, 'omitnan'))) + "</code>";
end
if isfield(replayData, 'ampAwareCrbAuditTable') && height(replayData.ampAwareCrbAuditTable) > 0
  metricLineList(end + 1, 1) = "• other RMSE/ampCRB median=<code>" + ...
    string(sprintf('%.3f', median(replayData.ampAwareCrbAuditTable.otherFullRmseOverAmpCrb, 'omitnan'))) + "</code>";
end
if isfield(replayData, 'crbCompareTable') && height(replayData.crbCompareTable) > 0
  metricLineList(end + 1, 1) = "• median DD/Det CRB ratio=<code>" + ...
    string(sprintf('%.3f', median(replayData.crbCompareTable.ddToDetJointRatio, 'omitnan'))) + "</code>";
end
metricLineList(end + 1, 1) = "• snapshot stores case, CRB, gain, Hessian and top-tail tables";
end

function localPrintReplayConfig(config, truth)
%LOCALPRINTREPLAYCONFIG Print compact replay configuration.

localPrintHeaderLine('Replay scale', sprintf('%d seeds x %d SNR points', ...
  config.numRepeat, numel(config.snrDbList)));
localPrintHeaderLine('SNR list (dB)', localFormatRow(config.snrDbList));
localPrintHeaderLine('Base seed', sprintf('%d', config.baseSeed));
localPrintHeaderLine('Probe seed count per SNR', sprintf('%d', config.probeSeedCountPerSnr));
localPrintHeaderLine('Hessian mean repeats', sprintf('%d', config.hessianMeanRepeat));
localPrintHeaderLine('Alpha sweep (probe only)', localFormatRow(config.alphaSat2List));
localPrintHeaderLine('Truth-init search range', 'full DoA grid');
localPrintHeaderLine('Selected satellites (global)', localFormatRow(truth.selectedSatIdxGlobal));
end

function localPrintReplaySummaryConfig(replayData)
%LOCALPRINTREPLAYSUMMARYCONFIG Print load-safe summary configuration lines.

config = localGetFieldOrDefault(replayData, 'config', struct());
truth = localGetFieldOrDefault(replayData, 'truth', struct());
if isfield(config, 'snrDbList')
  localPrintHeaderLine('SNR list (dB)', localFormatRow(config.snrDbList));
end
if isfield(config, 'numRepeat')
  localPrintHeaderLine('Repeats per SNR', sprintf('%d', config.numRepeat));
end
if isfield(config, 'hessianMeanRepeat')
  localPrintHeaderLine('Hessian mean repeats', sprintf('%d', config.hessianMeanRepeat));
end
if isfield(config, 'alphaSat2List')
  localPrintHeaderLine('Alpha sweep (probe only)', localFormatRow(config.alphaSat2List));
end
localPrintHeaderLine('Truth-init search range', 'full DoA grid');
if isfield(truth, 'selectedSatIdxGlobal')
  localPrintHeaderLine('Selected satellites (global)', localFormatRow(truth.selectedSatIdxGlobal));
end
end

function localPrintHeaderLine(labelText, valueText)
%LOCALPRINTHEADERLINE Print one aligned replay header line.

fprintf('  %-32s : %s\n', char(string(labelText)), char(string(valueText)));
end
