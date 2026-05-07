%SCANMFCPIPINTOOTHPERFMAP Scan controlled in-tooth CP/IP performance.
% Dev scan script. Edit the configuration section directly before running.
% It compares CP-K / CP-U / IP-K / IP-U using the same multi-sat static
% DoA seed while constraining fdRef to a truth-centered half-tooth range.
% The scan isolates phase-tying behavior from full-flow tooth selection.

clear; close all; clc;

%% Scan configuration

scanName = "scanMfCpIpInToothPerfMap";
saveSnapshot = true;
notifyTelegramEnable = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
optVerbose = false;
baseSeed = 253;
numRepeat = 3;
snrDbList = [0; 5; 10; 15];
frameCountList = [8; 10; 20];
oracleFdHalfToothFraction = 0.49;
oracleFdRateHalfWidthHzPerSec = 1000;
staticLocalDoaHalfWidthDeg = [0.002; 0.002];
toothResidualTolHz = 50;

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
frameCountList = reshape(double(frameCountList), [], 1);
staticLocalDoaHalfWidthDeg = reshape(double(staticLocalDoaHalfWidthDeg), [], 1);
numRepeat = numel(seedList);

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.baseSeed = baseSeed;
scanConfig.numRepeat = numRepeat;
scanConfig.seedList = seedList;
scanConfig.snrDbList = snrDbList;
scanConfig.frameCountList = frameCountList;
scanConfig.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
scanConfig.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
scanConfig.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg;
scanConfig.toothResidualTolHz = toothResidualTolHz;
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.checkpointResume = logical(checkpointResume);
scanConfig.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);
scanConfig.optVerbose = logical(optVerbose);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
checkpointDir = "";
checkpointRunState = struct();
checkpointCleanupReport = struct([]);
scanData = struct();

try
  %% Build context and scan tasks

  taskList = localBuildTaskList(frameCountList, snrDbList, seedList);
  scanConfig.numTask = numel(taskList);
  useParfor = localCanUseParfor(numel(taskList));
  if scanConfig.checkpointEnable
    checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, taskList, useParfor);
    checkpointDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
  else
    checkpointOpt = struct();
  end
  printMfScanHeader(char(scanName), scanConfig, checkpointDir);
  fprintf('Preparing runtime cache...\n');
  scanRuntime = localBuildScanRuntime(scanConfig);

  %% Run scan batch

  if scanConfig.checkpointEnable
    pendingCount = localCountPendingCheckpointTasks(checkpointDir, numel(taskList));
    progressTracker = localCreateScanProgressTracker(pendingCount);
    progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>
    if progressTracker.enabled
      checkpointOpt.progressFcn = @(~) localAdvanceScanProgress(progressTracker, false);
    end
    checkpointRunState = runPerfTaskGridWithCheckpoint(taskList, scanRuntime, ...
      @localRunCheckpointTask, checkpointOpt);
    taskOutCell = checkpointRunState.resultCell;
    localCloseScanProgressTracker(progressTracker);
    clear progressCleanup;
  else
    taskOutCell = cell(numel(taskList), 1);
    progressTracker = localCreateScanProgressTracker(numel(taskList));
    progressCleanup = onCleanup(@() localCloseScanProgressTracker(progressTracker)); %#ok<NASGU>
    if useParfor && progressTracker.enabled && ~isempty(progressTracker.queue)
      progressQueue = progressTracker.queue;
      parfor iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
        send(progressQueue, 1);
      end
    elseif useParfor
      parfor iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
      end
    else
      for iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), scanRuntime);
        localAdvanceScanProgress(progressTracker, false);
      end
    end
    localCloseScanProgressTracker(progressTracker);
    clear progressCleanup;
  end

  [perfTable, repeatOutCell] = localCollectTaskOutput(taskOutCell);
  aggregateTable = localBuildAggregateTable(perfTable, toothResidualTolHz);
  if scanConfig.checkpointEnable && localGetFieldOrDefault(checkpointRunState, 'isComplete', false) && ...
      scanConfig.checkpointCleanupOnSuccess
    checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(checkpointRunState, ...
      struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
  end
  checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunState, ...
    checkpointDir, checkpointCleanupReport);

  %% Data storage

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.perfTable = perfTable;
  scanData.aggregateTable = aggregateTable;
  scanData.repeatOutCell = repeatOutCell;
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupReport = checkpointCleanupReport;
  scanData.plotData = localBuildPlotData(aggregateTable);
  scanData.elapsedSec = toc(runTic);

  if ~scanConfig.checkpointEnable
    scanData = finalizeMfScanResult(scanData, "");
  end

  if scanConfig.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output and plotting

  printMfScanSection('In-tooth CP/IP performance-map aggregate', scanData.aggregateTable);
  fprintf('Frames                          : %s\n', localFormatRow(scanData.config.frameCountList));
  fprintf('SNR list (dB)                   : %s\n', localFormatRow(scanData.config.snrDbList));
  fprintf('Repeats per config              : %d\n', scanData.config.numRepeat);
  fprintf('Runtime cache                   : context per frame.\n');
  fprintf('fdRef oracle half-tooth fraction: %.3f\n', scanData.config.oracleFdHalfToothFraction);
  fprintf('fdRate oracle half-width        : %.2f Hz/s\n', scanData.config.oracleFdRateHalfWidthHzPerSec);
  fprintf('Truth-tooth rule                : toothIdx == 0 and |residual| <= %.2f Hz\n', ...
    scanData.config.toothResidualTolHz);
  if height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  scanData.plotData = localPlotScan(scanData);

  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "DONE", ...
    'config', scanData.config, ...
    'snapshotFile', scanData.snapshotFile, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Controlled in-tooth CP/IP scan completed."; ...
      "Use this result as the cleaner CP/IP trade-off reference, not a full-flow proof."]));

catch ME
  if exist('progressTracker', 'var')
    localCloseScanProgressTracker(progressTracker);
  end
  if exist('progressCleanup', 'var')
    clear progressCleanup;
  end
  localPrintCheckpointFailureHint(checkpointDir);
  fprintf('Scan failed while building controlled in-tooth CP/IP performance data.\n');
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function taskOut = localRunCheckpointTask(task, scanRuntime)
%LOCALRUNCHECKPOINTTASK Run one checkpointed scan task.

taskOut = localRunOneTask(task, scanRuntime);
end

function checkpointOpt = localBuildCheckpointOpt(scanName, scanConfig, taskList, useParfor)
%LOCALBUILDCHECKPOINTOPT Build task-level checkpoint runner options.

checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, taskList);
checkpointOpt.outputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.useParfor = logical(useParfor);
checkpointOpt.resume = logical(scanConfig.checkpointResume);
checkpointOpt.meta = localBuildCheckpointMeta(scanConfig, taskList);
checkpointOpt.progressFcn = [];
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end

function runKey = localBuildCheckpointRunKey(scanConfig, taskList)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable key for interrupted in-tooth scans.

frameList = unique([taskList.numFrame], 'stable');
snrList = unique([taskList.snrDb], 'stable');
seedList = unique([taskList.taskSeed], 'stable');
runKey = sprintf('frames%s_snr%s_seed%dto%d_rep%d_fdfrac%.3f_ratehw%.0f_doahw%s', ...
  localJoinNumericToken(frameList), localJoinNumericToken(snrList), seedList(1), ...
  seedList(end), numel(seedList), scanConfig.oracleFdHalfToothFraction, ...
  scanConfig.oracleFdRateHalfWidthHzPerSec, ...
  localJoinNumericToken(scanConfig.staticLocalDoaHalfWidthDeg(:).'));
runKey = string(runKey);
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
runKey = replace(runKey, ' ', '');
runKey = replace(runKey, '+', '');
end

function meta = localBuildCheckpointMeta(scanConfig, taskList)
%LOCALBUILDCHECKPOINTMETA Store the task-defining scan signature.

meta = struct();
meta.scanName = string(scanConfig.scanName);
meta.frameCountList = reshape(unique([taskList.numFrame], 'stable'), 1, []);
meta.snrDbList = reshape(unique([taskList.snrDb], 'stable'), 1, []);
meta.seedList = reshape(unique([taskList.taskSeed], 'stable'), 1, []);
meta.oracleFdHalfToothFraction = scanConfig.oracleFdHalfToothFraction;
meta.oracleFdRateHalfWidthHzPerSec = scanConfig.oracleFdRateHalfWidthHzPerSec;
meta.staticLocalDoaHalfWidthDeg = reshape(scanConfig.staticLocalDoaHalfWidthDeg, 1, []);
meta.toothResidualTolHz = scanConfig.toothResidualTolHz;
meta.methodList = ["CP-K", "CP-U", "IP-K", "IP-U"];
end

function pendingCount = localCountPendingCheckpointTasks(checkpointDir, numTask)
%LOCALCOUNTPENDINGCHECKPOINTTASKS Count unfinished task checkpoint files.

taskDir = fullfile(checkpointDir, 'task');
if ~isfolder(taskDir)
  pendingCount = numTask;
  return;
end
doneCount = 0;
for iTask = 1:numTask
  taskFile = fullfile(taskDir, sprintf('task_%06d.mat', iTask));
  if isfile(taskFile)
    doneCount = doneCount + 1;
  end
end
pendingCount = max(0, numTask - doneCount);
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runState, checkpointDir, cleanupReport)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Summarize task checkpoint resume and cleanup state.

if ~isstruct(runState) || ~isfield(runState, 'numTask')
  checkpointSummaryTable = table();
  return;
end
cleanupApplied = ~isempty(cleanupReport);
checkpointSummaryTable = table(string(checkpointDir), runState.numTask, ...
  runState.numDone, logical(runState.isComplete), cleanupApplied, ...
  'VariableNames', {'expectedRunDir', 'numTask', 'numDone', 'isComplete', 'cleanupApplied'});
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print preserved checkpoint path after failures.

checkpointDir = string(checkpointDir);
if strlength(checkpointDir) == 0
  fprintf('Scan failed. No checkpoint directory was created.\n');
  return;
end
fprintf('Scan failed. Completed checkpoint tasks were kept under tmp for resume.\n');
fprintf('  checkpoint dir: %s\n', char(checkpointDir));
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve the repository root from this scan script location.

scanDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scanDir)));
end

function token = localJoinNumericToken(value)
%LOCALJOINNUMERICTOKEN Format numeric values for a file-system-safe run key.

value = reshape(double(value), 1, []);
token = char(strjoin(compose('%.6g', value), 'x'));
end

function taskList = localBuildTaskList(frameCountList, snrDbList, seedList)
%LOCALBUILDTASKLIST Build independent scan tasks for the outer batch loop.

taskTemplate = struct('numFrame', NaN, 'snrDb', NaN, 'taskSeed', NaN);
taskList = repmat(taskTemplate, numel(frameCountList) * numel(snrDbList) * numel(seedList), 1);
idxTask = 0;
for iP = 1:numel(frameCountList)
  for iSnr = 1:numel(snrDbList)
    for iSeed = 1:numel(seedList)
      idxTask = idxTask + 1;
      taskList(idxTask).numFrame = frameCountList(iP);
      taskList(idxTask).snrDb = snrDbList(iSnr);
      taskList(idxTask).taskSeed = seedList(iSeed);
    end
  end
end
end

function taskOut = localRunOneTask(task, scanRuntime)
%LOCALRUNONETASK Run one frame-count, SNR, and seed configuration.

scanConfig = scanRuntime.config;
context = localGetRuntimeContext(scanRuntime, task.numFrame);
flowOpt = localGetRuntimeFlowOpt(scanRuntime, task.numFrame);
methodList = scanRuntime.methodList;
repeatOut = localRunOneRepeat(context, flowOpt, methodList, scanConfig, task.snrDb, task.taskSeed);

taskOut = struct();
taskOut.caseRowList = repeatOut.caseRowList;
taskOut.repeatSlim = localStripRepeatOut(repeatOut);
end

function scanRuntime = localBuildScanRuntime(scanConfig)
%LOCALBUILDSCANRUNTIME Cache frame-level context outside task workers.

frameCountList = reshape(unique(scanConfig.frameCountList, 'stable'), [], 1);
contextBank = repmat(struct('numFrame', NaN, 'context', struct(), 'flowOpt', struct()), numel(frameCountList), 1);
parallelOpt = localBuildSerialInnerParallelOpt();
for iFrame = 1:numel(frameCountList)
  numFrame = frameCountList(iFrame);
  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = parallelOpt;
  context = buildDynamicDualSatEciContext(contextOpt);
  context = localDisableSubsetBankForInTooth(context);
  contextBank(iFrame).numFrame = numFrame;
  contextBank(iFrame).context = context;
  contextBank(iFrame).flowOpt = localBuildFlowOpt(scanConfig.staticLocalDoaHalfWidthDeg, context.parallelOpt);
end

scanRuntime = struct();
scanRuntime.config = scanConfig;
scanRuntime.contextBank = contextBank;
scanRuntime.methodList = localBuildMethodList(scanConfig.staticLocalDoaHalfWidthDeg);
end

function context = localGetRuntimeContext(scanRuntime, numFrame)
%LOCALGETRUNTIMECONTEXT Return the cached context for one frame count.

idx = find([scanRuntime.contextBank.numFrame] == numFrame, 1, 'first');
if isempty(idx)
  error('scanMfCpIpInToothPerfMap:MissingRuntimeContext', ...
    'No cached runtime context found for numFrame=%g.', numFrame);
end
context = scanRuntime.contextBank(idx).context;
end

function flowOpt = localGetRuntimeFlowOpt(scanRuntime, numFrame)
%LOCALGETRUNTIMEFLOWOPT Return the cached flow options for one frame count.

idx = find([scanRuntime.contextBank.numFrame] == numFrame, 1, 'first');
if isempty(idx)
  error('scanMfCpIpInToothPerfMap:MissingRuntimeFlowOpt', ...
    'No cached runtime flow options found for numFrame=%g.', numFrame);
end
flowOpt = scanRuntime.contextBank(idx).flowOpt;
end

function parallelOpt = localBuildSerialInnerParallelOpt()
%LOCALBUILDSERIALINNERPARALLELOPT Disable inner parallelism for outer parfor scans.

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
end

function context = localDisableSubsetBankForInTooth(context)
%LOCALDISABLESUBSETBANKFORINTOOTH Remove subset fixtures from this controlled scan.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function flowOpt = localBuildFlowOpt(staticLocalDoaHalfWidthDeg, parallelOpt)
%LOCALBUILDFLOWOPT Build static and dynamic options shared by in-tooth methods.

flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false, ...
  'staticMsHalfWidth', staticLocalDoaHalfWidthDeg(:)));
flowOpt = applyDynamicTransitionFlowDefaults(flowOpt);
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if numFrame < 2
  error('scanMfCpIpInToothPerfMap:InvalidFrameCount', ...
    'Frame count must be at least two for multi-frame CP/IP scanning.');
end
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function methodList = localBuildMethodList(staticLocalDoaHalfWidthDeg)
%LOCALBUILDMETHODLIST Define the controlled CP/IP in-tooth method bank.

methodList = repmat(localMakeMethod("", "", "", false, "", "", []), 0, 1);
methodList(end + 1, 1) = localMakeMethod("CP-K", "continuous", "known", true, ...
  "truthHalfTooth", "known", staticLocalDoaHalfWidthDeg);
methodList(end + 1, 1) = localMakeMethod("CP-U", "continuous", "unknown", false, ...
  "truthHalfTooth", "truthLocal", staticLocalDoaHalfWidthDeg);
methodList(end + 1, 1) = localMakeMethod("IP-K", "independent", "known", true, ...
  "truthHalfTooth", "known", staticLocalDoaHalfWidthDeg);
methodList(end + 1, 1) = localMakeMethod("IP-U", "independent", "unknown", false, ...
  "truthHalfTooth", "truthLocal", staticLocalDoaHalfWidthDeg);
end

function method = localMakeMethod(displayName, phaseMode, fdRateMode, isKnownRate, fdRangeMode, fdRateRangeMode, doaHalfWidthDeg)
%LOCALMAKEMETHOD Construct one in-tooth method descriptor.

method = struct();
method.displayName = string(displayName);
method.phaseMode = string(phaseMode);
method.fdRateMode = string(fdRateMode);
method.isKnownRate = logical(isKnownRate);
method.fdRangeMode = string(fdRangeMode);
method.fdRateRangeMode = string(fdRateRangeMode);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
end

function repeatOut = localRunOneRepeat(context, flowOpt, methodList, scanConfig, snrDb, taskSeed)
%LOCALRUNONEREPEAT Run one noisy repeat and collect CP/IP method rows.

lastwarn('', '');
tRepeat = tic;
repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('scanMfCpIpInToothPerfMap:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for the controlled in-tooth scan.');
end
fdHalfWidthHz = scanConfig.oracleFdHalfToothFraction * toothStepHz;
fdRangeOracle = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeOracle = truthFdRateHzPerSec + scanConfig.oracleFdRateHalfWidthHzPerSec * [-1, 1];

staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fixture.fdRange, truth, context.otherSatIdxGlobal, scanConfig.optVerbose, ...
  flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);

caseRowList = repmat(localEmptyCaseRow(), numel(methodList), 1);
for iMethod = 1:numel(methodList)
  method = methodList(iMethod);
  tMethod = tic;
  caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, context, flowOpt, ...
    fdRangeOracle, fdRateRangeOracle, toothStepHz, scanConfig.optVerbose);
  wallTimeMs = 1000 * toc(tMethod);
  caseRowList(iMethod) = localBuildCaseRow(method, caseUse, truth, snrDb, taskSeed, ...
    fixture, toothStepHz, wallTimeMs);
end
repeatTotalMs = 1000 * toc(tRepeat);
[warningMessage, warningId] = lastwarn;

repeatOut = struct();
repeatOut.taskSeed = taskSeed;
repeatOut.snrDb = snrDb;
repeatOut.truth = truth;
repeatOut.toothStepHz = toothStepHz;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.fdRangeOracle = fdRangeOracle;
repeatOut.fdRateRangeOracle = fdRateRangeOracle;
repeatOut.caseRowList = caseRowList;
repeatOut.repeatTotalMs = repeatTotalMs;
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, context, flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz, optVerbose)
%LOCALRUNDYNAMICMETHOD Run one controlled multi-sat dynamic method.

fdRangeUse = fixture.fdRange;
fdRateRangeUse = fixture.fdRateRange;
if method.fdRangeMode == "truthHalfTooth"
  fdRangeUse = fdRangeOracle;
end
if method.fdRateRangeMode == "truthLocal"
  fdRateRangeUse = fdRateRangeOracle;
end

truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
seedDoaParam = staticBundle.caseStaticMs.estResult.doaParamEst(:);

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
initParam = localBuildTruthFreqInitParam(seedDoaParam, truthFdRefHz, truthFdRateHzPerSec, method.isKnownRate);
initOverride = struct('startTag', lower(method.displayName) + "-truth-freq", ...
  'initParam', initParam, 'initDoaParam', seedDoaParam(:), ...
  'initDoaHalfWidth', method.doaHalfWidthDeg(:), 'freezeDoa', false);

caseUse = runDynamicDoaDopplerCase("MS-MF-" + method.displayName + "-in-tooth", "multi", ...
  fixture.viewMs, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
  localGetFieldOrDefault(fixture, 'debugTruthMs', struct()), initOverride);
caseUse.inToothMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, ...
  'fdRateRangeUse', fdRateRangeUse);
end

function initParam = localBuildTruthFreqInitParam(doaParam, fdRefHz, fdRateHzPerSec, isKnownRate)
%LOCALBUILDTRUTHFREQINITPARAM Build an oracle frequency initializer.

if isKnownRate
  initParam = [doaParam(:); fdRefHz];
else
  initParam = [doaParam(:); fdRefHz; fdRateHzPerSec];
end
end

function row = localBuildCaseRow(method, caseUse, truth, snrDb, taskSeed, fixture, toothStepHz, wallTimeMs)
%LOCALBUILDCASEMETRICROW Convert one estimator case into a compact scan row.

row = localEmptyCaseRow();
row.displayName = method.displayName;
row.phaseMode = method.phaseMode;
row.fdRateMode = method.fdRateMode;
row.snrDb = snrDb;
row.taskSeed = taskSeed;
row.numFrame = fixture.sceneSeq.numFrame;
row.frameIntvlSec = median(diff(fixture.sceneSeq.timeOffsetSec));
row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
row.toothStepHz = toothStepHz;
row.inToothMode = method.fdRangeMode;
row.fdRateRangeMode = method.fdRateRangeMode;
row.doaSeedMode = "static-ms";
row.initMode = "truth-freq";
info = summarizeDoaDopplerCase(caseUse, truth);
dynSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth);
row.isResolved = logical(localGetFieldOrDefault(info, 'isResolved', false));
row.angleErrDeg = localGetFieldOrDefault(info, 'angleErrDeg', NaN);
row.fdRefAbsErrHz = abs(localGetFieldOrDefault(info, 'fdRefErrHz', NaN));
row.fdRateAbsErrHzPerSec = abs(localGetFieldOrDefault(info, 'fdRateErrHzPerSec', NaN));
row.fdLineRmseHz = localGetFieldOrDefault(info, 'fdLineRmseHz', NaN);
row.fdSatRmseHz = localGetFieldOrDefault(info, 'fdSatRmseHz', NaN);
row.deltaFdRmseHz = localGetFieldOrDefault(info, 'deltaFdRmseHz', NaN);
row.toothIdx = localGetFieldOrDefault(dynSummary, 'toothIdx', NaN);
row.toothResidualHz = localGetFieldOrDefault(dynSummary, 'toothResidualHz', NaN);
row.nonRefCoherenceFloor = localGetFieldOrDefault(dynSummary, 'nonRefCoherenceFloor', NaN);
row.nonRefRmsPhaseResidRad = localGetFieldOrDefault(dynSummary, 'nonRefRmsPhaseResidRad', NaN);
row.solveVariant = string(localGetFieldOrDefault(dynSummary, 'solveVariant', ""));
row.runTimeMs = localGetFieldOrDefault(info, 'runTimeMs', NaN);
row.wallTimeMs = wallTimeMs;
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for per-repeat case metrics.

row = struct('displayName', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'taskSeed', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, ...
  'windowSec', NaN, 'toothStepHz', NaN, 'inToothMode', "", ...
  'fdRateRangeMode', "", 'doaSeedMode', "", 'initMode', "", ...
  'angleErrDeg', NaN, 'fdRefAbsErrHz', NaN, 'fdRateAbsErrHzPerSec', NaN, ...
  'fdLineRmseHz', NaN, 'fdSatRmseHz', NaN, 'deltaFdRmseHz', NaN, ...
  'toothIdx', NaN, 'toothResidualHz', NaN, 'nonRefCoherenceFloor', NaN, ...
  'nonRefRmsPhaseResidRad', NaN, 'isResolved', false, 'solveVariant', "", ...
  'runTimeMs', NaN, 'wallTimeMs', NaN);
end

function [perfTable, repeatOutCell] = localCollectTaskOutput(taskOutCell)
%LOCALCOLLECTTASKOUTPUT Merge per-task rows and slim repeat diagnostics.

rowList = repmat(localEmptyCaseRow(), 0, 1);
repeatOutCell = cell(numel(taskOutCell), 1);
for iTask = 1:numel(taskOutCell)
  taskOut = taskOutCell{iTask};
  if isempty(taskOut)
    continue;
  end
  rowList = [rowList; taskOut.caseRowList(:)]; %#ok<AGROW>
  repeatOutCell{iTask} = taskOut.repeatSlim;
end
perfTable = struct2table(rowList(:));
end

function aggregateTable = localBuildAggregateTable(perfTable, toothResidualTolHz)
%LOCALBUILDAGGREGATETABLE Aggregate in-tooth performance metrics by case and grid point.

groupKey = unique(perfTable(:, {'displayName', 'phaseMode', 'fdRateMode', 'snrDb', 'numFrame', 'frameIntvlSec'}), 'rows', 'stable');
rowList = repmat(localEmptyAggregateRow(), height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = perfTable.displayName == groupKey.displayName(iRow) & ...
    perfTable.snrDb == groupKey.snrDb(iRow) & ...
    perfTable.numFrame == groupKey.numFrame(iRow) & ...
    abs(perfTable.frameIntvlSec - groupKey.frameIntvlSec(iRow)) < eps;
  row = localEmptyAggregateRow();
  row.displayName = groupKey.displayName(iRow);
  row.phaseMode = groupKey.phaseMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.numFrame = groupKey.numFrame(iRow);
  row.frameIntvlSec = groupKey.frameIntvlSec(iRow);
  row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
  row.numRepeat = nnz(mask);
  row.resolvedRate = mean(double(perfTable.isResolved(mask)), 'omitnan');
  row.truthToothHitRate = mean(double(abs(perfTable.toothIdx(mask)) == 0 & ...
    abs(perfTable.toothResidualHz(mask)) <= toothResidualTolHz), 'omitnan');
  row.angleRmseDeg = sqrt(mean(perfTable.angleErrDeg(mask).^2, 'omitnan'));
  row.angleMedianDeg = median(perfTable.angleErrDeg(mask), 'omitnan');
  row.angleP95Deg = prctile(perfTable.angleErrDeg(mask), 95);
  row.fdRefRmseHz = sqrt(mean(perfTable.fdRefAbsErrHz(mask).^2, 'omitnan'));
  row.fdRefMedianHz = median(perfTable.fdRefAbsErrHz(mask), 'omitnan');
  row.fdRefP95Hz = prctile(perfTable.fdRefAbsErrHz(mask), 95);
  row.fdRateRmseHzPerSec = sqrt(mean(perfTable.fdRateAbsErrHzPerSec(mask).^2, 'omitnan'));
  row.nonRefCoherenceFloorMedian = median(perfTable.nonRefCoherenceFloor(mask), 'omitnan');
  row.meanRunTimeMs = mean(perfTable.wallTimeMs(mask), 'omitnan');
  rowList(iRow) = row;
end
aggregateTable = struct2table(rowList);
aggregateTable = localAddCpIpRatios(aggregateTable);
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return a typed empty row for aggregate metrics.

row = struct('displayName', "", 'phaseMode', "", 'fdRateMode', "", ...
  'snrDb', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, 'windowSec', NaN, ...
  'numRepeat', NaN, 'resolvedRate', NaN, 'truthToothHitRate', NaN, ...
  'angleRmseDeg', NaN, 'angleMedianDeg', NaN, 'angleP95Deg', NaN, ...
  'fdRefRmseHz', NaN, 'fdRefMedianHz', NaN, 'fdRefP95Hz', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'nonRefCoherenceFloorMedian', NaN, ...
  'meanRunTimeMs', NaN, 'angleRatioToCp', NaN, 'fdRatioToCp', NaN);
end

function aggregateTable = localAddCpIpRatios(aggregateTable)
%LOCALADDCPIPRATIOS Add IP-to-CP ratios for matching SNR and frame-count rows.

aggregateTable.angleRatioToCp = nan(height(aggregateTable), 1);
aggregateTable.fdRatioToCp = nan(height(aggregateTable), 1);
for iRow = 1:height(aggregateTable)
  if aggregateTable.phaseMode(iRow) ~= "independent"
    continue;
  end
  cpName = "CP-K";
  if aggregateTable.fdRateMode(iRow) == "unknown"
    cpName = "CP-U";
  end
  maskCp = aggregateTable.displayName == cpName & aggregateTable.snrDb == aggregateTable.snrDb(iRow) & ...
    aggregateTable.numFrame == aggregateTable.numFrame(iRow) & ...
    abs(aggregateTable.frameIntvlSec - aggregateTable.frameIntvlSec(iRow)) < eps;
  idxCp = find(maskCp, 1, 'first');
  if ~isempty(idxCp)
    aggregateTable.angleRatioToCp(iRow) = aggregateTable.angleRmseDeg(iRow) / aggregateTable.angleRmseDeg(idxCp);
    aggregateTable.fdRatioToCp(iRow) = aggregateTable.fdRefRmseHz(iRow) / aggregateTable.fdRefRmseHz(idxCp);
  end
end
end

function plotData = localBuildPlotData(aggregateTable)
%LOCALBUILDPLOTDATA Build lightweight tables used to redraw scan figures.

plotData = struct();
plotData.figureFiles = strings(0, 1);
plotData.angleVsSnrTable = localBuildAngleVsSnrPlotTable(aggregateTable);
plotData.fdRefVsSnrTable = localBuildFdRefVsSnrPlotTable(aggregateTable);
plotData.angleRatioTable = aggregateTable(aggregateTable.phaseMode == "independent", :);
plotData.angleVsFrameTable = localBuildAngleVsFramePlotTable(aggregateTable);
end

function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw the default controlled in-tooth CP/IP figures.

plotData = localBuildPlotData(scanData.aggregateTable);
angleTable = plotData.angleVsSnrTable;
caseList = unique(angleTable.displayName, 'stable');
figure('Name', 'In-tooth CP/IP angle RMSE vs SNR'); hold on;
for iCase = 1:numel(caseList)
  mask = angleTable.displayName == caseList(iCase);
  plot(angleTable.snrDb(mask), angleTable.angleRmseDeg(mask), '-o');
end
grid on; xlabel('SNR (dB)'); ylabel('angle RMSE (deg)'); legend(cellstr(caseList), 'Location', 'best');
drawnow;

ratioTable = plotData.angleRatioTable;
figure('Name', 'In-tooth IP/CP angle ratio');
scatter(ratioTable.windowSec * 1e3, ratioTable.angleRatioToCp, 60, ratioTable.snrDb, 'filled');
grid on; xlabel('window span (ms)'); ylabel('IP / CP angle RMSE');
drawnow;

fdTable = plotData.fdRefVsSnrTable;
caseList = unique(fdTable.displayName, 'stable');
figure('Name', 'In-tooth CP/IP fdRef RMSE vs SNR'); hold on;
for iCase = 1:numel(caseList)
  mask = fdTable.displayName == caseList(iCase);
  plot(fdTable.snrDb(mask), fdTable.fdRefRmseHz(mask), '-o');
end
grid on; xlabel('SNR (dB)'); ylabel('fdRef RMSE (Hz)'); legend(cellstr(caseList), 'Location', 'best');
drawnow;

frameTable = plotData.angleVsFrameTable;
caseList = unique(frameTable.displayName, 'stable');
figure('Name', 'In-tooth CP/IP angle RMSE vs frame count'); hold on;
for iCase = 1:numel(caseList)
  mask = frameTable.displayName == caseList(iCase);
  plot(frameTable.numFrame(mask), frameTable.angleRmseDeg(mask), '-o');
end
grid on; xlabel('number of frames'); ylabel('angle RMSE (deg)'); legend(cellstr(caseList), 'Location', 'best');
drawnow;
end

function plotTable = localBuildAngleVsSnrPlotTable(aggregateTable)
%LOCALBUILDANGLEVSSNRPLOTTABLE Select max-frame rows for the angle-vs-SNR figure.

if isempty(aggregateTable)
  plotTable = aggregateTable;
  return;
end
maxFrame = max(aggregateTable.numFrame);
plotTable = aggregateTable(aggregateTable.numFrame == maxFrame, :);
end

function plotTable = localBuildFdRefVsSnrPlotTable(aggregateTable)
%LOCALBUILDFDREFVSSNRPLOTTABLE Select max-frame rows for the fdRef-vs-SNR figure.

if isempty(aggregateTable)
  plotTable = aggregateTable;
  return;
end
maxFrame = max(aggregateTable.numFrame);
plotTable = aggregateTable(aggregateTable.numFrame == maxFrame, :);
end

function plotTable = localBuildAngleVsFramePlotTable(aggregateTable)
%LOCALBUILDANGLEVSFRAMEPLOTTABLE Select max-SNR rows for the frame-count figure.

if isempty(aggregateTable)
  plotTable = aggregateTable;
  return;
end
maxSnr = max(aggregateTable.snrDb);
plotTable = aggregateTable(aggregateTable.snrDb == maxSnr, :);
end

function repeatSlim = localStripRepeatOut(repeatOut)
%LOCALSTRIPREPEATOUT Keep only lightweight repeat diagnostics for scanData.

repeatSlim = struct();
repeatSlim.taskSeed = repeatOut.taskSeed;
repeatSlim.snrDb = repeatOut.snrDb;
repeatSlim.truthFdRefHz = repeatOut.truthFdRefHz;
repeatSlim.truthFdRateHzPerSec = repeatOut.truthFdRateHzPerSec;
repeatSlim.toothStepHz = repeatOut.toothStepHz;
repeatSlim.fdRangeOracle = repeatOut.fdRangeOracle;
repeatSlim.fdRateRangeOracle = repeatOut.fdRateRangeOracle;
repeatSlim.caseRowList = repeatOut.caseRowList;
repeatSlim.repeatTotalMs = repeatOut.repeatTotalMs;
repeatSlim.warningSeen = repeatOut.warningSeen;
repeatSlim.warningId = repeatOut.warningId;
repeatSlim.warningMessage = repeatOut.warningMessage;
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve the Doppler comb spacing from fixture time offsets.

timeOffsetSec = reshape(fixture.sceneSeq.timeOffsetSec, [], 1);
if numel(timeOffsetSec) < 2
  toothStepHz = NaN;
  return;
end
frameStepSec = median(diff(timeOffsetSec));
toothStepHz = 1 / frameStepSec;
end

function value = localResolveTruthScalar(truth, fieldNameList)
%LOCALRESOLVETRUTHSCALAR Read the first finite scalar truth field from a list.

value = NaN;
for iField = 1:numel(fieldNameList)
  fieldName = char(fieldNameList{iField});
  if isfield(truth, fieldName) && ~isempty(truth.(fieldName))
    cand = truth.(fieldName);
    cand = cand(1);
    if isfinite(cand)
      value = cand;
      return;
    end
  end
end
end

function txt = localFormatRow(x)
%LOCALFORMATROW Format a numeric vector as one compact comma-separated string.

txt = strjoin(compose('%.6g', x(:).'), ', ');
end

function lineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build compact HTML metric lines for notifications.

config = scanData.config;
aggregateTable = scanData.aggregateTable;
lineList = strings(0, 1);
lineList(end + 1, 1) = "• Tasks: <code>" + string(localGetFieldOrDefault(config, 'numTask', height(scanData.perfTable))) + "</code>";
lineList(end + 1, 1) = "• Grid: <code>P=" + string(localFormatRow(config.frameCountList)) + ", SNR=" + ...
  string(localFormatRow(config.snrDbList)) + " dB</code>";
cpU = localFindMetricRow(aggregateTable, "CP-U", max(config.snrDbList), max(config.frameCountList));
ipU = localFindMetricRow(aggregateTable, "IP-U", max(config.snrDbList), max(config.frameCountList));
if ~isempty(cpU)
  lineList(end + 1, 1) = "• CP-U angle RMSE: <code>" + ...
    string(sprintf('%.4g deg', cpU.angleRmseDeg)) + "</code>";
  lineList(end + 1, 1) = "• CP-U fdRef RMSE: <code>" + ...
    string(sprintf('%.4g Hz', cpU.fdRefRmseHz)) + "</code>";
  lineList(end + 1, 1) = "• CP-U non-ref coh floor: <code>" + ...
    string(sprintf('%.3f', cpU.nonRefCoherenceFloorMedian)) + "</code>";
end
if ~isempty(ipU)
  lineList(end + 1, 1) = "• IP-U / CP-U angle ratio: <code>" + ...
    string(sprintf('%.3f', ipU.angleRatioToCp)) + "</code>";
  lineList(end + 1, 1) = "• IP-U / CP-U fd ratio: <code>" + ...
    string(sprintf('%.3f', ipU.fdRatioToCp)) + "</code>";
end
end

function row = localFindMetricRow(aggregateTable, displayName, snrDb, numFrame)
%LOCALFINDMETRICROW Find one aggregate row for compact notification metrics.

row = [];
if isempty(aggregateTable) || ~istable(aggregateTable)
  return;
end
mask = aggregateTable.displayName == string(displayName) & ...
  aggregateTable.snrDb == snrDb & aggregateTable.numFrame == numFrame;
idx = find(mask, 1, 'first');
if ~isempty(idx)
  row = aggregateTable(idx, :);
end
end

function progressTracker = localCreateScanProgressTracker(totalCount)
%LOCALCREATESCANPROGRESSTRACKER Create a best-effort scan progress tracker.

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

function localAdvanceScanProgress(progressTracker, isParforWorker)
%LOCALADVANCESCANPROGRESS Advance progress from serial or parfor execution.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
if isParforWorker
  if isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
    send(progressTracker.queue, 1);
  end
else
  try
    progressbar('advance');
  catch
  end
end
end

function localCloseScanProgressTracker(progressTracker)
%LOCALCLOSESCANPROGRESSTRACKER Close a best-effort scan progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Return true when an outer scan loop can use parfor.

useParfor = false;
if numTask <= 1
  return;
end
try
  useParfor = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && ...
    license('test', 'Distrib_Computing_Toolbox');
catch
  useParfor = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct field with a fallback value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
