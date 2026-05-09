% scanTemplate
% Copy this file to test/dev/scan/scanYourTopic.m before editing.
% The template standardizes scan sections, task-grid checkpointing, snapshot
% storage, compact summary output, and best-effort Telegram notification. It
% is not a shared scan framework and should not be used as a real result.

clear; close all; clc;

%% Scan configuration

scanName = "scanTemplate";
baseSeed = 253;
numRepeat = 3;
snrDbList = [0; 10];
saveSnapshot = false;
checkpointEnable = false;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = false;

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
numRepeat = numel(seedList);

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.baseSeed = baseSeed;
scanConfig.numRepeat = numRepeat;
scanConfig.seedList = seedList;
scanConfig.snrDbList = snrDbList;
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.checkpointResume = logical(checkpointResume);
scanConfig.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
repoRoot = localFindRepoRoot();
checkpointDir = "";
checkpointRunState = struct();
checkpointCleanupReport = struct([]);
scanData = struct();

try
  %% Build context and scan tasks

  taskList = localBuildTaskList(seedList, snrDbList);
  useParfor = localCanUseParfor(numel(taskList));
  if scanConfig.checkpointEnable
    checkpointOpt = localBuildCheckpointOpt(repoRoot, scanName, scanConfig, taskList, useParfor);
    checkpointDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
  else
    checkpointOpt = struct();
  end
  printMfScanHeader(char(scanName), scanConfig, checkpointDir);

  % TODO: Build scene, truth, model options, fixture options, or scan-specific
  % static context here. Keep heavy raw data out of scanData and checkpoint files.
  sharedData = struct();
  sharedData.scanName = scanName;

  %% Run scan batch

  if scanConfig.checkpointEnable
    checkpointRunState = runPerfTaskGridWithCheckpoint(taskList, sharedData, ...
      @localRunCheckpointTask, checkpointOpt);
    taskOutCell = checkpointRunState.resultCell;
  else
    taskOutCell = cell(numel(taskList), 1);
    if useParfor
      parfor iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), sharedData);
      end
    else
      for iTask = 1:numel(taskList)
        taskOutCell{iTask} = localRunOneTask(taskList(iTask), sharedData);
      end
    end
  end

  %% Data storage

  scanTable = localBuildScanTable(taskOutCell);
  aggregateTable = localBuildAggregateTable(scanTable);
  if scanConfig.checkpointEnable && localGetFieldOrDefault(checkpointRunState, 'isComplete', false) && ...
      scanConfig.checkpointCleanupOnSuccess
    checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(checkpointRunState, ...
      struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
  end
  checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunState, ...
    checkpointDir, checkpointCleanupReport);

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.scanTable = scanTable;
  scanData.aggregateTable = aggregateTable;
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupReport = checkpointCleanupReport;
  scanData.checkpointDir = checkpointDir;
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

  if ~exist('scanData', 'var') || ~isstruct(scanData)
    error('scanTemplate:MissingScanData', ...
      'Scan data is missing. Load a snapshot containing scanData before running this section.');
  end
  scanData = localValidateScanDataForSummary(scanData);
  scanNameForReport = string(localGetFieldOrDefault(scanData, 'scanName', "scanTemplate"));
  scanConfigForReport = localGetFieldOrDefault(scanData, 'config', struct());
  scanSnapshotFile = localGetFieldOrDefault(scanData, 'snapshotFile', "");
  scanCheckpointDir = string(localGetFieldOrDefault(scanData, 'checkpointDir', ""));
  scanElapsedSec = localGetFieldOrDefault(scanData, 'elapsedSec', NaN);

  printMfScanSection('Template aggregate summary', scanData.aggregateTable);
  if isfield(scanData, 'checkpointSummaryTable') && istable(scanData.checkpointSummaryTable) && ...
      height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  scanData.plotData = localPlotScan(scanData);

  metricLineList = localBuildTelegramMetricLines(scanData);
  notifyMfScanStatus(struct( ...
    'scanName', scanNameForReport, ...
    'statusText', "DONE", ...
    'config', scanConfigForReport, ...
    'snapshotFile', scanSnapshotFile, ...
    'checkpointDir', scanCheckpointDir, ...
    'elapsedSec', scanElapsedSec, ...
    'metricLineList', metricLineList, ...
    'commentLineList', [ ...
      "Template completed with placeholder task output."; ...
      "Replace local task, aggregate, and plot logic before using as a real scan."]));

catch ME
  localPrintCheckpointFailureHint(checkpointDir);
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


function scanData = localValidateScanDataForSummary(scanData)
%LOCALVALIDATESCANDATAFORSUMMARY Validate scanData contents for summary reruns.

if ~isfield(scanData, 'aggregateTable') || ~istable(scanData.aggregateTable)
  error('scanTemplate:MissingAggregateTable', ...
    'scanData.aggregateTable is missing. Store summary and plot inputs inside scanData before saving.');
end
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  scanData.plotData = localBuildPlotData(scanData.aggregateTable);
end
end

function taskOut = localRunCheckpointTask(task, sharedData)
%LOCALRUNCHECKPOINTTASK Run one checkpointed scan task.

taskOut = localRunOneTask(task, sharedData);
end

function taskOut = localRunOneTask(task, sharedData)
%LOCALRUNONETASK Placeholder for one independent scan task.
% Replace this helper with scan-specific orchestration. Do not put formal
% estimator logic here; call estimator / performance / flow helpers instead.

taskOut = struct();
taskOut.taskId = task.taskId;
taskOut.taskSeed = task.taskSeed;
taskOut.snrDb = task.snrDb;
taskOut.metricValue = NaN;
taskOut.status = "template-placeholder";
taskOut.scanName = string(sharedData.scanName);
end

function taskList = localBuildTaskList(seedList, snrDbList)
%LOCALBUILDTASKLIST Build one flat task list from explicit scan dimensions.

taskList = repmat(struct('taskId', 0, 'taskSeed', 0, 'snrDb', NaN), 0, 1);
taskId = 0;
for iSeed = 1:numel(seedList)
  for iSnr = 1:numel(snrDbList)
    taskId = taskId + 1;
    taskList(taskId, 1).taskId = taskId; %#ok<AGROW>
    taskList(taskId, 1).taskSeed = seedList(iSeed); %#ok<AGROW>
    taskList(taskId, 1).snrDb = snrDbList(iSnr); %#ok<AGROW>
  end
end
end

function checkpointOpt = localBuildCheckpointOpt(repoRoot, scanName, scanConfig, taskList, useParfor)
%LOCALBUILDCHECKPOINTOPT Build task-level checkpoint runner options.

checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, taskList);
checkpointOpt.outputRoot = fullfile(repoRoot, 'tmp');
checkpointOpt.useParfor = logical(useParfor);
checkpointOpt.resume = logical(scanConfig.checkpointResume);
checkpointOpt.meta = localBuildCheckpointMeta(scanConfig, taskList);
checkpointOpt.progressFcn = [];
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end

function runKey = localBuildCheckpointRunKey(scanConfig, taskList)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable run key from task-defining config.

seedList = unique([taskList.taskSeed], 'stable');
snrList = unique([taskList.snrDb], 'stable');
runKey = sprintf('seed%dto%d_rep%d_snr%s', seedList(1), seedList(end), ...
  scanConfig.numRepeat, localJoinNumericToken(snrList));
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
meta.seedList = reshape(unique([taskList.taskSeed], 'stable'), 1, []);
meta.snrDbList = reshape(unique([taskList.snrDb], 'stable'), 1, []);
meta.numTask = numel(taskList);
end

function scanTable = localBuildScanTable(taskOutCell)
%LOCALBUILDSCANTABLE Collect lightweight task output into one scan table.

numTask = numel(taskOutCell);
taskId = zeros(numTask, 1);
taskSeed = zeros(numTask, 1);
snrDb = NaN(numTask, 1);
metricValue = NaN(numTask, 1);
status = strings(numTask, 1);
for iTask = 1:numTask
  taskOut = taskOutCell{iTask};
  taskId(iTask) = taskOut.taskId;
  taskSeed(iTask) = taskOut.taskSeed;
  snrDb(iTask) = taskOut.snrDb;
  metricValue(iTask) = taskOut.metricValue;
  status(iTask) = taskOut.status;
end
scanTable = table(taskId, taskSeed, snrDb, metricValue, status);
end

function aggregateTable = localBuildAggregateTable(scanTable)
%LOCALBUILDAGGREGATETABLE Build a compact aggregate table for command output.

if isempty(scanTable)
  aggregateTable = table();
  return;
end
snrList = unique(scanTable.snrDb, 'stable');
numRow = numel(snrList);
snrDb = reshape(snrList, [], 1);
numTask = zeros(numRow, 1);
numFiniteMetric = zeros(numRow, 1);
metricMean = NaN(numRow, 1);
for iRow = 1:numRow
  mask = scanTable.snrDb == snrDb(iRow);
  metricVec = scanTable.metricValue(mask);
  numTask(iRow) = nnz(mask);
  numFiniteMetric(iRow) = nnz(isfinite(metricVec));
  if any(isfinite(metricVec))
    metricMean(iRow) = mean(metricVec(isfinite(metricVec)));
  end
end
aggregateTable = table(snrDb, numTask, numFiniteMetric, metricMean);
end

function plotData = localBuildPlotData(aggregateTable)
%LOCALBUILDPLOTDATA Keep only lightweight arrays needed to redraw figures.

plotData = struct();
if isempty(aggregateTable)
  plotData.snrDb = [];
  plotData.metricMean = [];
  return;
end
plotData.snrDb = aggregateTable.snrDb;
plotData.metricMean = aggregateTable.metricMean;
end

function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw a placeholder scan figure from lightweight plot data.
% Replace this block with scan-specific figures. Do not save image files.

plotData = scanData.plotData;
if isempty(plotData.snrDb)
  return;
end
figure('Name', char(scanData.scanName));
plot(plotData.snrDb, plotData.metricMean, '-o', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Template metric');
title(char(scanData.scanName), 'Interpreter', 'none');
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runState, checkpointDir, cleanupReport)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Summarize checkpoint status for scanData.

if strlength(string(checkpointDir)) == 0 && ~isfield(runState, 'numTask')
  checkpointSummaryTable = table();
  return;
end
runName = string(localGetFieldOrDefault(runState, 'runName', ""));
runKey = string(localGetFieldOrDefault(runState, 'runKey', ""));
runDir = string(localGetFieldOrDefault(runState, 'runDir', checkpointDir));
numTask = double(localGetFieldOrDefault(runState, 'numTask', 0));
numDone = double(localGetFieldOrDefault(runState, 'numDone', 0));
isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
cleanupItemCount = numel(cleanupReport);
checkpointSummaryTable = table(runName, runKey, runDir, numTask, numDone, ...
  isComplete, cleanupItemCount);
end

function tf = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Return true when an outer scan loop can use parfor.

tf = false;
if numTask <= 1
  return;
end
try
  tf = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && ...
    license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print retained checkpoint location on failure.

checkpointDir = string(checkpointDir);
if strlength(checkpointDir) == 0
  return;
end
if isfolder(checkpointDir)
  fprintf('Checkpoint directory retained for resume: %s\n', char(checkpointDir));
end
end

function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build scan-specific HTML-ready metric lines.

metricLineList = strings(0, 1);
if isfield(scanData, 'aggregateTable') && istable(scanData.aggregateTable) && ...
    height(scanData.aggregateTable) > 0
  metricLineList(end + 1, 1) = "• tasks=<code>" + ...
    string(sum(scanData.aggregateTable.numTask)) + "</code>, SNR points=<code>" + ...
    string(height(scanData.aggregateTable)) + "</code>";
end
if isfield(scanData, 'checkpointSummaryTable') && istable(scanData.checkpointSummaryTable) && ...
    height(scanData.checkpointSummaryTable) > 0
  metricLineList(end + 1, 1) = "• checkpoint done=<code>" + ...
    string(scanData.checkpointSummaryTable.numDone(1)) + "/" + ...
    string(scanData.checkpointSummaryTable.numTask(1)) + "</code>";
end
metricLineList(end + 1, 1) = "• report style: local metrics, lightweight scanData snapshot";
end

function token = localJoinNumericToken(value)
%LOCALJOINNUMERICTOKEN Join numeric values into a filesystem-friendly token.

value = reshape(double(value), 1, []);
if isempty(value)
  token = 'empty';
  return;
end
token = strjoin(compose('%g', value), '_');
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function repoRoot = localFindRepoRoot()
%LOCALFINDREPOROOT Find repository root from this template location.

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
repoRoot = thisDir;
while ~isempty(repoRoot) && ~isfile(fullfile(repoRoot, 'AGENTS.md'))
  parentDir = fileparts(repoRoot);
  if strcmp(parentDir, repoRoot)
    break;
  end
  repoRoot = parentDir;
end
if isempty(repoRoot) || ~isfile(fullfile(repoRoot, 'AGENTS.md'))
  error('scanTemplate:RepoRootNotFound', 'Cannot find repository root from template path.');
end
end
