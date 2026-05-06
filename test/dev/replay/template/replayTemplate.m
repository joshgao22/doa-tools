% replayTemplate
% Copy this file to test/dev/replay/replayYourTopic.m before editing.
% The template standardizes replay sections, snapshot storage, compact
% summary output, and best-effort Telegram notification. It is not a shared
% replay framework and should not be run as a real experiment.

clear; close all; clc;

%% Replay configuration

replayName = "replayTemplate";
baseSeed = 253;
numRepeat = 40;
snrDb = 10;
saveSnapshot = true;
checkpointEnable = false;
notifyTelegramEnable = false;

config = struct();
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.snrDb = snrDb;
config.saveSnapshot = saveSnapshot;
config.checkpointEnable = checkpointEnable;
config.notifyTelegramEnable = notifyTelegramEnable;

runTic = tic;
replayData = struct();
runDir = '';

try
  %% Build context and flow options

  repoRoot = localFindRepoRoot();
  runKey = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
  config.runKey = runKey;

  if config.checkpointEnable
    runDir = fullfile(repoRoot, 'tmp', char(replayName), runKey);
  end

  printMfReplayHeader(char(replayName), config, runDir);

  % TODO: Build context, flow options, and replay-specific candidate sets here.
  % Example:
  % contextOpt = struct('baseSeed', config.baseSeed, 'snrDb', config.snrDb);
  % context = buildDynamicDualSatEciContext(contextOpt);

  %% Run replay batch

  repeatResultCell = cell(config.numRepeat, 1);
  for iRepeat = 1:config.numRepeat
    repeatSeed = config.baseSeed + iRepeat - 1;
    repeatResultCell{iRepeat} = localRunSingleRepeat(repeatSeed, config);
  end

  %% Data storage

  replayData = struct();
  replayData.replayName = replayName;
  replayData.config = config;
  replayData.repeatResultCell = repeatResultCell;
  replayData.summaryTable = localBuildSummaryTable(repeatResultCell);
  replayData.elapsedSec = toc(runTic);
  replayData = finalizeMfReplayResult(replayData, runDir);

  if config.saveSnapshot
    saveOpt = struct();
    saveOpt.includeVars = {'replayData'};
    saveOpt.extraMeta = struct('replayName', char(replayName));
    saveOpt.verbose = true;
    replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
  else
    replayData.snapshotFile = '';
  end

  %% Summary output and plotting

  printMfReplaySection('Template summary', replayData.summaryTable);

  metricLineList = localBuildTelegramMetricLines(replayData);
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "DONE", ...
    'config', config, ...
    'snapshotFile', replayData.snapshotFile, ...
    'elapsedSec', replayData.elapsedSec, ...
    'metricLineList', metricLineList, ...
    'commentLineList', "Template completed. Replace placeholder logic before use."));

catch ME
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'checkpointDir', runDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function result = localRunSingleRepeat(repeatSeed, config)
%LOCALRUNSINGLEREPEAT Placeholder for one replay repeat.
% Replace this helper with replay-specific orchestration. Do not put formal
% estimator logic here; call estimator / flow / common helpers instead.

result = struct();
result.seed = repeatSeed;
result.snrDb = config.snrDb;
result.angleErrDeg = NaN;
result.status = "template-placeholder";
end

function summaryTable = localBuildSummaryTable(repeatResultCell)
%LOCALBUILDSUMMARYTABLE Build a lightweight replay summary table.

numRepeat = numel(repeatResultCell);
seed = zeros(numRepeat, 1);
angleErrDeg = NaN(numRepeat, 1);
status = strings(numRepeat, 1);
for iRepeat = 1:numRepeat
  repeatResult = repeatResultCell{iRepeat};
  seed(iRepeat) = repeatResult.seed;
  angleErrDeg(iRepeat) = repeatResult.angleErrDeg;
  status(iRepeat) = repeatResult.status;
end
summaryTable = table(seed, angleErrDeg, status);
end

function metricLineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build replay-specific HTML-ready metric lines.
% Keep replay-specific metrics local. The common notification helper only
% wraps these preformatted lines in the standard Telegram shell.

metricLineList = strings(0, 1);
if isfield(replayData, 'summaryTable') && ~isempty(replayData.summaryTable)
  metricLineList(end + 1, 1) = "• repeats: <code>" + string(height(replayData.summaryTable)) + "</code>";
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
  error('replayTemplate:RepoRootNotFound', 'Cannot find repository root from template path.');
end
end
