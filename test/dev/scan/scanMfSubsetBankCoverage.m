% scanMfSubsetBankCoverage
% Purpose: compare curated, rescue, and random subset banks for dynamic
% tooth selection coverage, no-truth consensus diagnostics, runtime cost,
% and curated12-easy damage.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves scanData in the workspace; saveSnapshot=true saves only
% scanData via saveExpSnapshot.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfSubsetBankCoverage";
saveSnapshot = true;
notifyTelegramEnable = true;
optVerbose = false;
strategyPreset = "scheduleComboQuick";
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
baseSeed = 253;
snrDb = 10;
numRepeat = 8;
toothResidualTolHz = 50;
nearToothIdxTol = 1;
baselineStrategyName = "curated12";
easyAngleTolDeg = 0.005;
damageAngleTolDeg = 1e-3;
randomTrialList = [1, 2, 4];
toothHistogramBinCount = 31;
toothConsensusResidualPenaltyHz = 75;

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

config = struct();
config.scanName = string(scanName);
config.saveSnapshot = logical(saveSnapshot);
config.notifyTelegramEnable = logical(notifyTelegramEnable);
config.optVerbose = logical(optVerbose);
config.strategyPreset = string(strategyPreset);
config.checkpointEnable = logical(checkpointEnable);
config.checkpointResume = logical(checkpointResume);
config.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
config.baseSeed = baseSeed;
config.snrDb = snrDb;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.toothResidualTolHz = toothResidualTolHz;
config.nearToothIdxTol = nearToothIdxTol;
config.baselineStrategyName = string(baselineStrategyName);
config.easyAngleTolDeg = easyAngleTolDeg;
config.damageAngleTolDeg = damageAngleTolDeg;
config.randomTrialList = randomTrialList;
config.toothHistogramBinCount = toothHistogramBinCount;
config.toothConsensusResidualPenaltyHz = toothConsensusResidualPenaltyHz;

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
config.runKey = string(runKey);
checkpointDir = "";
checkpointDirCell = strings(0, 1);
scanData = struct();

try
  %% Build context and scan tasks
  strategyList = localBuildStrategyList(config);
  numStrategy = numel(strategyList);
  config.numStrategy = numStrategy;
  config.numTask = numStrategy * config.numRepeat;

  tableCell = cell(numStrategy, 1);
  candidateCell = cell(numStrategy, 1);
  contextCell = cell(numStrategy, 1);
  flowOptCell = cell(numStrategy, 1);
  repeatCellByStrategy = cell(numStrategy, 1);
  runStateCell = cell(numStrategy, 1);
  batchOptCell = cell(numStrategy, 1);
  checkpointDirCell = strings(numStrategy, 1);
  for iStrategy = 1:numStrategy
    batchOptCell{iStrategy} = localBuildBatchOpt(config, strategyList(iStrategy));
    checkpointDirCell(iStrategy) = localGetBatchCheckpointDir(batchOptCell{iStrategy});
  end
  checkpointDir = localBuildCheckpointDirSummary(checkpointDirCell);

  printMfScanHeader(char(scanName), config, checkpointDir);
  localPrintScanExtraConfig(config);

  %% Run scan batch
  localLog(sprintf('Run %d subset-bank strategies.', numStrategy));
  for iStrategy = 1:numStrategy
    strategy = strategyList(iStrategy);
    localLog(sprintf('Run strategy %d/%d: %s.', iStrategy, numStrategy, strategy.strategyName));
    batchOpt = batchOptCell{iStrategy};
    [t, repeatCell, contextCell{iStrategy}, flowOptCell{iStrategy}, runStateCell{iStrategy}] = runSimpleDynamicFlowReplayBatch(batchOpt);
    t.strategyName = repmat(strategy.strategyName, height(t), 1);
    t.bankSize = repmat(strategy.bankSize, height(t), 1);
    t.numRandomTrial = repmat(strategy.numRandomTrial, height(t), 1);
    t = localAnnotateToothMetrics(t, config);
    tableCell{iStrategy} = t;
    candidateCell{iStrategy} = localBuildCandidateTable(repeatCell, strategy.strategyName);
    repeatCellByStrategy{iStrategy} = repeatCell;
  end

  config = localAttachFrameTimingConfig(config, contextCell);
  scanTable = vertcat(tableCell{:});
  scanTable = localAnnotateBaselineComparison(scanTable, config);
  candidateTable = localAnnotateCandidateToothMetrics(localVertcatTableCell(candidateCell), config);
  candidateTable = localAnnotateCandidateAnchorGroups(candidateTable, repeatCellByStrategy, config);
  candidateSeedCoverageTable = localBuildCandidateSeedCoverageTable(scanTable, candidateTable);
  consensusTable = localBuildConsensusDiagnosticTable(scanTable, candidateTable, config);
  consensusAggregateTable = localBuildConsensusAggregateTable(consensusTable);
  scheduleFeatureTable = localBuildScheduleFeatureTable(strategyList);
  transitionTable = localBuildBaselineTransitionTable(scanTable, config);
  aggregateTable = localBuildAggregateTable(scanTable, candidateTable, candidateSeedCoverageTable, config);
  toothHistogramTable = localBuildToothHistogramTable(scanTable, config);
  representative = localSelectRepresentative(scanTable, candidateTable, repeatCellByStrategy);
  checkpointCleanupTable = localCleanupCheckpointRunStates(runStateCell, config);
  checkpointSummaryTable = localBuildCheckpointSummaryTable(strategyList, runStateCell, checkpointDirCell, checkpointCleanupTable);

  %% Data storage
  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = config;
  scanData.strategyList = strategyList;
  scanData.contextSummary = localBuildContextSummary(contextCell, flowOptCell, strategyList);
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupTable = checkpointCleanupTable;
  scanData.scanTable = scanTable;
  scanData.candidateTable = candidateTable;
  scanData.candidateSeedCoverageTable = candidateSeedCoverageTable;
  scanData.consensusTable = consensusTable;
  scanData.consensusAggregateTable = consensusAggregateTable;
  scanData.scheduleFeatureTable = scheduleFeatureTable;
  scanData.transitionTable = transitionTable;
  scanData.aggregateTable = aggregateTable;
  scanData.toothHistogramTable = toothHistogramTable;
  scanData.representative = representative;
  scanData.plotData = localBuildPlotData(scanTable, aggregateTable, toothHistogramTable, candidateSeedCoverageTable, consensusAggregateTable);
  scanData.repeatCellByStrategy = localStripRepeatCellByStrategy(repeatCellByStrategy);
  scanData.elapsedSec = toc(runTic);
  if ~config.checkpointEnable
    scanData = finalizeMfScanResult(scanData, "");
  end

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output and plotting
  config = localEnsureConsensusConfig(scanData.config, scanData);
  scanData.config = config;
  scanTable = localAnnotateToothMetrics(scanData.scanTable, config);
  scanData.scanTable = localAnnotateBaselineComparison(scanTable, config);
  scanTable = scanData.scanTable;
  if ~isfield(scanData, 'candidateTable')
    scanData.candidateTable = table();
  end
  scanData.candidateTable = localAnnotateCandidateToothMetrics(scanData.candidateTable, config);
  if ~all(ismember({ 'anchorFdRefHz', 'anchorToothGroupIdx' }, scanData.candidateTable.Properties.VariableNames))
    scanData.candidateTable = localAnnotateCandidateAnchorGroups(scanData.candidateTable, {}, config);
  end
  if ~isfield(scanData, 'candidateSeedCoverageTable')
    scanData.candidateSeedCoverageTable = localBuildCandidateSeedCoverageTable(scanTable, scanData.candidateTable);
  end
  if ~isfield(scanData, 'consensusTable')
    scanData.consensusTable = localBuildConsensusDiagnosticTable(scanTable, scanData.candidateTable, config);
  end
  scanData.consensusAggregateTable = localBuildConsensusAggregateTable(scanData.consensusTable);
  if ~isfield(scanData, 'transitionTable')
    scanData.transitionTable = localBuildBaselineTransitionTable(scanTable, config);
  end
  aggregateTable = localBuildAggregateTable(scanTable, scanData.candidateTable, scanData.candidateSeedCoverageTable, config);
  consensusAggregateTable = scanData.consensusAggregateTable;
  scanData.aggregateTable = aggregateTable;
  if ~isfield(scanData, 'scheduleFeatureTable') || ~istable(scanData.scheduleFeatureTable)
    scanData.scheduleFeatureTable = localBuildScheduleFeatureTable(scanData.strategyList);
  end
  if isfield(scanData, 'toothHistogramTable')
    toothHistogramTable = scanData.toothHistogramTable;
  else
    toothHistogramTable = localBuildToothHistogramTable(scanTable, config);
    scanData.toothHistogramTable = toothHistogramTable;
  end

  fprintf('  tooth histogram rows            : %d (saved in scanData.toothHistogramTable)\n', height(toothHistogramTable));
  fprintf('  consensus diagnostic            : anchor-group soft residual penalty scale %.2f Hz (no-truth diagnostic)\n', ...
    config.toothConsensusResidualPenaltyHz);
  if isfield(scanData, 'checkpointSummaryTable') && ~isempty(scanData.checkpointSummaryTable)
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  printMfScanSection('Subset-bank aggregate', aggregateTable);
  localPrintTablePreview('Schedule feature table', localGetScheduleFeatureTableForPrint(scanData), 24);
  localPrintTablePreview('No-truth consensus diagnostic aggregate', consensusAggregateTable, 16);
  localPrintTablePreview('Strategy repeat preview', localSelectDisplayTable(scanTable), 16);
  localPrintTablePreview('Baseline transition preview', scanData.transitionTable, 16);
  localPrintTablePreview('Candidate seed coverage preview', scanData.candidateSeedCoverageTable, 16);
  localPrintTablePreview('No-truth consensus diagnostic preview', scanData.consensusTable, 16);
  localPrintTablePreview('Candidate preview', scanData.candidateTable, 16);

  scanData.plotData = localPlotScan(scanTable, aggregateTable, toothHistogramTable, ...
    scanData.candidateSeedCoverageTable, consensusAggregateTable);

  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "DONE", ...
    'config', scanData.config, ...
    'snapshotFile', scanData.snapshotFile, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Subset-bank coverage scan completed for offline schedule evaluation."; ...
      "Truth metrics are evaluation-only and must not enter runtime selector or adoption."]));

catch ME
  localPrintCheckpointFailureHint(checkpointDirCell);
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function localPrintScanExtraConfig(scanConfig)
%LOCALPRINTSCANEXTRACONFIG Print subset-bank-specific settings.
fprintf('  %-32s : %d\n', 'strategy count', scanConfig.numStrategy);
fprintf('  %-32s : %s\n', 'strategy preset', char(scanConfig.strategyPreset));
fprintf('  %-32s : %s\n', 'baseline strategy', char(scanConfig.baselineStrategyName));
fprintf('  %-32s : %s\n', 'random trial list', mat2str(scanConfig.randomTrialList));
fprintf('  %-32s : |toothIdx| <= %d, residual <= %.2f Hz\n', 'near-tooth tolerance', ...
  scanConfig.nearToothIdxTol, scanConfig.toothResidualTolHz);
fprintf('  %-32s : baseline angle <= %.4f deg, damage > %.4f deg or lost truth tooth\n', ...
  'easy damage rule', scanConfig.easyAngleTolDeg, scanConfig.damageAngleTolDeg);
end

function checkpointDir = localBuildCheckpointDirSummary(checkpointDirCell)
%LOCALBUILDCHECKPOINTDIRSUMMARY Return a compact checkpoint root for logs.
checkpointDirCell = string(checkpointDirCell(:));
checkpointDirCell = checkpointDirCell(strlength(checkpointDirCell) > 0);
if isempty(checkpointDirCell)
  checkpointDir = "";
  return;
end
parentDirList = strings(numel(checkpointDirCell), 1);
for iDir = 1:numel(checkpointDirCell)
  parentDirList(iDir) = string(fileparts(char(checkpointDirCell(iDir))));
end
if numel(unique(parentDirList)) == 1
  checkpointDir = parentDirList(1);
else
  checkpointDir = strjoin(checkpointDirCell, newline);
end
end

function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build scan-specific compact HTML metric lines.
metricLineList = strings(0, 1);
if ~isstruct(scanData) || ~isfield(scanData, 'aggregateTable') || ~istable(scanData.aggregateTable)
  return;
end
config = localGetFieldOrDefault(scanData, 'config', struct());
aggregateTable = scanData.aggregateTable;
metricLineList(end + 1, 1) = sprintf('• preset: <code>%s</code>, strategies: <code>%d</code>', ...
  char(localHtmlEscape(localGetFieldOrDefault(config, 'strategyPreset', ""))), height(aggregateTable));

baselineName = string(localGetFieldOrDefault(config, 'baselineStrategyName', ""));
baseRow = localFindStrategyRow(aggregateTable, baselineName);
if ~isempty(baseRow)
  metricLineList(end + 1, 1) = sprintf('• baseline truth/near hit: <code>%.3f / %.3f</code>', ...
    localTableScalar(baseRow, 'truthToothHitRate', NaN), localTableScalar(baseRow, 'nearToothHitRate', NaN));
end

[bestTruthHit, bestTruthIdx] = localBestTableValue(aggregateTable, 'truthToothHitRate', 'max');
if isfinite(bestTruthHit) && bestTruthIdx > 0
  bestName = aggregateTable.strategyName(bestTruthIdx);
  metricLineList(end + 1, 1) = sprintf('• best selected truth hit: <code>%s %.3f</code>', ...
    char(localHtmlEscape(bestName)), bestTruthHit);
end

maxEasyDamage = localTableColumnMax(aggregateTable, 'easyDamageRate');
maxSelectedMiss = localTableColumnMax(aggregateTable, 'selectedMissDespiteTruthCandidateRate');
metricLineList(end + 1, 1) = sprintf('• max easy damage / selected-miss: <code>%.3f / %.3f</code>', ...
  maxEasyDamage, maxSelectedMiss);
end

function row = localFindStrategyRow(t, strategyName)
%LOCALFINDSTRATEGYROW Return one aggregate row by strategy name.
row = table();
if isempty(t) || ~istable(t) || ~any(strcmp(t.Properties.VariableNames, 'strategyName'))
  return;
end
idx = find(t.strategyName == string(strategyName), 1, 'first');
if ~isempty(idx)
  row = t(idx, :);
end
end

function value = localTableScalar(t, varName, defaultValue)
%LOCALTABLESCALAR Read a scalar table value from the first row.
value = defaultValue;
if istable(t) && height(t) >= 1 && any(strcmp(t.Properties.VariableNames, varName)) && ~isempty(t.(varName))
  value = t.(varName)(1);
end
end

function [value, idx] = localBestTableValue(t, varName, modeText)
%LOCALBESTTABLEVALUE Return the best finite value in a numeric table column.
value = NaN;
idx = 0;
if isempty(t) || ~istable(t) || ~any(strcmp(t.Properties.VariableNames, varName))
  return;
end
x = double(t.(varName));
validIdx = find(isfinite(x));
if isempty(validIdx)
  return;
end
if strcmpi(modeText, 'min')
  [value, relIdx] = min(x(validIdx));
else
  [value, relIdx] = max(x(validIdx));
end
idx = validIdx(relIdx);
end

function value = localTableColumnMax(t, varName)
%LOCALTABLECOLUMNMAX Return max of a numeric table column with NaN fallback.
value = NaN;
if isempty(t) || ~istable(t) || ~any(strcmp(t.Properties.VariableNames, varName))
  return;
end
value = max(double(t.(varName)), [], 'omitnan');
end

function textValue = localHtmlEscape(textValue)
%LOCALHTMLESCAPE Escape text for scan-specific Telegram metric lines.
textValue = string(textValue);
textValue = replace(textValue, "&", "&amp;");
textValue = replace(textValue, "<", "&lt;");
textValue = replace(textValue, ">", "&gt;");
end

function strategyList = localBuildStrategyList(scanConfig)
%LOCALBUILDSTRATEGYLIST Build deterministic and random subset-bank variants.
[primaryOffsetCell, primaryLabelList, rescueOffsetCell, rescueLabelList] = getDynamicCuratedSubsetBank();
strategyList = repmat(localEmptyStrategy(), 0, 1);
strategyPreset = lower(string(scanConfig.strategyPreset));
firstRandomTrial = localFirstRandomTrial(scanConfig.randomTrialList);
[structuredOffsetCell, structuredLabelList] = localBuildStructuredScheduleBank(primaryOffsetCell, rescueOffsetCell);

if strategyPreset == "schedulecomboquick"
  strategyList(end + 1) = localMakeStrategy("curated12", primaryOffsetCell, primaryLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated124", [primaryOffsetCell, rescueOffsetCell(2)], ...
    [primaryLabelList; rescueLabelList(2)], 0);
  strategyList(end + 1) = localMakeStrategy("structuredCombo", ...
    structuredOffsetCell, structuredLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated12_structured", ...
    [primaryOffsetCell, structuredOffsetCell], [primaryLabelList; structuredLabelList], 0);
  return;
end

if strategyPreset == "schedulecomboconfirm"
  strategyList(end + 1) = localMakeStrategy("curated12", primaryOffsetCell, primaryLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated124", [primaryOffsetCell, rescueOffsetCell(2)], ...
    [primaryLabelList; rescueLabelList(2)], 0);
  strategyList(end + 1) = localMakeStrategy("fullRescue", ...
    [primaryOffsetCell, rescueOffsetCell], [primaryLabelList; rescueLabelList], firstRandomTrial);
  strategyList(end + 1) = localMakeStrategy("structuredCombo", ...
    structuredOffsetCell, structuredLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated12_structured", ...
    [primaryOffsetCell, structuredOffsetCell], [primaryLabelList; structuredLabelList], 0);
  strategyList(end + 1) = localMakeStrategy("curated124_structured", ...
    [primaryOffsetCell, rescueOffsetCell(2), structuredOffsetCell], ...
    [primaryLabelList; rescueLabelList(2); structuredLabelList], 0);
  return;
end

if strategyPreset == "scheduledesign"
  strategyList(end + 1) = localMakeStrategy("curated12", primaryOffsetCell, primaryLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated124", [primaryOffsetCell, rescueOffsetCell(2)], ...
    [primaryLabelList; rescueLabelList(2)], 0);
  strategyList(end + 1) = localMakeStrategy("curated12_random" + string(firstRandomTrial), ...
    primaryOffsetCell, primaryLabelList, firstRandomTrial);
  strategyList(end + 1) = localMakeStrategy("fullRescue", ...
    [primaryOffsetCell, rescueOffsetCell], [primaryLabelList; rescueLabelList], firstRandomTrial);
  for iSchedule = 1:numel(structuredOffsetCell)
    strategyList(end + 1) = localMakeStrategy(structuredLabelList(iSchedule), ...
      structuredOffsetCell(iSchedule), structuredLabelList(iSchedule), 0); %#ok<AGROW>
  end
  return;
end

if strategyPreset == "diagnosticcore"
  strategyList(end + 1) = localMakeStrategy("curated12", primaryOffsetCell, primaryLabelList, 0);
  strategyList(end + 1) = localMakeStrategy("curated124", [primaryOffsetCell, rescueOffsetCell(2)], ...
    [primaryLabelList; rescueLabelList(2)], 0);
  strategyList(end + 1) = localMakeStrategy("curated12_random" + string(firstRandomTrial), ...
    primaryOffsetCell, primaryLabelList, firstRandomTrial);
  strategyList(end + 1) = localMakeStrategy("fullRescue", ...
    [primaryOffsetCell, rescueOffsetCell], [primaryLabelList; rescueLabelList], firstRandomTrial);
  return;
end

strategyList(end + 1) = localMakeStrategy("curated12", primaryOffsetCell, primaryLabelList, 0);
strategyList(end + 1) = localMakeStrategy("curated123", [primaryOffsetCell, rescueOffsetCell(1)], ...
  [primaryLabelList; rescueLabelList(1)], 0);
strategyList(end + 1) = localMakeStrategy("curated124", [primaryOffsetCell, rescueOffsetCell(2)], ...
  [primaryLabelList; rescueLabelList(2)], 0);
strategyList(end + 1) = localMakeStrategy("curated1234", [primaryOffsetCell, rescueOffsetCell], ...
  [primaryLabelList; rescueLabelList], 0);

if strategyPreset == "cheapscreen"
  strategyList(end + 1) = localMakeStrategy("curated12_random" + string(firstRandomTrial), ...
    primaryOffsetCell, primaryLabelList, firstRandomTrial);
  strategyList(end + 1) = localMakeStrategy("fullRescue", ...
    [primaryOffsetCell, rescueOffsetCell], [primaryLabelList; rescueLabelList], firstRandomTrial);
  return;
end

if strategyPreset ~= "full"
  error('scanMfSubsetBankCoverage:UnknownStrategyPreset', ...
    'Unknown strategyPreset: %s. Use "scheduleComboQuick", "scheduleComboConfirm", "scheduleDesign", "diagnosticCore", "cheapScreen" or "full".', ...
    char(scanConfig.strategyPreset));
end

for iRandom = 1:numel(scanConfig.randomTrialList)
  numRandom = scanConfig.randomTrialList(iRandom);
  strategyList(end + 1) = localMakeStrategy("curated12_random" + string(numRandom), ...
    primaryOffsetCell, primaryLabelList, numRandom); %#ok<AGROW>
end
strategyList(end + 1) = localMakeStrategy("curated123_random1", ...
  [primaryOffsetCell, rescueOffsetCell(1)], [primaryLabelList; rescueLabelList(1)], 1);
strategyList(end + 1) = localMakeStrategy("curated124_random1", ...
  [primaryOffsetCell, rescueOffsetCell(2)], [primaryLabelList; rescueLabelList(2)], 1);
strategyList(end + 1) = localMakeStrategy("fullRescue", ...
  [primaryOffsetCell, rescueOffsetCell], [primaryLabelList; rescueLabelList], firstRandomTrial);
end
function numRandom = localFirstRandomTrial(randomTrialList)
%LOCALFIRSTRANDOMTRIAL Pick the smallest random-trial count for cheap screening.
randomTrialList = reshape(double(randomTrialList), [], 1);
randomTrialList = randomTrialList(isfinite(randomTrialList) & randomTrialList > 0);
if isempty(randomTrialList)
  numRandom = 1;
else
  numRandom = min(randomTrialList);
end
end



function [structuredOffsetCell, structuredLabelList] = localBuildStructuredScheduleBank(primaryOffsetCell, rescueOffsetCell)
%LOCALBUILDSTRUCTUREDSCHEDULEBANK Build deterministic nonuniform schedules for tooth screening.
offsetUniverse = localCollectOffsetUniverse([primaryOffsetCell, rescueOffsetCell]);
structuredLabelList = ["staggeredA"; "sparseRulerA"; "coprime34A"];
rawScheduleCell = { ...
  [-9, -8, -6, -3, 0, 1, 3, 6, 8, 10], ...
  [-9, -8, -5, -1, 0, 1, 4, 7, 9, 10], ...
  [-9, -8, -6, -4, -3, 0, 3, 4, 8, 9]};
structuredOffsetCell = cell(1, numel(rawScheduleCell));
for iSchedule = 1:numel(rawScheduleCell)
  structuredOffsetCell{iSchedule} = localNormalizeScheduleOffsets(rawScheduleCell{iSchedule}, offsetUniverse, numel(primaryOffsetCell{1}));
end
end

function offsetUniverse = localCollectOffsetUniverse(offsetCell)
%LOCALCOLLECTOFFSETUNIVERSE Collect the available legacy offset support.
valueList = [];
for iCell = 1:numel(offsetCell)
  valueList = [valueList, reshape(double(offsetCell{iCell}), 1, [])]; %#ok<AGROW>
end
valueList = unique(valueList(isfinite(valueList)));
if isempty(valueList)
  offsetUniverse = -9:10;
else
  offsetUniverse = min(valueList):max(valueList);
end
end

function schedule = localNormalizeScheduleOffsets(rawSchedule, offsetUniverse, targetCount)
%LOCALNORMALIZESCHEDULEOFFSETS Keep a structured schedule inside the fixture support.
schedule = unique(round(double(rawSchedule(:))).', 'stable');
schedule = schedule(ismember(schedule, offsetUniverse));
if ~ismember(0, schedule) && ismember(0, offsetUniverse)
  schedule = [schedule, 0];
end
candidate = reshape(double(offsetUniverse), 1, []);
[~, order] = sort(abs(candidate), 'ascend');
for iCand = 1:numel(order)
  if numel(schedule) >= targetCount
    break;
  end
  value = candidate(order(iCand));
  if ~ismember(value, schedule)
    schedule(end + 1) = value; %#ok<AGROW>
  end
end
schedule = unique(schedule, 'stable');
if numel(schedule) > targetCount
  keepMask = false(1, numel(schedule));
  keepMask(1:min(numel(schedule), targetCount)) = true;
  if ismember(0, schedule)
    keepMask(schedule == 0) = true;
  end
  if nnz(keepMask) > targetCount
    idx = find(keepMask & schedule ~= 0, 1, 'last');
    keepMask(idx) = false;
  end
  schedule = schedule(keepMask);
end
schedule = sort(schedule);
end

function featureTable = localBuildScheduleFeatureTable(strategyList)
%LOCALBUILDSCHEDULEFEATURETABLE Summarize frame-lag structure for every strategy schedule.
rowCell = {};
for iStrategy = 1:numel(strategyList)
  strategy = strategyList(iStrategy);
  for iSchedule = 1:numel(strategy.subsetOffsetCell)
    rowCell{end + 1, 1} = localBuildScheduleFeatureRow(strategy.strategyName, ...
      strategy.subsetLabelList(iSchedule), strategy.subsetOffsetCell{iSchedule}); %#ok<AGROW>
  end
end
if isempty(rowCell)
  featureTable = table();
else
  featureTable = struct2table([rowCell{:}].');
end
end

function row = localBuildScheduleFeatureRow(strategyName, scheduleLabel, offsetList)
%LOCALBUILDSCHEDULEFEATUREROW Build one lag-structure summary row.
offsetList = sort(reshape(double(offsetList), 1, []));
gapList = diff(offsetList);
lagList = localUniquePairwiseLagList(offsetList);
row = struct();
row.strategyName = string(strategyName);
row.scheduleLabel = string(scheduleLabel);
row.scheduleFamily = localInferScheduleFamily(scheduleLabel);
row.numFrame = numel(offsetList);
row.containsZero = any(offsetList == 0);
row.aperture = max(offsetList) - min(offsetList);
row.minGap = localMinOrNaN(gapList);
row.maxGap = localMaxOrNaN(gapList);
row.numUniqueLag = numel(lagList);
row.maxLag = localMaxOrNaN(lagList);
row.gcdLag = localGcdVector(lagList);
if isempty(lagList)
  row.lagRedundancy = NaN;
else
  row.lagRedundancy = nchoosek(numel(offsetList), 2) / numel(lagList);
end
row.offsetList = localJoinNumericList(offsetList);
row.consecutiveGapList = localJoinNumericList(gapList);
end

function family = localInferScheduleFamily(scheduleLabel)
%LOCALINFERSCHEDULEFAMILY Infer a compact family name from a schedule label.
label = lower(string(scheduleLabel));
if startsWith(label, "curated")
  family = "legacyCurated";
elseif startsWith(label, "random")
  family = "random";
elseif startsWith(label, "staggered")
  family = "staggered";
elseif startsWith(label, "sparse")
  family = "sparseRuler";
elseif startsWith(label, "coprime")
  family = "coprime";
else
  family = "custom";
end
end

function lagList = localUniquePairwiseLagList(offsetList)
%LOCALUNIQUEPAIRWISELAGLIST Return sorted positive pairwise frame lags.
lagList = [];
for iA = 1:numel(offsetList)
  for iB = (iA + 1):numel(offsetList)
    lagList(end + 1) = abs(offsetList(iB) - offsetList(iA)); %#ok<AGROW>
  end
end
lagList = unique(lagList(lagList > 0));
end

function value = localGcdVector(valueList)
%LOCALGCDVECTOR Compute gcd over an integer vector with NaN fallback.
valueList = abs(round(double(valueList(:))));
valueList = valueList(isfinite(valueList) & valueList > 0);
if isempty(valueList)
  value = NaN;
  return;
end
value = valueList(1);
for iValue = 2:numel(valueList)
  value = gcd(value, valueList(iValue));
end
end

function txt = localJoinNumericList(valueList)
%LOCALJOINNUMERICLIST Format numeric vectors for compact tables.
valueList = reshape(double(valueList), 1, []);
if isempty(valueList)
  txt = "";
else
  txt = strjoin(string(valueList), ",");
end
end

function value = localMinOrNaN(valueList)
%LOCALMINORNAN Return min(valueList) or NaN for an empty list.
if isempty(valueList)
  value = NaN;
else
  value = min(valueList);
end
end

function value = localMaxOrNaN(valueList)
%LOCALMAXORNAN Return max(valueList) or NaN for an empty list.
if isempty(valueList)
  value = NaN;
else
  value = max(valueList);
end
end

function featureTable = localGetScheduleFeatureTableForPrint(scanData)
%LOCALGETSCHEDULEFEATURETABLEFORPRINT Return schedule feature rows for summary output.
featureTable = table();
if isstruct(scanData) && isfield(scanData, 'scheduleFeatureTable') && istable(scanData.scheduleFeatureTable)
  featureTable = scanData.scheduleFeatureTable;
elseif isstruct(scanData) && isfield(scanData, 'strategyList')
  featureTable = localBuildScheduleFeatureTable(scanData.strategyList);
end
if isempty(featureTable) || ~istable(featureTable)
  return;
end
keepNames = {'strategyName', 'scheduleLabel', 'scheduleFamily', 'numFrame', 'aperture', ...
  'gcdLag', 'numUniqueLag', 'lagRedundancy', 'offsetList'};
keepNames = keepNames(ismember(keepNames, featureTable.Properties.VariableNames));
featureTable = featureTable(:, keepNames);
end

function strategy = localEmptyStrategy()
%LOCALEMPTYSTRATEGY Return a scalar template for subset-bank strategies.
strategy = struct('strategyName', "", 'subsetOffsetCell', {{}}, 'subsetLabelList', strings(0, 1), ...
  'numRandomTrial', 0, 'bankSize', 0);
end


function strategy = localMakeStrategy(name, subsetOffsetCell, subsetLabelList, numRandomTrial)
%LOCALMAKESTRATEGY Package one subset-bank strategy configuration.
strategy = localEmptyStrategy();
strategy.strategyName = string(name);
strategy.subsetOffsetCell = subsetOffsetCell;
strategy.subsetLabelList = reshape(string(subsetLabelList), [], 1);
strategy.numRandomTrial = numRandomTrial;
strategy.bankSize = numel(subsetOffsetCell) + numRandomTrial;
end


function batchOpt = localBuildBatchOpt(scanConfig, strategy)
%LOCALBUILDBATCHOPT Build one simple-flow batch option for a scan strategy.
parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
contextOpt = struct('baseSeed', scanConfig.baseSeed, ...
  'subsetOffsetCell', {strategy.subsetOffsetCell}, ...
  'subsetLabelList', {strategy.subsetLabelList}, ...
  'numSubsetRandomTrial', strategy.numRandomTrial, ...
  'parallelOpt', parallelOpt);
defaultLabels = reshape(strategy.subsetLabelList, 1, []);
if strategy.numRandomTrial > 0
  defaultLabels = [defaultLabels, "random" + string(1:strategy.numRandomTrial)];
end
flowOptOverride = struct();
flowOptOverride.subsetDefaultBankLabelList = defaultLabels;
flowOptOverride.subsetMarginFallbackBankLabelList = strings(1, 0);
flowOptOverride.subsetEnableMarginFallback = false;
flowOptOverride.periodicRefineFdHalfWidthHz = 50;
flowOptOverride.periodicRefineFdRateHalfWidthHzPerSec = 100;
flowOptOverride.periodicRefineDoaHalfWidthDeg = [1e-8; 1e-8];
flowOptOverride.periodicRefineFreezeDoa = true;
flowOptOverride.periodicRefineDoaSeedMode = "dualWhenMulti";
flowOptOverride.periodicPolishEnableWhenMulti = true;
flowOptOverride.periodicPolishDoaHalfWidthDeg = [0.002; 0.002];
flowOpt = buildSimpleDynamicFlowOpt(flowOptOverride);
checkpointOpt = localBuildCheckpointOpt(scanConfig, strategy);
batchOpt = struct('replayName', "MS-MF-Dyn-SubsetBankCoverage-" + strategy.strategyName, ...
  'snrDb', scanConfig.snrDb, 'baseSeed', scanConfig.baseSeed, 'seedList', scanConfig.seedList, ...
  'numRepeat', scanConfig.numRepeat, 'optVerbose', scanConfig.optVerbose, ...
  'contextOpt', contextOpt, 'parallelOpt', parallelOpt, 'flowOpt', flowOpt, ...
  'checkpointOpt', checkpointOpt);
end

function checkpointOpt = localBuildCheckpointOpt(scanConfig, strategy)
%LOCALBUILDCHECKPOINTOPT Build per-strategy repeat checkpoint options.
checkpointOpt = struct('enable', false);
if ~isfield(scanConfig, 'checkpointEnable') || ~scanConfig.checkpointEnable
  return;
end
checkpointOpt.enable = true;
checkpointOpt.runName = string(scanConfig.scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(scanConfig, strategy);
checkpointOpt.outputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.resume = logical(scanConfig.checkpointResume);
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end


function runKey = localBuildCheckpointRunKey(scanConfig, strategy)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable key so interrupted scans can resume.
seedList = reshape(double(scanConfig.seedList), 1, []);
runKey = sprintf('seed%dto%d_rep%d_snr%.2f_%s', seedList(1), seedList(end), ...
  numel(seedList), scanConfig.snrDb, char(strategy.strategyName));
runKey = string(runKey);
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
runKey = replace(runKey, ' ', '');
end


function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve the repository root from this scan script location.
scanDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scanDir)));
end


function runDir = localGetBatchCheckpointDir(batchOpt)
%LOCALGETBATCHCHECKPOINTDIR Return the expected checkpoint directory for logs.
runDir = "";
if ~isfield(batchOpt, 'checkpointOpt') || ~isstruct(batchOpt.checkpointOpt) || ...
    ~isfield(batchOpt.checkpointOpt, 'enable') || ~batchOpt.checkpointOpt.enable
  return;
end
checkpointOpt = batchOpt.checkpointOpt;
runDir = string(fullfile(checkpointOpt.outputRoot, checkpointOpt.runName, checkpointOpt.runKey));
end


function candidateTable = localBuildCandidateTable(repeatCell, strategyName)
%LOCALBUILDCANDIDATETABLE Collect evaluated subset candidates for one strategy.
rowCell = {};
for iRepeat = 1:numel(repeatCell)
  t = repeatCell{iRepeat}.msFlow.subsetCandidateTable;
  if isempty(t) || ~istable(t)
    continue;
  end
  t.strategyName = repmat(strategyName, height(t), 1);
  t.taskSeed = repmat(repeatCell{iRepeat}.taskSeed, height(t), 1);
  rowCell{end + 1, 1} = t; %#ok<AGROW>
end
if isempty(rowCell)
  candidateTable = table();
else
  candidateTable = vertcat(rowCell{:});
  candidateTable = movevars(candidateTable, {'strategyName', 'taskSeed'}, 'Before', 1);
end
end


function outTable = localVertcatTableCell(tableCell)
%LOCALVERTCATTABLECELL Concatenate a cell array of possibly empty tables.
nonEmptyMask = ~cellfun(@isempty, tableCell);
if ~any(nonEmptyMask)
  outTable = table();
else
  outTable = vertcat(tableCell{nonEmptyMask});
end
end


function scanTable = localAnnotateToothMetrics(scanTable, scanConfig)
%LOCALANNOTATETOOTHMETRICS Split integer-tooth and residual-aware hit metrics.
if isempty(scanTable) || ~istable(scanTable) || ~any(strcmp(scanTable.Properties.VariableNames, 'toothIdx'))
  return;
end
numRow = height(scanTable);
if any(strcmp(scanTable.Properties.VariableNames, 'toothResidualHz'))
  residualOk = abs(scanTable.toothResidualHz) <= scanConfig.toothResidualTolHz;
else
  residualOk = true(numRow, 1);
end
scanTable.truthToothIndexHit = scanTable.toothIdx == 0;
scanTable.nearToothIndexHit = abs(scanTable.toothIdx) <= scanConfig.nearToothIdxTol;
scanTable.truthToothStrictHit = scanTable.truthToothIndexHit & residualOk;
scanTable.nearToothStrictHit = scanTable.nearToothIndexHit & residualOk;
scanTable.sameToothResidualFail = scanTable.truthToothIndexHit & ~residualOk;
scanTable.nearToothResidualFail = scanTable.nearToothIndexHit & ~residualOk;
scanTable.truthToothHit = scanTable.truthToothStrictHit;
scanTable.nearToothHit = scanTable.nearToothStrictHit;
end


function candidateTable = localAnnotateCandidateToothMetrics(candidateTable, scanConfig)
%LOCALANNOTATECANDIDATETOOTHMETRICS Add candidate tooth coverage flags when possible.
if isempty(candidateTable) || ~istable(candidateTable) || ~any(strcmp(candidateTable.Properties.VariableNames, 'toothIdx'))
  return;
end
numRow = height(candidateTable);
residualAvailable = any(strcmp(candidateTable.Properties.VariableNames, 'toothResidualHz'));
if residualAvailable
  residualOk = abs(candidateTable.toothResidualHz) <= scanConfig.toothResidualTolHz;
else
  residualOk = true(numRow, 1);
end
candidateTable.truthToothIndexHit = candidateTable.toothIdx == 0;
candidateTable.nearToothIndexHit = abs(candidateTable.toothIdx) <= scanConfig.nearToothIdxTol;
candidateTable.truthToothHit = candidateTable.truthToothIndexHit & residualOk;
candidateTable.nearToothHit = candidateTable.nearToothIndexHit & residualOk;
candidateTable.toothResidualMetricAvailable = repmat(residualAvailable, numRow, 1);
end


function scanTable = localAnnotateBaselineComparison(scanTable, scanConfig)
%LOCALANNOTATEBASELINECOMPARISON Add offline baseline transition and damage flags.
scanTable = localAnnotateToothMetrics(scanTable, scanConfig);
numRow = height(scanTable);
scanTable.baselineEasy = false(numRow, 1);
scanTable.easyDamageFromBaseline = false(numRow, 1);
scanTable.angleDeltaFromBaselineDeg = nan(numRow, 1);
scanTable.fdRefAbsErrDeltaFromBaselineHz = nan(numRow, 1);
scanTable.truthToothGainFromBaseline = false(numRow, 1);
scanTable.nearToothGainFromBaseline = false(numRow, 1);
scanTable.truthToothIndexGainFromBaseline = false(numRow, 1);
scanTable.nearToothIndexGainFromBaseline = false(numRow, 1);
scanTable.sameToothResidualFixFromBaseline = false(numRow, 1);
scanTable.transitionClass = strings(numRow, 1);

baselineMask = scanTable.strategyName == scanConfig.baselineStrategyName;
if ~any(baselineMask)
  scanTable.transitionClass(:) = "noBaseline";
  return;
end
baselineTable = scanTable(baselineMask, :);
for iRow = 1:numRow
  seed = scanTable.taskSeed(iRow);
  idxBase = find(baselineTable.taskSeed == seed, 1, 'first');
  if isempty(idxBase)
    scanTable.transitionClass(iRow) = "missingBaselineSeed";
    continue;
  end
  baseRow = baselineTable(idxBase, :);
  baseTruthStrictHit = baseRow.truthToothStrictHit;
  baseNearStrictHit = baseRow.nearToothStrictHit;
  baseTruthIndexHit = baseRow.truthToothIndexHit;
  baseNearIndexHit = baseRow.nearToothIndexHit;
  baseSameToothResidualFail = baseRow.sameToothResidualFail;
  baseAngle = baseRow.angleErrDeg;
  isEasy = baseTruthStrictHit && baseAngle <= scanConfig.easyAngleTolDeg;
  scanTable.baselineEasy(iRow) = isEasy;
  scanTable.angleDeltaFromBaselineDeg(iRow) = scanTable.angleErrDeg(iRow) - baseAngle;
  scanTable.fdRefAbsErrDeltaFromBaselineHz(iRow) = abs(scanTable.fdRefErrHz(iRow)) - abs(baseRow.fdRefErrHz);
  scanTable.truthToothGainFromBaseline(iRow) = scanTable.truthToothStrictHit(iRow) && ~baseTruthStrictHit;
  scanTable.nearToothGainFromBaseline(iRow) = scanTable.nearToothStrictHit(iRow) && ~baseNearStrictHit;
  scanTable.truthToothIndexGainFromBaseline(iRow) = scanTable.truthToothIndexHit(iRow) && ~baseTruthIndexHit;
  scanTable.nearToothIndexGainFromBaseline(iRow) = scanTable.nearToothIndexHit(iRow) && ~baseNearIndexHit;
  scanTable.sameToothResidualFixFromBaseline(iRow) = baseSameToothResidualFail && scanTable.truthToothStrictHit(iRow);
  scanTable.easyDamageFromBaseline(iRow) = isEasy && ...
    (~scanTable.truthToothStrictHit(iRow) || scanTable.angleDeltaFromBaselineDeg(iRow) > scanConfig.damageAngleTolDeg);
  scanTable.transitionClass(iRow) = localClassifyBaselineTransition(baseRow, scanTable(iRow, :), scanConfig.baselineStrategyName);
end
end


function transitionClass = localClassifyBaselineTransition(baseRow, row, baselineStrategyName)
%LOCALCLASSIFYBASELINETRANSITION Classify one selected-result transition versus baseline.
if row.strategyName == baselineStrategyName
  transitionClass = "baseline";
elseif ~baseRow.truthToothStrictHit && row.truthToothStrictHit
  transitionClass = "rescueToTruthStrict";
elseif ~baseRow.nearToothStrictHit && row.nearToothStrictHit
  transitionClass = "rescueToNearStrict";
elseif ~baseRow.truthToothIndexHit && row.truthToothIndexHit && ~row.truthToothStrictHit
  transitionClass = "rescueToTruthIndexResidualFail";
elseif baseRow.truthToothIndexHit && baseRow.sameToothResidualFail && row.truthToothStrictHit
  transitionClass = "sameToothResidualFix";
elseif baseRow.truthToothStrictHit && ~row.truthToothStrictHit
  transitionClass = "truthStrictDamage";
elseif baseRow.nearToothStrictHit && ~row.nearToothStrictHit
  transitionClass = "nearStrictDamage";
elseif row.easyDamageFromBaseline
  transitionClass = "easyAngleDamage";
else
  transitionClass = "noStrictChange";
end
end


function transitionTable = localBuildBaselineTransitionTable(scanTable, scanConfig)
%LOCALBUILDBASELINETRANSITIONTABLE Build per-seed baseline-to-strategy transition rows.
if isempty(scanTable) || ~istable(scanTable)
  transitionTable = table();
  return;
end
scanTable = localAnnotateBaselineComparison(scanTable, scanConfig);
baselineMask = scanTable.strategyName == scanConfig.baselineStrategyName;
if ~any(baselineMask)
  transitionTable = table();
  return;
end
baselineTable = scanTable(baselineMask, :);
rowList = repmat(localEmptyTransitionRow(), height(scanTable), 1);
for iRow = 1:height(scanTable)
  row = localEmptyTransitionRow();
  seed = scanTable.taskSeed(iRow);
  idxBase = find(baselineTable.taskSeed == seed, 1, 'first');
  if isempty(idxBase)
    continue;
  end
  baseRow = baselineTable(idxBase, :);
  row.strategyName = scanTable.strategyName(iRow);
  row.taskSeed = seed;
  row.baselineSubsetLabel = string(baseRow.selectedSubsetLabel);
  row.selectedSubsetLabel = string(scanTable.selectedSubsetLabel(iRow));
  row.baselineToothIdx = baseRow.toothIdx;
  row.strategyToothIdx = scanTable.toothIdx(iRow);
  row.baselineResidualHz = baseRow.toothResidualHz;
  row.strategyResidualHz = scanTable.toothResidualHz(iRow);
  row.baselineTruthIndexHit = baseRow.truthToothIndexHit;
  row.strategyTruthIndexHit = scanTable.truthToothIndexHit(iRow);
  row.baselineTruthStrictHit = baseRow.truthToothStrictHit;
  row.strategyTruthStrictHit = scanTable.truthToothStrictHit(iRow);
  row.baselineNearStrictHit = baseRow.nearToothStrictHit;
  row.strategyNearStrictHit = scanTable.nearToothStrictHit(iRow);
  row.baselineSameToothResidualFail = baseRow.sameToothResidualFail;
  row.strategySameToothResidualFail = scanTable.sameToothResidualFail(iRow);
  row.transitionClass = scanTable.transitionClass(iRow);
  row.baselineAngleErrDeg = baseRow.angleErrDeg;
  row.strategyAngleErrDeg = scanTable.angleErrDeg(iRow);
  row.angleDeltaDeg = scanTable.angleDeltaFromBaselineDeg(iRow);
  row.fdRefAbsErrDeltaHz = scanTable.fdRefAbsErrDeltaFromBaselineHz(iRow);
  rowList(iRow) = row;
end
transitionTable = struct2table(rowList);
end


function row = localEmptyTransitionRow()
%LOCALEMPTYTRANSITIONROW Return a scalar template for baseline transition rows.
row = struct('strategyName', "", 'taskSeed', NaN, 'baselineSubsetLabel', "", ...
  'selectedSubsetLabel', "", 'baselineToothIdx', NaN, 'strategyToothIdx', NaN, ...
  'baselineResidualHz', NaN, 'strategyResidualHz', NaN, ...
  'baselineTruthIndexHit', false, 'strategyTruthIndexHit', false, ...
  'baselineTruthStrictHit', false, 'strategyTruthStrictHit', false, ...
  'baselineNearStrictHit', false, 'strategyNearStrictHit', false, ...
  'baselineSameToothResidualFail', false, 'strategySameToothResidualFail', false, ...
  'transitionClass', "", 'baselineAngleErrDeg', NaN, 'strategyAngleErrDeg', NaN, ...
  'angleDeltaDeg', NaN, 'fdRefAbsErrDeltaHz', NaN);
end


function coverageTable = localBuildCandidateSeedCoverageTable(scanTable, candidateTable)
%LOCALBUILDCANDIDATESEEDCOVERAGETABLE Summarize per-seed candidate coverage and selected misses.
if isempty(scanTable) || ~istable(scanTable)
  coverageTable = table();
  return;
end
rowList = repmat(localEmptyCandidateSeedCoverageRow(), height(scanTable), 1);
for iRow = 1:height(scanTable)
  row = localEmptyCandidateSeedCoverageRow();
  row.strategyName = scanTable.strategyName(iRow);
  row.taskSeed = scanTable.taskSeed(iRow);
  row.selectedTruthIndexHit = scanTable.truthToothIndexHit(iRow);
  row.selectedTruthStrictHit = scanTable.truthToothStrictHit(iRow);
  row.selectedNearStrictHit = scanTable.nearToothStrictHit(iRow);
  if ~isempty(candidateTable) && istable(candidateTable) && ...
      all(ismember({'strategyName', 'taskSeed'}, candidateTable.Properties.VariableNames))
    cmask = candidateTable.strategyName == row.strategyName & candidateTable.taskSeed == row.taskSeed;
    row.candidateEvalCountPerSeed = nnz(cmask);
    if row.candidateEvalCountPerSeed > 0 && any(strcmp(candidateTable.Properties.VariableNames, 'truthToothIndexHit'))
      c = candidateTable(cmask, :);
      row.anyCandidateTruthIndexHit = any(c.truthToothIndexHit);
      row.anyCandidateNearIndexHit = any(c.nearToothIndexHit);
      row.anyCandidateTruthStrictHit = any(c.truthToothHit);
      row.anyCandidateNearStrictHit = any(c.nearToothHit);
      row.bestCandidateAbsToothIdx = min(abs(c.toothIdx), [], 'omitnan');
      row.selectedMissDespiteTruthIndexCandidate = row.anyCandidateTruthIndexHit && ~row.selectedTruthIndexHit;
      row.selectedMissDespiteTruthStrictCandidate = row.anyCandidateTruthStrictHit && ~row.selectedTruthStrictHit;
      row.selectedMissDespiteNearStrictCandidate = row.anyCandidateNearStrictHit && ~row.selectedNearStrictHit;
    end
  end
  rowList(iRow) = row;
end
coverageTable = struct2table(rowList);
end


function row = localEmptyCandidateSeedCoverageRow()
%LOCALEMPTYCANDIDATESEEDCOVERAGEROW Return a scalar candidate coverage template.
row = struct('strategyName', "", 'taskSeed', NaN, 'candidateEvalCountPerSeed', 0, ...
  'anyCandidateTruthIndexHit', false, 'anyCandidateNearIndexHit', false, ...
  'anyCandidateTruthStrictHit', false, 'anyCandidateNearStrictHit', false, ...
  'bestCandidateAbsToothIdx', NaN, 'selectedTruthIndexHit', false, ...
  'selectedTruthStrictHit', false, 'selectedNearStrictHit', false, ...
  'selectedMissDespiteTruthIndexCandidate', false, ...
  'selectedMissDespiteTruthStrictCandidate', false, ...
  'selectedMissDespiteNearStrictCandidate', false);
end


function scanConfig = localAttachFrameTimingConfig(scanConfig, contextCell)
%LOCALATTACHFRAMETIMINGCONFIG Attach frame interval and tooth step for diagnostics.
scanConfig = localEnsureConsensusConfig(scanConfig, struct());
frameIntvlSec = localGetFieldOrDefault(scanConfig, 'frameIntvlSec', NaN);
for iContext = 1:numel(contextCell)
  contextUse = contextCell{iContext};
  if isstruct(contextUse) && isfield(contextUse, 'frameIntvlSec') && ...
      isfinite(contextUse.frameIntvlSec) && contextUse.frameIntvlSec > 0
    frameIntvlSec = contextUse.frameIntvlSec;
    break;
  end
end
scanConfig.frameIntvlSec = frameIntvlSec;
if isfinite(frameIntvlSec) && frameIntvlSec > 0
  scanConfig.toothStepHz = 1 / frameIntvlSec;
else
  scanConfig.toothStepHz = localGetFieldOrDefault(scanConfig, 'toothStepHz', NaN);
end
end


function scanConfig = localEnsureConsensusConfig(scanConfig, scanData)
%LOCALENSURECONSENSUSCONFIG Fill missing diagnostic fields for old snapshots.
if ~isfield(scanConfig, 'toothConsensusResidualPenaltyHz') || isempty(scanConfig.toothConsensusResidualPenaltyHz)
  if isfield(scanConfig, 'toothConsensusResidualTolHz') && ~isempty(scanConfig.toothConsensusResidualTolHz)
    scanConfig.toothConsensusResidualPenaltyHz = scanConfig.toothConsensusResidualTolHz;
  else
    scanConfig.toothConsensusResidualPenaltyHz = 75;
  end
end
if ~isfield(scanConfig, 'frameIntvlSec') || isempty(scanConfig.frameIntvlSec)
  scanConfig.frameIntvlSec = NaN;
end
if ~isfield(scanConfig, 'toothStepHz') || isempty(scanConfig.toothStepHz)
  scanConfig.toothStepHz = NaN;
end
if isstruct(scanData) && isfield(scanData, 'contextSummary') && isstruct(scanData.contextSummary) && ...
    isfield(scanData.contextSummary, 'frameIntvlSec') && isfinite(scanData.contextSummary.frameIntvlSec)
  scanConfig.frameIntvlSec = scanData.contextSummary.frameIntvlSec;
end
if ~isfinite(scanConfig.toothStepHz) && isfinite(scanConfig.frameIntvlSec) && scanConfig.frameIntvlSec > 0
  scanConfig.toothStepHz = 1 / scanConfig.frameIntvlSec;
end
end


function candidateTable = localAnnotateCandidateAnchorGroups(candidateTable, repeatCellByStrategy, scanConfig)
%LOCALANNOTATECANDIDATEANCHORGROUPS Add no-truth static-anchor tooth groups.
if isempty(candidateTable) || ~istable(candidateTable)
  return;
end
numRow = height(candidateTable);
anchorFdRefHz = nan(numRow, 1);
anchorToothGroupIdx = nan(numRow, 1);
anchorToothResidualHz = nan(numRow, 1);
toothStepHz = localGetFieldOrDefault(scanConfig, 'toothStepHz', NaN);
anchorTable = localBuildAnchorTable(repeatCellByStrategy);
for iRow = 1:numRow
  anchor = localLookupAnchorFdRef(anchorTable, candidateTable.strategyName(iRow), candidateTable.taskSeed(iRow));
  if ~isfinite(anchor)
    mask = candidateTable.strategyName == candidateTable.strategyName(iRow) & candidateTable.taskSeed == candidateTable.taskSeed(iRow);
    anchor = median(candidateTable.fdRefEstHz(mask), 'omitnan');
  end
  anchorFdRefHz(iRow) = anchor;
  if isfinite(anchor) && isfinite(toothStepHz) && toothStepHz > 0 && isfinite(candidateTable.fdRefEstHz(iRow))
    fdDrift = candidateTable.fdRefEstHz(iRow) - anchor;
    anchorToothGroupIdx(iRow) = round(fdDrift / toothStepHz);
    anchorToothResidualHz(iRow) = fdDrift - anchorToothGroupIdx(iRow) * toothStepHz;
  end
end
candidateTable.anchorFdRefHz = anchorFdRefHz;
candidateTable.anchorToothGroupIdx = anchorToothGroupIdx;
candidateTable.anchorToothResidualHz = anchorToothResidualHz;
end


function anchorTable = localBuildAnchorTable(repeatCellByStrategy)
%LOCALBUILDANCHORTABLE Collect static-seed fdRef anchors from repeat cells.
rowCell = {};
for iStrategy = 1:numel(repeatCellByStrategy)
  repeatCell = repeatCellByStrategy{iStrategy};
  if isempty(repeatCell)
    continue;
  end
  for iRepeat = 1:numel(repeatCell)
    repeatUse = repeatCell{iRepeat};
    if ~isstruct(repeatUse) || ~isfield(repeatUse, 'msFlow') || ~isfield(repeatUse, 'taskSeed')
      continue;
    end
    row = struct();
    row.strategyName = string(localGetFieldOrDefault(repeatUse.msFlow, 'selectedSubsetLabel', ""));
    row.taskSeed = repeatUse.taskSeed;
    row.anchorFdRefHz = localExtractStaticAnchorFdRef(repeatUse);
    rowCell{end + 1, 1} = row; %#ok<AGROW>
  end
end
if isempty(rowCell)
  anchorTable = table();
else
  anchorTable = struct2table([rowCell{:}].');
end
end


function anchorFdRefHz = localExtractStaticAnchorFdRef(repeatUse)
%LOCALEXTRACTSTATICANCHORFDREF Read the static MS fdRef estimate without truth.
anchorFdRefHz = NaN;
if ~isstruct(repeatUse) || ~isfield(repeatUse, 'staticBundle') || ~isstruct(repeatUse.staticBundle)
  return;
end
staticBundle = repeatUse.staticBundle;
if ~isfield(staticBundle, 'caseStaticMs') || ~isstruct(staticBundle.caseStaticMs)
  return;
end
estResult = localGetFieldOrDefault(staticBundle.caseStaticMs, 'estResult', struct());
anchorFdRefHz = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
end


function anchor = localLookupAnchorFdRef(anchorTable, strategyName, taskSeed)
%LOCALLOOKUPANCHORFDREF Look up a no-truth anchor for one strategy/seed pair.
anchor = NaN;
if isempty(anchorTable) || ~istable(anchorTable) || ~all(ismember({'strategyName', 'taskSeed'}, anchorTable.Properties.VariableNames))
  return;
end
mask = anchorTable.taskSeed == taskSeed;
idx = find(mask & isfinite(anchorTable.anchorFdRefHz), 1, 'first');
if isempty(idx)
  return;
end
anchor = anchorTable.anchorFdRefHz(idx);
end


function consensusTable = localBuildConsensusDiagnosticTable(scanTable, candidateTable, scanConfig)
%LOCALBUILDCONSENSUSDIAGNOSTICTABLE Evaluate a no-truth anchor-group selector offline.
if isempty(scanTable) || ~istable(scanTable)
  consensusTable = table();
  return;
end
rowList = repmat(localEmptyConsensusRow(), height(scanTable), 1);
for iRow = 1:height(scanTable)
  row = localEmptyConsensusRow();
  row.strategyName = scanTable.strategyName(iRow);
  row.taskSeed = scanTable.taskSeed(iRow);
  row.currentTruthToothIndexHit = scanTable.truthToothIndexHit(iRow);
  row.currentTruthToothHit = scanTable.truthToothStrictHit(iRow);
  row.currentNearToothHit = scanTable.nearToothStrictHit(iRow);
  row.currentSelectedSubsetLabel = string(scanTable.selectedSubsetLabel(iRow));
  if isempty(candidateTable) || ~istable(candidateTable)
    rowList(iRow) = row;
    continue;
  end
  mask = candidateTable.strategyName == row.strategyName & candidateTable.taskSeed == row.taskSeed;
  c = candidateTable(mask, :);
  row.numCandidate = height(c);
  if height(c) == 0
    rowList(iRow) = row;
    continue;
  end
  [bestCandidateIdx, groupDiag] = localSelectConsensusCandidate(c, scanConfig);
  row.numAnchorGroup = groupDiag.numAnchorGroup;
  row.selectedAnchorGroupIdx = groupDiag.selectedAnchorGroupIdx;
  row.selectedAnchorGroupSupport = groupDiag.selectedAnchorGroupSupport;
  row.selectedAnchorGroupBestRank = groupDiag.selectedAnchorGroupBestRank;
  row.selectedAnchorGroupBestHealthBucket = groupDiag.selectedAnchorGroupBestHealthBucket;
  row.selectedAnchorGroupBestObj = groupDiag.selectedAnchorGroupBestObj;
  row.selectedAnchorGroupBestResidualNorm = groupDiag.selectedAnchorGroupBestResidualNorm;
  if isnan(bestCandidateIdx)
    rowList(iRow) = row;
    continue;
  end
  best = c(bestCandidateIdx, :);
  row.consensusSubsetLabel = string(localGetTableScalar(best, 'subsetLabel', ""));
  row.consensusToothIdx = localGetTableScalar(best, 'toothIdx', NaN);
  row.consensusToothResidualHz = localGetTableScalar(best, 'toothResidualHz', NaN);
  row.consensusTruthToothIndexHit = logical(localGetTableScalar(best, 'truthToothIndexHit', false));
  row.consensusNearToothIndexHit = logical(localGetTableScalar(best, 'nearToothIndexHit', false));
  row.consensusTruthToothHit = logical(localGetTableScalar(best, 'truthToothHit', false));
  row.consensusNearToothHit = logical(localGetTableScalar(best, 'nearToothHit', false));
  row.consensusTruthGainVsSelected = row.consensusTruthToothHit && ~row.currentTruthToothHit;
  row.consensusNearGainVsSelected = row.consensusNearToothHit && ~row.currentNearToothHit;
  row.consensusTruthDamageVsSelected = row.currentTruthToothHit && ~row.consensusTruthToothHit;
  row.consensusNearDamageVsSelected = row.currentNearToothHit && ~row.consensusNearToothHit;
  rowList(iRow) = row;
end
consensusTable = struct2table(rowList);
end


function [bestCandidateIdx, groupDiag] = localSelectConsensusCandidate(candidateRows, scanConfig)
%LOCALSELECTCONSENSUSCANDIDATE Select one candidate through no-truth group consensus.
groupDiag = localEmptyConsensusGroupDiag();
bestCandidateIdx = NaN;
if isempty(candidateRows)
  return;
end
groupIdx = candidateRows.anchorToothGroupIdx;
finiteMask = isfinite(groupIdx);
if ~any(finiteMask)
  [~, bestCandidateIdx] = min(candidateRows.rank);
  return;
end
groupList = unique(groupIdx(finiteMask), 'stable');
groupDiag.numAnchorGroup = numel(groupList);
bestGroupScore = [];
for iGroup = 1:numel(groupList)
  g = groupList(iGroup);
  gMask = finiteMask & (groupIdx == g);
  idxList = find(gMask);
  [bestInGroupLocal, groupScore] = localBuildConsensusGroupScore(candidateRows, idxList, g, scanConfig);
  if isempty(bestGroupScore) || localIsLexicographicScoreBetter(groupScore, bestGroupScore)
    bestGroupScore = groupScore;
    bestCandidateIdx = bestInGroupLocal;
    groupDiag.selectedAnchorGroupIdx = g;
    groupDiag.selectedAnchorGroupSupport = nnz(gMask);
    groupDiag.selectedAnchorGroupBestRank = candidateRows.rank(bestInGroupLocal);
    groupDiag.selectedAnchorGroupBestHealthBucket = candidateRows.healthBucket(bestInGroupLocal);
    groupDiag.selectedAnchorGroupBestObj = candidateRows.finalObj(bestInGroupLocal);
    groupDiag.selectedAnchorGroupBestResidualNorm = candidateRows.finalResidualNorm(bestInGroupLocal);
  end
end
end


function [bestIdx, scoreVec] = localBuildConsensusGroupScore(candidateRows, idxList, groupIdx, scanConfig)
%LOCALBUILDCONSENSUSGROUPSCORE Build one no-truth group score.
bestIdx = idxList(1);
bestScore = [];
penaltyScaleHz = localGetFieldOrDefault(scanConfig, 'toothConsensusResidualPenaltyHz', 75);
if ~isfinite(penaltyScaleHz) || penaltyScaleHz <= 0
  penaltyScaleHz = 75;
end
for iIdx = 1:numel(idxList)
  idx = idxList(iIdx);
  residualPenalty = abs(candidateRows.anchorToothResidualHz(idx)) / penaltyScaleHz;
  rowScore = [candidateRows.healthBucket(idx), residualPenalty, candidateRows.rank(idx), ...
    candidateRows.finalObj(idx), candidateRows.finalResidualNorm(idx)];
  if isempty(bestScore) || localIsLexicographicScoreBetter(rowScore, bestScore)
    bestScore = rowScore;
    bestIdx = idx;
  end
end
supportCount = numel(idxList);
scoreVec = [bestScore(1), -supportCount, bestScore(2), abs(groupIdx), bestScore(3:end)];
end


function tf = localIsLexicographicScoreBetter(scoreA, scoreB)
%LOCALISLEXICOGRAPHICSCOREBETTER Compare two row score vectors with finite fallback.
tf = false;
for iScore = 1:min(numel(scoreA), numel(scoreB))
  a = scoreA(iScore);
  b = scoreB(iScore);
  if ~isfinite(a), a = inf; end
  if ~isfinite(b), b = inf; end
  if a < b
    tf = true;
    return;
  elseif a > b
    return;
  end
end
end


function diag = localEmptyConsensusGroupDiag()
%LOCALEMPTYCONSENSUSGROUPDIAG Return a scalar selected-group diagnostic template.
diag = struct('numAnchorGroup', 0, 'selectedAnchorGroupIdx', NaN, ...
  'selectedAnchorGroupSupport', 0, 'selectedAnchorGroupBestRank', NaN, ...
  'selectedAnchorGroupBestHealthBucket', NaN, 'selectedAnchorGroupBestObj', NaN, ...
  'selectedAnchorGroupBestResidualNorm', NaN);
end


function row = localEmptyConsensusRow()
%LOCALEMPTYCONSENSUSROW Return a scalar consensus diagnostic row template.
row = struct('strategyName', "", 'taskSeed', NaN, 'numCandidate', 0, ...
  'numAnchorGroup', 0, 'selectedAnchorGroupIdx', NaN, 'selectedAnchorGroupSupport', 0, ...
  'selectedAnchorGroupBestRank', NaN, 'selectedAnchorGroupBestHealthBucket', NaN, ...
  'selectedAnchorGroupBestObj', NaN, 'selectedAnchorGroupBestResidualNorm', NaN, ...
  'currentSelectedSubsetLabel', "", 'consensusSubsetLabel', "", ...
  'currentTruthToothIndexHit', false, 'currentTruthToothHit', false, ...
  'currentNearToothHit', false, 'consensusToothIdx', NaN, ...
  'consensusToothResidualHz', NaN, 'consensusTruthToothIndexHit', false, ...
  'consensusNearToothIndexHit', false, 'consensusTruthToothHit', false, ...
  'consensusNearToothHit', false, 'consensusTruthGainVsSelected', false, ...
  'consensusNearGainVsSelected', false, 'consensusTruthDamageVsSelected', false, ...
  'consensusNearDamageVsSelected', false);
end


function aggregateTable = localBuildConsensusAggregateTable(consensusTable)
%LOCALBUILDCONSENSUSAGGREGATETABLE Summarize diagnostic consensus hit rates.
if isempty(consensusTable) || ~istable(consensusTable)
  aggregateTable = table();
  return;
end
strategyList = unique(consensusTable.strategyName, 'stable');
rowList = repmat(localEmptyConsensusAggregateRow(), numel(strategyList), 1);
for iStrategy = 1:numel(strategyList)
  name = strategyList(iStrategy);
  mask = consensusTable.strategyName == name;
  row = localEmptyConsensusAggregateRow();
  row.strategyName = name;
  row.numRepeat = nnz(mask);
  row.meanCandidateCount = mean(consensusTable.numCandidate(mask), 'omitnan');
  row.meanAnchorGroupCount = mean(consensusTable.numAnchorGroup(mask), 'omitnan');
  row.meanSelectedGroupSupport = mean(consensusTable.selectedAnchorGroupSupport(mask), 'omitnan');
  row.currentTruthToothIndexHitRate = mean(double(consensusTable.currentTruthToothIndexHit(mask)), 'omitnan');
  row.consensusTruthToothIndexHitRate = mean(double(consensusTable.consensusTruthToothIndexHit(mask)), 'omitnan');
  row.currentTruthToothHitRate = mean(double(consensusTable.currentTruthToothHit(mask)), 'omitnan');
  row.consensusTruthToothHitRate = mean(double(consensusTable.consensusTruthToothHit(mask)), 'omitnan');
  row.currentNearToothHitRate = mean(double(consensusTable.currentNearToothHit(mask)), 'omitnan');
  row.consensusNearToothHitRate = mean(double(consensusTable.consensusNearToothHit(mask)), 'omitnan');
  row.consensusTruthGainRate = mean(double(consensusTable.consensusTruthGainVsSelected(mask)), 'omitnan');
  row.consensusNearGainRate = mean(double(consensusTable.consensusNearGainVsSelected(mask)), 'omitnan');
  row.consensusTruthDamageRate = mean(double(consensusTable.consensusTruthDamageVsSelected(mask)), 'omitnan');
  row.consensusNearDamageRate = mean(double(consensusTable.consensusNearDamageVsSelected(mask)), 'omitnan');
  row.consensusLabelSummary = localBuildSelectedLabelSummary(consensusTable.consensusSubsetLabel(mask));
  rowList(iStrategy) = row;
end
aggregateTable = struct2table(rowList);
end


function row = localEmptyConsensusAggregateRow()
%LOCALEMPTYCONSENSUSAGGREGATEROW Return a scalar consensus aggregate template.
row = struct('strategyName', "", 'numRepeat', 0, 'meanCandidateCount', NaN, ...
  'meanAnchorGroupCount', NaN, 'meanSelectedGroupSupport', NaN, ...
  'currentTruthToothIndexHitRate', NaN, 'consensusTruthToothIndexHitRate', NaN, ...
  'currentTruthToothHitRate', NaN, 'consensusTruthToothHitRate', NaN, ...
  'currentNearToothHitRate', NaN, 'consensusNearToothHitRate', NaN, ...
  'consensusTruthGainRate', NaN, 'consensusNearGainRate', NaN, ...
  'consensusTruthDamageRate', NaN, 'consensusNearDamageRate', NaN, ...
  'consensusLabelSummary', "");
end

function aggregateTable = localBuildAggregateTable(scanTable, candidateTable, candidateSeedCoverageTable, scanConfig)
%LOCALBUILDAGGREGATETABLE Summarize hit rate, damage, accuracy, runtime, and coverage.
strategyNameList = unique(scanTable.strategyName, 'stable');
rowList = repmat(localEmptyAggregateRow(), numel(strategyNameList), 1);
for iStrategy = 1:numel(strategyNameList)
  name = strategyNameList(iStrategy);
  mask = scanTable.strategyName == name;
  row = localEmptyAggregateRow();
  row.strategyName = name;
  row.bankSize = scanTable.bankSize(find(mask, 1, 'first'));
  row.numRandomTrial = scanTable.numRandomTrial(find(mask, 1, 'first'));
  row.numRepeat = nnz(mask);
  row.resolvedRate = mean(double(scanTable.isResolved(mask)), 'omitnan');
  row.truthToothIndexHitRate = mean(double(scanTable.truthToothIndexHit(mask)), 'omitnan');
  row.nearToothIndexHitRate = mean(double(scanTable.nearToothIndexHit(mask)), 'omitnan');
  row.truthToothHitRate = mean(double(scanTable.truthToothStrictHit(mask)), 'omitnan');
  row.nearToothHitRate = mean(double(scanTable.nearToothStrictHit(mask)), 'omitnan');
  row.sameToothResidualFailRate = mean(double(scanTable.sameToothResidualFail(mask)), 'omitnan');
  row.nearToothResidualFailRate = mean(double(scanTable.nearToothResidualFail(mask)), 'omitnan');
  row.truthToothGainRate = mean(double(scanTable.truthToothGainFromBaseline(mask)), 'omitnan');
  row.nearToothGainRate = mean(double(scanTable.nearToothGainFromBaseline(mask)), 'omitnan');
  row.truthToothIndexGainRate = mean(double(scanTable.truthToothIndexGainFromBaseline(mask)), 'omitnan');
  row.nearToothIndexGainRate = mean(double(scanTable.nearToothIndexGainFromBaseline(mask)), 'omitnan');
  row.sameToothResidualFixRate = mean(double(scanTable.sameToothResidualFixFromBaseline(mask)), 'omitnan');
  row.easyDamageRate = mean(double(scanTable.easyDamageFromBaseline(mask)), 'omitnan');
  row.easyDamageCount = nnz(scanTable.easyDamageFromBaseline(mask));
  row.angleRmseDeg = sqrt(mean(scanTable.angleErrDeg(mask).^2, 'omitnan'));
  row.angleP95Deg = prctile(scanTable.angleErrDeg(mask), 95);
  row.angleMaxDeg = max(scanTable.angleErrDeg(mask), [], 'omitnan');
  row.fdRefRmseHz = sqrt(mean(scanTable.fdRefErrHz(mask).^2, 'omitnan'));
  row.meanRunTimeMs = mean(scanTable.runTimeMs(mask), 'omitnan');
  row.selectedLabelSummary = localBuildSelectedLabelSummary(scanTable.selectedSubsetLabel(mask));
  if ~isempty(candidateSeedCoverageTable) && istable(candidateSeedCoverageTable)
    smask = candidateSeedCoverageTable.strategyName == name;
    row.candidateAnyTruthIndexRate = mean(double(candidateSeedCoverageTable.anyCandidateTruthIndexHit(smask)), 'omitnan');
    row.candidateAnyNearIndexRate = mean(double(candidateSeedCoverageTable.anyCandidateNearIndexHit(smask)), 'omitnan');
    row.candidateAnyTruthToothRate = mean(double(candidateSeedCoverageTable.anyCandidateTruthStrictHit(smask)), 'omitnan');
    row.candidateAnyNearToothRate = mean(double(candidateSeedCoverageTable.anyCandidateNearStrictHit(smask)), 'omitnan');
    row.selectedMissDespiteTruthCandidateRate = mean(double(candidateSeedCoverageTable.selectedMissDespiteTruthStrictCandidate(smask)), 'omitnan');
    row.selectedMissDespiteNearCandidateRate = mean(double(candidateSeedCoverageTable.selectedMissDespiteNearStrictCandidate(smask)), 'omitnan');
  end
  if ~isempty(candidateTable)
    cmask = candidateTable.strategyName == name;
    row.candidateEvalCount = nnz(cmask);
    row.meanCandidateEvalPerRepeat = row.candidateEvalCount / max(1, row.numRepeat);
    if any(strcmp(candidateTable.Properties.VariableNames, 'toothIdx'))
      row.candidateTruthToothIndexCount = nnz(candidateTable.truthToothIndexHit(cmask));
      row.candidateNearToothIndexCount = nnz(candidateTable.nearToothIndexHit(cmask));
      row.candidateTruthToothCount = localCountCandidateTooth(candidateTable(cmask, :), 0, scanConfig.toothResidualTolHz);
      row.candidateNearToothCount = localCountCandidateNearTooth(candidateTable(cmask, :), ...
        scanConfig.nearToothIdxTol, scanConfig.toothResidualTolHz);
      row.candidateTruthToothRate = row.candidateTruthToothCount / max(1, row.candidateEvalCount);
      row.candidateNearToothRate = row.candidateNearToothCount / max(1, row.candidateEvalCount);
      row.candidateTruthToothIndexRate = row.candidateTruthToothIndexCount / max(1, row.candidateEvalCount);
      row.candidateNearToothIndexRate = row.candidateNearToothIndexCount / max(1, row.candidateEvalCount);
    end
    if any(strcmp(candidateTable.Properties.VariableNames, 'subsetLabel'))
      row.candidateLabelSummary = localBuildSelectedLabelSummary(candidateTable.subsetLabel(cmask));
    end
  end
  rowList(iStrategy) = row;
end
aggregateTable = struct2table(rowList);
end


function count = localCountCandidateTooth(candidateTable, toothIdx, residualTolHz)
%LOCALCOUNTCANDIDATETOOTH Count candidates on a target truth tooth.
if isempty(candidateTable)
  count = 0;
  return;
end
mask = candidateTable.toothIdx == toothIdx;
if any(strcmp(candidateTable.Properties.VariableNames, 'toothResidualHz'))
  mask = mask & abs(candidateTable.toothResidualHz) <= residualTolHz;
end
count = nnz(mask);
end


function count = localCountCandidateNearTooth(candidateTable, toothIdxTol, residualTolHz)
%LOCALCOUNTCANDIDATENEARTOOTH Count candidates within a near-tooth band.
if isempty(candidateTable)
  count = 0;
  return;
end
mask = abs(candidateTable.toothIdx) <= toothIdxTol;
if any(strcmp(candidateTable.Properties.VariableNames, 'toothResidualHz'))
  mask = mask & abs(candidateTable.toothResidualHz) <= residualTolHz;
end
count = nnz(mask);
end


function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return a scalar template for aggregate output rows.
row = struct('strategyName', "", 'bankSize', NaN, 'numRandomTrial', NaN, 'numRepeat', 0, ...
  'resolvedRate', NaN, 'truthToothIndexHitRate', NaN, 'nearToothIndexHitRate', NaN, ...
  'truthToothHitRate', NaN, 'nearToothHitRate', NaN, ...
  'sameToothResidualFailRate', NaN, 'nearToothResidualFailRate', NaN, ...
  'truthToothGainRate', NaN, 'nearToothGainRate', NaN, ...
  'truthToothIndexGainRate', NaN, 'nearToothIndexGainRate', NaN, ...
  'sameToothResidualFixRate', NaN, 'easyDamageRate', NaN, 'easyDamageCount', 0, ...
  'angleRmseDeg', NaN, 'angleP95Deg', NaN, 'angleMaxDeg', NaN, ...
  'fdRefRmseHz', NaN, 'meanRunTimeMs', NaN, 'candidateEvalCount', 0, ...
  'meanCandidateEvalPerRepeat', NaN, 'candidateTruthToothIndexCount', 0, ...
  'candidateNearToothIndexCount', 0, 'candidateTruthToothCount', 0, ...
  'candidateNearToothCount', 0, 'candidateTruthToothIndexRate', NaN, ...
  'candidateNearToothIndexRate', NaN, 'candidateTruthToothRate', NaN, ...
  'candidateNearToothRate', NaN, 'candidateAnyTruthIndexRate', NaN, ...
  'candidateAnyNearIndexRate', NaN, 'candidateAnyTruthToothRate', NaN, ...
  'candidateAnyNearToothRate', NaN, 'selectedMissDespiteTruthCandidateRate', NaN, ...
  'selectedMissDespiteNearCandidateRate', NaN, 'selectedLabelSummary', "", ...
  'candidateLabelSummary', "");
end

function txt = localBuildSelectedLabelSummary(labelVec)
%LOCALBUILDSELECTEDLABELSUMMARY Format label counts as a compact string.
labelList = unique(string(labelVec), 'stable');
parts = strings(numel(labelList), 1);
for iLabel = 1:numel(labelList)
  parts(iLabel) = labelList(iLabel) + ":" + string(nnz(string(labelVec) == labelList(iLabel)));
end
txt = strjoin(parts, ', ');
end


function representative = localSelectRepresentative(scanTable, candidateTable, repeatCellByStrategy)
%LOCALSELECTREPRESENTATIVE Pick one hard repeat and attach its candidate tables.
[~, idx] = max(abs(scanTable.toothIdx) + scanTable.angleErrDeg / max(1, max(scanTable.angleErrDeg, [], 'omitnan')));
representative = table2struct(scanTable(idx, :));
representative.candidateTable = table();
if ~isempty(candidateTable)
  mask = candidateTable.strategyName == representative.strategyName & candidateTable.taskSeed == representative.taskSeed;
  representative.candidateTable = candidateTable(mask, :);
end
strategyList = unique(scanTable.strategyName, 'stable');
iStrategy = find(strategyList == representative.strategyName, 1, 'first');
if ~isempty(iStrategy)
  repeatCell = repeatCellByStrategy{iStrategy};
  repIdx = find(cellfun(@(x) x.taskSeed == representative.taskSeed, repeatCell), 1, 'first');
  if ~isempty(repIdx)
    representative.periodicCandidateTable = repeatCell{repIdx}.msFlow.periodicCandidateTable;
    representative.periodicDoaSeed = repeatCell{repIdx}.msFlow.periodicDoaSeed;
  end
end
end


function displayTable = localSelectDisplayTable(scanTable)
%LOCALSELECTDISPLAYTABLE Select compact repeat-level columns for printing.
varList = {'strategyName', 'taskSeed', 'selectedSubsetLabel', 'truthToothIndexHit', ...
  'truthToothHit', 'nearToothHit', 'sameToothResidualFail', 'toothIdx', ...
  'toothResidualHz', 'transitionClass', 'angleErrDeg', 'fdRefErrHz', ...
  'fdRateErrHzPerSec', 'selectedFinalTag', 'usedPolish', ...
  'easyDamageFromBaseline', 'angleDeltaFromBaselineDeg', 'runTimeMs'};
displayTable = scanTable(:, varList(ismember(varList, scanTable.Properties.VariableNames)));
end


function toothHistogramTable = localBuildToothHistogramTable(scanTable, scanConfig)
%LOCALBUILDTOOTHHISTOGRAMTABLE Build a compact tooth-index distribution table.
if isempty(scanTable) || ~any(strcmp(scanTable.Properties.VariableNames, 'toothIdx'))
  toothHistogramTable = table();
  return;
end
strategyNameList = unique(scanTable.strategyName, 'stable');
rowCell = cell(numel(strategyNameList), 1);
for iStrategy = 1:numel(strategyNameList)
  name = strategyNameList(iStrategy);
  toothIdxVec = scanTable.toothIdx(scanTable.strategyName == name);
  toothIdxVec = toothIdxVec(~isnan(toothIdxVec));
  if isempty(toothIdxVec)
    rowCell{iStrategy} = table();
    continue;
  end
  edges = localBuildToothHistogramEdges(toothIdxVec, scanConfig.toothHistogramBinCount);
  count = histcounts(toothIdxVec, edges);
  count = count(:);
  binLeft = edges(1:end-1).';
  binRight = edges(2:end).';
  toothIdxCenter = (binLeft + binRight) / 2;
  binIdx = (1:numel(count)).';
  strategyName = repmat(name, numel(count), 1);
  rate = count / max(1, numel(toothIdxVec));
  rowCell{iStrategy} = table(strategyName, binIdx, binLeft, binRight, toothIdxCenter, count, rate);
end
toothHistogramTable = localVertcatTableCell(rowCell);
end


function edges = localBuildToothHistogramEdges(toothIdxVec, maxBinCount)
%LOCALBUILDTOOTHHISTOGRAMEDGES Build exact or coarsened tooth histogram bins.
lo = floor(min(toothIdxVec)) - 0.5;
hi = ceil(max(toothIdxVec)) + 0.5;
numExactBin = max(1, round(hi - lo));
if numExactBin <= maxBinCount
  edges = lo:hi;
else
  edges = linspace(lo, hi, maxBinCount + 1);
end
end


function plotData = localBuildPlotData(scanTable, aggregateTable, toothHistogramTable, candidateSeedCoverageTable, consensusAggregateTable)
%LOCALBUILDPLOTDATA Store lightweight table data needed to redraw figures.
plotData = struct();
plotData.scanTable = scanTable;
plotData.aggregateTable = aggregateTable;
plotData.toothHistogramTable = toothHistogramTable;
plotData.candidateSeedCoverageTable = candidateSeedCoverageTable;
plotData.consensusAggregateTable = consensusAggregateTable;
end


function plotData = localPlotScan(scanTable, aggregateTable, toothHistogramTable, candidateSeedCoverageTable, consensusAggregateTable)
%LOCALPLOTSCAN Draw scan figures without saving image files.
plotData = localBuildPlotData(scanTable, aggregateTable, toothHistogramTable, candidateSeedCoverageTable, consensusAggregateTable);
strategyCat = categorical(aggregateTable.strategyName);
strategyCat = reordercats(strategyCat, cellstr(aggregateTable.strategyName));

figure('Name', 'Subset-bank coverage vs candidate cost');
hold on;
xCost = aggregateTable.meanCandidateEvalPerRepeat;
if all(isnan(xCost))
  xCost = aggregateTable.bankSize;
end
plot(xCost, aggregateTable.truthToothHitRate, '-o', 'DisplayName', 'truth strict');
plot(xCost, aggregateTable.nearToothHitRate, '-s', 'DisplayName', 'near strict');
for iRow = 1:height(aggregateTable)
  text(xCost(iRow), aggregateTable.truthToothHitRate(iRow), " " + aggregateTable.strategyName(iRow));
end
hold off;
grid on; xlabel('mean candidate eval per repeat'); ylabel('selected strict hit rate');
title('Selected tooth coverage vs candidate cost');
legend('Location', 'best');

figure('Name', 'Subset-bank selected tooth decomposition');
bar(strategyCat, [aggregateTable.truthToothIndexHitRate, aggregateTable.truthToothHitRate, ...
  aggregateTable.nearToothHitRate, aggregateTable.sameToothResidualFailRate]);
grid on; xlabel('strategy'); ylabel('rate');
legend({'tooth index = 0', 'truth strict', 'near strict', 'same-tooth residual fail'}, 'Location', 'best');
title('Selected integer-tooth and residual-aware hits');
xtickangle(25);

figure('Name', 'Subset-bank gain and damage');
bar(strategyCat, [aggregateTable.truthToothGainRate, aggregateTable.nearToothGainRate, ...
  aggregateTable.sameToothResidualFixRate, aggregateTable.easyDamageRate]);
grid on; xlabel('strategy'); ylabel('rate vs curated12');
legend({'truth strict gain', 'near strict gain', 'same-tooth residual fix', 'easy damage'}, 'Location', 'best');
title('Baseline transition summary');
xtickangle(25);

figure('Name', 'Subset-bank angle and fd tail');
xIdx = 1:height(aggregateTable);
yyaxis left;
bar(xIdx, aggregateTable.angleP95Deg);
ylabel('angle p95 (deg)');
yyaxis right;
plot(xIdx, aggregateTable.fdRefRmseHz, '-o');
ylabel('fdRef RMSE (Hz)');
set(gca, 'XTick', xIdx, 'XTickLabel', cellstr(aggregateTable.strategyName));
grid on; xlabel('strategy'); title('Angle tail and fdRef tail');
xtickangle(25);

if ~isempty(candidateSeedCoverageTable)
  figure('Name', 'Subset-bank candidate coverage');
  bar(strategyCat, [aggregateTable.candidateAnyTruthIndexRate, aggregateTable.candidateAnyTruthToothRate, ...
    aggregateTable.selectedMissDespiteTruthCandidateRate]);
  grid on; xlabel('strategy'); ylabel('per-seed rate');
  legend({'any candidate tooth index = 0', 'any candidate truth strict', ...
    'selected missed truth candidate'}, 'Location', 'best');
  title('Candidate availability vs selected winner');
  xtickangle(25);
end

if ~isempty(consensusAggregateTable)
  consensusCat = categorical(consensusAggregateTable.strategyName);
  consensusCat = reordercats(consensusCat, cellstr(consensusAggregateTable.strategyName));
  figure('Name', 'No-truth consensus diagnostic');
  bar(consensusCat, [consensusAggregateTable.currentTruthToothHitRate, ...
    consensusAggregateTable.consensusTruthToothHitRate, consensusAggregateTable.currentNearToothHitRate, ...
    consensusAggregateTable.consensusNearToothHitRate]);
  grid on; xlabel('strategy'); ylabel('strict hit rate');
  legend({'current truth strict', 'consensus truth strict', ...
    'current near strict', 'consensus near strict'}, 'Location', 'best');
  title('No-truth anchor-group consensus diagnostic (soft residual penalty)');
  xtickangle(25);
end

if ~isempty(toothHistogramTable)
  figure('Name', 'Subset-bank tooth distribution');
  hold on;
  strategyNameList = unique(toothHistogramTable.strategyName, 'stable');
  for iStrategy = 1:numel(strategyNameList)
    name = strategyNameList(iStrategy);
    mask = toothHistogramTable.strategyName == name;
    plot(toothHistogramTable.toothIdxCenter(mask), toothHistogramTable.rate(mask), '-o', ...
      'DisplayName', char(name));
  end
  hold off;
  grid on; xlabel('tooth index'); ylabel('repeat rate');
  title('Selected tooth-index distribution');
  legend('Location', 'best');
end

drawnow;
end


function cleanupTable = localCleanupCheckpointRunStates(runStateCell, scanConfig)
%LOCALCLEANUPCHECKPOINTRUNSTATES Clean completed per-strategy checkpoint dirs.
cleanupTable = table();
if ~isfield(scanConfig, 'checkpointEnable') || ~scanConfig.checkpointEnable || ...
    ~isfield(scanConfig, 'checkpointCleanupOnSuccess') || ~scanConfig.checkpointCleanupOnSuccess
  return;
end
rowCell = {};
for iRun = 1:numel(runStateCell)
  runState = runStateCell{iRun};
  if isempty(runState) || ~isstruct(runState) || ...
      ~isfield(runState, 'runDir') || strlength(string(runState.runDir)) == 0
    continue;
  end
  if ~isfield(runState, 'isComplete') || ~runState.isComplete
    continue;
  end
  runDir = string(runState.runDir);
  cleanupPerfTaskGridCheckpoint(runState, struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
  row = struct();
  row.runName = string(runState.runName);
  row.runKey = string(runState.runKey);
  row.runDir = runDir;
  row.cleanedOnSuccess = ~isfolder(runDir);
  rowCell{end + 1, 1} = row; %#ok<AGROW>
end
if ~isempty(rowCell)
  cleanupTable = struct2table([rowCell{:}].');
end
end


function summaryTable = localBuildCheckpointSummaryTable(strategyList, runStateCell, checkpointDirCell, cleanupTable)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Summarize checkpoint resume and cleanup state.
rowList = repmat(localEmptyCheckpointSummaryRow(), numel(strategyList), 1);
for iStrategy = 1:numel(strategyList)
  row = localEmptyCheckpointSummaryRow();
  row.strategyName = strategyList(iStrategy).strategyName;
  row.expectedRunDir = string(checkpointDirCell(iStrategy));
  runState = runStateCell{iStrategy};
  if ~isempty(runState) && isstruct(runState)
    row.runName = string(localGetFieldOrDefault(runState, 'runName', ""));
    row.runKey = string(localGetFieldOrDefault(runState, 'runKey', ""));
    row.runDir = string(localGetFieldOrDefault(runState, 'runDir', row.expectedRunDir));
    row.numTask = localGetFieldOrDefault(runState, 'numTask', 0);
    row.numDone = localGetFieldOrDefault(runState, 'numDone', 0);
    row.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
  end
  if ~isempty(cleanupTable) && any(strcmp(cleanupTable.Properties.VariableNames, 'runDir'))
    cleanupIdx = find(string(cleanupTable.runDir) == row.runDir, 1, 'first');
    if ~isempty(cleanupIdx)
      row.cleanedOnSuccess = logical(cleanupTable.cleanedOnSuccess(cleanupIdx));
    end
  end
  rowList(iStrategy) = row;
end
summaryTable = struct2table(rowList);
end


function row = localEmptyCheckpointSummaryRow()
%LOCALEMPTYCHECKPOINTSUMMARYROW Return a scalar checkpoint summary template.
row = struct('strategyName', "", 'runName', "", 'runKey', "", 'runDir', "", ...
  'expectedRunDir', "", 'numTask', 0, 'numDone', 0, 'isComplete', false, ...
  'cleanedOnSuccess', false);
end


function localPrintCheckpointFailureHint(checkpointDirCell)
%LOCALPRINTCHECKPOINTFAILUREHINT Print preserved checkpoint dirs after failures.
checkpointDirCell = checkpointDirCell(strlength(string(checkpointDirCell)) > 0);
if isempty(checkpointDirCell)
  fprintf('Scan failed. No checkpoint directory was created.\n');
  return;
end
fprintf('Scan failed. Completed checkpoint tasks were kept under tmp for resume.\n');
for iDir = 1:numel(checkpointDirCell)
  fprintf('  checkpoint dir: %s\n', char(checkpointDirCell(iDir)));
end
end


function contextSummary = localBuildContextSummary(contextCell, flowOptCell, strategyList)
%LOCALBUILDCONTEXTSUMMARY Store lightweight shared context metadata.
contextSummary = struct();
if isempty(contextCell) || isempty(contextCell{1})
  return;
end
context = contextCell{1};
contextSummary.frameIntvlSec = context.frameIntvlSec;
contextSummary.periodicOffsetIdx = context.periodicOffsetIdx;
contextSummary.masterOffsetIdx = context.masterOffsetIdx;
contextSummary.selectedSatIdxGlobal = context.selectedSatIdxGlobal;
contextSummary.refSatIdxGlobal = context.refSatIdxGlobal;
contextSummary.flowSignatureTable = localBuildFlowSignatureTable(flowOptCell, strategyList);
end


function flowSignatureTable = localBuildFlowSignatureTable(flowOptCell, strategyList)
%LOCALBUILDFLOWSIGNATURETABLE Keep compact flow fields that affect bank behavior.
rowList = repmat(localEmptyFlowSignatureRow(), numel(flowOptCell), 1);
for iFlow = 1:numel(flowOptCell)
  row = localEmptyFlowSignatureRow();
  row.strategyName = strategyList(iFlow).strategyName;
  flowOpt = flowOptCell{iFlow};
  if isempty(flowOpt)
    rowList(iFlow) = row;
    continue;
  end
  row.subsetDefaultBankLabelList = localJoinStringList(localGetFieldOrDefault(flowOpt, ...
    'subsetDefaultBankLabelList', strings(1, 0)));
  row.subsetMarginFallbackBankLabelList = localJoinStringList(localGetFieldOrDefault(flowOpt, ...
    'subsetMarginFallbackBankLabelList', strings(1, 0)));
  row.subsetEnableMarginFallback = logical(localGetFieldOrDefault(flowOpt, 'subsetEnableMarginFallback', false));
  row.periodicRefineFdHalfWidthHz = localGetFieldOrDefault(flowOpt, 'periodicRefineFdHalfWidthHz', NaN);
  row.periodicRefineFdRateHalfWidthHzPerSec = localGetFieldOrDefault(flowOpt, 'periodicRefineFdRateHalfWidthHzPerSec', NaN);
  row.periodicRefineFreezeDoa = logical(localGetFieldOrDefault(flowOpt, 'periodicRefineFreezeDoa', false));
  row.periodicPolishEnableWhenMulti = logical(localGetFieldOrDefault(flowOpt, 'periodicPolishEnableWhenMulti', false));
  rowList(iFlow) = row;
end
flowSignatureTable = struct2table(rowList);
end


function row = localEmptyFlowSignatureRow()
%LOCALEMPTYFLOWSIGNATUREROW Return a scalar template for compact flow metadata.
row = struct('strategyName', "", 'subsetDefaultBankLabelList', "", ...
  'subsetMarginFallbackBankLabelList', "", 'subsetEnableMarginFallback', false, ...
  'periodicRefineFdHalfWidthHz', NaN, 'periodicRefineFdRateHalfWidthHzPerSec', NaN, ...
  'periodicRefineFreezeDoa', false, 'periodicPolishEnableWhenMulti', false);
end


function txt = localJoinStringList(value)
%LOCALJOINSTRINGLIST Format a string vector for compact metadata tables.
value = reshape(string(value), 1, []);
if isempty(value)
  txt = "";
else
  txt = strjoin(value, ', ');
end
end


function value = localGetTableScalar(t, varName, defaultValue)
%LOCALGETTABLESCALAR Read one scalar table variable with a fallback.
if istable(t) && height(t) >= 1 && any(strcmp(t.Properties.VariableNames, varName))
  value = t.(varName)(1);
else
  value = defaultValue;
end
end


function value = localGetFieldOrDefault(s, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a struct field with a scalar fallback.
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
  value = s.(fieldName);
else
  value = defaultValue;
end
end


function slimCellByStrategy = localStripRepeatCellByStrategy(repeatCellByStrategy)
%LOCALSTRIPREPEATCELLBYSTRATEGY Keep only lightweight repeat diagnostics.
slimCellByStrategy = cell(size(repeatCellByStrategy));
for iStrategy = 1:numel(repeatCellByStrategy)
  repeatCell = repeatCellByStrategy{iStrategy};
  slim = cell(size(repeatCell));
  for iCell = 1:numel(repeatCell)
    item = repeatCell{iCell};
    slim{iCell} = struct('taskSeed', item.taskSeed, 'summary', item.summary, ...
      'subsetCandidateTable', item.msFlow.subsetCandidateTable, ...
      'periodicCandidateTable', item.msFlow.periodicCandidateTable, ...
      'periodicDoaSeed', item.msFlow.periodicDoaSeed);
  end
  slimCellByStrategy{iStrategy} = slim;
end
end


function localLog(messageText)
%LOCALLOG Print one timestamped scan orchestration message.
fprintf('[%s] %s\n', localNowString('HH:mm:ss'), char(messageText));
end


function txt = localNowString(formatText)
%LOCALNOWSTRING Format the current local time for compact logs.
txt = char(datetime('now', 'Format', formatText));
end


function localPrintTablePreview(titleText, t, maxRows)
%LOCALPRINTTABLEPREVIEW Print a compact preview while keeping full tables in scanData.
fprintf('\n========== %s ==========%s', char(titleText), newline);
if isempty(t) || ~istable(t) || height(t) == 0
  fprintf('(empty)\n');
  return;
end
numShow = min(maxRows, height(t));
disp(t(1:numShow, :));
if height(t) > numShow
  fprintf('... %d more rows saved in scanData.\n', height(t) - numShow);
end
end

