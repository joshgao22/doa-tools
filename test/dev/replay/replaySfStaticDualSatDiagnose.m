% replaySfStaticDualSatDiagnose
% Purpose: replay fixed-SNR single-frame static dual-satellite DoA-Doppler
% diagnostics on the canonical 4154+1165 pair. This replay is the
% mechanism-level companion of scanSfStaticMleCrbConsistency: it runs a
% small fixed-SNR repeat batch, keeps the SS/MS static transition path shared
% through common fixture helpers, and inspects weight-sweep, representative
% repeat, per-satellite, and reference-satellite invariant diagnostics.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted repeat runs from per-task files under the repo-root tmp. The
% outer repeat loop uses parfor automatically when available. saveSnapshot=true
% saves only replayData via saveExpSnapshot. Telegram notice is best-effort only.

clear; close all; clc;

%% Replay configuration

replayName = "replaySfStaticDualSatDiagnose";
saveSnapshot = false;
checkpointEnable = false;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = true;

snrDb = 10;
baseSeed = 253;
numRepeat = 40;
weightSweepAlpha = [0; 0.25; 0.5; 1];
staticMsHalfWidthDeg = [0.002; 0.002];
topWorstCount = 10;
representativePreviewRows = 6;

% buildSfStaticSingleSnrRepeatFixture currently defines the canonical static
% replay seed chain as baseSeed=253. Keep this replay as a thin wrapper over
% that shared fixture instead of duplicating static scene construction here.
if baseSeed ~= 253
  error('replaySfStaticDualSatDiagnose:UnsupportedBaseSeed', ...
    'This replay reuses buildSfStaticSingleSnrRepeatFixture, whose canonical baseSeed is 253.');
end

weightSweepAlpha = reshape(double(weightSweepAlpha), [], 1);
staticMsHalfWidthDeg = reshape(double(staticMsHalfWidthDeg), [], 1);
seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);

runTic = tic;
replayData = struct();
config = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and flow options

  config.replayName = string(replayName);
  config.snrDb = snrDb;
  config.baseSeed = baseSeed;
  config.numRepeat = numRepeat;
  config.seedList = seedList;
  config.weightSweepAlpha = weightSweepAlpha;
  config.staticMsHalfWidthDeg = staticMsHalfWidthDeg;
  config.topWorstCount = topWorstCount;
  config.representativePreviewRows = representativePreviewRows;
  config.saveSnapshot = logical(saveSnapshot);
  config.checkpointEnable = logical(checkpointEnable);
  config.checkpointResume = logical(checkpointResume);
  config.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
  config.notifyTelegramEnable = logical(notifyTelegramEnable);
  config.numTask = numRepeat;
  config.useParfor = localCanUseParfor(numRepeat);
  config.runKey = localBuildRunKey(config);

  taskList = localBuildTaskList(config.numRepeat, config.seedList);
  checkpointOpt = buildMfReplayCheckpointOpt(replayName, config, struct( ...
    'runName', replayName, ...
    'runKey', config.runKey, ...
    'useParfor', config.useParfor, ...
    'meta', localBuildCheckpointMeta(config)));
  if config.checkpointEnable
    checkpointRunDir = string(checkpointOpt.runDir);
    config.checkpointRunDir = checkpointRunDir;
  end

  printMfReplayHeader(char(replayName), config, char(checkpointRunDir));
  localPrintReplayConfig(config);

  %% Run replay batch

  try
    [repeatCell, runState] = localRunTaskBatch(config, taskList, checkpointOpt);
  catch ME
    if strlength(string(checkpointRunDir)) > 0
      fprintf('Replay failed. Checkpoint artifacts kept at: %s\n', char(checkpointRunDir));
    end
    rethrow(ME);
  end

  caseTable = localCollectCaseTable(repeatCell);
  weightRepeatTable = localCollectWeightTable(repeatCell);
  repeatCompareTable = localCollectRepeatCompareTable(repeatCell);
  caseSummaryTable = localBuildCaseSummaryTable(caseTable);
  staticGainSummaryTable = localBuildStaticGainSummaryTable(repeatCompareTable);
  weightPreferenceTable = localBuildWeightPreferenceTable(weightRepeatTable, repeatCompareTable, config.weightSweepAlpha);
  worstMsAngleLossTable = localBuildWorstMsAngleLossTable(repeatCompareTable, config.topWorstCount);
  representativeTable = localSelectRepresentativeRepeats(repeatCompareTable);
  representativeDetail = localBuildRepresentativeDetail(representativeTable, config);
  staticCrbSummaryTable = localBuildStaticCrbSummary(config);
  firstFixture = buildSfStaticSingleSnrRepeatFixture(1, config.snrDb);
  selectedSatelliteTable = localBuildSelectedSatelliteTable(firstFixture.truth);
  truthTable = localBuildTruthTable(firstFixture.truth);
  referenceInvariantTable = localBuildReferenceInvariantTable( ...
    [firstFixture.caseBundle.caseStaticRefOnly, firstFixture.caseBundle.caseStaticMs, firstFixture.caseBundle.weightCase], ...
    firstFixture.truth);
  plotData = localBuildPlotData(repeatCompareTable, weightRepeatTable, caseSummaryTable);

  %% Data storage

  replayData = struct();
  replayData.replayName = string(replayName);
  replayData.runKey = config.runKey;
  replayData.utcRun = datetime('now', 'TimeZone', 'local');
  replayData.config = config;
  replayData.truth = firstFixture.truth;
  replayData.selectedSatelliteTable = selectedSatelliteTable;
  replayData.truthTable = truthTable;
  replayData.caseTable = caseTable;
  replayData.weightRepeatTable = weightRepeatTable;
  replayData.repeatCompareTable = repeatCompareTable;
  replayData.caseSummaryTable = caseSummaryTable;
  replayData.staticGainSummaryTable = staticGainSummaryTable;
  replayData.weightPreferenceTable = weightPreferenceTable;
  replayData.worstMsAngleLossTable = worstMsAngleLossTable;
  replayData.representativeTable = representativeTable;
  replayData.representativeDetail = representativeDetail;
  replayData.staticCrbSummaryTable = staticCrbSummaryTable;
  replayData.referenceInvariantTable = referenceInvariantTable;
  if config.checkpointEnable
    replayData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  end
  replayData.plotData = plotData;
  replayData.elapsedSec = toc(runTic);
  replayData = finalizeMfReplayResult(replayData, "");

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'replayData'}}, ...
      'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
    replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
  else
    replayData.snapshotFile = "";
  end

  if config.checkpointEnable
    if config.checkpointCleanupOnSuccess
      replayData.checkpointSummary.cleanupReport = cleanupPerfTaskGridCheckpoint(runState, struct('verbose', false));
      replayData.checkpointSummary.cleanedOnSuccess = true;
    else
      replayData.checkpointSummary.cleanedOnSuccess = false;
    end
  end

  %% Summary output and plotting

  printMfReplaySection('Selected satellites', replayData.selectedSatelliteTable);
  printMfReplaySection('Truth and reference satellite', replayData.truthTable);
  printMfReplaySection('Single-SNR static case summary', replayData.caseSummaryTable);
  printMfReplaySection('MS vs SS static gain summary', replayData.staticGainSummaryTable);
  printMfReplaySection('Weight preference summary', replayData.weightPreferenceTable);
  printMfReplaySection('Worst MS angle-loss repeats', replayData.worstMsAngleLossTable);
  printMfReplaySection('Static CRB summary', replayData.staticCrbSummaryTable);
  printMfReplaySection('Reference-satellite invariant check', replayData.referenceInvariantTable);
  printMfReplaySection('Representative repeat index', replayData.representativeTable);

  for iRep = 1:numel(replayData.representativeDetail)
    repInfo = replayData.representativeDetail(iRep);
    fprintf('\n========== Representative repeat: %s (seed %d) ==========%s', ...
      char(repInfo.label), repInfo.taskSeed, newline);
    fprintf('---------- Main static cases ----------%s', newline);
    disp(repInfo.caseTable);
    fprintf('---------- Static ablation ----------%s', newline);
    disp(repInfo.ablationTable);
    fprintf('---------- Weight-sweep solve-path diagnostics ----------%s', newline);
    disp(repInfo.weightPathDiagTable);
    fprintf('---------- MS-SF-Static per-sat diagnostics ----------%s', newline);
    disp(repInfo.msSatDiagTable);
    if ~isempty(repInfo.bestWeightSatDiagTable) && height(repInfo.bestWeightSatDiagTable) > 0
      fprintf('---------- Best-weight MS-SF-Static per-sat diagnostics ----------%s', newline);
      disp(repInfo.bestWeightSatDiagTable);
    end
    fprintf('---------- Representative reference invariants ----------%s', newline);
    disp(repInfo.referenceInvariantTable);
  end

  if height(replayData.caseTable) > 2 * config.representativePreviewRows
    fprintf('\nCase table preview:\n');
    dispMfReplayTablePreview(replayData.caseTable, config.representativePreviewRows);
  end
  if height(replayData.weightRepeatTable) > 2 * config.representativePreviewRows
    fprintf('\nWeight repeat table preview:\n');
    dispMfReplayTablePreview(replayData.weightRepeatTable, config.representativePreviewRows);
  end

  replayData.plotData = localPlotReplay(replayData);

  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "DONE", ...
    'subtitleText', "Single-frame static dual-satellite replay diagnostics", ...
    'config', replayData.config, ...
    'snapshotFile', replayData.snapshotFile, ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', replayData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(replayData), ...
    'commentLineList', [ ...
      "Static replay completed through the shared SF static fixture path."; ...
      "Use scanSfStaticMleCrbConsistency for paper-facing MC curves and this replay for representative weight/per-sat diagnostics."]));

catch ME
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'subtitleText', "Single-frame static dual-satellite replay diagnostics", ...
    'config', config, ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function taskList = localBuildTaskList(numRepeat, seedList)
%LOCALBUILDTASKLIST Build one static repeat task per seed.

taskList = repmat(struct('taskId', 0, 'repeatIdx', 0, 'taskSeed', 0), numRepeat, 1);
for iRepeat = 1:numRepeat
  taskList(iRepeat).taskId = iRepeat;
  taskList(iRepeat).repeatIdx = iRepeat;
  taskList(iRepeat).taskSeed = seedList(iRepeat);
end
end

function [repeatCell, runState] = localRunTaskBatch(config, taskList, checkpointOpt)
%LOCALRUNTASKBATCH Run the static repeat batch with optional checkpointing.

sharedData = struct('config', config);
runState = struct();
if config.checkpointEnable
  runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt);
  runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, @localRunCheckpointTask, runnerOpt);
  repeatCell = runState.resultCell;
  return;
end

repeatCell = cell(numel(taskList), 1);
if config.useParfor
  parfor iTask = 1:numel(taskList)
    repeatCell{iTask} = localRunOneTask(taskList(iTask), sharedData);
  end
else
  progressbar('reset', numel(taskList));
  cleanupObj = onCleanup(@() progressbar('end')); %#ok<NASGU>
  for iTask = 1:numel(taskList)
    repeatCell{iTask} = localRunOneTask(taskList(iTask), sharedData);
    progressbar('advance');
  end
  clear cleanupObj;
end
end

function repeatOut = localRunCheckpointTask(taskInfo, sharedData)
%LOCALRUNCHECKPOINTTASK Run one checkpointed static repeat task.

repeatOut = localRunOneTask(taskInfo, sharedData);
end

function repeatOut = localRunOneTask(taskInfo, sharedData)
%LOCALRUNONETASK Run one static repeat using the shared canonical fixture.

fixture = buildSfStaticSingleSnrRepeatFixture(taskInfo.repeatIdx, sharedData.config.snrDb);
repeatOut = localPackRepeatOutput(fixture, taskInfo, sharedData.config);
end

function repeatOut = localPackRepeatOutput(fixture, taskInfo, config)
%LOCALPACKREPEATOUTPUT Strip one fixture into lightweight replay tables.

caseTable = localAddRepeatColumns(fixture.caseTable, taskInfo, config.snrDb);
weightTable = buildDoaDopplerWeightSweepTable( ...
  fixture.weightSweepAlpha, fixture.caseBundle.weightCase, fixture.truth);
weightTable = localAddRepeatColumns(weightTable, taskInfo, config.snrDb);
compareRow = localBuildRepeatCompareRow(taskInfo, config.snrDb, fixture.caseTable, weightTable);

repeatOut = struct();
repeatOut.repeatIdx = taskInfo.repeatIdx;
repeatOut.taskSeed = taskInfo.taskSeed;
repeatOut.snrDb = config.snrDb;
repeatOut.caseTable = caseTable;
repeatOut.weightSummaryTable = weightTable;
repeatOut.compareRow = compareRow;
end

function dataTable = localAddRepeatColumns(dataTable, taskInfo, snrDb)
%LOCALADDREPEATCOLUMNS Add repeat metadata to a compact table.

if isempty(dataTable) || ~istable(dataTable)
  return;
end
numRow = height(dataTable);
repeatIdx = repmat(taskInfo.repeatIdx, numRow, 1);
taskSeed = repmat(taskInfo.taskSeed, numRow, 1);
snrDbCol = repmat(snrDb, numRow, 1);
dataTable = addvars(dataTable, repeatIdx, taskSeed, snrDbCol, 'Before', 1, ...
  'NewVariableNames', {'repeatIdx', 'taskSeed', 'snrDb'});
if ismember('snrDb_1', dataTable.Properties.VariableNames)
  dataTable.snrDb_1 = [];
end
end

function compareRow = localBuildRepeatCompareRow(taskInfo, snrDb, caseTable, weightTable)
%LOCALBUILDREPEATCOMPAREROW Build one repeat-level SS/MS compare row.

ssRow = localSelectCaseRow(caseTable, "SS-SF-Static");
msRow = localSelectCaseRow(caseTable, "MS-SF-Static");
[bestWeightRow, bestAlpha] = localSelectBestWeightRow(weightTable);

ssAngleErrDeg = localTableValue(ssRow, 'angleErrDeg', NaN);
msAngleErrDeg = localTableValue(msRow, 'angleErrDeg', NaN);
ssFdErrHz = abs(localTableValue(ssRow, 'fdRefErrHz', NaN));
msFdErrHz = abs(localTableValue(msRow, 'fdRefErrHz', NaN));
bestWeightAngleErrDeg = localTableValue(bestWeightRow, 'angleErrDeg', NaN);
bestWeightFdErrHz = abs(localTableValue(bestWeightRow, 'fdRefErrHz', NaN));

compareRow = table(taskInfo.repeatIdx, taskInfo.taskSeed, snrDb, ...
  ssAngleErrDeg, msAngleErrDeg, msAngleErrDeg - ssAngleErrDeg, ...
  ssFdErrHz, msFdErrHz, msFdErrHz - ssFdErrHz, ...
  logical(localTableValue(ssRow, 'isResolved', false)), ...
  logical(localTableValue(msRow, 'isResolved', false)), ...
  bestAlpha, bestWeightAngleErrDeg, bestWeightFdErrHz, ...
  bestWeightAngleErrDeg - ssAngleErrDeg, bestWeightAngleErrDeg - msAngleErrDeg, ...
  'VariableNames', {'repeatIdx', 'taskSeed', 'snrDb', ...
  'ssAngleErrDeg', 'msAngleErrDeg', 'msMinusSsAngleDeg', ...
  'ssFdErrHz', 'msFdErrHz', 'msMinusSsFdHz', 'ssResolved', 'msResolved', ...
  'bestAlphaSat2', 'bestWeightAngleErrDeg', 'bestWeightFdErrHz', ...
  'bestWeightMinusSsAngleDeg', 'bestWeightMinusMsAngleDeg'});
end

function caseRow = localSelectCaseRow(caseTable, displayName)
%LOCALSELECTCASEROW Select one displayName row from a case table.

caseRow = table();
if isempty(caseTable) || ~istable(caseTable) || ~ismember('displayName', caseTable.Properties.VariableNames)
  return;
end
matchIdx = find(string(caseTable.displayName) == string(displayName), 1, 'first');
if ~isempty(matchIdx)
  caseRow = caseTable(matchIdx, :);
end
end

function [bestRow, bestAlpha] = localSelectBestWeightRow(weightTable)
%LOCALSELECTBESTWEIGHTROW Select the lowest-angle resolved weight row.

bestRow = table();
bestAlpha = NaN;
if isempty(weightTable) || ~istable(weightTable) || height(weightTable) == 0
  return;
end
mask = true(height(weightTable), 1);
if ismember('isResolved', weightTable.Properties.VariableNames)
  mask = mask & logical(weightTable.isResolved);
end
if ismember('angleErrDeg', weightTable.Properties.VariableNames)
  mask = mask & isfinite(weightTable.angleErrDeg);
else
  mask(:) = false;
end
if ~any(mask)
  return;
end
validIdx = find(mask);
[~, bestLocal] = min(weightTable.angleErrDeg(mask));
bestIdx = validIdx(bestLocal);
bestRow = weightTable(bestIdx, :);
bestAlpha = localTableValue(bestRow, 'alphaSat2', NaN);
end

function value = localTableValue(dataTable, varName, defaultValue)
%LOCALTABLEVALUE Read the first scalar table value with a default.

value = defaultValue;
if isempty(dataTable) || ~istable(dataTable) || ~ismember(varName, dataTable.Properties.VariableNames) || height(dataTable) == 0
  return;
end
rawValue = dataTable.(varName);
if isempty(rawValue)
  return;
end
if iscell(rawValue)
  rawValue = rawValue{1};
else
  rawValue = rawValue(1);
end
if isnumeric(rawValue) || islogical(rawValue)
  value = rawValue;
else
  value = string(rawValue);
end
end

function caseTable = localCollectCaseTable(repeatCell)
%LOCALCOLLECTCASETABLE Collect compact per-case rows from all repeats.

caseTable = table();
for iRepeat = 1:numel(repeatCell)
  if isempty(repeatCell{iRepeat}) || ~isfield(repeatCell{iRepeat}, 'caseTable')
    continue;
  end
  caseTable = [caseTable; repeatCell{iRepeat}.caseTable]; %#ok<AGROW>
end
end

function weightTable = localCollectWeightTable(repeatCell)
%LOCALCOLLECTWEIGHTTABLE Collect compact per-weight rows from all repeats.

weightTable = table();
for iRepeat = 1:numel(repeatCell)
  if isempty(repeatCell{iRepeat}) || ~isfield(repeatCell{iRepeat}, 'weightSummaryTable')
    continue;
  end
  weightTable = [weightTable; repeatCell{iRepeat}.weightSummaryTable]; %#ok<AGROW>
end
end

function repeatTable = localCollectRepeatCompareTable(repeatCell)
%LOCALCOLLECTREPEATCOMPARETABLE Collect one compare row per repeat.

repeatTable = table();
for iRepeat = 1:numel(repeatCell)
  if isempty(repeatCell{iRepeat}) || ~isfield(repeatCell{iRepeat}, 'compareRow')
    continue;
  end
  repeatTable = [repeatTable; repeatCell{iRepeat}.compareRow]; %#ok<AGROW>
end
end

function summaryTable = localBuildCaseSummaryTable(caseTable)
%LOCALBUILDCASESUMMARYTABLE Build a single-SNR case aggregate summary.

if isempty(caseTable) || height(caseTable) == 0
  summaryTable = table();
  return;
end
methodNameList = unique(string(caseTable.displayName), 'stable');
numCase = numel(methodNameList);
snrDb = nan(numCase, 1);
numRepeat = zeros(numCase, 1);
angleRmseDeg = nan(numCase, 1);
angleP95Deg = nan(numCase, 1);
fdRmseHz = nan(numCase, 1);
fdP95Hz = nan(numCase, 1);
resolveRate = nan(numCase, 1);

for iCase = 1:numCase
  mask = string(caseTable.displayName) == methodNameList(iCase);
  subTable = caseTable(mask, :);
  numRepeat(iCase) = height(subTable);
  snrDb(iCase) = median(subTable.snrDb, 'omitnan');
  angleFail = ~logical(subTable.isResolved) | ~isfinite(subTable.angleErrDeg);
  angleStat = summarizeMonteCarloStat(subTable.angleErrDeg, angleFail);
  angleRmseDeg(iCase) = angleStat.rmse;
  angleP95Deg(iCase) = angleStat.p95;
  resolveRate(iCase) = 1 - angleStat.failRate;
  if ismember('fdRefErrHz', subTable.Properties.VariableNames) && any(isfinite(subTable.fdRefErrHz))
    fdErrAbs = abs(subTable.fdRefErrHz);
    fdFail = ~logical(subTable.isResolved) | ~isfinite(fdErrAbs);
    fdStat = summarizeMonteCarloStat(fdErrAbs, fdFail);
    fdRmseHz(iCase) = fdStat.rmse;
    fdP95Hz(iCase) = fdStat.p95;
  end
end

summaryTable = table(methodNameList(:), snrDb, numRepeat, angleRmseDeg, angleP95Deg, ...
  fdRmseHz, fdP95Hz, resolveRate, ...
  'VariableNames', {'displayName', 'snrDb', 'numRepeat', 'angleRmseDeg', ...
  'angleP95Deg', 'fdRmseHz', 'fdP95Hz', 'resolveRate'});
end

function gainTable = localBuildStaticGainSummaryTable(repeatTable)
%LOCALBUILDSTATICGAINSUMMARYTABLE Build one compact SS/MS static gain table.

if isempty(repeatTable) || height(repeatTable) == 0
  gainTable = table();
  return;
end
maskStatic = repeatTable.ssResolved & repeatTable.msResolved & ...
  isfinite(repeatTable.ssAngleErrDeg) & isfinite(repeatTable.msAngleErrDeg);
maskFd = maskStatic & isfinite(repeatTable.ssFdErrHz) & isfinite(repeatTable.msFdErrHz);
maskBest = maskStatic & isfinite(repeatTable.bestWeightAngleErrDeg);

metricName = [ ...
  "numRepeat"; ...
  "msBetterAngleRate"; ...
  "msBetterFdRate"; ...
  "medianMsMinusSsAngleDeg"; ...
  "medianMsMinusSsFdHz"; ...
  "bestWeightBetterThanSsRate"; ...
  "bestWeightBetterThanMsRate"; ...
  "medianBestWeightMinusSsAngleDeg"; ...
  "medianBestWeightMinusMsAngleDeg"];
metricValue = [ ...
  height(repeatTable); ...
  localMeanLogical(repeatTable.msMinusSsAngleDeg(maskStatic) < 0); ...
  localMeanLogical(repeatTable.msMinusSsFdHz(maskFd) < 0); ...
  median(repeatTable.msMinusSsAngleDeg(maskStatic), 'omitnan'); ...
  median(repeatTable.msMinusSsFdHz(maskFd), 'omitnan'); ...
  localMeanLogical(repeatTable.bestWeightMinusSsAngleDeg(maskBest) < 0); ...
  localMeanLogical(repeatTable.bestWeightMinusMsAngleDeg(maskBest) < 0); ...
  median(repeatTable.bestWeightMinusSsAngleDeg(maskBest), 'omitnan'); ...
  median(repeatTable.bestWeightMinusMsAngleDeg(maskBest), 'omitnan')];

gainTable = table(metricName, metricValue, 'VariableNames', {'metricName', 'metricValue'});
end

function weightPrefTable = localBuildWeightPreferenceTable(weightTable, repeatTable, weightAlpha)
%LOCALBUILDWEIGHTPREFERENCETABLE Build compact weight preference summary.

numWeight = numel(weightAlpha);
numBestRepeat = zeros(numWeight, 1);
angleRmseDeg = nan(numWeight, 1);
fdRmseHz = nan(numWeight, 1);
resolveRate = nan(numWeight, 1);

for iWeight = 1:numWeight
  alpha = weightAlpha(iWeight);
  numBestRepeat(iWeight) = sum(abs(repeatTable.bestAlphaSat2 - alpha) < 1e-12);
  mask = abs(weightTable.alphaSat2 - alpha) < 1e-12;
  subTable = weightTable(mask, :);
  if isempty(subTable) || height(subTable) == 0
    continue;
  end
  angleFail = ~logical(subTable.isResolved) | ~isfinite(subTable.angleErrDeg);
  angleStat = summarizeMonteCarloStat(subTable.angleErrDeg, angleFail);
  angleRmseDeg(iWeight) = angleStat.rmse;
  resolveRate(iWeight) = 1 - angleStat.failRate;
  fdErrAbs = abs(subTable.fdRefErrHz);
  fdFail = ~logical(subTable.isResolved) | ~isfinite(fdErrAbs);
  fdStat = summarizeMonteCarloStat(fdErrAbs, fdFail);
  fdRmseHz(iWeight) = fdStat.rmse;
end

weightPrefTable = table(weightAlpha(:), numBestRepeat, angleRmseDeg, fdRmseHz, resolveRate, ...
  'VariableNames', {'alphaSat2', 'numBestRepeat', 'angleRmseDeg', 'fdRmseHz', 'resolveRate'});
end

function worstTable = localBuildWorstMsAngleLossTable(repeatTable, topWorstCount)
%LOCALBUILDWORSTMSANGLELOSSTABLE Select the largest MS-vs-SS angle-loss repeats.

if isempty(repeatTable) || height(repeatTable) == 0
  worstTable = table();
  return;
end
mask = isfinite(repeatTable.msMinusSsAngleDeg);
worstTable = repeatTable(mask, :);
worstTable = sortrows(worstTable, 'msMinusSsAngleDeg', 'descend');
worstTable = worstTable(1:min(topWorstCount, height(worstTable)), :);
end

function repTable = localSelectRepresentativeRepeats(repeatTable)
%LOCALSELECTREPRESENTATIVEREPEATS Select a few representative repeats.

labelList = strings(0, 1);
idxList = zeros(0, 1);

[bestGainIdx, hasBestGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'min');
if hasBestGain
  labelList(end + 1, 1) = "bestMsAngleGain";
  idxList(end + 1, 1) = repeatTable.repeatIdx(bestGainIdx);
end

[worstGainIdx, hasWorstGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'max');
if hasWorstGain
  labelList(end + 1, 1) = "worstMsAngleGain";
  idxList(end + 1, 1) = repeatTable.repeatIdx(worstGainIdx);
end

[medianMsIdx, hasMedianMs] = localSelectMedianResolved(repeatTable.msAngleErrDeg, repeatTable.msResolved);
if hasMedianMs
  labelList(end + 1, 1) = "medianMsStatic";
  idxList(end + 1, 1) = repeatTable.repeatIdx(medianMsIdx);
end

[bestWeightLiftIdx, hasWeightLift] = localSelectByMetric(repeatTable.bestWeightMinusMsAngleDeg, 'min');
if hasWeightLift
  labelList(end + 1, 1) = "bestWeightLift";
  idxList(end + 1, 1) = repeatTable.repeatIdx(bestWeightLiftIdx);
end

[idxUnique, keepIdx] = unique(idxList, 'stable');
taskSeed = nan(numel(idxUnique), 1);
for iRow = 1:numel(idxUnique)
  matchIdx = find(repeatTable.repeatIdx == idxUnique(iRow), 1, 'first');
  if ~isempty(matchIdx)
    taskSeed(iRow) = repeatTable.taskSeed(matchIdx);
  end
end
repTable = table(labelList(keepIdx), idxUnique, taskSeed, ...
  'VariableNames', {'label', 'repeatIdx', 'taskSeed'});
end

function repDetail = localBuildRepresentativeDetail(repTable, config)
%LOCALBUILDREPRESENTATIVEDETAIL Rerun representative repeats and keep diagnostics.

repTemplate = struct('label', "", 'repeatIdx', NaN, 'taskSeed', NaN, ...
  'caseTable', table(), 'ablationTable', table(), 'weightSummaryTable', table(), ...
  'weightPathDiagTable', table(), 'msSatDiagTable', table(), ...
  'bestWeightSatDiagTable', table(), 'referenceInvariantTable', table());
repDetail = repmat(repTemplate, height(repTable), 1);
for iRep = 1:height(repTable)
  repeatIdx = repTable.repeatIdx(iRep);
  fixture = buildSfStaticSingleSnrRepeatFixture(repeatIdx, config.snrDb);
  weightSummaryTable = buildDoaDopplerWeightSweepTable( ...
    fixture.weightSweepAlpha, fixture.caseBundle.weightCase, fixture.truth);
  [bestWeightCase, ~] = localSelectBestWeightCase(fixture.caseBundle.weightCase, fixture.truth);

  repDetail(iRep).label = string(repTable.label(iRep));
  repDetail(iRep).repeatIdx = repeatIdx;
  repDetail(iRep).taskSeed = repTable.taskSeed(iRep);
  repDetail(iRep).caseTable = fixture.caseTable;
  repDetail(iRep).ablationTable = fixture.ablationTable;
  repDetail(iRep).weightSummaryTable = weightSummaryTable;
  repDetail(iRep).weightPathDiagTable = buildDynamicWeightSolvePathDiagTable( ...
    fixture.caseBundle.caseMsDoa, fixture.caseBundle.caseStaticRefOnly, ...
    fixture.caseBundle.weightCase, fixture.weightSweepAlpha, fixture.truth);
  repDetail(iRep).msSatDiagTable = localBuildStaticSatDiagTable( ...
    fixture.caseBundle.caseStaticMs.estResult, fixture.truth, fixture.scene);
  if ~isempty(bestWeightCase)
    repDetail(iRep).bestWeightSatDiagTable = localBuildStaticSatDiagTable( ...
      bestWeightCase.estResult, fixture.truth, fixture.scene);
  end
  repDetail(iRep).referenceInvariantTable = localBuildReferenceInvariantTable( ...
    [fixture.caseBundle.caseStaticRefOnly, fixture.caseBundle.caseStaticMs, fixture.caseBundle.weightCase], ...
    fixture.truth);
end
end

function [bestCase, bestIdx] = localSelectBestWeightCase(weightCase, truth)
%LOCALSELECTBESTWEIGHTCASE Select the lowest-angle resolved weight case.

bestCase = [];
bestIdx = NaN;
bestAngle = inf;
truthLatlon = truth.latlonTrueDeg(:);
for iCase = 1:numel(weightCase)
  caseInfo = weightCase(iCase);
  if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
    continue;
  end
  estResult = caseInfo.estResult;
  isResolved = getDoaDopplerFieldOrDefault(estResult, 'isResolved', true);
  if ~isResolved
    continue;
  end
  latlonEst = getDoaDopplerLatlonEst(estResult);
  if any(~isfinite(latlonEst))
    continue;
  end
  angleErrDeg = calcLatlonAngleError(latlonEst(:), truthLatlon);
  if angleErrDeg < bestAngle
    bestAngle = angleErrDeg;
    bestCase = caseInfo;
    bestIdx = iCase;
  end
end
end

function crbTable = localBuildStaticCrbSummary(config)
%LOCALBUILDSTATICCRBSUMMARY Build one single-SNR static CRB anchor table.

fixture = buildSfStaticSingleSnrRepeatFixture(1, config.snrDb);
crbOpt = struct('doaType', 'latlon');
[crbSs, auxCrbSs] = crbPilotSfDoaDoppler( ...
  fixture.sceneRefOnly, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdRefTrueHz, 1, fixture.pwrNoise, crbOpt);
[crbMs, auxCrbMs] = crbPilotSfDoaDoppler( ...
  fixture.scene, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.truth.latlonTrueDeg, fixture.truth.fdRefTrueHz, 1, fixture.pwrNoise, crbOpt);

caseName = ["SS-SF-Static"; "MS-SF-Static"];
satMode = ["single"; "multi"];
angleCrbStdDeg = [ ...
  projectCrbToAngleMetric(crbSs(1:2, 1:2), fixture.truth.latlonTrueDeg, 'latlon'); ...
  projectCrbToAngleMetric(crbMs(1:2, 1:2), fixture.truth.latlonTrueDeg, 'latlon')];
fdRefCrbStdHz = [sqrt(max(real(crbSs(3, 3)), 0)); sqrt(max(real(crbMs(3, 3)), 0))];
fdRateMode = [string(auxCrbSs.fdRateMode); string(auxCrbMs.fdRateMode)];
crbTable = table(repmat(config.snrDb, 2, 1), caseName, satMode, angleCrbStdDeg, fdRefCrbStdHz, fdRateMode, ...
  'VariableNames', {'snrDb', 'displayName', 'satMode', 'angleCrbStdDeg', 'fdRefCrbStdHz', 'fdRateMode'});
end

function satTable = localBuildSelectedSatelliteTable(truth)
%LOCALBUILDSELECTEDSATELLITETABLE Build compact selected-satellite truth table.

numSat = numel(truth.selectedSatIdxGlobal);
localSatIdx = (1:numSat).';
satTable = table(localSatIdx, truth.selectedSatIdxGlobal(:), truth.usrElevationDeg(:), ...
  truth.fdSatTrueHz(:), truth.deltaFdTrueHz(:), localSatIdx == truth.refSatIdxLocal, ...
  'VariableNames', {'localSatIdx', 'globalSatIdx', 'usrElevationDeg', ...
  'fdSatTrueHz', 'deltaFdTrueHz', 'isRefSat'});
end

function truthTable = localBuildTruthTable(truth)
%LOCALBUILDTRUTHTABLE Build compact truth/reference table.

truthTable = table(truth.latlonTrueDeg(1), truth.latlonTrueDeg(2), truth.fdRefTrueHz, ...
  truth.refSatIdxLocal, truth.refSatIdxGlobal, string(truth.refStateSource), ...
  'VariableNames', {'latTrueDeg', 'lonTrueDeg', 'fdRefTrueHz', ...
  'refSatIdxLocal', 'refSatIdxGlobal', 'refStateSource'});
end

function invariantTable = localBuildReferenceInvariantTable(caseList, truth)
%LOCALBUILDREFERENCEINVARIANTTABLE Check fdRef / deltaFd / fdSat closure.

numCase = numel(caseList) + 1;
displayName = strings(numCase, 1);
deltaFdAtRefHz = nan(numCase, 1);
fdSatRefMinusFdRefHz = nan(numCase, 1);
maxFdClosureAbsHz = nan(numCase, 1);

refIdx = truth.refSatIdxLocal;
displayName(1) = "truth";
deltaFdAtRefHz(1) = truth.deltaFdTrueHz(refIdx);
fdSatRefMinusFdRefHz(1) = truth.fdSatTrueHz(refIdx) - truth.fdRefTrueHz;
maxFdClosureAbsHz(1) = max(abs(truth.fdSatTrueHz(:) - (truth.fdRefTrueHz + truth.deltaFdTrueHz(:))), [], 'omitnan');

for iCase = 1:numel(caseList)
  rowIdx = iCase + 1;
  caseInfo = caseList(iCase);
  displayName(rowIdx) = string(caseInfo.displayName);
  if ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
    continue;
  end
  estResult = caseInfo.estResult;
  aux = getDoaDopplerFieldOrDefault(estResult, 'aux', struct());
  fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
  fdSatEst = reshape(getDoaDopplerFieldOrDefault(aux, 'fdSatEst', []), [], 1);
  deltaFdEst = reshape(getDoaDopplerFieldOrDefault(aux, 'deltaFdRefEst', []), [], 1);
  if numel(fdSatEst) < refIdx || numel(deltaFdEst) < refIdx || ~isfinite(fdRefEst)
    continue;
  end
  deltaFdAtRefHz(rowIdx) = deltaFdEst(refIdx);
  fdSatRefMinusFdRefHz(rowIdx) = fdSatEst(refIdx) - fdRefEst;
  maxFdClosureAbsHz(rowIdx) = max(abs(fdSatEst(:) - (fdRefEst + deltaFdEst(:))), [], 'omitnan');
end

invariantTable = table(displayName, deltaFdAtRefHz, fdSatRefMinusFdRefHz, maxFdClosureAbsHz);
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

function plotData = localBuildPlotData(repeatTable, weightTable, caseSummaryTable)
%LOCALBUILDPLOTDATA Store lightweight arrays for replottable figures.

plotData = struct();
plotData.repeatTable = repeatTable;
plotData.weightTable = weightTable;
plotData.caseSummaryTable = caseSummaryTable;
end

function plotData = localPlotReplay(replayData)
%LOCALPLOTREPLAY Draw compact replay diagnostics from replayData only.

plotData = replayData.plotData;
localPlotStaticRepeatScatter(replayData.repeatCompareTable);
localPlotWeightSweep(replayData.weightPreferenceTable);
localPlotRepresentativeGain(replayData.repeatCompareTable);
plotData.figureGenerated = true;
end

function localPlotStaticRepeatScatter(repeatTable)
%LOCALPLOTSTATICREPEATSCATTER Plot one compact SS/MS repeat scatter view.

if isempty(repeatTable) || height(repeatTable) == 0
  return;
end
maskAngle = repeatTable.ssResolved & repeatTable.msResolved & ...
  isfinite(repeatTable.ssAngleErrDeg) & isfinite(repeatTable.msAngleErrDeg);
maskFd = maskAngle & isfinite(repeatTable.ssFdErrHz) & isfinite(repeatTable.msFdErrHz);

figure();
subplot(1, 2, 1);
scatter(repeatTable.ssAngleErrDeg(maskAngle), repeatTable.msAngleErrDeg(maskAngle), 28, 'filled');
hold on;
xyMax = max([repeatTable.ssAngleErrDeg(maskAngle); repeatTable.msAngleErrDeg(maskAngle)], [], 'omitnan');
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
xyMax = max([repeatTable.ssFdErrHz(maskFd); repeatTable.msFdErrHz(maskFd)], [], 'omitnan');
if isempty(xyMax) || ~isfinite(xyMax)
  xyMax = 1;
end
plot([0, xyMax], [0, xyMax], '--k', 'LineWidth', 1.1);
grid on;
xlabel('SS-SF-Static fd error (Hz)');
ylabel('MS-SF-Static fd error (Hz)');
title('Repeat-level static fd compare');
end

function localPlotWeightSweep(weightPreferenceTable)
%LOCALPLOTWEIGHTSWEEP Plot one compact fixed-SNR weight sweep summary.

if isempty(weightPreferenceTable) || height(weightPreferenceTable) == 0
  return;
end
figure();
subplot(1, 2, 1);
plot(weightPreferenceTable.alphaSat2, weightPreferenceTable.angleRmseDeg, 'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Angle RMSE (deg)');
title('Fixed-SNR static angle vs sat2 weight');

subplot(1, 2, 2);
plot(weightPreferenceTable.alphaSat2, weightPreferenceTable.fdRmseHz, 'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
grid on;
xlabel('sat2 weight \alpha');
ylabel('Reference Doppler RMSE (Hz)');
title('Fixed-SNR static fd vs sat2 weight');
end

function localPlotRepresentativeGain(repeatTable)
%LOCALPLOTREPRESENTATIVEGAIN Plot repeat-wise MS-SS gain traces.

if isempty(repeatTable) || height(repeatTable) == 0
  return;
end
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

function metricLineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build replay-specific HTML-ready metric lines.

metricLineList = strings(0, 1);
caseSummary = replayData.caseSummaryTable;
ssRow = localSelectCaseRow(caseSummary, "SS-SF-Static");
msRow = localSelectCaseRow(caseSummary, "MS-SF-Static");
gainTable = replayData.staticGainSummaryTable;
metricLineList(end + 1, 1) = "• static anchor: SNR=<code>" + string(replayData.config.snrDb) + ...
  " dB</code>, repeats=<code>" + string(replayData.config.numRepeat) + "</code>";
metricLineList(end + 1, 1) = "• SS/MS angle RMSE: <code>" + ...
  localFormatNumber(localTableValue(ssRow, 'angleRmseDeg', NaN)) + " / " + ...
  localFormatNumber(localTableValue(msRow, 'angleRmseDeg', NaN)) + " deg</code>";
metricLineList(end + 1, 1) = "• SS/MS fd RMSE: <code>" + ...
  localFormatNumber(localTableValue(ssRow, 'fdRmseHz', NaN)) + " / " + ...
  localFormatNumber(localTableValue(msRow, 'fdRmseHz', NaN)) + " Hz</code>";
if ~isempty(gainTable)
  msBetterRate = localMetricValue(gainTable, "msBetterAngleRate", NaN);
  metricLineList(end + 1, 1) = "• MS angle-better rate: <code>" + localFormatNumber(msBetterRate) + "</code>";
end
end

function value = localMetricValue(metricTable, metricName, defaultValue)
%LOCALMETRICVALUE Read one metric value from a metricName/metricValue table.

value = defaultValue;
if isempty(metricTable) || ~istable(metricTable) || ~all(ismember({'metricName', 'metricValue'}, metricTable.Properties.VariableNames))
  return;
end
idx = find(string(metricTable.metricName) == string(metricName), 1, 'first');
if ~isempty(idx)
  value = metricTable.metricValue(idx);
end
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
    error('replaySfStaticDualSatDiagnose:InvalidSelectMode', ...
      'Unsupported select mode: %s.', modeName);
end
finiteIdx = find(finiteMask);
selIdx = finiteIdx(idxLocal);
end

function [selIdx, hasValue] = localSelectMedianResolved(metricVec, resolvedMask)
%LOCALSELECTMEDIANRESOLVED Select one repeat closest to resolved median.

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

function runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only fields accepted by checkpoint runner.

runnerOpt = struct();
runnerOpt.runName = checkpointOpt.runName;
runnerOpt.runKey = checkpointOpt.runKey;
runnerOpt.outputRoot = checkpointOpt.outputRoot;
runnerOpt.useParfor = logical(getDoaDopplerFieldOrDefault(checkpointOpt, 'useParfor', false));
runnerOpt.resume = logical(getDoaDopplerFieldOrDefault(checkpointOpt, 'resume', true));
runnerOpt.meta = getDoaDopplerFieldOrDefault(checkpointOpt, 'meta', struct());
runnerOpt.cleanupOnSuccess = false;
runnerOpt.cleanupOpt = struct();
end

function meta = localBuildCheckpointMeta(config)
%LOCALBUILDCHECKPOINTMETA Build lightweight checkpoint manifest metadata.

meta = struct();
meta.snrDb = config.snrDb;
meta.seedList = reshape(config.seedList, 1, []);
meta.weightSweepAlpha = reshape(config.weightSweepAlpha, 1, []);
meta.staticMsHalfWidthDeg = reshape(config.staticMsHalfWidthDeg, 1, []);
end

function runKey = localBuildRunKey(config)
%LOCALBUILDRUNKEY Build a short stable checkpoint run key.

seedList = reshape(config.seedList, [], 1);
runKey = sprintf('snr%s_seed%dto%d_rep%d_static', ...
  localCompactNumber(config.snrDb), seedList(1), seedList(end), numel(seedList));
runKey = string(runKey);
end

function textValue = localCompactNumber(value)
%LOCALCOMPACTNUMBER Format one numeric value for run keys.

textValue = sprintf('%.6g', double(value));
textValue = strrep(textValue, '-', 'm');
textValue = strrep(textValue, '.', 'p');
end

function localPrintReplayConfig(config)
%LOCALPRINTREPLAYCONFIG Print replay-specific axes and static settings.

fprintf('  %-32s : %s\n', 'seed list', localFormatRow(config.seedList));
fprintf('  %-32s : %.2f dB\n', 'static SNR', config.snrDb);
fprintf('  %-32s : %s\n', 'weight alpha list', localFormatRow(config.weightSweepAlpha));
fprintf('  %-32s : [%s] deg\n', 'MS static DoA half-width', localFormatRow(config.staticMsHalfWidthDeg));
fprintf('  %-32s : %d / %d\n', 'task count / outer parfor', config.numTask, config.useParfor);
fprintf('  %-32s : %s\n', 'fixture helper', 'buildSfStaticSingleSnrRepeatFixture');
end

function textValue = localFormatRow(valueList)
%LOCALFORMATROW Format numeric values as one compact row.

valueList = reshape(double(valueList), 1, []);
if isempty(valueList)
  textValue = '';
  return;
end
if numel(valueList) > 2 && all(abs(diff(valueList) - 1) < 1e-12)
  textValue = sprintf('%.6g:%.6g (%d values)', valueList(1), valueList(end), numel(valueList));
else
  textValue = strjoin(compose('%.6g', valueList), ', ');
end
end

function value = localMeanLogical(mask)
%LOCALMEANLOGICAL Return mean(mask) with empty handling.

if isempty(mask)
  value = NaN;
  return;
end
value = mean(double(mask));
end

function textValue = localFormatNumber(value)
%LOCALFORMATNUMBER Format one number for compact text.

if ~(isnumeric(value) || islogical(value)) || isempty(value) || ~isfinite(double(value))
  textValue = "NaN";
  return;
end
textValue = string(sprintf('%.4g', double(value)));
end

function tf = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Check whether outer parfor can be used safely.

if nargin < 1 || isempty(numTask)
  numTask = 1;
end
tf = false;
if numTask < 2
  return;
end
try
  if exist('getCurrentTask', 'file') && ~isempty(getCurrentTask())
    return;
  end
catch
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  tf = ~isempty(ver('parallel'));
catch
  tf = false;
end
end
