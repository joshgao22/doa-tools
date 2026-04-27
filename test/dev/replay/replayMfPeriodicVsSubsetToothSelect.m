% replayMfPeriodicVsSubsetToothSelect
% Purpose: compare three stages in the dynamic flow: periodic-wide solving,
% selected subset tooth selection, and final periodic in-tooth refinement.
% The replay highlights whether non-periodic subsets improve tooth selection
% and whether the final periodic stage mainly refines within the selected tooth.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted repeat runs from per-repeat files under the repo-root tmp/.

clear; close all; clc;

%% Replay configuration
% The representative scenario, satellite pair, waveform, frame offsets, and
% default search boxes are defined by buildDynamicDualSatEciContext through
% runSimpleDynamicFlowReplayBatch. This replay only controls the small-MC
% repeat batch and the subset-periodic flow comparison.

replayName = "replayMfPeriodicVsSubsetToothSelect";

snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % First seed of the small replay batch.
numRepeat = 100;                    % Number of consecutive seeds to replay.
saveSnapshot = true;                % true saves the lightweight replayData only.
optVerbose = false;                 % true enables compact estimator/flow trace.
checkpointEnable = true;            % true enables per-repeat checkpoint/resume under repo-root tmp/.
toothHistogramBinCount = 8;         % Target number of |toothIdx| histogram bins; fewer bins merge wider tooth ranges.

% Checkpoint resume and success cleanup are fixed for this replay: an interrupted
% run resumes from existing task files, and a successful run removes them after
% replayData is assembled.

seedList = baseSeed + (0:(numRepeat - 1));  % Consecutive task seeds for the batch.
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

%% Build context and flow options

config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.saveSnapshot = saveSnapshot;
config.optVerbose = optVerbose;
config.checkpointEnable = checkpointEnable;
config.toothHistogramBinCount = toothHistogramBinCount;
config.checkpointResume = true;
config.checkpointCleanupOnSuccess = true;

fprintf('Running %s ...\n', char(replayName));
fprintf('  repeats                         : %d\n', config.numRepeat);
fprintf('  snr (dB)                        : %.2f\n', config.snrDb);
fprintf('  base seed                       : %d\n', config.baseSeed);
fprintf('  repeat mode                     : %s\n', 'parfor-auto');
fprintf('  save snapshot                   : %d\n', config.saveSnapshot);
fprintf('  checkpoint                      : %d\n', config.checkpointEnable);

% Build compact flow options here so the replay documents exactly which stages
% are being compared. The shared helper still owns fixture and signal creation.
replayBatchOpt = localBuildBatchOpt(replayName, config);
checkpointRunDir = "";
if config.checkpointEnable
  checkpointRunDir = replayBatchOpt.checkpointOpt.runDir;
  config.checkpointRunDir = string(checkpointRunDir);
  fprintf('  checkpoint resume               : %d\n', config.checkpointResume);
  fprintf('  checkpoint dir                  : %s\n', checkpointRunDir);
end

%% Run replay batch

% The batch helper owns the progressbar for the outer repeat loop. It builds
% the shared context, generates one repeat per seed, runs periodic-wide solving,
% then runs the current subset-periodic flow for each seed.
if config.checkpointEnable
  try
    [repeatTable, repeatCell, context, flowOpt, runState] = runSimpleDynamicFlowReplayBatch(replayBatchOpt);
  catch ME
    fprintf('Replay failed. Checkpoint artifacts kept at: %s\n', checkpointRunDir);
    rethrow(ME);
  end
else
  [repeatTable, repeatCell, context, flowOpt, runState] = runSimpleDynamicFlowReplayBatch(replayBatchOpt);
end
compareTable = localBuildCompareTable(repeatTable, repeatCell);
aggregateTable = localBuildAggregateTable(compareTable);
warningTable = localBuildWarningTable(repeatTable);
toothHistogramTable = localBuildToothHistogramTable(compareTable, config.toothHistogramBinCount);

%% Data storage

replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, flowOpt);
replayData.compareTable = compareTable;
replayData.aggregateTable = aggregateTable;
replayData.warningTable = warningTable;
replayData.toothHistogramTable = toothHistogramTable;
replayData.repeatTable = repeatTable;
replayData.repeatCell = localStripRepeatCell(repeatCell);
if config.checkpointEnable
  replayData.checkpointSummary = localBuildCheckpointSummary(runState);
end

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

if ~exist('replayData', 'var') || ~isstruct(replayData)
  error('Replay data is missing. Run the replay batch sections or load a snapshot containing replayData.');
end
config = replayData.config;
compareTable = replayData.compareTable;
aggregateTable = replayData.aggregateTable;
if isfield(replayData, 'warningTable')
  warningTable = replayData.warningTable;
else
  warningTable = table();
end
toothHistogramTable = replayData.toothHistogramTable;

fprintf('Running replayMfPeriodicVsSubsetToothSelect ...\n');
fprintf('  repeats                         : %d\n', height(compareTable));
fprintf('  snr (dB)                        : %.2f\n', config.snrDb);
fprintf('  base seed                       : %d\n', config.baseSeed);
if isfield(config, 'checkpointEnable')
  fprintf('  checkpoint                      : %d\n', config.checkpointEnable);
end
if isfield(config, 'checkpointRunDir')
  fprintf('  checkpoint dir                  : %s\n', char(config.checkpointRunDir));
end
fprintf('\n========== Periodic vs subset compare ==========\n');
localDisplayCompareTable(compareTable);
fprintf('\n========== Aggregate summary ==========\n');
disp(aggregateTable);
fprintf('\n========== Warning summary ==========\n');
localDisplayWarningTable(warningTable);
% The tooth histogram is saved in replayData and visualized below, but it is
% not printed because the long per-bin table is noisy for routine replay runs.
localPlotReplay(compareTable, toothHistogramTable);

%% Local helpers

function replayBatchOpt = localBuildBatchOpt(replayName, config)
% Build replay batch options while keeping estimator internals deterministic.
parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);

% The final periodic stage is intentionally narrow in DoA so this replay focuses
% on tooth selection first and in-tooth fd/fdRate refinement second.
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti", ...
  'periodicPolishEnableWhenMulti', true, ...
  'periodicPolishDoaHalfWidthDeg', [0.002; 0.002]));

contextOpt = struct( ...
  'baseSeed', config.baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt);

checkpointOpt = struct('enable', config.checkpointEnable);
if config.checkpointEnable
  checkpointRunKey = localBuildCheckpointRunKey(config);
  checkpointOutputRoot = fullfile(localGetRepoRoot(), 'tmp');
  checkpointOpt.resume = true;
  checkpointOpt.runName = string(replayName);
  checkpointOpt.runKey = checkpointRunKey;
  checkpointOpt.outputRoot = checkpointOutputRoot;
  checkpointOpt.runDir = fullfile(checkpointOutputRoot, char(replayName), char(checkpointRunKey));
end

replayBatchOpt = struct( ...
  'replayName', string(replayName), ...
  'snrDb', config.snrDb, ...
  'baseSeed', config.baseSeed, ...
  'seedList', config.seedList, ...
  'numRepeat', config.numRepeat, ...
  'optVerbose', config.optVerbose, ...
  'parallelOpt', parallelOpt, ...
  'contextOpt', contextOpt, ...
  'flowOpt', flowOpt, ...
  'checkpointOpt', checkpointOpt, ...
  'runPeriodicWide', true, ...
  'progressTitle', sprintf('%s repeat batch', char(replayName)));
end


function runKey = localBuildCheckpointRunKey(config)
% Build a stable checkpoint key for this replay configuration.
seedList = reshape(double(config.seedList), 1, []);
runKey = string(sprintf('seed%dto%d_rep%d_snr%.2f', ...
  seedList(1), seedList(end), numel(seedList), config.snrDb));
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
end


function repoRoot = localGetRepoRoot()
% Resolve the repository root from the replay script location.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
end


function checkpointSummary = localBuildCheckpointSummary(runState)
% Store lightweight checkpoint status without keeping checkpoint task results twice.
checkpointSummary = struct();
checkpointSummary.runDir = string(localGetFieldOrDefault(runState, 'runDir', ""));
checkpointSummary.numTask = localGetFieldOrDefault(runState, 'numTask', 0);
checkpointSummary.numDone = localGetFieldOrDefault(runState, 'numDone', 0);
checkpointSummary.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
checkpointSummary.cleanedOnSuccess = false;
end


function compareTable = localBuildCompareTable(repeatTable, repeatCell)
% Assemble per-repeat stage metrics into one compact comparison table.
numRepeat = height(repeatTable);
taskSeed = repeatTable.taskSeed;
selectedSubsetLabel = repeatTable.selectedSubsetLabel;
periodicWideToothIdx = repeatTable.periodicWideToothIdx;
subsetToothIdx = nan(numRepeat, 1);
finalToothIdx = repeatTable.toothIdx;
periodicWideAngleErrDeg = repeatTable.periodicWideAngleErrDeg;
subsetAngleErrDeg = nan(numRepeat, 1);
finalAngleErrDeg = repeatTable.angleErrDeg;
periodicWideFdRefErrHz = repeatTable.periodicWideFdRefErrHz;
subsetFdRefErrHz = nan(numRepeat, 1);
finalFdRefErrHz = repeatTable.fdRefErrHz;
usedPolish = repeatTable.usedPolish;
warningSeen = localGetTableColumnOrDefault(repeatTable, 'warningSeen', false(numRepeat, 1));
warningId = localGetTableColumnOrDefault(repeatTable, 'warningId', strings(numRepeat, 1));
for iRepeat = 1:numRepeat
  selectedSummary = repeatCell{iRepeat}.msFlow.selectedSubsetSummary;
  subsetToothIdx(iRepeat) = localGetFieldOrDefault(selectedSummary, 'toothIdx', NaN);
  subsetAngleErrDeg(iRepeat) = localGetFieldOrDefault(selectedSummary, 'angleErrDeg', NaN);
  subsetFdRefErrHz(iRepeat) = localGetFieldOrDefault(selectedSummary, 'fdRefErrHz', NaN);
end
compareTable = table(taskSeed, selectedSubsetLabel, periodicWideToothIdx, subsetToothIdx, finalToothIdx, ...
  periodicWideAngleErrDeg, subsetAngleErrDeg, finalAngleErrDeg, ...
  periodicWideFdRefErrHz, subsetFdRefErrHz, finalFdRefErrHz, usedPolish, warningSeen, warningId);
end


function aggregateTable = localBuildAggregateTable(compareTable)
% Summarize whether tooth estimates become more concentrated across stages.
stageName = ["periodicWide"; "selectedSubset"; "finalPeriodicRefine"];
[periodicHit, periodicWithinOne, periodicMeanAbs, periodicMedianAbs] = ...
  localBuildToothStats(compareTable.periodicWideToothIdx);
[subsetHit, subsetWithinOne, subsetMeanAbs, subsetMedianAbs] = ...
  localBuildToothStats(compareTable.subsetToothIdx);
[finalHit, finalWithinOne, finalMeanAbs, finalMedianAbs] = ...
  localBuildToothStats(compareTable.finalToothIdx);
strictToothHitRate = [periodicHit; subsetHit; finalHit];
withinOneToothRate = [periodicWithinOne; subsetWithinOne; finalWithinOne];
meanAbsToothIdx = [periodicMeanAbs; subsetMeanAbs; finalMeanAbs];
medianAbsToothIdx = [periodicMedianAbs; subsetMedianAbs; finalMedianAbs];
angleRmseDeg = [localRmse(compareTable.periodicWideAngleErrDeg); ...
  localRmse(compareTable.subsetAngleErrDeg); ...
  localRmse(compareTable.finalAngleErrDeg)];
fdRefRmseHz = [localRmse(compareTable.periodicWideFdRefErrHz); ...
  localRmse(compareTable.subsetFdRefErrHz); ...
  localRmse(compareTable.finalFdRefErrHz)];
aggregateTable = table(stageName, strictToothHitRate, withinOneToothRate, ...
  meanAbsToothIdx, medianAbsToothIdx, angleRmseDeg, fdRefRmseHz);
end



function warningTable = localBuildWarningTable(repeatTable)
% Summarize non-fatal MATLAB warnings observed during replay repeats.
numRepeat = height(repeatTable);
if numRepeat == 0 || ~ismember('warningSeen', repeatTable.Properties.VariableNames)
  warningTable = table();
  return;
end
warningSeen = logical(repeatTable.warningSeen);
if ~any(warningSeen)
  warningTable = table();
  return;
end
warningIdList = string(localGetTableColumnOrDefault(repeatTable, 'warningId', strings(numRepeat, 1)));
warningMessageList = string(localGetTableColumnOrDefault(repeatTable, 'warningMessage', strings(numRepeat, 1)));
taskSeedList = repeatTable.taskSeed;
warningKeyList = warningIdList;
emptyIdMask = strlength(warningKeyList) == 0;
warningKeyList(emptyIdMask) = localTruncateText(warningMessageList(emptyIdMask), 80);
warningKeyList(strlength(warningKeyList) == 0) = "unknown-warning";
uniqueKeyList = unique(warningKeyList(warningSeen), 'stable');
numWarning = numel(uniqueKeyList);
warningId = strings(numWarning, 1);
repeatCount = zeros(numWarning, 1);
taskSeedSummary = strings(numWarning, 1);
exampleMessage = strings(numWarning, 1);
for iWarning = 1:numWarning
  rowMask = warningSeen & warningKeyList == uniqueKeyList(iWarning);
  rowIdx = find(rowMask);
  warningId(iWarning) = warningIdList(rowIdx(1));
  if strlength(warningId(iWarning)) == 0
    warningId(iWarning) = uniqueKeyList(iWarning);
  end
  repeatCount(iWarning) = numel(rowIdx);
  taskSeedSummary(iWarning) = localFormatSeedList(taskSeedList(rowIdx), 12);
  exampleMessage(iWarning) = localTruncateText(warningMessageList(rowIdx(1)), 180);
end
warningTable = table(warningId, repeatCount, taskSeedSummary, exampleMessage);
end


function localDisplayCompareTable(compareTable)
% Print a bounded preview while keeping the full table in replayData.
maxFullRow = 20;
previewEdgeRow = 8;
numRow = height(compareTable);
if numRow <= maxFullRow
  disp(compareTable);
  return;
end
headIdx = 1:min(previewEdgeRow, numRow);
tailIdx = max(headIdx(end) + 1, numRow - previewEdgeRow + 1):numRow;
previewTable = [compareTable(headIdx, :); compareTable(tailIdx, :)];
fprintf('Full compare table has %d rows; showing first %d and last %d rows. Full table is saved in replayData.compareTable.\n', ...
  numRow, numel(headIdx), numel(tailIdx));
disp(previewTable);
end


function localDisplayWarningTable(warningTable)
% Print a compact warning summary without expanding every repeat row.
if isempty(warningTable)
  fprintf('No non-fatal MATLAB warnings were recorded by the repeat runner.\n');
else
  disp(warningTable);
end
end


function toothHistogramTable = localBuildToothHistogramTable(compareTable, binCount)
% Count absolute-tooth concentration using a target number of bins.
stageNameList = ["periodicWide"; "selectedSubset"; "finalPeriodicRefine"];
toothCell = {compareTable.periodicWideToothIdx; compareTable.subsetToothIdx; compareTable.finalToothIdx};
binCount = max(1, round(double(binCount)));
allAbsTooth = [];
for iStage = 1:numel(toothCell)
  allAbsTooth = [allAbsTooth; localAbsFiniteTooth(toothCell{iStage})]; %#ok<AGROW>
end
maxAbsTooth = max([0; allAbsTooth]);
numToothValue = maxAbsTooth + 1;
numRequestedBin = min(binCount, numToothValue);
binEdgeList = round(linspace(0, numToothValue, numRequestedBin + 1)).';
binLeftList = binEdgeList(1:end - 1);
binRightList = binEdgeList(2:end) - 1;
numStage = numel(stageNameList);
numBin = numel(binLeftList);
rowStageName = strings(numStage * numBin, 1);
rowBinLeftAbsToothIdx = zeros(numStage * numBin, 1);
rowBinRightAbsToothIdx = zeros(numStage * numBin, 1);
rowBinCenterAbsToothIdx = zeros(numStage * numBin, 1);
rowBinLabel = strings(numStage * numBin, 1);
rowRepeatCount = zeros(numStage * numBin, 1);
rowRepeatRate = nan(numStage * numBin, 1);
rowIdx = 0;
for iStage = 1:numStage
  absTooth = localAbsFiniteTooth(toothCell{iStage});
  totalCount = numel(absTooth);
  for iBin = 1:numBin
    rowIdx = rowIdx + 1;
    binLeft = binLeftList(iBin);
    binRight = binRightList(iBin);
    countUse = sum(absTooth >= binLeft & absTooth <= binRight);
    rowStageName(rowIdx) = stageNameList(iStage);
    rowBinLeftAbsToothIdx(rowIdx) = binLeft;
    rowBinRightAbsToothIdx(rowIdx) = binRight;
    rowBinCenterAbsToothIdx(rowIdx) = 0.5 * (binLeft + binRight);
    rowBinLabel(rowIdx) = localFormatToothBinLabel(binLeft, binRight);
    rowRepeatCount(rowIdx) = countUse;
    if totalCount > 0
      rowRepeatRate(rowIdx) = countUse / totalCount;
    end
  end
end
toothHistogramTable = table(rowStageName, rowBinLeftAbsToothIdx, rowBinRightAbsToothIdx, ...
  rowBinCenterAbsToothIdx, rowBinLabel, rowRepeatCount, rowRepeatRate, ...
  'VariableNames', {'stageName', 'binLeftAbsToothIdx', 'binRightAbsToothIdx', ...
  'binCenterAbsToothIdx', 'binLabel', 'repeatCount', 'repeatRate'});
end

function localPlotReplay(compareTable, toothHistogramTable)
% Plot stage traces and the absolute-tooth concentration histogram.
figure('Name', 'Periodic vs subset flow stages');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

stageLegend = {'periodic wide', 'selected subset', 'final periodic refine'};

nexttile;
plot(compareTable.taskSeed, compareTable.periodicWideToothIdx, '-o', 'DisplayName', stageLegend{1}); hold on;
plot(compareTable.taskSeed, compareTable.subsetToothIdx, '-s', 'DisplayName', stageLegend{2});
plot(compareTable.taskSeed, compareTable.finalToothIdx, '-^', 'DisplayName', stageLegend{3});
grid on;
xlabel('task seed');
ylabel('tooth index');
title('stage tooth index');
legend('Location', 'best');

nexttile;
[countMat, binLabelList, stageNameList] = localBuildHistogramMatrix(toothHistogramTable);
if isempty(binLabelList)
  text(0.5, 0.5, 'no finite tooth data', 'HorizontalAlignment', 'center');
  axis off;
else
  binX = 1:numel(binLabelList);
  bar(binX, countMat, 'grouped');
  grid on;
  xticks(binX);
  xticklabels(cellstr(binLabelList));
  xlabel('|tooth index| bin');
  ylabel('repeat count');
  title('tooth concentration');
  legend(stageLegend(1:numel(stageNameList)), 'Location', 'best');
end

nexttile;
plot(compareTable.taskSeed, compareTable.periodicWideAngleErrDeg, '-o', 'DisplayName', stageLegend{1}); hold on;
plot(compareTable.taskSeed, compareTable.subsetAngleErrDeg, '-s', 'DisplayName', stageLegend{2});
plot(compareTable.taskSeed, compareTable.finalAngleErrDeg, '-^', 'DisplayName', stageLegend{3});
grid on;
xlabel('task seed');
ylabel('angle err (deg)');
title('angle error');
legend('Location', 'best');

nexttile;
plot(compareTable.taskSeed, compareTable.periodicWideFdRefErrHz, '-o', 'DisplayName', stageLegend{1}); hold on;
plot(compareTable.taskSeed, compareTable.subsetFdRefErrHz, '-s', 'DisplayName', stageLegend{2});
plot(compareTable.taskSeed, compareTable.finalFdRefErrHz, '-^', 'DisplayName', stageLegend{3});
grid on;
xlabel('task seed');
ylabel('fdRef err (Hz)');
title('fdRef error');
legend('Location', 'best');
end


function contextSummary = localBuildContextSummary(context, flowOpt)
% Store the resolved shared context settings needed to interpret the replay.
contextSummary = struct();
contextSummary.frameIntvlSec = context.frameIntvlSec;
contextSummary.periodicOffsetIdx = context.periodicOffsetIdx;
contextSummary.masterOffsetIdx = context.masterOffsetIdx;
contextSummary.subsetDefaultBankLabelList = flowOpt.subsetDefaultBankLabelList;
contextSummary.subsetMarginFallbackBankLabelList = flowOpt.subsetMarginFallbackBankLabelList;
if isfinite(context.frameIntvlSec) && context.frameIntvlSec > 0
  contextSummary.toothStepHz = 1 / context.frameIntvlSec;
else
  contextSummary.toothStepHz = NaN;
end
end


function repeatCellSlim = localStripRepeatCell(repeatCell)
% Remove heavy fields while keeping the candidate tables needed for inspection.
repeatCellSlim = cell(size(repeatCell));
for iCell = 1:numel(repeatCell)
  item = repeatCell{iCell};
  repeatCellSlim{iCell} = struct('taskSeed', item.taskSeed, 'summary', item.summary, ...
    'periodicWideSummary', item.periodicWideSummary, ...
    'selectedSubsetSummary', item.msFlow.selectedSubsetSummary, ...
    'subsetCandidateTable', item.msFlow.subsetCandidateTable, ...
    'periodicCandidateTable', item.msFlow.periodicCandidateTable, ...
    'periodicDoaSeed', item.msFlow.periodicDoaSeed, ...
    'warning', localGetFieldOrDefault(item, 'warning', struct()));
end
end


function [hitRate, withinOneRate, meanAbsTooth, medianAbsTooth] = localBuildToothStats(toothIdx)
% Compute compact concentration metrics from finite tooth-index samples.
validTooth = toothIdx(isfinite(toothIdx));
if isempty(validTooth)
  hitRate = NaN;
  withinOneRate = NaN;
  meanAbsTooth = NaN;
  medianAbsTooth = NaN;
  return;
end
absTooth = abs(validTooth(:));
hitRate = mean(absTooth == 0);
withinOneRate = mean(absTooth <= 1);
meanAbsTooth = mean(absTooth);
medianAbsTooth = median(absTooth);
end


function absTooth = localAbsFiniteTooth(toothIdx)
% Convert finite tooth indices to nonnegative integer offset bins.
validTooth = toothIdx(isfinite(toothIdx));
absTooth = abs(round(validTooth(:)));
end


function value = localRmse(metricVec)
% Compute RMSE over finite entries only.
metricVec = metricVec(isfinite(metricVec));
if isempty(metricVec)
  value = NaN;
else
  value = sqrt(mean(metricVec(:).^2));
end
end


function [countMat, binLabelList, stageNameList] = localBuildHistogramMatrix(toothHistogramTable)
% Convert the long histogram table to a grouped-bar matrix.
if isempty(toothHistogramTable)
  countMat = [];
  binLabelList = strings(0, 1);
  stageNameList = strings(0, 1);
  return;
end
stageNameList = unique(toothHistogramTable.stageName, 'stable');
binLabelList = unique(toothHistogramTable.binLabel, 'stable');
countMat = zeros(numel(binLabelList), numel(stageNameList));
for iStage = 1:numel(stageNameList)
  for iBin = 1:numel(binLabelList)
    isRow = toothHistogramTable.stageName == stageNameList(iStage) & ...
      toothHistogramTable.binLabel == binLabelList(iBin);
    if any(isRow)
      countMat(iBin, iStage) = toothHistogramTable.repeatCount(find(isRow, 1, 'first'));
    end
  end
end
end


function binLabel = localFormatToothBinLabel(binLeft, binRight)
% Format one absolute-tooth histogram bin label.
if binLeft == binRight
  binLabel = string(sprintf('%d', binLeft));
else
  binLabel = string(sprintf('%d-%d', binLeft, binRight));
end
end



function value = localGetTableColumnOrDefault(dataTable, columnName, defaultValue)
% Read a table column when present and nonempty, otherwise return a default.
value = defaultValue;
if istable(dataTable) && ismember(columnName, dataTable.Properties.VariableNames)
  rawValue = dataTable.(columnName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end


function textOut = localFormatSeedList(seedList, maxCount)
% Format a bounded seed list for warning summaries.
seedList = reshape(double(seedList), 1, []);
if isempty(seedList)
  textOut = "";
  return;
end
shownSeed = seedList(1:min(maxCount, numel(seedList)));
seedText = join(string(shownSeed), ", ");
if numel(seedList) > numel(shownSeed)
  seedText = seedText + sprintf(", ... +%d more", numel(seedList) - numel(shownSeed));
end
textOut = string(seedText);
end


function textOut = localTruncateText(textIn, maxChar)
% Truncate long warning text for compact command-window output.
textOut = string(textIn);
for iText = 1:numel(textOut)
  if strlength(textOut(iText)) > maxChar
    textOut(iText) = extractBefore(textOut(iText), maxChar + 1) + "...";
  end
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
% Read a struct field when present and nonempty, otherwise return a default.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
