% replayMfFlowLikeGatedBasinEntryEffectiveness
% Purpose: validate the no-truth same-tooth basin-entry rescue inside the
% current subset-periodic flow. The replay compares the default flow against
% a gated wide+single-MF rescue bank without changing estimator internals.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; checkpointEnable=true resumes
% interrupted repeat runs from per-repeat files under the repo-root tmp/.
% saveSnapshot=true saves only replayData via saveExpSnapshot. Telegram notice is best-effort only.

clear; close all; clc;

%% Replay configuration
replayName = "replayMfFlowLikeGatedBasinEntryEffectiveness";

snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % First seed of the small replay batch.
numRepeat = 100;                    % Number of consecutive seeds to replay when focusedSeedList is empty.
focusedSeedList = [];               % Focused seeds for in-tooth/gate-miss triage; set [] for consecutive seeds.
saveSnapshot = true;                % true saves replayData via saveExpSnapshot.
notifyTelegramEnable = true;        % true sends best-effort HTML Telegram notice on completion or failure.
saveRepeatDetail = false;           % false keeps only diagnostic tables in replayData to avoid large snapshots.
optVerbose = false;                 % true enables compact estimator / flow trace.
checkpointEnable = true;            % true enables per-repeat checkpoint/resume under repo-root tmp/.
hardAngleThresholdDeg = 2e-3;       % Offline hard-case label based on disabled-flow angle error.
easyAngleThresholdDeg = 1e-3;       % Offline easy-case label based on disabled-flow angle error.
angleGainThresholdDeg = 5e-4;       % Offline rescue gain threshold.
damageThresholdDeg = 5e-4;          % Offline damage threshold.
gatedRescueCoherenceThreshold = 0.20; % No-truth trigger threshold for non-ref coherence collapse.
gatedRescuePhaseResidThresholdRad = 1.0; % No-truth trigger threshold for non-ref phase residual.

runTic = tic;

% Checkpoint resume and success cleanup are fixed for this replay: interrupted
% runs resume from existing task files, and successful runs remove them after
% replayData is assembled.

if isempty(focusedSeedList)
  seedList = baseSeed + (0:(numRepeat - 1));
else
  seedList = focusedSeedList;
end
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

%% Build context and flow options
config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.focusedSeedList = reshape(double(focusedSeedList), 1, []);
config.saveSnapshot = saveSnapshot;
config.notifyTelegramEnable = notifyTelegramEnable;
config.saveRepeatDetail = saveRepeatDetail;
config.optVerbose = optVerbose;
config.checkpointEnable = checkpointEnable;
config.checkpointResume = true;
config.checkpointCleanupOnSuccess = true;
config.hardAngleThresholdDeg = hardAngleThresholdDeg;
config.easyAngleThresholdDeg = easyAngleThresholdDeg;
config.angleGainThresholdDeg = angleGainThresholdDeg;
config.damageThresholdDeg = damageThresholdDeg;
config.gatedRescueCoherenceThreshold = gatedRescueCoherenceThreshold;
config.gatedRescuePhaseResidThresholdRad = gatedRescuePhaseResidThresholdRad;

printMfReplayHeader(char(replayName), config, '');
fprintf('  seed list                       : %s\n', mat2str(reshape(config.seedList, 1, [])));
fprintf('  save repeat detail              : %d\n', config.saveRepeatDetail);
fprintf('  gated rescue coherence threshold: %.3f\n', config.gatedRescueCoherenceThreshold);
fprintf('  gated rescue phase threshold    : %.3f rad\n', config.gatedRescuePhaseResidThresholdRad);

%% Run replay batch
disabledBatchOpt = localBuildBatchOpt(replayName, config, "disabled");
gatedBatchOpt = localBuildBatchOpt(replayName, config, "gatedWideSingle");

if config.checkpointEnable
  config.checkpointDisabledRunDir = string(disabledBatchOpt.checkpointOpt.runDir);
  config.checkpointGatedRunDir = string(gatedBatchOpt.checkpointOpt.runDir);
  fprintf('  checkpoint resume               : %d\n', config.checkpointResume);
  fprintf('  disabled checkpoint dir         : %s\n', char(config.checkpointDisabledRunDir));
  fprintf('  gated checkpoint dir            : %s\n', char(config.checkpointGatedRunDir));
end

try
  fprintf('[%s] Run disabled same-tooth rescue flow.\n', localNowString('HH:mm:ss'));
  [disabledTable, disabledCell, context, disabledFlowOpt, disabledRunState] = runSimpleDynamicFlowReplayBatch(disabledBatchOpt);

  fprintf('[%s] Run gated wide+single-MF same-tooth rescue flow.\n', localNowString('HH:mm:ss'));
  [gatedTable, gatedCell, ~, gatedFlowOpt, gatedRunState] = runSimpleDynamicFlowReplayBatch(gatedBatchOpt);

compareTable = localBuildCompareTable(disabledTable, gatedTable, config);
aggregateTable = localBuildAggregateTable(compareTable);
candidateTable = localBuildCandidateTable(gatedCell);
warningTable = localBuildWarningTable(disabledTable, gatedTable);

%% Data storage
replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, disabledFlowOpt, gatedFlowOpt);
replayData.compareTable = compareTable;
replayData.inToothCompareTable = localBuildInToothCompareTable(compareTable);
replayData.gateMissTable = localBuildGateMissTable(compareTable);
replayData.aggregateTable = aggregateTable;
replayData.inToothAggregateTable = localBuildAggregateTable(replayData.inToothCompareTable);
replayData.caseRoleSummaryTable = localBuildCaseRoleSummaryTable(compareTable);
replayData.gateReasonTable = localBuildGateReasonTable(compareTable);
replayData.inToothGateReasonTable = localBuildGateReasonTable(replayData.inToothCompareTable);
replayData.finalTagTransitionTable = localBuildTransitionCountTable(compareTable, ...
  'disabledFinalTag', 'gatedFinalTag', 'disabledFinalTag', 'gatedFinalTag');
replayData.toothTransitionTable = localBuildTransitionCountTable(compareTable, ...
  'disabledToothIdx', 'gatedToothIdx', 'disabledToothIdx', 'gatedToothIdx');
replayData.candidateTable = candidateTable;
replayData.warningTable = warningTable;
replayData.disabledRepeatTable = disabledTable;
replayData.gatedRepeatTable = gatedTable;
if config.saveRepeatDetail
  replayData.disabledRepeatCell = localStripRepeatCell(disabledCell);
  replayData.gatedRepeatCell = localStripRepeatCell(gatedCell);
end
if config.checkpointEnable
  replayData.checkpointSummary = localBuildCheckpointSummary(disabledRunState, gatedRunState);
end

if config.checkpointEnable && config.checkpointCleanupOnSuccess
  replayData.checkpointSummary.disabled.cleanupReport = cleanupPerfTaskGridCheckpoint(disabledRunState, struct('verbose', false));
  replayData.checkpointSummary.disabled.cleanedOnSuccess = true;
  replayData.checkpointSummary.gated.cleanupReport = cleanupPerfTaskGridCheckpoint(gatedRunState, struct('verbose', false));
  replayData.checkpointSummary.gated.cleanedOnSuccess = true;
end
replayData.storageSummaryTable = localBuildStorageSummaryTable(replayData);

if config.saveSnapshot
  saveOpt = struct('includeVars', {{'replayData'}}, ...
    'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
  replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
else
  replayData.snapshotFile = "";
end

replayData.elapsedSec = toc(runTic);
notifyMfReplayStatus(struct( ...
  'replayName', replayName, ...
  'statusText', "DONE", ...
  'config', config, ...
  'snapshotFile', replayData.snapshotFile, ...
  'checkpointDir', localBuildCheckpointDirText(config), ...
  'elapsedSec', replayData.elapsedSec, ...
  'metricLineList', localBuildTelegramMetricLines(replayData), ...
  'commentLineList', "Flow-like gated basin-entry stress replay completed."));
catch ME
  if config.checkpointEnable
    fprintf('Replay failed. Checkpoint artifacts were kept for resume.\n');
    fprintf('  disabled checkpoint dir: %s\n', char(config.checkpointDisabledRunDir));
    fprintf('  gated checkpoint dir   : %s\n', char(config.checkpointGatedRunDir));
  end
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'checkpointDir', localBuildCheckpointDirText(config), ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Summary output and plotting
if ~exist('replayData', 'var') || ~isstruct(replayData)
  error('Replay data is missing. Run the replay batch sections or load a snapshot containing replayData.');
end
config = replayData.config;
compareTable = replayData.compareTable;
aggregateTable = replayData.aggregateTable;
candidateTable = replayData.candidateTable;
warningTable = replayData.warningTable;
if isfield(replayData, 'inToothCompareTable')
  inToothCompareTable = replayData.inToothCompareTable;
else
  inToothCompareTable = localBuildInToothCompareTable(compareTable);
end
if isfield(replayData, 'gateMissTable')
  gateMissTable = replayData.gateMissTable;
else
  gateMissTable = localBuildGateMissTable(compareTable);
end
if isfield(replayData, 'caseRoleSummaryTable')
  caseRoleSummaryTable = replayData.caseRoleSummaryTable;
else
  caseRoleSummaryTable = localBuildCaseRoleSummaryTable(compareTable);
end
if isfield(replayData, 'gateReasonTable')
  gateReasonTable = replayData.gateReasonTable;
else
  gateReasonTable = localBuildGateReasonTable(compareTable);
end
if isfield(replayData, 'inToothGateReasonTable')
  inToothGateReasonTable = replayData.inToothGateReasonTable;
else
  inToothGateReasonTable = localBuildGateReasonTable(inToothCompareTable);
end
if isfield(replayData, 'finalTagTransitionTable')
  finalTagTransitionTable = replayData.finalTagTransitionTable;
else
  finalTagTransitionTable = localBuildTransitionCountTable(compareTable, ...
    'disabledFinalTag', 'gatedFinalTag', 'disabledFinalTag', 'gatedFinalTag');
end
if isfield(replayData, 'toothTransitionTable')
  toothTransitionTable = replayData.toothTransitionTable;
else
  toothTransitionTable = localBuildTransitionCountTable(compareTable, ...
    'disabledToothIdx', 'gatedToothIdx', 'disabledToothIdx', 'gatedToothIdx');
end
if isfield(replayData, 'storageSummaryTable')
  storageSummaryTable = replayData.storageSummaryTable;
else
  storageSummaryTable = localBuildStorageSummaryTable(replayData);
end

fprintf('Running replayMfFlowLikeGatedBasinEntryEffectiveness ...\n');
fprintf('  repeats                         : %d\n', height(compareTable));
fprintf('  snr (dB)                        : %.2f\n', config.snrDb);
fprintf('  base seed                       : %d\n', config.baseSeed);
if isfield(config, 'seedList')
  fprintf('  seed list                       : %s\n', mat2str(reshape(config.seedList, 1, [])));
end
if isfield(config, 'notifyTelegramEnable')
  fprintf('  telegram notify                 : %d\n', config.notifyTelegramEnable);
end
if isfield(config, 'saveRepeatDetail')
  fprintf('  save repeat detail              : %d\n', config.saveRepeatDetail);
end
if isfield(config, 'checkpointEnable')
  fprintf('  checkpoint                      : %d\n', config.checkpointEnable);
end
if isfield(config, 'checkpointDisabledRunDir')
  fprintf('  disabled checkpoint dir         : %s\n', char(config.checkpointDisabledRunDir));
end
if isfield(config, 'checkpointGatedRunDir')
  fprintf('  gated checkpoint dir            : %s\n', char(config.checkpointGatedRunDir));
end
localPrintTablePreview('Disabled vs gated flow compare', compareTable, 16);
fprintf('\n========== Flow-like gated rescue aggregate ==========\n');
disp(aggregateTable);
localPrintTablePreview('Same-tooth in-tooth compare', inToothCompareTable, 16);
if isfield(replayData, 'inToothAggregateTable')
  fprintf('\n========== Same-tooth in-tooth aggregate ==========\n');
  disp(replayData.inToothAggregateTable);
end
localPrintTablePreview('Gate-miss hard cases', gateMissTable, 16);
localPrintTablePreview('Case-role summary', caseRoleSummaryTable, 16);
localPrintTablePreview('Gate reason counts', gateReasonTable, 16);
localPrintTablePreview('Same-tooth gate reason counts', inToothGateReasonTable, 16);
localPrintTablePreview('Final tag transition counts', finalTagTransitionTable, 16);
localPrintTablePreview('Tooth transition counts', toothTransitionTable, 16);
localPrintTablePreview('Gated periodic candidates', candidateTable, 16);
localPrintTablePreview('Warning summary', warningTable, 12);
localPrintTablePreview('ReplayData storage summary', storageSummaryTable, 32);
localPlotReplay(compareTable);

%% Local helpers

function replayBatchOpt = localBuildBatchOpt(replayName, config, modeName)
%LOCALBUILDBATCHOPT Build one disabled or gated flow replay batch.
parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);

sameToothRescue = struct( ...
  'enable', false, ...
  'coherenceThreshold', config.gatedRescueCoherenceThreshold, ...
  'phaseResidThresholdRad', config.gatedRescuePhaseResidThresholdRad, ...
  'bankMode', "wide-single");
replayNameUse = string(replayName) + "-disabled";
if modeName == "gatedWideSingle"
  sameToothRescue.enable = true;
  replayNameUse = string(replayName) + "-gatedWideSingle";
end

flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti", ...
  'periodicPolishEnableWhenMulti', false, ...
  'sameToothRescue', sameToothRescue));

contextOpt = struct( ...
  'baseSeed', config.baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt);

checkpointOpt = localBuildCheckpointOpt(replayNameUse, config);

replayBatchOpt = struct( ...
  'replayName', replayNameUse, ...
  'snrDb', config.snrDb, ...
  'baseSeed', config.baseSeed, ...
  'seedList', config.seedList, ...
  'numRepeat', config.numRepeat, ...
  'optVerbose', config.optVerbose, ...
  'parallelOpt', parallelOpt, ...
  'contextOpt', contextOpt, ...
  'flowOpt', flowOpt, ...
  'checkpointOpt', checkpointOpt, ...
  'runPeriodicWide', false, ...
  'progressTitle', sprintf('%s %s repeat batch', char(replayName), char(modeName)));
end

function checkpointOpt = localBuildCheckpointOpt(replayNameUse, config)
%LOCALBUILDCHECKPOINTOPT Resolve per-repeat checkpoint storage for one flow mode.
checkpointOpt = struct('enable', config.checkpointEnable);
if ~config.checkpointEnable
  checkpointOpt.runDir = "";
  return;
end
checkpointRunKey = localBuildCheckpointRunKey(config);
checkpointOutputRoot = fullfile(localGetRepoRoot(), 'tmp');
checkpointOpt.resume = logical(config.checkpointResume);
checkpointOpt.runName = string(replayNameUse);
checkpointOpt.runKey = checkpointRunKey;
checkpointOpt.outputRoot = checkpointOutputRoot;
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
checkpointOpt.runDir = fullfile(checkpointOutputRoot, char(replayNameUse), char(checkpointRunKey));
end

function runKey = localBuildCheckpointRunKey(config)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable checkpoint key for the seed and SNR batch.
seedList = reshape(double(config.seedList), 1, []);
runKey = string(sprintf('seed%dto%d_rep%d_snr%.2f', ...
  seedList(1), seedList(end), numel(seedList), config.snrDb));
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
end

function repoRoot = localGetRepoRoot()
%LOCALGETREPOROOT Resolve the repository root from the replay script location.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
end

function checkpointSummary = localBuildCheckpointSummary(disabledRunState, gatedRunState)
%LOCALBUILDCHECKPOINTSUMMARY Store checkpoint status for both replay batches.
checkpointSummary = struct();
checkpointSummary.disabled = localBuildOneCheckpointSummary(disabledRunState);
checkpointSummary.gated = localBuildOneCheckpointSummary(gatedRunState);
end

function checkpointSummary = localBuildOneCheckpointSummary(runState)
%LOCALBUILDONECHECKPOINTSUMMARY Store lightweight checkpoint status for one batch.
checkpointSummary = struct();
checkpointSummary.runName = string(localGetFieldOrDefault(runState, 'runName', ""));
checkpointSummary.runKey = string(localGetFieldOrDefault(runState, 'runKey', ""));
checkpointSummary.runDir = string(localGetFieldOrDefault(runState, 'runDir', ""));
checkpointSummary.numTask = localGetFieldOrDefault(runState, 'numTask', 0);
checkpointSummary.numDone = localGetFieldOrDefault(runState, 'numDone', 0);
checkpointSummary.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
checkpointSummary.cleanedOnSuccess = false;
end

function compareTable = localBuildCompareTable(disabledTable, gatedTable, config)
%LOCALBUILDCOMPARETABLE Align disabled and gated flow rows by task seed.
seedList = disabledTable.taskSeed(:);
numSeed = numel(seedList);
disabledAngleErrDeg = nan(numSeed, 1);
gatedAngleErrDeg = nan(numSeed, 1);
angleGainDeg = nan(numSeed, 1);
disabledToothIdx = nan(numSeed, 1);
gatedToothIdx = nan(numSeed, 1);
disabledFdRefErrHz = nan(numSeed, 1);
gatedFdRefErrHz = nan(numSeed, 1);
disabledCoherenceFloor = nan(numSeed, 1);
gatedCoherenceFloor = nan(numSeed, 1);
rescueTriggered = false(numSeed, 1);
rescueSelected = false(numSeed, 1);
rescueGateReason = strings(numSeed, 1);
disabledFinalTag = strings(numSeed, 1);
gatedFinalTag = strings(numSeed, 1);
caseRole = strings(numSeed, 1);
isHardRescued = false(numSeed, 1);
isEasyDamaged = false(numSeed, 1);
isFdDamaged = false(numSeed, 1);

for iSeed = 1:numSeed
  seed = seedList(iSeed);
  idxDisabled = find(disabledTable.taskSeed == seed, 1, 'first');
  idxGated = find(gatedTable.taskSeed == seed, 1, 'first');
  disabledAngleErrDeg(iSeed) = disabledTable.angleErrDeg(idxDisabled);
  gatedAngleErrDeg(iSeed) = gatedTable.angleErrDeg(idxGated);
  angleGainDeg(iSeed) = disabledAngleErrDeg(iSeed) - gatedAngleErrDeg(iSeed);
  disabledToothIdx(iSeed) = disabledTable.toothIdx(idxDisabled);
  gatedToothIdx(iSeed) = gatedTable.toothIdx(idxGated);
  disabledFdRefErrHz(iSeed) = disabledTable.fdRefErrHz(idxDisabled);
  gatedFdRefErrHz(iSeed) = gatedTable.fdRefErrHz(idxGated);
  disabledCoherenceFloor(iSeed) = disabledTable.nonRefCoherenceFloor(idxDisabled);
  gatedCoherenceFloor(iSeed) = gatedTable.nonRefCoherenceFloor(idxGated);
  rescueTriggered(iSeed) = logical(gatedTable.basinEntryRescueTriggered(idxGated));
  rescueSelected(iSeed) = logical(gatedTable.basinEntryRescueSelected(idxGated));
  rescueGateReason(iSeed) = string(gatedTable.basinEntryGateReason(idxGated));
  disabledFinalTag(iSeed) = string(disabledTable.selectedFinalTag(idxDisabled));
  gatedFinalTag(iSeed) = string(gatedTable.selectedFinalTag(idxGated));
  if disabledAngleErrDeg(iSeed) >= config.hardAngleThresholdDeg
    caseRole(iSeed) = "hard-angle-tail";
  elseif disabledAngleErrDeg(iSeed) <= config.easyAngleThresholdDeg
    caseRole(iSeed) = "easy-angle";
  else
    caseRole(iSeed) = "middle-angle";
  end
  isHardRescued(iSeed) = (caseRole(iSeed) == "hard-angle-tail") && (angleGainDeg(iSeed) >= config.angleGainThresholdDeg);
  isEasyDamaged(iSeed) = (caseRole(iSeed) == "easy-angle") && (-angleGainDeg(iSeed) >= config.damageThresholdDeg);
  isFdDamaged(iSeed) = abs(gatedFdRefErrHz(iSeed)) > max(50, abs(disabledFdRefErrHz(iSeed)) + 50);
end

compareTable = table(seedList, caseRole, disabledAngleErrDeg, gatedAngleErrDeg, angleGainDeg, ...
  disabledToothIdx, gatedToothIdx, disabledFdRefErrHz, gatedFdRefErrHz, ...
  disabledCoherenceFloor, gatedCoherenceFloor, rescueTriggered, rescueSelected, rescueGateReason, ...
  disabledFinalTag, gatedFinalTag, isHardRescued, isEasyDamaged, isFdDamaged, ...
  'VariableNames', {'taskSeed', 'caseRole', 'disabledAngleErrDeg', 'gatedAngleErrDeg', 'angleGainDeg', ...
  'disabledToothIdx', 'gatedToothIdx', 'disabledFdRefErrHz', 'gatedFdRefErrHz', ...
  'disabledCoherenceFloor', 'gatedCoherenceFloor', 'rescueTriggered', 'rescueSelected', 'rescueGateReason', ...
  'disabledFinalTag', 'gatedFinalTag', 'isHardRescued', 'isEasyDamaged', 'isFdDamaged'});
end

function inToothCompareTable = localBuildInToothCompareTable(compareTable)
%LOCALBUILDINTOOTHCOMPARETABLE Keep seeds whose disabled and gated winners stay on the center tooth.
if isempty(compareTable) || height(compareTable) == 0
  inToothCompareTable = compareTable;
  return;
end
mask = (compareTable.disabledToothIdx == 0) & (compareTable.gatedToothIdx == 0);
inToothCompareTable = compareTable(mask, :);
end

function gateMissTable = localBuildGateMissTable(compareTable)
%LOCALBUILDGATEMISSTABLE Collect hard same-tooth misses where the rescue gate did not trigger.
if isempty(compareTable) || height(compareTable) == 0
  gateMissTable = compareTable;
  return;
end
hardMask = compareTable.caseRole == "hard-angle-tail";
sameToothMask = (compareTable.disabledToothIdx == 0) & (compareTable.gatedToothIdx == 0);
missMask = hardMask & sameToothMask & ~compareTable.rescueTriggered;
gateMissTable = compareTable(missMask, :);
end

function aggregateTable = localBuildAggregateTable(compareTable)
%LOCALBUILDAGGREGATETABLE Summarize flow-like rescue trigger, gain, and damage.
numSeed = height(compareTable);
hardMask = compareTable.caseRole == "hard-angle-tail";
easyMask = compareTable.caseRole == "easy-angle";
triggerRate = mean(compareTable.rescueTriggered);
selectedRate = mean(compareTable.rescueSelected);
hardRescueRate = localMeanLogical(compareTable.isHardRescued(hardMask));
easyDamageRate = localMeanLogical(compareTable.isEasyDamaged(easyMask));
fdDamageRate = mean(compareTable.isFdDamaged);
disabledMedianAngleDeg = median(compareTable.disabledAngleErrDeg, 'omitnan');
gatedMedianAngleDeg = median(compareTable.gatedAngleErrDeg, 'omitnan');
disabledP95AngleDeg = prctile(compareTable.disabledAngleErrDeg, 95);
gatedP95AngleDeg = prctile(compareTable.gatedAngleErrDeg, 95);
gatedMaxAngleDeg = max(compareTable.gatedAngleErrDeg, [], 'omitnan');
medianAngleGainDeg = median(compareTable.angleGainDeg, 'omitnan');
aggregateTable = table(numSeed, nnz(hardMask), nnz(easyMask), triggerRate, selectedRate, ...
  hardRescueRate, easyDamageRate, fdDamageRate, disabledMedianAngleDeg, gatedMedianAngleDeg, ...
  disabledP95AngleDeg, gatedP95AngleDeg, gatedMaxAngleDeg, medianAngleGainDeg);
end

function value = localMeanLogical(mask)
%LOCALMEANLOGICAL Return NaN for empty logical groups instead of zero.
if isempty(mask)
  value = NaN;
else
  value = mean(mask);
end
end


function caseRoleSummaryTable = localBuildCaseRoleSummaryTable(compareTable)
%LOCALBUILDCASEROLESUMMARYTABLE Summarize rescue behavior within each offline case role.
if isempty(compareTable) || height(compareTable) == 0
  caseRoleSummaryTable = table();
  return;
end
caseRoleList = unique(compareTable.caseRole, 'stable');
numRole = numel(caseRoleList);
numSeed = zeros(numRole, 1);
sameToothRate = nan(numRole, 1);
triggerRate = nan(numRole, 1);
selectedRate = nan(numRole, 1);
rescueRate = nan(numRole, 1);
damageRate = nan(numRole, 1);
disabledMedianAngleDeg = nan(numRole, 1);
gatedMedianAngleDeg = nan(numRole, 1);
medianAngleGainDeg = nan(numRole, 1);
for iRole = 1:numRole
  mask = compareTable.caseRole == caseRoleList(iRole);
  numSeed(iRole) = nnz(mask);
  sameToothRate(iRole) = mean(compareTable.disabledToothIdx(mask) == 0 & compareTable.gatedToothIdx(mask) == 0);
  triggerRate(iRole) = mean(compareTable.rescueTriggered(mask));
  selectedRate(iRole) = mean(compareTable.rescueSelected(mask));
  rescueRate(iRole) = mean(compareTable.isHardRescued(mask));
  damageRate(iRole) = mean(compareTable.isEasyDamaged(mask) | compareTable.isFdDamaged(mask));
  disabledMedianAngleDeg(iRole) = median(compareTable.disabledAngleErrDeg(mask), 'omitnan');
  gatedMedianAngleDeg(iRole) = median(compareTable.gatedAngleErrDeg(mask), 'omitnan');
  medianAngleGainDeg(iRole) = median(compareTable.angleGainDeg(mask), 'omitnan');
end
caseRoleSummaryTable = table(caseRoleList, numSeed, sameToothRate, triggerRate, selectedRate, ...
  rescueRate, damageRate, disabledMedianAngleDeg, gatedMedianAngleDeg, medianAngleGainDeg, ...
  'VariableNames', {'caseRole', 'numSeed', 'sameToothRate', 'triggerRate', 'selectedRate', ...
  'rescueRate', 'damageRate', 'disabledMedianAngleDeg', 'gatedMedianAngleDeg', 'medianAngleGainDeg'});
end

function gateReasonTable = localBuildGateReasonTable(compareTable)
%LOCALBUILDGATEREASONTABLE Count rescue gate reasons for compact post-load diagnosis.
if isempty(compareTable) || height(compareTable) == 0
  gateReasonTable = table();
  return;
end
gateReasonTable = localBuildCountTable(compareTable, {'rescueGateReason'});
end

function transitionTable = localBuildTransitionCountTable(compareTable, fromField, toField, fromName, toName)
%LOCALBUILDTRANSITIONCOUNTTABLE Count disabled-to-gated transitions for one pair of fields.
if isempty(compareTable) || height(compareTable) == 0
  transitionTable = table();
  return;
end
transitionTable = localBuildCountTable(compareTable, {fromField, toField});
transitionTable.Properties.VariableNames{1} = char(fromName);
transitionTable.Properties.VariableNames{2} = char(toName);
end

function countTable = localBuildCountTable(tableUse, fieldList)
%LOCALBUILDCOUNTTABLE Count unique combinations of one or more table variables.
if isempty(tableUse) || height(tableUse) == 0
  countTable = table();
  return;
end
keyTable = tableUse(:, fieldList);
[countGroup, keyValues] = findgroups(keyTable);
numSeed = splitapply(@numel, ones(height(tableUse), 1), countGroup);
countTable = [keyValues table(numSeed)];
countTable = sortrows(countTable, 'numSeed', 'descend');
end

function storageSummaryTable = localBuildStorageSummaryTable(replayData)
%LOCALBUILDSTORAGESUMMARYTABLE Summarize saved replayData fields for snapshot-size decisions.
fieldList = string(fieldnames(replayData));
numField = numel(fieldList);
fieldClass = strings(numField, 1);
numRow = nan(numField, 1);
numCol = nan(numField, 1);
containsRepeatDetail = false(numField, 1);
for iField = 1:numField
  value = replayData.(fieldList(iField));
  fieldClass(iField) = string(class(value));
  if istable(value)
    numRow(iField) = height(value);
    numCol(iField) = width(value);
  elseif iscell(value)
    numRow(iField) = numel(value);
    numCol(iField) = 1;
  elseif isstruct(value)
    numRow(iField) = numel(value);
    numCol(iField) = numel(fieldnames(value));
  else
    valueSize = size(value);
    if ~isempty(valueSize)
      numRow(iField) = valueSize(1);
      if numel(valueSize) >= 2
        numCol(iField) = valueSize(2);
      end
    end
  end
  containsRepeatDetail(iField) = contains(fieldList(iField), "RepeatCell");
end
storageSummaryTable = table(fieldList, fieldClass, numRow, numCol, containsRepeatDetail, ...
  'VariableNames', {'fieldName', 'fieldClass', 'numRow', 'numCol', 'containsRepeatDetail'});
end

function candidateTable = localBuildCandidateTable(repeatCell)
%LOCALBUILDCANDIDATETABLE Stack gated-flow periodic candidate tables.
rowCell = {};
for iRepeat = 1:numel(repeatCell)
  if ~isfield(repeatCell{iRepeat}, 'msFlow') || ~isfield(repeatCell{iRepeat}.msFlow, 'periodicCandidateTable')
    continue;
  end
  t = repeatCell{iRepeat}.msFlow.periodicCandidateTable;
  if isempty(t)
    continue;
  end
  t.taskSeed = repmat(repeatCell{iRepeat}.taskSeed, height(t), 1);
  t.snrDb = repmat(repeatCell{iRepeat}.snrDb, height(t), 1);
  rowCell{end + 1, 1} = movevars(t, {'taskSeed', 'snrDb'}, 'Before', 1); %#ok<AGROW>
end
if isempty(rowCell)
  candidateTable = table();
else
  candidateTable = vertcat(rowCell{:});
end
end

function warningTable = localBuildWarningTable(disabledTable, gatedTable)
%LOCALBUILDWARNINGTABLE Collect repeats that emitted warnings in either flow.
mask = disabledTable.warningSeen | gatedTable.warningSeen;
if ~any(mask)
  warningTable = table();
  return;
end
warningTable = table(disabledTable.taskSeed(mask), disabledTable.warningId(mask), gatedTable.warningId(mask), ...
  disabledTable.warningMessage(mask), gatedTable.warningMessage(mask), ...
  'VariableNames', {'taskSeed', 'disabledWarningId', 'gatedWarningId', 'disabledWarningMessage', 'gatedWarningMessage'});
end

function contextSummary = localBuildContextSummary(context, disabledFlowOpt, gatedFlowOpt)
%LOCALBUILDCONTEXTSUMMARY Store compact context and flow metadata.
contextSummary = struct();
contextSummary.refSatIdxGlobal = localGetFieldOrDefault(context, 'refSatIdxGlobal', NaN);
contextSummary.otherSatIdxGlobal = localGetFieldOrDefault(context, 'otherSatIdxGlobal', NaN);
contextSummary.disabledSameToothRescueEnable = localGetFieldOrDefault(disabledFlowOpt, 'sameToothRescueEnable', false);
contextSummary.gatedSameToothRescueEnable = localGetFieldOrDefault(gatedFlowOpt, 'sameToothRescueEnable', false);
contextSummary.gatedSameToothRescueBankMode = string(localGetFieldOrDefault(gatedFlowOpt, 'sameToothRescueBankMode', "unknown"));
end

function repeatCell = localStripRepeatCell(repeatCell)
%LOCALSTRIPREPEATCELL Drop heavy static bundles from saved repeat cells.
for iRepeat = 1:numel(repeatCell)
  if isfield(repeatCell{iRepeat}, 'staticBundle')
    repeatCell{iRepeat}.staticBundle = [];
  end
end
end


function checkpointDirText = localBuildCheckpointDirText(config)
%LOCALBUILDCHECKPOINTDIRTEXT Build a compact checkpoint directory summary.
checkpointDirText = "";
if isfield(config, 'checkpointDisabledRunDir') && strlength(string(config.checkpointDisabledRunDir)) > 0
  checkpointDirText = string(config.checkpointDisabledRunDir);
end
if isfield(config, 'checkpointGatedRunDir') && strlength(string(config.checkpointGatedRunDir)) > 0
  if strlength(checkpointDirText) > 0
    checkpointDirText = checkpointDirText + newline + string(config.checkpointGatedRunDir);
  else
    checkpointDirText = string(config.checkpointGatedRunDir);
  end
end
end

function lineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Extract compact key metrics from replayData.
lineList = strings(0, 1);
if isfield(replayData, 'aggregateTable') && height(replayData.aggregateTable) > 0
  t = replayData.aggregateTable;
  lineList(end + 1, 1) = sprintf(['All: trig/sel/hard/easy/fd = <code>%.3f/%.3f/%.3f/%.3f/%.3f</code>, ' ...
    'med <code>%.4g→%.4g deg</code>, p95 <code>%.4g→%.4g deg</code>, max <code>%.4g deg</code>'], ...
    localTableValue(t, 'triggerRate'), localTableValue(t, 'selectedRate'), ...
    localTableValue(t, 'hardRescueRate'), localTableValue(t, 'easyDamageRate'), ...
    localTableValue(t, 'fdDamageRate'), localTableValue(t, 'disabledMedianAngleDeg'), ...
    localTableValue(t, 'gatedMedianAngleDeg'), localTableValue(t, 'disabledP95AngleDeg'), ...
    localTableValue(t, 'gatedP95AngleDeg'), localTableValue(t, 'gatedMaxAngleDeg'));
end
if isfield(replayData, 'inToothAggregateTable') && height(replayData.inToothAggregateTable) > 0
  t = replayData.inToothAggregateTable;
  lineList(end + 1, 1) = sprintf(['In-tooth: n=<code>%d</code>, hard=<code>%d</code>, ' ...
    'trig/sel/hard = <code>%.3f/%.3f/%.3f</code>, med <code>%.4g→%.4g deg</code>'], ...
    round(localTableValue(t, 'numSeed')), round(localTableValue(t, 'Var2')), ...
    localTableValue(t, 'triggerRate'), localTableValue(t, 'selectedRate'), ...
    localTableValue(t, 'hardRescueRate'), localTableValue(t, 'disabledMedianAngleDeg'), ...
    localTableValue(t, 'gatedMedianAngleDeg'));
end
if isfield(replayData, 'gateMissTable') && height(replayData.gateMissTable) > 0
  reasonText = localFormatCountTable(localBuildGateReasonTable(replayData.gateMissTable), 'rescueGateReason');
  lineList(end + 1, 1) = sprintf('Gate-miss hard: <code>%d</code> (%s)', ...
    height(replayData.gateMissTable), reasonText);
else
  lineList(end + 1, 1) = 'Gate-miss hard: <code>0</code>';
end
if isfield(replayData, 'warningTable')
  lineList(end + 1, 1) = sprintf('Warnings: <code>%d</code>', height(replayData.warningTable));
end
end

function value = localTableValue(tableUse, fieldName)
%LOCALTABLEVALUE Read a scalar table value by name with NaN fallback.
value = NaN;
if istable(tableUse) && height(tableUse) > 0 && any(strcmp(tableUse.Properties.VariableNames, fieldName))
  rawValue = tableUse.(fieldName);
  if ~isempty(rawValue)
    value = rawValue(1);
  end
end
end

function reasonText = localFormatCountTable(countTable, reasonField)
%LOCALFORMATCOUNTTABLE Format count table entries for one-line notification.
if isempty(countTable) || height(countTable) == 0
  reasonText = '<code>none</code>';
  return;
end
maxItem = min(height(countTable), 4);
partList = strings(maxItem, 1);
for iItem = 1:maxItem
  partList(iItem) = sprintf('<code>%s=%d</code>', ...
    localHtmlEscape(string(countTable.(reasonField)(iItem))), countTable.numSeed(iItem));
end
if height(countTable) > maxItem
  reasonText = strjoin([partList; sprintf('<code>+%d more</code>', height(countTable) - maxItem)], ', ');
else
  reasonText = strjoin(partList, ', ');
end
end




function safeText = localHtmlEscape(textValue)
%LOCALHTMLESCAPE Escape text for Telegram HTML parse mode.
safeText = string(textValue);
safeText = replace(safeText, '&', '&amp;');
safeText = replace(safeText, '<', '&lt;');
safeText = replace(safeText, '>', '&gt;');
safeText = replace(safeText, '"', '&quot;');
safeText = char(safeText);
end

function localPrintTablePreview(titleText, tableUse, maxRows)
%LOCALPRINTTABLEPREVIEW Print a bounded table preview.
fprintf('\n========== %s ==========\n', titleText);
if isempty(tableUse) || height(tableUse) == 0
  fprintf('  <empty>\n');
  return;
end
if height(tableUse) <= maxRows
  disp(tableUse);
else
  disp(tableUse(1:maxRows, :));
  fprintf('  ... (%d rows total)\n', height(tableUse));
end
end

function localPlotReplay(compareTable)
%LOCALPLOTREPLAY Plot disabled vs gated angle and coherence summaries.
if isempty(compareTable) || height(compareTable) == 0
  return;
end
figure('Name', 'Flow-like gated basin-entry rescue', 'Color', 'w');
subplot(2, 1, 1);
plot(compareTable.taskSeed, compareTable.disabledAngleErrDeg, '-o'); hold on;
plot(compareTable.taskSeed, compareTable.gatedAngleErrDeg, '-s'); grid on;
xlabel('task seed'); ylabel('angle error (deg)');
legend({'disabled', 'gated'}, 'Location', 'best');
title('Final angle error');
subplot(2, 1, 2);
plot(compareTable.taskSeed, compareTable.disabledCoherenceFloor, '-o'); hold on;
plot(compareTable.taskSeed, compareTable.gatedCoherenceFloor, '-s'); grid on;
xlabel('task seed'); ylabel('non-ref coherence floor');
legend({'disabled', 'gated'}, 'Location', 'best');
title('Non-reference coherence floor');
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function text = localNowString(formatText)
%LOCALNOWSTRING Return a compact wall-clock time string.
text = char(datetime('now', 'Format', formatText));
end
