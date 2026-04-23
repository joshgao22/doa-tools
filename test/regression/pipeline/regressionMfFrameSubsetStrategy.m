% Regression report for deterministic frame-subset strategy selection.
% The goal here is not to keep adding more rescue branches into the main
% flow. Instead, we compare a small deterministic subset-strategy bank and
% recommend one default strategy that is both:
%   1) good at tooth selection under the no-truth-leak simple flow;
%   2) simple enough to use as a stable default without repeated manual
%      schedule tweaking.
%
% This regression now uses two practical improvements:
%   1) relaxed hit metrics: being within one tooth period is also tracked,
%      because later periodic in-tooth refinement may still recover the
%      final solution;
%   2) atomic-candidate caching: each repeat only runs the atomic subset
%      schedules once, and bank strategies are reconstructed by replaying
%      the same subset-selection score rule offline.
%
% In addition, results are checkpointed to disk after each completed repeat,
% so long runs can still preserve intermediate outputs.
clear(); close all;

localAddProjectPath();

%% Configuration
snrDbList = [10];
numRepeat = 24;
baseSeed = 253;
optVerbose = false;
windowHitFactor = 1.0;

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
context = buildDynamicDualSatEciContext(struct( ...
  'baseSeed', baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt));
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti", ...
  'periodicPolishEnableWhenMulti', true, ...
  'periodicPolishDoaHalfWidthDeg', [0.002; 0.002]));

atomicStrategyList = localBuildAtomicStrategyList(context.periodicOffsetIdx);
strategyList = localBuildStrategyList(atomicStrategyList);
numAtomic = numel(atomicStrategyList);
numStrategy = numel(strategyList);
numSnr = numel(snrDbList);
numJob = numRepeat * numSnr;

truthAnchor = buildDynamicRepeatData(context, snrDbList(1), baseSeed);
toothStepHz = localResolveToothStepHz(truthAnchor.periodicFixture);
windowTolHz = windowHitFactor * toothStepHz;
saveOpt = localBuildSaveOpt();
runInfo = struct();
runInfo.scriptName = mfilename();
runInfo.runStamp = saveOpt.runStamp;
runInfo.resultDir = saveOpt.resultDir;
runInfo.checkpointFile = saveOpt.checkpointFile;
runInfo.finalResultFile = saveOpt.finalResultFile;
runInfo.toothStepHz = toothStepHz;
runInfo.windowTolHz = windowTolHz;
runInfo.windowHitFactor = windowHitFactor;

fprintf('Running regressionMfFrameSubsetStrategy ...\n');
fprintf('  evaluated strategies          : %d\n', numStrategy);
fprintf('  atomic strategies             : %d\n', numAtomic);
fprintf('  snr list (dB)                 : %s\n', localFormatNumericRow(snrDbList));
fprintf('  repeats per SNR               : %d\n', numRepeat);
fprintf('  tooth step (Hz)               : %.6f\n', toothStepHz);
fprintf('  relaxed window tol (Hz)       : %.6f\n', windowTolHz);
fprintf('  result dir                    : %s\n', saveOpt.resultDir);

resultCellAtomic = cell(numJob, numAtomic);
resultCell = cell(numJob, numStrategy);
jobDoneMask = false(numJob, 1);
tracker = localCreateProgressTracker('Evaluating frame-subset strategies', numJob);
cleanupTracker = onCleanup(@() localCloseProgressTracker(tracker));
for iJob = 1:numJob
  [iSnr, iRepeat] = localDecodeJobIndex(iJob, numRepeat, numSnr);
  snrDb = snrDbList(iSnr);
  taskSeed = baseSeed + (iRepeat - 1);

  repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);
  staticBundle = buildDoaDopplerStaticTransitionBundle( ...
    repeatData.periodicFixture.viewRefOnly, repeatData.periodicFixture.viewOtherOnly, repeatData.periodicFixture.viewMs, ...
    context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
    repeatData.periodicFixture.fdRange, repeatData.truth, context.otherSatIdxGlobal, ...
    optVerbose, flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);

  atomicJobCell = cell(1, numAtomic);
  for iAtomic = 1:numAtomic
    subsetFixtureCell = localBuildSubsetFixtureCellForStrategy(repeatData.periodicFixture, atomicStrategyList(iAtomic), taskSeed);
    msFlow = runSimpleDynamicSubsetPeriodicFlow("MS-MF-Dyn", "multi", ...
      repeatData.periodicFixture, subsetFixtureCell, staticBundle.caseStaticMs, ...
      context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
      optVerbose, flowOpt);
    atomicJobCell{iAtomic} = localBuildStrategyResult( ...
      msFlow, repeatData.truth, atomicStrategyList(iAtomic), snrDb, taskSeed, toothStepHz, windowTolHz);
  end

  strategyJobCell = cell(1, numStrategy);
  for iStrategy = 1:numStrategy
    strategyJobCell{iStrategy} = localAssembleStrategyResultFromAtomic( ...
      atomicJobCell, atomicStrategyList, strategyList(iStrategy), toothStepHz, windowTolHz);
  end

  resultCellAtomic(iJob, :) = atomicJobCell;
  resultCell(iJob, :) = strategyJobCell;
  jobDoneMask(iJob) = true;
  localSaveCheckpoint(saveOpt, runInfo, snrDbList, numRepeat, baseSeed, ...
    atomicStrategyList, strategyList, resultCellAtomic, resultCell, jobDoneMask);
  progressbar('advance');
end
clear cleanupTracker;

summaryTable = localBuildStrategySummaryTable(resultCell, strategyList);
[bestIdx, rankingMat] = localSelectBestStrategy(summaryTable);
summaryTable.rank = rankingMat(:, 1);
summaryTable = movevars(summaryTable, 'rank', 'Before', 'strategyName');
summaryTable = sortrows(summaryTable, 'rank');

bestStrategy = strategyList(bestIdx);
bestSummary = summaryTable(summaryTable.strategyName == bestStrategy.strategyName, :);
periodicMask = summaryTable.strategyName == "periodic10";
if any(periodicMask)
  periodicSummary = summaryTable(periodicMask, :);
else
  periodicSummary = table();
end

fprintf('\n========== Deterministic frame-subset strategy summary ==========%s', newline);
disp(summaryTable);

fprintf('\n========== Recommended default strategy ==========%s', newline);
fprintf('  strategy name                   : %s\n', bestStrategy.strategyName);
fprintf('  strategy type                   : %s\n', bestStrategy.strategyType);
fprintf('  candidate count                 : %d\n', bestStrategy.numCandidate);
fprintf('  subset labels                   : %s\n', strjoin(cellstr(bestStrategy.labelList), ', '));
for iCand = 1:bestStrategy.numCandidate
  fprintf('  offsets %-24s : %s\n', char(bestStrategy.labelList(iCand)), ...
    localFormatIntegerRow(bestStrategy.offsetCell{iCand}));
end
fprintf('  strict subset-hit rate          : %.6f\n', bestSummary.subsetStrictToothHitRate);
fprintf('  relaxed subset-hit rate         : %.6f\n', bestSummary.subsetWindowHitRate);
fprintf('  relaxed subset-coverage rate    : %.6f\n', bestSummary.subsetCoverageWindowHitRate);
fprintf('  strict final-hit rate           : %.6f\n', bestSummary.finalStrictToothHitRate);
fprintf('  relaxed final-hit rate          : %.6f\n', bestSummary.finalWindowHitRate);
fprintf('  relaxed final-coverage rate     : %.6f\n', bestSummary.finalCoverageWindowHitRate);
fprintf('  final angle RMSE (deg)          : %.6f\n', bestSummary.finalAngleRmseDeg);
fprintf('  final angle P95 (deg)           : %.6f\n', bestSummary.finalAngleP95Deg);
fprintf('  final fdRef RMSE (Hz)           : %.6f\n', bestSummary.finalFdRmseHz);
fprintf('  final fdRate RMSE (Hz/s)        : %.6f\n', bestSummary.finalFdRateRmseHzPerSec);
if ~isempty(periodicSummary)
  fprintf('  vs periodic10 relaxed-hit delta : %.6f\n', ...
    bestSummary.finalWindowHitRate - periodicSummary.finalWindowHitRate);
  fprintf('  vs periodic10 strict-hit delta  : %.6f\n', ...
    bestSummary.finalStrictToothHitRate - periodicSummary.finalStrictToothHitRate);
  fprintf('  vs periodic10 angle-RMSE delta  : %.6f deg\n', ...
    periodicSummary.finalAngleRmseDeg - bestSummary.finalAngleRmseDeg);
end

finalPayload = struct();
finalPayload.runInfo = runInfo;
finalPayload.summaryTable = summaryTable;
finalPayload.bestStrategy = bestStrategy;
finalPayload.bestSummary = bestSummary;
finalPayload.periodicSummary = periodicSummary;
finalPayload.atomicStrategyList = atomicStrategyList;
finalPayload.strategyList = strategyList;
finalPayload.resultCellAtomic = resultCellAtomic;
finalPayload.resultCell = resultCell;
finalPayload.jobDoneMask = jobDoneMask;
finalPayload.config = struct( ...
  'snrDbList', snrDbList, ...
  'numRepeat', numRepeat, ...
  'baseSeed', baseSeed, ...
  'windowHitFactor', windowHitFactor, ...
  'parallelOpt', parallelOpt, ...
  'flowOpt', flowOpt);
localSaveFinalResult(saveOpt, finalPayload);

if bestSummary.finalResolvedRate <= 0
  error('regressionMfFrameSubsetStrategy:NoResolvedWinner', ...
    'The recommended strategy did not produce any resolved final cases.');
end
if ~isempty(periodicSummary) && bestSummary.subsetWindowHitRate < periodicSummary.subsetWindowHitRate
  error('regressionMfFrameSubsetStrategy:RecommendedStrategyWorseThanPeriodic', ...
    'The recommended strategy is worse than the periodic baseline on relaxed subset hit rate.');
end

fprintf('  saved final result              : %s\n', saveOpt.finalResultFile);
fprintf('PASS: regressionMfFrameSubsetStrategy\n');


function atomicStrategyList = localBuildAtomicStrategyList(periodicOffsetIdx)
[primaryOffsetCell, primaryLabelList, rescueOffsetCell, rescueLabelList] = getDynamicCuratedSubsetBank();
allOffsetCell = [primaryOffsetCell, rescueOffsetCell];
allLabelList = [primaryLabelList; rescueLabelList];

atomicStrategyList = repmat(struct( ...
  'strategyName', "", ...
  'strategyType', "", ...
  'offsetCell', {{}}, ...
  'labelList', strings(0, 1), ...
  'memberIdx', zeros(0, 1), ...
  'numCandidate', 0), 0, 1);

atomicStrategyList(end + 1, 1) = localMakeStrategy("periodic10", "single", ...
  {reshape(periodicOffsetIdx, 1, [])}, "periodic10"); %#ok<AGROW>
for iCurated = 1:numel(allOffsetCell)
  atomicStrategyList(end + 1, 1) = localMakeStrategy(allLabelList(iCurated), "single", ...
    {reshape(allOffsetCell{iCurated}, 1, [])}, allLabelList(iCurated)); %#ok<AGROW>
end
end

function strategyList = localBuildStrategyList(atomicStrategyList)
strategyList = atomicStrategyList;
labelMap = containers.Map();
for iAtomic = 1:numel(atomicStrategyList)
  labelMap(char(atomicStrategyList(iAtomic).strategyName)) = iAtomic;
end

strategyList(end + 1, 1) = localBuildBankStrategy("curated12", atomicStrategyList, labelMap, ["curated1", "curated2"]); %#ok<AGROW>
strategyList(end + 1, 1) = localBuildBankStrategy("curated123", atomicStrategyList, labelMap, ["curated1", "curated2", "curated3"]); %#ok<AGROW>
strategyList(end + 1, 1) = localBuildBankStrategy("curated124", atomicStrategyList, labelMap, ["curated1", "curated2", "curated4"]); %#ok<AGROW>
strategyList(end + 1, 1) = localBuildBankStrategy("curated34", atomicStrategyList, labelMap, ["curated3", "curated4"]); %#ok<AGROW>
strategyList(end + 1, 1) = localBuildBankStrategy("curated1234", atomicStrategyList, labelMap, ["curated1", "curated2", "curated3", "curated4"]); %#ok<AGROW>
end

function strategyUse = localBuildBankStrategy(strategyName, atomicStrategyList, labelMap, memberNames)
memberIdx = nan(numel(memberNames), 1);
offsetCell = cell(numel(memberNames), 1);
labelList = strings(numel(memberNames), 1);
for iMember = 1:numel(memberNames)
  memberIdx(iMember) = labelMap(char(memberNames(iMember)));
  offsetCell{iMember} = atomicStrategyList(memberIdx(iMember)).offsetCell{1};
  labelList(iMember) = atomicStrategyList(memberIdx(iMember)).labelList(1);
end
strategyUse = struct();
strategyUse.strategyName = string(strategyName);
strategyUse.strategyType = "bank";
strategyUse.offsetCell = offsetCell;
strategyUse.labelList = labelList;
strategyUse.memberIdx = memberIdx;
strategyUse.numCandidate = numel(memberIdx);
end

function strategyUse = localMakeStrategy(strategyName, strategyType, offsetCell, labelList)
strategyUse = struct();
strategyUse.strategyName = string(strategyName);
strategyUse.strategyType = string(strategyType);
strategyUse.offsetCell = offsetCell;
strategyUse.labelList = reshape(string(labelList), [], 1);
strategyUse.memberIdx = []; %#ok<STRNU>
strategyUse.numCandidate = numel(offsetCell);
end

function subsetFixtureCell = localBuildSubsetFixtureCellForStrategy(periodicFixture, strategyUse, taskSeed)
ctx = periodicFixture.subsetBuildContext;
parallelOpt = localMergeStruct(localGetFieldOrDefault(ctx, 'parallelOpt', struct()), struct( ...
  'enableSubsetEvalParfor', false, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false));
subsetFixtureCell = buildDynamicSubsetFixturesByLabel( ...
  ctx.sceneSeqMaster, ctx.linkParamCellMaster, ctx.rxSigCellMaster, ctx.masterOffsetIdx, ...
  strategyUse.offsetCell, strategyUse.labelList, 0, taskSeed, ctx.randomBaseOffsetIdx, ...
  ctx.gridSize, ctx.searchRange, ctx.E, ctx.wavelen, ctx.sampleRate, ...
  ctx.fdRangeDefault, ctx.fdRateRangeDefault, parallelOpt);
end

function resultUse = localBuildStrategyResult(msFlow, truth, strategyUse, snrDb, taskSeed, toothStepHz, windowTolHz)
selectedSummary = msFlow.selectedSubsetSummary;
finalSummary = msFlow.finalSummary;
resultUse = struct();
resultUse.strategyName = strategyUse.strategyName;
resultUse.strategyType = strategyUse.strategyType;
resultUse.numCandidate = strategyUse.numCandidate;
resultUse.snrDb = snrDb;
resultUse.taskSeed = taskSeed;
resultUse.selectedSubsetLabel = string(msFlow.selectedSubsetLabel);
resultUse.subsetToothIdx = localGetFieldOrDefault(selectedSummary, 'toothIdx', NaN);
resultUse.subsetToothResidualHz = abs(localGetFieldOrDefault(selectedSummary, 'toothResidualHz', NaN));
resultUse.subsetAngleErrDeg = localGetFieldOrDefault(selectedSummary, 'angleErrDeg', NaN);
resultUse.subsetFdRefErrHz = localGetFieldOrDefault(selectedSummary, 'fdRefErrHz', NaN);
resultUse.subsetFdRateErrHzPerSec = localGetFieldOrDefault(selectedSummary, 'fdRateErrHzPerSec', NaN);
resultUse.subsetResolved = logical(localGetFieldOrDefault(selectedSummary, 'isResolved', false));
resultUse.subsetStrictHit = localIsCorrectTooth(resultUse.subsetToothIdx, resultUse.subsetResolved);
resultUse.subsetWindowHit = localIsWithinWindow(resultUse.subsetFdRefErrHz, resultUse.subsetResolved, windowTolHz);
resultUse.finalToothIdx = localGetFieldOrDefault(finalSummary, 'toothIdx', NaN);
resultUse.finalToothResidualHz = abs(localGetFieldOrDefault(finalSummary, 'toothResidualHz', NaN));
resultUse.finalAngleErrDeg = localGetFieldOrDefault(finalSummary, 'angleErrDeg', NaN);
resultUse.finalFdRefErrHz = localGetFieldOrDefault(finalSummary, 'fdRefErrHz', NaN);
resultUse.finalFdRateErrHzPerSec = localGetFieldOrDefault(finalSummary, 'fdRateErrHzPerSec', NaN);
resultUse.finalResolved = logical(localGetFieldOrDefault(finalSummary, 'isResolved', false));
resultUse.finalStrictHit = localIsCorrectTooth(resultUse.finalToothIdx, resultUse.finalResolved);
resultUse.finalWindowHit = localIsWithinWindow(resultUse.finalFdRefErrHz, resultUse.finalResolved, windowTolHz);
resultUse.finalObj = localGetFieldOrDefault(finalSummary, 'finalObj', NaN);
resultUse.finalResidualNorm = localGetFieldOrDefault(finalSummary, 'finalResidualNorm', NaN);
resultUse.periodicReplaySeedSource = string(localGetFieldOrDefault(msFlow.periodicDoaSeed, 'selectedReplaySeedSource', "unknown"));
resultUse.periodicReplaySelectionReason = string(localGetFieldOrDefault(msFlow.periodicDoaSeed, 'selectedReplaySelectionReason', "unknown"));
resultUse.periodicPolishGateReason = string(localGetFieldOrDefault(msFlow.periodicDoaSeed, 'polishGateReason', "unknown"));
resultUse.toothStepHz = toothStepHz;
resultUse.windowTolHz = windowTolHz;
resultUse.truthLatlonDeg = reshape(truth.latlonTrueDeg, 1, []);
resultUse.subsetScoreVec = buildDynamicSelectionScoreVector(msFlow.selectedSubsetSummary);
resultUse.atomicWinnerLabel = strategyUse.strategyName;
resultUse.coverageStrictSubset = resultUse.subsetStrictHit;
resultUse.coverageWindowSubset = resultUse.subsetWindowHit;
resultUse.coverageStrictFinal = resultUse.finalStrictHit;
resultUse.coverageWindowFinal = resultUse.finalWindowHit;
resultUse.selectionReason = "single-candidate";
end

function resultUse = localAssembleStrategyResultFromAtomic(atomicJobCell, atomicStrategyList, strategyUse, toothStepHz, windowTolHz)
if strategyUse.strategyType == "single"
  memberNames = strategyUse.labelList(1);
  memberIdx = localFindAtomicMemberIdx(atomicStrategyList, memberNames);
else
  memberIdx = strategyUse.memberIdx;
end
memberResultCell = atomicJobCell(memberIdx);
scoreMat = nan(numel(memberIdx), numel(memberResultCell{1}.subsetScoreVec));
for iMember = 1:numel(memberIdx)
  scoreMat(iMember, :) = reshape(memberResultCell{iMember}.subsetScoreVec, 1, []);
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
winnerLocalIdx = orderIdx(1);
winnerResult = memberResultCell{winnerLocalIdx};

resultUse = winnerResult;
resultUse.strategyName = strategyUse.strategyName;
resultUse.strategyType = strategyUse.strategyType;
resultUse.numCandidate = strategyUse.numCandidate;
resultUse.toothStepHz = toothStepHz;
resultUse.windowTolHz = windowTolHz;
resultUse.atomicWinnerLabel = string(winnerResult.strategyName);
resultUse.selectionReason = localResolveSelectionReason(strategyUse, winnerResult);
resultUse.coverageStrictSubset = any(cellfun(@(r) logical(localGetFieldOrDefault(r, 'subsetStrictHit', false)), memberResultCell));
resultUse.coverageWindowSubset = any(cellfun(@(r) logical(localGetFieldOrDefault(r, 'subsetWindowHit', false)), memberResultCell));
resultUse.coverageStrictFinal = any(cellfun(@(r) logical(localGetFieldOrDefault(r, 'finalStrictHit', false)), memberResultCell));
resultUse.coverageWindowFinal = any(cellfun(@(r) logical(localGetFieldOrDefault(r, 'finalWindowHit', false)), memberResultCell));
end

function summaryTable = localBuildStrategySummaryTable(resultCell, strategyList)
numStrategy = numel(strategyList);
strategyName = strings(numStrategy, 1);
strategyType = strings(numStrategy, 1);
numCandidate = nan(numStrategy, 1);
subsetStrictToothHitRate = nan(numStrategy, 1);
subsetWindowHitRate = nan(numStrategy, 1);
subsetCoverageStrictHitRate = nan(numStrategy, 1);
subsetCoverageWindowHitRate = nan(numStrategy, 1);
subsetResolvedRate = nan(numStrategy, 1);
subsetAngleRmseDeg = nan(numStrategy, 1);
subsetFdRmseHz = nan(numStrategy, 1);
subsetFdRateRmseHzPerSec = nan(numStrategy, 1);
finalStrictToothHitRate = nan(numStrategy, 1);
finalWindowHitRate = nan(numStrategy, 1);
finalCoverageStrictHitRate = nan(numStrategy, 1);
finalCoverageWindowHitRate = nan(numStrategy, 1);
finalResolvedRate = nan(numStrategy, 1);
finalAngleRmseDeg = nan(numStrategy, 1);
finalAngleP95Deg = nan(numStrategy, 1);
finalFdRmseHz = nan(numStrategy, 1);
finalFdP95Hz = nan(numStrategy, 1);
finalFdRateRmseHzPerSec = nan(numStrategy, 1);
finalFdRateP95HzPerSec = nan(numStrategy, 1);
meanFinalObj = nan(numStrategy, 1);
rankingLossSubsetWindowRate = nan(numStrategy, 1);
rankingLossFinalWindowRate = nan(numStrategy, 1);

for iStrategy = 1:numStrategy
  strategyName(iStrategy) = strategyList(iStrategy).strategyName;
  strategyType(iStrategy) = strategyList(iStrategy).strategyType;
  numCandidate(iStrategy) = strategyList(iStrategy).numCandidate;
  resultVec = resultCell(:, iStrategy);
  subsetStrictHit = cellfun(@(r) logical(localGetFieldOrDefault(r, 'subsetStrictHit', false)), resultVec);
  subsetWindowHit = cellfun(@(r) logical(localGetFieldOrDefault(r, 'subsetWindowHit', false)), resultVec);
  subsetCoverageStrict = cellfun(@(r) logical(localGetFieldOrDefault(r, 'coverageStrictSubset', false)), resultVec);
  subsetCoverageWindow = cellfun(@(r) logical(localGetFieldOrDefault(r, 'coverageWindowSubset', false)), resultVec);
  finalStrictHit = cellfun(@(r) logical(localGetFieldOrDefault(r, 'finalStrictHit', false)), resultVec);
  finalWindowHit = cellfun(@(r) logical(localGetFieldOrDefault(r, 'finalWindowHit', false)), resultVec);
  finalCoverageStrict = cellfun(@(r) logical(localGetFieldOrDefault(r, 'coverageStrictFinal', false)), resultVec);
  finalCoverageWindow = cellfun(@(r) logical(localGetFieldOrDefault(r, 'coverageWindowFinal', false)), resultVec);
  subsetResolved = cellfun(@(r) logical(localGetFieldOrDefault(r, 'subsetResolved', false)), resultVec);
  finalResolved = cellfun(@(r) logical(localGetFieldOrDefault(r, 'finalResolved', false)), resultVec);

  subsetAngleVec = cellfun(@(r) localGetFieldOrDefault(r, 'subsetAngleErrDeg', NaN), resultVec);
  subsetFdVec = cellfun(@(r) localGetFieldOrDefault(r, 'subsetFdRefErrHz', NaN), resultVec);
  subsetFdRateVec = cellfun(@(r) localGetFieldOrDefault(r, 'subsetFdRateErrHzPerSec', NaN), resultVec);
  finalAngleVec = cellfun(@(r) localGetFieldOrDefault(r, 'finalAngleErrDeg', NaN), resultVec);
  finalFdVec = cellfun(@(r) localGetFieldOrDefault(r, 'finalFdRefErrHz', NaN), resultVec);
  finalFdRateVec = cellfun(@(r) localGetFieldOrDefault(r, 'finalFdRateErrHzPerSec', NaN), resultVec);
  finalObjVec = cellfun(@(r) localGetFieldOrDefault(r, 'finalObj', NaN), resultVec);

  subsetStrictToothHitRate(iStrategy) = mean(subsetStrictHit);
  subsetWindowHitRate(iStrategy) = mean(subsetWindowHit);
  subsetCoverageStrictHitRate(iStrategy) = mean(subsetCoverageStrict);
  subsetCoverageWindowHitRate(iStrategy) = mean(subsetCoverageWindow);
  subsetResolvedRate(iStrategy) = mean(subsetResolved);
  subsetAngleRmseDeg(iStrategy) = localCalcRmse(subsetAngleVec, ~subsetResolved);
  subsetFdRmseHz(iStrategy) = localCalcRmse(subsetFdVec, ~subsetResolved);
  subsetFdRateRmseHzPerSec(iStrategy) = localCalcRmse(subsetFdRateVec, ~subsetResolved);
  finalStrictToothHitRate(iStrategy) = mean(finalStrictHit);
  finalWindowHitRate(iStrategy) = mean(finalWindowHit);
  finalCoverageStrictHitRate(iStrategy) = mean(finalCoverageStrict);
  finalCoverageWindowHitRate(iStrategy) = mean(finalCoverageWindow);
  finalResolvedRate(iStrategy) = mean(finalResolved);
  finalAngleRmseDeg(iStrategy) = localCalcRmse(finalAngleVec, ~finalResolved);
  finalAngleP95Deg(iStrategy) = localCalcP95(finalAngleVec, ~finalResolved);
  finalFdRmseHz(iStrategy) = localCalcRmse(finalFdVec, ~finalResolved);
  finalFdP95Hz(iStrategy) = localCalcP95(finalFdVec, ~finalResolved);
  finalFdRateRmseHzPerSec(iStrategy) = localCalcRmse(finalFdRateVec, ~finalResolved);
  finalFdRateP95HzPerSec(iStrategy) = localCalcP95(finalFdRateVec, ~finalResolved);
  meanFinalObj(iStrategy) = mean(finalObjVec(finalResolved & isfinite(finalObjVec)), 'omitnan');
  rankingLossSubsetWindowRate(iStrategy) = mean(subsetCoverageWindow & ~subsetWindowHit);
  rankingLossFinalWindowRate(iStrategy) = mean(finalCoverageWindow & ~finalWindowHit);
end

summaryTable = table(strategyName, strategyType, numCandidate, ...
  subsetStrictToothHitRate, subsetWindowHitRate, subsetCoverageStrictHitRate, subsetCoverageWindowHitRate, ...
  rankingLossSubsetWindowRate, subsetResolvedRate, subsetAngleRmseDeg, subsetFdRmseHz, subsetFdRateRmseHzPerSec, ...
  finalStrictToothHitRate, finalWindowHitRate, finalCoverageStrictHitRate, finalCoverageWindowHitRate, ...
  rankingLossFinalWindowRate, finalResolvedRate, finalAngleRmseDeg, finalAngleP95Deg, ...
  finalFdRmseHz, finalFdP95Hz, finalFdRateRmseHzPerSec, finalFdRateP95HzPerSec, meanFinalObj);
end

function [bestIdx, rankingMat] = localSelectBestStrategy(summaryTable)
numStrategy = height(summaryTable);
rankingMat = nan(numStrategy, 6);
rankingMat(:, 1) = (1:numStrategy).';
scoreMat = [ ...
  -summaryTable.finalWindowHitRate, ...
  -summaryTable.finalStrictToothHitRate, ...
  summaryTable.finalFdP95Hz, ...
  summaryTable.finalAngleP95Deg, ...
  summaryTable.numCandidate, ...
  summaryTable.finalFdRateP95HzPerSec];
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
rankingMat(orderIdx, 1) = (1:numStrategy).';
bestIdx = orderIdx(1);
end

function tf = localIsCorrectTooth(toothIdx, isResolved)
tf = false;
if ~logical(isResolved)
  return;
end
if isempty(toothIdx) || ~isfinite(toothIdx)
  return;
end
tf = (abs(round(toothIdx)) == 0);
end

function tf = localIsWithinWindow(fdRefErrHz, isResolved, windowTolHz)
tf = false;
if ~logical(isResolved)
  return;
end
if isempty(fdRefErrHz) || ~isfinite(fdRefErrHz) || isempty(windowTolHz) || ~isfinite(windowTolHz)
  return;
end
tf = (abs(fdRefErrHz) <= windowTolHz);
end

function rmseValue = localCalcRmse(valueVec, invalidMask)
valueUse = reshape(valueVec, [], 1);
maskUse = isfinite(valueUse);
if nargin >= 2 && ~isempty(invalidMask)
  maskUse = maskUse & ~reshape(logical(invalidMask), [], 1);
end
if ~any(maskUse)
  rmseValue = NaN;
  return;
end
valueUse = valueUse(maskUse);
rmseValue = sqrt(mean(valueUse .^ 2));
end

function p95Value = localCalcP95(valueVec, invalidMask)
valueUse = reshape(valueVec, [], 1);
maskUse = isfinite(valueUse);
if nargin >= 2 && ~isempty(invalidMask)
  maskUse = maskUse & ~reshape(logical(invalidMask), [], 1);
end
if ~any(maskUse)
  p95Value = NaN;
  return;
end
valueUse = abs(valueUse(maskUse));
p95Value = prctile(valueUse, 95);
end

function [iSnr, iRepeat] = localDecodeJobIndex(iJob, numRepeat, numSnr)
iSnr = ceil(iJob / numRepeat);
iRepeat = iJob - (iSnr - 1) * numRepeat;
if iSnr < 1 || iSnr > numSnr || iRepeat < 1 || iRepeat > numRepeat
  error('regressionMfFrameSubsetStrategy:InvalidJobIndex', ...
    'Decoded job index is out of range.');
end
end

function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end

function tracker = localCreateProgressTracker(titleText, totalCount)
tracker = struct();
tracker.titleText = char(titleText);
tracker.totalCount = totalCount;

fprintf('%s\n', tracker.titleText);
progressbar('displaymode', 'replace');
progressbar('minimalupdateinterval', 0.2);
progressbar('reset', tracker.totalCount);
end

function localCloseProgressTracker(~)
pause(0.05);
progressbar('end');
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function outStruct = localMergeStruct(baseStruct, overrideStruct)
outStruct = baseStruct;
if nargin < 2 || isempty(overrideStruct)
  return;
end
fieldList = fieldnames(overrideStruct);
for iField = 1:numel(fieldList)
  outStruct.(fieldList{iField}) = overrideStruct.(fieldList{iField});
end
end

function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
if isempty(valueVec)
  textOut = '[]';
  return;
end
cellText = arrayfun(@(x) sprintf('%.6f', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end

function textOut = localFormatIntegerRow(valueVec)
valueVec = reshape(valueVec, 1, []);
if isempty(valueVec)
  textOut = '[]';
  return;
end
cellText = arrayfun(@(x) sprintf('%d', round(x)), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ' '), ']'];
end

function toothStepHz = localResolveToothStepHz(periodicFixture)
toothStepHz = NaN;
if isstruct(periodicFixture) && isfield(periodicFixture, 'frameIntvlSec')
  rawValue = periodicFixture.frameIntvlSec;
  if ~isempty(rawValue) && isfinite(rawValue) && (rawValue > 0)
    toothStepHz = 1 / rawValue;
    return;
  end
end
if isstruct(periodicFixture) && isfield(periodicFixture, 'sceneSeq')
  sceneSeq = periodicFixture.sceneSeq;
  if isstruct(sceneSeq) && isfield(sceneSeq, 'frameIntvlSec')
    rawValue = sceneSeq.frameIntvlSec;
    if ~isempty(rawValue) && isfinite(rawValue) && (rawValue > 0)
      toothStepHz = 1 / rawValue;
      return;
    end
  end
  if isstruct(sceneSeq) && isfield(sceneSeq, 'timeOffsetSec')
    timeOffsetSec = reshape(sceneSeq.timeOffsetSec, 1, []);
    dt = diff(timeOffsetSec);
    dt = dt(isfinite(dt) & (dt > 0));
    if ~isempty(dt)
      rawValue = median(dt);
      if isfinite(rawValue) && (rawValue > 0)
        toothStepHz = 1 / rawValue;
        return;
      end
    end
  end
end
error('regressionMfFrameSubsetStrategy:MissingFrameInterval', ...
  'Cannot resolve frame interval from periodicFixture.');
end

function saveOpt = localBuildSaveOpt()
scriptDir = fileparts(mfilename('fullpath'));
runStamp = datestr(now, 'yyyymmdd-HHMMSS');
resultDir = fullfile(scriptDir, 'results', 'regressionMfFrameSubsetStrategy', runStamp);
if ~exist(resultDir, 'dir')
  mkdir(resultDir);
end
saveOpt = struct();
saveOpt.runStamp = runStamp;
saveOpt.resultDir = resultDir;
saveOpt.checkpointFile = fullfile(resultDir, 'checkpoint.mat');
saveOpt.finalResultFile = fullfile(resultDir, 'result.mat');
saveOpt.summaryCsvFile = fullfile(resultDir, 'summary.csv');
end

function localSaveCheckpoint(saveOpt, runInfo, snrDbList, numRepeat, baseSeed, atomicStrategyList, strategyList, resultCellAtomic, resultCell, jobDoneMask)
checkpoint = struct();
checkpoint.runInfo = runInfo;
checkpoint.config = struct('snrDbList', snrDbList, 'numRepeat', numRepeat, 'baseSeed', baseSeed);
checkpoint.atomicStrategyList = atomicStrategyList;
checkpoint.strategyList = strategyList;
checkpoint.resultCellAtomic = resultCellAtomic;
checkpoint.resultCell = resultCell;
checkpoint.jobDoneMask = jobDoneMask;
save(saveOpt.checkpointFile, '-struct', 'checkpoint');
end

function localSaveFinalResult(saveOpt, finalPayload)
save(saveOpt.finalResultFile, '-struct', 'finalPayload');
try
  writetable(finalPayload.summaryTable, saveOpt.summaryCsvFile);
catch
end
end

function memberIdx = localFindAtomicMemberIdx(atomicStrategyList, memberNames)
memberIdx = nan(numel(memberNames), 1);
for iMember = 1:numel(memberNames)
  matchMask = arrayfun(@(s) s.strategyName == string(memberNames(iMember)), atomicStrategyList);
  if ~any(matchMask)
    error('regressionMfFrameSubsetStrategy:MissingAtomicMember', ...
      'Cannot find atomic strategy "%s".', string(memberNames(iMember)));
  end
  memberIdx(iMember) = find(matchMask, 1, 'first');
end
end

function textOut = localResolveSelectionReason(strategyUse, winnerResult)
if strategyUse.strategyType == "single"
  textOut = "single-candidate";
elseif string(winnerResult.strategyName) == string(winnerResult.atomicWinnerLabel)
  textOut = "best-atomic-score";
else
  textOut = "bank-best-score";
end
end
