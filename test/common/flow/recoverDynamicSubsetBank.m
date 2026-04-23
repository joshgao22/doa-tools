function [subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
  bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, ...
  selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable, ...
  subsetEscalationReason, subsetEscalationDiag, conditionalRandomReason, conditionalRandomDiag] = ...
  recoverDynamicSubsetBank(periodicFixture, subsetFixtureCell, subsetCaseCell, subsetSummaryCell, ...
  subsetRunTimeSec, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  subsetSeedInfo, truth, toothStepHz, runWideBranches, caseDynMsUnknownFdWide, ...
  fdWideUnknownSummary, caseDynMsKnown, knownSummaryForRecovery, ...
  runUnknownCaseFn, buildUnknownSummaryFn)
%RECOVERDYNAMICSUBSETBANK Run one fast subset bank plus lightweight rescue rounds.
% Keep subset ranking / random rescue scheduling out of the bundle entry so
% the main orchestration only wires stages together.

arguments
  periodicFixture (1,1) struct
  subsetFixtureCell cell
  subsetCaseCell cell
  subsetSummaryCell cell
  subsetRunTimeSec (:,1) double
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
  subsetSeedInfo (1,1) struct
  truth (1,1) struct
  toothStepHz (1,1) double
  runWideBranches (1,1) logical
  caseDynMsUnknownFdWide (1,1) struct
  fdWideUnknownSummary (1,1) struct
  caseDynMsKnown (1,1) struct
  knownSummaryForRecovery (1,1) struct
  runUnknownCaseFn (1,1) function_handle
  buildUnknownSummaryFn (1,1) function_handle
end

[bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable] = ...
  localRankSubsetBank(subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
  runWideBranches, caseDynMsUnknownFdWide, fdWideUnknownSummary, caseDynMsKnown, truth, toothStepHz, flowOpt, buildUnknownSummaryFn);

subsetEscalationReason = "not-needed";
subsetEscalationDiag = struct();
conditionalRandomReason = "not-needed";
conditionalRandomDiag = struct();

if logical(localGetFieldOrDefault(flowOpt, 'enableFastSubsetEscalation', false))
  subsetEscalationOpt = struct( ...
    'enableEscalation', true, ...
    'escalateOnToothDisagreement', localGetFieldOrDefault(flowOpt, 'fastEscalateOnToothDisagreement', true), ...
    'maxTrustedToothIdx', localGetFieldOrDefault(flowOpt, 'fastToothIdxThreshold', 0), ...
    'maxTrustedToothResidualHz', localGetFieldOrDefault(flowOpt, 'fastToothResidualHzThreshold', 50), ...
    'maxFdDriftFromKnownHz', localGetFieldOrDefault(flowOpt, 'fastUnknownFdDriftHzThreshold', 500), ...
    'maxDoaDriftFromKnownDeg', localGetFieldOrDefault(flowOpt, 'fastUnknownDoaDriftDegThreshold', 0.003));
  [needSubsetEscalation, subsetEscalationReason, subsetEscalationDiag] = ...
    shouldEscalateDynamicSubsetRecovery(subsetSummaryCell, selectedSubsetSummary, knownSummaryForRecovery, subsetEscalationOpt);
  if needSubsetEscalation
    [subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, didAppendRecovery] = ...
      localAppendRandomSubsetCandidates(periodicFixture, subsetFixtureCell, subsetCaseCell, ...
      subsetSummaryCell, subsetRunTimeSec, pilotWave, carrierFreq, sampleRate, ...
      optVerbose, flowOpt, subsetSeedInfo, truth, toothStepHz, runUnknownCaseFn, buildUnknownSummaryFn);
    if didAppendRecovery
      [bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable] = ...
        localRankSubsetBank(subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
        runWideBranches, caseDynMsUnknownFdWide, fdWideUnknownSummary, caseDynMsKnown, truth, toothStepHz, flowOpt, buildUnknownSummaryFn);
      subsetSelectReason = string(subsetSelectReason) + " | recovery fallback";
    end
  end
end

if logical(localGetFieldOrDefault(flowOpt, 'enableConditionalRandomSubsetRescue', false))
  [needConditionalRandomRescue, conditionalRandomReason, conditionalRandomDiag] = ...
    localShouldAppendConditionalRandomSubsetRescue(subsetFixtureCell, subsetSummaryCell, ...
    selectedSubsetSummary, knownSummaryForRecovery, flowOpt);
  if needConditionalRandomRescue
    randomRescueFlowOpt = flowOpt;
    randomRescueFlowOpt.fastRescueSubsetOffsetCell = {};
    randomRescueFlowOpt.fastRescueSubsetLabelList = strings(0, 1);
    randomRescueFlowOpt.fastNumRandomSubsetTrialFallback = ...
      max(1, round(localGetFieldOrDefault(flowOpt, 'conditionalRandomSubsetTrialCount', 1)));
    [subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, didAppendRandomRescue] = ...
      localAppendRandomSubsetCandidates(periodicFixture, subsetFixtureCell, subsetCaseCell, ...
      subsetSummaryCell, subsetRunTimeSec, pilotWave, carrierFreq, sampleRate, ...
      optVerbose, randomRescueFlowOpt, subsetSeedInfo, truth, toothStepHz, runUnknownCaseFn, buildUnknownSummaryFn);
    if didAppendRandomRescue
      [bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable] = ...
        localRankSubsetBank(subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
        runWideBranches, caseDynMsUnknownFdWide, fdWideUnknownSummary, caseDynMsKnown, truth, toothStepHz, flowOpt, buildUnknownSummaryFn);
      subsetSelectReason = string(subsetSelectReason) + " | conditional random rescue";
    else
      conditionalRandomReason = "random rescue requested but no new candidate was appended";
    end
  end
end
end


function [bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable] = ...
  localRankSubsetBank(subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
  runWideBranches, caseDynMsUnknownFdWide, fdWideUnknownSummary, caseDynMsKnown, truth, toothStepHz, flowOpt, buildUnknownSummaryFn)
%LOCALRANKSUBSETBANK Rank the current subset bank and pick one winner.

if isempty(subsetFixtureCell)
  bestSubsetIdx = NaN;
  subsetIsTrusted = false;
  subsetSelectReason = "no subset candidate";
  selectedSubsetFixture = struct('subsetLabel', "none", 'subsetOffsetIdx', []);
  subsetCandidateTable = table();
  if runWideBranches && localCaseHasUsableEstimate(caseDynMsUnknownFdWide)
    selectedSubsetSummary = fdWideUnknownSummary;
    selectedSubsetCase = caseDynMsUnknownFdWide;
  else
    selectedSubsetSummary = buildUnknownSummaryFn(caseDynMsKnown, truth, toothStepHz);
    selectedSubsetCase = caseDynMsKnown;
  end
  return;
end

[bestSubsetIdx, subsetIsTrusted, subsetSelectReason] = selectBestDynamicSubsetSummary(subsetSummaryCell);
[promotedSubsetIdx, promoteReason] = localMaybePromoteCentralRandomSubset( ...
  subsetFixtureCell, subsetSummaryCell, bestSubsetIdx, flowOpt);
if isfinite(promotedSubsetIdx) && promotedSubsetIdx ~= bestSubsetIdx
  bestSubsetIdx = promotedSubsetIdx;
  subsetSelectReason = string(subsetSelectReason) + " | " + string(promoteReason);
end
selectedSubsetFixture = subsetFixtureCell{bestSubsetIdx};
selectedSubsetSummary = subsetSummaryCell{bestSubsetIdx};
selectedSubsetCase = subsetCaseCell{bestSubsetIdx};
subsetCandidateTable = localBuildSubsetCandidateTable(subsetFixtureCell, subsetSummaryCell, subsetRunTimeSec);
end


function [bestIdx, reasonText] = localMaybePromoteCentralRandomSubset(subsetFixtureCell, subsetSummaryCell, selectedIdx, flowOpt)
%LOCALMAYBEPROMOTECENTRALRANDOMSUBSET Prefer one central random rescue over a weak curated winner.

reasonText = "not-needed";
bestIdx = selectedIdx;
if isempty(subsetSummaryCell) || selectedIdx < 1 || selectedIdx > numel(subsetSummaryCell)
  return;
end
selectedLabel = lower(strtrim(char(string(localGetFieldOrDefault( ...
  subsetFixtureCell{selectedIdx}, 'subsetLabel', "")))));
if startsWith(selectedLabel, 'random') || startsWith(selectedLabel, 'rescue-random')
  return;
end
selectedSummary = subsetSummaryCell{selectedIdx};
promoteResidualHz = localGetFieldOrDefault(flowOpt, 'subsetRandomRescuePromoteResidualHz', 20);
selectedNeedsRescue = localSummaryNeedsTrustedCentralRescue(selectedSummary, promoteResidualHz, inf);
selectedResidualHz = abs(localGetFieldOrDefault(selectedSummary, 'toothResidualHz', inf));
selectedHealthBucket = localBuildSubsetHealthBucket(selectedSummary);

randomIdxList = [];
randomSummaryCell = {};
for iCase = 1:numel(subsetFixtureCell)
  labelUse = lower(strtrim(char(string(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetLabel', "")))));
  if ~(startsWith(labelUse, 'random') || startsWith(labelUse, 'rescue-random'))
    continue;
  end
  summaryUse = subsetSummaryCell{iCase};
  if ~logical(localGetFieldOrDefault(summaryUse, 'isResolved', false))
    continue;
  end
  toothIdxUse = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
  toothResidualHzUse = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
  if ~(isfinite(toothIdxUse) && abs(toothIdxUse) == 0 && ...
      isfinite(toothResidualHzUse) && toothResidualHzUse <= promoteResidualHz)
    continue;
  end
  randomIdxList(end + 1, 1) = iCase; %#ok<AGROW>
  randomSummaryCell{end + 1, 1} = summaryUse; %#ok<AGROW>
end
if isempty(randomIdxList)
  return;
end

[bestRandomLocalIdx, ~] = selectBestDynamicSubsetSummary(randomSummaryCell);
bestRandomSummary = randomSummaryCell{bestRandomLocalIdx};
bestRandomResidualHz = abs(localGetFieldOrDefault(bestRandomSummary, 'toothResidualHz', inf));
bestRandomHealthBucket = localBuildSubsetHealthBucket(bestRandomSummary);
canUseHealthOverride = isfinite(selectedResidualHz) && isfinite(bestRandomResidualHz) && ...
  bestRandomResidualHz <= selectedResidualHz + promoteResidualHz && ...
  bestRandomHealthBucket + 1 < selectedHealthBucket;
if ~(selectedNeedsRescue || canUseHealthOverride)
  return;
end

bestIdx = randomIdxList(bestRandomLocalIdx);
if selectedNeedsRescue
  reasonText = "promote central random rescue subset";
else
  reasonText = "promote central random rescue subset by healthier non-reference fit";
end
end


function [subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, didAppendRandom] = ...
  localAppendRandomSubsetCandidates(periodicFixture, subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, subsetSeedInfo, truth, toothStepHz, runUnknownCaseFn, buildUnknownSummaryFn)
%LOCALAPPENDRANDOMSUBSETCANDIDATES Append one deterministic rescue bank or random fallback.

didAppendRandom = false;
ctx = localGetFieldOrDefault(periodicFixture, 'subsetBuildContext', struct());
if isempty(ctx) || ~isstruct(ctx)
  return;
end
parallelOpt = localGetFieldOrDefault(ctx, 'parallelOpt', struct());

rescueOffsetCell = localGetFieldOrDefault(flowOpt, 'fastRescueSubsetOffsetCell', {});
rescueLabelList = localGetFieldOrDefault(flowOpt, 'fastRescueSubsetLabelList', strings(0, 1));
[rescueOffsetCell, rescueLabelList] = localFilterNewRecoveryLabels(subsetFixtureCell, rescueOffsetCell, rescueLabelList);
if ~isempty(rescueOffsetCell)
  rescueFixtureCell = buildDynamicSubsetFixturesByLabel( ...
    ctx.sceneSeqMaster, ctx.linkParamCellMaster, ctx.rxSigCellMaster, ctx.masterOffsetIdx, ...
    rescueOffsetCell, rescueLabelList, 0, localGetFieldOrDefault(ctx, 'taskSeed', 253), ctx.randomBaseOffsetIdx, ...
    ctx.gridSize, ctx.searchRange, ctx.E, ctx.wavelen, ctx.sampleRate, ...
    ctx.fdRangeDefault, ctx.fdRateRangeDefault, parallelOpt);
  if ~isempty(rescueFixtureCell)
    [newCaseCell, newSummaryCell, newRunTimeSec] = evaluateDynamicSubsetBank( ...
      rescueFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
      periodicFixture.fdRange, periodicFixture.fdRateRange, subsetSeedInfo, truth, toothStepHz, ...
      runUnknownCaseFn, buildUnknownSummaryFn);
    subsetFixtureCell = [subsetFixtureCell; rescueFixtureCell];
    subsetCaseCell = [subsetCaseCell; newCaseCell];
    subsetSummaryCell = [subsetSummaryCell; newSummaryCell];
    subsetRunTimeSec = [subsetRunTimeSec; newRunTimeSec];
    didAppendRandom = true;
    return;
  end
end

numRandomTrial = localGetFieldOrDefault(flowOpt, 'fastNumRandomSubsetTrialFallback', 0);
if ~(isscalar(numRandomTrial) && isfinite(numRandomTrial) && numRandomTrial > 0)
  return;
end
randomSeed = localGetFieldOrDefault(ctx, 'taskSeed', 253) + localGetFieldOrDefault(flowOpt, 'fastRandomSeedOffset', 0);
randomFixtureCell = buildDynamicSubsetFixturesByLabel( ...
  ctx.sceneSeqMaster, ctx.linkParamCellMaster, ctx.rxSigCellMaster, ctx.masterOffsetIdx, ...
  {}, strings(0, 1), numRandomTrial, randomSeed, ctx.randomBaseOffsetIdx, ...
  ctx.gridSize, ctx.searchRange, ctx.E, ctx.wavelen, ctx.sampleRate, ...
  ctx.fdRangeDefault, ctx.fdRateRangeDefault, parallelOpt);
randomFixtureCell = localFilterDuplicateSubsetFixtures(subsetFixtureCell, randomFixtureCell);
if isempty(randomFixtureCell)
  return;
end
[newCaseCell, newSummaryCell, newRunTimeSec] = evaluateDynamicSubsetBank( ...
  randomFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, subsetSeedInfo, truth, toothStepHz, ...
  runUnknownCaseFn, buildUnknownSummaryFn);
subsetFixtureCell = [subsetFixtureCell; randomFixtureCell];
subsetCaseCell = [subsetCaseCell; newCaseCell];
subsetSummaryCell = [subsetSummaryCell; newSummaryCell];
subsetRunTimeSec = [subsetRunTimeSec; newRunTimeSec];
didAppendRandom = true;
end


function [offsetCellOut, labelListOut] = localFilterNewRecoveryLabels(subsetFixtureCell, offsetCellIn, labelListIn)
%LOCALFILTERNEWRECOVERYLABELS Drop recovery candidates whose labels or schedules already exist.

if isempty(offsetCellIn)
  offsetCellOut = {};
  labelListOut = strings(0, 1);
  return;
end
labelListIn = reshape(string(labelListIn), [], 1);
existingLabelList = strings(numel(subsetFixtureCell), 1);
existingOffsetCell = cell(numel(subsetFixtureCell), 1);
for iCase = 1:numel(subsetFixtureCell)
  existingLabelList(iCase) = string(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetLabel', ""));
  existingOffsetCell{iCase} = reshape(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetOffsetIdx', []), 1, []);
end
keepMask = true(numel(offsetCellIn), 1);
for iCase = 1:numel(offsetCellIn)
  if iCase > numel(labelListIn)
    keepMask(iCase) = false;
    continue;
  end
  offsetUse = reshape(offsetCellIn{iCase}, 1, []);
  keepMask(iCase) = ~any(existingLabelList == labelListIn(iCase)) && ...
    ~localHasMatchingSubsetOffsets(existingOffsetCell, offsetUse);
end
offsetCellOut = offsetCellIn(keepMask);
labelListOut = labelListIn(keepMask);
end


function tf = localHasMatchingSubsetOffsets(offsetCell, offsetUse)
%LOCALHASMATCHINGSUBSETOFFSETS Return true when one schedule already exists.

tf = false;
offsetUse = sort(reshape(offsetUse, 1, []));
for iCase = 1:numel(offsetCell)
  existingOffset = reshape(offsetCell{iCase}, 1, []);
  if numel(existingOffset) ~= numel(offsetUse)
    continue;
  end
  if isequal(sort(existingOffset), offsetUse)
    tf = true;
    return;
  end
end
end


function fixtureCellOut = localFilterDuplicateSubsetFixtures(existingFixtureCell, candidateFixtureCell)
%LOCALFILTERDUPLICATESUBSETFIXTURES Drop candidate fixtures whose schedules already exist.

if isempty(candidateFixtureCell)
  fixtureCellOut = candidateFixtureCell;
  return;
end
existingOffsetCell = cell(numel(existingFixtureCell), 1);
for iCase = 1:numel(existingFixtureCell)
  existingOffsetCell{iCase} = reshape(localGetFieldOrDefault(existingFixtureCell{iCase}, 'subsetOffsetIdx', []), 1, []);
end
fixtureCellOut = cell(0, 1);
for iCase = 1:numel(candidateFixtureCell)
  fixtureUse = candidateFixtureCell{iCase};
  offsetUse = reshape(localGetFieldOrDefault(fixtureUse, 'subsetOffsetIdx', []), 1, []);
  if localHasMatchingSubsetOffsets(existingOffsetCell, offsetUse)
    continue;
  end
  existingOffsetCell{end + 1, 1} = offsetUse; %#ok<AGROW>
  fixtureCellOut{end + 1, 1} = fixtureUse; %#ok<AGROW>
end
end


function bucket = localBuildSubsetHealthBucket(summary)
%LOCALBUILDSUBSETHEALTHBUCKET Count coarse non-reference quality failures.

bucket = 0;
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefSupportRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefFitRatioFloor', 1) < 0.95);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefConsistencyRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefCoherenceFloor', 1) < 0.995);
rmsPhaseResid = localGetFieldOrDefault(summary, 'nonRefRmsPhaseResidRad', NaN);
if isfinite(rmsPhaseResid)
  bucket = bucket + double(rmsPhaseResid > 0.003);
end
maxPhaseResid = localGetFieldOrDefault(summary, 'nonRefMaxAbsPhaseResidRad', NaN);
if isfinite(maxPhaseResid)
  bucket = bucket + double(maxPhaseResid > 0.005);
end
negativeRatio = localGetFieldOrDefault(summary, 'maxNonRefNegativeProjectionRatio', NaN);
if isfinite(negativeRatio)
  bucket = bucket + double(negativeRatio > 0.05);
end
end


function [tf, reasonText, diag] = localShouldAppendConditionalRandomSubsetRescue( ...
  subsetFixtureCell, summaryCell, selectedSummary, knownSummary, flowOpt)
%LOCALSHOULDAPPENDCONDITIONALRANDOMSUBSETRESCUE Trigger one cheap random rescue only when curated bank still looks suspicious.

reasonText = "not-needed";
diag = struct();
diag.numCandidate = numel(summaryCell);
diag.numResolved = 0;
diag.numCentralTrusted = 0;
diag.numResolvedOffTooth = 0;
diag.minCentralHealthBucket = inf;
diag.uniqueResolvedToothIdx = [];
diag.hasResolvedToothDisagreement = false;
diag.selectedToothIdx = localGetFieldOrDefault(selectedSummary, 'toothIdx', NaN);
diag.selectedToothResidualHz = abs(localGetFieldOrDefault(selectedSummary, 'toothResidualHz', inf));
diag.selectedHealthBucket = localBuildSubsetHealthBucket(selectedSummary);
diag.fdDriftFromKnownHz = localCalcSummaryFdDrift(selectedSummary, knownSummary);
diag.doaDriftFromKnownDeg = localCalcSummaryDoaDriftDeg(selectedSummary, knownSummary);

tf = false;
if ~logical(localGetFieldOrDefault(flowOpt, 'enableConditionalRandomSubsetRescue', false))
  reasonText = "disabled";
  return;
end
if isempty(summaryCell)
  reasonText = "empty subset bank";
  return;
end
if localHasRandomSubsetLabel(subsetFixtureCell)
  reasonText = "random subset already present";
  return;
end

residualTolHz = localGetFieldOrDefault(flowOpt, 'conditionalRandomToothResidualHzThreshold', 20);
maxHealthBucket = localGetFieldOrDefault(flowOpt, 'conditionalRandomMaxHealthBucket', 1);
resolvedToothIdxList = nan(numel(summaryCell), 1);
for iCase = 1:numel(summaryCell)
  summaryUse = summaryCell{iCase};
  if ~logical(localGetFieldOrDefault(summaryUse, 'isResolved', false))
    continue;
  end
  diag.numResolved = diag.numResolved + 1;
  toothIdx = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
  toothResidualHz = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
  if isfinite(toothIdx)
    resolvedToothIdxList(diag.numResolved) = toothIdx;
  end
  if isfinite(toothIdx) && abs(toothIdx) > 0
    diag.numResolvedOffTooth = diag.numResolvedOffTooth + 1;
  end
  if ~(isfinite(toothIdx) && abs(toothIdx) == 0 && isfinite(toothResidualHz) && toothResidualHz <= residualTolHz)
    continue;
  end
  diag.numCentralTrusted = diag.numCentralTrusted + 1;
  diag.minCentralHealthBucket = min(diag.minCentralHealthBucket, localBuildSubsetHealthBucket(summaryUse));
end
resolvedToothIdxList = resolvedToothIdxList(isfinite(resolvedToothIdxList));
if ~isempty(resolvedToothIdxList)
  diag.uniqueResolvedToothIdx = unique(resolvedToothIdxList).';
  diag.hasResolvedToothDisagreement = numel(diag.uniqueResolvedToothIdx) > 1;
end

selectedNeedsCentralRescue = localSummaryNeedsTrustedCentralRescue(selectedSummary, residualTolHz, maxHealthBucket);
if selectedNeedsCentralRescue
  tf = true;
  reasonText = "selected curated candidate is not a trusted central tooth";
  return;
end
if diag.numCentralTrusted == 0
  tf = true;
  reasonText = "no trusted central-tooth curated candidate";
  return;
end
minResolvedOffToothCount = localGetFieldOrDefault(flowOpt, 'conditionalRandomMinResolvedOffToothCount', 2);
if logical(localGetFieldOrDefault(flowOpt, 'conditionalRandomRescueOnToothDisagreement', true)) && ...
    diag.hasResolvedToothDisagreement && diag.numResolvedOffTooth >= minResolvedOffToothCount
  tf = true;
  reasonText = "resolved curated candidates disagree on tooth";
  return;
end
if diag.selectedHealthBucket > maxHealthBucket
  tf = true;
  reasonText = "selected curated candidate has weak non-reference phase health";
  return;
end
maxFdDriftHz = localGetFieldOrDefault(flowOpt, 'conditionalRandomMaxFdDriftHz', inf);
if isfinite(diag.fdDriftFromKnownHz) && diag.fdDriftFromKnownHz > maxFdDriftHz
  tf = true;
  reasonText = "selected curated candidate drifts too far from CP-K";
  return;
end
maxDoaDriftDeg = localGetFieldOrDefault(flowOpt, 'conditionalRandomMaxDoaDriftDeg', inf);
if isfinite(diag.doaDriftFromKnownDeg) && diag.doaDriftFromKnownDeg > maxDoaDriftDeg
  tf = true;
  reasonText = "selected curated candidate DoA drifts too far from CP-K";
  return;
end
end


function tf = localHasRandomSubsetLabel(subsetFixtureCell)
%LOCALHASRANDOMSUBSETLABEL Return true when one random/rescue-random label already exists.

tf = false;
for iCase = 1:numel(subsetFixtureCell)
  labelUse = lower(strtrim(char(string(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetLabel', "")))));
  if startsWith(labelUse, 'random') || startsWith(labelUse, 'rescue-random')
    tf = true;
    return;
  end
end
end


function tf = localSummaryNeedsTrustedCentralRescue(summaryUse, residualTolHz, maxHealthBucket)
%LOCALSUMMARYNEEDSTRUSTEDCENTRALRESCUE Return true when one summary still looks unsafe.

tf = true;
if ~(isstruct(summaryUse) && ~isempty(summaryUse) && ...
    logical(localGetFieldOrDefault(summaryUse, 'isResolved', false)))
  return;
end
toothIdx = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
if ~(isfinite(toothIdx) && abs(toothIdx) == 0)
  return;
end
toothResidualHz = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
if ~(isfinite(toothResidualHz) && toothResidualHz <= residualTolHz)
  return;
end
if nargin >= 3 && isfinite(maxHealthBucket)
  if localBuildSubsetHealthBucket(summaryUse) > maxHealthBucket
    return;
  end
end
tf = false;
end


function fdDriftHz = localCalcSummaryFdDrift(summaryA, summaryB)
%LOCALCALCSUMMARYFDDRIFT Build one fdRef drift metric between compact summaries.

fdDriftHz = NaN;
fdA = localGetFieldOrDefault(summaryA, 'fdRefEst', NaN);
fdB = localGetFieldOrDefault(summaryB, 'fdRefEst', NaN);
if isfinite(fdA) && isfinite(fdB)
  fdDriftHz = abs(fdA - fdB);
end
end


function doaDriftDeg = localCalcSummaryDoaDriftDeg(summaryA, summaryB)
%LOCALCALCSUMMARYDOADRIFTDEG Build one DoA drift metric between compact summaries.

doaDriftDeg = NaN;
doaA = localGetFieldOrDefault(summaryA, 'doaParamEst', []);
doaB = localGetFieldOrDefault(summaryB, 'doaParamEst', []);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = NaN;
end
end


function subsetTable = localBuildSubsetCandidateTable(subsetFixtureCell, summaryCell, subsetRunTimeSec)
%LOCALBUILDSUBSETCANDIDATETABLE Build one compact subset ranking table.

numCase = numel(summaryCell);
labelOut = strings(numCase, 1);
offsetText = strings(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
runTimeSec = nan(numCase, 1);
isResolved = false(numCase, 1);
for iCase = 1:numCase
  fixtureUse = subsetFixtureCell{iCase};
  s = summaryCell{iCase};
  labelOut(iCase) = string(localGetFieldOrDefault(fixtureUse, 'subsetLabel', sprintf('cand%d', iCase)));
  offsetText(iCase) = string(localFormatIntegerRow(localGetFieldOrDefault(fixtureUse, 'subsetOffsetIdx', [])));
  angleErrDeg(iCase) = localGetFieldOrDefault(s, 'angleErrDeg', NaN);
  fdRefErrHz(iCase) = localGetFieldOrDefault(s, 'fdRefErrHz', NaN);
  fdRateErrHzPerSec(iCase) = localGetFieldOrDefault(s, 'fdRateErrHzPerSec', NaN);
  toothIdx(iCase) = localGetFieldOrDefault(s, 'toothIdx', NaN);
  toothResidualHz(iCase) = localGetFieldOrDefault(s, 'toothResidualHz', NaN);
  runTimeSec(iCase) = subsetRunTimeSec(iCase);
  isResolved(iCase) = localGetFieldOrDefault(s, 'isResolved', false);
end
subsetTable = table(labelOut, offsetText, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeSec, ...
  'VariableNames', {'label', 'offsets', 'isResolved', 'angleErrDeg', ...
  'fdRefErrHz', 'fdRateErrHzPerSec', 'toothIdx', 'toothResidualHz', 'runTimeSec'});
end


function tf = localCaseHasUsableEstimate(caseUse)
%LOCALCASEHASUSABLEESTIMATE Return true when one case contains finite key estimates.

tf = false;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult) || ~isstruct(estResult)
  return;
end
fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
doaParamEst = localGetFieldOrDefault(estResult, 'doaParamEst', nan(2, 1));
tf = isfinite(fdRefEst) && isfinite(fdRateEst) && all(isfinite(doaParamEst(:)));
end


function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row compactly.

if isempty(valueRow)
  textOut = '[]';
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '), ']'];
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
try
  fieldValue = dataStruct.(fieldName);
catch
end
end
