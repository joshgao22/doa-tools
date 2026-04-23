function [bestIdx, isTrusted, reasonText] = selectBestDynamicSubsetSummary(summaryCell)
%SELECTBESTDYNAMICSUBSETSUMMARY Select one subset using tooth-first, truth-free scores.

arguments
  summaryCell cell
end

numCase = numel(summaryCell);
if numCase <= 0
  error('selectBestDynamicSubsetSummary:EmptySummaryCell', ...
    'summaryCell must contain at least one candidate summary.');
end

firstScoreVec = localBuildSubsetSelectionScore(summaryCell{1});
scoreMat = inf(numCase, numel(firstScoreVec));
scoreMat(1, :) = firstScoreVec;
for iCase = 2:numCase
  scoreMat(iCase, :) = localBuildSubsetSelectionScore(summaryCell{iCase});
end

[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
bestIdx = orderIdx(1);
bestSummary = summaryCell{bestIdx};
isTrusted = logical(localGetFieldOrDefault(bestSummary, 'isResolved', false)) && ...
  isfinite(localGetFieldOrDefault(bestSummary, 'finalObj', NaN)) && ...
  isfinite(localGetFieldOrDefault(bestSummary, 'finalResidualNorm', NaN));
reasonText = "selected by tooth-first objective ranking";
end


function scoreVec = localBuildSubsetSelectionScore(summary)
%LOCALBUILDSUBSETSELECTIONSCORE Build one truth-free tooth-selection score.
% Subset schedules primarily serve as tooth selectors. Keep tooth-first
% ordering at the cross-tooth level only. Once multiple candidates already
% land on the same tooth, ranking must avoid truth-aware residual proxies
% and should first use coarse phase-health gates, then objective / residual,
% and only then the lighter quality penalties.

isResolved = logical(localGetFieldOrDefault(summary, 'isResolved', false));
resolvedPenalty = double(~isResolved);

toothIdx = localGetFieldOrDefault(summary, 'toothIdx', NaN);
if isempty(toothIdx)
  toothIdx = NaN;
else
  toothIdx = toothIdx(1);
end
if ~isfinite(toothIdx)
  toothIdxPenalty = inf;
  absToothIdx = inf;
else
  toothIdxPenalty = double(abs(toothIdx) > 0);
  absToothIdx = abs(toothIdx);
end

finalObj = localGetFieldOrDefault(summary, 'finalObj', inf);
if isempty(finalObj)
  finalObj = inf;
else
  finalObj = finalObj(1);
  if ~isfinite(finalObj)
    finalObj = inf;
  end
end
finalResidualNorm = localGetFieldOrDefault(summary, 'finalResidualNorm', inf);
if isempty(finalResidualNorm)
  finalResidualNorm = inf;
else
  finalResidualNorm = finalResidualNorm(1);
  if ~isfinite(finalResidualNorm)
    finalResidualNorm = inf;
  end
end
runTimeMs = localGetFieldOrDefault(summary, 'runTimeMs', inf);
if isempty(runTimeMs)
  runTimeMs = inf;
else
  runTimeMs = runTimeMs(1);
  if ~isfinite(runTimeMs)
    runTimeMs = inf;
  end
end

supportBucket = localBuildFloorBucket(summary, 'nonRefSupportRatioFloor', 0.90);
fitBucket = localBuildFloorBucket(summary, 'nonRefFitRatioFloor', 0.80);
consistencyBucket = localBuildFloorBucket(summary, 'nonRefConsistencyRatioFloor', 0.98);
coherenceBucket = localBuildFloorBucket(summary, 'nonRefCoherenceFloor', 0.98);
phaseRmsBucket = localBuildCeilBucket(summary, 'nonRefRmsPhaseResidRad', 0.01);
phaseMaxBucket = localBuildCeilBucket(summary, 'nonRefMaxAbsPhaseResidRad', 0.02);
negativeBucket = localBuildCeilBucket(summary, 'maxNonRefNegativeProjectionRatio', 0.10);
phaseHealthBucket = supportBucket + fitBucket + consistencyBucket + ...
  coherenceBucket + phaseRmsBucket + phaseMaxBucket + negativeBucket;

supportPenalty = localGetQualityPenalty(summary, 'nonRefSupportRatioFloor');
fitPenalty = localGetQualityPenalty(summary, 'nonRefFitRatioFloor');
consistencyPenalty = localGetQualityPenalty(summary, 'nonRefConsistencyRatioFloor');
coherencePenalty = localGetQualityPenalty(summary, 'nonRefCoherenceFloor');
nonRefRmsPhaseResidRad = localGetFinitePenalty(summary, 'nonRefRmsPhaseResidRad', inf);
nonRefMaxAbsPhaseResidRad = localGetFinitePenalty(summary, 'nonRefMaxAbsPhaseResidRad', inf);
negativePenalty = localGetFinitePenalty(summary, 'maxNonRefNegativeProjectionRatio', 0);
negativePenalty = min(max(negativePenalty, 0), 1);

% Keep tooth-first ordering only at the cross-tooth level. Once several
% candidates have already landed on the same tooth, do not reintroduce
% truth-aware tooth residuals. Use coarse phase-health gates first, then
% objective / residual, and keep runtime only as the final tiebreaker.
scoreVec = [resolvedPenalty, toothIdxPenalty, absToothIdx, phaseHealthBucket, ...
  finalObj, finalResidualNorm, consistencyPenalty, coherencePenalty, ...
  nonRefRmsPhaseResidRad, nonRefMaxAbsPhaseResidRad, supportPenalty, ...
  fitPenalty, negativePenalty, runTimeMs];
end

function value = localBuildFloorBucket(summary, fieldName, minValue)
%LOCALBUILDFLOORBUCKET Build one coarse penalty bucket for lower-bounded metrics.

metricValue = localGetFinitePenalty(summary, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
metricValue = min(max(metricValue, 0), 1);
value = double(metricValue < minValue);
end


function value = localBuildCeilBucket(summary, fieldName, maxValue)
%LOCALBUILDCEILBUCKET Build one coarse penalty bucket for upper-bounded metrics.

metricValue = localGetFinitePenalty(summary, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
value = double(metricValue > maxValue);
end


function value = localGetQualityPenalty(summary, fieldName)
%LOCALGETQUALITYPENALTY Convert one unit-interval quality metric into a score penalty.

metricValue = localGetFieldOrDefault(summary, fieldName, NaN);
if isempty(metricValue)
  metricValue = NaN;
else
  metricValue = metricValue(1);
end
if ~isfinite(metricValue)
  value = 0;
  return;
end
metricValue = min(max(metricValue, 0), 1);
value = 1 - metricValue;
end


function value = localGetFinitePenalty(summary, fieldName, defaultValue)
%LOCALGETFINITEPENALTY Read one finite scalar penalty metric.

value = defaultValue;
if ~isstruct(summary) || ~isfield(summary, fieldName)
  return;
end
rawValue = summary.(fieldName);
if isempty(rawValue)
  return;
end
rawValue = rawValue(1);
if isfinite(rawValue)
  value = rawValue;
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field with default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue(1);
  end
end
end
