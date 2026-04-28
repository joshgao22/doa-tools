function [bestIdx, isTrusted, reasonText] = selectBestDynamicSubsetSummary(summaryCell)
%SELECTBESTDYNAMICSUBSETSUMMARY Select one subset using no-truth scores.

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
  isfinite(localGetFieldOrDefault(bestSummary, 'finalResidualNorm', NaN)) && ...
  localBuildPhaseHealthBucket(bestSummary) <= 0;
reasonText = "selected by no-truth health/objective ranking";
end


function scoreVec = localBuildSubsetSelectionScore(summary)
%LOCALBUILDSUBSETSELECTIONSCORE Build one truth-free subset-selection score.
% Subset schedules are compared with estimator-internal diagnostics only.
% Truth-relative tooth indices or residuals are intentionally ignored here;
% those fields are reserved for replay / scan evaluation tables.

isResolved = logical(localGetFieldOrDefault(summary, 'isResolved', false));
resolvedPenalty = double(~isResolved);

finalObj = localGetFiniteScalar(summary, 'finalObj', inf);
finalResidualNorm = localGetFiniteScalar(summary, 'finalResidualNorm', inf);
runTimeMs = localGetFiniteScalar(summary, 'runTimeMs', inf);

phaseHealthBucket = localBuildPhaseHealthBucket(summary);
supportPenalty = localGetQualityPenalty(summary, 'nonRefSupportRatioFloor');
fitPenalty = localGetQualityPenalty(summary, 'nonRefFitRatioFloor');
consistencyPenalty = localGetQualityPenalty(summary, 'nonRefConsistencyRatioFloor');
coherencePenalty = localGetQualityPenalty(summary, 'nonRefCoherenceFloor');
nonRefRmsPhaseResidRad = localGetFinitePenalty(summary, 'nonRefRmsPhaseResidRad', inf);
nonRefMaxAbsPhaseResidRad = localGetFinitePenalty(summary, 'nonRefMaxAbsPhaseResidRad', inf);
negativePenalty = localGetFinitePenalty(summary, 'maxNonRefNegativeProjectionRatio', 0);
negativePenalty = min(max(negativePenalty, 0), 1);

scoreVec = [resolvedPenalty, phaseHealthBucket, finalObj, finalResidualNorm, ...
  consistencyPenalty, coherencePenalty, nonRefRmsPhaseResidRad, ...
  nonRefMaxAbsPhaseResidRad, supportPenalty, fitPenalty, negativePenalty, ...
  runTimeMs];
end


function bucket = localBuildPhaseHealthBucket(summary)
%LOCALBUILDPHASEHEALTHBUCKET Count coarse non-reference phase-health failures.

supportBucket = localBuildFloorBucket(summary, 'nonRefSupportRatioFloor', 0.90);
fitBucket = localBuildFloorBucket(summary, 'nonRefFitRatioFloor', 0.80);
consistencyBucket = localBuildFloorBucket(summary, 'nonRefConsistencyRatioFloor', 0.98);
coherenceBucket = localBuildFloorBucket(summary, 'nonRefCoherenceFloor', 0.98);
phaseRmsBucket = localBuildCeilBucket(summary, 'nonRefRmsPhaseResidRad', 0.01);
phaseMaxBucket = localBuildCeilBucket(summary, 'nonRefMaxAbsPhaseResidRad', 0.02);
negativeBucket = localBuildCeilBucket(summary, 'maxNonRefNegativeProjectionRatio', 0.10);
bucket = supportBucket + fitBucket + consistencyBucket + coherenceBucket + ...
  phaseRmsBucket + phaseMaxBucket + negativeBucket;
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


function value = localGetFiniteScalar(summary, fieldName, defaultValue)
%LOCALGETFINITESCALAR Read one finite scalar used directly in the score.

value = localGetFieldOrDefault(summary, fieldName, defaultValue);
if isempty(value)
  value = defaultValue;
else
  value = value(1);
end
if ~isfinite(value)
  value = defaultValue;
end
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
    value = rawValue;
  end
end
end
