function [scoreVec, componentNameList] = buildDynamicSelectionScoreVector(summary)
%BUILDDYNAMICSELECTIONSCOREVECTOR Build truth-free dynamic ranking scores.
% This helper is used by the dynamic flow to rank final unknown candidates
% without leaking truth-aware metrics into the selection path.
%
% scoreVec is lexicographically sorted in ascending order:
%   1) resolved-status penalty (resolved first)
%   2) final objective value (smaller is better)
%   3) final residual norm (smaller is better)
%   4) non-reference support-floor bucket (coarse unhealthy gate)
%   5) non-reference fit-floor bucket (coarse unhealthy gate)
%   6) non-reference consistency-floor bucket
%   7) non-reference coherence-floor bucket
%   8) non-reference RMS phase-residual bucket
%   9) non-reference max |phase residual|-bucket
%  10) non-reference negative-projection bucket
%  11) non-reference support-floor penalty
%  12) non-reference fit-floor penalty
%  13) non-reference consistency-floor penalty
%  14) non-reference coherence-floor penalty
%  15) non-reference RMS phase residual
%  16) non-reference max |phase residual|
%  17) non-reference negative-projection penalty
%  18) runtime in milliseconds (smaller is better)

arguments
  summary (1, 1) struct
end

componentNameList = [ ...
  "resolved status"; ...
  "final objective"; ...
  "final residual norm"; ...
  "non-reference support floor bucket"; ...
  "non-reference fit floor bucket"; ...
  "non-reference consistency floor bucket"; ...
  "non-reference coherence floor bucket"; ...
  "non-reference RMS phase residual bucket"; ...
  "non-reference max phase residual bucket"; ...
  "non-reference negative projection bucket"; ...
  "non-reference support floor"; ...
  "non-reference fit floor"; ...
  "non-reference consistency floor"; ...
  "non-reference coherence floor"; ...
  "non-reference RMS phase residual"; ...
  "non-reference max phase residual"; ...
  "non-reference negative projection"; ...
  "run time"
  ];

resolvedPenalty = double(~localGetLogicalField(summary, 'isResolved', false));
finalObj = localGetFiniteField(summary, 'finalObj', inf);
finalResidualNorm = localGetFiniteField(summary, 'finalResidualNorm', inf);
supportBucket = localBuildFloorBucket(summary, 'nonRefSupportRatioFloor', 0.90);
fitBucket = localBuildFloorBucket(summary, 'nonRefFitRatioFloor', 0.80);
consistencyBucket = localBuildFloorBucket(summary, 'nonRefConsistencyRatioFloor', 0.98);
coherenceBucket = localBuildFloorBucket(summary, 'nonRefCoherenceFloor', 0.98);
phaseRmsBucket = localBuildCeilBucket(summary, 'nonRefRmsPhaseResidRad', 0.01);
phaseMaxBucket = localBuildCeilBucket(summary, 'nonRefMaxAbsPhaseResidRad', 0.02);
negativeBucket = localBuildCeilBucket(summary, 'maxNonRefNegativeProjectionRatio', 0.10);
supportPenalty = localGetQualityPenalty(summary, 'nonRefSupportRatioFloor');
fitPenalty = localGetQualityPenalty(summary, 'nonRefFitRatioFloor');
consistencyPenalty = localGetQualityPenalty(summary, 'nonRefConsistencyRatioFloor');
coherencePenalty = localGetQualityPenalty(summary, 'nonRefCoherenceFloor');
nonRefRmsPhaseResidRad = localGetFiniteField(summary, 'nonRefRmsPhaseResidRad', inf);
nonRefMaxAbsPhaseResidRad = localGetFiniteField(summary, 'nonRefMaxAbsPhaseResidRad', inf);
negativePenalty = localGetFiniteField(summary, 'maxNonRefNegativeProjectionRatio', 0);
negativePenalty = min(max(negativePenalty, 0), 1);
runTimeMs = localGetFiniteField(summary, 'runTimeMs', inf);
scoreVec = [resolvedPenalty, finalObj, finalResidualNorm, supportBucket, ...
  fitBucket, consistencyBucket, coherenceBucket, phaseRmsBucket, ...
  phaseMaxBucket, negativeBucket, supportPenalty, fitPenalty, ...
  consistencyPenalty, coherencePenalty, nonRefRmsPhaseResidRad, ...
  nonRefMaxAbsPhaseResidRad, negativePenalty, runTimeMs];
end

function value = localBuildFloorBucket(dataStruct, fieldName, minValue)
%LOCALBUILDFLOORBUCKET Build one coarse penalty bucket for lower-bounded metrics.

metricValue = localGetFiniteField(dataStruct, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
metricValue = min(max(metricValue, 0), 1);
value = double(metricValue < minValue);
end


function value = localBuildCeilBucket(dataStruct, fieldName, maxValue)
%LOCALBUILDCEILBUCKET Build one coarse penalty bucket for upper-bounded metrics.

metricValue = localGetFiniteField(dataStruct, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
value = double(metricValue > maxValue);
end


function value = localGetQualityPenalty(dataStruct, fieldName)
%LOCALGETQUALITYPENALTY Convert one unit-interval quality metric into a score penalty.

metricValue = localGetFiniteField(dataStruct, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
metricValue = min(max(metricValue, 0), 1);
value = 1 - metricValue;
end


function value = localGetFiniteField(dataStruct, fieldName, defaultValue)
%LOCALGETFINITEFIELD Read one scalar field and sanitize non-finite values.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end

rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
rawValue = rawValue(1);
if isfinite(rawValue)
  value = rawValue;
end
end


function value = localGetLogicalField(dataStruct, fieldName, defaultValue)
%LOCALGETLOGICALFIELD Read one logical-like scalar field.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end

rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
value = logical(rawValue(1));
end
