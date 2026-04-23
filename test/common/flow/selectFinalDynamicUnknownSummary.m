function [selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary(summaryCell, candidateTag, opt)
%SELECTFINALDYNAMICUNKNOWNSUMMARY Select one final unknown candidate.
% candidateTag and summaryCell must have the same number of elements.

arguments
  summaryCell cell
  candidateTag (:, 1) string
  opt (1, 1) struct = struct()
end

numCase = numel(summaryCell);
if numCase <= 0
  error('selectFinalDynamicUnknownSummary:EmptySummaryCell', ...
    'summaryCell must contain at least one candidate summary.');
end
if numCase ~= numel(candidateTag)
  error('selectFinalDynamicUnknownSummary:SizeMismatch', ...
    'summaryCell and candidateTag must have the same length.');
end

opt = localApplySelectionDefaults(opt);
[firstScoreVec, componentNameList] = buildDynamicSelectionScoreVector(summaryCell{1});
numComponent = numel(firstScoreVec);
scoreMat = inf(numCase, numComponent);
scoreMat(1, :) = firstScoreVec;
for iCase = 2:numCase
  [scoreMat(iCase, :), componentNameList] = buildDynamicSelectionScoreVector(summaryCell{iCase});
end

[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
selectedIdx = orderIdx(1);
selectedTag = candidateTag(selectedIdx);
reasonText = localBuildReasonText(orderIdx, scoreMat, componentNameList, candidateTag);
[selectedIdx, selectedTag, reasonText] = localApplySubsetAnchorGuard( ...
  summaryCell, candidateTag, scoreMat, selectedIdx, selectedTag, reasonText, opt);
[selectedIdx, selectedTag, reasonText] = localApplySameToothReleasePreference( ...
  summaryCell, candidateTag, scoreMat, selectedIdx, selectedTag, reasonText, opt);
end


function opt = localApplySelectionDefaults(opt)
%LOCALAPPLYSELECTIONDEFAULTS Fill one final-selection option struct.

if nargin < 1 || isempty(opt)
  opt = struct();
end
if ~isfield(opt, 'enableSubsetAnchorGuard') || isempty(opt.enableSubsetAnchorGuard)
  opt.enableSubsetAnchorGuard = false;
end
if ~isfield(opt, 'subsetAnchorIdx') || isempty(opt.subsetAnchorIdx)
  opt.subsetAnchorIdx = NaN;
end
if ~isfield(opt, 'subsetAnchorTag') || isempty(opt.subsetAnchorTag)
  opt.subsetAnchorTag = "subset-anchor";
end
if ~isfield(opt, 'maxFdRefDriftHz') || isempty(opt.maxFdRefDriftHz)
  opt.maxFdRefDriftHz = inf;
end
if ~isfield(opt, 'maxFdRateDriftHzPerSec') || isempty(opt.maxFdRateDriftHzPerSec)
  opt.maxFdRateDriftHzPerSec = inf;
end
if ~isfield(opt, 'maxDoaDriftDeg') || isempty(opt.maxDoaDriftDeg)
  opt.maxDoaDriftDeg = inf;
end
if ~isfield(opt, 'minObjGainToLeaveAnchor') || isempty(opt.minObjGainToLeaveAnchor)
  opt.minObjGainToLeaveAnchor = 0;
end
if ~isfield(opt, 'minResidualGainToLeaveAnchor') || isempty(opt.minResidualGainToLeaveAnchor)
  opt.minResidualGainToLeaveAnchor = 0;
end
if ~isfield(opt, 'sameToothBypassGuard') || isempty(opt.sameToothBypassGuard)
  opt.sameToothBypassGuard = true;
end
if ~isfield(opt, 'closerToZeroToothBypassGuard') || isempty(opt.closerToZeroToothBypassGuard)
  opt.closerToZeroToothBypassGuard = true;
end
if ~isfield(opt, 'preferReleasedSameTooth') || isempty(opt.preferReleasedSameTooth)
  opt.preferReleasedSameTooth = false;
end
if ~isfield(opt, 'sameToothObjectiveToleranceAbs') || isempty(opt.sameToothObjectiveToleranceAbs)
  opt.sameToothObjectiveToleranceAbs = 0;
end
if ~isfield(opt, 'sameToothResidualToleranceAbs') || isempty(opt.sameToothResidualToleranceAbs)
  opt.sameToothResidualToleranceAbs = 0;
end
if ~isfield(opt, 'sameToothObjectiveToleranceRel') || isempty(opt.sameToothObjectiveToleranceRel)
  opt.sameToothObjectiveToleranceRel = 0;
end
if ~isfield(opt, 'sameToothResidualToleranceRel') || isempty(opt.sameToothResidualToleranceRel)
  opt.sameToothResidualToleranceRel = 0;
end
end


function reasonText = localBuildReasonText(orderIdx, scoreMat, componentNameList, candidateTag)
%LOCALBUILDREASONTEXT Build one compact objective-only selection reason.

selectedTag = candidateTag(orderIdx(1));
reasonText = "tie";
if numel(orderIdx) < 2
  return;
end
bestScore = scoreMat(orderIdx(1), :);
secondScore = scoreMat(orderIdx(2), :);
diffIdx = find(bestScore ~= secondScore, 1, 'first');
if ~isempty(diffIdx)
  reasonText = selectedTag + " selected by smaller " + componentNameList(diffIdx);
end
end


function [selectedIdx, selectedTag, reasonText] = localApplySubsetAnchorGuard( ...
  summaryCell, candidateTag, scoreMat, selectedIdx, selectedTag, reasonText, opt)
%LOCALAPPLYSUBSETANCHORGUARD Keep one trusted subset-anchor when needed.

if ~logical(opt.enableSubsetAnchorGuard)
  return;
end
anchorIdx = opt.subsetAnchorIdx;
if ~isscalar(anchorIdx) || ~isfinite(anchorIdx) || anchorIdx < 1 || anchorIdx > numel(summaryCell)
  return;
end
if selectedIdx == anchorIdx
  return;
end
anchorSummary = summaryCell{anchorIdx};
selectedSummary = summaryCell{selectedIdx};
if ~logical(localGetLogicalField(anchorSummary, 'isResolved', false))
  return;
end
if ~logical(localGetLogicalField(selectedSummary, 'isResolved', false))
  return;
end

fdRefDriftHz = abs(localGetFiniteField(selectedSummary, 'fdRefEst', NaN) - ...
  localGetFiniteField(anchorSummary, 'fdRefEst', NaN));
fdRateDriftHzPerSec = abs(localGetFiniteField(selectedSummary, 'fdRateEst', NaN) - ...
  localGetFiniteField(anchorSummary, 'fdRateEst', NaN));
doaDriftDeg = localCalcDoaDriftDeg(selectedSummary, anchorSummary);
anchorObj = localGetFiniteField(anchorSummary, 'finalObj', inf);
selectedObj = localGetFiniteField(selectedSummary, 'finalObj', inf);
anchorResidual = localGetFiniteField(anchorSummary, 'finalResidualNorm', inf);
selectedResidual = localGetFiniteField(selectedSummary, 'finalResidualNorm', inf);
objGain = anchorObj - selectedObj;
residualGain = anchorResidual - selectedResidual;
anchorToothIdx = localGetFiniteField(anchorSummary, 'toothIdx', NaN);
selectedToothIdx = localGetFiniteField(selectedSummary, 'toothIdx', NaN);

if logical(opt.sameToothBypassGuard) && isfinite(anchorToothIdx) && isfinite(selectedToothIdx) && ...
    selectedToothIdx == anchorToothIdx
  return;
end
if logical(opt.closerToZeroToothBypassGuard) && isfinite(anchorToothIdx) && isfinite(selectedToothIdx) && ...
    abs(selectedToothIdx) < abs(anchorToothIdx)
  return;
end

hasLargeDrift = (fdRefDriftHz > opt.maxFdRefDriftHz) || ...
  (fdRateDriftHzPerSec > opt.maxFdRateDriftHzPerSec) || ...
  (doaDriftDeg > opt.maxDoaDriftDeg);
hasSmallScoreGain = (objGain < opt.minObjGainToLeaveAnchor) && ...
  (residualGain < opt.minResidualGainToLeaveAnchor);
if ~(hasLargeDrift && hasSmallScoreGain)
  return;
end

selectedIdx = anchorIdx;
selectedTag = candidateTag(anchorIdx);
reasonText = string(localGetFieldOrDefault(opt, 'subsetAnchorTag', "subset-anchor")) + ...
  " kept by subset-anchor guard";
end


function [selectedIdx, selectedTag, reasonText] = localApplySameToothReleasePreference( ...
  summaryCell, candidateTag, scoreMat, selectedIdx, selectedTag, reasonText, opt)
%LOCALAPPLYSAMETOOTHRELEASEPREFERENCE Prefer released same-tooth candidates.
% When the best candidates already agree on the same tooth and objective /
% residual gaps are numerically tiny, prefer the candidate that is not a
% fixed-DoA anchor and that keeps better non-reference phase quality.

if ~logical(localGetFieldOrDefault(opt, 'preferReleasedSameTooth', false))
  return;
end
bestSummary = summaryCell{selectedIdx};
if ~logical(localGetLogicalField(bestSummary, 'isResolved', false))
  return;
end
bestToothIdx = localGetFiniteField(bestSummary, 'toothIdx', NaN);
if ~isfinite(bestToothIdx)
  return;
end
bestObj = localGetFiniteField(bestSummary, 'finalObj', inf);
bestResidual = localGetFiniteField(bestSummary, 'finalResidualNorm', inf);
objTol = localBuildSameToothTolerance(bestObj, ...
  localGetFiniteField(opt, 'sameToothObjectiveToleranceAbs', 0), ...
  localGetFiniteField(opt, 'sameToothObjectiveToleranceRel', 0));
residualTol = localBuildSameToothTolerance(bestResidual, ...
  localGetFiniteField(opt, 'sameToothResidualToleranceAbs', 0), ...
  localGetFiniteField(opt, 'sameToothResidualToleranceRel', 0));

candidateIdx = zeros(0, 1);
releaseScoreMat = inf(0, 9);
for iCase = 1:numel(summaryCell)
  summaryUse = summaryCell{iCase};
  if ~logical(localGetLogicalField(summaryUse, 'isResolved', false))
    continue;
  end
  toothIdx = localGetFiniteField(summaryUse, 'toothIdx', NaN);
  if ~(isfinite(toothIdx) && toothIdx == bestToothIdx)
    continue;
  end
  finalObj = localGetFiniteField(summaryUse, 'finalObj', inf);
  finalResidual = localGetFiniteField(summaryUse, 'finalResidualNorm', inf);
  if finalObj > bestObj + objTol || finalResidual > bestResidual + residualTol
    continue;
  end
  candidateIdx(end + 1, 1) = iCase; %#ok<AGROW>
  releaseScoreMat(end + 1, :) = localBuildSameToothReleaseScore(summaryUse, candidateTag(iCase)); %#ok<AGROW>
end

if numel(candidateIdx) <= 1
  return;
end
[~, orderIdx] = sortrows(releaseScoreMat, 1:size(releaseScoreMat, 2));
selectedReleaseIdx = candidateIdx(orderIdx(1));
if selectedReleaseIdx == selectedIdx
  return;
end
selectedIdx = selectedReleaseIdx;
selectedTag = candidateTag(selectedIdx);
reasonText = selectedTag + " kept by same-tooth release preference";
end


function scoreVec = localBuildSameToothReleaseScore(summaryUse, candidateTag)
%LOCALBUILDSAMETOOTHRELEASESCORE Build one same-tooth tie-break score.

isFrozenLike = localGetLogicalField(summaryUse, 'isDoaFrozenLike', false);
if contains(lower(char(candidateTag)), 'anchor') && ~contains(lower(char(candidateTag)), 'polish')
  isFrozenLike = true;
end
supportPenalty = localQualityPenalty(summaryUse, 'nonRefSupportRatioFloor');
fitPenalty = localQualityPenalty(summaryUse, 'nonRefFitRatioFloor');
consistencyPenalty = localQualityPenalty(summaryUse, 'nonRefConsistencyRatioFloor');
coherencePenalty = localQualityPenalty(summaryUse, 'nonRefCoherenceFloor');
phaseRms = localGetFiniteField(summaryUse, 'nonRefRmsPhaseResidRad', inf);
phaseMax = localGetFiniteField(summaryUse, 'nonRefMaxAbsPhaseResidRad', inf);
finalObj = localGetFiniteField(summaryUse, 'finalObj', inf);
finalResidual = localGetFiniteField(summaryUse, 'finalResidualNorm', inf);
runTimeMs = localGetFiniteField(summaryUse, 'runTimeMs', inf);
scoreVec = [double(isFrozenLike), supportPenalty, fitPenalty, consistencyPenalty, ...
  coherencePenalty, phaseRms, phaseMax, finalObj, finalResidual + 1e-6 * runTimeMs];
end


function tol = localBuildSameToothTolerance(baseValue, absTol, relTol)
%LOCALBUILDSAMETOOTHTOLERANCE Build one mixed absolute/relative tolerance.

tol = absTol;
if ~isfinite(tol)
  tol = 0;
end
if isfinite(baseValue) && isfinite(relTol) && relTol > 0
  tol = max(tol, relTol * max(abs(baseValue), 1));
end
end


function value = localQualityPenalty(dataStruct, fieldName)
%LOCALQUALITYPENALTY Convert one unit-interval quality metric to a penalty.

metricValue = localGetFiniteField(dataStruct, fieldName, NaN);
if ~isfinite(metricValue)
  value = 0;
  return;
end
metricValue = min(max(metricValue, 0), 1);
value = 1 - metricValue;
end


function doaDriftDeg = localCalcDoaDriftDeg(summaryA, summaryB)
%LOCALCALCDOADRIFTDEG Compute one symmetric DoA drift metric.

doaDriftDeg = inf;
doaA = localGetFieldOrDefault(summaryA, 'doaParamEst', []);
doaB = localGetFieldOrDefault(summaryB, 'doaParamEst', []);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = inf;
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


function value = localGetFiniteField(dataStruct, fieldName, defaultValue)
%LOCALGETFINITEFIELD Read one scalar field with finite-value fallback.

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
