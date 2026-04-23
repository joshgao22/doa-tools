function repTable = selectDynamicRepresentativeRepeats(repeatTable, representativeLabelList)
%SELECTDYNAMICREPRESENTATIVEREPEATS Select one configurable set of dynamic repeats.

arguments
  repeatTable table
  representativeLabelList (:, 1) string = ["bestMsUnknownAngleGain"; "worstMsUnknownAngleGain"; "medianMsUnknown"]
end

labelList = strings(0, 1);
idxList = zeros(0, 1);
representativeLabelList = reshape(string(representativeLabelList), [], 1);

if any(representativeLabelList == "bestMsUnknownAngleGain")
  [bestGainIdx, hasBestGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'min');
  if hasBestGain
    labelList(end+1, 1) = "bestMsUnknownAngleGain";
    idxList(end+1, 1) = repeatTable.repeatIdx(bestGainIdx);
  end
end

if any(representativeLabelList == "worstMsUnknownAngleGain")
  [worstGainIdx, hasWorstGain] = localSelectByMetric(repeatTable.msMinusSsAngleDeg, 'max');
  if hasWorstGain
    labelList(end+1, 1) = "worstMsUnknownAngleGain";
    idxList(end+1, 1) = repeatTable.repeatIdx(worstGainIdx);
  end
end

if any(representativeLabelList == "medianMsUnknown")
  [medianMsIdx, hasMedianMs] = localSelectMedianResolved(repeatTable.msAngleErrDeg, repeatTable.msResolved);
  if hasMedianMs
    labelList(end+1, 1) = "medianMsUnknown";
    idxList(end+1, 1) = repeatTable.repeatIdx(medianMsIdx);
  end
end

[idxUnique, keepIdx] = unique(idxList, 'stable');
taskSeedList = nan(numel(idxUnique), 1);
for iRow = 1:numel(idxUnique)
  matchIdx = find(repeatTable.repeatIdx == idxUnique(iRow), 1, 'first');
  if ~isempty(matchIdx)
    taskSeedList(iRow) = repeatTable.taskSeed(matchIdx);
  end
end
repTable = table(labelList(keepIdx), idxUnique, taskSeedList, ...
  'VariableNames', {'label', 'repeatIdx', 'taskSeed'});
end


function [selIdx, hasValue] = localSelectByMetric(metricVec, modeName)
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
    error('selectDynamicRepresentativeRepeats:InvalidSelectMode', ...
      'Unsupported select mode: %s.', modeName);
end
finiteIdx = find(finiteMask);
selIdx = finiteIdx(idxLocal);
end

function [selIdx, hasValue] = localSelectMedianResolved(metricVec, resolvedMask)
metricVec = reshape(metricVec, [], 1);
resolvedMask = reshape(resolvedMask, [], 1) & isfinite(metricVec);
selIdx = NaN;
hasValue = any(resolvedMask);
if ~hasValue
  return;
end
metricResolved = metricVec(resolvedMask);
medianVal = median(metricResolved, 'omitnan');
resolvedIdx = find(resolvedMask);
[~, idxLocal] = min(abs(metricResolved - medianVal));
selIdx = resolvedIdx(idxLocal);
end
