function runtimeAggregateTable = buildMfRuntimeAggregateTable(runtimeTable)
%BUILDMFRUNTIMEAGGREGATETABLE Aggregate replay/scan runtime rows by stage.
% The input table must contain methodName, stageName, stageGroup, and
% wallTimeMs columns. The helper is experiment-side reporting only and does
% not affect estimator numerical paths.

rowTemplate = localEmptyAggregateRow();
if isempty(runtimeTable) || height(runtimeTable) == 0
  runtimeAggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end

keyTable = unique(runtimeTable(:, {'methodName', 'stageName', 'stageGroup'}), 'rows');
rowList = repmat(rowTemplate, height(keyTable), 1);
for iRow = 1:height(keyTable)
  mask = runtimeTable.methodName == keyTable.methodName(iRow) & ...
    runtimeTable.stageName == keyTable.stageName(iRow) & ...
    runtimeTable.stageGroup == keyTable.stageGroup(iRow);
  valueSec = runtimeTable.wallTimeMs(mask) / 1000;
  valueSec = valueSec(isfinite(valueSec));
  row = rowTemplate;
  row.methodName = keyTable.methodName(iRow);
  row.stageName = keyTable.stageName(iRow);
  row.stageGroup = keyTable.stageGroup(iRow);
  row.numCase = numel(valueSec);
  row.totalSec = sum(valueSec);
  row.meanSec = mean(valueSec, 'omitnan');
  row.medianSec = median(valueSec, 'omitnan');
  row.maxSec = max(valueSec, [], 'omitnan');
  rowList(iRow) = row;
end
runtimeAggregateTable = sortrows(struct2table(rowList), 'totalSec', 'descend');
end

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return one typed runtime aggregate row.
row = struct('methodName', "", 'stageName', "", 'stageGroup', "", ...
  'numCase', NaN, 'totalSec', NaN, 'meanSec', NaN, 'medianSec', NaN, ...
  'maxSec', NaN);
end
