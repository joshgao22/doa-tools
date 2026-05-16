function runtimeTable = buildMfRuntimeTable(repeatCell)
%BUILDMFRUNTIMETABLE Collect per-task stage timing rows.

numRepeat = numel(repeatCell);
rowCountList = zeros(numRepeat, 1);
for iRepeat = 1:numRepeat
  repeatOut = repeatCell{iRepeat};
  if isempty(repeatOut) || ~isfield(repeatOut, 'runtimeRowList') || isempty(repeatOut.runtimeRowList)
    continue;
  end
  rowCountList(iRepeat) = numel(repeatOut.runtimeRowList);
end

totalRow = sum(rowCountList);
if totalRow == 0
  runtimeTable = struct2table(repmat(emptyMfRuntimeRow(), 0, 1));
  return;
end

rowList = repmat(emptyMfRuntimeRow(), totalRow, 1);
rowOffset = 0;
for iRepeat = 1:numRepeat
  numRow = rowCountList(iRepeat);
  if numRow == 0
    continue;
  end
  rowIdx = rowOffset + (1:numRow);
  rowList(rowIdx) = repeatCell{iRepeat}.runtimeRowList(:);
  rowOffset = rowOffset + numRow;
end
runtimeTable = struct2table(rowList(:));
end
