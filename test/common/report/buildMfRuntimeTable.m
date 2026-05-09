function runtimeTable = buildMfRuntimeTable(repeatCell)
%BUILDMFRUNTIMETABLE Collect per-task stage timing rows.

rowList = repmat(emptyMfRuntimeRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isfield(repeatOut, 'runtimeRowList') && ~isempty(repeatOut.runtimeRowList)
    rowList = [rowList; repeatOut.runtimeRowList(:)]; %#ok<AGROW>
  end
end
runtimeTable = struct2table(rowList(:));
end
