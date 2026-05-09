function topTable = buildMfTopSlowRuntimeTable(runtimeTable, maxRow)
%BUILDMFTOPSLOWRUNTIMETABLE Select the slowest replay/scan runtime rows.
% This generic reporting helper only filters already-collected runtime rows.

if nargin < 2 || isempty(maxRow)
  maxRow = 12;
end
if isempty(runtimeTable) || height(runtimeTable) == 0
  topTable = runtimeTable;
  return;
end

runtimeTable = sortrows(runtimeTable, 'wallTimeMs', 'descend');
maxRow = min(maxRow, height(runtimeTable));
topTable = runtimeTable(1:maxRow, :);
end
