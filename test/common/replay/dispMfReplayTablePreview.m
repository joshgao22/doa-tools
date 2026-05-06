function dispMfReplayTablePreview(dataTable, edgeCount)
%DISPMFREPLAYTABLEPREVIEW Print compact first/last rows of long replay tables.

arguments
  dataTable
  edgeCount (1,1) double {mustBeInteger, mustBeNonnegative} = 8
end

if isempty(dataTable) || ~istable(dataTable)
  disp(dataTable);
  return;
end
if height(dataTable) <= 2 * edgeCount
  disp(dataTable);
  return;
end
fprintf('  showing first %d and last %d of %d rows\n', edgeCount, edgeCount, height(dataTable));
disp(dataTable(1:edgeCount, :));
fprintf('  ...\n');
disp(dataTable((height(dataTable) - edgeCount + 1):height(dataTable), :));
end
