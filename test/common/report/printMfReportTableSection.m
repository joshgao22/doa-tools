function printMfReportTableSection(sectionTitle, dataTable, opt)
%PRINTMFREPORTTABLESECTION Print compact replay/scan report tables.
% This helper is shared by replay and scan scripts. It only controls command
% line verbosity; full tables should remain in replayData/scanData.

arguments
  sectionTitle (1,:) char
  dataTable
  opt.Level (1,1) string = "compact"
  opt.RequiredLevel (1,1) string = "compact"
  opt.MaxPreviewRow (1,1) double {mustBeInteger, mustBeNonnegative} = 8
  opt.ColumnNames = string.empty(1, 0)
  opt.SectionPrinter = []
end

if ~localShouldPrint(opt.Level, opt.RequiredLevel)
  return;
end

printer = opt.SectionPrinter;
if isempty(printer)
  printer = @localDefaultSectionPrinter;
end

tableUse = dataTable;
if istable(tableUse) && ~isempty(opt.ColumnNames)
  keepCols = intersect(string(opt.ColumnNames(:)).', string(tableUse.Properties.VariableNames), 'stable');
  if ~isempty(keepCols)
    tableUse = tableUse(:, cellstr(keepCols));
  end
end

if istable(tableUse) && height(tableUse) > opt.MaxPreviewRow
  summaryTable = table(height(tableUse), 'VariableNames', {'numRows'});
  printer(sectionTitle, summaryTable);
  fprintf('%s preview:\n', sectionTitle);
  localDispTablePreview(tableUse, opt.MaxPreviewRow);
else
  printer(sectionTitle, tableUse);
end
end

function tf = localShouldPrint(level, requiredLevel)
%LOCALSHOULDPRINT Return whether one section belongs to the current level.
order = ["compact", "diagnostic", "full"];
levelIdx = find(order == string(level), 1, 'first');
requiredIdx = find(order == string(requiredLevel), 1, 'first');
if isempty(levelIdx)
  levelIdx = 1;
end
if isempty(requiredIdx)
  requiredIdx = 1;
end
tf = levelIdx >= requiredIdx;
end

function localDefaultSectionPrinter(sectionTitle, dataTable)
%LOCALDEFAULTSECTIONPRINTER Minimal section printer fallback.
fprintf('\n%s\n', sectionTitle);
disp(dataTable);
end

function localDispTablePreview(dataTable, edgeCount)
%LOCALDISPTABLEPREVIEW Print first/last rows without replay-only dependency.
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
