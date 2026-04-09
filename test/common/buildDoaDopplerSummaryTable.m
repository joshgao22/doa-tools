function estTable = buildDoaDopplerSummaryTable(caseResult, truth, summaryOpt)
%BUILDDOADOPPLERSUMMARYTABLE Build a displayed summary table for dev cases.
%
%summaryOpt.mode:
%  'static'  - keep the compact single-frame static columns
%  'dynamic' - keep the dynamic multi-frame columns
%  'full'    - keep all available columns

if nargin < 3 || isempty(summaryOpt)
  summaryOpt = struct();
end
modeName = string(getDoaDopplerFieldOrDefault(summaryOpt, 'mode', "full"));

numCase = numel(caseResult);
infoCell = cell(1, numCase);
for iCase = 1:numCase
  infoCell{iCase} = summarizeDoaDopplerCase(caseResult(iCase), truth);
end
estTable = struct2table([infoCell{:}], 'AsArray', true);

if modeName == "static"
  keepName = {'displayName', 'satMode', 'paramMode', 'dynamicMode', ...
    'latEstDeg', 'lonEstDeg', 'latErrDeg', 'lonErrDeg', 'angleErrDeg', ...
    'fdRefEstHz', 'fdRefErrHz', 'fdRateEstHzPerSec', 'fdSatRmseHz', ...
    'deltaFdRmseHz', 'exitflag', 'isResolved', 'runTimeMs', ...
    'funcCount', 'iterations'};
  estTable = localKeepColumns(estTable, keepName);
elseif modeName == "dynamic"
  keepName = {'displayName', 'satMode', 'frameMode', 'paramMode', ...
    'dynamicMode', 'latEstDeg', 'lonEstDeg', 'latErrDeg', 'lonErrDeg', ...
    'angleErrDeg', 'fdRefEstHz', 'fdRefErrHz', 'fdRateEstHzPerSec', ...
    'fdRateErrHzPerSec', 'fdLineRmseHz', 'fdSatRmseHz', 'exitflag', ...
    'isResolved', 'runTimeMs', 'funcCount', 'iterations'};
  estTable = localKeepColumns(estTable, keepName);
end
end

function outTable = localKeepColumns(inTable, keepName)
%LOCALKEEPCOLUMNS Keep one ordered subset of columns when they exist.

hasName = ismember(keepName, inTable.Properties.VariableNames);
outTable = inTable(:, keepName(hasName));
end
