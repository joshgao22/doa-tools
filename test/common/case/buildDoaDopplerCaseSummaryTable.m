function summaryTable = buildDoaDopplerCaseSummaryTable(caseList, truthList)
%BUILDDOADOPPLERCASESUMMARYTABLE Build one per-case summary table.
% truthList may be one truth struct, a struct array with one truth per
% case, or a cell array of truth structs.

arguments
  caseList (1, :) struct
  truthList
end

numCase = numel(caseList);
truthCell = localNormalizeTruthList(truthList, numCase);
caseCell = num2cell(caseList);
infoCell = cell(1, numCase);
for iCase = 1:numCase
  infoCell{iCase} = summarizeDoaDopplerCase(caseCell{iCase}, truthCell{iCase});
end
summaryTable = struct2table([infoCell{:}], 'AsArray', true);
end


function truthCell = localNormalizeTruthList(truthList, numCase)
%LOCALNORMALIZETRUTHLIST Normalize truth input to one cell per case.

if iscell(truthList)
  truthCell = reshape(truthList, 1, []);
elseif isstruct(truthList)
  if isscalar(truthList)
    truthCell = repmat({truthList}, 1, numCase);
  else
    truthCell = num2cell(reshape(truthList, 1, []));
  end
else
  error('buildDoaDopplerCaseSummaryTable:InvalidTruthList', ...
    'truthList must be a struct, struct array, or cell array.');
end

if numel(truthCell) ~= numCase
  error('buildDoaDopplerCaseSummaryTable:SizeMismatch', ...
    'truthList must have either one entry or the same number of entries as caseList.');
end
end
