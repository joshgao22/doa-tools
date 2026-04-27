function routeSummary = buildDynamicRouteSummaryTable(snrDbList, selectedSubsetLabel, selectedFinalTag)
%BUILDDYNAMICROUTESUMMARYTABLE Summarize selected subset/final route counts by SNR.
% Inputs are numRepeat-by-numSnr string-like arrays from dynamic perf output.

snrDbList = snrDbList(:).';
selectedSubsetLabel = string(selectedSubsetLabel);
selectedFinalTag = string(selectedFinalTag);

routeSummary = table();
for iSnr = 1:numel(snrDbList)
  routeSummary = [routeSummary; localBuildOneRouteTable(snrDbList(iSnr), "subset", selectedSubsetLabel(:, iSnr))]; %#ok<AGROW>
  routeSummary = [routeSummary; localBuildOneRouteTable(snrDbList(iSnr), "final", selectedFinalTag(:, iSnr))]; %#ok<AGROW>
end
end

function oneTable = localBuildOneRouteTable(snrDb, routeType, routeLabelVec)
routeLabelVec = routeLabelVec(:);
routeLabelVec(routeLabelVec == "") = "none";
[labelList, ~, labelIdx] = unique(routeLabelVec, 'stable');
countVec = accumarray(labelIdx, 1, [numel(labelList), 1]);
fractionVec = countVec ./ max(1, numel(routeLabelVec));
oneTable = table(repmat(snrDb, numel(labelList), 1), repmat(routeType, numel(labelList), 1), ...
  labelList, countVec, fractionVec, ...
  'VariableNames', {'snrDb', 'routeType', 'routeLabel', 'count', 'fraction'});
end
