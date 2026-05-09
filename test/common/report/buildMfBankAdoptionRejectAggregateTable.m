function rejectAggregateTable = buildMfBankAdoptionRejectAggregateTable(shadowTable)
%BUILDMFBANKADOPTIONREJECTAGGREGATETABLE Count shadow adoption outcomes by reason.
%
% This helper only aggregates an existing shadow-decision table. It does not
% evaluate objectives, run solver calls, or change estimator output.

emptyRows = repmat(localEmptyRejectAggregateRow(), 0, 1);
if isempty(shadowTable) || ~istable(shadowTable) || height(shadowTable) == 0 || ...
    ~all(ismember({'displayName', 'snrDb', 'rejectReason'}, shadowTable.Properties.VariableNames))
  rejectAggregateTable = struct2table(emptyRows(:));
  return;
end
keyTable = unique(shadowTable(:, {'displayName', 'snrDb', 'rejectReason'}), 'rows', 'stable');
rowList = repmat(localEmptyRejectAggregateRow(), 0, 1);
for iKey = 1:height(keyTable)
  mask = shadowTable.displayName == keyTable.displayName(iKey) & ...
    shadowTable.snrDb == keyTable.snrDb(iKey) & ...
    shadowTable.rejectReason == keyTable.rejectReason(iKey);
  rows = shadowTable(mask, :);
  row = localEmptyRejectAggregateRow();
  row.displayName = keyTable.displayName(iKey);
  row.snrDb = keyTable.snrDb(iKey);
  row.rejectReason = keyTable.rejectReason(iKey);
  row.numCase = height(rows);
  row.numCandidateMedian = localMedianIfPresent(rows, 'numCandidate');
  row.numEligibleMedian = localMedianIfPresent(rows, 'numEligible');
  row.baselineAngleMedianDeg = localMedianIfPresent(rows, 'baselineAngleErrDeg');
  if ismember('shadowAngleImproveDeg', rows.Properties.VariableNames) && ...
      ismember('shadowAdoptFlag', rows.Properties.VariableNames)
    adoptedRows = rows(rows.shadowAdoptFlag, :);
    row.adoptedAngleImproveMedianDeg = localMedianIfPresent(adoptedRows, 'shadowAngleImproveDeg');
  end
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
rejectAggregateTable = sortrows(struct2table(rowList(:)), {'displayName', 'snrDb', 'rejectReason'});
end

function value = localMedianIfPresent(rowTable, varName)
value = NaN;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames) && height(rowTable) > 0
  value = median(rowTable.(varName), 'omitnan');
end
end

function row = localEmptyRejectAggregateRow()
row = struct('displayName', "", 'snrDb', NaN, 'rejectReason', "", ...
  'numCase', NaN, 'numCandidateMedian', NaN, 'numEligibleMedian', NaN, ...
  'baselineAngleMedianDeg', NaN, 'adoptedAngleImproveMedianDeg', NaN);
end
