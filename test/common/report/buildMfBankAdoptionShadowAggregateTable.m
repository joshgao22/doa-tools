function aggregateTable = buildMfBankAdoptionShadowAggregateTable(shadowTable, varargin)
%BUILDMFBANKADOPTIONSHADOWAGGREGATETABLE Aggregate bank adoption shadow outcomes.
%
% This helper summarizes an offline shadow-selection table. It does not define
% a formal estimator winner rule and does not change any runtime estimate.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyAggregateRow(), 0, 1);
if isempty(shadowTable) || ~istable(shadowTable) || height(shadowTable) == 0
  aggregateTable = struct2table(rowList(:));
  return;
end
requiredVars = {'displayName', 'snrDb', 'shadowAdoptFlag', 'shadowAngleImproveFlag', ...
  'shadowDamageFlag', 'shadowObjMinusBaseline', 'shadowAngleImproveDeg', ...
  'shadowAngleImproveOverCrb'};
if ~all(ismember(requiredVars, shadowTable.Properties.VariableNames))
  aggregateTable = struct2table(rowList(:));
  return;
end
if ~ismember('primaryDiagnosis', shadowTable.Properties.VariableNames)
  shadowTable.primaryDiagnosis = repmat("all", height(shadowTable), 1);
end
groupBy = opt.GroupBy;
groupBy = groupBy(ismember(cellstr(groupBy), shadowTable.Properties.VariableNames));
if isempty(groupBy)
  groupBy = ["displayName", "snrDb"];
end
keyTable = unique(shadowTable(:, cellstr(groupBy)), 'rows', 'stable');
for iKey = 1:height(keyTable)
  mask = true(height(shadowTable), 1);
  for iVar = 1:numel(groupBy)
    varName = char(groupBy(iVar));
    mask = mask & shadowTable.(varName) == keyTable.(varName)(iKey);
  end
  rows = shadowTable(mask, :);
  row = localEmptyAggregateRow();
  row.displayName = localKeyValue(keyTable, iKey, 'displayName', "all");
  row.snrDb = localKeyValue(keyTable, iKey, 'snrDb', NaN);
  row.primaryDiagnosis = localKeyValue(keyTable, iKey, 'primaryDiagnosis', "all");
  row.numCase = height(rows);
  row.adoptCount = nnz(rows.shadowAdoptFlag);
  row.adoptRate = localSafeRatio(row.adoptCount, row.numCase);
  row.angleImproveCount = nnz(rows.shadowAdoptFlag & rows.shadowAngleImproveFlag);
  row.angleImproveRate = localSafeRatio(row.angleImproveCount, row.numCase);
  row.damageCount = nnz(rows.shadowAdoptFlag & rows.shadowDamageFlag);
  row.damageRate = localSafeRatio(row.damageCount, row.numCase);
  adoptRows = rows(rows.shadowAdoptFlag, :);
  row.medianObjDelta = median(adoptRows.shadowObjMinusBaseline, 'omitnan');
  row.medianAngleImproveDeg = median(adoptRows.shadowAngleImproveDeg, 'omitnan');
  row.medianAngleImproveOverCrb = median(adoptRows.shadowAngleImproveOverCrb, 'omitnan');
  row.medianFdRefErrOverTooth = localMedianIfPresent(adoptRows, 'shadowFdRefErrOverTooth');
  row.medianFdToothShift = localMedianIfPresent(adoptRows, 'shadowFdToothShift');
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
aggregateTable = sortrows(struct2table(rowList(:)), {'displayName', 'snrDb', 'primaryDiagnosis'});
end

function opt = localParseOpt(varargin)
opt = struct('GroupBy', ["displayName", "snrDb"]);
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankAdoptionShadowAggregateTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = lower(string(varargin{iArg}));
  value = varargin{iArg + 1};
  switch name
    case "groupby"
      opt.GroupBy = reshape(string(value), 1, []);
    otherwise
      error('buildMfBankAdoptionShadowAggregateTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function value = localKeyValue(keyTable, iKey, varName, defaultValue)
value = defaultValue;
if istable(keyTable) && ismember(varName, keyTable.Properties.VariableNames)
  value = keyTable.(varName)(iKey);
end
end

function value = localMedianIfPresent(rowTable, varName)
value = NaN;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames) && height(rowTable) > 0
  value = median(rowTable.(varName), 'omitnan');
end
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function row = localEmptyAggregateRow()
row = struct('displayName', "", 'snrDb', NaN, 'primaryDiagnosis', "", 'numCase', NaN, ...
  'adoptCount', NaN, 'adoptRate', NaN, ...
  'angleImproveCount', NaN, 'angleImproveRate', NaN, ...
  'damageCount', NaN, 'damageRate', NaN, ...
  'medianObjDelta', NaN, 'medianAngleImproveDeg', NaN, ...
  'medianAngleImproveOverCrb', NaN, 'medianFdRefErrOverTooth', NaN, ...
  'medianFdToothShift', NaN);
end
