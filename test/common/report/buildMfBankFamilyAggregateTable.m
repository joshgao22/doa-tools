function aggregateTable = buildMfBankFamilyAggregateTable(candidateTraceTable, varargin)
%BUILDMFBANKFAMILYAGGREGATETABLE Aggregate replay-only MS bank candidates by family.
%
% The helper summarizes existing candidate trace rows by candidate source,
% center source, axis direction, step size, and optional primary diagnosis.
% It does not run solvers, choose winners, or define an estimator adoption rule.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyAggregateRow(), 0, 1);
if isempty(candidateTraceTable) || ~istable(candidateTraceTable) || height(candidateTraceTable) == 0
  aggregateTable = struct2table(rowList(:));
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed', 'candidateSource', ...
  'bankCenterSource', 'bankAxisDirection', 'bankStepDeg', 'wouldBeatBaselineFlag', ...
  'angleImproveDeg', 'angleImproveOverCrb', 'damageFlag', 'objMinusBaseline', ...
  'finalAngleErrDeg'};
if ~all(ismember(requiredVars, candidateTraceTable.Properties.VariableNames))
  aggregateTable = struct2table(rowList(:));
  return;
end
rowsAll = candidateTraceTable(ismember(candidateTraceTable.candidateSource, opt.CandidateSourceList), :);
rowsAll = rowsAll(strlength(string(rowsAll.bankCenterSource)) > 0, :);
if isempty(rowsAll) || height(rowsAll) == 0
  aggregateTable = struct2table(rowList(:));
  return;
end
rowsAll.primaryDiagnosis = localLookupPrimaryDiagnosis(rowsAll, opt.PrimaryDiagnosisTable);
keyVars = {'displayName', 'snrDb', 'candidateSource', 'bankCenterSource', ...
  'bankAxisDirection', 'bankStepDeg', 'primaryDiagnosis'};
keyTable = unique(rowsAll(:, keyVars), 'rows', 'stable');
for iKey = 1:height(keyTable)
  mask = rowsAll.displayName == keyTable.displayName(iKey) & ...
    rowsAll.snrDb == keyTable.snrDb(iKey) & ...
    rowsAll.candidateSource == keyTable.candidateSource(iKey) & ...
    rowsAll.bankCenterSource == keyTable.bankCenterSource(iKey) & ...
    rowsAll.bankAxisDirection == keyTable.bankAxisDirection(iKey) & ...
    rowsAll.bankStepDeg == keyTable.bankStepDeg(iKey) & ...
    rowsAll.primaryDiagnosis == keyTable.primaryDiagnosis(iKey);
  rows = rowsAll(mask, :);
  row = localEmptyAggregateRow();
  row.displayName = keyTable.displayName(iKey);
  row.snrDb = keyTable.snrDb(iKey);
  row.candidateSource = keyTable.candidateSource(iKey);
  row.bankCenterSource = keyTable.bankCenterSource(iKey);
  row.bankAxisDirection = keyTable.bankAxisDirection(iKey);
  row.bankStepDeg = keyTable.bankStepDeg(iKey);
  row.primaryDiagnosis = keyTable.primaryDiagnosis(iKey);
  row.numCandidate = height(rows);
  row.numUniqueCase = height(unique(rows(:, {'displayName', 'snrDb', 'taskSeed'}), 'rows'));
  row.beatBaselineCount = nnz(rows.wouldBeatBaselineFlag);
  row.beatBaselineRate = localSafeRatio(row.beatBaselineCount, row.numCandidate);
  row.angleImproveCount = nnz(rows.angleImproveDeg > 0);
  row.angleImproveRate = localSafeRatio(row.angleImproveCount, row.numCandidate);
  row.damageCount = nnz(rows.damageFlag);
  row.damageRate = localSafeRatio(row.damageCount, row.numCandidate);
  row.medianObjDelta = median(rows.objMinusBaseline, 'omitnan');
  row.medianAngleImproveDeg = median(rows.angleImproveDeg, 'omitnan');
  row.medianAngleImproveOverCrb = median(rows.angleImproveOverCrb, 'omitnan');
  row.medianCandidateAngleDeg = median(rows.finalAngleErrDeg, 'omitnan');
  row.medianNonRefCoherenceFloor = localMedianIfPresent(rows, 'nonRefCoherenceFloor');
  row.medianFdRefErrOverTooth = localMedianIfPresent(rows, 'finalFdRefErrOverTooth');
  row.medianFdToothShift = localMedianIfPresent(rows, 'candidateToBaselineFdToothShift');
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
aggregateTable = sortrows(struct2table(rowList(:)), ...
  {'displayName', 'snrDb', 'candidateSource', 'primaryDiagnosis', ...
  'bankCenterSource', 'bankAxisDirection', 'bankStepDeg'});
end

function opt = localParseOpt(varargin)
opt = struct('CandidateSourceList', ["ms-bank-cpu-release"], 'PrimaryDiagnosisTable', table());
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankFamilyAggregateTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = lower(string(varargin{iArg}));
  value = varargin{iArg + 1};
  switch name
    case "candidatesourcelist"
      opt.CandidateSourceList = reshape(string(value), 1, []);
    case "primarydiagnosistable"
      opt.PrimaryDiagnosisTable = value;
    otherwise
      error('buildMfBankFamilyAggregateTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function primaryDiagnosis = localLookupPrimaryDiagnosis(rowTable, primaryDiagnosisTable)
primaryDiagnosis = repmat("all", height(rowTable), 1);
if isempty(primaryDiagnosisTable) || ~istable(primaryDiagnosisTable) || height(primaryDiagnosisTable) == 0
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed', 'primaryDiagnosis'};
if ~all(ismember(requiredVars, primaryDiagnosisTable.Properties.VariableNames))
  return;
end
for iRow = 1:height(rowTable)
  idx = find(primaryDiagnosisTable.displayName == rowTable.displayName(iRow) & ...
    primaryDiagnosisTable.snrDb == rowTable.snrDb(iRow) & ...
    primaryDiagnosisTable.taskSeed == rowTable.taskSeed(iRow), 1, 'first');
  if ~isempty(idx)
    primaryDiagnosis(iRow) = primaryDiagnosisTable.primaryDiagnosis(idx);
  end
end
end

function value = localMedianIfPresent(rowTable, varName)
value = NaN;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames)
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
row = struct('displayName', "", 'snrDb', NaN, 'candidateSource', "", ...
  'bankCenterSource', "", 'bankAxisDirection', "", 'bankStepDeg', NaN, ...
  'primaryDiagnosis', "", 'numCandidate', NaN, 'numUniqueCase', NaN, ...
  'beatBaselineCount', NaN, 'beatBaselineRate', NaN, ...
  'angleImproveCount', NaN, 'angleImproveRate', NaN, ...
  'damageCount', NaN, 'damageRate', NaN, ...
  'medianObjDelta', NaN, 'medianAngleImproveDeg', NaN, ...
  'medianAngleImproveOverCrb', NaN, 'medianCandidateAngleDeg', NaN, ...
  'medianNonRefCoherenceFloor', NaN, 'medianFdRefErrOverTooth', NaN, ...
  'medianFdToothShift', NaN);
end
