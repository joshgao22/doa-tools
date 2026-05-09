function aggregateTable = buildMfWideCoverageAggregateTable(candidateTraceTable, varargin)
%BUILDMFWIDECOVERAGEAGGREGATETABLE Aggregate static-wide or bank candidate utility.
%
% The function only summarizes existing candidate trace rows. It does not
% select a final winner or define an estimator adoption policy.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyWideCoverageAggregateRow(), 0, 1);
if isempty(candidateTraceTable) || ~istable(candidateTraceTable) || height(candidateTraceTable) == 0
  aggregateTable = struct2table(rowList(:));
  return;
end
candidateSource = string(opt.CandidateSource);
wideMask = candidateTraceTable.candidateSource == candidateSource & isfinite(candidateTraceTable.doaHalfWidthDeg);
wideRows = candidateTraceTable(wideMask, :);
if isempty(wideRows) || height(wideRows) == 0
  aggregateTable = struct2table(rowList(:));
  return;
end
groupKey = unique(wideRows(:, {'displayName', 'snrDb', 'doaHalfWidthDeg'}), 'rows', 'stable');
for iGroup = 1:height(groupKey)
  mask = wideRows.displayName == groupKey.displayName(iGroup) & ...
    wideRows.snrDb == groupKey.snrDb(iGroup) & ...
    wideRows.doaHalfWidthDeg == groupKey.doaHalfWidthDeg(iGroup);
  rows = wideRows(mask, :);
  row = localEmptyWideCoverageAggregateRow();
  row.displayName = groupKey.displayName(iGroup);
  row.snrDb = groupKey.snrDb(iGroup);
  row.halfWidthDeg = groupKey.doaHalfWidthDeg(iGroup);
  row.numCase = height(rows);
  row.beatBaselineCount = nnz(rows.wouldBeatBaselineFlag);
  row.beatBaselineRate = localSafeRatio(row.beatBaselineCount, row.numCase);
  row.angleImproveCount = nnz(rows.angleImproveDeg > 0);
  row.angleImproveRate = localSafeRatio(row.angleImproveCount, row.numCase);
  row.damageCount = nnz(rows.damageFlag);
  row.damageRate = localSafeRatio(row.damageCount, row.numCase);
  row.medianObjDelta = median(rows.objMinusBaseline, 'omitnan');
  row.medianAngleImproveDeg = median(rows.angleImproveDeg, 'omitnan');
  row.medianAngleImproveOverCrb = median(rows.angleImproveOverCrb, 'omitnan');
  row.medianCandidateAngleDeg = median(rows.finalAngleErrDeg, 'omitnan');
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
aggregateTable = sortrows(struct2table(rowList(:)), {'displayName', 'snrDb', 'halfWidthDeg'});
end

function opt = localParseOpt(varargin)
opt = struct('CandidateSource', "static-wide-solve");
if mod(numel(varargin), 2) ~= 0
  error('buildMfWideCoverageAggregateTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = string(varargin{iArg});
  value = varargin{iArg + 1};
  switch lower(name)
    case "candidatesource"
      opt.CandidateSource = string(value);
    otherwise
      error('buildMfWideCoverageAggregateTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function row = localEmptyWideCoverageAggregateRow()
row = struct('displayName', "", 'snrDb', NaN, 'halfWidthDeg', NaN, ...
  'numCase', NaN, 'beatBaselineCount', NaN, 'beatBaselineRate', NaN, ...
  'angleImproveCount', NaN, 'angleImproveRate', NaN, ...
  'damageCount', NaN, 'damageRate', NaN, 'medianObjDelta', NaN, ...
  'medianAngleImproveDeg', NaN, 'medianAngleImproveOverCrb', NaN, ...
  'medianCandidateAngleDeg', NaN);
end
