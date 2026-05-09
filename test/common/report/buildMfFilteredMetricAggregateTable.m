function filteredTable = buildMfFilteredMetricAggregateTable(caseTable, filterName, keepMask)
%BUILDMFFILTEREDMETRICAGGREGATETABLE Aggregate MF metrics after a row filter.
% Used by replay and scan scripts for health/resolved/trimmed summaries.

rowTemplate = localEmptyFilteredAggregateRow();
if isempty(caseTable) || height(caseTable) == 0
  filteredTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end

keepMask = logical(keepMask(:));
groupKey = unique(caseTable(:, {'displayName', 'frameMode', 'fdRateMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iRow = 1:height(groupKey)
  groupMask = caseTable.displayName == groupKey.displayName(iRow) & caseTable.snrDb == groupKey.snrDb(iRow);
  mask = groupMask & keepMask;
  row = rowTemplate;
  row.filterName = string(filterName);
  row.displayName = groupKey.displayName(iRow);
  row.frameMode = groupKey.frameMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.rawCount = nnz(groupMask);
  row.keepCount = nnz(mask);
  row.keepRate = localSafeRatio(row.keepCount, row.rawCount);
  row.rejectRate = 1 - row.keepRate;
  row.numRepeat = row.keepCount;
  if row.keepCount > 0
    row = localFillMetricFields(row, caseTable, mask);
  end
  rowList(iRow) = row;
end
filteredTable = struct2table(rowList);
end

function row = localFillMetricFields(row, caseTable, mask)
%LOCALFILLMETRICFIELDS Fill common metric fields for one filtered group.
row.angleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(mask));
row.angleMedianDeg = median(caseTable.sphericalAngleErrDeg(mask), 'omitnan');
row.angleMaxDeg = max(caseTable.sphericalAngleErrDeg(mask), [], 'omitnan');
row.angleCrbSphericalMedianDeg = median(caseTable.crbSphericalApproxStdDeg(mask), 'omitnan');
row.angleCrbTraceMedianDeg = median(caseTable.crbTraceStdDeg(mask), 'omitnan');
row.angleRmseOverSphericalCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbSphericalMedianDeg);
row.angleRmseOverTraceCrb = localSafeRatio(row.angleRmseDeg, row.angleCrbTraceMedianDeg);
row.fdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(mask));
row.fdRefCrbMedianHz = median(caseTable.fdRefCrbStdHz(mask), 'omitnan');
row.fdRefRmseOverCrb = localSafeRatio(row.fdRefRmseHz, row.fdRefCrbMedianHz);
row.fdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(mask));
row.initAngleMedianDeg = median(caseTable.initAngleErrDeg(mask), 'omitnan');
row.finalMinusInitMedianDeg = median(caseTable.finalMinusInitAngleDeg(mask), 'omitnan');
row.fdRefInitMoveMedianHz = median(abs(caseTable.fdRefInitMoveHz(mask)), 'omitnan');
row.fdRateInitMoveMedianHzPerSec = median(abs(caseTable.fdRateInitMoveHzPerSec(mask)), 'omitnan');
row.objectiveImproveMedian = median(caseTable.objectiveImprove(mask), 'omitnan');
row.boundaryHitRate = mean(double(caseTable.fdRefBoundaryHit(mask) | caseTable.fdRateBoundaryHit(mask)), 'omitnan');
row.candidateCountMedian = median(caseTable.candidateCount(mask), 'omitnan');
row.iterationsMedian = median(caseTable.iterations(mask), 'omitnan');
row.firstOrderOptMedian = median(caseTable.firstOrderOpt(mask), 'omitnan');
end

function row = localEmptyFilteredAggregateRow()
%LOCALEMPTYFILTEREDAGGREGATEROW Return one typed filtered aggregate row.
row = struct('displayName', "", 'frameMode', "", 'fdRateMode', "", 'snrDb', NaN, ...
  'numRepeat', NaN, 'angleRmseDeg', NaN, 'angleMedianDeg', NaN, ...
  'angleMaxDeg', NaN, 'angleCrbSphericalMedianDeg', NaN, ...
  'angleCrbTraceMedianDeg', NaN, 'angleRmseOverSphericalCrb', NaN, ...
  'angleRmseOverTraceCrb', NaN, 'fdRefRmseHz', NaN, ...
  'fdRefCrbMedianHz', NaN, 'fdRefRmseOverCrb', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'initAngleMedianDeg', NaN, ...
  'finalMinusInitMedianDeg', NaN, 'fdRefInitMoveMedianHz', NaN, ...
  'fdRateInitMoveMedianHzPerSec', NaN, 'objectiveImproveMedian', NaN, ...
  'boundaryHitRate', NaN, 'candidateCountMedian', NaN, ...
  'iterationsMedian', NaN, 'firstOrderOptMedian', NaN, 'filterName', "", ...
  'rawCount', NaN, 'keepCount', NaN, 'keepRate', NaN, 'rejectRate', NaN);
end

function value = localRmse(x)
%LOCALRMSE Root-mean-square of finite entries.
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = sqrt(mean(x(:).^2));
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide while protecting empty/zero denominators.
if isfinite(numerator) && isfinite(denominator) && denominator > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end
