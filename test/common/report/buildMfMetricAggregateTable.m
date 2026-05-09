function aggregateTable = buildMfMetricAggregateTable(caseTable)
%BUILDMFMETRICAGGREGATETABLE Aggregate common MF replay/scan metric rows.
% The helper expects the compact case-table columns produced by MF replay or
% scan scripts. It is reporting-only and does not change estimator behavior.

rowTemplate = localEmptyAggregateRow();
if isempty(caseTable) || height(caseTable) == 0
  aggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end

groupKey = unique(caseTable(:, {'displayName', 'frameMode', 'fdRateMode', 'snrDb'}), 'rows', 'stable');
rowList = repmat(rowTemplate, height(groupKey), 1);
for iRow = 1:height(groupKey)
  mask = caseTable.displayName == groupKey.displayName(iRow) & caseTable.snrDb == groupKey.snrDb(iRow);
  row = localBuildOneMetricAggregateRow(caseTable, mask, rowTemplate);
  row.displayName = groupKey.displayName(iRow);
  row.frameMode = groupKey.frameMode(iRow);
  row.fdRateMode = groupKey.fdRateMode(iRow);
  row.snrDb = groupKey.snrDb(iRow);
  row.numRepeat = nnz(mask);
  rowList(iRow) = row;
end
aggregateTable = struct2table(rowList);
end

function row = localBuildOneMetricAggregateRow(caseTable, mask, row)
%LOCALBUILDONEMETRICAGGREGATEROW Fill common metric fields for one group.
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

function row = localEmptyAggregateRow()
%LOCALEMPTYAGGREGATEROW Return one typed metric aggregate row.
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
  'iterationsMedian', NaN, 'firstOrderOptMedian', NaN);
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
