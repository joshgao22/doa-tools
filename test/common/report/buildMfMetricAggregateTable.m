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
row.angleMseDeg2 = localMse(caseTable.sphericalAngleErrDeg(mask));
row.angleRmseDeg = sqrt(row.angleMseDeg2);
row.angleMedianDeg = median(caseTable.sphericalAngleErrDeg(mask), 'omitnan');
row.angleMaxDeg = max(caseTable.sphericalAngleErrDeg(mask), [], 'omitnan');
row.angleCrbSphericalMedianDeg = median(caseTable.crbSphericalApproxStdDeg(mask), 'omitnan');
row.angleCrbTraceMedianDeg = median(caseTable.crbTraceStdDeg(mask), 'omitnan');
row.angleMseOverSphericalCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(mask), ...
  caseTable.crbSphericalApproxStdDeg(mask));
row.angleMseOverTraceCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(mask), ...
  caseTable.crbTraceStdDeg(mask));
row.angleRmseOverSphericalCrb = sqrt(row.angleMseOverSphericalCrb);
row.angleRmseOverTraceCrb = sqrt(row.angleMseOverTraceCrb);
row.fdRefMseHz2 = localMse(caseTable.fdRefAbsErrHz(mask));
row.fdRefRmseHz = sqrt(row.fdRefMseHz2);
row.fdRefCrbMedianHz = median(caseTable.fdRefCrbStdHz(mask), 'omitnan');
row.fdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(mask), ...
  caseTable.fdRefCrbStdHz(mask));
row.fdRefRmseOverCrb = sqrt(row.fdRefMseOverCrb);
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
  'numRepeat', NaN, 'angleMseDeg2', NaN, 'angleRmseDeg', NaN, ...
  'angleMedianDeg', NaN, 'angleMaxDeg', NaN, ...
  'angleCrbSphericalMedianDeg', NaN, 'angleCrbTraceMedianDeg', NaN, ...
  'angleRmseOverSphericalCrb', NaN, 'angleRmseOverTraceCrb', NaN, ...
  'angleMseOverSphericalCrb', NaN, 'angleMseOverTraceCrb', NaN, ...
  'fdRefMseHz2', NaN, 'fdRefRmseHz', NaN, ...
  'fdRefCrbMedianHz', NaN, 'fdRefRmseOverCrb', NaN, ...
  'fdRefMseOverCrb', NaN, ...
  'fdRateRmseHzPerSec', NaN, 'initAngleMedianDeg', NaN, ...
  'finalMinusInitMedianDeg', NaN, 'fdRefInitMoveMedianHz', NaN, ...
  'fdRateInitMoveMedianHzPerSec', NaN, 'objectiveImproveMedian', NaN, ...
  'boundaryHitRate', NaN, 'candidateCountMedian', NaN, ...
  'iterationsMedian', NaN, 'firstOrderOptMedian', NaN);
end

function value = localMse(x)
%LOCALMSE Mean-square of finite entries.
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = mean(x(:).^2);
end
end

function value = localRmse(x)
%LOCALRMSE Root-mean-square of finite entries.
value = sqrt(localMse(x));
end

function value = localNormalizedMse(err, crbStd)
%LOCALNORMALIZEDMSE Mean squared CRB-normalized error over finite entries.
err = reshape(double(err), [], 1);
crbStd = reshape(double(crbStd), [], 1);
mask = isfinite(err) & isfinite(crbStd) & crbStd > 0;
if ~any(mask)
  value = NaN;
else
  value = mean((err(mask) ./ crbStd(mask)).^2);
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
