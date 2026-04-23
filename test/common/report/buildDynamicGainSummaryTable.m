function gainTable = buildDynamicGainSummaryTable(angleErrDeg, fdErrHz, fdRateErrHzPerSec, isResolved, ...
  idxSsDynKnown, idxMsDynKnown, idxSsDynUnknown, idxMsDynUnknown)
%BUILDDYNAMICGAINSUMMARYTABLE Build one compact SS/MS dynamic gain table.

arguments
  angleErrDeg (:, :) double
  fdErrHz (:, :) double
  fdRateErrHzPerSec (:, :) double
  isResolved (:, :) logical
  idxSsDynKnown (1, 1) double {mustBeInteger, mustBePositive}
  idxMsDynKnown (1, 1) double {mustBeInteger, mustBePositive}
  idxSsDynUnknown (1, 1) double {mustBeInteger, mustBePositive}
  idxMsDynUnknown (1, 1) double {mustBeInteger, mustBePositive}
end

metricName = strings(0, 1);
metricValue = zeros(0, 1);

[metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
  "numRepeat", size(angleErrDeg, 1));
[metricName, metricValue] = localAppendPairGain(metricName, metricValue, ...
  "known", angleErrDeg(:, idxSsDynKnown), angleErrDeg(:, idxMsDynKnown), ...
  fdErrHz(:, idxSsDynKnown), fdErrHz(:, idxMsDynKnown), ...
  nan(size(fdRateErrHzPerSec, 1), 1), nan(size(fdRateErrHzPerSec, 1), 1), ...
  isResolved(:, idxSsDynKnown), isResolved(:, idxMsDynKnown));
[metricName, metricValue] = localAppendPairGain(metricName, metricValue, ...
  "unknown", angleErrDeg(:, idxSsDynUnknown), angleErrDeg(:, idxMsDynUnknown), ...
  fdErrHz(:, idxSsDynUnknown), fdErrHz(:, idxMsDynUnknown), ...
  fdRateErrHzPerSec(:, idxSsDynUnknown), fdRateErrHzPerSec(:, idxMsDynUnknown), ...
  isResolved(:, idxSsDynUnknown), isResolved(:, idxMsDynUnknown));

gainTable = table(metricName, metricValue, ...
  'VariableNames', {'metricName', 'metricValue'});
end


function [metricName, metricValue] = localAppendPairGain(metricName, metricValue, tagName, ...
  ssAngle, msAngle, ssFd, msFd, ssFdRate, msFdRate, ssResolved, msResolved)
%LOCALAPPENDPAIRGAIN Append one SS/MS compare block.

maskBase = ssResolved & msResolved & isfinite(ssAngle) & isfinite(msAngle);
maskFd = maskBase & isfinite(ssFd) & isfinite(msFd);
maskFdRate = maskBase & isfinite(ssFdRate) & isfinite(msFdRate);
angleDiff = msAngle - ssAngle;
fdDiff = msFd - ssFd;
fdRateDiff = msFdRate - ssFdRate;

[metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
  "msBetterAngleRate_" + tagName, localMeanLogical(angleDiff(maskBase) < 0));
[metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
  "msBetterFdRate_" + tagName, localMeanLogical(fdDiff(maskFd) < 0));
[metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
  "medianMsMinusSsAngleDeg_" + tagName, median(angleDiff(maskBase), 'omitnan'));
[metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
  "medianMsMinusSsFdHz_" + tagName, median(fdDiff(maskFd), 'omitnan'));
if any(maskFdRate)
  [metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
    "msBetterFdRateParamRate_" + tagName, localMeanLogical(fdRateDiff(maskFdRate) < 0));
  [metricName, metricValue] = localAppendGainMetric(metricName, metricValue, ...
    "medianMsMinusSsFdRateHzPerSec_" + tagName, median(fdRateDiff(maskFdRate), 'omitnan'));
end
end

function [metricName, metricValue] = localAppendGainMetric(metricName, metricValue, nameUse, valueUse)
metricName(end + 1, 1) = string(nameUse);
metricValue(end + 1, 1) = valueUse;
end

function value = localMeanLogical(mask)
%LOCALMEANLOGICAL Return mean(mask) with empty handling.

mask = logical(mask);
if isempty(mask)
  value = NaN;
else
  value = mean(mask);
end
end
