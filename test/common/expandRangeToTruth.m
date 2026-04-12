function rangeOut = expandRangeToTruth(rangeIn, truthVal, marginFrac, minMarginVal)
%EXPANDRANGETOTRUTH Expand one search range so truth stays inside.

arguments
  rangeIn (1, 2) double
  truthVal (:, 1) double
  marginFrac (1, 1) double = 0.1
  minMarginVal (1, 1) double = 0
end

rangeOut = rangeIn;
truthVal = truthVal(isfinite(truthVal));
if isempty(truthVal)
  return;
end

marginVal = max(minMarginVal, marginFrac * max(abs(truthVal)));
rangeOut(1) = min(rangeOut(1), min(truthVal) - marginVal);
rangeOut(2) = max(rangeOut(2), max(truthVal) + marginVal);
end
