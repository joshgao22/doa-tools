function textValue = formatMfReplaySeedList(seedList)
%FORMATMFREPLAYSEEDLIST Format seed lists for compact replay logs.
%
% Consecutive integer seed lists are shown as a colon range. Non-consecutive
% lists are printed in full so manually selected seeds remain visible.

seedList = reshape(double(seedList), 1, []);
if isempty(seedList)
  textValue = "";
  return;
end

if localIsConsecutiveIntegerList(seedList)
  textValue = sprintf('%.0f:%.0f (%d seeds)', seedList(1), seedList(end), numel(seedList));
  return;
end

textValue = strjoin(compose('%.6g', seedList), ', ');
end

function tf = localIsConsecutiveIntegerList(valueList)
%LOCALISCONSECUTIVEINTEGERLIST Return true for integer unit-stride vectors.

valueList = reshape(double(valueList), 1, []);
tf = numel(valueList) > 1 ...
  && all(isfinite(valueList)) ...
  && all(abs(valueList - round(valueList)) < sqrt(eps)) ...
  && all(diff(round(valueList)) == 1);
end
