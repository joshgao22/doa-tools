function plotDoaDopplerComparison(timeOffsetSec, truth, caseResult)
%PLOTDOADOPPLERCOMPARISON Plot reference-Doppler comparisons.
% This compact helper centralizes the reference-Doppler plot used by the
% static/dynamic development scripts.

timeMs = 1e3 * reshape(timeOffsetSec, 1, []);
caseName = strings(1, numel(caseResult));
for iCase = 1:numel(caseResult)
  caseName(iCase) = string(caseResult(iCase).displayName);
end

figure();
plot(timeMs, truth.fdRefSeries, '-ko', 'LineWidth', 1.3); hold on;

plotCaseList = ["SS-SF-Static", "MS-SF-Static", "SS-MF-CP-K", ...
  "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"];
legendCell = {'truth'};
for iCase = 1:numel(plotCaseList)
  idx = find(caseName == plotCaseList(iCase), 1, 'first');
  if isempty(idx)
    continue;
  end
  estResult = caseResult(idx).estResult;
  fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
  fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
  if ~isfinite(fdRefEst)
    continue;
  end
  if ~isfinite(fdRateEst)
    fdRateEst = 0;
  end
  fdLine = fdRefEst + fdRateEst * reshape(timeOffsetSec, 1, []);
  plot(timeMs, fdLine, 'LineWidth', 1.1);
  legendCell{end + 1} = char(plotCaseList(iCase)); %#ok<AGROW>
end

grid on;
xlabel('Time offset (ms)');
ylabel('Reference Doppler (Hz)');
legend(legendCell, 'Location', 'best');
title('Reference-Doppler comparison');
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct/object field with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end
