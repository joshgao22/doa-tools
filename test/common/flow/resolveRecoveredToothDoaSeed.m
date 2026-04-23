function [doaInitParam, seedSource] = resolveRecoveredToothDoaSeed(flowOpt, selectedSubsetCase, selectedSubsetSummary, knownSummary)
%RESOLVERECOVEREDTOOTHDOASEED Choose one DoA seed for recovered-tooth routes.
% Random rescue schedules only select the tooth. When one random subset wins
% tooth 0 but still drifts materially from the CP-K DoA anchor, keep the
% rescued fd/fdRate neighborhood while recentering the DoA seed back to the
% known-rate basin. This helper is shared by subset-anchor and in-tooth
% routes so recovered-tooth seeding semantics stay aligned.

arguments
  flowOpt (1,1) struct
  selectedSubsetCase (1,1) struct
  selectedSubsetSummary (1,1) struct
  knownSummary (1,1) struct
end

doaInitParam = reshape(selectedSubsetCase.estResult.doaParamEst(:), [], 1);
seedSource = "selected-subset";
if ~logical(localGetFieldOrDefault(flowOpt, 'enableRandomSubsetKnownDoaRecentering', true))
  return;
end
subsetLabel = lower(strtrim(char(string(localGetFieldOrDefault(selectedSubsetSummary, 'subsetLabel', "")))));
randomSelected = startsWith(subsetLabel, 'random') || startsWith(subsetLabel, 'rescue-random');
if ~randomSelected
  return;
end
knownDoaParam = reshape(localGetFieldOrDefault(knownSummary, 'doaParamEst', []), [], 1);
if isempty(knownDoaParam) || any(~isfinite(knownDoaParam)) || ...
    ~logical(localGetFieldOrDefault(knownSummary, 'isResolved', false))
  return;
end
minDriftDeg = localGetFieldOrDefault(flowOpt, 'randomSubsetKnownDoaRecenteringMinDriftDeg', 0.003);
doaDriftDeg = localCalcSummaryDoaDriftDeg(selectedSubsetSummary, knownSummary);
if ~(isfinite(doaDriftDeg) && doaDriftDeg >= minDriftDeg)
  return;
end
doaInitParam = knownDoaParam;
seedSource = "known-summary";
end


function doaDriftDeg = localCalcSummaryDoaDriftDeg(summaryA, summaryB)
%LOCALCALCSUMMARYDOADRIFTDEG Build one DoA drift metric between compact summaries.

doaDriftDeg = NaN;
doaA = localGetFieldOrDefault(summaryA, 'doaParamEst', []);
doaB = localGetFieldOrDefault(summaryB, 'doaParamEst', []);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = NaN;
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or a default value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
