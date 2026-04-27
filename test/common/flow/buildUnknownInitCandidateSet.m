function initCand = buildUnknownInitCandidateSet(cpKnownCase, staticCase, ...
  initParamCpKnown, initParamStatic, cpKnownHalfWidth, staticHalfWidth)
%BUILDUNKNOWNINITCANDIDATESET Build dual warm starts for CP-U release.
% The CP-K continuation should stay inside a tight local DoA box, while
% the static-seed release can reuse the slightly wider static box.

arguments
  cpKnownCase (1, 1) struct
  staticCase (1, 1) struct
  initParamCpKnown = []
  initParamStatic = []
  cpKnownHalfWidth (:, 1) double = []
  staticHalfWidth (:, 1) double = []
end

if isempty(staticHalfWidth)
  staticHalfWidth = cpKnownHalfWidth;
end

initCand = {};
initCand = localAppendInitCandidate(initCand, cpKnownCase, initParamCpKnown, cpKnownHalfWidth, "fromCpK");
initCand = localAppendInitCandidate(initCand, staticCase, initParamStatic, staticHalfWidth, "fromStatic");
if isempty(initCand)
  initCand = [];
end
end


function initCand = localAppendInitCandidate(initCand, sourceCase, initParamUse, initDoaHalfWidth, startTag)
%LOCALAPPENDINITCANDIDATE Append one well-formed outer warm-start candidate.

if isempty(initParamUse) || isempty(initDoaHalfWidth)
  return;
end
estResult = localGetFieldOrDefault(sourceCase, 'estResult', struct());
if isempty(estResult) || ~isstruct(estResult)
  return;
end
initDoaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if isempty(initDoaParam) || ~all(isfinite(initDoaParam))
  return;
end
initParamUse = reshape(initParamUse, [], 1);
initDoaHalfWidth = reshape(initDoaHalfWidth, [], 1);
if ~all(isfinite(initParamUse)) || ~all(isfinite(initDoaHalfWidth)) || any(initDoaHalfWidth <= 0)
  return;
end

initCand{end + 1} = struct( ...
  'initParam', initParamUse, ...
  'initDoaParam', initDoaParam, ...
  'initDoaHalfWidth', initDoaHalfWidth, ...
  'startTag', string(startTag), ...
  'sourceSolveVariant', string(localGetFieldOrDefault(estResult, 'solveVariant', "unknown")), ...
  'sourceIsResolved', logical(localGetFieldOrDefault(estResult, 'isResolved', false)));
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct field with a default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
