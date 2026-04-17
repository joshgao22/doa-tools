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
if ~isempty(initParamCpKnown)
  initCand{end + 1} = struct( ...
    'initParam', initParamCpKnown, ...
    'initDoaParam', reshape(cpKnownCase.estResult.doaParamEst, [], 1), ...
    'initDoaHalfWidth', cpKnownHalfWidth, ...
    'startTag', "fromCpK");
end
if ~isempty(initParamStatic)
  initCand{end + 1} = struct( ...
    'initParam', initParamStatic, ...
    'initDoaParam', reshape(staticCase.estResult.doaParamEst, [], 1), ...
    'initDoaHalfWidth', staticHalfWidth, ...
    'startTag', "fromStatic");
end
if isempty(initCand)
  initCand = [];
end
end
