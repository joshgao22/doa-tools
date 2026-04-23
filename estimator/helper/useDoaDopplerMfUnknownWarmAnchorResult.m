function useAlt = useDoaDopplerMfUnknownWarmAnchorResult(model, solveAlt, solveBase)
%USEDOADOPPLERMFUNKNOWNWARMANCHORRESULT Prefer the explicit CP-U warm-anchor path.

arguments
  model (1,1) struct
  solveAlt struct
  solveBase struct
end

useAlt = preferDoaDopplerMfSolveResult(solveAlt, solveBase);
if useAlt
  return;
end
if ~isfield(model, 'freezeDoa') || ~logical(model.freezeDoa)
  return;
end
if ~isstruct(solveAlt) || isempty(solveAlt)
  return;
end
if ~isstruct(solveBase) || isempty(solveBase)
  useAlt = true;
  return;
end
if solveAlt.isResolved && solveBase.isResolved
  useAlt = true;
end
end
