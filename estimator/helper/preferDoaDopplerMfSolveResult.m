function useAlt = preferDoaDopplerMfSolveResult(solveAlt, solveBase)
%PREFERDOADOPPLERMFSOLVERESULT Select the better MF solve result conservatively.

arguments
  solveAlt struct
  solveBase struct
end

useAlt = false;
if ~isstruct(solveAlt) || isempty(solveAlt)
  return;
end
if ~isstruct(solveBase) || isempty(solveBase)
  useAlt = true;
  return;
end

if solveAlt.isResolved && ~solveBase.isResolved
  useAlt = true;
  return;
end
if ~solveAlt.isResolved && solveBase.isResolved
  return;
end

baseObj = solveBase.fval;
altObj = solveAlt.fval;
if ~(isfinite(altObj) && isfinite(baseObj))
  useAlt = isfinite(altObj) && ~isfinite(baseObj);
  return;
end

objTol = 1e-9 * max([1, abs(baseObj), abs(altObj)]);
if altObj < baseObj - objTol
  useAlt = true;
end
end
