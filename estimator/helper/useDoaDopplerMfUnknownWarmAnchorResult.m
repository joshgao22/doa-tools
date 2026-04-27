function useAlt = useDoaDopplerMfUnknownWarmAnchorResult(model, solveAlt, solveBase)
%USEDOADOPPLERMFUNKNOWNWARMANCHORRESULT Prefer the explicit CP-U warm-anchor path.

arguments
  model (1,1) struct
  solveAlt struct
  solveBase struct
end

verbose = localResolveVerbose(model);
freezeDoa = localGetModelLogical(model, 'freezeDoa', false);
preferAlt = preferDoaDopplerMfSolveResult(solveAlt, solveBase);
useAlt = preferAlt;

if verbose
  fprintf('[WarmAnchorAdopt] freezeDoa=%d preferRule=%d\n', freezeDoa, preferAlt);
  localPrintCandidateSummary('base', solveBase);
  localPrintCandidateSummary('alt', solveAlt);
end

if useAlt
  if verbose
    fprintf('  decision: keep alt because preferDoaDopplerMfSolveResult returned true.\n');
  end
  return;
end
if ~freezeDoa
  if verbose
    fprintf('  decision: keep base because freezeDoa is off.\n');
  end
  return;
end

isAltFixedDoa = localIsFixedDoaWarmAnchorResult(solveAlt);
if ~isAltFixedDoa
  if verbose
    fprintf('  decision: keep base because alt is not fixed-DoA warm-anchor.\n');
  end
  return;
end
if ~isstruct(solveBase) || isempty(solveBase)
  useAlt = true;
  if verbose
    fprintf('  decision: take alt because base is empty.\n');
  end
  return;
end
if solveAlt.isResolved && ~solveBase.isResolved
  useAlt = true;
  if verbose
    fprintf('  decision: take alt because base is unresolved.\n');
  end
  return;
end
if ~(solveAlt.isResolved && solveBase.isResolved)
  if verbose
    fprintf('  decision: keep base because both candidates are not resolved.\n');
  end
  return;
end

isMateriallyWorse = localIsMateriallyWorse(solveAlt, solveBase);
useAlt = ~isMateriallyWorse;
if verbose
  fprintf('  flags: isAltFixedDoa=%d materiallyWorse=%d => takeAlt=%d\n', ...
    isAltFixedDoa, isMateriallyWorse, useAlt);
end
end


function tf = localIsFixedDoaWarmAnchorResult(solveResult)
%LOCALISFIXEDDOAWARMANCHORRESULT True only for explicit fixed-DoA warm-anchor solves.

tf = false;
if ~isstruct(solveResult) || isempty(solveResult) || ~isfield(solveResult, 'solveVariant')
  return;
end
solveVariant = string(solveResult.solveVariant);
if isempty(solveVariant)
  return;
end

tf = any(solveVariant == ["mainWarmAnchorFixedDoa", "mainWarmFixedDoa"]);
end


function tf = localIsMateriallyWorse(solveAlt, solveBase)
%LOCALISMATERIALLYWORSE True when one fixed-DoA warm-anchor solve is clearly worse.

tf = false;

altObj = localGetScalarField(solveAlt, 'fval', NaN);
baseObj = localGetScalarField(solveBase, 'fval', NaN);
altResidual = localGetEvalDiagScalar(solveAlt, 'residualNorm', NaN);
baseResidual = localGetEvalDiagScalar(solveBase, 'residualNorm', NaN);
altFitFloor = localGetEvalDiagScalar(solveAlt, 'nonRefFitRatioFloor', NaN);
baseFitFloor = localGetEvalDiagScalar(solveBase, 'nonRefFitRatioFloor', NaN);
altSupportFloor = localGetEvalDiagScalar(solveAlt, 'nonRefSupportRatioFloor', NaN);
baseSupportFloor = localGetEvalDiagScalar(solveBase, 'nonRefSupportRatioFloor', NaN);

if isfinite(altObj) && isfinite(baseObj)
  objTol = max(5, 1e-6 * max([1, abs(baseObj), abs(altObj)]));
  if altObj > baseObj + objTol
    tf = true;
    return;
  end
end
if isfinite(altResidual) && isfinite(baseResidual)
  residualTol = max(1, 1e-6 * max([1, abs(baseResidual), abs(altResidual)]));
  if altResidual > baseResidual + residualTol
    tf = true;
    return;
  end
end
qualityTol = 1e-3;
if isfinite(altFitFloor) && isfinite(baseFitFloor) && altFitFloor + qualityTol < baseFitFloor
  tf = true;
  return;
end
if isfinite(altSupportFloor) && isfinite(baseSupportFloor) && altSupportFloor + qualityTol < baseSupportFloor
  tf = true;
end
end


function value = localGetEvalDiagScalar(solveInfo, fieldName, defaultValue)
%LOCALGETEVALDIAGSCALAR Read one scalar from finalEvalDiag with default fallback.

value = defaultValue;
if ~isstruct(solveInfo) || ~isfield(solveInfo, 'finalEvalDiag') || ~isstruct(solveInfo.finalEvalDiag)
  return;
end
value = localGetScalarField(solveInfo.finalEvalDiag, fieldName, defaultValue);
end


function value = localGetScalarField(dataStruct, fieldName, defaultValue)
%LOCALGETSCALARFIELD Read one scalar field with default fallback.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
rawValue = dataStruct.(fieldName);
if isempty(rawValue) || ~isscalar(rawValue)
  return;
end
if isnumeric(rawValue) || islogical(rawValue)
  if ~isfinite(double(rawValue))
    return;
  end
end
value = rawValue;
end


function verbose = localResolveVerbose(model)
%LOCALRESOLVEVERBOSE Resolve one helper-local verbose flag.

verbose = false;
if isstruct(model) && isfield(model, 'verbose') && ~isempty(model.verbose)
  verbose = logical(model.verbose);
end
verbose = verbose || localResolveInheritedVerbose();
end


function value = localGetModelLogical(model, fieldName, defaultValue)
%LOCALGETMODELLOGICAL Read one logical model flag with default fallback.

value = defaultValue;
if ~isstruct(model) || ~isfield(model, fieldName) || isempty(model.(fieldName))
  return;
end
value = logical(model.(fieldName));
end


function localPrintCandidateSummary(tag, solveInfo)
%LOCALPRINTCANDIDATESUMMARY Print one compact adoption-candidate summary.

if isempty(solveInfo) || ~isstruct(solveInfo)
  fprintf('  %-8s : <empty>\n', tag);
  return;
end
fprintf(['  %-8s : variant=%s resolved=%d fval=%.6e ', ...
         'resid=%.6e fit=%.6e support=%.6e\n'], ...
  tag, ...
  string(localGetScalarField(solveInfo, 'solveVariant', "")), ...
  localGetScalarField(solveInfo, 'isResolved', false), ...
  localGetScalarField(solveInfo, 'fval', NaN), ...
  localGetEvalDiagScalar(solveInfo, 'residualNorm', NaN), ...
  localGetEvalDiagScalar(solveInfo, 'nonRefFitRatioFloor', NaN), ...
  localGetEvalDiagScalar(solveInfo, 'nonRefSupportRatioFloor', NaN));
end

function verbose = localResolveInheritedVerbose()
%LOCALRESOLVEINHERITEDVERBOSE Inherit one repo-level verbose flag.

verbose = false;
try
  if isappdata(0, 'doaToolsVerbose')
    verbose = logical(getappdata(0, 'doaToolsVerbose'));
  end
catch
  verbose = false;
end

if verbose
  return;
end

try
  hasVerbose = evalin('base', "exist(''doaToolsVerbose'', ''var'')");
  if isequal(hasVerbose, 1)
    verbose = logical(evalin('base', 'doaToolsVerbose'));
  end
catch
  verbose = false;
end
end

