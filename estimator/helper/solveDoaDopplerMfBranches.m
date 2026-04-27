function [solveUse, optimInfo, initEvalDiag] = solveDoaDopplerMfBranches(model, initParam, verbose)
%SOLVEDOADOPPLERMFBRANCHES Solve the multi-frame estimator branches.
% Runs the standard solve and, when requested by fdRateMode, the known-rate
% embedded branch or the unknown-rate continuation candidates. The returned
% solveUse is the selected best branch using the same preference rule as
% estimatorDoaDopplerMlePilotMfOpt.

arguments
  model (1,1) struct
  initParam (:,1) double
  verbose (1,1) logical = false
end

verbose = localResolveVerbose(model, verbose);
model.verbose = verbose;

[lb, ub] = buildDoaDopplerMfBounds(model);
[initObj, ~, initNoiseVar, initEvalAux] = evaluateDoaDopplerMfObjective(model, initParam);
initEvalDiag = buildDoaDopplerMfEvalDiag(model, initParam, initObj, initNoiseVar, initEvalAux);

displayMode = 'off';
if verbose && localGetModelLogical(model, 'verboseSolverIterations', false)
  displayMode = 'iter';
end

baseOptimOpt = optimoptions('fmincon', ...
  'Display', displayMode, ...
  'Algorithm', 'interior-point');
baseOptimField = fieldnames(model.optimOpt);
for iField = 1:numel(baseOptimField)
  baseOptimOpt.(baseOptimField{iField}) = model.optimOpt.(baseOptimField{iField});
end

solveStdOpt = struct('useScaledSolve', false, 'enableFallbackSqp', false);
solveStd = runDoaDopplerMfOptimization(model, initParam, lb, ub, [], [], baseOptimOpt, solveStdOpt);
solveUse = solveStd;
if verbose
  fprintf('[MfBranchTrace] mode=%s/%s initObj=%.6e freezeDoa=%d\n', ...
    string(model.phaseMode), string(model.fdRateMode), initObj, localGetModelLogical(model, 'freezeDoa', false));
  localPrintSolveSummary('std', solveStd);
end

if strcmp(model.fdRateMode, 'known')
  solveEmbed = runDoaDopplerMfKnownEmbeddedOptimization(model, initParam, baseOptimOpt);
  if verbose
    localPrintSolveSummary('knownEmbed', solveEmbed);
  end
  if preferDoaDopplerMfSolveResult(solveEmbed, solveStd)
    solveUse = solveEmbed;
    if verbose
      fprintf('[MfBranchTrace] selected known embedded branch.\n');
    end
  elseif verbose
    fprintf('[MfBranchTrace] kept standard known branch.\n');
  end
elseif strcmp(model.fdRateMode, 'unknown')
  if isfield(model, 'disableUnknownWarmAnchor') && logical(model.disableUnknownWarmAnchor)
    solveCand = struct([]);
    if verbose
      fprintf('[MfBranchTrace] unknown warm-anchor disabled.\n');
    end
  else
    solveCand = runDoaDopplerMfUnknownWarmAnchor(model, initParam, lb, ub, baseOptimOpt);
  end
  if ~isempty(solveCand)
    candidateObjective = solveCand.fval;
    candidateExitflag = solveCand.exitflag;
    candidateVariant = string(solveCand.solveVariant);
    if verbose
      localPrintSolveSummary('unknownCand', solveCand);
    end
    if useDoaDopplerMfUnknownWarmAnchorResult(model, solveCand, solveUse)
      solveUse = solveCand;
      if verbose
        fprintf('[MfBranchTrace] adopted warm-anchor candidate.\n');
      end
    elseif verbose
      fprintf('[MfBranchTrace] kept base solve over warm-anchor candidate.\n');
    end
    solveUse.candidateObjective = candidateObjective;
    solveUse.candidateExitflag = candidateExitflag;
    solveUse.candidateVariant = candidateVariant;
  elseif verbose
    fprintf('[MfBranchTrace] no warm-anchor candidate generated.\n');
  end
end

optimInfo = localBuildOptimInfo(solveUse.output, solveUse.exitflag, ...
  solveUse.runTimeSec, baseOptimOpt, model);
optimInfo.solveVariant = solveUse.solveVariant;
optimInfo.verbose = verbose;
if isfield(solveUse, 'candidateObjective')
  optimInfo.candidateObjective = solveUse.candidateObjective;
  optimInfo.candidateExitflag = solveUse.candidateExitflag;
  optimInfo.candidateVariant = solveUse.candidateVariant;
end
if verbose
  fprintf('[WarmAnchorExit] freezeDoa=%d final variant=%s fval=%.6e resolved=%d\n', ...
    localGetModelLogical(model, 'freezeDoa', false), ...
    string(localGetStructField(solveUse, 'solveVariant', "")), ...
    localGetStructField(solveUse, 'fval', NaN), ...
    localGetStructField(solveUse, 'isResolved', false));
end
end


function modelType = localGetModelType(model)
%LOCALGETMODELTYPE Get the estimator model tag for script-level dispatch.

phaseTagMap = struct();
phaseTagMap.continuous = 'Cp';
phaseTagMap.relaxed = 'Relaxed';
phaseTagMap.independent = 'Ip';

fdRateTagMap = struct();
fdRateTagMap.unknown = 'UnknownRate';
fdRateTagMap.known = 'KnownRate';
fdRateTagMap.zero = 'ZeroRate';

modelType = ['mf', phaseTagMap.(model.phaseMode), fdRateTagMap.(model.fdRateMode)];
end


function optimInfo = localBuildOptimInfo(output, exitflag, runTimeSec, optimOpt, model)
%LOCALBUILDOPTIMINFO Build a compact optimization summary for logging.

optimInfo = struct();
optimInfo.iterations = localGetStructField(output, 'iterations', NaN);
optimInfo.funcCount = localGetStructField(output, 'funcCount', NaN);
optimInfo.constrviolation = localGetStructField(output, 'constrviolation', NaN);
optimInfo.stepsize = localGetStructField(output, 'stepsize', NaN);
optimInfo.firstorderopt = localGetStructField(output, 'firstorderopt', NaN);
optimInfo.algorithm = localGetStructField(output, 'algorithm', '');
optimInfo.message = localGetStructField(output, 'message', '');
optimInfo.exitflag = exitflag;
optimInfo.runTimeSec = runTimeSec;
optimInfo.phaseMode = model.phaseMode;
optimInfo.fdRateMode = model.fdRateMode;
optimInfo.fdRateKnown = model.fdRateKnown;
optimInfo.numVar = localGetNumOptVar(model);
optimInfo.display = localGetStructField(optimOpt, 'Display', '');
optimInfo.modelType = localGetModelType(model);
optimInfo.solver = 'fmincon';
optimInfo.debugEnable = model.debugEnable;
optimInfo.debugStoreEvalTrace = model.debugStoreEvalTrace;
optimInfo.debugMaxEvalTrace = model.debugMaxEvalTrace;
optimInfo.enableFdAliasUnwrap = model.enableFdAliasUnwrap;
optimInfo.fdAliasStepHz = model.fdAliasStepHz;
optimInfo.usedScaledSolve = localGetStructField(output, 'usedScaledSolve', false);
optimInfo.usedFallbackAlgorithm = localGetStructField(output, 'usedFallbackAlgorithm', '');
optimInfo.fallbackTriggered = localGetStructField(output, 'fallbackTriggered', false);
end


function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct/object field with a default value.

fieldValue = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end
if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end


function numVar = localGetNumOptVar(model)
%LOCALGETNUMOPTVAR Get the current optimization vector size.

switch model.fdRateMode
  case 'unknown'
    numVar = 4;
  case {'known', 'zero'}
    numVar = 3;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end
end


function verbose = localResolveVerbose(model, verboseArg)
%LOCALRESOLVEVERBOSE Resolve one consistent branch-level verbose flag.

verbose = verboseArg;
if isstruct(model) && isfield(model, 'verbose') && ~isempty(model.verbose)
  verbose = verbose || logical(model.verbose);
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


function localPrintSolveSummary(tag, solveInfo)
%LOCALPRINTSOLVESUMMARY Print one compact branch solve summary.

if isempty(solveInfo) || ~isstruct(solveInfo)
  fprintf('  %-12s : <empty>\n', tag);
  return;
end
fprintf(['  %-12s : variant=%s resolved=%d exit=%d fval=%.6e ', ...
         'resid=%.6e fit=%.6e support=%.6e iter=%d\n'], ...
  tag, ...
  string(localGetStructField(solveInfo, 'solveVariant', "")), ...
  localGetStructField(solveInfo, 'isResolved', false), ...
  localGetStructField(solveInfo, 'exitflag', NaN), ...
  localGetStructField(solveInfo, 'fval', NaN), ...
  localGetEvalDiagField(localGetStructField(solveInfo, 'finalEvalDiag', struct()), 'residualNorm', NaN), ...
  localGetEvalDiagField(localGetStructField(solveInfo, 'finalEvalDiag', struct()), 'nonRefFitRatioFloor', NaN), ...
  localGetEvalDiagField(localGetStructField(solveInfo, 'finalEvalDiag', struct()), 'nonRefSupportRatioFloor', NaN), ...
  localGetOutputField(localGetStructField(solveInfo, 'output', struct()), 'iterations', -1));
end


function value = localGetEvalDiagField(evalDiag, fieldName, defaultValue)
%LOCALGETEVALDIAGFIELD Read one scalar evaluation diagnostic field.

value = defaultValue;
if ~isstruct(evalDiag) || ~isfield(evalDiag, fieldName)
  return;
end
rawValue = evalDiag.(fieldName);
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


function value = localGetOutputField(output, fieldName, defaultValue)
%LOCALGETOUTPUTFIELD Read one optimization-output scalar field.

value = defaultValue;
if ~isstruct(output) || ~isfield(output, fieldName)
  return;
end
rawValue = output.(fieldName);
if isempty(rawValue) || ~isscalar(rawValue)
  return;
end
value = rawValue;
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

