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

[lb, ub] = buildDoaDopplerMfBounds(model);
[initObj, ~, initNoiseVar, initEvalAux] = evaluateDoaDopplerMfObjective(model, initParam);
initEvalDiag = buildDoaDopplerMfEvalDiag(model, initParam, initObj, initNoiseVar, initEvalAux);

if verbose
  displayMode = 'iter';
else
  displayMode = 'off';
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

if strcmp(model.fdRateMode, 'known')
  solveEmbed = runDoaDopplerMfKnownEmbeddedOptimization(model, initParam, baseOptimOpt);
  if preferDoaDopplerMfSolveResult(solveEmbed, solveStd)
    solveUse = solveEmbed;
  end
elseif strcmp(model.fdRateMode, 'unknown')
  if isfield(model, 'disableUnknownWarmAnchor') && logical(model.disableUnknownWarmAnchor)
    solveCand = struct([]);
  else
    solveCand = runDoaDopplerMfUnknownWarmAnchor(model, initParam, lb, ub, baseOptimOpt);
  end
  if ~isempty(solveCand)
    candidateObjective = solveCand.fval;
    candidateExitflag = solveCand.exitflag;
    candidateVariant = string(solveCand.solveVariant);
    if useDoaDopplerMfUnknownWarmAnchorResult(model, solveCand, solveUse)
      solveUse = solveCand;
    end
    solveUse.candidateObjective = candidateObjective;
    solveUse.candidateExitflag = candidateExitflag;
    solveUse.candidateVariant = candidateVariant;
  end
end

optimInfo = localBuildOptimInfo(solveUse.output, solveUse.exitflag, ...
  solveUse.runTimeSec, baseOptimOpt, model);
optimInfo.solveVariant = solveUse.solveVariant;
if isfield(solveUse, 'candidateObjective')
  optimInfo.candidateObjective = solveUse.candidateObjective;
  optimInfo.candidateExitflag = solveUse.candidateExitflag;
  optimInfo.candidateVariant = solveUse.candidateVariant;
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
