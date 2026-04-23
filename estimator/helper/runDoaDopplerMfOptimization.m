function solveInfo = runDoaDopplerMfOptimization(model, initParam, lb, ub, Aeq, beq, baseOptimOpt, solveOpt)
%RUNDOADOPPLERMFOPTIMIZATION Run one MF optimization with optional warm-anchor safeguards.
% This helper owns the numerical solve plumbing: scaling, fallback SQP, DoA
% freeze equalities and compact debug traces. Branch selection stays in
% solveDoaDopplerMfBranches, while objective evaluation lives in the
% evaluation helpers.

arguments
  model (1,1) struct
  initParam (:,1) double
  lb (:,1) double
  ub (:,1) double
  Aeq double
  beq double
  baseOptimOpt
  solveOpt (1,1) struct = struct()
end

useScaledSolve = localGetStructField(solveOpt, 'useScaledSolve', false);
enableFallbackSqp = localGetStructField(solveOpt, 'enableFallbackSqp', false);
debugTrace = localInitDebugTrace(model);
[AeqUse, beqUse] = localApplyFreezeDoaConstraint(model, initParam, Aeq, beq);

if useScaledSolve
  solveInfo = localRunOptimizationScaled(model, initParam, lb, ub, AeqUse, beqUse, baseOptimOpt, debugTrace);
else
  solveInfo = localRunOptimization(model, initParam, lb, ub, AeqUse, beqUse, baseOptimOpt, debugTrace);
end

if ~enableFallbackSqp || ~localNeedUnknownWarmFallback(solveInfo)
  return;
end

fallbackOpt = baseOptimOpt;
fallbackOpt.Algorithm = 'sqp';
if useScaledSolve
  solveAlt = localRunOptimizationScaled(model, initParam, lb, ub, AeqUse, beqUse, fallbackOpt, localInitDebugTrace(model));
else
  solveAlt = localRunOptimization(model, initParam, lb, ub, AeqUse, beqUse, fallbackOpt, localInitDebugTrace(model));
end
solveAlt.output.usedFallbackAlgorithm = 'sqp';
solveAlt.output.fallbackTriggered = true;

if preferDoaDopplerMfSolveResult(solveAlt, solveInfo)
  solveInfo = solveAlt;
else
  solveInfo.output.fallbackTriggered = true;
end
end


function solveInfo = localRunOptimization(model, initParam, lb, ub, AeqUse, beqUse, optimOpt, debugTrace)
%LOCALRUNOPTIMIZATION Run one constrained optimization in original coordinates.

runTimer = tic;
objFun = @localObjFunTrace;
[optVar, fval, exitflag, output] = fmincon( ...
  objFun, initParam, [], [], AeqUse, beqUse, lb, ub, [], optimOpt);
runTimeSec = toc(runTimer);
output.usedScaledSolve = false;
output.usedFallbackAlgorithm = '';
output.fallbackTriggered = false;
[~, pathGain, noiseVar, estAux] = evaluateDoaDopplerMfObjective(model, optVar);
finalEvalDiag = buildDoaDopplerMfEvalDiag(model, optVar, fval, noiseVar, estAux);
solveInfo = localBuildSolveInfo(model, optVar, fval, exitflag, output, runTimeSec, ...
  pathGain, noiseVar, estAux, finalEvalDiag, debugTrace);

  function obj = localObjFunTrace(optVarLocal)
    [obj, ~, noiseVarEval, auxEval] = evaluateDoaDopplerMfObjective(model, optVarLocal);
    if model.debugStoreEvalTrace
      debugTrace = localAppendEvalTrace(debugTrace, model, optVarLocal, obj, noiseVarEval, auxEval, 1); %#ok<NASGU>
    end
  end
end


function solveInfo = localRunOptimizationScaled(model, initParam, lb, ub, AeqUse, beqUse, optimOpt, debugTrace)
%LOCALRUNOPTIMIZATIONSCALED Run one bounded solve in normalized coordinates.

[scaleCenter, scaleWidth] = localBuildOptimizationScaling(initParam, lb, ub);
scaledInit = zeros(size(initParam(:)));
scaledLb = (lb(:) - scaleCenter) ./ scaleWidth;
scaledUb = (ub(:) - scaleCenter) ./ scaleWidth;
if isempty(AeqUse)
  AeqScaled = [];
  beqScaled = [];
else
  AeqScaled = AeqUse .* reshape(scaleWidth, 1, []);
  beqScaled = beqUse(:) - AeqUse * scaleCenter;
end

runTimer = tic;
objFun = @localObjFunScaled;
[scaledVar, fval, exitflag, output] = fmincon( ...
  objFun, scaledInit, [], [], AeqScaled, beqScaled, scaledLb, scaledUb, [], optimOpt);
runTimeSec = toc(runTimer);
optVar = scaleCenter + scaleWidth .* scaledVar(:);
[~, pathGain, noiseVar, estAux] = evaluateDoaDopplerMfObjective(model, optVar);
finalEvalDiag = buildDoaDopplerMfEvalDiag(model, optVar, fval, noiseVar, estAux);
output.usedScaledSolve = true;
output.usedFallbackAlgorithm = '';
output.fallbackTriggered = false;
solveInfo = localBuildSolveInfo(model, optVar, fval, exitflag, output, runTimeSec, ...
  pathGain, noiseVar, estAux, finalEvalDiag, debugTrace);

  function obj = localObjFunScaled(optVarScaled)
    optVarUse = scaleCenter + scaleWidth .* optVarScaled(:);
    [obj, ~, noiseVarEval, auxEval] = evaluateDoaDopplerMfObjective(model, optVarUse);
    if model.debugStoreEvalTrace
      debugTrace = localAppendEvalTrace(debugTrace, model, optVarUse, obj, noiseVarEval, auxEval, 1); %#ok<NASGU>
    end
  end
end


function solveInfo = localBuildSolveInfo(model, optVar, fval, exitflag, output, runTimeSec, ...
  pathGain, noiseVar, estAux, finalEvalDiag, debugTrace)
%LOCALBUILDSOLVEINFO Package one solve result.

solveInfo = struct();
solveInfo.optVar = optVar(:);
solveInfo.fval = fval;
solveInfo.exitflag = exitflag;
solveInfo.output = output;
solveInfo.runTimeSec = runTimeSec;
solveInfo.pathGain = pathGain;
solveInfo.noiseVar = noiseVar;
solveInfo.estAux = estAux;
solveInfo.finalEvalDiag = finalEvalDiag;
solveInfo.debugTrace = debugTrace;
solveInfo.isResolved = localCheckResolved(model, optVar, pathGain, noiseVar, fval, exitflag);
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  solveInfo.solveVariant = "mainFixedDoa";
else
  solveInfo.solveVariant = "main";
end
end


function [scaleCenter, scaleWidth] = localBuildOptimizationScaling(initParam, lb, ub)
%LOCALBUILDOPTIMIZATIONSCALING Build one affine variable normalization.

scaleCenter = reshape(initParam, [], 1);
lb = reshape(lb, [], 1);
ub = reshape(ub, [], 1);
if numel(scaleCenter) ~= numel(lb) || numel(lb) ~= numel(ub)
  scaleWidth = ones(size(scaleCenter));
  return;
end

scaleCenter = min(max(scaleCenter, lb), ub);
scaleWidth = 0.5 * abs(ub - lb);
scaleFloor = 1e-6 * max(abs(scaleCenter), 1);
scaleWidth = max(scaleWidth, scaleFloor);
scaleWidth(~isfinite(scaleWidth) | scaleWidth <= 0) = 1;
end


function tf = localNeedUnknownWarmFallback(solveInfo)
%LOCALNEEDUNKNOWNWARMBACK Decide whether the warm solve needs SQP retry.

tf = false;
if ~isstruct(solveInfo) || isempty(solveInfo)
  tf = true;
  return;
end
if ~(isfinite(solveInfo.fval) && all(isfinite(solveInfo.optVar)))
  tf = true;
  return;
end
if solveInfo.isResolved && solveInfo.exitflag > 0
  return;
end
if solveInfo.exitflag <= 0
  tf = true;
  return;
end
firstOrderOpt = localGetStructField(solveInfo.output, 'firstorderopt', NaN);
constrViolation = localGetStructField(solveInfo.output, 'constrviolation', NaN);
if ~isfinite(firstOrderOpt) || ~isfinite(constrViolation)
  tf = true;
end
end


function [AeqUse, beqUse] = localApplyFreezeDoaConstraint(model, initParam, Aeq, beq)
%LOCALAPPLYFREEZEDOACONSTRAINT Convert one frozen DoA seed to equalities.

AeqUse = Aeq;
beqUse = beq;
if ~isfield(model, 'freezeDoa') || ~logical(model.freezeDoa)
  return;
end
freezeCenter = localResolveDoaFreezeCenter(model, initParam);
if isempty(freezeCenter)
  return;
end
numDoa = numel(freezeCenter);
AeqFreeze = zeros(numDoa, numel(initParam));
AeqFreeze(:, 1:numDoa) = eye(numDoa);
beqFreeze = freezeCenter(:);
if isempty(AeqUse)
  AeqUse = AeqFreeze;
  beqUse = beqFreeze;
else
  AeqUse = [AeqUse; AeqFreeze];
  beqUse = [beqUse(:); beqFreeze];
end
end


function freezeCenter = localResolveDoaFreezeCenter(model, initParam)
%LOCALRESOLVEDOAFREEZECENTER Resolve the exact DoA anchor for one solve.

if isfield(model, 'initDoaParam') && ~isempty(model.initDoaParam)
  freezeCenter = reshape(model.initDoaParam, [], 1);
else
  numDoa = min(2, numel(initParam));
  freezeCenter = reshape(initParam(1:numDoa), [], 1);
end
if isempty(freezeCenter)
  return;
end
if isfield(model, 'doaType') && strcmp(model.doaType, 'angle')
  freezeCenter(1) = mod(freezeCenter(1), 2 * pi);
end
end


function debugTrace = localInitDebugTrace(model)
%LOCALINITDEBUGTRACE Initialize the optional objective trace container.

debugTrace = struct();
debugTrace.isEnabled = model.debugStoreEvalTrace;
debugTrace.maxCount = model.debugMaxEvalTrace;
debugTrace.numEval = 0;
debugTrace.numStored = 0;
debugTrace.numDropped = 0;
debugTrace.optVarMat = zeros(localGetNumOptVar(model), 0);
debugTrace.obj = zeros(1, 0);
debugTrace.residualNorm = zeros(1, 0);
debugTrace.fdRef = zeros(1, 0);
debugTrace.fdRate = zeros(1, 0);
debugTrace.angle1 = zeros(1, 0);
debugTrace.angle2 = zeros(1, 0);
debugTrace.startIdx = zeros(1, 0);
end


function debugTrace = localAppendEvalTrace(debugTrace, model, optVar, obj, noiseVar, auxEval, startIdx)
%LOCALAPPENDEVALTRACE Append one compact objective-evaluation record.

if ~debugTrace.isEnabled
  return;
end
debugTrace.numEval = debugTrace.numEval + 1;
if debugTrace.numStored >= debugTrace.maxCount
  debugTrace.numDropped = debugTrace.numDropped + 1;
  return;
end
storeIdx = debugTrace.numStored + 1;
debugTrace.numStored = storeIdx;
debugTrace.optVarMat(:, storeIdx) = optVar(:);
debugTrace.obj(storeIdx) = obj;
debugTrace.residualNorm(storeIdx) = auxEval.residualNorm;
debugTrace.fdRef(storeIdx) = auxEval.fdRef;
debugTrace.fdRate(storeIdx) = auxEval.fdRate;
debugTrace.angle1(storeIdx) = optVar(1);
debugTrace.angle2(storeIdx) = optVar(2);
debugTrace.startIdx(storeIdx) = startIdx;
if nargin >= 5 && ~isempty(noiseVar) %#ok<INUSD>
end
if strcmp(model.doaType, 'angle')
  debugTrace.angle1(storeIdx) = auxEval.eciAngle(1);
  debugTrace.angle2(storeIdx) = auxEval.eciAngle(2);
end
end


function isResolved = localCheckResolved(model, optVar, pathGain, noiseVar, fval, exitflag)
%LOCALCHECKRESOLVED Check whether the optimization result is usable.

validGainMask = model.frameMask;
validNoiseMask = any(model.frameMask, 2);
isResolved = exitflag > 0 && ...
  isfinite(fval) && ...
  all(isfinite(optVar(:))) && ...
  all(isfinite(pathGain(validGainMask))) && ...
  all(isfinite(noiseVar(validNoiseMask)));
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
