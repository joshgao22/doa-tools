function solveInfo = runDoaDopplerMfKnownEmbeddedOptimization(model, initParamKnown, baseOptimOpt)
%RUNDOADOPPLERMFKNOWNEMBEDDEDOPTIMIZATION Solve known-rate mode on unknown layout.

arguments
  model (1,1) struct
  initParamKnown (:,1) double
  baseOptimOpt
end

solveModel = localBuildKnownEmbeddedModel(model);
[lbSolve, ubSolve] = buildDoaDopplerMfBounds(solveModel);
initParamSolve = [initParamKnown(:); model.fdRateKnown];
initParamSolve = min(max(initParamSolve, lbSolve), ubSolve);
Aeq = zeros(1, numel(initParamSolve));
Aeq(end) = 1;
beq = model.fdRateKnown;
solveOpt = struct('useScaledSolve', false, 'enableFallbackSqp', false);
solveInfoRaw = runDoaDopplerMfOptimization(solveModel, initParamSolve, lbSolve, ubSolve, Aeq, beq, baseOptimOpt, solveOpt);

optVarKnown = solveInfoRaw.optVar(1:3);
[fvalKnown, pathGainKnown, noiseVarKnown, estAuxKnown] = evaluateDoaDopplerMfObjective(model, optVarKnown);
finalEvalDiagKnown = buildDoaDopplerMfEvalDiag(model, optVarKnown, fvalKnown, noiseVarKnown, estAuxKnown);

solveInfo = solveInfoRaw;
solveInfo.optVar = optVarKnown(:);
solveInfo.fval = fvalKnown;
solveInfo.noiseVar = noiseVarKnown;
solveInfo.pathGain = pathGainKnown;
solveInfo.estAux = estAuxKnown;
solveInfo.finalEvalDiag = finalEvalDiagKnown;
solveInfo.isResolved = localCheckResolved(model, solveInfo.optVar, ...
  pathGainKnown, noiseVarKnown, fvalKnown, solveInfoRaw.exitflag);
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  solveInfo.solveVariant = "knownEmbeddedUnknownLayoutFixedDoa";
else
  solveInfo.solveVariant = "knownEmbeddedUnknownLayout";
end
end


function solveModel = localBuildKnownEmbeddedModel(model)
%LOCALBUILDKNOWNEMBEDDEDMODEL Clone one known-rate model onto unknown layout.

solveModel = model;
solveModel.fdRateMode = 'unknown';
solveModel.fdRateKnown = [];
solveModel.cachedBoundsReady = false;
solveModel.lb = [];
solveModel.ub = [];
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
