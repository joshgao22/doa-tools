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

[lb, ub] = localBuildBounds(model);

[initObj, ~, initNoiseVar, initEvalAux] = localNllPilotDynCon(model, initParam);
initEvalDiag = localBuildEvalDiag(model, initParam, initObj, initNoiseVar, initEvalAux);

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

solveStd = localRunOptimization(model, initParam, lb, ub, [], [], baseOptimOpt);
solveUse = solveStd;

if strcmp(model.fdRateMode, 'known')
  solveEmbed = localRunKnownEmbeddedOptimization(model, initParam, baseOptimOpt);
  if localPreferSolveResult(solveEmbed, solveStd)
    solveUse = solveEmbed;
  end
elseif strcmp(model.fdRateMode, 'unknown')
  if isfield(model, 'disableUnknownWarmAnchor') && logical(model.disableUnknownWarmAnchor)
    solveCand = struct([]);
  else
    solveCand = localRunUnknownMainFromWarmAnchor(model, initParam, lb, ub, baseOptimOpt);
  end
  if ~isempty(solveCand)
    candidateObjective = solveCand.fval;
    candidateExitflag = solveCand.exitflag;
    candidateVariant = string(solveCand.solveVariant);
    if localPreferSolveResult(solveCand, solveUse)
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

function solveInfo = localRunOptimization(model, initParam, lb, ub, Aeq, beq, baseOptimOpt)
%LOCALRUNOPTIMIZATION Run one constrained optimization on the given model.

if nargin < 5 || isempty(Aeq)
  Aeq = [];
end
if nargin < 6 || isempty(beq)
  beq = [];
end

debugTrace = localInitDebugTrace(model);
runTimer = tic;
objFun = @localObjFunTrace;
[optVar, fval, exitflag, output] = fmincon( ...
  objFun, initParam, [], [], Aeq, beq, lb, ub, [], baseOptimOpt);
runTimeSec = toc(runTimer);
[~, pathGain, noiseVar, estAux] = localNllPilotDynCon(model, optVar);
finalEvalDiag = localBuildEvalDiag(model, optVar, fval, noiseVar, estAux);

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
solveInfo.solveVariant = "main";

  function obj = localObjFunTrace(optVarLocal)
  %LOCALOBJFUNTRACE Wrapped objective with optional debug trace storage.

  [obj, ~, noiseVarEval, auxEval] = localNllPilotDynCon(model, optVarLocal);
  if model.debugStoreEvalTrace
    debugTrace = localAppendEvalTrace(debugTrace, model, optVarLocal, ...
      obj, noiseVarEval, auxEval, 1);
  end
  end
end

function solveCand = localRunUnknownMainFromWarmAnchor(model, initParam, lb, ub, baseOptimOpt)
%LOCALRUNUNKNOWNMAINFROMWARMANCHOR Stabilize CP-U with one warm-rate anchor.
% Keep the active CP-U logic compact: first hold fdRate at the current seed
% long enough to pull DoA into the better continuous-phase basin, then
% release the full unknown branch again while recentering the internal DoA
% floor around that warm anchor. To avoid getting numerically stuck at the
% fixed-rate warm point, probe a compact local set of fdRate release seeds
% around the current anchor, but still report the whole path as the single
% mainWarmAnchor branch.

solveCand = struct([]);
if ~strcmp(model.fdRateMode, 'unknown')
  return;
end
if model.numSat <= 1 || ~strcmp(model.phaseMode, 'continuous')
  return;
end
if numel(initParam) < 4
  return;
end

fdRateSeed = initParam(end);
Aeq = zeros(1, numel(initParam));
Aeq(end) = 1;
beq = fdRateSeed;
solveWarm = localRunOptimization(model, initParam, lb, ub, Aeq, beq, baseOptimOpt);

warmInit = solveWarm.optVar(:);
[lbAnchor, ubAnchor] = localBuildBounds(model);
numDoa = numel(model.initDoaParam);
if numDoa >= 1 && numel(lbAnchor) >= numDoa && numel(ubAnchor) >= numDoa
  [doaLbAnchor, doaUbAnchor] = localApplyUnknownDoaReleaseFloor(model, ...
    lbAnchor(1:numDoa), ubAnchor(1:numDoa), warmInit(1:numDoa));
  lbAnchor(1:numDoa) = doaLbAnchor;
  ubAnchor(1:numDoa) = doaUbAnchor;
end

fdRateReleaseSeed = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed);
if isempty(fdRateReleaseSeed)
  fdRateReleaseSeed = fdRateSeed;
end
fdRateReleaseSeed = reshape(fdRateReleaseSeed, 1, []);

solveBest = struct([]);
totalRunTimeSec = solveWarm.runTimeSec;
for iSeed = 1:numel(fdRateReleaseSeed)
  initRelease = warmInit;
  initRelease(end) = fdRateReleaseSeed(iSeed);

  [lbLocal, ubLocal, releaseHalfWidth] = localBuildUnknownFdRateReleaseBounds( ...
    model, lbAnchor, ubAnchor, fdRateReleaseSeed(iSeed), fdRateReleaseSeed);
  solveLocal = localRunOptimization(model, initRelease, lbLocal, ubLocal, [], [], baseOptimOpt);
  solveFull = localRunOptimization(model, solveLocal.optVar(:), lbAnchor, ubAnchor, [], [], baseOptimOpt);

  totalRunTimeSec = totalRunTimeSec + solveLocal.runTimeSec + solveFull.runTimeSec;
  solveFull.output = localMergeOptimOutput(solveLocal.output, solveFull.output);
  solveFull.warmStart = struct( ...
    'fdRateFixed', fdRateSeed, ...
    'warmFval', solveWarm.fval, ...
    'warmExitflag', solveWarm.exitflag, ...
    'warmIsResolved', solveWarm.isResolved, ...
    'warmVariant', "mainWarmFixedFdRate", ...
    'releaseSeedFdRate', fdRateReleaseSeed(iSeed), ...
    'releaseHalfWidth', releaseHalfWidth, ...
    'releaseLb', lbLocal(end), ...
    'releaseUb', ubLocal(end));

  if isempty(solveBest) || localPreferSolveResult(solveFull, solveBest)
    solveBest = solveFull;
  end
end

if isempty(solveBest)
  return;
end

solveBest.runTimeSec = totalRunTimeSec;
solveBest.output = localMergeOptimOutput(solveWarm.output, solveBest.output);
solveBest.solveVariant = "mainWarmAnchor";
solveBest.warmStart.fdRateReleaseSeedList = fdRateReleaseSeed(:);
solveCand = solveBest;
end

function solveCand = localRunUnknownLegacyWarmAnchors(model, initParam, lb, ub, baseOptimOpt)
%LOCALRUNUNKNOWNCONTINUATIONCANDIDATES Build robust CP-U continuation starts.
% The unknown-rate branch is more sensitive to noise because fdRate couples
% with the continuous-phase trajectory. For the multi-satellite continuous
% mainline, first solve several equality-constrained fixed-rate subproblems,
% then release fdRate only inside a local box around each fixed-rate anchor.
% This keeps CP-U close to the same basin as the corresponding CP-K anchor
% instead of inviting the optimizer to jump to the zero-rate branch.

solveCand = struct([]);
if ~strcmp(model.fdRateMode, 'unknown')
  return;
end
if model.numSat <= 1 || ~strcmp(model.phaseMode, 'continuous')
  return;
end
if numel(initParam) < 4
  return;
end

fdRateSeed = initParam(end);
fdRateCand = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed);
if isempty(fdRateCand)
  return;
end

for iCand = 1:numel(fdRateCand)
  rateVal = fdRateCand(iCand);
  initWarm = initParam(:);
  initWarm(end) = rateVal;

  Aeq = zeros(1, numel(initWarm));
  Aeq(end) = 1;
  beq = rateVal;
  solveWarm = localRunOptimization(model, initWarm, lb, ub, Aeq, beq, baseOptimOpt);

  warmInit = solveWarm.optVar(:);
  warmInit(end) = rateVal;
  [lbLocal, ubLocal, releaseHalfWidth] = localBuildUnknownFdRateReleaseBounds(model, lb, ub, rateVal, fdRateCand);
  solveRelease = localRunOptimization(model, warmInit, lbLocal, ubLocal, [], [], baseOptimOpt);
  solveRelease.runTimeSec = solveWarm.runTimeSec + solveRelease.runTimeSec;
  solveRelease.output = localMergeOptimOutput(solveWarm.output, solveRelease.output);
  solveRelease.solveVariant = sprintf('legacyWarmAnchor_fdRate%.6g', rateVal);
  solveRelease.warmStart = struct( ...
    'fdRateFixed', rateVal, ...
    'warmFval', solveWarm.fval, ...
    'warmExitflag', solveWarm.exitflag, ...
    'warmIsResolved', solveWarm.isResolved, ...
    'warmVariant', sprintf('legacyWarmFixed_fdRate%.6g', rateVal), ...
    'releaseHalfWidth', releaseHalfWidth, ...
    'releaseLb', lbLocal(end), ...
    'releaseUb', ubLocal(end));
  if isempty(solveCand)
    solveCand = solveRelease;
  else
    solveCand(end + 1) = solveRelease;
  end
end
end

function [lbLocal, ubLocal, releaseHalfWidth] = localBuildUnknownFdRateReleaseBounds(model, lb, ub, rateVal, fdRateCand)
%LOCALBUILDUNKNOWNFDRATERELEASEBOUNDS Tighten the unknown-rate release box.

lbLocal = lb(:);
ubLocal = ub(:);
releaseHalfWidth = 0;
if numel(lbLocal) ~= numel(ubLocal) || isempty(lbLocal)
  return;
end

fdRateIdx = numel(lbLocal);
rangeSpan = ubLocal(fdRateIdx) - lbLocal(fdRateIdx);
if ~(isfinite(rangeSpan) && rangeSpan > 0)
  return;
end

if nargin >= 5 && ~isempty(fdRateCand)
  fdRateCand = sort(fdRateCand(:));
  if numel(fdRateCand) >= 2
    gridStep = median(diff(fdRateCand));
  else
    gridStep = rangeSpan;
  end
else
  gridStep = rangeSpan;
end

releaseHalfWidth = min(0.1 * rangeSpan, max(0.5 * abs(gridStep), 1e-6 * max(abs(rateVal), 1)));
if ~isfinite(releaseHalfWidth) || releaseHalfWidth <= 0
  releaseHalfWidth = 0.25 * rangeSpan;
end

lbLocal(fdRateIdx) = max(lbLocal(fdRateIdx), rateVal - releaseHalfWidth);
ubLocal(fdRateIdx) = min(ubLocal(fdRateIdx), rateVal + releaseHalfWidth);
if lbLocal(fdRateIdx) > ubLocal(fdRateIdx)
  lbLocal(fdRateIdx) = lb(fdRateIdx);
  ubLocal(fdRateIdx) = ub(fdRateIdx);
  releaseHalfWidth = 0.5 * (ubLocal(fdRateIdx) - lbLocal(fdRateIdx));
end
end

function outputMerge = localMergeOptimOutput(outputWarm, outputRelease)
%LOCALMERGEOPTIMOUTPUT Merge warm-start and release optimization summaries.

outputMerge = outputRelease;
outputMerge.iterations = localGetStructField(outputWarm, 'iterations', 0) + ...
  localGetStructField(outputRelease, 'iterations', 0);
outputMerge.funcCount = localGetStructField(outputWarm, 'funcCount', 0) + ...
  localGetStructField(outputRelease, 'funcCount', 0);
outputMerge.constrviolation = localGetStructField(outputRelease, 'constrviolation', NaN);
outputMerge.stepsize = localGetStructField(outputRelease, 'stepsize', NaN);
outputMerge.firstorderopt = localGetStructField(outputRelease, 'firstorderopt', NaN);
outputMerge.algorithm = localGetStructField(outputRelease, 'algorithm', '');
outputMerge.message = localGetStructField(outputRelease, 'message', '');
end

function fdRateCand = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed)
%LOCALBUILDUNKNOWNFDRATEWARMSTARTSET Build local fdRate anchors for CP-U.
% CP-U should stay near the current dynamic seed instead of adding a coarse
% interior anchor such as -750 Hz/s when the seed lies at a range boundary.
% Such coarse anchors repeatedly trapped the optimizer in a stable but worse
% multi-satellite basin. Here we keep only a compact local set around the
% current seed and rely on the refined reference-satellite fdRate init to
% place that seed near the correct basin.

fdRateCand = [];
if isempty(model.fdRateRange) || numel(model.fdRateRange) ~= 2
  return;
end

fdMin = model.fdRateRange(1);
fdMax = model.fdRateRange(2);
if ~(isfinite(fdMin) && isfinite(fdMax) && fdMin <= fdMax)
  return;
end

rangeSpan = fdMax - fdMin;
if rangeSpan <= 0
  fdRateCand = fdMin;
  return;
end
if ~isfinite(fdRateSeed)
  fdRateSeed = 0.5 * (fdMin + fdMax);
end

localHalfWidth = min(0.08 * rangeSpan, max(25, 0.03 * max(abs(fdRateSeed), rangeSpan)));
candidateVec = [fdRateSeed - localHalfWidth, fdRateSeed, fdRateSeed + localHalfWidth];
fdRateCand = candidateVec(isfinite(candidateVec));
if isempty(fdRateCand)
  return;
end
fdRateCand = min(max(fdRateCand, fdMin), fdMax);
fdRateCand = localUniqueStableTol(fdRateCand(:).', 1e-9);
end

function uniqueVec = localUniqueStableTol(valVec, tol)
%LOCALUNIQUESTABLETOL Stable unique with tolerance for warm-start seeds.

if nargin < 2 || isempty(tol)
  tol = 0;
end
valVec = reshape(valVec, 1, []);
if isempty(valVec)
  uniqueVec = valVec;
  return;
end

keepMask = false(size(valVec));
uniqueList = zeros(1, 0);
for iVal = 1:numel(valVec)
  if isempty(uniqueList) || all(abs(valVec(iVal) - uniqueList) > tol)
    uniqueList(end + 1) = valVec(iVal); %#ok<AGROW>
    keepMask(iVal) = true;
  end
end
uniqueVec = valVec(keepMask);
end

function solveInfo = localRunKnownEmbeddedOptimization(model, initParamKnown, baseOptimOpt)
%LOCALRUNKNOWNEMBEDDEDOPTIMIZATION Solve known-rate mode on the unknown layout.
% The CP-known branch should follow the same parameter layout and evaluator
% path as CP-unknown, with fdRate fixed to the known value, so that the
% only model difference is whether fdRate is free or constrained.

solveModel = localBuildKnownEmbeddedModel(model);
[lbSolve, ubSolve] = localBuildBounds(solveModel);
initParamSolve = [initParamKnown(:); model.fdRateKnown];
initParamSolve = min(max(initParamSolve, lbSolve), ubSolve);
Aeq = zeros(1, numel(initParamSolve));
Aeq(end) = 1;
beq = model.fdRateKnown;

solveInfoRaw = localRunOptimization(solveModel, initParamSolve, lbSolve, ubSolve, Aeq, beq, baseOptimOpt);
optVarKnown = solveInfoRaw.optVar(1:3);
[fvalKnown, ~, noiseVarKnown, estAuxKnown] = localNllPilotDynCon(model, optVarKnown);
finalEvalDiagKnown = localBuildEvalDiag(model, optVarKnown, fvalKnown, noiseVarKnown, estAuxKnown);
pathGainKnown = solveInfoRaw.pathGain;

solveInfo = solveInfoRaw;
solveInfo.optVar = optVarKnown(:);
solveInfo.fval = fvalKnown;
solveInfo.noiseVar = noiseVarKnown;
solveInfo.pathGain = pathGainKnown;
solveInfo.estAux = estAuxKnown;
solveInfo.finalEvalDiag = finalEvalDiagKnown;
solveInfo.isResolved = localCheckResolved(model, solveInfo.optVar, ...
  pathGainKnown, noiseVarKnown, fvalKnown, solveInfoRaw.exitflag);
solveInfo.solveVariant = "knownEmbeddedUnknownLayout";
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

function useAlt = localPreferSolveResult(solveAlt, solveBase)
%LOCALPREFERSOLVERESULT Select the better solve result conservatively.

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

function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build the DoA search box used by all MF refinements.
% The local DoA refinement must be enforced by the optimizer bounds rather
% than only stored as a nominal debug setting. Precomputing the DoA box here
% makes every path reuse the same local basin definition.

if isfield(model, 'cachedBoundsReady') && logical(model.cachedBoundsReady) && ...
    isfield(model, 'doaLb') && isfield(model, 'doaUb') && ...
    ~isempty(model.doaLb) && ~isempty(model.doaUb)
  doaLb = model.doaLb;
  doaUb = model.doaUb;
  return;
end

baseRange = model.doaGrid{1}.range;
doaLb = baseRange(:, 1);
doaUb = baseRange(:, 2);

if isempty(model.initDoaParam) || isempty(model.initDoaHalfWidth)
  return;
end

doaCenter = model.initDoaParam(:);
doaHalfWidth = model.initDoaHalfWidth(:);

doaLbLocal = doaCenter - doaHalfWidth;
doaUbLocal = doaCenter + doaHalfWidth;

if strcmp(model.doaType, 'angle')
  doaCenter(1) = mod(doaCenter(1), 2 * pi);
  doaLbLocal(1) = doaCenter(1) - doaHalfWidth(1);
  doaUbLocal(1) = doaCenter(1) + doaHalfWidth(1);
  if doaLbLocal(1) < baseRange(1, 1) || doaUbLocal(1) > baseRange(1, 2)
    doaLbLocal(1) = baseRange(1, 1);
    doaUbLocal(1) = baseRange(1, 2);
  end
end

doaLb = max(doaLb, doaLbLocal);
doaUb = min(doaUb, doaUbLocal);
invalidMask = doaLb > doaUb;
doaLb(invalidMask) = baseRange(invalidMask, 1);
doaUb(invalidMask) = baseRange(invalidMask, 2);
end

function [lb, ub] = localBuildBounds(model)
%LOCALBUILDBOUNDS Build box constraints for fmincon.

if isfield(model, 'lb') && isfield(model, 'ub') && ...
    ~isempty(model.lb) && ~isempty(model.ub)
  lb = model.lb;
  ub = model.ub;
  return;
end

[doaLb, doaUb] = localBuildDoaBounds(model);
[doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb);

lb = [doaLb; model.fdRange(1)];
ub = [doaUb; model.fdRange(2)];

if strcmp(model.fdRateMode, 'unknown')
  lb = [lb; model.fdRateRange(1)];
  ub = [ub; model.fdRateRange(2)];
end
if isfield(model, 'cachedBoundsReady') && ~logical(model.cachedBoundsReady)
  model.lb = lb;
  model.ub = ub;
end
end

function [doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb, centerOverride)
%LOCALAPPLYUNKNOWNDOARELEASEFLOOR Enforce a minimum CP-U DoA release box.
% The default MF unknown branch previously inherited the very tight local
% DoA box from the CP-K/static seed. In the multi-satellite continuous
% branch this often froze DoA, so CP-U only released fdRate and could not
% continue descending along the joint DoA-Doppler basin. Here we keep the
% user-provided local box when it is already wide enough, but apply a small
% minimum half-width floor for the default CP-U solve. When a warm anchor is
% provided, the floor is also allowed to recenter around that anchor.

if nargin < 4
  centerOverride = [];
end
if ~strcmp(model.fdRateMode, 'unknown')
  return;
end
if model.numSat <= 1 || ~strcmp(model.phaseMode, 'continuous')
  return;
end
if isempty(model.initDoaParam) && isempty(centerOverride)
  return;
end

baseRange = model.doaGrid{1}.range;
if isempty(baseRange) || any(size(baseRange) ~= [2, 2])
  return;
end

releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model);
if isempty(releaseHalfWidth)
  return;
end

currentWidth = doaUb - doaLb;
targetWidth = 2 * releaseHalfWidth;
if ~isempty(centerOverride)
  center = reshape(centerOverride, [], 1);
else
  center = reshape(model.initDoaParam, [], 1);
end

needRecenter = false;
if ~isempty(centerOverride)
  centerTol = 1e-12;
  needRecenter = any(center < doaLb + releaseHalfWidth - centerTol) || ...
    any(center > doaUb - releaseHalfWidth + centerTol);
end
if all(currentWidth >= targetWidth - 1e-12) && ~needRecenter
  return;
end

doaLbFloor = center - releaseHalfWidth;
doaUbFloor = center + releaseHalfWidth;

if strcmp(model.doaType, 'angle')
  center(1) = mod(center(1), 2 * pi);
  doaLbFloor(1) = center(1) - releaseHalfWidth(1);
  doaUbFloor(1) = center(1) + releaseHalfWidth(1);
  if doaLbFloor(1) < baseRange(1, 1) || doaUbFloor(1) > baseRange(1, 2)
    doaLbFloor(1) = baseRange(1, 1);
    doaUbFloor(1) = baseRange(1, 2);
  end
end

for iDim = 1:numel(doaLb)
  if needRecenter || currentWidth(iDim) < targetWidth(iDim)
    doaLb(iDim) = max(baseRange(iDim, 1), doaLbFloor(iDim));
    doaUb(iDim) = min(baseRange(iDim, 2), doaUbFloor(iDim));
  end
end

invalidMask = doaLb > doaUb;
doAResetLb = baseRange(invalidMask, 1);
doAResetUb = baseRange(invalidMask, 2);
doAResetLb = reshape(doAResetLb, [], 1);
doAResetUb = reshape(doAResetUb, [], 1);
doAIdx = find(invalidMask);
for iIdx = 1:numel(doAIdx)
  doaLb(doAIdx(iIdx)) = doAResetLb(iIdx);
  doaUb(doAIdx(iIdx)) = doAResetUb(iIdx);
end
end

function releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model)
%LOCALRESOLVEUNKNOWNDOARELEASEHALFWIDTH Resolve the CP-U DoA box floor.

if isfield(model, 'unknownDoaReleaseHalfWidth') && ...
    ~isempty(model.unknownDoaReleaseHalfWidth)
  releaseHalfWidth = reshape(model.unknownDoaReleaseHalfWidth, [], 1);
else
  if strcmp(model.doaType, 'angle')
    defaultHalfWidth = deg2rad(0.03);
  else
    defaultHalfWidth = 0.03;
  end
  releaseHalfWidth = repmat(defaultHalfWidth, 2, 1);
end

if numel(releaseHalfWidth) == 1
  releaseHalfWidth = repmat(releaseHalfWidth, 2, 1);
end
releaseHalfWidth = reshape(releaseHalfWidth, [], 1);
if numel(releaseHalfWidth) ~= 2 || any(~isfinite(releaseHalfWidth)) || ...
    any(releaseHalfWidth < 0)
  releaseHalfWidth = [];
end
end

function [fdRefPerSat, fdRatePerSat] = localFitSatFdLine(timeOffsetSec, fdSeq)
%LOCALFITSATFDLINE Fit one affine Doppler trend per satellite.

numSat = size(fdSeq, 1);
fdRefPerSat = zeros(numSat, 1);
fdRatePerSat = zeros(numSat, 1);

timeVec = reshape(timeOffsetSec, [], 1);
if numel(timeVec) < 2 || max(timeVec) - min(timeVec) <= 0
  fdRefPerSat = fdSeq(:, 1);
  return;
end

designMat = [ones(numel(timeVec), 1), timeVec];

for iSat = 1:numSat
  fdVec = reshape(fdSeq(iSat, :), [], 1);
  validMask = isfinite(fdVec) & isfinite(timeVec);

  if nnz(validMask) == 0
    continue;
  elseif nnz(validMask) == 1
    fdRefPerSat(iSat) = fdVec(validMask);
    continue;
  end

  coef = designMat(validMask, :) \ fdVec(validMask);
  fdRefPerSat(iSat) = coef(1);
  fdRatePerSat(iSat) = coef(2);
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
  % Keep the signature aligned with localNllPilotDynCon outputs.
end
if strcmp(model.doaType, 'angle')
  debugTrace.angle1(storeIdx) = auxEval.eciAngle(1);
  debugTrace.angle2(storeIdx) = auxEval.eciAngle(2);
end
end

function evalDiag = localBuildEvalDiag(model, optVar, obj, noiseVar, evalAux)
%LOCALBUILDEVALDIAG Pack one compact evaluation summary for debugging.

evalDiag = struct();
evalDiag.optVar = optVar(:);
evalDiag.obj = obj;
evalDiag.doaType = model.doaType;
evalDiag.phaseMode = model.phaseMode;
evalDiag.fdRateMode = model.fdRateMode;
evalDiag.refSatIdxLocal = localGetStructField(evalAux, 'refSatIdxLocal', NaN);
evalDiag.doaParam = evalAux.doaParam;
evalDiag.eciAngle = evalAux.eciAngle;
evalDiag.latlon = evalAux.latlon;
evalDiag.fdRef = evalAux.fdRef;
evalDiag.fdRate = evalAux.fdRate;
evalDiag.deltaFdRef = evalAux.deltaFdRef(:);
evalDiag.deltaFdRefRaw = localGetStructField(evalAux, 'deltaFdRefRaw', evalAux.deltaFdRef(:));
evalDiag.deltaFdRefEval = localGetStructField(evalAux, 'deltaFdRefEval', evalAux.deltaFdRef(:));
evalDiag.deltaFdRate = evalAux.deltaFdRate(:);
evalDiag.deltaFdRateRaw = localGetStructField(evalAux, 'deltaFdRateRaw', evalAux.deltaFdRate(:));
evalDiag.deltaFdRateEval = localGetStructField(evalAux, 'deltaFdRateEval', evalAux.deltaFdRate(:));
evalDiag.fdSat = evalAux.fdSat;
evalDiag.fdSatRaw = localGetStructField(evalAux, 'fdSatRaw', evalAux.fdSat);
evalDiag.fdSatEval = evalAux.fdSatEval;
evalDiag.fdRateSat = evalAux.fdRateSat;
evalDiag.fdLocal = evalAux.fdLocal;
evalDiag.fdLocalRaw = localGetStructField(evalAux, 'fdLocalRaw', evalAux.fdLocal);
evalDiag.fdLocalEval = evalAux.fdLocalEval;
evalDiag.localDoaArrUsed = localGetStructField(evalAux, 'localDoaArrUsed', []);
evalDiag.fdAliasStepHz = evalAux.fdAliasStepHz;
evalDiag.fdAliasIndex = evalAux.fdAliasIndex;
evalDiag.fdAliasShiftHz = evalAux.fdAliasShiftHz;
evalDiag.noiseVar = noiseVar(:);
evalDiag.countPerSat = evalAux.countPerSat(:);
evalDiag.residualSat = localGetStructField(evalAux, 'residualSat', noiseVar(:) .* evalAux.countPerSat(:));
evalDiag.fitValueSat = localGetStructField(evalAux, 'fitValueSat', nan(size(evalDiag.residualSat)));
evalDiag.objectiveSat = localGetStructField(evalAux, 'objectiveSat', nan(size(evalDiag.residualSat)));
evalDiag.residualNorm = evalAux.residualNorm;
evalDiag.noiseVarGlobal = evalAux.noiseVarGlobal;
evalDiag.coherenceSat = localCalcSatCoherence(evalAux.zMat);
evalDiag.zMat = localGetStructField(evalAux, 'zMat', []);
evalDiag.etaMat = localGetStructField(evalAux, 'etaMat', []);
evalDiag.blockValueMat = localGetStructField(evalAux, 'blockValue', []);
evalDiag.blockNorm2Mat = localGetStructField(evalAux, 'blockNorm2', []);
evalDiag.blockResidualMat = evalDiag.blockNorm2Mat - evalDiag.blockValueMat;
evalDiag.countPerBlock = localGetStructField(evalAux, 'countPerBlock', []);
evalDiag.blockAbsZMat = abs(evalDiag.zMat);
evalDiag.blockPhaseResidMat = localCalcBlockPhaseResid(evalDiag.zMat, evalAux.framePhase);
end

function coherenceSat = localCalcSatCoherence(zMat)
%LOCALCALCSATCOHERENCE Compute |sum(z)| / sum(|z|) per satellite.

if isempty(zMat)
  coherenceSat = [];
  return;
end

magSum = sum(abs(zMat), 2);
coherenceSat = nan(size(magSum));
validMask = magSum > 0;
coherenceSat(validMask) = abs(sum(zMat(validMask, :), 2)) ./ magSum(validMask);
coherenceSat = coherenceSat(:);
end

function blockPhaseResidMat = localCalcBlockPhaseResid(zMat, framePhaseMat)
%LOCALCALCBLOCKPHASERESID Compute angle(z)-phase for each valid block.

blockPhaseResidMat = nan(size(zMat));
if isempty(zMat) || isempty(framePhaseMat)
  return;
end

validMask = isfinite(zMat) & isfinite(framePhaseMat);
blockPhaseResidMat(validMask) = angle(zMat(validMask)) - framePhaseMat(validMask);
blockPhaseResidMat(validMask) = mod(blockPhaseResidMat(validMask) + pi, 2 * pi) - pi;
end

function [obj, pathGain, noiseVar, aux] = localNllPilotDynCon(model, optVar)
%LOCALNLLPILOTDYNCON Concentrated dynamic pilot negative log-likelihood.

hyp = localBuildHypothesis(model, optVar);
[obj, prof, profAux] = evalDoaDopplerDynProfileLike(model, hyp);

if nargout <= 1
  return;
end

pathGain = prof.pathGainEst;
noiseVar = prof.noiseVarEst;

aux = struct();
aux.doaParam = hyp.doaParam;
aux.eciAngle = hyp.eciAngle;
aux.eciDirection = hyp.eciDirection;
aux.latlon = hyp.latlon;
aux.userPosEci = hyp.userPosEci;
aux.userVelEci = hyp.userVelEci;

aux.refSatIdxLocal = localGetStructField(profAux, 'refSatIdxLocal', NaN);
aux.deltaFdSeq = hyp.deltaFdSeq;
aux.deltaFdRef = hyp.deltaFdRef;
aux.deltaFdRefRaw = localGetStructField(profAux, 'deltaFdRefRaw', hyp.deltaFdRef);
aux.deltaFdRefEval = localGetStructField(profAux, 'deltaFdRefEval', hyp.deltaFdRef);
aux.deltaFdRate = hyp.deltaFdRate;
aux.deltaFdRateRaw = localGetStructField(profAux, 'deltaFdRateRaw', hyp.deltaFdRate);
aux.deltaFdRateEval = localGetStructField(profAux, 'deltaFdRateEval', hyp.deltaFdRate);
aux.fdRef = hyp.fdRef;
aux.fdRate = hyp.fdRate;
aux.fdSat = hyp.fdSat;
aux.fdSatRaw = localGetStructField(profAux, 'fdSatRaw', hyp.fdSat);
aux.fdRateSat = hyp.fdRateSat;
aux.fdLocal = profAux.fdLocal;
aux.fdLocalRaw = localGetStructField(profAux, 'fdLocalRaw', profAux.fdLocal);
aux.fdSatEval = profAux.fdSatEval;
aux.fdLocalEval = profAux.fdLocalEval;
aux.fdAliasStepHz = profAux.fdAliasStepHz;
aux.fdAliasIndex = profAux.fdAliasIndex;
aux.fdAliasShiftHz = profAux.fdAliasShiftHz;

aux.localDoaArrUsed = profAux.localDoaArrUsed;
aux.localDoaArr = profAux.localDoaArrUsed;
aux.phaseSat = prof.phaseSatEst;
aux.framePhase = prof.framePhaseEst;
aux.ampEst = prof.ampEst;
aux.noiseVarGlobal = prof.noiseVarGlobal;
aux.residualNorm = prof.residualNorm;
aux.countPerSat = prof.countPerSat;
aux.residualSat = localGetStructField(prof, 'residualSat', []);
aux.fitValueSat = localGetStructField(prof, 'fitValueSat', []);
aux.objectiveSat = localGetStructField(prof, 'objectiveSat', []);
aux.blockValue = prof.blockValue;
aux.blockNorm2 = prof.blockNorm2;
aux.countPerBlock = localGetStructField(prof, 'countPerBlock', []);
aux.zMat = prof.zMat;
aux.etaMat = prof.etaMat;
end

function hyp = localBuildHypothesis(model, optVar)
%LOCALBUILDHYPOTHESIS Map one optimization point to geometry and Doppler.

optVar = optVar(:);
numVar = localGetNumOptVar(model);
if numel(optVar) ~= numVar
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidOptVarSize', ...
    'The optimization vector must contain %d entries for fdRateMode = ''%s''.', ...
    numVar, model.fdRateMode);
end

doaParam = optVar(1:2);
fdRef = optVar(3);

switch model.fdRateMode
  case 'unknown'
    fdRate = optVar(4);
  case 'known'
    fdRate = model.fdRateKnown;
  case 'zero'
    fdRate = 0;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end

doaState = localBuildDoaState(model, doaParam);
fdSat = doaState.deltaFdRef + fdRef;
fdRateSat = doaState.deltaFdRate + fdRate;

hyp = doaState;
hyp.doaParam = doaParam(:);
hyp.fdRef = fdRef;
hyp.fdRate = fdRate;
hyp.fdSat = fdSat;
hyp.fdRateSat = fdRateSat;
end

function doaState = localBuildDoaState(model, doaParam)
%LOCALBUILDDOASTATE Map continuous DoA parameters to dynamic geometry.

doaParam = reshape(doaParam, 2, 1);

switch model.doaType
  case 'angle'
    doaState = localBuildAngleState(model, doaParam);

  case 'latlon'
    doaState = localBuildLatlonState(model, doaParam);

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end
end

function doaState = localBuildAngleState(model, doaParam)
%LOCALBUILDANGLESTATE Build frame-wise states for angle-mode parameterization.

eciAngle = doaParam;
eciAngle(1) = mod(eciAngle(1), 2*pi);
eciDirection = doa2dir(eciAngle);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refRotMat = model.rotMatCell{model.steeringRefFrameIdx};
    localDoaCell = eciToAngleGrid(eciDirection, refRotMat);
    if ~iscell(localDoaCell)
      localDoaCell = {localDoaCell};
    end

    for iFrame = 1:model.numFrame
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end

  case 'framewise'
    for iFrame = 1:model.numFrame
      localDoaCell = eciToAngleGrid(eciDirection, model.rotMatCell{iFrame});
      if ~iscell(localDoaCell)
        localDoaCell = {localDoaCell};
      end

      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

deltaFdSeq = zeros(model.numSat, model.numFrame);
for iFrame = 1:model.numFrame
  relSatVel = model.satVelEci(:, :, iFrame) - model.refVelEci(:, iFrame);
  deltaFdSeq(:, iFrame) = ((eciDirection.' * relSatVel) / model.wavelength).';
end

[deltaFdRef, deltaFdRate] = localFitSatFdLine(model.timeOffsetSec, deltaFdSeq);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = [];
doaState.userPosEci = [];
doaState.userVelEci = [];
doaState.deltaFdSeq = deltaFdSeq;
doaState.deltaFdRef = deltaFdRef;
doaState.deltaFdRate = deltaFdRate;
doaState.localDoaArr = localDoaArr;
end

function doaState = localBuildLatlonState(model, doaParam)
%LOCALBUILDLATLONSTATE Build frame-wise states for latlon parameterization.

latlon = doaParam;
[userPosEci, userVelEci] = localLatlonToUserState(latlon, model.userStateRef);

eciDirection = userPosEci ./ vecnorm(userPosEci, 2, 1);
eciAngle = localDirectionToAngle(eciDirection);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refFrameIdx = model.steeringRefFrameIdx;
    refUserPos = userPosEci(:, refFrameIdx);

    for iSat = 1:model.numSat
      localDoaArr(:, iSat, :) = repmat( ...
        globalToLocalDoa(refUserPos, model.satPosEci(:, iSat, refFrameIdx), ...
        model.rotMatCell{refFrameIdx}{iSat}), ...
        1, 1, model.numFrame);
    end

  case 'framewise'
    for iFrame = 1:model.numFrame
      currentUserPos = userPosEci(:, iFrame);
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = globalToLocalDoa( ...
          currentUserPos, model.satPosEci(:, iSat, iFrame), ...
          model.rotMatCell{iFrame}{iSat});
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

refDopplerState = buildReferenceDopplerState( ...
  model.refScene, model.satPosEci, model.satVelEci, ...
  userPosEci, userVelEci, model.wavelength, model.timeOffsetSec, model.refSatIdxLocal);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = latlon;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFdSeq = refDopplerState.deltaFd;
doaState.deltaFdRef = refDopplerState.deltaFdRef;
doaState.deltaFdRate = refDopplerState.deltaFdRate;
doaState.fdSatSeq = refDopplerState.fdSat;
doaState.fdRefSeq = refDopplerState.fdRef;
doaState.fdRefGeom = refDopplerState.fdRefRefFrame;
doaState.fdSatGeom = refDopplerState.fdSatRefFrame;
doaState.deltaFdGeom = refDopplerState.deltaFdRefFrame;
doaState.localDoaArr = localDoaArr;
end

function [userPosEci, userVelEci] = localLatlonToUserState(latlon, userStateRef)
%LOCALLATLONTOUSERSTATE Convert one ground point [lat; lon] into motion-aware ECI states.

[userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, userStateRef);
end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2*pi); ...
         asin(max(min(direction(3, :), 1), -1))];
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
