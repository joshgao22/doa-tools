function solveCand = runDoaDopplerMfUnknownWarmAnchor(model, initParam, lb, ub, baseOptimOpt)
%RUNDOADOPPLERMFUNKNOWNWARMANCHOR Stabilize CP-U with one warm-rate anchor.
% Keep the active CP-U logic compact: first hold fdRate at the current
% seed long enough to pull DoA into a better basin, then release fdRate in
% a compact local neighborhood. This helper only owns the unknown-route
% family; the top-level branch solver decides whether the family result wins.

arguments
  model (1,1) struct
  initParam (:,1) double
  lb (:,1) double
  ub (:,1) double
  baseOptimOpt
end

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
warmSolveOpt = struct( ...
  'useScaledSolve', localGetModelScalar(model, 'unknownWarmAnchorUseScaledSolve', false), ...
  'enableFallbackSqp', localGetModelScalar(model, 'unknownWarmAnchorFallbackSqp', false));
solveWarm = runDoaDopplerMfOptimization(model, initParam, lb, ub, Aeq, beq, baseOptimOpt, warmSolveOpt);

warmInit = solveWarm.optVar(:);
[lbAnchor, ubAnchor] = buildDoaDopplerMfBounds(model);
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
  solveLocal = runDoaDopplerMfOptimization(model, initRelease, lbLocal, ubLocal, [], [], baseOptimOpt, warmSolveOpt);
  solveFull = runDoaDopplerMfOptimization(model, solveLocal.optVar(:), lbAnchor, ubAnchor, [], [], baseOptimOpt, warmSolveOpt);

  totalRunTimeSec = totalRunTimeSec + solveLocal.runTimeSec + solveFull.runTimeSec;
  solveLocal.output = localMergeOptimOutput(solveWarm.output, solveLocal.output);
  solveLocal.warmStart = struct( ...
    'fdRateFixed', fdRateSeed, ...
    'warmFval', solveWarm.fval, ...
    'warmExitflag', solveWarm.exitflag, ...
    'warmIsResolved', solveWarm.isResolved, ...
    'warmVariant', "mainWarmFixedFdRate", ...
    'releaseSeedFdRate', fdRateReleaseSeed(iSeed), ...
    'releaseHalfWidth', releaseHalfWidth, ...
    'releaseLb', lbLocal(end), ...
    'releaseUb', ubLocal(end), ...
    'selectedReleaseStage', "local");
  solveLocal.runTimeSec = solveWarm.runTimeSec + solveLocal.runTimeSec;

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
    'releaseUb', ubLocal(end), ...
    'selectedReleaseStage', "full");

  solveCandUse = localSelectUnknownWarmReleaseSolve(model, solveLocal, solveFull);
  if isempty(solveBest) || preferDoaDopplerMfSolveResult(solveCandUse, solveBest)
    solveBest = solveCandUse;
  end
end

if isempty(solveBest)
  return;
end

solveBest.runTimeSec = totalRunTimeSec;
solveBest.output = localMergeOptimOutput(solveWarm.output, solveBest.output);
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  solveBest.solveVariant = "mainWarmAnchorFixedDoa";
else
  solveBest.solveVariant = "mainWarmAnchor";
end
solveBest.warmStart.fdRateReleaseSeedList = fdRateReleaseSeed(:);
solveCand = solveBest;
end


function solveUse = localSelectUnknownWarmReleaseSolve(model, solveLocal, solveFull)
%LOCALSELECTUNKNOWNWARMRELEASESOLVE Keep the anchor-preserving release when warranted.

solveUse = solveFull;
if ~isstruct(solveLocal) || isempty(solveLocal)
  return;
end
if ~isstruct(solveFull) || isempty(solveFull)
  solveUse = solveLocal;
  return;
end
if solveLocal.isResolved && ~solveFull.isResolved
  solveUse = solveLocal;
  return;
end
if ~solveLocal.isResolved || ~solveFull.isResolved
  return;
end

objLocal = solveLocal.fval;
objFull = solveFull.fval;
residualLocal = localGetEvalDiagField(solveLocal.finalEvalDiag, 'residualNorm', NaN);
residualFull = localGetEvalDiagField(solveFull.finalEvalDiag, 'residualNorm', NaN);
objGain = objLocal - objFull;
residualGain = residualLocal - residualFull;
doaMoveDeg = localCalcDoaMoveDeg(model, solveLocal.optVar, solveFull.optVar);
fitFloorLocal = localGetEvalDiagField(solveLocal.finalEvalDiag, 'nonRefFitRatioFloor', NaN);
fitFloorFull = localGetEvalDiagField(solveFull.finalEvalDiag, 'nonRefFitRatioFloor', NaN);
supportFloorLocal = localGetEvalDiagField(solveLocal.finalEvalDiag, 'nonRefSupportRatioFloor', NaN);
supportFloorFull = localGetEvalDiagField(solveFull.finalEvalDiag, 'nonRefSupportRatioFloor', NaN);

minObjGain = max(localGetModelScalar(model, 'unknownWarmAnchorMinObjGain', 5), ...
  1e-6 * max([1, abs(objLocal), abs(objFull)]));
minResidualGain = max(localGetModelScalar(model, 'unknownWarmAnchorMinResidualGain', 1), ...
  1e-6 * max([1, abs(residualLocal), abs(residualFull)]));
maxDoaMoveDeg = localGetModelScalar(model, 'unknownWarmAnchorMaxDoaMoveDeg', 1e-3);
qualityTol = 1e-3;
qualityDegraded = false;
if isfinite(fitFloorLocal) && isfinite(fitFloorFull) && fitFloorFull + qualityTol < fitFloorLocal
  qualityDegraded = true;
end
if isfinite(supportFloorLocal) && isfinite(supportFloorFull) && supportFloorFull + qualityTol < supportFloorLocal
  qualityDegraded = true;
end
smallGain = (objGain < minObjGain) && (residualGain < minResidualGain);
largeDoaMove = isfinite(doaMoveDeg) && doaMoveDeg > maxDoaMoveDeg;
if smallGain && (largeDoaMove || qualityDegraded)
  solveLocal.solveVariant = solveFull.solveVariant;
  solveLocal.warmStart = localGetStructField(solveLocal, 'warmStart', struct());
  solveLocal.warmStart.selectedReleaseStage = "local";
  solveLocal.warmStart.objGainToFull = objGain;
  solveLocal.warmStart.residualGainToFull = residualGain;
  solveLocal.warmStart.doaMoveDegToFull = doaMoveDeg;
  solveLocal.warmStart.fullNonRefFitRatioFloor = fitFloorFull;
  solveLocal.warmStart.fullNonRefSupportRatioFloor = supportFloorFull;
  solveUse = solveLocal;
else
  solveFull.warmStart = localGetStructField(solveFull, 'warmStart', struct());
  solveFull.warmStart.selectedReleaseStage = "full";
  solveFull.warmStart.objGainFromLocal = objGain;
  solveFull.warmStart.residualGainFromLocal = residualGain;
  solveFull.warmStart.doaMoveDegFromLocal = doaMoveDeg;
  solveFull.warmStart.localNonRefFitRatioFloor = fitFloorLocal;
  solveFull.warmStart.localNonRefSupportRatioFloor = supportFloorLocal;
  solveUse = solveFull;
end
end


function value = localGetEvalDiagField(evalDiag, fieldName, defaultValue)
%LOCALGETEVALDIAGFIELD Read one scalar evaluation-diagnostic field.

value = defaultValue;
if ~isstruct(evalDiag) || ~isfield(evalDiag, fieldName)
  return;
end
rawValue = evalDiag.(fieldName);
if isempty(rawValue)
  return;
end
rawValue = rawValue(1);
if isfinite(rawValue)
  value = rawValue;
end
end


function value = localGetModelScalar(model, fieldName, defaultValue)
%LOCALGETMODELSCALAR Read one finite scalar model option with fallback.

value = defaultValue;
if ~isstruct(model) || ~isfield(model, fieldName)
  return;
end
rawValue = model.(fieldName);
if isempty(rawValue) || ~isscalar(rawValue) || ~isfinite(rawValue)
  return;
end
value = rawValue;
end


function doaMoveDeg = localCalcDoaMoveDeg(model, optVarA, optVarB)
%LOCALCALCDOAMOVEDEG Measure DoA drift between two solve candidates.

doaMoveDeg = NaN;
try
  if strcmp(model.doaType, 'latlon')
    doaMoveDeg = calcLatlonAngleError(optVarA(1:2), optVarB(1:2));
  else
    doaMoveDeg = rad2deg(norm(optVarA(1:2) - optVarB(1:2)));
  end
catch
  doaMoveDeg = NaN;
end
end


function [lbLocal, ubLocal, releaseHalfWidth] = localBuildUnknownFdRateReleaseBounds(model, lb, ub, fdRateCand, fdRateCandList)
%LOCALBUILDUNKNOWNFDRATERELEASEBOUNDS Build one compact local fdRate box.

lbLocal = lb;
ubLocal = ub;
releaseHalfWidth = localGetModelScalar(model, 'unknownWarmAnchorFdRateReleaseHalfWidth', 0);
if numel(lbLocal) < 4 || numel(ubLocal) < 4
  return;
end
if releaseHalfWidth <= 0
  return;
end
fdRateIdx = 4;
fdRateCenter = fdRateCand;
fdRateLb = max(lb(fdRateIdx), fdRateCenter - releaseHalfWidth);
fdRateUb = min(ub(fdRateIdx), fdRateCenter + releaseHalfWidth);
if fdRateLb > fdRateUb
  fdRateLb = lb(fdRateIdx);
  fdRateUb = ub(fdRateIdx);
end
if nargin >= 5 && ~isempty(fdRateCandList)
  nearestGap = min(abs(fdRateCenter - fdRateCandList(abs(fdRateCenter - fdRateCandList) > 0)));
  if isfinite(nearestGap)
    halfGap = 0.45 * nearestGap;
    fdRateLb = max(fdRateLb, fdRateCenter - halfGap);
    fdRateUb = min(fdRateUb, fdRateCenter + halfGap);
  end
end
if fdRateLb <= fdRateUb
  lbLocal(fdRateIdx) = fdRateLb;
  ubLocal(fdRateIdx) = fdRateUb;
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
outputMerge.usedScaledSolve = logical(localGetStructField(outputWarm, 'usedScaledSolve', false) || ...
  localGetStructField(outputRelease, 'usedScaledSolve', false));
outputMerge.fallbackTriggered = logical(localGetStructField(outputWarm, 'fallbackTriggered', false) || ...
  localGetStructField(outputRelease, 'fallbackTriggered', false));
fallbackAlgRelease = localGetStructField(outputRelease, 'usedFallbackAlgorithm', '');
fallbackAlgWarm = localGetStructField(outputWarm, 'usedFallbackAlgorithm', '');
if ~isempty(fallbackAlgRelease)
  outputMerge.usedFallbackAlgorithm = fallbackAlgRelease;
else
  outputMerge.usedFallbackAlgorithm = fallbackAlgWarm;
end
end


function fdRateCand = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed)
%LOCALBUILDUNKNOWNFDRATEWARMSTARTSET Build local fdRate anchors for CP-U.

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


function [doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb, centerOverride)
%LOCALAPPLYUNKNOWNDOARELEASEFLOOR Local copy for warm-anchor recentering.

if nargin < 4
  centerOverride = [];
end
if ~strcmp(model.fdRateMode, 'unknown')
  return;
end
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  return;
end
if isfield(model, 'disableUnknownDoaReleaseFloor') && logical(model.disableUnknownDoaReleaseFloor)
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
doAResetLb = reshape(baseRange(invalidMask, 1), [], 1);
doAResetUb = reshape(baseRange(invalidMask, 2), [], 1);
doAIdx = find(invalidMask);
for iIdx = 1:numel(doAIdx)
  doaLb(doAIdx(iIdx)) = doAResetLb(iIdx);
  doaUb(doAIdx(iIdx)) = doAResetUb(iIdx);
end
end


function releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model)
%LOCALRESOLVEUNKNOWNDOARELEASEHALFWIDTH Resolve the CP-U DoA box floor.

if isfield(model, 'unknownDoaReleaseHalfWidth') && ~isempty(model.unknownDoaReleaseHalfWidth)
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
if numel(releaseHalfWidth) ~= 2 || any(~isfinite(releaseHalfWidth)) || any(releaseHalfWidth < 0)
  releaseHalfWidth = [];
end
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
