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
verbose = localResolveVerbose(model);

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

fdRateReleaseHalfWidth = localResolveUnknownWarmAnchorFdRateReleaseHalfWidth(model, lbAnchor(end), ubAnchor(end));
fdRateReleaseSeed = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed, lbAnchor(end), ubAnchor(end), fdRateReleaseHalfWidth);
if isempty(fdRateReleaseSeed)
  fdRateReleaseSeed = fdRateSeed;
end
fdRateReleaseSeed = reshape(fdRateReleaseSeed, 1, []);

if verbose
  fprintf('[WarmAnchorTrace] freezeDoa=%d seedFdRate=%.6f releaseSeedCount=%d\n', ...
    localGetModelLogical(model, 'freezeDoa', false), fdRateSeed, numel(fdRateReleaseSeed));
  localPrintSolveSummary('warm', solveWarm);
end

numReleaseSeed = numel(fdRateReleaseSeed);
seedSolveCell = cell(numReleaseSeed, 1);
useReleaseParfor = localShouldUseReleaseParfor(model, numReleaseSeed, verbose);

if useReleaseParfor
  parfor iSeed = 1:numReleaseSeed
    seedSolveCell{iSeed} = localRunUnknownWarmReleaseSeed(model, warmInit, ...
      lbAnchor, ubAnchor, fdRateSeed, fdRateReleaseSeed, iSeed, ...
      solveWarm, baseOptimOpt, warmSolveOpt);
  end
else
  for iSeed = 1:numReleaseSeed
    seedSolveCell{iSeed} = localRunUnknownWarmReleaseSeed(model, warmInit, ...
      lbAnchor, ubAnchor, fdRateSeed, fdRateReleaseSeed, iSeed, ...
      solveWarm, baseOptimOpt, warmSolveOpt);

    if verbose
      seedSolve = seedSolveCell{iSeed};
      fprintf('[WarmAnchorTrace] releaseSeed[%d/%d]=%.6f halfWidth=%.6f\n', ...
        iSeed, numReleaseSeed, fdRateReleaseSeed(iSeed), seedSolve.releaseHalfWidth);
      localPrintSolveSummary('local', seedSolve.solveLocal);
      localPrintSolveSummary('full', seedSolve.solveFull);
      localPrintSolveSummary('pick', seedSolve.solveCandUse);
    end
  end
end

solveBest = struct([]);
totalRunTimeSec = solveWarm.runTimeSec;
for iSeed = 1:numReleaseSeed
  seedSolve = seedSolveCell{iSeed};
  if isempty(seedSolve) || ~isstruct(seedSolve)
    continue;
  end
  totalRunTimeSec = totalRunTimeSec + seedSolve.runTimeSec;
  solveCandUse = seedSolve.solveCandUse;
  if isempty(solveBest) || preferDoaDopplerMfSolveResult(solveCandUse, solveBest)
    solveBest = solveCandUse;
  end
end

if isempty(solveBest)
  return;
end

solveBest.runTimeSec = totalRunTimeSec;
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  solveBest.solveVariant = "mainWarmAnchorFixedDoa";
else
  solveBest.solveVariant = "mainWarmAnchor";
end
solveBest.warmStart.fdRateReleaseSeedList = fdRateReleaseSeed(:);
if verbose
  localPrintSolveSummary('familyFinal', solveBest);
end
solveCand = solveBest;
end


function seedSolve = localRunUnknownWarmReleaseSeed(model, warmInit, lbAnchor, ubAnchor, ...
  fdRateSeed, fdRateReleaseSeed, iSeed, solveWarm, baseOptimOpt, warmSolveOpt)
%LOCALRUNUNKNOWNWARMRELEASESEED Solve one independent warm-anchor release seed.

initRelease = warmInit;
initRelease(end) = fdRateReleaseSeed(iSeed);

[lbLocal, ubLocal, releaseHalfWidth] = localBuildUnknownFdRateReleaseBounds( ...
  model, lbAnchor, ubAnchor, fdRateReleaseSeed(iSeed), fdRateReleaseSeed);
solveLocal = runDoaDopplerMfOptimization(model, initRelease, lbLocal, ubLocal, [], [], ...
  baseOptimOpt, warmSolveOpt);
solveFull = runDoaDopplerMfOptimization(model, solveLocal.optVar(:), lbAnchor, ubAnchor, [], [], ...
  baseOptimOpt, warmSolveOpt);
releaseRunTimeSec = solveLocal.runTimeSec + solveFull.runTimeSec;

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

seedSolve = struct();
seedSolve.solveLocal = solveLocal;
seedSolve.solveFull = solveFull;
seedSolve.solveCandUse = solveCandUse;
seedSolve.releaseHalfWidth = releaseHalfWidth;
seedSolve.runTimeSec = releaseRunTimeSec;
end


function tf = localShouldUseReleaseParfor(model, numReleaseSeed, verbose)
%LOCALSHOULDUSERELEASEPARFOR Decide whether the warm-anchor family should use parfor.

tf = shouldUseDoaDopplerMfWarmAnchorParfor(model, numReleaseSeed, verbose);
end


function fdRateSeedList = localBuildUnknownFdRateWarmStartSet(model, fdRateSeed, fdRateLb, fdRateUb, releaseHalfWidth)
%LOCALBUILDUNKNOWNFDRATEWARMSTARTSET Build one compact multi-seed release set.
%
% The default CP-U family should not silently collapse to one warm-rate seed.
% When the caller does not provide explicit offsets, build a compact
% symmetric multi-start set around the warm seed so the local release can
% actually leave the CP-K tooth. Explicit offsets still override this
% default and preserve the same data flow through modelOpt.

fdRateSeedList = [];
if ~isfinite(fdRateSeed)
  return;
end

offsetList = reshape(localGetStructField(model, 'unknownWarmAnchorFdRateReleaseOffsetList', []), [], 1);
if isempty(offsetList)
  if ~(isfinite(releaseHalfWidth) && releaseHalfWidth > 0)
    fdRateSeedList = fdRateSeed;
    return;
  end
  offsetList = releaseHalfWidth * [-1.0; -0.5; 0.0; 0.5; 1.0];
end

offsetList = offsetList(isfinite(offsetList));
if isempty(offsetList)
  fdRateSeedList = fdRateSeed;
  return;
end

fdRateSeedList = fdRateSeed + offsetList;
fdRateSeedList = [fdRateSeed; fdRateSeedList];
if isfinite(fdRateLb)
  fdRateSeedList = max(fdRateSeedList, fdRateLb);
end
if isfinite(fdRateUb)
  fdRateSeedList = min(fdRateSeedList, fdRateUb);
end
fdRateSeedList = unique(fdRateSeedList(isfinite(fdRateSeedList)), 'stable');
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
if numel(lbLocal) < 4 || numel(ubLocal) < 4
  releaseHalfWidth = 0;
  return;
end

fdRateIdx = 4;
releaseHalfWidth = localResolveUnknownWarmAnchorFdRateReleaseHalfWidth(model, lb(fdRateIdx), ub(fdRateIdx));
if releaseHalfWidth <= 0
  return;
end
fdRateCenter = fdRateCand;
fdRateLb = max(lb(fdRateIdx), fdRateCenter - releaseHalfWidth);
fdRateUb = min(ub(fdRateIdx), fdRateCenter + releaseHalfWidth);
if fdRateLb > fdRateUb
  fdRateLb = lb(fdRateIdx);
  fdRateUb = ub(fdRateIdx);
end
% Keep overlapping local boxes. The CP-U family is intentionally tiny, and
% allowing overlap is preferable to shrinking every release box back into a
% near-singleton neighborhood around the original seed.
if fdRateLb <= fdRateUb
  lbLocal(fdRateIdx) = fdRateLb;
  ubLocal(fdRateIdx) = fdRateUb;
  releaseHalfWidth = 0.5 * (fdRateUb - fdRateLb);
end
end


function releaseHalfWidth = localResolveUnknownWarmAnchorFdRateReleaseHalfWidth(model, fdRateLb, fdRateUb)
%LOCALRESOLVEUNKNOWNWARMANCHORFDRATERELEASEHALFWIDTH Resolve one nondegenerate CP-U rate-release box.
%
% The default must stay local, but it must not silently collapse to zero.
% A compact half width of about 250 Hz/s is enough to release the biased
% CP-K seeds used by the regression guards while still staying well inside
% the full fdRate search range.

releaseHalfWidth = localGetModelScalar(model, 'unknownWarmAnchorFdRateReleaseHalfWidth', NaN);
if isfinite(releaseHalfWidth)
  releaseHalfWidth = max(releaseHalfWidth, 0);
else
  fdRateSpan = fdRateUb - fdRateLb;
  if ~(isfinite(fdRateSpan) && fdRateSpan > 0)
    releaseHalfWidth = 0;
    return;
  end
  releaseHalfWidth = min(max(250, 0.01 * fdRateSpan), 0.45 * fdRateSpan);
end
end


function [doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb, doaCenter)
%LOCALAPPLYUNKNOWNDOARELEASEFLOOR Apply the minimum CP-U DoA release box.

if localGetModelLogical(model, 'freezeDoa', false)
  return;
end
if localGetModelLogical(model, 'disableUnknownDoaReleaseFloor', false)
  return;
end
if isempty(doaLb) || isempty(doaUb)
  return;
end

releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model, doaLb, doaUb);
if isempty(releaseHalfWidth)
  return;
end

doaCenter = reshape(doaCenter, [], 1);
numDoa = min([numel(doaLb), numel(doaUb), numel(doaCenter), numel(releaseHalfWidth)]);
if numDoa <= 0
  return;
end

for iDoa = 1:numDoa
  centerUse = doaCenter(iDoa);
  lbTarget = centerUse - releaseHalfWidth(iDoa);
  ubTarget = centerUse + releaseHalfWidth(iDoa);
  doaLb(iDoa) = min(doaLb(iDoa), lbTarget);
  doaUb(iDoa) = max(doaUb(iDoa), ubTarget);
end
end


function releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model, doaLb, doaUb)
%LOCALRESOLVEUNKNOWNDOARELEASEHALFWIDTH Resolve the minimum DoA release box.

releaseHalfWidth = reshape(localGetStructField(model, 'unknownDoaReleaseHalfWidth', []), [], 1);
if isempty(releaseHalfWidth)
  if strcmp(model.doaType, 'latlon')
    releaseHalfWidth = [0.03; 0.03];
  else
    releaseHalfWidth = deg2rad([0.03; 0.03]);
  end
end
if isscalar(releaseHalfWidth)
  releaseHalfWidth = repmat(releaseHalfWidth, numel(doaLb), 1);
end
numDoa = min([numel(releaseHalfWidth), numel(doaLb), numel(doaUb)]);
releaseHalfWidth = releaseHalfWidth(1:numDoa);
end


function outputUse = localMergeOptimOutput(outputA, outputB)
%LOCALMERGEOPTIMOUTPUT Merge one warm-start stage into the final solver stats.

outputUse = outputB;
if ~isstruct(outputUse)
  outputUse = struct();
end
for fieldName = ["iterations", "funcCount"]
  fieldChar = char(fieldName);
  valueA = localGetOutputScalar(outputA, fieldChar, 0);
  valueB = localGetOutputScalar(outputUse, fieldChar, 0);
  outputUse.(fieldChar) = valueA + valueB;
end
for fieldName = ["constrviolation", "stepsize", "firstorderopt"]
  fieldChar = char(fieldName);
  valueB = localGetOutputScalar(outputUse, fieldChar, NaN);
  if ~isfinite(valueB)
    outputUse.(fieldChar) = localGetOutputScalar(outputA, fieldChar, NaN);
  end
end
if ~isfield(outputUse, 'algorithm') || isempty(outputUse.algorithm)
  outputUse.algorithm = localGetOutputString(outputA, 'algorithm', '');
end
if ~isfield(outputUse, 'message') || isempty(outputUse.message)
  outputUse.message = localGetOutputString(outputA, 'message', '');
end
outputUse.usedScaledSolve = localGetOutputLogical(outputUse, 'usedScaledSolve', false) || ...
  localGetOutputLogical(outputA, 'usedScaledSolve', false);
outputUse.fallbackTriggered = localGetOutputLogical(outputUse, 'fallbackTriggered', false) || ...
  localGetOutputLogical(outputA, 'fallbackTriggered', false);
if ~isfield(outputUse, 'usedFallbackAlgorithm') || isempty(outputUse.usedFallbackAlgorithm)
  outputUse.usedFallbackAlgorithm = localGetOutputString(outputA, 'usedFallbackAlgorithm', '');
end
end


function value = localGetOutputScalar(outputStruct, fieldName, defaultValue)
%LOCALGETOUTPUTSCALAR Read one optimization-output scalar.

value = defaultValue;
if ~isstruct(outputStruct) || ~isfield(outputStruct, fieldName)
  return;
end
rawValue = outputStruct.(fieldName);
if isempty(rawValue) || ~isscalar(rawValue) || ~isfinite(rawValue)
  return;
end
value = rawValue;
end


function value = localGetOutputString(outputStruct, fieldName, defaultValue)
%LOCALGETOUTPUTSTRING Read one optimization-output string.

value = defaultValue;
if ~isstruct(outputStruct) || ~isfield(outputStruct, fieldName)
  return;
end
rawValue = outputStruct.(fieldName);
if isempty(rawValue)
  return;
end
value = rawValue;
end


function value = localGetOutputLogical(outputStruct, fieldName, defaultValue)
%LOCALGETOUTPUTLOGICAL Read one optimization-output logical flag.

value = defaultValue;
if ~isstruct(outputStruct) || ~isfield(outputStruct, fieldName)
  return;
end
rawValue = outputStruct.(fieldName);
if isempty(rawValue) || ~isscalar(rawValue)
  return;
end
value = logical(rawValue);
end


function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct field with default fallback.

fieldValue = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    fieldValue = rawValue;
  end
end
end


function value = localGetModelLogical(model, fieldName, defaultValue)
%LOCALGETMODELLOGICAL Read one logical model flag with default fallback.

value = defaultValue;
if ~isstruct(model) || ~isfield(model, fieldName) || isempty(model.(fieldName))
  return;
end
value = logical(model.(fieldName));
end


function verbose = localResolveVerbose(model)
%LOCALRESOLVEVERBOSE Resolve one helper-local verbose flag.

verbose = false;
if isstruct(model) && isfield(model, 'verbose') && ~isempty(model.verbose)
  verbose = logical(model.verbose);
end
verbose = verbose || localResolveInheritedVerbose();
end


function localPrintSolveSummary(tag, solveInfo)
%LOCALPRINTSOLVESUMMARY Print one compact solve summary for warm-anchor tracing.

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
  localGetOutputScalar(localGetStructField(solveInfo, 'output', struct()), 'iterations', -1));
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

