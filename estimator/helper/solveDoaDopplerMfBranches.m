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
end

[solveUse, doaBasinEntry] = localRunDoaBasinEntry( ...
  model, initParam, solveUse, baseOptimOpt, verbose);

if strcmp(model.fdRateMode, 'unknown')
  if isfield(model, 'disableUnknownWarmAnchor') && logical(model.disableUnknownWarmAnchor)
    solveCand = struct([]);
    if verbose
      fprintf('[MfBranchTrace] unknown warm-anchor disabled.\n');
    end
  else
    warmModel = localBuildWarmAnchorModel(model, solveUse);
    [warmLb, warmUb] = buildDoaDopplerMfBounds(warmModel);
    warmInit = localClipParamToBounds(solveUse.optVar, warmLb, warmUb, warmModel);
    solveCand = runDoaDopplerMfUnknownWarmAnchor(warmModel, warmInit, warmLb, warmUb, baseOptimOpt);
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
solveUse.doaBasinEntry = doaBasinEntry;

optimInfo = localBuildOptimInfo(solveUse.output, solveUse.exitflag, ...
  solveUse.runTimeSec, baseOptimOpt, model);
optimInfo.solveVariant = solveUse.solveVariant;
optimInfo.verbose = verbose;
if isfield(solveUse, 'candidateObjective')
  optimInfo.candidateObjective = solveUse.candidateObjective;
  optimInfo.candidateExitflag = solveUse.candidateExitflag;
  optimInfo.candidateVariant = solveUse.candidateVariant;
end
if isfield(solveUse, 'doaBasinEntry')
  optimInfo.doaBasinEntry = solveUse.doaBasinEntry;
end
if verbose
  fprintf('[WarmAnchorExit] freezeDoa=%d final variant=%s fval=%.6e resolved=%d\n', ...
    localGetModelLogical(model, 'freezeDoa', false), ...
    string(localGetStructField(solveUse, 'solveVariant', "")), ...
    localGetStructField(solveUse, 'fval', NaN), ...
    localGetStructField(solveUse, 'isResolved', false));
end
end


function [solveUse, basinDiag] = localRunDoaBasinEntry(model, initParam, solveUse, baseOptimOpt, verbose)
%LOCALRUNDOABASINENTRY Run truth-free wide DoA acquisition before final polish.
% The acquisition uses the same MF objective and keeps fd/fdRate ranges
% unchanged. Only the DoA search box is widened temporarily, then the winner
% is polished again inside the normal compact local box.

basinDiag = localInitDoaBasinEntryDiag(model);
[shouldRunEntry, skipReason] = localShouldRunDoaBasinEntry(model);
if ~shouldRunEntry
  basinDiag.reason = skipReason;
  return;
end

halfWidthMat = localResolveDoaBasinEntryHalfWidths(model);
if isempty(halfWidthMat)
  basinDiag.reason = "no-entry-half-width";
  return;
end

baseCenter = reshape(model.initDoaParam, [], 1);
compactHalfWidth = reshape(model.initDoaHalfWidth, [], 1);
solveOpt = struct('useScaledSolve', false, 'enableFallbackSqp', false);
entryRows = repmat(localEmptyDoaBasinEntryRow(), 0, 1);
bestSolve = solveUse;
bestTag = "baseline";
bestEntryIdx = 0;

for iEntry = 1:size(halfWidthMat, 2)
  entryHalfWidth = halfWidthMat(:, iEntry);
  entryTag = sprintf('doaBasinEntry%.4g', localHalfWidthDisplayValue(model, entryHalfWidth));
  entryModel = localBuildDoaBoxModel(model, baseCenter, entryHalfWidth);
  [entryLb, entryUb] = buildDoaDopplerMfBounds(entryModel);
  entryInit = localClipParamToBounds(initParam, entryLb, entryUb, entryModel);

  entrySolve = runDoaDopplerMfOptimization( ...
    entryModel, entryInit, entryLb, entryUb, [], [], baseOptimOpt, solveOpt);
  entrySolve.solveVariant = string(entryTag);

  adopted = preferDoaDopplerMfSolveResult(entrySolve, bestSolve);
  entryRows(end + 1, 1) = localBuildDoaBasinEntryRow( ... %#ok<AGROW>
    string(entryTag), entryHalfWidth, entrySolve, entrySolve, adopted);
  if adopted
    bestSolve = entrySolve;
    bestTag = string(entrySolve.solveVariant);
    bestEntryIdx = numel(entryRows);
  end

  if verbose
    localPrintSolveSummary(char(entryTag), entrySolve);
  end
end

if bestEntryIdx > 0
  polishModel = localBuildDoaBoxModel(model, bestSolve.optVar(1:2), compactHalfWidth);
  [polishLb, polishUb] = buildDoaDopplerMfBounds(polishModel);
  polishInit = localClipParamToBounds(bestSolve.optVar, polishLb, polishUb, polishModel);
  polishSolve = runDoaDopplerMfOptimization( ...
    polishModel, polishInit, polishLb, polishUb, [], [], baseOptimOpt, solveOpt);
  polishSolve.solveVariant = bestTag + "Polish";

  entryRows(bestEntryIdx).polishFval = localGetStructField(polishSolve, 'fval', NaN);
  entryRows(bestEntryIdx).polishExitflag = localGetStructField(polishSolve, 'exitflag', NaN);
  entryRows(bestEntryIdx).polishIterations = ...
    localGetOutputField(localGetStructField(polishSolve, 'output', struct()), 'iterations', NaN);
  if preferDoaDopplerMfSolveResult(polishSolve, bestSolve)
    bestSolve = polishSolve;
    bestTag = string(polishSolve.solveVariant);
  end
  entryRows(bestEntryIdx).selectedFval = localGetStructField(bestSolve, 'fval', NaN);
  entryRows(bestEntryIdx).selectedVariant = string(localGetStructField(bestSolve, 'solveVariant', ""));
end

solveUse = bestSolve;
basinDiag.enabled = true;
basinDiag.reason = "evaluated";
basinDiag.bestTag = bestTag;
basinDiag.numCandidate = numel(entryRows);
basinDiag.entryTable = entryRows;
end

function basinDiag = localInitDoaBasinEntryDiag(model)
%LOCALINITDOABASINENTRYDIAG Initialize one compact basin-entry diagnostic.

basinDiag = struct();
basinDiag.enabled = false;
basinDiag.reason = "disabled";
basinDiag.bestTag = "baseline";
basinDiag.numCandidate = 0;
basinDiag.entryTable = repmat(localEmptyDoaBasinEntryRow(), 0, 1);
basinDiag.phaseMode = string(localGetStructField(model, 'phaseMode', ""));
basinDiag.fdRateMode = string(localGetStructField(model, 'fdRateMode', ""));
basinDiag.scope = string(localGetStructField(model, 'doaBasinEntryScope', "single-sat"));
basinDiag.numSat = localGetStructField(model, 'numSat', NaN);
end

function [tf, reason] = localShouldRunDoaBasinEntry(model)
%LOCALSHOULDRUNDOABASINENTRY Decide whether the formal MF local solve needs DoA acquisition.

tf = false;
reason = "disabled";
if localGetModelLogical(model, 'disableDoaBasinEntry', false)
  reason = "disabled-by-option";
  return;
end
scope = string(localGetStructField(model, 'doaBasinEntryScope', "single-sat"));
if scope == "off"
  reason = "scope-off";
  return;
end
if scope == "single-sat" && localGetStructField(model, 'numSat', Inf) ~= 1
  reason = "scope-multi-sat-disabled";
  return;
end
if localGetModelLogical(model, 'freezeDoa', false)
  reason = "freeze-doa";
  return;
end
if ~strcmp(localGetStructField(model, 'phaseMode', ''), 'continuous')
  reason = "non-continuous-phase";
  return;
end
if isempty(localGetStructField(model, 'initDoaParam', [])) || ...
    isempty(localGetStructField(model, 'initDoaHalfWidth', []))
  reason = "missing-init-doa-box";
  return;
end
if numel(model.initDoaParam) ~= 2 || numel(model.initDoaHalfWidth) ~= 2
  reason = "invalid-init-doa-size";
  return;
end
if any(~isfinite(model.initDoaParam(:))) || ...
    any(~isfinite(model.initDoaHalfWidth(:))) || any(model.initDoaHalfWidth(:) <= 0)
  reason = "invalid-init-doa-box";
  return;
end
tf = true;
reason = "ready";
end

function halfWidthMat = localResolveDoaBasinEntryHalfWidths(model)
%LOCALRESOLVEDOABASINENTRYHALFWIDTHS Resolve wide acquisition boxes.

rawList = localGetStructField(model, 'doaBasinEntryHalfWidthList', []);
if isempty(rawList)
  if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
    rawList = deg2rad([0.012, 0.024, 0.048]);
  else
    rawList = [0.012, 0.024, 0.048];
  end
end

rawList = double(rawList);
if isvector(rawList)
  rawList = reshape(rawList, 1, []);
  halfWidthMat = repmat(rawList, 2, 1);
elseif size(rawList, 1) == 2
  halfWidthMat = rawList;
elseif size(rawList, 2) == 2
  halfWidthMat = rawList.';
else
  halfWidthMat = zeros(2, 0);
end

if isempty(halfWidthMat)
  return;
end
compactHalfWidth = reshape(model.initDoaHalfWidth, [], 1);
keepMask = all(isfinite(halfWidthMat), 1) & all(halfWidthMat > 0, 1) & ...
  any(halfWidthMat > compactHalfWidth + 1e-12, 1);
halfWidthMat = halfWidthMat(:, keepMask);
[~, keepIdx] = unique(round(halfWidthMat.' / 1e-12) * 1e-12, 'rows', 'stable');
halfWidthMat = halfWidthMat(:, keepIdx);
end

function modelBox = localBuildDoaBoxModel(model, doaCenter, doaHalfWidth)
%LOCALBUILDDOABOXMODEL Build one model copy with a different DoA box.

modelBox = model;
modelBox.initDoaParam = reshape(doaCenter, [], 1);
modelBox.initDoaHalfWidth = reshape(doaHalfWidth, [], 1);
if strcmp(localGetStructField(modelBox, 'doaType', ''), 'angle') && ~isempty(modelBox.initDoaParam)
  modelBox.initDoaParam(1) = mod(modelBox.initDoaParam(1), 2 * pi);
end
modelBox.doaLb = [];
modelBox.doaUb = [];
modelBox.lb = [];
modelBox.ub = [];
modelBox.cachedBoundsReady = false;
end

function modelWarm = localBuildWarmAnchorModel(model, solveUse)
%LOCALBUILDWARMANCHORMODEL Anchor CP-U release to the best DoA basin.

modelWarm = model;
if isstruct(solveUse) && isfield(solveUse, 'optVar') && numel(solveUse.optVar) >= 2 && ...
    ~isempty(model.initDoaHalfWidth)
  modelWarm = localBuildDoaBoxModel(model, solveUse.optVar(1:2), model.initDoaHalfWidth);
end
end

function optVar = localClipParamToBounds(optVar, lb, ub, model)
%LOCALCLIPPARAMTOBOUNDS Clip one parameter vector to a target box.

optVar = reshape(optVar, [], 1);
lb = reshape(lb, [], 1);
ub = reshape(ub, [], 1);
optVar = min(max(optVar, lb), ub);
if strcmp(localGetStructField(model, 'doaType', ''), 'angle') && ~isempty(optVar)
  optVar(1) = mod(optVar(1), 2 * pi);
end
end

function row = localEmptyDoaBasinEntryRow()
%LOCALEMPTYDOABASINENTRYROW Build one empty basin-entry diagnostic row.

row = struct( ...
  'tag', "", ...
  'halfWidth1', NaN, ...
  'halfWidth2', NaN, ...
  'entryFval', NaN, ...
  'polishFval', NaN, ...
  'selectedFval', NaN, ...
  'selectedVariant', "", ...
  'entryExitflag', NaN, ...
  'polishExitflag', NaN, ...
  'entryIterations', NaN, ...
  'polishIterations', NaN, ...
  'adoptedOverPreviousBest', false);
end

function row = localBuildDoaBasinEntryRow(tag, halfWidth, entrySolve, selectedSolve, adopted)
%LOCALBUILDDOABASINENTRYROW Build one basin-entry diagnostic row.

row = localEmptyDoaBasinEntryRow();
row.tag = tag;
row.halfWidth1 = halfWidth(1);
row.halfWidth2 = halfWidth(2);
row.entryFval = localGetStructField(entrySolve, 'fval', NaN);
row.selectedFval = localGetStructField(selectedSolve, 'fval', NaN);
row.selectedVariant = string(localGetStructField(selectedSolve, 'solveVariant', ""));
row.entryExitflag = localGetStructField(entrySolve, 'exitflag', NaN);
row.entryIterations = localGetOutputField(localGetStructField(entrySolve, 'output', struct()), 'iterations', NaN);
row.adoptedOverPreviousBest = logical(adopted);
end

function displayValue = localHalfWidthDisplayValue(model, halfWidth)
%LOCALHALFWIDTHDISPLAYVALUE Return a compact half-width value for tags.

if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
  displayValue = rad2deg(max(abs(halfWidth(:))));
else
  displayValue = max(abs(halfWidth(:)));
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
optimInfo.doaBasinEntryScope = string(localGetStructField(model, 'doaBasinEntryScope', "single-sat"));
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

