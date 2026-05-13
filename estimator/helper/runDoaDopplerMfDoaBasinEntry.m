function [solveUse, basinDiag] = runDoaDopplerMfDoaBasinEntry(model, initParam, solveUse, baseOptimOpt, verbose)
%RUNDOADOPPLERMFDOABASINENTRY Run truth-free wide DoA acquisition before final polish.
% The acquisition uses the same MF objective and keeps fd/fdRate ranges
% unchanged. Only the DoA search box is widened temporarily, then the winner
% is polished again inside the normal compact local box.

basinDiag = localInitDoaBasinEntryDiag(model);
basinDiag.baselineFval = localGetStructField(solveUse, 'fval', NaN);
basinDiag.baselineVariant = string(localGetStructField(solveUse, 'solveVariant', ""));
adoptionMode = localResolveDoaBasinEntryAdoptionMode(model);
basinDiag.adoptionMode = adoptionMode;
basinDiag.selectedFval = basinDiag.baselineFval;
basinDiag.selectedVariant = basinDiag.baselineVariant;
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
entrySpecList = localBuildDoaBasinEntrySpecList(model, baseCenter, compactHalfWidth, halfWidthMat);
if isempty(entrySpecList)
  basinDiag.reason = "no-entry-candidate";
  return;
end

numCandidateProposed = numel(entrySpecList);
[entrySpecList, dedupInfo] = localDeduplicateDoaBasinEntrySpecsByActualPath( ...
  model, entrySpecList, initParam);
basinDiag.numCandidateProposed = numCandidateProposed;
basinDiag.numCandidateAfterDedup = numel(entrySpecList);
basinDiag.numCandidateDeduped = dedupInfo.numDeduped;
basinDiag.numCandidateSkippedNoOverlap = dedupInfo.numNoOverlap;
if isempty(entrySpecList)
  basinDiag.reason = "no-entry-parent-overlap";
  return;
end

solveOpt = struct('useScaledSolve', false, 'enableFallbackSqp', false);
entryRows = repmat(localEmptyDoaBasinEntryRow(), 0, 1);
bestSolve = solveUse;
bestTag = "baseline";
bestEntryIdx = 0;

for iEntry = 1:numel(entrySpecList)
  entrySpec = entrySpecList(iEntry);
  if ~localDoaBoxOverlapsActiveEnvelope(model, entrySpec.center, entrySpec.halfWidth)
    continue;
  end
  entryModel = localBuildDoaBoxModel(model, entrySpec.center, entrySpec.halfWidth);
  [entryLb, entryUb] = buildDoaDopplerMfBounds(entryModel);
  entryInit = localBuildDoaBasinEntryInit(initParam, entrySpec.center, entryLb, entryUb, entryModel);

  entrySolve = runDoaDopplerMfOptimization( ...
    entryModel, entryInit, entryLb, entryUb, [], [], baseOptimOpt, solveOpt);
  entrySolve.solveVariant = string(entrySpec.tag);

  previousBestFval = localGetStructField(bestSolve, 'fval', NaN);
  baselineFval = localGetStructField(solveUse, 'fval', NaN);
  adopted = adoptionMode == "objective" && preferDoaDopplerMfSolveResult(entrySolve, bestSolve);
  entryRows(end + 1, 1) = localBuildDoaBasinEntryRow( ... %#ok<AGROW>
    entrySpec, entrySolve, entrySolve, adopted, previousBestFval, baselineFval, model);
  if adopted
    bestSolve = entrySolve;
    bestTag = string(entrySolve.solveVariant);
    bestEntryIdx = numel(entryRows);
  end

  if verbose
    localPrintSolveSummary(char(entrySpec.tag), entrySolve);
  end
end

if isempty(entryRows)
  basinDiag.reason = "no-entry-parent-overlap";
  return;
end

if bestEntryIdx > 0
  polishModel = localBuildDoaBoxModel(model, bestSolve.optVar(1:2), compactHalfWidth);
  [polishLb, polishUb] = buildDoaDopplerMfBounds(polishModel);
  polishInit = localClipParamToBounds(bestSolve.optVar, polishLb, polishUb, polishModel);
  polishSolve = runDoaDopplerMfOptimization( ...
    polishModel, polishInit, polishLb, polishUb, [], [], baseOptimOpt, solveOpt);
  polishSolve.solveVariant = bestTag + "Polish";

  entryRows(bestEntryIdx).polishFval = localGetStructField(polishSolve, 'fval', NaN);
  [entryRows(bestEntryIdx).polishLatDeg, entryRows(bestEntryIdx).polishLonDeg, ...
    entryRows(bestEntryIdx).polishFdRefHz, entryRows(bestEntryIdx).polishFdRateHzPerSec] = ...
    localExtractDoaBasinOptVarParam(polishSolve, model);
  entryRows(bestEntryIdx).polishObjDeltaToBaseline = ...
    entryRows(bestEntryIdx).polishFval - entryRows(bestEntryIdx).baselineFval;
  entryRows(bestEntryIdx).polishExitflag = localGetStructField(polishSolve, 'exitflag', NaN);
  entryRows(bestEntryIdx).polishIterations = ...
    localGetOutputField(localGetStructField(polishSolve, 'output', struct()), 'iterations', NaN);
  if preferDoaDopplerMfSolveResult(polishSolve, bestSolve)
    bestSolve = polishSolve;
    bestTag = string(polishSolve.solveVariant);
  end
  entryRows(bestEntryIdx).selectedFval = localGetStructField(bestSolve, 'fval', NaN);
  [entryRows(bestEntryIdx).selectedLatDeg, entryRows(bestEntryIdx).selectedLonDeg, ...
    entryRows(bestEntryIdx).selectedFdRefHz, entryRows(bestEntryIdx).selectedFdRateHzPerSec] = ...
    localExtractDoaBasinOptVarParam(bestSolve, model);
  entryRows(bestEntryIdx).selectedObjDeltaToBaseline = ...
    entryRows(bestEntryIdx).selectedFval - entryRows(bestEntryIdx).baselineFval;
  entryRows(bestEntryIdx).selectedVariant = string(localGetStructField(bestSolve, 'solveVariant', ""));
  entryRows(bestEntryIdx).selectedOverBaseline = localIsObjectiveBetter( ...
    entryRows(bestEntryIdx).selectedFval, entryRows(bestEntryIdx).baselineFval);
end

solveUse = bestSolve;
basinDiag.enabled = true;
basinDiag.reason = "evaluated";
basinDiag.bestTag = bestTag;
basinDiag.numCandidate = numel(entryRows);
basinDiag.selectedFval = localGetStructField(bestSolve, 'fval', NaN);
basinDiag.selectedVariant = string(localGetStructField(bestSolve, 'solveVariant', ""));
basinDiag.selectedOverBaseline = localIsObjectiveBetter(basinDiag.selectedFval, basinDiag.baselineFval);
basinDiag.entryTable = entryRows;
end

function adoptionMode = localResolveDoaBasinEntryAdoptionMode(model)
%LOCALRESOLVEDOABASINENTRYADOPTIONMODE Resolve candidate adoption behavior.

adoptionMode = string(localGetStructField(model, 'doaBasinEntryAdoptionMode', 'objective'));
adoptionMode = lower(strtrim(adoptionMode));
if ~ismember(adoptionMode, ["objective", "diagnostic-only"])
  adoptionMode = "objective";
end
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
basinDiag.numCandidateProposed = 0;
basinDiag.numCandidateAfterDedup = 0;
basinDiag.numCandidateDeduped = 0;
basinDiag.numCandidateSkippedNoOverlap = 0;
basinDiag.numExternalCenter = size(localNormalizeDoaMatrix( ...
  localGetStructField(model, 'doaBasinEntryCenterList', [])), 2);
basinDiag.numOffset = size(localNormalizeDoaMatrix( ...
  localGetStructField(model, 'doaBasinEntryOffsetList', [])), 2);
basinDiag.adoptionMode = localResolveDoaBasinEntryAdoptionMode(model);
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


function entrySpecList = localBuildDoaBasinEntrySpecList(model, baseCenter, compactHalfWidth, halfWidthMat)
%LOCALBUILDDOABASINENTRYSPECLIST Build truth-free DoA entry boxes.
% Default specs preserve the original wide-box behavior. External center and
% offset specs are opt-in and are used by MS replay/scan diagnostics to test
% truth-free family-aware basin entry without changing the SS default path.

entrySpecList = repmat(localEmptyDoaBasinEntrySpec(), 0, 1);
for iEntry = 1:size(halfWidthMat, 2)
  spec = localEmptyDoaBasinEntrySpec();
  spec.tag = sprintf('doaBasinEntry%.4g', localHalfWidthDisplayValue(model, halfWidthMat(:, iEntry)));
  spec.center = baseCenter(:);
  spec.halfWidth = halfWidthMat(:, iEntry);
  spec.offset = zeros(2, 1);
  spec.centerSource = "base";
  entrySpecList(end + 1, 1) = spec; %#ok<AGROW>
end

centerMat = localNormalizeDoaMatrix(localGetStructField(model, 'doaBasinEntryCenterList', []));
if isempty(centerMat)
  return;
end
centerSourceList = localResolveDoaBasinEntryCenterSourceList(model, size(centerMat, 2));
offsetMat = localNormalizeDoaMatrix(localGetStructField(model, 'doaBasinEntryOffsetList', []));
if isempty(offsetMat)
  offsetMat = zeros(2, 1);
end

for iCenter = 1:size(centerMat, 2)
  for iOffset = 1:size(offsetMat, 2)
    spec = localEmptyDoaBasinEntrySpec();
    spec.center = centerMat(:, iCenter) + offsetMat(:, iOffset);
    spec.center = localWrapDoaCenter(model, spec.center);
    spec.halfWidth = compactHalfWidth(:);
    spec.offset = offsetMat(:, iOffset);
    spec.centerSource = centerSourceList(iCenter);
    spec.tag = localBuildExternalDoaBasinEntryTag(model, iCenter, offsetMat(:, iOffset));
    entrySpecList(end + 1, 1) = spec; %#ok<AGROW>
  end
end

entrySpecList = localDeduplicateDoaBasinEntrySpecs(entrySpecList);
end

function spec = localEmptyDoaBasinEntrySpec()
%LOCALEMPTYDOABASINENTRYSPEC Return one typed DoA entry spec.

spec = struct('tag', "", 'center', nan(2, 1), 'halfWidth', nan(2, 1), ...
  'offset', nan(2, 1), 'centerSource', "");
end

function mat = localNormalizeDoaMatrix(rawMat)
%LOCALNORMALIZEDOAMATRIX Normalize a DoA 2-vector list to 2-by-N.

mat = zeros(2, 0);
if isempty(rawMat) || ~isnumeric(rawMat)
  return;
end
rawMat = double(rawMat);
if isvector(rawMat)
  if numel(rawMat) ~= 2
    return;
  end
  mat = reshape(rawMat, 2, 1);
elseif size(rawMat, 1) == 2
  mat = rawMat;
elseif size(rawMat, 2) == 2
  mat = rawMat.';
end
if isempty(mat)
  return;
end
keepMask = all(isfinite(mat), 1);
mat = mat(:, keepMask);
end

function centerSourceList = localResolveDoaBasinEntryCenterSourceList(model, numCenter)
%LOCALRESOLVEDOABASINENTRYCENTERSOURCELIST Resolve external center labels.

rawList = localGetStructField(model, 'doaBasinEntryCenterSourceList', strings(0, 1));
centerSourceList = strings(numCenter, 1);
for iCenter = 1:numCenter
  centerSourceList(iCenter) = sprintf('external%d', iCenter);
end
if isempty(rawList)
  return;
end
rawList = reshape(string(rawList), [], 1);
numUse = min(numel(rawList), numCenter);
for iCenter = 1:numUse
  if strlength(rawList(iCenter)) > 0
    centerSourceList(iCenter) = rawList(iCenter);
  end
end
end


function entrySpecList = localDeduplicateDoaBasinEntrySpecs(entrySpecList)
%LOCALDEDUPLICATEDOABASINENTRYSPECS Remove numerically duplicated specs.

if isempty(entrySpecList)
  return;
end
keyMat = zeros(numel(entrySpecList), 4);
for iSpec = 1:numel(entrySpecList)
  keyMat(iSpec, :) = [entrySpecList(iSpec).center(:).', entrySpecList(iSpec).halfWidth(:).'];
end
[~, keepIdx] = unique(round(keyMat / 1e-12) * 1e-12, 'rows', 'stable');
entrySpecList = entrySpecList(keepIdx);
end

function [entrySpecList, info] = localDeduplicateDoaBasinEntrySpecsByActualPath( ...
  model, entrySpecList, initParam)
%LOCALDEDUPLICATEDOABASINENTRYSPECSBYACTUALPATH Remove duplicated actual solves.
% Proposed wide boxes can collapse to the same effective path after the
% parent envelope is applied.  Keep only the first candidate with identical
% actual bounds and clipped init so the acquisition order remains stable.

info = struct('numDeduped', 0, 'numNoOverlap', 0);
if isempty(entrySpecList)
  return;
end

actualSpecList = repmat(localEmptyDoaBasinEntrySpec(), 0, 1);
keyMat = [];
for iSpec = 1:numel(entrySpecList)
  spec = entrySpecList(iSpec);
  if ~localDoaBoxOverlapsActiveEnvelope(model, spec.center, spec.halfWidth)
    info.numNoOverlap = info.numNoOverlap + 1;
    continue;
  end
  entryModel = localBuildDoaBoxModel(model, spec.center, spec.halfWidth);
  [entryLb, entryUb] = buildDoaDopplerMfBounds(entryModel);
  entryInit = localBuildDoaBasinEntryInit(initParam, spec.center, entryLb, entryUb, entryModel);
  keyMat(end + 1, :) = localBuildActualDoaBasinEntryKey(entryLb, entryUb, entryInit); %#ok<AGROW>
  actualSpecList(end + 1, 1) = spec; %#ok<AGROW>
end

if isempty(actualSpecList)
  entrySpecList = actualSpecList;
  info.numDeduped = 0;
  return;
end

[~, keepIdx] = unique(round(keyMat / 1e-10) * 1e-10, 'rows', 'stable');
entrySpecList = actualSpecList(keepIdx);
info.numDeduped = numel(actualSpecList) - numel(keepIdx);
end

function key = localBuildActualDoaBasinEntryKey(lb, ub, initParam)
%LOCALBUILDACTUALDOABASINENTRYKEY Build a stable key for one actual solve.

key = [reshape(lb, 1, []), reshape(ub, 1, []), reshape(initParam, 1, [])];
end

function tag = localBuildExternalDoaBasinEntryTag(model, centerIdx, offset)
%LOCALBUILDEXTERNALDOABASINENTRYTAG Build a compact external-center tag.

if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
  offsetDisplay = rad2deg(offset(:));
else
  offsetDisplay = offset(:);
end
tag = sprintf('doaBasinExt%d_%s_%s', centerIdx, ...
  char(localFormatDoaOffsetTag('lat', offsetDisplay(1))), ...
  char(localFormatDoaOffsetTag('lon', offsetDisplay(2))));
end

function tag = localFormatDoaOffsetTag(axisName, value)
%LOCALFORMATDOAOFFSETTAG Format one DoA offset into a tag component.

if abs(value) < 1e-12
  tag = string(axisName) + "0";
  return;
end
if value > 0
  signTag = "p";
else
  signTag = "m";
end
magTag = regexprep(sprintf('%.4g', abs(value)), '\\.', 'p');
tag = string(axisName) + signTag + string(magTag);
end

function optVar = localBuildDoaBasinEntryInit(initParam, doaCenter, lb, ub, model)
%LOCALBUILDDOABASINENTRYINIT Build one candidate init inside the target box.

optVar = reshape(initParam, [], 1);
if numel(optVar) >= 2
  optVar(1:2) = reshape(doaCenter(1:2), [], 1);
end
optVar = localClipParamToBounds(optVar, lb, ub, model);
end

function doaCenter = localWrapDoaCenter(model, doaCenter)
%LOCALWRAPDOACENTER Wrap azimuth-like DoA center when needed.

doaCenter = reshape(doaCenter, [], 1);
if strcmp(localGetStructField(model, 'doaType', ''), 'angle') && ~isempty(doaCenter)
  doaCenter(1) = mod(doaCenter(1), 2 * pi);
end
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
  'center1', NaN, ...
  'center2', NaN, ...
  'offset1', NaN, ...
  'offset2', NaN, ...
  'centerSource', "", ...
  'adoptionMode', "", ...
  'baselineFval', NaN, ...
  'previousBestFval', NaN, ...
  'entryFval', NaN, ...
  'entryLatDeg', NaN, ...
  'entryLonDeg', NaN, ...
  'entryFdRefHz', NaN, ...
  'entryFdRateHzPerSec', NaN, ...
  'entryObjDeltaToPreviousBest', NaN, ...
  'entryObjDeltaToBaseline', NaN, ...
  'polishFval', NaN, ...
  'polishLatDeg', NaN, ...
  'polishLonDeg', NaN, ...
  'polishFdRefHz', NaN, ...
  'polishFdRateHzPerSec', NaN, ...
  'polishObjDeltaToBaseline', NaN, ...
  'selectedFval', NaN, ...
  'selectedLatDeg', NaN, ...
  'selectedLonDeg', NaN, ...
  'selectedFdRefHz', NaN, ...
  'selectedFdRateHzPerSec', NaN, ...
  'selectedObjDeltaToBaseline', NaN, ...
  'selectedVariant', "", ...
  'entryExitflag', NaN, ...
  'polishExitflag', NaN, ...
  'entryIterations', NaN, ...
  'polishIterations', NaN, ...
  'adoptedOverPreviousBest', false, ...
  'selectedOverBaseline', false);
end

function row = localBuildDoaBasinEntryRow(entrySpec, entrySolve, selectedSolve, adopted, previousBestFval, baselineFval, model)
%LOCALBUILDDOABASINENTRYROW Build one basin-entry diagnostic row.

row = localEmptyDoaBasinEntryRow();
row.tag = string(entrySpec.tag);
row.halfWidth1 = entrySpec.halfWidth(1);
row.halfWidth2 = entrySpec.halfWidth(2);
row.center1 = entrySpec.center(1);
row.center2 = entrySpec.center(2);
row.offset1 = entrySpec.offset(1);
row.offset2 = entrySpec.offset(2);
row.centerSource = string(entrySpec.centerSource);
row.adoptionMode = localResolveDoaBasinEntryAdoptionMode(model);
row.baselineFval = baselineFval;
row.previousBestFval = previousBestFval;
row.entryFval = localGetStructField(entrySolve, 'fval', NaN);
[row.entryLatDeg, row.entryLonDeg, row.entryFdRefHz, row.entryFdRateHzPerSec] = ...
  localExtractDoaBasinOptVarParam(entrySolve, model);
row.entryObjDeltaToPreviousBest = row.entryFval - previousBestFval;
row.entryObjDeltaToBaseline = row.entryFval - baselineFval;
row.selectedFval = localGetStructField(selectedSolve, 'fval', NaN);
[row.selectedLatDeg, row.selectedLonDeg, row.selectedFdRefHz, row.selectedFdRateHzPerSec] = ...
  localExtractDoaBasinOptVarParam(selectedSolve, model);
row.selectedObjDeltaToBaseline = row.selectedFval - baselineFval;
row.selectedVariant = string(localGetStructField(selectedSolve, 'solveVariant', ""));
row.entryExitflag = localGetStructField(entrySolve, 'exitflag', NaN);
row.entryIterations = localGetOutputField(localGetStructField(entrySolve, 'output', struct()), 'iterations', NaN);
row.adoptedOverPreviousBest = logical(adopted);
row.selectedOverBaseline = localIsObjectiveBetter(row.selectedFval, baselineFval);
end

function [latDeg, lonDeg, fdRefHz, fdRateHzPerSec] = localExtractDoaBasinOptVarParam(solveInfo, model)
%LOCALEXTRACTDOABASINOPTVARPARAM Extract user-facing parameters from a solve.

latDeg = NaN;
lonDeg = NaN;
fdRefHz = NaN;
fdRateHzPerSec = NaN;
optVar = reshape(localGetStructField(solveInfo, 'optVar', []), [], 1);
if numel(optVar) < 2
  return;
end
if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
  latDeg = rad2deg(optVar(1));
  lonDeg = rad2deg(optVar(2));
else
  latDeg = optVar(1);
  lonDeg = optVar(2);
end
if numel(optVar) >= 3
  fdRefHz = optVar(3);
end
if numel(optVar) >= 4
  fdRateHzPerSec = optVar(4);
elseif strcmp(localGetStructField(model, 'fdRateMode', ''), 'known')
  fdRateHzPerSec = localGetStructField(model, 'fdRateKnown', NaN);
elseif strcmp(localGetStructField(model, 'fdRateMode', ''), 'zero')
  fdRateHzPerSec = 0;
end
end

function tf = localIsObjectiveBetter(candidateFval, referenceFval)
%LOCALISOBJECTIVEBETTER Compare objective values with a small tolerance.

tf = false;
if ~(isfinite(candidateFval) && isfinite(referenceFval))
  return;
end
objTol = max(1e-6, 1e-9 * max(abs(referenceFval), 1));
tf = candidateFval < referenceFval - objTol;
end

function displayValue = localHalfWidthDisplayValue(model, halfWidth)
%LOCALHALFWIDTHDISPLAYVALUE Return a compact half-width value for tags.

if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
  displayValue = rad2deg(max(abs(halfWidth(:))));
else
  displayValue = max(abs(halfWidth(:)));
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

function value = localGetModelLogical(model, fieldName, defaultValue)
%LOCALGETMODELLOGICAL Read one logical model flag with default fallback.

value = defaultValue;
if ~isstruct(model) || ~isfield(model, fieldName) || isempty(model.(fieldName))
  return;
end
value = logical(model.(fieldName));
end


function tf = localDoaBoxOverlapsActiveEnvelope(model, doaCenter, doaHalfWidth)
%LOCALDOABOXOVERLAPSACTIVEENVELOPE Check child-box overlap before solving.

baseRange = model.doaGrid{1}.range;
if isempty(baseRange) || any(size(baseRange) ~= [2, 2])
  tf = true;
  return;
end
center = reshape(doaCenter, [], 1);
halfWidth = reshape(doaHalfWidth, [], 1);
if numel(center) ~= 2 || numel(halfWidth) ~= 2 || ...
    any(~isfinite(center)) || any(~isfinite(halfWidth)) || any(halfWidth < 0)
  tf = false;
  return;
end
childLb = center - halfWidth;
childUb = center + halfWidth;
if strcmp(localGetStructField(model, 'doaType', ''), 'angle')
  center(1) = mod(center(1), 2 * pi);
  childLb(1) = center(1) - halfWidth(1);
  childUb(1) = center(1) + halfWidth(1);
  if childLb(1) < baseRange(1, 1) || childUb(1) > baseRange(1, 2)
    childLb(1) = baseRange(1, 1);
    childUb(1) = baseRange(1, 2);
  end
end
activeLb = baseRange(:, 1);
activeUb = baseRange(:, 2);
if isfield(model, 'parentDoaLb') && isfield(model, 'parentDoaUb') && ...
    ~isempty(model.parentDoaLb) && ~isempty(model.parentDoaUb)
  parentLb = reshape(model.parentDoaLb, [], 1);
  parentUb = reshape(model.parentDoaUb, [], 1);
  if numel(parentLb) == 2 && numel(parentUb) == 2 && ...
      all(isfinite(parentLb)) && all(isfinite(parentUb))
    activeLb = max(activeLb, parentLb);
    activeUb = min(activeUb, parentUb);
  end
end
finalLb = max(activeLb, childLb);
finalUb = min(activeUb, childUb);
tf = all(finalLb <= finalUb);
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
