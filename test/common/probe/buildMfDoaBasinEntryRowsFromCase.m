function rowList = buildMfDoaBasinEntryRowsFromCase(method, caseUse, truth, task, probeTag, probeGroup, baselineObj)
%BUILDMFDOABASINENTRYROWSFROMCASE Expose estimator DoA entry diagnostics.

rowList = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
basinDiag = localGetFieldOrDefault(optimInfo, 'doaBasinEntry', struct());
info = summarizeDoaDopplerCase(caseUse, truth);
finalObj = localGetFieldOrDefault(estResult, 'fval', NaN);
if ~isfinite(baselineObj)
  baselineObj = finalObj;
end
selectedVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));

if isempty(basinDiag) || ~isstruct(basinDiag)
  row = emptyMfDoaBasinEntryRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = string(probeTag);
  row.probeGroup = string(probeGroup);
  row.entryReason = "missing-diagnostic";
  row.baselineAngleErrDeg = info.angleErrDeg;
  row.selectedAngleErrDeg = info.angleErrDeg;
  row.message = "optimInfo.doaBasinEntry is missing";
  rowList = row;
  return;
end

entryRows = localGetFieldOrDefault(basinDiag, 'entryTable', repmat(struct(), 0, 1));
entryRows = reshape(entryRows, [], 1);
entryEnabled = logical(localGetFieldOrDefault(basinDiag, 'enabled', false));
entryReason = string(localGetFieldOrDefault(basinDiag, 'reason', ""));
entryScope = string(localGetFieldOrDefault(basinDiag, 'scope', ""));
entryNumSat = localScalarOrDefault(localGetFieldOrDefault(basinDiag, 'numSat', NaN), NaN);
entryBestTag = string(localGetFieldOrDefault(basinDiag, 'bestTag', ""));
entryNumCandidate = localScalarOrDefault(localGetFieldOrDefault(basinDiag, 'numCandidate', numel(entryRows)), numel(entryRows));

if isempty(entryRows)
  row = emptyMfDoaBasinEntryRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = string(probeTag);
  row.probeGroup = string(probeGroup);
  row.entryEnabled = entryEnabled;
  row.entryReason = entryReason;
  row.entryScope = entryScope;
  row.entryNumSat = entryNumSat;
  row.entryNumCandidate = entryNumCandidate;
  row.entryBestTag = entryBestTag;
  row.baselineAngleErrDeg = info.angleErrDeg;
  row.selectedAngleErrDeg = info.angleErrDeg;
  row.selectedVariant = selectedVariant;
  row.selectedFval = finalObj;
  row.selectedObjDeltaToBaseline = finalObj - baselineObj;
  row.selectedFlag = true;
  row.evalOk = entryEnabled;
  rowList = row;
  return;
end

for iEntry = 1:numel(entryRows)
  entry = entryRows(iEntry);
  row = emptyMfDoaBasinEntryRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = string(probeTag);
  row.probeGroup = string(probeGroup);
  row.entryEnabled = entryEnabled;
  row.entryReason = entryReason;
  row.entryScope = entryScope;
  row.entryNumSat = entryNumSat;
  row.entryNumCandidate = entryNumCandidate;
  row.entryBestTag = entryBestTag;
  row.entryTag = string(localGetFieldOrDefault(entry, 'tag', ""));
  row.entryHalfWidthLatDeg = localGetFieldOrDefault(entry, 'halfWidth1', NaN);
  row.entryHalfWidthLonDeg = localGetFieldOrDefault(entry, 'halfWidth2', NaN);
  row.entryFval = localGetFieldOrDefault(entry, 'entryFval', NaN);
  row.entryObjDeltaToBaseline = row.entryFval - baselineObj;
  row.entryIterations = localGetFieldOrDefault(entry, 'entryIterations', NaN);
  row.entryExitflag = localGetFieldOrDefault(entry, 'entryExitflag', NaN);
  row.adoptedOverPreviousBest = logical(localGetFieldOrDefault(entry, 'adoptedOverPreviousBest', false));
  row.polishFval = localGetFieldOrDefault(entry, 'polishFval', NaN);
  row.polishObjDeltaToBaseline = row.polishFval - baselineObj;
  row.polishIterations = localGetFieldOrDefault(entry, 'polishIterations', NaN);
  row.polishExitflag = localGetFieldOrDefault(entry, 'polishExitflag', NaN);
  row.selectedVariant = string(localGetFieldOrDefault(entry, 'selectedVariant', ""));
  if row.selectedVariant == ""
    row.selectedVariant = selectedVariant;
  end
  row.selectedFval = localGetFieldOrDefault(entry, 'selectedFval', NaN);
  row.selectedObjDeltaToBaseline = row.selectedFval - baselineObj;
  row.baselineAngleErrDeg = info.angleErrDeg;
  row.selectedAngleErrDeg = info.angleErrDeg;
  row.selectedFlag = selectedVariant == row.selectedVariant || selectedVariant == row.entryTag || ...
    (strlength(row.entryTag) > 0 && startsWith(selectedVariant, row.entryTag));
  row.evalOk = isfinite(row.entryFval) || isfinite(row.polishFval) || isfinite(row.selectedFval);
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
end

function value = localScalarOrDefault(value, defaultValue)
%LOCALSCALARORDEFAULT Return first scalar value or default.

if isempty(value)
  value = defaultValue;
  return;
end
value = double(value(1));
if ~isfinite(value)
  value = defaultValue;
end
end

function val = localGetFieldOrDefault(s, fieldName, defaultVal)
%LOCALGETFIELDORDEFAULT Get a field with a default fallback.

if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
