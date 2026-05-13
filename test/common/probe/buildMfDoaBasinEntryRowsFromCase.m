function rowList = buildMfDoaBasinEntryRowsFromCase(method, caseUse, truth, task, probeTag, probeGroup, baselineObj)
%BUILDMFDOABASINENTRYROWSFROMCASE Expose estimator DoA entry diagnostics.

rowList = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
basinDiag = localGetFieldOrDefault(optimInfo, 'doaBasinEntry', struct());
info = summarizeDoaDopplerCase(caseUse, truth);
finalObj = localGetFieldOrDefault(estResult, 'fval', NaN);
if ~isfinite(baselineObj)
  baselineObj = localGetFieldOrDefault(basinDiag, 'baselineFval', NaN);
end
if ~isfinite(baselineObj)
  baselineObj = finalObj;
end
baselineVariant = string(localGetFieldOrDefault(basinDiag, 'baselineVariant', ""));
selectedVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));

if isempty(basinDiag) || ~isstruct(basinDiag)
  row = emptyMfDoaBasinEntryRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = string(probeTag);
  row.probeGroup = string(probeGroup);
  row.entryReason = "missing-diagnostic";
  row.baselineFval = finalObj;
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
entryAdoptionMode = string(localGetFieldOrDefault(basinDiag, 'adoptionMode', ""));
entryNumSat = localScalarOrDefault(localGetFieldOrDefault(basinDiag, 'numSat', NaN), NaN);
entryBestTag = string(localGetFieldOrDefault(basinDiag, 'bestTag', ""));
basinSelectedVariant = string(localGetFieldOrDefault(basinDiag, 'selectedVariant', ""));
if strlength(basinSelectedVariant) == 0
  basinSelectedVariant = entryBestTag;
end
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
  row.entryAdoptionMode = entryAdoptionMode;
  row.entryNumSat = entryNumSat;
  row.entryNumCandidate = entryNumCandidate;
  row.entryBestTag = entryBestTag;
  row.baselineFval = baselineObj;
  row.baselineVariant = baselineVariant;
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
  row.entryAdoptionMode = entryAdoptionMode;
  row.entryNumSat = entryNumSat;
  row.entryNumCandidate = entryNumCandidate;
  row.entryBestTag = entryBestTag;
  row.baselineFval = localGetFieldOrDefault(entry, 'baselineFval', baselineObj);
  if ~isfinite(row.baselineFval)
    row.baselineFval = baselineObj;
  end
  row.baselineVariant = baselineVariant;
  row.entryTag = string(localGetFieldOrDefault(entry, 'tag', ""));
  row.entryHalfWidthLatDeg = localGetFieldOrDefault(entry, 'halfWidth1', NaN);
  row.entryHalfWidthLonDeg = localGetFieldOrDefault(entry, 'halfWidth2', NaN);
  row.entryCenterLatDeg = localGetFieldOrDefault(entry, 'center1', NaN);
  row.entryCenterLonDeg = localGetFieldOrDefault(entry, 'center2', NaN);
  row.entryOffsetLatDeg = localGetFieldOrDefault(entry, 'offset1', NaN);
  row.entryOffsetLonDeg = localGetFieldOrDefault(entry, 'offset2', NaN);
  row.entryCenterSource = string(localGetFieldOrDefault(entry, 'centerSource', ""));
  row.previousBestFval = localGetFieldOrDefault(entry, 'previousBestFval', NaN);
  row.entryFval = localGetFieldOrDefault(entry, 'entryFval', NaN);
  row.entryLatDeg = localGetFieldOrDefault(entry, 'entryLatDeg', NaN);
  row.entryLonDeg = localGetFieldOrDefault(entry, 'entryLonDeg', NaN);
  row.entryFdRefHz = localGetFieldOrDefault(entry, 'entryFdRefHz', NaN);
  row.entryFdRateHzPerSec = localGetFieldOrDefault(entry, 'entryFdRateHzPerSec', NaN);
  row.entryObjDeltaToPreviousBest = localGetFieldOrDefault(entry, 'entryObjDeltaToPreviousBest', row.entryFval - row.previousBestFval);
  row.entryObjDeltaToBaseline = localGetFieldOrDefault(entry, 'entryObjDeltaToBaseline', row.entryFval - row.baselineFval);
  row.entryIterations = localGetFieldOrDefault(entry, 'entryIterations', NaN);
  row.entryExitflag = localGetFieldOrDefault(entry, 'entryExitflag', NaN);
  row.adoptedOverPreviousBest = logical(localGetFieldOrDefault(entry, 'adoptedOverPreviousBest', false));
  row.polishFval = localGetFieldOrDefault(entry, 'polishFval', NaN);
  row.polishLatDeg = localGetFieldOrDefault(entry, 'polishLatDeg', NaN);
  row.polishLonDeg = localGetFieldOrDefault(entry, 'polishLonDeg', NaN);
  row.polishFdRefHz = localGetFieldOrDefault(entry, 'polishFdRefHz', NaN);
  row.polishFdRateHzPerSec = localGetFieldOrDefault(entry, 'polishFdRateHzPerSec', NaN);
  row.polishObjDeltaToBaseline = localGetFieldOrDefault(entry, 'polishObjDeltaToBaseline', row.polishFval - row.baselineFval);
  row.polishIterations = localGetFieldOrDefault(entry, 'polishIterations', NaN);
  row.polishExitflag = localGetFieldOrDefault(entry, 'polishExitflag', NaN);
  row.selectedVariant = string(localGetFieldOrDefault(entry, 'selectedVariant', ""));
  if strlength(row.selectedVariant) == 0
    row.selectedVariant = basinSelectedVariant;
  end
  if strlength(row.selectedVariant) == 0
    row.selectedVariant = selectedVariant;
  end
  row.selectedFval = localGetFieldOrDefault(entry, 'selectedFval', NaN);
  row.selectedLatDeg = localGetFieldOrDefault(entry, 'selectedLatDeg', NaN);
  row.selectedLonDeg = localGetFieldOrDefault(entry, 'selectedLonDeg', NaN);
  row.selectedFdRefHz = localGetFieldOrDefault(entry, 'selectedFdRefHz', NaN);
  row.selectedFdRateHzPerSec = localGetFieldOrDefault(entry, 'selectedFdRateHzPerSec', NaN);
  row.selectedObjDeltaToBaseline = localGetFieldOrDefault(entry, 'selectedObjDeltaToBaseline', row.selectedFval - row.baselineFval);
  row.selectedOverBaseline = logical(localGetFieldOrDefault(entry, 'selectedOverBaseline', row.selectedFval < row.baselineFval));
  row.baselineAngleErrDeg = info.angleErrDeg;
  row.selectedAngleErrDeg = info.angleErrDeg;
  row.selectedFlag = localIsSelectedBasinEntry(row.entryTag, row.selectedVariant, ...
    basinSelectedVariant, entryBestTag, selectedVariant);
  row.evalOk = isfinite(row.entryFval) || isfinite(row.polishFval) || isfinite(row.selectedFval);
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
end

function tf = localIsSelectedBasinEntry(entryTag, rowSelectedVariant, basinSelectedVariant, entryBestTag, finalSolveVariant)
%LOCALISSELECTEDBASINENTRY Match one row to the basin-entry winner.
% The final CP-U warm-anchor may replace the solveVariant after DoA entry.
% For family attribution, use the basin-entry selected variant / best tag first
% and fall back to the final solve variant only when the basin diagnostic is
% incomplete.

entryTag = string(entryTag);
candidateList = string.empty(0, 1);
if strlength(string(basinSelectedVariant)) > 0
  candidateList(end + 1, 1) = string(basinSelectedVariant);
end
if strlength(string(entryBestTag)) > 0
  candidateList(end + 1, 1) = string(entryBestTag);
end
if isempty(candidateList) && strlength(string(rowSelectedVariant)) > 0
  candidateList(end + 1, 1) = string(rowSelectedVariant);
end
if isempty(candidateList) && strlength(string(finalSolveVariant)) > 0
  candidateList(end + 1, 1) = string(finalSolveVariant);
end

tf = false;
for iCand = 1:numel(candidateList)
  tagUse = candidateList(iCand);
  if strlength(tagUse) == 0
    continue;
  end
  if tagUse == entryTag || (strlength(entryTag) > 0 && startsWith(tagUse, entryTag))
    tf = true;
    return;
  end
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
