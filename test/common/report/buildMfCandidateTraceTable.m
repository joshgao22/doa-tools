function candidateTraceTable = buildMfCandidateTraceTable(caseTable, objectiveProbeTable, solveProbeTable, doaBasinEntryTable, varargin)
%BUILDMFCANDIDATETRACETABLE Merge baseline, probe, objective, and entry candidates.
%
% This report helper is intentionally metric-only. It does not select a
% winner, change estimator output, or encode method-specific adoption rules.

opt = localParseOpt(varargin{:});
if nargin < 4
  doaBasinEntryTable = table();
end
if nargin < 3 || isempty(solveProbeTable)
  solveProbeTable = table();
end
if nargin < 2 || isempty(objectiveProbeTable)
  objectiveProbeTable = table();
end

rowList = repmat(localEmptyCandidateTraceRow(), 0, 1);
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  candidateTraceTable = struct2table(rowList(:));
  return;
end

methodPrefix = string(opt.MethodPrefix);
if strlength(methodPrefix) > 0
  caseMask = startsWith(caseTable.displayName, methodPrefix);
else
  caseMask = true(height(caseTable), 1);
end
caseRows = caseTable(caseMask, :);
for iCase = 1:height(caseRows)
  baseRow = caseRows(iCase, :);
  methodName = baseRow.displayName(1);
  snrDb = baseRow.snrDb(1);
  taskSeed = baseRow.taskSeed(1);
  solveRows = localFilterRows(solveProbeTable, methodName, taskSeed, snrDb);
  objRows = localFilterRows(objectiveProbeTable, methodName, taskSeed, snrDb);
  basinRows = localFilterRows(doaBasinEntryTable, methodName, taskSeed, snrDb);
  rowList = [rowList; localBuildSolveCandidateTraceRows(baseRow, solveRows)]; %#ok<AGROW>
  rowList = [rowList; localBuildBasinEntryCandidateTraceRows(baseRow, basinRows)]; %#ok<AGROW>
  rowList = [rowList; localBuildObjectiveCandidateTraceRows(baseRow, objRows)]; %#ok<AGROW>
end
candidateTraceTable = struct2table(rowList(:));
end

function opt = localParseOpt(varargin)
opt = struct('MethodPrefix', "");
if mod(numel(varargin), 2) ~= 0
  error('buildMfCandidateTraceTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = string(varargin{iArg});
  value = varargin{iArg + 1};
  switch lower(name)
    case "methodprefix"
      opt.MethodPrefix = string(value);
    otherwise
      error('buildMfCandidateTraceTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function rows = localFilterRows(rowTable, methodName, taskSeed, snrDb)
rows = table();
if isempty(rowTable) || ~istable(rowTable) || height(rowTable) == 0 || ...
    ~all(ismember({'displayName', 'taskSeed', 'snrDb'}, rowTable.Properties.VariableNames))
  return;
end
mask = rowTable.displayName == string(methodName) & ...
  rowTable.taskSeed == taskSeed & rowTable.snrDb == snrDb;
rows = rowTable(mask, :);
end

function rowList = localBuildSolveCandidateTraceRows(baseRow, solveRows)
rowList = repmat(localEmptyCandidateTraceRow(), 0, 1);
if isempty(solveRows) || ~istable(solveRows) || height(solveRows) == 0
  return;
end
baselineAngleErrDeg = baseRow.finalAngleErrDeg(1);
angleCrbStdDeg = baseRow.crbSphericalApproxStdDeg(1);
baselineObj = baseRow.finalObj(1);
for iRow = 1:height(solveRows)
  row = localEmptyCandidateTraceRow();
  row.displayName = solveRows.displayName(iRow);
  row.snrDb = solveRows.snrDb(iRow);
  row.taskSeed = solveRows.taskSeed(iRow);
  row.candidateTag = solveRows.probeTag(iRow);
  row.candidateSource = localCandidateSourceFromTag(row.candidateTag, "solve");
  row.initDoaSource = localCandidateDoaSourceFromTag(row.candidateTag);
  row.initFdSource = localCandidateFdSourceFromTag(row.candidateTag);
  row.doaHalfWidthDeg = max([solveRows.doaHalfWidthLatDeg(iRow), solveRows.doaHalfWidthLonDeg(iRow)], [], 'omitnan');
  row.fdRefHalfWidthHz = baseRow.fdRangeHalfWidthHz(1);
  row.fdRateHalfWidthHzPerSec = 0.5 * abs(baseRow.fdRateRangeHighHzPerSec(1) - baseRow.fdRateRangeLowHzPerSec(1));
  row.solveVariant = solveRows.solveVariant(iRow);
  row.candidateCount = solveRows.candidateCount(iRow);
  row.iterations = solveRows.iterations(iRow);
  if ismember('firstOrderOpt', solveRows.Properties.VariableNames)
    row.firstOrderOpt = solveRows.firstOrderOpt(iRow);
  end
  if ismember('nonRefCoherenceFloor', solveRows.Properties.VariableNames)
    row.nonRefCoherenceFloor = solveRows.nonRefCoherenceFloor(iRow);
  end
  if row.candidateTag == "baseline"
    row.firstOrderOpt = baseRow.firstOrderOpt(1);
    row.nonRefCoherenceFloor = baseRow.nonRefCoherenceFloor(1);
  end
  row.baselineAngleErrDeg = baselineAngleErrDeg;
  row.finalAngleErrDeg = solveRows.finalAngleErrDeg(iRow);
  row.finalAngleErrOverCrb = localSafeRatio(row.finalAngleErrDeg, angleCrbStdDeg);
  row.angleImproveDeg = baselineAngleErrDeg - row.finalAngleErrDeg;
  row.angleImproveOverCrb = localSafeRatio(row.angleImproveDeg, angleCrbStdDeg);
  row.finalFdRefErrHz = solveRows.finalFdRefErrHz(iRow);
  row.finalFdRateErrHzPerSec = solveRows.finalFdRateErrHzPerSec(iRow);
  row.objective = solveRows.finalObj(iRow);
  row.objMinusBaseline = solveRows.objMinusBaseline(iRow);
  if row.candidateTag == "baseline"
    row.objMinusBaseline = 0;
    row.objective = baselineObj;
  end
  row.selectedFlag = row.candidateTag == "baseline";
  row.wouldBeatBaselineFlag = isfinite(row.objMinusBaseline) && row.objMinusBaseline < -localObjectiveTol(baselineObj);
  row.damageFlag = ((isfinite(row.objMinusBaseline) && row.objMinusBaseline > localObjectiveTol(baselineObj)) || ...
    (isfinite(row.angleImproveDeg) && row.angleImproveDeg < -1e-6)) && ~row.wouldBeatBaselineFlag;
  row = localAnnotateBankAndToothDiagnostics(row, baseRow);
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
end

function rowList = localBuildBasinEntryCandidateTraceRows(baseRow, basinRows)
rowList = repmat(localEmptyCandidateTraceRow(), 0, 1);
if isempty(basinRows) || ~istable(basinRows) || height(basinRows) == 0
  return;
end
baselineAngleErrDeg = baseRow.finalAngleErrDeg(1);
angleCrbStdDeg = baseRow.crbSphericalApproxStdDeg(1);
baselineObj = baseRow.finalObj(1);
for iRow = 1:height(basinRows)
  if ~logical(basinRows.entryEnabled(iRow)) || strlength(string(basinRows.entryTag(iRow))) == 0
    continue;
  end
  row = localEmptyCandidateTraceRow();
  row.displayName = basinRows.displayName(iRow);
  row.snrDb = basinRows.snrDb(iRow);
  row.taskSeed = basinRows.taskSeed(iRow);
  row.candidateTag = "estimator-" + string(basinRows.probeTag(iRow)) + ":" + string(basinRows.entryTag(iRow));
  row.candidateSource = "estimator-doa-basin-entry";
  row.initDoaSource = "estimator-entry";
  row.initFdSource = localBasinEntryFdSourceFromProbeTag(basinRows.probeTag(iRow));
  row.doaHalfWidthDeg = max([basinRows.entryHalfWidthLatDeg(iRow), basinRows.entryHalfWidthLonDeg(iRow)], [], 'omitnan');
  row.fdRefHalfWidthHz = baseRow.fdRangeHalfWidthHz(1);
  row.fdRateHalfWidthHzPerSec = 0.5 * abs(baseRow.fdRateRangeHighHzPerSec(1) - baseRow.fdRateRangeLowHzPerSec(1));
  row.solveVariant = basinRows.selectedVariant(iRow);
  row.candidateCount = basinRows.entryNumCandidate(iRow);
  row.iterations = basinRows.entryIterations(iRow);
  row.baselineAngleErrDeg = baselineAngleErrDeg;
  if logical(basinRows.selectedFlag(iRow))
    row.finalAngleErrDeg = basinRows.selectedAngleErrDeg(iRow);
    row.finalAngleErrOverCrb = localSafeRatio(row.finalAngleErrDeg, angleCrbStdDeg);
    row.angleImproveDeg = baselineAngleErrDeg - row.finalAngleErrDeg;
    row.angleImproveOverCrb = localSafeRatio(row.angleImproveDeg, angleCrbStdDeg);
  end
  row.objective = basinRows.entryFval(iRow);
  row.objMinusBaseline = basinRows.entryObjDeltaToBaseline(iRow);
  row.selectedFlag = logical(basinRows.selectedFlag(iRow));
  row.wouldBeatBaselineFlag = isfinite(row.objMinusBaseline) && row.objMinusBaseline < -localObjectiveTol(baselineObj);
  row.damageFlag = isfinite(row.objMinusBaseline) && row.objMinusBaseline > localObjectiveTol(baselineObj);
  row = localAnnotateBankAndToothDiagnostics(row, baseRow);
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
end

function rowList = localBuildObjectiveCandidateTraceRows(baseRow, objRows)
rowList = repmat(localEmptyCandidateTraceRow(), 0, 1);
if isempty(objRows) || ~istable(objRows) || height(objRows) == 0
  return;
end
baselineAngleErrDeg = baseRow.finalAngleErrDeg(1);
angleCrbStdDeg = baseRow.crbSphericalApproxStdDeg(1);
baselineObj = localProbeObj(objRows, "finalPoint");
if ~isfinite(baselineObj)
  baselineObj = baseRow.finalObj(1);
end
for iRow = 1:height(objRows)
  row = localEmptyCandidateTraceRow();
  row.displayName = objRows.displayName(iRow);
  row.snrDb = objRows.snrDb(iRow);
  row.taskSeed = objRows.taskSeed(iRow);
  row.candidateTag = objRows.probeTag(iRow);
  row.candidateSource = localCandidateSourceFromTag(row.candidateTag, "objective");
  row.initDoaSource = localCandidateDoaSourceFromTag(row.candidateTag);
  row.initFdSource = localCandidateFdSourceFromTag(row.candidateTag);
  row.doaHalfWidthDeg = NaN;
  row.fdRefHalfWidthHz = baseRow.fdRangeHalfWidthHz(1);
  row.fdRateHalfWidthHzPerSec = 0.5 * abs(baseRow.fdRateRangeHighHzPerSec(1) - baseRow.fdRateRangeLowHzPerSec(1));
  row.solveVariant = "objective-eval";
  row.candidateCount = 0;
  row.iterations = 0;
  row.baselineAngleErrDeg = baselineAngleErrDeg;
  row.finalAngleErrDeg = objRows.angleErrDeg(iRow);
  row.finalAngleErrOverCrb = localSafeRatio(row.finalAngleErrDeg, angleCrbStdDeg);
  row.angleImproveDeg = baselineAngleErrDeg - row.finalAngleErrDeg;
  row.angleImproveOverCrb = localSafeRatio(row.angleImproveDeg, angleCrbStdDeg);
  row.finalFdRefErrHz = objRows.fdRefErrHz(iRow);
  row.finalFdRateErrHzPerSec = objRows.fdRateErrHzPerSec(iRow);
  row.objective = objRows.objective(iRow);
  row.objMinusBaseline = objRows.objMinusFinal(iRow);
  row.selectedFlag = row.candidateTag == "finalPoint";
  row.wouldBeatBaselineFlag = isfinite(row.objMinusBaseline) && row.objMinusBaseline < -localObjectiveTol(baselineObj);
  row.damageFlag = ((isfinite(row.objMinusBaseline) && row.objMinusBaseline > localObjectiveTol(baselineObj)) || ...
    (isfinite(row.angleImproveDeg) && row.angleImproveDeg < -1e-6)) && ~row.wouldBeatBaselineFlag;
  row = localAnnotateBankAndToothDiagnostics(row, baseRow);
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
end

function row = localAnnotateBankAndToothDiagnostics(row, baseRow)
%LOCALANNOTATEBANKANDTOOTHDIAGNOSTICS Add source-family and tooth diagnostics.

[centerSource, axisDirection, stepDeg] = localParseBankTag(row.candidateTag);
row.bankCenterSource = centerSource;
row.bankAxisDirection = axisDirection;
row.bankStepDeg = stepDeg;
row.toothStepHz = localTableScalar(baseRow, 'toothStepHz', NaN);
baselineFdRefErrHz = localTableScalar(baseRow, 'fdRefErrHz', NaN);
if isfinite(row.finalFdRefErrHz) && isfinite(row.toothStepHz) && abs(row.toothStepHz) > eps
  row.finalFdRefErrOverTooth = row.finalFdRefErrHz ./ row.toothStepHz;
  row.fdRefNearestToothIdx = round(row.finalFdRefErrHz ./ row.toothStepHz);
  row.fdRefDistanceToNearestToothHz = row.finalFdRefErrHz - row.fdRefNearestToothIdx .* row.toothStepHz;
end
if isfinite(row.finalFdRefErrHz) && isfinite(baselineFdRefErrHz) && ...
    isfinite(row.toothStepHz) && abs(row.toothStepHz) > eps
  row.candidateToBaselineFdToothShift = round((row.finalFdRefErrHz - baselineFdRefErrHz) ./ row.toothStepHz);
end
end

function [centerSource, axisDirection, stepDeg] = localParseBankTag(tagText)
%LOCALPARSEBANKTAG Parse replay-only MS bank tag into source, axis, and step.

centerSource = "";
axisDirection = "";
stepDeg = NaN;
tagText = string(tagText);
if startsWith(tagText, "ms-bank-cpk-")
  tail = extractAfter(tagText, strlength("ms-bank-cpk-"));
elseif startsWith(tagText, "ms-bank-cpu-release-")
  tail = extractAfter(tagText, strlength("ms-bank-cpu-release-"));
else
  return;
end
partList = split(tail, "-");
if numel(partList) < 2
  return;
end
centerSource = partList(1);
axisToken = partList(2);
if axisToken == "center"
  axisDirection = "center";
  stepDeg = 0;
  return;
end
expr = '^(latp|latm|lonp|lonm)(.*)$';
tok = regexp(char(axisToken), expr, 'tokens', 'once');
if isempty(tok)
  axisDirection = axisToken;
  return;
end
axisDirection = string(tok{1});
stepToken = string(tok{2});
stepToken = replace(stepToken, "p", ".");
stepToken = replace(stepToken, "m", "-");
stepDeg = str2double(stepToken);
end

function value = localTableScalar(rowTable, varName, defaultValue)
%LOCALTABLESCALAR Read one scalar value from a one-row table.

value = defaultValue;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames)
  valueUse = rowTable.(varName);
  if ~isempty(valueUse)
    value = valueUse(1);
  end
end
end

function source = localCandidateSourceFromTag(tagText, sourceKind)
tagText = string(tagText);
sourceKind = string(sourceKind);
if tagText == "baseline" || tagText == "finalPoint"
  source = "baseline";
elseif startsWith(tagText, "truth")
  source = "truth-oracle-" + sourceKind;
elseif startsWith(tagText, "static-doa-wide")
  source = "static-wide-" + sourceKind;
elseif startsWith(tagText, "ss-mf-seed")
  source = "ss-mf-seed-" + sourceKind;
elseif startsWith(tagText, "ms-bank-cpk")
  source = "ms-bank-cpk-preselect";
elseif startsWith(tagText, "ms-bank-cpu")
  source = "ms-bank-cpu-release";
elseif tagText == "initPoint" || tagText == "staticSeedPoint" || tagText == "truthDoaFinalFdPoint" || tagText == "finalDoaTruthFdPoint"
  source = "objective-point";
else
  source = sourceKind + "-probe";
end
end

function doaSource = localCandidateDoaSourceFromTag(tagText)
tagText = string(tagText);
if startsWith(tagText, "truth") || contains(tagText, "truthDoa", 'IgnoreCase', true)
  doaSource = "truth-oracle";
elseif startsWith(tagText, "static-doa") || tagText == "staticSeedPoint" || startsWith(tagText, "ms-bank-cpk-static") || startsWith(tagText, "ms-bank-cpu-static")
  doaSource = "static-seed";
elseif startsWith(tagText, "ss-mf-seed") || startsWith(tagText, "ms-bank-cpk-ssmf") || startsWith(tagText, "ms-bank-cpu-ssmf")
  doaSource = "ss-mf-seed";
elseif tagText == "initPoint"
  doaSource = "init";
elseif tagText == "baseline" || tagText == "finalPoint" || contains(tagText, "FinalDoa")
  doaSource = "final";
else
  doaSource = "unknown";
end
end

function fdSource = localCandidateFdSourceFromTag(tagText)
tagText = string(tagText);
if contains(tagText, "TruthFd") || startsWith(tagText, "truth") || startsWith(tagText, "static-doa-wide") || startsWith(tagText, "ss-mf-seed")
  fdSource = "truth-oracle";
elseif startsWith(tagText, "ms-bank-cpk")
  fdSource = "truth-fd-known-rate";
elseif startsWith(tagText, "ms-bank-cpu")
  fdSource = "cp-k-preselect";
elseif contains(tagText, "FinalFd") || tagText == "baseline" || tagText == "finalPoint"
  fdSource = "final";
elseif tagText == "initPoint" || tagText == "staticSeedPoint"
  fdSource = "init";
else
  fdSource = "unknown";
end
end

function fdSource = localBasinEntryFdSourceFromProbeTag(probeTag)
probeTag = string(probeTag);
if probeTag == "baseline"
  fdSource = "model-init";
elseif startsWith(probeTag, "truth") || startsWith(probeTag, "static-doa-wide") || startsWith(probeTag, "ss-mf-seed")
  fdSource = "truth-oracle";
else
  fdSource = "unknown";
end
end

function obj = localProbeObj(objRows, probeTag)
obj = NaN;
if isempty(objRows) || ~istable(objRows) || height(objRows) == 0 || ~ismember('probeTag', objRows.Properties.VariableNames)
  return;
end
idx = find(objRows.probeTag == string(probeTag), 1, 'first');
if ~isempty(idx) && ismember('objective', objRows.Properties.VariableNames)
  obj = objRows.objective(idx);
end
end

function objTol = localObjectiveTol(objValue)
objTol = max(1e-6, 1e-9 * max(abs(objValue), 1));
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function row = localEmptyCandidateTraceRow()
row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'candidateTag', "", 'candidateSource', "", 'initDoaSource', "", 'initFdSource', "", ...
  'doaHalfWidthDeg', NaN, 'fdRefHalfWidthHz', NaN, 'fdRateHalfWidthHzPerSec', NaN, ...
  'solveVariant', "", 'candidateCount', NaN, 'iterations', NaN, 'firstOrderOpt', NaN, ...
  'baselineAngleErrDeg', NaN, 'finalAngleErrDeg', NaN, 'finalAngleErrOverCrb', NaN, ...
  'angleImproveDeg', NaN, 'angleImproveOverCrb', NaN, ...
  'finalFdRefErrHz', NaN, 'finalFdRateErrHzPerSec', NaN, ...
  'toothStepHz', NaN, 'finalFdRefErrOverTooth', NaN, ...
  'fdRefNearestToothIdx', NaN, 'fdRefDistanceToNearestToothHz', NaN, ...
  'candidateToBaselineFdToothShift', NaN, ...
  'bankCenterSource', "", 'bankAxisDirection', "", 'bankStepDeg', NaN, ...
  'nonRefCoherenceFloor', NaN, 'objective', NaN, 'objMinusBaseline', NaN, ...
  'selectedFlag', false, 'wouldBeatBaselineFlag', false, 'damageFlag', false);
end
