function [adoptedCaseTable, detailTable] = buildMfBankAdoptedCaseTable(caseTable, candidateTraceTable, shadowTable, varargin)
%BUILDMFBANKADOPTEDCASETABLE Build replay-only adopted metric rows from shadow decisions.
%
% This report helper does not change estimator output. It reuses an offline
% no-truth shadow decision table plus the already computed candidate trace to
% construct a diagnostic case table: baseline rows are kept, and an extra
% adopted row per selected method copies the selected bank candidate metrics
% whenever the shadow decision adopted that candidate.

opt = localParseOpt(varargin{:});
emptyDetail = repmat(localEmptyDetailRow(), 0, 1);
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  adoptedCaseTable = caseTable;
  detailTable = struct2table(emptyDetail(:));
  return;
end
methodList = reshape(string(opt.MethodList), [], 1);
methodList = methodList(strlength(methodList) > 0);
if isempty(methodList)
  adoptedCaseTable = caseTable([], :);
  detailTable = struct2table(emptyDetail(:));
  return;
end
baseMask = ismember(caseTable.displayName, methodList);
baseRows = caseTable(baseMask, :);
if isempty(baseRows) || height(baseRows) == 0
  adoptedCaseTable = caseTable([], :);
  detailTable = struct2table(emptyDetail(:));
  return;
end
adoptRows = baseRows;
detailRows = repmat(localEmptyDetailRow(), height(baseRows), 1);
for iRow = 1:height(baseRows)
  baseName = string(baseRows.displayName(iRow));
  adoptRows.displayName(iRow) = baseName + "-" + string(opt.DiagnosticSuffix);
  if ismember('mfInitMode', adoptRows.Properties.VariableNames)
    adoptRows.mfInitMode(iRow) = string(baseRows.mfInitMode(iRow)) + "+" + lower(string(opt.DiagnosticSuffix));
  end
  detail = localEmptyDetailRow();
  detail.displayName = adoptRows.displayName(iRow);
  detail.baselineDisplayName = baseName;
  detail.snrDb = baseRows.snrDb(iRow);
  detail.taskSeed = baseRows.taskSeed(iRow);
  detail.baselineAngleErrDeg = localTableValue(baseRows(iRow, :), 'finalAngleErrDeg', NaN);
  detail.baselineAngleErrOverCrb = localTableValue(baseRows(iRow, :), 'angleErrOverSphericalCrb', NaN);
  detail.baselineFdRefErrOverCrb = localTableValue(baseRows(iRow, :), 'fdRefErrOverCrb', NaN);
  detail.baselineFirstOrderOpt = localTableValue(baseRows(iRow, :), 'firstOrderOpt', NaN);
  detail.baselineIterations = localTableValue(baseRows(iRow, :), 'iterations', NaN);
  detail.baselineNonRefCoherenceFloor = localTableValue(baseRows(iRow, :), 'nonRefCoherenceFloor', NaN);
  detail.baselineBoundaryHit = logical(localTableValue(baseRows(iRow, :), 'fdRefBoundaryHit', false)) || ...
    logical(localTableValue(baseRows(iRow, :), 'fdRateBoundaryHit', false));
  detail.baselineNoSolveFlag = ~logical(localTableValue(baseRows(iRow, :), 'notNoSolve', true));
  detail.baselineObj = localTableValue(baseRows(iRow, :), 'finalObj', NaN);

  [baselineBadEnough, baselineBadReason] = localIsBaselineBadEnough(detail, opt);
  if ~baselineBadEnough
    detail.baselineBadGateReason = baselineBadReason;
    detail.detailReason = "baseline-kept:baseline-not-bad";
    detailRows(iRow) = detail;
    continue;
  end
  detail.baselineBadGateReason = baselineBadReason;

  shadowRow = localFindShadowRow(shadowTable, baseName, baseRows.snrDb(iRow), baseRows.taskSeed(iRow));
  if isempty(shadowRow)
    detail.detailReason = "no-shadow-row";
    detailRows(iRow) = detail;
    continue;
  end
  detail.shadowAdoptFlag = logical(shadowRow.shadowAdoptFlag(1));
  detail.shadowRejectReason = string(shadowRow.rejectReason(1));
  if ~detail.shadowAdoptFlag
    detail.detailReason = "baseline-kept:" + detail.shadowRejectReason;
    detailRows(iRow) = detail;
    continue;
  end
  candidateRow = localFindShadowSelectedCandidate(candidateTraceTable, shadowRow);
  if isempty(candidateRow)
    detail.detailReason = "selected-candidate-missing";
    detailRows(iRow) = detail;
    continue;
  end

  adoptRows = localApplyCandidateToCaseRow(adoptRows, iRow, candidateRow, opt.BoundaryTolFraction);
  detail.selectedCandidateTag = string(candidateRow.candidateTag(1));
  detail.selectedCandidateSource = string(candidateRow.candidateSource(1));
  detail.selectedCenterSource = localTableValue(candidateRow, 'bankCenterSource', "");
  detail.selectedAxisDirection = localTableValue(candidateRow, 'bankAxisDirection', "");
  detail.selectedStepDeg = localTableValue(candidateRow, 'bankStepDeg', NaN);
  detail.angleImproveDeg = localTableValue(candidateRow, 'angleImproveDeg', NaN);
  detail.angleImproveOverCrb = localTableValue(candidateRow, 'angleImproveOverCrb', NaN);
  detail.fdToothShift = localTableValue(candidateRow, 'candidateToBaselineFdToothShift', NaN);
  detail.selectedFdRefErrOverTooth = localTableValue(candidateRow, 'finalFdRefErrOverTooth', NaN);
  detail.selectedFirstOrderOpt = localTableValue(candidateRow, 'firstOrderOpt', NaN);
  detail.selectedNonRefCoherenceFloor = localTableValue(candidateRow, 'nonRefCoherenceFloor', NaN);
  detail.selectedObjMinusBaseline = localTableValue(candidateRow, 'objMinusBaseline', NaN);
  detail.damageFlag = logical(localTableValue(candidateRow, 'damageFlag', false));
  detail.detailReason = "shadow-adopted";
  detailRows(iRow) = detail;
end
adoptedCaseTable = [baseRows; adoptRows];
detailTable = struct2table(detailRows(:));
end

function opt = localParseOpt(varargin)
opt = struct('MethodList', "MS-MF-CP-U", 'DiagnosticSuffix', "CenterRescue", ...
  'BoundaryTolFraction', 0.02, 'BaselineBadAngleNormMin', -Inf, ...
  'BaselineBadFdRefNormMin', Inf, 'BaselineFirstOrderOptMin', Inf, ...
  'BaselineIterationMin', Inf, 'BaselineNonRefCoherenceMax', -Inf, ...
  'BaselineRequireFirstOrderAndIteration', false, ...
  'BaselineBoundaryGateEnable', false, 'BaselineNoSolveGateEnable', false);
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankAdoptedCaseTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = lower(string(varargin{iArg}));
  value = varargin{iArg + 1};
  switch name
    case "methodlist"
      opt.MethodList = reshape(string(value), [], 1);
    case "diagnosticsuffix"
      opt.DiagnosticSuffix = string(value);
    case "boundarytolfraction"
      opt.BoundaryTolFraction = double(value);
    case "baselinebadanglenormmin"
      opt.BaselineBadAngleNormMin = double(value);
    case "baselinebadfdrefnormmin"
      opt.BaselineBadFdRefNormMin = double(value);
    case "baselinefirstorderoptmin"
      opt.BaselineFirstOrderOptMin = double(value);
    case "baselineiterationmin"
      opt.BaselineIterationMin = double(value);
    case "baselinenonrefcoherencemax"
      opt.BaselineNonRefCoherenceMax = double(value);
    case "baselinerequirefirstorderanditeration"
      opt.BaselineRequireFirstOrderAndIteration = logical(value);
    case "baselineboundarygateenable"
      opt.BaselineBoundaryGateEnable = logical(value);
    case "baselinenosolvegateenable"
      opt.BaselineNoSolveGateEnable = logical(value);
    otherwise
      error('buildMfBankAdoptedCaseTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end


function [tf, reason] = localIsBaselineBadEnough(detail, opt)
%LOCALISBASELINEBADENOUGH Gate adopted diagnostic rows by oracle or health badness.

reasonList = strings(0, 1);
angleGate = double(opt.BaselineBadAngleNormMin);
fdRefGate = double(opt.BaselineBadFdRefNormMin);
angleGateActive = isfinite(angleGate);
fdRefGateActive = isfinite(fdRefGate);

if angleGateActive && isfinite(detail.baselineAngleErrOverCrb) && detail.baselineAngleErrOverCrb > angleGate
  reasonList(end + 1, 1) = "angle"; %#ok<AGROW>
end
if fdRefGateActive && isfinite(detail.baselineFdRefErrOverCrb) && abs(detail.baselineFdRefErrOverCrb) > fdRefGate
  reasonList(end + 1, 1) = "fdRef"; %#ok<AGROW>
end
firstOrderOptMin = double(opt.BaselineFirstOrderOptMin);
firstOrderPass = isfinite(firstOrderOptMin) && isfinite(detail.baselineFirstOrderOpt) && detail.baselineFirstOrderOpt > firstOrderOptMin;
iterationMin = double(opt.BaselineIterationMin);
iterationPass = isfinite(iterationMin) && isfinite(detail.baselineIterations) && detail.baselineIterations >= iterationMin;
if logical(opt.BaselineRequireFirstOrderAndIteration)
  if firstOrderPass && iterationPass
    reasonList(end + 1, 1) = "firstOrderAndIteration"; %#ok<AGROW>
  end
else
  if firstOrderPass
    reasonList(end + 1, 1) = "firstOrderOpt"; %#ok<AGROW>
  end
  if iterationPass
    reasonList(end + 1, 1) = "iterations"; %#ok<AGROW>
  end
end
coherenceMax = double(opt.BaselineNonRefCoherenceMax);
if isfinite(coherenceMax) && isfinite(detail.baselineNonRefCoherenceFloor) && detail.baselineNonRefCoherenceFloor < coherenceMax
  reasonList(end + 1, 1) = "nonRefCoh"; %#ok<AGROW>
end
if logical(opt.BaselineBoundaryGateEnable) && logical(detail.baselineBoundaryHit)
  reasonList(end + 1, 1) = "boundary"; %#ok<AGROW>
end
if logical(opt.BaselineNoSolveGateEnable) && logical(detail.baselineNoSolveFlag)
  reasonList(end + 1, 1) = "noSolve"; %#ok<AGROW>
end

noGateConfigured = ~angleGateActive && ~fdRefGateActive && ...
  ~isfinite(firstOrderOptMin) && ~isfinite(iterationMin) && ~isfinite(coherenceMax) && ...
  ~logical(opt.BaselineBoundaryGateEnable) && ~logical(opt.BaselineNoSolveGateEnable);
if noGateConfigured
  tf = true;
  reason = "not-gated";
  return;
end

tf = ~isempty(reasonList);
if tf
  reason = strjoin(cellstr(reasonList), "+");
else
  reason = "healthy";
end
end

function shadowRow = localFindShadowRow(shadowTable, displayName, snrDb, taskSeed)
shadowRow = table();
if isempty(shadowTable) || ~istable(shadowTable) || height(shadowTable) == 0
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed'};
if ~all(ismember(requiredVars, shadowTable.Properties.VariableNames))
  return;
end
idx = find(shadowTable.displayName == string(displayName) & ...
  shadowTable.snrDb == snrDb & shadowTable.taskSeed == taskSeed, 1, 'first');
if ~isempty(idx)
  shadowRow = shadowTable(idx, :);
end
end

function candidateRow = localFindShadowSelectedCandidate(candidateTraceTable, shadowRow)
candidateRow = table();
if isempty(candidateTraceTable) || ~istable(candidateTraceTable) || height(candidateTraceTable) == 0 || ...
    isempty(shadowRow) || height(shadowRow) == 0
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed', 'candidateSource', 'candidateTag', 'objective'};
if ~all(ismember(requiredVars, candidateTraceTable.Properties.VariableNames)) || ...
    ~all(ismember({'shadowSelectedSource', 'shadowSelectedTag'}, shadowRow.Properties.VariableNames))
  return;
end
mask = candidateTraceTable.displayName == shadowRow.displayName(1) & ...
  candidateTraceTable.snrDb == shadowRow.snrDb(1) & ...
  candidateTraceTable.taskSeed == shadowRow.taskSeed(1) & ...
  candidateTraceTable.candidateSource == shadowRow.shadowSelectedSource(1) & ...
  candidateTraceTable.candidateTag == shadowRow.shadowSelectedTag(1);
rows = candidateTraceTable(mask, :);
if isempty(rows) || height(rows) == 0
  return;
end
if ismember('shadowObjective', shadowRow.Properties.VariableNames) && isfinite(shadowRow.shadowObjective(1))
  objDiff = abs(rows.objective - shadowRow.shadowObjective(1));
  [~, idx] = min(objDiff);
else
  [~, idx] = min(rows.objective);
end
candidateRow = rows(idx, :);
end

function caseRows = localApplyCandidateToCaseRow(caseRows, iRow, candidateRow, boundaryTolFraction)
angleErrDeg = localTableValue(candidateRow, 'finalAngleErrDeg', NaN);
fdRefErrHz = localTableValue(candidateRow, 'finalFdRefErrHz', NaN);
fdRateErrHzPerSec = localTableValue(candidateRow, 'finalFdRateErrHzPerSec', NaN);
if isfinite(angleErrDeg)
  caseRows = localSetIfPresent(caseRows, 'sphericalAngleErrDeg', iRow, angleErrDeg);
  caseRows = localSetIfPresent(caseRows, 'finalAngleErrDeg', iRow, angleErrDeg);
  if ismember('crbSphericalApproxStdDeg', caseRows.Properties.VariableNames)
    caseRows = localSetIfPresent(caseRows, 'angleErrOverSphericalCrb', iRow, ...
      localSafeRatio(angleErrDeg, caseRows.crbSphericalApproxStdDeg(iRow)));
  end
  if ismember('crbTraceStdDeg', caseRows.Properties.VariableNames)
    caseRows = localSetIfPresent(caseRows, 'angleErrOverTraceCrb', iRow, ...
      localSafeRatio(angleErrDeg, caseRows.crbTraceStdDeg(iRow)));
  end
  caseRows = localSetIfPresent(caseRows, 'finalLatDeg', iRow, NaN);
  caseRows = localSetIfPresent(caseRows, 'finalLonDeg', iRow, NaN);
  caseRows = localSetIfPresent(caseRows, 'latErrDeg', iRow, NaN);
  caseRows = localSetIfPresent(caseRows, 'lonErrDeg', iRow, NaN);
  caseRows = localSetIfPresent(caseRows, 'finalMinusInitAngleDeg', iRow, NaN);
end
if isfinite(fdRefErrHz)
  caseRows = localSetIfPresent(caseRows, 'fdRefErrHz', iRow, fdRefErrHz);
  caseRows = localSetIfPresent(caseRows, 'fdRefAbsErrHz', iRow, abs(fdRefErrHz));
  if all(ismember({'fdRefTrueHz', 'finalFdRefHz'}, caseRows.Properties.VariableNames))
    finalFdRefHz = caseRows.fdRefTrueHz(iRow) + fdRefErrHz;
    caseRows.finalFdRefHz(iRow) = finalFdRefHz;
    if ismember('fdRefCrbStdHz', caseRows.Properties.VariableNames)
      caseRows = localSetIfPresent(caseRows, 'fdRefErrOverCrb', iRow, ...
        localSafeRatio(abs(fdRefErrHz), caseRows.fdRefCrbStdHz(iRow)));
    end
    if ismember('initFdRefHz', caseRows.Properties.VariableNames)
      caseRows = localSetIfPresent(caseRows, 'fdRefInitMoveHz', iRow, finalFdRefHz - caseRows.initFdRefHz(iRow));
    end
    if all(ismember({'fdRangeLowHz', 'fdRangeHighHz', 'fdRefBoundaryHit'}, caseRows.Properties.VariableNames))
      caseRows.fdRefBoundaryHit(iRow) = localIsNearRangeBoundary(finalFdRefHz, ...
        [caseRows.fdRangeLowHz(iRow), caseRows.fdRangeHighHz(iRow)], boundaryTolFraction);
    end
  end
end
if isfinite(fdRateErrHzPerSec)
  caseRows = localSetIfPresent(caseRows, 'fdRateErrHzPerSec', iRow, fdRateErrHzPerSec);
  caseRows = localSetIfPresent(caseRows, 'fdRateAbsErrHzPerSec', iRow, abs(fdRateErrHzPerSec));
  if all(ismember({'fdRateTrueHzPerSec', 'finalFdRateHzPerSec'}, caseRows.Properties.VariableNames))
    finalFdRateHzPerSec = caseRows.fdRateTrueHzPerSec(iRow) + fdRateErrHzPerSec;
    caseRows.finalFdRateHzPerSec(iRow) = finalFdRateHzPerSec;
    if ismember('initFdRateHzPerSec', caseRows.Properties.VariableNames)
      caseRows = localSetIfPresent(caseRows, 'fdRateInitMoveHzPerSec', iRow, ...
        finalFdRateHzPerSec - caseRows.initFdRateHzPerSec(iRow));
    end
    if all(ismember({'fdRateRangeLowHzPerSec', 'fdRateRangeHighHzPerSec', 'fdRateBoundaryHit'}, caseRows.Properties.VariableNames))
      caseRows.fdRateBoundaryHit(iRow) = localIsNearRangeBoundary(finalFdRateHzPerSec, ...
        [caseRows.fdRateRangeLowHzPerSec(iRow), caseRows.fdRateRangeHighHzPerSec(iRow)], boundaryTolFraction);
    end
  end
end
caseRows = localSetIfPresent(caseRows, 'finalObj', iRow, localTableValue(candidateRow, 'objective', localTableValue(caseRows(iRow, :), 'finalObj', NaN)));
if all(ismember({'initObj', 'finalObj', 'objectiveImprove'}, caseRows.Properties.VariableNames))
  caseRows.objectiveImprove(iRow) = caseRows.initObj(iRow) - caseRows.finalObj(iRow);
end
caseRows = localSetIfPresent(caseRows, 'iterations', iRow, localTableValue(candidateRow, 'iterations', localTableValue(caseRows(iRow, :), 'iterations', NaN)));
caseRows = localSetIfPresent(caseRows, 'firstOrderOpt', iRow, localTableValue(candidateRow, 'firstOrderOpt', localTableValue(caseRows(iRow, :), 'firstOrderOpt', NaN)));
caseRows = localSetIfPresent(caseRows, 'nonRefCoherenceFloor', iRow, localTableValue(candidateRow, 'nonRefCoherenceFloor', localTableValue(caseRows(iRow, :), 'nonRefCoherenceFloor', NaN)));
caseRows = localSetIfPresent(caseRows, 'candidateCount', iRow, localTableValue(candidateRow, 'candidateCount', localTableValue(caseRows(iRow, :), 'candidateCount', NaN)));
caseRows = localSetIfPresent(caseRows, 'solveVariant', iRow, localTableValue(candidateRow, 'solveVariant', localTableValue(caseRows(iRow, :), 'solveVariant', "")));
caseRows = localSetIfPresent(caseRows, 'candidateVariantText', iRow, string(localTableValue(candidateRow, 'candidateTag', "")));
end

function tableOut = localSetIfPresent(tableOut, varName, iRow, value)
if istable(tableOut) && ismember(varName, tableOut.Properties.VariableNames)
  tableOut.(varName)(iRow) = value;
end
end

function value = localTableValue(rowTable, varName, defaultValue)
value = defaultValue;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames)
  values = rowTable.(varName);
  if ~isempty(values)
    value = values(1);
  end
end
end

function tf = localIsNearRangeBoundary(value, rangeValue, boundaryTolFraction)
tf = false;
rangeValue = reshape(double(rangeValue), [], 1);
if numel(rangeValue) < 2 || ~isfinite(value) || any(~isfinite(rangeValue(1:2)))
  return;
end
width = abs(rangeValue(2) - rangeValue(1));
if width <= 0
  return;
end
tol = max(double(boundaryTolFraction), 0) .* width;
tf = abs(value - rangeValue(1)) <= tol || abs(value - rangeValue(2)) <= tol;
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function row = localEmptyDetailRow()
row = struct('displayName', "", 'baselineDisplayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'shadowAdoptFlag', false, 'shadowRejectReason', "", 'baselineBadGateReason', "", ...
  'selectedCandidateTag', "", 'selectedCandidateSource', "", ...
  'selectedCenterSource', "", 'selectedAxisDirection', "", 'selectedStepDeg', NaN, ...
  'baselineAngleErrDeg', NaN, 'baselineAngleErrOverCrb', NaN, ...
  'angleImproveDeg', NaN, 'angleImproveOverCrb', NaN, ...
  'baselineFdRefErrOverCrb', NaN, 'baselineFirstOrderOpt', NaN, ...
  'baselineIterations', NaN, 'baselineNonRefCoherenceFloor', NaN, ...
  'baselineBoundaryHit', false, 'baselineNoSolveFlag', false, ...
  'selectedFdRefErrOverTooth', NaN, 'fdToothShift', NaN, ...
  'baselineObj', NaN, 'selectedObjMinusBaseline', NaN, 'selectedFirstOrderOpt', NaN, ...
  'selectedNonRefCoherenceFloor', NaN, 'damageFlag', false, 'detailReason', "");
end
