function shadowTable = buildMfBankAdoptionShadowTable(candidateTraceTable, varargin)
%BUILDMFBANKADOPTIONSHADOWTABLE Evaluate no-truth bank adoption guards offline.
%
% The helper selects at most one candidate per method/SNR/seed using only
% objective, first-order-opt, and coherence guard quantities that are already
% present in the candidate trace. Truth-based angle columns are only used to
% evaluate the shadow decision after the candidate has been selected.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyShadowRow(), 0, 1);
if isempty(candidateTraceTable) || ~istable(candidateTraceTable) || height(candidateTraceTable) == 0
  shadowTable = struct2table(rowList(:));
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed', 'candidateSource', 'candidateTag', ...
  'objective', 'objMinusBaseline', 'finalAngleErrDeg', 'baselineAngleErrDeg', ...
  'angleImproveDeg', 'angleImproveOverCrb', 'damageFlag'};
if ~all(ismember(requiredVars, candidateTraceTable.Properties.VariableNames))
  shadowTable = struct2table(rowList(:));
  return;
end
keyTable = unique(candidateTraceTable(:, {'displayName', 'snrDb', 'taskSeed'}), 'rows', 'stable');
for iKey = 1:height(keyTable)
  keyMask = candidateTraceTable.displayName == keyTable.displayName(iKey) & ...
    candidateTraceTable.snrDb == keyTable.snrDb(iKey) & ...
    candidateTraceTable.taskSeed == keyTable.taskSeed(iKey);
  rows = candidateTraceTable(keyMask, :);
  baselineRows = rows(rows.candidateSource == "baseline" | rows.candidateTag == "baseline", :);
  if isempty(baselineRows) || height(baselineRows) == 0
    continue;
  end
  baseline = baselineRows(1, :);
  candidateRows = rows(ismember(rows.candidateSource, opt.CandidateSourceList), :);
  row = localEmptyShadowRow();
  row.displayName = keyTable.displayName(iKey);
  row.snrDb = keyTable.snrDb(iKey);
  row.taskSeed = keyTable.taskSeed(iKey);
  row.primaryDiagnosis = localLookupPrimaryDiagnosis(row.displayName, row.snrDb, row.taskSeed, opt.PrimaryDiagnosisTable);
  row.baselineAngleErrDeg = baseline.finalAngleErrDeg(1);
  row.baselineObjective = baseline.objective(1);
  row.baselineNonRefCoherenceFloor = localTableValue(baseline, 'nonRefCoherenceFloor', NaN);
  row.baselineFirstOrderOpt = localTableValue(baseline, 'firstOrderOpt', NaN);
  row.numCandidate = height(candidateRows);
  if isempty(candidateRows) || height(candidateRows) == 0
    row.rejectReason = "no-candidate";
    rowList(end + 1, 1) = row; %#ok<AGROW>
    continue;
  end
  objTol = max(double(opt.ObjectiveAbsMargin), ...
    double(opt.ObjectiveRelMargin) * max(abs(row.baselineObjective), 1));
  eligibleMask = isfinite(candidateRows.objective) & isfinite(candidateRows.objMinusBaseline) & ...
    candidateRows.objMinusBaseline < -objTol;
  if ismember('firstOrderOpt', candidateRows.Properties.VariableNames) && isfinite(opt.MaxFirstOrderOpt)
    firstOrder = candidateRows.firstOrderOpt;
    eligibleMask = eligibleMask & (~isfinite(firstOrder) | firstOrder <= opt.MaxFirstOrderOpt);
  end
  if ismember('nonRefCoherenceFloor', candidateRows.Properties.VariableNames) && ...
      isfinite(row.baselineNonRefCoherenceFloor) && isfinite(opt.MaxCoherenceDrop)
    coherenceFloor = candidateRows.nonRefCoherenceFloor;
    eligibleMask = eligibleMask & (~isfinite(coherenceFloor) | ...
      coherenceFloor >= row.baselineNonRefCoherenceFloor - opt.MaxCoherenceDrop);
  end
  row.numEligible = nnz(eligibleMask);
  if ~any(eligibleMask)
    row.rejectReason = "guard-rejected";
    rowList(end + 1, 1) = row; %#ok<AGROW>
    continue;
  end
  eligibleIdx = find(eligibleMask);
  [~, bestLocalIdx] = min(candidateRows.objective(eligibleMask));
  selected = candidateRows(eligibleIdx(bestLocalIdx), :);
  row.shadowAdoptFlag = true;
  row.shadowSelectedTag = selected.candidateTag(1);
  row.shadowSelectedSource = selected.candidateSource(1);
  row.shadowBankCenterSource = localTableValue(selected, 'bankCenterSource', "");
  row.shadowBankAxisDirection = localTableValue(selected, 'bankAxisDirection', "");
  row.shadowBankStepDeg = localTableValue(selected, 'bankStepDeg', NaN);
  row.shadowFdRefErrOverTooth = localTableValue(selected, 'finalFdRefErrOverTooth', NaN);
  row.shadowFdToothShift = localTableValue(selected, 'candidateToBaselineFdToothShift', NaN);
  row.shadowObjective = selected.objective(1);
  row.shadowObjMinusBaseline = selected.objMinusBaseline(1);
  row.shadowFirstOrderOpt = localTableValue(selected, 'firstOrderOpt', NaN);
  row.shadowNonRefCoherenceFloor = localTableValue(selected, 'nonRefCoherenceFloor', NaN);
  row.shadowAngleErrDeg = selected.finalAngleErrDeg(1);
  row.shadowAngleImproveDeg = selected.angleImproveDeg(1);
  row.shadowAngleImproveOverCrb = selected.angleImproveOverCrb(1);
  row.shadowObjectiveBeatFlag = true;
  row.shadowAngleImproveFlag = isfinite(row.shadowAngleImproveDeg) && row.shadowAngleImproveDeg > 0;
  row.shadowDamageFlag = logical(selected.damageFlag(1));
  row.rejectReason = "adopt";
  rowList(end + 1, 1) = row; %#ok<AGROW>
end
shadowTable = struct2table(rowList(:));
end

function opt = localParseOpt(varargin)
opt = struct('CandidateSourceList', "ms-bank-cpu-release", ...
  'PrimaryDiagnosisTable', table(), ...
  'ObjectiveRelMargin', 1e-9, 'ObjectiveAbsMargin', 1e-6, ...
  'MaxFirstOrderOpt', Inf, 'MaxCoherenceDrop', Inf);
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankAdoptionShadowTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = lower(string(varargin{iArg}));
  value = varargin{iArg + 1};
  switch name
    case "candidatesourcelist"
      opt.CandidateSourceList = reshape(string(value), 1, []);
    case "primarydiagnosistable"
      opt.PrimaryDiagnosisTable = value;
    case "objectiverelmargin"
      opt.ObjectiveRelMargin = double(value);
    case "objectiveabsmargin"
      opt.ObjectiveAbsMargin = double(value);
    case "maxfirstorderopt"
      opt.MaxFirstOrderOpt = double(value);
    case "maxcoherencedrop"
      opt.MaxCoherenceDrop = double(value);
    otherwise
      error('buildMfBankAdoptionShadowTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function primaryDiagnosis = localLookupPrimaryDiagnosis(displayName, snrDb, taskSeed, primaryDiagnosisTable)
primaryDiagnosis = "all";
if isempty(primaryDiagnosisTable) || ~istable(primaryDiagnosisTable) || height(primaryDiagnosisTable) == 0
  return;
end
requiredVars = {'displayName', 'snrDb', 'taskSeed', 'primaryDiagnosis'};
if ~all(ismember(requiredVars, primaryDiagnosisTable.Properties.VariableNames))
  return;
end
idx = find(primaryDiagnosisTable.displayName == string(displayName) & ...
  primaryDiagnosisTable.snrDb == snrDb & primaryDiagnosisTable.taskSeed == taskSeed, 1, 'first');
if ~isempty(idx)
  primaryDiagnosis = primaryDiagnosisTable.primaryDiagnosis(idx);
end
end

function value = localTableValue(rowTable, varName, defaultValue)
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames)
  value = rowTable.(varName)(1);
else
  value = defaultValue;
end
end

function row = localEmptyShadowRow()
row = struct('displayName', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'primaryDiagnosis', "", 'numCandidate', NaN, 'numEligible', NaN, ...
  'baselineAngleErrDeg', NaN, 'baselineObjective', NaN, ...
  'baselineFirstOrderOpt', NaN, 'baselineNonRefCoherenceFloor', NaN, ...
  'shadowAdoptFlag', false, 'shadowSelectedTag', "", 'shadowSelectedSource', "", ...
  'shadowBankCenterSource', "", 'shadowBankAxisDirection', "", 'shadowBankStepDeg', NaN, ...
  'shadowFdRefErrOverTooth', NaN, 'shadowFdToothShift', NaN, ...
  'shadowObjective', NaN, 'shadowObjMinusBaseline', NaN, ...
  'shadowFirstOrderOpt', NaN, 'shadowNonRefCoherenceFloor', NaN, ...
  'shadowAngleErrDeg', NaN, 'shadowAngleImproveDeg', NaN, ...
  'shadowAngleImproveOverCrb', NaN, 'shadowObjectiveBeatFlag', false, ...
  'shadowAngleImproveFlag', false, 'shadowDamageFlag', false, 'rejectReason', "");
end
