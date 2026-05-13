function triggerTable = buildMfBankHealthGateTriggerOutcomeTable(caseTable, rescueCaseTable, rescueDetailTable, varargin)
%BUILDMFBANKHEALTHGATETRIGGEROUTCOMETABLE Summarize health-gated rescue outcomes by trigger family.
%
% This report helper only groups existing replay diagnostic rows. It does not
% run solvers, choose estimator winners, or define a formal rescue strategy.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyTriggerRow(), 0, 1);
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0 || ...
    isempty(rescueDetailTable) || ~istable(rescueDetailTable) || height(rescueDetailTable) == 0
  triggerTable = struct2table(rowList(:));
  return;
end
methodList = reshape(string(opt.MethodList), [], 1);
methodList = methodList(strlength(methodList) > 0);
if isempty(methodList)
  triggerTable = struct2table(rowList(:));
  return;
end
snrList = unique(caseTable.snrDb, 'stable');
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  rescueName = methodName + "-" + string(opt.DiagnosticSuffix);
  for iSnr = 1:numel(snrList)
    snrDb = snrList(iSnr);
    baseRows = localFilterRows(caseTable, methodName, snrDb);
    rescueRows = localFilterRows(rescueCaseTable, rescueName, snrDb);
    detailRows = localFilterRows(rescueDetailTable, rescueName, snrDb);
    if isempty(baseRows) || height(baseRows) == 0 || isempty(detailRows) || height(detailRows) == 0 || ...
        ~ismember('baselineBadGateReason', detailRows.Properties.VariableNames)
      continue;
    end
    triggerGroupList = strings(height(detailRows), 1);
    for iRow = 1:height(detailRows)
      triggerGroupList(iRow) = localTriggerGroup(detailRows.baselineBadGateReason(iRow));
    end
    groupList = unique(triggerGroupList(triggerGroupList ~= "notTriggered"), 'stable');
    for iGroup = 1:numel(groupList)
      groupName = groupList(iGroup);
      groupRows = detailRows(triggerGroupList == groupName, :);
      if isempty(groupRows) || height(groupRows) == 0
        continue;
      end
      row = localBuildGroupRow(methodName, rescueName, snrDb, groupName, groupRows, baseRows, rescueRows, opt);
      rowList(end + 1, 1) = row; %#ok<AGROW>
    end
  end
end
triggerTable = struct2table(rowList(:));
if height(triggerTable) > 0
  triggerTable = sortrows(triggerTable, {'displayName', 'snrDb', 'triggerGroup'});
end
end

function opt = localParseOpt(varargin)
opt = struct('MethodList', string.empty(0, 1), 'DiagnosticSuffix', "HealthGatedBankRescue", ...
  'TrimAngleNormMax', 5, 'SoftDamageAngleNormMargin', 0.1, 'LargeImproveAngleNormMargin', 0.5);
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankHealthGateTriggerOutcomeTable:InvalidNameValue', 'Name-value arguments must be paired.');
end
for iArg = 1:2:numel(varargin)
  name = lower(string(varargin{iArg}));
  value = varargin{iArg + 1};
  switch name
    case "methodlist"
      opt.MethodList = reshape(string(value), [], 1);
    case "diagnosticsuffix"
      opt.DiagnosticSuffix = string(value);
    case "trimanglenormmax"
      opt.TrimAngleNormMax = double(value);
    case "softdamageanglenormmargin"
      opt.SoftDamageAngleNormMargin = double(value);
    case "largeimproveanglenormmargin"
      opt.LargeImproveAngleNormMargin = double(value);
    otherwise
      error('buildMfBankHealthGateTriggerOutcomeTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
end
end

function row = localBuildGroupRow(methodName, rescueName, snrDb, groupName, detailRows, baseRows, rescueRows, opt)
numCase = height(detailRows);
baseAngle = NaN(numCase, 1);
rescueAngle = NaN(numCase, 1);
adoptFlag = false(numCase, 1);
baseTrimKeep = false(numCase, 1);
rescueTrimKeep = false(numCase, 1);
for iRow = 1:numCase
  taskSeed = detailRows.taskSeed(iRow);
  baseRow = localFindSeedRow(baseRows, taskSeed);
  rescueRow = localFindSeedRow(rescueRows, taskSeed);
  baseAngle(iRow) = localTableValue(baseRow, 'angleErrOverSphericalCrb', ...
    localTableValue(detailRows(iRow, :), 'baselineAngleErrOverCrb', NaN));
  rescueAngle(iRow) = localTableValue(rescueRow, 'angleErrOverSphericalCrb', baseAngle(iRow));
  adoptFlag(iRow) = logical(localTableValue(detailRows(iRow, :), 'shadowAdoptFlag', false));
  baseTrimKeep(iRow) = logical(localTableValue(baseRow, 'trimKeep', isfinite(baseAngle(iRow)) && baseAngle(iRow) <= opt.TrimAngleNormMax));
  rescueTrimKeep(iRow) = logical(localTableValue(rescueRow, 'trimKeep', isfinite(rescueAngle(iRow)) && rescueAngle(iRow) <= opt.TrimAngleNormMax));
end
adoptedDelta = rescueAngle - baseAngle;
softDamageMask = adoptFlag & isfinite(adoptedDelta) & adoptedDelta > double(opt.SoftDamageAngleNormMargin);
largeImproveMask = adoptFlag & isfinite(adoptedDelta) & -adoptedDelta >= double(opt.LargeImproveAngleNormMargin);
row = localEmptyTriggerRow();
row.displayName = methodName;
row.rescueDisplayName = rescueName;
row.snrDb = snrDb;
row.triggerGroup = groupName;
row.numCase = numCase;
row.adoptCount = nnz(adoptFlag);
row.adoptRate = localSafeRatio(row.adoptCount, numCase);
row.baselineTrimKeepCount = nnz(baseTrimKeep);
row.baselineTrimKeepRate = localSafeRatio(row.baselineTrimKeepCount, numCase);
row.rescueTrimKeepCount = nnz(rescueTrimKeep);
row.rescueTrimKeepRate = localSafeRatio(row.rescueTrimKeepCount, numCase);
row.softDamageCount = nnz(softDamageMask);
row.softDamageRate = localSafeRatio(row.softDamageCount, max(row.adoptCount, 1));
row.largeImproveCount = nnz(largeImproveMask);
row.largeImproveRate = localSafeRatio(row.largeImproveCount, max(row.adoptCount, 1));
row.medianBaselineAngleOverCrb = median(baseAngle(isfinite(baseAngle)), 'omitnan');
row.medianRescueAngleOverCrb = median(rescueAngle(isfinite(rescueAngle)), 'omitnan');
row.medianAdoptedDeltaOverCrb = median(adoptedDelta(adoptFlag & isfinite(adoptedDelta)), 'omitnan');
end

function triggerGroup = localTriggerGroup(reasonValue)
reason = string(reasonValue);
if strlength(reason) == 0 || reason == "healthy" || reason == "not-gated"
  triggerGroup = "notTriggered";
  return;
end
partList = split(reason, "+");
partList = partList(strlength(partList) > 0);
partSet = unique(partList, 'stable');
if numel(partSet) == 1
  switch partSet(1)
    case "firstOrderOpt"
      triggerGroup = "firstOrderOnly";
    case "iterations"
      triggerGroup = "iterationOnly";
    case {"firstOrderAndIteration", "firstOrderOptAndIteration"}
      triggerGroup = "firstOrderAndIteration";
    case "nonRefCoh"
      triggerGroup = "coherence";
    case "boundary"
      triggerGroup = "boundary";
    case "noSolve"
      triggerGroup = "noSolve";
    case "candidateObjective"
      triggerGroup = "candidateObjective";
    otherwise
      triggerGroup = partSet(1);
  end
  return;
end
if numel(partSet) == 2 && all(ismember(["firstOrderOpt"; "iterations"], partSet))
  triggerGroup = "firstOrderAndIteration";
else
  triggerGroup = "multiTrigger";
end
end

function rows = localFilterRows(tableIn, displayName, snrDb)
if isempty(tableIn) || ~istable(tableIn) || height(tableIn) == 0 || ...
    ~all(ismember({'displayName', 'snrDb'}, tableIn.Properties.VariableNames))
  rows = table();
  return;
end
rows = tableIn(tableIn.displayName == string(displayName) & tableIn.snrDb == snrDb, :);
end

function row = localFindSeedRow(rowTable, taskSeed)
row = table();
if isempty(rowTable) || ~istable(rowTable) || height(rowTable) == 0 || ~ismember('taskSeed', rowTable.Properties.VariableNames)
  return;
end
idx = find(rowTable.taskSeed == taskSeed, 1, 'first');
if ~isempty(idx)
  row = rowTable(idx, :);
end
end

function value = localTableValue(rowTable, varName, defaultValue)
value = defaultValue;
if istable(rowTable) && height(rowTable) > 0 && ismember(varName, rowTable.Properties.VariableNames)
  rawValue = rowTable.(varName);
  if ~isempty(rawValue)
    value = rawValue(1);
  end
end
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function row = localEmptyTriggerRow()
row = struct('displayName', "", 'rescueDisplayName', "", 'snrDb', NaN, 'triggerGroup', "", ...
  'numCase', NaN, 'adoptCount', NaN, 'adoptRate', NaN, ...
  'baselineTrimKeepCount', NaN, 'baselineTrimKeepRate', NaN, ...
  'rescueTrimKeepCount', NaN, 'rescueTrimKeepRate', NaN, ...
  'softDamageCount', NaN, 'softDamageRate', NaN, ...
  'largeImproveCount', NaN, 'largeImproveRate', NaN, ...
  'medianBaselineAngleOverCrb', NaN, 'medianRescueAngleOverCrb', NaN, ...
  'medianAdoptedDeltaOverCrb', NaN);
end
