function outcomeTable = buildMfBankRescueOutcomeTable(caseTable, rescueCaseTable, rescueDetailTable, shadowTable, varargin)
%BUILDMFBANKRESCUEOUTCOMETABLE Summarize replay-only bank rescue rates.
%
% This report helper only compares already-computed baseline rows, offline
% shadow decisions, and diagnostic adopted rows. It does not run solvers,
% define estimator winner adoption, or change numerical paths.

opt = localParseOpt(varargin{:});
rowList = repmat(localEmptyOutcomeRow(), 0, 1);
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  outcomeTable = struct2table(rowList(:));
  return;
end
methodList = reshape(string(opt.MethodList), [], 1);
methodList = methodList(strlength(methodList) > 0);
if isempty(methodList)
  outcomeTable = struct2table(rowList(:));
  return;
end
snrList = unique(caseTable.snrDb, 'stable');
for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  rescueName = methodName + "-" + string(opt.DiagnosticSuffix);
  for iSnr = 1:numel(snrList)
    snrDb = snrList(iSnr);
    baseRows = caseTable(caseTable.displayName == methodName & caseTable.snrDb == snrDb, :);
    if isempty(baseRows) || height(baseRows) == 0
      continue;
    end
    rescueRows = localFilterRows(rescueCaseTable, rescueName, snrDb);
    detailRows = localFilterRows(rescueDetailTable, rescueName, snrDb);
    shadowRows = localFilterRows(shadowTable, methodName, snrDb);

    row = localEmptyOutcomeRow();
    row.displayName = methodName;
    row.rescueDisplayName = rescueName;
    row.snrDb = snrDb;
    row.numCase = height(baseRows);
    row.baselineBadCount = nnz(localBadMask(baseRows, opt.TrimAngleNormMax));
    row.baselineBadRate = localSafeRatio(row.baselineBadCount, row.numCase);
    row.baselineHealthKeepCount = nnz(localLogicalColumn(baseRows, 'healthResolved', true(height(baseRows), 1)));
    row.baselineHealthKeepRate = localSafeRatio(row.baselineHealthKeepCount, row.numCase);
    row.baselineTrimKeepCount = nnz(localLogicalColumn(baseRows, 'trimKeep', true(height(baseRows), 1)));
    row.baselineTrimKeepRate = localSafeRatio(row.baselineTrimKeepCount, row.numCase);
    row.candidateCaseCount = localCountCandidateCases(shadowRows);
    row.candidateCaseRate = localSafeRatio(row.candidateCaseCount, row.numCase);
    row.adoptCount = nnz(localLogicalColumn(detailRows, 'shadowAdoptFlag', false(height(detailRows), 1)));
    row.adoptRate = localSafeRatio(row.adoptCount, row.numCase);
    row.adoptImproveCount = localCountDetailImprove(detailRows);
    row.adoptImproveRate = localSafeRatio(row.adoptImproveCount, row.numCase);
    row.damageCount = localCountDetailFlag(detailRows, 'damageFlag');
    row.damageRate = localSafeRatio(row.damageCount, row.numCase);
    row.sameToothAdoptCount = localCountSameTooth(detailRows);
    row.sameToothAdoptRate = localSafeRatio(row.sameToothAdoptCount, max(row.adoptCount, 1));
    row.rescueBadCount = nnz(localBadMask(rescueRows, opt.TrimAngleNormMax));
    row.rescueBadRate = localSafeRatio(row.rescueBadCount, row.numCase);
    row.rescueTrimKeepCount = nnz(localLogicalColumn(rescueRows, 'trimKeep', true(height(rescueRows), 1)));
    row.rescueTrimKeepRate = localSafeRatio(row.rescueTrimKeepCount, row.numCase);
    row.medianBaselineAngleOverCrb = localMedianIfPresent(baseRows, 'angleErrOverSphericalCrb');
    row.medianRescueAngleOverCrb = localMedianIfPresent(rescueRows, 'angleErrOverSphericalCrb');
    row.medianBaselineAngleDeg = localMedianIfPresent(baseRows, 'sphericalAngleErrDeg');
    row.medianRescueAngleDeg = localMedianIfPresent(rescueRows, 'sphericalAngleErrDeg');
    adoptedDetailRows = detailRows(localLogicalColumn(detailRows, 'shadowAdoptFlag', false(height(detailRows), 1)), :);
    row.medianAdoptAngleImproveDeg = localMedianIfPresent(adoptedDetailRows, 'angleImproveDeg');
    row.medianFdToothShift = localMedianIfPresent(adoptedDetailRows, 'fdToothShift');
    row = localAnnotateBaselineRescueOutcome(row, baseRows, rescueRows, detailRows, opt);
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
outcomeTable = struct2table(rowList(:));
if height(outcomeTable) > 0
  outcomeTable = sortrows(outcomeTable, {'displayName', 'snrDb'});
end
end

function opt = localParseOpt(varargin)
opt = struct('MethodList', string.empty(0, 1), 'DiagnosticSuffix', "BankRescue", ...
  'TrimAngleNormMax', 5, 'BadGateAngleNormMin', 1.5, ...
  'BadGateFdRefNormMin', Inf, 'SoftDamageAngleNormMargin', 0.1, ...
  'LargeImproveAngleNormMargin', 0.5, 'HealthGateFirstOrderOptMin', Inf, ...
  'HealthGateIterationMin', Inf, 'HealthGateNonRefCoherenceMax', -Inf, ...
  'HealthGateBoundaryEnable', false, 'HealthGateNoSolveEnable', false);
if mod(numel(varargin), 2) ~= 0
  error('buildMfBankRescueOutcomeTable:InvalidNameValue', 'Name-value arguments must be paired.');
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
    case "badgateanglenormmin"
      opt.BadGateAngleNormMin = double(value);
    case "badgatefdrefnormmin"
      opt.BadGateFdRefNormMin = double(value);
    case "softdamageanglenormmargin"
      opt.SoftDamageAngleNormMargin = double(value);
    case "largeimproveanglenormmargin"
      opt.LargeImproveAngleNormMargin = double(value);
    case "healthgatefirstorderoptmin"
      opt.HealthGateFirstOrderOptMin = double(value);
    case "healthgateiterationmin"
      opt.HealthGateIterationMin = double(value);
    case "healthgatenonrefcoherencemax"
      opt.HealthGateNonRefCoherenceMax = double(value);
    case "healthgateboundaryenable"
      opt.HealthGateBoundaryEnable = logical(value);
    case "healthgatenosolveenable"
      opt.HealthGateNoSolveEnable = logical(value);
    otherwise
      error('buildMfBankRescueOutcomeTable:UnknownOption', 'Unknown option "%s".', char(name));
  end
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

function mask = localBadMask(rowTable, trimAngleNormMax)
mask = false(height(rowTable), 1);
if isempty(rowTable) || ~istable(rowTable) || height(rowTable) == 0
  return;
end
if ismember('trimKeep', rowTable.Properties.VariableNames)
  mask = ~rowTable.trimKeep;
  return;
end
if ismember('angleErrOverSphericalCrb', rowTable.Properties.VariableNames)
  mask = rowTable.angleErrOverSphericalCrb > trimAngleNormMax;
end
end

function values = localLogicalColumn(rowTable, varName, defaultValues)
values = defaultValues;
if isempty(rowTable) || ~istable(rowTable) || height(rowTable) == 0
  values = false(0, 1);
  return;
end
if ismember(varName, rowTable.Properties.VariableNames)
  values = logical(rowTable.(varName));
else
  values = defaultValues;
end
values = reshape(values, [], 1);
end

function count = localCountCandidateCases(shadowRows)
count = 0;
if isempty(shadowRows) || ~istable(shadowRows) || height(shadowRows) == 0 || ...
    ~ismember('numCandidate', shadowRows.Properties.VariableNames)
  return;
end
count = nnz(shadowRows.numCandidate > 0);
end

function count = localCountShadowFlag(shadowRows, varName)
count = 0;
if isempty(shadowRows) || ~istable(shadowRows) || height(shadowRows) == 0 || ...
    ~all(ismember({'shadowAdoptFlag', varName}, shadowRows.Properties.VariableNames))
  return;
end
count = nnz(shadowRows.shadowAdoptFlag & logical(shadowRows.(varName)));
end

function count = localCountDetailFlag(detailRows, varName)
count = 0;
if isempty(detailRows) || ~istable(detailRows) || height(detailRows) == 0 || ...
    ~all(ismember({'shadowAdoptFlag', varName}, detailRows.Properties.VariableNames))
  return;
end
count = nnz(detailRows.shadowAdoptFlag & logical(detailRows.(varName)));
end

function count = localCountDetailImprove(detailRows)
count = 0;
if isempty(detailRows) || ~istable(detailRows) || height(detailRows) == 0 || ...
    ~all(ismember({'shadowAdoptFlag', 'angleImproveDeg'}, detailRows.Properties.VariableNames))
  return;
end
count = nnz(detailRows.shadowAdoptFlag & isfinite(detailRows.angleImproveDeg) & detailRows.angleImproveDeg > 0);
end

function count = localCountSameTooth(detailRows)
count = 0;
if isempty(detailRows) || ~istable(detailRows) || height(detailRows) == 0 || ...
    ~all(ismember({'shadowAdoptFlag', 'fdToothShift'}, detailRows.Properties.VariableNames))
  return;
end
count = nnz(detailRows.shadowAdoptFlag & detailRows.fdToothShift == 0);
end

function value = localMedianIfPresent(rowTable, varName)
value = NaN;
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames) && height(rowTable) > 0
  value = median(rowTable.(varName), 'omitnan');
end
end

function ratio = localSafeRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end


function row = localAnnotateBaselineRescueOutcome(row, baseRows, rescueRows, detailRows, opt)
%LOCALANNOTATEBASELINERESCUEOUTCOME Add CRB-normalized soft-damage rates.

numCase = max(row.numCase, 0);
baseAngle = localColumnOrNaN(baseRows, 'angleErrOverSphericalCrb');
baseFdRef = abs(localColumnOrNaN(baseRows, 'fdRefErrOverCrb'));
rescueAngle = NaN(size(baseAngle));
rescueFdRef = NaN(size(baseFdRef));
adoptFlag = false(size(baseAngle));
for iBase = 1:height(baseRows)
  taskSeed = baseRows.taskSeed(iBase);
  rescueRow = localFindSeedRow(rescueRows, taskSeed);
  detailRow = localFindSeedRow(detailRows, taskSeed);
  rescueAngle(iBase) = localTableValue(rescueRow, 'angleErrOverSphericalCrb', baseAngle(iBase));
  rescueFdRef(iBase) = abs(localTableValue(rescueRow, 'fdRefErrOverCrb', baseFdRef(iBase)));
  adoptFlag(iBase) = logical(localTableValue(detailRow, 'shadowAdoptFlag', false));
end
angleGate = double(opt.BadGateAngleNormMin);
fdRefGate = double(opt.BadGateFdRefNormMin);
angleBadGateMask = isfinite(angleGate) & isfinite(baseAngle) & baseAngle > angleGate;
fdRefBadGateMask = isfinite(fdRefGate) & isfinite(baseFdRef) & baseFdRef > fdRefGate;
badGateMask = angleBadGateMask | fdRefBadGateMask;
firstOrderHealthMask = localNumericGateMask(baseRows, 'firstOrderOpt', opt.HealthGateFirstOrderOptMin, 'greater');
iterationHealthMask = localNumericGateMask(baseRows, 'iterations', opt.HealthGateIterationMin, 'greaterEqual');
coherenceHealthMask = localNumericGateMask(baseRows, 'nonRefCoherenceFloor', opt.HealthGateNonRefCoherenceMax, 'less');
boundaryHealthMask = false(size(baseAngle));
if logical(opt.HealthGateBoundaryEnable)
  boundaryHealthMask = localLogicalColumn(baseRows, 'fdRefBoundaryHit', false(size(baseAngle))) | ...
    localLogicalColumn(baseRows, 'fdRateBoundaryHit', false(size(baseAngle)));
end
noSolveHealthMask = false(size(baseAngle));
if logical(opt.HealthGateNoSolveEnable)
  noSolveHealthMask = ~localLogicalColumn(baseRows, 'notNoSolve', true(size(baseAngle)));
end
healthGateMask = firstOrderHealthMask | iterationHealthMask | coherenceHealthMask | ...
  boundaryHealthMask | noSolveHealthMask;
healthyMask = isfinite(baseAngle) & ~badGateMask;
adoptedDelta = rescueAngle - baseAngle;
softDamageMask = adoptFlag & isfinite(adoptedDelta) & adoptedDelta > double(opt.SoftDamageAngleNormMargin);
largeImproveMask = adoptFlag & isfinite(adoptedDelta) & -adoptedDelta >= double(opt.LargeImproveAngleNormMargin);
angleBadRescueMask = angleBadGateMask & isfinite(rescueAngle) & rescueAngle <= angleGate;
fdRefBadRescueMask = fdRefBadGateMask & isfinite(rescueFdRef) & rescueFdRef <= fdRefGate;
badGateRescueMask = badGateMask & (~angleBadGateMask | angleBadRescueMask) & (~fdRefBadGateMask | fdRefBadRescueMask);

row.angleBadGateCaseCount = nnz(angleBadGateMask);
row.angleBadGateCaseRate = localSafeRatio(row.angleBadGateCaseCount, numCase);
row.fdRefBadGateCaseCount = nnz(fdRefBadGateMask);
row.fdRefBadGateCaseRate = localSafeRatio(row.fdRefBadGateCaseCount, numCase);
row.badGateCaseCount = nnz(badGateMask);
row.badGateCaseRate = localSafeRatio(row.badGateCaseCount, numCase);
row.healthGateCaseCount = nnz(healthGateMask);
row.healthGateCaseRate = localSafeRatio(row.healthGateCaseCount, numCase);
row.firstOrderHealthGateCaseCount = nnz(firstOrderHealthMask);
row.firstOrderHealthGateCaseRate = localSafeRatio(row.firstOrderHealthGateCaseCount, numCase);
row.iterationHealthGateCaseCount = nnz(iterationHealthMask);
row.iterationHealthGateCaseRate = localSafeRatio(row.iterationHealthGateCaseCount, numCase);
row.coherenceHealthGateCaseCount = nnz(coherenceHealthMask);
row.coherenceHealthGateCaseRate = localSafeRatio(row.coherenceHealthGateCaseCount, numCase);
row.boundaryHealthGateCaseCount = nnz(boundaryHealthMask);
row.boundaryHealthGateCaseRate = localSafeRatio(row.boundaryHealthGateCaseCount, numCase);
row.noSolveHealthGateCaseCount = nnz(noSolveHealthMask);
row.noSolveHealthGateCaseRate = localSafeRatio(row.noSolveHealthGateCaseCount, numCase);
row.healthGateAdoptCount = nnz(healthGateMask & adoptFlag);
row.healthGateAdoptRate = localSafeRatio(row.healthGateAdoptCount, max(row.healthGateCaseCount, 1));
row.healthyCaseCount = nnz(healthyMask);
row.healthyCaseRate = localSafeRatio(row.healthyCaseCount, numCase);
row.healthyAdoptCount = nnz(healthyMask & adoptFlag);
row.healthyAdoptRate = localSafeRatio(row.healthyAdoptCount, max(row.healthyCaseCount, 1));
row.badAdoptCount = nnz(badGateMask & adoptFlag);
row.badAdoptRate = localSafeRatio(row.badAdoptCount, max(row.badGateCaseCount, 1));
row.angleBadRescueCount = nnz(angleBadRescueMask);
row.angleBadRescueRate = localSafeRatio(row.angleBadRescueCount, max(row.angleBadGateCaseCount, 1));
row.fdRefBadRescueCount = nnz(fdRefBadRescueMask);
row.fdRefBadRescueRate = localSafeRatio(row.fdRefBadRescueCount, max(row.fdRefBadGateCaseCount, 1));
row.badGateRescueCount = nnz(badGateRescueMask);
row.badGateRescueRate = localSafeRatio(row.badGateRescueCount, max(row.badGateCaseCount, 1));
row.healthGateRescueCount = nnz(healthGateMask & isfinite(rescueAngle) & rescueAngle <= double(opt.TrimAngleNormMax));
row.healthGateRescueRate = localSafeRatio(row.healthGateRescueCount, max(row.healthGateCaseCount, 1));
row.softDamageCount = nnz(softDamageMask);
row.softDamageRate = localSafeRatio(row.softDamageCount, max(row.adoptCount, 1));
row.largeImproveCount = nnz(largeImproveMask);
row.largeImproveRate = localSafeRatio(row.largeImproveCount, max(row.adoptCount, 1));
row.medianAdoptedDeltaOverCrb = median(adoptedDelta(adoptFlag & isfinite(adoptedDelta)), 'omitnan');
end

function mask = localNumericGateMask(rowTable, varName, threshold, relation)
%LOCALNUMERICGATEMASK Build a finite-threshold mask from an existing numeric column.

mask = false(height(rowTable), 1);
if ~isfinite(threshold) || isempty(rowTable) || ~istable(rowTable) || ...
    height(rowTable) == 0 || ~ismember(varName, rowTable.Properties.VariableNames)
  return;
end
values = double(rowTable.(varName));
switch string(relation)
  case "greater"
    mask = isfinite(values) & values > threshold;
  case "greaterEqual"
    mask = isfinite(values) & values >= threshold;
  case "less"
    mask = isfinite(values) & values < threshold;
  otherwise
    error('buildMfBankRescueOutcomeTable:UnknownRelation', 'Unknown numeric gate relation "%s".', char(string(relation)));
end
mask = reshape(mask, [], 1);
end

function values = localColumnOrNaN(rowTable, varName)
values = NaN(height(rowTable), 1);
if istable(rowTable) && ismember(varName, rowTable.Properties.VariableNames)
  values = rowTable.(varName);
  values = reshape(double(values), [], 1);
end
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

function row = localEmptyOutcomeRow()
row = struct('displayName', "", 'rescueDisplayName', "", 'snrDb', NaN, 'numCase', NaN, ...
  'baselineBadCount', NaN, 'baselineBadRate', NaN, ...
  'baselineHealthKeepCount', NaN, 'baselineHealthKeepRate', NaN, ...
  'baselineTrimKeepCount', NaN, 'baselineTrimKeepRate', NaN, ...
  'candidateCaseCount', NaN, 'candidateCaseRate', NaN, ...
  'adoptCount', NaN, 'adoptRate', NaN, ...
  'adoptImproveCount', NaN, 'adoptImproveRate', NaN, ...
  'damageCount', NaN, 'damageRate', NaN, ...
  'sameToothAdoptCount', NaN, 'sameToothAdoptRate', NaN, ...
  'rescueBadCount', NaN, 'rescueBadRate', NaN, ...
  'rescueTrimKeepCount', NaN, 'rescueTrimKeepRate', NaN, ...
  'medianBaselineAngleOverCrb', NaN, 'medianRescueAngleOverCrb', NaN, ...
  'medianBaselineAngleDeg', NaN, 'medianRescueAngleDeg', NaN, ...
  'medianAdoptAngleImproveDeg', NaN, 'medianFdToothShift', NaN, ...
  'angleBadGateCaseCount', NaN, 'angleBadGateCaseRate', NaN, ...
  'fdRefBadGateCaseCount', NaN, 'fdRefBadGateCaseRate', NaN, ...
  'badGateCaseCount', NaN, 'badGateCaseRate', NaN, ...
  'healthGateCaseCount', NaN, 'healthGateCaseRate', NaN, ...
  'firstOrderHealthGateCaseCount', NaN, 'firstOrderHealthGateCaseRate', NaN, ...
  'iterationHealthGateCaseCount', NaN, 'iterationHealthGateCaseRate', NaN, ...
  'coherenceHealthGateCaseCount', NaN, 'coherenceHealthGateCaseRate', NaN, ...
  'boundaryHealthGateCaseCount', NaN, 'boundaryHealthGateCaseRate', NaN, ...
  'noSolveHealthGateCaseCount', NaN, 'noSolveHealthGateCaseRate', NaN, ...
  'healthGateAdoptCount', NaN, 'healthGateAdoptRate', NaN, ...
  'healthGateRescueCount', NaN, 'healthGateRescueRate', NaN, ...
  'healthyCaseCount', NaN, 'healthyCaseRate', NaN, ...
  'healthyAdoptCount', NaN, 'healthyAdoptRate', NaN, ...
  'badAdoptCount', NaN, 'badAdoptRate', NaN, ...
  'angleBadRescueCount', NaN, 'angleBadRescueRate', NaN, ...
  'fdRefBadRescueCount', NaN, 'fdRefBadRescueRate', NaN, ...
  'badGateRescueCount', NaN, 'badGateRescueRate', NaN, ...
  'softDamageCount', NaN, 'softDamageRate', NaN, ...
  'largeImproveCount', NaN, 'largeImproveRate', NaN, ...
  'medianAdoptedDeltaOverCrb', NaN);
end
