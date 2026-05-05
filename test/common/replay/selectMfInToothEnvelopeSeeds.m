function selectedSeedTable = selectMfInToothEnvelopeSeeds(tailDiagnosisTable, opt)
%SELECTMFINTOOTHENVELOPESEEDS Choose representative seeds for heavy in-tooth envelope surfaces.
%
% The selector uses only replay summary tables and offline labels. It is meant
% for replay seed coverage, not runtime estimator gating or candidate adoption.

if nargin < 2 || ~isstruct(opt)
  opt = struct();
end
if isempty(tailDiagnosisTable)
  selectedSeedTable = table();
  return;
end

maxEnvelopeSeedCount = localGetOpt(opt, 'maxEnvelopeSeedCount', 12);
maxHardSeedCount = localGetOpt(opt, 'maxHardSeedCount', 8);
maxGateMissSeedCount = localGetOpt(opt, 'maxGateMissSeedCount', 3);
maxEasyNegativeSeedCount = localGetOpt(opt, 'maxEasyNegativeSeedCount', 2);
maxFdNegativeSeedCount = localGetOpt(opt, 'maxFdNegativeSeedCount', 2);
gateCoherenceThreshold = localGetOpt(opt, 'gatedRescueCoherenceThreshold', 0.20);
truthDoaGapThresholdDeg = localGetOpt(opt, 'truthDoaGapThresholdDeg', 1e-3);

taskSeed = localTableColumn(tailDiagnosisTable, 'taskSeed', nan(height(tailDiagnosisTable), 1));
tailClass = string(localTableColumn(tailDiagnosisTable, 'tailClass', strings(height(tailDiagnosisTable), 1)));
defaultAngleErrDeg = localTableColumn(tailDiagnosisTable, 'multiAngleErrDeg', nan(height(tailDiagnosisTable), 1));
truthDoaAngleErrDeg = localTableColumn(tailDiagnosisTable, 'truthDoaAngleErrDeg', nan(height(tailDiagnosisTable), 1));
defaultCoherence = localTableColumn(tailDiagnosisTable, 'multiNonRefCoherenceFloor', nan(height(tailDiagnosisTable), 1));
gapToTruthDoaDeg = localTableColumn(tailDiagnosisTable, 'multiGapToTruthDoaDeg', defaultAngleErrDeg - truthDoaAngleErrDeg);

isFdNegative = contains(tailClass, "fd not healthy");
isCollapsed = contains(tailClass, "non-ref coherence collapsed");
isDoaBasin = contains(tailClass, "DoA/local-state basin");
isGapTail = isfinite(gapToTruthDoaDeg) & gapToTruthDoaDeg >= truthDoaGapThresholdDeg;
isHard = (isCollapsed | isDoaBasin | isGapTail) & ~isFdNegative;
isGateMiss = isHard & ~(isfinite(defaultCoherence) & defaultCoherence < gateCoherenceThreshold);
isEasyNegative = ~isHard & ~isFdNegative;

selected = false(height(tailDiagnosisTable), 1);
rowCell = cell(0, 1);
[rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, ...
  truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, isHard & ~isGateMiss, ...
  "hard-coherence-collapsed", maxHardSeedCount);
[rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, ...
  truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, isGateMiss, ...
  "hard-gate-miss", maxGateMissSeedCount);
[rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, ...
  truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, isFdNegative, ...
  "fd-not-healthy-negative", maxFdNegativeSeedCount);
[rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, ...
  truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, isEasyNegative, ...
  "easy-negative", maxEasyNegativeSeedCount);

numRemaining = max(0, maxEnvelopeSeedCount - numel(rowCell));
if numRemaining > 0
  fillMask = ~selected;
  [rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, ... %#ok<ASGLU>
    truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, fillMask, ...
    "top-gap-fill", numRemaining);
end

if isempty(rowCell)
  selectedSeedTable = table();
else
  selectedSeedTable = struct2table([rowCell{:}].');
  if height(selectedSeedTable) > maxEnvelopeSeedCount
    selectedSeedTable = selectedSeedTable(1:maxEnvelopeSeedCount, :);
  end
end
end

function [rowCell, selected] = localAppendTopRows(rowCell, selected, taskSeed, tailClass, defaultAngleErrDeg, truthDoaAngleErrDeg, defaultCoherence, gapToTruthDoaDeg, mask, role, maxCount)
%LOCALAPPENDTOPROWS Append the largest-gap rows from one seed class.

if maxCount <= 0 || ~any(mask)
  return;
end
candidateIdx = find(mask(:) & ~selected(:));
if isempty(candidateIdx)
  return;
end
metric = gapToTruthDoaDeg(candidateIdx);
metric(~isfinite(metric)) = defaultAngleErrDeg(candidateIdx(~isfinite(metric)));
metric(~isfinite(metric)) = -inf;
[~, order] = sort(metric, 'descend');
selectIdx = candidateIdx(order(1:min(maxCount, numel(order))));
for iSelect = 1:numel(selectIdx)
  idx = selectIdx(iSelect);
  row = struct();
  row.taskSeed = taskSeed(idx);
  row.selectionRole = string(role);
  row.selectionMetric = gapToTruthDoaDeg(idx);
  row.tailClass = tailClass(idx);
  row.defaultAngleErrDeg = defaultAngleErrDeg(idx);
  row.truthDoaOracleAngleErrDeg = truthDoaAngleErrDeg(idx);
  row.defaultNonRefCoherenceFloor = defaultCoherence(idx);
  rowCell{end + 1, 1} = row; %#ok<AGROW>
  selected(idx) = true;
end
end

function value = localTableColumn(dataTable, fieldName, defaultValue)
%LOCALTABLECOLUMN Read a table column if present, otherwise return a default vector.

if istable(dataTable) && any(strcmp(dataTable.Properties.VariableNames, fieldName))
  value = dataTable.(fieldName);
else
  value = defaultValue;
end
end

function value = localGetOpt(opt, fieldName, defaultValue)
%LOCALGETOPT Read an option field with a scalar fallback.

if isfield(opt, fieldName)
  value = opt.(fieldName);
else
  value = defaultValue;
end
end
