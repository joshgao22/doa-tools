function seedChainTable = buildMfSeedChainTable(caseTable)
%BUILDMFSEEDCHAINTABLE Compare static seeds, SS-MF, and MS-MF per seed.

seedChainTable = table();
if isempty(caseTable) || height(caseTable) == 0
  return;
end
keyTable = unique(caseTable(:, {'taskSeed', 'snrDb'}), 'rows', 'stable');
numRow = height(keyTable);
taskSeed = keyTable.taskSeed;
snrDb = keyTable.snrDb;
ssSfAngleDeg = NaN(numRow, 1);
msSfAngleDeg = NaN(numRow, 1);
ssMfKAngleDeg = NaN(numRow, 1);
ssMfUAngleDeg = NaN(numRow, 1);
msMfKAngleDeg = NaN(numRow, 1);
msMfUAngleDeg = NaN(numRow, 1);
msSfMinusSsSfDeg = NaN(numRow, 1);
msMfKMinusMsSfDeg = NaN(numRow, 1);
msMfUMinusMsSfDeg = NaN(numRow, 1);
nonRefCoherenceFloorK = NaN(numRow, 1);
nonRefCoherenceFloorU = NaN(numRow, 1);
for iRow = 1:numRow
  ssSfAngleDeg(iRow) = localCaseValue(caseTable, "SS-SF-Static", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  msSfAngleDeg(iRow) = localCaseValue(caseTable, "MS-SF-Static", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  ssMfKAngleDeg(iRow) = localCaseValue(caseTable, "SS-MF-CP-K", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  ssMfUAngleDeg(iRow) = localCaseValue(caseTable, "SS-MF-CP-U", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  msMfKAngleDeg(iRow) = localCaseValue(caseTable, "MS-MF-CP-K", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  msMfUAngleDeg(iRow) = localCaseValue(caseTable, "MS-MF-CP-U", taskSeed(iRow), snrDb(iRow), 'sphericalAngleErrDeg');
  msSfMinusSsSfDeg(iRow) = msSfAngleDeg(iRow) - ssSfAngleDeg(iRow);
  msMfKMinusMsSfDeg(iRow) = msMfKAngleDeg(iRow) - msSfAngleDeg(iRow);
  msMfUMinusMsSfDeg(iRow) = msMfUAngleDeg(iRow) - msSfAngleDeg(iRow);
  nonRefCoherenceFloorK(iRow) = localCaseValue(caseTable, "MS-MF-CP-K", taskSeed(iRow), snrDb(iRow), 'nonRefCoherenceFloor');
  nonRefCoherenceFloorU(iRow) = localCaseValue(caseTable, "MS-MF-CP-U", taskSeed(iRow), snrDb(iRow), 'nonRefCoherenceFloor');
end
seedChainTable = table(taskSeed, snrDb, ssSfAngleDeg, msSfAngleDeg, ...
  ssMfKAngleDeg, ssMfUAngleDeg, msMfKAngleDeg, msMfUAngleDeg, ...
  msSfMinusSsSfDeg, msMfKMinusMsSfDeg, msMfUMinusMsSfDeg, ...
  nonRefCoherenceFloorK, nonRefCoherenceFloorU);
end

function value = localCaseValue(caseTable, displayName, taskSeed, snrDb, fieldName)
%LOCALCASEVALUE Return one scalar from caseTable.

value = NaN;
if ~ismember(fieldName, caseTable.Properties.VariableNames)
  return;
end
idx = find(caseTable.displayName == string(displayName) & ...
  caseTable.taskSeed == taskSeed & caseTable.snrDb == snrDb, 1, 'first');
if ~isempty(idx)
  value = caseTable.(fieldName)(idx);
end
end
