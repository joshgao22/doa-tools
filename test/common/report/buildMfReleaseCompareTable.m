function releaseCompareTable = buildMfReleaseCompareTable(caseTable)
%BUILDMFRELEASECOMPARETABLE Compare CP-U against CP-K for the same seed and SNR.

releaseCompareTable = table();
if isempty(caseTable) || height(caseTable) == 0
  return;
end
kMaskAll = endsWith(caseTable.displayName, "CP-K") & contains(caseTable.displayName, "-MF-");
uMaskAll = endsWith(caseTable.displayName, "CP-U") & contains(caseTable.displayName, "-MF-");
keyTable = unique(caseTable(uMaskAll, {'displayName', 'satMode', 'taskSeed', 'snrDb'}), 'rows', 'stable');
if isempty(keyTable) || height(keyTable) == 0
  return;
end
displayName = keyTable.displayName;
satMode = keyTable.satMode;
seed = keyTable.taskSeed;
snrDb = keyTable.snrDb;
angleDeltaDeg = NaN(height(keyTable), 1);
fdRefAbsErrDeltaHz = NaN(height(keyTable), 1);
fdRateAbsErrHzPerSec = NaN(height(keyTable), 1);
finalObjDelta = NaN(height(keyTable), 1);
fdRateMoveAbsHzPerSec = NaN(height(keyTable), 1);
cpUCandidateCount = NaN(height(keyTable), 1);
cpUSolveVariant = strings(height(keyTable), 1);
cpUBoundaryHit = false(height(keyTable), 1);
for iRow = 1:height(keyTable)
  kName = replace(displayName(iRow), "CP-U", "CP-K");
  kIdx = find(kMaskAll & caseTable.displayName == kName & caseTable.taskSeed == seed(iRow) & caseTable.snrDb == snrDb(iRow), 1, 'first');
  uIdx = find(uMaskAll & caseTable.displayName == displayName(iRow) & caseTable.taskSeed == seed(iRow) & caseTable.snrDb == snrDb(iRow), 1, 'first');
  if isempty(kIdx) || isempty(uIdx)
    continue;
  end
  angleDeltaDeg(iRow) = caseTable.sphericalAngleErrDeg(uIdx) - caseTable.sphericalAngleErrDeg(kIdx);
  fdRefAbsErrDeltaHz(iRow) = caseTable.fdRefAbsErrHz(uIdx) - caseTable.fdRefAbsErrHz(kIdx);
  fdRateAbsErrHzPerSec(iRow) = caseTable.fdRateAbsErrHzPerSec(uIdx);
  finalObjDelta(iRow) = caseTable.finalObj(uIdx) - caseTable.finalObj(kIdx);
  fdRateMoveAbsHzPerSec(iRow) = abs(caseTable.fdRateInitMoveHzPerSec(uIdx));
  cpUCandidateCount(iRow) = caseTable.candidateCount(uIdx);
  cpUSolveVariant(iRow) = caseTable.solveVariant(uIdx);
  cpUBoundaryHit(iRow) = caseTable.fdRateBoundaryHit(uIdx);
end
releaseCompareTable = table(displayName, satMode, seed, snrDb, angleDeltaDeg, fdRefAbsErrDeltaHz, ...
  fdRateAbsErrHzPerSec, finalObjDelta, fdRateMoveAbsHzPerSec, cpUCandidateCount, ...
  cpUSolveVariant, cpUBoundaryHit);
end
