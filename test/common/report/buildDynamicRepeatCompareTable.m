function repeatTable = buildDynamicRepeatCompareTable(repeatIdx, taskSeed, snrDb, angleErrDeg, fdErrHz, ...
  fdRateErrHzPerSec, isResolved, selectedSubsetLabel, selectedFinalTag, ...
  bestWeightAlpha, bestWeightAngleErrDeg, bestWeightFdErrHz, idxSsDynUnknown, idxMsDynUnknown, includeWeightSweep)
%BUILDDYNAMICREPEATCOMPARETABLE Build one compact repeat-level compare table.

arguments
  repeatIdx (:, 1) double
  taskSeed (:, 1) double
  snrDb (1, 1) double
  angleErrDeg (:, :) double
  fdErrHz (:, :) double
  fdRateErrHzPerSec (:, :) double
  isResolved (:, :) logical
  selectedSubsetLabel (:, 1) string
  selectedFinalTag (:, 1) string
  bestWeightAlpha (:, 1) double
  bestWeightAngleErrDeg (:, 1) double
  bestWeightFdErrHz (:, 1) double
  idxSsDynUnknown (1, 1) double {mustBeInteger, mustBePositive}
  idxMsDynUnknown (1, 1) double {mustBeInteger, mustBePositive}
  includeWeightSweep (1, 1) logical = false
end

repeatTable = table();
repeatTable.repeatIdx = repeatIdx;
repeatTable.taskSeed = taskSeed;
repeatTable.snrDb = repmat(snrDb, numel(repeatIdx), 1);
repeatTable.ssAngleErrDeg = angleErrDeg(:, idxSsDynUnknown);
repeatTable.msAngleErrDeg = angleErrDeg(:, idxMsDynUnknown);
repeatTable.msMinusSsAngleDeg = angleErrDeg(:, idxMsDynUnknown) - angleErrDeg(:, idxSsDynUnknown);
repeatTable.ssFdErrHz = fdErrHz(:, idxSsDynUnknown);
repeatTable.msFdErrHz = fdErrHz(:, idxMsDynUnknown);
repeatTable.msMinusSsFdHz = fdErrHz(:, idxMsDynUnknown) - fdErrHz(:, idxSsDynUnknown);
repeatTable.ssFdRateErrHzPerSec = fdRateErrHzPerSec(:, idxSsDynUnknown);
repeatTable.msFdRateErrHzPerSec = fdRateErrHzPerSec(:, idxMsDynUnknown);
repeatTable.msMinusSsFdRateHzPerSec = fdRateErrHzPerSec(:, idxMsDynUnknown) - fdRateErrHzPerSec(:, idxSsDynUnknown);
repeatTable.ssResolved = isResolved(:, idxSsDynUnknown);
repeatTable.msResolved = isResolved(:, idxMsDynUnknown);
repeatTable.selectedSubsetLabel = selectedSubsetLabel;
repeatTable.selectedFinalTag = selectedFinalTag;
if includeWeightSweep
  repeatTable.bestAlphaSat2 = bestWeightAlpha;
  repeatTable.bestWeightAngleErrDeg = bestWeightAngleErrDeg;
  repeatTable.bestWeightFdErrHz = bestWeightFdErrHz;
  repeatTable.bestWeightMinusSsAngleDeg = bestWeightAngleErrDeg - angleErrDeg(:, idxSsDynUnknown);
  repeatTable.bestWeightMinusMsAngleDeg = bestWeightAngleErrDeg - angleErrDeg(:, idxMsDynUnknown);
end
repeatTable = sortrows(repeatTable, 'msMinusSsAngleDeg', 'descend');
end
