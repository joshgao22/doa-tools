function tailTable = buildMfTailTable(caseTable, opt)
%BUILDMFTAILTABLE Keep rejected or trimmed MF diagnostic rows.
% This table is report-only and is shared by replay/scan summaries.

arguments
  caseTable
  opt.SortColumns = {'snrDb', 'displayName', 'taskSeed'}
end

if isempty(caseTable) || height(caseTable) == 0
  tailTable = table();
  return;
end

tailMask = caseTable.isDynamicMethod & (~caseTable.healthResolved | ~caseTable.trimKeep);
columnList = {'displayName', 'snrDb', 'taskSeed', 'mfInitMode', 'sphericalAngleErrDeg', ...
  'angleErrOverSphericalCrb', 'fdRefAbsErrHz', 'fdRefErrOverCrb', ...
  'fdRateAbsErrHzPerSec', 'iterations', 'objectiveImprove', 'fdRefBoundaryHit', ...
  'fdRateBoundaryHit', 'healthResolved', 'trimKeep', 'rejectReason'};
columnList = columnList(ismember(columnList, caseTable.Properties.VariableNames));
tailTable = caseTable(tailMask, columnList);
if height(tailTable) > 0
  sortCols = opt.SortColumns(ismember(opt.SortColumns, tailTable.Properties.VariableNames));
  if ~isempty(sortCols)
    tailTable = sortrows(tailTable, sortCols);
  end
end
end
