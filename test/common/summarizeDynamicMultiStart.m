function multiStartTable = summarizeDynamicMultiStart(caseList, truth)
%SUMMARIZEDYNAMICMULTISTART Summarize CP-U multi-start release diagnostics.
% Returns one row per outer warm-start candidate so dev scripts can verify
% whether CP-U truly released from the CP-K basin or only copied it.

arguments
  caseList (1, :) struct
  truth (1, 1) struct = struct()
end

rowCell = cell(0, 1);
for iCase = 1:numel(caseList)
  caseInfo = caseList(iCase);
  estResult = localGetFieldOrDefault(caseInfo, 'estResult', struct());
  aux = localGetFieldOrDefault(estResult, 'aux', struct());
  multiStart = localGetFieldOrDefault(aux, 'multiStart', struct());
  if isempty(multiStart) || ~localGetFieldOrDefault(multiStart, 'isEnabled', false)
    continue;
  end

  tagList = reshape(string(localGetFieldOrDefault(multiStart, 'tagList', strings(0, 1))), [], 1);
  numRow = numel(tagList);
  if numRow == 0
    continue;
  end

  selectedIdx = localGetFieldOrDefault(multiStart, 'selectedIdx', NaN);
  selectedTag = string(localGetFieldOrDefault(multiStart, 'selectedTag', ""));
  fvalList = localExpandColumn(localGetFieldOrDefault(multiStart, 'fvalList', []), numRow);
  isResolvedList = logical(localExpandColumn(localGetFieldOrDefault(multiStart, 'isResolvedList', false(numRow, 1)), numRow));
  iterList = localExpandColumn(localGetFieldOrDefault(multiStart, 'iterationList', []), numRow);
  objImproveList = localExpandColumn(localGetFieldOrDefault(multiStart, 'objImproveList', []), numRow);
  initObjList = localExpandColumn(localGetFieldOrDefault(multiStart, 'initObjList', []), numRow);
  finalObjList = localExpandColumn(localGetFieldOrDefault(multiStart, 'finalObjList', []), numRow);
  moveNormList = localExpandColumn(localGetFieldOrDefault(multiStart, 'moveNormList', []), numRow);
  angleErrDegList = localExpandColumn(localGetFieldOrDefault(multiStart, 'angleErrDegList', []), numRow);
  fdRefErrHzList = localExpandColumn(localGetFieldOrDefault(multiStart, 'fdRefErrHzList', []), numRow);
  fdRateErrList = localExpandColumn(localGetFieldOrDefault(multiStart, 'fdRateErrHzPerSecList', []), numRow);
  solveVariantList = reshape(string(localGetFieldOrDefault(multiStart, 'solveVariantList', strings(numRow, 1))), [], 1);

  if numel(solveVariantList) ~= numRow
    solveVariantList = repmat("", numRow, 1);
  end

  bestFval = min(fvalList(isfinite(fvalList)));
  if isempty(bestFval)
    bestFval = NaN;
  end

  for iRow = 1:numRow
    rowInfo = struct();
    rowInfo.displayName = string(localGetFieldOrDefault(caseInfo, 'displayName', sprintf('case%d', iCase)));
    rowInfo.dynamicMode = string(localGetFieldOrDefault(caseInfo, 'dynamicMode', ""));
    rowInfo.tag = tagList(iRow);
    rowInfo.isSelected = iRow == selectedIdx || tagList(iRow) == selectedTag;
    rowInfo.isResolved = isResolvedList(iRow);
    rowInfo.iterations = iterList(iRow);
    rowInfo.initObj = initObjList(iRow);
    rowInfo.finalObj = finalObjList(iRow);
    rowInfo.objImprove = objImproveList(iRow);
    rowInfo.fvalGapToBest = fvalList(iRow) - bestFval;
    rowInfo.moveNorm = moveNormList(iRow);
    rowInfo.angleErrDeg = angleErrDegList(iRow);
    rowInfo.fdRefErrHz = fdRefErrHzList(iRow);
    rowInfo.fdRateErrHzPerSec = fdRateErrList(iRow);
    rowInfo.solveVariant = solveVariantList(iRow);
    rowCell{end + 1, 1} = rowInfo; %#ok<AGROW>
  end
end

if isempty(rowCell)
  multiStartTable = table();
  return;
end

multiStartTable = struct2table(vertcat(rowCell{:}));

if ~isempty(truth)
  multiStartTable.truthFdRefHz = repmat(localGetFieldOrDefault(truth, 'fdRefTrueHz', NaN), height(multiStartTable), 1);
  multiStartTable.truthFdRateHzPerSec = repmat(localGetFieldOrDefault(truth, 'fdRateTrueHzPerSec', NaN), height(multiStartTable), 1);
end
end

function valueCol = localExpandColumn(valueIn, numRow)
%LOCALEXPANDCOLUMN Expand one vector-like field into a column of length numRow.

if isempty(valueIn)
  valueCol = nan(numRow, 1);
  return;
end

valueCol = reshape(valueIn, [], 1);
if numel(valueCol) == numRow
  return;
end
if numel(valueCol) == 1
  valueCol = repmat(valueCol, numRow, 1);
  return;
end

valueCol = nan(numRow, 1);
end

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
