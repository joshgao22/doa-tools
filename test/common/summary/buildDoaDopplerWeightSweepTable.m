function weightTable = buildDoaDopplerWeightSweepTable(alphaList, caseList, truth)
%BUILDDOADOPPLERWEIGHTSWEEPTABLE Build one compact sat-weight sweep table.

arguments
  alphaList (:, 1) double
  caseList (1, :) struct
  truth
end

numCase = numel(caseList);
if numel(alphaList) ~= numCase
  error('buildDoaDopplerWeightSweepTable:SizeMismatch', ...
    'alphaList and caseList must have the same number of entries.');
end

truthList = repmat(localBuildBaseTruth(truth), 1, numCase);
summaryTable = buildDoaDopplerCaseSummaryTable(caseList, truthList);
weightTable = table(alphaList(:), summaryTable.displayName, summaryTable.angleErrDeg, ...
  summaryTable.fdRefErrHz, summaryTable.runTimeMs, summaryTable.isResolved, ...
  'VariableNames', {'alphaSat2', 'displayName', 'angleErrDeg', ...
  'fdRefErrHz', 'runTimeMs', 'isResolved'});
end


function caseTruth = localBuildBaseTruth(truth)
%LOCALBUILDBASETRUTH Build one reusable compact truth struct.

caseTruth = struct();
caseTruth.latlonTrueDeg = reshape(getDoaDopplerFieldOrDefault(truth, 'latlonTrueDeg', nan(2, 1)), [], 1);
caseTruth.fdRefTrueHz = getDoaDopplerFieldOrDefault(truth, 'fdRefTrueHz', NaN);
caseTruth.fdRateTrueHzPerSec = getDoaDopplerFieldOrDefault(truth, 'fdRateTrueHzPerSec', NaN);
caseTruth.refSatIdxGlobal = getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN);
caseTruth.selectedSatIdxGlobal = reshape(getDoaDopplerFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
caseTruth.fdSatTrueHz = reshape(getDoaDopplerFieldOrDefault(truth, 'fdSatTrueHz', []), [], 1);
if isfield(truth, 'fdRateSatTrueHzPerSec')
  caseTruth.fdRateSatTrueHzPerSec = reshape(truth.fdRateSatTrueHzPerSec, [], 1);
end
end
