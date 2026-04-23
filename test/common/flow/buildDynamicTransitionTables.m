function tableBundle = buildDynamicTransitionTables(stageTiming, subsetFixtureCell, subsetSummaryCell, ...
  subsetRunTimeSec, selectedSummaryCell, selectedTagList, selectedFinalTag, finalScoreMat, ...
  fdWideUnknownSummary, wideUnknownSummary, selectedSubsetSummary, anchorUnknownSummary, ...
  anchorDoaPolishSummary, refinedUnknownSummary, hasAnchorDoaPolish, didRunInToothRefine)
%BUILDDYNAMICTRANSITIONTABLES Build compact tables for the dynamic bundle.
% Keep reporting/table assembly out of the main transition orchestrator.

arguments
  stageTiming (1,1) struct
  subsetFixtureCell cell
  subsetSummaryCell cell
  subsetRunTimeSec (:,1) double
  selectedSummaryCell cell
  selectedTagList (:,1) string
  selectedFinalTag (1,1) string
  finalScoreMat double
  fdWideUnknownSummary (1,1) struct
  wideUnknownSummary (1,1) struct
  selectedSubsetSummary (1,1) struct
  anchorUnknownSummary (1,1) struct
  anchorDoaPolishSummary (1,1) struct
  refinedUnknownSummary (1,1) struct
  hasAnchorDoaPolish (1,1) logical
  didRunInToothRefine (1,1) logical
end

toothStageName = ["periodic-fd-wide"; "periodic-wide"; "subset-select"; "periodic-fd-anchor"];
toothSummaryCell = {fdWideUnknownSummary; wideUnknownSummary; selectedSubsetSummary; anchorUnknownSummary};
if hasAnchorDoaPolish
  toothStageName = [toothStageName; "periodic-fd-anchor-doa-polish"];
  toothSummaryCell = [toothSummaryCell; {anchorDoaPolishSummary}];
end
if didRunInToothRefine
  toothStageName = [toothStageName; "periodic-in-tooth"];
  toothSummaryCell = [toothSummaryCell; {refinedUnknownSummary}];
end

tableBundle = struct();
tableBundle.toothFlowTable = localBuildToothFlowSummaryTable(toothStageName, toothSummaryCell);
tableBundle.finalSelectTable = localBuildFinalSelectTable(selectedSummaryCell, selectedTagList, ...
  selectedFinalTag, finalScoreMat);
tableBundle.timingTable = localBuildDynamicTimingTable(stageTiming);
tableBundle.subsetTimingTable = localBuildSubsetTimingTable(subsetFixtureCell, subsetSummaryCell, subsetRunTimeSec);
end


function toothTable = localBuildToothFlowSummaryTable(stageName, summaryCell)
%LOCALBUILDTOOTHFLOWSUMMARYTABLE Build one compact tooth-selection summary.

numStage = numel(summaryCell);
solveVariant = strings(numStage, 1);
isResolved = false(numStage, 1);
angleErrDeg = nan(numStage, 1);
fdRefErrHz = nan(numStage, 1);
fdRateErrHzPerSec = nan(numStage, 1);
toothIdx = nan(numStage, 1);
toothResidualHz = nan(numStage, 1);
runTimeMs = nan(numStage, 1);
for iStage = 1:numStage
  s = summaryCell{iStage};
  solveVariant(iStage) = string(localGetFieldOrDefault(s, 'solveVariant', ""));
  isResolved(iStage) = localGetFieldOrDefault(s, 'isResolved', false);
  angleErrDeg(iStage) = localGetFieldOrDefault(s, 'angleErrDeg', NaN);
  fdRefErrHz(iStage) = localGetFieldOrDefault(s, 'fdRefErrHz', NaN);
  fdRateErrHzPerSec(iStage) = localGetFieldOrDefault(s, 'fdRateErrHzPerSec', NaN);
  toothIdx(iStage) = localGetFieldOrDefault(s, 'toothIdx', NaN);
  toothResidualHz(iStage) = localGetFieldOrDefault(s, 'toothResidualHz', NaN);
  runTimeMs(iStage) = localGetFieldOrDefault(s, 'runTimeMs', NaN);
end
toothTable = table(stageName, solveVariant, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs, ...
  'VariableNames', {'stageName', 'solveVariant', 'isResolved', 'angleErrDeg', ...
  'fdRefErrHz', 'fdRateErrHzPerSec', 'toothIdx', 'toothResidualHz', 'runTimeMs'});
end


function timingTable = localBuildDynamicTimingTable(stageTiming)
%LOCALBUILDDYNAMICTIMINGTABLE Build one compact stage-timing table.

stageName = ["staticTransitionSec"; "refKnownSec"; "refUnknownSec"; "msKnownSec"; ...
  "msUnknownFdWideSec"; "msUnknownWideSec"; "subsetSelectTotalSec"; ...
  "subsetRankingSec"; "subsetAnchorSec"; "anchorDoaPolishSec"; ...
  "inToothRefineSec"; "finalSelectSec"; "summaryBuildSec"; "totalBundleSec"];
runTimeSec = nan(numel(stageName), 1);
for iStage = 1:numel(stageName)
  runTimeSec(iStage) = localGetFieldOrDefault(stageTiming, char(stageName(iStage)), NaN);
end
timingTable = table(stageName, runTimeSec, 'VariableNames', {'stageName', 'runTimeSec'});
end


function subsetTimingTable = localBuildSubsetTimingTable(subsetFixtureCell, summaryCell, subsetRunTimeSec)
%LOCALBUILDSUBSETTIMINGTABLE Build one compact subset timing table.

numCase = numel(subsetFixtureCell);
if numCase == 0
  subsetTimingTable = table();
  return;
end

labelOut = strings(numCase, 1);
runTimeSec = nan(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
for iCase = 1:numCase
  labelOut(iCase) = string(localGetFieldOrDefault(subsetFixtureCell{iCase}, 'subsetLabel', sprintf('cand%d', iCase)));
  runTimeSec(iCase) = subsetRunTimeSec(iCase);
  summaryUse = summaryCell{iCase};
  angleErrDeg(iCase) = localGetFieldOrDefault(summaryUse, 'angleErrDeg', NaN);
  fdRefErrHz(iCase) = localGetFieldOrDefault(summaryUse, 'fdRefErrHz', NaN);
  fdRateErrHzPerSec(iCase) = localGetFieldOrDefault(summaryUse, 'fdRateErrHzPerSec', NaN);
end
subsetTimingTable = table(labelOut, runTimeSec, angleErrDeg, fdRefErrHz, fdRateErrHzPerSec, ...
  'VariableNames', {'label', 'runTimeSec', 'angleErrDeg', 'fdRefErrHz', 'fdRateErrHzPerSec'});
end


function finalSelectTable = localBuildFinalSelectTable(summaryCell, candidateTag, selectedTag, scoreMat)
%LOCALBUILDFINALSELECTTABLE Build one compact final-select table.

numCase = numel(summaryCell);
isSelected = (candidateTag == string(selectedTag));
isResolved = false(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
finalObj = nan(numCase, 1);
finalResidualNorm = nan(numCase, 1);
for iCase = 1:numCase
  s = summaryCell{iCase};
  isResolved(iCase) = s.isResolved;
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
  finalObj(iCase) = localGetFieldOrDefault(s, 'finalObj', NaN);
  finalResidualNorm(iCase) = localGetFieldOrDefault(s, 'finalResidualNorm', NaN);
end

baseVarName = {'candidateTag', 'isSelected', 'isResolved', 'angleErrDeg', 'fdRefErrHz', ...
  'fdRateErrHzPerSec', 'toothIdx', 'toothResidualHz', 'finalObj', 'finalResidualNorm'};
baseValue = {candidateTag, isSelected, isResolved, angleErrDeg, fdRefErrHz, ...
  fdRateErrHzPerSec, toothIdx, toothResidualHz, finalObj, finalResidualNorm};
scoreCell = cell(1, size(scoreMat, 2));
scoreVarName = localBuildScoreVarNames(size(scoreMat, 2));
for iCol = 1:size(scoreMat, 2)
  scoreCell{iCol} = scoreMat(:, iCol);
end
finalSelectTable = table(baseValue{:}, scoreCell{:}, 'VariableNames', [baseVarName, scoreVarName]);
end


function scoreVarName = localBuildScoreVarNames(numScore)
%LOCALBUILDSCOREVARNAMES Build score-column names for the final-select table.

scoreVarName = cellstr("score" + string(1:numScore));
templateName = {'scoreResolvedPenalty', 'scoreNonRefSupportPenalty', 'scoreNonRefFitPenalty', ...
  'scoreNonRefConsistencyPenalty', 'scoreNonRefCoherencePenalty', ...
  'scoreNonRefRmsPhaseResidRad', 'scoreNonRefMaxAbsPhaseResidRad', ...
  'scoreNonRefNegativePenalty', 'scoreFinalObj', 'scoreFinalResidualNorm', 'scoreRunTimeMs'};
numTemplate = min(numScore, numel(templateName));
scoreVarName(1:numTemplate) = templateName(1:numTemplate);
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or a default value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
