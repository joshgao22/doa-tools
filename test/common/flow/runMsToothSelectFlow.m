function flow = runMsToothSelectFlow(modeTag, periodicFixture, subsetFixtureInput, caseWide, truth, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, weightSweepAlpha, ...
  doaOnlyOpt, staticBaseOpt, dynBaseOpt, staticMsOpt, subsetDoaHalfWidthDeg, ...
  inToothFdHalfWidthHz, inToothFdRateHalfWidthHzPerSec, inToothDoaHalfWidthDeg, ...
  toothStepHz, parallelOpt)
%RUNMSTOOTHSELECTFLOW Run one reusable MS-MF tooth-selection flow.
%
% The dev flow is:
%   1) Keep one caller-provided periodic wide result as the baseline.
%   2) Run one bank of anti-periodic subset solves to pick one reliable
%      fdRef tooth / branch.
%   3) Return to the periodic full-data window and refine inside one narrow
%      in-tooth fdRef box while keeping DoA nearly frozen.
%
% modeTag:
%   "known"   -> MS-MF-CP-K on every subset and on the final in-tooth replay.
%   "unknown" -> MS-MF-CP-U on every subset and on the final in-tooth replay.

arguments
  modeTag (1, 1) string
  periodicFixture (1, 1) struct
  subsetFixtureInput
  caseWide (1, 1) struct
  truth (1, 1) struct
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical
  weightSweepAlpha (:, 1) double
  doaOnlyOpt (1, 1) struct
  staticBaseOpt (1, 1) struct
  dynBaseOpt (1, 1) struct
  staticMsOpt (1, 1) struct
  subsetDoaHalfWidthDeg (:, 1) double
  inToothFdHalfWidthHz (1, 1) double
  inToothFdRateHalfWidthHzPerSec (1, 1) double
  inToothDoaHalfWidthDeg (:, 1) double
  toothStepHz (1, 1) double
  parallelOpt (1, 1) struct
end

modeTag = lower(string(modeTag));
if modeTag ~= "known" && modeTag ~= "unknown"
  error('runMsToothSelectFlow:InvalidMode', ...
    'modeTag must be "known" or "unknown".');
end

subsetFixtureCell = localNormalizeSubsetFixtureBank(subsetFixtureInput);
numSubset = numel(subsetFixtureCell);
subsetSummaryCell = cell(numSubset, 1);
subsetCaseCell = cell(numSubset, 1);
useParfor = localShouldUseParfor(parallelOpt, numSubset, 'minSubsetBankForParfor');
if useParfor
  parfor iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    subsetCaseCell{iSubset} = localRunMsCaseFromFixture(modeTag, ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
      subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, subsetDoaHalfWidthDeg(:));
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
else
  for iSubset = 1:numSubset
    subsetFixtureUse = subsetFixtureCell{iSubset};
    subsetCaseCell{iSubset} = localRunMsCaseFromFixture(modeTag, ...
      subsetFixtureUse, pilotWave, carrierFreq, sampleRate, optVerbose, ...
      weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
      localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
      subsetFixtureUse.fdRange, subsetFixtureUse.fdRateRange, subsetDoaHalfWidthDeg(:));
    subsetSummaryCell{iSubset} = localBuildCaseSummary(subsetCaseCell{iSubset}, subsetFixtureUse.truth, toothStepHz);
    subsetSummaryCell{iSubset}.subsetLabel = string(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetLabel', "subset" + string(iSubset)));
    subsetSummaryCell{iSubset}.subsetOffsetIdx = reshape(localGetFieldOrDefault(subsetFixtureUse, ...
      'subsetOffsetIdx', []), 1, []);
  end
end

bestSubsetIdx = localSelectBestSubsetSummary(subsetSummaryCell);
bestSubsetCase = subsetCaseCell{bestSubsetIdx};
bestSubsetSummary = subsetSummaryCell{bestSubsetIdx};
wideSummary = localBuildCaseSummary(caseWide, truth, toothStepHz);

fdRangeInTooth = [bestSubsetCase.estResult.fdRefEst - inToothFdHalfWidthHz, ...
  bestSubsetCase.estResult.fdRefEst + inToothFdHalfWidthHz];
if modeTag == "unknown"
  fdRateRangeInTooth = [bestSubsetCase.estResult.fdRateEst - inToothFdRateHalfWidthHzPerSec, ...
    bestSubsetCase.estResult.fdRateEst + inToothFdRateHalfWidthHzPerSec];
else
  fdRateRangeInTooth = [];
end

refineCase = localRunMsCaseFromFixture(modeTag, ...
  periodicFixture, pilotWave, carrierFreq, sampleRate, optVerbose, ...
  weightSweepAlpha(:), doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]), ...
  fdRangeInTooth, fdRateRangeInTooth, inToothDoaHalfWidthDeg(:), ...
  bestSubsetCase.estResult.doaParamEst(:));
refineSummary = localBuildCaseSummary(refineCase, truth, toothStepHz);

toothFlowTable = localBuildToothFlowSummaryTable(wideSummary, bestSubsetSummary, refineSummary);
subsetCandidateTable = localBuildSubsetCandidateSummaryTable(subsetSummaryCell);

flow = struct();
flow.modeTag = modeTag;
flow.caseWide = caseWide;
flow.caseSubset = bestSubsetCase;
flow.caseFinal = refineCase;
flow.wideSummary = wideSummary;
flow.subsetSummary = bestSubsetSummary;
flow.refineSummary = refineSummary;
flow.toothFlowTable = toothFlowTable;
flow.subsetCandidateTable = subsetCandidateTable;
flow.bestSubsetIdx = bestSubsetIdx;
flow.selectedSubsetLabel = bestSubsetSummary.subsetLabel;
flow.selectedSubsetOffsets = bestSubsetSummary.subsetOffsetIdx;
flow.fdRangeInTooth = fdRangeInTooth;
flow.fdRateRangeInTooth = fdRateRangeInTooth;
flow.usedParfor = useParfor;
end

function caseResult = localRunMsCaseFromFixture(modeTag, fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, weightSweepAlpha, doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
  staticInitDoaHalfWidth, fdRangeUse, fdRateRangeUse, initDoaHalfWidth, initDoaParam, knownSeedCase)
%LOCALRUNMSCASEFROMFIXTURE Build one fixture-local MS-MF CP-K or CP-U case.

if nargin < 15 || isempty(initDoaHalfWidth)
  initDoaHalfWidth = [0.01; 0.01];
end
if nargin < 16
  initDoaParam = [];
end
if nargin < 17
  knownSeedCase = struct();
end

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  fixture.wavelen, pilotWave, carrierFreq, sampleRate, fdRangeUse, ...
  fixture.truth, fixture.otherSatIdxGlobal, optVerbose, doaOnlyOpt, ...
  staticBaseOpt, weightSweepAlpha, staticInitDoaHalfWidth);
bestStaticMsCase = caseBundle.bestStaticMsCase;
truthDebug = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());

if isempty(initDoaParam)
  initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
end

if modeTag == "known"
  dynMsKnownOpt = dynBaseOpt;
  dynMsKnownOpt.initDoaParam = initDoaParam(:);
  dynMsKnownOpt.initDoaHalfWidth = initDoaHalfWidth(:);
  dynMsKnownOpt.enableFdAliasUnwrap = true;
  caseResult = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
    fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
    fdRangeUse, [], optVerbose, dynMsKnownOpt, true, ...
    truthDebug, buildDynamicInitParamFromCase(bestStaticMsCase, true, fixture.truth.fdRateFit));
  return;
end

% unknown-rate: always seed from one known-rate solve built on the same fixture.
if isempty(knownSeedCase)
  knownSeedCase = localRunMsCaseFromFixture("known", fixture, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, weightSweepAlpha, doaOnlyOpt, staticBaseOpt, dynBaseOpt, ...
    staticInitDoaHalfWidth, fdRangeUse, [], initDoaHalfWidth, initDoaParam);
end

fdRateSeedKnown = localGetCaseFdRateEst(knownSeedCase, fixture.truth.fdRateFit);
initParamMsUnknownCpK = buildDynamicInitParamFromCase(knownSeedCase, false, fdRateSeedKnown);
initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = initDoaParam(:);
dynMsUnknownOpt.initDoaHalfWidth = initDoaHalfWidth(:);
dynMsUnknownOpt.enableFdAliasUnwrap = true;
msUnknownCand = buildUnknownInitCandidateSet( ...
  knownSeedCase, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  dynMsUnknownOpt.initDoaHalfWidth, caseBundle.staticMsOpt.initDoaHalfWidth);
caseResult = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, ...
  truthDebug, msUnknownCand);
end

function summary = localBuildCaseSummary(caseUse, truth, toothStepHz)
%LOCALBUILDCASESUMMARY Build one compact tooth-selection summary.

summary = struct();
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
summary.subsetLabel = string.empty(1, 0);
summary.subsetOffsetIdx = [];
end

function bestIdx = localSelectBestSubsetSummary(summaryCell)
%LOCALSELECTBESTSUBSETSUMMARY Select one subset with tooth-first ranking.

numCase = numel(summaryCell);
scoreMat = inf(numCase, 6);
for iCase = 1:numCase
  s = summaryCell{iCase};
  if ~logical(s.isResolved)
    continue;
  end
  scoreMat(iCase, :) = [0, abs(s.toothIdx), abs(s.toothResidualHz), ...
    s.angleErrDeg, abs(s.fdRateErrHzPerSec), abs(s.fdRefErrHz)];
end
[~, orderIdx] = sortrows(scoreMat, 1:size(scoreMat, 2));
bestIdx = orderIdx(1);
end

function toothTable = localBuildToothFlowSummaryTable(wideSummary, subsetSummary, refineSummary)
%LOCALBUILDTOOTHFLOWSUMMARYTABLE Build one compact anti-periodic flow table.

stageName = ["periodic-wide"; "subset-select"; "periodic-in-tooth"];
solveVariant = [wideSummary.solveVariant; subsetSummary.solveVariant; refineSummary.solveVariant];
isResolved = [wideSummary.isResolved; subsetSummary.isResolved; refineSummary.isResolved];
angleErrDeg = [wideSummary.angleErrDeg; subsetSummary.angleErrDeg; refineSummary.angleErrDeg];
fdRefErrHz = [wideSummary.fdRefErrHz; subsetSummary.fdRefErrHz; refineSummary.fdRefErrHz];
fdRateErrHzPerSec = [wideSummary.fdRateErrHzPerSec; subsetSummary.fdRateErrHzPerSec; refineSummary.fdRateErrHzPerSec];
toothIdx = [wideSummary.toothIdx; subsetSummary.toothIdx; refineSummary.toothIdx];
toothResidualHz = [wideSummary.toothResidualHz; subsetSummary.toothResidualHz; refineSummary.toothResidualHz];
runTimeMs = [wideSummary.runTimeMs; subsetSummary.runTimeMs; refineSummary.runTimeMs];

toothTable = table(stageName, solveVariant, isResolved, angleErrDeg, ...
  fdRefErrHz, fdRateErrHzPerSec, toothIdx, toothResidualHz, runTimeMs);
end

function subsetTable = localBuildSubsetCandidateSummaryTable(summaryCell)
%LOCALBUILDSUBSETCANDIDATESUMMARYTABLE Build one compact subset-bank table.

numCase = numel(summaryCell);
subsetLabel = strings(numCase, 1);
subsetOffsetStr = strings(numCase, 1);
isResolved = false(numCase, 1);
angleErrDeg = nan(numCase, 1);
fdRefErrHz = nan(numCase, 1);
fdRateErrHzPerSec = nan(numCase, 1);
toothIdx = nan(numCase, 1);
toothResidualHz = nan(numCase, 1);
for iCase = 1:numCase
  s = summaryCell{iCase};
  subsetLabel(iCase) = s.subsetLabel;
  subsetOffsetStr(iCase) = localFormatIntegerRow(s.subsetOffsetIdx);
  isResolved(iCase) = logical(s.isResolved);
  angleErrDeg(iCase) = s.angleErrDeg;
  fdRefErrHz(iCase) = s.fdRefErrHz;
  fdRateErrHzPerSec(iCase) = s.fdRateErrHzPerSec;
  toothIdx(iCase) = s.toothIdx;
  toothResidualHz(iCase) = s.toothResidualHz;
end
subsetTable = table(subsetLabel, subsetOffsetStr, isResolved, angleErrDeg, ...
  fdRefErrHz, fdRateErrHzPerSec, toothIdx, toothResidualHz);
end

function subsetFixtureCell = localNormalizeSubsetFixtureBank(subsetFixtureInput)
%LOCALNORMALIZESUBSETFIXTUREBANK Convert one struct or cell bank to one cell bank.

if iscell(subsetFixtureInput)
  subsetFixtureCell = reshape(subsetFixtureInput, [], 1);
else
  subsetFixtureCell = {subsetFixtureInput};
end
if isempty(subsetFixtureCell)
  error('runMsToothSelectFlow:EmptySubsetBank', ...
    'The subset fixture bank must contain at least one candidate.');
end
end

function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row for compact console printing.

if isempty(valueRow)
  textOut = "[]";
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '));
textOut = string(textOut);
end

function useParfor = localShouldUseParfor(parallelOpt, numCase, minFieldName)
%LOCALSHOULDUSEPARFOR Decide whether the current loop should use parfor.

useParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableSubsetSolveParfor', ...
  localGetFieldOrDefault(parallelOpt, 'enableParfor', true)));
if ~useParfor
  return;
end
minCase = localGetFieldOrDefault(parallelOpt, minFieldName, 12);
useParfor = useParfor && (numCase >= minCase) && localCanUseParfor();
end

function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether parfor can be used safely here.

tf = false;
if ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  tf = ~isempty(ver('parallel'));
catch
  tf = false;
end
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

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end

function fdRateEst = localGetCaseFdRateEst(caseUse, defaultValue)
%LOCALGETCASEFDRATEEST Safely read fdRateEst from one case struct.

if nargin < 2 || isempty(defaultValue)
  defaultValue = 0;
end
fdRateEst = defaultValue;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult)
  return;
end
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', defaultValue);
if isempty(fdRateEst) || ~isscalar(fdRateEst) || ~isfinite(fdRateEst)
  fdRateEst = defaultValue;
end
end
