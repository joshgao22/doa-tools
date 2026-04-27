function batchOut = runDynamicFixedSnrRepeatBatch(context, snrDb, numRepeat, baseSeed, ...
  flowOpt, caseName, idxSsDynUnknown, idxMsDynUnknown, representativeLabelList, runtimeOpt, optVerbose)
%RUNDYNAMICFIXEDSNRREPEATBATCH Run one fixed-SNR dynamic repeat batch.
% This helper keeps the dev entry focused on orchestration while reusing the
% shared repeat-data / bundle path for fast repeat statistics.

arguments
  context (1,1) struct
  snrDb (1,1) double
  numRepeat (1,1) double {mustBeInteger, mustBePositive}
  baseSeed (1,1) double
  flowOpt (1,1) struct
  caseName (:,1) string
  idxSsDynUnknown (1,1) double {mustBeInteger, mustBePositive}
  idxMsDynUnknown (1,1) double {mustBeInteger, mustBePositive}
  representativeLabelList (:,1) string = ["bestMsUnknownAngleGain"; "worstMsUnknownAngleGain"; "medianMsUnknown"]
  runtimeOpt (1,1) struct = struct()
  optVerbose (1,1) logical = false
end

numCase = numel(caseName);
enableWeightSweep = logical(localGetFieldOrDefault(flowOpt, 'enableWeightSweep', false));
weightSweepAlpha = reshape(localGetFieldOrDefault(flowOpt, 'weightSweepAlpha', zeros(0, 1)), [], 1);
idxWeightStart = 9;

angleErrDeg = nan(numRepeat, numCase);
fdErrHz = nan(numRepeat, numCase);
fdRateErrHzPerSec = nan(numRepeat, numCase);
isResolved = false(numRepeat, numCase);
bestWeightIdx = nan(numRepeat, 1);
bestWeightAlpha = nan(numRepeat, 1);
bestWeightAngleErrDeg = nan(numRepeat, 1);
bestWeightFdErrHz = nan(numRepeat, 1);
selectedSubsetLabel = strings(numRepeat, 1);
selectedFinalTag = strings(numRepeat, 1);
selectedFinalSolveVariant = strings(numRepeat, 1);

progressbar('reset', numRepeat);
repeatProgressQueue = parallel.pool.DataQueue;
afterEach(repeatProgressQueue, @(~) progressbar('advance'));

repeatTimer = tic;
repeatOutCell = cell(numRepeat, 1);
parfor iRepeat = 1:numRepeat
  taskSeed = baseSeed + iRepeat - 1;
  repeatOutCell{iRepeat} = localRunRepeatTask( ...
    context, snrDb, taskSeed, flowOpt, optVerbose, numCase, idxMsDynUnknown, weightSweepAlpha);
  send(repeatProgressQueue, iRepeat);
end
progressbar('end');
repeatLoopSec = toc(repeatTimer);

for iRepeat = 1:numRepeat
  repeatOut = repeatOutCell{iRepeat};
  angleErrDeg(iRepeat, :) = repeatOut.angleErrDeg;
  fdErrHz(iRepeat, :) = repeatOut.fdErrHz;
  fdRateErrHzPerSec(iRepeat, :) = repeatOut.fdRateErrHzPerSec;
  isResolved(iRepeat, :) = repeatOut.isResolved;
  bestWeightIdx(iRepeat) = repeatOut.bestWeightIdx;
  bestWeightAlpha(iRepeat) = repeatOut.bestWeightAlpha;
  bestWeightAngleErrDeg(iRepeat) = repeatOut.bestWeightAngleErrDeg;
  bestWeightFdErrHz(iRepeat) = repeatOut.bestWeightFdErrHz;
  selectedSubsetLabel(iRepeat) = repeatOut.selectedSubsetLabel;
  selectedFinalTag(iRepeat) = repeatOut.selectedFinalTag;
  selectedFinalSolveVariant(iRepeat) = repeatOut.selectedFinalSolveVariant;
end

repeatIdx = (1:numRepeat).';
taskSeed = baseSeed + repeatIdx - 1;
repeatTable = buildDynamicRepeatCompareTable(repeatIdx, taskSeed, snrDb, angleErrDeg, fdErrHz, ...
  fdRateErrHzPerSec, isResolved, selectedSubsetLabel, selectedFinalTag, ...
  bestWeightAlpha, bestWeightAngleErrDeg, bestWeightFdErrHz, idxSsDynUnknown, idxMsDynUnknown, enableWeightSweep);
repTable = selectDynamicRepresentativeRepeats(repeatTable, representativeLabelList);
repeatTimingTable = localBuildRepeatTimingSummaryTable(repeatOutCell);
numSlowTimingRepeat = localGetFieldOrDefault(runtimeOpt, 'numSlowTimingRepeat', 5);
slowRepeatTimingTable = localBuildSlowRepeatTimingTable(repeatOutCell, baseSeed, numSlowTimingRepeat);

batchOut = struct();
batchOut.angleErrDeg = angleErrDeg;
batchOut.fdErrHz = fdErrHz;
batchOut.fdRateErrHzPerSec = fdRateErrHzPerSec;
batchOut.isResolved = isResolved;
batchOut.bestWeightIdx = bestWeightIdx;
batchOut.bestWeightAlpha = bestWeightAlpha;
batchOut.bestWeightAngleErrDeg = bestWeightAngleErrDeg;
batchOut.bestWeightFdErrHz = bestWeightFdErrHz;
batchOut.selectedSubsetLabel = selectedSubsetLabel;
batchOut.selectedFinalTag = selectedFinalTag;
batchOut.selectedFinalSolveVariant = selectedFinalSolveVariant;
batchOut.repeatTable = repeatTable;
batchOut.repTable = repTable;
batchOut.repeatTimingTable = repeatTimingTable;
batchOut.slowRepeatTimingTable = slowRepeatTimingTable;
batchOut.repeatOutCell = repeatOutCell;
batchOut.repeatLoopSec = repeatLoopSec;
batchOut.enableWeightSweep = enableWeightSweep;
batchOut.weightSweepAlpha = weightSweepAlpha;
batchOut.idxWeightStart = idxWeightStart;
end

function repeatOut = localRunRepeatTask(context, snrDb, taskSeed, flowOpt, optVerbose, numCase, idxMsDynUnknown, weightSweepAlpha)
%LOCALRUNREPEATTASK Run one repeat through the shared repeat-data path.

repeatTimer = tic;
repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);

bundleTimer = tic;
bundle = buildDoaDopplerDynamicTransitionBundle( ...
  repeatData.periodicFixture, repeatData.subsetFixtureCell, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  optVerbose, flowOpt);
bundleSec = toc(bundleTimer);

caseList = bundle.caseResult;
if isfield(bundle, 'caseDynMsUnknownFinal') && numel(caseList) >= idxMsDynUnknown
  caseList(idxMsDynUnknown) = bundle.caseDynMsUnknownFinal;
end

angleVec = nan(1, numCase);
fdVec = nan(1, numCase);
fdRateVec = nan(1, numCase);
resolvedVec = false(1, numCase);
for iCase = 1:min(numel(caseList), numCase)
  [angleVec(iCase), fdVec(iCase), fdRateVec(iCase), resolvedVec(iCase)] = extractDynamicCaseMetric(caseList(iCase), repeatData.truth);
end

if numCase >= idxMsDynUnknown
  finalSummary = localGetFieldOrDefault(bundle, 'selectedFinalSummary', struct());
  angleOverride = localGetFieldOrDefault(finalSummary, 'angleErrDeg', NaN);
  fdOverride = abs(localGetFieldOrDefault(finalSummary, 'fdRefErrHz', NaN));
  fdRateOverride = abs(localGetFieldOrDefault(finalSummary, 'fdRateErrHzPerSec', NaN));
  resolvedOverride = logical(localGetFieldOrDefault(finalSummary, 'isResolved', resolvedVec(idxMsDynUnknown)));
  if isfinite(angleOverride)
    angleVec(idxMsDynUnknown) = angleOverride;
  end
  if isfinite(fdOverride)
    fdVec(idxMsDynUnknown) = fdOverride;
  end
  if isfinite(fdRateOverride)
    fdRateVec(idxMsDynUnknown) = fdRateOverride;
  end
  resolvedVec(idxMsDynUnknown) = resolvedOverride;
end

weightAngle = zeros(0, 1);
weightFd = zeros(0, 1);
weightResolved = false(0, 1);
if numCase >= 9
  weightAngle = angleVec(9:end).';
  weightFd = fdVec(9:end).';
  weightResolved = resolvedVec(9:end).' & isfinite(weightAngle);
end
bestWeightIdx = NaN;
bestWeightAlpha = NaN;
bestWeightAngleErrDeg = NaN;
bestWeightFdErrHz = NaN;
if ~isempty(weightSweepAlpha) && any(weightResolved)
  validIdx = find(weightResolved);
  [~, bestLocal] = min(weightAngle(weightResolved));
  bestWeightIdx = validIdx(bestLocal);
  if bestWeightIdx <= numel(weightSweepAlpha)
    bestWeightAlpha = weightSweepAlpha(bestWeightIdx);
  end
  bestWeightAngleErrDeg = weightAngle(bestWeightIdx);
  bestWeightFdErrHz = weightFd(bestWeightIdx);
end

repeatOut = struct();
repeatOut.angleErrDeg = angleVec;
repeatOut.fdErrHz = fdVec;
repeatOut.fdRateErrHzPerSec = fdRateVec;
repeatOut.isResolved = resolvedVec;
repeatOut.bestWeightIdx = bestWeightIdx;
repeatOut.bestWeightAlpha = bestWeightAlpha;
repeatOut.bestWeightAngleErrDeg = bestWeightAngleErrDeg;
repeatOut.bestWeightFdErrHz = bestWeightFdErrHz;
repeatOut.selectedSubsetLabel = string(localGetFieldOrDefault(bundle, 'selectedSubsetLabel', "none"));
repeatOut.selectedFinalTag = string(localGetFieldOrDefault(bundle, 'selectedFinalTag', "none"));
repeatOut.selectedFinalSolveVariant = string(localGetFieldOrDefault(localGetFieldOrDefault(bundle, 'selectedFinalSummary', struct()), 'solveVariant', ""));
repeatOut.timing = struct( ...
  'genSnapshotSec', repeatData.genSnapshotSec, ...
  'buildPeriodicFixtureSec', repeatData.buildPeriodicFixtureSec, ...
  'buildSubsetFixtureSec', repeatData.buildSubsetFixtureSec, ...
  'buildTransitionBundleSec', bundleSec, ...
  'totalRepeatSec', toc(repeatTimer));
if isfield(bundle, 'stageTiming')
  repeatOut.timing.stageTiming = bundle.stageTiming;
end
end

function repeatTimingTable = localBuildRepeatTimingSummaryTable(repeatOutCell)
%LOCALBUILDREPEATTIMINGSUMMARYTABLE Summarize repeat-level stage timings.

fieldList = { ...
  'genSnapshotSec', ...
  'buildPeriodicFixtureSec', ...
  'buildSubsetFixtureSec', ...
  'buildTransitionBundleSec', ...
  'totalRepeatSec'};
statName = ["mean"; "median"; "p95"; "max"];
valueMat = nan(numel(fieldList), numel(statName));
numRepeat = numel(repeatOutCell);
for iField = 1:numel(fieldList)
  valueVec = nan(numRepeat, 1);
  for iRepeat = 1:numRepeat
    timingUse = localGetRepeatTimingStruct(repeatOutCell{iRepeat});
    valueVec(iRepeat) = localGetTimingField(timingUse, fieldList{iField});
  end
  valueMat(iField, 1) = mean(valueVec, 'omitnan');
  valueMat(iField, 2) = median(valueVec, 'omitnan');
  valueMat(iField, 3) = localComputePercentile(valueVec, 95);
  valueMat(iField, 4) = max(valueVec, [], 'omitnan');
end
repeatTimingTable = table(string(fieldList(:)), valueMat(:, 1), valueMat(:, 2), valueMat(:, 3), valueMat(:, 4), ...
  'VariableNames', {'stageName', 'meanSec', 'medianSec', 'p95Sec', 'maxSec'});
end

function slowTable = localBuildSlowRepeatTimingTable(repeatOutCell, baseSeed, maxRow)
%LOCALBUILDSLOWREPEATTIMINGTABLE Build one slow-repeat timing table.

numRepeat = numel(repeatOutCell);
repeatIdx = (1:numRepeat).';
taskSeed = baseSeed + repeatIdx - 1;
totalRepeatSec = nan(numRepeat, 1);
genSnapshotSec = nan(numRepeat, 1);
buildPeriodicFixtureSec = nan(numRepeat, 1);
buildSubsetFixtureSec = nan(numRepeat, 1);
buildTransitionBundleSec = nan(numRepeat, 1);
for iRepeat = 1:numRepeat
  timingUse = localGetRepeatTimingStruct(repeatOutCell{iRepeat});
  totalRepeatSec(iRepeat) = localGetTimingField(timingUse, 'totalRepeatSec');
  genSnapshotSec(iRepeat) = localGetTimingField(timingUse, 'genSnapshotSec');
  buildPeriodicFixtureSec(iRepeat) = localGetTimingField(timingUse, 'buildPeriodicFixtureSec');
  buildSubsetFixtureSec(iRepeat) = localGetTimingField(timingUse, 'buildSubsetFixtureSec');
  buildTransitionBundleSec(iRepeat) = localGetTimingField(timingUse, 'buildTransitionBundleSec');
end
slowTable = table(repeatIdx, taskSeed, totalRepeatSec, genSnapshotSec, ...
  buildPeriodicFixtureSec, buildSubsetFixtureSec, buildTransitionBundleSec, ...
  'VariableNames', {'repeatIdx', 'taskSeed', 'totalRepeatSec', 'genSnapshotSec', ...
  'buildPeriodicFixtureSec', 'buildSubsetFixtureSec', 'buildTransitionBundleSec'});
slowTable = sortrows(slowTable, 'totalRepeatSec', 'descend');
if nargin >= 3 && ~isempty(maxRow)
  slowTable = slowTable(1:min(height(slowTable), maxRow), :);
end
end

function timingUse = localGetRepeatTimingStruct(repeatOut)
%LOCALGETREPEATTIMINGSTRUCT Return one repeat-level timing struct.

timingUse = struct();
if isempty(repeatOut) || ~isstruct(repeatOut)
  return;
end
if isfield(repeatOut, 'timing') && isstruct(repeatOut.timing)
  timingUse = repeatOut.timing;
end
end

function fieldValue = localGetTimingField(timingUse, fieldName)
%LOCALGETTIMINGFIELD Read one timing field with NaN fallback.

fieldValue = NaN;
if isempty(timingUse) || ~isstruct(timingUse)
  return;
end
if isfield(timingUse, fieldName)
  fieldValue = timingUse.(fieldName);
end
end

function pValue = localComputePercentile(valueVec, pct)
%LOCALCOMPUTEPERCENTILE Compute one percentile without toolbox dependence.

pValue = NaN;
valueVec = sort(valueVec(isfinite(valueVec)));
if isempty(valueVec)
  return;
end
if numel(valueVec) == 1
  pValue = valueVec;
  return;
end
rankPos = 1 + (numel(valueVec) - 1) * pct / 100;
rankLo = floor(rankPos);
rankHi = ceil(rankPos);
weightHi = rankPos - rankLo;
if rankLo == rankHi
  pValue = valueVec(rankLo);
else
  pValue = (1 - weightHi) * valueVec(rankLo) + weightHi * valueVec(rankHi);
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct field with a default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end
