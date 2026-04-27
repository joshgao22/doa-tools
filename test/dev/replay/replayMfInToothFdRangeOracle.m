% replayMfInToothFdRangeOracle
% Purpose: compare single-sat / multi-sat and single-frame / multi-frame
% DoA-Doppler baselines under a truth-centered half-tooth Doppler oracle.
% This replay tests whether the MS-MF dynamic estimator has a healthy
% in-tooth upper bound after wrong-tooth selection is factored out.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; saveSnapshot=true saves only
% replayData via saveExpSnapshot.

clear; close all; clc;

%% Replay configuration
replayName = "replayMfInToothFdRangeOracle";

snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % First seed of the small replay batch.
numRepeat = 50;                     % Number of consecutive seeds to replay.
saveSnapshot = true;                % true saves the lightweight replayData only.
optVerbose = false;                 % true enables compact estimator / branch trace.
oracleFdHalfToothFraction = 0.49;   % Half-width fraction of the 1/T_f tooth step.
oracleFdRateHalfWidthHzPerSec = 1000; % Local truth-centered fdRate half-width for oracle-rate variants.
staticLocalDoaHalfWidthDeg = [0.002; 0.002]; % Local DoA box around the static MS seed.
staticWideDoaHalfWidthDeg = [0.010; 0.010];  % Wider same-tooth DoA box for stress comparison.
truthLocalDoaHalfWidthDeg = [0.002; 0.002];  % Local DoA box around the truth DoA seed.
histogramBinCount = 12;             % Bin count for compact distribution plots.

seedList = baseSeed + (0:(numRepeat - 1));
seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

%% Build context and flow options
config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.numRepeat = numRepeat;
config.seedList = seedList;
config.saveSnapshot = saveSnapshot;
config.optVerbose = optVerbose;
config.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
config.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
config.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg(:);
config.staticWideDoaHalfWidthDeg = staticWideDoaHalfWidthDeg(:);
config.truthLocalDoaHalfWidthDeg = truthLocalDoaHalfWidthDeg(:);
config.histogramBinCount = histogramBinCount;

fprintf('Running %s ...\n', char(replayName));
fprintf('  repeats                         : %d\n', config.numRepeat);
fprintf('  snr (dB)                        : %.2f\n', config.snrDb);
fprintf('  base seed                       : %d\n', config.baseSeed);
fprintf('  repeat mode                     : %s\n', 'parfor-auto');
fprintf('  save snapshot                   : %d\n', config.saveSnapshot);
fprintf('  fd oracle half-tooth fraction   : %.3f\n', config.oracleFdHalfToothFraction);
fprintf('  fdRate oracle half-width        : %.2f Hz/s\n', config.oracleFdRateHalfWidthHzPerSec);

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
contextOpt = struct( ...
  'baseSeed', config.baseSeed, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', parallelOpt);
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false));
methodList = localBuildMethodList(config);

fprintf('[%s] Build shared dynamic context.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
context = buildDynamicDualSatEciContext(contextOpt);
context = localDisableSubsetBankForOracle(context);
fprintf('[%s] Shared dynamic context built.\n', char(datetime('now', 'Format', 'HH:mm:ss')));

%% Run replay batch
numTask = config.numRepeat;
repeatCell = cell(numTask, 1);
useParfor = localCanUseParfor() && numTask > 1;
tracker = localCreateProgressTracker(sprintf('%s repeat batch (%s)', char(replayName), localModeText(useParfor)), numTask, useParfor);
try
  if useParfor
    progressQueue = tracker.queue;
    parfor iRepeat = 1:numTask
      repeatCell{iRepeat} = localRunOneRepeat(iRepeat, config.seedList(iRepeat), context, flowOpt, methodList, config);
      if ~isempty(progressQueue)
        send(progressQueue, iRepeat);
      end
    end
  else
    for iRepeat = 1:numTask
      repeatCell{iRepeat} = localRunOneRepeat(iRepeat, config.seedList(iRepeat), context, flowOpt, methodList, config);
      localAdvanceProgressTracker(tracker);
    end
  end
  localCloseProgressTracker(tracker);
catch ME
  localCloseProgressTracker(tracker);
  rethrow(ME);
end

methodTable = localBuildMethodTable(repeatCell);
aggregateTable = localBuildAggregateTable(methodTable);
pairCompareTable = localBuildPairCompareTable(methodTable);
paperCompareTable = localBuildPaperCompareTable(methodTable);
singleMultiCompareTable = localBuildSingleMultiCompareTable(methodTable);
tailCaseTable = localBuildTailCaseTable(singleMultiCompareTable);
rangeTable = localBuildRangeTable(repeatCell);
timingTable = localBuildTimingTable(repeatCell);
timingAggregateTable = localBuildTimingAggregateTable(timingTable);
representative = localSelectRepresentative(methodTable, repeatCell);
plotData = localBuildPlotData(methodTable, singleMultiCompareTable, config.histogramBinCount);

%% Data storage
replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, methodList);
replayData.methodTable = methodTable;
replayData.aggregateTable = aggregateTable;
replayData.pairCompareTable = pairCompareTable;
replayData.paperCompareTable = paperCompareTable;
replayData.singleMultiCompareTable = singleMultiCompareTable;
replayData.tailCaseTable = tailCaseTable;
replayData.rangeTable = rangeTable;
replayData.timingTable = timingTable;
replayData.timingAggregateTable = timingAggregateTable;
replayData.representative = representative;
replayData.plotData = plotData;
replayData.methodList = localStripMethodList(methodList);

if config.saveSnapshot
  saveOpt = struct('includeVars', {{'replayData'}}, ...
    'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
  replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
else
  replayData.snapshotFile = "";
end

%% Summary output and plotting
if ~exist('replayData', 'var') || ~isstruct(replayData)
  error('Replay data is missing. Run the replay batch sections or load a snapshot containing replayData.');
end
methodTable = replayData.methodTable;
aggregateTable = replayData.aggregateTable;
pairCompareTable = replayData.pairCompareTable;
if isfield(replayData, 'paperCompareTable')
  paperCompareTable = replayData.paperCompareTable;
else
  paperCompareTable = localBuildPaperCompareTable(methodTable);
end
if isfield(replayData, 'singleMultiCompareTable')
  singleMultiCompareTable = replayData.singleMultiCompareTable;
else
  singleMultiCompareTable = localBuildSingleMultiCompareTable(methodTable);
end
if isfield(replayData, 'tailCaseTable')
  tailCaseTable = replayData.tailCaseTable;
else
  tailCaseTable = localBuildTailCaseTable(singleMultiCompareTable);
end
rangeTable = replayData.rangeTable;
if isfield(replayData, 'timingTable')
  timingTable = replayData.timingTable;
else
  timingTable = table();
end
if isfield(replayData, 'timingAggregateTable')
  timingAggregateTable = replayData.timingAggregateTable;
else
  timingAggregateTable = localBuildTimingAggregateTable(timingTable);
end
representative = replayData.representative;

fprintf('Running replayMfInToothFdRangeOracle ...\n');
fprintf('  repeats                         : %d\n', numel(unique(methodTable.taskSeed)));
fprintf('  methods                         : %d\n', numel(unique(methodTable.methodLabel)));
fprintf('  snr (dB)                        : %.2f\n', replayData.config.snrDb);
fprintf('  base seed                       : %d\n', replayData.config.baseSeed);
fprintf('\n========== Oracle method aggregate ==========\n');
disp(aggregateTable);
fprintf('\n========== MS-MF wide-vs-oracle pair compare ==========\n');
disp(pairCompareTable);
fprintf('\n========== Paper-claim upper-bound compare ==========\n');
disp(paperCompareTable);
if ~isempty(timingAggregateTable)
  fprintf('\n========== Runtime timing summary ==========\n');
  disp(timingAggregateTable);
end
fprintf('\n========== Oracle range summary ==========\n');
localDispTablePreview(rangeTable, 4);
fprintf('\n========== Single-vs-multi in-tooth tail summary ==========\n');
fprintf('  per-seed rows stored            : %d\n', height(singleMultiCompareTable));
fprintf('  MS-MF better angle rate         : %.3f\n', mean(double(singleMultiCompareTable.multiBetterAngle), 'omitnan'));
fprintf('  median angle gain SS-MF - MS-MF : %.6g deg\n', median(singleMultiCompareTable.angleGainDeg, 'omitnan'));
fprintf('\n========== Worst single-vs-multi tail cases ==========\n');
disp(tailCaseTable);
fprintf('  representative seed             : %d\n', representative.taskSeed);
fprintf('  representative method           : %s\n', representative.methodLabel);
fprintf('  representative angle err        : %.6g deg\n', representative.angleErrDeg);
fprintf('  representative fdRef err        : %.6g Hz\n', representative.fdRefErrHz);
localPlotReplay(replayData.plotData);

%% Local helpers

function context = localDisableSubsetBankForOracle(context)
%LOCALDISABLESUBSETBANKFORORACLE Skip curated subset fixtures in this oracle replay.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function methodList = localBuildMethodList(config)
%LOCALBUILDMETHODLIST Define the replay-only paper-claim method bank.

methodList = repmat(localMakeMethod("", "", "", "", "", false, false, ...
  "", "", "", [], false, "", ""), 0, 1);
methodList(end + 1, 1) = localMakeMethod("ss-mf-cp-u-wide", "baseline", "single", ...
  "multi", "dynamic-cp", false, false, "ref", "static", "staticZeroRate", ...
  config.staticWideDoaHalfWidthDeg, false, "default", "default");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-wide", "baseline", "multi", ...
  "multi", "dynamic-cp", false, false, "ms", "static", "staticZeroRate", ...
  config.staticWideDoaHalfWidthDeg, false, "default", "default");
methodList(end + 1, 1) = localMakeMethod("ss-mf-cp-u-in-tooth", "oracle", "single", ...
  "multi", "dynamic-cp", true, false, "ref", "static", "truthFreq", ...
  config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-in-tooth", "oracle", "multi", ...
  "multi", "dynamic-cp", true, false, "ms", "static", "truthFreq", ...
  config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-k-in-tooth", "oracle", "multi", ...
  "multi", "dynamic-cp", true, true, "ms", "static", "truthFreq", ...
  config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "known");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-truth-doa-oracle", "oracle-upper", "multi", ...
  "multi", "dynamic-cp", true, false, "ms", "truth", "truthFreq", ...
  config.truthLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
end

function method = localMakeMethod(label, family, satMode, frameMode, modelClass, isOracle, isKnownRate, viewMode, doaSeedMode, initMode, doaHalfWidthDeg, freezeDoa, fdRangeMode, fdRateRangeMode)
%LOCALMAKEMETHOD Construct one method descriptor with a stable field set.

method = struct();
method.label = string(label);
method.family = string(family);
method.satMode = string(satMode);
method.frameMode = string(frameMode);
method.modelClass = string(modelClass);
method.isOracle = logical(isOracle);
method.isKnownRate = logical(isKnownRate);
method.viewMode = string(viewMode);
method.doaSeedMode = string(doaSeedMode);
method.initMode = string(initMode);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
method.freezeDoa = logical(freezeDoa);
method.fdRangeMode = string(fdRangeMode);
method.fdRateRangeMode = string(fdRateRangeMode);
end

function repeatOut = localRunOneRepeat(iRepeat, taskSeed, context, flowOpt, methodList, config)
%LOCALRUNONEREPEAT Build one noisy repeat and run all oracle comparison methods.

lastwarn('', '');
tRepeat = tic;
tBuildData = tic;
repeatData = buildDynamicRepeatData(context, config.snrDb, taskSeed);
buildDataMs = 1000 * toc(tBuildData);
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('replayMfInToothFdRangeOracle:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for the in-tooth oracle replay.');
end
fdHalfWidthHz = config.oracleFdHalfToothFraction * toothStepHz;
fdRangeOracle = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
fdRateRangeOracle = truthFdRateHzPerSec + config.oracleFdRateHalfWidthHzPerSec * [-1, 1];

tStaticBundle = tic;
staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fixture.fdRange, truth, context.otherSatIdxGlobal, config.optVerbose, ...
  flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);
staticBundleMs = 1000 * toc(tStaticBundle);

rowCell = cell(numel(methodList) + 2, 1);
rowCell{1} = localBuildSummaryRow(iRepeat, taskSeed, "ss-sf-static", "baseline", ...
  "single", "single", "static", false, false, "static", "static", ...
  "default", "none", NaN, false, staticBundle.caseStaticRefOnly, truth, toothStepHz, NaN);
rowCell{2} = localBuildSummaryRow(iRepeat, taskSeed, "ms-sf-static", "baseline", ...
  "multi", "single", "static", false, false, "static", "static", ...
  "default", "none", NaN, false, staticBundle.caseStaticMs, truth, toothStepHz, NaN);

dynamicMethodWallTimeMs = nan(numel(methodList), 1);
for iMethod = 1:numel(methodList)
  method = methodList(iMethod);
  tMethod = tic;
  caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, ...
    context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, config.optVerbose, ...
    flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz);
  dynamicMethodWallTimeMs(iMethod) = 1000 * toc(tMethod);
  rowCell{iMethod + 2} = localBuildSummaryRow(iRepeat, taskSeed, method.label, method.family, ...
    method.satMode, method.frameMode, method.modelClass, method.isOracle, method.isKnownRate, ...
    method.doaSeedMode, method.initMode, method.fdRangeMode, method.fdRateRangeMode, ...
    localMaxAbs(method.doaHalfWidthDeg), method.freezeDoa, caseUse, truth, toothStepHz, ...
    dynamicMethodWallTimeMs(iMethod));
end
dynamicMethodTotalMs = sum(dynamicMethodWallTimeMs, 'omitnan');
dynamicMethodMaxMs = max(dynamicMethodWallTimeMs, [], 'omitnan');
repeatTotalMs = 1000 * toc(tRepeat);

[warningMessage, warningId] = lastwarn;
repeatOut = struct();
repeatOut.iRepeat = iRepeat;
repeatOut.taskSeed = taskSeed;
repeatOut.snrDb = config.snrDb;
repeatOut.toothStepHz = toothStepHz;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.fdRangeOracle = fdRangeOracle;
repeatOut.fdRateRangeOracle = fdRateRangeOracle;
repeatOut.methodSummary = struct2table([rowCell{:}].');
repeatOut.timing = struct( ...
  'buildDataMs', buildDataMs, ...
  'staticBundleMs', staticBundleMs, ...
  'dynamicMethodTotalMs', dynamicMethodTotalMs, ...
  'dynamicMethodMaxMs', dynamicMethodMaxMs, ...
  'repeatTotalMs', repeatTotalMs);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz)
%LOCALRUNDYNAMICMETHOD Run one dynamic method descriptor on the selected fixture view.

[viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle);
fdRangeUse = fixture.fdRange;
fdRateRangeUse = fixture.fdRateRange;
if method.fdRangeMode == "truthHalfTooth"
  fdRangeUse = fdRangeOracle;
end
if method.fdRateRangeMode == "truthLocal"
  fdRateRangeUse = fdRateRangeOracle;
end

truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
seedDoaParam = localResolveDoaSeed(method, staticCaseUse, truth);

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = method.freezeDoa;
if method.isOracle
  dynOpt.disableUnknownDoaReleaseFloor = true;
  dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
end

initOverride = [];
switch method.initMode
  case "auto"
    initOverride = [];
  case "staticZeroRate"
    initParam = buildDynamicInitParamFromCase(staticCaseUse, method.isKnownRate, 0);
    initOverride = struct('startTag', method.label, 'initParam', initParam, ...
      'initDoaParam', seedDoaParam(:), 'initDoaHalfWidth', method.doaHalfWidthDeg(:), ...
      'freezeDoa', method.freezeDoa);
  case "truthFreq"
    initParam = localBuildTruthFreqInitParam(seedDoaParam, truthFdRefHz, truthFdRateHzPerSec, method.isKnownRate);
    initOverride = struct('startTag', method.label, 'initParam', initParam, ...
      'initDoaParam', seedDoaParam(:), 'initDoaHalfWidth', method.doaHalfWidthDeg(:), ...
      'freezeDoa', method.freezeDoa);
  otherwise
    error('replayMfInToothFdRangeOracle:UnknownInitMode', ...
      'Unknown method initMode "%s".', method.initMode);
end

displayPrefix = "SS";
if method.satMode == "multi"
  displayPrefix = "MS";
end
displayName = displayPrefix + "-MF-" + method.label;
caseUse = runDynamicDoaDopplerCase(displayName, method.satMode, viewUse, truth, ...
  pilotWave, carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, optVerbose, ...
  dynOpt, method.isKnownRate, debugTruthUse, initOverride);
caseUse.oracleMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, 'fdRateRangeUse', fdRateRangeUse);
end

function [viewUse, debugTruthUse, staticCaseUse] = localResolveMethodView(method, fixture, staticBundle)
%LOCALRESOLVEMETHODVIEW Resolve view, debug truth, and static seed for one method.

switch method.viewMode
  case "ref"
    viewUse = fixture.viewRefOnly;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());
    staticCaseUse = staticBundle.caseStaticRefOnly;
  case "ms"
    viewUse = fixture.viewMs;
    debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());
    staticCaseUse = staticBundle.caseStaticMs;
  otherwise
    error('replayMfInToothFdRangeOracle:UnknownViewMode', ...
      'Unknown method viewMode "%s".', method.viewMode);
end
end

function seedDoaParam = localResolveDoaSeed(method, staticCase, truth)
%LOCALRESOLVEDOASEED Choose the DoA seed for one replay-only method.

switch method.doaSeedMode
  case "static"
    seedDoaParam = staticCase.estResult.doaParamEst(:);
  case "truth"
    seedDoaParam = reshape(truth.latlonTrueDeg, [], 1);
  otherwise
    error('replayMfInToothFdRangeOracle:UnknownDoaSeedMode', ...
      'Unknown DoA seed mode "%s".', method.doaSeedMode);
end
end

function initParam = localBuildTruthFreqInitParam(doaParam, fdRefHz, fdRateHzPerSec, isKnownRate)
%LOCALBUILDTRUTHFREQINITPARAM Build an oracle frequency initializer.

if isKnownRate
  initParam = [doaParam(:); fdRefHz];
else
  initParam = [doaParam(:); fdRefHz; fdRateHzPerSec];
end
end

function row = localBuildSummaryRow(iRepeat, taskSeed, methodLabel, methodFamily, satMode, frameMode, modelClass, isOracle, isKnownRate, doaSeedMode, initMode, fdRangeMode, fdRateRangeMode, doaHalfWidthMaxDeg, freezeDoa, caseUse, truth, toothStepHz, wallTimeMs)
%LOCALBUILDSUMMARYROW Convert one method result into a stable table row.

info = summarizeDoaDopplerCase(caseUse, truth);
dynSummary = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truth);
row = struct();
row.iRepeat = iRepeat;
row.taskSeed = taskSeed;
row.methodLabel = string(methodLabel);
row.methodFamily = string(methodFamily);
row.satMode = string(satMode);
row.frameMode = string(frameMode);
row.modelClass = string(modelClass);
row.rateMode = localInferRateMode(modelClass, isKnownRate);
row.oracleLevel = localInferOracleLevel(isOracle, fdRangeMode, fdRateRangeMode, doaSeedMode);
row.isOracle = logical(isOracle);
row.isKnownRate = logical(isKnownRate);
row.doaSeedMode = string(doaSeedMode);
row.initMode = string(initMode);
row.fdRangeMode = string(fdRangeMode);
row.fdRateRangeMode = string(fdRateRangeMode);
row.doaHalfWidthMaxDeg = doaHalfWidthMaxDeg;
row.freezeDoa = logical(freezeDoa);
row.solveVariant = string(localGetFieldOrDefault(dynSummary, 'solveVariant', ""));
row.isResolved = logical(localGetFieldOrDefault(info, 'isResolved', false));
row.angleErrDeg = localGetFieldOrDefault(info, 'angleErrDeg', NaN);
row.fdRefErrHz = localGetFieldOrDefault(info, 'fdRefErrHz', NaN);
row.fdRateErrHzPerSec = localGetFieldOrDefault(info, 'fdRateErrHzPerSec', NaN);
row.fdLineRmseHz = localGetFieldOrDefault(info, 'fdLineRmseHz', NaN);
row.fdSatRmseHz = localGetFieldOrDefault(info, 'fdSatRmseHz', NaN);
row.deltaFdRmseHz = localGetFieldOrDefault(info, 'deltaFdRmseHz', NaN);
row.toothIdx = localGetFieldOrDefault(dynSummary, 'toothIdx', NaN);
row.toothResidualHz = localGetFieldOrDefault(dynSummary, 'toothResidualHz', NaN);
row.finalObj = localGetFieldOrDefault(dynSummary, 'finalObj', NaN);
row.finalResidualNorm = localGetFieldOrDefault(dynSummary, 'finalResidualNorm', NaN);
row.nonRefCoherenceFloor = localGetFieldOrDefault(dynSummary, 'nonRefCoherenceFloor', NaN);
row.nonRefRmsPhaseResidRad = localGetFieldOrDefault(dynSummary, 'nonRefRmsPhaseResidRad', NaN);
row.runTimeMs = localGetFieldOrDefault(info, 'runTimeMs', NaN);
row.wallTimeMs = wallTimeMs;
end

function rateMode = localInferRateMode(modelClass, isKnownRate)
%LOCALINFERRATEMODE Return a compact rate-mode label for summary tables.

if string(modelClass) == "static"
  rateMode = "static";
elseif isKnownRate
  rateMode = "known-rate";
else
  rateMode = "unknown-rate";
end
end

function oracleLevel = localInferOracleLevel(isOracle, fdRangeMode, fdRateRangeMode, doaSeedMode)
%LOCALINFERORACLELEVEL Return the oracle scope used by one replay row.

if ~isOracle
  oracleLevel = "none";
  return;
end
labelPart = strings(0, 1);
if string(fdRangeMode) == "truthHalfTooth"
  labelPart(end + 1, 1) = "fd";
end
if string(fdRateRangeMode) == "truthLocal"
  labelPart(end + 1, 1) = "fdRate";
elseif string(fdRateRangeMode) == "known"
  labelPart(end + 1, 1) = "knownRate";
end
if string(doaSeedMode) == "truth"
  labelPart(end + 1, 1) = "truthDoA";
end
if isempty(labelPart)
  oracleLevel = "oracle";
else
  oracleLevel = strjoin(labelPart, "+");
end
end

function methodTable = localBuildMethodTable(repeatCell)
%LOCALBUILDMETHODTABLE Stack per-repeat method summaries into one table.

summaryCell = cell(numel(repeatCell), 1);
for iRepeat = 1:numel(repeatCell)
  summaryCell{iRepeat} = repeatCell{iRepeat}.methodSummary;
end
methodTable = vertcat(summaryCell{:});
end

function aggregateTable = localBuildAggregateTable(methodTable)
%LOCALBUILDAGGREGATETABLE Build compact method-level distribution metrics.

methodList = unique(methodTable.methodLabel, 'stable');
numMethod = numel(methodList);
methodLabel = strings(numMethod, 1);
methodFamily = strings(numMethod, 1);
satMode = strings(numMethod, 1);
frameMode = strings(numMethod, 1);
modelClass = strings(numMethod, 1);
rateMode = strings(numMethod, 1);
oracleLevel = strings(numMethod, 1);
numRepeat = zeros(numMethod, 1);
resolvedRate = nan(numMethod, 1);
toothHitRate = nan(numMethod, 1);
angleRmseDeg = nan(numMethod, 1);
angleMedianDeg = nan(numMethod, 1);
angleP95Deg = nan(numMethod, 1);
fdRefAbsMedianHz = nan(numMethod, 1);
fdRefAbsP95Hz = nan(numMethod, 1);
fdRateAbsMedianHzPerSec = nan(numMethod, 1);
fdRateAbsP95HzPerSec = nan(numMethod, 1);
fdLineRmseMedianHz = nan(numMethod, 1);
nonRefCoherenceFloorMedian = nan(numMethod, 1);
finalObjMedian = nan(numMethod, 1);
runTimeMedianMs = nan(numMethod, 1);
wallTimeMedianMs = nan(numMethod, 1);

for iMethod = 1:numMethod
  methodLabel(iMethod) = methodList(iMethod);
  mask = methodTable.methodLabel == methodList(iMethod);
  idxFirst = find(mask, 1, 'first');
  methodFamily(iMethod) = methodTable.methodFamily(idxFirst);
  satMode(iMethod) = methodTable.satMode(idxFirst);
  frameMode(iMethod) = methodTable.frameMode(idxFirst);
  modelClass(iMethod) = methodTable.modelClass(idxFirst);
  rateMode(iMethod) = methodTable.rateMode(idxFirst);
  oracleLevel(iMethod) = methodTable.oracleLevel(idxFirst);
  numRepeat(iMethod) = nnz(mask);
  resolvedRate(iMethod) = mean(double(methodTable.isResolved(mask)), 'omitnan');
  toothHitRate(iMethod) = mean(double(methodTable.toothIdx(mask) == 0), 'omitnan');
  angleVal = methodTable.angleErrDeg(mask);
  fdAbs = abs(methodTable.fdRefErrHz(mask));
  fdRateAbs = abs(methodTable.fdRateErrHzPerSec(mask));
  fdLineVal = methodTable.fdLineRmseHz(mask);
  nonRefCohVal = methodTable.nonRefCoherenceFloor(mask);
  finalObjVal = methodTable.finalObj(mask);
  runTimeVal = methodTable.runTimeMs(mask);
  if ismember('wallTimeMs', methodTable.Properties.VariableNames)
    wallTimeVal = methodTable.wallTimeMs(mask);
  else
    wallTimeVal = NaN;
  end
  angleRmseDeg(iMethod) = sqrt(mean(angleVal .^ 2, 'omitnan'));
  angleMedianDeg(iMethod) = median(angleVal, 'omitnan');
  angleP95Deg(iMethod) = localPercentile(angleVal, 95);
  fdRefAbsMedianHz(iMethod) = median(fdAbs, 'omitnan');
  fdRefAbsP95Hz(iMethod) = localPercentile(fdAbs, 95);
  fdRateAbsMedianHzPerSec(iMethod) = median(fdRateAbs, 'omitnan');
  fdRateAbsP95HzPerSec(iMethod) = localPercentile(fdRateAbs, 95);
  fdLineRmseMedianHz(iMethod) = median(fdLineVal, 'omitnan');
  nonRefCoherenceFloorMedian(iMethod) = median(nonRefCohVal, 'omitnan');
  finalObjMedian(iMethod) = median(finalObjVal, 'omitnan');
  runTimeMedianMs(iMethod) = median(runTimeVal, 'omitnan');
  wallTimeMedianMs(iMethod) = median(wallTimeVal, 'omitnan');
end
aggregateTable = table(methodLabel, methodFamily, satMode, frameMode, modelClass, rateMode, oracleLevel, ...
  numRepeat, resolvedRate, toothHitRate, angleRmseDeg, angleMedianDeg, angleP95Deg, ...
  fdRefAbsMedianHz, fdRefAbsP95Hz, fdRateAbsMedianHzPerSec, fdRateAbsP95HzPerSec, ...
  fdLineRmseMedianHz, nonRefCoherenceFloorMedian, finalObjMedian, runTimeMedianMs, wallTimeMedianMs);
end

function pairCompareTable = localBuildPairCompareTable(methodTable)
%LOCALBUILDPAIRCOMPARETABLE Compare the main wide baseline against oracle methods.

baseLabel = "ms-mf-cp-u-wide";
baseMask = methodTable.methodLabel == baseLabel;
baseTable = methodTable(baseMask, :);
oracleLabelList = unique(methodTable.methodLabel(methodTable.isOracle), 'stable');
numOracle = numel(oracleLabelList);
methodLabel = strings(numOracle, 1);
numAlignedSeed = zeros(numOracle, 1);
centralToothGain = nan(numOracle, 1);
angleRmseGainDeg = nan(numOracle, 1);
angleMedianGainDeg = nan(numOracle, 1);
fdRefAbsMedianGainHz = nan(numOracle, 1);

for iMethod = 1:numOracle
  methodLabel(iMethod) = oracleLabelList(iMethod);
  oracleTable = methodTable(methodTable.methodLabel == oracleLabelList(iMethod), :);
  [seedCommon, idxBase, idxOracle] = intersect(baseTable.taskSeed, oracleTable.taskSeed, 'stable');
  numAlignedSeed(iMethod) = numel(seedCommon);
  if isempty(seedCommon)
    continue;
  end
  baseToothHit = double(baseTable.toothIdx(idxBase) == 0);
  oracleToothHit = double(oracleTable.toothIdx(idxOracle) == 0);
  centralToothGain(iMethod) = mean(oracleToothHit - baseToothHit, 'omitnan');
  baseAngle = baseTable.angleErrDeg(idxBase);
  oracleAngle = oracleTable.angleErrDeg(idxOracle);
  angleRmseGainDeg(iMethod) = sqrt(mean(baseAngle .^ 2, 'omitnan')) - sqrt(mean(oracleAngle .^ 2, 'omitnan'));
  angleMedianGainDeg(iMethod) = median(baseAngle, 'omitnan') - median(oracleAngle, 'omitnan');
  fdRefAbsMedianGainHz(iMethod) = median(abs(baseTable.fdRefErrHz(idxBase)), 'omitnan') - ...
    median(abs(oracleTable.fdRefErrHz(idxOracle)), 'omitnan');
end
pairCompareTable = table(methodLabel, numAlignedSeed, centralToothGain, angleRmseGainDeg, ...
  angleMedianGainDeg, fdRefAbsMedianGainHz);
end

function paperCompareTable = localBuildPaperCompareTable(methodTable)
%LOCALBUILDPAPERCOMPARETABLE Compare MS-MF in-tooth CP-U against paper baselines.

targetLabel = "ms-mf-cp-u-in-tooth";
baselineLabelList = ["ss-sf-static"; "ms-sf-static"; "ss-mf-cp-u-in-tooth"; ...
  "ms-mf-cp-u-wide"; "ms-mf-cp-k-in-tooth"; "ms-mf-cp-u-truth-doa-oracle"];
targetTable = methodTable(methodTable.methodLabel == targetLabel, :);
numBaseline = numel(baselineLabelList);
baselineLabel = strings(numBaseline, 1);
numAlignedSeed = zeros(numBaseline, 1);
angleRmseGainDeg = nan(numBaseline, 1);
angleMedianGainDeg = nan(numBaseline, 1);
fdRefAbsMedianGainHz = nan(numBaseline, 1);
targetBetterAngleRate = nan(numBaseline, 1);

for iBase = 1:numBaseline
  baselineLabel(iBase) = baselineLabelList(iBase);
  baseTable = methodTable(methodTable.methodLabel == baselineLabelList(iBase), :);
  [seedCommon, idxTarget, idxBase] = intersect(targetTable.taskSeed, baseTable.taskSeed, 'stable');
  numAlignedSeed(iBase) = numel(seedCommon);
  if isempty(seedCommon)
    continue;
  end
  targetAngle = targetTable.angleErrDeg(idxTarget);
  baseAngle = baseTable.angleErrDeg(idxBase);
  angleRmseGainDeg(iBase) = sqrt(mean(baseAngle .^ 2, 'omitnan')) - sqrt(mean(targetAngle .^ 2, 'omitnan'));
  angleMedianGainDeg(iBase) = median(baseAngle, 'omitnan') - median(targetAngle, 'omitnan');
  fdRefAbsMedianGainHz(iBase) = median(abs(baseTable.fdRefErrHz(idxBase)), 'omitnan') - ...
    median(abs(targetTable.fdRefErrHz(idxTarget)), 'omitnan');
  targetBetterAngleRate(iBase) = mean(double(targetAngle < baseAngle), 'omitnan');
end
paperCompareTable = table(repmat(targetLabel, numBaseline, 1), baselineLabel, numAlignedSeed, ...
  angleRmseGainDeg, angleMedianGainDeg, fdRefAbsMedianGainHz, targetBetterAngleRate, ...
  'VariableNames', {'targetLabel', 'baselineLabel', 'numAlignedSeed', 'angleRmseGainDeg', ...
  'angleMedianGainDeg', 'fdRefAbsMedianGainHz', 'targetBetterAngleRate'});
end

function rangeTable = localBuildRangeTable(repeatCell)
%LOCALBUILDRANGETABLE Summarize the oracle fd / fdRate boxes used per repeat.

numRepeat = numel(repeatCell);
taskSeed = nan(numRepeat, 1);
toothStepHz = nan(numRepeat, 1);
truthFdRefHz = nan(numRepeat, 1);
truthFdRateHzPerSec = nan(numRepeat, 1);
fdRangeLowerHz = nan(numRepeat, 1);
fdRangeUpperHz = nan(numRepeat, 1);
fdRateRangeLowerHzPerSec = nan(numRepeat, 1);
fdRateRangeUpperHzPerSec = nan(numRepeat, 1);
for iRepeat = 1:numRepeat
  item = repeatCell{iRepeat};
  taskSeed(iRepeat) = item.taskSeed;
  toothStepHz(iRepeat) = item.toothStepHz;
  truthFdRefHz(iRepeat) = item.truthFdRefHz;
  truthFdRateHzPerSec(iRepeat) = item.truthFdRateHzPerSec;
  fdRangeLowerHz(iRepeat) = item.fdRangeOracle(1);
  fdRangeUpperHz(iRepeat) = item.fdRangeOracle(2);
  fdRateRangeLowerHzPerSec(iRepeat) = item.fdRateRangeOracle(1);
  fdRateRangeUpperHzPerSec(iRepeat) = item.fdRateRangeOracle(2);
end
rangeTable = table(taskSeed, toothStepHz, truthFdRefHz, truthFdRateHzPerSec, ...
  fdRangeLowerHz, fdRangeUpperHz, fdRateRangeLowerHzPerSec, fdRateRangeUpperHzPerSec);
end


function timingTable = localBuildTimingTable(repeatCell)
%LOCALBUILDTIMINGTABLE Stack per-repeat wall-clock timing summaries.

numRepeat = numel(repeatCell);
taskSeed = nan(numRepeat, 1);
buildDataMs = nan(numRepeat, 1);
staticBundleMs = nan(numRepeat, 1);
dynamicMethodTotalMs = nan(numRepeat, 1);
dynamicMethodMaxMs = nan(numRepeat, 1);
repeatTotalMs = nan(numRepeat, 1);
staticBundleShare = nan(numRepeat, 1);
dynamicMethodShare = nan(numRepeat, 1);
for iRepeat = 1:numRepeat
  item = repeatCell{iRepeat};
  taskSeed(iRepeat) = item.taskSeed;
  timing = localGetFieldOrDefault(item, 'timing', struct());
  buildDataMs(iRepeat) = localGetFieldOrDefault(timing, 'buildDataMs', NaN);
  staticBundleMs(iRepeat) = localGetFieldOrDefault(timing, 'staticBundleMs', NaN);
  dynamicMethodTotalMs(iRepeat) = localGetFieldOrDefault(timing, 'dynamicMethodTotalMs', NaN);
  dynamicMethodMaxMs(iRepeat) = localGetFieldOrDefault(timing, 'dynamicMethodMaxMs', NaN);
  repeatTotalMs(iRepeat) = localGetFieldOrDefault(timing, 'repeatTotalMs', NaN);
  if isfinite(repeatTotalMs(iRepeat)) && repeatTotalMs(iRepeat) > 0
    staticBundleShare(iRepeat) = staticBundleMs(iRepeat) / repeatTotalMs(iRepeat);
    dynamicMethodShare(iRepeat) = dynamicMethodTotalMs(iRepeat) / repeatTotalMs(iRepeat);
  end
end
timingTable = table(taskSeed, buildDataMs, staticBundleMs, dynamicMethodTotalMs, ...
  dynamicMethodMaxMs, repeatTotalMs, staticBundleShare, dynamicMethodShare);
end

function timingAggregateTable = localBuildTimingAggregateTable(timingTable)
%LOCALBUILDTIMINGAGGREGATETABLE Build compact timing summary across repeats.

if isempty(timingTable)
  timingAggregateTable = table();
  return;
end
stageLabel = ["build-repeat-data"; "static-transition-bundle"; "dynamic-methods-total"; ...
  "slowest-dynamic-method"; "repeat-total"];
valueCell = {timingTable.buildDataMs; timingTable.staticBundleMs; timingTable.dynamicMethodTotalMs; ...
  timingTable.dynamicMethodMaxMs; timingTable.repeatTotalMs};
shareCell = {NaN(height(timingTable), 1); timingTable.staticBundleShare; timingTable.dynamicMethodShare; ...
  timingTable.dynamicMethodMaxMs ./ timingTable.repeatTotalMs; ones(height(timingTable), 1)};
numStage = numel(stageLabel);
medianMs = nan(numStage, 1);
p95Ms = nan(numStage, 1);
totalSec = nan(numStage, 1);
medianRepeatShare = nan(numStage, 1);
for iStage = 1:numStage
  value = valueCell{iStage};
  share = shareCell{iStage};
  medianMs(iStage) = median(value, 'omitnan');
  p95Ms(iStage) = localPercentile(value, 95);
  totalSec(iStage) = sum(value, 'omitnan') / 1000;
  medianRepeatShare(iStage) = median(share, 'omitnan');
end
timingAggregateTable = table(stageLabel, medianMs, p95Ms, totalSec, medianRepeatShare);
end

function representative = localSelectRepresentative(methodTable, repeatCell)
%LOCALSELECTREPRESENTATIVE Pick one repeat where oracle helps most or fails most.

baseLabel = "ms-mf-cp-u-wide";
oracleLabel = "ms-mf-cp-u-in-tooth";
baseTable = methodTable(methodTable.methodLabel == baseLabel, :);
oracleTable = methodTable(methodTable.methodLabel == oracleLabel, :);
[seedCommon, idxBase, idxOracle] = intersect(baseTable.taskSeed, oracleTable.taskSeed, 'stable');
representative = struct('taskSeed', NaN, 'methodLabel', "", 'angleErrDeg', NaN, 'fdRefErrHz', NaN, ...
  'periodicWideToothIdx', NaN, 'oracleToothIdx', NaN, 'methodRows', table());
if isempty(seedCommon)
  return;
end
score = abs(baseTable.toothIdx(idxBase)) - abs(oracleTable.toothIdx(idxOracle));
score = score + (baseTable.angleErrDeg(idxBase) - oracleTable.angleErrDeg(idxOracle));
score(~isfinite(score)) = -inf;
[~, bestRel] = max(score);
seedUse = seedCommon(bestRel);
repeatIdx = find(arrayfun(@(x) repeatCell{x}.taskSeed == seedUse, 1:numel(repeatCell)), 1, 'first');
methodRows = methodTable(methodTable.taskSeed == seedUse, :);
[~, bestMethodRel] = min(methodRows.angleErrDeg);
representative.taskSeed = seedUse;
representative.methodLabel = methodRows.methodLabel(bestMethodRel);
representative.angleErrDeg = methodRows.angleErrDeg(bestMethodRel);
representative.fdRefErrHz = methodRows.fdRefErrHz(bestMethodRel);
representative.periodicWideToothIdx = baseTable.toothIdx(idxBase(bestRel));
representative.oracleToothIdx = oracleTable.toothIdx(idxOracle(bestRel));
representative.methodRows = methodRows;
if ~isempty(repeatIdx)
  representative.toothStepHz = repeatCell{repeatIdx}.toothStepHz;
end
end

function contextSummary = localBuildContextSummary(context, methodList)
%LOCALBUILDCONTEXTSUMMARY Keep only lightweight context metadata.

contextSummary = struct();
contextSummary.frameIntvlSec = context.frameIntvlSec;
contextSummary.periodicOffsetIdx = reshape(context.periodicOffsetIdx, 1, []);
contextSummary.masterOffsetIdx = reshape(context.masterOffsetIdx, 1, []);
contextSummary.selectedSatIdxGlobal = reshape(context.selectedSatIdxGlobal, 1, []);
contextSummary.refSatIdxGlobal = context.refSatIdxGlobal;
contextSummary.otherSatIdxGlobal = context.otherSatIdxGlobal;
contextSummary.baseSeed = context.baseSeed;
contextSummary.tleFileName = string(context.tleFileName);
contextSummary.methodLabelList = reshape([methodList.label], [], 1);
end

function methodListOut = localStripMethodList(methodList)
%LOCALSTRIPMETHODLIST Store method descriptors without function handles or heavy fields.

methodListOut = methodList;
end

function singleMultiCompareTable = localBuildSingleMultiCompareTable(methodTable)
%LOCALBUILDSINGLEMULTICOMPARETABLE Build per-seed SS-MF versus MS-MF in-tooth diagnostics.

singleLabel = "ss-mf-cp-u-in-tooth";
multiLabel = "ms-mf-cp-u-in-tooth";
truthDoaLabel = "ms-mf-cp-u-truth-doa-oracle";
singleTable = methodTable(methodTable.methodLabel == singleLabel, :);
multiTable = methodTable(methodTable.methodLabel == multiLabel, :);
truthDoaTable = methodTable(methodTable.methodLabel == truthDoaLabel, :);
[seedCommon, idxSingle, idxMulti] = intersect(singleTable.taskSeed, multiTable.taskSeed, 'stable');
numSeed = numel(seedCommon);
taskSeed = reshape(seedCommon, [], 1);
angleErrDegSingle = nan(numSeed, 1);
angleErrDegMulti = nan(numSeed, 1);
angleErrDegTruthDoa = nan(numSeed, 1);
angleGainDeg = nan(numSeed, 1);
multiGapToTruthDoaDeg = nan(numSeed, 1);
multiBetterAngle = false(numSeed, 1);
fdRefAbsErrSingleHz = nan(numSeed, 1);
fdRefAbsErrMultiHz = nan(numSeed, 1);
fdRateAbsErrSingleHzPerSec = nan(numSeed, 1);
fdRateAbsErrMultiHzPerSec = nan(numSeed, 1);
singleToothIdx = nan(numSeed, 1);
multiToothIdx = nan(numSeed, 1);
multiToothResidualHz = nan(numSeed, 1);
multiNonRefCoherenceFloor = nan(numSeed, 1);
multiFinalObj = nan(numSeed, 1);

for iSeed = 1:numSeed
  singleRow = singleTable(idxSingle(iSeed), :);
  multiRow = multiTable(idxMulti(iSeed), :);
  angleErrDegSingle(iSeed) = singleRow.angleErrDeg;
  angleErrDegMulti(iSeed) = multiRow.angleErrDeg;
  angleGainDeg(iSeed) = angleErrDegSingle(iSeed) - angleErrDegMulti(iSeed);
  multiBetterAngle(iSeed) = angleGainDeg(iSeed) > 0;
  fdRefAbsErrSingleHz(iSeed) = abs(singleRow.fdRefErrHz);
  fdRefAbsErrMultiHz(iSeed) = abs(multiRow.fdRefErrHz);
  fdRateAbsErrSingleHzPerSec(iSeed) = abs(singleRow.fdRateErrHzPerSec);
  fdRateAbsErrMultiHzPerSec(iSeed) = abs(multiRow.fdRateErrHzPerSec);
  singleToothIdx(iSeed) = singleRow.toothIdx;
  multiToothIdx(iSeed) = multiRow.toothIdx;
  multiToothResidualHz(iSeed) = multiRow.toothResidualHz;
  multiNonRefCoherenceFloor(iSeed) = multiRow.nonRefCoherenceFloor;
  multiFinalObj(iSeed) = multiRow.finalObj;
  idxTruth = find(truthDoaTable.taskSeed == taskSeed(iSeed), 1, 'first');
  if ~isempty(idxTruth)
    angleErrDegTruthDoa(iSeed) = truthDoaTable.angleErrDeg(idxTruth);
    multiGapToTruthDoaDeg(iSeed) = angleErrDegMulti(iSeed) - angleErrDegTruthDoa(iSeed);
  end
end
singleMultiCompareTable = table(taskSeed, angleErrDegSingle, angleErrDegMulti, angleErrDegTruthDoa, ...
  angleGainDeg, multiGapToTruthDoaDeg, multiBetterAngle, fdRefAbsErrSingleHz, fdRefAbsErrMultiHz, ...
  fdRateAbsErrSingleHzPerSec, fdRateAbsErrMultiHzPerSec, singleToothIdx, multiToothIdx, ...
  multiToothResidualHz, multiNonRefCoherenceFloor, multiFinalObj);
end

function tailCaseTable = localBuildTailCaseTable(singleMultiCompareTable)
%LOCALBUILDTAILCASETABLE Keep the worst per-seed single-vs-multi tail cases.

if isempty(singleMultiCompareTable)
  tailCaseTable = singleMultiCompareTable;
  return;
end
[~, orderIdx] = sort(singleMultiCompareTable.angleGainDeg, 'ascend');
numTail = min(8, numel(orderIdx));
tailCaseTable = singleMultiCompareTable(orderIdx(1:numTail), :);
tailRank = (1:height(tailCaseTable)).';
tailCaseTable = addvars(tailCaseTable, tailRank, 'Before', 1);
end

function plotData = localBuildPlotData(methodTable, singleMultiCompareTable, histogramBinCount)
%LOCALBUILDPLOTDATA Pack lightweight distribution data used by the plotting section.

coreMethodLabelList = ["ss-sf-static"; "ms-sf-static"; "ss-mf-cp-u-in-tooth"; ...
  "ms-mf-cp-u-in-tooth"; "ms-mf-cp-u-truth-doa-oracle"];
coreMask = ismember(methodTable.methodLabel, coreMethodLabelList);
plotData = struct();
plotData.histogramBinCount = histogramBinCount;
plotData.coreMethodLabelList = coreMethodLabelList;
plotData.distributionTable = methodTable(coreMask, {'taskSeed', 'methodLabel', 'angleErrDeg', ...
  'fdRefErrHz', 'fdRateErrHzPerSec', 'toothResidualHz'});
plotData.singleMultiCompareTable = singleMultiCompareTable;
end


function localDispTablePreview(dataTable, edgeCount)
%LOCALDISPTABLEPREVIEW Display short first/last preview for long tables.

if nargin < 2 || isempty(edgeCount)
  edgeCount = 4;
end
if isempty(dataTable) || height(dataTable) <= 2 * edgeCount
  disp(dataTable);
  return;
end
fprintf('  showing first %d and last %d of %d rows\n', edgeCount, edgeCount, height(dataTable));
disp(dataTable(1:edgeCount, :));
fprintf('  ...\n');
disp(dataTable((height(dataTable) - edgeCount + 1):height(dataTable), :));
end

function localPlotReplay(plotData)
%LOCALPLOTREPLAY Draw compact CDF and tail diagnostics for the replay.

if isempty(plotData) || ~isstruct(plotData)
  return;
end
histogramBinCount = localGetFieldOrDefault(plotData, 'histogramBinCount', 12);
distributionTable = localGetFieldOrDefault(plotData, 'distributionTable', table());
methodLabelList = localGetFieldOrDefault(plotData, 'coreMethodLabelList', strings(0, 1));
if ~isempty(distributionTable)
  figure('Name', 'In-tooth oracle empirical CDF');
  tiledlayout(3, 1, 'TileSpacing', 'compact');
  nexttile;
  localPlotGroupedEmpiricalCdf(distributionTable, methodLabelList, 'angleErrDeg', true, true);
  xlabel('angle error (deg)');
  ylabel('empirical CDF');
  title('Angle-error CDF');
  grid on;
  nexttile;
  localPlotGroupedEmpiricalCdf(distributionTable, methodLabelList, 'fdRefErrHz', true, true);
  xlabel('|fdRef error| (Hz)');
  ylabel('empirical CDF');
  title('Reference-Doppler center-deviation CDF');
  grid on;
  nexttile;
  localPlotGroupedEmpiricalCdf(distributionTable, methodLabelList, 'fdRateErrHzPerSec', true, true);
  xlabel('|fdRate error| (Hz/s)');
  ylabel('empirical CDF');
  title('Doppler-rate center-deviation CDF');
  grid on;
end

singleMultiCompareTable = localGetFieldOrDefault(plotData, 'singleMultiCompareTable', table());
if ~isempty(singleMultiCompareTable)
  figure('Name', 'Single-vs-multi in-tooth diagnostics');
  tiledlayout(2, 1, 'TileSpacing', 'compact');
  nexttile;
  gainVal = singleMultiCompareTable.angleGainDeg;
  gainVal = gainVal(isfinite(gainVal));
  if ~isempty(gainVal)
    histogram(gainVal, histogramBinCount, 'Normalization', 'probability');
  end
  xline(0, '--');
  xlabel('angle gain: SS-MF - MS-MF (deg)');
  ylabel('probability');
  title('MS-MF gain distribution');
  grid on;
  nexttile;
  localPlotSingleMultiScatter(singleMultiCompareTable);
  xlabel('SS-MF in-tooth angle error (deg)');
  ylabel('MS-MF in-tooth angle error (deg)');
  title('Paired seed scatter');
  grid on;
end
end

function localPlotGroupedEmpiricalCdf(dataTable, methodLabelList, fieldName, useAbsValue, useLogX)
%LOCALPLOTGROUPEDEMPIRICALCDF Plot one empirical CDF curve per method.

if nargin < 4
  useAbsValue = false;
end
if nargin < 5
  useLogX = false;
end
hold on;
legendLabel = strings(0, 1);
for iMethod = 1:numel(methodLabelList)
  methodLabel = methodLabelList(iMethod);
  mask = dataTable.methodLabel == methodLabel;
  if ~any(mask)
    continue;
  end
  value = dataTable.(fieldName)(mask);
  value = reshape(value, [], 1);
  if useAbsValue
    value = abs(value);
  end
  value = value(isfinite(value));
  if isempty(value)
    continue;
  end
  if useLogX
    value = localFloorForLogAxis(value);
  end
  value = sort(value);
  prob = (1:numel(value)).' ./ numel(value);
  stairs(value, prob, 'LineWidth', 1.2);
  legendLabel(end + 1, 1) = localShortMethodLabel(methodLabel); %#ok<AGROW>
end
hold off;
if useLogX
  set(gca, 'XScale', 'log');
end
ylim([0, 1]);
if ~isempty(legendLabel)
  legend(cellstr(legendLabel), 'Location', 'best', 'Interpreter', 'none');
end
end

function localPlotSingleMultiScatter(singleMultiCompareTable)
%LOCALPLOTSINGLEMULTISCATTER Plot paired SS-MF and MS-MF errors with a y=x guide.

xVal = singleMultiCompareTable.angleErrDegSingle;
yVal = singleMultiCompareTable.angleErrDegMulti;
validMask = isfinite(xVal) & isfinite(yVal);
if any(validMask)
  plot(xVal(validMask), yVal(validMask), 'o', 'LineStyle', 'none');
  hold on;
  axisVal = [xVal(validMask); yVal(validMask)];
  axisVal = axisVal(isfinite(axisVal));
  if ~isempty(axisVal)
    limMax = max(axisVal);
    limMin = min(axisVal);
    if limMax > limMin
      pad = 0.05 * (limMax - limMin);
      lim = [max(0, limMin - pad), limMax + pad];
    else
      lim = [0, limMax + eps];
    end
    plot(lim, lim, '--');
    xlim(lim);
    ylim(lim);
  end
  hold off;
end
end

function value = localFloorForLogAxis(value)
%LOCALFLOORFORLOGAXIS Replace zero values by a small positive floor for log plots.

value = abs(value(:));
positiveValue = value(value > 0 & isfinite(value));
if isempty(positiveValue)
  value(:) = eps;
  return;
end
floorValue = min(positiveValue) / 10;
value(value <= 0) = floorValue;
end

function shortLabel = localShortMethodLabel(methodLabel)
%LOCALSHORTMETHODLABEL Keep plot legends compact while preserving method identity.

methodLabel = string(methodLabel);
switch methodLabel
  case "ss-sf-static"
    shortLabel = "SS-SF static";
  case "ms-sf-static"
    shortLabel = "MS-SF static";
  case "ss-mf-cp-u-in-tooth"
    shortLabel = "SS-MF CP-U";
  case "ms-mf-cp-u-in-tooth"
    shortLabel = "MS-MF CP-U";
  case "ms-mf-cp-u-truth-doa-oracle"
    shortLabel = "MS-MF truth-DoA";
  otherwise
    shortLabel = methodLabel;
end
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve the fd tooth spacing from fixture timing.

toothStepHz = NaN;
if isfield(fixture, 'frameIntvlSec') && isfinite(fixture.frameIntvlSec) && fixture.frameIntvlSec > 0
  toothStepHz = 1 / fixture.frameIntvlSec;
  return;
end
if isfield(fixture, 'sceneSeq') && isfield(fixture.sceneSeq, 'timeOffsetSec')
  dt = diff(reshape(fixture.sceneSeq.timeOffsetSec, 1, []));
  dt = dt(isfinite(dt) & dt > 0);
  if ~isempty(dt)
    toothStepHz = 1 / median(dt);
  end
end
if ~(isfinite(toothStepHz) && toothStepHz > 0)
  error('replayMfInToothFdRangeOracle:MissingToothStep', ...
    'Cannot resolve tooth step from fixture timing.');
end
end

function value = localResolveTruthScalar(truth, fieldNameList)
%LOCALRESOLVETRUTHSCALAR Resolve a scalar truth field using ordered fallbacks.

value = NaN;
for iField = 1:numel(fieldNameList)
  fieldName = fieldNameList{iField};
  rawValue = localGetFieldOrDefault(truth, fieldName, NaN);
  if isscalar(rawValue) && isfinite(rawValue)
    value = rawValue;
    return;
  end
end
end

function maxVal = localMaxAbs(value)
%LOCALMAXABS Return max abs for a possibly empty numeric vector.

if isempty(value)
  maxVal = NaN;
else
  maxVal = max(abs(value(:)));
end
end

function percentileValue = localPercentile(valueVec, percentile)
%LOCALPERCENTILE Compute one percentile without requiring Statistics Toolbox.

valueVec = reshape(valueVec, [], 1);
valueVec = valueVec(isfinite(valueVec));
if isempty(valueVec)
  percentileValue = NaN;
  return;
end
valueVec = sort(valueVec);
if numel(valueVec) == 1
  percentileValue = valueVec;
  return;
end
rankPos = 1 + (numel(valueVec) - 1) * percentile / 100;
lo = floor(rankPos);
hi = ceil(rankPos);
if lo == hi
  percentileValue = valueVec(lo);
else
  weightHi = rankPos - lo;
  percentileValue = (1 - weightHi) * valueVec(lo) + weightHi * valueVec(hi);
end
end

function tracker = localCreateProgressTracker(titleText, totalCount, useParfor)
%LOCALCREATEPROGRESSTRACKER Create a client-side progressbar tracker.

tracker = struct('active', false, 'queue', [], 'totalCount', totalCount);
if totalCount <= 0
  return;
end
fprintf('%s\n', titleText);
if exist('progressbar', 'file') ~= 2
  fprintf('[%s] progressbar.m not found; continuing without progress bar.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  tracker.active = true;
  if useParfor
    tracker.queue = parallel.pool.DataQueue;
    afterEach(tracker.queue, @(~) progressbar('advance'));
  end
catch ME
  tracker.active = false;
  tracker.queue = [];
  fprintf('[%s] Progress bar disabled: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), ME.message);
end
end

function localAdvanceProgressTracker(tracker)
%LOCALADVANCEPROGRESSTRACKER Advance the progressbar in a serial loop.

if tracker.active && isempty(tracker.queue)
  progressbar('advance');
end
end

function localCloseProgressTracker(tracker)
%LOCALCLOSEPROGRESSTRACKER Close the progressbar if it is active.

if tracker.active
  pause(0.05);
  progressbar('end');
end
end

function textOut = localModeText(useParfor)
%LOCALMODETEXT Return a short run-mode tag.

if useParfor
  textOut = 'parfor';
else
  textOut = 'serial';
end
end

function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether an outer parfor repeat loop can be used.

tf = false;
try
  tf = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or the supplied default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
