% replayMfInToothTailCaseDiagnose
% Purpose: replay oracle in-tooth tail seeds and classify whether each tail is
% a wrong-tooth case, an fd/rate case, or a same-tooth DoA/local-state case.
% This replay deliberately uses only periodic all-frame fixtures and static or
% truth DoA seeds; in-tooth diagnosis does not use curated subset initialization.
% Fixed-DoA variants isolate whether fd/fdRate release alone can recover the
% non-reference coherence on priority tail seeds. Candidate objective probes
% then compare final-centered, wide-centered, implementable-center, and
% truth-line objective probes without solver adoption. Rescue-bank summaries
% compare disabled, wide-centered, single-MF-centered, combined wide+single-MF,
% and gated wide+single-MF implementable candidates on hard and negative seeds.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; saveSnapshot=true saves only
% replayData via saveExpSnapshot.

clear; close all; clc;

%% Replay configuration
replayName = "replayMfInToothTailCaseDiagnose";

snrDb = 10;                         % Snapshot SNR used to generate rx signals.
seedList = [277; 283; 298; 256; 293; 280; 268; 284]; % Tail and negative-check seeds exposed by the source in-tooth oracle replay.
contextBaseSeed = 253;              % Keep the source oracle context fixed while replaying tail seeds.
saveSnapshot = true;                % true saves the lightweight replayData only.
optVerbose = false;                 % true enables compact estimator / branch trace.
oracleFdHalfToothFraction = 0.49;   % Half-width fraction of the 1/T_f tooth step.
oracleFdRateHalfWidthHzPerSec = 1000; % Local truth-centered fdRate half-width for oracle-rate variants.
staticLocalDoaHalfWidthDeg = [0.002; 0.002]; % Local DoA box around the static seed.
staticWideDoaHalfWidthDeg = [0.010; 0.010];  % Wider same-tooth DoA box for stress comparison.
truthLocalDoaHalfWidthDeg = [0.002; 0.002];  % Local DoA box around the truth DoA seed.
fdRefHealthyAbsHz = 1;              % Same-tooth fdRef health threshold.
fdRateHealthyAbsHzPerSec = 50;      % Same-tooth fdRate health threshold.
coherenceCollapseThreshold = 0.95;  % Non-reference coherence floor below this value is a collapse.
truthDoaGapThresholdDeg = 1e-3;     % Gap to truth-DoA oracle used to flag DoA/local-state tails.
probeDoaStepDegList = [-0.002; -0.001; 0; 0.001; 0.002]; % Small local DoA offsets for objective probes.
probeFdRefStepHzList = [-0.05; -0.02; 0; 0.02; 0.05];    % Small fdRef offsets for joint local probes.
probeLineAlphaList = [0; 0.25; 0.5; 0.75; 0.85; 0.90; 0.95; 0.98; 0.99; 1]; % Line probe interpolation factors toward truth DoA.
probeCoarseDoaStepDegList = [-0.006; -0.004; -0.002; 0; 0.002; 0.004; 0.006]; % Implementable non-truth center sweep offsets.
probeAngleGainThresholdDeg = 5e-4; % Minimum angle gain used to label a probe as rescuing a tail.
probeDamageThresholdDeg = 5e-4;    % Minimum angle loss used to label easy/fd-not-healthy damage.
gatedRescueCoherenceThreshold = 0.20; % Trigger rescue only when default non-ref coherence has clearly collapsed.
gatedRescuePhaseResidThresholdRad = 1.0; % Additional phase-residual gate for collapse-triggered rescue.

seedList = reshape(double(seedList), [], 1);
numRepeat = numel(seedList);

%% Build context and flow options
config = struct();
config.snrDb = snrDb;
config.seedList = seedList;
config.numRepeat = numRepeat;
config.contextBaseSeed = contextBaseSeed;
config.saveSnapshot = saveSnapshot;
config.optVerbose = optVerbose;
config.oracleFdHalfToothFraction = oracleFdHalfToothFraction;
config.oracleFdRateHalfWidthHzPerSec = oracleFdRateHalfWidthHzPerSec;
config.staticLocalDoaHalfWidthDeg = staticLocalDoaHalfWidthDeg(:);
config.staticWideDoaHalfWidthDeg = staticWideDoaHalfWidthDeg(:);
config.truthLocalDoaHalfWidthDeg = truthLocalDoaHalfWidthDeg(:);
config.fdRefHealthyAbsHz = fdRefHealthyAbsHz;
config.fdRateHealthyAbsHzPerSec = fdRateHealthyAbsHzPerSec;
config.coherenceCollapseThreshold = coherenceCollapseThreshold;
config.truthDoaGapThresholdDeg = truthDoaGapThresholdDeg;
config.probeDoaStepDegList = reshape(double(probeDoaStepDegList), [], 1);
config.probeFdRefStepHzList = reshape(double(probeFdRefStepHzList), [], 1);
config.probeLineAlphaList = reshape(double(probeLineAlphaList), [], 1);
config.probeCoarseDoaStepDegList = reshape(double(probeCoarseDoaStepDegList), [], 1);
config.probeAngleGainThresholdDeg = probeAngleGainThresholdDeg;
config.probeDamageThresholdDeg = probeDamageThresholdDeg;
config.gatedRescueCoherenceThreshold = gatedRescueCoherenceThreshold;
config.gatedRescuePhaseResidThresholdRad = gatedRescuePhaseResidThresholdRad;

fprintf('Running %s ...\n', char(replayName));
fprintf('  tail seeds                       : %s\n', mat2str(config.seedList(:).'));
fprintf('  repeats                          : %d\n', config.numRepeat);
fprintf('  context base seed                : %d\n', config.contextBaseSeed);
fprintf('  snr (dB)                         : %.2f\n', config.snrDb);
fprintf('  repeat mode                      : %s\n', 'parfor-auto');
fprintf('  save snapshot                    : %d\n', config.saveSnapshot);
fprintf('  fd oracle half-tooth fraction    : %.3f\n', config.oracleFdHalfToothFraction);
fprintf('  fdRate oracle half-width         : %.2f Hz/s\n', config.oracleFdRateHalfWidthHzPerSec);
fprintf('  in-tooth DoA initialization      : static/truth only, no curated subset seed\n');
fprintf('  candidate probe DoA steps        : %s deg\n', mat2str(config.probeDoaStepDegList(:).'));
fprintf('  candidate probe fdRef steps      : %s Hz\n', mat2str(config.probeFdRefStepHzList(:).'));
fprintf('  line probe alpha values          : %s\n', mat2str(config.probeLineAlphaList(:).'));
fprintf('  implementable center DoA steps   : %s deg\n', mat2str(config.probeCoarseDoaStepDegList(:).'));
fprintf('  gated rescue coherence threshold : %.3f\n', config.gatedRescueCoherenceThreshold);
fprintf('  gated rescue phase threshold     : %.3f rad\n', config.gatedRescuePhaseResidThresholdRad);

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
contextOpt = struct( ...
  'baseSeed', config.contextBaseSeed, ...
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
context = localDisableSubsetBankForInTooth(context);
fprintf('[%s] Shared dynamic context built.\n', char(datetime('now', 'Format', 'HH:mm:ss')));

%% Run replay batch
repeatCell = cell(numRepeat, 1);
useParfor = localCanUseParfor() && numRepeat > 1;
tracker = localCreateProgressTracker(sprintf('%s repeat batch (%s)', char(replayName), localModeText(useParfor)), numRepeat, useParfor);
try
  if useParfor
    progressQueue = tracker.queue;
    parfor iRepeat = 1:numRepeat
      repeatCell{iRepeat} = localRunOneRepeat(iRepeat, config.seedList(iRepeat), context, flowOpt, methodList, config);
      if ~isempty(progressQueue)
        send(progressQueue, iRepeat);
      end
    end
  else
    for iRepeat = 1:numRepeat
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
tailDiagnosisTable = localBuildTailDiagnosisTable(methodTable, config);
candidateProbeTable = localBuildCandidateProbeTable(repeatCell);
candidateWinnerTable = localBuildCandidateWinnerTable(candidateProbeTable, tailDiagnosisTable, config);
lineProbeSummaryTable = localBuildLineProbeSummaryTable(candidateProbeTable, config);
rescueBankSummaryTable = localBuildRescueBankSummaryTable(candidateProbeTable, tailDiagnosisTable, config);
rescueBankAggregateTable = localBuildRescueBankAggregateTable(rescueBankSummaryTable);
rangeTable = localBuildRangeTable(repeatCell);
timingTable = localBuildTimingTable(repeatCell);
timingAggregateTable = localBuildTimingAggregateTable(timingTable);
identityTable = localBuildIdentityTable(repeatCell);
detailDiagnosticTable = localBuildDetailDiagnosticTable(methodTable, tailDiagnosisTable);
plotData = localBuildPlotData(methodTable, tailDiagnosisTable);
plotData.candidateWinnerTable = candidateWinnerTable;
plotData.lineProbeSummaryTable = lineProbeSummaryTable;
plotData.rescueBankSummaryTable = rescueBankSummaryTable;
plotData.rescueBankAggregateTable = rescueBankAggregateTable;

%% Data storage
replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, methodList);
replayData.methodTable = methodTable;
replayData.aggregateTable = aggregateTable;
replayData.tailDiagnosisTable = tailDiagnosisTable;
replayData.candidateProbeTable = candidateProbeTable;
replayData.candidateWinnerTable = candidateWinnerTable;
replayData.lineProbeSummaryTable = lineProbeSummaryTable;
replayData.rescueBankSummaryTable = rescueBankSummaryTable;
replayData.rescueBankAggregateTable = rescueBankAggregateTable;
replayData.identityTable = identityTable;
replayData.detailDiagnosticTable = detailDiagnosticTable;
replayData.rangeTable = rangeTable;
replayData.timingTable = timingTable;
replayData.timingAggregateTable = timingAggregateTable;
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
tailDiagnosisTable = replayData.tailDiagnosisTable;
if isfield(replayData, 'candidateProbeTable')
  candidateProbeTable = replayData.candidateProbeTable;
else
  candidateProbeTable = table();
end
if isfield(replayData, 'candidateWinnerTable')
  candidateWinnerTable = replayData.candidateWinnerTable;
else
  candidateWinnerTable = localBuildCandidateWinnerTable(candidateProbeTable, tailDiagnosisTable, replayData.config);
end
if isfield(replayData, 'lineProbeSummaryTable')
  lineProbeSummaryTable = replayData.lineProbeSummaryTable;
else
  lineProbeSummaryTable = localBuildLineProbeSummaryTable(candidateProbeTable, replayData.config);
end
if isfield(replayData, 'rescueBankSummaryTable')
  rescueBankSummaryTable = replayData.rescueBankSummaryTable;
else
  rescueBankSummaryTable = localBuildRescueBankSummaryTable(candidateProbeTable, tailDiagnosisTable, replayData.config);
end
if isfield(replayData, 'rescueBankAggregateTable')
  rescueBankAggregateTable = replayData.rescueBankAggregateTable;
else
  rescueBankAggregateTable = localBuildRescueBankAggregateTable(rescueBankSummaryTable);
end
if isfield(replayData, 'identityTable')
  identityTable = replayData.identityTable;
else
  identityTable = table();
end
if isfield(replayData, 'detailDiagnosticTable')
  detailDiagnosticTable = replayData.detailDiagnosticTable;
else
  detailDiagnosticTable = localBuildDetailDiagnosticTable(methodTable, tailDiagnosisTable);
end
rangeTable = replayData.rangeTable;
if isfield(replayData, 'timingAggregateTable')
  timingAggregateTable = replayData.timingAggregateTable;
else
  timingAggregateTable = localBuildTimingAggregateTable(replayData.timingTable);
end

fprintf('Running replayMfInToothTailCaseDiagnose ...\n');
fprintf('  repeats                          : %d\n', numel(unique(methodTable.taskSeed)));
if isfield(replayData.config, 'contextBaseSeed')
  fprintf('  context base seed                : %d\n', replayData.config.contextBaseSeed);
end
fprintf('  methods                          : %d\n', numel(unique(methodTable.methodLabel)));
fprintf('  snr (dB)                         : %.2f\n', replayData.config.snrDb);
fprintf('  in-tooth DoA initialization      : static/truth only, no curated subset seed\n');
if isfield(replayData.config, 'probeDoaStepDegList')
  fprintf('  candidate probe DoA steps        : %s deg\n', mat2str(replayData.config.probeDoaStepDegList(:).'));
end
if isfield(replayData.config, 'probeFdRefStepHzList')
  fprintf('  candidate probe fdRef steps      : %s Hz\n', mat2str(replayData.config.probeFdRefStepHzList(:).'));
end
if isfield(replayData.config, 'probeLineAlphaList')
  fprintf('  line probe alpha values          : %s\n', mat2str(replayData.config.probeLineAlphaList(:).'));
end
if isfield(replayData.config, 'probeCoarseDoaStepDegList')
  fprintf('  implementable center DoA steps   : %s deg\n', mat2str(replayData.config.probeCoarseDoaStepDegList(:).'));
end
if isfield(replayData.config, 'gatedRescueCoherenceThreshold')
  fprintf('  gated rescue coherence threshold : %.3f\n', replayData.config.gatedRescueCoherenceThreshold);
end
if isfield(replayData.config, 'gatedRescuePhaseResidThresholdRad')
  fprintf('  gated rescue phase threshold     : %.3f rad\n', replayData.config.gatedRescuePhaseResidThresholdRad);
end
if ~isempty(identityTable)
  fprintf('\n========== Replay identity check ==========\n');
  localDispTablePreview(identityTable, 4);
end
fprintf('\n========== Tail method aggregate ==========\n');
disp(aggregateTable);
fprintf('\n========== Tail classification ==========\n');
disp(tailDiagnosisTable);
if ~isempty(detailDiagnosticTable)
  fprintf('\n========== Tail per-seed detailed diagnostics ==========\n');
  disp(detailDiagnosticTable);
end
if ~isempty(candidateWinnerTable)
  fprintf('\n========== Tail candidate winner summary ==========\n');
  disp(candidateWinnerTable);
end
if ~isempty(lineProbeSummaryTable)
  fprintf('\n========== Tail line probe summary ==========\n');
  disp(lineProbeSummaryTable);
end
if ~isempty(rescueBankAggregateTable)
  fprintf('\n========== Tail rescue bank aggregate ==========\n');
  disp(rescueBankAggregateTable);
end
if ~isempty(rescueBankSummaryTable)
  fprintf('\n========== Tail rescue bank per-seed summary ==========\n');
  disp(rescueBankSummaryTable);
end
if ~isempty(candidateProbeTable)
  fprintf('\n========== Tail candidate objective probe preview ==========\n');
  localDispTablePreview(candidateProbeTable, 8);
end
if ~isempty(timingAggregateTable)
  fprintf('\n========== Runtime timing summary ==========\n');
  disp(timingAggregateTable);
end
fprintf('\n========== Oracle range summary ==========\n');
localDispTablePreview(rangeTable, 4);
localPlotReplay(replayData.plotData);

%% Local helpers

function context = localDisableSubsetBankForInTooth(context)
%LOCALDISABLESUBSETBANKFORINTOOTH Skip curated and random subset fixtures.

context.subsetOffsetCell = {};
context.subsetLabelList = strings(0, 1);
context.numSubsetRandomTrial = 0;
end

function methodList = localBuildMethodList(config)
%LOCALBUILDMETHODLIST Define compact tail-diagnosis method variants.

methodList = repmat(localMakeMethod("", "", "", false, "", "", [], false, "", ""), 0, 1);
methodList(end + 1, 1) = localMakeMethod("ss-mf-cp-u-in-tooth", "single-baseline", "ref", ...
  false, "static", "truthFreq", config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-in-tooth", "tail-target", "ms", ...
  false, "static", "truthFreq", config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-k-in-tooth", "known-rate-check", "ms", ...
  true, "static", "truthFreq", config.staticLocalDoaHalfWidthDeg, false, "truthHalfTooth", "known");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-static-doa-fixed-in-tooth", "fixed-doa-check", "ms", ...
  false, "static", "truthFreq", [0; 0], true, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-truth-doa-fixed-in-tooth", "truth-fixed-doa-check", "ms", ...
  false, "truth", "truthFreq", [0; 0], true, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-wide-doa-in-tooth", "wide-doa-check", "ms", ...
  false, "static", "truthFreq", config.staticWideDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
methodList(end + 1, 1) = localMakeMethod("ms-mf-cp-u-truth-doa-oracle", "truth-doa-upper", "ms", ...
  false, "truth", "truthFreq", config.truthLocalDoaHalfWidthDeg, false, "truthHalfTooth", "truthLocal");
end

function method = localMakeMethod(label, family, viewMode, isKnownRate, doaSeedMode, initMode, doaHalfWidthDeg, freezeDoa, fdRangeMode, fdRateRangeMode)
%LOCALMAKEMETHOD Construct one method descriptor with a stable field set.

method = struct();
method.label = string(label);
method.family = string(family);
method.viewMode = string(viewMode);
method.isKnownRate = logical(isKnownRate);
method.doaSeedMode = string(doaSeedMode);
method.initMode = string(initMode);
method.doaHalfWidthDeg = reshape(double(doaHalfWidthDeg), [], 1);
method.freezeDoa = logical(freezeDoa);
method.fdRangeMode = string(fdRangeMode);
method.fdRateRangeMode = string(fdRateRangeMode);
end

function repeatOut = localRunOneRepeat(iRepeat, taskSeed, context, flowOpt, methodList, config)
%LOCALRUNONEREPEAT Build one tail repeat and run all in-tooth variants.

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
  error('replayMfInToothTailCaseDiagnose:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for in-tooth tail diagnosis.');
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
identitySummary = localBuildIdentitySummary(taskSeed, truth, toothStepHz, fdRangeOracle, ...
  fdRateRangeOracle, staticBundle, numel(methodList) + 2);

rowCell = cell(numel(methodList) + 2, 1);
rowCell{1} = localBuildSummaryRow(iRepeat, taskSeed, "ss-sf-static", "static-baseline", ...
  "single", "single", false, "static", "static", "default", "none", NaN, false, ...
  staticBundle.caseStaticRefOnly, truth, toothStepHz, NaN);
rowCell{2} = localBuildSummaryRow(iRepeat, taskSeed, "ms-sf-static", "static-baseline", ...
  "multi", "single", false, "static", "static", "default", "none", NaN, false, ...
  staticBundle.caseStaticMs, truth, toothStepHz, NaN);

dynamicMethodWallTimeMs = nan(numel(methodList), 1);
methodCaseCell = cell(numel(methodList), 1);
for iMethod = 1:numel(methodList)
  method = methodList(iMethod);
  tMethod = tic;
  caseUse = localRunDynamicMethod(method, fixture, staticBundle, truth, ...
    context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, config.optVerbose, ...
    flowOpt, fdRangeOracle, fdRateRangeOracle, toothStepHz);
  dynamicMethodWallTimeMs(iMethod) = 1000 * toc(tMethod);
  methodCaseCell{iMethod} = caseUse;
  rowCell{iMethod + 2} = localBuildSummaryRow(iRepeat, taskSeed, method.label, method.family, ...
    localInferSatMode(method.viewMode), "multi", method.isKnownRate, method.doaSeedMode, method.initMode, ...
    method.fdRangeMode, method.fdRateRangeMode, localMaxAbs(method.doaHalfWidthDeg), method.freezeDoa, ...
    caseUse, truth, toothStepHz, dynamicMethodWallTimeMs(iMethod));
end
dynamicMethodTotalMs = sum(dynamicMethodWallTimeMs, 'omitnan');
dynamicMethodMaxMs = max(dynamicMethodWallTimeMs, [], 'omitnan');
candidateProbeTable = localBuildCandidateProbeForRepeat(iRepeat, taskSeed, fixture, staticBundle, truth, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, flowOpt, methodList, ...
  methodCaseCell, fdRangeOracle, fdRateRangeOracle, toothStepHz, config);
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
repeatOut.candidateProbeTable = candidateProbeTable;
repeatOut.identitySummary = identitySummary;
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
%LOCALRUNDYNAMICMETHOD Run one dynamic method descriptor on the periodic view.

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

dynOpt = localBuildDynOptForMethod(method, fixture, flowOpt, seedDoaParam);

switch method.initMode
  case "truthFreq"
    initParam = localBuildTruthFreqInitParam(seedDoaParam, truthFdRefHz, truthFdRateHzPerSec, method.isKnownRate);
    initOverride = struct('startTag', method.label, 'initParam', initParam, ...
      'initDoaParam', seedDoaParam(:), 'initDoaHalfWidth', method.doaHalfWidthDeg(:), ...
      'freezeDoa', method.freezeDoa);
  otherwise
    error('replayMfInToothTailCaseDiagnose:UnknownInitMode', ...
      'Unknown method initMode "%s".', method.initMode);
end

caseUse = runDynamicDoaDopplerCase("tail-" + method.label, localInferSatMode(method.viewMode), viewUse, truth, ...
  pilotWave, carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, optVerbose, ...
  dynOpt, method.isKnownRate, debugTruthUse, initOverride);
caseUse.oracleMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, 'fdRateRangeUse', fdRateRangeUse);
end

function dynOpt = localBuildDynOptForMethod(method, fixture, flowOpt, seedDoaParam)
%LOCALBUILDDYNOPTFORMETHOD Build dynamic estimator options for one replay method.

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = method.freezeDoa;
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
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
    error('replayMfInToothTailCaseDiagnose:UnknownViewMode', ...
      'Unknown method viewMode "%s".', method.viewMode);
end
end

function seedDoaParam = localResolveDoaSeed(method, staticCase, truth)
%LOCALRESOLVEDOASEED Choose the DoA seed without curated subset input.

switch method.doaSeedMode
  case "static"
    seedDoaParam = staticCase.estResult.doaParamEst(:);
  case "truth"
    seedDoaParam = reshape(truth.latlonTrueDeg, [], 1);
  otherwise
    error('replayMfInToothTailCaseDiagnose:UnknownDoaSeedMode', ...
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

function satMode = localInferSatMode(viewMode)
%LOCALINFERSATMODE Convert view mode into the common sat-mode label.

if string(viewMode) == "ref"
  satMode = "single";
else
  satMode = "multi";
end
end

function row = localBuildSummaryRow(iRepeat, taskSeed, methodLabel, methodFamily, satMode, frameMode, isKnownRate, doaSeedMode, initMode, fdRangeMode, fdRateRangeMode, doaHalfWidthMaxDeg, freezeDoa, caseUse, truth, toothStepHz, wallTimeMs)
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
row.rateMode = localInferRateMode(frameMode, isKnownRate);
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

initEval = localGetEvalDiag(caseUse, 'initEval');
finalEval = localGetEvalDiag(caseUse, 'finalEval');
[truthDeltaFdRefHz, truthDeltaFdRateHzPerSec] = localBuildTruthDeltaFdFit(truth);
row.initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
row.initResidualNorm = localGetFieldOrDefault(initEval, 'residualNorm', NaN);
initCohSat = localGetFieldOrDefault(initEval, 'coherenceSat', []);
finalCohSat = localGetFieldOrDefault(finalEval, 'coherenceSat', []);
initResidualSat = localGetFieldOrDefault(initEval, 'residualSat', []);
finalResidualSat = localGetFieldOrDefault(finalEval, 'residualSat', []);
initDeltaFdRef = localGetEvalDeltaFdRef(initEval);
finalDeltaFdRef = localGetEvalDeltaFdRef(finalEval);
initDeltaFdRate = localGetEvalDeltaFdRate(initEval);
finalDeltaFdRate = localGetEvalDeltaFdRate(finalEval);
for iSat = 1:2
  truthDeltaFdRef = localVectorElem(truthDeltaFdRefHz, iSat, NaN);
  truthDeltaFdRate = localVectorElem(truthDeltaFdRateHzPerSec, iSat, NaN);
  row.(sprintf('initCohSat%d', iSat)) = localVectorElem(initCohSat, iSat, NaN);
  row.(sprintf('finalCohSat%d', iSat)) = localVectorElem(finalCohSat, iSat, NaN);
  row.(sprintf('initResidualSat%d', iSat)) = localVectorElem(initResidualSat, iSat, NaN);
  row.(sprintf('finalResidualSat%d', iSat)) = localVectorElem(finalResidualSat, iSat, NaN);
  row.(sprintf('initDeltaFdRefErrSat%dHz', iSat)) = ...
    localVectorElem(initDeltaFdRef, iSat, NaN) - truthDeltaFdRef;
  row.(sprintf('finalDeltaFdRefErrSat%dHz', iSat)) = ...
    localVectorElem(finalDeltaFdRef, iSat, NaN) - truthDeltaFdRef;
  row.(sprintf('initDeltaFdRateErrSat%dHzPerSec', iSat)) = ...
    localVectorElem(initDeltaFdRate, iSat, NaN) - truthDeltaFdRate;
  row.(sprintf('finalDeltaFdRateErrSat%dHzPerSec', iSat)) = ...
    localVectorElem(finalDeltaFdRate, iSat, NaN) - truthDeltaFdRate;
end
end

function rateMode = localInferRateMode(frameMode, isKnownRate)
%LOCALINFERRATEMODE Return a compact rate-mode label for summary tables.

if string(frameMode) == "single"
  rateMode = "static";
elseif isKnownRate
  rateMode = "known-rate";
else
  rateMode = "unknown-rate";
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

function identitySummary = localBuildIdentitySummary(taskSeed, truth, toothStepHz, fdRangeOracle, fdRateRangeOracle, staticBundle, methodCount)
%LOCALBUILDIDENTITYSUMMARY Capture seed identity fields for reproducibility checks.

truthDoaDeg = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', [NaN; NaN]), [], 1);
staticRefDoaDeg = localResolveCaseDoa(staticBundle.caseStaticRefOnly);
staticMsDoaDeg = localResolveCaseDoa(staticBundle.caseStaticMs);
staticRefInfo = summarizeDoaDopplerCase(staticBundle.caseStaticRefOnly, truth);
staticMsInfo = summarizeDoaDopplerCase(staticBundle.caseStaticMs, truth);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
identitySummary = struct();
identitySummary.taskSeed = taskSeed;
identitySummary.truthDoa1Deg = localVectorElem(truthDoaDeg, 1, NaN);
identitySummary.truthDoa2Deg = localVectorElem(truthDoaDeg, 2, NaN);
identitySummary.staticRefDoa1Deg = localVectorElem(staticRefDoaDeg, 1, NaN);
identitySummary.staticRefDoa2Deg = localVectorElem(staticRefDoaDeg, 2, NaN);
identitySummary.staticMsDoa1Deg = localVectorElem(staticMsDoaDeg, 1, NaN);
identitySummary.staticMsDoa2Deg = localVectorElem(staticMsDoaDeg, 2, NaN);
identitySummary.staticRefAngleErrDeg = localGetFieldOrDefault(staticRefInfo, 'angleErrDeg', NaN);
identitySummary.staticMsAngleErrDeg = localGetFieldOrDefault(staticMsInfo, 'angleErrDeg', NaN);
identitySummary.truthFdRefHz = truthFdRefHz;
identitySummary.truthFdRateHzPerSec = truthFdRateHzPerSec;
identitySummary.toothStepHz = toothStepHz;
identitySummary.fdRangeLowerHz = fdRangeOracle(1);
identitySummary.fdRangeUpperHz = fdRangeOracle(2);
identitySummary.fdRateRangeLowerHzPerSec = fdRateRangeOracle(1);
identitySummary.fdRateRangeUpperHzPerSec = fdRateRangeOracle(2);
identitySummary.methodCount = methodCount;
identitySummary.identityKey = string(sprintf('seed=%d|truthDoa=%.12g,%.12g|staticMs=%.12g,%.12g|fd=%.12g|rate=%.12g', ...
  taskSeed, identitySummary.truthDoa1Deg, identitySummary.truthDoa2Deg, ...
  identitySummary.staticMsDoa1Deg, identitySummary.staticMsDoa2Deg, ...
  identitySummary.truthFdRefHz, identitySummary.truthFdRateHzPerSec));
end

function identityTable = localBuildIdentityTable(repeatCell)
%LOCALBUILDIDENTITYTABLE Stack per-repeat identity summaries into one table.

identityCell = cell(numel(repeatCell), 1);
for iRepeat = 1:numel(repeatCell)
  identityCell{iRepeat} = repeatCell{iRepeat}.identitySummary;
end
identityTable = struct2table([identityCell{:}].');
end

function doaParam = localResolveCaseDoa(caseUse)
%LOCALRESOLVECASEDOA Read estimated DoA parameters from one case result.

doaParam = [NaN; NaN];
if ~isstruct(caseUse) || ~isfield(caseUse, 'estResult') || ~isstruct(caseUse.estResult)
  return;
end
rawDoa = localGetFieldOrDefault(caseUse.estResult, 'doaParamEst', doaParam);
rawDoa = reshape(rawDoa, [], 1);
if numel(rawDoa) >= 2
  doaParam = rawDoa(1:2);
end
end

function aggregateTable = localBuildAggregateTable(methodTable)
%LOCALBUILDAGGREGATETABLE Build compact method-level distribution metrics.

methodList = unique(methodTable.methodLabel, 'stable');
numMethod = numel(methodList);
methodLabel = strings(numMethod, 1);
methodFamily = strings(numMethod, 1);
satMode = strings(numMethod, 1);
rateMode = strings(numMethod, 1);
numRepeat = zeros(numMethod, 1);
resolvedRate = nan(numMethod, 1);
toothHitRate = nan(numMethod, 1);
angleRmseDeg = nan(numMethod, 1);
angleMedianDeg = nan(numMethod, 1);
angleMaxDeg = nan(numMethod, 1);
fdRefAbsMedianHz = nan(numMethod, 1);
fdRateAbsMedianHzPerSec = nan(numMethod, 1);
nonRefCoherenceFloorMedian = nan(numMethod, 1);
wallTimeMedianMs = nan(numMethod, 1);

for iMethod = 1:numMethod
  methodLabel(iMethod) = methodList(iMethod);
  mask = methodTable.methodLabel == methodList(iMethod);
  idxFirst = find(mask, 1, 'first');
  methodFamily(iMethod) = methodTable.methodFamily(idxFirst);
  satMode(iMethod) = methodTable.satMode(idxFirst);
  rateMode(iMethod) = methodTable.rateMode(idxFirst);
  numRepeat(iMethod) = nnz(mask);
  resolvedRate(iMethod) = mean(double(methodTable.isResolved(mask)), 'omitnan');
  toothHitRate(iMethod) = mean(double(methodTable.toothIdx(mask) == 0), 'omitnan');
  angleVal = methodTable.angleErrDeg(mask);
  angleRmseDeg(iMethod) = sqrt(mean(angleVal .^ 2, 'omitnan'));
  angleMedianDeg(iMethod) = median(angleVal, 'omitnan');
  angleMaxDeg(iMethod) = max(angleVal, [], 'omitnan');
  fdRefAbsMedianHz(iMethod) = median(abs(methodTable.fdRefErrHz(mask)), 'omitnan');
  fdRateAbsMedianHzPerSec(iMethod) = median(abs(methodTable.fdRateErrHzPerSec(mask)), 'omitnan');
  nonRefCoherenceFloorMedian(iMethod) = median(methodTable.nonRefCoherenceFloor(mask), 'omitnan');
  wallTimeMedianMs(iMethod) = median(methodTable.wallTimeMs(mask), 'omitnan');
end
aggregateTable = table(methodLabel, methodFamily, satMode, rateMode, numRepeat, resolvedRate, toothHitRate, ...
  angleRmseDeg, angleMedianDeg, angleMaxDeg, fdRefAbsMedianHz, fdRateAbsMedianHzPerSec, ...
  nonRefCoherenceFloorMedian, wallTimeMedianMs);
end

function tailDiagnosisTable = localBuildTailDiagnosisTable(methodTable, config)
%LOCALBUILDTAILDIAGNOSISTABLE Classify each tail seed using in-tooth diagnostics.

seedList = unique(methodTable.taskSeed, 'stable');
numSeed = numel(seedList);
taskSeed = reshape(seedList, [], 1);
staticAngleErrDeg = nan(numSeed, 1);
singleAngleErrDeg = nan(numSeed, 1);
multiAngleErrDeg = nan(numSeed, 1);
wideDoaAngleErrDeg = nan(numSeed, 1);
knownAngleErrDeg = nan(numSeed, 1);
staticDoaFixedAngleErrDeg = nan(numSeed, 1);
truthDoaFixedAngleErrDeg = nan(numSeed, 1);
truthDoaAngleErrDeg = nan(numSeed, 1);
angleGainSingleMinusMultiDeg = nan(numSeed, 1);
multiGapToTruthDoaDeg = nan(numSeed, 1);
wideDoaGainDeg = nan(numSeed, 1);
knownGainDeg = nan(numSeed, 1);
multiToothIdx = nan(numSeed, 1);
multiToothResidualHz = nan(numSeed, 1);
multiFdRefAbsErrHz = nan(numSeed, 1);
multiFdRateAbsErrHzPerSec = nan(numSeed, 1);
multiFdLineRmseHz = nan(numSeed, 1);
multiNonRefCoherenceFloor = nan(numSeed, 1);
multiNonRefRmsPhaseResidRad = nan(numSeed, 1);
tailClass = strings(numSeed, 1);

for iSeed = 1:numSeed
  seedNow = taskSeed(iSeed);
  staticRow = localFindMethodRow(methodTable, seedNow, "ms-sf-static");
  singleRow = localFindMethodRow(methodTable, seedNow, "ss-mf-cp-u-in-tooth");
  multiRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-u-in-tooth");
  wideRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-u-wide-doa-in-tooth");
  knownRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-k-in-tooth");
  staticFixedRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-u-static-doa-fixed-in-tooth");
  truthFixedRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-u-truth-doa-fixed-in-tooth");
  truthDoaRow = localFindMethodRow(methodTable, seedNow, "ms-mf-cp-u-truth-doa-oracle");
  staticAngleErrDeg(iSeed) = localTableValue(staticRow, 'angleErrDeg', NaN);
  singleAngleErrDeg(iSeed) = localTableValue(singleRow, 'angleErrDeg', NaN);
  multiAngleErrDeg(iSeed) = localTableValue(multiRow, 'angleErrDeg', NaN);
  wideDoaAngleErrDeg(iSeed) = localTableValue(wideRow, 'angleErrDeg', NaN);
  knownAngleErrDeg(iSeed) = localTableValue(knownRow, 'angleErrDeg', NaN);
  staticDoaFixedAngleErrDeg(iSeed) = localTableValue(staticFixedRow, 'angleErrDeg', NaN);
  truthDoaFixedAngleErrDeg(iSeed) = localTableValue(truthFixedRow, 'angleErrDeg', NaN);
  truthDoaAngleErrDeg(iSeed) = localTableValue(truthDoaRow, 'angleErrDeg', NaN);
  angleGainSingleMinusMultiDeg(iSeed) = singleAngleErrDeg(iSeed) - multiAngleErrDeg(iSeed);
  multiGapToTruthDoaDeg(iSeed) = multiAngleErrDeg(iSeed) - truthDoaAngleErrDeg(iSeed);
  wideDoaGainDeg(iSeed) = multiAngleErrDeg(iSeed) - wideDoaAngleErrDeg(iSeed);
  knownGainDeg(iSeed) = multiAngleErrDeg(iSeed) - knownAngleErrDeg(iSeed);
  multiToothIdx(iSeed) = localTableValue(multiRow, 'toothIdx', NaN);
  multiToothResidualHz(iSeed) = localTableValue(multiRow, 'toothResidualHz', NaN);
  multiFdRefAbsErrHz(iSeed) = abs(localTableValue(multiRow, 'fdRefErrHz', NaN));
  multiFdRateAbsErrHzPerSec(iSeed) = abs(localTableValue(multiRow, 'fdRateErrHzPerSec', NaN));
  multiFdLineRmseHz(iSeed) = localTableValue(multiRow, 'fdLineRmseHz', NaN);
  multiNonRefCoherenceFloor(iSeed) = localTableValue(multiRow, 'nonRefCoherenceFloor', NaN);
  multiNonRefRmsPhaseResidRad(iSeed) = localTableValue(multiRow, 'nonRefRmsPhaseResidRad', NaN);
  tailClass(iSeed) = localClassifyTail(multiToothIdx(iSeed), multiFdRefAbsErrHz(iSeed), ...
    multiFdRateAbsErrHzPerSec(iSeed), multiNonRefCoherenceFloor(iSeed), ...
    multiGapToTruthDoaDeg(iSeed), config);
end

tailDiagnosisTable = table(taskSeed, staticAngleErrDeg, singleAngleErrDeg, multiAngleErrDeg, ...
  wideDoaAngleErrDeg, knownAngleErrDeg, staticDoaFixedAngleErrDeg, truthDoaFixedAngleErrDeg, ...
  truthDoaAngleErrDeg, angleGainSingleMinusMultiDeg, ...
  multiGapToTruthDoaDeg, wideDoaGainDeg, knownGainDeg, multiToothIdx, multiToothResidualHz, ...
  multiFdRefAbsErrHz, multiFdRateAbsErrHzPerSec, multiFdLineRmseHz, multiNonRefCoherenceFloor, ...
  multiNonRefRmsPhaseResidRad, tailClass);
end

function row = localFindMethodRow(methodTable, taskSeed, methodLabel)
%LOCALFINDMETHODROW Return one method row for one seed or an empty table.

mask = methodTable.taskSeed == taskSeed & methodTable.methodLabel == string(methodLabel);
row = methodTable(mask, :);
if height(row) > 1
  row = row(1, :);
end
end

function value = localTableValue(row, fieldName, defaultValue)
%LOCALTABLEVALUE Read a scalar table value with fallback.

value = defaultValue;
if isempty(row) || height(row) == 0 || ~ismember(fieldName, row.Properties.VariableNames)
  return;
end
rawValue = row.(fieldName);
if ~isempty(rawValue)
  value = rawValue(1);
end
end

function tailClass = localClassifyTail(toothIdx, fdRefAbsErrHz, fdRateAbsErrHzPerSec, nonRefCoherenceFloor, multiGapToTruthDoaDeg, config)
%LOCALCLASSIFYTAIL Assign one stable text class to a tail seed.

if ~(isfinite(toothIdx) && toothIdx == 0)
  tailClass = "wrong-tooth";
  return;
end
fdHealthy = isfinite(fdRefAbsErrHz) && fdRefAbsErrHz <= config.fdRefHealthyAbsHz && ...
  isfinite(fdRateAbsErrHzPerSec) && fdRateAbsErrHzPerSec <= config.fdRateHealthyAbsHzPerSec;
if ~fdHealthy
  tailClass = "same-tooth + fd not healthy";
  return;
end
if isfinite(nonRefCoherenceFloor) && nonRefCoherenceFloor < config.coherenceCollapseThreshold
  tailClass = "same-tooth + fd healthy + non-ref coherence collapsed";
  return;
end
if isfinite(multiGapToTruthDoaDeg) && multiGapToTruthDoaDeg >= config.truthDoaGapThresholdDeg
  tailClass = "same-tooth + fd healthy + DoA/local-state basin";
  return;
end
tailClass = "same-tooth light/unclear";
end


function candidateProbeTable = localBuildCandidateProbeForRepeat(iRepeat, taskSeed, fixture, staticBundle, truth, pilotWave, carrierFreq, sampleRate, flowOpt, methodList, methodCaseCell, fdRangeOracle, fdRateRangeOracle, toothStepHz, config)
%LOCALBUILDCANDIDATEPROBEFORREPEAT Re-evaluate local candidate points without changing solver adoption.

candidateProbeTable = table();
targetLabel = "ms-mf-cp-u-in-tooth";
targetIdx = find([methodList.label] == targetLabel, 1, 'first');
if isempty(targetIdx) || targetIdx > numel(methodCaseCell)
  return;
end
targetMethod = methodList(targetIdx);
targetCase = methodCaseCell{targetIdx};
defaultOptVar = localResolveOptVarFromCase(targetCase);
if numel(defaultOptVar) < 4
  return;
end
[viewUse, ~, staticCaseUse] = localResolveMethodView(targetMethod, fixture, staticBundle);
seedDoaParam = localResolveDoaSeed(targetMethod, staticCaseUse, truth);
dynOpt = localBuildDynOptForMethod(targetMethod, fixture, flowOpt, seedDoaParam);
dynOpt.verbose = false;
[model, ~, ~] = buildDoaDopplerMfModel(viewUse.sceneSeq, viewUse.rxSigMf, pilotWave, ...
  carrierFreq, sampleRate, viewUse.doaGrid, fdRangeOracle, fdRateRangeOracle, dynOpt);
[~, ~, model] = buildDoaDopplerMfInit(model, defaultOptVar);

rowCell = cell(0, 1);
rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
  "default-final", "default-final", targetLabel, defaultOptVar, [0; 0], 0, 0, ...
  model, truth, toothStepHz);

baselineLabelList = [
  "ms-mf-cp-u-static-doa-fixed-in-tooth";
  "ms-mf-cp-u-truth-doa-fixed-in-tooth";
  "ms-mf-cp-u-wide-doa-in-tooth";
  "ms-mf-cp-u-truth-doa-oracle"];
baselineFamilyList = [
  "static-doa-fixed-final";
  "truth-doa-fixed-final";
  "wide-doa-final";
  "truth-doa-oracle-final"];
for iBase = 1:numel(baselineLabelList)
  baseIdx = find([methodList.label] == baselineLabelList(iBase), 1, 'first');
  if isempty(baseIdx) || baseIdx > numel(methodCaseCell)
    continue;
  end
  optVar = localResolveOptVarFromCase(methodCaseCell{baseIdx});
  if numel(optVar) ~= numel(defaultOptVar)
    continue;
  end
  rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
    baselineFamilyList(iBase), baselineFamilyList(iBase), baselineLabelList(iBase), optVar, ...
    optVar(1:2) - defaultOptVar(1:2), optVar(3) - defaultOptVar(3), optVar(4) - defaultOptVar(4), ...
    model, truth, toothStepHz);
end

probeDoaStepDegList = reshape(localGetFieldOrDefault(config, 'probeDoaStepDegList', 0), [], 1);
probeFdRefStepHzList = reshape(localGetFieldOrDefault(config, 'probeFdRefStepHzList', 0), [], 1);
for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    doaStep = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
    optVar = defaultOptVar;
    optVar(1:2) = optVar(1:2) + doaStep;
    tag = sprintf('pure-doa-dlat%+.4g-dlon%+.4g', doaStep(1), doaStep(2));
    rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
      "pure-doa-grid", string(tag), targetLabel, optVar, doaStep, 0, 0, ...
      model, truth, toothStepHz);
  end
end

for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    for iFd = 1:numel(probeFdRefStepHzList)
      doaStep = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
      fdRefStep = probeFdRefStepHzList(iFd);
      optVar = defaultOptVar;
      optVar(1:2) = optVar(1:2) + doaStep;
      optVar(3) = optVar(3) + fdRefStep;
      tag = sprintf('joint-dlat%+.4g-dlon%+.4g-dfd%+.4g', doaStep(1), doaStep(2), fdRefStep);
      rowCell{end + 1, 1} = localEvaluateCandidateProbePoint(iRepeat, taskSeed, ...
        "joint-doa-fdref-grid", string(tag), targetLabel, optVar, doaStep, fdRefStep, 0, ...
        model, truth, toothStepHz);
    end
  end
end

wideOptVar = [];
wideIdx = find([methodList.label] == "ms-mf-cp-u-wide-doa-in-tooth", 1, 'first');
if ~isempty(wideIdx) && wideIdx <= numel(methodCaseCell)
  wideOptVar = localResolveOptVarFromCase(methodCaseCell{wideIdx});
  if numel(wideOptVar) == numel(defaultOptVar)
    rowCell = localAppendCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
      "wide-centered-pure-doa-grid", "wide-centered-joint-doa-fdref-grid", ...
      "ms-mf-cp-u-wide-doa-in-tooth", wideOptVar, defaultOptVar, ...
      probeDoaStepDegList, probeFdRefStepHzList, model, truth, toothStepHz);
  end
end

coarseDoaStepDegList = reshape(localGetFieldOrDefault(config, 'probeCoarseDoaStepDegList', [-0.006; -0.004; -0.002; 0; 0.002; 0.004; 0.006]), [], 1);
staticDoaParam = localResolveCaseDoa(staticBundle.caseStaticMs);
rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
  "default-coarse-doa-grid", "ms-mf-cp-u-in-tooth", defaultOptVar, defaultOptVar, ...
  coarseDoaStepDegList, model, truth, toothStepHz);
staticCenterOptVar = defaultOptVar;
if numel(staticDoaParam) >= 2 && all(isfinite(staticDoaParam(1:2)))
  staticCenterOptVar(1:2) = staticDoaParam(1:2);
  rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
    "static-coarse-doa-grid", "ms-sf-static-doa-center", staticCenterOptVar, defaultOptVar, ...
    coarseDoaStepDegList, model, truth, toothStepHz);
end
singleIdx = find([methodList.label] == "ss-mf-cp-u-in-tooth", 1, 'first');
if ~isempty(singleIdx) && singleIdx <= numel(methodCaseCell)
  singleOptVar = localResolveOptVarFromCase(methodCaseCell{singleIdx});
  if numel(singleOptVar) == numel(defaultOptVar)
    rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
      "single-mf-coarse-doa-grid", "ss-mf-cp-u-in-tooth", singleOptVar, defaultOptVar, ...
      coarseDoaStepDegList, model, truth, toothStepHz);
  end
end
if numel(wideOptVar) == numel(defaultOptVar)
  rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, ...
    "wide-coarse-doa-grid", "ms-mf-cp-u-wide-doa-in-tooth", wideOptVar, defaultOptVar, ...
    coarseDoaStepDegList, model, truth, toothStepHz);
end

lineAlphaList = reshape(localGetFieldOrDefault(config, 'probeLineAlphaList', [0; 0.25; 0.5; 0.75; 0.85; 0.90; 0.95; 0.98; 0.99; 1]), [], 1);
truthDoaParam = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', defaultOptVar(1:2)), [], 1);
truthDoaParam = truthDoaParam(1:2);
staticDoaParam = localResolveCaseDoa(staticBundle.caseStaticMs);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
truthOracleIdx = find([methodList.label] == "ms-mf-cp-u-truth-doa-oracle", 1, 'first');
truthLineFd = [truthFdRefHz; truthFdRateHzPerSec];
if ~isempty(truthOracleIdx) && truthOracleIdx <= numel(methodCaseCell)
  truthOracleOptVar = localResolveOptVarFromCase(methodCaseCell{truthOracleIdx});
  if numel(truthOracleOptVar) == numel(defaultOptVar)
    truthLineFd = truthOracleOptVar(3:4);
  end
end
rowCell = localAppendDoaLineProbe(rowCell, iRepeat, taskSeed, ...
  "default-to-truth-line-default-fd", "default", "default-fd", ...
  defaultOptVar(1:2), truthDoaParam, defaultOptVar(3:4), defaultOptVar, lineAlphaList, ...
  model, truth, toothStepHz);
rowCell = localAppendDoaLineProbe(rowCell, iRepeat, taskSeed, ...
  "default-to-truth-line-truth-fd", "default", "truth-fd", ...
  defaultOptVar(1:2), truthDoaParam, truthLineFd, defaultOptVar, lineAlphaList, ...
  model, truth, toothStepHz);
rowCell = localAppendDoaLineProbe(rowCell, iRepeat, taskSeed, ...
  "static-to-truth-line-default-fd", "static", "default-fd", ...
  staticDoaParam, truthDoaParam, defaultOptVar(3:4), defaultOptVar, lineAlphaList, ...
  model, truth, toothStepHz);
rowCell = localAppendDoaLineProbe(rowCell, iRepeat, taskSeed, ...
  "static-to-truth-line-truth-fd", "static", "truth-fd", ...
  staticDoaParam, truthDoaParam, truthLineFd, defaultOptVar, lineAlphaList, ...
  model, truth, toothStepHz);

candidateProbeTable = struct2table([rowCell{:}].');
defaultMask = candidateProbeTable.candidateFamily == "default-final";
if any(defaultMask)
  defaultObj = candidateProbeTable.obj(find(defaultMask, 1, 'first'));
  candidateProbeTable.objGainFromDefault = defaultObj - candidateProbeTable.obj;
end
end

function rowCell = localAppendCenteredProbeGrid(rowCell, iRepeat, taskSeed, pureFamily, jointFamily, sourceMethodLabel, centerOptVar, defaultOptVar, probeDoaStepDegList, probeFdRefStepHzList, model, truth, toothStepHz)
%LOCALAPPENDCENTEREDPROBEGRID Add pure-DoA and DoA-fdRef probes around a non-default center.

centerOptVar = reshape(double(centerOptVar), [], 1);
defaultOptVar = reshape(double(defaultOptVar), [], 1);
if numel(centerOptVar) ~= numel(defaultOptVar) || numel(centerOptVar) < 4
  return;
end
centerLabel = string(sourceMethodLabel) + "-center";
for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    doaStepLocal = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
    optVar = centerOptVar;
    optVar(1:2) = optVar(1:2) + doaStepLocal;
    doaStepFromDefault = optVar(1:2) - defaultOptVar(1:2);
    tag = sprintf('wide-center-dlat%+.4g-dlon%+.4g', doaStepLocal(1), doaStepLocal(2));
    row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, pureFamily, string(tag), ...
      sourceMethodLabel, optVar, doaStepFromDefault, optVar(3) - defaultOptVar(3), ...
      optVar(4) - defaultOptVar(4), model, truth, toothStepHz);
    row.probeCenterFamily = centerLabel;
    rowCell{end + 1, 1} = row; %#ok<AGROW>
  end
end
for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    for iFd = 1:numel(probeFdRefStepHzList)
      doaStepLocal = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
      fdRefStepLocal = probeFdRefStepHzList(iFd);
      optVar = centerOptVar;
      optVar(1:2) = optVar(1:2) + doaStepLocal;
      optVar(3) = optVar(3) + fdRefStepLocal;
      doaStepFromDefault = optVar(1:2) - defaultOptVar(1:2);
      tag = sprintf('wide-center-dlat%+.4g-dlon%+.4g-dfd%+.4g', ...
        doaStepLocal(1), doaStepLocal(2), fdRefStepLocal);
      row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, jointFamily, string(tag), ...
        sourceMethodLabel, optVar, doaStepFromDefault, optVar(3) - defaultOptVar(3), ...
        optVar(4) - defaultOptVar(4), model, truth, toothStepHz);
      row.probeCenterFamily = centerLabel;
      rowCell{end + 1, 1} = row; %#ok<AGROW>
    end
  end
end
end

function rowCell = localAppendPureCenteredProbeGrid(rowCell, iRepeat, taskSeed, candidateFamily, sourceMethodLabel, centerOptVar, defaultOptVar, probeDoaStepDegList, model, truth, toothStepHz)
%LOCALAPPENDPURECENTEREDPROBEGRID Add a wider DoA-only grid around one implementable center.

centerOptVar = reshape(double(centerOptVar), [], 1);
defaultOptVar = reshape(double(defaultOptVar), [], 1);
probeDoaStepDegList = reshape(double(probeDoaStepDegList), [], 1);
if numel(centerOptVar) ~= numel(defaultOptVar) || numel(centerOptVar) < 4
  return;
end
if any(~isfinite(centerOptVar(1:4))) || any(~isfinite(defaultOptVar(1:4)))
  return;
end
centerLabel = string(sourceMethodLabel) + "-coarse-center";
for iLat = 1:numel(probeDoaStepDegList)
  for iLon = 1:numel(probeDoaStepDegList)
    doaStepLocal = [probeDoaStepDegList(iLat); probeDoaStepDegList(iLon)];
    optVar = centerOptVar;
    optVar(1:2) = optVar(1:2) + doaStepLocal;
    doaStepFromDefault = optVar(1:2) - defaultOptVar(1:2);
    tag = sprintf('coarse-center-dlat%+.4g-dlon%+.4g', doaStepLocal(1), doaStepLocal(2));
    row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, candidateFamily, string(tag), ...
      sourceMethodLabel, optVar, doaStepFromDefault, optVar(3) - defaultOptVar(3), ...
      optVar(4) - defaultOptVar(4), model, truth, toothStepHz);
    row.probeCenterFamily = centerLabel;
    rowCell{end + 1, 1} = row; %#ok<AGROW>
  end
end
end


function rowCell = localAppendDoaLineProbe(rowCell, iRepeat, taskSeed, candidateFamily, lineStartFamily, fdLineSource, startDoaParam, truthDoaParam, fdLine, defaultOptVar, lineAlphaList, model, truth, toothStepHz)
%LOCALAPPENDDOALINEPROBE Add objective probes on a DoA line segment toward truth.

startDoaParam = reshape(double(startDoaParam), [], 1);
truthDoaParam = reshape(double(truthDoaParam), [], 1);
fdLine = reshape(double(fdLine), [], 1);
defaultOptVar = reshape(double(defaultOptVar), [], 1);
lineAlphaList = reshape(double(lineAlphaList), [], 1);
if numel(startDoaParam) < 2 || numel(truthDoaParam) < 2 || numel(fdLine) < 2 || numel(defaultOptVar) < 4
  return;
end
if any(~isfinite(startDoaParam(1:2))) || any(~isfinite(truthDoaParam(1:2))) || any(~isfinite(fdLine(1:2)))
  return;
end
for iAlpha = 1:numel(lineAlphaList)
  alpha = lineAlphaList(iAlpha);
  if ~isfinite(alpha)
    continue;
  end
  optVar = defaultOptVar;
  optVar(1:2) = (1 - alpha) * startDoaParam(1:2) + alpha * truthDoaParam(1:2);
  optVar(3:4) = fdLine(1:2);
  doaStepFromDefault = optVar(1:2) - defaultOptVar(1:2);
  tag = sprintf('%s-to-truth-alpha%.4g-%s', char(lineStartFamily), alpha, char(fdLineSource));
  row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, candidateFamily, string(tag), ...
    "line-probe", optVar, doaStepFromDefault, optVar(3) - defaultOptVar(3), ...
    optVar(4) - defaultOptVar(4), model, truth, toothStepHz);
  row.probeCenterFamily = string(lineStartFamily) + "-to-truth-line";
  row.lineAlpha = alpha;
  row.lineStartFamily = string(lineStartFamily);
  row.fdLineSource = string(fdLineSource);
  rowCell{end + 1, 1} = row; %#ok<AGROW>
end
end

function candidateProbeTable = localBuildCandidateProbeTable(repeatCell)
%LOCALBUILDCANDIDATEPROBETABLE Stack per-repeat candidate objective probes.

probeCell = cell(0, 1);
for iRepeat = 1:numel(repeatCell)
  if isfield(repeatCell{iRepeat}, 'candidateProbeTable') && ~isempty(repeatCell{iRepeat}.candidateProbeTable)
    probeCell{end + 1, 1} = repeatCell{iRepeat}.candidateProbeTable; %#ok<AGROW>
  end
end
if isempty(probeCell)
  candidateProbeTable = table();
else
  candidateProbeTable = vertcat(probeCell{:});
end
end

function row = localEvaluateCandidateProbePoint(iRepeat, taskSeed, candidateFamily, candidateTag, sourceMethodLabel, optVar, doaStepDeg, fdRefStepHz, fdRateStepHzPerSec, model, truth, toothStepHz)
%LOCALEVALUATECANDIDATEPROBEPOINT Evaluate one candidate optVar and pack compact metrics.

optVar = reshape(double(optVar), [], 1);
[obj, ~, noiseVar, evalAux] = evaluateDoaDopplerMfObjective(model, optVar);
evalDiag = buildDoaDopplerMfEvalDiag(model, optVar, obj, noiseVar, evalAux);
truthLatlon = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', [NaN; NaN]), [], 1);
truthFdRefHz = localResolveTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localResolveTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
[nonRefCoherenceFloor, nonRefRmsPhaseResidRad] = localExtractProbeNonRefMetrics(evalDiag);
fdRefErrHz = optVar(3) - truthFdRefHz;
fdRateErrHzPerSec = optVar(4) - truthFdRateHzPerSec;
toothIdx = NaN;
toothResidualHz = NaN;
if isfinite(toothStepHz) && toothStepHz > 0 && isfinite(fdRefErrHz)
  toothIdx = round(fdRefErrHz / toothStepHz);
  toothResidualHz = fdRefErrHz - toothIdx * toothStepHz;
end
row = struct();
row.iRepeat = iRepeat;
row.taskSeed = taskSeed;
row.candidateFamily = string(candidateFamily);
row.candidateTag = string(candidateTag);
row.sourceMethodLabel = string(sourceMethodLabel);
row.doaStep1Deg = localVectorElem(doaStepDeg, 1, NaN);
row.doaStep2Deg = localVectorElem(doaStepDeg, 2, NaN);
row.fdRefStepHz = fdRefStepHz;
row.fdRateStepHzPerSec = fdRateStepHzPerSec;
row.probeCenterFamily = "";
row.lineAlpha = NaN;
row.lineStartFamily = "";
row.fdLineSource = "";
row.doa1Deg = optVar(1);
row.doa2Deg = optVar(2);
row.fdRefHz = optVar(3);
row.fdRateHzPerSec = optVar(4);
row.obj = obj;
row.objGainFromDefault = NaN;
row.residualNorm = localGetFieldOrDefault(evalDiag, 'residualNorm', NaN);
row.angleErrDeg = calcLatlonAngleError(optVar(1:2), truthLatlon(1:2));
row.fdRefErrHz = fdRefErrHz;
row.fdRateErrHzPerSec = fdRateErrHzPerSec;
row.toothIdx = toothIdx;
row.toothResidualHz = toothResidualHz;
row.nonRefCoherenceFloor = nonRefCoherenceFloor;
row.nonRefRmsPhaseResidRad = nonRefRmsPhaseResidRad;
row.nonRefFitRatioFloor = localGetFieldOrDefault(evalDiag, 'nonRefFitRatioFloor', NaN);
row.nonRefSupportRatioFloor = localGetFieldOrDefault(evalDiag, 'nonRefSupportRatioFloor', NaN);
end

function [nonRefCoherenceFloor, nonRefRmsPhaseResidRad] = localExtractProbeNonRefMetrics(evalDiag)
%LOCALEXTRACTPROBENONREFMETRICS Extract non-reference coherence and phase residual metrics.

nonRefCoherenceFloor = NaN;
nonRefRmsPhaseResidRad = NaN;
refSatIdxLocal = localGetFieldOrDefault(evalDiag, 'refSatIdxLocal', NaN);
coherenceSat = reshape(localGetFieldOrDefault(evalDiag, 'coherenceSat', []), [], 1);
if isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= numel(coherenceSat)
  nonRefMask = true(size(coherenceSat));
  nonRefMask(refSatIdxLocal) = false;
  nonRefVals = coherenceSat(nonRefMask & isfinite(coherenceSat));
  if ~isempty(nonRefVals)
    nonRefCoherenceFloor = min(nonRefVals);
  end
end
blockPhaseResidMat = localGetFieldOrDefault(evalDiag, 'blockPhaseResidMat', []);
if isempty(blockPhaseResidMat) || ~(isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= size(blockPhaseResidMat, 1))
  return;
end
nonRefMask = true(size(blockPhaseResidMat, 1), 1);
nonRefMask(refSatIdxLocal) = false;
nonRefPhaseResid = abs(blockPhaseResidMat(nonRefMask, :));
nonRefPhaseResid = nonRefPhaseResid(isfinite(nonRefPhaseResid));
if ~isempty(nonRefPhaseResid)
  nonRefRmsPhaseResidRad = sqrt(mean(nonRefPhaseResid .^ 2));
end
end

function optVar = localResolveOptVarFromCase(caseUse)
%LOCALRESOLVEOPTVARFROMCASE Read a final optimizer vector from one case result.

optVar = [];
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if ~isstruct(estResult)
  return;
end
optVar = reshape(localGetFieldOrDefault(estResult, 'optVarEst', []), [], 1);
if ~isempty(optVar)
  return;
end
doaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
fdRef = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRate = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
if numel(doaParam) >= 2 && isfinite(fdRef) && isfinite(fdRate)
  optVar = [doaParam(1:2); fdRef; fdRate];
end
end

function candidateWinnerTable = localBuildCandidateWinnerTable(candidateProbeTable, tailDiagnosisTable, config)
%LOCALBUILDCANDIDATEWINNERTABLE Summarize the best pure-DoA and joint local probes per seed.

if isempty(candidateProbeTable)
  candidateWinnerTable = table();
  return;
end
seedList = unique(candidateProbeTable.taskSeed, 'stable');
numSeed = numel(seedList);
rowCell = cell(numSeed, 1);
for iSeed = 1:numSeed
  seedNow = seedList(iSeed);
  seedProbe = candidateProbeTable(candidateProbeTable.taskSeed == seedNow, :);
  defaultRow = localBestProbeRow(seedProbe, "default-final");
  pureRow = localBestProbeRow(seedProbe, "pure-doa-grid");
  jointRow = localBestProbeRow(seedProbe, "joint-doa-fdref-grid");
  wideCenteredRow = localBestProbeRowAny(seedProbe, ["wide-centered-pure-doa-grid"; "wide-centered-joint-doa-fdref-grid"]);
  implementableRow = localBestProbeRowAny(seedProbe, ["default-coarse-doa-grid"; "static-coarse-doa-grid"; ...
    "single-mf-coarse-doa-grid"; "wide-coarse-doa-grid"]);
  defaultLineTruthFdRow = localBestProbeRow(seedProbe, "default-to-truth-line-truth-fd");
  staticLineTruthFdRow = localBestProbeRow(seedProbe, "static-to-truth-line-truth-fd");
  truthRow = localBestProbeRow(seedProbe, "truth-doa-oracle-final");
  tailRow = tailDiagnosisTable(tailDiagnosisTable.taskSeed == seedNow, :);

  row = struct();
  row.taskSeed = seedNow;
  row.tailClass = string(localTableValue(tailRow, 'tailClass', ""));
  row.defaultAngleErrDeg = localProbeValue(defaultRow, 'angleErrDeg', NaN);
  row.defaultObj = localProbeValue(defaultRow, 'obj', NaN);
  row.defaultNonRefCoherenceFloor = localProbeValue(defaultRow, 'nonRefCoherenceFloor', NaN);
  row.bestPureDoaObjGain = localProbeValue(pureRow, 'objGainFromDefault', NaN);
  row.bestPureDoaAngleErrDeg = localProbeValue(pureRow, 'angleErrDeg', NaN);
  row.bestPureDoaNonRefCoherenceFloor = localProbeValue(pureRow, 'nonRefCoherenceFloor', NaN);
  row.bestPureDoaStep1Deg = localProbeValue(pureRow, 'doaStep1Deg', NaN);
  row.bestPureDoaStep2Deg = localProbeValue(pureRow, 'doaStep2Deg', NaN);
  row.bestJointObjGain = localProbeValue(jointRow, 'objGainFromDefault', NaN);
  row.bestJointAngleErrDeg = localProbeValue(jointRow, 'angleErrDeg', NaN);
  row.bestJointNonRefCoherenceFloor = localProbeValue(jointRow, 'nonRefCoherenceFloor', NaN);
  row.bestJointDoaStep1Deg = localProbeValue(jointRow, 'doaStep1Deg', NaN);
  row.bestJointDoaStep2Deg = localProbeValue(jointRow, 'doaStep2Deg', NaN);
  row.bestJointFdRefStepHz = localProbeValue(jointRow, 'fdRefStepHz', NaN);
  row.bestWideCenteredFamily = string(localProbeValue(wideCenteredRow, 'candidateFamily', ""));
  row.bestWideCenteredObjGain = localProbeValue(wideCenteredRow, 'objGainFromDefault', NaN);
  row.bestWideCenteredAngleErrDeg = localProbeValue(wideCenteredRow, 'angleErrDeg', NaN);
  row.bestWideCenteredNonRefCoherenceFloor = localProbeValue(wideCenteredRow, 'nonRefCoherenceFloor', NaN);
  row.bestWideCenteredDoaStep1Deg = localProbeValue(wideCenteredRow, 'doaStep1Deg', NaN);
  row.bestWideCenteredDoaStep2Deg = localProbeValue(wideCenteredRow, 'doaStep2Deg', NaN);
  row.bestWideCenteredFdRefStepHz = localProbeValue(wideCenteredRow, 'fdRefStepHz', NaN);
  row.bestImplementableCenterFamily = string(localProbeValue(implementableRow, 'candidateFamily', ""));
  row.bestImplementableCenterObjGain = localProbeValue(implementableRow, 'objGainFromDefault', NaN);
  row.bestImplementableCenterAngleErrDeg = localProbeValue(implementableRow, 'angleErrDeg', NaN);
  row.bestImplementableCenterNonRefCoherenceFloor = localProbeValue(implementableRow, 'nonRefCoherenceFloor', NaN);
  row.bestImplementableCenterDoaStep1Deg = localProbeValue(implementableRow, 'doaStep1Deg', NaN);
  row.bestImplementableCenterDoaStep2Deg = localProbeValue(implementableRow, 'doaStep2Deg', NaN);
  row.bestDefaultTruthLineTruthFdAlpha = localProbeValue(defaultLineTruthFdRow, 'lineAlpha', NaN);
  row.bestDefaultTruthLineTruthFdAngleErrDeg = localProbeValue(defaultLineTruthFdRow, 'angleErrDeg', NaN);
  row.bestDefaultTruthLineTruthFdCoherenceFloor = localProbeValue(defaultLineTruthFdRow, 'nonRefCoherenceFloor', NaN);
  row.bestStaticTruthLineTruthFdAlpha = localProbeValue(staticLineTruthFdRow, 'lineAlpha', NaN);
  row.bestStaticTruthLineTruthFdAngleErrDeg = localProbeValue(staticLineTruthFdRow, 'angleErrDeg', NaN);
  row.bestStaticTruthLineTruthFdCoherenceFloor = localProbeValue(staticLineTruthFdRow, 'nonRefCoherenceFloor', NaN);
  row.truthDoaOracleAngleErrDeg = localProbeValue(truthRow, 'angleErrDeg', NaN);
  row.truthDoaOracleObjGain = localProbeValue(truthRow, 'objGainFromDefault', NaN);
  row.candidateClass = localClassifyCandidateRescue(row, config);
  rowCell{iSeed} = row;
end
candidateWinnerTable = struct2table([rowCell{:}].');
end

function bestRow = localBestProbeRow(seedProbe, candidateFamily)
%LOCALBESTPROBEROW Return the lowest-objective row for one candidate family.

bestRow = seedProbe([], :);
if isempty(seedProbe) || ~ismember('candidateFamily', seedProbe.Properties.VariableNames)
  return;
end
mask = seedProbe.candidateFamily == string(candidateFamily) & isfinite(seedProbe.obj);
if ~any(mask)
  return;
end
idx = find(mask);
[~, idxRel] = min(seedProbe.obj(idx));
bestRow = seedProbe(idx(idxRel), :);
end

function bestRow = localBestProbeRowAny(seedProbe, candidateFamilyList)
%LOCALBESTPROBEROWANY Return the lowest-objective row across several candidate families.

bestRow = seedProbe([], :);
if isempty(seedProbe) || ~ismember('candidateFamily', seedProbe.Properties.VariableNames)
  return;
end
candidateFamilyList = string(candidateFamilyList(:));
mask = ismember(seedProbe.candidateFamily, candidateFamilyList) & isfinite(seedProbe.obj);
if ~any(mask)
  return;
end
idx = find(mask);
[~, idxRel] = min(seedProbe.obj(idx));
bestRow = seedProbe(idx(idxRel), :);
end

function value = localProbeValue(row, fieldName, defaultValue)
%LOCALPROBEVALUE Read a scalar field from one probe table row.

value = defaultValue;
if isempty(row) || height(row) == 0 || ~ismember(fieldName, row.Properties.VariableNames)
  return;
end
rawValue = row.(fieldName);
if ~isempty(rawValue)
  value = rawValue(1);
end
end

function candidateClass = localClassifyCandidateRescue(row, config)
%LOCALCLASSIFYCANDIDATERESCUE Label whether local probes can rescue one tail.

angleGainThreshold = localGetFieldOrDefault(config, 'probeAngleGainThresholdDeg', 5e-4);
cohThreshold = localGetFieldOrDefault(config, 'coherenceCollapseThreshold', 0.95);
defaultAngle = row.defaultAngleErrDeg;
pureAngleGain = defaultAngle - row.bestPureDoaAngleErrDeg;
jointAngleGain = defaultAngle - row.bestJointAngleErrDeg;
wideAngleGain = defaultAngle - row.bestWideCenteredAngleErrDeg;
implementableAngleGain = defaultAngle - row.bestImplementableCenterAngleErrDeg;
lineAngleGain = defaultAngle - min([row.bestDefaultTruthLineTruthFdAngleErrDeg; row.bestStaticTruthLineTruthFdAngleErrDeg], [], 'omitnan');
truthAngleGain = defaultAngle - row.truthDoaOracleAngleErrDeg;
pureRescued = isfinite(pureAngleGain) && pureAngleGain >= angleGainThreshold && ...
  isfinite(row.bestPureDoaNonRefCoherenceFloor) && row.bestPureDoaNonRefCoherenceFloor >= cohThreshold;
jointRescued = isfinite(jointAngleGain) && jointAngleGain >= angleGainThreshold && ...
  isfinite(row.bestJointNonRefCoherenceFloor) && row.bestJointNonRefCoherenceFloor >= cohThreshold;
wideRescued = isfinite(wideAngleGain) && wideAngleGain >= angleGainThreshold && ...
  isfinite(row.bestWideCenteredNonRefCoherenceFloor) && row.bestWideCenteredNonRefCoherenceFloor >= cohThreshold;
implementableRescued = isfinite(implementableAngleGain) && implementableAngleGain >= angleGainThreshold && ...
  isfinite(row.bestImplementableCenterNonRefCoherenceFloor) && row.bestImplementableCenterNonRefCoherenceFloor >= cohThreshold;
lineRescued = isfinite(lineAngleGain) && lineAngleGain >= angleGainThreshold;
truthRescued = isfinite(truthAngleGain) && truthAngleGain >= angleGainThreshold;
if jointRescued && (~pureRescued || row.bestJointAngleErrDeg <= row.bestPureDoaAngleErrDeg - 1e-4)
  candidateClass = "joint-doa-fdref-rescuable";
elseif pureRescued
  candidateClass = "pure-doa-rescuable";
elseif wideRescued
  candidateClass = "wide-centered-rescuable";
elseif implementableRescued
  candidateClass = "implementable-center-rescuable";
elseif lineRescued
  candidateClass = "line-to-truth-rescuable";
elseif truthRescued
  candidateClass = "truth-only-rescuable";
else
  candidateClass = "not-rescued-by-local-probe";
end
end


function rescueBankSummaryTable = localBuildRescueBankSummaryTable(candidateProbeTable, tailDiagnosisTable, config)
%LOCALBUILDRESCUEBANKSUMMARYTABLE Compare implementable rescue banks against the disabled default.

rescueBankSummaryTable = table();
if isempty(candidateProbeTable)
  return;
end
seedList = unique(candidateProbeTable.taskSeed, 'stable');
bankLabelList = ["disabled"; "wide-centered-coarse"; "single-mf-centered-coarse"; ...
  "wide-single-bank"; "gated-wide-single-bank"];
numRow = numel(seedList) * numel(bankLabelList);
rowCell = cell(numRow, 1);
rowCount = 0;
angleGainThreshold = localGetFieldOrDefault(config, 'probeAngleGainThresholdDeg', 5e-4);
damageThreshold = localGetFieldOrDefault(config, 'probeDamageThresholdDeg', 5e-4);
cohThreshold = localGetFieldOrDefault(config, 'coherenceCollapseThreshold', 0.95);

for iSeed = 1:numel(seedList)
  seedNow = seedList(iSeed);
  seedProbe = candidateProbeTable(candidateProbeTable.taskSeed == seedNow, :);
  defaultRow = localBestProbeRow(seedProbe, "default-final");
  if isempty(defaultRow) || height(defaultRow) == 0
    continue;
  end
  tailRow = tailDiagnosisTable(tailDiagnosisTable.taskSeed == seedNow, :);
  tailClass = string(localTableValue(tailRow, 'tailClass', ""));
  truthRow = localBestProbeRow(seedProbe, "truth-doa-oracle-final");
  truthAngle = localProbeValue(truthRow, 'angleErrDeg', NaN);
  defaultAngle = localProbeValue(defaultRow, 'angleErrDeg', NaN);
  defaultObj = localProbeValue(defaultRow, 'obj', NaN);

  rowCount = rowCount + 1;
  rowCell{rowCount} = localMakeRescueBankRow(seedNow, tailClass, truthAngle, defaultAngle, defaultObj, ...
    "disabled", defaultRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold);

  wideRow = localBestProbeRow(seedProbe, "wide-coarse-doa-grid");
  rowCount = rowCount + 1;
  rowCell{rowCount} = localMakeRescueBankRow(seedNow, tailClass, truthAngle, defaultAngle, defaultObj, ...
    "wide-centered-coarse", wideRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold);

  singleRow = localBestProbeRow(seedProbe, "single-mf-coarse-doa-grid");
  rowCount = rowCount + 1;
  rowCell{rowCount} = localMakeRescueBankRow(seedNow, tailClass, truthAngle, defaultAngle, defaultObj, ...
    "single-mf-centered-coarse", singleRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold);

  bankRow = localBestProbeRowAny(seedProbe, ["default-final"; "wide-coarse-doa-grid"; "single-mf-coarse-doa-grid"]);
  rowCount = rowCount + 1;
  rowCell{rowCount} = localMakeRescueBankRow(seedNow, tailClass, truthAngle, defaultAngle, defaultObj, ...
    "wide-single-bank", bankRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold);

  defaultCoherenceFloor = localProbeValue(defaultRow, 'nonRefCoherenceFloor', NaN);
  defaultRmsPhaseResidRad = localProbeValue(defaultRow, 'nonRefRmsPhaseResidRad', NaN);
  [triggerRescue, gateReason] = localShouldTriggerGatedRescue(defaultCoherenceFloor, ...
    defaultRmsPhaseResidRad, config);
  if triggerRescue
    gatedRow = bankRow;
  else
    gatedRow = defaultRow;
  end
  rowCount = rowCount + 1;
  rowCell{rowCount} = localMakeRescueBankRow(seedNow, tailClass, truthAngle, defaultAngle, defaultObj, ...
    "gated-wide-single-bank", gatedRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold, ...
    triggerRescue, gateReason, defaultCoherenceFloor, defaultRmsPhaseResidRad);
end

rowCell = rowCell(1:rowCount);
if isempty(rowCell)
  return;
end
rescueBankSummaryTable = struct2table([rowCell{:}].');
end

function row = localMakeRescueBankRow(taskSeed, tailClass, truthAngleErrDeg, defaultAngleErrDeg, defaultObj, bankLabel, selectedRow, defaultAngle, angleGainThreshold, damageThreshold, cohThreshold, rescueTriggered, gateReason, gateCoherenceFloor, gateRmsPhaseResidRad)
%LOCALMAKERESCUEBANKROW Convert one selected rescue candidate into a summary row.

if nargin < 12
  rescueTriggered = false;
end
if nargin < 13
  gateReason = "not-a-gated-bank";
end
if nargin < 14
  gateCoherenceFloor = NaN;
end
if nargin < 15
  gateRmsPhaseResidRad = NaN;
end

row = struct();
row.taskSeed = taskSeed;
row.tailClass = string(tailClass);
row.bankLabel = string(bankLabel);
row.rescueTriggered = logical(rescueTriggered);
row.gateReason = string(gateReason);
row.gateNonRefCoherenceFloor = gateCoherenceFloor;
row.gateNonRefRmsPhaseResidRad = gateRmsPhaseResidRad;
row.defaultAngleErrDeg = defaultAngleErrDeg;
row.defaultObj = defaultObj;
row.truthDoaOracleAngleErrDeg = truthAngleErrDeg;
row.selectedCandidateFamily = string(localProbeValue(selectedRow, 'candidateFamily', ""));
row.selectedCandidateTag = string(localProbeValue(selectedRow, 'candidateTag', ""));
row.selectedSourceMethodLabel = string(localProbeValue(selectedRow, 'sourceMethodLabel', ""));
row.selectedObj = localProbeValue(selectedRow, 'obj', NaN);
row.selectedObjGainFromDefault = localProbeValue(selectedRow, 'objGainFromDefault', NaN);
row.selectedAngleErrDeg = localProbeValue(selectedRow, 'angleErrDeg', NaN);
row.selectedNonRefCoherenceFloor = localProbeValue(selectedRow, 'nonRefCoherenceFloor', NaN);
row.selectedFdRefErrHz = localProbeValue(selectedRow, 'fdRefErrHz', NaN);
row.selectedFdRateErrHzPerSec = localProbeValue(selectedRow, 'fdRateErrHzPerSec', NaN);
row.selectedToothIdx = localProbeValue(selectedRow, 'toothIdx', NaN);
row.selectedDoaStep1Deg = localProbeValue(selectedRow, 'doaStep1Deg', NaN);
row.selectedDoaStep2Deg = localProbeValue(selectedRow, 'doaStep2Deg', NaN);
row.angleGainFromDefaultDeg = defaultAngle - row.selectedAngleErrDeg;
row.gapToTruthDoaOracleDeg = row.selectedAngleErrDeg - truthAngleErrDeg;
row.isCoherenceRecovered = isfinite(row.selectedNonRefCoherenceFloor) && row.selectedNonRefCoherenceFloor >= cohThreshold;
row.isAngleRescued = isfinite(row.angleGainFromDefaultDeg) && row.angleGainFromDefaultDeg >= angleGainThreshold;
row.isHardRescued = row.isAngleRescued && row.isCoherenceRecovered;
row.isDamaged = isfinite(row.angleGainFromDefaultDeg) && row.angleGainFromDefaultDeg <= -damageThreshold;
row.caseRole = localClassifyRescueCaseRole(row.tailClass);
end

function [triggerRescue, gateReason] = localShouldTriggerGatedRescue(defaultCoherenceFloor, defaultRmsPhaseResidRad, config)
%LOCALSHOULDTRIGGERGATEDRESCUE Gate the implementable bank to collapse-like defaults.

cohThreshold = localGetFieldOrDefault(config, 'gatedRescueCoherenceThreshold', 0.20);
phaseThreshold = localGetFieldOrDefault(config, 'gatedRescuePhaseResidThresholdRad', 1.0);
cohCollapsed = isfinite(defaultCoherenceFloor) && isfinite(cohThreshold) && defaultCoherenceFloor < cohThreshold;
phaseBad = ~isfinite(phaseThreshold) || phaseThreshold <= 0 || ...
  (isfinite(defaultRmsPhaseResidRad) && defaultRmsPhaseResidRad >= phaseThreshold);
triggerRescue = cohCollapsed && phaseBad;
if triggerRescue
  gateReason = "non-ref-coherence-collapse";
elseif ~cohCollapsed
  gateReason = "coherence-not-collapsed";
else
  gateReason = "phase-residual-not-large";
end
end

function caseRole = localClassifyRescueCaseRole(tailClass)
%LOCALCLASSIFYRESCUECASEROLE Group tail classes for rescue-bank aggregates.

tailClass = string(tailClass);
if contains(tailClass, "non-ref coherence collapsed") || contains(tailClass, "DoA/local-state basin")
  caseRole = "hard-collapse";
elseif contains(tailClass, "fd not healthy")
  caseRole = "fd-not-healthy-negative";
elseif contains(tailClass, "light") || contains(tailClass, "unclear")
  caseRole = "easy-negative";
elseif contains(tailClass, "wrong-tooth")
  caseRole = "wrong-tooth-negative";
else
  caseRole = "other";
end
end

function rescueBankAggregateTable = localBuildRescueBankAggregateTable(rescueBankSummaryTable)
%LOCALBUILDRESCUEBANKAGGREGATETABLE Build rescue and damage rates for each implementable bank.

rescueBankAggregateTable = table();
if isempty(rescueBankSummaryTable)
  return;
end
bankList = unique(rescueBankSummaryTable.bankLabel, 'stable');
rowCell = cell(numel(bankList), 1);
for iBank = 1:numel(bankList)
  bankNow = bankList(iBank);
  mask = rescueBankSummaryTable.bankLabel == bankNow;
  hardMask = mask & rescueBankSummaryTable.caseRole == "hard-collapse";
  easyMask = mask & rescueBankSummaryTable.caseRole == "easy-negative";
  fdMask = mask & rescueBankSummaryTable.caseRole == "fd-not-healthy-negative";
  row = struct();
  row.bankLabel = bankNow;
  row.numSeed = nnz(mask);
  row.numHardSeed = nnz(hardMask);
  row.numEasyNegativeSeed = nnz(easyMask);
  row.numFdNegativeSeed = nnz(fdMask);
  row.triggerRate = localMeanLogical(rescueBankSummaryTable.rescueTriggered(mask));
  row.hardTriggerRate = localMeanLogical(rescueBankSummaryTable.rescueTriggered(hardMask));
  row.easyTriggerRate = localMeanLogical(rescueBankSummaryTable.rescueTriggered(easyMask));
  row.fdNegativeTriggerRate = localMeanLogical(rescueBankSummaryTable.rescueTriggered(fdMask));
  row.hardRescueRate = localMeanLogical(rescueBankSummaryTable.isHardRescued(hardMask));
  row.hardMedianAngleGainDeg = median(rescueBankSummaryTable.angleGainFromDefaultDeg(hardMask), 'omitnan');
  row.hardMedianSelectedAngleDeg = median(rescueBankSummaryTable.selectedAngleErrDeg(hardMask), 'omitnan');
  row.hardCoherenceRecoveredRate = localMeanLogical(rescueBankSummaryTable.isCoherenceRecovered(hardMask));
  row.easyDamageRate = localMeanLogical(rescueBankSummaryTable.isDamaged(easyMask));
  row.fdNegativeDamageRate = localMeanLogical(rescueBankSummaryTable.isDamaged(fdMask));
  row.overallMedianSelectedAngleDeg = median(rescueBankSummaryTable.selectedAngleErrDeg(mask), 'omitnan');
  row.overallMaxSelectedAngleDeg = max(rescueBankSummaryTable.selectedAngleErrDeg(mask), [], 'omitnan');
  row.overallMedianCoherenceFloor = median(rescueBankSummaryTable.selectedNonRefCoherenceFloor(mask), 'omitnan');
  rowCell{iBank} = row;
end
rescueBankAggregateTable = struct2table([rowCell{:}].');
end

function value = localMeanLogical(maskValue)
%LOCALMEANLOGICAL Return the mean of a logical vector with empty fallback.

if isempty(maskValue)
  value = NaN;
else
  value = mean(double(maskValue), 'omitnan');
end
end

function lineProbeSummaryTable = localBuildLineProbeSummaryTable(candidateProbeTable, config)
%LOCALBUILDLINEPROBESUMMARYTABLE Summarize DoA line probes toward truth for each tail seed.

lineProbeSummaryTable = table();
if isempty(candidateProbeTable) || ~ismember('lineAlpha', candidateProbeTable.Properties.VariableNames)
  return;
end
lineMaskAll = isfinite(candidateProbeTable.lineAlpha) & contains(candidateProbeTable.candidateFamily, "-to-truth-line-");
if ~any(lineMaskAll)
  return;
end
cohThreshold = localGetFieldOrDefault(config, 'coherenceCollapseThreshold', 0.95);
angleGainThreshold = localGetFieldOrDefault(config, 'probeAngleGainThresholdDeg', 5e-4);
seedList = unique(candidateProbeTable.taskSeed(lineMaskAll), 'stable');
lineFamilyList = unique(candidateProbeTable.candidateFamily(lineMaskAll), 'stable');
rowCell = cell(numel(seedList) * numel(lineFamilyList), 1);
rowCount = 0;
for iSeed = 1:numel(seedList)
  seedNow = seedList(iSeed);
  seedProbe = candidateProbeTable(candidateProbeTable.taskSeed == seedNow, :);
  defaultRow = localBestProbeRow(seedProbe, "default-final");
  defaultAngle = localProbeValue(defaultRow, 'angleErrDeg', NaN);
  defaultCoh = localProbeValue(defaultRow, 'nonRefCoherenceFloor', NaN);
  for iFamily = 1:numel(lineFamilyList)
    familyNow = lineFamilyList(iFamily);
    familyProbe = seedProbe(seedProbe.candidateFamily == familyNow & isfinite(seedProbe.lineAlpha), :);
    if isempty(familyProbe)
      continue;
    end
    bestRow = localBestProbeRow(seedProbe, familyNow);
    cohRow = localFirstLineProbeRow(familyProbe, @(rowNow) localProbeValue(rowNow, 'nonRefCoherenceFloor', NaN) >= cohThreshold);
    rescueRow = localFirstLineProbeRow(familyProbe, @(rowNow) ...
      defaultAngle - localProbeValue(rowNow, 'angleErrDeg', NaN) >= angleGainThreshold && ...
      localProbeValue(rowNow, 'nonRefCoherenceFloor', NaN) >= cohThreshold);
    alphaOneRow = localBestLineAlphaRow(familyProbe, 1);

    row = struct();
    row.taskSeed = seedNow;
    row.lineFamily = string(familyNow);
    row.lineStartFamily = string(localProbeValue(bestRow, 'lineStartFamily', ""));
    row.fdLineSource = string(localProbeValue(bestRow, 'fdLineSource', ""));
    row.defaultAngleErrDeg = defaultAngle;
    row.defaultNonRefCoherenceFloor = defaultCoh;
    row.bestAlpha = localProbeValue(bestRow, 'lineAlpha', NaN);
    row.bestObjGain = localProbeValue(bestRow, 'objGainFromDefault', NaN);
    row.bestAngleErrDeg = localProbeValue(bestRow, 'angleErrDeg', NaN);
    row.bestNonRefCoherenceFloor = localProbeValue(bestRow, 'nonRefCoherenceFloor', NaN);
    row.firstCoherenceAlpha = localProbeValue(cohRow, 'lineAlpha', NaN);
    row.firstCoherenceAngleErrDeg = localProbeValue(cohRow, 'angleErrDeg', NaN);
    row.firstCoherenceObjGain = localProbeValue(cohRow, 'objGainFromDefault', NaN);
    row.firstRescueAlpha = localProbeValue(rescueRow, 'lineAlpha', NaN);
    row.firstRescueAngleErrDeg = localProbeValue(rescueRow, 'angleErrDeg', NaN);
    row.firstRescueObjGain = localProbeValue(rescueRow, 'objGainFromDefault', NaN);
    row.alphaOneObjGain = localProbeValue(alphaOneRow, 'objGainFromDefault', NaN);
    row.alphaOneAngleErrDeg = localProbeValue(alphaOneRow, 'angleErrDeg', NaN);
    row.alphaOneNonRefCoherenceFloor = localProbeValue(alphaOneRow, 'nonRefCoherenceFloor', NaN);
    rowCount = rowCount + 1;
    rowCell{rowCount} = row;
  end
end
rowCell = rowCell(1:rowCount);
if ~isempty(rowCell)
  lineProbeSummaryTable = struct2table([rowCell{:}].');
end
end

function bestRow = localFirstLineProbeRow(familyProbe, predicateFcn)
%LOCALFIRSTLINEPROBEROW Return the earliest alpha row that satisfies a predicate.

bestRow = familyProbe([], :);
if isempty(familyProbe)
  return;
end
[~, order] = sort(familyProbe.lineAlpha);
for iOrder = 1:numel(order)
  rowNow = familyProbe(order(iOrder), :);
  if predicateFcn(rowNow)
    bestRow = rowNow;
    return;
  end
end
end

function bestRow = localBestLineAlphaRow(familyProbe, alphaTarget)
%LOCALBESTLINEALPHAROW Return the line-probe row closest to a target alpha.

bestRow = familyProbe([], :);
if isempty(familyProbe) || ~ismember('lineAlpha', familyProbe.Properties.VariableNames)
  return;
end
alphaErr = abs(familyProbe.lineAlpha - alphaTarget);
alphaErr(~isfinite(alphaErr)) = inf;
[bestErr, idx] = min(alphaErr);
if isfinite(bestErr)
  bestRow = familyProbe(idx, :);
end
end

function detailDiagnosticTable = localBuildDetailDiagnosticTable(methodTable, tailDiagnosisTable)
%LOCALBUILDDETAILDIAGNOSTICTABLE Build per-method diagnostics, falling back to all seeds when no tail repeats.

detailDiagnosticTable = table();
if isempty(methodTable) || isempty(tailDiagnosisTable)
  return;
end
if ~ismember('tailClass', tailDiagnosisTable.Properties.VariableNames)
  return;
end
seedMask = tailDiagnosisTable.tailClass ~= "same-tooth light/unclear";
detailSeedList = tailDiagnosisTable.taskSeed(seedMask);
if isempty(detailSeedList)
  detailSeedList = tailDiagnosisTable.taskSeed;
end
mask = ismember(methodTable.taskSeed, detailSeedList);
if ~any(mask)
  return;
end
baseFields = {'taskSeed', 'methodLabel', 'angleErrDeg', 'fdRefErrHz', 'fdRateErrHzPerSec', ...
  'toothIdx', 'toothResidualHz', 'initObj', 'finalObj', 'initResidualNorm', 'finalResidualNorm', ...
  'initCohSat1', 'initCohSat2', 'finalCohSat1', 'finalCohSat2', ...
  'initResidualSat1', 'initResidualSat2', 'finalResidualSat1', 'finalResidualSat2', ...
  'initDeltaFdRefErrSat1Hz', 'initDeltaFdRefErrSat2Hz', ...
  'finalDeltaFdRefErrSat1Hz', 'finalDeltaFdRefErrSat2Hz', ...
  'initDeltaFdRateErrSat1HzPerSec', 'initDeltaFdRateErrSat2HzPerSec', ...
  'finalDeltaFdRateErrSat1HzPerSec', 'finalDeltaFdRateErrSat2HzPerSec', ...
  'nonRefCoherenceFloor', 'nonRefRmsPhaseResidRad'};
keepFields = baseFields(ismember(baseFields, methodTable.Properties.VariableNames));
detailDiagnosticTable = methodTable(mask, keepFields);
methodOrderList = unique(methodTable.methodLabel, 'stable');
methodOrder = nan(height(detailDiagnosticTable), 1);
for iRow = 1:height(detailDiagnosticTable)
  methodOrder(iRow) = find(methodOrderList == detailDiagnosticTable.methodLabel(iRow), 1, 'first');
end
[~, order] = sortrows([detailDiagnosticTable.taskSeed, methodOrder]);
detailDiagnosticTable = detailDiagnosticTable(order, :);
end

function evalDiag = localGetEvalDiag(caseUse, evalName)
%LOCALGETEVALDIAG Read one estimator debug eval block with fallback.

evalDiag = struct();
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
aux = localGetFieldOrDefault(estResult, 'aux', struct());
debugAux = localGetFieldOrDefault(aux, 'debug', struct());
evalDiag = localGetFieldOrDefault(debugAux, evalName, struct());
end

function deltaFdRef = localGetEvalDeltaFdRef(evalDiag)
%LOCALGETEVALDELTAFDREF Read the eval-domain differential Doppler vector.

deltaFdRef = localGetFieldOrDefault(evalDiag, 'deltaFdRefEval', ...
  localGetFieldOrDefault(evalDiag, 'deltaFdRef', []));
end

function deltaFdRate = localGetEvalDeltaFdRate(evalDiag)
%LOCALGETEVALDELTAFDRATE Read the eval-domain differential Doppler-rate vector.

deltaFdRate = localGetFieldOrDefault(evalDiag, 'deltaFdRateEval', ...
  localGetFieldOrDefault(evalDiag, 'deltaFdRate', []));
end

function [truthDeltaFdRefHz, truthDeltaFdRateHzPerSec] = localBuildTruthDeltaFdFit(truth)
%LOCALBUILDTRUTHDELTAFDFIT Resolve truth differential Doppler line fits.

truthDeltaFdRefHz = reshape(localGetFieldOrDefault(truth, 'deltaFdTrueHz', []), [], 1);
truthDeltaFdRateHzPerSec = reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1);
deltaFdSeries = localGetFieldOrDefault(truth, 'deltaFdSeries', []);
timeOffsetSec = reshape(localGetFieldOrDefault(truth, 'timeOffsetSec', []), [], 1);
if ~isempty(deltaFdSeries) && ~isempty(timeOffsetSec) && size(deltaFdSeries, 2) == numel(timeOffsetSec)
  numSat = size(deltaFdSeries, 1);
  truthDeltaFdRefHz = nan(numSat, 1);
  truthDeltaFdRateHzPerSec = nan(numSat, 1);
  for iSat = 1:numSat
    [truthDeltaFdRefHz(iSat), truthDeltaFdRateHzPerSec(iSat)] = ...
      localFitFdLine(timeOffsetSec, deltaFdSeries(iSat, :));
  end
end
end

function [fit0, fitRate] = localFitFdLine(timeOffsetSec, fdSeries)
%LOCALFITFDLINE Fit fdSeries ~= fit0 + fitRate * timeOffsetSec.

timeOffsetSec = reshape(timeOffsetSec, [], 1);
fdSeries = reshape(fdSeries, [], 1);
validMask = isfinite(timeOffsetSec) & isfinite(fdSeries);
if nnz(validMask) < 2
  fit0 = NaN;
  fitRate = NaN;
  return;
end
A = [ones(nnz(validMask), 1), timeOffsetSec(validMask)];
coeff = A \ fdSeries(validMask);
fit0 = coeff(1);
fitRate = coeff(2);
end

function value = localVectorElem(vec, idx, defaultValue)
%LOCALVECTORELEM Read one vector element with fallback.

value = defaultValue;
if isempty(vec)
  return;
end
vec = reshape(vec, [], 1);
if idx >= 1 && idx <= numel(vec) && isfinite(vec(idx))
  value = vec(idx);
end
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

function contextSummary = localBuildContextSummary(context, methodList)
%LOCALBUILDCONTEXTSUMMARY Keep only lightweight context metadata.

contextSummary = struct();
contextSummary.frameIntvlSec = context.frameIntvlSec;
contextSummary.periodicOffsetIdx = reshape(context.periodicOffsetIdx, 1, []);
contextSummary.masterOffsetIdx = reshape(context.masterOffsetIdx, 1, []);
contextSummary.selectedSatIdxGlobal = reshape(context.selectedSatIdxGlobal, 1, []);
contextSummary.refSatIdxGlobal = context.refSatIdxGlobal;
contextSummary.otherSatIdxGlobal = context.otherSatIdxGlobal;
contextSummary.contextBaseSeed = context.baseSeed;
contextSummary.tleFileName = string(context.tleFileName);
contextSummary.methodLabelList = reshape([methodList.label], [], 1);
contextSummary.subsetInitialization = "disabled: static/truth DoA seeds only";
end

function methodListOut = localStripMethodList(methodList)
%LOCALSTRIPMETHODLIST Store method descriptors without heavy fields.

methodListOut = methodList;
end

function plotData = localBuildPlotData(methodTable, tailDiagnosisTable)
%LOCALBUILDPLOTDATA Pack lightweight plotting inputs.

plotData = struct();
plotData.methodTable = methodTable(:, {'taskSeed', 'methodLabel', 'angleErrDeg', 'fdRefErrHz', ...
  'fdRateErrHzPerSec', 'toothIdx', 'nonRefCoherenceFloor'});
plotData.tailDiagnosisTable = tailDiagnosisTable;
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
%LOCALPLOTREPLAY Draw compact before/after rescue diagnostics.

if isempty(plotData) || ~isstruct(plotData)
  return;
end
rescueBankTable = localGetFieldOrDefault(plotData, 'rescueBankSummaryTable', table());
tailTable = localGetFieldOrDefault(plotData, 'tailDiagnosisTable', table());
if localPlotGatedRescueBeforeAfter(rescueBankTable, tailTable)
  return;
end
localPlotTailMethodBarFallback(tailTable);
end

function didPlot = localPlotGatedRescueBeforeAfter(rescueBankTable, tailTable)
%LOCALPLOTGATEDRESCUEBEFOREAFTER Plot default-vs-gated rescue results without line links.

didPlot = false;
if isempty(rescueBankTable) || ~ismember('bankLabel', rescueBankTable.Properties.VariableNames)
  return;
end
gatedMask = string(rescueBankTable.bankLabel) == "gated-wide-single-bank";
if ~any(gatedMask)
  return;
end
gatedTable = rescueBankTable(gatedMask, :);
if isempty(gatedTable)
  return;
end
seedLabel = string(gatedTable.taskSeed);
x = 1:height(gatedTable);
defaultAngle = gatedTable.defaultAngleErrDeg;
selectedAngle = gatedTable.selectedAngleErrDeg;
defaultCoherence = gatedTable.gateNonRefCoherenceFloor;
if any(~isfinite(defaultCoherence)) && ~isempty(tailTable) && ismember('multiNonRefCoherenceFloor', tailTable.Properties.VariableNames)
  defaultCoherence = localLookupTailValue(tailTable, gatedTable.taskSeed, 'multiNonRefCoherenceFloor', defaultCoherence);
end
selectedCoherence = gatedTable.selectedNonRefCoherenceFloor;

figure('Name', 'In-tooth gated rescue before-after diagnostics');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
bar(x, max([defaultAngle, selectedAngle], eps));
set(gca, 'YScale', 'log');
xticks(x);
xticklabels(seedLabel);
ylabel('angle error (deg)');
title('Default MS-MF vs gated rescue selected angle');
grid on;
legend({'default MS-MF', 'gated selected'}, 'Location', 'best', 'Interpreter', 'none');

nexttile;
bar(x, gatedTable.angleGainFromDefaultDeg);
yline(0, '--');
xticks(x);
xticklabels(seedLabel);
ylabel('angle gain (deg)');
title('Angle gain from gated rescue');
grid on;

nexttile;
bar(x, [defaultCoherence, selectedCoherence]);
yline(0.20, '--', 'gate');
yline(0.95, ':', 'healthy');
xticks(x);
xticklabels(seedLabel);
xlabel('seed');
ylabel('non-ref coherence floor');
title('Default MS-MF vs gated selected coherence');
grid on;
legend({'default MS-MF', 'gated selected'}, 'Location', 'best', 'Interpreter', 'none');
didPlot = true;
end

function valueVec = localLookupTailValue(tailTable, taskSeedVec, varName, defaultVec)
%LOCALLOOKUPTAILVALUE Look up one tail table column by task seed.

valueVec = defaultVec;
if isempty(tailTable) || ~ismember(varName, tailTable.Properties.VariableNames)
  return;
end
for iSeed = 1:numel(taskSeedVec)
  rowIdx = find(tailTable.taskSeed == taskSeedVec(iSeed), 1, 'first');
  if ~isempty(rowIdx)
    valueVec(iSeed) = tailTable.(varName)(rowIdx);
  end
end
end

function localPlotTailMethodBarFallback(tailTable)
%LOCALPLOTTAILMETHODBARFALLBACK Plot method comparison bars for older snapshots.

if isempty(tailTable)
  return;
end
seedLabel = string(tailTable.taskSeed);
x = 1:height(tailTable);
figure('Name', 'In-tooth tail method diagnostics');
tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile;
bar(x, max([tailTable.singleAngleErrDeg, tailTable.multiAngleErrDeg, ...
  tailTable.wideDoaAngleErrDeg, tailTable.truthDoaAngleErrDeg], eps));
set(gca, 'YScale', 'log');
xticks(x);
xticklabels(seedLabel);
ylabel('angle error (deg)');
title('Tail angle comparison');
grid on;
legend({'SS-MF CP-U', 'MS-MF CP-U', 'MS-MF wide-DoA', 'MS-MF truth-DoA'}, ...
  'Location', 'best', 'Interpreter', 'none');
nexttile;
bar(x, tailTable.multiNonRefCoherenceFloor);
yline(0.95, '--');
xticks(x);
xticklabels(seedLabel);
ylabel('non-ref coherence floor');
title('Non-reference coherence health');
grid on;
nexttile;
bar(x, tailTable.multiGapToTruthDoaDeg);
yline(0, '--');
xticks(x);
xticklabels(seedLabel);
xlabel('seed');
ylabel('MS gap to truth-DoA (deg)');
title('Gap to truth-DoA oracle');
grid on;
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
  error('replayMfInToothTailCaseDiagnose:MissingToothStep', ...
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
