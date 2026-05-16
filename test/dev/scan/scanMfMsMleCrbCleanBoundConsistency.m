% scanMfMsMleCrbCleanBoundConsistency
% Purpose: scan clean truth-centered local search bounds for SS/MS,
% SF/MF, CP/IP, and K/U MLE-vs-CRB consistency from the formal estimator
% path only.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves scanData in the workspace; checkpointEnable=true resumes
% interrupted task runs from per-task files under the repo-root tmp.
% saveSnapshot=true saves only scanData via saveExpSnapshot. Telegram notice
% is best-effort only and never affects numerical results.

clear; close all; clc;

%% Scan configuration

scanName = "scanMfMsMleCrbCleanBoundConsistency";
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = false;
notifyTelegramEnable = true;
optVerbose = false;

% Dynamic MF estimator dialect. The default is the clean MLE core: no CP
% engineering penalties, no fd alias unwrap, and no estimator-internal
% basin-entry / warm-anchor route. Keep this scan on core-only unless
% temporarily running a manual route diagnostic.
dynamicObjectiveMode = "pure-mle";
dynamicRouteMode = "core-only";
% Paper-model override shared with replayMfMsMleCrbCleanTrim: fixed steering
% at the reference frame and affine reference-tied Doppler over absolute time.
cleanMfSignalTimeModel = "linearRef";
cleanMfSteeringMode = "frozenRef";

numFrame = 10;
contextUsrLla = [55; 36.59; 0];
contextUtcRef = datetime([2026, 03, 18, 17, 08, 00], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
contextTleFileName = "statlink_20260318.tle";
contextSelectedSatIdxGlobal = [5259, 1243, 348, 5652, 14, 4437];
contextRefSatIdxGlobal = 5259;
% This scan fixes the CleanTrim numerical path and uses one explicit
% method-level local search envelope table. Tune rows in methodRangeTable when
% SS/MS, CP/IP, K/U, DoA and fd ranges need different local bounds.
snrDbList = (-15:5:10).';
baseSeed = 253;
numRepeat = 500;
seedList = baseSeed + (0:(numRepeat - 1));
methodNameList = [
  "SS-SF-DoA"
  "MS-SF-DoA"
  "SS-MF-Static"
  "MS-MF-Static"
  "SS-MF-CP-K"
  "SS-MF-CP-U"
  "MS-MF-CP-K"
  "MS-MF-CP-U"
  "SS-MF-IP-K"
  "SS-MF-IP-U"
  "MS-MF-IP-K"
  "MS-MF-IP-U"
];

% Per-method truth-centered local search envelope. This table is the intended
% tuning surface for the scan: each method can use a different DoA CRB-scale,
% fdRef range mode, fdRef scale/tooth width, and fdRate half-width.
% fdRefRangeMode="crb-scale" matches scanMfMsMleCrbCleanBoundScale and uses
% fdRefCrbScale as the actual half-width in CRB std units, with a tooth cap.
defaultDoaCrbScale = 2.0;
defaultFdRefRangeMode = "crb-scale";
defaultFdHalfToothFraction = 0.3;
defaultFdRefCrbScale = 2.5;
defaultFdRefMaxHalfToothFraction = 0.45;
defaultFdRateHalfWidthHzPerSec = 1000;
methodRangeTable = localBuildDefaultMethodRangeTable(methodNameList, ...
  defaultDoaCrbScale, defaultFdRefRangeMode, defaultFdHalfToothFraction, ...
  defaultFdRefCrbScale, defaultFdRateHalfWidthHzPerSec);
% Example: tune one route without adding an extra scan dimension.
% methodRangeTable.doaCrbScale(methodRangeTable.methodName == "MS-MF-CP-U") = 2.5;
% methodRangeTable.fdRefRangeMode(methodRangeTable.methodName == "MS-MF-IP-U") = "tooth-floor";
% methodRangeTable.fdRefHalfToothFraction(methodRangeTable.methodName == "MS-MF-IP-U") = 0.2;
% methodRangeTable.doaCrbScale(methodRangeTable.methodName == "SS-SF-Static") = 2.55;
% methodRangeTable.fdRefHalfToothFraction(methodRangeTable.methodName == "SS-SF-Static") = 0.1;
% methodRangeTable.doaCrbScale(methodRangeTable.methodName == "MS-SF-Static") = 100;

boundaryTolFraction = 0.02;
healthFdRefToothFraction = 0.25;
healthFdRateAbsMaxHzPerSec = 250;
trimAngleNormMax = 5;
trimFdRefNormMax = 5;
crbTargetMseRatio = 1.1;

seedList = reshape(double(seedList), [], 1);
snrDbList = reshape(double(snrDbList), [], 1);
methodNameList = reshape(string(methodNameList), [], 1);
methodRangeTable = localNormalizeMethodRangeTable(methodRangeTable, methodNameList, defaultFdRefRangeMode);

runTic = tic;
scanData = struct();
config = struct();
checkpointRunDir = "";
runState = struct();

try
  %% Build context and scan options

  config.scanName = string(scanName);
  config.numFrame = numFrame;
  config.contextUsrLla = reshape(double(contextUsrLla), [], 1);
  config.contextUtcRef = contextUtcRef;
  config.contextTleFileName = string(contextTleFileName);
  config.contextSelectedSatIdxGlobal = reshape(double(contextSelectedSatIdxGlobal), 1, []);
  config.contextRefSatIdxGlobal = double(contextRefSatIdxGlobal);
  config.snrDbList = snrDbList;
  config.seedList = seedList;
  config.baseSeed = seedList(1);
  config.numRepeat = numel(seedList);
  config.methodNameList = methodNameList;
  config.methodRangeTable = methodRangeTable;
  config.defaultDoaCrbScale = defaultDoaCrbScale;
  config.defaultFdRefRangeMode = string(defaultFdRefRangeMode);
  config.defaultFdHalfToothFraction = defaultFdHalfToothFraction;
  config.defaultFdRefCrbScale = defaultFdRefCrbScale;
  config.defaultFdRefMaxHalfToothFraction = defaultFdRefMaxHalfToothFraction;
  config.defaultFdRateHalfWidthHzPerSec = defaultFdRateHalfWidthHzPerSec;
  config.boundaryTolFraction = boundaryTolFraction;
  config.healthFdRefToothFraction = healthFdRefToothFraction;
  config.healthFdRateAbsMaxHzPerSec = healthFdRateAbsMaxHzPerSec;
  config.trimAngleNormMax = trimAngleNormMax;
  config.trimFdRefNormMax = trimFdRefNormMax;
  config.crbTargetMseRatio = crbTargetMseRatio;
  config.saveSnapshot = logical(saveSnapshot);
  config.checkpointEnable = logical(checkpointEnable);
  config.checkpointResume = logical(checkpointResume);
  config.checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
  config.notifyTelegramEnable = logical(notifyTelegramEnable);
  config.optVerbose = logical(optVerbose);
  config.dynamicObjectiveMode = localNormalizeDynamicObjectiveMode(dynamicObjectiveMode);
  config.dynamicRouteMode = localNormalizeDynamicRouteMode(dynamicRouteMode);
  config.cleanMfSignalTimeModel = localNormalizeCleanMfSignalTimeModel(cleanMfSignalTimeModel);
  config.cleanMfSteeringMode = localNormalizeCleanMfSteeringMode(cleanMfSteeringMode);
  config.initMode = "truth-crb-scaled-static-seed-fixed-bound";
  config.sweepMode = "snr-method-range";
  config.routeTag = "clean-bound-ms-ss-cp-ip-l" + ...
    string(numel(config.contextSelectedSatIdxGlobal)) + "-ref" + ...
    string(config.contextRefSatIdxGlobal) + "-" + ...
    config.dynamicObjectiveMode + "-" + config.dynamicRouteMode + "-" + ...
    config.cleanMfSignalTimeModel + "-" + config.cleanMfSteeringMode;

  taskList = localBuildTaskList(config.snrDbList, config.seedList);
  config.numTask = numel(taskList);
  config.useParfor = localCanUseParfor(config.numTask);
  config.runKey = localBuildRunKey(config);
  checkpointOpt = localBuildCheckpointOpt(scanName, config, config.useParfor);
  if config.checkpointEnable
    checkpointRunDir = string(checkpointOpt.runDir);
    config.checkpointRunDir = checkpointRunDir;
  end

  printMfScanHeader(char(scanName), config, char(checkpointRunDir));
  localPrintReplayConfig(config);

  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(config.numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.usrLla = config.contextUsrLla;
  contextOpt.utcRef = config.contextUtcRef;
  contextOpt.tleFileName = config.contextTleFileName;
  contextOpt.selectedSatIdxGlobal = config.contextSelectedSatIdxGlobal;
  contextOpt.refSatIdxGlobal = config.contextRefSatIdxGlobal;
  contextOpt.satPickCount = numel(config.contextSelectedSatIdxGlobal);
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = localBuildSerialInnerParallelOpt();
  context = buildDynamicMultiSatEciContext(contextOpt);
  context = localApplyCleanSignalModelOverride(context, config);
  context.subsetOffsetCell = {};
  context.subsetLabelList = strings(0, 1);
  context.numSubsetRandomTrial = 0;
  flowOpt = localBuildFlowOpt(context.parallelOpt);

  %% Run clean scan batch

  try
    [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, scanName);
  catch ME
    if strlength(string(checkpointRunDir)) > 0
      fprintf('Scan failed. Checkpoint artifacts kept at: %s\n', char(checkpointRunDir));
    end
    rethrow(ME);
  end

  caseTable = localBuildCaseTable(repeatCell);
  caseTable = localAnnotateHealthFlags(caseTable, config);
  cleanAggregateTable = buildMfMetricAggregateTable(caseTable);
  cleanHealthAggregateTable = buildMfFilteredMetricAggregateTable(caseTable, "health", caseTable.healthResolved);
  cleanTrimAggregateTable = buildMfFilteredMetricAggregateTable(caseTable, "joint-trim", caseTable.jointTrimKeep);
  cleanSsMsCompareTable = localBuildCleanSsMsCompareTable(caseTable, config);
  cleanKnownUnknownCompareTable = localBuildCleanKnownUnknownCompareTable(caseTable);
  cleanSnrCurveTable = localBuildCleanSnrCurveTable(cleanSsMsCompareTable);
  cleanEstimateCrbCurveTable = localBuildCleanEstimateCrbCurveTable( ...
    cleanTrimAggregateTable, cleanHealthAggregateTable, cleanAggregateTable);
  cleanSearchRangeAuditTable = localBuildCleanSearchRangeAuditTable(caseTable);
  cleanOutlierTable = localBuildCleanOutlierTable(caseTable);
  runtimeTable = buildMfRuntimeTable(repeatCell);
  runtimeAggregateTable = buildMfRuntimeAggregateTable(runtimeTable);
  topSlowRuntimeTable = buildMfTopSlowRuntimeTable(runtimeTable, 12);

  %% Data storage

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = config.runKey;
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = config;
  scanData.contextSummary = localBuildContextSummary(context);
  scanData.caseTable = caseTable;
  scanData.cleanAggregateTable = cleanAggregateTable;
  scanData.cleanHealthAggregateTable = cleanHealthAggregateTable;
  scanData.cleanTrimAggregateTable = cleanTrimAggregateTable;
  scanData.cleanSsMsCompareTable = cleanSsMsCompareTable;
  scanData.cleanKnownUnknownCompareTable = cleanKnownUnknownCompareTable;
  scanData.cleanSnrCurveTable = cleanSnrCurveTable;
  scanData.cleanEstimateCrbCurveTable = cleanEstimateCrbCurveTable;
  scanData.cleanSearchRangeAuditTable = cleanSearchRangeAuditTable;
  scanData.cleanOutlierTable = cleanOutlierTable;
  scanData.runtimeTable = runtimeTable;
  scanData.runtimeAggregateTable = runtimeAggregateTable;
  scanData.topSlowRuntimeTable = topSlowRuntimeTable;
  scanData.repeatCell = localStripRepeatCell(repeatCell);
  scanData.checkpointDir = string(checkpointRunDir);
  scanData.checkpointCleanupReport = struct([]);
  scanData.checkpointSummary = buildMfReplayCheckpointSummary(runState);
  if config.checkpointEnable && localGetFieldOrDefault(runState, 'isComplete', false) && ...
      config.checkpointCleanupOnSuccess
    scanData.checkpointCleanupReport = cleanupPerfTaskGridCheckpoint(runState, struct('verbose', false));
    scanData.checkpointSummary.cleanupReport = scanData.checkpointCleanupReport;
    scanData.checkpointSummary.cleanedOnSuccess = true;
  elseif isfield(scanData.checkpointSummary, 'cleanedOnSuccess')
    scanData.checkpointSummary.cleanedOnSuccess = false;
  end
  scanData.checkpointSummaryTable = localBuildCheckpointSummaryTable(scanData.checkpointSummary, ...
    scanData.checkpointDir);
  scanData.plotData = localBuildCleanPlotData(cleanEstimateCrbCurveTable);
  scanData.elapsedSec = toc(runTic);
  scanData = finalizeMfScanResult(scanData, "");

  if config.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output

  if ~exist('scanData', 'var') || ~isstruct(scanData)
    error('scanMfMsMleCrbCleanBoundConsistency:MissingScanData', ...
      'Load a snapshot containing scanData before running the summary section.');
  end
  scanData = localValidateScanDataForSummary(scanData);
  scanConfigForReport = localGetFieldOrDefault(scanData, 'config', struct());
  scanNameForReport = string(localGetFieldOrDefault(scanData, 'scanName', "scanMfMsMleCrbCleanBoundConsistency"));
  scanSnapshotFile = string(localGetFieldOrDefault(scanData, 'snapshotFile', ""));
  scanElapsedSec = localGetFieldOrDefault(scanData, 'elapsedSec', NaN);
  checkpointDirForReport = string(localGetFieldOrDefault(scanData, 'checkpointDir', ""));

  printMfScanSection('Clean raw aggregate', scanData.cleanAggregateTable);
  printMfScanSection('Clean health aggregate', scanData.cleanHealthAggregateTable);
  printMfScanSection('Clean joint-trim aggregate', scanData.cleanTrimAggregateTable);
  printMfScanSection('Clean MS-vs-SS comparison', scanData.cleanSsMsCompareTable);
  printMfScanSection('Clean known-vs-unknown comparison', scanData.cleanKnownUnknownCompareTable);
  printMfScanSection('Clean SNR curve summary', scanData.cleanSnrCurveTable);
  printMfScanSection('Clean estimate-vs-CRB curve summary', scanData.cleanEstimateCrbCurveTable);
  printMfScanSection('Clean search-range audit', scanData.cleanSearchRangeAuditTable);
  localPrintCompactTable('Clean outlier table', scanData.cleanOutlierTable, 12);
  printMfScanSection('Clean runtime aggregate', scanData.runtimeAggregateTable);
  localPrintCompactTable('Top slow clean tasks', scanData.topSlowRuntimeTable, 12);
  if isfield(scanData, 'checkpointSummaryTable') && istable(scanData.checkpointSummaryTable) && ...
      height(scanData.checkpointSummaryTable) > 0
    printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  end
  fprintf('Scan purpose                    : clean truth-centered CRB-scaled estimator + SNR sweep / trim statistics.\n');
  trimAngleNormMaxForReport = localGetFieldOrDefault(scanConfigForReport, 'trimAngleNormMax', NaN);
  trimFdRefNormMaxForReport = localGetFieldOrDefault(scanConfigForReport, 'trimFdRefNormMax', NaN);
  fprintf('Trim definition                   : health + angle <= %.2f CRB + fdRef <= %.2f CRB.\n', ...
    trimAngleNormMaxForReport, trimFdRefNormMaxForReport);
  fprintf('Main decision table               : cleanSsMsCompareTable.\n');
  fprintf('Known-rate baseline table         : cleanKnownUnknownCompareTable.\n');
  localPlotCleanEstimateVsCrbCurves(scanData.plotData.cleanEstimateCrbCurveTable);

  notifyMfScanStatus(struct( ...
    'scanName', scanNameForReport, ...
    'statusText', "DONE", ...
    'config', scanConfigForReport, ...
    'snapshotFile', scanSnapshotFile, ...
    'checkpointDir', checkpointDirForReport, ...
    'elapsedSec', scanElapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Clean-bound MS/SS CP/IP scan completed; inspect method-specific range settings and aggregate tables before changing estimator logic."]));

catch ME
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', config, ...
    'snapshotFile', "", ...
    'checkpointDir', checkpointRunDir, ...
    'elapsedSec', toc(runTic), ...
    'metricLineList', string(ME.identifier), ...
    'commentLineList', string(ME.message), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers

function scanData = localValidateScanDataForSummary(scanData)
%LOCALVALIDATESCANDATAFORSUMMARY Validate and rebuild stored summary data.

if ~isfield(scanData, 'config') || ~isstruct(scanData.config)
  scanData.config = struct();
end
scanData.config = localCompleteSummaryConfig(scanData.config);
if ~isfield(scanData, 'elapsedSec')
  scanData.elapsedSec = NaN;
end
if ~isfield(scanData, 'snapshotFile')
  scanData.snapshotFile = "";
end
if ~isfield(scanData, 'checkpointDir')
  scanData.checkpointDir = "";
end

hasCaseTable = isfield(scanData, 'caseTable') && istable(scanData.caseTable);
if hasCaseTable
  scanData.caseTable = localEnsureCaseTableSummaryFlags(scanData.caseTable, scanData.config);
  if ~isfield(scanData, 'cleanAggregateTable') || ~istable(scanData.cleanAggregateTable)
    scanData.cleanAggregateTable = buildMfMetricAggregateTable(scanData.caseTable);
  end
  if ~isfield(scanData, 'cleanHealthAggregateTable') || ~istable(scanData.cleanHealthAggregateTable)
    scanData.cleanHealthAggregateTable = buildMfFilteredMetricAggregateTable(scanData.caseTable, ...
      "health", scanData.caseTable.healthResolved);
  end
  if ~isfield(scanData, 'cleanTrimAggregateTable') || ~istable(scanData.cleanTrimAggregateTable)
    scanData.cleanTrimAggregateTable = buildMfFilteredMetricAggregateTable(scanData.caseTable, ...
      "joint-trim", scanData.caseTable.jointTrimKeep);
  end
  if ~isfield(scanData, 'cleanSsMsCompareTable') || ~istable(scanData.cleanSsMsCompareTable)
    scanData.cleanSsMsCompareTable = localBuildCleanSsMsCompareTable(scanData.caseTable, scanData.config);
  end
  if ~isfield(scanData, 'cleanKnownUnknownCompareTable') || ~istable(scanData.cleanKnownUnknownCompareTable)
    scanData.cleanKnownUnknownCompareTable = localBuildCleanKnownUnknownCompareTable(scanData.caseTable);
  end
  if ~isfield(scanData, 'cleanSearchRangeAuditTable') || ~istable(scanData.cleanSearchRangeAuditTable)
    scanData.cleanSearchRangeAuditTable = localBuildCleanSearchRangeAuditTable(scanData.caseTable);
  end
  if ~isfield(scanData, 'cleanOutlierTable') || ~istable(scanData.cleanOutlierTable)
    scanData.cleanOutlierTable = localBuildCleanOutlierTable(scanData.caseTable);
  end
end

if ~isfield(scanData, 'cleanAggregateTable') || ~istable(scanData.cleanAggregateTable)
  error('scanMfMsMleCrbCleanBoundConsistency:MissingAggregateTable', ...
    'scanData.cleanAggregateTable is missing. Store caseTable or aggregate tables inside scanData before saving.');
end
scanData = localEnsureCleanSummaryTable(scanData, 'cleanHealthAggregateTable');
scanData = localEnsureCleanSummaryTable(scanData, 'cleanTrimAggregateTable');
scanData = localEnsureCleanSummaryTable(scanData, 'cleanSsMsCompareTable');
scanData = localEnsureCleanSummaryTable(scanData, 'cleanKnownUnknownCompareTable');
scanData = localEnsureCleanSummaryTable(scanData, 'cleanSearchRangeAuditTable');
scanData = localEnsureCleanSummaryTable(scanData, 'cleanOutlierTable');

if ~isfield(scanData, 'cleanSnrCurveTable') || ~istable(scanData.cleanSnrCurveTable)
  scanData.cleanSnrCurveTable = localBuildCleanSnrCurveTable(scanData.cleanSsMsCompareTable);
end
if ~isfield(scanData, 'cleanEstimateCrbCurveTable') || ~istable(scanData.cleanEstimateCrbCurveTable)
  scanData.cleanEstimateCrbCurveTable = localBuildCleanEstimateCrbCurveTable( ...
    scanData.cleanTrimAggregateTable, scanData.cleanHealthAggregateTable, scanData.cleanAggregateTable);
end

if ~isfield(scanData, 'runtimeAggregateTable') || ~istable(scanData.runtimeAggregateTable)
  if isfield(scanData, 'runtimeTable') && istable(scanData.runtimeTable)
    scanData.runtimeAggregateTable = buildMfRuntimeAggregateTable(scanData.runtimeTable);
  else
    scanData.runtimeAggregateTable = table();
  end
end
if ~isfield(scanData, 'topSlowRuntimeTable') || ~istable(scanData.topSlowRuntimeTable)
  if isfield(scanData, 'runtimeTable') && istable(scanData.runtimeTable)
    scanData.topSlowRuntimeTable = buildMfTopSlowRuntimeTable(scanData.runtimeTable, 12);
  else
    scanData.topSlowRuntimeTable = table();
  end
end

if ~isfield(scanData, 'checkpointSummaryTable') || ~istable(scanData.checkpointSummaryTable)
  checkpointSummary = localGetFieldOrDefault(scanData, 'checkpointSummary', struct());
  scanData.checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointSummary, scanData.checkpointDir);
end
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  scanData.plotData = localBuildCleanPlotData(scanData.cleanEstimateCrbCurveTable);
end
if ~isfield(scanData.plotData, 'cleanEstimateCrbCurveTable') || ...
    ~istable(scanData.plotData.cleanEstimateCrbCurveTable)
  scanData.plotData.cleanEstimateCrbCurveTable = scanData.cleanEstimateCrbCurveTable;
end
end

function caseTable = localEnsureCaseTableSummaryFlags(caseTable, config)
%LOCALENSURECASETABLESUMMARYFLAGS Rebuild summary filters for loaded snapshots.

requiredFlagList = {'healthResolved', 'jointTrimKeep', 'rejectReason'};
hasAllFlags = all(ismember(requiredFlagList, caseTable.Properties.VariableNames));
if hasAllFlags
  return;
end
if localHasCleanHealthConfig(config)
  caseTable = localAnnotateHealthFlags(caseTable, config);
  return;
end
numRow = height(caseTable);
if ~ismember('healthResolved', caseTable.Properties.VariableNames)
  caseTable.healthResolved = true(numRow, 1);
end
if ~ismember('jointTrimKeep', caseTable.Properties.VariableNames)
  caseTable.jointTrimKeep = caseTable.healthResolved;
end
if ~ismember('rejectReason', caseTable.Properties.VariableNames)
  caseTable.rejectReason = strings(numRow, 1);
end
end

function tf = localHasCleanHealthConfig(config)
%LOCALHASCLEANHEALTHCONFIG True when loaded config can rebuild health filters.

fieldList = {'healthFdRefToothFraction', 'healthFdRateAbsMaxHzPerSec', ...
  'trimAngleNormMax', 'trimFdRefNormMax'};
tf = isstruct(config) && all(isfield(config, fieldList));
end

function config = localCompleteSummaryConfig(config)
%LOCALCOMPLETESUMMARYCONFIG Fill summary-only defaults for loaded snapshots.

if ~isfield(config, 'healthFdRefToothFraction')
  config.healthFdRefToothFraction = 0.25;
end
if ~isfield(config, 'healthFdRateAbsMaxHzPerSec')
  config.healthFdRateAbsMaxHzPerSec = 250;
end
if ~isfield(config, 'trimAngleNormMax')
  config.trimAngleNormMax = 5;
end
if ~isfield(config, 'trimFdRefNormMax')
  config.trimFdRefNormMax = 5;
end
if ~isfield(config, 'crbTargetMseRatio')
  config.crbTargetMseRatio = 1.1;
end
end

function scanData = localEnsureCleanSummaryTable(scanData, fieldName)
%LOCALENSURECLEANSUMMARYTABLE Create empty missing optional summary table.

if ~isfield(scanData, fieldName) || ~istable(scanData.(fieldName))
  scanData.(fieldName) = table();
end
end

function taskList = localBuildTaskList(snrDbList, seedList)
%LOCALBUILDTASKLIST Build one task per SNR and seed.

taskList = repmat(struct('taskId', 0, 'snrDb', NaN, 'taskSeed', NaN), 0, 1);
taskId = 0;
for iSnr = 1:numel(snrDbList)
  for iSeed = 1:numel(seedList)
    taskId = taskId + 1;
    taskList(taskId, 1).taskId = taskId; %#ok<AGROW>
    taskList(taskId, 1).snrDb = snrDbList(iSnr); %#ok<AGROW>
    taskList(taskId, 1).taskSeed = seedList(iSeed); %#ok<AGROW>
  end
end
end

function methodRangeTable = localBuildDefaultMethodRangeTable(methodNameList, ...
  doaCrbScale, fdRefRangeMode, fdHalfToothFraction, fdRefCrbScale, fdRateHalfWidthHzPerSec)
%LOCALBUILDDEFAULTMETHODRANGETABLE Build explicit per-method range rows.

methodNameList = reshape(string(methodNameList), [], 1);
rowList = repmat(struct('methodName', "", 'doaCrbScale', NaN, ...
  'fdRefRangeMode', "", 'fdRefHalfToothFraction', NaN, ...
  'fdRefCrbScale', NaN, 'fdRateHalfWidthHzPerSec', NaN), numel(methodNameList), 1);
for iMethod = 1:numel(methodNameList)
  rowList(iMethod).methodName = methodNameList(iMethod);
  rowList(iMethod).doaCrbScale = double(doaCrbScale);
  rowList(iMethod).fdRefRangeMode = localNormalizeFdRefRangeMode(fdRefRangeMode);
  rowList(iMethod).fdRefHalfToothFraction = double(fdHalfToothFraction);
  rowList(iMethod).fdRefCrbScale = double(fdRefCrbScale);
  rowList(iMethod).fdRateHalfWidthHzPerSec = double(fdRateHalfWidthHzPerSec);
end
methodRangeTable = struct2table(rowList(:));
end

function methodRangeTable = localNormalizeMethodRangeTable(methodRangeTable, methodNameList, defaultFdRefRangeMode)
%LOCALNORMALIZEMETHODRANGETABLE Validate and type method-level range columns.

if ~ismember('fdRefRangeMode', methodRangeTable.Properties.VariableNames)
  methodRangeTable.fdRefRangeMode = repmat(localNormalizeFdRefRangeMode(defaultFdRefRangeMode), height(methodRangeTable), 1);
end
requiredNameList = {'methodName', 'doaCrbScale', 'fdRefRangeMode', ...
  'fdRefHalfToothFraction', 'fdRefCrbScale', 'fdRateHalfWidthHzPerSec'};
for iName = 1:numel(requiredNameList)
  if ~ismember(requiredNameList{iName}, methodRangeTable.Properties.VariableNames)
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidMethodRangeTable', ...
      'methodRangeTable is missing required column "%s".', requiredNameList{iName});
  end
end
methodRangeTable.methodName = string(methodRangeTable.methodName);
for iRow = 1:height(methodRangeTable)
  if ~(methodRangeTable.methodName(iRow) == "*" || any(methodNameList == methodRangeTable.methodName(iRow)))
    error('scanMfMsMleCrbCleanBoundConsistency:UnknownMethodRange', ...
      'methodRangeTable contains unknown method "%s".', char(methodRangeTable.methodName(iRow)));
  end
end
for iMethod = 1:numel(methodNameList)
  methodName = methodNameList(iMethod);
  if ~any(methodRangeTable.methodName == methodName) && ~any(methodRangeTable.methodName == "*")
    error('scanMfMsMleCrbCleanBoundConsistency:MissingMethodRange', ...
      'methodRangeTable has no row for method "%s" and no wildcard fallback.', char(methodName));
  end
end
methodRangeTable.doaCrbScale = double(methodRangeTable.doaCrbScale);
fdRefRangeModeList = strings(height(methodRangeTable), 1);
for iRow = 1:height(methodRangeTable)
  fdRefRangeModeList(iRow) = localNormalizeFdRefRangeMode(methodRangeTable.fdRefRangeMode(iRow));
end
methodRangeTable.fdRefRangeMode = fdRefRangeModeList;
methodRangeTable.fdRefHalfToothFraction = double(methodRangeTable.fdRefHalfToothFraction);
methodRangeTable.fdRefCrbScale = double(methodRangeTable.fdRefCrbScale);
methodRangeTable.fdRateHalfWidthHzPerSec = double(methodRangeTable.fdRateHalfWidthHzPerSec);
end

function mode = localNormalizeFdRefRangeMode(mode)
%LOCALNORMALIZEFDREFRANGEMODE Normalize fdRef local range mode.

mode = lower(strtrim(string(mode)));
if mode == "tooth" || mode == "tooth-crb-floor"
  mode = "tooth-floor";
end
if ~(mode == "tooth-floor" || mode == "crb-scale")
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidFdRefRangeMode', ...
    'fdRefRangeMode must be "tooth-floor" or "crb-scale".');
end
end

function rangeCfg = localResolveMethodRangeConfig(config, methodName, crbBundle, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz)
%LOCALRESOLVEMETHODRANGECONFIG Resolve one method's local search range.

methodName = string(methodName);
methodRangeTable = config.methodRangeTable;
maskExact = methodRangeTable.methodName == methodName;
maskFallback = methodRangeTable.methodName == "*";
if any(maskExact)
  row = methodRangeTable(find(maskExact, 1, 'last'), :);
elseif any(maskFallback)
  row = methodRangeTable(find(maskFallback, 1, 'last'), :);
else
  error('scanMfMsMleCrbCleanBoundConsistency:MissingMethodRange', ...
    'No methodRangeTable row for method "%s".', char(methodName));
end
fdRefRangeMode = localNormalizeFdRefRangeMode(row.fdRefRangeMode);
if ~(isfinite(row.doaCrbScale) && row.doaCrbScale > 0 && ...
    isfinite(row.fdRefHalfToothFraction) && row.fdRefHalfToothFraction > 0 && ...
    isfinite(row.fdRefCrbScale) && row.fdRefCrbScale > 0 && ...
    isfinite(row.fdRateHalfWidthHzPerSec) && row.fdRateHalfWidthHzPerSec > 0)
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidMethodRangeValue', ...
    'methodRangeTable row for "%s" has nonpositive or nonfinite bound values.', ...
    char(methodName));
end
fdHalfWidthHz = double(row.fdRefHalfToothFraction) * toothStepHz;
fdRangeRequested = truthFdRefHz + [-fdHalfWidthHz, fdHalfWidthHz];
if fdRefRangeMode == "crb-scale"
  fdRangeUse = localBuildFdRefCrbScaledRange(truthFdRefHz, crbBundle, methodName, ...
    double(row.fdRefCrbScale), toothStepHz, ...
    localGetFieldOrDefault(config, 'defaultFdRefMaxHalfToothFraction', 0.45));
else
  fdRangeUse = localApplyFdRefCrbScaledFloor(fdRangeRequested, truthFdRefHz, ...
    crbBundle, methodName, double(row.fdRefCrbScale));
end
rangeCfg = struct();
rangeCfg.methodName = methodName;
rangeCfg.doaCrbScale = double(row.doaCrbScale);
rangeCfg.fdRefRangeMode = fdRefRangeMode;
rangeCfg.fdRefHalfToothFraction = localSafeRatio(0.5 * abs(fdRangeUse(2) - fdRangeUse(1)), toothStepHz);
rangeCfg.fdRefRequestedHalfToothFraction = double(row.fdRefHalfToothFraction);
rangeCfg.fdRefCrbScale = double(row.fdRefCrbScale);
rangeCfg.fdRateHalfWidthHzPerSec = double(row.fdRateHalfWidthHzPerSec);
rangeCfg.fdRangeRequested = fdRangeRequested;
rangeCfg.fdRangeUse = fdRangeUse;
rangeCfg.fdRateRangeUse = truthFdRateHzPerSec + double(row.fdRateHalfWidthHzPerSec) * [-1, 1];
end

function [repeatCell, runState] = localRunTaskBatch(context, flowOpt, config, taskList, scanName)
%LOCALRUNTASKBATCH Run clean tasks with optional checkpoint/resume.

numTask = numel(taskList);
runState = struct();
progressTracker = localCreateProgressTracker(numTask);
progressCleanup = onCleanup(@() localCloseProgressTracker(progressTracker)); %#ok<NASGU>
try
  if logical(localGetFieldOrDefault(config, 'checkpointEnable', false))
    checkpointOpt = localBuildCheckpointOpt(scanName, config, localGetFieldOrDefault(config, 'useParfor', false));
    taskDir = fullfile(checkpointOpt.runDir, 'task');
    numDoneTask = localCountCheckpointTaskFile(taskDir, numTask);
    localAdvanceProgressByStep(progressTracker, numDoneTask);
    sharedData = struct('context', context, 'flowOpt', flowOpt, 'config', config);
    runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker);
    runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, @localCheckpointTaskRunner, runnerOpt);
    repeatCell = runState.resultCell;
    return;
  end

  repeatCell = cell(numTask, 1);
  useParfor = logical(localGetFieldOrDefault(config, 'useParfor', false));
  if useParfor && isfield(progressTracker, 'queue') && ~isempty(progressTracker.queue)
    progressQueue = progressTracker.queue;
    parfor iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
      send(progressQueue, 1);
    end
  elseif useParfor
    parfor iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
    end
  else
    for iTask = 1:numTask
      repeatCell{iTask} = localRunOneTask(context, flowOpt, config, taskList(iTask));
      localAdvanceProgress(progressTracker);
    end
  end
catch ME
  localCloseProgressTracker(progressTracker);
  rethrow(ME);
end
end

function runnerOpt = localBuildCheckpointRunnerOpt(checkpointOpt, progressTracker)
%LOCALBUILDCHECKPOINTRUNNEROPT Keep only fields accepted by checkpoint runner.

runnerOpt = struct();
runnerOpt.runName = checkpointOpt.runName;
runnerOpt.runKey = checkpointOpt.runKey;
runnerOpt.outputRoot = checkpointOpt.outputRoot;
runnerOpt.useParfor = logical(localGetFieldOrDefault(checkpointOpt, 'useParfor', false));
runnerOpt.resume = checkpointOpt.resume;
runnerOpt.meta = checkpointOpt.meta;
runnerOpt.progressFcn = @(step) localAdvanceProgressByStep(progressTracker, step);
runnerOpt.cleanupOnSuccess = false;
runnerOpt.cleanupOpt = struct();
end

function numDone = localCountCheckpointTaskFile(taskDir, numTask)
%LOCALCOUNTCHECKPOINTTASKFILE Count completed checkpoint task files before progress reset.

numDone = 0;
if ~isfolder(taskDir)
  return;
end
for iTask = 1:numTask
  if isfile(fullfile(taskDir, sprintf('task_%06d.mat', iTask)))
    numDone = numDone + 1;
  end
end
end

function localAdvanceProgressByStep(progressTracker, step)
%LOCALADVANCEPROGRESSBYSTEP Advance progressbar by a task count.

for iStep = 1:max(0, round(step))
  localAdvanceProgress(progressTracker);
end
end

function repeatOut = localCheckpointTaskRunner(taskInfo, sharedData)
%LOCALCHECKPOINTTASKRUNNER Run one checkpointed clean task.

repeatOut = localRunOneTask(sharedData.context, sharedData.flowOpt, sharedData.config, taskInfo);
end

function repeatOut = localRunOneTask(context, flowOpt, config, task)
%LOCALRUNONETASK Run one SNR/seed clean truth-centered task.

lastwarn('', '');
tRepeat = tic;
runtimeRows = repmat(emptyMfRuntimeRow(), 0, 1);

tStage = tic;
repeatData = buildDynamicRepeatData(context, task.snrDb, task.taskSeed);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "repeat-data", "setup", 1000 * toc(tStage)); %#ok<AGROW>
fixture = repeatData.periodicFixture;
truth = fixture.truth;
toothStepHz = localResolveToothStepHz(fixture);
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  error('scanMfMsMleCrbCleanBoundConsistency:MissingTruthFdLine', ...
    'Truth fdRef/fdRate values are required for this truth-centered clean scan.');
end

noiseVar = 1 / (10^(task.snrDb / 10));
tStage = tic;
crbBundle = localBuildCrbBundleQuiet(fixture, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, noiseVar, config);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "crb-bundle", "crb", 1000 * toc(tStage)); %#ok<AGROW>

tStage = tic;
staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz, config);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, "", "static-bundle", "static", 1000 * toc(tStage)); %#ok<AGROW>

caseRowList = repmat(localEmptyCaseRow(), numel(config.methodNameList), 1);
caseMap = struct();
for iMethod = 1:numel(config.methodNameList)
  methodName = config.methodNameList(iMethod);
  rangeCfg = localResolveMethodRangeConfig(config, methodName, crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  tMethod = tic;
  lastwarn('', '');
  methodRuntimeRows = repmat(emptyMfRuntimeRow(), 0, 1);
  if methodName == "SS-SF-DoA"
    caseUse = staticBundle.caseSfDoaRefOnly;
    initMeta = staticBundle.initMetaSfDoaRefOnly;
  elseif methodName == "MS-SF-DoA"
    caseUse = staticBundle.caseSfDoaMs;
    initMeta = staticBundle.initMetaSfDoaMs;
  elseif methodName == "SS-SF-Static"
    caseUse = staticBundle.caseSfStaticRefOnly;
    initMeta = staticBundle.initMetaSfStaticRefOnly;
  elseif methodName == "MS-SF-Static"
    caseUse = staticBundle.caseSfStaticMs;
    initMeta = staticBundle.initMetaSfStaticMs;
  elseif methodName == "SS-MF-Static"
    caseUse = staticBundle.caseMfStaticRefOnly;
    initMeta = staticBundle.initMetaMfStaticRefOnly;
  elseif methodName == "MS-MF-Static"
    caseUse = staticBundle.caseMfStaticMs;
    initMeta = staticBundle.initMetaMfStaticMs;
  else
    method = localBuildDynamicMethod(methodName, rangeCfg.doaCrbScale, config.dynamicRouteMode);
    method = localApplyDoaCrbScaledWidth(method, crbBundle);
    staticSeedCase = localSelectStaticSeedCase(method, staticBundle);
    [caseUse, initMeta, methodRuntimeRows] = localRunDynamicMethod(method, fixture, truth, context, flowOpt, ...
      rangeCfg.fdRangeUse, rangeCfg.fdRateRangeUse, toothStepHz, config.optVerbose, staticSeedCase, config, task);
  end
  runtimeRows = [runtimeRows; methodRuntimeRows(:)]; %#ok<AGROW>
  if endsWith(methodName, "SF-DoA")
    fdRangeForRow = [NaN, NaN];
  else
    fdRangeForRow = rangeCfg.fdRangeUse;
  end
  row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
    fdRangeForRow, rangeCfg.fdRateRangeUse, crbBundle, config, initMeta, rangeCfg);
  [warningMessage, warningId] = lastwarn;
  row.warningFlag = ~(isempty(warningMessage) && isempty(warningId));
  caseRowList(iMethod) = row;
  caseMap.(localMethodKey(methodName)) = caseUse; %#ok<NASGU>
  runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, methodName, "method-total", "method", 1000 * toc(tMethod)); %#ok<AGROW>
end

[warningMessage, warningId] = lastwarn;
repeatOut = struct();
repeatOut.taskSeed = task.taskSeed;
repeatOut.snrDb = task.snrDb;
repeatOut.truthFdRefHz = truthFdRefHz;
repeatOut.truthFdRateHzPerSec = truthFdRateHzPerSec;
repeatOut.toothStepHz = toothStepHz;
repeatOut.caseRowList = caseRowList;
repeatOut.runtimeRowList = runtimeRows;
repeatOut.repeatTotalMs = 1000 * toc(tRepeat);
repeatOut.warningSeen = ~(isempty(warningMessage) && isempty(warningId));
repeatOut.warningId = string(warningId);
repeatOut.warningMessage = string(warningMessage);
end

function staticBundle = localBuildCleanStaticBundle(fixture, context, flowOpt, crbBundle, truth, ...
  truthFdRefHz, truthFdRateHzPerSec, toothStepHz, config)
%LOCALBUILDCLEANSTATICBUNDLE Run SF DoA-only and MF static seeds for clean init.

truthDoaParam = reshape(truth.latlonTrueDeg, [], 1);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundConsistency:MissingTruthDoa', ...
    'Truth lat/lon is required to build the clean truth-centered search box.');
end
truthDoaParam = truthDoaParam(1:2);

methodNameList = reshape(string(config.methodNameList), [], 1);
needRef = any(methodNameList == "SS-SF-DoA" | methodNameList == "SS-MF-Static" | startsWith(methodNameList, "SS-MF-"));
needMs = any(methodNameList == "MS-SF-DoA" | methodNameList == "MS-MF-Static" | startsWith(methodNameList, "MS-MF-"));

staticBundle = struct();
if needRef
  rangeCfgSfDoaRef = localResolveMethodRangeConfig(config, localFirstExistingMethod(config, ["SS-SF-DoA"; "SS-MF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "SS-MF-IP-K"; "SS-MF-IP-U"]), ...
    crbBundle, truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  [staticBundle.caseSfDoaRefOnly, staticBundle.initMetaSfDoaRefOnly] = localRunCleanSfDoaCase( ...
    "SS-SF-DoA", "single", fixture.viewRefOnly, context, crbBundle, truthDoaParam, ...
    rangeCfgSfDoaRef.doaCrbScale, config);
end
if needMs
  rangeCfgSfDoaMs = localResolveMethodRangeConfig(config, localFirstExistingMethod(config, ["MS-SF-DoA"; "MS-MF-Static"; "MS-MF-CP-K"; "MS-MF-CP-U"; "MS-MF-IP-K"; "MS-MF-IP-U"]), ...
    crbBundle, truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  [staticBundle.caseSfDoaMs, staticBundle.initMetaSfDoaMs] = localRunCleanSfDoaCase( ...
    "MS-SF-DoA", "multi", fixture.viewMs, context, crbBundle, truthDoaParam, ...
    rangeCfgSfDoaMs.doaCrbScale, config);
end

if any(methodNameList == "SS-SF-Static")
  rangeCfgRefSfStatic = localResolveMethodRangeConfig(config, "SS-SF-Static", crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  [staticBundle.caseSfStaticRefOnly, staticBundle.initMetaSfStaticRefOnly] = localRunCleanStaticCase( ...
    "SS-SF-Static", "single", fixture.viewRefOnly, context, flowOpt.staticBaseOpt, ...
    crbBundle, truthDoaParam, rangeCfgRefSfStatic.fdRangeUse, rangeCfgRefSfStatic.doaCrbScale, config);
end
if any(methodNameList == "MS-SF-Static")
  rangeCfgMsSfStatic = localResolveMethodRangeConfig(config, "MS-SF-Static", crbBundle, ...
    truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  [staticBundle.caseSfStaticMs, staticBundle.initMetaSfStaticMs] = localRunCleanStaticCase( ...
    "MS-SF-Static", "multi", fixture.viewMs, context, flowOpt.staticBaseOpt, ...
    crbBundle, truthDoaParam, rangeCfgMsSfStatic.fdRangeUse, rangeCfgMsSfStatic.doaCrbScale, config);
end

if any(methodNameList == "SS-MF-Static" | startsWith(methodNameList, "SS-MF-"))
  rangeCfgRefMfStatic = localResolveMethodRangeConfig(config, localFirstExistingMethod(config, ["SS-MF-Static"; "SS-MF-CP-K"; "SS-MF-CP-U"; "SS-MF-IP-K"; "SS-MF-IP-U"]), ...
    crbBundle, truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  methodRefStatic = localBuildDynamicMethod("SS-MF-Static", rangeCfgRefMfStatic.doaCrbScale, config.dynamicRouteMode);
  methodRefStatic = localApplyDoaCrbScaledWidth(methodRefStatic, crbBundle);
  [staticBundle.caseMfStaticRefOnly, staticBundle.initMetaMfStaticRefOnly] = localRunCleanMfStaticCase( ...
    methodRefStatic, fixture.viewRefOnly, truth, context, flowOpt.dynBaseOpt, rangeCfgRefMfStatic.fdRangeUse, ...
    staticBundle.caseSfDoaRefOnly, config);
end
if any(methodNameList == "MS-MF-Static" | startsWith(methodNameList, "MS-MF-"))
  rangeCfgMsMfStatic = localResolveMethodRangeConfig(config, localFirstExistingMethod(config, ["MS-MF-Static"; "MS-MF-CP-K"; "MS-MF-CP-U"; "MS-MF-IP-K"; "MS-MF-IP-U"]), ...
    crbBundle, truthFdRefHz, truthFdRateHzPerSec, toothStepHz);
  methodMsStatic = localBuildDynamicMethod("MS-MF-Static", rangeCfgMsMfStatic.doaCrbScale, config.dynamicRouteMode);
  methodMsStatic = localApplyDoaCrbScaledWidth(methodMsStatic, crbBundle);
  [staticBundle.caseMfStaticMs, staticBundle.initMetaMfStaticMs] = localRunCleanMfStaticCase( ...
    methodMsStatic, fixture.viewMs, truth, context, flowOpt.dynBaseOpt, rangeCfgMsMfStatic.fdRangeUse, ...
    staticBundle.caseSfDoaMs, config);
end

if isfield(staticBundle, 'caseSfDoaRefOnly')
  staticBundle.caseStaticRefOnly = staticBundle.caseSfDoaRefOnly;
end
if isfield(staticBundle, 'caseSfDoaMs')
  staticBundle.caseStaticMs = staticBundle.caseSfDoaMs;
end
if isfield(staticBundle, 'caseMfStaticRefOnly')
  staticBundle.bestStaticRefOnlyCase = staticBundle.caseMfStaticRefOnly;
end
if isfield(staticBundle, 'caseMfStaticMs')
  staticBundle.bestStaticMsCase = staticBundle.caseMfStaticMs;
end
end

function methodName = localFirstExistingMethod(config, candidateList)
%LOCALFIRSTEXISTINGMETHOD Return the first configured method from a candidate list.

candidateList = reshape(string(candidateList), [], 1);
methodNameList = reshape(string(config.methodNameList), [], 1);
for iCandidate = 1:numel(candidateList)
  if any(methodNameList == candidateList(iCandidate))
    methodName = candidateList(iCandidate);
    return;
  end
end
methodName = candidateList(1);
end

function [caseInfo, initMeta] = localRunCleanSfDoaCase(displayName, satMode, viewUse, context, crbBundle, truthDoaParam, doaCrbScale, config)
%LOCALRUNCLEANSFDOACASE Run one SF DoA-only estimator inside the clean DoA box.

modelOpt = struct();
modelOpt.initMethod = 'midpoint';
doaHalfWidth = localResolveDoaCrbScaledHalfWidth(displayName, doaCrbScale, crbBundle);
doaGridUse = localBuildTruthCenteredDoaGrid(viewUse.doaGrid, truthDoaParam, doaHalfWidth);
initParam = localBuildCrbOffsetDoaInit(truthDoaParam, doaHalfWidth);
[estResult, ~, ~] = estimatorDoaMlePilotOpt( ...
  viewUse.sceneRef.array, context.wavelen, viewUse.rxSigSf, context.pilotWave, ...
  doaGridUse, initParam, config.optVerbose, modelOpt);
estResult = localWrapCleanSfDoaOnlyResult(estResult);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa", "no-doppler", estResult);
initMeta = localBuildSfDoaInitMeta(estResult, truthDoaParam, doaHalfWidth, initParam);
end

function initParam = localBuildCrbOffsetDoaInit(truthDoaParam, doaHalfWidth)
%LOCALBUILDCRBOFFSETDOAINIT Build a deterministic non-truth CRB-box DoA seed.

truthDoaParam = reshape(double(truthDoaParam), [], 1);
doaHalfWidth = localNormalizeDoaHalfWidth(doaHalfWidth);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2))) || ...
    ~localIsValidPositiveDoaHalfWidth(doaHalfWidth)
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidSfDoaCrbInit', ...
    'A positive finite CRB-scaled DoA box is required for SF DoA-only initialization.');
end
offsetSign = [1; -1];
initParam = truthDoaParam(1:2) + 0.5 * offsetSign .* doaHalfWidth(1:2);
end

function estResult = localWrapCleanSfDoaOnlyResult(estResult)
%LOCALWRAPCLEANSFDOAONLYRESULT Add DoA-only fields expected by shared summary code.

estResult.modelType = 'sf-doa-only';
estResult.phaseMode = 'single-frame';
estResult.fdRateMode = 'none';
estResult.fdRefEst = NaN;
estResult.fdRateEst = NaN;
estResult.fdLineRmseHz = NaN;
estResult.fdSatRmseHz = NaN;
if ~isfield(estResult, 'aux') || ~isstruct(estResult.aux)
  estResult.aux = struct();
end
if isfield(estResult, 'latlonEst') && ~isempty(estResult.latlonEst)
  estResult.aux.latlonEst = estResult.latlonEst;
end
if isfield(estResult, 'localDoaEst') && ~isempty(estResult.localDoaEst)
  estResult.aux.localDoaEst = estResult.localDoaEst;
end
end

function initMeta = localBuildSfDoaInitMeta(estResult, truthDoaParam, doaHalfWidth, requestedInitParam)
%LOCALBUILDSFDOAINITMETA Record DoA-only initializer and local box.

initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', requestedInitParam(:)), [], 1);
requestedInitParam = reshape(double(requestedInitParam), [], 1);
initMeta = struct();
initMeta.requestedInitMode = "truth-crb-scaled-sf-doa-crb-offset";
initMeta.requestedInitParam = requestedInitParam(:);
initMeta.requestedInitDoaParam = truthDoaParam(:);
initMeta.initParam = initParam(:);
initMeta.initDoaParam = requestedInitParam(:);
if numel(initParam) >= 2 && all(isfinite(initParam(1:2)))
  initMeta.initDoaParam = initParam(1:2);
end
initMeta.initFdRefHz = NaN;
initMeta.initFdRateHzPerSec = NaN;
initMeta.searchHalfWidthDeg = reshape(doaHalfWidth, [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
initMeta.seedDoaParam = initMeta.initDoaParam(:);
initMeta.boundCenterDoaParam = truthDoaParam(:);
end

function doaGridUse = localBuildTruthCenteredDoaGrid(doaGridIn, truthDoaParam, doaHalfWidth)
%LOCALBUILDTRUTHCENTEREDDOAGRID Restrict an existing lat/lon grid to the clean truth box.

if iscell(doaGridIn)
  doaGridUse = cell(size(doaGridIn));
  for iGrid = 1:numel(doaGridIn)
    doaGridUse{iGrid} = localBuildTruthCenteredDoaGrid(doaGridIn{iGrid}, truthDoaParam, doaHalfWidth);
  end
  return;
end
truthDoaParam = reshape(double(truthDoaParam), [], 1);
doaHalfWidth = localNormalizeDoaHalfWidth(doaHalfWidth);
if ~localIsValidPositiveDoaHalfWidth(doaHalfWidth)
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidDoaGridHalfWidth', ...
    'A positive finite DoA half-width is required to build the clean truth-centered grid.');
end
localRange = [truthDoaParam(1) - doaHalfWidth(1), truthDoaParam(1) + doaHalfWidth(1); ...
  truthDoaParam(2) - doaHalfWidth(2), truthDoaParam(2) + doaHalfWidth(2)];
doaGridUse = doaGridIn;
if isfield(doaGridUse, 'range')
  doaGridUse.range = localRange;
end
end

function [caseInfo, initMeta] = localRunCleanMfStaticCase(method, viewUse, truth, context, dynBaseOpt, fdRangeUse, sfSeedCase, config)
%LOCALRUNCLEANMFSTATICCASE Run a zero-rate MF DoA-Doppler estimator seed.

dynOpt = dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.fdRateMode = 'zero';
dynOpt.steeringMode = char(config.cleanMfSteeringMode);
dynOpt.steeringRefFrameIdx = viewUse.sceneSeq.refFrameIdx;
dynOpt = localApplyDynamicObjectiveMode(dynOpt, config);
dynOpt = localApplyDynamicRouteMode(dynOpt, method.dynamicRouteMode);
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
[initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, sfSeedCase, truth, config);
dynOpt.initDoaParam = initRequestMeta.initDoaParam(:);
dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
initRequestMeta.initDoaHalfWidth = method.doaHalfWidthDeg(:);
[estResult, ~, ~] = estimatorDoaDopplerMlePilotMfStatOpt( ...
  viewUse.sceneSeq, viewUse.rxSigMf, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, viewUse.doaGrid, fdRangeUse, initParamOverride, ...
  config.optVerbose, dynOpt);
caseInfo = buildDoaDopplerCaseResult(method.displayName, method.satMode, "multi", ...
  "doa-doppler", "mf-static", estResult);
caseInfo.inToothMeta = struct('fdRangeUse', fdRangeUse, 'fdRateRangeUse', []);
initMeta = localBuildDynamicInitMeta(caseInfo, method, initRequestMeta);
end

function [caseInfo, initMeta] = localRunCleanStaticCase(displayName, satMode, viewUse, context, staticBaseOpt, ...
  crbBundle, truthDoaParam, fdRangeUse, doaCrbScale, config)
%LOCALRUNCLEANSTATICCASE Run one SF static estimator inside the CRB-scaled truth-centered box.

modelOpt = staticBaseOpt;
modelOpt.initDoaParam = truthDoaParam(:);
modelOpt.initDoaHalfWidth = localResolveDoaCrbScaledHalfWidth(displayName, doaCrbScale, crbBundle);
[estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  viewUse.sceneRef, viewUse.rxSigSf, context.pilotWave, context.carrierFreq, ...
  context.waveInfo.sampleRate, viewUse.doaGrid, fdRangeUse, 1, [], ...
  config.optVerbose, modelOpt);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa-doppler", "static", estResult);
initMeta = localBuildStaticInitMeta(estResult, modelOpt);
end

function initMeta = localBuildStaticInitMeta(estResult, modelOpt)
%LOCALBUILDSTATICINITMETA Record static estimator initial state and local box.

initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
initMeta = struct();
initMeta.requestedInitMode = "truth-crb-scaled-static";
initMeta.requestedInitParam = initParam(:);
initMeta.requestedInitDoaParam = reshape(modelOpt.initDoaParam, [], 1);
initMeta.initParam = initParam(:);
initMeta.initDoaParam = reshape(localGetFieldOrDefault(modelOpt, 'initDoaParam', NaN(2, 1)), [], 1);
if numel(initParam) >= 2 && all(isfinite(initParam(1:2)))
  initMeta.initDoaParam = initParam(1:2);
end
initMeta.initFdRefHz = NaN;
if numel(initParam) >= 3
  initMeta.initFdRefHz = initParam(3);
end
initMeta.initFdRateHzPerSec = NaN;
initMeta.searchHalfWidthDeg = reshape(localGetFieldOrDefault(modelOpt, 'initDoaHalfWidth', NaN(2, 1)), [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
end

function [caseUse, initMeta, runtimeRows] = localRunDynamicMethod(method, fixture, truth, context, flowOpt, ...
  fdRangeUse, fdRateRangeUse, toothStepHz, optVerbose, staticSeedCase, config, task)
%LOCALRUNDYNAMICMETHOD Run one clean MF estimator path.

runtimeRows = repmat(emptyMfRuntimeRow(), 0, 1);
if string(method.satMode) == "multi"
  viewUse = fixture.viewMs;
  debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthMs', struct());
else
  viewUse = fixture.viewRefOnly;
  debugTruthUse = localGetFieldOrDefault(fixture, 'debugTruthRef', struct());
end

fdRateRangeForMethod = fdRateRangeUse;
if method.isKnownRate || string(method.fdRateMode) == "zero"
  fdRateRangeForMethod = [];
end

dynOpt = flowOpt.dynBaseOpt;
dynOpt.phaseMode = char(method.phaseMode);
dynOpt.fdRateMode = char(method.fdRateMode);
if string(method.fdRateMode) == "known"
  dynOpt.fdRateKnown = truth.fdRateFit;
end
dynOpt.steeringMode = char(config.cleanMfSteeringMode);
dynOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynOpt = localApplyDynamicObjectiveMode(dynOpt, config);
dynOpt = localApplyDynamicRouteMode(dynOpt, method.dynamicRouteMode);
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = method.doaHalfWidthDeg(:);
[initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, staticSeedCase, truth, config);
if isfield(initRequestMeta, 'initDoaParam') && ~isempty(initRequestMeta.initDoaParam)
  dynOpt.initDoaParam = initRequestMeta.initDoaParam(:);
  dynOpt.initDoaHalfWidth = method.doaHalfWidthDeg(:);
  initRequestMeta.initDoaHalfWidth = method.doaHalfWidthDeg(:);
end

tStage = tic;
caseUse = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
  viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  fdRangeUse, fdRateRangeForMethod, optVerbose, dynOpt, method.isKnownRate, ...
  debugTruthUse, initParamOverride);
runtimeRows(end + 1, 1) = buildMfRuntimeRow(task, method.displayName, "baseline-solve", "estimator", 1000 * toc(tStage)); %#ok<AGROW>
caseUse.inToothMeta = struct('toothStepHz', toothStepHz, 'fdRangeUse', fdRangeUse, ...
  'fdRateRangeUse', fdRateRangeForMethod);
initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta);
end

function mode = localNormalizeDynamicObjectiveMode(mode)
%LOCALNORMALIZEDYNAMICOBJECTIVEMODE Normalize the dynamic objective dialect tag.

mode = lower(strtrim(string(mode)));
switch mode
  case {"pure-mle", "puremle", "pure"}
    mode = "pure-mle";
  case {"flow-default", "flowdefault", "robust", "engineering"}
    mode = "flow-default";
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidDynamicObjectiveMode', ...
      'Unsupported dynamicObjectiveMode: %s.', char(mode));
end
end

function dynOpt = localApplyDynamicObjectiveMode(dynOpt, config)
%LOCALAPPLYDYNAMICOBJECTIVEMODE Apply pure-MLE or robust objective options.

switch string(config.dynamicObjectiveMode)
  case "pure-mle"
    dynOpt.continuousPhaseConsistencyWeight = 0;
    dynOpt.continuousPhaseCollapsePenaltyWeight = 0;
    dynOpt.continuousPhaseNegativeProjectionPenaltyWeight = 0;
    dynOpt.continuousPhaseNonRefFitFloorWeight = 0;
    dynOpt.continuousPhaseNonRefSupportFloorWeight = 0;
    dynOpt.enableFdAliasUnwrap = false;
  case "flow-default"
    dynOpt.enableFdAliasUnwrap = true;
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidDynamicObjectiveMode', ...
      'Unsupported dynamic objective mode: %s.', char(config.dynamicObjectiveMode));
end
end

function mode = localNormalizeDynamicRouteMode(mode)
%LOCALNORMALIZEDYNAMICROUTEMODE Normalize the estimator route tag.

mode = lower(strtrim(string(mode)));
switch mode
  case {"core-only", "coreonly", "core"}
    mode = "core-only";
  case {"robust-estimator", "robustestimator", "robust", "default"}
    mode = "robust-estimator";
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidDynamicRouteMode', ...
      'Unsupported dynamicRouteMode: %s.', char(mode));
end
end

function dynOpt = localApplyDynamicRouteMode(dynOpt, routeMode)
%LOCALAPPLYDYNAMICROUTEMODE Apply one-shot core or robust branch route.

routeMode = localNormalizeDynamicRouteMode(routeMode);
switch string(routeMode)
  case "core-only"
    dynOpt.disableDoaBasinEntry = true;
    dynOpt.doaBasinEntryScope = 'off';
    dynOpt.disableKnownEmbedded = true;
    dynOpt.disableUnknownWarmAnchor = true;
  case "robust-estimator"
    dynOpt.disableDoaBasinEntry = false;
    dynOpt.disableKnownEmbedded = false;
    dynOpt.disableUnknownWarmAnchor = false;
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidDynamicRouteMode', ...
      'Unsupported dynamic route mode: %s.', char(routeMode));
end
end


function mode = localNormalizeCleanMfSignalTimeModel(mode)
%LOCALNORMALIZECLEANMFSIGNALTIMEMODE Normalize the scan signal phase model.

mode = lower(strtrim(string(mode)));
switch mode
  case {"linearref", "linear-ref", "linear_ref"}
    mode = "linearRef";
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidSignalTimeModel', ...
      'Unsupported cleanMfSignalTimeModel: %s.', char(mode));
end
end

function mode = localNormalizeCleanMfSteeringMode(mode)
%LOCALNORMALIZECLEANMFSTEERINGMODE Normalize the scan steering model.

mode = lower(strtrim(string(mode)));
switch mode
  case {"frozenref", "frozen-ref", "frozen_ref"}
    mode = "frozenRef";
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:InvalidSteeringMode', ...
      'Unsupported cleanMfSteeringMode: %s.', char(mode));
end
end

function context = localApplyCleanSignalModelOverride(context, config)
%LOCALAPPLYCLEANSIGNALMODELOVERRIDE Force the scan paper signal model.

if ~isfield(context, 'simOpt') || isempty(context.simOpt)
  context.simOpt = struct();
end
if isfield(context.simOpt, 'spatial')
  context.simOpt = rmfield(context.simOpt, 'spatial');
end
if ~isfield(context.simOpt, 'phase') || isempty(context.simOpt.phase)
  context.simOpt.phase = struct();
end
context.simOpt.phase.timeModel = char(config.cleanMfSignalTimeModel);
context.simOpt.phase.frameModel = 'shared';
context.simOpt.steeringMode = char(config.cleanMfSteeringMode);
context.simOpt.accessMode = 'framewise';
context.simOpt.steeringRefFrameIdx = context.sceneSeqMaster.refFrameIdx;
if ~isfield(context.simOpt, 'wave') || isempty(context.simOpt.wave)
  context.simOpt.wave = struct();
end
context.simOpt.wave.delayModel = 'phaseOnly';
context.simOpt.wave.timeRef = 'zero';
context.simOpt.wave.carrierPhaseModel = 'none';
end

function method = localBuildDynamicMethod(methodName, doaCrbScale, defaultRouteMode)
%LOCALBUILDDYNAMICMETHOD Build one SS/MS MF method descriptor.

if nargin < 3 || isempty(defaultRouteMode)
  defaultRouteMode = "core-only";
end
methodName = string(methodName);
baseMethodName = localBaseDynamicMethodName(methodName);
routeMode = localNormalizeDynamicRouteMode(defaultRouteMode);
if endsWith(methodName, "-Robust")
  routeMode = "robust-estimator";
elseif endsWith(methodName, "-Core")
  routeMode = "core-only";
end

if endsWith(baseMethodName, "MF-Static")
  fdRateMode = "zero";
  isKnownRate = true;
elseif contains(baseMethodName, "-K")
  fdRateMode = "known";
  isKnownRate = true;
elseif contains(baseMethodName, "-U")
  fdRateMode = "unknown";
  isKnownRate = false;
else
  error('scanMfMsMleCrbCleanBoundConsistency:UnknownMethod', 'Unsupported method name "%s".', char(methodName));
end
if endsWith(baseMethodName, "MF-Static") || contains(baseMethodName, "-CP-")
  phaseMode = "continuous";
elseif contains(baseMethodName, "-IP-")
  phaseMode = "independent";
else
  error('scanMfMsMleCrbCleanBoundConsistency:UnknownPhaseMethod', 'Unsupported method phase in "%s".', char(methodName));
end
method = struct();
method.displayName = methodName;
method.baseDisplayName = baseMethodName;
method.dynamicRouteMode = routeMode;
if startsWith(baseMethodName, "MS-")
  method.satMode = "multi";
else
  method.satMode = "single";
end
method.phaseMode = phaseMode;
method.fdRateMode = fdRateMode;
method.isKnownRate = logical(isKnownRate);
method.doaHalfWidthDeg = reshape(double(doaCrbScale), [], 1);
end

function baseMethodName = localBaseDynamicMethodName(methodName)
%LOCALBASEDYNAMICMETHODNAME Strip replay-only route suffix from a method name.

baseMethodName = string(methodName);
baseMethodName = erase(baseMethodName, "-Robust");
baseMethodName = erase(baseMethodName, "-Core");
end

function staticSeedCase = localSelectStaticSeedCase(method, staticBundle)
%LOCALSELECTSTATICSEEDCASE Select the matching static seed for SS/MS MF.

if string(method.satMode) == "multi"
  if string(method.fdRateMode) == "zero"
    staticSeedCase = staticBundle.caseSfDoaMs;
  else
    staticSeedCase = staticBundle.caseMfStaticMs;
  end
else
  if string(method.fdRateMode) == "zero"
    staticSeedCase = staticBundle.caseSfDoaRefOnly;
  else
    staticSeedCase = staticBundle.caseMfStaticRefOnly;
  end
end
end

function [initParamOverride, initRequestMeta] = localBuildDynamicInitOverride(method, staticSeedCase, truth, config)
%LOCALBUILDDYNAMICINITOVERRIDE Build a static seed with a fixed CRB-scaled MF box.

initMode = "truth-crb-scaled-static-seed-fixed-bound";
initParamOverride = [];
initRequestMeta = struct('initMode', initMode, 'initParam', [], ...
  'initDoaParam', [], 'initDoaHalfWidth', [], 'seedDoaParam', [], ...
  'seedFdRefHz', NaN, 'boundCenterDoaParam', []);

truthDoaParam = reshape(truth.latlonTrueDeg, [], 1);
if numel(truthDoaParam) < 2 || any(~isfinite(truthDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundConsistency:MissingTruthDoa', ...
    'Truth DoA is required to build the fixed clean MF search box.');
end
truthDoaParam = truthDoaParam(1:2);

estResult = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
seedDoaParam = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(seedDoaParam) < 2 || any(~isfinite(seedDoaParam(1:2)))
  error('scanMfMsMleCrbCleanBoundConsistency:MissingStaticDoaSeed', ...
    'A finite truth-centered static DoA seed is required for the MF refinement.');
end
seedDoaParam = seedDoaParam(1:2);
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
seedFdRefHz = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
if ~(isfinite(seedFdRefHz))
  seedFdRefHz = truthFdRefHz;
end

% The override is the solver seed. The DoA bound center below intentionally
% stays at truth so that static seed errors cannot expand the clean envelope.
if string(method.fdRateMode) == "unknown"
  initParamOverride = [seedDoaParam(:); seedFdRefHz; truthFdRateHzPerSec];
else
  initParamOverride = [seedDoaParam(:); seedFdRefHz];
end
initRequestMeta.initParam = initParamOverride(:);
initRequestMeta.initDoaParam = truthDoaParam(:);
initRequestMeta.initDoaHalfWidth = method.doaHalfWidthDeg(:);
initRequestMeta.seedDoaParam = seedDoaParam(:);
initRequestMeta.seedFdRefHz = seedFdRefHz;
initRequestMeta.boundCenterDoaParam = truthDoaParam(:);
end

function initMeta = localBuildDynamicInitMeta(caseUse, method, initRequestMeta)
%LOCALBUILDDYNAMICINITMETA Extract requested and internal initializer.

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
initParam = reshape(localGetFieldOrDefault(estResult, 'initParam', []), [], 1);
initMeta = struct();
initMeta.requestedInitMode = string(localGetFieldOrDefault(initRequestMeta, 'initMode', ""));
initMeta.requestedInitParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initParam', []), [], 1);
initMeta.requestedInitDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaParam', []), [], 1);
initMeta.seedDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'seedDoaParam', []), [], 1);
initMeta.seedFdRefHz = localGetFieldOrDefault(initRequestMeta, 'seedFdRefHz', NaN);
initMeta.boundCenterDoaParam = reshape(localGetFieldOrDefault(initRequestMeta, 'boundCenterDoaParam', []), [], 1);
initMeta.initParam = initParam(:);
initMeta.initDoaParam = NaN(2, 1);
initMeta.initFdRefHz = NaN;
initMeta.initFdRateHzPerSec = NaN;
if numel(initParam) >= 2
  initMeta.initDoaParam = initParam(1:2);
end
if numel(initParam) >= 3
  initMeta.initFdRefHz = initParam(3);
end
if string(method.fdRateMode) == "unknown" && numel(initParam) >= 4
  initMeta.initFdRateHzPerSec = initParam(4);
end
initMeta.dynamicRouteMode = string(localGetFieldOrDefault(method, 'dynamicRouteMode', ""));
initMeta.baseDisplayName = string(localGetFieldOrDefault(method, 'baseDisplayName', method.displayName));
initMeta.searchHalfWidthDeg = reshape(localGetFieldOrDefault(initRequestMeta, 'initDoaHalfWidth', NaN(2, 1)), [], 1);
if numel(initMeta.searchHalfWidthDeg) == 1
  initMeta.searchHalfWidthDeg = repmat(initMeta.searchHalfWidthDeg, 2, 1);
end
end

function row = localBuildCaseRow(methodName, caseUse, truth, task, fixture, toothStepHz, ...
  fdRangeUse, fdRateRangeUse, crbBundle, config, initMeta, rangeCfg)
%LOCALBUILDCASEMETRICROW Convert one clean case into a compact metric row.

row = localEmptyCaseRow();
methodName = string(methodName);
info = summarizeDoaDopplerCase(caseUse, truth);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
initEval = localGetFieldOrDefault(estDebug, 'initEval', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());

row.displayName = methodName;
row.dynamicObjectiveMode = string(config.dynamicObjectiveMode);
row.dynamicRouteMode = string(localGetFieldOrDefault(initMeta, 'dynamicRouteMode', config.dynamicRouteMode));
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.doaCrbScale = localGetFieldOrDefault(rangeCfg, 'doaCrbScale', NaN);
row.fdRefRangeMode = string(localGetFieldOrDefault(rangeCfg, 'fdRefRangeMode', ""));
row.fdRefHalfToothFraction = localGetFieldOrDefault(rangeCfg, 'fdRefHalfToothFraction', NaN);
row.fdRefCrbScale = localGetFieldOrDefault(rangeCfg, 'fdRefCrbScale', NaN);
row.fdRateHalfWidthHzPerSec = localGetFieldOrDefault(rangeCfg, 'fdRateHalfWidthHzPerSec', NaN);
row.numFrame = fixture.sceneSeq.numFrame;
row.frameIntvlSec = median(diff(fixture.sceneSeq.timeOffsetSec));
row.windowSec = (row.numFrame - 1) * row.frameIntvlSec;
row.toothStepHz = toothStepHz;
row.truthLatDeg = truth.latlonTrueDeg(1);
row.truthLonDeg = truth.latlonTrueDeg(2);
row.finalLatDeg = info.latEstDeg;
row.finalLonDeg = info.lonEstDeg;
row.latErrDeg = info.latErrDeg;
row.lonErrDeg = info.lonErrDeg;
row.sphericalAngleErrDeg = info.angleErrDeg;
row.finalAngleErrDeg = info.angleErrDeg;
row.fdRefTrueHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
row.fdRateTrueHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
row.finalFdRefHz = info.fdRefEstHz;
row.fdRefErrHz = info.fdRefErrHz;
row.fdRefAbsErrHz = abs(info.fdRefErrHz);
row.finalFdRateHzPerSec = info.fdRateEstHzPerSec;
row.fdRateErrHzPerSec = info.fdRateErrHzPerSec;
row.fdRateAbsErrHzPerSec = abs(info.fdRateErrHzPerSec);
row.fdLineRmseHz = info.fdLineRmseHz;
row.solverResolved = logical(info.isResolved);
row.exitflag = info.exitflag;
row.iterations = info.iterations;
row.funcCount = info.funcCount;
row.runTimeMs = info.runTimeMs;

if endsWith(methodName, "SF-DoA")
  row.satMode = string(ternary(startsWith(methodName, "MS-"), "multi", "single"));
  row.frameMode = "single";
  row.phaseMode = "single-frame";
  row.fdRateMode = "none";
  row.initLatDeg = initMeta.initDoaParam(1);
  row.initLonDeg = initMeta.initDoaParam(2);
  row.initAngleErrDeg = calcLatlonAngleError(initMeta.initDoaParam(:), truth.latlonTrueDeg(:));
  row.finalMinusInitAngleDeg = calcLatlonAngleError([row.finalLatDeg; row.finalLonDeg], initMeta.initDoaParam(:));
  row.initFdRefHz = NaN;
  row.initFdRateHzPerSec = NaN;
  row.searchHalfWidthLatDeg = initMeta.searchHalfWidthDeg(1);
  row.searchHalfWidthLonDeg = initMeta.searchHalfWidthDeg(min(2, numel(initMeta.searchHalfWidthDeg)));
  row.fdRateRangeLowHzPerSec = NaN;
  row.fdRateRangeHighHzPerSec = NaN;
  row.mfInitMode = initMeta.requestedInitMode;
elseif endsWith(methodName, "SF-Static")
  row.satMode = string(ternary(startsWith(methodName, "MS-"), "multi", "single"));
  row.frameMode = "single";
  row.phaseMode = "single-frame";
  row.fdRateMode = "zero";
  row.initLatDeg = initMeta.initDoaParam(1);
  row.initLonDeg = initMeta.initDoaParam(2);
  row.initAngleErrDeg = calcLatlonAngleError(initMeta.initDoaParam(:), truth.latlonTrueDeg(:));
  row.finalMinusInitAngleDeg = calcLatlonAngleError([row.finalLatDeg; row.finalLonDeg], initMeta.initDoaParam(:));
  row.initFdRefHz = initMeta.initFdRefHz;
  row.initFdRateHzPerSec = NaN;
  row.searchHalfWidthLatDeg = initMeta.searchHalfWidthDeg(1);
  row.searchHalfWidthLonDeg = initMeta.searchHalfWidthDeg(min(2, numel(initMeta.searchHalfWidthDeg)));
  row.fdRateRangeLowHzPerSec = NaN;
  row.fdRateRangeHighHzPerSec = NaN;
  row.mfInitMode = initMeta.requestedInitMode;
else
  baseMethodName = localBaseDynamicMethodName(methodName);
  row.satMode = string(ternary(startsWith(baseMethodName, "MS-"), "multi", "single"));
  row.frameMode = "multi";
  if contains(baseMethodName, "-IP-")
    row.phaseMode = "independent";
  else
    row.phaseMode = "continuous";
  end
  if endsWith(baseMethodName, "MF-Static")
    row.fdRateMode = "zero";
  else
    row.fdRateMode = string(ternary(endsWith(baseMethodName, "-K"), "known", "unknown"));
  end
  row.initLatDeg = initMeta.initDoaParam(1);
  row.initLonDeg = initMeta.initDoaParam(2);
  row.initAngleErrDeg = calcLatlonAngleError(initMeta.initDoaParam(:), truth.latlonTrueDeg(:));
  row.finalMinusInitAngleDeg = calcLatlonAngleError([row.finalLatDeg; row.finalLonDeg], initMeta.initDoaParam(:));
  row.initFdRefHz = initMeta.initFdRefHz;
  row.initFdRateHzPerSec = initMeta.initFdRateHzPerSec;
  row.searchHalfWidthLatDeg = initMeta.searchHalfWidthDeg(1);
  row.searchHalfWidthLonDeg = initMeta.searchHalfWidthDeg(min(2, numel(initMeta.searchHalfWidthDeg)));
  row.mfInitMode = initMeta.requestedInitMode;
end

row.seedLatDeg = NaN;
row.seedLonDeg = NaN;
seedDoaParam = reshape(localGetFieldOrDefault(initMeta, 'seedDoaParam', []), [], 1);
if numel(seedDoaParam) >= 2 && all(isfinite(seedDoaParam(1:2)))
  row.seedLatDeg = seedDoaParam(1);
  row.seedLonDeg = seedDoaParam(2);
end
row.boundCenterLatDeg = row.truthLatDeg;
row.boundCenterLonDeg = row.truthLonDeg;
boundCenterDoaParam = reshape(localGetFieldOrDefault(initMeta, 'boundCenterDoaParam', []), [], 1);
if numel(boundCenterDoaParam) >= 2 && all(isfinite(boundCenterDoaParam(1:2)))
  row.boundCenterLatDeg = boundCenterDoaParam(1);
  row.boundCenterLonDeg = boundCenterDoaParam(2);
end
boxTolDeg = 1e-10;
row.doaLatBoxViolation = abs(row.finalLatDeg - row.boundCenterLatDeg) > row.searchHalfWidthLatDeg + boxTolDeg;
row.doaLonBoxViolation = abs(row.finalLonDeg - row.boundCenterLonDeg) > row.searchHalfWidthLonDeg + boxTolDeg;
row.doaTruthBoxViolation = row.doaLatBoxViolation || row.doaLonBoxViolation;

[crbDoa, crbFull] = localResolveFiniteCrbForMethod(crbBundle, methodName);
[angleCrbStdDeg, angleAux] = projectCrbToAngleMetric(crbDoa, truth.latlonTrueDeg, 'latlon');
row.crbLatStdDeg = sqrt(max(real(crbDoa(1, 1)), 0));
row.crbLonStdDeg = sqrt(max(real(crbDoa(2, 2)), 0));
row.crbTraceStdDeg = sqrt(max(real(trace(crbDoa)), 0));
row.crbSphericalApproxStdDeg = angleCrbStdDeg;
row.crbAngleMetricVarRad2 = angleAux.angleVarRad2;
row.crbMetricTraceRatio = localSafeRatio(row.crbTraceStdDeg, row.crbSphericalApproxStdDeg);
row.angleErrOverSphericalCrb = localSafeRatio(row.sphericalAngleErrDeg, row.crbSphericalApproxStdDeg);
row.angleErrOverTraceCrb = localSafeRatio(row.sphericalAngleErrDeg, row.crbTraceStdDeg);
row.fdRefCrbStdHz = NaN;
if size(crbFull, 1) >= 3 && size(crbFull, 2) >= 3
  row.fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
end
row.fdRefErrOverCrb = localSafeRatio(row.fdRefAbsErrHz, row.fdRefCrbStdHz);

row.fdRangeLowHz = fdRangeUse(1);
row.fdRangeHighHz = fdRangeUse(2);
row.fdRangeHalfWidthHz = 0.5 * abs(fdRangeUse(2) - fdRangeUse(1));
if row.frameMode == "multi" && row.fdRateMode == "unknown"
  row.fdRateRangeLowHzPerSec = fdRateRangeUse(1);
  row.fdRateRangeHighHzPerSec = fdRateRangeUse(2);
end
row.fdRefBoundaryHit = localIsNearRangeBoundary(row.finalFdRefHz, fdRangeUse, config.boundaryTolFraction);
row.fdRateBoundaryHit = row.frameMode == "multi" && row.fdRateMode == "unknown" && ...
  localIsNearRangeBoundary(row.finalFdRateHzPerSec, fdRateRangeUse, config.boundaryTolFraction);

row.initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objectiveImprove = row.initObj - row.finalObj;
row.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
[row.refCoherence, row.nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval);
row.stepSize = localGetFieldOrDefault(optimInfo, 'stepsize', NaN);
row.firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
row.constrViolation = localGetFieldOrDefault(optimInfo, 'constrviolation', NaN);
row.algorithm = string(localGetFieldOrDefault(optimInfo, 'algorithm', ""));
row.usedScaledSolve = logical(localGetFieldOrDefault(optimInfo, 'usedScaledSolve', false));
row.fallbackTriggered = logical(localGetFieldOrDefault(optimInfo, 'fallbackTriggered', false));
row.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));
row.candidateCount = numel(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.candidateVariantText = strjoin(cellstr(reshape(string(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1))), 1, [])), ',');
row.fdRefInitMoveHz = row.finalFdRefHz - row.initFdRefHz;
row.fdRateInitMoveHzPerSec = row.finalFdRateHzPerSec - row.initFdRateHzPerSec;
end

function caseTable = localBuildCaseTable(repeatCell)
%LOCALBUILDCASETABLE Collect all per-method rows into one table.

numRepeat = numel(repeatCell);
rowCountList = zeros(numRepeat, 1);
for iRepeat = 1:numRepeat
  repeatOut = repeatCell{iRepeat};
  if isempty(repeatOut) || ~isfield(repeatOut, 'caseRowList') || isempty(repeatOut.caseRowList)
    continue;
  end
  rowCountList(iRepeat) = numel(repeatOut.caseRowList);
end

totalRow = sum(rowCountList);
if totalRow == 0
  caseTable = struct2table(repmat(localEmptyCaseRow(), 0, 1));
  return;
end

rowList = repmat(localEmptyCaseRow(), totalRow, 1);
rowOffset = 0;
for iRepeat = 1:numRepeat
  numRow = rowCountList(iRepeat);
  if numRow == 0
    continue;
  end
  rowIdx = rowOffset + (1:numRow);
  rowList(rowIdx) = repeatCell{iRepeat}.caseRowList(:);
  rowOffset = rowOffset + numRow;
end
caseTable = struct2table(rowList(:));
end

function caseTable = localAnnotateHealthFlags(caseTable, config)
%LOCALANNOTATEHEALTHFLAGS Add clean health and joint-trim flags.

if isempty(caseTable) || height(caseTable) == 0
  return;
end
numRow = height(caseTable);
isDynamic = contains(caseTable.displayName, "-MF-");
isUnknown = endsWith(caseTable.displayName, "-U") & isDynamic;

solverHealthy = true(numRow, 1);
solverHealthy(isDynamic) = isfinite(caseTable.finalObj(isDynamic)) & ...
  isfinite(caseTable.iterations(isDynamic)) & caseTable.iterations(isDynamic) > 0 & ...
  isfinite(caseTable.objectiveImprove(isDynamic));

freqHealthy = true(numRow, 1);
fdRefThresholdHz = config.healthFdRefToothFraction .* caseTable.toothStepHz;
freqHealthy(isDynamic) = caseTable.fdRefAbsErrHz(isDynamic) <= fdRefThresholdHz(isDynamic);

rateHealthy = true(numRow, 1);
rateHealthy(isUnknown) = caseTable.fdRateAbsErrHzPerSec(isUnknown) <= config.healthFdRateAbsMaxHzPerSec;

boundaryHealthy = true(numRow, 1);
boundaryHealthy(isDynamic) = ~(caseTable.fdRefBoundaryHit(isDynamic) | caseTable.fdRateBoundaryHit(isDynamic));

notNoSolve = true(numRow, 1);
noSolve = isDynamic & (caseTable.iterations <= 0 | ~isfinite(caseTable.iterations)) & ...
  (~isfinite(caseTable.objectiveImprove) | abs(caseTable.objectiveImprove) <= eps);
notNoSolve(noSolve) = false;

healthResolved = solverHealthy & freqHealthy & rateHealthy & boundaryHealthy & notNoSolve;
angleNormOk = ~isDynamic | caseTable.angleErrOverSphericalCrb <= config.trimAngleNormMax;
fdRefNormOk = ~isDynamic | caseTable.fdRefErrOverCrb <= config.trimFdRefNormMax;
trimKeep = healthResolved & angleNormOk;
jointTrimKeep = healthResolved & angleNormOk & fdRefNormOk;
fdRefNormTail = isDynamic & healthResolved & angleNormOk & ~fdRefNormOk;

caseTable.isDynamicMethod = isDynamic;
caseTable.solverHealthy = solverHealthy;
caseTable.freqHealthy = freqHealthy;
caseTable.rateHealthy = rateHealthy;
caseTable.boundaryHealthy = boundaryHealthy;
caseTable.notNoSolve = notNoSolve;
caseTable.healthResolved = healthResolved;
caseTable.trimKeep = trimKeep;
caseTable.jointTrimKeep = jointTrimKeep;
caseTable.fdRefNormTail = fdRefNormTail;
caseTable.rejectReason = localBuildRejectReason(caseTable, angleNormOk, fdRefNormOk);
end

function rejectReason = localBuildRejectReason(caseTable, angleNormOk, fdRefNormOk)
%LOCALBUILDREJECTREASON Build compact health / trim reason tags.

numRow = height(caseTable);
rejectReason = strings(numRow, 1);
for iRow = 1:numRow
  if ~caseTable.isDynamicMethod(iRow)
    rejectReason(iRow) = "baseline";
    continue;
  end
  reasonList = strings(0, 1);
  if ~caseTable.solverHealthy(iRow)
    reasonList(end + 1, 1) = "solver"; %#ok<AGROW>
  end
  if ~caseTable.freqHealthy(iRow)
    reasonList(end + 1, 1) = "fdRef"; %#ok<AGROW>
  end
  if ~caseTable.rateHealthy(iRow)
    reasonList(end + 1, 1) = "fdRate"; %#ok<AGROW>
  end
  if ~caseTable.boundaryHealthy(iRow)
    reasonList(end + 1, 1) = "boundary"; %#ok<AGROW>
  end
  if ~caseTable.notNoSolve(iRow)
    reasonList(end + 1, 1) = "no-solve"; %#ok<AGROW>
  end
  if caseTable.healthResolved(iRow) && ~angleNormOk(iRow)
    reasonList(end + 1, 1) = "trim-angle"; %#ok<AGROW>
  end
  if caseTable.healthResolved(iRow) && ~fdRefNormOk(iRow)
    reasonList(end + 1, 1) = "trim-fdRef"; %#ok<AGROW>
  end
  if isempty(reasonList)
    rejectReason(iRow) = "kept";
  else
    rejectReason(iRow) = strjoin(cellstr(reasonList), "+");
  end
end
end

function compareTable = localBuildCleanSsMsCompareTable(caseTable, config)
%LOCALBUILDCLEANSSMSCOMPARETABLE Compare configured SS-MF and MS-MF pairs.

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'phaseMode', "", ...
  'fdRateMode', "", 'ssName', "", 'msName', "", ...
  'ssKeepRate', NaN, 'msKeepRate', NaN, ...
  'ssAngleRmseDeg', NaN, 'msAngleRmseDeg', NaN, 'msAbsGainOverSs', NaN, ...
  'ssAngleMseOverOwnCrb', NaN, 'msAngleMseOverOwnCrb', NaN, 'msAngleMseOverSsCrb', NaN, ...
  'ssFdRefMseOverOwnCrb', NaN, 'msFdRefMseOverOwnCrb', NaN, ...
  'ssFdRateRmseHzPerSec', NaN, 'msFdRateRmseHzPerSec', NaN, ...
  'msBeatsSsAbs', false, 'msOwnCrbTargetPass', false, 'targetClass', "");
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  compareTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
methodNameList = unique(string(caseTable.displayName), 'stable');
ssNameList = methodNameList(startsWith(methodNameList, "SS-MF-"));
pairList = repmat(struct('ssName', "", 'msName', "", 'phaseMode', "", 'fdRateMode', ""), 0, 1);
for iPair = 1:numel(ssNameList)
  ssName = ssNameList(iPair);
  msName = "MS-" + extractAfter(ssName, "SS-");
  if ~any(methodNameList == msName)
    continue;
  end
  ssFirst = find(caseTable.displayName == ssName, 1, 'first');
  rowPair = struct('ssName', ssName, 'msName', msName, ...
    'phaseMode', string(caseTable.phaseMode(ssFirst)), ...
    'fdRateMode', string(caseTable.fdRateMode(ssFirst)));
  pairList(end + 1, 1) = rowPair; %#ok<AGROW>
end
if isempty(pairList)
  compareTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
filterNameList = ["raw"; "health"; "joint-trim"];
snrList = unique(caseTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList) * numel(pairList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iPair = 1:numel(pairList)
    pair = pairList(iPair);
    ssRaw = caseTable.displayName == pair.ssName & caseTable.snrDb == snrDb;
    msRaw = caseTable.displayName == pair.msName & caseTable.snrDb == snrDb;
    for iFilter = 1:numel(filterNameList)
      filterName = filterNameList(iFilter);
      idx = idx + 1;
      row = rowTemplate;
      row.filterName = filterName;
      row.snrDb = snrDb;
      row.phaseMode = pair.phaseMode;
      row.fdRateMode = pair.fdRateMode;
      row.ssName = pair.ssName;
      row.msName = pair.msName;
      ssMask = ssRaw & localFilterMask(caseTable, filterName);
      msMask = msRaw & localFilterMask(caseTable, filterName);
      row.ssKeepRate = localSafeRatio(nnz(ssMask), nnz(ssRaw));
      row.msKeepRate = localSafeRatio(nnz(msMask), nnz(msRaw));
      row.ssAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(ssMask));
      row.msAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(msMask));
      row.msAbsGainOverSs = localSafeRatio(row.ssAngleRmseDeg - row.msAngleRmseDeg, row.ssAngleRmseDeg);
      row.ssAngleMseOverOwnCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(ssMask), caseTable.crbSphericalApproxStdDeg(ssMask));
      row.msAngleMseOverOwnCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(msMask), caseTable.crbSphericalApproxStdDeg(msMask));
      row.msAngleMseOverSsCrb = localMsErrorOverMatchedSsCrb(caseTable, msMask, ssRaw);
      row.ssFdRefMseOverOwnCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(ssMask), caseTable.fdRefCrbStdHz(ssMask));
      row.msFdRefMseOverOwnCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(msMask), caseTable.fdRefCrbStdHz(msMask));
      row.ssFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(ssMask));
      row.msFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(msMask));
      row.msBeatsSsAbs = isfinite(row.msAbsGainOverSs) && row.msAbsGainOverSs > 0;
      row.msOwnCrbTargetPass = row.msAngleMseOverOwnCrb <= config.crbTargetMseRatio && ...
        row.msFdRefMseOverOwnCrb <= config.crbTargetMseRatio;
      if row.msOwnCrbTargetPass
        row.targetClass = "crb-level";
      elseif row.msBeatsSsAbs
        row.targetClass = "beats-ss-only";
      else
        row.targetClass = "not-yet";
      end
      rowList(idx) = row;
    end
  end
end
compareTable = struct2table(rowList(:));
end

function compareTable = localBuildCleanKnownUnknownCompareTable(caseTable)
%LOCALBUILDCLEANKNOWNUNKNOWNCOMPARETABLE Compare configured K/U clean MF pairs.

rowTemplate = struct('filterName', "", 'snrDb', NaN, 'satMode', "", ...
  'phaseMode', "", 'knownName', "", 'unknownName', "", ...
  'knownKeepRate', NaN, 'unknownKeepRate', NaN, ...
  'knownAngleRmseDeg', NaN, 'unknownAngleRmseDeg', NaN, ...
  'unknownRmseGainOverKnown', NaN, ...
  'knownAngleMseOverCrb', NaN, 'unknownAngleMseOverCrb', NaN, ...
  'unknownMinusKnownAngleMseOverCrb', NaN, ...
  'knownFdRefMseOverCrb', NaN, 'unknownFdRefMseOverCrb', NaN, ...
  'unknownMinusKnownFdRefMseOverCrb', NaN, ...
  'knownFdRefRmseHz', NaN, 'unknownFdRefRmseHz', NaN, ...
  'knownFdRateRmseHzPerSec', NaN, 'unknownFdRateRmseHzPerSec', NaN, ...
  'comparisonClass', "");
if isempty(caseTable) || ~istable(caseTable) || height(caseTable) == 0
  compareTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
methodNameList = unique(string(caseTable.displayName), 'stable');
knownNameList = methodNameList(endsWith(methodNameList, "-K") & contains(methodNameList, "-MF-"));
pairList = repmat(struct('knownName', "", 'unknownName', "", 'satMode', "", 'phaseMode', ""), 0, 1);
for iPair = 1:numel(knownNameList)
  knownName = knownNameList(iPair);
  unknownName = replace(knownName, "-K", "-U");
  if ~any(methodNameList == unknownName)
    continue;
  end
  knownFirst = find(caseTable.displayName == knownName, 1, 'first');
  rowPair = struct('knownName', knownName, 'unknownName', unknownName, ...
    'satMode', string(caseTable.satMode(knownFirst)), ...
    'phaseMode', string(caseTable.phaseMode(knownFirst)));
  pairList(end + 1, 1) = rowPair; %#ok<AGROW>
end
if isempty(pairList)
  compareTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
filterNameList = ["raw"; "health"; "joint-trim"];
snrList = unique(caseTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(filterNameList) * numel(snrList) * numel(pairList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iPair = 1:numel(pairList)
    pair = pairList(iPair);
    knownRaw = caseTable.displayName == pair.knownName & caseTable.snrDb == snrDb;
    unknownRaw = caseTable.displayName == pair.unknownName & caseTable.snrDb == snrDb;
    for iFilter = 1:numel(filterNameList)
      filterName = filterNameList(iFilter);
      idx = idx + 1;
      row = rowTemplate;
      row.filterName = filterName;
      row.snrDb = snrDb;
      row.satMode = pair.satMode;
      row.phaseMode = pair.phaseMode;
      row.knownName = pair.knownName;
      row.unknownName = pair.unknownName;
      knownMask = knownRaw & localFilterMask(caseTable, filterName);
      unknownMask = unknownRaw & localFilterMask(caseTable, filterName);
      row.knownKeepRate = localSafeRatio(nnz(knownMask), nnz(knownRaw));
      row.unknownKeepRate = localSafeRatio(nnz(unknownMask), nnz(unknownRaw));
      row.knownAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(knownMask));
      row.unknownAngleRmseDeg = localRmse(caseTable.sphericalAngleErrDeg(unknownMask));
      row.unknownRmseGainOverKnown = localSafeRatio(row.knownAngleRmseDeg - row.unknownAngleRmseDeg, row.knownAngleRmseDeg);
      row.knownAngleMseOverCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(knownMask), caseTable.crbSphericalApproxStdDeg(knownMask));
      row.unknownAngleMseOverCrb = localNormalizedMse(caseTable.sphericalAngleErrDeg(unknownMask), caseTable.crbSphericalApproxStdDeg(unknownMask));
      row.unknownMinusKnownAngleMseOverCrb = row.unknownAngleMseOverCrb - row.knownAngleMseOverCrb;
      row.knownFdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(knownMask), caseTable.fdRefCrbStdHz(knownMask));
      row.unknownFdRefMseOverCrb = localNormalizedMse(caseTable.fdRefAbsErrHz(unknownMask), caseTable.fdRefCrbStdHz(unknownMask));
      row.unknownMinusKnownFdRefMseOverCrb = row.unknownFdRefMseOverCrb - row.knownFdRefMseOverCrb;
      row.knownFdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(knownMask));
      row.unknownFdRefRmseHz = localRmse(caseTable.fdRefAbsErrHz(unknownMask));
      row.knownFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(knownMask));
      row.unknownFdRateRmseHzPerSec = localRmse(caseTable.fdRateAbsErrHzPerSec(unknownMask));
      row.comparisonClass = localClassifyKnownUnknownComparison(row);
      rowList(idx) = row;
    end
  end
end
compareTable = struct2table(rowList(:));
end

function classText = localClassifyKnownUnknownComparison(row)
%LOCALCLASSIFYKNOWNUNKNOWNCOMPARISON Classify finite-MC CP-K/CP-U comparison.

if ~(isfinite(row.knownAngleRmseDeg) && isfinite(row.unknownAngleRmseDeg))
  classText = "missing";
  return;
end
relTol = 0.02;
if row.unknownAngleRmseDeg < (1 - relTol) * row.knownAngleRmseDeg
  classText = "unknown-smaller-rmse";
elseif row.knownAngleRmseDeg < (1 - relTol) * row.unknownAngleRmseDeg
  classText = "known-smaller-rmse";
else
  classText = "similar";
end
end

function curveTable = localBuildCleanSnrCurveTable(compareTable)
%LOCALBUILDCLEANSNRCURVETABLE Keep compact joint-trim SNR metrics for summary.

rowTemplate = struct('snrDb', NaN, 'phaseMode', "", 'fdRateMode', "", ...
  'ssName', "", 'msName', "", ...
  'ssAngleMseOverCrb', NaN, 'msAngleMseOverCrb', NaN, ...
  'ssFdRefMseOverCrb', NaN, 'msFdRefMseOverCrb', NaN, 'ssAngleRmseDeg', NaN, ...
  'msAngleRmseDeg', NaN, 'ssKeepRate', NaN, 'msKeepRate', NaN, ...
  'msAbsGainOverSs', NaN, 'targetClass', "");
if isempty(compareTable) || ~istable(compareTable) || height(compareTable) == 0
  curveTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
mask = compareTable.filterName == "joint-trim";
if ~any(mask)
  mask = compareTable.filterName == "raw";
end
sortNameList = {'snrDb'};
if ismember('phaseMode', compareTable.Properties.VariableNames)
  sortNameList{end + 1} = 'phaseMode';
end
if ismember('fdRateMode', compareTable.Properties.VariableNames)
  sortNameList{end + 1} = 'fdRateMode';
end
curveSource = sortrows(compareTable(mask, :), sortNameList);
rowList = repmat(rowTemplate, height(curveSource), 1);
for iRow = 1:height(curveSource)
  row = rowTemplate;
  row.snrDb = curveSource.snrDb(iRow);
  if ismember('phaseMode', curveSource.Properties.VariableNames)
    row.phaseMode = curveSource.phaseMode(iRow);
  end
  if ismember('fdRateMode', curveSource.Properties.VariableNames)
    row.fdRateMode = curveSource.fdRateMode(iRow);
  end
  if ismember('ssName', curveSource.Properties.VariableNames)
    row.ssName = curveSource.ssName(iRow);
  end
  if ismember('msName', curveSource.Properties.VariableNames)
    row.msName = curveSource.msName(iRow);
  end
  row.ssAngleMseOverCrb = curveSource.ssAngleMseOverOwnCrb(iRow);
  row.msAngleMseOverCrb = curveSource.msAngleMseOverOwnCrb(iRow);
  row.ssFdRefMseOverCrb = curveSource.ssFdRefMseOverOwnCrb(iRow);
  row.msFdRefMseOverCrb = curveSource.msFdRefMseOverOwnCrb(iRow);
  row.ssAngleRmseDeg = curveSource.ssAngleRmseDeg(iRow);
  row.msAngleRmseDeg = curveSource.msAngleRmseDeg(iRow);
  row.ssKeepRate = curveSource.ssKeepRate(iRow);
  row.msKeepRate = curveSource.msKeepRate(iRow);
  row.msAbsGainOverSs = curveSource.msAbsGainOverSs(iRow);
  row.targetClass = curveSource.targetClass(iRow);
  rowList(iRow) = row;
end
curveTable = struct2table(rowList(:));
end

function curveTable = localBuildCleanEstimateCrbCurveTable(trimAggregateTable, healthAggregateTable, rawAggregateTable)
%LOCALBUILDCLEANESTIMATECRBCURVETABLE Build absolute RMSE/CRB curves for plotting.

rowTemplate = struct('filterName', "", 'displayName', "", 'snrDb', NaN, ...
  'angleRmseDeg', NaN, 'angleCrbSphericalMedianDeg', NaN, ...
  'fdRefRmseHz', NaN, 'fdRefCrbMedianHz', NaN);
sourceTable = trimAggregateTable;
filterName = "joint-trim";
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  sourceTable = healthAggregateTable;
  filterName = "health";
end
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  sourceTable = rawAggregateTable;
  filterName = "raw";
end
if isempty(sourceTable) || ~istable(sourceTable) || height(sourceTable) == 0
  curveTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
sourceTable = sortrows(sourceTable, {'snrDb', 'displayName'});
rowList = repmat(rowTemplate, height(sourceTable), 1);
for iRow = 1:height(sourceTable)
  row = rowTemplate;
  row.filterName = filterName;
  row.displayName = string(sourceTable.displayName(iRow));
  row.snrDb = sourceTable.snrDb(iRow);
  row.angleRmseDeg = sourceTable.angleRmseDeg(iRow);
  row.angleCrbSphericalMedianDeg = sourceTable.angleCrbSphericalMedianDeg(iRow);
  row.fdRefRmseHz = sourceTable.fdRefRmseHz(iRow);
  row.fdRefCrbMedianHz = sourceTable.fdRefCrbMedianHz(iRow);
  rowList(iRow) = row;
end
curveTable = struct2table(rowList(:));
end

function plotData = localBuildCleanPlotData(cleanEstimateCrbCurveTable)
%LOCALBUILDCLEANPLOTDATA Store all summary plot inputs in scanData.

plotData = struct();
plotData.cleanEstimateCrbCurveTable = cleanEstimateCrbCurveTable;
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointSummary, checkpointDir)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Convert stored checkpoint summary to a table.

checkpointSummaryTable = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
  zeros(0, 1), zeros(0, 1), false(0, 1), false(0, 1), ...
  'VariableNames', {'runName', 'runKey', 'runDir', 'numTask', 'numDone', ...
  'isComplete', 'cleanedOnSuccess'});
if nargin < 1 || ~isstruct(checkpointSummary) || isempty(checkpointSummary)
  return;
end
runNameList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runName', strings(0, 1))), [], 1);
runKeyList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runKey', strings(0, 1))), [], 1);
runDirList = reshape(string(localGetFieldOrDefault(checkpointSummary, 'runDir', strings(0, 1))), [], 1);
numTaskList = reshape(double(localGetFieldOrDefault(checkpointSummary, 'numTask', zeros(0, 1))), [], 1);
numDoneList = reshape(double(localGetFieldOrDefault(checkpointSummary, 'numDone', zeros(0, 1))), [], 1);
isCompleteList = reshape(logical(localGetFieldOrDefault(checkpointSummary, 'isComplete', false(0, 1))), [], 1);
if isempty(runNameList) && nargin >= 2 && strlength(string(checkpointDir)) > 0
  runNameList = "";
  runKeyList = "";
  runDirList = string(checkpointDir);
  numTaskList = 0;
  numDoneList = 0;
  isCompleteList = false;
end
numRow = max([numel(runNameList), numel(runKeyList), numel(runDirList), ...
  numel(numTaskList), numel(numDoneList), numel(isCompleteList)]);
if numRow == 0
  return;
end
runNameList = localPadColumn(runNameList, numRow, "");
runKeyList = localPadColumn(runKeyList, numRow, "");
runDirList = localPadColumn(runDirList, numRow, "");
numTaskList = localPadColumn(numTaskList, numRow, 0);
numDoneList = localPadColumn(numDoneList, numRow, 0);
isCompleteList = localPadColumn(isCompleteList, numRow, false);
cleanedOnSuccess = false(numRow, 1);
if isfield(checkpointSummary, 'cleanedOnSuccess') && ~isempty(checkpointSummary.cleanedOnSuccess)
  cleanedOnSuccess(:) = logical(checkpointSummary.cleanedOnSuccess);
end
checkpointSummaryTable = table(runNameList, runKeyList, runDirList, numTaskList, ...
  numDoneList, isCompleteList, cleanedOnSuccess, 'VariableNames', ...
  {'runName', 'runKey', 'runDir', 'numTask', 'numDone', 'isComplete', 'cleanedOnSuccess'});
end

function value = localPadColumn(value, numRow, fillValue)
%LOCALPADCOLUMN Pad or trim a column vector for lightweight summary tables.

value = reshape(value, [], 1);
if numel(value) >= numRow
  value = value(1:numRow);
  return;
end
value(end + 1:numRow, 1) = fillValue;
end

function localPlotCleanEstimateVsCrbCurves(curveTable)
%LOCALPLOTCLEANESTIMATEVSCRBCURVES Plot static/MF estimates and CRB on log axes.

if isempty(curveTable) || ~istable(curveTable) || height(curveTable) == 0
  return;
end
methodNameList = unique(string(curveTable.displayName), 'stable');
methodLabelList = strrep(methodNameList, "-", " ");
markerList = {'o', 's', '^', 'v', 'd', 'p', '>', '<', 'x', '+', '*', 'h'};
filterLabel = string(curveTable.filterName(1));
try
  figure('Name', 'Clean static/MF estimate vs CRB');
  tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

  axAngle = nexttile;
  hold(axAngle, 'on');
  for iMethod = 1:numel(methodNameList)
    methodMask = curveTable.displayName == methodNameList(iMethod);
    methodTable = sortrows(curveTable(methodMask, :), 'snrDb');
    if height(methodTable) == 0
      continue;
    end
    marker = markerList{1 + mod(iMethod - 1, numel(markerList))};
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.angleRmseDeg), ...
      ['-' marker], 'DisplayName', methodLabelList(iMethod) + " est");
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.angleCrbSphericalMedianDeg), ...
      ['--' marker], 'DisplayName', methodLabelList(iMethod) + " CRB");
  end
  xlabel('SNR (dB)');
  ylabel('Spherical angle RMSE and CRB std (deg)');
  title(sprintf('Angle, %s', char(filterLabel)));
  set(axAngle, 'YScale', 'log');
  grid(axAngle, 'on');
  enableLegendToggle(axAngle);

  axFdRef = nexttile;
  hold(axFdRef, 'on');
  for iMethod = 1:numel(methodNameList)
    methodMask = curveTable.displayName == methodNameList(iMethod);
    methodTable = sortrows(curveTable(methodMask, :), 'snrDb');
    if height(methodTable) == 0
      continue;
    end
    marker = markerList{1 + mod(iMethod - 1, numel(markerList))};
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.fdRefRmseHz), ...
      ['-' marker], 'DisplayName', methodLabelList(iMethod) + " est");
    semilogy(methodTable.snrDb, localPositiveForLog(methodTable.fdRefCrbMedianHz), ...
      ['--' marker], 'DisplayName', methodLabelList(iMethod) + " CRB");
  end
  xlabel('SNR (dB)');
  ylabel('fdRef RMSE and CRB std (Hz)');
  title(sprintf('Reference Doppler, %s', char(filterLabel)));
  set(axFdRef, 'YScale', 'log');
  grid(axFdRef, 'on');
  enableLegendToggle(axFdRef);
catch ME
  warning('scanMfMsMleCrbCleanBoundConsistency:PlotFailed', ...
    'Clean estimate-vs-CRB plot failed: %s', ME.message);
end
end

function value = localPositiveForLog(value)
%LOCALPOSITIVEFORLOG Replace nonpositive values with NaN for log-scale plots.

value(value <= 0 | ~isfinite(value)) = NaN;
end

function mask = localFilterMask(caseTable, filterName)
%LOCALFILTERMASK Return row filter for clean comparison.

filterName = string(filterName);
switch filterName
  case "raw"
    mask = true(height(caseTable), 1);
  case "health"
    mask = caseTable.healthResolved;
  case "joint-trim"
    mask = caseTable.jointTrimKeep;
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:UnknownFilter', 'Unsupported filter "%s".', char(filterName));
end
end

function value = localMsErrorOverMatchedSsCrb(caseTable, msMask, ssRaw)
%LOCALMSERROROVERMATCHEDSSCRB Normalize MS angle error with same-seed SS-MF CRB.

msIdx = find(msMask(:));
ratioList = NaN(numel(msIdx), 1);
for i = 1:numel(msIdx)
  idx = msIdx(i);
  ssIdx = find(ssRaw & caseTable.taskSeed == caseTable.taskSeed(idx), 1, 'first');
  if isempty(ssIdx)
    continue;
  end
  ssCrb = caseTable.crbSphericalApproxStdDeg(ssIdx);
  if isfinite(ssCrb) && ssCrb > 0
    ratioList(i) = (caseTable.sphericalAngleErrDeg(idx) ./ ssCrb).^2;
  end
end
value = mean(ratioList(isfinite(ratioList)), 'omitnan');
end

function auditTable = localBuildCleanSearchRangeAuditTable(caseTable)
%LOCALBUILDCLEANSEARCHRANGEAUDITTABLE Check clean local boxes against CRB scale.

if isempty(caseTable) || height(caseTable) == 0
  auditTable = table();
  return;
end
methodList = unique(caseTable.displayName, 'stable');
rowTemplate = struct('displayName', "", 'snrDb', NaN, 'numRepeat', NaN, ...
  'minDoaLatHalfOverCrb', NaN, 'minDoaLonHalfOverCrb', NaN, ...
  'minFdRefHalfOverCrb', NaN, 'doaBelowCrb', false, 'fdRefBelowCrb', false, ...
  'doaTruthBoxViolationRate', NaN);
snrList = unique(caseTable.snrDb, 'stable');
rowList = repmat(rowTemplate, numel(methodList) * numel(snrList), 1);
idx = 0;
for iSnr = 1:numel(snrList)
  snrDb = snrList(iSnr);
  for iMethod = 1:numel(methodList)
    idx = idx + 1;
    mask = caseTable.displayName == methodList(iMethod) & caseTable.snrDb == snrDb;
    row = rowTemplate;
    row.displayName = methodList(iMethod);
    row.snrDb = snrDb;
    row.numRepeat = nnz(mask);
    row.minDoaLatHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLatDeg(mask), caseTable.crbLatStdDeg(mask)), [], 'omitnan');
    row.minDoaLonHalfOverCrb = min(localSafeVectorRatio(caseTable.searchHalfWidthLonDeg(mask), caseTable.crbLonStdDeg(mask)), [], 'omitnan');
    row.minFdRefHalfOverCrb = min(localSafeVectorRatio(caseTable.fdRangeHalfWidthHz(mask), caseTable.fdRefCrbStdHz(mask)), [], 'omitnan');
    row.doaBelowCrb = row.minDoaLatHalfOverCrb < 1 || row.minDoaLonHalfOverCrb < 1;
    row.fdRefBelowCrb = row.minFdRefHalfOverCrb < 1;
    if ismember('doaTruthBoxViolation', caseTable.Properties.VariableNames)
      row.doaTruthBoxViolationRate = mean(double(caseTable.doaTruthBoxViolation(mask)), 'omitnan');
    end
    rowList(idx) = row;
  end
end
auditTable = struct2table(rowList(:));
end

function outlierTable = localBuildCleanOutlierTable(caseTable)
%LOCALBUILDCLEANOUTLIERTABLE Keep MS-MF unknown-rate rows rejected by clean trim.

if isempty(caseTable) || height(caseTable) == 0
  outlierTable = table();
  return;
end
mask = startsWith(caseTable.displayName, "MS-MF-") & endsWith(caseTable.displayName, "-U") & ~caseTable.jointTrimKeep;
columnList = {'displayName', 'snrDb', 'taskSeed', 'sphericalAngleErrDeg', ...
  'angleErrOverSphericalCrb', 'fdRefAbsErrHz', 'fdRefErrOverCrb', ...
  'fdRateAbsErrHzPerSec', 'iterations', 'firstOrderOpt', 'objectiveImprove', ...
  'fdRefBoundaryHit', 'fdRateBoundaryHit', 'healthResolved', 'jointTrimKeep', 'rejectReason'};
columnList = columnList(ismember(columnList, caseTable.Properties.VariableNames));
outlierTable = caseTable(mask, columnList);
if height(outlierTable) > 0
  outlierTable = sortrows(outlierTable, {'angleErrOverSphericalCrb'}, {'descend'});
end
end

function localPrintReplayConfig(config)
%LOCALPRINTREPLAYCONFIG Print clean scan axes.

fprintf('  %-32s : %s\n', 'SNR list (dB)', localFormatRow(config.snrDbList));
fprintf('  %-32s : %d\n', 'frame count', config.numFrame);
fprintf('  %-32s : %s\n', 'user LLA', localFormatRow(config.contextUsrLla));
fprintf('  %-32s : %s\n', 'selected satellites', localFormatRow(config.contextSelectedSatIdxGlobal));
fprintf('  %-32s : %.0f\n', 'reference satellite', config.contextRefSatIdxGlobal);
fprintf('  %-32s : %s\n', 'method list', strjoin(cellstr(config.methodNameList), ', '));
methodRangePreview = config.methodRangeTable(:, {'methodName', ...
  'doaCrbScale', 'fdRefRangeMode', 'fdRefHalfToothFraction', ...
  'fdRefCrbScale', 'fdRateHalfWidthHzPerSec'});
fprintf('  %-32s : %d rows\n', 'method range table', height(methodRangePreview));
fprintf('  %-32s : %.2f tooth\n', 'fdRef CRB-scale tooth cap', ...
  localGetFieldOrDefault(config, 'defaultFdRefMaxHalfToothFraction', NaN));
if height(methodRangePreview) <= 12
  disp(methodRangePreview);
else
  disp(methodRangePreview(1:12, :));
  fprintf('  ... %d more method range rows\n', height(methodRangePreview) - 12);
end
fprintf('  %-32s : %s\n', 'estimation chain', char(config.initMode));
fprintf('  %-32s : %s / %s\n', 'MF objective / route', ...
  char(config.dynamicObjectiveMode), char(config.dynamicRouteMode));
fprintf('  %-32s : %s / %s\n', 'MF signal / steering', ...
  char(config.cleanMfSignalTimeModel), char(config.cleanMfSteeringMode));
fprintf('  %-32s : %.2f / %.2f CRB\n', 'trim angle/fdRef cap', config.trimAngleNormMax, config.trimFdRefNormMax);
fprintf('  %-32s : %d / %d\n', 'task count / outer parfor', config.numTask, logical(config.useParfor));
end

function localPrintCompactTable(sectionTitle, dataTable, edgeCount)
%LOCALPRINTCOMPACTTABLE Print a compact scan table through common report helper.

if nargin < 3
  edgeCount = 12;
end
printMfReportTableSection(char(sectionTitle), dataTable, ...
  'MaxPreviewRow', edgeCount, 'SectionPrinter', @printMfScanSection);
end

function lineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build short notification metrics.

lineList = strings(0, 1);
compareTable = scanData.cleanSsMsCompareTable;
knownUnknownTable = localGetFieldOrDefault(scanData, 'cleanKnownUnknownCompareTable', table());
if istable(knownUnknownTable) && height(knownUnknownTable) > 0
  trimKuRow = knownUnknownTable(knownUnknownTable.filterName == "joint-trim" & knownUnknownTable.satMode == "multi" & knownUnknownTable.phaseMode == "continuous", :);
  if ~isempty(trimKuRow)
    row = trimKuRow(1, :);
    lineList(end + 1, 1) = sprintf('MS K/U: dMSE/CRB=%.3g, U gain=%.3g, class=%s', ...
      row.unknownMinusKnownAngleMseOverCrb, row.unknownRmseGainOverKnown, char(row.comparisonClass)); %#ok<AGROW>
  end
end
compareTable = scanData.cleanSsMsCompareTable;
if istable(compareTable) && height(compareTable) > 0
  trimRow = compareTable(compareTable.filterName == "joint-trim" & compareTable.phaseMode == "continuous" & compareTable.fdRateMode == "unknown", :);
  if ~isempty(trimRow)
    row = trimRow(1, :);
    lineList(end + 1, 1) = sprintf('joint-trim: ms/ownCRB=%.3g, ms/ssCRB=%.3g, msAbsGain=%.3g, class=%s', ...
      row.msAngleMseOverOwnCrb, row.msAngleMseOverSsCrb, row.msAbsGainOverSs, char(row.targetClass)); %#ok<AGROW>
  end
end
lineList(end + 1, 1) = sprintf('elapsed %.1f min', scanData.elapsedSec / 60); %#ok<AGROW>
end

function checkpointOpt = localBuildCheckpointOpt(scanName, config, useParfor)
%LOCALBUILDCHECKPOINTOPT Build common checkpoint options with a stable run key.

meta = struct('snrDbList', reshape(config.snrDbList, 1, []), ...
  'seedList', reshape(config.seedList, 1, []), ...
  'methodRangeTable', config.methodRangeTable, ...
  'contextUsrLla', reshape(config.contextUsrLla, 1, []), ...
  'contextSelectedSatIdxGlobal', reshape(config.contextSelectedSatIdxGlobal, 1, []), ...
  'contextRefSatIdxGlobal', config.contextRefSatIdxGlobal, ...
  'contextUtcRef', char(string(config.contextUtcRef)), ...
  'contextTleFileName', char(config.contextTleFileName), ...
  'dynamicObjectiveMode', char(config.dynamicObjectiveMode), ...
  'dynamicRouteMode', char(config.dynamicRouteMode), ...
  'cleanMfSignalTimeModel', char(config.cleanMfSignalTimeModel), ...
  'cleanMfSteeringMode', char(config.cleanMfSteeringMode), ...
  'routeTag', char(config.routeTag));
checkpointOpt = buildMfReplayCheckpointOpt(scanName, config, struct( ...
  'runKey', config.runKey, 'meta', meta, 'useParfor', logical(useParfor)));
end

function runKey = localBuildRunKey(config)
%LOCALBUILDRUNKEY Build path-safe run key.

runKey = sprintf('frames%d_snr%s_seed%dto%d_rep%d_%s', config.numFrame, ...
  localCompactNumberList(config.snrDbList), config.seedList(1), config.seedList(end), ...
  numel(config.seedList), localCompactStringHash(localBuildRunSignature(config)));
runKey = string(runKey);
end

function signature = localBuildRunSignature(config)
%LOCALBUILDRUNSIGNATURE Build deterministic signature for checkpoint compatibility.

signature = strjoin(string({ ...
  sprintf('frames%d', config.numFrame), ...
  sprintf('snr%s', mat2str(reshape(double(config.snrDbList), 1, []), 8)), ...
  sprintf('seed%s', mat2str(reshape(double(config.seedList), 1, []), 8)), ...
  sprintf('methods%s', char(strjoin(reshape(string(config.methodNameList), 1, []), ','))), ...
  sprintf('range%s', char(localMethodRangeSignature(config.methodRangeTable))), ...
  sprintf('usr%s', mat2str(reshape(double(config.contextUsrLla), 1, []), 8)), ...
  sprintf('sat%s', mat2str(reshape(double(config.contextSelectedSatIdxGlobal), 1, []), 8)), ...
  sprintf('ref%.0f', config.contextRefSatIdxGlobal), ...
  sprintf('utc%s', char(string(config.contextUtcRef))), ...
  sprintf('tle%s', char(config.contextTleFileName)), ...
  sprintf('init%s', char(string(config.initMode))), ...
  sprintf('objective%s', char(string(config.dynamicObjectiveMode))), ...
  sprintf('route%s', char(string(config.dynamicRouteMode))), ...
  sprintf('signal%s', char(string(config.cleanMfSignalTimeModel))), ...
  sprintf('steering%s', char(string(config.cleanMfSteeringMode))), ...
  sprintf('trimAngle%.8g', config.trimAngleNormMax), ...
  sprintf('trimFdRef%.8g', config.trimFdRefNormMax), ...
  char(config.routeTag) ...
  }), '|');
end

function signature = localMethodRangeSignature(methodRangeTable)
%LOCALMETHODRANGESIGNATURE Build compact text for method-range hashing.

if isempty(methodRangeTable) || ~istable(methodRangeTable)
  signature = "empty";
  return;
end
rowText = strings(height(methodRangeTable), 1);
for iRow = 1:height(methodRangeTable)
  rowText(iRow) = sprintf('%s:doa%.8g:fdmode%s:fd%.8g:fdcrb%.8g:rate%.8g', ...
    char(string(methodRangeTable.methodName(iRow))), ...
    methodRangeTable.doaCrbScale(iRow), char(string(methodRangeTable.fdRefRangeMode(iRow))), ...
    methodRangeTable.fdRefHalfToothFraction(iRow), methodRangeTable.fdRefCrbScale(iRow), ...
    methodRangeTable.fdRateHalfWidthHzPerSec(iRow));
end
signature = strjoin(rowText, ';');
end

function hashText = localCompactStringHash(textValue)
%LOCALCOMPACTSTRINGHASH Return a short deterministic hash for path-safe names.

byteValues = double(char(textValue));
hashValue = 5381;
modValue = 2147483647;
for iByte = 1:numel(byteValues)
  hashValue = mod(hashValue * 33 + byteValues(iByte), modValue);
end
hashText = sprintf('%08x', uint32(hashValue));
end

function textValue = localCompactNumberList(valueList)
%LOCALCOMPACTNUMBERLIST Format numeric axes for a checkpoint run key.

valueList = reshape(double(valueList), 1, []);
textValue = strjoin(compose('%.6g', valueList), '_');
textValue = replace(textValue, '-', 'm');
textValue = replace(textValue, '.', 'p');
end

function useParfor = localCanUseParfor(numTask)
%LOCALCANUSEPARFOR Check whether the outer scan task loop can use parfor.

useParfor = false;
if numTask < 2 || ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  useParfor = ~isempty(ver('parallel'));
catch
  useParfor = false;
end
end

function flowOpt = localBuildFlowOpt(parallelOpt)
%LOCALBUILDFLOWOPT Build static and dynamic options shared by clean methods.

flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicPolishEnableWhenMulti', false, ...
  'staticMsHalfWidth', [1e-8; 1e-8]));
flowOpt = applyDynamicTransitionFlowDefaults(flowOpt);
end

function parallelOpt = localBuildSerialInnerParallelOpt()
%LOCALBUILDSERIALINNERPARALLELOPT Disable estimator-inner parallelism.

parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(numFrame);
if numFrame < 2
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidFrameCount', ...
    'Frame count must be at least two for multi-frame diagnostics.');
end
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function crbBundle = localBuildCrbBundleQuiet(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar, config)
%LOCALBUILDCRBBUNDLEQUIET Build CP/IP CRB while suppressing expected warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbOptOverride = struct();
crbOptOverride.steeringMode = char(config.cleanMfSteeringMode);
crbOptOverride.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar, crbOptOverride);
crbBundle = localAddIndependentPhaseCrb(crbBundle, periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar, config);
end

function crbBundle = localAddIndependentPhaseCrb(crbBundle, periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar, config)
%LOCALADDINDEPENDENTPHASECRB Add local IP-K/IP-U CRB entries for this scan.

truth = periodicFixture.truth;
pathGainRef = ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
noiseVarRef = noiseVar * ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
pathGainMs = ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);
noiseVarMs = noiseVar * ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);

ipKnownOpt = struct();
ipKnownOpt.doaType = 'latlon';
ipKnownOpt.phaseMode = 'independent';
ipKnownOpt.fdRateMode = 'known';
ipKnownOpt.steeringMode = char(config.cleanMfSteeringMode);
ipKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
[crbBundle.crbMfRefIpKnown, crbBundle.auxCrbMfRefIpKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, ipKnownOpt);
[crbBundle.crbMfMsIpKnown, crbBundle.auxCrbMfMsIpKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, ipKnownOpt);

ipUnknownOpt = ipKnownOpt;
ipUnknownOpt.fdRateMode = 'unknown';
[crbBundle.crbMfRefIpUnknown, crbBundle.auxCrbMfRefIpUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, ipUnknownOpt);
[crbBundle.crbMfMsIpUnknown, crbBundle.auxCrbMfMsIpUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, ipUnknownOpt);
end

function [crb, aux] = localTryBuildDynamicCrb(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt)
%LOCALTRYBUILDDYNAMICCRB Build one dynamic CRB and return NaNs on failure.

try
  [crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = getDoaDopplerFieldOrDefault(crbOpt, 'fdRateMode', "unknown");
  aux.phaseMode = getDoaDopplerFieldOrDefault(crbOpt, 'phaseMode', "independent");
  aux.errorMessage = string(ME.message);
end
end

function [crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName)
%LOCALRESOLVECRBFORMETHOD Return the CRB block matching one method.

methodName = string(methodName);
methodName = localBaseDynamicMethodName(methodName);
switch methodName
  case "SS-SF-DoA"
    crbFull = crbBundle.crbSfDoaRef;
  case "MS-SF-DoA"
    crbFull = crbBundle.crbSfDoaMs;
  case "SS-SF-Static"
    crbFull = crbBundle.crbSfRef;
  case "MS-SF-Static"
    crbFull = crbBundle.crbSfMs;
  case "SS-MF-Static"
    crbFull = crbBundle.crbMfRefStatic;
  case "MS-MF-Static"
    crbFull = crbBundle.crbMfMsStatic;
  case "SS-MF-CP-K"
    crbFull = crbBundle.crbMfRefKnown;
  case "SS-MF-CP-U"
    crbFull = crbBundle.crbMfRefUnknown;
  case "MS-MF-CP-K"
    crbFull = crbBundle.crbMfMsKnown;
  case "MS-MF-CP-U"
    crbFull = crbBundle.crbMfMsUnknown;
  case "SS-MF-IP-K"
    crbFull = crbBundle.crbMfRefIpKnown;
  case "SS-MF-IP-U"
    crbFull = crbBundle.crbMfRefIpUnknown;
  case "MS-MF-IP-K"
    crbFull = crbBundle.crbMfMsIpKnown;
  case "MS-MF-IP-U"
    crbFull = crbBundle.crbMfMsIpUnknown;
  otherwise
    error('scanMfMsMleCrbCleanBoundConsistency:UnknownCrbMethod', 'Unsupported CRB method "%s".', char(methodName));
end
crbDoa = real(crbFull(1:2, 1:2));
end

function [crbDoa, crbFull, crbSource] = localResolveFiniteCrbForMethod(crbBundle, methodName)
%LOCALRESOLVEFINITECRBFORMETHOD Resolve a finite positive DoA CRB for clean metrics.

methodName = string(methodName);
[crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
crbSource = methodName;
if localIsValidPositiveDoaCrb(crbDoa)
  return;
end

fallbackName = localFallbackDoaCrbMethodName(methodName);
if fallbackName == ""
  return;
end
[crbDoaFallback, crbFullFallback] = localResolveCrbForMethod(crbBundle, fallbackName);
if localIsValidPositiveDoaCrb(crbDoaFallback)
  crbDoa = crbDoaFallback;
  crbFull = crbFullFallback;
  crbSource = fallbackName;
end
end

function fallbackName = localFallbackDoaCrbMethodName(methodName)
%LOCALFALLBACKDOACRBMETHODNAME Choose a matched static CRB when DoA-only CRB degenerates.

methodName = string(methodName);
switch methodName
  case "SS-SF-DoA"
    fallbackName = "SS-SF-Static";
  case "MS-SF-DoA"
    fallbackName = "MS-SF-Static";
  otherwise
    fallbackName = "";
end
end

function tf = localIsValidPositiveDoaCrb(crbDoa)
%LOCALISVALIDPOSITIVEDOACRB Return true for finite positive 2-D DoA CRB blocks.

tf = isnumeric(crbDoa) && all(size(crbDoa) >= [2 2]);
if ~tf
  return;
end
crbBlock = real(crbDoa(1:2, 1:2));
diagVal = diag(crbBlock);
tf = all(isfinite(crbBlock(:))) && all(diagVal > 0) && ...
  all(eig(0.5 * (crbBlock + crbBlock.')) > 0);
end

function method = localApplyDoaCrbScaledWidth(method, crbBundle)
%LOCALAPPLYDOACRBSCALEDWIDTH Set the local DoA box from the method CRB.

method.doaHalfWidthDeg = localResolveDoaCrbScaledHalfWidth( ...
  method.baseDisplayName, method.doaHalfWidthDeg, crbBundle);
end

function doaHalfWidthDeg = localResolveDoaCrbScaledHalfWidth(methodName, doaCrbScale, crbBundle)
%LOCALRESOLVEDOACRBSCALEDHALFWIDTH Build a method-specific CRB-scaled DoA box.

doaCrbScale = double(doaCrbScale);
if isempty(doaCrbScale) || ~isfinite(doaCrbScale(1)) || doaCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidDoaCrbScale', ...
    'The clean DoA CRB scale must be a positive finite scalar.');
end
[crbDoa, ~] = localResolveFiniteCrbForMethod(crbBundle, methodName);
doaStdDeg = sqrt(max(real(diag(crbDoa)), 0));
if numel(doaStdDeg) < 2 || any(~isfinite(doaStdDeg(1:2)))
  doaHalfWidthDeg = NaN(2, 1);
  return;
end
doaHalfWidthDeg = doaCrbScale(1) * reshape(doaStdDeg(1:2), [], 1);
if ~localIsValidPositiveDoaHalfWidth(doaHalfWidthDeg)
  doaHalfWidthDeg = NaN(2, 1);
end
end

function fdRangeUse = localBuildFdRefCrbScaledRange(truthFdRefHz, crbBundle, ...
  methodName, fdRefCrbScale, toothStepHz, maxHalfToothFraction)
%LOCALBUILDFDREFCRBSCALEDRANGE Build a pure CRB-scaled fdRef box.

fdRefCrbScale = double(fdRefCrbScale);
if isempty(fdRefCrbScale) || ~isfinite(fdRefCrbScale(1)) || fdRefCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidFdRefCrbScale', ...
    'The fdRef CRB scale must be a positive finite scalar.');
end
[~, crbFull] = localResolveFiniteFdRefCrbForMethod(crbBundle, methodName);
fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
if ~(isfinite(fdRefCrbStdHz) && fdRefCrbStdHz > 0)
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidFdRefCrbStd', ...
    'Cannot build a CRB-scaled fdRef box for method "%s".', char(methodName));
end
fdHalfWidthHz = fdRefCrbScale(1) * fdRefCrbStdHz;
maxHalfToothFraction = double(maxHalfToothFraction);
if isfinite(maxHalfToothFraction) && maxHalfToothFraction > 0 && ...
    isfinite(toothStepHz) && toothStepHz > 0
  fdHalfWidthHz = min(fdHalfWidthHz, maxHalfToothFraction * toothStepHz);
end
fdRangeUse = truthFdRefHz + fdHalfWidthHz * [-1, 1];
end

function [crbDoa, crbFull, crbSource] = localResolveFiniteFdRefCrbForMethod(crbBundle, methodName)
%LOCALRESOLVEFINITEFDREFCRBFORMETHOD Resolve a finite fdRef CRB for range building.

methodName = string(methodName);
[crbDoa, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
crbSource = methodName;
if localIsValidFdRefCrb(crbFull)
  return;
end
fallbackName = localFallbackDoaCrbMethodName(methodName);
if fallbackName == ""
  return;
end
[crbDoaFallback, crbFullFallback] = localResolveCrbForMethod(crbBundle, fallbackName);
if localIsValidFdRefCrb(crbFullFallback)
  crbDoa = crbDoaFallback;
  crbFull = crbFullFallback;
  crbSource = fallbackName;
end
end

function tf = localIsValidFdRefCrb(crbFull)
%LOCALISVALIDFDREFCRB Return true when a CRB matrix has a finite fdRef entry.

tf = isnumeric(crbFull) && all(size(crbFull) >= [3 3]) && ...
  isfinite(real(crbFull(3, 3))) && real(crbFull(3, 3)) > 0;
end

function fdRangeUse = localApplyFdRefCrbScaledFloor(fdRangeUse, truthFdRefHz, crbBundle, methodNameList, fdRefCrbScale)
%LOCALAPPLYFDREFCRBSCALEDFLOOR Ensure the fdRef box is not a one-sigma CRB truncation.

fdRefCrbScale = double(fdRefCrbScale);
if isempty(fdRefCrbScale) || ~isfinite(fdRefCrbScale(1)) || fdRefCrbScale(1) <= 0
  error('scanMfMsMleCrbCleanBoundConsistency:InvalidFdRefCrbScale', ...
    'The clean fdRef CRB floor scale must be a positive finite scalar.');
end
fdHalfWidthHz = 0.5 * abs(fdRangeUse(2) - fdRangeUse(1));
fdRefCrbFloorHz = NaN;
for iMethod = 1:numel(methodNameList)
  methodName = string(methodNameList(iMethod));
  try
    [~, crbFull] = localResolveCrbForMethod(crbBundle, methodName);
    fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
    if isfinite(fdRefCrbStdHz)
      fdRefCrbFloorHz = max([fdRefCrbFloorHz; fdRefCrbScale(1) * fdRefCrbStdHz], [], 'omitnan');
    end
  catch
  end
end
if isfinite(fdRefCrbFloorHz) && fdHalfWidthHz < fdRefCrbFloorHz
  fdRangeUse = truthFdRefHz + fdRefCrbFloorHz * [-1, 1];
end
end

function toothStepHz = localResolveToothStepHz(fixture)
%LOCALRESOLVETOOTHSTEPHZ Resolve Doppler tooth spacing from frame offsets.

timeOffsetSec = [];
if isfield(fixture, 'sceneSeq') && isfield(fixture.sceneSeq, 'timeOffsetSec')
  timeOffsetSec = fixture.sceneSeq.timeOffsetSec(:);
end
if numel(timeOffsetSec) < 2
  toothStepHz = 750;
  return;
end
dt = median(diff(sort(timeOffsetSec)));
toothStepHz = 1 / dt;
end

function doaHalfWidth = localNormalizeDoaHalfWidth(doaHalfWidth)
%LOCALNORMALIZEDOAHALFWIDTH Normalize scalar or vector DoA half-width input.

doaHalfWidth = reshape(double(doaHalfWidth), [], 1);
if isscalar(doaHalfWidth)
  doaHalfWidth = repmat(doaHalfWidth, 2, 1);
end
if numel(doaHalfWidth) < 2
  doaHalfWidth = NaN(2, 1);
else
  doaHalfWidth = doaHalfWidth(1:2);
end
end

function tf = localIsValidPositiveDoaHalfWidth(doaHalfWidth)
%LOCALISVALIDPOSITIVEDOAHALFWIDTH Return true for finite positive lat/lon half-width.

doaHalfWidth = localNormalizeDoaHalfWidth(doaHalfWidth);
tf = numel(doaHalfWidth) >= 2 && all(isfinite(doaHalfWidth(1:2))) && all(doaHalfWidth(1:2) > 0);
end

function tf = localIsNearRangeBoundary(value, rangeValue, boundaryTolFraction)
%LOCALISNEARRANGEBOUNDARY Check if one value is too close to a frequency bound.

tf = false;
rangeValue = reshape(double(rangeValue), 1, []);
if numel(rangeValue) ~= 2 || any(~isfinite(rangeValue)) || ~isfinite(value)
  return;
end
rangeWidth = abs(rangeValue(2) - rangeValue(1));
if rangeWidth <= 0
  return;
end
boundaryTol = max(boundaryTolFraction * rangeWidth, eps(rangeWidth));
tf = abs(value - min(rangeValue)) <= boundaryTol || abs(value - max(rangeValue)) <= boundaryTol;
end

function row = localEmptyCaseRow()
%LOCALEMPTYCASEROW Return a typed empty row for clean case metrics.

row = struct('displayName', "", 'satMode', "", 'frameMode', "", 'phaseMode', "", ...
  'fdRateMode', "", 'mfInitMode', "", 'dynamicObjectiveMode', "", ...
  'dynamicRouteMode', "", 'snrDb', NaN, 'taskSeed', NaN, ...
  'doaCrbScale', NaN, 'fdRefRangeMode', "", ...
  'fdRefHalfToothFraction', NaN, 'fdRefCrbScale', NaN, ...
  'fdRateHalfWidthHzPerSec', NaN, 'numFrame', NaN, ...
  'frameIntvlSec', NaN, 'windowSec', NaN, 'toothStepHz', NaN, ...
  'truthLatDeg', NaN, 'truthLonDeg', NaN, 'initLatDeg', NaN, 'initLonDeg', NaN, ...
  'finalLatDeg', NaN, 'finalLonDeg', NaN, 'latErrDeg', NaN, 'lonErrDeg', NaN, ...
  'initAngleErrDeg', NaN, 'finalAngleErrDeg', NaN, 'sphericalAngleErrDeg', NaN, ...
  'finalMinusInitAngleDeg', NaN, 'searchHalfWidthLatDeg', NaN, ...
  'searchHalfWidthLonDeg', NaN, 'seedLatDeg', NaN, 'seedLonDeg', NaN, ...
  'boundCenterLatDeg', NaN, 'boundCenterLonDeg', NaN, ...
  'doaLatBoxViolation', false, 'doaLonBoxViolation', false, ...
  'doaTruthBoxViolation', false, 'crbLatStdDeg', NaN, 'crbLonStdDeg', NaN, ...
  'crbTraceStdDeg', NaN, 'crbSphericalApproxStdDeg', NaN, ...
  'crbAngleMetricVarRad2', NaN, 'crbMetricTraceRatio', NaN, ...
  'angleErrOverSphericalCrb', NaN, 'angleErrOverTraceCrb', NaN, ...
  'fdRefTrueHz', NaN, 'initFdRefHz', NaN, 'finalFdRefHz', NaN, ...
  'fdRefErrHz', NaN, 'fdRefAbsErrHz', NaN, 'fdRefCrbStdHz', NaN, ...
  'fdRefErrOverCrb', NaN, 'fdRefInitMoveHz', NaN, 'fdRangeLowHz', NaN, ...
  'fdRangeHighHz', NaN, 'fdRangeHalfWidthHz', NaN, 'fdRefBoundaryHit', false, ...
  'fdRateTrueHzPerSec', NaN, 'initFdRateHzPerSec', NaN, ...
  'finalFdRateHzPerSec', NaN, 'fdRateErrHzPerSec', NaN, ...
  'fdRateAbsErrHzPerSec', NaN, 'fdRateInitMoveHzPerSec', NaN, ...
  'fdRateRangeLowHzPerSec', NaN, 'fdRateRangeHighHzPerSec', NaN, ...
  'fdRateBoundaryHit', false, 'fdLineRmseHz', NaN, 'solverResolved', false, ...
  'exitflag', NaN, 'iterations', NaN, 'funcCount', NaN, 'stepSize', NaN, ...
  'firstOrderOpt', NaN, 'constrViolation', NaN, 'algorithm', "", ...
  'usedScaledSolve', false, 'fallbackTriggered', false, 'solveVariant', "", ...
  'candidateCount', NaN, 'candidateVariantText', "", 'initObj', NaN, ...
  'finalObj', NaN, 'objectiveImprove', NaN, 'finalResidualNorm', NaN, ...
  'refCoherence', NaN, 'nonRefCoherenceFloor', NaN, 'runTimeMs', NaN, ...
  'isDynamicMethod', false, 'solverHealthy', false, 'freqHealthy', false, ...
  'rateHealthy', false, 'boundaryHealthy', false, 'notNoSolve', false, ...
  'healthResolved', false, 'trimKeep', false, 'jointTrimKeep', false, ...
  'fdRefNormTail', false, 'rejectReason', "", 'warningFlag', false);
end

function progressTracker = localCreateProgressTracker(totalCount)
%LOCALCREATEPROGRESSTRACKER Create a best-effort scan progress tracker.

progressTracker = struct('enabled', false, 'queue', []);
if totalCount <= 0
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', totalCount);
  progressTracker.enabled = true;
  try
    progressTracker.queue = parallel.pool.DataQueue;
    afterEach(progressTracker.queue, @(~) progressbar('advance'));
  catch
    progressTracker.queue = [];
  end
catch
  progressTracker.enabled = false;
  progressTracker.queue = [];
end
end

function localAdvanceProgress(progressTracker)
%LOCALADVANCEPROGRESS Advance progress when progressbar is available.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('advance');
catch
end
end

function localCloseProgressTracker(progressTracker)
%LOCALCLOSEPROGRESSTRACKER Close a best-effort scan progress tracker.

if ~(isstruct(progressTracker) && isfield(progressTracker, 'enabled') && progressTracker.enabled)
  return;
end
try
  progressbar('end');
catch
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a struct field or return a fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Divide only when denominator is finite and positive.

if isfinite(numerator) && isfinite(denominator) && denominator > 0
  ratio = numerator ./ denominator;
else
  ratio = NaN;
end
end

function ratio = localSafeVectorRatio(numerator, denominator)
%LOCALSAFEVECTORRATIO Element-wise safe ratio for table columns.

numerator = reshape(double(numerator), [], 1);
denominator = reshape(double(denominator), [], 1);
ratio = NaN(size(numerator));
mask = isfinite(numerator) & isfinite(denominator) & denominator > 0;
ratio(mask) = numerator(mask) ./ denominator(mask);
end

function value = localMse(x)
%LOCALMSE Mean-square of finite entries.

x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = mean(x(:).^2);
end
end

function value = localRmse(x)
%LOCALRMSE Root-mean-square of finite entries.

value = sqrt(localMse(x));
end

function value = localNormalizedMse(err, crbStd)
%LOCALNORMALIZEDMSE Mean squared CRB-normalized error over finite entries.

err = reshape(double(err), [], 1);
crbStd = reshape(double(crbStd), [], 1);
mask = isfinite(err) & isfinite(crbStd) & crbStd > 0;
if ~any(mask)
  value = NaN;
else
  value = mean((err(mask) ./ crbStd(mask)).^2);
end
end

function key = localMethodKey(methodName)
%LOCALMETHODKEY Convert a display name into a struct field key.

key = matlab.lang.makeValidName(char(regexprep(string(methodName), '[^A-Za-z0-9]', '_')));
end

function out = ternary(conditionValue, trueValue, falseValue)
%TERNARY Return one of two values based on a scalar condition.

if conditionValue
  out = trueValue;
else
  out = falseValue;
end
end

function textValue = localFormatRow(valueList)
%LOCALFORMATROW Format numeric/string row values.

if isstring(valueList) || ischar(valueList)
  textValue = strjoin(cellstr(reshape(string(valueList), 1, [])), ', ');
  return;
end
valueList = reshape(double(valueList), 1, []);
if isempty(valueList)
  textValue = '';
else
  textValue = strjoin(compose('%.6g', valueList), ', ');
end
end

function contextSummary = localBuildContextSummary(context)
%LOCALBUILDCONTEXTSUMMARY Keep a lightweight context summary.

contextSummary = struct();
contextSummary.carrierFreq = localGetFieldOrDefault(context, 'carrierFreq', NaN);
contextSummary.wavelen = localGetFieldOrDefault(context, 'wavelen', NaN);
contextSummary.usrLla = localGetFieldOrDefault(context, 'usrLla', NaN);
contextSummary.utcRef = localGetFieldOrDefault(context, 'utcRef', NaT);
contextSummary.selectedSatIdxGlobal = localGetFieldOrDefault(context, 'selectedSatIdxGlobal', NaN);
contextSummary.refSatIdxGlobal = localGetFieldOrDefault(context, 'refSatIdxGlobal', NaN);
contextSummary.refSatIdxLocal = localGetFieldOrDefault(context, 'refSatIdxLocal', NaN);
contextSummary.nonRefSatIdxGlobal = localGetFieldOrDefault(context, 'nonRefSatIdxGlobal', NaN);
contextSummary.otherSatIdxGlobal = localGetFieldOrDefault(context, 'otherSatIdxGlobal', NaN);
if isfield(context, 'waveInfo')
  contextSummary.sampleRate = localGetFieldOrDefault(context.waveInfo, 'sampleRate', NaN);
end
end

function repeatCell = localStripRepeatCell(repeatCell)
%LOCALSTRIPREPEATCELL Drop heavyweight task data before snapshot saving.

for iRepeat = 1:numel(repeatCell)
  if isempty(repeatCell{iRepeat})
    continue;
  end
  if isfield(repeatCell{iRepeat}, 'truth')
    repeatCell{iRepeat} = rmfield(repeatCell{iRepeat}, 'truth');
  end
end
end
