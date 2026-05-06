%SCANMFREGIMEMAPBYWINDOW Scan physical DoA-static / Doppler-dynamic regimes.
% Dev scan only. It maps window length, frame interval, and frame count to
% paper-level regime metrics without running Monte Carlo estimators.
% The script leaves scanData in the workspace; saveSnapshot=true saves only
% lightweight scanData through saveExpSnapshot.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfRegimeMapByWindow";
saveSnapshot = false;
notifyTelegramEnable = true;
checkpointEnable = false;
frameCountList = (1:100).';
frameIntvlMsList = [0.5; 0.75; 1; 1e3/750; 1.5; 2; 2.5; 3];
curveWindowMsList = linspace(0, 130, 600).';
representativeFrameCountList = [5; 8; 10; 15; 20];
representativeFrameIntvlMsList = [1; 1e3/750; 2];
carrierFreqHz = 11.7e9;
minSlantRangeM = 550e3;
transverseVelocityMps = 13.6e3;
primaryDoaSlowTolDeg = 0.1;
doaSlowTolDegList = [0.02; 0.05; 0.1; 0.2];
fdDynamicTolHz = 50;
staticDopplerPhaseTolRad = 1;
firstOrderFdResidualTolHz = 1;
firstOrderPhaseResidualTolRad = 0.1;

frameCountList = reshape(double(frameCountList), [], 1);
frameIntvlSecList = reshape(double(frameIntvlMsList), [], 1) * 1e-3;
curveWindowSecList = reshape(double(curveWindowMsList), [], 1) * 1e-3;
representativeFrameCountList = reshape(double(representativeFrameCountList), [], 1);
representativeFrameIntvlSecList = reshape(double(representativeFrameIntvlMsList), [], 1) * 1e-3;
doaSlowTolDegList = unique(reshape(double(doaSlowTolDegList), [], 1), 'stable');
if ~any(abs(doaSlowTolDegList - primaryDoaSlowTolDeg) < 1e-12)
  doaSlowTolDegList = sort([doaSlowTolDegList; primaryDoaSlowTolDeg]);
end

scanConfig = struct();
scanConfig.scanName = string(scanName);
scanConfig.saveSnapshot = logical(saveSnapshot);
scanConfig.notifyTelegramEnable = logical(notifyTelegramEnable);
scanConfig.checkpointEnable = logical(checkpointEnable);
scanConfig.frameCountList = frameCountList;
scanConfig.frameIntvlSecList = frameIntvlSecList;
scanConfig.curveWindowSecList = curveWindowSecList;
scanConfig.representativeFrameCountList = representativeFrameCountList;
scanConfig.representativeFrameIntvlSecList = representativeFrameIntvlSecList;
scanConfig.carrierFreqHz = carrierFreqHz;
scanConfig.minSlantRangeM = minSlantRangeM;
scanConfig.transverseVelocityMps = transverseVelocityMps;
scanConfig.doaSlowTolDeg = primaryDoaSlowTolDeg;
scanConfig.primaryDoaSlowTolDeg = primaryDoaSlowTolDeg;
scanConfig.doaSlowTolDegList = doaSlowTolDegList;
scanConfig.fdDynamicTolHz = fdDynamicTolHz;
scanConfig.staticDopplerPhaseTolRad = staticDopplerPhaseTolRad;
scanConfig.firstOrderFdResidualTolHz = firstOrderFdResidualTolHz;
scanConfig.firstOrderPhaseResidualTolRad = firstOrderPhaseResidualTolRad;

scanConfig.numTask = numel(frameCountList) * numel(frameIntvlSecList);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
checkpointDir = "";
scanJustRan = false;
printMfScanHeader(char(scanName), scanConfig, checkpointDir);

%% Build context and scan tasks
% This regime scan is deterministic and does not build estimator context or
% flow options. It only uses paper-level physical order metrics.
numConfig = scanConfig.numTask;
rowCell = cell(numConfig, 1);
configIdx = 0;

%% Run scan batch
try
  fprintf('[%s] Build %d deterministic window-regime rows.\n', char(datetime('now', 'Format', 'HH:mm:ss')), numConfig);
  for iP = 1:numel(scanConfig.frameCountList)
    numFrame = scanConfig.frameCountList(iP);
    for iTf = 1:numel(scanConfig.frameIntvlSecList)
      frameIntvlSec = scanConfig.frameIntvlSecList(iTf);
      configIdx = configIdx + 1;
      rowCell{configIdx} = localBuildRegimeRow(scanConfig, numFrame, frameIntvlSec);
    end
  end

  scanTable = struct2table([rowCell{:}].');
  aggregateTable = scanTable;
  curveTable = localBuildCurveTable(scanConfig);
  targetSummaryTable = localBuildTargetSummaryTable(scanTable, scanConfig);
  representativeTable = localBuildRepresentativeTable(scanTable, scanConfig);
  modelBoundaryTable = localBuildModelBoundaryTable(curveTable, scanConfig);
  doaToleranceSummaryTable = localBuildDoaToleranceSummaryTable(curveTable, scanConfig);

  %% Data storage
  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.config = scanConfig;
  scanData.scanTable = scanTable;
  scanData.aggregateTable = aggregateTable;
  scanData.curveTable = curveTable;
  scanData.targetSummaryTable = targetSummaryTable;
  scanData.doaToleranceSummaryTable = doaToleranceSummaryTable;
  scanData.representativeTable = representativeTable;
  scanData.modelBoundaryTable = modelBoundaryTable;
  scanData.plotData = localBuildPlotData(scanTable, aggregateTable, curveTable, doaToleranceSummaryTable);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.elapsedSec = toc(runTic);
  scanData = finalizeMfScanResult(scanData, checkpointDir);

  if scanConfig.saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end
  scanJustRan = true;
catch ME
  fprintf('Scan failed while building the deterministic regime table.\n');
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "FAILED", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Summary output and plotting
if ~exist('scanData', 'var') || ~isstruct(scanData)
  error('Scan data is missing. Run the scan batch sections or load a snapshot containing scanData.');
end

scanConfig = localNormalizeConfig(scanData.config);
scanData.config = scanConfig;
scanTable = scanData.scanTable;
aggregateTable = scanData.aggregateTable;
curveTable = scanData.curveTable;
scanData.modelBoundaryTable = localBuildModelBoundaryTable(curveTable, scanConfig);
if ~isfield(scanData, 'doaToleranceSummaryTable')
  scanData.doaToleranceSummaryTable = localBuildDoaToleranceSummaryTable(curveTable, scanConfig);
end

printMfScanSection('Regime summary by frame interval', scanData.targetSummaryTable);
printMfScanSection('DoA tolerance sensitivity', scanData.doaToleranceSummaryTable);
printMfScanSection('Model boundary summary', scanData.modelBoundaryTable);
printMfScanSection('Representative rows', scanData.representativeTable);

scanData.plotData = localPlotScan(scanTable, aggregateTable, curveTable, scanConfig, scanData.doaToleranceSummaryTable);

if exist('scanJustRan', 'var') && scanJustRan
  notifyMfScanStatus(struct( ...
    'scanName', scanName, ...
    'statusText', "DONE", ...
    'config', scanConfig, ...
    'snapshotFile', scanData.snapshotFile, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', localBuildTelegramMetricLines(scanData), ...
    'commentLineList', [ ...
      "Regime window map completed."; ...
      "Detailed tables remain in scanData and scan results docs."]));
end

%% Local helpers

function row = localBuildRegimeRow(scanConfig, numFrame, frameIntvlSec)
%LOCALBUILDREGIMEROW Compute one deterministic window-regime row.
windowSec = numFrame * frameIntvlSec;
timeOffsetSpanSec = (numFrame - 1) * frameIntvlSec;
metrics = localComputeWindowMetrics(scanConfig, windowSec);
row = metrics;
row.numFrame = numFrame;
row.frameIntvlSec = frameIntvlSec;
row.frameIntvlMs = frameIntvlSec * 1e3;
row.windowSec = windowSec;
row.windowMs = windowSec * 1e3;
row.timeOffsetSpanSec = timeOffsetSpanSec;
row.timeOffsetSpanMs = timeOffsetSpanSec * 1e3;
end


function metrics = localComputeWindowMetrics(scanConfig, windowSec)
%LOCALCOMPUTEWINDOWMETRICS Compute deterministic physical order metrics.
lightSpeedMps = 299792458;
fc = scanConfig.carrierFreqHz;
rho = scanConfig.minSlantRangeM;
vPerp = scanConfig.transverseVelocityMps;
fdRateApproxHzPerSec = fc / lightSpeedMps * vPerp^2 / rho;

rangeEndM = sqrt(rho^2 + (vPerp * windowSec)^2);
doaDriftRadExact = atan2(vPerp * windowSec, rho);
fdDriftHzExact = fc / lightSpeedMps * vPerp^2 * windowSec / rangeEndM;
phaseStaticDopplerErrRad = 2 * pi * fc / lightSpeedMps * (rangeEndM - rho);
phaseFirstOrderModelRad = pi * fdRateApproxHzPerSec * windowSec^2;
fdFirstOrderResidualHz = fdDriftHzExact - fdRateApproxHzPerSec * windowSec;
phaseFirstOrderResidualRad = phaseStaticDopplerErrRad - phaseFirstOrderModelRad;

doaDriftRadApprox = vPerp / rho * windowSec;
fdDriftHzApprox = fdRateApproxHzPerSec * windowSec;
quadPhaseRadApprox = phaseFirstOrderModelRad;

metrics = struct();
metrics.doaDriftRadExact = doaDriftRadExact;
metrics.doaDriftDegExact = doaDriftRadExact * 180 / pi;
metrics.doaDriftRadApprox = doaDriftRadApprox;
metrics.doaDriftDegApprox = doaDriftRadApprox * 180 / pi;
metrics.doaApproxResidualDeg = abs(metrics.doaDriftDegExact - metrics.doaDriftDegApprox);
metrics.fdRateApproxHzPerSec = fdRateApproxHzPerSec;
metrics.fdDriftHzExact = fdDriftHzExact;
metrics.fdDriftHzApprox = fdDriftHzApprox;
metrics.fdFirstOrderResidualHz = fdFirstOrderResidualHz;
metrics.phaseStaticDopplerErrRad = phaseStaticDopplerErrRad;
metrics.quadPhaseRadApprox = quadPhaseRadApprox;
metrics.quadPhaseCycleApprox = quadPhaseRadApprox / (2 * pi);
metrics.phaseFirstOrderResidualRad = phaseFirstOrderResidualRad;
metrics.phaseFirstOrderResidualCycle = phaseFirstOrderResidualRad / (2 * pi);
metrics.fdDriftPerDoaRadHz = fdRateApproxHzPerSec / (vPerp / rho);
metrics.quadPhasePerDoaRad2 = pi * fdRateApproxHzPerSec / (vPerp / rho)^2;
metrics.isDoaStaticValid = metrics.doaDriftDegExact <= scanConfig.primaryDoaSlowTolDeg;
metrics.isStaticDopplerInvalid = abs(metrics.fdDriftHzExact) >= scanConfig.fdDynamicTolHz || ...
  abs(metrics.phaseStaticDopplerErrRad) >= scanConfig.staticDopplerPhaseTolRad;
metrics.isFirstOrderDopplerValid = abs(metrics.fdFirstOrderResidualHz) <= scanConfig.firstOrderFdResidualTolHz && ...
  abs(metrics.phaseFirstOrderResidualRad) <= scanConfig.firstOrderPhaseResidualTolRad;
metrics.isTargetRegime = metrics.isDoaStaticValid && metrics.isStaticDopplerInvalid && ...
  metrics.isFirstOrderDopplerValid;
end


function curveTable = localBuildCurveTable(scanConfig)
%LOCALBUILDCURVETABLE Build continuous window-length curves for plotting.
windowSecList = scanConfig.curveWindowSecList;
curveCell = cell(numel(windowSecList), 1);
for iWin = 1:numel(windowSecList)
  windowSec = windowSecList(iWin);
  row = localComputeWindowMetrics(scanConfig, windowSec);
  row.windowSec = windowSec;
  row.windowMs = windowSec * 1e3;
  curveCell{iWin} = row;
end
curveTable = struct2table([curveCell{:}].');
end


function summaryTable = localBuildTargetSummaryTable(scanTable, scanConfig)
%LOCALBUILDTARGETSUMMARYTABLE Summarize target-regime coverage by frame interval.
frameIntvlSecList = scanConfig.frameIntvlSecList;
rowCell = cell(numel(frameIntvlSecList), 1);
for iTf = 1:numel(frameIntvlSecList)
  tf = frameIntvlSecList(iTf);
  idx = abs(scanTable.frameIntvlSec - tf) < 1e-12;
  subTable = scanTable(idx, :);
  targetP = subTable.numFrame(subTable.isTargetRegime);
  row = struct();
  row.frameIntvlMs = tf * 1e3;
  row.numConfig = height(subTable);
  row.numTarget = sum(subTable.isTargetRegime);
  row.minTargetFrame = localMinOrNan(targetP);
  row.maxTargetFrame = localMaxOrNan(targetP);
  row.maxTargetWindowMs = localMaxOrNan(subTable.windowMs(subTable.isTargetRegime));
  row.firstDoaStaticFailFrame = localFirstFrame(subTable, ~subTable.isDoaStaticValid);
  row.firstStaticDopplerInvalidFrame = localFirstFrame(subTable, subTable.isStaticDopplerInvalid);
  row.firstFirstOrderFailFrame = localFirstFrame(subTable, ~subTable.isFirstOrderDopplerValid);
  rowCell{iTf} = row;
end
summaryTable = struct2table([rowCell{:}].');
end


function representativeTable = localBuildRepresentativeTable(scanTable, scanConfig)
%LOCALBUILDREPRESENTATIVETABLE Extract compact rows for paper-scale windows.
keepIdx = false(height(scanTable), 1);
for iP = 1:numel(scanConfig.representativeFrameCountList)
  for iTf = 1:numel(scanConfig.representativeFrameIntvlSecList)
    keepIdx = keepIdx | (scanTable.numFrame == scanConfig.representativeFrameCountList(iP) & ...
      abs(scanTable.frameIntvlSec - scanConfig.representativeFrameIntvlSecList(iTf)) < 1e-12);
  end
end
fieldList = {'numFrame', 'frameIntvlMs', 'windowMs', 'timeOffsetSpanMs', ...
  'doaDriftDegExact', 'doaApproxResidualDeg', 'fdDriftHzExact', ...
  'fdRateApproxHzPerSec', 'phaseStaticDopplerErrRad', ...
  'fdFirstOrderResidualHz', 'phaseFirstOrderResidualRad', ...
  'isDoaStaticValid', 'isStaticDopplerInvalid', ...
  'isFirstOrderDopplerValid', 'isTargetRegime'};
representativeTable = scanTable(keepIdx, fieldList);
end


function boundaryTable = localBuildModelBoundaryTable(curveTable, scanConfig)
%LOCALBUILDMODELBOUNDARYTABLE Report continuous pass/fail boundary windows.
conditionCell = {curveTable.isDoaStaticValid; curveTable.isStaticDopplerInvalid; ...
  curveTable.isFirstOrderDopplerValid; curveTable.isTargetRegime};
metricList = ["DoA static validity"; "static Doppler invalid"; ...
  "first-order Doppler validity"; "target regime"];
ruleList = [string(sprintf('DoA drift <= %.4f deg', scanConfig.primaryDoaSlowTolDeg)); ...
  string(sprintf('|fd drift| >= %.2f Hz or static phase >= %.2f rad', ...
  scanConfig.fdDynamicTolHz, scanConfig.staticDopplerPhaseTolRad)); ...
  string(sprintf('|fd residual| <= %.2f Hz and |phase residual| <= %.2f rad', ...
  scanConfig.firstOrderFdResidualTolHz, scanConfig.firstOrderPhaseResidualTolRad)); ...
  "all three conditions"];
rowCell = cell(numel(conditionCell), 1);
for iRow = 1:numel(conditionCell)
  passIdx = conditionCell{iRow};
  row = struct();
  row.metric = metricList(iRow);
  row.passRule = ruleList(iRow);
  row.firstPassWindowMs = localFirstWindow(curveTable, passIdx);
  row.lastPassWindowMs = localLastWindow(curveTable, passIdx);
  row.firstFailAfterPassWindowMs = localFirstFailAfterPassWindow(curveTable, passIdx);
  rowCell{iRow} = row;
end
boundaryTable = struct2table([rowCell{:}].');
end


function summaryTable = localBuildDoaToleranceSummaryTable(curveTable, scanConfig)
%LOCALBUILDDOATOLERANCESUMMARYTABLE Summarize regime width for DoA tolerances.
tolList = scanConfig.doaSlowTolDegList;
staticDopplerLowerMs = localFirstWindow(curveTable, curveTable.isStaticDopplerInvalid);
firstOrderLastMs = localLastWindow(curveTable, curveTable.isFirstOrderDopplerValid);
rowCell = cell(numel(tolList), 1);
for iTol = 1:numel(tolList)
  tolDeg = tolList(iTol);
  doaValidIdx = curveTable.doaDriftDegExact <= tolDeg;
  targetIdx = doaValidIdx & curveTable.isStaticDopplerInvalid & curveTable.isFirstOrderDopplerValid;
  row = struct();
  row.doaTolDeg = tolDeg;
  row.doaLastValidWindowMs = localLastWindow(curveTable, doaValidIdx);
  row.doaFirstInvalidWindowMs = localFirstWindowAfterLastPass(curveTable, doaValidIdx);
  row.targetFirstWindowMs = localFirstWindow(curveTable, targetIdx);
  row.targetLastWindowMs = localLastWindow(curveTable, targetIdx);
  row.staticDopplerLowerWindowMs = staticDopplerLowerMs;
  row.firstOrderLastValidWindowMs = firstOrderLastMs;
  row.isPrimaryDoaTol = abs(tolDeg - scanConfig.primaryDoaSlowTolDeg) < 1e-12;
  rowCell{iTol} = row;
end
summaryTable = struct2table([rowCell{:}].');
end


function plotData = localBuildPlotData(scanTable, aggregateTable, curveTable, doaToleranceSummaryTable)
%LOCALBUILDPLOTDATA Store lightweight data needed to redraw figures.
plotData = struct();
plotData.scanTable = scanTable;
plotData.aggregateTable = aggregateTable;
plotData.curveTable = curveTable;
plotData.doaToleranceSummaryTable = doaToleranceSummaryTable;
end


function plotData = localPlotScan(scanTable, aggregateTable, curveTable, scanConfig, doaToleranceSummaryTable)
%LOCALPLOTSCAN Draw deterministic regime figures without saving image files.
plotData = localBuildPlotData(scanTable, aggregateTable, curveTable, doaToleranceSummaryTable);

figure('Name', 'Fixed-DoA validity curve');
plot(curveTable.windowMs, curveTable.doaDriftDegExact, 'LineWidth', 1.5); hold on;
for iTol = 1:numel(scanConfig.doaSlowTolDegList)
  tolDeg = scanConfig.doaSlowTolDegList(iTol);
  if abs(tolDeg - scanConfig.primaryDoaSlowTolDeg) < 1e-12
    yline(tolDeg, '--', sprintf('primary %.3g deg', tolDeg));
  else
    yline(tolDeg, ':', sprintf('%.3g deg', tolDeg));
  end
end
scatter(scanTable.windowMs, scanTable.doaDriftDegExact, 20, 'filled');
grid on; xlabel('paper window P T_f (ms)'); ylabel('DoA drift (deg)');
title('Fixed-DoA approximation check');

figure('Name', 'Static-Doppler mismatch curve');
plot(curveTable.windowMs, curveTable.fdDriftHzExact, 'LineWidth', 1.5); hold on;
yline(scanConfig.fdDynamicTolHz, '--', 'fd-drift threshold');
scatter(scanTable.windowMs, scanTable.fdDriftHzExact, 20, 'filled');
grid on; xlabel('paper window P T_f (ms)'); ylabel('Doppler drift (Hz)');
title('Constant-Doppler mismatch');

figure('Name', 'Doppler phase model residuals');
semilogy(curveTable.windowMs, max(abs(curveTable.phaseStaticDopplerErrRad), eps), 'LineWidth', 1.5); hold on;
semilogy(curveTable.windowMs, max(abs(curveTable.phaseFirstOrderResidualRad), eps), 'LineWidth', 1.5);
yline(scanConfig.staticDopplerPhaseTolRad, '--', 'static-Doppler phase threshold');
yline(scanConfig.firstOrderPhaseResidualTolRad, ':', 'first-order residual threshold');
grid on; xlabel('paper window P T_f (ms)'); ylabel('phase error (rad)');
legend({'constant-Doppler phase error', 'first-order Doppler residual'}, 'Location', 'best');
title('Static Doppler vs first-order Doppler model');

figure('Name', 'DoA-Doppler relation curve');
plot(curveTable.doaDriftDegExact, curveTable.fdDriftHzExact, 'LineWidth', 1.5); hold on;
scatter(scanTable.doaDriftDegExact, scanTable.fdDriftHzExact, 20, scanTable.frameIntvlMs, 'filled');
grid on; xlabel('DoA drift (deg)'); ylabel('Doppler drift (Hz)');
title('DoA drift and Doppler drift relation');
cb = colorbar; ylabel(cb, 'frame interval (ms)');

drawnow;
end


function val = localMinOrNan(x)
%LOCALMINORNAN Return min(x) or NaN for empty input.
if isempty(x)
  val = NaN;
else
  val = min(x);
end
end


function val = localMaxOrNan(x)
%LOCALMAXORNAN Return max(x) or NaN for empty input.
if isempty(x)
  val = NaN;
else
  val = max(x);
end
end


function frameVal = localFirstFrame(subTable, idx)
%LOCALFIRSTFRAME Return the first frame count satisfying a logical mask.
if ~any(idx)
  frameVal = NaN;
  return;
end
frameVal = min(subTable.numFrame(idx));
end


function windowMs = localFirstWindow(curveTable, idx)
%LOCALFIRSTWINDOW Return first window value satisfying a logical mask.
if ~any(idx)
  windowMs = NaN;
  return;
end
windowMs = min(curveTable.windowMs(idx));
end


function windowMs = localLastWindow(curveTable, idx)
%LOCALLASTWINDOW Return last window value satisfying a logical mask.
if ~any(idx)
  windowMs = NaN;
  return;
end
windowMs = max(curveTable.windowMs(idx));
end


function windowMs = localFirstFailAfterPassWindow(curveTable, passIdx)
%LOCALFIRSTFAILAFTERPASSWINDOW Return first failing window after any pass.
if ~any(passIdx)
  windowMs = localFirstWindow(curveTable, ~passIdx);
  return;
end
lastPassMs = localLastWindow(curveTable, passIdx);
failAfterIdx = ~passIdx & curveTable.windowMs > lastPassMs;
windowMs = localFirstWindow(curveTable, failAfterIdx);
end


function windowMs = localFirstWindowAfterLastPass(curveTable, passIdx)
%LOCALFIRSTWINDOWAFTERLASTPASS Return first sampled window after the pass region.
if ~any(passIdx)
  windowMs = localFirstWindow(curveTable, ~passIdx);
  return;
end
lastPassMs = localLastWindow(curveTable, passIdx);
windowMs = localFirstWindow(curveTable, ~passIdx & curveTable.windowMs > lastPassMs);
end


function scanConfig = localNormalizeConfig(scanConfig)
%LOCALNORMALIZECONFIG Add derived fields for old snapshots and current runs.
if ~isfield(scanConfig, 'primaryDoaSlowTolDeg')
  scanConfig.primaryDoaSlowTolDeg = scanConfig.doaSlowTolDeg;
end
if ~isfield(scanConfig, 'doaSlowTolDegList')
  scanConfig.doaSlowTolDegList = scanConfig.primaryDoaSlowTolDeg;
end
scanConfig.doaSlowTolDegList = reshape(double(scanConfig.doaSlowTolDegList), [], 1);
end


function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build compact scan-specific notification lines.

metricLineList = strings(0, 1);
config = scanData.config;
metricLineList(end + 1, 1) = sprintf('• Frame grid: <code>%d x %d</code>', ...
  numel(config.frameCountList), numel(config.frameIntvlSecList));
if isfield(scanData, 'targetSummaryTable') && ~isempty(scanData.targetSummaryTable)
  targetCount = sum(scanData.targetSummaryTable.numTarget);
  metricLineList(end + 1, 1) = sprintf('• Target-regime rows: <code>%d</code>', targetCount);
end
if isfield(config, 'primaryDoaSlowTolDeg')
  primaryDoaSlowTolDeg = config.primaryDoaSlowTolDeg;
elseif isfield(config, 'doaSlowTolDeg')
  primaryDoaSlowTolDeg = config.doaSlowTolDeg;
else
  primaryDoaSlowTolDeg = NaN;
end
if isfinite(primaryDoaSlowTolDeg)
  metricLineList(end + 1, 1) = sprintf('• Primary DoA tol: <code>%.4g deg</code>', ...
    primaryDoaSlowTolDeg);
end
end
