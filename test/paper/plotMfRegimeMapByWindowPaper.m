%PLOTMFREGIMEMAPBYWINDOWPAPER Redraw the paper-facing DoA/Doppler regime figure.
% Load the existing scanMfRegimeMapByWindow snapshot and combine the three
% main regime diagnostics into one compact figure for slides or papers.

clear; close all; clc;

%% Paper figure configuration

scriptName = "plotMfRegimeMapByWindowPaper";
repoRoot = localFindRepoRoot();
addpath(fullfile(repoRoot, 'utils', 'io'));
addpath(fullfile(repoRoot, 'test', 'common', 'util'));

preferredSnapshotFile = fullfile(repoRoot, 'test', 'data', 'cache', 'scan', ...
  'scanMfRegimeMapByWindow_20260427-213020.mat');
snapshotFile = resolveTestSnapshotFile(preferredSnapshotFile, ...
  'scanMfRegimeMapByWindow_*.mat');
mainFrameIntvlSec = 1 / 750;
highlightFrameCountList = [10; 20];
layoutMode = "horizontal"; % "horizontal" for slides, "vertical" for notes.
exportFigure = false;
exportDir = fullfile(repoRoot, 'tmp', char(scriptName));
exportBaseName = "mfRegimeMapByWindowPaper";

%% Load scan snapshot

snapshotData = loadExpSnapshot(snapshotFile, 'none', struct('varNames', {{'scanData'}}));
if ~isfield(snapshotData, 'scanData') || ~isstruct(snapshotData.scanData)
  error('plotMfRegimeMapByWindowPaper:MissingScanData', ...
    'The snapshot does not contain a valid scanData struct.');
end
scanData = localNormalizeScanData(snapshotData.scanData);
scanConfig = scanData.config;
scanData.scanTable = localEnsureWindowColumns(scanData.scanTable);
scanData.curveTable = localEnsureWindowColumns(scanData.curveTable);
curveTable = scanData.curveTable;
scanTable = scanData.scanTable;

highlightTable = localBuildHighlightTable(scanTable, mainFrameIntvlSec, highlightFrameCountList);
fprintf('\nPaper regime highlight rows:\n');
disp(highlightTable(:, {'numFrame', 'frameIntvlMs', 'windowMs', ...
  'doaDriftDegExact', 'fdDriftHzExact', 'phaseStaticDopplerErrRad', ...
  'phaseFirstOrderResidualRad'}));

summaryTable = localBuildDocumentSummaryTable(highlightTable);
fprintf('\nCompact document summary rows:\n');
disp(summaryTable);

%% Redraw compact paper figure

fig = localCreateFigure(layoutMode);
ax = gobjects(3, 1);

ax(1) = nexttile;
plot(curveTable.windowMs, curveTable.doaDriftDegExact, 'LineWidth', 1.5); hold on;
yline(scanConfig.primaryDoaSlowTolDeg, '--', ...
  sprintf('DoA tolerance %.3g deg', scanConfig.primaryDoaSlowTolDeg));
localMarkHighlightWindows(ax(1), highlightTable, 'doaDriftDegExact');
grid on;
xlabel('window length P T_f (ms)');
ylabel('DoA drift (deg)');
title('(a) Fixed-DoA approximation');

ax(2) = nexttile;
plot(curveTable.windowMs, curveTable.fdDriftHzExact, 'LineWidth', 1.5); hold on;
yline(scanConfig.fdDynamicTolHz, '--', ...
  sprintf('Doppler drift %.3g Hz', scanConfig.fdDynamicTolHz));
localMarkHighlightWindows(ax(2), highlightTable, 'fdDriftHzExact');
grid on;
xlabel('window length P T_f (ms)');
ylabel('Doppler drift (Hz)');
title('(b) Constant-Doppler mismatch');

ax(3) = nexttile;
staticPhase = max(abs(curveTable.phaseStaticDopplerErrRad), eps);
firstOrderResidual = max(abs(curveTable.phaseFirstOrderResidualRad), eps);
semilogy(curveTable.windowMs, staticPhase, 'LineWidth', 1.5); hold on;
semilogy(curveTable.windowMs, firstOrderResidual, 'LineWidth', 1.5);
yline(scanConfig.staticDopplerPhaseTolRad, '--', ...
  sprintf('static threshold %.3g rad', scanConfig.staticDopplerPhaseTolRad));
yline(scanConfig.firstOrderPhaseResidualTolRad, ':', ...
  sprintf('first-order threshold %.3g rad', scanConfig.firstOrderPhaseResidualTolRad));
localMarkPhaseHighlight(ax(3), highlightTable);
grid on;
xlabel('window length P T_f (ms)');
ylabel('phase error (rad)');
title('(c) Static Doppler vs first-order model');
legend({'constant-Doppler phase error', 'first-order Doppler residual'}, ...
  'Location', 'best');

localApplySharedXLimit(ax, curveTable.windowMs, highlightTable.windowMs);
sgtitle('DoA quasi-static / Doppler-dynamic multi-frame regime');

if exportFigure
  if ~isfolder(exportDir)
    mkdir(exportDir);
  end
  pngFile = fullfile(exportDir, exportBaseName + ".png");
  pdfFile = fullfile(exportDir, exportBaseName + ".pdf");
  exportgraphics(fig, pngFile, 'Resolution', 300);
  exportgraphics(fig, pdfFile, 'ContentType', 'vector');
  fprintf('Exported paper regime figure:\n  %s\n  %s\n', pngFile, pdfFile);
end

%% Local helpers


function tableIn = localEnsureWindowColumns(tableIn)
%LOCALENSUREWINDOWCOLUMNS Add window/time columns used by the paper script.

if ~ismember('windowMs', tableIn.Properties.VariableNames) && ...
    ismember('windowSec', tableIn.Properties.VariableNames)
  tableIn.windowMs = tableIn.windowSec * 1e3;
end
if ~ismember('frameIntvlMs', tableIn.Properties.VariableNames) && ...
    ismember('frameIntvlSec', tableIn.Properties.VariableNames)
  tableIn.frameIntvlMs = tableIn.frameIntvlSec * 1e3;
end
end

function summaryTable = localBuildDocumentSummaryTable(highlightTable)
%LOCALBUILDDOCUMENTSUMMARYTABLE Build a compact table suitable for notes.

summaryTable = table();
summaryTable.numFrame = highlightTable.numFrame;
summaryTable.windowMs = highlightTable.windowMs;
summaryTable.doaDriftDeg = highlightTable.doaDriftDegExact;
summaryTable.fdDriftHz = highlightTable.fdDriftHzExact;
summaryTable.staticPhaseErrRad = highlightTable.phaseStaticDopplerErrRad;
summaryTable.firstOrderResidualRad = highlightTable.phaseFirstOrderResidualRad;
summaryTable.isTargetRegime = highlightTable.isTargetRegime;
end

function scanData = localNormalizeScanData(scanData)
%LOCALNORMALIZESCANDATA Add fields needed by old snapshots and plotting.

if ~isfield(scanData, 'config') || ~isstruct(scanData.config)
  error('plotMfRegimeMapByWindowPaper:InvalidScanData', ...
    'scanData.config is missing or invalid.');
end
if ~isfield(scanData, 'curveTable') || ~istable(scanData.curveTable)
  error('plotMfRegimeMapByWindowPaper:InvalidScanData', ...
    'scanData.curveTable is missing or invalid.');
end
if ~isfield(scanData, 'scanTable') || ~istable(scanData.scanTable)
  error('plotMfRegimeMapByWindowPaper:InvalidScanData', ...
    'scanData.scanTable is missing or invalid.');
end

scanConfig = scanData.config;
scanData.scanTable = localEnsureWindowColumns(scanData.scanTable);
scanData.curveTable = localEnsureWindowColumns(scanData.curveTable);
if ~isfield(scanConfig, 'primaryDoaSlowTolDeg')
  if isfield(scanConfig, 'doaSlowTolDeg')
    scanConfig.primaryDoaSlowTolDeg = scanConfig.doaSlowTolDeg;
  else
    scanConfig.primaryDoaSlowTolDeg = 0.1;
  end
end
if ~isfield(scanConfig, 'fdDynamicTolHz')
  scanConfig.fdDynamicTolHz = 50;
end
if ~isfield(scanConfig, 'staticDopplerPhaseTolRad')
  scanConfig.staticDopplerPhaseTolRad = 1;
end
if ~isfield(scanConfig, 'firstOrderPhaseResidualTolRad')
  scanConfig.firstOrderPhaseResidualTolRad = 0.1;
end
scanData.config = scanConfig;
end

function highlightTable = localBuildHighlightTable(scanTable, mainFrameIntvlSec, frameCountList)
%LOCALBUILDHIGHLIGHTTABLE Pick representative rows for the main paper window.

frameCountList = reshape(double(frameCountList), [], 1);
rowIdx = false(height(scanTable), 1);
for iFrame = 1:numel(frameCountList)
  rowIdx = rowIdx | (scanTable.numFrame == frameCountList(iFrame) & ...
    abs(scanTable.frameIntvlSec - mainFrameIntvlSec) < 1e-12);
end
highlightTable = scanTable(rowIdx, :);
if height(highlightTable) ~= numel(frameCountList)
  error('plotMfRegimeMapByWindowPaper:MissingHighlightRows', ...
    'Could not find all requested P values at T_f = %.12g s.', mainFrameIntvlSec);
end
[~, orderIdx] = sort(highlightTable.numFrame);
highlightTable = highlightTable(orderIdx, :);
end

function fig = localCreateFigure(layoutMode)
%LOCALCREATEFIGURE Create the paper figure and tiled layout.

fig = figure('Name', 'Paper DoA-Doppler regime map', 'Color', 'w');
if layoutMode == "vertical"
  fig.Position(3:4) = [780, 900];
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
  fig.Position(3:4) = [1320, 430];
  tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end
end

function localMarkHighlightWindows(ax, highlightTable, yFieldName)
%LOCALMARKHIGHLIGHTWINDOWS Mark P=10/P=20 style representative windows.

axes(ax); %#ok<LAXES>
yValue = highlightTable.(yFieldName);
scatter(highlightTable.windowMs, yValue, 52, 'filled', 'MarkerEdgeColor', 'k');
for iRow = 1:height(highlightTable)
  xline(highlightTable.windowMs(iRow), ':', sprintf('P=%d', highlightTable.numFrame(iRow)), ...
    'LabelVerticalAlignment', 'bottom');
end
end

function localMarkPhaseHighlight(ax, highlightTable)
%LOCALMARKPHASEHIGHLIGHT Mark representative windows on both phase curves.

axes(ax); %#ok<LAXES>
scatter(highlightTable.windowMs, max(abs(highlightTable.phaseStaticDopplerErrRad), eps), ...
  52, 'filled', 'MarkerEdgeColor', 'k');
scatter(highlightTable.windowMs, max(abs(highlightTable.phaseFirstOrderResidualRad), eps), ...
  52, 'd', 'filled', 'MarkerEdgeColor', 'k');
for iRow = 1:height(highlightTable)
  xline(highlightTable.windowMs(iRow), ':', sprintf('P=%d', highlightTable.numFrame(iRow)), ...
    'LabelVerticalAlignment', 'bottom');
end
end

function localApplySharedXLimit(ax, windowMsList, highlightWindowMsList)
%LOCALAPPLYSHAREDXLIMIT Focus the x-axis on the useful paper-window range.

maxWindowMs = max([60; highlightWindowMsList(:) * 1.5]);
maxWindowMs = min(maxWindowMs, max(windowMsList));
for iAx = 1:numel(ax)
  xlim(ax(iAx), [0, maxWindowMs]);
end
end

function repoRoot = localFindRepoRoot()
%LOCALFINDREPOROOT Find repository root from this paper script location.

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
repoRoot = thisDir;
while ~isempty(repoRoot) && ~isfile(fullfile(repoRoot, 'AGENTS.md'))
  parentDir = fileparts(repoRoot);
  if strcmp(parentDir, repoRoot)
    break;
  end
  repoRoot = parentDir;
end
if isempty(repoRoot) || ~isfile(fullfile(repoRoot, 'AGENTS.md'))
  error('plotMfRegimeMapByWindowPaper:RepoRootNotFound', ...
    'Cannot find repository root from paper script path.');
end
end
