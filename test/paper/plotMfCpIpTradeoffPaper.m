%PLOTMFCPIPTRADEOFFPAPER Redraw the paper-facing controlled CP/IP trade-off figure.
% Load the existing scanMfCpIpInToothPerfMap snapshot and build a compact
% English figure for the continuous-phase versus independent-phase discussion.

clear; close all; clc;

%% Paper figure configuration

scriptName = "plotMfCpIpTradeoffPaper";
repoRoot = localFindRepoRoot();
addpath(fullfile(repoRoot, 'utils', 'io'));
addpath(fullfile(repoRoot, 'test', 'common', 'util'));

preferredSnapshotFile = fullfile(repoRoot, 'test', 'data', 'cache', 'scan', ...
  'scanMfCpIpInToothPerfMap_20260428-213853.mat');
snapshotFile = resolveTestSnapshotFile(preferredSnapshotFile, ...
  'scanMfCpIpInToothPerfMap_*.mat');
mainFrameIntvlSec = 1 / 750;
mainSnrDb = 10;
highlightFrameCountList = [10; 20];
layoutMode = "horizontal"; % "horizontal" for slides, "vertical" for notes.
exportFigure = false;
exportDir = fullfile(repoRoot, 'tmp', char(scriptName));
exportBaseName = "mfCpIpTradeoffPaper";

%% Load scan snapshot

snapshotData = loadExpSnapshot(snapshotFile, 'none', struct('varNames', {{'scanData'}}));
if ~isfield(snapshotData, 'scanData') || ~isstruct(snapshotData.scanData)
  error('plotMfCpIpTradeoffPaper:MissingScanData', ...
    'The snapshot does not contain a valid scanData struct.');
end
scanData = localNormalizeScanData(snapshotData.scanData, mainFrameIntvlSec, mainSnrDb);
mainSliceTable = scanData.mainSliceTable;
highlightTable = localBuildHighlightTable(mainSliceTable, highlightFrameCountList);

fprintf('\nPaper controlled CP/IP trade-off highlight rows:\n');
if ~isempty(highlightTable)
  disp(highlightTable(:, {'displayName', 'snrDb', 'numFrame', 'frameIntvlMs', ...
    'angleRmseDeg', 'fdRefRmseHz', 'nonRefCoherenceFloorMedian', ...
    'truthToothHitRate'}));
  summaryTable = localBuildDocumentSummaryTable(highlightTable);
  fprintf('\nCompact document summary rows:\n');
  disp(summaryTable);
end

%% Redraw compact paper figure

fig = localCreateFigure(layoutMode);
ax = gobjects(3, 1);

ax(1) = nexttile;
localPlotByMethod(mainSliceTable, 'angleRmseDeg', true);
localMarkHighlightRows(ax(1), highlightTable, 'angleRmseDeg');
grid on;
xlabel('Frame count P');
ylabel('Angle RMSE (deg)');
title('(a) DoA accuracy');
legend('Location', 'best');

ax(2) = nexttile;
localPlotByMethod(mainSliceTable, 'fdRefRmseHz', true);
localMarkHighlightRows(ax(2), highlightTable, 'fdRefRmseHz');
grid on;
xlabel('Frame count P');
ylabel('Reference-Doppler RMSE (Hz)');
title('(b) Doppler consistency');
legend('Location', 'best');

ax(3) = nexttile;
localPlotByMethod(mainSliceTable, 'nonRefCoherenceFloorMedian', false);
localMarkHighlightRows(ax(3), highlightTable, 'nonRefCoherenceFloorMedian');
grid on;
xlabel('Frame count P');
ylabel('Non-ref coherence floor');
ylim([0, 1.05]);
title('(c) Coherence consistency');
legend('Location', 'best');

sgtitle(sprintf('Controlled in-tooth CP/IP trade-off (SNR = %.0f dB, T_f = %.4g ms)', ...
  scanData.mainSnrDb, scanData.mainFrameIntvlSec * 1e3));

if exportFigure
  if ~isfolder(exportDir)
    mkdir(exportDir);
  end
  pngFile = fullfile(exportDir, exportBaseName + ".png");
  pdfFile = fullfile(exportDir, exportBaseName + ".pdf");
  exportgraphics(fig, pngFile, 'Resolution', 300);
  exportgraphics(fig, pdfFile, 'ContentType', 'vector');
  fprintf('Exported paper controlled CP/IP trade-off figure:\n  %s\n  %s\n', ...
    pngFile, pdfFile);
end

%% Local helpers


function summaryTable = localBuildDocumentSummaryTable(highlightTable)
%LOCALBUILDDOCUMENTSUMMARYTABLE Build a compact table suitable for notes.

summaryTable = table();
summaryTable.displayName = string(highlightTable.displayName);
summaryTable.numFrame = highlightTable.numFrame;
summaryTable.angleRmseDeg = highlightTable.angleRmseDeg;
summaryTable.fdRefRmseHz = highlightTable.fdRefRmseHz;
summaryTable.coherenceFloor = highlightTable.nonRefCoherenceFloorMedian;
summaryTable.truthToothHitRate = highlightTable.truthToothHitRate;
end

function scanData = localNormalizeScanData(scanData, mainFrameIntvlSec, mainSnrDb)
%LOCALNORMALIZESCANDATA Select the main paper slice from a controlled CP/IP scan.

if ~isfield(scanData, 'aggregateTable') || ~istable(scanData.aggregateTable)
  error('plotMfCpIpTradeoffPaper:InvalidScanData', ...
    'scanData.aggregateTable is missing or invalid.');
end

aggregateTable = scanData.aggregateTable;
aggregateTable = localEnsureAggregateColumns(aggregateTable);
aggregateTable.displayName = string(aggregateTable.displayName);
aggregateTable.phaseMode = string(aggregateTable.phaseMode);
aggregateTable.fdRateMode = string(aggregateTable.fdRateMode);
aggregateTable.frameIntvlMs = aggregateTable.frameIntvlSec * 1e3;
aggregateTable.windowMs = aggregateTable.windowSec * 1e3;

snrList = unique(aggregateTable.snrDb, 'stable');
[~, snrIdx] = min(abs(snrList - mainSnrDb));
mainSnrDb = snrList(snrIdx);
snrTable = aggregateTable(aggregateTable.snrDb == mainSnrDb, :);

frameIntvlList = unique(snrTable.frameIntvlSec, 'stable');
[~, tfIdx] = min(abs(frameIntvlList - mainFrameIntvlSec));
mainFrameIntvlSec = frameIntvlList(tfIdx);
mainSliceTable = snrTable(abs(snrTable.frameIntvlSec - mainFrameIntvlSec) < ...
  max(eps(mainFrameIntvlSec) * 16, 1e-15), :);

methodOrder = ["CP-K"; "CP-U"; "IP-K"; "IP-U"];
mainSliceTable.methodOrder = localMethodOrderIndex(mainSliceTable.displayName, methodOrder);
mainSliceTable = sortrows(mainSliceTable, {'methodOrder', 'numFrame'});

scanData.aggregateTable = aggregateTable;
scanData.mainSliceTable = mainSliceTable;
scanData.mainSnrDb = mainSnrDb;
scanData.mainFrameIntvlSec = mainFrameIntvlSec;
scanData.methodOrder = methodOrder;
end

function aggregateTable = localEnsureAggregateColumns(aggregateTable)
%LOCALENSUREAGGREGATECOLUMNS Validate required aggregate metrics.

requiredFieldList = ["displayName"; "phaseMode"; "fdRateMode"; "snrDb"; ...
  "numFrame"; "frameIntvlSec"; "angleRmseDeg"; ...
  "fdRefRmseHz"; "nonRefCoherenceFloorMedian"];
for iField = 1:numel(requiredFieldList)
  fieldName = char(requiredFieldList(iField));
  if ~ismember(fieldName, aggregateTable.Properties.VariableNames)
    error('plotMfCpIpTradeoffPaper:MissingColumn', ...
      'scanData.aggregateTable is missing required column: %s', fieldName);
  end
end
if ~ismember('windowSec', aggregateTable.Properties.VariableNames)
  aggregateTable.windowSec = (aggregateTable.numFrame - 1) .* aggregateTable.frameIntvlSec;
end
if ~ismember('truthToothHitRate', aggregateTable.Properties.VariableNames)
  aggregateTable.truthToothHitRate = nan(height(aggregateTable), 1);
end
end

function orderIdx = localMethodOrderIndex(displayName, methodOrder)
%LOCALMETHODORDERINDEX Return stable order keys for CP/IP methods.

orderIdx = nan(numel(displayName), 1);
for iMethod = 1:numel(methodOrder)
  orderIdx(displayName == methodOrder(iMethod)) = iMethod;
end
orderIdx(isnan(orderIdx)) = numel(methodOrder) + (1:nnz(isnan(orderIdx))).';
end

function highlightTable = localBuildHighlightTable(mainSliceTable, frameCountList)
%LOCALBUILDHIGHLIGHTTABLE Pick representative P values for command-window output.

frameCountList = reshape(double(frameCountList), [], 1);
rowIdx = false(height(mainSliceTable), 1);
for iFrame = 1:numel(frameCountList)
  rowIdx = rowIdx | (mainSliceTable.numFrame == frameCountList(iFrame));
end
highlightTable = mainSliceTable(rowIdx, :);
if isempty(highlightTable)
  warning('plotMfCpIpTradeoffPaper:NoHighlightRows', ...
    'No requested highlight frame counts were found in the selected slice.');
  return;
end
highlightTable = sortrows(highlightTable, {'methodOrder', 'numFrame'});
end

function fig = localCreateFigure(layoutMode)
%LOCALCREATEFIGURE Create the paper figure and tiled layout.

fig = figure('Name', 'Controlled CP/IP trade-off', 'Color', 'w');
if layoutMode == "vertical"
  fig.Position(3:4) = [780, 900];
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
  fig.Position(3:4) = [1320, 430];
  tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end
end

function localPlotByMethod(plotTable, yField, useLogY)
%LOCALPLOTBYMETHOD Plot one metric against frame count for CP/IP methods.

methodOrder = ["CP-K"; "CP-U"; "IP-K"; "IP-U"];
hold on;
for iMethod = 1:numel(methodOrder)
  mask = plotTable.displayName == methodOrder(iMethod);
  if ~any(mask)
    continue;
  end
  methodTable = sortrows(plotTable(mask, :), 'numFrame');
  yValue = methodTable.(yField);
  if useLogY
    yValue = max(yValue, eps);
    semilogy(methodTable.numFrame, yValue, '-o', 'LineWidth', 1.5, ...
      'DisplayName', localMethodLabel(methodOrder(iMethod)));
  else
    plot(methodTable.numFrame, yValue, '-o', 'LineWidth', 1.5, ...
      'DisplayName', localMethodLabel(methodOrder(iMethod)));
  end
end
end

function localMarkHighlightRows(ax, highlightTable, yFieldName)
%LOCALMARKHIGHLIGHTROWS Mark P=10/P=20 style representative points.

if isempty(highlightTable)
  return;
end
axes(ax); %#ok<LAXES>
if strcmp(get(ax, 'YScale'), 'log')
  yValue = max(highlightTable.(yFieldName), eps);
else
  yValue = highlightTable.(yFieldName);
end
scatter(highlightTable.numFrame, yValue, 42, 'filled', ...
  'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

function label = localMethodLabel(displayName)
%LOCALMETHODLABEL Format legend labels in English.

label = string(displayName);
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
  error('plotMfCpIpTradeoffPaper:RepoRootNotFound', ...
    'Cannot find repository root from paper script path.');
end
end
