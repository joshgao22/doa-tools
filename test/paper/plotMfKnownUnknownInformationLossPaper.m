%PLOTMFKNOWNUNKNOWNINFORMATIONLOSSPAPER Redraw the paper-facing known/unknown-rate EFIM loss figure.
% Load the existing scanMfKnownUnknownInformationLoss snapshot and build a compact
% English figure for the nuisance Doppler-rate information-loss discussion.

clear; close all; clc;

%% Paper figure configuration

scriptName = "plotMfKnownUnknownInformationLossPaper";
repoRoot = localFindRepoRoot();
addpath(fullfile(repoRoot, 'utils', 'io'));
addpath(fullfile(repoRoot, 'test', 'common', 'util'));

preferredSnapshotFile = fullfile(repoRoot, 'test', 'data', 'cache', 'scan', ...
  'scanMfKnownUnknownInformationLoss_20260428-151356.mat');
snapshotFile = resolveTestSnapshotFile(preferredSnapshotFile, ...
  'scanMfKnownUnknownInformationLoss_*.mat');
mainFrameIntvlSec = 1 / 750;
mainSnrDb = 10;
highlightFrameCountList = [10; 20];
layoutMode = "horizontal"; % "horizontal" for slides, "vertical" for notes.
exportFigure = false;
exportDir = fullfile(repoRoot, 'tmp', char(scriptName));
exportBaseName = "mfKnownUnknownInformationLossPaper";

%% Load scan snapshot

snapshotData = loadExpSnapshot(snapshotFile, 'none', struct('varNames', {{'scanData'}}));
if ~isfield(snapshotData, 'scanData') || ~isstruct(snapshotData.scanData)
  error('plotMfKnownUnknownInformationLossPaper:MissingScanData', ...
    'The snapshot does not contain a valid scanData struct.');
end
scanData = localNormalizeScanData(snapshotData.scanData, mainFrameIntvlSec, mainSnrDb);
mainSliceTable = scanData.mainSliceTable;
highlightTable = localBuildHighlightTable(mainSliceTable, highlightFrameCountList);

fprintf('\nPaper known/unknown-rate information-loss highlight rows:\n');
disp(highlightTable(:, {'satMode', 'snrDb', 'numFrame', 'frameIntvlMs', ...
  'windowMs', 'angleLossPct', 'fdRefLossPct', 'traceInfoLossPct'}));

summaryTable = localBuildDocumentSummaryTable(highlightTable);
fprintf('\nCompact document summary rows:\n');
disp(summaryTable);

%% Redraw compact paper figure

fig = localCreateFigure(layoutMode);
ax = gobjects(3, 1);

ax(1) = nexttile;
localPlotBySatMode(mainSliceTable, 'fdRefLossPct');
localMarkHighlightRows(ax(1), highlightTable, 'fdRefLossPct');
grid on;
xlabel('Frame count P');
ylabel('CRB std rollback (%)');
title('(a) Reference-Doppler loss');
legend('Location', 'northeast');

ax(2) = nexttile;
localPlotBySatMode(mainSliceTable, 'angleLossPct');
localMarkHighlightRows(ax(2), highlightTable, 'angleLossPct');
grid on;
xlabel('Frame count P');
ylabel('CRB std rollback (%)');
title('(b) DoA loss');
legend('Location', 'northwest');

ax(3) = nexttile;
localPlotBySatMode(mainSliceTable, 'traceInfoLossPct');
localMarkHighlightRows(ax(3), highlightTable, 'traceInfoLossPct');
grid on;
xlabel('Frame count P');
ylabel('EFIM trace loss (%)');
title('(c) Effective-information loss');
legend('Location', 'northeast');

sgtitle(sprintf('Known-rate vs unknown-rate EFIM loss (SNR = %.0f dB, T_f = %.4g ms)', ...
  scanData.mainSnrDb, scanData.mainFrameIntvlSec * 1e3));

if exportFigure
  if ~isfolder(exportDir)
    mkdir(exportDir);
  end
  pngFile = fullfile(exportDir, exportBaseName + ".png");
  pdfFile = fullfile(exportDir, exportBaseName + ".pdf");
  exportgraphics(fig, pngFile, 'Resolution', 300);
  exportgraphics(fig, pdfFile, 'ContentType', 'vector');
  fprintf('Exported paper known/unknown-rate information-loss figure:\n  %s\n  %s\n', ...
    pngFile, pdfFile);
end

%% Local helpers


function summaryTable = localBuildDocumentSummaryTable(highlightTable)
%LOCALBUILDDOCUMENTSUMMARYTABLE Build a compact table suitable for notes.

summaryTable = table();
summaryTable.satMode = string(highlightTable.satMode);
summaryTable.numFrame = highlightTable.numFrame;
summaryTable.windowMs = highlightTable.windowMs;
summaryTable.doaCrbRollbackPct = highlightTable.angleLossPct;
summaryTable.fdRefCrbRollbackPct = highlightTable.fdRefLossPct;
summaryTable.efimTraceLossPct = highlightTable.traceInfoLossPct;
end

function scanData = localNormalizeScanData(scanData, mainFrameIntvlSec, mainSnrDb)
%LOCALNORMALIZESCANDATA Select the main paper slice from a known/unknown-rate scan.

if ~isfield(scanData, 'lossTable') || ~istable(scanData.lossTable)
  error('plotMfKnownUnknownInformationLossPaper:InvalidScanData', ...
    'scanData.lossTable is missing or invalid.');
end

lossTable = scanData.lossTable;
lossTable = localEnsureLossColumns(lossTable);
lossTable.satMode = string(lossTable.satMode);
lossTable.timeOriginClass = string(lossTable.timeOriginClass);
lossTable.frameIntvlMs = lossTable.frameIntvlSec * 1e3;
lossTable.windowMs = lossTable.windowSec * 1e3;

snrList = unique(lossTable.snrDb, 'stable');
[~, snrIdx] = min(abs(snrList - mainSnrDb));
mainSnrDb = snrList(snrIdx);
snrTable = lossTable(lossTable.snrDb == mainSnrDb, :);

frameIntvlList = unique(snrTable.frameIntvlSec, 'stable');
[~, tfIdx] = min(abs(frameIntvlList - mainFrameIntvlSec));
mainFrameIntvlSec = frameIntvlList(tfIdx);
mainSliceTable = snrTable(abs(snrTable.frameIntvlSec - mainFrameIntvlSec) < ...
  max(eps(mainFrameIntvlSec) * 16, 1e-15), :);

originClass = localPickMainTimeOriginClass(mainSliceTable);
mainSliceTable = mainSliceTable(mainSliceTable.timeOriginClass == originClass, :);
mainSliceTable = sortrows(mainSliceTable, {'satMode', 'numFrame'});

scanData.lossTable = lossTable;
scanData.mainSliceTable = mainSliceTable;
scanData.mainTimeOriginClass = originClass;
scanData.mainSnrDb = mainSnrDb;
scanData.mainFrameIntvlSec = mainFrameIntvlSec;
end

function lossTable = localEnsureLossColumns(lossTable)
%LOCALENSURELOSSCOLUMNS Rebuild readable loss columns for older snapshots.

if ~ismember('angleLossPct', lossTable.Properties.VariableNames)
  lossTable.angleLossPct = 100 * (lossTable.angleStdLossRatio - 1);
end
if ~ismember('fdRefLossPct', lossTable.Properties.VariableNames)
  lossTable.fdRefLossPct = 100 * (lossTable.fdStdLossRatio - 1);
end
if ~ismember('traceInfoLossPct', lossTable.Properties.VariableNames)
  lossTable.traceInfoLossPct = 100 * (1 - lossTable.traceInfoRetention);
end
if ~ismember('minEigLossPct', lossTable.Properties.VariableNames)
  lossTable.minEigLossPct = 100 * (1 - lossTable.minEigRetention);
end
end

function originClass = localPickMainTimeOriginClass(sliceTable)
%LOCALPICKMAINTIMEORIGINCLASS Prefer the repository default right-biased window.

originList = unique(sliceTable.timeOriginClass, 'stable');
if any(originList == "rightBiasedRef")
  originClass = "rightBiasedRef";
elseif ~isempty(originList)
  originClass = originList(1);
else
  originClass = "";
end
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
  warning('plotMfKnownUnknownInformationLossPaper:NoHighlightRows', ...
    'No requested highlight frame counts were found in the selected slice.');
  return;
end
highlightTable = sortrows(highlightTable, {'satMode', 'numFrame'});
end

function fig = localCreateFigure(layoutMode)
%LOCALCREATEFIGURE Create the paper figure and tiled layout.

fig = figure('Name', 'Known/unknown-rate EFIM information loss', 'Color', 'w');
if layoutMode == "vertical"
  fig.Position(3:4) = [780, 900];
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
  fig.Position(3:4) = [1320, 430];
  tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end
end

function localPlotBySatMode(plotTable, yField)
%LOCALPLOTBYSATMODE Plot one metric against frame count for single and multi modes.

groupList = unique(plotTable.satMode, 'stable');
hold on;
for iGroup = 1:numel(groupList)
  mask = plotTable.satMode == groupList(iGroup);
  groupTable = sortrows(plotTable(mask, :), 'numFrame');
  plot(groupTable.numFrame, groupTable.(yField), '-o', 'LineWidth', 1.5, ...
    'DisplayName', localSatModeLabel(groupList(iGroup)));
end
end

function localMarkHighlightRows(ax, highlightTable, yFieldName)
%LOCALMARKHIGHLIGHTROWS Mark P=10/P=20 style representative points.

if isempty(highlightTable)
  return;
end
axes(ax); %#ok<LAXES>
scatter(highlightTable.numFrame, highlightTable.(yFieldName), 42, 'filled', ...
  'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

function label = localSatModeLabel(satMode)
%LOCALSATMODELABEL Format legend labels in English.

satMode = string(satMode);
if satMode == "single"
  label = "single-satellite";
elseif satMode == "multi"
  label = "multi-satellite";
else
  label = satMode;
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
  error('plotMfKnownUnknownInformationLossPaper:RepoRootNotFound', ...
    'Cannot find repository root from paper script path.');
end
end
