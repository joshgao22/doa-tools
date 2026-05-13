%PLOTMFMLECRBCONSISTENCYPAPER Redraw trimmed in-tooth MLE/CRB consistency figures.
% Load the existing SS and MS in-tooth scan snapshots and compare only the
% trimmed local-regime estimator error with the corresponding CRB curves.

clear; close all; clc;

%% Paper figure configuration

scriptName = "plotMfMleCrbConsistencyPaper";
repoRoot = localFindRepoRoot();
addpath(fullfile(repoRoot, 'utils', 'io'));

snapshotSpec = struct( ...
  'datasetName', {"SS-MF", "MS-MF"}, ...
  'snapshotFile', { ...
    fullfile(repoRoot, 'test', 'data', 'cache', 'scan', ...
      'scanMfMleCrbInToothConsistency_20260507-222606.mat'), ...
    fullfile(repoRoot, 'test', 'data', 'cache', 'scan', ...
      'scanMfMsMleCrbInToothConsistency_20260510-062927.mat')});

mainFrameCount = 10;
mainFrameIntvlSec = 1 / 750;
methodNameList = ["SS-MF-CP-K", "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"];
plotMetric = "rmse"; % "rmse", "mse", or "mseOverCrb".
layoutMode = "horizontal"; % "horizontal" for slides, "vertical" for notes.
exportFigure = false;
exportDir = fullfile(repoRoot, 'tmp', char(scriptName));
exportBaseName = "mfMleCrbConsistencyPaper";

%% Load scan snapshots

paperTable = table();
for iData = 1:numel(snapshotSpec)
  snapshotFile = snapshotSpec(iData).snapshotFile;
  if ~isfile(snapshotFile)
    error('plotMfMleCrbConsistencyPaper:MissingSnapshot', ...
      'Snapshot file is missing: %s', snapshotFile);
  end

  snapshotData = loadExpSnapshot(snapshotFile, 'none', struct('varNames', {{'scanData'}}));
  if ~isfield(snapshotData, 'scanData') || ~isstruct(snapshotData.scanData)
    error('plotMfMleCrbConsistencyPaper:MissingScanData', ...
      'The snapshot does not contain a valid scanData struct: %s', snapshotFile);
  end

  aggregateTable = localGetAggregateTable(snapshotData.scanData, snapshotSpec(iData).datasetName);
  datasetTable = localBuildPaperTable(aggregateTable, snapshotSpec(iData).datasetName, ...
    mainFrameCount, mainFrameIntvlSec, methodNameList);
  paperTable = [paperTable; datasetTable]; %#ok<AGROW>
end

if isempty(paperTable) || height(paperTable) == 0
  error('plotMfMleCrbConsistencyPaper:EmptyPaperTable', ...
    'No rows matched the requested dataset/method/frame configuration.');
end

fprintf('\nPaper trimmed MLE/CRB consistency compact table:\n');
disp(localCompactPrintTable(paperTable));

%% Redraw compact paper figure

fig = localCreateFigure(layoutMode);

ax = gobjects(4, 1);
plotSpec = [ ...
  struct('datasetName', "SS-MF", 'metricName', "angle", ...
    'titleText', "(a) SS-MF DoA", ...
    'subtitleText', "Trimmed local-regime error vs CRB"), ...
  struct('datasetName', "SS-MF", 'metricName', "fdRef", ...
    'titleText', "(b) SS-MF reference Doppler", ...
    'subtitleText', "Trimmed local-regime error vs CRB"), ...
  struct('datasetName', "MS-MF", 'metricName', "angle", ...
    'titleText', "(c) MS-MF DoA", ...
    'subtitleText', "Trimmed local-regime error vs CRB"), ...
  struct('datasetName', "MS-MF", 'metricName', "fdRef", ...
    'titleText', "(d) MS-MF reference Doppler", ...
    'subtitleText', "Trimmed local-regime error vs CRB")];

for iPlot = 1:numel(plotSpec)
  ax(iPlot) = nexttile;
  localPlotConsistencyPanel(ax(iPlot), paperTable, plotSpec(iPlot), plotMetric);
end

sgtitle(localBuildFigureTitle(plotMetric, mainFrameCount, mainFrameIntvlSec));

if exportFigure
  if ~isfolder(exportDir)
    mkdir(exportDir);
  end
  pngFile = fullfile(exportDir, exportBaseName + ".png");
  pdfFile = fullfile(exportDir, exportBaseName + ".pdf");
  exportgraphics(fig, pngFile, 'Resolution', 300);
  exportgraphics(fig, pdfFile, 'ContentType', 'vector');
  fprintf('Exported MLE/CRB consistency figure:\n  %s\n  %s\n', pngFile, pdfFile);
end

%% Local helpers

function aggregateTable = localGetAggregateTable(scanData, datasetName)
%LOCALGETAGGREGATETABLE Validate and fetch the scan aggregate table.

if ~isfield(scanData, 'aggregateTable') || ~istable(scanData.aggregateTable)
  error('plotMfMleCrbConsistencyPaper:InvalidScanData', ...
    '%s scanData.aggregateTable is missing or invalid.', datasetName);
end
aggregateTable = scanData.aggregateTable;
requiredFields = ["displayName", "snrDb", "numFrame", "frameIntvlSec", ...
  "trimAngleMseDeg2", "trimAngleRmseDeg", "angleCrbStdDeg", ...
  "trimAngleRmseOverCrb", "trimAngleMseOverCrb", ...
  "trimFdRefMseHz2", "trimFdRefRmseHz", "fdRefCrbStdHz", ...
  "trimFdRefRmseOverCrb", "trimFdRefMseOverCrb", ...
  "trimKeepRate", "crbLocalRate"];
missingFields = setdiff(requiredFields, string(aggregateTable.Properties.VariableNames));
if ~isempty(missingFields)
  error('plotMfMleCrbConsistencyPaper:MissingAggregateFields', ...
    '%s aggregateTable misses fields: %s', datasetName, strjoin(missingFields, ', '));
end
end

function paperTable = localBuildPaperTable(aggregateTable, datasetName, mainFrameCount, ...
  mainFrameIntvlSec, methodNameList)
%LOCALBUILDPAPERTABLE Convert aggregate rows into paper-plot rows.

aggregateTable.displayName = string(aggregateTable.displayName);
methodNameList = reshape(string(methodNameList), 1, []);

mask = aggregateTable.numFrame == mainFrameCount & ...
  abs(aggregateTable.frameIntvlSec - mainFrameIntvlSec) < 1e-12 & ...
  ismember(aggregateTable.displayName, methodNameList);
aggregateTable = aggregateTable(mask, :);

rowCell = {};
for iRow = 1:height(aggregateTable)
  sourceRow = aggregateTable(iRow, :);
  row = table();
  row.datasetName = datasetName;
  row.displayName = string(sourceRow.displayName);
  row.snrDb = sourceRow.snrDb;
  row.numFrame = sourceRow.numFrame;
  row.frameIntvlSec = sourceRow.frameIntvlSec;
  row.trimAngleMseDeg2 = sourceRow.trimAngleMseDeg2;
  row.trimAngleRmseDeg = sourceRow.trimAngleRmseDeg;
  row.angleCrbStdDeg = sourceRow.angleCrbStdDeg;
  row.angleCrbVarDeg2 = sourceRow.angleCrbStdDeg.^2;
  row.trimAngleRmseOverCrb = sourceRow.trimAngleRmseOverCrb;
  row.trimAngleMseOverCrb = sourceRow.trimAngleMseOverCrb;
  row.trimFdRefMseHz2 = sourceRow.trimFdRefMseHz2;
  row.trimFdRefRmseHz = sourceRow.trimFdRefRmseHz;
  row.fdRefCrbStdHz = sourceRow.fdRefCrbStdHz;
  row.fdRefCrbVarHz2 = sourceRow.fdRefCrbStdHz.^2;
  row.trimFdRefRmseOverCrb = sourceRow.trimFdRefRmseOverCrb;
  row.trimFdRefMseOverCrb = sourceRow.trimFdRefMseOverCrb;
  row.trimKeepRate = sourceRow.trimKeepRate;
  row.crbLocalRate = sourceRow.crbLocalRate;
  row.numRepeat = localReadOptionalScalar(sourceRow, 'numRepeat', NaN);
  row.trimmedCoreCount = localReadOptionalScalar(sourceRow, 'trimmedCoreCount', NaN);
  row.crbLocalCount = localReadOptionalScalar(sourceRow, 'crbLocalCount', NaN);
  rowCell{end + 1, 1} = row; %#ok<AGROW>
end

if isempty(rowCell)
  paperTable = table();
else
  paperTable = vertcat(rowCell{:});
  paperTable = sortrows(paperTable, {'datasetName', 'displayName', 'snrDb'});
end
end

function value = localReadOptionalScalar(sourceRow, fieldName, defaultValue)
%LOCALREADOPTIONALSCALAR Read a scalar table field if present.

if ismember(fieldName, sourceRow.Properties.VariableNames)
  value = sourceRow.(fieldName);
else
  value = defaultValue;
end
end

function fig = localCreateFigure(layoutMode)
%LOCALCREATEFIGURE Create the paper figure and tiled layout.

fig = figure('Name', 'Paper trimmed in-tooth MLE/CRB consistency', 'Color', 'w');
if layoutMode == "vertical"
  fig.Position(3:4) = [860, 920];
  tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
  fig.Position(3:4) = [1320, 760];
  tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end
end

function localPlotConsistencyPanel(ax, paperTable, plotSpec, plotMetric)
%LOCALPLOTCONSISTENCYPANEL Plot trimmed estimator error and CRB curves.

axes(ax); %#ok<LAXES>
hold on;
datasetMask = paperTable.datasetName == plotSpec.datasetName;
methodList = unique(paperTable.displayName(datasetMask), 'stable');
markerList = {'o', 's', '^', 'd'};
legendText = strings(0, 1);
legendHandle = gobjects(0, 1);

for iMethod = 1:numel(methodList)
  methodName = methodList(iMethod);
  mask = datasetMask & paperTable.displayName == methodName;
  curveTable = paperTable(mask, :);
  if isempty(curveTable) || height(curveTable) == 0
    continue;
  end
  curveTable = sortrows(curveTable, 'snrDb');
  [trimValue, crbValue] = localGetPlotValues(curveTable, plotSpec.metricName, plotMetric);
  trimLine = semilogy(curveTable.snrDb, trimValue, ...
    'LineWidth', 1.35, ...
    'LineStyle', '-', ...
    'Marker', markerList{mod(iMethod - 1, numel(markerList)) + 1});
  legendText(end + 1, 1) = methodName + " trim"; %#ok<AGROW>
  legendHandle(end + 1, 1) = trimLine; %#ok<AGROW>

  if plotMetric ~= "mseOverCrb"
    crbLine = semilogy(curveTable.snrDb, crbValue, ...
      'LineWidth', 1.1, ...
      'LineStyle', ':', ...
      'Marker', 'x', ...
      'Color', trimLine.Color);
    legendText(end + 1, 1) = methodName + " CRB"; %#ok<AGROW>
    legendHandle(end + 1, 1) = crbLine; %#ok<AGROW>
  end
end

if plotMetric == "mseOverCrb"
  yline(1, 'k:', 'CRB', 'LabelHorizontalAlignment', 'left');
end

grid on;
set(ax, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel(localBuildYLabel(plotSpec.metricName, plotMetric));
title({plotSpec.titleText; localBuildSubtitle(plotSpec.metricName, plotMetric)});
if ~isempty(legendHandle)
  legend(legendHandle, legendText, 'Location', 'best', 'Interpreter', 'none');
end
localApplyLogYLim(ax);
end

function [trimValue, crbValue] = localGetPlotValues(curveTable, metricName, plotMetric)
%LOCALGETPLOTVALUES Return trimmed estimator and CRB values for one panel.

switch metricName
  case "angle"
    switch plotMetric
      case "rmse"
        trimValue = curveTable.trimAngleRmseDeg;
        crbValue = curveTable.angleCrbStdDeg;
      case "mse"
        trimValue = curveTable.trimAngleMseDeg2;
        crbValue = curveTable.angleCrbVarDeg2;
      case "mseOverCrb"
        trimValue = curveTable.trimAngleMseOverCrb;
        crbValue = ones(height(curveTable), 1);
      otherwise
        error('plotMfMleCrbConsistencyPaper:UnknownPlotMetric', ...
          'Unknown plotMetric: %s', plotMetric);
    end
  case "fdRef"
    switch plotMetric
      case "rmse"
        trimValue = curveTable.trimFdRefRmseHz;
        crbValue = curveTable.fdRefCrbStdHz;
      case "mse"
        trimValue = curveTable.trimFdRefMseHz2;
        crbValue = curveTable.fdRefCrbVarHz2;
      case "mseOverCrb"
        trimValue = curveTable.trimFdRefMseOverCrb;
        crbValue = ones(height(curveTable), 1);
      otherwise
        error('plotMfMleCrbConsistencyPaper:UnknownPlotMetric', ...
          'Unknown plotMetric: %s', plotMetric);
    end
  otherwise
    error('plotMfMleCrbConsistencyPaper:UnknownMetricName', ...
      'Unknown metricName: %s', metricName);
end
end

function yLabelText = localBuildYLabel(metricName, plotMetric)
%LOCALBUILDYLABEL Build axis label for actual-unit or normalized plots.

if metricName == "angle"
  switch plotMetric
    case "rmse"
      yLabelText = 'DoA RMSE and CRB std (deg)';
    case "mse"
      yLabelText = 'DoA MSE and CRB variance (deg^2)';
    case "mseOverCrb"
      yLabelText = 'DoA MSE / CRB';
  end
else
  switch plotMetric
    case "rmse"
      yLabelText = 'fdRef RMSE and CRB std (Hz)';
    case "mse"
      yLabelText = 'fdRef MSE and CRB variance (Hz^2)';
    case "mseOverCrb"
      yLabelText = 'fdRef MSE / CRB';
  end
end
end

function subtitleText = localBuildSubtitle(metricName, plotMetric)
%LOCALBUILDSUBTITLE Describe the metric and trimming rule in the subplot title.

if plotMetric == "mseOverCrb"
  if metricName == "angle"
    subtitleText = "trimmed DoA MSE normalized by CRB variance";
  else
    subtitleText = "trimmed fdRef MSE normalized by CRB variance";
  end
elseif plotMetric == "mse"
  subtitleText = "trimmed MSE and CRB variance, log scale";
else
  subtitleText = "trimmed RMSE and CRB std, log scale";
end
end

function titleText = localBuildFigureTitle(plotMetric, mainFrameCount, mainFrameIntvlSec)
%LOCALBUILDFIGURETITLE Build figure-wide title.

frameIntvlText = sprintf('T_f = %.6g ms', mainFrameIntvlSec * 1e3);
switch plotMetric
  case "rmse"
    metricText = "Trimmed in-tooth RMSE vs CRB";
  case "mse"
    metricText = "Trimmed in-tooth MSE vs CRB";
  case "mseOverCrb"
    metricText = "Trimmed in-tooth MSE/CRB";
  otherwise
    metricText = "Trimmed in-tooth MLE/CRB consistency";
end
titleText = sprintf('%s, P = %d, %s', metricText, mainFrameCount, frameIntvlText);
end

function localApplyLogYLim(ax)
%LOCALAPPLYLOGYLIM Keep log-y axes readable for finite positive data.

yData = [];
lineObj = findobj(ax, 'Type', 'Line');
for iLine = 1:numel(lineObj)
  yData = [yData; lineObj(iLine).YData(:)]; %#ok<AGROW>
end
yData = yData(isfinite(yData) & yData > 0);
if isempty(yData)
  return;
end
yMin = min(yData);
yMax = max(yData);
if ~(isfinite(yMin) && isfinite(yMax) && yMax > 0)
  return;
end
if yMin == yMax
  yMin = yMin / sqrt(10);
  yMax = yMax * sqrt(10);
else
  yMin = 10 .^ (floor(log10(yMin)) - 0.1);
  yMax = 10 .^ (ceil(log10(yMax)) + 0.1);
end
ylim(ax, [yMin, yMax]);
end

function printTable = localCompactPrintTable(paperTable)
%LOCALCOMPACTPRINTTABLE Keep command output compact and useful.

printTable = paperTable(:, {'datasetName', 'displayName', 'snrDb', ...
  'trimAngleRmseDeg', 'angleCrbStdDeg', 'trimAngleRmseOverCrb', ...
  'trimAngleMseOverCrb', 'trimFdRefRmseHz', 'fdRefCrbStdHz', ...
  'trimFdRefRmseOverCrb', 'trimFdRefMseOverCrb', ...
  'trimKeepRate', 'crbLocalRate', 'numRepeat', 'trimmedCoreCount'});
printTable = sortrows(printTable, {'datasetName', 'displayName', 'snrDb'});
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
  error('plotMfMleCrbConsistencyPaper:RepoRootNotFound', ...
    'Could not find repository root containing AGENTS.md.');
end
end
