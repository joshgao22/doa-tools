%PLOTSYNCBLOCKOBSERVATIONWINDOWPAPER Draw the sync-block and multi-frame window schematic.
% This paper-facing script is intentionally data-free. It visualizes the
% Starlink-inspired frame-start known sync blocks, the multi-frame
% observation window, and the CP/IP phase-model distinction used in the
% paper positioning notes.

clear; close all; clc;

%% Paper figure configuration

scriptName = "plotSyncBlockObservationWindowPaper";
repoRoot = localFindRepoRoot();

frameCount = 10;
frameIntvlSec = 1 / 750;
fsHz = 240e6;
ofdmSize = 1024;
cpSize = 32;
numSyncOfdmSymbol = 2;
numSyncSample = numSyncOfdmSymbol * (ofdmSize + cpSize);

% The real sync block is much narrower than one frame. Enlarge it only for
% visual readability while keeping the numeric annotation exact.
syncWidthScale = 30;
maxFrameToLabel = 5;
labelLanguage = "en"; % "zh" or "en".
exportFigure = false;
exportDir = fullfile(repoRoot, 'tmp', char(scriptName));
exportBaseName = "syncBlockObservationWindowPaper";

%% Derived quantities

frameIntvlMs = frameIntvlSec * 1e3;
syncDurUs = numSyncSample / fsHz * 1e6;
windowMs = frameCount * frameIntvlMs;
syncWidthMsTrue = syncDurUs / 1e3;
syncWidthMsDraw = min(0.22 * frameIntvlMs, syncWidthMsTrue * syncWidthScale);
frameStartMs = (0:frameCount-1).' * frameIntvlMs;
syncCenterMsDraw = frameStartMs + 0.5 * syncWidthMsDraw;

fprintf('\nSync-block observation-window figure summary:\n');
fprintf('  frameCount              = %d\n', frameCount);
fprintf('  frameIntvlSec           = %.12g s\n', frameIntvlSec);
fprintf('  windowMs                = %.6g ms\n', windowMs);
fprintf('  fsHz                    = %.6g Hz\n', fsHz);
fprintf('  numSyncSample           = %d\n', numSyncSample);
fprintf('  syncBlockDurationUs     = %.6g us\n', syncDurUs);
fprintf('  drawnSyncWidthScale     = %.6g (for visualization only)\n', syncWidthScale);

%% Draw schematic

fig = figure('Name', 'Sync-block multi-frame observation window', ...
  'Color', 'w', 'Position', [80, 80, 1200, 680]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
localDrawFrameWindow(ax1, frameStartMs, frameIntvlMs, syncWidthMsDraw, ...
  windowMs, maxFrameToLabel, labelLanguage);

ax2 = nexttile(tl, 2);
localDrawContinuousPhase(ax2, frameStartMs, syncWidthMsDraw, frameIntvlMs, ...
  windowMs, labelLanguage);

ax3 = nexttile(tl, 3);
localDrawIndependentPhase(ax3, frameStartMs, syncWidthMsDraw, frameIntvlMs, ...
  windowMs, labelLanguage);

if labelLanguage == "zh"
  title(tl, '帧首同步块与多帧观测窗口示意图');
else
  title(tl, 'Frame-start sync blocks and multi-frame observation window');
end

annotationText = sprintf('T_f = %.4g ms,  P = %d,  N_p = %d samples = %.3g us at F_s = %.3g MHz', ...
  frameIntvlMs, frameCount, numSyncSample, syncDurUs, fsHz / 1e6);
annotation(fig, 'textbox', [0.12, 0.005, 0.78, 0.04], ...
  'String', annotationText, 'EdgeColor', 'none', ...
  'HorizontalAlignment', 'center', 'FontSize', 10);

if exportFigure
  if ~exist(exportDir, 'dir')
    mkdir(exportDir);
  end
  pngFile = fullfile(exportDir, exportBaseName + ".png");
  pdfFile = fullfile(exportDir, exportBaseName + ".pdf");
  exportgraphics(fig, pngFile, 'Resolution', 300);
  exportgraphics(fig, pdfFile, 'ContentType', 'vector');
  fprintf('Exported sync-block paper figure:\n  %s\n  %s\n', pngFile, pdfFile);
end

%% Local helpers

function localDrawFrameWindow(ax, frameStartMs, frameIntvlMs, syncWidthMsDraw, ...
  windowMs, maxFrameToLabel, labelLanguage)
  axes(ax); %#ok<LAXES>
  cla(ax); hold(ax, 'on');
  numFrame = numel(frameStartMs);
  frameHeight = 0.6;
  syncHeight = 0.6;

  for p = 1:numFrame
    x0 = frameStartMs(p);
    rectangle(ax, 'Position', [x0, 0.2, frameIntvlMs, frameHeight], ...
      'FaceColor', [0.94, 0.94, 0.94], 'EdgeColor', [0.65, 0.65, 0.65], ...
      'LineWidth', 0.8);
    rectangle(ax, 'Position', [x0, 0.2, syncWidthMsDraw, syncHeight], ...
      'FaceColor', [0.80, 0.88, 1.00], 'EdgeColor', [0.25, 0.35, 0.70], ...
      'LineWidth', 1.0);

    if p <= maxFrameToLabel || p == numFrame
      text(ax, x0 + 0.5 * frameIntvlMs, 0.92, sprintf('p=%d', p-1), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    end
  end

  plot(ax, [0, windowMs], [0.08, 0.08], 'k-', 'LineWidth', 1.0);
  plot(ax, [0, 0], [0.04, 0.12], 'k-', 'LineWidth', 1.0);
  plot(ax, [windowMs, windowMs], [0.04, 0.12], 'k-', 'LineWidth', 1.0);

  localStyleWindowAxis(ax, windowMs);
  ylim(ax, [0, 1.18]);
  xlabel(ax, localLabel(labelLanguage, 'time (ms)', '时间 (ms)'));
  ylabel(ax, localLabel(labelLanguage, 'frames', '帧'));
  title(ax, localLabel(labelLanguage, ...
    '(a) Known sync blocks are extracted at the beginning of each frame', ...
    '(a) 每帧帧首提取已知同步块'));
  text(ax, 0.01 * windowMs, 0.02, localLabel(labelLanguage, ...
    'multi-frame observation window P T_f', ...
    '多帧观测窗口 P T_f'), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 9);
  text(ax, frameStartMs(1) + syncWidthMsDraw, 0.5, localLabel(labelLanguage, ...
    'sync block width enlarged for display', ...
    '同步块宽度为示意性放大'), ...
    'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', [0.2, 0.2, 0.2]);
end

function localDrawContinuousPhase(ax, frameStartMs, syncWidthMsDraw, frameIntvlMs, ...
  windowMs, labelLanguage)
  axes(ax); %#ok<LAXES>
  cla(ax); hold(ax, 'on');
  numFrame = numel(frameStartMs);
  phaseFun = @(x) 0.18 + 0.45 * (x / windowMs) + 0.30 * (x / windowMs).^2;
  xDense = linspace(0, windowMs, 500);
  plot(ax, xDense, phaseFun(xDense), '-', 'LineWidth', 1.6, 'Color', [0.10, 0.30, 0.75]);

  for p = 1:numFrame
    xBlock = linspace(frameStartMs(p), frameStartMs(p) + syncWidthMsDraw, 20);
    yBlock = phaseFun(xBlock);
    plot(ax, xBlock, yBlock, '-', 'LineWidth', 5.0, 'Color', [0.10, 0.30, 0.75]);
    plot(ax, frameStartMs(p) + 0.5 * syncWidthMsDraw, phaseFun(frameStartMs(p) + 0.5 * syncWidthMsDraw), ...
      'o', 'MarkerSize', 4.5, 'MarkerFaceColor', [0.10, 0.30, 0.75], ...
      'MarkerEdgeColor', 'none');
  end

  localStyleWindowAxis(ax, windowMs);
  ylim(ax, [0, 1.05]);
  xlabel(ax, localLabel(labelLanguage, 'global time t_{p,t}', '全局时间 t_{p,t}'));
  ylabel(ax, localLabel(labelLanguage, 'carrier phase', '载波相位'));
  title(ax, localLabel(labelLanguage, ...
    '(b) CP model: one shared per-satellite phase and global time index', ...
    '(b) 连续相位模型：同一卫星共享公共相位并使用全局时间索引'));
  text(ax, 0.58 * windowMs, 0.78, ...
    '$\chi_\ell + \omega_\ell t_{p,t} + \frac{1}{2}\gamma_\ell t_{p,t}^2$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Color', [0.10, 0.30, 0.75]);
end

function localDrawIndependentPhase(ax, frameStartMs, syncWidthMsDraw, frameIntvlMs, ...
  windowMs, labelLanguage)
  axes(ax); %#ok<LAXES>
  cla(ax); hold(ax, 'on');
  numFrame = numel(frameStartMs);
  blockSlope = 0.16 / syncWidthMsDraw;
  phaseOffsets = 0.23 + 0.20 * sin((1:numFrame).' * 1.7) + 0.035 * (1:numFrame).';

  for p = 1:numFrame
    xBlock = linspace(frameStartMs(p), frameStartMs(p) + syncWidthMsDraw, 20);
    yBlock = phaseOffsets(p) + blockSlope * (xBlock - frameStartMs(p));
    plot(ax, xBlock, yBlock, '-', 'LineWidth', 5.0, 'Color', [0.78, 0.30, 0.10]);
    plot(ax, frameStartMs(p) + 0.5 * syncWidthMsDraw, phaseOffsets(p) + 0.5 * blockSlope * syncWidthMsDraw, ...
      'o', 'MarkerSize', 4.5, 'MarkerFaceColor', [0.78, 0.30, 0.10], ...
      'MarkerEdgeColor', 'none');
    if p <= 3
      text(ax, frameStartMs(p), yBlock(1) + 0.08, sprintf('$\\vartheta_{\\ell,%d}$', p-1), ...
        'Interpreter', 'latex', 'FontSize', 10, 'Color', [0.78, 0.30, 0.10]);
    end
  end

  localStyleWindowAxis(ax, windowMs);
  ylim(ax, [0, 1.05]);
  xlabel(ax, localLabel(labelLanguage, 'local time t in each frame', '每帧局部时间 t'));
  ylabel(ax, localLabel(labelLanguage, 'carrier phase', '载波相位'));
  title(ax, localLabel(labelLanguage, ...
    '(c) IP baseline: each satellite-frame block has an independent phase', ...
    '(c) 独立相位基线：每个卫星-帧观测块具有独立未知相位'));
  text(ax, 0.50 * windowMs, 0.78, ...
    '$\vartheta_{\ell,p} + \bar{\omega}_{\ell,p}t + \frac{1}{2}\gamma_\ell t^2$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Color', [0.78, 0.30, 0.10]);
end

function localStyleWindowAxis(ax, windowMs)
  xlim(ax, [-0.02 * windowMs, 1.02 * windowMs]);
  set(ax, 'YTick', [], 'Box', 'off', 'Layer', 'top');
  grid(ax, 'on');
end

function label = localLabel(labelLanguage, labelEn, labelZh)
  if labelLanguage == "zh"
    label = labelZh;
  else
    label = labelEn;
  end
end

function repoRoot = localFindRepoRoot()
  thisFile = mfilename('fullpath');
  thisDir = fileparts(thisFile);
  repoRoot = thisDir;
  while ~isfile(fullfile(repoRoot, 'README.md')) || ~isfile(fullfile(repoRoot, 'AGENTS.md'))
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
      error('plotSyncBlockObservationWindowPaper:RepoRootNotFound', ...
        'Cannot find repository root from %s.', thisDir);
    end
    repoRoot = parentDir;
  end
end
