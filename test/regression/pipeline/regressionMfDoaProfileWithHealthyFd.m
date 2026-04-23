% Regression diagnostic for the simplified dynamic flow.
% Find one MS-MF dynamic case whose fd chain is already healthy, then scan
% local DoA objective meshes for:
%   1) SS-SF-Static
%   2) MS-SF-Static
%   3) MS-MF-Subset-CP (selected non-periodic subset)
%   4) MS-MF-Periodic-CP
% using the original objective scale (no reciprocal transform).
clear(); close all;
localAddProjectPath();

context = buildDynamicDualSatEciContext(struct( ...
  'baseSeed', 253, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', struct('enableSubsetEvalParfor', true, 'minSubsetEvalParfor', 4)));
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', context.parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true));

snrDb = 10;
seedList = 253:264;
numSeed = numel(seedList);
healthyCaseCell = cell(numSeed, 1);
healthyScoreList = -inf(numSeed, 1);

seedTracker = localCreateProgressTracker('Searching healthy-fd seeds', numSeed);
cleanupSeedTracker = onCleanup(@() localCloseProgressTracker(seedTracker));
parfor iSeed = 1:numSeed
  taskSeed = seedList(iSeed);
  repeatData = buildDynamicRepeatData(context, snrDb, taskSeed);
  staticBundle = buildDoaDopplerStaticTransitionBundle( ...
    repeatData.periodicFixture.viewRefOnly, repeatData.periodicFixture.viewOtherOnly, repeatData.periodicFixture.viewMs, ...
    context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
    repeatData.periodicFixture.fdRange, repeatData.periodicFixture.truth, repeatData.periodicFixture.otherSatIdxGlobal, ...
    false, flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);
  msFlow = runSimpleDynamicSubsetPeriodicFlow("MS-MF-Dyn", "multi", ...
    repeatData.periodicFixture, repeatData.subsetFixtureCell, staticBundle.caseStaticMs, ...
    context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, false, flowOpt);
  summaryUse = msFlow.finalSummary;
  if logical(summaryUse.isResolved)
    fdHealthy = isfinite(summaryUse.fdRefErrHz) && isfinite(summaryUse.fdRateErrHzPerSec) && ...
      abs(summaryUse.fdRefErrHz) <= 5 && abs(summaryUse.fdRateErrHzPerSec) <= 50;
    if fdHealthy
      healthyScoreList(iSeed) = -localGetFieldOrDefault(summaryUse, 'angleErrDeg', inf);
      healthyCaseCell{iSeed} = struct('repeatData', repeatData, 'staticBundle', staticBundle, 'msFlow', msFlow);
    end
  end
  send(seedTracker.queue, 1);
end
clear cleanupSeedTracker;

[selectedScore, bestIdx] = max(healthyScoreList);
if ~isfinite(selectedScore) || isempty(healthyCaseCell{bestIdx})
  error('regressionMfDoaProfileWithHealthyFd:NoHealthyFdCase', ...
    'No healthy-fd MS-MF dynamic case was found in the tested seed list.');
end
selectedRep = healthyCaseCell{bestIdx};

truth = selectedRep.repeatData.truth;
ssStaticCase = selectedRep.staticBundle.caseStaticRefOnly;
msStaticCase = selectedRep.staticBundle.caseStaticMs;
subsetCase = selectedRep.msFlow.selectedSubsetCase;
periodicCase = selectedRep.msFlow.caseFinal;
periodicSummary = selectedRep.msFlow.finalSummary;
selectedSubsetFixture = selectedRep.repeatData.subsetFixtureCell{selectedRep.msFlow.bestSubsetIdx};
selectedSubsetLabel = string(localGetFieldOrDefault(selectedSubsetFixture, 'subsetLabel', selectedRep.msFlow.selectedSubsetLabel));
selectedSubsetOffsets = reshape(localGetFieldOrDefault(selectedSubsetFixture, 'subsetOffsetIdx', []), 1, []);

truthDoaDeg = truth.latlonTrueDeg(:);
ssStaticDoaDeg = ssStaticCase.estResult.doaParamEst(:);
msStaticDoaDeg = msStaticCase.estResult.doaParamEst(:);
subsetDoaDeg = subsetCase.estResult.doaParamEst(:);
periodicDoaDeg = periodicCase.estResult.doaParamEst(:);
[latGrid, lonGrid] = localBuildScanGrid(truthDoaDeg, ssStaticDoaDeg, msStaticDoaDeg, subsetDoaDeg, periodicDoaDeg);

ssStaticModel = localBuildSfProbeModel( ...
  selectedRep.repeatData.periodicFixture.viewRefOnly, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, ...
  selectedRep.repeatData.periodicFixture.fdRange, selectedRep.staticBundle.staticRefOpt);
msStaticModel = localBuildSfProbeModel( ...
  selectedRep.repeatData.periodicFixture.viewMs, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, ...
  selectedRep.repeatData.periodicFixture.fdRange, selectedRep.staticBundle.staticMsOpt);

probeOpt = struct();
probeOpt.debugEnable = false;
subsetModel = localBuildMfProbeModel(selectedSubsetFixture, subsetDoaDeg, flowOpt.dynBaseOpt, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate);
periodicModel = localBuildMfProbeModel(selectedRep.repeatData.periodicFixture, periodicDoaDeg, flowOpt.dynBaseOpt, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate);

ssStaticSurface = localScanSfDoaSurface( ...
  ssStaticModel, ssStaticCase.estResult.fdRefEst, latGrid, lonGrid, ...
  truthDoaDeg, ssStaticDoaDeg, "SS-SF-Static");
msStaticSurface = localScanSfDoaSurface( ...
  msStaticModel, msStaticCase.estResult.fdRefEst, latGrid, lonGrid, ...
  truthDoaDeg, msStaticDoaDeg, "MS-SF-Static");
subsetSurface = localScanMfDoaSurface( ...
  subsetModel, subsetCase.estResult.fdRefEst, subsetCase.estResult.fdRateEst, ...
  latGrid, lonGrid, truthDoaDeg, subsetDoaDeg, probeOpt, ...
  sprintf('MS-MF-Subset-CP (%s)', char(selectedSubsetLabel)));
periodicSurface = localScanMfDoaSurface( ...
  periodicModel, periodicCase.estResult.fdRefEst, periodicCase.estResult.fdRateEst, ...
  latGrid, lonGrid, truthDoaDeg, periodicDoaDeg, probeOpt, "MS-MF-Periodic-CP");

fprintf('Running regressionMfDoaProfileWithHealthyFd ...\n');
fprintf('  selected taskSeed              : %d\n', selectedRep.repeatData.taskSeed);
fprintf('  selected subset label          : %s\n', selectedSubsetLabel);
fprintf('  selected subset offsets        : %s\n', localFormatIntegerRow(selectedSubsetOffsets));
fprintf('  dyn angle err (deg)            : %.6f\n', periodicSummary.angleErrDeg);
fprintf('  dyn fdRef err (Hz)             : %.6f\n', periodicSummary.fdRefErrHz);
fprintf('  dyn fdRate err (Hz/s)          : %.6f\n', periodicSummary.fdRateErrHzPerSec);
localPrintSurfaceSummary(ssStaticSurface);
localPrintSurfaceSummary(msStaticSurface);
localPrintSurfaceSummary(subsetSurface);
localPrintSurfaceSummary(periodicSurface);
localPlotDoaMeshComparison(latGrid, lonGrid, [ssStaticSurface, msStaticSurface, subsetSurface, periodicSurface]);
fprintf('PASS: regressionMfDoaProfileWithHealthyFd\n');


function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function tracker = localCreateProgressTracker(titleText, totalCount)
tracker = struct();
tracker.queue = parallel.pool.DataQueue;
tracker.titleText = char(titleText);
tracker.totalCount = totalCount;

fprintf('%s\n', tracker.titleText);
progressbar('displaymode', 'replace');
progressbar('minimalupdateinterval', 0.2);
progressbar('reset', tracker.totalCount);
afterEach(tracker.queue, @(~) progressbar('advance'));
end

function localCloseProgressTracker(tracker)
if isempty(tracker)
  return;
end
pause(0.05);
progressbar('end');
end

function [latGrid, lonGrid] = localBuildScanGrid(varargin)
pointMat = nan(2, nargin);
for iArg = 1:nargin
  pointCur = reshape(varargin{iArg}, [], 1);
  pointMat(:, iArg) = pointCur(1:2);
end
latAll = pointMat(1, :);
lonAll = pointMat(2, :);
latSpan = max(latAll) - min(latAll);
lonSpan = max(lonAll) - min(lonAll);
latMargin = max(0.003, 1.5 * latSpan);
lonMargin = max(0.003, 1.5 * lonSpan);
latGrid = linspace(min(latAll) - latMargin, max(latAll) + latMargin, 25);
lonGrid = linspace(min(lonAll) - lonMargin, max(lonAll) + lonMargin, 25);
end

function model = localBuildSfProbeModel(view, pilotWave, carrierFreq, sampleRate, fdRange, modelOpt)
numSource = 1;
[model, ~, ~, ~] = buildDoaDopplerSfModel( ...
  view.sceneRef, view.rxSigSf, pilotWave, carrierFreq, sampleRate, ...
  view.doaGrid, fdRange, numSource, modelOpt);
end

function model = localBuildMfProbeModel(fixtureUse, initDoaDeg, dynBaseOpt, pilotWave, carrierFreq, sampleRate)
modelOpt = dynBaseOpt;
modelOpt.steeringRefFrameIdx = fixtureUse.sceneSeq.refFrameIdx;
modelOpt.fdRateMode = 'unknown';
modelOpt.initDoaParam = reshape(initDoaDeg, [], 1);
modelOpt.initDoaHalfWidth = [0.02; 0.02];
modelOpt.enableFdAliasUnwrap = true;
modelOpt.freezeDoa = false;
modelOpt.disableUnknownDoaReleaseFloor = true;
[model, ~, ~] = buildDoaDopplerMfModel( ...
  fixtureUse.sceneSeq, fixtureUse.rxSigCell, pilotWave, carrierFreq, sampleRate, ...
  fixtureUse.viewMs.doaGrid, fixtureUse.fdRange, fixtureUse.fdRateRange, modelOpt);
end

function surfaceInfo = localScanSfDoaSurface(model, fdRefFix, latGrid, lonGrid, truthDoaDeg, estDoaDeg, displayName)
objGrid = nan(numel(latGrid), numel(lonGrid));
tracker = localCreateProgressTracker(sprintf('Scanning %s mesh', char(displayName)), numel(latGrid));
cleanupTracker = onCleanup(@() localCloseProgressTracker(tracker));
parfor iLat = 1:numel(latGrid)
  objRow = nan(1, numel(lonGrid));
  for iLon = 1:numel(lonGrid)
    optVar = [latGrid(iLat); lonGrid(iLon); fdRefFix];
    objRow(iLon) = evalDoaDopplerSfProfileLike(model, optVar);
  end
  objGrid(iLat, :) = objRow;
  send(tracker.queue, 1);
end
clear cleanupTracker;

truthObj = evalDoaDopplerSfProfileLike(model, [truthDoaDeg(:); fdRefFix]);
estObj = evalDoaDopplerSfProfileLike(model, [estDoaDeg(:); fdRefFix]);
surfaceInfo = localBuildSurfaceInfo(displayName, objGrid, latGrid, lonGrid, truthDoaDeg, estDoaDeg, truthObj, estObj, ...
  sprintf('fdRef fixed = %.6f Hz', fdRefFix));
end

function surfaceInfo = localScanMfDoaSurface(model, fdRefFix, fdRateFix, latGrid, lonGrid, truthDoaDeg, estDoaDeg, probeOpt, displayName)
objGrid = nan(numel(latGrid), numel(lonGrid));
tracker = localCreateProgressTracker(sprintf('Scanning %s mesh', char(displayName)), numel(latGrid));
cleanupTracker = onCleanup(@() localCloseProgressTracker(tracker));
parfor iLat = 1:numel(latGrid)
  objRow = nan(1, numel(lonGrid));
  for iLon = 1:numel(lonGrid)
    pointCur = struct();
    pointCur.doaParam = [latGrid(iLat); lonGrid(iLon)];
    pointCur.fdRef = fdRefFix;
    pointCur.fdRate = fdRateFix;
    probeEval = evalDoaDopplerMfProbePoint(model, pointCur, probeOpt);
    objRow(iLon) = localGetFieldOrDefault(probeEval, 'obj', NaN);
  end
  objGrid(iLat, :) = objRow;
  send(tracker.queue, 1);
end
clear cleanupTracker;

truthPoint = struct('doaParam', truthDoaDeg(:), 'fdRef', fdRefFix, 'fdRate', fdRateFix);
estPoint = struct('doaParam', estDoaDeg(:), 'fdRef', fdRefFix, 'fdRate', fdRateFix);
truthEval = evalDoaDopplerMfProbePoint(model, truthPoint, probeOpt);
estEval = evalDoaDopplerMfProbePoint(model, estPoint, probeOpt);
truthObj = localGetFieldOrDefault(truthEval, 'obj', NaN);
estObj = localGetFieldOrDefault(estEval, 'obj', NaN);
surfaceInfo = localBuildSurfaceInfo(displayName, objGrid, latGrid, lonGrid, truthDoaDeg, estDoaDeg, truthObj, estObj, ...
  sprintf('fdRef fixed = %.6f Hz, fdRate fixed = %.6f Hz/s', fdRefFix, fdRateFix));
end

function surfaceInfo = localBuildSurfaceInfo(displayName, objGrid, latGrid, lonGrid, truthDoaDeg, estDoaDeg, truthObj, estObj, fixedLabel)
if any(~isfinite(objGrid), 'all')
  error('regressionMfDoaProfileWithHealthyFd:NonFiniteSurface', ...
    '%s surface contains non-finite values.', char(displayName));
end

[minObj, minLinearIdx] = min(objGrid(:));
[minLatIdx, minLonIdx] = ind2sub(size(objGrid), minLinearIdx);
if minLatIdx == 1 || minLatIdx == numel(latGrid) || minLonIdx == 1 || minLonIdx == numel(lonGrid)
  error('regressionMfDoaProfileWithHealthyFd:MinimumHitsBoundary', ...
    '%s surface minimum lies on the scan boundary. Enlarge the scan range.', char(displayName));
end

surfaceInfo = struct();
surfaceInfo.displayName = char(displayName);
surfaceInfo.fixedLabel = char(fixedLabel);
surfaceInfo.objGrid = objGrid;
surfaceInfo.truthDoaDeg = reshape(truthDoaDeg, [], 1);
surfaceInfo.estDoaDeg = reshape(estDoaDeg, [], 1);
surfaceInfo.truthObj = truthObj;
surfaceInfo.estObj = estObj;
surfaceInfo.gridMinObj = minObj;
surfaceInfo.gridMinDoaDeg = [latGrid(minLatIdx); lonGrid(minLonIdx)];
end

function localPrintSurfaceSummary(surfaceInfo)
fprintf('  %s\n', surfaceInfo.displayName);
fprintf('    %s\n', surfaceInfo.fixedLabel);
fprintf('    truth objective              : %.6f\n', surfaceInfo.truthObj);
fprintf('    estimate objective           : %.6f\n', surfaceInfo.estObj);
fprintf('    grid min objective           : %.6f\n', surfaceInfo.gridMinObj);
fprintf('    grid min DoA (deg)           : [%.6f, %.6f]\n', ...
  surfaceInfo.gridMinDoaDeg(1), surfaceInfo.gridMinDoaDeg(2));
fprintf('    truth DoA (deg)              : [%.6f, %.6f]\n', ...
  surfaceInfo.truthDoaDeg(1), surfaceInfo.truthDoaDeg(2));
fprintf('    estimate DoA (deg)           : [%.6f, %.6f]\n', ...
  surfaceInfo.estDoaDeg(1), surfaceInfo.estDoaDeg(2));
end

function localPlotDoaMeshComparison(latGrid, lonGrid, surfaceInfoList)
numSurface = numel(surfaceInfoList);
if numSurface ~= 4
  error('regressionMfDoaProfileWithHealthyFd:UnexpectedSurfaceCount', ...
    'Expected exactly 4 surfaces, but got %d.', numSurface);
end

[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
figure('Name', 'DoA mesh comparison near truth', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for iSurf = 1:numSurface
  surfaceInfo = surfaceInfoList(iSurf);
  nexttile;
  mesh(lonMesh, latMesh, surfaceInfo.objGrid);
  hold on;
  truthObj = localInterpSurfaceValue(lonGrid, latGrid, surfaceInfo.objGrid, ...
    surfaceInfo.truthDoaDeg(2), surfaceInfo.truthDoaDeg(1));
  estObj = localInterpSurfaceValue(lonGrid, latGrid, surfaceInfo.objGrid, ...
    surfaceInfo.estDoaDeg(2), surfaceInfo.estDoaDeg(1));
  plot3(surfaceInfo.truthDoaDeg(2), surfaceInfo.truthDoaDeg(1), truthObj, ...
    'wo', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
  plot3(surfaceInfo.estDoaDeg(2), surfaceInfo.estDoaDeg(1), estObj, ...
    'ks', 'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
  hold off;
  xlabel('Longitude (deg)');
  ylabel('Latitude (deg)');
  zlabel('Objective');
  title({surfaceInfo.displayName, surfaceInfo.fixedLabel}, 'Interpreter', 'none');
  view(45, 35);
  grid on;
  legend({'Objective mesh', 'Truth DoA', 'Estimate DoA'}, 'Location', 'best');
end
end

function value = localInterpSurfaceValue(lonGrid, latGrid, objGrid, lonDeg, latDeg)
[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
value = interp2(lonMesh, latMesh, objGrid, lonDeg, latDeg, 'linear');
if ~isfinite(value)
  [~, latIdx] = min(abs(latGrid - latDeg));
  [~, lonIdx] = min(abs(lonGrid - lonDeg));
  value = objGrid(latIdx, lonIdx);
end
end

function textOut = localFormatIntegerRow(valueVec)
valueVec = reshape(valueVec, 1, []);
if isempty(valueVec)
  textOut = '[]';
  return;
end
textOut = sprintf('%d ', valueVec);
textOut = ['[', strtrim(textOut), ']'];
end
