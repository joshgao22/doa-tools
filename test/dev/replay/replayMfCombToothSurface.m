% replayMfCombToothSurface
% Purpose: evaluate one representative seed with an fdRef comb scan and a
% DoA-fdRef objective surface. This replay visualizes the 1/T_f tooth
% structure and a local truth-centered DoA-fdRef coupling diagnostic.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; saveSnapshot=true saves only
% replayData via saveExpSnapshot.

clear; close all; clc;

%% Replay configuration
% The representative scenario, satellite pair, waveform, frame offsets, and
% default search boxes are defined by buildDynamicDualSatEciContext. This replay
% only controls the seed, SNR, estimator mode, and objective scan grids.

replayName = "replayMfCombToothSurface";

% Case controls.
snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % Fixed seed for the representative case.
saveSnapshot = true;                % true saves the lightweight replayData only.
optVerbose = false;                 % true enables compact estimator/flow trace.
modeTag = "unknown";                % "unknown" scans CP-U; use "known" for CP-K.

% Scan centers. "truth" shows the ideal tooth geometry; "finalEstimate" shows
% the tooth selected by the current flow. "staticSeed" is also supported.
centerNameList = ["truth", "finalEstimate"];

% fdRef comb line scan. scanHalfWidthHz=[] means numAliasSide * aliasStepHz,
% where aliasStepHz is resolved from the MF model time offsets, normally 1/T_f.
numAliasSide = 2;                   % Number of +/- alias teeth to include.
numFdGrid = 301;                    % Number of fdRef samples in each line scan.
scanHalfWidthHz = [];               % Manual fdRef half-width override in Hz.

% DoA-fdRef surface scan. The DoA axis is a 1-D cut along truth->finalEstimate,
% not a full latitude/longitude grid. Increase grid counts only for smoother plots.
surfaceFdGridCount = 121;           % Number of fdRef samples in the surface scan.
surfaceDoaGridCount = 61;           % Number of DoA-line samples in the surface scan.
surfaceDoaHalfWidthDeg = 0.004;     % DoA-line half-width around each center.

% Truth-centered DoA-fdRef coupling diagnostic. This local window is kept
% narrower than the comb surface so the fitted curvature reflects local coupling.
couplingFdGridCount = 101;          % Number of fdRef samples in the coupling scan.
couplingDoaGridCount = 101;         % Number of DoA-line samples in the coupling scan.
couplingFdHalfWidthHz = 300;        % Local fdRef half-width for coupling in Hz.
couplingDoaHalfWidthDeg = 0.006;    % Local DoA-line half-width for coupling in deg.

%% Build context and flow options

fprintf('Running %s ...\n', char(replayName));
fprintf('  seed                             : %d\n', baseSeed);
fprintf('  snr (dB)                         : %.2f\n', snrDb);
fprintf('  mode                             : %s\n', char(modeTag));
fprintf('  save snapshot                    : %d\n', saveSnapshot);

%% Run replay batch

fprintf('[%s] Build one fixed-seed replay case.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
[caseData, context, flowOpt] = localBuildReplayCase(baseSeed, snrDb, optVerbose);

% Build reusable center points once so every line/surface scan uses the same
% truth, static seed, and final-estimate parameter definitions.
centerRegistry = localBuildCenterRegistry(caseData.truth, caseData.staticSeedCase, caseData.msFlow.caseFinal);
centerNameList = reshape(string(centerNameList), [], 1);
if isempty(centerNameList)
  centerNameList = "truth";
end

% The probe model reuses the formal MF evaluator path. Only the probe point
% changes during scans; no signal or objective code is duplicated here.
fprintf('[%s] Build MF probe model for %s mode.\n', char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag));
model = localBuildProbeModel(caseData.periodicFixture, context, flowOpt, modeTag);
aliasStepHz = localResolveAliasStep(model);
scanHalfWidthHzResolved = scanHalfWidthHz;
if isempty(scanHalfWidthHzResolved)
  scanHalfWidthHzResolved = numAliasSide * aliasStepHz;
end

lineRowCell = cell(numel(centerNameList), 1);
surfaceRowCell = cell(numel(centerNameList), 1);
lineScanCell = cell(numel(centerNameList), 1);
surfaceScanCell = cell(numel(centerNameList), 1);
for iCenter = 1:numel(centerNameList)
  centerName = centerNameList(iCenter);
  basePoint = localResolveCenterPoint(centerRegistry, centerName);

  fprintf('[%s] Run fdRef comb line scan: center=%s.\n', char(datetime('now', 'Format', 'HH:mm:ss')), char(centerName));
  lineScan = localRunFdRefLineScan(model, basePoint, numFdGrid, scanHalfWidthHzResolved, aliasStepHz, char("fdRef line " + centerName));
  lineScan.centerName = centerName;
  lineScanCell{iCenter} = lineScan;

  fprintf('[%s] Run DoA-fdRef surface scan: center=%s.\n', char(datetime('now', 'Format', 'HH:mm:ss')), char(centerName));
  surfaceScan = localRunDoaFdRefSurfaceScan(model, basePoint, centerRegistry.truth, centerRegistry.finalEstimate, ...
    surfaceFdGridCount, surfaceDoaGridCount, surfaceDoaHalfWidthDeg, scanHalfWidthHzResolved, aliasStepHz, ...
    char("DoA-fdRef surface " + centerName));
  surfaceScan.centerName = centerName;
  surfaceScanCell{iCenter} = surfaceScan;

  lineRowCell{iCenter} = table(string(centerName), lineScan.aliasStepHz, lineScan.minDeltaFdRef, ...
    lineScan.minAliasIndex, lineScan.centerDeltaObj, lineScan.minObj, ...
    'VariableNames', {'centerName', 'aliasStepHz', 'minDeltaFdRefHz', 'minAliasIndex', 'centerDeltaObj', 'minObj'});
  surfaceRowCell{iCenter} = table(string(centerName), surfaceScan.minDoaOffsetDeg, surfaceScan.minFdOffsetHz, ...
    surfaceScan.minAliasIndex, surfaceScan.centerDeltaObj, ...
    'VariableNames', {'centerName', 'minDoaOffsetDeg', 'minFdOffsetHz', 'minAliasIndex', 'centerDeltaObj'});
end
lineTable = vertcat(lineRowCell{:});
surfaceTable = vertcat(surfaceRowCell{:});

fprintf('[%s] Run truth-centered DoA-fdRef coupling diagnostic.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
couplingScan = localRunDoaFdCouplingDiagnostic(model, centerRegistry.truth, centerRegistry.finalEstimate, ...
  couplingFdGridCount, couplingDoaGridCount, couplingDoaHalfWidthDeg, couplingFdHalfWidthHz, aliasStepHz, ...
  'truth-centered DoA-fdRef coupling');
couplingTable = table(string(couplingScan.centerName), couplingScan.rhoDoaFd, ...
  couplingScan.ridgeSlopeHzPerDeg, couplingScan.quadRidgeSlopeHzPerDeg, couplingScan.quadFitRmse, ...
  couplingScan.minDoaOffsetDeg, couplingScan.minFdOffsetHz, ...
  'VariableNames', {'centerName', 'rhoDoaFd', 'ridgeSlopeHzPerDeg', ...
  'quadRidgeSlopeHzPerDeg', 'quadFitRmse', 'minDoaOffsetDeg', 'minFdOffsetHz'});

%% Data storage

config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.saveSnapshot = saveSnapshot;
config.optVerbose = optVerbose;
config.modeTag = modeTag;
config.centerNameList = centerNameList;
config.numAliasSide = numAliasSide;
config.numFdGrid = numFdGrid;
config.scanHalfWidthHz = scanHalfWidthHz;
config.scanHalfWidthHzResolved = scanHalfWidthHzResolved;
config.surfaceFdGridCount = surfaceFdGridCount;
config.surfaceDoaGridCount = surfaceDoaGridCount;
config.surfaceDoaHalfWidthDeg = surfaceDoaHalfWidthDeg;
config.couplingFdGridCount = couplingFdGridCount;
config.couplingDoaGridCount = couplingDoaGridCount;
config.couplingFdHalfWidthHz = couplingFdHalfWidthHz;
config.couplingDoaHalfWidthDeg = couplingDoaHalfWidthDeg;

flowSummary = struct();
flowSummary.selectedSubsetLabel = string(localGetFieldOrDefault(caseData.msFlow, 'selectedSubsetLabel', "unknown"));
flowSummary.selectedSubsetTrusted = localGetFieldOrDefault(caseData.msFlow, 'selectedSubsetTrusted', false);
flowSummary.subsetSelectionReason = string(localGetFieldOrDefault(caseData.msFlow, 'subsetSelectionReason', "unknown"));
flowSummary.selectedFinalTag = "unknown";
if isfield(caseData.msFlow, 'periodicDoaSeed') && isfield(caseData.msFlow.periodicDoaSeed, 'selectedReplayTag')
  flowSummary.selectedFinalTag = string(caseData.msFlow.periodicDoaSeed.selectedReplayTag);
end
flowSummary.finalSummary = localGetFieldOrDefault(caseData.msFlow, 'finalSummary', []);

replayData = struct();
replayData.replayName = string(replayName);
replayData.config = config;
replayData.contextSummary = localBuildContextSummary(context, flowOpt, aliasStepHz);
replayData.flowSummary = flowSummary;
replayData.centerRegistry = centerRegistry;
replayData.lineTable = lineTable;
replayData.surfaceTable = surfaceTable;
replayData.couplingTable = couplingTable;
replayData.lineScanCell = lineScanCell;
replayData.surfaceScanCell = surfaceScanCell;
replayData.couplingScan = couplingScan;

if saveSnapshot
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
config = replayData.config;
lineTable = replayData.lineTable;
surfaceTable = replayData.surfaceTable;
couplingTable = localGetFieldOrDefault(replayData, 'couplingTable', table());
lineScanCell = replayData.lineScanCell;
surfaceScanCell = replayData.surfaceScanCell;
couplingScan = localGetFieldOrDefault(replayData, 'couplingScan', []);
flowSummary = replayData.flowSummary;
contextSummary = localGetFieldOrDefault(replayData, 'contextSummary', struct());
aliasStepHz = localGetFieldOrDefault(contextSummary, 'toothStepHz', NaN);
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz)) && ~isempty(lineScanCell)
  aliasStepHz = lineScanCell{1}.aliasStepHz;
end

fprintf('Running replayMfCombToothSurface ...\n');
fprintf('  seed                             : %d\n', config.baseSeed);
fprintf('  snr (dB)                         : %.2f\n', config.snrDb);
fprintf('  mode                             : %s\n', char(string(config.modeTag)));
fprintf('  alias step (Hz)                  : %.6f\n', aliasStepHz);
fprintf('  scan half width (Hz)             : %.6f\n', config.scanHalfWidthHzResolved);
fprintf('  selected subset                  : %s\n', char(string(flowSummary.selectedSubsetLabel)));
fprintf('  selected final tag               : %s\n', char(string(flowSummary.selectedFinalTag)));
disp('========== fdRef comb line summary ==========');
disp(lineTable);
disp('========== DoA-fdRef surface summary ==========');
disp(surfaceTable);
if ~isempty(couplingTable)
  disp('========== Truth-centered DoA-fdRef coupling summary ==========');
  disp(couplingTable);
end

localPlotReplay(lineScanCell, surfaceScanCell, couplingScan, aliasStepHz);

%% Local helpers

function [caseData, context, flowOpt] = localBuildReplayCase(baseSeed, snrDb, optVerbose)
% Build the shared fixture, generate one noisy repeat, and run the current
% subset-periodic flow to obtain the final estimate used as a scan center.
parallelOpt = struct( ...
  'enableSubsetEvalParfor', false, ...
  'minSubsetEvalParfor', 4, ...
  'enableFixtureBankParfor', false, ...
  'enableParfor', false);
context = buildDynamicDualSatEciContext(struct('baseSeed', baseSeed, 'parallelOpt', parallelOpt));

% Signal construction is inside buildDynamicRepeatData so this replay stays tied
% to the same fixture/signal path used by other dynamic dev entries.
repeatData = buildDynamicRepeatData(context, snrDb, baseSeed);

% Tight periodic refine settings make this replay focus on tooth/surface shape
% rather than changing the default flow search strategy.
flowOpt = buildSimpleDynamicFlowOpt(struct( ...
  'parallelOpt', parallelOpt, ...
  'periodicRefineFdHalfWidthHz', 50, ...
  'periodicRefineFdRateHalfWidthHzPerSec', 100, ...
  'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'periodicRefineFreezeDoa', true, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti", ...
  'periodicPolishEnableWhenMulti', true, ...
  'periodicPolishDoaHalfWidthDeg', [0.002; 0.002]));

% Static transition provides the MS-SF seed used before subset tooth selection
% and periodic in-tooth refinement.
staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  repeatData.periodicFixture.viewRefOnly, repeatData.periodicFixture.viewOtherOnly, repeatData.periodicFixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  repeatData.periodicFixture.fdRange, repeatData.truth, context.otherSatIdxGlobal, ...
  optVerbose, flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);

% Run the same high-level flow being debugged; the scans below only probe the
% final objective around selected centers.
msFlow = runSimpleDynamicSubsetPeriodicFlow("MS-MF-Dyn-CombSurfaceReplay", "multi", ...
  repeatData.periodicFixture, repeatData.subsetFixtureCell, staticBundle.caseStaticMs, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, optVerbose, flowOpt);

caseData = struct();
caseData.periodicFixture = repeatData.periodicFixture;
caseData.truth = repeatData.truth;
caseData.staticSeedCase = staticBundle.caseStaticMs;
caseData.msFlow = msFlow;
end

function model = localBuildProbeModel(periodicFixture, context, flowOpt, modeTag)
% Build an MF evaluator model for the selected known-rate or unknown-rate mode.
modelOpt = flowOpt.dynBaseOpt;
modelOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
modelOpt.enableFdAliasUnwrap = true;
switch lower(char(modeTag))
  case {'known', 'cp-k', 'k'}
    modelOpt.fdRateMode = 'known';
    modelOpt.fdRateKnown = periodicFixture.truth.fdRateFit;
    fdRateRangeUse = [];
  case {'unknown', 'cp-u', 'u'}
    modelOpt.fdRateMode = 'unknown';
    fdRateRangeUse = periodicFixture.fdRateRange;
  otherwise
    error('replayMfCombToothSurface:InvalidMode', 'Unsupported modeTag: %s', char(string(modeTag)));
end
[model, ~, ~, ~] = buildDoaDopplerMfModel( ...
  periodicFixture.viewMs.sceneSeq, periodicFixture.viewMs.rxSigMf, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  periodicFixture.viewMs.doaGrid, periodicFixture.fdRange, fdRateRangeUse, modelOpt);
end

function centerRegistry = localBuildCenterRegistry(truth, staticSeedCase, finalCase)
% Collect truth, static-seed, and final-estimate center points for probing.
centerRegistry = struct();
centerRegistry.truth = struct('latlon', reshape(truth.latlonTrueDeg(1:2), [], 1), ...
  'fdRef', truth.fdRefFit, 'fdRate', truth.fdRateFit);
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
centerRegistry.staticSeed = struct( ...
  'latlon', localResolveDoaParam(staticEst, truth.latlonTrueDeg(:)), ...
  'fdRef', localResolveScalarField(staticEst, 'fdRefEst', truth.fdRefFit), ...
  'fdRate', truth.fdRateFit);
finalEst = localGetFieldOrDefault(finalCase, 'estResult', struct());
centerRegistry.finalEstimate = struct( ...
  'latlon', localResolveDoaParam(finalEst, centerRegistry.staticSeed.latlon(:)), ...
  'fdRef', localResolveScalarField(finalEst, 'fdRefEst', centerRegistry.staticSeed.fdRef), ...
  'fdRate', localResolveScalarField(finalEst, 'fdRateEst', truth.fdRateFit));
end

function point = localResolveCenterPoint(centerRegistry, centerName)
% Resolve one named scan center from the prebuilt center registry.
switch lower(char(centerName))
  case 'truth'
    point = centerRegistry.truth;
  case 'staticseed'
    point = centerRegistry.staticSeed;
  case 'finalestimate'
    point = centerRegistry.finalEstimate;
  otherwise
    error('replayMfCombToothSurface:UnknownCenter', 'Unsupported centerName: %s', centerName);
end
end

function lineScan = localRunFdRefLineScan(model, basePoint, numFdGrid, scanHalfWidthHz, aliasStepHz, progressTitle)
% Hold DoA and fdRate fixed, then scan fdRef across several alias teeth.
baseLatlon = basePoint.latlon(:);
baseFdRef = basePoint.fdRef;
baseFdRate = basePoint.fdRate;
fdRefGrid = linspace(baseFdRef - scanHalfWidthHz, baseFdRef + scanHalfWidthHz, numFdGrid).';
fdRefEvalVec = fdRefGrid;
numGrid = numel(fdRefEvalVec);
objVec = nan(numGrid, 1);
useParfor = localUseParfor(numGrid);
tracker = localCreateProgressTracker(progressTitle, numGrid, useParfor);
if useParfor
  progressQueue = tracker.queue;
  parfor iGrid = 1:numGrid
    fdRef = fdRefEvalVec(iGrid);
    objVec(iGrid) = localEvalObjAtPoint(model, baseLatlon, fdRef, baseFdRate);
    if ~isempty(progressQueue)
      send(progressQueue, iGrid);
    end
  end
else
  for iGrid = 1:numGrid
    fdRef = fdRefEvalVec(iGrid);
    objVec(iGrid) = localEvalObjAtPoint(model, baseLatlon, fdRef, baseFdRate);
    localAdvanceProgressTracker(tracker);
  end
end
localCloseProgressTracker(tracker);
[minObj, minIdx] = min(objVec);
[~, centerIdx] = min(abs(fdRefGrid - baseFdRef));
lineScan = struct();
lineScan.basePoint = basePoint;
lineScan.aliasStepHz = aliasStepHz;
lineScan.fdRefGrid = fdRefGrid;
lineScan.deltaFdRef = fdRefGrid - baseFdRef;
lineScan.objVec = objVec;
lineScan.deltaObjVec = objVec - minObj;
lineScan.foldedOffset = localFoldOffset(lineScan.deltaFdRef, aliasStepHz);
lineScan.minObj = minObj;
lineScan.minIdx = minIdx;
lineScan.minDeltaFdRef = fdRefGrid(minIdx) - baseFdRef;
lineScan.minAliasIndex = round(lineScan.minDeltaFdRef / aliasStepHz);
lineScan.centerDeltaObj = objVec(centerIdx) - minObj;
end

function couplingScan = localRunDoaFdCouplingDiagnostic(model, truthPoint, finalPoint, ...
  couplingFdGridCount, couplingDoaGridCount, couplingDoaHalfWidthDeg, couplingFdHalfWidthHz, aliasStepHz, progressTitle)
% Run a local truth-centered DoA-fdRef surface and attach coupling metrics.
couplingScan = localRunDoaFdRefSurfaceScan(model, truthPoint, truthPoint, finalPoint, ...
  couplingFdGridCount, couplingDoaGridCount, couplingDoaHalfWidthDeg, couplingFdHalfWidthHz, aliasStepHz, progressTitle);
couplingScan.centerName = "truthCoupling";
couplingScan.localWindow = struct('doaHalfWidthDeg', couplingDoaHalfWidthDeg, ...
  'fdHalfWidthHz', couplingFdHalfWidthHz, 'doaGridCount', couplingDoaGridCount, ...
  'fdGridCount', couplingFdGridCount);
couplingScan = localAttachCouplingMetrics(couplingScan);
end

function surfaceScan = localRunDoaFdRefSurfaceScan(model, basePoint, truthPoint, finalPoint, ...
  surfaceFdGridCount, surfaceDoaGridCount, surfaceDoaHalfWidthDeg, scanHalfWidthHz, aliasStepHz, progressTitle)
% Scan fdRef jointly with a 1-D DoA cut. The DoA cut direction follows the
% observed truth-to-final drift so the surface highlights the active bad basin.
baseLatlon = basePoint.latlon(:);
baseFdRef = basePoint.fdRef;
baseFdRate = basePoint.fdRate;
fdOffsetGrid = linspace(-scanHalfWidthHz, scanHalfWidthHz, surfaceFdGridCount).';
doaOffsetGrid = linspace(-surfaceDoaHalfWidthDeg, surfaceDoaHalfWidthDeg, surfaceDoaGridCount).';
doaDir = localResolveDoaDirection(truthPoint.latlon, finalPoint.latlon);
doaDir = doaDir(:);
numDoa = numel(doaOffsetGrid);
numFd = numel(fdOffsetGrid);
numPoint = numDoa * numFd;

% Use sliced per-point vectors inside parfor. Indexing short grid vectors through
% ind2sub makes MATLAB classify them as broadcast variables and triggers needless
% communication warnings.
doaOffsetEvalVec = repmat(doaOffsetGrid, numFd, 1);
fdOffsetEvalVec = repelem(fdOffsetGrid, numDoa);
doaLatEvalVec = baseLatlon(1) + doaDir(1) * doaOffsetEvalVec;
doaLonEvalVec = baseLatlon(2) + doaDir(2) * doaOffsetEvalVec;
fdRefEvalVec = baseFdRef + fdOffsetEvalVec;

objVec = nan(numPoint, 1);
useParfor = localUseParfor(numPoint);
tracker = localCreateProgressTracker(progressTitle, numPoint, useParfor);
if useParfor
  progressQueue = tracker.queue;
  parfor iPoint = 1:numPoint
    doaParam = [doaLatEvalVec(iPoint); doaLonEvalVec(iPoint)];
    fdRef = fdRefEvalVec(iPoint);
    objVec(iPoint) = localEvalObjAtPoint(model, doaParam, fdRef, baseFdRate);
    if ~isempty(progressQueue)
      send(progressQueue, iPoint);
    end
  end
else
  for iPoint = 1:numPoint
    doaParam = [doaLatEvalVec(iPoint); doaLonEvalVec(iPoint)];
    fdRef = fdRefEvalVec(iPoint);
    objVec(iPoint) = localEvalObjAtPoint(model, doaParam, fdRef, baseFdRate);
    localAdvanceProgressTracker(tracker);
  end
end
localCloseProgressTracker(tracker);
objMat = reshape(objVec, numDoa, numFd);
[minObj, minLinearIdx] = min(objMat(:));
[minDoaIdx, minFdIdx] = ind2sub(size(objMat), minLinearIdx);
[~, centerDoaIdx] = min(abs(doaOffsetGrid));
[~, centerFdIdx] = min(abs(fdOffsetGrid));
surfaceScan = struct();
surfaceScan.basePoint = basePoint;
surfaceScan.doaDir = doaDir;
surfaceScan.doaOffsetGridDeg = doaOffsetGrid;
surfaceScan.fdOffsetGridHz = fdOffsetGrid;
surfaceScan.deltaObjMat = objMat - minObj;
surfaceScan.minObj = minObj;
surfaceScan.minDoaOffsetDeg = doaOffsetGrid(minDoaIdx);
surfaceScan.minFdOffsetHz = fdOffsetGrid(minFdIdx);
surfaceScan.minAliasIndex = round(surfaceScan.minFdOffsetHz / aliasStepHz);
surfaceScan.centerDeltaObj = objMat(centerDoaIdx, centerFdIdx) - minObj;
end

function surfaceScan = localAttachCouplingMetrics(surfaceScan)
% Fit local quadratic curvature and ridge slope metrics to one surface scan.
[doaGridMat, fdGridMat] = ndgrid(surfaceScan.doaOffsetGridDeg(:), surfaceScan.fdOffsetGridHz(:));
zMat = surfaceScan.deltaObjMat;
fitInfo = localFitQuadraticSurface(doaGridMat(:), fdGridMat(:), zMat(:));
surfaceScan.quadCoeff = fitInfo.coeff;
surfaceScan.quadHessian = fitInfo.hessian;
surfaceScan.rhoDoaFd = fitInfo.rhoDoaFd;
surfaceScan.quadFitRmse = fitInfo.rmse;
surfaceScan.quadRidgeSlopeHzPerDeg = fitInfo.quadRidgeSlopeHzPerDeg;
[~, bestFdIdx] = min(zMat, [], 2);
ridgeFdBestHz = nan(numel(surfaceScan.doaOffsetGridDeg), 1);
fdGrid = surfaceScan.fdOffsetGridHz(:);
for iDoa = 1:numel(ridgeFdBestHz)
  if isfinite(bestFdIdx(iDoa)) && bestFdIdx(iDoa) >= 1 && bestFdIdx(iDoa) <= numel(fdGrid)
    ridgeFdBestHz(iDoa) = fdGrid(bestFdIdx(iDoa));
  end
end
validRidge = isfinite(surfaceScan.doaOffsetGridDeg(:)) & isfinite(ridgeFdBestHz);
if nnz(validRidge) >= 2
  ridgeCoeff = polyfit(surfaceScan.doaOffsetGridDeg(validRidge), ridgeFdBestHz(validRidge), 1);
else
  ridgeCoeff = [NaN, NaN];
end
surfaceScan.ridgeFdBestHz = ridgeFdBestHz;
surfaceScan.ridgeSlopeHzPerDeg = ridgeCoeff(1);
surfaceScan.ridgeInterceptHz = ridgeCoeff(2);
end

function fitInfo = localFitQuadraticSurface(doaOffsetDeg, fdOffsetHz, deltaObj)
% Fit z = c0 + c1*x + c2*y + c3*x^2 + c4*x*y + c5*y^2 and report curvature coupling.
x = doaOffsetDeg(:);
y = fdOffsetHz(:);
z = deltaObj(:);
valid = isfinite(x) & isfinite(y) & isfinite(z);
x = x(valid);
y = y(valid);
z = z(valid);
fitInfo = struct('coeff', nan(6, 1), 'hessian', nan(2, 2), 'rhoDoaFd', NaN, ...
  'rmse', NaN, 'quadRidgeSlopeHzPerDeg', NaN);
if numel(z) < 6
  return;
end
X = [ones(size(x)), x, y, x.^2, x .* y, y.^2];
coeff = X \ z;
zFit = X * coeff;
hessian = [2 * coeff(4), coeff(5); coeff(5), 2 * coeff(6)];
diagProd = hessian(1, 1) * hessian(2, 2);
rhoDoaFd = NaN;
if isfinite(diagProd) && diagProd > 0
  rhoDoaFd = hessian(1, 2) / sqrt(diagProd);
end
quadRidgeSlopeHzPerDeg = NaN;
if isfinite(hessian(1, 2)) && isfinite(hessian(2, 2)) && abs(hessian(2, 2)) > eps
  quadRidgeSlopeHzPerDeg = -hessian(1, 2) / hessian(2, 2);
end
fitInfo.coeff = coeff;
fitInfo.hessian = hessian;
fitInfo.rhoDoaFd = rhoDoaFd;
fitInfo.rmse = sqrt(mean((z - zFit).^2));
fitInfo.quadRidgeSlopeHzPerDeg = quadRidgeSlopeHzPerDeg;
end

function objVal = localEvalObjAtPoint(model, doaParam, fdRef, fdRate)
% Evaluate the formal MF probe objective at one DoA, fdRef, and fdRate point.
pointCur = struct('doaParam', doaParam(:), 'fdRef', fdRef, 'fdRate', fdRate);
probeTemplate = struct('doaParam', doaParam(:), 'fdRef', fdRef, 'fdRate', fdRate, ...
  'obj', NaN, 'residualNorm', NaN);
probeEval = evalDoaDopplerMfProbePoint(model, pointCur, struct('probeTemplate', probeTemplate, 'debugEnable', false));
objVal = localGetFieldOrDefault(probeEval, 'obj', NaN);
end

function localPlotReplay(lineScanCell, surfaceScanCell, couplingScan, aliasStepHz)
% Plot line, folded-comb, heatmap, mesh, and local coupling diagnostic views.
for iCase = 1:numel(lineScanCell)
  lineScan = lineScanCell{iCase};
  surfaceScan = surfaceScanCell{iCase};
  centerName = char(lineScan.centerName);

  figure('Name', sprintf('comb and DoA-fdRef surface: %s', centerName));

  subplot(2, 2, 1);
  plot(lineScan.deltaFdRef, lineScan.deltaObjVec, 'LineWidth', 1.2);
  hold on;
  localOverlayAliasLines(aliasStepHz, lineScan.deltaFdRef);
  xline(0, ':');
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('\Delta objective');
  title(sprintf('fdRef line | %s', centerName), 'Interpreter', 'none');

  subplot(2, 2, 2);
  scatter(lineScan.foldedOffset, lineScan.deltaObjVec, 16, 'filled');
  grid on;
  xlabel(sprintf('wrapped \\Delta fdRef mod %.3f Hz', aliasStepHz));
  ylabel('\Delta objective');
  title('folded comb', 'Interpreter', 'none');

  subplot(2, 2, 3);
  imagesc(surfaceScan.fdOffsetGridHz, surfaceScan.doaOffsetGridDeg, surfaceScan.deltaObjMat);
  set(gca, 'YDir', 'normal');
  colorbar;
  hold on;
  xline(0, ':');
  yline(0, ':');
  localOverlayAliasLines(aliasStepHz, surfaceScan.fdOffsetGridHz);
  xlabel('\Delta fdRef (Hz)');
  ylabel('DoA-line offset (deg)');
  title('DoA-fdRef heatmap', 'Interpreter', 'none');

  subplot(2, 2, 4);
  mesh(surfaceScan.fdOffsetGridHz, surfaceScan.doaOffsetGridDeg, surfaceScan.deltaObjMat);
  grid on;
  view(38, 28);
  xlabel('\Delta fdRef (Hz)');
  ylabel('DoA-line offset (deg)');
  zlabel('\Delta objective');
  title('DoA-fdRef mesh', 'Interpreter', 'none');
end
if nargin >= 3 && isstruct(couplingScan) && ~isempty(couplingScan)
  localPlotCouplingDiagnostic(couplingScan, aliasStepHz);
end
end

function localPlotCouplingDiagnostic(couplingScan, aliasStepHz)
% Plot the truth-centered local surface and the best-fd ridge used for coupling.
figure('Name', 'truth-centered DoA-fdRef coupling');
imagesc(couplingScan.fdOffsetGridHz, couplingScan.doaOffsetGridDeg, couplingScan.deltaObjMat);
set(gca, 'YDir', 'normal');
colorbar;
hold on;
plot(couplingScan.ridgeFdBestHz, couplingScan.doaOffsetGridDeg, '.-', 'LineWidth', 1.2);
xline(0, ':');
yline(0, ':');
localOverlayAliasLines(aliasStepHz, couplingScan.fdOffsetGridHz);
grid on;
xlabel('\Delta fdRef (Hz)');
ylabel('DoA-line offset (deg)');
title(sprintf('truth local coupling: rho=%.3g, ridge=%.3g Hz/deg', ...
  couplingScan.rhoDoaFd, couplingScan.ridgeSlopeHzPerDeg), 'Interpreter', 'none');
end

function contextSummary = localBuildContextSummary(context, flowOpt, aliasStepHz)
% Store lightweight context metadata needed to interpret the replay snapshot.
contextSummary = struct();
contextSummary.usrLla = localGetFieldOrDefault(context, 'usrLla', []);
contextSummary.utcRef = string(localGetFieldOrDefault(context, 'utcRef', ""));
contextSummary.tleFileName = string(localGetFieldOrDefault(context, 'tleFileName', ""));
contextSummary.selectedSatIdxGlobal = localGetFieldOrDefault(context, 'selectedSatIdxGlobal', []);
contextSummary.refSatIdxGlobal = localGetFieldOrDefault(context, 'refSatIdxGlobal', []);
contextSummary.satElevationMaskDeg = localGetFieldOrDefault(context, 'satElevationMaskDeg', []);
contextSummary.arraySize = localGetFieldOrDefault(context, 'arraySize', []);
contextSummary.elemSpacingWavelength = localGetFieldOrDefault(context, 'elemSpacingWavelength', NaN);
contextSummary.satPickCount = localGetFieldOrDefault(context, 'satPickCount', NaN);
contextSummary.satPickAnchorIdx = localGetFieldOrDefault(context, 'satPickAnchorIdx', NaN);
contextSummary.satPickMode = localGetFieldOrDefault(context, 'satPickMode', "");
contextSummary.frameIntvlSec = localGetFieldOrDefault(context, 'frameIntvlSec', NaN);
contextSummary.sampleRate = localGetFieldOrDefault(context, 'sampleRate', NaN);
contextSummary.symbolRate = localGetFieldOrDefault(context, 'symbolRate', NaN);
contextSummary.numSym = localGetFieldOrDefault(context, 'numSym', NaN);
contextSummary.carrierFreq = localGetFieldOrDefault(context, 'carrierFreq', NaN);
contextSummary.gridSize = localGetFieldOrDefault(context, 'gridSize', []);
contextSummary.fdRangeDefault = localGetFieldOrDefault(context, 'fdRangeDefault', []);
contextSummary.fdRateRangeDefault = localGetFieldOrDefault(context, 'fdRateRangeDefault', []);
contextSummary.toothStepHz = aliasStepHz;
contextSummary.periodicOffsetIdx = localGetFieldOrDefault(context, 'periodicOffsetIdx', []);
contextSummary.masterOffsetIdx = localGetFieldOrDefault(context, 'masterOffsetIdx', []);
contextSummary.subsetDefaultBankLabelList = localGetFieldOrDefault(flowOpt, 'subsetDefaultBankLabelList', strings(0, 1));
contextSummary.subsetMarginFallbackBankLabelList = localGetFieldOrDefault(flowOpt, 'subsetMarginFallbackBankLabelList', strings(0, 1));
end

function aliasStepHz = localResolveAliasStep(model)
% Resolve the fdRef tooth spacing from the model or from its frame times.
aliasStepHz = localGetFieldOrDefault(model, 'fdAliasStepHz', NaN);
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  dt = diff(model.timeOffsetSec(:));
  dt = dt(isfinite(dt) & dt > 0);
  if isempty(dt)
    aliasStepHz = NaN;
  else
    aliasStepHz = 1 / median(dt);
  end
end
end

function foldedOffset = localFoldOffset(deltaFdRef, aliasStepHz)
% Fold fdRef offsets into one alias period for compact comb visualization.
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  foldedOffset = deltaFdRef(:);
  return;
end
foldedOffset = mod(deltaFdRef(:) + aliasStepHz / 2, aliasStepHz) - aliasStepHz / 2;
end

function doaDir = localResolveDoaDirection(truthDoa, finalDoa)
% Build the 1-D DoA scan direction from truth toward the final estimate.
doaDiff = reshape(finalDoa, [], 1) - reshape(truthDoa, [], 1);
if numel(doaDiff) < 2 || norm(doaDiff(1:2)) < 1e-10
  doaDir = [1; 0];
else
  doaDir = doaDiff(1:2) / norm(doaDiff(1:2));
end
end

function doaParam = localResolveDoaParam(estResult, fallbackDoa)
% Extract a two-element DoA estimate or fall back to the provided truth/seed DoA.
doaParam = reshape(fallbackDoa, [], 1);
cand = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(cand) >= 2 && all(isfinite(cand(1:2)))
  doaParam = cand(1:2);
end
end

function scalarValue = localResolveScalarField(estResult, fieldName, fallbackValue)
% Extract one scalar estimate field or keep a deterministic fallback value.
scalarValue = fallbackValue;
cand = localGetFieldOrDefault(estResult, fieldName, []);
if isscalar(cand) && isnumeric(cand) && isfinite(cand)
  scalarValue = cand;
end
end

function localOverlayAliasLines(aliasStepHz, xGrid)
% Overlay alias-tooth reference lines inside the current x-axis grid span.
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz) && aliasStepHz > 0)
  return;
end
xMin = min(xGrid);
xMax = max(xGrid);
numSide = ceil(max(abs([xMin, xMax])) / aliasStepHz) + 1;
for k = -numSide:numSide
  xVal = k * aliasStepHz;
  if xVal < xMin || xVal > xMax
    continue;
  end
  xline(xVal, ':');
end
end

function useParfor = localUseParfor(numPoint)
% Decide whether this independent scan loop can use an outer parfor.
if numPoint <= 1
  useParfor = false;
  return;
end
try
  useParfor = isempty(getCurrentTask()) && ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
catch
  useParfor = false;
end
end

function tracker = localCreateProgressTracker(titleText, totalCount, useParfor)
% Create a progressbar tracker and DataQueue callback for serial or parfor scans.
tracker = struct('active', false, 'queue', []);
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
% Advance the progressbar for serial loops while leaving parfor updates to DataQueue.
if tracker.active && isempty(tracker.queue)
  progressbar('advance');
end
end

function localCloseProgressTracker(tracker)
% Close an active progressbar without touching replay results.
if tracker.active
  pause(0.05);
  progressbar('end');
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
% Return a nonempty struct field when present; otherwise return the default.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
