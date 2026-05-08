% replayMfCombToothSurface
% Purpose: evaluate one representative seed with an fdRef comb scan and a
% DoA-fdRef objective surface. This replay visualizes the 1/T_f tooth
% structure and a local truth-centered DoA-fdRef coupling diagnostic.
%
% Usage: edit the configuration block below, then run this script directly.
% The script leaves replayData in the workspace; saveSnapshot=true saves only
% replayData via saveExpSnapshot. Telegram notice is best-effort only.

clear; close all; clc;

%% Replay configuration
% The representative scenario, satellite pair, waveform, frame offsets, and
% default search boxes are defined by buildDynamicDualSatEciContext. This replay
% only controls the seed, SNR, estimator modes, and objective scan grids.

replayName = "replayMfCombToothSurface";

% Case controls.
snrDb = 10;                         % Snapshot SNR used to generate rx signals.
baseSeed = 253;                     % Fixed seed for the representative case.
saveSnapshot = true;                % true saves the lightweight replayData only.
notifyTelegramEnable = true;        % true sends best-effort HTML Telegram notice on completion or failure.
checkpointEnable = false;           % fixed-seed surface replay does not create per-repeat checkpoint files.
optVerbose = false;                 % true enables compact estimator/flow trace.
modeTagList = ["unknown", "known"];   % "unknown" scans CP-U; "known" scans CP-K.

% Scan centers. "truth" shows the ideal tooth geometry. "finalByMode" maps to
% finalCpU for CP-U scans and finalCpK for CP-K scans. "finalCpU",
% "finalCpK", "finalEstimate" (legacy CP-U alias), and "staticSeed" are also supported.
centerNameList = ["truth", "finalByMode"];

% fdRef comb line scan. scanHalfWidthHz=[] means numAliasSide * aliasStepHz,
% where aliasStepHz is resolved from the MF model time offsets, normally 1/T_f.
numAliasSide = 2;                   % Number of +/- alias teeth to include.
numFdGrid = 301;                    % Number of fdRef samples in each line scan.
scanHalfWidthHz = [];               % Manual fdRef half-width override in Hz.

% DoA-fdRef surface scan. The DoA axis is a 1-D cut from truth to the active center,
% not a full latitude/longitude grid. Increase grid counts only for smoother plots.
surfaceFdGridCount = 121;           % Number of fdRef samples in the surface scan.
surfaceDoaGridCount = 61;           % Number of DoA-line samples in the surface scan.
surfaceDoaHalfWidthDeg = 0.004;     % DoA-line half-width around each center.
surfaceClipDeltaObj = 1e4;          % Clipped heatmap upper limit; [] keeps the full color scale.

% Optional CP-K local latitude/longitude surface. This checks whether the
% known-rate DoA shape is hidden by the 1-D cut direction or fdRef tooth scale.
latLonSurfaceModeList = "known";      % Empty disables the extra lat/lon surface scan.
latLonSurfaceCenterNameList = "finalByMode";
latLonSurfaceGridCount = 41;         % Number of samples per latitude/longitude axis.
latLonSurfaceHalfWidthDeg = 0.004;   % Local lat/lon half-width around each center.
latLonSurfaceFdOffsetTagList = ["centerFd", "surfaceMinFd"];  % center fd and surface-selected fd tooth.

% Truth-centered DoA-fdRef coupling diagnostic. This local window is kept
% narrower than the comb surface so the fitted curvature reflects local coupling.
couplingFdGridCount = 101;          % Number of fdRef samples in the coupling scan.
couplingDoaGridCount = 101;         % Number of DoA-line samples in the coupling scan.
couplingFdHalfWidthHz = 300;        % Local fdRef half-width for coupling in Hz.
couplingDoaHalfWidthDeg = 0.006;    % Local DoA-line half-width for coupling in deg.

runTic = tic;
notifyConfig = struct( ...
  'snrDb', snrDb, ...
  'baseSeed', baseSeed, ...
  'modeTagList', modeTagList, ...
  'saveSnapshot', saveSnapshot, ...
  'notifyTelegramEnable', notifyTelegramEnable, ...
  'checkpointEnable', checkpointEnable);

%% Build context and flow options

printMfReplayHeader(char(replayName), notifyConfig, '');
fprintf('  %-32s : %s\n', 'mode list', char(strjoin(string(modeTagList), ', ')));

%% Run replay batch

try
fprintf('[%s] Build one fixed-seed replay case.\n', char(datetime('now', 'Format', 'HH:mm:ss')));
[caseData, context, flowOpt] = localBuildReplayCase(baseSeed, snrDb, optVerbose);

% Build reusable center points once so every line/surface scan uses the same
% truth, static seed, and final-estimate parameter definitions.
centerRegistry = localBuildCenterRegistry(caseData.truth, caseData.staticSeedCase, caseData.msFlow.caseFinal, caseData.msKnownCase);
centerTable = localBuildCenterTable(centerRegistry);
centerNameList = reshape(string(centerNameList), [], 1);
if isempty(centerNameList)
  centerNameList = "truth";
end

% The probe model reuses the formal MF evaluator path. Only the probe point
% changes during scans; no signal or objective code is duplicated here.
modeTagList = reshape(string(modeTagList), [], 1);
centerNameList = reshape(string(centerNameList), [], 1);
if isempty(modeTagList)
  modeTagList = "unknown";
end
if isempty(centerNameList)
  centerNameList = "truth";
end

numMode = numel(modeTagList);
numCenter = numel(centerNameList);
numScan = numMode * numCenter;
lineRowCell = cell(numScan, 1);
surfaceRowCell = cell(numScan, 1);
sliceRowCell = cell(numScan, 1);
lineScanCell = cell(numScan, 1);
surfaceScanCell = cell(numScan, 1);
latLonSurfaceRowCell = {};
latLonSurfaceCell = {};
latLonSurfaceIdx = 0;
couplingRowCell = cell(numMode, 1);
couplingScanCell = cell(numMode, 1);
modeRowCell = cell(numMode, 1);
scanHalfWidthHzResolvedList = nan(numMode, 1);
aliasStepHzList = nan(numMode, 1);
scanIdx = 0;

for iMode = 1:numMode
  modeTag = modeTagList(iMode);
  fprintf('[%s] Build MF probe model for %s mode.\n', char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag));
  model = localBuildProbeModel(caseData.periodicFixture, context, flowOpt, modeTag);
  aliasStepHz = localResolveAliasStep(model);
  scanHalfWidthHzResolved = scanHalfWidthHz;
  if isempty(scanHalfWidthHzResolved)
    scanHalfWidthHzResolved = numAliasSide * aliasStepHz;
  end
  aliasStepHzList(iMode) = aliasStepHz;
  scanHalfWidthHzResolvedList(iMode) = scanHalfWidthHzResolved;
  modeRowCell{iMode} = table(string(modeTag), aliasStepHz, scanHalfWidthHzResolved, ...
    'VariableNames', {'modeTag', 'aliasStepHz', 'scanHalfWidthHzResolved'});

  for iCenter = 1:numCenter
    centerName = centerNameList(iCenter);
    basePoint = localResolveCenterPoint(centerRegistry, centerName, modeTag);
    scanIdx = scanIdx + 1;

    fprintf('[%s] Run fdRef comb line scan: mode=%s, center=%s.\n', ...
      char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag), char(basePoint.resolvedCenterName));
    lineScan = localRunFdRefLineScan(model, basePoint, numFdGrid, scanHalfWidthHzResolved, aliasStepHz, ...
      char("fdRef line " + modeTag + " " + basePoint.resolvedCenterName));
    lineScan.modeTag = modeTag;
    lineScan.centerName = centerName;
    lineScan.resolvedCenterName = basePoint.resolvedCenterName;
    lineScanCell{scanIdx} = lineScan;

    fprintf('[%s] Run DoA-fdRef surface scan: mode=%s, center=%s.\n', ...
      char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag), char(basePoint.resolvedCenterName));
    surfaceScan = localRunDoaFdRefSurfaceScan(model, basePoint, centerRegistry.truth, basePoint, ...
      surfaceFdGridCount, surfaceDoaGridCount, surfaceDoaHalfWidthDeg, scanHalfWidthHzResolved, aliasStepHz, ...
      char("DoA-fdRef surface " + modeTag + " " + basePoint.resolvedCenterName));
    surfaceScan.modeTag = modeTag;
    surfaceScan.centerName = centerName;
    surfaceScan.resolvedCenterName = basePoint.resolvedCenterName;
    surfaceScan = localAttachSurfaceSliceMetrics(surfaceScan);
    surfaceScanCell{scanIdx} = surfaceScan;

    if localShouldRunLatLonSurface(modeTag, centerName, basePoint.resolvedCenterName, ...
        latLonSurfaceModeList, latLonSurfaceCenterNameList)
      fdOffsetTable = localResolveLatLonSurfaceFdOffsetTable(latLonSurfaceFdOffsetTagList, surfaceScan);
      for iFdOffset = 1:height(fdOffsetTable)
        fixedFdTag = fdOffsetTable.fixedFdTag(iFdOffset);
        fixedFdOffsetHz = fdOffsetTable.fixedFdOffsetHz(iFdOffset);
        fprintf('[%s] Run local lat/lon surface scan: mode=%s, center=%s, fdTag=%s, fdOffset=%.3f Hz.\n', ...
          char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag), char(basePoint.resolvedCenterName), ...
          char(fixedFdTag), fixedFdOffsetHz);
        latLonSurfaceIdx = latLonSurfaceIdx + 1;
        latLonSurface = localRunLatLonSurfaceScan(model, basePoint, fixedFdOffsetHz, ...
          latLonSurfaceGridCount, latLonSurfaceHalfWidthDeg, aliasStepHz, ...
          char("lat/lon surface " + modeTag + " " + basePoint.resolvedCenterName + " " + fixedFdTag));
        latLonSurface.modeTag = modeTag;
        latLonSurface.centerName = centerName;
        latLonSurface.resolvedCenterName = basePoint.resolvedCenterName;
        latLonSurface.fixedFdTag = fixedFdTag;
        latLonSurfaceCell{latLonSurfaceIdx, 1} = latLonSurface;
        latLonSurfaceRowCell{latLonSurfaceIdx, 1} = localBuildLatLonSurfaceRow(latLonSurface);
      end
    end

    lineRowCell{scanIdx} = table(string(modeTag), string(centerName), string(lineScan.resolvedCenterName), ...
      lineScan.aliasStepHz, lineScan.minDeltaFdRef, lineScan.minAliasIndex, lineScan.centerDeltaObj, lineScan.minObj, ...
      'VariableNames', {'modeTag', 'centerName', 'resolvedCenterName', 'aliasStepHz', ...
      'minDeltaFdRefHz', 'minAliasIndex', 'centerDeltaObj', 'minObj'});
    surfaceRowCell{scanIdx} = table(string(modeTag), string(centerName), string(surfaceScan.resolvedCenterName), ...
      surfaceScan.minDoaOffsetDeg, surfaceScan.minFdOffsetHz, surfaceScan.minAliasIndex, surfaceScan.centerDeltaObj, ...
      'VariableNames', {'modeTag', 'centerName', 'resolvedCenterName', ...
      'minDoaOffsetDeg', 'minFdOffsetHz', 'minAliasIndex', 'centerDeltaObj'});
    sliceRowCell{scanIdx} = localBuildSurfaceSliceRow(surfaceScan);
  end

  fprintf('[%s] Run truth-centered DoA-fdRef coupling diagnostic: mode=%s.\n', ...
    char(datetime('now', 'Format', 'HH:mm:ss')), char(modeTag));
  finalByModePoint = localResolveCenterPoint(centerRegistry, "finalByMode", modeTag);
  couplingScan = localRunDoaFdCouplingDiagnostic(model, centerRegistry.truth, finalByModePoint, ...
    couplingFdGridCount, couplingDoaGridCount, couplingDoaHalfWidthDeg, couplingFdHalfWidthHz, aliasStepHz, ...
    char("truth-centered DoA-fdRef coupling " + modeTag));
  couplingScan.modeTag = modeTag;
  couplingScanCell{iMode} = couplingScan;
  couplingRowCell{iMode} = table(string(modeTag), string(couplingScan.centerName), couplingScan.rhoDoaFd, ...
    couplingScan.ridgeSlopeHzPerDeg, couplingScan.quadRidgeSlopeHzPerDeg, couplingScan.quadFitRmse, ...
    couplingScan.minDoaOffsetDeg, couplingScan.minFdOffsetHz, ...
    'VariableNames', {'modeTag', 'centerName', 'rhoDoaFd', 'ridgeSlopeHzPerDeg', ...
    'quadRidgeSlopeHzPerDeg', 'quadFitRmse', 'minDoaOffsetDeg', 'minFdOffsetHz'});
end
lineTable = vertcat(lineRowCell{:});
surfaceTable = vertcat(surfaceRowCell{:});
sliceTable = vertcat(sliceRowCell{:});
if isempty(latLonSurfaceRowCell)
  latLonSurfaceTable = table();
else
  latLonSurfaceTable = vertcat(latLonSurfaceRowCell{:});
end
couplingTable = vertcat(couplingRowCell{:});
modeTable = vertcat(modeRowCell{:});
aliasStepHz = aliasStepHzList(1);
scanHalfWidthHzResolved = scanHalfWidthHzResolvedList(1);

%% Data storage

config = struct();
config.snrDb = snrDb;
config.baseSeed = baseSeed;
config.saveSnapshot = saveSnapshot;
config.notifyTelegramEnable = notifyTelegramEnable;
config.checkpointEnable = checkpointEnable;
config.optVerbose = optVerbose;
config.modeTagList = modeTagList;
config.centerNameList = centerNameList;
config.numAliasSide = numAliasSide;
config.numFdGrid = numFdGrid;
config.scanHalfWidthHz = scanHalfWidthHz;
config.scanHalfWidthHzResolved = scanHalfWidthHzResolved;
config.scanHalfWidthHzResolvedList = scanHalfWidthHzResolvedList;
config.surfaceFdGridCount = surfaceFdGridCount;
config.surfaceDoaGridCount = surfaceDoaGridCount;
config.surfaceDoaHalfWidthDeg = surfaceDoaHalfWidthDeg;
config.surfaceClipDeltaObj = surfaceClipDeltaObj;
config.latLonSurfaceModeList = latLonSurfaceModeList;
config.latLonSurfaceCenterNameList = latLonSurfaceCenterNameList;
config.latLonSurfaceGridCount = latLonSurfaceGridCount;
config.latLonSurfaceHalfWidthDeg = latLonSurfaceHalfWidthDeg;
config.latLonSurfaceFdOffsetTagList = latLonSurfaceFdOffsetTagList;
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
replayData.centerTable = centerTable;
replayData.modeTable = modeTable;
replayData.lineTable = lineTable;
replayData.surfaceTable = surfaceTable;
replayData.sliceTable = sliceTable;
replayData.latLonSurfaceTable = latLonSurfaceTable;
replayData.couplingTable = couplingTable;
replayData.lineScanCell = lineScanCell;
replayData.surfaceScanCell = surfaceScanCell;
replayData.latLonSurfaceCell = latLonSurfaceCell;
replayData.couplingScanCell = couplingScanCell;
replayData.couplingScan = couplingScanCell{1};

if saveSnapshot
  saveOpt = struct('includeVars', {{'replayData'}}, ...
    'extraMeta', struct('replayName', char(replayName)), 'verbose', true);
  replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
else
  replayData.snapshotFile = "";
end
replayData.elapsedSec = toc(runTic);
notifyMfReplayStatus(struct( ...
  'replayName', replayName, ...
  'statusText', "DONE", ...
  'config', config, ...
  'snapshotFile', replayData.snapshotFile, ...
  'elapsedSec', replayData.elapsedSec, ...
  'commentLineList', "Fixed-seed comb/surface diagnostic completed."));
catch ME
  notifyMfReplayStatus(struct( ...
    'replayName', replayName, ...
    'statusText', "FAILED", ...
    'config', notifyConfig, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Summary output and plotting

if ~exist('replayData', 'var') || ~isstruct(replayData)
  error('Replay data is missing. Run the replay batch sections or load a snapshot containing replayData.');
end
config = replayData.config;
lineTable = replayData.lineTable;
surfaceTable = replayData.surfaceTable;
sliceTable = localGetFieldOrDefault(replayData, 'sliceTable', table());
latLonSurfaceTable = localGetFieldOrDefault(replayData, 'latLonSurfaceTable', table());
couplingTable = localGetFieldOrDefault(replayData, 'couplingTable', table());
centerTable = localGetFieldOrDefault(replayData, 'centerTable', table());
lineScanCell = replayData.lineScanCell;
surfaceScanCell = replayData.surfaceScanCell;
latLonSurfaceCell = localGetFieldOrDefault(replayData, 'latLonSurfaceCell', {});
couplingScanCell = localGetFieldOrDefault(replayData, 'couplingScanCell', {});
if isempty(couplingScanCell)
  couplingScanCell = {localGetFieldOrDefault(replayData, 'couplingScan', [])};
end
flowSummary = replayData.flowSummary;
contextSummary = localGetFieldOrDefault(replayData, 'contextSummary', struct());
aliasStepHz = localGetFieldOrDefault(contextSummary, 'toothStepHz', NaN);
if ~(isscalar(aliasStepHz) && isfinite(aliasStepHz)) && ~isempty(lineScanCell)
  aliasStepHz = lineScanCell{1}.aliasStepHz;
end

fprintf('Running replayMfCombToothSurface ...\n');
fprintf('  seed                             : %d\n', config.baseSeed);
fprintf('  snr (dB)                         : %.2f\n', config.snrDb);
modeTagText = strjoin(string(localGetFieldOrDefault(config, 'modeTagList', string(localGetFieldOrDefault(config, 'modeTag', "unknown")))), ', ');
fprintf('  mode list                        : %s\n', char(modeTagText));
fprintf('  alias step (Hz)                  : %.6f\n', aliasStepHz);
fprintf('  scan half width (Hz)             : %.6f\n', config.scanHalfWidthHzResolved);
fprintf('  selected subset                  : %s\n', char(string(flowSummary.selectedSubsetLabel)));
fprintf('  selected final tag               : %s\n', char(string(flowSummary.selectedFinalTag)));
if ~isempty(centerTable)
  disp('========== center summary ==========');
  disp(centerTable);
end
modeTable = localGetFieldOrDefault(replayData, 'modeTable', table());
if ~isempty(modeTable)
  disp('========== mode summary ==========');
  disp(modeTable);
end
disp('========== fdRef comb line summary ==========');
disp(lineTable);
disp('========== DoA-fdRef surface summary ==========');
disp(surfaceTable);
if ~isempty(sliceTable)
  disp('========== DoA-fdRef slice-scale summary ==========');
  disp(sliceTable);
end
if ~isempty(latLonSurfaceTable)
  disp('========== local lat/lon surface summary ==========');
  disp(latLonSurfaceTable);
end
if ~isempty(couplingTable)
  disp('========== Truth-centered DoA-fdRef coupling summary ==========');
  disp(couplingTable);
end

localPlotReplay(lineScanCell, surfaceScanCell, couplingScanCell, aliasStepHz, ...
  localGetFieldOrDefault(config, 'surfaceClipDeltaObj', []), latLonSurfaceCell);

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
caseData.msKnownCase = localRunMsKnownFinalCase(repeatData.periodicFixture, staticBundle.caseStaticMs, context, flowOpt, optVerbose);
end

function caseKnown = localRunMsKnownFinalCase(periodicFixture, staticSeedCase, context, flowOpt, optVerbose)
% Run one replay-local MS-MF-CP-K solve to define the CP-K final scan center.
viewUse = periodicFixture.viewMs;
truthUse = periodicFixture.truth;
dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.enableFdAliasUnwrap = true;
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
seedDoa = localResolveDoaParam(staticEst, truthUse.latlonTrueDeg(:));
dynOpt.initDoaParam = seedDoa(:);
dynOpt.initDoaHalfWidth = flowOpt.multiDoaHalfWidth(:);
initParam = buildDynamicInitParamFromCase(staticSeedCase, true, 0);
initCandidate = struct('startTag', "fromStaticCpK", 'initParam', initParam, ...
  'initDoaParam', seedDoa(:), 'initDoaHalfWidth', flowOpt.multiDoaHalfWidth(:));
debugTruth = localGetFieldOrDefault(periodicFixture, 'debugTruthMs', struct());
caseKnown = runDynamicDoaDopplerCase("MS-MF-CP-K-CombSurfaceReplay", "multi", ...
  viewUse, truthUse, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynOpt, true, debugTruth, initCandidate);
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

function centerRegistry = localBuildCenterRegistry(truth, staticSeedCase, finalCpUCase, finalCpKCase)
% Collect truth, static-seed, CP-U final, and CP-K final center points.
centerRegistry = struct();
centerRegistry.truth = localAttachCenterName(struct('latlon', reshape(truth.latlonTrueDeg(1:2), [], 1), ...
  'fdRef', truth.fdRefFit, 'fdRate', truth.fdRateFit), "truth");
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
centerRegistry.staticSeed = localAttachCenterName(struct( ...
  'latlon', localResolveDoaParam(staticEst, truth.latlonTrueDeg(:)), ...
  'fdRef', localResolveScalarField(staticEst, 'fdRefEst', truth.fdRefFit), ...
  'fdRate', truth.fdRateFit), "staticSeed");
finalCpUEst = localGetFieldOrDefault(finalCpUCase, 'estResult', struct());
centerRegistry.finalCpU = localAttachCenterName(struct( ...
  'latlon', localResolveDoaParam(finalCpUEst, centerRegistry.staticSeed.latlon(:)), ...
  'fdRef', localResolveScalarField(finalCpUEst, 'fdRefEst', centerRegistry.staticSeed.fdRef), ...
  'fdRate', localResolveScalarField(finalCpUEst, 'fdRateEst', truth.fdRateFit)), "finalCpU");
finalCpKEst = localGetFieldOrDefault(finalCpKCase, 'estResult', struct());
centerRegistry.finalCpK = localAttachCenterName(struct( ...
  'latlon', localResolveDoaParam(finalCpKEst, centerRegistry.staticSeed.latlon(:)), ...
  'fdRef', localResolveScalarField(finalCpKEst, 'fdRefEst', centerRegistry.staticSeed.fdRef), ...
  'fdRate', truth.fdRateFit), "finalCpK");
centerRegistry.finalEstimate = centerRegistry.finalCpU;
centerRegistry.finalEstimate.resolvedCenterName = "finalEstimate";
end

function centerTable = localBuildCenterTable(centerRegistry)
% Build a compact center table for CP-U / CP-K surface interpretation.
nameList = ["truth"; "staticSeed"; "finalCpU"; "finalCpK"];
rowCell = cell(numel(nameList), 1);
for iName = 1:numel(nameList)
  nameCur = nameList(iName);
  point = centerRegistry.(char(nameCur));
  latlon = reshape(point.latlon, [], 1);
  rowCell{iName} = table(nameCur, latlon(1), latlon(2), point.fdRef, point.fdRate, ...
    'VariableNames', {'centerName', 'latDeg', 'lonDeg', 'fdRefHz', 'fdRateHzPerSec'});
end
centerTable = vertcat(rowCell{:});
end

function point = localResolveCenterPoint(centerRegistry, centerName, modeTag)
% Resolve one named scan center from the prebuilt center registry.
centerKey = lower(char(centerName));
modeKey = lower(char(modeTag));
switch centerKey
  case 'truth'
    point = centerRegistry.truth;
  case 'staticseed'
    point = centerRegistry.staticSeed;
  case {'finalcpu', 'finalestimate'}
    point = centerRegistry.finalCpU;
    if strcmp(centerKey, 'finalestimate')
      point.resolvedCenterName = "finalEstimate";
    end
  case 'finalcpk'
    point = centerRegistry.finalCpK;
  case 'finalbymode'
    if ismember(modeKey, {'known', 'cp-k', 'k'})
      point = centerRegistry.finalCpK;
    else
      point = centerRegistry.finalCpU;
    end
  otherwise
    error('replayMfCombToothSurface:UnknownCenter', 'Unsupported centerName: %s', centerName);
end
end

function point = localAttachCenterName(point, resolvedCenterName)
% Attach a resolved center label used in tables and plot titles.
point.resolvedCenterName = string(resolvedCenterName);
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

function surfaceScan = localAttachSurfaceSliceMetrics(surfaceScan)
% Attach scale diagnostics along the best and center fd/DoA slices.
zMat = surfaceScan.deltaObjMat;
[~, minLinearIdx] = min(zMat(:));
[minDoaIdx, minFdIdx] = ind2sub(size(zMat), minLinearIdx);
[~, centerDoaIdx] = min(abs(surfaceScan.doaOffsetGridDeg));
[~, centerFdIdx] = min(abs(surfaceScan.fdOffsetGridHz));
surfaceScan.minDoaIdx = minDoaIdx;
surfaceScan.minFdIdx = minFdIdx;
surfaceScan.centerDoaIdx = centerDoaIdx;
surfaceScan.centerFdIdx = centerFdIdx;
surfaceScan.bestFdDoaSlice = zMat(:, minFdIdx);
surfaceScan.centerFdDoaSlice = zMat(:, centerFdIdx);
surfaceScan.bestDoaFdSlice = zMat(minDoaIdx, :).';
surfaceScan.centerDoaFdSlice = zMat(centerDoaIdx, :).';
surfaceScan.bestFdDoaSpan = localFiniteSpan(surfaceScan.bestFdDoaSlice);
surfaceScan.centerFdDoaSpan = localFiniteSpan(surfaceScan.centerFdDoaSlice);
surfaceScan.bestDoaFdSpan = localFiniteSpan(surfaceScan.bestDoaFdSlice);
surfaceScan.centerDoaFdSpan = localFiniteSpan(surfaceScan.centerDoaFdSlice);
surfaceScan.bestSliceFdToDoaSpanRatio = localSafeRatio(surfaceScan.bestDoaFdSpan, surfaceScan.bestFdDoaSpan);
surfaceScan.centerSliceFdToDoaSpanRatio = localSafeRatio(surfaceScan.centerDoaFdSpan, surfaceScan.centerFdDoaSpan);
end

function row = localBuildSurfaceSliceRow(surfaceScan)
% Build a compact table row that exposes DoA-vs-fdRef scale dominance.
row = table(string(localGetFieldOrDefault(surfaceScan, 'modeTag', "unknown")), ...
  string(localGetFieldOrDefault(surfaceScan, 'centerName', "unknown")), ...
  string(localGetFieldOrDefault(surfaceScan, 'resolvedCenterName', "unknown")), ...
  surfaceScan.minDoaOffsetDeg, surfaceScan.minFdOffsetHz, ...
  surfaceScan.bestFdDoaSpan, surfaceScan.centerFdDoaSpan, ...
  surfaceScan.bestDoaFdSpan, surfaceScan.centerDoaFdSpan, ...
  surfaceScan.bestSliceFdToDoaSpanRatio, surfaceScan.centerSliceFdToDoaSpanRatio, ...
  'VariableNames', {'modeTag', 'centerName', 'resolvedCenterName', ...
  'minDoaOffsetDeg', 'minFdOffsetHz', 'bestFdDoaSpan', 'centerFdDoaSpan', ...
  'bestDoaFdSpan', 'centerDoaFdSpan', 'bestSliceFdToDoaSpanRatio', 'centerSliceFdToDoaSpanRatio'});
end

function latLonSurface = localRunLatLonSurfaceScan(model, basePoint, fdOffsetHz, ...
  latLonSurfaceGridCount, latLonSurfaceHalfWidthDeg, aliasStepHz, progressTitle)
% Scan a local latitude/longitude surface at one fixed fdRef tooth.
baseLatlon = basePoint.latlon(:);
baseFdRef = basePoint.fdRef;
baseFdRate = basePoint.fdRate;
latOffsetGrid = linspace(-latLonSurfaceHalfWidthDeg, latLonSurfaceHalfWidthDeg, latLonSurfaceGridCount).';
lonOffsetGrid = linspace(-latLonSurfaceHalfWidthDeg, latLonSurfaceHalfWidthDeg, latLonSurfaceGridCount).';
numLat = numel(latOffsetGrid);
numLon = numel(lonOffsetGrid);
numPoint = numLat * numLon;
latOffsetEvalVec = repmat(latOffsetGrid, numLon, 1);
lonOffsetEvalVec = repelem(lonOffsetGrid, numLat);
latEvalVec = baseLatlon(1) + latOffsetEvalVec;
lonEvalVec = baseLatlon(2) + lonOffsetEvalVec;
fdRefUse = baseFdRef + fdOffsetHz;
objVec = nan(numPoint, 1);
useParfor = localUseParfor(numPoint);
tracker = localCreateProgressTracker(progressTitle, numPoint, useParfor);
if useParfor
  progressQueue = tracker.queue;
  parfor iPoint = 1:numPoint
    doaParam = [latEvalVec(iPoint); lonEvalVec(iPoint)];
    objVec(iPoint) = localEvalObjAtPoint(model, doaParam, fdRefUse, baseFdRate);
    if ~isempty(progressQueue)
      send(progressQueue, iPoint);
    end
  end
else
  for iPoint = 1:numPoint
    doaParam = [latEvalVec(iPoint); lonEvalVec(iPoint)];
    objVec(iPoint) = localEvalObjAtPoint(model, doaParam, fdRefUse, baseFdRate);
    localAdvanceProgressTracker(tracker);
  end
end
localCloseProgressTracker(tracker);
objMat = reshape(objVec, numLat, numLon);
[minObj, minLinearIdx] = min(objMat(:));
[minLatIdx, minLonIdx] = ind2sub(size(objMat), minLinearIdx);
[~, centerLatIdx] = min(abs(latOffsetGrid));
[~, centerLonIdx] = min(abs(lonOffsetGrid));
deltaObjMat = objMat - minObj;
latLonSurface = struct();
latLonSurface.basePoint = basePoint;
latLonSurface.fixedFdOffsetHz = fdOffsetHz;
latLonSurface.fixedFdAliasIndex = round(fdOffsetHz / aliasStepHz);
latLonSurface.fixedFdRefHz = fdRefUse;
latLonSurface.latOffsetGridDeg = latOffsetGrid;
latLonSurface.lonOffsetGridDeg = lonOffsetGrid;
latLonSurface.deltaObjMat = deltaObjMat;
latLonSurface.minObj = minObj;
latLonSurface.minLatOffsetDeg = latOffsetGrid(minLatIdx);
latLonSurface.minLonOffsetDeg = lonOffsetGrid(minLonIdx);
latLonSurface.minDoaRadialOffsetDeg = hypot(latLonSurface.minLatOffsetDeg, latLonSurface.minLonOffsetDeg);
latLonSurface.centerDeltaObj = objMat(centerLatIdx, centerLonIdx) - minObj;
latLonSurface.bestLonLatSlice = deltaObjMat(:, minLonIdx);
latLonSurface.bestLatLonSlice = deltaObjMat(minLatIdx, :).';
latLonSurface.centerLonLatSlice = deltaObjMat(:, centerLonIdx);
latLonSurface.centerLatLonSlice = deltaObjMat(centerLatIdx, :).';
latLonSurface.bestLonLatSpan = localFiniteSpan(latLonSurface.bestLonLatSlice);
latLonSurface.bestLatLonSpan = localFiniteSpan(latLonSurface.bestLatLonSlice);
latLonSurface.centerLonLatSpan = localFiniteSpan(latLonSurface.centerLonLatSlice);
latLonSurface.centerLatLonSpan = localFiniteSpan(latLonSurface.centerLatLonSlice);
end

function row = localBuildLatLonSurfaceRow(latLonSurface)
% Build a compact table row for the optional fixed-fdRef lat/lon surface.
row = table(string(localGetFieldOrDefault(latLonSurface, 'modeTag', "unknown")), ...
  string(localGetFieldOrDefault(latLonSurface, 'centerName', "unknown")), ...
  string(localGetFieldOrDefault(latLonSurface, 'resolvedCenterName', "unknown")), ...
  string(localGetFieldOrDefault(latLonSurface, 'fixedFdTag', "unknown")), ...
  latLonSurface.fixedFdOffsetHz, latLonSurface.fixedFdAliasIndex, ...
  latLonSurface.minLatOffsetDeg, latLonSurface.minLonOffsetDeg, latLonSurface.minDoaRadialOffsetDeg, ...
  latLonSurface.centerDeltaObj, latLonSurface.bestLonLatSpan, latLonSurface.bestLatLonSpan, ...
  latLonSurface.centerLonLatSpan, latLonSurface.centerLatLonSpan, ...
  'VariableNames', {'modeTag', 'centerName', 'resolvedCenterName', 'fixedFdTag', 'fixedFdOffsetHz', ...
  'fixedFdAliasIndex', 'minLatOffsetDeg', 'minLonOffsetDeg', 'minDoaRadialOffsetDeg', ...
  'centerDeltaObj', 'bestLonLatSpan', 'bestLatLonSpan', 'centerLonLatSpan', 'centerLatLonSpan'});
end

function shouldRun = localShouldRunLatLonSurface(modeTag, requestedCenterName, resolvedCenterName, ...
  latLonSurfaceModeList, latLonSurfaceCenterNameList)
% Decide whether the optional lat/lon surface applies to this mode and center.
shouldRun = localStringMatchesAny(modeTag, latLonSurfaceModeList) ...
  && (localStringMatchesAny(requestedCenterName, latLonSurfaceCenterNameList) ...
  || localStringMatchesAny(resolvedCenterName, latLonSurfaceCenterNameList));
end

function fdOffsetTable = localResolveLatLonSurfaceFdOffsetTable(fdOffsetTagList, surfaceScan)
% Resolve fixed-fdRef tags for local latitude/longitude surface scans.
fdOffsetTagList = reshape(string(fdOffsetTagList), [], 1);
fdOffsetTagList = fdOffsetTagList(strlength(fdOffsetTagList) > 0);
if isempty(fdOffsetTagList)
  fdOffsetTagList = "surfaceMinFd";
end
numTag = numel(fdOffsetTagList);
fixedFdTag = strings(numTag, 1);
fixedFdOffsetHz = nan(numTag, 1);
for iTag = 1:numTag
  tagKey = lower(char(fdOffsetTagList(iTag)));
  switch tagKey
    case {'centerfd', 'center', 'fd0'}
      fixedFdTag(iTag) = "centerFd";
      fixedFdOffsetHz(iTag) = 0;
    case {'surfaceminfd', 'bestfd', 'besttooth', 'minfd'}
      fixedFdTag(iTag) = "surfaceMinFd";
      fixedFdOffsetHz(iTag) = surfaceScan.minFdOffsetHz;
    otherwise
      error('replayMfCombToothSurface:InvalidLatLonFdTag', ...
        'Unsupported latLonSurfaceFdOffsetTagList entry: %s', char(fdOffsetTagList(iTag)));
  end
end
[~, uniqueIdx] = unique(fixedFdTag + "|" + string(fixedFdOffsetHz), 'stable');
fdOffsetTable = table(fixedFdTag(uniqueIdx), fixedFdOffsetHz(uniqueIdx), ...
  'VariableNames', {'fixedFdTag', 'fixedFdOffsetHz'});
end

function matches = localStringMatchesAny(value, listValue)
% Case-insensitive string-list match; an empty list disables the match.
listValue = reshape(string(listValue), [], 1);
listValue = listValue(strlength(listValue) > 0);
if isempty(listValue)
  matches = false;
  return;
end
matches = any(strcmpi(string(value), listValue));
end

function spanValue = localFiniteSpan(valueVec)
% Return the finite max-min span of a vector.
valueVec = valueVec(:);
valueVec = valueVec(isfinite(valueVec));
if isempty(valueVec)
  spanValue = NaN;
else
  spanValue = max(valueVec) - min(valueVec);
end
end

function ratioValue = localSafeRatio(numValue, denValue)
% Return a protected ratio used only for scale diagnostics.
if ~(isscalar(numValue) && isscalar(denValue) && isfinite(numValue) && isfinite(denValue)) || abs(denValue) < eps
  ratioValue = NaN;
else
  ratioValue = numValue / denValue;
end
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

function localPlotReplay(lineScanCell, surfaceScanCell, couplingScanCell, aliasStepHz, surfaceClipDeltaObj, latLonSurfaceCell)
% Plot comb, surface, clipped heatmap, slice, and optional lat/lon diagnostics.
for iCase = 1:numel(lineScanCell)
  lineScan = lineScanCell{iCase};
  surfaceScan = surfaceScanCell{iCase};
  centerName = char(localGetFieldOrDefault(lineScan, 'resolvedCenterName', string(lineScan.centerName)));
  modeTag = char(localGetFieldOrDefault(lineScan, 'modeTag', "unknown"));
  titleTag = sprintf('%s | %s', modeTag, centerName);

  figure('Name', sprintf('comb and DoA-fdRef surface: %s', titleTag));

  subplot(2, 2, 1);
  plot(lineScan.deltaFdRef, lineScan.deltaObjVec, 'LineWidth', 1.2);
  hold on;
  localOverlayAliasLines(aliasStepHz, lineScan.deltaFdRef);
  xline(0, ':');
  grid on;
  xlabel('\Delta fdRef (Hz)');
  ylabel('\Delta objective');
  title(sprintf('fdRef line | %s', titleTag), 'Interpreter', 'none');

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

  localPlotSurfaceSliceDiagnostic(surfaceScan, aliasStepHz, surfaceClipDeltaObj, titleTag);
end
if nargin >= 6 && ~isempty(latLonSurfaceCell)
  if isstruct(latLonSurfaceCell)
    latLonSurfaceCell = {latLonSurfaceCell};
  end
  for iLatLon = 1:numel(latLonSurfaceCell)
    latLonSurface = latLonSurfaceCell{iLatLon};
    if isstruct(latLonSurface) && ~isempty(latLonSurface)
      localPlotLatLonSurface(latLonSurface, surfaceClipDeltaObj);
    end
  end
end
if nargin >= 3 && ~isempty(couplingScanCell)
  if isstruct(couplingScanCell)
    couplingScanCell = {couplingScanCell};
  end
  for iCoupling = 1:numel(couplingScanCell)
    couplingScan = couplingScanCell{iCoupling};
    if isstruct(couplingScan) && ~isempty(couplingScan)
      localPlotCouplingDiagnostic(couplingScan, aliasStepHz);
    end
  end
end
end

function localPlotSurfaceSliceDiagnostic(surfaceScan, aliasStepHz, surfaceClipDeltaObj, titleTag)
% Plot clipped heatmap and 1-D slices so weak DoA curvature remains visible.
surfaceScan = localAttachSurfaceSliceMetrics(surfaceScan);
clipMax = localResolveClipMax(surfaceScan.deltaObjMat, surfaceClipDeltaObj);
figure('Name', sprintf('DoA-fdRef slice diagnostics: %s', titleTag));

subplot(1, 3, 1);
imagesc(surfaceScan.fdOffsetGridHz, surfaceScan.doaOffsetGridDeg, surfaceScan.deltaObjMat);
set(gca, 'YDir', 'normal');
colorbar;
if isfinite(clipMax)
  caxis([0, clipMax]);
end
hold on;
xline(0, ':');
yline(0, ':');
localOverlayAliasLines(aliasStepHz, surfaceScan.fdOffsetGridHz);
grid on;
xlabel('\Delta fdRef (Hz)');
ylabel('DoA-line offset (deg)');
title(sprintf('clipped heatmap | <= %.3g', clipMax), 'Interpreter', 'none');

subplot(1, 3, 2);
plot(surfaceScan.doaOffsetGridDeg, surfaceScan.bestFdDoaSlice, 'LineWidth', 1.2);
hold on;
plot(surfaceScan.doaOffsetGridDeg, surfaceScan.centerFdDoaSlice, '--', 'LineWidth', 1.2);
xline(surfaceScan.minDoaOffsetDeg, ':');
yline(0, ':');
grid on;
xlabel('DoA-line offset (deg)');
ylabel('\Delta objective');
legend({'at best fd tooth', 'at center fd'}, 'Location', 'best');
title(sprintf('DoA slice span %.3g / %.3g', surfaceScan.bestFdDoaSpan, surfaceScan.centerFdDoaSpan), ...
  'Interpreter', 'none');

subplot(1, 3, 3);
plot(surfaceScan.fdOffsetGridHz, surfaceScan.bestDoaFdSlice, 'LineWidth', 1.2);
hold on;
plot(surfaceScan.fdOffsetGridHz, surfaceScan.centerDoaFdSlice, '--', 'LineWidth', 1.2);
xline(surfaceScan.minFdOffsetHz, ':');
yline(0, ':');
localOverlayAliasLines(aliasStepHz, surfaceScan.fdOffsetGridHz);
grid on;
xlabel('\Delta fdRef (Hz)');
ylabel('\Delta objective');
legend({'at best DoA', 'at center DoA'}, 'Location', 'best');
title(sprintf('fd slice span %.3g / %.3g', surfaceScan.bestDoaFdSpan, surfaceScan.centerDoaFdSpan), ...
  'Interpreter', 'none');
end

function localPlotLatLonSurface(latLonSurface, surfaceClipDeltaObj)
% Plot the optional fixed-fdRef latitude/longitude local surface.
modeTag = char(localGetFieldOrDefault(latLonSurface, 'modeTag', "unknown"));
centerName = char(localGetFieldOrDefault(latLonSurface, 'resolvedCenterName', "unknown"));
fixedFdTag = char(localGetFieldOrDefault(latLonSurface, 'fixedFdTag', "unknown"));
clipMax = localResolveClipMax(latLonSurface.deltaObjMat, surfaceClipDeltaObj);
figure('Name', sprintf('local lat/lon surface: %s | %s | %s', modeTag, centerName, fixedFdTag));

subplot(1, 3, 1);
imagesc(latLonSurface.lonOffsetGridDeg, latLonSurface.latOffsetGridDeg, latLonSurface.deltaObjMat);
set(gca, 'YDir', 'normal');
colorbar;
if isfinite(clipMax)
  caxis([0, clipMax]);
end
hold on;
xline(0, ':');
yline(0, ':');
plot(latLonSurface.minLonOffsetDeg, latLonSurface.minLatOffsetDeg, 'kx', 'LineWidth', 1.2);
grid on;
xlabel('longitude offset (deg)');
ylabel('latitude offset (deg)');
title(sprintf('%s | %s | %s | fdOffset %.3f Hz', modeTag, centerName, fixedFdTag, latLonSurface.fixedFdOffsetHz), ...
  'Interpreter', 'none');

subplot(1, 3, 2);
plot(latLonSurface.latOffsetGridDeg, latLonSurface.bestLonLatSlice, 'LineWidth', 1.2);
hold on;
plot(latLonSurface.latOffsetGridDeg, latLonSurface.centerLonLatSlice, '--', 'LineWidth', 1.2);
xline(latLonSurface.minLatOffsetDeg, ':');
yline(0, ':');
grid on;
xlabel('latitude offset (deg)');
ylabel('\Delta objective');
legend({'at best longitude', 'at center longitude'}, 'Location', 'best');
title(sprintf('lat slice span %.3g / %.3g', latLonSurface.bestLonLatSpan, latLonSurface.centerLonLatSpan), ...
  'Interpreter', 'none');

subplot(1, 3, 3);
plot(latLonSurface.lonOffsetGridDeg, latLonSurface.bestLatLonSlice, 'LineWidth', 1.2);
hold on;
plot(latLonSurface.lonOffsetGridDeg, latLonSurface.centerLatLonSlice, '--', 'LineWidth', 1.2);
xline(latLonSurface.minLonOffsetDeg, ':');
yline(0, ':');
grid on;
xlabel('longitude offset (deg)');
ylabel('\Delta objective');
legend({'at best latitude', 'at center latitude'}, 'Location', 'best');
title(sprintf('lon slice span %.3g / %.3g', latLonSurface.bestLatLonSpan, latLonSurface.centerLatLonSpan), ...
  'Interpreter', 'none');
end

function clipMax = localResolveClipMax(deltaObjMat, surfaceClipDeltaObj)
% Resolve a positive color-axis upper bound for clipped surface plots.
clipMax = NaN;
if nargin >= 2 && isscalar(surfaceClipDeltaObj) && isnumeric(surfaceClipDeltaObj) ...
    && isfinite(surfaceClipDeltaObj) && surfaceClipDeltaObj > 0
  clipMax = surfaceClipDeltaObj;
  return;
end
finiteObj = deltaObjMat(isfinite(deltaObjMat));
if ~isempty(finiteObj)
  clipMax = max(finiteObj(:));
end
end

function localPlotCouplingDiagnostic(couplingScan, aliasStepHz)
% Plot the truth-centered local surface and the best-fd ridge used for coupling.
modeTag = char(localGetFieldOrDefault(couplingScan, 'modeTag', "unknown"));
figure('Name', sprintf('truth-centered DoA-fdRef coupling: %s', modeTag));
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
title(sprintf('truth local coupling | %s: rho=%.3g, ridge=%.3g Hz/deg', ...
  modeTag, couplingScan.rhoDoaFd, couplingScan.ridgeSlopeHzPerDeg), 'Interpreter', 'none');
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
