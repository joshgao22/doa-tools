% scanMfGeoSatGeometryCrbAudit
% Screen candidate user locations and nested cooperative Starlink satellite sets
% using multi-frame geometry, MF-CP CRB health, EFIM conditioning, and
% DoA-fdRef coupling diagnostics.
%
% This scan is the location-level coarse wrapper around the fixed-location
% satellite set audit. It does not run MLE, does not enumerate all satellite
% combinations, and by default only screens nested L = 1,2,4 sets for the top
% few user locations. Use scanMfPairGeometryCrbAudit for the full
% L = 1,2,4,6,8,10 deep audit after a location is selected. Lightweight
% task checkpoint files are written under tmp/ when enabled.

clear; close all; clc;

%% Scan configuration

scanName = "scanMfGeoSatGeometryCrbAudit";
contextSeed = 253;
snrDb = 0;
saveSnapshot = true;
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
notifyTelegramEnable = false;
useParfor = true;
minTaskForParfor = 4;

% First-pass Starlink-inspired latitude sweep. Keep the current location as an anchor.
usrLlaList = [ ...
  37.78, 36.59, 0; ...
  25.00, 36.59, 0; ...
  35.00, 36.59, 0; ...
  40.00, 36.59, 0; ...
  45.00, 36.59, 0; ...
  50.00, 36.59, 0; ...
  53.00, 36.59, 0; ...
  55.00, 36.59, 0; ...
  60.00, 36.59, 0];
geoSelectionMode = "cooperativeBestScore";
geoScanMode = "coarseL4";
numGeoGreedyKeep = 3;

tleFileName = "statlink_20260318.tle";
utc0 = datetime([2026, 03, 18, 17, 08, 00], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');

numFrame = 10;
frameIntvlSec = 1 / 750;
signalModelTag = "constant-doa-affine-doppler-cp";
refFrameIdx = round((numFrame + 1) / 2);
satCountList = [1 2 4];

numSym = 32;
sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
carrierFreq = 18e9;
lightSpeed = 299792458;
minUsrElevationDeg = 15;
maxSatOffAxisDeg = 55;
accessBatchSize = 500;

numElem = [4 4];
maxReferenceCandidateEval = 16;
refSelectionMode = "cooperativeBestScore";
manualRefSatIdxGlobal = 4154;
maxCooperativeReferenceEval = 12;
maxCandidateSecondSat = 8;
maxGreedyCandidatePerStep = 12;
forcedCandidateSatIdxGlobal = 1165;
numTopReport = 12;

% Score weights are diagnostic only. Tables keep the raw metrics so locations
% and satellite sets can be re-ranked offline without changing CRB results.
scoreOpt = struct();
scoreOpt.angleGainWeight = 1.0;
scoreOpt.fdGainWeight = 0.5;
scoreOpt.condPenaltyWeight = 0.25;
scoreOpt.couplingPenaltyWeight = 1.0;
scoreOpt.couplingFractionPenaltyWeight = 0.5;
scoreOpt.phaseResidualPenaltyWeight = 0.5;
scoreOpt.lowGainTolerance = 1.02;
scoreOpt.highCondThreshold = 1e8;
scoreOpt.highCouplingThreshold = 0.65;
scoreOpt.phaseResidualStressRad = 0.10;

checkpointEnable = logical(checkpointEnable);
checkpointResume = logical(checkpointResume);
checkpointCleanupOnSuccess = logical(checkpointCleanupOnSuccess);
saveSnapshot = logical(saveSnapshot);
notifyTelegramEnable = logical(notifyTelegramEnable);
useParfor = logical(useParfor);

runTic = tic;
runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
repoRoot = localFindRepoRoot();
checkpointDir = "";
checkpointRunStateList = struct([]);
checkpointCleanupReport = table();
scanData = struct();

try
  %% Build context and scan tasks

  scanConfig = struct();
  scanConfig.scanName = string(scanName);
  scanConfig.contextSeed = double(contextSeed);
  scanConfig.snrDb = double(snrDb);
  scanConfig.usrLlaList = double(usrLlaList);
  scanConfig.usrLla = reshape(double(usrLlaList(1, :)), 3, 1);
  scanConfig.geoSelectionMode = string(geoSelectionMode);
  scanConfig.geoScanMode = string(geoScanMode);
  scanConfig.numGeoGreedyKeep = double(numGeoGreedyKeep);
  scanConfig.tleFileName = string(tleFileName);
  scanConfig.utc0 = utc0;
  scanConfig.numFrame = double(numFrame);
  scanConfig.frameIntvlSec = double(frameIntvlSec);
  scanConfig.signalModelTag = string(signalModelTag);
  scanConfig.refFrameIdx = double(refFrameIdx);
  scanConfig.satCountList = reshape(double(satCountList), 1, []);
  scanConfig.numSym = double(numSym);
  scanConfig.sampleRate = double(sampleRate);
  scanConfig.symbolRate = double(symbolRate);
  scanConfig.osf = double(osf);
  scanConfig.carrierFreq = double(carrierFreq);
  scanConfig.lightSpeed = double(lightSpeed);
  scanConfig.minUsrElevationDeg = double(minUsrElevationDeg);
  scanConfig.maxSatOffAxisDeg = double(maxSatOffAxisDeg);
  scanConfig.accessBatchSize = double(accessBatchSize);
  scanConfig.numElem = reshape(double(numElem), 1, []);
  scanConfig.maxReferenceCandidateEval = double(maxReferenceCandidateEval);
  scanConfig.refSelectionMode = string(refSelectionMode);
  scanConfig.manualRefSatIdxGlobal = double(manualRefSatIdxGlobal);
  scanConfig.maxCooperativeReferenceEval = double(maxCooperativeReferenceEval);
  scanConfig.maxCandidateSecondSat = double(maxCandidateSecondSat);
  scanConfig.maxGreedyCandidatePerStep = double(maxGreedyCandidatePerStep);
  scanConfig.forcedCandidateSatIdxGlobal = reshape(double(forcedCandidateSatIdxGlobal), 1, []);
  scanConfig.numTopReport = double(numTopReport);
  scanConfig.scoreOpt = scoreOpt;
  scanConfig.saveSnapshot = saveSnapshot;
  scanConfig.checkpointEnable = checkpointEnable;
  scanConfig.checkpointResume = checkpointResume;
  scanConfig.checkpointCleanupOnSuccess = checkpointCleanupOnSuccess;
  scanConfig.notifyTelegramEnable = notifyTelegramEnable;
  scanConfig.useParfor = useParfor;
  scanConfig.minTaskForParfor = double(minTaskForParfor);
  scanConfig.runKey = char(runKey);
  scanConfig.geoCheckpointToken = "";
  scanConfig.geoIdx = NaN;
  scanConfig.numTask = size(scanConfig.usrLlaList, 1);

  rng(scanConfig.contextSeed, 'twister');
  localPrintRuntimeLog("context", "Loading TLE and building pilot/context objects.");
  tle = tleread(char(scanConfig.tleFileName));
  wavelen = scanConfig.lightSpeed / scanConfig.carrierFreq;
  elemSpace = wavelen / 2;
  arrUpa = createUpa(scanConfig.numElem, elemSpace);
  [pilotSym, ~] = genPilotSymbol(1, scanConfig.numSym, 'pn', 1);
  pulseOpt = struct('rolloff', 0.25, 'span', 8);
  [pilotWave, waveInfo] = genPilotWaveform(pilotSym, scanConfig.symbolRate, scanConfig.osf, 'rrc', pulseOpt);
  pwrNoise = 1 / (10 ^ (scanConfig.snrDb / 10));
  utcVec = scanConfig.utc0 + seconds(((1:scanConfig.numFrame) - scanConfig.refFrameIdx) * scanConfig.frameIntvlSec);

  baseSharedData = struct();
  baseSharedData.scanName = scanName;
  baseSharedData.config = scanConfig;
  baseSharedData.tle = tle;
  baseSharedData.utcVec = utcVec;
  baseSharedData.wavelen = wavelen;
  baseSharedData.arrUpa = arrUpa;
  baseSharedData.pilotWave = pilotWave;
  baseSharedData.waveInfo = waveInfo;
  baseSharedData.pwrNoise = pwrNoise;

  checkpointDir = localDefaultCheckpointDir(repoRoot, scanName, scanConfig);
  printMfScanHeader(char(scanName), scanConfig, checkpointDir);
  geoCandidateTable = localBuildGeoCandidateTable(scanConfig.usrLlaList);

  geoResultCell = cell(height(geoCandidateTable), 1);
  geoSummaryRows = repmat(localEmptyGeoReferenceSummaryRow(), height(geoCandidateTable), 1);
  referencePreselectCell = cell(height(geoCandidateTable), 1);
  referenceCandidateCell = cell(height(geoCandidateTable), 1);
  refCoopScoreCell = cell(height(geoCandidateTable), 1);
  allPairRankCell = cell(height(geoCandidateTable), 1);

  localPrintRuntimeLog("geo", sprintf('Evaluating %d candidate user location(s).', height(geoCandidateTable)));
  for iGeo = 1:height(geoCandidateTable)
    geoConfig = localConfigForGeo(scanConfig, iGeo, scanConfig.usrLlaList(iGeo, :).');
    [geoResult, checkpointRunStateList] = localEvaluateGeoCandidate( ...
      iGeo, geoConfig, baseSharedData, repoRoot, checkpointRunStateList);
    geoResultCell{iGeo} = geoResult;
    geoSummaryRows(iGeo) = geoResult.summaryRow;
    referencePreselectCell{iGeo} = geoResult.referencePreselectTable;
    referenceCandidateCell{iGeo} = geoResult.referenceCandidateTable;
    refCoopScoreCell{iGeo} = geoResult.refCooperativeScoreTable;
    allPairRankCell{iGeo} = geoResult.allReferencePairRankTable;
  end

  [geoReferenceSummaryTable, selectedGeoIdxList] = localRankGeoReferenceSummary( ...
    geoSummaryRows, scanConfig.numGeoGreedyKeep);
  referencePreselectTable = localVertcatNonEmptyTables(referencePreselectCell);
  referenceCandidateTable = localVertcatNonEmptyTables(referenceCandidateCell);
  refCooperativeScoreTable = localVertcatNonEmptyTables(refCoopScoreCell);
  allReferencePairRankTable = localVertcatNonEmptyTables(allPairRankCell);

  %% Run scan batch

  localPrintRuntimeLog("geoGreedy", sprintf('Running greedy satellite-set scan for top %d location(s).', numel(selectedGeoIdxList)));
  geoGreedyStepCell = cell(numel(selectedGeoIdxList), 1);
  geoGreedyCandidateCell = cell(numel(selectedGeoIdxList), 1);
  geoSatSetSummaryCell = cell(numel(selectedGeoIdxList), 1);
  for iRank = 1:numel(selectedGeoIdxList)
    geoIdx = selectedGeoIdxList(iRank);
    geoResult = geoResultCell{geoIdx};
    [greedyStepTableTmp, greedyCandidateTableTmp, satSetSummaryTableTmp, checkpointRunStateList] = ...
      localRunGeoGreedySelection(geoResult, scanConfig, repoRoot, checkpointRunStateList);
    if ~isempty(greedyStepTableTmp)
      greedyStepTableTmp.rankByGeoScore = repmat(iRank, height(greedyStepTableTmp), 1);
      greedyStepTableTmp = movevars(greedyStepTableTmp, 'rankByGeoScore', 'Before', 1);
    end
    if ~isempty(greedyCandidateTableTmp)
      greedyCandidateTableTmp.rankByGeoScore = repmat(iRank, height(greedyCandidateTableTmp), 1);
      greedyCandidateTableTmp = movevars(greedyCandidateTableTmp, 'rankByGeoScore', 'Before', 1);
    end
    if ~isempty(satSetSummaryTableTmp)
      satSetSummaryTableTmp.rankByGeoScore = repmat(iRank, height(satSetSummaryTableTmp), 1);
      satSetSummaryTableTmp = movevars(satSetSummaryTableTmp, 'rankByGeoScore', 'Before', 1);
    end
    geoGreedyStepCell{iRank} = greedyStepTableTmp;
    geoGreedyCandidateCell{iRank} = greedyCandidateTableTmp;
    geoSatSetSummaryCell{iRank} = satSetSummaryTableTmp;
  end
  geoGreedyStepTable = localVertcatNonEmptyTables(geoGreedyStepCell);
  geoGreedyCandidateTable = localVertcatNonEmptyTables(geoGreedyCandidateCell);
  geoSatSetSummaryTable = localVertcatNonEmptyTables(geoSatSetSummaryCell);
  paperCandidateSetTable = localBuildPaperCandidateSetTable(geoSatSetSummaryTable, geoReferenceSummaryTable);

  %% Data storage

  checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunStateList, checkpointDir);
  if scanConfig.checkpointEnable && scanConfig.checkpointCleanupOnSuccess
    localPrintRuntimeLog("checkpoint", "Cleaning completed checkpoint runs.");
    checkpointCleanupReport = localCleanupCompletedCheckpointRuns(checkpointRunStateList);
    checkpointSummaryTable = localBuildCheckpointSummaryTable(checkpointRunStateList, checkpointDir);
  end

  localPrintRuntimeLog("summary", "Building output tables and lightweight plot data.");
  contextSummaryTable = localBuildGeoContextSummaryTable(scanConfig, geoCandidateTable, geoReferenceSummaryTable);
  bestGeoSummaryTable = localPreviewRows(geoReferenceSummaryTable, min(scanConfig.numTopReport, height(geoReferenceSummaryTable)));

  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.contextSummaryTable = contextSummaryTable;
  scanData.geoCandidateTable = geoCandidateTable;
  scanData.geoReferenceSummaryTable = geoReferenceSummaryTable;
  scanData.referencePreselectTable = referencePreselectTable;
  scanData.referenceCandidateTable = referenceCandidateTable;
  scanData.refCooperativeScoreTable = refCooperativeScoreTable;
  scanData.allReferencePairRankTable = allReferencePairRankTable;
  scanData.geoGreedyStepTable = geoGreedyStepTable;
  scanData.geoGreedyCandidateTable = geoGreedyCandidateTable;
  scanData.geoSatSetSummaryTable = geoSatSetSummaryTable;
  scanData.paperCandidateSetTable = paperCandidateSetTable;
  scanData.bestGeoSummaryTable = bestGeoSummaryTable;
  scanData.checkpointSummaryTable = checkpointSummaryTable;
  scanData.checkpointCleanupReport = checkpointCleanupReport;
  scanData.checkpointDir = checkpointDir;
  scanData.plotData = localBuildGeoPlotData(geoReferenceSummaryTable, geoSatSetSummaryTable, paperCandidateSetTable);
  scanData.elapsedSec = toc(runTic);
  scanData = finalizeMfScanResult(scanData, "");

  if scanConfig.saveSnapshot
    saveOpt = struct();
    saveOpt.includeVars = {'scanData'};
    saveOpt.extraMeta = struct('scanName', char(scanName), ...
      'tleFileName', char(scanConfig.tleFileName), ...
      'signalModelTag', char(scanConfig.signalModelTag), ...
      'numGeoCandidate', height(geoCandidateTable));
    saveOpt.verbose = true;
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end

  %% Summary output and plotting

  printMfScanSection('Geo / satellite CRB context summary', scanData.contextSummaryTable);
  printMfScanSection('Geo candidates', scanData.geoCandidateTable);
  printMfScanSection('Geo cooperative reference scores', localPreviewRows(scanData.geoReferenceSummaryTable, scanConfig.numTopReport));
  printMfScanSection('Selected-location satellite-set summary', scanData.geoSatSetSummaryTable);
  printMfScanSection('Paper candidate location / satellite sets', localPreviewRows(scanData.paperCandidateSetTable, scanConfig.numTopReport));
  printMfScanSection('Checkpoint summary', scanData.checkpointSummaryTable);
  scanData.plotData = localPlotGeoScan(scanData);

  metricLineList = [ ...
    "• geo candidates: <code>" + string(height(scanData.geoCandidateTable)) + "</code>"; ...
    "• top location: <code>" + string(scanData.contextSummaryTable.topUsrLla(1)) + "</code>"; ...
    "• top reference: <code>" + string(scanData.contextSummaryTable.topRefSatIdxGlobal(1)) + "</code>"; ...
    "• paper candidates: <code>" + string(height(scanData.paperCandidateSetTable)) + "</code>" ...
    ];
  notifyMfScanStatus(struct( ...
    'scanName', string(scanName), ...
    'statusText', "DONE", ...
    'config', scanConfig, ...
    'snapshotFile', localGetFieldOrDefault(scanData, 'snapshotFile', ""), ...
    'checkpointDir', scanData.checkpointDir, ...
    'elapsedSec', scanData.elapsedSec, ...
    'metricLineList', metricLineList, ...
    'commentLineList', [ ...
      "Geo + satellite CRB-health scan completed."; ...
      "Default geo scan is coarse L=1/2/4; run scanMfPairGeometryCrbAudit for full L=1/2/4/6/8/10 deep audit."; ...
      "MLE is not run in this scan; validate chosen sets with paper-facing MLE/CRB scans."]));

catch ME
  if ~exist('scanConfig', 'var')
    scanConfig = struct();
  end
  localPrintCheckpointFailureHint(checkpointDir);
  notifyMfScanStatus(struct( ...
    'scanName', string(scanName), ...
    'statusText', "FAILED", ...
    'config', scanConfig, ...
    'checkpointDir', checkpointDir, ...
    'elapsedSec', toc(runTic), ...
    'errorObj', ME));
  rethrow(ME);
end

%% Local helpers


function tableOut = localBuildGeoCandidateTable(usrLlaList)
%LOCALBUILDGEOCANDIDATETABLE Build configured user-location candidates.
numGeo = size(usrLlaList, 1);
geoIdx = (1:numGeo).';
usrLatDeg = usrLlaList(:, 1);
usrLonDeg = usrLlaList(:, 2);
usrAltM = usrLlaList(:, 3);
usrLlaStr = strings(numGeo, 1);
for iGeo = 1:numGeo
  usrLlaStr(iGeo) = string(mat2str(usrLlaList(iGeo, :)));
end
tableOut = table(geoIdx, usrLatDeg, usrLonDeg, usrAltM, usrLlaStr);
end

function [geoResult, checkpointRunStateList] = localEvaluateGeoCandidate( ...
  geoIdx, geoConfig, baseSharedData, repoRoot, checkpointRunStateList)
%LOCALEVALUATEGEOCANDIDATE Run reference and cooperative-pair audit for one location.

geoResult = struct();
geoResult.geoIdx = double(geoIdx);
geoResult.usrLla = geoConfig.usrLla;
geoResult.usrLlaStr = string(mat2str(geoConfig.usrLla(:).'));
geoResult.availableSatIdxGlobal = [];
geoResult.referencePreselectTable = table();
geoResult.referenceCandidateTable = table();
geoResult.refCooperativeScoreTable = localEmptyRefCooperativeScoreTable();
geoResult.allReferencePairRankTable = table();
geoResult.pairRankTable = table();
geoResult.selectedRefRow = table();
geoResult.baselineMetrics = struct();
geoResult.sharedData = struct();
geoResult.summaryRow = localEmptyGeoReferenceSummaryRow();

stagePrefix = sprintf('geo%02d', geoIdx);
localPrintRuntimeLog(stagePrefix, sprintf('Evaluating user location %s.', char(geoResult.usrLlaStr)));

[satIdxRef, satAccessRef, satStateRef] = findVisibleSatFromTle( ...
  geoConfig.utc0, baseSharedData.tle, geoConfig.usrLla, ...
  geoConfig.minUsrElevationDeg, geoConfig.maxSatOffAxisDeg, geoConfig.accessBatchSize);
availableSatIdxGlobal = reshape(double(satIdxRef.available), 1, []);
geoResult.availableSatIdxGlobal = availableSatIdxGlobal;
localPrintRuntimeLog(stagePrefix, sprintf('Reference-epoch visible candidates: %d.', numel(availableSatIdxGlobal)));

sharedData = baseSharedData;
sharedData.config = geoConfig;
sharedData.satAccessRef = satAccessRef;
sharedData.satStateRef = satStateRef;
geoResult.sharedData = sharedData;

if isempty(availableSatIdxGlobal)
  geoResult.summaryRow = localBuildGeoReferenceSummaryRow(geoIdx, geoConfig.usrLla, ...
    0, 0, table(), localEmptyRefCooperativeScoreRow());
  geoResult.summaryRow.status = "no-visible-sat";
  geoResult.summaryRow.geoClass = "weak-coverage";
  return;
end

referencePreselectTable = localBuildReferencePreselectTable(availableSatIdxGlobal, satAccessRef, geoConfig);
referenceEvalTaskList = localBuildReferenceEvalTaskList(referencePreselectTable, geoConfig);
useParforResolved = localCanUseParfor(geoConfig.useParfor, geoConfig.minTaskForParfor, numel(referenceEvalTaskList));
[referenceOutCell, refRunState] = localRunTaskBatchWithState(referenceEvalTaskList, sharedData, ...
  string(stagePrefix) + "Ref", geoConfig, repoRoot, useParforResolved);
if ~isempty(refRunState)
  checkpointRunStateList = [checkpointRunStateList; refRunState(:)]; %#ok<AGROW>
end
referenceCandidateTable = localBuildSetResultTable(referenceOutCell);
referenceCandidateTable = localAttachReferenceScores(referenceCandidateTable, geoConfig);

if isempty(referenceCandidateTable)
  geoResult.summaryRow = localBuildGeoReferenceSummaryRow(geoIdx, geoConfig.usrLla, ...
    numel(availableSatIdxGlobal), 0, table(), localEmptyRefCooperativeScoreRow());
  geoResult.summaryRow.status = "no-reference-candidate";
  geoResult.summaryRow.geoClass = "bad-geometry";
  return;
end

[refCooperativeScoreTable, allReferencePairRankTable, checkpointRunStateList] = ...
  localBuildCooperativeReferenceAudit(referenceCandidateTable, availableSatIdxGlobal, ...
    sharedData, geoConfig, repoRoot, checkpointRunStateList);
selectedRefRow = localSelectReferenceRowFromCooperative(referenceCandidateTable, ...
  refCooperativeScoreTable, geoConfig);
refSatIdxGlobal = selectedRefRow.refSatIdxGlobal;
baselineMetrics = localFillSetResultDefaults(table2struct(selectedRefRow));
baselineMetrics.taskType = "baseline";

pairRankTable = allReferencePairRankTable(allReferencePairRankTable.refSatIdxGlobal == refSatIdxGlobal, :);
if ~isempty(pairRankTable)
  pairRankTable = localBuildPairRankTable(pairRankTable);
end

coopRow = localFindCoopScoreRow(refCooperativeScoreTable, refSatIdxGlobal);
summaryRow = localBuildGeoReferenceSummaryRow(geoIdx, geoConfig.usrLla, ...
  numel(availableSatIdxGlobal), height(referenceCandidateTable), selectedRefRow, coopRow);

geoResult.referencePreselectTable = localAddGeoColumns(referencePreselectTable, geoIdx, geoConfig.usrLla);
geoResult.referenceCandidateTable = localAddGeoColumns(referenceCandidateTable, geoIdx, geoConfig.usrLla);
geoResult.refCooperativeScoreTable = localAddGeoColumns(refCooperativeScoreTable, geoIdx, geoConfig.usrLla);
geoResult.allReferencePairRankTable = localAddGeoColumns(allReferencePairRankTable, geoIdx, geoConfig.usrLla);
geoResult.pairRankTable = localAddGeoColumns(pairRankTable, geoIdx, geoConfig.usrLla);
geoResult.selectedRefRow = selectedRefRow;
geoResult.baselineMetrics = baselineMetrics;
geoResult.summaryRow = summaryRow;
geoResult.sharedData.refSatIdxGlobal = refSatIdxGlobal;
geoResult.sharedData.baselineMetrics = baselineMetrics;
localPrintRuntimeLog(stagePrefix, sprintf('Selected reference %d; best pair %d; geo score %.3g.', ...
  double(refSatIdxGlobal), double(summaryRow.bestPairSatIdxGlobal), double(summaryRow.geoPreScore)));
end

function geoConfig = localConfigForGeo(scanConfig, geoIdx, usrLla)
%LOCALCONFIGFORGEO Return a per-location config with stable checkpoint token.
geoConfig = scanConfig;
geoConfig.geoIdx = double(geoIdx);
geoConfig.geoCheckpointToken = sprintf('geo%02d', geoIdx);
geoConfig.usrLla = reshape(double(usrLla), 3, 1);
end

function coopRow = localFindCoopScoreRow(refCooperativeScoreTable, refSatIdxGlobal)
%LOCALFINDCOOPSCOREROW Return one cooperative reference row or an empty row.
coopRow = localEmptyRefCooperativeScoreRow();
if isempty(refCooperativeScoreTable) || ~istable(refCooperativeScoreTable)
  return;
end
matchIdx = find(refCooperativeScoreTable.refSatIdxGlobal == refSatIdxGlobal, 1);
if isempty(matchIdx)
  return;
end
coopRow = table2struct(refCooperativeScoreTable(matchIdx, :));
end

function row = localEmptyGeoReferenceSummaryRow()
%LOCALEMPTYGEOREFERENCESUMMARYROW Return one empty geo-reference summary row.
row = struct();
row.geoIdx = NaN;
row.usrLlaStr = "";
row.numAvailableAtRef = NaN;
row.numReferenceCandidate = NaN;
row.selectedRefSatIdxGlobal = NaN;
row.selectedRefSatName = "";
row.selectedRefElevRefDeg = NaN;
row.selectedRefElevMinDeg = NaN;
row.selectedRefScore = NaN;
row.selectedRefCoopScore = NaN;
row.selectedRefCoopClass = "";
row.bestPairSatIdxGlobal = NaN;
row.bestPairScore = NaN;
row.bestPairClass = "";
row.bestPairAngleCrbDegCpK = NaN;
row.bestPairFdRefCrbHzCpU = NaN;
row.bestPairEfimCondCpU = NaN;
row.bestPairFdRefDoaCanonicalCorrCpU = NaN;
row.bestPairCouplingFractionCpU = NaN;
row.numPair = NaN;
row.numHealthyPair = NaN;
row.numHardCoupledPair = NaN;
row.healthyPairRate = NaN;
row.geoPreScore = NaN;
row.geoClass = "";
row.status = "";
end

function row = localBuildGeoReferenceSummaryRow(geoIdx, usrLla, numAvailable, numReferenceCandidate, selectedRefRow, coopRow)
%LOCALBUILDGEOREFERENCESUMMARYROW Build compact per-location summary.
row = localEmptyGeoReferenceSummaryRow();
row.geoIdx = double(geoIdx);
row.usrLlaStr = string(mat2str(usrLla(:).'));
row.numAvailableAtRef = double(numAvailable);
row.numReferenceCandidate = double(numReferenceCandidate);
if ~isempty(selectedRefRow) && istable(selectedRefRow)
  row.selectedRefSatIdxGlobal = double(selectedRefRow.refSatIdxGlobal(1));
  row.selectedRefSatName = string(selectedRefRow.selectedSatNameStr(1));
  row.selectedRefElevRefDeg = double(selectedRefRow.elevRefDeg(1));
  row.selectedRefElevMinDeg = double(selectedRefRow.elevMinDeg(1));
  row.selectedRefScore = double(selectedRefRow.refScore(1));
end
if isstruct(coopRow)
  row.selectedRefCoopScore = double(localGetFieldOrDefault(coopRow, 'refCoopScore', NaN));
  row.selectedRefCoopClass = string(localGetFieldOrDefault(coopRow, 'refCoopClass', ""));
  row.bestPairSatIdxGlobal = double(localGetFieldOrDefault(coopRow, 'bestPairSatIdxGlobal', NaN));
  row.bestPairScore = double(localGetFieldOrDefault(coopRow, 'bestPairScore', NaN));
  row.bestPairClass = string(localGetFieldOrDefault(coopRow, 'bestPairClass', ""));
  row.bestPairAngleCrbDegCpK = double(localGetFieldOrDefault(coopRow, 'bestPairAngleCrbDegCpK', NaN));
  row.bestPairFdRefCrbHzCpU = double(localGetFieldOrDefault(coopRow, 'bestPairFdRefCrbHzCpU', NaN));
  row.bestPairEfimCondCpU = double(localGetFieldOrDefault(coopRow, 'bestPairEfimCondCpU', NaN));
  row.bestPairFdRefDoaCanonicalCorrCpU = double(localGetFieldOrDefault(coopRow, 'bestPairFdRefDoaCanonicalCorrCpU', NaN));
  row.bestPairCouplingFractionCpU = double(localGetFieldOrDefault(coopRow, 'bestPairCouplingFractionCpU', NaN));
  row.numPair = double(localGetFieldOrDefault(coopRow, 'numPair', NaN));
  row.numHealthyPair = double(localGetFieldOrDefault(coopRow, 'numHealthyPair', NaN));
  row.numHardCoupledPair = double(localGetFieldOrDefault(coopRow, 'numHardCoupledPair', NaN));
  row.healthyPairRate = double(localGetFieldOrDefault(coopRow, 'healthyPairRate', NaN));
end
row.geoPreScore = localCalcGeoPreScore(row);
row.geoClass = localClassifyGeoCandidate(row);
row.status = "ok";
end

function score = localCalcGeoPreScore(row)
%LOCALCALCGEOPRESCORE Score a location by cooperative-reference health.
score = row.selectedRefCoopScore + 0.25 * log1p(max(row.numAvailableAtRef, 0)) ...
  - 0.5 * max(row.bestPairCouplingFractionCpU, 0);
if isfinite(row.bestPairEfimCondCpU)
  score = score - 0.10 * log10(max(row.bestPairEfimCondCpU, 1));
end
end

function geoClass = localClassifyGeoCandidate(row)
%LOCALCLASSIFYGEOCANDIDATE Classify location-level cooperative suitability.
if ~isfinite(row.numAvailableAtRef) || row.numAvailableAtRef < 10
  geoClass = "weak-coverage";
elseif row.numHealthyPair > 0
  geoClass = "good-geo-candidate";
elseif isfinite(row.bestPairCouplingFractionCpU) && row.bestPairCouplingFractionCpU <= 0.25
  geoClass = "multi-sat-promising";
elseif isfinite(row.bestPairScore)
  geoClass = "two-sat-hard-but-multisat-promising";
else
  geoClass = "bad-geometry";
end
end

function tableOut = localAddGeoColumns(tableIn, geoIdx, usrLla)
%LOCALADDGEOCOLUMNS Attach geo identity columns to an existing table.
tableOut = tableIn;
if isempty(tableOut) || ~istable(tableOut)
  return;
end
tableOut.geoIdx = repmat(double(geoIdx), height(tableOut), 1);
tableOut.usrLlaStr = repmat(string(mat2str(usrLla(:).')), height(tableOut), 1);
tableOut = movevars(tableOut, {'geoIdx', 'usrLlaStr'}, 'Before', 1);
end

function [geoReferenceSummaryTable, rankGeoIdx] = localRankGeoReferenceSummary(rowList, maxGeoGreedyKeep)
%LOCALRANKGEOREFERENCESUMMARY Build and rank per-location summary rows.
if isempty(rowList)
  geoReferenceSummaryTable = struct2table(localEmptyGeoReferenceSummaryRow());
  geoReferenceSummaryTable(1, :) = [];
  rankGeoIdx = [];
  return;
end
geoReferenceSummaryTable = struct2table(rowList);
geoReferenceSummaryTable = sortrows(geoReferenceSummaryTable, ...
  {'geoPreScore', 'selectedRefCoopScore', 'numAvailableAtRef'}, {'descend', 'descend', 'descend'});
geoReferenceSummaryTable.rankByGeoScore = (1:height(geoReferenceSummaryTable)).';
geoReferenceSummaryTable = movevars(geoReferenceSummaryTable, 'rankByGeoScore', 'Before', 'geoIdx');
validTable = geoReferenceSummaryTable(geoReferenceSummaryTable.status == "ok" & isfinite(geoReferenceSummaryTable.geoPreScore), :);
keepCount = min(height(validTable), max(1, double(maxGeoGreedyKeep)));
rankGeoIdx = double(validTable.geoIdx(1:keepCount));
end

function [greedyStepTable, greedyCandidateTable, satSetSummaryTable, checkpointRunStateList] = ...
  localRunGeoGreedySelection(geoResult, scanConfig, repoRoot, checkpointRunStateList)
%LOCALRUNGEOCOUNTSELECTION Run greedy satellite selection for one ranked geo candidate.
geoIdx = geoResult.geoIdx;
geoConfig = localConfigForGeo(scanConfig, geoIdx, geoResult.usrLla);
sharedData = geoResult.sharedData;
sharedData.config = geoConfig;
sharedData.refSatIdxGlobal = double(geoResult.selectedRefRow.refSatIdxGlobal(1));
sharedData.baselineMetrics = geoResult.baselineMetrics;
pairRankTable = geoResult.pairRankTable;
if isempty(pairRankTable)
  greedyStepTable = table();
  greedyCandidateTable = table();
  satSetSummaryTable = table();
  return;
end
localPrintRuntimeLog("geoGreedy", sprintf('Running greedy set selection for geo %d, ref %d.', ...
  geoIdx, double(sharedData.refSatIdxGlobal)));
[greedyStepTable, greedyCandidateTable, satSetSummaryTable, checkpointRunStateList] = ...
  localRunGreedySatelliteSelection(pairRankTable, sharedData.refSatIdxGlobal, sharedData, ...
    geoConfig, repoRoot, checkpointRunStateList);
greedyStepTable = localAddGeoColumns(greedyStepTable, geoIdx, geoResult.usrLla);
greedyCandidateTable = localAddGeoColumns(greedyCandidateTable, geoIdx, geoResult.usrLla);
satSetSummaryTable = localAddGeoColumns(satSetSummaryTable, geoIdx, geoResult.usrLla);
end

function paperTable = localBuildPaperCandidateSetTable(geoSatSetSummaryTable, geoReferenceSummaryTable)
%LOCALBUILDPAPERCANDIDATESETTABLE Build compact candidate location/set table.
if isempty(geoSatSetSummaryTable)
  paperTable = table();
  return;
end
paperTable = geoSatSetSummaryTable;
paperTable.paperSetScore = nan(height(paperTable), 1);
paperTable.paperCandidateClass = strings(height(paperTable), 1);
for iRow = 1:height(paperTable)
  geoIdx = paperTable.geoIdx(iRow);
  geoRows = paperTable(paperTable.geoIdx == geoIdx, :);
  baseRows = geoRows(geoRows.numSat == min(geoRows.numSat), :);
  if isempty(baseRows)
    continue;
  end
  baseAngle = baseRows.angleCrbDegCpU(1);
  baseFd = baseRows.fdRefCrbHzCpU(1);
  angleGain = localSafeRatio(baseAngle, paperTable.angleCrbDegCpU(iRow));
  fdGain = localSafeRatio(baseFd, paperTable.fdRefCrbHzCpU(iRow));
  paperTable.paperSetScore(iRow) = log(max(angleGain, eps)) + 0.7 * log(max(fdGain, eps)) ...
    - 0.25 * log10(max(paperTable.efimCondCpU(iRow), 1)) ...
    - paperTable.fdRefCouplingFractionCpU(iRow);
  paperTable.paperCandidateClass(iRow) = localClassifyPaperSetCandidate(paperTable(iRow, :), angleGain, fdGain);
end
paperTable = sortrows(paperTable, {'paperSetScore', 'numSat'}, {'descend', 'ascend'});
paperTable.rankByPaperSetScore = (1:height(paperTable)).';
paperTable = movevars(paperTable, 'rankByPaperSetScore', 'Before', 1);
if ~isempty(geoReferenceSummaryTable)
  [~, loc] = ismember(paperTable.geoIdx, geoReferenceSummaryTable.geoIdx);
  valid = loc > 0;
  selectedRefSatIdxGlobal = nan(height(paperTable), 1);
  selectedRefSatName = strings(height(paperTable), 1);
  geoClass = strings(height(paperTable), 1);
  selectedRefSatIdxGlobal(valid) = geoReferenceSummaryTable.selectedRefSatIdxGlobal(loc(valid));
  selectedRefSatName(valid) = geoReferenceSummaryTable.selectedRefSatName(loc(valid));
  geoClass(valid) = geoReferenceSummaryTable.geoClass(loc(valid));
  paperTable.selectedRefSatIdxGlobal = selectedRefSatIdxGlobal;
  paperTable.selectedRefSatName = selectedRefSatName;
  paperTable.geoClass = geoClass;
end
end

function className = localClassifyPaperSetCandidate(row, angleGain, fdGain)
%LOCALCLASSIFYPAPERSETCANDIDATE Classify final paper-candidate set.
if row.numSat <= 2
  className = "minimal-hard-coupled";
elseif row.fdRefCouplingFractionCpU <= 0.20 && fdGain >= 1.5 && angleGain >= 10
  className = "paper-main-candidate";
elseif row.fdRefCouplingFractionCpU <= 0.35 && fdGain >= 1.2
  className = "paper-promising-candidate";
else
  className = "stress-or-diagnostic-candidate";
end
end

function contextSummaryTable = localBuildGeoContextSummaryTable(config, geoCandidateTable, geoReferenceSummaryTable)
%LOCALBUILDGEOCONTEXTSUMMARYTABLE Build compact context table for geo scan.
if isempty(geoReferenceSummaryTable)
  selectedGeoIdx = NaN;
  selectedRef = NaN;
  selectedUsrLla = "";
else
  selectedGeoIdx = geoReferenceSummaryTable.geoIdx(1);
  selectedRef = geoReferenceSummaryTable.selectedRefSatIdxGlobal(1);
  selectedUsrLla = geoReferenceSummaryTable.usrLlaStr(1);
end
contextSummaryTable = table( ...
  string(config.tleFileName), ...
  string(config.signalModelTag), ...
  string(config.geoSelectionMode), ...
  string(localGetFieldOrDefault(config, 'geoScanMode', "")), ...
  string(config.refSelectionMode), ...
  string(config.utc0), ...
  config.numFrame, ...
  config.frameIntvlSec, ...
  height(geoCandidateTable), ...
  selectedGeoIdx, ...
  selectedUsrLla, ...
  selectedRef, ...
  'VariableNames', { ...
    'tleFileName', 'signalModelTag', 'geoSelectionMode', 'geoScanMode', 'refSelectionMode', 'utc0', ...
    'numFrame', 'frameIntvlSec', 'numGeoCandidate', 'topGeoIdx', 'topUsrLla', 'topRefSatIdxGlobal'});
end

function plotData = localBuildGeoPlotData(geoReferenceSummaryTable, geoSatSetSummaryTable, paperCandidateSetTable)
%LOCALBUILDGEOPLOTDATA Store lightweight plotting data for the geo scan.
plotData = struct();
plotData.geoReferenceSummaryTable = geoReferenceSummaryTable;
plotData.geoSatSetSummaryTable = geoSatSetSummaryTable;
plotData.paperCandidateSetTable = paperCandidateSetTable;
end

function scanData = localPlotGeoScan(scanData)
%LOCALPLOTGEOSCAN Draw quick-look geo and selected-set trends.
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  return;
end
geoTable = scanData.plotData.geoReferenceSummaryTable;
if istable(geoTable) && ~isempty(geoTable)
  topCount = min(12, height(geoTable));
  figure('Name', 'Geo cooperative score ranking');
  bar(geoTable.geoPreScore(1:topCount));
  grid on;
  xlabel('Geo CRB-health rank');
  ylabel('Geo pre-score');
  title('Candidate user locations by cooperative CRB health');
  set(gca, 'XTick', 1:topCount, 'XTickLabel', geoTable.usrLlaStr(1:topCount));
  xtickangle(30);
end
setTable = scanData.plotData.geoSatSetSummaryTable;
if istable(setTable) && ~isempty(setTable)
  figure('Name', 'Geo selected-set CRB trends');
  hold on;
  geoList = unique(setTable.geoIdx, 'stable');
  for iGeo = 1:numel(geoList)
    subTable = setTable(setTable.geoIdx == geoList(iGeo), :);
    plot(subTable.numSat, subTable.fdRefCrbHzCpU, '-o', 'DisplayName', "geo" + string(geoList(iGeo)));
  end
  hold off;
  grid on;
  xlabel('Number of satellites');
  ylabel('CP-U fdRef CRB std (Hz)');
  title('Selected-set fdRef CRB across top locations');
  legend('Location', 'best');
end
end

function rows = localPreviewByRows(tableIn, maxRows)
%LOCALPREVIEWBYROWS Preview helper with robust empty handling.
rows = localPreviewRows(tableIn, maxRows);
end

function tableOut = localVertcatNonEmptyTables(tableCell)
%LOCALVERTCATNONEMPTYTABLES Concatenate non-empty tables with safe empty fallback.
validMask = ~cellfun(@isempty, tableCell);
if ~any(validMask)
  tableOut = table();
  return;
end
tableOut = vertcat(tableCell{validMask});
end


function localPrintRuntimeLog(stageName, messageText)
%LOCALPRINTRUNTIMELOG Print a lightweight stage-level runtime message.
fprintf('[%s] %-12s %s\n', datestr(now, 'HH:MM:SS'), char(string(stageName) + ":"), char(string(messageText)));
end

function taskOut = localRunSetCheckpointTask(task, sharedData)
%LOCALRUNSETCHECKPOINTTASK Run one checkpointed CRB set task.
taskOut = localRunOneSetTask(task, sharedData);
end

function taskOut = localRunOneSetTask(task, sharedData)
%LOCALRUNONESETTASK Evaluate one satellite set CRB task.
taskOut = localEvaluateSatSetCrb(task.taskId, task.trialSatIdxGlobal, ...
  task.taskType, task.refSatIdxGlobal, sharedData);
taskOut.taskId = task.taskId;
taskOut.taskType = string(task.taskType);
taskOut.candSatIdxGlobal = localGetTaskCandidate(task);
taskOut.baseSelectedSatIdxGlobalStr = string(mat2str(task.baseSelectedSatIdxGlobal));
end

function taskOut = localEvaluateSatSetCrb(taskId, trialSatIdxGlobal, taskType, refSatIdxGlobal, sharedData)
%LOCALEVALUATESATSETCRB Evaluate geometry and MF CRB metrics for one satellite set.

config = sharedData.config;
trialSatIdxGlobal = reshape(double(trialSatIdxGlobal), 1, []);
numSat = numel(trialSatIdxGlobal);
row = localEmptySetResultRow();
row.taskId = double(taskId);
row.taskType = string(taskType);
row.refSatIdxGlobal = double(refSatIdxGlobal);
row.numSat = numSat;
row.selectedSatIdxGlobalStr = string(mat2str(trialSatIdxGlobal));
row.selectedSatNameStr = localJoinSatNames(sharedData.satAccessRef, trialSatIdxGlobal);
if numSat >= 2
  row.candSatIdxGlobal = trialSatIdxGlobal(end);
  row.candSatName = localExtractSatName(sharedData.satAccessRef, row.candSatIdxGlobal);
end

try
  sceneSeq = genMultiFrameScene(sharedData.utcVec, sharedData.tle, config.usrLla, ...
    trialSatIdxGlobal, [], sharedData.arrUpa, config.minUsrElevationDeg, ...
    config.maxSatOffAxisDeg, "satellite", refSatIdxGlobal, config.refFrameIdx);
  row.accessAllWindow = all(sceneSeq.access(:));
  row.refSatIdxLocal = localFindLocalSatIdx(sceneSeq.satIdx, refSatIdxGlobal);
  if ~row.accessAllWindow
    row.skipReason = "notAllAccessOverWindow";
    row.status = "skipped";
    row = localAttachSceneGeometry(row, sceneSeq, sharedData, trialSatIdxGlobal);
    taskOut = row;
    return;
  end

  row = localAttachSceneGeometry(row, sceneSeq, sharedData, trialSatIdxGlobal);
  linkParamCell = cell(1, sceneSeq.numFrame);
  for iFrame = 1:sceneSeq.numFrame
    linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, sharedData.wavelen);
  end
  truthDyn = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, ...
    sharedData.waveInfo.sampleRate, 1);
  row = localAttachDynamicTruth(row, truthDyn, sceneSeq);

  pathGain = ones(sceneSeq.numSat, sceneSeq.numFrame);
  noiseVar = sharedData.pwrNoise * ones(sceneSeq.numSat, sceneSeq.numFrame);
  crbOptBase = struct();
  crbOptBase.doaType = 'latlon';
  crbOptBase.phaseMode = 'continuous';
  crbOptBase.steeringMode = 'framewise';
  % The scan evaluates the paper-facing constant-global-DoA and affine-Doppler
  % CP model: fd(t)=fd0+fdRate*t with phase fd0*t+0.5*fdRate*t^2.
  % Framewise steering only maps the same global usrLla to each satellite frame.
  crbOptBase.steeringRefFrameIdx = sceneSeq.refFrameIdx;

  crbOptK = crbOptBase;
  crbOptK.fdRateMode = 'known';
  [crbK, auxK] = crbPilotMfDoaDoppler(sceneSeq, sharedData.pilotWave, ...
    config.carrierFreq, sharedData.waveInfo.sampleRate, config.usrLla(1:2, 1), ...
    truthDyn.fdRefFit, truthDyn.fdRateFit, pathGain, noiseVar, crbOptK);
  metricK = localExtractMfCrbMetrics(crbK, auxK);

  crbOptU = crbOptBase;
  crbOptU.fdRateMode = 'unknown';
  [crbU, auxU] = crbPilotMfDoaDoppler(sceneSeq, sharedData.pilotWave, ...
    config.carrierFreq, sharedData.waveInfo.sampleRate, config.usrLla(1:2, 1), ...
    truthDyn.fdRefFit, truthDyn.fdRateFit, pathGain, noiseVar, crbOptU);
  metricU = localExtractMfCrbMetrics(crbU, auxU);

  row.angleCrbDegCpK = metricK.angleTraceDeg;
  row.fdRefCrbHzCpK = metricK.fdRefStdHz;
  row.efimCondCpK = metricK.efimCond;
  row.doaBlockCondCpK = metricK.doaBlockCond;
  row.fdRefScalarStdHzCpK = metricK.fdRefScalarStdHz;
  row.fdRefCouplingFractionCpK = metricK.fdRefCouplingFraction;
  row.fdRefDoaCanonicalCorrCpK = metricK.fdRefDoaCanonicalCorr;
  row.fdRefStdDoaKnownHzCpK = metricK.fdRefStdDoaKnownHz;

  row.angleCrbDegCpU = metricU.angleTraceDeg;
  row.fdRefCrbHzCpU = metricU.fdRefStdHz;
  row.efimCondCpU = metricU.efimCond;
  row.doaBlockCondCpU = metricU.doaBlockCond;
  row.fdRefScalarStdHzCpU = metricU.fdRefScalarStdHz;
  row.fdRefCouplingFractionCpU = metricU.fdRefCouplingFraction;
  row.fdRefDoaCanonicalCorrCpU = metricU.fdRefDoaCanonicalCorr;
  row.fdRefStdDoaKnownHzCpU = metricU.fdRefStdDoaKnownHz;
  row.unknownOverKnownFdRefPct = 100 * (localSafeRatio(row.fdRefCrbHzCpU, row.fdRefCrbHzCpK) - 1);
  row.status = "ok";
catch ME
  row.status = "error";
  row.skipReason = "error";
  row.errorMessage = string(ME.message);
end

taskOut = row;
end

function row = localAttachSceneGeometry(row, sceneSeq, sharedData, trialSatIdxGlobal)
%LOCALATTACHSCENEGEOMETRY Attach access, elevation, separation and drift metrics.

elevMat = localExtractElevationMatrix(sceneSeq);
offAxisMat = localExtractOffAxisMatrix(sceneSeq);
row.elevRefDeg = localSafeMatrixElem(elevMat, row.refSatIdxLocal, sceneSeq.refFrameIdx);
row.elevMinDeg = min(elevMat(:), [], 'omitnan');
row.elevMeanDeg = mean(elevMat(:), 'omitnan');
row.offAxisMaxDeg = max(offAxisMat(:), [], 'omitnan');
row.angularSeparationMinDeg = localMinPairAngularSeparationDeg(sceneSeq);
if numel(trialSatIdxGlobal) >= 2
  row.candElevRefDeg = localSafeMatrixElem(elevMat, numel(trialSatIdxGlobal), sceneSeq.refFrameIdx);
  row.candElevMinDeg = min(elevMat(end, :), [], 'omitnan');
  row.candOffAxisMaxDeg = max(offAxisMat(end, :), [], 'omitnan');
  row.refCandAngularSeparationDeg = localRefCandidateSeparationDeg(sceneSeq, row.refSatIdxLocal, numel(trialSatIdxGlobal));
end
[row.localDirDriftMaxDeg, row.localDirDriftRmsDeg] = localCalcLocalDirDrift(sceneSeq);
row.baselineProxyKm = localBaselineProxyKm(sceneSeq);
end

function row = localAttachDynamicTruth(row, truthDyn, sceneSeq)
%LOCALATTACHDYNAMICTRUTH Attach Doppler trajectory and phase residual metrics.

refFrameIdx = sceneSeq.refFrameIdx;
row.fdRefTrueHz = truthDyn.fdRefSeries(refFrameIdx);
row.fdRefFitHz = truthDyn.fdRefFit;
row.fdRateFitHzPerSec = truthDyn.fdRateFit;
row.fdRefDeltaWinHz = truthDyn.deltaFdWin;
row.quadPhaseWinRad = truthDyn.quadPhaseWin;
row.fdLineResidualRmsHz = sqrt(mean(truthDyn.fdRefResidual(:) .^ 2));
row.fdLineResidualMaxHz = max(abs(truthDyn.fdRefResidual(:)));
row.deltaFdMaxAbsHz = max(abs(truthDyn.deltaFdSeries(:)));
row.deltaFdFitMaxAbsHz = max(abs(truthDyn.deltaFdFit(:)));
row.deltaFdRateMaxAbsHzPerSec = max(abs(truthDyn.deltaFdRate(:)));
if size(truthDyn.deltaFdSeries, 1) >= 2
  row.deltaFdRefCandHz = truthDyn.deltaFdSeries(end, refFrameIdx);
  row.deltaFdRateCandHzPerSec = truthDyn.deltaFdRate(end);
end
[phaseRms, phaseMax] = localCalcDopplerFitPhaseResidual(truthDyn);
row.phaseResidualRmsRad = phaseRms;
row.phaseResidualMaxRad = phaseMax;
end

function metric = localExtractMfCrbMetrics(crb, aux)
%LOCALEXTRACTMFCRBMETRICS Extract compact MF CRB, EFIM and coupling metrics.

metric = struct();
metric.angleTraceDeg = localCrbAngleTraceStdDeg(crb);
metric.fdRefStdHz = localCrbStd(crb, 3);
efim = localGetFieldOrDefault(aux, 'fimInterest', []);
if isempty(efim) || ~isnumeric(efim)
  metric.efimCond = NaN;
  metric.doaBlockCond = NaN;
  metric.fdRefScalarStdHz = NaN;
  metric.fdRefCouplingFraction = NaN;
  metric.fdRefDoaCanonicalCorr = NaN;
  metric.fdRefStdDoaKnownHz = NaN;
  return;
end
efim = (efim + efim') / 2;
metric.efimCond = localSafeCond(efim);
metric.doaBlockCond = localSafeCond(efim(1:2, 1:2));
fdInfo = localSafeMatrixElem(efim, 3, 3);
metric.fdRefScalarStdHz = localScalarInfoStd(fdInfo);
metric.fdRefStdDoaKnownHz = metric.fdRefScalarStdHz;
if size(efim, 1) >= 3
  doaBlock = efim(1:2, 1:2);
  couplingVec = efim(1:2, 3);
  doaInv = localSafeInv(doaBlock);
  if ~isempty(doaInv) && isfinite(fdInfo)
    couplingInfo = real(couplingVec' * doaInv * couplingVec);
    metric.fdRefCouplingFraction = localSafeRatio(couplingInfo, fdInfo);
    metric.fdRefDoaCanonicalCorr = sqrt(max(0, min(1, metric.fdRefCouplingFraction)));
  else
    metric.fdRefCouplingFraction = NaN;
    metric.fdRefDoaCanonicalCorr = NaN;
  end
else
  metric.fdRefCouplingFraction = NaN;
  metric.fdRefDoaCanonicalCorr = NaN;
end
end

function referencePreselectTable = localBuildReferencePreselectTable(availableSatIdxGlobal, satAccessRef, config)
%LOCALBUILDREFERENCEPRESELECTTABLE Build lightweight reference candidate metadata.

satIdx = reshape(double(availableSatIdxGlobal), [], 1);
numSat = numel(satIdx);
satName = strings(numSat, 1);
elevRefDeg = NaN(numSat, 1);
offAxisRefDeg = NaN(numSat, 1);
for iSat = 1:numSat
  idx = satIdx(iSat);
  satName(iSat) = localExtractSatName(satAccessRef, idx);
  elevRefDeg(iSat) = localSafeMatrixElem(satAccessRef.usrElevationDeg, idx, 1);
  offAxisRefDeg(iSat) = localSafeMatrixElem(satAccessRef.satOffAxisDeg, idx, 1);
end
referencePreselectTable = table(satIdx, satName, elevRefDeg, offAxisRefDeg, ...
  'VariableNames', {'refSatIdxGlobal', 'satName', 'elevRefDeg', 'offAxisRefDeg'});
referencePreselectTable = sortrows(referencePreselectTable, {'elevRefDeg'}, {'descend'});
referencePreselectTable.preselectRank = (1:height(referencePreselectTable)).';
referencePreselectTable = movevars(referencePreselectTable, 'preselectRank', 'Before', 'refSatIdxGlobal');
if height(referencePreselectTable) > config.maxReferenceCandidateEval
  referencePreselectTable = referencePreselectTable(1:config.maxReferenceCandidateEval, :);
end
end

function taskList = localBuildReferenceEvalTaskList(referencePreselectTable, config)
%LOCALBUILDREFERENCEEVALTASKLIST Build reference-evaluation satellite-set tasks.

if string(config.refSelectionMode) == "manual"
  satList = double(config.manualRefSatIdxGlobal);
else
  satList = double(referencePreselectTable.refSatIdxGlobal(:).');
end
satList = unique(satList, 'stable');
taskList = repmat(localEmptyTask(), numel(satList), 1);
for iTask = 1:numel(satList)
  taskList(iTask).taskId = iTask;
  taskList(iTask).taskType = "ref";
  taskList(iTask).refSatIdxGlobal = satList(iTask);
  taskList(iTask).baseSelectedSatIdxGlobal = satList(iTask);
  taskList(iTask).trialSatIdxGlobal = satList(iTask);
end
end

function selectedRow = localSelectReferenceRow(referenceCandidateTable, config)
%LOCALSELECTREFERENCEROW Select the reference satellite from evaluated rows.

if isempty(referenceCandidateTable)
  error('scanMfPairGeometryCrbAudit:EmptyReferenceCandidateTable', ...
    'No reference candidate was evaluated.');
end
validTable = referenceCandidateTable(referenceCandidateTable.status == "ok" & ...
  referenceCandidateTable.accessAllWindow, :);
if isempty(validTable)
  error('scanMfPairGeometryCrbAudit:NoAllAccessReference', ...
    'No reference candidate has all-window access.');
end
switch string(config.refSelectionMode)
  case "manual"
    matchIdx = find(validTable.refSatIdxGlobal == config.manualRefSatIdxGlobal, 1);
    if isempty(matchIdx)
      error('scanMfPairGeometryCrbAudit:ManualReferenceInvalid', ...
        'Manual reference satellite %d is not an all-access valid candidate.', ...
        config.manualRefSatIdxGlobal);
    end
    selectedRow = validTable(matchIdx, :);
  case "singleSatBestScore"
    validTable = sortrows(validTable, {'refScore', 'elevMinDeg'}, {'descend', 'descend'});
    selectedRow = validTable(1, :);
  case "bestScore"
    validTable = sortrows(validTable, {'refScore', 'elevMinDeg'}, {'descend', 'descend'});
    selectedRow = validTable(1, :);
  otherwise
    validTable = sortrows(validTable, {'elevRefDeg', 'elevMinDeg'}, {'descend', 'descend'});
    selectedRow = validTable(1, :);
end
end

function [refScoreTable, allPairRankTable, checkpointRunStateList] = ...
  localBuildCooperativeReferenceAudit(referenceCandidateTable, availableSatIdxGlobal, ...
    sharedData, config, repoRoot, checkpointRunStateList)
%LOCALBUILDCOOPERATIVEREFERENCEAUDIT Score references by their best cooperative pairs.

validRefTable = referenceCandidateTable(referenceCandidateTable.status == "ok" & ...
  referenceCandidateTable.accessAllWindow, :);
if isempty(validRefTable)
  refScoreTable = localEmptyRefCooperativeScoreTable();
  allPairRankTable = table();
  return;
end
validRefTable = sortrows(validRefTable, {'refScore', 'elevRefDeg'}, {'descend', 'descend'});
maxRef = min(height(validRefTable), max(1, double(config.maxCooperativeReferenceEval)));
validRefTable = validRefTable(1:maxRef, :);

[coopTaskList, refTaskCount] = localBuildCooperativePairTaskList( ...
  validRefTable, availableSatIdxGlobal, sharedData.satAccessRef, sharedData.satStateRef, config);
useParforStage = localCanUseParfor(config.useParfor, config.minTaskForParfor, numel(coopTaskList));
localPrintRuntimeLog("refCoop", sprintf('Auditing %d reference(s) with %d cooperative pair task(s).', ...
  height(validRefTable), numel(coopTaskList)));
[pairOutCell, stageRunState] = localRunTaskBatchWithState(coopTaskList, sharedData, ...
  "coopAll", config, repoRoot, useParforStage);
if ~isempty(stageRunState)
  checkpointRunStateList = [checkpointRunStateList; stageRunState(:)]; %#ok<AGROW>
end

scoreRows = repmat(localEmptyRefCooperativeScoreRow(), height(validRefTable), 1);
pairRankCell = cell(height(validRefTable), 1);
pairCandidateTableAll = localBuildSetResultTable(pairOutCell);
for iRef = 1:height(validRefTable)
  refSat = double(validRefTable.refSatIdxGlobal(iRef));
  localPrintRuntimeLog("refCoop", sprintf('Summarizing reference %d (%d/%d, %d task(s)).', ...
    refSat, iRef, height(validRefTable), refTaskCount(iRef)));
  baselineResult = localFillSetResultDefaults(table2struct(validRefTable(iRef, :)));
  pairCandidateTable = pairCandidateTableAll(pairCandidateTableAll.refSatIdxGlobal == refSat, :);
  pairCandidateTable = localAttachGainAndScore(pairCandidateTable, baselineResult, config);
  pairRankTable = localBuildPairRankTable(pairCandidateTable);
  pairRankCell{iRef} = pairRankTable;
  scoreRows(iRef) = localBuildRefCooperativeScoreRow(validRefTable(iRef, :), pairRankTable);
end
refScoreTable = struct2table(scoreRows);
refScoreTable = sortrows(refScoreTable, {'refCoopScore', 'bestPairScore', 'refElevRefDeg'}, ...
  {'descend', 'descend', 'descend'});
refScoreTable.rankByCoopScore = (1:height(refScoreTable)).';
refScoreTable = movevars(refScoreTable, 'rankByCoopScore', 'Before', 'refSatIdxGlobal');

if isempty(pairRankCell) || all(cellfun(@isempty, pairRankCell))
  allPairRankTable = table();
else
  allPairRankTable = vertcat(pairRankCell{~cellfun(@isempty, pairRankCell)});
end
end

function selectedRow = localSelectReferenceRowFromCooperative(referenceCandidateTable, refScoreTable, config)
%LOCALSELECTREFERENCEROWFROMCOOPERATIVE Select reference from cooperative score table.

if isempty(refScoreTable)
  selectedRow = localSelectReferenceRow(referenceCandidateTable, config);
  return;
end
validScoreTable = refScoreTable(refScoreTable.status == "ok" & isfinite(refScoreTable.refCoopScore), :);
if isempty(validScoreTable)
  selectedRow = localSelectReferenceRow(referenceCandidateTable, config);
  return;
end
validScoreTable = sortrows(validScoreTable, {'refCoopScore', 'bestPairScore', 'refElevRefDeg'}, ...
  {'descend', 'descend', 'descend'});
selectedRef = double(validScoreTable.refSatIdxGlobal(1));
matchIdx = find(referenceCandidateTable.refSatIdxGlobal == selectedRef, 1);
if isempty(matchIdx)
  error('scanMfPairGeometryCrbAudit:MissingCooperativeReference', ...
    'Selected cooperative reference %d is missing from referenceCandidateTable.', selectedRef);
end
selectedRow = referenceCandidateTable(matchIdx, :);
end

function row = localBuildRefCooperativeScoreRow(refRow, pairRankTable)
%LOCALBUILDREFCOOPERATIVESCOREROW Summarize one reference by its cooperative pair options.

row = localEmptyRefCooperativeScoreRow();
row.refSatIdxGlobal = double(refRow.refSatIdxGlobal);
row.refSatName = string(refRow.selectedSatNameStr);
row.refElevRefDeg = double(refRow.elevRefDeg);
row.refElevMinDeg = double(refRow.elevMinDeg);
row.refScore = double(refRow.refScore);
if isempty(pairRankTable)
  row.status = "no-pair";
  row.refCoopClass = "no-pair";
  return;
end
validPairTable = pairRankTable(pairRankTable.status == "ok" & pairRankTable.accessAllWindow & ...
  isfinite(pairRankTable.pairScore), :);
if isempty(validPairTable)
  row.status = "no-valid-pair";
  row.refCoopClass = "no-valid-pair";
  return;
end
validPairTable = sortrows(validPairTable, {'pairScore', 'angleCrbGainCpK'}, {'descend', 'descend'});
bestPair = validPairTable(1, :);
topCount = min(5, height(validPairTable));
topTable = validPairTable(1:topCount, :);
row.bestPairSatIdxGlobal = double(bestPair.candSatIdxGlobal);
row.bestPairScore = double(bestPair.pairScore);
row.bestPairClass = string(bestPair.pairClass);
row.bestPairAngleCrbDegCpK = double(bestPair.angleCrbDegCpK);
row.bestPairFdRefCrbHzCpU = double(bestPair.fdRefCrbHzCpU);
row.bestPairEfimCondCpU = double(bestPair.efimCondCpU);
row.bestPairFdRefDoaCanonicalCorrCpU = double(bestPair.fdRefDoaCanonicalCorrCpU);
row.bestPairCouplingFractionCpU = double(bestPair.fdRefCouplingFractionCpU);
row.medianTop5PairScore = median(double(topTable.pairScore), 'omitnan');
row.medianTop5EfimCondCpU = median(double(topTable.efimCondCpU), 'omitnan');
row.medianTop5CouplingFractionCpU = median(double(topTable.fdRefCouplingFractionCpU), 'omitnan');
row.numPair = height(validPairTable);
row.numHealthyPair = sum(validPairTable.pairClass == "healthy-gain");
row.numHardCoupledPair = sum(validPairTable.pairClass == "hard-coupled");
row.healthyPairRate = localSafeRatio(row.numHealthyPair, row.numPair);
row.refCoopScore = row.bestPairScore + 0.5 * row.medianTop5PairScore + ...
  0.2 * log1p(row.numHealthyPair) - 0.2 * log10(max(row.bestPairEfimCondCpU, 1));
row.refCoopClass = localClassifyReferenceCoop(row);
row.status = "ok";
end

function coopClass = localClassifyReferenceCoop(row)
%LOCALCLASSIFYREFERENCECOOP Classify a reference by cooperative pair health.
if row.status ~= "ok"
  coopClass = string(row.status);
elseif row.numHealthyPair > 0
  coopClass = "good-cooperative-ref";
elseif row.bestPairClass == "hard-coupled" && row.refScore > 0
  coopClass = "single-sat-good-but-coop-hard";
elseif row.bestPairClass == "hard-coupled"
  coopClass = "hard-coupled-ref";
else
  coopClass = "weak-reference";
end
end

function tableOut = localAttachReferenceScores(tableIn, config)
%LOCALATTACHREFERENCESCORES Attach a lightweight reference suitability score.

tableOut = tableIn;
if isempty(tableOut)
  return;
end
angleInfo = 1 ./ max(tableOut.angleCrbDegCpK, eps);
fdInfo = 1 ./ max(tableOut.fdRefCrbHzCpU, eps);
condPenalty = log10(max(tableOut.efimCondCpU, 1));
dynPenalty = abs(tableOut.quadPhaseWinRad - median(tableOut.quadPhaseWinRad, 'omitnan'));
tableOut.refScore = localSafeZScore(tableOut.elevMinDeg) + ...
  0.5 * localSafeZScore(log(angleInfo)) + ...
  0.25 * localSafeZScore(log(fdInfo)) - ...
  0.25 * localSafeZScore(condPenalty) - ...
  0.10 * localSafeZScore(dynPenalty);
end

function candidateTable = localBuildSecondSatPreselectTable(availableSatIdxGlobal, refSatIdxGlobal, satAccessRef, satStateRef, config)
%LOCALBUILDSECONDSATPRESELECTTABLE Build second-satellite preselection metadata.

candidateSatIdx = setdiff(reshape(double(availableSatIdxGlobal), 1, []), double(refSatIdxGlobal), 'stable');
forcedList = reshape(double(config.forcedCandidateSatIdxGlobal), 1, []);
forcedList = forcedList(isfinite(forcedList) & forcedList > 0 & forcedList ~= refSatIdxGlobal);
candidateSatIdx = unique([candidateSatIdx, forcedList], 'stable');
numCand = numel(candidateSatIdx);
rowList = repmat(localEmptyCandidatePreselectRow(), numCand, 1);
refLos = localLosUnitVectorAtRef(satStateRef, refSatIdxGlobal, config.usrLla);
for iCand = 1:numCand
  candIdx = candidateSatIdx(iCand);
  row = localEmptyCandidatePreselectRow();
  row.candSatIdxGlobal = candIdx;
  row.candSatName = localExtractSatName(satAccessRef, candIdx);
  row.candElevRefDeg = localSafeMatrixElem(satAccessRef.usrElevationDeg, candIdx, 1);
  row.candOffAxisRefDeg = localSafeMatrixElem(satAccessRef.satOffAxisDeg, candIdx, 1);
  candLos = localLosUnitVectorAtRef(satStateRef, candIdx, config.usrLla);
  row.refCandAngularSeparationDeg = localAngleBetweenUnitVecDeg(refLos, candLos);
  row.isForcedCandidate = ismember(candIdx, forcedList);
  row.geometryPreScore = row.candElevRefDeg + 0.25 * row.refCandAngularSeparationDeg - ...
    0.10 * row.candOffAxisRefDeg;
  rowList(iCand) = row;
end
candidateTable = struct2table(rowList);
if isempty(candidateTable)
  return;
end
candidateTable = sortrows(candidateTable, {'isForcedCandidate', 'geometryPreScore'}, {'descend', 'descend'});
candidateTable.preselectRank = (1:height(candidateTable)).';
candidateTable = movevars(candidateTable, 'preselectRank', 'Before', 'candSatIdxGlobal');
if height(candidateTable) > config.maxCandidateSecondSat
  forcedMask = candidateTable.isForcedCandidate;
  keepMask = false(height(candidateTable), 1);
  keepMask(1:min(config.maxCandidateSecondSat, height(candidateTable))) = true;
  keepMask = keepMask | forcedMask;
  candidateTable = candidateTable(keepMask, :);
  candidateTable = sortrows(candidateTable, {'isForcedCandidate', 'geometryPreScore'}, {'descend', 'descend'});
  candidateTable.preselectRank = (1:height(candidateTable)).';
end
end

function taskList = localBuildPairTaskList(candidatePreselectTable, refSatIdxGlobal, config)
%LOCALBUILDPAIRTASKLIST Build fixed-reference pair CRB tasks.

numTask = height(candidatePreselectTable);
taskList = repmat(localEmptyTask(), numTask, 1);
for iTask = 1:numTask
  candSat = double(candidatePreselectTable.candSatIdxGlobal(iTask));
  taskList(iTask).taskId = iTask;
  taskList(iTask).taskType = "pair";
  taskList(iTask).refSatIdxGlobal = refSatIdxGlobal;
  taskList(iTask).baseSelectedSatIdxGlobal = refSatIdxGlobal;
  taskList(iTask).trialSatIdxGlobal = [refSatIdxGlobal, candSat];
end
end

function [taskList, refTaskCount] = localBuildCooperativePairTaskList( ...
  validRefTable, availableSatIdxGlobal, satAccessRef, satStateRef, config)
%LOCALBUILDCOOPERATIVEPAIRTASKLIST Build one geo-level cooperative pair task batch.

numRef = height(validRefTable);
refTaskCell = cell(numRef, 1);
refTaskCount = zeros(numRef, 1);
for iRef = 1:numRef
  refSat = double(validRefTable.refSatIdxGlobal(iRef));
  candidatePreselectTable = localBuildSecondSatPreselectTable( ...
    availableSatIdxGlobal, refSat, satAccessRef, satStateRef, config);
  taskListRef = localBuildPairTaskList(candidatePreselectTable, refSat, config);
  refTaskCell{iRef} = taskListRef(:);
  refTaskCount(iRef) = numel(taskListRef);
end

if all(refTaskCount == 0)
  taskList = repmat(localEmptyTask(), 0, 1);
  return;
end
taskList = vertcat(refTaskCell{refTaskCount > 0});
for iTask = 1:numel(taskList)
  taskList(iTask).taskId = iTask;
end
end

function [greedyStepTable, greedyCandidateTable, satSetSummaryTable, checkpointRunStateList] = ...
  localRunGreedySatelliteSelection(pairRankTable, refSatIdxGlobal, sharedData, config, repoRoot, checkpointRunStateList)
%LOCALRUNGREEDYSATELLITESELECTION Greedily grow a nested satellite set.

maxSatCount = max(config.satCountList);
selectedSet = double(refSatIdxGlobal);
selectedResultList = localEmptySetResultRow();
selectedResultList(1) = sharedData.baselineMetrics;
selectedResultList(1).taskType = "selected";
selectedResultList(1).greedyStep = 1;
selectedResultList(1).addedSatIdxGlobal = refSatIdxGlobal;
selectedResultList(1).selectionClass = "reference";

allCandidateRowsCell = cell(0, 1);
stepRows = repmat(localEmptyGreedyStepRow(), maxSatCount, 1);
stepRows(1) = localBuildGreedyStepRow(1, refSatIdxGlobal, selectedSet, selectedResultList(1), NaN, "reference");

candidatePool = pairRankTable(pairRankTable.status == "ok" & pairRankTable.accessAllWindow, :);
if isempty(candidatePool)
  greedyStepTable = struct2table(stepRows(1));
  greedyCandidateTable = table();
  satSetSummaryTable = localBuildSatSetSummaryTable(selectedResultList, config);
  return;
end
candidatePool = sortrows(candidatePool, {'pairScore', 'angleCrbGainCpK'}, {'descend', 'descend'});
allCandidateSat = double(candidatePool.candSatIdxGlobal(:).');

for targetCount = 2:maxSatCount
  remainingSat = setdiff(allCandidateSat, selectedSet, 'stable');
  if isempty(remainingSat)
    break;
  end
  remainingSat = remainingSat(1:min(numel(remainingSat), config.maxGreedyCandidatePerStep));
  taskList = repmat(localEmptyTask(), numel(remainingSat), 1);
  for iTask = 1:numel(remainingSat)
    taskList(iTask).taskId = iTask;
    taskList(iTask).taskType = "greedy" + string(targetCount);
    taskList(iTask).refSatIdxGlobal = refSatIdxGlobal;
    taskList(iTask).baseSelectedSatIdxGlobal = selectedSet;
    taskList(iTask).trialSatIdxGlobal = [selectedSet, remainingSat(iTask)];
  end
  useParforStage = localCanUseParfor(config.useParfor, config.minTaskForParfor, numel(taskList));
  [outCell, stageRunState] = localRunTaskBatchWithState(taskList, sharedData, ...
    "greedy" + string(targetCount), config, repoRoot, useParforStage);
  if ~isempty(stageRunState)
    checkpointRunStateList = [checkpointRunStateList; stageRunState(:)]; %#ok<AGROW>
  end
  candidateTable = localBuildSetResultTable(outCell);
  candidateTable = localAttachGainAndScore(candidateTable, selectedResultList(end), config);
  candidateTable.greedyTargetCount = repmat(targetCount, height(candidateTable), 1);
  allCandidateRowsCell{end + 1, 1} = candidateTable; %#ok<AGROW>
  validTable = candidateTable(candidateTable.status == "ok" & candidateTable.accessAllWindow, :);
  if isempty(validTable)
    break;
  end
  validTable = sortrows(validTable, {'pairScore', 'angleCrbGainCpK'}, {'descend', 'descend'});
  chosen = validTable(1, :);
  addedSat = chosen.candSatIdxGlobal;
  selectedSet = [selectedSet, addedSat]; %#ok<AGROW>
  chosenResult = localFillSetResultDefaults(table2struct(chosen));
  chosenResult.taskType = "selected";
  chosenResult.greedyStep = targetCount;
  chosenResult.addedSatIdxGlobal = addedSat;
  chosenResult.selectionClass = chosen.pairClass;
  selectedResultList(end + 1, 1) = chosenResult; %#ok<AGROW>
  stepRows(targetCount) = localBuildGreedyStepRow(targetCount, addedSat, selectedSet, ...
    chosenResult, chosen.pairScore, chosen.pairClass);
end

greedyStepTable = struct2table(stepRows(1:numel(selectedResultList)));
if isempty(allCandidateRowsCell)
  greedyCandidateTable = table();
else
  greedyCandidateTable = vertcat(allCandidateRowsCell{:});
end
satSetSummaryTable = localBuildSatSetSummaryTable(selectedResultList, config);
end

function row = localBuildGreedyStepRow(step, addedSat, selectedSet, result, scoreValue, selectionClass)
%LOCALBUILDGREEDYSTEPROW Build one greedy selection summary row.

row = localEmptyGreedyStepRow();
row.step = double(step);
row.numSat = double(numel(selectedSet));
row.addedSatIdxGlobal = double(addedSat);
row.selectedSatIdxGlobalStr = string(mat2str(selectedSet));
row.angleCrbDegCpK = result.angleCrbDegCpK;
row.fdRefCrbHzCpK = result.fdRefCrbHzCpK;
row.angleCrbDegCpU = result.angleCrbDegCpU;
row.fdRefCrbHzCpU = result.fdRefCrbHzCpU;
row.efimCondCpU = result.efimCondCpU;
row.fdRefDoaCanonicalCorrCpU = result.fdRefDoaCanonicalCorrCpU;
row.fdRefCouplingFractionCpU = result.fdRefCouplingFractionCpU;
row.pairScore = scoreValue;
row.selectionClass = string(selectionClass);
end

function satSetSummaryTable = localBuildSatSetSummaryTable(resultList, config)
%LOCALBUILDSATSETSUMMARYTABLE Build final selected-set summary rows.

if isempty(resultList)
  satSetSummaryTable = struct2table(localEmptySetSummaryRow());
  satSetSummaryTable(1, :) = [];
  return;
end
rowList = repmat(localEmptySetSummaryRow(), numel(resultList), 1);
for iRow = 1:numel(resultList)
  result = resultList(iRow);
  row = localEmptySetSummaryRow();
  row.numSat = result.numSat;
  row.selectedSatIdxGlobalStr = result.selectedSatIdxGlobalStr;
  row.angleCrbDegCpK = result.angleCrbDegCpK;
  row.fdRefCrbHzCpK = result.fdRefCrbHzCpK;
  row.angleCrbDegCpU = result.angleCrbDegCpU;
  row.fdRefCrbHzCpU = result.fdRefCrbHzCpU;
  row.unknownOverKnownFdRefPct = result.unknownOverKnownFdRefPct;
  row.efimCondCpK = result.efimCondCpK;
  row.efimCondCpU = result.efimCondCpU;
  row.fdRefDoaCanonicalCorrCpU = result.fdRefDoaCanonicalCorrCpU;
  row.fdRefCouplingFractionCpU = result.fdRefCouplingFractionCpU;
  row.phaseResidualRmsRad = result.phaseResidualRmsRad;
  row.pairScore = localGetFieldOrDefault(result, 'pairScore', NaN);
  row.pairClass = localGetFieldOrDefault(result, 'pairClass', "");
  rowList(iRow) = row;
end
satSetSummaryTable = struct2table(rowList);
targetMask = ismember(satSetSummaryTable.numSat, reshape(double(config.satCountList), [], 1));
satSetSummaryTable = satSetSummaryTable(targetMask, :);
end

function row = localEmptyTask()
%LOCALEMPTYTASK Build an empty task struct.
row = struct();
row.taskId = 0;
row.taskType = "";
row.refSatIdxGlobal = NaN;
row.baseSelectedSatIdxGlobal = [];
row.trialSatIdxGlobal = [];
end

function row = localEmptySetResultRow()
%LOCALEMPTYSETRESULTROW Return one set-result row with scalar fields.
row = struct();
row.taskId = NaN;
row.taskType = "";
row.refSatIdxGlobal = NaN;
row.refSatIdxLocal = NaN;
row.candSatIdxGlobal = NaN;
row.candSatName = "";
row.numSat = NaN;
row.selectedSatIdxGlobalStr = "";
row.selectedSatNameStr = "";
row.baseSelectedSatIdxGlobalStr = "";
row.accessAllWindow = false;
row.elevRefDeg = NaN;
row.elevMinDeg = NaN;
row.elevMeanDeg = NaN;
row.offAxisMaxDeg = NaN;
row.candElevRefDeg = NaN;
row.candElevMinDeg = NaN;
row.candOffAxisMaxDeg = NaN;
row.angularSeparationMinDeg = NaN;
row.refCandAngularSeparationDeg = NaN;
row.baselineProxyKm = NaN;
row.localDirDriftMaxDeg = NaN;
row.localDirDriftRmsDeg = NaN;
row.fdRefTrueHz = NaN;
row.fdRefFitHz = NaN;
row.fdRateFitHzPerSec = NaN;
row.fdRefDeltaWinHz = NaN;
row.quadPhaseWinRad = NaN;
row.fdLineResidualRmsHz = NaN;
row.fdLineResidualMaxHz = NaN;
row.phaseResidualRmsRad = NaN;
row.phaseResidualMaxRad = NaN;
row.deltaFdRefCandHz = NaN;
row.deltaFdRateCandHzPerSec = NaN;
row.deltaFdMaxAbsHz = NaN;
row.deltaFdFitMaxAbsHz = NaN;
row.deltaFdRateMaxAbsHzPerSec = NaN;
row.angleCrbDegCpK = NaN;
row.fdRefCrbHzCpK = NaN;
row.efimCondCpK = NaN;
row.doaBlockCondCpK = NaN;
row.fdRefScalarStdHzCpK = NaN;
row.fdRefCouplingFractionCpK = NaN;
row.fdRefDoaCanonicalCorrCpK = NaN;
row.fdRefStdDoaKnownHzCpK = NaN;
row.angleCrbDegCpU = NaN;
row.fdRefCrbHzCpU = NaN;
row.efimCondCpU = NaN;
row.doaBlockCondCpU = NaN;
row.fdRefScalarStdHzCpU = NaN;
row.fdRefCouplingFractionCpU = NaN;
row.fdRefDoaCanonicalCorrCpU = NaN;
row.fdRefStdDoaKnownHzCpU = NaN;
row.unknownOverKnownFdRefPct = NaN;
row.angleCrbGainCpK = NaN;
row.fdRefCrbGainCpK = NaN;
row.angleCrbGainCpU = NaN;
row.fdRefCrbGainCpU = NaN;
row.pairScore = NaN;
row.pairClass = "";
row.greedyStep = NaN;
row.addedSatIdxGlobal = NaN;
row.selectionClass = "";
row.status = "pending";
row.skipReason = "";
row.errorMessage = "";
end

function row = localEmptyCandidatePreselectRow()
%LOCALEMPTYCANDIDATEPRESELECTROW Build second-satellite preselection row.
row = struct();
row.candSatIdxGlobal = NaN;
row.candSatName = "";
row.candElevRefDeg = NaN;
row.candOffAxisRefDeg = NaN;
row.refCandAngularSeparationDeg = NaN;
row.geometryPreScore = NaN;
row.isForcedCandidate = false;
end

function row = localEmptyGreedyStepRow()
%LOCALEMPTYGREEDYSTEPROW Build one greedy step row.
row = struct();
row.step = NaN;
row.numSat = NaN;
row.addedSatIdxGlobal = NaN;
row.selectedSatIdxGlobalStr = "";
row.angleCrbDegCpK = NaN;
row.fdRefCrbHzCpK = NaN;
row.angleCrbDegCpU = NaN;
row.fdRefCrbHzCpU = NaN;
row.efimCondCpU = NaN;
row.fdRefDoaCanonicalCorrCpU = NaN;
row.fdRefCouplingFractionCpU = NaN;
row.pairScore = NaN;
row.selectionClass = "";
end

function row = localEmptySetSummaryRow()
%LOCALEMPTYSETSUMMARYROW Build one selected-set summary row.
row = struct();
row.numSat = NaN;
row.selectedSatIdxGlobalStr = "";
row.angleCrbDegCpK = NaN;
row.fdRefCrbHzCpK = NaN;
row.angleCrbDegCpU = NaN;
row.fdRefCrbHzCpU = NaN;
row.unknownOverKnownFdRefPct = NaN;
row.efimCondCpK = NaN;
row.efimCondCpU = NaN;
row.fdRefDoaCanonicalCorrCpU = NaN;
row.fdRefCouplingFractionCpU = NaN;
row.phaseResidualRmsRad = NaN;
row.pairScore = NaN;
row.pairClass = "";
end

function tableOut = localEmptyRefCooperativeScoreTable()
%LOCALEMPTYREFCOOPERATIVESCORETABLE Return an empty cooperative-reference table.
row = localEmptyRefCooperativeScoreRow();
tableOut = struct2table(row);
tableOut(1, :) = [];
end

function row = localEmptyRefCooperativeScoreRow()
%LOCALEMPTYREFCOOPERATIVESCOREROW Build one cooperative-reference score row.
row = struct();
row.refSatIdxGlobal = NaN;
row.refSatName = "";
row.refElevRefDeg = NaN;
row.refElevMinDeg = NaN;
row.refScore = NaN;
row.bestPairSatIdxGlobal = NaN;
row.bestPairScore = NaN;
row.bestPairClass = "";
row.bestPairAngleCrbDegCpK = NaN;
row.bestPairFdRefCrbHzCpU = NaN;
row.bestPairEfimCondCpU = NaN;
row.bestPairFdRefDoaCanonicalCorrCpU = NaN;
row.bestPairCouplingFractionCpU = NaN;
row.medianTop5PairScore = NaN;
row.medianTop5EfimCondCpU = NaN;
row.medianTop5CouplingFractionCpU = NaN;
row.numPair = 0;
row.numHealthyPair = 0;
row.numHardCoupledPair = 0;
row.healthyPairRate = NaN;
row.refCoopScore = NaN;
row.refCoopClass = "";
row.status = "pending";
end

function [outCell, runState] = localRunTaskBatchWithState(taskList, sharedData, stageName, config, repoRoot, useParforStage)
%LOCALRUNTASKBATCHWITHSTATE Run a task batch and return checkpoint state.

runState = struct([]);
if isempty(taskList)
  outCell = cell(0, 1);
  localPrintRuntimeLog(stageName, "No tasks to run.");
  return;
end
stageTic = tic;
modeText = localParforModeText(useParforStage);
localPrintRuntimeLog(stageName, sprintf('Starting %d task(s), mode=%s, checkpoint=%d, resume=%d.', ...
  numel(taskList), char(modeText), logical(config.checkpointEnable), logical(config.checkpointResume)));
if config.checkpointEnable
  checkpointOpt = localBuildCheckpointOpt(repoRoot, config.scanName, config, taskList, stageName, useParforStage);
  pendingCount = localCountPendingCheckpointTasks(checkpointOpt, numel(taskList));
  localPrintRuntimeLog(stageName, sprintf('Checkpoint run key: %s; pending %d/%d task(s).', ...
    char(checkpointOpt.runKey), pendingCount, numel(taskList)));
  progressTracker = localStartTaskProgress(pendingCount, ...
    sprintf('%s tasks (%s)', char(stageName), char(modeText)));
  progressCleanup = onCleanup(@() localFinishTaskProgress(progressTracker)); %#ok<NASGU>
  checkpointOpt.progressFcn = @(step) localAdvanceTaskProgress(progressTracker, step);
  runState = runPerfTaskGridWithCheckpoint(taskList, sharedData, ...
    @localRunSetCheckpointTask, checkpointOpt);
  outCell = runState.resultCell;
else
  outCell = cell(numel(taskList), 1);
  progressTracker = localStartTaskProgress(numel(taskList), ...
    sprintf('%s tasks (%s)', char(stageName), char(modeText)));
  progressCleanup = onCleanup(@() localFinishTaskProgress(progressTracker)); %#ok<NASGU>
  if useParforStage
    progressQueue = [];
    if progressTracker.isActive
      progressQueue = parallel.pool.DataQueue;
      afterEach(progressQueue, @(~) progressbar('advance'));
    end
    parfor iTask = 1:numel(taskList)
      outCell{iTask} = localRunOneSetTask(taskList(iTask), sharedData);
      if ~isempty(progressQueue)
        send(progressQueue, 0);
      end
    end
  else
    for iTask = 1:numel(taskList)
      outCell{iTask} = localRunOneSetTask(taskList(iTask), sharedData);
      localAdvanceTaskProgress(progressTracker, 1);
    end
  end
end
localPrintRuntimeLog(stageName, sprintf('Finished in %.2f s.', toc(stageTic)));
end

function outCell = localRunTaskBatch(taskList, sharedData, stageName, config, repoRoot, useParforStage)
%LOCALRUNTASKBATCH Run one task batch and return only task outputs.
[outCell, ~] = localRunTaskBatchWithState(taskList, sharedData, stageName, config, repoRoot, useParforStage);
end

function checkpointOpt = localBuildCheckpointOpt(repoRoot, scanName, config, taskList, stageName, useParforStage)
%LOCALBUILDCHECKPOINTOPT Build checkpoint options for one scan stage.
checkpointOpt = struct();
checkpointOpt.runName = string(scanName);
checkpointOpt.runKey = localBuildCheckpointRunKey(config, taskList, stageName);
checkpointOpt.outputRoot = fullfile(repoRoot, 'tmp');
checkpointOpt.useParfor = logical(useParforStage);
checkpointOpt.resume = logical(config.checkpointResume);
checkpointOpt.meta = localBuildCheckpointMeta(config, taskList, stageName);
checkpointOpt.progressFcn = [];
checkpointOpt.cleanupOnSuccess = false;
checkpointOpt.cleanupOpt = struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true);
end

function pendingCount = localCountPendingCheckpointTasks(checkpointOpt, numTask)
%LOCALCOUNTPENDINGCHECKPOINTTASKS Count unfinished tasks before resetting progressbar.
pendingCount = double(numTask);
if ~isfield(checkpointOpt, 'resume') || ~checkpointOpt.resume
  return;
end
taskDir = fullfile(string(checkpointOpt.outputRoot), string(checkpointOpt.runName), ...
  string(checkpointOpt.runKey), "task");
if ~isfolder(taskDir)
  return;
end
numDone = 0;
for iTask = 1:numTask
  if isfile(fullfile(taskDir, sprintf('task_%06d.mat', iTask)))
    numDone = numDone + 1;
  end
end
pendingCount = max(0, double(numTask - numDone));
end

function tracker = localStartTaskProgress(totalCount, progressLabel)
%LOCALSTARTTASKPROGRESS Start a client-side progressbar when available.
tracker = struct('isActive', false, 'totalCount', double(totalCount), 'label', string(progressLabel));
fprintf('  %s: %d pending task(s)\n', char(progressLabel), double(totalCount));
if totalCount <= 0
  return;
end
if exist('progressbar', 'file') ~= 2
  fprintf('  progressbar.m not found; continuing without progress bar.\n');
  return;
end
try
  progressbar('displaymode', 'replace');
  progressbar('minimalupdateinterval', 0.2);
  progressbar('reset', double(totalCount));
  tracker.isActive = true;
catch ME
  tracker.isActive = false;
  fprintf('  progressbar disabled: %s\n', ME.message);
end
end

function localAdvanceTaskProgress(tracker, step)
%LOCALADVANCETASKPROGRESS Advance an active client-side progressbar.
if ~isstruct(tracker) || ~isfield(tracker, 'isActive') || ~tracker.isActive
  return;
end
for iStep = 1:max(1, double(step))
  progressbar('advance');
end
end

function localFinishTaskProgress(tracker)
%LOCALFINISHTASKPROGRESS Close an active progressbar without touching scan results.
if isstruct(tracker) && isfield(tracker, 'isActive') && tracker.isActive
  pause(0.05);
  progressbar('end');
end
end

function modeText = localParforModeText(useParforStage)
%LOCALPARFORMODETEXT Return a compact execution-mode label.
if useParforStage
  modeText = "parfor";
else
  modeText = "serial";
end
end

function runKey = localBuildCheckpointRunKey(config, taskList, stageName)
%LOCALBUILDCHECKPOINTRUNKEY Build a stable checkpoint key.
if isempty(taskList)
  firstSat = NaN;
  lastSat = NaN;
else
  firstSat = localGetTaskCandidate(taskList(1));
  lastSat = localGetTaskCandidate(taskList(end));
end
if isempty(taskList)
  refToken = config.manualRefSatIdxGlobal;
else
  refToken = taskList(1).refSatIdxGlobal;
end
geoToken = '';
if isfield(config, 'geoCheckpointToken') && strlength(string(config.geoCheckpointToken)) > 0
  geoToken = char(string(config.geoCheckpointToken) + "_");
end
runKey = sprintf('%s%s_ref%d_nf%d_snr%g_n%d_first%d_last%d', ...
  geoToken, char(stageName), refToken, config.numFrame, config.snrDb, ...
  numel(taskList), firstSat, lastSat);
runKey = string(runKey);
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
runKey = replace(runKey, ' ', '');
runKey = replace(runKey, '+', '');
end

function meta = localBuildCheckpointMeta(config, taskList, stageName)
%LOCALBUILDCHECKPOINTMETA Store the task-defining checkpoint signature.
meta = struct();
meta.scanName = string(config.scanName);
meta.stageName = string(stageName);
meta.tleFileName = string(config.tleFileName);
meta.utc0 = config.utc0;
meta.usrLla = config.usrLla;
if isfield(config, 'geoIdx')
  meta.geoIdx = config.geoIdx;
end
if isfield(config, 'geoCheckpointToken')
  meta.geoCheckpointToken = string(config.geoCheckpointToken);
end
meta.numFrame = config.numFrame;
meta.frameIntvlSec = config.frameIntvlSec;
meta.signalModelTag = string(config.signalModelTag);
meta.numTask = numel(taskList);
end

function dirText = localDefaultCheckpointDir(repoRoot, scanName, config)
%LOCALDEFAULTCHECKPOINTDIR Return the parent checkpoint directory shown in header.
dirText = string(fullfile(repoRoot, 'tmp', char(scanName)));
if ~config.checkpointEnable
  dirText = "";
end
end

function cleanupReport = localCleanupCompletedCheckpointRuns(runStateList)
%LOCALCLEANUPCOMPLETEDCHECKPOINTRUNS Clean completed checkpoint runs with fixed report fields.
cleanupReport = localEmptyCheckpointCleanupReportTable();
for iRun = 1:numel(runStateList)
  runState = runStateList(iRun);
  if ~isstruct(runState) || ~isfield(runState, 'isComplete') || ~runState.isComplete
    continue;
  end
  row = localBuildCheckpointCleanupRow(runState);
  try
    reportTmp = cleanupPerfTaskGridCheckpoint(runState, ...
      struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));
    row.cleanupStatus = "ok";
    row.cleanupMessage = localSummarizeCleanupReport(reportTmp);
  catch ME
    row.cleanupStatus = "failed";
    row.cleanupMessage = string(ME.message);
    warning('scanMfPairGeometryCrbAudit:CheckpointCleanupFailed', ...
      'Checkpoint cleanup failed for %s: %s', char(row.runDir), ME.message);
  end
  cleanupReport = [cleanupReport; struct2table(row)]; %#ok<AGROW>
end
end

function cleanupReport = localEmptyCheckpointCleanupReportTable()
%LOCALEMPTYCHECKPOINTCLEANUPREPORTTABLE Return an empty fixed-field cleanup table.
row = localBuildCheckpointCleanupRow(struct());
cleanupReport = struct2table(row);
cleanupReport(1, :) = [];
end

function row = localBuildCheckpointCleanupRow(runState)
%LOCALBUILDCHECKPOINTCLEANUPROW Build one fixed-field checkpoint cleanup row.
row = struct();
row.runName = string(localGetFieldOrDefault(runState, 'runName', ""));
row.runKey = string(localGetFieldOrDefault(runState, 'runKey', ""));
row.runDir = string(localGetFieldOrDefault(runState, 'runDir', ""));
row.numTask = double(localGetFieldOrDefault(runState, 'numTask', 0));
row.numDone = double(localGetFieldOrDefault(runState, 'numDone', 0));
row.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
row.cleanupStatus = "not-run";
row.cleanupMessage = "";
end

function message = localSummarizeCleanupReport(reportTmp)
%LOCALSUMMARIZECLEANUPREPORT Return a compact text summary for variable-shape cleanup output.
if isempty(reportTmp)
  message = "empty cleanup report";
  return;
end
if isstruct(reportTmp)
  message = sprintf('cleanup report rows=%d fields=%s', numel(reportTmp), strjoin(fieldnames(reportTmp), ','));
elseif istable(reportTmp)
  message = sprintf('cleanup report rows=%d variables=%s', height(reportTmp), strjoin(reportTmp.Properties.VariableNames, ','));
else
  message = sprintf('cleanup report class=%s', class(reportTmp));
end
message = string(message);
end

function resultTable = localBuildSetResultTable(outCell)
%LOCALBUILDSETRESULTTABLE Convert task output cell to a table.
if isempty(outCell)
  resultTable = struct2table(localEmptySetResultRow());
  resultTable(1, :) = [];
  return;
end
rowList = repmat(localEmptySetResultRow(), numel(outCell), 1);
for iOut = 1:numel(outCell)
  if isstruct(outCell{iOut})
    rowList(iOut) = localFillSetResultDefaults(outCell{iOut});
  end
end
resultTable = struct2table(rowList);
end

function row = localFillSetResultDefaults(rowIn)
%LOCALFILLSETRESULTDEFAULTS Fill known fields in a task output struct.
row = localEmptySetResultRow();
fieldList = fieldnames(rowIn);
validFieldList = fieldnames(row);
for iField = 1:numel(fieldList)
  fieldName = fieldList{iField};
  if any(strcmp(validFieldList, fieldName))
    row.(fieldName) = rowIn.(fieldName);
  end
end
end

function tableOut = localAttachGainAndScore(tableIn, baselineResult, config)
%LOCALATTACHGAINANDSCORE Attach baseline-relative gains and diagnostic score.

tableOut = tableIn;
if isempty(tableOut)
  return;
end
if isstruct(baselineResult)
  b = baselineResult;
else
  b = localTableRowToStruct(baselineResult);
end
tableOut.angleCrbGainCpK = localSafeRatio(b.angleCrbDegCpK, tableOut.angleCrbDegCpK);
tableOut.fdRefCrbGainCpK = localSafeRatio(b.fdRefCrbHzCpK, tableOut.fdRefCrbHzCpK);
tableOut.angleCrbGainCpU = localSafeRatio(b.angleCrbDegCpU, tableOut.angleCrbDegCpU);
tableOut.fdRefCrbGainCpU = localSafeRatio(b.fdRefCrbHzCpU, tableOut.fdRefCrbHzCpU);
score = NaN(height(tableOut), 1);
pairClass = strings(height(tableOut), 1);
for iRow = 1:height(tableOut)
  row = tableOut(iRow, :);
  [score(iRow), pairClass(iRow)] = localScoreCandidateRow(row, config.scoreOpt);
end
tableOut.pairScore = score;
tableOut.pairClass = pairClass;
end

function [score, pairClass] = localScoreCandidateRow(row, scoreOpt)
%LOCALSCORECANDIDATEROW Score one candidate row for diagnostic ranking.
if row.status ~= "ok" || ~row.accessAllWindow
  score = -Inf;
  pairClass = "bad-access";
  return;
end
angleGain = max(double(row.angleCrbGainCpK), eps);
fdGain = max(double(row.fdRefCrbGainCpU), eps);
condPenalty = log10(max(double(row.efimCondCpU), 1));
coupling = min(max(double(row.fdRefDoaCanonicalCorrCpU), 0), 1);
couplingFraction = min(max(double(row.fdRefCouplingFractionCpU), 0), 1);
phaseResidual = max(double(row.phaseResidualRmsRad), 0);
score = scoreOpt.angleGainWeight * log(angleGain) + ...
  scoreOpt.fdGainWeight * log(fdGain) - ...
  scoreOpt.condPenaltyWeight * condPenalty - ...
  scoreOpt.couplingPenaltyWeight * coupling - ...
  scoreOpt.couplingFractionPenaltyWeight * couplingFraction - ...
  scoreOpt.phaseResidualPenaltyWeight * phaseResidual;
if angleGain < scoreOpt.lowGainTolerance && fdGain < scoreOpt.lowGainTolerance
  pairClass = "trivial-gain";
elseif condPenalty >= log10(scoreOpt.highCondThreshold) || coupling >= scoreOpt.highCouplingThreshold
  pairClass = "hard-coupled";
elseif phaseResidual >= scoreOpt.phaseResidualStressRad
  pairClass = "dynamic-stress";
else
  pairClass = "healthy-gain";
end
end

function rankTable = localBuildPairRankTable(pairCandidateTable)
%LOCALBUILDPAIRRANKTABLE Build sorted pair candidate table.
rankVarList = {'rankByCrbHealth', 'rankByAngleGain', 'rankByLowCoupling'};
existingRankVar = intersect(rankVarList, pairCandidateTable.Properties.VariableNames, 'stable');
if ~isempty(existingRankVar)
  pairCandidateTable = removevars(pairCandidateTable, existingRankVar);
end
rankMask = pairCandidateTable.status == "ok" & pairCandidateTable.accessAllWindow & ...
  isfinite(pairCandidateTable.pairScore);
rankTable = pairCandidateTable(rankMask, :);
if isempty(rankTable)
  return;
end
rankTable = sortrows(rankTable, {'pairScore', 'angleCrbGainCpK', 'efimCondCpU'}, ...
  {'descend', 'descend', 'ascend'});
rankTable.rankByCrbHealth = (1:height(rankTable)).';
rankTable = movevars(rankTable, 'rankByCrbHealth', 'Before', 'candSatIdxGlobal');
angleGainForSort = rankTable.angleCrbGainCpK;
angleGainForSort(~isfinite(angleGainForSort)) = -Inf;
[~, angleOrder] = sort(angleGainForSort, 'descend');
rankTable.rankByAngleGain = NaN(height(rankTable), 1);
rankTable.rankByAngleGain(angleOrder) = (1:height(rankTable)).';
couplingForSort = rankTable.fdRefDoaCanonicalCorrCpU;
couplingForSort(~isfinite(couplingForSort)) = Inf;
[~, couplingOrder] = sort(couplingForSort, 'ascend');
rankTable.rankByLowCoupling = NaN(height(rankTable), 1);
rankTable.rankByLowCoupling(couplingOrder) = (1:height(rankTable)).';
end

function previewTable = localBuildRankPreviewTable(pairRankTable, satSetSummaryTable, config)
%LOCALBUILDRANKPREVIEWTABLE Build a compact table for printed summary.
if isempty(pairRankTable)
  previewTable = table();
  return;
end
topCount = min(config.numTopReport, height(pairRankTable));
previewTable = pairRankTable(1:topCount, {'rankByCrbHealth', 'candSatIdxGlobal', ...
  'candSatName', 'pairClass', 'pairScore', 'angleCrbGainCpK', 'fdRefCrbGainCpU', ...
  'efimCondCpU', 'fdRefDoaCanonicalCorrCpU', 'fdRefCouplingFractionCpU', ...
  'phaseResidualRmsRad'});
end

function contextSummaryTable = localBuildContextSummaryTable(config, availableSatIdxGlobal, refSatIdxGlobal, baselineResult)
%LOCALBUILDCONTEXTSUMMARYTABLE Build scan context summary.
contextSummaryTable = table( ...
  string(config.tleFileName), ...
  string(config.signalModelTag), ...
  string(config.refSelectionMode), ...
  string(mat2str(config.usrLla(:).')), ...
  string(config.utc0), ...
  config.numFrame, ...
  config.frameIntvlSec, ...
  numel(availableSatIdxGlobal), ...
  refSatIdxGlobal, ...
  baselineResult.elevRefDeg, ...
  baselineResult.elevMinDeg, ...
  baselineResult.angleCrbDegCpK, ...
  baselineResult.fdRefCrbHzCpU, ...
  'VariableNames', { ...
    'tleFileName', 'signalModelTag', 'refSelectionMode', 'usrLla', 'utc0', 'numFrame', 'frameIntvlSec', ...
    'numAvailableAtRef', 'selectedRefSatIdxGlobal', 'refElevRefDeg', ...
    'refElevMinDeg', 'ssAngleCrbDegCpK', 'ssFdRefCrbHzCpU'});
end

function checkpointSummaryTable = localBuildCheckpointSummaryTable(runStateList, checkpointDir)
%LOCALBUILDCHECKPOINTSUMMARYTABLE Build compact checkpoint summary rows.
if isempty(runStateList)
  checkpointSummaryTable = table(string(checkpointDir), 0, 0, false, "", ...
    'VariableNames', {'checkpointParentDir', 'numTask', 'numDone', 'isComplete', 'runDir'});
  checkpointSummaryTable(1, :) = [];
  return;
end
numRun = numel(runStateList);
checkpointParentDir = strings(numRun, 1);
numTask = zeros(numRun, 1);
numDone = zeros(numRun, 1);
isComplete = false(numRun, 1);
runDir = strings(numRun, 1);
for iRun = 1:numRun
  state = runStateList(iRun);
  checkpointParentDir(iRun) = string(checkpointDir);
  numTask(iRun) = localGetFieldOrDefault(state, 'numTask', 0);
  numDone(iRun) = localGetFieldOrDefault(state, 'numDone', 0);
  isComplete(iRun) = localGetFieldOrDefault(state, 'isComplete', false);
  runDir(iRun) = string(localGetFieldOrDefault(state, 'runDir', ""));
end
checkpointSummaryTable = table(checkpointParentDir, numTask, numDone, isComplete, runDir);
end

function scanData = localValidateScanDataForSummary(scanData)
%LOCALVALIDATESCANDATAFORSUMMARY Validate scanData for summary reruns.
requiredField = {'contextSummaryTable', 'referenceCandidateTable', 'pairRankTable', ...
  'greedyStepTable', 'satSetSummaryTable', 'checkpointSummaryTable'};
for iField = 1:numel(requiredField)
  fieldName = requiredField{iField};
  if ~isfield(scanData, fieldName) || ~istable(scanData.(fieldName))
    error('scanMfPairGeometryCrbAudit:MissingScanTable', ...
      'scanData.%s is missing or is not a table.', fieldName);
  end
end
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  scanData.plotData = localBuildPlotData(scanData.pairRankTable, scanData.satSetSummaryTable);
end
end

function plotData = localBuildPlotData(pairRankTable, satSetSummaryTable)
%LOCALBUILDPLOTDATA Store lightweight plotting tables.
plotData = struct();
plotData.pairRankTable = pairRankTable;
plotData.satSetSummaryTable = satSetSummaryTable;
end

function scanData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw quick-look plots from stored tables.
if ~isfield(scanData, 'plotData') || ~isstruct(scanData.plotData)
  return;
end
pairRankTable = scanData.plotData.pairRankTable;
if istable(pairRankTable) && ~isempty(pairRankTable)
  topCount = min(12, height(pairRankTable));
  topTable = pairRankTable(1:topCount, :);
  figure('Name', 'Pair CRB health ranking');
  bar(topTable.angleCrbGainCpK);
  grid on;
  xlabel('CRB health rank');
  ylabel('MS / SS angle CRB gain (CP-K)');
  title('Top second-satellite candidates by CRB health');
  set(gca, 'XTick', 1:topCount, 'XTickLabel', compose('%d', topTable.candSatIdxGlobal));
end
satSetSummaryTable = scanData.plotData.satSetSummaryTable;
if istable(satSetSummaryTable) && height(satSetSummaryTable) > 1
  figure('Name', 'Greedy selected-set CRB trend');
  yyaxis left;
  plot(satSetSummaryTable.numSat, satSetSummaryTable.angleCrbDegCpK, '-o');
  ylabel('CP-K angle CRB trace std (deg)');
  yyaxis right;
  plot(satSetSummaryTable.numSat, satSetSummaryTable.fdRefDoaCanonicalCorrCpU, '-s');
  ylabel('CP-U fdRef-DoA canonical corr');
  grid on;
  xlabel('Number of satellites');
  title('Nested greedy satellite-set CRB health');
end
end

function metricLineList = localBuildTelegramMetricLines(scanData)
%LOCALBUILDTELEGRAMMETRICLINES Build concise Telegram metrics.
numPair = height(scanData.pairRankTable);
selectedText = "";
if ~isempty(scanData.satSetSummaryTable)
  selectedText = scanData.satSetSummaryTable.selectedSatIdxGlobalStr(end);
end
metricLineList = [ ...
  "• ranked pairs: <code>" + string(numPair) + "</code>"; ...
  "• selected set: <code>" + selectedText + "</code>"; ...
  "• selected ref: <code>" + string(scanData.contextSummaryTable.selectedRefSatIdxGlobal(1)) + "</code>"; ...
  "• ref mode: <code>" + string(scanData.config.refSelectionMode) + "</code>" ...
];
end

function elevMat = localExtractElevationMatrix(sceneSeq)
%LOCALEXTRACTELEVATIONMATRIX Extract local satellite elevation over frames.
elevMat = NaN(sceneSeq.numSat, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  if isfield(sceneSeq.sceneCell{iFrame}, 'accessInfo') && ...
      isfield(sceneSeq.sceneCell{iFrame}.accessInfo, 'usrElevationDeg')
    elevMat(:, iFrame) = sceneSeq.sceneCell{iFrame}.accessInfo.usrElevationDeg(:, 1);
  end
end
end

function offAxisMat = localExtractOffAxisMatrix(sceneSeq)
%LOCALEXTRACTOFFAXISMATRIX Extract satellite off-axis angle over frames.
offAxisMat = NaN(sceneSeq.numSat, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  if isfield(sceneSeq.sceneCell{iFrame}, 'accessInfo') && ...
      isfield(sceneSeq.sceneCell{iFrame}.accessInfo, 'satOffAxisDeg')
    offAxisMat(:, iFrame) = sceneSeq.sceneCell{iFrame}.accessInfo.satOffAxisDeg(:, 1);
  end
end
end

function [maxDriftDeg, rmsDriftDeg] = localCalcLocalDirDrift(sceneSeq)
%LOCALCALCLOCALDIRDRIFT Build maximum local direction drift over all sats.
localDoaSeries = sceneSeq.localDoa;
if ndims(localDoaSeries) == 2
  localDoaSeries = reshape(localDoaSeries, 2, sceneSeq.numSat, 1);
end
driftAll = NaN(sceneSeq.numSat, sceneSeq.numFrame);
for iSat = 1:sceneSeq.numSat
  refLocalDoa = reshape(localDoaSeries(:, iSat, sceneSeq.refFrameIdx), 2, 1);
  refDir = localAngleToUnitVec(refLocalDoa);
  for iFrame = 1:sceneSeq.numFrame
    currentLocalDoa = reshape(localDoaSeries(:, iSat, iFrame), 2, 1);
    currentDir = localAngleToUnitVec(currentLocalDoa);
    driftAll(iSat, iFrame) = localAngleBetweenUnitVecDeg(refDir, currentDir);
  end
end
maxDriftDeg = max(driftAll(:), [], 'omitnan');
rmsDriftDeg = sqrt(mean(driftAll(:) .^ 2, 'omitnan'));
end

function dirVec = localAngleToUnitVec(localDoa)
%LOCALANGLETOTUNITVEC Convert [az; el] in radians to a 3x1 unit vector.
azRad = localDoa(1);
elRad = localDoa(2);
dirVec = [cos(elRad) * cos(azRad); cos(elRad) * sin(azRad); sin(elRad)];
end

function sepDeg = localMinPairAngularSeparationDeg(sceneSeq)
%LOCALMINPAIRANGULARSEPARATIONDEG Compute minimum LOS separation at ref frame.
sepDeg = NaN;
if sceneSeq.numSat < 2
  return;
end
usrPos = sceneSeq.usrPosEci(:, 1, sceneSeq.refFrameIdx);
satPos = sceneSeq.satPosEci(:, :, sceneSeq.refFrameIdx);
losMat = satPos - usrPos;
losMat = losMat ./ vecnorm(losMat, 2, 1);
sepList = [];
for iSat = 1:(sceneSeq.numSat - 1)
  for jSat = (iSat + 1):sceneSeq.numSat
    sepList(end + 1, 1) = localAngleBetweenUnitVecDeg(losMat(:, iSat), losMat(:, jSat)); %#ok<AGROW>
  end
end
sepDeg = min(sepList, [], 'omitnan');
end

function sepDeg = localRefCandidateSeparationDeg(sceneSeq, refLocalIdx, candLocalIdx)
%LOCALREFCANDIDATESEPARATIONDEG Compute ref-to-candidate LOS separation.
sepDeg = NaN;
if isnan(refLocalIdx) || candLocalIdx < 1 || candLocalIdx > sceneSeq.numSat
  return;
end
usrPos = sceneSeq.usrPosEci(:, 1, sceneSeq.refFrameIdx);
satPos = sceneSeq.satPosEci(:, :, sceneSeq.refFrameIdx);
refLos = satPos(:, refLocalIdx) - usrPos;
candLos = satPos(:, candLocalIdx) - usrPos;
sepDeg = localAngleBetweenUnitVecDeg(refLos ./ norm(refLos), candLos ./ norm(candLos));
end

function baselineKm = localBaselineProxyKm(sceneSeq)
%LOCALBASELINEPROXYKM Compute max inter-satellite distance at ref frame.
baselineKm = NaN;
if sceneSeq.numSat < 2
  baselineKm = 0;
  return;
end
satPos = sceneSeq.satPosEci(:, :, sceneSeq.refFrameIdx);
distList = [];
for iSat = 1:(sceneSeq.numSat - 1)
  for jSat = (iSat + 1):sceneSeq.numSat
    distList(end + 1, 1) = norm(satPos(:, iSat) - satPos(:, jSat)) / 1e3; %#ok<AGROW>
  end
end
baselineKm = max(distList, [], 'omitnan');
end

function losHat = localLosUnitVectorAtRef(satStateRef, satIdxGlobal, usrLla)
%LOCALLOSUNITVECTORATREF Return user-to-satellite unit vector at reference epoch.
losHat = [NaN; NaN; NaN];
if ~isstruct(satStateRef) || satIdxGlobal < 1 || satIdxGlobal > size(satStateRef.satPosEci, 2)
  return;
end
usrPos = satStateRef.usrPosEci(:, 1);
satPos = satStateRef.satPosEci(:, satIdxGlobal);
losVec = satPos - usrPos;
if norm(losVec) > 0
  losHat = losVec ./ norm(losVec);
end
end

function angleDeg = localAngleBetweenUnitVecDeg(leftVec, rightVec)
%LOCALANGLEBETWEENUNITVECDEG Compute robust vector angle in degrees.
if any(~isfinite(leftVec)) || any(~isfinite(rightVec))
  angleDeg = NaN;
  return;
end
leftVec = leftVec(:) ./ norm(leftVec(:));
rightVec = rightVec(:) ./ norm(rightVec(:));
cosVal = max(-1, min(1, real(leftVec' * rightVec)));
angleDeg = acosd(cosVal);
end

function [phaseRms, phaseMax] = localCalcDopplerFitPhaseResidual(truthDyn)
%LOCALCALCDOPPLERFITPHASERESIDUAL Approximate CP-projected phase residual from Doppler line fit.
phaseRms = NaN;
phaseMax = NaN;
if ~isfield(truthDyn, 'fdSatSeries') || ~isfield(truthDyn, 'fdSatFitHz') || ...
    ~isfield(truthDyn, 'fdRateSatTrueHzPerSec')
  return;
end
timeSec = reshape(truthDyn.timeOffsetSec, 1, []);
phaseResidAll = [];
for iSat = 1:size(truthDyn.fdSatSeries, 1)
  fdFit = truthDyn.fdSatFitHz(iSat) + truthDyn.fdRateSatTrueHzPerSec(iSat) * timeSec;
  fdResid = truthDyn.fdSatSeries(iSat, :) - fdFit;
  phaseResid = 2 * pi * cumtrapz(timeSec, fdResid);
  phaseResid = phaseResid - mean(phaseResid, 'omitnan');
  phaseResidAll = [phaseResidAll; phaseResid(:)]; %#ok<AGROW>
end
phaseRms = sqrt(mean(phaseResidAll .^ 2, 'omitnan'));
phaseMax = max(abs(phaseResidAll), [], 'omitnan');
end

function idx = localFindLocalSatIdx(satIdxGlobal, refSatIdxGlobal)
%LOCALFINDLOCALSATIDX Resolve reference local satellite index.
idx = find(double(satIdxGlobal(:)) == double(refSatIdxGlobal), 1, 'first');
if isempty(idx)
  idx = NaN;
end
end

function candSat = localGetTaskCandidate(task)
%LOCALGETTASKCANDIDATE Return the last satellite in a task set.
candSat = NaN;
if isstruct(task) && isfield(task, 'trialSatIdxGlobal') && ~isempty(task.trialSatIdxGlobal)
  candSat = double(task.trialSatIdxGlobal(end));
end
end

function satNameText = localJoinSatNames(satAccessRef, satIdxGlobal)
%LOCALJOINSATNAMES Join satellite names for one selected set.
nameList = strings(numel(satIdxGlobal), 1);
for iSat = 1:numel(satIdxGlobal)
  nameList(iSat) = localExtractSatName(satAccessRef, satIdxGlobal(iSat));
  if strlength(nameList(iSat)) == 0
    nameList(iSat) = "sat" + string(satIdxGlobal(iSat));
  end
end
satNameText = strjoin(cellstr(nameList(:).'), ',');
end

function satName = localExtractSatName(satAccessRef, satIdxGlobal)
%LOCALEXTRACTSATNAME Read one TLE satellite name when available.
satName = "";
if ~isfield(satAccessRef, 'satName') || isempty(satAccessRef.satName)
  return;
end
satNameList = satAccessRef.satName;
if iscell(satNameList)
  if satIdxGlobal >= 1 && satIdxGlobal <= numel(satNameList)
    satName = string(satNameList{satIdxGlobal});
  end
  return;
end
if isstring(satNameList) || ischar(satNameList)
  satNameList = string(satNameList);
  if satIdxGlobal >= 1 && satIdxGlobal <= numel(satNameList)
    satName = satNameList(satIdxGlobal);
  end
end
end

function value = localCrbAngleTraceStdDeg(crb)
%LOCALCRBANGLETRACESTDDEG Return sqrt(trace(2D CRB)) for lat/lon CRB in degrees.
value = NaN;
if isnumeric(crb) && size(crb, 1) >= 2 && size(crb, 2) >= 2
  traceVal = real(trace(crb(1:2, 1:2)));
  if isfinite(traceVal) && traceVal >= 0
    value = sqrt(traceVal);
  end
end
end

function value = localCrbStd(crb, idx)
%LOCALCRBSTD Return sqrt of one CRB diagonal entry.
value = NaN;
if isnumeric(crb) && size(crb, 1) >= idx && size(crb, 2) >= idx
  diagVal = real(crb(idx, idx));
  if isfinite(diagVal) && diagVal >= 0
    value = sqrt(diagVal);
  end
end
end

function value = localScalarInfoStd(infoValue)
%LOCALSCALARINFOSTD Return scalar standard deviation implied by information.
value = NaN;
if isfinite(infoValue) && infoValue > 0
  value = sqrt(1 / infoValue);
end
end

function value = localSafeMatrixElem(mat, rowIdx, colIdx)
%LOCALSAFEMATRIXELEM Read one matrix element safely.
value = NaN;
if isnumeric(mat) && size(mat, 1) >= rowIdx && size(mat, 2) >= colIdx
  value = mat(rowIdx, colIdx);
end
end

function value = localSafeCond(mat)
%LOCALSAFECOND Return condition number for finite square matrices.
value = NaN;
if isnumeric(mat) && ~isempty(mat) && size(mat, 1) == size(mat, 2) && all(isfinite(mat(:)))
  value = cond(mat);
end
end

function invMat = localSafeInv(mat)
%LOCALSAFEINV Invert or pseudo-invert one numeric matrix.
invMat = [];
if isnumeric(mat) && ~isempty(mat) && size(mat, 1) == size(mat, 2) && all(isfinite(mat(:)))
  invMat = pinv((mat + mat') / 2);
end
end

function ratio = localSafeRatio(numer, denom)
%LOCALSAFERATIO Safe element-wise ratio with scalar expansion.
try
  ratio = numer ./ denom;
catch
  ratio = NaN(size(denom));
  return;
end
mask = isfinite(numer) & isfinite(denom) & abs(denom) > 0;
try
  ratio(~mask) = NaN;
catch
  ratio = NaN(size(denom));
end
end

function z = localSafeZScore(x)
%LOCALSAFEZSCORE Robust z-score with finite fallback.
x = double(x);
finiteMask = isfinite(x);
if ~any(finiteMask(:))
  z = NaN(size(x));
  return;
end
finiteValue = x(finiteMask);
mu = mean(finiteValue);
sigma = std(finiteValue);
if ~isfinite(sigma) || sigma <= 0
  z = zeros(size(x));
  z(~finiteMask) = NaN;
else
  z = (x - mu) ./ sigma;
  z(~finiteMask) = NaN;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct or table field.
value = defaultValue;
if istable(dataStruct) && any(strcmp(dataStruct.Properties.VariableNames, fieldName))
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
elseif isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function rowStruct = localTableRowToStruct(rowTable)
%LOCALTABLEROWTOSTRUCT Convert one table row to a scalar struct.
rowStruct = table2struct(rowTable(1, :));
end

function preview = localPreviewRows(tableIn, maxRows)
%LOCALPREVIEWROWS Return the first rows of one table.
preview = tableIn;
if istable(tableIn) && height(tableIn) > maxRows
  preview = tableIn(1:maxRows, :);
end
end

function canUse = localCanUseParfor(requested, minTaskForParfor, numTask)
%LOCALCANUSEPARFOR Decide whether this task batch should use parallel execution.
canUse = logical(requested) && numTask >= minTaskForParfor && ...
  exist('gcp', 'file') == 2 && license('test', 'Distrib_Computing_Toolbox');
end

function repoRoot = localFindRepoRoot()
%LOCALFINDREPOROOT Find repository root by walking up from this script.
currentDir = fileparts(mfilename('fullpath'));
repoRoot = currentDir;
while true
  if exist(fullfile(repoRoot, 'AGENTS.md'), 'file') && exist(fullfile(repoRoot, 'README.md'), 'file')
    return;
  end
  parentDir = fileparts(repoRoot);
  if strcmp(parentDir, repoRoot)
    error('scanMfPairGeometryCrbAudit:RepoRootNotFound', ...
      'Unable to locate repository root containing AGENTS.md.');
  end
  repoRoot = parentDir;
end
end

function localPrintCheckpointFailureHint(checkpointDir)
%LOCALPRINTCHECKPOINTFAILUREHINT Print checkpoint recovery hint on failure.
if strlength(string(checkpointDir)) > 0
  fprintf('[checkpoint] Partial task results may remain under: %s\n', char(checkpointDir));
end
end
