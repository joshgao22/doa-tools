%SCANMFKNOWNUNKNOWNINFORMATIONLOSS Scan CRB information loss from unknown fdRate.
% Dev scan only. It compares known-rate and unknown-rate EFIM / CRB metrics
% across SNR, frame count, and frame interval for the paper nuisance-rate line.

clear; close all; clc;

%% Scan configuration
scanName = "scanMfKnownUnknownInformationLoss";
saveSnapshot = true;
contextSeed = 253;
snrDbList = [-5; 0; 5; 10];
frameCountList = 8:2:20;
frameIntvlSecList = [1/1000; 1/750; 1/500];
primaryPlotSnrDb = 10;
primaryFrameIntvlSec = 1/750;

snrDbList = reshape(double(snrDbList), [], 1);
frameCountList = reshape(double(frameCountList), [], 1);
frameIntvlSecList = reshape(double(frameIntvlSecList), [], 1);

scanConfig = struct();
scanConfig.saveSnapshot = saveSnapshot;
scanConfig.contextSeed = contextSeed;
scanConfig.snrDbList = snrDbList;
scanConfig.frameCountList = frameCountList;
scanConfig.frameIntvlSecList = frameIntvlSecList;
scanConfig.primaryPlotSnrDb = primaryPlotSnrDb;
scanConfig.primaryFrameIntvlSec = primaryFrameIntvlSec;

runKey = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
localPrintScanHeader(scanName, scanConfig);

try
  %% Run scan batch
  crbSummaryCell = cell(0, 1);
  lossRowList = repmat(localEmptyLossRow(), 0, 1);
  bundleCell = cell(0, 1);
  fimDiagRowList = repmat(localEmptyFimDiagRow(), 0, 1);

  for iP = 1:numel(frameCountList)
    numFrame = frameCountList(iP);
    offsetIdx = localCenteredOffsets(numFrame);
    for iTf = 1:numel(frameIntvlSecList)
      frameIntvlSec = frameIntvlSecList(iTf);
      contextOpt = struct( ...
        'frameIntvlSec', frameIntvlSec, ...
        'periodicOffsetIdx', offsetIdx, ...
        'masterOffsetIdx', offsetIdx);
      context = buildDynamicDualSatEciContext(contextOpt);
      periodicFixture = localBuildPeriodicFixture(context, primaryPlotSnrDb, contextSeed);
      for iSnr = 1:numel(snrDbList)
        snrDb = snrDbList(iSnr);
        noiseVar = 1 / (10^(snrDb / 10));
        crbBundle = localBuildCrbBundleQuiet(periodicFixture, context.pilotWave, ...
          context.carrierFreq, context.waveInfo.sampleRate, noiseVar);
        crbTable = buildDynamicCrbSummaryTable(crbBundle.truth, ...
          crbBundle.crbSfRef, crbBundle.auxCrbSfRef, ...
          crbBundle.crbSfMs, crbBundle.auxCrbSfMs, ...
          crbBundle.crbMfRefKnown, crbBundle.auxCrbMfRefKnown, ...
          crbBundle.crbMfMsKnown, crbBundle.auxCrbMfMsKnown, ...
          crbBundle.crbMfRefUnknown, crbBundle.auxCrbMfRefUnknown, ...
          crbBundle.crbMfMsUnknown, crbBundle.auxCrbMfMsUnknown, snrDb);
        crbTable.numFrame = repmat(numFrame, height(crbTable), 1);
        crbTable.frameIntvlSec = repmat(frameIntvlSec, height(crbTable), 1);
        crbTable.windowSec = repmat((numFrame - 1) * frameIntvlSec, height(crbTable), 1);
        crbSummaryCell{end + 1, 1} = crbTable; %#ok<AGROW>
        fimDiagRows = localBuildFimDiagRows(crbBundle, snrDb, numFrame, frameIntvlSec);
        for iDiag = 1:numel(fimDiagRows)
          fimDiagRowList(end + 1, 1) = fimDiagRows(iDiag); %#ok<AGROW>
        end
        lossRowList(end + 1) = localBuildLossRow("single", crbTable, crbBundle.auxCrbMfRefKnown, ...
          crbBundle.auxCrbMfRefUnknown, snrDb, numFrame, frameIntvlSec, offsetIdx); %#ok<AGROW>
        lossRowList(end + 1) = localBuildLossRow("multi", crbTable, crbBundle.auxCrbMfMsKnown, ...
          crbBundle.auxCrbMfMsUnknown, snrDb, numFrame, frameIntvlSec, offsetIdx); %#ok<AGROW>
        bundleCell{end + 1, 1} = localStripCrbBundle(crbBundle, snrDb, numFrame, frameIntvlSec); %#ok<AGROW>
      end
    end
  end

  crbSummaryTable = vertcat(crbSummaryCell{:});
  lossTable = struct2table(lossRowList(:));
  lossTable = localAppendLossPercentColumns(lossTable);
  fimDiagTable = struct2table(fimDiagRowList(:));
  fimDiagSummaryTable = localBuildFimDiagSummaryTable(fimDiagTable);
  primarySliceTable = localBuildPrimarySliceTable(lossTable, primaryPlotSnrDb);
  mainSliceTable = localBuildMainSliceTable(primarySliceTable, primaryFrameIntvlSec);
  mainTimeOriginClass = localGetMainTimeOriginClass(mainSliceTable);
  timeOriginSummaryTable = localBuildTimeOriginSummaryTable(primarySliceTable);
  primaryDisplayTable = localBuildPrimaryDisplayTable(mainSliceTable);
  tfSensitivityDisplayTable = localBuildTfSensitivityDisplayTable(primarySliceTable);
  snrSensitivityDisplayTable = localBuildSnrSensitivityDisplayTable(lossTable, primaryFrameIntvlSec);

  %% Data storage
  scanData = struct();
  scanData.scanName = string(scanName);
  scanData.runKey = string(runKey);
  scanData.utcRun = datetime('now', 'TimeZone', 'local');
  scanData.config = scanConfig;
  scanData.crbSummaryTable = crbSummaryTable;
  scanData.lossTable = lossTable;
  scanData.fimDiagTable = fimDiagTable;
  scanData.fimDiagSummaryTable = fimDiagSummaryTable;
  scanData.primarySliceTable = primarySliceTable;
  scanData.mainSliceTable = mainSliceTable;
  scanData.mainTimeOriginClass = mainTimeOriginClass;
  scanData.timeOriginSummaryTable = timeOriginSummaryTable;
  scanData.primaryDisplayTable = primaryDisplayTable;
  scanData.tfSensitivityDisplayTable = tfSensitivityDisplayTable;
  scanData.snrSensitivityDisplayTable = snrSensitivityDisplayTable;
  scanData.crbBundleSummaryCell = bundleCell;
  scanData.plotData = localBuildPlotData(primarySliceTable);

  if saveSnapshot
    saveOpt = struct('includeVars', {{'scanData'}}, ...
      'extraMeta', struct('scanName', char(scanName)), 'verbose', true);
    scanData.snapshotFile = saveExpSnapshot(char(scanName), saveOpt);
  else
    scanData.snapshotFile = "";
  end
catch ME
  fprintf('Scan failed while building known/unknown information-loss data.\n');
  rethrow(ME);
end

%% Summary output and plotting
if ~exist('scanData', 'var') || ~isstruct(scanData)
  error('scanMfKnownUnknownInformationLoss:MissingScanData', ...
    'Scan data is missing. Run the scan batch sections or load a snapshot containing scanData.');
end

fprintf('\n========== Main U/K information-loss slice ==========%s', newline);
fprintf('K = known fdRate, U = unknown fdRate nuisance. Loss percentages are 100*(U/K - 1).\n');
fprintf('Main table uses one consistent time-origin class: %s.\n', char(scanData.mainTimeOriginClass));
disp(scanData.primaryDisplayTable);
if height(scanData.timeOriginSummaryTable) > 1
  fprintf('\n========== Time-origin classes in primary SNR slice ==========\n');
  fprintf('Central-reference rows are kept in scanData but excluded from the main P-trend table.\n');
  disp(scanData.timeOriginSummaryTable);
end
fprintf('\n========== FIM condition summary ==========%s', newline);
fprintf('Expected full-FIM ill-conditioning warnings are suppressed in this scan; see the compact summary below.\n');
disp(scanData.fimDiagSummaryTable);
scanData.plotData = localPlotScan(scanData);

%% Local helpers

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.
numFrame = round(numFrame);
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end


function periodicFixture = localBuildPeriodicFixture(context, snrDb, taskSeed)
%LOCALBUILDPERIODICFIXTURE Build only the periodic fixture needed by CRB.
rng(taskSeed);
pwrNoise = 1 / (10^(snrDb / 10));
pathGainCellMaster = repmat({ones(context.sceneSeqMaster.numSat, context.sceneSeqMaster.numUser)}, ...
  1, context.sceneSeqMaster.numFrame);
rxSigCellMaster = genMultiFrameSnapshots(context.sceneSeqMaster, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, pwrNoise, pathGainCellMaster, context.simOpt);
periodicFixture = buildDynamicFrameSubsetFixture( ...
  context.sceneSeqMaster, context.linkParamCellMaster, rxSigCellMaster, ...
  context.masterOffsetIdx, context.periodicOffsetIdx, context.gridSize, context.searchRange, ...
  context.E, context.wavelen, context.waveInfo.sampleRate, ...
  context.fdRangeDefault, context.fdRateRangeDefault);
end


function crbBundle = localBuildCrbBundleQuiet(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%LOCALBUILDCRBBUNDLEQUIET Build CRB while suppressing expected full-FIM warnings.
warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar);
end


function rowList = localBuildFimDiagRows(crbBundle, snrDb, numFrame, frameIntvlSec)
%LOCALBUILDFIMDIAGROWS Build compact FIM conditioning diagnostics per MF CRB case.
caseList = { ...
  "single", "known", crbBundle.auxCrbMfRefKnown; ...
  "single", "unknown", crbBundle.auxCrbMfRefUnknown; ...
  "multi", "known", crbBundle.auxCrbMfMsKnown; ...
  "multi", "unknown", crbBundle.auxCrbMfMsUnknown};
rowList = repmat(localEmptyFimDiagRow(), size(caseList, 1), 1);
for iCase = 1:size(caseList, 1)
  row = localEmptyFimDiagRow();
  row.satMode = string(caseList{iCase, 1});
  row.fdRateMode = string(caseList{iCase, 2});
  row.snrDb = snrDb;
  row.numFrame = numFrame;
  row.frameIntvlSec = frameIntvlSec;
  row.windowSec = (numFrame - 1) * frameIntvlSec;
  aux = caseList{iCase, 3};
  fimFull = localGetFieldOrDefault(aux, 'fimFull', []);
  fimInterest = localGetFieldOrDefault(aux, 'fimInterest', []);
  row.rcondFull = localSymRcond(fimFull);
  row.rcondInterest = localSymRcond(fimInterest);
  row.minEigFull = localSymMinEig(fimFull);
  row.minEigInterest = localSymMinEig(fimInterest);
  row.fullFimNearSingular = isfinite(row.rcondFull) && row.rcondFull < 1e-12;
  row.interestFimNearSingular = isfinite(row.rcondInterest) && row.rcondInterest < 1e-12;
  rowList(iCase) = row;
end
end


function row = localEmptyFimDiagRow()
%LOCALEMPTYFIMDIAGROW Return a typed empty FIM diagnostic row.
row = struct('satMode', "", 'fdRateMode', "", 'snrDb', NaN, 'numFrame', NaN, ...
  'frameIntvlSec', NaN, 'windowSec', NaN, 'rcondFull', NaN, 'rcondInterest', NaN, ...
  'minEigFull', NaN, 'minEigInterest', NaN, 'fullFimNearSingular', false, ...
  'interestFimNearSingular', false);
end


function summaryTable = localBuildFimDiagSummaryTable(fimDiagTable)
%LOCALBUILDFIMDIAGSUMMARYTABLE Summarize full and interest FIM conditioning.
groupTable = unique(fimDiagTable(:, {'satMode', 'fdRateMode'}), 'rows', 'stable');
rowList = repmat(localEmptyFimDiagSummaryRow(), height(groupTable), 1);
for iGroup = 1:height(groupTable)
  mask = fimDiagTable.satMode == groupTable.satMode(iGroup) & ...
    fimDiagTable.fdRateMode == groupTable.fdRateMode(iGroup);
  row = localEmptyFimDiagSummaryRow();
  row.satMode = groupTable.satMode(iGroup);
  row.fdRateMode = groupTable.fdRateMode(iGroup);
  row.numCase = sum(mask);
  row.fullFimNearSingularCount = sum(fimDiagTable.fullFimNearSingular(mask));
  row.interestFimNearSingularCount = sum(fimDiagTable.interestFimNearSingular(mask));
  row.minRcondFull = localFiniteMin(fimDiagTable.rcondFull(mask));
  row.medianRcondFull = localFiniteMedian(fimDiagTable.rcondFull(mask));
  row.minRcondInterest = localFiniteMin(fimDiagTable.rcondInterest(mask));
  row.medianRcondInterest = localFiniteMedian(fimDiagTable.rcondInterest(mask));
  rowList(iGroup) = row;
end
summaryTable = struct2table(rowList(:));
end


function row = localEmptyFimDiagSummaryRow()
%LOCALEMPTYFIMDIAGSUMMARYROW Return a typed empty FIM summary row.
row = struct('satMode', "", 'fdRateMode', "", 'numCase', NaN, ...
  'fullFimNearSingularCount', NaN, 'interestFimNearSingularCount', NaN, ...
  'minRcondFull', NaN, 'medianRcondFull', NaN, ...
  'minRcondInterest', NaN, 'medianRcondInterest', NaN);
end


function row = localBuildLossRow(satMode, crbTable, auxKnown, auxUnknown, snrDb, numFrame, frameIntvlSec, offsetIdx)
%LOCALBUILDLOSSROW Summarize known-rate versus unknown-rate EFIM loss.
row = localEmptyLossRow();
row.satMode = string(satMode);
row.snrDb = snrDb;
row.numFrame = numFrame;
row.frameIntvlSec = frameIntvlSec;
row.windowSec = (numFrame - 1) * frameIntvlSec;
row.offsetMinIdx = min(offsetIdx);
row.offsetMaxIdx = max(offsetIdx);
row.offsetMeanIdx = mean(offsetIdx);
row.timeOriginClass = localClassifyTimeOrigin(offsetIdx);
if satMode == "single"
  knownName = "SS-MF-CP-K";
  unknownName = "SS-MF-CP-U";
else
  knownName = "MS-MF-CP-K";
  unknownName = "MS-MF-CP-U";
end
knownRow = crbTable(crbTable.displayName == knownName, :);
unknownRow = crbTable(crbTable.displayName == unknownName, :);
if ~isempty(knownRow) && ~isempty(unknownRow)
  row.angleCrbKnownStdDeg = knownRow.angleCrbStdDeg(1);
  row.angleCrbUnknownStdDeg = unknownRow.angleCrbStdDeg(1);
  row.fdCrbKnownStdHz = knownRow.fdRefCrbStdHz(1);
  row.fdCrbUnknownStdHz = unknownRow.fdRefCrbStdHz(1);
  row.angleStdLossRatio = localSafeRatio(row.angleCrbUnknownStdDeg, row.angleCrbKnownStdDeg);
  row.fdStdLossRatio = localSafeRatio(row.fdCrbUnknownStdHz, row.fdCrbKnownStdHz);
end
fimKnown = localGetFieldOrDefault(auxKnown, 'fimInterest', nan(3, 3));
fimUnknown = localGetFieldOrDefault(auxUnknown, 'fimInterest', nan(3, 3));
row.traceFimKnown = trace(fimKnown);
row.traceFimUnknown = trace(fimUnknown);
row.traceInfoRetention = localSafeRatio(row.traceFimUnknown, row.traceFimKnown);
row.traceInfoLossFraction = 1 - row.traceInfoRetention;
row.minEigKnown = localSymMinEig(fimKnown);
row.minEigUnknown = localSymMinEig(fimUnknown);
row.minEigRetention = localSafeRatio(row.minEigUnknown, row.minEigKnown);
end


function row = localEmptyLossRow()
%LOCALEMPTYLOSSROW Return a typed empty information-loss row.
row = struct('satMode', "", 'snrDb', NaN, 'numFrame', NaN, 'frameIntvlSec', NaN, ...
  'windowSec', NaN, 'offsetMinIdx', NaN, 'offsetMaxIdx', NaN, ...
  'offsetMeanIdx', NaN, 'timeOriginClass', "", ...
  'angleCrbKnownStdDeg', NaN, 'angleCrbUnknownStdDeg', NaN, ...
  'fdCrbKnownStdHz', NaN, 'fdCrbUnknownStdHz', NaN, 'angleStdLossRatio', NaN, ...
  'fdStdLossRatio', NaN, 'traceFimKnown', NaN, 'traceFimUnknown', NaN, ...
  'traceInfoRetention', NaN, 'traceInfoLossFraction', NaN, 'minEigKnown', NaN, ...
  'minEigUnknown', NaN, 'minEigRetention', NaN);
end


function bundleSummary = localStripCrbBundle(crbBundle, snrDb, numFrame, frameIntvlSec)
%LOCALSTRIPCRBBUNDLE Keep only lightweight CRB diagnostics for scanData.
bundleSummary = struct();
bundleSummary.snrDb = snrDb;
bundleSummary.numFrame = numFrame;
bundleSummary.frameIntvlSec = frameIntvlSec;
bundleSummary.truth = crbBundle.truth;
bundleSummary.auxKnownSingle = localStripAux(crbBundle.auxCrbMfRefKnown);
bundleSummary.auxUnknownSingle = localStripAux(crbBundle.auxCrbMfRefUnknown);
bundleSummary.auxKnownMulti = localStripAux(crbBundle.auxCrbMfMsKnown);
bundleSummary.auxUnknownMulti = localStripAux(crbBundle.auxCrbMfMsUnknown);
end


function auxSlim = localStripAux(aux)
%LOCALSTRIPAUX Remove heavy CRB auxiliary fields before snapshot storage.
auxSlim = struct();
auxSlim.phaseMode = localGetFieldOrDefault(aux, 'phaseMode', "");
auxSlim.fdRateMode = localGetFieldOrDefault(aux, 'fdRateMode', "");
auxSlim.numSat = localGetFieldOrDefault(aux, 'numSat', NaN);
auxSlim.numFrame = localGetFieldOrDefault(aux, 'numFrame', NaN);
auxSlim.paramName = localGetFieldOrDefault(aux, 'paramName', strings(0, 1));
auxSlim.fimInterest = localGetFieldOrDefault(aux, 'fimInterest', []);
auxSlim.rcondFull = localSymRcond(localGetFieldOrDefault(aux, 'fimFull', []));
auxSlim.rcondInterest = localSymRcond(auxSlim.fimInterest);
auxSlim.fullFimNearSingular = isfinite(auxSlim.rcondFull) && auxSlim.rcondFull < 1e-12;
auxSlim.interestFimNearSingular = isfinite(auxSlim.rcondInterest) && auxSlim.rcondInterest < 1e-12;
end


function lossTable = localAppendLossPercentColumns(lossTable)
%LOCALAPPENDLOSSPERCENTCOLUMNS Add readable U/K rollback percentage columns.
lossTable.angleLossPct = 100 * (lossTable.angleStdLossRatio - 1);
lossTable.fdRefLossPct = 100 * (lossTable.fdStdLossRatio - 1);
lossTable.traceInfoLossPct = 100 * (1 - lossTable.traceInfoRetention);
lossTable.minEigLossPct = 100 * (1 - lossTable.minEigRetention);
end


function primarySliceTable = localBuildPrimarySliceTable(lossTable, primaryPlotSnrDb)
%LOCALBUILDPRIMARYSLICETABLE Select one SNR slice for readable plots.
snrList = unique(lossTable.snrDb, 'stable');
[snrDist, snrIdx] = min(abs(snrList - primaryPlotSnrDb));
if isempty(snrIdx) || ~isfinite(snrDist)
  selectedSnrDb = lossTable.snrDb(1);
else
  selectedSnrDb = snrList(snrIdx);
end
primarySliceTable = lossTable(lossTable.snrDb == selectedSnrDb, :);
primarySliceTable = sortrows(primarySliceTable, {'satMode', 'frameIntvlSec', 'numFrame'});
end


function mainSliceTable = localBuildMainSliceTable(primarySliceTable, primaryFrameIntvlSec)
%LOCALBUILDMAINSLICETABLE Select one frame interval and one time-origin class.
mainSliceTable = localSelectNearestFrameIntvlTable(primarySliceTable, primaryFrameIntvlSec);
mainOriginClass = localPickMainTimeOriginClass(mainSliceTable);
mainSliceTable = mainSliceTable(mainSliceTable.timeOriginClass == mainOriginClass, :);
mainSliceTable = sortrows(mainSliceTable, {'satMode', 'numFrame'});
end


function originClass = localGetMainTimeOriginClass(mainSliceTable)
%LOCALGETMAINTIMEORIGINCLASS Return the single time-origin class used by main plots.
if isempty(mainSliceTable)
  originClass = "";
else
  originClass = string(mainSliceTable.timeOriginClass(1));
end
end


function originClass = localPickMainTimeOriginClass(sliceTable)
%LOCALPICKMAINTIMEORIGINCLASS Prefer the repository default right-biased reference window.
originList = unique(sliceTable.timeOriginClass, 'stable');
if any(originList == "rightBiasedRef")
  originClass = "rightBiasedRef";
elseif ~isempty(originList)
  originClass = originList(1);
else
  originClass = "";
end
end


function originClass = localClassifyTimeOrigin(offsetIdx)
%LOCALCLASSIFYTIMEORIGIN Classify whether the reference frame is central or biased.
offsetMean = mean(double(offsetIdx));
if abs(offsetMean) < 1e-12
  originClass = "centralRef";
elseif offsetMean > 0
  originClass = "rightBiasedRef";
else
  originClass = "leftBiasedRef";
end
end


function summaryTable = localBuildTimeOriginSummaryTable(primarySliceTable)
%LOCALBUILDTIMEORIGINSUMMARYTABLE Summarize frame-count classes in the primary SNR slice.
originList = unique(primarySliceTable.timeOriginClass, 'stable');
summaryTable = table();
for iOrigin = 1:numel(originList)
  mask = primarySliceTable.timeOriginClass == originList(iOrigin);
  row = table();
  row.timeOriginClass = originList(iOrigin);
  row.numRows = sum(mask);
  row.frameCountList = localFormatIntegerList(unique(primarySliceTable.numFrame(mask), 'stable'));
  row.offsetMeanIdxList = localFormatRow(unique(primarySliceTable.offsetMeanIdx(mask), 'stable'));
  summaryTable = [summaryTable; row]; %#ok<AGROW>
end
end


function selectedTable = localSelectNearestFrameIntvlTable(inputTable, targetFrameIntvlSec)
%LOCALSELECTNEARESTFRAMEINTVLTABLE Select rows at the nearest available frame interval.
frameIntvlList = unique(inputTable.frameIntvlSec, 'stable');
[~, idx] = min(abs(frameIntvlList - targetFrameIntvlSec));
selectedFrameIntvlSec = frameIntvlList(idx);
selectedTable = inputTable(abs(inputTable.frameIntvlSec - selectedFrameIntvlSec) < eps(selectedFrameIntvlSec) * 16, :);
end


function displayTable = localBuildPrimaryDisplayTable(mainSliceTable)
%LOCALBUILDPRIMARYDISPLAYTABLE Keep the command-window table focused on U/K rollback.
displayTable = table();
displayTable.satMode = mainSliceTable.satMode;
displayTable.snrDb = mainSliceTable.snrDb;
displayTable.numFrame = mainSliceTable.numFrame;
displayTable.frameIntvlMs = mainSliceTable.frameIntvlSec * 1e3;
displayTable.windowMs = mainSliceTable.windowSec * 1e3;
displayTable.offsetMeanIdx = mainSliceTable.offsetMeanIdx;
displayTable.angleLossPct = mainSliceTable.angleLossPct;
displayTable.fdRefLossPct = mainSliceTable.fdRefLossPct;
displayTable.traceInfoLossPct = mainSliceTable.traceInfoLossPct;
displayTable.minEigLossPct = mainSliceTable.minEigLossPct;
end


function displayTable = localBuildTfSensitivityDisplayTable(primarySliceTable)
%LOCALBUILDTFSENSITIVITYDISPLAYTABLE Keep compact frame-interval sensitivity data in scanData.
multiTable = primarySliceTable(primarySliceTable.satMode == "multi", :);
displayTable = table();
displayTable.numFrame = multiTable.numFrame;
displayTable.frameIntvlMs = multiTable.frameIntvlSec * 1e3;
displayTable.windowMs = multiTable.windowSec * 1e3;
displayTable.fdRefLossPct = multiTable.fdRefLossPct;
displayTable.traceInfoLossPct = multiTable.traceInfoLossPct;
end


function displayTable = localBuildSnrSensitivityDisplayTable(lossTable, primaryFrameIntvlSec)
%LOCALBUILDSNRSENSITIVITYDISPLAYTABLE Keep compact SNR sensitivity data in scanData.
snrSliceTable = localSelectNearestFrameIntvlTable(lossTable, primaryFrameIntvlSec);
displayTable = table();
displayTable.satMode = snrSliceTable.satMode;
displayTable.snrDb = snrSliceTable.snrDb;
displayTable.numFrame = snrSliceTable.numFrame;
displayTable.frameIntvlMs = snrSliceTable.frameIntvlSec * 1e3;
displayTable.fdRefLossPct = snrSliceTable.fdRefLossPct;
displayTable.traceInfoLossPct = snrSliceTable.traceInfoLossPct;
end


function plotData = localBuildPlotData(primarySliceTable)
%LOCALBUILDPLOTDATA Store lightweight table data needed to redraw figures.
plotData = struct();
plotData.primarySliceTable = primarySliceTable;
plotData.figureName = ["fdRef CRB rollback versus frame count"; ...
  "fdRef CRB rollback versus SNR"; "Frame-interval sensitivity of fdRef rollback"];
end


function plotData = localPlotScan(scanData)
%LOCALPLOTSCAN Draw focused known/unknown loss figures.
mainSliceTable = scanData.mainSliceTable;
primarySliceTable = scanData.primarySliceTable;
lossTable = scanData.lossTable;
primarySnrDb = mainSliceTable.snrDb(1);
primaryFrameIntvlMs = mainSliceTable.frameIntvlSec(1) * 1e3;

figure('Name', 'fdRef CRB rollback versus frame count');
localPlotBySatMode(mainSliceTable, 'fdRefLossPct');
grid on; xlabel('frame count P'); ylabel('fdRef CRB std rollback (%)');
title(sprintf('Unknown-rate fdRef rollback, SNR %.2f dB, Tf %.4g ms', ...
  primarySnrDb, primaryFrameIntvlMs));

figure('Name', 'fdRef CRB rollback versus SNR');
snrSliceTable = localSelectNearestFrameIntvlTable(lossTable, scanData.config.primaryFrameIntvlSec);
snrSliceTable = snrSliceTable(snrSliceTable.timeOriginClass == scanData.mainTimeOriginClass, :);
snrPlotFrameCountList = localPickRepresentativeFrameCounts(mainSliceTable.numFrame);
snrSliceTable = snrSliceTable(ismember(snrSliceTable.numFrame, snrPlotFrameCountList), :);
localPlotByFrameCount(snrSliceTable(snrSliceTable.satMode == "multi", :), ...
  'snrDb', 'fdRefLossPct');
grid on; xlabel('SNR (dB)'); ylabel('fdRef CRB std rollback (%)');
title(sprintf('SNR sensitivity for multi-sat, Tf %.4g ms', primaryFrameIntvlMs));

figure('Name', 'Frame-interval sensitivity of fdRef rollback');
tfPlotTable = primarySliceTable(primarySliceTable.satMode == "multi" & ...
  primarySliceTable.timeOriginClass == scanData.mainTimeOriginClass, :);
localPlotByFrameInterval(tfPlotTable, 'fdRefLossPct');
grid on; xlabel('frame count P'); ylabel('fdRef CRB std rollback (%)');
title(sprintf('Frame-interval sensitivity for multi-sat, SNR %.2f dB', primarySnrDb));

plotData = scanData.plotData;
plotData.mainSliceTable = mainSliceTable;
plotData.primarySliceTable = primarySliceTable;
plotData.snrSliceTable = snrSliceTable;
end


function localPlotBySatMode(plotTable, yField)
%LOCALPLOTBYSATMODE Plot a primary metric against frame count for each sat mode.
groupList = unique(plotTable.satMode, 'stable');
hold on;
legendText = strings(numel(groupList), 1);
for iGroup = 1:numel(groupList)
  mask = plotTable.satMode == groupList(iGroup);
  groupTable = sortrows(plotTable(mask, :), 'numFrame');
  plot(groupTable.numFrame, groupTable.(yField), '-o');
  legendText(iGroup) = groupList(iGroup);
end
legend(cellstr(legendText), 'Location', 'best');
end


function localPlotByFrameCount(plotTable, xField, yField)
%LOCALPLOTBYFRAMECOUNT Plot a metric versus x for each frame-count group.
groupList = unique(plotTable.numFrame, 'stable');
hold on;
legendText = strings(numel(groupList), 1);
for iGroup = 1:numel(groupList)
  mask = plotTable.numFrame == groupList(iGroup);
  groupTable = sortrows(plotTable(mask, :), xField);
  plot(groupTable.(xField), groupTable.(yField), '-o');
  legendText(iGroup) = sprintf('P=%d', groupList(iGroup));
end
legend(cellstr(legendText), 'Location', 'best');
end


function frameCountList = localPickRepresentativeFrameCounts(frameCountListIn)
%LOCALPICKREPRESENTATIVEFRAMECOUNTS Keep the SNR figure readable with low/mid/high P.
frameCountListIn = unique(frameCountListIn(:), 'stable');
if numel(frameCountListIn) <= 3
  frameCountList = frameCountListIn;
  return;
end
[~, midIdx] = min(abs(frameCountListIn - 10));
frameCountList = unique([frameCountListIn(1); frameCountListIn(midIdx); frameCountListIn(end)], 'stable');
end


function txt = localFormatIntegerList(x)
%LOCALFORMATINTEGERLIST Format integer vector values for compact diagnostic tables.
if isempty(x)
  txt = "";
else
  txt = strjoin(compose('%d', reshape(round(x), 1, [])), ', ');
end
end


function localPlotByFrameInterval(plotTable, yField)
%LOCALPLOTBYFRAMEINTERVAL Plot a metric versus frame count for each frame interval.
groupList = unique(plotTable.frameIntvlSec, 'stable');
hold on;
legendText = strings(numel(groupList), 1);
for iGroup = 1:numel(groupList)
  mask = plotTable.frameIntvlSec == groupList(iGroup);
  groupTable = sortrows(plotTable(mask, :), 'numFrame');
  plot(groupTable.numFrame, groupTable.(yField), '-o');
  legendText(iGroup) = sprintf('Tf=%.4g ms', groupList(iGroup) * 1e3);
end
legend(cellstr(legendText), 'Location', 'best');
end


function txt = localFormatRow(x)
%LOCALFORMATROW Format numeric vector for compact log output.
txt = strjoin(compose('%.6g', x(:).'), ', ');
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a struct field while preserving a fallback.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end


function ratio = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Compute a scalar ratio with NaN for invalid denominators.
if ~isscalar(numerator) || ~isscalar(denominator) || ~isfinite(denominator) || denominator == 0
  ratio = NaN;
else
  ratio = numerator / denominator;
end
end


function reciprocalCondition = localSymRcond(fim)
%LOCALSYMRCOND Compute a reciprocal condition estimate for a symmetrized matrix.
if isempty(fim) || any(~isfinite(fim(:)))
  reciprocalCondition = NaN;
  return;
end
fimSym = 0.5 * (fim + fim.');
reciprocalCondition = rcond(fimSym + eye(size(fimSym)) * eps);
end


function value = localFiniteMin(x)
%LOCALFINITEMIN Return the minimum finite value from a numeric vector.
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = min(x);
end
end


function value = localFiniteMedian(x)
%LOCALFINITEMEDIAN Return the median finite value from a numeric vector.
x = x(isfinite(x));
if isempty(x)
  value = NaN;
else
  value = median(x);
end
end


function minEig = localSymMinEig(fim)
%LOCALSYMMINEIG Compute the minimum eigenvalue of a symmetrized FIM block.
if isempty(fim) || any(~isfinite(fim(:)))
  minEig = NaN;
  return;
end
fimSym = 0.5 * (fim + fim.');
minEig = min(real(eig(fimSym)));
end


function localPrintScanHeader(scanName, scanConfig)
%LOCALPRINTSCANHEADER Print compact scan settings before heavy work starts.
fprintf('Running %s ...\n', char(scanName));
fprintf('  context seed                    : %d\n', scanConfig.contextSeed);
fprintf('  SNR list (dB)                   : %s\n', localFormatRow(scanConfig.snrDbList));
fprintf('  frame count list                : %s\n', localFormatRow(scanConfig.frameCountList));
fprintf('  frame interval list (ms)        : %s\n', localFormatRow(scanConfig.frameIntvlSecList * 1e3));
fprintf('  primary plot SNR (dB)           : %.2f\n', scanConfig.primaryPlotSnrDb);
fprintf('  primary frame interval (ms)     : %.6g\n', scanConfig.primaryFrameIntvlSec * 1e3);
fprintf('  CRB grid cases                  : %d\n', ...
  numel(scanConfig.snrDbList) * numel(scanConfig.frameCountList) * numel(scanConfig.frameIntvlSecList));
fprintf('  save snapshot                   : %d\n', logical(scanConfig.saveSnapshot));
end
