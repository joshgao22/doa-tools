% replayMfCrbTexConsistencyAudit
% Purpose: audit CRB model dispatch, parameter partition, FIM / EFIM blocks,
% and key case-pair differences against the paper-level SF-DoA and MF
% DoA-Doppler model hierarchy.
%
% Usage: edit the configuration block below, then run this script directly.
% The replay does not run MLE. It builds one fixed paper-facing multi-sat
% dynamic context, evaluates CRB helpers for the selected model-matched cases,
% and stores lightweight audit tables in replayData.

clear; close all; clc;

%% Replay configuration

replayName = "replayMfCrbTexConsistencyAudit";
saveSnapshot = true;
notifyTelegramEnable = false;

% This audit is model-matched: SF cases use the DoA-only effective CRB;
% MF cases use their corresponding zero-rate / known-rate / unknown-rate CRB.
auditMode = "model-matched";
numFrame = 10;
snrDb = 0;
baseSeed = 253;
usrLla = [55; 36.59; 0];
utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tleFileName = "statlink_20260318.tle";
selectedSatIdxGlobal = [5259 1243 348 5652 14 4437];
refSatIdxGlobal = 5259;
methodNameList = [
  "SS-SF-DoA"
  "MS-SF-DoA"
  "SS-MF-Static"
  "MS-MF-Static"
  "SS-MF-CP-K"
  "SS-MF-CP-U"
  "MS-MF-CP-K"
  "MS-MF-CP-U"
];

identicalRelTol = 1e-10;
summaryReuseRelTol = 1e-10;
saveSnapshot = logical(saveSnapshot);
notifyTelegramEnable = logical(notifyTelegramEnable);
methodNameList = reshape(string(methodNameList), [], 1);

runTic = tic;
replayData = struct();
config = struct();

try
  %% Build context and CRB audit options

  config.replayName = string(replayName);
  config.auditMode = string(auditMode);
  config.numFrame = double(numFrame);
  config.snrDb = double(snrDb);
  config.baseSeed = double(baseSeed);
  config.usrLla = reshape(double(usrLla), [], 1);
  config.utcRef = utcRef;
  config.tleFileName = string(tleFileName);
  config.selectedSatIdxGlobal = reshape(double(selectedSatIdxGlobal), 1, []);
  config.refSatIdxGlobal = double(refSatIdxGlobal);
  config.numSatSingle = 1;
  config.numSatMulti = NaN;
  config.methodNameList = methodNameList;
  config.identicalRelTol = double(identicalRelTol);
  config.summaryReuseRelTol = double(summaryReuseRelTol);
  config.saveSnapshot = saveSnapshot;
  config.notifyTelegramEnable = notifyTelegramEnable;
  config.runKey = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));

  printMfReplayHeader(char(replayName), config, '');
  localPrintReplayConfig(config);

  contextOpt = struct();
  contextOpt.usrLla = config.usrLla;
  contextOpt.utcRef = config.utcRef;
  contextOpt.tleFileName = config.tleFileName;
  contextOpt.selectedSatIdxGlobal = config.selectedSatIdxGlobal;
  contextOpt.refSatIdxGlobal = config.refSatIdxGlobal;
  contextOpt.periodicOffsetIdx = localCenteredOffsets(config.numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', inf);
  context = buildDynamicMultiSatEciContext(contextOpt);
  config.numSatMulti = context.sceneSeqMaster.numSat;
  context.subsetOffsetCell = {};
  context.subsetLabelList = strings(0, 1);
  context.numSubsetRandomTrial = 0;

  periodicFixture = localBuildPeriodicFixture(context, config.snrDb, config.baseSeed);

  %% Run CRB audit

  crbCaseList = localBuildCrbCaseList(periodicFixture, context, config);
  crbModelAuditTable = localBuildModelAuditTable(crbCaseList);
  crbFimBlockAuditTable = localBuildFimBlockAuditTable(crbCaseList);
  crbFdRefMarginalAuditTable = localBuildFdRefMarginalAuditTable(crbFimBlockAuditTable);
  crbEigenCouplingAuditTable = localBuildEigenCouplingAuditTable(crbCaseList);
  crbDerivativeNormTable = localBuildDerivativeNormTable(crbCaseList);
  crbSatContributionAuditTable = localBuildSatContributionAuditTable(crbDerivativeNormTable, context);
  crbPathAuditTable = localBuildCrbPathAuditTable(crbCaseList, periodicFixture, config);
  crbPathAuditAggregateTable = localBuildCrbPathAuditAggregateTable(crbPathAuditTable);
  crbCasePairAuditTable = localBuildCasePairAuditTable( ...
    crbCaseList, crbFimBlockAuditTable, config);
  crbSummaryReuseAuditTable = localBuildSummaryReuseAuditTable( ...
    crbCasePairAuditTable, config);
  crbContextSummaryTable = localBuildContextSummaryTable(context, periodicFixture);

  %% Data storage

  replayData = struct();
  replayData.replayName = string(replayName);
  replayData.runKey = string(config.runKey);
  replayData.utcRun = datetime('now', 'TimeZone', 'local');
  replayData.config = config;
  replayData.contextSummaryTable = crbContextSummaryTable;
  replayData.crbModelAuditTable = crbModelAuditTable;
  replayData.crbFimBlockAuditTable = crbFimBlockAuditTable;
  replayData.crbFdRefMarginalAuditTable = crbFdRefMarginalAuditTable;
  replayData.crbEigenCouplingAuditTable = crbEigenCouplingAuditTable;
  replayData.crbDerivativeNormTable = crbDerivativeNormTable;
  replayData.crbSatContributionAuditTable = crbSatContributionAuditTable;
  replayData.crbPathAuditTable = crbPathAuditTable;
  replayData.crbPathAuditAggregateTable = crbPathAuditAggregateTable;
  replayData.crbCasePairAuditTable = crbCasePairAuditTable;
  replayData.crbSummaryReuseAuditTable = crbSummaryReuseAuditTable;
  replayData.elapsedSec = toc(runTic);
  replayData = finalizeMfReplayResult(replayData, '');

  if config.saveSnapshot
    saveOpt = struct();
    saveOpt.includeVars = {'replayData'};
    saveOpt.extraMeta = struct('replayName', char(replayName), 'auditMode', char(config.auditMode));
    saveOpt.verbose = true;
    replayData.snapshotFile = saveExpSnapshot(char(replayName), saveOpt);
  else
    replayData.snapshotFile = '';
  end

  %% Summary output and plotting

  if ~exist('replayData', 'var') || ~isstruct(replayData)
    error('replayMfCrbTexConsistencyAudit:MissingReplayData', ...
      'replayData is missing. Load a snapshot containing replayData before running this section.');
  end
  replayData = localValidateReplayData(replayData);

  printMfReplaySection('CRB context summary', replayData.contextSummaryTable);
  printMfReplaySection('CRB model dispatch audit', replayData.crbModelAuditTable);
  printMfReplaySection('CRB K/U path audit aggregate', replayData.crbPathAuditAggregateTable);
  printMfReplaySection('CRB K/U path audit detail', replayData.crbPathAuditTable);
  printMfReplaySection('CRB FIM / EFIM block audit', replayData.crbFimBlockAuditTable);
  printMfReplaySection('CRB fdRef marginal audit', replayData.crbFdRefMarginalAuditTable);
  printMfReplaySection('CRB eigen / coupling audit', replayData.crbEigenCouplingAuditTable);
  printMfReplaySection('CRB case-pair audit', replayData.crbCasePairAuditTable);
  printMfReplaySection('CRB summary-reuse suspects', replayData.crbSummaryReuseAuditTable);
  printMfReplaySection('CRB per-sat contribution audit', replayData.crbSatContributionAuditTable);
  printMfReplaySection('CRB derivative norm preview', localPreviewRows(replayData.crbDerivativeNormTable, 24));

  if config.notifyTelegramEnable
    notifyMfReplayStatus(struct( ...
      'replayName', string(replayName), ...
      'statusText', "DONE", ...
      'config', config, ...
      'snapshotFile', localGetFieldOrDefault(replayData, 'snapshotFile', ''), ...
      'elapsedSec', replayData.elapsedSec, ...
      'metricLineList', localBuildTelegramMetricLines(replayData), ...
      'commentLineList', [ ...
        "CRB tex-consistency audit completed."; ...
        "Inspect unexpected-identical and summary-reuse-suspect pairs before changing CRB code."]));
  end

catch ME
  if notifyTelegramEnable
    notifyMfReplayStatus(struct( ...
      'replayName', string(replayName), ...
      'statusText', "FAILED", ...
      'config', config, ...
      'elapsedSec', toc(runTic), ...
      'errorObj', ME));
  end
  rethrow(ME);
end

%% Local helpers

function localPrintReplayConfig(config)
%LOCALPRINTREPLAYCONFIG Print the compact replay configuration.

printMfReplaySection('Replay configuration');
fprintf('auditMode: %s\n', char(config.auditMode));
fprintf('numFrame : %d\n', config.numFrame);
fprintf('snrDb    : %.1f dB\n', config.snrDb);
fprintf('baseSeed : %d\n', config.baseSeed);
fprintf('user LLA : %s\n', mat2str(config.usrLla(:).'));
fprintf('satellites: %s\n', mat2str(config.selectedSatIdxGlobal));
fprintf('ref sat  : %d\n', config.refSatIdxGlobal);
fprintf('methods  : %s\n', strjoin(cellstr(config.methodNameList.'), ', '));
end

function offsetIdx = localCenteredOffsets(numFrame)
%LOCALCENTEREDOFFSETS Build centered periodic frame offsets.

numFrame = round(double(numFrame));
if mod(numFrame, 2) == 0
  offsetIdx = (-(numFrame / 2 - 1)):(numFrame / 2);
else
  halfCount = floor(numFrame / 2);
  offsetIdx = -halfCount:halfCount;
end
offsetIdx = reshape(offsetIdx, 1, []);
end

function periodicFixture = localBuildPeriodicFixture(context, snrDb, taskSeed)
%LOCALBUILDPERIODICFIXTURE Build the periodic fixture used by the CRB audit.

rng(taskSeed, 'twister');
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
periodicFixture.pwrNoise = pwrNoise;
end

function crbCaseList = localBuildCrbCaseList(periodicFixture, context, config)
%LOCALBUILDCRBCASELIST Build the selected model-matched CRB cases.

numCase = numel(config.methodNameList);
crbCaseList = repmat(localEmptyCrbCase(), numCase, 1);
for iCase = 1:numCase
  methodName = config.methodNameList(iCase);
  crbCaseList(iCase) = localBuildSingleCrbCase(methodName, periodicFixture, context, config);
end
end

function caseInfo = localBuildSingleCrbCase(methodName, periodicFixture, context, config)
%LOCALBUILDSINGLECRBCASE Build one CRB case and record its expected model tags.

caseInfo = localEmptyCrbCase();
caseInfo.caseName = string(methodName);
caseInfo.expected = localExpectedModelSpec(methodName, config);
caseInfo.crb = [];
caseInfo.aux = struct();
caseInfo.errorMessage = "";

truth = periodicFixture.truth;
noiseVar = periodicFixture.pwrNoise;
pilotWave = context.pilotWave;
carrierFreq = context.carrierFreq;
sampleRate = context.waveInfo.sampleRate;

try
  switch string(methodName)
    case "SS-SF-DoA"
      sceneRefOnly = periodicFixture.sceneSeqRefOnly.sceneCell{periodicFixture.sceneSeqRefOnly.refFrameIdx};
      crbOpt = struct('doaType', 'latlon', 'effectiveGainMode', 'pilotModel');
      [caseInfo.crb, caseInfo.aux] = crbPilotSfDoaOnlyEffective( ...
        sceneRefOnly, pilotWave, carrierFreq, sampleRate, ...
        truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbOpt);

    case "MS-SF-DoA"
      crbOpt = struct('doaType', 'latlon', 'effectiveGainMode', 'pilotModel');
      [caseInfo.crb, caseInfo.aux] = crbPilotSfDoaOnlyEffective( ...
        periodicFixture.sceneRef, pilotWave, carrierFreq, sampleRate, ...
        truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbOpt);

    case {"SS-MF-Static", "SS-MF-CP-K", "SS-MF-CP-U"}
      sceneSeq = periodicFixture.sceneSeqRefOnly;
      [pathGain, noiseMat] = localOnesPathAndNoise(sceneSeq, noiseVar);
      crbOpt = localBuildMfCrbOpt(methodName, periodicFixture);
      [fdRef, fdRate] = localGetMfFrequencyTruth(methodName, truth);
      [caseInfo.crb, caseInfo.aux] = localBuildMfCrbQuiet( ...
        sceneSeq, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, ...
        fdRef, fdRate, pathGain, noiseMat, crbOpt);

    case {"MS-MF-Static", "MS-MF-CP-K", "MS-MF-CP-U"}
      sceneSeq = periodicFixture.sceneSeq;
      [pathGain, noiseMat] = localOnesPathAndNoise(sceneSeq, noiseVar);
      crbOpt = localBuildMfCrbOpt(methodName, periodicFixture);
      [fdRef, fdRate] = localGetMfFrequencyTruth(methodName, truth);
      [caseInfo.crb, caseInfo.aux] = localBuildMfCrbQuiet( ...
        sceneSeq, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, ...
        fdRef, fdRate, pathGain, noiseMat, crbOpt);

    otherwise
      error('replayMfCrbTexConsistencyAudit:UnsupportedMethod', ...
        'Unsupported CRB audit method: %s', methodName);
  end
catch ME
  caseInfo.crb = NaN(caseInfo.expected.numInterest);
  caseInfo.aux = struct();
  caseInfo.errorMessage = string(ME.message);
end

caseInfo.actual = localActualModelSpec(caseInfo);
end

function [pathGain, noiseMat] = localOnesPathAndNoise(sceneSeq, noiseVar)
%LOCALONESPATHANDNOISE Build unit deterministic gains and per-block noise.

pathGain = ones(sceneSeq.numSat, sceneSeq.numFrame);
noiseMat = noiseVar * ones(sceneSeq.numSat, sceneSeq.numFrame);
end

function crbOpt = localBuildMfCrbOpt(methodName, periodicFixture)
%LOCALBUILDMFCRBOPT Build the model-matched MF CRB options.

crbOpt = struct();
crbOpt.doaType = 'latlon';
crbOpt.phaseMode = 'continuous';
crbOpt.steeringMode = 'framewise';
crbOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;

switch string(methodName)
  case {"SS-MF-Static", "MS-MF-Static"}
    crbOpt.fdRateMode = 'zero';
  case {"SS-MF-CP-K", "MS-MF-CP-K"}
    crbOpt.fdRateMode = 'known';
  case {"SS-MF-CP-U", "MS-MF-CP-U"}
    crbOpt.fdRateMode = 'unknown';
  otherwise
    error('replayMfCrbTexConsistencyAudit:InvalidMfMethod', ...
      'Method %s is not an MF CRB method.', methodName);
end
end

function [fdRef, fdRate] = localGetMfFrequencyTruth(methodName, truth)
%LOCALGETMFFREQUENCYTRUTH Return model-matched reference frequency values.

switch string(methodName)
  case {"SS-MF-Static", "MS-MF-Static"}
    fdRef = truth.fdRefFit;
    fdRate = 0;
  case {"SS-MF-CP-K", "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"}
    fdRef = truth.fdRefFit;
    fdRate = truth.fdRateFit;
  otherwise
    error('replayMfCrbTexConsistencyAudit:InvalidMfTruthMethod', ...
      'Method %s is not an MF CRB method.', methodName);
end
end

function [crb, aux] = localBuildMfCrbQuiet(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRef, fdRate, pathGain, noiseMat, crbOpt)
%LOCALBUILDMFCRBQUIET Build MF CRB while suppressing expected conditioning warnings.

warnState = warning;
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
warning('off', 'crbPilotMfDoaDoppler:IllConditionedInterestFim');
[crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  doaParam, fdRef, fdRate, pathGain, noiseMat, crbOpt);
end

function modelAuditTable = localBuildModelAuditTable(crbCaseList)
%LOCALBUILDMODELAUDITTABLE Build the case dispatch and parameter partition table.

rowList = repmat(localEmptyModelAuditRow(), numel(crbCaseList), 1);
for iCase = 1:numel(crbCaseList)
  caseInfo = crbCaseList(iCase);
  expected = caseInfo.expected;
  actual = caseInfo.actual;
  row = localEmptyModelAuditRow();
  row.caseName = caseInfo.caseName;
  row.crbFunction = expected.crbFunction;
  row.expectedSignalModel = expected.signalModel;
  row.satMode = expected.satMode;
  row.frameMode = expected.frameMode;
  row.expectedTimeTemplate = expected.timeTemplate;
  row.actualTimeTemplate = actual.timeTemplate;
  row.expectedPhaseMode = expected.phaseMode;
  row.actualPhaseMode = actual.phaseMode;
  row.expectedFdRateMode = expected.fdRateMode;
  row.actualFdRateMode = actual.fdRateMode;
  row.expectedSteeringMode = expected.steeringMode;
  row.actualSteeringMode = actual.steeringMode;
  row.numSat = actual.numSat;
  row.numFrame = actual.numFrame;
  row.numInterest = actual.numInterest;
  row.numNuisance = actual.numNuisance;
  row.interestNames = actual.interestNames;
  row.nuisanceNames = actual.nuisanceNames;
  numFrameOk = (~isnan(expected.numFrame) && row.numFrame == expected.numFrame) || ...
    (isnan(expected.numFrame) && row.numFrame > 1);
  row.dispatchOk = localStringEq(row.expectedPhaseMode, row.actualPhaseMode) && ...
    localStringEq(row.expectedFdRateMode, row.actualFdRateMode) && ...
    row.numSat == expected.numSat && numFrameOk;
  row.errorMessage = caseInfo.errorMessage;
  rowList(iCase) = row;
end
modelAuditTable = struct2table(rowList);
end

function fimAuditTable = localBuildFimBlockAuditTable(crbCaseList)
%LOCALBUILDFIMBLOCKAUDITTABLE Build compact raw-FIM / EFIM diagnostics.

rowList = repmat(localEmptyFimBlockRow(), numel(crbCaseList), 1);
for iCase = 1:numel(crbCaseList)
  caseInfo = crbCaseList(iCase);
  row = localEmptyFimBlockRow();
  row.caseName = caseInfo.caseName;
  row.satMode = caseInfo.expected.satMode;
  row.frameMode = caseInfo.expected.frameMode;
  row.modelTag = caseInfo.expected.modelTag;
  row.crbDim = localMatrixDim(caseInfo.crb);
  [rawFim, efim, nuisanceFim, crossFim] = localExtractFimBlocks(caseInfo);
  row.rawFimDoaTrace = localSafeTrace(rawFim, 1:2);
  row.rawFimFdFd = localSafeMatrixElem(rawFim, 3, 3);
  row.rawFimRateRate = localSafeMatrixElem(rawFim, 4, 4);
  row.efimDoaTrace = localSafeTrace(efim, 1:2);
  row.efimFdFd = localSafeMatrixElem(efim, 3, 3);
  row.crossNorm = localSafeFroNorm(crossFim);
  row.nuisanceCond = localSafeCond(nuisanceFim);
  row.efimCond = localSafeCond(efim);
  row.crbAngleTraceDeg = localCrbAngleTraceStdDeg(caseInfo.crb);
  row.crbFdRefStdHz = localCrbStd(caseInfo.crb, 3);
  efimInv = localSafeInv(efim);
  rawInv = localSafeInv(rawFim);
  row.fdRefStdFullInvHz = localCrbStd(efimInv, 3);
  row.fdRefStdScalarEfimHz = localScalarInfoStd(row.efimFdFd);
  row.fdRefStdRawInvHz = localCrbStd(rawInv, 3);
  row.fdRefInfoFullInv = localScalarInfoFromStd(row.fdRefStdFullInvHz);
  row.fdRefScalarInfoEfim = row.efimFdFd;
  row.fdRefCouplingInflation = localSafeRatio(row.fdRefStdFullInvHz, row.fdRefStdScalarEfimHz);
  row.fdRefCouplingLossFraction = localCouplingLossFraction( ...
    row.fdRefStdFullInvHz, row.fdRefStdScalarEfimHz);
  row.crbFdRefReportedRelDiff = localRelDiff(row.crbFdRefStdHz, row.fdRefStdFullInvHz);
  row.status = localCaseStatus(caseInfo);
  rowList(iCase) = row;
end
fimAuditTable = struct2table(rowList);
end

function fdRefTable = localBuildFdRefMarginalAuditTable(fimAuditTable)
%LOCALBUILDFDREFMARGINALAUDITTABLE Compare fdRef CRB marginalization paths.

if isempty(fimAuditTable)
  fdRefTable = table();
  return;
end
fdRefTable = fimAuditTable(:, { ...
  'caseName', 'satMode', 'frameMode', 'modelTag', ...
  'efimFdFd', 'crbFdRefStdHz', 'fdRefStdFullInvHz', ...
  'fdRefStdScalarEfimHz', 'fdRefStdRawInvHz', ...
  'fdRefCouplingInflation', 'fdRefCouplingLossFraction', ...
  'crbFdRefReportedRelDiff', 'status'});
end

function eigenTable = localBuildEigenCouplingAuditTable(crbCaseList)
%LOCALBUILDEIGENCOUPLINGAUDITTABLE Explain EFIM eigen axes and DoA-fdRef coupling.

rowList = repmat(localEmptyEigenCouplingRow(), numel(crbCaseList), 1);
for iCase = 1:numel(crbCaseList)
  caseInfo = crbCaseList(iCase);
  row = localEmptyEigenCouplingRow();
  row.caseName = caseInfo.caseName;
  row.satMode = caseInfo.expected.satMode;
  row.frameMode = caseInfo.expected.frameMode;
  row.modelTag = caseInfo.expected.modelTag;
  row.status = localCaseStatus(caseInfo);
  [~, efim] = localExtractFimBlocks(caseInfo);
  if isnumeric(efim) && size(efim, 1) >= 3 && size(efim, 2) >= 3
    efim = (efim + efim') / 2;
    [eigValue, eigVector] = localSortedEigenSystem(efim(1:3, 1:3));
    if numel(eigValue) >= 3
      row.eigInfoMin = eigValue(1);
      row.eigInfoMid = eigValue(2);
      row.eigInfoMax = eigValue(3);
      row.efimCondEig = localSafeRatio(row.eigInfoMax, row.eigInfoMin);
      row.stdAlongMinInfoAxis = localScalarInfoStd(row.eigInfoMin);
      row.stdAlongMidInfoAxis = localScalarInfoStd(row.eigInfoMid);
      row.stdAlongMaxInfoAxis = localScalarInfoStd(row.eigInfoMax);
      row.fdRefWeightMinInfoAxis = abs(eigVector(3, 1));
      row.doaWeightMinInfoAxis = norm(eigVector(1:2, 1));
      row.fdRefWeightMaxInfoAxis = abs(eigVector(3, 3));
      row.doaWeightMaxInfoAxis = norm(eigVector(1:2, 3));
    end
    doaBlock = efim(1:2, 1:2);
    fdInfo = efim(3, 3);
    couplingVec = efim(1:2, 3);
    doaInv = localSafeInv(doaBlock);
    row.doaBlockCond = localSafeCond(doaBlock);
    row.corrDoa1FdRef = localNormalizedFimCoupling(efim, 1, 3);
    row.corrDoa2FdRef = localNormalizedFimCoupling(efim, 2, 3);
    row.fdRefInfoDoaKnown = fdInfo;
    row.fdRefStdDoaKnownHz = localScalarInfoStd(fdInfo);
    if ~isempty(doaInv) && isfinite(fdInfo)
      couplingInfo = real(couplingVec' * doaInv * couplingVec);
      marginalInfo = fdInfo - couplingInfo;
      row.fdRefDoaCouplingInfo = couplingInfo;
      row.fdRefInfoMarginalSchur = marginalInfo;
      row.fdRefDoaCouplingFraction = localSafeRatio(couplingInfo, fdInfo);
      row.fdRefDoaCanonicalCorr = sqrt(max(0, min(1, row.fdRefDoaCouplingFraction)));
      row.fdRefStdDoaUnknownHz = localScalarInfoStd(marginalInfo);
    end
    efimInv = localSafeInv(efim(1:3, 1:3));
    row.fdRefStdFullInvHz = localCrbStd(efimInv, 3);
    row.fdRefInfoMarginalFullInv = localScalarInfoFromStd(row.fdRefStdFullInvHz);
    row.fdRefMarginalInfoRelDiff = localRelDiff( ...
      row.fdRefInfoMarginalSchur, row.fdRefInfoMarginalFullInv);
  end
  rowList(iCase) = row;
end
eigenTable = struct2table(rowList);
end

function derivTable = localBuildDerivativeNormTable(crbCaseList)
%LOCALBUILDDERIVATIVENORMTABLE Build per-satellite / per-frame derivative diagnostics.

rowCell = cell(numel(crbCaseList), 1);
for iCase = 1:numel(crbCaseList)
  caseInfo = crbCaseList(iCase);
  if caseInfo.expected.isMf
    rowCell{iCase} = localBuildMfDerivativeRows(caseInfo);
  else
    rowCell{iCase} = localBuildSfContributionRows(caseInfo);
  end
end
if isempty(rowCell)
  derivTable = struct2table(localEmptyDerivativeRow());
  derivTable(1, :) = [];
else
  derivTable = vertcat(rowCell{:});
end
end

function pairTable = localBuildCasePairAuditTable(crbCaseList, fimAuditTable, config)
%LOCALBUILDCASEPAIRAUDITTABLE Compare key case pairs for unexpected identity.

pairSpec = localCasePairSpec();
rowList = repmat(localEmptyCasePairRow(), numel(pairSpec), 1);
for iPair = 1:numel(pairSpec)
  spec = pairSpec(iPair);
  leftIdx = find([crbCaseList.caseName] == spec.leftCase, 1);
  rightIdx = find([crbCaseList.caseName] == spec.rightCase, 1);
  row = localEmptyCasePairRow();
  row.leftCase = spec.leftCase;
  row.rightCase = spec.rightCase;
  row.expectedRelation = spec.expectedRelation;
  if isempty(leftIdx) || isempty(rightIdx)
    row.auditClass = "missing-case";
    rowList(iPair) = row;
    continue;
  end
  leftCase = crbCaseList(leftIdx);
  rightCase = crbCaseList(rightIdx);
  leftFim = fimAuditTable(leftIdx, :);
  rightFim = fimAuditTable(rightIdx, :);
  row.sameTimeTemplateFlag = localSameScalarString( ...
    leftCase.actual.timeTemplate, rightCase.actual.timeTemplate);
  row.sameInterestFlag = localSameNameList( ...
    leftCase.actual.interestNames, rightCase.actual.interestNames);
  row.sameNuisanceFlag = localSameNameList( ...
    leftCase.actual.nuisanceNames, rightCase.actual.nuisanceNames);
  row.relRawFimFdDiff = localRelDiff(leftFim.rawFimFdFd, rightFim.rawFimFdFd);
  row.relEfimFdDiff = localRelDiff(leftFim.efimFdFd, rightFim.efimFdFd);
  row.relCrbFdRefDiff = localRelDiff(leftFim.crbFdRefStdHz, rightFim.crbFdRefStdHz);
  row.relFullInvFdRefDiff = localRelDiff(leftFim.fdRefStdFullInvHz, rightFim.fdRefStdFullInvHz);
  row.relScalarEfimFdRefStdDiff = localRelDiff( ...
    leftFim.fdRefStdScalarEfimHz, rightFim.fdRefStdScalarEfimHz);
  row.relReportedFullInvLeft = localRelDiff(leftFim.crbFdRefStdHz, leftFim.fdRefStdFullInvHz);
  row.relReportedFullInvRight = localRelDiff(rightFim.crbFdRefStdHz, rightFim.fdRefStdFullInvHz);
  row.relCouplingInflationDiff = localRelDiff( ...
    leftFim.fdRefCouplingInflation, rightFim.fdRefCouplingInflation);
  row.relEfimFroDiff = localRelMatrixDiff(localExtractEfim(leftCase), localExtractEfim(rightCase));
  row.auditClass = localClassifyPairAudit(row, config);
  rowList(iPair) = row;
end
pairTable = struct2table(rowList);
end

function reuseTable = localBuildSummaryReuseAuditTable(pairTable, config)
%LOCALBUILDSUMMARYREUSEAUDITTABLE Keep only likely summary or identity suspects.

if isempty(pairTable)
  reuseTable = pairTable;
  return;
end
suspectMask = pairTable.auditClass == "summary-reuse-suspect" | ...
  pairTable.auditClass == "reported-fdref-crb-mismatch" | ...
  pairTable.auditClass == "unexpected-identical";
reuseTable = pairTable(suspectMask, :);
end

function contextSummaryTable = localBuildContextSummaryTable(context, periodicFixture)
%LOCALBUILDCONTEXTSUMMARYTABLE Build a compact scenario summary.

truth = periodicFixture.truth;
contextSummaryTable = table( ...
  string(context.tleFileName), ...
  string(context.utcRef), ...
  string(mat2str(context.usrLla(:).')), ...
  string(mat2str(context.selectedSatIdxGlobal)), ...
  string(mat2str(context.nonRefSatIdxGlobal)), ...
  truth.refSatIdxGlobal, ...
  truth.refSatIdxLocal, ...
  periodicFixture.sceneSeq.numSat, ...
  periodicFixture.sceneSeq.numFrame, ...
  string(mat2str(periodicFixture.sceneSeq.timeOffsetSec(:).')), ...
  truth.fdRefTrueHz, ...
  truth.fdRefFit, ...
  truth.fdRateFit, ...
  periodicFixture.pwrNoise, ...
  'VariableNames', { ...
    'tleFileName', 'utcRef', 'usrLla', 'selectedSatIdxGlobal', ...
    'nonRefSatIdxGlobal', 'refSatIdxGlobal', 'refSatIdxLocal', ...
    'numSat', 'numFrame', 'timeOffsetSec', 'fdRefTrueHz', 'fdRefFitHz', ...
    'fdRateFitHzPerSec', 'noiseVar'});
end


function satContributionTable = localBuildSatContributionAuditTable(derivTable, context)
%LOCALBUILDSATCONTRIBUTIONAUDITTABLE Summarize per-satellite CRB derivative contribution.

rowTemplate = localEmptySatContributionRow();
if isempty(derivTable) || ~istable(derivTable) || height(derivTable) == 0
  satContributionTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
caseList = unique(derivTable.caseName, 'stable');
rowList = repmat(rowTemplate, 0, 1);
for iCase = 1:numel(caseList)
  caseMask = derivTable.caseName == caseList(iCase);
  caseTable = derivTable(caseMask, :);
  satList = unique(caseTable.satIdxLocal(isfinite(caseTable.satIdxLocal)), 'stable');
  caseDoaEnergy = localDerivativeDoaEnergy(caseTable);
  caseFdRefEnergy = localDerivativeFdRefEnergy(caseTable);
  caseFdRateEnergy = localDerivativeFdRateEnergy(caseTable);
  for iSat = 1:numel(satList)
    satIdxLocal = satList(iSat);
    satMask = caseTable.satIdxLocal == satIdxLocal;
    row = rowTemplate;
    row.caseName = caseList(iCase);
    row.satMode = caseTable.satMode(find(satMask, 1, 'first'));
    row.frameMode = caseTable.frameMode(find(satMask, 1, 'first'));
    row.satIdxLocal = satIdxLocal;
    row.satIdxGlobal = localMapSatIdxGlobal(row.satMode, satIdxLocal, context);
    row.numRow = nnz(satMask);
    row.doaInfoProxy = localSumFinite(localDerivativeDoaEnergy(caseTable(satMask, :)));
    row.fdRefInfoProxy = localSumFinite(localDerivativeFdRefEnergy(caseTable(satMask, :)));
    row.fdRateInfoProxy = localSumFinite(localDerivativeFdRateEnergy(caseTable(satMask, :)));
    row.doaInfoProxyFraction = localSafeRatio(row.doaInfoProxy, localSumFinite(caseDoaEnergy));
    row.fdRefInfoProxyFraction = localSafeRatio(row.fdRefInfoProxy, localSumFinite(caseFdRefEnergy));
    row.fdRateInfoProxyFraction = localSafeRatio(row.fdRateInfoProxy, localSumFinite(caseFdRateEnergy));
    row.auditClass = localClassifySatContribution(row);
    rowList(end + 1, 1) = row; %#ok<AGROW>
  end
end
satContributionTable = struct2table(rowList(:));
end

function auditTable = localBuildCrbPathAuditTable(crbCaseList, periodicFixture, config)
%LOCALBUILDCRBPATHAUDITTABLE Build known-rate / unknown-rate CRB path rows.

rowList = repmat(localEmptyCrbPathAuditRow(), 0, 1);
rowList = [rowList; localBuildCrbPathPairRow(crbCaseList, periodicFixture.sceneSeqRefOnly, ...
  periodicFixture.truth, "single", "SS-MF-CP-K", "SS-MF-CP-U", config)]; %#ok<AGROW>
rowList = [rowList; localBuildCrbPathPairRow(crbCaseList, periodicFixture.sceneSeq, ...
  periodicFixture.truth, "multi", "MS-MF-CP-K", "MS-MF-CP-U", config)]; %#ok<AGROW>
auditTable = struct2table(rowList(:));
end

function aggregateTable = localBuildCrbPathAuditAggregateTable(auditTable)
%LOCALBUILDCRBPATHAUDITAGGREGATETABLE Summarize K/U CRB path health.

rowTemplate = struct('satMode', "", 'numRow', NaN, ...
  'crbMonotonicPassRate', NaN, 'efimMonotonicPassRate', NaN, ...
  'lossPsdPassRate', NaN, 'sharedBlockPassRate', NaN, ...
  'fdRefInvariantPassRate', NaN, 'fdRateInvariantPassRate', NaN, ...
  'gammaCouplingFiniteRate', NaN, 'angleCrbUnknownOverKnownMedian', NaN, ...
  'fdRefCrbUnknownOverKnownMedian', NaN, 'lossTraceRatioMedian', NaN, ...
  'gammaCouplingMedian', NaN, 'sharedBlockRelErrMax', NaN, ...
  'minEigEfimKnownMinusUnknownMin', NaN, 'minEigCrbUnknownMinusKnownMin', NaN, ...
  'primaryAuditClass', "");
if isempty(auditTable) || ~istable(auditTable) || height(auditTable) == 0
  aggregateTable = struct2table(repmat(rowTemplate, 0, 1));
  return;
end
satModeList = unique(auditTable.satMode, 'stable');
rowList = repmat(rowTemplate, numel(satModeList), 1);
for iMode = 1:numel(satModeList)
  mask = auditTable.satMode == satModeList(iMode);
  row = rowTemplate;
  row.satMode = satModeList(iMode);
  row.numRow = nnz(mask);
  row.crbMonotonicPassRate = mean(double(auditTable.crbMonotonicOk(mask)), 'omitnan');
  row.efimMonotonicPassRate = mean(double(auditTable.efimMonotonicOk(mask)), 'omitnan');
  row.lossPsdPassRate = mean(double(auditTable.lossPsdOk(mask)), 'omitnan');
  row.sharedBlockPassRate = mean(double(auditTable.sharedBlockOk(mask)), 'omitnan');
  row.fdRefInvariantPassRate = mean(double(auditTable.fdRefInvariantOk(mask)), 'omitnan');
  row.fdRateInvariantPassRate = mean(double(auditTable.fdRateInvariantOk(mask)), 'omitnan');
  row.gammaCouplingFiniteRate = mean(double(isfinite(auditTable.gammaCoupling(mask)) & ...
    auditTable.gammaCoupling(mask) > 0), 'omitnan');
  row.angleCrbUnknownOverKnownMedian = median(auditTable.angleCrbUnknownOverKnown(mask), 'omitnan');
  row.fdRefCrbUnknownOverKnownMedian = median(auditTable.fdRefCrbUnknownOverKnown(mask), 'omitnan');
  row.lossTraceRatioMedian = median(auditTable.lossTraceRatio(mask), 'omitnan');
  row.gammaCouplingMedian = median(auditTable.gammaCoupling(mask), 'omitnan');
  row.sharedBlockRelErrMax = max(auditTable.sharedBlockRelErr(mask), [], 'omitnan');
  row.minEigEfimKnownMinusUnknownMin = min(auditTable.minEigEfimKnownMinusUnknown(mask), [], 'omitnan');
  row.minEigCrbUnknownMinusKnownMin = min(auditTable.minEigCrbUnknownMinusKnown(mask), [], 'omitnan');
  row.primaryAuditClass = localClassifyCrbPathAggregate(row);
  rowList(iMode) = row;
end
aggregateTable = struct2table(rowList(:));
end

function row = localBuildCrbPathPairRow(crbCaseList, sceneSeq, truth, satMode, knownName, unknownName, config)
%LOCALBUILDCRBPATHPAIRROW Build one known / unknown CRB path audit row.

row = localEmptyCrbPathAuditRow();
knownIdx = find([crbCaseList.caseName] == string(knownName), 1);
unknownIdx = find([crbCaseList.caseName] == string(unknownName), 1);
row.snrDb = config.snrDb;
row.taskSeed = config.baseSeed;
row.satMode = string(satMode);
row.knownName = string(knownName);
row.unknownName = string(unknownName);
if isempty(knownIdx) || isempty(unknownIdx)
  row.auditClass = "missing-case";
  return;
end
knownCase = crbCaseList(knownIdx);
unknownCase = crbCaseList(unknownIdx);
row.refSatIdxLocal = localResolveRefSatIdxLocal(sceneSeq);
row.gammaParamIdx = localFindParamNameIndex(unknownCase.aux, "fdRate");

[row.angleCrbKnownDeg, row.fdRefCrbKnownHz] = localExtractCrbStdMetrics(knownCase.crb, truth);
[row.angleCrbUnknownDeg, row.fdRefCrbUnknownHz] = localExtractCrbStdMetrics(unknownCase.crb, truth);
row.angleCrbUnknownOverKnown = localSafeRatio(row.angleCrbUnknownDeg, row.angleCrbKnownDeg);
row.fdRefCrbUnknownOverKnown = localSafeRatio(row.fdRefCrbUnknownHz, row.fdRefCrbKnownHz);

fimKnown = localGetFiniteMatrixField(knownCase.aux, 'fimInterest', 3, 3);
fimUnknown = localGetFiniteMatrixField(unknownCase.aux, 'fimInterest', 3, 3);
crbKnown = localEnsureMatrixSize(knownCase.crb, 3, 3);
crbUnknown = localEnsureMatrixSize(unknownCase.crb, 3, 3);
row.minEigEfimKnownMinusUnknown = localMinSymEig(fimKnown - fimUnknown);
row.minEigCrbUnknownMinusKnown = localMinSymEig(crbUnknown - crbKnown);

schurAudit = localBuildGammaSchurAudit(unknownCase.aux, fimKnown);
row.lossMinEig = schurAudit.lossMinEig;
row.lossTraceRatio = schurAudit.lossTraceRatio;
row.gammaCoupling = schurAudit.gammaCoupling;
row.gammaInfo = schurAudit.gammaInfo;
row.sharedBlockRelErr = schurAudit.sharedBlockRelErr;

fdInvKnown = localCheckCrbReferenceStateInvariant(knownCase.aux, row.refSatIdxLocal, truth);
fdInvUnknown = localCheckCrbReferenceStateInvariant(unknownCase.aux, row.refSatIdxLocal, truth);
row.fdSatRefErrKnownHz = fdInvKnown.fdSatRefErrHz;
row.fdSatRefErrUnknownHz = fdInvUnknown.fdSatRefErrHz;
row.fdRateRefErrKnownHzPerSec = fdInvKnown.fdRateRefErrHzPerSec;
row.fdRateRefErrUnknownHzPerSec = fdInvUnknown.fdRateRefErrHzPerSec;
row.fdRefInvariantOk = fdInvKnown.fdRefInvariantOk && fdInvUnknown.fdRefInvariantOk;
row.fdRateInvariantOk = fdInvKnown.fdRateInvariantOk && fdInvUnknown.fdRateInvariantOk;

relEigTol = 1e-8;
row.crbMonotonicOk = localIsNonnegativeWithTol(row.minEigCrbUnknownMinusKnown, crbKnown, relEigTol);
row.efimMonotonicOk = localIsNonnegativeWithTol(row.minEigEfimKnownMinusUnknown, fimKnown, relEigTol);
row.lossPsdOk = localIsNonnegativeWithTol(row.lossMinEig, fimKnown, relEigTol);
row.sharedBlockOk = isfinite(row.sharedBlockRelErr) && row.sharedBlockRelErr <= 1e-6;
row.auditClass = localClassifyCrbPathRow(row);
end

function row = localEmptySatContributionRow()
%LOCALEMPTYSATCONTRIBUTIONROW Return one typed per-satellite contribution row.

row = struct('caseName', "", 'satMode', "", 'frameMode', "", ...
  'satIdxLocal', NaN, 'satIdxGlobal', NaN, 'numRow', NaN, ...
  'doaInfoProxy', NaN, 'fdRefInfoProxy', NaN, 'fdRateInfoProxy', NaN, ...
  'doaInfoProxyFraction', NaN, 'fdRefInfoProxyFraction', NaN, ...
  'fdRateInfoProxyFraction', NaN, 'auditClass', "");
end

function row = localEmptyCrbPathAuditRow()
%LOCALEMPTYCRBPATHAUDITROW Return one typed CRB path audit row.

row = struct('satMode', "", 'knownName', "", 'unknownName', "", ...
  'snrDb', NaN, 'taskSeed', NaN, 'refSatIdxLocal', NaN, 'gammaParamIdx', NaN, ...
  'angleCrbKnownDeg', NaN, 'angleCrbUnknownDeg', NaN, ...
  'angleCrbUnknownOverKnown', NaN, 'fdRefCrbKnownHz', NaN, ...
  'fdRefCrbUnknownHz', NaN, 'fdRefCrbUnknownOverKnown', NaN, ...
  'minEigEfimKnownMinusUnknown', NaN, 'minEigCrbUnknownMinusKnown', NaN, ...
  'lossMinEig', NaN, 'lossTraceRatio', NaN, 'gammaCoupling', NaN, ...
  'gammaInfo', NaN, 'sharedBlockRelErr', NaN, ...
  'fdSatRefErrKnownHz', NaN, 'fdSatRefErrUnknownHz', NaN, ...
  'fdRateRefErrKnownHzPerSec', NaN, 'fdRateRefErrUnknownHzPerSec', NaN, ...
  'crbMonotonicOk', false, 'efimMonotonicOk', false, 'lossPsdOk', false, ...
  'sharedBlockOk', false, 'fdRefInvariantOk', false, 'fdRateInvariantOk', false, ...
  'auditClass', "");
end

function rowTable = localBuildMfDerivativeRows(caseInfo)
%LOCALBUILDMFDERIVATIVEROWS Build derivative norm rows from MF CRB aux.

aux = caseInfo.aux;
numSat = localGetFieldOrDefault(aux, 'numSat', 0);
numFrame = localGetFieldOrDefault(aux, 'numFrame', 0);
rowList = repmat(localEmptyDerivativeRow(), max(numSat * numFrame, 0), 1);
rowCounter = 0;
for iFrame = 1:numFrame
  for iSat = 1:numSat
    rowCounter = rowCounter + 1;
    row = localEmptyDerivativeRow();
    row.caseName = caseInfo.caseName;
    row.satMode = caseInfo.expected.satMode;
    row.frameMode = caseInfo.expected.frameMode;
    row.satIdxLocal = iSat;
    row.frameIdx = iFrame;
    row.blockType = "mf-dAtom";
    dAtom = [];
    if isfield(aux, 'dAtomCell') && numel(aux.dAtomCell) >= sub2ind(size(aux.dAtomCell), iFrame, iSat)
      dAtom = aux.dAtomCell{iFrame, iSat};
    end
    row.normDMuDDoa1 = localColumnNorm(dAtom, 1);
    row.normDMuDDoa2 = localColumnNorm(dAtom, 2);
    row.normDMuDFdRef = localColumnNorm(dAtom, 3);
    row.normDMuDFdRate = localColumnNorm(dAtom, 4);
    row.blockFimDoa1 = NaN;
    row.blockFimDoa2 = NaN;
    row.blockFimFdRef = NaN;
    rowList(rowCounter) = row;
  end
end
if rowCounter == 0
  rowTable = struct2table(localEmptyDerivativeRow());
  rowTable(1, :) = [];
else
  rowTable = struct2table(rowList(1:rowCounter));
end
end

function rowTable = localBuildSfContributionRows(caseInfo)
%LOCALBUILDSFCONTRIBUTIONROWS Build SF DoA-only per-satellite contribution rows.

aux = caseInfo.aux;
numSat = localGetFieldOrDefault(aux, 'numArray', caseInfo.actual.numSat);
rowList = repmat(localEmptyDerivativeRow(), max(numSat, 0), 1);
for iSat = 1:numSat
  row = localEmptyDerivativeRow();
  row.caseName = caseInfo.caseName;
  row.satMode = caseInfo.expected.satMode;
  row.frameMode = caseInfo.expected.frameMode;
  row.satIdxLocal = iSat;
  row.frameIdx = 1;
  row.blockType = "sf-fimCell";
  fimCell = [];
  if isfield(aux, 'fimCell') && numel(aux.fimCell) >= iSat
    fimCell = aux.fimCell{iSat};
  end
  row.blockFimDoa1 = localSafeMatrixElem(fimCell, 1, 1);
  row.blockFimDoa2 = localSafeMatrixElem(fimCell, 2, 2);
  row.blockFimFdRef = NaN;
  rowList(iSat) = row;
end
if numSat == 0
  rowTable = struct2table(localEmptyDerivativeRow());
  rowTable(1, :) = [];
else
  rowTable = struct2table(rowList);
end
end

function [rawFim, efim, nuisanceFim, crossFim] = localExtractFimBlocks(caseInfo)
%LOCALEXTRACTFIMBLOCKS Extract raw interest FIM, EFIM, nuisance FIM and cross block.

aux = caseInfo.aux;
rawFim = [];
efim = [];
nuisanceFim = [];
crossFim = [];
if caseInfo.expected.isMf
  fimFull = localGetFieldOrDefault(aux, 'fimFull', []);
  if ~isempty(fimFull)
    interestIdx = 1:min(3, size(fimFull, 1));
    rawFim = fimFull(interestIdx, interestIdx);
    nuisanceIdx = localGetFieldOrDefault(aux, 'nuisanceParamIdx', []);
    nuisanceIdx = nuisanceIdx(nuisanceIdx >= 1 & nuisanceIdx <= size(fimFull, 1));
    if ~isempty(nuisanceIdx)
      nuisanceFim = fimFull(nuisanceIdx, nuisanceIdx);
      crossFim = fimFull(interestIdx, nuisanceIdx);
    end
  end
  efim = localGetFieldOrDefault(aux, 'fimInterest', []);
else
  rawFim = localGetFieldOrDefault(aux, 'fim', []);
  efim = rawFim;
end
end

function efim = localExtractEfim(caseInfo)
%LOCALEXTRACTEFIM Extract the effective FIM for pair comparison.

[~, efim] = localExtractFimBlocks(caseInfo);
end

function expected = localExpectedModelSpec(methodName, config)
%LOCALEXPECTEDMODELSPEC Return the paper-level expected model tags.

methodName = string(methodName);
expected = struct();
expected.caseName = methodName;
if startsWith(methodName, "MS-")
  expected.numSat = config.numSatMulti;
else
  expected.numSat = config.numSatSingle;
end
expected.isMf = contains(methodName, "-MF-");
expected.numFrame = 1;
expected.satMode = "single";
if startsWith(methodName, "MS-")
  expected.satMode = "multi";
end
expected.frameMode = "single";
expected.numInterest = 2;
expected.modelTag = "sfDoaOnly";
expected.signalModel = "single-frame DoA-only effective pilot model";
expected.crbFunction = "crbPilotSfDoaOnlyEffective";
expected.phaseMode = "doa-only";
expected.fdRateMode = "none";
expected.steeringMode = "single-frame";
expected.timeTemplate = "no-doppler-estimation";
expected.interestNames = "lat,lon";
expected.nuisanceNames = "deterministic-signal-power-profiled-in-crbDetDoa";

if expected.isMf
  expected.numFrame = NaN;
  expected.frameMode = "multi";
  expected.numInterest = 3;
  expected.signalModel = "multi-frame DoA-Doppler pilot model";
  expected.crbFunction = "crbPilotMfDoaDoppler";
  expected.phaseMode = "continuous";
  expected.steeringMode = "framewise";
  expected.timeTemplate = "absolute-time continuous phase";
  expected.interestNames = "lat,lon,fdRef";
  expected.nuisanceNames = "phase,amplitude";
  if contains(methodName, "Static")
    expected.modelTag = "mfStaticZeroRate";
    expected.fdRateMode = "zero";
    expected.signalModel = "multi-frame zero-rate continuous-phase DoA-Doppler model";
  elseif contains(methodName, "CP-K")
    expected.modelTag = "mfCpKnownRate";
    expected.fdRateMode = "known";
    expected.signalModel = "multi-frame affine-Doppler continuous-phase known-rate model";
  elseif contains(methodName, "CP-U")
    expected.modelTag = "mfCpUnknownRate";
    expected.fdRateMode = "unknown";
    expected.signalModel = "multi-frame affine-Doppler continuous-phase unknown-rate model";
    expected.nuisanceNames = "fdRateRef,phase,amplitude";
  end
end
end

function actual = localActualModelSpec(caseInfo)
%LOCALACTUALMODELSPEC Parse actual model tags from CRB aux fields.

aux = caseInfo.aux;
expected = caseInfo.expected;
actual = struct();
actual.phaseMode = localReadStringField(aux, 'phaseMode', expected.phaseMode);
actual.fdRateMode = localReadStringField(aux, 'fdRateMode', expected.fdRateMode);
actual.steeringMode = localReadStringField(aux, 'steeringMode', expected.steeringMode);
actual.numSat = localGetFieldOrDefault(aux, 'numSat', expected.numSat);
if ~expected.isMf && isfield(aux, 'numArray')
  actual.numSat = aux.numArray;
end
actual.numFrame = localGetFieldOrDefault(aux, 'numFrame', expected.numFrame);
if ~expected.isMf
  actual.numFrame = 1;
end
actual.timeTemplate = localInferTimeTemplate(actual.phaseMode, actual.fdRateMode, expected);
actual.paramNameFull = localStringList(localGetFieldOrDefault(aux, 'paramNameFull', []));
actual.interestNames = localStringList(localGetFieldOrDefault(aux, 'paramNameInterest', ...
  localGetFieldOrDefault(aux, 'paramName', [])));
if strlength(actual.interestNames) == 0
  actual.interestNames = expected.interestNames;
end
nuisanceIdx = localGetFieldOrDefault(aux, 'nuisanceParamIdx', []);
if ~isempty(nuisanceIdx) && isfield(aux, 'paramNameFull')
  paramNameFull = string(aux.paramNameFull);
  nuisanceIdx = nuisanceIdx(nuisanceIdx >= 1 & nuisanceIdx <= numel(paramNameFull));
  actual.nuisanceNames = localStringList(paramNameFull(nuisanceIdx));
else
  actual.nuisanceNames = "";
end
actual.numInterest = localCountDelimitedNames(actual.interestNames);
if expected.isMf && isfield(aux, 'fimInterest') && ~isempty(aux.fimInterest)
  actual.numInterest = size(aux.fimInterest, 1);
elseif ~expected.isMf && isfield(aux, 'fim') && ~isempty(aux.fim)
  actual.numInterest = size(aux.fim, 1);
end
actual.numNuisance = localCountDelimitedNames(actual.nuisanceNames);
end

function timeTemplate = localInferTimeTemplate(phaseMode, fdRateMode, expected)
%LOCALINFERTIMETEMPLATE Infer a compact time-template label for display.

if ~expected.isMf
  timeTemplate = expected.timeTemplate;
  return;
end
switch string(phaseMode)
  case {"continuous", "relaxed"}
    timeTemplate = "absolute-time continuous phase";
  case "independent"
    timeTemplate = "frame-local independent phase";
  otherwise
    timeTemplate = "unknown";
end
if string(fdRateMode) == "zero"
  timeTemplate = timeTemplate + " / zero rate";
end
end

function pairSpec = localCasePairSpec()
%LOCALCASEPAIRSPEC Return key case-pair comparisons for this audit.

pairSpec = [
  struct('leftCase', "SS-MF-Static", 'rightCase', "SS-MF-CP-K", 'expectedRelation', "zero-rate-degenerate-continuous-static")
  struct('leftCase', "MS-MF-Static", 'rightCase', "MS-MF-CP-K", 'expectedRelation', "zero-rate-degenerate-continuous-static")
  struct('leftCase', "SS-MF-CP-K", 'rightCase', "SS-MF-CP-U", 'expectedRelation', "known-vs-unknown-efim-loss")
  struct('leftCase', "MS-MF-CP-K", 'rightCase', "MS-MF-CP-U", 'expectedRelation', "known-vs-unknown-efim-loss")
  struct('leftCase', "SS-MF-CP-K", 'rightCase', "MS-MF-CP-K", 'expectedRelation', "single-vs-multi-sat")
  struct('leftCase', "SS-MF-CP-U", 'rightCase', "MS-MF-CP-U", 'expectedRelation', "single-vs-multi-sat")
  struct('leftCase', "SS-MF-Static", 'rightCase', "MS-MF-Static", 'expectedRelation', "single-vs-multi-sat")
  struct('leftCase', "SS-SF-DoA", 'rightCase', "MS-SF-DoA", 'expectedRelation', "single-vs-multi-sat-doa-only")
];
end

function auditClass = localClassifyPairAudit(row, config)
%LOCALCLASSIFYPAIRAUDIT Classify one pair comparison.

isIdentical = row.relRawFimFdDiff <= config.identicalRelTol && ...
  row.relEfimFdDiff <= config.identicalRelTol && ...
  row.relCrbFdRefDiff <= config.identicalRelTol && ...
  row.relEfimFroDiff <= config.identicalRelTol;
reportedMismatchTol = sqrt(config.summaryReuseRelTol);
if row.relReportedFullInvLeft > reportedMismatchTol || ...
    row.relReportedFullInvRight > reportedMismatchTol
  auditClass = "reported-fdref-crb-mismatch";
elseif contains(row.expectedRelation, "zero-rate-degenerate")
  if isIdentical
    auditClass = "expected-continuous-static-degenerate";
  else
    auditClass = "observed-static-cp-difference";
  end
elseif contains(row.expectedRelation, "known-vs-unknown")
  if row.relEfimFroDiff <= config.identicalRelTol
    auditClass = "unexpected-identical";
  else
    auditClass = "observed-known-unknown-difference";
  end
elseif contains(row.expectedRelation, "single-vs-multi")
  if row.relEfimFroDiff <= config.identicalRelTol
    auditClass = "unexpected-identical";
  elseif row.relScalarEfimFdRefStdDiff > reportedMismatchTol && ...
      row.relFullInvFdRefDiff <= reportedMismatchTol && ...
      row.relCrbFdRefDiff <= reportedMismatchTol
    auditClass = "fdref-coupling-cancelled-gain";
  else
    auditClass = "observed-ss-ms-difference";
  end
elseif isIdentical
  auditClass = "unexpected-identical";
elseif row.relCrbFdRefDiff <= config.summaryReuseRelTol && ...
    row.relFullInvFdRefDiff > reportedMismatchTol
  auditClass = "summary-reuse-suspect";
elseif row.relRawFimFdDiff > sqrt(config.identicalRelTol) && ...
    row.relEfimFdDiff <= config.identicalRelTol
  auditClass = "schur-collapse-suspect";
else
  auditClass = "observed-difference";
end
end

function replayData = localValidateReplayData(replayData)
%LOCALVALIDATEREPLAYDATA Validate summary tables for replay section reruns.

requiredField = {'crbModelAuditTable', 'crbFimBlockAuditTable', ...
  'crbFdRefMarginalAuditTable', 'crbEigenCouplingAuditTable', ...
  'crbDerivativeNormTable', 'crbSatContributionAuditTable', ...
  'crbPathAuditTable', 'crbPathAuditAggregateTable', 'crbCasePairAuditTable'};
for iField = 1:numel(requiredField)
  fieldName = requiredField{iField};
  if ~isfield(replayData, fieldName) || ~istable(replayData.(fieldName))
    error('replayMfCrbTexConsistencyAudit:MissingReplayTable', ...
      'replayData.%s is missing or is not a table.', fieldName);
  end
end
end

function metricLineList = localBuildTelegramMetricLines(replayData)
%LOCALBUILDTELEGRAMMETRICLINES Build short notification metrics.

pairTable = replayData.crbCasePairAuditTable;
numUnexpected = sum(pairTable.auditClass == "unexpected-identical");
numSummarySuspect = sum(pairTable.auditClass == "summary-reuse-suspect");
metricLineList = [
  "• CRB cases: <code>" + string(height(replayData.crbModelAuditTable)) + "</code>"
  "• unexpected-identical pairs: <code>" + string(numUnexpected) + "</code>"
  "• summary-reuse suspects: <code>" + string(numSummarySuspect) + "</code>"
  "• max EFIM cond: <code>" + string(localMaxFinite(replayData.crbEigenCouplingAuditTable.efimCondEig)) + "</code>"
];
end

function caseInfo = localEmptyCrbCase()
%LOCALEMPTYCRBCASE Build one empty CRB case struct.

caseInfo = struct();
caseInfo.caseName = "";
caseInfo.expected = struct();
caseInfo.actual = struct();
caseInfo.crb = [];
caseInfo.aux = struct();
caseInfo.errorMessage = "";
end

function row = localEmptyModelAuditRow()
%LOCALEMPTYMODELAUDITROW Return one empty model-audit row.

row = struct();
row.caseName = "";
row.crbFunction = "";
row.expectedSignalModel = "";
row.satMode = "";
row.frameMode = "";
row.expectedTimeTemplate = "";
row.actualTimeTemplate = "";
row.expectedPhaseMode = "";
row.actualPhaseMode = "";
row.expectedFdRateMode = "";
row.actualFdRateMode = "";
row.expectedSteeringMode = "";
row.actualSteeringMode = "";
row.numSat = NaN;
row.numFrame = NaN;
row.numInterest = NaN;
row.numNuisance = NaN;
row.interestNames = "";
row.nuisanceNames = "";
row.dispatchOk = false;
row.errorMessage = "";
end

function row = localEmptyFimBlockRow()
%LOCALEMPTYFIMBLOCKROW Return one empty FIM-block row.

row = struct();
row.caseName = "";
row.satMode = "";
row.frameMode = "";
row.modelTag = "";
row.crbDim = NaN;
row.rawFimDoaTrace = NaN;
row.rawFimFdFd = NaN;
row.rawFimRateRate = NaN;
row.efimDoaTrace = NaN;
row.efimFdFd = NaN;
row.crossNorm = NaN;
row.nuisanceCond = NaN;
row.efimCond = NaN;
row.crbAngleTraceDeg = NaN;
row.crbFdRefStdHz = NaN;
row.fdRefStdFullInvHz = NaN;
row.fdRefStdScalarEfimHz = NaN;
row.fdRefStdRawInvHz = NaN;
row.fdRefInfoFullInv = NaN;
row.fdRefScalarInfoEfim = NaN;
row.fdRefCouplingInflation = NaN;
row.fdRefCouplingLossFraction = NaN;
row.crbFdRefReportedRelDiff = NaN;
row.status = "";
end

function row = localEmptyEigenCouplingRow()
%LOCALEMPTYEIGENCOUPLINGROW Return one empty eigen/coupling audit row.

row = struct();
row.caseName = "";
row.satMode = "";
row.frameMode = "";
row.modelTag = "";
row.eigInfoMin = NaN;
row.eigInfoMid = NaN;
row.eigInfoMax = NaN;
row.efimCondEig = NaN;
row.stdAlongMinInfoAxis = NaN;
row.stdAlongMidInfoAxis = NaN;
row.stdAlongMaxInfoAxis = NaN;
row.fdRefWeightMinInfoAxis = NaN;
row.doaWeightMinInfoAxis = NaN;
row.fdRefWeightMaxInfoAxis = NaN;
row.doaWeightMaxInfoAxis = NaN;
row.doaBlockCond = NaN;
row.corrDoa1FdRef = NaN;
row.corrDoa2FdRef = NaN;
row.fdRefInfoDoaKnown = NaN;
row.fdRefInfoMarginalSchur = NaN;
row.fdRefInfoMarginalFullInv = NaN;
row.fdRefDoaCouplingInfo = NaN;
row.fdRefDoaCouplingFraction = NaN;
row.fdRefDoaCanonicalCorr = NaN;
row.fdRefStdDoaKnownHz = NaN;
row.fdRefStdDoaUnknownHz = NaN;
row.fdRefStdFullInvHz = NaN;
row.fdRefMarginalInfoRelDiff = NaN;
row.status = "";
end

function row = localEmptyDerivativeRow()
%LOCALEMPTYDERIVATIVEROW Return one empty derivative/contribution row.

row = struct();
row.caseName = "";
row.satMode = "";
row.frameMode = "";
row.satIdxLocal = NaN;
row.frameIdx = NaN;
row.blockType = "";
row.normDMuDDoa1 = NaN;
row.normDMuDDoa2 = NaN;
row.normDMuDFdRef = NaN;
row.normDMuDFdRate = NaN;
row.blockFimDoa1 = NaN;
row.blockFimDoa2 = NaN;
row.blockFimFdRef = NaN;
end

function row = localEmptyCasePairRow()
%LOCALEMPTYCASEPAIRROW Return one empty case-pair row.

row = struct();
row.leftCase = "";
row.rightCase = "";
row.expectedRelation = "";
row.sameTimeTemplateFlag = false;
row.sameInterestFlag = false;
row.sameNuisanceFlag = false;
row.relRawFimFdDiff = NaN;
row.relEfimFdDiff = NaN;
row.relCrbFdRefDiff = NaN;
row.relFullInvFdRefDiff = NaN;
row.relScalarEfimFdRefStdDiff = NaN;
row.relReportedFullInvLeft = NaN;
row.relReportedFullInvRight = NaN;
row.relCouplingInflationDiff = NaN;
row.relEfimFroDiff = NaN;
row.auditClass = "";
end


function energyValue = localDerivativeDoaEnergy(inputTable)
%LOCALDERIVATIVEDOAENERGY Build a DoA information proxy from derivative rows.

if isempty(inputTable) || height(inputTable) == 0
  energyValue = NaN(0, 1);
  return;
end
energyValue = inputTable.normDMuDDoa1.^2 + inputTable.normDMuDDoa2.^2;
sfProxy = inputTable.blockFimDoa1 + inputTable.blockFimDoa2;
useSf = ~isfinite(energyValue) & isfinite(sfProxy);
energyValue(useSf) = sfProxy(useSf);
end

function energyValue = localDerivativeFdRefEnergy(inputTable)
%LOCALDERIVATIVEFDREFENERGY Build an fdRef information proxy from derivative rows.

if isempty(inputTable) || height(inputTable) == 0
  energyValue = NaN(0, 1);
  return;
end
energyValue = inputTable.normDMuDFdRef.^2;
useSf = ~isfinite(energyValue) & isfinite(inputTable.blockFimFdRef);
energyValue(useSf) = inputTable.blockFimFdRef(useSf);
end

function energyValue = localDerivativeFdRateEnergy(inputTable)
%LOCALDERIVATIVEFDRATEENERGY Build an fdRate information proxy from derivative rows.

if isempty(inputTable) || height(inputTable) == 0
  energyValue = NaN(0, 1);
else
  energyValue = inputTable.normDMuDFdRate.^2;
end
end

function sumValue = localSumFinite(valueArray)
%LOCALSUMFINITE Sum finite nonnegative values and return NaN if none exist.

valueArray = valueArray(isfinite(valueArray) & valueArray >= 0);
if isempty(valueArray)
  sumValue = NaN;
else
  sumValue = sum(valueArray);
end
end

function satIdxGlobal = localMapSatIdxGlobal(satMode, satIdxLocal, context)
%LOCALMAPSATIDXGLOBAL Map local satellite index to global satellite index.

satIdxGlobal = NaN;
satIdxLocal = round(double(satIdxLocal));
if ~(isfinite(satIdxLocal) && satIdxLocal >= 1)
  return;
end
if string(satMode) == "single"
  satIdxGlobal = double(context.refSatIdxGlobal);
elseif numel(context.selectedSatIdxGlobal) >= satIdxLocal
  satIdxGlobal = double(context.selectedSatIdxGlobal(satIdxLocal));
end
end

function classText = localClassifySatContribution(row)
%LOCALCLASSIFYSATCONTRIBUTION Classify one per-satellite contribution row.

if ~(isfinite(row.numRow) && row.numRow > 0)
  classText = "missing-sat";
elseif ~(isfinite(row.doaInfoProxyFraction) || isfinite(row.fdRefInfoProxyFraction) || ...
    isfinite(row.fdRateInfoProxyFraction))
  classText = "missing-proxy";
elseif max([row.doaInfoProxy, row.fdRefInfoProxy, row.fdRateInfoProxy], [], 'omitnan') <= 0
  classText = "zero-contribution";
else
  classText = "sat-ok";
end
end

function [angleCrbStdDeg, fdRefCrbStdHz] = localExtractCrbStdMetrics(crbFull, truth)
%LOCALEXTRACTCRBSTDMETRICS Extract spherical angle and fdRef CRB std.

crbFull = localEnsureMatrixSize(crbFull, 3, 3);
angleCrbStdDeg = NaN;
fdRefCrbStdHz = NaN;
if all(isfinite(crbFull(:)))
  [angleCrbStdDeg, ~] = projectCrbToAngleMetric(real(crbFull(1:2, 1:2)), truth.latlonTrueDeg, 'latlon');
  fdRefCrbStdHz = sqrt(max(real(crbFull(3, 3)), 0));
end
end

function schurAudit = localBuildGammaSchurAudit(auxUnknown, fimKnown)
%LOCALBUILDGAMMASCHURAUDIT Rebuild gamma Schur loss after common nuisance removal.

schurAudit = struct('lossMinEig', NaN, 'lossTraceRatio', NaN, ...
  'gammaCoupling', NaN, 'gammaInfo', NaN, 'sharedBlockRelErr', NaN);
fimFull = localGetFieldOrDefault(auxUnknown, 'fimFull', []);
gammaIdx = localFindParamNameIndex(auxUnknown, "fdRate");
if isempty(fimFull) || ~isnumeric(fimFull) || ~isfinite(gammaIdx) || gammaIdx < 4
  return;
end
fimFull = real(fimFull);
numParam = size(fimFull, 1);
if size(fimFull, 2) ~= numParam || numParam < gammaIdx || numParam < 4
  return;
end
interestGammaIdx = [1, 2, 3, gammaIdx];
commonNuisanceIdx = setdiff(1:numParam, interestGammaIdx);
reducedFim = localSchurReduce(fimFull, interestGammaIdx, commonNuisanceIdx);
if size(reducedFim, 1) ~= 4 || any(~isfinite(reducedFim(:)))
  return;
end
etaEta = reducedFim(1:3, 1:3);
etaGamma = reducedFim(1:3, 4);
gammaGamma = reducedFim(4, 4);
if isfinite(gammaGamma) && abs(gammaGamma) > eps(max(1, norm(reducedFim, 'fro')))
  lossMat = (etaGamma * etaGamma') ./ gammaGamma;
  schurAudit.lossMinEig = localMinSymEig(lossMat);
  schurAudit.lossTraceRatio = localSafeRatio(real(trace(lossMat)), real(trace(etaEta)));
  schurAudit.gammaCoupling = norm(etaGamma) ./ sqrt(max(real(trace(etaEta)), 0) * abs(gammaGamma));
  schurAudit.gammaInfo = gammaGamma;
end
if all(isfinite(fimKnown(:))) && any(fimKnown(:) ~= 0)
  schurAudit.sharedBlockRelErr = norm(etaEta - fimKnown, 'fro') ./ max(norm(fimKnown, 'fro'), eps);
end
end

function reducedFim = localSchurReduce(fimFull, keepIdx, nuisanceIdx)
%LOCALSCHURREDUCE Eliminate nuisance indices from a symmetric FIM.

fimFull = localSymPart(real(fimFull));
keepIdx = reshape(double(keepIdx), 1, []);
nuisanceIdx = reshape(double(nuisanceIdx), 1, []);
if isempty(nuisanceIdx)
  reducedFim = fimFull(keepIdx, keepIdx);
  return;
end
fimKeepKeep = fimFull(keepIdx, keepIdx);
fimKeepNuisance = fimFull(keepIdx, nuisanceIdx);
fimNuisanceNuisance = fimFull(nuisanceIdx, nuisanceIdx);
reducedFim = fimKeepKeep - fimKeepNuisance * (pinv(fimNuisanceNuisance) * fimKeepNuisance');
reducedFim = localSymPart(reducedFim);
end

function invInfo = localCheckCrbReferenceStateInvariant(aux, refSatIdxLocal, truth)
%LOCALCHECKCRBREFERENCESTATEINVARIANT Check CRB-side reference Doppler state.

invInfo = struct('fdSatRefErrHz', NaN, 'fdRateRefErrHzPerSec', NaN, ...
  'fdRefInvariantOk', false, 'fdRateInvariantOk', false);
if ~isfinite(refSatIdxLocal) || refSatIdxLocal < 1
  return;
end
refIdx = round(refSatIdxLocal);
fdSat = reshape(localGetFieldOrDefault(aux, 'fdSat', []), [], 1);
fdRateSat = reshape(localGetFieldOrDefault(aux, 'fdRateSat', []), [], 1);
truthFdRefHz = localGetTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = localGetTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if numel(fdSat) >= refIdx && isfinite(fdSat(refIdx)) && isfinite(truthFdRefHz)
  invInfo.fdSatRefErrHz = fdSat(refIdx) - truthFdRefHz;
  fdTol = max(1e-5, 1e-10 * max(1, abs(truthFdRefHz)));
  invInfo.fdRefInvariantOk = abs(invInfo.fdSatRefErrHz) <= fdTol;
end
if numel(fdRateSat) >= refIdx && isfinite(fdRateSat(refIdx)) && isfinite(truthFdRateHzPerSec)
  invInfo.fdRateRefErrHzPerSec = fdRateSat(refIdx) - truthFdRateHzPerSec;
  fdRateTol = max(1e-5, 1e-10 * max(1, abs(truthFdRateHzPerSec)));
  invInfo.fdRateInvariantOk = abs(invInfo.fdRateRefErrHzPerSec) <= fdRateTol;
end
end

function scalarValue = localGetTruthScalar(truth, fieldNameList)
%LOCALGETTRUTHSCALAR Read the first finite scalar truth value from a field list.

scalarValue = NaN;
for iField = 1:numel(fieldNameList)
  value = localGetFieldOrDefault(truth, fieldNameList{iField}, []);
  if isnumeric(value) && ~isempty(value) && isfinite(value(1))
    scalarValue = double(value(1));
    return;
  end
end
end

function refSatIdxLocal = localResolveRefSatIdxLocal(sceneSeq)
%LOCALRESOLVEREFSATIDXLOCAL Resolve local reference-satellite index for an audit row.

refSatIdxLocal = NaN;
if isstruct(sceneSeq) && isfield(sceneSeq, 'refSatIdxLocal') && ~isempty(sceneSeq.refSatIdxLocal)
  refSatIdxLocal = double(sceneSeq.refSatIdxLocal(1));
elseif isstruct(sceneSeq) && isfield(sceneSeq, 'ref') && isfield(sceneSeq.ref, 'satIdxLocal') && ~isempty(sceneSeq.ref.satIdxLocal)
  refSatIdxLocal = double(sceneSeq.ref.satIdxLocal(1));
elseif isstruct(sceneSeq) && isfield(sceneSeq, 'sceneCell') && ~isempty(sceneSeq.sceneCell)
  refFrameIdx = round(localGetFieldOrDefault(sceneSeq, 'refFrameIdx', 1));
  refFrameIdx = min(max(refFrameIdx, 1), numel(sceneSeq.sceneCell));
  sceneRef = sceneSeq.sceneCell{refFrameIdx};
  if isstruct(sceneRef) && isfield(sceneRef, 'ref') && isfield(sceneRef.ref, 'satIdxLocal') && ~isempty(sceneRef.ref.satIdxLocal)
    refSatIdxLocal = double(sceneRef.ref.satIdxLocal(1));
  end
end
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && refSatIdxLocal >= 1)
  refSatIdxLocal = 1;
end
end

function idx = localFindParamNameIndex(aux, paramName)
%LOCALFINDPARAMNAMEINDEX Find a full-parameter index by exact or contains match.

idx = NaN;
nameList = string(localGetFieldOrDefault(aux, 'paramNameFull', strings(0, 1)));
if isempty(nameList)
  return;
end
matchIdx = find(nameList == string(paramName), 1, 'first');
if isempty(matchIdx)
  matchIdx = find(contains(lower(nameList), lower(string(paramName))), 1, 'first');
end
if ~isempty(matchIdx)
  idx = matchIdx;
end
end

function matrixValue = localGetFiniteMatrixField(dataStruct, fieldName, numRow, numCol)
%LOCALGETFINITEMATRIXFIELD Read one finite matrix field or return NaNs.

matrixValue = NaN(numRow, numCol);
value = localGetFieldOrDefault(dataStruct, fieldName, []);
if isnumeric(value) && isequal(size(value), [numRow, numCol])
  matrixValue = real(value);
end
end

function matrixValue = localEnsureMatrixSize(matrixValue, numRow, numCol)
%LOCALENSUREMATRIXSIZE Return a matrix with the requested size or NaNs.

if ~(isnumeric(matrixValue) && isequal(size(matrixValue), [numRow, numCol]))
  matrixValue = NaN(numRow, numCol);
else
  matrixValue = real(matrixValue);
end
end

function minEig = localMinSymEig(matrixValue)
%LOCALMINSYMEIG Return the smallest eigenvalue of a symmetric matrix.

if isempty(matrixValue) || any(~isfinite(matrixValue(:))) || size(matrixValue, 1) ~= size(matrixValue, 2)
  minEig = NaN;
  return;
end
minEig = min(real(eig(localSymPart(matrixValue))));
end

function matrixValue = localSymPart(matrixValue)
%LOCALSYMPART Force a matrix to its real symmetric part.

matrixValue = 0.5 * (real(matrixValue) + real(matrixValue)');
end

function tf = localIsNonnegativeWithTol(minEig, scaleMatrix, relTol)
%LOCALISNONNEGATIVEWITHTOL Check PSD with a relative numerical tolerance.

if ~(isfinite(minEig) && isnumeric(scaleMatrix))
  tf = false;
  return;
end
scaleValue = max(1, norm(real(scaleMatrix), 'fro'));
tf = minEig >= -relTol * scaleValue;
end

function classText = localClassifyCrbPathAggregate(row)
%LOCALCLASSIFYCRBPATHAGGREGATE Classify compact CRB path audit outcome.

if ~(isfinite(row.crbMonotonicPassRate) && isfinite(row.efimMonotonicPassRate) && ...
    isfinite(row.lossPsdPassRate) && isfinite(row.sharedBlockPassRate) && ...
    isfinite(row.fdRefInvariantPassRate) && isfinite(row.fdRateInvariantPassRate) && ...
    isfinite(row.gammaCouplingFiniteRate))
  classText = "audit-missing";
elseif row.crbMonotonicPassRate < 1
  classText = "crb-monotonic-fail";
elseif row.efimMonotonicPassRate < 1
  classText = "efim-monotonic-fail";
elseif row.lossPsdPassRate < 1
  classText = "schur-loss-fail";
elseif row.sharedBlockPassRate < 1
  classText = "shared-block-mismatch";
elseif row.fdRefInvariantPassRate < 1 || row.fdRateInvariantPassRate < 1
  classText = "reference-state-invariant-fail";
elseif row.gammaCouplingFiniteRate < 1
  classText = "gamma-coupling-weak-or-missing";
else
  classText = "path-ok";
end
end

function classText = localClassifyCrbPathRow(row)
%LOCALCLASSIFYCRBPATHROW Classify one K/U CRB path audit row.

if string(row.auditClass) == "missing-case"
  classText = "missing-case";
elseif ~(isfinite(row.minEigCrbUnknownMinusKnown) && isfinite(row.minEigEfimKnownMinusUnknown) && ...
    isfinite(row.lossMinEig) && isfinite(row.sharedBlockRelErr))
  classText = "audit-missing";
elseif ~row.crbMonotonicOk
  classText = "crb-monotonic-fail";
elseif ~row.efimMonotonicOk
  classText = "efim-monotonic-fail";
elseif ~row.lossPsdOk
  classText = "schur-loss-fail";
elseif ~row.sharedBlockOk
  classText = "shared-block-mismatch";
elseif ~(row.fdRefInvariantOk && row.fdRateInvariantOk)
  classText = "reference-state-invariant-fail";
elseif ~(isfinite(row.gammaCoupling) && row.gammaCoupling > 0)
  classText = "gamma-coupling-weak-or-missing";
else
  classText = "path-ok";
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a struct field or return a default.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end

function value = localReadStringField(dataStruct, fieldName, defaultValue)
%LOCALREADSTRINGFIELD Read a scalar string field with default.

value = string(defaultValue);
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = string(dataStruct.(fieldName));
  if ~isscalar(value)
    value = strjoin(cellstr(value(:).'), ',');
  end
end
end

function nameText = localStringList(nameValue)
%LOCALSTRINGLIST Convert a name vector to one comma-separated string.

if isempty(nameValue)
  nameText = "";
  return;
end
nameVec = string(nameValue);
nameText = strjoin(cellstr(nameVec(:).'), ',');
end

function countValue = localCountDelimitedNames(nameText)
%LOCALCOUNTDELIMITEDNAMES Count comma-delimited names in a string.

nameText = strtrim(string(nameText));
if strlength(nameText) == 0
  countValue = 0;
else
  countValue = numel(split(nameText, ','));
end
end

function isEqual = localStringEq(leftValue, rightValue)
%LOCALSTRINGEQ Compare two scalar strings.

isEqual = localSameScalarString(leftValue, rightValue);
end

function isEqual = localSameScalarString(leftValue, rightValue)
%LOCALSAMESCALARSTRING Compare scalar strings after trimming.

leftValue = strtrim(string(leftValue));
rightValue = strtrim(string(rightValue));
isEqual = isscalar(leftValue) && isscalar(rightValue) && leftValue == rightValue;
end

function isEqual = localSameNameList(leftValue, rightValue)
%LOCALSAMENAMELIST Compare comma-delimited or string-vector name lists.

leftName = localNormalizeNameList(leftValue);
rightName = localNormalizeNameList(rightValue);
isEqual = isequal(leftName, rightName);
end

function nameList = localNormalizeNameList(nameValue)
%LOCALNORMALIZENAMELIST Normalize comma-delimited names for exact comparison.

if isempty(nameValue)
  nameList = strings(0, 1);
  return;
end
nameValue = string(nameValue);
if numel(nameValue) == 1
  if strlength(strtrim(nameValue)) == 0
    nameList = strings(0, 1);
    return;
  end
  nameList = split(nameValue, ',');
else
  nameList = nameValue(:);
end
nameList = strtrim(nameList(:));
nameList = nameList(strlength(nameList) > 0);
end

function dimValue = localMatrixDim(matValue)
%LOCALMATRIXDIM Return matrix dimension for compact display.

if isempty(matValue) || ~isnumeric(matValue)
  dimValue = NaN;
else
  dimValue = size(matValue, 1);
end
end

function traceValue = localSafeTrace(matValue, idx)
%LOCALSAFETRACE Return trace of a selected principal block.

traceValue = NaN;
if isempty(matValue) || ~isnumeric(matValue)
  return;
end
idx = idx(idx >= 1 & idx <= size(matValue, 1) & idx <= size(matValue, 2));
if isempty(idx)
  return;
end
traceValue = trace(matValue(idx, idx));
end

function elemValue = localSafeMatrixElem(matValue, rowIdx, colIdx)
%LOCALSAFEMATRIXELEM Return one matrix element or NaN.

elemValue = NaN;
if isnumeric(matValue) && ~isempty(matValue) && ...
    size(matValue, 1) >= rowIdx && size(matValue, 2) >= colIdx
  elemValue = matValue(rowIdx, colIdx);
end
end

function normValue = localSafeFroNorm(matValue)
%LOCALSAFEFRONORM Return Frobenius norm or NaN.

normValue = NaN;
if isnumeric(matValue) && ~isempty(matValue)
  normValue = norm(matValue, 'fro');
end
end

function condValue = localSafeCond(matValue)
%LOCALSAFECOND Return condition number or NaN.

condValue = NaN;
if isnumeric(matValue) && ~isempty(matValue) && all(isfinite(matValue(:)))
  condValue = cond(matValue);
end
end

function [eigValue, eigVector] = localSortedEigenSystem(matValue)
%LOCALSORTEDEIGENSYSTEM Return ascending eigenvalues and aligned eigenvectors.

eigValue = [];
eigVector = [];
if isempty(matValue) || ~isnumeric(matValue) || size(matValue, 1) ~= size(matValue, 2) || ...
    any(~isfinite(matValue(:)))
  return;
end
matValue = (matValue + matValue') / 2;
[eigVector, eigDiag] = eig(matValue);
eigValue = real(diag(eigDiag));
[eigValue, sortIdx] = sort(eigValue, 'ascend');
eigVector = eigVector(:, sortIdx);
end

function invMat = localSafeInv(matValue)
%LOCALSAFEINV Return a symmetric pseudo-inverse for finite square matrices.

invMat = [];
if isempty(matValue) || ~isnumeric(matValue) || size(matValue, 1) ~= size(matValue, 2) || ...
    any(~isfinite(matValue(:)))
  return;
end
invMat = pinv((matValue + matValue') / 2);
invMat = (invMat + invMat') / 2;
end

function stdValue = localScalarInfoStd(infoValue)
%LOCALSCALARINFOSTD Convert positive scalar information to standard deviation.

stdValue = NaN;
if isfinite(infoValue) && infoValue > 0
  stdValue = sqrt(1 / infoValue);
end
end

function infoValue = localScalarInfoFromStd(stdValue)
%LOCALSCALARINFOFROMSTD Convert a positive standard deviation to information.

infoValue = NaN;
if isfinite(stdValue) && stdValue > 0
  infoValue = 1 / (stdValue ^ 2);
end
end

function couplingValue = localNormalizedFimCoupling(fimMat, rowIdx, colIdx)
%LOCALNORMALIZEDFIMCOUPLING Return a dimensionless FIM coupling coefficient.

couplingValue = NaN;
if isempty(fimMat) || ~isnumeric(fimMat) || ...
    size(fimMat, 1) < max(rowIdx, colIdx) || size(fimMat, 2) < max(rowIdx, colIdx)
  return;
end
diagProd = fimMat(rowIdx, rowIdx) * fimMat(colIdx, colIdx);
if isfinite(diagProd) && diagProd > 0
  couplingValue = fimMat(rowIdx, colIdx) / sqrt(diagProd);
end
end

function ratioValue = localSafeRatio(numerator, denominator)
%LOCALSAFERATIO Return numerator / denominator when both are finite.

ratioValue = NaN;
if isfinite(numerator) && isfinite(denominator) && denominator ~= 0
  ratioValue = numerator / denominator;
end
end

function lossValue = localCouplingLossFraction(fullStd, scalarStd)
%LOCALCOUPLINGLOSSFRACTION Estimate fdRef information lost to DoA coupling.

lossValue = NaN;
if isfinite(fullStd) && isfinite(scalarStd) && fullStd > 0 && scalarStd > 0
  lossValue = 1 - (scalarStd / fullStd) ^ 2;
end
end

function stdValue = localCrbStd(crbMat, idx)
%LOCALCRBSTD Return CRB standard deviation for one parameter.

stdValue = NaN;
if isnumeric(crbMat) && size(crbMat, 1) >= idx && size(crbMat, 2) >= idx && crbMat(idx, idx) >= 0
  stdValue = sqrt(crbMat(idx, idx));
end
end

function angleStdDeg = localCrbAngleTraceStdDeg(crbMat)
%LOCALCRBANGLETRACESTDDEG Return trace-based angle standard deviation in degrees.

angleStdDeg = NaN;
if isnumeric(crbMat) && size(crbMat, 1) >= 2 && all(diag(crbMat(1:2, 1:2)) >= 0)
  angleStdDeg = sqrt(trace(crbMat(1:2, 1:2)));
  angleStdDeg = localMaybeRadToDeg(angleStdDeg);
end
end

function valueOut = localMaybeRadToDeg(valueIn)
%LOCALMAYBERADTODEG Convert small angular standard deviations when already in radians.
% Lat/lon CRB in this repository is often expressed in degrees by the
% latlon Jacobian path. Keep the numeric value unchanged for audit display.

valueOut = valueIn;
end

function status = localCaseStatus(caseInfo)
%LOCALCASESTATUS Return one compact CRB case status.

if strlength(caseInfo.errorMessage) > 0
  status = "error";
elseif isempty(caseInfo.crb) || any(~isfinite(caseInfo.crb(:)))
  status = "nan-crb";
else
  status = "ok";
end
end

function colNorm = localColumnNorm(matValue, colIdx)
%LOCALCOLUMNNORM Return norm of one matrix column or NaN.

colNorm = NaN;
if isnumeric(matValue) && ~isempty(matValue) && size(matValue, 2) >= colIdx
  colNorm = norm(matValue(:, colIdx));
end
end

function relValue = localRelDiff(leftValue, rightValue)
%LOCALRELDIFF Return scalar relative difference with a stable denominator.

if ~(isfinite(leftValue) && isfinite(rightValue))
  relValue = NaN;
  return;
end
denom = max([abs(leftValue), abs(rightValue), eps]);
relValue = abs(leftValue - rightValue) / denom;
end

function relValue = localRelMatrixDiff(leftMat, rightMat)
%LOCALRELMATRIXDIFF Return relative Frobenius difference for two matrices.

if isempty(leftMat) || isempty(rightMat) || ~isnumeric(leftMat) || ~isnumeric(rightMat) || ...
    ~isequal(size(leftMat), size(rightMat))
  relValue = NaN;
  return;
end
denom = max([norm(leftMat, 'fro'), norm(rightMat, 'fro'), eps]);
relValue = norm(leftMat - rightMat, 'fro') / denom;
end

function previewTable = localPreviewRows(inputTable, maxRow)
%LOCALPREVIEWROWS Return the first rows of a table for compact display.

if height(inputTable) <= maxRow
  previewTable = inputTable;
else
  previewTable = inputTable(1:maxRow, :);
end
end

function maxValue = localMaxFinite(valueArray)
%LOCALMAXFINITE Return max finite value or NaN.

valueArray = valueArray(isfinite(valueArray));
if isempty(valueArray)
  maxValue = NaN;
else
  maxValue = max(valueArray);
end
end

