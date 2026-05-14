% replayMfCrbTexConsistencyAudit
% Purpose: audit CRB model dispatch, parameter partition, FIM / EFIM blocks,
% and key case-pair differences against the paper-level SF-DoA and MF
% DoA-Doppler model hierarchy.
%
% Usage: edit the configuration block below, then run this script directly.
% The replay does not run MLE. It builds one fixed dynamic context, evaluates
% CRB helpers for the selected model-matched cases, and stores lightweight
% audit tables in replayData.

clear; close all; clc;

%% Replay configuration

replayName = "replayMfCrbTexConsistencyAudit";
saveSnapshot = true;
notifyTelegramEnable = false;

% This first audit is model-matched: SF cases use the DoA-only effective CRB;
% MF cases use their corresponding zero-rate / known-rate / unknown-rate CRB.
auditMode = "model-matched";
numFrame = 10;
snrDb = 0;
baseSeed = 253;
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
  config.methodNameList = methodNameList;
  config.identicalRelTol = double(identicalRelTol);
  config.summaryReuseRelTol = double(summaryReuseRelTol);
  config.saveSnapshot = saveSnapshot;
  config.notifyTelegramEnable = notifyTelegramEnable;
  config.runKey = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));

  printMfReplayHeader(char(replayName), config, '');
  localPrintReplayConfig(config);

  contextOpt = struct();
  contextOpt.periodicOffsetIdx = localCenteredOffsets(config.numFrame);
  contextOpt.masterOffsetIdx = contextOpt.periodicOffsetIdx;
  contextOpt.numSubsetRandomTrial = 0;
  contextOpt.parallelOpt = struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', inf);
  context = buildDynamicDualSatEciContext(contextOpt);
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
  printMfReplaySection('CRB FIM / EFIM block audit', replayData.crbFimBlockAuditTable);
  printMfReplaySection('CRB fdRef marginal audit', replayData.crbFdRefMarginalAuditTable);
  printMfReplaySection('CRB eigen / coupling audit', replayData.crbEigenCouplingAuditTable);
  printMfReplaySection('CRB case-pair audit', replayData.crbCasePairAuditTable);
  printMfReplaySection('CRB summary-reuse suspects', replayData.crbSummaryReuseAuditTable);
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
  crbCaseList(iCase) = localBuildSingleCrbCase(methodName, periodicFixture, context);
end
end

function caseInfo = localBuildSingleCrbCase(methodName, periodicFixture, context)
%LOCALBUILDSINGLECRBCASE Build one CRB case and record its expected model tags.

caseInfo = localEmptyCrbCase();
caseInfo.caseName = string(methodName);
caseInfo.expected = localExpectedModelSpec(methodName);
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
  string(mat2str(context.selectedSatIdxGlobal)), ...
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
    'tleFileName', 'selectedSatIdxGlobal', 'refSatIdxGlobal', 'refSatIdxLocal', ...
    'numSat', 'numFrame', 'timeOffsetSec', 'fdRefTrueHz', 'fdRefFitHz', ...
    'fdRateFitHzPerSec', 'noiseVar'});
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

function expected = localExpectedModelSpec(methodName)
%LOCALEXPECTEDMODELSPEC Return the paper-level expected model tags.

methodName = string(methodName);
expected = struct();
expected.caseName = methodName;
expected.numSat = 1 + startsWith(methodName, "MS-");
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
  'crbDerivativeNormTable', 'crbCasePairAuditTable'};
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

