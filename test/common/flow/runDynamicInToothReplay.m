function [caseDynMsUnknown, runMeta] = runDynamicInToothReplay(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeInTooth, fdRateRangeInTooth, selectedSubsetCase, ...
  selectedSubsetSummary, anchorUnknownSummary, knownSummary)
%RUNDYNAMICINTOOTHREPLAY Run one narrow same-tooth periodic replay.
% This helper keeps the current fixed-tooth replay behavior but moves the
% recovered-tooth DoA seeding and same-tooth local-box resolution out of the
% main transition bundle file.

arguments
  periodicFixture (1, 1) struct
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical
  flowOpt (1, 1) struct
  fdRangeInTooth (1, 2) double
  fdRateRangeInTooth (1, 2) double
  selectedSubsetCase (1, 1) struct
  selectedSubsetSummary (1, 1) struct
  anchorUnknownSummary (1, 1) struct
  knownSummary (1, 1) struct
end

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
[inToothDoaInitParam, inToothDoaSeedSource] = resolveRecoveredToothDoaSeed( ...
  flowOpt, selectedSubsetCase, selectedSubsetSummary, knownSummary);
dynMsUnknownOpt.initDoaParam = inToothDoaInitParam(:);
[inToothDoaHalfWidthDeg, inToothOptimOpt] = localResolveInToothDoaConfig( ...
  flowOpt, selectedSubsetSummary, anchorUnknownSummary, knownSummary);
dynMsUnknownOpt.initDoaHalfWidth = inToothDoaHalfWidthDeg;

runMeta = struct();
runMeta.doaInitParam = reshape(inToothDoaInitParam(:), 1, []);
runMeta.doaSeedSource = string(inToothDoaSeedSource);
runMeta.doaHalfWidthDeg = reshape(inToothDoaHalfWidthDeg(:), 1, []);
runMeta.optimOpt = inToothOptimOpt;

dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.freezeDoa = logical(localGetFieldOrDefault(flowOpt, 'inToothFreezeDoa', true));
dynMsUnknownOpt.disableUnknownDoaReleaseFloor = logical(localGetFieldOrDefault(flowOpt, 'inToothDisableUnknownDoaReleaseFloor', true));
dynMsUnknownOpt.unknownDoaReleaseHalfWidth = inToothDoaHalfWidthDeg(:);
dynMsUnknownOpt.continuousPhaseConsistencyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseConsistencyWeight', ...
  localGetFieldOrDefault(dynMsUnknownOpt, 'continuousPhaseConsistencyWeight', 0));
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);
dynMsUnknownOpt.optimOpt = localMergeOptimOpt(localGetFieldOrDefault(dynMsUnknownOpt, 'optimOpt', struct()), ...
  inToothOptimOpt);

initParamMsUnknown = buildDynamicInitParamFromCase(selectedSubsetCase, false, selectedSubsetCase.estResult.fdRateEst);
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U-inTooth", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeInTooth, fdRateRangeInTooth, optVerbose, dynMsUnknownOpt, false, ...
  periodicFixture.debugTruthMs, initParamMsUnknown);
end

function [doaHalfWidthDeg, optimOpt] = localResolveInToothDoaConfig(flowOpt, selectedSubsetSummary, anchorUnknownSummary, knownSummary)
%LOCALRESOLVEINTOOTHDOACONFIG Choose one in-tooth DoA box and local solver budget.
% Keep the default tiny in-tooth replay for normal cases. Only when the
% selected tooth comes from a random rescue and the anchor result still
% drifts noticeably from the known-rate basin do we reopen a larger local
% DoA box.

doaHalfWidthDeg = reshape(localGetFieldOrDefault(flowOpt, 'inToothDoaHalfWidthDeg', [2e-4; 2e-4]), [], 1);
optimOpt = localGetFieldOrDefault(flowOpt, 'inToothOptimOpt', struct());
if ~logical(localGetFieldOrDefault(flowOpt, 'enableAdaptiveInToothDoaEscalation', false))
  return;
end
subsetLabel = lower(strtrim(char(string(localGetFieldOrDefault(selectedSubsetSummary, 'subsetLabel', "")))));
selectedToothIdx = localGetFieldOrDefault(selectedSubsetSummary, 'toothIdx', NaN);
if ~(isfinite(selectedToothIdx) && abs(selectedToothIdx) == 0)
  return;
end
selectedResidualHz = abs(localGetFieldOrDefault(selectedSubsetSummary, 'toothResidualHz', inf));
maxSelectedResidualHz = localGetFieldOrDefault(flowOpt, 'inToothEscalateMaxSelectedResidualHz', inf);
if ~(isfinite(selectedResidualHz) && selectedResidualHz <= maxSelectedResidualHz)
  return;
end
randomSelected = startsWith(subsetLabel, 'random') || startsWith(subsetLabel, 'rescue-random');
if randomSelected && ~logical(localGetFieldOrDefault(flowOpt, 'inToothEscalateWhenRandomSelected', true))
  return;
end
anchorDoaDriftFromKnownDeg = localCalcSummaryDoaDriftDeg(anchorUnknownSummary, knownSummary);
minAnchorDoaDriftFromKnownDeg = localGetFieldOrDefault(flowOpt, 'inToothEscalateMinAnchorDoaDriftFromKnownDeg', inf);
anchorHealthBucket = localBuildSubsetHealthBucket(anchorUnknownSummary);
maxAnchorHealthBucket = localGetFieldOrDefault(flowOpt, 'inToothEscalateMaxAnchorHealthBucket', 0);
needsEscalation = false;
if isfinite(anchorDoaDriftFromKnownDeg) && anchorDoaDriftFromKnownDeg >= minAnchorDoaDriftFromKnownDeg
  needsEscalation = true;
end
if anchorHealthBucket > maxAnchorHealthBucket
  needsEscalation = true;
end
if needsEscalation
  doaHalfWidthDeg = reshape(localGetFieldOrDefault(flowOpt, 'inToothEscalatedDoaHalfWidthDeg', doaHalfWidthDeg), [], 1);
  optimOpt = localMergeOptimOpt(optimOpt, localGetFieldOrDefault(flowOpt, 'inToothEscalatedOptimOpt', struct()));
end
end


function bucket = localBuildSubsetHealthBucket(summary)
%LOCALBUILDSUBSETHEALTHBUCKET Count coarse non-reference quality failures.

bucket = 0;
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefSupportRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefFitRatioFloor', 1) < 0.95);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefConsistencyRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefCoherenceFloor', 1) < 0.995);
rmsPhaseResid = localGetFieldOrDefault(summary, 'nonRefRmsPhaseResidRad', NaN);
if isfinite(rmsPhaseResid)
  bucket = bucket + double(rmsPhaseResid > 0.003);
end
maxPhaseResid = localGetFieldOrDefault(summary, 'nonRefMaxAbsPhaseResidRad', NaN);
if isfinite(maxPhaseResid)
  bucket = bucket + double(maxPhaseResid > 0.005);
end
negativeRatio = localGetFieldOrDefault(summary, 'maxNonRefNegativeProjectionRatio', NaN);
if isfinite(negativeRatio)
  bucket = bucket + double(negativeRatio > 0.05);
end
end


function doaDriftDeg = localCalcSummaryDoaDriftDeg(summaryA, summaryB)
%LOCALCALCSUMMARYDOADRIFTDEG Build one DoA drift metric between compact summaries.

doaDriftDeg = NaN;
doaA = localGetFieldOrDefault(summaryA, 'doaParamEst', []);
doaB = localGetFieldOrDefault(summaryB, 'doaParamEst', []);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = NaN;
end
end


function optimOpt = localMergeOptimOpt(baseOptimOpt, overrideOptimOpt)
%LOCALMERGEOPTIMOPT Merge one local optimizer override struct.

optimOpt = baseOptimOpt;
if nargin < 1 || isempty(optimOpt)
  optimOpt = struct();
end
if nargin < 2 || isempty(overrideOptimOpt)
  return;
end
if ~isstruct(overrideOptimOpt)
  error('runDynamicInToothReplay:InvalidLocalOptimOpt', ...
    'Optimizer overrides must be provided as a struct.');
end
fieldList = fieldnames(overrideOptimOpt);
for iField = 1:numel(fieldList)
  optimOpt.(fieldList{iField}) = overrideOptimOpt.(fieldList{iField});
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or a default value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
