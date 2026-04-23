function bundle = buildDoaDopplerDynamicTransitionBundle(periodicFixture, subsetFixtureCell, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt)
%BUILDDOADOPPLERDYNAMICTRANSITIONBUNDLE Build one shared static->dynamic MF bundle.
% This helper keeps the dynamic dev/perf scripts aligned with the static
% dual-satellite scripts: the entry script only prepares the scene/snapshot
% fixture and this helper runs the shared SF transition plus the current MF
% CP-K / CP-U flow.

arguments
  periodicFixture (1, 1) struct
  subsetFixtureCell cell
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical = false
  flowOpt (1, 1) struct = struct()
end

flowOpt = localApplyDynamicFlowDefaults(flowOpt);
bundleTimer = tic;
stageTiming = struct();
truth = periodicFixture.truth;
viewRefOnly = periodicFixture.viewRefOnly;
viewOtherOnly = periodicFixture.viewOtherOnly;
viewMs = periodicFixture.viewMs;

stageTimer = tic;
caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, periodicFixture.wavelen, pilotWave, ...
  carrierFreq, sampleRate, periodicFixture.fdRange, truth, ...
  periodicFixture.otherSatIdxGlobal, optVerbose, flowOpt.doaOnlyOpt, ...
  flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));
stageTiming.staticTransitionSec = toc(stageTimer);

caseRefDoa = caseBundle.caseRefDoa;
caseMsDoa = caseBundle.caseMsDoa;
caseStaticRefOnly = caseBundle.caseStaticRefOnly;
caseStaticRefAbl = caseBundle.caseStaticRefAbl;
caseStaticOtherOnly = caseBundle.caseStaticOtherOnly;
caseStaticMs = caseBundle.caseStaticMs;
weightCase = caseBundle.weightCase;
bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

initParamStaticRef = buildDynamicInitParamFromCase(caseStaticRefOnly, true, truth.fdRateFit);

periodicFixture.debugTruthRef = buildTruthDebugForView(viewRefOnly, truth);
periodicFixture.debugTruthMs = buildTruthDebugForView(viewMs, truth);
stageTimer = tic;
caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef);
stageTiming.refKnownSec = toc(stageTimer);
stageTimer = tic;
caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly);
stageTiming.refUnknownSec = toc(stageTimer);

stageTimer = tic;
caseDynMsKnown = localRunMsKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, bestStaticMsCase);
stageTiming.msKnownSec = toc(stageTimer);
knownSummaryForRecovery = localBuildUnknownCaseSummary(caseDynMsKnown, truth, ...
  1 / median(diff(periodicFixture.sceneSeq.timeOffsetSec)));

toothStepHz = 1 / median(diff(periodicFixture.sceneSeq.timeOffsetSec));
runWideBranches = logical(localGetFieldOrDefault(flowOpt, 'enableUnknownWideBranches', true));
runWideRefine = runWideBranches && logical(localGetFieldOrDefault(flowOpt, 'enableUnknownWideRefine', true));
deferWideBranches = runWideBranches && logical(localGetFieldOrDefault(flowOpt, 'deferUnknownWideBranches', false));
enableSubsetInToothRefine = logical(localGetFieldOrDefault(flowOpt, 'enableSubsetInToothRefine', true));
subsetSeedInfo = localBuildSubsetSeedReuseInfo(bestStaticMsCase, caseDynMsKnown, staticMsOpt, flowOpt);

caseDynMsUnknownFdWide = struct();
caseDynMsUnknownWide = struct();
fdWideUnknownSummary = localBuildSkippedUnknownSummary("periodic-fd-wide");
wideUnknownSummary = localBuildSkippedUnknownSummary("periodic-wide");
if runWideBranches && ~deferWideBranches
  stageTimer = tic;
  [caseDynMsUnknownFdWide, caseDynMsUnknownWide, fdWideUnknownSummary, wideUnknownSummary] = ...
    runDynamicUnknownWideRoutes(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
    optVerbose, flowOpt, bestStaticMsCase, truth, toothStepHz, runWideRefine, true, ...
    @localBuildUnknownCaseSummary);
  stageTiming.msUnknownFdWideSec = toc(stageTimer);
  stageTiming.msUnknownWideSec = localGetFieldOrDefault(localGetFieldOrDefault(caseDynMsUnknownWide, 'estResult', struct()), 'runTimeMs', NaN) / 1e3;
else
  stageTiming.msUnknownFdWideSec = 0;
  stageTiming.msUnknownWideSec = 0;
end

stageTimer = tic;
[subsetCaseCell, subsetSummaryCell, subsetRunTimeSec] = evaluateDynamicSubsetBank( ...
  subsetFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, subsetSeedInfo, truth, toothStepHz, ...
  @localRunMsUnknownCaseFromFixture, @localBuildUnknownCaseSummary);

[subsetFixtureCell, subsetCaseCell, subsetSummaryCell, subsetRunTimeSec, ...
  bestSubsetIdx, subsetIsTrusted, subsetSelectReason, selectedSubsetFixture, ...
  selectedSubsetSummary, selectedSubsetCase, subsetCandidateTable, ...
  subsetEscalationReason, subsetEscalationDiag, conditionalRandomReason, conditionalRandomDiag] = ...
  recoverDynamicSubsetBank(periodicFixture, subsetFixtureCell, subsetCaseCell, subsetSummaryCell, ...
  subsetRunTimeSec, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, ...
  subsetSeedInfo, truth, toothStepHz, runWideBranches, caseDynMsUnknownFdWide, ...
  fdWideUnknownSummary, caseDynMsKnown, knownSummaryForRecovery, ...
  @localRunMsUnknownCaseFromFixture, @localBuildUnknownCaseSummary);
stageTiming.subsetSelectTotalSec = toc(stageTimer);
stageTiming.subsetRankingSec = 0;

needWideFallback = false;
wideFallbackReason = "not-needed";
if runWideBranches && deferWideBranches
  wideResidualTolHz = localGetFieldOrDefault(flowOpt, 'fastWideToothResidualHzThreshold', 50);
  allowWideFallback = logical(localGetFieldOrDefault(flowOpt, 'enableFastWideFallback', true));
  if allowWideFallback && localHasResolvedCentralTooth(subsetSummaryCell, wideResidualTolHz)
    allowWideFallback = false;
    wideFallbackReason = "central-tooth subset already available";
  end
  wideFallbackOpt = struct( ...
    'enableEscalation', allowWideFallback, ...
    'escalateOnToothDisagreement', false, ...
    'maxTrustedToothIdx', localGetFieldOrDefault(flowOpt, 'fastWideToothIdxThreshold', 0), ...
    'maxTrustedToothResidualHz', wideResidualTolHz, ...
    'maxFdDriftFromKnownHz', localGetFieldOrDefault(flowOpt, 'fastWideFdDriftHzThreshold', 500), ...
    'maxDoaDriftFromKnownDeg', localGetFieldOrDefault(flowOpt, 'fastWideDoaDriftDegThreshold', 0.003));
  if allowWideFallback
    [needWideFallback, wideFallbackReason] = shouldEscalateDynamicSubsetRecovery( ...
      subsetSummaryCell, selectedSubsetSummary, knownSummaryForRecovery, wideFallbackOpt);
  else
    needWideFallback = false;
  end
  if needWideFallback
    stageTimer = tic;
    [caseDynMsUnknownFdWide, caseDynMsUnknownWide, fdWideUnknownSummary, wideUnknownSummary] = ...
      runDynamicUnknownWideRoutes(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
      optVerbose, flowOpt, selectedSubsetCase, truth, toothStepHz, runWideRefine, ...
      logical(localGetFieldOrDefault(flowOpt, 'enableFastFdWideFallback', false)), ...
      @localBuildUnknownCaseSummary);
    if localGetFieldOrDefault(flowOpt, 'enableFastFdWideFallback', false)
      stageTiming.msUnknownFdWideSec = toc(stageTimer);
      stageTiming.msUnknownWideSec = localGetFieldOrDefault(localGetFieldOrDefault(caseDynMsUnknownWide, 'estResult', struct()), 'runTimeMs', NaN) / 1e3;
    else
      stageTiming.msUnknownFdWideSec = 0;
      stageTiming.msUnknownWideSec = toc(stageTimer);
    end
  end
end

if ~localCaseHasUsableEstimate(selectedSubsetCase)
  if runWideBranches && localCaseHasUsableEstimate(caseDynMsUnknownFdWide)
    selectedSubsetCase = caseDynMsUnknownFdWide;
    selectedSubsetSummary = fdWideUnknownSummary;
    selectedSubsetFixture = struct('subsetLabel', "wide-fallback", 'subsetOffsetIdx', []);
  else
    selectedSubsetCase = caseDynMsKnown;
    selectedSubsetSummary = knownSummaryForRecovery;
    selectedSubsetFixture = struct('subsetLabel', "known-fallback", 'subsetOffsetIdx', []);
  end
end

fdRangeInTooth = [selectedSubsetCase.estResult.fdRefEst - flowOpt.inToothFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.inToothFdHalfWidthHz];
fdRateRangeInTooth = [selectedSubsetCase.estResult.fdRateEst - flowOpt.inToothFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.inToothFdRateHalfWidthHzPerSec];

[caseDynMsUnknownAnchor, anchorUnknownSummary, caseDynMsUnknownAnchorDoaPolish, ...
  anchorDoaPolishSummary, fdRangeSubsetAnchor, fdRateRangeSubsetAnchor, anchorRouteMeta] = ...
  runDynamicUnknownAnchorRoutes(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, selectedSubsetCase, selectedSubsetSummary, knownSummaryForRecovery, ...
  truth, toothStepHz, @localBuildUnknownCaseSummary);
stageTiming.subsetAnchorSec = anchorRouteMeta.subsetAnchorSec;
stageTiming.anchorDoaPolishSec = anchorRouteMeta.anchorDoaPolishSec;
hasAnchorDoaPolish = anchorRouteMeta.hasAnchorDoaPolish;

stageTimer = tic;
didRunInToothRefine = false;
inToothRunMeta = struct();
if enableSubsetInToothRefine && subsetIsTrusted && ...
    localShouldRunInToothRefine(selectedSubsetSummary, anchorUnknownSummary, knownSummaryForRecovery, flowOpt)
  [refinedCaseDynMsUnknown, inToothRunMeta] = runDynamicInToothReplay(periodicFixture, pilotWave, ...
    carrierFreq, sampleRate, optVerbose, flowOpt, ...
    fdRangeInTooth, fdRateRangeInTooth, selectedSubsetCase, ...
    selectedSubsetSummary, anchorUnknownSummary, knownSummaryForRecovery);
  refinedUnknownSummary = localBuildUnknownCaseSummary(refinedCaseDynMsUnknown, truth, toothStepHz);
  didRunInToothRefine = true;
else
  refinedCaseDynMsUnknown = caseDynMsUnknownAnchor;
  refinedUnknownSummary = anchorUnknownSummary;
end
stageTiming.inToothRefineSec = toc(stageTimer);

stageTimer = tic;
selectedSummaryCell = {};
selectedTagList = strings(0, 1);
caseFinalCell = {};
if runWideBranches && localCaseHasUsableEstimate(caseDynMsUnknownFdWide)
  selectedSummaryCell = [selectedSummaryCell; {wideUnknownSummary; fdWideUnknownSummary}];
  selectedTagList = [selectedTagList; "periodic-wide"; "periodic-fd-wide"];
  caseFinalCell = [caseFinalCell; {caseDynMsUnknownWide; caseDynMsUnknownFdWide}];
  subsetAnchorIdx = 3;
  enableSubsetAnchorGuard = subsetIsTrusted;
else
  subsetAnchorIdx = 1;
  enableSubsetAnchorGuard = false;
end
selectedSummaryCell = [selectedSummaryCell; {anchorUnknownSummary}];
selectedTagList = [selectedTagList; "periodic-fd-anchor"];
caseFinalCell = [caseFinalCell; {caseDynMsUnknownAnchor}];
if hasAnchorDoaPolish
  selectedSummaryCell = [selectedSummaryCell; {anchorDoaPolishSummary}];
  selectedTagList = [selectedTagList; "periodic-fd-anchor-doa-polish"];
  caseFinalCell = [caseFinalCell; {caseDynMsUnknownAnchorDoaPolish}];
end
if didRunInToothRefine
  selectedSummaryCell = [selectedSummaryCell; {refinedUnknownSummary}];
  selectedTagList = [selectedTagList; "periodic-in-tooth"];
  caseFinalCell = [caseFinalCell; {refinedCaseDynMsUnknown}];
end
finalSelectOpt = struct( ...
  'enableSubsetAnchorGuard', enableSubsetAnchorGuard, ...
  'subsetAnchorIdx', subsetAnchorIdx, ...
  'subsetAnchorTag', "periodic-fd-anchor", ...
  'maxFdRefDriftHz', flowOpt.finalSelectMaxFdRefDriftHz, ...
  'maxFdRateDriftHzPerSec', flowOpt.finalSelectMaxFdRateDriftHzPerSec, ...
  'maxDoaDriftDeg', flowOpt.finalSelectMaxDoaDriftDeg, ...
  'minObjGainToLeaveAnchor', flowOpt.finalSelectMinObjGainToLeaveAnchor, ...
  'minResidualGainToLeaveAnchor', flowOpt.finalSelectMinResidualGainToLeaveAnchor, ...
  'preferReleasedSameTooth', localGetFieldOrDefault(flowOpt, 'msFinalSelectPreferReleasedSameTooth', false), ...
  'sameToothObjectiveToleranceAbs', localGetFieldOrDefault(flowOpt, 'msFinalSelectSameToothObjectiveToleranceAbs', 0), ...
  'sameToothResidualToleranceAbs', localGetFieldOrDefault(flowOpt, 'msFinalSelectSameToothResidualToleranceAbs', 0), ...
  'sameToothObjectiveToleranceRel', localGetFieldOrDefault(flowOpt, 'msFinalSelectSameToothObjectiveToleranceRel', 0), ...
  'sameToothResidualToleranceRel', localGetFieldOrDefault(flowOpt, 'msFinalSelectSameToothResidualToleranceRel', 0));
[selectedFinalIdx, selectedFinalTag, finalSelectReason, finalScoreMat] = selectFinalDynamicUnknownSummary( ...
  selectedSummaryCell, selectedTagList, finalSelectOpt);
caseDynMsUnknown = caseFinalCell{selectedFinalIdx};
stageTiming.finalSelectSec = toc(stageTimer);

baseCase = [caseRefDoa, caseMsDoa, caseStaticRefOnly, caseStaticMs, ...
  caseDynRefKnown, caseDynRefUnknown, caseDynMsKnown, caseDynMsUnknown];
caseResult = [baseCase, weightCase];

stageTimer = tic;
dynSummary = summarizeDynamicEstimatorDebug(baseCase, periodicFixture.sceneSeq, truth);
dynMultiStartTable = summarizeDynamicMultiStart(baseCase, truth);
stageTiming.summaryBuildSec = toc(stageTimer);
stageTiming.totalBundleSec = toc(bundleTimer);

tableBundle = buildDynamicTransitionTables(stageTiming, subsetFixtureCell, subsetSummaryCell, ...
  subsetRunTimeSec, selectedSummaryCell, selectedTagList, string(selectedFinalTag), finalScoreMat, ...
  fdWideUnknownSummary, wideUnknownSummary, selectedSubsetSummary, anchorUnknownSummary, ...
  anchorDoaPolishSummary, refinedUnknownSummary, hasAnchorDoaPolish, didRunInToothRefine);
toothFlowTable = tableBundle.toothFlowTable;
finalSelectTable = tableBundle.finalSelectTable;

bundle = struct();
bundle.caseBundle = caseBundle;
bundle.caseRefDoa = caseRefDoa;
bundle.caseMsDoa = caseMsDoa;
bundle.caseStaticRefOnly = caseStaticRefOnly;
bundle.caseStaticRefAbl = caseStaticRefAbl;
bundle.caseStaticOtherOnly = caseStaticOtherOnly;
bundle.caseStaticMs = caseStaticMs;
bundle.weightCase = weightCase;
bundle.bestStaticMsCase = bestStaticMsCase;
bundle.caseDynRefKnown = caseDynRefKnown;
bundle.caseDynRefUnknown = caseDynRefUnknown;
bundle.caseDynMsKnown = caseDynMsKnown;
bundle.caseDynMsUnknownFdWide = caseDynMsUnknownFdWide;
bundle.caseDynMsUnknownWide = caseDynMsUnknownWide;
bundle.caseDynMsUnknownAnchor = caseDynMsUnknownAnchor;
bundle.caseDynMsUnknownAnchorDoaPolish = caseDynMsUnknownAnchorDoaPolish;
% Keep one explicit alias for the final selected CP-U case so repeat-level
% statistics and regressions can consume the post-final-select winner
% without having to infer it from tables or candidate tags.
bundle.caseDynMsUnknownFinal = caseDynMsUnknown;
bundle.caseDynMsUnknown = caseDynMsUnknown;
bundle.baseCase = baseCase;
bundle.caseResult = caseResult;
bundle.truth = truth;
bundle.toothStepHz = toothStepHz;
bundle.fdWideUnknownSummary = fdWideUnknownSummary;
bundle.wideUnknownSummary = wideUnknownSummary;
bundle.selectedSubsetSummary = selectedSubsetSummary;
bundle.anchorUnknownSummary = anchorUnknownSummary;
bundle.anchorDoaPolishSummary = anchorDoaPolishSummary;
bundle.refinedUnknownSummary = refinedUnknownSummary;
bundle.selectedFinalSummary = selectedSummaryCell{selectedFinalIdx};
bundle.selectedSubsetLabel = string(localGetFieldOrDefault(selectedSubsetFixture, 'subsetLabel', "none"));
bundle.selectedSubsetOffsetIdx = reshape(localGetFieldOrDefault(selectedSubsetFixture, 'subsetOffsetIdx', []), 1, []);
bundle.bestSubsetIdx = bestSubsetIdx;
bundle.subsetIsTrusted = subsetIsTrusted;
bundle.selectedFinalTag = string(selectedFinalTag);
bundle.finalSelectReason = string(finalSelectReason);
bundle.subsetSelectReason = string(subsetSelectReason);
bundle.subsetEscalationReason = string(subsetEscalationReason);
bundle.subsetEscalationDiag = subsetEscalationDiag;
bundle.conditionalRandomReason = string(conditionalRandomReason);
bundle.conditionalRandomDiag = conditionalRandomDiag;
bundle.wideFallbackReason = string(wideFallbackReason);
bundle.finalSelectTable = finalSelectTable;
bundle.fdRangeSubsetAnchor = fdRangeSubsetAnchor;
bundle.fdRateRangeSubsetAnchor = fdRateRangeSubsetAnchor;
bundle.fdRangeInTooth = fdRangeInTooth;
bundle.fdRateRangeInTooth = fdRateRangeInTooth;
bundle.inToothRunMeta = inToothRunMeta;
bundle.toothFlowTable = toothFlowTable;
bundle.subsetCandidateTable = subsetCandidateTable;
bundle.dynSummary = dynSummary;
bundle.dynMultiStartTable = dynMultiStartTable;
bundle.stageTiming = stageTiming;
bundle.timingTable = tableBundle.timingTable;
bundle.subsetTimingTable = tableBundle.subsetTimingTable;
end

function flowOpt = localApplyDynamicFlowDefaults(flowOpt)
%LOCALAPPLYDYNAMICFLOWDEFAULTS Fill one dynamic-flow option struct.

flowOpt = localMergeStruct(struct( ...
  'enableWeightSweep', true, ...
  'weightSweepAlpha', [0; 0.25; 0.5; 1], ...
  'staticMsHalfWidth', [0.002; 0.002], ...
  'doaOnlyOpt', struct('useLogObjective', true), ...
  'staticBaseOpt', struct('useLogObjective', true), ...
  'dynBaseOpt', struct( ...
    'useLogObjective', true, ...
    'initFdCount', 81, ...
    'useAccessMask', false, ...
    'phaseMode', 'continuous', ...
    'steeringMode', 'framewise', ...
    'debugEnable', true, ...
    'debugStoreEvalTrace', false, ...
    'debugMaxEvalTrace', 120), ...
  'msContinuousPhaseConsistencyWeight', 0.2, ...
  'msContinuousPhaseCollapsePenaltyWeight', 0.10, ...
  'msContinuousPhaseNegativeProjectionPenaltyWeight', 0.10, ...
  'msContinuousPhaseNonRefFitFloorWeight', 0.02, ...
  'msContinuousPhaseNonRefSupportFloorWeight', 0.02, ...
  'msFinalSelectPreferReleasedSameTooth', false, ...
  'msFinalSelectSameToothObjectiveToleranceAbs', 50, ...
  'msFinalSelectSameToothResidualToleranceAbs', 50, ...
  'msFinalSelectSameToothObjectiveToleranceRel', 1e-6, ...
  'msFinalSelectSameToothResidualToleranceRel', 1e-6, ...
  'enableUnknownWideBranches', true, ...
  'enableUnknownWideRefine', true, ...
  'deferUnknownWideBranches', false, ...
  'enableFastWideFallback', false, ...
  'fastWideToothIdxThreshold', 0, ...
  'fastWideToothResidualHzThreshold', 50, ...
  'fastWideFdDriftHzThreshold', 500, ...
  'fastWideDoaDriftDegThreshold', 0.003, ...
  'parallelOpt', struct('enableSubsetEvalParfor', true, 'minSubsetEvalParfor', 4), ...
  'enableSubsetInToothRefine', true, ...
  'enableSubsetCaseSeedReuse', true, ...
  'enableSubsetRandomEarlyStop', true, ...
  'subsetRandomEarlyStopMinEvaluated', 4, ...
  'subsetRandomEarlyStopMaxToothResidualHz', 20, ...
  'subsetRandomEarlyStopMaxHealthBucket', 0, ...
  'enableFastSubsetEscalation', false, ...
  'fastNumRandomSubsetTrialFallback', 0, ...
  'fastRandomSeedOffset', 0, ...
  'fastEscalateOnToothDisagreement', true, ...
  'fastToothIdxThreshold', 0, ...
  'fastToothResidualHzThreshold', 50, ...
  'fastUnknownFdDriftHzThreshold', 500, ...
  'fastUnknownDoaDriftDegThreshold', 0.003, ...
  'enableConditionalRandomSubsetRescue', false, ...
  'conditionalRandomSubsetTrialCount', 1, ...
  'conditionalRandomToothResidualHzThreshold', 20, ...
  'conditionalRandomMaxHealthBucket', 1, ...
  'conditionalRandomMaxFdDriftHz', 250, ...
  'conditionalRandomMaxDoaDriftDeg', 0.003, ...
  'conditionalRandomRescueOnToothDisagreement', true, ...
  'conditionalRandomMinResolvedOffToothCount', 2, ...
  'subsetRandomRescuePromoteResidualHz', 20, ...
  'refKnownDoaHalfWidth', [0.005; 0.005], ...
  'refUnknownDoaHalfWidth', [0.003; 0.003], ...
  'msKnownDoaHalfWidth', [0.003; 0.003], ...
  'msUnknownDoaHalfWidth', [0.002; 0.002], ...
  'subsetSelectDoaHalfWidthDeg', [0.01; 0.01], ...
  'subsetAnchorFdHalfWidthHz', 50, ...
  'subsetAnchorFdRateHalfWidthHzPerSec', 100, ...
  'subsetAnchorDoaHalfWidthDeg', [1e-8; 1e-8], ...
  'subsetAnchorFreezeDoa', true, ...
  'enableAnchorDoaPolish', false, ...
  'anchorDoaPolishDoaHalfWidthDeg', [0.005; 0.005], ...
  'anchorDoaPolishFdHalfWidthHz', 20, ...
  'anchorDoaPolishFdRateHalfWidthHzPerSec', 20, ...
  'anchorDoaPolishCentralResidualTolHz', 50, ...
  'anchorDoaPolishMinNonRefFitFloorToSkip', 0.9, ...
  'anchorDoaPolishMinNonRefSupportFloorToSkip', 0.8, ...
  'anchorDoaPolishDisableUnknownDoaReleaseFloor', true, ...
  'inToothFreezeDoa', false, ...
  'fdWideFreezeDoa', true, ...
  'wideRefineFdHalfWidthHz', 100, ...
  'wideRefineFdRateHalfWidthHzPerSec', 200, ...
  'wideRefineDoaHalfWidthDeg', [0.003; 0.003], ...
  'wideRefineDisableUnknownDoaReleaseFloor', true, ...
  'inToothFdHalfWidthHz', 50, ...
  'inToothFdRateHalfWidthHzPerSec', 100, ...
  'inToothDoaHalfWidthDeg', [2e-4; 2e-4], ...
  'inToothDisableUnknownDoaReleaseFloor', true, ...
  'inToothOptimOpt', struct('MaxIterations', 120, 'MaxFunctionEvaluations', 1200), ...
  'enableAdaptiveInToothDoaEscalation', false, ...
  'inToothEscalateWhenRandomSelected', true, ...
  'inToothEscalateMinAnchorDoaDriftFromKnownDeg', 0.006, ...
  'inToothEscalateMaxAnchorHealthBucket', 0, ...
  'inToothEscalateMaxSelectedResidualHz', 20, ...
  'inToothEscalatedDoaHalfWidthDeg', [0.01; 0.01], ...
  'inToothEscalatedOptimOpt', struct('MaxIterations', 180, 'MaxFunctionEvaluations', 1800), ...
  'enableInToothCentralSkip', true, ...
  'inToothCentralResidualTolHz', 5, ...
  'enableInToothHealthyAnchorSkip', true, ...
  'inToothSkipMinNonRefFitFloor', 0.95, ...
  'inToothSkipMinNonRefSupportFloor', 0.999, ...
  'inToothSkipMinNonRefConsistencyFloor', 0.999, ...
  'inToothSkipMinNonRefCoherenceFloor', 0.999, ...
  'inToothSkipMaxRmsPhaseResidRad', 0.003, ...
  'inToothSkipMaxAbsPhaseResidRad', 0.005, ...
  'inToothSkipMaxAnchorDoaDriftFromKnownDeg', 0.005, ...
  'finalSelectMaxFdRefDriftHz', 100, ...
  'finalSelectMaxFdRateDriftHzPerSec', 100, ...
  'finalSelectMaxDoaDriftDeg', 0.01, ...
  'finalSelectMinObjGainToLeaveAnchor', 1e4, ...
  'finalSelectMinResidualGainToLeaveAnchor', 1e3), flowOpt);

flowOpt.weightSweepAlpha = reshape(flowOpt.weightSweepAlpha, [], 1);
if ~logical(flowOpt.enableWeightSweep)
  flowOpt.weightSweepAlpha = zeros(0, 1);
end
flowOpt.staticMsHalfWidth = reshape(flowOpt.staticMsHalfWidth, [], 1);
flowOpt.refKnownDoaHalfWidth = reshape(flowOpt.refKnownDoaHalfWidth, [], 1);
flowOpt.refUnknownDoaHalfWidth = reshape(flowOpt.refUnknownDoaHalfWidth, [], 1);
flowOpt.msKnownDoaHalfWidth = reshape(flowOpt.msKnownDoaHalfWidth, [], 1);
flowOpt.msUnknownDoaHalfWidth = reshape(flowOpt.msUnknownDoaHalfWidth, [], 1);
flowOpt.subsetSelectDoaHalfWidthDeg = reshape(flowOpt.subsetSelectDoaHalfWidthDeg, [], 1);
flowOpt.subsetAnchorDoaHalfWidthDeg = reshape(flowOpt.subsetAnchorDoaHalfWidthDeg, [], 1);
flowOpt.anchorDoaPolishDoaHalfWidthDeg = reshape(flowOpt.anchorDoaPolishDoaHalfWidthDeg, [], 1);
flowOpt.wideRefineDoaHalfWidthDeg = reshape(flowOpt.wideRefineDoaHalfWidthDeg, [], 1);
flowOpt.inToothDoaHalfWidthDeg = reshape(flowOpt.inToothDoaHalfWidthDeg, [], 1);
flowOpt.inToothEscalatedDoaHalfWidthDeg = reshape(flowOpt.inToothEscalatedDoaHalfWidthDeg, [], 1);
if ~isstruct(flowOpt.inToothOptimOpt)
  error('buildDoaDopplerDynamicTransitionBundle:InvalidInToothOptimOpt', ...
    'flowOpt.inToothOptimOpt must be a struct.');
end
if ~isstruct(flowOpt.inToothEscalatedOptimOpt)
  error('buildDoaDopplerDynamicTransitionBundle:InvalidInToothEscalatedOptimOpt', ...
    'flowOpt.inToothEscalatedOptimOpt must be a struct.');
end
end


function caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef)
%LOCALRUNREFKNOWNCASE Run one SS-MF-CP-K case.

dynRefKnownOpt = flowOpt.dynBaseOpt;
dynRefKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefKnownOpt.initDoaHalfWidth = flowOpt.refKnownDoaHalfWidth;
dynRefKnownOpt.initDoaParam = reshape(initParamStaticRef(1:2), [], 1);

caseDynRefKnown = runDynamicDoaDopplerCase("SS-MF-CP-K", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefKnownOpt, true, ...
  periodicFixture.debugTruthRef, initParamStaticRef);
end


function caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly)
%LOCALRUNREFUNKNOWNCASE Run one SS-MF-CP-U case.

dynRefUnknownOpt = flowOpt.dynBaseOpt;
dynRefUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefUnknownOpt.initDoaParam = caseStaticRefOnly.estResult.doaParamEst(:);
dynRefUnknownOpt.initDoaHalfWidth = flowOpt.refUnknownDoaHalfWidth;

initParamRefUnknownCpK = buildDynamicInitParamFromCase(caseDynRefKnown, false, caseDynRefKnown.estResult.fdRateEst);
initParamRefUnknownStatic = buildDynamicInitParamFromCase(caseStaticRefOnly, false, periodicFixture.truth.fdRateFit);
refUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynRefKnown, caseStaticRefOnly, initParamRefUnknownCpK, initParamRefUnknownStatic, ...
  dynRefUnknownOpt.initDoaHalfWidth, flowOpt.refKnownDoaHalfWidth);
caseDynRefUnknown = runDynamicDoaDopplerCase("SS-MF-CP-U", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefUnknownOpt, false, ...
  periodicFixture.debugTruthRef, refUnknownCand);
end


function caseDynMsKnown = localRunMsKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, bestStaticMsCase)
%LOCALRUNMSKNOWNCASE Run one MS-MF-CP-K case.

dynMsKnownOpt = flowOpt.dynBaseOpt;
dynMsKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = flowOpt.msKnownDoaHalfWidth;
dynMsKnownOpt.enableFdAliasUnwrap = true;
dynMsKnownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsKnownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsKnownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsKnownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsKnownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsKnownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsKnownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

caseDynMsKnown = runDynamicDoaDopplerCase("MS-MF-CP-K", "multi", ...
  periodicFixture.viewMs, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynMsKnownOpt, true, ...
  periodicFixture.debugTruthMs, buildDynamicInitParamFromCase(bestStaticMsCase, true, periodicFixture.truth.fdRateFit));
end

function caseDynMsUnknown = localRunMsUnknownCaseFromFixture(fixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, fdRangeUse, fdRateRangeUse, subsetSeedInfo)
%LOCALRUNMSUNKNOWNCASEFROMFIXTURE Run one MS-MF-CP-U case on one frame subset.

if nargin < 9
  subsetSeedInfo = struct();
end
reusePeriodicSeed = logical(localGetFieldOrDefault(subsetSeedInfo, 'reusePeriodicSeeds', false));
if reusePeriodicSeed
  bestStaticMsCase = localGetFieldOrDefault(subsetSeedInfo, 'bestStaticMsCase', struct());
  knownSeedCase = localGetFieldOrDefault(subsetSeedInfo, 'knownSeedCase', struct());
  staticInitDoaHalfWidth = localGetFieldOrDefault(subsetSeedInfo, 'staticInitDoaHalfWidth', [0.01; 0.01]);
else
  caseBundle = buildDoaDopplerStaticTransitionBundle( ...
    fixture.viewRefOnly, fixture.viewOtherOnly, fixture.viewMs, fixture.wavelen, ...
    pilotWave, carrierFreq, sampleRate, fdRangeUse, fixture.truth, ...
    fixture.otherSatIdxGlobal, optVerbose, flowOpt.doaOnlyOpt, ...
    flowOpt.staticBaseOpt, flowOpt.weightSweepAlpha(:), flowOpt.staticMsHalfWidth(:));
  bestStaticMsCase = caseBundle.bestStaticMsCase;
  knownSeedCase = struct();
  staticInitDoaHalfWidth = localGetFieldOrDefault(caseBundle.staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]);
end

dynMsUnknownOpt = flowOpt.dynBaseOpt;
dynMsUnknownOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = flowOpt.subsetSelectDoaHalfWidthDeg;
dynMsUnknownOpt.enableFdAliasUnwrap = true;
dynMsUnknownOpt.continuousPhaseConsistencyWeight = flowOpt.msContinuousPhaseConsistencyWeight;
dynMsUnknownOpt.continuousPhaseCollapsePenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseCollapsePenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNegativeProjectionPenaltyWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNegativeProjectionPenaltyWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefFitFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefFitFloorWeight', 0);
dynMsUnknownOpt.continuousPhaseNonRefSupportFloorWeight = localGetFieldOrDefault(flowOpt, 'msContinuousPhaseNonRefSupportFloorWeight', 0);
dynMsUnknownOpt.unknownWarmAnchorUseScaledSolve = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorUseScaledSolve', true);
dynMsUnknownOpt.unknownWarmAnchorFallbackSqp = localGetFieldOrDefault(flowOpt, 'unknownWarmAnchorFallbackSqp', true);

initParamMsUnknownStatic = buildDynamicInitParamFromCase(bestStaticMsCase, false, fixture.truth.fdRateFit);
initSeed = initParamMsUnknownStatic;
if reusePeriodicSeed && localCaseHasUsableEstimate(knownSeedCase)
  fdRateSeedKnown = localGetCaseFdRateEst(knownSeedCase, fixture.truth.fdRateFit);
  initParamMsUnknownCpK = buildDynamicInitParamFromCase(knownSeedCase, false, fdRateSeedKnown);
  initSeed = buildUnknownInitCandidateSet( ...
    knownSeedCase, bestStaticMsCase, initParamMsUnknownCpK, initParamMsUnknownStatic, ...
    dynMsUnknownOpt.initDoaHalfWidth, staticInitDoaHalfWidth);
end
caseDynMsUnknown = runDynamicDoaDopplerCase("MS-MF-CP-U", "multi", ...
  fixture.viewMs, fixture.truth, pilotWave, carrierFreq, sampleRate, ...
  fdRangeUse, fdRateRangeUse, optVerbose, dynMsUnknownOpt, false, fixture.debugTruthMs, initSeed);
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

function tf = localSummaryNeedsTrustedCentralRescue(summaryUse, residualTolHz, maxHealthBucket)
%LOCALSUMMARYNEEDSTRUSTEDCENTRALRESCUE Return true when one summary still looks unsafe.
%
% Subset schedules are tooth selectors first. Treat the current winner as
% unsafe whenever it is unresolved, still off-tooth, or only weakly central.
% Optionally gate on one coarse non-reference health bucket so random rescue
% can also replace a fragile same-tooth curated winner.

tf = true;
if ~(isstruct(summaryUse) && ~isempty(summaryUse) && ...
    logical(localGetFieldOrDefault(summaryUse, 'isResolved', false)))
  return;
end
toothIdx = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
if ~(isfinite(toothIdx) && abs(toothIdx) == 0)
  return;
end
toothResidualHz = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
if ~(isfinite(toothResidualHz) && toothResidualHz <= residualTolHz)
  return;
end
if nargin >= 3 && isfinite(maxHealthBucket)
  if localBuildSubsetHealthBucket(summaryUse) > maxHealthBucket
    return;
  end
end
tf = false;
end


function tf = localHasResolvedCentralTooth(summaryCell, residualTolHz)
%LOCALHASRESOLVEDCENTRALTOOTH Return true when one resolved subset already sits on tooth 0.
% Fast statistics should not promote wide fallback once the subset bank
% already contains one trustworthy central-tooth candidate; at that point
% the remaining work is in-tooth refinement rather than wide tooth escape.

tf = false;
for iCase = 1:numel(summaryCell)
  summaryUse = summaryCell{iCase};
  if ~logical(localGetFieldOrDefault(summaryUse, 'isResolved', false))
    continue;
  end
  toothIdx = localGetFieldOrDefault(summaryUse, 'toothIdx', NaN);
  toothResidualHz = abs(localGetFieldOrDefault(summaryUse, 'toothResidualHz', inf));
  if isfinite(toothIdx) && abs(toothIdx) == 0 && isfinite(toothResidualHz) && toothResidualHz <= residualTolHz
    tf = true;
    return;
  end
end
end


function fdDriftHz = localCalcSummaryFdDrift(summaryA, summaryB)
%LOCALCALCSUMMARYFDDRIFT Build one fdRef drift metric between compact summaries.

fdDriftHz = NaN;
fdA = localGetFieldOrDefault(summaryA, 'fdRefEst', NaN);
fdB = localGetFieldOrDefault(summaryB, 'fdRefEst', NaN);
if isfinite(fdA) && isfinite(fdB)
  fdDriftHz = abs(fdA - fdB);
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


function tf = localShouldRunInToothRefine(selectedSubsetSummary, anchorUnknownSummary, knownSummary, flowOpt)
%LOCALSHOULDRUNINTOOTHREFINE Skip redundant in-tooth replay when anchor is already healthy.
% The released in-tooth replay is useful only when the subset winner still
% needs a local tooth-centered polish. Once the subset winner is already on
% tooth 0 and the periodic fd-anchor replay also stays on tooth 0 with strong
% non-reference phase quality, rerunning a tiny DoA release rarely changes the
% final winner but can dominate fast-statistics runtime. This skip applies to
% curated and random rescue subsets alike; random rescue only selects the tooth,
% so same-tooth replay should still prove that a local release is needed.

tf = true;
if ~isfield(flowOpt, 'enableInToothCentralSkip') || ...
    ~logical(flowOpt.enableInToothCentralSkip)
  return;
end

subsetToothIdx = localGetFieldOrDefault(selectedSubsetSummary, 'toothIdx', NaN);
if ~isfinite(subsetToothIdx) || abs(subsetToothIdx) > 0
  return;
end

subsetResidualHz = abs(localGetFieldOrDefault(selectedSubsetSummary, 'toothResidualHz', inf));
residualTolHz = localGetFieldOrDefault(flowOpt, 'inToothCentralResidualTolHz', 5);
if ~(isfinite(subsetResidualHz) && subsetResidualHz <= residualTolHz)
  return;
end

if ~(isstruct(anchorUnknownSummary) && ~isempty(anchorUnknownSummary) && ...
    logical(localGetFieldOrDefault(anchorUnknownSummary, 'isResolved', false)))
  if logical(localGetFieldOrDefault(flowOpt, 'inToothFreezeDoa', true))
    tf = false;
  end
  return;
end

anchorToothIdx = localGetFieldOrDefault(anchorUnknownSummary, 'toothIdx', NaN);
anchorResidualHz = abs(localGetFieldOrDefault(anchorUnknownSummary, 'toothResidualHz', inf));
if ~(isfinite(anchorToothIdx) && abs(anchorToothIdx) == 0 && ...
    isfinite(anchorResidualHz) && anchorResidualHz <= residualTolHz)
  return;
end

if ~logical(localGetFieldOrDefault(flowOpt, 'enableInToothHealthyAnchorSkip', true))
  return;
end

nonRefFitFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefFitRatioFloor', NaN);
nonRefSupportFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefSupportRatioFloor', NaN);
nonRefConsistencyFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefConsistencyRatioFloor', NaN);
nonRefCoherenceFloor = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefCoherenceFloor', NaN);
nonRefRmsPhaseResidRad = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefRmsPhaseResidRad', inf);
nonRefMaxAbsPhaseResidRad = localGetFieldOrDefault(anchorUnknownSummary, 'nonRefMaxAbsPhaseResidRad', inf);

minFitFloor = localGetFieldOrDefault(flowOpt, 'inToothSkipMinNonRefFitFloor', 0.95);
minSupportFloor = localGetFieldOrDefault(flowOpt, 'inToothSkipMinNonRefSupportFloor', 0.999);
minConsistencyFloor = localGetFieldOrDefault(flowOpt, 'inToothSkipMinNonRefConsistencyFloor', 0.999);
minCoherenceFloor = localGetFieldOrDefault(flowOpt, 'inToothSkipMinNonRefCoherenceFloor', 0.999);
maxRmsPhaseResidRad = localGetFieldOrDefault(flowOpt, 'inToothSkipMaxRmsPhaseResidRad', 0.003);
maxMaxAbsPhaseResidRad = localGetFieldOrDefault(flowOpt, 'inToothSkipMaxAbsPhaseResidRad', 0.005);

if isfinite(nonRefFitFloor) && nonRefFitFloor >= minFitFloor && ...
    isfinite(nonRefSupportFloor) && nonRefSupportFloor >= minSupportFloor && ...
    isfinite(nonRefConsistencyFloor) && nonRefConsistencyFloor >= minConsistencyFloor && ...
    isfinite(nonRefCoherenceFloor) && nonRefCoherenceFloor >= minCoherenceFloor && ...
    isfinite(nonRefRmsPhaseResidRad) && nonRefRmsPhaseResidRad <= maxRmsPhaseResidRad && ...
    isfinite(nonRefMaxAbsPhaseResidRad) && nonRefMaxAbsPhaseResidRad <= maxMaxAbsPhaseResidRad
  maxAnchorDoaDriftFromKnownDeg = localGetFieldOrDefault(flowOpt, 'inToothSkipMaxAnchorDoaDriftFromKnownDeg', inf);
  anchorDoaDriftFromKnownDeg = localCalcSummaryDoaDriftDeg(anchorUnknownSummary, knownSummary);
  if isfinite(anchorDoaDriftFromKnownDeg) && anchorDoaDriftFromKnownDeg > maxAnchorDoaDriftFromKnownDeg
    return;
  end
  tf = false;
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
  error('buildDoaDopplerDynamicTransitionBundle:InvalidLocalOptimOpt', ...
    'Optimizer overrides must be provided as a struct.');
end
fieldList = fieldnames(overrideOptimOpt);
for iField = 1:numel(fieldList)
  optimOpt.(fieldList{iField}) = overrideOptimOpt.(fieldList{iField});
end
end


function summary = localBuildUnknownCaseSummary(caseUse, truth, toothStepHz)
%LOCALBUILDUNKNOWNCASESUMMARY Build one compact tooth-selection summary.

summary = struct();
summary.solveVariant = string(localGetFieldOrDefault(caseUse.estResult, 'solveVariant', "unknown"));
summary.isResolved = logical(localGetFieldOrDefault(caseUse.estResult, 'isResolved', false));
summary.doaParamEst = reshape(caseUse.estResult.doaParamEst(:), 1, []);
summary.fdRefEst = caseUse.estResult.fdRefEst;
summary.fdRateEst = caseUse.estResult.fdRateEst;
summary.angleErrDeg = calcLatlonAngleError(caseUse.estResult.doaParamEst(:), truth.latlonTrueDeg(:));
summary.fdRefErrHz = caseUse.estResult.fdRefEst - truth.fdRefFit;
summary.fdRateErrHzPerSec = caseUse.estResult.fdRateEst - truth.fdRateFit;
summary.toothIdx = round(summary.fdRefErrHz / toothStepHz);
summary.toothResidualHz = summary.fdRefErrHz - summary.toothIdx * toothStepHz;
summary.runTimeMs = localGetFieldOrDefault(caseUse.estResult, 'runTimeMs', NaN);
debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseUse.estResult, 'aux', struct()), 'debug', struct());
finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
summary.finalObj = localGetFieldOrDefault(finalEval, 'obj', NaN);
summary.finalResidualNorm = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
summary.refFitRatio = localGetFieldOrDefault(finalEval, 'refFitRatio', NaN);
summary.nonRefFitRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefFitRatioFloor', NaN);
summary.refSupportRatio = localGetFieldOrDefault(finalEval, 'refSupportRatio', NaN);
summary.nonRefSupportRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefSupportRatioFloor', NaN);
summary.maxNonRefNegativeProjectionRatio = localGetFieldOrDefault(finalEval, 'maxNonRefNegativeProjectionRatio', NaN);
summary.nonRefFitFloorPenalty = localGetFieldOrDefault(finalEval, 'nonRefFitFloorPenalty', 0);
summary.nonRefSupportFloorPenalty = localGetFieldOrDefault(finalEval, 'nonRefSupportFloorPenalty', 0);
summary.additionalObjectivePenalty = localGetFieldOrDefault(finalEval, 'additionalObjectivePenalty', 0);
summary.refConsistencyNorm = localGetFieldOrDefault(finalEval, 'refConsistencyNorm', NaN);
summary.nonRefConsistencyRatioFloor = localGetFieldOrDefault(finalEval, 'nonRefConsistencyRatioFloor', NaN);
summary.isDoaFrozenLike = localIsFrozenLikeSummary(summary.solveVariant);
[summary.refCoherence, summary.nonRefCoherenceFloor, summary.nonRefMaxAbsPhaseResidRad, ...
  summary.nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval);
end

function subsetSeedInfo = localBuildSubsetSeedReuseInfo(bestStaticMsCase, knownSeedCase, staticMsOpt, flowOpt)
%LOCALBUILDSUBSETSEEDREUSEINFO Reuse periodic static/known seeds for subset ranking.

subsetSeedInfo = struct();
subsetSeedInfo.reusePeriodicSeeds = logical(localGetFieldOrDefault(flowOpt, 'enableSubsetCaseSeedReuse', true));
subsetSeedInfo.bestStaticMsCase = bestStaticMsCase;
subsetSeedInfo.knownSeedCase = knownSeedCase;
subsetSeedInfo.staticInitDoaHalfWidth = localGetFieldOrDefault(staticMsOpt, 'initDoaHalfWidth', [0.01; 0.01]);
end


function summary = localBuildSkippedUnknownSummary(stageTag)
%LOCALBUILDSKIPPEDUNKNOWNSUMMARY Build one placeholder summary for skipped branches.

summary = struct();
summary.solveVariant = "skipped";
summary.isResolved = false;
summary.doaParamEst = nan(1, 2);
summary.fdRefEst = NaN;
summary.fdRateEst = NaN;
summary.angleErrDeg = NaN;
summary.fdRefErrHz = NaN;
summary.fdRateErrHzPerSec = NaN;
summary.toothIdx = NaN;
summary.toothResidualHz = NaN;
summary.runTimeMs = 0;
summary.finalObj = NaN;
summary.finalResidualNorm = NaN;
summary.stageTag = string(stageTag);
end


function tf = localCaseHasUsableEstimate(caseUse)
%LOCALCASEHASUSABLEESTIMATE Return true when one case contains finite key estimates.

tf = false;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult) || ~isstruct(estResult)
  return;
end
fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
doaParamEst = localGetFieldOrDefault(estResult, 'doaParamEst', nan(2, 1));
tf = isfinite(fdRefEst) && isfinite(fdRateEst) && all(isfinite(doaParamEst(:)));
end


function fdRateEst = localGetCaseFdRateEst(caseUse, defaultValue)
%LOCALGETCASEFDRATEEST Safely read fdRateEst from one case struct.

if nargin < 2 || isempty(defaultValue)
  defaultValue = 0;
end
fdRateEst = defaultValue;
if isempty(caseUse) || ~isstruct(caseUse)
  return;
end
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
if isempty(estResult)
  return;
end
fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', defaultValue);
if isempty(fdRateEst) || ~isscalar(fdRateEst) || ~isfinite(fdRateEst)
  fdRateEst = defaultValue;
end
end

function [refCoherence, nonRefCoherenceFloor, nonRefMaxAbsPhaseResidRad, nonRefRmsPhaseResidRad] = localExtractNonRefEvalMetrics(finalEval)
%LOCALEXTRACTNONREFEVALMETRICS Build compact truth-free non-reference metrics.

refCoherence = NaN;
nonRefCoherenceFloor = NaN;
nonRefMaxAbsPhaseResidRad = NaN;
nonRefRmsPhaseResidRad = NaN;
if ~isstruct(finalEval)
  return;
end
refSatIdxLocal = localGetFieldOrDefault(finalEval, 'refSatIdxLocal', NaN);
coherenceSat = reshape(localGetFieldOrDefault(finalEval, 'coherenceSat', []), [], 1);
if isfinite(refSatIdxLocal) && ~isempty(coherenceSat) && refSatIdxLocal >= 1 && refSatIdxLocal <= numel(coherenceSat)
  refCoherence = coherenceSat(refSatIdxLocal);
  nonRefMask = true(size(coherenceSat));
  nonRefMask(refSatIdxLocal) = false;
  nonRefVals = coherenceSat(nonRefMask & isfinite(coherenceSat));
  if ~isempty(nonRefVals)
    nonRefCoherenceFloor = min(nonRefVals);
  end
end
blockPhaseResidMat = localGetFieldOrDefault(finalEval, 'blockPhaseResidMat', []);
if isempty(blockPhaseResidMat)
  return;
end
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= size(blockPhaseResidMat, 1))
  return;
end
nonRefMask = true(size(blockPhaseResidMat, 1), 1);
nonRefMask(refSatIdxLocal) = false;
nonRefPhaseResid = abs(blockPhaseResidMat(nonRefMask, :));
nonRefPhaseResid = nonRefPhaseResid(isfinite(nonRefPhaseResid));
if isempty(nonRefPhaseResid)
  return;
end
nonRefMaxAbsPhaseResidRad = max(nonRefPhaseResid);
nonRefRmsPhaseResidRad = sqrt(mean(nonRefPhaseResid .^ 2));
end


function tf = localIsFrozenLikeSummary(solveVariant)
%LOCALISFROZENLIKESUMMARY Return true for fixed-DoA-style summary variants.

solveVariant = lower(strtrim(string(solveVariant)));
tf = contains(solveVariant, "fixeddoa");
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
%LOCALBUILDTRUTHFDLOCALFORVIEW Build one single-frame truth local-frequency vector.

numSat = sceneRef.numSat;
fdLocalTrue = nan(numSat, 1);
selectedSatIdxGlobal = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
fdSatTrueHz = reshape(localGetFieldOrDefault(truth, 'fdSatTrueHz', []), [], 1);
if isempty(selectedSatIdxGlobal) || isempty(fdSatTrueHz)
  return;
end

sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
for iSat = 1:min(numSat, numel(sceneSatIdx))
  matchIdx = find(selectedSatIdxGlobal == sceneSatIdx(iSat), 1, 'first');
  if ~isempty(matchIdx) && matchIdx <= numel(fdSatTrueHz)
    fdLocalTrue(iSat) = fdSatTrueHz(matchIdx);
  end
end
end


function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
%LOCALEXTRACTSCENELOCALDOA Stack one scene sequence local-DoA tensor.

numSat = sceneSeq.numSat;
numFrame = sceneSeq.numFrame;
truthLocalDoa = nan(2, numSat, numFrame);
for iFrame = 1:numFrame
  sceneUse = sceneSeq.sceneCell{iFrame};
  if isfield(sceneUse, 'localDoa') && ~isempty(sceneUse.localDoa)
    localDoaUse = sceneUse.localDoa;
    truthLocalDoa(:, 1:size(localDoaUse, 2), iFrame) = localDoaUse;
  end
end
end


function out = localMergeStruct(base, override)
%LOCALMERGESTRUCT Merge two scalar structs with override precedence.

out = base;
if nargin < 2 || isempty(override)
  return;
end
if ~isstruct(override)
  error('buildDoaDopplerDynamicTransitionBundle:InvalidOverride', ...
    'Override data must be a struct.');
end
fieldList = fieldnames(override);
for iField = 1:numel(fieldList)
  out.(fieldList{iField}) = override.(fieldList{iField});
end
end


function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row compactly.

if isempty(valueRow)
  textOut = '[]';
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '), ']'];
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end
try
  fieldValue = dataStruct.(fieldName);
catch
end
end
