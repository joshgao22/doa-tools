function regressionSimpleFlowPolishTrigger(varargin)
% Regression check for the simple-flow very-small DoA polish gate.
% The simplified subset-periodic flow should only trigger the extra DoA
% polish on trusted same-tooth hard cases. Easy or clearly unsafe cases
% must stay on the bypass path.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
localPrint = @(varargin) fprintf(varargin{:});
flowOpt = buildSimpleDynamicFlowOpt();
flowOpt.verbose = verbose;
flowOpt.periodicRefineEnableVerySmallDoaPolish = true;
flowOpt.periodicPolishEnableWhenMulti = true;

hardSeed = localBuildSeedDiag("selected-subset", 4e-4, 8e-4, 2e-4);
hardSummary = localBuildSummary(0, 2, true, 1.0);
[runHard, reasonHard, diagHard] = shouldRunSimpleVerySmallDoaPolish("multi", hardSeed, hardSummary, flowOpt);

if ~runHard || string(reasonHard) ~= "triggered"
  error('regressionSimpleFlowPolishTrigger:HardCaseDidNotTrigger', ...
    'Trusted same-tooth hard case should trigger polish. reason=%s', char(string(reasonHard)));
end
if diagHard.selectedReplayHealthBucket ~= 0
  error('regressionSimpleFlowPolishTrigger:UnexpectedHealthBucket', ...
    'The healthy hard case should stay in bucket 0 before polish.');
end

% Easy case: the two frozen seeds nearly agree, so there is nothing to polish.
easySeed = localBuildSeedDiag("selected-subset", 4e-4, 1e-4, 2e-4);
[runEasy, reasonEasy] = shouldRunSimpleVerySmallDoaPolish("multi", easySeed, hardSummary, flowOpt);
if runEasy || string(reasonEasy) ~= "frozen-doa-gap-too-small"
  error('regressionSimpleFlowPolishTrigger:EasyCaseTriggered', ...
    'Easy same-tooth case should stay bypassed. reason=%s', char(string(reasonEasy)));
end

% Unsafe case: fd is not yet healthy, so the polish must remain disabled.
unsafeSummary = localBuildSummary(0, 2, true, 0.96);
[runUnsafe, reasonUnsafe] = shouldRunSimpleVerySmallDoaPolish("multi", hardSeed, unsafeSummary, flowOpt);
if runUnsafe || string(reasonUnsafe) ~= "fd-health-not-ready"
  error('regressionSimpleFlowPolishTrigger:UnsafeCaseTriggered', ...
    'Unhealthy same-tooth case should not trigger polish. reason=%s', char(string(reasonUnsafe)));
end

% Single-satellite simple flow must never use the multi-satellite polish path.
[runSingle, reasonSingle] = shouldRunSimpleVerySmallDoaPolish("single", hardSeed, hardSummary, flowOpt);
if runSingle || string(reasonSingle) ~= "single-mode"
  error('regressionSimpleFlowPolishTrigger:SingleModeTriggered', ...
    'Single-satellite simple flow should bypass polish. reason=%s', char(string(reasonSingle)));
end

localPrint('Running regressionSimpleFlowPolishTrigger\n');
localPrint('  hard-case reason              : %s\n', char(string(reasonHard)));
localPrint('  easy-case bypass reason       : %s\n', char(string(reasonEasy)));
localPrint('  unhealthy-case bypass reason  : %s\n', char(string(reasonUnsafe)));
localPrint('  single-mode bypass reason     : %s\n', char(string(reasonSingle)));
localPrint('PASS: regressionSimpleFlowPolishTrigger\n');

end

function seedDiag = localBuildSeedDiag(seedSource, subsetDriftDeg, frozenDoaGapDeg, frozenObjGap)
seedDiag = struct();
seedDiag.selectedReplaySeedSource = string(seedSource);
seedDiag.subsetDriftFromStaticDeg = subsetDriftDeg;
seedDiag.subsetRankMarginRelative = 1e-3;
seedDiag.frozenDoaDisagreementDeg = frozenDoaGapDeg;
seedDiag.frozenRelativeObjGap = frozenObjGap;
end

function summary = localBuildSummary(toothIdx, toothResidualHz, isResolved, qualityFloor)
summary = struct();
summary.isResolved = logical(isResolved);
summary.toothIdx = toothIdx;
summary.toothResidualHz = toothResidualHz;
summary.nonRefSupportRatioFloor = qualityFloor;
summary.nonRefFitRatioFloor = qualityFloor;
summary.nonRefConsistencyRatioFloor = qualityFloor;
summary.nonRefCoherenceFloor = qualityFloor;
summary.nonRefRmsPhaseResidRad = 1e-3;
summary.nonRefMaxAbsPhaseResidRad = 2e-3;
summary.maxNonRefNegativeProjectionRatio = 1e-2;
end
