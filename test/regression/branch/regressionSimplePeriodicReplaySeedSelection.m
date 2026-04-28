function regressionSimplePeriodicReplaySeedSelection(varargin)
% Regression check for the simple-flow periodic replay seed/fallback selection.
% The simplified subset-periodic flow should:
%   1) trust the subset replay when the subset seed is close and trusted;
%   2) fall back to the static replay when the subset drifts too far;
%   3) fall back to the static replay on same-tooth near-ties with a large
%      frozen DoA disagreement.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
localPrint = @(varargin) fprintf(varargin{:});
baseFlowOpt = buildSimpleDynamicFlowOpt();
baseFlowOpt.verbose = verbose;
baseFlowOpt.periodicRefineEnableVerySmallDoaPolish = false;
baseFlowOpt.periodicRefineDoaSeedMode = "dualWhenMulti";
baseFlowOpt.periodicRefineFreezeDoa = true;

staticDoa = [37.7800; 36.5900];
closeSubsetDoa = [37.7808; 36.5903];
farSubsetDoa = [37.7865; 36.5935];

dummyFixture = struct('sceneSeq', struct('refFrameIdx', 1), 'viewMs', struct('numSat', 2));
dummyTruth = struct();
dummySubsetTrust = struct('isTrusted', true, 'selectionReason', "objective-only", ...
  'relativeMargin', 1e-3, 'runnerUpLabel', "runner-up");

localPrint('Running regressionSimplePeriodicReplaySeedSelection\n');

% Case 1: trusted close subset replay should win.
caseStatic = localMakeCase(staticDoa, -1200, -700, true);
caseSubsetClose = localMakeCase(closeSubsetDoa, -1200, -700, true);
flowOptTrusted = baseFlowOpt;
seedList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOptTrusted, "multi", caseStatic, caseSubsetClose);
[finalCaseTrusted, periodicDiagTrusted] = selectSimpleDynamicFinalCase( ...
  "MS-MF-Dyn", "multi", dummyFixture, caseSubsetClose, dummySubsetTrust, ...
  750, seedList, {caseSubsetClose; caseStatic}, ...
  {localMakeSummary(100, 10, true), localMakeSummary(101, 10, true)}, ...
  {localMakeSummary(100, 10, true), localMakeSummary(101, 10, true)}, ...
  [], 1, 1, false, flowOptTrusted, struct(), dummyTruth);
if ~isequal(finalCaseTrusted.estResult.doaParamEst(:), caseSubsetClose.estResult.doaParamEst(:))
  error('regressionSimplePeriodicReplaySeedSelection:SubsetReplayDidNotWin', ...
    'Trusted close subset replay should win over the static replay.');
end
if string(periodicDiagTrusted.selectedReplaySeedSource) ~= "selected-subset" || ...
    string(periodicDiagTrusted.selectedReplaySelectionReason) ~= "best-frozen-subset-trusted"
  error('regressionSimplePeriodicReplaySeedSelection:UnexpectedTrustedSelectionReason', ...
    'Unexpected trusted selection reason: %s', char(string(periodicDiagTrusted.selectedReplaySelectionReason)));
end
localPrint('  trusted subset winner          : %s (%s)\n', ...
  char(string(periodicDiagTrusted.selectedReplaySeedSource)), ...
  char(string(periodicDiagTrusted.selectedReplaySelectionReason)));

% Case 2: subset drift is too large, so static replay must be used.
caseSubsetFar = localMakeCase(farSubsetDoa, -1200, -700, true);
flowOptLargeDrift = baseFlowOpt;
flowOptLargeDrift.periodicRefineMaxSubsetDoaDriftDeg = 0.003;
seedList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOptLargeDrift, "multi", caseStatic, caseSubsetFar);
[finalCaseDrift, periodicDiagDrift] = selectSimpleDynamicFinalCase( ...
  "MS-MF-Dyn", "multi", dummyFixture, caseSubsetFar, dummySubsetTrust, ...
  750, seedList, {caseSubsetFar; caseStatic}, ...
  {localMakeSummary(100, 10, true), localMakeSummary(101, 10, true)}, ...
  {localMakeSummary(100, 10, true), localMakeSummary(101, 10, true)}, ...
  [], 1, 1, false, flowOptLargeDrift, struct(), dummyTruth);
if ~isequal(finalCaseDrift.estResult.doaParamEst(:), caseStatic.estResult.doaParamEst(:))
  error('regressionSimplePeriodicReplaySeedSelection:LargeSubsetDriftDidNotFallback', ...
    'Large subset drift should force a fallback to the static replay.');
end
if string(periodicDiagDrift.selectedReplaySeedSource) ~= "static-seed" || ...
    string(periodicDiagDrift.selectedReplaySelectionReason) ~= "subset-static-drift-too-large"
  error('regressionSimplePeriodicReplaySeedSelection:UnexpectedLargeDriftReason', ...
    'Unexpected large-drift selection reason: %s', char(string(periodicDiagDrift.selectedReplaySelectionReason)));
end
localPrint('  large-drift fallback           : %s (%s)\n', ...
  char(string(periodicDiagDrift.selectedReplaySeedSource)), ...
  char(string(periodicDiagDrift.selectedReplaySelectionReason)));

% Case 3: same-tooth near tie with a large DoA disagreement should prefer the static replay.
flowOptTie = baseFlowOpt;
flowOptTie.periodicRefineMaxSubsetDoaDriftDeg = 0.02;
flowOptTie.periodicRefineMaxFrozenDoaDisagreementDeg = 0.002;
flowOptTie.periodicRefineMaxFrozenRelativeObjGapForTie = 5e-4;
seedList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOptTie, "multi", caseStatic, caseSubsetFar);
[finalCaseTie, periodicDiagTie] = selectSimpleDynamicFinalCase( ...
  "MS-MF-Dyn", "multi", dummyFixture, caseSubsetFar, dummySubsetTrust, ...
  750, seedList, {caseSubsetFar; caseStatic}, ...
  {localMakeSummary(100.0000, 10, true), localMakeSummary(100.00001, 10, true)}, ...
  {localMakeSummary(100.0000, 10, true), localMakeSummary(100.00001, 10, true)}, ...
  [], 1, 1, false, flowOptTie, struct(), dummyTruth);
if ~isequal(finalCaseTie.estResult.doaParamEst(:), caseStatic.estResult.doaParamEst(:))
  error('regressionSimplePeriodicReplaySeedSelection:SameToothTieDidNotFallback', ...
    'Same-tooth near tie with a large DoA disagreement should pick the static replay.');
end
if string(periodicDiagTie.selectedReplaySeedSource) ~= "static-seed" || ...
    string(periodicDiagTie.selectedReplaySelectionReason) ~= "frozen-subset-static-tie-with-large-doa-gap"
  error('regressionSimplePeriodicReplaySeedSelection:UnexpectedTieReason', ...
    'Unexpected same-tooth tie selection reason: %s', char(string(periodicDiagTie.selectedReplaySelectionReason)));
end
localPrint('  same-tooth tie fallback        : %s (%s)\n', ...
  char(string(periodicDiagTie.selectedReplaySeedSource)), ...
  char(string(periodicDiagTie.selectedReplaySelectionReason)));

localPrint('PASS: regressionSimplePeriodicReplaySeedSelection\n');

end

function caseUse = localMakeCase(doaParam, fdRefEst, fdRateEst, isResolved)
caseUse = struct();
caseUse.estResult = struct();
caseUse.estResult.doaParamEst = reshape(doaParam(:), [], 1);
caseUse.estResult.fdRefEst = fdRefEst;
caseUse.estResult.fdRateEst = fdRateEst;
caseUse.estResult.isResolved = logical(isResolved);
end

function summary = localMakeSummary(finalObj, finalResidualNorm, isResolved)
summary = struct();
summary.isResolved = logical(isResolved);
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.nonRefSupportRatioFloor = 1;
summary.nonRefFitRatioFloor = 1;
summary.nonRefConsistencyRatioFloor = 1;
summary.nonRefCoherenceFloor = 1;
summary.nonRefRmsPhaseResidRad = 1e-3;
summary.nonRefMaxAbsPhaseResidRad = 2e-3;
summary.maxNonRefNegativeProjectionRatio = 1e-3;
summary.toothIdx = 0;
summary.toothResidualHz = 5;
summary.runTimeMs = 1;
end
