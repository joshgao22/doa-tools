function regressionMfSubsetSelectNoTruthLeak(varargin)
% Regression check for tooth-first, truth-free subset selection.
% This script verifies one contract only:
%   1) subset selection must not depend on truth-aware error fields;
%   2) among candidates on the same tooth, ranking must first use phase
%      health gates, then objective / residual, and never truth-aware errors;
%   3) the selected subset is trusted only when the chosen candidate is a
%      resolved fit with finite objective diagnostics.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fprintf('Running regressionMfSubsetSelectNoTruthLeak ...\n');

summaryCell = { ...
  localBuildSummary(true,  -1.00e6, 6.55e4,  1.0e-4,   0.01,    0, 12.0), ...
  localBuildSummary(true,  -1.10e6, 6.60e4,  9.9e+1, 5000.0,    0, 11.0), ...
  localBuildSummary(false, -1.20e6, 1.00e4,  0.0,      0.0,    0, 10.0) ...
  };
[bestIdx, isTrusted, reasonText] = selectBestDynamicSubsetSummary(summaryCell);

if bestIdx ~= 2
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedBestIdx', ...
    'Expected same-tooth selection to ignore truth-aware errors, but bestIdx=%d.', bestIdx);
end
if ~isTrusted
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedTrustFlag', ...
    'The selected resolved candidate with finite diagnostics must be trusted.');
end

healthCell = { ...
  localBuildSummary(true, -1.00e6, 4.0e4, 10.0, 1.0e4, 0, 20.0, 1.00), ...
  localBuildSummary(true, -1.50e6, 3.0e4,  0.0, 0.0,   0, 10.0, 0.50) ...
  };
[healthBestIdx, healthTrusted] = selectBestDynamicSubsetSummary(healthCell);
if healthBestIdx ~= 1 || ~healthTrusted
  error('regressionMfSubsetSelectNoTruthLeak:HealthGateBypassed', ...
    ['A same-tooth candidate with bad phase-health diagnostics must not ', ...
     'beat a healthy candidate only because its objective is lower.']);
end

bestScore = buildDynamicSelectionScoreVector(summaryCell{bestIdx});
firstScore = buildDynamicSelectionScoreVector(summaryCell{1});
thirdScore = buildDynamicSelectionScoreVector(summaryCell{3});

fprintf('  truth-free selected idx     : %d\n', bestIdx);
fprintf('  truth-free selected reason  : %s\n', reasonText);
fprintf('  health-gated selected idx   : %d\n', healthBestIdx);
fprintf('  selected subset trusted     : %d\n', isTrusted);
fprintf('  score cand1                 : %s\n', localFormatNumericRow(firstScore));
fprintf('  score cand2                 : %s\n', localFormatNumericRow(bestScore));
fprintf('  score cand3                 : %s\n', localFormatNumericRow(thirdScore));
fprintf('PASS: regressionMfSubsetSelectNoTruthLeak\n');
end

function summary = localBuildSummary(isResolved, finalObj, finalResidualNorm, angleErrDeg, fdRefErrHz, toothIdx, runTimeMs, qualityFloor)
%LOCALBUILDSUMMARY Build one synthetic subset summary.

if nargin < 8
  qualityFloor = 1.0;
end
summary = struct();
summary.isResolved = isResolved;
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.angleErrDeg = angleErrDeg;
summary.fdRefErrHz = fdRefErrHz;
summary.fdRateErrHzPerSec = fdRefErrHz;
summary.toothIdx = toothIdx;
summary.toothResidualHz = fdRefErrHz;
summary.runTimeMs = runTimeMs;
summary.nonRefSupportRatioFloor = qualityFloor;
summary.nonRefFitRatioFloor = qualityFloor;
summary.nonRefConsistencyRatioFloor = qualityFloor;
summary.nonRefCoherenceFloor = qualityFloor;
summary.nonRefRmsPhaseResidRad = 1e-3 / max(qualityFloor, 0.1);
summary.nonRefMaxAbsPhaseResidRad = 2e-3 / max(qualityFloor, 0.1);
summary.maxNonRefNegativeProjectionRatio = max(0, 1 - qualityFloor) * 0.2;
end

function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end
