function regressionMfSubsetSelectNoTruthLeak(varargin)
% Regression check for truth-free subset selection.
% This script verifies one contract only:
%   1) subset selection must not depend on truth-aware error or tooth fields;
%   2) ranking must first use phase-health gates, then objective / residual;
%   3) the selected subset is trusted only when the chosen candidate is a
%      resolved fit with finite objective diagnostics and healthy fit metrics.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fprintf('Running regressionMfSubsetSelectNoTruthLeak ...\n');

summaryCell = { ...
  localBuildSummary(true,  -1.00e6, 6.55e4,  1.0e-4,      0.0,   0, 12.0), ...
  localBuildSummary(true,  -1.10e6, 6.60e4,  9.9e+1,  5000.0,  11, 11.0), ...
  localBuildSummary(false, -1.20e6, 1.00e4,  0.0,         0.0,   0, 10.0) ...
  };
[bestIdx, isTrusted, reasonText] = selectBestDynamicSubsetSummary(summaryCell);

if bestIdx ~= 2
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedBestIdx', ...
    'Expected selection to ignore truth-aware errors, but bestIdx=%d.', bestIdx);
end
if ~isTrusted
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedTrustFlag', ...
    'The selected resolved candidate with finite diagnostics must be trusted.');
end

shiftedTruthCell = localShiftTruthOnly(summaryCell);
[shiftedBestIdx, shiftedTrusted] = selectBestDynamicSubsetSummary(shiftedTruthCell);
if shiftedBestIdx ~= bestIdx || shiftedTrusted ~= isTrusted
  error('regressionMfSubsetSelectNoTruthLeak:TruthShiftChangedDecision', ...
    'Changing only truth-aware fields must not change subset selection or trust.');
end

healthCell = { ...
  localBuildSummary(true, -1.00e6, 4.0e4, 10.0, 1.0e4, -7, 20.0, 1.00), ...
  localBuildSummary(true, -1.50e6, 3.0e4,  0.0, 0.0,    0, 10.0, 0.50) ...
  };
[healthBestIdx, healthTrusted] = selectBestDynamicSubsetSummary(healthCell);
if healthBestIdx ~= 1 || ~healthTrusted
  error('regressionMfSubsetSelectNoTruthLeak:HealthGateBypassed', ...
    ['A candidate with bad phase-health diagnostics must not beat a healthy ', ...
     'candidate only because its objective is lower.']);
end

bestScore = buildDynamicSelectionScoreVector(summaryCell{bestIdx});
firstScore = buildDynamicSelectionScoreVector(summaryCell{1});
thirdScore = buildDynamicSelectionScoreVector(summaryCell{3});

fprintf('  truth-free selected idx     : %d\n', bestIdx);
fprintf('  shifted-truth selected idx  : %d\n', shiftedBestIdx);
fprintf('  truth-free selected reason  : %s\n', reasonText);
fprintf('  health-gated selected idx   : %d\n', healthBestIdx);
fprintf('  selected subset trusted     : %d\n', isTrusted);
if verbose
  fprintf('  score cand1                 : %s\n', localFormatNumericRow(firstScore));
  fprintf('  score cand2                 : %s\n', localFormatNumericRow(bestScore));
  fprintf('  score cand3                 : %s\n', localFormatNumericRow(thirdScore));
end
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
summary.runTimeMs = runTimeMs;
summary.nonRefSupportRatioFloor = qualityFloor;
summary.nonRefFitRatioFloor = qualityFloor;
summary.nonRefConsistencyRatioFloor = qualityFloor;
summary.nonRefCoherenceFloor = qualityFloor;
summary.nonRefRmsPhaseResidRad = 1e-3 / max(qualityFloor, 0.1);
summary.nonRefMaxAbsPhaseResidRad = 2e-3 / max(qualityFloor, 0.1);
summary.maxNonRefNegativeProjectionRatio = max(0, 1 - qualityFloor) * 0.2;
summary.angleErrDeg = angleErrDeg;
summary.fdRefErrHz = fdRefErrHz;
summary.fdRateErrHzPerSec = fdRefErrHz;
summary.toothIdx = toothIdx;
summary.toothResidualHz = fdRefErrHz;
summary.truthFdRefErrHz = fdRefErrHz;
summary.truthFdRateErrHzPerSec = fdRefErrHz;
summary.truthToothIdx = toothIdx;
summary.truthToothResidualHz = fdRefErrHz;
end

function shiftedCell = localShiftTruthOnly(summaryCell)
%LOCALSHIFTTRUTHONLY Perturb truth-aware fields while keeping decision fields fixed.

shiftedCell = summaryCell;
for iCell = 1:numel(shiftedCell)
  shiftedCell{iCell}.angleErrDeg = 1e3 - iCell;
  shiftedCell{iCell}.fdRefErrHz = 1e5 * iCell;
  shiftedCell{iCell}.fdRateErrHzPerSec = -1e5 * iCell;
  shiftedCell{iCell}.toothIdx = 100 - iCell;
  shiftedCell{iCell}.toothResidualHz = 250 + iCell;
  shiftedCell{iCell}.truthFdRefErrHz = shiftedCell{iCell}.fdRefErrHz;
  shiftedCell{iCell}.truthFdRateErrHzPerSec = shiftedCell{iCell}.fdRateErrHzPerSec;
  shiftedCell{iCell}.truthToothIdx = shiftedCell{iCell}.toothIdx;
  shiftedCell{iCell}.truthToothResidualHz = shiftedCell{iCell}.toothResidualHz;
end
end

function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one score vector for verbose output.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end
