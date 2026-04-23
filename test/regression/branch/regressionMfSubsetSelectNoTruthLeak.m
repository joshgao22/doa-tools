% Regression check for tooth-first, truth-free subset selection.
% This script verifies one contract only:
%   1) subset selection must not depend on truth-aware error fields;
%   2) among candidates on the same tooth, ranking still follows final
%      objective / residual rather than truth-aware error fields;
%   3) the selected subset is trusted only when the chosen candidate is a
%      resolved fit with finite objective diagnostics.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfSubsetSelectNoTruthLeak ...\n');

summaryCell = { ...
  localBuildSummary(true,  -1.00e6, 6.55e4,  1.0e-4,   0.01,    0, 12.0), ...
  localBuildSummary(true,  -1.10e6, 6.60e4,  9.9e+1, 5000.0,    0, 11.0), ...
  localBuildSummary(false, -1.20e6, 1.00e4,  0.0,      0.0,    0, 10.0) ...
  };
[bestIdx, isTrusted, reasonText] = selectBestDynamicSubsetSummary(summaryCell);

if bestIdx ~= 2
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedBestIdx', ...
    'Expected the selection to ignore truth-aware errors within the same tooth, but bestIdx=%d.', bestIdx);
end
if ~isTrusted
  error('regressionMfSubsetSelectNoTruthLeak:UnexpectedTrustFlag', ...
    'The selected resolved candidate with finite diagnostics must be trusted.');
end

bestScore = buildDynamicSelectionScoreVector(summaryCell{bestIdx});
firstScore = buildDynamicSelectionScoreVector(summaryCell{1});
thirdScore = buildDynamicSelectionScoreVector(summaryCell{3});

fprintf('  selected subset idx        : %d\n', bestIdx);
fprintf('  selected subset reason     : %s\n', reasonText);
fprintf('  selected subset trusted    : %d\n', isTrusted);
fprintf('  score cand1                : %s\n', localFormatNumericRow(firstScore));
fprintf('  score cand2                : %s\n', localFormatNumericRow(bestScore));
fprintf('  score cand3                : %s\n', localFormatNumericRow(thirdScore));
fprintf('PASS: regressionMfSubsetSelectNoTruthLeak\n');


function summary = localBuildSummary(isResolved, finalObj, finalResidualNorm, angleErrDeg, fdRefErrHz, toothIdx, runTimeMs)
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
end


function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
