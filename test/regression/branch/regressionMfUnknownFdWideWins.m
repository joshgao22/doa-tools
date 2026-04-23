% Regression check for one fd-only wide candidate winning the final selection.
% This script verifies one contract only:
%   1) periodic-fd-wide is allowed to enter the final candidate pool;
%   2) when periodic-fd-wide clearly beats both the subset anchor and the
%      post-wide refine by objective/residual, final selection may keep the
%      fixed-DoA tooth-escape result;
%   3) the subset-anchor guard must not force the final pick back to anchor.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownFdWideWins ...\n');

wideSummary = localBuildSummary(true, -1.5079e6, 6.5660e4, [37.7810, 36.5920], -28939.9, -3833.5, 0, 12.0);
fdWideSummary = localBuildSummary(true, -1.5084e6, 6.5600e4, [37.7795, 36.5905], -28940.0, -3833.6, 0, 11.0);
anchorSummary = localBuildSummary(true, -1.5081e6, 6.5635e4, [37.7890, 36.5990], -28190.0, -3834.0, 1, 10.0);
refineSummary = localBuildSummary(true, -1.5080e6, 6.5640e4, [37.7808, 36.5917], -28939.8, -3834.1, 0, 18.0);
finalSelectOpt = struct( ...
  'enableSubsetAnchorGuard', true, ...
  'subsetAnchorIdx', 3, ...
  'subsetAnchorTag', "periodic-fd-anchor", ...
  'maxFdRefDriftHz', 100, ...
  'maxFdRateDriftHzPerSec', 100, ...
  'maxDoaDriftDeg', 0.01, ...
  'minObjGainToLeaveAnchor', 1e4, ...
  'minResidualGainToLeaveAnchor', 1e3);
[selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary( ...
  {wideSummary; fdWideSummary; anchorSummary; refineSummary}, ...
  ["periodic-wide"; "periodic-fd-wide"; "periodic-fd-anchor"; "periodic-in-tooth"], finalSelectOpt);

if selectedIdx ~= 2 || selectedTag ~= "periodic-fd-wide"
  error('regressionMfUnknownFdWideWins:UnexpectedSelection', ...
    'Expected periodic-fd-wide to win, but got %s.', selectedTag);
end
if contains(reasonText, "subset-anchor guard")
  error('regressionMfUnknownFdWideWins:UnexpectedGuardReason', ...
    'fd-only wide win should not be blocked by the subset-anchor guard.');
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-fd-wide     : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('  score periodic-fd-anchor   : %s\n', localFormatNumericRow(scoreMat(3, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(4, :)));
fprintf('PASS: regressionMfUnknownFdWideWins\n');


function summary = localBuildSummary(isResolved, finalObj, finalResidualNorm, doaParamEst, fdRefEst, fdRateEst, toothIdx, runTimeMs)
summary = struct();
summary.isResolved = isResolved;
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.doaParamEst = reshape(doaParamEst, 1, []);
summary.fdRefEst = fdRefEst;
summary.fdRateEst = fdRateEst;
summary.angleErrDeg = 0;
summary.fdRefErrHz = 0;
summary.fdRateErrHzPerSec = 0;
summary.toothIdx = toothIdx;
summary.toothResidualHz = 0;
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
