% Regression check for one closer-to-zero-tooth wide escape.
% This script verifies one contract only:
%   1) the subset anchor may stay on a wrong nonzero tooth;
%   2) periodic-wide may move back toward the central tooth with a better
%      objective/residual pair;
%   3) the subset-anchor guard must not force the final pick back to anchor.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownCloserToZeroToothWins ...\n');

wideSummary = localBuildSummary(true, -1.5081e6, 6.5632e4, [37.7800, 36.5900], -28940.01, -3844.6, 0, 12.0);
anchorSummary = localBuildSummary(true, -1.5052e6, 6.5920e4, [37.7490, 36.5870], -30440.0, -3883.5, -2, 11.0);
refineSummary = localBuildSummary(true, -1.5052e6, 6.5922e4, [37.7492, 36.5872], -30440.0, -3883.6, -2, 18.0);
finalSelectOpt = struct( ...
  'enableSubsetAnchorGuard', true, ...
  'subsetAnchorIdx', 2, ...
  'subsetAnchorTag', "periodic-fd-anchor", ...
  'maxFdRefDriftHz', 100, ...
  'maxFdRateDriftHzPerSec', 100, ...
  'maxDoaDriftDeg', 0.01, ...
  'minObjGainToLeaveAnchor', 1e4, ...
  'minResidualGainToLeaveAnchor', 1e3);
[selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary( ...
  {wideSummary; anchorSummary; refineSummary}, ...
  ["periodic-wide"; "periodic-fd-anchor"; "periodic-in-tooth"], finalSelectOpt);

if selectedIdx ~= 1 || selectedTag ~= "periodic-wide"
  error('regressionMfUnknownCloserToZeroToothWins:UnexpectedSelection', ...
    'Expected periodic-wide to win when it moves closer to the central tooth, but got %s.', selectedTag);
end
if contains(reasonText, "subset-anchor guard")
  error('regressionMfUnknownCloserToZeroToothWins:UnexpectedGuardReason', ...
    'Closer-to-zero-tooth wide escape should not be blocked by the subset-anchor guard.');
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-fd-anchor   : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(3, :)));
fprintf('PASS: regressionMfUnknownCloserToZeroToothWins\n');


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
