% Regression check for one small-DoA wide refine after fd-only wide escape.
% This script verifies one contract only:
%   1) periodic-fd-wide may first escape to the correct tooth with DoA frozen;
%   2) periodic-wide may then make a small same-tooth objective improvement;
%   3) final selection should prefer the refined periodic-wide candidate.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownWideRefineFromFdWide ...\n');

wideSummary = localBuildSummary(true, -1.50845e6, 6.5595e4, [37.7798, 36.5902], -28940.0, -3833.5, 0, 13.0);
fdWideSummary = localBuildSummary(true, -1.50840e6, 6.5600e4, [37.7795, 36.5905], -28940.1, -3833.6, 0, 11.0);
anchorSummary = localBuildSummary(true, -1.50810e6, 6.5630e4, [37.7890, 36.5990], -28190.0, -3834.0, 1, 10.0);
refineSummary = localBuildSummary(true, -1.50835e6, 6.5605e4, [37.7797, 36.5904], -28939.9, -3833.7, 0, 18.0);
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

if selectedIdx ~= 1 || selectedTag ~= "periodic-wide"
  error('regressionMfUnknownWideRefineFromFdWide:UnexpectedSelection', ...
    'Expected periodic-wide to win after improving periodic-fd-wide, but got %s.', selectedTag);
end
if contains(reasonText, "subset-anchor guard")
  error('regressionMfUnknownWideRefineFromFdWide:UnexpectedGuardReason', ...
    'The small-DoA wide refine should not be blocked by the subset-anchor guard.');
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-fd-wide     : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('  score periodic-fd-anchor   : %s\n', localFormatNumericRow(scoreMat(3, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(4, :)));
fprintf('PASS: regressionMfUnknownWideRefineFromFdWide\n');


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
