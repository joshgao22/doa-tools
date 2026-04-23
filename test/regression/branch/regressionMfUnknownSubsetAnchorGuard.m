% Regression check for final unknown subset-anchor guard.
% This script verifies one contract only:
%   1) final unknown selection may still use objective-only ranking;
%   2) when one candidate leaves the subset-anchor tooth by a large amount
%      but only wins by a very small objective/residual margin, the guard
%      should keep the subset-anchor replay.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownSubsetAnchorGuard ...\n');

wideSummary = localBuildSummary(true, -1.5086e6, 6.5576e4, [37.7900, 36.6000], -29190.0, -3826.0, 1, 12.0);
anchorSummary = localBuildSummary(true, -1.50855e6, 6.5580e4, [37.7800, 36.5900], -28940.0, -3833.5, 0, 11.0);
refineSummary = localBuildSummary(true, -1.5084e6, 6.5590e4, [37.7810, 36.5910], -28939.8, -3832.0, 0, 18.0);
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

if selectedIdx ~= 2 || selectedTag ~= "periodic-fd-anchor"
  error('regressionMfUnknownSubsetAnchorGuard:UnexpectedSelection', ...
    'Expected the subset-anchor guard to keep periodic-fd-anchor, but got %s.', selectedTag);
end
if ~contains(reasonText, "subset-anchor guard")
  error('regressionMfUnknownSubsetAnchorGuard:UnexpectedReason', ...
    'Expected the selection reason to mention the subset-anchor guard, but got %s.', reasonText);
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-fd-anchor   : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(3, :)));
fprintf('PASS: regressionMfUnknownSubsetAnchorGuard\n');


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
