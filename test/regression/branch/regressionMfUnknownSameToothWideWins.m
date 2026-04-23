% Regression check for one same-tooth wide escape.
% This script verifies one contract only:
%   1) when periodic-wide stays on the same tooth as the subset anchor;
%   2) and wide wins by objective-only ranking;
%   3) the subset-anchor guard must not force the final pick back to anchor.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownSameToothWideWins ...\n');

wideSummary = localBuildSummary(true, -1.5083e6, 6.5607e4, [37.7790, 36.5910], -28940.01, -3833.5, 0, 12.0);
anchorSummary = localBuildSummary(true, -1.5080e6, 6.5643e4, [37.7890, 36.5990], -28940.01, -3837.2, 0, 11.0);
refineSummary = localBuildSummary(true, -1.5080e6, 6.5640e4, [37.7890, 36.5990], -28940.03, -3823.3, 0, 18.0);
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
  error('regressionMfUnknownSameToothWideWins:UnexpectedSelection', ...
    'Expected periodic-wide to win when it stays on the anchor tooth, but got %s.', selectedTag);
end
if contains(reasonText, "subset-anchor guard")
  error('regressionMfUnknownSameToothWideWins:UnexpectedGuardReason', ...
    'Same-tooth wide escape should not be blocked by the subset-anchor guard.');
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-fd-anchor   : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(3, :)));
fprintf('PASS: regressionMfUnknownSameToothWideWins\n');


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
