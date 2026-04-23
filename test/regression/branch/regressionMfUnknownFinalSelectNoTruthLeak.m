% Regression check for objective-only final unknown selection.
% This script verifies one contract only:
%   1) final unknown candidate selection must not depend on truth-aware
%      angle/fdRef/fdRate errors;
%   2) the selected candidate must follow objective-only ranking;
%   3) the reported selection reason must point to the decisive objective
%      metric when candidates differ.
clear(); close all;

localAddProjectPath();

fprintf('Running regressionMfUnknownFinalSelectNoTruthLeak ...\n');

wideSummary = localBuildSummary(true,  -1.5082e6, 6.5616e4, 1.0e-3,   1.0e-2,  10.0, 0, 12.0);
refineSummary = localBuildSummary(true, -1.5086e6, 6.5576e4, 5.0e-2, 500.0, 5000.0, 0, 18.0);
[selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary( ...
  {wideSummary; refineSummary}, ["periodic-wide"; "periodic-in-tooth"]);

if selectedIdx ~= 2 || selectedTag ~= "periodic-in-tooth"
  error('regressionMfUnknownFinalSelectNoTruthLeak:UnexpectedSelection', ...
    'Expected the refine candidate to win by final objective, but got %s.', selectedTag);
end
if ~contains(reasonText, "final objective")
  error('regressionMfUnknownFinalSelectNoTruthLeak:UnexpectedReason', ...
    'Expected the selection reason to mention final objective, but got %s.', reasonText);
end

fprintf('  selected final idx         : %d\n', selectedIdx);
fprintf('  selected final tag         : %s\n', selectedTag);
fprintf('  final-select reason        : %s\n', reasonText);
fprintf('  score periodic-wide        : %s\n', localFormatNumericRow(scoreMat(1, :)));
fprintf('  score periodic-in-tooth    : %s\n', localFormatNumericRow(scoreMat(2, :)));
fprintf('PASS: regressionMfUnknownFinalSelectNoTruthLeak\n');


function summary = localBuildSummary(isResolved, finalObj, finalResidualNorm, angleErrDeg, fdRefErrHz, fdRateErrHzPerSec, toothIdx, runTimeMs)
summary = struct();
summary.isResolved = isResolved;
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.angleErrDeg = angleErrDeg;
summary.fdRefErrHz = fdRefErrHz;
summary.fdRateErrHzPerSec = fdRateErrHzPerSec;
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
