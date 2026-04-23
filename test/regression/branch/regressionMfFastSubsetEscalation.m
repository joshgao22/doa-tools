% Regression checks for fast curated->random/wide escalation heuristics.
clear(); close all;
localAddProjectPath();

fprintf('Running regressionMfFastSubsetEscalation ...\n');

knownSummary = localBuildSummary(true, 0, 0.02, [37.78; 36.59]);

curatedStable = { ...
  localBuildSummary(true, 0, 0.01, [37.7801; 36.5902]), ...
  localBuildSummary(true, 0, 0.02, [37.7800; 36.5901])};
[selectedStable, ~, ~] = selectBestDynamicSubsetSummary(curatedStable);
[needStable, reasonStable] = shouldEscalateDynamicSubsetRecovery( ...
  curatedStable, curatedStable{selectedStable}, knownSummary, struct());
if needStable
  error('regressionMfFastSubsetEscalation:FalsePositiveEscalation', ...
    'A central curated winner should not trigger fallback escalation.');
end

curatedDisagree = { ...
  localBuildSummary(true, 0, 0.01, [37.7801; 36.5902]), ...
  localBuildSummary(true, 4, 0.02, [37.7800; 36.5901])};
[selectedDisagree, ~, ~] = selectBestDynamicSubsetSummary(curatedDisagree);
[needDisagree, reasonDisagree] = shouldEscalateDynamicSubsetRecovery( ...
  curatedDisagree, curatedDisagree{selectedDisagree}, knownSummary, struct());
if ~needDisagree
  error('regressionMfFastSubsetEscalation:MissingToothDisagreementEscalation', ...
    'Resolved curated candidates that disagree on tooth should trigger fallback escalation.');
end

curatedDrift = { ...
  localBuildSummary(true, 0, 0.01, [37.79; 36.60], 1200), ...
  localBuildSummary(true, 0, 0.02, [37.7901; 36.6001], 1300)};
[selectedDrift, ~, ~] = selectBestDynamicSubsetSummary(curatedDrift);
[needDrift, reasonDrift] = shouldEscalateDynamicSubsetRecovery( ...
  curatedDrift, curatedDrift{selectedDrift}, knownSummary, struct());
if ~needDrift
  error('regressionMfFastSubsetEscalation:MissingKnownDriftEscalation', ...
    'A curated winner far from the CP-K anchor should trigger fallback escalation.');
end

fprintf('  stable reason   : %s\n', reasonStable);
fprintf('  disagree reason : %s\n', reasonDisagree);
fprintf('  drift reason    : %s\n', reasonDrift);
fprintf('PASS: regressionMfFastSubsetEscalation\n');


function summary = localBuildSummary(isResolved, toothIdx, toothResidualHz, doaParamEst, fdRefDriftHz)
%LOCALBUILDSUMMARY Build one synthetic subset summary.

if nargin < 5
  fdRefDriftHz = 0;
end
summary = struct();
summary.isResolved = logical(isResolved);
summary.toothIdx = toothIdx;
summary.toothResidualHz = toothResidualHz;
summary.doaParamEst = reshape(doaParamEst, 1, []);
summary.fdRefEst = -28940 + toothIdx * 750 + fdRefDriftHz;
summary.fdRateEst = -3833.5;
summary.finalObj = -1e6 - 10 * abs(toothIdx) - abs(fdRefDriftHz);
summary.finalResidualNorm = 3.2e4 + 10 * abs(toothIdx);
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
