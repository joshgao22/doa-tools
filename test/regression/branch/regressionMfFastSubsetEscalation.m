function regressionMfFastSubsetEscalation(varargin)
% Regression checks for fast curated->random/wide escalation gates.
% This script keeps the lightweight recovery-gate contracts in one place so
% random/wide rescue does not become a duplicate blanket branch elsewhere.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

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
if ~needDisagree || string(reasonDisagree) ~= "resolved candidates disagree on tooth"
  error('regressionMfFastSubsetEscalation:MissingToothDisagreementEscalation', ...
    'Resolved curated candidates that disagree on tooth should trigger fallback escalation.');
end

curatedDrift = { ...
  localBuildSummary(true, 0, 0.01, [37.79; 36.60], 1200), ...
  localBuildSummary(true, 0, 0.02, [37.7901; 36.6001], 1300)};
[selectedDrift, ~, ~] = selectBestDynamicSubsetSummary(curatedDrift);
[needDrift, reasonDrift] = shouldEscalateDynamicSubsetRecovery( ...
  curatedDrift, curatedDrift{selectedDrift}, knownSummary, struct());
if ~needDrift || string(reasonDrift) ~= "selected fdRef drifts too far from CP-K"
  error('regressionMfFastSubsetEscalation:MissingKnownDriftEscalation', ...
    'A curated winner far from the CP-K anchor should trigger fallback escalation.');
end

curatedUntrusted = {localBuildSummary(false, 0, 0.01, [37.7801; 36.5902])};
[needUntrusted, reasonUntrusted] = shouldEscalateDynamicSubsetRecovery( ...
  curatedUntrusted, curatedUntrusted{1}, knownSummary, struct());
if ~needUntrusted || string(reasonUntrusted) ~= "selected subset is not trusted"
  error('regressionMfFastSubsetEscalation:MissingUntrustedEscalation', ...
    'An unresolved selected subset must trigger fallback escalation.');
end

curatedOffTooth = {localBuildSummary(true, -1, 0.01, [37.7801; 36.5902])};
[needOffTooth, reasonOffTooth] = shouldEscalateDynamicSubsetRecovery( ...
  curatedOffTooth, curatedOffTooth{1}, knownSummary, struct());
if ~needOffTooth || string(reasonOffTooth) ~= "selected tooth is not central"
  error('regressionMfFastSubsetEscalation:MissingOffToothEscalation', ...
    'A non-central selected tooth must trigger fallback escalation.');
end

[needDisabled, reasonDisabled] = shouldEscalateDynamicSubsetRecovery( ...
  curatedOffTooth, curatedOffTooth{1}, knownSummary, struct('enableEscalation', false));
if needDisabled || string(reasonDisabled) ~= "disabled"
  error('regressionMfFastSubsetEscalation:DisabledEscalated', ...
    'The disabled escalation path must never request random/wide rescue.');
end

fprintf('  stable reason     : %s\n', reasonStable);
fprintf('  disagree reason   : %s\n', reasonDisagree);
fprintf('  drift reason      : %s\n', reasonDrift);
fprintf('  untrusted reason  : %s\n', reasonUntrusted);
fprintf('  off-tooth reason  : %s\n', reasonOffTooth);
fprintf('  disabled reason   : %s\n', reasonDisabled);
fprintf('PASS: regressionMfFastSubsetEscalation\n');
end

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
