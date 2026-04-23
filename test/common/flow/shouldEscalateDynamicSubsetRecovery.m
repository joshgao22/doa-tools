function [tf, reasonText, diag] = shouldEscalateDynamicSubsetRecovery(summaryCell, selectedSummary, knownSummary, opt)
%SHOULDESCALATEDYNAMICSUBSETRECOVERY Decide whether fast subset recovery should escalate.
% The fast curated-only route is intentionally lightweight. This helper
% keeps it honest by triggering extra random or wide candidates only when the
% current subset winner looks suspicious under objective-only, truth-free
% diagnostics.

arguments
  summaryCell cell
  selectedSummary (1, 1) struct = struct()
  knownSummary (1, 1) struct = struct()
  opt (1, 1) struct = struct()
end

opt = localApplyDefaults(opt);
reasonText = "confident";
diag = struct();
diag.numCandidate = numel(summaryCell);
diag.numResolved = 0;
diag.uniqueResolvedTooth = [];
diag.selectedToothIdx = localGetFiniteField(selectedSummary, 'toothIdx', NaN);
diag.selectedToothResidualHz = abs(localGetFiniteField(selectedSummary, 'toothResidualHz', NaN));
diag.fdDriftFromKnownHz = localCalcFdDrift(selectedSummary, knownSummary);
diag.doaDriftFromKnownDeg = localCalcDoaDriftDeg(selectedSummary, knownSummary);

tf = false;
if ~logical(localGetLogicalField(opt, 'enableEscalation', true))
  reasonText = "disabled";
  return;
end

isResolved = localGetLogicalField(selectedSummary, 'isResolved', false);
hasFiniteScore = isfinite(localGetFiniteField(selectedSummary, 'finalObj', NaN)) && ...
  isfinite(localGetFiniteField(selectedSummary, 'finalResidualNorm', NaN));
if ~(isResolved && hasFiniteScore)
  tf = true;
  reasonText = "selected subset is not trusted";
  return;
end

resolvedTooth = nan(numel(summaryCell), 1);
numResolved = 0;
for iCase = 1:numel(summaryCell)
  summaryUse = summaryCell{iCase};
  if localGetLogicalField(summaryUse, 'isResolved', false)
    numResolved = numResolved + 1;
    resolvedTooth(numResolved) = localGetFiniteField(summaryUse, 'toothIdx', NaN);
  end
end
resolvedTooth = resolvedTooth(1:numResolved);
resolvedTooth = resolvedTooth(isfinite(resolvedTooth));
diag.numResolved = numel(resolvedTooth);
diag.uniqueResolvedTooth = unique(resolvedTooth);
if logical(localGetLogicalField(opt, 'escalateOnToothDisagreement', true)) && ...
    numel(diag.uniqueResolvedTooth) >= 2
  tf = true;
  reasonText = "resolved candidates disagree on tooth";
  return;
end

maxTrustedToothIdx = localGetFiniteField(opt, 'maxTrustedToothIdx', 0);
if isfinite(diag.selectedToothIdx) && abs(diag.selectedToothIdx) > maxTrustedToothIdx
  tf = true;
  reasonText = "selected tooth is not central";
  return;
end

maxTrustedToothResidualHz = localGetFiniteField(opt, 'maxTrustedToothResidualHz', inf);
if isfinite(diag.selectedToothResidualHz) && diag.selectedToothResidualHz > maxTrustedToothResidualHz
  tf = true;
  reasonText = "selected tooth residual is too large";
  return;
end

maxFdDriftFromKnownHz = localGetFiniteField(opt, 'maxFdDriftFromKnownHz', inf);
if isfinite(diag.fdDriftFromKnownHz) && diag.fdDriftFromKnownHz > maxFdDriftFromKnownHz
  tf = true;
  reasonText = "selected fdRef drifts too far from CP-K";
  return;
end

maxDoaDriftFromKnownDeg = localGetFiniteField(opt, 'maxDoaDriftFromKnownDeg', inf);
if isfinite(diag.doaDriftFromKnownDeg) && diag.doaDriftFromKnownDeg > maxDoaDriftFromKnownDeg
  tf = true;
  reasonText = "selected DoA drifts too far from CP-K";
  return;
end
end


function opt = localApplyDefaults(opt)
%LOCALAPPLYDEFAULTS Fill the escalation option struct.

if ~isfield(opt, 'enableEscalation') || isempty(opt.enableEscalation)
  opt.enableEscalation = true;
end
if ~isfield(opt, 'escalateOnToothDisagreement') || isempty(opt.escalateOnToothDisagreement)
  opt.escalateOnToothDisagreement = true;
end
if ~isfield(opt, 'maxTrustedToothIdx') || isempty(opt.maxTrustedToothIdx)
  opt.maxTrustedToothIdx = 0;
end
if ~isfield(opt, 'maxTrustedToothResidualHz') || isempty(opt.maxTrustedToothResidualHz)
  opt.maxTrustedToothResidualHz = 50;
end
if ~isfield(opt, 'maxFdDriftFromKnownHz') || isempty(opt.maxFdDriftFromKnownHz)
  opt.maxFdDriftFromKnownHz = 500;
end
if ~isfield(opt, 'maxDoaDriftFromKnownDeg') || isempty(opt.maxDoaDriftFromKnownDeg)
  opt.maxDoaDriftFromKnownDeg = 0.003;
end
end


function fdDriftHz = localCalcFdDrift(summaryA, summaryB)
%LOCALCALCFDDRIFT Compute one fdRef drift metric.

fdDriftHz = NaN;
fdA = localGetFiniteField(summaryA, 'fdRefEst', NaN);
fdB = localGetFiniteField(summaryB, 'fdRefEst', NaN);
if isfinite(fdA) && isfinite(fdB)
  fdDriftHz = abs(fdA - fdB);
end
end


function doaDriftDeg = localCalcDoaDriftDeg(summaryA, summaryB)
%LOCALCALCDOADRIFTDEG Compute one symmetric DoA drift metric.

doaDriftDeg = NaN;
doaA = localGetFieldOrDefault(summaryA, 'doaParamEst', []);
doaB = localGetFieldOrDefault(summaryB, 'doaParamEst', []);
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = NaN;
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field with default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end


function value = localGetFiniteField(dataStruct, fieldName, defaultValue)
%LOCALGETFINITEFIELD Read one scalar field with finite fallback.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
rawValue = rawValue(1);
if isfinite(rawValue)
  value = rawValue;
end
end


function value = localGetLogicalField(dataStruct, fieldName, defaultValue)
%LOCALGETLOGICALFIELD Read one logical-like scalar field.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
value = logical(rawValue(1));
end
