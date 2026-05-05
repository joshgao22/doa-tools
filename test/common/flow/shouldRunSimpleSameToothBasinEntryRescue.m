function [runRescue, gateReason, gateDiag] = shouldRunSimpleSameToothBasinEntryRescue(satMode, selectedDecisionSummary, subsetTrustDiag, flowOpt)
%SHOULDRUNSIMPLESAMETOOTHBASINENTRYRESCUE Gate same-tooth basin-entry candidates.
% The gate is intentionally truth-free. It only uses the selected periodic
% decision summary and subset trust diagnostics produced by the flow.

arguments
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  selectedDecisionSummary (1,1) struct
  subsetTrustDiag (1,1) struct
  flowOpt (1,1) struct
end

runRescue = false;
gateReason = "disabled";

gateDiag = struct();
gateDiag.enable = logical(localGetFieldOrDefault(flowOpt, 'sameToothRescueEnable', false));
gateDiag.enableWhenMulti = logical(localGetFieldOrDefault(flowOpt, 'sameToothRescueEnableWhenMulti', true));
gateDiag.subsetTrusted = logical(localGetFieldOrDefault(subsetTrustDiag, 'isTrusted', false));
gateDiag.nonRefCoherenceFloor = localGetFieldOrDefault(selectedDecisionSummary, 'nonRefCoherenceFloor', NaN);
gateDiag.nonRefRmsPhaseResidRad = localGetFieldOrDefault(selectedDecisionSummary, 'nonRefRmsPhaseResidRad', NaN);
gateDiag.coherenceThreshold = localGetFieldOrDefault(flowOpt, 'sameToothRescueCoherenceThreshold', 0.20);
gateDiag.phaseResidThresholdRad = localGetFieldOrDefault(flowOpt, 'sameToothRescuePhaseResidThresholdRad', 1.00);
gateDiag.coherenceCollapsed = isfinite(gateDiag.nonRefCoherenceFloor) && ...
  (gateDiag.nonRefCoherenceFloor < gateDiag.coherenceThreshold);
gateDiag.phaseResidualLarge = isfinite(gateDiag.nonRefRmsPhaseResidRad) && ...
  (gateDiag.nonRefRmsPhaseResidRad >= gateDiag.phaseResidThresholdRad);

if ~gateDiag.enable
  return;
end
if satMode ~= "multi" || ~gateDiag.enableWhenMulti
  gateReason = "not-multi";
  return;
end
if ~gateDiag.subsetTrusted
  gateReason = "subset-not-trusted";
  return;
end
if gateDiag.coherenceCollapsed
  runRescue = true;
  gateReason = "non-ref-coherence-collapse";
  return;
end
if gateDiag.phaseResidualLarge
  runRescue = true;
  gateReason = "non-ref-phase-residual-large";
  return;
end

gateReason = "coherence-not-collapsed";
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
