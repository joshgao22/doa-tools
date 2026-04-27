function [runPolish, gateReason, gateDiag] = shouldRunSimpleVerySmallDoaPolish(satMode, periodicDoaSeed, selectedReplaySummary, flowOpt)
%SHOULDRUNSIMPLEVERYSMALLDOAPOLISH Gate the simple-flow hard-case DoA polish.
% Keep the very-small DoA polish on a strict bypass. It should only run on
% same-tooth hard cases where fd has already stabilized but a small DoA
% disagreement still remains between trusted frozen seeds.

arguments
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  periodicDoaSeed (1,1) struct
  selectedReplaySummary (1,1) struct
  flowOpt (1,1) struct
end

runPolish = false;
gateReason = "disabled";
gateDiag = struct();
gateDiag.selectedReplaySeedSource = string(localGetFieldOrDefault(periodicDoaSeed, 'selectedReplaySeedSource', "unknown"));
gateDiag.subsetDriftFromStaticDeg = localGetFieldOrDefault(periodicDoaSeed, 'subsetDriftFromStaticDeg', NaN);
gateDiag.subsetRankMarginRelative = localGetFieldOrDefault(periodicDoaSeed, 'subsetRankMarginRelative', NaN);
gateDiag.frozenDoaDisagreementDeg = localGetFieldOrDefault(periodicDoaSeed, 'frozenDoaDisagreementDeg', NaN);
gateDiag.frozenRelativeObjGap = localGetFieldOrDefault(periodicDoaSeed, 'frozenRelativeObjGap', NaN);
gateDiag.selectedReplayHealthBucket = localBuildPeriodicHealthBucket(selectedReplaySummary);
gateDiag.selectedReplayToothIdx = localGetFieldOrDefault(selectedReplaySummary, 'toothIdx', NaN);
gateDiag.selectedReplayToothResidualHz = abs(localGetFieldOrDefault(selectedReplaySummary, 'toothResidualHz', inf));
if satMode ~= "multi"
  gateReason = "single-mode";
  return;
end
if ~localGetFieldOrDefault(flowOpt, 'periodicRefineEnableVerySmallDoaPolish', false)
  gateReason = "polish-disabled";
  return;
end
if string(gateDiag.selectedReplaySeedSource) ~= "selected-subset"
  gateReason = "replay-not-from-subset";
  return;
end
if ~logical(localGetFieldOrDefault(selectedReplaySummary, 'isResolved', false))
  gateReason = "selected-replay-unresolved";
  return;
end
if ~(isfinite(gateDiag.selectedReplayToothIdx) && abs(gateDiag.selectedReplayToothIdx) == 0)
  gateReason = "selected-replay-not-central-tooth";
  return;
end
maxToothResidualHz = localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMaxSelectedToothResidualHz', 50);
if ~(isfinite(gateDiag.selectedReplayToothResidualHz) && gateDiag.selectedReplayToothResidualHz <= maxToothResidualHz)
  gateReason = "selected-replay-tooth-residual-too-large";
  return;
end
if localBuildPeriodicHealthBucket(selectedReplaySummary) > localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMaxHealthBucket', 0)
  gateReason = "fd-health-not-ready";
  return;
end
subsetDriftDeg = gateDiag.subsetDriftFromStaticDeg;
if isfinite(subsetDriftDeg) && subsetDriftDeg > flowOpt.periodicRefinePolishTriggerDoaDriftDeg
  gateReason = "subset-static-drift-too-large";
  return;
end
subsetMargin = gateDiag.subsetRankMarginRelative;
if ~(isfinite(subsetMargin) && subsetMargin >= flowOpt.periodicRefineSubsetTrustMinRelativeMargin)
  gateReason = "subset-rank-margin-too-small";
  return;
end
frozenDoaGap = gateDiag.frozenDoaDisagreementDeg;
minFrozenDoaGap = localGetFieldOrDefault(flowOpt, 'periodicRefinePolishMinFrozenDoaDisagreementDeg', 5e-4);
if ~isfinite(frozenDoaGap)
  gateReason = "missing-frozen-doa-gap";
  return;
end
if frozenDoaGap < minFrozenDoaGap
  gateReason = "frozen-doa-gap-too-small";
  return;
end
if frozenDoaGap > flowOpt.periodicRefineMaxFrozenDoaDisagreementDeg
  gateReason = "frozen-doa-disagreement-too-large";
  return;
end
frozenObjGap = gateDiag.frozenRelativeObjGap;
if ~(isfinite(frozenObjGap) && frozenObjGap <= flowOpt.periodicRefinePolishMaxFrozenRelativeObjGap)
  gateReason = "frozen-objective-gap-too-large";
  return;
end
runPolish = true;
gateReason = "triggered";
end

function bucket = localBuildPeriodicHealthBucket(summary)
bucket = 0;
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefSupportRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefFitRatioFloor', 1) < 0.95);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefConsistencyRatioFloor', 1) < 0.995);
bucket = bucket + double(localGetFieldOrDefault(summary, 'nonRefCoherenceFloor', 1) < 0.995);
rmsPhaseResid = localGetFieldOrDefault(summary, 'nonRefRmsPhaseResidRad', NaN);
if isfinite(rmsPhaseResid)
  bucket = bucket + double(rmsPhaseResid > 0.003);
end
maxPhaseResid = localGetFieldOrDefault(summary, 'nonRefMaxAbsPhaseResidRad', NaN);
if isfinite(maxPhaseResid)
  bucket = bucket + double(maxPhaseResid > 0.005);
end
negativeRatio = localGetFieldOrDefault(summary, 'maxNonRefNegativeProjectionRatio', NaN);
if isfinite(negativeRatio)
  bucket = bucket + double(negativeRatio > 0.05);
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
