function [doaInitParam, seedSource, doaDriftDeg] = resolveSimplePeriodicRefineDoaSeed(flowOpt, staticSeedCase, selectedSubsetCase)
%RESOLVESIMPLEPERIODICREFINEDOASEED Compatibility wrapper for one replay seed.
% This wrapper preserves the old single-seed interface for callers that do
% not yet use the multi-candidate periodic replay path. The current simple
% flow should prefer resolveSimplePeriodicRefineDoaSeedCandidates and let
% the flow-level replay arbitration decide whether to trust the subset seed
% or recenter to the static seed.

arguments
  flowOpt (1,1) struct
  staticSeedCase (1,1) struct
  selectedSubsetCase (1,1) struct
end

seedCandidateList = resolveSimplePeriodicRefineDoaSeedCandidates( ...
  flowOpt, "multi", staticSeedCase, selectedSubsetCase);
if isempty(seedCandidateList)
  doaInitParam = nan(2, 1);
  seedSource = "missing";
  doaDriftDeg = NaN;
  return;
end

doAEntry = seedCandidateList(1);
doaInitParam = reshape(doAEntry.doaInitParam, [], 1);
seedSource = string(doAEntry.seedSource);
doaDriftDeg = doAEntry.subsetDriftFromStaticDeg;
end
