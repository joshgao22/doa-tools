function [refCoherence, nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval)
%EXTRACTMFPROBECOHERENCEFLOOR Extract ref/non-ref coherence from final eval.

refCoherence = NaN;
nonRefCoherenceFloor = NaN;
if ~isstruct(finalEval)
  return;
end
coherenceSat = reshape(localGetFieldOrDefault(finalEval, 'coherenceSat', []), [], 1);
if isempty(coherenceSat)
  return;
end
refSatIdxLocal = localScalarOrDefault(localGetFieldOrDefault(finalEval, 'refSatIdxLocal', []), 1);
if isfinite(refSatIdxLocal) && refSatIdxLocal >= 1 && refSatIdxLocal <= numel(coherenceSat)
  refCoherence = coherenceSat(refSatIdxLocal);
  nonRefMask = true(size(coherenceSat));
  nonRefMask(refSatIdxLocal) = false;
  nonRefVals = coherenceSat(nonRefMask & isfinite(coherenceSat));
  if ~isempty(nonRefVals)
    nonRefCoherenceFloor = min(nonRefVals);
  end
else
  refCoherence = coherenceSat(1);
end
end

function value = localScalarOrDefault(value, defaultValue)
%LOCALSCALARORDEFAULT Return first scalar value or default.

if isempty(value)
  value = defaultValue;
  return;
end
value = double(value(1));
if ~isfinite(value)
  value = defaultValue;
end
end

function val = localGetFieldOrDefault(s, fieldName, defaultVal)
%LOCALGETFIELDORDEFAULT Get a field with a default fallback.

if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
