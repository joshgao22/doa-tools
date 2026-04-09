function value = getDoaDopplerOptimField(estResult, fieldName, defaultValue)
%GETDOADOPPLEROPTIMFIELD Read one field from estResult.optimInfo.

value = defaultValue;
optimInfo = getDoaDopplerFieldOrDefault(estResult, 'optimInfo', struct());
if isstruct(optimInfo) && isfield(optimInfo, fieldName)
  value = optimInfo.(fieldName);
end
end
