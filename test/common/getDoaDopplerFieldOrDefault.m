function fieldValue = getDoaDopplerFieldOrDefault(dataStruct, fieldName, defaultValue)
%GETDOADOPPLERFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end
