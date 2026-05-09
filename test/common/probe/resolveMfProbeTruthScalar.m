function scalarValue = resolveMfProbeTruthScalar(truth, fieldNameList)
%RESOLVEMFPROBETRUTHSCALAR Resolve one scalar truth field from aliases.

scalarValue = NaN;
for iField = 1:numel(fieldNameList)
  fieldName = char(fieldNameList{iField});
  if isstruct(truth) && isfield(truth, fieldName)
    value = truth.(fieldName);
    if ~isempty(value)
      scalarValue = double(value(1));
      return;
    end
  end
end
end
