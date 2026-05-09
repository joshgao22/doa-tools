function twoValue = takeFirstTwoValues(value)
%TAKEFIRSTTWOVALUES Return a numeric two-vector or NaNs.

twoValue = NaN(2, 1);
value = reshape(double(value), [], 1);
if numel(value) >= 2
  twoValue = value(1:2);
end
end
