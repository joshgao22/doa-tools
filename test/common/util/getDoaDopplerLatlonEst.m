function latlonEst = getDoaDopplerLatlonEst(estResult)
%GETDOADOPPLERLATLONEST Extract a lat-lon estimate from estimator output.

latlonEst = [NaN; NaN];
if isempty(estResult)
  return;
end

aux = getDoaDopplerFieldOrDefault(estResult, 'aux', struct());
latlonAux = getDoaDopplerFieldOrDefault(aux, 'latlonEst', []);
if isnumeric(latlonAux) && numel(latlonAux) >= 2
  latlonEst = reshape(latlonAux(1:2), [], 1);
  return;
end

latlonDirect = getDoaDopplerFieldOrDefault(estResult, 'latlonEst', []);
if isnumeric(latlonDirect) && numel(latlonDirect) >= 2
  latlonEst = reshape(latlonDirect(1:2), [], 1);
  return;
end

doaType = getDoaDopplerFieldOrDefault(estResult, 'doaType', '');
doaParamEst = getDoaDopplerFieldOrDefault(estResult, 'doaParamEst', []);
if (ischar(doaType) || isstring(doaType)) && ...
    strcmpi(string(doaType), "latlon") && isnumeric(doaParamEst) && ...
    numel(doaParamEst) >= 2
  latlonEst = reshape(doaParamEst(1:2), [], 1);
end
end
