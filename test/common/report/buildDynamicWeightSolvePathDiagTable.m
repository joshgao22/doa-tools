function pathDiagTable = buildDynamicWeightSolvePathDiagTable(caseMsDoa, caseStaticRefOnly, weightCase, weightAlpha, truth)
%BUILDDYNAMICWEIGHTSOLVEPATHDIAGTABLE Build one compact static weight solve-path table.

arguments
  caseMsDoa (1, 1) struct
  caseStaticRefOnly (1, 1) struct
  weightCase (1, :) struct
  weightAlpha (:, 1) double
  truth (1, 1) struct
end

numWeight = numel(weightAlpha);
angleErrDeg = nan(numWeight, 1);
fdRefEstHz = nan(numWeight, 1);
fdRefShiftVsW0Hz = nan(numWeight, 1);
angleToMsDoaAnchorDeg = nan(numWeight, 1);
angleToW0Deg = nan(numWeight, 1);
funcCount = nan(numWeight, 1);
iterations = nan(numWeight, 1);

if numWeight == 0
  pathDiagTable = table();
  return;
end

msDoaLatlon = getDoaDopplerLatlonEst(caseMsDoa.estResult);
w0Latlon = getDoaDopplerLatlonEst(weightCase(1).estResult);
fdRefW0 = getDoaDopplerFieldOrDefault(weightCase(1).estResult, 'fdRefEst', NaN);

for iWeight = 1:numWeight
  estResult = weightCase(iWeight).estResult;
  latlonEst = getDoaDopplerLatlonEst(estResult);
  angleErrDeg(iWeight) = calcLatlonAngleError(latlonEst, truth.latlonTrueDeg(:));
  fdRefEstHz(iWeight) = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
  fdRefShiftVsW0Hz(iWeight) = fdRefEstHz(iWeight) - fdRefW0;
  angleToMsDoaAnchorDeg(iWeight) = calcLatlonAngleError(latlonEst, msDoaLatlon);
  angleToW0Deg(iWeight) = calcLatlonAngleError(latlonEst, w0Latlon);
  optimInfo = getDoaDopplerFieldOrDefault(estResult, 'optimInfo', struct());
  funcCount(iWeight) = getDoaDopplerFieldOrDefault(optimInfo, 'funcCount', NaN);
  iterations(iWeight) = getDoaDopplerFieldOrDefault(optimInfo, 'iterations', NaN);
end

pathDiagTable = table(weightAlpha(:), angleErrDeg, fdRefEstHz, fdRefShiftVsW0Hz, ...
  angleToMsDoaAnchorDeg, angleToW0Deg, funcCount, iterations, ...
  'VariableNames', {'alphaSat2', 'angleErrDeg', 'fdRefEstHz', 'fdRefShiftVsW0Hz', ...
  'angleToMsDoaAnchorDeg', 'angleToW0Deg', 'funcCount', 'iterations'});
end
