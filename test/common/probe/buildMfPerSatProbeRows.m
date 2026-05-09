function perSatRows = buildMfPerSatProbeRows(method, caseUse, truth, task, fixture)
%BUILDMFPERSATPROBEROWS Build per-satellite final-eval diagnostics for one MF case.

estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
coherenceSat = reshape(localGetFieldOrDefault(finalEval, 'coherenceSat', []), [], 1);
objectiveSat = reshape(localGetFieldOrDefault(finalEval, 'objectiveSat', []), [], 1);
residualSat = reshape(localGetFieldOrDefault(finalEval, 'residualSat', []), [], 1);
fitValueSat = reshape(localGetFieldOrDefault(finalEval, 'fitValueSat', []), [], 1);
localDoaArrUsed = localGetFieldOrDefault(finalEval, 'localDoaArrUsed', []);
truthLocalDoa = localExtractSceneLocalDoa(fixture.sceneSeq);
numSat = max([numel(coherenceSat), numel(objectiveSat), numel(residualSat), size(truthLocalDoa, 2), 1]);
refSatIdxLocal = localScalarOrDefault(localGetFieldOrDefault(finalEval, 'refSatIdxLocal', []), ...
  localScalarOrDefault(localGetFieldOrDefault(fixture.sceneSeq, 'refSatIdxLocal', []), 1));
perSatRows = repmat(emptyMfPerSatProbeRow(), numSat, 1);
for iSat = 1:numSat
  row = emptyMfPerSatProbeRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.satIdxLocal = iSat;
  row.isRefSat = isfinite(refSatIdxLocal) && iSat == refSatIdxLocal;
  row.coherence = localVectorValue(coherenceSat, iSat);
  row.objectiveSat = localVectorValue(objectiveSat, iSat);
  row.residualSat = localVectorValue(residualSat, iSat);
  row.fitValueSat = localVectorValue(fitValueSat, iSat);
  row.localDirErrMedianDeg = localLocalDoaMedianErr(localDoaArrUsed, truthLocalDoa, iSat);
  row.localDirErrMaxDeg = localLocalDoaMaxErr(localDoaArrUsed, truthLocalDoa, iSat);
  row.globalAngleErrDeg = localProbeAngleErr(localGetFieldOrDefault(estResult, 'doaParamEst', []), truth);
  row.fdRefErrHz = localGetFieldOrDefault(estResult, 'fdRefEst', NaN) - resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
  row.fdRateErrHzPerSec = localGetFieldOrDefault(estResult, 'fdRateEst', NaN) - resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
  perSatRows(iSat) = row;
end
end

function angleErrDeg = localProbeAngleErr(optVar, truth)
%LOCALPROBEANGLEERR Compute spherical angle error for an optimization vector.

angleErrDeg = NaN;
optVar = reshape(double(optVar), [], 1);
truthDoa = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
if numel(optVar) >= 2 && numel(truthDoa) >= 2
  angleErrDeg = calcLatlonAngleError(optVar(1:2), truthDoa(1:2));
end
end

function value = localScalarOrDefault(valueIn, defaultValue)
%LOCALSCALARORDEFAULT Return a finite scalar or default.

value = defaultValue;
valueIn = reshape(double(valueIn), [], 1);
if ~isempty(valueIn) && isfinite(valueIn(1))
  value = valueIn(1);
end
end

function value = localVectorValue(valueVec, idx)
%LOCALVECTORVALUE Safe vector index.

value = NaN;
if idx >= 1 && idx <= numel(valueVec)
  value = valueVec(idx);
end
end

function localDoa = localExtractSceneLocalDoa(sceneSeq)
%LOCALEXTRACTSCENELOCALDOA Extract truth local DoA as 2 x Ns x Nf.

localDoa = [];
if ~(isstruct(sceneSeq) && isfield(sceneSeq, 'localDoa') && ~isempty(sceneSeq.localDoa))
  return;
end
localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 4
  localDoa = squeeze(localDoaRaw(:, :, 1, :));
elseif ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
elseif ndims(localDoaRaw) == 2
  localDoa = reshape(localDoaRaw, size(localDoaRaw, 1), size(localDoaRaw, 2), 1);
end
end

function errMedian = localLocalDoaMedianErr(localDoaArrUsed, truthLocalDoa, satIdx)
%LOCALLOCALDOAMEDIANERR Median per-frame local direction error.

errVec = localLocalDoaErrVec(localDoaArrUsed, truthLocalDoa, satIdx);
if isempty(errVec)
  errMedian = NaN;
else
  errMedian = median(errVec, 'omitnan');
end
end

function errMax = localLocalDoaMaxErr(localDoaArrUsed, truthLocalDoa, satIdx)
%LOCALLOCALDOAMAXERR Max per-frame local direction error.

errVec = localLocalDoaErrVec(localDoaArrUsed, truthLocalDoa, satIdx);
if isempty(errVec)
  errMax = NaN;
else
  errMax = max(errVec, [], 'omitnan');
end
end

function errVec = localLocalDoaErrVec(localDoaArrUsed, truthLocalDoa, satIdx)
%LOCALLOCALDOAERRVEC Per-frame local direction angular error.

errVec = NaN(0, 1);
if isempty(localDoaArrUsed) || isempty(truthLocalDoa) || satIdx > size(truthLocalDoa, 2)
  return;
end
if ndims(localDoaArrUsed) == 2
  localDoaArrUsed = reshape(localDoaArrUsed, size(localDoaArrUsed, 1), size(localDoaArrUsed, 2), 1);
end
numFrame = min(size(localDoaArrUsed, 3), size(truthLocalDoa, 3));
if satIdx > size(localDoaArrUsed, 2) || numFrame <= 0
  return;
end
errVec = NaN(numFrame, 1);
for iFrame = 1:numFrame
  estDoa = reshape(localDoaArrUsed(:, satIdx, iFrame), [], 1);
  truthDoa = reshape(truthLocalDoa(:, satIdx, iFrame), [], 1);
  if numel(estDoa) >= 2 && numel(truthDoa) >= 2 && all(isfinite(estDoa(1:2))) && all(isfinite(truthDoa(1:2)))
    errVec(iFrame) = localAnglePairErrDeg(estDoa(1:2), truthDoa(1:2));
  end
end
end

function errDeg = localAnglePairErrDeg(estDoaDeg, truthDoaDeg)
%LOCALANGLEPAIRERRDEG Compute a compact two-angle local-direction error.

estDoaDeg = reshape(double(estDoaDeg), [], 1);
truthDoaDeg = reshape(double(truthDoaDeg), [], 1);
if numel(estDoaDeg) < 2 || numel(truthDoaDeg) < 2
  errDeg = NaN;
  return;
end
azErr = mod(estDoaDeg(1) - truthDoaDeg(1) + 180, 360) - 180;
elErr = estDoaDeg(2) - truthDoaDeg(2);
errDeg = hypot(azErr, elErr);
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
