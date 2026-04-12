function info = summarizeDoaDopplerCase(caseInfo, truth)
%SUMMARIZEDOADOPPLERCASE Summarize one estimate against truth.

estResult = getDoaDopplerFieldOrDefault(caseInfo, 'estResult', struct());
latlonEst = getDoaDopplerLatlonEst(estResult);
[angleErrDeg, latErrDeg, lonErrDeg] = calcLatlonAngleError(latlonEst, truth.latlonTrueDeg);

truthFdRef = localResolveTruthFdRef(truth);
truthFdRate = localResolveTruthFdRate(truth);

fdRefEst = getDoaDopplerFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRateEst = getDoaDopplerFieldOrDefault(estResult, 'fdRateEst', NaN);
fdLineRmseHz = localBuildFdLineRmse(fdRefEst, fdRateEst, truth);

fdSatRmseHz = NaN;
deltaFdRmseHz = NaN;
aux = getDoaDopplerFieldOrDefault(estResult, 'aux', struct());
fdSatEst = getDoaDopplerFieldOrDefault(aux, 'fdSatEst', []);
deltaFdEst = getDoaDopplerFieldOrDefault(aux, 'deltaFdRefEst', []);

truthFdSatSeries = getDoaDopplerFieldOrDefault(truth, 'fdSatSeries', []);
truthFdSatVec = getDoaDopplerFieldOrDefault(truth, 'fdSatTrueHz', []);
if ~isempty(fdSatEst)
  if isequal(size(fdSatEst), size(truthFdSatSeries))
    fdSatRmseHz = sqrt(mean((fdSatEst(:) - truthFdSatSeries(:)).^2));
  elseif numel(fdSatEst) == numel(truthFdSatVec)
    fdSatRmseHz = sqrt(mean((fdSatEst(:) - truthFdSatVec(:)).^2));
  end
end

truthDeltaFd = getDoaDopplerFieldOrDefault(truth, 'deltaFdTrueHz', []);
if ~isempty(deltaFdEst) && numel(deltaFdEst) == numel(truthDeltaFd)
  deltaFdRmseHz = sqrt(mean((deltaFdEst(:) - truthDeltaFd(:)).^2));
end

info = struct();
info.displayName = getDoaDopplerFieldOrDefault(caseInfo, 'displayName', string.empty(1, 0));
info.satMode = getDoaDopplerFieldOrDefault(caseInfo, 'satMode', string.empty(1, 0));
info.frameMode = getDoaDopplerFieldOrDefault(caseInfo, 'frameMode', string.empty(1, 0));
info.paramMode = getDoaDopplerFieldOrDefault(caseInfo, 'paramMode', string.empty(1, 0));
info.dynamicMode = getDoaDopplerFieldOrDefault(caseInfo, 'dynamicMode', string.empty(1, 0));
info.latEstDeg = latlonEst(1);
info.lonEstDeg = latlonEst(2);
info.latErrDeg = latErrDeg;
info.lonErrDeg = lonErrDeg;
info.angleErrDeg = angleErrDeg;
info.fdRefEstHz = fdRefEst;
info.fdRefErrHz = fdRefEst - truthFdRef;
info.fdRateEstHzPerSec = fdRateEst;
info.fdRateErrHzPerSec = fdRateEst - truthFdRate;
info.fdLineRmseHz = fdLineRmseHz;
info.fdSatRmseHz = fdSatRmseHz;
info.deltaFdRmseHz = deltaFdRmseHz;
info.exitflag = getDoaDopplerFieldOrDefault(estResult, 'exitflag', NaN);
info.isResolved = logical(getDoaDopplerFieldOrDefault(estResult, 'isResolved', false));
info.runTimeMs = 1e3 * getDoaDopplerOptimField(estResult, 'runTimeSec', NaN);
info.funcCount = getDoaDopplerOptimField(estResult, 'funcCount', NaN);
info.iterations = getDoaDopplerOptimField(estResult, 'iterations', NaN);
end

function fdLineRmseHz = localBuildFdLineRmse(fdRefEst, fdRateEst, truth)
%LOCALBUILDFDLINERMSE Build the reference-Doppler line RMSE when available.

fdLineRmseHz = NaN;
truthTime = getDoaDopplerFieldOrDefault(truth, 'timeOffsetSec', []);
truthFdLine = getDoaDopplerFieldOrDefault(truth, 'fdRefSeries', []);
if isempty(truthTime) || isempty(truthFdLine) || ~isfinite(fdRefEst)
  return;
end

if ~isfinite(fdRateEst)
  fdRateEst = 0;
end

fdLineEst = fdRefEst + fdRateEst * truthTime;
if numel(fdLineEst) == numel(truthFdLine)
  fdLineRmseHz = sqrt(mean((fdLineEst(:) - truthFdLine(:)).^2));
end
end


function truthFdRef = localResolveTruthFdRef(truth)
%LOCALRESOLVETRUTHFDREF Resolve reference-Doppler truth with subset awareness.
% Prefer subset-specific reference metadata when present so single-satellite
% ablation cases do not accidentally compare against the full-scene ref.

truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefTrueHz', NaN);

refSatIdxGlobal = getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN);
truthSatIdx = reshape(getDoaDopplerFieldOrDefault(truth, 'selectedSatIdxGlobal', []), [], 1);
truthFdSat = reshape(getDoaDopplerFieldOrDefault(truth, 'fdSatTrueHz', []), [], 1);
if isfinite(refSatIdxGlobal) && ~isempty(truthSatIdx) && ~isempty(truthFdSat)
  matchIdx = find(truthSatIdx == refSatIdxGlobal, 1, 'first');
  if ~isempty(matchIdx) && numel(truthFdSat) >= matchIdx && isfinite(truthFdSat(matchIdx))
    truthFdRef = truthFdSat(matchIdx);
    return;
  end
end

if ~isfinite(truthFdRef)
  truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefFit', NaN);
end
end


function truthFdRate = localResolveTruthFdRate(truth)
%LOCALRESOLVETRUTHFDRATE Resolve reference-rate truth with fallback.

truthFdRate = getDoaDopplerFieldOrDefault(truth, 'fdRateTrueHzPerSec', NaN);
if ~isfinite(truthFdRate)
  truthFdRate = getDoaDopplerFieldOrDefault(truth, 'fdRateFit', NaN);
end
end
