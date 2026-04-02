function stat = summarizeMonteCarloStat(sampleVal, failMask, statOpt)
%SUMMARIZEMONTECARLOSTAT Summarize Monte Carlo scalar error samples.
% Returns a compact statistic structure for Monte Carlo error sequences.
% The function is intended for unified post-processing of angle error,
% fdRef error, and other scalar metrics used in the paper experiments.
%
%Syntax:
%  stat = summarizeMonteCarloStat(sampleVal)
%  stat = summarizeMonteCarloStat(sampleVal, failMask)
%  stat = summarizeMonteCarloStat(sampleVal, failMask, statOpt)
%
%Inputs:
%  sampleVal         - scalar error samples. The summary is taken along one
%                      dimension specified by statOpt.dim.
%
%  failMask          - optional failure mask
%                      - []            : no additional failure mask
%                      - scalar logical: applied to all samples
%                      - same size as sampleVal
%
%  statOpt           - optional options structure
%    .dim            - summary dimension. Default: first non-singleton
%    .trimRatio      - two-sided trim ratio for trimmed RMSE. Default: 0.05
%    .winsorRatio    - two-sided winsor ratio. Default: trimRatio
%    .maxAbs         - samples with abs(error) larger than this value are
%                      treated as failures. Default: inf
%
%Output:
%  stat              - summary structure with fields
%    .dim
%    .numTrial
%    .numValid
%    .numFinite
%    .failCount
%    .failRate
%    .rmse
%    .meanAbs
%    .median
%    .p95
%    .maxAbs
%    .trimmedRmse
%    .winsorRmse
%    .trimRatio
%    .winsorRatio
%    .maxAbsThreshold
%
%Notes:
%  - All summary metrics are computed from abs(sampleVal).
%  - Non-finite samples are always counted as failures.
%  - failMask and maxAbs are merged into one unified valid-sample mask.
%  - For multi-dimensional input, output fields keep the input size with
%    the summary dimension collapsed to singleton.
%
%See also:
%  calcDoaAngleError, calcLatlonAngleError

arguments
  sampleVal {mustBeNumeric}
  failMask = []
  statOpt.dim = []
  statOpt.trimRatio (1,1) double = 0.05
  statOpt.winsorRatio = []
  statOpt.maxAbs (1,1) double = inf
end

localValidateRatio(statOpt.trimRatio, 'statOpt.trimRatio');
if isempty(statOpt.winsorRatio)
  statOpt.winsorRatio = statOpt.trimRatio;
end
localValidateRatio(statOpt.winsorRatio, 'statOpt.winsorRatio');
if ~isscalar(statOpt.maxAbs) || ~isnumeric(statOpt.maxAbs) || ...
    ~(isfinite(statOpt.maxAbs) || isinf(statOpt.maxAbs)) || statOpt.maxAbs < 0
  error('summarizeMonteCarloStat:InvalidMaxAbs', ...
    'statOpt.maxAbs must be a nonnegative scalar or inf.');
end

sumDim = statOpt.dim;
if isempty(sumDim)
  sumDim = localFirstNonSingletonDim(sampleVal);
end
localValidateDim(sumDim, ndims(sampleVal));

failMask = localExpandFailMask(failMask, size(sampleVal));
absVal = abs(sampleVal);

finiteMask = isfinite(absVal);
inRangeMask = absVal <= statOpt.maxAbs;
validMask = finiteMask & inRangeMask & ~failMask;

permOrder = localBuildPermOrder(ndims(absVal), sumDim);
absValPerm = permute(absVal, permOrder);
validMaskPerm = permute(validMask, permOrder);
finiteMaskPerm = permute(finiteMask, permOrder);

sizePerm = size(absValPerm);
numTrial = sizePerm(1);
numCase = prod(sizePerm(2:end));

absValMat = reshape(absValPerm, numTrial, numCase);
validMaskMat = reshape(validMaskPerm, numTrial, numCase);
finiteMaskMat = reshape(finiteMaskPerm, numTrial, numCase);

numValidVec = zeros(1, numCase);
numFiniteVec = zeros(1, numCase);
failCountVec = zeros(1, numCase);
rmseVec = nan(1, numCase);
meanAbsVec = nan(1, numCase);
medianVec = nan(1, numCase);
p95Vec = nan(1, numCase);
maxAbsVec = nan(1, numCase);
trimmedRmseVec = nan(1, numCase);
winsorRmseVec = nan(1, numCase);

for iCase = 1:numCase
  sampleCur = absValMat(:, iCase);
  validCur = validMaskMat(:, iCase);
  finiteCur = finiteMaskMat(:, iCase);

  numValidVec(iCase) = sum(validCur);
  numFiniteVec(iCase) = sum(finiteCur);
  failCountVec(iCase) = numTrial - numValidVec(iCase);

  if ~any(validCur)
    continue;
  end

  val = sampleCur(validCur);
  rmseVec(iCase) = sqrt(mean(val .^ 2));
  meanAbsVec(iCase) = mean(val);
  medianVec(iCase) = localPercentile(val, 50);
  p95Vec(iCase) = localPercentile(val, 95);
  maxAbsVec(iCase) = max(val);
  trimmedRmseVec(iCase) = localTrimmedRmse(val, statOpt.trimRatio);
  winsorRmseVec(iCase) = localWinsorRmse(val, statOpt.winsorRatio);
end

metricSizePerm = [1, sizePerm(2:end)];

stat = struct();
stat.dim = sumDim;
stat.numTrial = localRestoreMetric(repmat(numTrial, 1, numCase), metricSizePerm, permOrder);
stat.numValid = localRestoreMetric(numValidVec, metricSizePerm, permOrder);
stat.numFinite = localRestoreMetric(numFiniteVec, metricSizePerm, permOrder);
stat.failCount = localRestoreMetric(failCountVec, metricSizePerm, permOrder);
stat.failRate = stat.failCount ./ max(stat.numTrial, 1);
stat.rmse = localRestoreMetric(rmseVec, metricSizePerm, permOrder);
stat.meanAbs = localRestoreMetric(meanAbsVec, metricSizePerm, permOrder);
stat.median = localRestoreMetric(medianVec, metricSizePerm, permOrder);
stat.p95 = localRestoreMetric(p95Vec, metricSizePerm, permOrder);
stat.maxAbs = localRestoreMetric(maxAbsVec, metricSizePerm, permOrder);
stat.trimmedRmse = localRestoreMetric(trimmedRmseVec, metricSizePerm, permOrder);
stat.winsorRmse = localRestoreMetric(winsorRmseVec, metricSizePerm, permOrder);
stat.trimRatio = statOpt.trimRatio;
stat.winsorRatio = statOpt.winsorRatio;
stat.maxAbsThreshold = statOpt.maxAbs;
end

function failMask = localExpandFailMask(failMask, targetSize)
%LOCALEXPANDFAILMASK Expand optional failure mask to target size.

if isempty(failMask)
  failMask = false(targetSize);
  return;
end

if isscalar(failMask)
  failMask = repmat(logical(failMask), targetSize);
  return;
end

if ~islogical(failMask)
  failMask = logical(failMask);
end

if ~isequal(size(failMask), targetSize)
  error('summarizeMonteCarloStat:FailMaskSizeMismatch', ...
    'failMask must be empty, scalar, or have the same size as sampleVal.');
end
end

function dim = localFirstNonSingletonDim(x)
%LOCALFIRSTNONSINGLETONDIM Return first non-singleton dimension.

sz = size(x);
dim = find(sz > 1, 1, 'first');
if isempty(dim)
  dim = 1;
end
end

function permOrder = localBuildPermOrder(numDim, sumDim)
%LOCALBUILDPERMORDER Move summary dimension to the front.

permOrder = [sumDim, 1:sumDim-1, sumDim+1:numDim];
end

function localValidateDim(sumDim, numDim)
%LOCALVALIDATEDIM Validate summary dimension.

if ~isscalar(sumDim) || ~isnumeric(sumDim) || ~isfinite(sumDim) || ...
    sumDim < 1 || mod(sumDim, 1) ~= 0 || sumDim > numDim
  error('summarizeMonteCarloStat:InvalidDim', ...
    'statOpt.dim must be an integer within [1, ndims(sampleVal)].');
end
end

function localValidateRatio(value, nameStr)
%LOCALVALIDATERATIO Validate two-sided trimming or winsor ratio.

if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
    value < 0 || value >= 0.5
  error('summarizeMonteCarloStat:InvalidRatio', ...
    '%s must be a scalar within [0, 0.5).', nameStr);
end
end

function metric = localRestoreMetric(metricVec, metricSizePerm, permOrder)
%LOCALRESTOREMETRIC Restore metric array to the original dimension order.

metricPerm = reshape(metricVec, metricSizePerm);
metric = ipermute(metricPerm, permOrder);
end

function value = localPercentile(x, pct)
%LOCALPERCENTILE Compute percentile by linear interpolation.

x = sort(x(:));
numVal = numel(x);

if numVal == 0
  value = NaN;
  return;
end

if numVal == 1
  value = x;
  return;
end

pos = 1 + (numVal - 1) * pct / 100;
idxLo = floor(pos);
idxHi = ceil(pos);
frac = pos - idxLo;

value = (1 - frac) * x(idxLo) + frac * x(idxHi);
end

function value = localTrimmedRmse(x, trimRatio)
%LOCALTRIMMEDRMSE Compute two-sided trimmed RMSE.

x = sort(x(:));
numVal = numel(x);
numTrim = floor(trimRatio * numVal);

idxLo = numTrim + 1;
idxHi = numVal - numTrim;
if idxLo > idxHi
  value = NaN;
  return;
end

xTrim = x(idxLo:idxHi);
value = sqrt(mean(xTrim .^ 2));
end

function value = localWinsorRmse(x, winsorRatio)
%LOCALWINSORRMSE Compute two-sided winsorized RMSE.

x = sort(x(:));
numVal = numel(x);
if numVal == 0
  value = NaN;
  return;
end

numWinsor = floor(winsorRatio * numVal);
if numWinsor > 0
  idxLo = numWinsor + 1;
  idxHi = numVal - numWinsor;
  x(1:numWinsor) = x(idxLo);
  x(idxHi+1:end) = x(idxHi);
end

value = sqrt(mean(x .^ 2));
end
