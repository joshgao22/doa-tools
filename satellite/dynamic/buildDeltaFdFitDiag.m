function fdDiag = buildDeltaFdFitDiag(truth)
%BUILDDELTAFDFITDIAG Summarize linear-fit residuals for Doppler truth.
% This helper builds the compact reference-Doppler and differential-Doppler
% tables used by the development scripts.

timeOffsetSec = reshape(truth.timeOffsetSec, 1, []);
numSat = size(truth.deltaFdSeries, 1);

refRmsResidualHz = sqrt(mean(truth.fdRefResidual .^ 2));
refMaxAbsResidualHz = max(abs(truth.fdRefResidual));
quadPhaseWinRad = localGetFieldOrDefault(truth, 'quadPhaseWin', NaN);
refTable = table(truth.fdRefFit, truth.fdRateFit, refRmsResidualHz, ...
  refMaxAbsResidualHz, truth.deltaFdWin, quadPhaseWinRad, 'VariableNames', ...
  {'fdRefFitHz', 'fdRateFitHzPerSec', 'rmsResidualHz', ...
   'maxAbsResidualHz', 'fdWinHz', 'quadPhaseWinRad'});

localSatIdx = (1:numSat).';
globalSatIdx = nan(numSat, 1);
if isfield(truth, 'selectedSatIdxGlobal') && ~isempty(truth.selectedSatIdxGlobal)
  globalSatIdx = reshape(truth.selectedSatIdxGlobal, [], 1);
  if numel(globalSatIdx) ~= numSat
    globalSatIdx = nan(numSat, 1);
  end
end

deltaFdRefHz = zeros(numSat, 1);
deltaFdRateHzPerSec = zeros(numSat, 1);
rmsResidualHz = zeros(numSat, 1);
maxAbsResidualHz = zeros(numSat, 1);
refFrameIdx = find(abs(timeOffsetSec) == min(abs(timeOffsetSec)), 1, 'first');
deltaFdRefFrameHz = truth.deltaFdSeries(:, refFrameIdx);

for iSat = 1:numSat
  currentSeries = reshape(truth.deltaFdSeries(iSat, :), 1, []);
  [fit0, fitRate] = localFitFdLine(timeOffsetSec, currentSeries);
  fitSeries = fit0 + fitRate * timeOffsetSec;
  residual = currentSeries - fitSeries;
  deltaFdRefHz(iSat) = fit0;
  deltaFdRateHzPerSec(iSat) = fitRate;
  rmsResidualHz(iSat) = sqrt(mean(residual .^ 2));
  maxAbsResidualHz(iSat) = max(abs(residual));
end

deltaTable = table(localSatIdx, globalSatIdx, deltaFdRefHz, deltaFdRateHzPerSec, ...
  deltaFdRefFrameHz, rmsResidualHz, maxAbsResidualHz, 'VariableNames', ...
  {'localSatIdx', 'globalSatIdx', 'deltaFdFitHz', 'deltaFdRateFitHzPerSec', ...
   'deltaFdRefFrameHz', 'rmsResidualHz', 'maxAbsResidualHz'});

fdDiag = struct();
fdDiag.refTable = refTable;
fdDiag.deltaTable = deltaTable;
end


function [fit0, fitRate] = localFitFdLine(timeOffsetSec, fdSeries)
%LOCALFITFDLINE Fit fd(t) = fit0 + fitRate * t by least squares.

timeOffsetSec = reshape(timeOffsetSec, [], 1);
fdSeries = reshape(fdSeries, [], 1);
regMat = [ones(numel(timeOffsetSec), 1), timeOffsetSec];
coef = regMat \ fdSeries;
fit0 = coef(1);
fitRate = coef(2);
end


function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one struct/object field with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end
