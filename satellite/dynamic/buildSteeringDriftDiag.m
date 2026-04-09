function steerDiag = buildSteeringDriftDiag(sceneSeq, refFrameIdx, satIdxGlobal)
%BUILDSTEERINGDRIFTDIAG Summarize frame-wise local-DoA drift.
% This helper keeps the dynamic development scripts light by moving the
% steering-drift summary into a reusable utility.
%
%Syntax:
%  steerDiag = buildSteeringDriftDiag(sceneSeq, refFrameIdx)
%  steerDiag = buildSteeringDriftDiag(sceneSeq, refFrameIdx, satIdxGlobal)
%
%Output:
%  steerDiag - table with one row per satellite.

if ~isfield(sceneSeq, 'localDoa') || isempty(sceneSeq.localDoa)
  steerDiag = table();
  return;
end

localDoa = sceneSeq.localDoa;
if ndims(localDoa) == 3
  localDoaUse = localDoa;
elseif ndims(localDoa) == 4
  localDoaUse = localDoa(:, :, 1, :);
else
  error('buildSteeringDriftDiag:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end

numSat = size(localDoaUse, 2);
numFrame = size(localDoaUse, 3);
refLocalDoa = localDoaUse(:, :, refFrameIdx);
maxAbsAzDriftDeg = zeros(numSat, 1);
maxAbsElDriftDeg = zeros(numSat, 1);
maxNormDriftDeg = zeros(numSat, 1);
rmsNormDriftDeg = zeros(numSat, 1);

for iSat = 1:numSat
  deltaRad = zeros(2, numFrame);
  for iFrame = 1:numFrame
    deltaRad(:, iFrame) = localWrapDoaDiff(localDoaUse(:, iSat, iFrame) - refLocalDoa(:, iSat));
  end
  deltaDeg = rad2deg(deltaRad);
  normDriftDeg = sqrt(sum(deltaDeg .^ 2, 1));
  maxAbsAzDriftDeg(iSat) = max(abs(deltaDeg(1, :)));
  maxAbsElDriftDeg(iSat) = max(abs(deltaDeg(2, :)));
  maxNormDriftDeg(iSat) = max(normDriftDeg);
  rmsNormDriftDeg(iSat) = sqrt(mean(normDriftDeg .^ 2));
end

localSatIdx = (1:numSat).';
if nargin < 3 || isempty(satIdxGlobal)
  satIdxGlobal = nan(numSat, 1);
else
  satIdxGlobal = reshape(satIdxGlobal, [], 1);
  if numel(satIdxGlobal) ~= numSat
    satIdxGlobal = nan(numSat, 1);
  end
end

steerDiag = table(localSatIdx, satIdxGlobal, maxAbsAzDriftDeg, maxAbsElDriftDeg, ...
  maxNormDriftDeg, rmsNormDriftDeg, 'VariableNames', ...
  {'localSatIdx', 'globalSatIdx', 'maxAbsAzDriftDeg', 'maxAbsElDriftDeg', ...
   'maxNormDriftDeg', 'rmsNormDriftDeg'});
end


function deltaRad = localWrapDoaDiff(deltaRad)
%LOCALWRAPDOADIFF Wrap azimuth/elevation deltas to principal values.

deltaRad = reshape(deltaRad, 2, []);
deltaRad(1, :) = angle(exp(1j * deltaRad(1, :)));
deltaRad(2, :) = angle(exp(1j * deltaRad(2, :)));
end
