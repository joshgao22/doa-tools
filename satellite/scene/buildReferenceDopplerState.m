function dopplerState = buildReferenceDopplerState(scene, satPosEci, satVelEci, ...
  userPosEci, userVelEci, wavelength, timeOffsetSec, refSatIdxLocalOverride)
%BUILDREFERENCEDOPPLERSTATE Build reference-satellite Doppler states.
% This helper fixes the reference-satellite Doppler semantics shared by the
% static and dynamic estimators. Given one selected reference satellite, it
% computes the exact geometric Doppler of every satellite link and packs the
% reference-relative state
%
%   deltaFd = fdSat - fdRef
%
% together with optional line-fit summaries across a multi-frame window.
%
%Syntax:
%  dopplerState = buildReferenceDopplerState(scene, satPosEci, satVelEci, ...
%    userPosEci, userVelEci, wavelength)
%  dopplerState = buildReferenceDopplerState(scene, satPosEci, satVelEci, ...
%    userPosEci, userVelEci, wavelength, timeOffsetSec, refSatIdxLocalOverride)
%
%Inputs:
%  scene          - scene struct used to resolve the reference satellite.
%                   For multi-frame mode this should normally be the
%                   reference-frame scene.
%  satPosEci      - 3xNs or 3xNsxNf satellite ECI positions
%  satVelEci      - 3xNs or 3xNsxNf satellite ECI velocities
%  userPosEci     - 3xK user ECI positions
%  userVelEci     - 3xK user ECI velocities
%  wavelength     - carrier wavelength in meters
%  timeOffsetSec  - (optional) 1xK / Kx1 frame offsets in seconds. When
%                   provided, the helper also returns line-fit summaries for
%                   fdRef and deltaFd.
%  refSatIdxLocalOverride - (optional) explicit local reference-satellite
%                   index. When provided, the helper bypasses scene-based
%                   reference resolution and locks all Doppler states to the
%                   estimator-selected reference satellite.
%
%Output:
%  dopplerState   - struct with fields:
%                     .refState
%                     .refSatIdxLocal
%                     .refSatIdxGlobal
%                     .fdRef
%                     .fdSat
%                     .deltaFd
%                     .fdRefFit
%                     .fdRateFit
%                     .fdRefFitSeries
%                     .fdRefResidual
%                     .deltaFdRef
%                     .deltaFdRate
%                     .timeOffsetSec
%
%Invariant checks:
%  The helper enforces the project-level reference semantics directly:
%    1) deltaFd(refSat,:) = 0
%    2) fdSat(refSat,:)   = fdRef
%    3) fdSat             = fdRef + deltaFd
%
%See also:
%  resolveReferenceSatState, computeRelativeDopplerGeom

arguments
  scene
  satPosEci {mustBeNumeric, mustBeFinite}
  satVelEci {mustBeNumeric, mustBeFinite}
  userPosEci (3,:) {mustBeNumeric, mustBeFinite}
  userVelEci (3,:) {mustBeNumeric, mustBeFinite}
  wavelength (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  timeOffsetSec = []
  refSatIdxLocalOverride = []
end

if ~isempty(refSatIdxLocalOverride)
  refSatIdxLocal = localValidateReferenceSatIdx(refSatIdxLocalOverride, satPosEci);
  refState = localBuildReferenceStateFromIdx(scene, satPosEci, satVelEci, refSatIdxLocal, "override");
else
  [refSatPosResolve, refSatVelResolve] = localGetReferenceResolveState(scene, satPosEci, satVelEci);
  [refState, refSatIdxLocal] = resolveReferenceSatState(scene, refSatPosResolve, refSatVelResolve);
end
[deltaFd, satFd, refFd] = computeRelativeDopplerGeom( ...
  satPosEci, satVelEci, userPosEci, userVelEci, wavelength, refSatIdxLocal);

localCheckReferenceInvariant(deltaFd, satFd, refFd, refSatIdxLocal);

numState = size(deltaFd, 2);
if isempty(timeOffsetSec)
  timeOffsetSecUse = [];
else
  timeOffsetSecUse = reshape(timeOffsetSec, 1, []);
  if numel(timeOffsetSecUse) ~= numState
    error('buildReferenceDopplerState:TimeOffsetSizeMismatch', ...
      'numel(timeOffsetSec) must match the number of Doppler states.');
  end
end

if isempty(timeOffsetSecUse)
  fdRefFit = refFd(1);
  fdRateFit = 0;
  fdRefFitSeries = refFd;
  fdRefResidual = refFd - fdRefFitSeries;
  deltaFdRef = deltaFd(:, 1);
  deltaFdRate = zeros(size(deltaFd, 1), 1);
else
  [fdRefFit, fdRateFit] = localFitFdLine(timeOffsetSecUse, refFd);
  fdRefFitSeries = fdRefFit + fdRateFit * timeOffsetSecUse;
  fdRefResidual = refFd - fdRefFitSeries;
  [deltaFdRef, deltaFdRate] = localFitSatFdLine(timeOffsetSecUse, deltaFd);
end

refStateIdx = localResolveReferenceStateIdx(timeOffsetSecUse, numState);

dopplerState = struct();
dopplerState.refState = refState;
dopplerState.refSatIdxLocal = refSatIdxLocal;
dopplerState.refSatIdxGlobal = localGetStructField(refState, 'satIdxGlobal', NaN);
dopplerState.fdRef = refFd;
dopplerState.fdSat = satFd;
dopplerState.deltaFd = deltaFd;
dopplerState.fdRefFit = fdRefFit;
dopplerState.fdRateFit = fdRateFit;
dopplerState.fdRefFitSeries = fdRefFitSeries;
dopplerState.fdRefResidual = fdRefResidual;
dopplerState.deltaFdRef = deltaFdRef;
dopplerState.deltaFdRate = deltaFdRate;
dopplerState.timeOffsetSec = timeOffsetSecUse;
dopplerState.refStateIdx = refStateIdx;
dopplerState.fdRefRefFrame = refFd(refStateIdx);
dopplerState.fdSatRefFrame = satFd(:, refStateIdx);
dopplerState.deltaFdRefFrame = deltaFd(:, refStateIdx);
end


function [refSatPosResolve, refSatVelResolve] = localGetReferenceResolveState(scene, satPosEci, satVelEci)
%LOCALGETREFERENCERESOLVESTATE Get the 2-D satellite state used by resolver.

if isfield(scene, 'satPosEci') && isfield(scene, 'satVelEci') && ...
    isnumeric(scene.satPosEci) && isnumeric(scene.satVelEci) && ...
    isequal(size(scene.satPosEci), size(scene.satVelEci)) && ...
    size(scene.satPosEci, 1) == 3 && ndims(scene.satPosEci) == 2
  refSatPosResolve = scene.satPosEci;
  refSatVelResolve = scene.satVelEci;
  return;
end

if ndims(satPosEci) == 2
  refSatPosResolve = satPosEci;
  refSatVelResolve = satVelEci;
  return;
end

refSatPosResolve = satPosEci(:, :, 1);
refSatVelResolve = satVelEci(:, :, 1);
end


function localCheckReferenceInvariant(deltaFd, satFd, refFd, refSatIdxLocal)
%LOCALCHECKREFERENCEINVARIANT Check the shared reference-Doppler semantics.

tol = 1e-8 * max(1, max(abs(satFd(:))));
refFdMat = repmat(refFd, size(satFd, 1), 1);

if any(abs(deltaFd(refSatIdxLocal, :)) > tol)
  error('buildReferenceDopplerState:InvalidReferenceDeltaFd', ...
    'Reference satellite must satisfy deltaFd(refSat,:) = 0.');
end

if any(abs(satFd(refSatIdxLocal, :) - refFd) > tol)
  error('buildReferenceDopplerState:InvalidReferenceFd', ...
    'Reference satellite must satisfy fdSat(refSat,:) = fdRef.');
end

if any(abs(satFd(:) - (refFdMat(:) + deltaFd(:))) > tol)
  error('buildReferenceDopplerState:InvalidFdComposition', ...
    'All satellites must satisfy fdSat = fdRef + deltaFd.');
end
end


function [fd0, fdRate] = localFitSatFdLine(timeOffsetSec, fdSeriesMat)
%LOCALFITSATFDLINE Fit one Doppler line per satellite.

numSat = size(fdSeriesMat, 1);
fd0 = zeros(numSat, 1);
fdRate = zeros(numSat, 1);
for iSat = 1:numSat
  [fd0(iSat), fdRate(iSat)] = localFitFdLine(timeOffsetSec, fdSeriesMat(iSat, :));
end
end


function [fd0, fdRate] = localFitFdLine(timeOffsetSec, fdSeries)
%LOCALFITFDLINE Fit fd(t) = fd0 + fdRate * t by least squares.

timeOffsetSec = reshape(timeOffsetSec, [], 1);
fdSeries = reshape(fdSeries, [], 1);
regMat = [ones(numel(timeOffsetSec), 1), timeOffsetSec];
coef = regMat \ fdSeries;
fd0 = coef(1);
fdRate = coef(2);
end



function refStateIdx = localResolveReferenceStateIdx(timeOffsetSec, numState)
%LOCALRESOLVEREFERENCESTATEIDX Resolve which Doppler state is the reference frame.

if isempty(timeOffsetSec)
  refStateIdx = 1;
  return;
end

timeOffsetSec = reshape(timeOffsetSec, [], 1);
[~, refStateIdx] = min(abs(timeOffsetSec));
if isempty(refStateIdx) || ~isfinite(refStateIdx)
  refStateIdx = 1;
end
refStateIdx = min(max(round(refStateIdx), 1), numState);
end



function refSatIdxLocal = localValidateReferenceSatIdx(refSatIdxLocalOverride, satPosEci)
%LOCALVALIDATEREFERENCESATIDX Validate one explicit local ref-satellite index.

numSat = size(satPosEci, 2);
refSatIdxLocal = round(refSatIdxLocalOverride(1));
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && ...
    refSatIdxLocal >= 1 && refSatIdxLocal <= numSat)
  error('buildReferenceDopplerState:InvalidReferenceIndexOverride', ...
    'refSatIdxLocalOverride must be a valid local satellite index.');
end
end


function refState = localBuildReferenceStateFromIdx(scene, satPosEci, satVelEci, refSatIdxLocal, sourceTag)
%LOCALBUILDREFERENCESTATEFROMIDX Build one-hot reference metadata from an index.

refWeight = zeros(size(satPosEci, 2), 1);
refWeight(refSatIdxLocal) = 1;

refState = struct();
refState.type = 'selectedSat';
refState.weight = refWeight;
refState.posEci = reshape(localSelectSatState(satPosEci, refSatIdxLocal), 3, 1);
refState.velEci = reshape(localSelectSatState(satVelEci, refSatIdxLocal), 3, 1);
refState.satIdxLocal = refSatIdxLocal;
refState.satIdxGlobal = localResolveGlobalSatIdx(scene, refSatIdxLocal);
refState.source = char(sourceTag);
end


function satState = localSelectSatState(stateArray, satIdx)
%LOCALSELECTSATSTATE Select one satellite state from 2-D or 3-D arrays.

if ndims(stateArray) == 2
  satState = stateArray(:, satIdx);
else
  satState = stateArray(:, satIdx, 1);
end
end


function satIdxGlobal = localResolveGlobalSatIdx(scene, refSatIdxLocal)
%LOCALRESOLVEGLOBALSATIDX Resolve the global satellite index when available.

satIdxGlobal = NaN;
if isstruct(scene) && isfield(scene, 'satIdx') && ~isempty(scene.satIdx)
  satIdxVec = reshape(scene.satIdx, [], 1);
  if numel(satIdxVec) >= refSatIdxLocal
    satIdxGlobal = satIdxVec(refSatIdxLocal);
  end
end
end

function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct/object field with a default value.

fieldValue = defaultValue;

if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
