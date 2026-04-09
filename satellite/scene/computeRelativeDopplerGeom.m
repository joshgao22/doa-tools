function [deltaFd, satFd, refFd] = computeRelativeDopplerGeom( ...
  satPosEci, satVelEci, userPosEci, userVelEci, wavelength, refSatIdxLocal)
%COMPUTERELATIVEDOPPLERGEOM Compute exact reference-relative Doppler.
% The main parameterization in this project uses one selected reference
% satellite. This helper computes the exact geometric Doppler of every
% satellite link and then converts it into the reference-relative form
%
%   deltaFd = fdSat - fdRef
%
% so that the reference satellite always satisfies deltaFd = 0.
%
%Syntax:
%  [deltaFd, satFd, refFd] = computeRelativeDopplerGeom( ...
%    satPosEci, satVelEci, userPosEci, userVelEci, wavelength, refSatIdxLocal)
%
%Inputs:
%  satPosEci      - 3xNs or 3xNsxNf satellite ECI positions
%  satVelEci      - 3xNs or 3xNsxNf satellite ECI velocities
%  userPosEci     - 3xK candidate user positions. K is the number of static
%                   candidates or the number of frames for dynamic mode.
%  userVelEci     - 3xK candidate user velocities aligned with userPosEci
%  wavelength     - carrier wavelength in meters
%  refSatIdxLocal - local reference-satellite index
%
%Outputs:
%  deltaFd        - Ns x K reference-relative Doppler in Hz
%  satFd          - Ns x K absolute per-satellite Doppler in Hz
%  refFd          - 1 x K reference-satellite Doppler in Hz
%
%See also:
%  resolveReferenceSatState, estimatorDoaDopplerMlePilotSfOpt,
%  estimatorDoaDopplerMlePilotMfOpt

arguments
  satPosEci {mustBeNumeric, mustBeFinite}
  satVelEci {mustBeNumeric, mustBeFinite}
  userPosEci (3,:) {mustBeNumeric, mustBeFinite}
  userVelEci (3,:) {mustBeNumeric, mustBeFinite}
  wavelength (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  refSatIdxLocal (1,1) {mustBePositive, mustBeInteger}
end

if ~isequal(size(userPosEci), size(userVelEci))
  error('computeRelativeDopplerGeom:UserStateSizeMismatch', ...
    'userPosEci and userVelEci must have the same size.');
end

satDim = ndims(satPosEci);
if ~isequal(size(satPosEci), size(satVelEci))
  error('computeRelativeDopplerGeom:SatStateSizeMismatch', ...
    'satPosEci and satVelEci must have the same size.');
end
if size(satPosEci, 1) ~= 3
  error('computeRelativeDopplerGeom:InvalidSatStateSize', ...
    'satPosEci and satVelEci must start with size 3 x Ns.');
end

numSat = size(satPosEci, 2);
if refSatIdxLocal < 1 || refSatIdxLocal > numSat
  error('computeRelativeDopplerGeom:InvalidReferenceIndex', ...
    'refSatIdxLocal must be within [1, numSat].');
end

numState = size(userPosEci, 2);

deltaFd = zeros(numSat, numState);
satFd = zeros(numSat, numState);
refFd = zeros(1, numState);

switch satDim
  case 2
    for iState = 1:numState
      [deltaFd(:, iState), satFd(:, iState), refFd(iState)] = localComputeOne( ...
        satPosEci, satVelEci, userPosEci(:, iState), userVelEci(:, iState), ...
        wavelength, refSatIdxLocal);
    end

  case 3
    numFrame = size(satPosEci, 3);
    if ~(numState == 1 || numState == numFrame)
      error('computeRelativeDopplerGeom:FrameCountMismatch', ...
        ['When satPosEci is 3-D, the number of user states must be either ', ...
         '1 or equal to the number of frames.']);
    end

    for iFrame = 1:numFrame
      stateIdx = localResolveStateIdx(numState, iFrame);
      [deltaFd(:, iFrame), satFd(:, iFrame), refFd(iFrame)] = localComputeOne( ...
        satPosEci(:, :, iFrame), satVelEci(:, :, iFrame), ...
        userPosEci(:, stateIdx), userVelEci(:, stateIdx), wavelength, refSatIdxLocal);
    end

  otherwise
    error('computeRelativeDopplerGeom:InvalidSatStateDim', ...
      'satPosEci and satVelEci must be 2-D or 3-D arrays.');
end
end

function stateIdx = localResolveStateIdx(numState, frameIdx)
%LOCALRESOLVESTATEIDX Resolve which user state should pair with one frame.

if numState == 1
  stateIdx = 1;
else
  stateIdx = frameIdx;
end
end


function [deltaFd, satFd, refFd] = localComputeOne( ...
  satPosEci, satVelEci, userPosEci, userVelEci, wavelength, refSatIdxLocal)
%LOCALCOMPUTEONE Compute one exact reference-relative Doppler slice.

numSat = size(satPosEci, 2);
satFd = zeros(numSat, 1);

for iSat = 1:numSat
  satFd(iSat) = localComputeFd( ...
    satPosEci(:, iSat), satVelEci(:, iSat), userPosEci, userVelEci, wavelength);
end

refFd = satFd(refSatIdxLocal);
deltaFd = satFd - refFd;
deltaFd(refSatIdxLocal) = 0;
end


function fdHz = localComputeFd(rxPosEci, rxVelEci, userPosEci, userVelEci, wavelength)
%LOCALCOMPUTEFD Compute one exact geometric Doppler value.

losVec = rxPosEci - userPosEci;
rangeVal = norm(losVec);
if rangeVal <= 0
  error('computeRelativeDopplerGeom:InvalidRange', ...
    'Receiver and user positions must not coincide.');
end
losUnit = losVec / rangeVal;
relVel = rxVelEci - userVelEci;
rangeRate = losUnit' * relVel;
fdHz = -rangeRate / wavelength;
end
