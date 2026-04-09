function [userPosEci, userVelEci, aux] = buildUserStateFromLatlon(latlon, stateRef)
%BUILDUSERSTATEFROMLATLON Build motion-aware user ECI states from [lat; lon].
% This helper converts candidate ground coordinates into ECI position and
% velocity states while preserving any extra user-motion perturbation stored
% in the current scene/sceneSeq. The nominal Earth-fixed motion is computed
% from the candidate lat/lon itself, and then the scene-derived perturbation
% is added on top so that static and dynamic estimators use the same user
% motion semantics as the generated data.
%
%Syntax:
%  [userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, stateRef)
%  [userPosEci, userVelEci, aux] = buildUserStateFromLatlon(latlon, stateRef)
%
%Inputs:
%  latlon    - 2xK candidate [lat; lon] in degrees.
%  stateRef  - one of the following:
%                1) scene struct with field .utc or .utcVec and optional
%                   user-motion fields .usrPosEci, .usrVelEci,
%                   .usrPosNominalEci, .usrVelNominalEci
%                2) sceneSeq-like struct with .utcVec and optional user
%                   motion series fields matching genMultiFrameScene output
%                3) datetime scalar/vector, in which case no extra motion
%                   perturbation is applied.
%
%Outputs:
%  userPosEci - 3xN candidate ECI positions.
%  userVelEci - 3xN candidate ECI velocities.
%  aux        - diagnostic struct containing the nominal states and the
%               motion offsets added to them.
%
%Pairing rule:
%  - scalar time + K lat/lon candidates  -> evaluate all K points at the
%    same epoch (single-frame grid search).
%  - time vector + one lat/lon candidate -> evaluate one candidate across
%    the whole time grid (multi-frame dynamic mode).
%  - equal-sized time and lat/lon sets   -> pair them elementwise.
%
%See also:
%  estimatorDoaDopplerMlePilotSfOpt, estimatorDoaDopplerMlePilotMfOpt,
%  genMultiFrameScene

arguments
  latlon (2,:) {mustBeNumeric, mustBeFinite}
  stateRef
end

[utcVec, stateTag] = localResolveUtcVec(stateRef);
[llaMat, utcMat] = localBuildLlaInput(latlon, utcVec);

userPosNominal = lla2eci(llaMat, utcMat).';
userVelNominal = localComputeEarthFixedVel(userPosNominal);
[posOffsetEci, velOffsetEci, offsetTag] = localResolveMotionOffset(stateRef, size(userPosNominal, 2));

userPosEci = userPosNominal + posOffsetEci;
userVelEci = userVelNominal + velOffsetEci;

if nargout < 3
  return;
end

aux = struct();
aux.utcVec = utcVec;
aux.utcSource = stateTag;
aux.offsetSource = offsetTag;
aux.userPosNominalEci = userPosNominal;
aux.userVelNominalEci = userVelNominal;
aux.userPosOffsetEci = posOffsetEci;
aux.userVelOffsetEci = velOffsetEci;
end


function [utcData, stateTag] = localResolveUtcVec(stateRef)
%LOCALRESOLVEUTCVEC Resolve the UTC samples from scene metadata.

stateTag = "datetime";

if isa(stateRef, 'datetime')
  utcData = stateRef;
  return;
end

if ~isstruct(stateRef)
  error('buildUserStateFromLatlon:InvalidStateRef', ...
    'stateRef must be a datetime or a scene/sceneSeq-like struct.');
end

if isfield(stateRef, 'utcVec') && ~isempty(stateRef.utcVec)
  utcData = stateRef.utcVec;
  stateTag = "utcVec";
  return;
end

if isfield(stateRef, 'utc') && ~isempty(stateRef.utc)
  utcData = stateRef.utc;
  stateTag = "utc";
  return;
end

error('buildUserStateFromLatlon:MissingUtc', ...
  'stateRef must provide utc or utcVec.');
end


function [llaMat, utcMat] = localBuildLlaInput(latlon, utcData)
%LOCALBUILDLLAINPUT Build one LLA/UTC table for lla2eci.

[utcSampleMat, numUtc] = localNormalizeUtcSamples(utcData);
numDoa = size(latlon, 2);

if numUtc == 1
  llaMat = [latlon.', zeros(numDoa, 1)];
  utcMat = repmat(utcSampleMat(1, :), numDoa, 1);
  return;
end

if numDoa == 1
  llaMat = [repmat(latlon(1).', numUtc, 1), ...
            repmat(latlon(2).', numUtc, 1), ...
            zeros(numUtc, 1)];
  utcMat = utcSampleMat;
  return;
end

if numDoa == numUtc
  llaMat = [latlon.', zeros(numDoa, 1)];
  utcMat = utcSampleMat;
  return;
end

error('buildUserStateFromLatlon:LatlonUtcSizeMismatch', ...
  ['The lat/lon candidates and UTC samples must satisfy one of: ', ...
   'scalar UTC, scalar lat/lon, or equal counts.']);
end


function [utcMat, numUtc] = localNormalizeUtcSamples(utcData)
%LOCALNORMALIZEUTCSAMPLES Normalize UTC input to an Nx6 date-vector table.
%
%Supported UTC formats:
%  - datetime scalar/vector
%  - scalar datenum
%  - 1x6 date vector
%  - Nx6 date-vector table

if isa(utcData, 'datetime')
  utcMat = datevec(utcData(:));
  numUtc = size(utcMat, 1);
  return;
end

if ~isnumeric(utcData) || isempty(utcData)
  error('buildUserStateFromLatlon:InvalidUtcType', ...
    'UTC input must be datetime, scalar datenum, or a 1x6/Nx6 date vector.');
end

if isvector(utcData) && numel(utcData) == 6
  utcMat = reshape(utcData, 1, 6);
  numUtc = 1;
  return;
end

if ismatrix(utcData) && size(utcData, 2) == 6
  utcMat = reshape(utcData, [], 6);
  numUtc = size(utcMat, 1);
  return;
end

if isscalar(utcData)
  utcMat = datevec(utcData);
  numUtc = 1;
  return;
end

error('buildUserStateFromLatlon:InvalidUtcShape', ...
  ['UTC input must be datetime, scalar datenum, or a 1x6/Nx6 date-vector ', ...
   'array.']);
end


function [posOffsetEci, velOffsetEci, offsetTag] = localResolveMotionOffset(stateRef, numState)
%LOCALRESOLVEMOTIONOFFSET Resolve scene-derived user-motion perturbations.

posOffsetEci = zeros(3, numState);
velOffsetEci = zeros(3, numState);
offsetTag = "none";

if ~isstruct(stateRef)
  return;
end

[posAct, hasPosAct] = localExtractUserSeries(stateRef, 'usrPosEci');
[velAct, hasVelAct] = localExtractUserSeries(stateRef, 'usrVelEci');
[posNom, hasPosNom] = localExtractUserSeries(stateRef, 'usrPosNominalEci');
[velNom, hasVelNom] = localExtractUserSeries(stateRef, 'usrVelNominalEci');

if hasPosAct && hasPosNom
  posOffsetEci = localExpandStateSeries(posAct - posNom, numState, 'usrPos offset');
  offsetTag = "positionNominal";
end

if hasVelAct && hasVelNom
  velOffsetEci = localExpandStateSeries(velAct - velNom, numState, 'usrVel offset');
  if offsetTag == "none"
    offsetTag = "velocityNominal";
  else
    offsetTag = "positionVelocityNominal";
  end
  return;
end

if hasVelAct && hasPosAct
  velOffsetBase = localComputeEarthFixedVel(posAct);
  velOffsetEci = localExpandStateSeries(velAct - velOffsetBase, numState, 'usrVel offset');
  if offsetTag == "none"
    offsetTag = "velocityInferred";
  else
    offsetTag = "positionVelocityInferred";
  end
end
end


function [stateMat, isFound] = localExtractUserSeries(stateRef, fieldName)
%LOCALEXTRACTUSERSERIES Extract a 3xN user-state series for user 1.

stateMat = zeros(3, 0);
isFound = false;

if ~isfield(stateRef, fieldName) || isempty(stateRef.(fieldName))
  return;
end

fieldVal = stateRef.(fieldName);
if ~isnumeric(fieldVal) || size(fieldVal, 1) ~= 3
  return;
end

switch ndims(fieldVal)
  case 2
    stateMat = reshape(fieldVal(:, 1), 3, 1);
    isFound = true;

  case 3
    stateMat = squeeze(fieldVal(:, 1, :));
    stateMat = reshape(stateMat, 3, []);
    isFound = true;

  otherwise
    % Keep not found.
end
end


function stateSeries = localExpandStateSeries(stateSeriesIn, numState, seriesName)
%LOCALEXPANDSTATESERIES Expand one 3xN state series to 3xnumState.

stateSeriesIn = reshape(stateSeriesIn, 3, []);
numAvail = size(stateSeriesIn, 2);

if numAvail == numState
  stateSeries = stateSeriesIn;
  return;
end

if numAvail == 1
  stateSeries = repmat(stateSeriesIn, 1, numState);
  return;
end

error('buildUserStateFromLatlon:StateSeriesSizeMismatch', ...
  '%s must contain either 1 or %d states.', seriesName, numState);
end


function velEci = localComputeEarthFixedVel(posEci)
%LOCALCOMPUTEEARTHFIXEDVEL Compute the ECI velocity of Earth-fixed points.

omegaEarth = 7.2921150e-5;
omegaVec = repmat([0; 0; omegaEarth], 1, size(posEci, 2));
velEci = cross(omegaVec, posEci, 1);
end
