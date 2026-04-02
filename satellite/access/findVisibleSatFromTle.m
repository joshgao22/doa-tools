function [satIdx, satAccess, satState] = findVisibleSatFromTle( ...
  utc, tleInput, usrLla, minUsrElevationDeg, maxSatOffAxisDeg, batchSize)
%FINDVISIBLESATFROMTLE Find visible / available satellites from a TLE catalog.
% Propagates all satellites in a TLE catalog to a specified UTC epoch,
% converts ground-user geodetic coordinates to ECI, evaluates access by
% checkSatAccess, and returns the satellite indices that are visible,
% in-beam, and available for each user.
%
%Syntax:
%  [satIdx, satAccess, satState] = findVisibleSatFromTle(utc, tleInput, usrLla)
%
%  ... = findVisibleSatFromTle(utc, tleInput, usrLla, ...
%    minUsrElevationDeg, maxSatOffAxisDeg)
%
%  ... = findVisibleSatFromTle(utc, tleInput, usrLla, ...
%    minUsrElevationDeg, maxSatOffAxisDeg, batchSize)
%
%Inputs:
%  utc                - 1x1 datetime, evaluation epoch in UTC
%
%  tleInput           - TLE source, supported forms:
%                       - char / string scalar : TLE file path
%                       - TLE object / table accepted by propagateOrbit
%
%  usrLla             - 3xNu user geodetic coordinates [lat; lon; alt]
%                       latitude / longitude in degrees, altitude in meters
%
%  minUsrElevationDeg - (optional) minimum user elevation angle in degrees
%                       default: 15
%
%  maxSatOffAxisDeg   - (optional) maximum satellite nadir off-axis angle
%                       in degrees
%                       default: 55
%
%  batchSize          - (optional) number of satellites processed in each
%                       access-check batch
%                       - [] / inf : process all satellites at once
%                       - scalar K : process satellites in batches of K
%                       default: []
%
%Outputs:
%  satIdx             - structure with per-user satellite indices:
%                       .visible   : satisfies elevation constraint only
%                       .inBeam    : satisfies off-axis constraint only
%                       .available : satisfies both constraints
%                       Single-user outputs are row vectors. Multi-user
%                       outputs are 1xNu cell arrays.
%
%  satAccess          - access-result structure with fields:
%                       .isAvailable
%                       .isVisible
%                       .isInBeam
%                       .usrElevationDeg
%                       .satOffAxisDeg
%                       .slantRange
%                       .numSat
%                       .numUser
%                       .tleIdx
%                       .satName
%
%  satState           - propagated state structure with fields:
%                       .utc
%                       .usrLla
%                       .usrPosEci
%                       .usrVelEci
%                       .satPosEci
%                       .satVelEci
%                       .numSat
%                       .numUser
%
%Notes:
%  - The returned satellite indices always refer to the propagated TLE
%    order, regardless of whether batch processing is used.
%  - Ground-user ECI velocity is approximated by omegaEarth x rEci, which
%    is consistent with genMultiSatScene.
%  - For large TLE catalogs, set batchSize to reduce the temporary memory
%    usage of the access evaluation.
%
%Example:
%  tle = tleread("./tle/starlink_all.tle");
%  utc = datetime(2026, 3, 18, 12, 0, 0, 'TimeZone', 'UTC');
%  usrLla = [39.9042; 116.4074; 0];
%  [satIdx, satAccess] = findVisibleSatFromTle(utc, tle, usrLla, 20, 55, 500);
%
%See also:
%  checkSatAccess, propagateOrbit, tleread, lla2eci

arguments
  utc (1,1) datetime
  tleInput
  usrLla (3,:) {mustBeNumeric, mustBeReal, mustBeFinite}
  minUsrElevationDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 15
  maxSatOffAxisDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 55
  batchSize = []
end

% -------------------------------------------------------------------------
% Basic input check
% -------------------------------------------------------------------------
if minUsrElevationDeg < 0 || minUsrElevationDeg > 90
  error('findVisibleSatFromTle:InvalidMinUsrElevation', ...
    'minUsrElevationDeg must lie within [0, 90].');
end

if maxSatOffAxisDeg < 0 || maxSatOffAxisDeg > 90
  error('findVisibleSatFromTle:InvalidMaxSatOffAxis', ...
    'maxSatOffAxisDeg must lie within [0, 90].');
end

if ~isempty(batchSize)
  if ~isscalar(batchSize) || ~isnumeric(batchSize) || ~isreal(batchSize) || ...
      ~isfinite(batchSize) || batchSize <= 0 || mod(batchSize, 1) ~= 0
    error('findVisibleSatFromTle:InvalidBatchSize', ...
      'batchSize must be [] or a positive integer scalar.');
  end
end

if isempty(utc.TimeZone)
  utc.TimeZone = 'UTC';
end

% -------------------------------------------------------------------------
% Parse TLE input and propagate satellite states
% -------------------------------------------------------------------------
parsedTle = localParseTleInput(tleInput);
[satPosAll, satVelAll] = propagateOrbit(utc, parsedTle);

satPosEci = reshape(satPosAll, 3, []);
satVelEci = reshape(satVelAll, 3, []);
numSat = size(satPosEci, 2);
numUser = size(usrLla, 2);

% -------------------------------------------------------------------------
% User ECI states
% -------------------------------------------------------------------------
utcVec = datevec(utc);
usrPosEci = zeros(3, numUser);
usrVelEci = zeros(3, numUser);

for iUser = 1:numUser
  usrPosTmp = lla2eci(usrLla(:, iUser).', utcVec);
  usrPosEci(:, iUser) = usrPosTmp(:);
  usrVelEci(:, iUser) = localComputeEarthFixedVel(usrPosEci(:, iUser));
end

% -------------------------------------------------------------------------
% Access check for all user-satellite pairs
% -------------------------------------------------------------------------
satAccess = localCheckAccessInBatch(usrPosEci, satPosEci, ...
  minUsrElevationDeg, maxSatOffAxisDeg, batchSize);

satAccess.tleIdx = 1:numSat;
satAccess.satName = localExtractSatName(parsedTle, numSat);

% -------------------------------------------------------------------------
% Build per-user satellite index lists
% -------------------------------------------------------------------------
satIdx = struct();
satIdx.visible = localBuildIndexList(satAccess.isVisible);
satIdx.inBeam = localBuildIndexList(satAccess.isInBeam);
satIdx.available = localBuildIndexList(satAccess.isAvailable);

% -------------------------------------------------------------------------
% Pack propagated states
% -------------------------------------------------------------------------
satState = struct();
satState.utc = utc;
satState.usrLla = usrLla;
satState.usrPosEci = usrPosEci;
satState.usrVelEci = usrVelEci;
satState.satPosEci = satPosEci;
satState.satVelEci = satVelEci;
satState.numSat = numSat;
satState.numUser = numUser;
end

function parsedTle = localParseTleInput(tleInput)
%LOCALPARSETLEINPUT Parse TLE file path or pass through TLE object/table.

if ischar(tleInput) || (isstring(tleInput) && isscalar(tleInput))
  tlePath = char(tleInput);
  if exist(tlePath, 'file') ~= 2
    error('findVisibleSatFromTle:TleFileNotFound', ...
      'The TLE file does not exist: %s', tlePath);
  end
  parsedTle = tleread(tlePath);
  return;
end

parsedTle = tleInput;
end

function usrVelEci = localComputeEarthFixedVel(usrPosEci)
%LOCALCOMPUTEEARTHFIXEDVEL Approximate ECI velocity of an Earth-fixed user.

omegaEarth = 7.2921150e-5;             % rad/s
omegaVec = [0; 0; omegaEarth];
usrVelEci = cross(omegaVec, usrPosEci);
end

function satAccess = localCheckAccessInBatch(usrPosEci, satPosEci, ...
  minUsrElevationDeg, maxSatOffAxisDeg, batchSize)
%LOCALCHECKACCESSINBATCH Run checkSatAccess in one shot or in batches.

numSat = size(satPosEci, 2);
numUser = size(usrPosEci, 2);

if isempty(batchSize) || isinf(batchSize) || batchSize >= numSat
  satAccess = checkSatAccess(usrPosEci, satPosEci, ...
    minUsrElevationDeg, maxSatOffAxisDeg);
  return;
end

isAvailable = false(numSat, numUser);
isVisible = false(numSat, numUser);
isInBeam = false(numSat, numUser);
usrElevationDeg = zeros(numSat, numUser);
satOffAxisDeg = zeros(numSat, numUser);
slantRange = zeros(numSat, numUser);

for startIdx = 1:batchSize:numSat
  stopIdx = min(startIdx + batchSize - 1, numSat);
  satSel = startIdx:stopIdx;

  accessTmp = checkSatAccess(usrPosEci, satPosEci(:, satSel), ...
    minUsrElevationDeg, maxSatOffAxisDeg);

  isAvailable(satSel, :) = accessTmp.isAvailable;
  isVisible(satSel, :) = accessTmp.isVisible;
  isInBeam(satSel, :) = accessTmp.isInBeam;
  usrElevationDeg(satSel, :) = accessTmp.usrElevationDeg;
  satOffAxisDeg(satSel, :) = accessTmp.satOffAxisDeg;
  slantRange(satSel, :) = accessTmp.slantRange;
end

satAccess = struct();
satAccess.isAvailable = isAvailable;
satAccess.isVisible = isVisible;
satAccess.isInBeam = isInBeam;
satAccess.usrElevationDeg = usrElevationDeg;
satAccess.satOffAxisDeg = satOffAxisDeg;
satAccess.slantRange = slantRange;
satAccess.numSat = numSat;
satAccess.numUser = numUser;
end

function satIdx = localBuildIndexList(maskMat)
%LOCALBUILDINDEXLIST Convert an NsxNu logical mask to index lists.

numUser = size(maskMat, 2);

if numUser == 1
  satIdx = find(maskMat(:, 1)).';
  return;
end

satIdx = cell(1, numUser);
for iUser = 1:numUser
  satIdx{iUser} = find(maskMat(:, iUser)).';
end
end

function satName = localExtractSatName(tleData, numSat)
%LOCALEXTRACTSATNAME Try to extract satellite names from the parsed TLE data.

satName = strings(1, numSat);

if istable(tleData)
  varName = localFindNameVariable(tleData.Properties.VariableNames);
  if ~isempty(varName) && height(tleData) == numSat
    satName = string(tleData.(varName));
    satName = reshape(satName, 1, []);
  end
  return;
end

if isstruct(tleData) && numel(tleData) == numSat
  if isfield(tleData, 'Name')
    satName = string({tleData.Name});
    return;
  end
  if isfield(tleData, 'name')
    satName = string({tleData.name});
    return;
  end
end
end

function varName = localFindNameVariable(varNameList)
%LOCALFINDNAMEVARIABLE Find a likely satellite-name variable in a table.

candidateList = {'Name', 'name', 'SatelliteName', 'satelliteName'};
varName = '';

for iName = 1:numel(candidateList)
  hit = strcmp(varNameList, candidateList{iName});
  if any(hit)
    varName = varNameList{find(hit, 1, 'first')};
    return;
  end
end
end
