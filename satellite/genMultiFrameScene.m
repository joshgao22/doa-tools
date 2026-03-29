function sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
  minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg, refFrameIdx)
%GENMULTIFRAMESCENE Generate a multi-frame satellite / user geometry sequence.
% Repeatedly calls genMultiSatScene on a sequence of UTC epochs and packs
% all frame-dependent geometry into one structure for dynamic simulation and
% multi-frame joint estimation.
%
%Syntax:
%  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate)
%  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
%    minUsrElevationDeg, maxSatOffAxisDeg)
%  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
%    minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg)
%  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
%    minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg, refFrameIdx)
%
%Inputs:
%  utcVec             - 1xNf or Nfx1 datetime vector, frame epochs in UTC
%
%  tle                - TLE object / table accepted by propagateOrbit
%
%  usrLla             - 3xNu user geodetic coordinates [lat; lon; alt]
%                       latitude / longitude in degrees, altitude in meters
%
%  satIdx             - selected satellite indices
%                       - []      : use all propagated satellites
%                       - vector  : subset indices into propagated satellites
%
%  usrIdx             - selected user indices
%                       - []      : use all input users
%                       - vector  : subset indices into usrLla
%
%  arrayTemplate      - array model(s) attached to each satellite
%                       - scalar struct      : replicated to all satellites
%                       - struct array       : one array per satellite
%                       - cell array         : one array per satellite
%                       - []                 : leave sceneSeq.array empty
%
%  minUsrElevationDeg - (optional) minimum user elevation angle in degrees
%                       default: 15
%
%  maxSatOffAxisDeg   - (optional) maximum satellite off-axis angle in degrees
%                       default: 55
%
%  refType            - (optional) reference receiver type
%                       - "centroid"         : equal-weight satellite centroid
%                       - "weightedCentroid" : weighted satellite centroid
%                       - "satellite"        : use one selected satellite
%                       default: "centroid"
%
%  refArg             - (optional) reference parameter passed to
%                       genMultiSatScene
%
%  refFrameIdx        - (optional) reference frame index used to define
%                       sceneSeq.utcRef, sceneSeq.timeOffsetSec, and the
%                       frozen local DOA / rotation snapshot
%                       default: round((numFrame + 1) / 2)
%
%Output:
%  sceneSeq           - structure with fields:
%                       .utcVec
%                       .numFrame
%                       .refFrameIdx
%                       .utcRef
%                       .timeOffsetSec
%                       .frameDtSec
%                       .numUser
%                       .usrIdx
%                       .usrLla
%                       .usrPosEci
%                       .usrVelEci
%                       .numSat
%                       .satIdx
%                       .satPosEci
%                       .satVelEci
%                       .access
%                       .allAccess
%                       .anyAccess
%                       .frameAllAccess
%                       .rotMat
%                       .localDoa
%                       .localDoaRef
%                       .rotMatRef
%                       .array
%                       .ref
%                       .refPosEci
%                       .refVelEci
%                       .sceneCell
%
%  sceneSeq.sceneCell - 1xNf cell array, each entry is the output of
%                       genMultiSatScene at the corresponding utcVec frame
%
%  sceneSeq.access    - access indicators for all frames
%                       size: NsxNuxNf
%
%  sceneSeq.localDoa  - local DOAs for all frames
%                       - single-user : 2xNsxNf
%                       - multi-user  : 2xNsxNuxNf
%
%  sceneSeq.rotMat    - NfxNs cell array of satellite local-frame rotations
%
%Notes:
%  - This function does not discard frames with blocked links. Access
%    information is stored in sceneSeq.access and can be used later by the
%    signal generator or estimator.
%  - sceneSeq.localDoaRef and sceneSeq.rotMatRef are taken from the
%    reference frame refFrameIdx, which is convenient for frozen-steering
%    approximations in dynamic DoA-Doppler estimation.
%  - The selected user / satellite sets are fixed across all frames.
%
%See also:
%  genMultiSatScene, getLinkParam, getSceneSteering

arguments
  utcVec datetime
  tle
  usrLla (3,:) {mustBeNumeric, mustBeReal, mustBeFinite}
  satIdx {mustBePositive, mustBeFinite} = []
  usrIdx {mustBePositive, mustBeFinite} = []
  arrayTemplate = []
  minUsrElevationDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 15
  maxSatOffAxisDeg (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 55
  refType (1,1) string {mustBeMember(refType, ["centroid", "weightedCentroid", "satellite"])} = "centroid"
  refArg = []
  refFrameIdx = []
end

% -------------------------------------------------------------------------
% Parse frame epochs
% -------------------------------------------------------------------------
numUserAll = size(usrLla, 2);
usrIdx = localResolveIndexInput(usrIdx, numUserAll);
satIdxOut = satIdx;

if ~isvector(utcVec) || isempty(utcVec)
  error('genMultiFrameScene:InvalidUtcVec', ...
    'utcVec must be a nonempty datetime vector.');
end

utcVec = reshape(utcVec, 1, []);

if any(isnat(utcVec))
  error('genMultiFrameScene:NaTUtc', ...
    'utcVec must not contain NaT values.');
end

numFrame = numel(utcVec);

if numFrame > 1
  frameDtSec = seconds(diff(utcVec));
  if any(frameDtSec <= 0)
    error('genMultiFrameScene:NonIncreasingUtc', ...
      'utcVec must be strictly increasing.');
  end
else
  frameDtSec = [];
end

if isempty(refFrameIdx)
  refFrameIdx = round((numFrame + 1) / 2);
end

if ~isscalar(refFrameIdx) || ~isnumeric(refFrameIdx) || ~isfinite(refFrameIdx) || ...
    mod(refFrameIdx, 1) ~= 0 || refFrameIdx < 1 || refFrameIdx > numFrame
  error('genMultiFrameScene:InvalidRefFrameIdx', ...
    'refFrameIdx must be an integer within [1, numFrame].');
end

utcRef = utcVec(refFrameIdx);
timeOffsetSec = seconds(utcVec - utcRef);

% -------------------------------------------------------------------------
% Generate per-frame scenes
% -------------------------------------------------------------------------
sceneCell = cell(1, numFrame);

for iFrame = 1:numFrame
  sceneCell{iFrame} = genMultiSatScene( ...
    utcVec(iFrame), tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
    minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg);
end

% -------------------------------------------------------------------------
% Basic consistency checks
% -------------------------------------------------------------------------
sceneRef = sceneCell{refFrameIdx};
numSat = sceneRef.numSat;
numUser = sceneRef.numUser;

for iFrame = 1:numFrame
  localCheckFrameScene(sceneCell{iFrame}, numSat, numUser, iFrame);
end

% -------------------------------------------------------------------------
% Aggregate frame-dependent states
% -------------------------------------------------------------------------
satPosEci = zeros(3, numSat, numFrame);
satVelEci = zeros(3, numSat, numFrame);
usrPosEci = zeros(3, numUser, numFrame);
usrVelEci = zeros(3, numUser, numFrame);

access = false(numSat, numUser, numFrame);
rotMat = cell(numFrame, numSat);

refPosEci = zeros(3, numFrame);
refVelEci = zeros(3, numFrame);

if numUser == 1
  localDoa = zeros(2, numSat, numFrame);
else
  localDoa = zeros(2, numSat, numUser, numFrame);
end

for iFrame = 1:numFrame
  sceneTmp = sceneCell{iFrame};

  satPosEci(:, :, iFrame) = sceneTmp.satPosEci;
  satVelEci(:, :, iFrame) = sceneTmp.satVelEci;
  usrPosEci(:, :, iFrame) = sceneTmp.usrPosEci;
  usrVelEci(:, :, iFrame) = sceneTmp.usrVelEci;

  access(:, :, iFrame) = localExtractAccessMask(sceneTmp.access, numSat, numUser, iFrame);
  rotMat(iFrame, :) = reshape(sceneTmp.rotMat, 1, []);

  refPosEci(:, iFrame) = sceneTmp.ref.posEci;
  refVelEci(:, iFrame) = sceneTmp.ref.velEci;

  if numUser == 1
    localDoa(:, :, iFrame) = sceneTmp.localDoa;
  else
    localDoa(:, :, :, iFrame) = sceneTmp.localDoa;
  end
end

% -------------------------------------------------------------------------
% Aggregate access indicators
% -------------------------------------------------------------------------
allAccess = all(access, 3);
anyAccess = any(access, 3);
frameAllAccess = squeeze(all(reshape(access, [], numFrame), 1));

% -------------------------------------------------------------------------
% Pack outputs
% -------------------------------------------------------------------------
sceneSeq = struct();
sceneSeq.utcVec = utcVec;
sceneSeq.numFrame = numFrame;
sceneSeq.refFrameIdx = refFrameIdx;
sceneSeq.utcRef = utcRef;
sceneSeq.timeOffsetSec = timeOffsetSec;
sceneSeq.frameDtSec = frameDtSec;

sceneSeq.numUser = numUser;
sceneSeq.usrIdx = usrIdx;
sceneSeq.usrLla = sceneRef.usrLla;
sceneSeq.usrPosEci = usrPosEci;
sceneSeq.usrVelEci = usrVelEci;

sceneSeq.numSat = numSat;

if isempty(satIdxOut)
  sceneSeq.satIdx = 1:numSat;
else
  sceneSeq.satIdx = satIdxOut(:).';
end

sceneSeq.satPosEci = satPosEci;
sceneSeq.satVelEci = satVelEci;

sceneSeq.access = access;
sceneSeq.allAccess = allAccess;
sceneSeq.anyAccess = anyAccess;
sceneSeq.frameAllAccess = reshape(frameAllAccess, 1, []);

sceneSeq.rotMat = rotMat;
sceneSeq.localDoa = localDoa;
sceneSeq.localDoaRef = sceneRef.localDoa;
sceneSeq.rotMatRef = sceneRef.rotMat;

sceneSeq.array = sceneRef.array;

sceneSeq.ref = sceneRef.ref;
sceneSeq.refPosEci = refPosEci;
sceneSeq.refVelEci = refVelEci;

sceneSeq.sceneCell = sceneCell;
end

function localCheckFrameScene(sceneTmp, numSat, numUser, iFrame)
%LOCALCHECKFRAMESCENE Check size consistency of one frame scene.

if ~isstruct(sceneTmp)
  error('genMultiFrameScene:InvalidFrameScene', ...
    'sceneCell{%d} must be a struct returned by genMultiSatScene.', iFrame);
end

if ~isfield(sceneTmp, 'numSat') || ~isfield(sceneTmp, 'numUser')
  error('genMultiFrameScene:MissingFrameMeta', ...
    'sceneCell{%d} is missing numSat or numUser.', iFrame);
end

if sceneTmp.numSat ~= numSat || sceneTmp.numUser ~= numUser
  error('genMultiFrameScene:InconsistentFrameSize', ...
    'Frame %d has inconsistent numSat / numUser.', iFrame);
end

if ~isfield(sceneTmp, 'satPosEci') || ~isequal(size(sceneTmp.satPosEci), [3, numSat])
  error('genMultiFrameScene:InvalidSatPosSize', ...
    'sceneCell{%d}.satPosEci must have size 3xNs.', iFrame);
end

if ~isfield(sceneTmp, 'satVelEci') || ~isequal(size(sceneTmp.satVelEci), [3, numSat])
  error('genMultiFrameScene:InvalidSatVelSize', ...
    'sceneCell{%d}.satVelEci must have size 3xNs.', iFrame);
end

if ~isfield(sceneTmp, 'usrPosEci') || ~isequal(size(sceneTmp.usrPosEci), [3, numUser])
  error('genMultiFrameScene:InvalidUsrPosSize', ...
    'sceneCell{%d}.usrPosEci must have size 3xNu.', iFrame);
end

if ~isfield(sceneTmp, 'usrVelEci') || ~isequal(size(sceneTmp.usrVelEci), [3, numUser])
  error('genMultiFrameScene:InvalidUsrVelSize', ...
    'sceneCell{%d}.usrVelEci must have size 3xNu.', iFrame);
end

if ~isfield(sceneTmp, 'access') || isempty(sceneTmp.access)
  error('genMultiFrameScene:InvalidAccessSize', ...
    'sceneCell{%d}.access must not be empty.', iFrame);
end

localExtractAccessMask(sceneTmp.access, numSat, numUser, iFrame);

if ~isfield(sceneTmp, 'rotMat') || ~iscell(sceneTmp.rotMat) || numel(sceneTmp.rotMat) ~= numSat
  error('genMultiFrameScene:InvalidRotMatSize', ...
    'sceneCell{%d}.rotMat must be a cell array with one entry per satellite.', iFrame);
end

if ~isfield(sceneTmp, 'ref') || ~isstruct(sceneTmp.ref) || ...
    ~isfield(sceneTmp.ref, 'posEci') || ~isfield(sceneTmp.ref, 'velEci') || ...
    ~isequal(size(sceneTmp.ref.posEci), [3, 1]) || ~isequal(size(sceneTmp.ref.velEci), [3, 1])
  error('genMultiFrameScene:InvalidRefState', ...
    'sceneCell{%d}.ref must contain 3x1 posEci and velEci.', iFrame);
end
end

function accessMask = localExtractAccessMask(accessInput, numSat, numUser, iFrame)
%LOCALEXTRACTACCESSMASK Extract logical access mask from one frame scene.

if islogical(accessInput) || isnumeric(accessInput)
  accessMask = logical(accessInput);

elseif isstruct(accessInput)
  if isfield(accessInput, 'isAvailable') && ~isempty(accessInput.isAvailable)
    accessMask = logical(accessInput.isAvailable);
  elseif isfield(accessInput, 'isVisible') && ~isempty(accessInput.isVisible)
    accessMask = logical(accessInput.isVisible);
  else
    error('genMultiFrameScene:InvalidAccessStruct', ...
      ['sceneCell{%d}.access must contain field isAvailable or isVisible ', ...
       'when it is a struct.'], iFrame);
  end

else
  error('genMultiFrameScene:UnsupportedAccessType', ...
    'sceneCell{%d}.access must be a logical matrix, numeric matrix, or struct.', iFrame);
end

if ~isequal(size(accessMask), [numSat, numUser])
  error('genMultiFrameScene:InvalidAccessSize', ...
    'sceneCell{%d}.access must have size NsxNu.', iFrame);
end
end

function idx = localResolveIndexInput(idxIn, numAll)
%LOCALRESOLVEINDEXINPUT Resolve empty index input to full index range.

if isempty(idxIn)
  idx = 1:numAll;
else
  idx = idxIn(:).';
end
end
