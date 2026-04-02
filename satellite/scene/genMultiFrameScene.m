function sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
  minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg, refFrameIdx, motionOpt)
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
%  sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, usrIdx, arrayTemplate, ...
%    minUsrElevationDeg, maxSatOffAxisDeg, refType, refArg, refFrameIdx, motionOpt)
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
%  motionOpt          - (optional) user-motion perturbation options in ECI
%                       coordinates. This is intended for dynamic stress
%                       tests while keeping the estimator parameterization
%                       unchanged.
%                       .enableUsrMotion  - logical flag, default false
%                       .usrVelEciExtra   - 3x1 or 3xNu extra user velocity
%                                           in m/s, default 0
%                       .usrAccEciExtra   - 3x1 or 3xNu extra user acceleration
%                                           in m/s^2, default 0
%                       .motionRefFrameIdx
%                                         - reference frame index defining
%                                           zero displacement time, default
%                                           refFrameIdx
%                       .recomputeAccess  - whether to recompute access and
%                                           localDoa after motion update,
%                                           default true
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
%                       .usrPosNominalEci
%                       .usrVelNominalEci
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
%                       .userMotion
%                       .sceneCell
%
%  sceneSeq.sceneCell - 1xNf cell array, each entry is the output of
%                       genMultiSatScene at the corresponding utcVec frame.
%                       When motionOpt.enableUsrMotion=true, the user state,
%                       access mask, and localDoa fields are overwritten by
%                       the motion-perturbed values.
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
  motionOpt (1,1) struct = struct()
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
% Optional user-motion perturbation in ECI
% -------------------------------------------------------------------------
motionInfo = localParseMotionOpt(motionOpt, numUser, numFrame, refFrameIdx);
usrPosNominalEci = [];
usrVelNominalEci = [];

if motionInfo.enableUsrMotion
  usrPosNominalEci = zeros(3, numUser, numFrame);
  usrVelNominalEci = zeros(3, numUser, numFrame);
  motionDtSec = seconds(utcVec - utcVec(motionInfo.motionRefFrameIdx));

  for iFrame = 1:numFrame
    [sceneCell{iFrame}, usrPosNominalEci(:, :, iFrame), usrVelNominalEci(:, :, iFrame)] = ...
      localApplyUserMotion(sceneCell{iFrame}, motionDtSec(iFrame), motionInfo, ...
      minUsrElevationDeg, maxSatOffAxisDeg);
  end

  sceneRef = sceneCell{refFrameIdx};
else
  motionInfo.motionRefFrameIdx = refFrameIdx;
end

% -------------------------------------------------------------------------
% Aggregate frame-dependent states
% -------------------------------------------------------------------------
satPosEci = zeros(3, numSat, numFrame);
satVelEci = zeros(3, numSat, numFrame);
usrPosEci = zeros(3, numUser, numFrame);
usrVelEci = zeros(3, numUser, numFrame);

if isempty(usrPosNominalEci)
  usrPosNominalEci = zeros(3, numUser, numFrame);
  usrVelNominalEci = zeros(3, numUser, numFrame);
end

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

  if isfield(sceneTmp, 'usrPosNominalEci') && ~isempty(sceneTmp.usrPosNominalEci)
    usrPosNominalEci(:, :, iFrame) = sceneTmp.usrPosNominalEci;
  else
    usrPosNominalEci(:, :, iFrame) = sceneTmp.usrPosEci;
  end

  if isfield(sceneTmp, 'usrVelNominalEci') && ~isempty(sceneTmp.usrVelNominalEci)
    usrVelNominalEci(:, :, iFrame) = sceneTmp.usrVelNominalEci;
  else
    usrVelNominalEci(:, :, iFrame) = sceneTmp.usrVelEci;
  end

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
sceneSeq.usrPosNominalEci = usrPosNominalEci;
sceneSeq.usrVelNominalEci = usrVelNominalEci;

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

sceneSeq.userMotion = motionInfo;
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

function motionInfo = localParseMotionOpt(motionOpt, numUser, numFrame, refFrameIdx)
%LOCALPARSEMOTIONOPT Parse optional user-motion perturbation options.

motionInfo = struct();
motionInfo.enableUsrMotion = false;
motionInfo.usrVelEciExtra = zeros(3, numUser);
motionInfo.usrAccEciExtra = zeros(3, numUser);
motionInfo.motionRefFrameIdx = refFrameIdx;
motionInfo.recomputeAccess = true;

if isempty(fieldnames(motionOpt))
  return;
end

if isfield(motionOpt, 'usrVelEciExtra') && ~isempty(motionOpt.usrVelEciExtra)
  motionInfo.usrVelEciExtra = localExpandMotionState(motionOpt.usrVelEciExtra, numUser, 'usrVelEciExtra');
end

if isfield(motionOpt, 'usrAccEciExtra') && ~isempty(motionOpt.usrAccEciExtra)
  motionInfo.usrAccEciExtra = localExpandMotionState(motionOpt.usrAccEciExtra, numUser, 'usrAccEciExtra');
end

if isfield(motionOpt, 'recomputeAccess') && ~isempty(motionOpt.recomputeAccess)
  motionInfo.recomputeAccess = logical(motionOpt.recomputeAccess);
end

if isfield(motionOpt, 'motionRefFrameIdx') && ~isempty(motionOpt.motionRefFrameIdx)
  motionRefFrameIdx = motionOpt.motionRefFrameIdx;
  if ~isscalar(motionRefFrameIdx) || ~isnumeric(motionRefFrameIdx) || ...
      ~isfinite(motionRefFrameIdx) || mod(motionRefFrameIdx, 1) ~= 0 || ...
      motionRefFrameIdx < 1 || motionRefFrameIdx > numFrame
    error('genMultiFrameScene:InvalidMotionRefFrameIdx', ...
      'motionOpt.motionRefFrameIdx must be an integer within [1, numFrame].');
  end
  motionInfo.motionRefFrameIdx = motionRefFrameIdx;
end

if isfield(motionOpt, 'enableUsrMotion') && ~isempty(motionOpt.enableUsrMotion)
  motionInfo.enableUsrMotion = logical(motionOpt.enableUsrMotion);
else
  motionInfo.enableUsrMotion = any(abs(motionInfo.usrVelEciExtra(:)) > 0) || ...
    any(abs(motionInfo.usrAccEciExtra(:)) > 0);
end
end

function stateMat = localExpandMotionState(stateIn, numUser, fieldName)
%LOCALEXPANDMOTIONSTATE Expand 3x1 motion state to 3xNu if needed.

validateattributes(stateIn, {'numeric'}, {'real', 'finite', 'nonnan'});

if isequal(size(stateIn), [3, 1])
  stateMat = repmat(stateIn, 1, numUser);
elseif isequal(size(stateIn), [3, numUser])
  stateMat = stateIn;
else
  error('genMultiFrameScene:InvalidMotionStateSize', ...
    '%s must have size 3x1 or 3xNu.', fieldName);
end
end

function [sceneOut, usrPosNominalEci, usrVelNominalEci] = localApplyUserMotion( ...
  sceneIn, dtSec, motionInfo, minUsrElevationDeg, maxSatOffAxisDeg)
%LOCALAPPLYUSERMOTION Apply extra ECI motion to one frame scene.

sceneOut = sceneIn;
usrPosNominalEci = sceneIn.usrPosEci;
usrVelNominalEci = sceneIn.usrVelEci;

usrPosEci = usrPosNominalEci + motionInfo.usrVelEciExtra * dtSec + ...
  0.5 * motionInfo.usrAccEciExtra * (dtSec ^ 2);
usrVelEci = usrVelNominalEci + motionInfo.usrVelEciExtra + ...
  motionInfo.usrAccEciExtra * dtSec;

sceneOut.usrPosNominalEci = usrPosNominalEci;
sceneOut.usrVelNominalEci = usrVelNominalEci;
sceneOut.usrPosEci = usrPosEci;
sceneOut.usrVelEci = usrVelEci;

if ~motionInfo.recomputeAccess
  return;
end

accessInfo = checkSatAccess(usrPosEci, sceneIn.satPosEci, ...
  minUsrElevationDeg, maxSatOffAxisDeg);
sceneOut.accessInfo = accessInfo;
sceneOut.access = accessInfo.isAvailable;

numSat = sceneIn.numSat;
numUser = sceneIn.numUser;

if numUser == 1
  localDoa = zeros(2, numSat);
else
  localDoa = zeros(2, numSat, numUser);
end

for iSat = 1:numSat
  doaMat = globalToLocalDoa(usrPosEci, sceneIn.satPosEci(:, iSat), sceneIn.rotMat{iSat});
  if numUser == 1
    localDoa(:, iSat) = doaMat;
  else
    localDoa(:, iSat, :) = reshape(doaMat, 2, 1, numUser);
  end
end

sceneOut.localDoa = localDoa;
end
