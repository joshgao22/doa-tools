function [satIdx, aux] = pickVisibleSatByElevation(accessInput, numPick, userIdx, accessMode, frameIdx)
%PICKVISIBLESATBYELEVATION Pick highest-elevation satellites under access constraints.
% Returns satellite indices sorted by user-side elevation angle in
% descending order. The function is intended to unify the common
% "pick-the-best visible satellites" logic used by performance scripts.
%
%Syntax:
%  satIdx = pickVisibleSatByElevation(accessInput)
%
%  satIdx = pickVisibleSatByElevation(accessInput, numPick)
%
%  satIdx = pickVisibleSatByElevation(accessInput, numPick, userIdx)
%
%  satIdx = pickVisibleSatByElevation(accessInput, numPick, userIdx, accessMode)
%
%  satIdx = pickVisibleSatByElevation(accessInput, numPick, userIdx, accessMode, frameIdx)
%
%  [satIdx, aux] = pickVisibleSatByElevation(...)
%
%Inputs:
%  accessInput       - supported access container:
%                      1) satAccess / accessInfo struct returned by
%                         checkSatAccess or findVisibleSatFromTle
%                      2) single-frame scene struct containing
%                         .accessInfo or .access
%                      3) multi-frame sceneSeq struct containing .sceneCell
%
%  numPick           - (optional) number of satellites to pick
%                      - inf / [] : keep all eligible satellites
%                      - scalar K : pick the top-K eligible satellites
%                      default: inf
%
%  userIdx           - (optional) selected user index
%                      default: 1
%
%  accessMode        - (optional) access constraint used for candidate
%                      selection
%                      - "available" : require both visible and in-beam
%                      - "visible"   : require elevation constraint only
%                      - "inBeam"    : require beam constraint only
%                      default: "available"
%
%  frameIdx          - (optional) frame index when accessInput is a
%                      multi-frame scene sequence
%                      - [] : use sceneSeq.refFrameIdx if available,
%                             otherwise use frame 1
%                      default: []
%
%Outputs:
%  satIdx            - selected satellite indices in the local indexing of
%                      the input access container, sorted by elevation from
%                      high to low
%
%  aux               - diagnostic structure with fields:
%                      .userIdx
%                      .frameIdx
%                      .accessMode
%                      .numPick
%                      .numCandidate
%                      .candidateMask
%                      .candidateSatIdx
%                      .candidateElevationDeg
%                      .usrElevationDeg
%                      .globalSatIdx
%                      .sourceType
%
%Notes:
%  - For scene / sceneSeq input, satIdx always refers to the local
%    satellite indexing of that scene object. If the object also carries
%    scene.satIdx, the corresponding propagated-global indices are returned
%    in aux.globalSatIdx.
%  - For sceneSeq input, the access information is taken from one frame
%    sceneCell{frameIdx}. The default frame is refFrameIdx when available.
%  - If numPick is finite and there are fewer than numPick eligible
%    satellites, the function throws an error instead of silently returning
%    a partial result.
%
%See also:
%  checkSatAccess, findVisibleSatFromTle, selectSatScene, selectSatSceneSeq

arguments
  accessInput (1,1) struct
  numPick = inf
  userIdx (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 1
  accessMode (1,1) string {mustBeMember(accessMode, ["available", "visible", "inBeam"])} = "available"
  frameIdx = []
end

% -------------------------------------------------------------------------
% Basic input check
% -------------------------------------------------------------------------
if isempty(numPick)
  numPick = inf;
end

if ~(isscalar(numPick) && isnumeric(numPick) && isreal(numPick) && ...
    ((isfinite(numPick) && numPick >= 1 && mod(numPick, 1) == 0) || isinf(numPick)))
  error('pickVisibleSatByElevation:InvalidNumPick', ...
    'numPick must be inf, [], or a positive integer scalar.');
end

if mod(userIdx, 1) ~= 0 || userIdx < 1
  error('pickVisibleSatByElevation:InvalidUserIdx', ...
    'userIdx must be a positive integer scalar.');
end

if ~isempty(frameIdx)
  if ~(isscalar(frameIdx) && isnumeric(frameIdx) && isreal(frameIdx) && ...
      isfinite(frameIdx) && frameIdx >= 1 && mod(frameIdx, 1) == 0)
    error('pickVisibleSatByElevation:InvalidFrameIdx', ...
      'frameIdx must be [] or a positive integer scalar.');
  end
end

[accessInfo, satIdxGlobal, sourceType, frameIdxUsed] = localResolveAccessInfo(accessInput, frameIdx);
numSat = size(accessInfo.usrElevationDeg, 1);
numUser = size(accessInfo.usrElevationDeg, 2);

if userIdx > numUser
  error('pickVisibleSatByElevation:UserIdxOutOfRange', ...
    'userIdx exceeds the number of users in the access input.');
end

candidateMaskAll = localExtractCandidateMask(accessInfo, accessMode);
candidateMask = candidateMaskAll(:, userIdx);
usrElevationDeg = accessInfo.usrElevationDeg(:, userIdx);

if ~isequal(size(candidateMask), [numSat, 1])
  candidateMask = reshape(candidateMask, numSat, 1);
end

candidateSatIdx = find(candidateMask).';
candidateElevationDeg = usrElevationDeg(candidateSatIdx);

if isempty(candidateSatIdx)
  if isfinite(numPick)
    error('pickVisibleSatByElevation:NoEligibleSatellite', ...
      'No satellite satisfies the requested accessMode for user %d.', userIdx);
  end

  satIdx = zeros(1, 0);
  aux = localBuildAux(userIdx, frameIdxUsed, accessMode, numPick, candidateMask, ...
    candidateSatIdx, candidateElevationDeg, usrElevationDeg, satIdxGlobal, sourceType);
  aux.globalSatIdx = zeros(1, 0);
  return;
end

[~, sortIdx] = sort(candidateElevationDeg, 'descend');
candidateSatIdx = candidateSatIdx(sortIdx);
candidateElevationDeg = candidateElevationDeg(sortIdx);

if isfinite(numPick)
  if numel(candidateSatIdx) < numPick
    error('pickVisibleSatByElevation:NotEnoughEligibleSatellite', ...
      ['Only %d satellite(s) satisfy accessMode="%s" for user %d, ', ...
       'but numPick=%d was requested.'], ...
      numel(candidateSatIdx), accessMode, userIdx, numPick);
  end
  satIdx = candidateSatIdx(1:numPick);
else
  satIdx = candidateSatIdx;
end

aux = localBuildAux(userIdx, frameIdxUsed, accessMode, numPick, candidateMask, ...
  candidateSatIdx, candidateElevationDeg, usrElevationDeg, satIdxGlobal, sourceType);
aux.globalSatIdx = satIdxGlobal(satIdx);
end

function [accessInfo, satIdxGlobal, sourceType, frameIdxUsed] = localResolveAccessInfo(accessInput, frameIdx)
%LOCALRESOLVEACCESSINFO Resolve supported access containers to one accessInfo struct.

frameIdxUsed = [];
satIdxGlobal = [];

if isfield(accessInput, 'usrElevationDeg')
  accessInfo = localNormalizeAccessInfo(accessInput, 'accessInput');
  sourceType = "accessInfo";
  satIdxGlobal = localBuildIdentitySatIdx(accessInfo);
  return;
end

if isfield(accessInput, 'numFrame') || isfield(accessInput, 'sceneCell')
  [accessInfo, satIdxGlobal, frameIdxUsed] = localResolveFromSceneSeq(accessInput, frameIdx);
  sourceType = "sceneSeq";
  return;
end

[accessInfo, satIdxGlobal] = localResolveFromScene(accessInput, 'accessInput');
sourceType = "scene";
end

function [accessInfo, satIdxGlobal] = localResolveFromScene(scene, sceneName)
%LOCALRESOLVEFROMSCENE Resolve access information from one single-frame scene.

if isfield(scene, 'accessInfo') && ~isempty(scene.accessInfo)
  accessInfo = localNormalizeAccessInfo(scene.accessInfo, sprintf('%s.accessInfo', sceneName));
elseif isfield(scene, 'access') && isstruct(scene.access) && ~isempty(scene.access)
  accessInfo = localNormalizeAccessInfo(scene.access, sprintf('%s.access', sceneName));
else
  error('pickVisibleSatByElevation:MissingAccessInfo', ...
    ['%s must provide a satAccess/accessInfo struct with field ', ...
     'usrElevationDeg, or a sceneSeq.sceneCell entry containing it.'], sceneName);
end

satIdxGlobal = localExtractGlobalSatIdx(scene, size(accessInfo.usrElevationDeg, 1));
end

function [accessInfo, satIdxGlobal, frameIdxUsed] = localResolveFromSceneSeq(sceneSeq, frameIdx)
%LOCALRESOLVEFROMSCENESEQ Resolve frame-wise access information from sceneSeq.

if ~isfield(sceneSeq, 'sceneCell') || isempty(sceneSeq.sceneCell) || ~iscell(sceneSeq.sceneCell)
  error('pickVisibleSatByElevation:MissingSceneCell', ...
    'sceneSeq input must contain a non-empty sceneCell to pick by elevation.');
end

numFrame = numel(sceneSeq.sceneCell);
if isempty(frameIdx)
  if isfield(sceneSeq, 'refFrameIdx') && ~isempty(sceneSeq.refFrameIdx)
    frameIdxUsed = sceneSeq.refFrameIdx;
  else
    frameIdxUsed = 1;
  end
else
  frameIdxUsed = frameIdx;
end

if frameIdxUsed < 1 || frameIdxUsed > numFrame
  error('pickVisibleSatByElevation:FrameIdxOutOfRange', ...
    'frameIdx must lie within [1, sceneSeq.numFrame].');
end

sceneFrame = sceneSeq.sceneCell{frameIdxUsed};
if ~isstruct(sceneFrame)
  error('pickVisibleSatByElevation:InvalidFrameScene', ...
    'sceneSeq.sceneCell{%d} must be a scene struct.', frameIdxUsed);
end

[accessInfo, satIdxGlobal] = localResolveFromScene(sceneFrame, sprintf('sceneSeq.sceneCell{%d}', frameIdxUsed));
end

function accessInfo = localNormalizeAccessInfo(accessInfoIn, accessName)
%LOCALNORMALIZEACCESSINFO Validate access-info structure.

if ~isstruct(accessInfoIn)
  error('pickVisibleSatByElevation:InvalidAccessType', ...
    '%s must be a struct.', accessName);
end

if ~isfield(accessInfoIn, 'usrElevationDeg') || isempty(accessInfoIn.usrElevationDeg)
  error('pickVisibleSatByElevation:MissingUsrElevation', ...
    '%s must contain usrElevationDeg.', accessName);
end

usrElevationDeg = accessInfoIn.usrElevationDeg;
if ~isnumeric(usrElevationDeg) || ~isreal(usrElevationDeg) || ndims(usrElevationDeg) ~= 2
  error('pickVisibleSatByElevation:InvalidUsrElevation', ...
    '%s.usrElevationDeg must be a real NsxNu matrix.', accessName);
end

[numSat, numUser] = size(usrElevationDeg);
if numSat < 1 || numUser < 1
  error('pickVisibleSatByElevation:EmptyUsrElevation', ...
    '%s.usrElevationDeg must not be empty.', accessName);
end

accessInfo = accessInfoIn;
accessInfo.usrElevationDeg = usrElevationDeg;

% Build missing access masks from elevation when possible.
if ~isfield(accessInfo, 'isVisible') || isempty(accessInfo.isVisible)
  accessInfo.isVisible = true(numSat, numUser);
end
if ~isfield(accessInfo, 'isInBeam') || isempty(accessInfo.isInBeam)
  accessInfo.isInBeam = true(numSat, numUser);
end
if ~isfield(accessInfo, 'isAvailable') || isempty(accessInfo.isAvailable)
  accessInfo.isAvailable = logical(accessInfo.isVisible) & logical(accessInfo.isInBeam);
end

accessInfo.isVisible = localValidateAccessMask(accessInfo.isVisible, numSat, numUser, sprintf('%s.isVisible', accessName));
accessInfo.isInBeam = localValidateAccessMask(accessInfo.isInBeam, numSat, numUser, sprintf('%s.isInBeam', accessName));
accessInfo.isAvailable = localValidateAccessMask(accessInfo.isAvailable, numSat, numUser, sprintf('%s.isAvailable', accessName));
end

function accessMask = localValidateAccessMask(accessMaskIn, numSat, numUser, maskName)
%LOCALVALIDATEACCESSMASK Validate one logical access mask.

if ~(islogical(accessMaskIn) || isnumeric(accessMaskIn))
  error('pickVisibleSatByElevation:InvalidAccessMaskType', ...
    '%s must be a logical or numeric matrix.', maskName);
end

if ~isequal(size(accessMaskIn), [numSat, numUser])
  error('pickVisibleSatByElevation:InvalidAccessMaskSize', ...
    '%s must have size NsxNu.', maskName);
end

accessMask = logical(accessMaskIn);
end

function satIdxGlobal = localBuildIdentitySatIdx(accessInfo)
%LOCALBUILDIDENTITYSATIDX Build default or propagated-global satellite indices.

numSat = size(accessInfo.usrElevationDeg, 1);

if isfield(accessInfo, 'satIdx') && ~isempty(accessInfo.satIdx)
  satIdxGlobal = accessInfo.satIdx;
elseif isfield(accessInfo, 'tleIdx') && ~isempty(accessInfo.tleIdx)
  satIdxGlobal = accessInfo.tleIdx;
else
  satIdxGlobal = 1:numSat;
end

if ~isvector(satIdxGlobal) || numel(satIdxGlobal) ~= numSat
  error('pickVisibleSatByElevation:InvalidGlobalSatIdx', ...
    'The propagated-global satellite index vector must contain one entry per satellite.');
end

satIdxGlobal = reshape(satIdxGlobal, 1, []);
end

function satIdxGlobal = localExtractGlobalSatIdx(scene, numSat)
%LOCALEXTRACTGLOBALSATIDX Extract propagated-global satellite indices when available.

if isfield(scene, 'satIdx') && ~isempty(scene.satIdx)
  satIdxGlobal = scene.satIdx;
  if ~isvector(satIdxGlobal) || numel(satIdxGlobal) ~= numSat
    error('pickVisibleSatByElevation:InvalidSceneSatIdx', ...
      'scene.satIdx must contain one entry per satellite.');
  end
  satIdxGlobal = reshape(satIdxGlobal, 1, []);
else
  satIdxGlobal = 1:numSat;
end
end

function candidateMask = localExtractCandidateMask(accessInfo, accessMode)
%LOCALEXTRACTCANDIDATEMASK Extract candidate mask for one access mode.

switch accessMode
  case "available"
    candidateMask = logical(accessInfo.isAvailable);
  case "visible"
    candidateMask = logical(accessInfo.isVisible);
  case "inBeam"
    candidateMask = logical(accessInfo.isInBeam);
  otherwise
    error('pickVisibleSatByElevation:UnsupportedAccessMode', ...
      'Unsupported accessMode: %s', accessMode);
end
end

function aux = localBuildAux(userIdx, frameIdxUsed, accessMode, numPick, candidateMask, ...
  candidateSatIdx, candidateElevationDeg, usrElevationDeg, satIdxGlobal, sourceType)
%LOCALBUILDAUX Build diagnostic output structure.

aux = struct();
aux.userIdx = userIdx;
aux.frameIdx = frameIdxUsed;
aux.accessMode = accessMode;
aux.numPick = numPick;
aux.numCandidate = numel(candidateSatIdx);
aux.candidateMask = reshape(candidateMask, [], 1);
aux.candidateSatIdx = reshape(candidateSatIdx, 1, []);
aux.candidateElevationDeg = reshape(candidateElevationDeg, 1, []);
aux.usrElevationDeg = reshape(usrElevationDeg, [], 1);
aux.globalSatIdx = reshape(satIdxGlobal, 1, []);
aux.sourceType = sourceType;
end
