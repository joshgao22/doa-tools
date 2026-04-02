function truth = buildDynTruthFromLinkParam(linkParamCell, timeOffsetSec, sampleRate, userIdx)
%BUILDDYNTRUTHFROMLINKPARAM Build dynamic Doppler truth summaries from link parameters.
% Collects the frame-wise reference Doppler, per-satellite Doppler, and
% satellite-to-reference Doppler offsets from a multi-frame link-parameter
% sequence. A first-order fit of the reference Doppler trajectory is also
% returned, which is convenient for nuisance-rate simulations and dynamic
% truth visualization.
%
%Syntax:
%  truth = buildDynTruthFromLinkParam(linkParamCell, timeOffsetSec)
%  truth = buildDynTruthFromLinkParam(linkParamCell, timeOffsetSec, sampleRate)
%  truth = buildDynTruthFromLinkParam(linkParamCell, timeOffsetSec, sampleRate, userIdx)
%
%Inputs:
%  linkParamCell    - frame-wise link-parameter sequence, supported forms:
%                     - 1xNf / Nfx1 cell, each cell is one linkParam struct
%                     - 1xNf / Nfx1 struct array
%                     - scalar struct for a single frame
%
%  timeOffsetSec    - frame time offsets in seconds, relative to the
%                     reference epoch used by the dynamic model
%                     size must match the number of frames
%
%  sampleRate       - (optional) sampling rate in Hz. When provided, the
%                     output also contains normalized Doppler quantities in
%                     rad/sample and rad/sample^2
%                     default: []
%
%  userIdx          - (optional) user index to extract from each linkParam
%                     default: 1
%
%Output:
%  truth            - structure with fields:
%                     .numFrame
%                     .numSat
%                     .numUser
%                     .userIdx
%                     .timeOffsetSec
%                     .timeSpanSec
%                     .timeRadiusSec
%                     .fdRefSeries      : 1xNf reference Doppler in Hz
%                     .fdSatSeries      : NsxNf satellite Doppler in Hz
%                     .deltaFdSeries    : NsxNf differential Doppler in Hz
%                     .fdRefFit         : fitted reference Doppler at t = 0
%                     .fdRateFit        : fitted Doppler rate in Hz/s
%                     .fdRefFitSeries   : 1xNf fitted reference Doppler line
%                     .fdRefResidual    : 1xNf fit residual in Hz
%                     .deltaFdWin       : scalar window Doppler span in Hz
%
%  Optional fields when sampleRate is provided:
%                     .sampleRate
%                     .fdRefSeriesNorm  : 1xNf normalized Doppler in rad/sample
%                     .fdSatSeriesNorm  : NsxNf normalized Doppler in rad/sample
%                     .deltaFdSeriesNorm: NsxNf normalized differential Doppler
%                     .fdRefFitNorm     : fitted normalized Doppler in rad/sample
%                     .fdRateFitNorm    : fitted normalized Doppler rate in
%                                         rad/sample^2
%                     .deltaFdWinNorm   : scalar normalized window span in
%                                         rad/sample
%                     .quadPhaseWin     : scalar quadratic phase magnitude
%                                         at the window edge in radians
%
%Notes:
%  - deltaFdWin is defined as max(fdRefSeries) - min(fdRefSeries), i.e. the
%    actual Doppler span across the whole observation window.
%  - quadPhaseWin is defined as
%
%        pi * abs(fdRateFit) * timeRadiusSec^2,
%
%    which corresponds to the quadratic phase magnitude at the farthest
%    frame from the reference epoch.
%
%See also:
%  getLinkParam

arguments
  linkParamCell
  timeOffsetSec {mustBeNumeric, mustBeReal, mustBeFinite}
  sampleRate = []
  userIdx (1,1) {mustBePositive, mustBeInteger} = 1
end

% -------------------------------------------------------------------------
% Parse and validate frame sequence
% -------------------------------------------------------------------------
linkParamSeq = localNormalizeLinkParamSeq(linkParamCell);
numFrame = numel(linkParamSeq);
timeOffsetSec = reshape(timeOffsetSec, 1, []);

if numel(timeOffsetSec) ~= numFrame
  error('buildDynTruthFromLinkParam:TimeOffsetSizeMismatch', ...
    'numel(timeOffsetSec) must equal the number of frames in linkParamCell.');
end

[numSat, numUser] = localCheckFrameConsistency(linkParamSeq);
if userIdx > numUser
  error('buildDynTruthFromLinkParam:InvalidUserIdx', ...
    'userIdx exceeds the number of users in linkParamCell.');
end

if ~isempty(sampleRate)
  validateattributes(sampleRate, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'sampleRate');
end

% -------------------------------------------------------------------------
% Collect frame-wise Doppler truth
% -------------------------------------------------------------------------
fdRefSeries = zeros(1, numFrame);
fdSatSeries = zeros(numSat, numFrame);
deltaFdSeries = zeros(numSat, numFrame);

for iFrame = 1:numFrame
  linkParam = linkParamSeq{iFrame};
  fdRefSeries(iFrame) = linkParam.ref.fdGeom(userIdx);
  fdSatSeries(:, iFrame) = linkParam.fdGeom(:, userIdx);
  deltaFdSeries(:, iFrame) = linkParam.ref.deltaFdGeom(:, userIdx);
end

% -------------------------------------------------------------------------
% Fit reference Doppler trajectory: fdRef(t) ~= fdRefFit + fdRateFit * t
% -------------------------------------------------------------------------
[fdRefFit, fdRateFit] = localFitFdLine(timeOffsetSec, fdRefSeries);
fdRefFitSeries = fdRefFit + fdRateFit * timeOffsetSec;
fdRefResidual = fdRefSeries - fdRefFitSeries;

% -------------------------------------------------------------------------
% Window-level summaries
% -------------------------------------------------------------------------
timeSpanSec = max(timeOffsetSec) - min(timeOffsetSec);
timeRadiusSec = max(abs(timeOffsetSec));
deltaFdWin = max(fdRefSeries) - min(fdRefSeries);

% -------------------------------------------------------------------------
% Pack outputs
% -------------------------------------------------------------------------
truth = struct();
truth.numFrame = numFrame;
truth.numSat = numSat;
truth.numUser = numUser;
truth.userIdx = userIdx;
truth.timeOffsetSec = timeOffsetSec;
truth.timeSpanSec = timeSpanSec;
truth.timeRadiusSec = timeRadiusSec;
truth.fdRefSeries = fdRefSeries;
truth.fdSatSeries = fdSatSeries;
truth.deltaFdSeries = deltaFdSeries;
truth.fdRefFit = fdRefFit;
truth.fdRateFit = fdRateFit;
truth.fdRefFitSeries = fdRefFitSeries;
truth.fdRefResidual = fdRefResidual;
truth.deltaFdWin = deltaFdWin;

% Backward-compatible aliases for existing temporary scripts.
truth.refFd = fdRefSeries;
truth.satFd = fdSatSeries;
truth.deltaFd = deltaFdSeries;

% -------------------------------------------------------------------------
% Optional normalized quantities
% -------------------------------------------------------------------------
if isempty(sampleRate)
  return;
end

omegaScale = 2 * pi / sampleRate;
gammaScale = 2 * pi / (sampleRate ^ 2);

truth.sampleRate = sampleRate;
truth.fdRefSeriesNorm = omegaScale * fdRefSeries;
truth.fdSatSeriesNorm = omegaScale * fdSatSeries;
truth.deltaFdSeriesNorm = omegaScale * deltaFdSeries;
truth.fdRefFitNorm = omegaScale * fdRefFit;
truth.fdRateFitNorm = gammaScale * fdRateFit;
truth.deltaFdWinNorm = omegaScale * deltaFdWin;
truth.quadPhaseWin = pi * abs(fdRateFit) * (timeRadiusSec ^ 2);
end


function linkParamSeq = localNormalizeLinkParamSeq(linkParamCell)
%LOCALNORMALIZELINKPARAMSEQ Convert supported frame inputs to a row cell array.

if iscell(linkParamCell)
  if isempty(linkParamCell)
    error('buildDynTruthFromLinkParam:EmptyLinkParamCell', ...
      'linkParamCell must be non-empty.');
  end

  linkParamSeq = reshape(linkParamCell, 1, []);
else
  if ~isstruct(linkParamCell)
    error('buildDynTruthFromLinkParam:InvalidLinkParamType', ...
      'linkParamCell must be a struct, a struct array, or a cell array of structs.');
  end

  if isempty(linkParamCell)
    error('buildDynTruthFromLinkParam:EmptyLinkParamStruct', ...
      'linkParamCell must be non-empty.');
  end

  linkParamSeq = num2cell(reshape(linkParamCell, 1, []));
end

for iFrame = 1:numel(linkParamSeq)
  if ~isstruct(linkParamSeq{iFrame})
    error('buildDynTruthFromLinkParam:InvalidFrameType', ...
      'Each frame entry in linkParamCell must be a struct.');
  end
end
end


function [numSat, numUser] = localCheckFrameConsistency(linkParamSeq)
%LOCALCHECKFRAMECONSISTENCY Validate required fields and frame dimensions.

linkParamRef = linkParamSeq{1};
localCheckRequiredFields(linkParamRef);

numSat = localGetNumSat(linkParamRef);
numUser = localGetNumUser(linkParamRef);

for iFrame = 2:numel(linkParamSeq)
  linkParam = linkParamSeq{iFrame};
  localCheckRequiredFields(linkParam);

  if localGetNumSat(linkParam) ~= numSat
    error('buildDynTruthFromLinkParam:NumSatMismatch', ...
      'All frames in linkParamCell must have the same number of satellites.');
  end

  if localGetNumUser(linkParam) ~= numUser
    error('buildDynTruthFromLinkParam:NumUserMismatch', ...
      'All frames in linkParamCell must have the same number of users.');
  end
  
  if ~isequal(size(linkParam.fdGeom), [numSat, numUser]) || ...
      ~isequal(size(linkParam.ref.deltaFdGeom), [numSat, numUser]) || ...
      ~isequal(size(linkParam.ref.fdGeom), [1, numUser])
    error('buildDynTruthFromLinkParam:FieldSizeMismatch', ...
      'The Doppler field sizes in linkParamCell are inconsistent across frames.');
  end
end
end


function localCheckRequiredFields(linkParam)
%LOCALCHECKREQUIREDFIELDS Ensure Doppler fields exist in one linkParam frame.

if ~isfield(linkParam, 'fdGeom') || ~isnumeric(linkParam.fdGeom)
  error('buildDynTruthFromLinkParam:MissingFdGeom', ...
    'Each linkParam frame must contain numeric field fdGeom.');
end

if ~isfield(linkParam, 'ref') || ~isstruct(linkParam.ref)
  error('buildDynTruthFromLinkParam:MissingRef', ...
    'Each linkParam frame must contain struct field ref.');
end

if ~isfield(linkParam.ref, 'fdGeom') || ~isnumeric(linkParam.ref.fdGeom)
  error('buildDynTruthFromLinkParam:MissingRefFdGeom', ...
    'Each linkParam.ref must contain numeric field fdGeom.');
end

if ~isfield(linkParam.ref, 'deltaFdGeom') || ~isnumeric(linkParam.ref.deltaFdGeom)
  error('buildDynTruthFromLinkParam:MissingDeltaFdGeom', ...
    'Each linkParam.ref must contain numeric field deltaFdGeom.');
end
end


function numSat = localGetNumSat(linkParam)
%LOCALGETNUMSAT Get the number of satellites from one frame.

numSat = size(linkParam.fdGeom, 1);
end


function numUser = localGetNumUser(linkParam)
%LOCALGETNUMUSER Get the number of users from one frame.

numUser = size(linkParam.fdGeom, 2);
end


function [fdRefFit, fdRateFit] = localFitFdLine(timeOffsetSec, fdRefSeries)
%LOCALFITFDLINE Fit fdRefSeries ~= fdRefFit + fdRateFit * timeOffsetSec.

if numel(timeOffsetSec) <= 1 || max(timeOffsetSec) - min(timeOffsetSec) <= 0
  fdRefFit = fdRefSeries(1);
  fdRateFit = 0;
  return;
end

designMat = [ones(numel(timeOffsetSec), 1), timeOffsetSec(:)];
coef = designMat \ fdRefSeries(:);

fdRefFit = coef(1);
fdRateFit = coef(2);
end
