function [linkParamScaled, aux] = scaleLinkParamDoppler( ...
  linkParamCell, timeOffsetSec, dopplerScale, timeRefSec)
%SCALELINKPARAMDOPPLER Scale Doppler dynamics in a link-parameter sequence.
% Scales the frame-to-frame Doppler evolution while keeping the spatial
% geometry, delay terms, and reference-time Doppler anchors unchanged. This
% is mainly used to scan dynamic strength in multi-frame simulations
% without rebuilding the scene geometry.
%
%Syntax:
%  linkParamScaled = scaleLinkParamDoppler(linkParamCell, timeOffsetSec, dopplerScale)
%
%  linkParamScaled = scaleLinkParamDoppler( ...
%    linkParamCell, timeOffsetSec, dopplerScale, timeRefSec)
%
%  [linkParamScaled, aux] = scaleLinkParamDoppler(...)
%
%Inputs:
%  linkParamCell   - frame-wise link-parameter sequence, supported forms:
%                    - 1xNf / Nfx1 cell, each cell is one linkParam struct
%                    - 1xNf / Nfx1 struct array
%                    - scalar struct for a single frame
%
%  timeOffsetSec   - frame time offsets in seconds, relative to the dynamic
%                    reference epoch
%                    numel(timeOffsetSec) must equal the number of frames
%
%  dopplerScale    - nonnegative scalar Doppler-dynamics scale factor
%                    - 0 : freeze Doppler to the fitted reference-time value
%                    - 1 : keep the original Doppler sequence unchanged
%                    - >1: exaggerate the frame-to-frame Doppler evolution
%
%  timeRefSec      - (optional) anchor time in seconds at which the Doppler
%                    values are kept unchanged after scaling
%                    default: 0
%
%Outputs:
%  linkParamScaled - scaled link-parameter sequence with the same container
%                    type and shape as the input. The following fields are
%                    updated when present:
%                    .fdGeom
%                    .rangeRate
%                    .ref.fdGeom
%                    .ref.rangeRate
%                    .ref.deltaFdGeom
%                    .ref.deltaRangeRate
%
%  aux             - auxiliary information with fields:
%                    .numFrame
%                    .numSat
%                    .numUser
%                    .timeOffsetSec
%                    .timeRefSec
%                    .dopplerScale
%                    .fdAnchor
%                    .fdRefAnchor
%                    .fdRefSeries
%                    .fdRefSeriesScaled
%                    .fdRateFit
%                    .fdRateFitScaled
%
%Notes:
%  - The scaling rule is applied as
%
%        yScaled(t) = yAnchor + dopplerScale * (y(t) - yAnchor),
%
%    where yAnchor is the affine-fit value at timeRefSec.
%
%  - Delay and direction related fields such as tauGeom, losEci, and
%    steering geometry are left unchanged.
%
%  - ref.deltaFdGeom and ref.deltaRangeRate are recomputed from the scaled
%    absolute quantities to keep consistency.
%
%See also:
%  getLinkParam, buildDynTruthFromLinkParam

arguments
  linkParamCell
  timeOffsetSec {mustBeNumeric, mustBeReal, mustBeFinite}
  dopplerScale (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonnegative}
  timeRefSec (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
end

% -------------------------------------------------------------------------
% Normalize frame container
% -------------------------------------------------------------------------
[linkParamSeq, restoreInfo] = localNormalizeLinkParamSeq(linkParamCell);
numFrame = numel(linkParamSeq);
timeOffsetSec = reshape(timeOffsetSec, 1, []);

if numel(timeOffsetSec) ~= numFrame
  error('scaleLinkParamDoppler:TimeOffsetSizeMismatch', ...
    'numel(timeOffsetSec) must equal the number of frames in linkParamCell.');
end

[numSat, numUser] = localCheckFrameConsistency(linkParamSeq);

% Early return for trivial sequence.
if numFrame == 0
  linkParamScaled = linkParamCell;
  aux = struct();
  return;
end

% -------------------------------------------------------------------------
% Collect Doppler-related sequences
% -------------------------------------------------------------------------
fdSeries = localCollectFieldSeries(linkParamSeq, {'fdGeom'}, [numSat, numUser]);
fdRefSeries = localCollectFieldSeries(linkParamSeq, {'ref', 'fdGeom'}, [1, numUser]);

hasRangeRate = localHasField(linkParamSeq{1}, {'rangeRate'});
hasRefRangeRate = localHasField(linkParamSeq{1}, {'ref', 'rangeRate'});

if hasRangeRate ~= hasRefRangeRate
  error('scaleLinkParamDoppler:InconsistentRangeRateFields', ...
    'rangeRate and ref.rangeRate must either both exist or both be absent.');
end

if hasRangeRate
  rangeRateSeries = localCollectFieldSeries(linkParamSeq, {'rangeRate'}, [numSat, numUser]);
  refRangeRateSeries = localCollectFieldSeries(linkParamSeq, {'ref', 'rangeRate'}, [1, numUser]);
else
  rangeRateSeries = [];
  refRangeRateSeries = [];
end

% -------------------------------------------------------------------------
% Build reference-time anchors and apply scaling
% -------------------------------------------------------------------------
fdAnchor = localAffineAnchor(fdSeries, timeOffsetSec, timeRefSec);
fdRefAnchor = localAffineAnchor(fdRefSeries, timeOffsetSec, timeRefSec);

fdSeriesScaled = fdAnchor + dopplerScale * (fdSeries - fdAnchor);
fdRefSeriesScaled = fdRefAnchor + dopplerScale * (fdRefSeries - fdRefAnchor);

if hasRangeRate
  rangeRateAnchor = localAffineAnchor(rangeRateSeries, timeOffsetSec, timeRefSec);
  refRangeRateAnchor = localAffineAnchor(refRangeRateSeries, timeOffsetSec, timeRefSec);

  rangeRateSeriesScaled = rangeRateAnchor + dopplerScale * (rangeRateSeries - rangeRateAnchor);
  refRangeRateSeriesScaled = refRangeRateAnchor + dopplerScale * (refRangeRateSeries - refRangeRateAnchor);
else
  rangeRateSeriesScaled = [];
  refRangeRateSeriesScaled = [];
end

% -------------------------------------------------------------------------
% Write back scaled quantities
% -------------------------------------------------------------------------
linkParamScaledSeq = linkParamSeq;
for iFrame = 1:numFrame
  fdMat = reshape(fdSeriesScaled(:, iFrame), numSat, numUser);
  fdRefVec = reshape(fdRefSeriesScaled(:, iFrame), 1, numUser);

  linkParamScaledSeq{iFrame}.fdGeom = fdMat;
  linkParamScaledSeq{iFrame}.ref.fdGeom = fdRefVec;
  linkParamScaledSeq{iFrame}.ref.deltaFdGeom = fdMat - fdRefVec;

  if hasRangeRate
    rangeRateMat = reshape(rangeRateSeriesScaled(:, iFrame), numSat, numUser);
    refRangeRateVec = reshape(refRangeRateSeriesScaled(:, iFrame), 1, numUser);

    linkParamScaledSeq{iFrame}.rangeRate = rangeRateMat;
    linkParamScaledSeq{iFrame}.ref.rangeRate = refRangeRateVec;
    linkParamScaledSeq{iFrame}.ref.deltaRangeRate = rangeRateMat - refRangeRateVec;
  end
end

linkParamScaled = localRestoreLinkParamSeq(linkParamScaledSeq, restoreInfo);

% -------------------------------------------------------------------------
% Pack auxiliary information
% -------------------------------------------------------------------------
[fdRefFit, fdRateFit] = localFitAffineSeries(fdRefSeries, timeOffsetSec, timeRefSec);
[fdRefFitScaled, fdRateFitScaled] = localFitAffineSeries(fdRefSeriesScaled, timeOffsetSec, timeRefSec);

aux = struct();
aux.numFrame = numFrame;
aux.numSat = numSat;
aux.numUser = numUser;
aux.timeOffsetSec = timeOffsetSec;
aux.timeRefSec = timeRefSec;
aux.dopplerScale = dopplerScale;
aux.fdAnchor = reshape(fdAnchor, numSat, numUser);
aux.fdRefAnchor = reshape(fdRefAnchor, 1, numUser);
aux.fdRefSeries = reshape(fdRefSeries, numUser, numFrame);
aux.fdRefSeriesScaled = reshape(fdRefSeriesScaled, numUser, numFrame);
aux.fdRefFit = reshape(fdRefFit, 1, numUser);
aux.fdRateFit = reshape(fdRateFit, 1, numUser);
aux.fdRefFitScaled = reshape(fdRefFitScaled, 1, numUser);
aux.fdRateFitScaled = reshape(fdRateFitScaled, 1, numUser);
end


function [linkParamSeq, restoreInfo] = localNormalizeLinkParamSeq(linkParamCell)
%LOCALNORMALIZELINKPARAMSEQ Convert supported inputs to a row cell sequence.

if iscell(linkParamCell)
  if isempty(linkParamCell)
    error('scaleLinkParamDoppler:EmptyLinkParamCell', ...
      'linkParamCell must be non-empty.');
  end

  restoreInfo.containerType = 'cell';
  restoreInfo.containerSize = size(linkParamCell);
  linkParamSeq = reshape(linkParamCell, 1, []);

elseif isstruct(linkParamCell)
  if isempty(linkParamCell)
    error('scaleLinkParamDoppler:EmptyLinkParamStruct', ...
      'linkParamCell must be non-empty.');
  end

  restoreInfo.containerType = 'struct';
  restoreInfo.containerSize = size(linkParamCell);
  linkParamSeq = num2cell(reshape(linkParamCell, 1, []));

else
  error('scaleLinkParamDoppler:InvalidLinkParamType', ...
    'linkParamCell must be a struct, a struct array, or a cell array of structs.');
end

for iFrame = 1:numel(linkParamSeq)
  if ~isstruct(linkParamSeq{iFrame})
    error('scaleLinkParamDoppler:InvalidFrameType', ...
      'Each frame entry in linkParamCell must be a struct.');
  end
end
end


function linkParamOut = localRestoreLinkParamSeq(linkParamSeq, restoreInfo)
%LOCALRESTORELINKPARAMSEQ Restore the original container type and shape.

switch restoreInfo.containerType
  case 'cell'
    linkParamOut = reshape(linkParamSeq, restoreInfo.containerSize);
  case 'struct'
    linkParamOut = reshape([linkParamSeq{:}], restoreInfo.containerSize);
  otherwise
    error('scaleLinkParamDoppler:InvalidRestoreInfo', ...
      'Unsupported restore container type.');
end
end


function [numSat, numUser] = localCheckFrameConsistency(linkParamSeq)
%LOCALCHECKFRAMECONSISTENCY Validate required fields and frame sizes.

linkParamRef = linkParamSeq{1};
localCheckRequiredFields(linkParamRef);

numSat = size(linkParamRef.fdGeom, 1);
numUser = size(linkParamRef.fdGeom, 2);

for iFrame = 2:numel(linkParamSeq)
  linkParam = linkParamSeq{iFrame};
  localCheckRequiredFields(linkParam);

  if ~isequal(size(linkParam.fdGeom), [numSat, numUser])
    error('scaleLinkParamDoppler:FdGeomSizeMismatch', ...
      'All frames must have the same fdGeom size.');
  end

  if ~isequal(size(linkParam.ref.fdGeom), [1, numUser])
    error('scaleLinkParamDoppler:RefFdGeomSizeMismatch', ...
      'All frames must have ref.fdGeom with size 1xnumUser.');
  end

  if localHasField(linkParamRef, {'rangeRate'})
    if ~localHasField(linkParam, {'rangeRate'}) || ...
        ~isequal(size(linkParam.rangeRate), [numSat, numUser])
      error('scaleLinkParamDoppler:RangeRateSizeMismatch', ...
        'All frames must have rangeRate with size numSat x numUser.');
    end
  end

  if localHasField(linkParamRef, {'ref', 'rangeRate'})
    if ~localHasField(linkParam, {'ref', 'rangeRate'}) || ...
        ~isequal(size(linkParam.ref.rangeRate), [1, numUser])
      error('scaleLinkParamDoppler:RefRangeRateSizeMismatch', ...
        'All frames must have ref.rangeRate with size 1xnumUser.');
    end
  end
end
end


function localCheckRequiredFields(linkParam)
%LOCALCHECKREQUIREDFIELDS Ensure required Doppler fields exist.

if ~isfield(linkParam, 'fdGeom') || ~isnumeric(linkParam.fdGeom)
  error('scaleLinkParamDoppler:MissingFdGeom', ...
    'Each linkParam frame must contain numeric field fdGeom.');
end

if ~isfield(linkParam, 'ref') || ~isstruct(linkParam.ref)
  error('scaleLinkParamDoppler:MissingRef', ...
    'Each linkParam frame must contain struct field ref.');
end

if ~isfield(linkParam.ref, 'fdGeom') || ~isnumeric(linkParam.ref.fdGeom)
  error('scaleLinkParamDoppler:MissingRefFdGeom', ...
    'Each linkParam.ref must contain numeric field fdGeom.');
end
end


function tf = localHasField(s, fieldPath)
%LOCALHASFIELD Check whether a nested field path exists.

tf = true;
for iField = 1:numel(fieldPath)
  fieldName = fieldPath{iField};
  if ~isstruct(s) || ~isfield(s, fieldName)
    tf = false;
    return;
  end
  s = s.(fieldName);
end
end


function seriesMat = localCollectFieldSeries(linkParamSeq, fieldPath, fieldSize)
%LOCALCOLLECTFIELDSERIES Stack one field over frames as [numElem x numFrame].

numFrame = numel(linkParamSeq);
numElem = prod(fieldSize);
seriesMat = zeros(numElem, numFrame);

for iFrame = 1:numFrame
  fieldVal = localGetField(linkParamSeq{iFrame}, fieldPath);
  if ~isnumeric(fieldVal) || ~isequal(size(fieldVal), fieldSize)
    error('scaleLinkParamDoppler:FieldSizeMismatch', ...
      'Field %s must have size %s in every frame.', ...
      strjoin(fieldPath, '.'), mat2str(fieldSize));
  end
  seriesMat(:, iFrame) = fieldVal(:);
end
end


function fieldVal = localGetField(s, fieldPath)
%LOCALGETFIELD Read one nested field path.

fieldVal = s;
for iField = 1:numel(fieldPath)
  fieldVal = fieldVal.(fieldPath{iField});
end
end


function anchorMat = localAffineAnchor(seriesMat, timeOffsetSec, timeRefSec)
%LOCALAFFINEANCHOR Get affine-fit value at the reference time.

[anchorMat, ~] = localFitAffineSeries(seriesMat, timeOffsetSec, timeRefSec);
end


function [anchorMat, slopeMat] = localFitAffineSeries(seriesMat, timeOffsetSec, timeRefSec)
%LOCALFITAFFINESERIES Fit y(t) ~= anchor + slope * (t - timeRefSec).

numFrame = numel(timeOffsetSec);
if numFrame <= 1 || max(timeOffsetSec) - min(timeOffsetSec) <= 0
  anchorMat = seriesMat(:, 1);
  slopeMat = zeros(size(anchorMat));
  return;
end

timeShift = reshape(timeOffsetSec - timeRefSec, [], 1);
designMat = [ones(numFrame, 1), timeShift];
coef = designMat \ seriesMat.';     % 2 x numSeries
anchorMat = coef(1, :).';
slopeMat = coef(2, :).';
end
