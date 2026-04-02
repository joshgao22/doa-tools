function rxSigSub = selectRxSigBySat(rxSig, satIdx, containerType)
%SELECTRXSIGBYSAT Extract a satellite subset from snapshot containers.
% Clips satellite-dependent received-signal containers while preserving the
% frame grouping used by the estimators and test scripts.
%
%Syntax:
%  rxSigSub = selectRxSigBySat(rxSig, satIdx)
%
%  rxSigSub = selectRxSigBySat(rxSig, satIdx, containerType)
%
%Inputs:
%  rxSig            - received snapshot container. Supported forms:
%                     1) single-frame multi-satellite
%                        - 1xNs cell, one MxT block per satellite
%
%                     2) multi-frame
%                        - 1xNf cell, each frame is either
%                          * 1xNs cell, one MxT block per satellite
%                          * MxT matrix when Ns = 1
%
%                     3) single-satellite single-frame
%                        - MxT matrix
%
%  satIdx           - selected satellite indices
%                     - logical mask : one entry per input satellite
%                     - vector       : subset indices in 1...Ns
%
%  containerType    - optional container interpretation
%                     - 'auto'        : infer from rxSig (default)
%                     - 'singleFrame' : treat rxSig as one frame only
%                     - 'multiFrame'  : treat rxSig as a frame sequence
%
%Output:
%  rxSigSub         - clipped snapshot container
%                     - single-frame input:
%                         * matrix when one satellite is selected
%                         * 1xNsSel cell otherwise
%                     - multi-frame input:
%                         * 1xNf cell, one clipped frame per cell
%
%Auto inference:
%  - Nested cell input is interpreted as multi-frame.
%  - Numeric matrix input is interpreted as single-satellite single-frame.
%  - A flat cell array of numeric matrices is ambiguous when satIdx is the
%    scalar 1, because it may represent either
%      * one multi-satellite frame, or
%      * multiple single-satellite frames.
%    In that case, specify containerType explicitly.
%
%Notes:
%  - This function only clips by satellite. It does not clip frames or
%    modify signal values.
%  - Satellite indices are interpreted in the local ordering of the input
%    container.
%
%See also:
%  selectSatScene, selectSatSceneSeq, selectFrameSceneSeq

if nargin < 3 || isempty(containerType)
  containerType = 'auto';
end
containerType = localNormalizeContainerType(containerType);

if isnumeric(rxSig) || islogical(rxSig)
  satIdx = localNormalizeSatIdx(satIdx, 1, 'rxSig');
  if ~isequal(satIdx, 1)
    error('selectRxSigBySat:SingleSatOnly', ...
      'Numeric rxSig input represents one satellite only, so satIdx must be 1.');
  end
  rxSigSub = rxSig;
  return;
end

if ~iscell(rxSig)
  error('selectRxSigBySat:InvalidRxSigType', ...
    'rxSig must be a numeric matrix or a cell container.');
end

rxSig = reshape(rxSig, 1, []);
containerType = localResolveContainerType(rxSig, satIdx, containerType);

switch containerType
  case 'singleframe'
    rxSigSub = localSubsetSingleFrame(rxSig, satIdx, 'rxSig');

  case 'multiframe'
    numFrame = numel(rxSig);
    rxSigSub = cell(1, numFrame);
    for iFrame = 1:numFrame
      frameSig = rxSig{iFrame};
      frameName = sprintf('rxSig{%d}', iFrame);
      rxSigSub{iFrame} = localSubsetFrame(frameSig, satIdx, frameName);
    end

  otherwise
    error('selectRxSigBySat:InternalContainerType', ...
      'Unsupported containerType ''%s''.', containerType);
end
end


function containerType = localNormalizeContainerType(containerType)
%LOCALNORMALIZECONTAINERTYPE Normalize the container-type string.

if isstring(containerType)
  if numel(containerType) ~= 1
    error('selectRxSigBySat:InvalidContainerType', ...
      'containerType must be a character vector or a string scalar.');
  end
  containerType = char(containerType);
end

if ~ischar(containerType)
  error('selectRxSigBySat:InvalidContainerType', ...
    'containerType must be a character vector or a string scalar.');
end

containerType = lower(strtrim(containerType));
validType = {'auto', 'singleframe', 'multiframe'};
if ~ismember(containerType, validType)
  error('selectRxSigBySat:InvalidContainerType', ...
    'Unsupported containerType ''%s''.', containerType);
end
end


function containerType = localResolveContainerType(rxSig, satIdx, containerType)
%LOCALRESOLVECONTAINERTYPE Resolve the rxSig container interpretation.

if ~strcmp(containerType, 'auto')
  return;
end

isNestedCell = any(cellfun(@iscell, rxSig));
if isNestedCell
  containerType = 'multiframe';
  return;
end

isNumericCell = all(cellfun(@(x) isnumeric(x) || islogical(x), rxSig));
if ~isNumericCell
  error('selectRxSigBySat:UnsupportedCellContent', ...
    'Flat rxSig cell arrays must contain only numeric matrices.');
end

numEntry = numel(rxSig);
if numEntry <= 1
  containerType = 'singleframe';
  return;
end

satIdxFlat = satIdx(:).';
if isnumeric(satIdxFlat) && all(isfinite(satIdxFlat)) && any(satIdxFlat ~= 1)
  containerType = 'singleframe';
  return;
end

if islogical(satIdx)
  if numel(satIdx) ~= 1
    containerType = 'singleframe';
    return;
  end
end

error('selectRxSigBySat:AmbiguousFlatCell', ...
  ['Ambiguous flat rxSig cell array: it may represent a single multi-' ...
   'satellite frame or a multi-frame single-satellite sequence. ' ...
   'Specify containerType as ''singleFrame'' or ''multiFrame'' explicitly.']);
end


function rxSigSub = localSubsetSingleFrame(rxSig, satIdx, varName)
%LOCALSUBSETSINGLEFRAME Subset one single-frame rxSig container.

if iscell(rxSig)
  numSat = numel(rxSig);
  satIdx = localNormalizeSatIdx(satIdx, numSat, varName);
  if numel(satIdx) == 1
    rxSigSub = rxSig{satIdx};
  else
    rxSigSub = reshape(rxSig(satIdx), 1, []);
  end
  return;
end

satIdx = localNormalizeSatIdx(satIdx, 1, varName);
if ~isequal(satIdx, 1)
  error('selectRxSigBySat:SingleSatFrameOnly', ...
    '%s represents one satellite only, so satIdx must be 1.', varName);
end
rxSigSub = rxSig;
end


function frameSub = localSubsetFrame(frameSig, satIdx, varName)
%LOCALSUBSETFRAME Subset one frame inside a multi-frame container.

if iscell(frameSig)
  numSat = numel(frameSig);
  satIdx = localNormalizeSatIdx(satIdx, numSat, varName);
  if numel(satIdx) == 1
    frameSub = frameSig{satIdx};
  else
    frameSub = reshape(frameSig(satIdx), 1, []);
  end
  return;
end

if isnumeric(frameSig) || islogical(frameSig)
  satIdx = localNormalizeSatIdx(satIdx, 1, varName);
  if ~isequal(satIdx, 1)
    error('selectRxSigBySat:SingleSatFrameOnly', ...
      '%s represents one satellite only, so satIdx must be 1.', varName);
  end
  frameSub = frameSig;
  return;
end

error('selectRxSigBySat:UnsupportedFrameType', ...
  'Unsupported frame container type in %s.', varName);
end


function satIdx = localNormalizeSatIdx(satIdx, numSat, varName)
%LOCALNORMALIZESATIDX Validate and normalize satellite indices.

if islogical(satIdx)
  satIdx = reshape(satIdx, 1, []);
  if numel(satIdx) ~= numSat
    error('selectRxSigBySat:InvalidSatMask', ...
      '%s logical satIdx must have exactly %d entries.', varName, numSat);
  end
  satIdx = find(satIdx);
else
  validateattributes(satIdx, {'numeric'}, {'real', 'finite', 'vector'}, mfilename, 'satIdx');
  satIdx = reshape(satIdx, 1, []);
  if isempty(satIdx)
    error('selectRxSigBySat:EmptySatIdx', 'satIdx must not be empty.');
  end
  if any(satIdx < 1) || any(satIdx > numSat) || any(abs(satIdx - round(satIdx)) > 0)
    error('selectRxSigBySat:SatIdxOutOfRange', ...
      'satIdx must contain integer indices in the range 1...%d.', numSat);
  end
  satIdx = round(satIdx);
end

satIdx = unique(satIdx, 'stable');
if isempty(satIdx)
  error('selectRxSigBySat:EmptySelection', 'The selected satellite subset is empty.');
end
end
