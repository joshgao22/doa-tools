function rxSigSub = selectFrameRxSig(rxSigIn, frameIdx)
%SELECTFRAMERXSIG Extract one frame subset from multi-frame received signals.
% Clips the outer frame axis of one multi-frame rxSig container while
% preserving the per-frame payload shape. The function is intended to pair
% with selectFrameSceneSeq so dynamic dev / regression scripts can build
% periodic and non-periodic frame subsets from one common master window.
%
%Syntax:
%  rxSigSub = selectFrameRxSig(rxSigIn, frameIdx)
%
%Inputs:
%  rxSigIn   - multi-frame received-signal container
%              * 1xNf or Nfx1 cell array
%              * each cell can be one per-frame matrix or one per-frame
%                cell-of-satellites payload
%
%  frameIdx  - selected frame indices inside rxSigIn
%              * logical mask with one entry per input frame
%              * numeric subset indices in 1...Nf
%
%Output:
%  rxSigSub  - clipped frame subset, always returned as 1xNsub cell array
%
%See also:
%  selectFrameSceneSeq, selectSatSceneSeq, selectRxSigBySat

arguments
  rxSigIn
  frameIdx
end

if ~iscell(rxSigIn)
  error('selectFrameRxSig:InvalidInputType', ...
    'rxSigIn must be a multi-frame cell array.');
end

numFrameAll = numel(rxSigIn);
frameIdx = localValidateFrameIdx(frameIdx, numFrameAll);
rxSigSub = reshape(rxSigIn(frameIdx), 1, []);
end


function frameIdx = localValidateFrameIdx(frameIdx, numFrameAll)
%LOCALVALIDATEFRAMEIDX Validate and normalize one frame-index selector.

if islogical(frameIdx)
  if numel(frameIdx) ~= numFrameAll
    error('selectFrameRxSig:InvalidLogicalLength', ...
      'Logical frameIdx must contain exactly one entry per input frame.');
  end
  frameIdx = find(frameIdx);
else
  validateattributes(frameIdx, {'numeric'}, {'vector', 'integer', 'positive'}, ...
    mfilename, 'frameIdx');
  frameIdx = unique(frameIdx(:).', 'stable');
end

if isempty(frameIdx)
  error('selectFrameRxSig:EmptySelection', ...
    'At least one frame must be selected.');
end
if any(frameIdx < 1 | frameIdx > numFrameAll)
  error('selectFrameRxSig:FrameIdxOutOfRange', ...
    'frameIdx must stay inside 1...numel(rxSigIn).');
end
frameIdx = reshape(frameIdx, 1, []);
end
