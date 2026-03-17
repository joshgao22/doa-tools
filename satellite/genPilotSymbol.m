function [pilotSym, pilotInfo] = genPilotSymbol(numUser, numSnap, pilotType, pilotPower, pilotOpt)
%GENPILOTSYMBOL Generate user pilot / training symbols.
% Generates one pilot sequence per user in row-wise form, suitable for
% multi-user and multi-satellite signal simulation.
%
%Syntax:
%  pilotSym = genPilotSymbol(numUser, numSnap)
%  pilotSym = genPilotSymbol(numUser, numSnap, pilotType)
%  pilotSym = genPilotSymbol(numUser, numSnap, pilotType, pilotPower)
%  [pilotSym, pilotInfo] = genPilotSymbol(numUser, numSnap, pilotType, pilotPower, pilotOpt)
%
%Inputs:
%  numUser    - number of users
%
%  numSnap    - number of time samples
%
%  pilotType  - pilot type
%               - 'ones'       : all-one pilot
%               - 'bpsk'       : random BPSK pilot
%               - 'pn'         : alias of 'bpsk'
%               - 'qpsk'       : random QPSK pilot
%               - 'zadoffchu'  : Zadoff-Chu pilot
%               - 'zc'         : alias of 'zadoffchu'
%               - 'custom'     : user-provided pilot sequence
%               default: 'ones'
%
%  pilotPower - pilot power specification
%               - scalar       : shared power for all users
%               - vector Nu x 1: per-user power
%               default: 1
%
%  pilotOpt   - optional struct with fields depending on pilotType
%
%               Common fields:
%               .seed          : RNG seed used for random pilots
%
%               For 'zadoffchu':
%               .rootIdx       : scalar or Nu x 1 root indices
%               .cyclicShift   : scalar or Nu x 1 cyclic shifts
%
%               For 'custom':
%               .sequence      : custom sequence
%                                - 1 x N / N x 1 : shared by all users
%                                - Nu x N        : one row per user
%                                - scalar        : constant pilot
%
%Outputs:
%  pilotSig   - Nu x N pilot matrix, one row per user
%
%  pilotInfo  - auxiliary information with fields:
%               .type
%               .numUser
%               .numSnap
%               .pilotPower
%               .avgPower
%               .seed
%               .rootIdx
%               .cyclicShift
%
%Notes:
%  - The generated pilot is normalized to unit average power before applying
%    the requested pilotPower scaling.
%  - For 'zadoffchu', each user can use a different root and cyclic shift.
%  - Output pilotSig is directly compatible with genMultiSatSnapshots.
%
%See also:
%  applyFracDelay

arguments
  numUser (1,1) {mustBePositive, mustBeInteger}
  numSnap (1,1) {mustBePositive, mustBeInteger}
  pilotType = 'ones'
  pilotPower {mustBeNumeric, mustBeNonnegative} = 1
  pilotOpt (1,1) struct = struct()
end

pilotType = localParsePilotType(pilotType);
powerVec = localParsePilotPower(pilotPower, numUser);

seed = [];
if isfield(pilotOpt, 'seed') && ~isempty(pilotOpt.seed)
  seed = pilotOpt.seed;
end

restoreRng = false;
if ~isempty(seed)
  oldRngState = rng;
  rng(seed);
  restoreRng = true;
end

switch pilotType
  case 'ones'
    pilotSym = ones(numUser, numSnap);

  case {'bpsk', 'pn'}
    symbolIdx = randi([0 1], numUser, numSnap);
    pilotSym = 2 * symbolIdx - 1;

  case 'qpsk'
    realPart = 2 * randi([0 1], numUser, numSnap) - 1;
    imagPart = 2 * randi([0 1], numUser, numSnap) - 1;
    pilotSym = (realPart + 1j * imagPart) / sqrt(2);

  case {'zadoffchu', 'zc'}
    rootIdx = localParseZcRootIdx(pilotOpt, numUser, numSnap);
    cyclicShift = localParseZcShift(pilotOpt, numUser);
    pilotSym = zeros(numUser, numSnap);

    for iUser = 1:numUser
      baseSeq = localGenZcSequence(numSnap, rootIdx(iUser));
      if cyclicShift(iUser) ~= 0
        baseSeq = circshift(baseSeq, [0, cyclicShift(iUser)]);
      end
      pilotSym(iUser, :) = baseSeq;
    end

  case 'custom'
    if ~isfield(pilotOpt, 'sequence') || isempty(pilotOpt.sequence)
      error('genPilotSignal:MissingCustomSequence', ...
        'pilotOpt.sequence is required for pilotType="custom".');
    end
    pilotSym = localParseCustomSequence(pilotOpt.sequence, numUser, numSnap);

  otherwise
    error('genPilotSignal:UnsupportedPilotType', ...
      'Unsupported pilotType: %s.', pilotType);
end

if restoreRng
  rng(oldRngState);
end

% -------------------------------------------------------------------------
% Normalize row power and apply requested pilot power
% -------------------------------------------------------------------------
rowPower = mean(abs(pilotSym).^2, 2);
if any(rowPower <= 0)
  error('genPilotSignal:InvalidPilotPower', ...
    'Generated pilot contains zero-power rows.');
end

pilotSym = pilotSym ./ sqrt(rowPower);
pilotSym = pilotSym .* sqrt(powerVec);

pilotInfo = struct();
pilotInfo.type = pilotType;
pilotInfo.numUser = numUser;
pilotInfo.numSnap = numSnap;
pilotInfo.pilotPower = powerVec;
pilotInfo.avgPower = mean(abs(pilotSym).^2, 2);
pilotInfo.seed = seed;

if ismember(pilotType, {'zadoffchu', 'zc'})
  pilotInfo.rootIdx = rootIdx;
  pilotInfo.cyclicShift = cyclicShift;
else
  pilotInfo.rootIdx = [];
  pilotInfo.cyclicShift = [];
end
end

function pilotType = localParsePilotType(pilotType)
%LOCALPARSEPILOTTYPE Normalize pilot type input.

if isstring(pilotType)
  if ~isscalar(pilotType)
    error('genPilotSignal:InvalidPilotType', ...
      'pilotType must be a character vector or string scalar.');
  end
  pilotType = char(pilotType);
end

if ~ischar(pilotType)
  error('genPilotSignal:InvalidPilotType', ...
    'pilotType must be a character vector or string scalar.');
end

pilotType = lower(strtrim(pilotType));
end

function powerVec = localParsePilotPower(pilotPower, numUser)
%LOCALPARSEPILOTPOWER Convert pilot power input to an Nu x 1 vector.

if isempty(pilotPower)
  powerVec = ones(numUser, 1);
  return;
end

if isscalar(pilotPower)
  powerVec = repmat(pilotPower, numUser, 1);
  return;
end

pilotPower = pilotPower(:);
if numel(pilotPower) ~= numUser
  error('genPilotSignal:InvalidPilotPowerSize', ...
    'pilotPower must be a scalar or a vector of length numUser.');
end

powerVec = pilotPower;
end

function rootIdx = localParseZcRootIdx(pilotOpt, numUser, seqLen)
%LOCALPARSEZCROOTIDX Parse Zadoff-Chu root indices.

if isfield(pilotOpt, 'rootIdx') && ~isempty(pilotOpt.rootIdx)
  rootIdx = pilotOpt.rootIdx(:);
else
  rootIdx = localDefaultZcRootIdx(numUser, seqLen);
end

if isscalar(rootIdx)
  rootIdx = repmat(rootIdx, numUser, 1);
end

if numel(rootIdx) ~= numUser
  error('genPilotSignal:InvalidZcRootIdxSize', ...
    'pilotOpt.rootIdx must be a scalar or a vector of length numUser.');
end

for iUser = 1:numUser
  if ~isfinite(rootIdx(iUser)) || mod(rootIdx(iUser), 1) ~= 0
    error('genPilotSignal:InvalidZcRootIdx', ...
      'Each Zadoff-Chu root index must be a finite integer.');
  end
  if gcd(abs(rootIdx(iUser)), seqLen) ~= 1
    error('genPilotSignal:ZcRootNotCoprime', ...
      'Each Zadoff-Chu root index must be coprime with numSnap.');
  end
end
end

function cyclicShift = localParseZcShift(pilotOpt, numUser)
%LOCALPARSEZCSHIFT Parse Zadoff-Chu cyclic shifts.

if isfield(pilotOpt, 'cyclicShift') && ~isempty(pilotOpt.cyclicShift)
  cyclicShift = pilotOpt.cyclicShift(:);
else
  cyclicShift = zeros(numUser, 1);
end

if isscalar(cyclicShift)
  cyclicShift = repmat(cyclicShift, numUser, 1);
end

if numel(cyclicShift) ~= numUser
  error('genPilotSignal:InvalidZcShiftSize', ...
    'pilotOpt.cyclicShift must be a scalar or a vector of length numUser.');
end

if any(~isfinite(cyclicShift)) || any(mod(cyclicShift, 1) ~= 0)
  error('genPilotSignal:InvalidZcShift', ...
    'Each cyclic shift must be a finite integer.');
end
end

function rootIdx = localDefaultZcRootIdx(numUser, seqLen)
%LOCALDEFAULTZCROOTIDX Generate default coprime roots.

candidate = [];
for rr = 1:max(seqLen - 1, 1)
  if gcd(rr, seqLen) == 1
    candidate(end+1, 1) = rr; %#ok<AGROW>
  end
end

if isempty(candidate)
  error('genPilotSignal:NoValidZcRoot', ...
    'No valid Zadoff-Chu root index exists for the given numSnap.');
end

selectIdx = 1 + mod((0:numUser-1), numel(candidate));
rootIdx = candidate(selectIdx);
rootIdx = rootIdx(:);
end

function seq = localGenZcSequence(seqLen, rootIdx)
%LOCALGENZCSEQUENCE Generate one Zadoff-Chu sequence row.

n = 0:seqLen-1;

if mod(seqLen, 2) == 0
  seq = exp(-1j * pi * rootIdx * (n.^2) / seqLen);
else
  seq = exp(-1j * pi * rootIdx * n .* (n + 1) / seqLen);
end
end

function pilotSig = localParseCustomSequence(sequence, numUser, numSnap)
%LOCALPARSECUSTOMSEQUENCE Parse custom pilot sequence input.

if ~isnumeric(sequence) || isempty(sequence)
  error('genPilotSignal:InvalidCustomSequence', ...
    'pilotOpt.sequence must be a non-empty numeric array.');
end

if isscalar(sequence)
  pilotSig = repmat(sequence, numUser, numSnap);
  return;
end

if isvector(sequence)
  sequence = reshape(sequence, 1, []);
  if numel(sequence) ~= numSnap
    error('genPilotSignal:CustomSequenceLengthMismatch', ...
      'Vector custom sequence must have length numSnap.');
  end
  pilotSig = repmat(sequence, numUser, 1);
  return;
end

if isequal(size(sequence), [numUser, numSnap])
  pilotSig = sequence;
  return;
end

if isequal(size(sequence), [1, numSnap])
  pilotSig = repmat(sequence, numUser, 1);
  return;
end

error('genPilotSignal:InvalidCustomSequenceSize', ...
  'pilotOpt.sequence must be scalar, 1xnumSnap, or numUser x numSnap.');
end