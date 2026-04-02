function delayedSig = applyFracDelay(sig, delaySamples, filterHalfLen, windowType)
%APPLYFRACDELAY Apply row-wise fractional delay to signals.
% Applies a fractional sample delay to each row of the input signal matrix
% using windowed sinc interpolation. Positive delay means later arrival,
% i.e., the waveform is shifted to the right in time.
%
%Syntax:
%  delayedSig = applyFracDelay(sig, delaySamples)
%  delayedSig = applyFracDelay(sig, delaySamples, filterHalfLen)
%  delayedSig = applyFracDelay(sig, delaySamples, filterHalfLen, windowType)
%
%Inputs:
%  sig           - input signal matrix, size Ns x N
%                  each row is one signal sequence
%
%  delaySamples  - delay in samples
%                  - scalar     : shared by all rows
%                  - Ns x 1     : one delay per row
%
%  filterHalfLen - half length of sinc interpolation kernel
%                  actual FIR length is 2*filterHalfLen + 1
%                  default: 8
%
%  windowType    - window type for sinc truncation
%                  - 'rect'
%                  - 'hann'
%                  - 'hamming'
%                  - 'blackman'
%                  default: 'hann'
%
%Output:
%  delayedSig    - delayed signal matrix, same size as sig
%
%Definitions:
%  Let the total delay be
%      d = dInt + dFrac
%  where
%      dInt  = floor(d)
%      dFrac = d - dInt,   0 <= dFrac < 1
%
%  The function first applies the fractional part dFrac using a windowed
%  sinc kernel, and then applies the integer part dInt using zero-padded
%  sample shifting.
%
%Notes:
%  - Positive delay corresponds to x[n - d].
%  - The output length is kept equal to the input length.
%  - Boundary samples are zero-extended outside the observation interval.
%  - This function is suitable for pilot / training sequence simulation
%    before Doppler modulation is injected.
%
%See also:
%  conv, sinc

arguments
  sig {mustBeNumeric, mustBeFinite}
  delaySamples {mustBeNumeric, mustBeFinite}
  filterHalfLen (1,1) {mustBePositive, mustBeInteger} = 8
  windowType = 'hann'
end

if ~ismatrix(sig)
  error('applyFracDelay:InvalidSignalSize', ...
    'sig must be a 2D matrix.');
end

[numSig, numSnap] = size(sig);
delayVec = localParseDelay(delaySamples, numSig);
windowType = localParseWindowType(windowType);

delayedSig = zeros(numSig, numSnap, 'like', sig);

for iSig = 1:numSig
  currentSig = sig(iSig, :);
  currentDelay = delayVec(iSig);

  intDelay = floor(currentDelay);
  fracDelay = currentDelay - intDelay;

  if abs(fracDelay) < 1e-12
    fracSig = currentSig;
  else
    fracSig = localApplyFracKernel(currentSig, fracDelay, filterHalfLen, windowType);
  end

  delayedSig(iSig, :) = localApplyIntegerShift(fracSig, intDelay);
end
end

function delayVec = localParseDelay(delaySamples, numSig)
%LOCALPARSEDELAY Convert delay input to an Ns x 1 vector.

if isscalar(delaySamples)
  delayVec = repmat(delaySamples, numSig, 1);
  return;
end

delaySamples = delaySamples(:);
if numel(delaySamples) ~= numSig
  error('applyFracDelay:DelaySizeMismatch', ...
    'delaySamples must be a scalar or a vector of length size(sig,1).');
end

delayVec = delaySamples;
end

function windowType = localParseWindowType(windowType)
%LOCALPARSEWINDOWTYPE Normalize window type input.

if isstring(windowType)
  if ~isscalar(windowType)
    error('applyFracDelay:InvalidWindowType', ...
      'windowType must be a character vector or string scalar.');
  end
  windowType = char(windowType);
end

if ~ischar(windowType)
  error('applyFracDelay:InvalidWindowType', ...
    'windowType must be a character vector or string scalar.');
end

windowType = lower(strtrim(windowType));

validType = {'rect', 'hann', 'hamming', 'blackman'};
if ~ismember(windowType, validType)
  error('applyFracDelay:UnsupportedWindowType', ...
    'Unsupported windowType: %s.', windowType);
end
end

function fracSig = localApplyFracKernel(sigRow, fracDelay, filterHalfLen, windowType)
%LOCALAPPLYFRACKERNEL Apply fractional delay in [0,1) using sinc interpolation.

tapIndex = -filterHalfLen:filterHalfLen;
kernel = sinc(tapIndex - fracDelay);
window = localGenWindow(numel(tapIndex), windowType);

kernel = kernel .* window;
kernel = kernel / sum(kernel);

fracSig = conv(sigRow, kernel, 'same');
end

function shiftedSig = localApplyIntegerShift(sigRow, intDelay)
%LOCALAPPLYINTEGERSHIFT Apply zero-padded integer sample shift.

numSnap = numel(sigRow);
shiftedSig = zeros(1, numSnap, 'like', sigRow);

if intDelay >= 0
  if intDelay < numSnap
    shiftedSig((intDelay + 1):end) = sigRow(1:(end - intDelay));
  end
else
  advance = -intDelay;
  if advance < numSnap
    shiftedSig(1:(end - advance)) = sigRow((advance + 1):end);
  end
end
end

function window = localGenWindow(winLen, windowType)
%LOCALGENWINDOW Generate one window row vector.

if winLen == 1
  window = 1;
  return;
end

n = 0:winLen-1;

switch windowType
  case 'rect'
    window = ones(1, winLen);

  case 'hann'
    window = 0.5 - 0.5 * cos(2*pi*n/(winLen - 1));

  case 'hamming'
    window = 0.54 - 0.46 * cos(2*pi*n/(winLen - 1));

  case 'blackman'
    window = 0.42 ...
      - 0.5 * cos(2*pi*n/(winLen - 1)) ...
      + 0.08 * cos(4*pi*n/(winLen - 1));

  otherwise
    error('applyFracDelay:UnsupportedWindowType', ...
      'Unsupported windowType: %s.', windowType);
end
end