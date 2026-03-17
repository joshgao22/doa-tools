function [pilotWave, waveInfo] = genPilotWaveform(pilotSym, symRate, osf, pulseShape, pulseOpt)
%GENPILOTWAVEFORM Generate oversampled pilot waveform from pilot symbols.
% Converts row-wise pilot symbols to waveform samples by upsampling and
% pulse shaping. The output is suitable for fractional-delay simulation and
% subsequent Doppler injection.
%
%Syntax:
%  pilotWave = genPilotWaveform(pilotSym, symRate)
%  pilotWave = genPilotWaveform(pilotSym, symRate, osf)
%  pilotWave = genPilotWaveform(pilotSym, symRate, osf, pulseShape)
%  [pilotWave, waveInfo] = genPilotWaveform(pilotSym, symRate, osf, pulseShape, pulseOpt)
%
%Inputs:
%  pilotSym   - pilot symbol matrix, size Nu x Ns
%               each row corresponds to one user
%
%  symRate    - symbol rate in Hz
%
%  osf        - oversampling factor (samples per symbol)
%               default: 4
%
%  pulseShape - pulse shaping type
%               - 'rect'    : rectangular pulse
%               - 'rrc'     : root raised cosine pulse
%               - 'impulse' : zero-insertion only, no pulse shaping
%               default: 'rect'
%
%  pulseOpt   - optional struct with fields depending on pulseShape
%
%               Common fields:
%               .matchInputPower   : whether to match row average power of
%                                    waveform to input symbol power
%                                    default: true
%               .removeFilterDelay : whether to remove symmetric FIR group
%                                    delay after pulse shaping
%                                    default: true
%
%               For 'rect':
%               .span              : pulse span in symbols
%                                    default: 1
%
%               For 'rrc':
%               .rolloff           : roll-off factor in [0, 1]
%                                    default: 0.25
%               .span              : filter span in symbols
%                                    default: 8
%
%Outputs:
%  pilotWave  - oversampled waveform, size Nu x (Ns * osf)
%
%  waveInfo   - auxiliary waveform information with fields:
%               .numUser
%               .numSym
%               .numSample
%               .symRate
%               .sampleRate
%               .osf
%               .pulseShape
%               .pulse
%               .pulseLen
%               .groupDelay
%               .matchInputPower
%               .removeFilterDelay
%               .symbolPower
%               .wavePower
%
%Notes:
%  - pilotSym is interpreted as symbol-rate input, not waveform samples.
%  - The waveform sample rate is sampleRate = symRate * osf.
%  - The pulse is normalized to unit energy before filtering.
%  - If matchInputPower=true, each row is rescaled so that its average
%    waveform power matches the corresponding input symbol power.
%  - The output length is kept equal to numSym * osf when
%    removeFilterDelay=true.
%
%Example - Rectangular shaping:
%  [pilotSym, pilotInfo] = genPilotSymbol(3, 128, 'qpsk', 1);
%  
%  [pilotWave, waveInfo] = genPilotWaveform( ...
%    pilotSym, 1e6, 4, 'rect');
%
%Example - Raised-cosine shaping:
%  [pilotSym, pilotInfo] = genPilotSymbol(3, 128, 'zadoffchu', 1);
%  
%  pulseOpt = struct();
%  pulseOpt.rolloff = 0.25;
%  pulseOpt.span = 8;
%  
%  [pilotWave, waveInfo] = genPilotWaveform( ...
%    pilotSym, 1e6, 8, 'rrc', pulseOpt);
%
%See also:
%  genPilotSignal, applyFracDelay

arguments
  pilotSym {mustBeNumeric, mustBeFinite}
  symRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  osf (1,1) {mustBePositive, mustBeInteger} = 4
  pulseShape = 'rect'
  pulseOpt (1,1) struct = struct()
end

if ~ismatrix(pilotSym) || isempty(pilotSym)
  error('genPilotWaveform:InvalidPilotSym', ...
    'pilotSym must be a non-empty 2D matrix.');
end

[numUser, numSym] = size(pilotSym);
pulseShape = localParsePulseShape(pulseShape);

matchInputPower = true;
if isfield(pulseOpt, 'matchInputPower') && ~isempty(pulseOpt.matchInputPower)
  matchInputPower = logical(pulseOpt.matchInputPower);
end

removeFilterDelay = true;
if isfield(pulseOpt, 'removeFilterDelay') && ~isempty(pulseOpt.removeFilterDelay)
  removeFilterDelay = logical(pulseOpt.removeFilterDelay);
end

% -------------------------------------------------------------------------
% Build pulse shaping filter
% -------------------------------------------------------------------------
[pulse, pulseMeta] = localBuildPulse(pulseShape, osf, pulseOpt);
pulseLen = numel(pulse);
groupDelay = floor((pulseLen - 1) / 2);

% -------------------------------------------------------------------------
% Upsample symbol sequence
% -------------------------------------------------------------------------
pilotUp = zeros(numUser, numSym * osf, 'like', pilotSym);
pilotUp(:, 1:osf:end) = pilotSym;

% -------------------------------------------------------------------------
% Pulse shaping
% -------------------------------------------------------------------------
if strcmp(pulseShape, 'impulse')
  pilotWave = pilotUp;
else
  fullLen = size(pilotUp, 2) + pulseLen - 1;
  pilotWaveFull = zeros(numUser, fullLen, 'like', pilotSym);

  for iUser = 1:numUser
    pilotWaveFull(iUser, :) = conv(pilotUp(iUser, :), pulse, 'full');
  end

  if removeFilterDelay
    startIdx = groupDelay + 1;
    stopIdx = startIdx + size(pilotUp, 2) - 1;
    pilotWave = pilotWaveFull(:, startIdx:stopIdx);
  else
    pilotWave = pilotWaveFull;
  end
end

% -------------------------------------------------------------------------
% Match row average power to input symbol power
% -------------------------------------------------------------------------
symbolPower = mean(abs(pilotSym).^2, 2);
wavePower = mean(abs(pilotWave).^2, 2);

if any(wavePower <= 0)
  error('genPilotWaveform:ZeroWavePower', ...
    'Generated waveform contains zero-power rows.');
end

if matchInputPower
  scaleVec = sqrt(symbolPower ./ wavePower);
  pilotWave = pilotWave .* scaleVec;
  wavePower = mean(abs(pilotWave).^2, 2);
end

% -------------------------------------------------------------------------
% Pack outputs
% -------------------------------------------------------------------------
waveInfo = struct();
waveInfo.numUser = numUser;
waveInfo.numSym = numSym;
waveInfo.numSample = size(pilotWave, 2);
waveInfo.symRate = symRate;
waveInfo.sampleRate = symRate * osf;
waveInfo.osf = osf;
waveInfo.pulseShape = pulseShape;
waveInfo.pulse = pulse;
waveInfo.pulseLen = pulseLen;
waveInfo.groupDelay = groupDelay;
waveInfo.matchInputPower = matchInputPower;
waveInfo.removeFilterDelay = removeFilterDelay;
waveInfo.symbolPower = symbolPower;
waveInfo.wavePower = wavePower;
waveInfo.pulseMeta = pulseMeta;
end

function pulseShape = localParsePulseShape(pulseShape)
%LOCALPARSEPULSESHAPE Normalize pulse shape input.

if isstring(pulseShape)
  if ~isscalar(pulseShape)
    error('genPilotWaveform:InvalidPulseShape', ...
      'pulseShape must be a character vector or string scalar.');
  end
  pulseShape = char(pulseShape);
end

if ~ischar(pulseShape)
  error('genPilotWaveform:InvalidPulseShape', ...
    'pulseShape must be a character vector or string scalar.');
end

pulseShape = lower(strtrim(pulseShape));

validShape = {'rect', 'rrc', 'impulse'};
if ~ismember(pulseShape, validShape)
  error('genPilotWaveform:UnsupportedPulseShape', ...
    'Unsupported pulseShape: %s.', pulseShape);
end
end

function [pulse, pulseMeta] = localBuildPulse(pulseShape, osf, pulseOpt)
%LOCALBUILDPULSE Build pulse shaping FIR.

switch pulseShape
  case 'impulse'
    pulse = 1;
    pulseMeta = struct();
    pulseMeta.type = 'impulse';

  case 'rect'
    span = 1;
    if isfield(pulseOpt, 'span') && ~isempty(pulseOpt.span)
      span = pulseOpt.span;
    end

    if ~isscalar(span) || ~isfinite(span) || span <= 0 || mod(span, 1) ~= 0
      error('genPilotWaveform:InvalidRectSpan', ...
        'pulseOpt.span for rect pulse must be a positive integer scalar.');
    end

    pulse = ones(1, span * osf);
    pulse = pulse / norm(pulse);

    pulseMeta = struct();
    pulseMeta.type = 'rect';
    pulseMeta.span = span;

  case 'rrc'
    rolloff = 0.25;
    if isfield(pulseOpt, 'rolloff') && ~isempty(pulseOpt.rolloff)
      rolloff = pulseOpt.rolloff;
    end

    span = 8;
    if isfield(pulseOpt, 'span') && ~isempty(pulseOpt.span)
      span = pulseOpt.span;
    end

    if ~isscalar(rolloff) || ~isfinite(rolloff) || rolloff < 0 || rolloff > 1
      error('genPilotWaveform:InvalidRolloff', ...
        'pulseOpt.rolloff must be a scalar in [0, 1].');
    end

    if ~isscalar(span) || ~isfinite(span) || span <= 0 || mod(span, 1) ~= 0
      error('genPilotWaveform:InvalidRrcSpan', ...
        'pulseOpt.span for rrc pulse must be a positive integer scalar.');
    end

    pulse = localGenRrcPulse(osf, rolloff, span);
    pulse = pulse / norm(pulse);

    pulseMeta = struct();
    pulseMeta.type = 'rrc';
    pulseMeta.rolloff = rolloff;
    pulseMeta.span = span;

  otherwise
    error('genPilotWaveform:UnsupportedPulseShape', ...
      'Unsupported pulseShape: %s.', pulseShape);
end
end

function pulse = localGenRrcPulse(osf, rolloff, span)
%LOCALGENRRCPULSE Generate root raised cosine pulse.
% Time axis is in symbol durations with T = 1.

t = (-span/2 : 1/osf : span/2);
pulse = zeros(size(t));

if rolloff == 0
  pulse = sinc(t);
  return;
end

for ii = 1:numel(t)
  ti = t(ii);

  if abs(ti) < 1e-12
    pulse(ii) = 1 + rolloff * (4/pi - 1);

  elseif abs(abs(ti) - 1/(4*rolloff)) < 1e-12
    pulse(ii) = (rolloff / sqrt(2)) * ...
      ((1 + 2/pi) * sin(pi / (4*rolloff)) + ...
       (1 - 2/pi) * cos(pi / (4*rolloff)));

  else
    numerator = sin(pi * ti * (1 - rolloff)) + ...
      4 * rolloff * ti * cos(pi * ti * (1 + rolloff));
    denominator = pi * ti * (1 - (4 * rolloff * ti)^2);
    pulse(ii) = numerator / denominator;
  end
end
end