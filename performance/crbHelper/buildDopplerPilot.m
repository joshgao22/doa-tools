function [pilotFd, timeSec] = buildDopplerPilot(pilotWave, fd, timeInput)
%BUILDDOPPLERPILOT Apply Doppler shift to pilot waveform.
%
%Syntax:
%  pilotFd = buildDopplerPilot(pilotWave, fd, sampleRate)
%  pilotFd = buildDopplerPilot(pilotWave, fd, timeSec)
%  [pilotFd, timeSec] = buildDopplerPilot(...)
%
%Inputs:
%  pilotWave - row-wise baseband pilot waveform, size Nw x Ns
%              each row corresponds to one waveform sequence
%
%  fd        - Doppler frequency in Hz
%              - scalar : common Doppler applied to all rows
%              - vector : one Doppler per row, length Nw
%
%  timeInput - time specification
%              - scalar : sample rate in Hz, time starts from zero
%              - vector : sample time axis in seconds, length Ns
%
%Outputs:
%  pilotFd   - Doppler-shifted waveform, size Nw x Ns
%
%  timeSec   - row time axis in seconds, size 1 x Ns
%
%Description:
%  Builds the Doppler-shifted waveform
%
%    pilotFd = pilotWave .* exp(1j * 2 * pi * fd * t)
%
%  using either a sample-rate input or an explicit time axis.
%
%Notes:
%  - The Doppler phase is applied row-wise using implicit expansion.
%  - When timeInput is a scalar sample rate, the generated time axis is
%    timeSec = (0:Ns-1) / sampleRate.
%  - This function is the non-derivative counterpart of
%    buildDopplerPilotJacobian.
%
%Example:
%  [pilotFd, timeSec] = buildDopplerPilot(pilotPad, 250, sampleRate);
%
%See also:
%  buildDopplerPilotJacobian, genPilotWaveform

arguments
  pilotWave {mustBeNumeric, mustBeFinite}
  fd {mustBeNumeric, mustBeFinite}
  timeInput {mustBeNumeric, mustBeFinite}
end

if ~ismatrix(pilotWave) || isempty(pilotWave)
  error('buildDopplerPilot:InvalidPilotWave', ...
    'pilotWave must be a non-empty 2D matrix.');
end

[numWave, numSample] = size(pilotWave);
fdVec = localExpandFd(fd, numWave);
timeSec = localParseTimeInput(timeInput, numSample);

phaseArg = 2 * pi * (fdVec * timeSec);
phaseTerm = exp(1j * phaseArg);

pilotFd = pilotWave .* phaseTerm;

end

function fdVec = localExpandFd(fd, numWave)
%LOCALEXPANDFD Expand Doppler input to a column vector.

if isscalar(fd)
  fdVec = repmat(fd, numWave, 1);
else
  if ~isvector(fd) || numel(fd) ~= numWave
    error('buildDopplerPilot:InvalidFdSize', ...
      'fd must be a scalar or a vector with one entry per waveform row.');
  end
  fdVec = reshape(fd, [], 1);
end

end

function timeSec = localParseTimeInput(timeInput, numSample)
%LOCALPARSETIMEINPUT Convert time specification to a row vector in seconds.

if isscalar(timeInput)
  sampleRate = timeInput;
  if ~isfinite(sampleRate) || sampleRate <= 0
    error('buildDopplerPilot:InvalidSampleRate', ...
      'sampleRate must be a positive finite scalar.');
  end

  timeSec = (0:numSample-1) / sampleRate;

else
  if ~isvector(timeInput) || numel(timeInput) ~= numSample
    error('buildDopplerPilot:InvalidTimeAxis', ...
      'timeSec must be a vector whose length matches the number of samples.');
  end

  timeSec = reshape(timeInput, 1, []);
end

end
