function [pilotFd, dPilotDfd, timeSec] = buildDopplerPilotJacobian(pilotWave, fd, timeInput)
%BUILDDOPPLERPILOTJACOBIAN Apply Doppler shift to pilot waveform and
%compute its derivative.
%
%Syntax:
%  [pilotFd, dPilotDfd] = buildDopplerPilotJacobian(pilotWave, fd, sampleRate)
%  [pilotFd, dPilotDfd] = buildDopplerPilotJacobian(pilotWave, fd, timeSec)
%  [pilotFd, dPilotDfd, timeSec] = buildDopplerPilotJacobian(...)
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
%  dPilotDfd - derivative of pilotFd with respect to Doppler frequency
%              - if fd is scalar, this is the derivative with respect to
%                the common Doppler
%              - if fd is a vector, row i is the derivative with respect to
%                fd(i)
%
%  timeSec   - row time axis in seconds, size 1 x Ns
%
%Description:
%  Builds the Doppler-shifted waveform
%
%    pilotFd = pilotWave .* exp(1j * 2 * pi * fd * t)
%
%  and its first-order derivative with respect to Doppler frequency
%
%    dPilotDfd = pilotFd .* (1j * 2 * pi * t).
%
%Notes:
%  - The Doppler phase is applied row-wise using implicit expansion.
%  - When timeInput is a scalar sample rate, the generated time axis is
%    timeSec = (0:Ns-1) / sampleRate.
%  - This function is intended for Doppler-aware pilot modeling in joint
%    DoA-Doppler estimation and CRB derivations.
%
%Example:
%  [pilotFd, dPilotDfd, timeSec] = buildDopplerPilotJacobian( ...
%    pilotPad, 250, sampleRate);
%
%See also:
%  genPilotWaveform

arguments
  pilotWave {mustBeNumeric, mustBeFinite}
  fd {mustBeNumeric, mustBeFinite}
  timeInput {mustBeNumeric, mustBeFinite}
end

if ~ismatrix(pilotWave) || isempty(pilotWave)
  error('buildDopplerPilotJacobian:InvalidPilotWave', ...
    'pilotWave must be a non-empty 2D matrix.');
end

[numWave, numSample] = size(pilotWave);
fdVec = localExpandFd(fd, numWave);
timeSec = localParseTimeInput(timeInput, numSample);

phaseArg = 2 * pi * (fdVec * timeSec);
phaseTerm = exp(1j * phaseArg);

pilotFd = pilotWave .* phaseTerm;
dPilotDfd = pilotFd .* (1j * 2 * pi * timeSec);

end

function fdVec = localExpandFd(fd, numWave)
%LOCALEXPANDFD Expand Doppler input to a column vector.

if isscalar(fd)
  fdVec = repmat(fd, numWave, 1);
else
  if ~isvector(fd) || numel(fd) ~= numWave
    error('buildDopplerPilotJacobian:InvalidFdSize', ...
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
    error('buildDopplerPilotJacobian:InvalidSampleRate', ...
      'sampleRate must be a positive finite scalar.');
  end

  timeSec = (0:numSample-1) / sampleRate;

else
  if ~isvector(timeInput) || numel(timeInput) ~= numSample
    error('buildDopplerPilotJacobian:InvalidTimeAxis', ...
      'timeSec must be a vector whose length matches the number of samples.');
  end

  timeSec = reshape(timeInput, 1, []);
end

end
