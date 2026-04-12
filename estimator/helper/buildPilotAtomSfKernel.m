function [atomMat, atomVec, zVal, etaVal] = buildPilotAtomSfKernel(array, ...
  wavelength, localDoa, pilotVec, timeSec, fdHz, yMat)
%BUILDPILOTATOMSFKERNEL Build the single-frame static pilot atom.
%
% Builds the same pilot atom used by the single-frame static
% DoA-Doppler estimator for one source and one satellite block.
%
%Syntax:
%  [atomMat, atomVec] = buildPilotAtomSfKernel(array, wavelength, ...
%    localDoa, pilotVec, timeSec, fdHz)
%
%  [atomMat, atomVec, zVal, etaVal] = buildPilotAtomSfKernel(array, ...
%    wavelength, localDoa, pilotVec, timeSec, fdHz, yMat)
%
%Inputs:
%  array          - array struct used by steeringMatrix
%  wavelength     - carrier wavelength in m
%  localDoa       - 2x1 local DoA parameter vector
%  pilotVec       - 1xT or Tx1 known pilot waveform
%  timeSec        - 1xT or Tx1 local sample times in s
%  fdHz           - scalar Doppler in Hz
%  yMat           - optional MxT received block
%
%Outputs:
%  atomMat        - MxT pilot atom matrix
%  atomVec        - vectorized atomMat
%  zVal           - matched-filter output atomVec' * yVec when yMat is given
%  etaVal         - atom energy real(atomVec' * atomVec)
%
%See also:
%  estimatorDoaDopplerMlePilotSfOpt, evalDoaDopplerDynProfileLike,
%  steeringMatrix

arguments
  array (1,1) struct
  wavelength (1,1) double {mustBeFinite, mustBePositive}
  localDoa (:,1) double {mustBeFinite}
  pilotVec {mustBeNumeric, mustBeFinite}
  timeSec {mustBeNumeric, mustBeFinite}
  fdHz (1,1) double {mustBeFinite}
  yMat {mustBeNumeric, mustBeFinite} = []
end

currentSteering = steeringMatrix(array, wavelength, localDoa);
currentSteering = currentSteering(:, 1);
sourceSig = reshape(pilotVec, 1, []) .* exp(1j * 2 * pi * fdHz * reshape(timeSec, 1, []));
atomMat = currentSteering * sourceSig;
atomVec = atomMat(:);
etaVal = real(atomVec' * atomVec);

if isempty(yMat)
  zVal = [];
  return;
end

yVec = yMat(:);
zVal = atomVec' * yVec;
end
