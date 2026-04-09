function [crb, aux] = crbPilotMfStatDoaDoppler( ...
  sceneSeq, pilotWave, carrierFreq, sampleRate, doaParam, fdRef, fdRate, pathGain, noiseVar, modelOpt)
%CRBPILOTMFSTATDOADOPPLER Wrapper for the zero-rate MF CRB baseline.
% The static multi-frame baseline is reduced to the unified multi-frame CRB
% with continuous phase and fdRateMode = 'zero'.
%
% See also: crbPilotMfDoaDoppler

if nargin < 10 || isempty(modelOpt)
  modelOpt = struct();
end
modelOpt.phaseMode = 'continuous';
modelOpt.fdRateMode = 'zero';
[crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, ...
  sampleRate, doaParam, fdRef, fdRate, pathGain, noiseVar, modelOpt);
end
