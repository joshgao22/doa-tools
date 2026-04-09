function [estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotMfStatOpt( ...
  sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, fdRange, initParam, verbose, modelOpt)
%ESTIMATORDOADOPPLERMLEPILOTMFSTATOPT Wrapper for the zero-rate MF baseline.
% The static multi-frame baseline is no longer maintained as a separate
% implementation. It is reduced to the unified multi-frame estimator with
% continuous phase and fdRateMode = 'zero'.
%
% See also: estimatorDoaDopplerMlePilotMfOpt

if nargin < 8 || isempty(initParam)
  initParam = [];
end
if nargin < 9 || isempty(verbose)
  verbose = false;
end
if nargin < 10 || isempty(modelOpt)
  modelOpt = struct();
end
modelOpt.phaseMode = 'continuous';
modelOpt.fdRateMode = 'zero';
[estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotMfOpt( ...
  sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, ...
  fdRange, [], initParam, verbose, modelOpt);
end
