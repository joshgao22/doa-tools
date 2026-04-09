function [estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotOpt(varargin)
%ESTIMATORDOADOPPLERMLEPILOTOPT Compatibility wrapper for the SF solver.
% This legacy entry is retained only for backward compatibility.
% The formal technical route uses estimatorDoaDopplerMlePilotSfOpt as the
% single-frame static DoA-Doppler baseline, and estimatorDoaDopplerMlePilotMfOpt
% as the unified multi-frame solver.
%
% See also: estimatorDoaDopplerMlePilotSfOpt, estimatorDoaDopplerMlePilotMfOpt

[estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotSfOpt(varargin{:});
end
