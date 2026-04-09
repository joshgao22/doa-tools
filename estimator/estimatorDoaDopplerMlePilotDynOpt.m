function [estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotDynOpt(varargin)
%ESTIMATORDOADOPPLERMLEPILOTDYNOPT Compatibility wrapper for the MF solver.
% This legacy entry is retained only for backward compatibility.
% The formal technical route is consolidated in estimatorDoaDopplerMlePilotMfOpt,
% which unifies continuous, relaxed, and independent multi-frame likelihoods.
%
% See also: estimatorDoaDopplerMlePilotMfOpt

[estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotMfOpt(varargin{:});
end
