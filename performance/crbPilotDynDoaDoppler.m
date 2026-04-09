function [crb, aux] = crbPilotDynDoaDoppler(varargin)
%CRBPILOTDYNDOADOPPLER Compatibility wrapper for the MF pilot CRB.
% This legacy entry is retained only for backward compatibility.
% The formal multi-frame CRB route is consolidated in crbPilotMfDoaDoppler.
%
% See also: crbPilotMfDoaDoppler

[crb, aux] = crbPilotMfDoaDoppler(varargin{:});
end
