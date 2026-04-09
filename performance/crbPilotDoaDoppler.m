function [crb, aux] = crbPilotDoaDoppler(varargin)
%CRBPILOTDOADOPPLER Compatibility wrapper for the SF pilot CRB.
% This legacy entry is retained only for backward compatibility.
% The formal static single-frame CRB route uses crbPilotSfDoaDoppler.
%
% See also: crbPilotSfDoaDoppler, crbPilotMfDoaDoppler

[crb, aux] = crbPilotSfDoaDoppler(varargin{:});
end
