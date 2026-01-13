function array = createUla(numElements, spacing, positionErrors, gainErrors, phaseErrors)
%CREATEULA Create a uniform linear array (ULA) as a struct.
%
%Syntax:
%   array = createUla(numElements, spacing)
%   array = createUla(numElements, spacing, positionErrors, gainErrors, phaseErrors)
%
%Inputs:
%   numElements    - Number of array elements (positive integer)
%   spacing        - Inter-element spacing (positive scalar)
%   positionErrors - (Optional) 1×M numeric perturbations (meters) for ULA
%   gainErrors     - (Optional) 1×M numeric gain perturbations (linear)
%   phaseErrors    - (Optional) 1×M numeric phase perturbations (radians)
%
%Output:
%   array - Struct with fields (BaseArray compatible):
%     positions       : 1×M (meters)
%     dim             : 1
%     count           : M
%     indices         : 1×M (optional)
%     positionErrors  : 1×M (optional)
%     gainErrors      : 1×M (optional)
%     phaseErrors     : 1×M (optional)
%     spacing         : scalar
%
%See also: steeringMatrix

arguments
  numElements (1,1) {mustBePositive, mustBeInteger}
  spacing (1,1) {mustBePositive, mustBeNumeric}
  positionErrors {mustBeNumeric} = []
  gainErrors {mustBeNumeric} = []
  phaseErrors {mustBeNumeric} = []
end

% -------------------------------------------------------------------------
% Generate ULA positions (1×M)
% -------------------------------------------------------------------------
m = double(numElements);
idx = (0:m-1);

positions = double(idx) * double(spacing);   % 1×M
positions = reshape(positions, 1, []);       % enforce row

dim = 1;
count = m;
indices = 0:(numElements - 1);

% -------------------------------------------------------------------------
% Errors (optional): same size as positions
% -------------------------------------------------------------------------
if ~isempty(positionErrors)
  positionErrors = double(positionErrors);
  positionErrors = reshape(positionErrors, 1, []);
  if ~isequal(size(positionErrors), size(positions))
    error(['Size mismatch: positionErrors ', mat2str(size(positionErrors)), ...
      ' vs positions ', mat2str(size(positions))]);
  end
end

if ~isempty(gainErrors)
  gainErrors = double(gainErrors);
  gainErrors = reshape(gainErrors, 1, []);
  if ~isequal(size(gainErrors), size(positions))
    error(['Size mismatch: gainErrors ', mat2str(size(gainErrors)), ...
      ' vs positions ', mat2str(size(positions))]);
  end
end

if ~isempty(phaseErrors)
  phaseErrors = double(phaseErrors);
  phaseErrors = reshape(phaseErrors, 1, []);
  if ~isequal(size(phaseErrors), size(positions))
    error(['Size mismatch: phaseErrors ', mat2str(size(phaseErrors)), ...
      ' vs positions ', mat2str(size(positions))]);
  end
end

% -------------------------------------------------------------------------
% Assemble output struct
% -------------------------------------------------------------------------
array = struct();
array.positions      = positions;
array.dim            = dim;
array.count          = count;

array.indices        = indices;
array.positionErrors = positionErrors;
array.gainErrors     = gainErrors;
array.phaseErrors    = phaseErrors;

array.spacing        = double(spacing);

end
