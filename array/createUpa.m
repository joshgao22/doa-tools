function array = createUpa(arraySize, spacing, positionErrors, gainErrors, phaseErrors)
%CREATEUPA Create a uniform planar array (UPA) as a struct.
% Generates a 2D uniform rectangular array on the xy-plane.
%
%Syntax:
%   array = createUpa(arraySize, spacing)
%   array = createUpa(arraySize, spacing, positionErrors, gainErrors, phaseErrors)
%
%Inputs:
%   arraySize      - [numX, numY], number of elements along x/y axes
%   spacing        - Inter-element spacing (meters)
%   positionErrors - (Optional) 2×M numeric position perturbations (meters)
%   gainErrors     - (Optional) 2×M numeric gain perturbations (linear)
%   phaseErrors    - (Optional) 2×M numeric phase perturbations (radians)
%
%Output:
%   array - Struct with fields (BaseArray compatible):
%     positions       : 2×M element positions (meters)
%     dim             : 2
%     count           : M
%     indices         : 2×M (optional)
%     positionErrors  : 2×M (optional)
%     gainErrors      : 2×M (optional)
%     phaseErrors     : 2×M (optional)
%     spacing         : scalar
%
%Notes:
%   - Elements lie on the xy-plane (z=0).
%   - Indexing order: x varies fastest, then y (MATLAB column-major style).
%
%See also: createUla, steeringMatrix

arguments
  arraySize (1,2) {mustBePositive, mustBeInteger}
  spacing (1,1) {mustBePositive, mustBeNumeric}
  positionErrors {mustBeNumeric} = []
  gainErrors {mustBeNumeric} = []
  phaseErrors {mustBeNumeric} = []
end

numX = double(arraySize(1));
numY = double(arraySize(2));
numElements = numX * numY;

% -------------------------------------------------------------------------
% Generate grid indices (x-fastest ordering)
% -------------------------------------------------------------------------
[xIdx, yIdx] = meshgrid(0:numX-1, 0:numY-1);
xIdx = xIdx(:).';
yIdx = yIdx(:).';
indices = [xIdx; yIdx];

% -------------------------------------------------------------------------
% Element positions: 2×M (xy-plane)
% -------------------------------------------------------------------------
positions = [xIdx * spacing;
             yIdx * spacing];

dim   = 2;
count = numElements;

% -------------------------------------------------------------------------
% Errors (optional): must match positions size
% -------------------------------------------------------------------------
if ~isempty(positionErrors)
  if ~isequal(size(positionErrors), size(positions))
    error(['Size mismatch: positionErrors ', mat2str(size(positionErrors)), ...
      ' vs positions ', mat2str(size(positions))]);
  end
end

if ~isempty(gainErrors)
  if ~isequal(size(gainErrors), size(positions))
    error(['Size mismatch: gainErrors ', mat2str(size(gainErrors)), ...
      ' vs positions ', mat2str(size(positions))]);
  end
end

if ~isempty(phaseErrors)
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
