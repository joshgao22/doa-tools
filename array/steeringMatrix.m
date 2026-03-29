function [A, dA] = steeringMatrix(array, wavelength, doas, useNormalizedAngles)
%STEERINGMATRIX Create the steering matrix for arrays.
%
%Syntax:
%   A = steeringMatrix(array, wavelength, doas)
%   [A, dA] = steeringMatrix(array, wavelength, doas)
%   ... = steeringMatrix(array, wavelength, doas, useNormalizedAngles)
%
%Inputs:
%   array - Struct with fields:
%           positions       : pDim×M numeric
%           type            : arrayDesign.ArrayType
%           spacing         : scalar (required if useNormalizedAngles=true)
%           indices         : (optional)
%           positionErrors  : 1×M or 3×M
%           gainErrors      : M×1
%           phaseErrors     : M×1
%
%   wavelength - Positive scalar
%   doas       - 1×N (broadside) or 2×N ([az; el])
%   useNormalizedAngles - logical (default: false)
%
%Outputs:
%   A  - Steering matrix (M×N)
%   dA - Struct of derivatives
%        1D / broadside-only: dA.theta  (M×N)
%        az/el form         : dA.az, dA.el (each M×N)

arguments
  array (1,1) struct
  wavelength (1,1) {mustBePositive, mustBeNumeric}
  doas {mustBeNumeric}
  useNormalizedAngles (1,1) logical = false
end

needDerivative = (nargout == 2);

if needDerivative
  dA = struct();
else
  dA = [];
end

elementPositions = array.positions;
arrayDim = size(elementPositions, 1);
numElements = size(elementPositions, 2);

% -------------------------------------------------------------------------
% Apply position errors
% -------------------------------------------------------------------------
if isfield(array, 'positionErrors') && ~isempty(array.positionErrors)
  elementPositions = applyPositionErrors(elementPositions, array.positionErrors);
  arrayDim = size(elementPositions, 1);
end

% -------------------------------------------------------------------------
% Steering matrix computation
% -------------------------------------------------------------------------
if arrayDim == 1
  if ~isvector(doas)
    broadside = ae2broad(doas(1,:), doas(2,:));
  else
    broadside = reshape(doas, 1, []);
  end

  if useNormalizedAngles
    posNorm = elementPositions(:) / array.spacing;
    A = exp(1j * (posNorm * broadside));
    if needDerivative
      dA.theta = (1j * posNorm) .* A;
    end
  else
    k = 2*pi / wavelength;
    x = elementPositions(:);
    A = exp(1j * k * (x * sin(broadside)));
    if needDerivative
      dA.theta = A .* (1j * k * (x * cos(broadside)));
    end
  end

else
  if arrayDim ~= 2 && arrayDim ~= 3
    error('steeringMatrix:InvalidArrayDim', ...
      'array.positions must be 1D, 2D, or 3D.');
  end
  if useNormalizedAngles
    error('steeringMatrix:NormalizedAnglesNotSupported', ...
      'Normalized angles are not supported for 2D/3D arrays.');
  end

  k = 2*pi / wavelength;

  x = elementPositions(1,:).';
  y = elementPositions(2,:).';

  if size(doas,1) == 1
    % broadside only (assume source on xy-plane, angle relative to y-axis)
    broadside = reshape(doas, 1, []);
    A = exp(1j * k * (x * sin(broadside) + y * cos(broadside)));
    if needDerivative
      dA.theta = A .* (1j * k * (x * cos(broadside) - y * sin(broadside)));
    end

  else
    % azimuth / elevation
    az = doas(1,:);
    el = doas(2,:);

    cosAz = cos(az);
    sinAz = sin(az);
    cosEl = cos(el);
    sinEl = sin(el);

    % direction cosines
    cc = cosEl .* cosAz;   % x coefficient
    cs = cosEl .* sinAz;   % y coefficient

    if arrayDim == 2
      A = exp(1j * k * (x * cc + y * cs));

      if needDerivative
        % d/d az
        dphi_daz = k * (-x * (cosEl .* sinAz) + y * (cosEl .* cosAz));
        dA.az = 1j * dphi_daz .* A;

        % d/d el   (z = 0)
        dphi_del = k * (-x * (sinEl .* cosAz) - y * (sinEl .* sinAz));
        dA.el = 1j * dphi_del .* A;
      end

    else
      z = elementPositions(3,:).';
      A = exp(1j * k * (x * cc + y * cs + z * sinEl));

      if needDerivative
        % d/d az
        dphi_daz = k * (-x * (cosEl .* sinAz) + y * (cosEl .* cosAz));
        dA.az = 1j * dphi_daz .* A;

        % d/d el
        dphi_del = k * (-x * (sinEl .* cosAz) ...
                      - y * (sinEl .* sinAz) ...
                      + z * cosEl);
        dA.el = 1j * dphi_del .* A;
      end
    end
  end
end

% -------------------------------------------------------------------------
% Apply gain / phase errors
% -------------------------------------------------------------------------
if isfield(array, 'gainErrors') && ~isempty(array.gainErrors)
  g = array.gainErrors(:);
  if numel(g) ~= numElements
    error('steeringMatrix:GainErrorSizeMismatch', ...
      'array.gainErrors length must match number of elements.');
  end
  A = g .* A;

  if needDerivative
    fn = fieldnames(dA);
    for ii = 1:numel(fn)
      dA.(fn{ii}) = g .* dA.(fn{ii});
    end
  end
end

if isfield(array, 'phaseErrors') && ~isempty(array.phaseErrors)
  p = array.phaseErrors(:);
  if numel(p) ~= numElements
    error('steeringMatrix:PhaseErrorSizeMismatch', ...
      'array.phaseErrors length must match number of elements.');
  end
  ph = exp(1j * p);
  A = ph .* A;

  if needDerivative
    fn = fieldnames(dA);
    for ii = 1:numel(fn)
      dA.(fn{ii}) = ph .* dA.(fn{ii});
    end
  end
end

end