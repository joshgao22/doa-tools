function angleErrDeg = calcDoaAngleError(doaEst, doaTrue, doaType)
%CALCDOAANGLEERROR Compute equivalent global DoA angle error in degrees.
%
% Computes the angular separation between estimated and true DoA parameters
% after mapping them to a common global unit-direction representation.
% The function is intended for unified Monte Carlo post-processing, so both
% angle-mode and lat-lon-mode parameterizations share one output metric.
%
%Syntax:
%  angleErrDeg = calcDoaAngleError(doaEst, doaTrue)
%  angleErrDeg = calcDoaAngleError(doaEst, doaTrue, doaType)
%
%Inputs:
%  doaEst          - estimated DoA parameters, one estimate per column
%                    - doaType = 'angle'  : 2xN [az; el] in radians
%                    - doaType = 'latlon' : 2xN [lat; lon] in degrees
%
%  doaTrue         - true DoA parameters, one truth value per column
%                    Supports either one column or the same number of
%                    columns as doaEst. A single column will be expanded to
%                    match the other input.
%
%  doaType         - (optional) global DoA parameterization
%                    - 'angle'  : global [az; el] in radians
%                    - 'latlon' : [lat; lon] in degrees
%                    default: 'angle'
%
%Output:
%  angleErrDeg     - 1xN equivalent angular error in degrees
%
%Notes:
%  - For doaType = 'angle', the function assumes [az; el] follows the same
%    convention as doa2dir.
%  - For doaType = 'latlon', the metric is the great-circle angle between
%    the corresponding Earth-centered unit vectors. Therefore the result is
%    independent of longitude wrapping and does not require a UTC epoch.
%
%See also:
%  doa2dir

arguments
  doaEst {mustBeNumeric, mustBeReal}
  doaTrue {mustBeNumeric, mustBeReal}
  doaType = 'angle'
end

doaType = validatestring(doaType, {'angle', 'latlon'}, mfilename, 'doaType');
doaEst = localNormalizeDoaParam(doaEst, 'doaEst');
doaTrue = localNormalizeDoaParam(doaTrue, 'doaTrue');
[doaEst, doaTrue] = localMatchColumnCount(doaEst, doaTrue);

switch doaType
  case 'angle'
    dirEst = doa2dir(doaEst);
    dirTrue = doa2dir(doaTrue);

  case 'latlon'
    dirEst = localLatlonToDir(doaEst);
    dirTrue = localLatlonToDir(doaTrue);

  otherwise
    error('calcDoaAngleError:InvalidDoaType', ...
      'Unsupported doaType: %s.', doaType);
end

cosVal = sum(conj(dirEst) .* dirTrue, 1);
cosVal = real(cosVal);
cosVal = min(max(cosVal, -1), 1);
angleErrDeg = rad2deg(acos(cosVal));
end


function doaParam = localNormalizeDoaParam(doaParam, inputName)
%LOCALNORMALIZEDOAPARAM Convert input DoA parameters to a 2xN matrix.

if isempty(doaParam)
  error('calcDoaAngleError:EmptyInput', ...
    '%s must be non-empty.', inputName);
end

if size(doaParam, 1) == 2
  return;
end

if size(doaParam, 2) == 2
  doaParam = doaParam.';
  return;
end

error('calcDoaAngleError:InvalidInputSize', ...
  '%s must have size 2xN or Nx2.', inputName);
end


function [lhs, rhs] = localMatchColumnCount(lhs, rhs)
%LOCALMATCHCOLUMNCOUNT Expand singleton columns if needed.

numColLhs = size(lhs, 2);
numColRhs = size(rhs, 2);

if numColLhs == numColRhs
  return;
end

if numColLhs == 1
  lhs = repmat(lhs, 1, numColRhs);
  return;
end

if numColRhs == 1
  rhs = repmat(rhs, 1, numColLhs);
  return;
end

error('calcDoaAngleError:ColumnCountMismatch', ...
  'doaEst and doaTrue must have the same number of columns, or one input must have one column.');
end


function dirVec = localLatlonToDir(latlonDeg)
%LOCALLATLONTODIR Convert [lat; lon] in degrees to Earth-centered unit vectors.

latRad = deg2rad(latlonDeg(1, :));
lonRad = deg2rad(latlonDeg(2, :));

dirVec = [cos(latRad) .* cos(lonRad);
          cos(latRad) .* sin(lonRad);
          sin(latRad)];
end
