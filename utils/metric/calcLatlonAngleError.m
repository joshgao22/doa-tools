function [angleErrDeg, latErrDeg, lonErrDeg] = calcLatlonAngleError(latlonEst, latlonTrue)
%CALCLATLONANGLEERROR Compute angular and signed lat-lon errors in degrees.
%
% Computes the equivalent global angular error between estimated and true
% [lat; lon] parameters, while also returning the signed latitude and
% wrapped longitude errors for diagnostic use.
%
%Syntax:
%  angleErrDeg = calcLatlonAngleError(latlonEst, latlonTrue)
%
%  [angleErrDeg, latErrDeg, lonErrDeg] = calcLatlonAngleError(latlonEst, latlonTrue)
%
%Inputs:
%  latlonEst       - estimated [lat; lon] parameters in degrees, one
%                    estimate per column
%
%  latlonTrue      - true [lat; lon] parameters in degrees
%                    Supports either one column or the same number of
%                    columns as latlonEst. A single column will be
%                    expanded to match the other input.
%
%Outputs:
%  angleErrDeg     - 1xN equivalent global angular error in degrees
%
%  latErrDeg       - 1xN signed latitude error in degrees
%
%  lonErrDeg       - 1xN signed longitude error in degrees after wrapping
%                    to the interval [-180, 180)
%
%Description:
%  This function is a lat-lon convenience wrapper for unified metric
%  evaluation in Monte Carlo scripts. The main output angleErrDeg is the
%  great-circle separation angle, while latErrDeg and lonErrDeg are useful
%  for debugging parameter-wise estimator bias.
%
%Notes:
%  - The angular error is rotation-invariant, so no UTC epoch is required.
%  - Longitude differences are wrapped to avoid artificial jumps across the
%    +/-180 degree discontinuity.
%
%See also:
%  calcDoaAngleError

arguments
  latlonEst {mustBeNumeric, mustBeReal}
  latlonTrue {mustBeNumeric, mustBeReal}
end

latlonEst = localNormalizeLatlon(latlonEst, 'latlonEst');
latlonTrue = localNormalizeLatlon(latlonTrue, 'latlonTrue');
[latlonEst, latlonTrue] = localMatchColumnCount(latlonEst, latlonTrue);

latErrDeg = latlonEst(1, :) - latlonTrue(1, :);
lonErrDeg = localWrapTo180(latlonEst(2, :) - latlonTrue(2, :));
angleErrDeg = calcDoaAngleError(latlonEst, latlonTrue, 'latlon');
end


function latlon = localNormalizeLatlon(latlon, inputName)
%LOCALNORMALIZELATLON Convert input latitude-longitude parameters to 2xN.

if isempty(latlon)
  error('calcLatlonAngleError:EmptyInput', ...
    '%s must be non-empty.', inputName);
end

if size(latlon, 1) == 2
  return;
end

if size(latlon, 2) == 2
  latlon = latlon.';
  return;
end

error('calcLatlonAngleError:InvalidInputSize', ...
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

error('calcLatlonAngleError:ColumnCountMismatch', ...
  'latlonEst and latlonTrue must have the same number of columns, or one input must have one column.');
end


function lonDeg = localWrapTo180(lonDeg)
%LOCALWRAPTO180 Wrap angles in degrees to the interval [-180, 180).

lonDeg = mod(lonDeg + 180, 360) - 180;
end
