function [angleStdDeg, aux] = projectCrbToAngleMetric(crbDoa, doaParam, doaType)
%PROJECTCRBTOANGLEMETRIC Project a DoA CRB onto the equivalent angle metric.
%
% Converts a 2-parameter DoA covariance matrix into a scalar lower bound on
% the equivalent global angular RMSE used in the paper plots.
%
%Syntax:
%  angleStdDeg = projectCrbToAngleMetric(crbDoa, doaParam)
%  angleStdDeg = projectCrbToAngleMetric(crbDoa, doaParam, doaType)
%  [angleStdDeg, aux] = projectCrbToAngleMetric(...)
%
%Inputs:
%  crbDoa          - DoA covariance matrix or covariance pages.
%                    Supported layouts:
%                    - 2x2      : one DoA covariance matrix
%                    - KxK      : one full covariance matrix; the leading
%                                 2x2 block is interpreted as the DoA CRB
%                    - 2x2xN    : one DoA covariance page per parameter set
%                    - KxKxN    : one full covariance page per parameter
%                                 set; the leading 2x2 block is used
%
%  doaParam        - 2x1 or 2xN DoA parameter vector(s)
%                    angle mode : [az; el] in radians
%                    latlon mode: [lat; lon] in degrees
%
%  doaType         - (optional) DoA parameterization
%                    - 'angle'  : global [az; el] in radians
%                    - 'latlon' : [lat; lon] in degrees
%                    default: 'angle'
%
%Outputs:
%  angleStdDeg     - 1xN projected scalar angle lower bound in degrees
%
%  aux             - Auxiliary structure with fields:
%                    .doaType
%                    .dirVec
%                    .dirJacobian
%                    .metricMat
%                    .crbDoaUsed
%                    .angleVarRad2
%
%Description:
%  For a small DoA perturbation dp around the true direction, the squared
%  equivalent angle error satisfies
%
%      angle^2 \approx ||J_u dp||_2^2,
%
%  where J_u is the Jacobian of the unit-direction vector with respect to
%  the two DoA parameters. Therefore the scalar CRB in the angle metric is
%
%      E[angle^2] >= trace(J_u * C_p * J_u^T),
%
%  with C_p being the 2x2 DoA covariance block.
%
%Notes:
%  - The returned quantity is a first-order local approximation, which is
%    the natural way to project a parameter-space CRB onto the global angle
%    metric.
%  - For lat-lon inputs, the Jacobian is taken with respect to degrees, so
%    the covariance input should also be in degree-squared units.
%
%See also:
%  calcDoaAngleError, doa2dirJacobian

arguments
  crbDoa {mustBeNumeric, mustBeReal}
  doaParam {mustBeNumeric, mustBeReal}
  doaType = 'angle'
end

doaType = validatestring(doaType, {'angle', 'latlon'}, mfilename, 'doaType');
doaParam = localNormalizeDoaParam(doaParam, 'doaParam');
[crbPage, numPage] = localParseCrb(crbDoa);
doaParam = localMatchPageCount(doaParam, numPage);

angleVarRad2 = zeros(1, numPage);
dirVec = zeros(3, numPage);
dirJacobian = zeros(3, 2, numPage);
metricMat = zeros(2, 2, numPage);
crbUsed = zeros(2, 2, numPage);

for iPage = 1:numPage
  currentParam = doaParam(:, iPage);
  currentCrb = localSymmetrize(crbPage(:, :, iPage));

  switch doaType
    case 'angle'
      [currentDir, du] = doa2dirJacobian(currentParam);
      jacDir = [du.az, du.el];

    case 'latlon'
      [currentDir, jacDir] = localLatlonDirJacobian(currentParam);

    otherwise
      error('projectCrbToAngleMetric:InvalidDoaType', ...
        'Unsupported doaType: %s.', doaType);
  end

  currentMetric = jacDir.' * jacDir;
  currentVar = real(trace(jacDir * currentCrb * jacDir.'));
  angleVarRad2(iPage) = max(currentVar, 0);

  dirVec(:, iPage) = currentDir;
  dirJacobian(:, :, iPage) = jacDir;
  metricMat(:, :, iPage) = currentMetric;
  crbUsed(:, :, iPage) = currentCrb;
end

angleStdDeg = rad2deg(sqrt(angleVarRad2));

if nargout < 2
  return;
end

aux = struct();
aux.doaType = doaType;
aux.dirVec = dirVec;
aux.dirJacobian = dirJacobian;
aux.metricMat = metricMat;
aux.crbDoaUsed = crbUsed;
aux.angleVarRad2 = angleVarRad2;
end


function doaParam = localNormalizeDoaParam(doaParam, inputName)
%LOCALNORMALIZEDOAPARAM Convert input DoA parameters to a 2xN matrix.

if isempty(doaParam)
  error('projectCrbToAngleMetric:EmptyDoaParam', ...
    '%s must be non-empty.', inputName);
end

if size(doaParam, 1) == 2
  return;
end

if size(doaParam, 2) == 2
  doaParam = doaParam.';
  return;
end

error('projectCrbToAngleMetric:InvalidDoaParamSize', ...
  '%s must have size 2xN or Nx2.', inputName);
end


function [crbPage, numPage] = localParseCrb(crbIn)
%LOCALPARSECRB Normalize covariance input to a 2x2xN array.

if ndims(crbIn) > 3
  error('projectCrbToAngleMetric:InvalidCrbDim', ...
    'crbDoa must be a 2-D matrix or a 3-D covariance array.');
end

if size(crbIn, 1) < 2 || size(crbIn, 2) < 2 || size(crbIn, 1) ~= size(crbIn, 2)
  error('projectCrbToAngleMetric:InvalidCrbSize', ...
    'crbDoa must be a square matrix with size at least 2x2.');
end

numPage = size(crbIn, 3);
crbPage = crbIn(1:2, 1:2, :);
end


function doaParam = localMatchPageCount(doaParam, numPage)
%LOCALMATCHPAGECOUNT Expand or validate DoA parameter pages.

numParam = size(doaParam, 2);
if numParam == numPage
  return;
end

if numParam == 1
  doaParam = repmat(doaParam, 1, numPage);
  return;
end

if numPage == 1
  return;
end

error('projectCrbToAngleMetric:PageCountMismatch', ...
  'The number of doaParam columns must match the covariance pages, or doaParam must contain one column.');
end


function covMat = localSymmetrize(covMat)
%LOCALSYMMETRIZE Force Hermitian symmetry in the covariance block.

covMat = real(0.5 * (covMat + covMat.'));
end


function [dirVec, jacDir] = localLatlonDirJacobian(latlonDeg)
%LOCALLATLONDIRJACOBIAN Unit vector and Jacobian for [lat; lon] in degrees.

latRad = deg2rad(latlonDeg(1));
lonRad = deg2rad(latlonDeg(2));

cosLat = cos(latRad);
sinLat = sin(latRad);
cosLon = cos(lonRad);
sinLon = sin(lonRad);

dirVec = [cosLat * cosLon;
          cosLat * sinLon;
          sinLat];

scale = pi / 180;
dLat = scale * [-sinLat * cosLon;
                -sinLat * sinLon;
                 cosLat];
dLon = scale * [-cosLat * sinLon;
                 cosLat * cosLon;
                 0];

jacDir = [dLat, dLon];
end
