function tf = shouldUseDoaDopplerMfWarmAnchorParfor(model, numReleaseSeed, verbose)
%SHOULDUSEDOADOPPLERMFWARMANCHORPARFOR Gate warm-anchor release-seed parfor.
% The warm-anchor release seeds are independent solves, but this gate keeps
% the default estimator path serial.  In particular, verbose traces and
% fixed-DoA tooth-guard branches must never switch to parfor because small
% worker-side optimizer differences can change the final tooth winner.

if nargin < 3 || isempty(verbose)
  verbose = false;
end

% Keep the conservative default unless every opt-in guard passes.
tf = false;
if verbose
  return;
end
if localGetLogicalField(model, 'freezeDoa', false)
  return;
end
if ~localGetLogicalField(model, 'unknownWarmAnchorUseParfor', false)
  return;
end
minSeed = localGetFiniteField(model, 'unknownWarmAnchorMinParforSeed', 3);
if ~(isfinite(minSeed) && numReleaseSeed >= minSeed)
  return;
end

tf = localCanUseParfor();
end

function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether parfor can be used safely here.

tf = false;
if exist('getCurrentTask', 'file') == 2 && ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  poolObj = gcp('nocreate');
  tf = ~isempty(poolObj);
catch
  tf = false;
end
end

function value = localGetLogicalField(dataStruct, fieldName, defaultValue)
%LOCALGETLOGICALFIELD Read one logical-like scalar field.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
value = logical(rawValue(1));
end

function value = localGetFiniteField(dataStruct, fieldName, defaultValue)
%LOCALGETFINITEFIELD Read one finite scalar field.

value = defaultValue;
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
rawValue = dataStruct.(fieldName);
if isempty(rawValue)
  return;
end
rawValue = rawValue(1);
if isfinite(rawValue)
  value = rawValue;
end
end
