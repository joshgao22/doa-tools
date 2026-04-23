function [lb, ub] = buildDoaDopplerMfBounds(model)
%BUILDDOADOPPLERMFBOUNDS Build MF optimizer box constraints.
% Keep DoA local-box construction and CP-U release-floor logic in one
% helper so the branch solver only orchestrates known/unknown paths.

arguments
  model (1,1) struct
end

if isfield(model, 'lb') && isfield(model, 'ub') && ...
    ~isempty(model.lb) && ~isempty(model.ub)
  lb = model.lb;
  ub = model.ub;
  return;
end

[doaLb, doaUb] = localBuildDoaBounds(model);
[doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb);

lb = [doaLb; model.fdRange(1)];
ub = [doaUb; model.fdRange(2)];

if strcmp(model.fdRateMode, 'unknown')
  lb = [lb; model.fdRateRange(1)];
  ub = [ub; model.fdRateRange(2)];
end
end


function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build the DoA search box used by MF refinements.

baseRange = model.doaGrid{1}.range;
doaLb = baseRange(:, 1);
doaUb = baseRange(:, 2);

if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  return;
end
if isempty(model.initDoaParam) || isempty(model.initDoaHalfWidth)
  return;
end

doaCenter = model.initDoaParam(:);
doaHalfWidth = model.initDoaHalfWidth(:);
doaLbLocal = doaCenter - doaHalfWidth;
doaUbLocal = doaCenter + doaHalfWidth;

if strcmp(model.doaType, 'angle')
  doaCenter(1) = mod(doaCenter(1), 2 * pi);
  doaLbLocal(1) = doaCenter(1) - doaHalfWidth(1);
  doaUbLocal(1) = doaCenter(1) + doaHalfWidth(1);
  if doaLbLocal(1) < baseRange(1, 1) || doaUbLocal(1) > baseRange(1, 2)
    doaLbLocal(1) = baseRange(1, 1);
    doaUbLocal(1) = baseRange(1, 2);
  end
end

doaLb = max(doaLb, doaLbLocal);
doaUb = min(doaUb, doaUbLocal);
invalidMask = doaLb > doaUb;
doaLb(invalidMask) = baseRange(invalidMask, 1);
doaUb(invalidMask) = baseRange(invalidMask, 2);
end


function [doaLb, doaUb] = localApplyUnknownDoaReleaseFloor(model, doaLb, doaUb, centerOverride)
%LOCALAPPLYUNKNOWNDOARELEASEFLOOR Enforce one minimum CP-U DoA release box.

if nargin < 4
  centerOverride = [];
end
if ~strcmp(model.fdRateMode, 'unknown')
  return;
end
if isfield(model, 'freezeDoa') && logical(model.freezeDoa)
  return;
end
if isfield(model, 'disableUnknownDoaReleaseFloor') && logical(model.disableUnknownDoaReleaseFloor)
  return;
end
if model.numSat <= 1 || ~strcmp(model.phaseMode, 'continuous')
  return;
end
if isempty(model.initDoaParam) && isempty(centerOverride)
  return;
end

baseRange = model.doaGrid{1}.range;
if isempty(baseRange) || any(size(baseRange) ~= [2, 2])
  return;
end

releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model);
if isempty(releaseHalfWidth)
  return;
end

currentWidth = doaUb - doaLb;
targetWidth = 2 * releaseHalfWidth;
if ~isempty(centerOverride)
  center = reshape(centerOverride, [], 1);
else
  center = reshape(model.initDoaParam, [], 1);
end

needRecenter = false;
if ~isempty(centerOverride)
  centerTol = 1e-12;
  needRecenter = any(center < doaLb + releaseHalfWidth - centerTol) || ...
    any(center > doaUb - releaseHalfWidth + centerTol);
end
if all(currentWidth >= targetWidth - 1e-12) && ~needRecenter
  return;
end

doaLbFloor = center - releaseHalfWidth;
doaUbFloor = center + releaseHalfWidth;

if strcmp(model.doaType, 'angle')
  center(1) = mod(center(1), 2 * pi);
  doaLbFloor(1) = center(1) - releaseHalfWidth(1);
  doaUbFloor(1) = center(1) + releaseHalfWidth(1);
  if doaLbFloor(1) < baseRange(1, 1) || doaUbFloor(1) > baseRange(1, 2)
    doaLbFloor(1) = baseRange(1, 1);
    doaUbFloor(1) = baseRange(1, 2);
  end
end

for iDim = 1:numel(doaLb)
  if needRecenter || currentWidth(iDim) < targetWidth(iDim)
    doaLb(iDim) = max(baseRange(iDim, 1), doaLbFloor(iDim));
    doaUb(iDim) = min(baseRange(iDim, 2), doaUbFloor(iDim));
  end
end

invalidMask = doaLb > doaUb;
doAResetLb = reshape(baseRange(invalidMask, 1), [], 1);
doAResetUb = reshape(baseRange(invalidMask, 2), [], 1);
doAIdx = find(invalidMask);
for iIdx = 1:numel(doAIdx)
  doaLb(doAIdx(iIdx)) = doAResetLb(iIdx);
  doaUb(doAIdx(iIdx)) = doAResetUb(iIdx);
end
end


function releaseHalfWidth = localResolveUnknownDoaReleaseHalfWidth(model)
%LOCALRESOLVEUNKNOWNDOARELEASEHALFWIDTH Resolve the CP-U DoA box floor.

if isfield(model, 'unknownDoaReleaseHalfWidth') && ~isempty(model.unknownDoaReleaseHalfWidth)
  releaseHalfWidth = reshape(model.unknownDoaReleaseHalfWidth, [], 1);
else
  if strcmp(model.doaType, 'angle')
    defaultHalfWidth = deg2rad(0.03);
  else
    defaultHalfWidth = 0.03;
  end
  releaseHalfWidth = repmat(defaultHalfWidth, 2, 1);
end

if numel(releaseHalfWidth) == 1
  releaseHalfWidth = repmat(releaseHalfWidth, 2, 1);
end
releaseHalfWidth = reshape(releaseHalfWidth, [], 1);
if numel(releaseHalfWidth) ~= 2 || any(~isfinite(releaseHalfWidth)) || any(releaseHalfWidth < 0)
  releaseHalfWidth = [];
end
end
