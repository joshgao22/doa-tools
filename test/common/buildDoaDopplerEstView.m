function view = buildDoaDopplerEstView(sceneRef, rxSigIn, gridSize, searchRange, E, viewOpt)
%BUILDDOADOPPLERESTVIEW Build a unified estimator view for dev scripts.
% The returned view uses one common layout for both single-frame and
% multi-frame DoA/DoA-Doppler test cases.
%
%Syntax:
%  view = buildDoaDopplerEstView(sceneRef, rxSigIn, gridSize, searchRange, E)
%  view = buildDoaDopplerEstView(..., viewOpt)
%
%Inputs:
%  sceneRef         - reference-frame scene used by the estimator
%  rxSigIn          - single-frame array snapshot or multi-frame cell input
%  gridSize         - DoA grid size
%  searchRange      - search range used by genDoaGrid
%  E                - reference ellipsoid
%  viewOpt          - optional options structure
%    .sceneSeq      - multi-frame scene sequence. Default: []
%
%Output:
%  view             - unified estimator view with fields
%    .sceneRef
%    .sceneSeq
%    .rxSigSf
%    .rxSigMf
%    .doaGrid
%    .frameMode
%    .timeOffsetSec

if nargin < 6 || isempty(viewOpt)
  viewOpt = struct();
end

sceneSeq = localGetFieldOrDefault(viewOpt, 'sceneSeq', []);
doaGridCell = localBuildDoaGridCell(sceneRef, gridSize, searchRange, E);
if numel(doaGridCell) == 1
  doaGridOut = doaGridCell{1};
else
  doaGridOut = doaGridCell;
end

view = struct();
view.sceneRef = sceneRef;
view.sceneSeq = sceneSeq;
view.doaGrid = doaGridOut;

if isempty(sceneSeq)
  view.frameMode = "single";
  view.rxSigSf = rxSigIn;
  view.rxSigMf = [];
  view.timeOffsetSec = 0;
else
  view.frameMode = "multi";
  view.rxSigMf = rxSigIn;
  view.rxSigSf = rxSigIn{sceneSeq.refFrameIdx};
  view.timeOffsetSec = reshape(sceneSeq.timeOffsetSec, [], 1);
end
end

function doaGridCell = localBuildDoaGridCell(sceneRef, gridSize, searchRange, E)
%LOCALBUILDDOAGRIDCELL Build one lat-lon search grid per satellite.

numSat = sceneRef.numSat;
doaGridCell = cell(1, numSat);
for iSat = 1:numSat
  doaGridCell{iSat} = genDoaGrid("latlon", 2, gridSize, searchRange, ...
    'eci', datevec(sceneRef.utc), sceneRef.satPosEci(:, iSat), ...
    sceneRef.rotMat{iSat}, E);
end
end

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end
