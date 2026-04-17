function [model, lb, ub, modelOpt] = buildDoaDopplerSfModel(scene, rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, fdRange, numSource, modelOpt)
%BUILDDOADOPPLERSFMODEL Build the single-frame DoA-Doppler estimator model.
% This helper extracts the model-construction stage from
% estimatorDoaDopplerMlePilotSfOpt and returns the parsed model together
% with the optimization bounds.

arguments
  scene (1,1) struct
  rxSig
  pilotWave
  carrierFreq (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  sampleRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  doaGrid
  fdRange = []
  numSource (1,1) {mustBePositive, mustBeInteger} = 1
  modelOpt (1,1) struct = struct()
end

modelOpt = localParseModelOpt(modelOpt);
model = localBuildModel(scene, rxSig, pilotWave, carrierFreq, sampleRate, ...
  doaGrid, fdRange, numSource, modelOpt);
[lb, ub] = localBuildBounds(model);
end

function modelOpt = localParseModelOpt(modelOpt)
%LOCALPARSEMODELOPT Parse optional estimator settings.

if ~isfield(modelOpt, 'lightSpeed') || isempty(modelOpt.lightSpeed)
  modelOpt.lightSpeed = 299792458;
end
if ~isfield(modelOpt, 'useLogObjective') || isempty(modelOpt.useLogObjective)
  modelOpt.useLogObjective = true;
end
if ~isfield(modelOpt, 'initFdCount') || isempty(modelOpt.initFdCount)
  modelOpt.initFdCount = 65;
end
if ~isfield(modelOpt, 'optimOpt') || isempty(modelOpt.optimOpt)
  modelOpt.optimOpt = struct();
end
if ~isfield(modelOpt, 'initDoaParam') || isempty(modelOpt.initDoaParam)
  modelOpt.initDoaParam = [];
end
if ~isfield(modelOpt, 'initDoaHalfWidth') || isempty(modelOpt.initDoaHalfWidth)
  modelOpt.initDoaHalfWidth = [];
end
if ~isfield(modelOpt, 'satWeight') || isempty(modelOpt.satWeight)
  modelOpt.satWeight = [];
end
if ~isfield(modelOpt, 'evalOnly') || isempty(modelOpt.evalOnly)
  modelOpt.evalOnly = false;
end
if ~isfield(modelOpt, 'useScaledSolver') || isempty(modelOpt.useScaledSolver)
  modelOpt.useScaledSolver = true;
end
if ~isfield(modelOpt, 'enableDoaAnchorFallback') || isempty(modelOpt.enableDoaAnchorFallback)
  modelOpt.enableDoaAnchorFallback = true;
end
if ~isfield(modelOpt, 'doaAnchorFallbackObjTol') || isempty(modelOpt.doaAnchorFallbackObjTol)
  modelOpt.doaAnchorFallbackObjTol = 2.0;
end

if ~isscalar(modelOpt.lightSpeed) || ~isfinite(modelOpt.lightSpeed) || modelOpt.lightSpeed <= 0
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidLightSpeed', ...
    'modelOpt.lightSpeed must be a positive finite scalar.');
end
if ~isscalar(modelOpt.initFdCount) || modelOpt.initFdCount <= 0 || mod(modelOpt.initFdCount, 1) ~= 0
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidInitFdCount', ...
    'modelOpt.initFdCount must be a positive integer scalar.');
end
if ~isstruct(modelOpt.optimOpt)
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidOptimOpt', ...
    'modelOpt.optimOpt must be a struct.');
end
if ~isempty(modelOpt.initDoaParam)
  if ~isnumeric(modelOpt.initDoaParam) || size(modelOpt.initDoaParam, 1) ~= 2 || ...
      any(~isfinite(modelOpt.initDoaParam(:)))
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidInitDoaParam', ...
      'modelOpt.initDoaParam must be empty or a finite 2xK array.');
  end
end
if ~isempty(modelOpt.initDoaHalfWidth)
  if ~isnumeric(modelOpt.initDoaHalfWidth) || numel(modelOpt.initDoaHalfWidth) ~= 2 || ...
      any(~isfinite(modelOpt.initDoaHalfWidth(:))) || any(modelOpt.initDoaHalfWidth(:) < 0)
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidInitDoaHalfWidth', ...
      'modelOpt.initDoaHalfWidth must be empty or a nonnegative 2x1 vector.');
  end
  if isempty(modelOpt.initDoaParam)
    error('estimatorDoaDopplerMlePilotSfOpt:MissingInitDoaParam', ...
      'modelOpt.initDoaHalfWidth requires modelOpt.initDoaParam.');
  end
end
if ~isempty(modelOpt.satWeight)
  if ~isnumeric(modelOpt.satWeight) || ~isvector(modelOpt.satWeight) || ...
      any(~isfinite(modelOpt.satWeight(:))) || any(modelOpt.satWeight(:) < 0)
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidSatWeight', ...
      'modelOpt.satWeight must be empty or a nonnegative finite vector.');
  end
end
if ~isscalar(modelOpt.evalOnly) || (~islogical(modelOpt.evalOnly) && ~ismember(modelOpt.evalOnly, [0, 1]))
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidEvalOnly', ...
    'modelOpt.evalOnly must be a logical scalar.');
end
if ~isscalar(modelOpt.useScaledSolver) || ...
    (~islogical(modelOpt.useScaledSolver) && ~ismember(modelOpt.useScaledSolver, [0, 1]))
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidScaledSolver', ...
    'modelOpt.useScaledSolver must be a logical scalar.');
end
modelOpt.evalOnly = logical(modelOpt.evalOnly);
end

function model = localBuildModel(scene, rxSig, pilotWave, carrierFreq, sampleRate, ...
  doaGridIn, fdRange, numSource, modelOpt)
%LOCALBUILDMODEL Build fixed quantities used by the estimator.

[rxCell, numSat, numSample] = localParseRxSig(rxSig);
arrayCell = localParseSceneArray(scene);
rotMatCell = localParseRotMat(scene);
[satPosEci, satVelEci] = localParseSatState(scene);
sceneUtc = localParseSceneUtc(scene, doaGridIn);

doagrid = localParseDoaGrid(doaGridIn, numSat, sceneUtc);

[arrayCell, rotMatCell] = localBroadcastSceneContainers(arrayCell, rotMatCell, numSat);
[satPosEci, satVelEci] = localBroadcastSatState(satPosEci, satVelEci, numSat);
[refState, refSatIdxLocal] = resolveReferenceSatState(scene, satPosEci, satVelEci);

if numel(arrayCell) ~= numSat
  error('estimatorDoaDopplerMlePilotSfOpt:ArrayCountMismatch', ...
    'Number of scene arrays must match the number of rxSig entries.');
end
if numel(rotMatCell) ~= numSat
  error('estimatorDoaDopplerMlePilotSfOpt:RotMatCountMismatch', ...
    'Number of scene rotation matrices must match the number of rxSig entries.');
end
if size(satPosEci, 2) ~= numSat || size(satVelEci, 2) ~= numSat
  error('estimatorDoaDopplerMlePilotSfOpt:SatCountMismatch', ...
    'scene.satPosEci and scene.satVelEci must match the number of rxSig entries.');
end
if numel(doagrid) ~= numSat
  error('estimatorDoaDopplerMlePilotSfOpt:DoaGridCountMismatch', ...
    'Number of doaGrid entries must match the number of rxSig entries.');
end

[pilotPad, pilotEnergy, pilotLen] = localParsePilotWave(pilotWave, numSource, numSample);
[numElement, sampleCovCell] = localBuildSampleCov(rxCell, arrayCell, numSample);
[doaType, globalFrame, eciAngleGrid, eciDirectionGrid, latlonGrid] = localParseDoaMetadata(doagrid);
fdRange = localParseFdRange(fdRange, carrierFreq, modelOpt.lightSpeed, satVelEci, refState.velEci);

model = struct();
model.array = arrayCell;
model.rotMat = rotMatCell;
model.rxSig = rxCell;
model.sampleCov = sampleCovCell;
model.pilotPad = pilotPad;
model.pilotLen = pilotLen;
model.pilotEnergy = pilotEnergy;
model.numSat = numSat;
model.numSource = numSource;
model.numSample = numSample;
model.numElement = numElement;
model.sampleRate = sampleRate;
model.carrierFreq = carrierFreq;
model.lightSpeed = modelOpt.lightSpeed;
model.wavelength = modelOpt.lightSpeed / carrierFreq;
model.timeSec = (0:numSample-1) / sampleRate;
model.useLogObjective = logical(modelOpt.useLogObjective);
model.doaGrid = doagrid;
model.doaType = doaType;
model.globalFrame = globalFrame;
model.eciAngleGrid = eciAngleGrid;
model.eciDirectionGrid = eciDirectionGrid;
model.latlonGrid = latlonGrid;
model.fdRange = fdRange;
model.sceneUtc = sceneUtc;
model.userStateRef = localBuildUserStateRef(scene, sceneUtc);
model.satPosEci = satPosEci;
model.satVelEci = satVelEci;
model.sceneRef = scene;
model.ref = refState;
model.refSatIdxLocal = refSatIdxLocal;
model.refSatIdxGlobal = localGetStructField(refState, 'satIdxGlobal', NaN);
model.refStateSource = char(localGetStructField(refState, 'source', ""));
model.relSatVelEci = satVelEci - refState.velEci;
model.initFdCount = modelOpt.initFdCount;
model.initDoaParam = modelOpt.initDoaParam;
model.initDoaHalfWidth = reshape(modelOpt.initDoaHalfWidth, [], 1);
model.satWeight = localParseSatWeight(modelOpt.satWeight, numSat);
model.evalOnly = logical(modelOpt.evalOnly);
model.useScaledSolver = logical(modelOpt.useScaledSolver);
model.enableDoaAnchorFallback = logical(modelOpt.enableDoaAnchorFallback);
model.doaAnchorFallbackObjTol = double(modelOpt.doaAnchorFallbackObjTol);
model.optimOpt = modelOpt.optimOpt;
[doaLbCache, doaUbCache] = localBuildDoaBounds(model);
model.doaLb = doaLbCache;
model.doaUb = doaUbCache;
model.lb = [doaLbCache(:); repmat(model.fdRange(1), model.numSource, 1)];
model.ub = [doaUbCache(:); repmat(model.fdRange(2), model.numSource, 1)];
model.cachedBoundsReady = true;
end

function satWeight = localParseSatWeight(satWeightIn, numSat)
%LOCALPARSESATWEIGHT Parse one optional per-satellite objective weight vector.

if isempty(satWeightIn)
  satWeight = ones(numSat, 1);
  return;
end

satWeight = reshape(satWeightIn, [], 1);
if numel(satWeight) ~= numSat
  error('estimatorDoaDopplerMlePilotSfOpt:SatWeightCountMismatch', ...
    'modelOpt.satWeight must contain one entry per satellite.');
end
if all(satWeight <= 0)
  error('estimatorDoaDopplerMlePilotSfOpt:ZeroSatWeight', ...
    'At least one entry in modelOpt.satWeight must be positive.');
end
end

function [numElement, sampleCovCell] = localBuildSampleCov(rxCell, arrayCell, numSample)
%LOCALBUILDSAMPLECOV Build per-satellite sample covariance matrices.

numSat = numel(rxCell);
numElement = zeros(numSat, 1);
sampleCovCell = cell(1, numSat);

for iSat = 1:numSat
  currentRx = rxCell{iSat};
  if size(currentRx, 2) ~= numSample
    error('estimatorDoaDopplerMlePilotSfOpt:RxLengthMismatch', ...
      'All rxSig entries must share the same snapshot length.');
  end

  currentArray = arrayCell{iSat};
  if ~isstruct(currentArray) || ~isfield(currentArray, 'positions')
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidArray', ...
      'Each scene.array entry must be a valid array struct.');
  end

  numElement(iSat) = size(currentArray.positions, 2);
  if size(currentRx, 1) ~= numElement(iSat)
    error('estimatorDoaDopplerMlePilotSfOpt:RxElementMismatch', ...
      'rxSig{%d} row size must match the number of elements in scene.array{%d}.', iSat, iSat);
  end

  sampleCovCell{iSat} = (currentRx * currentRx') / numSample;
end
end

function [lb, ub] = localBuildBounds(model)
%LOCALBUILDBOUNDS Build box constraints for fmincon.

if isfield(model, 'lb') && isfield(model, 'ub') && ...
    ~isempty(model.lb) && ~isempty(model.ub)
  lb = model.lb;
  ub = model.ub;
  return;
end

[doaLb, doaUb] = localBuildDoaBounds(model);
lb = [doaLb(:); repmat(model.fdRange(1), model.numSource, 1)];
ub = [doaUb(:); repmat(model.fdRange(2), model.numSource, 1)];
if isfield(model, 'cachedBoundsReady') && ~logical(model.cachedBoundsReady)
  model.lb = lb;
  model.ub = ub;
end
end

function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build a per-source DoA search box.

if isfield(model, 'cachedBoundsReady') && logical(model.cachedBoundsReady) && ...
    isfield(model, 'doaLb') && isfield(model, 'doaUb') && ...
    ~isempty(model.doaLb) && ~isempty(model.doaUb)
  doaLb = model.doaLb;
  doaUb = model.doaUb;
  return;
end

baseRange = model.doaGrid{1}.range;
doaLb = repmat(baseRange(:, 1), 1, model.numSource);
doaUb = repmat(baseRange(:, 2), 1, model.numSource);

if isempty(model.initDoaParam) || isempty(model.initDoaHalfWidth)
  return;
end

doaCenter = reshape(model.initDoaParam, 2, []);
if size(doaCenter, 2) == 1 && model.numSource > 1
  doaCenter = repmat(doaCenter, 1, model.numSource);
end
if size(doaCenter, 2) ~= model.numSource
  error('estimatorDoaDopplerMlePilotSfOpt:InitDoaSourceMismatch', ...
    'modelOpt.initDoaParam must provide one 2x1 initializer per source.');
end

doaHalfWidth = reshape(model.initDoaHalfWidth, 2, 1);
for iSrc = 1:model.numSource
  currentCenter = doaCenter(:, iSrc);
  currentLb = currentCenter - doaHalfWidth;
  currentUb = currentCenter + doaHalfWidth;

  if strcmp(model.doaType, 'angle')
    currentCenter(1) = mod(currentCenter(1), 2 * pi);
    currentLb(1) = currentCenter(1) - doaHalfWidth(1);
    currentUb(1) = currentCenter(1) + doaHalfWidth(1);
    if currentLb(1) < baseRange(1, 1) || currentUb(1) > baseRange(1, 2)
      currentLb(1) = baseRange(1, 1);
      currentUb(1) = baseRange(1, 2);
    end
  end

  doaLb(:, iSrc) = max(doaLb(:, iSrc), currentLb);
  doaUb(:, iSrc) = min(doaUb(:, iSrc), currentUb);
end

invalidMask = doaLb > doaUb;
for iSrc = 1:model.numSource
  badIdx = invalidMask(:, iSrc);
  doaLb(badIdx, iSrc) = baseRange(badIdx, 1);
  doaUb(badIdx, iSrc) = baseRange(badIdx, 2);
end
end

function userStateRef = localBuildUserStateRef(scene, sceneUtc)
%LOCALBUILDUSERSTATEREF Build the user-state reference used by lat/lon mode.

userStateRef = struct();
userStateRef.utc = sceneUtc;

copyField = {'usrPosEci', 'usrVelEci', 'usrPosNominalEci', 'usrVelNominalEci'};
for iField = 1:numel(copyField)
  fieldName = copyField{iField};
  if isfield(scene, fieldName) && ~isempty(scene.(fieldName))
    userStateRef.(fieldName) = scene.(fieldName);
  end
end
end

function [rxCell, numSat, numSample] = localParseRxSig(rxSig)
%LOCALPARSERXSIG Normalize rxSig to a 1xNs cell array.

if iscell(rxSig)
  rxCell = reshape(rxSig, 1, []);
elseif isnumeric(rxSig)
  rxCell = {rxSig};
else
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidRxSigType', ...
    'rxSig must be a numeric matrix or a cell array of matrices.');
end

numSat = numel(rxCell);
if numSat == 0
  error('estimatorDoaDopplerMlePilotSfOpt:EmptyRxSig', ...
    'rxSig must not be empty.');
end

numSample = size(rxCell{1}, 2);
for iSat = 1:numSat
  currentRx = rxCell{iSat};
  if ~ismatrix(currentRx) || isempty(currentRx) || ~isnumeric(currentRx)
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidRxSigEntry', ...
      'Each rxSig entry must be a non-empty numeric matrix.');
  end
  if any(~isfinite(real(currentRx(:)))) || any(~isfinite(imag(currentRx(:))))
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidRxSigValue', ...
      'Each rxSig entry must contain finite values.');
  end
end
end

function [pilotPad, pilotEnergy, pilotLen] = localParsePilotWave(pilotWave, numSource, numSample)
%LOCALPARSEPILOTWAVE Convert pilotWave to a numSource x N template matrix.

if isempty(pilotWave) || ~isnumeric(pilotWave)
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidPilotWave', ...
    'pilotWave must be a non-empty numeric vector or matrix.');
end

if isvector(pilotWave)
  pilotMat = reshape(pilotWave, 1, []);
else
  if size(pilotWave, 1) == numSource
    pilotMat = pilotWave;
  elseif size(pilotWave, 2) == numSource
    pilotMat = pilotWave.';
  else
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidPilotWaveSize', ...
      'pilotWave must be a vector or a numSource x T matrix.');
  end
end

if size(pilotMat, 1) == 1 && numSource > 1
  pilotMat = repmat(pilotMat, numSource, 1);
end
if size(pilotMat, 1) ~= numSource
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidPilotRowCount', ...
    'pilotWave must contain one row per source.');
end
if any(~isfinite(real(pilotMat(:)))) || any(~isfinite(imag(pilotMat(:))))
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidPilotWaveValue', ...
    'pilotWave must contain finite values.');
end

pilotLen = size(pilotMat, 2);
if pilotLen > numSample
  error('estimatorDoaDopplerMlePilotSfOpt:PilotTooLong', ...
    'pilotWave length must not exceed the received snapshot length.');
end

pilotPad = zeros(numSource, numSample, 'like', pilotMat);
pilotPad(:, 1:pilotLen) = pilotMat;
pilotEnergy = real(sum(conj(pilotPad) .* pilotPad, 2));
if any(pilotEnergy <= eps)
  error('estimatorDoaDopplerMlePilotSfOpt:ZeroPilotEnergy', ...
    'Each pilotWave row must have non-zero energy.');
end
end

function doagrid = localParseDoaGrid(doaGrid, numSat, sceneUtc)
%LOCALPARSEDOAGRID Normalize doaGrid to a 1xNs cell array.

if iscell(doaGrid)
  doagrid = reshape(doaGrid, 1, []);
elseif isstruct(doaGrid)
  if isscalar(doaGrid)
    doagrid = {doaGrid};
  else
    doagrid = num2cell(reshape(doaGrid, 1, []));
  end
else
  doagrid = {doaGrid};
end

if isempty(doagrid)
  error('estimatorDoaDopplerMlePilotSfOpt:EmptyDoaGrid', ...
    'doaGrid must not be empty.');
end

if numel(doagrid) == 1 && numSat > 1
  doagrid = repmat(doagrid, 1, numSat);
elseif numel(doagrid) < numSat
  doagrid = [doagrid, repmat(doagrid(1), 1, numSat - numel(doagrid))];
elseif numel(doagrid) > numSat
  doagrid = doagrid(1:numSat);
end

for iSat = 1:numSat
  grid = doagrid{iSat};
  if ~isstruct(grid) || ~isfield(grid, 'type') || ~isfield(grid, 'dimension')
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidDoaGrid', ...
      'Each doaGrid entry must be a grid struct generated by genDoaGrid.');
  end
  if grid.dimension ~= 2
    error('estimatorDoaDopplerMlePilotSfOpt:UnsupportedGridDimension', ...
      'Only 2D doaGrid inputs are supported.');
  end

  grid.type = char(lower(string(grid.type)));
  switch grid.type
    case 'angle'
      if ~isfield(grid, 'angleGrid') || isempty(grid.angleGrid)
        error('estimatorDoaDopplerMlePilotSfOpt:MissingAngleGrid', ...
          'Each angle-mode doaGrid entry must contain a non-empty angleGrid.');
      end
      if ~isfield(grid, 'range') || ~isequal(size(grid.range), [2, 2])
        error('estimatorDoaDopplerMlePilotSfOpt:InvalidAngleRange', ...
          'Angle-mode doaGrid.range must have size 2x2.');
      end
      if ~isfield(grid, 'resolution') || isempty(grid.resolution)
        error('estimatorDoaDopplerMlePilotSfOpt:MissingAngleResolution', ...
          'Angle-mode doaGrid must contain resolution for ECI-grid reconstruction.');
      end
      grid.globalFrame = 'eci';

    case 'latlon'
      if ~isfield(grid, 'latlonGrid') || isempty(grid.latlonGrid)
        error('estimatorDoaDopplerMlePilotSfOpt:MissingLatlonGrid', ...
          'Each latlon-mode doaGrid entry must contain latlonGrid.');
      end
      if ~isfield(grid, 'angleGrid') || isempty(grid.angleGrid)
        error('estimatorDoaDopplerMlePilotSfOpt:MissingAngleGrid', ...
          'Each latlon-mode doaGrid entry must contain angleGrid.');
      end
      if ~isfield(grid, 'globalFrame') || isempty(grid.globalFrame)
        grid.globalFrame = 'ecef';
      else
        grid.globalFrame = char(lower(string(grid.globalFrame)));
      end
      if ~ismember(grid.globalFrame, {'ecef', 'eci'})
        error('estimatorDoaDopplerMlePilotSfOpt:InvalidGlobalFrame', ...
          'Latlon-mode doaGrid.globalFrame must be ''ecef'' or ''eci''.');
      end
      if (~isfield(grid, 'utc') || isempty(grid.utc)) && ~isempty(sceneUtc)
        grid.utc = sceneUtc;
      end

    otherwise
      error('estimatorDoaDopplerMlePilotSfOpt:InvalidDoaGridType', ...
        'doaGrid.type must be ''angle'' or ''latlon''.');
  end

  doagrid{iSat} = grid;
end
end

function [arrayCell, rotMatCell] = localBroadcastSceneContainers(arrayCell, rotMatCell, numSat)
%LOCALBROADCASTSCENECONTAINERS Broadcast singleton scene containers.
%
% This keeps duplicate-reference diagnostic cases usable even when helper
% builders collapse repeated scene metadata to one array / rotation entry.

if numSat <= 1
  return;
end

if numel(arrayCell) == 1
  arrayCell = repmat(arrayCell, 1, numSat);
end

if numel(rotMatCell) == 1
  rotMatCell = repmat(rotMatCell, 1, numSat);
end
end

function [satPosEci, satVelEci] = localBroadcastSatState(satPosEci, satVelEci, numSat)
%LOCALBROADCASTSATSTATE Broadcast singleton satellite states.
%
% Duplicate-reference diagnostic cases may intentionally reuse one physical
% satellite multiple times. In that case, scene.satPosEci / satVelEci can
% arrive with one column while rxSig carries multiple repeated entries.

if numSat <= 1
  return;
end

if size(satPosEci, 2) == 1 && size(satVelEci, 2) == 1
  satPosEci = repmat(satPosEci, 1, numSat);
  satVelEci = repmat(satVelEci, 1, numSat);
end
end

function [doaType, globalFrame, eciAngleGrid, eciDirectionGrid, latlonGrid] = localParseDoaMetadata(doaGrid)
%LOCALPARSEDOAMETADATA Validate shared grid metadata.

baseGrid = doaGrid{1};
doaType = baseGrid.type;
eciAngleGrid = [];
eciDirectionGrid = [];
latlonGrid = [];

for iSat = 2:numel(doaGrid)
  grid = doaGrid{iSat};
  if ~strcmpi(grid.type, doaType)
    error('estimatorDoaDopplerMlePilotSfOpt:GridTypeMismatch', ...
      'All doaGrid entries must share the same grid type.');
  end
  if size(grid.angleGrid, 2) ~= size(baseGrid.angleGrid, 2)
    error('estimatorDoaDopplerMlePilotSfOpt:GridSizeMismatch', ...
      'All doaGrid entries must share the same number of grid points.');
  end
end

switch doaType
  case 'angle'
    globalFrame = 'eci';
    [eciAngleGrid, eciDirectionGrid] = genEciDirectionGrid(baseGrid.resolution, baseGrid.range);
    if size(eciAngleGrid, 2) ~= size(baseGrid.angleGrid, 2)
      error('estimatorDoaDopplerMlePilotSfOpt:AngleGridOrderMismatch', ...
        ['Angle-mode doaGrid must be built from one common ECI direction grid ', ...
         'so that the number of hypotheses matches genEciDirectionGrid(resolution, range).']);
    end

  case 'latlon'
    globalFrame = baseGrid.globalFrame;
    latlonGrid = baseGrid.latlonGrid;
    for iSat = 2:numel(doaGrid)
      grid = doaGrid{iSat};
      if ~strcmpi(grid.globalFrame, globalFrame)
        error('estimatorDoaDopplerMlePilotSfOpt:GlobalFrameMismatch', ...
          'All latlon doaGrid entries must share the same globalFrame.');
      end
      if any(size(grid.latlonGrid) ~= size(latlonGrid)) || any(abs(grid.latlonGrid(:) - latlonGrid(:)) > 1e-10)
        error('estimatorDoaDopplerMlePilotSfOpt:LatlonGridMismatch', ...
          'All latlon doaGrid entries must share the same latlonGrid ordering.');
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidDoaGridType', ...
      'doaGrid.type must be ''angle'' or ''latlon''.');
end
end

function fdRange = localParseFdRange(fdRange, carrierFreq, lightSpeed, satVelEci, refVelEci)
%LOCALPARSEFDRANGE Fill and validate Doppler search range.

if isempty(fdRange)
  wavelength = lightSpeed / carrierFreq;
  satSpeed = max(vecnorm(satVelEci, 2, 1));
  refSpeed = norm(refVelEci);
  fdMax = (satSpeed + refSpeed) / wavelength;
  fdRange = [-fdMax, fdMax];
end

if ~isnumeric(fdRange) || numel(fdRange) ~= 2
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidFdRange', ...
    'fdRange must contain two numeric entries.');
end

fdRange = reshape(fdRange, 1, 2);
if any(~isfinite(fdRange))
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidFdRangeValue', ...
    'fdRange must contain finite values.');
end
if fdRange(1) > fdRange(2)
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidFdRangeOrder', ...
    'fdRange must satisfy fdMin <= fdMax.');
end
end

function arrayCell = localParseSceneArray(scene)
%LOCALPARSESCENEARRAY Convert scene.array to a 1xNs cell array.

if ~isfield(scene, 'array') || isempty(scene.array)
  error('estimatorDoaDopplerMlePilotSfOpt:MissingArrayField', ...
    'scene.array is required.');
end

arrayInput = scene.array;
if iscell(arrayInput)
  arrayCell = reshape(arrayInput, 1, []);
  return;
end

if isstruct(arrayInput)
  if isscalar(arrayInput)
    arrayCell = {arrayInput};
  else
    arrayCell = num2cell(reshape(arrayInput, 1, []));
  end
  return;
end

error('estimatorDoaDopplerMlePilotSfOpt:InvalidArrayField', ...
  'scene.array must be a struct, struct array, or cell array.');
end

function rotMatCell = localParseRotMat(scene)
%LOCALPARSEROTMAT Convert scene.rotMat to a 1xNs cell array.

if ~isfield(scene, 'rotMat') || isempty(scene.rotMat)
  error('estimatorDoaDopplerMlePilotSfOpt:MissingRotMatField', ...
    'scene.rotMat is required.');
end

rotMat = scene.rotMat;
if isnumeric(rotMat)
  if ~isequal(size(rotMat), [3, 3])
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidRotMatSize', ...
      'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
  end
  rotMatCell = {rotMat};
  return;
end

if ~iscell(rotMat) || isempty(rotMat)
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidRotMatType', ...
    'scene.rotMat must be a 3x3 matrix or a cell array of 3x3 matrices.');
end

rotMatCell = reshape(rotMat, 1, []);
for iSat = 1:numel(rotMatCell)
  currentRotMat = rotMatCell{iSat};
  if ~isnumeric(currentRotMat) || ~isequal(size(currentRotMat), [3, 3])
    error('estimatorDoaDopplerMlePilotSfOpt:InvalidRotMatCell', ...
      'Each scene.rotMat entry must be a 3x3 numeric matrix.');
  end
end
end

function [satPosEci, satVelEci] = localParseSatState(scene)
%LOCALPARSESATSTATE Validate satellite ECI states.

if ~isfield(scene, 'satPosEci') || ~isfield(scene, 'satVelEci') || ...
    isempty(scene.satPosEci) || isempty(scene.satVelEci)
  error('estimatorDoaDopplerMlePilotSfOpt:MissingSatState', ...
    'scene.satPosEci and scene.satVelEci are required.');
end

satPosEci = scene.satPosEci;
satVelEci = scene.satVelEci;

if ~isnumeric(satPosEci) || ~isnumeric(satVelEci) || ...
    ~isequal(size(satPosEci), size(satVelEci)) || size(satPosEci, 1) ~= 3
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidSatStateSize', ...
    'scene.satPosEci and scene.satVelEci must both have size 3xNs.');
end
end

function sceneUtc = localParseSceneUtc(scene, doaGrid)
%LOCALPARSESCENEUTC Resolve the scene epoch as a 1x6 UTC vector.

sceneUtc = [];
if isfield(scene, 'utc') && ~isempty(scene.utc)
  sceneUtc = localNormalizeUtc(scene.utc);
  return;
end

if iscell(doaGrid)
  baseGrid = doaGrid{1};
elseif isstruct(doaGrid)
  baseGrid = doaGrid;
else
  baseGrid = struct();
end

if isfield(baseGrid, 'utc') && ~isempty(baseGrid.utc)
  sceneUtc = localNormalizeUtc(baseGrid.utc);
end
end

function utcVec = localNormalizeUtc(utcIn)
%LOCALNORMALIZEUTC Normalize utc to a 1x6 numeric date vector.

if isa(utcIn, 'datetime')
  utcVec = datevec(utcIn);
elseif isnumeric(utcIn) && isequal(size(utcIn), [1, 6])
  utcVec = utcIn;
elseif isnumeric(utcIn) && size(utcIn, 2) == 6 && size(utcIn, 1) >= 1
  utcVec = utcIn(1, :);
else
  error('estimatorDoaDopplerMlePilotSfOpt:InvalidUtc', ...
    'utc must be a datetime scalar, 1x6 date vector, or N-by-6 date matrix.');
end
end

function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct/object field with a default value.

fieldValue = defaultValue;

if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
  return;
end

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
