function [model, lb, ub, modelOpt] = buildDoaDopplerMfModel( ...
  sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, ...
  fdRange, fdRateRange, modelOpt)
%BUILDDOADOPPLERMFMODEL Build the multi-frame DoA-Doppler estimator model.
% This helper extracts the model-construction stage from
% estimatorDoaDopplerMlePilotMfOpt so that scripts and wrappers can reuse
% the same parsing logic, inferred ranges, and DoA bounds.

arguments
  sceneSeq (1,1) struct
  rxSig
  pilotWave
  carrierFreq (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  sampleRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  fdRange = []
  fdRateRange = []
  modelOpt (1,1) struct = struct()
end

modelOpt = localParseModelOpt(modelOpt, sceneSeq);
model = localBuildModel(sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, ...
  doaGrid, fdRange, fdRateRange, modelOpt);
[lb, ub] = localBuildBounds(model);
end

function modelOpt = localParseModelOpt(modelOpt, sceneSeq)
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
if ~isfield(modelOpt, 'steeringMode') || isempty(modelOpt.steeringMode)
  modelOpt.steeringMode = 'framewise';
end
if ~isfield(modelOpt, 'steeringRefFrameIdx') || isempty(modelOpt.steeringRefFrameIdx)
  modelOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
end
if ~isfield(modelOpt, 'useAccessMask') || isempty(modelOpt.useAccessMask)
  modelOpt.useAccessMask = true;
end
if ~isfield(modelOpt, 'frameMask') || isempty(modelOpt.frameMask)
  modelOpt.frameMask = [];
end
if ~isfield(modelOpt, 'optimOpt') || isempty(modelOpt.optimOpt)
  modelOpt.optimOpt = struct();
end
if ~isfield(modelOpt, 'phaseMode') || isempty(modelOpt.phaseMode)
  modelOpt.phaseMode = 'continuous';
end
if ~isfield(modelOpt, 'phaseGridCount') || isempty(modelOpt.phaseGridCount)
  modelOpt.phaseGridCount = 72;
end
if ~isfield(modelOpt, 'phaseRefine') || isempty(modelOpt.phaseRefine)
  modelOpt.phaseRefine = true;
end
if ~isfield(modelOpt, 'fdRateMode') || isempty(modelOpt.fdRateMode)
  modelOpt.fdRateMode = 'unknown';
end
if ~isfield(modelOpt, 'fdRateKnown')
  modelOpt.fdRateKnown = [];
end
if ~isfield(modelOpt, 'enableFdAliasUnwrap') || isempty(modelOpt.enableFdAliasUnwrap)
  modelOpt.enableFdAliasUnwrap = false;
end
if ~isfield(modelOpt, 'fdAliasStepHz')
  modelOpt.fdAliasStepHz = [];
end
if ~isfield(modelOpt, 'fdSatPriorHz')
  modelOpt.fdSatPriorHz = [];
end
if ~isfield(modelOpt, 'timeSecCell')
  modelOpt.timeSecCell = [];
end
if ~isfield(modelOpt, 'timeOffsetSec')
  modelOpt.timeOffsetSec = [];
end
if ~isfield(modelOpt, 'initDoaParam') || isempty(modelOpt.initDoaParam)
  modelOpt.initDoaParam = [];
end
if ~isfield(modelOpt, 'initDoaHalfWidth') || isempty(modelOpt.initDoaHalfWidth)
  modelOpt.initDoaHalfWidth = [];
end
if ~isfield(modelOpt, 'debugEnable') || isempty(modelOpt.debugEnable)
  modelOpt.debugEnable = false;
end
if ~isfield(modelOpt, 'debugTruth') || isempty(modelOpt.debugTruth)
  modelOpt.debugTruth = struct();
end
if ~isfield(modelOpt, 'debugStoreEvalTrace') || isempty(modelOpt.debugStoreEvalTrace)
  modelOpt.debugStoreEvalTrace = false;
end
if ~isfield(modelOpt, 'debugMaxEvalTrace') || isempty(modelOpt.debugMaxEvalTrace)
  modelOpt.debugMaxEvalTrace = 400;
end
if ~isfield(modelOpt, 'disableUnknownWarmAnchor') || isempty(modelOpt.disableUnknownWarmAnchor)
  modelOpt.disableUnknownWarmAnchor = false;
end
if ~isfield(modelOpt, 'freezeDoa') || isempty(modelOpt.freezeDoa)
  modelOpt.freezeDoa = false;
end
if ~isfield(modelOpt, 'disableUnknownDoaReleaseFloor') || isempty(modelOpt.disableUnknownDoaReleaseFloor)
  modelOpt.disableUnknownDoaReleaseFloor = false;
end
if ~isfield(modelOpt, 'unknownDoaReleaseHalfWidth') || isempty(modelOpt.unknownDoaReleaseHalfWidth)
  modelOpt.unknownDoaReleaseHalfWidth = [];
end
if ~isfield(modelOpt, 'continuousPhaseConsistencyWeight') || isempty(modelOpt.continuousPhaseConsistencyWeight)
  modelOpt.continuousPhaseConsistencyWeight = 0;
end
if ~isfield(modelOpt, 'continuousPhaseCollapsePenaltyWeight') || isempty(modelOpt.continuousPhaseCollapsePenaltyWeight)
  modelOpt.continuousPhaseCollapsePenaltyWeight = 0;
end
if ~isfield(modelOpt, 'continuousPhaseNegativeProjectionPenaltyWeight') || isempty(modelOpt.continuousPhaseNegativeProjectionPenaltyWeight)
  modelOpt.continuousPhaseNegativeProjectionPenaltyWeight = 0;
end
if ~isfield(modelOpt, 'continuousPhaseNonRefFitFloorWeight') || isempty(modelOpt.continuousPhaseNonRefFitFloorWeight)
  modelOpt.continuousPhaseNonRefFitFloorWeight = 0;
end
if ~isfield(modelOpt, 'continuousPhaseNonRefSupportFloorWeight') || isempty(modelOpt.continuousPhaseNonRefSupportFloorWeight)
  modelOpt.continuousPhaseNonRefSupportFloorWeight = 0;
end
if ~isfield(modelOpt, 'unknownWarmAnchorUseScaledSolve') || isempty(modelOpt.unknownWarmAnchorUseScaledSolve)
  modelOpt.unknownWarmAnchorUseScaledSolve = true;
end
if ~isfield(modelOpt, 'unknownWarmAnchorFallbackSqp') || isempty(modelOpt.unknownWarmAnchorFallbackSqp)
  modelOpt.unknownWarmAnchorFallbackSqp = true;
end

if ~isscalar(modelOpt.lightSpeed) || ~isfinite(modelOpt.lightSpeed) || modelOpt.lightSpeed <= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidLightSpeed', ...
    'modelOpt.lightSpeed must be a positive finite scalar.');
end
if ~isscalar(modelOpt.initFdCount) || modelOpt.initFdCount <= 0 || ...
    mod(modelOpt.initFdCount, 1) ~= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidInitFdCount', ...
    'modelOpt.initFdCount must be a positive integer scalar.');
end
if ~isstruct(modelOpt.optimOpt)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidOptimOpt', ...
    'modelOpt.optimOpt must be a struct.');
end

if isstring(modelOpt.steeringMode)
  modelOpt.steeringMode = char(modelOpt.steeringMode);
end
modelOpt.steeringMode = lower(strtrim(modelOpt.steeringMode));
if ~ismember(modelOpt.steeringMode, {'frozenref', 'framewise'})
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
    'modelOpt.steeringMode must be ''frozenRef'' or ''framewise''.');
end

if isstring(modelOpt.phaseMode)
  modelOpt.phaseMode = char(modelOpt.phaseMode);
end
modelOpt.phaseMode = lower(strtrim(modelOpt.phaseMode));
if ~ismember(modelOpt.phaseMode, {'continuous', 'relaxed', 'independent'})
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidPhaseMode', ...
    ['modelOpt.phaseMode must be ''continuous'', ''relaxed'', ', ...
     'or ''independent''.']);
end
if ~isscalar(modelOpt.phaseGridCount) || modelOpt.phaseGridCount < 8 || ...
    mod(modelOpt.phaseGridCount, 1) ~= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidPhaseGridCount', ...
    'modelOpt.phaseGridCount must be an integer scalar no smaller than 8.');
end
if ~isscalar(modelOpt.phaseRefine) || ~islogical(modelOpt.phaseRefine)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidPhaseRefine', ...
    'modelOpt.phaseRefine must be a logical scalar.');
end
if ~isscalar(modelOpt.debugEnable) || ~islogical(modelOpt.debugEnable)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDebugEnable', ...
    'modelOpt.debugEnable must be a logical scalar.');
end
if ~isstruct(modelOpt.debugTruth)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDebugTruth', ...
    'modelOpt.debugTruth must be a struct.');
end
if ~isscalar(modelOpt.debugStoreEvalTrace) || ~islogical(modelOpt.debugStoreEvalTrace)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDebugStoreEvalTrace', ...
    'modelOpt.debugStoreEvalTrace must be a logical scalar.');
end
if ~isscalar(modelOpt.debugMaxEvalTrace) || modelOpt.debugMaxEvalTrace <= 0 || ...
    mod(modelOpt.debugMaxEvalTrace, 1) ~= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDebugMaxEvalTrace', ...
    'modelOpt.debugMaxEvalTrace must be a positive integer scalar.');
end
if ~isscalar(modelOpt.disableUnknownWarmAnchor) || ~islogical(modelOpt.disableUnknownWarmAnchor)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDisableUnknownWarmAnchor', ...
    'modelOpt.disableUnknownWarmAnchor must be a logical scalar.');
end
if ~isscalar(modelOpt.freezeDoa) || ~islogical(modelOpt.freezeDoa)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFreezeDoa', ...
    'modelOpt.freezeDoa must be a logical scalar.');
end
if ~isscalar(modelOpt.disableUnknownDoaReleaseFloor) || ~islogical(modelOpt.disableUnknownDoaReleaseFloor)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidDisableUnknownDoaReleaseFloor', ...
    'modelOpt.disableUnknownDoaReleaseFloor must be a logical scalar.');
end
if ~isempty(modelOpt.unknownDoaReleaseHalfWidth)
  if ~isnumeric(modelOpt.unknownDoaReleaseHalfWidth) || ...
      any(~isfinite(modelOpt.unknownDoaReleaseHalfWidth)) || ...
      any(modelOpt.unknownDoaReleaseHalfWidth < 0)
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidUnknownDoaReleaseHalfWidth', ...
      ['modelOpt.unknownDoaReleaseHalfWidth must be empty or a ', ...
       'nonnegative finite 2x1 half-width vector.']);
  end
  if numel(modelOpt.unknownDoaReleaseHalfWidth) == 1
    modelOpt.unknownDoaReleaseHalfWidth = repmat(modelOpt.unknownDoaReleaseHalfWidth, 2, 1);
  end
  if numel(modelOpt.unknownDoaReleaseHalfWidth) ~= 2
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidUnknownDoaReleaseHalfWidthSize', ...
      'modelOpt.unknownDoaReleaseHalfWidth must contain one value per DoA dimension.');
  end
end
if ~isscalar(modelOpt.continuousPhaseConsistencyWeight) || ...
    ~isfinite(modelOpt.continuousPhaseConsistencyWeight) || ...
    modelOpt.continuousPhaseConsistencyWeight < 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidContinuousPhaseConsistencyWeight', ...
    'modelOpt.continuousPhaseConsistencyWeight must be a nonnegative finite scalar.');
end
if ~isscalar(modelOpt.continuousPhaseCollapsePenaltyWeight) || ...
    ~isfinite(modelOpt.continuousPhaseCollapsePenaltyWeight) || ...
    modelOpt.continuousPhaseCollapsePenaltyWeight < 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidContinuousPhaseCollapsePenaltyWeight', ...
    'modelOpt.continuousPhaseCollapsePenaltyWeight must be a nonnegative finite scalar.');
end
if ~isscalar(modelOpt.continuousPhaseNegativeProjectionPenaltyWeight) || ...
    ~isfinite(modelOpt.continuousPhaseNegativeProjectionPenaltyWeight) || ...
    modelOpt.continuousPhaseNegativeProjectionPenaltyWeight < 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidContinuousPhaseNegativeProjectionPenaltyWeight', ...
    ['modelOpt.continuousPhaseNegativeProjectionPenaltyWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~isscalar(modelOpt.continuousPhaseNonRefFitFloorWeight) || ...
    ~isfinite(modelOpt.continuousPhaseNonRefFitFloorWeight) || ...
    modelOpt.continuousPhaseNonRefFitFloorWeight < 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidContinuousPhaseNonRefFitFloorWeight', ...
    ['modelOpt.continuousPhaseNonRefFitFloorWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~isscalar(modelOpt.continuousPhaseNonRefSupportFloorWeight) || ...
    ~isfinite(modelOpt.continuousPhaseNonRefSupportFloorWeight) || ...
    modelOpt.continuousPhaseNonRefSupportFloorWeight < 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidContinuousPhaseNonRefSupportFloorWeight', ...
    ['modelOpt.continuousPhaseNonRefSupportFloorWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~isscalar(modelOpt.unknownWarmAnchorUseScaledSolve) || ...
    ~islogical(modelOpt.unknownWarmAnchorUseScaledSolve)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidUnknownWarmAnchorUseScaledSolve', ...
    'modelOpt.unknownWarmAnchorUseScaledSolve must be a logical scalar.');
end
if ~isscalar(modelOpt.unknownWarmAnchorFallbackSqp) || ...
    ~islogical(modelOpt.unknownWarmAnchorFallbackSqp)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidUnknownWarmAnchorFallbackSqp', ...
    'modelOpt.unknownWarmAnchorFallbackSqp must be a logical scalar.');
end
if ~modelOpt.debugEnable
  modelOpt.debugStoreEvalTrace = false;
end

if isstring(modelOpt.fdRateMode)
  modelOpt.fdRateMode = char(modelOpt.fdRateMode);
end
modelOpt.fdRateMode = lower(strtrim(modelOpt.fdRateMode));
if ~ismember(modelOpt.fdRateMode, {'unknown', 'known', 'zero'})
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
    'modelOpt.fdRateMode must be ''unknown'', ''known'', or ''zero''.');
end

switch modelOpt.fdRateMode
  case 'unknown'
    if ~isempty(modelOpt.fdRateKnown) && ...
        (~isscalar(modelOpt.fdRateKnown) || ~isfinite(modelOpt.fdRateKnown))
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateKnown', ...
        'modelOpt.fdRateKnown must be empty or a finite scalar.');
    end

  case 'known'
    if ~isscalar(modelOpt.fdRateKnown) || ~isfinite(modelOpt.fdRateKnown)
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateKnown', ...
        'modelOpt.fdRateKnown must be a finite scalar when fdRateMode = ''known''.');
    end

  case 'zero'
    modelOpt.fdRateKnown = 0;
end

if ~isscalar(modelOpt.enableFdAliasUnwrap) || ~islogical(modelOpt.enableFdAliasUnwrap)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidEnableFdAliasUnwrap', ...
    'modelOpt.enableFdAliasUnwrap must be a logical scalar.');
end
if ~(isempty(modelOpt.fdAliasStepHz) || ...
    (isscalar(modelOpt.fdAliasStepHz) && isfinite(modelOpt.fdAliasStepHz) && ...
     modelOpt.fdAliasStepHz > 0))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdAliasStepHz', ...
    'modelOpt.fdAliasStepHz must be empty or a positive finite scalar.');
end
if ~(isempty(modelOpt.fdSatPriorHz) || isnumeric(modelOpt.fdSatPriorHz))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdSatPriorHz', ...
    'modelOpt.fdSatPriorHz must be empty or numeric.');
end
if ~isempty(modelOpt.fdSatPriorHz) && numel(modelOpt.fdSatPriorHz) ~= sceneSeq.numSat
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdSatPriorHzSize', ...
    'modelOpt.fdSatPriorHz must contain one Doppler prior per satellite.');
end
if ~(isempty(modelOpt.timeSecCell) || iscell(modelOpt.timeSecCell))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeSecCell', ...
    'modelOpt.timeSecCell must be empty or a cell array.');
end
if ~(isempty(modelOpt.timeOffsetSec) || isnumeric(modelOpt.timeOffsetSec))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeOffsetSec', ...
    'modelOpt.timeOffsetSec must be empty or numeric.');
end
if ~isempty(modelOpt.initDoaParam)
  if ~isnumeric(modelOpt.initDoaParam) || numel(modelOpt.initDoaParam) ~= 2 || ...
      any(~isfinite(modelOpt.initDoaParam))
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidInitDoaParam', ...
      'modelOpt.initDoaParam must be empty or a finite 2x1 parameter vector.');
  end
end
if ~isempty(modelOpt.initDoaHalfWidth)
  if ~isnumeric(modelOpt.initDoaHalfWidth) || numel(modelOpt.initDoaHalfWidth) ~= 2 || ...
      any(~isfinite(modelOpt.initDoaHalfWidth)) || any(modelOpt.initDoaHalfWidth < 0)
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidInitDoaHalfWidth', ...
      ['modelOpt.initDoaHalfWidth must be empty or a nonnegative ', ...
       'finite 2x1 half-width vector.']);
  end
  if isempty(modelOpt.initDoaParam)
    error('estimatorDoaDopplerMlePilotMfOpt:MissingInitDoaParam', ...
      ['modelOpt.initDoaHalfWidth requires modelOpt.initDoaParam so ', ...
       'that the local multi-frame DoA box is anchored to a real coarse ', ...
       'initializer.']);
  end
end

if ~isscalar(modelOpt.steeringRefFrameIdx) || ...
    ~isfinite(modelOpt.steeringRefFrameIdx) || ...
    mod(modelOpt.steeringRefFrameIdx, 1) ~= 0 || ...
    modelOpt.steeringRefFrameIdx < 1 || ...
    modelOpt.steeringRefFrameIdx > sceneSeq.numFrame
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringRefFrameIdx', ...
    'modelOpt.steeringRefFrameIdx must be an integer within [1, numFrame].');
end

if ~isempty(modelOpt.frameMask)
  modelOpt.frameMask = localNormalizeFrameMask(modelOpt.frameMask, sceneSeq.numSat, sceneSeq.numFrame, ...
    'modelOpt.frameMask');
end
end

function model = localBuildModel(sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, ...
  doaGridIn, fdRange, fdRateRange, modelOpt)
%LOCALBUILDMODEL Build fixed quantities used by the dynamic estimator.

localCheckSceneSeq(sceneSeq);

numFrame = sceneSeq.numFrame;
numSat = sceneSeq.numSat;

arrayCell = localParseSceneArray(sceneSeq);
[numElement, rotMatCell, satPosEci, satVelEci] = ...
  localParseSceneKinematics(sceneSeq, arrayCell);
refScene = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[refState, refSatIdxLocal] = resolveReferenceSatState( ...
  refScene, satPosEci(:, :, sceneSeq.refFrameIdx), satVelEci(:, :, sceneSeq.refFrameIdx));
[refSatPosEci, refSatVelEci] = localExtractRefSatKinematics(satPosEci, satVelEci, refSatIdxLocal);

[rxCell, numSampleFrame] = localParseRxSig(rxSig, numFrame, numSat, numElement);
pilotPadCell = localParsePilotWave(pilotWave, numFrame, numSampleFrame);

doaGrid = localParseDoaGrid(doaGridIn, numSat);
[doaType, globalFrame, eciAngleGrid, eciDirectionGrid, latlonGrid] = localParseDoaMetadata(doaGrid);

frameMask = localParseFrameMask(sceneSeq, modelOpt, numSat, numFrame);

sampleCovAgg = localBuildAggSampleCov(rxCell, frameMask);
fdRange = localParseFdRange(fdRange, carrierFreq, modelOpt.lightSpeed, satVelEci, refSatVelEci);
fdRateRange = localParseFdRateRange(fdRateRange, carrierFreq, modelOpt.lightSpeed, ...
  satVelEci, refSatVelEci, sceneSeq.timeOffsetSec);

timeSecCell = localParseTimeSecCell(modelOpt.timeSecCell, numFrame, numSampleFrame, sampleRate);
timeOffsetSec = localParseTimeOffsetSec(modelOpt.timeOffsetSec, sceneSeq.timeOffsetSec, numFrame);

model = struct();
model.array = arrayCell;
model.rotMatCell = rotMatCell;
model.rxSig = rxCell;
model.sampleCovAgg = sampleCovAgg;
model.pilotPadCell = pilotPadCell;

model.numSat = numSat;
model.numFrame = numFrame;
model.wavelength = modelOpt.lightSpeed / carrierFreq;

model.timeSecCell = timeSecCell;
model.timeOffsetSec = timeOffsetSec;
model.userStateRef = localBuildUserStateRef(sceneSeq, timeOffsetSec);

model.useLogObjective = logical(modelOpt.useLogObjective);
model.steeringMode = modelOpt.steeringMode;
model.steeringRefFrameIdx = modelOpt.steeringRefFrameIdx;
model.frameMask = frameMask;

model.phaseMode = modelOpt.phaseMode;
model.phaseGridCount = modelOpt.phaseGridCount;
model.phaseRefine = modelOpt.phaseRefine;

model.doaGrid = doaGrid;
model.doaType = doaType;
model.globalFrame = globalFrame;
model.eciAngleGrid = eciAngleGrid;
model.eciDirectionGrid = eciDirectionGrid;
model.latlonGrid = latlonGrid;

model.fdRange = fdRange;
model.fdRateRange = fdRateRange;
model.fdRateMode = modelOpt.fdRateMode;
if ismember(model.fdRateMode, {'known', 'zero'})
  model.fdRateKnown = modelOpt.fdRateKnown;
else
  model.fdRateKnown = [];
end
model.enableFdAliasUnwrap = logical(modelOpt.enableFdAliasUnwrap);
model.fdAliasStepHz = modelOpt.fdAliasStepHz;
model.fdSatPriorHz = modelOpt.fdSatPriorHz;

model.utcVec = sceneSeq.utcVec(:);
model.sceneCell = sceneSeq.sceneCell;
model.refScene = refScene;
model.refFrameIdx = sceneSeq.refFrameIdx;
model.refSatIdxLocal = refSatIdxLocal;
model.refSatIdxGlobal = localGetStructField(refState, 'satIdxGlobal', NaN);
model.refStateSource = char(localGetStructField(refState, 'source', ""));
model.satPosEci = satPosEci;
model.satVelEci = satVelEci;
model.refPosEci = refSatPosEci;
model.refVelEci = refSatVelEci;

model.initFdCount = modelOpt.initFdCount;
model.initDoaParam = reshape(modelOpt.initDoaParam, [], 1);
model.initDoaHalfWidth = reshape(modelOpt.initDoaHalfWidth, [], 1);
model.optimOpt = modelOpt.optimOpt;
model.debugEnable = modelOpt.debugEnable;
model.debugTruth = modelOpt.debugTruth;
model.debugStoreEvalTrace = modelOpt.debugStoreEvalTrace;
model.debugMaxEvalTrace = modelOpt.debugMaxEvalTrace;
model.disableUnknownWarmAnchor = logical(modelOpt.disableUnknownWarmAnchor);
model.freezeDoa = logical(modelOpt.freezeDoa);
model.disableUnknownDoaReleaseFloor = logical(modelOpt.disableUnknownDoaReleaseFloor);
model.unknownDoaReleaseHalfWidth = reshape(modelOpt.unknownDoaReleaseHalfWidth, [], 1);
model.continuousPhaseConsistencyWeight = modelOpt.continuousPhaseConsistencyWeight;
model.continuousPhaseCollapsePenaltyWeight = modelOpt.continuousPhaseCollapsePenaltyWeight;
model.continuousPhaseNegativeProjectionPenaltyWeight = modelOpt.continuousPhaseNegativeProjectionPenaltyWeight;
model.continuousPhaseNonRefFitFloorWeight = modelOpt.continuousPhaseNonRefFitFloorWeight;
model.continuousPhaseNonRefSupportFloorWeight = modelOpt.continuousPhaseNonRefSupportFloorWeight;
model.unknownWarmAnchorUseScaledSolve = logical(modelOpt.unknownWarmAnchorUseScaledSolve);
model.unknownWarmAnchorFallbackSqp = logical(modelOpt.unknownWarmAnchorFallbackSqp);
[doaLbCache, doaUbCache] = localBuildDoaBounds(model);
model.doaLb = doaLbCache;
model.doaUb = doaUbCache;
model.lb = [doaLbCache; model.fdRange(1)];
model.ub = [doaUbCache; model.fdRange(2)];
if strcmp(model.fdRateMode, 'unknown')
  model.lb = [model.lb; model.fdRateRange(1)];
  model.ub = [model.ub; model.fdRateRange(2)];
end
model.cachedBoundsReady = true;
end

function localCheckSceneSeq(sceneSeq)
%LOCALCHECKSCENESEQ Check the basic validity of sceneSeq.

requiredField = {'numFrame', 'numSat', 'sceneCell', 'utcVec', 'timeOffsetSec', 'refFrameIdx'};
for iField = 1:numel(requiredField)
  if ~isfield(sceneSeq, requiredField{iField})
    error('estimatorDoaDopplerMlePilotMfOpt:MissingSceneSeqField', ...
      'sceneSeq.%s is required.', requiredField{iField});
  end
end

if ~isscalar(sceneSeq.numFrame) || sceneSeq.numFrame < 1 || mod(sceneSeq.numFrame, 1) ~= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidNumFrame', ...
    'sceneSeq.numFrame must be a positive integer scalar.');
end
if ~isscalar(sceneSeq.numSat) || sceneSeq.numSat < 1 || mod(sceneSeq.numSat, 1) ~= 0
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidNumSat', ...
    'sceneSeq.numSat must be a positive integer scalar.');
end
if ~iscell(sceneSeq.sceneCell) || numel(sceneSeq.sceneCell) ~= sceneSeq.numFrame
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidSceneCell', ...
    'sceneSeq.sceneCell must contain one scene per frame.');
end
if ~isa(sceneSeq.utcVec, 'datetime') || numel(sceneSeq.utcVec) ~= sceneSeq.numFrame
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidUtcVec', ...
    'sceneSeq.utcVec must be a datetime vector with one entry per frame.');
end
if numel(sceneSeq.timeOffsetSec) ~= sceneSeq.numFrame
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeOffset', ...
    'sceneSeq.timeOffsetSec must contain one value per frame.');
end
end

function arrayCell = localParseSceneArray(sceneSeq)
%LOCALPARSESCENEARRAY Parse array geometry from the reference frame scene.

refScene = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

if ~isfield(refScene, 'array') || isempty(refScene.array)
  error('estimatorDoaDopplerMlePilotMfOpt:MissingArrayField', ...
    'sceneSeq.sceneCell{refFrameIdx}.array is required.');
end

arrayInput = refScene.array;
if iscell(arrayInput)
  arrayCell = reshape(arrayInput, 1, []);
elseif isstruct(arrayInput)
  if isscalar(arrayInput)
    arrayCell = {arrayInput};
  else
    arrayCell = num2cell(reshape(arrayInput, 1, []));
  end
else
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidArrayField', ...
    'scene.array must be a struct, struct array, or cell array.');
end

if numel(arrayCell) ~= sceneSeq.numSat
  error('estimatorDoaDopplerMlePilotMfOpt:ArrayCountMismatch', ...
    'The number of arrays must match sceneSeq.numSat.');
end
end

function [numElement, rotMatCell, satPosEci, satVelEci] = ...
  localParseSceneKinematics(sceneSeq, arrayCell)
%LOCALPARSESCENEKINEMATICS Parse frame-wise rotations and ECI kinematics.

numSat = sceneSeq.numSat;
numFrame = sceneSeq.numFrame;

numElement = zeros(numSat, 1);
for iSat = 1:numSat
  if ~isstruct(arrayCell{iSat}) || ~isfield(arrayCell{iSat}, 'positions')
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidArray', ...
      'Each scene.array entry must be a valid array struct.');
  end
  numElement(iSat) = size(arrayCell{iSat}.positions, 2);
end

rotMatCell = cell(1, numFrame);
satPosEci = zeros(3, numSat, numFrame);
satVelEci = zeros(3, numSat, numFrame);
for iFrame = 1:numFrame
  sceneTmp = sceneSeq.sceneCell{iFrame};

  if ~isfield(sceneTmp, 'rotMat') || ~iscell(sceneTmp.rotMat) || numel(sceneTmp.rotMat) ~= numSat
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidRotMat', ...
      'Each frame scene must contain a 1xNs cell array scene.rotMat.');
  end
  rotMatCell{iFrame} = reshape(sceneTmp.rotMat, 1, []);

  if ~isfield(sceneTmp, 'satPosEci') || ~isequal(size(sceneTmp.satPosEci), [3, numSat])
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSatPos', ...
      'Each frame scene must contain scene.satPosEci with size 3xNs.');
  end
  if ~isfield(sceneTmp, 'satVelEci') || ~isequal(size(sceneTmp.satVelEci), [3, numSat])
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSatVel', ...
      'Each frame scene must contain scene.satVelEci with size 3xNs.');
  end
  satPosEci(:, :, iFrame) = sceneTmp.satPosEci;
  satVelEci(:, :, iFrame) = sceneTmp.satVelEci;
end
end

function [refSatPosEci, refSatVelEci] = localExtractRefSatKinematics(satPosEci, satVelEci, refSatIdxLocal)
%LOCALEXTRACTREFSATKINEMATICS Extract one frame-wise reference-satellite state.

refSatPosEci = reshape(satPosEci(:, refSatIdxLocal, :), 3, []);
refSatVelEci = reshape(satVelEci(:, refSatIdxLocal, :), 3, []);
end

function [rxCell, numSampleFrame] = localParseRxSig(rxSig, numFrame, numSat, numElement)
%LOCALPARSERXSIG Normalize rxSig to a frame-by-frame cell container.

if iscell(rxSig)
  if numFrame == 1 && ~(iscell(rxSig{1}) || isnumeric(rxSig{1}))
    frameInput = {rxSig};
  else
    if numel(rxSig) ~= numFrame
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidRxSigFrameCount', ...
        'rxSig must contain exactly one entry per frame.');
    end
    frameInput = reshape(rxSig, 1, []);
  end
else
  if numFrame ~= 1
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidRxSigType', ...
      'For multi-frame input, rxSig must be a 1xNf cell array.');
  end
  frameInput = {rxSig};
end

rxCell = cell(1, numFrame);
numSampleFrame = zeros(1, numFrame);

for iFrame = 1:numFrame
  currentFrame = frameInput{iFrame};

  if iscell(currentFrame)
    currentCell = reshape(currentFrame, 1, []);
  elseif isnumeric(currentFrame)
    currentCell = {currentFrame};
  else
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidRxSigEntry', ...
      'Each frame rxSig entry must be a numeric matrix or a cell array of matrices.');
  end

  if numel(currentCell) ~= numSat
    error('estimatorDoaDopplerMlePilotMfOpt:RxSigSatCountMismatch', ...
      'Frame %d rxSig entry count must match sceneSeq.numSat.', iFrame);
  end

  currentLen = [];

  for iSat = 1:numSat
    currentRx = currentCell{iSat};

    if ~isnumeric(currentRx) || ~ismatrix(currentRx) || isempty(currentRx)
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidRxSigMatrix', ...
        'Each rxSig block must be a non-empty numeric matrix.');
    end
    if any(~isfinite(real(currentRx(:)))) || any(~isfinite(imag(currentRx(:))))
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidRxSigValue', ...
        'Each rxSig block must contain finite values.');
    end
    if size(currentRx, 1) ~= numElement(iSat)
      error('estimatorDoaDopplerMlePilotMfOpt:RxElementMismatch', ...
        'Frame %d, satellite %d row size must match the array element count.', ...
        iFrame, iSat);
    end

    if isempty(currentLen)
      currentLen = size(currentRx, 2);
    elseif size(currentRx, 2) ~= currentLen
      error('estimatorDoaDopplerMlePilotMfOpt:RxLengthMismatch', ...
        'Within one frame, all satellite blocks must share the same sample length.');
    end

  end

  rxCell{iFrame} = currentCell;
  numSampleFrame(iFrame) = currentLen;
end
end

function pilotPadCell = localParsePilotWave(pilotWave, numFrame, numSampleFrame)
%LOCALPARSEPILOTWAVE Parse shared or frame-wise pilot templates.

if iscell(pilotWave)
  if numel(pilotWave) ~= numFrame
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidPilotFrameCount', ...
      'When pilotWave is a cell array, it must contain exactly one entry per frame.');
  end
  pilotInput = reshape(pilotWave, 1, []);
else
  pilotInput = repmat({pilotWave}, 1, numFrame);
end

pilotPadCell = cell(1, numFrame);

for iFrame = 1:numFrame
  currentPilot = pilotInput{iFrame};

  if ~isnumeric(currentPilot) || isempty(currentPilot)
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidPilotWave', ...
      'Each pilotWave entry must be a non-empty numeric vector or matrix.');
  end

  if isvector(currentPilot)
    pilotRow = reshape(currentPilot, 1, []);
  else
    if size(currentPilot, 1) == 1
      pilotRow = currentPilot;
    elseif size(currentPilot, 2) == 1
      pilotRow = currentPilot.';
    else
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidPilotWaveSize', ...
        'This dynamic estimator only supports single-source pilotWave inputs.');
    end
  end

  if any(~isfinite(real(pilotRow(:)))) || any(~isfinite(imag(pilotRow(:))))
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidPilotWaveValue', ...
      'pilotWave must contain finite values.');
  end

  pilotLen = numel(pilotRow);
  if pilotLen > numSampleFrame(iFrame)
    error('estimatorDoaDopplerMlePilotMfOpt:PilotTooLong', ...
      'pilotWave length must not exceed the received snapshot length.');
  end

  pilotPad = zeros(1, numSampleFrame(iFrame), 'like', pilotRow);
  pilotPad(1:pilotLen) = pilotRow;

  pilotPadCell{iFrame} = pilotPad;
  pilotEnergy = real(pilotPad * pilotPad');
  if pilotEnergy <= eps
    error('estimatorDoaDopplerMlePilotMfOpt:ZeroPilotEnergy', ...
      'Each pilotWave entry must have non-zero energy.');
  end
end
end

function doaGrid = localParseDoaGrid(doaGrid, numSat)
%LOCALPARSEDOAGRID Normalize doaGrid to a 1xNs cell array.

if iscell(doaGrid)
  doaGrid = reshape(doaGrid, 1, []);
else
  if numSat ~= 1
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaGridContainer', ...
      'For multi-satellite input, doaGrid must be a cell array with one grid per satellite.');
  end
  doaGrid = {doaGrid};
end

if numel(doaGrid) ~= numSat
  error('estimatorDoaDopplerMlePilotMfOpt:DoaGridCountMismatch', ...
    'The number of doaGrid entries must match the number of satellites.');
end

for iSat = 1:numSat
  grid = doaGrid{iSat};
  if ~isstruct(grid) || ~isfield(grid, 'type') || ~isfield(grid, 'dimension')
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaGrid', ...
      'Each doaGrid entry must be a grid struct generated by genDoaGrid.');
  end
  if grid.dimension ~= 2
    error('estimatorDoaDopplerMlePilotMfOpt:UnsupportedGridDimension', ...
      'Only 2D doaGrid inputs are supported.');
  end

  grid.type = char(lower(string(grid.type)));
  switch grid.type
    case 'angle'
      if ~isfield(grid, 'angleGrid') || isempty(grid.angleGrid)
        error('estimatorDoaDopplerMlePilotMfOpt:MissingAngleGrid', ...
          'Each angle-mode doaGrid entry must contain angleGrid.');
      end
      if ~isfield(grid, 'range') || ~isequal(size(grid.range), [2, 2])
        error('estimatorDoaDopplerMlePilotMfOpt:InvalidAngleRange', ...
          'Angle-mode doaGrid.range must have size 2x2.');
      end
      if ~isfield(grid, 'resolution') || isempty(grid.resolution)
        error('estimatorDoaDopplerMlePilotMfOpt:MissingAngleResolution', ...
          'Angle-mode doaGrid must contain resolution.');
      end
      grid.globalFrame = 'eci';

    case 'latlon'
      if ~isfield(grid, 'latlonGrid') || isempty(grid.latlonGrid)
        error('estimatorDoaDopplerMlePilotMfOpt:MissingLatlonGrid', ...
          'Each latlon-mode doaGrid entry must contain latlonGrid.');
      end
      if ~isfield(grid, 'angleGrid') || isempty(grid.angleGrid)
        error('estimatorDoaDopplerMlePilotMfOpt:MissingAngleGrid', ...
          'Each latlon-mode doaGrid entry must contain angleGrid.');
      end
      if ~isfield(grid, 'range') || ~isequal(size(grid.range), [2, 2])
        error('estimatorDoaDopplerMlePilotMfOpt:InvalidLatlonRange', ...
          'Latlon-mode doaGrid.range must have size 2x2.');
      end
      if ~isfield(grid, 'globalFrame') || isempty(grid.globalFrame)
        grid.globalFrame = 'ecef';
      else
        grid.globalFrame = char(lower(string(grid.globalFrame)));
      end
      if ~ismember(grid.globalFrame, {'ecef', 'eci'})
        error('estimatorDoaDopplerMlePilotMfOpt:InvalidGlobalFrame', ...
          'Latlon-mode doaGrid.globalFrame must be ''ecef'' or ''eci''.');
      end

    otherwise
      error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaGridType', ...
        'doaGrid.type must be ''angle'' or ''latlon''.');
  end

  doaGrid{iSat} = grid;
end
end

function [doaType, globalFrame, eciAngleGrid, eciDirectionGrid, latlonGrid] = ...
  localParseDoaMetadata(doaGrid)
%LOCALPARSEDOAMETADATA Validate shared grid metadata.

baseGrid = doaGrid{1};
doaType = baseGrid.type;
eciAngleGrid = [];
eciDirectionGrid = [];
latlonGrid = [];

for iSat = 2:numel(doaGrid)
  grid = doaGrid{iSat};
  if ~strcmpi(grid.type, doaType)
    error('estimatorDoaDopplerMlePilotMfOpt:GridTypeMismatch', ...
      'All doaGrid entries must share the same grid type.');
  end
  if size(grid.angleGrid, 2) ~= size(baseGrid.angleGrid, 2)
    error('estimatorDoaDopplerMlePilotMfOpt:GridSizeMismatch', ...
      'All doaGrid entries must share the same number of grid points.');
  end
end

switch doaType
  case 'angle'
    globalFrame = 'eci';
    [eciAngleGrid, eciDirectionGrid] = genEciDirectionGrid(baseGrid.resolution, baseGrid.range);
    if size(eciAngleGrid, 2) ~= size(baseGrid.angleGrid, 2)
      error('estimatorDoaDopplerMlePilotMfOpt:AngleGridOrderMismatch', ...
        ['Angle-mode doaGrid must be built from one common ECI direction grid ', ...
         'so that its hypothesis count matches genEciDirectionGrid(resolution, range).']);
    end

  case 'latlon'
    globalFrame = baseGrid.globalFrame;
    latlonGrid = baseGrid.latlonGrid;
    for iSat = 2:numel(doaGrid)
      grid = doaGrid{iSat};
      if ~strcmpi(grid.globalFrame, globalFrame)
        error('estimatorDoaDopplerMlePilotMfOpt:GlobalFrameMismatch', ...
          'All latlon doaGrid entries must share the same globalFrame.');
      end
      if any(size(grid.latlonGrid) ~= size(latlonGrid)) || ...
          any(abs(grid.latlonGrid(:) - latlonGrid(:)) > 1e-10)
        error('estimatorDoaDopplerMlePilotMfOpt:LatlonGridMismatch', ...
          'All latlon doaGrid entries must share the same latlonGrid ordering.');
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', doaType);
end
end

function frameMask = localParseFrameMask(sceneSeq, modelOpt, numSat, numFrame)
%LOCALPARSEFRAMEMASK Parse the frame-satellite validity mask.

if ~isempty(modelOpt.frameMask)
  frameMask = localNormalizeFrameMask(modelOpt.frameMask, numSat, numFrame, 'modelOpt.frameMask');
  return;
end

frameMask = true(numSat, numFrame);

if ~modelOpt.useAccessMask
  return;
end

if ~isfield(sceneSeq, 'access') || isempty(sceneSeq.access)
  return;
end

access = sceneSeq.access;
accessSize = size(access);

if isequal(accessSize, [numSat, 1, numFrame]) || ...
    (numel(accessSize) == 3 && accessSize(1) == numSat && accessSize(2) == 1 && accessSize(3) == numFrame)
  frameMask = reshape(access(:, 1, :), numSat, numFrame);
elseif isequal(accessSize, [numSat, numFrame])
  frameMask = access;
elseif isvector(access) && numSat == 1 && numel(access) == numFrame
  frameMask = reshape(access, 1, numFrame);
elseif numel(accessSize) == 3 && accessSize(1) == numSat && accessSize(3) == numFrame
  if accessSize(2) ~= 1
    error('estimatorDoaDopplerMlePilotMfOpt:UnsupportedMultiUserAccess', ...
      'This dynamic estimator is single-source only and expects one user in sceneSeq.access.');
  end
  frameMask = reshape(access(:, 1, :), numSat, numFrame);
else
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidAccessSize', ...
    'sceneSeq.access has an unsupported size.');
end

frameMask = localNormalizeFrameMask(frameMask, numSat, numFrame, 'sceneSeq.access');
end

function frameMask = localNormalizeFrameMask(frameMaskIn, numSat, numFrame, srcName)
%LOCALNORMALIZEFRAMEMASK Normalize one frame mask to size Ns x Nf.

frameMask = logical(frameMaskIn);

if isequal(size(frameMask), [numSat, numFrame])
  return;
end

if numSat == 1 && isvector(frameMask) && numel(frameMask) == numFrame
  frameMask = reshape(frameMask, 1, numFrame);
  return;
end

if numFrame == 1 && isvector(frameMask) && numel(frameMask) == numSat
  frameMask = reshape(frameMask, numSat, 1);
  return;
end

error('estimatorDoaDopplerMlePilotMfOpt:InvalidFrameMaskSize', ...
  '%s must have size Ns x Nf.', srcName);
end

function sampleCovAgg = localBuildAggSampleCov(rxCell, frameMask)
%LOCALBUILDAGGSAMPLECOV Build frame-aggregated sample covariance matrices.

numFrame = numel(rxCell);
numSat = numel(rxCell{1});
sampleCovAgg = cell(1, numSat);

for iSat = 1:numSat
  rSum = [];
  count = 0;

  for iFrame = 1:numFrame
    if ~frameMask(iSat, iFrame)
      continue;
    end

    currentRx = rxCell{iFrame}{iSat};
    if isempty(rSum)
      rSum = currentRx * currentRx';
    else
      rSum = rSum + currentRx * currentRx';
    end
    count = count + size(currentRx, 2);
  end

  if isempty(rSum) || count == 0
    currentM = size(rxCell{1}{iSat}, 1);
    sampleCovAgg{iSat} = zeros(currentM, currentM);
  else
    sampleCovAgg{iSat} = 0.5 * ((rSum / count) + (rSum / count)');
  end
end
end

function fdRange = localParseFdRange(fdRange, carrierFreq, lightSpeed, satVelEci, refVelEci)
%LOCALPARSEFDRANGE Fill and validate Doppler intercept search range.

if isempty(fdRange)
  wavelength = lightSpeed / carrierFreq;
  satSpeed = max(reshape(vecnorm(satVelEci, 2, 1), 1, []), [], 'omitnan');
  refSpeed = max(vecnorm(refVelEci, 2, 1));
  fdMax = (satSpeed + refSpeed) / wavelength;
  fdRange = [-fdMax, fdMax];
end

if ~isnumeric(fdRange) || numel(fdRange) ~= 2
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRange', ...
    'fdRange must contain two numeric entries.');
end

fdRange = reshape(fdRange, 1, 2);
if any(~isfinite(fdRange))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRangeValue', ...
    'fdRange must contain finite values.');
end
if fdRange(1) > fdRange(2)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRangeOrder', ...
    'fdRange must satisfy fdMin <= fdMax.');
end
end

function fdRateRange = localParseFdRateRange(fdRateRange, carrierFreq, lightSpeed, ...
  satVelEci, refVelEci, timeOffsetSec)
%LOCALPARSEFDRATERANGE Fill and validate Doppler-rate search range.

if isempty(fdRateRange)
  if numel(timeOffsetSec) < 2 || max(timeOffsetSec) - min(timeOffsetSec) <= 0
    fdRateRange = [0, 0];
  else
    wavelength = lightSpeed / carrierFreq;
    dt = diff(timeOffsetSec(:).');
    dt(dt == 0) = [];

    if isempty(dt)
      fdRateRange = [0, 0];
    else
      satAcc = diff(satVelEci, 1, 3);
      satAcc = satAcc ./ reshape(dt, 1, 1, []);
      satAccMax = max(reshape(vecnorm(satAcc, 2, 1), 1, []), [], 'omitnan');

      refAcc = diff(refVelEci, 1, 2);
      refAcc = refAcc ./ reshape(dt, 1, []);
      refAccMax = max(vecnorm(refAcc, 2, 1), [], 'omitnan');

      fdRateMax = 2 * (satAccMax + refAccMax) / wavelength;
      if ~isfinite(fdRateMax) || fdRateMax <= 0
        fdRateMax = 0;
      end
      fdRateRange = [-fdRateMax, fdRateMax];
    end
  end
end

if ~isnumeric(fdRateRange) || numel(fdRateRange) ~= 2
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateRange', ...
    'fdRateRange must contain two numeric entries.');
end

fdRateRange = reshape(fdRateRange, 1, 2);
if any(~isfinite(fdRateRange))
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateRangeValue', ...
    'fdRateRange must contain finite values.');
end
if fdRateRange(1) > fdRateRange(2)
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateRangeOrder', ...
    'fdRateRange must satisfy fdRateMin <= fdRateMax.');
end
end

function timeSecCell = localParseTimeSecCell(timeSecCellIn, numFrame, numSampleFrame, sampleRate)
%LOCALPARSETIMESECCELL Parse optional per-frame local time vectors.

if isempty(timeSecCellIn)
  timeSecCell = cell(1, numFrame);
  for iFrame = 1:numFrame
    timeSecCell{iFrame} = (0:numSampleFrame(iFrame)-1) / sampleRate;
  end
  return;
end

if ~iscell(timeSecCellIn) || numel(timeSecCellIn) ~= numFrame
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeSecCell', ...
    'modelOpt.timeSecCell must contain one cell per frame.');
end

timeSecCell = cell(1, numFrame);
for iFrame = 1:numFrame
  timeVec = timeSecCellIn{iFrame};
  if ~isnumeric(timeVec) || numel(timeVec) ~= numSampleFrame(iFrame) || ...
      any(~isfinite(timeVec))
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeSecCellEntry', ...
      'modelOpt.timeSecCell{%d} must match the pilot length and be finite.', iFrame);
  end
  timeSecCell{iFrame} = reshape(timeVec, 1, []);
end
end

function timeOffsetSec = localParseTimeOffsetSec(timeOffsetSecIn, defaultTimeOffsetSec, numFrame)
%LOCALPARSETIMEOFFSETSEC Parse optional frame time offsets.

if isempty(timeOffsetSecIn)
  timeOffsetSec = reshape(defaultTimeOffsetSec, 1, []);
else
  if ~isnumeric(timeOffsetSecIn) || numel(timeOffsetSecIn) ~= numFrame || ...
      any(~isfinite(timeOffsetSecIn))
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidTimeOffsetSec', ...
      'modelOpt.timeOffsetSec must contain one finite entry per frame.');
  end
  timeOffsetSec = reshape(timeOffsetSecIn, 1, []);
end
end

function [doaLb, doaUb] = localBuildDoaBounds(model)
%LOCALBUILDDOABOUNDS Build the DoA search box used by all MF refinements.
% The local DoA refinement must be enforced by the optimizer bounds rather
% than only stored as a nominal debug setting. Precomputing the DoA box here
% makes every path reuse the same local basin definition.

if isfield(model, 'cachedBoundsReady') && logical(model.cachedBoundsReady) && ...
    isfield(model, 'doaLb') && isfield(model, 'doaUb') && ...
    ~isempty(model.doaLb) && ~isempty(model.doaUb)
  doaLb = model.doaLb;
  doaUb = model.doaUb;
  return;
end

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

function [lb, ub] = localBuildBounds(model)
%LOCALBUILDBOUNDS Build box constraints for fmincon.

if isfield(model, 'lb') && isfield(model, 'ub') && ...
    ~isempty(model.lb) && ~isempty(model.ub)
  lb = model.lb;
  ub = model.ub;
  return;
end

[doaLb, doaUb] = localBuildDoaBounds(model);

lb = [doaLb; model.fdRange(1)];
ub = [doaUb; model.fdRange(2)];

if strcmp(model.fdRateMode, 'unknown')
  lb = [lb; model.fdRateRange(1)];
  ub = [ub; model.fdRateRange(2)];
end
if isfield(model, 'cachedBoundsReady') && ~logical(model.cachedBoundsReady)
  model.lb = lb;
  model.ub = ub;
end
end

function userStateRef = localBuildUserStateRef(sceneSeq, timeOffsetSec)
%LOCALBUILDUSERSTATEREF Build the motion reference used by dynamic lat/lon mode.

userStateRef = struct();
userStateRef.utcVec = sceneSeq.utcVec(:);
userStateRef.timeOffsetSec = reshape(timeOffsetSec, 1, []);

copyField = {'usrPosEci', 'usrVelEci', 'usrPosNominalEci', 'usrVelNominalEci'};
for iField = 1:numel(copyField)
  fieldName = copyField{iField};
  if isfield(sceneSeq, fieldName) && ~isempty(sceneSeq.(fieldName))
    userStateRef.(fieldName) = sceneSeq.(fieldName);
  end
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
