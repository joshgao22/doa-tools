function [obj, prof, aux] = evalDoaDopplerDynProfileLike(model, hyp)
%EVALDOADOPPLERDYNPROFILELIKE Evaluate the multi-frame profile likelihood.
%
% Evaluates the single-source profile likelihood used by the multi-frame
% pilot-based DoA-Doppler estimator. The formal mainline is the
% continuous-phase model around a reliable coarse DoA initializer, while the
% relaxed and independent modes are retained as controlled baselines. The
% function supports three phase modes:
%
%   1) phaseMode = 'continuous'
%      Shared satellite phase across frames. The per-frame amplitude stays
%      frame-dependent, while one common phase is profiled for each
%      satellite.
%
%   2) phaseMode = 'relaxed'
%      Uses the same absolute-time phase template as the continuous mode,
%      but relaxes each valid (satellite, frame) block to an independent
%      complex gain.
%
%   3) phaseMode = 'independent'
%      Uses a frame-local phase template and profiles one independent
%      complex gain for each valid (satellite, frame) block.
%
% Doppler-rate handling is unified through hyp.fdRateSat. Therefore:
%
%   - unknown-rate : pass the candidate fdRateSat built from the optimizer
%   - known-rate   : pass the fixed fdRateSat
%   - zero-rate    : pass fdRateSat = 0
%
%Syntax:
%  obj = evalDoaDopplerDynProfileLike(model, hyp)
%
%  [obj, prof] = evalDoaDopplerDynProfileLike(model, hyp)
%
%  [obj, prof, aux] = evalDoaDopplerDynProfileLike(model, hyp)
%
%Inputs:
%  model           - fixed model structure. Required fields:
%                    .array           : 1xNs cell of array structs
%                    .rxSig           : 1xNf cell, each frame is 1xNs cell
%                                       of MxT snapshot matrices
%                    .pilotPadCell    : 1xNf cell of row pilot vectors
%                    .timeSecCell     : 1xNf cell of local time vectors in s
%                    .timeOffsetSec   : 1xNf frame start times in s
%                    .frameMask       : Ns x Nf logical mask
%                    .wavelength      : wavelength in m
%                    .phaseMode       : 'continuous', 'relaxed',
%                                       or 'independent'
%                    .useLogObjective : logical
%
%                    Optional fields:
%                    .doaType         : 'angle' or 'latlon'
%                    .fdRateMode      : 'unknown', 'known', or 'zero'
%                    .phaseGridCount  : coarse common-phase grid size
%                                       default: 72
%                    .phaseRefine     : whether to refine the common phase
%                                       with fminbnd
%                                       default: true
%                    .continuousPhaseConsistencyWeight
%                                     : nonnegative weight for the shared
%                                       cross-frame continuous-phase bonus
%                                       used only in continuous mode
%                                       default: 0
%                    .continuousPhaseCollapsePenaltyWeight
%                                     : nonnegative weight for penalizing
%                                       support collapse onto only a few
%                                       surviving frames in continuous mode
%                                       default: 0
%                    .continuousPhaseNegativeProjectionPenaltyWeight
%                                     : nonnegative weight for penalizing
%                                       large negative shared-phase
%                                       projections in continuous mode
%                                       default: 0
%                    .continuousPhaseNonRefFitFloorWeight
%                                     : nonnegative weight for penalizing
%                                       multi-satellite solutions whose
%                                       non-reference explained-energy floor
%                                       collapses far below the reference
%                                       satellite level in continuous mode
%                                       default: 0
%                    .continuousPhaseNonRefSupportFloorWeight
%                                     : nonnegative weight for penalizing
%                                       multi-satellite solutions whose
%                                       non-reference shared-phase support
%                                       floor collapses far below the
%                                       reference satellite level in
%                                       continuous mode
%                                       default: 0
%                    .fdAliasStepHz    : optional Doppler ambiguity step in
%                                       Hz. When empty, a uniform frame-step
%                                       based value is inferred from
%                                       timeOffsetSec.
%                    .enableFdAliasUnwrap
%                                     : whether to unwrap per-satellite
%                                       Doppler candidates to the nearest
%                                       prior modulo fdAliasStepHz
%                                       default: false
%                    .fdSatPriorHz     : optional Ns x 1 Doppler prior used
%                                       only when enableFdAliasUnwrap is
%                                       true
%
%  hyp             - candidate hypothesis. Required fields:
%                    .localDoaArr     : 2xNs or 2xNsxNf local DoA array
%                    .fdSat           : Ns x 1 reference Doppler in Hz
%
%                    Optional fields:
%                    .fdRateSat       : Ns x 1 Doppler rate in Hz/s
%                                       default: zeros(Ns,1)
%
%Outputs:
%  obj             - negative log-likelihood surrogate after profiling the
%                    nuisance amplitude and phase terms. When
%                    model.useLogObjective is true, the objective is
%                    sum(count .* log(noiseVar)); otherwise it is the total
%                    residual energy.
%
%  prof            - profiled nuisance result structure with fields:
%                    .modelType
%                    .doaType
%                    .phaseMode
%                    .fdRateMode
%                    .frameMask
%                    .timeOffsetSec
%                    .pathGainEst     : Ns x Nf complex path gains
%                    .ampEst          : Ns x Nf nonnegative amplitudes
%                    .phaseSatEst     : Ns x 1 common phases for
%                                       continuous mode
%                    .framePhaseEst   : Ns x Nf frame phases for
%                                       relaxed or independent modes
%                    .noiseVarEst     : Ns x 1 residual noise powers
%                    .noiseVarGlobal  : scalar residual noise power
%                    .residualNorm    : total squared residual norm
%                    .countPerSat     : Ns x 1 valid sample counts
%                    .blockValue      : Ns x Nf profiled fit values
%                    .blockNorm2      : Ns x Nf block energies
%                    .countPerBlock   : Ns x Nf block sample counts
%                    .zMat            : Ns x Nf matched-filter outputs
%                    .etaMat          : Ns x Nf atom energies
%                    .fdSat           : Ns x 1 reference Doppler in Hz
%                    .fdLocal         : Ns x Nf frame-start Doppler in Hz
%                    .fdSatEval       : Ns x 1 Doppler values actually used
%                                       in the phase template after optional
%                                       ambiguity unwrapping
%                    .fdLocalEval     : Ns x Nf frame-start Doppler values
%                                       after optional ambiguity unwrapping
%                    .fdAliasStepHz   : ambiguity step used for optional
%                                       unwrapping
%                    .fdAliasIndex    : Ns x 1 integer ambiguity indices
%                    .fdAliasShiftHz  : Ns x 1 applied ambiguity shifts in Hz
%                    .fdRateSat       : Ns x 1 Doppler rates in Hz/s
%                    .effectiveFrameSupportSat
%                                     : Ns x 1 effective shared-phase frame
%                                       support count in continuous mode
%                    .effectiveFrameSupportRatioSat
%                                     : Ns x 1 normalized effective support
%                    .activeFrameSupportRatioSat
%                                     : Ns x 1 ratio of materially active
%                                       shared-phase frames
%                    .dominantFrameRatioSat
%                                     : Ns x 1 largest per-frame fit share
%                    .positiveAlignmentRatioSat
%                                     : Ns x 1 ratio of positive shared-phase
%                                       projection mass
%                    .negativeProjectionRatioSat
%                                     : Ns x 1 ratio of negative shared-phase
%                                       projection energy
%                    .collapsedFrameCountSat
%                                     : Ns x 1 number of frames clipped by
%                                       the shared-phase nonnegative profile
%                    .satFitRatio     : Ns x 1 per-satellite explained-energy
%                                       ratios
%                    .refFitRatio     : reference-satellite fit ratio
%                    .nonRefFitRatioFloor
%                                     : smallest non-reference fit ratio
%                    .refSupportRatio : reference-satellite shared-phase
%                                       support ratio
%                    .nonRefSupportRatioFloor
%                                     : smallest non-reference support ratio
%                    .nonRefFitFloorPenalty
%                                     : global CP penalty against ref-only
%                                       fit dominance
%                    .nonRefSupportFloorPenalty
%                                     : global CP penalty against non-ref
%                                       support collapse
%                    .localDoaArrUsed : 2xNsxNf steering angles used for
%                                       likelihood evaluation
%
%  aux             - auxiliary per-hypothesis diagnostics with fields:
%                    .modelType
%                    .doaType
%                    .phaseMode
%                    .fdRateMode
%                    .localDoaArrUsed
%                    .fdSat
%                    .fdRateSat
%                    .fdLocal
%                    .fdSatEval
%                    .fdLocalEval
%                    .fdAliasStepHz
%                    .fdAliasIndex
%                    .fdAliasShiftHz
%                    .countPerBlock
%                    .blockNorm2
%                    .effectiveFrameSupportSat
%                    .effectiveFrameSupportRatioSat
%                    .activeFrameSupportRatioSat
%                    .dominantFrameRatioSat
%                    .positiveAlignmentRatioSat
%                    .negativeProjectionRatioSat
%                    .collapsedFrameCountSat
%                    .satFitRatio
%                    .refFitRatio
%                    .nonRefFitRatioFloor
%                    .refSupportRatio
%                    .nonRefSupportRatioFloor
%                    .nonRefFitFloorPenalty
%                    .nonRefSupportFloorPenalty
%
%Notes:
%  - This function assumes one source.
%  - The returned pathGainEst always has the complex form b*exp(1j*phase),
%    even in the continuous mode where phase is shared across frames.
%
%See also:
%  estimatorDoaDopplerMlePilotMfOpt, steeringMatrix

arguments
  model (1,1) struct
  hyp (1,1) struct
end

model = localParseModel(model);
hyp = localParseHyp(model, hyp);
[blockStat, localDoaArrUsed, hypEval] = localBuildBlockStat(model, hyp);
prof = localProfileBlocks(model, hypEval, blockStat, localDoaArrUsed);

validSat = prof.countPerSat > 0 & isfinite(prof.noiseVarEst);
if ~any(validSat)
  obj = inf;
else
  if model.useLogObjective
    obj = sum(prof.countPerSat(validSat) .* log(max(prof.noiseVarEst(validSat), eps)));
  else
    obj = prof.residualNorm;
  end
  obj = obj + localGetStructField(prof, 'additionalObjectivePenalty', 0);
end

if nargout < 3
  return;
end

aux = struct();
aux.modelType = model.modelType;
aux.doaType = model.doaType;
aux.phaseMode = model.phaseMode;
aux.fdRateMode = model.fdRateMode;
aux.localDoaArrUsed = localDoaArrUsed;
aux.refSatIdxLocal = prof.refSatIdxLocal;
aux.fdSat = hypEval.fdSat;
aux.fdSatRaw = prof.fdSatRaw;
aux.fdRateSat = hypEval.fdRateSat;
aux.fdLocal = prof.fdLocal;
aux.fdLocalRaw = prof.fdLocalRaw;
aux.fdSatEval = prof.fdSatEval;
aux.fdLocalEval = prof.fdLocalEval;
aux.fdRefRaw = prof.fdRefRaw;
aux.fdRefEval = prof.fdRefEval;
aux.fdRateRef = prof.fdRateRef;
aux.deltaFdRefRaw = prof.deltaFdRefRaw;
aux.deltaFdRefEval = prof.deltaFdRefEval;
aux.deltaFdRateRaw = prof.deltaFdRateRaw;
aux.deltaFdRateEval = prof.deltaFdRateEval;
aux.residualSat = prof.residualSat;
aux.fitValueSat = prof.fitValueSat;
aux.objectiveSat = prof.objectiveSat;
aux.fdAliasStepHz = prof.fdAliasStepHz;
aux.fdAliasIndex = prof.fdAliasIndex;
aux.fdAliasShiftHz = prof.fdAliasShiftHz;
aux.countPerBlock = prof.countPerBlock;
aux.blockNorm2 = prof.blockNorm2;
aux.effectiveFrameSupportSat = prof.effectiveFrameSupportSat;
aux.effectiveFrameSupportRatioSat = prof.effectiveFrameSupportRatioSat;
aux.negativeProjectionRatioSat = prof.negativeProjectionRatioSat;
aux.activeFrameSupportRatioSat = prof.activeFrameSupportRatioSat;
aux.dominantFrameRatioSat = prof.dominantFrameRatioSat;
aux.positiveAlignmentRatioSat = prof.positiveAlignmentRatioSat;
aux.collapsedFrameCountSat = prof.collapsedFrameCountSat;
aux.supportBonusSat = prof.supportBonusSat;
aux.collapsePenaltySat = prof.collapsePenaltySat;
aux.negativePenaltySat = prof.negativePenaltySat;
aux.consistencyScoreSat = prof.consistencyScoreSat;
aux.consistencyNormSat = prof.consistencyNormSat;
aux.refConsistencyNorm = prof.refConsistencyNorm;
aux.nonRefConsistencyRatioFloor = prof.nonRefConsistencyRatioFloor;
end


function model = localParseModel(model)
%LOCALPARSEMODEL Validate the fixed model structure.

requiredField = {'array', 'rxSig', 'pilotPadCell', 'timeSecCell', ...
  'timeOffsetSec', 'frameMask', 'wavelength', 'phaseMode', 'useLogObjective'};
for iField = 1:numel(requiredField)
  if ~isfield(model, requiredField{iField})
    error('evalDoaDopplerDynProfileLike:MissingModelField', ...
      'model.%s is required.', requiredField{iField});
  end
end

if ~iscell(model.array) || isempty(model.array)
  error('evalDoaDopplerDynProfileLike:InvalidArray', ...
    'model.array must be a non-empty cell array.');
end
if ~iscell(model.rxSig) || isempty(model.rxSig)
  error('evalDoaDopplerDynProfileLike:InvalidRxSig', ...
    'model.rxSig must be a non-empty frame cell array.');
end
if ~iscell(model.pilotPadCell) || numel(model.pilotPadCell) ~= numel(model.rxSig)
  error('evalDoaDopplerDynProfileLike:InvalidPilotCell', ...
    'model.pilotPadCell must contain one entry per frame.');
end
if ~iscell(model.timeSecCell) || numel(model.timeSecCell) ~= numel(model.rxSig)
  error('evalDoaDopplerDynProfileLike:InvalidTimeSecCell', ...
    'model.timeSecCell must contain one entry per frame.');
end
if numel(model.timeOffsetSec) ~= numel(model.rxSig)
  error('evalDoaDopplerDynProfileLike:InvalidTimeOffset', ...
    'model.timeOffsetSec must contain one entry per frame.');
end
if ~isscalar(model.wavelength) || ~isfinite(model.wavelength) || model.wavelength <= 0
  error('evalDoaDopplerDynProfileLike:InvalidWavelength', ...
    'model.wavelength must be a positive finite scalar.');
end
if ~isscalar(model.useLogObjective) || ~islogical(model.useLogObjective)
  error('evalDoaDopplerDynProfileLike:InvalidUseLogObjective', ...
    'model.useLogObjective must be a logical scalar.');
end

model.numSat = numel(model.array);
model.numFrame = numel(model.rxSig);
model.frameMask = localNormalizeFrameMask(model.frameMask, model.numSat, model.numFrame);
model.timeOffsetSec = reshape(model.timeOffsetSec, 1, []);

if isstring(model.phaseMode)
  model.phaseMode = char(model.phaseMode);
end
model.phaseMode = lower(strtrim(model.phaseMode));
if ~ismember(model.phaseMode, {'continuous', 'relaxed', 'independent'})
  error('evalDoaDopplerDynProfileLike:InvalidPhaseMode', ...
    ['model.phaseMode must be ''continuous'', ''relaxed'' ', ...
     'or ''independent''.']);
end

if ~isfield(model, 'doaType') || isempty(model.doaType)
  model.doaType = '';
elseif isstring(model.doaType)
  model.doaType = char(model.doaType);
end
if ~isempty(model.doaType)
  model.doaType = lower(strtrim(model.doaType));
end

if ~isfield(model, 'fdRateMode') || isempty(model.fdRateMode)
  model.fdRateMode = 'unknown';
elseif isstring(model.fdRateMode)
  model.fdRateMode = char(model.fdRateMode);
end
model.fdRateMode = lower(strtrim(model.fdRateMode));
if ~ismember(model.fdRateMode, {'unknown', 'known', 'zero'})
  error('evalDoaDopplerDynProfileLike:InvalidFdRateMode', ...
    'model.fdRateMode must be ''unknown'', ''known'', or ''zero''.');
end

if ~isfield(model, 'phaseGridCount') || isempty(model.phaseGridCount)
  model.phaseGridCount = 72;
end
if ~isfield(model, 'phaseRefine') || isempty(model.phaseRefine)
  model.phaseRefine = true;
end
if ~isfield(model, 'continuousPhaseConsistencyWeight') || isempty(model.continuousPhaseConsistencyWeight)
  model.continuousPhaseConsistencyWeight = 0;
end
if ~isfield(model, 'continuousPhaseCollapsePenaltyWeight') || isempty(model.continuousPhaseCollapsePenaltyWeight)
  model.continuousPhaseCollapsePenaltyWeight = 0;
end
if ~isfield(model, 'continuousPhaseNegativeProjectionPenaltyWeight') || isempty(model.continuousPhaseNegativeProjectionPenaltyWeight)
  model.continuousPhaseNegativeProjectionPenaltyWeight = 0;
end
if ~isfield(model, 'continuousPhaseNonRefFitFloorWeight') || isempty(model.continuousPhaseNonRefFitFloorWeight)
  model.continuousPhaseNonRefFitFloorWeight = 0;
end
if ~isfield(model, 'continuousPhaseNonRefSupportFloorWeight') || isempty(model.continuousPhaseNonRefSupportFloorWeight)
  model.continuousPhaseNonRefSupportFloorWeight = 0;
end
if ~isfield(model, 'fdAliasStepHz') || isempty(model.fdAliasStepHz)
  model.fdAliasStepHz = localInferFdAliasStep(model.timeOffsetSec);
end
if ~isfield(model, 'enableFdAliasUnwrap') || isempty(model.enableFdAliasUnwrap)
  model.enableFdAliasUnwrap = false;
end
if ~isfield(model, 'fdSatPriorHz') || isempty(model.fdSatPriorHz)
  model.fdSatPriorHz = nan(model.numSat, 1);
end
if ~isscalar(model.phaseGridCount) || model.phaseGridCount < 8 || ...
    mod(model.phaseGridCount, 1) ~= 0
  error('evalDoaDopplerDynProfileLike:InvalidPhaseGridCount', ...
    'model.phaseGridCount must be an integer scalar no smaller than 8.');
end
if ~isscalar(model.phaseRefine) || ~islogical(model.phaseRefine)
  error('evalDoaDopplerDynProfileLike:InvalidPhaseRefine', ...
    'model.phaseRefine must be a logical scalar.');
end
if ~isscalar(model.continuousPhaseConsistencyWeight) || ...
    ~isfinite(model.continuousPhaseConsistencyWeight) || ...
    model.continuousPhaseConsistencyWeight < 0
  error('evalDoaDopplerDynProfileLike:InvalidContinuousPhaseConsistencyWeight', ...
    'model.continuousPhaseConsistencyWeight must be a nonnegative finite scalar.');
end
if ~isscalar(model.continuousPhaseCollapsePenaltyWeight) || ...
    ~isfinite(model.continuousPhaseCollapsePenaltyWeight) || ...
    model.continuousPhaseCollapsePenaltyWeight < 0
  error('evalDoaDopplerDynProfileLike:InvalidContinuousPhaseCollapsePenaltyWeight', ...
    'model.continuousPhaseCollapsePenaltyWeight must be a nonnegative finite scalar.');
end
if ~isscalar(model.continuousPhaseNegativeProjectionPenaltyWeight) || ...
    ~isfinite(model.continuousPhaseNegativeProjectionPenaltyWeight) || ...
    model.continuousPhaseNegativeProjectionPenaltyWeight < 0
  error('evalDoaDopplerDynProfileLike:InvalidContinuousPhaseNegativeProjectionPenaltyWeight', ...
    ['model.continuousPhaseNegativeProjectionPenaltyWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~isscalar(model.continuousPhaseNonRefFitFloorWeight) || ...
    ~isfinite(model.continuousPhaseNonRefFitFloorWeight) || ...
    model.continuousPhaseNonRefFitFloorWeight < 0
  error('evalDoaDopplerDynProfileLike:InvalidContinuousPhaseNonRefFitFloorWeight', ...
    ['model.continuousPhaseNonRefFitFloorWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~isscalar(model.continuousPhaseNonRefSupportFloorWeight) || ...
    ~isfinite(model.continuousPhaseNonRefSupportFloorWeight) || ...
    model.continuousPhaseNonRefSupportFloorWeight < 0
  error('evalDoaDopplerDynProfileLike:InvalidContinuousPhaseNonRefSupportFloorWeight', ...
    ['model.continuousPhaseNonRefSupportFloorWeight must be a ', ...
     'nonnegative finite scalar.']);
end
if ~(isscalar(model.enableFdAliasUnwrap) && islogical(model.enableFdAliasUnwrap))
  error('evalDoaDopplerDynProfileLike:InvalidEnableFdAliasUnwrap', ...
    'model.enableFdAliasUnwrap must be a logical scalar.');
end
if ~(isscalar(model.fdAliasStepHz) && (isnan(model.fdAliasStepHz) || ...
    (isfinite(model.fdAliasStepHz) && model.fdAliasStepHz > 0)))
  error('evalDoaDopplerDynProfileLike:InvalidFdAliasStepHz', ...
    'model.fdAliasStepHz must be empty, NaN, or a positive finite scalar.');
end
if ~isnumeric(model.fdSatPriorHz) || numel(model.fdSatPriorHz) ~= model.numSat
  error('evalDoaDopplerDynProfileLike:InvalidFdSatPriorHz', ...
    'model.fdSatPriorHz must contain one Doppler prior per satellite.');
end
model.fdSatPriorHz = reshape(model.fdSatPriorHz, [], 1);

model.modelType = localGetModelType(model.phaseMode, model.fdRateMode);

for iFrame = 1:model.numFrame
  currentFrame = model.rxSig{iFrame};
  if ~iscell(currentFrame)
    error('evalDoaDopplerDynProfileLike:InvalidRxFrame', ...
      'model.rxSig{%d} must be a 1xNs cell array.', iFrame);
  end
  if numel(currentFrame) ~= model.numSat
    error('evalDoaDopplerDynProfileLike:RxSatCountMismatch', ...
      'model.rxSig{%d} must contain one block per satellite.', iFrame);
  end

  pilotVec = model.pilotPadCell{iFrame};
  if ~isnumeric(pilotVec) || isempty(pilotVec)
    error('evalDoaDopplerDynProfileLike:InvalidPilotVector', ...
      'model.pilotPadCell{%d} must be a non-empty numeric vector.', iFrame);
  end
  pilotVec = pilotVec(:);
  model.pilotPadCell{iFrame} = pilotVec;

  timeVec = model.timeSecCell{iFrame};
  if ~isnumeric(timeVec) || isempty(timeVec) || numel(timeVec) ~= numel(pilotVec)
    error('evalDoaDopplerDynProfileLike:InvalidTimeVector', ...
      'model.timeSecCell{%d} must match the pilot length.', iFrame);
  end
  model.timeSecCell{iFrame} = timeVec(:);

  for iSat = 1:model.numSat
    yMat = currentFrame{iSat};
    if ~isnumeric(yMat) || ~ismatrix(yMat) || isempty(yMat)
      error('evalDoaDopplerDynProfileLike:InvalidRxBlock', ...
        'model.rxSig{%d}{%d} must be a non-empty numeric matrix.', iFrame, iSat);
    end
    if size(yMat, 2) ~= numel(pilotVec)
      error('evalDoaDopplerDynProfileLike:RxLengthMismatch', ...
        ['The sample length of model.rxSig{%d}{%d} must match the length ', ...
         'of model.pilotPadCell{%d}.'], iFrame, iSat, iFrame);
    end
  end
end
end


function frameMask = localNormalizeFrameMask(frameMask, numSat, numFrame)
%LOCALNORMALIZEFRAMEMASK Normalize the frame validity mask.

if isvector(frameMask) && numSat == 1 && numel(frameMask) == numFrame
  frameMask = reshape(logical(frameMask), 1, []);
elseif isvector(frameMask) && numFrame == 1 && numel(frameMask) == numSat
  frameMask = reshape(logical(frameMask), [], 1);
elseif islogical(frameMask) && isequal(size(frameMask), [numSat, numFrame])
  % Keep as is.
elseif isnumeric(frameMask) && isequal(size(frameMask), [numSat, numFrame])
  frameMask = logical(frameMask);
else
  error('evalDoaDopplerDynProfileLike:InvalidFrameMask', ...
    'model.frameMask must have size Ns x Nf.');
end

frameMask = logical(frameMask);
end


function hyp = localParseHyp(model, hyp)
%LOCALPARSEHYP Validate the candidate hypothesis.

requiredField = {'localDoaArr', 'fdSat'};
for iField = 1:numel(requiredField)
  if ~isfield(hyp, requiredField{iField})
    error('evalDoaDopplerDynProfileLike:MissingHypField', ...
      'hyp.%s is required.', requiredField{iField});
  end
end

if ~isfield(hyp, 'fdRateSat') || isempty(hyp.fdRateSat)
  hyp.fdRateSat = zeros(model.numSat, 1);
end

if ~isnumeric(hyp.fdSat) || numel(hyp.fdSat) ~= model.numSat
  error('evalDoaDopplerDynProfileLike:InvalidFdSat', ...
    'hyp.fdSat must contain one Doppler value per satellite.');
end
if ~isnumeric(hyp.fdRateSat) || numel(hyp.fdRateSat) ~= model.numSat
  error('evalDoaDopplerDynProfileLike:InvalidFdRateSat', ...
    'hyp.fdRateSat must contain one Doppler-rate value per satellite.');
end

hyp.fdSat = reshape(hyp.fdSat, [], 1);
hyp.fdRateSat = reshape(hyp.fdRateSat, [], 1);

if any(~isfinite(hyp.fdSat)) || any(~isfinite(hyp.fdRateSat))
  error('evalDoaDopplerDynProfileLike:InvalidDopplerValue', ...
    'hyp.fdSat and hyp.fdRateSat must contain finite values.');
end

localDoaArr = hyp.localDoaArr;
if ~isnumeric(localDoaArr) || size(localDoaArr, 1) ~= 2 || ...
    ~(isequal(size(localDoaArr), [2, model.numSat]) || ...
      isequal(size(localDoaArr), [2, model.numSat, model.numFrame]))
  error('evalDoaDopplerDynProfileLike:InvalidLocalDoaArr', ...
    ['hyp.localDoaArr must have size 2xNs or 2xNsxNf. ', ...
     'Use the 2xNsxNf form when steering is frame-dependent.']);
end
end


function [blockStat, localDoaArrUsed, hypEval] = localBuildBlockStat(model, hyp)
%LOCALBUILDBLOCKSTAT Build matched-filter statistics for all valid blocks.

nanComplex = complex(nan(model.numSat, model.numFrame), nan(model.numSat, model.numFrame));
blockStat = struct();
blockStat.zMat = nanComplex;
blockStat.etaMat = nan(model.numSat, model.numFrame);
blockStat.blockNorm2 = nan(model.numSat, model.numFrame);
blockStat.fdLocal = nan(model.numSat, model.numFrame);
blockStat.fdSatEval = nan(model.numSat, 1);
blockStat.fdLocalEval = nan(model.numSat, model.numFrame);
blockStat.fdAliasIndex = zeros(model.numSat, 1);
blockStat.fdAliasShiftHz = zeros(model.numSat, 1);
blockStat.fdAliasStepHz = model.fdAliasStepHz;
blockStat.countPerBlock = zeros(model.numSat, model.numFrame);

hypEval = hyp;
localDoaArrUsed = localNormalizeLocalDoaArr(hyp.localDoaArr, model.numSat, model.numFrame, ...
  'hyp.localDoaArr');
for iFrame = 1:model.numFrame
  pilotVec = model.pilotPadCell{iFrame};
  timeLoc = model.timeSecCell{iFrame};
  timeAbs = model.timeOffsetSec(iFrame) + timeLoc;

  for iSat = 1:model.numSat
    if ~model.frameMask(iSat, iFrame)
      continue;
    end

    yMat = model.rxSig{iFrame}{iSat};
    localDoa = localDoaArrUsed(:, iSat, iFrame);
    aVec = steeringMatrix(model.array{iSat}, model.wavelength, localDoa);
    aVec = aVec(:);

    fdSat = hypEval.fdSat(iSat);
    fdRateSat = hypEval.fdRateSat(iSat);
    [fdSatEval, fdAliasIndex, fdAliasShiftHz] = localResolveFdAlias(model, fdSat, iSat);
    fdLocal = fdSat + fdRateSat * model.timeOffsetSec(iFrame);
    fdLocalEval = fdSatEval + fdRateSat * model.timeOffsetSec(iFrame);

    pilotCol = reshape(pilotVec, [], 1);
    timeLocCol = reshape(timeLoc, [], 1);
    timeAbsCol = model.timeOffsetSec(iFrame) + timeLocCol;

    switch model.phaseMode
      case {'continuous', 'relaxed'}
        phase = 2 * pi * (fdSatEval * timeAbsCol + 0.5 * fdRateSat * (timeAbsCol .^ 2));

      case 'independent'
        phase = 2 * pi * (fdLocalEval * timeLocCol + 0.5 * fdRateSat * (timeLocCol .^ 2));

      otherwise
        error('evalDoaDopplerDynProfileLike:InvalidPhaseMode', ...
          'Unsupported phase mode: %s.', model.phaseMode);
    end

    qCol = pilotCol .* exp(1j * phase);
    zVal = aVec' * yMat * conj(qCol);
    etaVal = real(aVec' * aVec) * real(qCol' * qCol);

    if ~isfinite(etaVal) || etaVal <= eps
      error('evalDoaDopplerDynProfileLike:InvalidAtomEnergy', ...
        'The pilot atom energy is non-positive for satellite %d, frame %d.', ...
        iSat, iFrame);
    end

    blockStat.zMat(iSat, iFrame) = zVal;
    blockStat.etaMat(iSat, iFrame) = etaVal;
    blockStat.blockNorm2(iSat, iFrame) = real(sum(abs(yMat(:)) .^ 2));
    blockStat.fdLocal(iSat, iFrame) = fdLocal;
    blockStat.fdSatEval(iSat) = fdSatEval;
    blockStat.fdLocalEval(iSat, iFrame) = fdLocalEval;
    blockStat.fdAliasIndex(iSat) = fdAliasIndex;
    blockStat.fdAliasShiftHz(iSat) = fdAliasShiftHz;
    blockStat.countPerBlock(iSat, iFrame) = numel(yMat);
  end
end
end


function localDoaArr = localNormalizeLocalDoaArr(localDoaArrIn, numSat, numFrame, fieldName)
%LOCALNORMALIZELOCALDOAARR Normalize a local-DoA array to 2xNsxNf.

if isempty(localDoaArrIn)
  error('evalDoaDopplerDynProfileLike:MissingLocalDoaArr', ...
    '%s must not be empty.', fieldName);
end

if ndims(localDoaArrIn) == 2
  if ~isequal(size(localDoaArrIn), [2, numSat])
    error('evalDoaDopplerDynProfileLike:InvalidLocalDoaArr', ...
      '%s must have size 2xNs or 2xNsxNf.', fieldName);
  end
  localDoaArr = repmat(localDoaArrIn, 1, 1, numFrame);
  return;
end

if ~isequal(size(localDoaArrIn), [2, numSat, numFrame])
  error('evalDoaDopplerDynProfileLike:InvalidLocalDoaArr', ...
    '%s must have size 2xNs or 2xNsxNf.', fieldName);
end

localDoaArr = localDoaArrIn;
end


function satVec = localNormalizeSatVec(satVecIn, numSat, fieldName)
%LOCALNORMALIZESATVEC Normalize one per-satellite vector.

satVec = reshape(satVecIn, [], 1);
if numel(satVec) ~= numSat
  error('evalDoaDopplerDynProfileLike:InvalidSatVector', ...
    '%s must contain exactly one entry per satellite.', fieldName);
end
end


function prof = localProfileBlocks(model, hyp, blockStat, localDoaArrUsed)
%LOCALPROFILEBLOCKS Profile nuisance amplitude and phase parameters.

nanComplex = complex(nan(model.numSat, model.numFrame), nan(model.numSat, model.numFrame));
refSatIdxLocal = localResolveRefSatIdxForEval(model, hyp);
deltaFdRefRaw = localBuildRelativeSatValue(hyp.fdSat, refSatIdxLocal);
deltaFdRefEval = localBuildRelativeSatValue(blockStat.fdSatEval, refSatIdxLocal);
deltaFdRate = localBuildRelativeSatValue(hyp.fdRateSat, refSatIdxLocal);

prof = struct();
prof.modelType = model.modelType;
prof.doaType = model.doaType;
prof.phaseMode = model.phaseMode;
prof.fdRateMode = model.fdRateMode;
prof.frameMask = model.frameMask;
prof.timeOffsetSec = model.timeOffsetSec;
prof.refSatIdxLocal = refSatIdxLocal;
prof.pathGainEst = nanComplex;
prof.ampEst = nan(model.numSat, model.numFrame);
prof.phaseSatEst = nan(model.numSat, 1);
prof.framePhaseEst = nan(model.numSat, model.numFrame);
prof.noiseVarEst = nan(model.numSat, 1);
prof.noiseVarGlobal = nan;
prof.residualNorm = 0;
prof.countPerSat = zeros(model.numSat, 1);
prof.blockValue = nan(model.numSat, model.numFrame);
prof.blockNorm2 = blockStat.blockNorm2;
prof.countPerBlock = blockStat.countPerBlock;
prof.zMat = blockStat.zMat;
prof.etaMat = blockStat.etaMat;
prof.fdSat = hyp.fdSat;
prof.fdSatRaw = hyp.fdSat;
prof.fdLocal = blockStat.fdLocal;
prof.fdLocalRaw = blockStat.fdLocal;
prof.fdSatEval = blockStat.fdSatEval;
prof.fdLocalEval = blockStat.fdLocalEval;
prof.fdAliasStepHz = blockStat.fdAliasStepHz;
prof.fdAliasIndex = blockStat.fdAliasIndex;
prof.fdAliasShiftHz = blockStat.fdAliasShiftHz;
prof.fdRateSat = hyp.fdRateSat;
prof.fdRefRaw = hyp.fdSat(refSatIdxLocal);
prof.fdRefEval = blockStat.fdSatEval(refSatIdxLocal);
prof.fdRateRef = hyp.fdRateSat(refSatIdxLocal);
prof.deltaFdRefRaw = deltaFdRefRaw;
prof.deltaFdRefEval = deltaFdRefEval;
prof.deltaFdRateRaw = deltaFdRate;
prof.deltaFdRateEval = deltaFdRate;
prof.residualSat = nan(model.numSat, 1);
prof.fitValueSat = nan(model.numSat, 1);
prof.objectiveSat = nan(model.numSat, 1);
prof.effectiveFrameSupportSat = nan(model.numSat, 1);
prof.effectiveFrameSupportRatioSat = nan(model.numSat, 1);
prof.negativeProjectionRatioSat = nan(model.numSat, 1);
prof.activeFrameSupportRatioSat = nan(model.numSat, 1);
prof.dominantFrameRatioSat = nan(model.numSat, 1);
prof.positiveAlignmentRatioSat = nan(model.numSat, 1);
prof.collapsedFrameCountSat = nan(model.numSat, 1);
prof.supportBonusSat = nan(model.numSat, 1);
prof.collapsePenaltySat = nan(model.numSat, 1);
prof.negativePenaltySat = nan(model.numSat, 1);
prof.consistencyScoreSat = nan(model.numSat, 1);
prof.consistencyNormSat = nan(model.numSat, 1);
prof.satFitRatio = nan(model.numSat, 1);
prof.refFitRatio = nan;
prof.nonRefFitRatioFloor = nan;
prof.refSupportRatio = nan;
prof.nonRefSupportRatioFloor = nan;
prof.refConsistencyNorm = nan;
prof.nonRefConsistencyRatioFloor = nan;
prof.maxNonRefNegativeProjectionRatio = nan;
prof.nonRefFitFloorPenalty = 0;
prof.nonRefSupportFloorPenalty = 0;
prof.additionalObjectivePenalty = 0;
prof.localDoaArrUsed = localDoaArrUsed;

for iSat = 1:model.numSat
  validFrameIdx = find(model.frameMask(iSat, :));
  if isempty(validFrameIdx)
    continue;
  end

  zVec = blockStat.zMat(iSat, validFrameIdx).';
  etaVec = blockStat.etaMat(iSat, validFrameIdx).';
  blockNorm2 = blockStat.blockNorm2(iSat, validFrameIdx).';
  countVec = blockStat.countPerBlock(iSat, validFrameIdx).';

  switch model.phaseMode
    case 'continuous'
      [phaseSat, scoreVec, consistencyScore, consistencyDiag] = localEstimateCommonPhase(zVec, etaVec, model);
      projVec = real(exp(-1j * phaseSat) .* zVec);
      projPos = max(projVec, 0);
      ampVec = projPos ./ etaVec;
      gainVec = ampVec .* exp(1j * phaseSat);

      prof.phaseSatEst(iSat) = phaseSat;
      prof.framePhaseEst(iSat, validFrameIdx) = phaseSat;
      prof.effectiveFrameSupportSat(iSat) = consistencyDiag.effectiveFrameSupport;
      prof.effectiveFrameSupportRatioSat(iSat) = consistencyDiag.effectiveFrameSupportRatio;
      prof.negativeProjectionRatioSat(iSat) = consistencyDiag.negativeProjectionRatio;
      prof.activeFrameSupportRatioSat(iSat) = consistencyDiag.activeFrameSupportRatio;
      prof.dominantFrameRatioSat(iSat) = consistencyDiag.dominantFrameRatio;
      prof.positiveAlignmentRatioSat(iSat) = consistencyDiag.positiveAlignmentRatio;
      prof.collapsedFrameCountSat(iSat) = consistencyDiag.collapsedFrameCount;
      prof.supportBonusSat(iSat) = consistencyDiag.supportBonus;
      prof.collapsePenaltySat(iSat) = consistencyDiag.collapsePenalty;
      prof.negativePenaltySat(iSat) = consistencyDiag.negativePenalty;
      prof.consistencyScoreSat(iSat) = consistencyDiag.consistencyScore;

    case {'relaxed', 'independent'}
      consistencyScore = 0;
      consistencyDiag = localBuildEmptyConsistencyDiag(numel(zVec));
      prof.activeFrameSupportRatioSat(iSat) = consistencyDiag.activeFrameSupportRatio;
      prof.dominantFrameRatioSat(iSat) = consistencyDiag.dominantFrameRatio;
      prof.positiveAlignmentRatioSat(iSat) = consistencyDiag.positiveAlignmentRatio;
      prof.supportBonusSat(iSat) = consistencyDiag.supportBonus;
      prof.collapsePenaltySat(iSat) = consistencyDiag.collapsePenalty;
      prof.negativePenaltySat(iSat) = consistencyDiag.negativePenalty;
      prof.consistencyScoreSat(iSat) = consistencyDiag.consistencyScore;
      scoreVec = (abs(zVec) .^ 2) ./ etaVec;
      gainVec = zVec ./ etaVec;
      ampVec = abs(gainVec);
      phaseVec = angle(gainVec);

      prof.framePhaseEst(iSat, validFrameIdx) = phaseVec.';
      prof.phaseSatEst(iSat) = nan;

    otherwise
      error('evalDoaDopplerDynProfileLike:InvalidPhaseMode', ...
        'Unsupported phase mode: %s.', model.phaseMode);
  end

  rawFitValueSat = real(sum(scoreVec));
  fitValueSat = max(real(rawFitValueSat + consistencyScore), 0);
  residualSat = max(real(sum(blockNorm2) - fitValueSat), 0);
  countSat = sum(countVec);

  prof.pathGainEst(iSat, validFrameIdx) = gainVec.';
  prof.ampEst(iSat, validFrameIdx) = ampVec.';
  prof.blockValue(iSat, validFrameIdx) = scoreVec.';
  prof.countPerSat(iSat) = countSat;
  prof.fitValueSat(iSat) = fitValueSat;
  prof.residualSat(iSat) = residualSat;

  satEnergy = sum(blockNorm2);
  if isfinite(satEnergy) && satEnergy > 0
    prof.satFitRatio(iSat) = max(min(fitValueSat / satEnergy, 1), 0);
  end
  if isfinite(rawFitValueSat) && rawFitValueSat > 0 && isfinite(fitValueSat)
    prof.consistencyNormSat(iSat) = min(max(fitValueSat / rawFitValueSat, 0), 1);
  end

  if countSat > 0
    prof.noiseVarEst(iSat) = max(residualSat / countSat, eps);
    if model.useLogObjective
      prof.objectiveSat(iSat) = countSat * log(max(prof.noiseVarEst(iSat), eps));
    else
      prof.objectiveSat(iSat) = residualSat;
    end
    prof.residualNorm = prof.residualNorm + residualSat;
  end
end

totalCount = sum(prof.countPerSat);
if totalCount > 0
  prof.noiseVarGlobal = max(prof.residualNorm / totalCount, eps);
end

validSat = prof.countPerSat > 0 & isfinite(prof.objectiveSat);
if any(validSat)
  prof.objective = sum(prof.objectiveSat(validSat));
else
  prof.objective = inf;
end

[prof.nonRefFitFloorPenalty, prof.nonRefSupportFloorPenalty, prof.additionalObjectivePenalty, ...
  prof.refFitRatio, prof.nonRefFitRatioFloor, prof.refSupportRatio, ...
  prof.nonRefSupportRatioFloor, prof.refConsistencyNorm, ...
  prof.nonRefConsistencyRatioFloor, prof.maxNonRefNegativeProjectionRatio] = ...
  localBuildNonRefFloorPenalty(model, prof);
prof.objective = prof.objective + prof.additionalObjectivePenalty;
end


function [fitPenalty, supportPenalty, objectivePenalty, refFitRatio, ...
  nonRefFitFloor, refSupportRatio, nonRefSupportFloor, refConsistencyNorm, ...
  nonRefConsistencyRatioFloor, maxNonRefNegativeRatio] = ...
  localBuildNonRefFloorPenalty(model, prof)
%LOCALBUILDNONREFFLOORPENALTY Penalize CP solutions that collapse onto ref-only support.
% Wrong-tooth CP candidates can still look attractive when the reference
% satellite explains most energy while the non-reference satellites survive
% only through a much weaker support floor. Build one compact truth-free
% penalty from the final profiled per-satellite diagnostics so the main
% objective distinguishes those ref-dominant branches more clearly.

fitPenalty = 0;
supportPenalty = 0;
objectivePenalty = 0;
refFitRatio = nan;
nonRefFitFloor = nan;
refSupportRatio = nan;
nonRefSupportFloor = nan;
refConsistencyNorm = nan;
nonRefConsistencyRatioFloor = nan;
maxNonRefNegativeRatio = nan;

if ~strcmp(model.phaseMode, 'continuous') || model.numSat <= 1
  return;
end
refSatIdxLocal = localGetStructField(prof, 'refSatIdxLocal', NaN);
if ~(isscalar(refSatIdxLocal) && isfinite(refSatIdxLocal) && ...
    refSatIdxLocal >= 1 && refSatIdxLocal <= model.numSat)
  return;
end

fitRatioSat = reshape(localGetStructField(prof, 'satFitRatio', nan(model.numSat, 1)), [], 1);
supportRatioSat = reshape(localGetStructField(prof, 'effectiveFrameSupportRatioSat', nan(model.numSat, 1)), [], 1);
activeSupportRatioSat = reshape(localGetStructField(prof, 'activeFrameSupportRatioSat', nan(model.numSat, 1)), [], 1);
activeMask = isfinite(activeSupportRatioSat);
supportRatioSat(activeMask) = min(supportRatioSat(activeMask), activeSupportRatioSat(activeMask));
negativeRatioSat = reshape(localGetStructField(prof, 'negativeProjectionRatioSat', nan(model.numSat, 1)), [], 1);
consistencyNormSat = reshape(localGetStructField(prof, 'consistencyNormSat', nan(model.numSat, 1)), [], 1);
fitValueSat = reshape(localGetStructField(prof, 'fitValueSat', nan(model.numSat, 1)), [], 1);
blockNorm2Mat = localGetStructField(prof, 'blockNorm2', []);

if isempty(blockNorm2Mat)
  satEnergy = nan(model.numSat, 1);
else
  satEnergy = sum(blockNorm2Mat, 2, 'omitnan');
  satEnergy = reshape(satEnergy, [], 1);
end

refFitRatio = localSanitizeUnitInterval(fitRatioSat(refSatIdxLocal), nan);
refSupportRatio = localSanitizeUnitInterval(supportRatioSat(refSatIdxLocal), nan);
refConsistencyNorm = localSanitizeNonnegative(consistencyNormSat(refSatIdxLocal), nan);
nonRefMask = true(model.numSat, 1);
nonRefMask(refSatIdxLocal) = false;
nonRefFitFloor = localFiniteReducer(fitRatioSat(nonRefMask), @min, nan);
nonRefSupportFloor = localFiniteReducer(supportRatioSat(nonRefMask), @min, nan);
nonRefConsistencyFloor = localFiniteReducer(consistencyNormSat(nonRefMask), @min, nan);
if isfinite(refConsistencyNorm) && refConsistencyNorm > 0 && isfinite(nonRefConsistencyFloor)
  nonRefConsistencyRatioFloor = min(max(nonRefConsistencyFloor / refConsistencyNorm, 0), 1);
end
maxNonRefNegativeRatio = localFiniteReducer(negativeRatioSat(nonRefMask), @max, nan);

fitGap = 0;
if isfinite(refFitRatio) && isfinite(nonRefFitFloor)
  fitGap = max(refFitRatio - nonRefFitFloor, 0);
end
supportGap = 0;
if isfinite(refSupportRatio) && isfinite(nonRefSupportFloor)
  supportGap = max(refSupportRatio - nonRefSupportFloor, 0);
end

fitWeight = localGetStructField(model, 'continuousPhaseNonRefFitFloorWeight', 0);
supportWeight = localGetStructField(model, 'continuousPhaseNonRefSupportFloorWeight', 0);
if ~(isfinite(fitWeight) && fitWeight >= 0)
  fitWeight = 0;
end
if ~(isfinite(supportWeight) && supportWeight >= 0)
  supportWeight = 0;
end
if fitWeight == 0 && supportWeight == 0
  return;
end

energyScale = localFiniteReducer(satEnergy(isfinite(fitValueSat) & isfinite(satEnergy) & satEnergy > 0), @mean, nan);
countScale = localFiniteReducer(localGetStructField(prof, 'countPerSat', zeros(model.numSat, 1)), @sum, 0);
if model.useLogObjective
  penaltyScale = max(countScale, 1);
else
  penaltyScale = max(energyScale, 1);
end
if ~(isfinite(penaltyScale) && penaltyScale > 0)
  penaltyScale = 1;
end

fitPenalty = fitWeight * penaltyScale * fitGap .^ 2;
supportPenalty = supportWeight * penaltyScale * supportGap .^ 2;
objectivePenalty = fitPenalty + supportPenalty;
end


function value = localSanitizeUnitInterval(valueIn, defaultValue)
%LOCALSANITIZEUNITINTERVAL Clip one scalar metric to [0, 1].

value = defaultValue;
if isempty(valueIn)
  return;
end
valueIn = valueIn(1);
if ~isfinite(valueIn)
  return;
end
value = min(max(real(valueIn), 0), 1);
end


function value = localSanitizeNonnegative(valueIn, defaultValue)
%LOCALSANITIZENONNEGATIVE Clip one scalar metric to [0, inf).

value = defaultValue;
if isempty(valueIn)
  return;
end
valueIn = valueIn(1);
if ~isfinite(valueIn)
  return;
end
value = max(real(valueIn), 0);
end


function value = localFiniteReducer(valueVec, reducer, defaultValue)
%LOCALFINITEREDUCER Reduce one finite-valued vector with fallback.

value = defaultValue;
if isempty(valueVec)
  return;
end
valueVec = reshape(valueVec, [], 1);
valueVec = valueVec(isfinite(valueVec));
if isempty(valueVec)
  return;
end
value = reducer(valueVec);
end


function refSatIdxLocal = localResolveRefSatIdxForEval(model, hyp)
%LOCALRESOLVEREFSATIDXFOREVAL Resolve the local reference-satellite index.

refSatIdxLocal = [];

if isfield(model, 'refSatIdxLocal') && isscalar(model.refSatIdxLocal) && ...
    isfinite(model.refSatIdxLocal)
  refSatIdxLocal = round(model.refSatIdxLocal);
end

if isempty(refSatIdxLocal) && isfield(hyp, 'deltaFdRef') && ~isempty(hyp.deltaFdRef)
  deltaFdRef = reshape(hyp.deltaFdRef, [], 1);
  validMask = isfinite(deltaFdRef);
  if any(validMask)
    validIdx = find(validMask);
    [~, idxRel] = min(abs(deltaFdRef(validMask)));
    refSatIdxLocal = validIdx(idxRel);
  end
end

if isempty(refSatIdxLocal)
  refSatIdxLocal = 1;
end

refSatIdxLocal = max(1, min(model.numSat, refSatIdxLocal));
end


function deltaVal = localBuildRelativeSatValue(valueVec, refSatIdxLocal)
%LOCALBUILDRELATIVESATVALUE Build satellite values relative to the reference.

valueVec = reshape(valueVec, [], 1);
if isempty(valueVec) || refSatIdxLocal < 1 || refSatIdxLocal > numel(valueVec)
  deltaVal = nan(size(valueVec));
  return;
end

refVal = valueVec(refSatIdxLocal);
deltaVal = valueVec - refVal;
end


function [phaseSat, scoreVec, consistencyScore, consistencyDiag] = localEstimateCommonPhase(zVec, etaVec, model)
%LOCALESTIMATECOMMONPHASE Estimate the shared satellite phase.

if numel(zVec) == 1
  phaseSat = angle(zVec);
else
  coarsePhase = linspace(-pi, pi, model.phaseGridCount + 1);
  coarsePhase(end) = [];
  coarseScore = localEvalCommonPhaseScore(coarsePhase, zVec, etaVec, model);
  [~, bestIdx] = max(coarseScore);
  phaseSat = coarsePhase(bestIdx);

  if model.phaseRefine && localHasFminbnd()
    searchHalfWidth = 2 * pi / model.phaseGridCount;
    searchLeft = phaseSat - searchHalfWidth;
    searchRight = phaseSat + searchHalfWidth;
    refineFun = @(x) -localEvalCommonPhaseScore(x, zVec, etaVec, model);
    phaseRefined = fminbnd(refineFun, searchLeft, searchRight);

    if isfinite(phaseRefined)
      refinedScore = localEvalCommonPhaseScore(phaseRefined, zVec, etaVec, model);
      if refinedScore >= coarseScore(bestIdx)
        phaseSat = phaseRefined;
      end
    end
  end
end

phaseSat = localWrapToPi(phaseSat);
scoreVec = localEvalCommonPhaseTerm(phaseSat, zVec, etaVec);
scoreVec = scoreVec(:);
[consistencyScore, consistencyDiag] = localEvalCommonPhaseConsistency(phaseSat, zVec, etaVec, model);
end


function score = localEvalCommonPhaseScore(phaseVal, zVec, etaVec, model)
%LOCALEVALCOMMONPHASESCORE Evaluate the common-phase profile objective.

scoreTerm = localEvalCommonPhaseTerm(phaseVal, zVec, etaVec);
score = sum(scoreTerm, 2) + localEvalCommonPhaseConsistency(phaseVal, zVec, etaVec, model);
end


function scoreTerm = localEvalCommonPhaseTerm(phaseVal, zVec, etaVec)
%LOCALEVALCOMMONPHASETERM Evaluate the per-frame common-phase objective.

phaseVal = phaseVal(:);
rotProj = bsxfun(@times, exp(-1j * phaseVal), reshape(zVec, 1, []));
rotProj = real(rotProj);
rotProj(rotProj < 0) = 0;
scoreTerm = bsxfun(@rdivide, rotProj .^ 2, reshape(etaVec, 1, []));
end



function [consistencyScore, detail] = localEvalCommonPhaseConsistency(phaseVal, zVec, etaVec, model)
%LOCALEVALCOMMONPHASECONSISTENCY Add one stronger CP tying term.
% The formal per-frame nonnegative-amplitude profile can zero out frames
% that disagree with the shared phase, which leaves a strong 1/T_f comb in
% hard multi-frame cases. Keep the original cross-frame coherence bonus,
% but also penalize two wrong-tooth signatures:
%   1) support collapse onto only a few surviving frames; and
%   2) large negative shared-phase projections that are being clipped away.

if ~isfield(model, 'continuousPhaseConsistencyWeight') || ...
    isempty(model.continuousPhaseConsistencyWeight)
  bonusWeight = 0;
else
  bonusWeight = model.continuousPhaseConsistencyWeight;
end
if ~isfield(model, 'continuousPhaseCollapsePenaltyWeight') || ...
    isempty(model.continuousPhaseCollapsePenaltyWeight)
  collapseWeight = 0;
else
  collapseWeight = model.continuousPhaseCollapsePenaltyWeight;
end
if ~isfield(model, 'continuousPhaseNegativeProjectionPenaltyWeight') || ...
    isempty(model.continuousPhaseNegativeProjectionPenaltyWeight)
  negativeWeight = 0;
else
  negativeWeight = model.continuousPhaseNegativeProjectionPenaltyWeight;
end

phaseVal = phaseVal(:);
numFrame = numel(etaVec);
if ~(isfinite(bonusWeight) && bonusWeight >= 0 && ...
    isfinite(collapseWeight) && collapseWeight >= 0 && ...
    isfinite(negativeWeight) && negativeWeight >= 0) || ...
    (bonusWeight == 0 && collapseWeight == 0 && negativeWeight == 0)
  consistencyScore = zeros(numel(phaseVal), 1);
  if numel(phaseVal) == 1
    detail = localBuildEmptyConsistencyDiag(numFrame);
  else
    detail = struct();
  end
  return;
end

phaseVal = phaseVal(:);
rotVec = bsxfun(@times, exp(-1j * phaseVal), reshape(zVec, 1, []));
projRaw = real(rotVec);
projPos = max(projRaw, 0);
projNeg = max(-projRaw, 0);
etaRow = reshape(etaVec, 1, []);
etaSqrt = sqrt(max(etaRow, eps));
scoreRaw = projPos .^ 2 ./ max(etaRow, eps);
fitBase = sum(scoreRaw, 2);

sharedDenom = max(sum(etaRow, 2), eps);
coherentBase = abs(sum(rotVec, 2)) .^ 2 ./ sharedDenom;
penaltyBase = max(fitBase, coherentBase);

projNorm = projPos ./ etaSqrt;
effSupport = (sum(projNorm, 2) .^ 2) ./ max(sum(projNorm .^ 2, 2), eps);
effSupport = min(max(effSupport, 0), numFrame);
effSupportRatio = effSupport / max(numFrame, 1);

scorePeak = max(scoreRaw, [], 2);
activeMask = scoreRaw > bsxfun(@times, 0.05, scorePeak);
activeFrameCount = sum(activeMask, 2);
activeFrameSupportRatio = activeFrameCount / max(numFrame, 1);

dominantFrameRatio = max(scoreRaw, [], 2) ./ max(fitBase, eps);
dominantFrameRatio = min(max(dominantFrameRatio, 0), 1);

absProjSum = max(sum(abs(projRaw), 2), eps);
positiveAlignmentRatio = sum(projPos, 2) ./ absProjSum;
positiveAlignmentRatio = min(max(positiveAlignmentRatio, 0), 1);

negDenom = max(sum(abs(rotVec) .^ 2, 2), eps);
negativeRatio = sum(projNeg .^ 2, 2) ./ negDenom;
negativeRatio = min(max(negativeRatio, 0), 1);

supportRatio = min(effSupportRatio, activeFrameSupportRatio);
collapseMass = max(1 - supportRatio, 0) + max(dominantFrameRatio - 1 / max(numFrame, 1), 0);
negativeMass = 0.5 * negativeRatio + 0.5 * (1 - positiveAlignmentRatio);

supportBonus = bonusWeight * supportRatio .* coherentBase;
collapsePenalty = collapseWeight * collapseMass .* penaltyBase;
negativePenalty = negativeWeight * negativeMass .* penaltyBase;
consistencyScore = supportBonus - collapsePenalty - negativePenalty;

if numel(phaseVal) == 1
  detail = struct();
  detail.effectiveFrameSupport = effSupport;
  detail.effectiveFrameSupportRatio = effSupportRatio;
  detail.activeFrameSupportRatio = activeFrameSupportRatio;
  detail.dominantFrameRatio = dominantFrameRatio;
  detail.positiveAlignmentRatio = positiveAlignmentRatio;
  detail.negativeProjectionRatio = negativeRatio;
  detail.collapsedFrameCount = numFrame - activeFrameCount;
  detail.coherentBase = coherentBase;
  detail.fitBase = fitBase;
  detail.penaltyBase = penaltyBase;
  detail.supportBonus = supportBonus;
  detail.collapsePenalty = collapsePenalty;
  detail.negativePenalty = negativePenalty;
  detail.consistencyScore = consistencyScore;
else
  detail = struct();
end
end


function detail = localBuildEmptyConsistencyDiag(numFrame)
%LOCALBUILDEMPTYCONSISTENCYDIAG Build the default CP consistency diagnostics.

detail = struct();
detail.effectiveFrameSupport = nan;
detail.effectiveFrameSupportRatio = nan;
detail.activeFrameSupportRatio = nan;
detail.dominantFrameRatio = nan;
detail.positiveAlignmentRatio = nan;
detail.negativeProjectionRatio = nan;
detail.collapsedFrameCount = nan;
detail.coherentBase = 0;
detail.fitBase = 0;
detail.penaltyBase = 0;
detail.supportBonus = 0;
detail.collapsePenalty = 0;
detail.negativePenalty = 0;
detail.consistencyScore = 0;
if nargin >= 1 && isfinite(numFrame) && numFrame >= 0
  detail.collapsedFrameCount = numFrame;
end
end


function modelType = localGetModelType(phaseMode, fdRateMode)
%LOCALGETMODELTYPE Build a compact model tag.

phaseTagMap = struct();
phaseTagMap.continuous = 'Cp';
phaseTagMap.relaxed = 'Relaxed';
phaseTagMap.independent = 'Ip';

fdRateTagMap = struct();
fdRateTagMap.unknown = 'UnknownRate';
fdRateTagMap.known = 'KnownRate';
fdRateTagMap.zero = 'ZeroRate';

modelType = ['mf', phaseTagMap.(phaseMode), fdRateTagMap.(fdRateMode)];
end


function tf = localHasFminbnd()
%LOCALHASFMINBND Check whether fminbnd is available.

tf = exist('fminbnd', 'file') == 2;
end


function fdAliasStepHz = localInferFdAliasStep(timeOffsetSec)
%LOCALINFERFDALIASSTEP Infer the frame-to-frame Doppler ambiguity step.

timeOffsetSec = reshape(timeOffsetSec, 1, []);
if numel(timeOffsetSec) < 2
  fdAliasStepHz = nan;
  return;
end

dt = diff(timeOffsetSec);
dt = dt(isfinite(dt) & abs(dt) > eps);
if isempty(dt)
  fdAliasStepHz = nan;
  return;
end

dtMed = median(dt);
if ~isfinite(dtMed) || abs(dtMed) <= eps
  fdAliasStepHz = nan;
  return;
end

if max(abs(dt - dtMed)) > 1e-6 * max(abs(dtMed), eps)
  fdAliasStepHz = nan;
  return;
end

fdAliasStepHz = 1 / abs(dtMed);
end


function [fdSatEval, fdAliasIndex, fdAliasShiftHz] = localResolveFdAlias(model, fdSat, satIdx)
%LOCALRESOLVEFDALIAS Optionally unwrap one satellite Doppler by a prior.

fdSatEval = fdSat;
fdAliasIndex = 0;
fdAliasShiftHz = 0;

if ~model.enableFdAliasUnwrap
  return;
end
if ~(isfinite(model.fdAliasStepHz) && model.fdAliasStepHz > 0)
  return;
end
if satIdx < 1 || satIdx > numel(model.fdSatPriorHz)
  return;
end

fdPrior = model.fdSatPriorHz(satIdx);
if ~isfinite(fdPrior)
  return;
end

fdAliasIndex = round((fdPrior - fdSat) / model.fdAliasStepHz);
fdAliasShiftHz = fdAliasIndex * model.fdAliasStepHz;
fdSatEval = fdSat + fdAliasShiftHz;
end


function fieldValue = localGetStructField(dataStruct, fieldName, defaultValue)
%LOCALGETSTRUCTFIELD Read one struct field with a default fallback.

fieldValue = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    fieldValue = rawValue;
  end
end
end


function angleOut = localWrapToPi(angleIn)
%LOCALWRAPTOPI Wrap angles to [-pi, pi).

angleOut = mod(angleIn + pi, 2 * pi) - pi;
end
