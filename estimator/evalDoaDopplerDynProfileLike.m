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
      [phaseSat, scoreVec] = localEstimateCommonPhase(zVec, etaVec, model);
      projVec = real(exp(-1j * phaseSat) .* zVec);
      projPos = max(projVec, 0);
      ampVec = projPos ./ etaVec;
      gainVec = ampVec .* exp(1j * phaseSat);

      prof.phaseSatEst(iSat) = phaseSat;
      prof.framePhaseEst(iSat, validFrameIdx) = phaseSat;

    case {'relaxed', 'independent'}
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

  fitValueSat = real(sum(scoreVec));
  residualSat = max(real(sum(blockNorm2 - scoreVec)), 0);
  countSat = sum(countVec);

  prof.pathGainEst(iSat, validFrameIdx) = gainVec.';
  prof.ampEst(iSat, validFrameIdx) = ampVec.';
  prof.blockValue(iSat, validFrameIdx) = scoreVec.';
  prof.countPerSat(iSat) = countSat;
  prof.fitValueSat(iSat) = fitValueSat;
  prof.residualSat(iSat) = residualSat;

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


function [phaseSat, scoreVec] = localEstimateCommonPhase(zVec, etaVec, model)
%LOCALESTIMATECOMMONPHASE Estimate the shared satellite phase.

if numel(zVec) == 1
  phaseSat = angle(zVec);
else
  coarsePhase = linspace(-pi, pi, model.phaseGridCount + 1);
  coarsePhase(end) = [];
  coarseScore = localEvalCommonPhaseScore(coarsePhase, zVec, etaVec);
  [~, bestIdx] = max(coarseScore);
  phaseSat = coarsePhase(bestIdx);

  if model.phaseRefine && localHasFminbnd()
    searchHalfWidth = 2 * pi / model.phaseGridCount;
    searchLeft = phaseSat - searchHalfWidth;
    searchRight = phaseSat + searchHalfWidth;
    refineFun = @(x) -localEvalCommonPhaseScore(x, zVec, etaVec);
    phaseRefined = fminbnd(refineFun, searchLeft, searchRight);

    if isfinite(phaseRefined)
      refinedScore = localEvalCommonPhaseScore(phaseRefined, zVec, etaVec);
      if refinedScore >= coarseScore(bestIdx)
        phaseSat = phaseRefined;
      end
    end
  end
end

phaseSat = localWrapToPi(phaseSat);
scoreVec = localEvalCommonPhaseTerm(phaseSat, zVec, etaVec);
scoreVec = scoreVec(:);
end


function score = localEvalCommonPhaseScore(phaseVal, zVec, etaVec)
%LOCALEVALCOMMONPHASESCORE Evaluate the common-phase profile objective.

scoreTerm = localEvalCommonPhaseTerm(phaseVal, zVec, etaVec);
score = sum(scoreTerm, 2);
end


function scoreTerm = localEvalCommonPhaseTerm(phaseVal, zVec, etaVec)
%LOCALEVALCOMMONPHASETERM Evaluate the per-frame common-phase objective.

phaseVal = phaseVal(:);
rotProj = bsxfun(@times, exp(-1j * phaseVal), reshape(zVec, 1, []));
rotProj = real(rotProj);
rotProj(rotProj < 0) = 0;
scoreTerm = bsxfun(@rdivide, rotProj .^ 2, reshape(etaVec, 1, []));
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


function angleOut = localWrapToPi(angleIn)
%LOCALWRAPTOPI Wrap angles to [-pi, pi).

angleOut = mod(angleIn + pi, 2 * pi) - pi;
end
