function debugAux = buildDoaDopplerMfDebug(model, initDiag, initEvalDiag, finalEvalDiag, debugTrace)
%BUILDDOADOPPLERMFDEBUG Build compact multi-frame debug summaries.
% This helper centralizes the compact debug payload used by the
% multi-frame estimator, including init/final summaries, optional probe
% diagnostics, truth-state checks, and stored evaluation trace.

arguments
  model (1,1) struct
  initDiag (1,1) struct = struct()
  initEvalDiag (1,1) struct = struct()
  finalEvalDiag (1,1) struct = struct()
  debugTrace (1,1) struct = struct()
end

debugAux = localBuildDebugAux(model, initDiag, initEvalDiag, finalEvalDiag, debugTrace);
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

function refSatIdx = localResolveRefSatIdx(model)
%LOCALRESOLVEREFSATIDX Resolve the local reference-satellite index.

refSatIdx = [];

if isfield(model, 'refSatIdxLocal') && ~isempty(model.refSatIdxLocal) && ...
    isscalar(model.refSatIdxLocal) && isfinite(model.refSatIdxLocal)
  refSatIdx = model.refSatIdxLocal;
end

if isempty(refSatIdx) && isfield(model, 'sceneCell') && ~isempty(model.sceneCell) && ...
    isfield(model, 'refFrameIdx') && isscalar(model.refFrameIdx) && ...
    model.refFrameIdx >= 1 && model.refFrameIdx <= numel(model.sceneCell)
  refScene = model.sceneCell{model.refFrameIdx};
  if isstruct(refScene) && isfield(refScene, 'ref') && isstruct(refScene.ref) && ...
      isfield(refScene.ref, 'weight') && ~isempty(refScene.ref.weight)
    refWeight = refScene.ref.weight(:);
    refSatIdx = find(refWeight > 0.5, 1, 'first');
  end
end

if isempty(refSatIdx)
  refSatIdx = localGetStructField(model.debugTruth, 'refSatIdxLocal', []);
end

if isempty(refSatIdx) || ~isscalar(refSatIdx) || ~isfinite(refSatIdx)
  refSatIdx = 1;
end

refSatIdx = max(1, min(model.numSat, round(refSatIdx)));
end

function [fdRefPerSat, fdRatePerSat] = localFitSatFdLine(timeOffsetSec, fdSeq)
%LOCALFITSATFDLINE Fit one affine Doppler trend per satellite.

numSat = size(fdSeq, 1);
fdRefPerSat = zeros(numSat, 1);
fdRatePerSat = zeros(numSat, 1);

timeVec = reshape(timeOffsetSec, [], 1);
if numel(timeVec) < 2 || max(timeVec) - min(timeVec) <= 0
  fdRefPerSat = fdSeq(:, 1);
  return;
end

designMat = [ones(numel(timeVec), 1), timeVec];

for iSat = 1:numSat
  fdVec = reshape(fdSeq(iSat, :), [], 1);
  validMask = isfinite(fdVec) & isfinite(timeVec);

  if nnz(validMask) == 0
    continue;
  elseif nnz(validMask) == 1
    fdRefPerSat(iSat) = fdVec(validMask);
    continue;
  end

  coef = designMat(validMask, :) \ fdVec(validMask);
  fdRefPerSat(iSat) = coef(1);
  fdRatePerSat(iSat) = coef(2);
end
end

function evalDiag = localBuildEvalDiag(model, optVar, obj, noiseVar, evalAux)
%LOCALBUILDEVALDIAG Pack one compact evaluation summary for debugging.

evalDiag = struct();
evalDiag.optVar = optVar(:);
evalDiag.obj = obj;
evalDiag.doaType = model.doaType;
evalDiag.phaseMode = model.phaseMode;
evalDiag.fdRateMode = model.fdRateMode;
evalDiag.refSatIdxLocal = localGetStructField(evalAux, 'refSatIdxLocal', NaN);
evalDiag.doaParam = evalAux.doaParam;
evalDiag.eciAngle = evalAux.eciAngle;
evalDiag.latlon = evalAux.latlon;
evalDiag.fdRef = evalAux.fdRef;
evalDiag.fdRate = evalAux.fdRate;
evalDiag.deltaFdRef = evalAux.deltaFdRef(:);
evalDiag.deltaFdRefRaw = localGetStructField(evalAux, 'deltaFdRefRaw', evalAux.deltaFdRef(:));
evalDiag.deltaFdRefEval = localGetStructField(evalAux, 'deltaFdRefEval', evalAux.deltaFdRef(:));
evalDiag.deltaFdRate = evalAux.deltaFdRate(:);
evalDiag.deltaFdRateRaw = localGetStructField(evalAux, 'deltaFdRateRaw', evalAux.deltaFdRate(:));
evalDiag.deltaFdRateEval = localGetStructField(evalAux, 'deltaFdRateEval', evalAux.deltaFdRate(:));
evalDiag.fdSat = evalAux.fdSat;
evalDiag.fdSatRaw = localGetStructField(evalAux, 'fdSatRaw', evalAux.fdSat);
evalDiag.fdSatEval = evalAux.fdSatEval;
evalDiag.fdRateSat = evalAux.fdRateSat;
evalDiag.fdLocal = evalAux.fdLocal;
evalDiag.fdLocalRaw = localGetStructField(evalAux, 'fdLocalRaw', evalAux.fdLocal);
evalDiag.fdLocalEval = evalAux.fdLocalEval;
evalDiag.localDoaArrUsed = localGetStructField(evalAux, 'localDoaArrUsed', []);
evalDiag.fdAliasStepHz = evalAux.fdAliasStepHz;
evalDiag.fdAliasIndex = evalAux.fdAliasIndex;
evalDiag.fdAliasShiftHz = evalAux.fdAliasShiftHz;
evalDiag.noiseVar = noiseVar(:);
evalDiag.countPerSat = evalAux.countPerSat(:);
evalDiag.residualSat = localGetStructField(evalAux, 'residualSat', noiseVar(:) .* evalAux.countPerSat(:));
evalDiag.fitValueSat = localGetStructField(evalAux, 'fitValueSat', nan(size(evalDiag.residualSat)));
evalDiag.objectiveSat = localGetStructField(evalAux, 'objectiveSat', nan(size(evalDiag.residualSat)));
evalDiag.residualNorm = evalAux.residualNorm;
evalDiag.noiseVarGlobal = evalAux.noiseVarGlobal;
evalDiag.coherenceSat = localCalcSatCoherence(evalAux.zMat);
evalDiag.zMat = localGetStructField(evalAux, 'zMat', []);
evalDiag.etaMat = localGetStructField(evalAux, 'etaMat', []);
evalDiag.blockValueMat = localGetStructField(evalAux, 'blockValue', []);
evalDiag.blockNorm2Mat = localGetStructField(evalAux, 'blockNorm2', []);
evalDiag.blockResidualMat = evalDiag.blockNorm2Mat - evalDiag.blockValueMat;
evalDiag.countPerBlock = localGetStructField(evalAux, 'countPerBlock', []);
evalDiag.blockAbsZMat = abs(evalDiag.zMat);
evalDiag.blockPhaseResidMat = localCalcBlockPhaseResid(evalDiag.zMat, evalAux.framePhase);
end

function coherenceSat = localCalcSatCoherence(zMat)
%LOCALCALCSATCOHERENCE Compute |sum(z)| / sum(|z|) per satellite.

if isempty(zMat)
  coherenceSat = [];
  return;
end

magSum = sum(abs(zMat), 2);
coherenceSat = nan(size(magSum));
validMask = magSum > 0;
coherenceSat(validMask) = abs(sum(zMat(validMask, :), 2)) ./ magSum(validMask);
coherenceSat = coherenceSat(:);
end

function blockPhaseResidMat = localCalcBlockPhaseResid(zMat, framePhaseMat)
%LOCALCALCBLOCKPHASERESID Compute angle(z)-phase for each valid block.

blockPhaseResidMat = nan(size(zMat));
if isempty(zMat) || isempty(framePhaseMat)
  return;
end

validMask = isfinite(zMat) & isfinite(framePhaseMat);
blockPhaseResidMat(validMask) = angle(zMat(validMask)) - framePhaseMat(validMask);
blockPhaseResidMat(validMask) = mod(blockPhaseResidMat(validMask) + pi, 2 * pi) - pi;
end

function probeDiag = localBuildProbeDiag(model, initEvalDiag)
%LOCALBUILDPROBEDIAG Evaluate truth and mixed debug probe points.

probeDiag = struct();
probeDiag.isEnabled = false;

debugTruth = model.debugTruth;
if ~isstruct(debugTruth) || isempty(fieldnames(debugTruth))
  return;
end

truthDoa = localGetStructField(debugTruth, 'doaParam', []);
truthFdRef = localGetStructField(debugTruth, 'fdRef', NaN);
truthFdRate = localGetStructField(debugTruth, 'fdRate', NaN);

if isempty(truthDoa) || ~isfinite(truthFdRef)
  return;
end

switch model.fdRateMode
  case 'known'
    truthFdRate = model.fdRateKnown;
  case 'zero'
    truthFdRate = 0;
  otherwise
    if ~isfinite(truthFdRate)
      return;
    end
end

probeDiag.isEnabled = true;
probeDiag.truth = localEvalProbePoint(model, truthDoa, truthFdRef, truthFdRate, false);
probeDiag.mixDoaTrue = localEvalProbePoint(model, truthDoa, initEvalDiag.fdRef, initEvalDiag.fdRate, false);
probeDiag.mixFdTrue = localEvalProbePoint(model, initEvalDiag.doaParam, truthFdRef, truthFdRate, false);
end

function probeEval = localEvalProbePoint(model, doaParam, fdRef, fdRate, clipToBounds)
%LOCALEVALPROBEPOINT Evaluate one compact debug probe point.

if nargin < 5
  clipToBounds = true;
end

optVar = localBuildProbeOptVar(model, doaParam, fdRef, fdRate, clipToBounds);
[obj, ~, noiseVar, evalAux] = localNllPilotDynCon(model, optVar);
probeEval = localBuildEvalDiag(model, optVar, obj, noiseVar, evalAux);
end

function optVar = localBuildProbeOptVar(model, doaParam, fdRef, fdRate, clipToBounds)
%LOCALBUILDPROBEOPTVAR Build one optimization vector for probing.
% Truth and mixed-point probes are diagnostic only. They should reflect the
% actual objective at the requested point rather than the current local DoA
% box, otherwise CP-U continuation runs can report a fake "truth" point
% that was silently clipped back into a bad local basin.

if nargin < 5
  clipToBounds = true;
end

switch model.fdRateMode
  case 'unknown'
    optVar = [doaParam(:); fdRef; fdRate];
  case {'known', 'zero'}
    optVar = [doaParam(:); fdRef];
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end

if clipToBounds
  [lb, ub] = localBuildBounds(model);
  optVar = min(max(optVar(:), lb), ub);
else
  optVar = optVar(:);
end
if strcmp(model.doaType, 'angle')
  optVar(1) = mod(optVar(1), 2*pi);
end
end

function debugAux = localBuildDebugAux(model, initDiag, initEvalDiag, finalEvalDiag, debugTrace)
%LOCALBUILDDEBUGAUX Build compact debug outputs for diagnosis.

debugAux = struct();
debugAux.isEnabled = model.debugEnable;
debugAux.init = initDiag;
debugAux.initEval = initEvalDiag;
debugAux.finalEval = finalEvalDiag;
debugAux.probe = localBuildProbeDiag(model, initEvalDiag);
debugAux.truthState = localBuildTruthStateDiag(model);
debugAux.satOrder = localBuildSatOrderDiag(model);
debugAux.objImprove = initEvalDiag.obj - finalEvalDiag.obj;
debugAux.residualImprove = initEvalDiag.residualNorm - finalEvalDiag.residualNorm;
debugAux.trace = debugTrace;
debugAux.fdRange = model.fdRange(:);
debugAux.fdRateRange = model.fdRateRange(:);
debugAux.frameMask = model.frameMask;
debugAux.timeOffsetSec = model.timeOffsetSec(:);
debugAux.steeringMode = model.steeringMode;
debugAux.phaseMode = model.phaseMode;
debugAux.fdRateMode = model.fdRateMode;
debugAux.enableFdAliasUnwrap = model.enableFdAliasUnwrap;
debugAux.fdAliasStepHz = model.fdAliasStepHz;
debugAux.fdSatPriorHz = model.fdSatPriorHz;
end

function truthStateDiag = localBuildTruthStateDiag(model)
%LOCALBUILDTRUTHSTATEDIAG Evaluate truth-point mapping consistency.

truthStateDiag = struct();
truthStateDiag.isEnabled = false;
truthStateDiag.summaryTable = table();
truthStateDiag.maxLocalDirErrDeg = NaN;
truthStateDiag.maxFdLocalErrHz = NaN;
truthStateDiag.refFrameIdx = model.refFrameIdx;

if ~model.debugEnable
  return;
end

debugTruth = model.debugTruth;
if ~isstruct(debugTruth) || isempty(fieldnames(debugTruth))
  return;
end

truthDoa = localGetStructField(debugTruth, 'doaParam', []);
truthFdRef = localGetStructField(debugTruth, 'fdRef', NaN);
truthFdRate = localGetStructField(debugTruth, 'fdRate', NaN);
truthLocalDoa = localGetStructField(debugTruth, 'localDoaTrue', []);
truthFdLocal = localGetStructField(debugTruth, 'fdLocalTrue', []);
satIdxGlobal = reshape(localGetStructField(debugTruth, 'satIdxGlobal', []), [], 1);

if isempty(truthDoa) || ~isfinite(truthFdRef)
  return;
end

switch model.fdRateMode
  case 'known'
    truthFdRate = model.fdRateKnown;
  case 'zero'
    truthFdRate = 0;
  otherwise
    if ~isfinite(truthFdRate)
      return;
    end
end

modelTruth = localBuildTruthStateModel(model);
truthEval = localEvalProbePoint(modelTruth, truthDoa, truthFdRef, truthFdRate, false);
numSat = max([size(truthEval.localDoaArrUsed, 2), size(truthLocalDoa, 2), numel(satIdxGlobal), model.numSat]);
numFrame = max([size(truthEval.localDoaArrUsed, 3), size(truthLocalDoa, 3), size(truthFdLocal, 2), model.numFrame]);
refFrameIdx = min(max(model.refFrameIdx, 1), max(numFrame, 1));
rowCell = {};
maxLocalDirErrDeg = NaN;
maxFdLocalErrHz = NaN;

for iSat = 1:numSat
  [refDirErrDeg, rmsDirErrDeg, maxDirErrDeg, refAzErrDeg, refElErrDeg] = ...
    localCalcTruthStateLocalDoaErr(truthEval.localDoaArrUsed, truthLocalDoa, iSat, refFrameIdx);
  [fdRefErrHz, fdRmseHz] = localCalcTruthStateFdErr(truthEval.fdLocalRaw, truthFdLocal, iSat, refFrameIdx);
  rowCell{end + 1, 1} = struct( ...
    'localSatIdx', iSat, ...
    'globalSatIdx', localGetVectorElem(satIdxGlobal, iSat, NaN), ...
    'truthLocalDirRefErrDeg', refDirErrDeg, ...
    'truthRmsLocalDirErrDeg', rmsDirErrDeg, ...
    'truthMaxLocalDirErrDeg', maxDirErrDeg, ...
    'truthAzRefErrDeg', refAzErrDeg, ...
    'truthElRefErrDeg', refElErrDeg, ...
    'truthFdLocalRefErrHz', fdRefErrHz, ...
    'truthFdLocalRmseHz', fdRmseHz); %#ok<AGROW>
  maxLocalDirErrDeg = localMaxFinite(maxLocalDirErrDeg, maxDirErrDeg);
  maxFdLocalErrHz = localMaxFinite(maxFdLocalErrHz, abs(fdRefErrHz));
end

truthStateDiag.isEnabled = true;
truthStateDiag.refFrameIdx = refFrameIdx;
truthStateDiag.localDoaEval = truthEval.localDoaArrUsed;
truthStateDiag.fdLocalEval = truthEval.fdLocalRaw;
truthStateDiag.aliasBypass = true;
truthStateDiag.summaryTable = struct2table([rowCell{:}], 'AsArray', true);
truthStateDiag.maxLocalDirErrDeg = maxLocalDirErrDeg;
truthStateDiag.maxFdLocalErrHz = maxFdLocalErrHz;
end

function modelTruth = localBuildTruthStateModel(model)
%LOCALBUILDTRUTHSTATEMODEL Build a mapping-only model for truth diagnostics.
%
% The truth-state check is intended to validate the internal geometry and
% local-state mapping chain. It must therefore bypass ambiguity-branch
% priors so that the reported fdLocal error reflects only the raw mapping
% consistency, not a branch choice tied to the initializer.

modelTruth = model;
modelTruth.enableFdAliasUnwrap = false;
modelTruth.fdSatPriorHz = nan(model.numSat, 1);
end

function satOrderDiag = localBuildSatOrderDiag(model)
%LOCALBUILDSATORDERDIAG Pack compact satellite-order diagnostics.

numSat = model.numSat;
refSatIdxLocal = localResolveRefSatIdx(model);
satIdxGlobal = nan(numSat, 1);
refWeight = nan(numSat, 1);
frameCount = sum(model.frameMask, 2);

if ~isempty(model.sceneCell) && model.refFrameIdx >= 1 && model.refFrameIdx <= numel(model.sceneCell)
  sceneRef = model.sceneCell{model.refFrameIdx};
  if isfield(sceneRef, 'satIdx') && ~isempty(sceneRef.satIdx)
    satIdxGlobal(1:min(numSat, numel(sceneRef.satIdx))) = reshape(sceneRef.satIdx(1:min(numSat, numel(sceneRef.satIdx))), [], 1);
  end
  if isfield(sceneRef, 'ref') && isstruct(sceneRef.ref) && isfield(sceneRef.ref, 'weight') && ~isempty(sceneRef.ref.weight)
    refWeight(1:min(numSat, numel(sceneRef.ref.weight))) = reshape(sceneRef.ref.weight(1:min(numSat, numel(sceneRef.ref.weight))), [], 1);
  end
end

debugTruth = localGetStructField(model, 'debugTruth', struct());
satIdxGlobalDbg = reshape(localGetStructField(debugTruth, 'satIdxGlobal', []), [], 1);
refWeightDbg = reshape(localGetStructField(debugTruth, 'refWeight', []), [], 1);

numCopySat = min(numSat, numel(satIdxGlobalDbg));
fillMaskSat = ~isfinite(satIdxGlobal) & (1:numSat).' <= numCopySat;
satIdxGlobal(fillMaskSat) = satIdxGlobalDbg(fillMaskSat);

numCopyWeight = min(numSat, numel(refWeightDbg));
fillMaskWeight = ~isfinite(refWeight) & (1:numSat).' <= numCopyWeight;
refWeight(fillMaskWeight) = refWeightDbg(fillMaskWeight);

satOrderDiag = table((1:numSat).', satIdxGlobal, refWeight, frameCount(:), ...
  (1:numSat).' == refSatIdxLocal, 'VariableNames', ...
  {'localSatIdx', 'globalSatIdx', 'refWeight', 'numValidFrame', 'isRefSat'});
end

function [refDirErrDeg, rmsDirErrDeg, maxDirErrDeg, refAzErrDeg, refElErrDeg] = ...
  localCalcTruthStateLocalDoaErr(evalLocalDoa, truthLocalDoa, satIdx, refFrameIdx)
%LOCALCALCTRUTHSTATELOCALDOAERR Compare truth-point local steering states.

refDirErrDeg = NaN;
rmsDirErrDeg = NaN;
maxDirErrDeg = NaN;
refAzErrDeg = NaN;
refElErrDeg = NaN;

if isempty(evalLocalDoa) || isempty(truthLocalDoa)
  return;
end
if ndims(evalLocalDoa) ~= 3 || ndims(truthLocalDoa) ~= 3
  return;
end
if size(evalLocalDoa, 2) < satIdx || size(truthLocalDoa, 2) < satIdx
  return;
end
numFrame = min(size(evalLocalDoa, 3), size(truthLocalDoa, 3));
if numFrame < 1
  return;
end

errVec = nan(1, numFrame);
for iFrame = 1:numFrame
  evalDoa = evalLocalDoa(:, satIdx, iFrame);
  truthDoa = truthLocalDoa(:, satIdx, iFrame);
  if all(isfinite(evalDoa)) && all(isfinite(truthDoa))
    errVec(iFrame) = localCalcDirectionErrFromLocalDoa(evalDoa, truthDoa);
  end
end

refFrameIdx = min(max(refFrameIdx, 1), numFrame);
refDirErrDeg = errVec(refFrameIdx);
if any(isfinite(errVec))
  rmsDirErrDeg = sqrt(mean(errVec(isfinite(errVec)).^2));
  maxDirErrDeg = max(errVec(isfinite(errVec)));
end

refAzErrDeg = localWrapTo180Deg(rad2deg(evalLocalDoa(1, satIdx, refFrameIdx) - truthLocalDoa(1, satIdx, refFrameIdx)));
refElErrDeg = rad2deg(evalLocalDoa(2, satIdx, refFrameIdx) - truthLocalDoa(2, satIdx, refFrameIdx));
end

function [fdRefErrHz, fdRmseHz] = localCalcTruthStateFdErr(evalFdLocal, truthFdLocal, satIdx, refFrameIdx)
%LOCALCALCTRUTHSTATEFDERR Compare truth-point local Doppler states.

fdRefErrHz = NaN;
fdRmseHz = NaN;

if isempty(evalFdLocal) || isempty(truthFdLocal)
  return;
end
if size(evalFdLocal, 1) < satIdx || size(truthFdLocal, 1) < satIdx
  return;
end
numFrame = min(size(evalFdLocal, 2), size(truthFdLocal, 2));
if numFrame < 1
  return;
end

refFrameIdx = min(max(refFrameIdx, 1), numFrame);
fdRefErrHz = evalFdLocal(satIdx, refFrameIdx) - truthFdLocal(satIdx, refFrameIdx);
fdErrVec = evalFdLocal(satIdx, 1:numFrame) - truthFdLocal(satIdx, 1:numFrame);
validMask = isfinite(fdErrVec);
if any(validMask)
  fdRmseHz = sqrt(mean(fdErrVec(validMask).^2));
end
end

function dirErrDeg = localCalcDirectionErrFromLocalDoa(doaEval, doaTruth)
%LOCALCALCDIRECTIONERRFROMLOCALDOA Compute local unit-vector angle error.

vecEval = doa2dir(doaEval(:));
vecTruth = doa2dir(doaTruth(:));
cosVal = max(min(real(vecEval' * vecTruth), 1), -1);
dirErrDeg = rad2deg(acos(cosVal));
end

function valueOut = localMaxFinite(valueA, valueB)
%LOCALMAXFINITE Max helper that tolerates NaN inputs.

if ~isfinite(valueA)
  valueOut = valueB;
elseif ~isfinite(valueB)
  valueOut = valueA;
else
  valueOut = max(valueA, valueB);
end
end

function wrapDeg = localWrapTo180Deg(angleDeg)
%LOCALWRAPTO180DEG Wrap angles to [-180, 180].

wrapDeg = mod(angleDeg + 180, 360) - 180;
end

function [obj, pathGain, noiseVar, aux] = localNllPilotDynCon(model, optVar)
%LOCALNLLPILOTDYNCON Concentrated dynamic pilot negative log-likelihood.

hyp = localBuildHypothesis(model, optVar);
[obj, prof, profAux] = evalDoaDopplerDynProfileLike(model, hyp);

if nargout <= 1
  return;
end

pathGain = prof.pathGainEst;
noiseVar = prof.noiseVarEst;

aux = struct();
aux.doaParam = hyp.doaParam;
aux.eciAngle = hyp.eciAngle;
aux.eciDirection = hyp.eciDirection;
aux.latlon = hyp.latlon;
aux.userPosEci = hyp.userPosEci;
aux.userVelEci = hyp.userVelEci;

aux.refSatIdxLocal = localGetStructField(profAux, 'refSatIdxLocal', NaN);
aux.deltaFdSeq = hyp.deltaFdSeq;
aux.deltaFdRef = hyp.deltaFdRef;
aux.deltaFdRefRaw = localGetStructField(profAux, 'deltaFdRefRaw', hyp.deltaFdRef);
aux.deltaFdRefEval = localGetStructField(profAux, 'deltaFdRefEval', hyp.deltaFdRef);
aux.deltaFdRate = hyp.deltaFdRate;
aux.deltaFdRateRaw = localGetStructField(profAux, 'deltaFdRateRaw', hyp.deltaFdRate);
aux.deltaFdRateEval = localGetStructField(profAux, 'deltaFdRateEval', hyp.deltaFdRate);
aux.fdRef = hyp.fdRef;
aux.fdRate = hyp.fdRate;
aux.fdSat = hyp.fdSat;
aux.fdSatRaw = localGetStructField(profAux, 'fdSatRaw', hyp.fdSat);
aux.fdRateSat = hyp.fdRateSat;
aux.fdLocal = profAux.fdLocal;
aux.fdLocalRaw = localGetStructField(profAux, 'fdLocalRaw', profAux.fdLocal);
aux.fdSatEval = profAux.fdSatEval;
aux.fdLocalEval = profAux.fdLocalEval;
aux.fdAliasStepHz = profAux.fdAliasStepHz;
aux.fdAliasIndex = profAux.fdAliasIndex;
aux.fdAliasShiftHz = profAux.fdAliasShiftHz;

aux.localDoaArrUsed = profAux.localDoaArrUsed;
aux.localDoaArr = profAux.localDoaArrUsed;
aux.phaseSat = prof.phaseSatEst;
aux.framePhase = prof.framePhaseEst;
aux.ampEst = prof.ampEst;
aux.noiseVarGlobal = prof.noiseVarGlobal;
aux.residualNorm = prof.residualNorm;
aux.countPerSat = prof.countPerSat;
aux.residualSat = localGetStructField(prof, 'residualSat', []);
aux.fitValueSat = localGetStructField(prof, 'fitValueSat', []);
aux.objectiveSat = localGetStructField(prof, 'objectiveSat', []);
aux.blockValue = prof.blockValue;
aux.blockNorm2 = prof.blockNorm2;
aux.countPerBlock = localGetStructField(prof, 'countPerBlock', []);
aux.zMat = prof.zMat;
aux.etaMat = prof.etaMat;
end

function hyp = localBuildHypothesis(model, optVar)
%LOCALBUILDHYPOTHESIS Map one optimization point to geometry and Doppler.

optVar = optVar(:);
numVar = localGetNumOptVar(model);
if numel(optVar) ~= numVar
  error('estimatorDoaDopplerMlePilotMfOpt:InvalidOptVarSize', ...
    'The optimization vector must contain %d entries for fdRateMode = ''%s''.', ...
    numVar, model.fdRateMode);
end

doaParam = optVar(1:2);
fdRef = optVar(3);

switch model.fdRateMode
  case 'unknown'
    fdRate = optVar(4);
  case 'known'
    fdRate = model.fdRateKnown;
  case 'zero'
    fdRate = 0;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end

doaState = localBuildDoaState(model, doaParam);
fdSat = doaState.deltaFdRef + fdRef;
fdRateSat = doaState.deltaFdRate + fdRate;

hyp = doaState;
hyp.doaParam = doaParam(:);
hyp.fdRef = fdRef;
hyp.fdRate = fdRate;
hyp.fdSat = fdSat;
hyp.fdRateSat = fdRateSat;
end

function doaState = localBuildDoaState(model, doaParam)
%LOCALBUILDDOASTATE Map continuous DoA parameters to dynamic geometry.

doaParam = reshape(doaParam, 2, 1);

switch model.doaType
  case 'angle'
    doaState = localBuildAngleState(model, doaParam);

  case 'latlon'
    doaState = localBuildLatlonState(model, doaParam);

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidDoaType', ...
      'Unsupported doa type: %s.', model.doaType);
end
end

function doaState = localBuildAngleState(model, doaParam)
%LOCALBUILDANGLESTATE Build frame-wise states for angle-mode parameterization.

eciAngle = doaParam;
eciAngle(1) = mod(eciAngle(1), 2*pi);
eciDirection = doa2dir(eciAngle);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refRotMat = model.rotMatCell{model.steeringRefFrameIdx};
    localDoaCell = eciToAngleGrid(eciDirection, refRotMat);
    if ~iscell(localDoaCell)
      localDoaCell = {localDoaCell};
    end

    for iFrame = 1:model.numFrame
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end

  case 'framewise'
    for iFrame = 1:model.numFrame
      localDoaCell = eciToAngleGrid(eciDirection, model.rotMatCell{iFrame});
      if ~iscell(localDoaCell)
        localDoaCell = {localDoaCell};
      end

      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = localDoaCell{iSat};
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

deltaFdSeq = zeros(model.numSat, model.numFrame);
for iFrame = 1:model.numFrame
  relSatVel = model.satVelEci(:, :, iFrame) - model.refVelEci(:, iFrame);
  deltaFdSeq(:, iFrame) = ((eciDirection.' * relSatVel) / model.wavelength).';
end

[deltaFdRef, deltaFdRate] = localFitSatFdLine(model.timeOffsetSec, deltaFdSeq);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = [];
doaState.userPosEci = [];
doaState.userVelEci = [];
doaState.deltaFdSeq = deltaFdSeq;
doaState.deltaFdRef = deltaFdRef;
doaState.deltaFdRate = deltaFdRate;
doaState.localDoaArr = localDoaArr;
end

function doaState = localBuildLatlonState(model, doaParam)
%LOCALBUILDLATLONSTATE Build frame-wise states for latlon parameterization.

latlon = doaParam;
[userPosEci, userVelEci] = localLatlonToUserState(latlon, model.userStateRef);

eciDirection = userPosEci ./ vecnorm(userPosEci, 2, 1);
eciAngle = localDirectionToAngle(eciDirection);

localDoaArr = zeros(2, model.numSat, model.numFrame);

switch model.steeringMode
  case 'frozenref'
    refFrameIdx = model.steeringRefFrameIdx;
    refUserPos = userPosEci(:, refFrameIdx);

    for iSat = 1:model.numSat
      localDoaArr(:, iSat, :) = repmat( ...
        globalToLocalDoa(refUserPos, model.satPosEci(:, iSat, refFrameIdx), ...
        model.rotMatCell{refFrameIdx}{iSat}), ...
        1, 1, model.numFrame);
    end

  case 'framewise'
    for iFrame = 1:model.numFrame
      currentUserPos = userPosEci(:, iFrame);
      for iSat = 1:model.numSat
        localDoaArr(:, iSat, iFrame) = globalToLocalDoa( ...
          currentUserPos, model.satPosEci(:, iSat, iFrame), ...
          model.rotMatCell{iFrame}{iSat});
      end
    end

  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidSteeringMode', ...
      'Unsupported steering mode: %s.', model.steeringMode);
end

refDopplerState = buildReferenceDopplerState( ...
  model.refScene, model.satPosEci, model.satVelEci, ...
  userPosEci, userVelEci, model.wavelength, model.timeOffsetSec, model.refSatIdxLocal);

doaState = struct();
doaState.eciAngle = eciAngle;
doaState.eciDirection = eciDirection;
doaState.latlon = latlon;
doaState.userPosEci = userPosEci;
doaState.userVelEci = userVelEci;
doaState.deltaFdSeq = refDopplerState.deltaFd;
doaState.deltaFdRef = refDopplerState.deltaFdRef;
doaState.deltaFdRate = refDopplerState.deltaFdRate;
doaState.fdSatSeq = refDopplerState.fdSat;
doaState.fdRefSeq = refDopplerState.fdRef;
doaState.fdRefGeom = refDopplerState.fdRefRefFrame;
doaState.fdSatGeom = refDopplerState.fdSatRefFrame;
doaState.deltaFdGeom = refDopplerState.deltaFdRefFrame;
doaState.localDoaArr = localDoaArr;
end

function [userPosEci, userVelEci] = localLatlonToUserState(latlon, userStateRef)
%LOCALLATLONTOUSERSTATE Convert one ground point [lat; lon] into motion-aware ECI states.

[userPosEci, userVelEci] = buildUserStateFromLatlon(latlon, userStateRef);
end

function angle = localDirectionToAngle(direction)
%LOCALDIRECTIONTOANGLE Convert ECI unit directions to [az; el].

direction = direction ./ vecnorm(direction, 2, 1);
angle = [mod(atan2(direction(2, :), direction(1, :)), 2*pi); ...
         asin(max(min(direction(3, :), 1), -1))];
end

function elemValue = localGetVectorElem(vec, idx, defaultValue)
%LOCALGETVECTORELEM Read one vector element with a default fallback.

elemValue = defaultValue;
if isempty(vec)
  return;
end
vec = vec(:);
if idx >= 1 && idx <= numel(vec)
  elemValue = vec(idx);
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

function numVar = localGetNumOptVar(model)
%LOCALGETNUMOPTVAR Get the current optimization vector size.

switch model.fdRateMode
  case 'unknown'
    numVar = 4;
  case {'known', 'zero'}
    numVar = 3;
  otherwise
    error('estimatorDoaDopplerMlePilotMfOpt:InvalidFdRateMode', ...
      'Unsupported fdRate mode: %s.', model.fdRateMode);
end
end
