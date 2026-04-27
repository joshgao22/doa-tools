function [estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotMfOpt( ...
  sceneSeq, rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, ...
  fdRange, fdRateRange, initParam, verbose, modelOpt)
%ESTIMATORDOADOPPLERMLEPILOTMFOPT Multi-frame pilot-based DoA-Doppler MLE.
% Keeps the main entry focused on orchestration:
%   1) build one unified MF model;
%   2) build one robust initializer;
%   3) solve CP/IP branches through the shared MF solver helper;
%   4) assemble one compact unified result structure.
%
% The heavy model/init/solver/debug logic is delegated to:
%   - buildDoaDopplerMfModel
%   - buildDoaDopplerMfInit
%   - solveDoaDopplerMfBranches
%   - buildDoaDopplerMfDebug
%
% This file intentionally keeps only entry-layer glue and result packing so
% that future algorithm changes happen in the shared helpers instead of
% being duplicated between entry code and legacy local functions.
%
%See also:
%  buildDoaDopplerMfModel, buildDoaDopplerMfInit,
%  solveDoaDopplerMfBranches, buildDoaDopplerMfDebug,
%  evalDoaDopplerDynProfileLike

arguments (Input)
  sceneSeq (1,1) struct
  rxSig
  pilotWave
  carrierFreq (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  sampleRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  fdRange = []
  fdRateRange = []
  initParam = []
  verbose (1,1) logical = false
  modelOpt (1,1) struct = struct()
end

if isfield(modelOpt, 'verbose') && ~isempty(modelOpt.verbose)
  verbose = verbose || logical(modelOpt.verbose);
end
verbose = verbose || localResolveInheritedVerbose();
modelOpt.verbose = verbose;

[model, ~, ~, modelOpt] = buildDoaDopplerMfModel(sceneSeq, rxSig, pilotWave, ...
  carrierFreq, sampleRate, doaGrid, fdRange, fdRateRange, modelOpt);
[initParam, initDiag, model] = buildDoaDopplerMfInit(model, initParam);
[solveUse, optimInfo, initEvalDiag] = solveDoaDopplerMfBranches(model, initParam, verbose);
debugAux = buildDoaDopplerMfDebug(model, initDiag, initEvalDiag, ...
  solveUse.finalEvalDiag, solveUse.debugTrace);

estResult = localBuildEstResult(model, solveUse, initParam, optimInfo, ...
  initDiag, initEvalDiag, debugAux);
pathGain = solveUse.pathGain;
noiseVar = solveUse.noiseVar;
end


function estResult = localBuildEstResult(model, solveUse, initParam, optimInfo, ...
  initDiag, initEvalDiag, debugAux)
%LOCALBUILDESTRESULT Build the unified estimator output structure.

estAux = solveUse.estAux;

estResult = struct();
estResult.modelType = localGetModelType(model);
estResult.phaseMode = model.phaseMode;
estResult.fdRateMode = model.fdRateMode;
estResult.doaType = model.doaType;
estResult.numSource = 1;
estResult.numSat = model.numSat;
estResult.numFrame = model.numFrame;
estResult.timeOffsetSec = model.timeOffsetSec;
estResult.frameMask = model.frameMask;

estResult.doaParamEst = estAux.doaParam;
estResult.fdRefEst = estAux.fdRef(:);
estResult.fdRateEst = estAux.fdRate(:);
estResult.pathGainEst = solveUse.pathGain;
estResult.noiseVarEst = solveUse.noiseVar;

estResult.initParam = initParam(:);
estResult.optVarEst = solveUse.optVar(:);
estResult.fval = solveUse.fval;
estResult.exitflag = solveUse.exitflag;
estResult.isResolved = logical(localGetStructField(solveUse, 'isResolved', false));
% Expose the final branch identity directly on estResult so downstream
% flow-level summary / regression code does not need to re-open optimInfo
% just to learn which CP-U branch won the final selection.
estResult.solveVariant = string(localGetStructField(optimInfo, 'solveVariant', ""));
estResult.candidateVariant = localGetStructField(optimInfo, 'candidateVariant', strings(0, 1));
estResult.candidateObjective = localGetStructField(optimInfo, 'candidateObjective', []);
estResult.optimInfo = optimInfo;
estResult.aux = localBuildAux(model, estAux, solveUse.pathGain, solveUse.noiseVar, ...
  initDiag, initEvalDiag, solveUse.finalEvalDiag, solveUse.debugTrace, debugAux);
end


function modelType = localGetModelType(model)
%LOCALGETMODELTYPE Get the estimator model tag for script-level dispatch.

phaseTagMap = struct();
phaseTagMap.continuous = 'Cp';
phaseTagMap.relaxed = 'Relaxed';
phaseTagMap.independent = 'Ip';

fdRateTagMap = struct();
fdRateTagMap.unknown = 'UnknownRate';
fdRateTagMap.known = 'KnownRate';
fdRateTagMap.zero = 'ZeroRate';

modelType = ['mf', phaseTagMap.(model.phaseMode), fdRateTagMap.(model.fdRateMode)];
end


function aux = localBuildAux(model, estAux, pathGain, noiseVar, ...
  initDiag, initEvalDiag, finalEvalDiag, debugTrace, debugAux)
%LOCALBUILDAUX Build auxiliary geometry and compact diagnostics.

aux = struct();

if model.numSat == 1
  aux.doaGrid = model.doaGrid{1};
else
  aux.doaGrid = model.doaGrid;
end

aux.modelType = localGetModelType(model);
aux.doaType = model.doaType;
aux.globalFrame = model.globalFrame;
aux.phaseMode = model.phaseMode;
aux.phaseGridCount = model.phaseGridCount;
aux.phaseRefine = model.phaseRefine;
aux.fdRateMode = model.fdRateMode;
aux.fdRateKnown = model.fdRateKnown;
aux.steeringMode = model.steeringMode;
aux.enableFdAliasUnwrap = model.enableFdAliasUnwrap;
aux.fdSatPriorHz = model.fdSatPriorHz;
aux.refSatIdxGlobal = localGetStructField(model, 'refSatIdxGlobal', NaN);

aux.eciAngleEst = estAux.eciAngle;
aux.eciDirectionEst = estAux.eciDirection;
aux.latlonEst = estAux.latlon;
aux.userPosEciEst = estAux.userPosEci;
aux.userVelEciEst = estAux.userVelEci;

aux.refSatIdxLocal = localGetStructField(estAux, 'refSatIdxLocal', NaN);
aux.deltaFdSeqEst = estAux.deltaFdSeq;
aux.deltaFdRefEst = estAux.deltaFdRef;
aux.deltaFdRateEst = estAux.deltaFdRate;
aux.fdSatEst = estAux.fdSat;
aux.fdRateSatEst = estAux.fdRateSat;
aux.fdLocalEst = estAux.fdLocal;
aux.localDoaEst = estAux.localDoaArr;

aux.phaseSatEst = estAux.phaseSat;
aux.framePhaseEst = estAux.framePhase;
aux.ampEst = estAux.ampEst;
aux.blockValue = estAux.blockValue;
aux.blockNorm2 = localGetStructField(estAux, 'blockNorm2', []);
aux.countPerBlock = localGetStructField(estAux, 'countPerBlock', []);
aux.zMat = estAux.zMat;
aux.etaMat = estAux.etaMat;
aux.residualNorm = estAux.residualNorm;
aux.residualSatEst = localGetStructField(estAux, 'residualSat', []);
aux.fitValueSatEst = localGetStructField(estAux, 'fitValueSat', []);
aux.objectiveSatEst = localGetStructField(estAux, 'objectiveSat', []);
aux.effectiveFrameSupportSat = localGetStructField(estAux, 'effectiveFrameSupportSat', []);
aux.effectiveFrameSupportRatioSat = localGetStructField(estAux, 'effectiveFrameSupportRatioSat', []);
aux.negativeProjectionRatioSat = localGetStructField(estAux, 'negativeProjectionRatioSat', []);
aux.collapsedFrameCountSat = localGetStructField(estAux, 'collapsedFrameCountSat', []);
aux.satFitRatio = localGetStructField(estAux, 'satFitRatio', []);
aux.refFitRatio = localGetStructField(estAux, 'refFitRatio', NaN);
aux.nonRefFitRatioFloor = localGetStructField(estAux, 'nonRefFitRatioFloor', NaN);
aux.refSupportRatio = localGetStructField(estAux, 'refSupportRatio', NaN);
aux.nonRefSupportRatioFloor = localGetStructField(estAux, 'nonRefSupportRatioFloor', NaN);
aux.refConsistencyNorm = localGetStructField(estAux, 'refConsistencyNorm', NaN);
aux.nonRefConsistencyRatioFloor = localGetStructField(estAux, 'nonRefConsistencyRatioFloor', NaN);
aux.maxNonRefNegativeProjectionRatio = localGetStructField(estAux, 'maxNonRefNegativeProjectionRatio', NaN);
aux.nonRefFitFloorPenalty = localGetStructField(estAux, 'nonRefFitFloorPenalty', 0);
aux.nonRefSupportFloorPenalty = localGetStructField(estAux, 'nonRefSupportFloorPenalty', 0);
aux.additionalObjectivePenalty = localGetStructField(estAux, 'additionalObjectivePenalty', 0);
aux.noiseVarGlobal = estAux.noiseVarGlobal;
aux.countPerSat = estAux.countPerSat;
aux.pathGain = pathGain;
aux.noiseVar = noiseVar;

aux.isStaticModel = strcmp(model.fdRateMode, 'zero');
aux.isDiscrete = false;
aux.debug = debugAux;
aux.debug.initDiag = initDiag;
aux.debug.initEval = initEvalDiag;
aux.debug.finalEval = finalEvalDiag;
aux.debug.evalTrace = localGetStructField(debugTrace, 'traceTable', []);
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

function verbose = localResolveInheritedVerbose()
%LOCALRESOLVEINHERITEDVERBOSE Inherit one repo-level verbose flag.

verbose = false;
try
  if isappdata(0, 'doaToolsVerbose')
    verbose = logical(getappdata(0, 'doaToolsVerbose'));
  end
catch
  verbose = false;
end

if verbose
  return;
end

try
  hasVerbose = evalin('base', "exist(''doaToolsVerbose'', ''var'')");
  if isequal(hasVerbose, 1)
    verbose = logical(evalin('base', 'doaToolsVerbose'));
  end
catch
  verbose = false;
end
end

