function [estResult, pathGain, noiseVar] = estimatorDoaDopplerMlePilotSfOpt(scene, ...
  rxSig, pilotWave, carrierFreq, sampleRate, doaGrid, fdRange, ...
  numSource, initParam, verbose, modelOpt)
%ESTIMATORDOADOPPLERMLEPILOTSFOPT Single-frame pilot-based DoA-Doppler MLE.
% Keeps the main entry focused on stage orchestration:
%   1) build one unified SF model;
%   2) build one DoA/fd initializer;
%   3) evaluate or optimize one concentrated profile objective;
%   4) assemble one unified estimator result.
%
% The heavy parsing, model-building, and profile-likelihood logic lives in:
%   - buildDoaDopplerSfModel
%   - buildDoaDopplerSfInit
%   - evalDoaDopplerSfProfileLike
%
% This keeps the entry file short and reduces the maintenance cost of
% future static-estimator changes.
%
%See also:
%  buildDoaDopplerSfModel, buildDoaDopplerSfInit,
%  evalDoaDopplerSfProfileLike, estimatorDoaDopplerMlePilotMfOpt

arguments (Input)
  scene (1,1) struct
  rxSig
  pilotWave {mustBeNumeric, mustBeFinite}
  carrierFreq (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  sampleRate (1,1) {mustBePositive, mustBeNumeric, mustBeFinite}
  doaGrid (1,:) {mustBeA(doaGrid, ["struct", "cell"])}
  fdRange = []
  numSource (1,1) {mustBePositive, mustBeInteger} = 1
  initParam = []
  verbose (1,1) logical = false
  modelOpt (1,1) struct = struct()
end

[model, lb, ub, modelOpt] = buildDoaDopplerSfModel(scene, rxSig, pilotWave, ...
  carrierFreq, sampleRate, doaGrid, fdRange, numSource, modelOpt);
initParam = buildDoaDopplerSfInit(model, initParam);

if verbose
  displayMode = 'iter';
else
  displayMode = 'off';
end

baseOptimOpt = optimoptions('fmincon', ...
  'Display', displayMode, ...
  'Algorithm', 'interior-point');

baseOptimField = fieldnames(model.optimOpt);
for iField = 1:numel(baseOptimField)
  baseOptimOpt.(baseOptimField{iField}) = model.optimOpt.(baseOptimField{iField});
end

objFun = @(optVar) evalDoaDopplerSfProfileLike(model, optVar);
runTimer = tic;
if model.evalOnly
  output = struct();
  optVar = localSortOptVar(model, initParam);
  [fval, pathGain, noiseVar, estAux] = evalDoaDopplerSfProfileLike(model, optVar);
  runTimeSec = toc(runTimer);
  exitflag = 99;
  optimInfo = localBuildEvalOnlyOptimInfo(runTimeSec, model);
else
  [optVar, fval, exitflag, output] = fmincon( ...
    objFun, initParam, [], [], [], [], lb, ub, [], baseOptimOpt);
  runTimeSec = toc(runTimer);
  optVar = localSortOptVar(model, optVar);
  [~, pathGain, noiseVar, estAux] = evalDoaDopplerSfProfileLike(model, optVar);
  optimInfo = localBuildOptimInfo(output, exitflag, runTimeSec, baseOptimOpt, model);
end

estResult = localBuildEstResult(model, estAux, pathGain, noiseVar, ...
  initParam, optVar, fval, exitflag, optimInfo);
end


function optVar = localSortOptVar(model, optVar)
%LOCALSORTOPTVAR Sort sources by DoA so output ordering stays stable.

if model.numSource <= 1
  optVar = optVar(:);
  return;
end

numDoaVar = 2 * model.numSource;
doaParam = reshape(optVar(1:numDoaVar), 2, model.numSource);
fdRef = optVar(numDoaVar + 1:end);
[~, sortIdx] = sortrows(doaParam.', [1, 2]);
doaParam = doaParam(:, sortIdx);
fdRef = fdRef(sortIdx);
optVar = [doaParam(:); fdRef(:)];
end


function estResult = localBuildEstResult(model, estAux, pathGain, noiseVar, ...
  initParam, optVar, fval, exitflag, optimInfo)
%LOCALBUILDESTRESULT Build the unified estimator output structure.

estResult = struct();
estResult.modelType = 'sfStatic';
estResult.phaseMode = 'singleFrame';
estResult.fdRateMode = 'zero';
estResult.doaType = model.doaType;
estResult.numSource = model.numSource;
estResult.numSat = model.numSat;
estResult.numFrame = 1;
estResult.timeOffsetSec = 0;
estResult.frameMask = true(model.numSat, 1);

estResult.doaParamEst = estAux.doaParam;
estResult.fdRefEst = estAux.fdRef(:);
estResult.fdRateEst = zeros(model.numSource, 1);
estResult.pathGainEst = pathGain;
estResult.noiseVarEst = noiseVar;

estResult.initParam = initParam(:);
estResult.optVarEst = optVar(:);
estResult.fval = fval;
estResult.exitflag = exitflag;
estResult.isResolved = localCheckResolved(optVar, pathGain, noiseVar, fval, exitflag);
estResult.optimInfo = optimInfo;
estResult.aux = localBuildAux(model, estAux, pathGain);
end


function aux = localBuildAux(model, estAux, pathGain)
%LOCALBUILDAUX Build auxiliary geometry and diagnostic outputs.

aux = struct();

if model.numSat == 1
  aux.doaGrid = model.doaGrid{1};
else
  aux.doaGrid = model.doaGrid;
end

aux.globalFrame = model.globalFrame;
aux.phaseMode = 'singleFrame';
aux.fdRateMode = 'zero';
aux.fdRateKnown = 0;
aux.steeringMode = 'singleFrame';
aux.refSatIdxLocal = model.refSatIdxLocal;
aux.refSatIdxGlobal = model.refSatIdxGlobal;
aux.eciAngleEst = estAux.eciAngle;
aux.eciDirectionEst = estAux.eciDirection;
aux.latlonEst = estAux.latlon;
aux.userPosEciEst = estAux.userPosEci;
aux.userVelEciEst = estAux.userVelEci;

aux.deltaFdSeqEst = estAux.deltaFd;
aux.deltaFdRefEst = estAux.deltaFd;
aux.deltaFdRateEst = zeros(size(estAux.deltaFd));
aux.fdSatEst = estAux.fd;
aux.fdLocalEst = estAux.fd;
aux.localDoaEst = estAux.localDoaArr;

aux.phaseSatEst = angle(pathGain);
aux.framePhaseEst = angle(pathGain);
aux.ampEst = abs(pathGain);
aux.blockValue = nan(size(pathGain));
aux.zMat = nan(size(pathGain));
aux.etaMat = nan(size(pathGain));
aux.residualNorm = estAux.residualNorm;
aux.residualNormSat = estAux.residualNormSat;
aux.objectiveSat = estAux.objectiveSat;
aux.noiseVarSat = estAux.noiseVarSat;
aux.noiseVarGlobal = estAux.noiseVarGlobal;
aux.countPerSat = estAux.countPerSat;
aux.satWeight = estAux.satWeight;

aux.isStaticModel = true;
aux.isDiscrete = false;
end


function isResolved = localCheckResolved(optVar, pathGain, noiseVar, fval, exitflag)
%LOCALCHECKRESOLVED Check whether the optimization result is usable.

isResolved = exitflag > 0 && ...
  isfinite(fval) && ...
  all(isfinite(optVar(:))) && ...
  all(isfinite(pathGain(:))) && ...
  all(isfinite(noiseVar(:)));
end


function optimInfo = localBuildEvalOnlyOptimInfo(runTimeSec, model)
%LOCALBUILDEVALONLYOPTIMINFO Build one compact eval-only summary.

optimInfo = struct();
optimInfo.iterations = 0;
optimInfo.funcCount = 1;
optimInfo.constrviolation = 0;
optimInfo.stepsize = 0;
optimInfo.firstorderopt = NaN;
optimInfo.algorithm = 'evalOnly';
optimInfo.message = 'Objective evaluated at one fixed point.';
optimInfo.exitflag = 99;
optimInfo.runTimeSec = runTimeSec;
optimInfo.phaseMode = 'singleFrame';
optimInfo.fdRateMode = 'zero';
optimInfo.fdRateKnown = 0;
optimInfo.numVar = 3 * model.numSource;
optimInfo.display = 'off';
optimInfo.modelType = 'sfStatic';
optimInfo.solver = 'none';
end


function optimInfo = localBuildOptimInfo(output, exitflag, runTimeSec, optimOpt, model)
%LOCALBUILDOPTIMINFO Build a compact optimization summary for logging.

optimInfo = struct();
optimInfo.iterations = localGetStructField(output, 'iterations', NaN);
optimInfo.funcCount = localGetStructField(output, 'funcCount', NaN);
optimInfo.constrviolation = localGetStructField(output, 'constrviolation', NaN);
optimInfo.stepsize = localGetStructField(output, 'stepsize', NaN);
optimInfo.firstorderopt = localGetStructField(output, 'firstorderopt', NaN);
optimInfo.algorithm = localGetStructField(output, 'algorithm', '');
optimInfo.message = localGetStructField(output, 'message', '');
optimInfo.exitflag = exitflag;
optimInfo.runTimeSec = runTimeSec;
optimInfo.phaseMode = 'singleFrame';
optimInfo.fdRateMode = 'zero';
optimInfo.fdRateKnown = 0;
optimInfo.numVar = 3 * model.numSource;
optimInfo.display = localGetStructField(optimOpt, 'Display', '');
optimInfo.modelType = 'sfStatic';
optimInfo.solver = 'fmincon';
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
