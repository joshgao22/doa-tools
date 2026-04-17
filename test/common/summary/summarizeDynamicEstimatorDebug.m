function dynSummary = summarizeDynamicEstimatorDebug(caseResult, sceneSeq, truth)
%SUMMARIZEDYNAMICESTIMATORDEBUG Build compact dynamic debug summaries.
% This helper collects the dynamic-estimator tables that were previously
% assembled inside the development script. Keeping them here makes the dev
% script shorter while preserving the existing diagnostics.
%
%Syntax:
%  dynSummary = summarizeDynamicEstimatorDebug(caseResult, sceneSeq, truth)
%
%Output:
%  dynSummary - struct with fields:
%                 .diagTable
%                 .objectiveTable
%                 .perSatTable
%                 .blockTable
%                 .localStateTable
%                 .truthStateTable
%                 .satOrderTable

arguments
  caseResult
  sceneSeq struct
  truth struct
end

dynSummary = struct();
dynSummary.diagTable = localBuildDynamicDiagTable(caseResult, truth);
dynSummary.objectiveTable = localBuildDynamicObjectiveTable(caseResult);
dynSummary.perSatTable = localBuildDynamicPerSatTable(caseResult, truth);
dynSummary.blockTable = localBuildDynamicBlockTable(caseResult);
dynSummary.localStateTable = localBuildDynamicLocalStateTable(caseResult, sceneSeq, truth);
dynSummary.truthStateTable = localBuildDynamicTruthStateTable(caseResult);
dynSummary.satOrderTable = localBuildDynamicSatOrderTable(caseResult);
end

function dynDiagTable = localBuildDynamicDiagTable(caseResult, truth)
%LOCALBUILDDYNAMICDIAGTABLE Build diagnostics for dynamic CP estimators.

keepMask = false(1, numel(caseResult));
for iCase = 1:numel(caseResult)
  keepMask(iCase) = contains(caseResult(iCase).displayName, "MF-CP");
end

caseDyn = caseResult(keepMask);
numCase = numel(caseDyn);
infoCell = cell(1, numCase);
for iCase = 1:numCase
  caseInfo = caseDyn(iCase);
  estResult = caseInfo.estResult;
  initParam = localGetFieldOrDefault(estResult, 'initParam', []);
  initLatlon = nan(2, 1);
  initFdRef = NaN;
  initFdRate = NaN;
  if ~isempty(initParam)
    initParam = initParam(:);
    if numel(initParam) >= 2
      initLatlon = initParam(1:2);
    end
    if numel(initParam) >= 3
      initFdRef = initParam(3);
    end
    if numel(initParam) >= 4
      initFdRate = initParam(4);
    end
  end

  [initAngleErrDeg, ~, ~] = calcLatlonAngleError(initLatlon, truth.latlonTrueDeg);
  [finalAngleErrDeg, ~, ~] = calcLatlonAngleError(localGetLatlonEst(estResult), truth.latlonTrueDeg);

  fdSatRmseHz = NaN;
  if isfield(estResult, 'aux') && isfield(estResult.aux, 'fdSatEst') && ...
      isequal(size(estResult.aux.fdSatEst), size(truth.fdSatSeries))
    fdSatRmseHz = sqrt(mean((estResult.aux.fdSatEst(:) - truth.fdSatSeries(:)) .^ 2));
  end

  residualNorm = NaN;
  if isfield(estResult, 'aux') && isfield(estResult.aux, 'residualNorm')
    residualNorm = estResult.aux.residualNorm;
  end

  infoCell{iCase} = struct( ...
    'displayName', caseInfo.displayName, ...
    'satMode', caseInfo.satMode, ...
    'dynamicMode', caseInfo.dynamicMode, ...
    'initAngleErrDeg', initAngleErrDeg, ...
    'finalAngleErrDeg', finalAngleErrDeg, ...
    'initFdRefErrHz', initFdRef - truth.fdRefFit, ...
    'finalFdRefErrHz', localGetFieldOrDefault(estResult, 'fdRefEst', NaN) - truth.fdRefFit, ...
    'initFdRateErrHzPerSec', initFdRate - truth.fdRateFit, ...
    'finalFdRateErrHzPerSec', localGetFieldOrDefault(estResult, 'fdRateEst', NaN) - truth.fdRateFit, ...
    'fdSatRmseHz', fdSatRmseHz, ...
    'residualNorm', residualNorm, ...
    'exitflag', localGetFieldOrDefault(estResult, 'exitflag', NaN), ...
    'isResolved', logical(localGetFieldOrDefault(estResult, 'isResolved', false)) ...
    );
end

dynDiagTable = struct2table([infoCell{:}], 'AsArray', true);
end




function dynObjTable = localBuildDynamicObjectiveTable(caseResult)
%LOCALBUILDDYNAMICOBJECTIVETABLE Summarize objective reduction for dynamic cases.

caseDyn = localSelectDynamicCase(caseResult);
numCase = numel(caseDyn);
infoCell = cell(1, numCase);
for iCase = 1:numCase
  caseInfo = caseDyn(iCase);
  estResult = caseInfo.estResult;
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(estResult, 'aux', struct()), 'debug', struct());
  initEval = localGetFieldOrDefault(debugAux, 'initEval', struct());
  finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
  probe = localGetFieldOrDefault(debugAux, 'probe', struct());

  initObj = localGetFieldOrDefault(initEval, 'obj', NaN);
  finalObj = localGetFieldOrDefault(finalEval, 'obj', NaN);
  initResidual = localGetFieldOrDefault(initEval, 'residualNorm', NaN);
  finalResidual = localGetFieldOrDefault(finalEval, 'residualNorm', NaN);
  truthObj = localGetFieldOrDefault(localGetFieldOrDefault(probe, 'truth', struct()), 'obj', NaN);
  mixDoaTrueObj = localGetFieldOrDefault(localGetFieldOrDefault(probe, 'mixDoaTrue', struct()), 'obj', NaN);
  mixFdTrueObj = localGetFieldOrDefault(localGetFieldOrDefault(probe, 'mixFdTrue', struct()), 'obj', NaN);

  infoCell{iCase} = struct( ...
    'displayName', caseInfo.displayName, ...
    'satMode', caseInfo.satMode, ...
    'dynamicMode', caseInfo.dynamicMode, ...
    'initObj', initObj, ...
    'finalObj', finalObj, ...
    'truthObj', truthObj, ...
    'mixDoaTrueObj', mixDoaTrueObj, ...
    'mixFdTrueObj', mixFdTrueObj, ...
    'objImprove', initObj - finalObj, ...
    'finalGapToTruth', finalObj - truthObj, ...
    'mixDoaGapToTruth', mixDoaTrueObj - truthObj, ...
    'mixFdGapToTruth', mixFdTrueObj - truthObj, ...
    'initResidualNorm', initResidual, ...
    'finalResidualNorm', finalResidual, ...
    'residualImprove', initResidual - finalResidual ...
    );
end

dynObjTable = struct2table([infoCell{:}], 'AsArray', true);
end



function dynSatTable = localBuildDynamicPerSatTable(caseResult, truth)
%LOCALBUILDDYNAMICPERSATTABLE Summarize residual and coherence per satellite.

caseDyn = localSelectDynamicCase(caseResult);
numCase = numel(caseDyn);
numSat = size(truth.fdSatSeries, 1);
[truthDeltaFdRefHz, truthDeltaFdRateHzPerSec] = localBuildTruthDeltaFdFit(truth);
infoCell = cell(1, numCase);

for iCase = 1:numCase
  caseInfo = caseDyn(iCase);
  estResult = caseInfo.estResult;
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(estResult, 'aux', struct()), 'debug', struct());
  initEval = localGetFieldOrDefault(debugAux, 'initEval', struct());
  finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());

  info = struct();
  info.displayName = caseInfo.displayName;
  info.satMode = caseInfo.satMode;
  info.dynamicMode = caseInfo.dynamicMode;
  for iSat = 1:numSat
    info.(sprintf('initResidualSat%d', iSat)) = localGetVectorElem( ...
      localGetFieldOrDefault(initEval, 'residualSat', []), iSat, NaN);
    info.(sprintf('finalResidualSat%d', iSat)) = localGetVectorElem( ...
      localGetFieldOrDefault(finalEval, 'residualSat', []), iSat, NaN);
    info.(sprintf('initCohSat%d', iSat)) = localGetVectorElem( ...
      localGetFieldOrDefault(initEval, 'coherenceSat', []), iSat, NaN);
    info.(sprintf('finalCohSat%d', iSat)) = localGetVectorElem( ...
      localGetFieldOrDefault(finalEval, 'coherenceSat', []), iSat, NaN);
    initDeltaFdRefVec = localGetFieldOrDefault(initEval, 'deltaFdRefEval', ...
      localGetFieldOrDefault(initEval, 'deltaFdRef', []));
    finalDeltaFdRefVec = localGetFieldOrDefault(finalEval, 'deltaFdRefEval', ...
      localGetFieldOrDefault(finalEval, 'deltaFdRef', []));
    initDeltaFdRateVec = localGetFieldOrDefault(initEval, 'deltaFdRateEval', ...
      localGetFieldOrDefault(initEval, 'deltaFdRate', []));
    finalDeltaFdRateVec = localGetFieldOrDefault(finalEval, 'deltaFdRateEval', ...
      localGetFieldOrDefault(finalEval, 'deltaFdRate', []));
    info.(sprintf('initDeltaFdRefErrSat%dHz', iSat)) = localGetVectorElem(initDeltaFdRefVec, iSat, NaN) - truthDeltaFdRefHz(iSat);
    info.(sprintf('finalDeltaFdRefErrSat%dHz', iSat)) = localGetVectorElem(finalDeltaFdRefVec, iSat, NaN) - truthDeltaFdRefHz(iSat);
    info.(sprintf('initDeltaFdRateErrSat%dHzPerSec', iSat)) = localGetVectorElem(initDeltaFdRateVec, iSat, NaN) - truthDeltaFdRateHzPerSec(iSat);
    info.(sprintf('finalDeltaFdRateErrSat%dHzPerSec', iSat)) = localGetVectorElem(finalDeltaFdRateVec, iSat, NaN) - truthDeltaFdRateHzPerSec(iSat);
  end
  infoCell{iCase} = info;
end

dynSatTable = struct2table([infoCell{:}], 'AsArray', true);
end




function dynBlockTable = localBuildDynamicBlockTable(caseResult)
%LOCALBUILDDYNAMICBLOCKTABLE Summarize compact per-block diagnostics.

caseDyn = localSelectDynamicCase(caseResult);
numCase = numel(caseDyn);
infoCell = cell(1, numCase);
for iCase = 1:numCase
  caseInfo = caseDyn(iCase);
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseInfo.estResult, 'aux', struct()), 'debug', struct());
  finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
  blockResidual = localGetFieldOrDefault(finalEval, 'blockResidualMat', []);
  phaseResid = localGetFieldOrDefault(finalEval, 'blockPhaseResidMat', []);

  info = struct();
  info.displayName = caseInfo.displayName;
  info.satMode = caseInfo.satMode;
  info.dynamicMode = caseInfo.dynamicMode;
  if ~isempty(blockResidual)
    info.maxBlockResidual = max(blockResidual(:), [], 'omitnan');
    info.rmsBlockResidual = sqrt(mean(blockResidual(:).^2, 'omitnan'));
  else
    info.maxBlockResidual = NaN;
    info.rmsBlockResidual = NaN;
  end
  if ~isempty(phaseResid)
    info.maxAbsPhaseResidRad = max(abs(phaseResid(:)), [], 'omitnan');
    info.rmsPhaseResidRad = sqrt(mean(phaseResid(:).^2, 'omitnan'));
  else
    info.maxAbsPhaseResidRad = NaN;
    info.rmsPhaseResidRad = NaN;
  end
  infoCell{iCase} = info;
end

dynBlockTable = struct2table([infoCell{:}], 'AsArray', true);
end



function dynLocalStateTable = localBuildDynamicLocalStateTable(caseResult, sceneSeq, truth)
%LOCALBUILDDYNAMICLOCALSTATETABLE Summarize local DoA/fd state mapping.

caseDyn = localSelectDynamicCase(caseResult);
truthLocalDoa = localExtractSceneLocalDoa(sceneSeq);
if nargin < 3
  truth = struct();
end
numCase = numel(caseDyn);
infoCell = {};
for iCase = 1:numCase
  caseInfo = caseDyn(iCase);
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseInfo.estResult, 'aux', struct()), 'debug', struct());
  initEval = localGetFieldOrDefault(debugAux, 'initEval', struct());
  finalEval = localGetFieldOrDefault(debugAux, 'finalEval', struct());
  initLocalDoa = localAlignLocalDoaForCase(localGetFieldOrDefault(initEval, 'localDoaArrUsed', []), truthLocalDoa, caseInfo);
  finalLocalDoa = localAlignLocalDoaForCase(localGetFieldOrDefault(finalEval, 'localDoaArrUsed', []), truthLocalDoa, caseInfo);
  truthLocalDoaCase = localAlignLocalDoaForCase([], truthLocalDoa, caseInfo);
  initFdLocal = localGetFieldOrDefault(initEval, 'fdLocalEval', localGetFieldOrDefault(initEval, 'fdLocal', []));
  finalFdLocal = localGetFieldOrDefault(finalEval, 'fdLocalEval', localGetFieldOrDefault(finalEval, 'fdLocal', []));
  numSatCase = localInferNumSatForCase(initLocalDoa, finalLocalDoa, [], caseInfo);
  numFrameCase = localInferNumFrameForCase(initLocalDoa, finalLocalDoa, truthLocalDoaCase, caseInfo);
  refFrameIdxCase = min(max(sceneSeq.refFrameIdx, 1), numFrameCase);
  truthFdLocalSeriesMat = localBuildTruthFdLocalSeries(sceneSeq, truth, numSatCase, numFrameCase);
  for iSat = 1:numSatCase
    [initRefErr, initRmsErr, ~] = localCalcLocalDoaErrDeg(initLocalDoa, truthLocalDoaCase, iSat, refFrameIdxCase);
    [finalRefErr, finalRmsErr, finalMaxErr] = localCalcLocalDoaErrDeg(finalLocalDoa, truthLocalDoaCase, iSat, refFrameIdxCase);
    [~, truthRmsErr, ~] = localCalcLocalDoaErrDeg(truthLocalDoaCase, truthLocalDoaCase, iSat, refFrameIdxCase);
    [initAzRefErrDeg, initElRefErrDeg] = localGetLocalDoaComponentErr(initLocalDoa, truthLocalDoaCase, iSat, refFrameIdxCase);
    [finalAzRefErrDeg, finalElRefErrDeg] = localGetLocalDoaComponentErr(finalLocalDoa, truthLocalDoaCase, iSat, refFrameIdxCase);
    truthFdLocalSeries = reshape(truthFdLocalSeriesMat(iSat, :), 1, []);
    initFdLocalRefErrHz = localGetMatrixElem(initFdLocal, iSat, refFrameIdxCase, NaN) - localGetVectorElem(truthFdLocalSeries, refFrameIdxCase, NaN);
    finalFdLocalRefErrHz = localGetMatrixElem(finalFdLocal, iSat, refFrameIdxCase, NaN) - localGetVectorElem(truthFdLocalSeries, refFrameIdxCase, NaN);
    if ~isempty(finalFdLocal) && size(finalFdLocal, 1) >= iSat
      finalFdLocalRow = reshape(finalFdLocal(iSat, 1:min(size(finalFdLocal, 2), numFrameCase)), 1, []);
      truthFdLocalRow = truthFdLocalSeries(1:numel(finalFdLocalRow));
      fdLocalRmseHz = sqrt(mean((finalFdLocalRow - truthFdLocalRow).^2, 'omitnan'));
    else
      fdLocalRmseHz = NaN;
    end
    infoCell{end+1} = struct( ...
      'displayName', caseInfo.displayName, ...
      'satMode', caseInfo.satMode, ...
      'dynamicMode', caseInfo.dynamicMode, ...
      'localSatIdx', iSat, ...
      'initLocalDirRefErrDeg', initRefErr, ...
      'finalLocalDirRefErrDeg', finalRefErr, ...
      'rmsInitLocalDirErrDeg', initRmsErr, ...
      'rmsFinalLocalDirErrDeg', finalRmsErr, ...
      'rmsTruthLocalDirSelfErrDeg', truthRmsErr, ...
      'maxFinalLocalDirErrDeg', finalMaxErr, ...
      'initAzRefErrDeg', initAzRefErrDeg, ...
      'initElRefErrDeg', initElRefErrDeg, ...
      'finalAzRefErrDeg', finalAzRefErrDeg, ...
      'finalElRefErrDeg', finalElRefErrDeg, ...
      'initFdLocalRefErrHz', initFdLocalRefErrHz, ...
      'finalFdLocalRefErrHz', finalFdLocalRefErrHz, ...
      'fdLocalRmseHz', fdLocalRmseHz ...
      );
  end
end
if isempty(infoCell)
  dynLocalStateTable = table();
else
  dynLocalStateTable = struct2table([infoCell{:}], 'AsArray', true);
end
end



function dynTruthStateTable = localBuildDynamicTruthStateTable(caseResult)
%LOCALBUILDDYNAMICTRUTHTABLE Summarize truth-point internal mapping errors.

caseDyn = localSelectDynamicCase(caseResult);
infoCell = {};
for iCase = 1:numel(caseDyn)
  caseInfo = caseDyn(iCase);
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseInfo.estResult, 'aux', struct()), 'debug', struct());
  truthState = localGetFieldOrDefault(debugAux, 'truthState', struct());
  truthTable = localGetFieldOrDefault(truthState, 'summaryTable', table());
  if ~istable(truthTable) || isempty(truthTable)
    continue;
  end
  for iRow = 1:height(truthTable)
    infoCell{end + 1, 1} = struct( ...
      'displayName', caseInfo.displayName, ...
      'satMode', caseInfo.satMode, ...
      'dynamicMode', caseInfo.dynamicMode, ...
      'localSatIdx', localTableElemOrDefault(truthTable, 'localSatIdx', iRow, NaN), ...
      'globalSatIdx', localTableElemOrDefault(truthTable, 'globalSatIdx', iRow, NaN), ...
      'truthLocalDirRefErrDeg', localTableElemOrDefault(truthTable, 'truthLocalDirRefErrDeg', iRow, NaN), ...
      'truthRmsLocalDirErrDeg', localTableElemOrDefault(truthTable, 'truthRmsLocalDirErrDeg', iRow, NaN), ...
      'truthMaxLocalDirErrDeg', localTableElemOrDefault(truthTable, 'truthMaxLocalDirErrDeg', iRow, NaN), ...
      'truthAzRefErrDeg', localTableElemOrDefault(truthTable, 'truthAzRefErrDeg', iRow, NaN), ...
      'truthElRefErrDeg', localTableElemOrDefault(truthTable, 'truthElRefErrDeg', iRow, NaN), ...
      'truthFdLocalRefErrHz', localTableElemOrDefault(truthTable, 'truthFdLocalRefErrHz', iRow, NaN), ...
      'truthFdLocalRmseHz', localTableElemOrDefault(truthTable, 'truthFdLocalRmseHz', iRow, NaN)); %#ok<AGROW>
  end
end
if isempty(infoCell)
  dynTruthStateTable = table();
else
  dynTruthStateTable = struct2table([infoCell{:}], 'AsArray', true);
end
end



function dynSatOrderTable = localBuildDynamicSatOrderTable(caseResult)
%LOCALBUILDDYNAMICSATORDERTABLE Summarize satellite-order diagnostics.

caseDyn = localSelectDynamicCase(caseResult);
infoCell = {};
for iCase = 1:numel(caseDyn)
  caseInfo = caseDyn(iCase);
  debugAux = localGetFieldOrDefault(localGetFieldOrDefault(caseInfo.estResult, 'aux', struct()), 'debug', struct());
  satOrder = localGetFieldOrDefault(debugAux, 'satOrder', table());
  if ~istable(satOrder) || isempty(satOrder)
    continue;
  end
  for iRow = 1:height(satOrder)
    infoCell{end + 1, 1} = struct( ...
      'displayName', caseInfo.displayName, ...
      'satMode', caseInfo.satMode, ...
      'dynamicMode', caseInfo.dynamicMode, ...
      'localSatIdx', localTableElemOrDefault(satOrder, 'localSatIdx', iRow, NaN), ...
      'globalSatIdx', localTableElemOrDefault(satOrder, 'globalSatIdx', iRow, NaN), ...
      'refWeight', localTableElemOrDefault(satOrder, 'refWeight', iRow, NaN), ...
      'numValidFrame', localTableElemOrDefault(satOrder, 'numValidFrame', iRow, NaN), ...
      'isRefSat', logical(localTableElemOrDefault(satOrder, 'isRefSat', iRow, false))); %#ok<AGROW>
  end
end
if isempty(infoCell)
  dynSatOrderTable = table();
else
  dynSatOrderTable = struct2table([infoCell{:}], 'AsArray', true);
end
end


function truthLocalDoa = localExtractSceneLocalDoa(sceneSeq)
%LOCALEXTRACTSCENELOCALDOA Extract 2xNsxNf local DoA tensor from sceneSeq.
%
% Prefer the consolidated sceneSeq.localDoa field, which is the same source
% used by the steering-drift diagnostic. Fall back to per-frame sceneCell
% fields only when the consolidated tensor is unavailable.

truthLocalDoa = [];
if isfield(sceneSeq, 'localDoa') && ~isempty(sceneSeq.localDoa)
  localDoa = sceneSeq.localDoa;
  if ndims(localDoa) == 3
    truthLocalDoa = localDoa;
    return;
  elseif ndims(localDoa) == 4
    truthLocalDoa = localDoa(:, :, 1, :);
    return;
  end
end

numFrame = numel(sceneSeq.sceneCell);
numSat = localGetFieldOrDefault(sceneSeq, 'numSat', 0);
truthLocalDoa = nan(2, numSat, numFrame);
for iFrame = 1:numFrame
  sceneNow = sceneSeq.sceneCell{iFrame};
  localDoa = localGetFieldOrDefault(sceneNow, 'localDoa', []);
  if isempty(localDoa)
    localDoa = localGetFieldOrDefault(localGetFieldOrDefault(sceneNow, 'ref', struct()), 'localDoa', []);
  end
  if isempty(localDoa)
    localDoa = localGetFieldOrDefault(localGetFieldOrDefault(sceneNow, 'ref', struct()), 'doa', []);
  end
  if isempty(localDoa)
    continue;
  end
  if ndims(localDoa) == 3
    truthLocalDoa(:, 1:min(numSat, size(localDoa, 2)), iFrame) = localDoa(:, 1:min(numSat, size(localDoa, 2)), 1);
  elseif ismatrix(localDoa)
    truthLocalDoa(:, 1:min(numSat, size(localDoa, 2)), iFrame) = localDoa(:, 1:min(numSat, size(localDoa, 2)));
  end
end
end



function truthFdLocalSeriesMat = localBuildTruthFdLocalSeries(sceneSeq, truth, numSatCase, numFrameCase)
%LOCALBUILDTRUTHFDLOCALSERIES Build per-satellite truth fdLocal matrix.
%
% Prefer truth.fdSatSeries, which is already aligned with the selected
% satellite set and reference-star parameterization used in this script.

truthFdLocalSeriesMat = nan(numSatCase, numFrameCase);
fdTruth = localGetFieldOrDefault(truth, 'fdSatSeries', []);
if ~isempty(fdTruth)
  ns = min(numSatCase, size(fdTruth, 1));
  nf = min(numFrameCase, size(fdTruth, 2));
  truthFdLocalSeriesMat(1:ns, 1:nf) = fdTruth(1:ns, 1:nf);
  return;
end

for iFrame = 1:min(numFrameCase, numel(sceneSeq.sceneCell))
  sceneNow = sceneSeq.sceneCell{iFrame};
  fdVec = localGetFieldOrDefault(sceneNow, 'fd', []);
  if isempty(fdVec)
    fdVec = localGetFieldOrDefault(localGetFieldOrDefault(sceneNow, 'ref', struct()), 'fd', []);
  end
  for iSat = 1:numSatCase
    truthFdLocalSeriesMat(iSat, iFrame) = localGetMatrixElem(fdVec, iSat, 1, NaN);
  end
end
end



function localDoaCase = localAlignLocalDoaForCase(localDoaEval, truthLocalDoa, caseInfo)
%LOCALALIGNLOCALDOAFORCASE Align one local DoA tensor to current case size.

numSatCase = localInferNumSatForCase(localDoaEval, [], [], caseInfo);
numFrameCase = localInferNumFrameForCase(localDoaEval, truthLocalDoa, [], caseInfo);
localDoaCase = nan(2, numSatCase, numFrameCase);
source = localDoaEval;
if isempty(source)
  source = truthLocalDoa;
end
if isempty(source)
  return;
end
numSatSrc = min(numSatCase, size(source, 2));
if ndims(source) < 3
  source = repmat(source(:, 1:numSatSrc), 1, 1, numFrameCase);
end
numFrameSrc = min(numFrameCase, size(source, 3));
localDoaCase(:, 1:numSatSrc, 1:numFrameSrc) = source(:, 1:numSatSrc, 1:numFrameSrc);
end


function numSatCase = localInferNumSatForCase(varargin)
%LOCALINFERNUMSATFORCASE Infer satellite count from available tensors.

numSatCase = 1;
for iArg = 1:nargin
  value = varargin{iArg};
  if isnumeric(value) && ~isempty(value)
    numSatCase = max(numSatCase, size(value, 2));
  elseif isstruct(value) && isfield(value, 'satMode')
    if strcmp(string(value.satMode), "multi")
      numSatCase = max(numSatCase, 2);
    else
      numSatCase = max(numSatCase, 1);
    end
  end
end
end


function numFrameCase = localInferNumFrameForCase(varargin)
%LOCALINFERNUMFRAMEFORCASE Infer frame count from available tensors.

numFrameCase = 1;
for iArg = 1:nargin
  value = varargin{iArg};
  if isnumeric(value) && ndims(value) >= 3 && ~isempty(value)
    numFrameCase = max(numFrameCase, size(value, 3));
  end
end
end


function [refErrDeg, rmsErrDeg, maxErrDeg] = localCalcLocalDoaErrDeg(localDoaEval, truthLocalDoa, satIdx, refFrameIdx)
%LOCALCALCLOCALDOAERRDEG Compute local-direction error via unit-vector angle.

if nargin < 4 || isempty(refFrameIdx)
  refFrameIdx = 1;
end
refErrDeg = NaN;
rmsErrDeg = NaN;
maxErrDeg = NaN;
if isempty(localDoaEval) || isempty(truthLocalDoa)
  return;
end
numFrame = min(size(localDoaEval, 3), size(truthLocalDoa, 3));
if size(localDoaEval, 2) < satIdx || size(truthLocalDoa, 2) < satIdx
  return;
end
errDeg = nan(1, numFrame);
for iFrame = 1:numFrame
  evalNow = localDoaEval(:, satIdx, iFrame);
  truthNow = truthLocalDoa(:, satIdx, iFrame);
  if any(~isfinite(evalNow)) || any(~isfinite(truthNow))
    continue;
  end
  uEval = localDoaToUnitVec(evalNow);
  uTruth = localDoaToUnitVec(truthNow);
  cosVal = max(min(dot(uEval, uTruth), 1), -1);
  errDeg(iFrame) = acosd(cosVal);
end
validErr = errDeg(isfinite(errDeg));
if isempty(validErr)
  return;
end
refFrameIdx = min(max(refFrameIdx, 1), numFrame);
refErrDeg = errDeg(refFrameIdx);
if ~isfinite(refErrDeg)
  refErrDeg = validErr(1);
end
rmsErrDeg = sqrt(mean(validErr .^ 2));
maxErrDeg = max(validErr);
end


function [azErrDeg, elErrDeg] = localGetLocalDoaComponentErr(localDoaEval, truthLocalDoa, satIdx, frameIdx)
%LOCALGETLOCALDOACOMPONENTERR Get wrapped az/el component error at one frame.

azErrDeg = NaN;
elErrDeg = NaN;
if isempty(localDoaEval) || isempty(truthLocalDoa)
  return;
end
if size(localDoaEval, 2) < satIdx || size(truthLocalDoa, 2) < satIdx
  return;
end
if size(localDoaEval, 3) < frameIdx || size(truthLocalDoa, 3) < frameIdx
  return;
end
evalNow = localDoaEval(:, satIdx, frameIdx);
truthNow = truthLocalDoa(:, satIdx, frameIdx);
if any(~isfinite(evalNow)) || any(~isfinite(truthNow))
  return;
end
azErrDeg = rad2deg(localWrapDoaDiff(evalNow(1) - truthNow(1)));
elErrDeg = rad2deg(evalNow(2) - truthNow(2));
end



function u = localDoaToUnitVec(localDoa)
%LOCALDOATOUNITVEC Convert [az; el] local DoA to unit vector.

az = localDoa(1);
el = localDoa(2);
ce = cos(el);
u = [ce * cos(az); ce * sin(az); sin(el)];
end



function caseDyn = localSelectDynamicCase(caseResult)
%LOCALSELECTDYNAMICCASE Select MF-CP cases only.

keepMask = false(1, numel(caseResult));
for iCase = 1:numel(caseResult)
  keepMask(iCase) = contains(caseResult(iCase).displayName, "MF-CP");
end
caseDyn = caseResult(keepMask);
end



function [truthDeltaFdRefHz, truthDeltaFdRateHzPerSec] = localBuildTruthDeltaFdFit(truth)
%LOCALBUILDTRUTHDELTAFDFIT Build per-satellite truth line fits.

numSat = size(truth.deltaFdSeries, 1);
truthDeltaFdRefHz = zeros(numSat, 1);
truthDeltaFdRateHzPerSec = zeros(numSat, 1);
for iSat = 1:numSat
  [truthDeltaFdRefHz(iSat), truthDeltaFdRateHzPerSec(iSat)] = ...
    localFitFdLine(truth.timeOffsetSec, truth.deltaFdSeries(iSat, :));
end
end



function value = localGetMatrixElem(mat, rowIdx, colIdx, defaultValue)
%LOCALGETMATRIXELEM Get one matrix element with bounds checking.
%
% Supports both matrices and vectors. When MAT is a vector, COLIDX is
% ignored and the helper falls back to linear indexing by ROWIDX so that it
% can be safely used on per-satellite series extracted from scene structs.

value = defaultValue;
if isempty(mat)
  return;
end
if isvector(mat)
  vec = mat(:);
  if rowIdx >= 1 && rowIdx <= numel(vec)
    value = vec(rowIdx);
  end
  return;
end
if rowIdx >= 1 && rowIdx <= size(mat, 1) && colIdx >= 1 && colIdx <= size(mat, 2)
  value = mat(rowIdx, colIdx);
end
end



function value = localGetVectorElem(vec, idx, defaultValue)
%LOCALGETVECTORELEM Get one vector element with bounds checking.

value = defaultValue;
if isempty(vec)
  return;
end
vec = vec(:);
if idx >= 1 && idx <= numel(vec)
  value = vec(idx);
end
end




function [fit0, fitRate] = localFitFdLine(timeOffsetSec, fdSeries)
%LOCALFITFDLINE Fit fdSeries ~= fit0 + fitRate * timeOffsetSec.

timeOffsetSec = reshape(timeOffsetSec, [], 1);
fdSeries = reshape(fdSeries, [], 1);
A = [ones(numel(timeOffsetSec), 1), timeOffsetSec];
coeff = A \ fdSeries;
fit0 = coeff(1);
fitRate = coeff(2);
end



function deltaRad = localWrapDoaDiff(deltaRad)
%LOCALWRAPDOADIFF Wrap angular difference to [-pi, pi].

deltaRad = mod(deltaRad + pi, 2 * pi) - pi;
end




function latlonEst = localGetLatlonEst(estResult)
%LOCALGETLATLONEST Extract lat/lon estimate from estimator result.

latlonEst = [NaN; NaN];
if isempty(estResult)
  return;
end

aux = localGetFieldOrDefault(estResult, 'aux', struct());
latlonAux = localGetFieldOrDefault(aux, 'latlonEst', []);
if isnumeric(latlonAux) && numel(latlonAux) >= 2
  latlonEst = latlonAux(1:2);
  latlonEst = latlonEst(:);
  return;
end

latlonDirect = localGetFieldOrDefault(estResult, 'latlonEst', []);
if isnumeric(latlonDirect) && numel(latlonDirect) >= 2
  latlonEst = latlonDirect(1:2);
  latlonEst = latlonEst(:);
  return;
end

doaType = localGetFieldOrDefault(estResult, 'doaType', '');
doaParamEst = localGetFieldOrDefault(estResult, 'doaParamEst', []);
if ischar(doaType) || isstring(doaType)
  if strcmpi(string(doaType), "latlon") && isnumeric(doaParamEst) && numel(doaParamEst) >= 2
    latlonEst = doaParamEst(1:2);
    latlonEst = latlonEst(:);
  end
end
end



function value = localTableElemOrDefault(dataTable, varName, rowIdx, defaultValue)
%LOCALTABLEELEMORDEFAULT Read one table element with default fallback.

value = defaultValue;

if nargin < 4
  defaultValue = [];
  value = defaultValue;
end

if isempty(dataTable) || ~istable(dataTable)
  return;
end

if ~ismember(varName, dataTable.Properties.VariableNames)
  return;
end

if rowIdx < 1 || rowIdx > height(dataTable)
  return;
end

rawValue = dataTable{rowIdx, varName};
if iscell(rawValue)
  if isempty(rawValue)
    return;
  end
  value = rawValue{1};
  return;
end

if isempty(rawValue)
  return;
end

if isscalar(rawValue)
  value = rawValue;
else
  value = rawValue(1);
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
