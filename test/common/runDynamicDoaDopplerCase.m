function caseInfo = runDynamicDoaDopplerCase(displayName, satMode, view, ...
  truth, pilotWave, carrierFreq, sampleRate, fdRange, fdRateRange, ...
  verbose, dynBaseOpt, isKnownRate, debugTruth, initParamOverride)
%RUNDYNAMICDOADOPPLERCASE Run one multi-frame continuous-phase case.

if nargin < 13 || isempty(debugTruth)
  debugTruth = [];
end
if nargin < 14
  initParamOverride = [];
end

dynOpt = dynBaseOpt;
if isKnownRate
  dynOpt.fdRateMode = 'known';
  dynOpt.fdRateKnown = truth.fdRateFit;
  fdRateRangeUse = [];
  dynamicMode = "cp-known";
else
  dynOpt.fdRateMode = 'unknown';
  fdRateRangeUse = fdRateRange;
  dynamicMode = "cp-unknown";
end

if ~isempty(debugTruth)
  dynOpt.debugTruth = debugTruth;
end

[estResult, multiStartDiag] = localRunDynamicEstimatorWithCandidates( ...
  view, truth, pilotWave, carrierFreq, sampleRate, fdRange, fdRateRangeUse, ...
  verbose, dynOpt, initParamOverride);

caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "multi", ...
  "doa-doppler", dynamicMode, estResult);
if multiStartDiag.isEnabled
  aux = localGetFieldOrDefault(caseInfo.estResult, 'aux', struct());
  aux.multiStart = multiStartDiag;
  caseInfo.estResult.aux = aux;
end
end


function [estResultBest, multiStartDiag] = localRunDynamicEstimatorWithCandidates( ...
  view, truth, pilotWave, carrierFreq, sampleRate, fdRange, fdRateRange, ...
  verbose, dynOptBase, initParamOverride)
%LOCALRUNDYNAMICESTIMATORWITHCANDIDATES Evaluate one or multiple init seeds.

multiStartDiag = struct();
multiStartDiag.isEnabled = false;
multiStartDiag.selectedIdx = 1;
multiStartDiag.selectedTag = "single";
multiStartDiag.tagList = strings(0, 1);
multiStartDiag.fvalList = [];
multiStartDiag.isResolvedList = [];
multiStartDiag.iterationList = [];
multiStartDiag.objImproveList = [];
multiStartDiag.initObjList = [];
multiStartDiag.finalObjList = [];
multiStartDiag.moveNormList = [];
multiStartDiag.angleErrDegList = [];
multiStartDiag.fdRefErrHzList = [];
multiStartDiag.fdRateErrHzPerSecList = [];
multiStartDiag.solveVariantList = strings(0, 1);
multiStartDiag.innerCandidateObjectiveList = {};
multiStartDiag.innerCandidateVariantList = {};

if isempty(initParamOverride) || (~iscell(initParamOverride) && ~isstruct(initParamOverride))
  [estResultBest, ~, ~] = estimatorDoaDopplerMlePilotMfOpt( ...
    view.sceneSeq, view.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
    view.doaGrid, fdRange, fdRateRange, initParamOverride, verbose, dynOptBase);
  return;
end

if iscell(initParamOverride)
  candidateCell = reshape(initParamOverride, 1, []);
elseif isstruct(initParamOverride)
  candidateCell = num2cell(reshape(initParamOverride, 1, []));
else
  candidateCell = {initParamOverride};
end

numCandidate = numel(candidateCell);
resultCell = cell(1, numCandidate);
tagList = strings(numCandidate, 1);
fvalList = nan(numCandidate, 1);
isResolvedList = false(numCandidate, 1);
iterationList = nan(numCandidate, 1);
objImproveList = nan(numCandidate, 1);
initObjList = nan(numCandidate, 1);
finalObjList = nan(numCandidate, 1);
moveNormList = nan(numCandidate, 1);
angleErrDegList = nan(numCandidate, 1);
fdRefErrHzList = nan(numCandidate, 1);
fdRateErrHzPerSecList = nan(numCandidate, 1);
solveVariantList = strings(numCandidate, 1);
innerCandidateObjectiveList = cell(numCandidate, 1);
innerCandidateVariantList = cell(numCandidate, 1);

for iCand = 1:numCandidate
  candidateItem = candidateCell{iCand};
  dynOptCur = dynOptBase;
  initParamCur = [];
  tagList(iCand) = sprintf("cand%d", iCand);

  if isnumeric(candidateItem)
    initParamCur = candidateItem;
  elseif isstruct(candidateItem)
    if isfield(candidateItem, 'initParam')
      initParamCur = candidateItem.initParam;
    end
    if isfield(candidateItem, 'initDoaParam') && ~isempty(candidateItem.initDoaParam)
      dynOptCur.initDoaParam = candidateItem.initDoaParam;
    end
    if isfield(candidateItem, 'initDoaHalfWidth') && ~isempty(candidateItem.initDoaHalfWidth)
      dynOptCur.initDoaHalfWidth = candidateItem.initDoaHalfWidth;
    end
    if isfield(candidateItem, 'tag') && ~isempty(candidateItem.tag)
      tagList(iCand) = string(candidateItem.tag);
    end
  else
    error('runDynamicDoaDopplerCase:InvalidInitCandidate', ...
      'Each init candidate must be numeric or a struct.');
  end

  [resultCell{iCand}, ~, ~] = estimatorDoaDopplerMlePilotMfOpt( ...
    view.sceneSeq, view.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
    view.doaGrid, fdRange, fdRateRange, initParamCur, verbose, dynOptCur);

  resultCur = resultCell{iCand};
  fvalList(iCand) = localGetFieldOrDefault(resultCur, 'fval', NaN);
  isResolvedList(iCand) = logical(localGetFieldOrDefault(resultCur, 'isResolved', false));
  optimInfoCur = localGetFieldOrDefault(resultCur, 'optimInfo', struct());
  iterationList(iCand) = localGetFieldOrDefault(optimInfoCur, 'iterations', NaN);
  auxCur = localGetFieldOrDefault(resultCur, 'aux', struct());
  debugCur = localGetFieldOrDefault(auxCur, 'debug', struct());
  initEvalCur = localGetFieldOrDefault(debugCur, 'initEval', struct());
  finalEvalCur = localGetFieldOrDefault(debugCur, 'finalEval', struct());
  initObjList(iCand) = localGetFieldOrDefault(initEvalCur, 'obj', NaN);
  finalObjList(iCand) = localGetFieldOrDefault(finalEvalCur, 'obj', NaN);
  objImproveList(iCand) = initObjList(iCand) - finalObjList(iCand);
  initParamEst = reshape(localGetFieldOrDefault(resultCur, 'initParam', []), [], 1);
  optVarEst = reshape(localGetFieldOrDefault(resultCur, 'optVarEst', []), [], 1);
  if ~isempty(initParamEst) && numel(initParamEst) == numel(optVarEst)
    moveNormList(iCand) = norm(optVarEst - initParamEst);
  end

  [angleErrDegList(iCand), fdRefErrHzList(iCand), fdRateErrHzPerSecList(iCand)] = ...
    localComputeDynamicCandidateMetric(resultCur, truth);
  solveVariantList(iCand) = string(localGetFieldOrDefault(optimInfoCur, 'solveVariant', ""));
  innerCandidateObjectiveList{iCand} = localGetFieldOrDefault(optimInfoCur, 'candidateObjective', []);
  innerCandidateVariantList{iCand} = localGetFieldOrDefault(optimInfoCur, 'candidateVariant', strings(0, 1));
end

selectedIdx = localSelectBestDynamicCandidate(fvalList, isResolvedList, tagList, iterationList, objImproveList, moveNormList);
estResultBest = resultCell{selectedIdx};

multiStartDiag.isEnabled = true;
multiStartDiag.selectedIdx = selectedIdx;
multiStartDiag.selectedTag = tagList(selectedIdx);
multiStartDiag.tagList = tagList;
multiStartDiag.fvalList = fvalList;
multiStartDiag.isResolvedList = isResolvedList;
multiStartDiag.iterationList = iterationList;
multiStartDiag.objImproveList = objImproveList;
multiStartDiag.initObjList = initObjList;
multiStartDiag.finalObjList = finalObjList;
multiStartDiag.moveNormList = moveNormList;
multiStartDiag.angleErrDegList = angleErrDegList;
multiStartDiag.fdRefErrHzList = fdRefErrHzList;
multiStartDiag.fdRateErrHzPerSecList = fdRateErrHzPerSecList;
multiStartDiag.solveVariantList = solveVariantList;
multiStartDiag.innerCandidateObjectiveList = innerCandidateObjectiveList;
multiStartDiag.innerCandidateVariantList = innerCandidateVariantList;
end


function selectedIdx = localSelectBestDynamicCandidate(fvalList, isResolvedList, tagList, iterationList, objImproveList, moveNormList)
%LOCALSELECTBESTDYNAMICCANDIDATE Prefer resolved candidates with a real release.
% When CP-U is launched from multiple warm starts, a CP-K continuation can
% sometimes return immediately with one iteration and effectively copy the
% CP-K solution. If a static-seed release reaches nearly the same objective,
% prefer the candidate that actually moved and improved.

selectedIdx = 1;
validResolved = find(isResolvedList & isfinite(fvalList));
if isempty(validResolved)
  validAny = find(isfinite(fvalList));
  if ~isempty(validAny)
    [~, idxRel] = min(fvalList(validAny));
    selectedIdx = validAny(idxRel);
  end
  return;
end

bestFval = min(fvalList(validResolved));
% Allow a modest objective tolerance so CP-U can prefer a real release from
% the static seed when it lands in essentially the same basin as the CP-K
% continuation. Tiny pure-fval tolerances repeatedly selected the one-step
% CP-K copy and masked a usable static-seed release.
objTol = max([1e-3, 1e-7 * max(abs(bestFval), 1), 250]);
nearBest = validResolved(abs(fvalList(validResolved) - bestFval) <= objTol);
if numel(nearBest) == 1
  selectedIdx = nearBest(1);
  return;
end

% Prefer candidates that actually improved or moved away from the warm seed.
scoreImprove = objImproveList(nearBest);
scoreImprove(~isfinite(scoreImprove)) = -inf;
scoreMove = moveNormList(nearBest);
scoreMove(~isfinite(scoreMove)) = -inf;
scoreIter = iterationList(nearBest);
scoreIter(~isfinite(scoreIter)) = -inf;
preferStatic = contains(lower(string(tagList(nearBest))), "static");

scoreMat = [scoreImprove(:), scoreMove(:), scoreIter(:), double(preferStatic(:))];
[~, idxRel] = sortrows(scoreMat, [-1, -2, -3, -4]);
selectedIdx = nearBest(idxRel(1));
end



function [angleErrDeg, fdRefErrHz, fdRateErrHzPerSec] = localComputeDynamicCandidateMetric(estResult, truth)
%LOCALCOMPUTEDYNAMICCANDIDATEMETRIC Compute compact candidate errors.

angleErrDeg = NaN;
fdRefErrHz = NaN;
fdRateErrHzPerSec = NaN;

if isempty(estResult) || ~isstruct(estResult)
  return;
end

truthLatlon = reshape(localGetFieldOrDefault(truth, 'latlonTrueDeg', []), [], 1);
estLatlon = reshape(localGetFieldOrDefault(estResult, 'doaParamEst', []), [], 1);
if numel(estLatlon) >= 2 && numel(truthLatlon) >= 2
  angleErrDeg = calcLatlonAngleError(estLatlon(1:2), truthLatlon(1:2));
end

fdRefEst = localGetFieldOrDefault(estResult, 'fdRefEst', NaN);
fdRefTrueHz = localGetFieldOrDefault(truth, 'fdRefTrueHz', localGetFieldOrDefault(truth, 'fdRefFit', NaN));
fdRefErrHz = localResolveScalarMetric(fdRefEst, fdRefTrueHz);

fdRateEst = localGetFieldOrDefault(estResult, 'fdRateEst', NaN);
fdRateTrueHzPerSec = localGetFieldOrDefault(truth, 'fdRateTrueHzPerSec', localGetFieldOrDefault(truth, 'fdRateFit', NaN));
fdRateErrHzPerSec = localResolveScalarMetric(fdRateEst, fdRateTrueHzPerSec);
end

function scalarErr = localResolveScalarMetric(estValue, truthValue)
%LOCALRESOLVESCALARMETRIC Convert possibly vector-valued fields into one error.

scalarErr = NaN;
if isempty(estValue) || isempty(truthValue)
  return;
end

estValue = estValue(1);
truthValue = truthValue(1);
if isfinite(estValue) && isfinite(truthValue)
  scalarErr = estValue - truthValue;
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

if isobject(dataStruct) && isprop(dataStruct, fieldName)
  fieldValue = dataStruct.(fieldName);
end
end
