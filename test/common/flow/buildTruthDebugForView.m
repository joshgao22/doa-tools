function debugTruth = buildTruthDebugForView(view, truth)
%BUILDTRUTHDEBUGFORVIEW Build one truth-debug struct for a scene view.
% This helper centralizes the truth-state packing used by the dynamic
% transition flow and representative replay scripts, so the parfor workers
% do not depend on file-local glue hidden inside one caller.

arguments
  view (1, 1) struct
  truth (1, 1) struct
end

sceneRef = view.sceneRef;
sceneSeq = view.sceneSeq;

if ~isfield(truth, 'latlonTrueDeg') || isempty(truth.latlonTrueDeg)
  error('buildTruthDebugForView:MissingTruthDoa', ...
    'truth.latlonTrueDeg is required to build the debug-truth view.');
end
if ~isfield(truth, 'fdRefFit') || isempty(truth.fdRefFit)
  error('buildTruthDebugForView:MissingTruthFdRef', ...
    'truth.fdRefFit is required to build the debug-truth view.');
end
if ~isfield(truth, 'fdRateFit') || isempty(truth.fdRateFit)
  error('buildTruthDebugForView:MissingTruthFdRate', ...
    'truth.fdRateFit is required to build the debug-truth view.');
end

debugTruth = struct();
debugTruth.doaParam = truth.latlonTrueDeg(:);
debugTruth.fdRef = truth.fdRefFit;
debugTruth.fdRate = truth.fdRateFit;
debugTruth.localDoaTrue = localExtractSceneLocalDoa(sceneSeq);
debugTruth.fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth);
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
debugTruth.refSatIdxLocal = refSatIdxLocal;
if isfield(sceneRef, 'satIdx') && ~isempty(sceneRef.satIdx)
  debugTruth.satIdxGlobal = reshape(sceneRef.satIdx, [], 1);
else
  debugTruth.satIdxGlobal = reshape(truth.selectedSatIdxGlobal(1:sceneRef.numSat), [], 1);
end
if isfield(sceneRef, 'ref') && isstruct(sceneRef.ref) && isfield(sceneRef.ref, 'weight') ...
    && ~isempty(sceneRef.ref.weight)
  debugTruth.refWeight = reshape(sceneRef.ref.weight, [], 1);
else
  debugTruth.refWeight = nan(sceneRef.numSat, 1);
end
end


function fdLocalTrue = localBuildTruthFdLocalForView(sceneRef, truth)
numSat = sceneRef.numSat;
numFrame = size(truth.fdSatSeries, 2);
fdLocalTrue = nan(numSat, numFrame);
sceneSatIdx = reshape(localGetFieldOrDefault(sceneRef, 'satIdx', []), 1, []);
truthSatIdx = reshape(localGetFieldOrDefault(truth, 'selectedSatIdxGlobal', []), 1, []);
if ~isempty(sceneSatIdx) && ~isempty(truthSatIdx)
  [isFound, loc] = ismember(sceneSatIdx, truthSatIdx);
  if any(isFound)
    fdLocalTrue(isFound, :) = truth.fdSatSeries(loc(isFound), :);
  end
  return;
end
numCopy = min(numSat, size(truth.fdSatSeries, 1));
fdLocalTrue(1:numCopy, :) = truth.fdSatSeries(1:numCopy, :);
end


function localDoa = localExtractSceneLocalDoa(sceneSeq)
localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
elseif ndims(localDoaRaw) == 4
  localDoa = localDoaRaw(:, :, 1, :);
else
  error('buildTruthDebugForView:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isempty(dataStruct)
  return;
end
if isfield(dataStruct, fieldName)
  value = dataStruct.(fieldName);
end
end
