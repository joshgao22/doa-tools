function caseTruth = buildDoaDopplerCaseTruthFromScene(truth, sceneView)
%BUILDDOADOPPLERCASETRUTHFROMSCENE Build one compact case-truth struct.
% The returned truth uses the subset-scene reference satellite when it can
% be resolved, so single-satellite ablation cases can reuse the same
% summary and metric helpers as the full multi-satellite case.

arguments
  truth
  sceneView
end

sceneInfo = localResolveSceneStruct(sceneView);
subsetSatIdxGlobal = localResolveSceneSatIdxGlobal(sceneInfo);
if isempty(subsetSatIdxGlobal)
  subsetSatIdxGlobal = reshape(getDoaDopplerFieldOrDefault(truth, 'selectedSatIdxGlobal', []), [], 1);
end

refSatIdxGlobal = localResolveSceneRefSatIdxGlobal(sceneInfo);
if ~isfinite(refSatIdxGlobal) && numel(subsetSatIdxGlobal) == 1
  refSatIdxGlobal = subsetSatIdxGlobal(1);
end

subsetFdSatTrueHz = localSliceTruthSatValue(truth, 'fdSatTrueHz', subsetSatIdxGlobal);
subsetFdRateSatTrueHzPerSec = localSliceTruthSatValue(truth, 'fdRateSatTrueHzPerSec', subsetSatIdxGlobal);

fdRefTrueHz = localResolveSceneTruthFdRef(sceneInfo, truth);
if isfinite(refSatIdxGlobal)
  fdRefMatch = localResolveTruthSatValue(truth, 'fdSatTrueHz', refSatIdxGlobal);
  if isfinite(fdRefMatch)
    fdRefTrueHz = fdRefMatch;
  end
end

fdRateTrueHzPerSec = localResolveSceneTruthFdRate(sceneInfo, truth);
if isfinite(refSatIdxGlobal)
  fdRateMatch = localResolveTruthSatValue(truth, 'fdRateSatTrueHzPerSec', refSatIdxGlobal);
  if isfinite(fdRateMatch)
    fdRateTrueHzPerSec = fdRateMatch;
  end
end

caseTruth = localBuildCaseTruth( ...
  getDoaDopplerFieldOrDefault(truth, 'latlonTrueDeg', nan(2, 1)), ...
  fdRefTrueHz, fdRateTrueHzPerSec, refSatIdxGlobal, ...
  subsetSatIdxGlobal(:).', subsetFdSatTrueHz);
if ~isempty(subsetFdRateSatTrueHzPerSec)
  caseTruth.fdRateSatTrueHzPerSec = reshape(subsetFdRateSatTrueHzPerSec, [], 1);
end
end


function sceneInfo = localResolveSceneStruct(sceneView)
%LOCALRESOLVESCENESTRUCT Resolve one scene struct from scene or view input.

sceneInfo = sceneView;
if isstruct(sceneView) && isfield(sceneView, 'sceneRef') && isstruct(sceneView.sceneRef)
  sceneInfo = sceneView.sceneRef;
end
end


function subsetSatIdxGlobal = localResolveSceneSatIdxGlobal(sceneInfo)
%LOCALRESOLVESCENESATIDXGLOBAL Resolve subset-scene global satellite list.

subsetSatIdxGlobal = reshape(getDoaDopplerFieldOrDefault(sceneInfo, 'satIdx', []), [], 1);
end


function refSatIdxGlobal = localResolveSceneRefSatIdxGlobal(sceneInfo)
%LOCALRESOLVESCENEREFSATIDXGLOBAL Resolve one subset-scene reference index.

refSatIdxGlobal = NaN;
refInfo = getDoaDopplerFieldOrDefault(sceneInfo, 'ref', struct());
if ~isstruct(refInfo) || isempty(refInfo)
  return;
end

refSatIdxGlobal = getDoaDopplerFieldOrDefault(refInfo, 'satIdxGlobal', NaN);
if isfinite(refSatIdxGlobal)
  return;
end

sceneSatIdx = reshape(getDoaDopplerFieldOrDefault(sceneInfo, 'satIdx', []), [], 1);
refSatIdxLocal = getDoaDopplerFieldOrDefault(refInfo, 'satIdxLocal', NaN);
if isempty(sceneSatIdx) || ~isscalar(refSatIdxLocal) || ~isfinite(refSatIdxLocal)
  refSatIdxGlobal = NaN;
  return;
end
if refSatIdxLocal < 1 || refSatIdxLocal > numel(sceneSatIdx)
  refSatIdxGlobal = NaN;
  return;
end
refSatIdxGlobal = sceneSatIdx(refSatIdxLocal);
end


function fdRefTrueHz = localResolveSceneTruthFdRef(sceneInfo, truth)
%LOCALRESOLVESCENETRUTHFDREF Resolve subset-aware reference Doppler truth.

fdRefTrueHz = getDoaDopplerFieldOrDefault(truth, 'fdRefTrueHz', NaN);
if ~isfinite(fdRefTrueHz)
  fdRefTrueHz = getDoaDopplerFieldOrDefault(truth, 'fdRefFit', NaN);
end

truthSatIdx = reshape(getDoaDopplerFieldOrDefault(truth, 'selectedSatIdxGlobal', []), [], 1);
truthFdSat = reshape(getDoaDopplerFieldOrDefault(truth, 'fdSatTrueHz', []), [], 1);
if isempty(truthSatIdx) || isempty(truthFdSat)
  return;
end

refSatIdxGlobal = localResolveSceneRefSatIdxGlobal(sceneInfo);
if ~isfinite(refSatIdxGlobal)
  sceneSatIdx = localResolveSceneSatIdxGlobal(sceneInfo);
  if numel(sceneSatIdx) == 1
    refSatIdxGlobal = sceneSatIdx(1);
  end
end
if ~isfinite(refSatIdxGlobal)
  return;
end

matchIdx = find(truthSatIdx == refSatIdxGlobal, 1, 'first');
if ~isempty(matchIdx) && numel(truthFdSat) >= matchIdx && isfinite(truthFdSat(matchIdx))
  fdRefTrueHz = truthFdSat(matchIdx);
end
end


function fdRateTrueHzPerSec = localResolveSceneTruthFdRate(sceneInfo, truth)
%LOCALRESOLVESCENETRUTHFDRATE Resolve subset-aware reference rate truth.

fdRateTrueHzPerSec = getDoaDopplerFieldOrDefault(truth, 'fdRateTrueHzPerSec', NaN);
if ~isfinite(fdRateTrueHzPerSec)
  fdRateTrueHzPerSec = getDoaDopplerFieldOrDefault(truth, 'fdRateFit', NaN);
end

refSatIdxGlobal = localResolveSceneRefSatIdxGlobal(sceneInfo);
if ~isfinite(refSatIdxGlobal)
  sceneSatIdx = localResolveSceneSatIdxGlobal(sceneInfo);
  if numel(sceneSatIdx) == 1
    refSatIdxGlobal = sceneSatIdx(1);
  end
end
if ~isfinite(refSatIdxGlobal)
  return;
end

truthFdRateSat = localResolveTruthSatValue(truth, 'fdRateSatTrueHzPerSec', refSatIdxGlobal);
if isfinite(truthFdRateSat)
  fdRateTrueHzPerSec = truthFdRateSat;
  return;
end

truthRefSatIdxGlobal = getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN);
if isfinite(truthRefSatIdxGlobal) && refSatIdxGlobal ~= truthRefSatIdxGlobal
  fdRateTrueHzPerSec = NaN;
end
end


function satValue = localResolveTruthSatValue(truth, fieldName, satIdxGlobal)
%LOCALRESOLVETRUTHSATVALUE Resolve one per-satellite truth value.

satValue = NaN;
if ~isscalar(satIdxGlobal) || ~isfinite(satIdxGlobal)
  return;
end

truthSatIdx = reshape(getDoaDopplerFieldOrDefault(truth, 'selectedSatIdxGlobal', []), [], 1);
truthVal = reshape(getDoaDopplerFieldOrDefault(truth, fieldName, []), [], 1);
if isempty(truthSatIdx) || isempty(truthVal)
  return;
end

matchIdx = find(truthSatIdx == satIdxGlobal, 1, 'first');
if ~isempty(matchIdx) && numel(truthVal) >= matchIdx && isfinite(truthVal(matchIdx))
  satValue = truthVal(matchIdx);
end
end


function subsetVal = localSliceTruthSatValue(truth, fieldName, satIdxGlobalList)
%LOCALSLICETRUTHSATVALUE Slice one truth vector by global satellite index.

subsetVal = [];
satIdxGlobalList = reshape(satIdxGlobalList, [], 1);
if isempty(satIdxGlobalList)
  return;
end

subsetVal = nan(numel(satIdxGlobalList), 1);
for iSat = 1:numel(satIdxGlobalList)
  subsetVal(iSat) = localResolveTruthSatValue(truth, fieldName, satIdxGlobalList(iSat));
end
if all(~isfinite(subsetVal))
  subsetVal = [];
end
end


function caseTruth = localBuildCaseTruth(latlonTrueDeg, fdRefTrueHz, fdRateTrueHzPerSec, ...
  refSatIdxGlobal, selectedSatIdxGlobal, fdSatTrueHz)
%LOCALBUILDCASETRUTH Build one compact truth struct.

if nargin < 3 || isempty(fdRateTrueHzPerSec)
  fdRateTrueHzPerSec = NaN;
end
if nargin < 4 || isempty(refSatIdxGlobal)
  refSatIdxGlobal = NaN;
end
if nargin < 5 || isempty(selectedSatIdxGlobal)
  selectedSatIdxGlobal = [];
end
if nargin < 6 || isempty(fdSatTrueHz)
  fdSatTrueHz = [];
end

caseTruth = struct();
caseTruth.latlonTrueDeg = reshape(latlonTrueDeg, [], 1);
caseTruth.fdRefTrueHz = fdRefTrueHz;
caseTruth.fdRateTrueHzPerSec = fdRateTrueHzPerSec;
caseTruth.refSatIdxGlobal = refSatIdxGlobal;
caseTruth.selectedSatIdxGlobal = reshape(selectedSatIdxGlobal, 1, []);
caseTruth.fdSatTrueHz = reshape(fdSatTrueHz, [], 1);
end
