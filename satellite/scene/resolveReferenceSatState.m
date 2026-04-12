function [refState, refSatIdxLocal] = resolveReferenceSatState(scene, satPosEci, satVelEci)
%RESOLVEREFERENCESATSTATE Resolve one reference-satellite state consistently.
% The project-level parameterization assumes that DoA and Doppler are
% defined with respect to one selected reference satellite, rather than an
% arbitrary centroid reference. This helper resolves that reference
% satellite from the scene metadata and returns a one-hot reference state
% aligned with the selected local satellite index.
%
%Syntax:
%  [refState, refSatIdxLocal] = resolveReferenceSatState(scene, satPosEci, satVelEci)
%
%Inputs:
%  scene           - single-frame scene structure
%  satPosEci       - 3xNs satellite ECI position matrix
%  satVelEci       - 3xNs satellite ECI velocity matrix
%
%Outputs:
%  refState        - reference state struct with fields:
%                      .type
%                      .weight
%                      .posEci
%                      .velEci
%                      .satIdxLocal
%                      .satIdxGlobal
%                      .source
%  refSatIdxLocal  - resolved local reference-satellite index
%
%Resolution rule:
%  1) Prefer explicit scene.ref.satIdxLocal / satIdxGlobal tags.
%  2) Otherwise prefer a one-hot or dominant scene.ref.weight entry.
%  3) Otherwise match scene.ref.posEci to one satellite state.
%  4) Fall back to the single available satellite, or satellite 1.
%
%See also:
%  selectSatScene, estimatorDoaDopplerMlePilotSfOpt,
%  estimatorDoaDopplerMlePilotMfOpt

numSat = size(satPosEci, 2);
if size(satVelEci, 2) ~= numSat
  error('resolveReferenceSatState:SatCountMismatch', ...
    'satPosEci and satVelEci must have the same satellite count.');
end
if numSat < 1
  error('resolveReferenceSatState:EmptySatelliteSet', ...
    'At least one satellite is required.');
end

[refSatIdxLocal, sourceTag] = localResolveReferenceSatIdx(scene, satPosEci);

refWeight = zeros(numSat, 1);
refWeight(refSatIdxLocal) = 1;

refState = struct();
refState.type = 'selectedSat';
refState.weight = refWeight;
refState.posEci = reshape(satPosEci(:, refSatIdxLocal), 3, 1);
refState.velEci = reshape(satVelEci(:, refSatIdxLocal), 3, 1);
refState.satIdxLocal = refSatIdxLocal;
refState.satIdxGlobal = localResolveGlobalSatIdx(scene, refSatIdxLocal);
refState.source = sourceTag;
end


function [refSatIdxLocal, sourceTag] = localResolveReferenceSatIdx(scene, satPosEci)
%LOCALRESOLVEREFERENCESATIDX Resolve one local reference-satellite index.

numSat = size(satPosEci, 2);
refSatIdxLocal = [];
sourceTag = "fallback";

if isfield(scene, 'ref') && isstruct(scene.ref)
  if isfield(scene.ref, 'satIdxLocal') && isscalar(scene.ref.satIdxLocal) && ...
      isnumeric(scene.ref.satIdxLocal) && isfinite(scene.ref.satIdxLocal) && ...
      scene.ref.satIdxLocal >= 1 && scene.ref.satIdxLocal <= numSat && ...
      mod(scene.ref.satIdxLocal, 1) == 0
    refSatIdxLocal = scene.ref.satIdxLocal;
    sourceTag = "satIdxLocal";
  end

  if isempty(refSatIdxLocal) && isfield(scene.ref, 'satIdxGlobal') && ...
      isscalar(scene.ref.satIdxGlobal) && isnumeric(scene.ref.satIdxGlobal) && ...
      isfinite(scene.ref.satIdxGlobal) && mod(scene.ref.satIdxGlobal, 1) == 0 && ...
      isfield(scene, 'satIdx') && ~isempty(scene.satIdx)
    satIdxVec = reshape(scene.satIdx, [], 1);
    matchIdx = find(satIdxVec == scene.ref.satIdxGlobal, 1, 'first');
    if ~isempty(matchIdx)
      refSatIdxLocal = matchIdx;
      sourceTag = "satIdxGlobal";
    end
  end

  if isempty(refSatIdxLocal) && isfield(scene.ref, 'weight') && ~isempty(scene.ref.weight)
    refWeight = scene.ref.weight(:);
    if numel(refWeight) == numSat && all(isfinite(refWeight))
      strongIdx = find(refWeight > 0.5, 1, 'first');
      if ~isempty(strongIdx)
        refSatIdxLocal = strongIdx;
        sourceTag = "weight";
      else
        [weightMax, weightIdx] = max(refWeight);
        if isfinite(weightMax) && weightMax > 0
          refSatIdxLocal = weightIdx;
          sourceTag = "weightMax";
        end
      end
    end
  end

  if isempty(refSatIdxLocal) && isfield(scene.ref, 'posEci') && ...
      isnumeric(scene.ref.posEci) && isequal(size(scene.ref.posEci), [3, 1])
    refPos = reshape(scene.ref.posEci, 3, 1);
    errVec = vecnorm(satPosEci - refPos, 2, 1);
    [errMin, idxMin] = min(errVec);
    posScale = max(vecnorm(satPosEci, 2, 1));
    if isempty(posScale) || ~isfinite(posScale) || posScale <= 0
      posScale = 1;
    end
    if isfinite(errMin) && errMin <= 1e-9 * posScale
      refSatIdxLocal = idxMin;
      sourceTag = "position";
    end
  end
end

if isempty(refSatIdxLocal)
  if numSat == 1
    refSatIdxLocal = 1;
    sourceTag = "singleSat";
  else
    refSatIdxLocal = 1;
    sourceTag = "defaultFirstSat";
  end
end
end


function satIdxGlobal = localResolveGlobalSatIdx(scene, refSatIdxLocal)
%LOCALRESOLVEGLOBALSATIDX Resolve the global satellite index when available.

satIdxGlobal = NaN;

if isstruct(scene) && isfield(scene, 'satIdx') && ~isempty(scene.satIdx)
  satIdxVec = reshape(scene.satIdx, [], 1);
  if numel(satIdxVec) >= refSatIdxLocal
    satIdxGlobal = satIdxVec(refSatIdxLocal);
    return;
  end
end

if isstruct(scene) && isfield(scene, 'ref') && isstruct(scene.ref) && ...
    isfield(scene.ref, 'satIdxLocal') && isfield(scene.ref, 'satIdxGlobal')
  refSatIdxLocalScene = scene.ref.satIdxLocal;
  refSatIdxGlobalScene = scene.ref.satIdxGlobal;
  if isscalar(refSatIdxLocalScene) && isnumeric(refSatIdxLocalScene) && ...
      isfinite(refSatIdxLocalScene) && mod(refSatIdxLocalScene, 1) == 0 && ...
      isscalar(refSatIdxGlobalScene) && isnumeric(refSatIdxGlobalScene) && ...
      isfinite(refSatIdxGlobalScene) && mod(refSatIdxGlobalScene, 1) == 0 && ...
      round(refSatIdxLocalScene) == refSatIdxLocal
    satIdxGlobal = refSatIdxGlobalScene;
  end
end
end
