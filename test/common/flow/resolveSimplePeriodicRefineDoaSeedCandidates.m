function seedCandidateList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOpt, satMode, staticSeedCase, selectedSubsetCase)
%RESOLVESIMPLEPERIODICREFINEDOASEEDCANDIDATES Build periodic replay DoA seed candidates.
% The simple subset-periodic flow uses the subset branch to pick one fd
% tooth, but the periodic replay should not blindly trust either the static
% seed or the subset DoA. Instead, the flow can compare multiple frozen DoA
% seed basins inside the same narrow [fdRef, fdRate] box and only then run a
% very-small same-tooth DoA polish when warranted.
%
% Supported seed modes via flowOpt.periodicRefineDoaSeedMode:
%   - "staticSeed"      : only static DoA
%   - "selectedSubset"  : only selected-subset DoA
%   - "subsetWhenClose" : subset when it stays close to static, otherwise static
%   - "dualWhenMulti"   : for multi-sat replay, compare subset and static
%                          frozen replays; single-sat still uses subset only
%   - "dualSeed"        : always compare subset and static when both exist

arguments
  flowOpt (1,1) struct
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  staticSeedCase (1,1) struct
  selectedSubsetCase (1,1) struct
end

staticDoa = localGetCaseDoa(staticSeedCase);
subsetDoa = localGetCaseDoa(selectedSubsetCase);
doAFreeze = logical(localGetFieldOrDefault(flowOpt, 'periodicRefineFreezeDoa', true));
doAHalfWidthDeg = reshape(localGetFieldOrDefault(flowOpt, 'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8]), [], 1);
maxSubsetDriftDeg = localGetFieldOrDefault(flowOpt, 'periodicRefineMaxSubsetDoaDriftDeg', 0.003);
modeUse = lower(strtrim(char(string(localGetFieldOrDefault(flowOpt, ...
  'periodicRefineDoaSeedMode', "dualWhenMulti")))));

driftDeg = localCalcDoaDriftDeg(subsetDoa, staticDoa);
hasStatic = ~isempty(staticDoa) && all(isfinite(staticDoa));
hasSubset = ~isempty(subsetDoa) && all(isfinite(subsetDoa));

seedCandidateList = struct('doaInitParam', {}, 'seedSource', {}, ...
  'subsetDriftFromStaticDeg', {}, 'freezeDoa', {}, 'doaHalfWidthDeg', {}, ...
  'startTag', {});

switch modeUse
  case 'staticseed'
    seedCandidateList = localAppendIfValid(seedCandidateList, staticDoa, "static-seed", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-static-seed");
  case 'selectedsubset'
    seedCandidateList = localAppendIfValid(seedCandidateList, subsetDoa, "selected-subset", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-subset-seed");
  case 'subsetwhenclose'
    if hasSubset && (~hasStatic || (~isfinite(driftDeg)) || (driftDeg <= maxSubsetDriftDeg))
      seedCandidateList = localAppendIfValid(seedCandidateList, subsetDoa, "selected-subset", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-subset-seed");
    else
      seedCandidateList = localAppendIfValid(seedCandidateList, staticDoa, "static-seed", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-static-seed");
    end
  case 'dualwhenmulti'
    if satMode == "multi"
      seedCandidateList = localAppendDual(seedCandidateList, subsetDoa, staticDoa, driftDeg, doAFreeze, doAHalfWidthDeg);
    else
      seedCandidateList = localAppendIfValid(seedCandidateList, subsetDoa, "selected-subset", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-subset-seed");
    end
  case 'dualseed'
    seedCandidateList = localAppendDual(seedCandidateList, subsetDoa, staticDoa, driftDeg, doAFreeze, doAHalfWidthDeg);
  otherwise
    error('resolveSimplePeriodicRefineDoaSeedCandidates:InvalidSeedMode', ...
      'Unsupported periodicRefineDoaSeedMode: %s.', modeUse);
end

if isempty(seedCandidateList)
  if hasSubset
    seedCandidateList = localAppendIfValid(seedCandidateList, subsetDoa, "selected-subset", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-subset-seed");
  elseif hasStatic
    seedCandidateList = localAppendIfValid(seedCandidateList, staticDoa, "static-seed", driftDeg, doAFreeze, doAHalfWidthDeg, "periodic-static-seed");
  else
    seedCandidateList = localAppendIfValid(seedCandidateList, nan(2, 1), "missing", NaN, doAFreeze, doAHalfWidthDeg, "periodic-missing-seed");
  end
end
end

function seedCandidateList = localAppendDual(seedCandidateList, subsetDoa, staticDoa, driftDeg, freezeDoa, doaHalfWidthDeg)
seedCandidateList = localAppendIfValid(seedCandidateList, subsetDoa, "selected-subset", driftDeg, freezeDoa, doaHalfWidthDeg, "periodic-subset-seed");
if ~localDoaAlreadyPresent(seedCandidateList, staticDoa)
  seedCandidateList = localAppendIfValid(seedCandidateList, staticDoa, "static-seed", driftDeg, freezeDoa, doaHalfWidthDeg, "periodic-static-seed");
end
end

function seedCandidateList = localAppendIfValid(seedCandidateList, doaParam, seedSource, driftDeg, freezeDoa, doaHalfWidthDeg, startTag)
if isempty(doaParam) || any(~isfinite(doaParam(:)))
  return;
end
entry = struct();
entry.doaInitParam = reshape(doaParam(:), [], 1);
entry.seedSource = string(seedSource);
entry.subsetDriftFromStaticDeg = driftDeg;
entry.freezeDoa = logical(freezeDoa);
entry.doaHalfWidthDeg = reshape(doaHalfWidthDeg(:), [], 1);
entry.startTag = string(startTag);
seedCandidateList(end + 1, 1) = entry; %#ok<AGROW>
end

function tf = localDoaAlreadyPresent(seedCandidateList, doaParam)
tf = false;
if isempty(seedCandidateList) || isempty(doaParam) || any(~isfinite(doaParam(:)))
  return;
end
for iCand = 1:numel(seedCandidateList)
  refDoa = reshape(seedCandidateList(iCand).doaInitParam, [], 1);
  if numel(refDoa) ~= numel(doaParam)
    continue;
  end
  if all(abs(refDoa(:) - doaParam(:)) <= 1e-12)
    tf = true;
    return;
  end
end
end

function doaParam = localGetCaseDoa(caseUse)
doaParam = [];
if ~isstruct(caseUse) || ~isfield(caseUse, 'estResult') || isempty(caseUse.estResult)
  return;
end
estResult = caseUse.estResult;
if ~isfield(estResult, 'doaParamEst') || isempty(estResult.doaParamEst)
  return;
end
rawValue = reshape(estResult.doaParamEst, [], 1);
if numel(rawValue) ~= 2 || any(~isfinite(rawValue))
  return;
end
doaParam = rawValue;
end

function doaDriftDeg = localCalcDoaDriftDeg(doaA, doaB)
doaDriftDeg = NaN;
if isempty(doaA) || isempty(doaB)
  return;
end
try
  doaDriftDeg = calcLatlonAngleError(doaA(:), doaB(:));
catch
  doaDriftDeg = NaN;
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
