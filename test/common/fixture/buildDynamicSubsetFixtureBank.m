function subsetFixtureCell = buildDynamicSubsetFixtureBank(sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, ...
  masterOffsetIdx, deterministicOffsetIdx, numRandomTrial, scheduleSeed, ...
  gridSize, searchRange, E, wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault, parallelOpt)
%BUILDDYNAMICSUBSETFIXTUREBANK Build one deterministic-plus-random subset bank.
%
% The first subset is always the caller-provided deterministic schedule.
% The remaining subsets are random schedules with the same frame count and a
% fixed zero-offset reference frame. This builder stays serial by default
% because per-fixture construction is usually too light to amortize parfor
% startup and scheduling overhead.

if nargin < 15 || isempty(parallelOpt)
  parallelOpt = struct();
end
if nargin < 7 || isempty(scheduleSeed)
  scheduleSeed = 253;
end
if nargin < 6 || isempty(numRandomTrial)
  numRandomTrial = 0;
end

subsetSize = numel(deterministicOffsetIdx);
if subsetSize < 2
  error('buildDynamicSubsetFixtureBank:InvalidSubsetSize', ...
    'The subset schedule must contain at least two frames.');
end
if ~any(deterministicOffsetIdx == 0)
  error('buildDynamicSubsetFixtureBank:MissingReferenceFrame', ...
    'The deterministic subset schedule must contain one zero-offset frame.');
end

masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
numPickOther = subsetSize - 1;
if numPickOther > numel(masterNoZero)
  error('buildDynamicSubsetFixtureBank:MasterWindowTooShort', ...
    'The master window does not contain enough non-reference frames.');
end

scheduleCell = cell(numRandomTrial + 1, 1);
labelCell = strings(numRandomTrial + 1, 1);
scheduleCell{1} = reshape(deterministicOffsetIdx, 1, []);
labelCell(1) = "deterministic";

rngPrev = rng();
cleanupObj = onCleanup(@() rng(rngPrev));
rng(scheduleSeed);
for iTrial = 1:numRandomTrial
  pickUse = sort([0, masterNoZero(randperm(numel(masterNoZero), numPickOther))]);
  scheduleCell{iTrial + 1} = reshape(pickUse, 1, []);
  labelCell(iTrial + 1) = "random" + string(iTrial);
end
clear cleanupObj;
rng(rngPrev);

subsetFixtureCell = cell(numel(scheduleCell), 1);
useParfor = localShouldUseParfor(parallelOpt, numel(scheduleCell), 'minFixtureBankForParfor');
if useParfor
  parfor iSubset = 1:numel(scheduleCell)
    fixtureUse = buildDynamicFrameSubsetFixture( ...
      sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, ...
      masterOffsetIdx, scheduleCell{iSubset}, gridSize, searchRange, E, ...
      wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault);
    fixtureUse.subsetLabel = labelCell(iSubset);
    fixtureUse.subsetOffsetIdx = reshape(scheduleCell{iSubset}, 1, []);
    subsetFixtureCell{iSubset} = fixtureUse;
  end
else
  for iSubset = 1:numel(scheduleCell)
    fixtureUse = buildDynamicFrameSubsetFixture( ...
      sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, ...
      masterOffsetIdx, scheduleCell{iSubset}, gridSize, searchRange, E, ...
      wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault);
    fixtureUse.subsetLabel = labelCell(iSubset);
    fixtureUse.subsetOffsetIdx = reshape(scheduleCell{iSubset}, 1, []);
    subsetFixtureCell{iSubset} = fixtureUse;
  end
end
end

function useParfor = localShouldUseParfor(parallelOpt, numCase, minFieldName)
%LOCALSHOULDUSEPARFOR Decide whether the current loop should use parfor.

useParfor = logical(localGetFieldOrDefault(parallelOpt, 'enableFixtureBankParfor', ...
  localGetFieldOrDefault(parallelOpt, 'enableParfor', false)));
if ~useParfor
  return;
end
minCase = localGetFieldOrDefault(parallelOpt, minFieldName, 12);
useParfor = useParfor && (numCase >= minCase) && localCanUseParfor();
end

function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether parfor can be used safely here.

tf = false;
if ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  tf = ~isempty(ver('parallel'));
catch
  tf = false;
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
