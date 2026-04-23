function subsetFixtureCell = buildDynamicSubsetFixturesByLabel( ...
  sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, masterOffsetIdx, ...
  explicitOffsetCell, explicitLabelList, numRandomTrial, scheduleSeed, randomBaseOffsetIdx, ...
  gridSize, searchRange, E, wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault, parallelOpt)
%BUILDDYNAMICSUBSETFIXTURESBYLABEL Build curated and random dynamic subset fixtures.
% Explicit subset schedules are built first and keep their caller-provided
% labels. Optional random schedules are appended afterwards using the same
% subset size while always keeping the zero-offset reference frame.

arguments
  sceneSeqMaster (1, 1) struct
  linkParamCellMaster (1, :) cell
  rxSigCellMaster (1, :) cell
  masterOffsetIdx (1, :) double
  explicitOffsetCell (1, :) cell
  explicitLabelList (:, 1) string
  numRandomTrial (1, 1) double {mustBeNonnegative, mustBeInteger}
  scheduleSeed (1, 1) double
  randomBaseOffsetIdx (1, :) double
  gridSize (1, 2) double
  searchRange (2, 2) double
  E (1, 1)
  wavelen (1, 1) double {mustBePositive}
  sampleRate (1, 1) double {mustBePositive}
  fdRangeDefault (1, 2) double
  fdRateRangeDefault (1, 2) double
  parallelOpt (1, 1) struct = struct()
end

if isempty(explicitOffsetCell)
  explicitOffsetCell = {};
end
if isempty(explicitLabelList)
  explicitLabelList = strings(0, 1);
else
  explicitLabelList = reshape(string(explicitLabelList), [], 1);
end
if isempty(numRandomTrial)
  numRandomTrial = 0;
end
if isempty(scheduleSeed)
  scheduleSeed = 253;
end
if isempty(randomBaseOffsetIdx)
  randomBaseOffsetIdx = [0, 1];
end

numExplicit = numel(explicitOffsetCell);
if numExplicit ~= numel(explicitLabelList)
  error('buildDynamicSubsetFixturesByLabel:SizeMismatch', ...
    'explicitOffsetCell and explicitLabelList must have the same length.');
end

if isempty(explicitOffsetCell) && numRandomTrial <= 0
  subsetFixtureCell = cell(0, 1);
  return;
end

randomBaseOffsetIdx = reshape(randomBaseOffsetIdx, 1, []);
if ~any(randomBaseOffsetIdx == 0)
  error('buildDynamicSubsetFixturesByLabel:MissingReferenceFrame', ...
    'randomBaseOffsetIdx must contain one zero-offset frame.');
end

scheduleCell = cell(numExplicit + numRandomTrial, 1);
labelCell = strings(numExplicit + numRandomTrial, 1);
for iSubset = 1:numExplicit
  scheduleCell{iSubset} = reshape(explicitOffsetCell{iSubset}, 1, []);
  labelCell(iSubset) = string(explicitLabelList(iSubset));
end

if numRandomTrial > 0
  randomScheduleCell = localBuildRandomScheduleCell(masterOffsetIdx, randomBaseOffsetIdx, numRandomTrial, scheduleSeed);
  for iTrial = 1:numRandomTrial
    outIdx = numExplicit + iTrial;
    scheduleCell{outIdx} = randomScheduleCell{iTrial};
    labelCell(outIdx) = "random" + string(iTrial);
  end
end

[scheduleCell, labelCell] = localDeduplicateScheduleCell(scheduleCell, labelCell);

subsetFixtureCell = cell(numel(scheduleCell), 1);
useParfor = localShouldUseParfor(parallelOpt, numel(scheduleCell), 'minFixtureBankForParfor');
if useParfor
  parfor iSubset = 1:numel(scheduleCell)
    subsetFixtureCell{iSubset} = localBuildOneFixture( ...
      sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, masterOffsetIdx, ...
      scheduleCell{iSubset}, labelCell(iSubset), gridSize, searchRange, E, ...
      wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault);
  end
else
  for iSubset = 1:numel(scheduleCell)
    subsetFixtureCell{iSubset} = localBuildOneFixture( ...
      sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, masterOffsetIdx, ...
      scheduleCell{iSubset}, labelCell(iSubset), gridSize, searchRange, E, ...
      wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault);
  end
end
end


function fixtureUse = localBuildOneFixture(sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, masterOffsetIdx, ...
  subsetOffsetIdx, subsetLabel, gridSize, searchRange, E, wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault)
%LOCALBUILDONEFIXTURE Build and label one subset fixture.

fixtureUse = buildDynamicFrameSubsetFixture( ...
  sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, ...
  masterOffsetIdx, subsetOffsetIdx, gridSize, searchRange, E, ...
  wavelen, sampleRate, fdRangeDefault, fdRateRangeDefault);
fixtureUse.subsetLabel = string(subsetLabel);
fixtureUse.subsetOffsetIdx = reshape(subsetOffsetIdx, 1, []);
end


function scheduleCell = localBuildRandomScheduleCell(masterOffsetIdx, randomBaseOffsetIdx, numRandomTrial, scheduleSeed)
%LOCALBUILDRANDOMSCHEDULECELL Build random subset schedules with one zero frame.

subsetSize = numel(randomBaseOffsetIdx);
if subsetSize < 2
  error('buildDynamicSubsetFixturesByLabel:InvalidSubsetSize', ...
    'The subset schedule must contain at least two frames.');
end
masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
numPickOther = subsetSize - 1;
if numPickOther > numel(masterNoZero)
  error('buildDynamicSubsetFixturesByLabel:MasterWindowTooShort', ...
    'The master window does not contain enough non-reference frames.');
end

rngPrev = rng();
cleanupObj = onCleanup(@() rng(rngPrev));
rng(scheduleSeed);
scheduleCell = cell(numRandomTrial, 1);
for iTrial = 1:numRandomTrial
  pickUse = sort([0, masterNoZero(randperm(numel(masterNoZero), numPickOther))]);
  scheduleCell{iTrial} = reshape(pickUse, 1, []);
end
clear cleanupObj;
rng(rngPrev);
end


function [scheduleCellOut, labelCellOut] = localDeduplicateScheduleCell(scheduleCellIn, labelCellIn)
%LOCALDEDUPLICATESCHEDULECELL Drop duplicate subset schedules while preserving order.
% Random rescue schedules may occasionally reproduce one curated schedule or
% another already-generated random schedule. Building and solving duplicate
% fixtures only adds runtime and does not change the final tooth-selection
% pool, so keep the first occurrence and discard later duplicates.

numSchedule = numel(scheduleCellIn);
keepMask = true(numSchedule, 1);
seenScheduleCell = cell(0, 1);
for iSchedule = 1:numSchedule
  scheduleUse = sort(reshape(scheduleCellIn{iSchedule}, 1, []));
  if localHasMatchingSchedule(seenScheduleCell, scheduleUse)
    keepMask(iSchedule) = false;
    continue;
  end
  seenScheduleCell{end + 1, 1} = scheduleUse; %#ok<AGROW>
end
scheduleCellOut = scheduleCellIn(keepMask);
labelCellOut = labelCellIn(keepMask);
end


function tf = localHasMatchingSchedule(scheduleCell, scheduleUse)
%LOCALHASMATCHINGSCHEDULE Return true when one identical schedule already exists.

tf = false;
scheduleUse = sort(reshape(scheduleUse, 1, []));
for iSchedule = 1:numel(scheduleCell)
  scheduleKnown = reshape(scheduleCell{iSchedule}, 1, []);
  if numel(scheduleKnown) ~= numel(scheduleUse)
    continue;
  end
  if isequal(sort(scheduleKnown), scheduleUse)
    tf = true;
    return;
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
