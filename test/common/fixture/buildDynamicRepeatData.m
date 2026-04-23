function repeatData = buildDynamicRepeatData(context, snrDb, taskSeed)
%BUILDDYNAMICREPEATDATA Build one repeat-level snapshot and fixture bundle.
% The simplified dynamic dev/regression entry points share the same repeat
% preparation path so that subset fixtures, periodic fixtures, and truth are
% constructed identically.

arguments
  context (1,1) struct
  snrDb (1,1) double
  taskSeed (1,1) double
end

rng(taskSeed);

pwrNoise = 1 / (10^(snrDb / 10));
pathGainCellMaster = repmat({ones(context.sceneSeqMaster.numSat, context.sceneSeqMaster.numUser)}, ...
  1, context.sceneSeqMaster.numFrame);

snapshotTimer = tic;
[rxSigCellMaster, ~, ~, ~, ~] = genMultiFrameSnapshots( ...
  context.sceneSeqMaster, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  pwrNoise, pathGainCellMaster, context.simOpt);
genSnapshotSec = toc(snapshotTimer);

[periodicFixture, subsetFixtureCell, fixtureTiming] = buildDynamicRepeatFixtures( ...
  context.sceneSeqMaster, context.linkParamCellMaster, rxSigCellMaster, ...
  context.masterOffsetIdx, context.periodicOffsetIdx, ...
  context.subsetOffsetCell, context.subsetLabelList, context.numSubsetRandomTrial, taskSeed, ...
  context.gridSize, context.searchRange, context.E, context.wavelen, context.waveInfo.sampleRate, ...
  context.fdRangeDefault, context.fdRateRangeDefault, context.parallelOpt);

repeatData = struct();
repeatData.taskSeed = taskSeed;
repeatData.snrDb = snrDb;
repeatData.rxSigCellMaster = rxSigCellMaster;
repeatData.periodicFixture = periodicFixture;
repeatData.subsetFixtureCell = subsetFixtureCell;
repeatData.truth = periodicFixture.truth;
repeatData.genSnapshotSec = genSnapshotSec;
repeatData.buildPeriodicFixtureSec = fixtureTiming.buildPeriodicFixtureSec;
repeatData.buildSubsetFixtureSec = fixtureTiming.buildSubsetFixtureSec;
end
