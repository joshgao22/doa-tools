function row = buildMfRuntimeRow(task, methodName, stageName, stageGroup, wallTimeMs)
%BUILDMFRUNTIMEROW Build one timing row.

row = emptyMfRuntimeRow();
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.methodName = string(methodName);
row.stageName = string(stageName);
row.stageGroup = string(stageGroup);
row.wallTimeMs = wallTimeMs;
row.wallTimeSec = wallTimeMs / 1000;
end
