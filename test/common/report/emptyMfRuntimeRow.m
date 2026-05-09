function row = emptyMfRuntimeRow()
%EMPTYMFRUNTIMEROW Return one typed runtime row.

row = struct('snrDb', NaN, 'taskSeed', NaN, 'methodName', "", ...
  'stageName', "", 'stageGroup', "", 'wallTimeMs', NaN, 'wallTimeSec', NaN);
end
