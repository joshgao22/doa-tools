function checkpointSummaryTable = buildMfScanCheckpointSummaryTable(runStateList, checkpointDir)
%BUILDMFSCANCHECKPOINTSUMMARYTABLE Build compact scan checkpoint summary rows.
% Accepts a scalar or struct-array checkpoint run state returned by
% runPerfTaskGridWithCheckpoint. The helper stores only lightweight fields so
% scanData can be reloaded and summarized without task payloads.

arguments
  runStateList struct = struct([])
  checkpointDir = ""
end

if isstring(checkpointDir)
  checkpointDir = char(checkpointDir);
end
if isempty(runStateList)
  checkpointSummaryTable = table(string(checkpointDir), 0, 0, false, "", ...
    'VariableNames', {'checkpointParentDir', 'numTask', 'numDone', 'isComplete', 'runDir'});
  checkpointSummaryTable(1, :) = [];
  return;
end

runStateList = runStateList(:);
numRun = numel(runStateList);
checkpointParentDir = strings(numRun, 1);
numTask = zeros(numRun, 1);
numDone = zeros(numRun, 1);
isComplete = false(numRun, 1);
runDir = strings(numRun, 1);
for iRun = 1:numRun
  state = runStateList(iRun);
  checkpointParentDir(iRun) = string(checkpointDir);
  numTask(iRun) = localGetFieldOrDefault(state, 'numTask', 0);
  numDone(iRun) = localGetFieldOrDefault(state, 'numDone', 0);
  isComplete(iRun) = logical(localGetFieldOrDefault(state, 'isComplete', false));
  runDir(iRun) = string(localGetFieldOrDefault(state, 'runDir', ""));
end
checkpointSummaryTable = table(checkpointParentDir, numTask, numDone, isComplete, runDir);
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
% Return a nonempty struct field value or a fallback.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
