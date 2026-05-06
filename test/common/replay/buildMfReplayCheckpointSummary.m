function checkpointSummary = buildMfReplayCheckpointSummary(runState)
%BUILDMFREPLAYCHECKPOINTSUMMARY Store lightweight replay checkpoint status.
% Accepts either one runState or a struct array for multi-stage replay scripts.

arguments
  runState struct = struct()
end

checkpointSummary = struct('runName', strings(0, 1), 'runKey', strings(0, 1), ...
  'runDir', strings(0, 1), 'numTask', zeros(0, 1), 'numDone', zeros(0, 1), ...
  'isComplete', false(0, 1), 'cleanedOnSuccess', false, 'cleanupReport', struct([]));
if isempty(runState)
  return;
end

numState = numel(runState);
checkpointSummary.runName = strings(numState, 1);
checkpointSummary.runKey = strings(numState, 1);
checkpointSummary.runDir = strings(numState, 1);
checkpointSummary.numTask = zeros(numState, 1);
checkpointSummary.numDone = zeros(numState, 1);
checkpointSummary.isComplete = false(numState, 1);
for iState = 1:numState
  checkpointSummary.runName(iState) = string(localGetFieldOrDefault(runState(iState), 'runName', ""));
  checkpointSummary.runKey(iState) = string(localGetFieldOrDefault(runState(iState), 'runKey', ""));
  checkpointSummary.runDir(iState) = string(localGetFieldOrDefault(runState(iState), 'runDir', ""));
  checkpointSummary.numTask(iState) = localGetFieldOrDefault(runState(iState), 'numTask', 0);
  checkpointSummary.numDone(iState) = localGetFieldOrDefault(runState(iState), 'numDone', 0);
  checkpointSummary.isComplete(iState) = logical(localGetFieldOrDefault(runState(iState), 'isComplete', false));
end
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
