function cleanupReport = cleanupMfScanCheckpointRuns(runStateList, opt)
%CLEANUPMFSCANCHECKPOINTRUNS Clean completed scan checkpoint run directories.
% The caller owns when cleanup is safe. This helper only loops over completed
% run states, calls cleanupPerfTaskGridCheckpoint, and returns a fixed-shape
% lightweight report for scanData.

arguments
  runStateList struct = struct([])
  opt (1,1) struct = struct()
end

warningId = string(localGetFieldOrDefault(opt, 'warningId', 'cleanupMfScanCheckpointRuns:CheckpointCleanupFailed'));
cleanupOpt = localGetFieldOrDefault(opt, 'cleanupOpt', struct('verbose', false, 'logEnable', true, 'removeEmptyParent', true));

cleanupReport = localEmptyCheckpointCleanupReportTable();
for iRun = 1:numel(runStateList)
  runState = runStateList(iRun);
  if ~isstruct(runState) || ~isfield(runState, 'isComplete') || ~runState.isComplete
    continue;
  end
  row = localBuildCheckpointCleanupRow(runState);
  try
    reportTmp = cleanupPerfTaskGridCheckpoint(runState, cleanupOpt);
    row.cleanupStatus = "ok";
    row.cleanupMessage = localSummarizeCleanupReport(reportTmp);
  catch ME
    row.cleanupStatus = "failed";
    row.cleanupMessage = string(ME.message);
    warning(char(warningId), 'Checkpoint cleanup failed for %s: %s', char(row.runDir), ME.message);
  end
  cleanupReport = [cleanupReport; struct2table(row)]; %#ok<AGROW>
end
end

function cleanupReport = localEmptyCheckpointCleanupReportTable()
% Return an empty fixed-field cleanup table.
row = localBuildCheckpointCleanupRow(struct());
cleanupReport = struct2table(row);
cleanupReport(1, :) = [];
end

function row = localBuildCheckpointCleanupRow(runState)
% Build one fixed-field checkpoint cleanup row.
row = struct();
row.runName = string(localGetFieldOrDefault(runState, 'runName', ""));
row.runKey = string(localGetFieldOrDefault(runState, 'runKey', ""));
row.runDir = string(localGetFieldOrDefault(runState, 'runDir', ""));
row.numTask = double(localGetFieldOrDefault(runState, 'numTask', 0));
row.numDone = double(localGetFieldOrDefault(runState, 'numDone', 0));
row.isComplete = logical(localGetFieldOrDefault(runState, 'isComplete', false));
row.cleanupStatus = "not-run";
row.cleanupMessage = "";
end

function message = localSummarizeCleanupReport(reportTmp)
% Return a compact text summary for variable-shape cleanup output.
if isempty(reportTmp)
  message = "empty cleanup report";
  return;
end
if isstruct(reportTmp)
  message = sprintf('cleanup report rows=%d fields=%s', numel(reportTmp), strjoin(fieldnames(reportTmp), ','));
elseif istable(reportTmp)
  message = sprintf('cleanup report rows=%d variables=%s', height(reportTmp), strjoin(reportTmp.Properties.VariableNames, ','));
else
  message = sprintf('cleanup report class=%s', class(reportTmp));
end
message = string(message);
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
