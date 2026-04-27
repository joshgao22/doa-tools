function replayData = finalizeMfReplayResult(replayData, runDir)
%FINALIZEMFREPLAYRESULT Clean temporary replay artifacts after successful runs.
% Replay scripts save only lightweight replayData. Temporary folders are
% removed after successful completion. Fixed single-case replays may pass an
% empty runDir when they do not write intermediate artifacts.

arguments
  replayData (1,1) struct
  runDir = ''
end

if isstring(runDir)
  runDir = char(runDir);
end
if isempty(runDir)
  replayData.tmpDirUsed = false;
  replayData.tmpDirCleaned = true;
  return;
end

replayData.tmpDirUsed = isfolder(runDir);
if isfolder(runDir)
  fprintf('[%s] Clean replay temporary artifacts: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), runDir);
  cleanupRunArtifacts(runDir, struct('requiredPrefix', '', 'verbose', false));
end
replayData.tmpDirCleaned = true;
end
