function scanData = finalizeMfScanResult(scanData, runDir)
%FINALIZEMFSCANRESULT Clean temporary scan artifacts after successful runs.
% Scan scripts save only lightweight scanData. Short scans may pass an empty
% runDir when they do not write checkpoint or intermediate artifacts.

arguments
  scanData (1,1) struct
  runDir = ''
end

if isstring(runDir)
  runDir = char(runDir);
end
if isempty(runDir)
  scanData.tmpDirUsed = false;
  scanData.tmpDirCleaned = true;
  return;
end

scanData.tmpDirUsed = isfolder(runDir);
if isfolder(runDir)
  fprintf('[%s] Clean scan temporary artifacts: %s\n', ...
    char(datetime('now', 'Format', 'HH:mm:ss')), runDir);
  cleanupRunArtifacts(runDir, struct('requiredPrefix', '', 'verbose', false));
end
scanData.tmpDirCleaned = true;
end
