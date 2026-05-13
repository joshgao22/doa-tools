function snapshotFile = resolveTestSnapshotFile(preferredSnapshotFile, fallbackPattern)
%RESOLVETESTSNAPSHOTFILE Resolve an experiment snapshot with newest-pattern fallback.
%
%Syntax:
%  snapshotFile = resolveTestSnapshotFile(preferredSnapshotFile, fallbackPattern)
%
%The preferred snapshot is returned when it exists. If it is missing, the
%newest file matching fallbackPattern in the same directory is returned.
%This helper is intended for paper/replay/scan redraw scripts where a
%representative snapshot is documented but a local rerun may have a newer
%timestamp.

if nargin < 2 || isempty(fallbackPattern)
  fallbackPattern = "";
end

preferredSnapshotFile = char(string(preferredSnapshotFile));
fallbackPattern = char(string(fallbackPattern));

if exist(preferredSnapshotFile, 'file') == 2
  snapshotFile = preferredSnapshotFile;
  return;
end

cacheDir = fileparts(preferredSnapshotFile);
if isempty(fallbackPattern)
  error('resolveTestSnapshotFile:MissingSnapshot', ...
    'Snapshot file is missing: %s', preferredSnapshotFile);
end

matchInfo = dir(fullfile(cacheDir, fallbackPattern));
if isempty(matchInfo)
  error('resolveTestSnapshotFile:MissingSnapshot', ...
    ['Snapshot file is missing: %s\n' ...
     'No fallback snapshot matches pattern %s in %s.'], ...
    preferredSnapshotFile, fallbackPattern, cacheDir);
end

[~, newestIdx] = max([matchInfo.datenum]);
snapshotFile = fullfile(matchInfo(newestIdx).folder, matchInfo(newestIdx).name);
warning('resolveTestSnapshotFile:UsingFallbackSnapshot', ...
  'Preferred snapshot is missing. Using newest fallback snapshot: %s', snapshotFile);
end
