function dataPath = resolveTestDataPath(fileName, subDir, extraDirList)
%RESOLVETESTDATAPATH Resolve one file from common test-data locations.
%
%Syntax:
%  dataPath = resolveTestDataPath(fileName)
%  dataPath = resolveTestDataPath(fileName, subDir)
%  dataPath = resolveTestDataPath(fileName, subDir, extraDirList)
%
%Inputs:
%  fileName      - file name or relative path to resolve
%  subDir        - optional test/data subdirectory, e.g. "tle" or "cache"
%  extraDirList  - optional string/cellstr list of extra directories to try

arguments
  fileName {mustBeTextScalar}
  subDir {mustBeTextScalar} = ""
  extraDirList = {}
end

fileName = char(string(fileName));
subDir = char(string(subDir));
extraDirList = localNormalizeExtraDirList(extraDirList);

[thisDir, ~, ~] = fileparts(mfilename('fullpath'));
testDir = fileparts(thisDir);
dataRoot = fullfile(testDir, 'data');

candidateRoot = {pwd, dataRoot, fullfile(pwd, 'test', 'data'), fullfile('test', 'data'), fullfile('..', 'data')};
candidateRoot = [extraDirList(:).', candidateRoot];

candidateList = {};
if isempty(subDir)
  candidateList{end + 1} = fileName;
  for iRoot = 1:numel(candidateRoot)
    candidateList{end + 1} = fullfile(candidateRoot{iRoot}, fileName);
  end
else
  candidateList{end + 1} = fullfile(subDir, fileName);
  for iRoot = 1:numel(candidateRoot)
    candidateList{end + 1} = fullfile(candidateRoot{iRoot}, subDir, fileName);
  end
end

candidateList = localUniquePathList(candidateList);
for iPath = 1:numel(candidateList)
  if exist(candidateList{iPath}, 'file')
    dataPath = candidateList{iPath};
    return;
  end
end

error('resolveTestDataPath:FileNotFound', ...
  'Unable to locate "%s" under the common test-data paths.', fileName);
end


function extraDirList = localNormalizeExtraDirList(extraDirList)
%LOCALNORMALIZEEXTRADIRLIST Normalize extra directory input to cellstr.

if isempty(extraDirList)
  extraDirList = {};
  return;
end
if ischar(extraDirList) || isstring(extraDirList)
  extraDirList = cellstr(string(extraDirList));
  return;
end
if iscell(extraDirList)
  extraDirList = cellfun(@char, cellstr(string(extraDirList)), 'UniformOutput', false);
  return;
end
error('resolveTestDataPath:InvalidExtraDirList', ...
  'extraDirList must be text or a cell array of text.');
end


function pathList = localUniquePathList(pathList)
%LOCALUNIQUEPATHLIST Keep the first occurrence of each candidate path.

isKeep = true(size(pathList));
for iPath = 1:numel(pathList)
  if ~isKeep(iPath)
    continue;
  end
  matchMask = strcmp(pathList{iPath}, pathList);
  matchIdx = find(matchMask);
  if numel(matchIdx) > 1
    isKeep(matchIdx(2:end)) = false;
  end
end
pathList = pathList(isKeep);
end
