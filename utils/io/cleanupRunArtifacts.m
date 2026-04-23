function cleanupReport = cleanupRunArtifacts(targetPathList, opt)
%CLEANUPRUNARTIFACTS Safely delete run artifacts under repo tmp/cache roots.
%
%   cleanupRunArtifacts(targetPathList)
%   cleanupRunArtifacts(targetPathList, opt)
%
% targetPathList can be a char/string path or a cell array of paths.
%
% Optional fields in opt:
%   allowedRoots    - allowed deletion roots, default {<repoRoot>/tmp,
%                     <repoRoot>/test/data/cache}
%   requiredPrefix  - basename must start with this prefix or one of these
%   verbose         - print concise delete summary (default true)
%   dryRun          - only validate, do not delete (default false)
%
% The helper is intentionally narrow: it only deletes files/directories that
% are strictly inside allowed roots, and it refuses to delete the allowed
% roots themselves.

if nargin < 2 || isempty(opt)
  opt = struct();
end

opt = localResolveOpt(opt);
targetPathList = localNormalizePathList(targetPathList);
cleanupReport = repmat(localEmptyReportRow(), numel(targetPathList), 1);

for iPath = 1:numel(targetPathList)
  rawPath = targetPathList{iPath};
  cleanupReport(iPath).inputPath = rawPath;

  if ~(isfile(rawPath) || isfolder(rawPath))
    cleanupReport(iPath).status = 'missing';
    cleanupReport(iPath).message = 'Target does not exist.';
    continue;
  end

  canonicalPath = localCanonicalPath(rawPath);
  cleanupReport(iPath).canonicalPath = canonicalPath;

  [isAllowed, allowedRoot] = localResolveAllowedRoot(canonicalPath, opt.allowedRoots);
  if ~isAllowed
    error('cleanupRunArtifacts:OutsideAllowedRoots', ...
      'Refuse to delete outside allowed roots: %s', canonicalPath);
  end
  cleanupReport(iPath).allowedRoot = allowedRoot;

  if strcmp(canonicalPath, allowedRoot)
    error('cleanupRunArtifacts:RefuseDeleteRoot', ...
      'Refuse to delete allowed root itself: %s', canonicalPath);
  end

  [~, baseName, ext] = fileparts(canonicalPath);
  itemName = [baseName ext];
  cleanupReport(iPath).itemName = itemName;

  if ~isempty(opt.requiredPrefix)
    hasPrefix = any(cellfun(@(prefix) startsWith(itemName, prefix), opt.requiredPrefix));
    if ~hasPrefix
      error('cleanupRunArtifacts:PrefixMismatch', ...
        'Refuse to delete item without required prefix: %s', canonicalPath);
    end
  end

  if isfolder(canonicalPath)
    cleanupReport(iPath).itemType = 'dir';
  else
    cleanupReport(iPath).itemType = 'file';
  end

  if opt.dryRun
    cleanupReport(iPath).status = 'dryRun';
    cleanupReport(iPath).message = 'Validation passed. No deletion performed.';
    continue;
  end

  if strcmp(cleanupReport(iPath).itemType, 'dir')
    rmdir(canonicalPath, 's');
  else
    delete(canonicalPath);
  end

  cleanupReport(iPath).status = 'deleted';
  cleanupReport(iPath).message = 'Deleted successfully.';
end

if opt.verbose
  deletedCount = nnz(strcmp({cleanupReport.status}, 'deleted'));
  dryRunCount = nnz(strcmp({cleanupReport.status}, 'dryRun'));
  missingCount = nnz(strcmp({cleanupReport.status}, 'missing'));
  fprintf('cleanupRunArtifacts: deleted=%d, dryRun=%d, missing=%d\n', ...
    deletedCount, dryRunCount, missingCount);
end
end


function opt = localResolveOpt(opt)
allowedField = {'allowedRoots', 'requiredPrefix', 'verbose', 'dryRun'};
optField = fieldnames(opt);
extraField = setdiff(optField, allowedField);
if ~isempty(extraField)
  error('cleanupRunArtifacts:UnknownOption', ...
    'Unknown cleanupRunArtifacts option field: %s', strjoin(extraField, ', '));
end

repoRoot = localGetProjectRoot();
defaultAllowedRoots = { ...
  fullfile(repoRoot, 'tmp'), ...
  fullfile(repoRoot, 'test', 'data', 'cache')};

if ~isfield(opt, 'allowedRoots') || isempty(opt.allowedRoots)
  opt.allowedRoots = defaultAllowedRoots;
else
  opt.allowedRoots = localNormalizePathList(opt.allowedRoots);
end
if ~isfield(opt, 'requiredPrefix') || isempty(opt.requiredPrefix)
  opt.requiredPrefix = {};
else
  opt.requiredPrefix = cellstr(string(opt.requiredPrefix(:)'));
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose)
  opt.verbose = true;
end
if ~isfield(opt, 'dryRun') || isempty(opt.dryRun)
  opt.dryRun = false;
end

opt.verbose = logical(opt.verbose);
opt.dryRun = logical(opt.dryRun);
opt.allowedRoots = cellfun(@localCanonicalPath, opt.allowedRoots, 'UniformOutput', false);
end


function targetPathList = localNormalizePathList(targetPathList)
if ischar(targetPathList) || isstring(targetPathList)
  targetPathList = cellstr(string(targetPathList(:)'));
elseif iscell(targetPathList)
  targetPathList = cellfun(@char, cellstr(string(targetPathList(:)')), 'UniformOutput', false);
else
  error('cleanupRunArtifacts:InvalidPathList', ...
    'targetPathList must be a path or a cell array of paths.');
end
end


function [isAllowed, allowedRoot] = localResolveAllowedRoot(canonicalPath, allowedRoots)
isAllowed = false;
allowedRoot = '';
for iRoot = 1:numel(allowedRoots)
  rootPath = allowedRoots{iRoot};
  if strcmp(canonicalPath, rootPath) || startsWith(canonicalPath, [rootPath filesep])
    isAllowed = true;
    allowedRoot = rootPath;
    return;
  end
end
end


function pathStr = localCanonicalPath(pathStr)
if isstring(pathStr)
  pathStr = char(pathStr);
end

pathStr = char(pathStr);
if isempty(pathStr)
  error('cleanupRunArtifacts:EmptyPath', 'Path cannot be empty.');
end

try
  pathStr = char(java.io.File(pathStr).getCanonicalPath());
catch
  pathStr = char(java.io.File(pathStr).getAbsolutePath());
end

pathStr = strrep(pathStr, '/', filesep);
while numel(pathStr) > 1 && any(pathStr(end) == ['/' '\'])
  pathStr(end) = [];
end
end


function repoRoot = localGetProjectRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end


function row = localEmptyReportRow()
row = struct( ...
  'inputPath', '', ...
  'canonicalPath', '', ...
  'allowedRoot', '', ...
  'itemName', '', ...
  'itemType', '', ...
  'status', '', ...
  'message', '');
end
