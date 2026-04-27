function resultFile = saveExpSnapshot(prefix, opt)
%SAVEEXPSNAPSHOT Save caller workspace snapshot for experiment recovery.
%
%   resultFile = saveExpSnapshot()
%   resultFile = saveExpSnapshot(prefix)
%   resultFile = saveExpSnapshot(prefix, opt)
%
% Default behavior matches the old implementation: save all variables from
% the caller workspace so loadExpSnapshot can restore them later.
%
% Optional fields in opt:
%   includeVars  - caller variable names to save explicitly
%   excludeVars  - caller variable names to skip from full workspace save
%   outputDir    - snapshot directory, default <repoRoot>/test/data/cache
%                  or a task-type subdirectory such as replay/scan/perf when
%                  the caller location identifies one
%   maxVarBytes  - skip auto-collected variables larger than this threshold
%                  includeVars has higher priority and bypasses this filter
%   extraMeta    - extra metadata struct to store in meta.extraMeta
%   verbose      - print one concise save summary (default true)

if nargin < 1
  prefix = '';
end
if nargin < 2 || isempty(opt)
  opt = struct();
end

callerWhos = evalin('caller', 'whos');
callerEntry = evalin('caller', 'mfilename(''fullpath'')');
callerEntryName = evalin('caller', 'mfilename');
callerInfo = localBuildCallerInfo(callerEntry, callerEntryName);
prefix = localResolvePrefix(prefix, callerInfo.entryName);
opt = localResolveOpt(opt, callerInfo, prefix);

if ~isfolder(opt.outputDir)
  mkdir(opt.outputDir);
end

timestamp = datestr(now, 'yyyymmdd-HHMMSS');
resultFile = fullfile(opt.outputDir, sprintf('%s_%s.mat', prefix, timestamp));

[inventory, saveMask] = localPlanSnapshotSelection(callerWhos, opt);
data = struct();
saveIdx = find(saveMask);
for iSave = 1:numel(saveIdx)
  iVar = saveIdx(iSave);
  varName = callerWhos(iVar).name;
  try
    data.(varName) = evalin('caller', varName);
  catch readErr
    inventory(iVar).isSaved = false;
    inventory(iVar).skipReason = sprintf('evalFailed:%s', readErr.identifier);
    if isfield(data, varName)
      data = rmfield(data, varName);
    end
  end
end

meta = localBuildMeta(callerInfo, opt, data, inventory, resultFile);

saveFlag = localChooseSaveFlag(inventory);
try
  localSaveSnapshot(resultFile, data, meta, inventory, saveFlag);
catch firstErr
  if strcmp(saveFlag, '-v7.3')
    rethrow(firstErr);
  end
  if isfile(resultFile)
    delete(resultFile);
  end
  meta.matVersionUsed = '-v7.3';
  localSaveSnapshot(resultFile, data, meta, inventory, '-v7.3');
end

if opt.verbose
  savedMask = [inventory.isSaved];
  fprintf('Saved experiment snapshot to %s (%d/%d vars saved).\n', ...
    resultFile, nnz(savedMask), numel(savedMask));
end
end


function opt = localResolveOpt(opt, callerInfo, prefix)
allowedField = {'includeVars', 'excludeVars', 'outputDir', 'maxVarBytes', 'extraMeta', 'verbose'};
optField = fieldnames(opt);
extraField = setdiff(optField, allowedField);
if ~isempty(extraField)
  error('saveExpSnapshot:UnknownOption', ...
    'Unknown saveExpSnapshot option field: %s', strjoin(extraField, ', '));
end

if ~isfield(opt, 'includeVars') || isempty(opt.includeVars)
  opt.includeVars = {};
else
  opt.includeVars = cellstr(string(opt.includeVars(:)'));
end
if ~isfield(opt, 'excludeVars') || isempty(opt.excludeVars)
  opt.excludeVars = {};
else
  opt.excludeVars = cellstr(string(opt.excludeVars(:)'));
end
if ~isempty(opt.includeVars) && ~isempty(opt.excludeVars)
  error('saveExpSnapshot:ConflictingSelection', ...
    'includeVars and excludeVars cannot be used at the same time.');
end
if ~isfield(opt, 'outputDir') || isempty(opt.outputDir)
  opt.outputDir = localResolveDefaultOutputDir(callerInfo, prefix);
end
if ~isfield(opt, 'maxVarBytes') || isempty(opt.maxVarBytes)
  opt.maxVarBytes = [];
end
if ~isfield(opt, 'extraMeta') || isempty(opt.extraMeta)
  opt.extraMeta = struct();
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose)
  opt.verbose = true;
end

opt.outputDir = char(string(opt.outputDir));
opt.verbose = logical(opt.verbose);
if ~isempty(opt.maxVarBytes)
  opt.maxVarBytes = double(opt.maxVarBytes);
  if ~isscalar(opt.maxVarBytes) || ~isfinite(opt.maxVarBytes) || opt.maxVarBytes <= 0
    error('saveExpSnapshot:InvalidMaxVarBytes', ...
      'maxVarBytes must be one positive finite scalar.');
  end
end
if ~isstruct(opt.extraMeta)
  error('saveExpSnapshot:InvalidExtraMeta', ...
    'extraMeta must be a struct.');
end
end


function outputDir = localResolveDefaultOutputDir(callerInfo, prefix)
% Resolve the default snapshot directory from caller path or filename prefix.
cacheRoot = fullfile(localGetProjectRoot(), 'test', 'data', 'cache');
entryPath = lower(strrep(char(string(localGetStructField(callerInfo, 'entry', ''))), '\', '/'));
entryName = lower(char(string(localGetStructField(callerInfo, 'entryName', ''))));
prefixName = lower(char(string(prefix)));

if ~isempty(strfind(entryPath, '/test/dev/replay/')) || strncmp(prefixName, 'replay', numel('replay')) %#ok<STREMP>
  outputDir = fullfile(cacheRoot, 'replay');
elseif ~isempty(strfind(entryPath, '/test/dev/scan/')) || strncmp(prefixName, 'scan', numel('scan')) %#ok<STREMP>
  outputDir = fullfile(cacheRoot, 'scan');
elseif ~isempty(strfind(entryPath, '/test/regression/')) || strncmp(prefixName, 'regression', numel('regression')) %#ok<STREMP>
  outputDir = fullfile(cacheRoot, 'regression');
elseif (~isempty(strfind(entryPath, '/test/dev/')) && ~isempty(strfind(entryName, 'perf'))) ...
    || ~isempty(strfind(prefixName, 'perf')) %#ok<STREMP>
  outputDir = fullfile(cacheRoot, 'perf');
else
  outputDir = cacheRoot;
end
end


function value = localGetStructField(dataStruct, fieldName, defaultValue)
% Return one struct field when it exists; otherwise keep the default value.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  value = dataStruct.(fieldName);
end
end


function callerInfo = localBuildCallerInfo(callerEntry, callerEntryName)
callerInfo = struct();
callerInfo.entry = char(string(callerEntry));
callerInfo.entryName = char(string(callerEntryName));
callerInfo.code = '';

if isempty(callerInfo.entry) || isempty(callerInfo.entryName)
  callerFrame = localFindExternalCallerFrame();
  if ~isempty(callerFrame)
    [stackDir, stackName] = fileparts(callerFrame.file);
    if isempty(callerInfo.entryName)
      callerInfo.entryName = stackName;
    end
    if isempty(callerInfo.entry)
      callerInfo.entry = fullfile(stackDir, stackName);
    end
  end
end

if ~isempty(callerInfo.entry)
  entryFile = callerInfo.entry;
  if ~endsWith(entryFile, '.m')
    entryFile = [entryFile '.m'];
  end
  if isfile(entryFile)
    try
      callerInfo.code = fileread(entryFile);
    catch
      callerInfo.code = '';
    end
  end
end
end


function callerFrame = localFindExternalCallerFrame()
% Return the first stack frame outside this snapshot helper file.
callerFrame = [];
stack = dbstack('-completenames');
if isempty(stack)
  return;
end
thisFile = localNormalizeMFilePath(stack(1).file);
for iFrame = 1:numel(stack)
  stackFile = localNormalizeMFilePath(stack(iFrame).file);
  if ~strcmp(stackFile, thisFile)
    callerFrame = stack(iFrame);
    return;
  end
end
end


function pathText = localNormalizeMFilePath(pathText)
% Normalize one m-file path for stack-frame comparisons.
pathText = char(string(pathText));
if isempty(pathText)
  return;
end
if ~endsWith(lower(pathText), '.m')
  pathText = [pathText '.m'];
end
pathText = lower(strrep(pathText, '\', '/'));
end


function prefix = localResolvePrefix(prefix, entryName)
if nargin < 1 || isempty(prefix)
  prefix = entryName;
end
if isempty(prefix)
  prefix = 'exp_snapshot';
end
prefix = char(string(prefix));
prefix = regexprep(prefix, '[\\/:*?"<>|\s]+', '_');
if isempty(prefix)
  prefix = 'exp_snapshot';
end
end


function [inventory, saveMask] = localPlanSnapshotSelection(callerWhos, opt)
varNameList = {callerWhos.name};
inventory = repmat(localEmptyInventoryRow(), numel(callerWhos), 1);
saveMask = false(numel(callerWhos), 1);

includeMode = ~isempty(opt.includeVars);
if includeMode
  missingVar = setdiff(opt.includeVars, varNameList);
  if ~isempty(missingVar)
    error('saveExpSnapshot:MissingIncludeVar', ...
      'Included caller variable not found: %s', strjoin(missingVar, ', '));
  end
end

for iVar = 1:numel(callerWhos)
  varInfo = callerWhos(iVar);
  varName = varInfo.name;

  inventory(iVar).name = varName;
  inventory(iVar).class = varInfo.class;
  inventory(iVar).size = varInfo.size;
  inventory(iVar).bytes = varInfo.bytes;

  [shouldSave, skipReason] = localSelectVariable(varName, varInfo.bytes, opt, includeMode);
  inventory(iVar).isSaved = shouldSave;
  inventory(iVar).skipReason = skipReason;
  saveMask(iVar) = shouldSave;
end
end


function [shouldSave, skipReason] = localSelectVariable(varName, varBytes, opt, includeMode)
if includeMode
  shouldSave = any(strcmp(varName, opt.includeVars));
  if shouldSave
    skipReason = '';
  else
    skipReason = 'notIncluded';
  end
  return;
end

if any(strcmp(varName, opt.excludeVars))
  shouldSave = false;
  skipReason = 'excluded';
  return;
end

if ~isempty(opt.maxVarBytes) && varBytes > opt.maxVarBytes
  shouldSave = false;
  skipReason = 'tooLarge';
  return;
end

shouldSave = true;
skipReason = '';
end


function meta = localBuildMeta(callerInfo, opt, data, inventory, resultFile)
savedVarNames = fieldnames(data);
skipMask = ~[inventory.isSaved];
meta = struct();
meta.time = datetime('now');
meta.entry = callerInfo.entry;
meta.entryName = callerInfo.entryName;
meta.code = callerInfo.code;
meta.outputDir = fileparts(resultFile);
meta.resultFile = resultFile;
meta.savedVarNames = savedVarNames;
meta.skippedVarNames = {inventory(skipMask).name};
meta.extraMeta = opt.extraMeta;
meta.matVersionUsed = 'default';

if ~isempty(opt.includeVars)
  meta.selectionMode = 'include';
elseif ~isempty(opt.excludeVars)
  meta.selectionMode = 'exclude';
elseif ~isempty(opt.maxVarBytes)
  meta.selectionMode = 'allWithSizeLimit';
else
  meta.selectionMode = 'all';
end

meta.maxVarBytes = opt.maxVarBytes;
end


function saveFlag = localChooseSaveFlag(inventory)
savedMask = [inventory.isSaved];
savedBytes = [inventory(savedMask).bytes];
if isempty(savedBytes)
  saveFlag = 'default';
  return;
end

if any(savedBytes >= 2^31) || sum(double(savedBytes)) >= 2^31
  saveFlag = '-v7.3';
else
  saveFlag = 'default';
end
end


function localSaveSnapshot(resultFile, data, meta, inventory, saveFlag)
if strcmp(saveFlag, 'default')
  save(resultFile, 'data', 'meta', 'inventory');
else
  save(resultFile, 'data', 'meta', 'inventory', saveFlag);
end
end


function row = localEmptyInventoryRow()
row = struct( ...
  'name', '', ...
  'class', '', ...
  'size', [], ...
  'bytes', 0, ...
  'isSaved', false, ...
  'skipReason', '');
end


function projectRoot = localGetProjectRoot()
thisFile = mfilename('fullpath');
projectRoot = fileparts(fileparts(fileparts(thisFile)));
end
