function [data, meta, loadedVar, inventory] = loadExpSnapshot(snapshotFile, targetWs, opt)
%LOADEXPSNAPSHOT Load one experiment snapshot and optionally restore variables.
%
%   loadExpSnapshot(snapshotFile)
%   loadExpSnapshot(snapshotFile, targetWs)
%   loadExpSnapshot(snapshotFile, targetWs, opt)
%   [data, meta, loadedVar, inventory] = loadExpSnapshot(...)
%
% targetWs:
%   'caller' - restore to caller workspace (default)
%   'base'   - restore to base workspace
%   'none'   - only load snapshot content, do not assign variables
%
% Optional fields in opt:
%   varNames - restore only these saved variables
%   strict   - error when requested varNames are missing (default true)
%   metaOnly - only load meta/inventory, skip variable restore (default false)

if nargin < 2 || isempty(targetWs)
  targetWs = 'caller';
end
if nargin < 3 || isempty(opt)
  opt = struct();
end

if ~any(strcmp(targetWs, {'caller', 'base', 'none'}))
  error('loadExpSnapshot:InvalidTargetWs', ...
    'targetWs must be ''caller'', ''base'', or ''none''.');
end

opt = localResolveOpt(opt);
snapshot = load(snapshotFile);

if ~isfield(snapshot, 'data')
  error('loadExpSnapshot:MissingData', ...
    'The snapshot file does not contain the field "data".');
end

data = snapshot.data;
if ~isstruct(data)
  error('loadExpSnapshot:InvalidData', ...
    'The field "data" must be a struct.');
end

if isfield(snapshot, 'meta')
  meta = snapshot.meta;
else
  meta = struct();
end

if isfield(snapshot, 'inventory')
  inventory = snapshot.inventory;
else
  inventory = localBuildLegacyInventory(data);
end

loadedVar = {};
if opt.metaOnly
  targetWs = 'none';
else
  loadedVar = localResolveLoadedVar(data, opt);
end

if ~strcmp(targetWs, 'none')
  for iVar = 1:numel(loadedVar)
    varName = loadedVar{iVar};
    assignin(targetWs, varName, data.(varName));
  end
end

fprintf('Loaded snapshot from %s\n', snapshotFile);
if opt.metaOnly
  fprintf('Loaded metadata only.\n');
elseif strcmp(targetWs, 'none')
  fprintf('Loaded %d saved variables without workspace restore.\n', numel(loadedVar));
else
  fprintf('Restored %d variables to workspace: %s\n', numel(loadedVar), targetWs);
end

if isfield(meta, 'selectionMode')
  fprintf('Snapshot selection mode: %s\n', meta.selectionMode);
end
end


function opt = localResolveOpt(opt)
allowedField = {'varNames', 'strict', 'metaOnly'};
optField = fieldnames(opt);
extraField = setdiff(optField, allowedField);
if ~isempty(extraField)
  error('loadExpSnapshot:UnknownOption', ...
    'Unknown loadExpSnapshot option field: %s', strjoin(extraField, ', '));
end

if ~isfield(opt, 'varNames') || isempty(opt.varNames)
  opt.varNames = {};
else
  opt.varNames = cellstr(string(opt.varNames(:)'));
end
if ~isfield(opt, 'strict') || isempty(opt.strict)
  opt.strict = true;
end
if ~isfield(opt, 'metaOnly') || isempty(opt.metaOnly)
  opt.metaOnly = false;
end
opt.strict = logical(opt.strict);
opt.metaOnly = logical(opt.metaOnly);
end


function loadedVar = localResolveLoadedVar(data, opt)
savedVar = fieldnames(data);
if isempty(opt.varNames)
  loadedVar = savedVar;
  return;
end

missingVar = setdiff(opt.varNames, savedVar);
if opt.strict && ~isempty(missingVar)
  error('loadExpSnapshot:MissingRequestedVar', ...
    'Requested saved variable not found: %s', strjoin(missingVar, ', '));
end

loadedVar = intersect(savedVar, opt.varNames, 'stable');
end


function inventory = localBuildLegacyInventory(data)
varNameList = fieldnames(data);
inventory = repmat(localEmptyInventoryRow(), numel(varNameList), 1);
for iVar = 1:numel(varNameList)
  varName = varNameList{iVar};
  varValue = data.(varName);
  varInfo = whos('varValue');
  inventory(iVar).name = varName;
  inventory(iVar).class = class(varValue);
  inventory(iVar).size = size(varValue);
  inventory(iVar).bytes = varInfo.bytes;
  inventory(iVar).isSaved = true;
  inventory(iVar).skipReason = '';
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
