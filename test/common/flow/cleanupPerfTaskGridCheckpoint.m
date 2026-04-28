function cleanupReport = cleanupPerfTaskGridCheckpoint(runState, opt)
%CLEANUPPERFTASKGRIDCHECKPOINT Delete one completed checkpoint run directory.
% The helper also removes the script-level checkpoint parent directory when it
% becomes empty, keeping tmp cleanup policy next to the checkpoint runner.

arguments
  runState (1,1) struct
  opt (1,1) struct = struct()
end

opt = localResolveOpt(opt);
cleanupReport = struct([]);
runDir = localGetFieldOrDefault(runState, 'runDir', "");
if strlength(string(runDir)) == 0 || ~isfolder(runDir)
  return;
end

if opt.logEnable
  fprintf('[%s] Clean checkpoint artifacts: %s\n', ...
    char(datetime('now', 'Format', 'HH:mm:ss')), char(runDir));
end

runPrefix = localResolveRunPrefix(runState, runDir);
cleanupReport = cleanupRunArtifacts(runDir, struct( ...
  'requiredPrefix', runPrefix, ...
  'verbose', opt.verbose));

if opt.removeEmptyParent
  parentReport = localCleanupEmptyRunParent(runState, runDir, opt);
  if ~isempty(parentReport)
    if isempty(cleanupReport)
      cleanupReport = parentReport(:);
    else
      cleanupReport = [cleanupReport(:); parentReport(:)];
    end
  end
end
end


function opt = localResolveOpt(opt)
%LOCALRESOLVEOPT Resolve cleanup defaults for checkpoint artifacts.
allowedField = {'verbose', 'logEnable', 'removeEmptyParent'};
optField = fieldnames(opt);
extraField = setdiff(optField, allowedField);
if ~isempty(extraField)
  error('cleanupPerfTaskGridCheckpoint:UnknownOption', ...
    'Unknown checkpoint cleanup option field: %s', strjoin(extraField, ', '));
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose)
  opt.verbose = false;
end
if ~isfield(opt, 'logEnable') || isempty(opt.logEnable)
  opt.logEnable = true;
end
if ~isfield(opt, 'removeEmptyParent') || isempty(opt.removeEmptyParent)
  opt.removeEmptyParent = true;
end
opt.verbose = logical(opt.verbose);
opt.logEnable = logical(opt.logEnable);
opt.removeEmptyParent = logical(opt.removeEmptyParent);
end


function runPrefix = localResolveRunPrefix(runState, runDir)
%LOCALRESOLVERUNPREFIX Choose the required prefix for deleting the run folder.
runKey = localGetFieldOrDefault(runState, 'runKey', "");
if strlength(string(runKey)) > 0
  runPrefix = char(runKey);
  return;
end
[~, runPrefix] = fileparts(char(runDir));
end


function cleanupReport = localCleanupEmptyRunParent(runState, runDir, opt)
%LOCALCLEANUPEMPTYRUNPARENT Remove an empty script-level checkpoint parent.
cleanupReport = struct([]);
parentDir = fileparts(char(runDir));
if isempty(parentDir) || ~isfolder(parentDir)
  return;
end
entryInfo = dir(parentDir);
entryName = {entryInfo.name};
hasContent = any(~ismember(entryName, {'.', '..'}));
if hasContent
  return;
end
parentPrefix = localResolveParentPrefix(runState, parentDir);
cleanupReport = cleanupRunArtifacts(parentDir, struct( ...
  'requiredPrefix', parentPrefix, ...
  'verbose', opt.verbose));
end


function parentPrefix = localResolveParentPrefix(runState, parentDir)
%LOCALRESOLVEPARENTPREFIX Choose the required prefix for the run parent.
runName = localGetFieldOrDefault(runState, 'runName', "");
if strlength(string(runName)) > 0
  parentPrefix = char(runName);
  return;
end
[~, parentPrefix] = fileparts(char(parentDir));
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read a nonempty struct field or return a default.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
