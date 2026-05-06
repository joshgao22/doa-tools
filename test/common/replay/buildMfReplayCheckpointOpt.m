function checkpointOpt = buildMfReplayCheckpointOpt(replayName, config, opt)
%BUILDMFREPLAYCHECKPOINTOPT Build common replay checkpoint options.
% The caller owns the replay-specific run-key signature and checkpoint meta.
% This helper only resolves the common output root, resume flag, run directory,
% and disabled-checkpoint shape used by replay scripts.

arguments
  replayName
  config (1,1) struct
  opt (1,1) struct = struct()
end

checkpointEnable = localGetLogicalField(config, 'checkpointEnable', false);
checkpointOpt = struct('enable', checkpointEnable);
if ~checkpointEnable
  checkpointOpt.runDir = "";
  return;
end

runName = string(localGetFieldOrDefault(opt, 'runName', replayName));
checkpointOpt.resume = localGetLogicalField(config, 'checkpointResume', true);
checkpointOpt.runName = runName;
checkpointOpt.runKey = string(localGetFieldOrDefault(opt, 'runKey', localBuildDefaultRunKey(config)));
checkpointOpt.outputRoot = char(localGetFieldOrDefault(opt, 'outputRoot', fullfile(localGetRepoRoot(), 'tmp')));
checkpointOpt.runDir = fullfile(checkpointOpt.outputRoot, char(checkpointOpt.runName), char(checkpointOpt.runKey));
checkpointOpt.meta = localGetFieldOrDefault(opt, 'meta', localBuildDefaultMeta(config));
if isfield(opt, 'useParfor') && ~isempty(opt.useParfor)
  checkpointOpt.useParfor = logical(opt.useParfor);
end
end

function runKey = localBuildDefaultRunKey(config)
% Build a conservative fallback run key when the caller does not pass one.
seedList = reshape(double(localGetFieldOrDefault(config, 'seedList', [])), 1, []);
if isempty(seedList)
  baseSeed = localGetFieldOrDefault(config, 'baseSeed', 0);
  numRepeat = localGetFieldOrDefault(config, 'numRepeat', 1);
  seedList = baseSeed + (0:(numRepeat - 1));
end
snrDb = localGetFieldOrDefault(config, 'snrDb', NaN);
runKey = string(sprintf('seed%dto%d_rep%d_snr%.2f', seedList(1), seedList(end), numel(seedList), snrDb));
runKey = replace(runKey, '.', 'p');
runKey = replace(runKey, '-', 'm');
end

function meta = localBuildDefaultMeta(config)
% Keep only stable fields that identify a replay batch.
seedList = localGetFieldOrDefault(config, 'seedList', []);
meta = struct( ...
  'snrDb', localGetFieldOrDefault(config, 'snrDb', NaN), ...
  'seedList', reshape(seedList, 1, []));
end

function repoRoot = localGetRepoRoot()
% Resolve the repository root from this common helper location.
helperDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(helperDir)));
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

function value = localGetLogicalField(dataStruct, fieldName, defaultValue)
% Return a logical struct field value or a fallback.
value = logical(defaultValue);
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = logical(dataStruct.(fieldName));
end
end
