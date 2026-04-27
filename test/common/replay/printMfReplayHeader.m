function printMfReplayHeader(replayName, replayConfig, runDir)
%PRINTMFREPLAYHEADER Print compact runtime settings for replay functions.

arguments
  replayName (1,:) char
  replayConfig (1,1) struct
  runDir = ''
end

fprintf('Running %s ...\n', replayName);
if isfield(replayConfig, 'numRepeat')
  fprintf('  repeats                         : %d\n', localGetFieldOrDefault(replayConfig, 'numRepeat', NaN));
end
fprintf('  snr (dB)                        : %.2f\n', localGetFieldOrDefault(replayConfig, 'snrDb', NaN));
fprintf('  base seed                       : %d\n', localGetFieldOrDefault(replayConfig, 'baseSeed', NaN));
if isfield(replayConfig, 'numRepeat')
  fprintf('  repeat mode                     : %s\n', 'parfor-auto');
end
fprintf('  save snapshot                   : %d\n', logical(localGetFieldOrDefault(replayConfig, 'saveSnapshot', false)));
if isstring(runDir)
  runDir = char(runDir);
end
if ~isempty(runDir)
  fprintf('  run dir                         : %s\n', runDir);
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
