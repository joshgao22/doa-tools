function printMfReplayHeader(replayName, replayConfig, runDir)
%PRINTMFREPLAYHEADER Print compact runtime settings for replay scripts.

arguments
  replayName (1,:) char
  replayConfig (1,1) struct
  runDir = ''
end

fprintf('Running %s ...\n', replayName);
localPrintField(replayConfig, 'numRepeat', 'repeats', '%d');
localPrintField(replayConfig, 'numSearchRepeat', 'search repeats', '%d');
localPrintField(replayConfig, 'snrDb', 'snr (dB)', '%.2f');
localPrintField(replayConfig, 'baseSeed', 'base seed', '%d');
if isfield(replayConfig, 'contextBaseSeed')
  localPrintField(replayConfig, 'contextBaseSeed', 'context base seed', '%d');
end
if isfield(replayConfig, 'numRepeat') || isfield(replayConfig, 'numSearchRepeat')
  fprintf('  %-32s : %s\n', 'repeat mode', 'parfor-auto');
end
localPrintLogicalField(replayConfig, 'saveSnapshot', 'save snapshot');
localPrintLogicalField(replayConfig, 'checkpointEnable', 'checkpoint');
localPrintLogicalField(replayConfig, 'notifyTelegramEnable', 'telegram notify');
localPrintTextField(replayConfig, 'runKey', 'run key');

if isstring(runDir)
  runDir = char(runDir);
end
if ~isempty(runDir)
  fprintf('  %-32s : %s\n', 'run dir', runDir);
end
end

function localPrintField(dataStruct, fieldName, labelText, formatText)
% Print one scalar field when present.
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
value = dataStruct.(fieldName);
if isempty(value)
  return;
end
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  fprintf('  %-32s : %s\n', labelText, char(string(value)));
  return;
end
fprintf(['  %-32s : ' formatText '\n'], labelText, value);
end

function localPrintLogicalField(dataStruct, fieldName, labelText)
% Print one logical field when present.
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
value = dataStruct.(fieldName);
if isempty(value)
  return;
end
fprintf('  %-32s : %d\n', labelText, logical(value));
end

function localPrintTextField(dataStruct, fieldName, labelText)
% Print one text field when present.
if ~isstruct(dataStruct) || ~isfield(dataStruct, fieldName)
  return;
end
value = dataStruct.(fieldName);
if isempty(value)
  return;
end
fprintf('  %-32s : %s\n', labelText, char(string(value)));
end
