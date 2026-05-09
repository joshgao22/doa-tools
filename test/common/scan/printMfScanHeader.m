function printMfScanHeader(scanName, scanConfig, runDir)
%PRINTMFSCANHEADER Print compact runtime settings for scan scripts.

arguments
  scanName (1,:) char
  scanConfig (1,1) struct
  runDir = ''
end

fprintf('Running %s ...\n', scanName);
localPrintField(scanConfig, 'numRepeat', 'repeats', '%d');
localPrintField(scanConfig, 'contextSeed', 'context seed', '%d');
localPrintField(scanConfig, 'baseSeed', 'base seed', '%d');
localPrintNumericList(scanConfig, 'seedList', 'seed list');
localPrintNumericList(scanConfig, 'snrDbList', 'SNR list (dB)');
localPrintNumericList(scanConfig, 'frameCountList', 'frame count list');
if localHasField(scanConfig, 'frameIntvlSecList')
  localPrintText('frame interval list (ms)', localFormatNumericList(scanConfig.frameIntvlSecList * 1e3, 'frameIntvlSecList'));
end
localPrintField(scanConfig, 'snrDb', 'snr (dB)', '%.2f');
localPrintField(scanConfig, 'primaryPlotSnrDb', 'primary plot SNR (dB)', '%.2f');
if localHasField(scanConfig, 'primaryFrameIntvlSec')
  localPrintText('primary frame interval (ms)', sprintf('%.6g', scanConfig.primaryFrameIntvlSec * 1e3));
end
localPrintField(scanConfig, 'numTask', 'task count', '%d');
localPrintLogicalField(scanConfig, 'saveSnapshot', 'save snapshot');
localPrintLogicalField(scanConfig, 'checkpointEnable', 'checkpoint');
localPrintLogicalField(scanConfig, 'notifyTelegramEnable', 'telegram notify');
localPrintTextField(scanConfig, 'runKey', 'run key');

if isstring(runDir)
  runDir = char(runDir);
end
if ~isempty(runDir)
  localPrintText('run dir', runDir);
end
end

function localPrintField(dataStruct, fieldName, labelText, formatText)
% Print one scalar field when present.
if ~localHasField(dataStruct, fieldName)
  return;
end
value = dataStruct.(fieldName);
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  localPrintText(labelText, char(string(value)));
  return;
end
fprintf(['  %-32s : ' formatText '\n'], labelText, value);
end

function localPrintLogicalField(dataStruct, fieldName, labelText)
% Print one logical field when present.
if ~localHasField(dataStruct, fieldName)
  return;
end
fprintf('  %-32s : %d\n', labelText, logical(dataStruct.(fieldName)));
end

function localPrintTextField(dataStruct, fieldName, labelText)
% Print one text field when present.
if ~localHasField(dataStruct, fieldName)
  return;
end
localPrintText(labelText, char(string(dataStruct.(fieldName))));
end

function localPrintNumericList(dataStruct, fieldName, labelText)
% Print one numeric vector field when present.
if ~localHasField(dataStruct, fieldName)
  return;
end
value = dataStruct.(fieldName);
if ~(isnumeric(value) || islogical(value))
  localPrintText(labelText, char(string(value)));
  return;
end
localPrintText(labelText, localFormatNumericList(value, fieldName));
end

function localPrintText(labelText, textValue)
% Print one aligned text line.
fprintf('  %-32s : %s\n', labelText, char(string(textValue)));
end

function tf = localHasField(dataStruct, fieldName)
% Return true for a nonempty struct field.
tf = isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName));
end

function textValue = localFormatNumericList(value, fieldName)
% Format numeric vectors for compact logs.
value = reshape(double(value), 1, []);
if isempty(value)
  textValue = '';
  return;
end

if strcmp(fieldName, 'seedList') && localIsConsecutiveSeedList(value)
  textValue = sprintf('%.0f:%.0f (%d seeds)', value(1), value(end), numel(value));
  return;
end

maxFullCount = 16;
maxEdgeCount = 6;
if numel(value) <= maxFullCount
  textValue = localJoinNumericValues(value);
  return;
end

headText = localJoinNumericValues(value(1:maxEdgeCount));
tailText = localJoinNumericValues(value((end - maxEdgeCount + 1):end));
textValue = sprintf('%s, ..., %s (%d values)', headText, tailText, numel(value));
end

function tf = localIsConsecutiveSeedList(value)
% Return true when seeds are an integer unit-stride range.
tf = numel(value) > 1 ...
  && all(isfinite(value)) ...
  && all(abs(value - round(value)) < sqrt(eps)) ...
  && all(diff(round(value)) == 1);
end

function textValue = localJoinNumericValues(value)
% Join numeric values using the same compact scalar format as before.
textValue = strjoin(compose('%.6g', value), ', ');
end
