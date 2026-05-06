function ok = notifyMfReplayStatus(statusOpt)
%NOTIFYMFREPLAYSTATUS Send a compact best-effort HTML replay notification.
% This helper only formats the common replay status shell. Replay-specific
% metric lines must be passed in by the caller as HTML-ready strings.

arguments
  statusOpt (1,1) struct
end

ok = false;
config = localGetFieldOrDefault(statusOpt, 'config', struct());
if ~localGetLogicalField(config, 'notifyTelegramEnable', false)
  return;
end

try
  replayName = string(localGetFieldOrDefault(statusOpt, 'replayName', 'replay'));
  statusText = upper(string(localGetFieldOrDefault(statusOpt, 'statusText', 'DONE')));
  messageText = localBuildReplayStatusHtml(statusOpt, replayName, statusText);
  titleText = sprintf('MATLAB REPLAY - %s', char(statusText));
  ok = notifyTelegram(titleText, messageText, "HTML");
catch notifyME
  warning('Replay Telegram notification wrapper failed: %s', notifyME.message);
end
end

function messageText = localBuildReplayStatusHtml(statusOpt, replayName, statusText)
% Build the common HTML message shell.
config = localGetFieldOrDefault(statusOpt, 'config', struct());
elapsedSec = localGetFieldOrDefault(statusOpt, 'elapsedSec', NaN);
snapshotFile = string(localGetFieldOrDefault(statusOpt, 'snapshotFile', ''));
checkpointDir = string(localGetFieldOrDefault(statusOpt, 'checkpointDir', ''));
subtitleText = string(localGetFieldOrDefault(statusOpt, 'subtitleText', ...
  'Continuous-phase multi-frame DoA-Doppler replay'));
metricLineList = localAsStringColumn(localGetFieldOrDefault(statusOpt, 'metricLineList', strings(0, 1)));
commentLineList = localAsStringColumn(localGetFieldOrDefault(statusOpt, 'commentLineList', strings(0, 1)));
errorObj = localGetFieldOrDefault(statusOpt, 'errorObj', []);

statusIcon = localStatusIcon(statusText);
lineList = strings(0, 1);
lineList(end + 1, 1) = "<b>Simulation report</b>";
lineList(end + 1, 1) = "<i>" + localHtmlEscape(subtitleText) + "</i>";
lineList(end + 1, 1) = "";
lineList(end + 1, 1) = "<b>Status</b>: " + statusIcon + " " + localHtmlEscape(statusText);
lineList(end + 1, 1) = "<b>Script</b>: <code>" + localHtmlEscape(replayName) + "</code>";
if isfinite(double(elapsedSec))
  lineList(end + 1, 1) = "<b>Elapsed</b>: <code>" + localHtmlEscape(localFormatElapsed(elapsedSec)) + "</code>";
end
lineList = [lineList; localBuildConfigLines(config)]; %#ok<AGROW>

if strlength(snapshotFile) > 0
  lineList(end + 1, 1) = "";
  lineList(end + 1, 1) = "<b>Snapshot</b>";
  lineList(end + 1, 1) = "<code>" + localHtmlEscape(snapshotFile) + "</code>";
end
if strlength(checkpointDir) > 0
  lineList(end + 1, 1) = "";
  lineList(end + 1, 1) = "<b>Checkpoint</b>";
  lineList(end + 1, 1) = "<code>" + localHtmlEscape(checkpointDir) + "</code>";
end
if ~isempty(metricLineList)
  lineList(end + 1, 1) = "";
  lineList(end + 1, 1) = "<b>Key metrics</b>";
  lineList = [lineList; metricLineList]; %#ok<AGROW>
end
if ~isempty(commentLineList)
  lineList(end + 1, 1) = "";
  lineList(end + 1, 1) = "<b>Comment</b>";
  for iLine = 1:numel(commentLineList)
    lineList(end + 1, 1) = localHtmlEscape(commentLineList(iLine)); %#ok<AGROW>
  end
end
if isa(errorObj, 'MException')
  lineList(end + 1, 1) = "";
  lineList(end + 1, 1) = "<b>Error</b>";
  lineList(end + 1, 1) = "<code>" + localHtmlEscape(errorObj.identifier) + "</code>";
  lineList(end + 1, 1) = localHtmlEscape(errorObj.message);
end

messageText = strjoin(lineList, newline);
end

function lineList = localBuildConfigLines(config)
% Build compact common configuration lines.
lineList = strings(0, 1);
if ~isstruct(config)
  return;
end
if isfield(config, 'baseSeed')
  lineList(end + 1, 1) = "<b>Base seed</b>: <code>" + localHtmlEscape(localScalarToString(config.baseSeed)) + "</code>";
end
if isfield(config, 'numRepeat')
  lineList(end + 1, 1) = "<b>Repeats</b>: <code>" + localHtmlEscape(localScalarToString(config.numRepeat)) + "</code>";
end
if isfield(config, 'snrDb')
  lineList(end + 1, 1) = "<b>SNR</b>: <code>" + localHtmlEscape(localScalarToString(config.snrDb)) + " dB</code>";
end
end

function textValue = localScalarToString(value)
% Convert a scalar-ish value into compact display text.
if isstring(value) || ischar(value)
  textValue = char(string(value));
elseif isnumeric(value) || islogical(value)
  if isempty(value)
    textValue = '';
  elseif isscalar(value)
    if islogical(value)
      textValue = sprintf('%d', logical(value));
    elseif isfinite(double(value)) && abs(double(value) - round(double(value))) < eps(max(1, abs(double(value))))
      textValue = sprintf('%d', round(double(value)));
    else
      textValue = sprintf('%.6g', double(value));
    end
  else
    textValue = sprintf('[%s]', strjoin(cellstr(string(value(:)')), ', '));
  end
else
  textValue = '<unprintable>';
end
end

function statusIcon = localStatusIcon(statusText)
% Return a compact unicode icon for the status field.
if strcmpi(char(statusText), 'DONE') || strcmpi(char(statusText), 'OK')
  statusIcon = "✅";
elseif strcmpi(char(statusText), 'FAILED') || strcmpi(char(statusText), 'ERROR')
  statusIcon = "❌";
else
  statusIcon = "ℹ️";
end
end

function textValue = localFormatElapsed(elapsedSec)
% Format elapsed seconds as HH:MM:SS.
elapsedSec = max(0, double(elapsedSec));
hourValue = floor(elapsedSec / 3600);
minuteValue = floor(mod(elapsedSec, 3600) / 60);
secondValue = floor(mod(elapsedSec, 60));
textValue = sprintf('%02d:%02d:%02d', hourValue, minuteValue, secondValue);
end

function stringList = localAsStringColumn(value)
% Convert input text into a string column.
if isempty(value)
  stringList = strings(0, 1);
  return;
end
stringList = string(value(:));
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
% Return struct field value when present.
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function value = localGetLogicalField(dataStruct, fieldName, defaultValue)
% Return a logical struct field value when present.
value = logical(defaultValue);
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = logical(dataStruct.(fieldName));
end
end

function textValue = localHtmlEscape(textValue)
% Escape text for Telegram HTML parse mode.
textValue = string(textValue);
textValue = replace(textValue, "&", "&amp;");
textValue = replace(textValue, "<", "&lt;");
textValue = replace(textValue, ">", "&gt;");
end
