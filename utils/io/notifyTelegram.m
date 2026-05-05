function ok = notifyTelegram(titleText, messageText, parseMode)
% Send a Telegram notification from MATLAB.
%
% parseMode:
%   ""       : plain text
%   "HTML"   : Telegram HTML formatting
%   "MarkdownV2" : Telegram MarkdownV2 formatting

if nargin < 3
    parseMode = "";
end

ok = false;

cfg = localLoadTelegramConfig();

if strlength(cfg.BotToken) == 0 || strlength(cfg.ChatId) == 0
    warning("Telegram notification skipped: missing bot token or chat id.");
    return;
end

fullText = sprintf("[%s]\n%s", string(titleText), string(messageText));
chunkList = localSplitText(string(fullText), 3500);

url = "https://api.telegram.org/bot" + cfg.BotToken + "/sendMessage";

options = weboptions( ...
    "MediaType", "application/json", ...
    "Timeout", 15);

try
    for idx = 1:numel(chunkList)
        payload = struct();
        payload.chat_id = char(cfg.ChatId);
        payload.text = char(chunkList(idx));
        payload.disable_web_page_preview = true;

        if strlength(string(parseMode)) > 0
            payload.parse_mode = char(parseMode);
        end

        webwrite(url, payload, options);
    end

    ok = true;
catch ME
    warning("Telegram notification failed: %s", ME.message);
end
end

function cfg = localLoadTelegramConfig()
% Load Telegram configuration from environment or local preference file.

cfg = struct();
cfg.BotToken = string(getenv("TG_BOT_TOKEN"));
cfg.ChatId = string(getenv("TG_CHAT_ID"));

if strlength(cfg.BotToken) > 0 && strlength(cfg.ChatId) > 0
    return;
end

configFile = fullfile(prefdir, "telegramNotifyConfig.mat");

if ~isfile(configFile)
    return;
end

loadedData = load(configFile);

if isfield(loadedData, "tg")
    tg = loadedData.tg;

    if strlength(cfg.BotToken) == 0 && isfield(tg, "BotToken")
        cfg.BotToken = string(tg.BotToken);
    end

    if strlength(cfg.ChatId) == 0 && isfield(tg, "ChatId")
        cfg.ChatId = string(tg.ChatId);
    end
end
end

function chunkList = localSplitText(textValue, maxLen)
% Split long text into Telegram-safe chunks.

textValue = string(textValue);

if strlength(textValue) <= maxLen
    chunkList = textValue;
    return;
end

chunkList = strings(0, 1);
remainingText = textValue;

while strlength(remainingText) > maxLen
    chunkList(end + 1, 1) = extractBefore(remainingText, maxLen + 1);
    remainingText = extractAfter(remainingText, maxLen);
end

if strlength(remainingText) > 0
    chunkList(end + 1, 1) = remainingText;
end
end

