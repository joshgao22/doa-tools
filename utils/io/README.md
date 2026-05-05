# utils/io

本目录放通用 I/O 与运行辅助工具，不包含 DoA-Doppler 模型、estimator 数值路径或 replay / scan 专属 summary 逻辑。

## 文件职责

- `saveExpSnapshot.m`：保存轻量实验 snapshot。
- `loadExpSnapshot.m`：恢复轻量实验 snapshot。
- `cleanupRunArtifacts.m`：清理 cache / tmp 运行产物。
- `notifyTelegram.m`：通过 Telegram bot 发送可选运行通知。

## notifyTelegram

`notifyTelegram` 用于从 MATLAB 向 Telegram bot 发送运行状态消息，主要服务长时间 replay、scan、Monte Carlo 或 perf smoke。

推荐用途：

- 脚本正常结束后发送 `DONE` 通知；
- `catch ME` 后发送 `FAILED` 通知并继续 `rethrow(ME)`；
- snapshot 保存后附带 snapshot 路径；
- 只发送少量关键 summary，不替代 `results/` 文档。

使用边界：

- 通知失败只应 `warning`，不得影响仿真主流程；
- 不得在仓库中提交 bot token、chat_id 或本机私有配置；
- 不得在 `estimator/`、objective、residual、CRB 或正式算法 helper 中调用；
- replay / scan 可在顶层脚本或 orchestration wrapper 中调用，保持通知逻辑旁路；
- 不在 `notifyTelegram` 内维护特定 replay / scan 的指标解释、results 解析或 winner adoption 逻辑。

配置方式：

1. 优先读取环境变量 `TG_BOT_TOKEN` 与 `TG_CHAT_ID`。
2. GUI MATLAB 可使用 MATLAB `prefdir` 下的私有 `telegramNotifyConfig.mat`。
3. 私有配置文件只保留在本机，不进入 git。

GUI MATLAB 的本机配置示例：

```matlab
tg = struct();
tg.BotToken = "your-bot-token";
tg.ChatId = "your-chat-id";
save(fullfile(prefdir, "telegramNotifyConfig.mat"), "tg");
```

最小测试：

```matlab
notifyTelegram("MATLAB TEST", "Telegram notification works.", "HTML");
```

HTML 模式适合成功通知；失败通知若直接拼接任意 error message，应先确认文本已经 HTML escape，或改用 plain text，避免通知本身因为特殊字符发送失败。
