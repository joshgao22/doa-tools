# replay template

本目录只放 copy-only replay 模板，不放正式可运行 replay 入口。

新增 replay 时，先复制：

```text
test/dev/replay/template/replayTemplate.m
```

到：

```text
test/dev/replay/replayYourTopic.m
```

然后再改文件名、`replayName`、头部配置、batch 逻辑、summary 和 plot。模板只规范 section、checkpoint/snapshot/Telegram 壳和打印风格，不提供通用实验框架，也不承载 replay-specific 结果解释。

使用边界：

- `.m` 文件头注释与 inline comment 使用英文；
- 默认参数直接写在 `Replay configuration`；
- snapshot 只保存轻量 `replayData`；
- Telegram 指标行由具体 replay 本地构造，common helper 只负责 HTML 通知壳；
- 不把 truth-aware selector、gate、adoption 或 replay-specific policy 解析抽到 common helper。

## Telegram 通知格式模板

replay 模板采用“固定外壳 + 弹性指标行”的 HTML 通知格式。固定外壳由 `test/common/replay/notifyMfReplayStatus.m` 生成，包含状态、脚本名、耗时、常用配置、snapshot 文件/目录、失败错误和 location。各 replay 只在本地构造 `metricLineList` 与 `commentLineList`。

指标行要求：

- 只放 3–6 行最关键 summary，不替代 results 文档；
- 每行尽量在手机上一屏内可读，避免长句堆叠；
- `<code>` 只用于短数字、短 tag 或短状态，不包裹完整 snapshot 路径、长 candidate family 名或长错误信息；
- 动态文本若不是刻意写入 HTML tag，必须先做 HTML escape；
- `commentLineList` 写本次 replay 的结论、recommendation 或下一步判断，不只写 `completed`。

通知仍是 best-effort 旁路：失败只 `warning`，不得改变 replay 数值路径、`replayData`、snapshot、表格或 results 口径。
