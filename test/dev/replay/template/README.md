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
