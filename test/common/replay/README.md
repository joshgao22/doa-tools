# test/common/replay

本目录放 replay / scan 顶层 orchestration 可复用的薄 helper，不放正式 estimator 数值路径，也不维护具体 replay 的结果解释逻辑。

## 文件职责

- `runSimpleDynamicFlowReplayBatch.m`：小 MC replay / scan batch runner，支持外层 parfor、progressbar、可选 per-repeat checkpoint / resume。
- `printMfReplayHeader.m`：统一打印 replay 头部配置，包括 repeat、SNR、seed、snapshot、checkpoint、Telegram、run dir。
- `printMfReplaySection.m`：统一打印 summary section banner，可选直接 `disp` 一个轻量 table / struct。
- `notifyMfReplayStatus.m`：统一 HTML Telegram 状态壳，内部调用 `utils/io/notifyTelegram.m`，只负责 DONE / FAILED、耗时、snapshot、checkpoint 和调用者传入的 HTML-ready metric lines。
- `finalizeMfReplayResult.m`：成功构造轻量 `replayData` 后清理 replay tmp 目录。
- `selectMfInToothEnvelopeSeeds.m`：in-tooth envelope replay 的 seed selection 复用逻辑。

## 边界

可以放进本目录：

- header / section banner 打印；
- best-effort 通知壳；
- 外层 batch / checkpoint / resume glue；
- 不改变 estimator 数值路径的 replay 组织 helper。

不要放进本目录：

- formal objective / residual / branch solve / winner adoption；
- truth-aware runtime selector、gate 或 adoption；
- replay-specific metric parser、policy recommendation、candidate family ranking；
- 需要具体 replay 内部状态的一次性 plot helper。

具体 replay 的指标行应在对应 replay 本地构造，再传给 `notifyMfReplayStatus`。common helper 不解析 `replayData` 的具体表格字段，避免形成第二套 results parser。
