# test/common/replay

本目录放 replay / scan 顶层 orchestration 可复用的薄 helper，不放正式 estimator 数值路径，也不维护具体 replay 的结果解释逻辑。

## 文件职责

- `runSimpleDynamicFlowReplayBatch.m`：小 MC replay / scan batch runner，支持外层 parfor、progressbar、可选 per-repeat checkpoint / resume。
- `printMfReplayHeader.m`：统一打印 replay 头部配置，包括 repeat、SNR、seed、snapshot、checkpoint、Telegram、run dir。连续 `seedList` 使用 `a:b (N seeds)`，非连续 seed list 完整打印，便于手动 seed bank 复查。
- `formatMfReplaySeedList.m`：统一 seed list 日志格式化；只负责 replay 日志展示，不参与 run-key、checkpoint 或数值路径。
- `printMfReplaySection.m`：统一打印 summary section banner，可选直接 `disp` 一个轻量 table / struct。
- `notifyMfReplayStatus.m`：统一 HTML Telegram 状态壳，内部调用 `utils/io/notifyTelegram.m`，只负责 DONE / FAILED、耗时、snapshot、checkpoint 和调用者传入的 HTML-ready metric lines。
- `buildMfReplayCheckpointOpt.m`：统一 replay checkpoint 的 `runName / runKey / outputRoot / resume / runDir / meta` 外壳；支持 replay 本地传入 stage-specific `runName` 与 `useParfor`，具体 run-key signature 与 meta 字段仍由 replay 本地决定。
- `buildMfReplayCheckpointSummary.m`：统一从单阶段或多阶段 checkpoint runner 返回值提取轻量 `runName / runKey / runDir / numTask / numDone / isComplete / cleanedOnSuccess` summary。
- `dispMfReplayTablePreview.m`：统一长表命令行预览，只打印首尾若干行，完整表仍保存在 `replayData`。
- `finalizeMfReplayResult.m`：成功构造轻量 `replayData` 后清理 replay tmp 目录。
- `selectMfInToothEnvelopeSeeds.m`：in-tooth envelope replay 的 seed selection 复用逻辑。

## 边界

可以放进本目录：

- header / section banner 打印；
- best-effort 通知壳；
- checkpoint option / summary / table preview glue；
- 外层 batch / checkpoint / resume glue；
- 不改变 estimator 数值路径的 replay 组织 helper。

不要放进本目录：

- formal objective / residual / branch solve / winner adoption；
- truth-aware runtime selector、gate 或 adoption；
- replay-specific metric parser、policy recommendation、candidate family ranking；
- 需要具体 replay 内部状态的一次性 plot helper。

具体 replay 的指标行应在对应 replay 本地构造，再传给 `notifyMfReplayStatus`。common helper 不解析 `replayData` 的具体表格字段，避免形成第二套 results parser。
