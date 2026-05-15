# test/common/scan

本目录放 scan 顶层 orchestration 可复用的薄 helper，不放正式 estimator 数值路径，也不维护具体 scan 的结果解释逻辑。

## 文件职责

- `printMfScanHeader.m`：统一打印 scan 头部配置，包括 seed / SNR / frame grid、snapshot、checkpoint、Telegram 和 run dir。
- `printMfScanSection.m`：统一打印 summary section banner，可选直接 `disp` 一个轻量 table / struct。
- `notifyMfScanStatus.m`：统一 HTML Telegram 状态壳，内部调用 `utils/io/notifyTelegram.m`，只负责 DONE / FAILED、耗时、snapshot、checkpoint 和调用者传入的 HTML-ready metric lines。
- `finalizeMfScanResult.m`：成功构造轻量 `scanData` 后清理 scan tmp 目录；短 scan 可传空 run dir，不创建也不清理 tmp。
- `buildMfScanCheckpointSummaryTable.m`：把一个或多个 checkpoint run state 转为轻量 summary table，供 `scanData` 保存和重载 summary 使用。
- `cleanupMfScanCheckpointRuns.m`：在 scan 顶层确认安全后批量清理已完成 checkpoint run，并返回固定字段 cleanup report。

## 边界

可以放进本目录：

- header / section banner 打印；
- best-effort 通知壳；
- checkpoint / tmp cleanup glue；
- 不改变 estimator 数值路径的 scan 组织 helper。

不要放进本目录：

- formal objective / residual / branch solve / winner adoption；
- truth-aware runtime selector、gate 或 adoption；
- scan-specific metric parser、policy recommendation、candidate family ranking；
- 需要具体 scan 内部状态的一次性 plot helper。

具体 scan 的指标行应在对应 scan 本地构造，再传给 `notifyMfScanStatus`。common helper 不解析 `scanData` 的具体表格字段，避免形成第二套 results parser。
