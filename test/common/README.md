# test/common 目录说明

`test/common/` 放实验侧复用 helper。这里的 helper 服务 regression、dev、replay、scan、perf orchestration，不应长期承载正式 estimator 主核或正式 scene 语义。

## 子目录职责

### `case/`

负责 case result、truth view、summary table。

典型文件：

- `buildDoaDopplerCaseResult.m`
- `buildDoaDopplerEstView.m`
- `buildDoaDopplerCaseSummaryTable.m`

不应包含：

- branch selection；
- profile likelihood；
- estimator 主逻辑。

### `fixture/`

负责可复用 fixture/context 构造。

典型文件：

- `buildDynamicDualSatEciContext.m`
- `buildDynamicSubsetFixtureBank.m`
- `buildMfUnknownReleaseFixture.m`

`buildDynamicDualSatEciContext.m` 负责动态双星 dev/replay/scan/regression 的统一代表性场景入口。它内置一组固定默认仿真参数；只有当某个 replay / scan 明确要改变场景或 waveform 时，才由上层脚本显式传入覆盖项。

不应包含：

- formal scene semantics 的替代实现；
- estimator 主核；
- strategy selection 规则。

### `flow/`

负责 test/dev 层 orchestration flow。

典型文件：

- `runSimpleDynamicSubsetPeriodicFlow.m`
- `evaluateDynamicSubsetBank.m`

不应包含：

- 正式 objective evaluator；
- profile likelihood 主核；
- satellite geometry 语义。

如果 flow 里的逻辑开始成为正式算法，应迁移到 `estimator/helper/`。

### `probe/`

负责 probe 点评估 helper。

典型文件：

- `evalDoaDopplerMfProbePoint.m`

不应包含：

- 长 scan orchestration；
- 正式 estimator branch selection。

### `report/`

负责报告表、CRB bundle、代表性 replay report。

典型文件：

- `buildDynamicCaseSummaryTable.m`
- `buildDynamicCrbSummaryTable.m`

不应包含：

- selection rule；
- winner adoption；
- objective 主核。

### `summary/`

负责 compact summary helper。

典型文件：

- `summarizeDynamicEstimatorDebug.m`
- `summarizeDynamicMultiStart.m`

不应包含：

- 大矩阵 dump；
- 全量 trace cache；
- 改变 estimator 结果的逻辑。

### `plot/`

负责多脚本复用绘图。

典型文件：

- `plotDoaDopplerGeometryComparison.m`

不应包含：

- 单脚本一次性图；
- 需要特定 replay 内部状态的临时图。

### `replay/`

负责 replay batch、header、finalize。

典型文件：

- `runSimpleDynamicFlowReplayBatch.m`
- `printMfReplayHeader.m`
- `finalizeMfReplayResult.m`
- `test/common/flow/runPerfTaskGridWithCheckpoint.m`
- `test/common/flow/cleanupPerfTaskGridCheckpoint.m`

使用边界：

- 顶层 replay 是脚本，文件头固定为英文说明 + `clear; close all; clc;` + `Replay configuration`，默认参数写在配置 section；
- common/replay helper 只负责 batch 执行、header、progressbar 和 checkpoint glue；checkpoint run 目录清理由 `test/common/flow/cleanupPerfTaskGridCheckpoint.m` 统一处理，replay 脚本不再维护私有 cleanup helper；固定单样本 replay 如果没有中间落盘需求，不创建 tmp，也不打印占位行；checkpoint/tmp 统一位于仓库根目录 `tmp/`；
- 数据落盘由 replay 脚本的 `Data storage` section 调用 `saveExpSnapshot` 完成；
- summary 与画图由 replay 脚本的 `Summary output and plotting` section 完成，且应只依赖 `replayData`；
- common/replay 不保存图片，也不维护具体 replay 的 plot helper。

不应包含：

- 具体 replay 的实验默认参数；
- regression contract；
- 正式 estimator path；
- 长期策略规则。

### `util/`

负责 test 通用工具。

典型文件：

- `runRegressionSuite.m`
- `parseRegressionCaseOpt.m`

不应包含：

- 模型语义 helper；
- estimator-specific objective helper。

## 什么时候应该迁移出 test/common

### 迁移到 `estimator/helper/`

当内容变成以下正式估计器逻辑时：

- profile likelihood；
- objective evaluation；
- 参数打包；
- branch solving；
- winner adoption；
- warm-anchor 主流程。

### 迁移到 `satellite/scene/`

当内容变成以下正式场景语义时：

- reference-sat；
- scene 裁剪；
- user-state；
- Doppler geometry；
- reference Doppler state。

### 迁移到 `satellite/signal/`

当内容变成以下正式信号逻辑时：

- pilot 生成；
- rx signal 生成；
- snapshot 构造；
- frame/sat 信号裁剪。

### 迁移到 `performance/`

当内容属于：

- CRB；
- FIM；
- EFIM；
- Jacobian；
- information loss。

### 迁移到 `utils/`

当内容是多处复用且与模型无关的通用工具。

## 公共落盘与清理规范索引

- 最终 snapshot 的保存、恢复和 cache 删除规范见 `test/data/cache/README.md`。
- 本目录只说明 replay / scan / perf orchestration 中的 checkpoint、batch 执行和中间目录清理。
- replay / scan 脚本不要复制 `saveExpSnapshot`、`loadExpSnapshot`、`cleanupRunArtifacts` 或 checkpoint cleanup 的使用说明；只在各自 `Data storage` / `Summary output and plotting` section 给出局部用法。
- checkpoint run 目录清理由 `cleanupPerfTaskGridCheckpoint` 负责；普通 tmp/cache 运行产物清理由 `cleanupRunArtifacts` 负责。

## replay helper

### `runSimpleDynamicFlowReplayBatch.m`

- 作用：运行 replay 的小 MC repeat batch。
- 特点：支持外层 parfor、progressbar，以及可选 per-repeat checkpoint / resume。`progressbar` 只覆盖外层 repeat loop；并行时必须通过 `parallel.pool.DataQueue` 在 client 侧更新，worker 内不直接调用 `progressbar`。checkpoint `runState` 会记录 `runName`、`runKey`、`runDir` 与 cleanup 状态。
- checkpoint：只在 replay 显式传入 `checkpointOpt.enable=true` 时启用；每个 repeat 保存一个 task 文件，manifest 校验 seed、SNR、contextOpt 和 flow signature；默认根目录为仓库根目录 `tmp/`。关闭 checkpoint 时不创建 tmp、不打印 checkpoint dir、不写空状态字段；恢复时只重跑缺失 task，成功后由 `cleanupPerfTaskGridCheckpoint` 清理 checkpoint run 目录；若脚本级 tmp 父目录已经为空，也由该 helper 删除该空父目录。
- 日志：进入 repeat loop 前默认打印紧凑 stage log 与耗时，覆盖 option resolve、shared context build、repeat mode、resume 进度和 enter parfor；这些不是 estimator verbose，不受 `optVerbose` 控制。
- 不做：硬性 pass/fail assert。

### `printMfReplayHeader.m`

- 作用：打印 replay 配置、progress、保存策略；只有确实创建临时目录时才打印 `run dir`，没有临时目录时不打印占位行。
- 目的：让长 replay 能看到运行到哪一步。

### `finalizeMfReplayResult.m`

- 注意：snapshot 保存由各 replay 脚本通过 `saveExpSnapshot(includeVars={'replayData'})` 完成。

### `runPerfTaskGridWithCheckpoint.m` 与 `cleanupPerfTaskGridCheckpoint.m`

- 作用：前者负责 per-task checkpoint/resume，后者负责删除已完成 run 目录，并在脚本级 tmp 父目录为空时同步删除父目录。
- 使用边界：checkpoint runner 可以通过 `cleanupOnSuccess=true` 自动清理；如果脚本需要先组装 `replayData` 或保存 snapshot，则在 Data storage 结束后显式调用 `cleanupPerfTaskGridCheckpoint(runState, ...)`。
- 放置原因：runDir、runKey、taskDir 和 manifest 结构由 checkpoint runner 定义，清理逻辑不应复制到具体 replay 脚本 local helper。

## 避免 common 膨胀的规则

- 一个 helper 只服务一个脚本时，先保留 local；
- 两个以上脚本复用时，再提升到 `test/common/`；
- 一旦 helper 变成正式算法逻辑，应迁移到正式目录；
- 不为小改动提前建立过重框架；
- common 不维护第二套 estimator。


## `common/scan/`

`common/scan/` 暂不维护统一配置 resolver。scan 参数差异较大，默认配置和临时调整都应写在对应 scan 脚本头部。若未来出现两个以上 scan 复用同一段非参数配置逻辑，再按职责新增具体 helper。


## MATLAB 实现细节约束

- common helper 和脚本 local helper 都必须有最少一句英文函数注释，说明职责；没有职责说明的短 helper 优先内联或删除。

- 公共 helper 不保留未使用参数、dummy `~` 参数、空 `runDir` 兼容壳或只为旧脚本服务的无效字段。旧入口若仍被正式引用，可以保留 wrapper；否则直接删。
- common helper 内部若使用 `parfor`，必须让 grid / case / candidate 输入以 sliced 形式进入 worker；不要让短 grid 因二级下标访问变成 broadcast。
- common helper 内部若维护 progressbar，必须在 loop 外 reset / end；`parfor` worker 只发送 `DataQueue` 消息，不能直接调用 `progressbar`。
- checkpointed progressbar 只 reset 未完成 task 数，不能对已完成 task 伪造 advance；resume 状态用日志和 checkpoint summary 表达。
- helper 中的日志时间使用 `datetime('now', 'Format', ...)`，不要新增 `datestr(now, ...)`。
- common 只抽真正复用的 batch、summary、fixture、plot、report 逻辑；一个脚本私有的短 glue 不进 common。
- 删除冗余 helper 或分支时，必须确认没有改变 estimator 主核、flow selection、reference-sat 语义、输出字段契约或 snapshot 轻量保存口径。
