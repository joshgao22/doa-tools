# test/common 目录说明

`test/common/` 放实验侧复用 helper。这里的 helper 服务 regression、dev、replay、scan、perf orchestration，不应长期承载正式 estimator 主核或正式 scene 语义。

新增 test common helper 前，先查本目录已有 fixture / flow / probe / report / summary / replay / scan helper，避免为单个脚本复制 checkpoint、snapshot、summary、notify 或表格预览壳。只有两个以上非 legacy 入口稳定复用、且不改变 estimator 数值路径的工程逻辑，才适合放到这里。

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
- `buildDoaDopplerStaticTransitionBundle.m`
- `evaluateDynamicSubsetBank.m`
- `shouldRunSimpleSameToothBasinEntryRescue.m`

当前约束：

- flow 决策必须使用 decision summary，只允许读取 objective、residual、health、candidate margin、DoA drift 等接收数据 / estimator 内部量；
- truth-aware evaluation summary 只用于 replay / scan 表格和离线结果分析，`truthTooth*`、`toothIdx`、`fdRefErrHz`、`angleErrDeg` 等字段不能进入 selector、rescue gate 或 polish gate；
- `runSimpleDynamicSubsetPeriodicFlow.m` 同时保存 `subsetDecisionSummaryCell / periodicDecisionSummaryCell` 与 truth-aware summary，前者用于决策，后者用于输出。

不应包含：

- 正式 objective evaluator；
- profile likelihood 主核；
- satellite geometry 语义。

如果 flow 里的逻辑开始成为正式算法，应迁移到 `estimator/helper/`。`buildDoaDopplerStaticTransitionBundle.m` 默认保持完整 SF static transition ladder；只在 replay / scan 显式传入 `bundleMode="ms-seed-only"` 时跳过 other-sat ablation 与 weight sweep，用于 MS dynamic init 诊断提速，不改变默认 flow 语义。

### `probe/`

负责 probe 点评估 helper。

典型文件：

- `evalDoaDopplerMfProbePoint.m`
- `buildMfObjectiveProbeRows.m`
- `buildMfTargetedSolveProbeRows.m`
- `buildMsCpuBankProbeRows.m`
- `buildMfSolveProbeRowFromCase.m`
- `buildMfDoaBasinEntryRowsFromCase.m`
- `buildMfPerSatProbeRows.m`
- `makeMfProbeOptVar.m` / `resolveMfProbeTruthScalar.m` / `extractMfProbeCoherenceFloor.m`
- `emptyMfObjectiveProbeRow.m` / `emptyMfSolveProbeRow.m` / `emptyMfDoaBasinEntryRow.m` / `emptyMfPerSatProbeRow.m`

`buildMfObjectiveProbeRows.m` 只负责构造 MF model 并重评估 static-seed / final / truth / mixed fixed-point objective，供 replay / scan 诊断复用；它不运行 solver、不选择 winner、不改变 estimator adoption。

`buildMfTargetedSolveProbeRows.m` 负责 replay / scan 侧的 targeted controlled MF solver rerun：每个 method 记录 baseline row，默认只有 MS bad condition 才触发 truth-DoA、static-wide 或 SS-MF seed probe，并把 MS CP-U bank 委托给 `buildMsCpuBankProbeRows.m`。`solveProbeRoute="ms-bank-only"` 对所有 MS CP-U case 跑 compact bank，用于 coverage-first 诊断和 adopted-row 上限确认；上一版 no-truth-health-only targeted bank 会漏 hidden angle-only tail，当前不再作为 replay 公开路线。上层显式传入 `solveProbeRoute="ss-parity"` 时，它只为 SS baseline parity 检查运行 SS truth / static-wide probe，用于对照 `replayMfSsMleCrbMetricDiagnose`；这些 probe 都不改变 estimator final winner。

`buildMsCpuBankProbeRows.m` 负责 MS-MF-CP-U 的 replay-only axis-cross bank probe：先用 CP-K cheap preselect 排序，再对 top candidate 做 CP-U release，用于定位 basin-entry / nuisance-rate release 风险；它不定义正式 rescue strategy，不进入 estimator 默认路径。

`buildMfSolveProbeRowFromCase.m`、`buildMfDoaBasinEntryRowsFromCase.m`、`buildMfPerSatProbeRows.m` 与 `emptyMf*Row.m` 是 probe / replay 共用的表行 primitive，只负责把已有 estimator 输出重组为诊断表，不运行新 solver，也不做 candidate adoption。`buildMfDoaBasinEntryRowsFromCase.m` 会透传 estimator 的 baseline / entry / polish / selected objective、entry / polish / selected 参数、entry center / offset / center-source 与 `entryAdoptionMode` 字段，用于 replay 侧判断 estimator-side DoA basin-entry 是未构造、未评估、外部 center 是否进入，评估后没有赢过 baseline，还是已有 truth-free candidate 本身就无法解释 trim 后 CRB gap；selected row 归因优先使用 basin-entry 诊断中的 `selectedVariant / bestTag`，避免 unknown warm-anchor 后续 `solveVariant` 覆盖 entry family。外部 center 可携带 source label，便于 replay 汇总 MS entry family 的 selected rate、selected-trim rate 与 CRB-normalized residual；当前 MS replay 默认采用 trim-focused family：只保留 `ssmf-center` 作为外部 truth-free center；`static-center`、`static-axis`、`static-diagonal` 与 `ssmf-axis` 不再默认维护。CRB/FIM variant 与 trim-kept residual 分母对照仍保留在 replay local helper 中，尚未稳定到 common report helper。

不应包含：

- 长 scan orchestration；
- 正式 estimator branch selection。

### `report/`

负责报告表、CRB bundle、代表性 replay report，以及 replay / scan 共用的 compact metric / runtime / tail 表格与命令行预览壳。

`buildMfReleaseCompareTable.m`、`buildMfSeedChainTable.m`、`buildMfPerSatProbeTable.m`、`buildMfRuntimeTable.m` 与 `buildMfRuntimeRow.m` 是 replay / scan 侧的稳定报告表 helper，只收集或比较已有 case / probe rows，不定义诊断标签、不触发新 estimator 调用。

典型文件：

- `buildDynamicCaseSummaryTable.m`
- `buildDynamicCrbSummaryTable.m`
- `buildMfMetricAggregateTable.m`

  - MF metric aggregate 的 `angleMseOver*SphericalCrb` / `fdRefMseOverCrb` 使用逐样本 `mean((error/CRB)^2)`；`RMSE/CRB` 为该值平方根。物理量 RMSE、MSE 与 CRB median 仍保留为诊断字段，但不再用 `MSE/median(CRB)^2` 作为主 MSE/CRB 口径。
- `buildMfFilteredMetricAggregateTable.m`
- `buildMfTailTable.m`
- `buildMfRuntimeAggregateTable.m`
- `buildMfTopSlowRuntimeTable.m`
- `buildMfCandidateTraceTable.m`
- `buildMfWideCoverageAggregateTable.m`
- `buildMfBankAdoptionShadowTable.m`
- `buildMfBankAdoptionShadowAggregateTable.m`
- `buildMfBankAdoptionRejectAggregateTable.m`
- `buildMfBankAdoptedCaseTable.m`
- `buildMfBankRescueOutcomeTable.m`
- `buildMfBankHealthGateTriggerOutcomeTable.m`
- `buildMfBankFamilyAggregateTable.m`
- `printMfReportTableSection.m`

不应包含：

- selection rule；
- winner adoption；
- objective 主核。

`buildMfCandidateTraceTable.m`、`buildMfWideCoverageAggregateTable.m`、`buildMfBankFamilyAggregateTable.m`、`buildMfBankAdoptionShadow*.m`、`buildMfBankAdoptionRejectAggregateTable.m`、`buildMfBankAdoptedCaseTable.m`、`buildMfBankRescueOutcomeTable.m` 与 `buildMfBankHealthGateTriggerOutcomeTable.m` 只重组已有 case / probe 表，允许 replay 与 scan 复用候选来源、objective delta、angle improvement、damage 统计、bank family 拆分、in-tooth tooth/comb 诊断字段、离线 shadow-adoption 评价、fixed-SNR rescue 统计比率、angle-or-fdRef bad-gated replay 诊断、runtime-health-gated replay 诊断、truth-free candidate-objective 健康触发、触发源分组诊断和 CRB-normalized soft-damage / healthy-adopt / bad-rescue 指标；它们不定义 MS rescue 策略，也不改变 estimator final winner。

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


### `scan/`

负责 scan 顶层 orchestration 的薄工程 helper。

典型文件：

- `printMfScanHeader.m`
- `printMfScanSection.m`
- `notifyMfScanStatus.m`
- `finalizeMfScanResult.m`

使用边界：

- 顶层 scan 是脚本，文件头固定为英文说明 + `clear; close all; clc;` + `Scan configuration`，默认参数写在配置 section；
- common/scan helper 只负责 header、section banner、best-effort Telegram 状态壳和 tmp cleanup glue；
- 数据落盘由 scan 脚本的 `Data storage` section 调用 `saveExpSnapshot` 完成；
- summary 与画图由 scan 脚本的 `Summary output and plotting` section 完成，且应只依赖 `scanData`；
- common/scan 不维护具体 scan 的 metric parser、results parser、plot helper、strategy recommendation 或 winner adoption 逻辑。

不应包含：

- 具体 scan 的实验默认参数；
- regression contract；
- 正式 estimator path；
- 长期策略规则。

### `replay/`

负责 replay batch、flow-like scan batch、header、finalize。

典型文件：

- `runSimpleDynamicFlowReplayBatch.m`
- `printMfReplayHeader.m`
- `printMfReplaySection.m`
- `notifyMfReplayStatus.m`
- `finalizeMfReplayResult.m`
- `test/common/flow/runPerfTaskGridWithCheckpoint.m`
- `test/common/flow/cleanupPerfTaskGridCheckpoint.m`

使用边界：

- 顶层 replay 是脚本，文件头固定为英文说明 + `clear; close all; clc;` + `Replay configuration`，默认参数写在配置 section；
- common/replay helper 只负责 replay / scan 的 batch 执行、header、section banner、best-effort Telegram 状态壳、progressbar 和 checkpoint glue；checkpoint run 目录清理由 `test/common/flow/cleanupPerfTaskGridCheckpoint.m` 统一处理，replay / scan 脚本不再维护私有 cleanup helper；固定单样本 replay 或短 scan 如果没有中间落盘需求，不创建 tmp，也不打印占位行；checkpoint/tmp 统一位于仓库根目录 `tmp/`；
- 数据落盘由 replay / scan 脚本的 `Data storage` section 调用 `saveExpSnapshot` 完成；
- summary 与画图由 replay / scan 脚本的 `Summary output and plotting` section 完成，且应只依赖 `replayData` 或 `scanData`；
- common/replay 不保存图片，也不维护具体 replay 的 plot helper、metric parser、policy recommendation 或 winner adoption 逻辑。

不应包含：

- 具体 replay / scan 的实验默认参数；
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

- 作用：运行 replay / scan 共用的小 MC repeat batch。
- 特点：支持外层 parfor、progressbar，以及可选 per-repeat checkpoint / resume。`progressbar` 只覆盖外层 repeat loop；并行时必须通过 `parallel.pool.DataQueue` 在 client 侧更新，worker 内不直接调用 `progressbar`。checkpoint `runState` 会记录 `runName`、`runKey`、`runDir` 与 cleanup 状态。
- checkpoint：只在 replay / scan 显式传入 `checkpointOpt.enable=true` 时启用；每个 repeat 保存一个 task 文件，manifest 校验 seed、SNR、contextOpt 和 flow signature；默认根目录为仓库根目录 `tmp/`。关闭 checkpoint 时不创建 tmp、不打印 checkpoint dir、不写空状态字段；恢复时只重跑缺失 task，成功后由 `cleanupPerfTaskGridCheckpoint` 清理 checkpoint run 目录；若脚本级 tmp 父目录已经为空，也由该 helper 删除该空父目录。
- 日志：进入 repeat loop 前默认打印紧凑 stage log 与耗时，覆盖 option resolve、shared context build、repeat mode、resume 进度和 enter parfor；这些不是 estimator verbose，不受 `optVerbose` 控制。
- 不做：硬性 pass/fail assert。


### `selectMfInToothEnvelopeSeeds.m`

- 作用：根据 envelope scout 阶段的 `tailDiagnosisTable` 选择 hard-collapse、gate-miss、easy-negative 与 fd-not-healthy 代表 seed，避免对所有 search seeds 做重 surface。
- 边界：只服务 replay seed coverage；使用离线 label 和 truth-aware summary，不能进入 runtime selector、gate、candidate adoption 或 final winner。

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
