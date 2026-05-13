# test 目录说明

`test/` 是当前仓库里最常用的实验、回归和排障入口。这个目录不直接承载正式算法主核；正式算法应在 `estimator/`、`satellite/`、`performance/` 中实现，`test/` 只负责组织、验证、回放和报告。

## 文档边界

- 本 README 负责说明 test 入口怎么选、各子目录职责和统一落盘边界。
- replay / scan 的脚本格式和文件索引分别放在 `test/dev/replay/README.md` 与 `test/dev/scan/README.md`。
- 具体运行结果、长表格、观察现象和 snapshot 绑定放在对应 `results/` 文档。
- 大 `.mat` snapshot 的保存、恢复、清理入口只以 `test/data/cache/README.md` 为单一索引；本 README 不复制调用细节。
- 排障主记录只摘取影响当前优先级的机制结论，不替代本 README 的入口说明。

## 如何选择入口

### 我想快速验证代码没回退

- 去：`test/regression/`
- 入口：`runRegressionQuick`
- 适用场景：改 estimator、flow、scene、summary、report 后。

### 我想验证 dynamic flow 分支规则

- 去：`test/regression/branch/`
- 典型入口：`regressionMfSubsetSelectNoTruthLeak`
- 适用场景：改 subset ranking、final selection、warm-anchor gate、polish trigger。

### 我想验证完整代表性机制

- 去：`test/regression/pipeline/`
- 典型入口：`regressionMfDoaProfileWithHealthyFd`
- 适用场景：改完整 fixture -> estimator -> summary 链路。

### 我想验证 checkpoint / snapshot / cleanup

- 去：`test/regression/perf/`
- 入口：`runRegressionPerfSmoke`
- 适用场景：改长任务工程链路。

### 我想复现一个固定问题或小 MC 现象

- 去：`test/dev/replay/`
- 典型入口：`replayMfSameToothHardCase`
- 特点：可画图、可保存 slim snapshot、默认不作为正确性护栏。

### 我想看 objective 曲线、profile 或局部曲面

- 去：`test/dev/probe/`
- 典型入口：`probeMfDoaProfileWithHealthyFdSurface`
- 特点：通常需要人工看图，不进 regression。

### 我想扫 frame schedule、comb、truth-neighborhood

- 去：`test/dev/scan/`
- 典型入口：`scanMfFdRefComb`、`scanMfTauSchedule`、`scanMfTruthNeighborhood`
- 特点：比 replay 更重，用于建立机制证据；入口文件统一以 `scan` 开头。

### 我想重画论文 / 汇报图

- 去：`test/paper/`
- 典型入口：`plotMfRegimeMapByWindowPaper`、`plotMfKnownUnknownInformationLossPaper`
- 特点：读取已有轻量 snapshot，整理成 paper-facing 图；不重新运行 estimator、scan 或 replay。

### 我想比较 flow 策略或记录不通方向

- 去：`test/dev/strategy/`
- 典型入口：`compareMfSubsetPeriodicRefineStrategies`
- 文档：`blockedDynamicDirections.md`

### 我想复用 fixture、summary、flow helper

- 去：`test/common/`
- 注意：这里的 helper 不应长期承载正式 estimator 主逻辑。

## 顶层 dev 文件

### `doaDopplerStatDualSatUraEci.m`

- 作用：static 双星代表性样本。
- 主要检查：SS/MS、DoA/static、权重 sweep、static 锚点。
- 主要输出：estimator summary、static summary、weight sweep。

### `doaDopplerDynDualSatUraEci.m`

- 作用：dynamic 双星代表性样本。
- 主要检查：CP-K/CP-U、known/unknown、transition、branch、per-sat 诊断。
- 主要输出：dynamic summary、objective summary、per-sat summary、branch trace。

### `doaDopplerDynDualSatUraEciPerf.m`

- 作用：dynamic 长 Monte Carlo / 正式数据生产入口。
- 主要输出：summary table、snapshot、checkpoint。
- 注意：长任务结果优先通过 snapshot 保存，不散落 ad-hoc save。

### `doaDopplerDynDualSatUraEciSimple.m`

- 作用：简化 dynamic flow 平台。
- 用途：轻量观察 subset bank、periodic refine、polish 机制。
- 定位：机制验证平台，不是正式 dynamic 主流程替代品。

### `doaDopplerStatDualSatUraEciPerf.m`

- 作用：static perf 对照。
- 主要输出：RMSE、CRB 对比、cache。

## 子目录职责

### `regression/`

自动 pass/fail 护栏。

- 只放已稳定的契约；
- 不放需要人工看图的脚本；
- 不放未修好的 hard case；
- 不放长 Monte Carlo 或策略探索。

### `dev/replay/`

固定 seed 或小 Monte Carlo 的问题回放。

- 可打印 trace；
- 可画图；
- 可保存 slim snapshot；
- 默认不 `error`，不作为 quick regression。

### `dev/probe/`

局部机制探针。

- 常用于 objective surface、profile curve、truth/final objective gap；
- 需要人工判断时不要迁移到 regression。

### `dev/scan/`

系统扫描。

- 适合 comb scan、schedule scan、tooth scan、参数 sweep；
- 比 replay 更重，通常不作为日常入口。

### `dev/strategy/`

策略比较与决策记录。

- 比较 subset bank、periodic/non-periodic selection、polish gate 等规则；
- `blockedDynamicDirections.md` 记录已证伪或不建议回流主路径的方向。

### `paper/`

论文 / 汇报图重绘入口。

- 读取 `test/data/cache/` 下的轻量 snapshot；
- 将 scan / replay 已有结果整理为 paper-facing 图或简表；
- 不重新运行 estimator、scan 或 replay，也不维护第二套数值逻辑。

### `dev/trace/`

文本 trace。

- 用于 branch、constraint、outer-start 的精简输出；
- 不替代 replay，也不替代 regression。

### `common/`

多脚本复用 helper。

- fixture/context 构造；
- flow orchestration；
- summary/report/plotData；
- replay batch helper；
- scan 运行 / 汇总 helper；不维护统一配置合并入口；
- regression runner 工具。

## replay / probe / scan / strategy 的区别

### replay

- 关注：固定现象是否可复现，在哪些 seed / 分支出现。
- 可做小 MC，通常 2–20 次。
- 可以画图，可以保存 snapshot。
- 不作为正确性护栏。

### probe

- 关注：一个局部机制，例如 profile 或 surface。
- 通常不做 MC。
- 经常画图。
- 不作为正确性护栏。

### scan

- 关注：某个参数、schedule、tooth、窗口长度如何系统影响结果。
- 通常较重。
- 可并行，可保存 cache。
- 不作为 quick regression。

### strategy

- 关注：规则比较和路线决策。
- 通常输出对比表。
- 不把策略优劣固化成 regression，除非已经变成明确契约。

## 临时文件与结果保存

### 临时文件

- 默认放：`tmp/<scriptName>/<runKey>/`
- 正常结束：默认 cleanup。
- 失败、显式调试或画图：可以保留现场。

### 最终结果

- 默认放：`test/data/cache/<taskType>/`，例如 `replay/`、`scan/`、`perf/`。
- 文件名建议：`脚本名_yyyymmdd-HHMMSS.mat`。
- 保存最终 snapshot 用 `saveExpSnapshot`；恢复用 `loadExpSnapshot`；清理 cache / tmp 运行产物用 `cleanupRunArtifacts`；checkpoint run 目录清理用 `cleanupPerfTaskGridCheckpoint`。
- 上述入口的路径、调用示例和安全边界统一维护在 `test/data/cache/README.md`。
- replay 优先保存 `replayData`，scan 优先保存 `scanData`；大 MAT 只放 cache，人工分析放对应 `results/` 文档。
- 优先保存 slim result / summary / config / plotData / meta。
- 不默认保存图片文件；图片由保存的数据重画。
- 不默认保存大中间量，例如 `rxSigCell`、`sceneSeq`、`fixtureCell`、全量 objective map。

## 避免重复验证

- subset no-truth-leak 只由 `regressionMfSubsetSelectNoTruthLeak` 验证。
- final winner 规则只由 `regressionMfUnknownFinalSelectionRules` 验证。
- fast subset escalation 只由 `regressionMfFastSubsetEscalation` 验证。
- warm-anchor parfor 默认 gate 只由 `regressionMfWarmAnchorParforGate` 验证。
- same-tooth hard case 在修复前放 replay，不放 pipeline regression。
- CP/IP time-axis 由 invariant regression 验证，不在 replay 里重复 assert。
