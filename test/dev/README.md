# dev 目录说明

`test/dev/` 用于机制探索、问题定位、策略比较、代表性样本回放和长实验 orchestration。这里的脚本可以画图、打印诊断表、保存中间结果、保留失败现场，但不能替代 regression。

## 顶层代表性入口

### `doaDopplerStatDualSatUraEci.m`

- 作用：static 双星代表性实验。
- 推荐用途：检查 SS/MS、DoA/static、权重 sweep、static 锚点。
- 主要输出：estimator summary、static summary、weight sweep、geometry plot。

### `doaDopplerDynDualSatUraEci.m`

- 作用：dynamic 双星代表性实验。
- 推荐用途：检查 CP-K/CP-U、known/unknown、transition、branch、per-sat 诊断。
- 主要输出：dynamic summary、objective summary、per-sat summary、branch trace。

### `doaDopplerDynDualSatUraEciSimple.m`

- 作用：简化 dynamic flow 平台。
- 推荐用途：快速观察 subset bank、periodic refine、polish 机制。
- 主要输出：simple-flow summary、subset/final table。
- 定位：轻量机制验证平台，不是正式 dynamic 主流程替代品。

### `doaDopplerDynDualSatUraEciPerf.m`

- 作用：dynamic 长 Monte Carlo / 正式结果生产。
- 推荐用途：跑 task grid、checkpoint/resume、snapshot。
- 主要输出：summary table、snapshot、cache。

### `doaDopplerStatDualSatUraEciPerf.m`

- 作用：static perf 对照。
- 推荐用途：static RMSE 与 CRB 对比。
- 主要输出：RMSE table、CRB table、snapshot。

### `doaDopplerRankSecondSatDynTle.m`

- 作用：第二星筛选 / 排序。
- 推荐用途：选代表性 pair。
- 主要输出：pair ranking。

### `doaDopplerScreenSingleSatDynTle.m`

- 作用：单星动态筛选。
- 推荐用途：找稳定参考星或代表性单星。
- 主要输出：sat ranking。

## 子目录职责

### `replay/`

固定 seed 或小 Monte Carlo 回放问题定位链。

- 典型文件：`replayMfCombToothSurface.m`、`replayMfPeriodicVsSubsetToothSelect.m`、`replayMfSameToothHardCase.m`、`replayMfUnknownReleaseRoute.m`；
- 顶层 replay 统一是脚本，文件头部显式写默认参数；
- replay 文件头固定为英文说明 + `clear; close all; clc;` + `Replay configuration`；
- 参数调整直接写在 replay 文件头部；不再使用 workspace override；
- 运行后在 workspace 中留下 `replayData`；
- `Data storage` section 只保存轻量 `replayData`，不画图；
- `Summary output and plotting` section 只依赖 `replayData`，因此从 snapshot 恢复后可直接运行这一节；
- `saveSnapshot=true` 时只用 `saveExpSnapshot` 保存 `replayData`，后续可用 `loadExpSnapshot` 恢复；
- 不默认保存图片文件，图像由 `replayData.plotData` 或 replay 的 plotting section 重画；
- replay 特有代码规范只维护在 `test/dev/replay/README.md`，详细运行结果和 snapshot 绑定维护在 `test/dev/replay/results/`；
- 默认不作为正确性护栏；
- random rescue、polish eligible 与 fastStats/fullDiag 差异都先在 replay 中观察。

### `probe/`

局部 objective/profile 探针。

- 典型文件：`probeMfDoaProfileWithHealthyFdSurface.m`；
- 通常需要人工看图；
- 不进入 regression。

### `scan/`

参数、schedule、tooth、objective surface 和论文数值图候选扫描。

- 顶层 scan 统一是 section 化脚本，文件名以 `scan` 开头；
- scan 文件头固定为英文说明 + `clear; close all; clc;` + `Scan configuration`；
- 文件头部显式写默认参数；参数调整也直接写在文件头部，不再使用 workspace override；
- `Data storage` section 只构造并保存轻量 `scanData`；
- `Summary output and plotting` section 只依赖 `scanData`，因此 load snapshot 后可直接重跑这一节；
- `saveSnapshot=true` 时通过 `saveExpSnapshot` 只保存 `scanData`；
- 不默认保存图片文件，图像由 `scanData` 中的曲线/曲面/table 数据重画；
- 典型入口：`scanMfRegimeMapByWindow.m`、`scanMfCpIpTying.m`、`scanMfFdRefComb.m`、`scanMfTruthNeighborhood.m`、`scanMfSubsetBankCoverage.m`；
- 用于建立机制证据，比 replay 更重，不作为 regression 护栏；详细运行结果和 snapshot 绑定维护在 `test/dev/scan/results/`。

具体文件用途、推荐观察对象和旧文件替换关系见 `test/dev/scan/README.md`。

### `strategy/`

策略比较和已证伪方向记录。

- 典型文件：`compareMfSubsetPeriodicRefineStrategies.m`；
- 决策文档：`blockedDynamicDirections.md`。

### `trace/`

branch、constraint、outer-start 的文本 trace。

- 典型文件：`traceMfUnknownDoaConstraint.m`；
- 用于轻量排查，不替代 replay。

## replay / probe / scan / strategy 的使用边界

### replay

- 关注点：一个固定现象是否可复现，在哪些 seed / 分支出现。
- 是否小 MC：可以，通常 2–20 次。
- 是否画图：可以。
- 是否作为正确性护栏：否。

### probe

- 关注点：一个局部机制，例如 profile 或 surface。
- 是否小 MC：通常不做。
- 是否画图：经常画图。
- 是否作为正确性护栏：否。

### scan

- 关注点：某个维度系统变化时机制如何变化。
- 是否小 MC：可选，但通常较重。
- 是否画图：通常画曲线或曲面。
- 是否作为正确性护栏：否。

### strategy

- 关注点：比较不同流程策略或记录已证伪路线。
- 是否小 MC：可以。
- 是否画图：可选。
- 是否作为正确性护栏：否，除非结论迁移为 regression。

## 什么时候迁移到 regression

一个 dev 结论只有满足以下条件，才迁移到 regression：

- 结论已经固定为代码契约；
- 不需要人工看图或看表；
- 可以写明确 assert；
- seed 固定且 runtime 可接受；
- 不依赖临时 cache；
- 不复制正式主逻辑；
- 失败能指向具体契约。

## 长任务与保存规则

- 运行时临时存储放 `tmp/<scriptName>/<runKey>/`；
- 最终结果放 `test/data/cache/<taskType>/`，例如 `replay/`、`scan/`、`perf/`；
- 正常结束默认 cleanup；checkpoint task grid 的 run 目录和空父目录由 common checkpoint cleanup helper 统一删除，不在 replay/scan 脚本里复制私有 cleanup helper；
- replay 保存 `replayData`，scan 保存 `scanData`；两者都应只保存 summary、table、plot data、config 和 meta；详细人工结果分别记录在 `test/dev/replay/results/` 和 `test/dev/scan/results/`；
- 长实验优先通过 snapshot 保存 slim result、summary、config、meta；
- 不默认保存 `rxSigCell`、完整 `sceneSeq`、完整 fixture cache 或全量 objective map。
