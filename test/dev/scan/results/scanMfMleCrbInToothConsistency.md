# scanMfMleCrbInToothConsistency 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-gate / SS paper-facing candidate / MS needs replay diagnosis` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-003909.mat` |
| 当前一句话结论 | `SS-MF` 在 auto init + in-tooth 条件下，local-consistent trimmed 子样本已经接近 CRB 尺度；`MS-MF` 当前不能稳定进入 multi-sat local basin，暂不适合作为正文 MS-vs-CRB 主图。 |
| 论文图定位 | `SS-MF`: main-figure candidate after metric wording check；`MS-MF`: internal diagnostic / replay input，不作为当前主图。 |
| 决策影响 | 当前 scan 版本可以固定为统计筛查入口；下一步转 replay 排查 MS 的 same-tooth tail、non-ref coherence / fdRef near-boundary branch 和 CP-U fdRate-unresolved。 |
| 下一步动作 | 新增或运行 `replayMfSsMleCrbMetricDiagnose` 与 `replayMfMsMleCrbMetricDiagnose`；SS 用于确认 angle CRB metric 与 trimmed conditional 口径，MS 用于定位 auto-init / warm-anchor / coherence / boundary 失败。 |
| 禁止误用 | 不要把 `trimmedCore` 曲线解释为 unconditional estimator efficiency；不要把 `trimmed RMSE/CRB < 1` 写成超过 CRB；不要用当前 MS 的极低 trimmed 样本率声称 MS-MF 已接近 CRB；不要把 in-tooth oracle range 当成 runtime tooth selector。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfMleCrbInToothConsistency.m`
- 结果文档：`test/dev/scan/results/scanMfMleCrbInToothConsistency.md`
- scan 类型：`paper-facing curve / local consistency screening / engineering diagnostic gate`
- 主要问题：在 coarse Doppler compensation 已将参考星 `fdRef` 限制到同一 Doppler tooth 内时，`SS/MS-MF-CP-K/U` 的 dynamic local MLE 是否能在 resolved/local 条件下与对应 CRB 对齐？
- 扫描对象：`SNR=-15:5:10 dB`，`P=10`，`SS/MS`，`CP-K/CP-U`，`auto` internal MF initializer，`full/resolved/core/trimmed` 四层统计。
- 不覆盖范围：不验证 full-flow tooth acquisition；不验证 subset tooth selection；不比较 rescue / ordinary-wide / flow-like gate；不下结论 estimator 默认路径已经修复；不形成 regression 契约。
- truth 使用口径：truth 只用于 offline evaluation、CRB、in-tooth oracle local range 与 error / tail label；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：`SS-MF` 可作为局部一致性论文候选；`MS-MF` 当前仅为 diagnostic，需 replay 修复或解释后再回到 scan。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `Doppler-aided / in-tooth` | 用 `truth ± 0.25 tooth` 的 `fdRef` local box 表示前端 coarse Doppler compensation 已消除全局 tooth ambiguity。 | Oracle / eval only | 验证 tooth-resolved local MLE 与 CRB 是否同尺度。 | 不能解释为真实 runtime tooth acquisition 已通过。 |
| `initMode=auto` | 不再把 static DoA seed 或 truth-frequency init 作为外部初值传入 dynamic case；由 MF estimator 内部 initializer 构造初值。 | No | 更真实暴露 estimator 内部 init / warm-anchor / basin 问题。 | 不能和旧 `staticTruthFreq` 结果直接混成同一性能曲线。 |
| `SS-MF-CP-K` | 单星、多帧、连续相位、known `fdRate`。 | Eval only | 单星 dynamic known-rate local MLE 基线。 | 不能代表 multi-sat 结果。 |
| `SS-MF-CP-U` | 单星、多帧、连续相位、unknown `fdRate` nuisance。 | Eval only | 观察 nuisance-rate release 对 DoA / `fdRef` 的影响。 | 低 SNR `resolvedRate` 下降不能直接解释为 CRB 失效。 |
| `MS-MF-CP-K` | 多星、多帧、连续相位、known `fdRate`。 | Eval only | 用于检查 MS local basin、non-ref coherence 和 fdRef tail。 | 当前不能作为 MS-vs-CRB 主图。 |
| `MS-MF-CP-U` | 多星、多帧、连续相位、unknown `fdRate` nuisance。 | Eval only | 用于暴露 MS unknown-rate warm-anchor / release 失败。 | 当前不能作为 unknown-rate information loss 的性能结论。 |
| `full-sample` | 所有 finite estimator 输出。 | Eval only | 显示 auto-init 工程风险和 tail 污染。 | 不要求 full RMSE 必须贴 CRB。 |
| `resolved` | solver-valid、in-tooth、无 frequency boundary hit；CP-U 还要求 `fdRate` 在预设健康范围内。 | Eval only | 用于排除明显非局部 / boundary / fdRate failure 样本。 | angle error 不参与 loose resolved，不能用它循环筛选 angle 性能。 |
| `coreResolved` | 对 MS 额外要求 `nonRefCoherenceFloor >= 0.8`；对 SS 等同 resolved。 | Eval only | 区分 MS non-ref coherence collapse 与健康 coherence 子样本。 | core 通过不代表 CRB-normalized error 一定健康。 |
| `trimmedCore` | 在 core 样本上固定剔除 `angleNormErr > 5` 或 `fdRefNormErr > 5` 的 tail。 | Eval only | 表示 local-consistent conditional subset，用于观察 tail 剔除后与 CRB 的尺度关系。 | 不是 unconditional efficiency；必须同时报告 `trimmedCoreRate` / `trimKeepRate`。 |
| `trimmedCoreRate` | `trimmedCoreCount / numRepeat`。 | Eval only | 表示所有 repeat 中真正可用于 local-consistent CRB 对比的比例。 | 不能用只剩极少样本的 trimmed 曲线支撑主结论。 |
| `trimKeepRate` | `trimmedCoreCount / coreResolvedCount`。 | Eval only | 表示 core 样本内部有多少没有被 5σ tail 剔除。 | 不等同于总 seed 正确率。 |
| `RMSE/CRB` | `RMSE / CRB standard deviation`。 | Eval only | 用于看误差标准差是否与 CRB 同尺度。 | 小于 1 尤其在 trimmed 条件样本中不代表超过 CRB。 |
| `MSE/CRB` | `mean(error^2) / CRB variance`。 | Eval only | 更直接回答 MSE 是否接近 CRB。 | 仍受条件样本选择影响。 |
| `KKT singular / ill-conditioned warning` | `fmincon` barrier / KKT system 中出现病态或奇异警告。 | No | 指示 MS/unknown 分支优化数值状态不稳定，应进入 replay 诊断。 | 不直接等价于该 seed 失败，需与输出和 failure reason 联合判断。 |

常见口径固定如下：

- `full-sample` 用于工程风险；`resolved/core/trimmed` 用于分层判断是否进入 local CRB 对比区间。
- `trimmedCore` 只能作为 local-consistent conditional curve；它必须和 keep rate / failure reason 一起出现。
- `truth` 只用于 offline scan 评价与 oracle in-tooth local range，不迁移到 estimator 或 flow 的 runtime 分支。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-003909.mat` | 2026-05-07 | `representative` | `SS/MS-MF-CP-K/U`，`initMode=auto`，`P=10`，`SNR=-15:5:10 dB`，`numRepeat=200`，`fdRef half-tooth=0.25`，`fdRate half-width=1000 Hz/s`，`coreCoherenceFloor=0.8`，`trim cap=5 sigma` | SS trimmed/local 子样本接近 CRB；MS trimmed 样本率极低，MS-K same-tooth tail 与 MS-U fdRate-unresolved 明显；scan 固定为统计筛查入口，MS 转 replay 排查。 | 本文档只保留该代表性结果，旧结果不再列出。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:452`
- `numRepeat = 200`
- `snrDbList = [-15, -10, -5, 0, 5, 10]`
- `frameCountList = 10`
- 关键 scan 轴：`SNR × methodNameList`
- active methods：`SS-MF-CP-K`, `SS-MF-CP-U`, `MS-MF-CP-K`, `MS-MF-CP-U`
- `initMode = auto`
- `oracleFdHalfToothFraction = 0.25`
- `oracleFdRateHalfWidthHzPerSec = 1000`
- `resolvedToothHalfWidthFraction = 0.25`
- `resolvedFdRateAbsTolHzPerSec = 250`
- `coreCoherenceFloor = 0.8`
- `trimNormCap = 5`
- `trimMinKeepRate = 0.8`
- task count：`1200`
- checkpoint：enabled；runKey 使用短 hash `f10_snrm15to10_seed253to452_rep200_050e3a5d`；正常完成后已清理 checkpoint artifacts。
- snapshot 保存变量：`scanData`
- 运行时间：约 `1 h 37 m 46 s`

### 4.2 存档数据检查

- 顶层 snapshot 内容：`scanData`
- `scanData` 主要字段：`scanName`, `runKey`, `utcRun`, `config`, `perfTable`, `aggregateTable`, `failureSummaryTable`, `topTailTable`, `repeatOutCell`, `checkpointSummaryTable`, `checkpointCleanupReport`, `plotData`, `elapsedSec`
- 未保存大体量数据：`rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace、图片。
- warning / fail 计数：运行过程中出现多次 `fmincon` KKT singular / ill-conditioned warning，调用链主要经过 `runDoaDopplerMfOptimization`、`runDoaDopplerMfUnknownWarmAnchor`、`solveDoaDopplerMfBranches`。本轮没有导致 scan 中断，但它们与 MS / CP-U 不稳定现象方向一致，后续应进入 replay 诊断。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

#### SS-MF local-consistent 统计

| case | SNR range | samples per SNR | resolved rate | trimmedCoreRate range | trimKeepRate range | trim angle RMSE/CRB | trim fdRef RMSE/CRB | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `SS-MF-CP-K` | `-15:5:10` | 200 | `1.00` all SNR | `0.80 -> 0.93` | `0.80 -> 0.93` | `0.964 -> 0.930` | `1.047 -> 1.032` | Known-rate SS 的 trimmed/local 子样本稳定接近 CRB 尺度；full/core fdRef 被 near-boundary tail 污染。 |
| `SS-MF-CP-U` | `-15:5:10` | 200 | `0.67 -> 0.91` | `0.665 -> 0.91` | `0.993, 0.994, 1.0, 1.0, 1.0, 1.0` | `0.993 -> 0.935` | `0.990 -> 1.089` | Unknown-rate 在低 SNR 有 boundary / fdRate unresolved；resolved 后 trimmed 子样本接近 CRB，但低于 1 的点只能按 conditional subset 解释。 |

#### MS-MF 统计筛查

| case | SNR range | samples per SNR | coreResolvedRate range | trimmedCoreRate range | trimKeepRate range | main failure / pollution | 备注 |
|---|---:|---:|---:|---:|---:|---|---|
| `MS-MF-CP-K` | `-15:5:10` | 200 | `0.695 -> 0.78` | `0.015 -> 0.12` | `0.0216 -> 0.162` | `resolved-coherence-low` 与 `fdRef` near-boundary tail；中高 SNR trimmed 样本率只有约 `1.5%~2.5%`。 | 当前不能作为 MS-vs-CRB 主图；需要 replay 查 same-tooth tail / local basin。 |
| `MS-MF-CP-U` | `-15:5:10` | 200 | `0.265 -> 0.52` | `0 -> 0.09` | `0 -> 0.264` | `fdRate-unresolved` 长期占主导，另有 frequency boundary hit。 | 当前 unknown-rate MS release 不稳；需要 replay 查 warm-anchor / fdRate branch。 |

### 5.2 按扫描轴汇总

#### `SS-MF-CP-K`

| SNR (dB) | resolvedRate | trimmedCoreRate | trimAngleRmseOverCrb | trimAngleMseOverCrb | trimFdRefRmseOverCrb | trimFdRefMseOverCrb |
|---:|---:|---:|---:|---:|---:|---:|
| -15 | 1.000 | 0.800 | 0.9636 | 0.9284 | 1.0473 | 1.0968 |
| -10 | 1.000 | 0.830 | 0.9603 | 0.9221 | 1.0422 | 1.0862 |
| -5 | 1.000 | 0.845 | 0.9543 | 0.9106 | 1.0368 | 1.0749 |
| 0 | 1.000 | 0.855 | 0.9476 | 0.8980 | 1.0468 | 1.0958 |
| 5 | 1.000 | 0.875 | 0.9418 | 0.8869 | 1.0432 | 1.0882 |
| 10 | 1.000 | 0.930 | 0.9298 | 0.8644 | 1.0322 | 1.0655 |

#### `SS-MF-CP-U`

| SNR (dB) | resolvedRate | trimmedCoreRate | trimKeepRate | trimAngleRmseOverCrb | trimAngleMseOverCrb | trimFdRefRmseOverCrb | trimFdRefMseOverCrb | fdRate RMSE (Hz/s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -15 | 0.670 | 0.665 | 0.993 | 0.9933 | 0.9867 | 0.9897 | 0.9795 | 492.37 |
| -10 | 0.785 | 0.780 | 0.994 | 0.9513 | 0.9049 | 1.0918 | 1.1919 | 458.43 |
| -5 | 0.855 | 0.855 | 1.000 | 0.9405 | 0.8845 | 1.1035 | 1.2177 | 380.13 |
| 0 | 0.850 | 0.850 | 1.000 | 0.9396 | 0.8829 | 1.1085 | 1.2288 | 369.64 |
| 5 | 0.905 | 0.905 | 1.000 | 0.9467 | 0.8962 | 1.0900 | 1.1881 | 308.72 |
| 10 | 0.910 | 0.910 | 1.000 | 0.9353 | 0.8748 | 1.0893 | 1.1865 | 292.25 |

#### `MS-MF-CP-K`

| SNR (dB) | resolvedRate | coreResolvedRate | trimmedCoreRate | trimKeepRate | coherenceCollapseRate | trimAngleRmseOverCrb | trimFdRefRmseOverCrb |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -15 | 0.995 | 0.780 | 0.105 | 0.135 | 0.215 | 2.2525 | 0.8850 |
| -10 | 1.000 | 0.740 | 0.120 | 0.162 | 0.260 | 2.5512 | 1.2123 |
| -5 | 1.000 | 0.770 | 0.115 | 0.149 | 0.230 | 3.4897 | 1.7186 |
| 0 | 1.000 | 0.760 | 0.025 | 0.0329 | 0.240 | 2.6336 | 1.4616 |
| 5 | 0.985 | 0.695 | 0.015 | 0.0216 | 0.290 | 2.2918 | 0.3669 |
| 10 | 0.980 | 0.715 | 0.020 | 0.0280 | 0.265 | 2.8741 | 3.3130 |

> 注意：`MS-MF-CP-K` 的 trimmed 样本率太低，上表中的 trimmed RMSE/CRB 不应作为性能曲线，只能作为“极小 local-consistent 子集”的诊断值。

#### `MS-MF-CP-U`

| SNR (dB) | resolvedRate | coreResolvedRate | trimmedCoreRate | trimKeepRate | fdRate-unresolved rate | frequency-boundary-hit rate | trimAngleRmseOverCrb | trimFdRefRmseOverCrb |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -15 | 0.265 | 0.265 | 0.070 | 0.264 | 0.545 | 0.190 | 2.7802 | 1.2246 |
| -10 | 0.350 | 0.350 | 0.090 | 0.257 | 0.500 | 0.150 | 3.2439 | 1.2382 |
| -5 | 0.390 | 0.375 | 0.040 | 0.107 | 0.585 | 0.025 | 3.7469 | 0.8981 |
| 0 | 0.440 | 0.435 | 0.000 | 0.000 | 0.545 | 0.015 | NaN | NaN |
| 5 | 0.500 | 0.495 | 0.010 | 0.020 | 0.495 | 0.005 | 1.3730 | 0.4241 |
| 10 | 0.540 | 0.520 | 0.010 | 0.019 | 0.435 | 0.025 | 2.6600 | 3.5919 |

> 注意：`MS-MF-CP-U` 在 0 dB 没有 trimmedCore 样本；5/10 dB 也只有 2 个 trimmed samples。该方法当前只能作为 failure diagnostic。

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Angle / fdRef RMSE with CRB | SNR | `RMSE` 与 CRB std | `SS/MS-MF-CP-K/U`, full/core/trimmed | SS: Yes after wording check；MS: No | MS trimmed 样本率太低；SS trimmed 曲线必须配 keep rate。 |
| MSE/CRB ratio | SNR | `MSE / CRB variance` | 同上 | SS: conditional curve；MS: diagnostic only | 低于 1 的点只说明 conditional subset，不说明超过 CRB。 |
| Resolved / core / trimmed rate | SNR | rate | loose/core/trimmed rates | Yes, required companion plot | 论文图必须与 RMSE/CRB 同时展示，避免 cherry-pick。 |
| Failure reason summary | SNR / method | failure rate | resolved / fdRate-unresolved / boundary / coherence-low | Appendix / diagnostic | 用于说明 MS 暂缓原因与 CP-U threshold 行为。 |
| Top-tail preview | SNR / seed | CRB-normalized tail score | worst seeds | diagnostic only | 供 replay 选 seed，不直接作为论文曲线。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `SS-MF-CP-K` 在所有 SNR 上 loose resolved rate 为 `1.0`，trimmedCoreRate 从 `0.80` 单调提升到 `0.93`；trimmed angle RMSE/CRB 在 `0.93~0.96`，trimmed fdRef RMSE/CRB 在 `1.03~1.05`。这说明 SS known-rate 的 local-consistent 子样本已经接近 CRB 尺度。
- `SS-MF-CP-U` 在中高 SNR 的 local-consistent 子样本也稳定：`SNR >= -5 dB` 时 trimmedCoreRate 为 `0.855, 0.850, 0.905, 0.910`，trimmed fdRef RMSE/CRB 约 `1.09~1.10`。unknown-rate nuisance 主要体现在低 SNR resolved rate 与 fdRate RMSE，而不是在 trimmed 子样本中完全破坏 DoA / fdRef。
- `MS-MF-CP-K` 的 loose resolved rate 虽高，但 trimmedCoreRate 极低，尤其 `0/5/10 dB` 只有 `0.025/0.015/0.020`。这说明 MS-K 不是单纯 tooth 失败，而是在 same-tooth / non-ref coherence / fdRef tail 处没有稳定进入 local basin。
- `MS-MF-CP-U` 的 `fdRate-unresolved` 率在全 SNR 范围都很高，最低也约 `43.5%`。这说明 MS unknown-rate release 当前不适合继续用大 MC 硬扫，而应进入 seed-level replay。
- top-tail 表中大量 `fdRefAbsErrHz` 聚集在 `173~187 Hz`，接近 `0.25 tooth = 187.5 Hz` 的 local box 边缘，说明许多 tail 是 in-tooth near-boundary branch，而不是全局 wrong-tooth。

### 6.2 反向、污染或未解决现象

- `SS-MF-CP-K` 的 full/core fdRef RMSE/CRB 很大，原因是 core 仍包含 near-boundary fdRef tail。trimmed 后贴 CRB 并不意味着 full estimator 分布已经高效。
- `SS-MF-CP-U` 的 trimmed RMSE/CRB 有低于 1 的点，尤其 -15 dB 的 fdRef 与 angle。该现象来自 resolved/trimmed 条件样本统计，不能写成 estimator 超过 CRB。
- `MS-MF-CP-K` 的 `nonRefCoherenceFloorMedian` 很高，但仍有 large angle/fdRef tail；这说明仅用 coherence floor 不能充分定义 MS local consistency。
- `MS-MF-CP-U` 在运行中伴随多次 KKT singular / ill-conditioned warning，且 failure reason 以 `fdRate-unresolved` 为主；unknown-rate MS branch 需要查 warm-anchor / solver / local box。
- `MS-MF` 的 trimmed 样本数过少，导致一些 trimmed RMSE/CRB 数值看似可用但统计意义不足。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `SS-MF-CP-K`, all SNR | fdRef tail | top-tail seed 反复出现 `fdRefAbsErrHz≈173~175 Hz`，但 failure reason 仍为 `resolved`。 | 说明 loose resolved 不足以定义 local CRB consistency；trimmed/local 口径必须保留。 |
| `SS-MF-CP-U`, low SNR | fdRate / boundary | -15 dB resolved rate `0.67`，frequency-boundary-hit `0.215`，fdRate-unresolved `0.115`。 | 说明 unknown-rate nuisance 在低 SNR 有 threshold 行为。 |
| `MS-MF-CP-K`, 0/5/10 dB | same-tooth tail | trimmedCoreRate 只有 `0.025/0.015/0.020`，coreResolvedRate 仍约 `0.695~0.76`。 | 说明 MS-K 的 core 条件仍混入大量 non-local tail，不能作为论文主图。 |
| `MS-MF-CP-U`, 0 dB | no trimmed sample | trimmedCoreRate `0`，trimmed RMSE/CRB 为 NaN。 | 直接阻止 MS-U 作为当前 MLE-vs-CRB 主结果。 |
| MS warning traces | solver warning | 多次 `fmincon` KKT singular / ill-conditioned warning 出现在 `runDoaDopplerMfOptimization` / `UnknownWarmAnchor` 链。 | 支持后续用 replay 查 solver / warm-anchor，而不是继续扩 scan。 |

## 7. 机制解释

### 7.1 当前解释

本次结果把 `scanMfMleCrbInToothConsistency` 的职责进一步固定为“统计筛查”和“是否值得进入论文图”的判定入口。它已经能清楚区分两类现象：

1. `SS-MF` 在 auto init 条件下，good / local-consistent 子样本已经可以达到 CRB 尺度附近；这说明 reference-sat 单星 continuous-phase MLE 主链基本健康。
2. `MS-MF` 在同样 auto init 与 in-tooth local box 条件下，仍无法稳定进入 multi-sat local basin。MS-K 的问题表现为 same-tooth tail、near-boundary `fdRef` branch 和 coherence / local consistency 判据不足；MS-U 的问题表现为 `fdRate-unresolved` 和 boundary / warm-anchor 不稳。

这里的关键不是“MS CRB 不存在”或“MS 理论不成立”。MS 有对应 CRB，但当前 estimator auto-init / solver / local basin 还不能稳定产生与该 CRB 对应的样本分布。因此，MS 暂不应继续靠扩大 repeat 形成论文主结果，而应该先用 replay 对 seed-level 路径进行分类。

### 7.2 这个 scan 支持什么

- 支持将 `SS-MF-CP-K/U` 作为后续 local-consistent MLE-vs-CRB 的先行统计基线。
- 支持论文仿真采用 `full/core/trimmed + keep rate + failure reason` 的联合报告口径，而不是只画一条 RMSE/CRB 曲线。
- 支持把 `trimmedCore` 明确命名为 `local-consistent conditional subset`，并把 full/core tail 作为工程风险或 threshold/outlier 风险报告。
- 支持将 MS 失败转入 replay，而不是继续在 scan 中扩大 repeat。
- 支持下一步优先排查 `MS-MF-CP-K` 的 same-tooth near-boundary fdRef branch 与 `MS-MF-CP-U` 的 fdRate release / warm-anchor branch。

### 7.3 这个 scan 不证明什么

- 不证明 `MS-MF` 已接近 MS CRB。
- 不证明 `trimmed RMSE/CRB < 1` 是 estimator 超过 CRB。
- 不证明 in-tooth oracle local box 等价于真实 coarse Doppler compensation 的全部实现。
- 不证明 full-flow tooth acquisition、subset selection 或 rescue flow 已经通过。
- 不证明当前结果可以写 regression 契约。
- 不证明 MS 理论模型错误；当前更像 estimator auto-init / local basin / solver 路径问题。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。当前 scan 只暴露统计现象；不直接修改 estimator 主核、objective、residual、reference-sat 语义或 search bound。 |
| flow 默认路径 | 不改。该 scan 是 in-tooth local consistency scan，不回流到 subset / rescue flow。 |
| replay 下一步 | 需要。建议新增或运行 `replayMfSsMleCrbMetricDiagnose` 与 `replayMfMsMleCrbMetricDiagnose`。SS replay 查 angle CRB metric、fdRef tail、CP-U release；MS replay 查 non-ref coherence / near-boundary fdRef branch / fdRate-unresolved / KKT warning seed。 |
| regression | 不写。只有 replay 证明稳定契约后，再考虑补最小 regression，例如 CP-U release、reference-sat 不变量、boundary / warm-anchor 不回退等。 |
| 论文图 | SS 可继续作为 main figure candidate；MS 暂时只能放 diagnostic / appendix，不能作为正文 MS-vs-CRB 主图。 |
| 排障记录 | 主记录可摘一句：`scanMfMleCrbInToothConsistency_20260507-003909` 固定为统计筛查入口，SS 可用，MS 转 replay 排查。机制归并版可补充 MS tail / fdRate-unresolved 证据。 |

### 后续仿真应围绕什么展开

后续 paper-facing 仿真应先围绕 **SS-MF local-consistent MLE-vs-CRB** 展开，配套报告：

- `full/core/trimmed` RMSE 与 MSE；
- CRB-normalized error；
- `resolvedRate`、`trimmedCoreRate`、`trimKeepRate`；
- failure reason summary；
- known / unknown-rate 信息损失；
- 后续再扩展到 `P` / SNR / CP-IP / known-unknown 的论文主图。

MS 不应继续作为大 MC 主线硬扫。MS 当前应围绕以下 replay 目标展开：

- `MS-MF-CP-K`：same-tooth fdRef near-boundary branch、core 通过但 trimmed 失败的 seed、non-ref coherence 高但 angle/fdRef 仍错的样本。
- `MS-MF-CP-U`：`fdRate-unresolved`、frequency boundary hit、warm-anchor selected seed、KKT singular warning 与 final `fdRate/fdRef` 的关系。

### 达到什么条件后可以回到本 scan

满足以下条件后，再回到 `scanMfMleCrbInToothConsistency` 做正式 rerun：

1. **SS 条件**：确认 angle error 与 angle CRB metric 同口径，或在结果文档中明确写清 conditional / transformed angular CRB 口径；`SS-MF-CP-K/U` 在目标 SNR 区间的 `trimmedCoreRate >= 0.85`，且 trimmed MSE/CRB 大体在 `0.8~1.3` 或有明确有限样本解释。
2. **MS-K 条件**：通过 replay 解释并修复或稳定 gate 掉 near-boundary fdRef tail；中高 SNR 下 `coreResolvedRate >= 0.8`，`trimmedCoreRate >= 0.7~0.8`，`coherenceCollapseRate <= 0.1~0.15`，且 trimmed 样本不是个位数。
3. **MS-U 条件**：通过 replay 降低 `fdRate-unresolved` 与 frequency-boundary-hit；中高 SNR 下 `resolvedRate >= 0.8`，`trimmedCoreRate >= 0.7`，`fdRate-unresolved <= 0.1~0.15`。
4. **solver 条件**：KKT singular / ill-conditioned warning 不再集中出现在同一类 successful-looking result 中；若 warning 仍存在，必须能在 failure summary 或 replay 分类中解释，不污染 paper-facing curves。
5. **文档条件**：scan 只保留统计和 compact top-tail，不把 replay 内部 trace 下沉到 scan；如果需要新的内部诊断，继续放 replay。

## 9. 限制与禁止解释

- 不要把 offline truth evaluation 字段迁移到 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 in-tooth oracle box 解释成 full-flow tooth acquisition 已解决。
- 不要把 trimmed-only 曲线单独作为 estimator efficiency 结论；必须同时报告 `trimmedCoreRate`、`trimKeepRate` 和 failure reason。
- 不要把 `trimmed RMSE/CRB < 1` 写成超过 CRB；这是条件样本统计效应或 metric 口径 warning。
- 不要把 MS 当前极低 trimmedCoreRate 下的 RMSE/CRB 值作为论文主图或正结论。
- 不要用单一 hit rate 代替 RMSE / P95 / P99 / CRB-normalized error / resolved rate / outlier rate。
- 不要把当前 scan 结果直接迁移为 regression；replay 分类和稳定契约形成前，regression 太早。
- 不要继续在本 scan 中堆 objective trace、block compare、warm-anchor trace 或 per-frame residual matrix；这些属于 replay。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-003909.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfMleCrbInToothConsistency.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact aggregate table、failure summary、top-tail preview 和图。

## 11. 历史备注

- 本文档只保留 `scanMfMleCrbInToothConsistency_20260507-003909.mat` 这一组代表性结果。
- 旧 SS-only、旧 MS/SS mixed 或 `staticTruthFreq` 结果不在本文档中保留，避免读者混淆 auto init 与 oracle/static init 口径。
- 当前版本的后续角色已经固定：scan 负责统计筛查和论文候选曲线；replay 负责 seed-level estimator / metric / solver 排障；regression 等 replay 形成稳定契约后再考虑。
