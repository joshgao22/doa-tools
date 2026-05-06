# scanMfBlockLength 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / appendix candidate` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfBlockLength_20260428-105157.mat` |
| 当前一句话结论 | 增大 pilot block length 会显著增强 truth-centered `fdRef` comb / alias tooth separation；但即使 4x block，块内 Doppler-rate 二次相位仍极小。 |
| 论文图定位 | `mechanism figure / appendix candidate`，可用于解释 block length 与 single-block Doppler / tooth 分辨力。 |
| 决策影响 | 固定当前机制结论；不进入 regression；不把结果迁移成 estimator 默认路径修改。 |
| 下一步动作 | 暂停扩跑；若论文需要图，可基于该 snapshot 重画 `aliasGap` vs `blockLenBaseRatio` 的 log-y 曲线。 |
| 禁止误用 | 不能把该机制 scan 写成“单个同步块内部 Doppler-rate 已显著”，也不能用 final-centered 诊断替代 truth-centered alias-gap 主结论。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfBlockLength.m`
- 结果文档：`test/dev/scan/results/scanMfBlockLength.md`
- scan 类型：`mechanism scan / block-length alias-tooth scan`
- 主要问题：pilot block length 对连续相位多帧 `fdRef` 一维 objective comb / alias tooth separation 有多大影响。
- 扫描对象：`blockLen=[528,1056,2112,4224,8448]`，`CP-K / CP-U`，truth / seed / final centers。
- 不覆盖范围：不比较 subset bank、rescue flow、same-tooth polish；不做 Monte Carlo estimator 性能统计；不验证 regression 契约。
- truth 使用口径：truth center 用于机制线扫和 alias-gap 评价；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Appendix candidate / mechanism figure only。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `blockLen` | 单帧参与估计的 pilot block 样本数。 | No | 横轴；越大表示单块相干观测更长。 | 不代表帧数 `P`。 |
| `blockLenBaseRatio` | 相对默认 `2112` samples 的比例。 | No | 便于比较 `0.25x / 0.5x / 1x / 2x / 4x`。 | 不要和 `blockLenRatio` 相对 full block 混淆。 |
| `aliasGap1 / aliasGap2` | truth center 下相邻 / 次邻 alias tooth 与 center tooth 的 objective gap。 | Oracle center | 衡量 wrong tooth 被压低的强度；主机制指标。 | 不代表真实 flow tooth hit-rate。 |
| `inBlockFdDriftHz` | 单个 block 内由 `fdRate` 引起的 Doppler drift。 | No | 说明块内 Doppler 动态量级。 | 不能解释跨帧 Doppler dynamic 是否重要。 |
| `inBlockQuadPhaseRad` | 单个 block 内二次相位量级。 | No | 用于证明 block 内 rate 仍可忽略。 | 不要写成单块 rate 已主导。 |
| `finalEstimate` center | 围绕 estimator final point 的线扫。 | Eval only | 诊断 final basin / wrong-tooth 吸引。 | 不应和 truth-centered alias-gap 机制主图混放。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：N/A；本 scan 不是 MC performance scan。
- `resolved-sample`：N/A；不做 CRB / MLE resolved 对比。
- `outlier rate`：N/A。
- `truth-tooth / oracle range`：truth center 只用于机制评价。
- `stress-test`：否；本 scan 是机制曲线。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfBlockLength_20260428-105157.mat` | 2026-04-28 | `representative` | `Fs=512 MHz`；`symbolRate=128 MHz`；`osf=4`；`baseBlockLen=2112`；`standardBlockLen=8448`；`numFrame=10`；`T_f=1/750 s`；`SNR=10 dB`；`numGrid=801`；`numAliasSide=2`；CP-K / CP-U 均运行；checkpoint resume enabled。 | `aliasGap1` 从 528 samples 的 `0.2998` 增至 8448 samples 的 `1.357e4`；4x block 内二次相位仍只有 `3.28e-6 rad`。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = N/A`；该 scan 使用固定机制样本和固定 shared phase。
- `numRepeat = N/A`；非 Monte Carlo。
- `snrDb = 10`
- `frameCount = 10`
- `frameIntvlSec = 1/750`
- `blockLenList = [528, 1056, 2112, 4224, 8448]`
- `numGrid = 801`
- `numAliasSide = 2`
- `scanHalfWidthHz = []`，由 alias side / grid 构造线扫范围。
- 关键 resolved / outlier 判据：N/A；本 scan 是 objective line mechanism scan。
- checkpoint：enabled；本次从部分 block task resume，最终 `numTask=5`、`numDone=5`、`isComplete=true`，成功后清理 tmp 目录。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录 `elapsedSec`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`truth`、`scanConfig`、`frameInfo`、`blockLenList`、`blockResult`、`summaryKnown`、`summaryUnknown`、`toothKnown`、`toothUnknown`、`aggregateTable`、`checkpointSummaryTable`、`checkpointCleanupReport`、`scanName`、`runKey`、`utcRun`、`snapshotFile`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：无 hard fail；checkpoint 完整结束并清理 tmp。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

主结果只取 `truth` center 的 `fdRef` line。当前 snapshot 中 CP-K 与 CP-U 在 truth-centered alias gap 上一致，因此表中合并展示。

| blockLen | 相对 2112 | block duration (us) | `|fdRate| * T_block` (Hz) | quad phase (rad) | aliasGap1 | aliasGap2 | 备注 |
|---:|---:|---:|---:|---:|---:|---:|---|
| 528 | 0.25x | 1.03125 | 0.003953 | `1.28e-8` | 0.2998 | 1.2366 | comb 很平，wrong tooth 区分弱。 |
| 1056 | 0.5x | 2.06250 | 0.007907 | `5.12e-8` | 2.7226 | 10.9158 | tooth gap 开始增强。 |
| 2112 | 1x | 4.12500 | 0.015813 | `2.05e-7` | 24.3386 | 98.7082 | 默认 block 已有明显 tooth separation。 |
| 4224 | 2x | 8.25000 | 0.031626 | `8.20e-7` | 275.3637 | 1118.0680 | wrong tooth 被明显压低。 |
| 8448 | 4x | 16.50000 | 0.063253 | `3.28e-6` | 13574.8601 | 53875.1723 | alias gap 极大，但块内二次相位仍极小。 |

### 5.2 按扫描轴汇总

| axis value | case | metric 1 | metric 2 | metric 3 | 解释 |
|---:|---|---:|---:|---:|---|
| `528 -> 2112` | truth-centered CP-K/U | `aliasGap1: 0.2998 -> 24.3386` | `aliasGap2: 1.2366 -> 98.7082` | quad phase `1.28e-8 -> 2.05e-7` | 从 0.25x 到默认 block，tooth separation 明显增强。 |
| `2112 -> 8448` | truth-centered CP-K/U | `aliasGap1: 24.3386 -> 13574.8601` | `aliasGap2: 98.7082 -> 53875.1723` | quad phase `2.05e-7 -> 3.28e-6` | 4x block 极大增强 wrong-tooth gap，但不是因为块内 rate 变显著。 |
| all block | CP-K vs CP-U | truth-centered alias gap identical | N/A | N/A | 当前机制主要由 block length 主导，不由 known/unknown branch 主导。 |

final-centered 诊断结果：

| mode | blockLen | aliasGap1 | aliasGap2 | centerDeltaFinal | minDeltaFinalHz | 解释 |
|---|---:|---:|---:|---:|---:|---|
| CP-K | 528 | 0.2998 | 1.2366 | 3.1360 | 1500 | 短 block final center 可落到远齿。 |
| CP-K | 1056 | 2.7226 | 10.9158 | 167235.0001 | -375 | final-centered 诊断异常，不作为主机制指标。 |
| CP-K | 2112 | 24.3386 | 98.7082 | 0 | 0 | 默认 block final center 正常。 |
| CP-K | 4224 | 275.3637 | 1118.0680 | 0 | 0 | 长 block final center 正常。 |
| CP-K | 8448 | 13574.8601 | 53875.1723 | 0 | 0 | 长 block final center 正常。 |
| CP-U | 528 | 0.2998 | 1.2366 | 3.1360 | 1500 | 短 block 不稳。 |
| CP-U | 1056 | 2.7226 | 10.9158 | 0 | 0 | final center 正常。 |
| CP-U | 2112 | 24.3386 | 98.7082 | 0.6718 | 750 | 可见邻齿吸引。 |
| CP-U | 4224 | 275.3637 | 1118.0680 | 8.5554 | 750 | final-centered 仍有邻齿诊断。 |
| CP-U | 8448 | 13574.8601 | 53875.1723 | 0 | 0 | 长 block 正常。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `Truth-centered alias gap vs block length` | `blockLenBaseRatio` 或 `blockLen` | `aliasGap1 / aliasGap2` | CP-K / CP-U 可合并 | Appendix / mechanism figure | 建议 log-y；只用 truth center。 |
| `In-block dynamic scale` | `blockLenBaseRatio` | `inBlockFdDriftHz`、`inBlockQuadPhaseRad` | single curve | Appendix / text support | 用于说明块内 rate 动态仍小。 |
| `Final-centered minDelta` | `blockLen` | `minDeltaFinalHz` | CP-K / CP-U | Diagnostic only | 不和 truth-centered alias-gap 主图混放。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `aliasGap1` 随 block length 快速增大：`528 -> 2112 -> 8448` 对应 `0.2998 -> 24.3386 -> 13574.8601`。
- `aliasGap2` 同样快速增大：`1.2366 -> 98.7082 -> 53875.1723`。
- 4x block 下 `inBlockFdDriftHz` 也只有 `0.063253 Hz`，`inBlockQuadPhaseRad` 只有 `3.28e-6 rad`，说明单块内部 Doppler-rate 动态仍然不是主导因素。
- CP-K 与 CP-U truth-centered alias gap 一致，说明 block length 是该机制 scan 的主要变量。
- checkpoint resume / cleanup 成功，说明该重 scan 的工程外壳可用。

### 6.2 反向、污染或未解决现象

- final-centered 诊断在短 block 下不稳定，例如 528 samples 的 `minDeltaFinalHz=1500`，1056 samples 的 CP-K `centerDeltaFinal` 很大。该现象说明 final basin 可受短 block 影响，但不是 truth-centered alias-gap 主结论。
- 该 scan 没有 Monte Carlo repeat，因此不能说明更长 block 一定提高 full-flow hit-rate。
- 该 scan 不比较不同 block length 的计算成本 / 实际接收可用性。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `blockLen=528, finalEstimate` | diagnostic | `minDeltaFinalHz=1500` | 暴露短 block final-centered wrong-tooth 风险。 |
| `blockLen=1056, CP-K finalEstimate` | diagnostic | `centerDeltaFinal≈1.67235e5`，`minDeltaFinalHz=-375` | 说明 final-centered 线扫受 basin 影响，不能替代 truth-centered机制指标。 |
| `blockLen=8448, truth` | mechanism upper end | `aliasGap1≈1.357e4`，quad phase `3.28e-6 rad` | 说明长 block tooth separation 强，但块内 rate 仍很小。 |

## 7. 机制解释

### 7.1 当前解释

更长的 pilot block 提供更长的单块相干积分时间，因此 `fdRef` 线扫中不同 alias tooth 的相位匹配差异被放大，truth-centered alias gap 快速增加。这解释了为什么默认 `2112` samples 已明显优于 528 / 1056，而 2x / 4x block 可进一步压低 wrong tooth。

但 block 内时间尺度仍是微秒级。当前 truth `fdRateFit≈-3833.5 Hz/s`，即使 `8448` samples 对应 `16.5 us`，块内 Doppler drift 也只有 `0.063 Hz`，二次相位约 `3.28e-6 rad`。所以本文 Doppler dynamic 的必要性仍来自跨帧观测窗口，而不是单个同步块内部的 rate 曲率。

### 7.2 这个 scan 支持什么

- 支持更长 known pilot block 可以增强单块 `fdRef` / alias tooth 分辨力。
- 支持默认 `2112` samples 已比 528 / 1056 有明显 stronger tooth gap。
- 支持 block-length 图可作为 mechanism / appendix 证据。
- 支持把跨帧 Doppler dynamic 与块内 Doppler-rate 动态区分开。

### 7.3 这个 scan 不证明什么

- 不证明 full-flow tooth hit-rate 已随 block length 提升。
- 不证明 estimator 默认路径需要改 block length。
- 不证明单块内部必须引入 Doppler-rate 模型。
- 不证明可以写 regression 契约。
- 不证明 4x block 在真实系统中一定可用或代价可接受。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 不要求修改 estimator。 |
| flow 默认路径 | 不改；没有验证 subset / rescue flow。 |
| replay 下一步 | 不需要；若未来验证 flow hit-rate，应另建 flow / MC scan。 |
| regression | 不写；机制趋势不是自动 pass/fail 契约。 |
| 论文图 | `appendix / mechanism figure`，主图建议画 truth-centered `aliasGap1/2` vs block length。 |
| 排障记录 | 可摘一句“longer block improves tooth separation, but in-block rate remains negligible”；不复制长表。 |

## 9. 限制与禁止解释

- 不要把 truth-centered alias gap 当作 runtime selected hit-rate。
- 不要把 final-centered `minDeltaFinalHz` 当作主机制指标。
- 不要写成单个 pilot block 内 Doppler-rate 已显著。
- 不要把 block-length mechanism scan 迁移成 regression。
- 不要在该 scan 内继续加入 subset bank、rescue、same-tooth polish 等正交问题。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfBlockLength_20260428-105157.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfBlockLength.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：当前脚本如已 template 化，checkpoint / notification 外壳可能和历史 snapshot 不完全一致；机制数值以本节 snapshot 为准。

## 11. 历史备注

- 当前只绑定 `scanMfBlockLength_20260428-105157.mat` 作为代表性 snapshot。
- 本 snapshot 已包含 checkpoint resume / cleanup 成功记录，后续没有必要为了工程链路重复跑大规模 block scan。
