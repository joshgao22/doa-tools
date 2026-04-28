# scanMfCpIpTying 结果记录

## 对应 scan

- `test/dev/scan/scanMfCpIpTying.m`

## 扫描目标

比较 `continuous / relaxed / independent` 三种 phase tying 在同一多帧样本上的 `fdRef` 一维 objective 结构，重点观察 `1/T_f` comb tooth、truth tooth 与 wrong tooth 的 objective gap。

该 scan 是 CP/IP phase-tying mechanism scan，不是 Monte Carlo 性能 scan，也不是 flow replay。它固定场景、固定中心点后沿 `fdRef` 网格扫描 objective，用于解释连续相位 tying 如何影响 comb 结构，以及 CP-K / CP-U 的最终中心是否落在正确 tooth。

## Snapshot index

| snapshot | 状态 | 配置 | 结论 |
|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpTying_20260428-143341.mat` | current | `numFrame=20`，`T_f=1/750 s`，`SNR=10 dB`，`phaseList=[continuous, relaxed, independent]`，`numGrid=801`，`numAliasSide=2`，CP-K / CP-U 均运行，只保存轻量 `scanData` | 当前代表性结果。CP continuous 在 truth center 上仍把 `deltaFdRef=0` 作为最低点，但 CP-K final 与 CP-U final 都落到 `-1 tooth` 附近；IP/relaxed 不是更好的 `fdRef` 估计，只是弱化跨帧相位约束后的 baseline。 |
| 未绑定 snapshot | superseded | 2026-04-28 控制台日志；旧版 `saveSnapshot=false`，`numFrame=10`，`T_f=1/750 s`，`SNR=10 dB`，`phaseList=[continuous, relaxed, independent]`，`numGrid=801`，CP-K / CP-U 均运行 | 旧日志显示 CP-U final 曾在该 10 帧样本回到 center tooth。该结论只作为历史观察保留，不再作为当前代表性结果。 |

## 存档数据检查

`scanMfCpIpTying_20260428-143341.mat` 为 `saveExpSnapshot` 生成的轻量 snapshot：

- 顶层变量为 `data / meta / inventory`。
- `data` 中只保存 `scanData`，未保存 `rxSigCell`、完整 `sceneSeq` 或其它大体量运行对象。
- `scanData` 包含 `truth`、`scanConfig`、`centerRegistry`、`frameInfo`、`cpKnown`、`cpUnknown`、`scanName`、`utcRun`、`snapshotFile`。
- `frameInfo` 显示本次实际运行为 `numFrame=20`，`frameIntvlSec=1/750 s`，`refFrameIdx=10`，`numPilotSample=512`，`refSatIdxLocal=1`。
- `cpKnown` 与 `cpUnknown` 均包含 `3 phase modes x 3 centers x 801 fdRef points` 的扫描数据，并保存 `summaryTable / toothTable` 以及可重画曲线的 `phaseScan`。

该 snapshot 与当前 scan 落盘规范一致：结果分析写在本 md，`.mat` 文件只作为可恢复的 scan 数据存放在 `test/data/cache/scan/`。

## 当前代表性结果

本次代表性 snapshot 配置为 20 帧、`T_f=1/750 s`、`SNR=10 dB`。KNOWN 与 UNKNOWN 各自扫描 `3 phase modes x 3 centers x 801 fdRef points`。

### CP-K phase / center summary

| phase mode | center | minDeltaFd (Hz) | aliasIdx | centerDeltaObj |
|---|---|---:|---:|---:|
| continuous | truth | 0.000000 | 0 | 0 |
| continuous | staticSeed | -348.750000 | 0 | 322793.798204 |
| continuous | finalEstimate | -750.000000 | -1 | 1.494901 |
| relaxed | truth | -251.250000 | 0 | 0.682014 |
| relaxed | staticSeed | -457.500000 | -1 | 2.213423 |
| relaxed | finalEstimate | -277.500000 | 0 | 0.817795 |
| independent | truth | -251.250000 | 0 | 0.682014 |
| independent | staticSeed | -457.500000 | -1 | 2.213423 |
| independent | finalEstimate | -277.500000 | 0 | 0.817795 |

### CP-U phase / center summary

| phase mode | center | minDeltaFd (Hz) | aliasIdx | centerDeltaObj |
|---|---|---:|---:|---:|
| continuous | truth | 0.000000 | 0 | 0 |
| continuous | cpKnownSeed | -750.000000 | -1 | 1.494901 |
| continuous | finalEstimate | -750.000000 | -1 | 4.077799 |
| relaxed | truth | -251.250000 | 0 | 0.682014 |
| relaxed | cpKnownSeed | -277.500000 | 0 | 0.817795 |
| relaxed | finalEstimate | -630.000000 | -1 | 4.237186 |
| independent | truth | -251.250000 | 0 | 0.682014 |
| independent | cpKnownSeed | -277.500000 | 0 | 0.817795 |
| independent | finalEstimate | -630.000000 | -1 | 4.237186 |

## 可观察现象

### 1. Continuous phase tying 不能消除 comb

在 truth center 下，`continuous` 的最低点仍是 `deltaFdRef=0`，说明 CP 模型在正确中心附近保留了区分正确 tooth 的信息。但相邻 alias tooth 仍然存在，且仍可能吸引实际 flow。

CP-K / CP-U 的 `continuous | truth` alias table 为：

| aliasIndex | deltaFdRef (Hz) | deltaObj |
|---:|---:|---:|
| -2 | -1500 | 15.890611 |
| -1 | -750 | 1.956850 |
| 0 | 0 | 0 |
| 1 | 750 | 10.020294 |
| 2 | 1500 | 32.016549 |

相较旧 10 帧日志，本次 20 帧结果中 truth tooth 与近邻 tooth 的 gap 变大，但 `-750 Hz` wrong tooth 仍然是明确可见的局部支路。因此 CP continuous 不能写成“消除 comb ambiguity”，只能写成“在正确中心附近增强 tooth 区分”。

### 2. CP-K final 落到 `-1 tooth` 附近

`CP-K | continuous | finalEstimate` 的最低点是 `-750 Hz`，`aliasIdx=-1`，且 `centerDeltaObj=1.494901`。该中心附近的 alias table 为：

| aliasIndex | deltaFdRef (Hz) | deltaObj |
|---:|---:|---:|
| -2 | -1500 | 4.784291 |
| -1 | -750 | 0 |
| 0 | 0 | 1.494901 |
| 1 | 750 | 9.268793 |
| 2 | 1500 | 23.320762 |

这说明 CP-K final 已经停在 wrong-tooth 分支附近。即使 CP continuous 在 truth center 上是健康的，flow / seed 一旦到 wrong tooth，局部 objective 也会支持该错齿。

### 3. CP-U final 本次没有把 final center 拉回 center tooth

旧 10 帧日志中，`CP-U | continuous | finalEstimate` 曾回到 `deltaFdRef=0`。但当前 20 帧 snapshot 显示，`CP-U | continuous | finalEstimate` 的最低点仍是 `-750 Hz`，`aliasIdx=-1`，`centerDeltaObj=4.077799`。

该中心附近的 alias table 为：

| aliasIndex | deltaFdRef (Hz) | deltaObj |
|---:|---:|---:|
| -2 | -1500 | 7.891140 |
| -1 | -750 | 0 |
| 0 | 0 | 4.077799 |
| 1 | 750 | 20.124055 |
| 2 | 1500 | 48.136865 |

因此，当前代表性结果更保守也更清楚：unknown-rate release 并不保证自动救回 comb wrong tooth。它可以改变 objective 与中心位置，但选齿仍需要 subset / non-periodic schedule / flow rescue 来承担。

### 4. relaxed 与 independent 完全一致是合理现象

`relaxed` 与 `independent` 在本次输出中仍然数值一致。原因是二者都允许更局部的 frame/block complex gain profile，absolute-time 与 local-time 之间的常数相位差会被独立相位吸收。因此它们在该 scan 中都表示“切断或弱化跨帧 continuous phase tying 后”的 IP-like baseline。

### 5. independent phase 不是更好的 `fdRef` 估计

本次结果不支持“独立相位的 `fdRef` 更好 / 没有 comb”这一解读：

- 在 truth center 下，`independent` 的最优点是 `-251.25 Hz`，不是 `0 Hz`；而 `continuous` 的最优点是 `0 Hz`。
- 在 CP-K final center 下，`independent` 的最低点是 `-277.5 Hz`，只是弱化/平滑了 `fdRef` 约束，并没有形成物理上更正确的 `fdRef`。
- 在 CP-U final center 下，`independent` 的最低点是 `-630 Hz`，仍接近 `-1 tooth`，也没有救回 center tooth。

independent 看起来 comb 结构不那么尖锐，更像是跨帧绝对相位约束被释放后 `fdRef` 可辨识性被削弱，而不是更好的 Doppler 估计。

## 当前结论

当前代表性结果应写成：**连续相位 tying 不能消除 `1/T_f` comb；它只是在正确中心附近把正确 tooth 作为最低点，并给出有限 alias gap。若 flow / seed 已经落到 wrong tooth，CP objective 仍可能稳定支持该错齿；CP-U release 也不能保证自动把 final center 拉回 truth tooth。**

因此，论文中不应把 CP 表述为“解决 comb ambiguity”的机制，而应表述为：

- CP 是主物理模型，因为它保留跨帧连续相位信息；
- IP / relaxed 是切断这部分信息后的对比基线；
- comb / wrong-tooth 是 CP objective 仍然保留的结构性风险；
- 实际算法仍需要 subset / non-periodic selection、conditional rescue 和 same-tooth refine 来稳定接住正确 tooth 与正确 DoA basin。

## 对代码和论文图的影响

- 该 scan 保留在 `test/dev/scan/` 是合适的，因为它扫描的是 objective curve / folded tooth surface，而不是固定 seed 的 flow replay。
- 图形应保留每个 subplot 的 phase-mode legend，但 center line 与 alias reference line 不应进入 legend；当前代码已把这些辅助线的 `HandleVisibility` 关闭。
- 新版代码默认保存轻量 `scanData`，其中包含 `summaryTable`、`toothTable` 和可重画曲线数据；不保存图片，也不保存大体量原始观测。
- 该结果可作为论文 CP/IP mechanism 图候选，但不能替代 `scanMfCpIpPerfMap` 的 SNR / frame-count 性能图。

## 后续建议

`scanMfCpIpTying` 当前不需要继续改主逻辑。下一步若围绕论文机制图推进，可以直接基于该 snapshot 重出表格和图；若围绕 estimator / flow 推进，应回到 subset tooth selection、same-tooth tail diagnose 和 conditional refine，而不是把 IP 改成主模型。

## 历史 / superseded snapshots

- 未绑定 snapshot：2026-04-28 控制台日志；`numFrame=10`，`saveSnapshot=false`。该旧结果中 CP-U final 曾回到 center tooth，但已被当前 20 帧 snapshot 取代为代表性结果。
