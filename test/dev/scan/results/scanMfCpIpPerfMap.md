# scanMfCpIpPerfMap 结果记录

## 对应 scan

- `test/dev/scan/scanMfCpIpPerfMap.m`

## 扫描目标

扫描 `CP-K / CP-U / IP-K / IP-U` 在不同 SNR 与联合帧数下的实际估计性能。该 scan 使用当前 simple dynamic flow 生成 `CP-K / CP-U` 结果，再从同一 static seed 上重跑 `IP-K / IP-U`，用于观察当前完整 flow 下 CP/IP 与 known/unknown-rate 的性能表现。

需要特别说明：该 scan 是 **full-flow performance map**，不是受控的 in-tooth oracle。它会同时受到 wrong-tooth、same-tooth bad basin、candidate adoption 和 unknown warm-anchor release 的影响。因此，当前结果只能用于判断“现有 flow 是否已经足以支撑 CP/IP 性能图”，不能直接解释为 CP/IP 统计模型本身的理论优劣。

## Snapshot index

| snapshot | 状态 | 配置 | 结论 |
|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpPerfMap_20260428-195242.mat` | current | `baseSeed=253`，`seedList=[253,254,255]`，`SNR=[0,5,10,15] dB`，`P=[8,10,20]`，`T_f=1/750 s`，每格 `numRepeat=3`，truth-tooth 判据为 `toothIdx==0 && abs(toothResidualHz)<=50 Hz`，只保存轻量 `scanData` | 当前代表性结果。完整 36 个 task 均跑通，但 CP-K / CP-U 的 truth-tooth hit-rate 明显低于 IP-K / IP-U；该结果不能作为“CP 优于 IP”的论文性能图，只能说明当前 full flow 仍被 wrong-tooth / same-tooth basin 污染。 |

## 存档数据检查

`scanMfCpIpPerfMap_20260428-195242.mat` 为 `saveExpSnapshot` 生成的轻量 scan snapshot：

- 顶层变量为 `data / meta / inventory`。
- `data` 中只保存 `scanData`。
- `scanData` 包含 `config`、`perfTable`、`aggregateTable`、`repeatOutCell` 与 `plotData`。
- `config` 中记录 `baseSeed=253`、`numRepeat=3`、`seedList=[253,254,255]`、`snrDbList=[0,5,10,15]`、`frameCountList=[8,10,20]`、`toothResidualTolHz=50`。
- snapshot 中未保存大体量 `rxSigCell`、完整 `sceneSeq`、完整 fixture cache 或图片文件，符合 scan cache 的轻量保存口径。

运行日志显示 36/36 个 task 完成，总耗时约 `22 min 11 s`。日志中出现 `84` 次 `矩阵接近奇异值` warning 和 `2` 次 `矩阵在工作精度内为奇异的` warning，主要来自 unknown warm-anchor / fmincon 内部 release seed 链路。这些 warning 没有导致 hard fail，但说明当前 full-flow scan 的 unknown-stage 局部优化条件数仍然较差，后续应在 scan 层做 compact warning 计数，而不是把 warning 压到 estimator 主核里。

## 当前代表性结果

### 全部格点汇总

下表把 `P x SNR x repeat` 的 36 个样本按 case 汇总。`truthToothHitRate` 使用当前 scan 的 strict 判据：`toothIdx==0 && abs(toothResidualHz)<=50 Hz`。

| case | samples | truth-tooth hit | wrong-tooth / residual fail | angle RMSE (deg) | angle median (deg) | angle p95 (deg) | fdRef RMSE (Hz) | fdRef median (Hz) | fdRef p95 (Hz) | fdRate RMSE (Hz/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `CP-K` | 36 | 0.333 | 0.667 | 0.004460 | 0.003567 | 0.008434 | 267 | 124 | 563 | 0 |
| `CP-U` | 36 | 0.222 | 0.778 | 0.007285 | 0.001214 | 0.014249 | 25800 | 1500 | 52500 | 714 |
| `IP-K` | 36 | 0.722 | 0.278 | 0.002254 | 0.000508 | 0.004486 | 66.8 | 30.6 | 150 | 0 |
| `IP-U` | 36 | 0.667 | 0.333 | 0.001715 | 0.000495 | 0.002960 | 54.5 | 33.3 | 108 | 4660 |

该表的关键现象非常直接：当前 full-flow 下，`IP-K / IP-U` 在 hit-rate、angle RMSE 与 `fdRef` RMSE 上都优于 `CP-K / CP-U`。这与论文主线中希望展示的“连续相位利用跨帧信息带来增益”并不一致，因此这版 full-flow 结果不能作为 CP/IP 论文性能图。

### 按帧数汇总

| P | case | truth-tooth hit | angle RMSE (deg) | fdRef RMSE (Hz) | fdRate RMSE (Hz/s) |
|---:|---|---:|---:|---:|---:|
| 8 | `CP-K` | 0.250 | 0.004937 | 219 | 0 |
| 8 | `CP-U` | 0.250 | 0.007983 | 5587 | 78.3 |
| 8 | `IP-K` | 0.750 | 0.002577 | 45.6 | 0 |
| 8 | `IP-U` | 0.667 | 0.002063 | 47.3 | 5629 |
| 10 | `CP-K` | 0.500 | 0.004123 | 148 | 0 |
| 10 | `CP-U` | 0.333 | 0.007815 | 16233 | 1079 |
| 10 | `IP-K` | 0.667 | 0.002471 | 56.9 | 0 |
| 10 | `IP-U` | 0.583 | 0.001933 | 61.1 | 4709 |
| 20 | `CP-K` | 0.250 | 0.004277 | 380 | 0 |
| 20 | `CP-U` | 0.083 | 0.005867 | 41202 | 600 |
| 20 | `IP-K` | 0.750 | 0.001581 | 89.9 | 0 |
| 20 | `IP-U` | 0.750 | 0.000913 | 54.2 | 3367 |

随着 `P` 增加到 20，`IP-U` 的 angle RMSE 继续下降到 `9.13e-4 deg`，但 `CP-U` 的 truth-tooth hit-rate 反而降到 `0.083`，`fdRef` RMSE 增大到 `4.12e4 Hz`。这说明当前 CP-U 不是在长窗口中稳定利用连续相位，而是更容易被 comb tooth / release basin 拉到远齿上。

### 高 SNR / 长窗口代表格点

`P=20, SNR=15 dB` 是最接近论文性能图直觉的高 SNR 长窗口格点，但当前结果仍不支持 CP 性能图：

| case | truth-tooth hit | angle RMSE (deg) | angle p95 (deg) | fdRef RMSE (Hz) | fdRate RMSE (Hz/s) |
|---|---:|---:|---:|---:|---:|
| `CP-K` | 0.667 | 0.003223 | 0.004465 | 203 | 0 |
| `CP-U` | 0.000 | 0.006823 | 0.010294 | 31235 | 1017 |
| `IP-K` | 1.000 | 0.000269 | 0.000291 | 12.7 | 0 |
| `IP-U` | 1.000 | 0.000269 | 0.000291 | 13.6 | 1831 |

如果直接把该格点画进论文，会得到“IP 明显优于 CP”的错误叙事。更合理的解释是：当前 `scanMfCpIpPerfMap` 测到的是完整 flow 失败模式，而不是纯 CP/IP 模型差异。

## 可观察现象

### 1. full-flow CP 结果仍被 wrong-tooth 主导

`CP-K` 的整体 truth-tooth hit-rate 只有 `0.333`，`CP-U` 只有 `0.222`。这说明 CP 路径仍经常没有停在 center tooth。对于 CP-K，`fdRate` 已知且 `fdRate RMSE=0`，但 `fdRef` RMSE 仍达到 `267 Hz`，说明问题不是 unknown-rate release 本身，而是 `fdRef` tooth / same-tooth basin 没有被 flow 稳定接住。

### 2. CP-U 的中位数和 RMSE 分离，说明存在远齿 tail

`CP-U` 的 angle median 只有 `0.001214 deg`，但 angle RMSE 为 `0.007285 deg`；`fdRef` median 为 `1500 Hz`，但 `fdRef` RMSE 达到 `25800 Hz`，p95 达到 `52500 Hz`。这不是平滑噪声导致的性能退化，而是少数远 tooth / bad-basin 样本把 RMSE 和 p95 拉坏。

典型远齿样本包括：

- `P=20, SNR=10 dB, seed=254`，`CP-U fdRefAbsErrHz=100443 Hz`，`toothIdx=-134`；
- `P=20, SNR=5 dB, seed=254`，`CP-U fdRefAbsErrHz=57000 Hz`，`toothIdx=76`；
- `P=10, SNR=15 dB, seed=253`，`CP-U fdRefAbsErrHz=50988 Hz`，`toothIdx=68`。

这些样本说明 CP-U 的 release 可以落到非常远的 comb tooth，而不是只在 `±1/T_f` 邻近齿间摇摆。

### 3. IP 在当前实现里看起来更好，但不能解释为 IP 理论上更优

`IP-K / IP-U` 的 hit-rate 和 RMSE 明显优于 CP，但这不应直接写成“IP 模型优于 CP 模型”。IP baseline 切断了跨帧公共相位 tying，它可能弱化了当前 flow 中的 comb / bad-basin 锁定，使局部优化更容易得到小 angle error；但它同时也不是本文主物理模型，不能替代 CP 的机制论证。

因此，这个结果更适合写成：当前 full-flow CP 路径仍有选齿和 same-tooth basin 风险；在该风险未隔离前，CP/IP 性能图会被算法 flow 失败污染。

### 4. SNR 提高不能自动修复 CP-U

`CP-U` 在 `SNR=15 dB` 下的整体表现并没有随 SNR 单调变好。尤其：

- `P=10, SNR=15 dB`：`truth-tooth hit=0`，`angle RMSE=0.0130 deg`，`fdRef RMSE=31365 Hz`；
- `P=20, SNR=15 dB`：`truth-tooth hit=0`，`angle RMSE=0.00682 deg`，`fdRef RMSE=31235 Hz`。

这说明当前瓶颈不是低 SNR 噪声主导，而是结构性 tooth / basin 问题。继续增加 SNR 或 repeat，不会把这个 scan 变成可用的 CP/IP 性能图。

### 5. 该结果与 CP/IP tying mechanism scan 一致但分工不同

`scanMfCpIpTying` 已经说明：CP continuous 在 truth center 上仍把 center tooth 作为最低点，但 flow / final center 一旦落到 wrong tooth，CP objective 也会稳定支持该错齿。因此，`scanMfCpIpPerfMap` 当前失败并不推翻 CP 物理模型；它只是说明当前完整 flow 仍没有可靠地把 CP 解带到正确 tooth / 正确 DoA basin。

## 当前结论

当前代表性结果应写成：**`scanMfCpIpPerfMap` 工程上已经跑通，但科学结论上不能作为 CP/IP 论文性能图。它显示当前完整 dynamic flow 下 CP-K / CP-U 仍被 wrong-tooth 与 bad-basin 主导，导致 IP baseline 在表观 RMSE 和 hit-rate 上反而更好。**

因此，本文档不把该 snapshot 作为“CP 优于 IP”的证据，而把它作为一个负结果 / gating evidence：在进入正式 CP/IP performance map 前，必须先完成受控 in-tooth 对比或先让 subset / rescue / same-tooth refine 更稳定。

## 对代码、replay 和论文图的影响

- `scanMfCpIpPerfMap.m` 可以继续保留为 full-flow CP/IP stress test，但当前结果不应进入论文主图。
- 下一步如果要服务论文 CP/IP 主线，应新增或运行受控入口，例如 `scanMfCpIpInToothPerfMap.m`：固定 correct tooth / truth-centered half-tooth `fdRef` range 后，再比较 `CP-K / CP-U / IP-K / IP-U` 的 SNR 与 frame-count 性能。
- 若目标是修 flow，应回到 `scanMfSubsetBankCoverage`、`replayMfPeriodicVsSubsetToothSelect`、`replayMfInToothTailCaseDiagnose` 与 conditional same-tooth refine，而不是继续扩大 `scanMfCpIpPerfMap` 的 repeat。
- 若后续仍使用 full-flow CP/IP 图，必须先保证 CP 路径的 truth-tooth hit-rate 与 same-tooth tail 不再主导 RMSE；否则图会把 flow 失败误写成 phase-model 差异。

## 后续建议

1. 不继续扩大本 snapshot 的 repeat。当前 `numRepeat=3` 已足以暴露 CP full-flow 失败模式；扩大 repeat 只会更稳定地证明该负结果。
2. 单开受控 in-tooth scan，而不是把 in-tooth 逻辑塞进本脚本。full-flow 与 controlled in-tooth 是两种不同口径，混在同一脚本会污染结果解释。
3. 后续如果修改 `scanMfCpIpPerfMap.m`，优先做 warning compact summary 和结果口径标注，不改 estimator 默认路径、不吞掉主核 warning。
4. 若 flow 稳定后重跑本 scan，应把当前 snapshot 保留为 superseded 负结果，新的 current snapshot 必须重新说明 truth-tooth hit-rate、angle RMSE 和 `fdRef` tail 是否已经恢复到可解释状态。

## 历史 / superseded snapshots

- 当前未保留更早的 superseded snapshot。前两次运行中的 frame subset 越界与 `msKnownDoaHalfWidth` 缺字段属于脚本合规清理过程中的工程错误，已修复后才产生本次 current snapshot。
