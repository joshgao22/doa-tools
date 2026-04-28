# scanMfCpIpInToothPerfMap 结果记录

## 对应 scan

- `test/dev/scan/scanMfCpIpInToothPerfMap.m`

## 扫描目标

在受控 in-tooth 条件下比较 `CP-K / CP-U / IP-K / IP-U` 随 SNR 与联合帧数变化的性能。该 scan 使用 truth-centered half-tooth `fdRef` 搜索范围、truth-local unknown-rate 范围，以及同一 `MS-SF-Static` DoA seed，刻意隔离 full-flow tooth selection、candidate adoption 和 same-tooth rescue 的影响。

因此，该结果只代表 **given-correct-tooth / controlled in-tooth** 条件下的模型与求解上限，不代表当前完整 dynamic flow 的系统性能。它应与 `scanMfCpIpPerfMap.md` 中的 full-flow 负结果分开解释。

## Snapshot index

| snapshot | 状态 | 配置 | 结论 |
|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpInToothPerfMap_20260428-213853.mat` | current | `baseSeed=253`，`seedList=253:262`，`SNR=[0,5,10,15] dB`，`P=[8,10,20]`，`T_f=1/750 s`，每格 `numRepeat=10`，`fdRef` oracle half-tooth fraction `0.49`，`fdRate` oracle half-width `1000 Hz/s`，static DoA local half-width `[0.002,0.002] deg`，truth-tooth 判据为 `toothIdx==0 && abs(toothResidualHz)<=50 Hz` | 当前代表性结果。CP-K / CP-U 在所有 `P x SNR` 格点中 strict truth-tooth hit-rate 均为 `1.0`，说明 in-tooth 口径成功隔离了 full-flow wrong-tooth 污染；CP 在 `fdRef` consistency 和 non-ref coherence 上显著优于 IP，但 IP 在当前 estimator / metric 下经常给出更小 angle RMSE。 |

## 存档数据检查

`scanMfCpIpInToothPerfMap_20260428-213853.mat` 为 `saveExpSnapshot` 生成的轻量 scan snapshot：

- 顶层变量为 `data / meta / inventory`。
- `data` 中只保存 `scanData`。
- `scanData` 包含 `config`、`perfTable`、`aggregateTable`、`repeatOutCell`、`checkpointSummaryTable`、`checkpointCleanupReport` 与 `plotData`。
- `config` 中记录 `baseSeed=253`、`numRepeat=10`、`seedList=253:262`、`snrDbList=[0,5,10,15]`、`frameCountList=[8,10,20]`、`oracleFdHalfToothFraction=0.49`、`oracleFdRateHalfWidthHzPerSec=1000`、`staticLocalDoaHalfWidthDeg=[0.002,0.002]`、`toothResidualTolHz=50`。
- `repeatOutCell` 包含 `120` 个 task，每个 task 对应一个 `(P, SNR, seed)`，每个 task 内有 `CP-K / CP-U / IP-K / IP-U` 四个 case；共 `480` 条 case 结果。
- snapshot 中未保存大体量 `rxSigCell`、完整 `sceneSeq`、完整 fixture cache、transition bundle 或图片文件，符合 scan cache 的轻量保存口径。

运行日志显示本次从 checkpoint 恢复，已有 `107/120` 个 task 完成，最终补跑 `13` 个 task 后清理 checkpoint 目录并保存 snapshot。日志中仍出现 unknown warm-anchor / `fmincon` 近奇异 warning；snapshot 的 task 级 `warningSeen` 标记显示有 `1` 个 task 出现 warning。该 warning 没有导致 hard fail，但后续若继续维护该 scan，可在 scan 层补 compact warning count，不应在 estimator 主核中吞 warning。

## 当前代表性结果

### 全部样本汇总

下表对全部 `120` 个 task、`480` 个 case 按 case 汇总。`truth-tooth hit` 使用当前 scan 的 strict 判据：`toothIdx==0 && abs(toothResidualHz)<=50 Hz`。

| case   |   samples |   truthToothHitRate |   angleRmseDeg |   angleMedianDeg |   angleP95Deg |   fdRefRmseHz |   fdRefMedianHz |   fdRefP95Hz |   fdRateRmseHzPerSec |   nonRefCoherenceFloorMedian |
|:-------|----------:|--------------------:|---------------:|-----------------:|--------------:|--------------:|----------------:|-------------:|---------------------:|-----------------------------:|
| CP-K   |       120 |                 1   |     0.00400572 |      0.00155925  |    0.00846097 |       2.13575 |       0.0233542 |       5.0714 |                0     |                    0.999246  |
| CP-U   |       120 |                 1   |     0.00364185 |      0.00101753  |    0.00819631 |       4.90411 |       0.0238246 |       5.1346 |              101.011 |                    0.999947  |
| IP-K   |       120 |                 0.6 |     0.00252741 |      0.000677988 |    0.00568758 |      96.3658  |      37.9316    |     200.888  |                0     |                    0.0749021 |
| IP-U   |       120 |                 0.6 |     0.00252741 |      0.000678    |    0.00568758 |      96.3886  |      37.5878    |     200.888  |              645.362 |                    0.0749021 |

该表的主要读法是：

- `CP-K / CP-U` 的 strict truth-tooth hit-rate 均为 `1.0`，说明受控 `fdRef` half-tooth 范围成功把 CP 解锁在 center tooth 内。
- `IP-K / IP-U` 的整体 truth-tooth hit-rate 只有 `0.6`，且 `fdRef` RMSE 约 `96 Hz`，明显差于 CP 的 `2.14 Hz / 4.90 Hz`。
- IP 的 angle RMSE 约 `0.00253 deg`，小于 CP 的 `0.00401 deg / 0.00364 deg`。因此不能把该结果写成“CP 在 angle 上全面优于 IP”；更合理的解释是 IP 用独立相位自由度换取局部角度拟合，但丢失 `fdRef` consistency 和跨星相干一致性。

### 按联合帧数汇总

|   P | case   |   truthToothHitRate |   angleRmseDeg |   fdRefRmseHz |   fdRateRmseHzPerSec |   nonRefCoherenceFloorMedian |
|----:|:-------|--------------------:|---------------:|--------------:|---------------------:|-----------------------------:|
|   8 | CP-K   |               1     |     0.00352573 |       2.33891 |               0      |                   0.99999    |
|   8 | CP-U   |               1     |     0.00306159 |       2.34551 |              38.6061 |                   0.999992   |
|   8 | IP-K   |               0.525 |     0.00215884 |     103.472   |               0      |                   0.0749021  |
|   8 | IP-U   |               0.525 |     0.00215885 |     103.46    |             418.318  |                   0.0749021  |
|  10 | CP-K   |               1     |     0.00436924 |       2.57335 |               0      |                   0.999933   |
|  10 | CP-U   |               1     |     0.00379811 |       2.47861 |              23.668  |                   0.999992   |
|  10 | IP-K   |               0.475 |     0.00282054 |     112.92    |               0      |                   0.111337   |
|  10 | IP-U   |               0.475 |     0.00282054 |     112.872   |             651.611  |                   0.111338   |
|  20 | CP-K   |               1     |     0.00407632 |       1.26161 |               0      |                   0.200113   |
|  20 | CP-U   |               1     |     0.0039988  |       7.77855 |             168.994  |                   0.993796   |
|  20 | IP-K   |               0.8   |     0.00255879 |      66.3457  |               0      |                   0.00889216 |
|  20 | IP-U   |               0.8   |     0.00255877 |      66.5447  |             806.158  |                   0.00888769 |

帧数维度上，CP 的 `fdRef` consistency 明显稳定：`P=8/10/20` 下 CP-K 的 `fdRef` RMSE 分别为 `2.34 / 2.57 / 1.26 Hz`，CP-U 分别为 `2.35 / 2.48 / 7.78 Hz`。IP 的 `fdRef` RMSE 则长期停留在 `66 ~ 113 Hz` 量级。

需要注意的是，`P=20` 时 CP-K 的 non-ref coherence median 降到 `0.200`，CP-U 仍为 `0.994`；这说明 known-rate 与 unknown-rate 在长窗口下的 local basin / coherence 行为仍有差异，不能只看 angle RMSE。

### 按 SNR 汇总

|   SNR | case   |   truthToothHitRate |   angleRmseDeg |   fdRefRmseHz |   fdRateRmseHzPerSec |   nonRefCoherenceFloorMedian |
|------:|:-------|--------------------:|---------------:|--------------:|---------------------:|-----------------------------:|
|     0 | CP-K   |            1        |    0.00680309  |      2.35351  |               0      |                    0.994837  |
|     0 | CP-U   |            1        |    0.0063882   |      2.34764  |              38.0618 |                    0.997861  |
|     0 | IP-K   |            0.266667 |    0.00464391  |    160.527    |               0      |                    0.0743753 |
|     0 | IP-U   |            0.266667 |    0.00464392  |    160.596    |             683.109  |                    0.0743794 |
|     5 | CP-K   |            1        |    0.0033168   |      0.60791  |               0      |                    0.999996  |
|     5 | CP-U   |            1        |    0.00301335  |      0.613142 |              18.0325 |                    0.999997  |
|     5 | IP-K   |            0.433333 |    0.0019019   |     93.2059   |               0      |                    0.0747064 |
|     5 | IP-U   |            0.433333 |    0.00190187  |     93.1629   |             707.034  |                    0.0747064 |
|    10 | CP-K   |            1        |    0.00241819  |      0.122415 |               0      |                    0.999999  |
|    10 | CP-U   |            1        |    0.00152417  |      0.117748 |              23.8696 |                    0.999999  |
|    10 | IP-K   |            0.733333 |    0.000547059 |     45.7942   |               0      |                    0.0748935 |
|    10 | IP-U   |            0.733333 |    0.000547089 |     45.8449   |             682.849  |                    0.0748894 |
|    15 | CP-K   |            1        |    0.00102586  |      3.51029  |               0      |                    0.995824  |
|    15 | CP-U   |            1        |    0.000916467 |      9.50262  |             196.136  |                    0.996741  |
|    15 | IP-K   |            0.966667 |    0.00026228  |     24.3346   |               0      |                    0.074999  |
|    15 | IP-U   |            0.966667 |    0.000262204 |     24.3101   |             482.859  |                    0.074999  |

SNR 维度上，CP 的 angle RMSE 随 SNR 提高整体下降；但高 SNR 下 `fdRef` RMSE 并不完全单调。例如 `SNR=15 dB` 时 CP-U 的 `fdRef` RMSE 为 `9.50 Hz`，主要由少数高 SNR / 长窗口 tail 拉大，而中位数仍只有 `2.62 Hz` 左右。这个现象说明 in-tooth 条件下仍可能存在 same-tooth tail，但已不再是 full-flow 里的远 tooth 失败。

### 高 SNR / 长窗口代表格点

`P=20, SNR=15 dB` 是最接近论文性能图直觉的高 SNR 长窗口格点。

| case   |   truthToothHitRate |   angleRmseDeg |   angleP95Deg |   fdRefRmseHz |   fdRefMedianHz |   fdRefP95Hz |   fdRateRmseHzPerSec |   nonRefCoherenceFloorMedian |
|:-------|--------------------:|---------------:|--------------:|--------------:|----------------:|-------------:|---------------------:|-----------------------------:|
| CP-K   |                   1 |    0.000993794 |   0.00139346  |       2.36651 |         2.62972 |      2.73826 |                0     |                   0.596063   |
| CP-U   |                   1 |    0.000891987 |   0.00125147  |      15.5325  |         1.0777  |     28.0758  |              336.466 |                   0.997785   |
| IP-K   |                   1 |    0.000225447 |   0.00029742  |       9.80542 |         8.53638 |     15.3997  |                0     |                   0.00884915 |
| IP-U   |                   1 |    0.000225206 |   0.000297193 |       9.87384 |         8.87453 |     15.7724  |              632.408 |                   0.00884019 |

该格点的可解释性比 full-flow `scanMfCpIpPerfMap` 明显更好：四个 case 都在 strict truth-tooth 内。但它仍显示 angle 与 Doppler/coherence 存在 trade-off：

- IP-K / IP-U 的 angle RMSE 约 `2.25e-4 deg`，小于 CP-K / CP-U；
- CP-K / CP-U 的 non-ref coherence floor median 分别为 `0.596 / 0.998`，远高于 IP 的 `0.00885`；
- CP-U 的 `fdRef` RMSE 为 `15.53 Hz`，但 median 只有 `1.08 Hz`，说明仍有少数 same-tooth tail。

### IP/CP trade-off 汇总

下表按每个 `(P,SNR)` 格点计算 `IP / CP` 比值后汇总。`angle IP/CP < 1` 表示 IP angle RMSE 更小；`fdRef IP/CP > 1` 表示 IP 的 `fdRef` RMSE 更差。

| rate    |   median angle IP/CP |   min angle IP/CP |   max angle IP/CP |   median fdRef IP/CP |   min fdRef IP/CP |   max fdRef IP/CP |
|:--------|---------------------:|------------------:|------------------:|---------------------:|------------------:|------------------:|
| known   |             0.410717 |          0.178603 |          0.707137 |              115.337 |          4.14341  |           4052.36 |
| unknown |             0.607466 |          0.243505 |          0.900124 |              115.139 |          0.635691 |           4444.71 |

可以看到，IP 的 angle RMSE 中位数约为 CP 的 `0.41`（known）或 `0.61`（unknown），但 IP 的 `fdRef` RMSE 中位数约为 CP 的 `115` 倍。该 trade-off 是本次 snapshot 最重要的解释口径。

### CP-U same-tooth tail 样本

下表列出 CP-U 按 `fdRefAbsErrHz` 排序的最大若干样本。它们全部仍满足 `toothIdx=0`，因此不是 full-flow 远 tooth 失败，而是受控 in-tooth 条件下的 residual / same-tooth tail。

|   P |   SNR |   seed |   angleErrDeg |   fdRefAbsErrHz |   fdRateAbsErrHzPerSec |   toothIdx |   toothResidualHz |   nonRefCoherenceFloor |
|----:|------:|-------:|--------------:|----------------:|-----------------------:|-----------:|------------------:|-----------------------:|
|  20 |    15 |    260 |   0.00139637  |        48.9535  |           742.978      |          0 |          48.9535  |               0.991165 |
|  10 |     0 |    261 |   0.010503    |        10.8127  |            42.2927     |          0 |          10.8127  |               0.966951 |
|   8 |    15 |    256 |   0.000909207 |         6.86029 |             0.00492727 |          0 |           6.86029 |               0.999855 |
|   8 |    15 |    259 |   0.00117159  |         6.73661 |             0.00161175 |          0 |           6.73661 |               0.999935 |
|   8 |    15 |    255 |   0.000302118 |         6.59584 |             0.00143982 |          0 |           6.59584 |               0.992537 |
|  10 |    15 |    261 |   0.00027836  |         5.28493 |             0.084297   |          0 |           5.28493 |               0.999658 |

最大样本为 `P=20, SNR=15 dB, seed=260`，`fdRefAbsErrHz=48.95 Hz`，接近 `50 Hz` strict residual 门限，同时 `fdRateAbsErrHzPerSec=742.98`。这类样本可以解释 CP-U 高 SNR 下 `fdRef` RMSE 被拉大的原因；但它和 full-flow 中几十 kHz 级远 tooth tail 不是同一类问题。

## 可观察现象

### 1. in-tooth 口径成功隔离了 full-flow wrong-tooth 污染

`CP-K / CP-U` 在全部 `120` 个 task 中 strict truth-tooth hit-rate 都为 `1.0`。与 `scanMfCpIpPerfMap` 的 full-flow 负结果相比，这说明 CP 路径此前表现差的主要污染源确实来自 tooth selection / adoption / bad basin，而不是 CP 模型在正确 tooth 内完全不可用。

### 2. CP 的优势主要体现在 `fdRef` consistency 与跨星相干一致性

CP 的 `fdRef` RMSE 为 `2.14 Hz / 4.90 Hz`，而 IP 为约 `96 Hz`。同时 CP 的 non-ref coherence floor median 接近 `1`，IP 只有约 `0.075`。这说明 independent phase baseline 切断跨帧公共相位后，可能仍能局部拟合 angle，但难以保持参考 Doppler 与非参考星相干链的一致性。

因此，论文中如果使用该图，应避免只画 angle RMSE。更合适的图组应同时包含：

- angle RMSE；
- `fdRef` RMSE；
- truth-tooth hit-rate；
- non-ref coherence floor 或相近的相干一致性指标。

### 3. IP angle 更小不能解释为 IP 模型物理上更好

IP 的 angle RMSE 在多数格点更小，尤其高 SNR 下非常明显。但 IP 同时出现明显 `fdRef` degradation 和 non-ref coherence collapse。这更像是独立相位自由度吸收了跨帧相位不一致，使 angle 局部拟合更容易，而不是 IP 更准确地恢复了本文主参数。

该结果与论文主线并不矛盾：本文强调 CP 是主物理模型，因为它保留跨帧连续相位信息；IP 是切断该信息的对比基线。当前 snapshot 的正确叙事应是 **angle-vs-Doppler consistency trade-off**，而不是“CP angle 全面优于 IP”。

### 4. CP-U 相比 CP-K 没有明显信息损失，受控 rate 范围下甚至略好

整体上 CP-U 的 angle RMSE 为 `0.00364 deg`，略好于 CP-K 的 `0.00401 deg`；多数 SNR / P 格点中 CP-U 与 CP-K 接近，甚至 CP-U 稍优。这说明在 truth-local `±1000 Hz/s` 的 `fdRate` 范围内，unknown-rate 对 angle 主参数的额外损失并不明显。

这不等价于“unknown Doppler-rate 没有信息损失”。它只说明当前 estimator + truth-local controlled range 下，信息损失被局部 oracle 条件压低了。正式的 known/unknown 信息损失解释仍应以 `scanMfKnownUnknownInformationLoss` 的 CRB / EFIM 口径为主。

### 5. same-tooth tail 仍存在，但已经不是远 tooth tail

CP-U 最大 `fdRef` tail 只有 `48.95 Hz`，且仍满足 strict truth-tooth 判据。这和 full-flow 负结果中的几十 kHz 远 tooth tail 不同。后续若要继续优化该 scan，不应再围绕 wrong-tooth rescue，而应关注 same-tooth residual、fdRate local profile 与 high-SNR / long-window 局部解的稳定性。

### 6. runtime 成本主要来自 unknown-rate 分支

从 `runTimeMs` 看，CP-U 和 IP-U 明显慢于 known-rate case。整体平均 run time 约为：

- `CP-K`：`4.51 s`;
- `CP-U`：`74.44 s`;
- `IP-K`：`4.10 s`;
- `IP-U`：`33.47 s`.

这符合 unknown-rate warm-anchor / release 更重的预期。当前 snapshot 可作为受控结果存档，不建议为追求更大 repeat 继续扩大运行成本。

## 当前结论

当前代表性结果应写成：**在 truth-centered in-tooth 控制下，CP-K / CP-U 已经摆脱 full-flow wrong-tooth 污染，并能稳定保持 center tooth、低 `fdRef` 误差和较高 non-ref coherence；IP baseline 虽然在 angle RMSE 上经常更小，但代价是 `fdRef` consistency 和跨星相干一致性明显变差。**

因此，这个 snapshot 可以作为 CP/IP 论文主线的受控候选结果，但不能被表述为“CP 在所有指标上全面优于 IP”。更稳妥的论文口径是：

- CP 是保持跨帧连续相位与参考 Doppler 物理一致性的主模型；
- IP 是切断跨帧相位 tying 后的对比基线；
- 在受控 correct-tooth 条件下，CP 的 Doppler consistency 明显强于 IP；
- 当前实现中 IP 可获得更小 angle RMSE，说明 angle-only 指标不足以评价 CP/IP 模型层级。

## 对代码、replay 和论文图的影响

- `scanMfCpIpInToothPerfMap.m` 应继续保留为 controlled CP/IP performance-map 入口，与 `scanMfCpIpPerfMap.m` 的 full-flow stress test 分工明确。
- 本 snapshot 可作为当前代表性 in-tooth 结果固定，不建议继续扩大 repeat；继续扩大只会增加 unknown-rate 成本，并不会改变本次最重要的 trade-off 结论。
- 论文图候选不应只画 angle RMSE surface。建议至少并排画 angle RMSE 与 `fdRef` RMSE，必要时补 truth-tooth hit-rate 或 non-ref coherence floor。
- 若后续目标是完整系统性能图，应回到 `scanMfSubsetBankCoverage`、`replayMfPeriodicVsSubsetToothSelect` 和 same-tooth tail / polish replay，而不是把 full-flow 逻辑塞回这个受控 scan。
- 若后续目标是 known/unknown-rate 信息损失，应继续使用 `scanMfKnownUnknownInformationLoss` 的 CRB / EFIM 口径；本 scan 只说明 estimator 在 truth-local rate 范围内的受控表现。

## 后续建议

1. 先固定该 snapshot 为 current，不再对本结果做 repeat 扩展。
2. 若写论文图，采用“angle-vs-fdRef consistency trade-off”口径，不写 CP angle 全面优于 IP。
3. 后续如果要增强该 scan，优先增加派生 summary：truth-tooth-conditioned summary、angle-vs-fdRef Pareto summary、same-tooth residual tail summary；不要改 estimator 默认数值路径。
4. 若继续追 CP-U 的 high-SNR / long-window tail，可单独做小 replay 固定 `P=20, SNR=15, seed=260`，但它属于 same-tooth residual / rate-local-profile 诊断，不应混入本结果文档继续膨胀。
5. 后续 flow 稳定后，可重跑 `scanMfCpIpPerfMap` 作为 full-flow performance map；届时本 in-tooth snapshot 仍作为 controlled upper-bound / trade-off 对照保留。

## 历史 / superseded snapshots

- 当前未保留更早的 superseded snapshot。本文件只绑定 `scanMfCpIpInToothPerfMap_20260428-213853.mat` 作为 current。
