# scanMfCpIpInToothPerfMap 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing / representative` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfCpIpInToothPerfMap_20260428-213853.mat` |
| 当前一句话结论 | 在 truth-centered in-tooth 控制下，CP-K / CP-U 稳定保持 center tooth、低 `fdRef` 误差和高 non-ref coherence；IP 的 angle RMSE 经常更小，但代价是 `fdRef` consistency 与跨星相干一致性明显变差。 |
| 论文图定位 | `main figure / controlled trade-off figure`，用于 CP/IP angle-vs-Doppler-consistency 对照。 |
| 决策影响 | 固定为当前 controlled CP/IP 代表性结果；full-flow 性能仍需另看 `scanMfCpIpPerfMap` 或后续 resolved full-flow。 |
| 下一步动作 | 论文图不要只画 angle RMSE；至少并排画 angle RMSE 与 `fdRef` RMSE，可补 non-ref coherence 或 truth-tooth hit-rate。 |
| 禁止误用 | 不能把受控 in-tooth / oracle range 结果解释为默认 full-flow 已通过；也不能写成 CP 在 angle 上全面优于 IP。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfCpIpInToothPerfMap.m`
- 结果文档：`test/dev/scan/results/scanMfCpIpInToothPerfMap.md`
- scan 类型：`controlled / oracle scan`，也是 CP/IP paper-facing trade-off scan。
- 主要问题：在剥离 full-flow tooth selection 与 candidate adoption 污染后，`CP-K / CP-U / IP-K / IP-U` 的 angle、`fdRef` 与 coherence trade-off 如何变化。
- 扫描对象：`P`、`SNR`、`seed`、`case=CP-K/CP-U/IP-K/IP-U`。
- 不覆盖范围：不验证 default full-flow；不验证 subset selector；不证明 same-tooth rescue 已进入默认；不替代 known/unknown-rate CRB information-loss scan。
- truth 使用口径：`oracle-controlled`。truth 用于构造 half-tooth `fdRef` 搜索范围、truth-local unknown-rate 范围和 strict truth-tooth offline evaluation；truth 不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Yes，但必须标注 controlled in-tooth 条件。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `CP-K` | continuous-phase，known `fdRate`。 | Oracle range | 保留跨帧公共相位 tying 的 known-rate 受控表现。 | 不代表 full-flow CP-K。 |
| `CP-U` | continuous-phase，unknown `fdRate`，truth-local rate range。 | Oracle range | 保留 CP tying，并在局部 rate range 内 release `fdRate`。 | 不证明 unknown-rate 没有理论信息损失。 |
| `IP-K / IP-U` | independent-phase baseline。 | Oracle range | 切断跨帧公共相位约束后的对比。 | IP angle 更小不能解释成物理模型更好。 |
| `truthToothHitRate` | `toothIdx==0 && abs(toothResidualHz)<=50 Hz` 的比例。 | Eval only | 判断结果是否仍在 center tooth 内。 | 不能迁移到 runtime selector。 |
| `nonRefCoherenceFloorMedian` | 非参考星 coherence floor 的中位数。 | Eval only | 反映跨星链路相干一致性。 | 不能单独作为 estimator correctness。 |
| `fdRefRmseHz` | 参考星 `fdRef` 误差 RMSE。 | Eval only | 本 scan 中评价 Doppler consistency 的核心指标。 | 不要只看 angle RMSE 后忽略 `fdRef`。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：本 scan 的全部 controlled in-tooth samples；不是 full-flow samples。
- `resolved-sample`：可以把 strict truth-tooth rows 作为 offline resolved 近似；但这是 oracle/evaluation 口径，不能进入 runtime。
- `outlier rate`：可由 `1 - truthToothHitRate` 或 tail rows 辅助解释；但不作为 full-flow outlier rate。
- `truth-tooth / oracle range`：本 scan 的核心控制条件，只用于受控上限和离线评价。
- `stress-test`：否；full-flow stress 由 `scanMfCpIpPerfMap.md` 承担。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpInToothPerfMap_20260428-213853.mat` | 2026-04-28 | `representative` | `baseSeed=253`；`seedList=253:262`；`numRepeat=10`；`SNR=[0,5,10,15] dB`；`P=[8,10,20]`；`T_f=1/750 s`；`fdRef` half-tooth fraction `0.49`；`fdRate` half-width `1000 Hz/s`；static DoA local half-width `[0.002,0.002] deg`；truth-tooth residual tol `50 Hz`。 | CP 保住 center tooth、`fdRef` 与 coherence；IP angle 更小但 `fdRef` 与 coherence 明显退化。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:262`
- `numRepeat = 10`
- `snrDbList = [0, 5, 10, 15]`
- `frameCountList = [8, 10, 20]`
- `oracleFdHalfToothFraction = 0.49`
- `oracleFdRateHalfWidthHzPerSec = 1000`
- `staticLocalDoaHalfWidthDeg = [0.002, 0.002]`
- `toothResidualTolHz = 50`
- 关键 scan 轴：`P x SNR x seed x case`
- 关键 resolved / outlier 判据：strict truth-tooth rule `toothIdx == 0 && |toothResidualHz| <= 50 Hz`
- checkpoint：本次从 checkpoint 恢复，已有 `107/120` 个 task 完成，补跑 `13` 个 task 后清理 checkpoint 目录。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录顶层 `elapsedSec`；旧结果文档记录 average runtime 主要由 unknown-rate 分支主导。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`scanName`、`runKey`、`utcRun`、`config`、`perfTable`、`aggregateTable`、`repeatOutCell`、`checkpointSummaryTable`、`checkpointCleanupReport`、`plotData`
- task / case 规模：`120` 个 task，每个 task 含 `CP-K / CP-U / IP-K / IP-U` 四个 case，共 `480` 条 case 结果。
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、完整 fixture cache、transition bundle、全量 objective map 或图片文件。
- warning / fail 计数：运行日志中出现 unknown warm-anchor / `fmincon` near-singular warning；旧 task-level `warningSeen` 标记显示 `1` 个 task 有 warning。该 warning 未导致 hard fail。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

下表对全部 `120` 个 task、`480` 个 case 按 case 汇总。`truthToothHitRate` 使用 strict rule：`toothIdx==0 && abs(toothResidualHz)<=50 Hz`。

| case | samples | truthToothHitRate | angleRmseDeg | angleMedianDeg | angleP95Deg | fdRefRmseHz | fdRefMedianHz | fdRefP95Hz | fdRateRmseHzPerSec | nonRefCoherenceFloorMedian |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CP-K | 120 | 1.000 | 0.00400572 | 0.00155925 | 0.00846097 | 2.13575 | 0.0233542 | 5.0714 | 0 | 0.999246 |
| CP-U | 120 | 1.000 | 0.00364185 | 0.00101753 | 0.00819631 | 4.90411 | 0.0238246 | 5.1346 | 101.011 | 0.999947 |
| IP-K | 120 | 0.600 | 0.00252741 | 0.000677988 | 0.00568758 | 96.3658 | 37.9316 | 200.888 | 0 | 0.0749021 |
| IP-U | 120 | 0.600 | 0.00252741 | 0.000678000 | 0.00568758 | 96.3886 | 37.5878 | 200.888 | 645.362 | 0.0749021 |

读法：CP-K / CP-U 的 strict truth-tooth hit-rate 均为 `1.0`，并且 `fdRef` RMSE 只有 `2.14 / 4.90 Hz`；IP 的 angle RMSE 更小，但 `fdRef` RMSE 约 `96 Hz`，non-ref coherence floor median 约 `0.075`，明显低于 CP。

### 5.2 按扫描轴汇总

#### 按联合帧数汇总

| P | case | truthToothHitRate | angleRmseDeg | fdRefRmseHz | fdRateRmseHzPerSec | nonRefCoherenceFloorMedian |
|---:|---|---:|---:|---:|---:|---:|
| 8 | CP-K | 1.000 | 0.00352573 | 2.33891 | 0 | 0.999990 |
| 8 | CP-U | 1.000 | 0.00306159 | 2.34551 | 38.6061 | 0.999992 |
| 8 | IP-K | 0.525 | 0.00215884 | 103.472 | 0 | 0.0749021 |
| 8 | IP-U | 0.525 | 0.00215885 | 103.460 | 418.318 | 0.0749021 |
| 10 | CP-K | 1.000 | 0.00436924 | 2.57335 | 0 | 0.999933 |
| 10 | CP-U | 1.000 | 0.00379811 | 2.47861 | 23.668 | 0.999992 |
| 10 | IP-K | 0.475 | 0.00282054 | 112.920 | 0 | 0.111337 |
| 10 | IP-U | 0.475 | 0.00282054 | 112.872 | 651.611 | 0.111338 |
| 20 | CP-K | 1.000 | 0.00407632 | 1.26161 | 0 | 0.200113 |
| 20 | CP-U | 1.000 | 0.00399880 | 7.77855 | 168.994 | 0.993796 |
| 20 | IP-K | 0.800 | 0.00255879 | 66.3457 | 0 | 0.00889216 |
| 20 | IP-U | 0.800 | 0.00255877 | 66.5447 | 806.158 | 0.00888769 |

#### 按 SNR 汇总

| SNR | case | truthToothHitRate | angleRmseDeg | fdRefRmseHz | fdRateRmseHzPerSec | nonRefCoherenceFloorMedian |
|---:|---|---:|---:|---:|---:|---:|
| 0 | CP-K | 1.000 | 0.00680309 | 2.35351 | 0 | 0.994837 |
| 0 | CP-U | 1.000 | 0.00638820 | 2.34764 | 38.0618 | 0.997861 |
| 0 | IP-K | 0.266667 | 0.00464391 | 160.527 | 0 | 0.0743753 |
| 0 | IP-U | 0.266667 | 0.00464392 | 160.596 | 683.109 | 0.0743794 |
| 5 | CP-K | 1.000 | 0.00331680 | 0.607910 | 0 | 0.999996 |
| 5 | CP-U | 1.000 | 0.00301335 | 0.613142 | 18.0325 | 0.999997 |
| 5 | IP-K | 0.433333 | 0.00190190 | 93.2059 | 0 | 0.0747064 |
| 5 | IP-U | 0.433333 | 0.00190187 | 93.1629 | 707.034 | 0.0747064 |
| 10 | CP-K | 1.000 | 0.00241819 | 0.122415 | 0 | 0.999999 |
| 10 | CP-U | 1.000 | 0.00152417 | 0.117748 | 23.8696 | 0.999999 |
| 10 | IP-K | 0.733333 | 0.000547059 | 45.7942 | 0 | 0.0748935 |
| 10 | IP-U | 0.733333 | 0.000547089 | 45.8449 | 682.849 | 0.0748894 |
| 15 | CP-K | 1.000 | 0.00102586 | 3.51029 | 0 | 0.995824 |
| 15 | CP-U | 1.000 | 0.000916467 | 9.50262 | 196.136 | 0.996741 |
| 15 | IP-K | 0.966667 | 0.000262280 | 24.3346 | 0 | 0.0749990 |
| 15 | IP-U | 0.966667 | 0.000262204 | 24.3101 | 482.859 | 0.0749990 |

#### IP/CP trade-off 汇总

| rate | median angle IP/CP | min angle IP/CP | max angle IP/CP | median fdRef IP/CP | min fdRef IP/CP | max fdRef IP/CP |
|---|---:|---:|---:|---:|---:|---:|
| known | 0.410717 | 0.178603 | 0.707137 | 115.337 | 4.14341 | 4052.36 |
| unknown | 0.607466 | 0.243505 | 0.900124 | 115.139 | 0.635691 | 4444.71 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Angle RMSE surface / slices | `SNR` 或 `P` | `angleRmseDeg` | `CP-K / CP-U / IP-K / IP-U` | Yes | 必须和 `fdRef` 或 coherence 图一起解释。 |
| `fdRef` RMSE surface / slices | `SNR` 或 `P` | `fdRefRmseHz` | `CP-K / CP-U / IP-K / IP-U` | Yes | CP/IP trade-off 的核心图，避免只看 angle。 |
| truth-tooth hit rate | `SNR` 或 `P` | `truthToothHitRate` | 四个 case | Appendix / diagnostic | 这是 offline truth evaluation，不进入 runtime selector。 |
| non-ref coherence | `SNR` 或 `P` | `nonRefCoherenceFloorMedian` | 四个 case | Mechanism figure | 用于解释 CP 保持相干一致性，IP 切断相位 tying。 |
| same-tooth tail table | `seed / P / SNR` | `fdRefAbsErrHz` | CP-U tail samples | Diagnostic | 不应膨胀成主性能图。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- in-tooth 口径隔离了 full-flow wrong-tooth 污染：CP-K / CP-U 的 strict truth-tooth hit-rate 均为 `1.0`。
- CP 的主要优势是 `fdRef` consistency 与跨星相干一致性：CP 的 `fdRef` RMSE 为 `2.14 / 4.90 Hz`，IP 约 `96 Hz`；CP non-ref coherence median 接近 `1`，IP 约 `0.075`。
- IP angle 更小但代价明显：IP angle RMSE 约为 CP 的 `0.41`（known）或 `0.61`（unknown），但 IP `fdRef` RMSE 中位数约为 CP 的 `115` 倍。
- 高 SNR / 长窗口格点 `P=20, SNR=15 dB` 中，四个 case 都在 strict truth-tooth 内，仍体现 angle-vs-Doppler consistency trade-off。

### 6.2 反向、污染或未解决现象

- IP 的 angle RMSE 经常更小，说明论文不能写成“CP angle 全面优于 IP”。
- CP-U 在 high-SNR / long-window 下仍有 same-tooth tail；例如 `P=20, SNR=15, seed=260` 的 `fdRefAbsErrHz=48.95 Hz`，接近 strict residual 门限。
- CP-U 相比 CP-K 没有明显 angle loss，不等价于 unknown-rate 没有理论信息损失；正式信息损失仍以 `scanMfKnownUnknownInformationLoss` 为准。
- unknown-rate 分支 runtime 成本明显高，当前 snapshot 不建议继续盲目扩大 repeat。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `P=20, SNR=15 dB, seed=260, CP-U` | same-tooth tail | `fdRefAbsErrHz=48.9535 Hz`，`fdRateAbsErrHzPerSec=742.978`，`toothIdx=0` | 解释 CP-U 高 SNR 下 `fdRef` RMSE 被少数 tail 拉大；不是远 tooth failure。 |
| `P=20, SNR=15 dB` | high-SNR long-window slice | IP angle RMSE 约 `2.25e-4 deg`，CP-U `fdRef` median 约 `1.08 Hz`，IP coherence floor 约 `0.00885` | 说明 angle-only 指标不足以评价 CP/IP 模型层级。 |
| IP overall | controlled baseline | truth-tooth hit-rate `0.6`，`fdRef` RMSE 约 `96 Hz` | 说明 IP 独立相位自由度会牺牲 Doppler consistency。 |

## 7. 机制解释

### 7.1 当前解释

CP 模型保留同一颗卫星跨帧公共相位 tying，因此它更强地约束 reference `fdRef` 与 non-ref-sat 相干链。受控 in-tooth 后，CP 不再被 full-flow wrong-tooth 污染，所以可以稳定保持 center tooth 和较低 `fdRef` 误差。

IP 模型为每个 `(sat, frame)` 引入独立相位，自由度更高，可能在 angle-only 局部拟合上更灵活，因此 angle RMSE 经常更小；但这种灵活性会切断跨帧相位一致性，导致 `fdRef` 与 non-ref coherence 明显退化。因此当前结果应解释为 angle-vs-Doppler-consistency trade-off，而不是简单胜负。

### 7.2 这个 scan 支持什么

- 支持 controlled in-tooth 条件下 CP 模型在 `fdRef` consistency 与 non-ref coherence 上显著优于 IP。
- 支持 full-flow CP/IP 负结果需要与受控 in-tooth 结果分开解释。
- 支持论文图不能只看 angle RMSE，应同时报告 `fdRef` 和 coherence / consistency 指标。
- 支持 CP/IP 的论文口径写成 angle-vs-Doppler-consistency trade-off。

### 7.3 这个 scan 不证明什么

- 不证明 default full-flow 已修复。
- 不证明 CP 在 angle 上全面优于 IP。
- 不证明 unknown-rate 没有理论信息损失。
- 不证明 truth-centered oracle range 可以进入 runtime selector、gate、candidate adoption 或 final winner。
- 不证明 same-tooth tail 已解决。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 只使用受控范围评价 CP/IP trade-off。 |
| flow 默认路径 | 不改；不能把 truth-centered in-tooth 结果迁移为 full-flow selector。 |
| replay 下一步 | 若继续追 tail，可固定 `P=20, SNR=15, seed=260` 做 same-tooth residual / rate-profile replay。 |
| regression | 不写；当前是 controlled scan 结果，不是自动契约。 |
| 论文图 | `main` 候选，但必须标注 controlled in-tooth / oracle range。 |
| 排障记录 | 主记录保留“controlled CP/IP trade-off 已固定；full-flow stress 另看”的结论即可。 |

## 9. 限制与禁止解释

- 不要把 truth-centered half-tooth range 用作 runtime fd range selector。
- 不要把 strict truth-tooth hit-rate 用作 default winner adoption gate。
- 不要只画 angle RMSE 后声称 IP 或 CP 全面更好。
- 不要用该 scan 替代 known/unknown-rate CRB / EFIM information-loss scan。
- 不要把 CP-U 的受控表现解释成 full-flow CP-U 没有 bad-basin / tail。
- 不要把 high-SNR same-tooth tail 当作 wrong-tooth 问题处理。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfCpIpInToothPerfMap_20260428-213853.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfCpIpInToothPerfMap.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：当前脚本默认 `numRepeat` 可能已经改小用于 smoke；复现实验时应以本 snapshot 的 `numRepeat=10` 配置为准。

## 11. 历史备注

- 当前只绑定 `scanMfCpIpInToothPerfMap_20260428-213853.mat` 作为 controlled in-tooth 代表性 snapshot。
- `scanMfCpIpPerfMap.md` 的 full-flow 结果应继续作为 stress-test 单独解释；不要把两者合并成一组无边界 CP/IP 结论。
