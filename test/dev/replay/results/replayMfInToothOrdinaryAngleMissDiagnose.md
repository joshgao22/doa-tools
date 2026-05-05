# replayMfInToothOrdinaryAngleMissDiagnose

## 对应 replay

- `test/dev/replay/replayMfInToothOrdinaryAngleMissDiagnose.m`

## 观察目标

在 truth-centered half-tooth `fdRef` 范围内，专门诊断 controlled in-tooth 条件下的 non-collapse ordinary DoA miss。该 replay 的目标不是继续扩大 hard-collapse rescue，而是回答：

1. raw `MS-MF-CP-U-in-tooth` 剩余 `angleErr > 0.002 deg` 的样本中，ordinary miss 占多少；
2. ordinary miss 是否只是 final 点附近的细小偏差，能否由 final-centered small polish 修复；
3. ordinary miss 是否仍需要从 `wide` 或 `single-MF` center 重新做 basin-entry；
4. 这些可实现 candidate 与 truth-DoA oracle 的差距有多大；
5. 是否可以据此继续放宽 hard-collapse gate，或需要另开 ordinary-wide gate / adoption 诊断。

该 replay 只验证 tooth 已正确时的 in-tooth DoA 机制，不做 subset tooth selection，不改变 estimator / flow / objective，也不进入 regression。truth 只用于 offline miss label、caseRole 与 truth-DoA oracle 上限；`default / wide / single-MF / final-small-polish` candidate 本身仍是可实现 probe。

## 最新代码口径

- 使用 truth-centered half-tooth `fdRef` 范围，剥离 wrong-tooth 因素。
- 默认 `numRepeat=100`、`baseSeed=253`、`snrDb=10`。
- `angleHitThresholdDeg=0.002`。
- `ordinaryCoherenceThreshold=0.950`，ordinary miss 要求 same tooth、fd healthy、non-ref coherence 不 collapse。
- `final-small-polish` 只做 final-centered DoA 小网格 probe，DoA steps 为 `[-0.002 -0.001 -0.0005 0 0.0005 0.001 0.002] deg`。
- basin-entry DoA steps 为 `[-0.006 -0.004 -0.002 0 0.002 0.004 0.006] deg`。
- 默认 `includeFamilySafeAdoptBank=true`、`includeJointSafeAdoptBank=false`、`includeLegacyFailedBanks=false`。
- Telegram 通知为旁路 I/O；成功保存 snapshot 后发送 HTML `DONE`，失败时发送 HTML `FAILED` 并 `rethrow`。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat` | `baseSeed=253`，`numRepeat=100`，`snrDb=10`，`ordinaryCoh=0.950`，`rescueCoh=0.200`，`safeStep=0.004 deg`，`familyStep=0.006 deg`，`joint=0`，`family=1`，`angleHitThreshold=0.002 deg` | ordinary miss 存在但占比不高；final-centered small polish 完全没有改善；ordinary miss 全部可由 `wide-basin-entry` 救回并接近 truth-DoA oracle；`single-MF` 对 ordinary miss 明显不稳；blanket wide 仍会放大整体 tail，因此下一步应找 truth-free ordinary-wide gate，而不是继续放宽 hard-collapse gate。 |

## 当前代表性结果

2026-05-05 的 100-repeat run 是当前唯一代表性结果。snapshot 保存为：

```text
 test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat
```

运行配置摘要：

| 配置项 | 值 |
|---|---:|
| repeats | 100 |
| SNR | 10 dB |
| base seed | 253 |
| seed range | 253:352 |
| repeat mode | parfor-auto |
| fd oracle half-tooth fraction | 0.490 |
| fdRate oracle half-width | 1000 Hz/s |
| ordinary coherence threshold | 0.950 |
| gated rescue coherence threshold | 0.200 |
| gated rescue phase gate | inactive |
| safe adopt max DoA step | 0.004 deg |
| family safe max DoA step | 0.006 deg |
| safe adopt coherence threshold | 0.950 |
| safe adopt max fdRef step | 300 Hz |
| angle hit threshold | 0.002 deg |

checkpoint 从空目录开始，100 个 repeat 完成后保存轻量 `replayData` 并清理 checkpoint run 目录。parfor batch wall time 约 24 min 38 s。

## Ordinary miss aggregate

| 指标 | 值 | 解释 |
|---|---:|---|
| `numSeed` | 100 | 当前小 MC 总 seed 数。 |
| `numAngleMiss` | 10 | raw `MS-MF-CP-U-in-tooth` 中 `angleErr > 0.002 deg` 的样本数。 |
| `numOrdinaryAngleMiss` | 4 | same tooth、fd healthy、non-collapse 的 ordinary angle miss。 |
| `ordinaryAngleMissRate` | 0.04 | ordinary miss 占全部 seed 的比例。 |
| `collapseHardCount` | 4 | coherence-collapse hard miss 数。 |
| `midCoherenceMissCount` | 0 | 本轮没有落在中间 coherence 档的 miss。 |
| `wrongToothCount` | 0 | controlled in-tooth 条件下没有 wrong-tooth 污染。 |
| `fdNotHealthyCount` | 5 | fd 不健康负样本数。 |
| `easyHitCount` | 87 | raw 已满足 hit threshold 的样本数。 |
| default ordinary median angle | 0.0030293 deg | ordinary miss raw median。 |
| default ordinary P95 angle | 0.0042986 deg | ordinary miss raw P95。 |
| default ordinary median coherence | 1 | ordinary miss 并不表现为 non-ref coherence collapse。 |
| best implementable rescue rate | 1 | ordinary miss 均可由可实现 candidate 救回。 |
| best implementable hit rate | 1 | best implementable 全部进入 `0.002 deg` hit threshold。 |
| best implementable damage rate | 0 | ordinary miss 上没有 damage。 |
| best implementable median angle | 0.00029545 deg | 接近 truth-DoA oracle。 |
| best implementable median gain | 0.0027666 deg | ordinary miss 的典型改善幅度。 |

这张表说明：剩余 raw angle miss 不是单一机制。`collapseHardCount=4` 对应已有 hard-collapse family-safe 线；`numOrdinaryAngleMiss=4` 是新的 non-collapse ordinary miss 线；`fdNotHealthyCount=5` 不应被 ordinary DoA rescue 混入统计。

## Ordinary candidate family aggregate

| candidate group | ordinary miss | hit rate | rescue rate | damage rate | median angle (deg) | P95 angle (deg) | median gain (deg) | median coherence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `default` | 4 | 0 | 0 | 0 | 0.0030293 | 0.0042986 | 0 | 1 |
| `final-small-polish` | 4 | 0 | 0 | 0 | 0.0030293 | 0.0042986 | 0 | 1 |
| `wide-basin-entry` | 4 | 1 | 1 | 0 | 0.00029545 | 0.00049975 | 0.0027666 | 1 |
| `single-mf-basin-entry` | 4 | 0.25 | 0.25 | 0.75 | 0.0075701 | 0.00847 | -0.0049796 | 0.99304 |
| `best-implementable` | 4 | 1 | 1 | 0 | 0.00029545 | 0.00049975 | 0.0027666 | 1 |
| `best-basin-entry` | 4 | 1 | 1 | 0 | 0.00029545 | 0.00049975 | 0.0027666 | 1 |
| `truth-doa-oracle` | 4 | 1 | 1 | 0 | 0.00028681 | 0.00047371 | 0.0027814 | 1 |

关键观察：

- `final-small-polish` 与 `default` 完全一致，说明 ordinary miss 不是 final 点周围的小盒 refine 问题。
- `wide-basin-entry` 全部救回 ordinary miss，且 median / P95 与 truth-DoA oracle 非常接近。
- `single-MF` 对 ordinary miss 不稳，damage rate 达到 `0.75`，不能作为 ordinary miss 默认修复方向。
- `best-implementable` 与 `wide-basin-entry` 一致，说明本轮 ordinary miss 的可实现救力主要来自 wide center。

## Baseline compare：SS-SF 与 SS-MF

当前脚本已经包含 `SS-SF static` 与 `SS-MF CP-U in-tooth` baseline。若只看单星 baseline 的单帧 / 多帧收益，本轮结果如下：

| baseline | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | fdRef median abs (Hz) | fdRate median abs (Hz/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 | 138.4 | 3833.5 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 | 0.0095962 | 6.8894 |

该对比非常关键：`SS-MF CP-U in-tooth` 相比 `SS-SF static` 有明显收益，说明多帧连续相位单星 baseline 本身是健康的。它也为论文中“多帧 CP 相对单帧 static 的收益”提供了受控证据。

如果论文图需要更完整 baseline，建议后续在当前脚本或 completeness smoke 中单独补一个 `baselineCompareTable`，显式列出 `SS-SF static / MS-SF static / SS-MF CP-U / MS-MF CP-U / family-safe / truth oracle`，避免每次从 method aggregate 和 ladder 表里手工提取。

## In-tooth method aggregate

| method | family | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | fdRef median abs (Hz) | fdRate median abs (Hz/s) | coherence median |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | static-baseline | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 | 138.4 | 3833.5 | — |
| `ms-sf-static` | static-baseline | 0.99 | 0.0020359 | 0.0015905 | 0.0034644 | 0.0045189 | 0.66 | 0.004673 | 88.803 | 3833.5 | — |
| `ss-mf-cp-u-in-tooth` | single-center | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 | 0.0095962 | 6.8894 | — |
| `ms-mf-cp-u-in-tooth` | default-target | 1.00 | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 | 0.0099094 | 4.3808 | 1 |
| `ms-mf-cp-u-wide-doa-in-tooth` | wide-center | 1.00 | 0.0021782 | 0.00044087 | 0.0037959 | 0.011068 | 0.93 | 0.014519 | 0.0099498 | 4.3288 | 1 |
| `ms-mf-cp-u-truth-doa-oracle` | truth-doa-label | 1.00 | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 | 0.009783 | 4.318 | 1 |

raw `MS-MF-CP-U-in-tooth` 的 median 最好，但 P95 / P99 / max 被少数 tail 拉坏。`wide-doa-in-tooth` 的 hit rate 从 `0.90` 提到 `0.93`，但 RMSE / P99 / max 明显变差，说明 blanket wide 不能作为默认路径。

## In-tooth DoA ladder

| stage | method | role | trigger | selected by bank | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | gap to truth oracle median (deg) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `ss-sf-static` | single-sat single-frame static baseline | — | — | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 | 0.0014795 |
| 2 | `ms-sf-static` | multi-sat single-frame static baseline | — | — | 0.0020359 | 0.0015905 | 0.0034644 | 0.0045189 | 0.66 | 0.004673 | 0.0012853 |
| 3 | `ss-mf-cp-u-in-tooth` | single-sat multi-frame baseline | — | — | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 | 0.0003028 |
| 4 | `ms-mf-cp-u-in-tooth` | raw multi-sat multi-frame target | — | — | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 | 1.6318e-05 |
| 5a | strict DoA safe-adopt | selected gated wide+single-MF DoA-only bank | 0.04 | 0.02 | 0.0012688 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.0054541 | 1.5279e-05 |
| 5b | family-safe-adopt | family-aware wide+single-MF bank | 0.04 | 0.04 | 0.0010127 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 0.0044189 | 1.5279e-05 |
| 6 | truth-DoA oracle | offline upper bound | — | — | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 | 0 |

这张 ladder 同时说明两点：

- hard-collapse 线仍由 family-safe-adopt 改善，overall hit rate 从 `0.90` 到 `0.94`。
- 即使 family-safe 之后，距离 truth oracle 仍有空间；这部分不能靠继续扩大 hard-collapse gate 解释，ordinary miss 需要单独 gate/adoption。

## Gated rescue aggregate

| bank | hard cases | easy negative | fd-negative | trigger rate | hard trigger | hard rescue | easy damage | fd-negative damage | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | hard hit rate | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `disabled` | 5 | 90 | 5 | 0 | 0 | 0 | 0 | 0 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.20 | 0.0054541 |
| strict safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.4 | 0 | 0 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.60 | 0.0054541 |
| family-safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.8 | 0 | 0 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 1.00 | 0.0044189 |

本轮 ordinary diagnose 同时复现了 gated rescue 结论：family-safe-adopt 仍是当前 hard-collapse controlled in-tooth 的最优候选。它没有解决全部 ordinary miss，但保持了 hard-collapse 线的稳定性。

## 关键 seed 与分类观察

### seed 343：ordinary angle miss 代表

日志中可见 seed `343` 被标为 `ordinary-angle-miss`：

| 字段 | 值 |
|---|---:|
| default angle error | 0.0020406 deg |
| truth-DoA oracle angle error | 0.00040927 deg |
| oracle gap | 0.0016313 deg |
| default tooth idx | 0 |
| fdRef error | -0.012933 Hz |
| fdRate error | 14.234 Hz/s |
| non-ref coherence floor | 1 |

该样本展示了 ordinary miss 的典型症状：tooth 正确、fd 基本健康、coherence 不 collapse，但 DoA 已略超 hit threshold。它不能被 coherence-collapse gate 捕捉。

### seed 256 / 345：hard-collapse 线仍成立

- seed `256`：default angle error `0.0054541 deg`，coherence floor `0.00067049`，family-safe 选择 `single-mf-coarse-doa-grid`，selected angle error `0.00034933 deg`，non-ref coherence 恢复到 `0.99845`。
- seed `345`：default angle error `0.003759 deg`，coherence floor `0.10043`，family-safe 选择 `wide-coarse-doa-grid`，selected angle error `7.5606e-05 deg`，non-ref coherence 恢复到 `1`。

这些结果继续支持 hard-collapse 使用 family-aware `wide + single-MF` basin-entry；但它们不是 ordinary miss 的 gate 证据。

### seed 346：fd-not-healthy 负样本

seed `346` 被标为 `fd-not-healthy-negative`，default angle error `0.0022914 deg`，fdRef error `7.5262 Hz`，non-ref coherence floor `0.98368`，truth-DoA oracle 有上限收益。它不应被 ordinary DoA gate 混入正样本，否则容易把 fd 问题误判成 DoA gate 问题。

## 可观察现象

1. `numOrdinaryAngleMiss=4`，ordinary miss 是真实存在的，但不是当前全部 miss。
2. ordinary miss 的 default non-ref coherence median 为 `1`，所以当前 `coherence-v1` hard-collapse gate 对它天然无效。
3. final-centered small polish 对 ordinary miss 完全无效，说明该类 miss 不是 final 点附近的小扰动问题。
4. `wide-basin-entry` 能把 4 个 ordinary miss 全部救回，且几乎达到 truth-DoA oracle 水平。
5. `single-MF` 对 ordinary miss 明显不稳定，不能和 hard-collapse 的 `wide + single-MF` 结论简单合并。
6. blanket wide 虽然改善 ordinary miss，但在全体 100 seed 上 RMSE / P99 / max 变差，因此不能作为默认路径。
7. family-safe hard-collapse rescue 仍保持有效：hard rescue rate `0.8`、hard angle hit rate `1.0`、easy / fd-negative damage 为 0。

## 当前结论

本轮结果把 controlled in-tooth DoA 问题拆得更清楚：

- hard-collapse miss：继续由 family-safe `wide + single-MF` basin-entry 处理，当前 replay 结果仍支持该候选。
- ordinary non-collapse miss：不是 small polish 问题，主要需要 `wide` center basin-entry；但当前没有 truth-free gate。
- fd-not-healthy miss：不应混入 ordinary DoA rescue 评价，应继续作为负样本或单独机制处理。

因此，当前不能继续通过以下方式冲 hit rate：

- 直接把 hard-collapse coherence threshold 调高；
- 直接扩大 `familySafeAdoptMaxAbsDoaStepDeg`；
- blanket 打开 wide center；
- 把 `single-MF` 当成 ordinary miss 默认 candidate；
- 把 final-small-polish 接入默认路径。

最值得继续推进的是：为 ordinary miss 设计 truth-free `wide` candidate gate / adoption rule。

## 对主流程的影响

- `replayMfInToothOrdinaryAngleMissDiagnose` 结果支持新增一个更窄的后续诊断入口，例如 `replayMfInToothOrdinaryWideGateDiagnose`。
- 后续 gate 不应只用 non-ref coherence；ordinary miss 的 coherence floor 往往为 `1`。
- 可优先检查以下 no-truth proxy：
  - raw final 与 wide-center candidate 的 objective gap；
  - raw final 与 wide-center candidate 的 coherence 是否都接近 1；
  - raw final 与 static / SS-MF / wide center 的 DoA disagreement；
  - wide candidate 的 DoA step 分布；
  - objective flatness / Hessian proxy；
  - fdRef / fdRate 是否保持健康。
- 当前仍不建议下沉 estimator 主核；所有 ordinary-wide gate 先留在 replay 层验证。
- 若要写论文或结果小节，当前可固定的表述是：`SS-MF CP-U` 相比 `SS-SF static` 有明显多帧收益；raw `MS-MF CP-U` 的 median 很好但 tail 被少数 same-tooth miss 拉坏；family-safe 能救 hard-collapse，ordinary miss 则需要另一个 truth-free wide gate。

## 当前建议下一步

1. 保留该 snapshot 作为 ordinary miss 第一版代表性结果，不需要立即扩大 repeat。
2. 新增或扩展 ordinary-wide gate 诊断，只比较 `default` 与 `wide`，不要再默认带入 `single-MF`。
3. 对 ordinary miss 与 easy-hit 样本同时输出 objective / coherence / DoA disagreement proxy，寻找 no-truth adoption 条件。
4. 如果 ordinary-wide gate 能救 4/4 ordinary miss 且 easy damage 为 0，再考虑把 hard-collapse gate 与 ordinary-wide gate 写成两档 replay-level policy。
5. 只有 controlled in-tooth gate 稳定后，才推进 flow-like 或 `runSimpleDynamicSubsetPeriodicFlow` 层验证。
