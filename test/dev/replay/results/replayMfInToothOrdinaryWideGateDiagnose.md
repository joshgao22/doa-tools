# replayMfInToothOrdinaryWideGateDiagnose 结果记录

## 对应 replay

`test/dev/replay/replayMfInToothOrdinaryWideGateDiagnose.m`

## 观察目标

在 controlled in-tooth 条件下，检查 ordinary non-collapse DoA miss 是否能由 truth-free 的 `wide` basin-entry gate 接住，同时避免 easy / fd-negative damage。

该 replay 不做 subset tooth selection，不改变 estimator / flow 默认路径。truth 只用于 offline label、rescue / damage 评价和 truth-DoA oracle 上限。

本轮 500-repeat 结果用于覆盖此前 100-repeat 的乐观结论：ordinary-wide gate 仍有诊断价值，但不足以作为稳定 replay-level preferred candidate。

## Snapshot index

| snapshot | 日期 | 状态 | 备注 |
|---|---:|---|---|
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-233931.mat` | 2026-05-05 | representative | 500-repeat 扩大验证；template cleanup 后默认旁路 final-small-polish、single-MF ordinary probe 与 hard-collapse auxiliary tables。结果将 ordinary-wide 从 preferred candidate 降级为 hold / partial diagnostic。 |
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-202052.mat` | 2026-05-05 | superseded | 100-repeat 中 `wide-obj10-minStep0.001-maxStep0.004` 曾表现为 preferred candidate；该结论未通过 500-repeat 扩大验证。 |
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-181942.mat` | 2026-05-05 | superseded | 首轮 wide-gate sweep；证明 `obj0` 可救 ordinary 但 easy trigger 过宽。 |

## 当前代表性配置

- `numRepeat = 500`
- `snrDb = 10`
- `baseSeed = 253`
- seed range: `253:752`
- `fd oracle half-tooth fraction = 0.49`
- `fdRate oracle half-width = 1000 Hz/s`
- `ordinaryCoherenceMinThreshold = 0.95`
- `angleHitThresholdDeg = 0.002`
- `includeFinalSmallPolishProbe = false`
- `includeSingleMfBasinEntryProbe = false`
- `includeHardCollapseAuxTables = false`
- `includeJointSafeAdoptBank = false`
- `includeFamilySafeAdoptBank = false`
- `wideGateObjGainThresholdList = [0 5 10 20 30 50 80 100 150 200]`
- `wideGateMinAbsDoaStepDegList = [0 0.001 0.0015 0.002]`
- `wideGateMaxAbsDoaStepDegList = [0.004 0.006]`
- `wideGateCandidateCoherenceThreshold = 0.95`
- `wideGateRequireDefaultNonCollapse = true`
- preferred policy preset: `objGain >= 10` and `0.001 <= wide DoA step <= 0.004 deg`

运行约 1 h 44 min，保存轻量 `replayData` snapshot，并成功清理 checkpoint run 目录。

## Ordinary miss 结构

| 指标 | 数值 |
|---|---:|
| `numSeed` | 500 |
| `numAngleMiss` | 65 |
| `numOrdinaryAngleMiss` | 22 |
| `ordinaryAngleMissRate` | 0.044 |
| `collapseHardCount` | 39 |
| `midCoherenceMissCount` | 0 |
| `wrongToothCount` | 0 |
| `fdNotHealthyCount` | 13 |
| `easyHitCount` | 426 |
| `defaultOrdinaryMedianAngleDeg` | 0.0022747 |
| `defaultOrdinaryP95AngleDeg` | 0.0036078 |
| `defaultOrdinaryMedianCoherence` | 1 |
| `bestImplementableRescueRate` | 0.59091 |
| `bestImplementableHitRate` | 0.59091 |
| `bestImplementableDamageRate` | 0 |
| `bestImplementableMedianAngleDeg` | 0.00063892 |
| `bestImplementableMedianGainDeg` | 0.0015971 |
| `angleHitThresholdDeg` | 0.002 |

观察：500-repeat 下 ordinary miss 数从 100-repeat 的 4 个扩展到 22 个，样本量更可信。它们的 median coherence 为 1，说明 ordinary miss 不表现为 non-ref coherence collapse，不能用 hard-collapse gate 直接捕获。

## Ordinary candidate family aggregate

| candidate group | ordinary miss | hit rate | rescue rate | damage rate | median angle err (deg) | P95 angle err (deg) | median gain (deg) | median coherence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `default` | 22 | 0 | 0 | 0 | 0.0022747 | 0.0036078 | 0 | 1 |
| `wide-basin-entry` | 22 | 0.63636 | 0.63636 | 0.31818 | 0.00046751 | 0.0070408 | 0.00164 | 1 |
| `best-implementable` | 22 | 0.59091 | 0.59091 | 0 | 0.00063892 | 0.0027479 | 0.0015971 | 1 |
| `best-basin-entry` | 22 | 0.63636 | 0.63636 | 0.31818 | 0.00046751 | 0.0070408 | 0.00164 | 1 |
| `truth-doa-oracle` | 22 | 1 | 1 | 0 | 0.00040962 | 0.00076736 | 0.0018871 | 1 |

关键观察：

1. `wide-basin-entry` 仍有明显救力，median angle error 接近 oracle，但并不稳定。
2. blanket wide 在 ordinary miss 上 damage rate 达到 `0.31818`，不能作为默认路径。
3. `best-implementable` 通过 guard 避开 damage，但 rescue / hit rate 只剩 `0.59091`，说明当前 truth-free guard 只能形成 partial rescue。
4. truth-DoA oracle 对 22 个 ordinary miss 全部有效，说明这些 tail 仍属于可解释的 basin / adoption 问题，而不是模型本身无解。

## Ordinary-wide recommended policy

当前自动推荐仍选中 preset policy：

```text
wide-obj10-minStep0.001-maxStep0.004
```

但 500-repeat 下推荐结论为 `hold`，不是 `preferred-candidate`。

| policy | ordinary rescue | ordinary hit | easy trigger | easy trigger count | easy damage | fd-negative damage | overall hit | P95 (deg) | P99 (deg) | max (deg) | recommendation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `wide-obj10-minStep0.001-maxStep0.004` | 0.59091 | 0.59091 | 0.023474 | 10 | 0 | 0 | 0.898 | 0.0037382 | 0.0051413 | 0.0054663 | `hold` |

与 disabled 对比：

| policy | overall hit | median (deg) | P95 (deg) | P99 (deg) | max (deg) |
|---|---:|---:|---:|---:|---:|
| `disabled` | 0.87 | 0.0004843 | 0.0037611 | 0.0051413 | 0.0054663 |
| `wide-obj10-minStep0.001-maxStep0.004` | 0.898 | 0.00043519 | 0.0037382 | 0.0051413 | 0.0054663 |

该 policy 有小幅收益，并保持 zero easy / fd-negative damage，但 ordinary rescue coverage 不足，不能推进到 combined / flow 默认路径。

## Gate policy sweep 关键观察

### 过宽 policy

| policy | ordinary rescue | easy trigger | easy trigger count | easy damage | overall hit | recommendation |
|---|---:|---:|---:|---:|---:|---|
| `wide-obj0-minStep0-maxStep0.004` | 0.59091 | 0.53521 | 228 | 0.0023474 | 0.898 | `unsafe` |
| `wide-obj0-minStep0-maxStep0.006` | 0.59091 | 0.53521 | 228 | 0.0023474 | 0.898 | `unsafe` |

`obj0/minStep0` 在 500-repeat 中不仅 easy trigger 过宽，还出现 easy damage，因此 blanket wide / near-blanket wide 必须明确否定。

### 保守但 coverage 不足的 policy

| policy | ordinary rescue | easy trigger | easy trigger count | easy damage | overall hit | recommendation |
|---|---:|---:|---:|---:|---:|---|
| `wide-obj0-minStep0.001-maxStep0.004` | 0.59091 | 0.042254 | 18 | 0 | 0.898 | `hold` |
| `wide-obj5-minStep0.001-maxStep0.004` | 0.59091 | 0.037559 | 16 | 0 | 0.898 | `hold` |
| `wide-obj10-minStep0.001-maxStep0.004` | 0.59091 | 0.023474 | 10 | 0 | 0.898 | `hold` |

这些 policy 能控制误伤，但都只救约 13/22 ordinary miss。

### 过紧 policy

| policy family | 典型 ordinary rescue | 现象 |
|---|---:|---|
| `minStep >= 0.0015` | 0.45455 | 开始漏掉大量 ordinary miss。 |
| `minStep >= 0.002` | 0.27273 | coverage 明显不足。 |
| `objGain >= 20` | 0.40909 | objective gate 开始过紧。 |
| `objGain >= 50` | 0.045455 | 几乎丢掉 ordinary miss。 |
| `objGain >= 150` | 0 | ordinary miss 完全不触发。 |

结论：单纯继续调 `objGain` 与 DoA-step lower guard 已经到头。阈值过松会误伤 / 触发过宽，阈值过紧会漏救；当前普通二维 gate 不足以稳定解决 ordinary miss。

## Wide candidate feature table

| miss type | num case | default hit | wide hit | wide rescue | wide damage | median obj gain | P95 obj gain | median angle gain (deg) | median wide angle (deg) | P95 wide angle (deg) | median default coherence | median wide coherence | median max abs DoA step (deg) | P05 step | P95 step |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `easy-hit` | 426 | 1 | 0.93427 | 0 | 0.077465 | 5.1487e-05 | 7.1051 | 1.5644e-06 | 0.00037938 | 0.0039108 | 1 | 1 | 1.5123e-05 | 1.2501e-06 | 0.0033126 |
| `collapse-hard` | 39 | 0 | 0.64103 | 0 | 0.12821 | 8.836e+05 | 9.0654e+05 | 0.0030895 | 0.00076606 | 0.010781 | 0.10011 | 1 | 0.004049 | 0.003045 | 0.0079204 |
| `fd-not-healthy-negative` | 13 | 0.69231 | 0.69231 | 0 | 0.23077 | 5.4092e+05 | 8.2655e+05 | 0.00093587 | 0.00032224 | 0.0090909 | 0.98268 | 1 | 0.002059 | 0.00071691 | 0.0081767 |
| `ordinary-angle-miss` | 22 | 0 | 0.63636 | 0.63636 | 0.31818 | 12.93 | 46.942 | 0.00164 | 0.00046751 | 0.0070408 | 1 | 1 | 0.0020616 | 0.0011673 | 0.0075163 |

该表是本轮最重要的新证据：ordinary miss 并非单一机制，wide 对其中一部分有效，但也会 damage。collapse-hard 的 objective gain 与 coherence pattern 与 ordinary 完全不同，应继续分支处理。

## In-tooth method aggregate

| method | family | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | fdRef median abs (Hz) | fdRate median abs (Hz/s) | coherence median |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | static-baseline | 0.968 | 0.0022082 | 0.0018705 | 0.0038167 | 0.0047003 | 0.548 | 0.0056668 | 116.96 | 3833.5 | — |
| `ms-sf-static` | static-baseline | 0.996 | 0.0019684 | 0.0016108 | 0.0034213 | 0.0039843 | 0.632 | 0.0047385 | 82.993 | 3833.5 | — |
| `ss-mf-cp-u-in-tooth` | single-center | 1.00 | 0.00095184 | 0.00071409 | 0.0017016 | 0.0028258 | 0.966 | 0.0034199 | 0.0099662 | 6.7737 | — |
| `ms-mf-cp-u-in-tooth` | default-target | 1.00 | 0.001487 | 0.0004843 | 0.0037611 | 0.0051413 | 0.87 | 0.0054663 | 0.010534 | 4.9728 | 1 |
| `ms-mf-cp-u-wide-doa-in-tooth` | wide-center | 1.00 | 0.0027142 | 0.00039533 | 0.005958 | 0.011036 | 0.906 | 0.01476 | 0.0104 | 4.699 | 1 |
| `ms-mf-cp-u-truth-doa-oracle` | truth-doa-label | 1.00 | 0.00047336 | 0.00032835 | 0.00093547 | 0.0011753 | 1.00 | 0.0016946 | 0.010086 | 4.7977 | 1 |

关键解释：

- `MS-MF CP-U raw` 的 median 优于 `SS-MF CP-U`，说明多星多帧局部精度潜力存在。
- `MS-MF CP-U raw` 的 P95/P99/max 明显差，说明当前 low-complexity / in-tooth realization 有 bad-basin / threshold tail。
- `MS-MF CP-U wide` 提升 hit rate 到 `0.906`，但 RMSE / P95 / P99 / max 明显变差，说明 blanket wide 不是可接受默认路径。
- truth-DoA oracle 达到 `angleHitRate=1` 且 tail 很小，说明当前 tail 更像 basin/adoption 问题，而非模型不可达。

## In-tooth DoA ladder

| stage | method | role | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | gap to truth oracle median (deg) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | `ss-sf-static` | single-sat single-frame static baseline | 0.0022082 | 0.0018705 | 0.0038167 | 0.0047003 | 0.548 | 0.0056668 | 0.0014568 |
| 2 | `ms-sf-static` | multi-sat single-frame static baseline | 0.0019684 | 0.0016108 | 0.0034213 | 0.0039843 | 0.632 | 0.0047385 | 0.0012819 |
| 3 | `ss-mf-cp-u-in-tooth` | single-sat multi-frame baseline | 0.00095184 | 0.00071409 | 0.0017016 | 0.0028258 | 0.966 | 0.0034199 | 0.00029354 |
| 4 | `ms-mf-cp-u-in-tooth` | raw multi-sat multi-frame target | 0.001487 | 0.0004843 | 0.0037611 | 0.0051413 | 0.87 | 0.0054663 | 1.7052e-05 |
| 6 | `ms-mf-cp-u-truth-doa-oracle` | offline upper bound | 0.00047336 | 0.00032835 | 0.00093547 | 0.0011753 | 1.00 | 0.0016946 | 0 |

该 ladder 对论文主线有启发：多星多帧 raw target 在中位数上最接近 truth oracle，但 tail 影响 RMSE / P95 / hit rate。后续论文主仿真应回到 CRB consistency、conditional RMSE 与 resolved/outlier rate，而不是把该 replay 当作 99% hit-rate 工程目标。

## 可观察现象

1. `wrongToothCount=0`，说明本轮确实剥离了 wrong-tooth 污染。
2. `numOrdinaryAngleMiss=22`，ordinary non-collapse miss 在 500-repeat 下稳定存在。
3. ordinary miss 的 median angle 仅为 `0.0022747 deg`，很多样本只是略超 `0.002 deg` 诊断阈值；不应把所有 ordinary miss 都解释成必须强修的 hard failure。
4. `wide-basin-entry` 对 ordinary miss 有救力，但 damage rate 达到 `0.31818`；ordinary miss 应继续拆成 wide-rescued / wide-damaged / wide-not-enough 子类型。
5. `objGain >= 10` 与 `0.001 <= wide DoA step <= 0.004 deg` 能把 easy damage 压到 0，但 ordinary rescue coverage 只有 `0.59091`。
6. blanket wide 在 easy-hit 上也有 `0.077465` damage rate，因此不能作为默认路径。
7. collapse-hard 与 ordinary miss 的 feature 分布明显不同，二者不能混用同一个 gate。
8. 该 replay 当前更适合作为 outlier / bad-basin diagnostic，而不是主论文精度目标。

## 当前结论

1. 100-repeat 中的 `preferred-candidate` 结论被 500-repeat 扩大验证推翻。
2. 当前 ordinary-wide gate 只能定性为 `hold / partial diagnostic candidate`。
3. `wide-obj10-minStep0.001-maxStep0.004` 是低误伤 partial rescue：它改善 overall hit rate，但只救约 59% ordinary miss，不能晋级 flow 或 estimator 默认路径。
4. blanket wide 与 `obj0/minStep0` 系列明确不安全。
5. 继续只调 `objGain` 与 DoA-step threshold 不太可能解决全部 ordinary miss；下一步若继续此线，应先做 outcome classification 与 baseline-relative analysis。
6. 对论文主线而言，该 replay 的最重要结论不是“如何把 hit rate 冲到 99%”，而是说明当前 low-complexity realization 存在 threshold / bad-basin outlier；主性能验证应回到 CRB consistency、conditional RMSE、resolved rate 与 SS/MS/SF/MF/CP 对比。

## 对主流程的影响

- 不把 ordinary-wide gate 推进到 flow-like / flow 默认路径。
- 不把该结果写成 regression 契约。
- 不继续围绕 ordinary-wide 二维阈值调参来追 99% hit rate。
- hard-collapse 的 family-safe `wide + single-MF` 分支仍可作为独立机制保留，但不能用本 replay 的 ordinary 结果为其背书。
- 论文-facing 仿真应优先报告：
  - RMSE / CRB 或 EFIM ratio；
  - conditional RMSE on resolved / local-converged samples；
  - resolved / outlier rate；
  - SS-SF、MS-SF、SS-MF、MS-MF、CP/IP、known/unknown-rate 对照。

## 建议下一步

若继续维护该 replay，建议只做诊断增强，不做默认 rescue 接入：

1. 新增 `ordinaryWideOutcomeTable`，把 22 个 ordinary miss 拆成 `wide-rescued / wide-damaged / wide-not-enough / gate-rejected`。
2. 新增 `angleExcessDeg = angleErrDeg - angleHitThresholdDeg`，区分 near-threshold 与真正 tail ordinary miss。
3. 新增 same-seed baseline-relative 字段，比较 `SS-SF / MS-SF / SS-MF / MS-MF / wide / truth oracle`，避免只用固定 `0.002 deg` 阈值解释精度。
4. 做 `fdHealthy` threshold sensitivity，把 `1 Hz / 50 Hz/s` strict 口径和更宽松口径分开。
5. 若论文结果需要，新增 CRB-normalized error 与 resolved-rate 表；这比继续打磨 ordinary-wide rescue 更贴近主线。

在此之前，本结果文档已足够支撑当前阶段判断，不需要继续修改代码才能写入排障记录或 results 索引。
