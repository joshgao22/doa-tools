# replayMfInToothGatedRescueEffectiveness

## 对应 replay

- `test/dev/replay/replayMfInToothGatedRescueEffectiveness.m`

## 观察目标

在 truth-centered half-tooth `fdRef` 范围内，批量验证 no-truth gated rescue 是否能够稳定救回 same-tooth collapse，并确认 easy / fd-not-healthy 负样本不会被 basin-entry bank 误伤。

该 replay 只验证 tooth 已正确时的 in-tooth rescue 上限，不做 subset tooth selection，不改变默认 flow，也不进入 regression。`caseRole / isHardRescued / isDamaged` 可使用 truth 做离线评价；`rescueTriggered / triggerReason / selectedCandidateFamily` 只能来自默认估计结果、non-ref coherence、phase residual、candidate objective 和卫星几何回代诊断。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260427-132413.mat` | `baseSeed=253, numRepeat=100, snrDb=10`，`coherence-or-phase-v2` gate，`checkpointEnable=true` | `gated-wide-single-bank` 是唯一通过 verdict 的 bank；hard rescue rate 为 `0.85714`，easy / fd-not-healthy 负样本无误伤，overall P95 / max 均优于 disabled；blanket reference 虽全救 hard case，但会误伤 easy / fd-negative 并拉坏 max。 |

## 当前代表性结果

2026-04-27 的 100-repeat v2 gate confirmation 是当前唯一保留的代表性结果。配置为 `baseSeed=253`、`numRepeat=100`、`snrDb=10`、`checkpointEnable=true`，snapshot 保存为 `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260427-132413.mat`。

运行耗时约 `23 min 22 s`。checkpoint 从空目录开始，完成后保存轻量 `replayData` 并清理 checkpoint run 目录。

### In-tooth estimator 上限

`MS-MF-CP-U-in-tooth` 在 half-tooth oracle 范围内保持 `toothHitRate=1`，`angleMedianDeg=0.00057009`，`angleP95Deg=0.0036204`，`fdRefAbsMedianHz=0.009968`，`fdRateAbsMedianHzPerSec=4.3077`。这说明 tooth 与频率链本身健康，剩余坏点主要是 same-tooth DoA / local-state basin。

`MS-MF-CP-U-truth-doa-oracle` 的 `angleP95Deg=0.00096105`，明显优于 default target，说明少数 tail seed 可以被更好的 DoA basin 救回。

### Gated rescue aggregate

| bank | triggerRate | hardTriggerRate | easyTriggerRate | fdNegativeTriggerRate | hardRescueRate | easyDamageRate | fdNegativeDamageRate | overallP95AngleDeg | overallMaxAngleDeg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `disabled` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0036204 | 0.0054539 |
| `gated-wide-only` | 0.06 | 0.85714 | 0 | 0 | 0.57143 | 0 | 0 | 0.0023635 | 0.0044191 |
| `gated-single-mf-only` | 0.06 | 0.85714 | 0 | 0 | 0.42857 | 0 | 0 | 0.0024993 | 0.0080813 |
| `gated-wide-single-bank` | 0.06 | 0.85714 | 0 | 0 | 0.85714 | 0 | 0 | 0.0020833 | 0.0044191 |
| `blanket-wide-single-bank` | 1 | 1 | 1 | 1 | 1 | 0.011364 | 0.2 | 0.00099125 | 0.010301 |

本轮共有 7 个 hard-collapse、88 个 easy-negative、5 个 fd-not-healthy-negative。`gated-wide-single-bank` 只触发 6 个 hard-collapse，不触发 easy / fd-negative，因此当前 v2 gate 的 no-truth 触发行为基本符合预期。相比 disabled，`gated-wide-single-bank` 将 overall P95 从 `0.0036204 deg` 降到 `0.0020833 deg`，并将 max 从 `0.0054539 deg` 降到 `0.0044191 deg`。

### Gated rescue verdict

| bank | hardRescueRate | easyDamageRate | fdNegativeDamageRate | p95ImprovementVsDisabledDeg | maxImprovementVsDisabledDeg | recommendedPass | recommendation |
|---|---:|---:|---:|---:|---:|:---:|---|
| `gated-wide-only` | 0.57143 | 0 | 0 | 0.0012569 | 0.0010348 | false | `hold` |
| `gated-single-mf-only` | 0.42857 | 0 | 0 | 0.0011211 | -0.0026274 | false | `hold` |
| `gated-wide-single-bank` | 0.85714 | 0 | 0 | 0.0015371 | 0.0010348 | true | `ready-for-flow-like-replay` |
| `blanket-wide-single-bank` | 1 | 0.011364 | 0.2 | 0.0026292 | -0.0048466 | false | `hold` |

`gated-wide-single-bank` 是唯一通过 hard rescue、damage、P95 不变差、max 不变差和 hard-case 数量检查的 bank。单独 `wide-only` 与 `single-mf-only` 都能提供一定改善，但 hard rescue rate 不足；其中 `single-mf-only` 还会拉坏 max，因此不应作为独立主路线。

### 关键 seed 解释

seed `256` 是代表性 hard-collapse：default angle error 为 `0.0054539 deg`，default non-ref coherence floor 为 `0.00066016`，truth-DoA oracle angle error 为 `2.1464e-05 deg`。v2 gate 由 `coherence-collapse-phase-unavailable` 触发，`gated-wide-single-bank` 选择 `single-mf-coarse-doa-grid` 中 `coarse-center-dlat+0-dlon+0.004`，selected angle error 降到 `0.00034936 deg`，non-ref coherence 恢复到 `0.99845`。

seed `345` 也是 hard-collapse：default angle error 为 `0.0037598 deg`，default non-ref coherence floor 为 `0.10043`，truth-DoA oracle angle error 为 `6.1079e-05 deg`。v2 gate 同样由 `coherence-collapse-phase-unavailable` 触发，`gated-wide-single-bank` 选择 `wide-coarse-doa-grid` 的中心候选，selected angle error 降到 `7.5596e-05 deg`，non-ref coherence 恢复到 `1`。

### blanket reference 的作用

`blanket-wide-single-bank` 只作为误伤参考保留。它说明 wide + single-MF candidate bank 本身有能力全救 hard case，并能把 overall P95 降到 `0.00099125 deg`；但它会误伤 easy 与 fd-not-healthy-negative 样本，`easyDamageRate=0.011364`、`fdNegativeDamageRate=0.2`，并把 `overallMaxAngleDeg` 拉到 `0.010301`。因此 blanket 不能进入默认流程，只能保留为 replay reference。

## 主要统计口径

- `rescueEffectAggregateTable`：按 bank 汇总 trigger rate、hard rescue rate、easy / fd-negative damage rate、overall angle 分布。
- `rescueEffectVerdictTable`：按固定 replay promotion checks 汇总各 bank 是否满足 hard rescue、damage、P95 / max 不变差等条件。
- `rescueEffectCaseTable`：每个 seed 保留 `gated-wide-single-bank` 的 default 指标、trigger、selected candidate、angle gain、hard rescue 与 damage 标志。
- `rescueBankDecisionTable`：逐 seed、逐 bank 保存 `disabled`、gated bank 与 blanket reference 的完整选择记录。
- `triggerReasonTable`：只保留 no-truth gate 诊断，包括 non-ref coherence floor、phase residual、`phaseResidAvailable` 和 trigger reason。

## 可观察现象

- `phaseResidAvailable` 在本轮均为 false，因此实际 gate 主要依赖 non-ref coherence collapse；phase residual 缺失不会 veto coherence trigger。
- `gated-wide-single-bank` 的 `easyTriggerRate=0`、`fdNegativeTriggerRate=0`，说明当前 no-truth gate 没有覆盖 easy / fd-not-healthy 负样本。
- 已触发的 hard case 基本都被救回；剩余未救 hard case主要是 gate 漏检，而不是 candidate bank 完全无效。
- blanket reference 的 P95 更好但 max 更差，说明它不能作为默认常驻 rescue。

## 当前结论

当前代码不需要改变 estimator、flow 或 rescue candidate objective。100-repeat confirmation 支持 `gated-wide-single-bank` 作为 in-tooth conditional rescue 的首选候选：它是唯一通过 verdict 的 bank，能显著降低 P95 / max，且 easy / fd-not-healthy 负样本无误伤。

该结论仍然是 in-tooth oracle 条件下的 replay 证据，不等价于已经可以直接接入 `runSimpleDynamicSubsetPeriodicFlow`。下一步应写 flow-like replay，验证 subset tooth selection、periodic same-tooth refine 与 gated rescue 是否能在真实 flow 条件下衔接。

## 对主流程的影响

- 不再继续为该 replay 增加重复次数；当前 100-repeat confirmation 已足够作为 in-tooth gated rescue 证据。
- 下一步进入 flow-like replay，而不是直接改默认 flow。
- 若 flow-like replay 中 damage 变高，再转入 gate threshold scan；threshold scan 应放到 scan，不要继续膨胀该 replay。
- 暂不把 `gated-wide-single-bank` 抽成正式 estimator helper；当前仍保留在 replay / flow-like 验证层，避免未经过真实 flow 验证的机制过早固化。
