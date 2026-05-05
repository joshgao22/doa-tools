# replayMfInToothGatedRescueEffectiveness

## 对应 replay

- `test/dev/replay/replayMfInToothGatedRescueEffectiveness.m`

## 观察目标

在 truth-centered half-tooth `fdRef` 范围内，批量验证 no-truth gated rescue 是否能够稳定救回 same-tooth collapse，并确认 easy / fd-not-healthy 负样本不会被 basin-entry bank 误伤。

该 replay 只验证 tooth 已正确时的 in-tooth rescue 上限，不做 subset tooth selection，不改变默认 flow，也不进入 regression。`caseRole / isHardRescued / isDamaged` 可使用 truth 做离线评价；`rescueTriggered / triggerReason / selectedCandidateFamily` 只能来自默认估计结果、non-ref coherence、phase residual、candidate objective 和卫星几何回代诊断。

## 最新代码口径

- gate 口径为 `coherence-v1`：只由 default non-ref coherence collapse 触发；phase residual 字段仅保留为诊断占位。
- 默认旁路已验证不稳或无增益的 `wide-only / single-mf-only / unsafe gated / blanket / joint` 对照。
- 默认保留三类 bank / row：`disabled`、严格 DoA-only safe-adopt、`gated-wide-single-bank-family-safe-adopt`。
- strict safe-adopt 使用统一 `safeAdoptMaxAbsDoaStepDeg=0.004`。
- family-safe-adopt 只对 `wide / single-MF` basin-entry family 放宽 DoA step guard 到 `0.006 deg`；default / final-centered candidate 仍按严格 guard 处理。
- `angleHitThresholdDeg=0.002`，method / ladder / rescue aggregate 输出 DoA hit rate 与 P99。
- `includeJointSafeAdoptBank=false`；joint bank 只在 envelope 结果支持后显式打开。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` | `baseSeed=253`，`numRepeat=100`，`snrDb=10`，`coherenceThreshold=0.20`，`safeStep=0.004 deg`，`familyStep=0.006 deg`，`safeCoh=0.95`，`safeFdRefStep=300 Hz`，`joint=0`，`family=1`，`angleHitThreshold=0.002 deg` | family-safe-adopt 是当前最好的 controlled in-tooth rescue 版本：RMSE / median / P95 / P99 / max / hit rate 全部优于 raw MS-MF 与 strict safe-adopt；hard rescue rate 从 strict 的 `0.4` 提升到 `0.8`，easy / fd-negative damage 仍为 0。 |

## 当前代表性结果

2026-05-04 的 100-repeat family-safe-adopt run 是当前唯一保留的代表性结果。snapshot 保存为 `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat`。

运行配置摘要：

| 配置项 | 值 |
|---|---:|
| repeats | 100 |
| SNR | 10 dB |
| base seed | 253 |
| fd oracle half-tooth fraction | 0.490 |
| fdRate oracle half-width | 1000 Hz/s |
| gated rescue coherence threshold | 0.200 |
| safe adopt max DoA step | 0.004 deg |
| family safe max DoA step | 0.006 deg |
| safe adopt coherence threshold | 0.950 |
| safe adopt max fdRef step | 300 Hz |
| angle hit threshold | 0.002 deg |

checkpoint 从空目录开始，完成后保存轻量 `replayData` 并清理 checkpoint run 目录。

## In-tooth method aggregate

| method | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 |
| `ms-sf-static` | 0.99 | 0.0020359 | 0.0015905 | 0.0034644 | 0.0045189 | 0.66 | 0.004673 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 |
| `ms-mf-cp-u-in-tooth` | 1.00 | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 |
| `ms-mf-cp-u-wide-doa-in-tooth` | 1.00 | 0.0021782 | 0.00044087 | 0.0037959 | 0.011068 | 0.93 | 0.014519 |
| `ms-mf-cp-u-truth-doa-oracle` | 1.00 | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 |

raw MS-MF 的 median 已经很好，但 P95 / P99 / max 被少数 hard seed 拉坏；truth-DoA oracle 说明这些 tail 仍存在可救上限。

## In-tooth DoA ladder

| stage | method | trigger | selected by bank | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | `ss-mf-cp-u-in-tooth` | — | — | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 |
| 4 | `ms-mf-cp-u-in-tooth` | — | — | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 |
| 5a | strict DoA safe-adopt | 0.04 | 0.02 | 0.0012688 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.0054541 |
| 5b | family-safe-adopt | 0.04 | 0.04 | 0.0010127 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 0.0044189 |
| 6 | truth-DoA oracle | — | — | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 |

family-safe-adopt 相比 raw MS-MF 同时改善 RMSE、median、P95、P99、max 与 hit rate；相比 strict safe-adopt，也明显压低 P99 / max，并把 selected-by-bank rate 从 `0.02` 提高到 `0.04`。

## Gated rescue aggregate

| bank | hard cases | easy negative | fd-negative | trigger rate | hard trigger | hard rescue | easy damage | fd-negative damage | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | hard hit rate | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `disabled` | 5 | 90 | 5 | 0 | 0 | 0 | 0 | 0 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.20 | 0.0054541 |
| strict safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.4 | 0 | 0 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.60 | 0.0054541 |
| family-safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.8 | 0 | 0 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 1.00 | 0.0044189 |

family-safe-adopt 的关键价值是：只触发 4% case，却把 hard rescue rate 从 strict 的 `0.4` 提到 `0.8`，并让 hard angle hit rate 达到 `1.0`；同时 easy / fd-negative damage 仍为 0。

## Gated rescue verdict

| bank | hardRescueRate | easyDamageRate | fdNegativeDamageRate | P95 improvement (deg) | max improvement (deg) | recommendedPass | recommendation |
|---|---:|---:|---:|---:|---:|:---:|---|
| strict safe-adopt | 0.4 | 0 | 0 | 0.0013185 | 0 | false | `hold` |
| family-safe-adopt | 0.8 | 0 | 0 | 0.0015752 | 0.0010353 | true | `candidate-but-not-target` |

`candidate-but-not-target` 的含义是：该 replay 已经支持 family-safe-adopt 作为 controlled in-tooth rescue candidate，但它还不是 full-flow / default estimator 结论。

## 关键 seed 解释

### seed 256：family-safe-adopt 解决 strict guard 过紧问题

seed `256` 是这轮最关键的 hard-collapse 样本：

| 字段 | 值 |
|---|---:|
| default angle error | 0.0054541 deg |
| truth-DoA oracle angle error | 3.2556e-05 deg |
| default non-ref coherence floor | 0.00067049 |
| selected family | `single-mf-coarse-doa-grid` |
| selected tag | `coarse-center-dlat+0-dlon+0.004` |
| selected angle error | 0.00034933 deg |
| selected non-ref coherence floor | 0.99845 |
| angle gain | 0.0051048 deg |
| safe adopt reason | `accepted` |

strict safe-adopt 会因为 `doa-step-too-large` 拒绝这个 candidate；family-safe-adopt 允许 `wide / single-MF` center 使用 `0.006 deg` guard，因此成功采纳该 candidate。

### seed 345：wide center 继续稳定救回

seed `345` 的 default angle error 为 `0.003759 deg`，default non-ref coherence floor 为 `0.10043`。family-safe-adopt 选择 `wide-coarse-doa-grid` 的中心候选，selected angle error 降到 `7.5606e-05 deg`，non-ref coherence 恢复到 `1`，angle gain 为 `0.0036834 deg`。

### seed 253：唯一 hard miss，但不影响 hit threshold

当前 hard miss summary 只剩 seed `253`：

| 字段 | 值 |
|---|---:|
| default angle error | 0.001526 deg |
| truth-DoA oracle angle error | 0.00048515 deg |
| default non-ref coherence floor | 0.20455 |
| rescueTriggered | false |
| miss reason | `gate-not-triggered-bank-not-confirmed` |

该 seed 没有触发，是因为 coherence floor `0.20455` 略高于 gate threshold `0.2`。不过它的 default angle error 已小于 `0.002 deg` hit threshold，所以不影响 hard angle hit rate。

## 可观察现象

- `phaseResidAvailable=false`，当前实际 gate 仍只依赖 non-ref coherence collapse。
- `triggerRate=0.04`、`hardTriggerRate=0.8`、`easyTriggerRate=0`、`fdNegativeTriggerRate=0`，说明 gate 仍很干净。
- family-safe-adopt 并没有 blanket 改所有样本，而是只多采纳了 strict guard 原本拒绝的可实现 wide / single-MF candidate。
- `MS-MF-CP-U raw in-tooth` 的 median 已优于 `SS-MF-CP-U-in-tooth`，但 tail 更差；family-safe-adopt 把 RMSE 拉到接近 SS-MF，并显著降低 P99 / max。
- overall hit rate 从 `0.90` 提升到 `0.94`，尚未达到 `0.99`；剩余 miss 不全是 coherence-collapse hard case，后续若要冲 99% 需要单独分析 non-collapse 普通 >0.002 deg miss。

## 当前结论

`gated-wide-single-bank-family-safe-adopt` 是当前 controlled in-tooth DoA rescue 的最优候选。它在不改变 estimator / flow / objective 的前提下，显著压低 raw MS-MF 的 DoA tail，并且保持 easy / fd-negative damage 为 0。

当前结论仍然是 replay-level：

- 可以固定为 controlled in-tooth rescue 结果；
- 暂不下沉到 estimator 主核；
- 暂不接入 full-flow 默认路径；
- 暂不继续默认打开 joint fdRef / fdRate bank；
- 若后续追求 `>99%` hit rate，应新开普通 non-collapse miss 诊断，而不是继续放宽当前 hard-collapse rescue。

## 对主流程的影响

- 结果支持把 `family-safe-adopt` 作为下一阶段 flow-like gated same-tooth basin-entry 的候选，但接入前仍需在真实 subset-periodic flow 中验证 gate、candidate adoption 和 wrong-tooth interaction。
- 不建议继续扩大 repeat 只为确认当前结论；100-repeat 已足够说明 family-safe-adopt 的方向。
- 不建议把 `coherenceThreshold` 直接调高或 blanket 放宽 DoA；seed `253` 虽是 gate miss，但已满足 `0.002 deg` hit threshold。
- 不建议默认打开 joint bank；对应 envelope 结果没有支持非零 fdRef joint step 的稳定收益。
