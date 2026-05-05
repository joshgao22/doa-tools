# replayMfFlowLikeGatedBasinEntryEffectiveness

## 对应 replay

- `test/dev/replay/replayMfFlowLikeGatedBasinEntryEffectiveness.m`

## 观察目标

在真实 `subset-periodic flow` 结构内，比较 disabled baseline 与 no-truth gated `wide + single-MF` same-tooth basin-entry rescue。该入口不使用 truth-centered half-tooth oracle，也不跳过 subset selection，因此会同时暴露：

1. basin-entry candidate family 在真实 flow 中是否仍有救力；
2. same-tooth in-tooth DoA hard case 是否能被 no-truth gate 捕捉；
3. gate / adoption 是否会误伤 easy case 或放大 tail；
4. wrong-tooth / subset selection 因素是否会污染该 rescue 的整体统计。

该 replay 只作为 flow-like stress check 和机制分流入口；不改变 estimator、objective、residual、reference-sat 语义，也不作为默认 flow 通过依据。

## 最新代码口径

- disabled 与 gated 两个 batch 使用同一组 `seedList`。
- `focusedSeedList=[]` 时使用 `baseSeed + (0:(numRepeat-1))` 的连续 seed 小 MC；`focusedSeedList` 非空时只跑指定 seed 做 triage。
- gated 分支通过 `buildSimpleDynamicFlowOpt.sameToothRescue` 打开 `wide + single-MF` same-tooth basin-entry candidate。
- gate 只读取 selected periodic decision summary 与 subset trust 诊断中的 non-ref coherence、phase residual 和 trust 状态；truth 只用于结果表里的 angle / tooth / fd 评价。
- snapshot 只保存轻量 `replayData`，默认 `saveRepeatDetail=false`，不保存完整 `disabledRepeatCell / gatedRepeatCell`。
- `replayData` 保留轻量表：`compareTable`、`inToothCompareTable`、`gateMissTable`、`aggregateTable`、`inToothAggregateTable`、`caseRoleSummaryTable`、`gateReasonTable`、`inToothGateReasonTable`、`finalTagTransitionTable`、`toothTransitionTable`、`candidateTable`、`warningTable`、`disabledRepeatTable`、`gatedRepeatTable`、`checkpointSummary`、`storageSummaryTable`。
- 可选 `notifyTelegramEnable=true`，仿真完成或失败后发送 best-effort HTML 通知；通知不影响 replay 结果或 snapshot 保存。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat` | `baseSeed=253`，`numRepeat=100`，`seedList=253:352`，`snrDb=10`，`focusedSeedList=[]`，`gatedRescueCoherenceThreshold=0.200`，`gatedRescuePhaseResidThresholdRad=1.000`，`saveRepeatDetail=false`，`checkpointEnable=true`，`notifyTelegramEnable=true` | gated `wide + single-MF` 能改善 median 并救回部分 hard tail，但当前 gate / adoption 不安全：overall 与 same-tooth 子集均出现 easy damage，P95 / max 变差。因此该结果是 negative / diagnostic，不应推广到默认 flow。 |

## 当前代表性结果

2026-05-05 的 100-repeat run 是当前代表性 flow-like 结果。snapshot 保存为：

```text
 test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat
```

运行配置摘要：

| 配置项 | 值 |
|---|---:|
| repeats | 100 |
| SNR | 10 dB |
| base seed | 253 |
| seed list | 253:352 |
| focused seed list | empty |
| repeat mode | `parfor-auto` |
| save snapshot | true |
| save repeat detail | false |
| checkpoint enable | true |
| checkpoint resume | true |
| Telegram notify | true |
| gated rescue coherence threshold | 0.200 |
| gated rescue phase threshold | 1.000 rad |

运行过程：disabled batch 命中已完成 checkpoint 并直接恢复 100/100 task；gated batch 从空 checkpoint 运行 100/100 task。完成后清理 disabled / gated checkpoint 目录并保存轻量 snapshot。运行中出现 1 个 seed 的 singular matrix warning，见“Warning”。

## Overall flow-like aggregate

| 指标 | disabled | gated | 变化 |
|---|---:|---:|---:|
| numSeed | 100 | 100 | — |
| triggerRate | — | 0.51 | — |
| selectedRate | — | 0.51 | — |
| hardRescueRate | — | 0.45946 | — |
| easyDamageRate | — | 0.11111 | 非零，不能默认推广 |
| fdDamageRate | — | 0 | 未观察到 fd damage |
| median angle (deg) | 0.0017468 | 0.00095966 | 改善 |
| P95 angle (deg) | 0.0035023 | 0.0050759 | 变差 |
| max angle (deg) | — | 0.010284 | 仍有大 tail |
| medianAngleGainDeg | — | 0 | 中位 gain 统计为 0 |

整体结果说明：`wide + single-MF` basin-entry family 仍有救力，能显著改善 median；但当前 gate / adoption 会放大 tail，且 easy damage 非零。因此这不是通过结果，而是一个高价值 negative / diagnostic 结果。

## Same-tooth in-tooth aggregate

same-tooth 子集定义为 `disabledToothIdx == 0 && gatedToothIdx == 0`，共 58 个 seed。该子集更接近当前 in-tooth DoA 主问题。

| 指标 | disabled | gated | 变化 |
|---|---:|---:|---:|
| numSeed | 58 | 58 | — |
| triggerRate | — | 0.5 | — |
| selectedRate | — | 0.5 | — |
| hardRescueRate | — | 0.38095 | — |
| easyDamageRate | — | 0.076923 | 非零，不能默认推广 |
| fdDamageRate | — | 0 | 未观察到 fd damage |
| median angle (deg) | 0.0015356 | 0.00091144 | 改善 |
| P95 angle (deg) | 0.0037733 | 0.0052297 | 变差 |
| max angle (deg) | — | 0.010284 | 仍有大 tail |

Same-tooth 子集也呈现同样模式：median 改善，但 P95 / max 变差，且 easy damage 非零。这说明当前问题不是单纯 wrong-tooth 污染；即使剥离到 same-tooth，当前 gate / adoption 仍不安全。

## Case-role summary

| caseRole | numSeed | sameToothRate | triggerRate | selectedRate | rescueRate | damageRate | disabledMedianAngleDeg | gatedMedianAngleDeg | medianAngleGainDeg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `hard-angle-tail` | 37 | 0.56757 | 0.56757 | 0.56757 | 0.45946 | 0 | 0.0028606 | 0.002058 | 0 |
| `middle-angle` | 45 | 0.53333 | 0.57778 | 0.57778 | 0 | 0 | 0.0014915 | 0.0011507 | 0 |
| `easy-angle` | 18 | 0.72222 | 0.22222 | 0.22222 | 0 | 0.11111 | 0.00077983 | 0.00074891 | 0 |

关键观察：

- hard tail 确实有 45.946% rescue rate，说明 candidate family 有救力。
- middle angle 没有 damage，但也不算 hard rescue。
- easy angle 触发率为 22.222%，damage rate 为 11.111%，这是当前 gate 不能默认推广的核心原因。

## Gate reason counts

### 全部 100 seed

| rescueGateReason | numSeed |
|---|---:|
| `coherence-not-collapsed` | 39 |
| `non-ref-coherence-collapse` | 33 |
| `non-ref-phase-residual-large` | 18 |
| `subset-not-trusted` | 10 |

### Same-tooth 58 seed

| rescueGateReason | numSeed |
|---|---:|
| `coherence-not-collapsed` | 22 |
| `non-ref-coherence-collapse` | 17 |
| `non-ref-phase-residual-large` | 12 |
| `subset-not-trusted` | 7 |

同齿内仍有较多 `coherence-not-collapsed` 样本，说明剩余 miss 不全是 non-ref coherence collapse hard case；后续必须单独诊断 ordinary / non-collapse angle miss。`subset-not-trusted` 则说明部分 collapse-like hard case 可能被 trust gate 拦截。

## Hard gate-miss cases

`gateMissTable` 只列出 same-tooth hard-angle-tail 且 rescue 未触发的 case，共 11 个 seed。

| seed | angle err (deg) | fdRef err (Hz) | coherence floor | gate reason | final tag |
|---:|---:|---:|---:|---|---|
| 256 | 0.004673 | 0.079397 | 0.07501 | `subset-not-trusted` | `periodic-static-seed` |
| 261 | 0.002058 | 7.7889 | 0.97856 | `coherence-not-collapsed` | `periodic-static-seed` |
| 285 | 0.0034602 | 0.14463 | 0.03552 | `subset-not-trusted` | `periodic-static-seed` |
| 298 | 0.0034128 | 0.38534 | 0.038124 | `subset-not-trusted` | `periodic-static-seed` |
| 305 | 0.0021977 | -0.011465 | 1 | `coherence-not-collapsed` | `periodic-subset-seed` |
| 315 | 0.0020672 | 0.13836 | 0.12249 | `subset-not-trusted` | `periodic-static-seed` |
| 317 | 0.0027852 | 1.0517 | 0.63926 | `coherence-not-collapsed` | `periodic-static-seed` |
| 322 | 0.0035444 | 5.2386 | 0.99221 | `coherence-not-collapsed` | `periodic-static-seed` |
| 323 | 0.0028606 | -52.94 | 0.93694 | `coherence-not-collapsed` | `periodic-static-seed` |
| 335 | 0.0020391 | -0.012703 | 1 | `coherence-not-collapsed` | `periodic-subset-seed` |
| 346 | 0.0032445 | 188.95 | 0.99999 | `coherence-not-collapsed` | `periodic-static-seed` |

这些 case 将后续问题拆成两类：

1. `subset-not-trusted`：例如 256 / 285 / 298 / 315。它们的 coherence floor 很低，更像 collapse-like hard case，但 trust gate 拦截了 rescue。下一步应检查 trust gate，而不是直接扩 candidate。
2. `coherence-not-collapsed`：例如 261 / 305 / 317 / 322 / 323 / 335 / 346。它们不是 coherence-collapse 型 hard case，应进入 ordinary angle miss replay，而不是继续调 collapse gate。

## Final tag transition

| disabledFinalTag | gatedFinalTag | numSeed |
|---|---|---:|
| `periodic-static-seed` | `periodic-static-seed` | 36 |
| `periodic-static-seed` | `same-tooth-rescue-single-mf` | 35 |
| `periodic-static-seed` | `same-tooth-rescue-wide` | 16 |
| `periodic-subset-seed` | `periodic-subset-seed` | 13 |

Gated 分支实际大量选择了 same-tooth rescue candidate：`single-MF` 35 次、`wide` 16 次。这说明 candidate family 的参与度足够高；当前主要问题不再是“candidate 完全没有被评估”，而是 gate / adoption 的安全性和 tail 控制。

## Tooth transition

Top transition：

| disabledToothIdx | gatedToothIdx | numSeed |
|---:|---:|---:|
| 0 | 0 | 58 |
| 1 | 1 | 7 |
| -2 | -2 | 3 |
| -1 | -1 | 3 |
| -19 | -19 | 2 |
| -5 | -5 | 2 |
| 4 | 4 | 2 |
| 5 | 5 | 2 |
| 75 | 75 | 2 |

没有观察到 gated rescue 改变 tooth index；它主要是在原 tooth 上改 final candidate。这也说明：本 replay 的好坏主要反映 same-tooth adoption / basin-entry 行为，而不是 rescue 直接修正 wrong-tooth。

## Warning

| seed | disabledWarningId | gatedWarningId | warning |
|---:|---|---|---|
| 273 | `MATLAB:singularMatrix` | `MATLAB:singularMatrix` | Matrix is singular to working precision. |

warning 同时出现在 disabled 和 gated 分支。当前将其作为诊断元数据保留，不单独解释为 gated rescue 造成的错误。

## Storage summary

当前 snapshot 仍只保存 `replayData` 一个 workspace 变量，但 `replayData` 内部包含足够重出 summary 的轻量表。

| field | class | size | repeat detail |
|---|---|---:|:---:|
| `compareTable` | table | 100 x 19 | false |
| `inToothCompareTable` | table | 58 x 19 | false |
| `gateMissTable` | table | 11 x 19 | false |
| `aggregateTable` | table | 1 x 14 | false |
| `inToothAggregateTable` | table | 1 x 14 | false |
| `caseRoleSummaryTable` | table | 3 x 10 | false |
| `gateReasonTable` | table | 4 x 2 | false |
| `inToothGateReasonTable` | table | 4 x 2 | false |
| `finalTagTransitionTable` | table | 4 x 3 | false |
| `toothTransitionTable` | table | 28 x 3 | false |
| `candidateTable` | table | 302 x 10 | false |
| `warningTable` | table | 1 x 5 | false |
| `disabledRepeatTable` | table | 100 x 36 | false |
| `gatedRepeatTable` | table | 100 x 36 | false |

该保存口径是合理的：文件不会因为 full repeat cell 变大，但 load 后仍能运行 `Summary output and plotting` 小节重打核心表。

## 当前结论

这次 100-repeat flow-like replay 是一个明确的 negative / diagnostic 结果：

- `wide + single-MF` basin-entry 仍然是有价值的 candidate family；它改善 median，并能救回部分 hard-angle-tail。
- 当前 gate / adoption 不够安全；overall 与 same-tooth 子集均有 easy damage，且 P95 / max 变差。
- 当前策略不能推进到默认 flow，也不能作为 full-flow 通过证据。
- 当前不应继续扩大 DoA step，也不应默认打开 joint fdRef / fdRate bank。
- 下一步应拆分排查：easy damage、`subset-not-trusted` hard miss、`coherence-not-collapsed` ordinary angle miss。

## 对主流程的影响

- 不接入默认 `runSimpleDynamicSubsetPeriodicFlow`。
- 不下沉到 `estimatorDoaDopplerMlePilotMfOpt`、objective 或 residual。
- 保留 `wide / single-MF` basin-entry 作为候选 family，但需要更严格或分层 gate。
- `non-ref-phase-residual-large` 相关触发需要重点复核，因为 easy-angle 中存在非零 damage。
- `subset-not-trusted` hard miss 需要单独 replay 判断：这些 case 是不是 gate 过严，而不是 candidate 不足。
- `coherence-not-collapsed` hard miss 应进入 `replayMfInToothOrdinaryAngleMissDiagnose`，不要继续用 collapse gate 硬修。

## 建议下一步

1. 先从本 snapshot 的 `compareTable` 中抽出 easy damage 样本，检查其 `rescueGateReason / disabledCoherenceFloor / gatedFinalTag / angleGainDeg`，判断误伤主要来自哪类 trigger。
2. 针对 `subset-not-trusted` hard miss seeds `256, 285, 298, 315`，做 controlled replay：固定同齿与 fd healthy 条件，比较放行 trust gate 后 `wide / single-MF` candidate 是否能救且是否误伤。
3. 针对 `coherence-not-collapsed` hard miss seeds `261, 305, 317, 322, 323, 335, 346`，新增或扩展 ordinary angle miss replay，比较 final-centered small polish、wide / single-MF basin-entry 与 truth-DoA oracle。
4. 在上述两类问题未分清前，不再扩大 `replayMfFlowLikeGatedBasinEntryEffectiveness` 的 repeat，也不推进默认 flow。
