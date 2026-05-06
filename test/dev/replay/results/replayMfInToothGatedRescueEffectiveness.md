# replayMfInToothGatedRescueEffectiveness 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `candidate-but-not-target` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` |
| 当前一句话结论 | 在 controlled in-tooth 条件下，family-safe `wide + single-MF` basin-entry 能显著压低 hard-collapse DoA tail，hard rescue rate 达到 `0.8`，easy / fd-negative damage 为 0。 |
| 决策影响 | 支持推进到 flow-like gated same-tooth basin-entry 验证；暂不进入 estimator 主核或 flow 默认路径。 |
| 下一步动作 | 保留该 snapshot；下一步只在真实 subset-periodic flow 中验证 gate / adoption / wrong-tooth interaction。 |
| 禁止误用 | 不能当作 full-flow 已通过；不能据此默认打开 joint fd bank；不能继续为了 99% hit rate 放宽 coherence gate 或 blanket wide。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothGatedRescueEffectiveness.m`
- 结果文档：`test/dev/replay/results/replayMfInToothGatedRescueEffectiveness.md`
- replay 类型：controlled in-tooth 小 MC replay。
- 主要问题：same-tooth、fd 已健康但 non-ref coherence collapse 的 hard case，是否能由 truth-free gated `wide + single-MF` basin-entry bank 救回，同时避免 easy / fd-negative damage。
- 观察范围：truth-centered half-tooth `fdRef` 范围；比较 disabled、strict safe-adopt、family-safe-adopt 与 truth-DoA oracle。
- 不覆盖范围：不做 subset tooth selection；不验证 full-flow；不改变 estimator / objective / residual；不证明论文 CRB consistency。
- truth 使用口径：`caseRole / isHardRescued / isDamaged` 可用 truth 做离线评价；`rescueTriggered / triggerReason / selectedCandidateFamily` 必须来自 default result、non-ref coherence、candidate objective 和 no-truth 诊断。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `coherence-v1` | 只由 default non-ref coherence collapse 触发的 gate；phase residual 暂为诊断占位。 | No | 决定是否启用 rescue bank。 | 当前 gate 很干净，easy / fd-negative trigger 为 0。 |
| `hard-collapse` | same-tooth / fd healthy，但 non-ref coherence floor 很低且 angle miss 的 hard case。 | Evaluation only | 只改离线分类。 | 与 ordinary non-collapse miss 分开处理。 |
| `strict safe-adopt` | DoA-only bank，统一 `safeAdoptMaxAbsDoaStepDeg=0.004`。 | No | 改 candidate adoption。 | 能救一部分 hard case，但会拒绝部分有效 wide / single-MF center。 |
| `family-safe-adopt` | 对 `wide / single-MF` basin-entry family 放宽 DoA step guard 到 `0.006 deg`；default / final-centered 仍用 strict guard。 | No | 改 family-aware adoption。 | 当前 controlled in-tooth 最优 replay-level candidate。 |
| `wide / single-MF basin-entry` | 从 wide DoA center 或 single-MF center 重新进入 basin 的候选 family。 | No | 改 DoA basin center。 | 对 hard-collapse 有救力，但只能 gated 使用，不能 blanket 常驻。 |
| `truth-DoA oracle` | 离线 truth-centered 上限。 | Oracle only | 只改评价上限。 | 说明 tail 可达，不是模型无解；不能进入 runtime。 |
| `damage` | 对 easy / fd-negative 样本造成角度变差或越过阈值。 | Evaluation only | 只改评价标签。 | damage 为 0 是推进到下一层验证的必要条件。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` | 2026-05-04 | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=100`，`coherenceThreshold=0.20`，`safeStep=0.004 deg`，`familyStep=0.006 deg`，`safeCoh=0.95`，`safeFdRefStep=300 Hz`，`joint=0`，`family=1`。 | family-safe-adopt 同时改善 RMSE / median / P95 / P99 / max / hit rate；hard rescue rate `0.8`，easy / fd-negative damage 为 0。 | 当前唯一代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 100`
- seed range：`253:352`
- `fd oracle half-tooth fraction = 0.49`
- `fdRate oracle half-width = 1000 Hz/s`
- `gated rescue coherence threshold = 0.20`
- `safe adopt max DoA step = 0.004 deg`
- `family safe max DoA step = 0.006 deg`
- `safe adopt coherence threshold = 0.95`
- `safe adopt max fdRef step = 300 Hz`
- `angleHitThresholdDeg = 0.002`
- `includeJointSafeAdoptBank = false`
- checkpoint：`100/100` 完成，成功保存 `replayData` 后清理 checkpoint run 目录。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| hard cases | 5 | non-ref coherence collapse hard 样本数。 |
| easy negative | 90 | raw 已健康、用于检查误触发 / damage 的样本。 |
| fd-negative | 5 | fd 不健康负样本。 |
| family-safe trigger rate | 0.04 | 只触发 4% case，未 blanket 改所有样本。 |
| family-safe hard trigger | 0.8 | 5 个 hard case 中触发 4 个。 |
| family-safe hard rescue | 0.8 | hard rescue rate 明显高于 strict 的 0.4。 |
| easy damage | 0 | 没有误伤 raw easy 样本。 |
| fd-negative damage | 0 | 没有误伤 fd-negative 样本。 |
| overall hit rate | 0.94 | raw `0.90` 提升到 `0.94`。 |
| hard hit rate | 1.00 | hard case 在 family-safe 后全部进入 `0.002 deg` hit threshold。 |

### 4.3 关键对比表

#### In-tooth method aggregate

| method | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 |
| `ms-sf-static` | 0.99 | 0.0020359 | 0.0015905 | 0.0034644 | 0.0045189 | 0.66 | 0.004673 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 |
| `ms-mf-cp-u-in-tooth` | 1.00 | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 |
| `ms-mf-cp-u-wide-doa-in-tooth` | 1.00 | 0.0021782 | 0.00044087 | 0.0037959 | 0.011068 | 0.93 | 0.014519 |
| `ms-mf-cp-u-truth-doa-oracle` | 1.00 | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 |

#### In-tooth DoA ladder

| stage | method | trigger | selected by bank | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | `ss-mf-cp-u-in-tooth` | — | — | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 |
| 4 | `ms-mf-cp-u-in-tooth` | — | — | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 |
| 5a | strict DoA safe-adopt | 0.04 | 0.02 | 0.0012688 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.0054541 |
| 5b | family-safe-adopt | 0.04 | 0.04 | 0.0010127 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 0.0044189 |
| 6 | truth-DoA oracle | — | — | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 |

#### Gated rescue aggregate

| bank | hard cases | easy negative | fd-negative | trigger rate | hard trigger | hard rescue | easy damage | fd-negative damage | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | hard hit rate | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `disabled` | 5 | 90 | 5 | 0 | 0 | 0 | 0 | 0 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.20 | 0.0054541 |
| strict safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.4 | 0 | 0 | 0.00048794 | 0.0022989 | 0.0053885 | 0.92 | 0.60 | 0.0054541 |
| family-safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.8 | 0 | 0 | 0.00047908 | 0.0020422 | 0.0036254 | 0.94 | 1.00 | 0.0044189 |

#### Verdict

| bank | hardRescueRate | easyDamageRate | fdNegativeDamageRate | P95 improvement (deg) | max improvement (deg) | recommendedPass | recommendation |
|---|---:|---:|---:|---:|---:|:---:|---|
| strict safe-adopt | 0.4 | 0 | 0 | 0.0013185 | 0 | false | `hold` |
| family-safe-adopt | 0.8 | 0 | 0 | 0.0015752 | 0.0010353 | true | `candidate-but-not-target` |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `phaseResidAvailable=false`，当前 gate 实际只依赖 non-ref coherence collapse；即便如此，trigger 仍很干净。
- family-safe-adopt 只触发 `4%` case，却把 hard rescue rate 从 strict 的 `0.4` 提到 `0.8`。
- easy / fd-negative damage 均为 0，说明 family-safe 不是 blanket 放宽。
- raw `MS-MF-CP-U-in-tooth` 的 median 已很好，但 P95 / P99 / max 被少数 hard seed 拉坏；family-safe 明显压低 tail。

### 5.2 仍未解决或反向的现象

- overall hit rate 只从 `0.90` 到 `0.94`，没有达到 `0.99`；剩余 miss 不全是 hard-collapse。
- seed `253` 没有触发 gate，因为 coherence floor `0.20455` 略高于阈值 `0.2`；但它本身已满足 `0.002 deg` hit threshold。
- 该 replay 没有验证真实 subset selection 下的 wrong-tooth interaction，因此不能直接推进 flow default。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 256 | hard-collapse | default angle `0.0054541 deg`，coherence floor `0.00067049`；family-safe 选择 `single-mf-coarse-doa-grid`，angle 降到 `0.00034933 deg`，coherence 恢复到 `0.99845`。 | 证明 strict `0.004 deg` guard 会拒绝有效 candidate，family `0.006 deg` guard 有必要。 |
| 345 | hard-collapse | default angle `0.003759 deg`，coherence floor `0.10043`；family-safe 选择 `wide-coarse-doa-grid`，angle 降到 `7.5606e-05 deg`，coherence 恢复到 `1`。 | 证明 wide center 对 hard-collapse 仍稳定有效。 |
| 253 | gate-miss but hit | default angle `0.001526 deg`，truth oracle `0.00048515 deg`，coherence floor `0.20455`，未触发。 | 说明不应为了该 seed 盲目调高 coherence threshold；它已经是 hit。 |

## 6. 机制解释

### 6.1 当前解释

hard-collapse 的核心症状不是 wrong-tooth，而是 same-tooth 内 non-ref coherence collapse。raw `MS-MF-CP-U-in-tooth` 的频率链和 tooth 已经健康，但 DoA / local-state 被带入坏 basin，导致非参考星相干性塌陷。`wide / single-MF` basin-entry 的作用是提供新的 DoA center，让 estimator 从另一个 basin 进入；family-safe-adopt 则只对这类 basin-entry family 放宽 DoA step guard，避免 strict guard 拒绝真正有效的 candidate。

这个结果与 ordinary-wide 线必须分开解释：hard-collapse 的 gate 依赖 non-ref coherence collapse；ordinary miss 往往 coherence 为 1，不能用同一个 gate 解决。

### 6.2 这个结果支持什么

- 支持 family-safe-adopt 作为 controlled in-tooth replay-level candidate。
- 支持下一步进入 flow-like gated same-tooth basin-entry 验证。
- 支持继续默认旁路 joint fdRef / fdRate bank，因为本轮收益来自 DoA basin-entry，不来自 nonzero fd joint offset。

### 6.3 这个结果不证明什么

- 不证明 estimator 默认路径已经修复。
- 不证明 full-flow 性能已经可用于论文主图。
- 不证明 ordinary non-collapse miss 也能由 hard-collapse gate 解决。
- 不证明可以写 regression 契约。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉 objective / residual / estimator 主核。 |
| flow 默认路径 | 暂不改；下一步放入 flow-like replay 验证。 |
| regression | 暂不写；只有真实 flow 稳定后再考虑 gate/no-truth 契约。 |
| replay / scan 下一步 | 用 flow-like replay 或 `runSimpleDynamicSubsetPeriodicFlow` 的 final same-tooth rescue 层比较 disabled / strict / family-safe。 |
| 论文图 / 论文口径 | diagnostic / optional robustness；不作为正文主精度指标。 |
| 排障记录 | 作为 controlled in-tooth hard-collapse candidate 证据保留。 |

## 8. 限制与禁止解释

- 不要把 `candidate-but-not-target` 解读成 default target。
- 不要因为 seed `253` gate miss 就调高 coherence threshold；该 seed 已满足 hit threshold。
- 不要 blanket 打开 wide 或 single-MF candidate。
- 不要默认打开 joint fdRef / fdRate bank。
- 不要把当前 hard-collapse 结果外推到 ordinary non-collapse miss。
- 不要把 truth label 放入 runtime selector、gate、candidate adoption 或 final winner。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothGatedRescueEffectiveness.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 snapshot 固定了 controlled in-tooth hard-collapse 线的当前最优候选：family-safe-adopt。
- 后续 ordinary-wide 500-repeat 已说明 non-collapse ordinary miss 不能并入本 gate；两条线应保持分开。
