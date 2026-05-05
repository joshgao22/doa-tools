# replayMfInToothDoaFdRangeEnvelope

## 对应 replay

- `test/dev/replay/replayMfInToothDoaFdRangeEnvelope.m`

## 观察目标

在 controlled in-tooth 条件下，从 `baseSeed` 自动搜索 hard / gate-miss / easy / fd-negative 代表 seed，再只对代表 seed 扫描 DoA-fdRef / DoA-fdRate 局部 objective envelope。该 replay 的目标不是验证最终 rescue 性能，而是回答三个问题：

1. hard seed 的可实现 candidate 是否存在；
2. 改善主要来自 DoA basin-entry center，还是来自非零 fdRef / fdRate joint offset；
3. 哪个 gate / safe-adopt policy 更值得进入 `replayMfInToothGatedRescueEffectiveness`。

该 replay 只做离线 surface / policy 诊断，不做 subset tooth selection，不改变 estimator / flow / objective，也不进入 regression。truth 只用于 seed 分类、oracle 上限和结果评价，不进入 runtime gate、candidate adoption 或 final winner。

## 最新代码口径

- 使用 truth-centered half-tooth `fdRef` 范围剥离 wrong-tooth。
- 默认 `seedMode="base-seed-search"`：先从 `baseSeed / numSearchRepeat` 跑轻量 scout，只对自动选出的 hard / gate-miss / easy / fd-negative 代表 seed 做重 envelope surface。
- 默认 `checkpointEnable=true`：scout 与 envelope 分阶段写 per-seed checkpoint，成功保存 `replayData` 后默认清理；中断或失败后重新运行可 resume。
- `static-ms / single-mf / wide` 为可实现 center；`truth` 只作 oracle center。
- 主要结果表：`selectedEnvelopeSeedTable`、`scoutTailDiagnosisTable`、`candidateCoverageTable`、`policyRecommendationTable`、`envelopeSummaryTable` 与 `rangeRecommendationTable`。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat` | `seedMode=base-seed-search`，`baseSeed=253`，`numSearchRepeat=100`，`snrDb=10`，`oracleFdHalfToothFraction=0.49`，`angleHitThresholdDeg=0.002` | 自动 seed selection 生效；hard seed 的可救方向主要来自 `wide / single-MF` DoA basin-entry center，而不是非零 fdRef joint step；policy sweep 支持 family-safe DoA step `0.006 deg`，不支持默认打开 joint fdRef / fdRate bank。 |

## 当前代表性结果

2026-05-04 的 base-seed envelope 是当前唯一保留的代表性结果。该结果先 scout `253:352` 共 100 个 seed，再自动选择 12 个 envelope seed 做 heavy surface。

### 运行与 checkpoint

| 阶段 | seed 数 | 运行时间 | checkpoint |
|---|---:|---:|---|
| scout pass | 100 | 约 28 min 51 s | 成功后清理 |
| envelope heavy pass | 12 | 约 4 min 26 s | 成功后清理 |

两阶段设计是有效的：它允许先从 100 个 seed 中发现代表性 hard / negative case，同时避免对所有 seed 做二维 / 三维 heavy surface。

### 自动选择的 envelope seed

| selection reason | seeds | 作用 |
|---|---|---|
| `hard-coherence-collapsed` | `256, 277, 345, 312, 298, 283` | 检查 coherence collapse hard case 是否存在可实现 basin-entry candidate。 |
| `hard-gate-miss` | `328, 319, 306` | 检查 hard 但 coherence gate 不触发的样本，评估是否需要新 gate proxy。 |
| `fd-not-healthy-negative` | `346, 285` | 检查 rescue / policy 是否误伤 fd 不健康负样本。 |
| `easy-negative` | `341` | 检查 rescue / policy 是否误伤 easy 样本。 |

这组 seed coverage 比固定手工 seed 更合理：它覆盖了已触发 hard、gate miss hard、fd-negative 和 easy-negative 四类样本。

## 主要统计与观察

### Candidate coverage

`candidateCoverageTable` 的核心现象是：hard seed 里确实存在可实现的好 candidate，但这些 candidate 主要来自 `wide-coarse-doa-grid` 或 `single-mf-coarse-doa-grid`，而不是非零 fdRef joint step。

| 观察项 | 结果 | 含义 |
|---|---|---|
| `jointAddsOverPure` | false | joint DoA-fdRef 没有比 pure DoA candidate 多救 seed。 |
| `jointUsesNonzeroFdRef` | false | best joint candidate 的 fdRef step 仍为 0。 |
| `bestJointFdRefStepHz` | 0 | 当前 evidence 不支持默认引入 fdRef offset。 |
| hard seed best center | 多为 `wide-coarse-doa-grid` / `single-mf-coarse-doa-grid` | 当前 tail 更像 DoA basin-entry center 问题。 |

代表性样本：

- seed `256` 的 best implementable center 是 `single-mf-coarse-doa-grid`，angle error 可降到约 `0.00034936 deg`，non-ref coherence 可恢复到约 `0.99845`。
- seed `277` 的 best implementable center 是 `wide-coarse-doa-grid`，angle error 可降到约 `3.6871e-05 deg`。
- seed `345` 的 best implementable center 是 `wide-coarse-doa-grid`，并与后续 `family-safe-adopt` 结果一致。

### DoA-fdRef / DoA-fdRate envelope

- DoA-fdRef surface 没有给出稳定的非零 fdRef offset envelope；有效 candidate 的 `bestFdRefOffsetHz` 通常仍为 0。
- DoA-fdRate surface 多处出现边界命中，说明 fdRate offset 方向尚未形成可下沉的稳定搜索范围。
- 因此，`replayMfCombToothSurface` 提示的 DoA-fdRef coupling 目前更适合作为 objective 机制解释，而不是直接进入默认 low-complexity rescue bank。

### Policy sweep

`policyRecommendationTable` 用 cached candidate 做离线 policy sweep，不重跑 estimator。当前最有价值的对比是固定 `coherenceThreshold=0.20` 时，把 family DoA step guard 从 `0.004 deg` 放宽到 `0.006 deg`。

| policy | selectedAngleHitRate | hardAngleHitRate | overallP95AngleDeg | overallMaxAngleDeg | easyDamageRate | fdNegativeDamageRate |
|---|---:|---:|---:|---:|---:|---:|
| step `0.004 deg` | 0.50 | 0.55556 | 0.0048848 | 0.0054539 | 0 | 0 |
| step `0.006 deg` | 0.58333 | 0.66667 | 0.0039798 | 0.0044191 | 0 | 0 |

该结果支持新增 `family-safe-adopt`：对 `wide / single-MF` basin-entry family 允许更大的 DoA step，而不是对所有 candidate blanket 放宽。

## 当前结论

1. `baseSeed` scout + 自动 seed selection 已经足够支撑后续分析，不需要继续用手工 seed list 作为默认入口。
2. hard seed 的好 candidate 主要来自 `wide / single-MF` DoA basin-entry center；当前结果不支持默认打开 DoA-fdRef joint bank。
3. fdRate surface 仍有边界命中，不适合进入默认 candidate bank。
4. policy sweep 支持 `familySafeAdoptMaxAbsDoaStepDeg=0.006`，并且在 selected envelope seeds 上未观察到 easy / fd-negative damage。
5. 下一步应把 `family-safe-adopt` 放进 `replayMfInToothGatedRescueEffectiveness` 验证 aggregate 表现；`Envelope` 本身暂时不需要扩大 seed 数或 fdRef / fdRate 范围。

## 对主流程的影响

- 保留 `replayMfInToothDoaFdRangeEnvelope` 作为搜索范围 / policy 发现器。
- 不把 joint fdRef / fdRate bank 默认接入 gated effectiveness 或 flow。
- 不修改 estimator 主 objective / residual。
- 若后续需要冲更高 DoA hit rate，应优先继续做 gate / adoption policy sweep，而不是扩大 fdRef / fdRate joint grid。
