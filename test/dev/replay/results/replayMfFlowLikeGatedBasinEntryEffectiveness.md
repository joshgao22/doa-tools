# replayMfFlowLikeGatedBasinEntryEffectiveness 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / negative` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat` |
| 当前一句话结论 | flow-like 真实 subset-periodic 结构中，`wide + single-MF` basin-entry family 仍有救力并改善 median，但当前 gate / adoption 会产生 easy damage，并使 P95 / max 变差，因此不能推进默认 flow。 |
| 决策影响 | 保留为 flow-like stress replay 和 negative diagnostic；不接入默认 `runSimpleDynamicSubsetPeriodicFlow`，不下沉 estimator / objective / residual。 |
| 下一步动作 | 拆分 easy damage、`subset-not-trusted` hard miss 与 `coherence-not-collapsed` ordinary miss；在这些子问题未分清前，不扩大 repeat，也不推进默认 gate。 |
| 禁止误用 | 不能把 controlled in-tooth positive result 直接外推到 full-flow；不能把当前 gate 写成 regression；不能因为 median 改善就忽略 P95 / max / easy damage。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfFlowLikeGatedBasinEntryEffectiveness.m`
- 结果文档：`test/dev/replay/results/replayMfFlowLikeGatedBasinEntryEffectiveness.md`
- replay 类型：小 MC / flow-like stress replay。
- 主要问题：在真实 `subset-periodic flow` 结构内，no-truth gated `wide + single-MF` same-tooth basin-entry rescue 是否仍能救 hard tail，并且不误伤 easy case。
- 观察范围：100-repeat 小 MC；disabled baseline 与 gated flow-like 分支使用相同 seed list；gated 分支通过 `buildSimpleDynamicFlowOpt.sameToothRescue` 打开 basin-entry candidate。
- 不覆盖范围：不提供 oracle 上限；不跳过 subset selection；不证明 full-flow 已通过；不改变 estimator 默认路径；不证明论文主图可直接使用 full-flow 统计。
- truth 使用口径：truth 只用于结果表中的 angle / tooth / fd 评价、case role 标注和 damage/rescue 离线统计；不进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `flow-like` | 保留 subset selection、periodic refine、same-tooth rescue gate 等真实 flow 结构的 replay。 | No | 改变 replay 运行环境，不跳过真实 flow 决策。 | 比 controlled in-tooth 更接近默认 flow，但也会混入 wrong-tooth、subset trust 和 adoption 风险。 |
| `disabled` | 不打开 same-tooth basin-entry rescue 的 baseline flow。 | No | 作为对照基线。 | 用于判断 gated 分支带来的 rescue、damage 和 tail 变化。 |
| `gated` | 打开 no-truth same-tooth rescue gate 的 flow-like 分支。 | No | 允许 flow 在同一 tooth 内选择 `wide / single-MF` basin-entry candidate。 | 当前能改善 median，但 gate/adoption 不安全。 |
| `wide + single-MF basin-entry` | 以较宽 DoA 或单星多帧 candidate 改变 DoA basin center 的 candidate family。 | No | 改变最终 DoA basin entry，而不是修改 objective。 | family 本身有救力；是否可进入默认 flow 取决于 gate 与 damage 控制。 |
| `non-ref-coherence-collapse` | 非参考星 coherence 明显塌陷的 no-truth 触发原因。 | No | 触发 same-tooth rescue。 | 可解释部分 hard-collapse case，但不是所有 ordinary miss。 |
| `non-ref-phase-residual-large` | 非参考星 phase residual 较大的 no-truth 触发原因。 | No | 触发 same-tooth rescue。 | 当前需要重点复核，因为 easy damage 可能来自过宽触发。 |
| `subset-not-trusted` | subset trust 诊断认为当前 selected subset 不可靠。 | No | 阻止或限制 rescue。 | 部分 hard-collapse case 可能被该 gate 过严拦截。 |
| `coherence-not-collapsed` | coherence 未塌陷，因此 collapse gate 不触发。 | No | 阻止 collapse 型 rescue。 | 这类 hard miss 更像 ordinary / non-collapse angle miss，应由 ordinary angle replay 继续拆分。 |
| `hard-angle-tail` | disabled baseline angle 已进入 tail 的离线 case role。 | Evaluation only | 只改结果分类。 | 用于观察 rescue 是否救 hard case；不能进入 runtime gate。 |
| `easy-angle` | disabled baseline 已较好的离线 easy case。 | Evaluation only | 只改 damage 评价。 | easy damage 是当前 gate 不能默认推广的核心风险。 |
| `toothIdx` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于评价和结果标注。 | 本 replay 中 gated 基本不改变 tooth，主要观察 same-tooth adoption。 |
| `damage` | gated 分支使原本 easy / baseline-hit 样本变差。 | Evaluation only | 只改结果评价。 | damage 优先级高于单纯 median gain 或 rescue rate。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat` | 2026-05-05 | representative | `baseSeed=253`，`numRepeat=100`，`snrDb=10`，`seedList=253:352`，`gatedRescueCoherenceThreshold=0.200`，`gatedRescuePhaseResidThresholdRad=1.000`，`checkpointEnable=true`，`saveRepeatDetail=false` | basin-entry family 有救力，但 gate/adoption 不安全：median 改善，easy damage 非零，P95 / max 变差。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 100`
- seed range：`253:352`
- `focusedSeedList = []`
- repeat mode：`parfor-auto`
- gated rescue coherence threshold：`0.200`
- gated rescue phase residual threshold：`1.000 rad`
- checkpoint：disabled batch 命中已完成 checkpoint 并恢复 `100/100` task；gated batch 从空 checkpoint 完成 `100/100` task；成功后清理 checkpoint run 目录。
- snapshot 保存变量：`replayData`
- repeat detail：`saveRepeatDetail=false`，仅保存轻量 table。
- Telegram：`notifyTelegramEnable=true`，best-effort completion / failure notification。

运行中出现 1 个 seed 的 singular matrix warning，且 disabled 与 gated 分支同时出现，当前仅作为诊断元数据保留。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `numSeed` | 100 | 小 MC flow-like stress 样本数。 |
| `triggerRate` / `selectedRate` | 0.51 / 0.51 | gated 分支约一半 seed 触发并选择 rescue candidate，说明 candidate family 参与度足够高。 |
| `hardRescueRate` | 0.45946 | hard-angle-tail 中约 46% 被救回，说明 basin-entry family 有实际救力。 |
| `easyDamageRate` | 0.11111 | easy case 出现非零误伤，是当前不能默认推广的主要原因。 |
| `fdDamageRate` | 0 | 未观察到 fd damage，但不足以抵消 angle easy damage。 |
| disabled / gated median angle | 0.0017468 / 0.00095966 deg | median 明显改善。 |
| disabled / gated P95 angle | 0.0035023 / 0.0050759 deg | tail 变差，说明 gate/adoption 放大了部分坏样本。 |
| gated max angle | 0.010284 deg | 仍有大 tail，不能作为 flow/default 通过证据。 |

### 4.3 关键对比表

#### Overall flow-like aggregate

| 方法 | numSeed | trigger / selected | hard rescue | easy damage | fd damage | median angle (deg) | P95 angle (deg) | max angle (deg) | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| disabled | 100 | — | — | — | — | 0.0017468 | 0.0035023 | — | baseline flow |
| gated | 100 | 0.51 / 0.51 | 0.45946 | 0.11111 | 0 | 0.00095966 | 0.0050759 | 0.010284 | median 改善，但 P95 / max / easy damage 变差 |

#### Same-tooth in-tooth aggregate

same-tooth 子集定义为 `disabledToothIdx == 0 && gatedToothIdx == 0`，共 58 个 seed。该子集更接近当前 in-tooth DoA 主问题。

| 方法 | numSeed | trigger / selected | hard rescue | easy damage | fd damage | median angle (deg) | P95 angle (deg) | max angle (deg) | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| disabled | 58 | — | — | — | — | 0.0015356 | 0.0037733 | — | same-tooth baseline |
| gated | 58 | 0.50 / 0.50 | 0.38095 | 0.076923 | 0 | 0.00091144 | 0.0052297 | 0.010284 | 剥离 wrong-tooth 后仍有 easy damage 与 tail 变差 |

#### Case-role summary

| caseRole | numSeed | sameToothRate | triggerRate | selectedRate | rescueRate | damageRate | disabled median (deg) | gated median (deg) | 解释 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `hard-angle-tail` | 37 | 0.56757 | 0.56757 | 0.56757 | 0.45946 | 0 | 0.0028606 | 0.002058 | hard tail 确有 rescue 价值 |
| `middle-angle` | 45 | 0.53333 | 0.57778 | 0.57778 | 0 | 0 | 0.0014915 | 0.0011507 | 有 median 改善，但不属于 hard rescue |
| `easy-angle` | 18 | 0.72222 | 0.22222 | 0.22222 | 0 | 0.11111 | 0.00077983 | 0.00074891 | easy trigger 与 damage 非零，是当前 blocker |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- Overall 与 same-tooth 子集都显示 median 改善，说明 `wide + single-MF` basin-entry family 不是无效候选。
- hard-angle-tail 的 rescue rate 为 `0.45946`，说明真实 flow-like 结构内仍能救回部分 hard tail。
- final tag transition 中，gated 分支大量选择 same-tooth rescue candidate：`same-tooth-rescue-single-mf` 35 次，`same-tooth-rescue-wide` 16 次；当前主要问题不是 candidate 完全没被评估。
- tooth transition 的 top row 是 `0 -> 0` 共 58 个 seed，且未观察到 gated rescue 改变 tooth index；本 replay 主要反映 same-tooth adoption / basin-entry 行为，而不是直接修 wrong-tooth。

### 5.2 仍未解决或反向的现象

- Overall easy damage rate 为 `0.11111`，same-tooth 子集 easy damage rate 为 `0.076923`，说明 gate / adoption 不安全。
- Overall P95 从 `0.0035023 deg` 变差到 `0.0050759 deg`，same-tooth P95 从 `0.0037733 deg` 变差到 `0.0052297 deg`，说明 median gain 掩盖了 tail 放大。
- `coherence-not-collapsed` 占全部 39 个 seed、same-tooth 22 个 seed，说明剩余 miss 不全是 non-ref coherence collapse hard case。
- `subset-not-trusted` 有全部 10 个 seed、same-tooth 7 个 seed；其中部分 hard-collapse-like case 可能被 trust gate 过严拦截。
- `non-ref-phase-residual-large` 触发了 18 个 seed，且 easy damage 非零，说明 phase-residual 触发条件需要单独复核。

### 5.3 代表性 seed / case

`gateMissTable` 列出 same-tooth hard-angle-tail 且 rescue 未触发的 11 个 seed。它们把后续问题拆成两类。

| seed | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 256 | `subset-not-trusted` hard miss | angle `0.004673 deg`，coherence floor `0.07501`，final tag `periodic-static-seed` | coherence 很低但被 trust gate 拦截；后续应检查 trust gate 是否过严。 |
| 285 | `subset-not-trusted` hard miss | angle `0.0034602 deg`，coherence floor `0.03552` | 与 seed 256 同类，支持单独做 trust gate 放行诊断。 |
| 298 | `subset-not-trusted` hard miss | angle `0.0034128 deg`，coherence floor `0.038124` | collapse-like hard case 未触发 rescue。 |
| 315 | `subset-not-trusted` hard miss | angle `0.0020672 deg`，coherence floor `0.12249` | 同类 trust-gate hard miss。 |
| 261 | `coherence-not-collapsed` hard miss | angle `0.002058 deg`，fdRef err `7.7889 Hz`，coherence floor `0.97856` | 不是 coherence-collapse 型 hard case，应进入 ordinary angle miss 分类。 |
| 305 | `coherence-not-collapsed` hard miss | angle `0.0021977 deg`，coherence floor `1`，final tag `periodic-subset-seed` | 高 coherence 仍有 angle miss，不能继续靠 collapse gate 硬修。 |
| 323 | `coherence-not-collapsed` hard miss | angle `0.0028606 deg`，fdRef err `-52.94 Hz`，coherence floor `0.93694` | ordinary / non-collapse miss 的代表之一。 |
| 346 | `coherence-not-collapsed` hard miss | angle `0.0032445 deg`，fdRef err `188.95 Hz`，coherence floor `0.99999` | 高 coherence + fd 偏差并存，需单独诊断。 |
| 273 | warning seed | disabled 与 gated 均出现 `MATLAB:singularMatrix` | warning 不应单独归因于 gated rescue。 |

## 6. 机制解释

### 6.1 当前解释

这轮 replay 与 controlled in-tooth positive result 的差异非常关键。controlled in-tooth replay 已经剥离 wrong-tooth，并在固定 half-tooth / oracle 约束下观察到 `wide + single-MF` basin-entry family 可以救 hard-collapse tail；而本 replay 把同类 family 放回真实 subset-periodic flow 后，会同时面对 subset trust、same-tooth gate、ordinary non-collapse miss 和 final adoption 之间的耦合。

结果显示，candidate family 的方向没有错：它能改善 median，并救回约一半 hard-angle-tail。但是当前 no-truth gate 太粗，无法只捕捉真正需要 basin-entry 的 hard-collapse case。一部分 easy case 被触发后遭到 damage；另一部分 hard case 因 `subset-not-trusted` 或 `coherence-not-collapsed` 没有触发。也就是说，当前 blocker 不是“candidate 不存在”，而是 **gate / adoption 不能稳定区分 hard-collapse、ordinary miss 与 easy case**。

same-tooth 子集也出现同样模式，说明问题不只是 wrong-tooth 污染。即使 tooth 不变，当前 gated candidate 仍可能改善 median 但放大 tail。这支持把后续工作拆成：easy damage 溯源、trust gate hard miss、ordinary angle miss，而不是继续扩大 DoA step 或 blanket 打开 rescue。

### 6.2 这个结果支持什么

- 支持 `wide + single-MF` basin-entry family 在真实 flow-like 结构内仍有救力。
- 支持把 flow-like replay 作为 controlled replay 与 default flow 之间的必要 stress check。
- 支持将 hard-collapse 与 ordinary non-collapse angle miss 分开处理。
- 支持暂缓扩大 repeat，先做 seed-level / role-level 拆分诊断。

### 6.3 这个结果不证明什么

- 不证明 `family-safe-adopt` 可以进入默认 flow。
- 不证明 same-tooth rescue gate 已经 no-truth 安全。
- 不证明 full-flow CP/IP 或论文主图可以使用当前 flow-like 结果。
- 不证明可以写 regression 契约；当前统计仍是 negative / diagnostic。
- 不证明需要扩大 DoA step、提高 coherence threshold，或默认打开 joint `fdRef/fdRate` bank。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉到 `estimatorDoaDopplerMlePilotMfOpt`、objective 或 residual。 |
| flow 默认路径 | 不接入默认 `runSimpleDynamicSubsetPeriodicFlow`；保留为 gated rescue 候选 family 的 flow-like stress 证据。 |
| regression | 不写；当前 gate/adoption 尚未稳定，damage 非零。 |
| replay / scan 下一步 | 从本 snapshot 抽 easy damage 样本、`subset-not-trusted` hard miss、`coherence-not-collapsed` ordinary miss 分别做 targeted replay。 |
| 论文图 / 论文口径 | diagnostic-only / stress test；不能作为 paper-facing 精度图。论文主图仍应优先用 CRB consistency、resolved-regime RMSE、resolved/outlier rate、CP/IP controlled trade-off。 |
| 排障记录 | 可在机制归并版中保留一句：flow-like gated basin-entry 当前为 negative diagnostic，说明 controlled in-tooth positive result 不能直接推广默认 flow。 |

## 8. 限制与禁止解释

- 不要把 controlled in-tooth `family-safe-adopt` 的正结果直接推广到真实 full-flow。
- 不要因为 median angle 改善就忽略 P95 / max / easy damage；当前 negative 结论主要来自 tail 与 damage。
- 不要把 truth label、truth tooth、truth DoA、truth `fdRef/fdRate` 放入 runtime selector、gate、candidate adoption 或 final winner。
- 不要把当前 gate 写成 regression；它还不能稳定避免 easy damage。
- 不要继续 blanket 扩大 `wide + single-MF` DoA step 或默认打开 joint `fdRef/fdRate` bank。
- 不要把 `coherence-not-collapsed` hard miss 继续强行归入 collapse-gated rescue；它应进入 ordinary angle miss / non-collapse 分类。
- 不要把 `subset-not-trusted` hard miss 简单解释为 candidate 不足；它可能是 trust gate 过严。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
test/dev/replay/replayMfFlowLikeGatedBasinEntryEffectiveness.m
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 snapshot 是当前唯一代表性 flow-like gated basin-entry 结果。
- 它与 controlled `replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` 的关系是：controlled replay 支持 family 有救力，本 flow-like replay 暴露 gate/adoption 在真实 flow 中不安全。
- 在 easy damage、trust gate hard miss 与 ordinary non-collapse miss 未拆清前，不建议继续扩大该 replay 的 repeat。
