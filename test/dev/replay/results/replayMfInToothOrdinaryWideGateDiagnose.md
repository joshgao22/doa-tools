# replayMfInToothOrdinaryWideGateDiagnose 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `hold / partial diagnostic` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-233931.mat` |
| 当前一句话结论 | 500-repeat 扩大验证推翻 100-repeat 的 preferred-candidate 乐观结论；ordinary-wide gate 低误伤但 coverage 不足，只能作为 ordinary non-collapse miss 的 bad-basin 诊断。 |
| 决策影响 | 不推进 estimator / flow 默认路径，不写 regression；保留为 threshold / outlier 机制证据。 |
| 下一步动作 | 若继续维护，只做 outcome classification、baseline-relative 与 CRB-normalized 诊断；不继续调二维 gate 冲 `99%` hit rate。 |
| 禁止误用 | 不能把 `wide-obj10-minStep0.001-maxStep0.004` 视为 default candidate；不能用 truth label 做 runtime gate；不能用本 replay 证明 full-flow 已修好。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothOrdinaryWideGateDiagnose.m`
- 结果文档：`test/dev/replay/results/replayMfInToothOrdinaryWideGateDiagnose.md`
- replay 类型：oracle / controlled replay + 小 MC gate diagnostic。
- 主要问题：controlled in-tooth 条件下，ordinary non-collapse DoA miss 是否能用 truth-free `wide` basin-entry gate 接住，并避免 easy / fd-negative damage。
- 观察范围：truth-centered half-tooth `fdRef` 范围内的 `MS-MF-CP-U-in-tooth` tail；重点看 ordinary miss、wide candidate feature、gate policy sweep 和 damage。
- 不覆盖范围：不做 subset tooth selection；不验证真实 subset-periodic full-flow；不改变 estimator / flow / objective；不证明 CRB consistency。
- truth 使用口径：truth 只用于 offline label、ordinary / collapse / fd-negative 分类、truth-DoA oracle 上限和结果评价；gate policy 的输入必须来自 default / wide candidate 的 truth-free feature。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `ordinary-angle-miss` | same tooth、fd healthy、non-ref coherence 不 collapse，但 raw angle 超过 `0.002 deg` 的普通 DoA miss。 | Evaluation only | 只改离线分类。 | 用于把 non-collapse miss 从 hard-collapse / fd-negative 中拆出来，不能作为 runtime gate。 |
| `wide-basin-entry` | 从较宽 DoA center 重新进入局部 basin 的候选。 | No | 改 DoA basin center。 | 对部分 ordinary miss 有救力，但 blanket 使用会 damage；不能直接下沉默认路径。 |
| `ordinary-wide gate` | 用 objective gain、wide DoA step、candidate coherence 等 truth-free proxy 判断是否采用 wide candidate。 | No | 改 candidate adoption。 | 当前只形成低误伤 partial rescue，recommendation 为 `hold`。 |
| `best-implementable` | 在当前可实现候选中按 replay 评价选出的上限，不包含 truth-DoA oracle。 | Evaluation only | 只改结果汇总。 | 说明候选 family 的可实现上限，不等于 runtime selector。 |
| `truth-doa-oracle` | 把 DoA 放到 truth 附近的离线上限候选。 | Oracle only | 只用于上限诊断。 | 用来判断模型 / local objective 是否可达；不能进入 runtime flow。 |
| `easy-hit damage` | wide gate 使原本 raw 已 hit 的样本变差并越过阈值。 | Evaluation only | 只改评价标签。 | damage 优先级高于 rescue rate；出现 damage 的 policy 不能推进默认路径。 |
| `fd-negative` | fdRef / fdRate 不健康的负样本。 | Evaluation only | 只改分类。 | 不应被 ordinary DoA gate 混入正样本，否则会把频率问题误判成 DoA gate 问题。 |
| `toothIdx` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于评价和结果标注。 | 本 replay 中 `wrongToothCount=0`，说明已剥离 wrong-tooth 污染。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-233931.mat` | 2026-05-05 | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=500`，`seed=253:752`，wide gate sweep，默认旁路 final-small-polish / single-MF ordinary probe / hard-collapse auxiliary tables。 | `wide-obj10-minStep0.001-maxStep0.004` 低误伤但只救约 `59%` ordinary miss，recommendation 为 `hold`。 | 覆盖 100-repeat preferred-candidate 结论。 |
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-202052.mat` | 2026-05-05 | superseded | 100-repeat，ordinary miss 只有 4 个。 | 曾显示 preset policy 像 preferred candidate。 | 被 500-repeat 推翻。 |
| `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-181942.mat` | 2026-05-05 | superseded | 首轮 wide-gate sweep。 | 证明 `obj0` 可救 ordinary 但 easy trigger 过宽。 | 被后续 gate sweep 覆盖。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 500`
- seed range：`253:752`
- `fd oracle half-tooth fraction = 0.49`
- `fdRate oracle half-width = 1000 Hz/s`
- `ordinaryCoherenceMinThreshold = 0.95`
- `angleHitThresholdDeg = 0.002`
- `includeFinalSmallPolishProbe = false`
- `includeSingleMfBasinEntryProbe = false`
- `includeHardCollapseAuxTables = false`
- `includeJointSafeAdoptBank = false`
- `includeFamilySafeAdoptBank = false`
- wide gate sweep：`objGain=[0 5 10 20 30 50 80 100 150 200]`，`minStep=[0 0.001 0.0015 0.002] deg`，`maxStep=[0.004 0.006] deg`
- checkpoint：`500/500` 完成，成功保存 `replayData` 后清理 checkpoint run 目录。
- 运行时间：约 `1 h 44 min`。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `numSeed` | 500 | 当前代表性小 MC 总数。 |
| `numAngleMiss` | 65 | raw `MS-MF-CP-U-in-tooth` 中超过 `0.002 deg` 的样本。 |
| `numOrdinaryAngleMiss` | 22 | same tooth、fd healthy、non-collapse 的 ordinary miss。 |
| `ordinaryAngleMissRate` | 0.044 | ordinary miss 占全部 seed 的比例。 |
| `collapseHardCount` | 39 | hard-collapse 与 ordinary miss 是不同机制。 |
| `wrongToothCount` | 0 | controlled in-tooth 条件剥离了 wrong-tooth 污染。 |
| `fdNotHealthyCount` | 13 | fd-negative 不应混入 ordinary DoA gate。 |
| `easyHitCount` | 426 | gate 的 damage 风险主要来自这些 raw 已 hit 样本。 |
| `defaultOrdinaryMedianAngleDeg` | 0.0022747 | ordinary miss 多数只是略超阈值。 |
| `defaultOrdinaryP95AngleDeg` | 0.0036078 | ordinary tail 仍需要记录 P95，而不是只看 hit rate。 |
| `bestImplementableRescueRate` | 0.59091 | 当前可实现 gate 只能救约 13/22 ordinary miss。 |
| `bestImplementableDamageRate` | 0 | 推荐 policy 的优势是低误伤。 |

### 4.3 关键对比表

#### Ordinary candidate family aggregate

| candidate group | ordinary miss | hit rate | rescue rate | damage rate | median angle err (deg) | P95 angle err (deg) | median gain (deg) | median coherence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `default` | 22 | 0 | 0 | 0 | 0.0022747 | 0.0036078 | 0 | 1 |
| `wide-basin-entry` | 22 | 0.63636 | 0.63636 | 0.31818 | 0.00046751 | 0.0070408 | 0.00164 | 1 |
| `best-implementable` | 22 | 0.59091 | 0.59091 | 0 | 0.00063892 | 0.0027479 | 0.0015971 | 1 |
| `truth-doa-oracle` | 22 | 1 | 1 | 0 | 0.00040962 | 0.00076736 | 0.0018871 | 1 |

#### Recommended policy

| policy | ordinary rescue | ordinary hit | easy trigger | easy trigger count | easy damage | fd-negative damage | overall hit | P95 (deg) | P99 (deg) | max (deg) | recommendation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `wide-obj10-minStep0.001-maxStep0.004` | 0.59091 | 0.59091 | 0.023474 | 10 | 0 | 0 | 0.898 | 0.0037382 | 0.0051413 | 0.0054663 | `hold` |

与 disabled 对比：overall hit rate 从 `0.87` 到 `0.898`，median 从 `0.0004843 deg` 到 `0.00043519 deg`，但 P99 / max 没变；这说明它有小幅收益，但不具备推进 default 的强证据。

#### Wide candidate feature table

| miss type | num case | default hit | wide hit | wide rescue | wide damage | median obj gain | median angle gain (deg) | median wide angle (deg) | P95 wide angle (deg) | median default coherence | median wide coherence | median max abs DoA step (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `easy-hit` | 426 | 1 | 0.93427 | 0 | 0.077465 | 5.1487e-05 | 1.5644e-06 | 0.00037938 | 0.0039108 | 1 | 1 | 1.5123e-05 |
| `collapse-hard` | 39 | 0 | 0.64103 | 0 | 0.12821 | 8.836e+05 | 0.0030895 | 0.00076606 | 0.010781 | 0.10011 | 1 | 0.004049 |
| `fd-not-healthy-negative` | 13 | 0.69231 | 0.69231 | 0 | 0.23077 | 5.4092e+05 | 0.00093587 | 0.00032224 | 0.0090909 | 0.98268 | 1 | 0.002059 |
| `ordinary-angle-miss` | 22 | 0 | 0.63636 | 0.63636 | 0.31818 | 12.93 | 0.00164 | 0.00046751 | 0.0070408 | 1 | 1 | 0.0020616 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- 500-repeat 中 `numOrdinaryAngleMiss=22`，ordinary non-collapse miss 稳定存在，但占比只有 `4.4%`；它不应继续主导论文主仿真。
- ordinary miss 的 median coherence 为 `1`，说明 hard-collapse 的 `coherence-v1` gate 天然抓不到它。
- `wide-basin-entry` 对 ordinary miss 有明显救力，median angle 可接近 truth-DoA oracle。
- 推荐 policy 能保持 easy / fd-negative damage 为 0，说明存在低误伤 truth-free proxy。

### 5.2 仍未解决或反向的现象

- blanket wide 在 ordinary miss 上自身也有 `0.31818` damage，并在 easy-hit 上有 `0.077465` damage，不能作为默认路径。
- 推荐 policy 只救约 `59%` ordinary miss；P99 / max 基本不改善。
- 过宽 `obj0/minStep0` 会触发 `228` 个 easy 样本并出现 easy damage；过紧 policy 又迅速漏掉 ordinary miss。
- 单纯继续调 `objGain + DoA step` 二维阈值已经到头，需要 outcome classification 而不是继续阈值搜索。

### 5.3 代表性 case

| case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| ordinary miss group | non-collapse ordinary | `toothIdx=0`、fd healthy、coherence median `1`，但 angle 略超 `0.002 deg`。 | 说明 ordinary miss 与 hard-collapse 是不同机制。 |
| easy-hit group | damage risk | blanket wide 造成 `7.7465%` wide damage。 | 说明不能 blanket 打开 wide。 |
| recommended policy | hold | `easyDamage=0`、`fdNegativeDamage=0`，但 ordinary rescue 只有 `0.59091`。 | 支持 partial diagnostic，不支持 default。 |

## 6. 机制解释

### 6.1 当前解释

ordinary miss 不是 wrong-tooth，也不是 non-ref coherence collapse。它更像是 same-tooth 内的 threshold / bad-basin 子类型：raw `MS-MF-CP-U-in-tooth` 的中位表现很好，但少数样本在 DoA basin / adoption 上略超 `0.002 deg` 阈值。`wide-basin-entry` 能把一部分样本带到更好的 basin，但 wide center 对不同样本的方向并不一致，因此 blanket 或近 blanket policy 会误伤 easy 和 fd-negative 样本。

这说明 ordinary miss 不适合直接并入 hard-collapse 的 family-safe 分支。hard-collapse 的特征是 non-ref coherence floor 明显塌陷、objective gain 极大；ordinary miss 的 coherence 通常为 1，objective gain 也小得多，必须单独作为 outcome diagnostic 处理。

### 6.2 这个结果支持什么

- 支持把 ordinary-wide 作为 threshold / bad-basin diagnostic 保留。
- 支持在结果分析中区分 `wide-rescued / wide-damaged / wide-not-enough / gate-rejected`。
- 支持论文主线回到 CRB consistency、resolved RMSE、resolved/outlier rate，而不是追单一 fixed hit-rate。

### 6.3 这个结果不证明什么

- 不证明 ordinary-wide gate 可以进入 estimator 或 flow 默认路径。
- 不证明 combined hard-collapse + ordinary-wide policy 已经稳定。
- 不证明 full-flow CP/IP 或 MLE-vs-CRB 已经可用。
- 不证明继续调 `objGain` 和 DoA step 阈值能达到 `99%` hit rate。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；ordinary-wide 不下沉 estimator 主核。 |
| flow 默认路径 | 不推进；只作为 hold / partial diagnostic。 |
| regression | 不写；当前不是稳定 contract。 |
| replay / scan 下一步 | 若继续此线，新增 outcome classification、baseline-relative、CRB-normalized error；不继续盲调 gate。 |
| 论文图 / 论文口径 | diagnostic-only；可用来解释 full-sample outlier / threshold tail，不作为正文主精度目标。 |
| 排障记录 | 主记录已摘取“500-repeat 将 ordinary-wide 降级为 hold / partial diagnostic”。 |

## 8. 限制与禁止解释

- 不要把 `wide-obj10-minStep0.001-maxStep0.004` 写成 preferred candidate。
- 不要把 0 easy damage 与 0 fd-negative damage 单独当作晋级条件；coverage 不足同样重要。
- 不要用 fixed `angleHitThresholdDeg=0.002` 代替 RMSE / P95 / P99 / outlier / CRB-normalized 的联合判断。
- 不要把 ordinary miss 与 hard-collapse miss 合并成一个 gate。
- 不要把 truth label、truth tooth、truth DoA 或 truth `fdRef/fdRate` 放入 runtime selector、gate、candidate adoption 或 final winner。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-233931.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothOrdinaryWideGateDiagnose.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 100-repeat 中 ordinary miss 只有 4 个，`wide-obj10-minStep0.001-maxStep0.004` 曾显得像 preferred candidate。
- 500-repeat 后普通二维 gate 的局限变清楚：低误伤但 coverage 不足，因此当前结论降级为 `hold / partial diagnostic`。
