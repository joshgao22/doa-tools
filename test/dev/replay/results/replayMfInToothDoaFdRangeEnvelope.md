# replayMfInToothDoaFdRangeEnvelope 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / policy-discovery` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat` |
| 当前一句话结论 | base-seed scout + envelope 结果支持 `wide / single-MF` DoA basin-entry 与 `0.006 deg` family-safe step，但不支持默认引入非零 `fdRef / fdRate` joint bank。 |
| 决策影响 | 保留为 controlled in-tooth candidate envelope 与 policy 发现器；推进结论交给 `replayMfInToothGatedRescueEffectiveness` / flow-like replay 验证，不直接改 estimator 或 flow 默认路径。 |
| 下一步动作 | 当前 snapshot 已足够支撑 envelope 结论；后续若继续推进，应优先验证 gate / adoption，而不是扩大 fdRef/fdRate joint grid。 |
| 禁止误用 | 不能把 truth-centered half-tooth / truth center 放入 runtime；不能把 nonzero joint fd offset 下沉默认路径；不能把本 replay 的 policy sweep 直接写成 regression 或默认 flow。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothDoaFdRangeEnvelope.m`
- 结果文档：`test/dev/replay/results/replayMfInToothDoaFdRangeEnvelope.md`
- replay 类型：oracle-controlled / envelope 机制诊断 / policy-discovery replay。
- 主要问题：在 controlled in-tooth 条件下，hard same-tooth seed 的可实现好 candidate 来自 DoA basin-entry，还是需要非零 `fdRef / fdRate` joint offset。
- 观察范围：`baseSeed=253` 起的 100 个 scout seed；自动选择 12 个 hard / gate-miss / easy / fd-negative 代表 seed 做 heavy envelope surface、candidate coverage 与 policy sweep。
- 不覆盖范围：不做真实 subset tooth selection；不验证 full-flow；不证明 estimator default 已修复；不证明 CRB consistency；不改变 objective / residual / reference-sat 语义。
- truth 使用口径：truth 用于 half-tooth oracle fd box、truth-center oracle、seed 分类和结果评价；不进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `base-seed scout` | 从 `baseSeed` 开始做轻量 100-seed scout，先定位 hard、gate-miss、easy 与 fd-negative 代表样本。 | Evaluation only | 只改变本 replay 的离线选样。 | 用来减少手工 seed bias；不能变成 runtime seed selector。 |
| `selected envelope seed` | 从 scout 结果中挑出的 12 个 heavy surface seed。 | Evaluation only | 只改变 envelope 分析样本。 | 覆盖 hard collapse、gate miss、fd-negative 与 easy negative 控制样本。 |
| `truth-centered half-tooth` | 以 truth `fdRef` 所在 tooth 为中心，给 `fdRef` 一个 `0.49 / T_f` 半齿盒。 | Oracle only | 固定 in-tooth 观察条件。 | 用于剥离 wrong-tooth 影响；不能进入 runtime flow。 |
| `wide-coarse-doa-grid` | 在较宽 DoA box 内生成的可实现 DoA basin-entry center。 | No | 改变 DoA basin center。 | 当前 hard seed 的好 candidate 常来自它；但是否可 default 取决于后续 no-truth gate。 |
| `single-mf-coarse-doa-grid` | 用单星多帧 candidate 作为另一类可实现 DoA basin-entry center。 | No | 改变 DoA basin center。 | 与 wide family 互补，常用于救 non-ref coherence collapse hard seed。 |
| `truth` center | 以 truth DoA 作为 oracle center。 | Oracle only | 只提供上限。 | 可说明“如果 center 够好是否能救回”，不能作为 runtime candidate。 |
| `joint DoA-fdRef envelope` | 在 DoA offset 与 `fdRef` offset 上扫局部 objective surface。 | Evaluation only / Oracle context | 只改变离线 surface 扫描。 | 用来判断非零 `fdRef` step 是否提供额外收益；当前 evidence 不支持默认 joint fd bank。 |
| `joint DoA-fdRate envelope` | 在 DoA offset 与 `fdRate` offset 上扫局部 objective surface。 | Evaluation only / Oracle context | 只改变离线 surface 扫描。 | 多处边界命中，说明 fdRate offset 尚不具备可下沉的稳定范围。 |
| `family-safe-adopt` | 对 `wide / single-MF` basin-entry family 使用较宽 DoA step guard，而不是对全部 candidate blanket 放宽。 | No | 改变离线 policy adoption。 | 当前支持 `0.006 deg` family step；仍需 aggregate / flow-like 验证。 |
| `jointAddsOverPure` | best joint candidate 是否比 pure DoA candidate 多救 seed。 | Evaluation only | 只改结果评价。 | 当前为 false，说明 joint fd offset 没有相对 pure DoA center 形成额外救援。 |
| `jointUsesNonzeroFdRef` | best joint candidate 是否实际使用非零 `fdRef` step。 | Evaluation only | 只改结果评价。 | 当前为 false，是不推进默认 joint fd bank 的核心证据。 |
| `damage` | policy 或 candidate 使 easy / fd-negative 负样本变差。 | Evaluation only | 只改离线评价。 | damage 优先级高于单纯 hard rescue rate；本 snapshot 未观察到 easy / fd-negative damage。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat` | 2026-05-04 | representative | `seedMode=base-seed-search`，`baseSeed=253`，`numSearchRepeat=100`，`snrDb=10`，`oracleFdHalfToothFraction=0.49`，`maxEnvelopeSeedCount=12`，`angleHitThresholdDeg=0.002` | hard seed 的可实现好 candidate 主要来自 `wide / single-MF` DoA basin-entry；policy sweep 支持 `0.006 deg` family-safe step；不支持默认 nonzero `fdRef/fdRate` joint bank。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `seedMode = "base-seed-search"`
- `baseSeed = 253`
- `numSearchRepeat = 100`
- scout seed range：`253:352`
- envelope seed count：12
- selected envelope seeds：`[256, 277, 345, 312, 298, 283, 328, 319, 306, 346, 285, 341]`
- `oracleFdHalfToothFraction = 0.49`
- `oracleFdRateHalfWidthHzPerSec = 1000`
- `staticLocalDoaHalfWidthDeg = [0.002; 0.002]`
- `staticWideDoaHalfWidthDeg = [0.010; 0.010]`
- `truthLocalDoaHalfWidthDeg = [0.002; 0.002]`
- `probeDoaStepDegList = [-0.002; -0.001; 0; 0.001; 0.002]`
- `probeFdRefStepHzList = [-300; -150; 0; 150; 300]`
- `ridgeDoaOffsetDegList = -0.012:0.002:0.012 deg`
- `ridgeFdRefOffsetHzList = -360:60:360 Hz`
- `ridgeFdRateOffsetHzPerSecList = -240:80:240 Hz/s`
- `angleHitThresholdDeg = 0.002`
- gated rescue proxy：`coherenceThreshold=0.20`，`phaseResidThreshold=1.0 rad`
- checkpoint：scout pass `100/100`，envelope pass `12/12`；配置启用 checkpoint / resume / cleanup-on-success。
- snapshot 保存变量：`replayData`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| scout seed 数 | 100 | 用于自动发现代表 hard / negative case。 |
| heavy envelope seed 数 | 12 | 只对代表 seed 做 heavy surface，避免全量二维 / 三维扫描。 |
| hard-coherence-collapsed seed | 6 | 主要目标样本，用于检查 collapse hard case 是否存在可实现 basin-entry candidate。 |
| hard-gate-miss seed | 3 | 用于检查当前 coherence gate 漏掉的 hard 样本。 |
| fd-negative control seed | 2 | 用于检查 policy 是否误伤 fd 不健康负样本。 |
| easy-negative control seed | 1 | 用于检查 policy 是否误伤 easy 样本。 |
| `jointAddsOverPure` | false | joint DoA-fdRef 没有比 pure DoA center 多救 seed。 |
| `jointUsesNonzeroFdRef` | false | best joint candidate 的 `fdRef` step 仍为 0。 |
| best joint `fdRef` step | 0 Hz | 当前不支持默认引入 `fdRef` offset。 |
| `0.006 deg` family-safe damage | 0 | 在 selected envelope seeds 上未观察到 easy / fd-negative damage。 |

### 4.3 关键对比表

#### 自动选择的 envelope seed

| selection reason | seeds | 作用 |
|---|---|---|
| `hard-coherence-collapsed` | `256, 277, 345, 312, 298, 283` | 检查 coherence collapse hard case 是否存在可实现 basin-entry candidate。 |
| `hard-gate-miss` | `328, 319, 306` | 检查 hard 但 coherence gate 不触发的样本，评估是否需要新 gate proxy。 |
| `fd-not-healthy-negative` | `346, 285` | 检查 rescue / policy 是否误伤 fd 不健康负样本。 |
| `easy-negative` | `341` | 检查 rescue / policy 是否误伤 easy 样本。 |

#### Policy sweep：family DoA step guard

`policyRecommendationTable` 使用 cached candidate 做离线 policy sweep，不重跑 estimator。当前最关键的是固定 `coherenceThreshold=0.20`，比较 family DoA step guard `0.004 deg` 与 `0.006 deg`。

| policy | selectedAngleHitRate | hardAngleHitRate | overallP95AngleDeg | overallMaxAngleDeg | easyDamageRate | fdNegativeDamageRate | 备注 |
|---|---:|---:|---:|---:|---:|---:|---|
| family step `0.004 deg` | 0.50 | 0.55556 | 0.0048848 | 0.0054539 | 0 | 0 | strict guard 会拒绝部分有效 basin-entry candidate。 |
| family step `0.006 deg` | 0.58333 | 0.66667 | 0.0039798 | 0.0044191 | 0 | 0 | 当前 selected envelope seeds 上更优，进入 gated effectiveness 验证。 |

#### 代表性 hard seed candidate

| seed | best implementable center | 结果现象 | 对结论的作用 |
|---:|---|---|---|
| 256 | `single-mf-coarse-doa-grid` | angle error 可降到约 `0.00034936 deg`，non-ref coherence 可恢复到约 `0.99845`。 | 说明 single-MF basin-entry 能救 collapse hard seed。 |
| 277 | `wide-coarse-doa-grid` | angle error 可降到约 `3.6871e-05 deg`。 | 说明 wide basin-entry 是可实现好 center。 |
| 345 | `wide-coarse-doa-grid` | 与后续 `family-safe-adopt` 结果一致。 | 支持 wide family 与 `0.006 deg` guard。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- base-seed scout 自动挑出了 hard collapse、hard gate miss、fd-negative 和 easy-negative 四类样本，比固定手工 seed 更适合做 envelope policy 发现。
- hard seed 中确实存在可实现的好 candidate，且 best implementable center 主要来自 `wide-coarse-doa-grid` 或 `single-mf-coarse-doa-grid`。
- `jointAddsOverPure=false`、`jointUsesNonzeroFdRef=false`、best joint `fdRef` step 为 0，说明当前收益不是来自非零 `fdRef` joint offset。
- `0.006 deg` family-safe step 相比 `0.004 deg` 同时提升 selected / hard angle hit，并降低 P95 / max；selected envelope seeds 上 easy / fd-negative damage 均为 0。
- 两阶段 checkpoint / scout + envelope 设计有效：先用 100 个 scout seed 找代表 case，再只对 12 个 seed 做 heavy surface。

### 5.2 仍未解决或反向的现象

- DoA-fdRate surface 多处边界命中，说明 fdRate offset 方向尚未形成稳定可用的搜索盒。
- `hard-gate-miss` seed 仍存在，说明单纯 coherence collapse proxy 不能覆盖全部 hard case。
- 本 replay 的 negative controls 很少：fd-negative 2 个、easy-negative 1 个；因此“无 damage”只能作为进入下一步验证的条件，不能当作默认安全证明。
- envelope 只在 truth-centered half-tooth 中进行，不能说明 wrong-tooth / subset selection 已解决。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 256 | hard-coherence-collapsed | `single-mf-coarse-doa-grid` 可把 angle error 降到约 `0.00034936 deg`，coherence 恢复到约 `0.99845`。 | 支持 single-MF basin-entry family。 |
| 277 | hard-coherence-collapsed | `wide-coarse-doa-grid` 可把 angle error 降到约 `3.6871e-05 deg`。 | 支持 wide basin-entry family。 |
| 345 | hard-coherence-collapsed | wide family 与后续 family-safe adoption 一致。 | 支持把 `0.006 deg` family step 送入 gated effectiveness。 |
| 328 / 319 / 306 | hard-gate-miss | hard 但当前 gate proxy 不一定触发。 | 提醒后续 gate 不能只看 collapse coherence。 |
| 346 / 285 | fd-negative control | 用于观察 policy 是否误伤 fd 不健康样本。 | 当前 snapshot 未观察到 damage，但样本数有限。 |
| 341 | easy-negative control | 用于观察 policy 是否误伤 easy 样本。 | 当前 snapshot 未观察到 damage，但不能替代 flow-like 验证。 |

## 6. 机制解释

### 6.1 当前解释

这轮 envelope replay 的关键价值是把 same-tooth hard case 的候选来源拆开。`replayMfCombToothSurface` 和 ridge replay 已提示 DoA-Doppler coupling 存在，但 coupling 存在并不等于应该把非零 `fdRef / fdRate` offset 下沉默认路径。本 replay 在 truth-centered half-tooth 条件下直接比较 implementable center、joint DoA-fdRef / DoA-fdRate envelope 和离线 policy，结果显示：可实现收益主要来自 **换一个更好的 DoA basin-entry center**，而不是给最终点附近加一个非零 `fdRef` 或 `fdRate` offset。

hard collapse seed 能被 `wide / single-MF` center 救回，说明当前 same-tooth tail 更像 DoA / local-state basin-entry 问题；非参考星 coherence collapse 是症状之一。相反，joint fd 方向没有稳定贡献：best joint `fdRef` step 仍为 0，fdRate surface 还存在边界命中。因此，继续堆 joint fd bank 的维护成本和误伤风险都高于当前证据能支持的收益。

policy sweep 的作用是把“candidate 存在”转换成“哪个 adoption guard 值得下一步验证”。`0.006 deg` family-safe step 只对 `wide / single-MF` basin-entry family 放宽，而不是 blanket 放宽全部 candidate，这与当前机制解释一致：hard seed 需要 basin-entry center，但 easy / fd-negative control 不能被随意改坏。

### 6.2 这个结果支持什么

- 支持 `wide / single-MF` DoA basin-entry 是当前 controlled in-tooth hard-collapse 的主要可实现 candidate family。
- 支持 `familySafeAdoptMaxAbsDoaStepDeg=0.006` 作为下一步 `replayMfInToothGatedRescueEffectiveness` 的 policy 候选。
- 支持继续把 joint DoA-fdRef / DoA-fdRate surface 保留为机制诊断，而不是默认救援 bank。
- 支持 base-seed scout + 自动 seed selection 作为 envelope replay 的默认入口，减少手工 seed bias。

### 6.3 这个结果不证明什么

- 不证明 estimator 默认路径已修复。
- 不证明 full-flow 或真实 subset-periodic flow 已经安全。
- 不证明可以写 regression 契约。
- 不证明可以把 truth-centered half-tooth、truth center 或 oracle label 放入 runtime。
- 不证明 `0.006 deg` family-safe step 在大 MC 或 flow-like 条件下无 damage。
- 不证明 nonzero `fdRef/fdRate` joint bank 有默认收益；当前证据恰好不支持默认打开它。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉到 estimator 主核、objective、residual 或 reference-sat 语义。 |
| flow 默认路径 | 不改；本 replay 只发现 candidate / policy，是否进入 flow 需由 gated effectiveness 和 flow-like replay 决定。 |
| regression | 不写；这是 oracle-controlled envelope 和 offline policy sweep，不是稳定 runtime 契约。 |
| replay / scan 下一步 | 已推进到 `replayMfInToothGatedRescueEffectiveness` 验证 aggregate 表现；若继续研究，应优先拆 gate-miss 与 damage，而不是扩 fd grid。 |
| 论文图 / 论文口径 | diagnostic-only / appendix 机制材料候选；不能作为 paper-facing 精度图。论文主线仍应优先 CRB consistency、resolved-regime RMSE、CP/IP controlled trade-off。 |
| 排障记录 | 可保留为机制证据：same-tooth tail 的主要收益来自 DoA basin-entry center，joint fd bank 暂不推进默认路径。 |

## 8. 限制与禁止解释

- 不要把 truth-centered half-tooth oracle 当作 runtime tooth selector。
- 不要把 truth DoA center 或 truth-oracle candidate 放入 runtime candidate bank。
- 不要把 `policyRecommendationTable` 的离线 sweep 直接写成 regression 或默认 flow。
- 不要因为 selected envelope seeds 上无 easy / fd-negative damage，就认为大 MC 或 flow-like gate 已经安全。
- 不要把 DoA-Doppler coupling 解释为必须默认引入 nonzero `fdRef/fdRate` joint offset；当前 candidate coverage 不支持这个结论。
- 不要继续扩大 fdRef/fdRate joint grid 作为第一反应；当前更值得处理的是 gate / adoption / flow-like damage。
- 不要把本 replay 的 angle hit rate 当作论文主指标；它只是 hard/tail candidate 发现指标。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
test/dev/replay/replayMfInToothDoaFdRangeEnvelope.m
```

只运行 `Summary output and plotting` 小节，即可重出 compact table、candidate coverage、policy recommendation 和 envelope 图。

## 10. 历史备注

- 本 snapshot 是当前唯一代表性 envelope 结果。
- 它承接 `replayMfInToothTailCaseDiagnose` 暴露的 same-tooth hard-collapse 问题，并把候选来源进一步收缩到 `wide / single-MF` DoA basin-entry family。
- 它也解释了为什么 `replayMfInToothGatedRescueEffectiveness` 应先验证 family-safe `0.006 deg` step，而不是默认打开 nonzero `fdRef/fdRate` joint bank。
- 若后续 flow-like replay 发现 easy damage，应回到 gate / adoption 拆分，而不是直接否定本 replay 中“candidate family 有救力”的结论。
