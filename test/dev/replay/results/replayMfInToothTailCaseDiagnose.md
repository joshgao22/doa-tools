# replayMfInToothTailCaseDiagnose 结果记录

## 对应 replay

- 脚本：`test/dev/replay/replayMfInToothTailCaseDiagnose.m`
- 目标：固定 `replayMfInToothFdRangeOracle` 暴露出的 in-tooth tail seed，定位 same-tooth hard case 是否来自 non-ref coherence collapse、fd/rate 不健康，或 DoA/local-state basin，并验证 gated rescue bank 是否只救 hard case、不误伤负样本。

## Snapshot index

| snapshot / log | 状态 | 配置摘要 | 备注 |
|---|---|---|---|
| `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat` | representative | `snrDb=10`，`contextBaseSeed=253`，`seedList=[277 283 298 256 293 280 268 284]`，`fdOracleHalfToothFraction=0.49`，`fdRateOracleHalfWidthHzPerSec=1000`，gated rescue coherence threshold `0.20`，phase threshold `1.00 rad` | 当前保留的 in-tooth tail diagnose / gated rescue bank 代表性结果；`saveSnapshot=1`，只保存轻量 `replayData` |

## 2026-04-27 10:20 结果分析

### 配置与总体状态

- `contextBaseSeed=253`，与源 half-tooth oracle replay 对齐。
- tail / negative-check seeds：`277, 283, 298, 256, 293, 280, 268, 284`。
- 所有 method 的 `toothHitRate=1`，说明本 replay 仍处于 in-tooth oracle 范围内；当前问题不是 wrong-tooth。
- `MS-MF-CP-U-in-tooth` 的 aggregate：angle RMSE `0.0032273 deg`，median `0.0027772 deg`，non-ref coherence floor median `0.50663`。
- `MS-MF-CP-U-truth-DoA-oracle` 的 aggregate：angle RMSE `0.00036907 deg`，median `0.00033505 deg`，non-ref coherence floor median `1`。
- 当前一句话结论：默认 MS in-tooth 路径仍被少数 same-tooth non-ref coherence collapse tail 拉坏；但 collapse-gated `wide + single-MF` rescue bank 在本组 seeds 上能全救 hard case，并避免 easy / fd-not-healthy 负样本误伤。本轮已保存 snapshot，可从 `replayData` 恢复表格与图。

### Tail classification

| seed | tail class | default MS angle (deg) | truth-DoA oracle angle (deg) | default non-ref coherence floor | default non-ref RMS phase residual (rad) | 备注 |
|---:|---|---:|---:|---:|---:|---|
| 277 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0047349 | 0.0000325 | 0.089376 | 1.7993 | wide center 可救 |
| 283 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0036180 | 0.0005192 | 0.099451 | 1.6232 | wide center 可救 |
| 298 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0035943 | 0.0001577 | 0.10056 | 1.6218 | single-MF center 可救 |
| 256 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0054539 | 0.0000215 | 0.000660 | 2.2202 | single-MF center 可救 |
| 293 | `same-tooth light/unclear` | 0.0006742 | 0.0006498 | 1.0000 | 0.000913 | easy negative |
| 280 | `same-tooth light/unclear` | 0.0002487 | 0.0003634 | 1.0000 | 0.001482 | easy negative |
| 268 | `same-tooth + fd not healthy` | 0.0019601 | 0.0003067 | 0.9127 | 0.4235 | fd-not-healthy negative |
| 284 | `same-tooth + fd not healthy` | 0.0008910 | 0.0003812 | 0.9938 | 0.1112 | fd-not-healthy negative；blanket wide 会误伤 |

### Candidate / line probe 观察

- final-centered small DoA 或 DoA+fdRef probe 仍不能救 `277 / 283 / 298 / 256`，说明不能把 final-centered very-small polish 接进默认 flow。
- 加密 line probe 显示 good basin 很窄：hard seeds 通常需要沿 default/static 到 truth DoA 方向移动到 `alpha=0.95~0.98` 才恢复 coherence。
- 不含 truth 的 implementable center sweep 给出可实现入口：
  - `277 / 283` 由 `wide-coarse-doa-grid` 救回；
  - `298 / 256` 由 `single-mf-coarse-doa-grid` 救回。
- `284` 是关键负样本：wide candidate objective 和 coherence 都很强，但 angle 会从 `0.000891 deg` 被拉坏到 `0.010301 deg`；因此不能只按 objective gain 或 coherence gain blanket adoption。

### Rescue bank aggregate

| bank | trigger rate | hard trigger rate | easy trigger rate | fd-negative trigger rate | hard rescue rate | hard median selected angle (deg) | easy damage rate | fd-negative damage rate | overall max selected angle (deg) | 结论 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `disabled` | 0 | 0 | 0 | 0 | 0 | 0.0041764 | 0 | 0 | 0.0054539 | baseline，不救 hard case |
| `wide-centered-coarse` | 0 | 0 | 0 | 0 | 0.5 | 0.0018424 | 0 | 0.5 | 0.010301 | 能救 277/283，但误伤 284 |
| `single-mf-centered-coarse` | 0 | 0 | 0 | 0 | 0.5 | 0.0042682 | 1 | 1 | 0.0080813 | 能救 298/256，但 blanket 使用误伤严重 |
| `wide-single-bank` | 0 | 0 | 0 | 0 | 1 | 0.0004135 | 0 | 0.5 | 0.010301 | hard case 全救，但仍误伤 284 |
| `gated-wide-single-bank` | 0.5 | 1 | 0 | 0 | 1 | 0.0004135 | 0 | 0 | 0.0019601 | 当前最稳：hard 全救，负样本不误伤 |

### Gated rescue per-seed 行为

| seed | gate triggered | gate reason | selected family | selected angle (deg) | selected coherence floor | 是否 hard rescued / damaged |
|---:|---|---|---|---:|---:|---|
| 277 | true | `non-ref-coherence-collapse` | `wide-coarse-doa-grid` | 0.0000369 | 1.0000 | rescued，不 damaged |
| 283 | true | `non-ref-coherence-collapse` | `wide-coarse-doa-grid` | 0.0005268 | 1.0000 | rescued，不 damaged |
| 298 | true | `non-ref-coherence-collapse` | `single-mf-coarse-doa-grid` | 0.0004776 | 0.9998 | rescued，不 damaged |
| 256 | true | `non-ref-coherence-collapse` | `single-mf-coarse-doa-grid` | 0.0003494 | 0.9985 | rescued，不 damaged |
| 293 | false | `coherence-not-collapsed` | `default-final` | 0.0006742 | 1.0000 | easy 不触发，不 damaged |
| 280 | false | `coherence-not-collapsed` | `default-final` | 0.0002487 | 1.0000 | easy 不触发，不 damaged |
| 268 | false | `coherence-not-collapsed` | `default-final` | 0.0019601 | 0.9127 | fd-negative 不触发，不 damaged |
| 284 | false | `coherence-not-collapsed` | `default-final` | 0.0008910 | 0.9938 | fd-negative 不触发，成功挡住 wide 误伤 |

### No-truth gate 可实现性说明

- 当前 `gated-wide-single-bank` 的触发条件只使用默认 `MS-MF-CP-U-in-tooth` 估计结果回代得到的 non-ref coherence floor 与 non-ref RMS phase residual。它们来自接收信号、已知训练块/阵列模型以及卫星轨道几何下的 profile residual / phase-consistency 评价，不使用 truth DoA、truth `fdRef/fdRate`、`tailClass` 或人工 seed 标签。
- `wide-coarse-doa-grid` 与 `single-mf-coarse-doa-grid` 的候选中心也应理解为 data-derived center：前者来自同一批接收信号上的 wide-DoA 多星候选，后者来自单星多帧估计结果。二者都不是由 truth oracle 给出。
- replay 中的 `truth-doa-fixed`、`truth-doa-oracle`、default/static-to-truth line probe 和 `tailClass` 只用于机制定位与结果标注，不能进入实际系统的 trigger、candidate 生成或 adoption 逻辑。
- 因此，按当前设计，gated rescue 符合“实际估计只能使用接收信号、已知训练结构、阵列/卫星轨道参数”的要求；后续若进入 flow 层实现，应继续保持这个边界，并在代码中避免读取 truth 字段或 oracle table 参与 gate。

### 图形口径

- 旧图把不同 seed 之间用折线连接，容易让读者误解为连续曲线；当前不再保留该画法。
- 默认图改为 `gated-wide-single-bank` 的前后对比：按 seed 分组显示 default MS-MF angle error 与 gated selected angle error、angle gain、default / selected non-ref coherence floor。
- 若恢复的是旧 snapshot，且没有 `rescueBankSummaryTable`，脚本才退回 method-level grouped bar fallback；fallback 也不再用 seed 间连线。

### 可观察现象

- same-tooth hard case 的共同特征是：`toothIdx=0`、`fdRef/fdRate` 误差小，但 non-ref coherence floor 极低，且 truth-DoA oracle 能把 coherence 恢复到 1。
- `CP-K` 与 `CP-U` 在 hard seeds 上几乎同坏，说明当前 tail 不是 unknown-rate release 的主问题。
- `truth-doa-fixed` 能恢复 coherence，说明 phase chain / evaluator 上限存在；默认路径的问题是没有进入好 DoA/local-state basin。
- `wide` 与 `single-MF` center 分别能救不同 hard seed：后续候选 bank 不能只保留其中一个。
- blanket `wide-single-bank` 不安全：它会采用 284 的 wide candidate 并造成明显角度误伤。
- collapse gate 当前非常干净：`hardTriggerRate=1`，`easyTriggerRate=0`，`fdNegativeTriggerRate=0`。

## 当前结论

- `gated-wide-single-bank` 是目前最有希望的 same-tooth hard-case rescue 候选。
- 当前 gate 口径 `nonRefCoherenceFloor < 0.20` 且 `nonRefRmsPhaseResidRad >= 1.0` 在本组 seeds 上能准确触发 4 个 hard-collapse seed，并避开 4 个负样本。
- 该结果支持继续推进 flow-like replay 验证，但还不支持直接修改 estimator 默认路径。
- 不建议把 `wide-centered`、`single-MF-centered` 或二者组合做成 blanket 常驻；必须先由 collapse gate 限制触发范围。

## 可能的修改方向

### 推荐继续验证

1. 新增或扩展一个更接近 `runSimpleDynamicSubsetPeriodicFlow` 的 gated rescue replay。
   - 比较 `disabled`、`gated-wide-only`、`gated-single-mf-only`、`gated-wide-single-bank`。
   - 小 MC 建议先用 `20~40` repeats，不直接跑大 MC。
   - 重点统计 `triggerRate`、`hardRescueRate`、`easyDamageRate`、`fdNegativeDamageRate`、angle median / p95 / max 与 selected candidate family。
2. 将 collapse gate 保持为 truth-free 指标。
   - 可先使用 `defaultNonRefCoherenceFloor < 0.20` 与 `defaultNonRefRmsPhaseResidRad >= 1.0`。
   - 后续若触发不足，可再测试阈值，不要引入 truth-aware gate。
3. 若 flow-like replay 成立，再考虑在 flow 层加入 conditional basin-entry bank。
   - 放置位置应在 flow orchestration / final same-tooth rescue 层，不下沉到 estimator 主核。
   - 不改变 `replayMfInToothFdRangeOracle` 的 oracle 上限职责。

### 暂不建议

- 不把 final-centered very-small DoA polish 或 final-centered DoA+fdRef joint probe 接进 flow；当前 hard seed 已证明它们救不动。
- 不做 blanket wide DoA；284 已经给出明确误伤。
- 不做 blanket single-MF coarse；293/280/268/284 的误伤风险更高。
- 不把 truth-line 或 truth-DoA oracle 结果用于 selector / adoption；这些只作为 replay 机制定位。
- 不在 estimator 默认主路径中直接加 rescue bank；先在 flow-like replay 中验证 no-truth gate 与负样本误伤。

## 与 `runSimpleDynamicSubsetPeriodicFlow` 的关系

`runSimpleDynamicSubsetPeriodicFlow` 的验证不应放进 `replayMfInToothFdRangeOracle`。

原因：

- `replayMfInToothFdRangeOracle` 是 truth-centered half-tooth oracle，用于回答“正确 tooth 内 estimator 上限是否足够好”，它绕开了真实 flow 的 subset tooth selection 和 no-truth adoption。
- `runSimpleDynamicSubsetPeriodicFlow` 要验证的是真实 flow 下的 subset 选齿、periodic refine、conditional rescue / adoption 是否能在 no-truth 条件下接住好盆地。
- 如果把 flow 验证塞进 oracle replay，会混淆“oracle 上限”和“可实现 flow”的证据，导致表格解释变得不清楚。

更合理的路线是：

1. `replayMfInToothFdRangeOracle` 继续只做 half-tooth oracle 上限和 tail seed 暴露；
2. `replayMfInToothTailCaseDiagnose` 继续只做 fixed tail 的机制定位和 gated rescue bank 预验证；
3. 新增一个 flow-like gated rescue replay，或在 `replayMfPeriodicVsSubsetToothSelect` 后续阶段中增加 gated rescue 对照，用于验证真实 flow 可实现性；
4. flow-like replay 通过后，再决定是否进入 `runSimpleDynamicSubsetPeriodicFlow` 默认 conditional rescue。

## 对主流程的影响

- 当前结果把 same-tooth tail 的候选修复方向从“泛化 joint refine”进一步收敛为 **collapse-gated wide+single-MF basin-entry bank**。
- `runSimpleDynamicSubsetPeriodicFlow` 后续应重点观察 gate 是否只在 non-ref coherence collapse 时触发，以及 rescue 后是否保持 selected tooth 不变。
- 如果后续小 MC 中 gated bank 仍满足 hard case 改善且 easy / fd-negative 不误伤，可以考虑把它升级为 flow 层 conditional rescue 候选。
- 在进入 regression 前，还需要至少一轮 flow-like replay 证明 no-truth gate 稳定；当前固定 8 seed 结果不宜直接固化成 regression。

## 恢复方式

本次 snapshot 已保存为 `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat`。恢复时建议使用：

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后运行 `replayMfInToothTailCaseDiagnose.m` 的 `Summary output and plotting` section。
