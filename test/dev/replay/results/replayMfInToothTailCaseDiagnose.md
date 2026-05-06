# replayMfInToothTailCaseDiagnose 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `candidate mechanism / fixed-tail diagnostic` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat` |
| 当前一句话结论 | 固定 half-tooth oracle 暴露出的 8 个 tail / negative seeds 后，collapse-gated `wide + single-MF` basin-entry bank 能全救 hard-collapse seeds，并避开 easy 与 fd-not-healthy 负样本。 |
| 决策影响 | 支持把候选方向收敛到 flow 层 conditional same-tooth basin-entry rescue；仍不能直接进入 estimator 主核或 regression。 |
| 下一步动作 | 进入 flow-like replay，验证在真实 subset-periodic flow 中 gate 是否仍只触发 hard case、是否保持 selected tooth 不变、是否无 easy / fd-negative damage。 |
| 禁止误用 | 不能把固定 8 seed 结果当成 MC 统计；不能把 tail class、truth DoA 或 truth fd 放入 runtime gate；不能 blanket 打开 wide / single-MF。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothTailCaseDiagnose.m`
- 结果文档：`test/dev/replay/results/replayMfInToothTailCaseDiagnose.md`
- replay 类型：oracle-controlled fixed-tail mechanism replay。
- 主要问题：`replayMfInToothFdRangeOracle` 暴露的 same-tooth tail 是否由 non-ref coherence collapse 主导，以及 data-derived `wide / single-MF` basin-entry center 能否在 no-truth gate 下救回 hard case。
- 观察范围：固定 seeds `277, 283, 298, 256, 293, 280, 268, 284`；仍处于 truth-centered half-tooth oracle 条件；比较 default、truth-DoA oracle、candidate probes、line probes 与 rescue bank。
- 不覆盖范围：不验证真实 subset tooth selection；不验证 full-flow adoption；不证明 estimator default 已修复；不形成 regression contract。
- truth 使用口径：truth 用于 tail classification、truth-DoA oracle、default-to-truth line probe 与结果评价；gated rescue 的触发条件与候选中心不得使用 truth。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `same-tooth + fd healthy + non-ref coherence collapsed` | `toothIdx=0` 且 `fdRef/fdRate` 健康，但非参考星 coherence floor 很低、phase residual 大。 | Evaluation only | 只改离线分类。 | 本 replay 的 hard case；说明问题是 DoA/local-state basin，不是 wrong-tooth。 |
| `fd-not-healthy negative` | 频率链不健康或不应被 DoA rescue 接管的负样本。 | Evaluation only | 只改离线分类。 | 用来测试 gate 是否误伤频率问题样本。 |
| `wide-coarse-doa-grid` | 从 wide DoA 多星候选中心重新进入 basin。 | No | 改 DoA basin center。 | 能救部分 hard seed，但 blanket 使用会误伤 seed 284。 |
| `single-mf-coarse-doa-grid` | 从单星多帧估计中心重新进入 basin。 | No | 改 DoA basin center。 | 能救另一部分 hard seed，但 blanket 使用误伤更明显。 |
| `wide-single-bank` | 同时评估 wide 与 single-MF basin-entry family。 | No | 扩大 candidate bank。 | hard seeds 全救，但没有 gate 时仍可能采用坏 wide candidate。 |
| `gated-wide-single-bank` | 只在 default non-ref coherence collapse 且 phase residual 足够大时启用 wide+single-MF bank。 | No | 改 rescue 触发与 candidate adoption。 | 当前 fixed-tail 最稳候选，但仍需 flow-like 小 MC 验证。 |
| `truth-DoA oracle` | DoA 放在 truth 附近的离线上限。 | Oracle only | 只改上限评价。 | 说明 hard seed 可达，不是模型上限不足；不能进入 runtime。 |
| `default/static-to-truth line probe` | 从 default / static center 沿 truth DoA 方向做离线线扫。 | Oracle only | 只改诊断。 | 用于说明好 basin 很窄；不能当作可实现 selector。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat` | 2026-04-27 | representative | `snrDb=10`，`contextBaseSeed=253`，`seedList=[277 283 298 256 293 280 268 284]`，`fdOracleHalfToothFraction=0.49`，`fdRateOracleHalfWidthHzPerSec=1000`，gate coherence threshold `0.20`，phase threshold `1.0 rad`。 | `gated-wide-single-bank` 触发 4 个 hard-collapse seed，hard rescue rate `1`，easy / fd-negative damage 为 0。 | 当前唯一代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `contextBaseSeed = 253`
- `numRepeat = 8`
- seed list：`277, 283, 298, 256, 293, 280, 268, 284`
- `fdRef` oracle 范围：truth-centered `0.49` tooth
- `fdRate` oracle 范围：truth-centered `±1000 Hz/s`
- fd healthy thresholds：`abs(fdRefErr) <= 1 Hz`，`abs(fdRateErr) <= 50 Hz/s`
- hard collapse thresholds：`coherenceCollapseThreshold = 0.95`，`truthDoaGapThresholdDeg = 0.001`
- gated rescue thresholds：`defaultNonRefCoherenceFloor < 0.20` 且 `defaultNonRefRmsPhaseResidRad >= 1.0`
- snapshot 保存变量：`replayData`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `numSeed` | 8 | 固定 tail / negative seeds，不是大 MC。 |
| `all method toothHitRate` | 1 | 当前问题已剥离 wrong-tooth。 |
| `default MS angle RMSE` | `0.0032273 deg` | same-tooth hard tail 明显。 |
| `default MS angle median` | `0.0027772 deg` | 8 seed 中 default 多数较差。 |
| `default non-ref coherence median` | `0.50663` | hard seed 存在 non-ref coherence collapse。 |
| `truth-DoA oracle RMSE` | `0.00036907 deg` | 上限存在，hard seed 可被救回。 |
| `gated-wide-single-bank hard rescue rate` | 1 | 4 个 hard seed 全救。 |
| `gated-wide-single-bank easy damage` | 0 | easy negative 未误伤。 |
| `gated-wide-single-bank fd-negative damage` | 0 | 频率负样本未误伤。 |

### 4.3 关键对比表

#### Tail classification

| seed | tail class | default MS angle (deg) | truth-DoA oracle angle (deg) | default non-ref coherence floor | default non-ref RMS phase residual (rad) | 备注 |
|---:|---|---:|---:|---:|---:|---|
| 277 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0047349 | 0.0000325 | 0.089376 | 1.7993 | wide center 可救。 |
| 283 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0036180 | 0.0005192 | 0.099451 | 1.6232 | wide center 可救。 |
| 298 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0035943 | 0.0001577 | 0.10056 | 1.6218 | single-MF center 可救。 |
| 256 | `same-tooth + fd healthy + non-ref coherence collapsed` | 0.0054539 | 0.0000215 | 0.000660 | 2.2202 | single-MF center 可救。 |
| 293 | `same-tooth light/unclear` | 0.0006742 | 0.0006498 | 1.0000 | 0.000913 | easy negative。 |
| 280 | `same-tooth light/unclear` | 0.0002487 | 0.0003634 | 1.0000 | 0.001482 | easy negative。 |
| 268 | `same-tooth + fd not healthy` | 0.0019601 | 0.0003067 | 0.9127 | 0.4235 | fd-not-healthy negative。 |
| 284 | `same-tooth + fd not healthy` | 0.0008910 | 0.0003812 | 0.9938 | 0.1112 | fd-not-healthy negative；blanket wide 会误伤。 |

#### Rescue bank aggregate

| bank | trigger rate | hard trigger rate | easy trigger rate | fd-negative trigger rate | hard rescue rate | hard median selected angle (deg) | easy damage | fd-negative damage | overall max selected angle (deg) | 结论 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `disabled` | 0 | 0 | 0 | 0 | 0 | 0.0041764 | 0 | 0 | 0.0054539 | baseline，不救 hard。 |
| `wide-centered-coarse` | 0 | 0 | 0 | 0 | 0.5 | 0.0018424 | 0 | 0.5 | 0.010301 | 能救 277/283，但误伤 284。 |
| `single-mf-centered-coarse` | 0 | 0 | 0 | 0 | 0.5 | 0.0042682 | 1 | 1 | 0.0080813 | 能救 298/256，但 blanket 误伤严重。 |
| `wide-single-bank` | 0 | 0 | 0 | 0 | 1 | 0.0004135 | 0 | 0.5 | 0.010301 | hard 全救，但仍误伤 284。 |
| `gated-wide-single-bank` | 0.5 | 1 | 0 | 0 | 1 | 0.0004135 | 0 | 0 | 0.0019601 | 当前最稳：hard 全救、负样本不误伤。 |

#### Gated rescue per-seed behavior

| seed | gate triggered | gate reason | selected family | selected angle (deg) | selected coherence floor | 评价 |
|---:|---|---|---|---:|---:|---|
| 277 | true | `non-ref-coherence-collapse` | `wide-coarse-doa-grid` | 0.0000369 | 1.0000 | rescued，不 damaged。 |
| 283 | true | `non-ref-coherence-collapse` | `wide-coarse-doa-grid` | 0.0005268 | 1.0000 | rescued，不 damaged。 |
| 298 | true | `non-ref-coherence-collapse` | `single-mf-coarse-doa-grid` | 0.0004776 | 0.9998 | rescued，不 damaged。 |
| 256 | true | `non-ref-coherence-collapse` | `single-mf-coarse-doa-grid` | 0.0003494 | 0.9985 | rescued，不 damaged。 |
| 293 | false | `coherence-not-collapsed` | `default-final` | 0.0006742 | 1.0000 | easy 不触发。 |
| 280 | false | `coherence-not-collapsed` | `default-final` | 0.0002487 | 1.0000 | easy 不触发。 |
| 268 | false | `coherence-not-collapsed` | `default-final` | 0.0019601 | 0.9127 | fd-negative 不触发。 |
| 284 | false | `coherence-not-collapsed` | `default-final` | 0.0008910 | 0.9938 | 成功挡住 wide 误伤。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- 4 个 hard seed 同时满足 same-tooth、fd healthy、non-ref coherence collapsed；这把 tail 从 wrong-tooth / fd failure 中分离出来。
- `truth-DoA oracle` 能把 hard seeds 的 angle 拉到 `1e-4 deg` 量级并恢复 coherence，说明 evaluator / phase chain 上限存在。
- `wide` 与 `single-MF` center 分别救不同 hard seed，单独保留任一 family 都不够。
- blanket `wide-single-bank` 能全救 hard，但会误伤 seed 284；必须用 no-truth collapse gate 限制触发。
- `gated-wide-single-bank` 的 `hardTriggerRate=1`、`easyTriggerRate=0`、`fdNegativeTriggerRate=0`，是当前最强正证据。

### 5.2 仍未解决或反向的现象

- 本 replay 只有固定 8 seed，不是稳定 MC；它证明机制可行，但不能证明 gate 的总体误伤率。
- 触发条件目前在 oracle in-tooth 条件下验证，尚未和真实 subset tooth selection、periodic refine、flow warning 交互。
- final-centered small DoA 或 DoA+fdRef probe 救不动 hard seeds，因此不能把 final-centered polish 作为默认修复。

### 5.3 代表性 seed / case

| seed | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 277 / 283 | hard-collapse | wide center 可救，selected angle 约 `3.69e-5 / 5.27e-4 deg`。 | 说明 wide basin-entry family 必须保留。 |
| 298 / 256 | hard-collapse | single-MF center 可救，selected angle约 `4.78e-4 / 3.49e-4 deg`。 | 说明 single-MF basin-entry family 必须保留。 |
| 284 | fd-negative damage sentinel | wide candidate objective / coherence 强，但 angle 可被拉坏到约 `0.010301 deg`。 | 说明不能 blanket adoption。 |
| 293 / 280 | easy negative | gate 不触发，保持 default。 | 说明当前 no-truth gate 对 easy 样本干净。 |

## 6. 机制解释

### 6.1 当前解释

same-tooth hard seed 的频率链已经足够健康，但 default DoA / local-state basin 没有同时接住非参考星。表现上是 reference link 可以拟合，但 non-ref coherence floor 很低、phase residual 很大。truth-DoA oracle 能恢复 coherence，说明不是模型上限问题，而是 basin-entry center 问题。

`wide` 与 `single-MF` center 是两类不同的 data-derived basin-entry：前者来自宽 DoA 多星候选，后者来自单星多帧中心。它们能救的 hard seed 不完全相同，所以当前候选需要是 family bank，而不是单一路线。但 candidate objective / coherence 本身仍可能在负样本上很诱人，因此必须由 default non-ref coherence collapse gate 限制触发。

### 6.2 这个结果支持什么

- 支持将修复方向从 generic joint refine 收敛到 collapse-gated `wide + single-MF` basin-entry bank。
- 支持把 no-truth non-ref coherence / phase residual 作为 flow-like gate 的首版 proxy。
- 支持暂缓 final-centered very-small polish 和 blanket bank。

### 6.3 这个结果不证明什么

- 不证明 full-flow 已通过；本 replay 绕开了真实 subset tooth selection。
- 不证明 `gated-wide-single-bank` 可以直接进入 estimator 默认路径。
- 不证明固定阈值 `0.20 / 1.0 rad` 对大 MC 或不同 SNR 都稳定。
- 不证明 tail label、truth DoA 或 truth line probe 可以进入 runtime gate。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；rescue bank 不下沉 estimator 主核。 |
| flow 默认路径 | 暂不改；推进到 flow-like replay 或 `runSimpleDynamicSubsetPeriodicFlow` final same-tooth rescue 层验证。 |
| regression | 不写；固定 8 seed 机制证据不是稳定契约。 |
| replay / scan 下一步 | 比较 disabled / gated-wide-only / gated-single-MF-only / gated-wide-single-bank；小 MC 先用 20–40 repeats。 |
| 论文图 / 论文口径 | diagnostic-only；可解释 same-tooth outlier 机制，不作为正文主性能图。 |
| 排障记录 | 机制归并版可摘取“collapse-gated wide+single-MF 是当前 fixed-tail 最稳候选”。 |

## 8. 限制与禁止解释

- 不要把 `gated-wide-single-bank` 从 fixed-tail replay 直接升级为 default flow。
- 不要把 `wide-centered`、`single-MF-centered` 或二者组合作为 blanket 常驻候选。
- 不要把 `tailClass`、truth line、truth DoA、truth `fdRef/fdRate` 放入 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 objective gain 或 coherence gain 单独当作 adoption 充分条件；seed 284 已证明会误伤。
- 不要把本 replay 写成 99% hit-rate 证据；它是机制定位。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothTailCaseDiagnose.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 replay 承接 `replayMfInToothFdRangeOracle_20260426-104249` 的 tail seeds。
- 当前结果把候选从“泛化 joint refine / final-centered polish”收缩为“collapse-gated wide+single-MF basin-entry”。
