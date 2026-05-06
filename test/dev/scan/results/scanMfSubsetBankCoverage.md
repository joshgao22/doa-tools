# scanMfSubsetBankCoverage 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `negative / diagnostic-only` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfSubsetBankCoverage_20260428-164620.mat` |
| 当前一句话结论 | `scheduleComboQuick` 说明 structured nonuniform schedule 的 lag 指标更好，但没有转化为 selected tooth hit；hybrid bank 只增加 candidate coverage 和成本，当前问题更像 candidate-to-final adoption，而不是 schedule 数量不足。 |
| 论文图定位 | `not for paper`；最多作为 internal negative evidence，用于解释为什么不把 subset schedule / tooth acquisition 作为正文主线。 |
| 决策影响 | 固定为 structured schedule quick screen 负结果；不运行 `scheduleComboConfirm`；不进入 regression；不修改默认 subset bank / flow。 |
| 下一步动作 | 暂停继续堆 structured schedule；若后续仍研究 flow，优先诊断 candidate-to-final carryover、periodic final validation 和 no-truth winner adoption。 |
| 禁止误用 | 不能把 candidate truth coverage 增加解释成默认 flow 已改善；truth 字段只用于 offline evaluation，不能进入 runtime selector、gate、candidate adoption 或 final winner。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfSubsetBankCoverage.m`
- 结果文档：`test/dev/scan/results/scanMfSubsetBankCoverage.md`
- scan 类型：`schedule / subset bank scan`，当前结果为 `negative / diagnostic-only`
- 主要问题：比较 legacy curated subset bank、structured nonuniform schedule 和 hybrid bank 对 simple dynamic flow 选齿覆盖、误伤风险、候选成本和 no-truth consensus 的影响。
- 扫描对象：`strategyPreset=scheduleComboQuick`；`curated12 / curated124 / structuredCombo / curated12_structured`；`seedList=253:260`；`SNR=10 dB`。
- 不覆盖范围：不验证 estimator 默认路径；不验证 full paper-facing MLE-vs-CRB；不证明 structured schedule 的一般最优性；不替代 flow-like replay。
- truth 使用口径：truth 只用于 offline evaluation，包括 `truthTooth*`、`nearTooth*`、`toothIdx`、`toothResidualHz`、`angleErrDeg`、`fdRefErrHz`、baseline damage / transition label；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：No。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `curated12` | legacy empirical baseline，包含 `curated1 + curated2` 两条 subset schedule。 | Eval only | 作为当前 subset bank quick screen 的 baseline。 | 不能由本 scan 证明它是全局最优。 |
| `curated124` | 在 `curated12` 基础上增加 `curated4`。 | Eval only | 检查 legacy rescue schedule 是否带来边际 selected gain。 | 候选数增加不等于 flow improvement。 |
| `structuredCombo` | `staggeredA + sparseRulerA + coprime34A` 的纯 structured nonuniform schedule 组合。 | Eval only | 检查较好 lag feature 是否可转化为 selected tooth hit。 | 不能只凭 lag feature 更好就进入默认 bank。 |
| `curated12_structured` | `curated12` 与 structured schedule 的 hybrid bank。 | Eval only | 检查新增 structured 候选是否能救 baseline miss。 | candidate coverage 增加不能解释为 selected winner 已改善。 |
| `truth index hit` | selected candidate 的 tooth index 是否等于 truth center tooth。 | Eval only | 粗看是否选中 center tooth。 | 不能作为 runtime selector。 |
| `near index hit` | `abs(toothIdx) <= 1` 的 relaxed index hit。 | Eval only | 允许邻齿附近候选，用于观察 near-tooth coverage。 | 不能代替 strict paper metric。 |
| `truth strict hit` | center tooth index hit 且 `abs(toothResidualHz) <= 50 Hz`。 | Eval only | 当前 quick screen 的主要 selected tooth 成功指标。 | 不能迁移为默认 flow 的 truth-aware gate。 |
| `easy damage` | baseline angle `<=0.005 deg`，但新策略 angle 增量 `>0.001 deg` 或 lost truth tooth。 | Eval only | 衡量新 bank 对 easy case 的误伤。 | 不是 proof of estimator failure；只是 strategy-level damage label。 |
| `candidateAnyTruthIndexRate` | 某策略候选池里是否至少存在 truth-index candidate。 | Eval only | 衡量候选覆盖上限。 | 不能单独代表 selected winner 可用。 |
| `selectedMissDespite*CandidateRate` | 候选池里已有 truth / near candidate，但最终 selected 仍 miss 的比例。 | Eval only | 指向 candidate-to-final adoption / ranking 问题。 | 不能用来做 runtime truth-aware adoption。 |
| `consensus diagnostic` | no-truth anchor-group residual consensus 诊断。 | No direct truth in decision; truth only used after-the-fact for scoring | 观察无 truth 的 group consensus 是否可改善 selected hit。 | 当前结果不支持下沉为默认 winner rule。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：每个策略 8 个 repeat，共 32 个 strategy-repeat 结果，用于 schedule quick screen。
- `resolved-sample`：N/A；本 scan 不做 CRB / MLE resolved 对比。
- `outlier rate`：N/A；用 `truth strict miss`、`near miss`、`fdRef RMSE / tail` 和 `easy damage` 作为 schedule-level failure signal。
- `truth-tooth / oracle range`：truth 只用于 offline evaluation / label。
- `stress-test`：不是 estimator stress；是 subset schedule quick screen / negative diagnostic。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfSubsetBankCoverage_20260428-164620.mat` | 2026-04-28 | `representative` | `strategyPreset=scheduleComboQuick`；`baseSeed=253`；`seedList=253:260`；`SNR=10 dB`；`numRepeat=8`；baseline `curated12`；near-tooth rule `abs(toothIdx)<=1 && abs(toothResidualHz)<=50 Hz`；easy damage rule baseline angle `<=0.005 deg` and angle delta `>0.001 deg` or lost truth tooth。 | structured schedule quick screen 是负结果；hybrid 只提高 candidate coverage，没有改善 selected hit 或 fd tail。 | none |

维护说明：当前只保留这一组 `scheduleComboQuick` 结果作为代表性 snapshot，不混放早期 `cheapScreen / scheduleDesign / diagnosticCore` 数据。

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = [253, 254, 255, 256, 257, 258, 259, 260]`
- `numRepeat = 8` per strategy
- `snrDb = 10`
- `frameIntvlSec = 1/750`
- `toothStepHz = 750`
- `strategyPreset = scheduleComboQuick`
- `baselineStrategyName = curated12`
- `toothResidualTolHz = 50`
- `nearToothIdxTol = 1`
- `easyAngleTolDeg = 0.005`
- `damageAngleTolDeg = 0.001`
- `toothConsensusResidualPenaltyHz = 75`
- 关键 scan 轴：subset bank strategy / schedule family / bank size。
- 关键 offline 判据：truth strict hit、near strict hit、candidate coverage、selected miss despite candidate、easy damage。
- checkpoint：enabled / resume enabled / cleanup on success enabled；每个 strategy 单独 checkpoint。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录 `elapsedSec`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`scanName`、`runKey`、`config`、`strategyList`、`contextSummary`、`checkpointSummaryTable`、`checkpointCleanupTable`、`scanTable`、`candidateTable`、`candidateSeedCoverageTable`、`consensusTable`、`consensusAggregateTable`、`scheduleFeatureTable`、`transitionTable`、`aggregateTable`、`toothHistogramTable`、`representative`、`plotData`、`repeatCellByStrategy`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：代表性结果中未见 hard fail；本结果主要用于 strategy-level quick screen。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

`truth strict hit` 与 `near strict hit` 是本 scan 的主要 selected 结果指标；`candidate` 相关字段只说明候选池上限，不能代替 selected winner 成功。

| strategy | bank size | truth index hit | near index hit | truth strict hit | near strict hit | easy damage | angle RMSE (deg) | angle P95 (deg) | fdRef RMSE (Hz) | candidates / repeat | selected miss despite truth cand. | selected miss despite near cand. | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `curated12` | 2 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 2 | 0 | 0.125 | legacy baseline。 |
| `curated124` | 3 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 3 | 0 | 0.125 | 增加 `curated4` 没有 selected gain。 |
| `structuredCombo` | 3 | 0.375 | 0.50 | 0.375 | 0.50 | 0.375 | 0.0026108 | 0.004673 | 10460 | 3 | 0.375 | 0.25 | 纯 structured bank 明显差，且误伤 easy。 |
| `curated12_structured` | 5 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 5 | 0.125 | 0.125 | hybrid 只增加成本，没有 selected gain。 |

### 5.2 按扫描轴汇总

| axis value | case | metric 1 | metric 2 | metric 3 | 解释 |
|---:|---|---:|---:|---:|---|
| `bankSize=2 -> 3` | `curated12 -> curated124` | truth strict `0.75 -> 0.75` | fdRef RMSE `8116.8 -> 8116.8 Hz` | candidates `2 -> 3` | legacy `curated4` 没有在本组 seed 中形成边际收益。 |
| structured-only | `structuredCombo` | truth strict `0.375` | near strict `0.50` | easy damage `0.375` | lag feature 更好的 pure structured bank 反而 selected 更差。 |
| hybrid | `curated12_structured` | truth strict `0.75` | fdRef RMSE `8116.8 Hz` | candidates `5` | 候选数和覆盖增加，但 selected 结果与 baseline 相同。 |
| candidate coverage | `curated12_structured` | candidate truth coverage `0.875` | selected strict `0.75` | selected miss despite truth cand. `0.125` | 主要问题转向 candidate-to-final adoption / validation。 |

### 5.3 Schedule feature summary

| schedule | family | num frame | aperture | gcd lag | unique lag | lag redundancy | offset list |
|---|---|---:|---:|---:|---:|---:|---|
| `curated1` | legacyCurated | 10 | 17 | 1 | 17 | 2.6471 | `-8,-7,-5,-4,-2,-1,0,4,7,9` |
| `curated2` | legacyCurated | 10 | 17 | 1 | 17 | 2.6471 | `-7,-4,-1,0,3,5,7,8,9,10` |
| `curated4` | legacyCurated | 10 | 16 | 1 | 15 | 3.0000 | `-6,-3,-2,-1,0,2,3,6,8,10` |
| `staggeredA` | staggered | 10 | 19 | 1 | 19 | 2.3684 | `-9,-8,-6,-3,0,1,3,6,8,10` |
| `sparseRulerA` | sparseRuler | 10 | 19 | 1 | 19 | 2.3684 | `-9,-8,-5,-1,0,1,4,7,9,10` |
| `coprime34A` | coprime | 10 | 18 | 1 | 18 | 2.5000 | `-9,-8,-6,-4,-3,0,3,4,8,9` |

`staggeredA / sparseRulerA / coprime34A` 的 `aperture / unique lag / lag redundancy` 指标比 legacy curated 更整齐，但本次 selected 结果没有跟着改善；这就是本 scan 最重要的负结果。

### 5.4 No-truth consensus diagnostic

| strategy | current truth strict | consensus truth strict | current near strict | consensus near strict | consensus truth gain | consensus damage | mean candidate count | mean anchor group count | 解释 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `curated12` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 2 | 1.375 | consensus 不改变 baseline。 |
| `curated124` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 3 | 1.625 | 无 selected gain。 |
| `structuredCombo` | 0.375 | 0.50 | 0.50 | 0.50 | 0.125 | 0 | 3 | 2.000 | consensus 对 structured-only 有小修复，但仍低于 baseline。 |
| `curated12_structured` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 5 | 2.625 | hybrid 的 consensus 没有收益。 |

### 5.5 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `Selected tooth hit by strategy` | strategy | truth / near strict hit | four strategies | No | 用于 internal negative screen；truth only offline。 |
| `Candidate coverage vs selected hit` | strategy | candidate coverage and selected hit | four strategies | No | 重点展示 candidate coverage 没转成 selected gain。 |
| `Easy damage by strategy` | strategy | easy damage rate | four strategies | No | 说明 pure structured bank 会误伤 easy case。 |
| `Schedule feature comparison` | schedule | aperture / unique lag / redundancy | six schedules | Diagnostic only | lag feature 只解释 schedule 结构，不解释 selected performance。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `curated124` 在 8-repeat quick screen 中没有超过 `curated12`：truth / near strict hit、angle tail、`fdRef RMSE` 都相同，只是候选数从 2 增到 3。
- `structuredCombo` 明显不适合作为单独 bank：truth strict hit 只有 `0.375`，near strict hit 只有 `0.50`，easy damage 达到 `0.375`，`fdRef RMSE=10460 Hz` 高于 baseline。
- `curated12_structured` 没有带来 selected 命中率提升，angle 与 fd tail 也没有改善；它只是把候选数从 2 增加到 5。
- `curated12_structured` 的 candidate coverage 提高到约 `0.875`，但 selected strict hit 仍停在 `0.75`，说明新增 candidate 没有被 final winner 正确采用。
- no-truth consensus diagnostic 对 `structuredCombo` 只有小修复，对 hybrid 没有收益，因此当前 consensus 也不足以下沉到默认 flow。

### 6.2 反向、污染或未解决现象

- 该 scan 是 quick screen，`numRepeat=8` 很小；负结果足以阻止 `scheduleComboConfirm`，但不等价于所有 structured schedule 都无价值。
- `candidateAny*` 提升说明候选池中确实出现过更接近 truth 的 candidate；未转成 selected gain 暗示 final validation / winner adoption 仍可能有研究价值。
- `fdRef RMSE` 对少数 wrong-tooth seed 很敏感；本 scan 不应被用于 paper-facing estimator 精度结论。
- consensus diagnostic 只是 no-truth 后验诊断，不是已经验证的 runtime selector。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `seed=254` under `curated12` | baseline wrong-tooth seed | baseline 主要 miss seed 之一。 | 用于说明 legacy baseline 仍有可救空间。 |
| `seed=259` under `curated12` | baseline wrong-tooth seed | baseline 主要 miss seed 之一。 | 用于后续 candidate-to-final adoption 诊断。 |
| `curated124` on baseline miss seeds | noStrictChange | 对 `254 / 259` 没有 strict rescue。 | 说明 `curated4` 本轮没有边际收益。 |
| `curated12_structured` | candidate-without-selected-gain | 候选覆盖增加，但 selected hit 不变。 | 指向 adoption / final validation，而不是继续堆 schedule。 |
| `structuredCombo` | easy damage | easy damage `0.375`。 | 说明 pure structured bank 的 winner selection 不稳定。 |
| representative damage sample | selected wrong tooth | `taskSeed=260`；`angleErrDeg=0.002612`；`fdRefErrHz=-21000.0`；`toothIdx=-28`；`toothResidualHz≈-0.0055`；`easyDamageFromBaseline=1`。 | 说明 bank 可在 easy baseline 上选到远齿，不能 blanket 下沉。 |

## 7. 机制解释

### 7.1 当前解释

这个 scan 的核心负结果是：structured nonuniform schedule 在 lag feature 上更“好看”，但 selected winner 结果没有随之改善。原因很可能不是 candidate 完全不存在，而是候选到最终 periodic refine / final winner 的承接不稳定。换句话说，schedule 设计改善了搜索空间覆盖，但当前 flow 的 candidate ranking、final validation 或 adoption 没有稳定把好 candidate 变成最终解。

因此，本结果不支持继续沿“新增更多 structured schedule”单线推进。若后续还研究 tooth acquisition / low-complexity flow，应先看 no-truth final validation、top-K carryover、candidate family retention 或 selected-to-final consistency，而不是把 structured schedule blanket 加进默认 bank。

从论文主线看，这个 scan 也支持当前阶段切换：global tooth acquisition / subset schedule 是工程诊断层问题，不应成为正文主精度目标。正文应继续转向 Doppler-aided / in-tooth local MLE、CRB consistency、CP/IP trade-off 和 known/unknown-rate 信息损失。

### 7.2 这个 scan 支持什么

- 支持 structured schedule quick screen 当前为负结果。
- 支持 `curated12_structured` 的 candidate coverage 增加没有转化为 selected hit。
- 支持后续若继续修 flow，应优先看 candidate-to-final adoption / validation，而不是继续堆 schedule。
- 支持 truth 只作为 offline evaluation / results label，不进入 runtime selector、gate、adoption 或 final winner。
- 支持把 subset schedule 相关工作降级为 engineering diagnostic，不作为论文主仿真线。

### 7.3 这个 scan 不证明什么

- 不证明所有 structured schedule 都必然无效。
- 不证明 `curated12` 是最终默认最优 bank。
- 不证明 estimator 默认路径有数值错误。
- 不证明 full-flow 已经可以用于 paper-facing CP/IP 或 CRB 图。
- 不证明 candidate-to-final adoption 修复无价值。
- 不证明可以写 regression 契约。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 只评价 subset schedule / flow candidate，不触碰 estimator 主核。 |
| flow 默认路径 | 不改；structured schedule 不进入默认 bank，`curated12_structured` 不进入默认 hybrid。 |
| replay 下一步 | 暂不新增 schedule replay；若继续研究，应针对 candidate-to-final adoption 建立 flow-like replay，而不是扩 `scheduleComboConfirm`。 |
| regression | 不写；当前是 negative quick screen，不是稳定 pass/fail 契约。 |
| 论文图 | 不作为论文图；最多作为内部负结果，解释为何论文主线不围绕 tooth acquisition / schedule 展开。 |
| 排障记录 | 保留一句“structured schedule quick screen 为负结果，后续优先看 candidate-to-final adoption”；不复制长表。 |

## 9. 限制与禁止解释

- 不要把 offline truth 评价字段用于 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 candidate coverage 增加当成 selected flow 已改善。
- 不要把 structured lag feature 更好当成进入默认 bank 的充分依据。
- 不要因为 `structuredCombo` 负结果就否定所有可能的 no-truth final validation / top-K carryover 方向。
- 不要把本 scan 的 truth strict hit / fdRef RMSE 用作 paper-facing estimator 精度结论。
- 不要把该 quick screen 迁移成 regression。
- 不要继续扩大 `staggeredB / sparseRulerB / coprimeB`，除非先有 candidate-to-final adoption 的新证据。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfSubsetBankCoverage_20260428-164620.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfSubsetBankCoverage.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：若当前脚本已经按 template 化，checkpoint / Telegram 外壳可能和该历史 snapshot 不完全一致；机制结论以本 snapshot 的 `scheduleComboQuick` 配置为准。

## 11. 历史备注

* 当前只绑定 `scanMfSubsetBankCoverage_20260428-164620.mat` 作为代表性 structured schedule quick screen 结果。
* 本文件不再混放早期 `cheapScreen / scheduleDesign / diagnosticCore` 数据，避免把旧实验口径和当前 `scheduleComboQuick` 结论混在一起。
* 当前结果已经足以暂停 `scheduleComboConfirm`；若未来 flow adoption 逻辑发生实质变化，再重新运行并把本 snapshot 标为 `superseded`。
*