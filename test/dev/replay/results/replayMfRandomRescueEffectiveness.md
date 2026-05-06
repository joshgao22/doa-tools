# replayMfRandomRescueEffectiveness 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / small-MC stress` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat` |
| 当前一句话结论 | rescue/random subset bank 明显提升 central tooth 命中率并压低 wrong-tooth tail，但已有 easy-case damage 且运行成本约 `1.74x`；因此只能作为 conditional tooth rescue 候选，不能 blanket 常驻默认 flow。 |
| 决策影响 | 支持保留 `curated4 / random1` 作为离线或条件触发候选；系统性边际贡献应交给 `scanMfSubsetBankCoverage`，不在本 replay 中继续膨胀 bank ablation。 |
| 下一步动作 | 若继续推进，应先做 no-truth gate / scan 验证，重点看 selected margin、tooth confidence、candidate score separation 与 easy damage；不要直接改默认候选池。 |
| 禁止误用 | 不能用 truth tooth 或人工 seed 标签设计 trigger；不能只看 tooth hit 提升而忽略 easy damage、成本和 angle RMSE；不能把 24-repeat 统计写成 regression contract。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfRandomRescueEffectiveness.m`
- 结果文档：`test/dev/replay/results/replayMfRandomRescueEffectiveness.md`
- replay 类型：小 Monte Carlo / subset bank rescue diagnostic
- 主要问题：比较 curated-only subset bank 与 curated+rescue/random subset bank 在真实小 MC 中是否能救回 wrong-tooth case，以及是否会误伤 easy case。
- 观察范围：`numRepeat=24` 的 curated-only vs rescue/random 对照、subset label 评估/选中次数、central tooth rescue、easy damage 与代表性 seed。
- 不覆盖范围：不系统搜索所有 subset bank 组合，不验证 final flow gate，不解决 same-tooth DoA tail，不证明 estimator 主核改变。
- truth 使用口径：truth 只用于 toothIdx、central tooth hit、rescue / damage offline label 和结果评价；不能进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `curated-only` | 只使用基础 curated subset bank 的 flow 对照。 | No | subset candidate 覆盖较少。 | baseline；用于观察默认 bank 的 wrong-tooth tail。 |
| `rescue/random bank` | 在 curated bank 外增加 `curated3 / curated4 / random1` 等扩展候选。 | No | 增加 subset candidate 覆盖和运行成本。 | 若命中 central tooth 增加，说明覆盖有价值；若有 damage，说明不能 blanket 常驻。 |
| `curated3 / curated4` | 手工设计的额外 deterministic subset。 | No | 增加可解释的非默认 subset 候选。 | 可用于 one-tooth rescue，但边际贡献需 scan 系统比较。 |
| `random1` | 条件生成的随机 subset 候选。 | No | 增加随机覆盖，成本较高。 | 对严重 wrong-tooth 有救力，但必须 gate，不能所有 seed 常驻。 |
| `toothIdx` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于结果评价。 | central tooth hit 是选齿指标，不能进入 runtime selector。 |
| `central tooth rescue` | curated-only 不在 tooth 0，而 rescue/random 结果回到 tooth 0。 | Evaluation only | 只改结果标签。 | 说明扩展 bank 有救齿价值，但不代表 angle / CRB 表现改善。 |
| `easy-case damage` | curated-only 已经健康，但 rescue/random 选中后变差。 | Evaluation only | 只改结果标签。 | damage 优先级高于单纯 rescue rate；出现 damage 后不能 blanket 打开。 |
| `random evaluated / selected rate` | random 候选被实际评估 / 被选为 winner 的比例。 | No | 反映成本和实际 winner 贡献。 | evaluated 高代表成本，selected 高代表有贡献；两者都需 gate 平衡。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat` | 2026-04-27 | `representative` | `snrDb=10`，`baseSeed=253`，`numRepeat=24`，`rescueRandomTrials=1`，`repeatMode=parfor-auto` | central tooth hit 从 `17/24` 提升到 `22/24`，mean abs toothIdx 从 `1.9583` 降到 `0.083333`，但出现 `1/24` easy damage 且耗时约 `1.74x`。 | 当前 24-repeat 代表性结果；覆盖早期 12-repeat smoke 的方向性结论。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 24`
- `rescueRandomTrials = 1`
- repeat mode：`parfor-auto`
- 对照策略：curated-only vs curated+rescue/random bank
- snapshot 保存变量：`replayData`
- curated-only batch 耗时约 `306.16 s`
- rescue/random batch 耗时约 `533.33 s`
- 总耗时约 `839.49 s`
- rescue/random 相对 curated-only 耗时约 `1.74x`

### 4.2 aggregate summary

| 指标 | curated-only | rescue/random bank | 变化 / 解释 |
|---|---:|---:|---|
| tooth=0 命中率 | `17/24 = 0.70833` | `22/24 = 0.91667` | 增加 `5/24`，救齿效果明确。 |
| mean abs toothIdx | `1.9583` | `0.083333` | wrong-tooth tail 被显著压缩。 |
| angle RMSE (deg) | `0.0019948` | `0.0020495` | 略差；该 replay 主要改善选齿，不解决 same-tooth DoA tail。 |
| random evaluated rate | — | `10/24 = 0.41667` | random 候选是条件评估，但仍有明显成本。 |
| random selected rate | — | `6/24 = 0.25` | `random1` 有实际 winner 贡献。 |
| central tooth rescue rate | — | `5/24 = 0.20833` | 与 tooth hit 提升一致。 |
| easy-case damage rate | — | `1/24 = 0.041667` | 已出现误伤，需要 no-truth gate。 |
| mean angle delta (deg) | — | `-6.7969e-05` | rescue angle 平均略大，不能把 angle 改善作为主结论。 |

### 4.3 subset bank coverage

| subset label | curated eval count | rescue eval count | curated selected count | rescue selected count | 观察 |
|---|---:|---:|---:|---:|---|
| `curated1` | `24` | `24` | `11` | `6` | primary bank 仍有贡献。 |
| `curated2` | `24` | `24` | `13` | `2` | curated-only 主要依赖 `curated1/2`。 |
| `curated3` | `0` | `24` | `0` | `3` | 24-repeat 下已有少量贡献，不应直接删除。 |
| `curated4` | `0` | `24` | `0` | `7` | 当前 rescue 中最强的 deterministic 扩展候选。 |
| `random1` | `0` | `10` | `0` | `6` | random 有明确 hard rescue 价值，但成本与误伤风险需要 gate。 |

### 4.4 代表性 seed 对比

| seed | curated subset | rescue subset | curated tooth | rescue tooth | 关键变化 | 解释 |
|---:|---|---|---:|---:|---|---|
| `254` | `curated2` | `random1` | `-19` | `0` | `fdRefErrHz` 从 `-14166` 变为 `252.5`，`fdRateErrHzPerSec` 从 `3233.5` 变为 `-14.615` | 典型 random1 严重 wrong-tooth rescue。 |
| `263` | `curated2` | `curated4` | `1` | `0` | `fdRefErrHz` 从 `749.97` 变为 `-0.029226` | `curated4` 可救 one-tooth error。 |
| `264` | `curated2` | `random1` | `1` | `0` | `fdRefErrHz` 从 `749.87` 变为 `-2.2474`，但 `fdRateErrHzPerSec` 变为 `3233.5` | tooth 被救回，但 rate 可能变差；不能只看 tooth。 |
| `259` | `curated1` | `curated1` | `1` | `1` | random 被评估但未选中 | rescue bank 未必能救所有 one-tooth case。 |
| `261` | `curated1` | `curated4` | `-1` | `-1` | 仍停在 one-tooth error | `curated4` 不是万能救齿，需要 scan / gate 继续比较。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- rescue/random bank 将 central tooth hit 从 `0.70833` 提高到 `0.91667`，并将 mean abs toothIdx 从 `1.9583` 降到 `0.083333`。
- `curated4` selected `7/24`，`random1` selected `6/24`，说明扩展 bank 不是纯冗余。
- seed `254` 证明 `random1` 对严重 wrong-tooth tail 有实际救力。
- 24-repeat 下 `curated3` 也有 `3/24` selected，说明不能仅凭早期 smoke 删除它。

### 5.2 仍未解决或反向的现象

- angle RMSE 从 `0.0019948 deg` 略升到 `0.0020495 deg`，说明该 replay 不支持“same-tooth DoA tail 已修复”。
- 已出现 `1/24` easy-case damage；这比单纯 rescue rate 更重要，说明 blanket rescue/random 不安全。
- rescue/random batch 耗时约为 curated-only 的 `1.74x`，默认常驻会增加长 MC 成本。
- seed `264` 显示 tooth 被救回时，`fdRate` 仍可能变差；后续 damage 分类不能只看 tooth。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `254` | severe wrong-tooth rescue | `curated2` 的 `toothIdx=-19` 被 `random1` 拉回 tooth 0 | 支持 random rescue 有覆盖价值。 |
| `263` | one-tooth rescue | `curated4` 将 tooth `1` 拉回 `0` | 支持 deterministic 扩展 subset 有边际贡献。 |
| `264` | tooth rescued but rate suspicious | tooth 回到 0，但 `fdRate` 变差 | 说明后续 gate / damage 不能只看 tooth hit。 |
| `259` | unrecovered one-tooth | random 被评估但未选中，tooth 仍为 `1` | 说明 rescue bank 不万能。 |
| `261` | unrecovered one-tooth | `curated4` 被选中但 tooth 仍为 `-1` | 支持后续需要 scan / gate，而不是继续堆 replay variants。 |

## 6. 机制解释

### 6.1 当前解释

该 replay 的主要信息是：默认 curated subset bank 对大多数 seed 已经有效，但仍存在少数 wrong-tooth tail；扩展 `curated3 / curated4 / random1` 后，候选覆盖更广，能把部分 severe wrong-tooth 或 one-tooth case 拉回 central tooth。由于改善主要体现在 tooth hit 和 mean abs toothIdx，而 angle RMSE 没有改善，当前机制应解释为 **tooth selection coverage 改善**，不是 same-tooth DoA basin-entry 已解决。

同时，easy-case damage 与运行成本说明扩展 bank 不能无条件常驻。后续如果要利用 `curated4 / random1`，应设计 no-truth 条件触发，例如基于 selected subset margin、candidate score separation、tooth confidence 或 consistency proxy，而不是使用 truth tooth 或人工 label。

### 6.2 这个结果支持什么

- 支持 `curated4 / random1` 对 wrong-tooth rescue 有实际价值。
- 支持把 random rescue 保留为 conditional tooth rescue 候选。
- 支持将系统性 bank 组合比较转到 `scanMfSubsetBankCoverage`。
- 支持继续保留 no-truth-leak regression 护栏，避免把 truth label 误放入 selector。

### 6.3 这个结果不证明什么

- 不证明 `curated4 / random1` 可以 blanket 加入默认 flow。
- 不证明 same-tooth DoA / local-state tail 已解决。
- 不证明 angle RMSE 会因 rescue bank 改善；本轮 angle RMSE 反而略差。
- 不证明当前 24-repeat 统计可写成 regression contract。
- 不证明可以用 truth tooth、truth DoA 或人工 seed 类型设计 runtime trigger。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；subset bank 属于 flow / replay 侧，不下沉 estimator 主核。 |
| flow 默认路径 | 不 blanket 增加 random rescue；若推进，需先通过 no-truth gate / scan 验证。 |
| regression | 不写 24-repeat 统计契约；继续保留 `regressionMfSubsetSelectNoTruthLeak` 与 `regressionMfFastSubsetEscalation`。 |
| replay / scan 下一步 | 转向 `scanMfSubsetBankCoverage` 系统比较 `curated123`、`curated124`、`curated1234`、`curated12_random1`、`curated124_random1`、`fullRescue` 等 bank variants。 |
| 论文图 / 论文口径 | diagnostic-only / stress；不作为 paper-facing MLE-vs-CRB 或 CP/IP 主图证据。 |
| 排障记录 | 若只整理格式，不必新增主记录；若后续 scan 证实可用，再更新机制归并版的 subset bank 条目。 |

## 8. 限制与禁止解释

- 不要把 `random1` 当作默认常驻候选；当前已有成本和 easy damage 信号。
- 不要把当前 replay 扩成 exhaustive bank scan；这会混淆 replay 与 scan 职责。
- 不要仅凭当前结果删除 `curated3`；24-repeat 已显示它有少量 selected 贡献。
- 不要用 truth tooth、truth DoA、truth `fdRef/fdRate` 或人工 seed 标签设计 runtime trigger。
- 不要把 angle RMSE 的轻微变化解释成 estimator 主核改善或退化；本 replay 的主要观测对象是 tooth selection。
- 不要把 central tooth hit rate 单独当作论文成功标准；paper-facing 口径仍应使用 RMSE / P95 / P99 / resolved rate / outlier rate 与 CRB-normalized 指标。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfRandomRescueEffectiveness.m`
```

只运行 `Summary output and plotting` 小节，即可重出 aggregate table、bank coverage、tooth histogram 和代表性 seed 表。

## 10. 历史备注

- 早期 12-repeat smoke 的方向与本轮一致：rescue/random bank 能提升 central tooth hit。
- 本轮 24-repeat 更重要的新增信息是：`curated3` 出现少量 selected 贡献，同时开始出现 `1/24` easy-case damage；因此当前结论从“rescue 有用”收敛为“有用但必须 gate，不能 blanket 常驻”。
