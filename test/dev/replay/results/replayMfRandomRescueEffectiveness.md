# replayMfRandomRescueEffectiveness 结果记录

## 对应 replay

- 脚本：`test/dev/replay/replayMfRandomRescueEffectiveness.m`
- 目标：比较 curated-only subset bank 与 curated+rescue/random subset bank 在小 Monte Carlo 中的选齿效果、wrong-tooth rescue、easy-case damage 与候选 bank 覆盖情况。
- 定位：该 replay 只回答“扩展 rescue bank 是否有救齿价值、是否有误伤风险”；系统性比较 `curated3 / curated4 / random1 / fullRescue` 的边际贡献应放到 `scanMfSubsetBankCoverage`，不在本 replay 中继续膨胀 bank ablation。

## Snapshot index

| snapshot / log | 状态 | 配置摘要 | 备注 |
|---|---|---|---|
| `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat` | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=24`，`rescueRandomTrials=1`，`repeatMode=parfor-auto`，`saveSnapshot=1` | 当前保留的 24-repeat curated rescue 代表性结果；只保存轻量 `replayData` |

## 2026-04-27 11:04 结果分析

### 配置与总体状态

- 本轮使用 `24` 个 repeat，base seed 为 `253`，SNR 为 `10 dB`。
- curated-only batch 耗时约 `306.16 s`；curated-rescue/random batch 耗时约 `533.33 s`；总运行时间约 `839.49 s`。
- rescue 侧相对 curated-only 的运行时间约为 `1.74x`，说明扩展 bank 有明显计算成本，不能在默认 flow 中无条件常驻。
- snapshot 已保存为 `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat`，保存变量为 `replayData`。
- `toothHistogramTable` 共 `40` 行，完整 tooth 分布保存在 `replayData.toothHistogramTable` 中。
- 当前一句话结论：rescue/random bank 明显提升 central tooth 命中率并显著压低 mean abs toothIdx，但已经出现 `1/24` 的 easy-case damage，因此后续应验证 no-truth 条件触发 gate，而不是把 rescue/random blanket 加入默认 flow。

### Aggregate summary

| 指标 | curated-only | rescue/random bank | 变化 / 解释 |
|---|---:|---:|---|
| tooth=0 命中率 | `17/24 = 0.70833` | `22/24 = 0.91667` | 增加 `5/24`，救齿效果明确 |
| mean abs toothIdx | `1.9583` | `0.083333` | wrong-tooth tail 被显著压缩 |
| angle RMSE (deg) | `0.0019948` | `0.0020495` | 略差，说明该 replay 主要改善选齿，不解决 same-tooth DoA tail |
| random evaluated rate | — | `10/24 = 0.41667` | random 候选是条件评估，不是所有 seed 都进入 |
| random selected rate | — | `6/24 = 0.25` | `random1` 有实际 winner 贡献 |
| central tooth rescue rate | — | `5/24 = 0.20833` | 与 tooth hit 提升一致 |
| easy-case damage rate | — | `1/24 = 0.041667` | 已出现误伤，需要 gate 约束 |
| mean angle delta (deg) | — | `-6.7969e-05` | rescue angle 平均略大，不能把 angle 改善作为本 replay 主结论 |

### Subset bank coverage

| subset label | curated eval count | rescue eval count | curated selected count | rescue selected count | 观察 |
|---|---:|---:|---:|---:|---|
| `curated1` | 24 | 24 | 11 | 6 | primary bank 仍有贡献 |
| `curated2` | 24 | 24 | 13 | 2 | curated-only 主要依赖 `curated1/2` |
| `curated3` | 0 | 24 | 0 | 3 | 24-repeat 下已有少量贡献，不应直接删除 |
| `curated4` | 0 | 24 | 0 | 7 | 当前 rescue 中最强的 curated 扩展候选 |
| `random1` | 0 | 10 | 0 | 6 | random 有明确 hard rescue 价值，但成本与误伤风险需要 gate |

### 代表性 seed 观察

| seed | curated subset | rescue subset | curated tooth | rescue tooth | 关键变化 | 解释 |
|---:|---|---|---:|---:|---|---|
| 254 | `curated2` | `random1` | -19 | 0 | `fdRefErrHz` 从 `-14166` 变为 `252.5`，`fdRateErrHzPerSec` 从 `3233.5` 变为 `-14.615` | 典型 random1 严重 wrong-tooth rescue |
| 263 | `curated2` | `curated4` | 1 | 0 | `fdRefErrHz` 从 `749.97` 变为 `-0.029226` | `curated4` 可救 one-tooth error |
| 264 | `curated2` | `random1` | 1 | 0 | `fdRefErrHz` 从 `749.87` 变为 `-2.2474`，但 `fdRateErrHzPerSec` 变为 `3233.5` | tooth 被救回，但 rate 可能变差；后续 damage 分类不能只看 tooth |
| 259 | `curated1` | `curated1` | 1 | 1 | random 被评估但未选中 | rescue bank 未必能救所有 one-tooth case |
| 261 | `curated1` | `curated4` | -1 | -1 | 仍停在 one-tooth error | `curated4` 不是万能救齿；需要 scan 继续比较 bank / gate |

### 与早期 12-repeat smoke 的关系

- 方向保持一致：rescue/random bank 继续显著提升 central tooth hit，并降低 wrong-tooth tail。
- 24-repeat 下 `curated3` 从“未被选中”变成 `3/24` selected，说明不能仅凭 12-repeat smoke 删除 `curated3`。
- 24-repeat 下 `curated4` selected `7/24`，`random1` selected `6/24`，二者仍是最值得系统比较的 rescue 候选。
- 24-repeat 下开始出现 `1/24` easy-case damage，因此后续重点应从“rescue 有没有用”转为“如何 no-truth 条件触发”。

### 可观察现象

- curated-only 已能在多数 seed 上命中 central tooth，但仍存在少数严重 wrong-tooth tail，例如 seed `254` 的 `toothIdx=-19`。
- 扩展 rescue bank 后，central tooth hit 从 `17/24` 提升到 `22/24`，mean abs toothIdx 从 `1.9583` 降到 `0.083333`，说明 rescue bank 的主要价值是选齿覆盖，而不是角度精修。
- `curated4` 与 `random1` 都有实际 winner 贡献；其中 `random1` 对严重 wrong-tooth rescue 更明显，`curated4` 对 one-tooth rescue 更可解释。
- angle RMSE 未改善，且略有恶化；因此该 replay 不能支持“rescue bank 已解决 same-tooth DoA/local-state hard case”的结论。
- rescue 侧耗时明显增加，且已经有 easy-case damage；因此 blanket rescue/random 不安全。
- `randomEvaluatedRate=0.41667` 与 `randomSelectedRate=0.25` 表明 random1 不只是冗余候选，但应被 gate 限制在 hard case 中。

## 当前结论

- `replayMfRandomRescueEffectiveness` 已完成当前 replay 职责：证明 curated+rescue/random bank 在真实小 MC 中确实能救回部分 wrong-tooth case。
- 当前结果支持继续保留 `random1` / rescue subset 作为 conditional tooth rescue 候选。
- 当前结果不支持把 `curated4` 或 `random1` blanket 加入所有默认候选池，因为成本约 `1.74x`，且出现 `1/24` easy-case damage。
- 当前结果也不支持把 same-tooth DoA tail 归因于 subset bank 覆盖不足；angle RMSE 基本不变，same-tooth hard case 仍应由 `replayMfInToothTailCaseDiagnose` / gated basin-entry replay 继续定位。
- 系统性比较 `curated3 / curated4 / random1 / fullRescue` 的边际贡献应转到 `scanMfSubsetBankCoverage`，本 replay 不继续扩成 bank ablation scan。

## 对主流程的影响

- `buildSimpleDynamicFlowOpt` 的默认 flow 不应因为本 replay 直接 blanket 增加 random rescue。
- 下一步更合理的是在 `scanMfSubsetBankCoverage` 中系统比较：`curated123`、`curated124`、`curated1234`、`curated12_random1`、`curated124_random1`、`fullRescue` 等 bank variants。
- 若 scan 仍显示 `curated4 / random1` 有稳定边际贡献，再回到 flow 层设计 no-truth gate，例如基于 selected subset margin、tooth confidence、candidate score separation 或 same-tooth consistency，而不是基于 truth tooth。
- 回归护栏仍应保持 `regressionMfSubsetSelectNoTruthLeak` 与 `regressionMfFastSubsetEscalation`；本 replay 的 24-repeat 统计不应直接固化成 regression contract。
- 后续若修改 flow adoption 或 rescue gate，建议至少运行：`regressionMfSubsetSelectNoTruthLeak`、`regressionMfFastSubsetEscalation`、`runRegressionDynamicFlow`，并重新跑小 repeat replay 确认 easy-case damage 没有扩大。

## 暂不建议

- 不在 `replayMfRandomRescueEffectiveness` 中继续加入 exhaustive bank variants；这会把 replay 变成 scan。
- 不把 `random1` 当作默认常驻候选；当前已有成本与误伤信号。
- 不仅凭当前 replay 删除 `curated3`；24-repeat 已显示 `curated3` 有少量 selected 贡献。
- 不用 truth tooth 或人工 seed 标签设计 trigger；这些只能用于结果标注和离线评价。
- 不把 angle RMSE 的轻微变化解释成 estimator 主核改善或退化；该 replay 的主要观测对象是 tooth selection。

## 恢复方式

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后运行 `replayMfRandomRescueEffectiveness.m` 的 `Summary output and plotting` section。
