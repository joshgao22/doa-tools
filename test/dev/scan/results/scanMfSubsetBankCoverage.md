# scanMfSubsetBankCoverage 结果记录

## 对应 scan

- `test/dev/scan/scanMfSubsetBankCoverage.m`

## 扫描目标

`scanMfSubsetBankCoverage` 用于离线比较 subset bank / structured nonuniform schedule 对真实 simple dynamic flow 的选齿覆盖、误伤风险和候选成本。该 scan 可以用 truth 做结果评价，但不能把 truth 反馈到 selector、gate、periodic refine 或默认 flow。

## Snapshot index

- `test/data/cache/scan/scanMfSubsetBankCoverage_20260428-164620.mat`

## 本次运行配置

| item | value |
|---|---|
| strategy preset | `scheduleComboQuick` |
| repeat seeds | `253:260` |
| repeats per strategy | 8 |
| SNR | 10 dB |
| baseline strategy | `curated12` |
| near-tooth rule | `abs(toothIdx) <= 1` and `abs(toothResidualHz) <= 50 Hz` |
| easy damage rule | baseline angle `<= 0.005 deg`，且 angle 增量 `> 0.001 deg` 或 lost truth tooth |
| snapshot | `scanMfSubsetBankCoverage_20260428-164620.mat` |

本次只保留这一组 `scheduleComboQuick` 结果作为当前结果记录，不再混放早期 `cheapScreen / scheduleDesign / diagnosticCore` 数据。

## Strategy definition

| strategy | bank size | schedule labels | 角色 |
|---|---:|---|---|
| `curated12` | 2 | `curated1 + curated2` | legacy empirical baseline |
| `curated124` | 3 | `curated1 + curated2 + curated4` | 当前历史结果中相对有价值的 legacy rescue 对照 |
| `structuredCombo` | 3 | `staggeredA + sparseRulerA + coprime34A` | 纯 structured nonuniform schedule 组合 |
| `curated12_structured` | 5 | `curated1 + curated2 + staggeredA + sparseRulerA + coprime34A` | legacy baseline + structured schedule 的 hybrid bank |

## 选择策略原则

1. **truth 只用于离线评价，不进入选择器。** 允许在 `scanData.aggregateTable / transitionTable / candidateSeedCoverageTable / consensusTable` 中使用 `truthTooth*`、`angleErrDeg`、`fdRefErrHz` 等字段评价 schedule，但真实 flow 不能使用这些字段判断是否 rescue、是否可信或选哪个 tooth。
2. **先看 selected 结果，再看 candidate coverage。** 一个 schedule 仅提高 `candidateAnyTruthIndexRate` 或 `candidateAnyNearIndexRate` 不足以进入默认 flow；必须同时降低 `selectedMissDespite*CandidateRate`，并带来 selected hit / fd tail 的实际改善。
3. **candidate 数量必须付出可解释收益。** 扩大 bank 会增加候选成本，也会增加误选风险。若 `meanCandidateEvalPerRepeat` 增加但 `truthToothHitRate / nearToothHitRate / fdRefRmseHz` 不改善，则不应进入下一轮确认。
4. **structured schedule 需要证明机制收益，而不是只证明 lag 指标更好。** `aperture` 更大、`numUniqueLag` 更多、`lagRedundancy` 更低，只能说明 schedule 设计更有结构；是否有效仍以 tooth hit、candidate miss、damage 和 fd tail 为准。
5. **不原地替换 legacy curated。** `curated1/2/4` 继续作为历史 baseline 和对照，不直接改写 `getDynamicCuratedSubsetBank.m`。structured schedule 只有在 scan / replay 中稳定超过 legacy 后，才考虑进入 flow 配置。
6. **不把 poor selected 结果交给更重 confirmation。** quick screen 若已经显示 structured bank 明显差或 hybrid 没有 selected gain，就不继续跑 `scheduleComboConfirm`。

## Aggregate result

| strategy | bank size | truth index hit | near index hit | truth strict hit | near strict hit | easy damage | angle RMSE (deg) | angle p95 (deg) | fdRef RMSE (Hz) | candidates / repeat | selected miss despite truth cand. | selected miss despite near cand. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `curated12` | 2 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 2 | 0 | 0.125 |
| `curated124` | 3 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 3 | 0 | 0.125 |
| `structuredCombo` | 3 | 0.375 | 0.50 | 0.375 | 0.50 | 0.375 | 0.0026108 | 0.004673 | 10460 | 3 | 0.375 | 0.25 |
| `curated12_structured` | 5 | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0.0022066 | 0.004673 | 8116.8 | 5 | 0.125 | 0.125 |

### 直接观察

- `curated124` 在本次 8-repeat quick screen 中没有超过 `curated12`。虽然 `curated4` 被选中 3 次，但 selected truth / near hit、angle tail 和 `fdRef RMSE` 都与 baseline 相同。
- `structuredCombo` 明显不适合作为单独 bank：truth strict hit 只有 `0.375`，near strict hit 只有 `0.50`，easy damage 达到 `0.375`，`fdRef RMSE` 也高于 baseline。
- `curated12_structured` 没有带来 selected 命中率提升，angle 与 fd tail 也没有改善；它只是把候选数从 2 增加到 5。
- `curated12_structured` 的 `candidateAnyTruthIndexRate` 和 `candidateAnyTruthToothRate` 提高到 `0.875`，但 selected hit 仍停在 `0.75`，说明新增 structured schedule 主要暴露了 candidate-to-final adoption 问题，而不是直接解决选齿问题。

## Schedule feature summary

| schedule | family | num frame | aperture | gcd lag | unique lag | lag redundancy | offset list |
|---|---|---:|---:|---:|---:|---:|---|
| `curated1` | legacyCurated | 10 | 17 | 1 | 17 | 2.6471 | `-8,-7,-5,-4,-2,-1,0,4,7,9` |
| `curated2` | legacyCurated | 10 | 17 | 1 | 17 | 2.6471 | `-7,-4,-1,0,3,5,7,8,9,10` |
| `curated4` | legacyCurated | 10 | 16 | 1 | 15 | 3.0000 | `-6,-3,-2,-1,0,2,3,6,8,10` |
| `staggeredA` | staggered | 10 | 19 | 1 | 19 | 2.3684 | `-9,-8,-6,-3,0,1,3,6,8,10` |
| `sparseRulerA` | sparseRuler | 10 | 19 | 1 | 19 | 2.3684 | `-9,-8,-5,-1,0,1,4,7,9,10` |
| `coprime34A` | coprime | 10 | 18 | 1 | 18 | 2.5000 | `-9,-8,-6,-4,-3,0,3,4,8,9` |

`staggeredA / sparseRulerA / coprime34A` 的 lag feature 比 legacy curated 更“结构化”，但本次结果说明：结构化 lag 指标没有直接转化成 selected tooth hit。后续不应继续只靠新增更多单条 structured schedule 试错。

## No-truth consensus diagnostic

| strategy | current truth strict | consensus truth strict | current near strict | consensus near strict | consensus truth gain | consensus damage | mean candidate count | mean anchor group count |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `curated12` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 2 | 1.375 |
| `curated124` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 3 | 1.625 |
| `structuredCombo` | 0.375 | 0.50 | 0.50 | 0.50 | 0.125 | 0 | 3 | 2.000 |
| `curated12_structured` | 0.75 | 0.75 | 0.75 | 0.75 | 0 | 0 | 5 | 2.625 |

consensus diagnostic 对 `structuredCombo` 有小幅修复，但修复后仍显著低于 legacy baseline；对 `curated12_structured` 没有收益。当前 no-truth anchor-group consensus 不能作为把 structured bank 下沉到默认 flow 的依据。

## Seed-level observation

- `curated12` 的主要 wrong-tooth seed 是 `254` 和 `259`。
- `curated124` 对这两个 seed 没有 strict rescue；transition 均为 `noStrictChange`。
- `curated12_structured` 虽然增加了候选覆盖，但没有把这些新增候选转成 selected gain。
- `structuredCombo` 在 easy seed 上有明显误伤，说明纯 structured bank 的 winner selection 不稳定。

## 当前结论

1. **本轮 structured schedule 是负结果。** `structuredCombo` selected hit 明显差，且有高 easy damage；不应继续进入 24-repeat confirmation。
2. **hybrid bank 暂时没有价值。** `curated12_structured` 提高了候选覆盖，但没有提升 selected hit，也没有降低 fd tail；候选成本从 2 增到 5，不值得作为默认或确认路线。
3. **`curated124` 本轮没有表现出边际收益。** 它仍可作为 legacy rescue baseline 保留，但本组 seed 下没有超过 `curated12`。
4. **当前问题更像 candidate-to-final adoption，而不是 schedule 数量不足。** 当新增 bank 只提高 candidate coverage、却没有降低 selected miss 或提升 selected hit 时，应停止继续堆 schedule，转向检查 candidate ranking、periodic final validation 和 winner adoption。
5. **本结果不支持继续扩展 `staggeredB / sparseRulerB / coprimeB`。** 若后续还要研究 structured schedule，应先设计 final-level validation，而不是继续单独或组合增加 schedule。

## 下一步建议

当前不建议运行 `scheduleComboConfirm`。更经济的后续路线是：

1. 保留本次结果作为 structured schedule quick screen 的负结果。
2. scan 默认后续可回到更短的 practical rescue 对照：`curated12 / curated124 / curated12_random1 / fullRescue`。
3. 针对 `candidateAny*` 已有但 selected miss 的 seed，新增或增强 candidate-to-final adoption 诊断，重点看 subset candidate 到 `periodic-static-seed / periodic-subset-seed` final 的承接链。
4. 若要进入 flow 修改，应优先验证 no-truth final validation / top-K carryover，而不是把 structured schedule blanket 加进默认 bank。

## 对 replay / regression / 论文口径的影响

- 本次结果只作为 scan 机制结论，不迁移为 regression。
- structured schedule 暂不作为论文低复杂度方法候选；最多作为“尝试过 deterministic nonuniform schedule，但单纯结构化选帧不足以解决 adoption 问题”的内部负结果。
- 若论文保留低复杂度实现，应弱化 schedule 本身，强调主线仍是 continuous-phase + nuisance-rate MLE；subset / schedule 只作为可选工程实现，不应喧宾夺主。

## 历史 / superseded snapshots

本文件只保留当前代表性 snapshot：

- `scanMfSubsetBankCoverage_20260428-164620.mat`
