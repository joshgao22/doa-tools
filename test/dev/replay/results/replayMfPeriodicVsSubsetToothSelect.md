# replayMfPeriodicVsSubsetToothSelect 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / mechanism diagnostic` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat` |
| 当前一句话结论 | subset 的主要作用是选齿，periodic full-data 的主要作用是同齿 refine；final periodic refine 几乎不改变 selected subset 的 tooth，因此不能指望它自动修 wrong-tooth。 |
| 决策影响 | 保持 `selected subset -> final periodic refine` 的 flow 分工；后续优先修 subset rescue / candidate-to-final adoption，而不是把 alias-aware 字段灌回 objective。 |
| 下一步动作 | 保留为 tooth-selection 机制样板；若重跑长 MC，继续保留 warning summary 与 selected/final tooth transition。 |
| 禁止误用 | 不能把 final periodic refine 当作 wrong-tooth rescue；不能用 truth tooth 做 selector；不能因为 warning seed 正常返回就剔除 repeat。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfPeriodicVsSubsetToothSelect.m`
- 结果文档：`test/dev/replay/results/replayMfPeriodicVsSubsetToothSelect.md`
- replay 类型：小 MC mechanism diagnostic。
- 主要问题：periodic wide、selected subset、final periodic refine 三阶段分别承担什么职责。
- 观察范围：真实小 MC flow 中的 tooth index、angle error、fdRef error、warning metadata 和 selected/final tooth transition。
- 不覆盖范围：不验证 in-tooth hard-collapse rescue；不验证 ordinary-wide gate；不改变 estimator 默认路径；不作为 CRB consistency 证据。
- truth 使用口径：truth 只用于计算 `toothIdx`、angle / fd error 和 aggregate 评价；subset selection 与 final refine 不能使用 truth tooth。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `periodic wide` | 在完整等间隔周期帧上做较宽 dynamic search。 | No | 同时影响 tooth 与 DoA basin。 | 暴露 full-data baseline 是否会落入 wrong-tooth；本轮 exact tooth 只有 `42%`。 |
| `selected subset` | 从非周期 subset bank 中按 objective 选出的候选。 | No | 主要改变 tooth selection。 | 核心价值是缩小 tooth error，不等价于最终 angle refine。 |
| `final periodic refine` | 回到完整 periodic data，在 selected subset 附近做 refine。 | No | 主要改善同齿内 angle / Doppler consistency。 | 几乎不改变 selected tooth，不能自动修 wrong-tooth。 |
| `toothIdx` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于评价和结果标注。 | 不进入 runtime selector / gate / final winner。 |
| `abs(toothIdx)<=1` | one-tooth-near 评价口径。 | Evaluation only | 只改评价阈值。 | 可作为辅助观察，但不能替代 exact tooth hit。 |
| `warningTable` | 正常返回 repeat 中记录的非致命 MATLAB warning。 | No | 只改诊断元数据。 | aggregate 不剔除 warning repeat；若 warning 集中在 tail，再单独排查。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat` | 2026-04-25 | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=100`，`checkpointEnable=true`，新增 `warningTable` / `warningSeen` 字段。 | subset 明显改善 tooth selection；final periodic refine 明显改善 angle，但不负责修齿。 | 当前代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 100`
- `toothStepHz = 750`
- checkpoint：`100/100` 个任务完成。
- snapshot 保存变量：`replayData`
- 关键字段：`warningTable`、`repeatTable.warningSeen`、`repeatTable.warningId`、`repeatTable.warningMessage`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| periodic wide tooth=0 | 42/100 | full periodic wide 很容易落 wrong-tooth。 |
| selected subset tooth=0 | 76/100 | subset 显著改善 exact tooth hit。 |
| final periodic refine tooth=0 | 75/100 | final refine 几乎继承 selected subset tooth。 |
| selected subset 修回 periodic wrong-tooth | 38/58 | periodic wrong-tooth 中多数被 subset 拉回 truth tooth。 |
| final / subset tooth 相同 | 99/100 | final refine 主要做同齿细化。 |
| final 改善 angle | 91/100 | final refine 的主要收益是 angle。 |
| warning repeat | 1/100 | 非致命 warning 进入诊断表，不剔除 aggregate。 |

### 4.3 关键对比表

| 阶段 | tooth=0 命中率 | abs(toothIdx)<=1 比例 | abs(toothIdx) median / p95 / max | angle error median / p95 | abs(fdRefErrHz) median / p95 |
|---|---:|---:|---:|---:|---:|
| periodic wide | 42/100 | 51/100 | 1 / 228.55 / 305 | 0.005657 deg / 0.012967 deg | 750.018 Hz / 171.413 kHz |
| selected subset | 76/100 | 87/100 | 0 / 19.15 / 47 | 0.011009 deg / 0.013998 deg | 0.0327 Hz / 14.363 kHz |
| final periodic refine | 75/100 | 87/100 | 0 / 19.15 / 47 | 0.001542 deg / 0.003465 deg | 1.983 Hz / 14.360 kHz |

补充统计：

- selected subset 相比 periodic wide 的 `abs(toothIdx)`：`56/100` 更小，`39/100` 相同，`5/100` 更大。
- final periodic refine 相比 selected subset 的 `abs(toothIdx)`：`99/100` 相同，`1/100` 更大，没有样本跨齿改善。
- selected / final 后仍有 `13/100` 超出 `abs(toothIdx)<=1`，其中 `9/100` 超过 `5` 个 tooth，`5/100` 超过 `20` 个 tooth；没有超过 `50` 个 tooth 的远端错误。
- 本数据中 `usedPolish=0`，只能说明当前流程没有误触发 polish，不能证明 conditional polish 触发后是否有效。

## 5. 可观察现象

### 5.1 支持当前结论的现象

- periodic full-data 可以给出可用 angle，但会落到很远的 wrong tooth；最大达到 `305` 个 tooth。
- selected subset 把 exact tooth hit 从 `42%` 提高到 `76%`，把 `abs(toothIdx)<=1` 从 `51%` 提高到 `87%`，并消除了 `>|50|` tooth 的远端错误。
- final periodic refine 在 `91/100` 样本中降低 angle error，final / subset 的 angle error 中位比值约为 `0.152`。
- final periodic refine 与 selected subset 的 tooth 在 `99/100` 样本中相同，说明它不是修齿步骤。

### 5.2 仍未解决或反向的现象

- selected subset / final periodic refine 后仍有 `13/100` 超出 one-tooth-near，说明 subset bank / rescue 仍不完备。
- selected subset 有 `5/100` 比 periodic wide 的 `abs(toothIdx)` 更大；subset 不是无风险 oracle。
- final periodic refine 在 seed `297` 从 tooth `0` 变到 `-1`，说明 refine 也可能轻微破坏 tooth。
- 当前没有 polish 触发，不能评价 conditional polish 的真实有效性。

### 5.3 warning 记录与代表 seed

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 329 | warning repeat | final periodic refine `toothIdx=-1`，`angleErrDeg=0.001833`，`fdRefErrHz=-750.083`。 | warning repeat 正常返回，aggregate 不剔除；它不是 angle tail，但仍是 one-tooth error。 |
| 297 | tooth transition | final periodic refine 将 tooth 从 `0` 变到 `-1`。 | 说明 final refine 不能被视作只会改善的安全修齿步骤。 |
| periodic wrong-tooth group | wrong-tooth | `58/100` 不在 truth tooth，其中 `38/58` 被 subset 拉回。 | 支持 subset 的核心职责是选齿。 |

## 6. 机制解释

### 6.1 当前解释

完整 periodic data 因等间隔帧结构保留强 `1/T_f` comb，容易在 `fdRef` tooth 上落错；非周期 subset 通过破坏周期结构，更适合在候选层选到正确 tooth 或至少缩小 tooth error。selected subset 的代价是 angle 不一定最细，因为它只使用 subset 数据或非最终配置；因此最终还需要回到完整 periodic data 做同齿内 refine。

final periodic refine 的价值是 angle / Doppler consistency 细化，而不是重新选齿。本轮 `99/100` 样本 final tooth 与 selected subset 相同，说明如果 selected subset 已经错齿，final refine 通常不会自动修回来。

### 6.2 这个结果支持什么

- 支持 flow 分工：subset 负责 tooth selection，periodic full-data 负责 same-tooth refine。
- 支持继续维护 no-truth subset selector / gate 护栏。
- 支持后续关注 subset rescue / random bank / candidate-to-final adoption，而不是只缩 periodic fd 搜索盒。

### 6.3 这个结果不证明什么

- 不证明 selected subset 已经足够稳定；仍有 `13/100` 超出 one-tooth-near。
- 不证明 final periodic refine 可修 wrong-tooth。
- 不证明 conditional polish 无效；本轮没有 polish 触发。
- 不证明可以把 truth tooth、truth DoA 或 truth `fdRef` 用于 runtime selection。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不触碰 estimator 主核。 |
| flow 默认路径 | 保持 subset -> final periodic refine 分工；后续只在 flow 层改 candidate bank / adoption。 |
| regression | 继续依赖 no-truth selection / same-tooth ranking 护栏；本 replay 不新增 regression。 |
| replay / scan 下一步 | 用 `replayMfRandomRescueEffectiveness`、`scanMfSubsetBankCoverage` 检查 rescue bank 与 candidate-to-final adoption。 |
| 论文图 / 论文口径 | diagnostic-only；可解释 full-flow outlier / wrong-tooth 风险，不作为 paper-facing 主图。 |
| 排障记录 | 保留“subset 选齿、periodic refine 同齿细化”的机制结论。 |

## 8. 限制与禁止解释

- 不要把 `abs(toothIdx)<=1` 当成 exact tooth hit 的替代。
- 不要在 runtime selector、gate 或 final winner 中使用 truth tooth。
- 不要把 warning repeat 自动剔除；除非 warning 与 tail 集中相关，否则只作为诊断元数据。
- 不要因为 final periodic refine 能改善 angle，就期望它修 wrong-tooth。
- 不要把本 replay 的机制结果写成 estimator numerical path 的修改依据。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfPeriodicVsSubsetToothSelect.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 snapshot 新增 warning 记录字段，当前处理口径是正常返回的 warning repeat 不剔除 aggregate。
- 该 replay 的长期作用是固定 flow 分工解释，而不是继续承担 same-tooth in-tooth rescue 或 ordinary-wide gate 验证。
