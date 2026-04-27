# replayMfPeriodicVsSubsetToothSelect 结果记录

## 对应 replay

- 脚本：`test/dev/replay/replayMfPeriodicVsSubsetToothSelect.m`
- 目标：对比 periodic wide、selected subset、final periodic refine 的选齿与同齿细化分工。

## Snapshot index

| snapshot | 状态 | 配置摘要 | 备注 |
|---|---|---|---|
| `test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat` | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=100`，`checkpointEnable=true` | 同配置重跑结果；新增 warning 记录字段，当前代表性结果 |

## 当前代表性结果

### 配置

- snapshot：`test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat`
- `snrDb=10`
- `baseSeed=253`
- `numRepeat=100`
- `toothStepHz=750`
- checkpoint 完成 `100/100` 个任务
- snapshot 保存变量：`replayData`
- 关键新增字段：`replayData.warningTable`、`replayData.repeatTable.warningSeen`、`replayData.repeatTable.warningId`、`replayData.repeatTable.warningMessage`

### 主要统计

| 阶段 | tooth=0 命中率 | abs(toothIdx)<=1 比例 | abs(toothIdx) median / p95 / max | angle error median / p95 | abs(fdRefErrHz) median / p95 |
|---|---:|---:|---:|---:|---:|
| periodic wide | 42/100 | 51/100 | 1 / 228.55 / 305 | 0.005657 deg / 0.012967 deg | 750.018 Hz / 171.413 kHz |
| selected subset | 76/100 | 87/100 | 0 / 19.15 / 47 | 0.011009 deg / 0.013998 deg | 0.0327 Hz / 14.363 kHz |
| final periodic refine | 75/100 | 87/100 | 0 / 19.15 / 47 | 0.001542 deg / 0.003465 deg | 1.983 Hz / 14.360 kHz |

补充统计：

- periodic wide 有 `58/100` 个样本不在 truth tooth；selected subset / final periodic refine 可把其中 `38/58` 个样本拉回 `toothIdx=0`。
- final periodic refine 与 selected subset 的 tooth 在 `99/100` 个样本中相同，只在 seed `297` 从 `0` 变到 `-1`。
- final periodic refine 在 `91/100` 个样本中降低了 angle error；final / subset 的 angle error 中位比值约为 `0.152`。
- selected subset 相比 periodic wide 的 `abs(toothIdx)`：`56/100` 更小，`39/100` 相同，`5/100` 更大。
- final periodic refine 相比 selected subset 的 `abs(toothIdx)`：`99/100` 相同，`1/100` 更大，没有样本跨齿改善。
- selected subset / final periodic refine 后仍有 `13/100` 个样本超出 `abs(toothIdx)<=1`，其中 `9/100` 超过 `5` 个 tooth，`5/100` 超过 `20` 个 tooth；但没有超过 `50` 个 tooth 的远端错误。
- 该数据中 `usedPolish=0`，因此这组数据只能说明当前流程没有误触发 polish，不能证明 conditional polish 触发后是否有效。

### warning 记录

- 本次 `100` 个 repeat 中有 `1` 个 repeat 记录到非致命 MATLAB warning。
- warning repeat：seed `329`。
- seed `329` 的 final periodic refine 结果为 `toothIdx=-1`，`angleErrDeg=0.001833`，`fdRefErrHz=-750.083`。
- 该 warning seed 的角度误差不属于 tail，但它仍处在 one-tooth error 上；因此不应因为 warning 自动丢弃该 repeat，也不应把它当成已修复样本。
- 当前处理口径：只要 repeat 正常返回结果，warning 作为诊断元数据进入 `warningTable` / `repeatTable`，aggregate 统计不剔除该样本。若后续 warning 在错齿 tail 中集中出现，再单独检查 solver 条件数、fixed-DoA warm-anchor 或候选 adoption，而不是先改主 objective。

### 可观察现象

- periodic full-data 可以给出很小的 angle error，但会落到很远的 wrong tooth；最大达到 `305` 个 tooth，对应 `fdRef` 偏差约 `228.75 kHz`。
- 非周期 selected subset 的核心价值是选齿：它把 exact tooth 命中率从 `42%` 提高到 `76%`，把 `abs(toothIdx)<=1` 的比例从 `51%` 提高到 `87%`，并消除了 `>|50|` tooth 的远端错误。
- final periodic refine 的核心价值是同齿角度细化：它显著降低 angle error，但几乎不改变 selected subset 的 tooth。因此它不能被当作 wrong-tooth rescue。
- 若后续允许“一个周期内可接受”，`abs(toothIdx)<=1` 应作为辅助观察口径；但仍有 `13/100` 超出一个 tooth，需要 conditional rescue / subset bank 覆盖继续处理。

### 当前结论

- 当前 flow 的合理分工仍应保持为：subset 负责 tooth selection，periodic full-data 负责同齿 refine。
- 不应把 final periodic refine 期望成自动修齿步骤；进入 final refine 前必须先有可信的 subset tooth。
- 本次只有 `1/100` 个非致命 warning，且已保存到轻量诊断字段；warning 不影响本次机制结论，但后续长 MC 应继续保留 warning summary。
- 下一步更值得优先改的是 subset rescue / random bank 的触发条件，以及 same-tooth 内的 DoA-local-state refine，而不是继续缩小 periodic fd 搜索盒或把 alias-aware 字段灌回主 objective。

## 恢复方式

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后运行 `replayMfPeriodicVsSubsetToothSelect.m` 的 `Summary output and plotting` section。
