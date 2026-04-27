# replayMfInToothFdRangeOracle 结果记录

## 对应 replay

- 脚本：`test/dev/replay/replayMfInToothFdRangeOracle.m`
- 目标：在 truth-centered 半齿 `fdRef` oracle 条件下，对比单星 / 多星、单帧 / 多帧、wide / in-tooth、known-rate / unknown-rate 的性能上限，并定位 `MS-MF-CP-U-in-tooth` 是否已经具备支撑论文主张的上限优势。

## Snapshot index

| snapshot | 状态 | 配置摘要 | 备注 |
|---|---|---|---|
| `test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat` | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=50`，`fdOracleHalfToothFraction=0.49`，`fdRateOracleHalfWidthHzPerSec=1000` | 当前保留的 half-tooth oracle / paper-claim upper-bound / same-tooth tail 代表性结果 |

## 2026-04-26 10:42 结果分析

### 配置与总体状态

- snapshot：`test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat`
- 配置：`snrDb=10`，`baseSeed=253`，`numRepeat=50`
- `toothStepHz=750`
- `fdRef` oracle 范围：truth-centered `0.49` tooth，即约 `±367.5 Hz`
- `fdRate` oracle 范围：truth-centered `±1000 Hz/s`
- replay 定位：论文主张上限 replay，只用于机制判断，不进入默认 flow / regression
- 当前一句话结论：`MS-MF-CP-U-in-tooth` 已经显著优于单帧 static 与 wide baseline，且在多数 seed / median 口径上优于 `SS-MF-CP-U-in-tooth`；但 RMSE 仍被少数 same-tooth multi-sat coherence hard case 拉坏，因此还不能写成“多星多帧稳定最好”。

### Oracle method aggregate

| method | tooth hit | angle RMSE (deg) | angle median (deg) | angle p95 (deg) | `abs(fdRef)` median / p95 (Hz) | `abs(fdRate)` median / p95 (Hz/s) | wall time median (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 1.00 | 0.002153 | 0.001691 | 0.003610 | 130.339 / 272.776 | 3833.5 / 3833.5 | — |
| `ms-sf-static` | 0.98 | 0.001871 | 0.001466 | 0.003340 | 84.315 / 271.973 | 3833.5 / 3833.5 | — |
| `ss-mf-cp-u-wide` | 0.66 | 0.003945 | 0.000838 | 0.010985 | 122.074 / 1500.012 | 12.262 / 3833.488 | 2109 |
| `ms-mf-cp-u-wide` | 0.46 | 0.007887 | 0.005122 | 0.012580 | 747.065 / 171000.014 | 6.736 / 3333.497 | 24648 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.000976 | 0.000691 | 0.001488 | 0.00949 / 0.03389 | 7.724 / 19.338 | 1797 |
| `ms-mf-cp-u-in-tooth` | 1.00 | 0.001407 | 0.000416 | 0.003607 | 0.00989 / 0.52708 | 4.306 / 14.295 | 46476 |
| `ms-mf-cp-k-in-tooth` | 1.00 | 0.002127 | 0.001358 | 0.004084 | 0.01346 / 0.41224 | 0 / 0 | 4190 |
| `ms-mf-cp-u-truth-doa-oracle` | 1.00 | 0.000493 | 0.000343 | 0.000873 | 0.00975 / 0.03244 | 4.504 / 12.886 | 41638 |

### Paper-claim upper-bound compare

`paperCompareTable` 以 `ms-mf-cp-u-in-tooth` 为 target。正的 `angleRmseGainDeg` 表示 target 优于 baseline。

| target | baseline | angle RMSE gain (deg) | angle median gain (deg) | `abs(fdRef)` median gain (Hz) | target better angle rate |
|---|---|---:|---:|---:|---:|
| `ms-mf-cp-u-in-tooth` | `ss-sf-static` | +0.000745 | +0.001275 | +130.329 | 0.86 |
| `ms-mf-cp-u-in-tooth` | `ms-sf-static` | +0.000463 | +0.001050 | +84.305 | 0.88 |
| `ms-mf-cp-u-in-tooth` | `ss-mf-cp-u-in-tooth` | -0.000431 | +0.000275 | -0.000402 | 0.76 |
| `ms-mf-cp-u-in-tooth` | `ms-mf-cp-u-wide` | +0.006479 | +0.004706 | +747.055 | 0.82 |
| `ms-mf-cp-u-in-tooth` | `ms-mf-cp-k-in-tooth` | +0.000719 | +0.000941 | +0.003566 | 0.72 |
| `ms-mf-cp-u-in-tooth` | `ms-mf-cp-u-truth-doa-oracle` | -0.000915 | -0.000073 | -0.000138 | 0.32 |

### Single-vs-multi in-tooth tail summary

- `singleMultiCompareTable` 共 `50` 行，逐 seed 记录 `SS-MF-CP-U-in-tooth` 与 `MS-MF-CP-U-in-tooth` 的角度、频率和 coherence 差异。
- `MS-MF-CP-U-in-tooth` 在 `38/50` 个 seed 中角度误差小于 `SS-MF-CP-U-in-tooth`，即 better angle rate 为 `0.76`。
- `MS-MF-CP-U-in-tooth` 的 angle median 更小：`0.000416 deg` vs `0.000691 deg`。
- 但 `MS-MF-CP-U-in-tooth` 的 angle RMSE 更大：`0.001407 deg` vs `0.000976 deg`，说明少数 tail case 主导 RMSE。

Worst tail cases：

| rank | seed | single angle (deg) | multi angle (deg) | truth-DoA angle (deg) | gain = single - multi (deg) | multi gap to truth-DoA (deg) | multi coherence floor | 解释 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 277 | 0.000556 | 0.004735 | 0.000033 | -0.004179 | 0.004702 | 0.0894 | 非参考星 coherence 崩，典型 same-tooth hard case |
| 2 | 283 | 0.000528 | 0.003618 | 0.000519 | -0.003090 | 0.003099 | 0.0995 | 非参考星 coherence 崩 |
| 3 | 298 | 0.001311 | 0.003594 | 0.000158 | -0.002283 | 0.003437 | 0.1006 | 非参考星 coherence 崩 |
| 4 | 256 | 0.003420 | 0.005454 | 0.000021 | -0.002034 | 0.005432 | 0.00066 | 最强 hard case，非参考星几乎完全没接住 |
| 5 | 268 | 0.001262 | 0.001960 | 0.000307 | -0.000698 | 0.001653 | 0.9127 | 部分 coherence / fd residual 问题 |
| 6 | 284 | 0.000646 | 0.000891 | 0.000381 | -0.000245 | 0.000510 | 0.9938 | 同齿内小 tail |
| 7 | 293 | 0.000468 | 0.000674 | 0.000650 | -0.000206 | 0.000024 | 1.0000 | 轻微波动，不是主要 hard case |
| 8 | 280 | 0.000080 | 0.000249 | 0.000363 | -0.000169 | -0.000115 | 1.0000 | 轻微波动，truth-DoA oracle 也不占优 |

### Runtime timing summary

| stage | median (ms) | p95 (ms) | total (s) | median repeat share |
|---|---:|---:|---:|---:|
| build repeat data | 1900 | 2306 | 96.7 | — |
| static transition bundle | 3414 | 4956 | 177.3 | 0.0265 |
| dynamic methods total | 122712 | 150208 | 6316.7 | 0.9594 |
| slowest dynamic method | 48778 | 67147 | 2607.3 | 0.3774 |
| repeat total | 127665 | 156421 | 6591.4 | 1.0000 |

- static transition bundle 只占 repeat 中位时间约 `2.65%`，不是瓶颈。
- dynamic methods 占 repeat 中位时间约 `95.94%`，主耗时来自 dynamic estimator 本身。
- 最重的 method 是 `ms-mf-cp-u-in-tooth` 与 `ms-mf-cp-u-truth-doa-oracle`。
- 后续若需要降运行时间，应优先拆 method bank 或分两轮 replay，而不是优化 bundle。

### 可观察现象

- `MS-MF-CP-U-wide` 的 tooth hit rate 只有 `0.46`，说明非 oracle flow 仍有明显 wrong-tooth 风险。
- `MS-MF-CP-U-in-tooth` 的 tooth hit rate 为 `1.00`，`fdRef` median error 约 `0.0099 Hz`，`fdRate` median error 约 `4.31 Hz/s`，说明正确 tooth 内频率链是健康的。
- `MS-MF-CP-U-in-tooth` 明显优于 `SS-SF-Static` 与 `MS-SF-Static`，说明多帧 continuous-phase dynamic 模型在正确 tooth 内相对单帧 static 有实际增益。
- `MS-MF-CP-U-in-tooth` 明显优于 `MS-MF-CP-U-wide`，说明 half-tooth oracle 有效拆掉了 wide baseline 中的 comb / wrong-tooth 问题。
- `MS-MF-CP-U-in-tooth` 在 median 与多数 seed 上优于 `SS-MF-CP-U-in-tooth`，但 RMSE 与 p95 仍被少数 same-tooth tail 拉坏。
- `MS-MF-CP-U-truth-doa-oracle` 明显优于所有非 truth-DoA 方法，说明多星多帧的信息上限存在；当前差距主要来自 DoA / non-ref coherence basin 承接不足，而不是多星信息本身无价值。

### 当前结论

- 该 replay 支持继续推进论文主线，但不能把当前结果写成“多星多帧已经稳定最好”。
- 更准确的阶段性结论是：在正确 tooth 内，`MS-MF-CP-U` 的频率链健康，并相对单帧 static 与 wide baseline 有明确增益；但当前实现仍存在少数 same-tooth non-ref coherence hard case，导致 RMSE 暂时不能稳定压过 `SS-MF-CP-U-in-tooth`。
- `truth-DoA oracle` 结果表明，多星多帧的上限仍然明显优于单星多帧；因此当前主问题应判为 same-tooth hard-case refinement / flow adoption，而不是模型上限不足。
- `fdRate` unknown-rate 分支不是当前拖累项；`MS-MF-CP-U-in-tooth` 优于 `MS-MF-CP-K-in-tooth`。
- timing 结果确认 bundle 不是运行瓶颈；后续性能优化优先拆 method bank，不优先动 static transition bundle。

### 对主流程的影响

- 仍需要运行 / 增强已有 `scanMfSubsetBankCoverage`，因为 wide baseline 的 wrong-tooth 风险仍明显；但它只能解决 tooth selection，不能解决 same-tooth tail，不再另建 `scanMfCuratedSubsetScheduleSearch`。
- 下一步应先用固定 hard seeds 验证 gated joint local refine：建议 seed `277, 283, 298, 256, 268`。
- replay 对照建议至少包含 disabled / DoA-only polish / joint DoA-fdRef local refine，并用 no-truth 指标触发，例如 `nonRefCoherenceFloor`、multi-sat block coherence、`fd` 已健康但 final objective 异常等。
- 只有 same-tooth tail 被压住后，`MS-MF-CP-U-in-tooth` 的 majority / median 增益才可能转化为稳定 RMSE 增益。
- 不建议 blanket 放宽 DoA，不建议把 truth-DoA oracle 灌回默认 flow，不建议把 alias-aware 或 probe 字段灌回主 residual / objective。

## 恢复方式

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后运行 `replayMfInToothFdRangeOracle.m` 的 `Summary output and plotting` section。
