# replayMfCombToothSurface 结果记录

## 对应 replay

- 脚本：`test/dev/replay/replayMfCombToothSurface.m`
- 目标：固定样本可视化 `1/T_f` fdRef comb、DoA-fdRef surface，以及 truth 附近的 DoA-fdRef coupling。

## Snapshot index

| snapshot | 状态 | 配置摘要 | 备注 |
|---|---|---|---|
| `test/data/cache/replay/replayMfCombToothSurface_20260425-204458.mat` | representative | `snrDb=10`，`baseSeed=253`，`toothStepHz=750`，coupling 窗口 `±0.006 deg / ±300 Hz` | 当前保留的 comb / wrong-tooth / DoA-fdRef coupling 机制代表性结果 |

## 2026-04-25 20:44 结果分析

### 配置与总体状态

- snapshot：`test/data/cache/replay/replayMfCombToothSurface_20260425-204458.mat`
- 配置：`snrDb=10`，`baseSeed=253`，`frameIntvlSec=1/750 s`，`toothStepHz=750 Hz`
- line scan：`numFdGrid=301`，覆盖 `±1500 Hz`
- surface scan：`surfaceDoaGridCount=61`，`surfaceFdGridCount=121`，DoA 半宽 `0.004 deg`
- coupling scan：`couplingDoaGridCount=101`，`couplingFdGridCount=101`，窗口 `±0.006 deg / ±300 Hz`
- flow final：`angleErrDeg=0.002142`，`fdRefErrHz=603.221`，`fdRateErrHzPerSec=-15.246`，`toothIdx=1`，`toothResidualHz=-146.779`
- flow 选择：`selectedSubsetLabel=curated1`，`selectedFinalTag=periodic-static-seed`

### line scan 结果

| center | min fd offset | min alias index | center `deltaObj` | `deltaObj(-1500)` | `deltaObj(-750)` | `deltaObj(0)` | `deltaObj(+750)` | `deltaObj(+1500)` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| truth | 0 Hz | 0 | 0 | 1754.65 | 470.59 | 0 | 344.69 | 1503.33 |
| final estimate | -750 Hz | -1 | 171.16 | 241.33 | 0 | 171.16 | 754.09 | 1746.69 |

- truth center 处正确 tooth 仍是 line scan 全局最优；但相邻 `±750 Hz` tooth 的 objective 差只有几百量级，而普通非 alias 区域的 `deltaObj` 中位数约为 `1.829e6`。
- final-estimate center 的一维 fdRef line scan 最优 tooth 在 `-750 Hz`，当前 final center 只比该 tooth 差 `171.16`，说明 final 仍停在相邻 comb tooth 附近。
- `1/T_f = 750 Hz` 的 comb 结构仍然是 objective 层真实结构，不是 summary 口径导致的假象。

### DoA-fdRef surface 结果

| center | min DoA offset | min fd offset | min alias index | center `deltaObj` | min objective |
|---|---:|---:|---:|---:|---:|
| truth | 0 deg | 0 Hz | 0 | 0 | -1.9648e6 |
| final estimate | -0.0024 deg | -600 Hz | -1 | 8.7031e5 | -1.9003e6 |

- truth surface 的最优点仍在 `(doaOffset=0, fdOffset=0)`，确认正确 tooth 是 objective 层真实峰值。
- final-estimate surface 的最优点不是严格的 `-750 Hz`，而是 `fdOffset=-600 Hz` 同时配合 `doaOffset=-0.0024 deg`；该点把 `fdRef` 从 `-28337.09 Hz` 拉到约 `-28937.09 Hz`，已经接近 truth `-28940.31 Hz`。
- 因此当前 hard case 不是单纯一维 fdRef alias 问题，而是 DoA 与 fdRef 可以互相补偿的局部曲面问题。只靠 fdRef 缩盒或按齿硬切不够，后续更适合在 flow 层做带 gate 的 joint local rescue / polish。

### truth-centered DoA-fdRef coupling 结果

| center | `rhoDoaFd` | ridge slope | quadratic ridge slope | quadratic fit RMSE | min DoA offset | min fd offset |
|---|---:|---:|---:|---:|---:|---:|
| truth coupling | 0.4934 | -48.92 Hz/deg | -1.2513e4 Hz/deg | 2.2415e5 | 0 deg | 0 Hz |

- `rhoDoaFd≈0.49` 明显不接近 0，说明 truth 附近的 profile objective 在 DoA 与 fdRef 两个方向上存在可见交叉曲率。
- coupling surface 的最优点仍在 truth，说明本次 coupling 诊断不是在说 truth 被移走，而是在说明 truth 附近局部等高线不是严格轴对齐。
- ridge 提取的 `ridgeSlopeHzPerDeg` 与二次拟合的 `quadRidgeSlopeHzPerDeg` 差异较大，且 `quadFitRmse≈2.24e5`，说明 `±0.006 deg / ±300 Hz` 窗口内的 profile surface 不是非常理想的纯二次面；因此不应过度解释 slope 的具体数值。
- 该结果更稳妥的用法是：把 `rhoDoaFd` 与二维 surface 的联合最优偏移一起作为 objective-level 证据，说明 DoA 与 fdRef 在 profile MLE / EFIM 层不宜先验视为完全可分离。
- 严格的 FIM 交叉项结论仍应由 FIM 直接计算或更小窗口的局部二次拟合验证；本 replay 只作为机制可视化证据。

### 本次结论

- 本次结果确认 `1/T_f` comb 是 objective surface 的真实结构。
- 正确 tooth 在 truth center 处仍是最优，但相邻 tooth 的 objective gap 过小，足以吸引 flow 落到 wrong-tooth 或 near-tooth 局部盆地。
- final center 附近的小范围 DoA-fdRef surface 内已经包含接近 truth 的更好点，并且该更好点需要 DoA 与 fdRef 联合移动；这支持后续优先检查 flow 层是否能用条件触发的 joint local refine 接住该盆地。
- truth 附近 coupling 诊断给出 `rhoDoaFd≈0.49`，支持“DoA-Doppler 在 profile objective / EFIM 层存在可见耦合”的解释，但不应写成“该 replay 严格证明 FIM 交叉项非零”。
- 后续若要把该点写入论文或 CRB/EFIM 说明，建议再用更小 coupling 窗口复跑，或直接计算 FIM / EFIM 交叉块。
- 该结果不支持把 truth override、alias-aware 字段或 blanket DoA release 重新灌回 estimator 默认主路径。

## 恢复方式

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后运行 `replayMfCombToothSurface.m` 的 `Summary output and plotting` section。
