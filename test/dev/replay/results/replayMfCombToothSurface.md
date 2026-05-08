# replayMfCombToothSurface 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / mechanism-visualization` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfCombToothSurface_20260508-141941.mat` |
| 当前一句话结论 | 在 MS 周期窗口中，`1/T_f=750 Hz` comb 是 objective 层真实结构；CP-K / known-rate 的 DoA 局部信息没有退化，主要错误是 `fdRef` 落在相邻 tooth；CP-U / unknown-rate 同时出现 DoA 偏移与 `fdRef/fdRate` coupling，不能与 CP-K wrong-tooth 混为一谈。 |
| 决策影响 | 保留为 comb / wrong-tooth / DoA-Doppler coupling 机制可视化；可作为 appendix / diagnostic 图候选，不作为 estimator 或 flow 默认改动依据。 |
| 下一步动作 | 当前 surface 证据已足够支撑 CP-K wrong-tooth 机制解释；若继续推进论文主图，应转向 in-tooth / CRB-local scan 与 known/unknown-rate 信息损失，而不是继续扩展本 replay。 |
| 禁止误用 | 不能把 surface 最小点当作 runtime truth oracle；不能将 ridge slope、surface-min fd、truth tooth 或 `finalByMode` 解析结果下沉到默认 refine / selector / adoption；不能据此恢复 truth override、alias-aware objective injection 或 blanket DoA release。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfCombToothSurface.m`
- 结果文档：`test/dev/replay/results/replayMfCombToothSurface.md`
- replay 类型：机制可视化 / fixed-seed surface replay
- 主要问题：固定一个代表性 MS hard sample，观察 `fdRef` comb、DoA-fdRef surface、truth 附近 coupling，以及 CP-K / CP-U 在 wrong-tooth 下的差异。
- 观察范围：单个 seed、固定 MS 周期窗口、`unknown` / `known` 两种 probe model、truth center 与 `finalByMode` center 的 line / surface / coupling 切片。
- 不覆盖范围：不验证 flow 默认策略，不做 Monte Carlo 统计，不证明 CRB consistency，不决定 candidate bank 是否进入默认路径。
- truth 使用口径：truth 只用于结果评价、truth-centered line / surface / coupling center 与 mechanism label；不能进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `1/T_f` comb | 多帧等间隔观测下，`fdRef` 方向以 `1/frameIntvlSec` 为步长出现的周期性 objective 支路。本 snapshot 中 `T_f=1/750 s`，对应 `750 Hz`。 | Evaluation only | 只改变 line / surface 观察坐标。 | 说明 wrong-tooth 是 objective 层真实结构，不是 summary 口径假象。 |
| `unknown` / CP-U | `fdRateMode='unknown'`，参考星 Doppler rate 作为 nuisance 参数释放。 | No | 使用 unknown-rate probe model。 | 用于观察 unknown-rate release 后的 DoA-`fdRef`-`fdRate` coupling / bad basin。 |
| `known` / CP-K | `fdRateMode='known'`，固定参考星 Doppler rate。 | No | 使用 known-rate probe model。 | 用于分离 known-rate 条件下的 `fdRef` comb tooth 与 DoA 局部信息。 |
| `truth center` | 以真值 DoA 与真值 `fdRef` 为中心做 line / surface scan。 | Oracle only | 只用于机制上限和正确 tooth 参照。 | 可判断 truth 附近是否仍是 objective peak；不能作为 runtime center。 |
| `finalByMode` | mode-aware final center：`unknown` 解析为 `finalCpU`，`known` 解析为 `finalCpK`。 | No | 只改变诊断中心。 | 避免把 CP-U final 与 CP-K final 混用。 |
| `centerFd` | local lat/lon surface 中固定 `fdOffset=0`，即固定在当前 center 的 `fdRef`。 | No | 只改变 lat/lon 小面固定的 `fdRef`。 | 对 `finalCpK` 来说，这是 wrong-tooth estimated fdRef。 |
| `surfaceMinFd` | local lat/lon surface 中固定在 DoA-fdRef surface 找到的最优 `fdOffset`。 | Surface evaluation only | 只改变 lat/lon 小面固定的 `fdRef`。 | 对本 snapshot 的 `known + finalCpK` 来说，是 `fdOffset=+750 Hz`，即回到 truth tooth。 |
| `fdRef line scan` | 固定 DoA center，仅沿 `fdRef` 方向扫描 objective。 | center 可用 truth / final | 只观察一维 `fdRef` comb。 | 能说明 tooth 周期性，但不能解释 DoA 与 `fdRef` 的联合补偿。 |
| `DoA-fdRef surface` | 在一个 DoA offset 与 `fdRef` offset 网格上扫描 objective。 | center 可用 truth / final | 同时扰动 DoA 与 `fdRef`。 | 用来观察 final 附近是否存在联合移动后更好的局部点。 |
| `slice-scale summary` | 分别统计 DoA 方向切片跨度和 `fdRef` 方向切片跨度。 | No | 只增加诊断表。 | 用于判断 heatmap 中 DoA 方向是否被 `fdRef` 大尺度变化压扁。 |
| `local lat/lon surface` | 固定一个 `fdRef` offset 后，在局部 lat/lon 二维网格上扫描 objective。 | No | 只观察 DoA 局部面。 | 用于判断 CP-K 的 DoA 信息是否真的退化，还是被 `fdRef` comb 色阶掩盖。 |
| `truth-centered coupling` | 在 truth 附近小窗口内估计 DoA-fdRef 交叉耦合和 ridge 方向。 | Oracle only | 只改评价切片，不改 estimator。 | 只能作为 profile objective 机制证据；严格 FIM / EFIM 结论仍需直接计算。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfCombToothSurface_20260508-141941.mat` | 2026-05-08 | `representative` | `snrDb=10`，`baseSeed=253`，MS 周期窗口，`modeTagList=[unknown, known]`，`aliasStepHz=750 Hz`，line scan `±1500 Hz`，DoA-fdRef surface `61 x 121`，CP-K local lat/lon surface 同时固定 `centerFd` 与 `surfaceMinFd`。 | CP-K 的 DoA 局部曲率健康，`finalCpK` 主要错在 `fdRef=-1 tooth`；CP-U 的 `finalCpU` 同时存在 DoA 偏移与 `fdRef/fdRate` coupling。 | 当前代表性机制结果，取代 20260425 单模式旧结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `modeTagList = ["unknown", "known"]`
- `aliasStepHz = 750 Hz`
- `scanHalfWidthHzResolved = 1500 Hz`
- `selectedSubsetLabel = curated1`
- `selectedFinalTag = periodic-static-seed`
- line scan：`numFdGrid = 301`，覆盖 `±1500 Hz`
- DoA-fdRef surface：`surfaceDoaGridCount = 61`，`surfaceFdGridCount = 121`，DoA 半宽 `0.004 deg`
- local lat/lon surface：`latLonSurfaceGridCount = 41`，`latLonSurfaceHalfWidthDeg = 0.004 deg`，仅对 `known + finalByMode` 跑 `centerFd` 与 `surfaceMinFd`
- truth-centered coupling：`couplingDoaGridCount = 101`，`couplingFdGridCount = 101`
- snapshot 保存变量：`replayData`

### 4.2 center summary

| center | latDeg | lonDeg | fdRefHz | fdRateHzPerSec | 解释 |
|---|---:|---:|---:|---:|---|
| `truth` | `37.78` | `36.59` | `-28940` | `-3833.5` | 真值中心。 |
| `staticSeed` | `37.778` | `36.588` | `-29184` | `-3833.5` | static seed，`fdRate` 固定为真值 / known-rate 口径。 |
| `finalCpU` | `37.778` | `36.588` | `-29690` | `-3851.6` | unknown-rate final，DoA 与 `fdRef/fdRate` 均偏离 truth。 |
| `finalCpK` | `37.78` | `36.59` | `-29690` | `-3833.5` | known-rate final，DoA 基本贴近 truth，但 `fdRef` 比 truth 偏约 `-750 Hz`。 |

`finalCpK.fdRef - truth.fdRef ≈ -750 Hz`，正好对应一个 comb tooth；因此 CP-K 的首要现象是 wrong-tooth，而不是 DoA 失效。

### 4.3 fdRef comb line summary

| mode | center | resolved center | min fd offset | min alias index | center `deltaObj` | min objective | 解释 |
|---|---|---|---:|---:|---:|---:|---|
| `unknown` | `truth` | `truth` | `0 Hz` | `0` | `0` | `-1.9648e6` | truth tooth 是 truth center 下的最优齿。 |
| `unknown` | `finalByMode` | `finalCpU` | `750 Hz` | `1` | `236.82` | `-1.0307e6` | CP-U final 位于相邻 tooth，沿 `fdRef` 加回 `+750 Hz` 更优。 |
| `known` | `truth` | `truth` | `0 Hz` | `0` | `0` | `-1.9648e6` | CP-K truth center 正常。 |
| `known` | `finalByMode` | `finalCpK` | `750 Hz` | `1` | `470.6` | `-1.9648e6` | CP-K final 的 DoA 已对，但 `fdRef` 在 wrong tooth；correct tooth 只比 wrong tooth 好约 `470.6`。 |

### 4.4 DoA-fdRef surface summary

| mode | center | resolved center | min DoA offset | min fd offset | min alias index | center `deltaObj` | 解释 |
|---|---|---|---:|---:|---:|---:|---|
| `unknown` | `truth` | `truth` | `0 deg` | `0 Hz` | `0` | `0` | truth surface 最优点在 truth。 |
| `unknown` | `finalByMode` | `finalCpU` | `-0.0024 deg` | `750 Hz` | `1` | `9.3245e5` | CP-U final 不是局部 surface 好点，需要同时修 DoA 与 `fdRef`。 |
| `known` | `truth` | `truth` | `0 deg` | `0 Hz` | `0` | `0` | CP-K truth surface 正常。 |
| `known` | `finalByMode` | `finalCpK` | `-0.00013333 deg` | `750 Hz` | `1` | `470.76` | CP-K 只需把 `fdRef` 加回一个 tooth，DoA 几乎不动。 |

这里最重要的对比是 `centerDeltaObj`：CP-U final 离 surface minimum 差 `9.3245e5`，而 CP-K final 只差 `470.76`。因此 CP-U 是更强的 coupling / bad basin 问题，CP-K 是 near-equivalent wrong-tooth 问题。

### 4.5 DoA-fdRef slice-scale summary

| mode | center | resolved center | min DoA offset | min fd offset | bestFdDoaSpan | centerFdDoaSpan | bestDoaFdSpan | centerDoaFdSpan | best ratio | center ratio |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `unknown` | `truth` | `truth` | `0 deg` | `0 Hz` | `9.3648e5` | `9.3648e5` | `1.8819e6` | `1.8819e6` | `2.0095` | `2.0095` |
| `unknown` | `finalByMode` | `finalCpU` | `-0.0024 deg` | `750 Hz` | `9.38e5` | `9.3776e5` | `1.8793e6` | `9.4399e5` | `2.0035` | `1.0066` |
| `known` | `truth` | `truth` | `0 deg` | `0 Hz` | `9.3648e5` | `9.3648e5` | `1.8819e6` | `1.8819e6` | `2.0095` | `2.0095` |
| `known` | `finalByMode` | `finalCpK` | `-0.00013333 deg` | `750 Hz` | `507.49` | `507.5` | `1.8819e6` | `1.8819e6` | `3708.3` | `3708.2` |

`known + finalCpK` 的 `fdRef` 方向跨度约为 DoA 方向跨度的 `3708` 倍。因此普通 heatmap 如果不做 clipped color scale，DoA 曲率会被 `fdRef` 方向的大动态范围压扁，看起来像“known-rate objective 和 DoA 无关”。这个视觉现象不等价于 DoA 信息退化。

### 4.6 local lat/lon surface summary

| mode | center | resolved center | fixed fd tag | fixed fd offset | alias index | min lat offset | min lon offset | radial offset | center `deltaObj` | lon-slice span | lat-slice span | 解释 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `known` | `finalByMode` | `finalCpK` | `centerFd` | `0 Hz` | `0` | `0 deg` | `0 deg` | `0 deg` | `0` | `9.3622e5` | `9.3847e5` | 固定在 CP-K 估计出来的 wrong-tooth `fdRef` 上，DoA 小面仍有强曲率，且 center DoA 是该小面最优。 |
| `known` | `finalByMode` | `finalCpK` | `surfaceMinFd` | `750 Hz` | `1` | `0 deg` | `0 deg` | `0 deg` | `0` | `9.3645e5` | `9.3871e5` | 固定在 correct tooth 上，DoA 小面同样有强曲率，center DoA 仍是该小面最优。 |

这两张 lat/lon 小面共同说明：CP-K 的 DoA 局部信息没有坏；wrong tooth 和 correct tooth 上的 DoA 曲率都很明显。`centerDeltaObj=0` 只表示各自固定 `fdRef` 的 lat/lon 小面内部 center DoA 是最优，不表示两个 tooth 的全局 objective 一样。tooth 之间的差异应看 DoA-fdRef surface 中的 `centerDeltaObj=470.76`。

### 4.7 truth-centered coupling summary

| mode | center | `rhoDoaFd` | ridge slope | quadratic ridge slope | quadratic fit RMSE | min DoA offset | min fd offset |
|---|---|---:|---:|---:|---:|---:|---:|
| `unknown` | `truthCoupling` | `0.49337` | `-48.923 Hz/deg` | `-12513 Hz/deg` | `2.2415e5` | `0 deg` | `0 Hz` |
| `known` | `truthCoupling` | `NaN` | `0 Hz/deg` | `-20.011 Hz/deg` | `3.2531e5` | `0.00012 deg` | `0 Hz` |

`unknown` 下 truth 附近存在可见 DoA-`fdRef` 交叉耦合。`known` 下 coupling 指标退化，不应解释成 DoA 无信息；结合 lat/lon 小面看，更合理的解释是 CP-K 已基本固定 rate 且 DoA center 正确，truth 附近的 1-D coupling 指标对当前切片不敏感。

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `fdRef` line scan 在 truth center 和 `finalByMode` center 都显示 `750 Hz` comb；这支持 wrong-tooth / near-tooth 是 objective 层结构。
- `known + finalCpK` 的 DoA-fdRef surface minimum 位于 `fdOffset=+750 Hz`、`doaOffset=-0.00013333 deg`；说明 CP-K 主要是 `fdRef` 错一个 tooth，DoA 基本正确。
- `known + finalCpK` 的 slice-scale ratio 约 `3708`；说明 heatmap 中 DoA 方向被 `fdRef` 大动态范围掩盖，而不是 DoA 信息消失。
- `known + finalCpK` 的 `centerFd` 与 `surfaceMinFd` 两张 lat/lon surface 都有约 `9.36e5` 量级 DoA 切片跨度，证明 CP-K local DoA 曲率健康。
- `unknown + finalCpU` 的 surface minimum 需要 `doaOffset=-0.0024 deg` 与 `fdOffset=+750 Hz`，且 `centerDeltaObj=9.3245e5`；说明 CP-U 的问题比 CP-K 更像 DoA-`fdRef`-`fdRate` coupling / bad basin。

### 5.2 仍未解决或反向的现象

- 该 replay 是单 seed 机制可视化，不能说明这种 surface 形态在大 MC 中的覆盖率。
- CP-K 的 correct tooth 与 wrong tooth objective gap 只有几百量级，说明 full-flow / acquisition 仍可能选择 wrong tooth；但这不应被写成局部 MLE / CRB 主目标失败。
- CP-U 的 bad basin 需要在 in-tooth / CRB-local 口径下另行隔离，不能用本 replay 直接判断 CP-U 理论模型失败。
- ridge slope 数值不稳定，尤其 known-rate coupling 指标出现退化；当前 replay 不适合给出可下沉策略的 slope 参数。
- line / surface 结果只能解释 local objective 结构，不能自动给出 safe gate 或 adoption rule。

### 5.3 代表性 case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `253 / CP-K` | known-rate wrong-tooth sample | `finalCpK` DoA 贴近 truth，但 `fdRef=-29690 Hz`，比 truth `-28940 Hz` 偏一个 `-750 Hz` tooth；lat/lon 小面显示 DoA 曲率健康。 | 证明 CP-K 的主要问题是 comb tooth ambiguity，不是 DoA objective 退化。 |
| `253 / CP-U` | unknown-rate coupling sample | `finalCpU` 同时有 DoA 偏移、`fdRef=-29690 Hz` 与 `fdRate=-3851.6 Hz/s` 偏移；surface minimum 在 `doaOffset=-0.0024 deg, fdOffset=+750 Hz`。 | 说明 CP-U 的问题需要按 unknown-rate nuisance coupling / bad basin 另行解释。 |

## 6. 机制解释

### 6.1 当前解释

本 replay 把当前 bad sample 拆成 CP-K 与 CP-U 两层来看。

对 CP-K / known-rate，`finalCpK` 的 DoA 已经几乎等于 truth，但 `fdRef` 落在相邻 `-1 tooth`。DoA-fdRef surface 显示，只要把 `fdRef` 加回 `+750 Hz`，几乎不需要移动 DoA 就能回到正确 tooth；local lat/lon surface 又显示，无论固定 wrong tooth 还是 correct tooth，DoA 小面都有明显曲率。因此，CP-K 的结果应解释为 **near-equivalent fdRef comb branch / acquisition ambiguity**，而不是 known-rate 局部 DoA 信息损失。

对 CP-U / unknown-rate，`finalCpU` 不仅落在相邻 tooth，还伴随 DoA 与 `fdRate` 偏移；surface minimum 与 final center 的 objective 差达到 `9.3245e5`，远大于 CP-K wrong tooth 的 `470.76`。因此 CP-U 更像是 unknown nuisance-rate release 后的 DoA-`fdRef`-`fdRate` coupling / bad basin，需要在 in-tooth / CRB-local 统计中隔离，而不应直接与 CP-K wrong-tooth 混写。

### 6.2 这个结果支持什么

- 支持 `1/T_f` comb 是 objective 层真实结构。
- 支持 CP-K wrong-tooth 是 acquisition / comb ambiguity，而不是 DoA 信息消失。
- 支持 CP-U 的坏点需要从 unknown-rate nuisance coupling / bad basin 角度解释。
- 支持把 `replayMfCombToothSurface` 保留为 mechanism visualization / diagnostic 图候选。
- 支持论文正文主精度验证优先使用 Doppler-aided / in-tooth local estimation 口径；full-flow wrong-tooth 作为工程边界或附录机制讨论。

### 6.3 这个结果不证明什么

- 不证明 estimator 默认路径已修复。
- 不证明 CP-K 能自动完成 global tooth acquisition。
- 不证明 CP-U 理论模型失败；当前看到的是单 seed full-flow / probe-level bad basin。
- 不证明可以把 `surfaceMinFd`、ridge slope、quadratic ridge slope 或 nonzero `fdRef` offset 下沉为默认 refine 策略。
- 不证明 joint fd bank 应进入默认 candidate family；后续 envelope / gated replay 已显示默认收益主要来自 `wide / single-MF` DoA basin-entry。
- 不证明可以使用 truth center、truth tooth 或 truth DoA 做 runtime selector。
- 不证明 FIM / EFIM 交叉项的严格数值；该结论仍需专门计算。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉 ridge slope、truth override、surface-min fd 或 alias-aware objective 字段。 |
| flow 默认路径 | 不直接改；仅作为后续 basin-entry / wrong-tooth / same-tooth rescue 的机制背景。 |
| regression | 不写；单 seed surface 不是自动 pass/fail 契约。 |
| replay / scan 下一步 | 本 replay 暂停继续扩图；下一步更应推进 MS in-tooth / CRB-local scan 与 known/unknown-rate 信息损失统计。 |
| 论文图 / 论文口径 | 可作为 appendix / diagnostic-only 机制图候选；正文主图仍优先 CRB consistency、resolved RMSE、CP/IP 与 known/unknown-rate 信息损失。 |
| 排障记录 | 若只整理本结果，不必新增主记录；若未来作为论文机制图使用，可在机制归并版补一条“CP-K DoA 未退化，wrong-tooth 属于 comb acquisition ambiguity”。 |

## 8. 限制与禁止解释

- 不要把本 replay 的离线 surface minimum 当作 runtime 可用 candidate。
- 不要把 truth-centered coupling 的正结果迁移到 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 ridge slope 数值公式化为默认路径参数。
- 不要把 `surfaceMinFd` 直接作为 estimator 的默认 tooth correction。
- 不要用单 seed surface 替代 Monte Carlo 指标、CRB / EFIM 结论或 flow-like 验证。
- 不要恢复 truth override、alias-aware objective injection、blanket DoA release 或 blanket joint fd bank。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfCombToothSurface_20260508-141941.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfCombToothSurface.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本文档只保留 `20260508-141941` 作为当前代表性机制结果。
- `20260425-204458` 旧结果已被本次 CP-U / CP-K 双模式、slice-scale 与 dual lat/lon surface 结果取代；旧结果仍可作为最早 `1/T_f` comb 证据回溯，但不再作为当前 representative snapshot。
