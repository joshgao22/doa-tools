# replayMfCombToothSurface 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / mechanism-visualization` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfCombToothSurface_20260425-204458.mat` |
| 当前一句话结论 | `1/T_f` comb 是 objective 层真实结构；final 附近存在 DoA-fdRef 联合补偿方向，但本 replay 只用于解释 surface / coupling，不支持把 ridge slope 或 joint fd offset 下沉默认路径。 |
| 决策影响 | 保留为 comb / wrong-tooth / DoA-Doppler coupling 机制可视化；可作为 appendix / diagnostic 图候选，不作为 estimator 或 flow 默认改动依据。 |
| 下一步动作 | 若需要论文图或更严格 EFIM 解释，可用更小 coupling 窗口复跑，或直接计算 FIM / EFIM 交叉块；当前 replay 不继续扩成策略搜索。 |
| 禁止误用 | 不能把 surface 最小点当作 runtime truth oracle；不能将 ridge slope 公式化进入默认 refine；不能据此恢复 truth override、alias-aware objective injection 或 blanket DoA release。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfCombToothSurface.m`
- 结果文档：`test/dev/replay/results/replayMfCombToothSurface.md`
- replay 类型：机制可视化 / fixed-seed surface replay
- 主要问题：固定一个代表性 hard sample，观察 `fdRef` comb、DoA-fdRef surface 与 truth 附近 coupling 是否真实存在。
- 观察范围：单个 seed、固定 flow final、truth center 与 final-estimate center 的 line / surface / coupling 切片。
- 不覆盖范围：不验证 flow 默认策略，不做 Monte Carlo 统计，不证明 CRB consistency，不决定 candidate bank 是否进入默认路径。
- truth 使用口径：truth 只用于结果评价、truth-centered line / surface / coupling center 与 mechanism label；不能进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `1/T_f` comb | 多帧等间隔观测下，`fdRef` 方向以 `1/frameIntvlSec` 为步长出现的周期性 objective 支路。本 snapshot 中 `T_f=1/750 s`，对应 `750 Hz`。 | Evaluation only | 只改变 line / surface 观察坐标。 | 说明 wrong-tooth 是 objective 层真实结构，不是 summary 口径假象。 |
| `truth center` | 以真值 DoA 与真值 `fdRef` 为中心做 line / surface scan。 | Oracle only | 只用于机制上限和正确 tooth 参照。 | 可判断 truth 附近是否仍是 objective peak；不能作为 runtime center。 |
| `final-estimate center` | 以当前 flow final estimate 为中心做 line / surface scan。 | No | 固定当前估计点附近的局部切片。 | 用来解释当前 final 为什么卡在 near-tooth / wrong-tooth 邻域。 |
| `fdRef line scan` | 固定 DoA center，仅沿 `fdRef` 方向扫描 objective。 | center 可用 truth / final | 只观察一维 `fdRef` comb。 | 能说明 tooth 周期性，但不能解释 DoA 与 `fdRef` 的联合补偿。 |
| `DoA-fdRef surface` | 在一个 DoA offset 与 `fdRef` offset 网格上扫描 objective。 | center 可用 truth / final | 同时扰动 DoA 与 `fdRef`。 | 用来观察 final 附近是否存在联合移动后更好的局部点。 |
| `truth-centered coupling` | 在 truth 附近小窗口内估计 DoA-fdRef 交叉耦合和 ridge 方向。 | Oracle only | 只改评价切片，不改 estimator。 | 只能作为 profile objective 机制证据；严格 FIM / EFIM 结论仍需直接计算。 |
| `ridge slope` | 从 surface / 二次拟合中估计的局部 ridge 方向。 | Evaluation only | 只用于解释曲面形状。 | 数值本身不稳定，不能直接变成默认 refine 斜率或策略参数。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfCombToothSurface_20260425-204458.mat` | 2026-04-25 | `representative` | `snrDb=10`，`baseSeed=253`，`frameIntvlSec=1/750 s`，`toothStepHz=750 Hz`，line scan `±1500 Hz`，surface `61 x 121`，coupling window `±0.006 deg / ±300 Hz` | comb 与 DoA-fdRef coupling 均可见；final 附近存在联合偏移更优点，但不支持把 ridge slope 下沉默认路径。 | 当前代表性机制结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `frameIntvlSec = 1/750 s`
- `toothStepHz = 750 Hz`
- line scan：`numFdGrid = 301`，覆盖 `±1500 Hz`
- DoA-fdRef surface：`surfaceDoaGridCount = 61`，`surfaceFdGridCount = 121`，DoA 半宽 `0.004 deg`
- truth-centered coupling：`couplingDoaGridCount = 101`，`couplingFdGridCount = 101`，窗口 `±0.006 deg / ±300 Hz`
- snapshot 保存变量：`replayData`

### 4.2 flow final summary

| 指标 | 数值 | 解释 |
|---|---:|---|
| `selectedSubsetLabel` | `curated1` | 当前 flow 选择的 subset。 |
| `selectedFinalTag` | `periodic-static-seed` | 当前最终解来源。 |
| `angleErrDeg` | `0.002142` | final DoA 已较小，但仍可能停在 near-tooth 局部盆地。 |
| `fdRefErrHz` | `603.221` | 接近 `+1` tooth，但 residual 不是 0。 |
| `fdRateErrHzPerSec` | `-15.246` | rate 基本健康，主现象不再是 rate 崩塌。 |
| `toothIdx` | `1` | final 落在相邻 tooth。 |
| `toothResidualHz` | `-146.779` | `fdRef` 偏差不是严格整数 tooth，还与 DoA / local-state 共同补偿。 |

### 4.3 fdRef line scan

| center | min fd offset | min alias index | center `deltaObj` | `deltaObj(-1500)` | `deltaObj(-750)` | `deltaObj(0)` | `deltaObj(+750)` | `deltaObj(+1500)` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| truth | `0 Hz` | `0` | `0` | `1754.65` | `470.59` | `0` | `344.69` | `1503.33` |
| final estimate | `-750 Hz` | `-1` | `171.16` | `241.33` | `0` | `171.16` | `754.09` | `1746.69` |

- truth center 处正确 tooth 仍是 line scan 全局最优。
- 相邻 `±750 Hz` tooth 的 objective gap 只有几百量级，而普通非 alias 区域的 `deltaObj` 中位数约为 `1.829e6`。
- final-estimate center 的一维最优 tooth 位于 `-750 Hz`，当前 final center 只比该 tooth 差 `171.16`，说明 final 仍停在相邻 comb tooth 附近。

### 4.4 DoA-fdRef surface

| center | min DoA offset | min fd offset | min alias index | center `deltaObj` | min objective |
|---|---:|---:|---:|---:|---:|
| truth | `0 deg` | `0 Hz` | `0` | `0` | `-1.9648e6` |
| final estimate | `-0.0024 deg` | `-600 Hz` | `-1` | `8.7031e5` | `-1.9003e6` |

- truth surface 的最优点仍在 `(0 deg, 0 Hz)`，说明正确 tooth / truth center 是真实 objective peak。
- final-estimate surface 的最优点不是纯 `-750 Hz`，而是 `fdOffset=-600 Hz` 且 `doaOffset=-0.0024 deg` 的联合偏移。
- 该联合偏移把 `fdRef` 从约 `-28337.09 Hz` 拉到约 `-28937.09 Hz`，接近 truth `-28940.31 Hz`。

### 4.5 truth-centered coupling

| center | `rhoDoaFd` | ridge slope | quadratic ridge slope | quadratic fit RMSE | min DoA offset | min fd offset |
|---|---:|---:|---:|---:|---:|---:|
| truth coupling | `0.4934` | `-48.92 Hz/deg` | `-1.2513e4 Hz/deg` | `2.2415e5` | `0 deg` | `0 Hz` |

- `rhoDoaFd≈0.49` 不接近 0，说明 truth 附近 DoA 与 `fdRef` 两个方向存在可见交叉曲率。
- coupling surface 的最优点仍在 truth，说明该诊断不是说 truth 被移走，而是说明局部等高线不是严格轴对齐。
- ridge slope 与二次拟合 slope 差异很大，且二次拟合 RMSE 较高；因此 slope 的具体数值不稳定，不应过度解释。

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `fdRef` line scan 在 truth center 和 final center 都显示 `750 Hz` 周期 comb；这支持 wrong-tooth / near-tooth 是 objective 层结构。
- truth-centered surface 最优点仍在 truth，说明模型/生成链不是把 truth 本身推走。
- final-centered surface 中更优点需要 DoA 与 `fdRef` 联合移动；这解释了为什么单纯按齿硬切或只缩 `fdRef` 盒不一定足够。
- `rhoDoaFd≈0.49` 说明 profile objective 中 DoA 与 `fdRef` 方向存在可见耦合。

### 5.2 仍未解决或反向的现象

- 该 replay 是单 seed 机制可视化，不能说明这种 surface 形态在大 MC 中的覆盖率。
- ridge slope 数值不稳定；当前窗口下 surface 不是理想二次面。
- final-centered 更优点仍来自离线扫描，不代表 runtime flow 已经具备找到该点的能力。
- line / surface 结果只能解释 local objective 结构，不能自动给出 safe gate 或 adoption rule。

### 5.3 代表性 case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `253` | near-tooth / comb hard sample | flow final `toothIdx=1`，`fdRefErrHz=603.221`，surface 更优点在 `doaOffset=-0.0024 deg`、`fdOffset=-600 Hz` | 说明 bad basin 可能需要 DoA-fdRef 联合解释，但仍只是机制诊断。 |

## 6. 机制解释

### 6.1 当前解释

本 replay 把当前 bad sample 拆成两层来看。第一层是 `fdRef` 方向的 `1/T_f` comb：在等间隔多帧连续相位模型下，相邻 tooth 的 objective gap 可以远小于普通非 alias 区域，因此 flow 有机会落到相邻 tooth 或 near-tooth 盆地。第二层是 DoA 与 `fdRef` 的联合补偿：final center 附近的二维 surface 显示，单独移动 `fdRef` 不是唯一更优方向，DoA 同时移动后能进入更接近 truth 的 basin。

因此，该 replay 更适合解释为什么 wrong-tooth / near-tooth 与 same-tooth DoA hard case 会彼此缠绕，而不是直接导出某个固定 slope 策略。若要把 DoA-Doppler coupling 写入论文机制解释，应优先引用 FIM / EFIM 交叉块或更小窗口二次拟合，而不是把当前 ridge slope 数值作为定量结论。

### 6.2 这个结果支持什么

- 支持 `1/T_f` comb 是 objective 层真实结构。
- 支持 DoA 与 `fdRef` 在 profile objective 中存在可见局部耦合。
- 支持把 `replayMfCombToothSurface` 保留为 mechanism visualization / diagnostic 图候选。
- 支持后续 flow 需要关注 same-tooth / near-tooth basin-entry，而不能只看 `fdRef` tooth correction。

### 6.3 这个结果不证明什么

- 不证明 estimator 默认路径已修复。
- 不证明可以把 ridge slope、quadratic ridge slope 或 nonzero `fdRef` offset 下沉为默认 refine 策略。
- 不证明 joint fd bank 应进入默认 candidate family；后续 envelope / gated replay 已显示默认收益主要来自 `wide / single-MF` DoA basin-entry。
- 不证明可以使用 truth center、truth tooth 或 truth DoA 做 runtime selector。
- 不证明 FIM / EFIM 交叉项的严格数值；该结论仍需专门计算。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉 ridge slope、truth override 或 alias-aware objective 字段。 |
| flow 默认路径 | 不直接改；仅作为后续 basin-entry / same-tooth rescue 的机制背景。 |
| regression | 不写；单 seed surface 不是自动 pass/fail 契约。 |
| replay / scan 下一步 | 若需要更强证据，可用更小 coupling 窗口复跑，或转向 CRB / EFIM 交叉块计算。 |
| 论文图 / 论文口径 | 可作为 appendix / diagnostic-only 机制图候选；正文主图仍优先 CRB consistency、resolved RMSE、CP/IP 与 known/unknown-rate 信息损失。 |
| 排障记录 | 若只整理格式，不必新增主记录；若未来作为论文机制图使用，可在机制归并版补一条。 |

## 8. 限制与禁止解释

- 不要把本 replay 的离线 surface minimum 当作 runtime 可用 candidate。
- 不要把 truth-centered coupling 的正结果迁移到 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 ridge slope 数值公式化为默认路径参数。
- 不要用单 seed surface 替代 Monte Carlo 指标、CRB / EFIM 结论或 flow-like 验证。
- 不要恢复 truth override、alias-aware objective injection、blanket DoA release 或 blanket joint fd bank。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfCombToothSurface_20260425-204458.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfCombToothSurface.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本文档只保留 `20260425-204458` 作为当前代表性机制结果。
- 旧 comb / fdRef scan 结果若只重复说明 `1/T_f` comb 存在，可不再单独保留；若未来有更小窗口 coupling 或 EFIM 交叉块结果，再作为新 representative 更新本文件。
