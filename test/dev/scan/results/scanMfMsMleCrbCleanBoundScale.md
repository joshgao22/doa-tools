# scanMfMsMleCrbCleanBoundScale 结果记录

> 本文档只归档当前代表性 snapshot：`scanMfMsMleCrbCleanBoundScale_20260515-223012.mat`。旧占位描述和早期粗网格结论已删除；后续若新增 snapshot，应先判断是否取代本次结果，再更新 Snapshot index。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic / controlled-local range surface / representative` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfMsMleCrbCleanBoundScale_20260515-223012.mat` |
| 当前一句话结论 | 在 6 星 CleanTrim truth-centered local-bound 口径下，MS-MF 主线推荐使用 `truthLocalDoaCrbScale = 2`、`truthLocalFdRefCrbScale = 2.5`：0 dB / 30 repeats 时 `MS-MF-CP-K/U` 的 DoA RMSE 分别约为 `1.0016× / 1.0003× CRB`，`fdRef` RMSE 约为 `1.0784× / 1.0979× CRB`，且 keep rate 为 1、无 boundary reject。 |
| 论文图定位 | 不作为正文主性能图；作为 full-SNR CleanBound / MLE-vs-CRB scan 的 bound 选择依据，可在附录或方法说明中作为 controlled-local sensitivity 证据。 |
| 决策影响 | 后续全 SNR 主扫描建议固定 `(DoA scale, fdRef scale) = (2, 2.5)`；可用 `(2.5, 2.25)` 做 sensitivity。 |
| 下一步动作 | 用 `(2, 2.5)` 跑全 SNR；若需要证明不依赖单一 bound，再用小 repeat 跑 `(2.5, 2.25)` sensitivity。 |
| 禁止误用 | 不能把本 scan 解释为 full-flow acquisition 成功；`1.75×CRB` 附近 RMSE/CRB 低于 1 不应写成 estimator 超过 CRB，而应解释为 truth-centered hard-box / oracle-local 截断效应。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfMsMleCrbCleanBoundScale.m`
- 结果文档：`test/dev/scan/results/scanMfMsMleCrbCleanBoundScale.md`
- scan 类型：controlled-local / oracle-bound range-surface scan。
- 主要问题：`truthLocalDoaCrbScale` 与 `truthLocalFdRefCrbScale` 作为 hard-bound scale 改变时，clean local MLE 的 DoA / `fdRef` RMSE/CRB、joint-trim keep rate、boundary hit 和 outlier reason 如何变化；并据此选择后续全 SNR scan 的固定 bound。
- 扫描对象：DoA hard box 的 CRB 倍数 × `fdRef` hard box 的 CRB 倍数。
- 不覆盖范围：不做 full-flow tooth acquisition，不做 subset selection，不做 candidate adoption，不验证默认 runtime selector，不证明 truth-free basin-entry 能力。
- truth 使用口径：truth 只用于构造 local hard bounds、计算 CRB 和离线评价；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：`Appendix / parameter-selection evidence only`；正文性能图应使用后续固定 bound 的 full-SNR scan。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `doaCrbScale` | DoA latitude / longitude hard box half-width 使用对应方法 CRB 标准差乘以该 scale。 | Oracle | 控制 DoA local-search box 宽度，用于检查 bound 截断和 RMSE/CRB knee。 | 不等价于真实 coarse acquisition 误差，也不能写成 runtime selector。 |
| `fdRefCrbScale` | `fdRef` hard box half-width 使用对应方法 `fdRef` CRB 标准差乘以该 scale，并由同 tooth cap 限制。 | Oracle | 控制参考星 Doppler local-search box 宽度，用于检查 `fdRef` boundary 和 fd RMSE/CRB knee。 | 不证明全局 Doppler tooth acquisition 已解决。 |
| `scaleReadableRmseCrbTable` | 主读表，比较 raw 与 joint-trim 后的 DoA / `fdRef` RMSE/CRB、keep rate 和 reject reason。 | Eval only | 直接用于选择后续固定 bound。 | 不能只看最小 RMSE；还要看是否低于 CRB、是否 boundary-limited。 |
| `scaleSurfaceSummaryTable` | 每个 method / SNR 按当前 score 选出的 best scale 点。 | Eval only | 可快速找到 surface minimum。 | best 点若 RMSE/CRB < 1，不应自动作为 paper-facing 主设置。 |
| `scaleSearchRangeAuditTable` | 实际 hard bound 相对 CRB 的 coverage、boundary hit 和 truth box violation 诊断。 | Oracle / Eval | 判断范围是否低于 CRB、truth 是否落在 hard box、是否贴边。 | `fdRefBelowCrb=false` 不等价于没有 boundary 截断，仍需看 `fdRefBoundaryHitRate`。 |
| `joint-trim` | health 且 angle / `fdRef` 误差均不超过 `5×CRB` 的 trimmed sample。 | Eval only | 用于 local MLE-vs-CRB 口径。 | 不能当作 full-sample 工程成功率。 |

常用口径：

- `raw`：全部样本统计。
- `health`：去除 solver / boundary 等 health 失败后的统计。
- `joint-trim`：health 且 angle / `fdRef` 均在 `5×CRB` 内的统计。
- `trimKeepRate`：进入 joint-trim 的比例。本 scan 中 MS-MF 主线在所有细化格点均为 1。
- `topTrimRejectReason = "boundary"`：主要表示被 health / boundary 判据排除，不表示理论 CRB 错误。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfMsMleCrbCleanBoundScale_20260515-223012.mat` | 2026-05-15 | `representative` | `SNR=0 dB`，`numRepeat=30`，`seed=253:282`，`P=10`，6 星 `[5259 1243 348 5652 14 4437]`，`refSat=5259`，DoA scale `[1.75 2 2.5 3]`，`fdRef` scale `[2.25 2.5 2.75 3 3.5]`，methods = `SS/MS-MF-CP-K/U`。 | MS-MF 主推荐 bound 为 `(2, 2.5)`；`(2.5, 2.25)` 可作为 sensitivity；`1.75` 的 DoA RMSE/CRB 偏低，不宜作为主设置。 | 覆盖早期“尚无代表性 snapshot”的占位文档；本文只保留本次结果。 |

维护规则：当前代表性结果只保留本次 snapshot。若未来新增 full-SNR scan，应另写对应结果文档或在本文中明确标为 fixed-bound follow-up，不要把 surface-scan 与 full-SNR 性能曲线混写。

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:282`
- `numRepeat = 30`
- `snrDbList = 0`
- `frameCount = 10`
- `user LLA = [55, 36.59, 0]`
- `selectedSatIdxGlobal = [5259 1243 348 5652 14 4437]`
- `refSatIdxGlobal = 5259`
- `methodNameList = ["SS-MF-CP-K", "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"]`
- `truthLocalDoaCrbScaleList = [1.75, 2, 2.5, 3]`
- `truthLocalFdRefCrbScaleList = [2.25, 2.5, 2.75, 3, 3.5]`
- `fdRef range mode = crb-scale / 0.45 tooth cap`
- trim definition：health + angle `<= 5×CRB` + `fdRef <= 5×CRB`
- checkpoint：`600/600` complete，success 后已清理
- snapshot 保存变量：`scanData`
- 总 wall time：未在当前 summary 中直接打印；runtime 表保存了各 stage 累计耗时。

### 4.2 存档数据检查

- 顶层变量：`scanData`
- 关键字段：`scaleRawAggregateTable`、`scaleHealthAggregateTable`、`scaleJointTrimAggregateTable`、`scaleReadableRmseCrbTable`、`scaleSurfaceSummaryTable`、`scaleSearchRangeAuditTable`、`scaleOutlierTable`、`runtimeAggregateTable`、`topSlowRuntimeTable`、`checkpointSummaryTable`
- 未保存大体量数据：未见 `rxSigCell`、完整 `sceneSeq`、fixture cache、全量 objective map 或图片进入 snapshot。
- warning / fail：summary 中未显示 run-level failure；outlier table 共 24 行，主要来自 SS-MF 的 boundary case。
- runtime 观察：MS-MF-CP-U 最重，`method-total` 累计约 `6340.6 s / 600 cases`，median `10.426 s`；MS-MF-CP-K 累计约 `4272.2 s`，median `6.9543 s`。

## 5. 主要统计与曲线结果

### 5.1 主推荐点

| bound | case | trim samples | DoA RMSE / CRB | fdRef RMSE / CRB | trim status | top reject | 解释 |
|---|---|---:|---:|---:|---|---|---|
| `(2, 2.5)` | `MS-MF-CP-K` | `30/30` | `1.0016` | `1.0784` | `near-crb` | `""` | 主推荐点；DoA 几乎正好贴 CRB 且略大于 CRB。 |
| `(2, 2.5)` | `MS-MF-CP-U` | `30/30` | `1.0003` | `1.0979` | `near-crb` | `""` | 主推荐点；unknown-rate 主参数仍贴近 CRB，`fdRef` 略高于 CRB。 |
| `(2.5, 2.25)` | `MS-MF-CP-K` | `30/30` | `1.0529` | `1.0653` | `near-crb` | `""` | sensitivity 候选；DoA 稍宽，两个维度均大于 CRB。 |
| `(2.5, 2.25)` | `MS-MF-CP-U` | `30/30` | `1.0344` | `1.0480` | `near-crb` | `""` | sensitivity 候选；`fdRef` 更贴 CRB，但 DoA 比主推荐点稍远。 |
| `(2, 3)` | `MS-MF-CP-K` | `30/30` | `1.0125` | `1.0834` | `near-crb` | `""` | 可用但没有 `(2,2.5)` 更贴 DoA CRB。 |
| `(2, 3)` | `MS-MF-CP-U` | `30/30` | `0.99582` | `1.1141` | `kept` | `""` | 可用但 CP-U 的 DoA 略低于 CRB，`fdRef` 更高。 |

### 5.2 surface best 点与人工选择差异

| case | surface best bound | best DoA RMSE / CRB | best fdRef RMSE / CRB | keep rate | 人工判断 |
|---|---|---:|---:|---:|---|
| `MS-MF-CP-K` | `(1.75, 3.5)` | `0.96902` | `1.0792` | `1` | DoA 低于 CRB，偏 oracle-tight，不作为主设置。 |
| `MS-MF-CP-U` | `(1.75, 3)` | `0.95095` | `1.0769` | `1` | DoA 低于 CRB，偏 oracle-tight，不作为主设置。 |
| `SS-MF-CP-K` | `(1.75, 2.75)` | `0.87282` | `0.88893` | `0.96667` | SS 受 bounded-local 影响更明显，不用于反推主 bound。 |
| `SS-MF-CP-U` | `(1.75, 2.25)` | `0.88876` | `0.97685` | `0.96667` | 存在 boundary reject，不作为主设置。 |

### 5.3 DoA scale 对 MS 主线的影响

固定 `fdRefCrbScale = 2.5` 时：

| `doaCrbScale` | `MS-MF-CP-K` DoA RMSE / CRB | `MS-MF-CP-U` DoA RMSE / CRB | 解释 |
|---:|---:|---:|---|
| `1.75` | `0.9776` | `0.97825` | 略低于 CRB，仍有 hard-box 截断嫌疑。 |
| `2` | `1.0016` | `1.0003` | 最贴近 CRB 且略大于 CRB。 |
| `2.5` | `1.0388` | `1.0541` | 稍宽，作为 sensitivity 合理。 |
| `3` | `1.0650` | `1.0571` | 更宽但 DoA RMSE 明显抬高。 |

### 5.4 fdRef scale 对 MS 主线的影响

固定 `doaCrbScale = 2` 时：

| `fdRefCrbScale` | `MS-MF-CP-K` fdRef RMSE / CRB | `MS-MF-CP-U` fdRef RMSE / CRB | MS keep rate | 解释 |
|---:|---:|---:|---:|---|
| `2.25` | `1.0530` | `1.0906` | `1` | 已无 MS boundary；K 的 fdRef 更贴 CRB，但 U 不比 `2.5` 明显好。 |
| `2.5` | `1.0784` | `1.0979` | `1` | 主推荐；DoA 与 fdRef 综合最均衡。 |
| `2.75` | `1.0833` | `1.1007` | `1` | 与 `2.5` 接近，但没有收益。 |
| `3` | `1.0834` | `1.1141` | `1` | 稳定但 `fdRef` 略高。 |
| `3.5` | `1.0881` | `1.1073` | `1` | 稳定但更宽，没有主设置价值。 |

### 5.5 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| DoA scale × fdRef scale heatmap | `doaCrbScale` / `fdRefCrbScale` | DoA RMSE / CRB、`fdRef` RMSE / CRB | `MS-MF-CP-K/U` 为主，SS 作辅助 | Appendix / internal diagnostic | 不要把 surface best 直接当主设置；低于 1 的格点要标注 oracle-bound。 |
| fixed-bound sensitivity table | bound pair | RMSE/CRB + keep rate | `(2,2.5)` vs `(2.5,2.25)` | Appendix candidate | 用于说明 full-SNR 主图的 bound 选择不是任意调参。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- **MS-MF 主线在细化网格内已经无 boundary 污染。** `scaleSearchRangeAuditTable` 中 `MS-MF-CP-K/U` 所有 `doaCrbScale × fdRefCrbScale` 组合均为 `fdRefBoundaryHitRate = 0`、`doaTruthBoxViolationRate = 0`、`doaBelowCrb=false`、`fdRefBelowCrb=false`。
- **`doaCrbScale=2` 是最干净的 DoA knee。** 固定 `fdRefCrbScale=2.5` 时，MS-MF-CP-K/U 的 DoA RMSE/CRB 从 `1.75` 的约 `0.978` 上升到 `2` 的约 `1.00`，再到 `2.5/3` 的 `1.03–1.06`；`2` 最符合“贴近 CRB 且略大于 CRB”的目标。
- **`fdRefCrbScale=2.5` 已足够稳定。** 固定 `doaCrbScale=2` 时，`fdRef` RMSE/CRB 在 `2.25–3.5` 之间都稳定在约 `1.05–1.11`，且 MS keep rate 均为 1；继续放宽没有收益。
- **SS-MF 不适合反向决定主 bound。** SS-MF 的 DoA RMSE/CRB 长期低于 1，且部分 `fdRef` scale 下存在 1 个 seed 的 boundary reject；这更像 bounded-local / SS-specific 截断现象，不应牵引 MS 主线 bound 选择。

### 6.2 反向、污染或未解决现象

- `scaleSurfaceSummaryTable` 的 best point 对 MS 落在 `doaCrbScale=1.75`，但这主要因为更窄 DoA hard box 压低了 DoA RMSE/CRB；它不是 paper-facing 主设置。
- `SS-MF-CP-K` 在 `fdRefCrbScale < 3.5` 时常有 `trimKeepRate=0.96667` 和 `topRejectReason="boundary"`；`SS-MF-CP-U` 在 `fdRefCrbScale < 2.75` 时也有类似单 seed boundary。该现象不影响 MS bound 选择，但说明 SS 在这个 surface 中仍有局部截断 / boundary tail。
- 本 scan 只在 `0 dB` 下做 30 repeats，不能推出低 SNR threshold 区、高 SNR branch tail 或 full-SNR curve 的最终形状。

### 6.3 代表性异常格点 / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `SS-MF-CP-K`, seed `255`, 多个 fdRef scale | boundary | `fdRefErrOverCrb` 基本贴在对应 scale 右侧，例如 `2.25 / 2.5 / 2.75 / 3.0` 附近；`fdRefBoundaryHit=true`。 | 说明 SS 的局部 `fdRef` box 仍可能截断，不能用 SS surface best 决定主设置。 |
| `SS-MF-CP-U`, seed `255`, `fdRefScale=2.25/2.5` | boundary | `fdRefBoundaryHit=true`，`fdRateAbsErrHzPerSec` 约 `19–22`。 | unknown-rate SS 有类似单 seed boundary；主结论仍以 MS-MF 为准。 |
| `MS-MF-CP-K/U`, 全部细化格点 | clean | `trimKeepRate=1`，`topRejectReason=""`。 | 支持 `(2,2.5)` 可作为 MS full-SNR 主 bound。 |

## 7. 机制解释

### 7.1 当前解释

这个 scan 是 truth-centered hard-box surface，不是 full-flow 搜索。小的 `doaCrbScale` 会把 DoA 求解限制在真值附近很窄的 CRB box 内，因此可以人为压低 DoA RMSE/CRB，甚至低于 1。随着 DoA box 从 `1.75×CRB` 放宽到 `2×CRB`，MS-MF 的 DoA RMSE/CRB 从略低于 1 过渡到约 1；继续放宽到 `2.5–3×CRB` 后，finite-sample local basin / noise-induced variation 开始反映出来，DoA RMSE/CRB 略高于 1。

`fdRefCrbScale` 的作用在本次 MS surface 中相对温和。早期粗网格中 `fdRefScale=1.5/2` 仍可能带来 boundary reject；本次细化从 `2.25` 起已经让 MS-MF 全部格点无 boundary hit，说明 `fdRef` local box 的 knee 已经被覆盖。综合 DoA 与 `fdRef` 两个指标，`(2,2.5)` 是当前最稳妥的折中点。

### 7.2 这个 scan 支持什么

- 支持后续 full-SNR CleanBound 主 scan 使用 `(truthLocalDoaCrbScale, truthLocalFdRefCrbScale) = (2, 2.5)`。
- 支持将 `(2.5, 2.25)` 作为 sensitivity bound：该点两个主参数都略大于 CRB，但 DoA 稍远。
- 支持把 `1.75×CRB` 解释为偏 oracle-tight，不作为 paper-facing 主设置。
- 支持 MS-MF 在 0 dB、controlled local-bound 条件下已经可以达到约 `1.00×` DoA CRB 和 `1.08–1.10×` `fdRef` CRB。

### 7.3 这个 scan 不证明什么

- 不证明 estimator 默认 full-flow 已经解决 tooth acquisition、same-tooth bad basin 或 candidate adoption。
- 不证明全 SNR 下每个 SNR 都能保持 `1.0–1.1×CRB`；低 SNR 仍可能出现 threshold / outlier，高 SNR 仍需检查 fd branch tail。
- 不证明可以把 truth-centered bounds 迁移到 runtime selector。
- 不证明 SS-MF 的 DoA RMSE/CRB 低于 1 是 CRB 错误；更合理解释是 bounded-local / oracle-box 截断。
- 不应写成 regression 契约；这是参数选择和机制诊断 scan。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。该结果只选择 dev/scan 的 controlled-local bound，不改变 estimator objective、residual、初始化、reference-sat 语义或默认路径。 |
| flow 默认路径 | 不改。不进入 full-flow selector / rescue / adoption。 |
| replay 下一步 | 不需要新增 replay；可用 `replayMfMsMleCrbCleanTrim` 在 `(2,2.5)` 下做单点 sanity check。 |
| scan 下一步 | full-SNR 主 scan 固定 `(2,2.5)`；可选 sensitivity 固定 `(2.5,2.25)`。 |
| regression | 不写。bound surface 不是稳定自动契约。 |
| 论文图 | 主文性能图不直接用该 surface；可作为 appendix / parameter-selection 说明。 |
| 排障记录 | 可在主记录或机制归并版中摘一句：CleanBoundScale 细化网格支持全 SNR 主 bound 取 `(2,2.5)`。 |

## 9. 限制与禁止解释

- 不要把 `doaCrbScale=1.75` 的较小 RMSE 写成 estimator 优于 CRB。
- 不要用 `scaleSurfaceSummaryTable` 的 best point 直接替代人工 bound 选择；best point可能是 oracle-box 压缩后的 minimum。
- 不要把本 scan 的 `trimKeepRate=1` 解释为 full-flow resolved rate = 1。
- 不要把 truth-centered local hard bounds 用于 runtime selector、gate、candidate adoption 或 final winner。
- 不要因为 SS-MF 的 boundary seed 去扩大 MS 主 bound；SS 在本 scan 中只作辅助 baseline。
- 不要把 0 dB / 30 repeats 的 surface 直接推广到所有 SNR；必须通过固定 bound 的 full-SNR scan 验证。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfMsMleCrbCleanBoundScale_20260515-223012.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfMsMleCrbCleanBoundScale.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table、surface summary、range audit、outlier table 和 runtime table。

## 11. 历史备注

- 早期占位文档只说明该 scan 尚无代表性 snapshot；已由本次 `20260515-223012` 结果取代。
- 粗网格曾提示 `(2,3)` 可作为保守候选；本次细化后，`(2,2.5)` 更符合“MLE 尽可能接近 CRB 且略大于 CRB”的目标，因此升级为主推荐 bound。
