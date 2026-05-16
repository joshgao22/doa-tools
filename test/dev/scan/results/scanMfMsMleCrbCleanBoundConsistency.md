# scanMfMsMleCrbCleanBoundConsistency 结果记录

> 本文档只归档当前代表性 snapshot：`scanMfMsMleCrbCleanBoundConsistency_20260516-135405.mat`。旧占位描述已删除；后续若新增 snapshot，应先判断是否取代本次 full-SNR 结果，再更新 Snapshot index。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing / controlled-local / representative` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfMsMleCrbCleanBoundConsistency_20260516-135405.mat` |
| 当前一句话结论 | 在 6 星、10 帧、truth-centered CRB-scaled local bound `(DoA, fdRef) = (2, 2.5)` 下，`MS-MF-CP-K/U` 在 joint-trim 口径下达到 CRB-level；MS 绝对 DoA 误差随对应更紧 CRB 显著低于 SS，CP 相比 IP 显著保留跨帧连续相位带来的 DoA / `fdRef` 信息量。 |
| 论文图定位 | 可作为正文或附录的 controlled-local MLE-vs-CRB / CP-IP / K-U 对照结果；图注必须标明 truth-centered local bound 和 joint-trim 口径。 |
| 决策影响 | 固定本次 `(2, 2.5)` bound 作为 CleanBound 全 SNR 主结果；后续论文图可从 `cleanEstimateCrbCurveTable` / `cleanSnrCurveTable` 出图，若需要 sensitivity 再单独跑 `(2.5, 2.25)`。 |
| 下一步动作 | 结果可进入论文图候选；若写入正文，应同时报告 keep rate / outlier 口径，并说明 MSE/CRB 略低于 1 来自 controlled local bound 和 trim。 |
| 禁止误用 | 不能把本结果解释为 full-flow acquisition、tooth selection、candidate adoption 或 same-tooth rescue 已通过；不能把低中 SNR 的 MSE/CRB 略低于 1 写成 estimator 超过无约束 CRB。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfMsMleCrbCleanBoundConsistency.m`
- 结果文档：`test/dev/scan/results/scanMfMsMleCrbCleanBoundConsistency.md`
- scan 类型：paper-facing controlled-local / oracle-bound SNR curve。
- 主要问题：在由 CleanBoundScale 选出的 fixed local bound 下，SS/MS、SF/MF、CP/IP、K/U 的 MLE 是否与对应 CRB 一致，并能否支撑 continuous-phase、multi-satellite 和 known/unknown-rate 的论文结论。
- 扫描对象：`SNR × seed × method`，其中 SNR 为 `[-15, -10, -5, 0, 5, 10] dB`，每个 SNR 500 seeds。
- 不覆盖范围：不做 full-flow tooth acquisition、subset schedule、candidate adoption、basin-entry rescue 或同齿 hard-tail 修复。
- truth 使用口径：truth 只用于构造 local hard bounds、计算误差 / CRB 和离线 trim；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Yes，但只能作为 controlled-local / resolved-regime 证据；不是 full-flow 工程成功率图。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `raw` | 全部 500 seeds 的统计。 | Eval only | 用于观察 full local-bound 下的原始 tail 和 boundary。 | 不能单独作为 CRB-local 结论。 |
| `health` | 去除 solver / boundary / health 失败后的统计。 | Eval only | 用于区分数值健康样本和边界 / fdRate 异常。 | 不能等同于真实捕获成功率。 |
| `joint-trim` | health 且 angle / `fdRef` 误差均不超过 `5×CRB` 的统计。 | Eval only | 主 MLE-vs-CRB 口径；本文主要用它判断 CRB-level。 | 不能解释为 full-flow resolved rate；trim 使用 truth 评价。 |
| `angleMseOverSphericalCrb` | angle MSE 相对 spherical angular CRB 方差的比值。 | Eval only | 主 CRB-normalized angle 指标；接近 1 表示 local MLE 与对应 CRB 同阶。 | 不要与 trace CRB 混用；也不要把略低于 1 写成超 CRB。 |
| `fdRefMseOverCrb` | `fdRef` MSE 相对对应 `fdRef` CRB 方差的比值。 | Eval only | 判断 reference Doppler 是否达到 CRB-level。 | 不证明 tooth acquisition 成功。 |
| `keepRate` | 进入 health / joint-trim 的样本比例。 | Eval only | 必须和 RMSE/CRB 同时报告，避免只看筛选后曲线。 | 不等同于 runtime hit rate。 |
| `CP-K / CP-U` | continuous-phase，known / unknown reference Doppler-rate。 | No | 论文主模型和 nuisance-rate 对照。 | 不要把 K/U 的小差异解释成没有 nuisance；还要看 EFIM loss 和 keep rate。 |
| `IP-K / IP-U` | independent-phase，切断跨帧公共相位。 | No | CP 的受控模型层级 baseline。 | IP 自身 CRB-level 不表示 IP 信息量与 CP 相同。 |
| `crb-scale` bound | DoA half-width = `2×CRB`，`fdRef` half-width = `2.5×CRB`，并受同 tooth cap 约束。 | Oracle | controlled local-bound 口径，用于验证 resolved/local MLE。 | 不能迁移为实际 coarse acquisition 或 runtime selector。 |

常用口径：

- `full-sample / raw`：所有 local-bound 任务的统计，用于暴露 tail。
- `joint-trim`：本文主读口径，用于 MLE-vs-CRB 曲线和 paper-facing controlled-local 结论。
- `outlier table`：只解释失败机制，不用于调整 estimator 默认路径。
- `truth-centered local bound`：本文主图必须标注的 oracle 条件。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfMsMleCrbCleanBoundConsistency_20260516-135405.mat` | 2026-05-16 | `representative` | `P=10`，`SNR=-15:5:10 dB`，`numRepeat=500`，`seed=253:752`，6 星 `[5259 1243 348 5652 14 4437]`，`refSat=5259`，`DoA scale=2`，`fdRef mode=crb-scale`，`fdRef scale=2.5`，methods = SS/MS、SF/MF、CP/IP、K/U。 | CP 主线达到 CRB-level；MS 显著降低绝对 DoA 误差；CP 相比 IP 保留跨帧连续相位信息；K/U 主参数在 CRB-local 口径下差异较小。 | 覆盖早期“pending first run”的占位文档；本文只保留本次最终结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:752`
- `numRepeat = 500`
- `snrDbList = [-15, -10, -5, 0, 5, 10]`
- `taskCount = 3000`
- `frameCount = 10`
- `user LLA = [55, 36.59, 0]`
- `selectedSatIdxGlobal = [5259, 1243, 348, 5652, 14, 4437]`
- `refSatIdxGlobal = 5259`
- `methodNameList = [SS-SF-DoA, MS-SF-DoA, SS-MF-Static, MS-MF-Static, SS-MF-CP-K, SS-MF-CP-U, MS-MF-CP-K, MS-MF-CP-U, SS-MF-IP-K, SS-MF-IP-U, MS-MF-IP-K, MS-MF-IP-U]`
- `defaultDoaCrbScale = 2`
- `defaultFdRefRangeMode = "crb-scale"`
- `defaultFdRefCrbScale = 2.5`
- `fdRef CRB-scale tooth cap = 0.45 tooth`
- `fdRateHalfWidthHzPerSec = 1000`
- trim definition：health + angle `<= 5×CRB` + `fdRef <= 5×CRB`
- checkpoint：`3000/3000` complete；`cleanedOnSuccess = false`，本次 summary 从已完成 checkpoint / snapshot 恢复。
- snapshot 保存变量：`scanData`
- summary 恢复耗时：snapshot 中 `elapsedSec = 157.35 s`；完整求解 wall time 不由该字段代表。

### 4.2 存档数据检查

- 顶层变量：`data.scanData`，并带 `meta` / `inventory`。
- 关键字段：`config`、`contextSummary`、`caseTable`、`cleanAggregateTable`、`cleanHealthAggregateTable`、`cleanTrimAggregateTable`、`cleanSsMsCompareTable`、`cleanKnownUnknownCompareTable`、`cleanSnrCurveTable`、`cleanEstimateCrbCurveTable`、`cleanSearchRangeAuditTable`、`cleanOutlierTable`、`runtimeAggregateTable`、`topSlowRuntimeTable`、`checkpointSummaryTable`、`plotData`。
- 未保存大体量数据：未见 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map 或图片进入 snapshot。
- run-level failure：summary 中显示 task 全部完成；outlier 表主要用于机制解释，非 run failure。
- 主要 runtime 观察：最重 method 为 `MS-MF-IP-U`，`method-total` 累计约 `42173 s / 3000 cases`，median `12.276 s`；其次为 `MS-MF-CP-U`，累计约 `32495 s`，median `10.407 s`。

## 5. 主要统计与曲线结果

### 5.1 MS-MF-CP 主线：joint-trim MLE-vs-CRB

| SNR (dB) | case | angle MSE / CRB | fdRef MSE / CRB | angle RMSE (deg) | fdRef RMSE (Hz) | keep rate | 结论 |
|---:|---|---:|---:|---:|---:|---:|---|
| -15 | `MS-MF-CP-K` | `0.92345` | `0.87409` | `2.3503e-05` | `0.11158` | `0.970` | CRB-level，低中 SNR 略低于 1，受 local bound / trim 影响。 |
| -15 | `MS-MF-CP-U` | `0.93198` | `0.92731` | `2.3611e-05` | `0.12342` | `0.984` | CRB-level；U 的 keep rate 高于 K。 |
| -10 | `MS-MF-CP-K` | `0.93393` | `0.88411` | `1.3291e-05` | `0.063104` | `0.958` | CRB-level。 |
| -10 | `MS-MF-CP-U` | `0.92107` | `0.93560` | `1.3199e-05` | `0.069713` | `0.986` | CRB-level。 |
| -5 | `MS-MF-CP-K` | `0.90856` | `0.91892` | `7.3720e-06` | `0.036178` | `0.960` | CRB-level。 |
| -5 | `MS-MF-CP-U` | `0.91691` | `0.92392` | `7.4058e-06` | `0.038957` | `0.986` | CRB-level。 |
| 0 | `MS-MF-CP-K` | `0.90141` | `0.89199` | `4.1292e-06` | `0.020044` | `0.956` | CRB-level，略低于 1。 |
| 0 | `MS-MF-CP-U` | `0.89782` | `0.91044` | `4.1210e-06` | `0.021747` | `0.986` | CRB-level。 |
| 5 | `MS-MF-CP-K` | `0.93857` | `0.90857` | `2.3694e-06` | `0.011376` | `0.946` | CRB-level。 |
| 5 | `MS-MF-CP-U` | `0.91470` | `0.93458` | `2.3391e-06` | `0.012390` | `0.992` | CRB-level。 |
| 10 | `MS-MF-CP-K` | `1.16680` | `0.90469` | `1.4856e-06` | `0.0063834` | `0.862` | 高 SNR 下 angle 略高于 CRB，且有 boundary seed。 |
| 10 | `MS-MF-CP-U` | `1.10170` | `0.95138` | `1.4436e-06` | `0.0070298` | `0.992` | 高 SNR 下 angle 略高于 CRB，仍为 CRB-level。 |

### 5.2 SS-MF vs MS-MF：绝对误差和多星信息量

| SNR (dB) | pair | SS angle RMSE (deg) | MS angle RMSE (deg) | SS angle MSE / CRB | MS angle MSE / CRB | 解释 |
|---:|---|---:|---:|---:|---:|---|
| -15 | `CP-K` | `0.017569` | `2.3503e-05` | `0.96238` | `0.92345` | MS 的绝对 DoA 误差远低于 SS，同时仍贴自己的更紧 CRB。 |
| 0 | `CP-K` | `0.0031161` | `4.1292e-06` | `0.95737` | `0.90141` | MS 误差随 MS CRB 同步下降，不是用 SS CRB 分母制造增益。 |
| 10 | `CP-K` | `0.00096439` | `1.4856e-06` | `0.91698` | `1.16680` | 高 SNR 下 MS angle 稍高于 CRB，但绝对误差仍显著小于 SS。 |
| -15 | `CP-U` | `0.017712` | `2.3611e-05` | `0.97810` | `0.93198` | unknown-rate 主参数仍为 CRB-level。 |
| 0 | `CP-U` | `0.0031211` | `4.1210e-06` | `0.96043` | `0.89782` | MS 形成数量级上的 DoA 绝对精度提升。 |
| 10 | `CP-U` | `0.00097538` | `1.4436e-06` | `0.93801` | `1.10170` | 高 SNR 下 MS 仍为 CRB-level，且 keep rate 高。 |

### 5.3 CP vs IP：连续相位信息量

| SNR (dB) | MS case | angle RMSE / CRB | angle RMSE (deg) | angle CRB (deg) | fdRef RMSE / CRB | fdRef RMSE (Hz) | fdRef CRB (Hz) | keep rate |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | `MS-MF-CP-K` | `0.94942` | `4.1292e-06` | `4.3492e-06` | `0.94445` | `0.020044` | `0.021223` | `0.956` |
| 0 | `MS-MF-IP-K` | `0.94851` | `0.0012482` | `0.001316` | `0.96974` | `67.476` | `69.581` | `0.984` |
| 0 | `MS-MF-CP-U` | `0.94753` | `4.1210e-06` | `4.3492e-06` | `0.95417` | `0.021747` | `0.022791` | `0.986` |
| 0 | `MS-MF-IP-U` | `0.94029` | `0.0012374` | `0.001316` | `0.90671` | `64.044` | `70.634` | `0.452` |
| 10 | `MS-MF-CP-K` | `1.0802` | `1.4856e-06` | `1.3753e-06` | `0.95115` | `0.0063834` | `0.0067112` | `0.862` |
| 10 | `MS-MF-IP-K` | `0.94812` | `0.00039455` | `0.00041614` | `0.96936` | `21.329` | `22.003` | `0.984` |
| 10 | `MS-MF-CP-U` | `1.0496` | `1.4436e-06` | `1.3753e-06` | `0.97539` | `0.0070298` | `0.0072072` | `0.992` |
| 10 | `MS-MF-IP-U` | `0.95028` | `0.00039545` | `0.00041614` | `0.93477` | `20.879` | `22.336` | `0.578` |

> 读法：IP 也能相对自己的 CRB 达到 CRB-level，但它的 angle CRB 和 `fdRef` CRB 远大于 CP；因此 CP 的优势不是“归一化误差一定更小”，而是保留了跨帧连续相位带来的绝对信息量。

### 5.4 K/U 对照：unknown-rate 的主参数损失较小

| SNR (dB) | sat / phase | known angle MSE / CRB | unknown angle MSE / CRB | Δ angle MSE / CRB | known fdRef MSE / CRB | unknown fdRef MSE / CRB | Δ fdRef MSE / CRB | 解释 |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| -15 | `MS / CP` | `0.92345` | `0.93198` | `+0.00853` | `0.87409` | `0.92731` | `+0.05322` | U 略损失但仍 CRB-level。 |
| 0 | `MS / CP` | `0.90141` | `0.89782` | `-0.00359` | `0.89199` | `0.91044` | `+0.01845` | angle 几乎无损失，fdRef 小幅膨胀。 |
| 10 | `MS / CP` | `1.16680` | `1.10170` | `-0.06510` | `0.90469` | `0.95138` | `+0.04669` | 高 SNR U 的 keep rate 更高；不能只看 K/U 单调性。 |
| 0 | `MS / IP` | `0.89967` | `0.88414` | `-0.01553` | `0.94040` | `0.82212` | `-0.11828` | IP-U 的 trimmed 样本有强筛选效应，需结合低 keep rate。 |
| 10 | `MS / IP` | `0.89894` | `0.90302` | `+0.00408` | `0.93966` | `0.87379` | `-0.06587` | IP-U 不能单独作为信息损失主证据。 |

### 5.5 搜索范围与异常机制

| 项目 | 观察 | 解释 |
|---|---|---|
| CP / static 的 DoA range | 所有 method / SNR 的 `minDoaLatHalfOverCrb = minDoaLonHalfOverCrb = 2`，`doaBelowCrb=false`，`doaTruthBoxViolationRate=0`。 | DoA hard box 按计划执行，没有 truth box violation。 |
| CP / static 的 `fdRef` range | `SS/MS-MF-CP-K/U` 和 static 的 `minFdRefHalfOverCrb = 2.5`，`fdRefBelowCrb=false`。 | CP 主线使用了 CleanBoundScale 选出的 `2.5×CRB` fdRef bound。 |
| IP 的低 SNR `fdRef` range | `SS/MS-MF-IP-*` 在 -15 dB 到部分 -10 dB 下可能出现 `fdRefBelowCrb=true`，因为同 tooth cap 限制后实际 half-width 小于 IP 自身更大的 CRB。 | 解释 IP 低 SNR boundary / keepRate 下降；不应混入 CP 主线。 |
| static 行 | `SS/MS-MF-Static` 在 health / joint-trim 下 `keepRate=0`。 | static 在动态 health / joint-trim 口径下不适合作为主判断，只作反例和 anchor。 |
| outlier preview | 主要由 `MS-MF-IP-U` 的 `fdRate+boundary` / `fdRef+boundary` 组成；另有 `MS-MF-CP-U, SNR=10, seed=347` 的 boundary case。 | IP-U 的 unknown-rate / boundary 敏感性强；CP-U 主线仍较稳。 |

### 5.6 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Estimate-vs-CRB angle curve | SNR (dB) | angle RMSE 与 spherical CRB | SS/MS、SF/MF、CP/IP、K/U | Yes | 建议主图突出 MS-MF-CP-K/U、SS-MF-CP-K/U；IP 可作对照或副图。 |
| Estimate-vs-CRB `fdRef` curve | SNR (dB) | `fdRef` RMSE 与 CRB | MF CP/IP K/U | Yes / appendix | CP 与 IP 的绝对 `fdRef` CRB 差异巨大，适合说明 continuous-phase 信息量。 |
| Clean SNR curve summary | SNR (dB) | MSE/CRB、keep rate | SS/MS pair | Yes | 论文中应同时报告 keep rate，避免 trimmed-only 误导。 |
| Outlier preview | SNR / seed | error-over-CRB、boundary reason | MS-MF-IP-U 和少量 CP-U | No | 只作机制说明，不作为主性能图。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- **CP 主线达到 CRB-level。** `MS-MF-CP-K/U` 在 -15 到 10 dB 的 joint-trim angle MSE/CRB 大致落在 `0.90–1.17`，`fdRef` MSE/CRB 大致落在 `0.87–0.95`；低中 SNR 略低于 1，高 SNR angle 略高于 1，整体符合 controlled-local CRB-level 结果。
- **MS 显著降低绝对 DoA 误差。** 0 dB 时 `SS-MF-CP-U` 的 angle RMSE 为 `0.0031211 deg`，而 `MS-MF-CP-U` 为 `4.1210e-06 deg`；二者都贴自己的 CRB，说明 MS 的绝对增益来自更强几何 / 多星信息量，而不是比较口径错配。
- **CP 相比 IP 保留了数量级更强的信息量。** 0 dB 时 `MS-MF-CP-K` angle CRB 为 `4.3492e-06 deg`，`MS-MF-IP-K` 为 `0.001316 deg`；`fdRef` CRB 分别为 `0.021223 Hz` 和 `69.581 Hz`。IP 虽然也能相对自己的 CRB-level，但绝对边界远松于 CP。
- **K/U 主参数差异较小。** 对 MS-MF-CP，K/U 的 angle MSE/CRB 曲线非常接近，U 的 `fdRef` MSE/CRB 多数只小幅变化；这说明当前设置下 unknown `fdRate` 的主要影响不是显著抬高 trimmed DoA CRB，而是进入 health / resolved rate 与 `fdRate` nuisance 稳定性。
- **IP-U 的 keep rate 明显较弱。** `MS-MF-IP-U` 的 joint-trim keep rate 约为 `0.35, 0.524, 0.586, 0.452, 0.456, 0.578`，明显低于 `MS-MF-CP-U` 的 `0.984–0.992` 主区间。

### 6.2 反向、污染或未解决现象

- **MSE/CRB 略低于 1 不能写成 estimator 超 CRB。** 低中 SNR 的 CP 主线多在 `0.90–0.96`，应解释为 truth-centered finite box + joint-trim 截尾 + median CRB 标准化共同作用，而不是无约束 unbiased CRB 被突破。
- **`MS-MF-CP-K` 在 10 dB 的 keep rate 下降到 `0.862`。** outlier preview 中存在 `MS-MF-CP-U, seed=347` boundary case；高 SNR 下 branch / boundary tail 仍需在图注或正文中保留 resolved-rate 说明。
- **static 在 dynamic health 下全空。** `SS/MS-MF-Static` 在 joint-trim 下 keep rate 为 0；raw 中的 static angle 不能作为动态主模型结论。
- **IP-U 的 fdRate boundary 污染明显。** outlier preview 显示多个 `MS-MF-IP-U` 样本在 `fdRate+boundary` 或 `fdRef+fdRate+boundary` 上失败；这支持 IP-U 作为 stress / trade-off 机制，而不是正文主性能路线。

### 6.3 代表性异常格点 / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `MS-MF-IP-U`, seeds `289 / 698`，SNR `0 / 5 / 10` | fdRate boundary | angle over CRB 可达约 `1.8`，`fdRateAbsErrHzPerSec` 约 `1000`，reject reason 多为 `fdRate+boundary`。 | 说明 IP-U 的 unknown-rate 路径在 local-bound 下仍容易被 fdRate 边界控制。 |
| `MS-MF-IP-U`, seed `394`, SNR `-10` | fdRef + fdRate + boundary | `fdRefErrOverCrb ≈ 1.511`，`fdRateAbsErrHzPerSec ≈ 998`。 | 说明 IP-U 的 trimmed 样本选择强，不能只看 trimmed RMSE/CRB。 |
| `MS-MF-CP-U`, seed `347`, SNR `10` | boundary | angle / fdRef 误差接近 0，但 `fdRefBoundaryHit=true`，first-order optimality 非常大。 | 高 SNR CP 仍有边界 / KKT 类 tail，需要 keep rate 伴随报告。 |

## 7. 机制解释

### 7.1 当前解释

本 scan 是 Doppler-aided / in-tooth / truth-centered local MLE 验证。DoA 和 `fdRef` 搜索范围由各自 CRB 倍数控制，因此它回答的是：在 coarse Doppler compensation 和 local bound 已成立时，MLE 是否能达到对应 CRB 量级；它不回答 full acquisition 或 truth-free basin-entry 成功率。

在这个口径下，CP-K/U 的 MLE 与 CRB 基本一致。低中 SNR 下 MSE/CRB 稳定略低于 1，主要是 fixed local box 和 joint-trim 对 tail 的截断；10 dB 下 angle MSE/CRB 回到略高于 1，说明并非所有点都系统性低于 CRB。该结果可以作为“controlled local MLE-vs-CRB consistency”的论文证据，但必须在图注中标明 local / joint-trim 条件。

CP/IP 差异是本次最强的论文机制证据。IP 模型切断跨帧公共相位，导致自身 CRB 大幅变松；因此 IP 在归一化 RMSE/CRB 上也能接近 1，但绝对 DoA / `fdRef` 精度远低于 CP。这个现象比单纯比较 angle RMSE 胜负更适合支撑 continuous-phase 必要性。

### 7.2 这个 scan 支持什么

- 支持在 6 星、10 帧、controlled local-bound 条件下，`MS-MF-CP-K/U` 达到 CRB-level。
- 支持多星 MS 相对 SS 显著降低绝对 DoA CRB 和估计误差，且估计误差随对应 CRB 同步下降。
- 支持 CP 相比 IP 显著保留跨帧连续相位信息量；IP 的 CRB-level 只是相对其更松 CRB 的一致性，不是与 CP 信息量相等。
- 支持 K/U 在当前 CRB-local 主参数上差异较小；unknown-rate 的论文解释应更多结合 EFIM loss、`fdRef` 小幅膨胀、keep rate 和 `fdRate` health，而不是宣称 DoA CRB 显著膨胀。
- 支持 static 路线不作为 dynamic 主模型性能结论：static 在 joint-trim 口径下全空，更适合作反例 / anchor。

### 7.3 这个 scan 不证明什么

- 不证明 full-flow tooth acquisition、subset selection、candidate adoption 或 same-tooth rescue 已经成功。
- 不证明 truth-free estimator 默认路径可以直接达到这些 local-bound 结果。
- 不证明 `MSE/CRB < 1` 表示 CRB 公式错或 estimator 超过 CRB；更合理解释是 bounded local / trimmed 口径。
- 不证明 IP 是合适主模型；IP 只是在自身更松 CRB 下达到一致性。
- 不证明 unknown-rate 没有信息损失；本 scan 只说明 trimmed 主参数差异小，信息损失仍应通过 EFIM / CRB loss 和 keep-rate 行为说明。
- 不应据此写 regression 契约；这是 paper-facing controlled-local scan，而不是自动 pass/fail guard。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。该结果不改变 objective、residual、初始化、reference-sat 语义、method list 或默认 solver。 |
| flow 默认路径 | 不改。不进入 full-flow selector、rescue、candidate adoption 或 final winner。 |
| replay 下一步 | 不需要为本结论新增 replay；如要解释个别 outlier，可用 fixed-seed replay 定点查 IP-U / 高 SNR boundary。 |
| scan 下一步 | 若论文需要 sensitivity，可补 `(2.5, 2.25)` 小 repeat；否则本 snapshot 可作为 CleanBoundConsistency 代表结果。 |
| regression | 不写。local-bound paper-facing curve 不是稳定自动契约。 |
| 论文图 | 可作为 main / appendix 候选：主文建议突出 CP-K/U 的 estimate-vs-CRB 和 CP/IP 的绝对 CRB 差异；IP-U keep-rate 可放副图或表。 |
| 排障记录 | 可摘入主记录：CleanBoundConsistency 全 SNR 固定 `(2,2.5)` 下，MS-MF-CP-K/U 达到 CRB-level，CP/IP 绝对信息量差异显著。 |

## 9. 限制与禁止解释

- 不要把本结果解释为 full-flow acquisition 成功。
- 不要把 truth-centered local hard bounds 用于 runtime selector、gate、candidate adoption 或 final winner。
- 不要只画 trimmed RMSE/CRB 而不报告 keep rate；尤其是 IP-U。
- 不要把低中 SNR 的 MSE/CRB 略低于 1 写成 estimator 超过 CRB。
- 不要用 trace CRB 作为主判断；本结果主口径使用 spherical angular CRB。
- 不要把 static raw 中的 angle 表现当作 dynamic 主模型结论；static 在 health / joint-trim 下全空。
- 不要把 IP 的 CRB-level 归一化结果解释为 IP 与 CP 信息量相当。
- 不要把 K/U angle CRB 相同解释为 `fdRate` nuisance 没有进入模型；应结合 `fdRef` CRB、EFIM loss 和 keep-rate 行为判断。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfMsMleCrbCleanBoundConsistency_20260516-135405.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfMsMleCrbCleanBoundConsistency.m`
```

只运行 `Summary output and plotting` 小节，即可重出 raw / health / joint-trim aggregate、SS-vs-MS、K-vs-U、estimate-vs-CRB、range audit、outlier preview、runtime 表和图。estimate-vs-CRB 图已支持点击 legend 项临时显示 / 隐藏曲线。

## 11. 历史备注

- 早期文档只是 `pending first run` 占位；已由本次 `20260516-135405` full-SNR snapshot 取代。
- `scanMfMsMleCrbCleanBoundScale_20260515-223012.mat` 用于选择本次固定 bound `(DoA, fdRef) = (2, 2.5)`；本结果是该 bound 下的 full-SNR follow-up。
- 本次结果比 CleanBoundScale surface 更适合作论文曲线，但仍是 controlled-local / truth-centered local-bound 口径，不替代 full-flow stress test。
