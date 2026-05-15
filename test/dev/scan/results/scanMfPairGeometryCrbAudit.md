# scanMfPairGeometryCrbAudit 结果记录

> 本文档按 `test/dev/scan/results/template/scanResultTemplate.md` 记录 `scanMfPairGeometryCrbAudit.m` 的代表性运行结果。本文只记录 scan 结果、snapshot 绑定、曲线口径、机制解释和当前决策影响；不复制完整运行日志，不替代脚本 README，也不把几何 / CRB-health 结论解释为 MLE 已贴近 CRB。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing scenario screening / fixed-location satellite-set audit / representative` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfPairGeometryCrbAudit_20260515-112711.mat` |
| 当前一句话结论 | 在 `[45,36.59,0]` 与 `[55,36.59,0]` 两个候选位置中，`[55,36.59,0]` 的 `L=6` 组合 `[5259 1243 348 5652 14 4437]` 是当前最适合作为论文主仿真的多星组合；它兼顾明显 CRB 增益、`healthy-gain` 分类、低 scaled condition 和低 DoA-`fdRef` coupling。 |
| 论文图定位 | setup justification / scenario screening；可支撑后续 MLE-vs-CRB 主仿真的场景选择，但本身不是性能曲线。 |
| 决策影响 | 后续 paper-facing MLE / CleanBound / in-tooth consistency 优先固定 `usrLla=[55,36.59,0]`、`refSat=5259`、`selectedSatIdxGlobal=[5259 1243 348 5652 14 4437]`。 |
| 下一步动作 | 用 `loc2-L6` 跑 MS-MF `CP-K/CP-U` 的 controlled local MLE-vs-CRB；用 `loc2-L4` 做轻量 smoke；用 `loc2-L8` 做多星数增强确认。 |
| 禁止误用 | 不要把本 scan 写成全球选址 / 星座调度优化；不要把 hard-coupled 的纯 CRB 最小组合直接当论文主场景；不要把 CRB-health 解释为 estimator 默认路径已通过。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfPairGeometryCrbAudit.m`
- 结果文档：`test/dev/scan/results/scanMfPairGeometryCrbAudit.md`
- scan 类型：`fixed-location geometry / CRB-health / satellite-set audit`
- 主要问题：在 `scanMfGeoSatGeometryCrbAudit` 粗筛出的地理位置上，进一步判断哪个 reference satellite 和 nested multi-satellite set 更适合后续论文仿真。
- 扫描对象：两个用户位置 `[45,36.59,0]` 和 `[55,36.59,0]`；每个位置用 cooperative reference selection 选择参考星，再按 `satCountList=[1 2 4 6 8 10 12 14 16]` 贪心扩展卫星集合。
- 不覆盖范围：不运行 MLE；不验证 full-flow；不比较 rescue / adoption 策略；不枚举所有卫星组合；不证明最终 estimator RMSE 达到 CRB。
- truth 使用口径：不使用 truth 做 runtime 选择；本 scan 只基于 TLE、几何、CRB / EFIM、condition、coupling 和相位残差等离线指标。
- 是否 paper-facing：Yes，限于场景与卫星组合选择；最终论文性能图仍需 MLE/CRB scan 或 replay 支撑。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---:|---:|---|---|
| `locationIdx` / `usrLlaStr` | 候选用户位置，本轮为 `[45,36.59,0]` 与 `[55,36.59,0]`。 | No | 用于比较固定位置下的可见星、reference 和 cooperative geometry。 | 不要解释为全球最优位置。 |
| `selectedRefSatIdxGlobal` | cooperative reference selection 选出的参考星 global index。 | No | 定义 reference-satellite Doppler 参数化 anchor。 | 不要和最大仰角星简单混同。 |
| `selectedSatIdxGlobalStr` | 贪心构造得到的 nested satellite set。 | No | 后续仿真可直接使用的多星集合。 | 不代表组合枚举全局最优。 |
| `angleCrbDegCpK / angleCrbDegCpU` | `CP-K/CP-U` 下 DoA CRB 标准差，单位 deg。 | No | 越小表示几何上 DoA 信息越强。 | 不等于实际 MLE RMSE。 |
| `fdRefCrbHzCpK / fdRefCrbHzCpU` | `CP-K/CP-U` 下参考星 `fdRef` CRB 标准差，单位 Hz。 | No | 越小表示 reference Doppler 信息越强。 | 不说明 full-flow tooth acquisition 已解决。 |
| `unknownOverKnownFdRefPct` | `100*(fdRefCrbU/fdRefCrbK-1)`。 | No | unknown-rate nuisance 对 `fdRef` CRB 的相对回退。 | 不是实际 CP-U estimator 的 RMSE 回退。 |
| `efimScaledCondCpU` | 主参数 EFIM 经尺度归一后的 condition。 | No | 用于判断 DoA 与 `fdRef` 信息是否均衡；越小越适合解释。 | 不能替代 raw CRB。 |
| `efimScaledMinEigDoaFracCpU / efimScaledMinEigFdFracCpU` | scaled 最弱信息轴中 DoA / `fdRef` 分量占比。 | No | 接近 `0.5/0.5` 表示最弱轴较均衡；极端接近 1 或 0 表示信息轴偏向单一参数。 | 不要单独作为好坏判据。 |
| `fdRefDoaCanonicalCorrCpU` | DoA 子空间与 `fdRef` 的 canonical coupling 指标。 | No | 越低表示 `fdRef` 与 DoA 的耦合越弱，CRB 增益更容易解释。 | 不等价于估计误差相关系数。 |
| `fdRefCouplingFractionCpU` | `fdRef` 信息被 DoA coupling 抵消的比例。 | No | 越低越说明 `fdRef` gain 没被 coupling 吃掉。 | 不要把高 coupling 直接解释成模型错误。 |
| `pairClass` | 根据 CRB gain、condition、coupling、phase residual 得到的诊断分类。 | No | `healthy-gain` 更适合作 paper-facing 主例；`hard-coupled` 更偏 stress / appendix。 | 不要把 `hard-coupled` 当作严格不可用。 |
| `phaseResidualRmsRad` | 几何 Doppler affine approximation 的相位残差 RMS。 | No | 当前所有候选均在 `1e-6 rad` 量级，说明 affine-Doppler phase 模型本身不是主要限制。 | 不要用它证明 MLE 已收敛。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：N/A，本 scan 不做 Monte Carlo repeat。
- `resolved-sample`：N/A，本 scan 不跑 MLE，因此没有 resolved / outlier 过滤。
- `truth-tooth / oracle range`：N/A。
- `stress-test`：`loc1` 系列与 `loc2-L10+` 更适合作 hard geometry / stress，不建议直接作主论文组合。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfPairGeometryCrbAudit_20260515-112711.mat` | 2026-05-15 | `representative` | `contextSeed=253`；`SNR=0 dB`；`P=10`；`T_f=1/750 s`；`usrLlaList=[45,36.59,0;55,36.59,0]`；`satCountList=[1 2 4 6 8 10 12 14 16]`；`refSelectionMode=cooperativeBestScore`。 | `[55,36.59,0]` 的 `L=6` 组合是当前主推荐；`L=4` 为轻量 smoke；`L=8` 为增强确认；`[45,36.59,0]` 系列整体 hard-coupled。 | 取代 coarse geo scan 中“先深挖 45 deg”的临时优先级；45 deg 降级为 hard-coupled stress 对照。 |

维护规则：当前代表性结果只保留本 snapshot。若后续用新版本脚本重跑并移除 checkpoint housekeeping 字段，只要核心表格和结论一致，可在本节追加新 snapshot 并把本行标为 `superseded-by-clean-snapshot`。

## 4. 最新代表性运行

### 4.1 配置

- `contextSeed = 253`
- `snrDb = 0`
- `numRepeat = N/A`，本 scan 是几何 / CRB-health scan，无 MC repeat。
- `usrLlaList = [45.00,36.59,0; 55.00,36.59,0]`
- `numFrame = 10`
- `frameIntvlSec = 1/750 = 0.0013333 s`
- `refFrameIdx = 6`
- `satCountList = [1 2 4 6 8 10 12 14 16]`
- `numElem = [4 4]`
- `carrierFreq = 11.7 GHz`
- `sampleRate = 512 MHz`
- `minUsrElevationDeg = 15`
- `maxSatOffAxisDeg = 55`
- `refSelectionMode = cooperativeBestScore`
- `maxReferenceCandidateEval = 40`
- `maxCooperativeReferenceEval = 12`
- `maxCandidateSecondSat = 200`
- `maxGreedyCandidatePerStep = 80`
- checkpoint：本 snapshot 来自 checkpoint-enabled run；当前 `.mat` 中仍包含旧版运行期 checkpoint housekeeping 字段，但本文档不使用它们解释结果。
- snapshot 保存变量：`scanData`
- 运行时间：`elapsedSec = 4734.22 s`，约 `78.9 min`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 主要结果字段：`config`、`contextSummaryTable`、`referenceCandidateTable`、`refCooperativeScoreTable`、`pairRankTable`、`greedyStepTable`、`satSetSummaryTable`、`rankPreviewTable`、`locationSummaryTable`、`locationSatSetSummaryTable`、`plotData`、`elapsedSec`。
- 当前 snapshot 中还存在旧版运行期字段：`checkpointSummaryTable`、`checkpointCleanupReport`、`checkpointDir`。这些字段只反映 run housekeeping，不作为结果分析依据；后续按已修正脚本重跑后不需要进入 snapshot。
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map 或图片文件。
- warning / fail 计数：上传 summary 中没有显示任务失败；所有打印出的 context / set summary 状态均可用于本轮结论。

## 5. 主要统计与曲线结果

### 5.1 两个位置的 context 对比

| location | `usrLla` | available sats | selected ref | ref elev min (deg) | SS angle CRB CP-K (deg) | SS fdRef CRB CP-U (Hz) | 解读 |
|---:|---|---:|---:|---:|---:|---:|---|
| 1 | `[45 36.59 0]` | 56 | 7791 | 55.279 | 0.0039294 | 0.055178 | 可见星更多，单星 DoA CRB 略好，但后续 cooperative set 全部 classified 为 `hard-coupled`。 |
| 2 | `[55 36.59 0]` | 44 | 5259 | 55.262 | 0.0043964 | 0.055178 | 单星 DoA CRB 略差、可见星略少，但 `L=4/6/8` 能形成 `healthy-gain` 多星集合。 |

### 5.2 当前推荐组合

| 用途 | location | L | selected satellites | class | angle CRB CP-U (deg) | fdRef CRB CP-U (Hz) | scaled cond CP-U | corr | coupling fraction | 推荐结论 |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 主论文仿真 | 2 | 6 | `[5259 1243 348 5652 14 4437]` | `healthy-gain` | `5.9066e-06` | `0.022773` | `1.3554` | `0.14688` | `0.021572` | **首选**。CRB 增益明显，scaled 最弱轴 DoA / fdRef 分量接近 `0.50/0.50`，DoA-`fdRef` coupling 很低，解释性最好。 |
| 轻量 smoke / debug | 2 | 4 | `[5259 1243 348 5652]` | `healthy-gain` | `2.4928e-05` | `0.033629` | `3.763` | `0.5718` | `0.32695` | 卫星数少、运行便宜；适合先跑 MLE smoke，但 coupling 明显高于 L6。 |
| 增强确认 / 多星数趋势 | 2 | 8 | `[5259 1243 348 5652 14 4437 2372 5878]` | `healthy-gain` | `4.2620e-06` | `0.020047` | `1.8658` | `0.23017` | `0.052976` | CRB 比 L6 进一步下降，但 scaled condition 与 coupling 稍差；适合确认 L 从 6 到 8 的增益是否能被 estimator 吃到。 |
| CRB 上界 / stress | 2 | 16 | `[5259 1243 348 5652 14 4437 2372 5878 1632 7122 9772 2356 9776 5150 8495 49]` | `hard-coupled` | `1.6358e-06` | `0.013917` | `1.3297` | `0.13234` | `0.017513` | 纯 CRB 很强，尤其 `fdRef` 最小；但 raw EFIM condition 触发 hard-coupled，不适合作主场景第一张图。 |
| 角度 CRB 极小 stress | 1 | 16 | `[7791 1211 4138 5675 8947 8173 7417 4306 4079 9003 81 4117 9727 3073 7273 8440]` | `hard-coupled` | `1.0037e-06` | `0.018837` | `5.4966` | `0.68096` | `0.46371` | 纯 DoA CRB 最小，但 coupling 极强；只适合 hard geometry / appendix，不建议做主论文组合。 |

### 5.3 `loc2=[55,36.59,0]` 的卫星数趋势

| L | selected satellites | class | angle CRB CP-U (deg) | fdRef CRB CP-U (Hz) | U/K fdRef rollback (%) | scaled cond CP-U | DoA frac | fdRef frac | coupling fraction | 解释 |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | `5259` | N/A | `4.3964e-03` | `5.5178e-02` | 7.4856 | 1.9933 | 1.0000 | `6.09e-08` | `1.04e-08` | 单星 anchor；`fdRef` 信息主要来自参考链路。 |
| 2 | `[5259 1243]` | `hard-coupled` | `1.6492e-03` | `5.3287e-02` | 3.8339 | 587.3 | 0.99999 | `1.14e-05` | 0.46390 | 二星已有 DoA 增益，但 coupling 和 scaled condition 太强。 |
| 4 | `[5259 1243 348 5652]` | `healthy-gain` | `2.4928e-05` | `3.3629e-02` | 4.8534 | 3.7630 | 0.60256 | 0.39744 | 0.32695 | 第一组可用 healthy set；适合轻量验证。 |
| 6 | `[5259 1243 348 5652 14 4437]` | `healthy-gain` | `5.9066e-06` | `2.2773e-02` | 7.3063 | 1.3554 | 0.50187 | 0.49813 | 0.02157 | 最均衡，当前主推荐。 |
| 8 | `[5259 1243 348 5652 14 4437 2372 5878]` | `healthy-gain` | `4.2620e-06` | `2.0047e-02` | 7.0471 | 1.8658 | 0.76346 | 0.23654 | 0.05298 | CRB 继续改善，但最弱轴更偏 DoA，解释性略弱于 L6。 |
| 10 | `[5259 1243 348 5652 14 4437 2372 5878 1632 7122]` | `hard-coupled` | `2.4121e-06` | `1.8150e-02` | 6.8580 | 2.0284 | 0.98920 | 0.01080 | 0.07580 | CRB 更低，但最弱轴几乎全是 DoA，已不如 L6 均衡。 |
| 12 | `[5259 1243 348 5652 14 4437 2372 5878 1632 7122 9772 2356]` | `hard-coupled` | `1.8257e-06` | `1.6157e-02` | 7.2519 | 1.4213 | 0.58110 | 0.41890 | 0.02803 | raw condition 仍高；可作高 L stress。 |
| 14 | `[5259 1243 348 5652 14 4437 2372 5878 1632 7122 9772 2356 9776 5150]` | `hard-coupled` | `1.7625e-06` | `1.4780e-02` | 7.4477 | 1.1445 | 0.49369 | 0.50631 | 0.00452 | scaled 指标很好，但 raw EFIM condition 触发 hard-coupled；可作上界对照。 |
| 16 | `[5259 1243 348 5652 14 4437 2372 5878 1632 7122 9772 2356 9776 5150 8495 49]` | `hard-coupled` | `1.6358e-06` | `1.3917e-02` | 7.3394 | 1.3297 | 0.50436 | 0.49564 | 0.01751 | 纯 CRB 最强，但不建议第一阶段作为主场景。 |

### 5.4 `loc1=[45,36.59,0]` 的作用

`loc1` 的全部 multi-sat set 均被标为 `hard-coupled`。它的 raw angle CRB 随 L 增大下降很快，例如 `L=16` 达到 `1.0037e-06 deg`，但 `fdRefDoaCanonicalCorrCpU` 长期在 `0.53–0.68`，`fdRefCouplingFractionCpU` 也保持在 `0.28–0.46` 区间。它更适合作 hard geometry / stress 或 appendix 对照，不适合作论文主仿真的首选组合。

### 5.5 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `Angle CRB vs L` | `numSat` | `angleCrbDegCpU` 或 `angleCrbDegCpK` | `loc1`、`loc2` | Appendix / setup | 可以说明增加卫星数带来几何 CRB 下降，但不能说明 estimator 已贴 CRB。 |
| `fdRef CRB vs L` | `numSat` | `fdRefCrbHzCpU` | `loc1`、`loc2` | Appendix / setup | 用于说明多星对参考 Doppler 也有信息增益。 |
| `scaled condition / coupling vs L` | `numSat` | `efimScaledCondCpU`、`fdRefCouplingFractionCpU` | `loc1`、`loc2` | Yes, setup justification | 推荐用于解释为什么 `loc2-L6` 比纯 CRB 最小组合更适合作主例。 |
| `classification vs L` | `numSat` | `pairClass` | `loc1`、`loc2` | Text / table | `healthy-gain` 只出现在 `loc2-L4/L6/L8`，这是主推荐的核心依据。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `loc2-L6` 是当前最均衡的组合：`angleCrbDegCpU=5.9066e-06 deg`、`fdRefCrbHzCpU=0.022773 Hz`、`efimScaledCondCpU=1.3554`、最弱轴 DoA / `fdRef` 分量约 `0.50187/0.49813`、`fdRefCouplingFractionCpU=0.021572`。这比只看 CRB 最小更适合论文主图解释。
- `loc2` 的 `L=4/6/8` 均为 `healthy-gain`，说明该位置存在一段可解释的 cooperative gain 区间；`L=6` 位于该区间中部，兼顾星数、CRB 和 condition。
- `loc1` 虽然可见星更多、部分高 L raw CRB 更小，但所有多星 set 均为 `hard-coupled`，且 DoA-`fdRef` coupling 长期偏强，不适合作为第一主场景。
- `phaseResidualRmsRad` 均在 `~2e-6–3e-6 rad` 量级，说明 affine-Doppler / continuous-phase 几何近似在这些候选上不是主要限制。
- unknown-rate 对 `fdRef` 的 U/K rollback 在 `loc2` 多星组合中约 `4.85%–7.45%`，是可解释的 nuisance 信息损失，不会改变 L6 作为主组合的判断。

### 6.2 反向、污染或未解决现象

- 纯 CRB 最小不是主推荐：`loc1-L16` 的 angle CRB 最小，`loc2-L16` 的 fdRef CRB 最小，但二者都标为 `hard-coupled`，不宜直接作为论文主图组合。
- `loc2-L14/L16` 的 scaled condition 和 coupling fraction 很好，但 raw EFIM condition 仍触发 hard-coupled 分类；它们可以做上界 / stress，不应第一阶段替代 L6。
- 本 scan 没有 estimator MC，因此不能判断 `MS-MF-CP-K/U` 在这些组合上是否进入正确 basin，也不能判断 same-tooth bad basin 是否缓解。
- 本轮只比较两个候选位置，不构成全球场景最优；如果论文需要更强的场景选择说法，必须另跑更大地理网格或明确写成 representative scenario。

### 6.3 代表性异常格点 / set

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `loc1-L16` | raw CRB optimum / hard geometry | angle CRB `1.0037e-06 deg`，但 `fdRefCouplingFractionCpU=0.46371`、`fdRefDoaCanonicalCorrCpU=0.68096`。 | 说明“CRB 最小”不等于“最适合论文主仿真”。 |
| `loc2-L2` | minimal cooperation hard pair | angle CRB 降到 `0.0016492 deg`，但 scaled condition `587.3`、coupling fraction `0.46390`。 | 说明二星协作在此场景仍偏 hard-coupled，主价值出现在 L=4 以后。 |
| `loc2-L10` | hard-coupled transition | angle / fdRef CRB 继续改善，但 scaled 最弱轴 DoA fraction `0.9892`、fdRef fraction `0.0108`。 | 说明继续加星后可能出现信息轴不均衡，L6/L8 更适合作主解释。 |
| `loc2-L14/L16` | high-L upper bound | scaled condition 接近 1，coupling fraction 很低，但 pairClass 仍为 `hard-coupled`。 | 可用于上界 / stress confirm，不宜第一阶段作为 main scenario。 |

## 7. 机制解释

### 7.1 当前解释

这个 scan 说明，多星集合选择不能只看 DoA CRB 或 `fdRef` CRB 的绝对最小值。对本文的 reference-satellite continuous-phase 模型来说，论文主场景还需要 EFIM 的最弱信息轴不要过度偏向单一参数，并且 DoA 与 reference `fdRef` 的 Schur-coupling 不要过强。否则，即使 raw CRB 很小，后续 MLE 也更容易表现为 hard geometry、参数耦合或局部 basin stress，而不是干净的 cooperative gain。

`loc2-L6` 的优势正好在于它不是单一指标极小，而是各项指标均衡：DoA CRB 相比单星下降约三数量级，`fdRef` CRB 也从 `0.055178 Hz` 降到 `0.022773 Hz`；同时 scaled condition 只有 `1.3554`，最弱轴几乎等分 DoA / `fdRef`，coupling fraction 仅 `0.021572`。这使它更适合论文中解释“多星几何带来主参数有效信息增益”，而不是把结果归因于病态耦合或极端几何。

相比之下，`loc1` 系列虽然 raw angle CRB 很强，但几乎所有二星和多星 set 都带有强 DoA-`fdRef` coupling 或高 condition。它适合作为 stress / appendix，用于说明不合适的 cooperative geometry 会把信息增益变成 hard-coupled 结构；但不适合作主文第一条 MLE-vs-CRB 曲线。

### 7.2 这个 scan 支持什么

- 支持将 `[55,36.59,0]` 作为下一步 fixed-location paper-facing 主场景。
- 支持将 `[5259 1243 348 5652 14 4437]` 作为当前最优的主仿真多星组合。
- 支持 `L=4 -> L=6 -> L=8` 作为轻量、主例、增强确认的 nested set 选择。
- 支持二星 minimal pair 不足以代表本文的多星 cooperative gain；至少需要 L=4，主推荐为 L=6。
- 支持用 condition / coupling / classification 约束“最优组合”的定义，而不是只按 CRB 最小排序。

### 7.3 这个 scan 不证明什么

- 不证明 `MS-MF-CP-K/U` estimator 在 `loc2-L6` 已经贴近 CRB。
- 不证明 full-flow acquisition、subset tooth selection 或 same-tooth rescue 已经稳定。
- 不证明 `[55,36.59,0]` 是全球最优地理位置。
- 不证明 `L=14/16` 不可用；它们只是目前不适合作为第一阶段主论文组合。
- 不适合写 regression 契约；它是场景筛选和结果记录，不是自动 pass/fail guardrail。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。该 scan 只做 geometry / CRB-health 场景筛选。 |
| flow 默认路径 | 不改。未涉及 tooth selection、rescue、candidate adoption 或 final winner。 |
| replay 下一步 | 建议新增或配置一个 `loc2-L6` 代表性 replay / CleanBound run，先验证 `MS-MF-CP-K/U` 的 local MLE-vs-CRB。 |
| regression | 不写。卫星集合选择和场景筛选不是稳定自动契约。 |
| 论文图 | 可作为 setup justification / scenario selection；最终 main performance figure 仍需 MLE/CRB scan。 |
| 排障记录 | 主记录可摘一句：fixed-location pair audit 将 `[55,36.59,0]` 的 L6 healthy set 作为下一步 MS-MF 主战场，`[45,36.59,0]` 降级为 hard-coupled stress。 |

## 9. 限制与禁止解释

- 不要把本 scan 直接写成“选出了 Starlink 最优多星组合”；更准确说法是“在当前两个候选位置和当前贪心规则下，选出了最适合论文仿真的 representative cooperative set”。
- 不要只用 `angleCrbDegCpU` 或 `fdRefCrbHzCpU` 的最小值定义最优；论文主组合必须同时考虑 `pairClass`、`efimScaledCondCpU`、最弱轴分量和 DoA-`fdRef` coupling。
- 不要把 `hard-coupled` 解释为模型错误或不能用；它表示该组合更偏 stress / appendix，需要后续 MLE 结果确认。
- 不要把 `loc2-L6` 的 CRB-health 解释成 estimator 默认路径已修好；它只说明该场景值得进入 MLE scan。
- 不要把 coarse `scanMfGeoSatGeometryCrbAudit` 中 `[45,36.59,0]` 的临时优先级继续当成最终结论；本轮 fixed-location deep scan 已经把主推荐转向 `[55,36.59,0]`。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfPairGeometryCrbAudit_20260515-112711.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfPairGeometryCrbAudit.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

后续推荐仿真配置可以先固定为：

```matlab
usrLla = [55; 36.59; 0];
selectedRefSatIdxGlobal = 5259;
selectedSatIdxGlobal = [5259 1243 348 5652 14 4437];
numFrame = 10;
frameIntvlSec = 1 / 750;
```

轻量 smoke：

```matlab
selectedSatIdxGlobal = [5259 1243 348 5652];
```

增强确认：

```matlab
selectedSatIdxGlobal = [5259 1243 348 5652 14 4437 2372 5878];
```

## 11. 历史备注

- `scanMfGeoSatGeometryCrbAudit_20260514-223422.mat` 的 coarse-L4 结果曾把 `[45,36.59,0]` 作为优先 deep-scan 候选，并将 `[55,36.59,0]` 作为备选。
- 本轮 fixed-location deep scan 显示 `[45,36.59,0]` 在 L=2 到 L=16 均为 `hard-coupled`，而 `[55,36.59,0]` 在 L=4/6/8 形成连续的 `healthy-gain` 区间，因此当前主优先级应切换到 `[55,36.59,0]`。
- 由于本 snapshot 仍含旧版 checkpoint housekeeping 字段，若后续需要提交仓库结果文档，可用修正后的脚本重跑得到更干净的 snapshot；若数值表不变，本结果结论无需重写。
