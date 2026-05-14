# scanMfGeoSatGeometryCrbAudit 结果记录

> 本文档按 `test/dev/scan/results/template/scanResultTemplate.md` 记录 `scanMfGeoSatGeometryCrbAudit.m` 的代表性 coarse-L4 结果。本文只记录 scan 结果、snapshot 绑定、曲线口径、机制解释和当前决策影响；不作为脚本使用手册，不复制完整运行日志，不替代固定地理位置下的 `scanMfPairGeometryCrbAudit.m`。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing scenario screening / coarse representative` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfGeoSatGeometryCrbAudit_20260514-223422.mat` |
| 当前一句话结论 | 在本轮固定经度 `36.59 deg`、纬度 coarse sweep 中，`[45,36.59,0]` 是当前最适合作为下一步 fixed-location paper-facing scan 的主候选；`[60,36.59,0]` 虽然 score 最高，但 coverage 少、EFIM condition 高且 coupling 仍强，更适合作 stress / diagnostic；原 `[37.78,36.59,0]` 降级为 anchor / stress 对照。 |
| 论文图定位 | setup justification / appendix scenario screening；不是最终 MLE-vs-CRB 性能图。 |
| 决策影响 | 推进到 `scanMfPairGeometryCrbAudit` 对 `[45;36.59;0]` 做完整 `L=1/2/4/6/8/10` 深挖；可选 `[55;36.59;0]` 作备选；不进入 regression。 |
| 下一步动作 | 固定地理位置后重跑 single-location satellite-set audit，再进入 MLE/CRB paper-facing scan；暂不扩大为全球经纬网格。 |
| 禁止误用 | 不要把 `geoPreScore` 第一名直接写成论文主场景；不要把 coarse-L4 CRB-health 当作 MLE 已贴 CRB；不要把二星 hard-coupled 解释成多星模型失败；不要把本结果写成选址 / 调度优化贡献。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfGeoSatGeometryCrbAudit.m`
- 结果文档：`test/dev/scan/results/scanMfGeoSatGeometryCrbAudit.md`
- scan 类型：`paper-facing scenario screening / multi-location geometry-CRB health scan`
- 主要问题：在一组 Starlink-inspired 候选地理位置中，哪个 `usrLla` 和哪个初始 cooperative satellite set 更适合作为后续论文仿真起点。
- 扫描对象：纬度候选 `[37.78,25,35,40,45,50,53,55,60] deg`，经度固定 `36.59 deg`，高度 `0 m`；每个位置做 cooperative reference / pair 粗筛，只对 top 3 位置贪心扩到 `L=4`。
- 不覆盖范围：不跑 MLE；不证明 estimator default；不做全球位置优化；不遍历全星座组合；不输出最终 `L=6/8/10` nested set；不比较 rescue / flow / candidate adoption 策略。
- truth 使用口径：未使用 truth 做选择；只使用 TLE、可见性、geometry、CRB-health、EFIM condition、DoA-fdRef coupling 等离线指标。
- 是否 paper-facing：Yes，限于 scenario-screening / setup justification；最终论文性能图需要 fixed-location scan 和 MLE/CRB scan 支撑。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `usrLla` | 候选用户地理位置 `[lat, lon, alt]`。本轮只扫纬度，经度固定 `36.59 deg`。 | No | 用于粗筛 Starlink 可见性和 cooperative geometry。 | 不要解释为全球最优位置搜索。 |
| `numAvailableAtRef` | 参考时刻满足可见性 / access 约束的候选卫星数量。 | No | 是 coverage 下限指标。数量过少的位置即使 score 高也更偏 stress。 | 不要只按可见星数量选主场景。 |
| `geoPreScore` | 地理位置 coarse score，综合 best cooperative reference / pair、可见星数量、EFIM condition 与 DoA-fdRef coupling。 | No | 只用于位置初筛排序。 | 不要把最高 `geoPreScore` 直接等同于论文主场景。 |
| `paperSetScore` | 对 top-geo 的 selected set 做 paper-candidate 排序的辅助分数。 | No | 用来组织 `paperCandidateSetTable`，需要结合 raw CRB / coupling / condition 判断。 | 不要忽略 `paperCandidateClass`。 |
| `selectedRefSatIdxGlobal` | 该位置下 cooperative-aware reference selection 选出的参考星 global index。 | No | 表示该位置的 reference-Doppler 参数化 anchor。 | 不要和最大仰角 reference 混同。 |
| `bestPairSatIdxGlobal` | selected reference 下粗筛得到的 best two-sat partner。 | No | 用于判断 minimal cooperation 是否可用。 | 不要把 best pair 当最终多星 set。 |
| `pairClass="hard-coupled"` | pair / set 的 EFIM condition 或 DoA-fdRef coupling 仍偏强。 | No | 表示可以作为 hard baseline 或 stress case。 | 不等价于 CRB / 模型错误。 |
| `paper-promising-candidate` | coarse scan 中较适合进入下一步 fixed-location deep scan 的候选。 | No | 可作为下一轮 `scanMfPairGeometryCrbAudit` 入口。 | 不表示 MLE 已经贴 CRB。 |
| `stress-or-diagnostic-candidate` | 分数可能高，但 coverage、coupling、EFIM condition 或解释性不够稳。 | No | 保留为 stress / appendix / robustness 候选。 | 不建议作为唯一 main scenario。 |
| `coarseL4` | Geo scan 默认只筛到 `L=4`。 | No | 低成本判断地理位置是否值得完整深挖。 | 不证明 `L=6/8/10` 趋势。 |
| `fdRefDoaCanonicalCorrCpU` / `fdRefCouplingFractionCpU` | `CP-U` EFIM 中 reference-Doppler 与 DoA 子空间的耦合指标。 | No | 越大说明 `fdRef` marginal gain 越容易被 DoA-fdRef coupling 抵消。 | 不要单独用它们否定多星协作。 |
| `efimCondCpU` | `CP-U` 主参数 EFIM 条件数。 | No | 大条件数表示几何 / 参数化下局部信息很尖锐、病态。 | 不等于严格奇异；不能单独判定不可用。 |

常见 scan 口径说明：

- `full-sample`：本 scan 不做 Monte Carlo repeat，因此没有 full-sample RMSE。
- `resolved-sample`：本 scan 不跑 MLE，因此没有 resolved / outlier 过滤。
- `truth-tooth / oracle range`：本 scan 不使用 truth tooth 或 oracle range。
- `stress-test`：`60 deg N` 这类结果可作为 stress / diagnostic candidate，但不应作为主场景直接写入论文性能图。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfGeoSatGeometryCrbAudit_20260514-223422.mat` | 2026-05-14 | `representative` | `geoScanMode="coarseL4"`；9 个纬度候选、经度固定 `36.59 deg`；每个 geo 评估 16 个 reference candidate、12 个 cooperative reference、每个 reference 8 个 second-sat candidate；top 3 geo 贪心到 `L=4`；`coopAll` 合并同一 geo 下的 cooperative pair task。 | `[45,36.59,0]` 是 main deep-scan candidate；`[60,36.59,0]` 是 stress / diagnostic；原 `[37.78,36.59,0]` 排名靠后，降级为 anchor / stress 对照。 | 覆盖早期未按 `coarseL4` / `coopAll` 固定的初版 geo-scan 结果。 |

维护说明：当前代表性结果只保留这一行；若后续加入经度 sweep 或固定位置 deep scan，应在对应 scan 结果文档另建代表性 snapshot，不用覆盖本 coarse-L4 记录。

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed / context seed = 253`
- `seedList = N/A`，本 scan 不跑 Monte Carlo repeat。
- `numRepeat = N/A`，本 scan 是 geometry / CRB-health task scan。
- `snrDbList = 0 dB`
- `frameCountList = 10`
- 关键 scan 轴：`usrLla`，纬度 `[37.78,25,35,40,45,50,53,55,60] deg`，经度固定 `36.59 deg`。
- 关键模型口径：`signalModelTag="constant-doa-affine-doppler-cp"`，即 DoA 恒定、Doppler 一阶变化、continuous phase。
- `tleFileName = "statlink_20260318.tle"`
- `utc0 = "2026-03-18 17:08:00.000000"`
- `numFrame = 10`
- `frameIntvlSec = 1/750`
- `geoSelectionMode = "cooperativeBestScore"`
- `geoScanMode = "coarseL4"`
- `refSelectionMode = "cooperativeBestScore"`
- 关键 resolved / outlier 判据：N/A；不跑 estimator，不做 resolved filtering。
- checkpoint：completed and cleaned；checkpoint summary 中各 run 均 `isComplete=true`。
- snapshot 保存变量：`scanData`
- 运行时间：约 54 分钟，从 `21:40:30` 到 `22:34:21`。

### 4.2 存档数据检查

- 顶层变量：snapshot 日志显示 `scanData` 被保存，`1/91 vars saved`。
- `scanData` 字段：`geoCandidateTable`、`geoReferenceSummaryTable`、`referenceCandidateTable`、`refCooperativeScoreTable`、`allReferencePairRankTable`、`geoGreedyStepTable`、`geoGreedyCandidateTable`、`geoSatSetSummaryTable`、`paperCandidateSetTable`、`bestGeoSummaryTable`、`checkpointSummaryTable`、`plotData`。
- 未保存大体量数据：`rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace、图片。
- warning / fail 计数：运行日志未显示影响结果的 warning 或 failed task；checkpoint summary 全部 complete。
- runtime 结构：每个 geo 的 `coopAll` 为一次约 96 task 的大 `parfor`，典型耗时约 4 分钟；top 3 geo 的 greedy `L=2/3/4` 仍是小 batch。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

本 scan 不产生 MC samples，因此 `samples` 表示该位置在本轮粗筛中进入 geo-level summary 的 location 数，`primary metric` 取 `geoPreScore` 或 `paperSetScore`。`median/P95/P99/resolved/outlier` 不适用于该 scan，用 `N/A` 标注。

| 条件 / case | samples | primary metric | median | P95 | P99 | resolved rate | outlier rate | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Geo coarse sweep | 9 locations | `geoPreScore` best = `-3.6287` at `[60 36.59 0]` | N/A | N/A | N/A | N/A | N/A | score 第一是 stress candidate，不直接作为 main scenario。 |
| Top paper candidate | 3 selected geos × `L=1/2/4` | best `paperSetScore = 4.9590` at `[60 36.59 0]`, L=4 | N/A | N/A | N/A | N/A | N/A | paper class 为 `stress-or-diagnostic-candidate`。 |
| Main recommended candidate | 1 selected geo/set | `[45 36.59 0]`, L=4, `paperSetScore = 4.1002` | N/A | N/A | N/A | N/A | N/A | 唯一 `paper-promising-candidate`，建议进入 fixed-location deep scan。 |
| Current anchor | 1 selected geo | `[37.78 36.59 0]`, `rankByGeoScore=8/9`, `geoPreScore=-4.3836` | N/A | N/A | N/A | N/A | N/A | 当前坐标不再作为唯一主场景。 |

### 5.2 按扫描轴汇总

| axis value (`usrLla`) | case | visible sats | selected ref / best pair | best-pair `angleCrbDegCpK` | best-pair `fdRefCrbHzCpU` | best-pair `efimCondCpU` | best-pair coupling | 解释 |
|---|---|---:|---|---:|---:|---:|---:|---|
| `[60 36.59 0]` | geo rank 1 | 16 | `9147 / 6493` | `0.0057009` | `0.075343` | `9.49e6` | `0.46412` | score 高但 coverage 少，stress 候选。 |
| `[55 36.59 0]` | geo rank 2 | 44 | `5259 / 1243` | `0.0023302` | `0.075355` | `1.8577e6` | `0.46428` | 备选，L=4 后 condition 较温和但 coupling 仍偏高。 |
| `[45 36.59 0]` | geo rank 3 | 56 | `7791 / 1211` | `0.0019300` | `0.075345` | `5.2778e6` | `0.46415` | 主候选，L=4 后最均衡。 |
| `[53 36.59 0]` | geo rank 4 | 52 | `4437 / 5150` | `0.0016614` | `0.075352` | `1.0159e7` | `0.46424` | hold。 |
| `[40 36.59 0]` | geo rank 5 | 56 | `4117 / 8859` | `0.0021474` | `0.075339` | `7.9082e6` | `0.46406` | hold。 |
| `[25 36.59 0]` | geo rank 6 | 36 | `6120 / 4826` | `0.0020734` | `0.075354` | `3.0134e7` | `0.46428` | hold。 |
| `[35 36.59 0]` | geo rank 7 | 35 | `3831 / 9062` | `0.0029558` | `0.075344` | `1.4328e7` | `0.46413` | hold。 |
| `[37.78 36.59 0]` | geo rank 8 | 52 | `4527 / 4520` | `0.0021687` | `0.075350` | `4.578e7` | `0.46422` | current anchor / stress 对照。 |
| `[50 36.59 0]` | geo rank 9 | 55 | `5652 / 7488` | `0.0016182` | `0.075340` | `1.7864e7` | `0.46408` | hold。 |

Top-3 selected-location `L=1/2/4` summary：

| rank | `usrLla` | L | selected set | `angleCrbDegCpU` | `fdRefCrbHzCpU` | `efimCondCpU` | `fdRefCouplingFractionCpU` | class / 定位 |
|---:|---|---:|---|---:|---:|---:|---:|---|
| 1 | `[60 36.59 0]` | 1 | `9147` | `0.023957` | `0.078000` | `186.91` | `8.15e-09` | single-sat baseline |
| 1 | `[60 36.59 0]` | 2 | `[9147 6493]` | `0.0057027` | `0.075343` | `9.49e6` | `0.46412` | hard-coupled |
| 1 | `[60 36.59 0]` | 4 | `[9147 6493 1388 8171]` | `1.1637e-05` | `0.053812` | `6.6594e9` | `0.47476` | stress / diagnostic |
| 2 | `[55 36.59 0]` | 1 | `5259` | `0.0062175` | `0.078000` | `648.12` | `1.04e-08` | single-sat baseline |
| 2 | `[55 36.59 0]` | 2 | `[5259 1243]` | `0.0023304` | `0.075355` | `1.8577e6` | `0.46428` | hard-coupled |
| 2 | `[55 36.59 0]` | 4 | `[5259 1243 348 5878]` | `2.3045e-05` | `0.050110` | `2.7439e8` | `0.39427` | backup / stress |
| 3 | `[45 36.59 0]` | 1 | `7791` | `0.0055570` | `0.078000` | `796.56` | `1.45e-09` | single-sat baseline |
| 3 | `[45 36.59 0]` | 2 | `[7791 1211]` | `0.0019302` | `0.075345` | `5.2778e6` | `0.46415` | hard-coupled |
| 3 | `[45 36.59 0]` | 4 | `[7791 1211 4138 8947]` | `1.0957e-05` | `0.045529` | `9.0602e8` | `0.26624` | main deep-scan candidate |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Geo coarse score preview | `usrLatDeg` | `geoPreScore` / rank | 9 个 `usrLla` | Appendix / setup | 只用于位置初筛；score 第一的 60 deg N 需按 stress 解释。 |
| Coverage / reference preview | `usrLatDeg` | `numAvailableAtRef`、selected ref elevation | 9 个 `usrLla` | Appendix / setup | coverage 少的位置不宜作为唯一主场景。 |
| L=4 CRB-health preview | selected top geos | `angleCrbDegCpU`、`fdRefCrbHzCpU`、`efimCondCpU`、coupling fraction | top 3 geos | Appendix / setup | 不是 MLE 性能图；不含 L=6/8/10。 |
| Paper candidate set preview | candidate set rank | `paperSetScore` 与 raw health columns | top 3 geos × `L=1/2/4` | Internal / setup | 必须同时显示 `paperCandidateClass`。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `[45,36.59,0]` 的 L=4 set `[7791 1211 4138 8947]` 是本轮唯一 `paper-promising-candidate`。它的 `fdRefCrbHzCpU=0.045529`、`fdRefCouplingFractionCpU=0.26624`，在 top 3 中比 60 deg N 和 55 deg N 更均衡。
- `[60,36.59,0]` 的 `geoPreScore` 和 `paperSetScore` 最高，但 `numAvailableAtRef=16`、L=4 `efimCondCpU=6.6594e9`、`fdRefCouplingFractionCpU=0.47476`，因此更适合作 stress / diagnostic，而不是 main scenario。
- `[37.78,36.59,0]` 的 `rankByGeoScore=8/9`，selected reference 仍为 `4527`、best pair 为 `4520`，但 `geoPreScore=-4.3836`，说明当前坐标不宜继续作为唯一 paper-facing 主场景。
- 所有 geo 的 best pair 都是 `hard-coupled`，`numHealthyPair=0`，说明 minimal two-sat cooperation 在当前 reference-Doppler 参数化下普遍存在 DoA-fdRef coupling，不应强行把二星写成最漂亮主结论。
- L=4 已经能区分地理位置质量：45 deg N 的 coupling fraction 降到 `0.26624`，而 60 deg N 仍为 `0.47476`，55 deg N 为 `0.39427`。
- `coopAll` 合并后，9 个地理点完整 coarse run 约 54 分钟完成；每个 geo 的 cooperative audit 约 4 分钟，说明当前 scan 已接近“位置初筛”而不是完整单位置深挖。

### 6.2 反向、污染或未解决现象

- `geoPreScore` 对 60 deg N 偏乐观，说明当前 score 对 angle CRB gain 的奖励仍可能压过 coverage / condition / coupling 惩罚；最终场景必须看 raw table 与 `paperCandidateClass`。
- `selectedRefCoopClass` 在 summary 中仍显示 `pending`，说明 reference class 字段尚未作为稳定分类口径沉淀；当前结论不依赖该字段。
- `coarseL4` 每个 reference 只看 8 个 second-sat candidate，可能漏掉低 coupling partner；该 scan 是初筛，不证明任何位置的全局最优 cooperative set。
- Top-3 greedy 只跑到 `L=4`，不能外推 `L=6/8/10` 是否继续改善；必须由 fixed-location `scanMfPairGeometryCrbAudit` 深挖。
- 本 scan 不含 MLE，所以不能判断 estimator 是否能达到 CRB。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `[60 36.59 0]`, L=4 | stress | `angleCrbDegCpU=1.1637e-05` 很强，但 `fdRefCrbHzCpU=0.053812`、`efimCondCpU=6.6594e9`、`couplingFraction=0.47476`。 | 说明 score 第一不等于 main scenario。 |
| 所有 geos, best L=2 pair | hard-coupled | `pairClass=hard-coupled`，coupling fraction 约 `0.464`，`numHealthyPair=0`。 | 支持二星只作为 minimal hard baseline。 |
| `[37.78 36.59 0]` | anchor degradation | `rankByGeoScore=8/9`。 | 支持当前坐标降级为 anchor / stress 对照。 |
| `[45 36.59 0]`, L=4 | promising | `fdRefCrbHzCpU=0.045529`，`couplingFraction=0.26624`，唯一 `paper-promising-candidate`。 | 支持将 45 deg N 推进到 fixed-location scan。 |

## 7. 机制解释

### 7.1 当前解释

当前 scan 支持一个与前期 CRB audit 一致的机制：在 reference-Doppler 参数化下，二星协作容易把新增信息放到 DoA-fdRef 混合方向上。也就是说，第二颗星能够改善 angle CRB，但 `fdRef` marginal CRB 和 DoA-fdRef coupling 未必同步改善。因此所有候选地理位置的 best pair 都是 `hard-coupled`，二星不适合作为“多星协作已经充分健康”的主证据。

L=4 后，多个非参考星的几何方向开始提供更丰富的差分 Doppler / local steering 约束。不同位置在 L=4 的 coupling、fdRef CRB 和 EFIM condition 之间出现可区分差异。`[45,36.59,0]` 虽然 `geoPreScore` 不是第一，但 L=4 指标更均衡：`fdRefCrbHzCpU` 明显优于 60 / 55 deg N，coupling fraction 也降到 `0.26624`。这使它比 score 第一但 stress 特征明显的 60 deg N 更适合作为主仿真入口。

60 deg N 的结果说明：高纬位置可能带来强 angle geometry，但 coverage 少、condition 极大、coupling 未降时，不能只看 angle CRB 或综合 score。它适合作为 robustness / stress case，用于解释场景选择对 CRB-health 的影响。

### 7.2 这个 scan 支持什么

- 支持将 `[45;36.59;0]` 作为下一步 fixed-location `scanMfPairGeometryCrbAudit` 的主候选。
- 支持将 `[60;36.59;0]` 保留为 stress / diagnostic candidate，而不是 main scenario。
- 支持将原 `[37.78;36.59;0]` 降级为 anchor / stress 对照。
- 支持二星 minimal cooperation 普遍 hard-coupled 的判断；后续论文主场景应优先看 `L=4` 及以上。
- 支持 geo scan 默认采用 `coarseL4`：L=4 已足以筛出更值得深挖的地理位置，同时避免把 geo scan 变成完整选星调度。

### 7.3 这个 scan 不证明什么

- 不证明 45 deg N 是全球最优位置。
- 不证明 `[7791 1211 4138 8947]` 是最终论文 set。
- 不证明 L=6/8/10 会延续 L=4 趋势。
- 不证明 MLE 能贴 CRB。
- 不证明二星协作理论失败。
- 不证明可以把位置 / 选星性能写成 regression 契约。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。该 scan 只做 geometry / CRB-health 场景筛选。 |
| flow 默认路径 | 不改。未涉及 estimator flow、rescue、candidate adoption 或 MLE objective。 |
| replay 下一步 | 暂不需要 replay；先对候选位置做 fixed-location scan。 |
| regression | 不写。地理位置和选星性能不是稳定自动契约。 |
| 论文图 | 可作为 setup justification / appendix scenario screening；最终 main performance figure 仍需 MLE/CRB scan。 |
| 排障记录 | 主记录可摘一句：当前 coarse-L4 geo scan 将 45 deg N 标为主候选，60 deg N 为 stress candidate，原 37.78 deg N 降级为 anchor / stress 对照。 |

## 9. 限制与禁止解释

- 不要用 `geoPreScore` 单独决定最终论文位置；必须结合可见星数、EFIM condition、coupling fraction、`paperCandidateClass` 和后续 fixed-location scan。
- 不要把 coarse-L4 selected set 写成最终 `L=6/8/10` nested satellite set。
- 不要把 CRB-health 选星结果解释为 estimator 已能贴 CRB。
- 不要因为二星 hard-coupled 就否定多星 continuous-phase 模型；多星主价值可能出现在 `L=4` 及以上。
- 不要把 60 deg N 当作 main scenario，除非后续 fixed-location scan 和 MLE/CRB scan 都证明它在更高 L 下健康。
- 不要把本 scan 迁移为 regression。
- 不要把该 scan 包装成 Starlink 全局覆盖 / 调度优化；它只是 paper-facing 场景初筛。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfGeoSatGeometryCrbAudit_20260514-223422.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfGeoSatGeometryCrbAudit.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 11. 历史备注

- 早期 geo scan 近似对每个地理位置跑完整 single-location audit，第一个地理点约 31 分钟，不适合作为地理位置初筛。
- 当前代表性结果切换为 `coarseL4`，并使用 `coopAll` 合并同一 geo 下的 cooperative pair task；9 个地理位置约 54 分钟完成。
- 本结果只用于选出下一步固定位置深挖入口；若 `[45;36.59;0]` 在完整 `L=1/2/4/6/8/10` scan 中不健康，应回到 geo scan 扩大 `maxCandidateSecondSat` 或增加经度 sweep。
