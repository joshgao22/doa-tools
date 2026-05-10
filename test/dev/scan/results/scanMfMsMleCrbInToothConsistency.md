# scanMfMsMleCrbInToothConsistency 结果记录

> 本文档记录 MS-only in-tooth MLE/CRB scan 的 1000-repeat 结果，用于固定各 SNR 下的错误类型分布和代表性 tail seed，指导 `replayMfMsMleCrbFlowDiagnose` 的后续定点排查。本文档只记录 scan 结果和 replay 选 seed 依据；不改变 estimator、flow selector、runtime adoption 或 regression 契约。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / replay-seed index` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfMsMleCrbInToothConsistency_20260510-062927.mat` |
| 当前一句话结论 | MS-MF 不是单一故障：`CP-K` 主要由 `fdref-branch-tail`、`doa-basin-limited`、`coherence-collapse-tail` 三类 tail 污染；`CP-U` 在低中 SNR 主要是 `fd-rate-tail`，在高 SNR 主要转为 `fdref-branch-tail`。 |
| 论文图定位 | 暂不作为 paper-facing MS-MF consistency 图；当前用于 replay 机制排查和 seed 索引。 |
| 决策影响 | 固定 1000-repeat 错误类型统计；下一步推进 `replayMfMsMleCrbFlowDiagnose` 的 scan-tail-driven targeted replay。 |
| 下一步动作 | 按本文档的 representative seed 列表分别排查 fdRef branch、DoA basin、coherence collapse、CP-U fdRate release。 |
| 禁止误用 | 不要把当前 MS-MF full / resolved aggregate 解释为已达到 CRB；不要把错误类型标签下沉到 estimator / flow runtime selector；不要用 truth label 做 candidate adoption。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfMsMleCrbInToothConsistency.m`
- 结果文档：`test/dev/scan/results/scanMfMsMleCrbInToothConsistency.md`
- scan 类型：`diagnostic scan / replay seed selection`
- 主要问题：统计 `MS-MF-CP-K/U` 在 `-15:5:10 dB` 下的 CRB-normalized tail 类型，为 replay 机制排查提供稳定 seed。
- 扫描对象：SNR、`MS-MF-CP-K`、`MS-MF-CP-U`、full / resolved / core / CRB-local 样本、错误类型和代表 tail seed。
- 不覆盖范围：不比较 SS；不运行 replay heavy bank / rescue solver；不验证 estimator default 已修复；不定义 regression pass/fail。
- truth 使用口径：truth 只用于 offline label、CRB-normalized error、resolved / tail 类型归类和 representative seed selection，不进入 runtime selector。
- 是否 paper-facing：No。后续若 MS-MF 经过 replay 修复或形成稳定 resolved 口径，可另开 paper-facing MS scan。

## 2. 术语与统计口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `resolvedRate` | solver-valid、in-tooth、无 frequency-boundary hit；`CP-U` 还要求 fdRate healthy。 | Eval only | 判断是否进入基本 local / resolved regime。 | 不等于 CRB-local，也不保证 angle / fdRef 已贴 CRB。 |
| `coreResolvedRate` | resolved 样本再要求 multi-sat non-ref coherence floor 达标。 | Eval only | 区分 non-ref coherence collapse 与其他 tail。 | 不代表 estimator 主路径已无坏盆地。 |
| `crbLocalRate` | 当前等价于 trimmed core；core 样本再按固定 CRB-normalized angle / fdRef cap 保留。 | Eval only | paper-facing local-regime 候选样本率；当前 MS 的该比例很低。 | 不要只看 crbLocal RMSE/CRB 而忽略 keep rate。 |
| `fdref-branch-tail` | fdRef CRB-normalized error 主导，常伴随 `outside-resolved-tooth` 或 same-tooth fdRef branch。 | Eval only | 优先查 fdRef branch / tooth 边界 / final selection。 | 不等同于 non-ref coherence collapse。 |
| `fd-rate-tail` | `CP-U` 的 fdRate unresolved 或 fdRate error 主导。 | Eval only | 优先查 unknown-rate warm-anchor release / candidate family / final winner。 | 不应解释为 known/unknown EFIM loss 本身。 |
| `coherence-collapse-tail` | non-ref coherence floor 很低，通常 `resolved-coherence-low`。 | Eval only | 优先查 non-ref sat 的 local DoA / fdLocal / phase chain 是否进入似然。 | 不要与 high-coherence 的 DoA basin tail 混合。 |
| `doa-basin-limited` | fdRef / fdRate 相对健康，coherence 通常接近 1，但 angle CRB-normalized error 很大。 | Eval only | 优先查 same-tooth DoA / local-state bad basin。 | 不要继续优先放宽 fd 范围。 |
| `representativeSeedTable` | 每个 `displayName × SNR × error type` 取 tail score 最大的 3 个 seed。 | Eval only | replay 优先从这里挑 case。 | 不定义 regression 固定样本契约。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfMsMleCrbInToothConsistency_20260510-062927.mat` | 2026-05-10 | `representative` | `P=10`，`SNR=-15:5:10 dB`，`numRepeat=1000`，`baseSeed=253`，`seedList=253:1252`，只跑 `MS-MF-CP-K/U`。 | 1000-repeat 固定 MS-MF 错误类型分布；可作为 replay seed 索引。 | 当前唯一保留的代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:1252`
- `numRepeat = 1000`
- `snrDbList = [-15, -10, -5, 0, 5, 10]`
- `frameCountList = 10`
- active methods：`MS-MF-CP-K`, `MS-MF-CP-U`
- `initMode = auto`
- `frameIntvlSec = 1/750`
- `windowSec = 0.012`
- `oracleFdHalfToothFraction = 0.40`
- `oracleFdRateHalfWidthHzPerSec = 1000`
- `resolvedToothHalfWidthFraction = 0.25`
- `resolvedFdRateAbsTolHzPerSec = 250`
- `coreCoherenceFloor = 0.8`
- `trimNormCap = 5`
- task count：`6000`
- checkpoint：enabled；runKey `f10_snrm15to10_seed253to1252_rep1000_6ea0c2f1`；正常完成后已清理 checkpoint artifacts。
- snapshot 保存变量：`scanData`
- 运行时间：约 `6 h 35 m 49 s`

### 4.2 存档数据检查

- 顶层 snapshot 内容：`scanData`
- `scanData` 主要字段：`config`, `perfTable`, `aggregateTable`, `crbLocalSummaryTable`, `failureSummaryTable`, `msErrorTypeSummaryTable`, `representativeSeedTable`, `topTailTable`, `topTailExportTable`, `repeatOutCell`, `checkpointSummaryTable`, `plotData`。
- 未保存大体量数据：`rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace、图片。
- warning / fail 计数：本轮命令行输出未显示导致 scan 中断的错误；后续 MATLAB idle timeout 和 software OpenGL warning 与 scan 数值无关。

## 5. CRB-normalized compact summary

### 5.1 `MS-MF-CP-K`

| SNR (dB) | resolvedRate | coreResolvedRate | crbLocalRate | outlierRate | freqBoundaryHitRate | angle RMSE/CRB (CRB-local) | fdRef RMSE/CRB (CRB-local) | 诊断 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| -15 | 0.963 | 0.798 | 0.108 | 0.037 | 0.002 | 2.2373 | 1.2218 | resolved 很高，但 CRB-local 只剩约 11%；主要被 DoA / fdRef / coherence tail 削掉。 |
| -10 | 0.917 | 0.740 | 0.086 | 0.083 | 0.015 | 2.6559 | 1.2859 | CRB-local 比例下降，fdRef branch 和 DoA basin 都明显。 |
| -5 | 0.901 | 0.727 | 0.091 | 0.099 | 0.024 | 3.7464 | 1.8898 | CRB-local ratio 已偏高，不适合作 paper-facing MS consistency。 |
| 0 | 0.894 | 0.716 | 0.030 | 0.106 | 0.033 | 2.9245 | 2.1319 | 进入 CRB-local 的样本很少，full/core 被 tail 主导。 |
| 5 | 0.908 | 0.679 | 0.007 | 0.092 | 0.033 | 2.5418 | 1.3054 | `crbLocalCount=7`，只作为好样本存在性参考。 |
| 10 | 0.889 | 0.646 | 0.014 | 0.111 | 0.036 | 2.5374 | 2.3316 | 高 SNR 下 CRB-local rate 仍低，tail 随 SNR 不消失。 |

### 5.2 `MS-MF-CP-U`

| SNR (dB) | resolvedRate | coreResolvedRate | crbLocalRate | outlierRate | freqBoundaryHitRate | fdRateResolvedRate | angle RMSE/CRB (CRB-local) | fdRef RMSE/CRB (CRB-local) | 诊断 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| -15 | 0.249 | 0.241 | 0.045 | 0.751 | 0.147 | 0.254 | 2.6617 | 1.0711 | 低 SNR 主要是 fdRate release / boundary 问题。 |
| -10 | 0.271 | 0.264 | 0.048 | 0.729 | 0.109 | 0.279 | 3.3582 | 1.2655 | fdRate unresolved 仍是主导。 |
| -5 | 0.309 | 0.305 | 0.024 | 0.691 | 0.049 | 0.319 | 3.9121 | 1.2685 | fdRate tail 下降但仍大，fdRef branch 开始上升。 |
| 0 | 0.365 | 0.356 | 0.002 | 0.635 | 0.028 | 0.378 | 2.5665 | 0.2298 | `crbLocalCount=2`，ratio 不稳定；不能解释为已贴 CRB。 |
| 5 | 0.465 | 0.456 | 0.004 | 0.535 | 0.019 | 0.480 | 1.8340 | 2.3018 | `crbLocalCount=4`；高 SNR 主要问题已转向 fdRef branch / DoA basin。 |
| 10 | 0.488 | 0.474 | 0.003 | 0.512 | 0.033 | 0.495 | 2.9343 | 0.9511 | `crbLocalCount=3`；只说明好样本存在，不能形成曲线结论。 |

## 6. 错误类型统计

### 6.1 `MS-MF-CP-K` 错误类型分布

| SNR (dB) | crb-local | fdref-branch-tail | doa-basin-limited | coherence-collapse-tail | solver-conditioning-tail | 主导机制 |
|---:|---:|---:|---:|---:|---:|---|
| -15 | 108 (10.8%) | 322 (32.2%) | 417 (41.7%) | 153 (15.3%) | 0 | DoA basin 最大，fdRef branch 次之。 |
| -10 | 86 (8.6%) | 306 (30.6%) | 449 (44.9%) | 159 (15.9%) | 0 | DoA basin 最大。 |
| -5 | 91 (9.1%) | 322 (32.2%) | 424 (42.4%) | 163 (16.3%) | 0 | DoA basin / fdRef branch 双主导。 |
| 0 | 30 (3.0%) | 365 (36.5%) | 448 (44.8%) | 157 (15.7%) | 0 | DoA basin 与 fdRef branch 同时严重。 |
| 5 | 7 (0.7%) | 422 (42.2%) | 372 (37.2%) | 199 (19.9%) | 0 | fdRef branch 成为第一主导。 |
| 10 | 14 (1.4%) | 488 (48.8%) | 291 (29.1%) | 206 (20.6%) | 1 (0.1%) | 高 SNR 下 fdRef branch 明显主导。 |

### 6.2 `MS-MF-CP-U` 错误类型分布

| SNR (dB) | crb-local | fd-rate-tail | fdref-branch-tail | doa-basin-limited | coherence-collapse-tail | 主导机制 |
|---:|---:|---:|---:|---:|---:|---|
| -15 | 45 (4.5%) | 667 (66.7%) | 88 (8.8%) | 192 (19.2%) | 8 (0.8%) | fdRate release 主导。 |
| -10 | 48 (4.8%) | 614 (61.4%) | 119 (11.9%) | 212 (21.2%) | 7 (0.7%) | fdRate release 主导。 |
| -5 | 24 (2.4%) | 443 (44.3%) | 254 (25.4%) | 275 (27.5%) | 4 (0.4%) | fdRate、fdRef branch、DoA basin 三者混合。 |
| 0 | 2 (0.2%) | 105 (10.5%) | 550 (55.0%) | 334 (33.4%) | 9 (0.9%) | fdRef branch 转为主导。 |
| 5 | 4 (0.4%) | 25 (2.5%) | 564 (56.4%) | 398 (39.8%) | 9 (0.9%) | fdRef branch + DoA basin 主导。 |
| 10 | 3 (0.3%) | 9 (0.9%) | 607 (60.7%) | 367 (36.7%) | 14 (1.4%) | fdRef branch 主导，DoA basin 次之。 |

### 6.3 Resolved failure reason 统计

#### `MS-MF-CP-K`

| SNR (dB) | resolved | outside-resolved-tooth | resolved-coherence-low | solver-invalid | 诊断 |
|---:|---:|---:|---:|---:|---|
| -15 | 798 (79.8%) | 37 (3.7%) | 165 (16.5%) | 0 | resolved 多，但 coherence tail 不能忽略。 |
| -10 | 740 (74.0%) | 83 (8.3%) | 177 (17.7%) | 0 | outside-tooth 与 coherence-low 同时存在。 |
| -5 | 727 (72.7%) | 99 (9.9%) | 174 (17.4%) | 0 | 与 -10 dB 类似。 |
| 0 | 716 (71.6%) | 106 (10.6%) | 178 (17.8%) | 0 | coherence-low 稳定存在。 |
| 5 | 679 (67.9%) | 92 (9.2%) | 229 (22.9%) | 0 | coherence collapse 比例升高。 |
| 10 | 646 (64.6%) | 110 (11.0%) | 243 (24.3%) | 1 (0.1%) | 高 SNR 下 coherence-low 和 outside-tooth 都没有自然消失。 |

#### `MS-MF-CP-U`

| SNR (dB) | resolved | fdRate-unresolved | frequency-boundary-hit | outside-resolved-tooth | resolved-coherence-low | 诊断 |
|---:|---:|---:|---:|---:|---:|---|
| -15 | 241 (24.1%) | 570 (57.0%) | 142 (14.2%) | 39 (3.9%) | 8 (0.8%) | 低 SNR 主要卡在 fdRate / boundary。 |
| -10 | 264 (26.4%) | 573 (57.3%) | 107 (10.7%) | 49 (4.9%) | 7 (0.7%) | fdRate unresolved 稳定为最大问题。 |
| -5 | 305 (30.5%) | 596 (59.6%) | 46 (4.6%) | 49 (4.9%) | 4 (0.4%) | fdRate unresolved 仍最大。 |
| 0 | 356 (35.6%) | 566 (56.6%) | 23 (2.3%) | 46 (4.6%) | 9 (0.9%) | 虽然 error type 转为 fdRef branch，但 resolved-rule 仍被 fdRate 卡住。 |
| 5 | 456 (45.6%) | 471 (47.1%) | 15 (1.5%) | 49 (4.9%) | 9 (0.9%) | resolved 与 fdRate-unresolved 接近对半。 |
| 10 | 474 (47.4%) | 455 (45.5%) | 28 (2.8%) | 29 (2.9%) | 14 (1.4%) | fdRate-unresolved 仍未消失。 |

## 7. Replay seed 索引

> 使用建议：`replayMfMsMleCrbFlowDiagnose` 先不要全量加重；按下表每类挑 1–3 个 seed 做 targeted replay。若同一类 seed 现象一致，再扩到整组。表中每格为 `taskSeed`，来自 `representativeSeedTable` 每类 tail score 最大的 3 个 seed。

### 7.1 `MS-MF-CP-K` 推荐 replay seed

| SNR (dB) | fdref-branch-tail | doa-basin-limited | coherence-collapse-tail | solver-conditioning-tail |
|---:|---|---|---|---|
| -15 | `574, 457, 610` | `295, 318, 1234` | `324, 444, 1240` | - |
| -10 | `308, 693, 849` | `372, 845, 1185` | `663, 437, 741` | - |
| -5 | `960, 575, 787` | `1172, 1220, 849` | `1038, 1059, 1148` | - |
| 0 | `1051, 675, 1037` | `297, 405, 498` | `621, 1075, 679` | - |
| 5 | `1085, 392, 557` | `1018, 667, 888` | `423, 1113, 1187` | - |
| 10 | `888, 955, 808` | `697, 466, 511` | `689, 295, 305` | `1097` |

### 7.2 `MS-MF-CP-U` 推荐 replay seed

| SNR (dB) | fd-rate-tail | fdref-branch-tail | doa-basin-limited | coherence-collapse-tail |
|---:|---|---|---|---|
| -15 | `820, 1220, 356` | `1224, 1018, 1131` | `1155, 1068, 438` | `627, 923, 287` |
| -10 | `265, 744, 259` | `976, 1240, 1134` | `1200, 1208, 449` | `899, 745, 518` |
| -5 | `1018, 1047, 396` | `1113, 1040, 761` | `946, 999, 526` | `1150, 505, 1210` |
| 0 | `1054, 574, 838` | `939, 266, 602` | `1230, 271, 423` | `722, 965, 957` |
| 5 | `1208, 1115, 776` | `954, 1039, 1173` | `1167, 585, 567` | `699, 444, 463` |
| 10 | `639, 290, 1232` | `941, 356, 1070` | `1181, 472, 1096` | `820, 677, 1045` |

## 8. 面向 `replayMfMsMleCrbFlowDiagnose` 的排查计划

### 8.1 第一组：fdRef branch / tooth-boundary tail

优先 case：

- `MS-MF-CP-K`: `(-10,308)`, `(0,1051)`, `(5,1085)`, `(10,888)`
- `MS-MF-CP-U`: `(-10,976)`, `(0,939)`, `(5,954)`, `(10,941)`

建议 replay 输出：

```text
fdRefTrue / fdRefInit / fdRefFinal
fdRefErrHz / fdRefErrOverTooth / fdRefInitMoveHz
objectiveFinal / objectiveAtTruthFdRef
objectiveAtFinalDoaTruthFd / objectiveAtTruthDoaFinalFd
localFdRefGridBestErrHz / localFdRefGridBestObjMinusFinal
failureReason / tailSubtype / boundary flag
```

判断目标：确认 `fdRef≈285–300 Hz` 的 branch 是 oracle box 边界、comb tooth branch、final selection 问题，还是 objective 层真实更优分支。

### 8.2 第二组：same-tooth DoA / local-state bad basin

优先 case：

- `MS-MF-CP-K`: `(-15,295)`, `(-10,372)`, `(-5,1172)`, `(0,297)`, `(5,1018)`, `(10,697)`
- `MS-MF-CP-U`: `(-15,1155)`, `(-10,1200)`, `(-5,946)`, `(0,1230)`, `(5,1167)`, `(10,1181)`

建议 replay 输出：

```text
baseline final
truth-doa-final-fd
final-doa-truth-fd
truth-doa-truth-fd
static-doa / wide-center candidate
finalMinusInitAngleDeg
nonRefCoherenceFloor
per-sat localDirErrDeg
```

判断目标：确认 high-coherence / healthy-fd 条件下 DoA 是否被局部 basin 拉偏；优先验证是否存在 truth-free wide/static center 能救回，而不是继续放宽 fd。

### 8.3 第三组：non-ref coherence collapse

优先 case：

- `MS-MF-CP-K`: `(-15,324)`, `(-10,663)`, `(-5,1038)`, `(0,621)`, `(5,423)`, `(10,689)`
- `MS-MF-CP-U`: `(-15,627)`, `(-10,899)`, `(-5,1150)`, `(0,722)`, `(5,699)`, `(10,820)`

建议 replay 输出：

```text
satIdx / isRefSat
finalObjectiveSat / truthObjectiveSat / finalMinusTruthObjectiveSat
finalCoherence / truthCoherence
truthDoaFinalFdCoherence / finalDoaTruthFdCoherence
finalLocalDirErrDeg / finalFdLocalErrHz
per-sat residual or profile score
```

判断目标：区分 coherence collapse 是 DoA/local-state 问题、non-ref Doppler chain 问题，还是 per-sat phase / signal construction 问题。

### 8.4 第四组：CP-U fdRate release

优先 case：

- `(-15,820)`, `(-10,265)`, `(-5,1018)`, `(0,1054)`, `(5,1208)`, `(10,639)`

建议 replay 输出：

```text
fdRateTrue / fdRateInit / fdRateFinal / fdRateErrHzPerSec
fromCpKObjective / fromStaticObjective
mainWarmAnchorObjective / mainWarmAnchorFixedDoaObjective
bestCandidateTag
candidateFdRateErrList / candidateObjectiveList
fdRateBoundaryHit / fdRateResolved flag
```

判断目标：确认 bad case 中好候选是否存在但未被 final winner 采用；若存在，优先做 replay-level adopted row；若不存在，再考虑扩 release family。

## 9. 当前机制解释

### 9.1 `CP-K` 结论

`CP-K` 的 `fdRateResolvedRate=1`，但 CRB-local rate 只有 `0.7%–10.8%`。这说明 known-rate 分支并不是 rate 问题，而是三类 MS-MF tail 混在一起：

1. `fdref-branch-tail`：高 SNR 下成为主导，代表 fdRef branch / tooth 边界 / final winner 问题；
2. `doa-basin-limited`：低中 SNR 与 0 dB 附近主导，代表 same-tooth DoA / local-state bad basin；
3. `coherence-collapse-tail`：稳定存在，尤其 5–10 dB 达到约 `20–24%`，代表 non-ref sat 没有被相干利用。

### 9.2 `CP-U` 结论

`CP-U` 的主问题随 SNR 转移：

- 低中 SNR：`fd-rate-tail` 占主导，`-15/-10/-5 dB` 分别约 `66.7% / 61.4% / 44.3%`；
- 0 dB 以后：`fdref-branch-tail` 成为主导，`0/5/10 dB` 分别约 `55.0% / 56.4% / 60.7%`；
- `doa-basin-limited` 在全 SNR 仍保持 `19.2%–39.8%`，说明即使 rate release 改善，DoA basin 仍会限制 MS-MF；
- `coherence-collapse-tail` 对 CP-U 比例较低，但代表 seed 仍值得抽查，以免被 CP-K 问题掩盖。

### 9.3 对论文口径的影响

当前结果不适合作为 MS-MF paper-facing MLE-vs-CRB 主图。更合理的使用方式是：

- 在 results 文档中固定 MS tail 分类和 representative seeds；
- 用 replay 机制排查和 candidate adoption 先证明好候选能否无 truth 地进入 final result；
- 稳定后再重跑 paper-facing MS scan；
- 如果 full-sample 仍受 tail 污染，论文中必须同时报告 resolved / CRB-local rate 和 outlier type，而不是只画筛选后 RMSE/CRB。

## 10. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。当前 scan 是分类和 seed index，不足以下沉新策略。 |
| flow 默认路径 | 不改。truth label 只用于 offline 分类，不进入 runtime selector。 |
| replay 下一步 | 修改 `replayMfMsMleCrbFlowDiagnose`，增加 scan-tail-driven targeted replay preset / summary 表；先跑代表 seed，不全量加 heavy bank。 |
| regression | 不写。错误类型分类和 representative seed 不是稳定 pass/fail 契约。 |
| 论文图 | 暂缓 MS-MF consistency 图；当前用于解释 MS outlier 机制和确定后续 replay 路线。 |
| 排障记录 | 主记录可摘一句：`scanMfMsMleCrbInToothConsistency_20260510-062927` 固定 1000-repeat MS-MF 错误类型分布，MS-MF 下一步转入 targeted replay。 |

## 11. 限制与禁止解释

- 不要把 `crbLocalAngleRmseOverCrb` 看起来接近 `2–4x` 的局部好样本直接当成 paper-facing 结论；当前 `crbLocalRate` 太低。
- 不要把 `CP-U` 的 `fd-rate-tail` 直接解释为 EFIM nuisance information loss；它更像 solver / release / candidate family 健康问题。
- 不要因为 `fdref-branch-tail` 很多就 blanket 放宽 fd 范围；本文档的 replay 目标是判断 branch 来源，而不是先改搜索域。
- 不要把 `representativeSeedTable` 迁移成 regression 固定样本；它只是机制 replay 的索引。
- 不要把 truth-derived `msErrorType` 用于 runtime selector、gate、candidate adoption 或 final winner。

## 12. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfMsMleCrbInToothConsistency_20260510-062927.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
test/dev/scan/scanMfMsMleCrbInToothConsistency.m
```

只运行 `Summary output and plotting` 小节，即可重出 compact table、错误类型统计、representative seed 和图。

## 13. 历史备注

- 本轮结果把该 scan 的定位固定为 MS-MF replay seed index；不是 paper-facing MS consistency 曲线。
- 本文档只保留 `scanMfMsMleCrbInToothConsistency_20260510-062927.mat` 的 1000-repeat 代表性结果；不再保留旧 snapshot 的统计或结论。
