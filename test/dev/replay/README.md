# replay 脚本说明

`test/dev/replay/` 中的文件是人工排障用的 section 化脚本，不是 regression，也不是可复用 API。每个脚本固定一个代表性问题或小 Monte Carlo 诊断，用于输出 summary、保存轻量 `replayData`，并根据 `replayData` 中的轻量曲线、曲面或表格数据重新画图。

详细运行结果、长表格、观察现象和 snapshot 绑定统一放在 `test/dev/replay/results/`；本 README 只保留脚本规范、入口索引和当前一句话结论。

## 模板入口

新增 replay 时先复制 `test/dev/replay/template/replayTemplate.m`，再改名到 `test/dev/replay/replayXxx.m`。模板只规范 section、checkpoint/snapshot/Telegram 壳和打印风格，不是可运行实验入口，也不是统一 replay 框架。

模板目录 `test/dev/replay/template/` 不放正式 replay；正式入口仍只放在 `test/dev/replay/` 根目录。

新增或改 replay 时，先复用模板和 `test/common/replay|report|probe|summary` 中已有工程壳；不要在单个 replay 里复制 checkpoint、snapshot、notify、compact table 或 solver-probe 表行逻辑。只把脚本特有配置、局部指标和一次性绘图保留 local。

## 统一脚本格式

每个 replay 文件头部采用固定形式：

```matlab
% English purpose comment.
% English usage / storage comment.

clear; close all; clc;

%% Replay configuration
```

默认参数直接写在 `Replay configuration` 小节中。不要使用 workspace override，也不要新增默认参数 helper。如果要改变 repeat 数、seed、SNR、搜索范围或是否保存 snapshot，直接修改该脚本头部参数。

## 固定 section 顺序

1. `Replay configuration`：显式写本 replay 的默认参数。
2. `Build context and flow options`：构造 context、flow option、run key 与 tmp 路径。
3. `Run replay batch`：执行固定样本或小 Monte Carlo，可在外层使用 `parfor`。
4. `Data storage`：构造轻量 `replayData`，用 `saveExpSnapshot` 只保存 `replayData`，并清理 tmp。
5. `Summary output and plotting`：只依赖 `replayData` 输出 summary 并画图。恢复 snapshot 后可直接运行这一节。
6. `Local helpers`：仅放本脚本私有 glue、summary 和绘图 helper。

## 可选运行通知

长时间 replay 可在顶层脚本结束、失败 `catch` 或 snapshot 保存后调用 `test/common/replay/notifyMfReplayStatus.m` 发送可选 HTML Telegram 通知；该 helper 内部只包装 `utils/io/notifyTelegram.m` 的通用 I/O 壳。通知必须是 best-effort：发送失败只 `warning`，不得改变 `replayData`、snapshot、图表数据、数值路径或任何 replay 结论口径。

通知内容只保留脚本名、状态、耗时、snapshot 文件/目录和少量关键 summary。具体 metric line 由各 replay 本地构造并传入 common helper；common helper 不解析 `replayData` 中的具体结果表，不维护第二套 replay results parser。详细结果仍写入 `test/dev/replay/results/`。

当前已验证但降频使用的第三批 diagnostic replay（`replayMfCombToothSurface`、`replayMfInToothDoaDopplerRidgeTrace`、`replayMfInToothPolishGate`、`replayMfSameToothHardCase`、`replayMfSubsetRankingTrace`、`replayMfRandomRescueEffectiveness`、`replayMfUnknownReleaseRoute`、`replayMfWarmAnchorParforSensitivity`、`replayMfFastStatsRepresentativeDivergence`，以及 flow-like stress replay）已统一接入 common header 与 `notifyMfReplayStatus`。其中固定单样本或轻量诊断默认 `checkpointEnable=false`，不为了格式统一额外创建 per-repeat checkpoint task；若脚本原本已有短期 runDir 用于绘图数据或 cleanup，则保持原样并在完成后清理。已有 per-repeat checkpoint 价值的 flow-like / 大 MC replay 继续保留 checkpoint / resume。

Telegram HTML 格式遵循 `template/replayTemplate.m` 的示例：固定外壳统一，正文指标弹性；`<code>` 只用于短数字、短 tag 或短状态，不包裹完整路径、长 candidate family 名或长错误文本；`commentLineList` 应写 replay-level 结论、recommendation 或下一步判断，而不只是 `completed`。

## 推荐推进顺序

这组 replay 按当前排障主线组织。后续整体仿真、推进和观察现象时，建议按下面顺序运行，而不是按文件名字母顺序运行。

## static anchor 诊断入口

#### `replaySfStaticDualSatDiagnose.m`

- 作用：对 canonical `4154 + 1165` 双星做固定 SNR、小 MC 的 single-frame static replay，作为 `scanSfStaticMleCrbConsistency.m` 的机制级伴随入口，也是当前唯一活跃的 SF static-anchor replay。
- 职责边界：scan 负责 paper-facing 的 SNR sweep 与大 repeat CRB consistency；本 replay 只看固定 SNR 下的 SS/MS static 对比、sat2 weight sweep、代表性 repeat、per-satellite static diagnostic、reference-sat 不变量，以及 compact DoA-only CRB-scale sanity。
- SF-Doppler 解释口径：`SF-Static` 只作为 single-frame constant-Doppler matched-model anchor，不承担 Doppler-dynamic 贡献叙事；Doppler evolution、continuous phase 和 unknown-rate information loss 仍由 MF CP-K/CP-U scan / replay 支撑。
- DoA-only 口径：`SS/MS-SF-DoA` 只保留 compact pilot-model effective-gain CRB sanity，用于防止旧 unit-gain CRB gap 被误读为 estimator failure；不再维护独立 heavy DoA-only CRB audit replay。
- 复用链路：默认通过 `buildSfStaticSingleSnrRepeatFixture` 复用 shared SF static scene / snapshot / estimator flow，通过 `printMfReplayHeader`、`printMfReplaySection`、`buildMfReplayCheckpointOpt`、`dispMfReplayTablePreview`、`notifyMfReplayStatus` 和 `finalizeMfReplayResult` 复用 replay 工程壳；不在 replay 内维护第二套 static estimator。
- 主要输出：`caseSummaryTable`、`staticGainSummaryTable`、`weightPreferenceTable`、`repeatCompareTable`、`worstMsAngleLossTable`、`representativeDetail`、`staticCrbSummaryTable`、`doaOnlyCrbScaleTable`、`referenceInvariantTable` 与三张 compact diagnostic 图。
- 当前结论口径：用于解释 static hard pair / weight / per-sat 机制和 DoA-only CRB signal-scale sanity；不替代 `scanSfStaticMleCrbConsistency` 的正文级 static anchor 曲线，也不用于改 estimator 默认路径。旧 `replaySfMsDoaCrbDiagnose.m` 已移除，历史结果保留在 `results/replaySfMsDoaCrbDiagnose.md`。

### 1. 先确认 comb / tooth 是真实 objective 结构

#### `replayMfCombToothSurface.m`

- 作用：固定样本可视化 `1/T_f` fdRef comb 和 DoA-fdRef surface；当前默认在同一 MS 周期窗口内同时扫描 CP-U (`unknown`) 与 CP-K (`known`) objective。
- 先看现象：wrong-tooth 是否沿 `1/T_f` 重复，truth / `finalByMode` 是否在不同 tooth 上；`finalByMode` 在 CP-U 图中解析为 `finalCpU`，在 CP-K 图中解析为 `finalCpK`。
- 配置口径：头部 `modeTagList` 控制要扫描的 evaluator，`centerNameList` 支持 `truth / finalByMode / finalCpU / finalCpK / finalEstimate / staticSeed`；`finalEstimate` 仅作为旧 CP-U final alias 保留。`surfaceClipDeltaObj` 只控制 clipped heatmap 色阶，不改变 objective 数据。
- 主要输出：`modeTable`、带 `modeTag / resolvedCenterName` 的 fdRef line / DoA-fdRef surface summary、`sliceTable`、folded comb、DoA-fdRef heatmap / mesh、clipped heatmap、best-tooth / center-tooth DoA 与 fdRef slice 图，以及按 mode 输出的 truth-centered DoA-fdRef coupling heatmap。
- CP-K DoA 形状诊断：头部 `latLonSurfaceModeList / latLonSurfaceCenterNameList / latLonSurfaceGridCount / latLonSurfaceHalfWidthDeg / latLonSurfaceFdOffsetTagList` 控制额外 lat/lon 小网格；默认只对 `known + finalByMode` 同时扫描 `centerFd`（固定在当前 center 的 `fdRef`，即 `fdOffset=0`）与 `surfaceMinFd`（固定在 DoA-fdRef surface 选出的最优 fd tooth）。这样可直接比较 wrong-tooth 与 corrected-tooth 上的 CP-K DoA 曲率，判断 known-rate surface 是否只是被 1-D DoA cut 或 fdRef comb 色阶掩盖。输出 `latLonSurfaceTable.fixedFdTag` 与对应 lat/lon heatmap / slice 图。
- 耦合诊断：头部 `couplingFdHalfWidthHz / couplingDoaHalfWidthDeg / couplingFdGridCount / couplingDoaGridCount` 可调；输出 `rhoDoaFd`、ridge slope 和 quadratic ridge slope，用于判断 truth 附近 profile surface 是否近似轴对齐。
- 当前结论：comb / wrong-tooth 是 objective 层真实结构；该 replay 只用于观察 CP-U / CP-K surface 与 final center，不改变 estimator、flow 或默认 rescue 策略。
- 详细结果：`results/replayMfCombToothSurface.md`。

### 2. 再确认 subset 负责选齿，periodic 负责同齿细化

#### `replayMfPeriodicVsSubsetToothSelect.m`

- 作用：对比 periodic wide、selected subset、final periodic refine 的选齿与细化分工。
- 配置：头部只保留 `snrDb / baseSeed / numRepeat / saveSnapshot / notifyTelegramEnable / optVerbose / checkpointEnable / toothHistogramBinCount`，并用英文注释说明；共享场景、卫星、TLE、waveform 与 frame offset 由 `runSimpleDynamicFlowReplayBatch` 调用 `buildDynamicDualSatEciContext` 构造。
- 临时存档：每个 repeat 运行时间较长，头部 `checkpointEnable=true` 时启用 per-repeat checkpoint / resume；设为 `false` 时不创建 tmp、不恢复、不清理。任务文件写入仓库根目录 `tmp/replayMfPeriodicVsSubsetToothSelect/<runKey>/task/task_*.mat`。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存与 checkpoint 清理后通过 `notifyMfReplayStatus` 发送 HTML `DONE` 通知，在失败时发送 `FAILED` 通知；通知只包含脚本名、耗时、snapshot / checkpoint 与少量 tooth-selection 指标，发送失败只 warning。
- 进度条：`Run replay batch` 阶段必须显示外层 repeat progressbar；启用 checkpoint 时，progressbar 只统计未完成的 todo repeat；并行时由 checkpoint runner 的 `parallel.pool.DataQueue` 在 client 侧更新。
- 先看现象：subset 是否缩小 tooth error，periodic refine 是否只在同齿内改善 angle / fd。
- 主要输出：三阶段 `toothIdx`、angle error、fdRef error 的 compare table、aggregate table 与合并 subplot 图。
- 当前结论：subset 负责 tooth selection，periodic full-data 负责同齿 refine；final periodic refine 不能被当作 wrong-tooth rescue。
- 详细结果：`results/replayMfPeriodicVsSubsetToothSelect.md`。

#### `replayMfSubsetRankingTrace.m`

- 作用：解释真实 repeat 下 subset candidate 的排序、trusted reason、selected subset 与 final tag。
- 先看现象：same subset bank 下是否同时出现 truth tooth 与 wrong tooth，ranking 是否有 no-truth-leak 风险。
- 主要输出：`selectedSubsetLabel`、`toothIdx`、candidate rank、subset bank 覆盖。
- 命令行输出：只打印 selection、aggregate、subset bank 和 representative candidate 的紧凑预览；完整表保存在 `replayData`。
- 当前结论：用于解释 subset selection 的 candidate 覆盖与 ranking，不直接固化为 regression。

### 3. 再看 conditional rescue 是否救齿、是否误伤 easy case

#### `replayMfRandomRescueEffectiveness.m`

- 作用：比较 curated-only 与 curated+rescue/random subset bank；当前 rescue 口径为 primary curated bank 加 curated rescue bank，`random1` 仍作为 margin fallback 候选观察。
- 先看现象：`curated3 / random1` 是否只在 hard case 中补齿，是否把 easy / median case 拉坏。
- 主要输出：wrong-tooth rescue、easy-case damage、各 subset label evaluated / selected 次数、`replayData.toothHistogramTable` 中的 `|toothIdx|` 分布表与对应 histogram 图。
- 命令行口径：逐 seed compare table 与 subset bank coverage 只打印预览，完整表保存在 `replayData`；长 histogram table 默认不在命令行打印。
- 当前结论：random rescue 只能作为条件触发的 tooth rescue，不应 blanket 常驻。
- 详细结果：`results/replayMfRandomRescueEffectiveness.md`。

### 4. paper-facing in-tooth MLE / CRB 口径诊断

#### `replayMfCrbTexConsistencyAudit.m`

- 作用：对照论文模型层级检查 CRB case dispatch、interest / nuisance 参数切分、time template、FIM / EFIM block、多星贡献、K/U Schur loss 与 summary 取值是否自洽；默认不跑 MLE，不替代性能 replay。
- 默认场景：使用 paper-facing 6 星 context，`usrLla=[55,36.59,0]`、`selectedSatIdxGlobal=[5259 1243 348 5652 14 4437]`、`refSatIdxGlobal=5259`。SS case 仍取 reference-only view，MS case 使用完整 6 星 view；`buildDynamicDualSatEciContext` 不再作为该 replay 的默认入口。
- 默认模型集合：`SS/MS-SF-DoA`、`SS/MS-MF-Static`、`SS/MS-MF-CP-K`、`SS/MS-MF-CP-U`。SF 只使用 DoA-only pilot-model effective CRB；MF 才进入 DoA-Doppler / Doppler-rate CRB。
- 主要输出：`contextSummaryTable`、`crbModelAuditTable`、`crbPathAuditTable` / `crbPathAuditAggregateTable`、`crbFimBlockAuditTable`、`crbFdRefMarginalAuditTable`、`crbSatContributionAuditTable`、`crbDerivativeNormTable`、`crbCasePairAuditTable`、`crbSummaryReuseAuditTable`。重点观察 `dispatchOk`、`path-ok`、`unexpected-identical`、`reported-fdref-crb-mismatch`、`summary-reuse-suspect` 和 static / CP-K、known / unknown、SS / MS pair 的 FIM 差异层级；`crbFdRefMarginalAuditTable` 专门并排显示 reported `fdRef` CRB、EFIM full inverse、scalar `1/EFIM_fd,fd` 与 DoA-fdRef coupling inflation，用于判断 fdRef CRB 不变是输出复用还是边际化耦合导致。
- 职责边界：这是 CRB 构造诊断 replay；只用于决定是否需要修改 `performance/` CRB helper 或后续固化 regression，不负责 MLE-vs-CRB 曲线、full-flow outlier 或 estimator tuning。
- 详细结果：`results/replayMfCrbTexConsistencyAudit.md`。当前固定结论为：CRB dispatch / reported fdRef / K-U Schur loss 正常，`MF-Static` 是 continuous zero-rate `CP-K` 退化，MS 的 fdRef scalar gain 可被 DoA-fdRef marginal coupling 抵消；该 replay 运行很快，结果文档不依赖固定 snapshot。


#### `replayMfSsMleCrbMetricDiagnose.m`

- 作用：对 `scanMfMleCrbInToothConsistency` 暴露出的单星 SS-SF / SS-MF MLE-vs-CRB gap 做小样本诊断；它复用 scan 的 truth-centered in-tooth `fdRef` / `fdRate` 盒，只跑少量 seed 与 SNR，不承担论文统计曲线。
- 诊断口径：同时输出 spherical angle error、`projectCrbToAngleMetric(latlon)` 投影 CRB、lat/lon 单轴 CRB、trace CRB、requested init / estimator init / dynamic final、`fdRef` init/final/CRB、legacy median-CRB ratio 与 seed-normalized `fdRefNormRmse` / bias ratio、CP-U `fdRate` init/final、objective improve、iteration、first-order optimality、boundary hit 与 warm-anchor candidate tag；默认 `mfInitMode="static-doa-truth-fd"`，即使用 static DoA seed 与 truth-centered frequency seed 观察 local in-tooth MLE；`mfInitMode="internal-estimator"` 只作为 stress 诊断。
- 主要输出：`caseTable`、raw `aggregateTable`、`healthAggregateTable`、`trimAggregateTable`（均保留 legacy median-denominator ratio 与 seed-normalized `fdRef` CRB ratio）、`tailTable`、`objectiveProbeTable`、`solveProbeTable`、`pathProbeTable`、`diagnosisTable`、`crbMetricTable`、`releaseCompareTable`、`rangeTable` 与三张诊断图；命令行打印 raw / health / trim 聚合表、objective probe、controlled solve probe、DoA path probe、诊断决策表、被 health / trim 排除的 tail 预览、CRB metric projection check、CP-K/CP-U release compare 和 oracle range summary，完整逐 seed 结果保存在 `replayData`。
- 存储口径：默认启用 per-task checkpoint / resume，并按 task 数与 Parallel Computing Toolbox 自动启用外层 `parfor`；estimator 内层并行仍显式关闭，任务文件写入仓库根目录 `tmp/replayMfSsMleCrbMetricDiagnose/<runKey>/task/`；`runKey` 使用短语义前缀加 8 位 hash，避免 checkpoint 路径过长；成功构造并保存轻量 `replayData` 后清理 checkpoint。
- 当前定位：作为 paper-facing scan 的伴随流程定位入口；raw aggregate 永远保留所有点，health / trim aggregate 只用于区分 solver/frequency/boundary/no-solve 与 CRB-normalized tail，不能替代 raw 结果。MF angle 未对齐 CRB 时，应优先看 `objectiveProbeTable` 判断 truth / mixed point 在 MF objective 中是否优于 final，再看 `solveProbeTable` 判断 truth-DoA 或 wide-static-start 是否能接住好盆地，最后用 `pathProbeTable` 观察 final/static 到 truth 方向上的 objective 走势；`diagnosisTable` 只给排查分流标签，不作为论文统计。它不修改 estimator、scan 或 regression 默认路径。

#### `replayMfMsMleCrbCleanTrim.m`

- 作用：作为 MS-vs-SS 的干净 truth-centered / in-tooth MLE-vs-CRB replay，也作为当前轻量 clean SNR sweep 入口；用于指导正式估计器在给定真值附近搜索范围内的局部估计、初始化和 trim 口径修改。
- 默认扫描：文件头默认按 `snrDbList=(-10:10:10).' × seedList` 构造 task grid；当前默认 `numRepeat=100`、seed `253:352`、`P=10`、DoA hard box `3×CRB`、MF `fdRef` truth-local `0.4 tooth` half-width。默认场景已切到 `buildDynamicMultiSatEciContext` 的 6 星配置：`usrLla=[55;36.59;0]`、`tleFileName="statlink_20260318.tle"`、`selectedSatIdxGlobal=[5259 1243 348 5652 14 4437]`、`refSatIdxGlobal=5259`。快速 smoke test 时可临时把 `numRepeat` 降到 10 或 20。该脚本已经在 aggregate、compare、search-range audit 与 CRB path audit 中按 SNR 分组；后续更系统的 DoA box scale、fdRef half-width、geometry、frame count 与 information-loss 曲线交给 scan，不再在 CleanTrim 内继续扩展扫描维度。
- 固定路线：默认跑 `SS-SF-DoA / MS-SF-DoA / SS-MF-Static / MS-MF-Static / SS-MF-CP-K / SS-MF-CP-U / MS-MF-CP-K / MS-MF-CP-U`。`SF-DoA` 是 no-Doppler 单帧 DoA-only 失配 baseline，不估计 `fdRef`，对应 `crbPilotSfDoaOnlyEffective` 的 pilot-model effective-gain CRB；`MF-Static` 是 zero-rate 多帧 DoA-Doppler 失配 baseline，仍估计 `fdRef`，并作为后续 `MF-CP-K/U` 的 DoA / `fdRef` seed。`CP-K` 是 known-rate baseline，`CP-U` 是 unknown-rate nuisance 口径；MF 路线均保持 `pure-mle / core-only`，不再默认跑 `SS-MF-CP-U-Robust`。MS-MF 只读取 MS-MF-Static 结果，不读取 SS-MF final。
- MF 口径：文件头 `dynamicObjectiveMode="pure-mle"` 与 `dynamicRouteMode="core-only"` 是当前默认，用于检查最干净的 one-shot MLE core：关闭 CP engineering penalty、fd alias unwrap、DoA basin-entry、known embedded branch 与 unknown warm-anchor。本 replay 还在本地覆盖生成端和 MF 估计 / CRB 为 `cleanMfSignalTimeModel="linearRef"`、`cleanMfSteeringMode="frozenRef"`，即生成端使用固定参考帧 steering 与 reference-tied affine Doppler phase，MF estimator / MF CRB 使用同一 fixed-DoA 口径；默认 context、正式 estimator 和其它 replay/scan 不随之改变。MF 的 DoA hard bound 固定为 truth-centered `truthLocalDoaCrbScale × 当前方法 CRB std`，不随 seed 漂移，也不使用固定绝对角度宽度；若 seed 超出该 hard bound，由 estimator init check 裁回 box 内。此前用于 SS route 对照的 robust 行已从默认 method list 删除；如需复查 basin-entry，只在临时诊断中显式改 route，不混入 CleanTrim 主输出。
- 搜索范围：DoA 不再使用固定绝对 half-width，而是由 `truthLocalDoaCrbScale` 对每个 SNR / 方法的 DoA CRB std 统一缩放得到；`fdRef` 先使用 truth-centered in-tooth requested 范围，再按 `truthLocalFdRefCrbScale × fdRef CRB std` 做 floor 保护。Static 与 MF 分别计算 floor：static 可以被放宽到避免 1σ 截尾，MF 仍优先保持用户给定的 in-tooth 范围，只有当 MF 自身 CRB floor 大于 requested 范围时才放宽；`fdRate` 使用 truth-centered 固定局部范围。`cleanSearchRangeAuditTable` 同时汇总 static 与 MF 的实际 half-width / CRB 最小比值，并报告 `doaTruthBoxViolationRate`，用于确认 final DoA 没有越过固定 truth-centered hard bound。
- 主要输出：`caseTable`、`cleanAggregateTable`、`cleanHealthAggregateTable`、`cleanTrimAggregateTable`、`cleanSsMsCompareTable`、`cleanKnownUnknownCompareTable`、`cleanSnrCurveTable`、`cleanEstimateCrbCurveTable`、`cleanSearchRangeAuditTable`、`cleanStaticCrbMetricAuditTable`、`cleanCrbPathAuditTable`、`cleanCrbPathAuditAggregateTable`、`cleanOutlierTable`、`runtimeTable / runtimeAggregateTable / topSlowRuntimeTable`。其中 `cleanSsMsCompareTable` 按 `MF-Static / MF-CP-K / MF-CP-U` 分族报告 MS/SS 绝对 angle RMSE、MS 相对自身 CRB 的 MSE、MS 相对同族 SS-MF CRB 的 MSE、fdRef MSE/CRB、keep rate 与 `targetClass`；`cleanKnownUnknownCompareTable` 在 raw / health / joint-trim 口径下分别比较 SS 与 MS 的 `CP-K` / `CP-U` angle RMSE、angle MSE/CRB、fdRef MSE/CRB、keep rate 和 K/U class；`cleanStaticCrbMetricAuditTable` 只针对显式启用的 matched `SS/MS-SF-Static` 行，默认 `SF-DoA` 不进入该 matched DoA-Doppler audit；`cleanCrbPathAuditAggregateTable` 汇总 K/U monotonicity、EFIM Schur loss PSD、shared block、reference `fdRef/fdRate` invariant、`gamma_ref` coupling、angle/fdRef CRB U/K ratio 与 `primaryAuditClass`。`cleanEstimateCrbCurveTable` 同时保存 joint-trim 与 raw 两种曲线口径；脚本末尾分别画 trimmed 与 raw 两张 log-y 图，每张图各含 DoA / fdRef 两个子图。DoA 子图按当前 `methodNameList` 显示有 finite angle 数据的方法；`fdRef` 子图强制跳过 `SS/MS-SF-DoA` 这类 DoA-only 行，只画有 `fdRef` 估计和 CRB 的方法，避免 DoA-only 的 fallback CRB 进入 Doppler 图例。
- trim 口径：raw 永远保留全样本；health 只剔除 solver/frequency/rate/boundary/no-solve 不健康样本；joint-trim 使用 `health + angleErr/CRB <= trimAngleNormMax + fdRefErr/CRB <= trimFdRefNormMax`。该 trim 是 replay 后统计，不进入 estimator selector。
- 新增 static CRB source audit：`cleanStaticCrbSourceAuditTable` 拆出 static matched CRB 的 raw interest FIM、Schur-eliminated EFIM、no-nuisance CRB/current CRB、nuisance loss、condition number 与 per-sat EFIM sum consistency；`cleanStaticCrbSatContributionTable` 按卫星拆分 raw/EFIM trace 与 fdRef 信息占比；`cleanStaticCrbSsMsPairAuditTable` 对齐同一 seed 下 SS/MS 的 fdRef CRB、等效信息和 sat2 贡献；`cleanStaticCrbDiffStepAuditTable` 重算多组 `jacOpt.diffStep`，用于判断 lat/lon finite-difference Jacobian 是否引入 CRB 步长敏感性。这些表只读 CRB aux 和 replay case summary，不改变 estimator、CRB helper 默认公式或 objective。
- 详细结果：`results/replayMfMsMleCrbCleanTrim.md`。当前代表性 snapshot 为 `test/data/cache/replay/replayMfMsMleCrbCleanTrim_20260513-181333.mat`；`20260513-172449` 的 `3×CRB` hard-box 版本只作为 sensitivity 证据，不作为代表性曲线。
- 当前定位：该 replay 已固定为 representative / diagnostic-only / fixed-bound local sanity 入口。`replayMfMsMleCrbCleanTrim_20260513-181333.mat` 显示 SS-MF 在 bounded local 口径下稳定贴近 CRB，MS-MF 在 `2×CRB` DoA hard box 内形成绝对 angle RMSE gain，并在部分 SNR 达到 own-CRB level；但 `3×CRB` sensitivity 已显示 MS 会明显远离 own-CRB，因此 `2×CRB` 结果必须解释为强 oracle local guard，而不是自然 basin-entry 能力。K/U CRB path audit 为 `path-ok`，说明 K/U angle CRB 接近不应再作为 CRB 主公式错误处理；unknown-rate information loss 后续应通过 Schur loss、fdRef inflation、主方向 loss 与专门 scan 展示。若只想检查 `MS-SF-Static` 是否真的低于 matched static CRB，可临时把 `methodNameList` 缩到 `"SS-SF-Static" / "MS-SF-Static"` 或只保留 `"MS-SF-Static"`。先看 `cleanStaticCrbMetricAuditTable` 的 `whitenedSecondTraceOverDim`、`fdRefMseOverCrb`、有限 MC 95% 下界字段和 `metricAuditClass`，再看 `cleanStaticCrbSourceAuditAggregateTable`、`cleanStaticCrbSatContributionAggregateTable`、`cleanStaticCrbSsMsPairAggregateTable` 与 `cleanStaticCrbDiffStepAuditAggregateTable`：前两张检查 nuisance Schur loss、per-sat FIM 贡献和 sat2 fdRef 信息是否被消掉，pair 表检查 MS 相对 SS 的 fdRef 等效信息增益是否过小，diff-step 表检查 lat/lon finite-difference Jacobian 是否对步长敏感。不要只看 spherical scalar MSE/CRB。当 `replayMfMsMleCrbFlowDiagnose` 路径过重或表格过多时，优先跑本脚本判断 clean estimator 是否能在统一 CRB 倍数的真值附近范围内贴近 CRB 或至少优于 SS-MF；若仍不达标，先看 `cleanKnownUnknownCompareTable`、`cleanStaticCrbMetricAuditTable` 与 `cleanCrbPathAuditAggregateTable`，再决定是否转向 scan 或回到 estimator 主流程修改，不继续往 replay 中添加候选路径。Summary output / plotting 小节只读取 `replayData` 中已保存的数据段，支持 load snapshot 后重跑。

#### `replayMfMsMleCrbFlowDiagnose.m`

- 作用：对 `scanMfMleCrbInToothConsistency` 暴露出的多星 MS-MF MLE-vs-CRB gap 做小样本流程定位；它不是论文统计入口，而是固定为 **in-tooth MS-MF-CP-U local-basin / nuisance-rate 诊断路线**，用于检查 MS static seed、SS-MF seed、MS-MF final、truth / mixed objective 和 per-satellite final-eval。
- 固定路线：文件头 `runScale` 支持 `"ss-sanity"`、`"ms-crb-audit"`、`"ms-fast"`、`"ms-fixed-snr"`、兼容别名 `"ms-sanity"` 与 `"full"`。当前默认 `"ms-crb-audit"`，用于按 trim MSE/CRB 目标排查 estimator-side MS DoA basin-entry 与 CRB 归一化口径：它启用 `estimatorDoaBasinEntryScope="all"`，关闭 replay-only MS bank，直接观察正式 MF estimator 的 angle/fdRef normalized MSE 是否能接近 1.1。上一轮 `ms-targeted` 会漏 hidden angle-only tail，`ms-tail-expanded` 的 0.024/top-2 axis bank 没有继续压低 residual angle gap 且运行明显变慢，二者不再作为本 replay 的公开路线保留。若要重新做固定 SNR coverage-first bank 上限检查，应把文件头改为 `runScale="ms-fast"`。`ms-fast` 按 `fixedSnrDb / baseSeed / numRepeat` 做固定 SNR 连续 seed 批量 replay，当前默认 `fixedSnrDb=-5`、`baseSeed=253`、`numRepeat=40`，只保留 `SS-SF-Static / MS-SF-Static / SS-MF-CP-U / MS-MF-CP-U` 方法集，并让 static transition bundle 使用 `ms-seed-only` 轻量模式。该路线设置 `solveProbeRoute="ms-bank-only"`、`msBankRoute="cpu-release-axis0p006-0p012-core"`，对每个 MS CP-U repeat 跑 compact axis-cross bank，仅用于 replay-only adoption 上限确认，不改变 estimator winner。`ms-fixed-snr` 保留全量固定 SNR 诊断，方法集为 `SS-SF-Static / MS-SF-Static / SS-MF-CP-K/U / MS-MF-CP-K/U`，并对所有 MS CP-U case 跑 compact bank，用于阶段性 coverage 确认。`ss-sanity` 仍只用于 SS-MF DoA basin-entry wiring 快速复查；`full` 才恢复代表性 seed/SNR 和完整 targeted MS 诊断。SS/MS 动态方法都会显式传入 `estimatorDoaBasinEntryScope`，由 scope 值决定哪些方法实际 evaluate entry；脚本保留 `ssDoaBasinEntrySummaryTable` 作为 compact 护栏，确保 MS-only route 不污染 SS。
- probe helper 边界：fixed-point objective probe 已抽到 `test/common/probe/buildMfObjectiveProbeRows.m`，targeted solve probe 已抽到 `test/common/probe/buildMfTargetedSolveProbeRows.m`，MS CP-U bank probe 已抽到 `test/common/probe/buildMsCpuBankProbeRows.m`。solve row、DoA basin-entry row、per-sat probe row、MF probe opt vector、truth scalar、coherence floor 与 typed empty rows 进一步抽到 `test/common/probe/` 的小 helper，避免 replay 与 bank / solve probe 维护重复构表逻辑。这些 helper 只负责 replay-side fixed-point evaluation / controlled rerun / CP-K preselect -> CP-U release probe 与表行重组，不定义 estimator winner adoption，也不改变 estimator 默认路径；MS 诊断标签、primary diagnosis 与 representative case 仍保留在本 replay 内，方便继续收敛分类口径。
- MS 专用输出：`perSatProbeTable` 记录每颗星的 final coherence、per-sat objective/residual、local direction median/max error、global angle error 与 reference Doppler error，用于定位 non-ref coherence collapse 或 local-state bad basin；`msDiagnosisTable / msDiagnosisAggregateTable` 将 MS-MF 样本拆成 `ms-doa-basin-limited`、`ms-nonref-coherence-collapse`、`ms-fdref-branch-suspect`、`ms-truth-init-rescuable`、`ms-wide-static-insufficient`、`ms-ss-mf-seed-not-helpful`、`ms-cp-u-worse-than-cpk`、`ms-cp-u-objective-better-angle-worse`、`ms-cp-u-rate-release-weak`、`ms-cp-u-rate-coupled-basin`、`ms-kkt-conditioning`、`ms-cp-u-kkt-conditioning` 等标签；`msPrimaryDiagnosisTable / msPrimaryDiagnosisAggregateTable` 在多标签基础上给出每个样本的第一分流类型，方便决定下一步是查 coherence collapse、fdRef branch、DoA basin，还是 solver conditioning；`diagnosisTable` 仍保留通用 MF 分流标签，不作为论文统计。
- 候选追踪口径：`msCandidateTraceTable` 将 MS baseline、controlled solve probe、estimator 内部 `doaBasinEntry` acquisition 与 objective fixed-point probe 统一整理为 candidate trace，包含 candidate source、DoA / fd seed 来源、DoA half-width、objective delta、angle / fd 误差、angle improvement / CRB-normalized improvement、是否优于 baseline 与潜在 damage 标记。对 MS bank rows 还会解析 `bankCenterSource / bankAxisDirection / bankStepDeg`，并补充 `fdRef` over-tooth、nearest-tooth distance 与 candidate-to-baseline tooth shift 诊断字段，用于判断改善是否仍在 in-tooth 同齿内；这些字段不参与 runtime selector。`damageFlag` 同时考虑 objective 变差和 angle 变差，只作为 replay-level 诊断，不改变 estimator adoption。
- estimator 内部 basin-entry 口径：`msDoaBasinEntryTable` 直接展开 `optimInfo.doaBasinEntry`，记录 baseline 与 controlled probe 中 `doaBasinEntry*` 是否触发、触发原因、候选数、best tag、baseline / entry / polish / selected objective、entry 相对 baseline / previous-best 的 objective delta、entry half-width、entry center / offset / center-source、iteration / exitflag、selected variant 与是否优于 baseline；该表会保留 SS/MS 动态方法的入口记录，便于横向核对 SS 侧 estimator-side `wide acquisition -> compact polish` 逻辑是否真实进入。`ms-crb-audit` 只对 `MS-MF-*` 显式传入 trim-focused truth-free external center；当前默认切到 `ssmf-axis-diagnostic`，在 `ssmf-center` 之外只加入小范围 SS-MF-centered axis-cross candidate，并把 estimator 侧 `doaBasinEntryAdoptionMode` 设为 `diagnostic-only`，因此这些新增 candidate 只进入 `msDoaBasinEntryTable` / joint-health / shadow audit，不改变 estimator final winner。SS-MF 不接收这些 external centers，继续走默认 single-sat wide acquisition，防止 SS 护栏被 MS 路线污染。`static-center / static-axis / static-diagonal` 仍仅能作为手动 `full-family` ablation 的历史对照；`static-axis` 已在 compact-gated ablation 中再次显示会拉高 trim-normalized angle MSE/CRB，不进入默认路线。`msEstimatorEntryFamilyAggregateTable` 使用 basin-entry 阶段的 `bestTag / selectedVariant` 归因，而不是 unknown warm-anchor 之后的最终 `solveVariant`，并输出 selected / selected-trim 的 angle 与 fdRef normalized MSE，用于判断 trimmed RMSE/CRB gap 是否仍由可保留 family 主导。`msEstimatorSelectedVsTruthResidualTable` 默认只保留 `trimKeep=true` 的 residual seed，把已被 angle-trim 剔除的 tail 留给 tail table，避免继续围绕可剔除 seed 过度调参；`msEstimatorResidualAggregateTable` 进一步按 residual type 汇总 trim-kept 样本的 angle/fdRef normalized MSE、truth-better rate 与 selected family 来源，用于判断下一步是搜索 gap、finite-sample residual 还是 CRB/FIM 口径排查。`msBestEntryOracleTable / AggregateTable` 使用 estimator 已经评估过的 truth-free basin-entry 候选，分别统计 objective-best candidate 与离线 angle-oracle upper bound 的 angle/fdRef MSE/CRB；它不运行新 solver、不用 truth 选择 runtime winner，只用于判断“已有候选族是否已经足够贴近 CRB”还是“候选族本身打不到 CRB-local”。`msJointHealthOracleTable / AggregateTable` 在同一批已评估候选上进一步要求 angle、`fdRef` 与 CP-U `fdRate` 同时健康，输出 strict target、loose joint-health、angle-good/fd-bad 三类 reachability；aggregate 同时报告 coverage-aware 口径，包括 eligible 子集、missing 子集以及“有 joint candidate 用 candidate、没有则回退 final”的 fallback normalized MSE，避免把 NaN-omit oracle 误读为全样本上限。`msAngleGoodFdBadTable` 专门列出 angle 已接近目标但 Doppler 状态不可采纳的候选，用来区分 candidate-family coverage、objective trade-off 与参数读取一致性问题。`msEntryShadowAdoptedCaseTable / DetailTable / OutcomeTable` 进一步构造 replay-only 的 truth-free entry-shadow adoption 诊断行：默认同时输出 strict 与 relaxed 两套 selector。strict 只采纳 objective 明显优于 final、来源为 `ssmf-center` 且相对 final 的 `fdRef` tooth step / `fdRate` step 均健康的候选；relaxed 允许 `ssmf-center/ssmf-axis` 候选在 objective 不明显劣于 final 的小容差内参与比较，并使用稍宽的 fd / rate step guard，用于判断“好 joint-health 候选只是 objective 略差”还是“selector proxy 本身无效”。两个 selector 都不读取 truth angle、truth `fdRef/fdRate` 或人工 tail label；truth 只在候选选定后用于评价 shadow row 的 RMSE/CRB、damage 与 target/loose pass rate。`DetailTable` 同时拆分 `source / objective / fdRef-step / fdRate-step` reject reason 和逐级候选计数，避免再把所有未采纳样本混成 `guard-rejected`。`msEntryShadowConfusionTable / AggregateTable` 将 strict / relaxed shadow 选择结果与 `msJointHealthOracleTable` 对齐，离线统计 `true-positive-joint`、`false-positive-damage`、`false-positive-no-joint`、`false-negative-joint-exists` 与 `true-negative-no-joint`，用于判断 selector 是误伤、漏救还是候选族缺失；该 confusion 只用于 replay 后评价，truth 不进入 shadow selector。这些表只属于 replay-level offline audit，不新增 solver，不改变 estimator winner、CRB 主公式、objective、search bounds 或默认 adoption。`ssDoaBasinEntrySummaryTable` 额外按 SS-MF method/SNR 汇总 baseline entry 是否启用、是否选中 entry/polish 以及 selected angle 中位数，作为 MS replay 中的 SS 护栏。`ssOnlyDebugMode=true` 时，compact 输出额外打印 SS solve-probe baseline/wide parity 与 SS entry trace，便于快速看到 `doaBasinEntry0.012 / doaBasinEntry0.048Polish` 是否重新出现。只有 MS 行会并入 `msCandidateTraceTable`。
- 搜索范围与 rescue 统计口径：`msSeedOffsetTable` 统计 static seed、SS-MF seed、baseline final 与 best solve 到 truth 的 lat/lon offset、球面距离和 fd offset；`msFdBranchProbeTable` 只对 MS-MF-CP-K 的 fdRef branch suspect 输出 final / truth mixed objective、fdRef over-tooth、first-order-opt 与 non-ref coherence，避免把频率 branch 和 DoA basin 混在一起解释；`msSearchRangeHintTable` 进一步给出 static-centered box 是否覆盖 truth、SS-MF center 是否更接近、best solve 与 seed 的距离、best static-wide candidate及 replay-only 推荐中心 / half-width；`msBankCoverageAggregateTable` 聚合 replay-only MS bank 中 CP-K preselect 后 top-K CP-U release 的收益；`msBankFamilyAggregateTable` 按 center source、axis direction、step 与 primary diagnosis 拆分 CP-U release 候选收益，`msBankCpkPreselectAggregateTable` 用已有 CP-K preselect rows 观察 cheap preselect 是否已经找对 DoA basin；`msBankAdoptionShadowTable / msBankAdoptionShadowAggregateTable / msBankAdoptionShadowDiagnosisAggregateTable` 用 objective / first-order-opt / non-ref coherence 这类 no-truth guard 离线评价“若采用 bank winner 是否会改善或误伤”，其中 diagnosis aggregate 只做分型统计；`msBankAdoptionRejectAggregateTable` 汇总 `no-candidate / guard-rejected / adopt` 等原因。`msBankRescueCaseTable / msBankRescueAggregateTable / msBankRescueDetailTable` 构造更直观的 `MS-MF-CP-U-BankRescue` replay 诊断行：若 shadow guard 采纳当前 bank candidate，就用该 candidate 的已计算 angle / fd / objective 指标替换 baseline，否则保持 baseline。`msBankBadGatedRescueCaseTable / Aggregate / Detail / Outcome / Tail` 构造 oracle 上限性质的 `MS-MF-CP-U-BadGatedBankRescue` 诊断行，只在 baseline 的 CRB-normalized angle 或 `fdRef` error 超过对应 bad-gate 阈值时允许采用，用于检查 fixed-SNR 下“只救 tail”能否让 trim / resolved 口径贴近 CRB。`msBankHealthGatedRescueCaseTable / Aggregate / TrimAggregate / Detail / Outcome / Tail` 进一步构造 `MS-MF-CP-U-HealthGatedBankRescue` 诊断行，只用运行时可观察的 solver / coherence / boundary / no-solve 健康信号，以及 shadow candidate 的 objective margin 触发 replay-only bank rescue；其中 first-order-opt 与 iteration 默认需要同时出现，candidate-objective trigger 只使用已计算 objective 差，不读取 truth / CRB error。它仍不改变 estimator 默认 winner。`msBankHealthGatedTriggerOutcomeTable` 按触发源统计 health gate 的 adopt、trim keep、soft damage 与收益，用于判断是否应继续收紧 first-order / iteration / candidate-objective 类触发。`msBankHealthGatedResidualPriorityTable / Aggregate` 只对 health-gated rescue 后仍未贴近 CRB 的 residual seed 做离线优先级排序，区分 `rescue-angle-tail`、`rescue-angle-residual`、`rescue-fdref-tail`、`guard-rejected` 和 `no-candidate`，并给出是否应查 fdRef branch、no-truth guard 或转向 estimator-side selected-vs-truth residual 的下一步提示；这些字段使用 CRB-normalized error 仅做 replay 后处理，不参与 runtime selector。`msBankRescueOutcomeTable`、`msBankBadGatedRescueOutcomeTable` 与 `msBankHealthGatedRescueOutcomeTable` 按 fixed-SNR repeat 批量统计 baseline bad rate、oracle bad-gate case rate、health-gate case rate、candidate coverage、adopt rate、healthy/bad adopt rate、soft damage rate、large improve rate、same-tooth adopt rate、rescue bad/trim keep rate 与 rescue 前后中位数；对应 tail table 列出仍未解决、发生 adoption 或 gate 漏救的 repeat，便于下一轮按 seed 批量发现问题。相关表只服务搜索范围、candidate-family、trim/resolved 口径和 replay-only bank-rescue 设计，不把 truth 覆盖判断放入 runtime selector。
- 统计 / CRB-floor 口径：`trimAggregateTable`、`msBank*TrimAggregateTable` 当前使用 angle-tail trim，不再用 `fdRef` error 决定 `fdRef` MSE，避免 trimmed `fdRef` MSE 被条件筛选压到 CRB 以下；`jointTrimAggregateTable`、`msBank*JointTrimAggregateTable` 保留 angle+fdRef cap，仅作为诊断对照。`crbFloorCheckTable` 标出 MSE/CRB 低于 1 的 baseline 行，`objectiveTruthAuditTable / Aggregate` 和 `estimatorBasinObjectiveAuditTable` 将真值代价、final 代价、estimator-side DoA basin-entry 候选代价放在同一口径下比较，方便判断搜索范围不足、候选未构造还是候选未赢；新增 `crbNormalizedAuditTable` 按 raw / health / angle-trim / joint-trim 直接报告 `mean((error/CRB)^2)`、1.1 目标是否通过，以及旧 `MSE/median(CRB)^2` 口径的对照，避免把 median-CRB 聚合误读成逐样本 CRB-normalized MSE；`crbMethodPairAuditTable` 按 SS/MS method pair 对比 CRB median、normalized MSE 和 MS/SS CRB 比值，`crbComponentAuditTable` 汇总 angle CRB 的 lat/lon 分量、trace/spherical ratio、相对 SS-MF 的 tightness、interest/full FIM condition、unit-scaled FIM condition 与 DoA-fdRef FIM coupling，用于判断 trim 后 gap 是 estimator residual、有限样本 objective variance，还是 CRB/FIM 口径异常；`crbVariantAuditTable` 进一步用 current MS CRB、SS-MF 等价 CRB、fdRef-known DoA block CRB 和 nuisance-known DoA block CRB 重新归一化同一批误差，避免把“MS CRB 分母过紧”误判成 estimator 唯一问题；`msTrimResidualCrbVariantTable` 只列 `trimKeep=true` 的 MS residual seed，比较这些 seed 在不同 CRB 分母下是否仍远离 CRB，不让 angle-trim 可剔除 seed 牵引候选调参；`fdRefNormTail` 会把 angle 合格但 `fdRef` 超 cap 的样本留在 tail 表中，便于专门回查 branch / summary bias。

- 代表样本口径：`msRepresentativeCaseTable` 按 primary diagnosis 自动挑选每类最多 `representativeMaxCasePerType` 个高 severity 样本，并附带下一轮 replay hint；它只用于后续手动 full deep-probe / flow-like replay 的 seed 列表，不在本脚本内二次 rerun，也不引入 truth-aware runtime selector。
- MS bank 口径：固定为 replay controlled probe 层的 `MS-MF-CP-U` axis-cross compact DoA candidates；全局默认 `msBankAxisStepDegList=[0.006;0.012;0.024]`，但固定 SNR 的 `ms-fast / ms-fixed-snr / ms-sanity` 会收缩为 `[0.006;0.012]`；0.024/top-2 residual-tail 路线已删除，`0.048 deg` 只在 `snrDb<=msBankWideStepMaxSnrDb` 时追加，避免高 SNR blanket-wide damage；先用 CP-K cheap preselect 排序，再只对 top `msBankReleaseTopK` 做 CP-U release。该 bank 仍使用同一个 truth-centered in-tooth `fdRef/fdRate` 诊断盒，不放宽 fd 范围，不新增非零 fdRef / fdRate offset bank，不进入 estimator 默认路径，不改变 objective / residual / reference-sat 语义；若运行时间仍过长，优先调小 axis step / top-K，而不是恢复旧路线开关。
- 运行时间与稳定报告表口径：`releaseCompareTable`、`seedChainTable`、`perSatProbeTable`、`runtimeTable / runtimeAggregateTable / topSlowRuntimeTable` 的稳定构表 helper 已抽到 `test/common/report/`；它们只收集或比较已有 case / probe rows，不定义 MS 诊断标签，也不触发新 estimator 调用。runtime 表记录 repeat-data、static bundle、CRB、baseline solve、objective probe、solve probe、per-sat probe 与 method-total 的 wall-clock 时间，用于判断瓶颈是在 fixture/CRB、MS-MF-CP-U baseline、solve probe 还是 objective probe，不据此改变 estimator 默认并行策略。
- 打印口径：默认 `summaryPrintLevel="compact"`，命令行只打印 raw / trim aggregate、MS primary diagnosis、wide coverage、bank-rescue / bad-gated / health-gated raw 与 trim aggregate、runtime top 与代表样本等短表；`diagnostic` 才打印 seed chain、objective / solve / candidate trace、search range 等长表。完整表始终保存在 `replayData`。
- 存储与并行口径：`runScale="ss-sanity"` 与 `"ms-sanity"` 默认关闭 per-task checkpoint / resume 和 snapshot，但保留 Telegram 通知与自动外层 task-level `parfor`，用于在低成本 wiring / 表结构检查中仍能及时收到完成/失败消息，并在 MS sanity 的少量独立 task 上利用 CPU；`runScale="full"` 时启用原代表性 full replay 的 per-task checkpoint / resume，任务文件写入仓库根目录 `tmp/replayMfMsMleCrbFlowDiagnose/<runKey>/task/`。外层 task loop 参考 `scanMfMleCrbInToothConsistency` 自动检测 Parallel Computing Toolbox，checkpoint 与非 checkpoint 路径均优先使用 task-level `parfor`，不可用或嵌套 worker 中自动退回串行；`runKey` 使用短语义前缀加 8 位 hash，且包含 `runScale` / route tag，避免不同规模或 wiring 版本复用旧 checkpoint；成功构造并保存轻量 `replayData` 后清理 checkpoint。
- 当前定位：作为 MS 估计流程剖面入口；若 scan 里 MS 不贴 CRB，应先从 `topTailExportTable` 选 seed，再用本 replay 看 seed chain、per-sat coherence/local-state、objective probe、bank-rescue aggregate 和 reject aggregate，不要把 scan 膨胀成内部 trace 脚本。`MS-MF-CP-U-BankRescue` 只回答“已有 replay-only bank candidate 若按 no-truth shadow guard 被采用，统计会怎样”；它仍是 replay-level 诊断行，不代表 estimator 默认 winner 已改变。

### 5. tooth 已对后，检查 same-tooth DoA 坏盆地

#### `replayMfInToothFdRangeOracle.m`

- 作用：把参考星 `fdRef` 搜索盒固定到 truth-centered 半齿范围，并用同一批 seed 对比 `SS-SF-Static / MS-SF-Static / SS-MF-CP-U / MS-MF-CP-U / MS-MF-CP-K`。
- 选帧口径：该 replay 已经固定 tooth，不再做 subset tooth selection；repeat fixture 只构造 periodic all-frame view，跳过 curated / random subset fixture bank。
- 先看现象：当 wrong-tooth 被 oracle 排除后，`MS-MF-CP-U-in-tooth` 是否优于单星、单帧与 wide baseline；若它仍不优于这些上限基线，说明论文主张上限不足，应先回查模型层级或 same-tooth refine，而不是继续扩大 subset bank。
- 主要输出：method aggregate、MS-MF wide-vs-oracle pair compare、paper-claim upper-bound compare、runtime timing summary、oracle fd/fdRate range preview、`singleMultiCompareTable`、`tailCaseTable` 与两张分布 / tail 诊断图。
- 表格口径：`aggregateTable` 额外保留 `satMode / frameMode / modelClass / rateMode / oracleLevel / wallTimeMedianMs`，`paperCompareTable` 以 `MS-MF-CP-U-in-tooth` 为 target，正 gain 表示 target 优于对应 baseline；`singleMultiCompareTable` 逐 seed 比较 `SS-MF-CP-U-in-tooth` 与 `MS-MF-CP-U-in-tooth`；`timingTable / timingAggregateTable` 只记录 repeat 级轻量 wall-clock 计时，用于判断 static bundle 与 dynamic method 的耗时占比。
- 命令行口径：`rangeTable` 仍完整保存在 `replayData`，但命令行通过 `dispMfReplayTablePreview` 只打印前四行和后四行，避免 repeat 较多时长表刷屏。
- 临时存档：默认 `checkpointEnable=true`，按 repeat 在仓库根目录 `tmp/replayMfInToothFdRangeOracle/<runKey>/task/` 保存轻量 checkpoint；成功保存 `replayData` 后默认清理，失败时保留并打印路径。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存与 checkpoint 清理后发送 HTML `DONE` 通知，在失败时发送 `FAILED` 通知；通知只包含脚本名、耗时、snapshot / checkpoint、MS-MF better angle rate、median angle gain 和 tail case 数量。
- 图形口径：默认不再画 method-level RMSE bar，也不再叠加多个方法的直方图；改为 angle / fdRef / fdRate 经验 CDF、单独的 MS-vs-SS angle gain 直方图与 paired SS-MF vs MS-MF seed scatter。
- 当前结论：只作为论文主张上限 replay，不进入默认 flow / regression，不改变 no-truth-leak selector。

#### `replayMfInToothTailCaseDiagnose.m`

- 作用：固定 `replayMfInToothFdRangeOracle` 暴露出的 tail seed，默认回放 coherence-collapse tail seed `277 / 283 / 298 / 256` 与负样本 `293 / 280 / 268 / 284`，重放 `SS-MF-in-tooth / MS-MF-in-tooth / MS-MF-truth-DoA-oracle` 等核心 case，并分类 same-tooth tail。
- 选帧与初始化口径：该 replay 已经在 in-tooth oracle 范围内定位 tail，不做 subset tooth selection，也不使用 curated subset 初始化；DoA seed 只来自 static seed 或 truth-DoA oracle；`contextBaseSeed` 固定为源 oracle replay 的 base seed，避免因为只改 seed 列表顺序而换掉共享场景上下文。
- 主要输出：`identityTable`、`aggregateTable`、`tailDiagnosisTable`、`detailDiagnosticTable`、`candidateProbeTable`、`candidateWinnerTable`、`lineProbeSummaryTable`、`rescueBankSummaryTable`、`rescueBankAggregateTable`、oracle range / runtime summary，以及 gated rescue 前后 angle / coherence 对比图。
- 临时存档：默认 `checkpointEnable=true`，按 tail seed 在仓库根目录 `tmp/replayMfInToothTailCaseDiagnose/<runKey>/task/` 保存轻量 checkpoint；成功保存 `replayData` 后默认清理，失败时保留并打印路径。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存与 checkpoint 清理后发送 HTML `DONE` 通知，在失败时发送 `FAILED` 通知；通知只包含脚本名、耗时、snapshot / checkpoint、tail seed 数量和 rescue aggregate 行数。
- 对照口径：除 static / truth DoA 可释放分支外，额外包含 `static-doa-fixed` 与 `truth-doa-fixed` 对照，用于判断只释放 `fdRef/fdRate` 时 non-ref coherence 是否可恢复。
- 诊断口径：`identityTable` 用于核对同一 seed 的 truth DoA、static DoA、oracle fd/fdRate 范围和 method count 是否与源 oracle replay一致；`detailDiagnosticTable` 优先展开非 light tail seed，若本轮没有非 light seed，则展开全部 seed 的 init/final per-sat coherence、residual 与差分 Doppler / rate 误差，方便判断 collapse 是初值已坏、优化后被拉坏，还是坏点已不复现；`candidateProbeTable` 只重评估 default / fixed-DoA / wide / truth-DoA、final 附近的小网格、wide-final 附近的小网格、default/static 到 truth DoA 的加密 line probe，以及 default final/static MS/SS-MF final/wide final 这些不含 truth 的 implementable center 粗 DoA 网格，不改变 solver adoption；`candidateWinnerTable` 总结每个 seed 是否可由 final-centered、wide-centered、implementable-center 或 line-to-truth 候选救回；`lineProbeSummaryTable` 记录 line probe 中 objective / coherence / angle 最早在什么 alpha 恢复，用于判断好 basin 离 default/static 有多远；`rescueBankSummaryTable` 与 `rescueBankAggregateTable` 比较 disabled、wide-centered coarse、single-MF-centered coarse 与 wide+single-MF bank，专门观察 hard seed 的救回率以及 easy / fd-not-healthy 负样本是否被误伤；默认图形只画 default MS-MF 与 `gated-wide-single-bank` selected result 的前后对比，不再用连线连接不同 seed。
- 分类口径：区分 `wrong-tooth`、`same-tooth + fd not healthy`、`same-tooth + fd healthy + non-ref coherence collapsed`、`same-tooth + fd healthy + DoA/local-state basin` 与轻微/不明确 tail。
- 当前结论：只用于 tail 定位和 gated refine 前置分类；候选 objective probe 用于决定下一步是 conditional DoA polish、wide-centered gated refine、single-MF-centered basin-entry，还是更宽的 basin 进入机制；rescue bank 只模拟候选选择，不改变 solver adoption；line probe 含 truth 只作 replay 机制定位，不进入 regression，不改变默认 flow。

#### `replayMfInToothDoaDopplerRidgeTrace.m`

- 作用：固定 in-tooth tail / negative seed，沿 final→static-MS、final→SS-MF、final→wide 与 final→truth 四类方向重评估 DoA-fdRef / DoA-fdRate 局部 surface，并做 center-to-center line trace。
- 选帧与初始化口径：继续使用 truth-centered half-tooth `fdRef` 范围，只做 objective 诊断，不做 subset tooth selection，不改变 solver adoption，也不把 truth 信息带入可实现方向。
- 主要输出：`ridgeSummaryTable`、`centerLineSummaryTable`、`candidateProbeTable`、`tailDiagnosisTable`、`aggregateTable`、oracle range / runtime summary，以及代表性 DoA-fdRef / DoA-fdRate ridge heatmap。
- 诊断口径：`static-ms / single-mf / wide` 是 truth-free implementable 方向，`truth` 只作 oracle 标记；重点看 ridge slope、minimum 是否贴边界、objective gain、non-ref coherence 是否恢复，以及 easy / fd-negative seed 是否也被 surface candidate 吸走。
- 当前结论：新增为 in-tooth DoA/local-state 坏盆地的 ridge 诊断入口；代表性结果显示 final-centered ridge 不适合直接公式化；当前只把 gated `wide+single-MF` basin-entry 作为 controlled in-tooth DoA rescue 候选，不把 ridge slope 或 flow 接入作为下一步。
- 详细结果：`results/replayMfInToothDoaDopplerRidgeTrace.md`。


#### `replayMfInToothDoaFdRangeEnvelope.m`

- 作用：用 `baseSeed / numSearchRepeat` 先跑轻量 scout 分类，再只对 hard / gate-miss / easy / fd-negative 代表 seed 扫描 DoA-fdRef 与 DoA-fdRate 局部 envelope，用于决定 controlled in-tooth rescue 的 DoA / fdRef 搜索半宽。
- 选帧与初始化口径：只做 replay-level objective surface 诊断，不做 subset tooth selection，不改变 solver adoption；`static-ms / single-mf / wide` 是可实现 center，`truth` 仍只作 oracle 方向。`fixed-list` 模式保留为手动复查入口。
- 主要输出：`selectedEnvelopeSeedTable`、`scoutTailDiagnosisTable`、`candidateCoverageTable`、`policyRecommendationTable`、`envelopeSummaryTable`、`rangeRecommendationTable`、`ridgeSummaryTable`、`centerLineSummaryTable`、`candidateProbeTable`、`tailDiagnosisTable`、oracle range / runtime summary。
- 临时存储：默认 `checkpointEnable=true`，scout 与 envelope 两个阶段分别在仓库根目录 `tmp/replayMfInToothDoaFdRangeEnvelope-scout/<runKey>/task/` 与 `tmp/replayMfInToothDoaFdRangeEnvelope-envelope/<runKey>/task/` 保存 per-seed checkpoint；成功保存 `replayData` 后默认清理，失败时保留并打印路径。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在两阶段 checkpoint summary 与 snapshot 保存后发送 HTML `DONE` 通知，在失败时发送 `FAILED` 通知；通知只包含脚本名、耗时、snapshot、scout/envelope checkpoint 路径和少量 seed / coverage / policy 指标。
- 诊断口径：重点看每个 center/surface 的 `minDoaOffsetDeg`、`minFdRefOffsetHz`、`boundaryHit`、`p95AbsDoaOffsetDeg` 与 `p95AbsFdOffsetHz`；同时用 `candidateCoverageTable` 判断 implementable center 是否覆盖 truth-DoA oracle 的 hit，用 `policyRecommendationTable` 离线比较 gate threshold、family DoA step 与 safe coherence 的组合。若 hard seed minimum 频繁贴边界，应先扩大 replay envelope，不把该范围直接写进 rescue bank。
- 当前结论：新增为搜索范围决策入口；默认两阶段运行，避免把全部 scout seed 都做重 surface。它不替代 `replayMfInToothGatedRescueEffectiveness` 的效果验证，只为后者是否需要 joint DoA-fdRef candidate bank 提供范围依据。
- 详细结果：`results/replayMfInToothDoaFdRangeEnvelope.md`。

#### `replayMfInToothGatedRescueEffectiveness.m`

- 作用：在更多 repeat 上验证 in-tooth 条件下的 no-truth gated rescue 是否稳定救回 same-tooth collapse，并确认 easy / fd-not-healthy 负样本是否被误伤。
- 选帧与初始化口径：使用 truth-centered half-tooth `fdRef` 范围，只验证 tooth 已正确时的 basin-entry rescue；不做 subset tooth selection，也不改变默认 flow。
- 候选口径：`caseRole / isHardRescued / isDamaged` 可用 truth 做离线评价；`rescueTriggered / triggerReason / selectedCandidateFamily` 只能由默认估计的 non-ref coherence、phase residual、candidate objective 和卫星几何回代诊断决定。当前 gate 采用 `coherence-v1`，只用 default non-ref coherence collapse 触发；phase residual 仅保留字段占位，不参与本 replay 的触发。
- 主要输出：`inToothDoaLadderTable`、`rescueEffectAggregateTable`、`rescueEffectVerdictTable`、`rescueEffectCaseTable`、`rescueBankDecisionTable`、`triggerReasonTable`、`hardMissSummaryTable`、`candidateProbeTable`、method / range / timing summary，以及按 task seed 画出的 default / gated / blanket / truth-DoA oracle angle / coherence 线图和 trigger / case-role 线图。
- 工程外壳：该脚本已按第一批 replayTemplate 化口径接入 common header、checkpoint option / summary、长表预览与 `notifyMfReplayStatus` 通知壳；run-key signature 与 gated-rescue metric lines 仍保留在本地，避免改变 checkpoint 兼容性和结果解释口径。
- 对照 bank：默认保留 `disabled`、严格 DoA safe-adopt 的 `gated-wide-single-bank-safe-adopt`，以及新增的 `gated-wide-single-bank-family-safe-adopt`。后者只放宽 wide / single-MF basin-entry center 的 DoA step guard，用于验证 seed 256 这类被统一阈值拒绝的好 candidate 是否可安全采纳。已验证无增益或不稳的 `wide-only / single-mf-only / unsafe gated / blanket` 默认旁路。如需复查旧路径，可显式打开 `includeLegacyFailedBanks`。
- joint bank 口径：DoA-fdRef joint candidate 仍保留为 opt-in 诊断分支，只有在文件头打开 `includeJointSafeAdoptBank` 时才生成 joint grid 并输出 `gated-wide-single-joint-bank-safe-adopt`，避免默认增加候选和表格噪声。
- safe adoption 口径：safe-adopt 只在 replay 层增加 truth-free adoption guard，要求 candidate 相对 default 有 objective gain、non-ref coherence 恢复、DoA step 和 fdRef step 均不超过头部阈值；family-safe 仅对 wide / single-MF center 使用更宽的 `familySafeAdoptMaxAbsDoaStepDeg`，default/final-centered candidate 仍用严格阈值。不改变 estimator、flow 或 objective。
- 当前结论：只作为 gated rescue 批量验证 replay；当前默认代码收敛为 coherence-only gate + DoA-only strict/family safe-adopt 对照，用于区分 gate 漏检、candidate bank 不足、统一 DoA-step guard 过紧与 triggered-selected damage。joint offset 先交给 envelope replay 判断，默认不进入效果验证主表。当前短期目标收回到 controlled in-tooth DoA 修复，不把该 replay 结论直接推进 flow 或默认 estimator。
- 详细结果：`results/replayMfInToothGatedRescueEffectiveness.md`。

#### `replayMfInToothOrdinaryAngleMissDiagnose.m`

- 作用：在 controlled in-tooth 条件下专门诊断 non-collapse ordinary DoA miss，区分普通 angle miss、collapse hard、mid-coherence miss、wrong-tooth 与 fd-not-healthy 负样本。
- 选帧与初始化口径：沿用 truth-centered half-tooth `fdRef` 范围，不做 subset tooth selection，不改变默认 flow；truth 只用于 offline miss label 与 truth-DoA oracle 上限。
- 候选口径：比较 default final、final-centered small DoA polish、wide-centered DoA-only basin-entry、single-MF-centered DoA-only basin-entry、best implementable 与 truth-DoA oracle；final-centered polish 只作为 replay probe，不进入 estimator / flow 默认路径。
- 主要输出：`ordinarySeedTable`、`ordinaryCandidateTable`、`ordinaryAggregateTable`、`ordinaryFamilyAggregateTable`、原有 method / gated rescue 辅助表、runtime timing summary，以及 ordinary miss 角误差 / gain / coherence 对比图。
- 工程外壳：该脚本已按第一批 replayTemplate 化口径接入 common header、checkpoint option / summary、长表预览与 `notifyMfReplayStatus` 通知壳；run-key signature 与 ordinary metric lines 仍保留在本地，避免改变 checkpoint 兼容性和结果解释口径。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存后发送 HTML `DONE` 通知，在 batch / storage 失败时发送 `FAILED` 通知并继续 `rethrow`；通知只包含脚本名、耗时、seed、snapshot 与少量 ordinary summary，发送失败只 warning。
- 当前结论：新增为 controlled in-tooth DoA 主线的 ordinary-miss 分叉诊断入口；用于判断剩余 `>0.002 deg` miss 是否可由 small polish 或 basin-entry truth-free 候选解释，不据此直接放宽 hard-collapse rescue gate。

#### `replayMfInToothOrdinaryWideGateDiagnose.m`

- 作用：在 controlled in-tooth 条件下专门诊断 ordinary miss 的 `wide` basin-entry 候选能否被 truth-free gate / adoption 接住。
- 选帧与初始化口径：沿用 truth-centered half-tooth `fdRef` 范围，不做 subset tooth selection，不改变默认 flow；truth 只用于 offline miss label、rescue / damage 评价和 truth-DoA oracle 上限。
- gate 口径：policy sweep 只使用 default / wide candidate 的 objective gain、non-ref coherence、DoA step 上下界、fdRef / fdRate step 与 default non-collapse 状态；不读取 truth tooth、truth DoA、truth `fdRef/fdRate` 或人工 case label。当前默认细扫 objective gain `0~200`，并增加 default-vs-wide DoA-disagreement 下限，用于观察能否压低 easy trigger 同时保留 ordinary rescue。
- 推荐 policy 口径：当前 replay-level 首选 ordinary-wide gate 固定为 `objGain >= 10`、`0.001 <= |wide DoA step| <= 0.004 deg`、wide-only、default non-collapse、candidate coherence `>=0.95`。脚本会额外生成 `wideGateRecommendedPolicyTable`，按 full ordinary rescue、zero damage、低 easy trigger、低整体 tail 的顺序自动选择最稳候选，避免命令行和 Telegram 误选第一个过宽 pass policy。
- 旁路口径：该脚本已按第一批 replayTemplate 化口径整理为 common header、checkpoint option / summary、snapshot、长表预览与 `notifyMfReplayStatus` 通知壳；默认旁路已证伪或非本入口主线的 `final-small-polish` probe、`single-MF` ordinary basin-entry probe 与 hard-collapse auxiliary tables。若要复核旧负结果，可在头部显式打开 `includeFinalSmallPolishProbe`、`includeSingleMfBasinEntryProbe` 或 `includeHardCollapseAuxTables`。run-key signature 与 wide-gate metric lines 仍保留在本地，避免改变 checkpoint 兼容性和结果解释口径。
- 主要输出：继承 ordinary replay 的 `ordinarySeedTable`、`ordinaryCandidateTable`、`ordinaryAggregateTable`，并新增 `wideGateFeatureTable`、`wideGateDecisionTable`、`wideGatePolicyTable`、`wideGateRecommendedPolicyTable`，用于比较不同 objective-gain / DoA-step lower/upper guard 下的 ordinary rescue、easy trigger / damage、fd-negative damage 和 overall tail。默认旁路的 legacy probe 不再出现在 ordinary family aggregate 中。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存后通过 `notifyMfReplayStatus` 发送 HTML `DONE` 通知，在 batch / storage 失败时发送 `FAILED` 通知并继续 `rethrow`；通知只包含脚本名、耗时、seed、snapshot、ordinary summary、自动选择的 wide-gate policy 摘要及 easy trigger/damage 口径，发送失败只 warning。
- 当前结论：作为 `replayMfInToothOrdinaryAngleMissDiagnose` 的后续 gate/adoption 诊断入口；它只判断 `wide` candidate 是否存在可用 truth-free proxy。即使某个 policy 标为 `preferred-candidate` 或 `candidate`，也只是 replay-level gate 候选，不据此直接打开 blanket wide 或修改 estimator 主核。

#### `replayMfFlowLikeGatedBasinEntryEffectiveness.m`

- 作用：在真实 subset-periodic flow 结构内比较 disabled 与 no-truth gated `wide + single-MF` same-tooth basin-entry rescue。
- 选帧与初始化口径：不使用 truth-centered 半齿 oracle，不跳过 subset selection；只通过 `buildSimpleDynamicFlowOpt.sameToothRescue` 显式打开 flow-level rescue，默认 flow 仍保持关闭。
- gate 口径：只读取 selected periodic decision summary 与 subset trust 诊断中的 non-ref coherence、phase residual 和 trust 状态；truth 只用于结果表里的 angle / tooth / fd 评价。
- 主要输出：disabled-vs-gated compare table、same-tooth in-tooth compare / aggregate、gate-miss hard-case table、case-role summary、gate-reason count、final-tag / tooth transition count、gated periodic candidate table、warning table、storage summary，以及 angle / coherence 对比图。
- 存储口径：snapshot 仍只保存 `replayData` 一个变量；`saveRepeatDetail=false` 时只保存轻量 table 与 checkpoint summary，适合扩大 repeat 后通过 `loadExpSnapshot` 重跑 summary 小节；若需要追单 seed 的完整 flow 结构，才临时打开 `saveRepeatDetail=true`。
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存后发送 HTML `DONE` 通知，在 batch / storage 失败时发送 `FAILED` 通知并继续 `rethrow`；通知只包含脚本名、耗时、seed、snapshot / checkpoint 和少量 aggregate / gate-miss 指标，发送失败只 warning。
- 临时存档：头部 `checkpointEnable=true` 时分别为 disabled 与 gated batch 建立 per-repeat checkpoint / resume，任务文件写入仓库根目录 `tmp/replayMfFlowLikeGatedBasinEntryEffectiveness-disabled/<runKey>/task/` 与 `tmp/replayMfFlowLikeGatedBasinEntryEffectiveness-gatedWideSingle/<runKey>/task/`；成功构造 `replayData` 并清理 checkpoint 后再保存 snapshot。
- 当前结论：该入口只作为真实 flow 条件下的 stress check；由于会混入 wrong-tooth / subset selection 因素，当前 in-tooth DoA 主线不继续围绕它调参，也不据此打开默认 flow。

#### `replayMfSameToothHardCase.m`

- 作用：追踪 `fdRef/fdRate` 健康且 same-tooth 已成立，但 DoA 仍停在坏盆地的样本。
- 先看现象：`fdHealthyRate` 与 `sameToothHardRate` 是否同时为非零。
- 主要输出：hard-case repeat、representative candidate 表、fd/DoA 对比。
- 当前结论：fd 健康不代表 DoA 已进入好盆地，same-tooth refine 仍需要单独观察。

#### `replayMfInToothPolishGate.m`

- 作用：搜索 polish-eligible seed，并比较 conditional polish 与 disabled polish。
- 先看现象：polish 是否只在 eligible hard case 触发，是否避免 blanket polish。
- 主要输出：eligible seed、trigger reason、angle improvement、误触发情况。
- 当前结论：没有 eligible 样本时只能说明未误触发，不能证明 polish 无效。

### 6. unknown-rate release 路线说明

#### `replayMfUnknownReleaseRoute.m`

- 作用：比较 CP-K seed、CP-U fd-only release、CP-U DoA+fdRate release。
- 先看现象：只 release fdRate 时 DoA 是否基本不动；放开很小 DoA 盒后，objective / angle 是否改善。
- 主要输出：route summary、DoA box sweep、DoA bound margin。
- 当前结论：unknown-rate release 要区分 fdRate release 与 DoA release，不应只看 final objective。

#### `replayMfWarmAnchorParforSensitivity.m`

- 作用：比较串行 warm-anchor 与显式内部 warm-anchor parfor 对 winner / tooth 的敏感性。
- 先看现象：内部 parfor 是否改变 winner、tooth 或 final tag。
- 主要输出：serial/parfor winner 对比、tooth 变化、angle/fd pair plot。
- 当前结论：warm-anchor 内层 parfor 只能显式 opt-in，不能进入 estimator 默认路径。

### 7. 最后检查 compact flow 与 representative/fullDiag 是否分叉

#### `replayMfFastStatsRepresentativeDivergence.m`

- 作用：比较 fast compact batch 与 representative/full-diagnostic rerun。
- 先看现象：compact flow 是否遗漏会改变 final winner 的候选，是否存在 adoption 差异。
- 主要输出：candidate coverage、winner adoption、selectedFinalTag 与 tooth 差异。
- 当前结论：若 compact 与 representative 分叉，优先收紧 flow candidate coverage，不膨胀 estimator branch。

## results 文档与 snapshot 绑定

- replay 的详细运行结果放 `test/dev/replay/results/<scriptName>.md`。
- 大 `.mat` snapshot 放 `test/data/cache/replay/`。
- 结果文档记录 snapshot 文件名、关键配置、统计表、观察现象和当前结论。
- 若某个 replay 的结果很多，可将结果文档升级为 `results/<scriptName>/README.md` 与多个子结果文档。
- 新结果先更新对应 results 文档；排障主记录只摘取影响当前优先级的结论。

## 存储与画图规范

- 默认画图，不再提供绘图开关。
- 多阶段或多策略曲线图必须给每个 subplot 配清楚图例；同一脚本内的 stage 命名要一致，避免图例和表格字段各叫一套。
- 不保存图片文件，只保存可重画图的数据，例如 line/surface scan 数据、summary table、candidate table、representative case。
- tmp 只作为运行时临时目录，且统一位于仓库根目录 `tmp/<scriptName>/<runKey>/`；不能写到 `test/tmp`、当前工作目录或脚本目录。正常完成后必须清理，失败时由 `catch` 打印现场路径并保留。
- 小 MC replay 若单个 repeat 较慢，可以在文件头保留一个 `checkpointEnable` 开关。开启时按 repeat 做 checkpoint / resume：每个 repeat 一个轻量 task 文件，manifest 记录 seed、SNR、contextOpt 和影响分支行为的 flow signature；关闭时不创建 tmp，仍直接运行完整 repeat batch。恢复时只重跑缺失 task。
- checkpoint 只保存 repeat 结果，不保存 `rxSigCell`、完整 `sceneSeq`、fixture cell、transition bundle 或全量 objective map；成功构造 `replayData` 后默认调用 `cleanupPerfTaskGridCheckpoint` 清理 checkpoint run 目录；如果脚本级 tmp 父目录已经为空，也由该 common helper 同步删除该空目录。replay 脚本不得为 checkpoint run 目录再写私有 cleanup helper。
- 不需要中间落盘的固定单样本 replay 不创建 tmp 目录，也不打印 temporary run dir disabled 这类无信息输出；只保留真正影响运行和结果解释的配置。
- snapshot 默认只保存 `replayData`，保存路径为 `test/data/cache/replay/`；不要保存 `rxSigCell`、`sceneSeq`、fixture cell、transition bundle 或全量 objective map。

## 恢复后重出结果

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开对应 replay 文件，直接运行 `Summary output and plotting` 小节即可重新输出表格和图。

## 头部开关收敛规则

- 不再写 `saveMaxVarBytes`。`saveExpSnapshot` 已经有 `maxVarBytes` 默认值；当前 replay 只通过 `includeVars={'replayData'}` 保存轻量结果，不需要在每个脚本重复定义。
- 不再写并行开关。小 Monte Carlo 默认优先走外层 `parfor`，并行工具箱不可用或 repeat 数不足时自动退回串行。
- 不再写进度条开关。固定样本 surface/grid scan 和小 MC repeat batch 默认显示 progressbar；找不到 progressbar 时只打印一次紧凑提示并继续运行。
- 不再写 context / parallel 的通用 override。需要改变 context 或内部并行策略时，只在具体脚本对应逻辑处显式修改，不做统一头部开关。
- 不再写 snapshot 目录或前缀。目录使用 `saveExpSnapshot` 的 replay task type 路径，前缀使用脚本名。
- figure 和 table 默认全部输出，不提供细碎开关。
- `tmp` 正常结束一定清理；失败时保留现场并打印路径。
- checkpoint / resume 只加在确实耗时的小 MC replay 上；固定单样本短 replay 不为了形式保留 checkpoint 外壳。

## MATLAB 实现细节规则

### `parfor` 广播变量

- 使用 `parfor` 前必须检查 MATLAB Code Analyzer 的广播变量提示。
- 对大数组、grid、offset、candidate list，不能在 `parfor` 内用派生下标间接访问，例如先用 `ind2sub` 得到 `iFd` 后再读 `fdOffsetGrid(iFd)`；这会让 MATLAB 把整个 grid 识别为 broadcast。
- 正确做法是在 `parfor` 前预先展开为按循环变量直接索引的 sliced vector / sliced matrix，例如 `fdRefEvalVec(iPoint)`、`doaLatEvalVec(iPoint)`、`candidateCell{iRepeat}`。
- 小标量、必要的只读 model、`parallel.pool.DataQueue` 这类确实必须广播的对象可以保留；除此之外，看到不必要 broadcast 就要立即清理。
- 不为消除警告改变 objective、candidate 顺序、随机数路径或 estimator 默认并行语义。

### 时间戳写法

- 新代码不要再用 `datestr(now, ...)`。
- 日志时间统一用 `datetime('now', 'Format', ...)`，需要传给 `fprintf` 时再转 `char(...)`。
- snapshot 文件名仍由 `saveExpSnapshot` 统一处理，脚本内不要自造 timestamp 前缀。

### 冗余变量、分支和 helper 零容忍

- 固定单样本 replay 不保留无效 `numRepeat`、空 `runDir`、空 tmp 状态字段、无用 `runKey`、temporary run dir disabled 打印或只为包装而存在的 `replayConfig`。
- 没有中间落盘需求时，不创建 tmp，不写 cleanup 分支，也不维护失败现场路径。
- 没有被实际使用的 helper、字段、summary 表、plotData wrapper、兼容分支应直接删除。
- 一次性短 glue 直接写在执行流程中；objective scan、surface scan、summary、plot、progress 这类块状逻辑才保留 local helper。

### local helper 注释规则

- replay 脚本中的每个 local helper 至少在函数行后保留一句英文注释，说明它负责什么。
- 注释只解释职责边界和关键语义，不重复逐行代码；短 glue 也要说明为什么存在。
- 新增或清理 helper 时同步检查是否真的需要拆分：一次性短 glue 直接写在执行流程中，objective scan、surface scan、summary、plot、progress 这类块状逻辑才保留 local helper。

### 小 MC 统计输出规则

- 小 MC replay 如果用于比较多个阶段或策略，除了逐 seed 表格，还应提供能体现分布收敛的统计视图。
- tooth selection 类 replay 默认保存 `|toothIdx|` histogram table 并画 histogram subplot，用 repeat count 或 rate 表示偏离中心 tooth 的样本数；如果 tooth 分布跨度较大导致柱子太细，可以在脚本头部提供一个明确的 histogram bin count 参数；bin count 越小，每个 bin 覆盖的相邻 tooth range 越宽。
- 长 histogram table 默认不在命令行打印；命令行只打印 compare table 与 aggregate table，分布细节通过图和 `replayData.toothHistogramTable` 查看。
- histogram / aggregate table 只进入轻量 `replayData`，不保存大中间量，不改变 selection、ranking 或 final result。

### progressbar 与运行日志规则

- 小 MC replay 的外层 repeat loop 必须默认显示 progressbar，不再提供 `progressEnable` 开关。
- `progressbar('reset', totalCount)`、`progressbar('advance')`、`progressbar('end')` 必须成对维护；异常 fallback 只能关闭进度条，不能影响 repeat 结果。
- `parfor` 中不能直接调用 `progressbar`；必须在 client 侧创建 `parallel.pool.DataQueue`，用 `afterEach(queue, @(~) progressbar('advance'))` 更新。
- 串行 `for` loop 可以在每个 repeat 结束后直接 `progressbar('advance')`。
- progressbar 是运行可视化，不写入 snapshot，不创建 tmp 文件，也不进入 `replayData`。
- checkpointed replay 的 progressbar 只对未完成 task reset；已经完成的 task 通过 resume 计入日志，不再伪 advance。
- progressbar 开始前的长耗时不能靠 `optVerbose` 解释；应在 common replay helper 中默认打印紧凑 stage log，例如 option resolve、context build、repeat mode、enter parfor。
- `optVerbose=true` 只用于 estimator / flow 内部 trace，不作为 replay orchestration 进度显示开关。
- 清理只能删除无效工程外壳，不能改变默认数值路径、reference-sat 语义、subset 顺序、candidate ranking 或 final selection。
