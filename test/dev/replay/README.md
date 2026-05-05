# replay 脚本说明

`test/dev/replay/` 中的文件是人工排障用的 section 化脚本，不是 regression，也不是可复用 API。每个脚本固定一个代表性问题或小 Monte Carlo 诊断，用于输出 summary、保存轻量 `replayData`，并根据 `replayData` 中的轻量曲线、曲面或表格数据重新画图。

详细运行结果、长表格、观察现象和 snapshot 绑定统一放在 `test/dev/replay/results/`；本 README 只保留脚本规范、入口索引和当前一句话结论。

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

长时间 replay 可在顶层脚本结束、失败 `catch` 或 snapshot 保存后调用 `utils/io/notifyTelegram.m` 发送可选 Telegram 通知。通知必须是 best-effort：发送失败只 `warning`，不得改变 `replayData`、snapshot、图表数据、数值路径或任何 replay 结论口径。

通知内容只保留脚本名、状态、耗时、snapshot 路径和少量关键 summary；详细结果仍写入 `test/dev/replay/results/`，不要在通知 helper 中维护第二套 replay results 解析逻辑。

## 推荐推进顺序

这组 replay 按当前排障主线组织。后续整体仿真、推进和观察现象时，建议按下面顺序运行，而不是按文件名字母顺序运行。

### 1. 先确认 comb / tooth 是真实 objective 结构

#### `replayMfCombToothSurface.m`

- 作用：固定样本可视化 `1/T_f` fdRef comb 和 DoA-fdRef surface。
- 先看现象：wrong-tooth 是否沿 `1/T_f` 重复，truth / final 是否在不同 tooth 上。
- 主要输出：fdRef line、folded comb、DoA-fdRef heatmap、DoA-fdRef mesh，以及 truth-centered DoA-fdRef coupling heatmap。
- 耦合诊断：头部 `couplingFdHalfWidthHz / couplingDoaHalfWidthDeg / couplingFdGridCount / couplingDoaGridCount` 可调；输出 `rhoDoaFd`、ridge slope 和 quadratic ridge slope，用于判断 truth 附近 profile surface 是否近似轴对齐。
- 当前结论：comb / wrong-tooth 是 objective 层真实结构；final periodic refine 更适合同齿 polish，不是跨齿 rescue。
- 详细结果：`results/replayMfCombToothSurface.md`。

### 2. 再确认 subset 负责选齿，periodic 负责同齿细化

#### `replayMfPeriodicVsSubsetToothSelect.m`

- 作用：对比 periodic wide、selected subset、final periodic refine 的选齿与细化分工。
- 配置：头部只保留 `snrDb / baseSeed / numRepeat / saveSnapshot / optVerbose / checkpointEnable / toothHistogramBinCount`，并用英文注释说明；共享场景、卫星、TLE、waveform 与 frame offset 由 `runSimpleDynamicFlowReplayBatch` 调用 `buildDynamicDualSatEciContext` 构造。
- 临时存档：每个 repeat 运行时间较长，头部 `checkpointEnable=true` 时启用 per-repeat checkpoint / resume；设为 `false` 时不创建 tmp、不恢复、不清理。任务文件写入仓库根目录 `tmp/replayMfPeriodicVsSubsetToothSelect/<runKey>/task/task_*.mat`。
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

### 4. tooth 已对后，检查 same-tooth DoA 坏盆地

#### `replayMfInToothFdRangeOracle.m`

- 作用：把参考星 `fdRef` 搜索盒固定到 truth-centered 半齿范围，并用同一批 seed 对比 `SS-SF-Static / MS-SF-Static / SS-MF-CP-U / MS-MF-CP-U / MS-MF-CP-K`。
- 选帧口径：该 replay 已经固定 tooth，不再做 subset tooth selection；repeat fixture 只构造 periodic all-frame view，跳过 curated / random subset fixture bank。
- 先看现象：当 wrong-tooth 被 oracle 排除后，`MS-MF-CP-U-in-tooth` 是否优于单星、单帧与 wide baseline；若它仍不优于这些上限基线，说明论文主张上限不足，应先回查模型层级或 same-tooth refine，而不是继续扩大 subset bank。
- 主要输出：method aggregate、MS-MF wide-vs-oracle pair compare、paper-claim upper-bound compare、runtime timing summary、oracle fd/fdRate range preview、`singleMultiCompareTable`、`tailCaseTable` 与两张分布 / tail 诊断图。
- 表格口径：`aggregateTable` 额外保留 `satMode / frameMode / modelClass / rateMode / oracleLevel / wallTimeMedianMs`，`paperCompareTable` 以 `MS-MF-CP-U-in-tooth` 为 target，正 gain 表示 target 优于对应 baseline；`singleMultiCompareTable` 逐 seed 比较 `SS-MF-CP-U-in-tooth` 与 `MS-MF-CP-U-in-tooth`；`timingTable / timingAggregateTable` 只记录 repeat 级轻量 wall-clock 计时，用于判断 static bundle 与 dynamic method 的耗时占比。
- 命令行口径：`rangeTable` 仍完整保存在 `replayData`，但命令行只打印前四行和后四行，避免 repeat 较多时长表刷屏。
- 图形口径：默认不再画 method-level RMSE bar，也不再叠加多个方法的直方图；改为 angle / fdRef / fdRate 经验 CDF、单独的 MS-vs-SS angle gain 直方图与 paired SS-MF vs MS-MF seed scatter。
- 当前结论：只作为论文主张上限 replay，不进入默认 flow / regression，不改变 no-truth-leak selector。

#### `replayMfInToothTailCaseDiagnose.m`

- 作用：固定 `replayMfInToothFdRangeOracle` 暴露出的 tail seed，默认回放 coherence-collapse tail seed `277 / 283 / 298 / 256` 与负样本 `293 / 280 / 268 / 284`，重放 `SS-MF-in-tooth / MS-MF-in-tooth / MS-MF-truth-DoA-oracle` 等核心 case，并分类 same-tooth tail。
- 选帧与初始化口径：该 replay 已经在 in-tooth oracle 范围内定位 tail，不做 subset tooth selection，也不使用 curated subset 初始化；DoA seed 只来自 static seed 或 truth-DoA oracle；`contextBaseSeed` 固定为源 oracle replay 的 base seed，避免因为只改 seed 列表顺序而换掉共享场景上下文。
- 主要输出：`identityTable`、`aggregateTable`、`tailDiagnosisTable`、`detailDiagnosticTable`、`candidateProbeTable`、`candidateWinnerTable`、`lineProbeSummaryTable`、`rescueBankSummaryTable`、`rescueBankAggregateTable`、oracle range / runtime summary，以及 gated rescue 前后 angle / coherence 对比图。
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
- 诊断口径：重点看每个 center/surface 的 `minDoaOffsetDeg`、`minFdRefOffsetHz`、`boundaryHit`、`p95AbsDoaOffsetDeg` 与 `p95AbsFdOffsetHz`；同时用 `candidateCoverageTable` 判断 implementable center 是否覆盖 truth-DoA oracle 的 hit，用 `policyRecommendationTable` 离线比较 gate threshold、family DoA step 与 safe coherence 的组合。若 hard seed minimum 频繁贴边界，应先扩大 replay envelope，不把该范围直接写进 rescue bank。
- 当前结论：新增为搜索范围决策入口；默认两阶段运行，避免把全部 scout seed 都做重 surface。它不替代 `replayMfInToothGatedRescueEffectiveness` 的效果验证，只为后者是否需要 joint DoA-fdRef candidate bank 提供范围依据。
- 详细结果：`results/replayMfInToothDoaFdRangeEnvelope.md`。

#### `replayMfInToothGatedRescueEffectiveness.m`

- 作用：在更多 repeat 上验证 in-tooth 条件下的 no-truth gated rescue 是否稳定救回 same-tooth collapse，并确认 easy / fd-not-healthy 负样本是否被误伤。
- 选帧与初始化口径：使用 truth-centered half-tooth `fdRef` 范围，只验证 tooth 已正确时的 basin-entry rescue；不做 subset tooth selection，也不改变默认 flow。
- 候选口径：`caseRole / isHardRescued / isDamaged` 可用 truth 做离线评价；`rescueTriggered / triggerReason / selectedCandidateFamily` 只能由默认估计的 non-ref coherence、phase residual、candidate objective 和卫星几何回代诊断决定。当前 gate 采用 `coherence-v1`，只用 default non-ref coherence collapse 触发；phase residual 仅保留字段占位，不参与本 replay 的触发。
- 主要输出：`inToothDoaLadderTable`、`rescueEffectAggregateTable`、`rescueEffectVerdictTable`、`rescueEffectCaseTable`、`rescueBankDecisionTable`、`triggerReasonTable`、`hardMissSummaryTable`、`candidateProbeTable`、method / range / timing summary，以及按 task seed 画出的 default / gated / blanket / truth-DoA oracle angle / coherence 线图和 trigger / case-role 线图。
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
- Telegram 通知：头部 `notifyTelegramEnable=true` 时在 snapshot 保存后发送 HTML `DONE` 通知，在 batch / storage 失败时发送 `FAILED` 通知并继续 `rethrow`；通知只包含脚本名、耗时、seed、snapshot 与少量 ordinary summary，发送失败只 warning。
- 当前结论：新增为 controlled in-tooth DoA 主线的 ordinary-miss 分叉诊断入口；用于判断剩余 `>0.002 deg` miss 是否可由 small polish 或 basin-entry truth-free 候选解释，不据此直接放宽 hard-collapse rescue gate。

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

### 5. unknown-rate release 路线说明

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

### 6. 最后检查 compact flow 与 representative/fullDiag 是否分叉

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
