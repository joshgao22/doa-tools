# regression 目录说明

`test/regression/` 只放自动 pass/fail 的护栏脚本。它的职责不是记录所有排障现象，而是固化当前代码已经承诺满足的契约。失败时应能指向一个明确的语义、分支或流程回退。

## 子目录职责

### `invariant/`

低层不变量、数学结构契约、reference 语义、CRB/FIM 基本关系。

适合放：

- reference-sat 语义；
- `fdRef / deltaFdRef / fdSat` 组合关系；
- user-state 与 local/global DoA 映射；
- CP/IP time-axis 语义；
- CRB known/unknown 基本一致性。

不放：

- 最终性能比较；
- 策略优劣；
- 长 Monte Carlo；
- 需要人工看图的脚本。

### `branch/`

分支选择、guard、fallback、winner adoption、subset ranking、polish trigger 等流程契约。

适合放：

- no-truth-leak；
- final selection rule；
- warm-anchor gate；
- random/wide rescue gate；
- polish 触发条件。

不放：

- replay；
- 曲面图；
- 大策略比较；
- 未修好的 hard case。

### `pipeline/`

少量代表性完整流程护栏。

适合放：

- fixture -> estimator -> summary 的代表性闭环；
- 已稳定且可自动判断的机制结论。

不放：

- 大量 hard case；
- 长任务；
- checkpoint 工程链路。

### `perf/`

长任务工程链路的 smoke-small 护栏。

适合放：

- checkpoint / resume；
- snapshot；
- cleanup；
- task grid；
- summary 保存。

不放：

- 常规 quick 必跑数值 regression；
- 长 Monte Carlo 正式实验。

## Runner

### `runRegressionQuick.m`

- 作用：常规快速护栏。
- 覆盖：invariant 与 branch 中最重要的自动契约。
- 何时跑：每次改 estimator、flow、scene、summary 后优先跑。

### `runRegressionDynamicFlow.m`

- 作用：dynamic-flow 相关组合 runner。
- 何时跑：改 subset bank、periodic refine、polish、recovery gate 后跑。

### `runRegressionPipeline.m`

- 作用：代表性完整流程护栏。
- 何时跑：改完整 fixture -> estimator -> summary 链路后跑。

### `perf/runRegressionPerfSmoke.m`

- 作用：perf 工程链路 smoke runner。
- 何时跑：改 checkpoint、resume、snapshot、cleanup 或 task grid 后跑。

## invariant regression

### `regressionRefStateInvariant.m`

- 唯一契约：reference-Doppler 状态不变量。
- 重点检查：`fdSat(ref)=fdRef`、`deltaFd(ref)=0`、`fdSat=fdRef+deltaFd`。
- 覆盖：full scene 与 subset scene。
- 失败含义：reference 重映射或 Doppler 组合被改坏。

### `regressionUserStateFromLatlon.m`

- 唯一契约：`lat/lon -> user ECI state` 的 motion-aware 语义。
- 重点检查：sceneSeq 与 nominal state 间的运动偏移一致。
- 失败含义：候选 user-state 与生成场景的运动语义不一致。

### `regressionLocalDoaMapping.m`

- 唯一契约：global DoA 与 local DoA / local direction 的闭环一致。
- 重点检查：az/el 映射、方向向量回转、ref-frame local DoA。
- 失败含义：几何坐标链被改坏。

### `regressionSfModelBuild.m`

- 唯一契约：single-frame model builder 的字段与语义稳定。
- 重点检查：reference identity、`satWeight`、DoA bounds。
- 失败含义：SF model 构造契约回退。

### `regressionSfObjectiveShape.m`

- 唯一契约：noiseless SF profile objective 在 truth 附近形状正确。
- 重点检查：truth objective 最优，DoA/fd 扰动变差。
- 失败含义：SF objective 主核或输入构造被改坏。

### `regressionSfStaticRefOnlyEquivalence.m`

- 唯一契约：SS-SF-Static 与 MS-SF-Static zero-weight ref-only 分支等价。
- 重点检查：zero-weight 不应改变 ref-only 结果。

### `regressionSfStaticReferenceDopplerInvariant.m`

- 唯一契约：static reference Doppler 组合关系在多个 case 中成立。
- 重点检查：`deltaFdRef(ref)=0`、`fdSat(ref)=fdRef`、composition error。

### `regressionSfStaticTruthReplayMultiSat.m`

- 唯一契约：multi-sat static truth replay 的 objective additivity 与 fd composition 正确。
- 重点检查：per-sat objective 与 total objective 一致。

### `regressionSfStaticFdInitZeroWeightInvariant.m`

- 唯一契约：zero-weight second sat 不应改变 ref-only fdRef 初始化。

### `regressionMsSfStaticObjectiveDecomposition.m`

- 唯一契约：MS-SF-Static objective 随 `satWeight` 的分解关系正确。

### `regressionMsSfStaticW0EqualsSsStatic.m`

- 唯一契约：MS-SF-Static 的 W0 与 SS-SF-Static 等价。

### `regressionMsSfStaticWeightInjection.m`

- 唯一契约：satWeight 只按预期注入对应卫星 objective。

### `regressionMfModelBuild.m`

- 唯一契约：multi-frame dynamic model builder 的 reference、bounds、unknown vars 结构稳定。

### `regressionMfInitFdLine.m`

- 唯一契约：MF fd line 初始化在 ref-only 与 multi-sat 输入下保持 reference-first 语义。

### `regressionSnapshotTruthReplayMf.m`

- 唯一契约：snapshot truth replay 下 truth residual 与扰动 residual 关系正确。

### `regressionOtherSatOnlyMf.m`

- 唯一契约：other-sat-only MF dynamic 链路单独成立。
- 作用：排除“非参考星天然坏”的误判。

### `regressionMfCpSupportCollapsePenalty.m`

- 唯一契约：CP support / consistency penalty 不应被错误支路轻易绕开。

### `regressionMfCpIpTimeAxisInvariant.m`

- 唯一契约：CP 使用 absolute/global time，IP 使用 frame-local time。
- 重点：只验证 time-axis 语义，不比较 CP/IP 性能优劣。

### `regressionCrbKnownUnknownConsistency.m`

- 唯一契约：known-rate 与 unknown-rate CRB/EFIM 的基本包含关系与量级一致。
- 注意：full-FIM singular / pinv warning 是该 case 的 expected warning。

## branch regression

### `regressionMfSubsetSelectNoTruthLeak.m`

- 唯一契约：subset selection 不能使用 truth leak。
- 重点检查：只改变 `angleErrDeg / fdRefErrHz / toothIdx / truthToothIdx` 等 truth-aware 字段时，selected subset 与 trust flag 必须完全不变。
- 额外覆盖：phase-health gate 优先于 objective 的排序契约。
- 不再另建重复的 same-tooth ranking regression。

### `regressionMfFastSubsetEscalation.m`

- 唯一契约：fast subset escalation / random-wide rescue gate 的触发条件。
- 覆盖：confident、tooth disagree、drift、untrusted、off-tooth、disabled。
- 不再另建重复的 random-rescue regression。

### `regressionMfUnknownFinalSelectionRules.m`

- 唯一契约：unknown final winner selection 的表驱动规则。
- 覆盖：objective-only、subset-anchor guard、fd-wide wins、wide refine、same-tooth wide、closer-to-zero tooth。

### `regressionMfUnknownFixedDoaWarmAnchor.m`

- 唯一契约：fixed-DoA warm-anchor 能实际 release fdRate，并正确报告 `solveVariant`。

### `regressionMfUnknownWarmStartSet.m`

- 唯一契约：unknown warm-start set 的 source tags、fdRate seed 与 DoA half-width 构造稳定。

### `regressionMfUnknownReleaseFromCpK.m`

- 唯一契约：CP-U 从 CP-K seed 出发时能完成 rate release，不退化成单点 seed。

### `regressionMfUnknownBestStartSelection.m`

- 唯一契约：unknown outer start 选择按 objective / improvement 规则稳定。

### `regressionMfBranchKnownUnknown.m`

- 唯一契约：CP-K 与 CP-U branch 主链的 known/unknown 分支语义稳定。

### `regressionMfWarmAnchorParforGate.m`

- 唯一契约：warm-anchor 内层 parfor 默认关闭，只允许显式 opt-in。
- 重点：verbose、fixed-DoA tooth guard、小 release family 都必须保持串行。

### `regressionMfUnknownNoisyWrongToothGuard.m`

- 唯一契约：noisy wrong-tooth guard 不被默认路径破坏。
- 重点：验证 warm-anchor 默认串行与 final adoption guard 的稳定性。

## pipeline regression

### `regressionMfDoaProfileWithHealthyFd.m`

- 唯一契约：在 fd 已健康前提下，DoA profile / final DoA 不应明显偏离 truth basin。
- 作用：当前 same-tooth hard case 的核心机制护栏。

### 其他 pipeline case

- 只放已经稳定的完整流程契约；
- 若还需要人工看图、看表或判断策略优劣，先放 replay 或 strategy。

## 不重复验证规则

### subset selection

- no-truth-leak 与 same-tooth ranking 契约归 `regressionMfSubsetSelectNoTruthLeak.m`。
- 不另建功能重复的 subset ranking regression。

### final selection

- final winner 规则归 `regressionMfUnknownFinalSelectionRules.m`。
- 新的 winner case 应合并进该文件的表驱动 case list。

### random / wide rescue

- fast subset escalation 与 rescue gate 归 `regressionMfFastSubsetEscalation.m`。
- 不另建只改名字的 random rescue regression。
- curated3 / curated4 / random1 的性能和成本比较属于 scan，归 `test/dev/scan/scanMfSubsetBankCoverage.m`，不进入 regression。

### warm-anchor parfor

- 默认 gate 归 `regressionMfWarmAnchorParforGate.m`。
- 敏感性证据放 `replayMfWarmAnchorParforSensitivity.m`，不放 regression。

### same-tooth hard case

- 未完全稳定前放 replay。
- 只有当结论可自动 assert 且 runtime 可接受时，才迁移到 pipeline regression。
