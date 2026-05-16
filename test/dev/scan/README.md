# scan 脚本说明

`test/dev/scan/` 中的文件用于较重的参数、schedule、曲线或曲面扫描。scan 不是 regression，也不是 replay：它通常扫描一个维度或二维网格，用于解释机制和生成可复查的 `scanData`。

## 模板入口

新增 scan 时，优先复制：

```text
test/dev/scan/template/scanTemplate.m
```

到：

```text
test/dev/scan/scanYourTopic.m
```

然后只修改 `Scan configuration`、task grid、summary 和 plot。模板只统一 section、checkpoint / snapshot / Telegram 壳和轻量 `scanData` 组织方式；不要把具体 strategy、schedule、candidate selection、resolved/outlier 分类或 metric parser 提前公共化。短 scan 保持 `checkpointEnable=false`，重 scan 才启用 `tmp/<scanName>/<stableRunKey>/` checkpoint。

新增或改 scan 时，优先复用 `test/common/scan|report|summary` 中已有 header、checkpoint、snapshot、notify 和 compact table 壳；不要复制 replay / scan 中已有的工程 helper。具体扫描指标、plotData 和尚未稳定的分类逻辑继续留在脚本 local 或 results 文档。

## 统一脚本格式

每个 scan 文件头部采用固定形式：

```matlab
% English purpose comment.
% English usage / storage comment.

clear; close all; clc;

%% Scan configuration
```

默认参数直接写在 `Scan configuration` 小节中。不要使用 workspace override，也不要新增默认参数 helper。不同 scan 的参数语义差异很大，统一 opt/override 只会增加误用风险。

## 固定 section 顺序

1. `Scan configuration`：显式写扫描参数、seed、grid、是否保存 snapshot。
2. `Build context and scan tasks`：构造场景、truth、fixture、model / flow option 和显式 task grid。
3. `Run scan batch`：执行扫描，可在外层使用 `parfor`，不可用时自动退回串行。
4. `Data storage`：构造轻量 `scanData`，用 `saveExpSnapshot` 只保存 `scanData`，并清理 tmp。
5. `Summary output and plotting`：只依赖 `scanData` 输出 summary 并画图。恢复 snapshot 后可直接运行这一节。
6. `Local helpers`：仅放本脚本私有 glue、summary 和绘图 helper。


## 公共工程 helper

scan 顶层工程外壳可使用 `test/common/scan/` 的薄 helper：

- `printMfScanHeader` / `printMfScanSection` 统一头部与 section 打印；
- `notifyMfScanStatus` 统一 best-effort HTML Telegram 状态壳；
- `finalizeMfScanResult` 在成功构造轻量 `scanData` 后清理普通 tmp 目录，短 scan 可传空 run dir。

这些 helper 不解析具体 `scanData` 表格、不构造论文指标、不处理 strategy / schedule / candidate / resolved 逻辑；具体 metric line、summary table 和 plot 仍保留在各 scan 本地。

## 可选运行通知

长时间 scan 可在顶层脚本结束、失败 `catch`、checkpoint 完成或 snapshot 保存后调用 `utils/io/notifyTelegram.m` 发送可选 Telegram 通知。通知必须是 best-effort：发送失败只 `warning`，不得改变 `scanData`、checkpoint / snapshot、图表数据、数值路径或任何 scan 结论口径。

通知内容只报告脚本名、状态、耗时、snapshot 路径和少量关键 summary；详细扫描结果仍写入 `test/dev/scan/results/`，不要在通知 helper 中维护第二套 scan results 解析逻辑。

## 推荐推进顺序

这组 scan 按“论文主线 + 当前 dynamic flow 排障”组织。整体仿真时建议先用较轻 scan 建立机制，再跑重的性能图。

### 1. 先确认论文 regime 与模型层级

#### `scanMfRegimeMapByWindow.m`

- 作用：扫描窗口长度、帧间隔和帧数下的 DoA-slow / Doppler-dynamic 区间。
- 用途：说明为什么当前论文路线不是 full dynamic state，也不是 static Doppler，并解释固定 DoA、静态 Doppler、一阶 Doppler 三种模型层级的适用边界。
- 主要输出：窗口量级、DoA drift、Doppler drift、静态 Doppler 相位失配、一阶 Doppler 残差、DoA 容忍度敏感性表、model boundary summary 与连续曲线图。
- 阈值口径：`primaryDoaSlowTolDeg` 只用于当前 target regime 判定；`doaSlowTolDegList` 用于同时给出更严格和更宽松的 DoA static 容忍度边界，避免把某个角度阈值误写成 estimator 精度指标。

#### `scanMfBlockLength.m`

- 作用：扫描 pilot block length 对 comb 和动态可辨性的影响。
- 用途：确认同步块长度变化是否改变 tooth 可辨性、known/unknown-rate 中心差异和块内 Doppler 动态量级；默认保留 `2112` 样本基准，并新增 `2x/4x` 长块用于验证更长已知块的影响。
- 默认波形：`sampleRate=512e6`、`symbolRate=128e6`、`osf=4`、`baseBlockLen=2112`，生成最长 `standardBlockLen=8448` 样本。
- 主要输出：block-length aggregate summary、alias-tooth table、in-block Doppler drift/quadratic-phase 指标、相对最长块与相对 `2112` 基准的 block ratio、fdRef comb 曲线和 known/unknown 对比。
- 存储口径：`scanData` 只保留轻量 summary、tooth table、center line plot data、checkpoint summary 和配置；不保存 `pilotWave`、`viewMs`、`rxSigCell` 或完整 fixture。
- checkpoint：默认开启 block-length 级 checkpoint，路径为仓库根目录 `tmp/scanMfBlockLength/<stableRunKey>/`；中断后用同一配置直接重跑可恢复已完成 block，成功构造 `scanData` 后默认清理 checkpoint 目录。

### 2. 再看 CP/IP 与 known/unknown 的论文主线图

#### `scanMfCpIpTying.m`

- 作用：比较 CP / relaxed / IP tying 对 fdRef comb 的影响。
- 用途：解释连续相位 tying 为什么是论文主模型，而 IP 只是对比基线。
- 主要输出：不同 phase tying 的 comb 曲线、folded tooth 图、phase/center summary table 与 alias-tooth table。
- 存储口径：默认保存轻量 `scanData` 到 `test/data/cache/scan/`，只包含配置、summary table、alias-tooth table 和可重画曲线数据，不保存 `rxSigCell`、完整 `sceneSeq` 或图片。

#### `scanMfCpIpPerfMap.m`

- 作用：扫描 CP/IP 在 SNR 和 frame count 上的性能。
- 用途：形成 CP/IP 性能图候选；只比较同一 static seed 下 CP-K / CP-U 与 IP-K / IP-U 的估计性能，不承担 dynamic flow 正确性验证；不同 frame count 使用本脚本内的 frame-count-aware subset schedule，避免复用 10 帧 curated bank 时越出当前 master window。
- 主要输出：`perfTable`、`aggregateTable`、轻量 `repeatOutCell` 与可重画的 `plotData`；命令行只打印 aggregate table，完整逐 repeat 表留在 `scanData`。
- 存储口径：默认不创建 tmp、不保存图片、不保存 snapshot；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。

#### `scanMfCpIpInToothPerfMap.m`

- 作用：扫描受控 in-tooth 条件下 CP/IP 在 SNR 和 frame count 上的性能。
- 用途：在 truth-centered half-tooth `fdRef` range、truth-local unknown-rate range 与同一 `MS-SF-Static` DoA seed 下隔离 full-flow tooth selection / same-tooth adoption 污染，形成 CP/IP 论文性能图的受控候选；该结果只代表 given-correct-tooth 的模型上限，不代表完整系统性能。
- 主要输出：`perfTable`、`aggregateTable`、轻量 `repeatOutCell`、`checkpointSummaryTable` 与可重画的 `plotData`；表中保留 `inToothMode`、`fdRateRangeMode`、truth-tooth hit、non-ref coherence floor 和 IP/CP ratio。命令行只打印 aggregate table 和 checkpoint summary，完整逐 repeat 表留在 `scanData`。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanMfCpIpInToothPerfMap/<stableRunKey>/`。中断后用同一配置直接重跑可恢复已完成 task；成功构造 `scanData` 后默认清理 checkpoint 目录，失败时 `catch` 打印保留路径。
- 存储口径：默认不保存图片；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。

#### `scanSfStaticMleCrbConsistency.m`

- 作用：扫描 single-frame static DoA-Doppler MLE 与 pilot/static CRB 的一致性。
- 用途：作为 static paper-facing 锚点，固定 SS/MS、DoA-only/static 四个 canonical case 的 full / resolved / CRB-local RMSE、P95、keep rate 与 CRB-normalized gap，为后续 dynamic MLE-vs-CRB 和 outlier 解释提供同口径对照。
- 主要输出：`repeatTable`、`crbTable`、`perfTable`、`aggregateTable`、`crbLocalSummaryTable`、`topTailTable`、`topTailExportTable`、`checkpointSummaryTable` 与可重画 `plotData`；`crbTable` 额外包含 `angleCrbDoaUnit*`、`angleCrbDoaPilot*` 与 pilot/unit FIM ratio 字段；命令行打印 aggregate、CRB-local compact summary、top-tail 预览和 checkpoint summary，完整逐 repeat 表留在 `scanData`。
- 统计口径：full 统计使用所有 finite estimate；resolved 统计使用 estimator `isResolved`；CRB-local 统计在 resolved 样本上按固定 normalized angle / fdRef cap 剔除极端 tail，并报告 keep rate / outlier rate。DoA-only case 默认使用 pilot-model effective-gain DoA-only CRB；`crbTable` 同时保留 unit-gain DoA-only CRB 与 matched DoA-Doppler CRB 对照。Static DoA-Doppler case 仍使用 matched `crbPilotSfDoaDoppler`。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanSfStaticMleCrbConsistency/<stableRunKey>/`。中断后用同一配置直接重跑可恢复；成功构造 `scanData` 后默认清理 checkpoint 目录。
- 存储口径：默认不保存图片；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。

#### `scanMfMleCrbInToothConsistency.m`

- 作用：扫描 Doppler-aided / in-tooth 条件下 CP-K / CP-U local MLE 与 CRB 的一致性。
- 用途：作为论文 MLE-vs-CRB resolved-regime 主图候选；在 truth-centered 单 tooth `fdRef` range 与 truth-local unknown-rate range 下，报告 full-sample、loose resolved、core resolved、trimmed core、CRB-local、CRB-normalized error、keep rate 与 outlier reason。
- truth 口径：truth 默认只用于构造 oracle in-tooth 频率盒、known-rate 真值条件、CRB 计算和离线 resolved / trimming / top-tail 评价；默认 `initMode="auto"` 不使用 static 结果或 truth-frequency `initParam`。`initMode="staticTruthFreq"` 仅作为诊断 oracle baseline，不能作为论文主估计口径。
- 主要输出：`perfTable`、`aggregateTable`、`crbLocalSummaryTable`、`failureSummaryTable`、`topTailTable`、`topTailExportTable`、轻量 `repeatOutCell`、`checkpointSummaryTable` 与可重画 `plotData`；命令行打印 aggregate、CRB-local compact summary、failure reason、top-tail 与 compact export 预览和 checkpoint summary，完整逐 repeat 表留在 `scanData`。`topTailTable` / `topTailExportTable` 增加 `tailSubtype`，用于区分 `same-tooth-fdref-tail`、`outside-tooth-tail`、`fd-boundary-tail`、`angle-local-tail` 等 seed-level 轻量诊断字段，不保存 estimator 内部大 trace。
- 方法选择：配置区保留完整 `methodNameList=["SS-MF-CP-K", "SS-MF-CP-U", "MS-MF-CP-K", "MS-MF-CP-U"]`；需要先看动态单星时可直接注释掉 MS 两行，需要只看多星时可注释掉 SS 两行。该列表只影响 scan 运行的 method bank，不改 estimator 默认路径。
- 初始化口径：默认 `initMode="auto"`，每个 MF 方法内部自行用 frame-aggregated MUSIC/Bartlett 形成 DoA seed，并用 reference-sat single-sat frame-line / continuous-phase refiner 形成 `fdRef/fdRate` seed；不共享 upstream static 估计结果。若需要复现早期 oracle-local 结果，可手动改为 `"staticTruthFreq"`。
- 统计口径：loose resolved 不使用 angle error；core resolved 在 multi-sat case 中进一步要求 non-ref coherence floor 达标；当前 `crbLocal*` / trimmed-core 使用 angle-tail trim，不再用 `fdRef` error 决定 `fdRef` MSE，以避免 trimmed `fdRef` MSE 被条件筛选压到 CRB 以下；MSE/CRB 聚合统一使用逐样本 `mean((error/CRB)^2)`，RMSE/CRB 为其平方根，不再用 `MSE/median(CRB)^2` 作为主口径；`fdRefNormTail` 与 tail subtype 单独标记 angle 合格但 `fdRef` 超 cap 的样本，`crbFloorSummaryTable` 标出 MSE/CRB 低于 1 的行，`crbTargetAuditTable` 标出 CRB-local angle/fdRef 是否通过 1.1 目标。
- 图形口径：默认画 static perf 风格的 angle/fdRef 并列 RMSE-vs-CRB 图、angle/fdRef 并列 MSE/CRB 图，以及 loose/core/CRB-local keep-rate 图；不保存图片。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanMfMleCrbInToothConsistency/<shortRunKey>/`。`shortRunKey` 只保留 frame / SNR range / seed range / repeat 和 8 位 hash，完整语义签名写入 checkpoint meta，避免 Windows 路径超过 260 字符；中断后用同一配置直接重跑可恢复已完成 task，成功构造 `scanData` 后默认清理 checkpoint 目录。
- 存储口径：默认不保存图片；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。

#### `scanMfMsMleCrbInToothConsistency.m`

- 作用：MS 专用的 Doppler-aided / in-tooth MLE-vs-CRB 跨 SNR scan，只跑 `MS-MF-CP-K` 与 `MS-MF-CP-U`。
- 用途：在不混入 SS 对照和 replay heavy bank/rescue solver 的前提下，统计 `-15:5:10 dB` 下 MS-MF 的 full / resolved / core / CRB-local 表现、outlier 类型分布和代表性 tail seed，用于反向指导 `replayMfMsMleCrbFlowDiagnose` 的定点机制排查。
- 默认配置：`snrDbList=(-15:5:10).'`、`numRepeat=40`、`baseSeed=253`、`frameCountList=10`、`dynamicObjectiveMode="pureMle"`、`localSearchBoxMode="truthCrbFloor"`、`localBoxCrbSigmaMultiplier=1`、`localBoxRequestedDoaHalfWidthListDeg=[0.003;0.006;0.012;0.024;0.048]`、`localBoxRequestedFdHalfToothFractionList=0.4`、`methodNameList=["MS-MF-CP-K"; "MS-MF-CP-U"]`。该版本默认 `numRepeat=40` 先做宽度失效点粗筛，确认趋势后再提高到 `100` 或更多。
- 目标函数口径：默认 `dynamicObjectiveMode="pureMle"`，仅在本 scan 调用 dynamic estimator 前将 CP consistency / collapse / negative-projection / non-ref floor 等工程惩罚置零，并关闭 `enableFdAliasUnwrap`，使 MLE-vs-CRB 与无惩罚 deterministic CRB 口径一致；`dynamicObjectiveMode="flowDefault"` 仅用于复现旧工程诊断口径，不作为 P0 paper-facing 主结果。
- 局部搜索盒口径：默认 `localSearchBoxMode="truthCrbFloor"`，这是 controlled oracle-box 诊断而不是 runtime selector。该模式把 DoA 与 `fdRef` 搜索盒中心放在 truth，先按文件头请求半宽构造 hard envelope，再与 `localBoxCrbSigmaMultiplier` 倍对应 MS-MF CRB 标准差取较大值；后续 core solve 只能在该 truth-centered hard bound 内。当前主用法是扫描 `localBoxRequestedDoaHalfWidthListDeg`，从真值附近逐步扩大 DoA 搜索范围，观察从哪个宽度开始出现 bad basin / CRB-local 失效；`localBoxRequestedFdHalfToothFractionList` 默认只有 `0.4 tooth`，如需单独排查 fdRef 宽度再扩成列表。若要复现纯 auto 初始化，把 `localSearchBoxMode="off"`。若要隔离 DoA 或 fdRef，可分别设置 `localSearchBoxMode="truthDoaCrbFloor"` 或 `"truthFdCrbFloor"`。
- 主要输出：继承通用 in-tooth scan 的 `perfTable`、`aggregateTable`、`crbLocalSummaryTable`、`failureSummaryTable`、`topTailTable`、`topTailExportTable`、`repeatOutCell`、`checkpointSummaryTable` 与 `plotData`，并额外输出 `msErrorTypeSummaryTable`、`representativeSeedTable`、`crbFloorSummaryTable` 和 `crbTargetAuditTable`。`crbLocalSummaryTable` 采用 CRB-normalized metric 优先口径，先报告 CRB-local / core / resolved / full 的 `RMSE/CRB` 与 CRB-local / core 的 `MSE/CRB`，再报告 `crbLocalRate` 等 keep-rate 诊断字段；其中 CRB-local trim 为 angle-tail trim，`fdRef` 超 cap 只进入 `fdRefNormTail` / CRB-floor 诊断；宽度扫描轴会进入 `requestedDoaHalfWidthDeg`、`requestedFdHalfToothFraction`、实际 CRB-floor 后 half-width 与 `doaTruthBoxViolationRate` 等字段。
- 错误类型口径：`msErrorTypeSummaryTable` 将逐 repeat 离线分类为 `crb-local`、`fdref-branch-tail`、`fd-rate-tail`、`coherence-collapse-tail`、`doa-basin-limited`、`angle-local-tail`、`solver-conditioning-tail` 或 `other-tail`。这些标签只用于 offline scan / replay seed selection，不进入 estimator、flow selector 或 runtime adoption。
- 代表 seed 口径：`representativeSeedTable` 每个 `displayName × SNR × requested width × error type` 默认选 3 个 tail score 最大的 seed，字段包含 `taskSeed`、CRB-normalized angle/fdRef、fdRate error、non-ref coherence、iterations、first-order optimality、failure reason 与 tail subtype。后续 replay 应优先从这里选 seed，而不是凭单次日志印象挑 case。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanMfMsMleCrbInToothConsistency/<shortRunKey>/`；中断后同配置重跑可恢复，成功构造 `scanData` 后默认清理 checkpoint 目录。
- 存储口径：默认不保存图片；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。该 scan 是 baseline 分类器，不运行 `BankRescue / HealthGated / 0.006+0.012` replay bank。

#### `scanMfMsMleCrbCleanBoundScale.m`

- 作用：在 `replayMfMsMleCrbCleanTrim` 的 6 星 truth-centered clean-bound 口径上，二维扫描 DoA hard box 的 CRB 倍数与 `fdRef` hard box 的 CRB 倍数，观察局部 MLE 的 RMSE/CRB、trim keep rate、range boundary 和 outlier 类型如何随搜索范围变化。
- 用途：回答 CleanTrim 中 `truthLocalDoaCrbScale` / `truthLocalFdRefCrbScale` 是否过紧、是否出现 oracle-bound 截断、以及 DoA/fdRef 范围是否存在耦合；它是 controlled / oracle local range-surface scan，不证明 full-flow acquisition。
- 默认配置：固定 `snrDbList=0`、`numRepeat=30`、`baseSeed=253`、`numFrame=10`、6 星组合 `usrLla=[55;36.59;0]`、`selectedSatIdxGlobal=[5259 1243 348 5652 14 4437]`、`refSatIdxGlobal=5259`，默认只跑 `SS/MS-MF-CP-K/U`。扫描轴为 `truthLocalDoaCrbScaleList=[1;1.5;2;3;5]` 和 `truthLocalFdRefCrbScaleList=[1;1.5;2;3;5;8]`；确认 knee 后再局部加密或扩展 SNR。
- 局部范围口径：`fdRefRangeMode="crb-scale"` 使用纯 CRB-scaled `fdRef` hard box，并以 `fdRefMaxHalfToothFraction` 限制在 tooth 内；这不同于 CleanTrim 中“请求 tooth box + CRB floor”的默认写法，目的是让 `fdRef` scale 列表真正成为扫描轴。`fdRate` 仍使用 truth-local 固定半宽，不参与本 scan 维度。
- 主要输出：`caseTable`、`scaleRawAggregateTable`、`scaleHealthAggregateTable`、`scaleJointTrimAggregateTable`、`scaleReadableRmseCrbTable`、`scaleSurfaceSummaryTable`、`scaleSearchRangeAuditTable`、`scaleOutlierTable`、`runtimeAggregateTable`、`topSlowRuntimeTable`、`checkpointSummaryTable` 和 `plotData`。summary section 只读取 `scanData` 内字段；从 snapshot 恢复后不依赖 workspace 临时变量。
- 解释口径：若小 scale 下 RMSE/CRB 小于 1，应优先解释为 truth-centered hard box 截断 / oracle-bound 条件；若放宽 DoA 或 `fdRef` 后 RMSE/CRB、boundary hit 或 reject reason 明显变化，才说明对应范围参与了 basin / branch 选择。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanMfMsMleCrbCleanBoundScale/<runKey>/`；成功构造 `scanData` 后默认清理 checkpoint 目录。
- 存储口径：`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache，不保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map 或图片。

#### `scanMfMsMleCrbCleanBoundConsistency.m`

- 作用：把 `replayMfMsMleCrbCleanTrim` 固定成 paper-facing 的 truth-centered clean-bound MLE-vs-CRB scan，同时纳入 `SS/MS`、`SF/MF`、`CP/IP`、`K/U` 方法组。
- 用途：服务论文中的局部 MLE 与对应 CRB 对照。在显式给定 DoA / `fdRef` 局部范围后，检查 estimator 是否达到对应模型的 CRB，并把 full / health / joint-trim、SS-vs-MS、K-vs-U、range audit 和 outlier 表一起保存到 `scanData`。
- 默认配置：`snrDbList=(-15:5:10).'`、`numRepeat=200`、`baseSeed=253`、`numFrame=10`、默认 method 包含 `SS/MS-SF-DoA`、`SS/MS-MF-Static`、`SS/MS-MF-CP-K/U`。默认 `methodRangeTable` 对所有 method 使用 `doaCrbScale=2`、`fdRefRangeMode="crb-scale"`、`fdRefCrbScale=2.5`、`fdRefHalfToothFraction=0.3`、`fdRateHalfWidthHzPerSec=1000`，并用 `defaultFdRefMaxHalfToothFraction=0.45` 限制 CRB-scaled fdRef box 仍在同 tooth 内。
- 局部范围口径：配置区的 `methodRangeTable` 是显式 method-level 表格，每行绑定一个 `methodName`，可分别指定 DoA half-width 的 CRB 倍数、`fdRefRangeMode`、`fdRef` half-tooth fraction、`fdRef` CRB 倍数和 `fdRate` truth-local half-width。`fdRefRangeMode="crb-scale"` 时，`fdRefCrbScale` 是实际 half-width 的 CRB std 倍数，`fdRefHalfToothFraction` 只保留为切回 `"tooth-floor"` 时的 tooth-box 宽度；`fdRefRangeMode="tooth-floor"` 时沿用旧口径：先用 `fdRefHalfToothFraction × toothStep` 给宽盒，再用 `fdRefCrbScale × CRB` 做最小 floor。需要精细调参时，直接修改对应 method 行；不要把不同范围作为额外 profile 维度重复跑所有 method，也不向 estimator helper 扩接口。
- IP 口径：脚本在本地 CRB bundle 中补充 independent-phase 的 `IP-K/IP-U` CRB，不改 `buildDynamicCrbBundle` 公共 helper；estimator 侧只通过 `dynOpt.phaseMode='independent'` 调用正式 dynamic estimator。
- 主要输出：`caseTable`、`cleanAggregateTable`、`cleanHealthAggregateTable`、`cleanTrimAggregateTable`、`cleanSsMsCompareTable`、`cleanKnownUnknownCompareTable`、`cleanSnrCurveTable`、`cleanEstimateCrbCurveTable`、`cleanSearchRangeAuditTable`、`cleanOutlierTable`、`runtimeAggregateTable`、`topSlowRuntimeTable`、`checkpointSummaryTable` 和 `plotData`。summary section 只读取 `scanData` 内字段；从 snapshot 恢复后不依赖 workspace 中的临时 aggregate 变量。estimate-vs-CRB 图调用 `enableLegendToggle`，图窗中点击 legend 项可显示 / 隐藏对应曲线，便于临时观察；该交互不影响 `scanData` 或保存口径。
- 限制：这是 controlled / oracle local scan，不证明默认 full-flow acquisition 已通过；`2×CRB` 或其它 truth-centered DoA hard box 必须在论文图注中标明为 local-bound 条件。DoA box、`fdRef` box、geometry、frame count 与 information-loss 曲线的系统扫描应由本 scan 或专门 information-loss scan 承担，不再继续堆到 CleanTrim replay。
- checkpoint：默认开启 per-task checkpoint，路径为仓库根目录 `tmp/scanMfMsMleCrbCleanBoundConsistency/<shortRunKey>/`；成功构造 `scanData` 后默认清理 checkpoint 目录。
- 存储口径：`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache，不保存 `rxSigCell`、完整 `sceneSeq` 或图片。

#### `scanMfKnownUnknownInformationLoss.m`

- 作用：扫描 known/unknown Doppler-rate 条件下的 CRB / EFIM 信息损失。
- 用途：支撑 nuisance-rate 的 EFIM / information-loss 解释，优先作为论文理论机制图候选，不承担 dynamic flow 正确性验证。
- 主要输出：`lossTable`、`primarySliceTable`、`mainSliceTable`、`timeOriginSummaryTable`、`primaryDisplayTable`、`tfSensitivityDisplayTable`、`snrSensitivityDisplayTable`、`crbSummaryTable`、`fimDiagTable`、`fimDiagSummaryTable`、轻量 `crbBundleSummaryCell`。命令行只打印主 SNR + 主帧间隔 + 单一 time-origin class 下的 U/K rollback 百分比表、time-origin 分类提示和 FIM condition summary，完整表留在 `scanData`。
- 图形口径：默认画 3 张聚焦图：主 SNR / 主帧间隔 / 单一 time-origin class 下 `fdRef CRB std rollback (%)` 随帧数变化、主帧间隔下 multi-sat `fdRef CRB std rollback (%)` 随 SNR 变化、主 SNR 下 multi-sat `fdRef CRB std rollback (%)` 对帧间隔的敏感性；不再画信息含义不清的混合 window-span 曲线。若手动加入奇数帧数，`centralRef` 样本会保留在 `scanData.lossTable`，但不会混入默认主趋势图，避免奇偶帧参考时刻不同导致的锯齿误读。
- warning 口径：该 scan 会在调用 CRB bundle 时抑制预期内的 full-FIM ill-conditioned / pinv warning，并用 `fimDiagSummaryTable` 汇总 full-FIM 与 interest-FIM 的条件数；interest-FIM warning 不作为预期噪声压制。
- 存储口径：默认只在 workspace 保留 `scanData`，不创建 tmp、不保存图片；`saveSnapshot=true` 时只保存轻量 `scanData` 到 scan cache。

### 3. 建立 comb / tooth 的 objective 证据

#### `scanMfFdRefComb.m`

- 作用：扫描 fdRef 一维 comb。
- 用途：确认 `1/T_f` wrong-tooth 是 objective 结构，不是单次优化偶然。
- 主要输出：fdRef line、reciprocal peak、folded comb。

#### `scanMfCombTeeth.m`

- 作用：比较相邻 tooth 的 objective 结构。
- 用途：观察 truth tooth 与 wrong tooth 的 objective gap。
- 主要输出：tooth table、tooth objective curve。

#### `scanMfTauSchedule.m`

- 作用：扫描 time-offset schedule 对 comb 的影响。
- 用途：解释 uniform / jittered / gap schedule 对 tooth 等价性的影响。
- 主要输出：tau schedule summary、folded tooth 图。

### 4. 再看 DoA-Doppler coupling 与 same-tooth basin

#### `scanMfTruthNeighborhood.m`

- 作用：扫描 truth、static seed、final estimate 附近的 objective neighborhood。
- 用途：判断最终点是在 truth basin、斜坡还是同齿坏盆地。
- 主要输出：neighborhood surface、center compare。

#### `scanMfDoaToothSlice.m`

- 作用：扫描 DoA 偏移与 tooth 选择的耦合。
- 用途：观察 DoA 错位是否会把 fdRef tooth 拉到错误齿。
- 主要输出：DoA-tooth slice、winner tooth map。

#### `scanMfFdRefFdRateCoupling.m`

- 作用：扫描 fdRef-fdRate 二维 coupling。
- 用途：观察 unknown-rate release 后 fdRef / fdRate 是否存在近等价 ridge。
- 主要输出：fdRef-fdRate surface、alias ridge。

#### `scanMfInToothDoaBasin.m`

- 作用：观察 same-tooth DoA basin 与 polish width。
- 用途：给 conditional very-small polish 的 width 和触发条件提供机制证据。
- 主要输出：in-tooth basin curve、polish-width 对比。

### 5. 再看 subset / flow 机制

#### `scanMfPeriodicVsRandomSubset.m`

- 作用：比较 periodic 与 non-periodic/random subset schedule。
- 用途：解释为什么 subset 选齿、periodic 同齿细化。
- 主要输出：schedule 对比表、tooth hit-rate。

#### `scanMfSubsetBankCoverage.m`

- 作用：观察 subset bank / structured nonuniform schedule 覆盖率；这是当前 legacy curated、random rescue 与 staggered / sparse-ruler / coprime-like schedule 的系统比较入口，不再另建 `scanMfCuratedSubsetScheduleSearch.m`。
- 用途：判断经验 `curated` 是否仍有价值，同时快速筛选更有信号处理解释的非均匀选帧 schedule 是否能生成更好的 Doppler tooth candidate。
- 默认策略：默认 `strategyPreset="scheduleComboQuick"`，只跑 `curated12`、`curated124`、`structuredCombo` 和 `curated12_structured`，用于快速验证 structured 多 schedule 组合是否比单 schedule 更有价值；需要保留 `fullRescue` 上限和 `curated124_structured` 时改为 `"scheduleComboConfirm"`，需要复查单条 structured schedule 时改为 `"scheduleDesign"`，只复查旧核心 bank 时改为 `"diagnosticCore"`，完整 random 二轮确认时改为 `"full"`。
- 主要输出：`aggregateTable`、`scanTable`、`candidateTable`、`candidateSeedCoverageTable`、`consensusTable`、`consensusAggregateTable`、`scheduleFeatureTable`、`transitionTable`、`toothHistogramTable`、subset label selected/evaluated 统计、integer-tooth / residual-aware strict 命中率、same-tooth residual fail、相对 `curated12` 的 rescue/damage transition、候选是否存在但未被 winner 接住、候选评估成本、angle RMSE/P95/max、schedule lag feature 和 tooth 分布图。
- 评价口径：truth 只用于离线评价 schedule/bank 覆盖，不进入 selector；`truthToothIndexHitRate` 只看整数齿，`truthToothHitRate` 继续表示 residual-aware strict 命中，`sameToothResidualFailRate` 单独记录“已回到 tooth=0 但 residual 未收好”的样本；`transitionTable` 按 seed 对比 `curated12` 到其它 strategy 的 rescue / damage 类型；`candidateSeedCoverageTable` 用于判断好候选是否存在但没被 selector/adoption 接住；`consensusTable` 只做 no-truth anchor-group diagnostic，使用 static-anchor residual soft penalty，不再用 residual hard gate 直接丢弃候选 group。该 scan 复用 `runSimpleDynamicFlowReplayBatch` 以保持 repeat 构造、static seed 和 simple-flow 执行与 replay 完全一致；命令行只打印 aggregate、schedule feature、checkpoint summary 与预览，完整候选表、transition、candidate coverage、repeat 表、checkpoint summary 和 histogram 保存在 `scanData`。
- 加速口径：默认使用 schedule combo quick screen，不再第一轮逐个跑所有 structured 单 schedule，也不默认跑 `fullRescue`；先看 `structuredCombo` 与 `curated12_structured` 是否改善 candidate coverage 和 consensus hit，只有有效时再用 `scheduleComboConfirm` 或 24 repeats 做确认。
- checkpoint：默认开启 per-strategy repeat checkpoint，路径为仓库根目录 `tmp/scanMfSubsetBankCoverage/<stableRunKey>/`。中断后可用同一配置直接重跑恢复；成功构造 `scanData` 后默认清理 checkpoint 目录，失败时 `catch` 打印保留路径。

#### `scanMfSubsetRankingLandscape.m`

- 作用：观察 subset ranking landscape。
- 用途：分析 ranking margin、trusted flag 和 selected subset 的稳定性。
- 主要输出：ranking landscape、candidate score table。

#### `scanMfPeriodicTruthNarrowBox.m`

- 作用：观察 truth 附近窄盒 periodic refine。
- 用途：判断 periodic refine 在同齿内是否能作为 polish / refine 的上界参考。
- 主要输出：narrow-box refine table、objective gap。

## results 文档与 snapshot 绑定

- scan 的详细运行结果放 `test/dev/scan/results/<scriptName>.md`。
- 大 `.mat` snapshot 放 `test/data/cache/scan/`。
- 结果文档记录 snapshot 文件名、关键配置、扫描维度、曲线/曲面现象和当前结论。
- 若某个 scan 的结果很多，可将结果文档升级为 `results/<scriptName>/README.md` 与多个子结果文档。
- 新结果先更新对应 results 文档；排障主记录只摘取影响当前优先级的结论。

## 存储与画图规范

- 默认画图，不再提供绘图开关。
- 多阶段、多策略或多 schedule 图必须给每个 subplot 配清楚图例；同一脚本内的 stage / strategy / bank 命名要和表格字段一致，避免图例和 summary 各叫一套。
- 不保存图片文件，只保存可重画图的数据，例如 grid、curve、surface、summary table、candidate table、histogram table 和 representative case。
- tmp 只作为运行时临时目录，且统一位于仓库根目录 `tmp/<scriptName>/<runKey>/`；不能写到 `test/tmp`、当前工作目录或脚本目录。正常完成后必须清理，失败时由 `catch` 打印现场路径并保留。
- 重 scan 若单个 grid / repeat / strategy 很慢，可以在文件头保留一个 `checkpointEnable` 开关。开启时按独立任务写轻量 task 文件，manifest 记录 seed、grid、schedule、strategy、contextOpt 和影响分支行为的 flow signature；关闭时不创建 tmp，直接运行完整 scan。
- checkpoint 只保存独立任务结果，不保存 `rxSigCell`、完整 `sceneSeq`、fixture cell、transition bundle、全量 objective map 或全量 debug trace；成功构造 `scanData` 后默认调用公共 cleanup 入口清理 checkpoint run 目录。scan 脚本不得为 checkpoint run 目录再写私有 cleanup helper。
- 不需要中间落盘的短 scan 不创建 tmp 目录，也不打印 temporary run dir disabled 这类无信息输出；只保留真正影响运行和结果解释的配置。
- snapshot 默认只保存 `scanData`，保存路径为 `test/data/cache/scan/`；不要保存大体量原始观测、fixture 全量缓存、完整 scene / transition bundle 或全量 objective map。
- 详细运行结果、曲线/曲面观察和 snapshot 绑定记录在 `test/dev/scan/results/`，不写进本 README。

## 恢复后重出结果

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开对应 scan 文件，直接运行 `Summary output and plotting` 小节即可重新输出表格和图。

## 头部开关收敛规则

- 不再写 `saveMaxVarBytes`。`saveExpSnapshot` 已有默认 `maxVarBytes`；scan 只通过 `includeVars={'scanData'}` 保存轻量结果，不需要每个脚本重复定义。
- 不再写 snapshot 输出目录或 snapshot 前缀；目录使用 `saveExpSnapshot` 的 scan task type 路径，前缀使用脚本名。
- 不再写并行开关、自动开池开关或最小并行网格阈值。scan 网格 / batch 默认优先用外层 `parfor`，不可用时自动串行。
- 不再写进度条开关。长 grid / surface / repeat / strategy loop 默认显示 progressbar；找不到 progressbar 时只打印一次紧凑提示并继续运行。
- figure 和 table 默认全部输出，不再维护 `showXXXFigure`、`showXXXTable` 或其它细碎显示开关。
- 不再使用通用 context / parallel override。scan 参数差异很大，真正需要改的 seed、grid、schedule、bank、context 或内部并行策略应直接写在本脚本 `Scan configuration` 或对应构造逻辑中。
- checkpoint / resume 只加在确实耗时的 scan 上；短 scan 不为了形式保留 checkpoint 外壳。
- `tmp` 正常结束一定清理；失败时保留现场并打印路径。

## MATLAB 实现细节规则

### `parfor` 广播变量

- scan 通常比 replay 更重，使用 `parfor` 时必须优先消除不必要 broadcast。
- 网格、surface、schedule、candidate bank 进入 `parfor` 前，应展开为按循环变量直接切片的 vector / matrix / cell；不要在 `parfor` 内用 `ind2sub` 产生的二级下标再访问短 grid。
- 必要的只读 model / fixture 可作为 broadcast，但不要把完整 fixture cache、sceneSeq、transition bundle 或全量 objective map 带进 worker。
- 若消除 broadcast 需要改变执行顺序、随机数路径或 candidate ranking，则不要改并行路径，先保持串行一致性。

### 时间戳写法

- 新 scan 不使用 `datestr(now, ...)`。
- 日志与轻量 meta 时间统一用 `datetime('now', 'Format', ...)`；结果文件名交给 `saveExpSnapshot`。

### 冗余变量、分支和 helper 零容忍

- scan 头部只保留真实会调整的实验参数；不维护 `scanOptOverride`、空 context override、空 parallel override、无效 show/plot/table 开关。
- 没有 checkpoint / resume / 中间缓存需求的 scan 不创建 tmp；有 tmp 的 scan 必须正常结束 cleanup，失败时才保留现场。
- 未使用字段、历史兼容分支、重复 table/plotData wrapper、只服务一处的短 helper 应直接删除。
- 重 scan 的公共执行逻辑只有在两个以上入口复用时才提升到 `test/common/scan/`；不要为了“统一”提前增加配置 resolver。

### local helper 注释规则

- scan 脚本中的每个 local helper 至少在函数行后保留一句英文注释，说明它负责什么。
- 注释只解释职责边界和关键语义，不重复逐行代码；短 glue 也要说明为什么存在。
- 新增或清理 helper 时同步检查是否真的需要拆分：一次性短 glue 直接写在执行流程中，grid / surface scan、summary、plot、progress、checkpoint 这类块状逻辑才保留 local helper。

### scan 统计输出规则

- 比较多个 strategy / schedule / bank 的 scan，除了逐 seed / 逐 grid 表格，还应提供能体现分布收敛的 aggregate table。
- tooth selection 类 scan 默认保存 `|toothIdx|` histogram table，并画 histogram 或 rate subplot，用 repeat count 或 rate 表示偏离中心 tooth 的样本数；若 tooth 分布跨度很大，可在脚本头部提供明确的 histogram bin count。
- subset / rescue bank 类 scan 应同时保存 selected label、evaluated label、candidate cost、truth-tooth / near-tooth hit rate、easy-case damage 与 runtime / candidate Pareto 视图。
- 长 histogram table、candidate table 和逐 grid 大表默认不在命令行完整打印；命令行只打印 aggregate table 与必要的紧凑预览，完整表进入 `scanData`。
- histogram / aggregate / candidate table 只进入轻量 `scanData`，不保存大中间量，不改变 selection、ranking 或 final result。

### progressbar 与运行日志规则

- 长 scan 的外层 grid / repeat / strategy loop 必须默认显示 progressbar，不再提供 `progressEnable` 开关。
- `progressbar('reset', totalCount)`、`progressbar('advance')`、`progressbar('end')` 必须成对维护；异常 fallback 只能关闭进度条，不能影响 scan 结果。
- `parfor` 中不能直接调用 `progressbar`；必须在 client 侧创建 `parallel.pool.DataQueue`，用 `afterEach(queue, @(~) progressbar('advance'))` 更新。
- 串行 `for` loop 可以在每个独立任务结束后直接 `progressbar('advance')`。
- progressbar 是运行可视化，不写入 snapshot，不创建 tmp 文件，也不进入 `scanData`。
- checkpointed scan 的 progressbar 只对未完成 task reset；已经完成的 task 通过 resume 计入日志，不再伪 advance。
- progressbar 开始前的长耗时不能靠 `optVerbose` 解释；应默认打印紧凑 stage log，例如 context build、strategy resolve、grid build、repeat mode、enter parfor。
- `optVerbose=true` 只用于 estimator / flow 内部 trace，不作为 scan orchestration 进度显示开关。
- 清理只能删除无效工程外壳，不能改变默认数值路径、reference-sat 语义、subset 顺序、candidate ranking 或 final selection。

- `methodRangeTable` uses method-level ranges. The static seed stage is lazy: SS-only runs do not require a `MS-SF-Static` range row.

#### `scanMfPairGeometryCrbAudit.m`

- 作用：在大 TLE 星座文件中做 reference / cooperative satellite set 的几何与 MF-CRB 健康度筛选；不跑 MLE，不遍历全组合。
- 用途：排查当前 `4154+1165` 这类 MS-MF hard pair 是否由选星几何病态、EFIM condition、DoA-fdRef coupling 或 Doppler 线性残差导致，并为后续 `replayMfPairChoiceCompare` 或 paper-facing MS-MF scan 选出更健康的代表性卫星集合。
- 默认场景：显式批处理 `usrLlaList=[45,36.59,0;55,36.59,0]`、`statlink_20260318.tle`、`numFrame=10`、`frameIntvlSec=1/750`。TLE 文件直接交给 `tleread`，要求用户已把 TLE 所在目录加入 MATLAB path，不再在脚本内维护私有 path resolver。脚本会按 `usrLlaList` 外层串行逐位置运行，位置内部继续使用 checkpoint / parfor；若只想跑单位置，保留 `usrLlaList` 的一行即可。每个位置先按参考时刻可见星筛 reference candidate，再用 `refSelectionMode="cooperativeBestScore"` 对 reference 的 second-sat 配合能力做 pair audit，选出 cooperative reference 后固定该 reference 枚举 second-sat candidates，最后用 greedy marginal CRB health 扩展到默认 `L=1,2,4,6,8,10,12,14,16` 的 nested sets。
- 模型口径：只评估 paper-facing 的 global DoA 恒定、Doppler 一阶变化、continuous-phase MF-CP 模型；代码中 `signalModelTag="constant-doa-affine-doppler-cp"`，CRB 使用 `phaseMode='continuous'`、`fdRateMode='known/unknown'` 与拟合得到的 `fdRefFit/fdRateFit`。`steeringMode='framewise'` 只表示同一个固定 `usrLla` 在不同卫星帧下重新映射为 local steering，不表示估计时引入时变 DoA 主参数。
- reference 选择口径：默认 `refSelectionMode="cooperativeBestScore"`，不是单纯选择最大仰角 reference。脚本先评估最多 `maxCooperativeReferenceEval` 个 single-sat reference candidate，并对每个 reference 的 second-sat 候选计算 pair CRB health，形成 `refCooperativeScoreTable`；该表报告 best/top-5 pair score、healthy/hard-coupled pair 数量、best-pair EFIM condition 和 DoA-fdRef coupling。`refSelectionMode="maxElevation"`、`"singleSatBestScore"` 和 `"manual"` 仍保留为对照。
- 排名口径：主 pair 排名使用 `pairScore`，它综合 `MS/SS` angle CRB gain、fdRef CRB gain、MS EFIM condition、DoA-fdRef canonical coupling、fdRef coupling fraction 与 Doppler fitted-phase residual。表中同时保留 `rankByAngleGain`、`rankByLowCoupling` 和 raw CRB/coupling 字段，避免把单一 score 写死成论文结论。
- 主要输出：单位置时保留 `contextSummaryTable`、`referencePreselectTable`、`referenceCandidateTable`、`refCooperativeScoreTable`、`allReferencePairRankTable`、`candidatePreselectTable`、`pairCandidateTable`、`pairRankTable`、`greedyCandidateTable`、`greedyStepTable`、`satSetSummaryTable`、`rankPreviewTable` 和轻量 `plotData`；多位置批处理时额外保留 `locationScanDataCell`，并在聚合表中加入 `locationIdx` 与 `usrLlaStr`，同时提供 `locationSummaryTable` 与 `locationSatSetSummaryTable` 便于比较不同位置的 L=1/2/4/.../16 结果。`checkpointSummaryTable` / `checkpointCleanupReport` 不进入 `scanData`，避免把运行期 housekeeping 写进 snapshot。命令行只打印 reference / selected-set summary 和必要 preview，不打印 9000+ TLE 的全量信息。
- checkpoint / progress / runtime log：默认开启 task-level checkpoint，路径为仓库根目录 `tmp/scanMfPairGeometryCrbAudit/<stageRunKey>/`，checkpoint key 带 `locXX` 前缀，避免批处理不同位置之间互相覆盖。reference、cooperative reference `coopAll`、selected-reference pair 和每个 greedy step 分 stage 保存；cooperative reference audit 不再按 reference 逐个开小 `parfor`，而是把同一位置下的 `maxCooperativeReferenceEval × maxCandidateSecondSat` 个 pair task 合并成一个 `coopAll` checkpoint / parfor 批次，以提高 worker 利用率并减少 stage overhead。中断后同配置重跑可恢复，progressbar 只对未完成 task reset，parfor 时通过 client-side `DataQueue` 或 checkpoint callback 更新；脚本会打印轻量 stage-level runtime log，说明当前处于 location、context、reference、coopAll、pair、greedy、checkpoint cleanup 或 summary 阶段。若运行中断或报错，checkpoint 会保留，下一次同配置可继续 resume；正常完成时只执行 cleanup，不把 cleanup report 或 checkpoint summary 写入 `scanData`。
- 最新结果文档：`test/dev/scan/results/scanMfPairGeometryCrbAudit.md` 记录 2026-05-15 fixed-location 深挖结果；当前推荐 `[55,36.59,0]`、`refSat=5259`、`L=6 [5259 1243 348 5652 14 4437]` 作为论文主仿真组合，`L=4` 作为轻量 smoke，`L=8` 作为增强确认，`[45,36.59,0]` 系列降级为 hard-coupled stress 对照。
- 限制：该 scan 只评估几何与 CRB 健康度，不证明 estimator 在所选 set 上一定贴 CRB；被选出的 healthy / hard-coupled set 需要再用小 MC replay 或 paper-facing scan 验证 MLE / CRB consistency。

#### `scanMfGeoSatGeometryCrbAudit.m`

- 作用：在一组候选地理位置上做 Starlink-inspired location + cooperative satellite-set 初筛；不跑 MLE，不遍历全组合，不替代固定位置的 `scanMfPairGeometryCrbAudit.m`。
- 用途：为论文仿真选择 representative `usrLla` 与候选卫星集合。该 scan 是上层场景筛选器，默认采用 `geoScanMode="coarseL4"`：每个候选 `usrLla` 只做轻量 cooperative reference / pair CRB-health 粗筛，并只对排名靠前的少数位置贪心到 `L=1,2,4`；完整 `L=1,2,4,6,8,10` 深挖交给固定位置的 `scanMfPairGeometryCrbAudit.m`。
- 默认地理位置：第一轮只做轻量纬度 sweep，并保留当前坐标作为 anchor：`[37.78,36.59,0]`、`25/35/40/45/50/53/55/60 deg` 纬度且经度固定 `36.59 deg`。若该轮显示中纬度更健康，再手动扩展经度列表；不要把第一版写成全球经纬网格搜索。
- 模型口径：与 `scanMfPairGeometryCrbAudit.m` 一致，只评估 global DoA 恒定、Doppler 一阶变化、continuous-phase MF-CP 模型；`signalModelTag="constant-doa-affine-doppler-cp"`，CRB 使用 `phaseMode='continuous'`、`fdRateMode='known/unknown'`、`steeringMode='framewise'`。
- 地理位置评分：`geoPreScore` 只作为排序辅助，综合 best cooperative reference score、可见星数量、best-pair EFIM condition 与 DoA-fdRef coupling。表中保留 raw `numAvailableAtRef`、best pair class、best pair CRB、coupling fraction、healthy/hard-coupled pair 数量，避免把单一 score 写成论文结论。
- coarse 口径：默认 `maxReferenceCandidateEval=16`、`maxCooperativeReferenceEval=12`、`maxCandidateSecondSat=8`、`maxGreedyCandidatePerStep=12`，避免每个地理点完整复制单位置 scan 的 `40 + 12×51` task 流程。cooperative reference audit 不再按 reference 逐个开小 `parfor`，而是把同一地理点下的 `maxCooperativeReferenceEval × maxCandidateSecondSat` 个 pair task 合并成一个 `geoXX_coopAll` checkpoint / parfor 批次，以提高 worker 利用率并减少 stage overhead。默认只对 `numGeoGreedyKeep=3` 个地理位置运行 `L=1,2,4` greedy satellite-set scan。输出 `paperCandidateSetTable` 用于挑后续单位置深挖 / MLE-CRB paper-facing scan 起点；它不是最终性能结论。
- 主要输出：`geoCandidateTable`、`geoReferenceSummaryTable`、`referenceCandidateTable`、`refCooperativeScoreTable`、`allReferencePairRankTable`、`geoGreedyStepTable`、`geoGreedyCandidateTable`、`geoSatSetSummaryTable`、`paperCandidateSetTable`、`bestGeoSummaryTable`、`checkpointSummaryTable` 和轻量 `plotData`。命令行只打印地理位置列表、top geo cooperative summary、top geo 的 selected-set summary 和 paper candidate preview。
- checkpoint / progress / runtime log：默认开启 task-level checkpoint，路径为仓库根目录 `tmp/scanMfGeoSatGeometryCrbAudit/<stageRunKey>/`。checkpoint key 带 `geoXX` 前缀，避免不同地理位置下相同 reference / pair stage 互相覆盖；成功构造 `scanData` 后默认清理已完成 checkpoint。
- 与 `scanMfPairGeometryCrbAudit.m` 的关系：本 scan 负责“多个地理位置中选场景”，默认只筛到 `L=4`；`scanMfPairGeometryCrbAudit.m` 负责“给定一个或少量显式地理位置后深挖 reference / pair / `L=1,2,4,6,8,10,12,14,16` greedy set”。两者都应保留；若要解释某个具体位置为什么被选中，或要得到最终 6/8/10 星集合，应回到固定位置 scan 继续看。
- 最新结果文档：`test/dev/scan/results/scanMfGeoSatGeometryCrbAudit.md` 记录 2026-05-14 coarse-L4 地理位置初筛结果；当前结论是 `[45,36.59,0]` 为主候选，`[60,36.59,0]` 为 stress candidate，原 `[37.78,36.59,0]` 降级为 anchor / stress 对照。

