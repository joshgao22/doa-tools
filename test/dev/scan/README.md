# scan 脚本说明

`test/dev/scan/` 中的文件用于较重的参数、schedule、曲线或曲面扫描。scan 不是 regression，也不是 replay：它通常扫描一个维度或二维网格，用于解释机制和生成可复查的 `scanData`。

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
2. `Build context and flow options`：构造场景、truth、fixture、model option 或 flow option。
3. `Run scan batch`：执行扫描，可在外层使用 `parfor`，不可用时自动退回串行。
4. `Data storage`：构造轻量 `scanData`，用 `saveExpSnapshot` 只保存 `scanData`，并清理 tmp。
5. `Summary output and plotting`：只依赖 `scanData` 输出 summary 并画图。恢复 snapshot 后可直接运行这一节。
6. `Local helpers`：仅放本脚本私有 glue、summary 和绘图 helper。

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
- 主要输出：不同 phase tying 的 comb 曲线、folded tooth 图。

#### `scanMfCpIpPerfMap.m`

- 作用：扫描 CP/IP 在 SNR 和 frame count 上的性能。
- 用途：形成 CP/IP 性能图候选。
- 主要输出：RMSE / hit-rate surface、summary table。

#### `scanMfKnownUnknownInformationLoss.m`

- 作用：扫描 known/unknown Doppler-rate 条件下的信息损失。
- 用途：支撑 nuisance-rate 的 EFIM / information-loss 解释。
- 主要输出：known/unknown performance gap、loss surface。

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

- 作用：观察 subset bank 覆盖率；这是当前 curated / random bank 系统比较入口，不再另建 `scanMfCuratedSubsetScheduleSearch.m`。
- 用途：判断 `curated3`、`curated4`、`random1` 与 full rescue bank 是否覆盖 hard case，random/rescue 是否有必要。
- 默认策略：默认 `strategyPreset="cheapScreen"`，只跑 `curated12`、`curated123`、`curated124`、`curated1234`、`curated12_random1`、`fullRescue`；需要完整二轮确认时把 `strategyPreset` 改为 `"full"`，再比较 `curated12_random{1,2,4}`、`curated123_random1`、`curated124_random1` 等更重组合。
- 主要输出：`aggregateTable`、`scanTable`、`candidateTable`、`candidateSeedCoverageTable`、`transitionTable`、`toothHistogramTable`、subset label selected/evaluated 统计、integer-tooth / residual-aware strict 命中率、same-tooth residual fail、相对 `curated12` 的 rescue/damage transition、候选是否存在但未被 winner 接住、候选评估成本、angle RMSE/P95/max 和 tooth 分布图。
- 评价口径：truth 只用于离线评价 schedule/bank 覆盖，不进入 selector；`truthToothIndexHitRate` 只看整数齿，`truthToothHitRate` 继续表示 residual-aware strict 命中，`sameToothResidualFailRate` 单独记录“已回到 tooth=0 但 residual 未收好”的样本；`transitionTable` 按 seed 对比 `curated12` 到其它 strategy 的 rescue / damage 类型；`candidateSeedCoverageTable` 用于判断好候选是否存在但没被 selector/adoption 接住。该 scan 复用 `runSimpleDynamicFlowReplayBatch` 以保持 repeat 构造、static seed 和 simple-flow 执行与 replay 完全一致；命令行只打印 aggregate、checkpoint summary 与预览，完整候选表、transition、candidate coverage、repeat 表、checkpoint summary 和 histogram 保存在 `scanData`。
- 加速口径：默认使用 cheap screen，不在第一轮直接跑完整 bank；若 cheap screen 中 `curated123` / `curated124` 已接近 `fullRescue`，后续只对 top 2-3 个 strategy 增大 repeat 做 confirmation，不继续扩大 random 数量。
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
