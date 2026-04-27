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
- 用途：说明为什么当前论文路线不是 full dynamic state，也不是 static Doppler。
- 主要输出：窗口量级、DoA drift、Doppler drift、二次相位强度。

#### `scanMfBlockLength.m`

- 作用：扫描 pilot block length 对 comb 和动态可辨性的影响。
- 用途：确认同步块长度变化是否改变 tooth 可辨性和 CP/IP 对比强度。
- 主要输出：block-length summary、fdRef comb 曲线、known/unknown 对比。

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

- 作用：观察 subset bank 覆盖率。
- 用途：判断 curated bank 是否覆盖 hard case，random/rescue 是否有必要。
- 主要输出：subset label coverage、selected/evaluated 统计。

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
- 不保存图片文件，只保存可重画图的数据，例如 grid、curve、surface、summary table。
- tmp 只作为运行时临时目录；正常完成后必须清理。失败时由 `catch` 打印现场路径并保留。
- snapshot 默认只保存 `scanData`，保存路径为 `test/data/cache/scan/`；不要保存大体量原始观测、fixture 全量缓存或 transition bundle。
- 详细运行结果、曲线/曲面观察和 snapshot 绑定记录在 `test/dev/scan/results/`，不写进本 README。

## 恢复后重出结果

```matlab
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开对应 scan 文件，直接运行 `Summary output and plotting` 小节即可重新输出表格和图。

## 头部开关收敛规则

- 不再写 `saveMaxVarBytes`、snapshot 输出目录或 snapshot 前缀；统一交给 `saveExpSnapshot` 的 scan task type 路径处理。
- 不再写并行开关、自动开池开关或最小并行网格阈值。scan 网格 / batch 默认优先用外层 `parfor`，不可用时自动串行。
- 不再写进度条开关。长循环默认显示 progressbar。
- figure 和 table 默认全部输出，不再维护 `showXXXFigure` 或 `showXXXTable`。
- 不再使用通用 context / parallel override。scan 参数差异很大，真正需要改的实验参数应直接写在本脚本 `Scan configuration` 中。
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
