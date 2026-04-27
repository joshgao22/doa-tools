# doa-tools

本仓库用于 DoA / Doppler / CRB 相关研究代码，当前主线是多星连续相位多帧 DoA-Doppler 估计。代码组织目标是：正式算法逻辑集中，测试与排障入口清楚，回归与 replay 职责分开，避免同一现象在多个脚本里重复验证。

## 规范入口

AI 或自动化工具修改代码前，先读 `AGENTS.md`。`AGENTS.md` 负责说明读文档顺序、规范职责和冲突裁决；根 README 负责全局代码风格、注释语言、文件存储和目录职责；更细的局部规则放在目标目录 README 中。

默认读取链条：

```text
AGENTS.md -> README.md -> 目标目录 README -> 更内层 README -> 必要时读 results / 排障记录
```

不维护单独的长 `CODING_RULES.md`。长期规则应放在根 README 或局部 README；具体运行结果应放在 `results/` 文档；当前优先级和机制结论放在排障记录。外部 coding prompt 只约束 AI 的执行风格和风险偏好，不替代仓库内 `AGENTS.md` 与 README 的文件职责。

## 文档分工与冲突裁决

| 文档 | 承担内容 | 不承担内容 |
|---|---|---|
| `AGENTS.md` | 最短读法、优先级、冲突裁决、默认修改原则 | 长实验结果、脚本使用手册 |
| 根 `README.md` | 仓库导航、目录职责、全局风格、落盘规范索引 | 单个 replay / scan 的长结果 |
| 目录 README | 本目录入口、文件职责、运行方式、局部规则 | 历史排障流水账 |
| `results/` 文档 | 具体运行结果、长表格、snapshot 绑定、观察现象 | 通用 coding 规范 |
| 排障主记录 | 当前主问题、优先级、必跑护栏、已证伪路线 | 脚本使用手册和完整结果表 |
| 机制归并版 / 历史归档 | 机制证据和历史回溯 | 日常开发入口 |
| 外部 coding prompt | AI 执行风格、风险偏好、交付说明格式 | 仓库文件索引和长期实验记录 |

若文档之间冲突，先按 `AGENTS.md` 的读法和优先级处理。若冲突会影响默认数值路径、reference-sat 语义、搜索边界、输出字段或 summary 口径，优先选择更小、更保守的修改，并在交付说明中写明取舍。

## 当前论文主线

当前研究围绕以下问题展开：

- 在多帧观测窗口内，DoA 近似慢变，但 Doppler 漂移已经不可忽略；
- 主参数是参考时刻的全局 DoA 与参考星链路上的 Doppler；
- 参考星 Doppler rate 作为 nuisance parameter 进入模型；
- 连续相位 CP 模型是主模型，跨帧独立相位 IP 是对比基线；
- CRB / EFIM 用于解释 known-rate 与 unknown-rate 条件下的主参数信息损失。

因此，代码修改默认围绕这个主线服务，不为与论文主线无关的泛化框架重构接口。

## 常用入口

### 快速检查代码是否回退

- 入口：`test/regression/runRegressionQuick.m`
- 用途：日常改 estimator、flow、scene、summary 后的第一道自动护栏。

### 集中检查 dynamic flow

- 入口：`test/regression/runRegressionDynamicFlow.m`
- 用途：改 subset bank、periodic refine、polish、recovery gate 后跑。

### 跑少量完整 pipeline 护栏

- 入口：`test/regression/runRegressionPipeline.m`
- 用途：改完整 fixture -> estimator -> summary 链路后跑。

### 跑 checkpoint / snapshot / cleanup smoke

- 入口：`test/regression/perf/runRegressionPerfSmoke.m`
- 用途：改长任务工程链路后跑。

### 单次代表性实验

- static：`test/dev/doaDopplerStatDualSatUraEci.m`
- dynamic：`test/dev/doaDopplerDynDualSatUraEci.m`

### 固定问题回放与小 MC 诊断

- 入口目录：`test/dev/replay/`
- 用途：复现 hard case、打印 trace、画图、保存 slim snapshot。
- 详细结果记录：`test/dev/replay/results/`

### 机制扫描与策略比较

- scan：`test/dev/scan/`
- scan 详细结果记录：`test/dev/scan/results/`
- strategy：`test/dev/strategy/`

## 仓库目录职责

### `array/`

通用阵列构造、阵列旋转、steering、array-only snapshot。

- 不依赖 satellite scene 或当前论文的 reference-sat 语义；
- 旧实现保留在 `legacy/`，非主线优先入口。

### `estimator/`

正式估计器入口与 estimator 侧 helper。

- profile likelihood、init、branch solve、winner adoption 等正式算法逻辑放这里；
- dev-only probe、replay、summary 不应长期留在这里；
- helper 分工见 `estimator/helper/README.md`。

### `performance/`

CRB、FIM、EFIM、理论性能界。

- 只放性能界计算；
- Monte Carlo orchestration、表格汇总、plot 数据构造放 `test/`。

### `satellite/`

卫星几何、场景、信号生成、动态 truth。

- reference-sat、scene slicing、Doppler geometry 放 `satellite/scene/`；
- pilot / rx snapshot 放 `satellite/signal/`；
- steering drift、line-fit truth 诊断放 `satellite/dynamic/`。

### `test/`

regression、dev、replay、scan、strategy、实验 helper。

- regression 是自动 pass/fail 护栏；
- replay 是固定问题回放，不是正确性契约；
- scan 是较系统的参数 / schedule / surface 扫描；
- common 是实验侧复用 helper，不长期维护第二套正式算法。

### `utils/` 与 `plottool/`

通用工具和通用绘图。

- 只有确实没有 DoA-Doppler 语义的工具才放 `utils/`；
- 单脚本一次性图先留 local，复用后再提升到 `test/common/plot/` 或 `plottool/`。

## README 索引

### 规范与结果入口

- `AGENTS.md`：AI / 自动化修改代码前的最短入口，说明读文档链条和优先级。
- `test/data/cache/README.md`：大 `.mat` snapshot 的集中存储规则。

### 测试与排障入口

- `test/README.md`：test 总入口，决定跑 regression、replay、scan 还是 strategy。
- `test/regression/README.md`：每个 regression 的唯一契约和 runner 说明。
- `test/dev/README.md`：dev 主入口、子目录职责和脚本选择。
- `test/dev/replay/README.md`：replay 脚本规范、入口索引和当前一句话结论。
- `test/dev/replay/results/README.md`：replay 详细结果记录与 snapshot 绑定。
- `test/dev/scan/README.md`：scan 的扫描维度、输出和相关 replay。
- `test/dev/scan/results/README.md`：scan 详细结果记录与 snapshot 绑定。
- `test/dev/strategy/README.md`：策略比较与已证伪方向记录。
- `test/common/README.md`：test helper 的放置边界。

### 正式实现与理论分析

- `estimator/README.md`：正式 estimator 入口和模型层级。
- `estimator/helper/README.md`：model/init/bounds/branch/warm-anchor helper 分工。
- `satellite/scene/README.md`：reference-sat、scene、Doppler 不变量。
- `performance/README.md`：CRB / FIM / EFIM 与论文分析层对应关系。

## 全局 coding 约束

### 行为保持

默认保持现有函数签名、默认参数、初始化策略、搜索边界、reference-sat 语义、sat 顺序、输出字段、summary 口径和 fallback 顺序。

如果修改可能改变数值行为，交付说明必须写明：

- 哪个模块变了；
- 哪个数值路径变了；
- 风险点在哪里；
- 为什么必要；
- 建议跑哪些 regression / replay。

### MATLAB 风格

- `.m` 文件中的函数头注释、脚本头注释和 inline comment 使用英文。
- README、排障记录和说明文档使用中文。
- 命名和格式沿用当前目录同类型、非 legacy 文件风格。
- 不为小改动新增过重 `arguments`、`validateattributes`、通用 opt resolver 或框架化封装。
- 主入口只做 orchestration；正式复用逻辑按职责放到 estimator / satellite / performance / test common。

### 文件与结果存储

- 运行时临时文件放仓库根目录 `tmp/<scriptName>/<runKey>/`。
- 大 `.mat` snapshot 统一放 `test/data/cache/<taskType>/`，例如 `replay/`、`scan/`、`perf/`。
- replay 保存轻量 `replayData`；scan 保存轻量 `scanData`；perf 保存 final result struct、summary table、关键配置和 meta。
- 不默认保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- 人工可读的结果分析不写进 cache 目录，放到对应任务目录的 `results/` 文档。

常用落盘入口只在 `test/data/cache/README.md` 维护详细说明：

| 操作 | 推荐入口 | 位置 |
|---|---|---|
| 保存最终 snapshot | `saveExpSnapshot` | `utils/io/saveExpSnapshot.m` |
| 恢复 snapshot | `loadExpSnapshot` | `utils/io/loadExpSnapshot.m` |
| 清理 cache / tmp 运行产物 | `cleanupRunArtifacts` | `utils/io/cleanupRunArtifacts.m` |
| 清理 checkpoint run 目录 | `cleanupPerfTaskGridCheckpoint` | `test/common/flow/cleanupPerfTaskGridCheckpoint.m` |

根 README 只保留入口索引；save / load / cleanup 的调用细节不要在其它 README 中重复维护。

### 文档同步

- 新增或重命名 regression / replay / scan / common helper 后，同步更新对应局部 README。
- 新增入口文件、移动 helper、改变结果保存位置、改变 replay / scan 输出字段或改变 summary/table 口径后，同步更新对应 README。
- 只改内部数值实现且入口、输出字段、保存位置不变时，不强制更新 README；若行为变化可能影响使用者判断，仍应在交付说明中说明。
- 新增重要 replay / scan 结果时，先更新对应 `results/*.md`；只把影响当前决策的结论摘到排障主记录。
- README 负责入口、职责和运行方式；排障记录负责当前优先级、机制结论和已证伪路线。

## 放置规则摘要

### 正式算法逻辑

放到：

- `estimator/helper/`
- `satellite/scene/`
- `satellite/signal/`
- `performance/`

不要长期放在：

- `test/dev/` local helper；
- `test/common/flow/` 的 orchestration glue；
- replay / probe 脚本内部。

### 回归护栏

放到 `test/regression/`，前提是：

- 能自动 pass/fail；
- 结论已经稳定为代码契约；
- 失败能指向具体契约；
- 不需要人工看图或比较表格。

### 问题定位

放到 `test/dev/replay/`、`test/dev/probe/`、`test/dev/scan/` 或 `test/dev/strategy/`。

- 固定问题或小 MC：`replay/`；
- 局部曲面 / profile：`probe/`；
- 参数和 schedule 扫描：`scan/`；
- 策略比较和已证伪路线：`strategy/`。
