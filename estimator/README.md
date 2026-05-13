# estimator 目录说明

`estimator/` 放正式估计器入口和估计器侧 helper。当前主线是 multi-frame continuous-phase DoA-Doppler 估计，重点是参考星参数化、CP/IP、known/unknown Doppler-rate、profile likelihood 与 branch solve。

## 当前主线入口

### `estimatorDoaMlePilotOpt.m`

- 模型层级：pilot DoA MLE。
- 主要用途：纯 DoA baseline。
- 调用者：static dev、regression。

### `estimatorDoaDopplerMlePilotSfOpt.m`

- 模型层级：single-frame static DoA-Doppler。
- 主要用途：SF static 主估计器。
- 调用者：static dev、regression、transition。

### `estimatorDoaDopplerMlePilotMfOpt.m`

- 模型层级：multi-frame dynamic DoA-Doppler。
- 主要用途：MF CP-K / CP-U 主估计器。
- 调用者：dynamic dev、replay、regression。

### `evalDoaDopplerDynProfileLike.m`

- 模型层级：MF dynamic profile likelihood evaluator。
- 主要用途：正式 dynamic objective 评估。
- 调用者：MF estimator 与 probe helper。

## compatibility / 旧入口

### `estimatorDoaDopplerMlePilotMfStatOpt.m`

- 定位：MF static wrapper / compatibility。
- 用途：兼容旧入口或固定层级调用。
- 原则：如仍被引用，保留 wrapper；不要复制正式主逻辑。

### `estimatorDoaDopplerMlePilotDynOpt.m`

- 定位：dynamic wrapper / compatibility。
- 用途：兼容旧动态入口。
- 原则：行为保持，尽量复用正式 helper。

### `estimatorDoaDopplerMusic.m`

- 定位：DoA-Doppler MUSIC 类方法。
- 用途：对照或旧实验。
- 当前状态：非当前论文主线。

## 子目录职责

### `helper/`

正式 estimator helper。

负责：

- model build；
- init；
- bounds；
- branch solve；
- profile helper；
- selection；
- warm-anchor。

### `doaGrid/`

DoA / ground / ECI grid 生成和邻域。

### `covariance/`

covariance 侧 estimator 工具。

### `legacy/`

旧实现和兼容参考。

- 不作为当前主线默认实现；
- 需要迁移逻辑时，优先抽正式 helper，而不是长期维护镜像分支。

## 修改放置规则

新增 estimator 侧逻辑前，先查 `estimator/helper/README.md` 中已有 model / init / bounds / branch / warm-anchor helper；不要在主入口或 compatibility wrapper 里复制正式逻辑。若只是 replay/probe 诊断，放到 `test/dev/` 或 `test/common/`，不要下沉到 estimator。

### 新的正式 MF model 字段构造

放到：

- `estimator/helper/buildDoaDopplerMfModel.m`；
- 或新建职责明确的 helper。

### 新的 MF init / seed family

放到：

- `estimator/helper/buildDoaDopplerMfInit.m`；
- 或与 warm-start 直接相关的 helper。

### known/unknown branch solve

放到：

- `estimator/helper/solveDoaDopplerMfBranches.m`；
- `estimator/helper/runDoaDopplerMfDoaBasinEntry.m`；
- `estimator/helper/runDoaDopplerMfUnknownWarmAnchor.m`。

其中 `solveDoaDopplerMfBranches.m` 只负责 branch orchestration；truth-free DoA basin-entry acquisition 与 compact polish 已独立到 `runDoaDopplerMfDoaBasinEntry.m`，以便区分 one-shot MLE core 与 estimator route enhancement。该 acquisition 使用同一个 objective 与 fd/fdRate 范围，只临时扩大 DoA 盒捕获盆地，并在最终出口回到 compact local polish。默认只作用于 single-sat MF local solve；multi-sat MS-MF 默认旁路，除非显式把 `modelOpt.doaBasinEntryScope` 设为 `all`。外部 center / offset family 仍只能通过 truth-free orchestration 输入，不允许传入 truth DoA。`modelOpt.disableKnownEmbedded`、`modelOpt.disableDoaBasinEntry` 与 `modelOpt.disableUnknownWarmAnchor` 仅用于显式 core-only replay / scan 口径，默认不改变正式 estimator 行为。不要把 replay 的 truth probe 或 path probe 下沉到 estimator。

### candidate preference / final adoption

放到：

- `estimator/helper/preferDoaDopplerMfSolveResult.m`；
- `estimator/helper/useDoaDopplerMfUnknownWarmAnchorResult.m`。

### 只为 dev/replay 服务的诊断

不要放进 estimator 主路径。

应放到：

- `test/dev/replay/`；
- `test/dev/probe/`；
- `test/common/report/`；
- `test/common/summary/`。

## 行为保持规则

修改 estimator 时默认保持：

- 函数签名；
- 默认参数；
- 初始化策略；
- 搜索边界；
- reference-sat 语义；
- sat 顺序；
- 输出字段；
- summary 口径；
- fallback 顺序。

如果必须改变行为，需要说明：

- 哪个模块变了；
- 哪个数值路径变了；
- 风险点在哪里；
- 为什么必要；
- 应跑哪些 regression / replay。
