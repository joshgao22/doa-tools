# estimator/helper 目录说明

`estimator/helper/` 放正式 estimator helper。这里的代码属于主实现，不是 dev-only glue。新增 helper 应直接反映功能边界，避免 `core`、`misc`、`tmp`、`util2` 这类含糊命名。

## model / init / bounds

### `buildDoaDopplerSfModel.m`

- 职责：SF static model 构造。
- 内容：reference-sat、satWeight、DoA/fd 参数域。

### `buildDoaDopplerMfModel.m`

- 职责：MF dynamic model 构造。
- 内容：CP/IP、known/unknown rate、timeOffset、reference state。

### `buildDoaDopplerSfInit.m`

- 职责：SF static 初始化。
- 内容：fdRef init、DoA anchor、static branch。

### `buildDoaDopplerMfInit.m`

- 职责：MF dynamic 初始化。
- 内容：fd line、fdRate seed、warm-start set。

### `buildDoaDopplerMfBounds.m`

- 职责：MF 参数 bounds。
- 内容：DoA/fd/fdRate 局部盒与默认范围。

## objective / profile helper

### `evalDoaDopplerSfProfileLike.m`

- 职责：SF profile likelihood evaluator。
- 定位：SF objective 主核。

### `evaluateDoaDopplerMfObjective.m`

- 职责：MF objective 包装评估。
- 用途：分支、selection、probe 中统一调用正式 objective。

### `buildPilotAtomSfKernel.m`

- 职责：SF block kernel helper。
- 作用：保持 SF/static block 语义一致。

## optimization / branch solve

### `runDoaDopplerMfOptimization.m`

- 职责：MF 主优化包装。
- 作用：统一调用优化器。

### `runDoaDopplerMfKnownEmbeddedOptimization.m`

- 职责：known-rate 嵌入式优化。
- 作用：known-rate 与 unknown layout 对齐。

### `solveDoaDopplerMfBranches.m`

- 职责：MF known/unknown 分支求解。
- 当前连续相位 MF local solve 在存在 `initDoaParam / initDoaHalfWidth` 时，会先运行 truth-free DoA basin-entry acquisition：用较宽 DoA 盒捕获好盆地，再回到原 compact local box 做 final polish；该逻辑只改变 DoA entry，不改变 fdRef / fdRate 范围、reference-sat 语义或 objective。
- CP-U warm-anchor 会锚定在当前最佳 DoA basin 后再释放 fdRate，避免 unknown-rate 分支重复继承旧 static basin。
- 注意：这是 branch orchestration，不应塞 dev-only 诊断或大体量 trace；truth / path probe 仍只能放在 replay。

### `runDoaDopplerMfUnknownWarmAnchor.m`

- 职责：unknown warm-anchor release。
- 内容：mainWarmAnchor / fixed-DoA warm-anchor。

## candidate preference / adoption

### `useDoaDopplerMfUnknownWarmAnchorResult.m`

- 职责：warm-anchor adoption。
- 内容：控制 solveVariant 和最终出口。

### `preferDoaDopplerMfSolveResult.m`

- 职责：candidate preference。
- 内容：objective / branch preference 规则。

### `shouldUseDoaDopplerMfWarmAnchorParfor.m`

- 职责：warm-anchor 内层 parfor gate。
- 规则：默认必须串行，显式 opt-in 才可能开启。
- 原因：内层 parfor 可能因 worker 侧微小优化差异改变敏感 winner。

## debug / evalDiag

### `buildDoaDopplerMfDebug.m`

- 职责：紧凑 debug 组织。
- 规则：只能组织旁路诊断，不改变 objective。

### `buildDoaDopplerMfEvalDiag.m`

- 职责：evalDiag 组织。
- 规则：只能组织旁路诊断，不改变 residual、参数语义或 branch selection。

## refinement helper

### `refineFdPeakFromGrid.m`

- 职责：fd grid peak refinement。
- 用途：局部峰值细化。

### `refineGridEstimate.m`

- 职责：通用网格估计 refinement。

### `gridIndexToDoaParam.m`

- 职责：grid index 到 DoA 参数转换。

## 并行化边界

### 允许优先考虑 parfor 的地方

- test/perf 脚本层；
- Monte Carlo 外层；
- subset candidate 独立评估；
- 多起点评估；
- comb/tooth/grid 扫描。

### 不默认开启 parfor 的地方

- estimator 主求解内核；
- warm-anchor 内层 release seed；
- nested parfor；
- 会改变 winner 的敏感 branch。

### warm-anchor parfor 特别规则

- estimator 默认路径保持串行；
- DoA basin-entry acquisition 也是串行候选求解；不要在 estimator 内层默认打开 parfor 或 nested parfor；
- dev/perf 可显式 opt-in；
- opt-in 后必须复核 wrong-tooth guard 和 quick regression；
- 敏感性用 `replayMfWarmAnchorParforSensitivity.m` 观察。

## 什么不应放在这里

### dev-only replay / probe

放到：

- `test/dev/replay/`；
- `test/dev/probe/`。

### summary / table / plotData

放到：

- `test/common/report/`；
- `test/common/summary/`；
- `test/common/plot/`。

### fixture / context 构造

放到：

- `test/common/fixture/`。

### scene / reference-sat / Doppler geometry

放到：

- `satellite/scene/`。

### CRB / FIM / EFIM

放到：

- `performance/`；
- `performance/crbHelper/`。

## helper 拆分原则

- 主入口只做 orchestration；
- 重复出现且职责独立的逻辑抽正式 helper；
- 不为了形式拆一串一次性 helper；
- 不保留“主入口 + 巨型 core/helper”的伪拆分；
- static 与 dynamic 共用逻辑优先抽共享 helper；
- debug/probe 默认旁路，不污染主 objective 或 residual。
