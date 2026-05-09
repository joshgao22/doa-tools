# performance 目录说明

`performance/` 放 CRB、FIM、EFIM 和相关 Jacobian。它服务论文理论分析层，不负责 Monte Carlo orchestration、summary 表格或实验图生产。

## 当前论文对应关系

### static single-frame 分析

对应内容：

- single-frame static DoA-Doppler；
- SS/MS static baseline；
- static CRB 对照。

主要入口：

- `crbPilotSfDoaDoppler.m`；
- `crbPilotSfDoaOnlyEffective.m`。

### multi-frame dynamic 分析

对应内容：

- continuous-phase multi-frame model；
- known Doppler-rate；
- unknown Doppler-rate；
- nuisance-rate 信息损失；
- CP/IP 对比。

主要入口：

- `crbPilotMfDoaDoppler.m`；
- `crbPilotDynDoaDoppler.m`。

### 通用 DoA CRB

对应内容：

- 通用 deterministic / stochastic DoA baseline；
- 早期 array-only 对照。

主要入口：

- `crbDetDoa.m`；
- `crbStoDoa.m`。

## 主要文件

### `crbPilotSfDoaDoppler.m`

- 模型：single-frame static DoA-Doppler。
- 作用：SF static CRB。
- 对应实验：static perf、SF regression。

### `crbPilotSfDoaOnlyEffective.m`

- 模型：single-frame DoA-only pilot-model deterministic CRB。
- 作用：把已知 pilot、几何诱导 Doppler、path gain 与噪声方差折算为 `crbDetDoa` 的 per-sat effective source power。
- 典型用途：DoA-only estimator 没有建模 Doppler，但接收 pilot 已经带 Doppler 相位时，用 pilot-model effective-gain CRB 对比 MLE；matched DoA-Doppler estimator 仍使用 `crbPilotSfDoaDoppler.m`。
- 输出口径：`effectiveGainMode="unit"` 给 array-only unit-gain 对照，`"pilotModel"` 给当前 no-Doppler DoA-only atom 的 coherent-projection 口径。

### `crbPilotMfDoaDoppler.m`

- 模型：multi-frame DoA-Doppler。
- 作用：MF CP/IP 或动态主模型 CRB 主体。
- 对应实验：dynamic perf、CP/IP 对比。

### `crbPilotDynDoaDoppler.m`

- 模型：dynamic wrapper / compatibility。
- 作用：动态 CRB 兼容入口。
- 对应实验：旧脚本或过渡入口。

### `crbPilotMfStatDoaDoppler.m`

- 模型：MF static wrapper / compatibility。
- 作用：多帧 static 兼容入口。
- 对应实验：旧脚本或过渡入口。

### `crbPilotDoaDoppler.m`

- 模型：通用 wrapper。
- 作用：旧入口兼容。
- 当前建议：不作为新主逻辑首选。

### `crbDetDoa.m` 与 `crbStoDoa.m`

- 模型：通用 DoA CRB。
- 作用：通用阵列 baseline。
- 对应实验：早期 / 通用实验。
- 注意：`crbDetDoa` 不直接接收 path gain 或 Doppler pilot；这些信号模型细节应在上层 wrapper 中折算成 `pwrSource`，例如 `crbPilotSfDoaOnlyEffective.m`。

## crbHelper

### `buildDoaDopplerJacobian.m`

- 职责：static DoA-Doppler Jacobian。

### `buildDynDoaDopplerJacobian.m`

- 职责：dynamic DoA-Doppler Jacobian。

### `buildDopplerPilot.m`

- 职责：Doppler pilot 相关构造。

### `buildDopplerPilotJacobian.m`

- 职责：Doppler pilot Jacobian。

### `doa2dirJacobian.m`

- 职责：DoA 到方向向量的 Jacobian。

### `projectCrbToAngleMetric.m`

- 职责：将 CRB 投影到 angle metric。

## 放置规则

### 放在 performance 的内容

- FIM；
- EFIM；
- CRB；
- Jacobian；
- Schur complement / nuisance 消元；
- 信息损失分析。

### 不放在 performance 的内容

Monte Carlo orchestration：

- 放 `test/dev/perf/` 或 dev 顶层 perf 脚本。

CRB summary table：

- 放 `test/common/report/`。

CRB plotData：

- 放 `test/common/report/` 或 `test/common/plot/`。

估计器 objective / profile likelihood：

- 放 `estimator/` 或 `estimator/helper/`。

卫星几何 / reference Doppler 状态：

- 放 `satellite/scene/`。

## 维护原则

### 理论层与实验层分开

- performance 只负责计算理论界；
- 实验脚本负责调用 CRB 并整理结果；
- 不在 CRB 函数中保存文件、画图、跑 task grid。

### known-rate / unknown-rate 语义要清楚

- known-rate 是理想基准；
- unknown-rate 通过 nuisance parameter 消元得到 EFIM；
- 不要把 unknown-rate 的退化只解释为“参数更多”，应保留信息损失项的结构解释。

### warning 处理

- 某些 full-FIM singular / pinv warning 可以是 expected warning；
- 局部 regression 可 suppress；
- `verbose=true` 时应能暴露 CRB/FIM 诊断。
