# scanMfRegimeMapByWindow 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfRegimeMapByWindow_20260427-213020.mat` |
| 当前一句话结论 | 在 Starlink-inspired `T_f=1/750 s`、`P=10~20` 窗口内，DoA drift 仍小，而 Doppler drift 与静态 Doppler 相位误差已不可忽略；fixed-DoA + first-order Doppler / nuisance-rate 的模型层级成立。 |
| 论文图定位 | `main figure / regime justification`，服务论文的 DoA quasi-static / Doppler-dynamic 区间识别。 |
| 决策影响 | 固定为当前代表性 order-analysis scan；不进入 regression；后续若需要更强证据，可补 exact-scene geometry scan。 |
| 下一步动作 | 论文图展示窗口边界与 `P=10~20, T_f=1/750 s` 的相对位置；可选补 scene-based exact geometry residual。 |
| 禁止误用 | `0.1 deg` 是窗口内 DoA drift 容许上限，不是 estimator 成功阈值；该 scan 不验证 estimator、CRB 或 flow。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfRegimeMapByWindow.m`
- 结果文档：`test/dev/scan/results/scanMfRegimeMapByWindow.md`
- scan 类型：`paper-facing curve / deterministic order-analysis scan`
- 主要问题：是否存在一个多帧窗口，使 DoA 可以近似静态，但 Doppler 漂移和静态 Doppler 相位误差已经必须建模，且一阶 Doppler 仍足够。
- 扫描对象：`P`、`T_f`、窗口长度、DoA drift、Doppler drift、静态 Doppler 相位误差、一阶 Doppler residual。
- 不覆盖范围：不运行 Monte Carlo estimator；不计算 CRB；不调用 subset bank；不验证 full-flow 或 rescue。
- truth 使用口径：未使用 estimator truth；这是解析量级 / deterministic model-boundary scan。
- 是否 paper-facing：Yes。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `windowMs` | 多帧观测窗口长度，近似由 `P*T_f` 决定。 | No | 横轴核心变量；target regime 主要由窗口总长度决定。 | 不要只按 `P` 或 `T_f` 单独解释。 |
| `DoA static validity` | `DoA drift <= doaSlowTolDeg`。 | No | 判断 fixed-DoA 模型是否仍可接受。 | `doaSlowTolDeg` 不是 estimator 精度门限。 |
| `static Doppler invalid` | `|fd drift| >= 50 Hz` 或静态 Doppler 相位误差 `>= 1 rad`。 | No | 判断常 Doppler / static Doppler 模型是否已经失配。 | 不代表 actual estimator failure。 |
| `first-order Doppler validity` | 一阶 Doppler residual 仍低于 `1 Hz` 与 `0.1 rad`。 | No | 支持 Doppler-rate 一阶模型足够。 | 不证明更高阶几何永远不需要。 |
| `target regime` | fixed-DoA 有效、static Doppler 失效、一阶 Doppler 有效三者同时成立。 | No | 论文主模型的量级合理区间。 | 不要写成算法通过区间。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：N/A，无 Monte Carlo repeat。
- `resolved-sample`：N/A，无 estimator 输出。
- `outlier rate`：N/A。
- `truth-tooth / oracle range`：N/A。
- `stress-test`：否。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfRegimeMapByWindow_20260427-213020.mat` | 2026-04-27 | `representative` | `frameCountList=1:100`；`frameIntvlMsList=[0.5,0.75,1,1.3333,1.5,2,2.5,3]`；`fc=11.7 GHz`；`rho=550 km`；`v_perp=13.6 km/s`；`doaSlowTolDeg=0.1`。 | target regime 约为 `3.91–70.53 ms`；主实验窗口 `P=10~20, T_f=1/750 s` 明确落在 target regime 内。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `frameCountList = 1:100`
- `frameIntvlSecList = [0.0005, 0.00075, 0.001, 1/750, 0.0015, 0.002, 0.0025, 0.003]`
- `curveWindowSecList = linspace(0, 0.13, 600)`
- `representativeFrameCountList = [5, 8, 10, 15, 20]`
- `representativeFrameIntvlSecList = [0.001, 1/750, 0.002]`
- `carrierFreqHz = 11.7e9`
- `minSlantRangeM = 550e3`
- `transverseVelocityMps = 13600`
- `doaSlowTolDeg = 0.1`
- `fdDynamicTolHz = 50`
- `staticDopplerPhaseTolRad = 1`
- `firstOrderFdResidualTolHz = 1`
- `firstOrderPhaseResidualTolRad = 0.1`
- checkpoint：N/A，deterministic light scan。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录 `elapsedSec`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`scanName`、`runKey`、`config`、`scanTable`、`aggregateTable`、`curveTable`、`targetSummaryTable`、`representativeTable`、`modelBoundaryTable`、`plotData`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：N/A；该 scan 是解析量级扫描。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

当前 snapshot 扫描了 `100 x 8 = 800` 个离散窗口配置。target regime 判据是：DoA fixed 模型有效、static Doppler 失效、一阶 Doppler 有效。

| frame interval (ms) | target P range | max target window (ms) | first DoA-static fail P | first static-Doppler invalid P | first first-order fail P |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 8–100 | 50.0 | — | 8 | — |
| 0.75 | 6–94 | 70.5 | 95 | 6 | — |
| 1.0 | 4–70 | 70.0 | 71 | 4 | — |
| 1.3333 | 3–52 | 69.333 | 53 | 3 | — |
| 1.5 | 3–47 | 70.5 | 48 | 3 | — |
| 2.0 | 2–35 | 70.0 | 36 | 2 | — |
| 2.5 | 2–28 | 70.0 | 29 | 2 | — |
| 3.0 | 2–23 | 69.0 | 24 | 2 | — |

### 5.2 按扫描轴汇总

| axis value | case | metric 1 | metric 2 | metric 3 | 解释 |
|---:|---|---:|---:|---:|---|
| `window≈3.91 ms` | lower boundary | static Doppler invalid starts | DoA still valid | first-order valid | 从这里开始常 Doppler 模型已经不够。 |
| `window≈70.53 ms` | upper boundary | DoA static last valid | static Doppler invalid | first-order valid | `0.1 deg` DoA drift 容许下的 target regime 上界。 |
| `P=10, T_f=1/750 s` | main setting | `window=13.333 ms` | `DoA drift=0.01889 deg` | `fd drift=174.99 Hz` | 明确在 target regime 内。 |
| `P=20, T_f=1/750 s` | sensitivity setting | `window=26.667 ms` | `DoA drift=0.03778 deg` | `fd drift=349.98 Hz` | 仍在 target regime 内，静态 Doppler 相位误差更大。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Regime boundary map | `windowMs` 或 `P` | DoA drift / Doppler drift / static phase / first-order residual | model-boundary curves | Yes | 主图建议标出 `P=10` 与 `P=20` 的主实验窗口。 |
| Target regime by frame interval | `T_f` | target `P` range / max target window | frame interval lines | Appendix | 用于说明 target regime 主要由窗口长度决定。 |
| Representative rows | discrete `(P,T_f)` | drift / residual values | selected rows | Yes / table | 可作为正文表或图注补充。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- target regime 主要由窗口总长度决定，不同 `T_f` 下的最大 target window 都集中在约 `69–70.5 ms`。
- 主配置 `P=10, T_f=1/750 s` 的 DoA drift 约 `0.01889 deg`，Doppler drift 约 `174.99 Hz`，静态相位误差约 `7.33 rad`。
- 扩展配置 `P=20, T_f=1/750 s` 的 DoA drift 约 `0.03778 deg`，Doppler drift 约 `349.98 Hz`，静态相位误差约 `29.32 rad`。
- 一阶 Doppler residual 在当前 `0–130 ms` 曲线范围内没有触发失效，说明该量级模型支持“一阶 Doppler 足够”的主建模层级。

### 6.2 反向、污染或未解决现象

- 当前 scan 是 order-analysis，不是 exact-scene geometry residual scan；如果审稿需要更强的几何证据，应新增 scene-based scan。
- `0.1 deg` 是模型层级容许的未建模 DoA drift，不是 estimator angle hit threshold。
- 该 scan 不说明 MLE 贴 CRB，也不说明 CP/IP 性能差异。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `window < 3.91 ms` | not target | static Doppler 还未明显失效 | 说明过短窗口没有必要引入 dynamic Doppler 模型。 |
| `window > 70.53 ms` with `0.1 deg` DoA tol | upper-bound violation | fixed-DoA 近似不再满足 primary tolerance | 说明本文主模型不是任意长窗口模型。 |
| `0–130 ms` curve range | first-order validity | 一阶 Doppler residual 未失效 | 支持在当前窗口内不升级为更高阶动态状态。 |

## 7. 机制解释

### 7.1 当前解释

LEO 几何下，短窗口内 DoA 累积变化随窗口长度近似线性增长，但在 `P=10~20, T_f=1/750 s` 时仍只有 `0.02–0.04 deg` 量级；因此把 DoA 作为参考时刻常值主参数是合理的。

与此同时，Doppler drift 也随窗口增长，但静态 Doppler 诱导的累计相位误差随窗口近似二次增长，因此十几毫秒窗口内已经达到数 rad。这个量级直接支持把 reference Doppler rate 作为 nuisance parameter 纳入连续相位多帧模型，而不是继续采用常 Doppler 或跨帧独立相位的简化模型。

### 7.2 这个 scan 支持什么

- 支持论文的 regime identification：DoA quasi-static / Doppler-dynamic 区间真实存在。
- 支持 `P=10~20, T_f=1/750 s` 作为主实验窗口。
- 支持 fixed-DoA + first-order Doppler / nuisance-rate 的模型层级。
- 支持把 static Doppler 模型作为不足的 baseline，而不是主模型。

### 7.3 这个 scan 不证明什么

- 不证明 estimator 默认路径已修复。
- 不证明 `CP-U` 或 `CP-K` 的 RMSE 接近 CRB。
- 不证明 full-flow 没有 wrong-tooth / same-tooth bad basin。
- 不证明 `doaSlowTolDeg=0.1` 是 estimator 的最终成功标准。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 只服务模型区间论证。 |
| flow 默认路径 | 不改；与 subset / rescue flow 无关。 |
| replay 下一步 | 不需要 replay；若要增强证据，可新增 exact-scene geometry scan。 |
| regression | 不写；该结果是论文 order-analysis，不是自动契约。 |
| 论文图 | `main`，用于 Section 2 / system model 的 regime justification。 |
| 排障记录 | 只需保留“paper-facing regime justification 已固定”的结论，不复制长表。 |

## 9. 限制与禁止解释

- 不要把 `0.1 deg` 写成 estimator 精度要求。
- 不要把该 deterministic scan 解释为 Monte Carlo 性能。
- 不要把 first-order valid 扩展到任意长窗口或任意几何。
- 不要把主配置落在 target regime 内解释为 low-complexity flow 已通过。
- 不要为了该 scan 增加 `numRepeat`、SNR 维度或 estimator branch；这会混淆 order-analysis 与 performance scan。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfRegimeMapByWindow_20260427-213020.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfRegimeMapByWindow.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：若当前脚本已新增多 tolerance summary，本 snapshot 仍以 `doaSlowTolDeg=0.1` 的历史代表性配置为准。

## 11. 历史备注

- 当前只绑定 `scanMfRegimeMapByWindow_20260427-213020.mat` 作为代表性 snapshot。
- 后续若使用新版脚本重跑多 tolerance summary，应追加新 snapshot，并标明是否覆盖本结果；不要在同一个 snapshot 下混写新字段。
