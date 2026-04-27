# scanMfRegimeMapByWindow

## 对应 scan

- `test/dev/scan/scanMfRegimeMapByWindow.m`

## 扫描目标

验证论文主线中的多帧观测区间是否存在：在同一窗口内，DoA drift 仍可作为固定 DoA 近似处理，但静态 Doppler 已出现明显失配，而一阶 Doppler / Doppler-rate 模型仍然足够。

该 scan 是确定性量级分析，不运行 Monte Carlo estimator，不调用 subset bank，也不改变 dynamic flow。它用于支撑 `DoA quasi-static / Doppler-dynamic` regime justification 与论文图候选。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/scan/scanMfRegimeMapByWindow_20260427-213020.mat` | `frameCountList=1:100`，`frameIntvlMsList=[0.5,0.75,1,1.3333,1.5,2,2.5,3]`，`fc=11.7 GHz`，`rho=550 km`，`v_perp=13.6 km/s`，`doaSlowTolDeg=0.1` | target regime 约为 `3.91 ms` 到 `70.53 ms`；主实验窗口 `P=10~20, T_f=1/750 s` 明确落在 target regime 内。 |

## 当前代表性结果

2026-04-27 的代表性结果扫描了 `100 x 8 = 800` 个离散窗口配置。当前 primary DoA-static tolerance 为 `0.1 deg`，target regime 判据为：固定 DoA 有效、静态 Doppler 失效、一阶 Doppler 有效。

### Regime summary by frame interval

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

该表说明 target regime 主要由窗口总长度 `P T_f` 决定，而不是由帧数或帧间隔单独决定。不同 `T_f` 下的最大 target window 都集中在约 `69–70.5 ms`。

### Model boundary summary

| metric | rule | boundary |
|---|---|---|
| DoA static validity | `DoA drift <= 0.1 deg` | 最后有效窗口约 `70.53 ms`，下一采样点约 `70.75 ms` 失效 |
| static Doppler invalid | `|fd drift| >= 50 Hz` 或静态相位误差 `>= 1 rad` | 约 `3.91 ms` 开始明显失配 |
| first-order Doppler validity | `|fd residual| <= 1 Hz` 且 phase residual `<= 0.1 rad` | 在当前 `0–130 ms` 曲线范围内没有失效 |
| target regime | 三个条件同时成立 | 约 `3.91–70.53 ms` |

旧表头 `firstWindowMs / lastWindowMs` 容易被误读。代码已改为 `firstPassWindowMs / lastPassWindowMs / firstFailAfterPassWindowMs`，明确每一行都是“该条件成立”的窗口范围。

### 主实验窗口解释

| P | `T_f` (ms) | window (ms) | DoA drift (deg) | fd drift (Hz) | static phase error (rad) | first-order residual |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1.3333 | 13.333 | 0.01889 | 174.99 | 7.3301 | `~1e-5 Hz / ~2e-7 rad` |
| 20 | 1.3333 | 26.667 | 0.03778 | 349.98 | 29.32 | `~8e-5 Hz / ~3e-6 rad` |

这两个主配置都满足 fixed-DoA、static-Doppler-invalid、first-order-Doppler-valid 三个条件。因此它们可以直接支撑论文中“DoA 在窗口内近似静态，但 Doppler 漂移已经必须建模”的主张。

## DoA tolerance 解释

`0.1 deg` 不应被写成 estimator 精度要求，而应解释为窗口内固定 DoA 模型允许的未建模 DoA drift 上限。若使用更保守的 `0.05 deg`，target regime 上界会收缩到约 `35 ms`；若使用 `0.1 deg`，上界约为 `70 ms`。两者都覆盖主实验的 `P=10~20, T_f=1/750 s`。

因此结果文档和图中应同时展示多个 DoA tolerance boundary，而不是只保留单一阈值。代码已新增 `doaSlowTolDegList` 与 `doaToleranceSummaryTable`，用于同时记录严格、保守和宽松标准下的 fixed-DoA 有效窗口。

## 可观察现象

- DoA drift 近似随 `P T_f` 线性增长。
- Doppler drift 也近似随 `P T_f` 线性增长，但静态 Doppler 诱导的累计相位误差随窗口平方增长，因此在十几毫秒窗口内已经达到数 rad。
- 一阶 Doppler residual 在当前几何近似下很小，说明该 order scan 内部支持“一阶 Doppler 足够”的建模层级。
- 当前 scan 是量级模型，不是 exact-scene 几何 residual scan；若后续需要更强证据，应另做 scene-based exact geometry scan，逐帧计算真实 `DoA(t)` 与 `fd(t)` 后拟合 fixed-DoA / constant-Doppler / first-order-Doppler residual。

## 当前结论

当前结果支持论文主线：在 Starlink-inspired `T_f=1/750 s` 且 `P=10~20` 的观测窗口内，DoA drift 约为 `0.0189–0.0378 deg`，Doppler drift 约为 `175–350 Hz`，静态 Doppler 相位误差约为 `7.3–29.3 rad`。因此静态 Doppler 模型不足，而 fixed-DoA + first-order Doppler / nuisance-rate 的模型层级是合理的。

## 对代码和论文图的影响

- scan 继续保持确定性量级分析，不引入 `numRepeat`、`snrDb` 或 estimator replay batch。
- 图应保留连续曲线，并叠加离散 `P,T_f` 扫描点；这样更容易看出主实验窗口与模型边界的相对位置。
- DoA static 图应显示多条 tolerance boundary，避免把 `0.05 deg` 或 `0.1 deg` 写成绝对算法性能标准。
- 结果可作为论文 Section 2 中 order analysis / regime justification 的图候选；不应迁移到 regression。

## 后续建议

该 scan 当前不需要继续增加 MC 维度。下一步若要增强论证，应新增或扩展 exact-scene scan，用真实轨道几何检查固定 DoA、一阶 Doppler 拟合的 residual，而不是继续膨胀当前 order scan。
