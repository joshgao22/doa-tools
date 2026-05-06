# scanMfKnownUnknownInformationLoss 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `paper-facing` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfKnownUnknownInformationLoss_20260428-151356.mat` |
| 当前一句话结论 | unknown `fdRate` 作为 nuisance parameter 时，对 DoA CRB 的影响几乎为零，主要信息损失集中在参考星 `fdRef`；该损失随联合帧数增加下降，多星相对单星更不敏感。 |
| 论文图定位 | `main figure / mechanism figure`，用于 unknown-rate nuisance information-loss 机制图。 |
| 决策影响 | 固定为当前代表性理论 scan；不进入 regression；后续 estimator RMSE-vs-CRB 需要另做 controlled MC。 |
| 下一步动作 | 论文图优先画 `fdRef rollback (%)` vs `P`；若需要实际 MLE 贴 CRB，另跑 resolved/local consistency scan。 |
| 禁止误用 | 不能把该纯 CRB / EFIM scan 解释成实际 `CP-U` estimator 已贴近 CRB，也不能用于判断 tooth selection、same-tooth rescue 或 full-flow 性能。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfKnownUnknownInformationLoss.m`
- 结果文档：`test/dev/scan/results/scanMfKnownUnknownInformationLoss.md`
- scan 类型：`paper-facing curve / CRB-EFIM mechanism scan`
- 主要问题：known `fdRate` 与 unknown `fdRate` nuisance 两种条件下，主参数 DoA 与参考星 `fdRef` 的 CRB / EFIM 信息损失如何变化。
- 扫描对象：`SNR`、联合帧数 `P`、帧间隔 `T_f`、`single / multi`、`known / unknown fdRate`。
- 不覆盖范围：不运行 estimator MC；不验证 MLE 收敛；不覆盖 subset tooth selection、full-flow、same-tooth polish、rescue gate 或 bad-basin 行为。
- truth 使用口径：不使用 truth 选择候选；只在 CRB 构造中使用理论真值点作为 Fisher information 的线性化点。
- 是否 paper-facing：Yes。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `K` | known `fdRate` 条件下的 CRB。 | No | 理想已知 Doppler-rate 条件。 | 不是实际 estimator 的 `CP-K` MC RMSE。 |
| `U` | unknown `fdRate` 作为 nuisance，经 EFIM 消元后的 CRB。 | No | 描述 unknown-rate 对主参数信息的理论削弱。 | 不是实际 estimator 的 `CP-U` MC RMSE。 |
| `rollback (%)` | `100 * (CRB_std_U / CRB_std_K - 1)`。 | No | 相对 known-rate 的标准差回退，适合比较信息损失比例。 | 不能当作绝对误差，也不能和 full-flow outlier rate 混用。 |
| `single / multi` | 单星参考链路 CRB 与多星联合 CRB。 | No | 比较多星几何对 unknown-rate loss 的缓解作用。 | 不代表当前 full-flow 多星 estimator 一定优于单星。 |
| `timeOriginClass` | 多帧 offset 相对于 reference time 的类别。 | No | `rightBiasedRef` 内的 P-trend 可比；`centralRef` 与 `rightBiasedRef` 不应混成同一单调曲线。 | 不要把奇偶 P 混画后解释锯齿。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：N/A，本 scan 没有 Monte Carlo repeat。
- `resolved-sample`：N/A，本 scan 是理论 CRB / EFIM。
- `outlier rate`：N/A，本 scan 不运行 estimator。
- `truth-tooth / oracle range`：N/A。
- `stress-test`：否。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfKnownUnknownInformationLoss_20260428-151356.mat` | 2026-04-28 | `representative` | `contextSeed=253`；`SNR=[-5,0,5,10] dB`；`P=[8,10,12,14,16,18,20]`；`T_f=[1,1.3333,2] ms`；主切片 `SNR=10 dB`、`T_f=1/750 s`、`rightBiasedRef`。 | unknown `fdRate` 对 DoA CRB 几乎无影响，主要损失集中在 `fdRef`；multi-sat 的 `fdRef` rollback 约为 single-sat 的一半。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `contextSeed = 253`
- `numRepeat = N/A`，理论 CRB / EFIM scan 无 MC repeat。
- `snrDbList = [-5, 0, 5, 10]`
- `frameCountList = [8, 10, 12, 14, 16, 18, 20]`
- `frameIntvlSecList = [1/1000, 1/750, 1/500]`
- 主切片：`primaryPlotSnrDb = 10`，`primaryFrameIntvlSec = 1/750`
- 关键 scan 轴：`P`、`T_f`、`SNR`、`mode=single/multi`、`fdRateMode=known/unknown`
- 关键 resolved / outlier 判据：N/A；本 scan 不运行 estimator。
- checkpoint：该 snapshot 未记录 checkpoint summary。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录 `elapsedSec`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`scanName`、`runKey`、`utcRun`、`config`、`crbSummaryTable`、`lossTable`、`fimDiagTable`、`fimDiagSummaryTable`、`primarySliceTable`、`mainSliceTable`、`mainTimeOriginClass`、`timeOriginSummaryTable`、`primaryDisplayTable`、`tfSensitivityDisplayTable`、`snrSensitivityDisplayTable`、`crbBundleSummaryCell`、`plotData`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：full FIM 的 expected ill-conditioning 在 scan 结果中以 compact condition summary 记录；主参数 interest FIM 未 near-singular。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

主切片固定 `SNR=10 dB`、`T_f=1/750 s`，并只使用统一的 `rightBiasedRef` time-origin class。

| mode | P | window (ms) | offset mean | DoA rollback (%) | fdRef rollback (%) | trace info loss (%) | min-eig loss (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| multi | 8 | 9.3333 | 0.5 | 0.001739 | 6.1017 | 10.048 | 11.171 |
| multi | 10 | 12.0000 | 0.5 | 0.004372 | 3.8559 | 6.7916 | 7.2876 |
| multi | 12 | 14.6667 | 0.5 | 0.009210 | 2.6599 | 4.8647 | 5.1149 |
| multi | 14 | 17.3333 | 0.5 | 0.017219 | 1.9467 | 3.6432 | 3.7826 |
| multi | 16 | 20.0000 | 0.5 | 0.029545 | 1.4870 | 2.8247 | 2.9090 |
| multi | 18 | 22.6667 | 0.5 | 0.047500 | 1.1733 | 2.2514 | 2.3059 |
| multi | 20 | 25.3333 | 0.5 | 0.072567 | 0.94968 | 1.8352 | 1.8726 |
| single | 8 | 9.3333 | 0.5 | ~0 | 11.8700 | 0.007531 | 20.096 |
| single | 10 | 12.0000 | 0.5 | 0 | 7.5727 | 0.007997 | 13.584 |
| single | 12 | 14.6667 | 0.5 | ~0 | 5.2514 | 0.008272 | 9.7299 |
| single | 14 | 17.3333 | 0.5 | ~0 | 3.8554 | 0.008445 | 7.2867 |
| single | 16 | 20.0000 | 0.5 | 0 | 2.9505 | 0.008560 | 5.6498 |
| single | 18 | 22.6667 | 0.5 | 0 | 2.3307 | 0.008639 | 4.5033 |
| single | 20 | 25.3333 | 0.5 | 0 | 1.8875 | 0.008695 | 3.6707 |

### 5.2 按扫描轴汇总

| axis value | case | metric 1 | metric 2 | metric 3 | 解释 |
|---:|---|---:|---:|---:|---|
| `P=8 -> 20` | multi | `fdRef rollback: 6.10% -> 0.95%` | `DoA rollback: 0.0017% -> 0.0726%` | `trace loss: 10.05% -> 1.84%` | 增加帧数会显著降低 unknown-rate 对 `fdRef` 的相对损失。 |
| `P=8 -> 20` | single | `fdRef rollback: 11.87% -> 1.89%` | `DoA rollback: ~0` | `min-eig loss: 20.10% -> 3.67%` | single-sat 的 `fdRef` rollback 始终高于 multi-sat。 |
| `SNR=-5..10 dB` | U/K ratio | 近似不变 | N/A | N/A | rollback 是相对 CRB 比值，主要由几何和时间窗口决定，对 SNR 不敏感。 |
| `T_f=[1,1.3333,2] ms` | primary sensitivity | 诊断用途 | N/A | N/A | 可用于检查窗口时长影响；论文主图建议先固定 `T_f=1/750 s`。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `fdRef CRB rollback versus frame count` | `P` | `100*(U/K-1)` for `fdRef` CRB std | `single`、`multi` | Yes | 主图候选；固定 `SNR=10 dB`、`T_f=1/750 s`、`rightBiasedRef`。 |
| `fdRef CRB rollback versus SNR` | `SNR` | `fdRef rollback (%)` | `single`、`multi` | No / appendix | 主要说明 U/K rollback 对 SNR 不敏感。 |
| `Frame-interval sensitivity` | `P` 或 `windowMs` | `fdRef rollback (%)` | 不同 `T_f` | Appendix / diagnostic | 若要解释窗口时长，建议横轴改成 `windowMs`，不要和默认 P-trend 混用。 |
| `DoA rollback table` | `P` | `DoA rollback (%)` | `single`、`multi` | No | DoA rollback 太小，适合文字说明，不必单独成主图。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- DoA rollback 几乎为零：single-sat 基本是数值零，multi-sat 在主切片中也只有 `0.0017% -> 0.0726%`。
- `fdRef` 是主要信息损失指标：multi-sat 从 `6.10%` 降到 `0.95%`，single-sat 从 `11.87%` 降到 `1.89%`。
- multi-sat 比 single-sat 更不敏感：同一 `P` 下 multi-sat 的 `fdRef` rollback 始终低于 single-sat，例如 `P=10` 时为 `3.86%` vs `7.57%`。
- SNR 扩展不是主要信息来源：由于 U/K ratio 近似抵消噪声缩放，继续扩大 SNR 范围预计主要得到近似水平线。
- interest FIM 未 near-singular：即使 multi unknown 的 full FIM near-singular，主参数 interest block 仍可用。

### 6.2 反向、污染或未解决现象

- 该 scan 不能验证实际 estimator 是否贴近 CRB；`CP-K / CP-U` 的 RMSE-vs-CRB 需要专门的 controlled MC。
- full FIM near-singular 是 nuisance 全参数化病态，不应当作主参数 EFIM 不可用。
- 奇偶帧混画会引入 time-origin class 锯齿；主趋势只在统一 `rightBiasedRef` 内解释。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `mode=multi, fdRate=unknown` | expected full-FIM ill-conditioning | `full near-singular = 84/84`，但 `interest near-singular = 0` | 说明应压制 full-FIM warning，但保留 interest FIM condition summary。 |
| odd/even mixed `P` | time-origin mixing | 奇偶 P 可能出现锯齿 | 说明主图必须固定 time-origin class。 |

## 7. 机制解释

### 7.1 当前解释

unknown `fdRate` 作为 nuisance parameter 后，`fdRef` 与 `fdRate` 在跨帧相位轨迹中存在直接耦合，因此 EFIM 消元会削弱参考时刻 `fdRef` 的有效信息。随着 `P` 增加，跨帧相位轨迹提供更多关于斜率和截距的联合约束，`fdRate` nuisance 被消元后的额外损失快速下降。

DoA 主参数在当前几何与参考星参数化下，与 `fdRate` nuisance 的 EFIM 耦合较弱，因此 DoA rollback 远小于 `fdRef` rollback。多星几何进一步提供额外约束，使 `fdRef` 对 unknown-rate 的相对损失低于 single-sat。

### 7.2 这个 scan 支持什么

- 支持 known/unknown-rate 信息损失主要体现在 `fdRef`，不是 DoA。
- 支持 `fdRef` rollback 随 `P` 增加快速下降。
- 支持多星协作不仅降低绝对 CRB，也能降低 unknown-rate nuisance 对 `fdRef` 的相对回退。
- 支持论文中把 unknown `gamma_ref` 写成 nuisance parameter，并用 EFIM loss 解释退化机制。

### 7.3 这个 scan 不证明什么

- 不证明实际 estimator 默认路径已修复。
- 不证明 `CP-U` Monte Carlo RMSE 一定贴近 EFIM。
- 不证明 full-flow tooth selection、same-tooth rescue 或 low-complexity flow 已通过。
- 不适合写 regression 契约；这是理论曲线结果，不是 pass/fail guardrail。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 不触碰 estimator 数值路径。 |
| flow 默认路径 | 不改；与 subset / rescue flow 无关。 |
| replay 下一步 | 若要闭环论文图，下一步需要 controlled MLE-vs-CRB local consistency scan / replay。 |
| regression | 不写；只作为结果文档和论文图候选。 |
| 论文图 | `main / mechanism figure`，优先画 `fdRef rollback (%)` vs `P`。 |
| 排障记录 | 主记录保留“known/unknown-rate 信息损失依赖 CRB/EFIM 口径”的结论即可，不复制长表。 |

## 9. 限制与禁止解释

- 不要把 CRB rollback 当作实际 RMSE rollback。
- 不要把该结果用于证明 default `CP-U` estimator 没有 tail。
- 不要把 full-FIM near-singular 解读成主参数 EFIM 不可用。
- 不要混合 `centralRef` 与 `rightBiasedRef` 后解释 P-trend。
- 不要把 SNR 曲线平坦解释成 noise 不重要；这里比较的是 known/unknown 的相对比值。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfKnownUnknownInformationLoss_20260428-151356.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfKnownUnknownInformationLoss.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：当前脚本的默认 `saveSnapshot` 可能已经与该历史 snapshot 不同，复现实验时应以本节 snapshot 配置为准。

## 11. 历史备注

- 当前只绑定 `scanMfKnownUnknownInformationLoss_20260428-151356.mat` 作为代表性 snapshot。
- 旧版若混入奇数 `P`，会引入 `centralRef / rightBiasedRef` time-origin 锯齿；默认主图不再混画该口径。
