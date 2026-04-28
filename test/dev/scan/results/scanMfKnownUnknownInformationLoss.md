# scanMfKnownUnknownInformationLoss 结果记录

## 对应 scan

- `test/dev/scan/scanMfKnownUnknownInformationLoss.m`

## 扫描目标

扫描 known / unknown Doppler-rate 条件下的多帧 CRB / EFIM 信息损失，量化把参考星 Doppler rate 从已知量改为 nuisance parameter 后，对主参数 DoA 与参考星 `fdRef` 的理论性能界回退。

该 scan 只比较理论界：

- `K`：known `fdRate`，对应理想已知 Doppler-rate 条件；
- `U`：unknown `fdRate`，对应 `fdRate` 作为 nuisance parameter 后用 EFIM 消元；
- `U/K rollback (%) = 100 * (CRB_std_U / CRB_std_K - 1)`。

它不是实际 estimator Monte Carlo 性能 scan，不验证 tooth selection、subset flow、same-tooth polish 或优化 basin。实际 `CP-K / CP-U` RMSE 是否贴近 CRB / EFIM，需要另由 performance / controlled replay 验证。

## Snapshot index

| snapshot | 状态 | 配置 | 结论 |
|---|---|---|---|
| `test/data/cache/scan/scanMfKnownUnknownInformationLoss_20260428-151356.mat` | current | `contextSeed=253`，`SNR=[-5,0,5,10] dB`，`P=[8,10,12,14,16,18,20]`，`T_f=[1,1.3333,2] ms`，主切片 `SNR=10 dB`、`T_f=1/750 s`，`timeOriginClass=rightBiasedRef`，`84` 个 CRB grid case，只保存轻量 `scanData` | 当前代表性结果。unknown `fdRate` 对 DoA CRB 几乎没有影响，主要损失集中在 `fdRef`；`fdRef` U/K 回退随帧数增加单调下降，multi-sat 的回退明显小于 single-sat；interest FIM 未出现 near-singular，multi unknown 的 full FIM ill-conditioning 属于 nuisance 全参数化的预期病态。 |

## 当前代表性结果

当前代表性结果只保留偶数帧数，并在主表中使用同一个 reference-time class：`rightBiasedRef`。该口径与仓库默认 `P=10` 的周期窗口一致，避免把奇数帧 `centralRef` 和偶数帧 `rightBiasedRef` 混画成一条 P-trend 曲线。

### 主切片：SNR = 10 dB，T_f = 1/750 s

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

## 可观察现象

### 1. unknown `fdRate` 对 DoA CRB 几乎没有负面影响

主切片中，single-sat 的 DoA rollback 基本为数值零；multi-sat 的 DoA rollback 也只有 `0.0017% -> 0.0726%`，远小于 `fdRef` rollback。因此，这个 scan 不支持“unknown rate 会显著损害 DoA 理论界”的说法。

更准确的说法是：在当前参考星参数化和观测窗口内，`fdRate` nuisance 与 DoA 主参数的 EFIM 耦合较弱；其主要影响不在 DoA，而在参考星 `fdRef`。

### 2. `fdRef` 是主要信息损失指标

`fdRef` rollback 随帧数增加单调下降：

- multi-sat：`6.10% -> 0.95%`；
- single-sat：`11.87% -> 1.89%`。

这说明 unknown Doppler-rate 的主要代价是削弱参考时刻 Doppler 的理论分辨能力。随着帧数增加，跨帧相位轨迹提供更多关于 `fdRef / fdRate` 的联合约束，`fdRate` nuisance 被 EFIM 消元后的额外损失快速下降。

### 3. multi-sat 比 single-sat 对 unknown-rate 更不敏感

同样的 `P` 和 `T_f` 下，multi-sat 的 `fdRef` rollback 始终低于 single-sat。例如：

- `P=8`：single `11.87%`，multi `6.10%`；
- `P=10`：single `7.57%`，multi `3.86%`；
- `P=20`：single `1.89%`，multi `0.95%`。

这支持一个对论文有用的解释：多星几何信息不仅能改善绝对 CRB，也能降低 unknown `fdRate` 对 `fdRef` 主参数的相对信息损失。

### 4. SNR 扩展不是该 scan 的主要信息来源

该 scan 关注的是 U/K ratio 或 rollback percentage。known / unknown FIM 都会随噪声方差近似整体缩放，因此相对回退对 SNR 不敏感。当前 `SNR=[-5,0,5,10] dB` 已足够说明该点；继续扩展到更宽 SNR 范围，预计主要得到近似水平线，不会显著增加机制信息。

若论文需要 RMSE vs SNR 曲线，应在 `scanMfCpIpPerfMap` 或专门 controlled MC 中比较实际 estimator 的 `CP-K / CP-U` RMSE，并叠加 CRB / EFIM，而不是把这个纯理论 information-loss scan 扩展成性能 MC。

### 5. 帧数趋势应只在同一 time-origin class 内解读

旧版 `P=8:20` 结果出现奇偶锯齿，是因为奇数帧的 offset 对称于 reference frame，属于 `centralRef`；偶数帧采用 `rightBiasedRef`，offset mean 为 `0.5`。两类时间原点下 `fdRef` 与 `fdRate` 的耦合不同，不能混成一条单调 P 曲线。

当前代表性结果只使用偶数 P，主趋势图使用统一 `rightBiasedRef`，因此 `fdRefLossPct` 随 P 单调下降，结果解释是自洽的。若后续想研究 reference-time placement 本身，应单独新增或扩展一个 sensitivity 图，而不是并入默认主图。

### 6. full FIM ill-conditioning 不影响当前主结论

FIM condition summary 为：

| mode | fdRate mode | cases | full near-singular | interest near-singular | min rcond full | median rcond full | min rcond interest | median rcond interest |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| single | known | 84 | 0 | 0 | 6.189e-05 | 3.4009e-04 | 3.4658e-04 | 1.9071e-03 |
| single | unknown | 84 | 0 | 0 | 6.2053e-11 | 1.9413e-09 | 2.7682e-04 | 1.7681e-03 |
| multi | known | 84 | 0 | 0 | 2.4225e-10 | 2.5047e-10 | 2.5540e-10 | 2.5540e-10 |
| multi | unknown | 84 | 84 | 0 | 4.8726e-16 | 2.8584e-15 | 2.5217e-10 | 2.5488e-10 |

`multi + unknown` 的 full FIM 全部 near-singular，但 interest FIM 的 near-singular 计数为 0。因此，本 scan 中 suppress full-FIM warning 是合理的：full nuisance 参数化病态不应刷屏；主参数 EFIM / interest block 的条件数仍在 summary 中保留。

## 图形解释

当前 scan 的图应按以下方式解读：

1. **`fdRef CRB rollback versus frame count`**：主图。横轴为 `P`，纵轴为 `100*(U/K-1)`。它回答“unknown `fdRate` 相对 known `fdRate` 会让 `fdRef` 理论标准差增加多少”。该图应作为本文 nuisance-rate information-loss 机制图的优先候选。
2. **`fdRef CRB rollback versus SNR`**：诊断图。若不同 SNR 曲线近似水平，说明 U/K rollback 主要由几何和时间窗口决定，而不是由 SNR 决定。该图可用于解释为什么不需要继续扩大 SNR 范围。
3. **`Frame-interval sensitivity of fdRef rollback`**：诊断图。比较同一 `P` 下不同 `T_f` 的相对回退。若曲线接近，说明在当前设置中 `P` / time-origin placement 比帧间隔本身更主导相对损失；若后续要深入研究窗口时长，可将横轴改成 `windowMs` 另做专门 sensitivity scan。

默认论文图建议只保留第 1 张；第 2 / 3 张放到补充结果或作为内部 sanity check。

## 当前结论

当前结果没有明显数值问题。它支持以下论文表述：

> 在连续相位多帧模型中，把参考星 Doppler rate 作为 unknown nuisance parameter 主要带来参考星 `fdRef` 的 EFIM 信息损失，而对 DoA 主参数的理论性能影响很小。该 `fdRef` 损失随联合帧数增加快速下降，并且多星协作相较单星能进一步降低 unknown-rate 的相对回退。

同时，该结果也给出边界：

- 这不是实际 estimator 性能图；
- 不能用于证明 `CP-U` 实际 RMSE 一定优于或接近 `CP-K`；
- 不能替代 flow / tooth-selection / same-tooth refine 的 replay 或 scan；
- 不应把 full FIM near-singular 解读成主参数 EFIM 不可用，因为 interest FIM 当前没有 near-singular。

## 对 replay / regression / 论文图的影响

- 该 scan 可作为论文 Section 4 / Section 5 中 unknown-rate nuisance information-loss 的理论机制图候选。
- 主文若只放一张图，建议画 multi / single 的 `fdRef rollback (%)` vs `P`，固定 `SNR=10 dB`、`T_f=1/750 s`、`rightBiasedRef`。
- DoA rollback 可以不单独成图，只在文字或表格中说明其量级接近 0。
- SNR 不建议继续扩大；若需要 SNR 曲线，应转到 estimator MC performance scan。
- 帧数不建议混入奇数 P 默认主图；若要讨论 reference-time placement，则另开小节或新增专门 sensitivity 结果。
- 该结果不应迁移为 regression。它是理论 scan 结果，不是 pass/fail 契约。

## 后续建议

1. 保留当前 scan 默认配置：偶数 `P`、主 `T_f=1/750 s`、主 `SNR=10 dB`。
2. 不继续扩大 SNR 范围；当前 `[-5,0,5,10] dB` 已足以说明相对 rollback 对 SNR 不敏感。
3. 若需要更完整的论文性能闭环，下一步应跑 controlled MC：比较实际 `CP-K / CP-U` RMSE 与 CRB / EFIM 的关系。
4. 若要解释旧版奇偶锯齿，可单独写 reference-time placement sensitivity，不要混入该 scan 的默认主趋势。

## 历史 / superseded snapshots

- 当前未保留 superseded snapshot 记录；已删除的旧 cache 文件不再在本文档中索引。
