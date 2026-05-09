# replaySfMsDoaCrbDiagnose 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / CRB-scale audit complete` |
| 最新代表性 snapshot | `test/data/cache/replay/replaySfMsDoaCrbDiagnose_20260509-155057.mat` |
| 当前一句话结论 | `MS-SF-DoA` 的 unit-gain CRB gap 主要来自 DoA-only CRB 没有按 second-sat effective pilot gain / path gain 缩放；amp-aware CRB 后 `Other-SF-DoA` 与 `MS-SF-DoA` 均回到 CRB-level 附近。 |
| 决策影响 | 不继续调 `MS-SF-DoA` estimator；下一步应把 replay-local amp-aware deterministic CRB 口径 formalize 到 `performance/` 或 scan wrapper，再回到 `scanSfStaticMleCrbConsistency` 验证四种 static case。 |
| 下一步动作 | 暂停本 replay 机制扩展；用本 snapshot 作为 CRB signal-scale 诊断代表，后续转向正式 CRB helper / `scanSfStaticMleCrbConsistency` 的 unit-vs-amp 对照字段。 |
| 禁止误用 | 不能把 unit-gain CRB 下的 `MS-SF-DoA` gap 写成 estimator failure；不能把 replay-local amp-aware CRB 直接当成已正式修改的 paper-facing CRB；不能把 alpha sweep 下沉为默认权重策略。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replaySfMsDoaCrbDiagnose.m`
- 结果文档：`test/dev/replay/results/replaySfMsDoaCrbDiagnose.md`
- replay 类型：小 MC / static anchor CRB 口径诊断 / post-probe representative run。
- 主要问题：解释 `scanSfStaticMleCrbConsistency.m` 中 `MS-SF-DoA` 相对 DoA-only CRB 稳定约 `1.16 x` 的 gap，到底来自 MLE / solver / fusion 路径，还是来自 CRB / metric / signal-scale 口径。
- 观察范围：single-frame DoA-only，global sat `1, 2`，`SS-SF-DoA`、`Other-SF-DoA`、`MS-SF-DoA` 和 `MS-SF-DoA-TruthInit`，同时输出 unit-gain deterministic CRB、single-frame DoA-Doppler CRB 与 replay-local pilot-model amp-aware CRB。
- 不覆盖范围：不覆盖 single-frame static DoA-Doppler；不覆盖 multi-frame dynamic；不验证 full-flow tooth acquisition；不改变 estimator 默认路径；不提供论文最终 MC 曲线。
- truth 使用口径：truth 只用于 truth-init solver probe、noise-free sanity、truth/final objective probe 和 amp-aware CRB 的 clean-signal effective gain 估计；不进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `SS-SF-DoA` | 参考星单星、单帧、DoA-only MLE。 | No | 作为 ref-only DoA baseline。 | 用于观察公共 DoA-only MLE / CRB gap，以及 multi-sat 是否比 ref-only 有增益。 |
| `Other-SF-DoA` | 第二星单独、单帧、DoA-only MLE。 | No | 作为 other-sat 单独信息量审计对象。 | 若 other 单星相对 unit CRB 大幅偏离，但相对 amp-aware CRB 接近 1，说明 unit CRB 对 second-sat 信息量过乐观。 |
| `MS-SF-DoA` | 两星联合、单帧、DoA-only MLE。 | No | 联合两星 DoA-only objective。 | 用于判断 multi-sat DoA-only 是否实现 CRB 预期的 joint gain。 |
| `MS-SF-DoA-TruthInit` | 与 `MS-SF-DoA` 相同，但初始 DoA seed 设为 truth；仍走完整 DoA 搜索域。 | Oracle only | 只改变初值，不收紧搜索盒。 | 若与 baseline 一致，则排除 init / local basin 作为主因。 |
| unit-gain deterministic CRB | 原 DoA-only deterministic CRB 口径，未按 replay 中估计到的 per-sat effective pilot gain 缩放。 | No | 只改评价 denominator。 | 适合作为旧 scan 口径对照；本 replay 显示它对 second sat 信息量过乐观。 |
| amp-aware CRB | replay-local pilot-model CRB 诊断，用 clean signal 在 truth DoA 处估计 effective pilot gain，并按增益平方缩放 FIM。 | Oracle / diagnostic | 只改评价 denominator，不改 estimator。 | 若 MLE/CRB 回到 1 附近，说明 gap 主要来自 signal-scale 口径；当前还不是正式 `performance/` helper。 |
| `empirical joint gain` | 用实测 `SS-SF-DoA` 与 `Other-SF-DoA` 方差预测 joint gain。 | Evaluation only | 只改离线解释。 | 用于验证 MS 实际 gain 是否符合 “other sat 实际较弱” 的经验模型。 |
| `mean-noisy Hessian` | 固定 probe seed，在 noisy data 上重复估计数值 Hessian 后取均值。 | Diagnostic only | 只做局部曲率审计。 | 用来辅助判断 objective curvature 与 analytic FIM 是否同阶；不作为 paper-facing CRB 证明。 |
| `alphaSat2` sweep | 在 probe seeds 上改变 second-sat objective 权重。 | Diagnostic only | 只做 probe-level objective weighting。 | 只用于确认 second sat 确实进入 objective；不能据此设默认权重。 |
| `noise-free sanity` | 用无噪声信号跑 SS / Other / MS DoA-only estimator。 | Oracle / diagnostic | 只改输入噪声。 | 三条 clean case 均回 truth 时，可排除 steering、sat-order、local/global DoA 映射硬错。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replaySfMsDoaCrbDiagnose_20260509-155057.mat` | 2026-05-09 | `representative` | `numRepeat=500`，seed `253:752`，SNR `-20:5:10 dB`，global sat `1,2`，checkpoint / progressbar 启用，snapshot 保存 `replayData`。 | unit-gain CRB 下 `Other-SF-DoA` 约 `2.92-3.40 x CRB`、`MS-SF-DoA` 约 `1.15-1.16 x CRB`；amp-aware 后分别回到约 `1.02-1.18 x` 与 `1.06-1.07 x`，证明 gap 主要是 second-sat signal-scale CRB 口径。 | 当前代表性结果；覆盖此前 200-repeat probe。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDbList = [-20, -15, -10, -5, 0, 5, 10]`
- `baseSeed = 253`
- `numRepeat = 500`
- seed range：`253:752`
- task 数：`3500 = 7 SNR x 500 seeds`
- global selected satellites：`1, 2`
- `probeSeedCountPerSnr = 4`
- `hessianMeanRepeat = 8`
- alpha sweep：`[0, 0.1, 0.25, 0.5, 1, 2, 5]`，只用于 probe seeds。
- truth-init：完整 DoA grid，不使用局部 truth box。
- checkpoint：从 `126/3500` resume，最终 `3500/3500` 完成，运行后清理 `tmp/replaySfMsDoaCrbDiagnose/snrm20to10_seed253to752_rep500_35a49c6e`。
- snapshot 保存变量：`replayData`
- 运行时间：约 `7 min 30 s`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `SS-SF-DoA` CRB-local RMSE / CRB | `1.0817-1.0908` | ref-only DoA-only MLE 本身有约 `8-9%` 公共 finite-sample / metric gap。`-20 dB` full raw 被 `0.4%` outlier 拉到 `4.3827 x`，不作为主口径。 |
| `Other-SF-DoA` / unit CRB | `2.9188-3.4014` | second sat 相对 unit-gain CRB 严重偏离，说明旧 CRB 对其信息量过强。 |
| `Other-SF-DoA` / amp-aware CRB | `1.0165-1.1846` | 按 effective pilot gain 缩放后，second sat 单独 DoA 结果回到 CRB-level 附近。 |
| `MS-SF-DoA` / unit joint CRB | `1.1532-1.1626` | 旧口径下看似有稳定 `~16%` gap。 |
| `MS-SF-DoA` / amp-aware joint CRB | `1.0577-1.0663` | amp-aware 后 MS joint 结果接近 CRB-level，说明 estimator / fusion 不是主因。 |
| analytic expected joint gain | `1.0784` | unit-gain CRB 预期 ref-only 到 joint 应有 `7.84%` 标准差增益。 |
| amp-aware expected joint gain | `1.0111` | 考虑实际 second-sat effective gain 后，预期 joint 增益只有约 `1.1%`。 |
| realized joint gain | `1.0034-1.0135` | 实际 MS 增益与 amp-aware / empirical 预期一致。 |
| other effective pilot gain | `~0.348` | second sat 有效增益明显低于 ref sat 的 `~0.978`。 |
| other amp / unit FIM trace ratio | `0.12129` | second-sat FIM 约为 unit-gain 假设下的 `12.1%`，解释 unit CRB 过乐观。 |
| truth-init diff | median diff `0`，max diff `~1.2e-6 deg` | MS gap 不是初值 / basin 问题。 |
| noise-free sanity | SS / Other / MS clean angle error 均为 `0` | 排除 steering、sat-order、local/global DoA 映射硬错。 |

### 4.3 关键对比表

#### Unit CRB vs amp-aware CRB

| SNR (dB) | Other / unit CRB | Other / amp CRB | MS / unit CRB | MS / amp CRB | expected gain amp | realized gain |
|---:|---:|---:|---:|---:|---:|---:|
| -20 | 3.4014 | 1.1846 | 1.1606 | 1.0645 | 1.0111 | 1.0135 |
| -15 | 3.0325 | 1.0561 | 1.1557 | 1.0599 | 1.0111 | 1.0135 |
| -10 | 2.9537 | 1.0287 | 1.1535 | 1.0580 | 1.0111 | 1.0133 |
| -5 | 2.9300 | 1.0204 | 1.1532 | 1.0577 | 1.0111 | 1.0125 |
| 0 | 2.9223 | 1.0177 | 1.1552 | 1.0595 | 1.0111 | 1.0102 |
| 5 | 2.9197 | 1.0168 | 1.1593 | 1.0633 | 1.0111 | 1.0064 |
| 10 | 2.9188 | 1.0165 | 1.1626 | 1.0663 | 1.0111 | 1.0034 |

#### Empirical joint-gain prediction

| SNR (dB) | expected analytic gain | expected empirical gain | realized gain | 解释 |
|---:|---:|---:|---:|---|
| -20 | 1.0784 | 1.0069 | 1.0135 | 实测 other 方差很大，经验模型只预期约 `0.7%` gain；MS 实际略高于经验预测。 |
| -15 | 1.0784 | 1.0067 | 1.0135 | 同上。 |
| -10 | 1.0784 | 1.0066 | 1.0133 | 同上。 |
| -5 | 1.0784 | 1.0065 | 1.0125 | 同上。 |
| 0 | 1.0784 | 1.0065 | 1.0102 | 同上。 |
| 5 | 1.0784 | 1.0066 | 1.0064 | 基本完全一致。 |
| 10 | 1.0784 | 1.0066 | 1.0034 | 高 SNR 下实际 gain 略低，但仍接近经验 / amp-aware 预期。 |

#### Other-sat covariance audit

| SNR (dB) | unit eig ratio 1 | unit eig ratio 2 | amp eig ratio 1 | amp eig ratio 2 | 解释 |
|---:|---:|---:|---:|---:|---|
| -20 | 16.481 | 10.942 | 1.999 | 1.327 | 低 SNR 仍有 threshold / tail，amp-aware 后显著收敛但不完全贴。 |
| -15 | 8.650 | 9.239 | 1.049 | 1.121 | amp-aware 后两个 tangent 方向均接近 1。 |
| -10 | 8.633 | 8.720 | 1.047 | 1.058 | amp-aware 后接近 1。 |
| -5 | 8.661 | 8.561 | 1.050 | 1.038 | amp-aware 后接近 1。 |
| 0 | 8.688 | 8.507 | 1.054 | 1.032 | amp-aware 后接近 1。 |
| 5 | 8.706 | 8.488 | 1.056 | 1.030 | amp-aware 后接近 1。 |
| 10 | 8.718 | 8.480 | 1.057 | 1.029 | amp-aware 后接近 1。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `MS-SF-DoA-TruthInit` 与 baseline 完全重合，`medianEstimateDiffDeg=0`，`maxEstimateDiffDeg≈1.2e-6 deg`，说明当前 gap 不是 init / local basin / solver 起点问题。
- `Other-SF-DoA` 在 unit-gain CRB 下约 `3 x CRB`，但在 amp-aware CRB 下除 `-20 dB` 外约 `1.02-1.06 x CRB`，说明 second-sat 大 gap 主要由 signal-scale denominator 引起。
- `MS-SF-DoA` 在 amp-aware joint CRB 下稳定约 `1.06 x CRB`，接近 ref-only 公共 DoA-only gap，说明 multi-sat fusion 没有明显额外失效。
- `normalizationAuditTable` 中 `effectivePilotGainAbs` 显示 ref sat 约 `0.978`、other sat 约 `0.348`，且 other `ampToDetFimTraceRatio=0.12129`，与 second-sat 信息量被 unit CRB 高估约 `1/0.121` 的现象一致。
- `empiricalGainTable` 中 expected empirical gain 约 `1.0065-1.0069`，与 realized gain `1.0034-1.0135` 基本吻合；这说明 MS 实际结果符合“other sat 实际较弱”的经验模型。
- noise-free sanity 中 SS / Other / MS 三条路径均回 truth，说明 sat-order、steering、global/local DoA 映射没有硬错误。
- checkpoint / progressbar 外壳稳定：本轮从已有 checkpoint 恢复，`3500/3500` task 完成后清理临时目录并保存 snapshot。

### 5.2 仍未解决或反向的现象

- `SS-SF-DoA` 在 `-20 dB` 的 full raw RMSE 被极少量 outlier 拉到 `4.3827 x CRB`，但 CRB-local RMSE 为 `1.0908 x`。因此低 SNR 下解释 gain 时应使用 CRB-local / empirical 口径，不应直接读 `ampAwareCrbAuditTable` 中可能被 full SS outlier 污染的 `realizedGainVsRef`。
- `Other-SF-DoA` 在 `-20 dB` 即使使用 amp-aware CRB，仍有 `otherFullRmseOverAmpCrb=1.1846`，且 covariance amp eig ratio 可到 `1.999 / 1.327`。这更像低 SNR threshold / finite-sample tail，不是 CRB scale 主因。
- `MS-SF-DoA` amp-aware 后仍约 `1.06 x CRB`，并非严格 `1.00 x`。这可作为 finite-sample / metric / profile MLE 小 gap 记录，不能过度宣称“完全贴合”。
- 数值 Hessian 在高 SNR 下的 trace ratio 不稳定，尤其 other / joint 的 noisy Hessian trace 会受 profile noise variance 或 objective scale 影响；Hessian probe 只作为辅助 scale 诊断，不作为正式 CRB 证明。
- alpha sweep 只覆盖 probe seeds，不能用来确定 default sat2 weight。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `SNR=-20, SS-SF-DoA full` | low-SNR outlier caveat | full raw `SS-SF-DoA` 为 `4.3827 x CRB`，但 CRB-local 为 `1.0908 x`，outlier rate `0.004`。 | 说明 `-20 dB` full RMSE 不适合用于 gain denominator；应使用 CRB-local / empirical 表。 |
| `SNR=-15..10, Other-SF-DoA` | signal-scale mismatch | unit eig ratio 约 `8.5-9.2`，amp eig ratio 约 `1.03-1.12`。 | 直接证明 unit CRB 的 second-sat FIM 约高估一个固定尺度。 |
| `SNR=-10..10, MS-SF-DoA` | amp-aware consistency | `MS / amp CRB` 约 `1.058-1.066`。 | 说明 MS estimator / fusion 在正确 signal-scale CRB 下接近 CRB-level。 |
| clean SS / Other / MS | model sanity | 三者 `angleErrDeg=0`。 | 排除 steering / sat-order / local-global DoA mapping 硬错。 |
| truth-init pair | solver sanity | baseline 与 truth-init RMSE 相同，最大估计差约 `1.2e-6 deg`。 | 排除初值和 local basin 作为主因。 |

## 6. 机制解释

### 6.1 当前解释

这轮结果把 `MS-SF-DoA` 的稳定 gap 从 estimator / fusion 问题，收缩到 deterministic DoA-only CRB 的 signal-scale 口径问题。旧 unit-gain CRB 假设 second sat 的有效信号尺度与其 analytic geometry FIM 匹配，但 replay 中 second sat 的 effective pilot gain 只有约 `0.348`，对应 FIM trace 缩放约 `0.121`。因此 unit CRB 把 second sat 看得过强，导致 `Other-SF-DoA` 看似有 `~3 x CRB` gap，也导致 joint CRB 预期的 `1.0784 x` gain 过高。

amp-aware CRB 用 clean pilot-model effective gain 对 per-sat FIM 缩放后，`Other-SF-DoA` 的 covariance eigen ratio 大多回到 `1.03-1.12`，`MS-SF-DoA` 的 RMSE/CRB 回到 `1.06 x` 附近。此时 realized joint gain 与 empirical / amp-aware expected gain 基本一致，说明 MS 路径不是“没有吃到 second sat”，而是 second sat 在当前信号模型中本来就只提供很小的有效额外信息。

这也解释了为什么 `MS-SF-DoA` 不应继续通过 alpha sweep 或 solver 修改来“救”。alpha sweep 能证明 second sat 的 objective contribution 存在，但不能改变 second sat 信号幅度较弱这一事实。真正应推进的是把 deterministic CRB 的 per-sat amplitude / effective-SNR scaling 正式化，并在 scan 中同时报告 unit / amp 两套口径，直到 paper-facing CRB 定稿。

### 6.2 这个结果支持什么

- 支持停止继续优化 `MS-SF-DoA` estimator 或 truth-init / local-box probe。
- 支持将 unit-gain DoA-only CRB 标为旧对照口径，而不是最终 paper-facing denominator。
- 支持 formalize amp-aware / pilot-model deterministic CRB：FIM 需要按每星 effective pilot gain、path gain 或接收 SNR 缩放。
- 支持回到 `scanSfStaticMleCrbConsistency`，增加 `rmseOverUnitCrb` 与 `rmseOverAmpCrb` compact 字段，验证 `SS-SF-DoA / MS-SF-DoA / SS-SF-Static / MS-SF-Static` 四条 static 曲线。
- 支持在论文中将 CRB 解释为 deterministic / conditional CRB，条件于实际链路幅度或等效接收 SNR，而不是隐含所有卫星 unit gain。

### 6.3 这个结果不证明什么

- 不证明正式 `performance/` CRB helper 已经修改；amp-aware CRB 当前仍是 replay-local 诊断口径。
- 不证明 stochastic CRB 是正确模型；当前论文主线仍应使用 deterministic / conditional CRB。
- 不证明 dynamic MS-MF 已贴 CRB；本 replay 是 single-frame DoA-only static anchor 诊断。
- 不证明 alpha sweep 应进入 estimator 默认权重。
- 不证明 full-flow acquisition / tooth selection / Doppler-aided runtime 条件已经通过。
- 不证明 `-20 dB` full-sample low-SNR tail 可以忽略；如果作为论文图，仍需报告 full / resolved / outlier rate。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。`TruthInit`、noise-free sanity 与 alpha sweep 都没有指向 solver / fusion 主路径问题。 |
| CRB / performance | 下一步应把 amp-aware deterministic CRB 口径 formalize 到 `performance/` 或 wrapper；先保留 unit / amp 两套字段做对照，不直接覆盖旧结果。 |
| flow 默认路径 | 不改。本 replay 不涉及 dynamic flow、subset selection 或 same-tooth rescue。 |
| regression | 暂不新增 regression。CRB 口径尚未从 replay-local helper 正式下沉，不应写 pass/fail 契约。 |
| replay 下一步 | 本 replay 冻结为 representative diagnostic；只有当 CRB helper 或 signal model 改动后才重跑。 |
| scan 下一步 | 修改 / 扩展 `scanSfStaticMleCrbConsistency`，输出 unit vs amp CRB 对照，确认四种 static case 的 paper-facing 口径。 |
| 论文图 / 论文口径 | 本结果只支持 static anchor CRB denominator 修正；论文主图仍应来自 scan，并同时标注 full / resolved / outlier rate。 |
| 排障记录 | 可在主记录 / 机制归并版补一句：`MS-SF-DoA` unit-CRB gap 已定位为 second-sat signal-scale CRB 口径差，不再继续调 DoA-only MLE。 |

## 8. 限制与禁止解释

- 不要把 `Other-SF-DoA ~3 x unit CRB` 写成 other-sat estimator 失败；应先说明 unit CRB 未按 effective gain 缩放。
- 不要把 replay-local amp-aware CRB 直接当成正式 paper-facing CRB；它仍需迁移到 `performance/` 或 scan wrapper 并统一命名。
- 不要用 `MS-SF-DoA` 的 unit-CRB gap 继续驱动 solver / estimator 主核修改。
- 不要把 `alphaSat2=2` 或 `alphaSat2=5` 的 probe 行解释成默认权重建议；probe seed 数太少，且该方向绕开了 CRB 口径主因。
- 不要把 clean sanity 写成 noisy MLE 一定贴 CRB；它只排除硬模型 / 映射错误。
- 不要把 `-20 dB` full raw outlier 与 CRB-local / amp-aware 结论混用；低 SNR 要同时报告 outlier / CRB-local rate。
- 不要把本 single-frame DoA-only 结论直接外推到 `SS/MS-SF-Static` 的 DoA-Doppler CRB 或 dynamic CP-K/U CRB。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replaySfMsDoaCrbDiagnose_20260509-155057.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
test/dev/replay/replaySfMsDoaCrbDiagnose.m
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。完整 per-seed / Hessian / alpha sweep 表保存在 `replayData` 中；命令行默认只打印 compact preview。

## 10. 历史备注

- 早期 200-repeat probe 中，`MS-SF-DoA-TruthInit` 与 baseline 已经完全一致，初步排除了 init / basin 问题。
- 后续新增 `Other-SF-DoA`、per-sat Hessian、normalization audit、empirical gain 与 amp-aware CRB 后，问题从 “MS 没吃到 second sat” 改判为 “second sat effective signal scale 被 unit-gain CRB 高估”。
- 当前 500-repeat representative snapshot 固化了该结论，并验证 checkpoint / progressbar / snapshot 外壳正常。后续不再继续向本 replay 加机制 probe。
