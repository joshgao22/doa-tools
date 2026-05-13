# replayMfMsMleCrbCleanTrim 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / fixed-bound local sanity` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfMsMleCrbCleanTrim_20260513-181333.mat` |
| 当前一句话结论 | CleanTrim 已可暂时固定为 truth-centered / CRB-scaled / in-tooth / core-only 的局部一致性 replay：SS 线稳定贴近 bounded local CRB，MS 在 `2×CRB` DoA hard box 内有绝对 RMSE 增益并在部分 SNR 达到 own-CRB level，但仍暴露 MS local basin / fd tail 边界。 |
| 决策影响 | 保留为 replay-level diagnostic 与 paper-facing scan 前的 sanity check；不继续在本 replay 中扩大扫描，不修改 estimator / CRB 主公式 / 默认 flow。 |
| 下一步动作 | 将更系统的 DoA box scale、fdRef half-tooth、geometry、frame count 与 information-loss 曲线交给 scan；主记录 / 机制归并版摘录本结论。 |
| 禁止误用 | 不能把 `2×CRB` hard-box 结果写成自然 local MLE 已完全贴 CRB；不能把本 replay 当作 full-flow acquisition 证明；不能据此打开 truth-aware selector 或修改 estimator 主核。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfMsMleCrbCleanTrim.m`
- 结果文档：`test/dev/replay/results/replayMfMsMleCrbCleanTrim.md`
- replay 类型：truth-centered / CRB-scaled / in-tooth / core-only 小 MC replay。
- 主要问题：在固定 hard envelope 和 in-tooth local range 下，检查 `SS/MS-SF-Static`、`SS/MS-MF-CP-K`、`SS/MS-MF-CP-U` 是否与对应 CRB 口径一致，并判断 MS-MF 是否能形成相对 SS-MF 的 local gain。
- 观察范围：`P=10`，SNR sweep `[-15,-10,-5,0,5,10] dB`，seed `253:452`，`200` repeats，`2×CRB` DoA hard box，MF `fdRef` truth-local `0.300 tooth` half-width。
- 不覆盖范围：不验证 full-flow tooth acquisition、subset selection、truth-free basin-entry、same-tooth rescue 或默认 estimator robustness；不用于证明 global / unbounded MLE 性能。
- truth 使用口径：truth 只用于构造 controlled hard envelope、in-tooth local range、CRB-scaled box 和结果评价；不进入 runtime selector / adoption / final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `SS-SF-Static` | 单星单帧 static DoA-Doppler anchor。 | Evaluation / controlled range | 固定 truth-centered static fdRef range 与 DoA hard box。 | static 地基 sanity；主要看 static fdRef 和 angle 是否回到 CRB-level。 |
| `MS-SF-Static` | 多星单帧 static DoA-Doppler anchor。 | Evaluation / controlled range | 相同 local range 下加入 second sat。 | 检查 static MS 是否形成绝对增益。 |
| `SS-MF-CP-K` | 单星多帧 continuous-phase known-rate baseline。 | Evaluation / controlled range | `fdRate` known，不估 `gamma_ref`。 | 已知 rate 的 MF local CRB sanity。 |
| `SS-MF-CP-U` | 单星多帧 continuous-phase unknown-rate baseline。 | Evaluation / controlled range | `gamma_ref` 作为 nuisance 估计。 | unknown-rate 信息损失与 fdRate tail 的单星对照。 |
| `MS-MF-CP-K` | 多星多帧 continuous-phase known-rate baseline。 | Evaluation / controlled range | 多星几何进入 local MLE，rate known。 | 判断 MS angle CRB 更紧时 estimator 是否仍能跟上。 |
| `MS-MF-CP-U` | 多星多帧 continuous-phase unknown-rate 主诊断行。 | Evaluation / controlled range | 多星 + unknown `gamma_ref` nuisance。 | 当前最关键行；看 MS 相对 SS 的绝对增益、own-CRB gap、fdRef/fdRate tail。 |
| `joint-trim` | health + `angle <= 5×CRB` + `fdRef <= 5×CRB` 的统计子集。 | Evaluation only | 只改统计口径。 | 用来剥离少数 fdRef / fdRate tail；不能替代 full-sample RMSE。 |
| `CRB path audit` | K/U shared block、Schur loss、monotonicity、reference-state invariant 的路径审计。 | Evaluation only | 只改诊断表。 | 用于确认 K/U CRB 接近是否为路径错误；通过时不应改 CRB 主公式。 |
| `2×CRB DoA hard box` | truth-centered DoA hard envelope。 | Controlled / oracle range | 限制 DoA local optimizer 可搜索范围。 | 强 local sanity guard；不能解释为自然 basin-entry 能力。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfMsMleCrbCleanTrim_20260513-181333.mat` | 2026-05-13 | representative | `numRepeat=200`，seed `253:452`，SNR `[-15,-10,-5,0,5,10] dB`，`P=10`，`fdRefHalfWidth=0.300 tooth`，DoA hard box `2×CRB`，core-only。 | SS 线稳定；MS 在 `2×CRB` box 内有绝对 angle gain，`-5/0 dB` 达到 own-CRB level，其余 SNR 为 `beats-ss-only`；K/U CRB path 已确认正常。 | 当前代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = [-15, -10, -5, 0, 5, 10]`
- `baseSeed = 253`
- `numRepeat = 200`
- seed range：`253:452`
- frame count：`10`
- method list：`SS-SF-Static, MS-SF-Static, SS-MF-CP-K, SS-MF-CP-U, MS-MF-CP-K, MS-MF-CP-U`
- `fdRef requested half-width = 0.300 tooth`
- `static fdRef CRB floor scale = 3.00 × CRB`
- `fdRate truth-local half-width = 1000 Hz/s`
- `DoA CRB half-width scale = 2.00 × CRB`
- estimation chain：`truth-crb-scaled-static-seed-fixed-bound`
- MF objective / route：`pure-mle / core-only`
- trim definition：`health + angle <= 5×CRB + fdRef <= 5×CRB`
- checkpoint：`completed 1200/1200`，成功后清理 checkpoint artifacts。
- snapshot 保存变量：`replayData`
- 运行时间：约 `24 min 40 s`。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| task count | `1200` | `6 SNR × 200 repeats`。 |
| representative snapshot | `replayMfMsMleCrbCleanTrim_20260513-181333.mat` | 当前固定结果。 |
| SS-MF-CP-U angle MSE/CRB | `0.945–0.960` | bounded local 口径下 SS-MF 稳定贴近 CRB；低于 1 不解释为 estimator 超 CRB。 |
| MS-MF-CP-U angle MSE/CRB | `0.942–1.877` | MS 有绝对 RMSE 增益，但 own-CRB 更紧，低/高 SNR 仍为 `beats-ss-only`。 |
| MS-MF-CP-U fdRef MSE/CRB | `0.983–1.124` | joint-trim 后 fdRef 基本回到 CRB-level。 |
| MS-MF-CP-U keep rate | `0.94–1.00` | 低 SNR / `-5 dB` 有少量 fdRef / fdRate tail。 |
| CRB path audit | `path-ok` | K/U monotonicity、Schur loss、shared block、reference-state invariant 通过；不改 CRB 主公式。 |
| outlier rows | `26` | 主要来自 `MS-MF-CP-U` 的 `fdRate`、`fdRef` 或 `fdRate+boundary`。 |

### 4.3 关键对比表

#### MS-MF-CP-U vs SS-MF-CP-U, joint-trim

| SNR (dB) | SS angle MSE/CRB | MS angle MSE/CRB | SS fdRef MSE/CRB | MS fdRef MSE/CRB | SS angle RMSE (deg) | MS angle RMSE (deg) | SS keep | MS keep | MS abs gain | target class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| -15 | 0.945 | 1.877 | 0.973 | 1.026 | 0.012243 | 0.011375 | 0.865 | 0.950 | 0.071 | `beats-ss-only` |
| -10 | 0.968 | 1.126 | 0.976 | 1.027 | 0.006966 | 0.004956 | 0.985 | 0.980 | 0.289 | `beats-ss-only` |
| -5 | 0.962 | 0.942 | 1.034 | 1.034 | 0.003906 | 0.002549 | 1.000 | 0.940 | 0.347 | `crb-level` |
| 0 | 0.961 | 0.983 | 1.034 | 0.983 | 0.002195 | 0.001464 | 1.000 | 1.000 | 0.333 | `crb-level` |
| 5 | 0.961 | 1.244 | 1.035 | 1.110 | 0.001234 | 0.000926 | 1.000 | 1.000 | 0.250 | `beats-ss-only` |
| 10 | 0.959 | 1.666 | 1.035 | 1.124 | 0.000694 | 0.000603 | 1.000 | 1.000 | 0.131 | `beats-ss-only` |

#### Estimate-vs-CRB representative points

| Method | SNR (dB) | angle RMSE (deg) | angle CRB (deg) | fdRef RMSE (Hz) | fdRef CRB (Hz) | 备注 |
|---|---:|---:|---:|---:|---:|---|
| `SS-MF-CP-U` | -15 | 0.012243 | 0.012593 | 0.30637 | 0.31054 | SS-MF local sanity。 |
| `MS-MF-CP-U` | -15 | 0.011375 | 0.0083034 | 0.30362 | 0.29981 | MS 绝对 RMSE 好于 SS，但 own-CRB gap 明显。 |
| `SS-MF-CP-U` | -5 | 0.0039061 | 0.0039822 | 0.099872 | 0.098201 | SS-MF CRB-level。 |
| `MS-MF-CP-U` | -5 | 0.0025491 | 0.0026258 | 0.096412 | 0.094808 | MS-MF 达到 own-CRB level。 |
| `SS-MF-CP-U` | 0 | 0.0021947 | 0.0022393 | 0.05615 | 0.055223 | SS-MF CRB-level。 |
| `MS-MF-CP-U` | 0 | 0.0014637 | 0.0014766 | 0.052852 | 0.053315 | MS-MF CRB-level。 |
| `SS-MF-CP-U` | 10 | 0.00069361 | 0.00070814 | 0.017768 | 0.017463 | SS-MF CRB-level。 |
| `MS-MF-CP-U` | 10 | 0.00060277 | 0.00046694 | 0.017871 | 0.01686 | MS 绝对 RMSE 好于 SS，但 own-CRB gap 存在。 |

#### K/U CRB path audit

| satMode | rows | monotonic / loss / shared / invariant | angle CRB U/K median | fdRef CRB U/K median | loss trace ratio | gamma coupling | audit class |
|---|---:|---|---:|---:|---:|---:|---|
| `single` | 1200 | all pass | 1.000 | 1.0757 | `7.9973e-05` | 0.00894 | `path-ok` |
| `multi` | 1200 | all pass | 1.000 | 1.0386 | 0.067916 | 0.26061 | `path-ok` |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `SS-MF-CP-K/U` 在所有 SNR 上稳定为 bounded local CRB-level，说明 SS-MF core path 可作为 sanity anchor。
- `MS-MF-CP-U` 在所有 SNR 上都有相对 `SS-MF-CP-U` 的绝对 angle RMSE 增益；但按 MS 自己更紧的 CRB 归一化，只有 `-5 dB` 和 `0 dB` 明确达到 `crb-level`。
- `MS-MF-CP-K/U` 的 `fdRef` 在 joint-trim 后基本回到 `~1×CRB`，说明 raw fdRef tail 不是主 CRB 分母错误。
- K/U CRB path audit 已确认 `CRB_U >= CRB_K` 单调性、shared block、Schur loss 与 reference-sat Doppler / rate invariant 均通过；K/U angle CRB 接近不应解释为 CRB path 漏接。
- DoA search range 没有 truth-box violation；MF fdRef range 相对 CRB 极宽，不再支持“MF fdRef 范围低于 CRB”这一解释。

### 5.2 仍未解决或反向的现象

- MS-MF 的 own-CRB gap 在低 SNR 和高 SNR 仍存在，尤其 `10 dB` 的 `MS-MF-CP-U` angle MSE/CRB 约 `1.666`。
- 之前 `3×CRB` DoA hard-box 运行显示 MS 会明显远离 own-CRB，说明当前 `2×CRB` 是强 oracle guard，不能代表自然 basin-entry。
- outlier 仍主要集中在 `MS-MF-CP-U`，原因多为 `fdRate`、`fdRef` 或 `fdRate+boundary`，不是 CRB path failure。
- 当前结果中 static `SS-SF-Static` 的 search-range audit 曾出现 `minFdRefHalfOverCrb≈1` 与少量 `fdRefBelowCrb` 标记；由于 static fdRef RMSE/CRB 已约为 `1`，该问题不影响主结论，但若要进入论文图，需要进一步检查 bound-hit / margin 口径。
- scalar angle CRB 中 K/U 差异很小；unknown-rate information loss 应通过 Schur loss、fdRef inflation、主方向 loss 或专门 scan 展示，而不是只看 spherical angle CRB ratio。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 330 / -10 dB | `trim-fdRef` | angle 到 `2×CRB` 附近，`fdRefErrOverCrb≈16`。 | 表明 joint-trim 主要剔除 fdRef tail。 |
| 450 / -15 dB | `fdRef` tooth tail | `fdRefAbsErr≈750 Hz`，`fdRefErrOverCrb≈2500`。 | 表明 raw tail 仍包含 tooth-scale fdRef miss。 |
| 347 / -5 dB | `fdRate+boundary` | `fdRateAbsErr≈1000 Hz/s`，hit rate boundary。 | 表明 unknown-rate tail 与 fdRate range 边界有关。 |
| 378 / -15 dB | `fdRef` tooth tail | angle 很小但 `fdRefAbsErr≈750 Hz`。 | 表明 angle health 和 fdRef health 必须分开统计。 |

## 6. 机制解释

### 6.1 当前解释

CleanTrim 现在的定位已经比较清楚：它不是 full-flow estimator，而是一个 truth-centered、CRB-scaled、in-tooth、core-only 的 local sanity replay。`2×CRB` DoA hard box 把 DoA 搜索限制在极强的 oracle-local 区域内，因此 SS-MF 的 angle MSE/CRB 略低于 1 是 bounded-local 统计口径，不应解释为 estimator 超越 CRB。

在这个 controlled local 条件下，SS-MF-CP-K/U 已经稳定，说明 single-sat MF continuous-phase core path、fdRef local range 和 known/unknown-route glue 基本健康。MS-MF 的绝对 RMSE 仍普遍好于 SS-MF，但 MS CRB 更紧，导致 MS own-CRB gap 在低 / 高 SNR 更明显。这说明当前 MS 问题不是“多星没有增益”，而是“多星 local objective 更容易暴露 DoA / fdRef / fdRate coupling 与 tail”。

K/U 的 angle CRB 几乎相同并不说明没有 information loss。CRB audit 显示 `gamma_ref` coupling 和 Schur loss 存在，且 multi-sat coupling 更强；只是这些 loss 在当前几何、窗口和 spherical angle scalar 上投影较弱。更适合用 scan 展示的是 `Δ_gamma`、主方向 loss、fdRef bound inflation、frame count / geometry / time aperture 对 loss 的影响。

### 6.2 这个结果支持什么

- 支持把 `replayMfMsMleCrbCleanTrim` 暂时固定为 local CRB sanity replay。
- 支持不再继续在 CleanTrim 里堆更多路线或 candidate bank；更多参数面交给 scan。
- 支持 `SS-MF-CP-K/U` 作为当前 MF local sanity anchor。
- 支持 `MS-MF` 的 paper-facing 描述采用“有绝对 RMSE gain，但 own-CRB gap / tail 仍需报告”的口径。
- 支持 K/U CRB path 暂不修改主公式，information-loss 展示应转向 dedicated scan。

### 6.3 这个结果不证明什么

- 不证明 full-flow acquisition 或 subset tooth selection 已解决。
- 不证明 MS-MF 自然 local optimizer 在无 oracle hard-box 下稳定贴 CRB。
- 不证明 `2×CRB` hard-box 可作为论文主 estimator 默认设置。
- 不证明 unknown-rate 对所有 scalar 指标都有明显退化。
- 不证明可以把 truth-centered range、truth label 或 oracle bound 放进 runtime selector / gate / adoption。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；CleanTrim 只证明 controlled local core sanity。 |
| flow 默认路径 | 不改；不把 CleanTrim 结果推进到 subset / adoption / rescue。 |
| CRB 主公式 | 不改；CRB path audit 为 `path-ok`，K/U angle CRB 接近需要解释，不是直接 bug。 |
| regression | 不写新的 regression；这是 fixed-bound local replay，不是 pass/fail 契约。 |
| replay / scan 下一步 | CleanTrim 暂停继续扩展；新增或使用 scan 做 DoA box scale、fdRef half-width、frame count、geometry、information-loss 曲线。 |
| 论文图 / 论文口径 | 可作为 paper-facing scan 前的 sanity / appendix evidence；正文主图应使用 scan，并标注 local / resolved / oracle-bound 口径。 |
| 排障记录 | 主记录摘一句“CleanTrim 固定为 local sanity，更多参数面交给 scan”；机制归并版保留 K/U CRB path 与 2×/3× hard-box sensitivity 结论。 |

## 8. 限制与禁止解释

- 不要把 `2×CRB` DoA hard box 下的结果解释为自然 basin-entry 能力。
- 不要把 `SS-MF` 的 angle MSE/CRB 略低于 1 写成 estimator 超 CRB；这是 bounded local / finite MC / hard-envelope 口径。
- 不要用 CleanTrim 替代 `scanMfMsMleCrbInToothConsistency` 或 information-loss scan。
- 不要因为 K/U scalar angle CRB 接近就说没有 information loss；应直接看 Schur loss、fdRef inflation 和主方向 loss。
- 不要把 truth-centered `fdRef/fdRate` range、truth DoA 或 hard seed label 进入 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 outlier table 中的 fdRef/fdRate tail 当成 CRB 主公式错误；应先区分 fdRef tooth tail、fdRate boundary 和 local basin tail。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfMsMleCrbCleanTrim_20260513-181333.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfMsMleCrbCleanTrim.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- `SS-MF-CP-U-Robust` 已在 fixed-bound CleanTrim 中降级：hard-envelope intersection 后与 core-only 基本一致，默认无价值。
- static fdRef hidden floor 曾导致 `fdRef` MSE/CRB 偏低；当前 replay 继续保留 static fdRef range / floor 审计，但 paper-facing 图前仍应增加 bound-hit / margin 检查。
- `3×CRB` DoA hard-box 版本显示 MS-MF 明显远离 own-CRB，因此当前 `2×CRB` 版本只能作为 oracle-local sanity，不应扩展解释。
- K/U CRB path audit 已把“MS K/U CRB 太接近所以公式错”的路线降级；information-loss 需要由 dedicated scan 展示。
