# replayMfSsMleCrbMetricDiagnose 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / post-fix validation` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfSsMleCrbMetricDiagnose_20260507-174134.mat` |
| 当前一句话结论 | SS-MF 的 DoA basin-entry 修复后，`SS-MF-CP-K/U` 在大部分 SNR 下稳定回到约 `1.17 x spherical CRB` 的局部 MLE 行为；低 SNR 仅剩少量 residual basin / fdRate nuisance tail。 |
| 决策影响 | SS estimator 主路径可阶段性固定；本 replay 结果用于诊断与结果解释，不作为 paper-facing 大 MC 统计。 |
| 下一步动作 | 可回到 `scanMfMleCrbInToothConsistency` 做 paper-facing / coverage scan；scan 中重点验证更大 seed 下 angle/CRB、fdRef/CRB 与 CP-U resolved rate。 |
| 禁止误用 | 不能宣称 CP-K `fdRef` 超越 CRB；不能把 5-seed replay 当成论文统计曲线；不能把 truth/oracle in-tooth 条件当成 runtime full-flow 条件。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfSsMleCrbMetricDiagnose.m`
- 结果文档：`test/dev/replay/results/replayMfSsMleCrbMetricDiagnose.md`
- replay 类型：小 MC / diagnostic replay / post-fix local MLE validation。
- 主要问题：在 SS 固定 in-tooth / local 条件下，检查 `SS-SF-Static`、`SS-MF-CP-K`、`SS-MF-CP-U` 的 DoA / `fdRef` / `fdRate` 是否与对应 CRB、truth/final objective probe 和 solver probe 口径一致。
- 观察范围：单星 SS，10 帧多帧 CP 模型，truth-centered in-tooth `fdRef` 范围，known-rate / unknown-rate 两个 MF 分支。
- 不覆盖范围：不覆盖 MS；不覆盖 full-flow tooth selection；不验证 subset-periodic flow；不提供论文最终 MC 曲线；不直接证明 estimator 在所有 seed / all-flow 下达到 CRB。
- truth 使用口径：truth 用于 in-tooth oracle range、objective probe、mixed-point probe、path probe 和离线评价；不进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `SS-SF-Static` | 单星单帧 static DoA-Doppler baseline。 | No | 作为 static 对照锚点。 | 用于判断 static baseline 是否健康，以及 MF 是否确实带来角度 CRB 收缩。 |
| `SS-MF-CP-K` | 单星多帧 continuous-phase known-rate local MLE。 | No | `fdRate` 固定为 known truth / known condition；估计 DoA 与 `fdRef`。 | 用于验证 known-rate MF 下 DoA 与 `fdRef` 的局部 MLE / CRB 行为。 |
| `SS-MF-CP-U` | 单星多帧 continuous-phase unknown-rate local MLE。 | No | 释放 `fdRate` nuisance。 | 用于验证 unknown-rate release 是否破坏 DoA / `fdRef`，以及低 SNR fdRate tail 是否需要 resolved 口径。 |
| `fdRef oracle half-width` | 以 truth `fdRef` 为中心的 in-tooth 搜索范围，本轮为 `0.400 tooth`。 | Oracle only | 限制 `fdRef` 搜索在同一 tooth 附近。 | 只用于局部 MLE / metric 诊断；不能外推为 full-flow tooth acquisition 已完成。 |
| `MF init mode = static-doa-truth-fd` | DoA 从 static seed 起步，频率中心使用 truth-local oracle。 | Oracle only | 降低 wrong-tooth / fd 初值污染。 | 适合诊断 DoA basin-entry 与 metric，不适合直接作为真实 runtime 初始化。 |
| `fdRefNormRmse` | `sqrt(mean((fdRefErrHz / fdRefCrbHz)^2))`。 | Evaluation only | 只改统计口径。 | 比 legacy `RMSE / median(CRB std)` 更适合判断 CRB-normalized fdRef 统计。 |
| `health` filter | 剔除 solver / frequency / boundary / no-solve 不健康样本。 | Evaluation only | 只改结果汇总。 | 用于 unknown-rate nuisance tail 的 resolved 解释，不改变 raw 结果。 |
| `trim` filter | health + CRB-normalized angle / fdRef cap。 | Evaluation only | 只改结果汇总。 | 用于检查 tail 对 RMSE 的影响；不等价于 runtime rejection。 |
| `truth/final objective probe` | 在 final、truth、truthDoA-finalFd、finalDoA-truthFd 等 mixed point 评价 objective。 | Oracle only | 只做离线 objective 对比。 | 用于判断是否仍存在系统性 DoA basin miss。 |
| `solve probe` | 重跑 baseline、truth-DoA 和 wide-static-start probe。 | Oracle / diagnostic | 只做诊断性重求解。 | 修复后 probe 可能走 estimator 默认 `doaBasinEntry*`，因此不应过度解释为纯 local truth oracle。 |
| `path probe` | 沿 final/static 到 truth DoA 方向做 objective line probe。 | Oracle only | 只观察目标函数形状。 | 用于解释是否有 truth-side objective slope；不能进入默认流程。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfSsMleCrbMetricDiagnose_20260507-174134.mat` | 2026-05-07 | `representative` | `numRepeat=5`，seed `[266,271,347,253,260]`，SNR `[-15:3:9] dB`，`fdRefOracleHalfWidth=0.400 tooth`，外层 parfor，objective/solve/path probe 全开。 | post-fix SS-MF 已恢复到约 `1.17 x CRB` 的 DoA 局部 MLE 行为；CP-K `fdRefNormRmse < 1` 作为 diagnostic caveat 记录。 | 当前代表性结果。 |
| `test/data/cache/replay/replayMfSsMleCrbMetricDiagnose_20260507-170849.mat` | 2026-05-07 | `superseded` | `numRepeat=5`，SNR `[-15,-10,-5,0,5,10] dB`，targeted basin probe。 | 证明 DoA basin-entry 修复有效，但 fdRef CRB-normalized 口径尚未补全。 | 被最新 snapshot 的更密 SNR 与 fdRef normalized metrics 覆盖。 |
| `test/data/cache/replay/replayMfSsMleCrbMetricDiagnose_20260507-123251.mat` | 2026-05-07 | `superseded` | compact-only / pre basin-entry fix。 | 低 SNR `SS-MF` angle/CRB 约 `3 x`，用于说明旧 compact local box 存在系统性 basin-entry failure。 | 被 estimator basin-entry fix 后结果覆盖。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDbList = [-15, -12, -9, -6, -3, 0, 3, 6, 9]`
- `baseSeed = 266`
- `numRepeat = 5`
- seed list：`266, 271, 347, 253, 260`
- `numFrame = 10`
- method list：`SS-SF-Static`，`SS-MF-CP-K`，`SS-MF-CP-U`
- `fdRef oracle half-width = 0.400 tooth`
- `fdRate oracle half-width = 1000 Hz/s`
- `DoA release half-width = [0.002, 0.002] deg`
- `MF init mode = static-doa-truth-fd`
- health threshold：`fdRef <= 0.250 tooth`，`fdRate <= 250 Hz/s`
- objective / solve / path probes：`1 / 1 / 1`
- checkpoint：`45/45` 完成，运行后清理 checkpoint artifacts。
- 外层 parfor：启用，`task count / outer parfor = 45 / 1`
- snapshot 保存变量：`replayData`
- 运行时间：约 `4 min 12 s`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| task 数 | 45 | `9 SNR x 5 seeds`。 |
| raw `SS-MF-CP-K` angle/CRB | `1.305 @ -15 dB`，`~1.173-1.176 @ -12..9 dB` | DoA 已从旧版本低 SNR `~3 x CRB` 回到局部 CRB-level。 |
| raw `SS-MF-CP-U` angle/CRB | `1.307 @ -15 dB`，`~1.174-1.176 @ -12..9 dB` | unknown-rate release 没有系统性破坏 DoA。 |
| health `SS-MF-CP-U` angle/CRB | `1.1338 @ -15 dB`，其余 SNR 基本与 raw 一致 | 低 SNR 剔除 fdRate nuisance tail 后，DoA 更接近 CRB-level。 |
| raw `SS-SF-Static` angle/CRB | `1.174-1.205` | static baseline 稳定贴近自身 spherical CRB。 |
| `SS-MF-CP-K` fdRefNormRmse | `0.809 @ -15 dB` 到 `0.791 @ 9 dB` | 5-seed diagnostic replay 中系统低于显示 CRB；记录为统计 / 条件化 caveat，不宣称 beat CRB。 |
| `SS-MF-CP-U` fdRefNormRmse | `0.984 @ -15 dB` 到 `0.961 @ 9 dB` | unknown-rate `fdRef` 接近 CRB；低 SNR health 后为 `1.023`。 |
| CP-U health reject | `2/5 @ -15 dB`，reject reason = `fdRate` | 低 SNR unknown-rate nuisance 有 outlier；不应解读为 DoA 主参数失败。 |
| boundary hit rate | 全部主表为 0 | 当前局部盒 / fd range 没有 boundary-hit 主导问题。 |
| CRB metric projection | `traceOverSphericalMedian = 1.1363` | spherical angle CRB 与 trace CRB 口径稳定，二者差异为常数因子。 |

### 4.3 关键对比表

#### Angle / CRB compact summary

| SNR (dB) | SS-SF-Static | SS-MF-CP-K | SS-MF-CP-U raw | SS-MF-CP-U health |
|---:|---:|---:|---:|---:|
| -15 | 1.174 | 1.305 | 1.307 | 1.134 |
| -12 | 1.184 | 1.176 | 1.176 | 1.176 |
| -9 | 1.191 | 1.175 | 1.175 | 1.175 |
| -6 | 1.196 | 1.173 | 1.174 | 1.174 |
| -3 | 1.199 | 1.174 | 1.174 | 1.174 |
| 0 | 1.201 | 1.174 | 1.174 | 1.174 |
| 3 | 1.203 | 1.174 | 1.174 | 1.174 |
| 6 | 1.204 | 1.174 | 1.174 | 1.174 |
| 9 | 1.205 | 1.173 | 1.173 | 1.173 |

#### fdRef normalized RMSE compact summary

| SNR (dB) | SS-SF-Static | SS-MF-CP-K | SS-MF-CP-U raw | 备注 |
|---:|---:|---:|---:|---|
| -15 | 1.120 | 0.809 | 0.984 | CP-U health 后为 1.023。 |
| -12 | 1.123 | 0.806 | 0.979 | CP-K 系统低于 1。 |
| -9 | 1.125 | 0.803 | 0.975 | 不是 legacy median denominator 问题。 |
| -6 | 1.126 | 0.801 | 0.972 | 同一 SNR 下 CRB std 对 seed 基本相同。 |
| -3 | 1.127 | 0.800 | 0.969 | 记录为 5-seed diagnostic caveat。 |
| 0 | 1.128 | 0.799 | 0.967 | 后续 scan 需用更多 seed 验证。 |
| 3 | 1.128 | 0.796 | 0.966 | 不能写成 CP-K 超越 CRB。 |
| 6 | 1.129 | 0.794 | 0.963 | 需要 paper-facing scan 给统计解释。 |
| 9 | 1.129 | 0.791 | 0.961 | 同上。 |

#### Tail / reject cases

| method | SNR | seed | angle err (deg) | fdRefErr/CRB | fdRate abs err (Hz/s) | reject reason | 解释 |
|---|---:|---:|---:|---:|---:|---|---|
| `SS-MF-CP-U` | -15 | 253 | 0.010661 | 0.963 | 302.90 | `fdRate` | DoA 不坏，低 SNR unknown-rate nuisance 出 tail。 |
| `SS-MF-CP-U` | -15 | 266 | 0.025073 | 0.882 | 269.53 | `fdRate` | 同上；health 后 CP-U angle/CRB 明显改善。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `SS-MF-CP-K/U` 在 `-12..9 dB` 的 angle RMSE 几乎稳定在 `1.17 x spherical CRB`，说明 basin-entry 修复后，旧的低中 SNR 系统性 local-box miss 已经基本消除。
- `SS-MF-CP-U` 与 `SS-MF-CP-K` 的 angle 结果几乎重合，且 CP-U 的 objective release compare 多数为负，说明 unknown-rate release 确实降低 objective，但没有破坏 DoA 主参数。
- `boundaryHitRate = 0`，说明当前结果不是搜索边界或 fd range 截断导致。
- `SS-SF-Static` 在所有 SNR 上维持约 `1.17-1.20 x spherical CRB`，作为 static 对照锚点是健康的。
- CRB projection table 中 `traceOverSphericalMedian = 1.1363` 稳定，说明 angle error 与 spherical CRB 投影口径已经明确；trace CRB 只是另一个固定尺度参考。

### 5.2 仍未解决或反向的现象

- `SS-MF-CP-K` 的 `fdRefNormRmse` 在 5-seed replay 中系统低于 1，范围约 `0.79-0.81`。由于本 replay 使用 oracle in-tooth range 和固定小 seed set，该现象暂记为统计 / 条件化 caveat；不能解释为真正超越 CRB。
- `SS-MF-CP-U` 在 `-15 dB` 有 `2/5` fdRate health reject，说明低 SNR unknown-rate nuisance 仍需要 resolved / health 口径。
- solve probe 中 `truth-doa-truth-fd` 仍会经过默认 `doaBasinEntry*` variant，因此它不再是纯 local truth-oracle；后续解释应优先依赖 aggregate、objective probe 与 diagnosis table。
- `-15 dB seed 347` 仍体现 residual basin/local issue：truth-DoA probe 可把 angle 从约 `0.0188 deg` 降到约 `0.00785 deg`，说明极低 SNR 个例仍可能有 basin tail。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `seed=347, SNR=-15` | residual basin / local-box tail | `SS-MF-CP-K/U` baseline angle 约 `0.0188 deg`，truth-DoA solve 可到约 `0.00785 deg`。 | 说明修复后仍有低 SNR 个例 tail，但不再是系统性问题。 |
| `seed=253, SNR=-15, CP-U` | fdRate nuisance tail | `fdRateAbsErr=302.9 Hz/s`，被 health reject，angle 仅 `0.010661 deg`。 | 表明 CP-U 低 SNR 的主要 tail 是 nuisance-rate，不是 DoA 主参数。 |
| `seed=266, SNR=-15, CP-U` | fdRate nuisance + high angle tail | `fdRateAbsErr=269.53 Hz/s`，raw angle/CRB 较高。 | health filter 后 CP-U angle/CRB 改善到 `1.1338` 的主要来源。 |
| `SNR=9, CP-K` | high-SNR diagnostic caveat | angle/CRB `1.1732`，`fdRefNormRmse=0.79068`。 | 说明 fdRef < CRB 不是低 SNR 偶然，而是当前小 seed 诊断 caveat。 |

## 6. 机制解释

### 6.1 当前解释

这轮结果说明，SS-MF 的主要 DoA 问题已经从旧版本的“static seed + compact `0.002 deg` local box 无法进入 MF basin”转变为“post-fix local MLE 在有限样本下略高于 spherical CRB”。`SS-MF-CP-K/U` 在多数 SNR 上稳定约 `1.17 x CRB`，而不是低 SNR `3 x CRB`，说明 DoA basin-entry acquisition 已经发挥作用。

CP-U 结果显示 unknown-rate release 并没有损坏 DoA：CP-U 和 CP-K 的 angle RMSE 基本一致，并且 CP-U objective 相对 CP-K 通常有小幅改善。低 SNR CP-U 的问题集中在 fdRate nuisance tail，适合用 health / resolved 口径解释，而不是继续修改 DoA estimator。

`fdRef` 方面，CP-U 接近 CRB，而 CP-K 在这个 5-seed diagnostic replay 中系统低于 CRB。由于本 replay 使用 oracle in-tooth `fdRef` 范围、truth-related initialization 条件和固定少量 seed，这个现象更适合记录为“CRB/statistical caveat”。后续需要用 scan 的更多 seed 和 paper-facing 统计口径确认，而不是在 estimator 主路径上人为调差。

### 6.2 这个结果支持什么

- 支持固定当前 SS-MF DoA basin-entry estimator 修改，不建议回退。
- 支持 `SS-MF-CP-K/U` 在 SS local / in-tooth 条件下已经进入 CRB-level 行为。
- 支持 CP-U unknown-rate release 在 DoA 主参数上基本安全，但低 SNR fdRate nuisance 需要 resolved / health 口径。
- 支持从 replay 诊断切回 scan：后续需要用 scan 做更大 seed 的 paper-facing 统计，而不是继续扩这个 replay。

### 6.3 这个结果不证明什么

- 不证明 CP-K `fdRef` 真正超越 CRB。
- 不证明 full-flow tooth acquisition / subset selection 已经通过。
- 不证明 MS-MF 会自动同样贴 CRB；MS 仍需单独 replay / scan。
- 不证明当前 5-seed RMSE 可直接用于论文主图。
- 不证明 truth-centered in-tooth oracle 条件可作为 runtime estimator 条件。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | SS DoA basin-entry 修改可阶段性固定；不建议继续修改 estimator 主核。 |
| flow 默认路径 | 不直接影响 full-flow / subset-periodic flow；本 replay 不验证 tooth selection。 |
| regression | 暂不新增 regression。`fdRef < CRB` 还不是稳定 contract；低 SNR fdRate tail 也不应固化成 pass/fail。 |
| replay 下一步 | 本 replay 可作为 representative diagnostic 固定；后续只在 estimator / CRB 口径再变时重跑。 |
| scan 下一步 | 可以回到 `scanMfMleCrbInToothConsistency`，用更多 seed 验证 SS static / SS-MF 的 angle/CRB、fdRef/CRB 和 CP-U resolved rate。 |
| 论文图 / 论文口径 | 本 replay 只支持“SS local estimator 修复有效”的机制证据；论文图应来自 scan 的 full / resolved / outlier 统计。 |
| 排障记录 | 可在主记录中补一句：SS-MF post-fix replay 已回到 CRB-level，下一步切回 scan 验证统计代表性。 |

## 8. 限制与禁止解释

- 不要把 `fdRefNormRmse < 1` 写成 estimator 超过 CRB；当前只能写成 small-repeat diagnostic caveat。
- 不要把 `health` / `trim` 结果写成 runtime filtering，二者只是离线统计口径。
- 不要把 `truth-doa-truth-fd` solve probe 当作纯 local oracle；当前 estimator 默认 basin-entry 会影响 solve variant。
- 不要用本 replay 替代 `scanMfMleCrbInToothConsistency` 的 paper-facing 大样本统计。
- 不要因为 SS 已修复就直接推断 MS 已修复；MS 有 multi-sat coherence、branch 和 conditioning 额外机制。
- 不要恢复 truth-dependent selector、truth DoA gate、truth fd adoption 或任何 oracle runtime 逻辑。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfSsMleCrbMetricDiagnose_20260507-174134.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
test/dev/replay/replayMfSsMleCrbMetricDiagnose.m
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 早期 compact-only 版本中，`SS-MF-CP-K/U` 在低中 SNR 出现约 `3 x CRB` 的 angle gap；truth/final probe 证明当时主要是 DoA basin-entry / local box 问题。
- 引入 DoA basin-entry 后，`SS-MF-CP-K/U` angle gap 明显下降，当前代表性 snapshot 稳定在约 `1.17 x CRB`。
- 本轮新增 `fdRefNormRmse / fdRefRmseOverMeanCrbStd / fdRefMseOverMeanCrb` 后，确认 CP-K `fdRef < CRB` 不是 legacy median denominator 造成；但由于样本小且条件带 oracle，仍不作为 estimator / CRB 错误结论。
