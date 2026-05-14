# replayMfCrbTexConsistencyAudit 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / no-snapshot-required` |
| 最新代表性 snapshot | 不固定 snapshot；该 replay 运行很快，按终端输出记录结论，后续直接重跑即可。 |
| 当前一句话结论 | CRB 主链 dispatch、reported fdRef CRB、K/U Schur loss 和 SS/MS FIM aggregation 均未发现明显错误；当前最关键机制是 MS 的 fdRef scalar information gain 会被 DoA-fdRef marginal coupling 抵消，且 MF-Static 是 CP-K 的 continuous zero-rate 退化。 |
| 决策影响 | 固定为 CRB / tex consistency 诊断入口；不修改 `performance/` CRB 主公式，不新增 regression，不用于 MLE 性能曲线。 |
| 下一步动作 | 后续若 MS-MF 仍有 gap，优先转向 pair geometry / CRB coupling scan、pure-MLE-vs-CRB 口径、generated-vs-fitted phase residual 和 eigen-axis error，而不是继续怀疑本 replay 已排除的 summary / dispatch 问题。 |
| 禁止误用 | 不能把本 replay 当作 estimator 已贴 CRB 的证据；不能把 `fdref-coupling-cancelled-gain` 解释成 CRB 输出 bug；不能把 continuous-static 与 CP-K 完全一致解释成 summary 复用。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfCrbTexConsistencyAudit.m`
- 结果文档：`test/dev/replay/results/replayMfCrbTexConsistencyAudit.md`
- replay 类型：CRB / tex consistency 结构诊断 replay。
- 主要问题：检查 `SS/MS-SF-DoA`、`SS/MS-MF-Static`、`SS/MS-MF-CP-K`、`SS/MS-MF-CP-U` 的 CRB case dispatch、interest / nuisance 切分、time template、FIM / EFIM block、fdRef marginal CRB 和 DoA-fdRef coupling 是否与论文模型层级一致。
- 观察范围：`auditMode="model-matched"`，`numFrame=10`，`snrDb=0 dB`，`baseSeed=253`，卫星对为 `starlink_pair_4154_1165_20260318_170800.tle` 中的 `[1 2]`，参考星 local index 为 `1`。
- 不覆盖范围：不跑 MLE；不验证 full-flow、candidate adoption、subset tooth selection、same-tooth rescue、pair selection 或 estimator local basin。
- truth 使用口径：只用于构造 model-matched CRB context 和结果解释，不进入 runtime selector / adoption / final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `SS/MS-SF-DoA` | 单帧 DoA-only pilot-model effective CRB；不估计 `fdRef`。 | Context / evaluation | 只用于 DoA-only CRB audit。 | SF baseline 现在只作为 DoA-only anchor，不再保留 SF-static Doppler 估计。 |
| `MF-Static` | 多帧 zero-rate continuous-phase DoA-Doppler CRB。 | Context / evaluation | `fdRateMode="zero"`，interest 为 `[lat,lon,fdRef]`。 | 当前语义是 continuous static；它是 CP-K 的 zero-rate / known-rate 退化。 |
| `MF-CP-K` | 多帧 affine-Doppler continuous-phase known-rate CRB。 | Context / evaluation | `fdRate` fixed，不进入 nuisance。 | 与 `MF-Static` 在主参数导数上可能完全一致。 |
| `MF-CP-U` | 多帧 affine-Doppler continuous-phase unknown-rate CRB。 | Context / evaluation | `fdRate` 作为 nuisance，经 Schur complement 消去。 | 用于检查 K/U EFIM loss 是否真实进入。 |
| `reported fdRef CRB` | CRB helper 返回的 `fdRef` 标准差。 | Evaluation only | 只改诊断表。 | 应与 `sqrt(inv(EFIM)(3,3))` 一致；本次已验证一致。 |
| `scalar EFIM fdRef std` | `sqrt(1 / EFIM(3,3))`，假设 DoA 已知或忽略 DoA-fdRef coupling 的标量口径。 | Evaluation only | 只改诊断表。 | 可显示 MS 的 raw fdRef information gain，但不是 full marginal CRB。 |
| `fdRef coupling inflation` | `fdRefStdFullInvHz / fdRefStdScalarEfimHz`。 | Evaluation only | 只改诊断表。 | 大于 1 表示 DoA-fdRef marginal coupling 抬高 fdRef CRB。 |
| `fdRefDoaCanonicalCorr` | fdRef 与 DoA block 的 canonical correlation。 | Evaluation only | 只改诊断表。 | 解释 scalar fdRef information gain 是否被 DoA marginalization 抵消。 |
| `expected-continuous-static-degenerate` | `MF-Static` 与 `CP-K` 在 continuous zero-rate 语义下完全一致。 | No | 只改 pair audit 分类。 | 合理退化，不再列为 summary-reuse suspect。 |
| `fdref-coupling-cancelled-gain` | MS 的 scalar fdRef 信息有增益，但 full marginal fdRef CRB 无改善。 | Evaluation only | 只改诊断标签。 | 不是 reported CRB bug，而是 DoA-fdRef coupling 机制。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| 无固定 snapshot | 2026-05-14 | representative terminal run | `auditMode=model-matched`，`snrDb=0 dB`，`baseSeed=253`，`numFrame=10`，8 个主链 case。 | CRB dispatch / reported fdRef / K-U Schur loss 正常；MS fdRef scalar gain 被 DoA-fdRef coupling 抵消；MF-Static 是 CP-K continuous zero-rate 退化。 | 后续直接重跑脚本即可复现。 |

维护说明：本 replay 运行很快，当前不将 `.mat` snapshot 作为结果文档的必要依赖。若本地运行仍保存 snapshot，只作为临时恢复便利，不作为代表性数据归档。

## 4. 最新代表性运行

### 4.1 配置

- `auditMode = "model-matched"`
- `snrDb = 0 dB`
- `baseSeed = 253`
- `numFrame = 10`
- `timeOffsetSec = [-0.0053333, -0.004, -0.0026667, -0.0013333, 0, 0.0013333, 0.0026667, 0.004, 0.0053333, 0.0066667]`
- selected local sats：`[1 2]`
- reference local sat：`1`
- `fdRefTrueHz = -28940`
- `fdRefFitHz = -28940`
- `fdRateFitHzPerSec = -3833.5`
- method list：`SS-SF-DoA, MS-SF-DoA, SS-MF-Static, MS-MF-Static, SS-MF-CP-K, SS-MF-CP-U, MS-MF-CP-K, MS-MF-CP-U`
- snapshot：不固定；按本次终端输出记录结果。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| Dispatch pass | `8/8` | 所有 case 的 `dispatchOk=true`，没有发现走错 CRB 入口。 |
| `SS-MF-Static` vs `SS-MF-CP-K` | identical | continuous zero-rate static 与 CP-K 在主参数 EFIM 上退化一致。 |
| `MS-MF-Static` vs `MS-MF-CP-K` | identical | 多星下同样为合理退化，而不是 summary-reuse。 |
| `SS-MF-CP-K -> SS-MF-CP-U` fdRef CRB ratio | `0.051335 -> 0.055223 Hz` | unknown rate nuisance 带来 fdRef CRB inflation。 |
| `MS-MF-CP-K -> MS-MF-CP-U` fdRef CRB ratio | `0.051335 -> 0.053315 Hz` | 多星下 K/U 差异仍存在，但受 coupling 改写。 |
| `SS-MF-CP-K -> MS-MF-CP-K` scalar fdRef std | `0.051335 -> 0.036299 Hz` | 如果只看 scalar EFIM，MS 给出约 `1/sqrt(2)` 改善。 |
| `SS-MF-CP-K -> MS-MF-CP-K` full marginal fdRef std | `0.051335 -> 0.051335 Hz` | full inverse 后 fdRef 增益被 DoA-fdRef coupling 抵消。 |
| `MS-MF-CP-K` EFIM condition | `3.7892e9` | MS EFIM 极度病态，后续 MLE/CRB 应看 whitened / eigen-axis error。 |
| Summary-reuse suspects | empty | 当前没有真正的 reported-CRB / summary 复用嫌疑。 |

### 4.3 关键对比表

#### CRB dispatch compact table

| case | CRB function | interest | nuisance count | phase / rate | numSat | numFrame | dispatch |
|---|---|---|---:|---|---:|---:|---|
| `SS-SF-DoA` | `crbPilotSfDoaOnlyEffective` | `lat,lon` | 0 | `doa-only / none` | 1 | 1 | ok |
| `MS-SF-DoA` | `crbPilotSfDoaOnlyEffective` | `lat,lon` | 0 | `doa-only / none` | 2 | 1 | ok |
| `SS-MF-Static` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 11 | `continuous / zero` | 1 | 10 | ok |
| `MS-MF-Static` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 22 | `continuous / zero` | 2 | 10 | ok |
| `SS-MF-CP-K` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 11 | `continuous / known` | 1 | 10 | ok |
| `SS-MF-CP-U` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 12 | `continuous / unknown` | 1 | 10 | ok |
| `MS-MF-CP-K` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 22 | `continuous / known` | 2 | 10 | ok |
| `MS-MF-CP-U` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 23 | `continuous / unknown` | 2 | 10 | ok |

#### FIM / EFIM block audit

| case | rawFimFdFd | efimFdFd | crbAngleTraceDeg | crbFdRefStdHz | EFIM cond | 备注 |
|---|---:|---:|---:|---:|---:|---|
| `SS-SF-DoA` | NaN | NaN | 0.0082264 | NaN | 1.5073 | DoA-only。 |
| `MS-SF-DoA` | NaN | NaN | 0.0081506 | NaN | 1.5542 | second-sat DoA-only gain 很小。 |
| `SS-MF-Static` | 391.03 | 379.46 | 0.0025446 | 0.051335 | 1020.5 | 与 SS CP-K 一致。 |
| `MS-MF-Static` | 782.06 | 758.93 | 0.0016514 | 0.051335 | 3.7892e9 | scalar fdRef 信息翻倍，marginal CRB 不变。 |
| `SS-MF-CP-K` | 391.03 | 379.46 | 0.0025446 | 0.051335 | 1020.5 | known-rate baseline。 |
| `SS-MF-CP-U` | 391.03 | 327.92 | 0.0025446 | 0.055223 | 1180.9 | unknown-rate Schur loss 生效。 |
| `MS-MF-CP-K` | 782.06 | 758.93 | 0.0016514 | 0.051335 | 3.7892e9 | fdRef gain 被 coupling 抵消。 |
| `MS-MF-CP-U` | 782.06 | 655.84 | 0.0016515 | 0.053315 | 3.8094e9 | K/U loss + coupling。 |

#### fdRef marginal audit

| case | efimFdFd | reported fdRef std | full-inv fdRef std | scalar EFIM fdRef std | coupling inflation | coupling loss fraction |
|---|---:|---:|---:|---:|---:|---:|
| `SS-MF-Static` | 379.46 | 0.051335 | 0.051335 | 0.051335 | 1.0000 | `~0` |
| `MS-MF-Static` | 758.93 | 0.051335 | 0.051335 | 0.036299 | 1.4142 | 0.5000 |
| `SS-MF-CP-K` | 379.46 | 0.051335 | 0.051335 | 0.051335 | 1.0000 | `~0` |
| `SS-MF-CP-U` | 327.92 | 0.055223 | 0.055223 | 0.055223 | 1.0000 | `~0` |
| `MS-MF-CP-K` | 758.93 | 0.051335 | 0.051335 | 0.036299 | 1.4142 | 0.5000 |
| `MS-MF-CP-U` | 655.84 | 0.053315 | 0.053315 | 0.039048 | 1.3653 | 0.46357 |

#### eigen / coupling audit

| case | eigInfoMin | eigInfoMid | eigInfoMax | cond | fdRefDoaCanonicalCorr | fdRefDoaCouplingFraction | fdRef std if DoA known | fdRef full marginal std | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `SS-MF-CP-K` | 379.46 | 2.5691e5 | 3.8724e5 | 1020.5 | 1.3153e-05 | `~0` | 0.051335 | 0.051335 | SS coupling 很弱。 |
| `MS-MF-CP-K` | 379.46 | 3.6668e5 | 1.4379e12 | 3.7892e9 | 0.70711 | 0.5000 | 0.036299 | 0.051335 | MS fdRef scalar gain 被 coupling 抵消。 |
| `SS-MF-CP-U` | 327.92 | 2.5691e5 | 3.8724e5 | 1180.9 | 1.2226e-05 | `~0` | 0.055223 | 0.055223 | unknown-rate loss 主要投到 fdRef。 |
| `MS-MF-CP-U` | 351.81 | 3.6665e5 | 1.3402e12 | 3.8094e9 | 0.68086 | 0.46357 | 0.039048 | 0.053315 | K/U loss 与 MS coupling 共同作用。 |

#### case-pair audit

| left | right | relation | audit class | 关键结论 |
|---|---|---|---|---|
| `SS-MF-Static` | `SS-MF-CP-K` | zero-rate-degenerate-continuous-static | `expected-continuous-static-degenerate` | 合理退化。 |
| `MS-MF-Static` | `MS-MF-CP-K` | zero-rate-degenerate-continuous-static | `expected-continuous-static-degenerate` | 合理退化。 |
| `SS-MF-CP-K` | `SS-MF-CP-U` | known-vs-unknown-efim-loss | `observed-known-unknown-difference` | K/U loss 正常。 |
| `MS-MF-CP-K` | `MS-MF-CP-U` | known-vs-unknown-efim-loss | `observed-known-unknown-difference` | 多星下 K/U loss 正常。 |
| `SS-MF-CP-K` | `MS-MF-CP-K` | single-vs-multi-sat | `fdref-coupling-cancelled-gain` | scalar fdRef gain 被 full marginal coupling 抵消。 |
| `SS-MF-CP-U` | `MS-MF-CP-U` | single-vs-multi-sat | `observed-ss-ms-difference` | full marginal fdRef 仍有小差异，angle CRB 明显更紧。 |
| `SS-MF-Static` | `MS-MF-Static` | single-vs-multi-sat | `fdref-coupling-cancelled-gain` | static 下同样存在 coupling cancellation。 |
| `SS-SF-DoA` | `MS-SF-DoA` | single-vs-multi-sat-doa-only | `observed-ss-ms-difference` | second-sat DoA-only contribution 很弱。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `crbModelAuditTable` 中 8 个 case 均 `dispatchOk=true`，且 CP-U 比 CP-K 正好多出 `fdRate` nuisance；因此当前没有 evidence 支持 case dispatch 错误。
- `crbFdRefReportedRelDiff` 在 MF case 中约为 `1e-16`，说明 reported `fdRef` CRB 与 `sqrt(inv(EFIM)(3,3))` 一致；因此 reported CRB 不是 summary 复用或字段取错。
- `MS-MF-CP-K` 的 `efimFdFd` 从 `379.46` 增至 `758.93`，但 full marginal `fdRef` std 不变；`fdRefDoaCanonicalCorr=0.70711`、`fdRefDoaCouplingFraction=0.5` 直接解释了这个抵消机制。
- `MS-MF-CP-K` 的 `crbAngleTraceDeg` 从 SS 的 `0.0025446` 下降到 `0.0016514`，说明 MS 确实提供了更紧 DoA CRB；但 EFIM 条件数增至 `3.7892e9`，提示后续 MLE/CRB 对比必须看 whitened / eigen-axis 指标。
- `MF-Static` 与 `CP-K` 在 raw FIM、EFIM、full inverse 和 pair audit 上完全一致，但分类为 `expected-continuous-static-degenerate`；这固定了当前 `MF-Static` 的模型语义。
- `SS/MS-SF-DoA` 的 MS gain 很小，且 derivative preview 中 sat2 的 DoA FIM block 远小于 sat1，继续支持 DoA-only effective CRB 的弱 second-sat contribution 解释。

### 5.2 仍未解决或反向的现象

- 当前 replay 只解释 CRB structure，不解释 MS-MF MLE 为什么不能稳定贴 own-CRB；后者仍需 estimator / pair / phase-mismatch replay 或 scan。
- `MS-MF-CP-K/U` 的 EFIM 极度病态，意味着 angle scalar MSE/CRB 容易误导；需要在 MLE 结果中补 whitened second moment 与 eigen-axis projected error。
- 如果论文需要一个与 CP-K 明显不同的 static baseline，当前 `MF-Static` 不够；需要另行定义 independent-phase static 或 `MF-IP-Static` baseline。
- 选星是否导致 MS performance 不好，本 replay 只能提示当前 pair 的 MS EFIM 病态和 coupling 强；不能替代 pair geometry / CRB scan。

### 5.3 代表性 case

| case | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `MS-MF-CP-K` | fdRef coupling cancellation | scalar fdRef std `0.036299 Hz`，full marginal `0.051335 Hz`。 | 证明 MS fdRef scalar gain 被 DoA-fdRef coupling 抵消。 |
| `MS-MF-CP-U` | K/U loss with coupling | full marginal fdRef `0.053315 Hz`，canonical corr `0.68086`。 | 证明 unknown-rate loss 存在，但被 MS coupling 改写。 |
| `MS-MF-Static` | continuous-static degeneracy | 与 `MS-MF-CP-K` raw FIM / EFIM / full CRB 完全一致。 | 固定 MF-Static 当前是 zero-rate CP-K 退化。 |
| `MS-SF-DoA` | weak second-sat DoA-only gain | angle trace `0.0082264 -> 0.0081506 deg`。 | 支持 SF-DoA 只作弱多星 anchor，不作为 dynamic failure 证据。 |

## 6. 机制解释

### 6.1 当前解释

当前 CRB 链路没有暴露出“主公式明显错”或“summary 取值错”。`fdRef` CRB 看起来在 SS/MS 间不变，是因为 reported value 采用 full marginal inverse，而不是 scalar `1/EFIM_fd,fd`。MS 让 `EFIM_fd,fd` 翻倍，但也引入强 DoA-fdRef coupling；边际化 DoA 后，fdRef 的 scalar information gain 被 coupling inflation 正好抵消。

`MF-Static` 与 `CP-K` 的完全一致也不是异常。当前 static 定义是 continuous zero-rate phase：`psi = 2*pi*fd*t_abs`。当 CP-K 的 rate 是 known fixed quantity 且不进入 interest / nuisance 时，它对 `[lat,lon,fdRef]` 的导数结构可与 zero-rate continuous static 完全一致。因此这类一致应写成模型退化，而不是 summary-reuse suspect。

MS 的 DoA CRB 明显更紧，但 EFIM 极度病态。这个结果解释了为什么后续 MS-MF 的 normalized angle gap 不能只用 scalar angle RMSE/CRB 判断：MS CRB 的强主轴 / 弱主轴结构很尖，估计误差若不沿 CRB 主轴匹配，就会出现很大的 normalized gap。

### 6.2 这个结果支持什么

- 支持暂不修改 `crbPilotMfDoaDoppler` 主公式。
- 支持把 `replayMfCrbTexConsistencyAudit` 固定为轻量 CRB structure sanity replay。
- 支持在结果和论文注释中明确：`MF-Static` 是 continuous zero-rate CP-K 退化，不是 independent-phase static baseline。
- 支持后续 MS-MF MLE/CRB 诊断必须加入 whitened / eigen-axis error，而不是只看 spherical angle scalar。
- 支持把 pair selection / geometry 病态作为下一阶段检查项：当前 pair 的 MS EFIM condition 和 DoA-fdRef coupling 已非常强。

### 6.3 这个结果不证明什么

- 不证明 MLE 已经达到 CRB。
- 不证明当前 satellite pair 是 paper-main 最优 pair。
- 不证明 MS-MF full-flow 或 local-basin 问题已解决。
- 不证明 CP-K / CP-U 的 generation model 与 estimator / CRB 在相位层完全 matched；这需要 generated-vs-fitted phase residual 另查。
- 不证明 independent-phase baseline 或 CP/IP trade-off 已覆盖；当前 replay 第一版没有纳入 IP。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；本 replay 不运行 estimator。 |
| flow 默认路径 | 不改；不涉及 subset / candidate / adoption。 |
| CRB 主公式 | 不改；当前 evidence 指向 coupling mechanism，而不是 CRB output bug。 |
| regression | 暂不新增；该 replay 仍是结构诊断，不是稳定 pass/fail 契约。未来若要固化，可只抽 case dispatch、K/U EFIM monotonicity、duplicate-sat FIM aggregation 这类小契约。 |
| replay / scan 下一步 | 若要解释 MS-MF gap，优先做 `scanMfPairGeometryCrbAudit`、pure-MLE-vs-CRB replay、generated-vs-fitted phase residual replay、MLE eigen-axis error audit。 |
| 论文图 / 论文口径 | 可作为 appendix / diagnostic 依据，帮助解释 CRB structure、continuous-static degeneracy 和 MS coupling；不直接作为主性能图。 |
| 排障记录 | 主记录可摘一句：CRB dispatch / fdRef reported / K-U Schur loss 正常，MS fdRef scalar gain 被 DoA-fdRef coupling 抵消。 |

## 8. 限制与禁止解释

- 不要把 `fdRefStdScalarEfimHz` 当成正式 `fdRef` CRB；正式 comparison 应使用 full marginal `sqrt(inv(EFIM)(3,3))`。
- 不要把 MS `fdRef` full marginal CRB 不变解释为 summary bug；本 replay 已验证 reported value 与 full inverse 一致。
- 不要继续把 `MF-Static == CP-K` 放入 summary-reuse suspect；当前它是 continuous-static 退化。
- 不要把 `MS-MF` angle CRB 变紧直接等价为 MLE 一定更容易贴 CRB；MS EFIM 同时极度病态。
- 不要把本 replay 的 CRB structure 结论外推到 pair selection、full-flow acquisition 或 estimator local-basin robustness。
- 不要把 truth 或 CRB-normalized posterior 标签迁移到 runtime selector、candidate adoption 或 final winner。

## 9. 恢复与复现

该 replay 运行很快，当前不固定 snapshot。直接运行：

```matlab
run('test/dev/replay/replayMfCrbTexConsistencyAudit.m')
```

若本地临时保存了 snapshot，可按常规方式恢复：

```matlab
snapshotFile = 'test/data/cache/replay/replayMfCrbTexConsistencyAudit_<runKey>.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfCrbTexConsistencyAudit.m`
```

只运行 summary output 小节即可重出 compact tables。

## 10. 历史备注

- 早期版本把 `MF-Static vs CP-K` 标成 `unexpected-identical`；加入 fdRef marginal / eigen coupling audit 后，已改判为 `expected-continuous-static-degenerate`。
- 早期怀疑 `crbFdRefStdHz` 可能取错字段；新增 `fdRefStdFullInvHz` 后确认 reported value 与 full inverse 一致。
- `MS fdRef CRB 不变` 的解释从“可能 summary 复用”更新为“scalar information gain 被 DoA-fdRef coupling cancellation 抵消”。
- 当前版本已不再打印 summary-reuse suspects；若未来重新出现，应优先查 reported-vs-full-inverse mismatch，而不是直接改 CRB 主公式。
