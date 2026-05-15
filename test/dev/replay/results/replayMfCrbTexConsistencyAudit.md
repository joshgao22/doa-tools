# replayMfCrbTexConsistencyAudit 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic-only / no-snapshot-required` |
| 最新代表性运行 | 2026-05-15 终端运行；该 replay 运行较快，不固定 `.mat` snapshot，后续直接重跑即可。 |
| 当前一句话结论 | 6 星场景下 CRB model dispatch、K/U EFIM loss、reported `fdRef` CRB、reference-sat invariant、summary-reuse 与 per-sat contribution 均未发现结构性错误；MS 的 `fdRef` 信息近似按 6 星叠加，DoA 信息主要由非参考星几何贡献。 |
| 决策影响 | 固定为 6 星 paper-facing CRB / TeX consistency 诊断入口；不修改 `performance/` CRB 主公式，不新增 regression，不用于 MLE 性能曲线。 |
| 下一步动作 | 后续若 `MS-MF` local MLE 仍有 own-CRB gap，优先转向 CleanBound / CleanTrim scan、bound-margin / eigen-axis error / estimator residual，而不是继续怀疑本 replay 已排除的 dispatch、reported CRB 或 summary 复用问题。 |
| 禁止误用 | 不能把本 replay 当作 estimator 已贴 CRB 的证据；不能把 `MF-Static == CP-K` 解释成 summary 复用；不能把 `CP-K/CP-U` angle CRB 接近解释成 CRB 错误。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfCrbTexConsistencyAudit.m`
- 结果文档：`test/dev/replay/results/replayMfCrbTexConsistencyAudit.md`
- replay 类型：CRB / TeX consistency 结构诊断 replay。
- 主要问题：检查 `SS/MS-SF-DoA`、`SS/MS-MF-Static`、`SS/MS-MF-CP-K`、`SS/MS-MF-CP-U` 的 CRB case dispatch、interest / nuisance 切分、time template、FIM / EFIM block、K/U Schur loss、reported `fdRef` marginal CRB、reference-sat invariant 和 per-sat contribution 是否与论文模型层级一致。
- 当前默认场景：paper-facing 6 星 multi-sat context。
- 不覆盖范围：不跑 MLE；不验证 full-flow、candidate adoption、subset tooth selection、same-tooth rescue、pair selection 或 estimator local basin。
- truth 使用口径：只用于构造 model-matched CRB context 和结果解释，不进入 runtime selector / adoption / final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `SS/MS-SF-DoA` | 单帧 DoA-only pilot-model effective CRB；不估计 `fdRef`。 | Context / evaluation | 只用于 DoA-only CRB audit。 | SF baseline 作为 DoA-only anchor；MS 的 `numSat` 应等于当前 multi-sat context 的星数。 |
| `MF-Static` | 多帧 zero-rate continuous-phase DoA-Doppler CRB。 | Context / evaluation | `fdRateMode="zero"`，interest 为 `[lat,lon,fdRef]`。 | 当前语义是 continuous static；在当前线性化点上是 CP-K 的 zero-rate / known-rate 退化。 |
| `MF-CP-K` | 多帧 affine-Doppler continuous-phase known-rate CRB。 | Context / evaluation | `fdRate` fixed，不进入 nuisance。 | known-rate benchmark；用于和 CP-U 比较 nuisance-rate loss。 |
| `MF-CP-U` | 多帧 affine-Doppler continuous-phase unknown-rate CRB。 | Context / evaluation | `fdRate` 作为 nuisance，经 Schur complement 消去。 | 用于检查 unknown `gamma_ref` 的 EFIM loss 是否真实进入。 |
| `reported fdRef CRB` | CRB helper 返回的 `fdRef` 标准差。 | Evaluation only | 只改诊断表。 | 应与 `sqrt(inv(EFIM)(3,3))` 一致；本次已验证一致。 |
| `K/U path audit` | known-rate 与 unknown-rate 的 EFIM monotonicity、loss PSD、shared block、reference invariant 检查。 | Evaluation only | 只改诊断表。 | 通过时说明 K/U Schur loss 路径自洽。 |
| `per-sat contribution` | 按卫星拆分 DoA / `fdRef` / `fdRate` 信息 proxy。 | Evaluation only | 只改诊断表。 | 用于确认 6 颗星都进入 CRB，并观察 DoA 信息由哪些卫星主导。 |
| `expected-continuous-static-degenerate` | `MF-Static` 与 `CP-K` 在 continuous zero-rate 语义下完全一致。 | No | 只改 pair audit 分类。 | 合理退化，不应列为 summary-reuse suspect。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| 无固定 snapshot | 2026-05-15 | representative terminal run | `auditMode=model-matched`，`snrDb=0 dB`，`baseSeed=253`，`numFrame=10`，6 星 `[5259 1243 348 5652 14 4437]`，参考星 `5259`。 | 6 星 CRB dispatch、K/U path、reported fdRef CRB、summary-reuse 与 per-sat contribution 均通过；MS `fdRef` CRB 约为 SS 的 `0.413x`，接近 `1/sqrt(6)`；unknown-rate 对 `fdRef` CRB 约有 `7.4%` inflation，对 angle CRB 几乎无影响。 | 取代旧 dual-sat 代表结论。 |

维护说明：本 replay 运行很快，当前不将 `.mat` snapshot 作为结果文档的必要依赖。若本地运行仍保存 snapshot，只作为临时恢复便利，不作为代表性数据归档。

## 4. 最新代表性运行

### 4.1 配置

- `auditMode = "model-matched"`
- `snrDb = 0 dB`
- `baseSeed = 253`
- `numFrame = 10`
- `tleFileName = "statlink_20260318.tle"`
- `utcRef = 2026-03-18 17:08:00 UTC`
- `usrLla = [55, 36.59, 0]`
- `selectedSatIdxGlobal = [5259 1243 348 5652 14 4437]`
- `nonRefSatIdxGlobal = [1243 348 5652 14 4437]`
- `refSatIdxGlobal = 5259`
- `refSatIdxLocal = 1`
- `numSat = 6`
- `timeOffsetSec = [-0.0053333, -0.004, -0.0026667, -0.0013333, 0, 0.0013333, 0.0026667, 0.004, 0.0053333, 0.0066667]`
- `fdRefTrueHz = -1.0435e5`
- `fdRefFitHz = -1.0435e5`
- `fdRateFitHzPerSec = -2910.7`
- method list：`SS-SF-DoA, MS-SF-DoA, SS-MF-Static, MS-MF-Static, SS-MF-CP-K, SS-MF-CP-U, MS-MF-CP-K, MS-MF-CP-U`
- snapshot：不固定；按本次终端输出记录结果。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| Dispatch pass | `8/8` | 所有 case 的 `dispatchOk=true`，没有发现走错 CRB 入口。 |
| `contextSummaryTable.numSat` | `6` | Tex audit 已切到 6 星场景，不再是 legacy dual-sat 审计。 |
| `MS-MF-CP-K` nuisance count | `66` | 6 颗星相位 nuisance + `6 × 10` frame amplitude nuisance。 |
| `MS-MF-CP-U` nuisance count | `67` | 相比 CP-K 多一个 `fdRate` nuisance。 |
| K/U path audit | `single=path-ok, multi=path-ok` | K/U monotonicity、loss PSD、shared block 和 reference invariant 均通过。 |
| `SS-MF-CP-K -> SS-MF-CP-U` fdRef CRB ratio | `0.051335 -> 0.055223 Hz`，ratio `1.0757` | unknown `gamma_ref` 对 single-sat `fdRef` CRB 有约 `7.6%` inflation。 |
| `MS-MF-CP-K -> MS-MF-CP-U` fdRef CRB ratio | `0.021223 -> 0.022791 Hz`，ratio `1.0739` | unknown `gamma_ref` 对 6 星 `fdRef` CRB 有约 `7.4%` inflation。 |
| `SS-MF-CP-K -> MS-MF-CP-K` fdRef CRB ratio | `0.051335 -> 0.021223 Hz`，ratio `0.4134` | MS `fdRef` CRB 接近 `1/sqrt(6)` 改善。 |
| `SS-MF-CP-U -> MS-MF-CP-U` fdRef CRB ratio | `0.055223 -> 0.022791 Hz`，ratio `0.4127` | unknown-rate 下同样形成近似 6 星信息增益。 |
| `angleCrbUnknownOverKnown` | `1.0` | 当前场景下 unknown-rate 对 angle CRB 几乎没有可见 inflation。 |
| Summary-reuse suspects | empty | 当前没有 reported-CRB / summary 复用嫌疑。 |
| Per-sat contribution | all `sat-ok` | 6 颗星都进入 CRB，没有 missing 或 unexpected-zero。 |

### 4.3 关键对比表

#### CRB dispatch compact table

| case | CRB function | interest | nuisance count | phase / rate | numSat | numFrame | dispatch |
|---|---|---|---:|---|---:|---:|---|
| `SS-SF-DoA` | `crbPilotSfDoaOnlyEffective` | `lat,lon` | 0 | `doa-only / none` | 1 | 1 | ok |
| `MS-SF-DoA` | `crbPilotSfDoaOnlyEffective` | `lat,lon` | 0 | `doa-only / none` | 6 | 1 | ok |
| `SS-MF-Static` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 11 | `continuous / zero` | 1 | 10 | ok |
| `MS-MF-Static` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 66 | `continuous / zero` | 6 | 10 | ok |
| `SS-MF-CP-K` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 11 | `continuous / known` | 1 | 10 | ok |
| `SS-MF-CP-U` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 12 | `continuous / unknown` | 1 | 10 | ok |
| `MS-MF-CP-K` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 66 | `continuous / known` | 6 | 10 | ok |
| `MS-MF-CP-U` | `crbPilotMfDoaDoppler` | `lat,lon,fdRef` | 67 | `continuous / unknown` | 6 | 10 | ok |

#### K/U path audit

| satMode | CRB monotonic | EFIM monotonic | loss PSD | shared block | fdRef invariant | fdRate invariant | gamma coupling finite | fdRef U/K CRB ratio | audit class |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `single` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1.0757 | `path-ok` |
| `multi` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1.0739 | `path-ok` |

#### FIM / EFIM block audit

| case | rawFimFdFd | efimFdFd | crbAngleTraceDeg | crbFdRefStdHz | EFIM cond | 备注 |
|---|---:|---:|---:|---:|---:|---|
| `SS-SF-DoA` | NaN | NaN | 0.018851 | NaN | 3.1181 | DoA-only。 |
| `MS-SF-DoA` | NaN | NaN | 0.0066702 | NaN | 1.9835 | 6 星 DoA-only CRB 明显更紧。 |
| `SS-MF-Static` | 391.03 | 379.46 | 0.0043965 | 0.051335 | 561.47 | 与 SS CP-K 一致。 |
| `MS-MF-Static` | 2346.2 | 2276.8 | 5.9066e-06 | 0.021223 | 4.1721e7 | 6 星 `fdRef` 信息接近 6 倍。 |
| `SS-MF-CP-K` | 391.03 | 379.46 | 0.0043965 | 0.051335 | 561.47 | known-rate baseline。 |
| `SS-MF-CP-U` | 391.03 | 327.92 | 0.0043965 | 0.055223 | 649.72 | unknown-rate Schur loss 生效。 |
| `MS-MF-CP-K` | 2346.2 | 2276.8 | 5.9066e-06 | 0.021223 | 4.1721e7 | 多星 known-rate baseline。 |
| `MS-MF-CP-U` | 2346.2 | 1967.5 | 5.9066e-06 | 0.022791 | 4.7958e7 | K/U loss + 6 星融合。 |

#### fdRef marginal audit

| case | efimFdFd | reported fdRef std | full-inv fdRef std | scalar EFIM fdRef std | reported rel diff | status |
|---|---:|---:|---:|---:|---:|---|
| `SS-MF-Static` | 379.46 | 0.051335 | 0.051335 | 0.051335 | 0 | ok |
| `MS-MF-Static` | 2276.8 | 0.021223 | 0.021223 | 0.020957 | `1.63e-16` | ok |
| `SS-MF-CP-K` | 379.46 | 0.051335 | 0.051335 | 0.051335 | 0 | ok |
| `SS-MF-CP-U` | 327.92 | 0.055223 | 0.055223 | 0.055223 | `1.26e-16` | ok |
| `MS-MF-CP-K` | 2276.8 | 0.021223 | 0.021223 | 0.020957 | `1.63e-16` | ok |
| `MS-MF-CP-U` | 1967.5 | 0.022791 | 0.022791 | 0.022545 | 0 | ok |

#### case-pair audit

| left | right | relation | audit class | 关键结论 |
|---|---|---|---|---|
| `SS-MF-Static` | `SS-MF-CP-K` | zero-rate-degenerate-continuous-static | `expected-continuous-static-degenerate` | 合理退化。 |
| `MS-MF-Static` | `MS-MF-CP-K` | zero-rate-degenerate-continuous-static | `expected-continuous-static-degenerate` | 多星下同样合理退化。 |
| `SS-MF-CP-K` | `SS-MF-CP-U` | known-vs-unknown-efim-loss | `observed-known-unknown-difference` | K/U loss 正常。 |
| `MS-MF-CP-K` | `MS-MF-CP-U` | known-vs-unknown-efim-loss | `observed-known-unknown-difference` | 多星下 K/U loss 正常。 |
| `SS-MF-CP-K` | `MS-MF-CP-K` | single-vs-multi-sat | `observed-ss-ms-difference` | MS `fdRef` 和 angle CRB 均显著更紧。 |
| `SS-MF-CP-U` | `MS-MF-CP-U` | single-vs-multi-sat | `observed-ss-ms-difference` | unknown-rate 下同样形成多星增益。 |
| `SS-MF-Static` | `MS-MF-Static` | single-vs-multi-sat | `observed-ss-ms-difference` | static/CP-K 退化下的多星增益一致。 |
| `SS-SF-DoA` | `MS-SF-DoA` | single-vs-multi-sat-doa-only | `observed-ss-ms-difference` | 6 星 DoA-only CRB 明显更紧。 |

#### per-sat contribution audit

| case family | 关键现象 | 解释 |
|---|---|---|
| `MS-SF-DoA` | 6 颗星均 `sat-ok`；DoA proxy fraction 约为 `0.1517, 0.2201, 0.2527, 0.0273, 0.0797, 0.2685`。 | 单帧 DoA-only 下各星均进入 FIM，`4437 / 348 / 1243` 贡献较大。 |
| `MS-MF-Static / CP-K / CP-U` | 6 颗星均 `sat-ok`；`fdRefInfoProxyFraction = 0.16667`。 | `fdRef` 公共参考星 Doppler 信息基本按 6 颗星等权叠加。 |
| `MS-MF-Static / CP-K / CP-U` | DoA proxy fraction 约为 `4.8e-06, 0.00218, 0.02214, 0.21993, 0.27192, 0.48382`。 | 多帧 DoA 信息主要由 `4437`、`14`、`5652` 提供，参考星 `5259` 不是 DoA 主贡献星。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `crbModelAuditTable` 中 8 个 case 均 `dispatchOk=true`；MS case 的 `numSat=6`，CP-U 比 CP-K 正好多一个 `fdRate` nuisance。
- `crbPathAuditAggregateTable` 中 single / multi 均为 `path-ok`，且 monotonicity、loss PSD、shared block 和 reference invariant 全部通过。
- `crbFdRefReportedRelDiff` 在 MF case 中约为 0，说明 reported `fdRef` CRB 与 `sqrt(inv(EFIM)(3,3))` 一致；当前没有 summary 复用或字段取错证据。
- `SS -> MS` 的 `fdRef` EFIM 约从 `379.46` 增至 `2276.8`，接近 6 倍；`fdRef` CRB 从 `0.051335 Hz` 降至 `0.021223 Hz`，接近 `1/sqrt(6)` 改善。
- `CP-U` 相对 `CP-K` 的 `fdRef` CRB inflation 在 single / multi 下分别为 `1.0757` 和 `1.0739`，说明 unknown-rate nuisance 信息损失稳定存在。
- `angleCrbUnknownOverKnownMedian = 1`，说明当前场景下 unknown-rate 对 angle CRB 的直接膨胀很小；这与 CleanTrim 中 K/U angle 性能接近的现象一致。
- `MF-Static` 与 `CP-K` 在 raw FIM、EFIM、full inverse 和 pair audit 上完全一致，但分类为 `expected-continuous-static-degenerate`；这固定了当前 `MF-Static` 的模型语义。
- `crbSatContributionAuditTable` 中 6 颗星均为 `sat-ok`，没有 missing 或 unexpected zero；6 星没有被静默裁掉。

### 5.2 仍未解决或反向的现象

- 当前 replay 只解释 CRB structure，不解释 MS-MF MLE 为什么在高 SNR 或 tight own-CRB 下仍可能有 `1.1x–1.4x` normalized angle gap；后者仍需 CleanBound / CleanTrim scan、bound-margin 和 eigen-axis error 诊断。
- `MS-MF` 的 angle CRB 极紧，绝对误差很小也可能对应明显 normalized gap；后续 MLE/CRB 对比不能只看 scalar RMSE/CRB。
- 如果论文需要一个与 CP-K 明显不同的 static baseline，当前 `MF-Static` 不够；需要另行定义 independent-phase static 或 `MF-IP-Static` baseline。
- 选星是否最优，本 replay 只能说明当前 L6 CRB path 自洽、per-sat contribution 健康；不能替代 `scanMfGeoSatGeometryCrbAudit` / `scanMfPairGeometryCrbAudit` 或正式 SNR scan。

### 5.3 代表性 case

| case | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `MS-MF-CP-K` | 6 星 known-rate CRB | angle trace `5.9066e-06 deg`，`fdRef` std `0.021223 Hz`。 | 证明多星 CRB 显著紧于 SS。 |
| `MS-MF-CP-U` | 6 星 unknown-rate CRB | angle trace不变，`fdRef` std 增至 `0.022791 Hz`。 | 证明 unknown-rate loss 主要落在 `fdRef`，对 angle 影响很小。 |
| `MS-MF-Static` | continuous-static degeneracy | 与 `MS-MF-CP-K` raw FIM / EFIM / full CRB 一致。 | 固定 MF-Static 当前是 zero-rate CP-K 退化。 |
| `MS-MF-CP-U` per-sat | per-sat contribution | `fdRefInfoProxyFraction=1/6`，DoA proxy 由 `4437/14/5652` 主导。 | 证明 6 颗星均进入 CRB，且 DoA 增益来自非参考星几何。 |

## 6. 机制解释

### 6.1 当前解释

当前 6 星 CRB 链路没有暴露“主公式明显错”或“summary 取值错”。`fdRef` CRB 在 SS/MS 间的变化符合 6 星信息叠加：`MS-MF-CP-K` 的 `fdRef` 标准差约为 `SS-MF-CP-K` 的 `0.413x`，接近 `1/sqrt(6)`。这与 per-sat contribution 中 `fdRefInfoProxyFraction=1/6` 一致。

`CP-K` 与 `CP-U` 的 angle CRB 接近 1 也不是异常。当前场景下 `gamma_ref` 与 angle 主参数的有效耦合较弱；消去 `gamma_ref` 后，主要损失反映在 `fdRef` 维度，约为 `7.4%` CRB inflation。这个现象已被 K/U path audit 和 fdRef marginal audit 同时确认。

`MF-Static` 与 `CP-K` 的完全一致也不是 summary 复用。当前 static 定义是 continuous zero-rate phase；当 CP-K 的 rate 是 known fixed quantity 且不进入 interest / nuisance 时，它对 `[lat,lon,fdRef]` 的导数结构可与 zero-rate continuous static 在当前线性化点完全一致。因此这类一致应写成模型退化，而不是 summary-reuse suspect。

多星 DoA 信息高度不均匀。参考星 `5259` 在多帧 DoA proxy 中几乎不是主贡献项，主要贡献来自 `4437`、`14`、`5652`。这说明当前 L6 组合的 DoA CRB 变紧确实来自多星几何，而不是参考星自身单独支撑。

### 6.2 这个结果支持什么

- 支持暂不修改 `crbPilotMfDoaDoppler` 主公式。
- 支持把 `replayMfCrbTexConsistencyAudit` 固定为 6 星 CRB / TeX structure sanity replay。
- 支持在结果和论文注释中明确：`MF-Static` 是 continuous zero-rate CP-K 退化，不是 independent-phase static baseline。
- 支持当前 6 星 L6 场景作为后续 paper-facing CleanBound / CleanTrim scan 的候选主场景。
- 支持把 unknown-rate 的主要 CRB 信息损失解释为 `fdRef` 维度损失，而不是 angle CRB 的明显膨胀。
- 支持后续 MS-MF MLE/CRB 诊断优先看 estimator/bound/eigen-axis，而不是 CRB dispatch 或 summary reuse。

### 6.3 这个结果不证明什么

- 不证明 MLE 已经达到 CRB。
- 不证明当前 6 星组合是全局最优组合。
- 不证明 MS-MF full-flow 或 local-basin 问题已解决。
- 不证明 CP/IP trade-off 已覆盖；当前 replay 没有纳入 independent-phase CRB。
- 不证明 static estimator 在 dynamic truth 下会贴 static CRB；CleanTrim 中 static `fdRef` 误差主要是 dynamic truth 下的 static model mismatch。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；本 replay 不运行 estimator。 |
| flow 默认路径 | 不改；不涉及 subset / candidate / adoption。 |
| CRB 主公式 | 不改；当前 evidence 指向 CRB path 自洽，而不是 CRB output bug。 |
| regression | 暂不新增；该 replay 仍是结构诊断，不是稳定 pass/fail 契约。未来若要固化，可只抽 case dispatch、K/U EFIM monotonicity、reference invariant、per-sat non-missing 这类小契约。 |
| replay / scan 下一步 | 进入 `scanMfMsMleCrbCleanBoundConsistency` 的 6 星 paper-facing SNR scan；必要时补 bound-margin / eigen-axis error audit。 |
| 论文图 / 论文口径 | 可作为 appendix / diagnostic 依据，帮助解释 CRB structure、continuous-static degeneracy、K/U loss 和 6 星 per-sat contribution；不直接作为主性能图。 |
| 排障记录 | 主记录可摘一句：6 星 CRB dispatch / K-U Schur loss / reported fdRef / per-sat contribution 均通过，unknown-rate loss 主要体现在 `fdRef` 而非 angle。 |

## 8. 限制与禁止解释

- 不要把 `crbSatContributionAuditTable` 的 proxy fraction 当作严格理论分解；它是诊断口径，用于发现 missing / zero / dominance，不替代完整 EFIM 分解。
- 不要把 `angleCrbUnknownOverKnown=1` 解释成 unknown-rate 没有任何影响；它说明当前场景下影响主要落在 `fdRef`，不是 angle trace。
- 不要继续把 `MF-Static == CP-K` 放入 summary-reuse suspect；当前它是 continuous-static 退化。
- 不要把 MS `fdRef` CRB 约 `1/sqrt(6)` 改善外推成 angle RMSE 一定同等改善；DoA 信息由几何主导，per-sat contribution 明显不均匀。
- 不要把本 replay 的 CRB structure 结论外推到 full-flow acquisition、estimator local-basin robustness 或 CP/IP trade-off。
- 不要把 truth 或 CRB-normalized标签迁移到 runtime selector、candidate adoption 或 final winner。

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

- 早期版本使用 legacy dual-sat context，代表结论是“MS scalar fdRef gain 被 DoA-fdRef coupling 抵消”。该结论不适用于当前 6 星 L6 场景；当前 6 星结果显示 `fdRef` CRB 接近 `1/sqrt(6)` 改善。
- 早期版本把 `MF-Static vs CP-K` 标成 `unexpected-identical`；加入 fdRef marginal / K-U path / pair audit 后，已改判为 `expected-continuous-static-degenerate`。
- 早期怀疑 `crbFdRefStdHz` 可能取错字段；新增 `fdRefStdFullInvHz` 后确认 reported value 与 full inverse 一致。
- 当前版本已不再打印 summary-reuse suspects；若未来重新出现，应优先查 reported-vs-full-inverse mismatch，而不是直接改 CRB 主公式。
