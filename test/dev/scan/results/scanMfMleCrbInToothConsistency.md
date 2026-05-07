# scanMfMleCrbInToothConsistency 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / SS paper-facing candidate / CRB-local main-figure candidate` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-222606.mat` |
| 当前一句话结论 | 500-repeat SS-only in-tooth scan 表明，`SS-MF-CP-K/U` 在 CRB-local 样本上 DoA 与 `fdRef` 均稳定接近 CRB；full / tooth-resolved 结果仍含 tooth / boundary tail，因此主图应使用 `crbLocal*` 口径并同步报告 resolved / CRB-local rate。 |
| 论文图定位 | `SS-MF`: main-figure candidate；`CP-U resolved-rate`: paper-facing auxiliary curve；full/resolved RMSE: outlier / tail diagnostic，不作为主 CRB-consistency 曲线。 |
| 决策影响 | 固定该 snapshot 作为当前 SS in-tooth MLE-vs-CRB 代表结果；不继续改 estimator，不再调 `trimNormCap=5`；后续可基于 `crbLocalAngleRmseOverCrb` / `crbLocalFdRefRmseOverCrb` 画初版论文图。 |
| 下一步动作 | 写/更新本文档与图脚本；若要扩展，下一步单独恢复 MS 或做 SS/MS 总览 scan，不在本 SS-only snapshot 上继续加诊断。 |
| 禁止误用 | 不要把 full/resolved RMSE 当作 CRB-local 主曲线；不要把 `crbLocal*` 解释成 unconditional efficiency；不要把 `CP-U` 的低 SNR RMSE 单独画出而不报告 resolved / boundary rate；不要把 in-tooth oracle range 当成 runtime tooth acquisition 已通过。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfMleCrbInToothConsistency.m`
- 结果文档：`test/dev/scan/results/scanMfMleCrbInToothConsistency.md`
- scan 类型：`paper-facing curve / local CRB consistency / resolved-regime diagnostic`
- 主要问题：在 coarse Doppler compensation 已将参考星 `fdRef` 限制到单个 Doppler tooth 内时，`SS-MF-CP-K/U` 的 in-tooth local MLE 是否能在 CRB-local 样本上与对应 CRB 对齐？
- 扫描对象：`SNR=-15:3:9 dB`，`P=10`，`SS-MF-CP-K/U`，`numRepeat=500`，`initMode=auto`，`fdRef oracle half-tooth=0.40`，`fdRate oracle half-width=1000 Hz/s`。
- 不覆盖范围：不验证 full-flow tooth acquisition；不验证 subset tooth selection；不比较 rescue / ordinary-wide / flow-like gate；不验证 MS；不形成 regression 契约。
- truth 使用口径：truth 只用于构造 oracle in-tooth 频率盒、known-rate 真值条件、CRB 计算和离线 resolved / CRB-local / top-tail 评价；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：`SS-MF` 的 `crbLocal*` 曲线可作为正文主图候选；full / resolved 曲线只用于解释 threshold / tooth / fdRate tail。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `Doppler-aided / in-tooth` | 用 truth-centered single-tooth local `fdRef` range 表示前端 coarse Doppler compensation 已消除全局 tooth ambiguity。 | Oracle / eval only | 验证 tooth-resolved local MLE 与 CRB 是否同尺度。 | 不能解释为真实 runtime tooth acquisition 已通过。 |
| `initMode=auto` | MF 方法内部自行构造 DoA 与 frequency seed，不使用 replay 中的 `static-doa-truth-fd`。 | No | 更接近论文主估计器自身初始化口径。 | 不与旧 oracle-init replay 结果直接混作同一曲线。 |
| `resolved` | solver-valid、in-tooth、无 frequency boundary hit；CP-U 还要求 `fdRate` 在健康范围内。 | Eval only | 排除明显非局部 / boundary / fdRate failure 样本。 | angle error 不参与 resolved，不能把 resolved 直接等同 CRB-local。 |
| `coreResolved` | 对 SS 等同于 resolved；对 MS 会额外要求 non-ref coherence。 | Eval only | 本轮 SS-only 中与 resolved 等价。 | 不应在 SS 结果里过度解释 coherence 字段。 |
| `trimmedCore` | 在 core 样本上用固定 CRB-normalized angle / fdRef cap 剔除极端 tail。 | Eval only | dev/statistical 旧字段。 | 不是 unconditional estimator efficiency。 |
| `crbLocal*` | 当前等价于 `trimmedCore` 的 paper-facing alias。 | Eval only | 表示 CRB-local resolved 样本，用于 MLE-vs-CRB 主曲线。 | 不是另一个新 estimator，也不是通过调参改善结果。 |
| `crbLocalRate` | `crbLocalCount / numRepeat`。 | Eval only | 表示所有 repeat 中进入 CRB-local 对比区间的比例。 | 不能忽略该比例只画 RMSE。 |
| `crbLocalKeepRate` | `crbLocalCount / coreResolvedCount`。 | Eval only | 表示 resolved/core 样本中有多少没有被 5σ CRB-normalized cap 剔除。 | 不等同于总 resolved rate。 |
| `full-sample` | 所有 finite estimator 输出。 | Eval only | 显示 threshold、tooth tail 与 boundary 风险。 | 不要求 full RMSE 必须贴 CRB。 |
| `tailSubtype` | top-tail seed 的轻量诊断标签，例如 `outside-tooth-tail`。 | Eval only | 用于解释污染来源。 | 不作为 estimator runtime selector。 |

常见口径固定如下：

- 论文主 CRB-consistency 曲线优先使用 `crbLocalAngleRmseOverCrb` 与 `crbLocalFdRefRmseOverCrb`。
- 每张 CRB-local RMSE 图必须配套 `resolvedRate` 或 `crbLocalRate`，尤其是 `CP-U`。
- `full-sample` / `resolved` 用于说明 outlier 和 tooth / boundary 风险，不作为证明 MLE 达到 CRB 的主曲线。
- `truth` 只用于 offline evaluation 与 oracle in-tooth local range，不迁移到 estimator 或 flow 的 runtime 分支。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-222606.mat` | 2026-05-07 | `representative` | `SS-MF-CP-K/U`，`initMode=auto`，`P=10`，`SNR=-15:3:9 dB`，`numRepeat=500`，`fdRef half-tooth=0.40`，`fdRate half-width=1000 Hz/s`，`trim/crbLocal cap=5 sigma` | CRB-local 样本上 DoA / `fdRef` 均稳定接近 CRB；CP-U 在低 SNR 主要损失 resolved rate，而非 CRB-local 精度。 | 当前唯一保留的 SS paper-facing 代表结果。 |
| `test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-003909.mat` | 2026-05-07 | `historical / diagnostic` | `SS/MS-MF-CP-K/U`，`P=10`，`SNR=-15:5:10 dB`，`numRepeat=200`，旧 trimmed-only 口径 | SS 可用、MS 需要 replay 诊断；用于历史说明，不作为当前 SS 主图口径。 | 被后续 SS-only `crbLocal*` 结果覆盖。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:752`
- `numRepeat = 500`
- `snrDbList = [-15, -12, -9, -6, -3, 0, 3, 6, 9]`
- `frameCountList = 10`
- active methods：`SS-MF-CP-K`, `SS-MF-CP-U`
- `initMode = auto`
- `oracleFdHalfToothFraction = 0.40`
- `oracleFdRateHalfWidthHzPerSec = 1000`
- `resolvedToothHalfWidthFraction = 0.25`
- `resolvedFdRateAbsTolHzPerSec = 250`
- `coreCoherenceFloor = 0.8`（本轮 SS-only 中不额外起作用）
- `trimNormCap = 5`
- `crbLocal* = trimmedCore*` paper-facing alias
- task count：`4500`
- checkpoint：enabled；runKey `f10_snrm15to9_seed253to752_rep500_5800bd69`；正常完成后已清理 checkpoint artifacts。
- snapshot 保存变量：`scanData`
- 运行时间：约 `1 h 7 m 0 s`

### 4.2 存档数据检查

- 顶层 snapshot 内容：`scanData`
- `scanData` 主要字段：`scanName`, `runKey`, `utcRun`, `config`, `perfTable`, `aggregateTable`, `crbLocalSummaryTable`, `failureSummaryTable`, `topTailTable`, `topTailExportTable`, `repeatOutCell`, `checkpointSummaryTable`, `checkpointCleanupReport`, `plotData`, `elapsedSec`
- 未保存大体量数据：`rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace、图片。
- warning / fail 计数：本轮命令行输出未显示导致 scan 中断的错误；checkpoint 正常清理，snapshot 正常保存。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

| case | samples / SNR | resolvedRate range | crbLocalRate range | crbLocalKeepRate range | crbLocal angle RMSE/CRB | crbLocal fdRef RMSE/CRB | 备注 |
|---|---:|---:|---:|---:|---:|---:|---|
| `SS-MF-CP-K` | 500 | `0.964–0.998` | `0.660–0.904` | `0.685–0.909` | `1.0336–1.0675` | `1.0203–1.0751` | Known-rate SS 在 CRB-local 样本上稳定贴近 CRB；低 SNR 仍有 outside-tooth tail，full/resolved 不宜作主曲线。 |
| `SS-MF-CP-U` | 500 | `0.558–0.900` | `0.554–0.894` | `0.9837–0.9970` | `1.0201–1.0553` | `0.9922–1.0443` | Unknown-rate 一旦 resolved，CRB-local 精度不差；低 SNR 主要损失来自 boundary / fdRate resolved rate。 |

### 5.2 按扫描轴汇总

#### `SS-MF-CP-K`

| SNR (dB) | resolvedRate | crbLocalRate | crbLocalKeepRate | crbLocalAngleRmseOverCrb | crbLocalFdRefRmseOverCrb | 备注 |
|---:|---:|---:|---:|---:|---:|---|
| -15 | 0.964 | 0.660 | 0.68465 | 1.0597 | 1.0315 | 低 SNR CRB-local 精度好，但总样本保留率最低。 |
| -12 | 0.970 | 0.746 | 0.76907 | 1.0408 | 1.0751 | fdRef ratio 是 CP-K 中较高点，但仍为 CRB-level。 |
| -9 | 0.976 | 0.766 | 0.78484 | 1.0440 | 1.0613 | 进入稳定区间。 |
| -6 | 0.976 | 0.714 | 0.73156 | 1.0675 | 1.0454 | angle ratio 是 CP-K 中较高点。 |
| -3 | 0.982 | 0.728 | 0.74134 | 1.0565 | 1.0310 | CRB-local 稳定。 |
| 0 | 0.986 | 0.786 | 0.79716 | 1.0414 | 1.0501 | full/resolved 仍有 tail，CRB-local 正常。 |
| 3 | 0.992 | 0.824 | 0.83065 | 1.0336 | 1.0374 | 中高 SNR 保留率提升。 |
| 6 | 0.998 | 0.856 | 0.85772 | 1.0437 | 1.0327 | resolved 几乎全通过。 |
| 9 | 0.994 | 0.904 | 0.90946 | 1.0450 | 1.0203 | 最高 SNR 下保留率最高。 |

#### `SS-MF-CP-U`

| SNR (dB) | resolvedRate | crbLocalRate | crbLocalKeepRate | crbLocalAngleRmseOverCrb | crbLocalFdRefRmseOverCrb | freqBoundaryHitRate | fdRateResolvedRate | 备注 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| -15 | 0.558 | 0.554 | 0.99283 | 1.0553 | 1.0334 | 0.308 | 0.578 | 低 SNR 主要损失是 boundary / rate 解析率，不是 CRB-local 精度。 |
| -12 | 0.630 | 0.628 | 0.99683 | 1.0201 | 0.9922 | 0.292 | 0.662 | fdRef ratio 略低于 1，但接近 CRB-level，不解释为超过 CRB。 |
| -9 | 0.676 | 0.674 | 0.99704 | 1.0469 | 1.0171 | 0.274 | 0.712 | resolved 后几乎都保留为 CRB-local。 |
| -6 | 0.652 | 0.646 | 0.99080 | 1.0458 | 1.0314 | 0.296 | 0.682 | 边界率仍高。 |
| -3 | 0.706 | 0.698 | 0.98867 | 1.0492 | 1.0092 | 0.262 | 0.728 | CRB-local 稳定。 |
| 0 | 0.736 | 0.724 | 0.98370 | 1.0475 | 1.0147 | 0.230 | 0.752 | 低中 SNR 过渡。 |
| 3 | 0.784 | 0.772 | 0.98469 | 1.0468 | 1.0443 | 0.194 | 0.798 | resolved rate 继续上升。 |
| 6 | 0.860 | 0.854 | 0.99302 | 1.0269 | 1.0432 | 0.128 | 0.862 | 高 SNR 下 CRB-local 与 resolved 基本一致。 |
| 9 | 0.900 | 0.894 | 0.99333 | 1.0376 | 1.0316 | 0.094 | 0.902 | 最高 SNR 下表现稳定。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| CRB-local angle RMSE/CRB | SNR (dB) | `crbLocalAngleRmseOverCrb` | `SS-MF-CP-K`, `SS-MF-CP-U` | Yes | 主图候选；建议加 `y=1` 参考线。 |
| CRB-local fdRef RMSE/CRB | SNR (dB) | `crbLocalFdRefRmseOverCrb` | `SS-MF-CP-K`, `SS-MF-CP-U` | Yes | `CP-U -12 dB` 略低于 1，按有限样本 / 条件样本 CRB-level 解读。 |
| Resolved / CRB-local rate | SNR (dB) | `resolvedRate`, `crbLocalRate` | `SS-MF-CP-K`, `SS-MF-CP-U` | Yes / companion | 必须配合主 RMSE 图，尤其说明 CP-U 低 SNR resolved-rate 损失。 |
| Full / resolved RMSE | SNR (dB) | full/resolved RMSE 或 RMSE/CRB | `SS-MF-CP-K`, `SS-MF-CP-U` | No / diagnostic | 用于说明 tail，不作为 estimator efficiency 主结论。 |
| Failure reason | SNR (dB) | failure reason rate | `CP-U` boundary / fdRate | Appendix / diagnostic | 解释 unknown-rate nuisance 的 resolved-rate 损失。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `SS-MF-CP-K` 在所有 SNR 的 `crbLocalAngleRmseOverCrb` 均处于 `1.0336–1.0675`，`crbLocalFdRefRmseOverCrb` 处于 `1.0203–1.0751`，说明 known-rate single-sat continuous-phase local MLE 已经达到 CRB-level。
- `SS-MF-CP-U` 在 CRB-local 样本上的 angle / fdRef 也接近 CRB：angle ratio 处于 `1.0201–1.0553`，fdRef ratio 处于 `0.9922–1.0443`；这说明 unknown-rate nuisance 一旦 resolved，并未破坏主参数局部精度。
- CP-U 的 `crbLocalKeepRate` 基本接近 1，说明 CRB-local cap 不是 CP-U 低 SNR 样本减少的主因；样本减少主要发生在 resolved 之前。
- CP-K 的 `crbLocalRate` 随 SNR 从 `0.660` 上升到 `0.904`，CP-U 从 `0.554` 上升到 `0.894`，可用于论文中同时报告精度与可解析率。

### 6.2 反向、污染或未解决现象

- full / resolved 结果仍受 tooth / boundary tail 污染；即使 CP-K 的 resolved rate 很高，仍不能把 resolved RMSE/CRB 直接当作 local CRB consistency。
- `CP-U` 在低 SNR 下 `freqBoundaryHitRate` 明显高，例如 `-15 dB` 为 `0.308`，`-12 dB` 为 `0.292`，说明 unknown-rate branch 的主要代价体现在 resolved-rate，而不是 CRB-local 样本精度。
- `CP-U -12 dB` 的 `crbLocalFdRefRmseOverCrb=0.9922` 略低于 1，不能写成超过 CRB；当前解释为 500-repeat 条件样本下的 CRB-level 波动。
- top-tail 的 compact export 显示大量 tail 被标为 `outside-tooth-tail`，说明 `tailSubtype` 已可帮助解释 full/resolved 污染来源。

### 6.3 代表性异常格点 / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `SS-MF-CP-K`, `-15 dB` | outside-tooth tail | top-tail 中 `taskSeed=293` 的 `fdRefNormErr=989.15`，failure reason 为 `outside-resolved-tooth`，tail subtype 为 `outside-tooth-tail`。 | 说明 full/top-tail 被 tooth tail 主导，不能用 full RMSE 做 CRB 主曲线。 |
| `SS-MF-CP-U`, `-15 dB` | boundary / fdRate | `freqBoundaryHitRate=0.308`，`fdRateResolvedRate=0.578`。 | 说明 unknown-rate nuisance 的低 SNR 代价主要是解析率下降。 |
| `SS-MF-CP-U`, `-12 dB` | near-CRB below-one caveat | `crbLocalFdRefRmseOverCrb=0.9922`。 | 作为 CRB-level 波动处理，不写成超过 CRB。 |
| `SS-MF-CP-K`, `9 dB` | high-SNR local consistency | `crbLocalRate=0.904`，angle/fdRef ratio 分别为 `1.0450 / 1.0203`。 | 支持高 SNR 下 local estimator 主体稳定贴 CRB。 |

## 7. 机制解释

### 7.1 当前解释

这轮 500-repeat scan 说明，前面 replay 中修复的 SS-MF DoA basin-entry 已经在更大 Monte Carlo 中表现稳定。过去低 SNR `SS-MF` 的 `angle/CRB≈3×` 问题已经不再是当前主现象；在 CRB-local 条件下，known-rate 与 unknown-rate 的 angle / `fdRef` 都稳定处于 `≈1.0–1.08× CRB` 范围内。

同时，full / tooth-resolved 分布仍然存在 tail。这说明 `tooth-resolved` 只表示没有明显跳出 coarse Doppler tooth，不等价于 CRB-local basin；真正用于 CRB consistency 的主口径应该是 `crbLocal*`，并用 `crbLocalRate` / `resolvedRate` 同步说明样本有效率。

`CP-U` 的主要代价不是 CRB-local 样本上的 DoA / `fdRef` 精度，而是 unknown-rate nuisance 使低 SNR 下更容易出现 frequency-boundary hit 或 fdRate-unresolved，从而降低 resolved / CRB-local rate。这与论文主线中 “unknown Doppler-rate 作为 nuisance 会带来信息损失和解析率损失” 的解释方向一致，但本文档不把该 scan 单独当作 known/unknown-rate EFIM 证明。

### 7.2 这个 scan 支持什么

- 支持将 `SS-MF-CP-K/U` 的 `crbLocal*` 曲线作为 paper-facing MLE-vs-CRB 主图候选。
- 支持论文图同时报告 RMSE/CRB 与 resolved / CRB-local rate，而不是只画一条误差曲线。
- 支持把 full / resolved tail 作为 threshold / tooth / boundary 风险单独解释。
- 支持当前 estimator 主路径在 SS in-tooth local regime 下已经健康，不需要继续围绕 SS 改 estimator。
- 支持后续把 MS 或 SS/MS/SF/MF 对比作为独立下一步，而不是继续修改本 SS-only scan。

### 7.3 这个 scan 不证明什么

- 不证明 full-flow tooth acquisition、subset selection 或 rescue flow 已通过。
- 不证明 full-sample / tooth-resolved RMSE 应该贴 CRB。
- 不证明 `CP-U` 在所有样本上都等价于 known-rate；它的 resolved-rate 明显低于 CP-K。
- 不证明 estimator 超过 CRB；`crbLocalFdRefRmseOverCrb < 1` 的单个点只作为 CRB-level 有限样本波动。
- 不证明 MS-MF 已接近 MS CRB；本 snapshot 只运行 SS 方法。
- 不形成 regression 契约；scan 用于 paper-facing 统计，regression 只守已稳定自动契约。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改。当前结果支持 SS estimator 主路径已经足够用于 local CRB-consistency 图。 |
| flow 默认路径 | 不改。该 scan 是 in-tooth local consistency scan，不回流到 subset / rescue flow。 |
| replay 下一步 | SS 不需要继续 replay；若恢复 MS，应使用 MS replay/scan 单独诊断，不混入本结果。 |
| regression | 不写。当前是统计结果，不是自动 pass/fail 契约。 |
| 论文图 | 可用作 SS in-tooth MLE-vs-CRB 主图候选；建议主图画 `crbLocalAngleRmseOverCrb` 与 `crbLocalFdRefRmseOverCrb`，副图/右轴画 `resolvedRate` 或 `crbLocalRate`。 |
| 排障记录 | 主记录可摘一句：`scanMfMleCrbInToothConsistency_20260507-222606` 固定 SS in-tooth CRB-local 结果，SS-MF-CP-K/U 已达到 CRB-level，后续转向 MS / CP-IP / known-unknown 信息损失。 |

### 后续仿真建议

1. 基于该 snapshot 先画 SS-only 初版图：angle RMSE/CRB、fdRef RMSE/CRB、resolved/CRB-local rate。
2. 若要进入 SS/MS 总览，另开配置恢复 MS 方法行；不要在本 SS-only 结果中补 MS 结论。
3. 若要解释 CP-U 的低 SNR rate 损失，优先结合 `failureSummaryTable` 与 known/unknown-rate CRB / EFIM scan，不在本 scan 中继续加 estimator trace。
4. 若论文图需要更平滑曲线，可后续只扩 SNR 轴或帧数轴，不必改变 `trimNormCap=5`。

## 9. 限制与禁止解释

- 不要把 offline truth 评价字段用于 runtime selector、gate、candidate adoption 或 final winner。
- 不要把 `crbLocal*` 曲线解释为所有样本的 unconditional efficiency。
- 不要只画 `CP-U` 的 CRB-local RMSE 而不报告 `resolvedRate / crbLocalRate`。
- 不要用 full / resolved tail 直接否定 estimator 主路径；tail 应单独作为 threshold / outlier 风险报告。
- 不要把 in-tooth oracle local range 直接解释为 full receiver 已完成 coarse Doppler / tooth acquisition。
- 不要把本结果迁移成 regression，除非后续形成稳定、自动、可重复的契约。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-222606.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfMleCrbInToothConsistency.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 11. 历史备注

- 早期 `20260507-003909` 同时跑 SS/MS，结论是 SS 可用、MS 需要 replay 诊断；该结果不再作为当前 SS 主图依据。
- 当前 500-repeat 结果把 SS 的 paper-facing 口径固定为 `CRB-local RMSE/CRB + resolved / CRB-local rate`；full/resolved tail 只作为 outlier 机制说明。
