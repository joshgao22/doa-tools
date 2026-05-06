# replayMfInToothFdRangeOracle 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / oracle-controlled` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat` |
| 当前一句话结论 | truth-centered half-tooth oracle 下，`MS-MF-CP-U-in-tooth` 的 tooth 与频率链已经健康，并相对单帧 static / wide baseline 有明确增益；但少数 same-tooth non-ref coherence tail 仍拉坏 RMSE。 |
| 决策影响 | 支持继续推进 same-tooth tail diagnose 与 gated basin-entry 机制；不进入 estimator / flow 默认路径，不写 regression。 |
| 下一步动作 | 用本 replay 暴露出的 tail seeds 驱动 `replayMfInToothTailCaseDiagnose`、`replayMfInToothDoaDopplerRidgeTrace` 和后续 flow-like gated rescue；paper-facing 只使用 resolved / controlled 口径，不把 full-flow wide baseline 当主图。 |
| 禁止误用 | 不能把 truth-centered half-tooth oracle 当作 runtime selector；不能用 truth-DoA oracle 进入 gate / adoption；不能写成“多星多帧已经在 full-flow 中稳定最好”。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothFdRangeOracle.m`
- 结果文档：`test/dev/replay/results/replayMfInToothFdRangeOracle.md`
- replay 类型：oracle / controlled replay。
- 主要问题：当 `fdRef` 被限制在 truth-centered 半齿内时，`MS-MF-CP-U` 是否已经具备正确 tooth 内的性能上限优势；若仍有 tail，tail 是频率链问题、wrong-tooth 问题，还是 same-tooth DoA / non-ref coherence basin 问题。
- 观察范围：`snrDb=10`、`baseSeed=253`、`numRepeat=50` 的 fixed-seed 小 MC；比较单星 / 多星、单帧 / 多帧、wide / in-tooth、known / unknown-rate 与 truth-DoA oracle。
- 不覆盖范围：不验证真实 subset-periodic flow；不验证 no-truth tooth selector；不改变 estimator 主核；不证明 full-sample MLE-vs-CRB 已完成。
- truth 使用口径：truth 用于构造 half-tooth oracle 范围、truth-DoA oracle 上限、tail label 和结果评价；这些信息不能进入 runtime selector、gate、candidate adoption 或 final winner。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `half-tooth oracle` | 将 `fdRef` 搜索范围限制在 truth-centered `0.49` tooth，即约 `±367.5 Hz`。 | Oracle only | 限制 `fdRef` tooth 范围。 | 用来剥离 wrong-tooth / comb 污染，判断正确 tooth 内 estimator 上限。 |
| `fdRate oracle box` | 将 `fdRate` 限制在 truth-centered `±1000 Hz/s`。 | Oracle only | 限制 unknown-rate 干扰参数范围。 | 用于 controlled 上限，不是 runtime 可用设置。 |
| `wide baseline` | 不使用 half-tooth oracle 的较宽 `MS-MF-CP-U` 搜索。 | No | 保留真实 wide search 风险。 | 用来暴露 wrong-tooth / comb 影响，不能作为当前 CP/IP 论文主性能口径。 |
| `in-tooth` | 使用 half-tooth oracle 后的同齿估计。 | Oracle only | 剥离 wrong-tooth，只观察 same-tooth basin。 | 若这里仍有 tail，说明问题不只是 tooth selection。 |
| `truth-DoA oracle` | 把 DoA 放在 truth 附近的离线上限候选。 | Oracle only | 只改上限评价。 | 说明多星多帧信息上限存在；不能进入 runtime flow。 |
| `single-vs-multi tail` | 同 seed 下 `MS-MF-CP-U-in-tooth` 比 `SS-MF-CP-U-in-tooth` 更差的样本。 | Evaluation only | 只改结果分类。 | 用来发现 same-tooth tail seed，后续进入 tail diagnose。 |
| `non-ref coherence floor` | 非参考星 profile coherence 的最小值。 | No for metric, truth only for tail interpretation | 只用于诊断。 | 低值说明 non-ref link 没有被同齿 DoA / local-state basin 接住。 |
| `toothIdx` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于评价和结果标注。 | 本 replay 中 in-tooth 方法 `toothHitRate=1`，不代表真实 flow 已能选齿。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat` | 2026-04-26 | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=50`，`fdOracleHalfToothFraction=0.49`，`fdRateOracleHalfWidthHzPerSec=1000`。 | in-tooth 频率链健康，MS 相对 static / wide 有明确增益，但 same-tooth coherence tail 拉坏 RMSE。 | 当前唯一代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 50`
- seed range：`253:302`
- `toothStepHz = 750`
- `fdRef` oracle 范围：truth-centered `0.49` tooth，约 `±367.5 Hz`
- `fdRate` oracle 范围：truth-centered `±1000 Hz/s`
- 关键方法：`ss-sf-static`、`ms-sf-static`、`ss/ms-mf-cp-u-wide`、`ss/ms-mf-cp-u-in-tooth`、`ms-mf-cp-k-in-tooth`、`ms-mf-cp-u-truth-doa-oracle`
- snapshot 保存变量：`replayData`
- 运行时间：repeat total 约 `6591.4 s`，dynamic methods total 约 `6316.7 s`。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `MS-MF-CP-U-in-tooth toothHitRate` | 1.00 | half-tooth oracle 成功剥离 wrong-tooth 污染。 |
| `MS-MF-CP-U-in-tooth angle median` | `0.000416 deg` | 多数 seed 上多星多帧同齿结果很好。 |
| `MS-MF-CP-U-in-tooth angle RMSE` | `0.001407 deg` | RMSE 仍被少数 tail 拉坏。 |
| `MS-MF-CP-U-in-tooth fdRef median / p95` | `0.00989 / 0.52708 Hz` | 正确 tooth 内频率链健康。 |
| `MS-MF-CP-U-in-tooth fdRate median / p95` | `4.306 / 14.295 Hz/s` | unknown-rate 估计在 in-tooth 条件下健康。 |
| `target better than SS-MF rate` | 0.76 | `MS-MF-CP-U-in-tooth` 在 38/50 seed 上角度优于 `SS-MF-CP-U-in-tooth`。 |
| `truth-DoA oracle angle RMSE` | `0.000493 deg` | 多星多帧信息上限仍明显更好。 |
| `wide MS-MF toothHitRate` | 0.46 | 非 oracle wide baseline 仍被 wrong-tooth / comb 污染。 |

### 4.3 关键对比表

#### Oracle method aggregate

| method | tooth hit | angle RMSE (deg) | angle median (deg) | angle P95 (deg) | `abs(fdRef)` median / P95 (Hz) | `abs(fdRate)` median / P95 (Hz/s) | wall time median (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 1.00 | 0.002153 | 0.001691 | 0.003610 | 130.339 / 272.776 | 3833.5 / 3833.5 | — |
| `ms-sf-static` | 0.98 | 0.001871 | 0.001466 | 0.003340 | 84.315 / 271.973 | 3833.5 / 3833.5 | — |
| `ss-mf-cp-u-wide` | 0.66 | 0.003945 | 0.000838 | 0.010985 | 122.074 / 1500.012 | 12.262 / 3833.488 | 2109 |
| `ms-mf-cp-u-wide` | 0.46 | 0.007887 | 0.005122 | 0.012580 | 747.065 / 171000.014 | 6.736 / 3333.497 | 24648 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.000976 | 0.000691 | 0.001488 | 0.00949 / 0.03389 | 7.724 / 19.338 | 1797 |
| `ms-mf-cp-u-in-tooth` | 1.00 | 0.001407 | 0.000416 | 0.003607 | 0.00989 / 0.52708 | 4.306 / 14.295 | 46476 |
| `ms-mf-cp-k-in-tooth` | 1.00 | 0.002127 | 0.001358 | 0.004084 | 0.01346 / 0.41224 | 0 / 0 | 4190 |
| `ms-mf-cp-u-truth-doa-oracle` | 1.00 | 0.000493 | 0.000343 | 0.000873 | 0.00975 / 0.03244 | 4.504 / 12.886 | 41638 |

#### Paper-claim upper-bound compare

`target = MS-MF-CP-U-in-tooth`，正的 gain 表示 target 优于 baseline。

| baseline | angle RMSE gain (deg) | angle median gain (deg) | `abs(fdRef)` median gain (Hz) | target better angle rate | 解释 |
|---|---:|---:|---:|---:|---|
| `ss-sf-static` | +0.000745 | +0.001275 | +130.329 | 0.86 | target 明显优于单星单帧 static。 |
| `ms-sf-static` | +0.000463 | +0.001050 | +84.305 | 0.88 | target 明显优于多星单帧 static。 |
| `ss-mf-cp-u-in-tooth` | -0.000431 | +0.000275 | -0.000402 | 0.76 | median / 多数 seed 更好，但 RMSE 被 tail 拉坏。 |
| `ms-mf-cp-u-wide` | +0.006479 | +0.004706 | +747.055 | 0.82 | half-tooth oracle 剥离 wide wrong-tooth 风险。 |
| `ms-mf-cp-k-in-tooth` | +0.000719 | +0.000941 | +0.003566 | 0.72 | unknown-rate 分支不是当前拖累项。 |
| `truth-doa-oracle` | -0.000915 | -0.000073 | -0.000138 | 0.32 | truth-DoA 上限仍明显更好，说明 tail 来自 DoA / local-state basin。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `MS-MF-CP-U-wide` 的 tooth hit rate 只有 `0.46`，而 `MS-MF-CP-U-in-tooth` 为 `1.00`，说明当前 full-flow / wide baseline 仍不能直接用于论文主性能图。
- `MS-MF-CP-U-in-tooth` 的 `fdRef` / `fdRate` 误差很小，说明正确 tooth 内频率链不是主瓶颈。
- `MS-MF-CP-U-in-tooth` 相对 `SS-SF-Static`、`MS-SF-Static` 与 `MS-MF-CP-U-wide` 有明显增益，支持连续相位多帧模型在 controlled / resolved 条件下有价值。
- `truth-DoA oracle` 明显优于非 oracle 方法，说明多星多帧上限存在；当前差距主要是 DoA basin / non-ref coherence 承接不足。

### 5.2 仍未解决或反向的现象

- `MS-MF-CP-U-in-tooth` 的 RMSE 高于 `SS-MF-CP-U-in-tooth`，尽管 median 和 76% seed 更好；这是 tail 主导，不应只看 median 写结论。
- 少数 tail seed 中 non-ref coherence floor 很低，说明 same-tooth 条件下仍可能出现 non-ref link collapse。
- `MS-MF-CP-U-truth-doa-oracle` 仍显著更好，说明当前 estimator / flow 还没稳定进入上限 basin。

### 5.3 代表性 seed / case

| rank | seed | single angle (deg) | multi angle (deg) | truth-DoA angle (deg) | multi coherence floor | 对结论的作用 |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 277 | 0.000556 | 0.004735 | 0.000033 | 0.0894 | 典型 same-tooth non-ref coherence collapse tail。 |
| 2 | 283 | 0.000528 | 0.003618 | 0.000519 | 0.0995 | hard tail，可由后续 wide center 诊断。 |
| 3 | 298 | 0.001311 | 0.003594 | 0.000158 | 0.1006 | hard tail，可由后续 single-MF center 诊断。 |
| 4 | 256 | 0.003420 | 0.005454 | 0.000021 | 0.00066 | 最强 collapse，non-ref coherence 几乎完全没接住。 |
| 5 | 268 | 0.001262 | 0.001960 | 0.000307 | 0.9127 | 较轻 tail / fd residual 相关负样本。 |

## 6. 机制解释

### 6.1 当前解释

这个 replay 把问题拆成了两层：第一层是 comb / wrong-tooth，第二层是 same-tooth DoA / non-ref coherence basin。half-tooth oracle 把第一层剥离后，`MS-MF-CP-U-in-tooth` 的 tooth 和频率链都已经健康，因此不能再把所有 bad seed 都归因于 `fdRef/fdRate` 或 reference-sat 语义。剩余 tail 更像是：DoA / local-state 没有进入能同时解释非参考星相干相位的好 basin。

这也解释了为什么 `MS-MF-CP-U-in-tooth` 的 median 与多数 seed 表现很好，但 RMSE 仍不好：少数 hard seed 的角度误差大约在 `0.003~0.005 deg`，足以支配 50-repeat 的 RMSE 和 P95。

### 6.2 这个结果支持什么

- 支持把 same-tooth tail 作为独立问题推进，而不是继续只修 tooth selection。
- 支持后续使用 fixed tail replay 定位 non-ref coherence collapse 与 basin-entry center。
- 支持 paper-facing 仿真采用 controlled / resolved 口径，同时报告 outlier / tail rate。

### 6.3 这个结果不证明什么

- 不证明真实 flow 已经能选对 tooth。
- 不证明 `MS-MF-CP-U` 在 full-sample / wide-search 下稳定优于单星多帧。
- 不证明可以用 truth-centered oracle 范围作为 estimator 默认搜索域。
- 不证明 truth-DoA oracle 或 tail label 可以进入 gate / adoption。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；oracle 结果只作上限诊断。 |
| flow 默认路径 | 不改；本 replay 暴露 tail，但不验证 no-truth flow。 |
| regression | 不写；tail 机制尚未形成稳定契约。 |
| replay / scan 下一步 | 固定 seeds `277, 283, 298, 256, 268` 进入 `TailCaseDiagnose` / `RidgeTrace` / gated rescue；wide wrong-tooth 另交给 subset scan / flow-like replay。 |
| 论文图 / 论文口径 | 可作为 controlled / resolved 机制证据；不作为 full-flow 主性能图。 |
| 排障记录 | 主记录和机制归并版可摘取“in-tooth 频率链健康，但 same-tooth tail 拉坏 RMSE”。 |

## 8. 限制与禁止解释

- 不要把 half-tooth oracle 当作可部署 estimator 的搜索策略。
- 不要把 `toothHitRate=1` 解读为真实 flow 已经解决 tooth selection。
- 不要只用 median 证明 MS-MF 已稳定优于 SS-MF；必须同时看 RMSE / P95 / tail seeds。
- 不要把 truth-DoA oracle 作为 runtime candidate center。
- 不要把本 replay 升级为 regression；它是机制定位和 tail seed 暴露入口。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothFdRangeOracle.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 snapshot 是后续 `TailCaseDiagnose`、`DoaDopplerRidgeTrace` 与 gated rescue replay 的 seed 来源。
- 当前结果不再支持“继续只修 comb tooth”作为唯一主线；正确 tooth 内仍有 same-tooth DoA / coherence tail。
