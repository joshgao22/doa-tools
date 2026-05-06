# replayMfInToothOrdinaryAngleMissDiagnose 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / diagnostic` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat` |
| 当前一句话结论 | 100-repeat 证明 ordinary non-collapse miss 确实存在，final-centered small polish 完全无效；ordinary miss 主要由 `wide-basin-entry` 救回，并接近 truth-DoA oracle。 |
| 决策影响 | 推动新增 / 扩展 `replayMfInToothOrdinaryWideGateDiagnose`；不改变 estimator / flow 默认路径。 |
| 下一步动作 | 该 snapshot 已完成第一版诊断使命；后续使用 500-repeat ordinary-wide gate 结果作为更强代表。 |
| 禁止误用 | 不能把 4/4 ordinary miss 被 wide 救回解读成 blanket wide 可用；不能用 single-MF 或 final-small-polish 修 ordinary；不能把 truth label 进入 runtime gate。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothOrdinaryAngleMissDiagnose.m`
- 结果文档：`test/dev/replay/results/replayMfInToothOrdinaryAngleMissDiagnose.md`
- replay 类型：oracle / controlled mechanism triage。
- 主要问题：controlled in-tooth 条件下，raw `MS-MF-CP-U-in-tooth` 剩余 angle miss 中的 ordinary non-collapse 子类是什么、能被哪些候选救回、不能被哪些候选救回。
- 观察范围：truth-centered half-tooth `fdRef` 范围；比较 default、final-small-polish、wide-basin-entry、single-MF basin-entry、family-safe hard-collapse rescue 和 truth-DoA oracle。
- 不覆盖范围：不做 subset tooth selection；不验证 ordinary-wide truth-free gate 稳定性；不改变 estimator / flow / objective；不进入 regression。
- truth 使用口径：truth 只用于 offline miss label、caseRole、truth-DoA oracle 上限和 rescue / damage 评价；`default / wide / single-MF / final-small-polish` 本身仍是可实现 probe。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `ordinary-angle-miss` | same tooth、fd healthy、non-ref coherence 不 collapse，但 angle 超过 `0.002 deg` 的 non-collapse miss。 | Evaluation only | 只改离线分类。 | 说明剩余 tail 不全是 hard-collapse；不能直接用 coherence gate 捕获。 |
| `final-small-polish` | 围绕 raw final 点做很小 DoA 网格扰动。 | No | 只改 final 附近局部 DoA。 | 本轮完全无效，说明 ordinary miss 不是 final 点周围小扰动问题。 |
| `wide-basin-entry` | 从 wide DoA center 重新进入 basin 的候选。 | No | 改 DoA basin center。 | 本轮 4/4 ordinary miss 被救回，但还缺 truth-free adoption gate。 |
| `single-mf-basin-entry` | 用单星多帧中心作为 basin-entry candidate。 | No | 改 DoA basin center。 | 对 ordinary miss 不稳，不能和 hard-collapse 的 `wide + single-MF` 结论混用。 |
| `family-safe-adopt` | hard-collapse 线的 gated `wide + single-MF` family-aware adoption。 | No | 改 hard-collapse candidate adoption。 | 本 replay 同时复现其有效性，但它不是 ordinary miss 的 gate 证据。 |
| `fd-not-healthy-negative` | fdRef / fdRate 不健康但 angle 也可能坏的负样本。 | Evaluation only | 只改离线分类。 | 不应混入 ordinary DoA rescue 评价。 |
| `truth-DoA oracle` | 离线 truth-centered DoA 上限。 | Oracle only | 只改评价上限。 | 说明 ordinary miss 可达，不等于 runtime candidate。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat` | 2026-05-05 | representative | `snrDb=10`，`baseSeed=253`，`numRepeat=100`，`ordinaryCoh=0.950`，`rescueCoh=0.200`，`safeStep=0.004 deg`，`familyStep=0.006 deg`，`joint=0`，`family=1`。 | ordinary miss 存在但占比不高；small polish 无效；wide center 能救回全部 ordinary miss；single-MF 不稳。 | 后续 ordinary-wide 500-repeat 对 gate 稳定性给出更强结论。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `baseSeed = 253`
- `numRepeat = 100`
- seed range：`253:352`
- `fd oracle half-tooth fraction = 0.49`
- `fdRate oracle half-width = 1000 Hz/s`
- `ordinaryCoherenceThreshold = 0.950`
- `gatedRescueCoherenceThreshold = 0.200`
- `gated rescue phase gate = inactive`
- `safe adopt max DoA step = 0.004 deg`
- `family safe max DoA step = 0.006 deg`
- `safe adopt coherence threshold = 0.950`
- `safe adopt max fdRef step = 300 Hz`
- `angleHitThresholdDeg = 0.002`
- checkpoint：`100/100` 完成，成功保存轻量 `replayData` 后清理 checkpoint run 目录。
- 运行时间：parfor batch wall time 约 `24 min 38 s`。

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `numSeed` | 100 | 当前小 MC 总 seed 数。 |
| `numAngleMiss` | 10 | raw `MS-MF-CP-U-in-tooth` 超过 hit threshold 的样本。 |
| `numOrdinaryAngleMiss` | 4 | same-tooth、fd healthy、non-collapse ordinary miss。 |
| `collapseHardCount` | 4 | hard-collapse 由 family-safe 线处理。 |
| `wrongToothCount` | 0 | controlled in-tooth 条件剥离 wrong-tooth。 |
| `fdNotHealthyCount` | 5 | fd-negative 不应混入 ordinary DoA gate。 |
| `easyHitCount` | 87 | raw 已 hit 样本，用于检查 blanket wide damage。 |
| default ordinary median angle | 0.0030293 deg | ordinary miss raw median。 |
| default ordinary median coherence | 1 | ordinary miss 不表现为 coherence collapse。 |
| best implementable rescue rate | 1 | 本轮 4 个 ordinary miss 均可由 wide candidate 救回。 |
| best implementable damage rate | 0 | ordinary miss 上没有 damage。 |

### 4.3 关键对比表

#### Ordinary candidate family aggregate

| candidate group | ordinary miss | hit rate | rescue rate | damage rate | median angle (deg) | P95 angle (deg) | median gain (deg) | median coherence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `default` | 4 | 0 | 0 | 0 | 0.0030293 | 0.0042986 | 0 | 1 |
| `final-small-polish` | 4 | 0 | 0 | 0 | 0.0030293 | 0.0042986 | 0 | 1 |
| `wide-basin-entry` | 4 | 1 | 1 | 0 | 0.00029545 | 0.00049975 | 0.0027666 | 1 |
| `single-mf-basin-entry` | 4 | 0.25 | 0.25 | 0.75 | 0.0075701 | 0.00847 | -0.0049796 | 0.99304 |
| `best-implementable` | 4 | 1 | 1 | 0 | 0.00029545 | 0.00049975 | 0.0027666 | 1 |
| `truth-doa-oracle` | 4 | 1 | 1 | 0 | 0.00028681 | 0.00047371 | 0.0027814 | 1 |

#### Baseline compare：SS-SF 与 SS-MF

| baseline | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | fdRef median abs (Hz) | fdRate median abs (Hz/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 | 138.4 | 3833.5 |
| `ss-mf-cp-u-in-tooth` | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 | 0.0095962 | 6.8894 |

#### In-tooth method aggregate

| method | family | tooth hit | RMSE (deg) | median (deg) | P95 (deg) | P99 (deg) | hit rate @0.002 deg | max (deg) | coherence median |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ss-sf-static` | static-baseline | 0.92 | 0.0022363 | 0.0018886 | 0.0038394 | 0.00418 | 0.54 | 0.0052371 | — |
| `ms-sf-static` | static-baseline | 0.99 | 0.0020359 | 0.0015905 | 0.0034644 | 0.0045189 | 0.66 | 0.004673 | — |
| `ss-mf-cp-u-in-tooth` | single-center | 1.00 | 0.0009994 | 0.00072022 | 0.0017592 | 0.0024848 | 0.96 | 0.0034199 | — |
| `ms-mf-cp-u-in-tooth` | default-target | 1.00 | 0.0013708 | 0.00049554 | 0.0036174 | 0.0053885 | 0.90 | 0.0054541 | 1 |
| `ms-mf-cp-u-wide-doa-in-tooth` | wide-center | 1.00 | 0.0021782 | 0.00044087 | 0.0037959 | 0.011068 | 0.93 | 0.014519 | 1 |
| `ms-mf-cp-u-truth-doa-oracle` | truth-doa-label | 1.00 | 0.0004864 | 0.00036235 | 0.00095348 | 0.0013319 | 1.00 | 0.0016946 | 1 |

#### Gated rescue aggregate

| bank | hard cases | easy negative | fd-negative | trigger rate | hard trigger | hard rescue | easy damage | fd-negative damage | hit rate @0.002 deg | hard hit rate | max (deg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `disabled` | 5 | 90 | 5 | 0 | 0 | 0 | 0 | 0 | 0.90 | 0.20 | 0.0054541 |
| strict safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.4 | 0 | 0 | 0.92 | 0.60 | 0.0054541 |
| family-safe-adopt | 5 | 90 | 5 | 0.04 | 0.8 | 0.8 | 0 | 0 | 0.94 | 1.00 | 0.0044189 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `numOrdinaryAngleMiss=4`，ordinary miss 真实存在，但不是全部 tail。
- ordinary miss 的 default non-ref coherence median 为 `1`，所以 `coherence-v1` hard-collapse gate 对它天然无效。
- final-centered small polish 与 default 完全一致，说明它不是 final 点附近小扰动问题。
- `wide-basin-entry` 4/4 救回 ordinary miss，且几乎达到 truth-DoA oracle 水平。
- `SS-MF CP-U in-tooth` 相比 `SS-SF static` 有明显收益，说明多帧 CP baseline 本身健康。

### 5.2 仍未解决或反向的现象

- `single-MF` 对 ordinary miss damage rate 达到 `0.75`，不能作为 ordinary 默认 candidate。
- blanket wide 虽然改善 ordinary miss，但在全体 100 seed 上 RMSE / P99 / max 变差。
- 本轮只有 4 个 ordinary miss，不能据此得出稳定 gate 结论；后续 500-repeat ordinary-wide 已更适合代表 gate 稳定性。

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 343 | ordinary angle miss | default angle `0.0020406 deg`，truth-DoA oracle `0.00040927 deg`，tooth `0`，fdRef error `-0.012933 Hz`，fdRate error `14.234 Hz/s`，coherence floor `1`。 | ordinary miss 的典型症状：fd healthy、same-tooth、coherence 不 collapse。 |
| 256 | hard-collapse | default angle `0.0054541 deg`，coherence floor `0.00067049`，family-safe selected angle `0.00034933 deg`。 | 证明 hard-collapse 线仍成立，但不是 ordinary gate 证据。 |
| 345 | hard-collapse | default angle `0.003759 deg`，coherence floor `0.10043`，family-safe selected angle `7.5606e-05 deg`。 | 继续支持 hard-collapse family-safe bank。 |
| 346 | fd-negative | default angle `0.0022914 deg`，fdRef error `7.5262 Hz`，coherence floor `0.98368`。 | 说明 fd-not-healthy 不应混入 ordinary DoA rescue。 |

## 6. 机制解释

### 6.1 当前解释

本 replay 把 controlled in-tooth DoA tail 拆成三类：hard-collapse、ordinary non-collapse、fd-not-healthy。hard-collapse 由 non-ref coherence collapse 暴露，适合由 gated `wide + single-MF` basin-entry 处理；ordinary miss 不 collapse，fd 也健康，说明它更像 same-tooth 内的 DoA basin / adoption 问题。final-small-polish 无效，说明 final 点附近的小盒不足；wide center 有效，说明需要 basin-entry，而不是 final-centered 微调。

### 6.2 这个结果支持什么

- 支持新增 ordinary-wide gate 诊断 replay。
- 支持把 hard-collapse 与 ordinary non-collapse 分开处理。
- 支持继续把 final-small-polish 与 single-MF ordinary 默认路径旁路。

### 6.3 这个结果不证明什么

- 不证明 ordinary-wide gate 已经稳定；本轮只有 4 个 ordinary miss。
- 不证明 blanket wide 可用。
- 不证明 family-safe hard-collapse gate 能解决 ordinary miss。
- 不证明 estimator 默认路径或 full-flow 已修复。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉 final-small-polish、wide 或 single-MF。 |
| flow 默认路径 | 不改；该 replay 只推动后续 ordinary-wide gate 诊断。 |
| regression | 不写；只是机制 triage。 |
| replay / scan 下一步 | 已由 `replayMfInToothOrdinaryWideGateDiagnose` 扩大到 500-repeat 验证。 |
| 论文图 / 论文口径 | 可作为 outlier / bad-basin 机制说明；不作为 paper-facing 主精度结果。 |
| 排障记录 | 只摘取 ordinary miss 分类和 hard-collapse / ordinary 分流结论。 |

## 8. 限制与禁止解释

- 不要把 4/4 wide rescue 外推为稳定 policy。
- 不要继续通过放宽 hard-collapse coherence threshold 来解决 ordinary miss。
- 不要把 `single-MF` 当成 ordinary miss 默认 candidate。
- 不要把 final-small-polish 接入默认路径。
- 不要把 truth label、truth DoA 或 truth `fdRef/fdRate` 放入 runtime selector、gate、candidate adoption 或 final winner。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothOrdinaryAngleMissDiagnose.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 replay 的主要作用是拆分 ordinary non-collapse miss，并直接促成后续 `replayMfInToothOrdinaryWideGateDiagnose`。
- 后续 500-repeat ordinary-wide 结果已经将 gate 结论从“可能 preferred”降级为 `hold / partial diagnostic`，因此本文档不再作为 ordinary-wide gate 稳定性的最终证据。
