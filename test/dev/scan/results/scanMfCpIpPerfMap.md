# scanMfCpIpPerfMap 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `stress-test` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfCpIpPerfMap_20260428-195242.mat` |
| 当前一句话结论 | 当前 full-flow CP/IP performance map 被 wrong-tooth 与 same-tooth bad basin 明显污染；IP 在表观 angle / `fdRef` RMSE 上更好，不能解释为 IP 理论上优于 CP。 |
| 论文图定位 | `not for paper main figure`；可作为 appendix / internal stress-test，用于说明 full-flow threshold 与 tooth-acquisition 风险。 |
| 决策影响 | 固定为 full-flow 负结果；不继续扩大 repeat；CP/IP 论文主图转向 controlled in-tooth 或 resolved/local 口径。 |
| 下一步动作 | 暂停扩跑；若未来 flow 已稳定，再重跑并把本 snapshot 标为 `superseded`。 |
| 禁止误用 | 不能把该 full-flow stress-test 负结果写成“CP 理论不如 IP”，也不能用它否定 continuous-phase 主模型。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfCpIpPerfMap.m`
- 结果文档：`test/dev/scan/results/scanMfCpIpPerfMap.md`
- scan 类型：`full-flow stress scan / CP-IP performance-map candidate rejection`
- 主要问题：当前 simple dynamic full flow 是否已经足够稳定，能够支撑 CP/IP、known/unknown-rate 的论文性能图。
- 扫描对象：`SNR`、联合帧数 `P`、`CP-K / CP-U / IP-K / IP-U`，每个 `(P,SNR)` 使用 3 个 seed。
- 不覆盖范围：不隔离 correct-tooth / in-tooth 条件；不做 CRB-normalized resolved comparison；不验证单独的 CP/IP 理论模型。
- truth 使用口径：truth 只用于 offline evaluation，包括 `toothIdx`、`toothResidualHz` 与 truth-tooth hit-rate；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：No，当前仅作为 full-flow stress / negative evidence。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `CP-K` | continuous-phase，known `fdRate`。 | Eval only | 观察当前 full-flow CP known-rate 路径能否停在正确 tooth / basin。 | 不能当作纯 CP 理论性能。 |
| `CP-U` | continuous-phase，unknown `fdRate` nuisance。 | Eval only | 观察 unknown-stage release 与 CP flow 是否稳定。 | 不能把 tail 直接解释为 EFIM 信息损失。 |
| `IP-K / IP-U` | independent-phase baseline，切断跨帧公共相位 tying。 | Eval only | 作为对比基线；当前实现中可能弱化 comb / basin 锁定。 | 不能写成 IP 是本文主物理模型。 |
| `truth-tooth hit` | `toothIdx==0 && abs(toothResidualHz)<=50 Hz`。 | Eval only | 衡量 full-flow 是否落在 center tooth 附近。 | 不是 runtime gate，也不是论文最终成功标准。 |
| `full-sample RMSE` | 全部 36 个 full-flow 样本统计。 | Eval only | 用于展示 stress-test tail 与 bad-basin 风险。 | 不能直接和 CRB 做 local asymptotic 对比。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：全部 `P x SNR x seed = 36` 个样本，作为 stress-test 统计。
- `resolved-sample`：未定义；本 snapshot 没有用于论文 CRB 对比的 resolved filtering。
- `outlier rate`：可由 wrong-tooth / residual fail 间接观察，但不是 paper-facing resolved/outlier 判据。
- `truth-tooth / oracle range`：truth 只用于 offline label。
- `stress-test`：是。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpPerfMap_20260428-195242.mat` | 2026-04-28 | `representative` | `baseSeed=253`；`seedList=[253,254,255]`；`SNR=[0,5,10,15] dB`；`P=[8,10,20]`；`T_f=1/750 s`；每格 `numRepeat=3`；truth-tooth 判据 `toothIdx==0 && abs(toothResidualHz)<=50 Hz`。 | 36/36 个 task 跑通，但 CP-K / CP-U 的 truth-tooth hit-rate 和 RMSE 被 wrong-tooth / bad-basin tail 拉坏；该结果拒绝作为 CP/IP 论文性能图。 | none |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = [253, 254, 255]`
- `numRepeat = 3` per `(P,SNR)` grid
- `snrDbList = [0, 5, 10, 15]`
- `frameCountList = [8, 10, 20]`
- `frameIntvlSec = 1/750`
- 关键 scan 轴：`P`、`SNR`、`phaseMode=CP/IP`、`fdRateMode=known/unknown`
- 关键 offline label：`toothIdx==0 && abs(toothResidualHz)<=50 Hz`
- checkpoint：disabled / not recorded；本 scan 直接保存轻量 `scanData`。
- snapshot 保存变量：`scanData`
- 运行时间：控制台记录约 `22 min 11 s`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`scanName`、`runKey`、`utcRun`、`config`、`perfTable`、`aggregateTable`、`repeatOutCell`、`plotData`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：运行日志出现 `84` 次 near-singular warning 和 `2` 次 singular warning，主要来自 unknown warm-anchor / local release 过程；未导致 hard fail，但说明该 full-flow stress scan 的 unknown-stage 局部条件数较差。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

全部 `P x SNR x repeat = 36` 个样本按 case 汇总。`truth-tooth hit` 使用 strict offline 判据：`toothIdx==0 && abs(toothResidualHz)<=50 Hz`。

| case | samples | truth-tooth hit | wrong-tooth / residual fail | angle RMSE (deg) | angle median (deg) | angle P95 (deg) | fdRef RMSE (Hz) | fdRef median (Hz) | fdRef P95 (Hz) | fdRate RMSE (Hz/s) | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `CP-K` | 36 | 0.333 | 0.667 | 0.004460 | 0.003567 | 0.008434 | 267 | 124 | 563 | 0 | known-rate 仍经常未停在 center tooth。 |
| `CP-U` | 36 | 0.222 | 0.778 | 0.007285 | 0.001214 | 0.014249 | 25800 | 1500 | 52500 | 714 | 中位数较小但远齿 tail 极重。 |
| `IP-K` | 36 | 0.722 | 0.278 | 0.002254 | 0.000508 | 0.004486 | 66.8 | 30.6 | 150 | 0 | 当前 full-flow 下表观更稳。 |
| `IP-U` | 36 | 0.667 | 0.333 | 0.001715 | 0.000495 | 0.002960 | 54.5 | 33.3 | 108 | 4660 | angle / fdRef 表观最好，但不能解释成物理主模型。 |

### 5.2 按扫描轴汇总

| P | case | truth-tooth hit | angle RMSE (deg) | fdRef RMSE (Hz) | fdRate RMSE (Hz/s) | 解释 |
|---:|---|---:|---:|---:|---:|---|
| 8 | `CP-K` | 0.250 | 0.004937 | 219 | 0 | CP-K 已受 wrong-tooth 影响。 |
| 8 | `CP-U` | 0.250 | 0.007983 | 5587 | 78.3 | CP-U tail 已出现。 |
| 8 | `IP-K` | 0.750 | 0.002577 | 45.6 | 0 | IP-K 表观稳定。 |
| 8 | `IP-U` | 0.667 | 0.002063 | 47.3 | 5629 | IP-U angle 小但 `fdRate` 不作为主指标。 |
| 10 | `CP-K` | 0.500 | 0.004123 | 148 | 0 | 仍不足以支撑论文图。 |
| 10 | `CP-U` | 0.333 | 0.007815 | 16233 | 1079 | 远齿 tail 加重。 |
| 10 | `IP-K` | 0.667 | 0.002471 | 56.9 | 0 | IP baseline 仍较稳。 |
| 10 | `IP-U` | 0.583 | 0.001933 | 61.1 | 4709 | 表观 RMSE 小。 |
| 20 | `CP-K` | 0.250 | 0.004277 | 380 | 0 | 长窗口没有自动修正 CP tooth。 |
| 20 | `CP-U` | 0.083 | 0.005867 | 41202 | 600 | CP-U 在长窗口下 tooth hit 最差。 |
| 20 | `IP-K` | 0.750 | 0.001581 | 89.9 | 0 | IP-K 表观继续改善。 |
| 20 | `IP-U` | 0.750 | 0.000913 | 54.2 | 3367 | 高 P 下 angle 最小。 |

代表高 SNR / 长窗口格点 `P=20, SNR=15 dB`：

| case | truth-tooth hit | angle RMSE (deg) | angle P95 (deg) | fdRef RMSE (Hz) | fdRate RMSE (Hz/s) | 解释 |
|---|---:|---:|---:|---:|---:|---|
| `CP-K` | 0.667 | 0.003223 | 0.004465 | 203 | 0 | 仍有 CP tooth / basin 风险。 |
| `CP-U` | 0.000 | 0.006823 | 0.010294 | 31235 | 1017 | 明显不能作为 CP-U 论文性能点。 |
| `IP-K` | 1.000 | 0.000269 | 0.000291 | 12.7 | 0 | full-flow 表观很好。 |
| `IP-U` | 1.000 | 0.000269 | 0.000291 | 13.6 | 1831 | 表观同样很好。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `Full-flow CP/IP angle RMSE` | `SNR` 或 `P` | full-sample angle RMSE | `CP-K / CP-U / IP-K / IP-U` | No | 当前会把 flow failure 误写成 phase-model 差异。 |
| `Truth-tooth hit-rate stress plot` | `P` 或 `SNR` | truth-tooth hit-rate | four cases | Appendix / diagnostic | 只作为 full-flow stress / tooth acquisition 风险说明。 |
| `fdRef tail plot` | `P` 或 `SNR` | `fdRef` RMSE / P95 | four cases | Diagnostic | 用于暴露 CP-U 远齿 tail，不直接和 CRB 比。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- CP-K / CP-U 的 truth-tooth hit-rate 分别只有 `0.333 / 0.222`，说明 full-flow CP 路径经常没有停在 center tooth。
- CP-U 的 angle median 为 `0.001214 deg`，但 RMSE 为 `0.007285 deg`；`fdRef` median 为 `1500 Hz`，但 RMSE 达 `25800 Hz`、P95 达 `52500 Hz`。这说明是远齿 / bad-basin tail，而不是平滑噪声退化。
- `P=20, SNR=15 dB` 下 IP-K / IP-U 的 truth-tooth hit-rate 为 `1.0`，而 CP-U 为 `0.0`。若直接画论文图，会得到错误叙事。
- SNR 提高不能自动修 CP-U：`P=10/20, SNR=15 dB` 下 CP-U 仍有大 `fdRef` tail。

### 6.2 反向、污染或未解决现象

- IP 在当前 full-flow 里看起来更好，但这更可能来自 independent-phase 弱化了当前 CP flow 的 comb / basin 锁定，而不是 IP 理论上更优。
- CP-K 已知 rate 仍有 tooth failure，说明问题不只在 unknown-rate nuisance release。
- 该 scan 没有 resolved/local filtering，因此不能拿 full-sample RMSE 对 CRB。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `P=20, SNR=10 dB, seed=254, CP-U` | far-tooth outlier | `fdRefAbsErrHz≈100443 Hz`，`toothIdx=-134` | 说明 CP-U 可落到非常远的 comb tooth。 |
| `P=20, SNR=5 dB, seed=254, CP-U` | far-tooth outlier | `fdRefAbsErrHz≈57000 Hz`，`toothIdx=76` | 说明 tail 不是单个格点偶然。 |
| `P=10, SNR=15 dB, seed=253, CP-U` | far-tooth outlier | `fdRefAbsErrHz≈50988 Hz`，`toothIdx=68` | 说明高 SNR 不能自动修正 tooth / basin。 |

## 7. 机制解释

### 7.1 当前解释

该 scan 测到的是完整 flow 结果，而不是纯 CP/IP 统计模型差异。CP 模型保留跨帧公共相位 tying，因此在正确局部 basin 内可提供 `fdRef` / coherence consistency；但 full-flow 若先落入 wrong tooth 或 same-tooth bad basin，CP objective 也会在该错误中心附近形成稳定支路。IP / relaxed baseline 切断跨帧公共相位后，可能弱化当前实现中的 comb 锁定，使 angle 表观更好，但这同时意味着它不再保留本文主模型的连续相位信息。

因此，本 snapshot 应作为 “full-flow 仍不稳定” 的负结果，而不是 CP/IP 理论优劣的证据。真正可用于论文 CP/IP trade-off 的结果应使用 controlled in-tooth 或 resolved/local 口径。

### 7.2 这个 scan 支持什么

- 支持当前 full-flow CP/IP performance map 被 wrong-tooth / same-tooth basin 污染。
- 支持 full-flow CP-U 的远齿 tail 是结构性风险，不能靠扩大 SNR 或 repeat 自然消失。
- 支持把 `scanMfCpIpInToothPerfMap` 作为 CP/IP 论文口径，而不是使用本 full-flow scan。
- 支持将 full-flow CP/IP 结果保留为 stress-test / engineering boundary。

### 7.3 这个 scan 不证明什么

- 不证明 IP 理论上优于 CP。
- 不证明 continuous-phase 主模型错误。
- 不证明 `CP-U` 的 EFIM 信息损失会导致这些远齿 outlier。
- 不证明 estimator 默认路径应该改成 IP。
- 不证明可以写 regression 契约。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 只暴露 full-flow stress failure，不触碰 estimator 主核。 |
| flow 默认路径 | 暂不改；若要修 flow，应回到 selected tooth / same-tooth rescue 的 replay 或 scan。 |
| replay 下一步 | 不因本 scan 单独新增 replay；已有 in-tooth / subset / tail replay 可解释该类污染。 |
| regression | 不写；full-flow stress 负结果不是稳定契约。 |
| 论文图 | 不作为主图；最多作为 appendix / engineering diagnostic。 |
| 排障记录 | 保留“full-flow CP/IP 被污染，controlled in-tooth 结果才可解释”的结论即可，不复制长表。 |

## 9. 限制与禁止解释

- 不要用该 full-flow stress-test 负结果直接否定 CP 理论模型。
- 不要把 IP 表观 RMSE 更小写成 IP 物理上更好。
- 不要把 truth-tooth hit-rate 迁移成 runtime selector / gate。
- 不要把 full-sample RMSE 与 CRB 做 resolved/local 对比。
- 不要继续扩大 repeat 试图把该 snapshot 变成论文主图；需要先隔离 tooth / basin 污染。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfCpIpPerfMap_20260428-195242.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfCpIpPerfMap.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：若后续脚本已按 template 修改，默认通知 / checkpoint 外壳可能与该历史 snapshot 不完全一致；复现结论以本节配置为准。

## 11. 历史备注

- 当前只绑定 `scanMfCpIpPerfMap_20260428-195242.mat` 作为代表性 full-flow stress snapshot。
- 早期 frame subset 越界与 `msKnownDoaHalfWidth` 缺字段属于脚本合规清理过程中的工程错误，不作为结果结论保留。
- 若未来 flow 稳定后重跑，应保留本 snapshot 作为 `superseded` negative baseline，并明确新结果如何消除了 wrong-tooth / same-tooth tail。
