# scanMfCpIpTying 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only` |
| 最新代表性 snapshot | `test/data/cache/scan/scanMfCpIpTying_20260428-143341.mat` |
| 当前一句话结论 | continuous phase tying 在 truth center 上仍把 center tooth 作为最低点，但不能消除 `1/T_f` comb；若 final / seed 已落到 wrong tooth，CP objective 也会支持该错齿。 |
| 论文图定位 | `mechanism figure / appendix candidate`；用于解释 CP / relaxed / IP phase tying 与 comb tooth 结构。 |
| 决策影响 | 固定为 CP/IP tying 机制证据；不作为 SNR/P 性能图；不进入 regression。 |
| 下一步动作 | 暂停扩跑；如论文需要机制图，可基于该 snapshot 重画 folded tooth curve / alias table。 |
| 禁止误用 | 不能写成 CP 消除了 comb ambiguity；也不能把 relaxed/IP 的 line-scan 最小点解释为更好的 `fdRef` estimator。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanMfCpIpTying.m`
- 结果文档：`test/dev/scan/results/scanMfCpIpTying.md`
- scan 类型：`mechanism scan / objective line scan / CP-IP phase tying`
- 主要问题：`continuous / relaxed / independent` 三种 phase tying 对同一多帧样本的 `fdRef` one-dimensional objective comb 结构有何影响。
- 扫描对象：phase mode、CP-K / CP-U、truth / seed / final centers、`fdRef` line grid。
- 不覆盖范围：不做 Monte Carlo 性能统计；不比较 SNR / P 曲线；不验证 full-flow；不证明 estimator default 已修复。
- truth 使用口径：truth center 只用于机制线扫；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Appendix / mechanism figure only。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `continuous` | CP 模型，同一卫星跨帧共享公共相位 tying。 | No | 本文主物理模型的 phase tying。 | 不等于自动消除 comb。 |
| `relaxed` | 弱化 continuous tying 的诊断模式。 | No | 用于理解绝对相位 tying 被释放后的 objective 变化。 | 不作为正式 estimator 模型。 |
| `independent` | 每个观测块独立相位的 IP baseline。 | No | 对比切断跨帧公共相位后的 `fdRef` 可辨识性。 | 不代表更好的 Doppler 估计。 |
| `center` | `truth / staticSeed / cpKnownSeed / finalEstimate` 等线扫中心。 | Oracle for truth center | 区分正确中心与 final bad-basin 中心。 | 不同 center 的最小点不能混成同一性能指标。 |
| `minDeltaFd` | line scan 最小点相对 center 的 `fdRef` offset。 | Eval only | 判断 objective 最低 tooth 位置。 | 不能当作实际 estimator error，除非 center 是 actual estimate。 |
| `aliasIdx` | `minDeltaFd` 相对 `1/T_f=750 Hz` 的 tooth index。 | Eval only | 用于识别 center tooth 或 `±1` tooth。 | 不是 runtime selector。 |
| `centerDeltaObj` | center point 相对 line minimum 的 objective gap。 | Eval only | 判断 center tooth 是否是当前 line 最低点。 | 不是 RMSE / CRB 指标。 |

常见 scan 口径在本文件中的取值：

- `full-sample`：N/A；非 MC performance scan。
- `resolved-sample`：N/A。
- `outlier rate`：N/A。
- `truth-tooth / oracle range`：truth center 只用于机制线扫。
- `stress-test`：否；本 scan 是 objective mechanism scan。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanMfCpIpTying_20260428-143341.mat` | 2026-04-28 | `representative` | `numFrame=20`；`T_f=1/750 s`；`SNR=10 dB`；`phaseList=[continuous, relaxed, independent]`；`numGrid=801`；`numAliasSide=2`；CP-K / CP-U 均运行；只保存轻量 `scanData`。 | CP continuous 在 truth center 上 center tooth 最低；但 CP-K final 与 CP-U final 均落到 `-1 tooth` 附近； relaxed/IP 不是更好的 `fdRef` 估计。 | 覆盖旧 10-frame console-only 观察。 |
| 未绑定 snapshot | 2026-04-28 | `superseded` | 旧控制台日志；`saveSnapshot=false`；`numFrame=10`；同类 phase list 与 grid。 | 旧日志中 CP-U final 曾回到 center tooth。 | 被当前 20-frame snapshot 取代，不再作为代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = N/A`；固定机制样本。
- `numRepeat = N/A`；非 Monte Carlo。
- `snrDb = 10`
- `numFrame = 20`
- `frameIntvlSec = 1/750`
- `numPilotSample = 512`
- `refFrameIdx = 10`
- `phaseList = [continuous, relaxed, independent]`
- `numGrid = 801`
- `numAliasSide = 2`
- 关键 scan 轴：`phaseMode`、`center`、`fdRateMode=known/unknown`、`fdRef` line offset。
- 关键 resolved / outlier 判据：N/A；本 scan 不运行 MC estimator 统计。
- checkpoint：disabled / not recorded。
- snapshot 保存变量：`scanData`
- 运行时间：snapshot 未记录 `elapsedSec`。

### 4.2 存档数据检查

- 顶层变量：`data / meta / inventory`
- `data.scanData` 字段：`truth`、`scanConfig`、`centerRegistry`、`frameInfo`、`cpKnown`、`cpUnknown`、`scanName`、`utcRun`、`snapshotFile`
- `cpKnown / cpUnknown` 字段：`fdRateMode`、`phaseList`、`centerNameList`、`phaseScan`、`summaryTable`、`toothTable`
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
- warning / fail 计数：无 hard fail 记录。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

#### CP-K phase / center summary

| phase mode | center | minDeltaFd (Hz) | aliasIdx | centerDeltaObj | 备注 |
|---|---|---:|---:|---:|---|
| continuous | truth | 0.000000 | 0 | 0 | 正确中心下 CP 保住 center tooth。 |
| continuous | staticSeed | -348.750000 | 0 | 322793.798204 | seed center 与 truth 差异极大，line minimum 不在 center。 |
| continuous | finalEstimate | -750.000000 | -1 | 1.494901 | final center 已支持 `-1 tooth`。 |
| relaxed | truth | -251.250000 | 0 | 0.682014 | 放松 tying 后 center 不是最低点。 |
| relaxed | staticSeed | -457.500000 | -1 | 2.213423 | 弱 tying 下可偏向邻齿。 |
| relaxed | finalEstimate | -277.500000 | 0 | 0.817795 | 不代表更准确 `fdRef`。 |
| independent | truth | -251.250000 | 0 | 0.682014 | 与 relaxed 一致。 |
| independent | staticSeed | -457.500000 | -1 | 2.213423 | 与 relaxed 一致。 |
| independent | finalEstimate | -277.500000 | 0 | 0.817795 | 与 relaxed 一致。 |

#### CP-U phase / center summary

| phase mode | center | minDeltaFd (Hz) | aliasIdx | centerDeltaObj | 备注 |
|---|---|---:|---:|---:|---|
| continuous | truth | 0.000000 | 0 | 0 | unknown-rate 下 truth center 也保住 center tooth。 |
| continuous | cpKnownSeed | -750.000000 | -1 | 1.494901 | CP-K seed 已在 wrong tooth 附近。 |
| continuous | finalEstimate | -750.000000 | -1 | 4.077799 | CP-U final 没有救回 center tooth。 |
| relaxed | truth | -251.250000 | 0 | 0.682014 | relaxed / IP 一致。 |
| relaxed | cpKnownSeed | -277.500000 | 0 | 0.817795 | 弱 tying 下偏移不等于正确 Doppler。 |
| relaxed | finalEstimate | -630.000000 | -1 | 4.237186 | 仍接近 `-1 tooth`。 |
| independent | truth | -251.250000 | 0 | 0.682014 | 与 relaxed 一致。 |
| independent | cpKnownSeed | -277.500000 | 0 | 0.817795 | 与 relaxed 一致。 |
| independent | finalEstimate | -630.000000 | -1 | 4.237186 | 与 relaxed 一致。 |

### 5.2 按扫描轴汇总

CP continuous truth center 的 alias table：

| aliasIndex | deltaFdRef (Hz) | deltaObj | 解释 |
|---:|---:|---:|---|
| -2 | -1500 | 15.890611 | 次邻 wrong tooth 仍存在。 |
| -1 | -750 | 1.956850 | 最近 wrong tooth gap 有限。 |
| 0 | 0 | 0 | center tooth 为最低点。 |
| 1 | 750 | 10.020294 | 另一侧 wrong tooth gap 更大。 |
| 2 | 1500 | 32.016549 | 远 tooth 被压低。 |

CP-K continuous finalEstimate center 的 alias table：

| aliasIndex | deltaFdRef (Hz) | deltaObj | 解释 |
|---:|---:|---:|---|
| -2 | -1500 | 4.784291 | final wrong-tooth basin 周围的次邻齿。 |
| -1 | -750 | 0 | 当前 line minimum。 |
| 0 | 0 | 1.494901 | center tooth 反而高于 wrong tooth。 |
| 1 | 750 | 9.268793 | 另一侧 tooth。 |
| 2 | 1500 | 23.320762 | 远 tooth。 |

CP-U continuous finalEstimate center 的 alias table：

| aliasIndex | deltaFdRef (Hz) | deltaObj | 解释 |
|---:|---:|---:|---|
| -2 | -1500 | 7.891140 | final wrong-tooth basin 周围的次邻齿。 |
| -1 | -750 | 0 | 当前 line minimum。 |
| 0 | 0 | 4.077799 | center tooth 未被 release 救回。 |
| 1 | 750 | 20.124055 | 另一侧 tooth。 |
| 2 | 1500 | 48.136865 | 远 tooth。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `fdRef line around truth center` | `deltaFdRef` or folded tooth index | `deltaObj` | continuous / relaxed / independent | Appendix / mechanism figure | 主图候选；强调 CP truth center 最低，但 gap 有限。 |
| `fdRef line around final center` | `deltaFdRef` | `deltaObj` | CP-K / CP-U continuous | Diagnostic only | 用于说明 final 已在 wrong tooth；不作为性能图。 |
| `alias table` | `aliasIndex` | `deltaObj` | truth / final centers | Mechanism table | 适合正文文字或附录表。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `continuous | truth` 的最低点是 `deltaFdRef=0`，说明在正确中心附近 CP 确实保留了区分 center tooth 的信息。
- 但 `continuous | truth` 的 `-750 Hz` tooth 只比 center 高 `1.956850`，说明 CP 不会消除 comb，只是给出有限 gap。
- CP-K final 与 CP-U final 的 continuous line minimum 都是 `-750 Hz`，分别对应 `centerDeltaObj=1.494901` 与 `4.077799`，说明 final / seed 一旦在 wrong tooth，CP objective 会支持该错齿。
- relaxed 与 independent 在该 snapshot 中数值一致，说明它们都相当于切断或弱化跨帧 absolute phase tying 后的 IP-like baseline。
- independent truth center 的最小点是 `-251.25 Hz`，不是 `0 Hz`，说明 IP 不是更好的 `fdRef` 估计，而是 `fdRef` 可辨识性变弱。

### 6.2 反向、污染或未解决现象

- 旧 10-frame 日志曾出现 CP-U final 回到 center tooth，但没有绑定 snapshot，且被当前 20-frame 代表性结果取代。
- 该 scan 不包含 SNR / P / MC 性能曲线，因此不能直接替代 `scanMfCpIpInToothPerfMap` 或未来 MLE-vs-CRB scan。
- 该 scan 不解释 full-flow selector 为什么会落到 wrong tooth，只解释落到 wrong center 后 objective 的局部结构。

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `CP-K continuous finalEstimate` | wrong-tooth center | minimum at `-750 Hz`，`aliasIdx=-1` | 说明 CP final center 可被 wrong tooth 稳定支持。 |
| `CP-U continuous finalEstimate` | wrong-tooth center | minimum at `-750 Hz`，`centerDeltaObj=4.077799` | 说明 unknown-rate release 不保证救回 center tooth。 |
| `relaxed / independent truth` | IP-like ambiguity | minimum at `-251.25 Hz` | 说明切断 tying 后 `fdRef` 可辨识性下降。 |

## 7. 机制解释

### 7.1 当前解释

continuous phase tying 保留跨帧公共相位，因此在 truth center 附近，它能让 center tooth 成为 line minimum，并对邻近 alias tooth 给出有限 objective gap。但由于帧间隔仍带来 `1/T_f=750 Hz` 的 comb 结构，该 gap 不是无限大；当 flow / seed / final estimate 已经落入 wrong tooth 附近时，CP objective 会在该错误中心附近形成自洽局部支路。

relaxed 与 independent 模式放松或切断跨帧公共相位，允许更多 frame/block-wise complex phase 被 profile 掉。这样 comb 看起来可能更平滑，angle optimization 有时更容易，但代价是 `fdRef` 物理可辨识性下降。因此该结果支持 CP/IP 的机制区分，而不是支持把 IP 改成主模型。

### 7.2 这个 scan 支持什么

- 支持 CP continuous 在正确中心附近保留 center-tooth 信息。
- 支持 CP 不能消除 `1/T_f` comb ambiguity。
- 支持 full-flow CP wrong-tooth 结果不能直接归因于理论模型错误，而是 seed / flow / basin 问题。
- 支持 relaxed / independent 是 IP-like baseline，不是更准确的 `fdRef` estimator。

### 7.3 这个 scan 不证明什么

- 不证明 CP/IP 的 Monte Carlo 性能优劣。
- 不证明 CP-U release 能自动救回 wrong tooth。
- 不证明 estimator 默认 flow 已修复。
- 不证明 IP 应作为主模型。
- 不证明可以写 regression 契约。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；该 scan 只解释 objective 机制。 |
| flow 默认路径 | 不改；wrong-tooth 仍需 subset selection / flow replay 处理。 |
| replay 下一步 | 不需要新增；若排查 wrong-tooth，应回到 subset / periodic-vs-subset / in-tooth tail replay。 |
| regression | 不写；机制曲线不是自动契约。 |
| 论文图 | `appendix / mechanism figure`，可辅助解释 CP/IP tying 与 comb。 |
| 排障记录 | 保留“CP 保留 center-tooth 信息但不能消除 comb；final wrong center 可自洽”的结论即可。 |

## 9. 限制与禁止解释

- 不要写成 CP 消除了 comb ambiguity。
- 不要把 relaxed / independent 的 line minimum 当作更好的 `fdRef` 估计。
- 不要把该 objective line scan 当作 Monte Carlo performance map。
- 不要用旧 10-frame console-only 结果覆盖当前 snapshot 结论。
- 不要把 truth center 线扫中的 oracle 信息迁移到 runtime selector 或 gate。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanMfCpIpTying_20260428-143341.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanMfCpIpTying.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：当前脚本如已 template 化，通知 / snapshot 外壳可能和历史 snapshot 不完全一致；机制表以本节 snapshot 为准。

## 11. 历史备注

- 未绑定 snapshot 的旧 10-frame 日志只作为历史观察保留，不作为当前代表性结果。
- 当前 20-frame snapshot 更保守：CP-U final 未救回 center tooth，因此本文档按“CP-U release 不保证救齿”来解释。
