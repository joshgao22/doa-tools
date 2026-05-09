# scanSfStaticMleCrbConsistency 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / paper-facing static anchor` |
| 最新代表性 snapshot | `test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-001642.mat` |
| 当前一句话结论 | 2000-repeat static SF scan 表明，`SS-SF-Static` 全 SNR 稳定贴近单星 static CRB，`MS-SF-Static` 在低 SNR 有轻微 threshold / finite-sample gap，但随 SNR 升高稳定收敛到 joint CRB，可作为 dynamic MLE-vs-CRB gap 的 static 健康锚点。 |
| 论文图定位 | `main / appendix static anchor`；可作为 SS/MS、DoA-only/static、static joint CRB consistency 的 paper-facing 对照。 |
| 决策影响 | 固定该 snapshot 作为当前 static representative anchor；不继续优化 static estimator；后续 dynamic 若在 resolved / CRB-local 口径仍有 gap，应优先排查 dynamic-specific 机制。 |
| 下一步动作 | 后续转向 dynamic in-tooth / MS / known-unknown / CP-IP；若要作图，优先画 `MS-SF-Static` 和 `SS-SF-Static` 的 angle / fdRef RMSE-over-CRB，并配套 outlier / keep rate。 |
| 禁止误用 | 不要把本 static 结果解释为 dynamic full-flow 已贴 CRB；不要把低 SNR full-sample tail 当成 static estimator 结构性错误；不要用 DoA-only case 的 CRB-normalized ratio 判断 static Doppler 链健康性。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanSfStaticMleCrbConsistency.m`
- 结果文档：`test/dev/scan/results/scanSfStaticMleCrbConsistency.md`
- scan 类型：`paper-facing curve / static anchor / SF MLE-vs-CRB consistency`
- 主要问题：single-frame static DoA-Doppler estimator 在 SS/MS、DoA-only/static 四个 canonical case 下是否与对应 CRB 同尺度，能否作为 dynamic scan 的健康锚点？
- 扫描对象：`SNR=-20:3:13 dB`，`numRepeat=2000`，`SS-SF-DoA`、`MS-SF-DoA`、`SS-SF-Static`、`MS-SF-Static`。
- 不覆盖范围：不验证 multi-frame dynamic；不验证 CP/IP；不验证 unknown `fdRate` nuisance；不覆盖 tooth acquisition、subset selection、same-tooth rescue 或 flow-level bad-basin 机制。
- truth 使用口径：truth 用于 CRB 计算、误差统计、CRB-normalized tail 标记和 CRB-local 离线评价；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Yes，作为 static anchor / appendix 或 main supporting figure 候选。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `SS-SF-DoA` | 单星、单帧、DoA-only MLE 与对应 DoA CRB。 | Eval only | 纯角度基线，检查 DoA-only 链路。 | 不评估 `fdRef`，不能说明 Doppler 链健康。 |
| `MS-SF-DoA` | 多星、单帧、DoA-only MLE 与对应 joint DoA CRB。 | Eval only | 纯角度多星基线。 | 不能直接与 static DoA-Doppler 的 CRB-normalized ratio 混作同一信息结构。 |
| `SS-SF-Static` | 单星、单帧、联合 DoA-`fdRef` static MLE 与单星 static CRB。 | Eval only | static Doppler 链和 angle / `fdRef` metric 的健康锚点。 | 不代表 multi-frame CP / IP 已通过。 |
| `MS-SF-Static` | 多星、单帧、联合 DoA-`fdRef` static MLE 与 joint static CRB。 | Eval only | multi-sat static anchor；中高 SNR 贴 joint CRB 说明 static 地基健康。 | 低 SNR gap 不能直接外推到 dynamic 机制。 |
| `full-sample` | 所有 finite estimator 输出的统计。 | Eval only | 暴露 low-SNR threshold / extreme tail 对 RMSE 的拉高。 | 不要求 full-sample 在最低 SNR 必须贴 CRB。 |
| `resolved` | 本 scan 中所有 case `resolvedRate=1`，等价于 finite resolved 输出。 | Eval only | 表示没有出现 solver-level unresolved 样本。 | 不能说明没有 CRB-normalized tail。 |
| `crbLocal*` | 在 resolved 样本上用 `trimNormCap=5` 的 CRB-normalized cap 剔除极少数 low-SNR tail 后的统计。 | Eval only | 用于 paper-facing MLE-vs-CRB 口径。 | 不是 unconditional efficiency；必须同时报告 keep / outlier rate。 |
| `crbLocalKeepRate` | CRB-local 保留率。 | Eval only | 低 SNR tail 是否影响 CRB-local 统计的比例指标。 | 不等同于实际 receiver 命中率。 |
| `outlierRate` | 被 CRB-local cap 标记为 tail 的比例。 | Eval only | 解释 low-SNR full-sample 与 CRB-local 差异。 | 不作为 runtime selector。 |

固定解释口径：

- static anchor 的主曲线优先使用 `crbLocalAngleRmseOverCrb` 和 `crbLocalFdRefRmseOverCrb`。
- `fullAngleRmseOverCrb` 在 `-20 dB` 会被极少数 angle tail 放大，不作为 static estimator 是否健康的主判断。
- `resolvedRate=1` 说明本 scan 没有 unresolved 失败；但 `crbLocalKeepRate` 仍需报告，用于区分低 SNR tail。
- `truth` 只用于 offline 评价与 CRB 线性化，不迁移到任何 runtime 选择逻辑。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-001642.mat` | 2026-05-09 | `representative` | `baseSeed=253`，`seedList=253:2252`，`numRepeat=2000`，`SNR=-20:3:13 dB`，4 canonical case，`trimNormCap=5`，`staticMsHalfWidthDeg=[0.002,0.002]`，selected global sats `[1,2]` | `SS-SF-Static` 全 SNR 约 `1.031–1.037 × CRB`；`MS-SF-Static` angle 从 `1.157 × CRB` 收敛到 `1.0209 × CRB`，`fdRef` 从 `1.037 × CRB` 收敛到 `1.0025 × CRB`。 | 当前唯一 static paper-facing anchor snapshot。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:2252`
- `numRepeat = 2000`
- `snrDbList = [-20, -17, -14, -11, -8, -5, -2, 1, 4, 7, 10, 13]`
- active methods：`SS-SF-DoA`, `MS-SF-DoA`, `SS-SF-Static`, `MS-SF-Static`
- selected satellites global index：`[1, 2]`
- `trimEnable = true`
- `trimNormCap = 5`
- `staticMsHalfWidthDeg = [0.002, 0.002] deg`
- `fdRangeBase = [-2e5, 2e5] Hz`
- task count：`24000`
- checkpoint：enabled；run dir `tmp/scanSfStaticMleCrbConsistency/seed253_sat1_2_snr-20to13_n12_rep2000_case4_fed8b80b`；`24000/24000` complete；success 后 cleanup requested。
- snapshot 保存变量：`scanData`
- 运行时间：task grid 约 `2 h 24 m 30 s`；总 elapsed 约 `2 h 39 m 47 s`。

### 4.2 存档数据检查

- 顶层 snapshot 内容：`data / meta / inventory`，其中 `data.scanData` 为主结果。
- `scanData` 主要字段：`scanName`, `runKey`, `utcRun`, `config`, `truth`, `repeatTable`, `crbTable`, `perfTable`, `aggregateTable`, `crbLocalSummaryTable`, `topTailTable`, `topTailExportTable`, `checkpointSummaryTable`, `checkpointCleanupReport`, `plotData`, `elapsedSec`。
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片。
- warning / fail 计数：命令行未显示导致 scan 中断的错误；checkpoint 正常完成并请求清理；snapshot 正常保存。
- 备注：本次运行 header 仍完整打印了 2000 个 seed，说明运行时尚未使用后续的 seed-list 压缩打印补丁；该问题只影响日志可读性，不影响 `scanData` 数值结果。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

| case | samples / SNR | resolvedRate | crbLocalKeepRate range | crbLocal angle RMSE/CRB | crbLocal fdRef RMSE/CRB | 备注 |
|---|---:|---:|---:|---:|---:|---|
| `SS-SF-DoA` | 2000 | `1.000` | `0.9975–1.0000` | `1.0515–1.0571` | N/A | DoA-only 单星稳定，`-20 dB` 有极少数 angle CRB tail。 |
| `MS-SF-DoA` | 2000 | `1.000` | `0.9995–1.0000` | `1.1663–1.1760` | N/A | DoA-only 多星相对 joint CRB 约 `1.17×`，用于对照，不作为 static Doppler 链判断。 |
| `SS-SF-Static` | 2000 | `1.000` | `0.9975–1.0000` | `1.0313–1.0373` | `0.98189–0.98552` | 单星 static angle / `fdRef` 全 SNR 贴近 CRB，可作为单星 static 地基。 |
| `MS-SF-Static` | 2000 | `1.000` | `0.9995–1.0000` | `1.0209–1.1570` | `1.0025–1.0370` | 多星 static 低 SNR 有轻微 angle gap，中高 SNR 收敛到 joint CRB。 |

### 5.2 按扫描轴汇总

#### Static DoA-Doppler 主切片

| SNR (dB) | `SS-SF-Static` angle/CRB | `SS-SF-Static` fdRef/CRB | `MS-SF-Static` angle/CRB | `MS-SF-Static` fdRef/CRB | keep / outlier 说明 |
|---:|---:|---:|---:|---:|---|
| -20 | 1.0373 | 0.98552 | 1.1570 | 1.0370 | `SS` keep `0.9975`、outlier `0.0025`；`MS` keep `0.9995`、outlier `0.0005`。 |
| -17 | 1.0344 | 0.98338 | 1.1486 | 1.0338 | keep `1.000`，outlier `0`。 |
| -14 | 1.0331 | 0.98259 | 1.1394 | 1.0310 | keep `1.000`，outlier `0`。 |
| -11 | 1.0324 | 0.98221 | 1.1282 | 1.0279 | keep `1.000`，outlier `0`。 |
| -8 | 1.0320 | 0.98203 | 1.1142 | 1.0243 | keep `1.000`，outlier `0`。 |
| -5 | 1.0317 | 0.98194 | 1.0976 | 1.0205 | keep `1.000`，outlier `0`。 |
| -2 | 1.0316 | 0.98191 | 1.0786 | 1.0163 | keep `1.000`，outlier `0`。 |
| 1 | 1.0315 | 0.98189 | 1.0583 | 1.0116 | keep `1.000`，outlier `0`。 |
| 4 | 1.0314 | 0.98189 | 1.0401 | 1.0069 | keep `1.000`，outlier `0`。 |
| 7 | 1.0314 | 0.98189 | 1.0272 | 1.0039 | keep `1.000`，outlier `0`。 |
| 10 | 1.0313 | 0.98189 | 1.0215 | 1.0027 | keep `1.000`，outlier `0`。 |
| 13 | 1.0313 | 0.98189 | 1.0209 | 1.0025 | keep `1.000`，outlier `0`。 |

#### DoA-only 对照

| case | angle RMSE/CRB range | 解释 |
|---|---:|---|
| `SS-SF-DoA` | `1.0515–1.0571` | 单星 DoA-only 贴近对应 CRB，说明基础 angle metric 正常。 |
| `MS-SF-DoA` | `1.1663–1.1760` | 多星 DoA-only 相对 joint CRB 存在稳定有限样本 / grid-local gap；不影响 static DoA-Doppler 链健康性判断。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| Static angle RMSE/CRB | SNR (dB) | `crbLocalAngleRmseOverCrb` | `SS-SF-Static`, `MS-SF-Static` | Yes | 主 static anchor 图；建议加 `y=1` 参考线。 |
| Static fdRef RMSE/CRB | SNR (dB) | `crbLocalFdRefRmseOverCrb` | `SS-SF-Static`, `MS-SF-Static` | Yes | `SS` 略低于 1 按有限样本 / 条件样本 CRB-level 解读，不写成超过 CRB。 |
| Static keep / outlier rate | SNR (dB) | `crbLocalKeepRate`, `outlierRate` | 四个 case | Companion / appendix | 说明 `-20 dB` 极少数 tail，其他 SNR 全保留。 |
| DoA-only angle RMSE/CRB | SNR (dB) | `crbLocalAngleRmseOverCrb` | `SS-SF-DoA`, `MS-SF-DoA` | Optional | 用作 DoA-only baseline；不要和 static-Doppler 信息结构硬比较。 |
| Full-sample RMSE/CRB | SNR (dB) | `fullAngleRmseOverCrb`, `fullFdRefRmseOverCrb` | 四个 case | No / diagnostic | `-20 dB` 被少数 angle tail 放大，只用于 threshold 解释。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `SS-SF-Static` 的 `crbLocalAngleRmseOverCrb` 在全 SNR 仅为 `1.0313–1.0373`，`crbLocalFdRefRmseOverCrb` 为 `0.98189–0.98552`，说明单星 static DoA-Doppler estimator 与 CRB metric 已经对齐。
- `MS-SF-Static` 的 angle ratio 随 SNR 从 `1.1570` 单调收敛到 `1.0209`，`fdRef` ratio 从 `1.0370` 收敛到 `1.0025`，说明多星 static joint CRB 口径在中高 SNR 下可被 MLE 接近。
- 除 `-20 dB` 外，四个 case 的 `crbLocalKeepRate=1`、`outlierRate=0`；`-20 dB` 也只有 `SS` 的 `0.25%` 和 `MS` 的 `0.05%` 极少数 CRB-normalized tail。
- `resolvedRate=1` 覆盖所有 SNR / case，说明本次 static scan 没有 solver-level unresolved failure，低 SNR full-sample gap 主要来自 CRB-normalized tail 而不是流程失败。
- `SS-SF-DoA` 约 `1.05× CRB`，说明基础 angle error 与 CRB metric 没有明显错配；static-Doppler 的 angle ratio 更低，说明引入 `fdRef` 联合估计没有破坏 angle 链路。

### 6.2 反向、污染或未解决现象

- `MS-SF-Static` 在 `-20 dB` 的 angle ratio 仍为 `1.1570`，低 SNR 不应写成完全 efficient；更合适解释为 threshold / finite-sample / hard-pair trade-off。
- `SS-SF-Static` 的 `fdRef` ratio 稳定低于 1，幅度约 `0.982`；这应按 Monte Carlo 条件样本、有限 repeat 和 CRB-level 统计波动解释，不应写成 estimator 违反 CRB。
- `MS-SF-DoA` DoA-only 相对 joint CRB 约 `1.17×`，提示纯 DoA-only 多星对照未完全贴 joint DoA CRB；但 static DoA-Doppler 主链在中高 SNR 已贴 CRB，因此不把它作为 static 地基回退证据。
- 本 scan 不涉及 dynamic CP phase tying、unknown `fdRate` nuisance、comb / tooth branch 或 same-tooth bad basin；这些 dynamic gap 不能由本结果直接解释完。

### 6.3 代表性异常格点 / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `-20 dB`, `SS-SF-DoA` / `SS-SF-Static` | angle CRB tail | top seed 包括 `1975`, `671`, `1611`, `867`, `617`，angle normalized error 可达几十倍 CRB。 | 解释 full-sample `-20 dB` angle RMSE/CRB 被少数 tail 放大。 |
| `-20 dB`, `MS-SF-Static` | single mixed tail | top seed `1611` 同时出现 angle / fd CRB tail；其余 top rows 多为 `crb-local`。 | 说明 MS static low-SNR full-sample gap 主要受极少数 tail 影响，CRB-local keep rate 仍为 `0.9995`。 |
| `SNR >= -17 dB` | no tail | top-tail export 中后续 rows 均为 `crb-local` preview，summary 中 outlier 为 0。 | 支持中高 SNR static anchor 稳定。 |

## 7. 机制解释

### 7.1 当前解释

本 scan 是 single-frame static anchor。`SS-SF-Static` 全 SNR 稳定贴近 CRB，说明 reference-sat Doppler state、angle error metric、static CRB 线性化和 single-frame static estimator 主链已经对齐。`MS-SF-Static` 在低 SNR 有轻微 angle gap，但随 SNR 升高快速收敛到 joint CRB，说明多星 static 融合层也已经恢复健康。

因此，后续若 dynamic resolved / CRB-local scan 仍出现明显 MLE-vs-CRB gap，应优先把原因归到 dynamic-specific 机制：continuous-phase 跨帧 tying、`fdRate` nuisance、comb / tooth branch、same-tooth DoA bad basin、non-ref coherence 或 local optimizer floor，而不是重新怀疑 static geometry、reference-sat 语义或 static CRB metric。

### 7.2 这个 scan 支持什么

- 支持 static 地基健康：`SS-SF-Static` 与 `MS-SF-Static` 已可作为 dynamic 结果的对照锚点。
- 支持 `4154 + 1165` 这组非平凡双星 pair 在 single-frame static 下能形成 joint CRB-level 对照。
- 支持后续 paper-facing 结果中将 static 与 dynamic 分层解释：static anchor 贴 CRB，dynamic gap 需要另由 dynamic 机制解释。
- 支持停止继续优化 static estimator；后续精力应转向 dynamic in-tooth / MS / known-unknown / CP-IP。

### 7.3 这个 scan 不证明什么

- 不证明 dynamic full-flow 已贴 CRB。
- 不证明 multi-frame CP / IP 估计器没有 bad basin 或 tooth / comb tail。
- 不证明 unknown `fdRate` nuisance 的实际 MLE 已达到 EFIM。
- 不证明 full-sample low-SNR RMSE 必须贴 CRB。
- 不适合迁移为 regression 契约；它是 paper-facing static scan 结果，不是自动 pass/fail 护栏。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；本结果只记录现有 SF static estimator 的 MC 统计。 |
| flow 默认路径 | 不改；与 dynamic flow / subset / rescue 无关。 |
| replay 下一步 | 不需要继续 static replay；dynamic 若 `crbLocal*` 不贴，再回 targeted dynamic replay。 |
| regression | 不写；static 地基已有其它 regression / quick suite 护栏，本 scan 只作为结果锚点。 |
| 论文图 | `main / appendix static anchor`；建议作为 MLE-vs-CRB consistency 的 static 对照图或文字表。 |
| 排障记录 | 可在主记录 / 机制归并版摘一句：“`scanSfStaticMleCrbConsistency_20260509-001642.mat` 固定为 static anchor，static 不再作为主战场。” |

## 9. 限制与禁止解释

- 不要用本 scan 证明 dynamic estimator 默认路径已经通过。
- 不要把 `crbLocal*` 条件样本统计写成 unconditional estimator efficiency。
- 不要把 `SS-SF-Static` 的 `fdRef` ratio 略低于 1 解释为超过 CRB；这里只能写作 CRB-level / Monte Carlo 统计波动。
- 不要用 DoA-only `MS-SF-DoA` 的 `1.17× CRB` 直接否定 multi-sat static DoA-Doppler，因为二者信息结构和 parameter coupling 不同。
- 不要把 `-20 dB` full-sample tail 当成中高 SNR static 锚点失败。
- 不要把本结果迁移到 runtime selector、gate、candidate adoption 或 final winner 逻辑。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-001642.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanSfStaticMleCrbConsistency.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。注意：复现实验时应确认当前 `printMfScanHeader` 已应用 seed-list 压缩补丁，否则仅日志头部会再次打印完整 seed list。

## 11. 历史备注

- 该 scan 是 `doaDopplerStatDualSatUraEciPerf` 的标准 scan 化结果；用途是固定 static paper-facing anchor，而不是新增 static rescue 或改变 estimator。
- 旧 static perf 经验显示中高 SNR `MS-SF-Static` 可贴 joint CRB；本次 2000-repeat 结果用标准 `scanData` / checkpoint / snapshot / results 文档口径固化这一结论。
- 本结果与 dynamic 排障记录中的当前优先级一致：static 只作为锚点和护栏，后续主线转向 dynamic resolved / CRB-local MLE-vs-CRB、known/unknown 信息损失与 controlled CP/IP trade-off。
