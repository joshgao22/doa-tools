# scanSfStaticMleCrbConsistency 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `representative / paper-facing static anchor / CRB口径已更新` |
| 最新代表性 snapshot | `test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-163111.mat` |
| 当前一句话结论 | 最新 200-repeat scan 证明 DoA-only 主口径切到 pilot-model effective-gain CRB 后，`SS-SF-DoA` 与 `MS-SF-DoA` 均稳定在约 `1.07–1.08 × CRB`；旧 unit-gain CRB 下的 `MS-SF-DoA 1.16–1.18×` gap 已被解释为 second-sat DoA-only signal-scale 口径差，而不是 estimator failure。 |
| 论文图定位 | `main / appendix static anchor`；可作为 single-frame static SS/MS、DoA-only/static、unit-vs-pilot CRB 口径说明的 paper-facing 对照。 |
| 决策影响 | 固定本 snapshot 为当前 static representative anchor；DoA-only case 主归一化口径改为 `angleCrbDoaPilot*`；unit-gain DoA CRB 只保留为 ideal / diagnostic 对照。 |
| 下一步动作 | 若要出论文图，建议用更大 repeat 重跑同一新口径确认平滑曲线；代码开发上不继续调 `SS/MS-SF-DoA` estimator，转向 dynamic / MF 的 pilot-model CRB 与 resolved MLE-vs-CRB。 |
| 禁止误用 | 不要继续用旧 unit-gain DoA CRB 判定 `MS-SF-DoA` 失败；不要把 DoA-only pilot-model CRB 用作 matched static DoA-Doppler CRB；不要用本 static scan 证明 dynamic full-flow 已通过。 |

## 1. Scan 身份

- 脚本：`test/dev/scan/scanSfStaticMleCrbConsistency.m`
- 结果文档：`test/dev/scan/results/scanSfStaticMleCrbConsistency.md`
- scan 类型：`paper-facing curve / static anchor / SF MLE-vs-CRB consistency`
- 主要问题：single-frame static 层中，SS/MS、DoA-only/static 四个 canonical case 是否与**对应模型口径**的 CRB 同尺度，能否作为 dynamic scan 的健康锚点？
- 扫描对象：`SNR=-20:5:15 dB`，`numRepeat=200`，`SS-SF-DoA`、`MS-SF-DoA`、`SS-SF-Static`、`MS-SF-Static`。
- 不覆盖范围：不验证 multi-frame dynamic；不验证 CP/IP；不验证 unknown `fdRate` nuisance；不覆盖 tooth acquisition、subset selection、same-tooth rescue 或 flow-level bad-basin 机制。
- truth 使用口径：truth 用于 CRB 计算、误差统计、CRB-normalized tail 标记和 CRB-local 离线评价；不进入 runtime selector、gate、candidate adoption 或 final winner。
- 是否 paper-facing：Yes，作为 static anchor / appendix 或 main supporting figure 候选。

## 2. 术语与曲线口径

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `SS-SF-DoA` | 单星、单帧、DoA-only MLE。 | Eval only | 主口径应看 `angleCrbDoaPilotStdDeg` 与 `*OverDoaPilotCrb`。 | 不评估 `fdRef`，不能说明 Doppler 链健康。 |
| `MS-SF-DoA` | 多星、单帧、DoA-only MLE。 | Eval only | 主口径应看 pilot-model effective-gain DoA-only CRB；unit CRB 只作理想相干增益对照。 | 不要用旧 unit-gain CRB 的 `1.16–1.18×` gap 判定 estimator 失败。 |
| `SS-SF-Static` | 单星、单帧、联合 DoA-`fdRef` static MLE。 | Eval only | 主口径仍是 matched `crbPilotSfDoaDoppler`，即 `angleCrbStdDeg` / `fdRefCrbStdHz`。 | 不要改用 DoA-only pilot CRB 来评价 matched DoA-Doppler estimator。 |
| `MS-SF-Static` | 多星、单帧、联合 DoA-`fdRef` static MLE。 | Eval only | 用 matched static joint CRB 评价 static DoA-Doppler 链；pilot DoA CRB 只作角度信息口径对照。 | 低 SNR gap 不能直接外推为 dynamic 机制。 |
| `angleCrbDoaUnitStdDeg` | unit-gain array-only DoA CRB。 | Eval only | 理想化上限：假设每颗星 DoA-only pilot 信息完整相干。 | 不代表当前 DoA-only no-Doppler estimator 的实际 signal model。 |
| `angleCrbDoaPilotStdDeg` | pilot-model effective-gain DoA-only CRB。 | Eval only | 对 DoA-only estimator 的主 CRB 口径；包含 Doppler-shifted pilot 到 no-Doppler atom 的有效相干增益。 | 不等于 matched DoA-Doppler CRB。 |
| `angleCrbStdDeg` / `angleCrbDdStdDeg` | matched static DoA-Doppler CRB 投影到角度的标量标准差。 | Eval only | static DoA-Doppler case 的主 CRB 口径。 | 不应用来评价 DoA-only estimator 是否贴 CRB。 |
| `crbLocal*` | 在 resolved 样本上用 `trimNormCap=5` 的 CRB-normalized cap 剔除极端 tail 后的统计。 | Eval only | 用于 paper-facing MLE-vs-CRB 口径。 | 不是 unconditional efficiency；必须同时报告 keep / outlier rate。 |
| `outlierRate` | 被 CRB-local cap 标记为 tail 的比例。 | Eval only | 解释 low-SNR full-sample 与 CRB-local 差异。 | 不作为 runtime selector。 |

固定解释口径：

- DoA-only case 的主曲线使用 `crbLocalAngleRmseOverDoaPilotCrb`。
- Static DoA-Doppler case 的主曲线使用 `crbLocalAngleRmseOverCrb` 与 `crbLocalFdRefRmseOverCrb`。
- `angleCrbDoaUnit*` 保留为 ideal/unit-gain 对照，用于说明旧 gap 来源，不作为当前 DoA-only paper-facing 主 CRB。
- `truth` 只用于 offline 评价与 CRB 线性化，不迁移到任何 runtime 选择逻辑。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-163111.mat` | 2026-05-09 | `representative` | `baseSeed=253`，`seedList=253:452`，`numRepeat=200`，`SNR=-20:5:15 dB`，4 canonical case，checkpoint enabled，selected global sats `[1,2]` | DoA-only pilot CRB 后 `SS-SF-DoA≈1.072–1.074×`、`MS-SF-DoA≈1.068–1.084×`；static DoA-Doppler matched CRB 下 `SS-SF-Static≈1.094×`，`MS-SF-Static` 从 `1.197×` 收敛到 `1.079×`。 | 当前代表性结果；取代旧文档中仅按 unit DoA CRB 解释 `MS-SF-DoA` gap 的结论。 |
| `test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-001642.mat` | 2026-05-09 | `superseded / old CRB口径` | `numRepeat=2000`，`SNR=-20:3:13 dB`，旧 DoA-only 口径主要看 unit CRB。 | 当时 `MS-SF-DoA≈1.17× unit CRB`，`SS/MS-SF-Static` 健康；后续 replay 证明该 DoA-only gap 主要来自 second-sat signal-scale 口径差。 | 被 `20260509-163111` 的 pilot-model DoA-only CRB 口径覆盖；旧 2000-repeat 仍可作为 static matched CRB 历史锚点。 |

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = 253`
- `seedList = 253:452`
- `numRepeat = 200`
- `snrDbList = [-20, -15, -10, -5, 0, 5, 10, 15]`
- active methods：`SS-SF-DoA`, `MS-SF-DoA`, `SS-SF-Static`, `MS-SF-Static`
- selected satellites global index：`[1, 2]`
- `trimNormCap = 5`
- `staticMsHalfWidthDeg = [0.002, 0.002] deg`
- task count：`1600`
- checkpoint：enabled；run dir `tmp/scanSfStaticMleCrbConsistency/seed253_sat1_2_snr-20to15_n8_rep200_case4_cdd82ed7`；`1600/1600` complete；success 后 cleanup requested。
- snapshot 保存变量：`scanData`
- 运行时间：task grid 约 `9 min 54 s`；static CRB SNR sweep 约 `1.47 s`。

### 4.2 存档数据检查

- snapshot：`test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-163111.mat`
- `scanData` 主要字段仍包含：`config`, `truth`, `repeatTable`, `crbTable`, `perfTable`, `aggregateTable`, `crbLocalSummaryTable`, `topTailTable`, `topTailExportTable`, `checkpointSummaryTable`, `plotData`。
- 未保存大体量数据：未保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片。
- warning / fail 计数：命令行未显示导致 scan 中断的错误；checkpoint 正常完成并请求清理；snapshot 正常保存。

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

| case | samples / SNR | 主 CRB 口径 | crbLocal keep | angle RMSE / 主 CRB | fdRef RMSE / CRB | 备注 |
|---|---:|---|---:|---:|---:|---|
| `SS-SF-DoA` | 200 | pilot-model DoA-only CRB | `1.000` | `1.0719–1.0739` | N/A | 单星 DoA-only 已与 pilot-model CRB 对齐；unit CRB 下为 `1.0958–1.0979`。 |
| `MS-SF-DoA` | 200 | pilot-model DoA-only CRB | `1.000` | `1.0681–1.0835` | N/A | 旧 unit CRB 下仍为 `1.1645–1.1814`，但这是 ideal/unit-gain 对照，不再作为主判断。 |
| `SS-SF-Static` | 200 | matched static DoA-Doppler CRB | `1.000` | `1.0938–1.0954` | `0.9014–0.9179` | 单星 static angle / `fdRef` 为 CRB-level；`fdRef<1` 只作有限 MC / 条件统计解释。 |
| `MS-SF-Static` | 200 | matched static DoA-Doppler CRB | `1.000` | `1.1973 → 1.0787` | `1.1116 → 1.0515` | 低 SNR 有轻微 gap，高 SNR 收敛到 CRB-level。 |

### 5.2 按扫描轴汇总

#### DoA-only：pilot-model CRB 修正后的主切片

| SNR (dB) | `SS-SF-DoA` / pilot CRB | `SS-SF-DoA` / unit CRB | `MS-SF-DoA` / pilot CRB | `MS-SF-DoA` / unit CRB | 解释 |
|---:|---:|---:|---:|---:|---|
| -20 | 1.0719 | 1.0958 | 1.0681 | 1.1645 | pilot CRB 后 MS 不再显示旧 `1.16×` gap。 |
| -15 | 1.0723 | 1.0962 | 1.0685 | 1.1649 | SS/MS 同量级。 |
| -10 | 1.0728 | 1.0968 | 1.0692 | 1.1658 | unit gap 只保留为 ideal 对照。 |
| -5 | 1.0733 | 1.0973 | 1.0705 | 1.1672 | pilot CRB 后 MS 约 `1.07×`。 |
| 0 | 1.0736 | 1.0976 | 1.0735 | 1.1704 | SS/MS 几乎一致。 |
| 5 | 1.0737 | 1.0977 | 1.0781 | 1.1755 | 高 SNR MS 略升但仍为 CRB-level。 |
| 10 | 1.0738 | 1.0978 | 1.0818 | 1.1795 | pilot CRB 后约 `1.08×`。 |
| 15 | 1.0739 | 1.0979 | 1.0835 | 1.1814 | unit CRB 不应再用于判定失败。 |

#### Static DoA-Doppler：matched CRB 主切片

| SNR (dB) | `SS-SF-Static` angle/CRB | `SS-SF-Static` fdRef/CRB | `MS-SF-Static` angle/CRB | `MS-SF-Static` fdRef/CRB | 解释 |
|---:|---:|---:|---:|---:|---|
| -20 | 1.0954 | 0.9179 | 1.1973 | 1.1116 | MS 低 SNR 有轻微 angle/fd gap。 |
| -15 | 1.0940 | 0.9098 | 1.1889 | 1.1061 | 仍为低 SNR finite-sample / hard-pair 区。 |
| -10 | 1.0938 | 0.9057 | 1.1740 | 1.0991 | gap 开始收敛。 |
| -5 | 1.0938 | 0.9036 | 1.1511 | 1.0898 | MS static 逐渐接近 joint CRB。 |
| 0 | 1.0939 | 0.9025 | 1.1219 | 1.0785 | 中 SNR 已明显收敛。 |
| 5 | 1.0939 | 0.9019 | 1.0909 | 1.0614 | angle 与 fdRef 均为 CRB-level。 |
| 10 | 1.0940 | 0.9016 | 1.0786 | 1.0526 | 高 SNR 稳定接近 CRB。 |
| 15 | 1.0940 | 0.9014 | 1.0787 | 1.0515 | 高 SNR plateau 约 `1.08×` angle、`1.05×` fdRef。 |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| DoA-only pilot CRB consistency | SNR (dB) | `crbLocalAngleRmseOverDoaPilotCrb` | `SS-SF-DoA`, `MS-SF-DoA` | Yes / appendix | 主曲线；建议同时用浅线或附表给出 unit CRB 对照。 |
| DoA-only unit-vs-pilot gap | SNR (dB) | `*OverDoaUnitCrb` 与 `*OverDoaPilotCrb` | `MS-SF-DoA` | Mechanism / appendix | 用于解释旧 `MS-SF-DoA` gap，不作为主性能判据。 |
| Static angle RMSE/CRB | SNR (dB) | `crbLocalAngleRmseOverCrb` | `SS-SF-Static`, `MS-SF-Static` | Yes | matched static DoA-Doppler CRB；建议加 `y=1` 参考线。 |
| Static fdRef RMSE/CRB | SNR (dB) | `crbLocalFdRefRmseOverCrb` | `SS-SF-Static`, `MS-SF-Static` | Yes | `SS` 略低于 1 按 CRB-level 统计波动解读。 |
| keep / outlier rate | SNR (dB) | `crbLocalKeepRate`, `outlierRate` | 四个 case | Companion / appendix | 本次所有 case / SNR keep=1、outlier=0。 |
| Full-sample RMSE/CRB | SNR (dB) | `fullAngleRmseOver*`, `fullFdRefRmseOverCrb` | 四个 case | No / diagnostic | 本次 full 与 CRB-local 基本一致；仍不作为唯一论文成功标准。 |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `MS-SF-DoA` 在 pilot-model CRB 下全 SNR 为 `1.0681–1.0835×`，而 unit CRB 下仍为 `1.1645–1.1814×`。这正好复现 replay 结论：旧 gap 来自 DoA-only CRB signal-scale 口径，而不是 MLE / fusion failure。
- `SS-SF-DoA` 在 pilot-model CRB 下为 `1.0719–1.0739×`，和 `MS-SF-DoA` 同量级，说明 pilot-model DoA-only CRB 口径已经把 SS/MS 对齐。
- `SS-SF-Static` 在 matched static CRB 下 angle 稳定约 `1.094×`，`fdRef` 约 `0.90–0.92×`，说明单星 static DoA-Doppler 链健康。
- `MS-SF-Static` 在 matched static CRB 下从 `1.1973×` / `1.1116×` 收敛到 `1.0787×` / `1.0515×`，说明多星 static DoA-Doppler 链整体健康，只在低 SNR 有轻微 gap。
- 本次 `crbLocalKeepRate=1`、`resolvedRate=1`、`outlierRate=0` 覆盖所有 case / SNR，说明这些结论不是靠 trim 掩盖大量错误样本得到的。

### 6.2 反向、污染或未解决现象

- `MS-SF-Static` 低 SNR angle 仍有 `1.19×` 左右 gap；若论文中放低 SNR static 曲线，应解释为 finite-sample / hard-pair / low-SNR threshold 区，而不是完全 asymptotic efficient。
- `SS-SF-Static` 的 `fdRef` ratio 稳定低于 1；只能写成 finite MC / profile / scalar CRB 口径下的 CRB-level，不写成“超过 CRB”。
- `MS-SF-DoA` 在 unit CRB 下仍有 `1.16–1.18×`，这个 gap 现在应作为口径对照，不再驱动 estimator 修改。
- 本 scan 不涉及 dynamic CP phase tying、unknown `fdRate` nuisance、comb / tooth branch 或 same-tooth bad basin；这些 dynamic gap 不能由本结果直接解释完。

### 6.3 代表性异常格点 / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| 全 SNR / 全 case | no outlier | `crbLocalKeepRate=1`、`outlierRate=0`。 | 说明 pilot CRB 修正后的 DoA-only gap 不是 trim 产物。 |
| `-20 dB`, `MS-SF-Static` | low-SNR tail-like top rows | top-tail preview 中 angle / fdRef normalized error 可到 `3×` 左右，但仍在 CRB-local cap 内。 | 支持低 SNR gap 是有限样本 tail / hard-pair 效应，而不是 unresolved failure。 |
| `10–15 dB`, `MS-SF-Static` | high-SNR plateau | angle 约 `1.079×`，fdRef 约 `1.05×`。 | 支持 static 多星 matched CRB 已形成稳定健康锚点。 |

## 7. 机制解释

### 7.1 当前解释

这次 scan 的关键变化是：DoA-only estimator 不再用理想 unit-gain array-only CRB 作为主判据，而使用 pilot-model effective-gain DoA-only CRB。`MS-SF-DoA` 不估 Doppler，因此第二星的 Doppler-shifted pilot 在 no-Doppler DoA-only atom 上只有有限相干投影。旧 unit CRB 默认每颗星都以完整相干 pilot power 贡献 DoA 信息，因此对 second-sat 信息量过乐观。pilot-model CRB 把这部分 effective gain 折算进 `pwrSource` 后，SS/MS DoA-only 的 RMSE/CRB 立即对齐到 `~1.07×`。

Static DoA-Doppler case 则不同：matched `crbPilotSfDoaDoppler` 显式包含 Doppler-shifted pilot atom，第二星的 Doppler 不再是 DoA-only mismatch，因此 static DoA-Doppler 仍应使用 matched static CRB，而不是 pilot-model DoA-only CRB。`SS-SF-Static` 与 `MS-SF-Static` 的结果说明 static 地基是健康的。

因此，后续若 dynamic resolved / CRB-local scan 仍出现明显 MLE-vs-CRB gap，应优先把原因归到 dynamic-specific 机制：continuous-phase 跨帧 tying、`fdRate` nuisance、comb / tooth branch、same-tooth DoA bad basin、non-ref coherence 或 local optimizer floor，而不是重新怀疑 static DoA-only / static DoA-Doppler CRB 口径。

### 7.2 这个 scan 支持什么

- 支持 `crbPilotSfDoaOnlyEffective` 的 pilot-model effective-gain DoA-only CRB 作为 `SS/MS-SF-DoA` 的主比较口径。
- 支持旧 `MS-SF-DoA` unit-CRB gap 是 CRB signal-scale mismatch，而不是 estimator / truth-init / MS fusion failure。
- 支持 static 地基健康：`SS-SF-Static` 与 `MS-SF-Static` 已可作为 dynamic 结果的对照锚点。
- 支持停止继续优化 single-frame DoA-only estimator；后续精力应转向 dynamic in-tooth / MS / known-unknown / CP-IP。

### 7.3 这个 scan 不证明什么

- 不证明 dynamic full-flow 已贴 CRB。
- 不证明 multi-frame CP / IP 估计器没有 bad basin 或 tooth / comb tail。
- 不证明 unknown `fdRate` nuisance 的实际 MLE 已达到 EFIM。
- 不证明 pilot-model DoA-only CRB 可以替代 matched DoA-Doppler CRB。
- 不适合直接写成 regression 契约；它是 paper-facing static scan 结果，不是自动 pass/fail 护栏。

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；本结果说明当前 DoA-only MLE 不需要继续“救”。 |
| CRB helper | 已支持 `crbPilotSfDoaOnlyEffective`；后续 DoA-only scan 默认使用 pilot-model CRB。 |
| flow 默认路径 | 不改；与 dynamic flow / subset / rescue 无关。 |
| replay 下一步 | 不需要继续 `replaySfMsDoaCrbDiagnose` 机制扩展；保留其作为 CRB-scale audit 证据。 |
| regression | 已有 `regressionSfDoaOnlyPilotModelCrb` 锁住 effective-gain scaling；本 scan 本身不迁移为 regression。 |
| 论文图 | `main / appendix static anchor`；建议同时放 DoA-only pilot CRB 与 static matched CRB 两组曲线，unit CRB 作为附表或机制说明。 |
| 排障记录 | 主记录可补一句：DoA-only unit-CRB gap 已被 pilot-model effective CRB 解释，static 地基继续视为健康。 |

## 9. 限制与禁止解释

- 不要用本 scan 证明 dynamic estimator 默认路径已经通过。
- 不要把 `crbLocal*` 条件样本统计写成 unconditional estimator efficiency。
- 不要把 `SS-SF-Static` 的 `fdRef` ratio 略低于 1 解释为超过 CRB；这里只能写作 CRB-level / Monte Carlo 统计波动。
- 不要继续用 `angleCrbDoaUnitStdDeg` 作为 `MS-SF-DoA` 主评价口径；它只是 ideal/unit-gain 对照。
- 不要把 DoA-only pilot-model CRB 与 matched DoA-Doppler CRB 混用。
- 不要把本结果迁移到 runtime selector、gate、candidate adoption 或 final winner 逻辑。

## 10. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/scan/scanSfStaticMleCrbConsistency_20260509-163111.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/scanSfStaticMleCrbConsistency.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 11. 历史备注

- `scanSfStaticMleCrbConsistency_20260509-001642.mat` 的 2000-repeat 旧结果证明 static matched DoA-Doppler 锚点健康，但当时 DoA-only `MS-SF-DoA` 主要按 unit-gain CRB 解读，因此留下约 `1.17×` gap。
- `replaySfMsDoaCrbDiagnose_20260509-155057.mat` 先定位到 second-sat effective pilot gain / signal-scale 口径差；本次 scan 把该结论回填到 paper-facing static scan。
- 当前代表性结果改为 `scanSfStaticMleCrbConsistency_20260509-163111.mat`，因为它同时保留 unit-gain 与 pilot-model DoA-only CRB 字段，并把 DoA-only 主口径切到 pilot-model CRB。
