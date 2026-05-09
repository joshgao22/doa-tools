# replay results 说明

本目录记录 replay 的详细运行结果、观察现象、机制解释和 snapshot 绑定。`test/dev/replay/README.md` 只保留 replay 入口、运行规范和当前一句话结论；具体统计表、长现象描述、代表性 seed、历史 snapshot 索引和 replay-level 决策影响放在这里。

## 与 cache snapshot 的关系

- 大 `.mat` 文件放在 `test/data/cache/replay/`。
- 本目录不保存 `.mat`，只记录 snapshot 文件名、关键配置、主要统计、观察现象和结论。
- 如果 snapshot 文件未随仓库同步，本目录中的摘要仍作为人工排障记录。
- 若本地有 snapshot，可恢复 `replayData` 后运行对应 replay 的 `Summary output and plotting` section 重出表格和图。

## 结果文档模板

新增或重写单个 replay 结果文档时，先复制：

```text
test/dev/replay/results/template/replayTemplate.md
```

到：

```text
test/dev/replay/results/<scriptName>.md
```

再替换占位内容。模板的目标不是把每个结果文档写成长文，而是保证读者第一屏就能恢复上下文：

- 这个 replay 问什么；
- 最新代表性 snapshot 是哪个；
- 当前结论是什么；
- 对 estimator / flow / regression / 论文图有什么影响；
- 不能用这个结果证明什么；
- 文中反复出现的方法名和现象名分别是什么意思。

旧结果文档不需要一次性全部迁移。只有在新增 snapshot、重写代表性结论、或当前文档已经难以快速恢复上下文时，再按模板顺手整理。

## 单个结果文档结构

每个 replay 默认一个结果文档：

```text
results/<scriptName>.md
```

推荐结构：

1. `状态摘要`：当前状态、最新 snapshot、一句话结论、决策影响、下一步动作和禁止误用。
2. `Replay 身份`：对应脚本、结果类型、主要问题、观察范围、不覆盖范围和 truth 使用口径。
3. `机制词典与方法地图`：解释方法名、标签和现象名，例如 `periodic wide`、`selected subset`、`family-safe-adopt`、`ordinary-wide gate`、`fd-negative`、`non-ref coherence collapse`。
4. `Snapshot index`：代表性 snapshot 与 superseded snapshot 的配置、结论和覆盖关系。
5. `最新代表性运行`：关键配置、compact 统计表、主要对比表、checkpoint / warning / snapshot 信息。
6. `可观察现象`：支持当前结论的现象、仍未解决的现象、代表性 seed / case。
7. `机制解释`：说明现象背后的 tooth selection、same-tooth DoA basin、coherence collapse、CP/IP 或 unknown-rate 机制。
8. `对主流程的影响`：是否影响 estimator 默认路径、flow 默认路径、regression、后续 replay / scan、论文图和排障记录。
9. `限制与禁止解释`：明确哪些结论不能从该 replay 推出，尤其是 replay-level candidate 不等于 default path。
10. `恢复与复现`：`loadExpSnapshot` 命令，以及恢复后运行哪个 section。
11. `历史备注`：只保留会改变结论或优先级的旧结果，不复制完整旧表。

如果某个 replay 的结果数量太多，可以升级为：

```text
results/<scriptName>/README.md
results/<scriptName>/<caseName>.md
```

## 按 replay 类型填写的重点

### 机制可视化 replay

适用于 `replayMfCombToothSurface`、`replayMfInToothDoaDopplerRidgeTrace` 等。

- 重点写清图或 surface 的坐标轴、固定量、扫描量和 objective / residual 口径。
- 必须写 `这个结果不证明什么`，例如 ridge / surface 只能解释 coupling，不等价于可以把 slope 或 joint offset 下沉默认路径。
- 若作为论文机制图候选，写清它是正文、附录还是 diagnostic-only。

### 小 MC replay

适用于 `replayMfPeriodicVsSubsetToothSelect`、`replayMfRandomRescueEffectiveness`、`replayMfInToothGatedRescueEffectiveness`、`replayMfFlowLikeGatedBasinEntryEffectiveness` 等。

- 重点写 primary metric、tail metric、damage metric 和 runtime / checkpoint 口径。
- 必须同时报告改善与误伤，不能只写 rescue rate 或 hit rate。
- 若结论只是 replay-level candidate，必须说明是否需要 flow-like replay、scan 或更大 repeat 才能推进。

### oracle / controlled replay

适用于 `replayMfInToothFdRangeOracle`、`replayMfInToothTailCaseDiagnose`、`replayMfInToothDoaFdRangeEnvelope`、`replayMfInToothOrdinaryAngleMissDiagnose`、`replayMfInToothOrdinaryWideGateDiagnose` 等。

- 必须写清 truth 使用口径：truth 只用于 oracle 上限、offline label 或结果评价，不能进入 runtime selector / gate / adoption / winner。
- 重点区分 upper bound、mechanism triage、gate candidate 和 flow/default candidate。
- 若结果降级为 `hold`、`partial diagnostic` 或 `blocked`，应在状态摘要和限制解释中同时写清。

## 当前代表性结果索引

| replay | 结果文档 | representative snapshot | 当前结论 |
|---|---|---|---|
| `replayMfCombToothSurface.m` | `replayMfCombToothSurface.md` | `test/data/cache/replay/replayMfCombToothSurface_20260508-141941.mat` | CP-K 的 DoA 局部信息未退化，主要问题是 `fdRef` wrong-tooth；CP-U 另有 DoA-`fdRef`-`fdRate` coupling / bad basin。 |
| `replayMfPeriodicVsSubsetToothSelect.m` | `replayMfPeriodicVsSubsetToothSelect.md` | `test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat` | subset 负责选齿，periodic 主要负责同齿 refine；warning 只作为诊断元数据保留。 |
| `replayMfRandomRescueEffectiveness.m` | `replayMfRandomRescueEffectiveness.md` | `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat` | rescue/random bank 明显提升 central tooth 命中率，但已有 easy-case damage，后续应转 scan/gate 验证而非 blanket 常驻。 |
| `replayMfInToothFdRangeOracle.m` | `replayMfInToothFdRangeOracle.md` | `test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat` | half-tooth oracle 下频率链健康，但少数 same-tooth non-ref coherence tail 拉坏 RMSE。 |
| `replayMfInToothTailCaseDiagnose.m` | `replayMfInToothTailCaseDiagnose.md` | `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat` | gated `wide+single-MF` rescue bank 能全救 fixed hard-collapse seeds，并避开 easy / fd-not-healthy 负样本。 |
| `replayMfInToothDoaDopplerRidgeTrace.m` | `replayMfInToothDoaDopplerRidgeTrace.md` | `test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat` | final-centered ridge 证明 DoA-Doppler coupling 存在，但 minimum 常贴边界；优先推进 gated `wide+single-MF` basin-entry，而不是公式化 ridge slope。 |
| `replayMfInToothGatedRescueEffectiveness.m` | `replayMfInToothGatedRescueEffectiveness.md` | `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` | family-safe-adopt 显著压低 controlled in-tooth MS-MF DoA tail，hard rescue rate 提升到 `0.8`，easy / fd-negative damage 为 0；仍是 replay-level candidate。 |
| `replayMfInToothDoaFdRangeEnvelope.m` | `replayMfInToothDoaFdRangeEnvelope.md` | `test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat` | base-seed envelope 支持 wide / single-MF DoA basin-entry 与 `0.006 deg` family-safe step，不支持默认引入非零 fdRef / fdRate joint bank。 |
| `replayMfInToothOrdinaryAngleMissDiagnose.m` | `replayMfInToothOrdinaryAngleMissDiagnose.md` | `test/data/cache/replay/replayMfInToothOrdinaryAngleMissDiagnose_20260505-171659.mat` | ordinary miss 不是 final-small-polish 问题，主要由 wide-basin-entry 救回且接近 truth-DoA oracle；后续应找 truth-free ordinary-wide gate，不 blanket 打开 wide。 |
| `replayMfInToothOrdinaryWideGateDiagnose.m` | `replayMfInToothOrdinaryWideGateDiagnose.md` | `test/data/cache/replay/replayMfInToothOrdinaryWideGateDiagnose_20260505-233931.mat` | 500-repeat 扩大验证将 ordinary-wide 从 preferred candidate 降级为 hold / partial diagnostic；`objGain >= 10` 且 `0.001 <= wide DoA step <= 0.004 deg` 低误伤但只救约 59% ordinary miss，不推进默认 flow。 |
| `replayMfFlowLikeGatedBasinEntryEffectiveness.m` | `replayMfFlowLikeGatedBasinEntryEffectiveness.md` | `test/data/cache/replay/replayMfFlowLikeGatedBasinEntryEffectiveness_20260505-160714.mat` | 100-repeat flow-like 结果显示 basin-entry family 有救力，但当前 gate/adoption 有 easy damage 且 P95/max 变差；作为 negative / diagnostic，不推进默认 flow。 |
| `replaySfMsDoaCrbDiagnose.m`（已移除） | `replaySfMsDoaCrbDiagnose.md` | `test/data/cache/replay/replaySfMsDoaCrbDiagnose_20260509-155057.mat` | 历史 CRB-scale audit：`MS-SF-DoA` 的 unit-CRB gap 已定位为 second-sat signal-scale CRB 口径差；当前 active 入口合并到 `replaySfStaticDualSatDiagnose.m` 的 compact `doaOnlyCrbScaleTable`。 |

## 维护规则

- 新 snapshot 先追加到对应结果文档的 Snapshot index。
- 当前代表性结果只保留一个；若新结果完全覆盖旧结果且结论一致，可删除旧 snapshot 记录，不必保留 `superseded` 行。
- 结果文档记录证据、长表、机制解释和 snapshot 绑定；排障主记录只摘取影响当前优先级的结论。
- replay-level candidate 进入默认 flow 或 regression 前，必须另有 flow-like replay、scan 或稳定契约支撑；不要只因为某个 results md 写了 `candidate` 就修改 estimator 默认路径。
- 结果文档中的 truth label、truth tooth、truth DoA、truth `fdRef/fdRate` 只能作为 oracle / evaluation / offline label；不能迁移到 runtime selector、gate、candidate adoption 或 final winner。
