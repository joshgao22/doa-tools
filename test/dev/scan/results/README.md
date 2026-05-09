# scan results 说明

本目录记录 scan 的详细运行结果、曲线 / 曲面观察、summary table、snapshot 绑定、机制解释和论文图候选口径。`test/dev/scan/README.md` 只保留 scan 脚本规范、入口索引和当前观察目标；具体统计表、长现象描述、历史 snapshot 索引和 scan-level 决策影响放在这里。

## 与 cache snapshot 的关系

- 大 `.mat` 文件放在 `test/data/cache/scan/`。
- 本目录不保存 `.mat`，只记录 snapshot 文件名、关键配置、主要 surface / table 结论、当前解释和论文图定位。
- 如果 snapshot 文件未随仓库同步，本目录中的摘要仍作为人工排障和论文图候选记录。
- 若本地有 snapshot，可恢复 `scanData` 后运行对应 scan 的 `Summary output and plotting` section 重出表格和图。
- scan snapshot 应保持轻量：保存 `scanData`，不要保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片。

## 结果文档模板

新增或重写单个 scan 结果文档时，先复制：

```text
test/dev/scan/results/template/scanResultTemplate.md
```

到：

```text
test/dev/scan/results/<scriptName>.md
```

再替换占位内容。模板的目标不是把每个结果文档写成长文，而是保证读者第一屏就能恢复上下文：

- 这个 scan 问什么；
- 最新代表性 snapshot 是哪个；
- 当前结论是什么；
- 横轴、纵轴、case、resolved / outlier 条件和 truth 字段分别是什么意思；
- 对 estimator / flow / replay / regression / 论文图有什么影响；
- 不能用这个结果证明什么。

旧结果文档不需要一次性全部迁移。只有在新增 snapshot、重写代表性结论、或当前文档已经难以快速恢复上下文时，再按模板顺手整理。

## 单个结果文档结构

每个 scan 默认一个结果文档：

```text
results/<scriptName>.md
```

推荐结构：

1. `状态摘要`：当前状态、最新 snapshot、一句话结论、论文图定位、决策影响、下一步动作和禁止误用。
2. `Scan 身份`：对应脚本、scan 类型、主要问题、扫描对象、不覆盖范围、truth 使用口径和 paper-facing 状态。
3. `术语与曲线口径`：解释横轴、纵轴、case、metric、resolved / outlier 条件、truth 字段和 stress-test 口径。
4. `Snapshot index`：代表性 snapshot 与 superseded snapshot 的配置、结论和覆盖关系。
5. `最新代表性运行`：关键配置、scanData 字段、checkpoint / warning / snapshot 信息。
6. `主要统计与曲线结果`：主表 / 主切片、扫描轴汇总和图形口径。
7. `可观察现象`：支持当前结论的现象、反向 / 污染 / 未解决现象、代表性异常格点 / strategy / seed。
8. `机制解释`：说明现象背后的 CP/IP phase tying、known/unknown-rate 信息损失、wrong-tooth / same-tooth basin、schedule adoption 或其他机制。
9. `对主流程的影响`：是否影响 estimator 默认路径、flow 默认路径、replay、regression、论文图和排障记录。
10. `限制与禁止解释`：明确哪些结论不能从该 scan 推出，尤其是 truth evaluation、oracle-controlled 和 full-flow stress-test 的边界。
11. `恢复与复现`：`loadExpSnapshot` 命令，以及恢复后运行哪个 section。
12. `历史备注`：只保留会改变结论或优先级的旧结果，不复制完整旧表。

如果某个 scan 有多个 sweep 或大量 snapshot，可以升级为：

```text
results/<scriptName>/README.md
results/<scriptName>/<caseName>.md
```

## 按 scan 类型填写的重点

### paper-facing curve / CRB scan

适用于 `scanMfKnownUnknownInformationLoss`、`scanMfRegimeMapByWindow` 等。

- 重点写清横轴、纵轴、主切片和图形口径。
- 若用于 CRB / EFIM 或 MLE-vs-CRB 对比，必须区分 `full-sample`、`loose/core resolved-sample`、`trimmed-sample`、`resolved/keep rate` 和 `outlier rate`。
- 结论要服务论文主线：regime identification、continuous-phase 必要性、unknown-rate 信息损失或 SS/MS/SF/MF 层级比较。
- 不要把单一 hit rate 写成论文主成功标准。

### controlled / oracle scan

适用于 `scanMfCpIpInToothPerfMap` 这类受控 in-tooth 或 oracle-range scan。

- 必须写清 truth 使用口径：truth 只用于 oracle 上限、离线评价或受控边界，不进入 runtime selector、gate、candidate adoption 或 final winner。
- 重点区分 upper bound、controlled mechanism、paper-facing trade-off 和 full-flow 真实性能。
- 正结果不能直接解释为默认 flow 已通过。

### full-flow stress scan

适用于 `scanMfCpIpPerfMap` 等暴露 wrong-tooth、same-tooth bad basin 或 threshold/outlier 的完整流程 scan。

- 必须写清当前结果是否被 flow 污染。
- full-flow 负结果可作为 stress test 或 outlier 机制说明，但不能直接否定 CP / IP 理论模型。
- 若要转为论文性能图，必须先有稳定的 selected-tooth、same-tooth rescue 或 resolved filtering 口径。

### mechanism / surface scan

适用于 `scanMfBlockLength`、`scanMfCpIpTying` 等。

- 重点写清图或 surface 的坐标轴、固定量、扫描量和 objective / residual 口径。
- 机制趋势可以作为论文机制图或 appendix 候选，但不等价于 regression 契约。
- 必须写“这个 scan 不证明什么”，避免把局部 surface、line scan 或 block-level 机制误迁移到 estimator 默认路径。

### schedule / subset bank scan

适用于 `scanMfSubsetBankCoverage` 等。

- truth 只用于 offline evaluation，不进入 runtime selector、gate、candidate adoption 或 final winner。
- 先看 selected 结果，再看 candidate coverage；candidate coverage 增加但 selected hit 不改善时，应优先检查 candidate-to-final adoption，而不是继续堆 schedule。
- 必须同时报告命中、误伤、tail、候选成本和 runtime / checkpoint。
- 负结果要保留，避免重复回到已证伪或收益不足的 schedule 方向。

## 当前结果文档索引

| scan | 结果文档 | 建议状态 | 当前结论 |
|---|---|---|---|
| `scanMfRegimeMapByWindow.m` | `scanMfRegimeMapByWindow.md` | `paper-facing` | 多帧窗口内 DoA 可近似静态，而 Doppler 漂移和静态 Doppler 相位误差已不可忽略；支持论文 regime justification。 |
| `scanMfKnownUnknownInformationLoss.m` | `scanMfKnownUnknownInformationLoss.md` | `paper-facing` | known / unknown Doppler-rate 的 CRB / EFIM 信息损失 scan，是后续 unknown-rate nuisance 退化解释的主入口。 |
| `scanMfCpIpInToothPerfMap.m` | `scanMfCpIpInToothPerfMap.md` | `paper-facing / representative` | controlled in-tooth 条件下，CP 稳定保住 center tooth、`fdRef` consistency 和 non-ref coherence；IP 在 angle RMSE 上可能更小，应解释为 angle-vs-Doppler-consistency trade-off。 |
| `scanSfStaticMleCrbConsistency.m` | `scanSfStaticMleCrbConsistency.md` | `paper-facing / representative` | Single-frame static SS/MS DoA-only 与 static DoA-Doppler MLE/CRB anchor；DoA-only 主口径已切到 pilot-model effective-gain CRB，旧 `MS-SF-DoA` unit-CRB gap 已解释为 signal-scale 口径差，static DoA-Doppler matched CRB 仍健康。 |
| `scanMfMleCrbInToothConsistency.m` | `scanMfMleCrbInToothConsistency.md` | `paper-facing / needs-rerun` | Doppler-aided / in-tooth local MLE 与 CRB 的一致性 scan；已升级为默认 auto MF initializer、loose/core/trimmed、MSE/CRB、top-tail 和 angle/fdRef 并列图，需按新口径重跑。 |
| `scanMfCpIpPerfMap.m` | `scanMfCpIpPerfMap.md` | `stress-test` | full-flow CP/IP 结果被 wrong-tooth / same-tooth bad basin 污染，当前不能直接作为 CP/IP 论文性能图。 |
| `scanMfBlockLength.m` | `scanMfBlockLength.md` | `diagnostic-only / appendix candidate` | 更长 pilot block 显著增强 `fdRef` comb / alias tooth separation，但块内 Doppler-rate 二次相位仍很小。 |
| `scanMfCpIpTying.m` | `scanMfCpIpTying.md` | `diagnostic-only` | 用于解释 CP / relaxed / IP phase tying 对 `fdRef` comb tooth 结构的影响，不直接作为性能图。 |
| `scanMfSubsetBankCoverage.m` | `scanMfSubsetBankCoverage.md` | `negative / diagnostic-only` | structured schedule quick screen 是负结果；candidate coverage 增加没有转化为 selected hit，后续应优先看 candidate-to-final adoption。 |

## 维护规则

- 新 snapshot 先追加到对应结果文档的 Snapshot index。
- 当前代表性结果只保留一个；旧结果若被完全覆盖且结论一致，可以删除旧行；若结论发生变化，标记为 `superseded` 并说明变化原因。
- scan 结果用于机制证据、参数曲线、论文图候选和负结果归档；只有稳定成自动契约后才迁移到 regression。
- paper-facing scan 应优先报告 RMSE / P95 / P99、CRB-normalized error、resolved rate 和 outlier rate；hit-rate 只作为工程捕获指标，不单独定义论文成功标准。
- full-flow stress-test 与 controlled / oracle scan 必须分开解释；不要用受控结果证明默认 flow 已修复，也不要用 full-flow 负结果否定理论模型。
- 结果文档中的 truth label、truth tooth、truth DoA、truth `fdRef/fdRate` 只能作为 oracle / evaluation / offline label；不能迁移到 runtime selector、gate、candidate adoption 或 final winner。
- README 只保留入口、结构和维护规则；长表、详细观察、snapshot 变更和具体结论放在单个 results 文档中。
