# <scanName> 结果记录模板

> 使用方式：复制本文件到 `test/dev/scan/results/<scanName>.md`，再替换尖括号占位符。本文档只记录 scan 结果、snapshot 绑定、曲线口径、机制解释和当前决策影响；不要把脚本使用手册、完整运行日志、README 内容或排障记录全文复制进来。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `<paper-facing / representative / diagnostic-only / stress-test / negative / hold / superseded / needs-rerun>` |
| 最新代表性 snapshot | `<test/data/cache/scan/<scanName>_yyyymmdd-HHMMSS.mat>` |
| 当前一句话结论 | `<用一句话说明本 scan 支持什么>` |
| 论文图定位 | `<main figure / appendix / mechanism figure / internal diagnostic / not for paper>` |
| 决策影响 | `<固定结果 / 需要重跑 / 暂缓 / 推进到 replay / 不进入 regression>` |
| 下一步动作 | `<补图、补 resolved 统计、重跑 confirm、更新主记录、暂停>` |
| 禁止误用 | `<例如：不能把 full-flow stress-test 负结果解释成 CP 理论不如 IP>` |

## 1. Scan 身份

- 脚本：`test/dev/scan/<scanName>.m`
- 结果文档：`test/dev/scan/results/<scanName>.md`
- scan 类型：`<paper-facing curve / mechanism scan / stress test / schedule scan / engineering smoke>`
- 主要问题：`<这个 scan 要回答的一个核心问题>`
- 扫描对象：`<SNR / P / T_f / blockLen / strategy / CP-IP / known-unknown / resolved condition>`
- 不覆盖范围：`<例如不验证 full-flow、不证明 estimator default、不比较 rescue 策略>`
- truth 使用口径：`<未使用 / 只用于 offline label / oracle-controlled / evaluation only>`
- 是否 paper-facing：`<Yes / No / Appendix only / 需要 resolved 口径后再说>`

## 2. 术语与曲线口径

> 本节必须保留。scan 结果最容易混淆的是“横轴、纵轴、case、resolved 条件和 truth 字段到底代表什么”。每个结果文档都应解释文中反复出现的字段、曲线和状态名，避免几周后只看到一串方法名却想不起含义。

| 名称 / 字段 | 含义 | 是否使用 truth | 如何解读 | 禁止解释 |
|---|---|---:|---|---|
| `<metric>` | `<1-2 行说明它是什么>` | `<No / Eval only / Oracle>` | `<它能说明什么>` | `<它不能说明什么>` |
| `<axis>` | `<扫描轴，例如 P, SNR, blockLen>` | No | `<横轴含义>` | `<不要跨不同 time-origin / search-domain 口径硬比较>` |
| `<case>` | `<CP-K / CP-U / IP-K / IP-U / strategy tag>` | No | `<模型、估计模式或策略>` | `<不要把 stress-test failure 解释成理论结论>` |

常见 scan 口径模板：

- `full-sample`：所有 repeat / task 的统计，用于说明 full-flow threshold、bad-basin 或 outlier 风险。
- `resolved-sample`：满足预定义 resolved 条件的样本统计，用于和 CRB / EFIM 对比。
- `outlier rate`：未进入 resolved / local regime 的比例，用于工程捕获能力说明，不作为主贡献。
- `truth-tooth / oracle range`：只用于受控上限或离线评价，不能进入 runtime selector、gate、candidate adoption 或 final winner。
- `stress-test`：用于暴露失败机制，不直接作为论文主性能图。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `<test/data/cache/scan/...mat>` | `<yyyy-mm-dd>` | `representative` | `<P、SNR、numRepeat、strategy、关键开关>` | `<本 snapshot 的一句话结论>` | `<none / 覆盖旧 snapshot>` |
| `<old snapshot>` | `<yyyy-mm-dd>` | `superseded` | `<旧配置>` | `<为什么被取代>` | `<被最新代表性结果覆盖>` |

维护规则：

- 当前代表性结果只保留一个 `representative`。
- 旧结果如果被新结果完全覆盖且结论一致，可以删除旧行；如果结论发生变化，必须保留变化原因。
- snapshot 未随仓库同步时，本节仍作为人工复现实验索引。

## 4. 最新代表性运行

### 4.1 配置

- `baseSeed = <...>`
- `seedList = <...>`
- `numRepeat = <...>`
- `snrDbList = <...>`
- `frameCountList = <...>`
- 关键 scan 轴：`<...>`
- 关键 resolved / outlier 判据：`<...>`
- checkpoint：`<disabled / completed N/N / cleaned / retained at tmp/...>`
- snapshot 保存变量：`scanData`
- 运行时间：`<...>`

### 4.2 存档数据检查

- 顶层变量：`<data / meta / inventory>`
- `scanData` 字段：`<config / taskTable / resultTable / aggregateTable / metricTable / plotData / checkpointSummaryTable>`
- 未保存大体量数据：`<rxSigCell / full sceneSeq / fixture cache / transition bundle / full objective map / images>`
- warning / fail 计数：`<若有，写 compact 数量和影响>`

## 5. 主要统计与曲线结果

### 5.1 主表 / 主切片

| 条件 / case | samples | primary metric | median | P95 | P99 | resolved rate | outlier rate | 备注 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `<case>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` |

### 5.2 按扫描轴汇总

| axis value | case | metric 1 | metric 2 | metric 3 | 解释 |
|---:|---|---:|---:|---:|---|
| `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` |

### 5.3 图形口径

| 图 | 横轴 | 纵轴 | 曲线 | 是否论文候选 | 注意事项 |
|---|---|---|---|---:|---|
| `<Figure A>` | `<P / SNR / blockLen / windowMs>` | `<RMSE / CRB inflation / aliasGap>` | `<case list>` | `<Yes / No>` | `<log-y / resolved only / full sample>` |

## 6. 可观察现象

### 6.1 支持当前结论的现象

- `<现象 1：写清楚来自哪张表、哪条曲线、哪个指标，以及为什么重要。>`
- `<现象 2>`

### 6.2 反向、污染或未解决现象

- `<例如 full-flow 被 wrong-tooth / same-tooth bad basin 污染。>`
- `<例如 candidate coverage 增加但 selected hit 没改善。>`
- `<例如 SNR 提高没有自动消除 tail。>`

### 6.3 代表性异常格点 / strategy / seed

| 条件 | 类型 | 现象 | 对结论的作用 |
|---|---|---|---|
| `<P=..., SNR=...>` | `<outlier / stress / warning / negative>` | `<关键数值>` | `<说明污染、tail、或机制边界>` |

## 7. 机制解释

### 7.1 当前解释

`<用 1-3 段解释为什么这些结果会出现。优先使用论文主线和当前排障机制语言，例如 continuous-phase tying、unknown-rate nuisance、EFIM loss、wrong-tooth / same-tooth bad basin、angle-vs-Doppler-consistency trade-off。不要只重复方法名。>`

### 7.2 这个 scan 支持什么

- `<例如：支持 known/unknown-rate 信息损失主要体现在 fdRef 或 DoA- fdRef coupling。>`
- `<例如：支持 controlled CP/IP 是 angle-vs-Doppler-consistency trade-off。>`
- `<例如：支持 block length 增强 alias tooth separation。>`

### 7.3 这个 scan 不证明什么

- `<例如：不证明 estimator 默认路径已修复。>`
- `<例如：不证明 full-flow CP/IP 结果可以作为论文主图。>`
- `<例如：不证明 structured schedule 可以进入默认 bank。>`
- `<例如：不证明可以写 regression 契约。>`

## 8. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | `<不改 / 无影响 / 需要另行验证>` |
| flow 默认路径 | `<不改 / 暂缓 / 进入 flow-like replay>` |
| replay 下一步 | `<需要 / 不需要 / 对应 replay 名称>` |
| regression | `<不写 / 只在契约稳定后写 / 已有护栏>` |
| 论文图 | `<main / appendix / stress test / diagnostic-only>` |
| 排障记录 | `<不更新 / 主记录摘一句 / 机制归并版补证据 / 历史归档追加 snapshot>` |

## 9. 限制与禁止解释

- 不要把 offline truth 评价字段用于 runtime selector、gate、candidate adoption 或 final winner。
- 不要用 full-flow stress-test 的负结果直接否定模型理论。
- 不要用单一 hit rate 代替 RMSE / P95 / P99 / resolved rate / outlier rate。
- 不要把 oracle-controlled 或 in-tooth scan 直接解释为 full-flow 已通过。
- 不要把机制 scan 直接迁移为 regression，除非已经形成稳定自动契约。
- `<按当前 scan 补充其他禁止解释。>`

## 10. 恢复与复现

```matlab
snapshotFile = '<test/data/cache/scan/<scanName>_yyyymmdd-HHMMSS.mat>';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

随后打开：

```text
`test/dev/scan/<scanName>.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 11. 历史备注

- `<只记录会改变结论或后续优先级的历史变化。不要复制完整旧表。>`
- `<例如：full-flow 结果被污染，因此转为 stress-test；controlled in-tooth 结果作为代表性 trade-off。>`
