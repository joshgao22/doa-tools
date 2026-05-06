# <replayName> 结果记录模板

> 使用方式：复制本文件到 `test/dev/replay/results/<scriptName>.md`，再替换尖括号占位符。本文档只记录 replay 结果、snapshot 绑定、机制解释和当前决策影响；不要把脚本使用手册、完整运行日志或排障记录全文复制进来。

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `<representative / superseded / diagnostic-only / hold / candidate / blocked>` |
| 最新代表性 snapshot | `<test/data/cache/replay/<scriptName>_yyyymmdd-HHMMSS.mat>` |
| 当前一句话结论 | `<用一句话说明这轮结果支持什么>` |
| 决策影响 | `<不影响默认路径 / 推进到 flow-like replay / 保留为 diagnostic / 可作为论文机制图候选>` |
| 下一步动作 | `<重跑、更换 repeat、迁移到 scan、更新主记录、暂缓>` |
| 禁止误用 | `<例如：不能作为 estimator default 已修复证据；不能写 regression；不能用 truth label 做 runtime gate>` |

## 1. Replay 身份

- 脚本：`test/dev/replay/<scriptName>.m`
- 结果文档：`test/dev/replay/results/<scriptName>.md`
- replay 类型：`<机制可视化 / 小 MC / oracle-controlled / flow-like stress / 工程链路 smoke>`
- 主要问题：`<这个 replay 要回答的一个核心问题>`
- 观察范围：`<例如 controlled in-tooth、真实 subset-periodic flow、固定 hard seed、surface/ridge 切片>`
- 不覆盖范围：`<例如不验证 full-flow、不改变 estimator 默认路径、不证明 CRB consistency>`
- truth 使用口径：`<未使用 / 只用于 offline label / 只用于 oracle 上限 / 只用于结果评价，不进入 runtime selection>`

## 2. 机制词典与方法地图

> 本节必须保留。每个结果文档都应解释本 replay 中反复出现的方法名、标签和现象名，避免几周后只看到一串名字但想不起含义。

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `<method/tag>` | `<1-2 行说明它是什么>` | `<No / Evaluation only / Oracle only>` | `<选齿 / DoA basin center / fdRef/fdRate box / adoption / 只改评价标签>` | `<它能说明什么，不能说明什么>` |
| `<toothIdx>` | `fdRef` 相对 truth tooth 的 `1/T_f` 周期编号。 | Evaluation only | 只用于评价和结果标注。 | 不能进入 runtime selector / gate / final winner。 |
| `<representative example>` | `<按当前 replay 补充>` | `<...>` | `<...>` | `<...>` |

常见解释模板：

- `periodic wide`：完整等间隔帧上的较宽搜索，主要暴露 full-data baseline 是否会落入 wrong-tooth 或坏 basin。
- `selected subset`：从非周期 subset bank 中按 objective 选出的候选，主要承担 tooth selection，不等价于最终角度细化。
- `final periodic refine`：回到完整周期帧，在 selected subset 附近做同齿 refine；它通常改善 angle / Doppler consistency，但不应被期待自动修 wrong-tooth。
- `truth-centered oracle`：只用于上限诊断，不能进入 runtime flow、selector、gate 或 regression 契约。
- `wide / single-MF basin-entry`：改变 DoA basin center 的候选 family；若只在 replay 中有效，不能直接下沉 estimator 主核。
- `damage`：候选或 gate 使原本 easy / fd-healthy / baseline-hit 样本变差；damage 指标优先级高于单纯 rescue rate。

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `<test/data/cache/replay/...mat>` | `<yyyy-mm-dd>` | `representative` | `<snrDb、numRepeat、baseSeed、关键开关>` | `<本 snapshot 的一句话结论>` | `<none / 覆盖某旧 snapshot>` |
| `<old snapshot>` | `<yyyy-mm-dd>` | `superseded` | `<旧配置>` | `<为什么被取代>` | `<被最新 representative 覆盖>` |

维护规则：

- 当前代表性结果只保留一个 `representative`。
- 旧结果若被新结果完全覆盖且结论一致，可以删除旧行；若结论发生变化，保留 `superseded` 行并说明变化原因。
- snapshot 未随仓库同步时，本节仍作为人工排障记录。

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = <...>`
- `baseSeed = <...>`
- `numRepeat = <...>`
- seed range：`<...>`
- 关键场景 / frame / tooth 设置：`<...>`
- 关键 gate / policy / oracle 设置：`<...>`
- checkpoint：`<disabled / completed N/N / cleaned / retained at tmp/...>`
- snapshot 保存变量：`replayData`
- 运行时间：`<...>`

### 4.2 主要统计

> 只放 compact table。长表、per-repeat 细节和大候选表留在 `replayData` 中。若必须展示长表，只保留影响当前结论的列。

| 指标 | 数值 | 解释 |
|---|---:|---|
| `<primaryMetric>` | `<...>` | `<为什么它是主指标>` |
| `<tailMetric>` | `<...>` | `<P95/P99/max/outlier 的含义>` |
| `<damageMetric>` | `<...>` | `<easy / fd-negative / baseline-hit damage>` |
| `<runtimeMetric>` | `<...>` | `<checkpoint / warning / elapsed>` |

### 4.3 关键对比表

| 方法 / 阶段 | role | RMSE | median | P95 | P99 | hit / resolved rate | damage | 备注 |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `<baseline>` | `<baseline role>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` |
| `<candidate>` | `<candidate role>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` | `<...>` |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- `<现象 1：写清楚对应哪张表、哪个指标、为什么重要>`
- `<现象 2>`

### 5.2 仍未解决或反向的现象

- `<现象 1：例如 coverage 不足、P99 变差、easy damage、warning 集中在 tail>`
- `<现象 2>`

### 5.3 代表性 seed / case

| seed / case | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| `<seed>` | `<hard / easy / fd-negative / gate-miss / warning>` | `<关键数值或标签>` | `<说明 gate 有效、误伤、漏救、或只是诊断>` |

## 6. 机制解释

### 6.1 当前解释

`<用 1-3 段解释为什么这些结果会出现。优先用“选齿、same-tooth DoA basin、non-ref coherence collapse、CP/IP phase tying、unknown-rate nuisance”等机制语言，而不是只重复方法名。>`

### 6.2 这个结果支持什么

- `<例如：支持 subset 主要负责选齿>`
- `<例如：支持 family-safe-adopt 作为 controlled in-tooth replay-level candidate>`
- `<例如：支持 ordinary-wide 只作为 partial diagnostic>`

### 6.3 这个结果不证明什么

- `<例如：不证明 estimator 默认路径已修复>`
- `<例如：不证明 full-flow 性能已经可用于论文主图>`
- `<例如：不证明可以写 regression 契约>`
- `<例如：不证明可以把 truth-aware oracle 逻辑进入 runtime>`

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | `<不改 / 暂不下沉 / 需要单独验证后再说>` |
| flow 默认路径 | `<不改 / 推进到 flow-like replay / 暂缓>` |
| regression | `<不写 / 只在契约稳定后补 / 已有护栏名称>` |
| replay / scan 下一步 | `<下一步入口>` |
| 论文图 / 论文口径 | `<paper-facing / appendix / stress test / diagnostic-only>` |
| 排障记录 | `<不更新 / 主记录摘一句 / 机制归并版补证据 / 历史归档追加 snapshot>` |

## 8. 限制与禁止解释

- 不要把本 replay 的 replay-level candidate 直接当作 estimator default。
- 不要把 truth label、truth tooth、truth DoA、truth `fdRef/fdRate` 放入 runtime selector、gate、candidate adoption 或 final winner。
- 不要用单一 hit rate 代替 RMSE / P95 / P99 / outlier / damage 的联合判断。
- 不要把 controlled / oracle replay 的正结果直接解释为 full-flow 已通过。
- `<按当前 replay 补充其他禁止解释>`

## 9. 恢复与复现

```matlab
snapshotFile = '<test/data/cache/replay/<scriptName>_yyyymmdd-HHMMSS.mat>';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/<scriptName>.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- `<只记录会改变结论或后续优先级的历史变化。不要复制完整旧表。>`
- `<例如：100-repeat 曾是 preferred-candidate；500-repeat 后降级为 hold / partial diagnostic。>`
