# strategy 目录说明

`test/dev/strategy/` 用于策略比较和决策记录。这里的脚本回答“哪种流程规则更值得采用”或“哪些方向不要再重复投入”，不作为 regression 护栏。

## 当前策略文件

### `compareMfSubsetPeriodicRefineStrategies.m`

- 比较内容：subset 选齿、periodic refine、random rescue 等 flow 策略。
- 主要问题：哪种两阶段策略更稳，哪些 case 仍错。
- 输出：策略对比表、summary。

### `compareMfUnknownDoaReleaseEffect.m`

- 比较内容：unknown-rate 分支 DoA release 效果。
- 主要问题：放开 DoA 是否改善 CP-U，而不是只 release fdRate。
- 输出：before/after objective、angle、fd summary。

### `probeMfFrameSubsetStrategy.m`

- 比较内容：frame subset 策略探针。
- 主要问题：不同 subset schedule 是否改善 tooth selection。
- 输出：subset/tooth summary。
- 当前定位：历史策略探针，不作为当前 curated bank 系统扫描入口；新的 bank coverage / random rescue 比较集中到 `test/dev/scan/scanMfSubsetBankCoverage.m`。

### `probeMfSubsetToothSelectThenPeriodicRefine.m`

- 比较内容：subset 选齿后再 periodic refine。
- 主要问题：先选 tooth 再同齿 refine 是否优于单一路径。
- 输出：tooth/refine 对比。

### `sweepMfInToothDoaRelease.m`

- 比较内容：in-tooth DoA release 范围。
- 主要问题：same-tooth hard case 是否需要 very-small DoA release。
- 输出：release width sweep。

### `sweepMfUnknownDoaBox.m`

- 比较内容：unknown DoA box 大小。
- 主要问题：CP-U 是否因 DoA box 过窄而退化。
- 输出：box sweep summary。

### `blockedDynamicDirections.md`

- 内容：已证伪或不建议重复投入的方向。
- 主要问题：哪些路线不要再回流主路径。
- 定位：策略黑名单 + 重试准入条件，不是 README，不是完整排障记录。

## results 目录

- `results/` 可选记录较长的策略对比结果、snapshot 绑定和详细表格。
- 短结论、已证伪方向和重试准入条件仍优先写入 `blockedDynamicDirections.md`。
- 排障主记录只摘取影响当前优先级的策略结论。

## `blockedDynamicDirections.md` 的维护方式

每个条目建议包含：

- 已尝试方向；
- 当时想解决的问题；
- 失败现象或副作用；
- 当前结论：放弃 / 暂缓 / 仅 dev 可用；
- 未来允许重试的条件；
- 相关 replay / scan / regression / 排障记录。

不要写入：

- 完整运行日志；
- 大段历史复述；
- 与 README 重复的目录说明；
- 与 regression README 重复的契约说明；
- 与 replay README 重复的参数说明。

## 当前不建议回流主路径的典型方向

### 粗 DoA 下直接做 `[fdRef, fdRate]` joint 初始化

- 当前结论：放弃。
- 原因：粗 DoA 误差、差分 Doppler 映射误差和 CP profile 敏感性会耦合，甚至可能改坏单星 dynamic。
- 可重试条件：只作为 dev/probe 对照，且 DoA 初值已锁定。

### alias-prior refine 作为主修复

- 当前结论：放弃。
- 原因：容易把诊断量、主 objective 和主参数语义混在一起。
- 可重试条件：只能旁路诊断，不能灌回主似然。

### CP-U 显式塞 `fdRate=0` 或全局 coarse anchor

- 当前结论：放弃。
- 原因：会把 CP-U 吸到错误但有吸引力的锚点。
- 可重试条件：只允许局部 seed family，不能回到 blanket coarse anchor。

### blanket DoA release / blanket anchorDoaPolish

- 当前结论：放弃。
- 原因：可能误伤 easy / median case。
- 可重试条件：只在 fd healthy + same-tooth + DoA suspicious 的 hard case 条件触发。

### truth override / alias-aware 字段灌回主 objective

- 当前结论：放弃。
- 原因：会污染主 residual、主参数语义和默认分支选择。
- 可重试条件：只允许 probe/replay 旁路，不进入正式 estimator。

### warm-anchor 内层 parfor 作为 estimator 默认路径

- 当前结论：放弃默认化。
- 原因：已观察到可能改变 wrong-tooth sentinel 的 winner。
- 可重试条件：仅 dev/perf 显式 opt-in，并配合 replay / regression 复核。

## 与 replay / scan / regression 的关系

### strategy 到 replay

当策略比较暴露出某个固定代表样本时，迁移到 replay：

- 固定 seed；
- 打印 trace；
- 可画图；
- 可保存 snapshot。

### strategy 到 scan

当需要系统扫参数或 schedule 时，迁移到 scan：

- sweep 维度明确；
- 输出曲线或曲面；
- 不作为日常 quick。

### strategy 到 regression

只有当策略已经变成稳定契约时，才迁移到 regression：

- 能自动 assert；
- 不需要人工判断；
- runtime 可接受；
- 不和已有 regression 重复。
