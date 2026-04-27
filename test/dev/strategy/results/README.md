# strategy results 说明

本目录可选用于记录策略比较的详细结果。策略结论如果只是当前路线选择，优先写入 `blockedDynamicDirections.md` 或排障主记录；只有当对比表较长、需要绑定 snapshot 或需要保留多组策略实验时，才在本目录新增结果文档。

## 推荐结构

```text
results/<strategyName>.md
```

每个文档建议记录：

- 对应 strategy 脚本；
- 比较目标；
- snapshot / cache 文件名；
- 策略配置；
- 主要对比表；
- 当前采用 / 放弃 / 暂缓的结论；
- 允许重试的条件。

## 与其他文档的关系

- `blockedDynamicDirections.md` 记录已证伪或不建议重复投入的方向。
- `results/*.md` 记录长表格和具体证据。
- 排障主记录只摘取影响当前优先级的策略结论。
