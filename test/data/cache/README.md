# test/data/cache 说明

`test/data/cache/` 只放可恢复的大数据 snapshot，不放人工分析文档。人工可读结果记录放在对应任务目录的 `results/` 下，例如 `test/dev/replay/results/` 或 `test/dev/scan/results/`。

本 README 是实验落盘、恢复和清理入口的单一索引。`utils/io/` 只放通用实现文件，不单独维护一份重复的 save / load / delete 规范；其它 README 只引用本文件，不复制调用细节。

## 推荐目录

```text
test/data/cache/
├── replay/
├── scan/
├── perf/
└── regression/
```

默认规则：

- replay snapshot 放 `test/data/cache/replay/`；
- scan snapshot 放 `test/data/cache/scan/`；
- perf / 长 Monte Carlo snapshot 放 `test/data/cache/perf/`；
- regression 一般不长期保存大结果，确需保存时放 `test/data/cache/regression/`。

## 文件命名

默认文件名保持：

```text
<scriptName>_yyyymmdd-HHMMSS.mat
```

即使移动到 task type 子目录，文件名也保留脚本名前缀，便于从文件名直接判断来源。

## 落盘、恢复与删除入口

| 操作 | 推荐入口 | 说明 |
|---|---|---|
| 保存最终 snapshot | `utils/io/saveExpSnapshot.m` | replay 保存 `replayData`，scan 保存 `scanData`，perf 保存 final result struct、summary table、关键配置和 meta |
| 恢复 snapshot | `utils/io/loadExpSnapshot.m` | 只恢复需要的变量，例如 `replayData` 或 `scanData`；恢复后运行脚本的 `Summary output and plotting` section |
| 清理 cache / tmp 运行产物 | `utils/io/cleanupRunArtifacts.m` | 通用安全清理入口，删除范围默认限制在仓库 `tmp/` 与 `test/data/cache/` 下 |
| 清理 checkpoint run 目录 | `test/common/flow/cleanupPerfTaskGridCheckpoint.m` | checkpoint runner 专用清理入口，成功组装结果并保存最终 snapshot 后调用 |

常用恢复示例：

```matlab
% Replay snapshot restore
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));

% Scan snapshot restore
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'scanData'}}));
```

使用原则：

- replay / scan 脚本的最终落盘统一放在 `Data storage` section。
- checkpoint 只保存中间 task 结果，不替代最终 snapshot。
- 成功构造 `replayData` / `scanData` 并保存最终 snapshot 后，再清理 tmp / checkpoint。
- 删除 `tmp/` 或 `test/data/cache/` 下运行产物时，优先使用公共清理入口；不要在 replay / scan 脚本里复制私有递归删除 helper。
- 确认只需要读取 metadata 时，使用 `loadExpSnapshot(..., 'none', struct('metaOnly', true))`，避免把大变量恢复到工作区。

## 何时按脚本再建子目录

默认不要按脚本建更深目录。若某个脚本积累了大量 snapshot，才升级为：

```text
test/data/cache/replay/<scriptName>/<scriptName>_yyyymmdd-HHMMSS.mat
```

对应的人工索引仍写在 `test/dev/replay/results/<scriptName>.md` 或该 replay 的结果子目录中。

## 保存内容

snapshot 默认只保存轻量数据结构：

- replay 保存 `replayData`；
- scan 保存 `scanData`；
- perf 保存 final result struct、summary table、关键配置和 meta。

不要默认保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片文件。
