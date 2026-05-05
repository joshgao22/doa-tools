# replay results 说明

本目录记录 replay 的详细运行结果、观察现象和 snapshot 绑定。`test/dev/replay/README.md` 只保留 replay 入口、运行规范和当前一句话结论；具体统计表、长现象描述和历史 snapshot 索引放在这里。

## 与 cache snapshot 的关系

- 大 `.mat` 文件放在 `test/data/cache/replay/`。
- 本目录不保存 `.mat`，只记录 snapshot 文件名、关键配置、主要统计、观察现象和结论。
- 如果 snapshot 文件未随仓库同步，本目录中的摘要仍作为人工排障记录。
- 若本地有 snapshot，可恢复 `replayData` 后运行对应 replay 的 `Summary output and plotting` section 重出表格和图。

## 单个结果文档结构

每个 replay 默认一个结果文档：

```text
results/<scriptName>.md
```

推荐结构：

1. 对应 replay。
2. 观察目标。
3. Snapshot index。
4. 当前代表性结果。
5. 主要统计。
6. 可观察现象。
7. 当前结论。
8. 对主流程的影响。
9. 历史 / superseded snapshots。

如果某个 replay 的结果数量太多，可以升级为：

```text
results/<scriptName>/README.md
results/<scriptName>/<caseName>.md
```

## 当前代表性结果索引

| replay | 结果文档 | representative snapshot | 当前结论 |
|---|---|---|---|
| `replayMfCombToothSurface.m` | `replayMfCombToothSurface.md` | `test/data/cache/replay/replayMfCombToothSurface_20260425-180232.mat` | `1/T_f` comb 是 objective 层真实结构。 |
| `replayMfPeriodicVsSubsetToothSelect.m` | `replayMfPeriodicVsSubsetToothSelect.md` | `test/data/cache/replay/replayMfPeriodicVsSubsetToothSelect_20260425-182612.mat` | subset 负责选齿，periodic 主要负责同齿 refine；warning 只作为诊断元数据保留。 |
| `replayMfRandomRescueEffectiveness.m` | `replayMfRandomRescueEffectiveness.md` | `test/data/cache/replay/replayMfRandomRescueEffectiveness_20260427-110427.mat` | rescue/random bank 明显提升 central tooth 命中率，但已有 easy-case damage，后续应转 scan/gate 验证而非 blanket 常驻。 |
| `replayMfInToothFdRangeOracle.m` | `replayMfInToothFdRangeOracle.md` | `test/data/cache/replay/replayMfInToothFdRangeOracle_20260426-104249.mat` | half-tooth oracle 下频率链健康，但少数 same-tooth non-ref coherence tail 拉坏 RMSE。 |
| `replayMfInToothTailCaseDiagnose.m` | `replayMfInToothTailCaseDiagnose.md` | `test/data/cache/replay/replayMfInToothTailCaseDiagnose_20260427-102008.mat` | gated `wide+single-MF` rescue bank 能全救 fixed hard-collapse seeds，并避开 easy / fd-not-healthy 负样本。 |
| `replayMfInToothDoaDopplerRidgeTrace.m` | `replayMfInToothDoaDopplerRidgeTrace.md` | `test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat` | final-centered ridge 证明 DoA-Doppler coupling 存在，但 minimum 常贴边界；优先推进 gated `wide+single-MF` basin-entry，而不是公式化 ridge slope。 |
| `replayMfInToothGatedRescueEffectiveness.m` | `replayMfInToothGatedRescueEffectiveness.md` | `test/data/cache/replay/replayMfInToothGatedRescueEffectiveness_20260504-215932.mat` | family-safe-adopt 显著压低 controlled in-tooth MS-MF DoA tail，hard rescue rate 提升到 `0.8`，easy / fd-negative damage 为 0；仍是 replay-level candidate。 |
| `replayMfInToothDoaFdRangeEnvelope.m` | `replayMfInToothDoaFdRangeEnvelope.md` | `test/data/cache/replay/replayMfInToothDoaFdRangeEnvelope_20260504-204938.mat` | base-seed envelope 支持 wide / single-MF DoA basin-entry 与 `0.006 deg` family-safe step，不支持默认引入非零 fdRef / fdRate joint bank。 |

## 维护规则

- 新 snapshot 先追加到对应结果文档的 Snapshot index。
- 当前代表性结果只保留一个；若新结果完全覆盖旧结果且结论一致，可删除旧 snapshot 记录，不必保留 `superseded` 行。
- 结果文档记录证据，排障主记录只摘取影响当前优先级的结论。
