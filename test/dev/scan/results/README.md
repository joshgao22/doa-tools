# scan results 说明

本目录记录 scan 的详细运行结果、曲线 / 曲面观察、summary table 和 snapshot 绑定。`test/dev/scan/README.md` 只保留 scan 脚本规范、入口索引和当前观察目标；具体结果放在这里。

## 与 cache snapshot 的关系

- 大 `.mat` 文件放在 `test/data/cache/scan/`。
- 本目录不保存 `.mat`，只记录 snapshot 文件名、关键配置、主要 surface / table 结论和当前解释。
- 如果 snapshot 文件未随仓库同步，本目录中的摘要仍作为人工排障和论文图候选记录。
- 若本地有 snapshot，可恢复 `scanData` 后运行对应 scan 的 `Summary output and plotting` section 重出表格和图。

## 当前结果文档

- `scanMfRegimeMapByWindow.md`：多帧窗口下 fixed-DoA / static-Doppler / first-order-Doppler regime 量级分析。
- `scanMfBlockLength.md`：pilot block length 对 `fdRef` comb / alias tooth separation 与块内 Doppler-rate 量级的影响。
- `scanMfCpIpTying.md`：CP/relaxed/IP phase tying 对 `fdRef` comb tooth 结构的影响。
- `scanMfCpIpPerfMap.md`：CP/IP 在 SNR 和 frame count 上的 full-flow 性能 scan；当前 snapshot 作为负结果，说明 full-flow 仍被 wrong-tooth / same-tooth basin 污染，不能直接作为 CP/IP 论文性能图。
- `scanMfCpIpInToothPerfMap.md`：CP/IP 在 truth-centered in-tooth 条件下的受控性能 scan；当前 snapshot 显示 CP 稳定保住 center tooth、`fdRef` consistency 和 non-ref coherence，但 IP 在 angle RMSE 上经常更小，应作为 angle-vs-Doppler-consistency trade-off 解释。
- `scanMfKnownUnknownInformationLoss.md`：known / unknown Doppler-rate 的 CRB / EFIM 信息损失 scan。
- `scanMfSubsetBankCoverage.md`：curated / random subset bank 的选齿覆盖与 runtime-cost scan。

## 单个结果文档结构

每个 scan 默认一个结果文档：

```text
results/<scriptName>.md
```

推荐结构：

1. 对应 scan。
2. 扫描目标。
3. Snapshot index。
4. 当前代表性结果。
5. 关键配置与扫描维度。
6. 主要表格 / 曲线 / 曲面现象。
7. 当前结论。
8. 对 replay / regression / 论文图的影响。
9. 历史 / superseded snapshots。

如果某个 scan 有多个 sweep 或大量 snapshot，可以升级为：

```text
results/<scriptName>/README.md
results/<scriptName>/<caseName>.md
```

## 维护规则

- 新 snapshot 先追加到对应结果文档的 Snapshot index。
- 当前代表性结果只保留一个；旧结果标记为 `superseded`，不要全文复制多份。
- scan 结果用于机制证据和论文图候选；只有稳定成自动契约后才迁移到 regression。
