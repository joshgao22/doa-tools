# paper 图脚本说明

`test/paper/` 用于放置论文 / 汇报图的轻量重绘脚本。这里的脚本默认**不重新运行 estimator、scan 或 replay**，只读取已经保存的轻量 `scanData` / `replayData` snapshot，并把已有结果整理成更适合论文或 slides 的图。

## 当前入口

### `plotMfRegimeMapByWindowPaper.m`

- 来源数据：`test/data/cache/scan/scanMfRegimeMapByWindow_20260427-213020.mat`。
- 作用：读取 `scanMfRegimeMapByWindow` 的 `scanData`，将 DoA drift、Doppler drift、常 Doppler 相位误差 / 一阶 Doppler residual 合并重画为一张三子图。
- 用途：支撑论文开头的 **DoA 准静态 / Doppler 动态显著** 观测区间说明。
- 默认行为：只画图并打印 `P=10`、`P=20`、`T_f=1/750 s` 的代表性数值；不保存图片。
- 若需要导出图片：将脚本头部 `exportFigure=false` 改为 `true`，图片导出到 `tmp/plotMfRegimeMapByWindowPaper/`，不要提交生成图片。


### `plotSyncBlockObservationWindowPaper.m`

- 来源数据：不读取 replay / scan snapshot，直接使用论文设定中的代表性帧结构参数。
- 作用：画出帧首同步块、多帧观测窗口、连续相位 CP 模型与独立相位 IP 基线的概念示意图。
- 用途：支撑论文模型部分中“帧首已知同步块 + 多帧连续相位观测窗口”的建模说明。
- 默认行为：只画图并打印关键参数；同步块真实持续时间按数值标注，图中宽度为便于阅读的示意性放大。
- 若需要导出图片：将脚本头部 `exportFigure=false` 改为 `true`，图片导出到 `tmp/plotSyncBlockObservationWindowPaper/`，不要提交生成图片。


### `plotMfKnownUnknownInformationLossPaper.m`

- 来源数据：`test/data/cache/scan/scanMfKnownUnknownInformationLoss_20260428-151356.mat`。
- 作用：读取 `scanMfKnownUnknownInformationLoss` 的 `scanData`，重画 known-rate / unknown-rate EFIM 信息损失机制图。
- 图中文字：全部使用英文，适合直接放入论文或 slides。
- 默认主切片：`SNR=10 dB`、`T_f=1/750 s`、`rightBiasedRef` 时间原点；曲线展示 single-satellite 与 multi-satellite。
- 默认输出：三子图，分别展示 reference-Doppler CRB 标准差回退、DoA CRB 标准差回退、EFIM trace loss。
- 默认行为：只画图并打印 `P=10`、`P=20` 的代表性数值；不保存图片。
- 若需要导出图片：将脚本头部 `exportFigure=false` 改为 `true`，图片导出到 `tmp/plotMfKnownUnknownInformationLossPaper/`，不要提交生成图片。


### `plotMfCpIpTradeoffPaper.m`

- 来源数据：`test/data/cache/scan/scanMfCpIpInToothPerfMap_20260428-213853.mat`。
- 作用：读取 `scanMfCpIpInToothPerfMap` 的 `scanData`，重画受控 in-tooth 条件下 CP/IP 折中关系图。
- 图中文字：全部使用英文，适合直接放入论文或 slides。
- 默认主切片：`SNR=10 dB`、`T_f=1/750 s`；曲线展示 `CP-K`、`CP-U`、`IP-K`、`IP-U`。
- 默认输出：三子图，分别展示 DoA angle RMSE、reference-Doppler RMSE 与 non-ref coherence floor。
- 使用口径：该图只表示 truth-centered in-tooth 条件下的模型折中关系，不代表完整 tooth acquisition / full-flow 性能。
- 默认行为：只画图并打印 `P=10`、`P=20` 的代表性数值；不保存图片。
- 若需要导出图片：将脚本头部 `exportFigure=false` 改为 `true`，图片导出到 `tmp/plotMfCpIpTradeoffPaper/`，不要提交生成图片。


### `plotMfMleCrbConsistencyPaper.m`

- 来源数据：`test/data/cache/scan/scanMfMleCrbInToothConsistency_20260507-222606.mat` 与 `test/data/cache/scan/scanMfMsMleCrbInToothConsistency_20260510-062927.mat`。
- 作用：读取 SS / MS 的 in-tooth MLE-vs-CRB scan snapshot，将 trim 后的 local-regime estimator 误差与对应 CRB 曲线重画到同一张 paper-facing consistency 图中。
- 默认主配置：`P=10`、`T_f=1/750 s`；SS 使用 `SS-MF-CP-K/U`，MS 使用 `MS-MF-CP-K/U`。
- 默认输出：四子图，分别展示 SS DoA、SS reference-Doppler、MS DoA、MS reference-Doppler 的 trim RMSE 与 CRB std，纵轴为对数坐标；命令行打印 trim RMSE/CRB、trim MSE/CRB 与 keep-rate 紧凑表。
- 可选口径：脚本头部 `plotMetric` 可设为 `"rmse"`、`"mse"` 或 `"mseOverCrb"`；默认 `"rmse"` 用于直接比较 RMSE 与 CRB 标准差。
- 解读边界：SS snapshot 可作为当前 paper-facing consistency 结果；MS snapshot 主要用于说明 trim 后仍存在的 MS 边界，不应直接包装成 MS 已贴近 CRB 的主结论。
- 默认行为：只画图并打印紧凑表；不保存图片。
- 若需要导出图片：将脚本头部 `exportFigure=false` 改为 `true`，图片导出到 `tmp/plotMfMleCrbConsistencyPaper/`，不要提交生成图片。

## 放置边界

- 本目录只做 paper-facing 重绘、轻量表格整理和图形排版。
- 不在这里维护第二套 estimator、CRB、scan 或 replay 逻辑。
- 若需要新增数值结果，先到 `test/dev/scan/` 或 `test/dev/replay/` 生成轻量 snapshot，再在本目录读取 snapshot 重画。
