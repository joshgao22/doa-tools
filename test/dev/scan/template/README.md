# scan template

本目录只放 copy-only scan 模板，不放正式可运行的论文或排障 scan 入口。

新增 scan 时，先复制：

```text
test/dev/scan/template/scanTemplate.m
```

到：

```text
test/dev/scan/scanYourTopic.m
```

然后再改文件名、`scanName`、头部配置、task grid、summary 和 plot。模板只规范 section、checkpoint / snapshot / Telegram 壳、轻量 `scanData` 与打印风格，并示范调用 `test/common/scan/` 的薄工程 helper；它不提供通用 scan 框架，也不承载 scan-specific 结果解释。

使用边界：

- `.m` 文件头注释与 inline comment 使用英文；
- 默认参数直接写在 `Scan configuration`；
- 短 scan 不创建 tmp，也不保留 checkpoint 外壳；
- 重 scan 才启用 `checkpointEnable`，checkpoint 位于 `tmp/<scanName>/<stableRunKey>/`；
- snapshot 只保存轻量 `scanData`，默认变量选择为 `includeVars={'scanData'}`；
- Telegram 指标行由具体 scan 本地构造，再传给 `notifyMfScanStatus`；
- 不把 strategy、schedule、candidate selection、resolved/outlier 分类或 scan-specific metric parser 抽成公共逻辑。
- `printMfScanHeader`、`printMfScanSection`、`notifyMfScanStatus` 和 `finalizeMfScanResult` 只负责工程外壳；具体表格、plot 与 metric 仍留在 scan 本地。

## 两类 scan 的推荐用法

### 轻 scan

适合 CRB / EFIM、regime map、单曲线或少量 grid。保留 `checkpointEnable=false`，不创建 tmp。若需要保存结果，打开 `saveSnapshot=true`。

### 重 scan

适合多 repeat、多 strategy、多 schedule 或耗时 task grid。打开：

```matlab
checkpointEnable = true;
checkpointResume = true;
checkpointCleanupOnSuccess = true;
```

重 scan 的 task 输出应是轻量结构；checkpoint 不保存 `rxSigCell`、完整 `sceneSeq`、fixture cache、transition bundle、全量 objective map、完整 debug trace 或图片。成功构造 `scanData` 后清理 tmp；失败时保留 checkpoint 路径并在命令行和可选 Telegram 通知中提示。

## Telegram 通知格式模板

scan 模板采用“固定外壳 + 弹性指标行”的 HTML 通知格式。固定外壳包含状态、脚本名、耗时、snapshot 文件/目录、checkpoint 目录、失败错误和少量配置。各 scan 只在本地构造 `metricLineList` 与 `commentLineList`。

指标行要求：

- 只放 3–6 行最关键 summary，不替代 results 文档；
- 每行尽量短，适合手机通知阅读；
- `<code>` 只用于短数字、短 tag 或短状态，不包裹完整 snapshot 路径、长 strategy 名或长错误信息；
- 动态文本若不是刻意写入 HTML tag，必须先做 HTML escape；
- `commentLineList` 写本次 scan 的结论、recommendation 或下一步判断，不只写 `completed`。

通知仍是 best-effort 旁路：失败只 `warning`，不得改变 scan 数值路径、`scanData`、checkpoint、snapshot、表格、图或 results 口径。
