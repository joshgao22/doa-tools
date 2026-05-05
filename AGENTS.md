# AGENTS.md

本文件是 AI / 自动化修改代码前的最短入口。它只说明读文档顺序、规范职责、冲突裁决和默认修改原则，不替代各目录 README，也不记录实验结果。

## 修改前读取顺序

1. 读本文件。
2. 读根目录 `README.md`，确认全局代码风格、注释语言、文件存储和目录职责。
3. 读目标文件所在目录 README。
4. 若目标文件位于更深层目录，继续读上一级和当前目录 README。
5. 若修改 replay / scan / strategy / regression，读对应局部 README。
6. 若修改 estimator helper，读 `estimator/README.md` 与 `estimator/helper/README.md`。
7. 若修改 scene、reference-sat、Doppler 几何或裁剪，读 `satellite/scene/README.md`。
8. 若修改 dynamic flow、subset tooth selection、conditional random rescue、same-tooth polish、warm-anchor 或相关 replay，读当前排障主记录；只有需要机制背景或旧路线回溯时，再读机制归并版和历史归档。

## 规范职责

- `AGENTS.md`：最短读法、规范优先级、冲突裁决和默认修改原则。
- 根 `README.md`：仓库级目录职责、代码风格、落盘规范索引和文档同步规则。
- 目标目录 README：该目录入口、文件职责、局部运行方式和局部代码约束。
- `results/` 文档：replay / scan / strategy 的具体运行结果、长表格、snapshot 绑定和观察现象。
- 排障主记录：当前主问题、当前优先级、必跑护栏和已证伪路线。
- 机制归并版 / 历史归档：机制证据和历史追溯，只在需要理解 why 或追旧线索时读取。
- 外部 coding prompt：约束 AI 的执行风格、风险偏好和交付说明，不作为仓库文件索引的替代品。

## 优先级与冲突裁决

1. 当前用户明确要求。
2. 本文件。
3. 根目录 `README.md`。
4. 目标目录 README。
5. 更内层 README。
6. 当前排障主记录。
7. 机制归并版 / 历史归档。

如果用户要求、外部 prompt 或 README 之间存在冲突，优先选择更小、更局部、更不改变默认数值路径的方案。若冲突涉及函数签名、默认参数、reference-sat 语义、搜索边界、sat 顺序、输出字段、summary 口径或 fallback 顺序，不能把风险藏在实现里，必须在交付说明中明确说明风险和取舍。

## 默认代码策略

- 默认保持行为，不改变数值路径。
- 不机械附和临时方案；如果某个改法会增加接口膨胀、重复文档、上下文负担、已证伪路线回流、默认数值风险或维护成本，应选择更小、更稳的替代实现。
- 不擅自改函数签名、默认参数、初始化策略、搜索边界、reference-sat 语义、sat 顺序、输出字段、summary 口径或 fallback 顺序。
- `.m` 文件中的函数头注释、脚本头注释和 inline comment 使用英文。
- README、排障记录和说明文档使用中文。
- 参考目标目录同类型、非 legacy 文件的风格。
- 不为小改动新增过重 opt、arguments、validateattributes、配置 resolver 或框架化封装。
- 能在上层 orchestration 内部消化的改动，不继续下沉扩展底层接口。
- 正式 estimator 算法逻辑放 `estimator/` 或 `estimator/helper/`；scene / reference-sat / Doppler 几何放 `satellite/scene/`；CRB / FIM / EFIM 放 `performance/`；实验组织和 summary 放 `test/`。

## test 子系统边界

- `test/regression/`：自动 pass/fail 契约护栏。
- `test/dev/replay/`：固定 seed 或小 Monte Carlo 回放，用于记录现象和机制证据。
- `test/dev/scan/`：较重的参数、schedule、曲线或曲面扫描。
- `test/dev/strategy/`：策略比较、已证伪路线和重试准入条件。
- `test/common/`：实验侧复用 helper，不长期维护第二套正式算法。


## 长任务通知边界

- 长时间运行的 replay / scan / Monte Carlo 脚本可在顶层 orchestration 层调用 `utils/io/notifyTelegram.m` 发送完成或失败通知。
- 通知逻辑必须是 best-effort：发送失败只 `warning`，不得影响仿真主结果、snapshot 保存、`replayData` / `scanData` 内容或 estimator 数值路径。
- `notifyTelegram` 只作为通用 I/O 工具使用，不下沉到 `estimator/`、`performance/`、objective、residual、CRB 或正式算法 helper。
- token、chat_id 与本机配置不得写入仓库；优先使用环境变量或 MATLAB `prefdir` 下的私有配置。
- 新增或修改长时间 replay / scan 时，可在脚本结束、失败 `catch`、snapshot 保存后追加简短通知；消息只报告脚本名、状态、耗时、snapshot 路径和少量关键 summary，不维护第二套 results 解析逻辑。

## 结果、落盘与文档边界

- 运行时临时文件统一放仓库根目录 `tmp/<scriptName>/<runKey>/`。
- 大 `.mat` snapshot 统一放 `test/data/cache/<taskType>/`。
- replay 详细运行结果写入 `test/dev/replay/results/`。
- scan 详细运行结果写入 `test/dev/scan/results/`。
- strategy 长结果写入 `test/dev/strategy/results/`。
- perf 详细运行结果写入对应 perf results 文档。
- 保存、恢复和清理入口以 `test/data/cache/README.md` 为单一索引；不要在多个 README 中复制 save/load/delete 细节。
- README 只保留入口、职责、运行方式和输出解释；不记录长结果。
- 排障记录只保留当前优先级、机制结论和已证伪路线，不复制 README 或完整运行日志。
