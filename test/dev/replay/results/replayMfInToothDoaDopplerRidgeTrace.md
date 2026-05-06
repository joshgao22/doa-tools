# replayMfInToothDoaDopplerRidgeTrace 结果记录

## 0. 状态摘要

| 项目 | 内容 |
|---|---|
| 当前状态 | `diagnostic-only / mechanism visualization` |
| 最新代表性 snapshot | `test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat` |
| 当前一句话结论 | final-centered DoA-Doppler ridge 证明 bad basin 附近存在 coupling，但 minimum 常贴当前 scan box 边界；当前更支持 gated `wide + single-MF` basin-entry，而不是固定 ridge slope 或 joint fd offset 默认化。 |
| 决策影响 | 保留为机制解释和候选设计依据；不改 estimator / flow 默认路径，不新增 regression。 |
| 下一步动作 | 若继续此线，只作为 flow-like gated basin-entry 的机制补充；不要继续把 ridge slope 公式化。 |
| 禁止误用 | 不能把 ridge surface 当作严格 FIM 证明；不能把 truth line / ridge best offset 下沉 runtime；不能默认引入 nonzero `fdRef/fdRate` joint bank。 |

## 1. Replay 身份

- 脚本：`test/dev/replay/replayMfInToothDoaDopplerRidgeTrace.m`
- 结果文档：`test/dev/replay/results/replayMfInToothDoaDopplerRidgeTrace.md`
- replay 类型：机制可视化 replay / fixed-seed ridge trace。
- 主要问题：same-tooth hard seed 的 bad basin 是否存在 DoA-fdRef / DoA-fdRate coupling；final-centered 小扰动是否足够，还是需要更换 basin-entry center。
- 观察范围：固定 seeds `277, 283, 298, 256, 293, 280, 268, 284`；half-tooth oracle in-tooth 条件；ridge surface、center-to-center line、candidate winner 和 rescue bank 对照。
- 不覆盖范围：不做 MC 统计；不验证 full-flow；不改变 estimator 主 objective / residual；不证明 joint fd bank 应进入默认路径。
- truth 使用口径：truth 用于离线 line probe、tail label 和结果评价；ridge / center / rescue 的 runtime 可实现候选必须来自 data-derived center，不允许使用 truth。

## 2. 机制词典与方法地图

| 名称 | 含义 | 是否使用 truth | 改变了什么 | 如何解读 |
|---|---|---:|---|---|
| `DoA-fdRef ridge` | 固定 seed 下沿 DoA offset 与 `fdRef` offset 扫描 objective 的局部 surface。 | Evaluation / diagnostic | 只改诊断扫描点。 | 用来观察 coupling，不等价于可部署的 slope rule。 |
| `DoA-fdRate ridge` | 固定 seed 下沿 DoA offset 与 `fdRate` offset 扫描 objective 的局部 surface。 | Evaluation / diagnostic | 只改诊断扫描点。 | 可解释 bad basin，但不能直接下沉 default optimizer。 |
| `boundary-hit minimum` | ridge best point 贴到当前 scan box 边界。 | Evaluation only | 只改结果标注。 | 表示小盒没包住稳定 minimum；不适合提炼固定扰动公式。 |
| `center-to-center line` | 从 default / static / wide / single-MF center 到其他 center 的线扫。 | Mixed：truth line 是 oracle，data-derived center 是 No | 只改诊断。 | 用来判断需要换 basin-entry center，还是 final 附近 polish。 |
| `wide-centered` | 从 wide DoA candidate center 重新进入 basin。 | No | 改 DoA basin center。 | 能救一部分 hard seed；应 gated 使用。 |
| `single-MF-centered` | 从单星多帧 result center 重新进入 basin。 | No | 改 DoA basin center。 | 能救另一部分 hard seed；blanket 使用会误伤。 |
| `gated-wide-single-bank` | 用 no-truth coherence collapse gate 触发 `wide + single-MF` bank。 | No | 改 rescue 触发与候选 adoption。 | 当前 ridge replay 里最稳候选，但仍只是 fixed-seed 机制证据。 |
| `truth-DoA oracle` | truth-centered DoA 上限。 | Oracle only | 只改上限评价。 | 说明上限存在，不能进入 runtime。 |

## 3. Snapshot index

| snapshot | 日期 | 状态 | 配置摘要 | 结论 | 覆盖 / 取代 |
|---|---:|---|---|---|---|
| `test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat` | 2026-04-29 | representative | `snrDb=10`，`contextBaseSeed=253`，`seedList=[277 283 298 256 293 280 268 284]`，ridge DoA offsets、`fdRef` offsets、`fdRate` offsets、center-line sweep、gated rescue bank。 | ridge 证明 coupling 存在但不稳定可公式化；`gated-wide-single-bank` 仍是更稳的可实现方向。 | 当前唯一代表性结果。 |

## 4. 最新代表性运行

### 4.1 配置

- `snrDb = 10`
- `contextBaseSeed = 253`
- `numRepeat = 8`
- seed list：`277, 283, 298, 256, 293, 280, 268, 284`
- hard-collapse seeds：`283, 298, 256`
- same-tooth light / unclear seeds：`277, 293, 280`
- fd-not-healthy negative seeds：`268, 284`
- `fdRef` oracle 范围：truth-centered `0.49` tooth
- `fdRate` oracle 范围：truth-centered `±1000 Hz/s`
- ridge scan：DoA offset、`fdRef` offset、`fdRate` offset 与 center-line alpha sweep。
- gated rescue thresholds：`defaultNonRefCoherenceFloor < 0.20` 且 `defaultNonRefRmsPhaseResidRad >= 1.0`
- snapshot 保存变量：`replayData`

### 4.2 主要统计

| 指标 | 数值 | 解释 |
|---|---:|---|
| `MS-MF-CP-U-in-tooth toothHitRate` | 1 | 当前 replay 已剥离 wrong-tooth。 |
| `MS-MF-CP-U-in-tooth angle RMSE` | `0.0031024 deg` | same-tooth DoA/local-state tail 仍明显。 |
| `SS-MF-CP-U-in-tooth angle RMSE` | `0.0014248 deg` | 多星默认路径在 hard seeds 上仍可差于单星。 |
| `truth-DoA oracle angle RMSE` | `0.00036916 deg` | MS-MF 上限存在，主要问题是 basin entry。 |
| `gated-wide-single-bank trigger rate` | 0.375 | 只触发 ridge snapshot 中的 hard-collapse seeds。 |
| `gated-wide-single-bank hard trigger rate` | 1 | hard seeds 全部被 no-truth gate 抓到。 |
| `gated-wide-single-bank easy / fd-negative trigger rate` | 0 / 0 | 负样本不触发。 |
| `gated-wide-single-bank hard coherence recovered rate` | 1 | hard seeds 的 non-ref coherence 被恢复。 |

### 4.3 关键对比表

#### Ridge / center observation

| 现象 | 观察 | 解读 |
|---|---|---|
| final-centered ridge | hard seeds 的 best point 经常贴近 DoA `+0.0045~+0.006 deg`、`fdRef≈-0.08 Hz`、`fdRate≈-160 Hz/s` 等 scan box 边界。 | coupling 存在，但小盒没有稳定包住 minimum，不能公式化 slope。 |
| center-to-center line | hard seeds 的好点更接近 wide 或 single-MF center，而不是 default final 附近。 | 需要换 basin-entry center，而不是 very-small polish。 |
| candidate winner | `283` 由 wide-centered 救，`298/256` 由 single-MF-centered 救。 | 候选 bank 不能只保留一种 center。 |
| rescue gate | gated bank 只触发 hard seeds，easy / fd-negative 不触发。 | no-truth collapse gate 是比 blanket bank 更安全的方向。 |

#### Implementable center examples

| seed | default angle (deg) | best implementable center | selected angle (deg) | 现象 |
|---:|---:|---|---:|---|
| 283 | 0.003618 | `wide-centered` | 0.00052846 | coherence 恢复到 1。 |
| 298 | 0.0053878 | `single-MF-centered` | 0.00047757 | coherence 恢复到 0.99975。 |
| 256 | 0.0054541 | `single-MF-centered` | 0.00034933 | coherence 恢复到 0.99845。 |

#### Rescue bank aggregate

| bank | hard rescue rate | easy damage | fd-negative damage | 结论 |
|---|---:|---:|---:|---|
| `disabled` | 0 | 0 | 0 | baseline。 |
| `wide-centered-coarse` | 1/3 | 0 | 0 | 只能救部分 hard。 |
| `single-mf-centered-coarse` | 2/3 | 1 | 1 | 会误伤 easy / fd-negative。 |
| `wide-single-bank` | 1 | 0 | 0 | hard 全救，但 blanket 常驻仍不安全。 |
| `gated-wide-single-bank` | 1 | 0 | 0 | 只触发 hard-collapse，当前最稳。 |

## 5. 可观察现象

### 5.1 支持当前结论的现象

- final-centered surface 的更优方向常同时包含 DoA 与 Doppler offset，说明 bad basin 附近确有 DoA-Doppler coupling。
- ridge minimum 贴边界，说明当前小扰动没有形成稳定可复用的局部公式。
- `wide` 和 `single-MF` center 能把 hard seeds 拉回好 basin，说明更换 basin-entry center 比 final-centered very-small polish 更有效。
- `gated-wide-single-bank` 在本 snapshot 中 hard rescue rate 为 1、easy / fd-negative damage 为 0，和 tail diagnose 的方向一致。

### 5.2 仍未解决或反向的现象

- ridge slope / best offset 不稳定，不能提炼成固定 `fdRef` 或 `fdRate` offset rule。
- 本 replay 是 8 seed 机制可视化，不提供大 MC damage 率。
- 若把 `single-MF` 或 `wide` center blanket 常驻，仍存在误伤风险；gate 才是重点。

### 5.3 代表性 seed / case

| seed | 类型 | 现象 | 对结论的作用 |
|---:|---|---|---|
| 283 | hard-collapse | wide-centered 可救到 `0.00052846 deg`。 | 说明 wide center 不能删。 |
| 298 | hard-collapse | single-MF-centered 可救到 `0.00047757 deg`。 | 说明 single-MF center 不能删。 |
| 256 | hard-collapse | single-MF-centered 可救到 `0.00034933 deg`。 | 说明 hard seed 更需要 basin-entry center。 |
| 268 / 284 | fd-negative | 不应被 DoA ridge / basin-entry 误伤。 | 说明必须保留 gate / damage 指标。 |

## 6. 机制解释

### 6.1 当前解释

ridge surface 说明同齿 bad basin 不是纯 DoA 或纯 Doppler 单变量问题：在 default final 附近，objective 的下降方向可能同时包含 DoA 与 `fdRef/fdRate` 的小偏移。但是这些局部 surface 的 minimum 经常贴在当前 scan 边界上，说明它们更多是“解释 bad basin 形状”的诊断，而不是一个可直接写入 flow 的稳定修正公式。

更稳的可实现现象来自 center 级别：把起点换到 wide 或 single-MF data-derived center 后，hard seed 的 non-ref coherence 可以恢复。这说明当前阶段应优先推进 gated basin-entry family，而不是把 `ridgeSlopeHzPerDeg`、truth line 或 nonzero fd offset 下沉默认路径。

### 6.2 这个结果支持什么

- 支持“same-tooth bad basin 存在 DoA-Doppler coupling”的机制解释。
- 支持“hard seed 需要 basin-entry center，而不是 final-centered very-small polish”。
- 支持继续把 `wide + single-MF` family 作为 flow-like gated rescue 候选。

### 6.3 这个结果不证明什么

- 不证明可以把 ridge slope 公式化。
- 不证明 nonzero `fdRef/fdRate` joint bank 应默认打开。
- 不证明 estimator objective / residual 需要修改。
- 不证明 full-flow 性能已经可用于论文主图。

## 7. 对主流程的影响

| 项目 | 影响 |
|---|---|
| estimator 默认路径 | 不改；不下沉 ridge / joint fd 规则。 |
| flow 默认路径 | 暂不改；只支持后续 flow-like gated basin-entry 验证。 |
| regression | 不写；surface / ridge 是机制诊断，不是 pass/fail contract。 |
| replay / scan 下一步 | 与 `TailCaseDiagnose` 一起作为 flow-like gated rescue 的机制依据；不继续扩大 ridge slope 搜索。 |
| 论文图 / 论文口径 | 可作为附录或 diagnostic 机制图候选；不作为主性能图。 |
| 排障记录 | 机制归并版可保留“ridge 解释 coupling，但不支持公式化 slope”。 |

## 8. 限制与禁止解释

- 不要把 ridge minimum 当成稳定 estimator 更新方向。
- 不要把 truth line probe、truth DoA 或 tail label 放入 runtime gate / adoption。
- 不要因为 surface 有 coupling 就默认启用 nonzero `fdRef/fdRate` joint bank。
- 不要把 8 seed 机制结果当作大 MC damage 率。
- 不要让本 replay 覆盖 `TailCaseDiagnose` 的 fixed-tail rescue 预验证职责；两者互补。

## 9. 恢复与复现

```matlab
snapshotFile = 'test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat';
loadExpSnapshot(snapshotFile, 'caller', struct('varNames', {{'replayData'}}));
```

随后打开：

```text
`test/dev/replay/replayMfInToothDoaDopplerRidgeTrace.m`
```

只运行 `Summary output and plotting` 小节，即可重出 compact table 和图。

## 10. 历史备注

- 本 replay 承接 `FdRangeOracle` / `TailCaseDiagnose` 的 same-tooth tail seeds。
- 当前结论不是“ridge 线可下沉”，而是“ridge 解释 coupling，真正可实现方向仍是 gated basin-entry center”。
