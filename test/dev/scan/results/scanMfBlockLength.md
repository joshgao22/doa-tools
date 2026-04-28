# scanMfBlockLength 结果记录

## 对应 scan

- `test/dev/scan/scanMfBlockLength.m`

## 扫描目标

验证 pilot block length 对连续相位多帧 MF `fdRef` comb / alias tooth separation 的影响，并同时记录单个同步块内部 Doppler-rate 动态量级。

该 scan 使用同一组场景、同一个 shared phase 和同一份最大长度 pilot / snapshot，裁剪出不同 block length 后比较 `CP-K` 与 `CP-U` 的 `fdRef` 一维 objective line。它是 block-length mechanism scan，不用于比较 subset bank、rescue flow 或 same-tooth polish，也不把结果迁移为 regression 契约。

## Snapshot index

| snapshot | 配置 | 结论 |
|---|---|---|
| `test/data/cache/scan/scanMfBlockLength_20260428-105157.mat` | `Fs=512 MHz`，`symbolRate=128 MHz`，`osf=4`，`baseBlockLen=2112`，`standardBlockLen=8448`，`numFrame=10`，`T_f=1/750 s`，`SNR=10 dB`，`blockLen=[528,1056,2112,4224,8448]`，`numGrid=801`，`numAliasSide=2`，`CP-K/CP-U` 均运行，checkpoint resume 开启 | block length 增长会显著增大 truth-centered alias tooth gap；`8448` samples 下块内二次相位仍只有 `3.28e-6 rad`，因此本结果支持“长 pilot block 增强单块 Doppler / tooth 分辨力”，不支持“单块内部 Doppler-rate 已显著”的解释。 |

## 当前代表性结果

2026-04-28 的代表性结果扫描了 5 个 block length，每个 block 分别评估 `CP-K` 与 `CP-U`，并在 `truth / staticSeed / finalEstimate` 或 `truth / cpKnownSeed / finalEstimate` 三类中心附近扫描 `fdRef` line。checkpoint 从 `3/5` 已完成任务恢复，最终 `numTask=5`、`numDone=5`、`isComplete=true`，成功构造 `scanData` 后清理 checkpoint run 目录。

### Truth-centered alias gap

下表只取 `truth` 中心的 `fdRef` line。由于 `CP-K` 与 `CP-U` 在 truth-centered line 上得到相同 gap，表中合并展示两种 rate mode 的共同结果。

| blockLen | 相对 2112 | block duration (us) | abs(fdRate) * T_block (Hz) | quad phase (rad) | aliasGap1 | aliasGap2 |
|---:|---:|---:|---:|---:|---:|---:|
| 528 | 0.25x | 1.03125 | 0.003953 | 1.28e-8 | 0.2998 | 1.2366 |
| 1056 | 0.5x | 2.06250 | 0.007907 | 5.12e-8 | 2.7226 | 10.9158 |
| 2112 | 1x | 4.12500 | 0.015813 | 2.05e-7 | 24.3386 | 98.7082 |
| 4224 | 2x | 8.25000 | 0.031626 | 8.20e-7 | 275.3637 | 1118.0680 |
| 8448 | 4x | 16.50000 | 0.063253 | 3.28e-6 | 13574.8601 | 53875.1723 |

该表是本 scan 的主结果。`528` samples 时相邻 tooth 的 objective gap 只有 `0.3` 量级，comb 很平；`2112` samples 时相邻 tooth gap 已增至 `24.34`；`4224` 与 `8448` samples 进一步把 wrong tooth 明显压低。说明 block length 对 `fdRef` tooth separation 的影响非常强。

### Known / unknown aggregate summary

| mode | blockLen | aliasGap1 | aliasGap2 | centerDeltaFinal | minDeltaFinalHz |
|---|---:|---:|---:|---:|---:|
| CP-K | 528 | 0.2998 | 1.2366 | 3.1360 | 1500 |
| CP-K | 1056 | 2.7226 | 10.9158 | 167235.0001 | -375 |
| CP-K | 2112 | 24.3386 | 98.7082 | 0 | 0 |
| CP-K | 4224 | 275.3637 | 1118.0680 | 0 | 0 |
| CP-K | 8448 | 13574.8601 | 53875.1723 | 0 | 0 |
| CP-U | 528 | 0.2998 | 1.2366 | 3.1360 | 1500 |
| CP-U | 1056 | 2.7226 | 10.9158 | 0 | 0 |
| CP-U | 2112 | 24.3386 | 98.7082 | 0.6718 | 750 |
| CP-U | 4224 | 275.3637 | 1118.0680 | 8.5554 | 750 |
| CP-U | 8448 | 13574.8601 | 53875.1723 | 0 | 0 |

`centerDeltaFinal / minDeltaFinalHz` 是 final-centered 诊断，不是本 scan 的主机制指标。短块下它会受到 estimator final basin 的影响，例如 `528` 可能落到 `+1500 Hz`，`1056` 的 `CP-K` final-centered line 还出现很大的 `centerDeltaFinal`。这些现象说明短 block 下 final estimate 更容易不稳，但不应和 truth-centered alias-gap 机制结论混为一谈。

### Block 内 Doppler-rate 量级

本次 truth `fdRateFit` 约为 `-3833.5 Hz/s`。即使在最长的 `8448` samples 下，单个 block 也只有 `16.5 us`，块内 Doppler drift 约 `0.063 Hz`，二次相位约 `3.28e-6 rad`。因此本 scan 不能用来证明“单个同步块内部必须建模 Doppler rate”；它证明的是 block length 增加带来了更强的单块相干积分与 Doppler / tooth 分辨能力。

## 可观察现象

- `aliasGap1` 和 `aliasGap2` 随 block length 单调、快速增大；从 `2112` 到 `8448`，`aliasGap1` 从 `24.34` 增至 `1.36e4`。
- `CP-K` 与 `CP-U` 的 truth-centered alias gap 在当前设置下完全一致，说明该 scan 中 tooth separation 的主导因素是 block length，而不是 known-rate / unknown-rate 分支差异。
- `truthGap=0` 且 `minDeltaTruthHz=0` 对全部 block 成立，说明在 truth center 下，当前 evaluator 的最小点仍在中心 tooth。
- final-centered 诊断在短 block 下不稳定，尤其 `528` 与 `1056`；这类现象更适合解释 estimator basin，而不是解释 block-length 物理机制。
- checkpoint 机制已能按 block task 恢复，成功后清理 tmp 目录，符合重 scan 的工程存储要求。

## 当前结论

当前结果足够支撑 `scanMfBlockLength` 的预期结论：**更长 pilot block 会显著增强 `fdRef` comb / alias tooth separation**。默认 `2112` samples 已经比 `528 / 1056` 明显更稳；`2x / 4x` 更长 block 进一步显著提高 tooth gap。

同时，结果也明确了边界：**4x block 的块内 Doppler-rate 动态仍可忽略**。因此，论文主线中 Doppler dynamic 的必要性仍应来自跨帧连续相位窗口，而不是单个同步块内部的二次相位。若正文使用该 scan，应写成“longer known block improves single-block Doppler/tooth discriminability”，不要写成“longer block makes in-block Doppler rate dominant”。

## 对 replay / regression / 论文图的影响

- 该 scan 可作为 block-length mechanism evidence 和论文图候选：主图建议画 `aliasGap1 / aliasGap2` 随 `blockLenBaseRatio` 的 log-y 曲线。
- 副图可画 `inBlockQuadPhaseRad` 或 `inBlockFdDriftHz`，用于说明单块内部 rate 动态仍很小。
- `finalEstimate` 中心的 `minDeltaFinalHz` 只作为诊断图，不应和 truth-centered alias-gap 主图混放。
- 该结果不直接进入 regression。它解释的是机制趋势，不是自动 pass/fail 契约。
- 若后续要验证“更长 block 是否提高真实 flow 选齿命中率”，应另做 flow / MC scan 或扩展 subset-flow 入口；不要继续把 subset bank、rescue、same-tooth polish 塞进 `scanMfBlockLength`。

## 后续建议

`scanMfBlockLength` 当前不需要继续增加 Monte Carlo、subset strategy 或 rescue 维度。建议只保留当前结果，并在后续代码中优化图形口径：主图突出 truth-centered alias gap，诊断图单独展示 final-centered min-tooth 偏移。

## 历史 / superseded snapshots

- 无。
