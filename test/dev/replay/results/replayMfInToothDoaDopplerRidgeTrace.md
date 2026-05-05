# replayMfInToothDoaDopplerRidgeTrace 结果记录

## 观察目标

本 replay 用固定 in-tooth tail / negative seeds 观察 DoA-Doppler 局部 ridge 和可实现 basin-entry center。重点不是给出性能 MC，而是判断 final-centered 小扰动是否值得下沉为规则，或者是否应继续优先使用 gated `wide + single-MF` basin-entry bank。

## Snapshot index

| snapshot | seeds | SNR | 说明 |
|---|---:|---:|---|
| `test/data/cache/replay/replayMfInToothDoaDopplerRidgeTrace_20260429-100106.mat` | 8 | 10 dB | 当前代表性结果；包含 hard-collapse、light/easy 与 fd-not-healthy negative seeds。 |

## 当前代表性结果

固定 seeds：`277 / 283 / 298 / 256 / 293 / 280 / 268 / 284`。

- hard-collapse：`283 / 298 / 256`。
- same-tooth light / unclear：`277 / 293 / 280`。
- fd-not-healthy negative：`268 / 284`。

`MS-MF-CP-U-in-tooth` 在 half-tooth oracle 下仍有 same-tooth DoA/local-state tail：tooth hit 为 1，但 angle RMSE 为 `0.0031024 deg`，差于 `SS-MF-CP-U-in-tooth` 的 `0.0014248 deg`。`MS-MF-CP-U-truth-doa-oracle` 可到 `0.00036916 deg`，说明 MS-MF 的上限存在，主要问题是 bad basin entry。

## 主要现象

### 1. final-centered ridge 不能直接公式化

hard seeds 的 DoA-fdRef / DoA-fdRate surface 经常在当前 scan box 边界取得 best candidate，例如 DoA offset 贴近 `+0.0045 ~ +0.006 deg`，fdRef offset 贴近 `-0.08 Hz`，fdRate offset 贴近 `-160 Hz/s`。这说明坏点附近确实存在联合方向，但当前小盒并没有稳定包住 minimum。

因此，本 replay 不支持把 `ridgeSlopeHzPerDeg` 或固定 fdRef/fdRate 扰动斜率直接写入 flow。

### 2. 可实现 center 比 final-centered 小扰动更关键

center-to-center line 与 candidate winner 结果显示，hard seeds 的主要救回来自 truth-free basin-entry center：

| seed | default angle (deg) | best implementable center | selected angle (deg) | 现象 |
|---:|---:|---|---:|---|
| 283 | 0.003618 | wide-centered | 0.00052846 | coherence 恢复到 1 |
| 298 | 0.0053878 | single-MF-centered | 0.00047757 | coherence 恢复到 0.99975 |
| 256 | 0.0054541 | single-MF-centered | 0.00034933 | coherence 恢复到 0.99845 |

这说明 hard seed 更像需要换 basin-entry center，而不是在 final 附近做 very-small polish。

### 3. gated-wide-single-bank 是当前最稳候选

rescue bank aggregate 中：

| bank | hard rescue rate | easy damage | fd-negative damage | 结论 |
|---|---:|---:|---:|---|
| disabled | 0 | 0 | 0 | baseline |
| wide-centered-coarse | 1/3 | 0 | 0 | 只能救部分 hard |
| single-mf-centered-coarse | 2/3 | 1 | 1 | 会误伤 easy / fd-negative |
| wide-single-bank | 1 | 0 | 0 | blanket 下有效但不应默认常驻 |
| gated-wide-single-bank | 1 | 0 | 0 | 只触发 hard-collapse，当前最稳 |

`gated-wide-single-bank` 的 trigger rate 为 `0.375`，hard trigger rate 为 `1`，easy / fd-negative trigger rate 为 `0`，hard coherence recovered rate 为 `1`。

## 当前结论

1. in-tooth bad basin 的确存在 DoA-Doppler coupling。
2. 但 ridge minimum 经常贴边界，暂不应提炼成固定扰动公式。
3. 当前最稳的可实现修复方向是 no-truth collapse-gated `wide + single-MF` basin-entry bank。
4. 下一步应进入 flow-like replay，验证该 bank 在真实 subset-periodic flow 中是否仍能救 hard seed 且不误伤 easy / fd-negative seed。

## 对主流程的影响

- 不改 estimator 主核、objective 或 residual。
- 不新增 regression。
- 不把 truth DoA、truth tooth、truth fdRef/fdRate 或 tail label 放入 gate / adoption。
- 先新增 / 运行 flow-like gated basin-entry replay；通过后再考虑在 `runSimpleDynamicSubsetPeriodicFlow` 的 final same-tooth rescue 层显式启用。
