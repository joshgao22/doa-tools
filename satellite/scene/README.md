# satellite/scene 目录说明

`satellite/scene/` 放场景构造、reference-sat 解析、用户状态、Doppler 几何和 scene 裁剪。这里是 static 与 dynamic 共用的几何语义层，不能在 estimator 或 test 中长期维护另一套等价实现。

## 核心职责

### 场景构造

- 多星单帧 scene；
- 多星多帧 sceneSeq；
- 卫星位置、速度、姿态；
- 用户位置、速度；
- local/global DoA 映射所需状态。

### reference-sat 语义

- 解析 reference satellite 的 local/global index；
- 构造 reference-sat Doppler state；
- scene / subset 裁剪后重映射 reference；
- 保持 static 与 dynamic 语义一致。

### Doppler 几何

- 构造 `fdRef`；
- 构造 `deltaFdRef`；
- 构造 `fdSat`；
- 计算相对 Doppler geometry。

### 用户状态

- 从 lat/lon 构造 motion-aware ECI user state；
- 保留 nominal state 与实际 scene 中额外用户运动偏移的一致性。

## 核心不变量

所有修改 reference-Doppler 或 subset scene 的代码，都应显式检查：

- `deltaFdRef(refSat) = 0`；
- `fdSat(refSat) = fdRef`；
- `fdSat = fdRef + deltaFdRef`。

对于 subset scene，还需要检查：

- 裁剪后 `refSatIdxLocal` 正确重映射；
- `refSatIdxGlobal` 不被偷换；
- sat 顺序与 estimator 看到的 local order 一致；
- `timeOffsetSec` 与 frame 裁剪保持一致。

## 主要文件

### `genMultiSatScene.m`

- 职责：构造多星单帧 scene。
- 调用者：static dev/regression、SF estimator。

### `genMultiFrameScene.m`

- 职责：构造多星多帧 sceneSeq。
- 调用者：dynamic dev/replay/regression。

### `getLinkParam.m`

- 职责：计算链路方向、距离、Doppler 等基础参数。
- 调用者：scene 构造、truth 构造。

### `getSceneSteering.m`

- 职责：从 scene 得到阵列 steering 相关量。
- 调用者：signal generation、estimator。

### `resolveReferenceSatState.m`

- 职责：解析 reference satellite local/global index 和 state。
- 调用者：static/dynamic estimator。

### `buildReferenceDopplerState.m`

- 职责：构造 `fdRef / deltaFdRef / fdSat` 等 reference-Doppler 状态。
- 调用者：estimator model build、regression。

### `computeRelativeDopplerGeom.m`

- 职责：计算几何差分 Doppler。
- 调用者：static/dynamic reference state。

### `buildUserStateFromLatlon.m`

- 职责：从 lat/lon 构造 motion-aware ECI user state。
- 调用者：estimator model、regression。

### `selectSatScene.m`

- 职责：单帧 scene 按卫星裁剪并重映射 reference。
- 调用者：ablation、subset regression。

### `selectSatSceneSeq.m`

- 职责：多帧 sceneSeq 按卫星裁剪并重映射 reference。
- 调用者：subset flow、replay。

### `selectFrameSceneSeq.m`

- 职责：多帧 sceneSeq 按 frame 裁剪。
- 调用者：periodic/random subset、scan。

## 修改风险点

### reference-sat 选择

修改 reference-sat 选择或 override 时，需要同时检查：

- local/global index；
- subset 后 reference remap；
- `fdRef` truth 口径；
- summary 是否拿的是 subset reference truth。

### `fdRef / deltaFdRef / fdSat` 构造

修改 Doppler 构造时，需要同时检查：

- ref sat 本身 delta 是否为 0；
- 非参考星差分 Doppler 是否与几何一致；
- static 与 dynamic 是否复用同一语义；
- summary 是否区分 raw-domain 与 eval-domain。

### frame 裁剪

修改 frame subset 时，需要同时检查：

- `timeOffsetSec` 是否保留原始时间；
- CP absolute-time 语义是否保留；
- IP frame-local time 是否未被误改；
- sceneSeq 和 rxSig 的 frame 顺序是否一致。

## 不应放在这里的内容

### estimator 参数打包和优化

放到：

- `estimator/helper/`。

### pilot / rx signal / snapshot 生成

放到：

- `satellite/signal/`。

### steering drift、Doppler line-fit truth 诊断

放到：

- `satellite/dynamic/`。

### regression fixture 组合

放到：

- `test/common/fixture/`。

### summary / report 表格

放到：

- `test/common/report/`；
- `test/common/summary/`。

## 相关 regression

- `regressionRefStateInvariant.m`：reference-Doppler 不变量；
- `regressionUserStateFromLatlon.m`：motion-aware user state；
- `regressionLocalDoaMapping.m`：local/global DoA 映射；
- `regressionSfStaticReferenceDopplerInvariant.m`：static reference Doppler 组合关系；
- `regressionMfCpIpTimeAxisInvariant.m`：CP/IP time-axis 语义。
