# doa-tools

A MATLAB toolbox for array processing, DoA estimation, and the current multi-satellite DoA-Doppler research code.

This repository contains two layers of code:

1. **Generic array/estimation utilities** such as array construction, steering generation, classical estimators, and CRB helpers.
2. **Project-specific multi-satellite DoA-Doppler code** for static and dynamic LEO uplink scenarios, including scene construction, signal generation, MLE/CRB implementations, regression scripts, and development entry scripts.

The project has gradually evolved from a generic DOA toolbox into a mixed repository with both reusable toolbox code and active research code. The folder rules below are intended to keep that boundary clear and reduce maintenance cost.

## Quick start

Run `startup.m` to add the repository to the MATLAB path.

For active project development, the usual entry points are:

- `test/dev/doaDopplerStatDualSatUraEci.m`
- `test/dev/doaDopplerDynDualSatUraEci.m`
- `test/regression/runRegressionQuick.m`
- `test/regression/runRegressionPipeline.m`

If the regression runners are not yet committed in your local branch, add them under `test/regression/` and keep their shared orchestration helper out of the estimator core.

## Top-level folder structure

### `array/`
Generic array construction and steering-related utilities.

Put code here when it is **array-topology specific** and does not depend on the satellite/scene model.

Typical contents:

- array geometry builders such as ULA/UPA creation
- coordinate-frame array rotation helpers
- generic snapshot generation for array-only experiments
- generic steering matrix helpers

Guideline:

- New reusable array helpers go here.
- Old prototype or compatibility implementations stay under `array/legacy/`.

### `estimator/`
Active estimators and estimator-side helper logic.

This folder contains both generic estimators and the project's current DoA-Doppler MLE implementations.

Typical contents:

- single-frame and multi-frame MLE estimators
- MUSIC / MVDR and related estimators
- shared profile-likelihood evaluators
- grid construction and local refinement helpers
- estimator-specific orchestration helpers in `estimator/helper/`

Use `estimator/helper/` for logic that is **part of the formal estimator pipeline** and may be reused by multiple estimators or branches, for example:

- model construction
- parameter packing/unpacking
- initialization construction
- branch solving
- profile-likelihood helper logic

Do **not** leave those pieces buried in `test/dev` local functions once they become stable.

Subfolders:

- `estimator/doaGrid/`: grid generation and neighborhood helpers
- `estimator/covariance/`: covariance-side estimator utilities
- `estimator/helper/`: active shared estimator helpers
- `estimator/legacy/`: frozen legacy implementations kept for reference or compatibility

### `performance/`
Performance analysis and CRB-related code.

Put code here when it computes **bounds, FIM/EFIM, or theoretical performance quantities**, not when it runs Monte Carlo experiments.

Typical contents:

- single-frame DoA / DoA-Doppler CRB
- multi-frame dynamic CRB
- helper logic for FIM blocks and shared CRB calculations

Guideline:

- If the code produces a bound from a model, it belongs here.
- If the code runs simulations to compare estimators, it belongs under `test/`.
- Compatibility wrappers for old CRB interfaces may remain, but the main logic should converge into the current unified entry points.

### `satellite/`
Satellite-specific geometry, scene, and signal-generation code.

This folder should contain code that depends on orbital geometry, reference-satellite semantics, or scene construction.

Recommended placement inside this folder:

#### `satellite/scene/`
Scene construction and reference-state logic.

Put code here when it defines or manipulates:

- satellite/user geometry
- reference satellite semantics
- scene slicing by satellite or by frame
- local/global DoA mapping tied to scene state
- Doppler and differential Doppler geometric construction
- motion-aware user-state construction

Examples already in the repo include:

- `resolveReferenceSatState`
- `buildReferenceDopplerState`
- `buildUserStateFromLatlon`
- `computeRelativeDopplerGeom`
- `selectSatScene`, `selectSatSceneSeq`, `selectFrameSceneSeq`

#### `satellite/signal/`
Signal and snapshot generation.

Put code here when it builds or modifies the received signal itself, for example:

- pilot/synchronization waveform generation
- multi-satellite or multi-frame snapshot generation
- path gain construction
- fractional delay application
- received-signal slicing helpers

Examples already in the repo include:

- `genPilotWaveform`
- `genMultiSatSnapshots`
- `genMultiFrameSnapshots`
- `buildFramePathGain`
- `selectRxSigBySat`

#### `satellite/dynamic/`
Dynamic diagnostics and dynamic truth construction.

Put code here when it is specific to **time evolution over the multi-frame window**, for example:

- steering drift diagnostics
- Doppler line-fit truth construction
- dynamic-strength scaling for experiments
- dynamic truth summary helpers

Examples already in the repo include:

- `buildDynTruthFromLinkParam`
- `scaleLinkParamDoppler`
- `buildSteeringDriftDiag`
- `buildDeltaFdFitDiag`

### `plottool/`
Reusable plotting helpers.

Put code here only when the plotting logic is generic enough to be shared across multiple scripts.

Do not move every one-off debug figure here. If a plot is only used by one development script, keep it local until it becomes genuinely reusable.

### `solvers/`
External solver interfaces and solver-side linear algebra helpers.

Use this folder for:

- solver capability checks
- SDP-related conversion helpers
- third-party solver glue that is not estimator-specific

### `utils/`
Low-level reusable utilities that do not belong to the signal model or estimator model.

Typical contents:

- math helpers
- structure utilities
- generic runners
- progress utilities
- metric helpers
- I/O helpers

Recommended placement:

- `utils/metric/`: scalar metrics such as DoA/lat-lon error calculations
- `utils/io/`: save/load helpers for experiment snapshots and artifacts

Rule of thumb:

- If a function has no DoA-Doppler model semantics and would still make sense in another project, it likely belongs in `utils/`.
- If it knows about reference satellites, frames, or profile likelihood, it probably belongs elsewhere.

### `test/`
Experiment, development, and regression code.

This folder is intentionally separated from `estimator/` and `performance/`: it should contain scripts and helpers that **exercise** the formal modules, not replace them.

Current subfolders and their intended roles are:

#### `test/common/`
Shared experiment-side helpers.

Put code here when it is reused by multiple test/dev/regression scripts but is **not** part of the formal estimator core.

Typical contents:

- case/result struct builders
- experiment-side summary tables
- static-to-dynamic transition helpers
- subset fixture builders
- multi-start summaries
- experiment-side plotting helpers
- regression/test path helpers

Examples already in the repo include:

- `buildDoaDopplerCaseResult`
- `buildDoaDopplerStaticTransitionBundle`
- `buildDynamicFrameSubsetFixture`
- `buildDynamicSubsetFixtureBank`
- `runDynamicDoaDopplerCase`
- `summarizeDynamicEstimatorDebug`
- `summarizeDynamicMultiStart`

What should **not** go here:

- core estimator mathematics that should live in `estimator/helper/`
- scene/reference-state construction that should live in `satellite/scene/`
- raw signal generation that should live in `satellite/signal/`

#### `test/dev/`
Development entry scripts and interactive diagnosis scripts.

Use this folder for scripts that you run manually while developing or debugging.

Typical contents:

- single-shot verification scripts
- scenario scripts used to inspect summaries and figures
- scan/sweep/trace style scripts for mechanism diagnosis
- ranking or screening scripts for satellite pairs or schedules

Examples already in the repo include:

- `doaDopplerStatDualSatUraEci`
- `doaDopplerDynDualSatUraEci`
- `doaDopplerDynFdRefCombScan`
- `doaDopplerDynTruthNeighborhoodScan`
- `doaDopplerDynCpIpTyingScan`

Rule of thumb:

- `test/dev/` scripts should mainly do orchestration, call helpers, and print results.
- Once a local function becomes long, reusable, or numerically important, move it out to a proper helper folder.

#### `test/regression/`
Focused regression scripts that guard stable behavior.

Only keep scripts here if they check a **clear contract** and are suitable for repeated regression runs.

Recommended categories, even if the folder is still flat in your branch:

- **Invariant** regressions: geometry, reference-state, model-build, truth-replay
- **Branch** regressions: known/unknown branch behavior, warm-start construction, candidate selection
- **Pipeline** regressions: the active end-to-end flow that should not silently regress

Typical examples in the repo:

- model-build regressions
- reference-state regressions
- unknown-rate warm-start / release regressions
- subset-tooth-selection pipeline regressions

Keep the runner scripts here as well:

- `runRegressionQuick.m`: fast invariant + branch checks
- `runRegressionPipeline.m`: active pipeline checks

What should **not** stay here:

- one-off trace scripts
- large exploratory sweeps without clear pass/fail semantics
- dev-only visualization scripts

Those belong in `test/dev/`.

#### `test/paper/`
Paper-figure generation scripts and experiment wrappers intended for final figures/tables.

Use this folder when a script is no longer only for diagnosis and is meant to generate stable paper-ready outputs.

#### `test/data/`
Small test data, saved fixtures, or local experiment assets needed by regression/dev scripts.

Keep only lightweight project data here. Large temporary dumps should stay outside the repo.

#### `test/archive/`
Retired or historical scripts.

Use this only for scripts that you want to preserve for reference but do not want to treat as active development or regression code.

Do not move active-but-messy code here just to avoid organizing it. Prefer a clear active folder whenever the script still has value.

### `examples/`
Generic examples inherited from the original toolbox.

These are useful for toolbox-level sanity checks and public-facing demos, but they are not the main home of the current multi-satellite project workflow.

## Placement rules for new code

Use the following rules when adding new code.

### Put code in `estimator/helper/` if it is:

- part of the formal estimator path
- reused across multiple estimator branches
- numerically important enough to deserve direct regression
- about model construction, initialization, branch solving, or profile likelihood

### Put code in `satellite/scene/` if it is about:

- scene construction
- reference-satellite semantics
- local/global DoA mapping
- user/satellite state construction
- Doppler geometry derived from scene state

### Put code in `satellite/signal/` if it is about:

- waveform generation
- snapshot generation
- path gains
- delay/Doppler application to signals
- received-signal slicing

### Put code in `performance/` if it is about:

- CRB / FIM / EFIM
- theoretical bounds
- analytic performance metrics derived from a model

### Put code in `test/common/` if it is:

- shared by multiple scripts
- experiment-side orchestration or reporting
- not part of the formal estimator kernel

### Keep code local in a script only if it is:

- short
- file-private
- glue only
- unlikely to be reused

Examples of acceptable script-local functions:

- compact formatting helpers
- tiny guards
- simple plot cosmetics
- default-value shims

Do **not** keep the following buried as long local functions once they stabilize:

- reference-Doppler semantics
- truth construction
- subset fixture generation
- path gain construction
- profile-likelihood logic
- unknown-rate release logic
- branch selection logic

## Legacy code policy

Folders named `legacy/` are retained for reference, compatibility, or historical comparison.

Guidelines:

- Do not add new mainline functionality there.
- Do not copy legacy code into new project code unless there is a deliberate compatibility reason.
- When a legacy entry point must remain callable, prefer a thin wrapper around the current implementation.

## Practical workflow

A practical development flow for this repo is:

1. Build or modify formal modules under `satellite/`, `estimator/`, or `performance/`.
2. Wire them into `test/common/` experiment helpers if needed.
3. Verify behavior with `test/dev/` scenario scripts.
4. Guard the stable behavior with `test/regression/`.
5. Only after the path is stable, move figure-generation wrappers to `test/paper/`.

## Notes for this project

For the current multi-satellite DoA-Doppler line, keep the following discipline:

- Use `test/dev/doaDopplerStatDualSatUraEci.m` and `test/dev/doaDopplerDynDualSatUraEci.m` as scenario entry scripts, not as long-term homes for core logic.
- Keep static and dynamic implementations aligned through shared helpers whenever the semantics are the same.
- Treat reference-satellite semantics as explicit invariants.
- Keep debug/probe logic on a bypass path so it does not alter the main objective or default outputs.
- Use regression scripts to guard invariants, branch behavior, and the active tooth-selection / in-tooth-refine pipeline.

## License

The source code is released under the [MIT](LICENSE.md) license.
