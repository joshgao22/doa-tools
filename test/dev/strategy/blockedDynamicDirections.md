# Blocked dynamic directions

This file records directions that have already been tried or diagnosed as high-risk for the current DoA-Doppler dynamic flow. It is not a regression file. Entries here should prevent repeated attempts unless new evidence changes the premise.

## Coarse-DoA joint `[fdRef, fdRate]` initialization

- Status: blocked for the default path.
- Affected layer: estimator initialization.
- Observation: using the final dynamic profile likelihood to search `[fdRef, fdRate]` while DoA is still coarse couples DoA error, differential-Doppler mapping error, and CP profile sensitivity. It can damage even the single-satellite dynamic baseline.
- Allowed future use: dev/probe only, after DoA is already tightly anchored and the default initialization path is unchanged.

## Alias-prior refine in the main path

- Status: blocked for the default path.
- Affected layer: dynamic initialization / branch selection.
- Observation: alias-aware quantities are useful diagnostics, but when alias-prior refine is allowed to alter the main initialization path it can push otherwise healthy known-rate cases onto the wrong branch.
- Allowed future use: read-only diagnostics or explicit dev replay/probe comparisons.

## Explicit `fdRate = 0` or global coarse anchors for CP-U

- Status: blocked for the default path.
- Affected layer: unknown-rate warm start.
- Observation: these anchors are attractive wrong basins in CP-U. They can keep the unknown-rate branch near zero rate or a coarse anchor instead of performing local release from a physically meaningful seed.
- Allowed future use: historical replay only, if needed to explain why the route was dropped.

## Blanket DoA release or blanket anchor-DoA polish

- Status: blocked for the default path.
- Affected layer: same-tooth periodic refinement.
- Observation: current evidence favors a conditional very-small DoA polish only when fd is healthy, the selected tooth is stable, and same-tooth DoA remains suspicious. Blanket release can move easy/median cases into worse basins.
- Allowed future use: dev strategy comparison with explicit before/after tables; not a default-flow change without a new regression contract.

## Truth override or alias-aware fields in the main objective/residual

- Status: blocked for the default path.
- Affected layer: evaluator and summary semantics.
- Observation: truth override and alias-aware fields are useful to localize failures, but feeding them back into the main objective/residual mixes diagnostics with model semantics and can corrupt the estimator path.
- Allowed future use: probe/replay-only diagnostic evaluation that never changes winner selection.

## Warm-anchor inner `parfor` as estimator default

- Status: blocked for the default path.
- Affected layer: unknown warm-anchor release seeds.
- Observation: release-seed candidates are independent in principle, but worker-side numerical differences can change sensitive tooth guards and final winners. The estimator default must remain serial.
- Allowed future use: explicit opt-in dev/perf runs, followed by replay/regression checks of the wrong-tooth sentinel.

## Treating `otherSat(g2)` bookkeeping as the first-priority estimator bug

- Status: deprioritized.
- Affected layer: summary/reporting.
- Observation: several `otherSat(g2)` large-error rows were traced to subset-reference truth bookkeeping and summary fallback paths, not necessarily to the estimator core.
- Allowed future use: fix reporting glue when it blocks interpretation, but do not let it displace current flow-layer same-tooth work.

## Treating CRB singular / pinv warnings as the current main fault

- Status: deprioritized.
- Affected layer: CRB diagnostics.
- Observation: known/unknown CRB consistency warnings can be expected for selected full-FIM configurations. They are useful diagnostics but are not the current blocker for dynamic flow correctness.
- Allowed future use: local CRB/FIM diagnostics and targeted invariant regression only.
