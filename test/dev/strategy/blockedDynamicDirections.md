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

## Truth-dependent tooth selection or truth-aware adoption

- Status: blocked for the default path.
- Affected layer: subset tooth selection, rescue gate, candidate adoption, final winner selection.
- Observation: truth tooth, truth DoA, truth `fdRef/fdRate`, oracle labels, and manually assigned tail classes are useful for oracle replay, offline scan evaluation, and result annotation. Letting them influence the runtime selector would invalidate the no-truth-leak contract and weaken the paper/system boundary.
- Allowed future use: oracle replay, offline schedule evaluation, result labeling, and mechanism explanation only. Runtime gates must use data-derived metrics such as profile residual, non-ref coherence, phase residual, candidate consensus, or objective comparisons computed from the received data and known geometry.

## Blanket wide or single-MF basin-entry bank

- Status: blocked for the default path.
- Affected layer: same-tooth rescue / final periodic refinement.
- Observation: in-tooth tail replay shows that wide-centered and single-MF-centered candidates can rescue different non-ref coherence collapse hard cases, but blanket use can damage easy or fd-negative cases. The viable route is a truth-free collapse gate, not unconditional adoption.
- Allowed future use: flow-like replay comparisons with disabled, gated-wide-only, gated-single-MF-only, gated-wide-single-bank, and blanket-damage controls. Default-flow adoption requires evidence that hard cases improve and easy / fd-negative cases are not damaged.

## Blanket structured nonuniform subset schedules

- Status: blocked for the default path until adoption evidence improves.
- Affected layer: subset bank / tooth selection.
- Observation: structured schedules with better lag features can improve candidate coverage without improving selected tooth hit, and can introduce easy-case damage. The current evidence points to candidate-to-final adoption and validation, not simply schedule count or lag design.
- Allowed future use: offline scan or strategy comparisons. A structured schedule can only move toward the default bank after it improves selected hit / fd tails under no-truth selection and does not increase easy damage at acceptable candidate cost.

## Expanding full-flow CP/IP performance maps before flow stabilization

- Status: deferred.
- Affected layer: paper performance scan production.
- Observation: current full-flow CP/IP performance map is polluted by wrong-tooth and same-tooth bad-basin behavior. Enlarging repeats would mostly stabilize a flow failure result rather than measure the CP/IP model hierarchy.
- Allowed future use: rerun after selected-tooth flow and same-tooth gated rescue are stable. Until then, use controlled in-tooth CP/IP scans for the angle-vs-Doppler-consistency trade-off and keep full-flow scans as stress-test evidence.
