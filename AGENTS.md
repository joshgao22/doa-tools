# AGENTS.md

This file is the shortest entry point for AI or automated code changes. It defines the required reading order, document responsibilities, conflict resolution rules, default modification principles, and delivery format. It does not replace directory-level README files and does not record experiment results.

## Pre-change reading order

1. Read this file.
2. Read the repo-root `README.md` to confirm global code style, comment language, storage conventions, and directory responsibilities.
3. Read the README in the target file's directory.
4. If the target file is in a deeper directory, continue reading the parent and current directory README files as needed.
5. If modifying replay / scan / strategy / regression files, read the corresponding local README.
6. If modifying estimator helpers, read `estimator/README.md` and `estimator/helper/README.md`.
7. If modifying scene logic, reference-satellite handling, Doppler geometry, or scene slicing, read `satellite/scene/README.md`.
8. If modifying dynamic flow, subset tooth selection, conditional random rescue, same-tooth polish, warm-anchor logic, or related replay files, read the current main troubleshooting record. Read the mechanism summary and historical archive only when mechanism background or old-line traceback is needed.

## Document responsibilities

- `AGENTS.md`: shortest reading path, priority rules, conflict resolution, default modification principles, and delivery format.
- Repo-root `README.md`: repository-level directory responsibilities, code style, storage convention index, and documentation synchronization rules.
- Target directory README: local entry points, file responsibilities, local run instructions, and local code constraints.
- `results/` documents: concrete replay / scan / strategy run results, long tables, snapshot bindings, and observations.
- Main troubleshooting record: current main problem, current priority, required guardrails, and disproven routes.
- Mechanism summary / historical archive: mechanism evidence and historical traceback; read only when understanding why or tracing old clues is necessary.
- External coding prompt: AI execution style, risk preference, and delivery explanation format; it is not a substitute for repository file indexes.

## Priority and conflict resolution

1. Current explicit user request.
2. This file.
3. Repo-root `README.md`.
4. Target directory README.
5. Deeper directory README.
6. Current main troubleshooting record.
7. Mechanism summary / historical archive.

If the user request, external prompt, or README files conflict, choose the smaller, more local, and more behavior-preserving option. If a conflict touches function signatures, default parameters, reference-satellite semantics, search bounds, satellite order, output fields, summary semantics, or fallback order, do not hide the risk in the implementation. State the risk and trade-off clearly in the delivery notes.

## Default code strategy

- Preserve behavior by default. Do not change numerical paths unless explicitly required.
- Do not mechanically follow a temporary proposal. If a proposed change would increase interface bloat, duplicated documentation, context burden, reintroduce a disproven route, increase default numerical risk, or raise maintenance cost, choose a smaller and safer alternative.
- Do not change function signatures, default parameters, initialization strategy, search bounds, reference-satellite semantics, satellite order, output fields, summary semantics, or fallback order without explicit need.
- Function headers, script headers, and inline comments in `.m` files must be written in English.
- README files, troubleshooting records, and explanatory documents are normally written in Chinese unless the user explicitly asks otherwise.
- Follow the style of non-legacy files of the same type in the target directory.
- Do not add heavyweight `opt` fields, `arguments` blocks, `validateattributes`, configuration resolvers, or framework-style wrappers for small changes.
- If a change can be absorbed in upper-level orchestration, do not push new interfaces down into lower-level helpers.
- Formal estimator logic belongs in `estimator/` or `estimator/helper/`; scene / reference-satellite / Doppler geometry belongs in `satellite/scene/`; CRB / FIM / EFIM belongs in `performance/`; experiment orchestration and summaries belong in `test/`.

## Test subsystem boundaries

- `test/regression/`: automatic pass/fail contract guardrails.
- `test/dev/replay/`: fixed-seed or small Monte Carlo replay for recording phenomena and mechanism evidence.
- `test/dev/scan/`: heavier parameter, schedule, curve, or surface scans.
- `test/dev/strategy/`: strategy comparison, disproven routes, and retry admission conditions.
- `test/common/`: experiment-side reusable helpers; do not maintain a second formal algorithm implementation here long term.

## Long-task notification boundary

- Long-running replay / scan / Monte Carlo scripts may call `utils/io/notifyTelegram.m` at the top-level orchestration layer to send completion or failure notifications.
- Notification logic must be best-effort: a send failure may only issue a `warning`; it must not affect simulation results, snapshot saving, `replayData` / `scanData`, or estimator numerical paths.
- `notifyTelegram` is a generic I/O utility only. Do not push it into `estimator/`, `performance/`, objective functions, residual functions, CRB code, or formal algorithm helpers.
- Tokens, chat IDs, and local machine configuration must not be committed to the repository. Prefer environment variables or private MATLAB `prefdir` configuration.
- When adding or modifying long-running replay / scan scripts, a short notification may be added at script completion, failure `catch`, or after snapshot saving. The message should report only the script name, status, elapsed time, snapshot path, and a small key summary. Do not maintain a second results parser in notifications.

## Results, storage, and documentation boundaries

- Runtime temporary files go under repo-root `tmp/<scriptName>/<runKey>/`.
- Large `.mat` snapshots go under `test/data/cache/<taskType>/`.
- Detailed replay run results go under `test/dev/replay/results/`.
- Detailed scan run results go under `test/dev/scan/results/`.
- Long strategy results go under `test/dev/strategy/results/`.
- Detailed performance results go under the corresponding performance results document.
- Save, load, and cleanup entry points are indexed only by `test/data/cache/README.md`. Do not duplicate save / load / delete details across multiple README files.
- README files should keep only entry points, responsibilities, run instructions, and output explanations. They should not record long results.
- Troubleshooting records should keep only current priorities, mechanism conclusions, and disproven routes. Do not copy README content or full run logs into troubleshooting records.

## Modification delivery format

When the user asks for direct code modification and asks to receive the changed files, the default delivery is a zip archive.

- The zip archive must preserve repository-relative paths.
- The archive must include only files added or modified in the current change.
- The archive must not include unchanged files.
- The archive must not include `.git/`, `tmp/`, cache directories, large `.mat` files, generated figures, or unrelated temporary files.
- The user should be able to unzip the archive at the repository root and overwrite files directly.
- The final response must list:
  - added / modified files;
  - why the files belong in those paths;
  - whether any interface changed;
  - whether any default behavior changed;
  - which local helpers were kept and why;
  - which logic, if any, was extracted into helpers and why;
  - recommended regression / replay cases to run.
