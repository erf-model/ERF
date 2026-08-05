# Validation Operators

These rules govern parity campaigns, retreat, ledgers, and acceptance decisions.

## Path and Lane Model
- `tag_frontier` answers implementation-parity questions at controlled scope.
- `plotfile_parity` answers evolved-input and physical-regression questions at milestone scope.
- Stop at the first materially failing frontier. Do not continue forward validation past an unexplained divergence.

## Canonical Validation Order
1. Stabilize crashes and gross runtime failures first.
2. Establish bridge-oracle behavior with tracked assets.
3. Use tag-frontier diagnostics to isolate the first failing block.
4. Escalate to line-level evidence only after block-level narrowing.
5. Move to plotfile milestones only after controlled-scope acceptance is in place.

## Retreat Rules
- Use Stop, Diff, Retreat:
  - stop on first material divergence
  - record the first divergent group/tag/step
  - mark downstream scope as retreat if earlier scope reopens
- Retreat is two-stage:
  - block-by-block narrowing
  - line-by-line proof with matched `PRE_<site>` and `POST_<site>` evidence
- Do not patch logic without first-line evidence at the scoped retreat site.

## Restart and Plotfile Rules
- Use shared-source one-step restart reproduction to classify first-fail plotfile divergence.
- Use the restart-state matrix before invoking later plotfile-to-tag linkage when source-state provenance is unclear.
- Serial-first scope remains the default until serial acceptance is complete.
- Treat milestone parity as a campaign with explicit stop/go gates, not as an ad hoc sequence of runs.
- Use variable-aware growth adjudication for epsilon-scale thermo/pressure residuals so numerical quantization is not over-classified.

## Diagnostic Contract
- Tier 1 canonical tags are the default forward-validation surface.
- Tier 1.5 block signatures sit between canonical tags and full line traces.
- Tier 2 forensic traces are required for proven line-level attribution.
- Host-side diagnostic prints after GPU work must synchronize first.
- Diagnostic field order is contractual once a trace schema is declared.

## Ledger Model
- `validation_manifest.tsv`: compact current frontier state, not append-only history.
- `validation_runs.tsv`: legacy combined run ledger.
- `executions.tsv`: execution facts only.
- `decisions.tsv`: manifest/frontier mutations citing execution evidence.
- `regression_cases.tsv`: scoped acceptance and failure obligations for known regressions.

## Evidence Discipline
- No formal manifest change without clean-SHA evidence.
- Dirty-SHA executions are acceptable for triage and debug only.
- Use tracked assets only for formal validation outcomes.
- Archive first-fail compare artifacts immediately and keep provenance attached to the decision record.
- Formal execution rows must use unique plot, checkpoint, and log roots per run so evidence is non-overwriting and restartable.
- Prefer dedicated ignored or stashed run directories for validation artifacts so clean-SHA reruns do not require destructive cleanup of prior evidence.
