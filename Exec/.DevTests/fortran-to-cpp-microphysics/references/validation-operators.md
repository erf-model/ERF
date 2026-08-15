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
- Tier 2 forensic traces are required for proven line-level attribution, except
  under the Tier 1-R reconstruction clause below.

### Tier 1-R: reconstruction as line-level proof

A closed-form reconstruction from printed Tier 1 inputs may stand in for a Tier 2
trace, but only when every one of these holds. Reconstruction is cheaper than
instrumentation and does not perturb the build, so the conditions are what keep
it from becoming a general escape hatch.

- Inputs to the retreat site are bitwise equal at `PRE_<site>` across both legs,
  established over the full printed `k` range, not at a single cell.
- The reconstruction reproduces **both** legs bitwise, to the full precision of
  the print contract, under the two competing hypotheses. Matching only the leg
  you believe is wrong proves nothing; it is matching the oracle *and* the port
  from one arithmetic model that excludes alternatives.
- Every cell is accounted for. Cells that do not reproduce must be explained by
  a named mechanism, such as a clamp or a gate, and that mechanism must itself
  be shown to make the two legs agree.
- The divergence is a constant factor or another closed form that can be
  evaluated exactly outside the code. State-dependent or iterative divergence is
  not eligible; use Tier 1.5 or Tier 2.
- The reconstruction is recorded in `decisions.tsv` with the cell count, the
  columns and steps covered, and `first_line_divergence=<file>:<line>`.

Prefer Tier 1-R over Tier 2 when a scheme has no Tier 2 plumbing, rather than
treating the absence of the instrument as permission to patch on block-level
evidence alone. Two confirming cases in the WDM6 Bubble campaign: the `slope_*`
local `t0c` shadow (G4) and the `qck1` literal-product evaluation (G13a); in both
the planned Tier 2 decomposition was rendered unnecessary by reconstruction.
- Host-side diagnostic prints after GPU work must synchronize first.
- Diagnostic field order is contractual once a trace schema is declared.
- Measure a diagnostic surface by coverage per tag, never by tag count. Coverage
  has three independent axes and none of them implies the others: the `k` range
  a tag spans, the field set it prints, and whether the two legs emit a matched
  shape that can actually be compared. A scheme can carry three times as many
  tags as a well-instrumented one and still be strictly weaker on every axis.
- A tag that reads zero where its inputs cannot activate the block it guards is a
  NULL RESULT, not a closure. Establish that a group's inputs are non-null at the
  step and the level being sampled before reading agreement as evidence. Two
  forms recur and both must be checked:
  - temporal, when the case starts from rest and a gated branch has not yet
    received nonzero input at the sampled step;
  - spatial, when a single-cell tag samples a level at which the species that
    group consumes are identically zero.
- Prefer converting the remaining single-cell tags over opening a retreat on a
  downstream group, whenever an unconverted group sits upstream of the candidate
  and produces anything the candidate consumes. Attribution to a group that
  merely displays an inherited divergence is the characteristic failure of an
  incomplete tag surface.

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
