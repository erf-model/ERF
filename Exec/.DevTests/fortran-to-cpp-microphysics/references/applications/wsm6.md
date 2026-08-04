# WSM6 Application Note

## Role in This Skill
WSM6 is the primary source corpus for this package. It provides:
- the long-form rulebook
- scheme-specific implementation notes
- the validation operator
- sidecar/operator campaign notes
- structured TSV ledgers and audits

Primary in-package source corpus:
- `Exec/.DevTests/fortran-to-cpp-microphysics/references/source-corpus/`

Historical origin:
- extracted from the `wsm6-mpi-cpu` worktree on August 3, 2026

## Source Corpus
- `references/source-corpus/fortran_to_cpp_microphysics_skill.md`
- `references/source-corpus/wsm6_implementation_notes.md`
- `references/source-corpus/wsm6_validation_operator.md`
- `references/source-corpus/wsm6_diagnostic_sidecar_notes.md`
- `references/source-corpus/rules_audit.tsv`
- `references/source-corpus/rules_audit_summary.md`
- `references/source-corpus/validation_manifest.tsv`
- `references/source-corpus/validation_runs.tsv`
- `references/source-corpus/executions.tsv`
- `references/source-corpus/decisions.tsv`
- `references/source-corpus/regression_cases.tsv`

Tracked validation assets referenced by the corpus:
- `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_WSM6`
- `Exec/CanonicalTests/SquallLine_2D/run_r600b_subtaskA.sh`
- `Source/Microphysics/WSM6/README`

## Extraction Guidance
- `fortran_to_cpp_microphysics_skill.md` is the main rule and vocabulary source.
- `wsm6_implementation_notes.md` contains both generic lessons and WSM6-only details; promote cautiously.
- `wsm6_validation_operator.md` is the main source for validation-lane procedures and ledger discipline.
- `wsm6_diagnostic_sidecar_notes.md` is a process/operator source, not a core-rule source.
- `rules_audit.tsv` is the row-level inventory for rule extraction.
- `rules_audit_summary.md` is interpretive guidance and should not override the source corpus silently.
- `validation_manifest.tsv`, `validation_runs.tsv`, `executions.tsv`, `decisions.tsv`, and `regression_cases.tsv` define validation schemas and decision semantics that must survive the split.

## Current WSM6 State Captured by the Corpus
- Path A bridge is complete and used as oracle.
- Path B native C++ is advanced enough to support detailed retreat and plotfile campaigns.
- Validation history is rich enough that many rules are earned rather than hypothetical.
- SquallLine is the canonical tracked validation lane in the current corpus.

## Promotion Policy
- Promote reusable bridge/native/validation rules to `core-rules.md` or `validation-operators.md`.
- Keep WSM6 execution ordering, tag tables, milestone step numbers, and variable-to-tag mappings here unless confirmed generically elsewhere.
- Keep branch-history facts as evidence pointers only.
