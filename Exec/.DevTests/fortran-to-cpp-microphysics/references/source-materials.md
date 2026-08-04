# Source Materials

Primary in-package source corpus:
- `Exec/.DevTests/fortran-to-cpp-microphysics/references/source-corpus/`

Historical origin:
- extracted from the `wsm6-mpi-cpu` worktree on August 3, 2026

## Markdown Sources
- `source-corpus/fortran_to_cpp_microphysics_skill.md`
  - long-form rulebook, vocabulary, path/phase model, and skill source
- `source-corpus/wsm6_implementation_notes.md`
  - WSM6 implementation notes, draft-to-final rule refinements, and scheme-local details
- `source-corpus/wsm6_validation_operator.md`
  - validation operator workflow, ledger discipline, and non-negotiable campaign rules
- `source-corpus/wsm6_diagnostic_sidecar_notes.md`
  - process/operator notes, debug discipline, harness pointers, and evidence practices
- `source-corpus/rules_audit_summary.md`
  - curated summary of what is generic, what should move, and what remains scheme-local

## TSV Sources
- `source-corpus/rules_audit.tsv`
  - structured row-level inventory of rulebook content and recommended homes
- `source-corpus/validation_manifest.tsv`
  - compact current-state frontier table
- `source-corpus/validation_runs.tsv`
  - legacy combined validation ledger
- `source-corpus/executions.tsv`
  - split execution-facts ledger
- `source-corpus/decisions.tsv`
  - split decision ledger for manifest/frontier mutations
- `source-corpus/regression_cases.tsv`
  - scoped regression cases, acceptance criteria, and linked-rule obligations

## Related Tracked Assets
- `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_WSM6`
- `Exec/CanonicalTests/SquallLine_2D/run_r600b_subtaskA.sh`
- `Source/Microphysics/WSM6/README`

## Materiality Rule
- Markdown files contribute prose rules, operators, process notes, and scheme context.
- TSV files contribute schemas, state models, acceptance logic, and decision discipline.
- Related tracked assets contribute canonical validation harness context.
- All three categories are part of the source-of-truth corpus for this skill package.
