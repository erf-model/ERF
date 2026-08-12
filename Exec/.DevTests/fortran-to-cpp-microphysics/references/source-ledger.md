# Source Ledger

This ledger records the source corpus and the major extracted items that must survive into the split package.

| source_id | source_file | source_loc | source_kind | category | status | destination |
| --- | --- | --- | --- | --- | --- | --- |
| SRC_FILE_RULEBOOK | `fortran_to_cpp_microphysics_skill.md` | file | markdown rulebook | vocabulary/rules | preserved | `core-rules.md`, `validation-operators.md`, `applications/wsm6.md` |
| SRC_FILE_IMPL | `wsm6_implementation_notes.md` | file | markdown notes | scheme note | preserved | `applications/wsm6.md`, `core-rules.md` |
| SRC_FILE_OPERATOR | `wsm6_validation_operator.md` | file | markdown operator | validation operator | preserved | `validation-operators.md` |
| SRC_FILE_SIDECAR | `wsm6_diagnostic_sidecar_notes.md` | file | markdown process | process operator | preserved | `process-operators.md` |
| SRC_FILE_AUDIT_SUMMARY | `rules_audit_summary.md` | file | markdown summary | evidence pointer | preserved | `applications/wsm6.md` |
| SRC_FILE_AUDIT_TSV | `rules_audit.tsv` | file | TSV inventory | schema | preserved | `source-materials.md`, `coverage-matrix.md` |
| SRC_FILE_MANIFEST | `validation_manifest.tsv` | file | TSV ledger | schema/state model | preserved | `validation-operators.md` |
| SRC_FILE_RUNS | `validation_runs.tsv` | file | TSV ledger | schema/state model | preserved | `validation-operators.md` |
| SRC_FILE_EXECUTIONS | `executions.tsv` | file | TSV ledger | schema/state model | preserved | `validation-operators.md` |
| SRC_FILE_DECISIONS | `decisions.tsv` | file | TSV ledger | schema/state model | preserved | `validation-operators.md` |
| SRC_FILE_REGRESSIONS | `regression_cases.tsv` | file | TSV ledger | schema/acceptance policy | preserved | `validation-operators.md` |
| SRC_FILE_INPUT_WSM6 | `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_WSM6` | file | tracked asset | evidence pointer | preserved | `applications/wsm6.md` |
| SRC_FILE_HARNESS_R600B | `Exec/CanonicalTests/SquallLine_2D/run_r600b_subtaskA.sh` | file | tracked asset | evidence pointer | preserved | `applications/wsm6.md` |
| SRC_FILE_WSM6_README | `Source/Microphysics/WSM6/README` | file | tracked asset | evidence pointer | preserved | `applications/wsm6.md` |
| RULE_PATH_A_MODEL | `fortran_to_cpp_microphysics_skill.md` | vocabulary/path model | markdown rulebook | vocabulary | preserved | `core-rules.md` |
| RULE_PHASE_MODEL | `fortran_to_cpp_microphysics_skill.md` | phase model | markdown rulebook | vocabulary | preserved | `core-rules.md` |
| RULE_LANE_MODEL | `fortran_to_cpp_microphysics_skill.md` | lane model | markdown rulebook | vocabulary | preserved | `validation-operators.md` |
| RULE_1 | `rules_audit.tsv` | `Rule_1` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_2 | `rules_audit.tsv` | `Rule_2` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_3 | `rules_audit.tsv` | `Rule_3` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_4 | `rules_audit.tsv` | `Rule_4` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_9 | `rules_audit.tsv` | `Rule_9` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_10 | `rules_audit.tsv` | `Rule_10` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_11 | `rules_audit.tsv` | `Rule_11` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_12 | `rules_audit.tsv` | `Rule_12` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_13 | `rules_audit.tsv` | `Rule_13` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_14 | `rules_audit.tsv` | `Rule_14` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_21 | `rules_audit.tsv` | `Rule_21` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_21A | `rules_audit.tsv` | `Rule_21_Addendum_ExactHelperSemantics` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_22 | `rules_audit.tsv` | `Rule_22` | markdown rule row | process operator | preserved | `process-operators.md` |
| RULE_23 | `rules_audit.tsv` | `Rule_23` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_24 | `rules_audit.tsv` | `Rule_24` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_24A | `rules_audit.tsv` | `Rule_24_Addendum_NISLFVBoundary` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_25 | `rules_audit.tsv` | `Rule_25` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_26 | `rules_audit.tsv` | `Rule_26` | markdown rule row | process operator | preserved | `process-operators.md` |
| RULE_27 | `rules_audit.tsv` | `Rule_27` | markdown rule row | scheme note | preserved | `applications/wsm6.md` |
| RULE_28 | `rules_audit.tsv` | `Rule_28` | markdown rule row | scheme note | preserved | `applications/wsm6.md` |
| RULE_29 | `rules_audit.tsv` | `Rule_29` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_30 | `rules_audit.tsv` | `Rule_30` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_30A | `rules_audit.tsv` | `Rule_30A` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_30B | `rules_audit.tsv` | `Rule_30B` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_30D | `rules_audit.tsv` | `Rule_30D` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_31 | `rules_audit.tsv` | `Rule_31` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_32 | `rules_audit.tsv` | `Rule_32` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_32A | `rules_audit.tsv` | `Rule_32A` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_33 | `rules_audit.tsv` | `Rule_33` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_34 | `rules_audit.tsv` | `Rule_34` | markdown rule row | rule | preserved | `core-rules.md` |
| RULE_35 | `rules_audit.tsv` | `Rule_35` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_36 | `rules_audit.tsv` | `Rule_36` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE36_MATRIX | `rules_audit.tsv` | `Rule36_RestartStateMatrixPolicy` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_GROWTH | `rules_audit.tsv` | `VariableAwareGrowthAdjudication` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_37 | `rules_audit.tsv` | `Rule_37` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| RULE_38 | `rules_audit.tsv` | `Rule_38` | markdown rule row | validation operator | preserved | `validation-operators.md` |
| IMPL_RULE9_SAVEVARS | `rules_audit_summary.md` | generic move candidates | markdown summary | rule | preserved | `core-rules.md` |
| IMPL_SCOP_GATE | `rules_audit_summary.md` | Appendix A scope gate | markdown summary | process operator | preserved | `process-operators.md` |
| IMPL_FCOMPARE | `rules_audit_summary.md` | Appendix A fcompare workflow | markdown summary | validation operator | preserved | `validation-operators.md` |
| IMPL_ARCHIVE_DIFF | `rules_audit_summary.md` | archive diff output | markdown summary | process operator | preserved | `process-operators.md` |
| IMPL_RESTART_REPRO | `rules_audit_summary.md` | restart reproducibility appendix | markdown summary | validation operator | preserved | `validation-operators.md` |
| IMPL_BOUNDED_REFINEMENT | `rules_audit_summary.md` | bounded per-step refinement | markdown summary | process operator | preserved | `process-operators.md` |
| IMPL_CAMPAIGN_NAMING | `wsm6_implementation_notes.md` | `:426` run root, `:444-452` build-type id convention | markdown notes | process operator | late extraction 2026-08-12 | `process-operators.md` |
