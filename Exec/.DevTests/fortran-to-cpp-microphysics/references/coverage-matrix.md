# Coverage Matrix

This matrix records how the source corpus was carried into the split package.

| source_id | destination_section | preservation_mode | orchestrator_mandatory | rationale |
| --- | --- | --- | --- | --- |
| SRC_FILE_RULEBOOK | `core-rules`, `validation-operators`, `applications/wsm6` | split across destinations | yes | The rulebook is too large for direct use but remains the primary source for rules and vocabulary. |
| SRC_FILE_IMPL | `applications/wsm6`, `core-rules` | split across destinations | yes | Generic lessons were promoted; WSM6-only structures remain scheme-local. |
| SRC_FILE_OPERATOR | `validation-operators` | near-verbatim | yes | Operator logic and ledger discipline are first-class. |
| SRC_FILE_SIDECAR | `process-operators` | normalized wording | yes | Converted branch/campaign notes into reusable process discipline. |
| SRC_FILE_AUDIT_SUMMARY | `applications/wsm6` | normalized wording | no | Used as interpretive support, not as sole authority. |
| SRC_FILE_AUDIT_TSV | `source-materials`, `coverage-matrix` | normalized wording | yes | The TSV inventory anchors provenance and extraction coverage. |
| SRC_FILE_MANIFEST | `validation-operators` | near-verbatim | yes | Frontier-state model is essential operator context. |
| SRC_FILE_RUNS | `validation-operators` | normalized wording | yes | Legacy ledger semantics preserved as migration context. |
| SRC_FILE_EXECUTIONS | `validation-operators` | near-verbatim | yes | Execution-facts split is part of the formal operator model. |
| SRC_FILE_DECISIONS | `validation-operators` | near-verbatim | yes | Decision evidence discipline must survive intact. |
| SRC_FILE_REGRESSIONS | `validation-operators` | normalized wording | yes | Acceptance/failure obligations were preserved as validation rules. |
| SRC_FILE_INPUT_WSM6 | `applications/wsm6` | evidence pointer | no | Tracked canonical asset for WSM6 campaign context. |
| SRC_FILE_HARNESS_R600B | `applications/wsm6` | evidence pointer | no | Harness provenance retained without overloading the core skill. |
| SRC_FILE_WSM6_README | `applications/wsm6` | evidence pointer | no | Run/compare workflow retained as tracked support material. |
| RULE_PATH_A_MODEL | `core-rules` | near-verbatim | yes | The path vocabulary is part of the operating model. |
| RULE_PHASE_MODEL | `core-rules` | near-verbatim | yes | The phase model remains central to orchestrated porting. |
| RULE_LANE_MODEL | `validation-operators` | near-verbatim | yes | Validation lane terminology is contractual. |
| RULE_1 | `core-rules` | near-verbatim | yes | Core reusable translation pattern. |
| RULE_2 | `core-rules` | near-verbatim | yes | Explicit pack/unpack remains required. |
| RULE_3 | `core-rules` | normalized wording | yes | Shortened but obligation strength preserved. |
| RULE_4 | `core-rules` | near-verbatim | yes | High-risk numerical drift guardrail. |
| RULE_9 | `core-rules` | near-verbatim | yes | Init/runtime separation preserved. |
| RULE_10 | `core-rules` | normalized wording | yes | Variant-flag principle preserved. |
| RULE_11 | `core-rules` | normalized wording | yes | Translation-strategy choices preserved. |
| RULE_12 | `validation-operators` | near-verbatim | yes | Oracle policy is non-negotiable. |
| RULE_13 | `core-rules` | normalized wording | yes | 2D slab allocation lesson preserved. |
| RULE_14 | `core-rules` | normalized wording | yes | Init contract retained in compact form. |
| RULE_21 | `core-rules` | near-verbatim | yes | Exact helper semantics remain required. |
| RULE_21A | `core-rules` | near-verbatim | yes | Runtime helper parity preserved. |
| RULE_22 | `process-operators` | normalized wording | yes | Convention survey is process discipline, not physics law. |
| RULE_23 | `core-rules` | normalized wording | yes | Translate intent, not dead scaffolding. |
| RULE_24 | `core-rules` | normalized wording | yes | Generic sedimentation concept preserved. |
| RULE_24A | `validation-operators` | normalized wording | yes | Boundary-proof retreat rule preserved. |
| RULE_25 | `core-rules` | near-verbatim | yes | The five-phase workflow remains central. |
| RULE_26 | `process-operators` | normalized wording | yes | Context-budget discipline is process-level guidance. |
| RULE_27 | `applications/wsm6` | near-verbatim | yes | WSM6 execution ordering remains scheme-local. |
| RULE_28 | `applications/wsm6` | near-verbatim | yes | WSM6 working-array layout remains scheme-local. |
| RULE_29 | `core-rules` | normalized wording | yes | Minor-timestep outer-loop concept preserved. |
| RULE_30 | `validation-operators` | near-verbatim | yes | Phase 3 diagnostic protocol is contractual. |
| RULE_30A | `validation-operators` | near-verbatim | yes | Runtime high-precision diagnostics preserved. |
| RULE_30B | `validation-operators` | near-verbatim | yes | Tier2 invocation criteria preserved. |
| RULE_30D | `validation-operators` | near-verbatim | yes | Tier1.5 block-signature layer preserved. |
| RULE_31 | `core-rules` | normalized wording | yes | Runtime toggle rule preserved without Morrison-only framing. |
| RULE_32 | `validation-operators` | normalized wording | yes | Debug-first validation order preserved. |
| RULE_32A | `validation-operators` | near-verbatim | yes | Stop, Diff, Retreat remains hard requirement. |
| RULE_33 | `core-rules` | near-verbatim | yes | Pre-compile alignment pass preserved. |
| RULE_34 | `core-rules` | normalized wording | yes | Control-flow hazard policy preserved. |
| RULE_35 | `validation-operators` | normalized wording | yes | Plotfile campaign structure preserved without WSM6 step numbers. |
| RULE_36 | `validation-operators` | near-verbatim | yes | Shared-source restart reproduction preserved. |
| RULE36_MATRIX | `validation-operators` | near-verbatim | yes | Restart-state matrix policy preserved. |
| RULE_GROWTH | `validation-operators` | near-verbatim | yes | Variable-aware growth adjudication preserved. |
| RULE_37 | `validation-operators` | normalized wording | yes | Plotfile-to-tag linkage preserved without WSM6-only mappings. |
| RULE_38 | `validation-operators` | near-verbatim | yes | TSV schema extension concept preserved. |
| IMPL_RULE9_SAVEVARS | `core-rules` | normalized wording | yes | Save-variable init lesson promoted to generic init-separation guidance. |
| IMPL_SCOP_GATE | `process-operators` | near-verbatim | yes | Serial-first scope gate preserved. |
| IMPL_FCOMPARE | `validation-operators` | normalized wording | yes | First-fail compare workflow preserved. |
| IMPL_ARCHIVE_DIFF | `process-operators` | near-verbatim | yes | Immediate diff-archive habit preserved. |
| IMPL_RESTART_REPRO | `validation-operators` | normalized wording | yes | Appendix-specific WSM6 phrasing collapsed into generic restart operator guidance. |
| IMPL_BOUNDED_REFINEMENT | `process-operators` | normalized wording | yes | Bounded refinement remains operator guidance. |
