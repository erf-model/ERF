# WSM6 Validation Corpus Audit Summary

## 1. True generic rules
The following rules apply to any Fortran-to-C++ microphysics port and should survive unchanged into a Morrison or Thompson conversion. They are selected from `GENERIC_METHOD_RULE`, `GENERIC_MICROPHYSICS_CONCEPT`, and `VALIDATION_POLICY` rows with `applies_to=any_scheme`.

Rules currently in `fortran_to_cpp_microphysics_skill.md`:
  - `Section_Status_PARTIAL` (`# Status: PARTIAL — validated on Morrison, pending MYNN 2.5 + WSM6`): Declares maturity state and outstanding scheme validation scope for the skill.
  - `Rule_1` (`## Rule 1: Two-Layer Enum Pattern [COMPLETE]`): Use a two-layer enum strategy to preserve Fortran semantic grouping while enabling C++ kernel structure.
  - `Rule_2` (`## Rule 2: Pack/Unpack Bridge Pattern [COMPLETE]`): Define explicit pack/unpack conventions at the Fortran bridge boundary so path parity remains testable.
  - `Rule_4` (`## Rule 4: Physical Constants Resolution Hierarchy [COMPLETE]`): Resolve physical constants via an explicit precedence hierarchy to avoid silent drift across paths.
  - `Rule_9` (`## Rule 9: Separate Coefficient Init [COMPLETE — confirmed Morrison, WSM6, MYNN]`): Keep coefficient initialization logically separated so runtime kernels stay deterministic and auditable.
  - `Rule_12` (`## Rule 12: Fortran Bridge as Ground Truth During Native Port [COMPLETE]`): Treat the Fortran bridge path as oracle while incrementally validating the native C++ implementation.
  - `Rule_22` (`## Rule 22: Name New Enum Entries by Convention Survey Before Writing [process rule]`): Survey existing naming conventions first, then add enum entries that preserve local codebase semantic patterns.
  - `Rule_23` (`## Rule 23: Dead and Simplified Fortran Constructs — Omit or Inline in C++ Port`): Identify dead/simplifiable Fortran constructs and translate intent, not unnecessary legacy structure, into C++.
  - `Rule_24` (`## Rule 24: nislfv Semi-Lagrangian Column Kernels`): Column-serial sedimentation kernel pattern for falling-species routines is generic across bulk schemes; WSM6 `nislfv` is the demonstrated instance.
  - `Rule_25` (`## Rule 25: Phase Structure for Large Routine Ports`): Mandates a staged five-phase workflow for large routine ports so implementation and validation remain tractable.
  - `Rule_29` (`## Rule 29: The loops/dtcld Minor Timestep Outer Loop`): Minor-timestep subdivision wrapping the full process group sequence is generic across bulk schemes; WSM6 `dtcld` is the demonstrated instance.
  - `Rule_30` (`## Rule 30: Phase 3 Diagnostic Instrumentation Protocol`): Defines canonical per-group diagnostic tags and gated target-column prints to create a Fortran-versus-C++ retreat oracle.
  - `Rule_30A` (`## Rule 30 Addendum: Runtime High-Precision Diagnostic Instrumentation`): Upgrades diagnostics to runtime-controlled high-precision output with synchronized Path A/Path B print contracts.
  - `Rule_30D` (`## Rule 30D: Tier1.5 Block-Signature Diagnostics`): Defines scoped block-signature probes between compact Tier1 tags and Tier2 line traces for axis classification before first-line escalation.
  - `Rule_31` (`## Rule 31: Runtime Fortran/C++ toggle — Morrison pattern`): Requires a runtime implementation-path toggle via ParmParse instead of brittle build-time switching.
  - `Rule_33` (`## Rule 33: Pre-compile reflexive alignment pass — group-by-group`): Requires a pre-compile side-by-side scan for dropped boundary statements and off-by-one loop-bound translation errors.
  - `Rule_35` (`## Rule 35: Plotfile Parity Lane — Campaign Structure`): Defines ordered milestone-based plotfile parity campaigns, serial-first scope, and first-fail narrowing before retreat.
  - `Rule_36` (`## Rule 36: Restart Reproducibility Check`): Requires one-step shared-source restart reproduction at first-fail step to separate physics divergence from restart-state issues.
  - `Rule_37` (`## Rule 37: Plotfile Divergence to Tag Retreat Linkage`): Maps first-fail plotfile metadata back into canonical tag retreat with escalating debug levels and value-level evidence.
  - `Rule_38` (`## Rule 38: Plotfile Campaign TSV Schema Extension`): Defines structured run-ledger schema fields and fail-artifact requirements for plotfile parity campaigns.

Rules currently in `wsm6_implementation_notes.md` that are generic and should move:
  - `ImplNote_Rule9_ModuleSave_New` (`## Rule 9: Module save Variables → Class Member Init [NEW - WSM6]`): Initial hypothesis for translating Fortran save-variable initialization into C++ class-member initialization lifecycle.
  - `ImplNote_Rule9_ModuleSave_Complete` (`## Rule 9: Module save Variables → Class Member Init [COMPLETE]`): Finalizes the method for resolving save-variable init inputs and one-time initialization capture into kernels.
  - `ImplNote_AppendixA_ScopeGate` (`### Scope gate`): States serial-first campaign scope and defers MPI/GPU parity lanes until serial acceptance is complete.
  - `ImplNote_AppendixA_FcompareWorkflow` (`### Per-milestone fcompare workflow`): Defines first-fail fcompare usage including `-z` and `-d` evidence capture plus diff artifact archival.
  - `ImplNote_AppendixA_ArchiveDiffOutput` (`# Archive default diff output immediately`): Requires immediate rename/archive of default diff output to preserve first-fail provenance.
  - `ImplNote_AppendixA_RestartRepro` (`### Restart reproducibility (Rule 36)`): Applies shared-source one-step restart reproducibility protocol before classifying divergence for retreat.
  - `ImplNote_AppendixA_BoundedRefinement` (`### Bounded per-step refinement (coarse milestone failure)`): Defines bounded per-step restart refinement to narrow coarse-cadence failures to earliest failing substep `N*`.

Operator-policy additions now captured in `wsm6_validation_operator.md`:
  - `Rule36_RestartStateMatrixPolicy` (`## Rule36 Restart-State Matrix Policy`): Formalizes the four-leg restart-state matrix (`Fpath_from_Fchk`, `Cpath_from_Fchk`, `Fpath_from_Cchk`, `Cpath_from_Cchk`) after one-step mismatch and blocks Rule37 until restart/source-state provenance is classified.
  - `VariableAwareGrowthAdjudication` (`## Variable-Aware Growth Adjudication`): Requires read-only growth adjudication over existing fcompare artifacts so epsilon-scale pressure/thermo residuals are not over-classified as material failures; includes diff-plotfile + fextract checks for quantized binary-quantum roundoff lattices.

## 2. The name-layer crossover
Many sections sit at a three-way crossover: generic physics concept, conventional name, and WSM6-specific implementation detail. Future porters should inherit concepts and workflow, reuse conventional names where they match, and re-derive only the truly scheme-specific layer.

### 2a. Generic physics concepts
These are `GENERIC_MICROPHYSICS_CONCEPT` rows and represent process patterns a bulk microphysics porter should inherit.

  - `Rule_24`: Column-serial sedimentation kernel pattern for falling-species routines is generic across bulk schemes; WSM6 `nislfv` is the demonstrated instance.
  - `Rule_29`: Minor-timestep subdivision wrapping the full process group sequence is generic across bulk schemes; WSM6 `dtcld` is the demonstrated instance.

### 2b. Conventional names (shared but not universal)
The rows below are `COMMON_MICROPHYSICS_CONVENTION` and document name-layer artifacts that sit above generic concepts.

Audit rows in this tier:
  - `ImplNote_MicVarWSM6_EnumDraft`
  - `ImplNote_WSM6Ind_EnumDraft`
  - `ImplNote_Rule6_SpeciesIndexedArrays`
  - `ImplNote_Rule7_1DArrays`
  - `ImplNote_WSM6Ind_EnumUpdated`
  - `ImplNote_ERFWSM6H_ClassMembers`
  - `ImplNote_AppendixA_VarTagMapping`

Concept/name crossover table:

```
concept | WSM6_name | Morrison_name | Thompson_name | notes
water vapor mixing ratio | qv (also q) | unknown in current corpus | unknown in current corpus | concept is generic; symbol likely conventional but not audited here
cloud water mixing ratio | qc | unknown in current corpus | unknown in current corpus | common bulk variable naming layer
cloud ice mixing ratio | qi | unknown in current corpus | unknown in current corpus | common bulk variable naming layer
rain mixing ratio | qr | unknown in current corpus | unknown in current corpus | common bulk variable naming layer
snow mixing ratio | qs | unknown in current corpus | unknown in current corpus | common bulk variable naming layer
graupel mixing ratio | qg | unknown in current corpus | unknown in current corpus | may be scheme-dependent; see Open Question 2
autoconversion cloud->rain | praut | unknown in current corpus | unknown in current corpus | generic concept; symbol may differ by scheme
accretion cloud by rain | pracw | unknown in current corpus | unknown in current corpus | generic concept; symbol may differ by scheme
rain evaporation | prevp | unknown in current corpus | unknown in current corpus | generic concept; symbol may differ by scheme
substep timestep variable | dtcld | unknown in current corpus | unknown in current corpus | concept is generic; name is convention-layer
```

### 2c. Truly WSM6-specific
These are `SCHEME_SPECIFIC_FACT` rows and must be re-derived for other schemes.

  - `Rule_27`: Canonical WSM6 execution group ordering and dependencies for incremental implementation and retreat.
  - `Rule_28`: WSM6 working-array allocation layout across `fab_box`, `box2d`, and column-local stack arrays.
  - `ImplNote_Rule5_StructuralDiff`: WSM6-specific structural deltas relative to Morrison baseline.
  - `ImplNote_Rule10_RuntimeVariantFlag`: WSM6 hail-versus-graupel runtime branch and exact `hail_opt`-controlled variable set.
  - `Caveat_WSM6_SG_DefaultGraupelOnly_HailDeferred`: The NISLFV_SG closure evidence is scoped to default graupel mode; hail-mode coefficient routing is deferred to a separate validation lane.

Lower-confidence flag:
  - `Rule_28` could partially be convention-layer because the FAB allocation pattern itself is reusable, but the concrete field layout remains WSM6-specific.

## 3. Accumulated debugging memory
These rows encode knowledge learned from real failures and should be treated as high-loss if removed. They are not obvious from first principles, and future ports would re-incur these failures.

  - `Rule_13`: Prevents slab allocation/indexing mismatch for 2D accumulators (`box2d` vs full `fab_box`) that can trigger parity drift or memory faults.
  - `Rule_21`: Prevents silent numerical drift from replacing exact special-function behavior with non-equivalent approximations at initialization time.
  - `Rule_21_Addendum_ExactHelperSemantics`: Extends exact-helper parity requirements into runtime local helper paths and records non-causal parity cleanups separately from causal fixes.
  - `Rule_24_Addendum_NISLFVBoundary`: Captures the high-loss boundary-proof method for sedimentation kernels: bracket working temporaries (for example `denqrs1`) before patching downstream update logic.
  - `RegCase_WSM6_SG_RULE37_NISLFV_SG_STEP506_TO_600`: Captures the SG mixed-phase closure chain (`patch_commit=4264cfb4`) from Tier1.5 boundary proof through Rule36 restart-state matrix and variable-aware growth adjudication.
  - `Evidence_WSM6_SG_PostFix_Pressure_ULP_Lattice`: Captures the postpatch ULP-lattice pattern (`pressure_quantum=2^-36`) used to adjudicate step-511 pressure residuals as non-material.
  - `Rule_32A`: Prevents uncontrolled forward implementation after first divergence by enforcing stop-diff-retreat and stepwise narrowing gates.
  - `Rule_34`: Prevents semantic mistranslation of Fortran control flow (`goto`, loop targets, `cycle`/`exit` mapping) that can pass smoke tests but fail later.

`Rule_32A` and `Rule_34` were correctly reclassified from scheme-framed interpretations to `ACCUMULATED_DEBUGGING_MEMORY` because their failure modes recur across any scheme port.

## 4. Test-case-specific behavior
The following rows are specific to SquallLine_2D and should not be treated as WSM6 scheme facts. A different WSM6 testcase would keep the same validation concept but use different numbers, schedules, and artifacts.

  - `ImplNote_AppendixA_PlotfileParity`: Generic concept is campaign orchestration; testcase-specific part is this SquallLine_2D campaign framing.
  - `ImplNote_AppendixA_ExecutablePolicy`: Generic concept is build/run policy; testcase-specific part is exact executable/run-location policy for this case.
  - `ImplNote_AppendixA_DtConvention`: Generic concept is fixed-time-step reproducibility; testcase-specific part is `fixed_dt=1.0` for this family.
  - `ImplNote_AppendixA_MilestoneTable`: Generic concept is milestone checkpoints; testcase-specific part is the exact steps/times (A-I, 300/500/600/3000/6000/9000/12000).
  - `ImplNote_AppendixA_RunStructure`: Generic concept is short/long pair decomposition; testcase-specific part is this case’s cadence settings.
  - `ImplNote_AppendixA_CommandTemplates`: Generic concept is paired command templates; testcase-specific part is concrete command lines and paths.
  - `ImplNote_AppendixA_SoftBaselines`: Generic concept is soft epsilon reference; testcase-specific part is numeric values from a specific WSM6 run.

Milestone step numbers, onset timing, and campaign artifact paths should live in `regression_cases.tsv`-linked records, not in generic entry rules.

Post-SG-fix testcase-specific evidence now recorded:
  - `Evidence_WSM6_SG_PostFix_Pressure_ULP_Lattice`: ULP-scan artifacts under `rule36_refine_plt00600_from_chk00500_post_sg_fix_clean_6782b1f6_20260506T094302Z/artifacts/ulpscan` show step-511 isolated pressure delta (`4.3655745685100555e-11`) and steps 558/600 quantization to `1.4551915228366852e-11` (`2^-36`), supporting variable-aware `EPSILON_OK` adjudication.

## 5. Short entry skill index
The lists below are lookup indices only.

### Skill 1 — Wire a Fortran physics scheme into ERF
Question: Can ERF call the original Fortran scheme as a reference path?
Phase: `phase_1`
Path: `path_a`
Rules:
  - `Rule_2`, `Rule_3`, `Rule_4`, `Rule_10`, `Rule_12`, `Rule_13`, `Rule_14`, `Rule_25`, `Rule_27`, `Rule_30`, `Rule_30A`, `Rule_31`, `Rule_32`, `Rule_32A`, `Rule_35`, `Rule_36`, `ImplNote_5b_OptionalArgs`, `ImplNote_5c_ErrmsgErrflg`, `ImplNote_PendingBeforeComplete`, `ImplNote_FilesToCreate`, `ImplNote_AppendixA_PlotfileParity`, `ImplNote_AppendixA_ExecutablePolicy`, `ImplNote_AppendixA_CommandTemplates`.

Wrong-file blockers before clean citation:
  - `Rule_27` should move scheme-specific content to `wsm6_implementation_notes` with a pointer from the generic rule.
  - `ImplNote_5b_OptionalArgs` and `ImplNote_5c_ErrmsgErrflg` are convention-level and should migrate to the generic rulebook/operator corpus.

### Skill 2 — Port a Fortran physics scheme to native AMReX C++
Question: Can I produce a native C++ implementation that mirrors the Fortran reference?
Phase: `phase_2` through `phase_4`
Path: `path_b`
Rules:
  - `Rule_1`, `Rule_2`, `Rule_4`, `Rule_9`, `Rule_10`, `Rule_11`, `Rule_12`, `Rule_13`, `Rule_14`, `Rule_21`, `Rule_22`, `Rule_23`, `Rule_24`, `Rule_25`, `Rule_27`, `Rule_28`, `Rule_29`, `Rule_30`, `Rule_30A`, `Rule_31`, `Rule_32`, `Rule_32A`, `Rule_33`, `Rule_34`, `ImplNote_Rule5_StructuralDiff`, `ImplNote_5a_2DIndexing`, `ImplNote_5d_VrecVsqrt`, `ImplNote_MicVarWSM6_EnumDraft`, `ImplNote_WSM6Ind_EnumDraft`, `ImplNote_PendingBeforeComplete`, `ImplNote_Rule6_SpeciesIndexedArrays`, `ImplNote_Rule7_1DArrays`, `ImplNote_Rule8_IntLogicalArrays`, `ImplNote_WSM6Ind_EnumUpdated`, `ImplNote_Rule9_ModuleSave_New`, `ImplNote_Rule9_ModuleSave_Complete`, `ImplNote_Rule10_RuntimeVariantFlag`, `ImplNote_ERFWSM6H_ClassMembers`, `ImplNote_AppendixA_VarTagMapping`.

Wrong-file blockers before clean citation:
  - `Rule_27`, `Rule_28` are currently in skill doc but carry heavy WSM6-specific payload.
  - `ImplNote_5a_2DIndexing`, `ImplNote_5d_VrecVsqrt`, `ImplNote_Rule8_IntLogicalArrays`, `ImplNote_Rule9_ModuleSave_New`, and `ImplNote_Rule9_ModuleSave_Complete` are generic enough to migrate into rulebook/operator-facing docs.

### Skill 3 — Validate and retreat a native port against the reference
Question: Where does native C++ first diverge from Fortran, and what evidence is sufficient to change state?
Phase: `phase_3` through `phase_5`
Path: `path_a vs path_b`
Rules:
  - `Section_Status_PARTIAL`, `Rule_12`, `Rule_30`, `Rule_30A`, `Rule_32`, `Rule_32A`, `Rule_35`, `Rule_36`, `Rule_37`, `Rule_38`, `ImplNote_AppendixA_PlotfileParity`, `ImplNote_AppendixA_ScopeGate`, `ImplNote_AppendixA_FcompareBinary`, `ImplNote_AppendixA_ExecutablePolicy`, `ImplNote_AppendixA_DtConvention`, `ImplNote_AppendixA_MilestoneTable`, `ImplNote_AppendixA_RunStructure`, `ImplNote_AppendixA_CommandTemplates`, `ImplNote_AppendixA_FcompareWorkflow`, `ImplNote_AppendixA_ArchiveDiffOutput`, `ImplNote_AppendixA_RestartRepro`, `ImplNote_AppendixA_BoundedRefinement`, `ImplNote_AppendixA_SoftBaselines`, `ImplNote_AppendixA_VarTagMapping`.

Wrong-file blockers before clean citation:
  - `ImplNote_AppendixA_ScopeGate`, `ImplNote_AppendixA_FcompareBinary`, `ImplNote_AppendixA_FcompareWorkflow`, `ImplNote_AppendixA_ArchiveDiffOutput`, `ImplNote_AppendixA_RestartRepro`, and `ImplNote_AppendixA_BoundedRefinement` are generic validation-operator content currently embedded in scheme notes.

## 6. Migration targets
The following moves are concrete relocation actions, not abstract suggestions.

### 6a. From wsm6_implementation_notes → rulebook
  - Source: `ImplNote_5a_2DIndexing`. Reason: AMReX mapping convention is reusable across schemes. Pointer to keep: one line in WSM6 notes stating “see rulebook indexing convention; WSM6 has no exception.”
  - Source: `ImplNote_5b_OptionalArgs`. Reason: bridge optional-argument handling is cross-scheme interop policy. Pointer to keep: WSM6-specific optional list only.
  - Source: `ImplNote_5c_ErrmsgErrflg`. Reason: error signaling translation policy is general. Pointer to keep: WSM6-specific error fields only.
  - Source: `ImplNote_5d_VrecVsqrt`. Reason: Fortran vector-helper substitution pattern is generic to GPU ports. Pointer to keep: WSM6 call-site examples.
  - Source: `ImplNote_Rule8_IntLogicalArrays`. Reason: kernel-local integer/logical handling is an AMReX convention. Pointer to keep: WSM6 variable names that instantiate it.
  - Source: `ImplNote_AppendixA_FcompareBinary`. Reason: repository-path convention belongs with operator/rulebook execution instructions. Pointer to keep: WSM6 appendix can reference canonical location.

### 6b. From wsm6_implementation_notes → validation_operator
  - Source: `ImplNote_AppendixA_ScopeGate`. Reason: serial-first scope gate is lane policy, not WSM6 content. Pointer to keep: WSM6 appendix reference to operator section.
  - Source: `ImplNote_AppendixA_FcompareWorkflow`. Reason: first-fail `fcompare` evidence flow is generic operator procedure. Pointer to keep: WSM6-specific example commands only.
  - Source: `ImplNote_AppendixA_ArchiveDiffOutput`. Reason: diff artifact archival policy is generic provenance discipline. Pointer to keep: WSM6 naming example.
  - Source: `ImplNote_AppendixA_RestartRepro`. Reason: shared-source restart reproducibility is generic Rule 36 operator logic. Pointer to keep: WSM6 restart path format example.
  - Source: `ImplNote_AppendixA_BoundedRefinement`. Reason: coarse-window to `N*` narrowing is generic Rule 35/36 bridge logic. Pointer to keep: WSM6 bounded window command example.

### 6c. From skill doc → wsm6_implementation_notes
  - Source: `Rule_30`. Move: WSM6 canonical tag-name list (`DENFAC`, `NISLFV_R`, etc.). Replacement: keep generic instrumentation protocol in `Rule_30` with pointer to WSM6 tag table in implementation notes.
  - Source: `Rule_30A`. Move: WSM6-specific print-field ordering contract details. Replacement: keep generic print precision/runtime-level contract in `Rule_30A` with pointer.
  - Source: `Rule_32A`. Move: WSM6-specific epsilon thresholds and dry-benchmark numbers. Replacement: keep generic stop-diff-retreat doctrine in `Rule_32A` with pointer to testcase baselines.
  - Source: `Rule_35`. Move: WSM6 milestone step table. Replacement: keep generic campaign structure and first-fail narrowing in `Rule_35` with pointer.
  - Source: `Rule_37`. Move: WSM6 var-to-tag mapping table. Replacement: keep generic divergence-to-retreat linkage in `Rule_37` with pointer to scheme table.

Enum placement check:
  - `MicVar_WSM6` and `WSM6Ind` drafts are correctly placed in implementation notes.
  - Fields `qv`, `qc`, `qi`, `qr`, `qs`, and likely `qg` are convention-layer names (`COMMON_MICROPHYSICS_CONVENTION`) unless scheme evidence proves otherwise.
  - Fields such as `denfac` and WSM6Ind layout-specific derived intermediates are scheme-layer (`SCHEME_SPECIFIC_FACT`) until mapped to cross-scheme equivalents.

### 6d. From skill doc or implementation notes → regression_cases.tsv
  - Regression case: `WSM6_G1b_DENFAC_EVOLVED_STEP254`. Required artifact: shared restart source `plotfile_parity_v1d/long_fortran/chk00253`. Pass criterion: `G1b/DENFAC EPSILON_OK` at step 254 with structural count parity and clean-SHA evidence. Status today: row exists and is structurally complete but still `OPEN`.
  - Regression case: `WSM6_SG_RULE37_NISLFV_SG_STEP506_TO_600`. Required artifacts: shared `F505` Rule37 Tier1.5 boundary run and postpatch stepwise `501..600` fcompare growth artifacts. Pass criterion: boundary exact at `k_dbg=21..24`, shared-source step-506 mixed-phase signal removed, step-508 pressure failure removed, and `501..600` variable-aware acceptance (`STEP511_EPSILON_OK_CONTINUE`).
  - Regression case family to add: milestone parity gates A-I as explicit cases. Required artifacts: per-milestone paired plotfiles plus first-fail `-z/-d` artifacts when failing. Pass criterion: `PLOTFILE_AGREE` or within agreed epsilon for each milestone.
  - Regression case family to add: bounded refinement `N*` cases for coarse cadence failures. Required artifacts: shared-checkpoint restart pair and per-step `fcompare` series from `K+1..N`. Pass criterion: earliest fail `N*` identified and reproducible under Rule 36.

## 7. Recommended next actions
The list below is ordered by agentic reusability impact first, then knowledge-loss risk, then cleanup.

1. Extract the WSM6 tag-name list from `Rule_30` into a named table in `wsm6_implementation_notes.md` and leave a pointer in `Rule_30`. Files: `fortran_to_cpp_microphysics_skill.md`, `wsm6_implementation_notes.md`. Commit message: "Rule30: move WSM6 canonical tag table to implementation notes; keep generic instrumentation protocol with pointer." This is first priority because `Rule_30` is the load-bearing instrumentation oracle in Phase 3/4, and cross-scheme readers currently cannot separate generic procedure from WSM6 content.

2. Extract the var-to-tag mapping from `Rule_37` into a named implementation-notes section and replace it with a pointer in `Rule_37`. Files: `fortran_to_cpp_microphysics_skill.md`, `wsm6_implementation_notes.md`. Commit message: "Rule37: move WSM6 var-to-tag map to implementation notes; retain generic divergence-to-retreat linkage." This is second because it is the bridge from plotfile failure to tag retreat and is essential for reproducible diagnostics.

3. Move Appendix A generic validation-policy subsections (`ScopeGate`, `FcompareWorkflow`, `ArchiveDiffOutput`, `RestartRepro`, `BoundedRefinement`) into `wsm6_validation_operator.md`, then leave one-line pointers in implementation notes. Files: `wsm6_implementation_notes.md`, `wsm6_validation_operator.md`. Commit message: "Operator: absorb generic plotfile policy from WSM6 Appendix A; leave scheme-note pointers." This is third because operator decisions should not depend on scheme-notes discovery.

4. Move generic `CODEBASE_CONVENTION` subsections (`5a`-`5d`, `Rule8`, `FcompareBinary`) into the generic rule corpus. Files: `wsm6_implementation_notes.md`, `fortran_to_cpp_microphysics_skill.md` (or operator doc for binary path convention). Commit message: "Rulebook: migrate generic AMReX/interop conventions out of WSM6 notes." This is fourth because it improves reuse and reduces duplicate rediscovery across future ports.

5. Write the three short entry skills as index files that reference rule IDs only, using Section 5 of this summary as source of truth. Files: new skill index files under docs/skill area (project choice), with no rule text duplication. Commit message: "Add three entry-skill rule indices (bridge/native/validation) sourced from rules_audit.tsv." This is fifth because it unlocks agent onboarding speed with minimal maintenance burden.

6. Add a named `COMMON_MICROPHYSICS_CONVENTION` subsection in `wsm6_implementation_notes.md` with columns `concept | WSM6_name | Morrison_name | Thompson_name | notes`; fully populate WSM6 and mark unknowns explicitly for others. Files: `wsm6_implementation_notes.md`. Commit message: "Add concept/name crossover table for common microphysics conventions (WSM6 filled, others TBD)." This is sixth because it creates the reusable concept-to-symbol bridge required for future scheme ports.

7. Mark superseded sections in the skill doc with explicit `SUPERSEDED` headers and pointers to replacements. Files: `fortran_to_cpp_microphysics_skill.md`. Commit message: "Skill doc hygiene: mark superseded rules explicitly with replacement pointers." This is seventh because it reduces citation ambiguity and prevents stale rule use by agents.

8. Add a `TODO` note to `Rule_14` indicating required future subdivision, without splitting now. Files: `fortran_to_cpp_microphysics_skill.md`. Commit message: "Rule14 backlog: add explicit TODO for planned subdivision." This is eighth because it is important but not currently blocking reusable operations.

## 8. Open questions
The following were not resolvable from the audit corpus alone and need human decision or additional source data.

1. Morrison and Thompson process-name equivalents for the Section 2b convention table remain unknown in this corpus. Resolution input needed: authoritative symbol/process inventories from active Morrison and Thompson ERF ports.

2. Graupel as a separate prognostic species (`qg`) may be common but was not proven universal in the audited files. Resolution input needed: cross-scheme species inventory comparison to decide whether `qg` should remain convention-layer or become conditional scheme-specific guidance.

3. `Rule_14` remains broad and was flagged for subdivision. Resolution input needed: a design decision on whether Init-signature patterns are generally reusable at rulebook level or tightly coupled to scheme-specific member composition.

4. `Rule 5` remains marked `PARTIAL` with unresolved optional-argument and `errmsg/errflg` policy questions. Resolution input needed: ERF-wide interop convention decision and confirmation whether these are Phase-5 blockers or post-parity cleanup items.

5. `SCHEME_SPECIFIC_FACT` count is lower than initially expected. Resolution input needed: targeted per-field review of `MicVar_WSM6` and `WSM6Ind` entries to confirm whether some rows currently in `COMMON_MICROPHYSICS_CONVENTION` should be recategorized as scheme-specific.
