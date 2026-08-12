Follow Exec/.DevTests/fortran-to-cpp-microphysics/SKILL.md exactly.

Target: close WDM6 Bubble Milestone A on the `plotfile_parity_serial` lane.
Milestone A currently FAILS. Milestones B and C are deliberately not attempted.

========================
READ FIRST (skill + sub-operators)
========================

Skill workflow order, all six steps. The extracted operators are compact BY
DESIGN and are the intended entry point — do not reach for the source corpus by
reflex. But follow SKILL.md's own routing: the prior session skipped step 1 and
missed `core-rules.md:37`, which names the exact bug class behind the heap
corruption fixed in `38a3f6ead`.

1. `references/core-rules.md`
2. `references/applications/wdm6.md`
3. `references/validation-operators.md`   <- sub-operator: lanes, ledgers, retreat
4. `references/process-operators.md`      <- sub-operator: campaign habits, provenance
5. `references/coverage-matrix.md`
6. `references/source-materials.md` / `references/source-ledger.md`

The two `*-operators.md` files are the skill's sub-operators. `validation-operators.md`
governs this task (Path/Lane model, Stop-Diff-Retreat, restart rules, ledger model,
evidence discipline). `process-operators.md` governs provenance repair and debug
discipline. Neither is optional here.

Then the WDM6/Bubble-specific set:

- `references/applications/wdm6/bubble_campaign.md` — the campaign definition
  (milestone table, run posture, fcompare binaries, variable-to-tag map,
  "sites outside every group tag")
- `references/applications/wdm6/group_map.md` — canonical file-order groups
- `references/applications/wdm6/campaign_wdm6_bubble_v2_report.md` — the Milestone A
  evidence chain
- `references/applications/wdm6/campaign_{executions,decisions}.tsv` — Rule 38
  lane ledgers

When to reach past the compact operators, and only then:

- Rule 35 CAMPAIGN PROCEDURE is the one operative gap. `validation-operators.md`
  establishes that milestone parity is a campaign with stop/go gates, but the
  milestone-ladder structure, the fcompare protocol (zero tolerance first,
  `-z`/`-d` archive on first fail), and the per-step narrowing loop live only in
  `references/source-corpus/fortran_to_cpp_microphysics_skill.md` (Rules 35-38).
  For Bubble specifically, `bubble_campaign.md` now carries all of it — prefer it.
- The compact operators DO carry the Rule 36 -> Rule 37 ordering ("use
  shared-source one-step restart reproduction to classify first-fail plotfile
  divergence"; "use the restart-state matrix before invoking later
  plotfile-to-tag linkage"). Read them closely rather than substituting the corpus.
- `references/source-corpus/wsm6_validation_operator.md` adds the 9-axis
  differential diagnosis, Tier 1.5 block-signature policy, and variable-aware
  growth adjudication. Use for a hard divergence, not for routine progression.

========================
HAZARD: STALE SKILL MIRROR
========================

`.claude/skills/fortran-to-cpp-microphysics/` is an UNTRACKED, PARTIAL mirror.
Its `references/applications/` holds only `wdm6.md` and `wsm6.md` — it has no
`wdm6/` subdirectory at all, so `group_map.md`, `bubble_campaign.md`, the
campaign report, and the lane ledgers are INVISIBLE when the skill loads from
there. Invoking the skill reports that directory as its base.

Treat `Exec/.DevTests/fortran-to-cpp-microphysics/` as the source of truth and
open the WDM6 files by explicit path. Consider syncing or symlinking the mirror.

========================
STATE AS OF 2026-08-12
========================

Commits on `wdm6-mpi-cpu` (unpushed at time of writing):
- `38a3f6ead` WDM6: fix bridge heap corruption in 2D precip buffer allocation
- `dfcc24231` WDM6: add serial plotfile-parity campaign lane for Bubble

Resolved: the Fortran bridge abort recorded open since G5a
("free(): invalid pointer downstream of G5a"). It was NOT a diagnostic-buffering
artifact — it reproduced at `microphysics_debug=0`, `GFORTRAN_UNBUFFERED_ALL=y`,
no `mpirun`, as a glibc `_int_malloc` assertion. valgrind attributed 28 invalid
reads + 15 invalid writes to `mp_wdm6_run_c`, all past 6400-byte (200x4
active-tile) FArrayBox blocks. The seven 2D precip buffers were allocated on the
ACTIVE tile box while the bridge declares them `dimension(ims:ime,jms:jme)` on
STORAGE bounds. Fixed per `ERF_AdvanceWSM6.cpp:963-994`. Post-fix valgrind:
0 errors from 0 contexts. The bridge now completes step 1 and writes `plt00001`.

Open: Milestone A FAIL at step 1, zero tolerance. `plt00000` pair AGREE, so the
divergence is entirely inside the first WDM6 call. `temp` rel 2.291e-2 (6.63 K),
`eq_pot_temp` rel 1.986e-2, `rhotheta` rel 8.859e-3, `pressure` rel 8.722e-3,
`qi` rel 2.647e-3, `qsnow`/`nc`/`nr` rel ~6.2e-4. `density`, `qc`, `qrain`,
`qgraup` and all velocities bitwise equal; `qv` rel 2.667e-11.
Zone `(i,j,k) = (199,3,90)`. Rule 36 `restart_repro_status=MATCH`.

========================
ORDERED NEXT ACTIONS
========================

A) Clean-SHA re-run of Milestone A.
   All existing lane rows are `dirty_flag=dirty`, triage-only, because the fix
   was uncommitted when they were produced. The fix is now committed. Re-run and
   append clean rows before any promotion. Check the embedded runtime hash first
   (`process-operators.md` provenance rule); if it disagrees with the worktree
   SHA, repair with the tiered build-info refresh — note tier 1 (object + dep)
   was INSUFFICIENT last time because the generated source carried the stale
   hash; tier 2 (remove `build_wdm6/erf_srclib/AMReX_buildInfo.cpp`) was required
   and triggers a CMake regeneration.

B) Restore the WSM6 all-k print contract in WDM6 group tags. PREREQUISITE for
   trusting the frontier.
   WDM6 group tags index every field at `kts` with no `do k` loop, e.g.
   `'WDM6-FORT_PRE_G7', kts, qrs(i_dbg_local,kts,1), ...`. Each tag therefore
   validates ONE cell. WSM6 on branch `wsm6-mpi-cpu` loops the full column:
   `if (mpdbg_level >= 1 .and. loop == 1) then / do k = kts, kte / ...`.
   Consequence: all 33 P0..G17 closures rest on ~1/100th of the evidence their
   WSM6 counterparts did, and NO group evidence exists at k=90 where Milestone A
   fails. A Rule 37 retreat at the failing column (199,3) returned
   `max_abs_diff=0.000e+00` for every group — inconclusive by construction.
   Also: `G10A/G10B/G10C` emit ~20000 Fortran lines (full i-k plane) against one
   C++ line per group, so those have no matched surface either.
   Do not reopen or promote any group until this is restored and re-run clean.

   THE HARNESS NEEDS THE SAME FIX. The tag-frontier harness is `tools/build_ci`
   (`./tools/build_ci <group> [--no-build] [--all-tags]`), untracked, and it is
   what every group closure in the ledgers was produced with
   (`command=./tools/build_ci G13e`). Two limitations to lift together with the
   print contract, or restored all-k prints will still report single-cell
   agreement:
   - it hardwires `erf.micro_diag_target_column="100 0"`; parameterize it, since
     Milestone A fails at (199,3)
   - its comparison takes only the FIRST matching marker line per tag
     (`grep ... | head -1`, and the same in its embedded python), so it cannot
     see per-k lines even when they are emitted. Match markers by the k field
     and reduce over the column.
   `tools/` is a throwaway local area holding nothing tracked, so treat build_ci
   results as triage-only for the same reason as the fcompare binaries.

C) Resolve the `theta` writeback design question.
   CANDIDATE site, NOT confirmed root cause. The Fortran branch converts updated
   absolute temperature back to potential temperature after the bridge call
   (`theta_arr(i,j,k) = t_arr(i,j,k) / exner`, `ERF_AdvanceWDM6.cpp` Fortran
   branch only). The native branch never writes `mic_fab_vars[theta]`; it ends at
   precipitation accumulation. `Copy_Micro_to_State` then forms
   `RhoTheta = rho * theta` from that field, so the native path appears to carry a
   pre-microphysics theta and drop microphysical latent heating. That matches the
   signature (moisture agrees, thermodynamics diverges at percent scale).
   Two open points before patching:
   - retreat discipline requires paired `PRE_<site>`/`POST_<site>` first-line
     evidence at the scoped site; structural absence is not proof
   - the design choice is deliberate: mirror the bridge conversion in the native
     branch, or have the native kernel update `theta` directly
   Residual `qi`/`qsnow`/`nc`/`nr` differences are NOT explained by writeback
   alone. Do not treat C) as settled by itself.

D) Fix `ni` plotfile emission.
   With `erf.plot_vars_1` overridden to append `nc ni nr`, only 31 variables
   appear: `nc` and `nr` are present, `ni` is not. `RhoQ8`/`nn` (CCN) is therefore
   uncompared. Confirmed an ERF emission gap, not a tooling artifact — both
   fcompare binaries agree to every printed digit. ERF gates `ni` on
   `n_qstate_moist >= 8` at `Source/IO/ERF_Plotfile.cpp:1263-1282`; `nr` requires
   `>= 9` and does appear, so the gate is not the obvious cause.
   Naming caveat: WDM6 `RhoQ8` is CCN/aerosol `nn` but ERF's plot label for that
   slot is `ni`. Do not rename; the mismatch is recorded in `bubble_campaign.md`.

E) Only after A-D: Milestones B (step 10) and C (step 100), in order, stopping at
   first failure. C is the direct analogue of the nightly `Bubble_WSM6` page.

========================
DEFERRED LANES
========================

- 4-rank MPI: needs an `ERF_ENABLE_MPI=ON` build (`build_wdm6` is MPI=OFF).
  Distinct `validation_lane`, own gates. Serial acceptance first.
- GPU: separate lane, also gated on serial acceptance.
- Short SquallLine_2D WDM6 campaign: the real cold-phase gate. Bubble does
  exercise cold-phase groups at upper levels (contrary to an earlier "warm-rain
  only" framing) but weakly. Use in-repo `inputs_squallline` (dt=0.25,
  plot_int_1=60, 192x4x128) with a `moisture_model` override; note
  `inputs_moisture_WSM6` is not on this branch. Keep it much shorter than WSM6's
  A-I ladder: WSM6 had all species active by t=600, and F/G/H at t=3000/6000/9000
  were paper comparisons, not coverage. Determine onset empirically.
- Blessing a `Bubble_WDM6` benchmark + nightly test: `Tests/CTestList.cmake` has
  no `Bubble_WSM6` entry, so the nightly lives in the external regression-suite
  config, not this repo.

========================
UNPLACED FILES
========================

Left uncommitted pending a decision on where they should live:
- `references/source-corpus/{decisions,executions}.tsv` (modified; G13e-G17 rows,
  several `formal_clean_rerun_required=true`)
- package-root `decisions.tsv` / `executions.tsv` (ad-hoc 8-column schema that
  cannot hold Rule 38 lane fields)
- `references/applications/wdm6/pre_post_variable_audit.md`
- `references/applications/wdm6/wsm6_group_mapping_comparison.md`
- `PHASE2_PARITY_FIX_LOG.txt`, `g11_post_handoff.md`, `skill_package_handoff.md`

Note the closure bookkeeping is split across two ledgers with different schemas
and neither is complete. Union of WDM6 rows: G4, G5a, G5b, G9, G10a, G10b, G10d,
G11, G12, G13a, G13b, G13c, G13e, G13f, G13g, G14, G15, G16b, G17. Absent from
both: P0, G1a, G1b, G1c, G2, G3, G5c, G6, G7, G8, G10c, G10e, G13d, G16a.

========================
RUN POSTURE (quick reference)
========================

Runner: `Exec/.DevTests/fortran-to-cpp-microphysics/references/applications/wdm6/campaign_bubble_wdm6 <campaign_id> <pair> <max_step> <plot_int>`
Serial, no `mpirun`, `microphysics_debug=0`, `GFORTRAN_UNBUFFERED_ALL=y` on the
Fortran leg only, unique `campaign_id` per execution (`campaign_*/` is gitignored).
Compare zero-tolerance first; never use tolerance flags to hide a divergence.
