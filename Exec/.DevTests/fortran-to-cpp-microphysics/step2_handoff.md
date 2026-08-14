Follow Exec/.DevTests/fortran-to-cpp-microphysics/SKILL.md exactly.

Target: close the WDM6 Bubble step-2 discontinuity on the `plotfile_parity_serial`
and `tag_frontier` lanes.

**Milestone A now PASSES bitwise (32/32 variables).** Milestone B fails. The
frontier is a step-2 discontinuity of a different character to anything closed
so far. Milestone C is not attempted.

Supersedes `milestone_a_handoff.md`, whose actions A, C and D are all closed.

========================
READ FIRST
========================

Skill workflow, all six steps, from `Exec/.DevTests/fortran-to-cpp-microphysics/`
by explicit path. The mirror at `.claude/skills/fortran-to-cpp-microphysics/`
still has no `references/applications/wdm6/` subdirectory, so the group map,
campaign definition, lane ledgers and manifest are invisible when the skill
auto-loads. This was true at the start of the last session and was not fixed.

New since the last handoff, read these too:

- `references/applications/wdm6/campaign_validation_manifest.tsv` — the compact
  frontier table. **Read this first**; it answers "where are we" in one line per
  milestone. It did not exist before and it has gone stale twice, so verify it
  against `campaign_decisions.tsv` if the two disagree.
- `references/applications/wdm6/group_map.md` — now carries a
  "Bridge-wrapper sites" section declaring the `W1..W4` namespace for sites
  outside the Fortran kernel span. These are NOT `G*` groups and must never join
  the `P0..G17` sequence.
- `Source/Microphysics/WDM6/ERF_WDM6.H` — the LITERAL PRECISION CONTRACT comment
  block. Load-bearing; read before touching any WDM6 constant.

========================
STATE AS OF 2026-08-14
========================

HEAD `0817610b7`, `26.06-318-g0817610b756b`, 2 commits unpushed.
Ledgers: 34 execution rows, 26 decision rows, 3 manifest rows.

Milestone status:

| | status |
| --- | --- |
| A (step 1) | **PASS bitwise, 32/32** |
| B (step 10) | FAIL, `rhoQ3` rel 1.385e+02 |
| C (step 100) | never run |

Tag frontier at step 1: no material divergence in any of 33 groups. The only
nonzero readings are four instances of `2.220e-16` on `denfac` at k=9 in G8,
G11, G13b and G13c — one ULP, consistent across four independent groups,
quantization not defect.

**The frontier**: `plt00001` AGREES bitwise on all 32 variables while `plt00002`
FAILS with 24 nonzero and `rhoQ3` at 148.4 relative. Bitwise-identical to
catastrophic in one step is a discrete event, not amplification, so a
step-dependent code path is implicated rather than a constant or an expression.

========================
WHAT CLOSED THIS SESSION
========================

Six defects, in the order found. All but the last were outside every group tag.

1. `W1_THETAWB` — the native path never wrote `mic_fab_vars[theta]`, so
   `Copy_Micro_to_State` formed `RhoTheta` from a pre-microphysics field.
   Confirmed by paired PRE/POST at the site: native PRE == POST at 100/100
   levels, divergence 6.625510053 at k=90 matching the Milestone A `theta`
   error to ten significant figures. Fixed by mirroring the bridge conversion.
2. Adam's legacy Steps 10 and 11 — untagged pre-group blocks that double-
   sedimented `qs`/`qg`/`qi` after `POST_G17` and applied `nc`/`nr` floors with
   no Fortran counterpart. `#if 0`'d. A controlled A/B proved they were
   pre-existing and independent of the theta fix.
3. Plotfile column rotation — `derived_names` in `ERF.H` declares
   `qt,qn,qp,qsat,nc,ni,nr,ns,ng` while the writer emitted the number
   concentrations first, so four columns reported each other's data. Shipped
   with PR #3218 (30 May 2026) and **also affects Morrison**. Fixed by moving
   the writer blocks to match the declared contract.
4. `nn` had no plot variable at all — WDM6's CCN slot was uncomparable. Added as
   a first-class variable with a `ccn_number` capability. Note `ni` means cloud
   ice number ERF-wide (`ERF_Morrison.H:47`) and WDM6 has none; `ERF_DataStruct.H`
   already set `ni = -1` for WDM6, so the old campaign note claiming ERF labels
   that slot `ni` was wrong and is retracted.
5. `dimax` literal precision — G13e `pidep`. See below.
6. `roqimax` literal precision — G13f `psaut`. See below.

========================
THE LITERAL PRECISION CONTRACT
========================

This is the most reusable finding and it is not finished.

Every parameter in `ERF_module_mp_wdm6.F90` is declared
`real(kind=kind_phys), parameter :: dimax = 500.e-6` with **no kind suffix**, so
Fortran evaluates the literal in single precision and only then widens it. All
22 are written this way; none is suffixed. `dimax` is therefore
`+5.00000023748725652695E-04`, bit-for-bit `float32(5e-4)` widened, against a
true double in C++.

`ERF_WDM6_F32_LITERALS` (CMake, **defaults ON**) routes all 39 constants in
`ERF_WDM6.H` through `wdm6_literal()`, which reproduces that rounding. Verified:
the Tier 2 `pidep` decomposition goes 14/14 bitwise equal with **no Fortran
change**. Turning it OFF without also setting `ERF_WDM6_FORTRAN_REAL8` puts the
paths back out of agreement — four constants (`actk`, `bvtr`, `bvts`, `satmax`)
previously carried hand-added `f` suffixes, an undocumented 4-of-22 partial fix,
now absorbed into the toggle.

**Unfinished**: `wdm6_literal` covers only the header constants. There are **98
further float32-lossy inline literals** in the WDM6 C++ sources — 87 in
`ERF_AdvanceWDM6.cpp`, 11 in `ERF_InitWDM6.cpp` — including Euler-Mascheroni
`0.577215664901532`. Each is a latent divergence. `roqimax`'s `2.08e22` was one
of them and it alone accounted for the whole G13f residual. Do not wrap them
blindly: some sit in `#if 0` legacy blocks and three header constants (`r0`,
`peaut`, `alpha_wdm6`) have no matching Fortran parameter, so correspondence
must be checked case by case.

Integer versus real exponents are a second, distinct hazard in the same code.
Fortran `dimax**8` has an integer exponent and gfortran expands it to repeated
multiplication; C++ used `std::pow(dimax, Real(8.0))`, a libm call. For that
value the two agree bitwise, measured, but that is luck and not a guarantee.
`ERF_InitWDM6.cpp` now writes both `**8` and `**2` as explicit multiplication.

========================
ORDERED NEXT ACTIONS
========================

A) **Instrument `W3_STATECOPY` and resolve the step-2 discontinuity.**

   Confirmed structural asymmetry: `mp_wdm6_run_c` calls `mp_wdm6_run` which
   calls `wdm62D` **directly**, and `mp_wdm6_run` contains zero references to
   `itimestep`. The `wdm6` driver's `if (itimestep .eq. 1) nn = ccn0` block at
   `ERF_module_mp_wdm6.F90:216-224` is therefore **unreachable through the
   bridge**. The native path emulates it with `m_nn_initialized`
   (`ERF_InitWDM6.cpp:69,103-129`, cleared in `ERF_UpdateWDM6.cpp:54`). Two
   different mechanisms for reseeding CCN; disagreement about when they fire
   would present exactly as a step-2 discontinuity.

   **A tension must be resolved before acting on that.** The G1c PRE tag reports
   `nn` as 1.0e8 on the bridge leg against 1.8405808145911936e+05 on the native
   leg at step-2 entry, a 9.98e+07 difference matching the reported group
   divergence — yet `plt00001` is bitwise equal on all 32 variables including
   `nn`, so the conserved state entering step 2 is identical. Under corpus axis 5
   resolve structure before interpreting values. The unconverted Tier 1 tags are
   single-cell at `kts`; the isohelper silently substitutes `its` for the target
   column when it falls outside the tile (`isohelper:110`) and suppresses debug
   on non-target `j` slices (`isohelper:156`). The tag and the plotfile may not
   be sampling the same cell.

   Do this with paired PRE/POST records around `Copy_State_to_Micro`, not by
   inferring from group tags.

B) **Finish the Tier 1 all-k contract.** 14 of 33 groups converted (G8, G10a-e,
   G11, G12, G13a-f). Still single-cell at `kts`: P0, G1a-G7, G13g, G14-G17.
   G9 is **deliberately excluded** — it is surface-scoped, its native inputs are
   `box2d` slabs and its Fortran payload is 1D surface arrays, so an all-k loop
   would manufacture ~100 false differences per step. Everything from G13g on is
   downstream of the current frontier and cannot be an origin, so P0-G7 is the
   priority half.

C) **Option B validation**, if the build in the other terminal completed. See
   "not recoverable" below.

D) The remaining 98 inline literals, case by case.

E) Only after the above: Milestone C.

========================
STATE NOT RECOVERABLE FROM THE REPO
========================

- **An Option B build was started in another terminal** into
  `build_wdm6_doubleliteral` with `-DERF_WDM6_F32_LITERALS=OFF
  -DERF_WDM6_FORTRAN_REAL8=ON`. Status unknown at handoff. Before trusting it,
  check the flag was scoped:
  `grep -rl "fdefault-real-8" build_wdm6_doubleliteral/Exec/CMakeFiles/erf_srclib.dir/`
  should hit only `ERF_module_mp_wdm6*`. `tools/build_ci` hardwires
  `BUILD_DIR="$ERF_ROOT/build_wdm6"` and has no `--build-dir` flag, so either add
  one or run the two legs directly.
  Expected: `diameter` bitwise equal on both legs but at the true-double value.
  The point of the test is whether anything *else* moved — `work1i` and `xmi`
  come from `sqrt` and `**1.31` chains and are where a lowering change would show.
- **`build_wdm6` is currently in F32 mode** (`ERF_WDM6_F32_LITERALS=ON` in the
  cache). Every measurement in this handoff was taken that way.
- **A dcg config change was requested and never applied.** `general.hook_timeout_ms`
  is at its 1000 ms default and `~/.config/dcg/config.toml` does not exist. The
  write was blocked by the auto-mode classifier as a security-relevant change.
  Jean must create it manually; the content is in the session transcript. Timeout
  exhaustion **blocks** the command, which caused repeated friction on long
  pipelines and heredocs.
- `tools/build_ci` is now **tracked** (`b96258e12`), so its results are no longer
  disqualified from formal use on the untracked-asset rule. Phase 1 rows remain
  triage for the independent reason that they came from a dirty worktree. The
  harness gained `--target-column`, `--max-step`, `--debug`, per-k comparison
  with the k of the maximum, and first/last block passes so growth across
  timesteps is visible in one run. `mpirun` was removed to match campaign posture.
- `plotfile-numconc-fix` is a local branch, 2 commits off `development`, holding
  the column-rotation fix and the core `nn` support for upstream. **Never
  compiled, stale base.** Morrison is exposed to that same rotation bug and no
  upstream issue has been filed.

========================
OPERATOR NOTES
========================

- **Builds.** The box reports 20 cores but ~11 GB RAM and Jean runs several
  worktrees concurrently from separate sessions. Two `-j4` builds OOM-killed each
  other twice. Check `free -g` and `ps -o pid,etime,rss,args -C cc1plus` first —
  the paths reveal which worktree owns each compile; never kill another session's.
  `-j1` when anything else is running. Touching `ERF.H` or `ERF_WDM6.H` triggers
  a near-full ~200 TU rebuild that exceeds the Bash timeout: run it in the
  background. A single scheme file is one TU and finishes in a minute.
- **Provenance.** `AMReXBuildInfo.cmake:181` captures `git describe` at CMake
  *configure* time, so build-info tier 1 can never repair a stale hash on a CMake
  build. Go straight to tier 2 (remove `build_wdm6/erf_srclib/AMReX_buildInfo.cpp`).
  Watch for the inverse hazard too: after a plain `cmake --build` on a dirty tree
  the binary reports a *clean* hash because buildinfo was not regenerated.
- **Milestone B is currently useless as a discriminator.** Its divergence is
  dominated by the step-2 event, so any microphysics change reads `1.000x` there.
  Step 1 and the tag frontier are the only meaningful signals until step 2 is clean.
- Commit messages: subject line only. Rationale belongs in the ledger rows.
