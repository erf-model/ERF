# WDM6 Bubble Plotfile-Parity Campaign

Scheme- and case-specific campaign definition required by Rule 35
("`fcompare` binary location ... must be defined in scheme notes";
"Milestones are scheme- and case-specific. Define them in implementation notes").

This is the WDM6/Bubble analogue of the WSM6/SquallLine campaign appendix in
`source-corpus/wsm6_implementation_notes.md:420-580`.

## Lane and scope

- `validation_lane = plotfile_parity_serial`
- serial, single-rank, CPU, non-MPI executable (Rule 35 scope gate)
- MPI and GPU parity are deferred to distinct lanes with their own gates

## Run root and assets

- run root: `Exec/RegTests/Bubble`
- input: `inputs_BF02_moist_bubble` with `erf.moisture_model=WDM6` override
- executable: `build_wdm6/Exec/erf_exec` (`CMAKE_BUILD_TYPE=Release`,
  `ERF_ENABLE_MPI=OFF`) — already satisfies the serial scope gate
- campaign runner: `Exec/.DevTests/fortran-to-cpp-microphysics/references/applications/wdm6/campaign_bubble_wdm6`
- **no `mpirun`** for any run in this campaign
- `GFORTRAN_UNBUFFERED_ALL=y` on the Fortran leg only

### fcompare binary

Rule 35 requires this note to name the binary. Either is acceptable:

- `tools/fcompare_serial`
- `Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex` — built from the tracked
  submodule, and the binary the WSM6 campaign used

Verified equivalent on the Milestone A pair: both report identical values to
every printed digit (`temp` 6.631607696 / 2.290772895e-2, `rhotheta` 3.01690336,
`pressure` 867.1722896, `qi` 5.347781098e-07). Milestone results are therefore
independent of the compare tool, and the missing `ni` column is an ERF plotfile
emission gap rather than a tooling artifact.

Both are build artifacts and neither is tracked; prefer the AMReX one when
provenance matters, since it is reproducible from tracked submodule source.

Do not edit the tracked input: it is shared with the nightly `Bubble_WSM6`
regression test. Apply all campaign settings as run-time overrides.

## dt convention

Keep the input's `erf.fixed_dt = 0.5`. Do not adopt WSM6's `fixed_dt=1.0` —
0.5 is what the nightly Bubble page runs, and step-100 equivalence is the point.

## Compared variable set

The input's `erf.plot_vars_1` carries 29 variables and stops at `rhoQ6`, so the
double-moment state is invisible to `fcompare` by default. Override it to append
number concentrations:

```
erf.plot_vars_1="<29 input vars> nc ni nr"
```

Use `nc nn nr`, not `nc ni nr`. WDM6 carries **no ice number**: `RhoQ8` is
CCN/total aerosol `nn`.

Earlier revisions of this note claimed "ERF's plot label for that slot is `ni`.
Do not rename anything." That was wrong on both counts and is retracted. ERF has
never labelled WDM6's `RhoQ8` as `ni` — `ERF_DataStruct.H` sets `ni = -1` for
WDM6 with the comment *"ice number (ni) not used - WDM6 uses nn for total
aerosol"*, so the selection gate drops `ni` and the slot went unnamed entirely.
`ni` is bound to cloud ice number ERF-wide by Morrison
(`ERF_Morrison.H:47`, `ERF_UpdateMorrison.cpp:55-59`), and the two are
physically different: every other number concentration pairs with a mass species
(`nc`/`qc`, `nr`/`qr`, `ni`/`qi`), while `nn` has none — it is an aerosol
reservoir that exchanges number with `nc` and `nr` through activation
(`ERF_module_mp_wdm6.F90:2790`) and evaporation (`:2835`, `:1855`), bounded to
1e8-2e10 m^-3 (`:602`). Magnitudes differ by three to five orders.

`nn` is now a first-class plot variable: `MoistureComponentIndices::nn`, a
`ccn_number` capability, and a writer block that reads the slot from the
per-scheme registry rather than hardcoding `RhoQ8`. Morrison is untouched.

**Resolved:** the earlier "31 variables, not 32" gap had two independent causes,
both fixed. `ni` was correctly absent (WDM6 has no ice number) — the missing
field was `nn`, which had no name at all. Separately, the columns after `qgraup`
were rotated because the writer emitted `nc..ng` before `qt/qn/qp/qsat` while the
header names came from `derived_names` in `ERF.H`, which declares the reverse;
`qt`, `qp`, `nc` and `nr` each reported another field's data. That defect shipped
with PR #3218 (30 May 2026) and also affects Morrison. Any `qt`/`qp`/`nc`/`nr`
value recorded before that fix is void.

## Milestone table (WDM6 Bubble)

`plot_int_1=100` makes plot steps multiples of 100; A and B use `plot_int_1=1`.

| Milestone | Stop | Plot step | t (s) | Role |
| --- | --- | --- | --- | --- |
| A | step 1 | 1 | 0.5 | copyback / initial-state check |
| B | step 10 | 10 | 5 | first sequential accumulation |
| C | step 100 | 100 | 50 | nightly `Bubble_WSM6` equivalence point; PR gate |

Run structure:
- short pair A-B: `max_step=10`, `plot_int_1=1`, `check_int=1`
- C pair: `max_step=100`, `plot_int_1=100`, `check_int=100`

## Coverage: what Bubble does and does not exercise

Bubble is **not** warm-phase-only. Measured at Milestone A, `qi` and `qsnow`
carry nonzero bridge-vs-native differences and the max-error zone is
`(i,j,k) = (199,3,90)` — roughly 9 km, well below freezing. Cold-phase groups
are active at upper levels even though the initial bubble state has no ice.

Do not claim Bubble is a warm-rain-only surface. It is, however, a weaker
cold-phase surface than SquallLine, where `qi/qs/qg` become active throughout
the domain; keep SquallLine as the cold-phase gate.

## Diagnostic surface caveat (blocks frontier trust)

WDM6 group tags print **one cell**: every field is indexed at `kts` with no
`do k` loop, e.g.

```fortran
'WDM6-FORT_PRE_G7', kts, qrs(i_dbg_local,kts,1), ... t(i_dbg_local,kts), ...
```

The WSM6 contract on `wsm6-mpi-cpu` instead loops the full column:

```fortran
if (mpdbg_level >= 1 .and. loop == 1) then
  do k = kts, kte
    write(*,'(A,I3,6E24.16)') 'WSM6-FORT_DENFAC ', k, denfac(its,k), ...
  enddo
endif
```

Consequence: each WDM6 group closure rests on ~1/100th of the evidence its WSM6
counterpart did, and no group evidence exists at the k where Milestone A fails.
Restore the WSM6 all-k print contract before trusting the group frontier or
attributing any plotfile residual to a group.

Exception: `G10A`/`G10B`/`G10C` emit ~20,000 Fortran lines (full i-k plane)
while the C++ side emits one line per group, so those groups have no matched
comparison surface either.

## Campaign directory discipline

Unique `campaign_id` root per execution, never reused. `.gitignore` carries
`plotfile_parity_*/` so run artifacts and logs are untracked and clean-SHA reruns
need no destructive cleanup.

Naming follows the WSM6 convention promoted into `process-operators.md`. The
WDM6/Bubble instantiation:

| Field | Value |
| --- | --- |
| campaign root | `plotfile_parity_wdm6_bubble_v<N><letter>` |
| `campaign_id` | same as the root |
| `run_pair_id` | `<campaign_id>_short<max_step>`, or `<campaign_id>_r36_chk<NNNNN>_step<N>` |
| `execution_id` | `<UTC>Z_plotparity_wdm6_bubble_v<N><letter>_<pair>_<leg>_dbg<level>` |
| `run_pair` | `short` / `r36` / `r37` — the runner's second argument |

Scheme and case stay in the name because `campaign_id` is a ledger column read
out of context. Omit the clean-SHA token WSM6 used on retreat sub-runs when a
whole campaign sits at one SHA; `git_sha` already carries it in full.

This supersedes the earlier `campaign_wdm6_bubble_v*` roots. Those were invented
without a corpus check and their numbering is not a reliable index into the
ledgers — v3/v4/v5 were run but never recorded. Rows and reports already written
against `campaign_wdm6_bubble_v1`/`v2` remain valid history under the naming in
force at the time; do not rewrite them.

## Command templates

Fortran (oracle) leg:

```bash
cd Exec/RegTests/Bubble
GFORTRAN_UNBUFFERED_ALL=y ../../../build_wdm6/Exec/erf_exec inputs_BF02_moist_bubble \
  erf.moisture_model=WDM6 erf.use_wdm6_cpp_answer=0 erf.microphysics_debug=0 \
  max_step=<N> erf.plot_int_1=<P> erf.check_int=<P> \
  erf.plot_vars_1="$PV" \
  erf.plot_file_1=<campaign_id>/<pair>_fortran/plt \
  erf.check_file=<campaign_id>/<pair>_fortran/chk \
  > <campaign_id>/log.<pair>_fortran 2>&1
```

Native leg: identical with `erf.use_wdm6_cpp_answer=1`, no `GFORTRAN_UNBUFFERED_ALL`.

Per-milestone compare (zero tolerance first):

```bash
<fcompare_binary> <campaign_id>/<pair>_fortran/plt<N> \
                      <campaign_id>/<pair>_cpp/plt<N>
```

On first fail, archive immediately:

```bash
<fcompare_binary> -z <var> <fort_plt_N> <cpp_plt_N>
<fcompare_binary> -d <var> <fort_plt_N> <cpp_plt_N>
mv diffs <campaign_id>/diffs_plt<N>_<var>
```

## Variable-to-tag mapping (WDM6 groups)

Adapted from the WSM6 table (`wsm6_implementation_notes.md:671`) onto WDM6
group ids:

- `qr` -> G5b (NISLFV_R), G13a (praut/pracw/prevp)
- `qs`/`qg` -> G5c (NISLFV_SG), G13b/G13c (psaci/psacw), G13d (pracs), G13f (psaut), G13g (psevp)
- `qi` -> G8 (VICE), G13b (praci), G13e (pidep), G13f (psaut)
- `qv` -> G13a (prevp), G16b (PCOND), G1c (QSAT), G15 (QSAT2)
- `t` -> G7 (MELT), G10a-G10d (PHASE), G14 (UPDATE), G16b (PCOND)
- `nc`/`nn`/`nr` -> G13a (nraut/nracw/nccol/nrcol), G16a (di82 collapse), G16b (activation), G17 (lambda-bound ncr recompute)
- `den` -> G4/G6/G11 (SLOPE1/2/3)

Retreat order is the canonical file order in `group_map.md`.

Sites **outside** every group tag, and therefore invisible to the frontier —
check these before group retreat. They now have a declared namespace and tag
surface; see "Bridge-wrapper sites" in `group_map.md`:
- `W1_THETAWB` — the T -> theta conversion after the bridge call
  (`ERF_AdvanceWDM6.cpp` Fortran branch only). Instrumented.
- `W2_PACKUNPACK` — pre-kernel pack / post-kernel unpack in `mp_wdm6_run_c`
- `W3_STATECOPY` — `Copy_Micro_to_State` / `Copy_State_to_Micro`
- `W4_LEGACY` — the retained `#if 0` pre-group regions in the native path
  (Steps 4-8, 9, and now 10 and 11). Compiled out, so eliminated by inspection.

These are not `G*` groups and must not be added to the `P0..G17` file-order
sequence: they have no Fortran counterpart at all.
