# Spec: driver-owned I/O — restart, output, and static inputs

> Status: living document · Owns: checkpoint/restart, the per-step land
> plotfile, and the static-input reads, as seen from the ERF driver. · Companion:
> [`spec-noahmp-api.md`](spec-noahmp-api.md),
> [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md). · Files: `ERF_NOAHMP.cpp`,
> `ERF_NOAHMP.H`. · User-facing docs:
> <https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html#checkpoint-and-restart>.

> **Where the bytes are written.** The NetCDF reading/writing itself lives on the
> Noah-MP side (`NoahmpWriteLandMod`, `NoahmpWrite/ReadRestartMod`, `ReadLand*`,
> `ReadNamelist`, `ReadTable`) and is specified under
> `Submodules/Noah-MP/drivers/erf/dev/`. This spec covers what the *driver*
> orchestrates: when those routines run and what ERF-side state accompanies them.

## 1. Checkpoint / restart

Noah-MP participates in ERF's standard checkpoint/restart (`erf.check_file`,
`erf.check_int`, `amr.restart`) with no extra input options. A restart is
**bitwise reproducible** — the restarted run matches an uninterrupted run to the
last bit — which requires persisting two pieces of state beyond the ordinary ERF
`MultiFab`s.

### 1a. The substep counter

When Noah-MP subcycles (`NOAH_TIMESTEP > dt`; see
[`spec-noahmp-api.md`](spec-noahmp-api.md) §5), the firing schedule is driven by
`m_itimestep`. ERF writes it to `<chkfile>/lsm_step` (one value per AMR level)
via `Get_LSM_Step()` and restores it with `Set_LSM_Step(step)`.

`Get_LSM_Step` returns the **class member** `m_itimestep`, not
`noahmpio_vect[0].itimestep` — the class member is valid on every rank, including
land-free ranks whose `noahmpio_vect` is empty. `Set_LSM_Step` writes
`m_itimestep` *and* fans the value out to every box, so the schedule resumes
identically across ranks.

### 1b. The full prognostic state

The complete land state — soil temperature/moisture, the snowpack (active layers
+ layer count), canopy/vegetation, aquifer, albedo history, phenology,
accumulators — is serialized at the model's working precision to
`<chkfile>/noahmp_restart/Level_<lev>.nc`. The driver exposes this through two
`NullSurf` overrides that simply fan out over the boxes:

- `Write_Lsm_Restart(dir)` → `noahmpio.WriteRestart(dir)` for each box.
- `Read_Lsm_Restart(dir)`  → `noahmpio.ReadRestart(dir)` for each box.

Each local box writes its own tile into the single global-domain NetCDF file, in a
collective operation that every rank must enter together.

### 1c. Restart ordering

On restart, ERF **first cold-initializes** the land state from the WRF input file
and tables (exactly as `Init()` does for a fresh run), **then overwrites** it
with the checkpointed state via `Read_Lsm_Restart`. Noah-MP's per-step input
transfer (stages 1–3 of [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)) pulls the
restored state into the physics on the first `Advance`.

> **Legacy checkpoints.** Checkpoints predating this capability lack `lsm_step`
> and `noahmp_restart/`. ERF detects their absence, resets the counter to zero,
> cold-initializes the state, and warns that the restarted trajectory will differ
> from a cold start (i.e. *not* bitwise reproducible).
>
> **Archiving.** `lsm_step` and `noahmp_restart/` are plain files inside the
> checkpoint directory, *not* AMReX `MultiFab`s. A manual copy/archive of a
> checkpoint must include them.

## 2. Per-step land output

`Plot_Landfile(nstep)` fans out `noahmpio.WriteLand(nstep)` over the boxes to emit
the per-timestep NetCDF land diagnostics. The cadence is owned by the caller, not
the driver: ERF's time loop calls `Plot_Landfile` for every level whenever it
writes the level-1 3-D plotfile (governed by `erf.plot_int_1`/`erf.plot_per_1`;
see `ERF.cpp`). `Init()` does `pp.query("plot_int_1", m_plot_int_1)`, but
`m_plot_int_1` is currently unused — the driver does not gate output itself. The
initial land file is written during `Init()` itself with tag `0` (`WriteLand(0)`).
Which fields appear in that file is controlled Fortran-side in
`NoahmpWriteLandMod.F90`.

## 3. Static inputs read during `Init()`

The driver requires two files in the run directory — `namelist.erf` (Noah-MP
parameters, incl. `NOAH_TIMESTEP`) and `NoahmpTable.TBL` (parameter tables) —
plus the NetCDF land file (Noah-MP currently requires `erf.init_type =
"WRFInput"`). Per box, `Init()` invokes, in order: `ReadNamelist()`,
`ReadLandHeader()`, `ReadTable()`, `ReadLandMain()`. These set up the
Fortran-owned arrays before `InitMain()` computes any remaining initial values.
The `DTBL` value parsed from the namelist is broadcast to `m_dtbl` so the firing
schedule is defined on every rank (see [`spec-noahmp-api.md`](spec-noahmp-api.md)
§4–5).

## 4. Invariants (do not regress)

1. **Counter from the class member.** `Get_LSM_Step` must return `m_itimestep`,
   never `noahmpio_vect[0]` — the latter is absent on land-free ranks.
2. **Cold-init then overwrite.** Restart relies on `Init()` running first and
   `Read_Lsm_Restart` overwriting it; do not skip the cold init on restart.
3. **Collective writes.** Restart and land output are written collectively over
   the boxes; every rank must reach them (they sit alongside the collective
   `FillBoundary` in the step).
4. **Restart files travel with the checkpoint.** Keep `lsm_step` and
   `noahmp_restart/` inside the checkpoint directory.
