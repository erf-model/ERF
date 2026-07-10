# Spec: driver source reorganization (developer-friendliness)

> Status: living document · Owns: the *file layout* of the ERF-side Noah-MP
> driver and the field-registry convention. · Companion:
> [`spec-noahmp-api.md`](spec-noahmp-api.md) (the class & lifecycle),
> [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md) (per-step data movement),
> [`spec-noahmp-io.md`](spec-noahmp-io.md) (restart & output).

## 1. Why this exists

The driver was, historically, a single 824-line `ERF_NOAHMP.cpp` plus a 494-line
`ERF_NOAHMP.H`. Every other large ERF class (`Microphysics/SAM`,
`TimeIntegration`, `Diffusion`) is split across concern-named translation units;
Noah-MP was the outlier. Two concrete problems followed:

- **Navigability.** `Init` (~270 lines) and `Advance_With_State` (~520 lines)
  were monoliths of 6–7 inline phases each.
- **The four-enum sync hazard** — the API spec's "most common coupling bug."
  Adding one coupled field meant editing 4–6 disjoint, index-aligned sites (the
  `NoahmpInputComp`/`NoahmpOutputComp`/`LsmData_NOAHMP` enums, two device kernels,
  two host copy loops, and a hand-maintained `fixed_names` table matched to the
  enum order only by a runtime assert).

This reorganization is a **pure refactor**: no physics, component order, buffer
sizing, or output changes — a Noah-MP run reproduces its prior trajectory
bitwise.

## 2. File layout

| File | Contents |
|------|----------|
| `ERF_NOAHMP.H` | The `NOAHMP` class: members, inline accessors, private-helper declarations. Includes `ERF_NOAHMP_Fields.H`. |
| `ERF_NOAHMP_Fields.H` | The coupling **contract**: the five field `enum`s and the X-macro field registry (§4). |
| `ERF_NOAHMP_Init.cpp` | `NOAHMP::Init`. |
| `ERF_NOAHMP_Advance.cpp` | `NOAHMP::Advance_With_State` + the per-step helpers (`interp_from_lev0`, `time_to_fire`, `stage_forcing`, `read_results`). |
| `ERF_NOAHMP_Precip.cpp` | precip source collection, snapshot alloc/seed/advance, and clamp/invariant reporting. |
| `ERF_NOAHMP_IO.cpp` | `Plot_Landfile`, `Write_Lsm_Restart`, `Read_Lsm_Restart`. |

Both build systems list these explicitly: `Make.package` (`CEXE_sources`) and
`CMake/BuildERFExe.cmake` (the `ERF_ENABLE_NOAHMP` block). Adding a `.cpp` means
editing **both**.

## 3. `Advance_With_State` structure

`Advance_With_State` is now a thin orchestrator; the phases described in
[`spec-noahmp-gpu.md`](spec-noahmp-gpu.md) map to private methods:

- `interp_from_lev0(...)` — the `!m_has_nc_file` fine-level interpolation branch.
- `time_to_fire(elapsed_time)` — the firing-schedule gate + lockstep counter
  advance (identical on every rank; see api spec §5).
- `prepare_precip_sources(...)` / `advance_precip_snapshots(...)` /
  `report_precip_diagnostics(...)` — precip bookkeeping (`_Precip.cpp`).
- per box: `stage_forcing(...)` (gpu-spec steps 1–3) then `DriverMain()` (step 4)
  then `read_results(...)` (steps 5–6).

The GPU invariants in gpu-spec §6 are unchanged (sync placement, one component
per field, pinned arena, k/j transpose, counter mirrored-not-owned).

## 4. The field registry (X-macro)

`ERF_NOAHMP_Fields.H` defines the mechanical, 1:1 coupled fields once as X-macro
lists, then expands them to generate the enum entries, the `LsmData_NOAHMP` name
strings, and the host↔`NoahmpIO` copy lines — a single source of truth that
cannot drift:

```c
// (name-token, noahmp-member-token)
#define NOAHMP_INPUT_3D(X)  X(u_phy,U_PHY) X(v_phy,V_PHY) X(t_phy,T_PHY) \
                            X(qv_curr,QV_CURR) X(p8w,P8W)
#define NOAHMP_INPUT_2D(X)  X(swdown,SWDOWN) X(glw,GLW) X(coszen,COSZEN)
#define NOAHMP_OUTPUT_2D(X) X(hfx,HFX) X(lh,LH) X(tau_ew,TAU_EW) X(tau_ns,TAU_NS) \
                            X(tsk,TSK) X(emiss,EMISS) X(o_grdflx,GRDFLX) \
                            X(o_fira,FIRAXY) X(o_sav,SAVXY) X(o_sag,SAGXY) \
                            X(o_albedo,ALBEDO) X(o_sfcrunoff,SFCRUNOFF) \
                            X(o_udrunoff,UDRUNOFF)
```

The 3-D lists copy with the k/j transpose (`noahmpio->MEMBER(i,1,j)`); the 2-D
lists use `(i,j)`.

### What the table does *not* cover (by necessity — stays explicit)

A fully data-driven runtime table is impossible here: Fortran-backed members are
named (not integer-indexable) and several fields carry per-field math. These stay
hand-written, sectioned, and commented:

- **Computed forcing:** winds (face→center average), `t_phy`/`p8w` (EOS helpers).
- **Precipitation:** the per-slot delta, the clamp guard, and the
  `MP_SNOW+MP_GRAUP ≤ MP_RAINNC` invariant rescale.
- **Banded albedos:** `ALBSFCDIRXY`/`ALBSFCDIFXY` at band 1 (VIS) / 2 (NIR).
- **Sentinels:** `smstav`/`smstot` emit `lsm_undefined` (not computed by this core).
- **Guarded outputs:** flux `÷rho` + the `-9999` fill-value guard.
- **Soil profile:** the per-`nsoil` `SMOIS`/`SH2O`/`TSLB` loops.

So: adding a *mechanical* field is one X-macro line + one enum position; adding a
*computed* field follows the explicit path (api spec §6, still authoritative for
the general workflow).

## 5. Invariants (do not regress)

1. **Component order & counts unchanged.** `NoahmpInputComp::NumComps`,
   `NoahmpOutputComp::NumComps`, and `LsmData_NOAHMP::NumVars` keep their values;
   the X-macro must not renumber components or the pinned buffers desync.
2. **Both build systems updated together** for any file add/remove.
3. **Bitwise reproducibility** — refactor only; verify against a pre-change
   baseline plus a restart run (io-spec §1).
