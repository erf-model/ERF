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
- `collect_precip_sources(...)` / `prepare_precip_snapshots(...)` /
  `advance_precip_snapshots(...)` / `report_precip_diagnostics(...)` — precip
  bookkeeping (`_Precip.cpp`).
- per box: `stage_forcing(...)` (gpu-spec steps 1–3) then `DriverMain()` (step 4)
  then `read_results(...)` (steps 5–6).

The GPU invariants in gpu-spec §6 are unchanged (sync placement, one component
per field, pinned arena, k/j transpose, counter mirrored-not-owned).

## 4. The field registry (X-macro)

`ERF_NOAHMP_Fields.H` defines the mechanical, 1:1 coupled fields once as
**field registries** — `NOAHMP_*_FIELDS(X)` lists of `X(...)` rows — then expands
each with a small **emitter** macro (the per-row action) to generate the enum
entries, the `LsmData_NOAHMP` name strings, and the host↔`NoahmpIO` copy lines.
A single source of truth that cannot drift:

```c
// Registries. Boundary rows are  X( ERF comp, NoahmpIO member ).
#define NOAHMP_INPUT_3D_FIELDS(X)  X(u_phy,U_PHY) X(v_phy,V_PHY) X(t_phy,T_PHY) \
                                   X(qv_curr,QV_CURR) X(p8w,P8W)
#define NOAHMP_INPUT_2D_FIELDS(X)  X(swdown,SWDOWN) X(glw,GLW) X(coszen,COSZEN)
// ... NOAHMP_OUTPUT_2D_FIELDS(X), NOAHMP_LSMDATA_FIELDS(X) (comp-only rows) ...

// Reusable emitters (defined once) -> no #define/#undef ceremony per site.
#define NOAHMP_ENUM(comp)          comp,    // (comp)         -> enum entry
#define NOAHMP_QUOTE(comp)         #comp,   // (comp)         -> "comp"
#define NOAHMP_COMP(comp, member)  comp,    // (comp, member) -> enum entry
```

An expansion site is then one clean line, e.g. `NOAHMP_INPUT_2D_FIELDS(NOAHMP_COMP)`
inside the enum, or `NOAHMP_LSMDATA_FIELDS(NOAHMP_QUOTE)` for the name table.
The host↔`NoahmpIO` copy loops reference local buffers, so their emitters
(`NOAHMP_STAGE_IN_3D/2D`, `NOAHMP_READ_OUT_2D`) are defined next to each loop in
`_Advance.cpp`: 3-D members copy with the k/j transpose (`noahmpio->MEMBER(i,1,j)`);
2-D members use `(i,j)`.

A fourth registry, `NOAHMP_RESULT_FIELDS`, ties the Noah-MP result scatter
together. Each row is a triple `X( ERF Array4 alias, LsmData_NOAHMP comp,
NoahmpOutputComp comp )` and drives the three index-aligned sites in
`_Advance.cpp` that previously drifted by hand: the `Array4` bind
(`NOAHMP_BIND_RESULT`), the valid-cell copy (`NOAHMP_COPY_RESULT`), and the
`-9999` fill-value → `lsm_undefined` branch (`NOAHMP_UNDEF_RESULT`). Only the
mechanical (pure 1:1) results are listed; the per-field math around them stays
explicit (below).

### What the table does *not* cover (by necessity — stays explicit)

A fully data-driven runtime table is impossible here: Fortran-backed members are
named (not integer-indexable) and several fields carry per-field math. These stay
hand-written, sectioned, and commented:

- **Computed forcing:** winds (face→center average), `t_phy`/`p8w` (EOS helpers).
- **Precipitation:** the per-slot delta, the clamp guard, and the
  `MP_SNOW+MP_GRAUP ≤ MP_RAINNC` invariant rescale.
- **Banded albedos:** the `NoahmpIO` *read* of `ALBSFCDIRXY`/`ALBSFCDIFXY` at
  band 1 (VIS) / 2 (NIR) stays explicit; their subsequent scatter into ERF *is*
  in `NOAHMP_RESULT_FIELDS` (that part is a plain 1:1 copy).
- **Sentinels:** `smstav`/`smstot` emit `lsm_undefined` (not computed by this
  core), so they are kept out of `NOAHMP_RESULT_FIELDS` and stay explicit.
- **Fluxes:** the surface-layer flux `÷rho` conversions behind the `-9999`
  guard (the guard *branch itself* now drives the table-driven result scatter).
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
