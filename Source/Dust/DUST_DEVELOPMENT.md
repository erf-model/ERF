# ERF-Dust Hazardous Particulate Dispersion Module — Development Log

## Module Purpose

The ERF-Dust module simulates wind-driven mineral dust emission, transport, and deposition
from mine tailings, haul roads, blast events, and chemically-altered surfaces in the
atmosphere. It is tightly integrated with the ERF mesoscale atmospheric model (MRF PBL,
MOST surface layer) and produces regulatory output for EPA NAAQS and MSHA compliance.

**Reference scenarios:**
- Marticorena & Bergametti (1995): Saltation-based dust emission from arid soils
- Shao & Lu (2000): Threshold friction velocity for wind erosion
- Zhang et al. (2001): Dry deposition resistance scheme for coarse particles
- Hong & Pan (1996): MRF PBL for turbulent vertical mixing of dust scalars

---

## Twenty-Three-Phase Implementation Roadmap

| Phase | Title | Key Deliverables | Status |
|-------|-------|------------------|--------|
| 1 | **Directory Scaffold & CMake** | Compilable stub; parameter reading; grid construction | ✅ COMPLETE |
| 2 | **2D Dust Grid Definition** | `ERF_DustGrid.H/cpp`; refined 2D slab from ATM level-0 | ✅ COMPLETE |
| 3 | **Surface Property Map Reader** | `ERF_DustSurfaceReader.H/cpp`; soil type, silt, crust, moisture, suppression rasters | ✅ COMPLETE |
| 4 | **PHREEQC Output File Reader & u\*t Mapper** | `ERF_PhreeqcReader.H/cpp`; offline geochemical coupling | ✅ COMPLETE |
| 5 | **Threshold Friction Velocity** | `ERF_DustThreshold.H`; Bagnold/Shao/Marticorena formulation | ✅ COMPLETE |
| 6 | **Saltation & Emission Flux** | `ERF_DustEmission.H/cpp`; MB95 vertical flux | ✅ COMPLETE |
| 7 | **Blast Schedule Emission Source** | `ERF_DustBlastSchedule.H`; time-indexed events + MPI broadcast | ✅ COMPLETE |
| 8 | **Suppression Agent Degradation** | `ERF_DustSuppression.H`; temperature/wind-driven half-life decay | ✅ COMPLETE |
| 9 | **Wind Extraction from ERF-Atm** | `ERF_DustWindExtract.H/cpp`; zref interpolation; SurfaceLayer u* mapping | ✅ COMPLETE |
| 10 | **ERF-Dust → ERF-Atm Aerosol Injection** | `ERF_DustAtmCoupling.H/cpp`; coarsen + inject at k=0 | ✅ COMPLETE |
| 11 | **Multi-Bin Settling** | `ERF_DustSettling.H`; Stokes-Cunningham; first-order upwind flux divergence | ✅ COMPLETE |
| 12 | **Dry Deposition BC** | `ERF_DustDeposition.H`; Zhang et al. resistance scheme; accumulator | ✅ COMPLETE |
| 13 | **Two-Way Coupling: ATM → Dust** | `ERF_DustAtmReturn.H`; loading feedback + moisture inhibition | ✅ COMPLETE |
| 14 | **MRF Countergradient Extension** | `EddyDiff::Scalar_v` in `ERF_ComputeDiffusivityMRF.cpp`; dust Schmidt number | ✅ COMPLETE |
| 15 | **MRF Scale-Aware Blending at Fine AMR Levels** | High-res diffusivity bounds; multi-level transport validation | 🔲 PLANNED |
| 16 | **Dust Diagnostics & Plotfile Output** | `ERF_DustPlotfile.H/cpp`; `ERF_DustStatsOutput.H`; `dust_diag.dat` | ✅ COMPLETE |
| 17 | **EPA NAAQS PM2.5/PM10 Compliance** | `ERF_DustPM.H`; 24h running average; exceedance flags; `dust_naaqs.csv` | ✅ COMPLETE |
| 18 | **MSHA Worker Exposure Tracking** | `ERF_DustMSHA.H`; `ERF_DustMSHAOutput.H`; 8h TWA; shift summary CSV | ✅ COMPLETE |
| 19 | **Lagrangian Super-Particle Tracking** | `ERF_DustPC.H/cpp`; source-receptor attribution; nearest-cell advection | ✅ COMPLETE |
| 20 | **Multi-Mine Domain / Per-Site PHREEQC** | `ERF_DustSiteRegistry.H`; site_id MultiFab; per-site deposition CSV | ✅ COMPLETE |
| 21 | **PHREEQC Deposition Feedback Writer** | `ERF_DustPHREEQCWriter.H`; periodic grid write + per-site summary CSV | ✅ COMPLETE |
| 22 | **Haul Road Vehicle Schedule Emission** | `ERF_DustRoadSchedule.H`; AP-42 Ch. 13.2.2 road resuspension | ✅ COMPLETE |
| 23 | **DOE Critical Materials Assessment** | `ERF_DustCriticalMaterials.H`; per-bin CM fractions; budget CSV | ✅ COMPLETE |
| 24 | **Regression Test Suite & Documentation** | `DustIntegration` end-to-end test; all phase tests audited | ✅ COMPLETE |

---

## Phase 1: Directory Scaffold & CMake Integration

### Overview

Phase 1 establishes a fully compilable, zero-physics foundation. The module:
- Reads all parameters from `erf.dust.*` ParmParse namespace (`ERF_DustParams.H`)
- Constructs a 2D dust computational grid refined from the atmospheric level-0 grid
- Allocates all MultiFabs for dust state variables
- Validates prerequisites (grid divisibility, parameter ranges)
- Emits structured debug output at every step

### All Parameters (`ERF_DustParams.H`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enable` | `false` | Enable dust module |
| `dust_debug` | `false` | Per-step debug output |
| `n_size_bins` | `3` | Number of particle size bins |
| `grid_ratio` | `1` | Dust grid refinement vs ATM level-0 |
| `particle_density` | `2650.0` | Bulk particle density [kg/m³] |
| `z0_dust` | `0.01` | Roughness length over bare tailings [m] |
| `silt_fraction` | `0.10` | Default surface silt mass fraction [-] |
| `threshold_A_coeff` | `0.0123` | Bagnold threshold coefficient A |
| `crust_index` | `0.0` | Default surface crust index [0,1] |
| `alpha_crust` | `0.5` | Crust reduction coefficient for u\*t |
| `alpha_efflor` | `0.3` | Efflorescence reduction coefficient for u\*t |
| `test_ustar` | `0.0` | Placeholder u\* before Phase 9 [m/s] |
| `rho_air` | `1.225` | Air density [kg/m³] |
| `zref` | `10.0` | Wind extraction height [m] |
| `supp_tau_base_s` | `3600.0` | Base suppression half-life [s] |
| `test_surf_temp_K` | `293.15` | Placeholder T\_sfc [K] |
| `test_wind_speed` | `0.0` | Placeholder 10m wind [m/s] |
| `atm_feedback` | `1.0` | Coupling strength [0,1] |
| `transport_bins_separately` | `false` | Separate 3D scalar per bin vs. single sum |
| `bin_diameters` | `{7e-6, 2.5e-6, 50e-6}` | Per-bin particle diameters [m] |
| `deposition_E0` | `3e-3` | Surface collection efficiency [-] |
| `loading_feedback_coeff` | `0.0` | Shao (2001) loading feedback |
| `use_dynamic_moisture` | `false` | Use Q1fx3 for moisture inhibition |
| `phreeqc_update_interval_s` | `86400.0` | PHREEQC re-read interval [s] |
| `blast_schedule_file` | `""` | Blast event CSV |
| `blast_reactivity` | `2.0` | Fresh blast surface reactivity multiplier |
| `road_schedule_file` | `""` | Haul road CSV |
| `cm_fractions` | `{}` | Per-bin critical material mass fractions |
| `dust_plot_int` | `-1` | Plotfile interval (-1 = disabled) |
| `dust_diag_file` | `"dust_diag.dat"` | Diagnostics CSV |
| `dust_naaqs_file` | `"dust_naaqs.csv"` | NAAQS compliance CSV |
| `msha_pel_mg_m3` | `5.0` | PEL threshold [mg/m³], 30 CFR 56.5001 |
| `msha_shift_duration_s` | `28800.0` | Shift duration [s] |
| `msha_exposure_file` | `"msha_exposure.csv"` | MSHA per-step CSV |
| `msha_shift_file` | `"msha_shift_summary.csv"` | End-of-shift summary |
| `phreeqc_feedback_interval_s` | `0.0` | PHREEQC deposition write interval [s] |
| `phreeqc_feedback_file` | `"dust_dep_feedback.dat"` | Deposition grid file |
| `phreeqc_site_summary_file` | `"dust_dep_site_summary.csv"` | Per-site CSV |
| `cm_budget_file` | `"dust_cm_budget.csv"` | CM budget CSV |
| `road_diag_file` | `"dust_road_diag.csv"` | Road diagnostics CSV |

### Dust Grid Structure (2D Slab)

```
atm level-0 grid: nx × ny × nz
  ↓ extract k=0 slab
  ↓ refine by grid_ratio in x,y
dust grid: (nx * grid_ratio) × (ny * grid_ratio) × 1
```

Same physical domain, finer surface resolution. `grid_ratio=1` means one dust cell per ATM column.

---

## Phase 2: 2D Dust Grid Definition

### File: `ERF_DustGrid.H/.cpp`

**`create_dust_grid(ba_atm, dm_atm, geom_atm, grid_ratio)`**:
1. Extract k=0 slice from each ATM box → `ba_2d`
2. Refine by `IntVect(grid_ratio, grid_ratio, 1)` → `ba_dust`
3. Reuse `dm_atm` (no box count change after refine)
4. Build refined 2D Geometry: `hi = old_hi * grid_ratio + (grid_ratio - 1)`; physical domain unchanged

**Key invariant:** Box `i` in `ba_dust` is always owned by the same rank as box `i` in `ba_atm`.
This eliminates inter-rank communication in all dust-to-ATM coupling operations.

---

## Phase 3: Surface Property Map Reader

### File: `ERF_DustSurfaceReader.H/.cpp`

Reads ESRI ASCII (`.asc`) or NetCDF (`.nc`) raster files and maps to the dust grid:
- `soil_type_file` → `dust_soil_type` (STATSGO codes; 100+ = mine-specific)
- `silt_fraction_file` → `dust_silt_fraction`
- `crust_index_file` → `dust_crust_index`
- `moisture_flag_file` → `dust_moisture_flag`
- `suppression_file` → `dust_suppression`

When a file path is empty, the corresponding MultiFab is filled with the uniform default from `DustParams`.

---

## Phase 4: PHREEQC Output File Reader

### File: `ERF_PhreeqcReader.H/.cpp`

Reads PHREEQC geochemical output and updates dust surface properties:
- `crust_index` → modifies u\*\_t via `alpha_crust`
- `silt_fraction` → controls emission flux magnitude
- `efflorescence` → additional u\*\_t reduction
- `suppression_mod` → modifies suppression half-life

Update triggered when `(time - m_last_phreeqc_update) >= phreeqc_update_interval_s`.
Pattern: rank-0 reads file; result stored in per-cell MultiFabs (no MPI needed — file read is 2D,
same on all ranks if they read the same file, or scatter is used for large files).

---

## Phase 5: Threshold Friction Velocity

### File: `ERF_DustThreshold.H`

**Formulation (Bagnold 1941 + Shao & Lu 2000 + Marticorena & Bergametti 1995):**
```
u*_t_base = A * sqrt(rho_p * g * d / rho_a)          [Bagnold, threshold velocity]
u*_t      = u*_t_base * (1 - alpha_crust * CI)        [crust reduction]
           * (1 - alpha_efflor * efflor)               [efflorescence reduction]
           * sqrt(1 + a_f * max(w - w', 0))            [Fecan moisture inhibition, Phase 13]
```

Called each timestep via `recompute_dust_ustar_t()` using current crust/moisture/suppression state.

---

## Phase 6: Saltation & Emission Flux

### File: `ERF_DustEmission.H/.cpp`

**Vertical emission flux (Marticorena & Bergametti 1995):**
```
F_e = C_d * rho_a / g * ustar_in^3 * (1 - (u*_t/ustar_in)^2)
    * (1 + u*_t/ustar_in) * silt_fraction   [kg/m^2/s]
```

Applied only where `ustar_in > u*_t`. Zero elsewhere.
Distributed equally across `n_size_bins` (bin weighting added in Phase 11+).

---

## Phase 7: Blast Schedule Emission Source

### File: `ERF_DustBlastSchedule.H`

**CSV format:**
```
time_s  cx_m  cy_m  radius_m  mass_kg_m2  [mineral_type]  [priority]
```

**MPI pattern:** Rank-0 reads CSV → `MPI_Bcast` of POD event array → all ranks apply identically.
Each event adds a Gaussian-weighted mass impulse to `dust_emission_flux` at the blast time.
`blast_reactivity` multiplier scales the fresh-surface emission relative to weathered tailings.

---

## Phase 8: Suppression Agent Degradation

### File: `ERF_DustSuppression.H`

**Exponential decay:**
```
suppression(t+dt) = suppression(t) * exp(-dt / tau_eff)
tau_eff = supp_tau_base_s * f(T_sfc, u_10m)
```

Temperature correction: `tau_eff` decreases at high surface temperature (faster evaporation).
Wind correction: `tau_eff` decreases at high wind (mechanical removal).
`dust_retreat_flag` marks cells where suppression drops below threshold (triggers re-emission eligibility).

**Note:** `amrex.fpe_trap_invalid = 0` is required in all test inputs because the suppression
exponential can produce a denormal float at the boundary of the treated/untreated interface.

---

## Phase 9: Wind Extraction from ERF-Atm to ERF-Dust

### File: `ERF_DustWindExtract.H/.cpp`

Three extraction functions (adapts fire module pattern from `Source/Fire/`):

**`fill_dust_ustar_from_surface_layer(dust_ustar_in, u_star, m_dg)`**
- Coarsens ATM `u_star` (2D, from `SurfLayer->get_u_star(0)`) to dust grid
- Pattern: `dust_field(i,j,0) = atm_field(i/C, j/C, 0)`

**`fill_dust_wind_from_interpolation(dust_wind_ref, xvel, yvel, z_phys_cc, m_dg, zref, nz)`**
- Vertical interpolation of face-staggered ATM velocity to `z_surf + zref`
- Finds bracket `z_cc(k_lo) <= z_target < z_cc(k_hi)`, linearly interpolates
- Maps ATM column `(i/C, j/C)` to each dust cell `(i, j)`

**`fill_dust_scalar_from_atm(dust_field, atm_field_2d, m_dg)`**
- Generic k=0 coarsening for T_sfc, PBLH scalars

**PBLH fix (Phase 9 post-merge):** MRF stores PBLH in `EddyDiff::Turb_lengthscale` at all k levels.
`SurfLayer->get_pblh()` is never populated by MRF (only MYNN25 calls `update_pblh()`).
Fix: after `ComputeDiffusivityMRF` MFI loop, copy `Turb_lengthscale(k=klo)` into `SurfLayer->pblh`.
Without this fix, `PBLH_max = 1e+150` (bogus_large_value).

---

## Phase 10: ERF-Dust → ERF-Atm Aerosol Injection

### File: `ERF_DustAtmCoupling.H/.cpp`

**`coarsen_dust_flux_to_atm(Q_dust_atm, Q_dust_or_tmp, geom_dust, geom_atm, grid_ratio)`**
- Coarsens per-bin flux from 2D dust grid to 2D ATM slab using `amrex::average_down`
- When `grid_ratio > 1`: averages from `dust_bin_tmp` (on dust grid) → `dust_flux_atm` (on ATM grid)
- When `grid_ratio = 1`: equivalent to direct copy

**`apply_dust_tendency_to_cc_source(cc_source, Q_dust_atm, z_phys_cc, geom_atm, comp, feedback)`**
- Injects at k=0 only: `src(i,j,k,comp) += feedback * flux(i,j,0) / dz_k0`
- Explicit one-step lag: flux from step n injected at step n+1

**`transport_bins_separately` flag:**
- `false` (default): all bins summed into `dust_bin_tmp`, injected into single component
- `true`: each bin `b` injected separately into `m_dust_scalar_comp + b`

**Grid ratio > 1 bug (fixed):** Old code called `MultiFab::Copy(*dust_flux_atm, *dust_emission_flux, b, 0, 1, ...)`.
This crashes when `grid_ratio > 1` because `dust_flux_atm` (ATM resolution) and
`dust_emission_flux` (dust grid resolution) have different BoxArrays.
Fix: accumulate into `dust_bin_tmp` (allocated on dust grid), then `average_down` to `dust_flux_atm`.

---

## Phase 11: Multi-Bin Settling

### File: `ERF_DustSettling.H`

**Stokes settling velocity with Cunningham correction (Allen & Raabe 1985):**
```
C_c = 1 + (2λ/d) * (A1 + A2 * exp(-A3 * d / λ))
v_s = (rho_p - rho_a) * g * d^2 * C_c / (18 * mu_a)
```

**First-order upwind flux divergence (applied to cc_source):**
```
tend(k) = -(v_s * rho_dust(k) - v_s * rho_dust(k-1)) / dz
```
At `k=klo`: flux from below = 0 (deposition handled by Phase 12).

Applied per bin with `dust_comp = m_dust_scalar_comp + b` when `transport_bins_separately = true`.

---

## Phase 12: Dry Deposition Lower BC

### File: `ERF_DustDeposition.H`

**Resistance scheme (Zhang et al. 2001):**
```
v_d = 1 / (Ra + Rb + Rs) + v_s
Ra = [ln(zref/z0) - psi_h] / (kappa * u*)     [aerodynamic resistance]
Rb = Sc^(2/3) / (kappa * u*)                   [quasi-laminar sublayer]
Rs = 1 / (E0 * u*)                              [surface collection efficiency]
```

Deposition flux extracted at k=klo only:
```
dep_flux = v_d * rho_dust(k=klo)   [kg/m^2/s]
cc_source(k=klo, comp) -= dep_flux / dz
```

Accumulated total stored in `dust_deposition_rate` MultiFab [kg/m²] for PHREEQC feedback.

---

## Phase 13: Two-Way Coupling (ATM → Dust)

### File: `ERF_DustAtmReturn.H`

**`fill_dust_conc_from_atm(dust_conc_sfc, S_cons, dust_comp, geom_atm, C)`**
- Maps `RhoDust(k=klo)` from 3D conserved state to 2D dust grid
- Pattern: `conc(i,j,0) = max(S_cons(i/C, j/C, klo, dust_comp), 0)`

**`fill_dust_moist_from_atm(dust_surf_moist, Q1fx3, geom_atm, C)`**
- Maps vertical moisture flux at z-face klo to 2D dust grid
- Null-safe: when `Q1fx3 == nullptr` (no moisture scheme), zeroes `dust_surf_moist`

Loading feedback (Shao 2001):
```
u*_t *= (1 + loading_feedback_coeff * conc_sfc)
```

Fecan moisture inhibition:
```
excess = max(w - w_prime, 0)
u*_t *= sqrt(1 + a_f * excess)
```

---

## Phase 14: MRF Countergradient Extension to Dust Scalars

### Files: `ERF_ComputeDiffusivityMRF.cpp`, `ERF_TurbStruct.H`

`EddyDiff::Scalar_v` set equal to `EddyDiff::Theta_v` in all three MRF mixing branches
(K-profile unstable, stable inside PBL, free atmosphere above PBL).

Optional Schmidt number scaling:
```cpp
K_turb(i, j, k, EddyDiff::Scalar_v) *= Pr_t / dust_mrf_Sc_t;
```

PBLH write-back (post-Phase-9 fix in `ERF_ComputeDiffusivityMRF.cpp`):
After the MFIter loop, copy `Turb_lengthscale(k=klo)` into `SurfLayer->get_pblh(level)`.
This is the only location where MRF computes PBLH; without this, `dust_pblh` reads `1e+150`.

---

## Phase 16: Dust Diagnostics & Plotfile Output

### Files: `ERF_DustPlotfile.H/.cpp`, `ERF_DustStatsOutput.H`, `ERF_DustPlotfileCatalog.H`

**Plotfile format** (VisMF::Write + WriteGenericPlotfileHeader on native dust grid):
```
plt_dust_NNNNN/
  ├── Header           (AMReX header)
  └── Level_0/Cell     (VisMF binary, all dust MultiFab components)
```

**`dust_diag.dat`** columns: `step, time_s, emission_flux_max, emission_flux_sum, deposition_rate_sum, ustar_max, conc_sfc_max`

Write triggered by `dust_plot_int > 0` or `is_final = true`.

---

## Phase 17: EPA NAAQS PM2.5/PM10 Compliance

### Files: `ERF_DustPM.H`, `ERF_DustNAAQSOutput.H`

**PM size classification:**
- PM2.5: bins with `d <= 2.5e-6 m`
- PM10: bins with `d <= 10e-6 m`

**24-hour running average (exponential moving average):**
```
C_24h(t) = C_24h(t-dt) * (T-dt)/T + C_now * dt/T,   T = 86400 s
```

**NAAQS thresholds:**
- PM2.5 24h: 35 µg/m³ (40 CFR Part 50)
- PM10 24h: 150 µg/m³

Output to `dust_naaqs.csv`: `step, time_s, PM25_max, PM25_24h_max, PM10_max, PM10_24h_max, PM25_exceed_cells, PM10_exceed_cells`

---

## Phase 18: MSHA Worker Exposure Tracking

### Files: `ERF_DustMSHA.H`, `ERF_DustMSHAOutput.H`

**8-hour time-weighted average:**
```
dose(t+dt) = dose(t) + conc_sfc * dt   [mg/m^3 * s]
TWA = dose / elapsed_shift_time         [mg/m^3]
```

Shift reset when `floor(t / shift_duration) > floor((t-dt) / shift_duration)`.
At reset: `dust_msha_shift_twa` saved, `dust_msha_dose` zeroed, `m_msha_shift_count++`.

Output:
- `msha_exposure.csv`: per-step TWA, exceedance flag, dose accumulator
- `msha_shift_summary.csv`: one row per shift with end-of-shift TWA and peak

MSHA PEL: 5.0 mg/m³ respirable dust (30 CFR 56.5001).

---

## Phase 19: Lagrangian Super-Particle Tracking

### Files: `ERF_DustPC.H/.cpp`

**Particle container:** `ParticleContainer<0,0,5,0>` — 5 real SoA attributes: `{d_m, rho_p, dep_x, dep_y, source_id}`.

**`ReleaseParticles(emission_flux, geom_atm, geom_dust, dt, d_m, rho_p)`:**
- Uses `LoopOnCpu` over dust grid to release particles proportional to emission flux
- Each emission cell releases one particle with position at cell center at z=klo
- Does NOT use `ParallelFor` (see `DUST_MPI_SKILLS.md` Rule D1)

**`AdvanceParticles(xvel, yvel, zvel, source_map, geom_atm, geom_dust, dt)`:**
- Nearest-cell velocity interpolation (face → cell center average at particle cell)
- Stokes settling: `z -= v_s * dt`
- Deposition: when `z < z_lo + 0.5*dz_atm`, particle deposits and increments `dust_source_map`

**Limitation:** Nearest-cell interpolation. Full trilinear (as in `ERFPCEvolve.cpp`) is a future improvement.

---

## Phase 20: Multi-Mine Domain / Per-Site PHREEQC

### File: `ERF_DustSiteRegistry.H`

**`populate_dust_site_id(dust_site_id, geom_dust, site_x_lo, site_y_lo, site_x_hi, site_y_hi)`:**
- Assigns integer site index (1..N) to each dust cell based on bounding box containment
- Cells outside all bounding boxes → site_id = 0 (unassigned)

**`count_site_cells(dust_site_id, n_sites)`:**
- Returns `Vector<int>` of size n_sites+1: counts[0] = unassigned, counts[s] = cells in site s

Used by Phase 21 for per-site deposition budget.

---

## Phase 21: PHREEQC Deposition Feedback Writer

### File: `ERF_DustPHREEQCWriter.H`

**`write_phreeqc_deposition_file(filename, dust_deposition_rate, geom_dust, cur_time, nstep)`:**
- Writes full 2D grid of accumulated deposition [kg/m²] in a 64-row binary text format
- File is **overwritten** each interval (not appended)
- Triggered when `(cur_time - m_last_phreeqc_write_time) >= phreeqc_feedback_interval_s`

**`append_phreeqc_site_summary(filename, dust_deposition_rate, dust_site_id, geom_dust, cur_time, site_names)`:**
- Appends one row per interval to per-site CSV: `time_s, site_name, dep_total_kg_m2`

**Duplicate-write guard:** `m_last_phreeqc_write_step` prevents double-write when called from both
`WriteAtIntermediateTime` and `WriteAtFinalTime` in the same step.

---

## Phase 22: Haul Road Vehicle Schedule Emission

### File: `ERF_DustRoadSchedule.H`

**CSV format:**
```
road_name, x_lo, y_lo, x_hi, y_hi, road_width_m, vehicle_weight_t, silt_pct, vmt_per_h, start_time_s, end_time_s
haul_main, 1000, 800, 2800, 900, 10.0, 130, 8.5, 20, 0.0, -1
```

**EPA AP-42 Ch. 13.2.2 unpaved road formula:**
```
E_road = k * (s/12)^a * (W/3)^b   [lb/VMT]
E_road → F_road [kg/m^2/s] via unit conversion and cell area
```

Active when `t >= start_time_s AND (end_time_s < 0 OR t <= end_time_s)`.
Emission added to `dust_emission_flux` in cells overlapping the road bounding box.

**MPI pattern:** Rank-0 reads CSV → broadcast via `MPI_Bcast` of POD struct array (same as blast schedule).

---

## Phase 23: DOE Critical Materials Assessment

### File: `ERF_DustCriticalMaterials.H`

**`compute_cm_flux(dust_cm_flux, dust_emission_flux, cm_fractions, n_active)`:**
- Per-cell: `cm_flux(i,j,0) = sum(cm_fractions[b] * emission_flux(i,j,0,b), b=0..n_active-1)`
- When `transport_bins_separately=false`: n_active=1, fraction applied to summed bin

**`append_cm_budget(filename, dust_cm_flux, dust_site_id, geom_dust, cur_time, nstep, site_names)`:**
- One row per step: domain total + one row per site (including unassigned)
- Columns: `step, time_s, site_id, site_name, cm_flux_kg_m2_s, cm_flux_sum_kg_s`

Reference: US DOE (2023) Critical Materials Assessment. REE mass fractions in mine tailings: 0.001–0.01 kg\_CM/kg\_dust.

---

## Phase 24: Regression Test Suite

### `DustIntegration` End-to-End Test

**Location:** `Exec/CanonicalTests/Dust/DustIntegration/`

**Configuration:**
- 8×8×64 ATM cells, neutral ABL, MRF PBL, 15 m/s geostrophic wind
- `max_step=10`, `dt=0.5 s`
- All phases exercised simultaneously
- `grid_ratio=1`, `amr.max_level=0`

**All 14 pass criteria verified by the test:**
1. Exit code 0 at step 10
2. `emission_flux_max > 0` every step (wind-driven + road)
3. `cm_flux_max > 0` every step (`cm_fractions=0.001`)
4. `active_roads=1` from step 2 (`haul_main` active at `t >= 1.0s`)
5. Phase 22 road debug printed from step 2
6. Phase 23 `cm_flux_max ~ 0.001 × emission_max`
7. Phase 18 shift reset at t=5.0s and t=10.0s
8. `dust_diag.dat`: 10 rows, `emission_total > 0`
9. `dust_naaqs.csv`: 10 rows, `PM10_max > 0`
10. `msha_exposure.csv`: 10 rows, TWA increasing
11. `dust_dep_feedback.dat`: written at t=3, 6, 9, 10
12. `dust_road_diag.csv`: rows for `haul_main` from step 2
13. `dust_cm_budget.csv`: 10 groups (1 domain + 1 unassigned per step)
14. `PBLH_max ~ 500–1000 m` (not `1e+150`) after MRF PBLH fix

---

## Post-Merge Bug Fixes

### Phase 10 Grid-Ratio Bug (Fixed)

| Commit | Bug | Symptom | Fix |
|--------|-----|---------|-----|
| `grid_ratio_fix` | `MultiFab::Copy` with mismatched BoxArrays | Crash/assertion when `grid_ratio > 1` in `apply_to_cc_source` | Allocate `dust_bin_tmp` on dust grid; use `average_down` to coarsen to ATM slab. Remove `coarsen_dust_flux_to_atm(dust_flux_atm, dust_emission_flux, ...)` call. |

### Phase 14 PBLH Bug (Fixed)

| Commit | Bug | Symptom | Fix |
|--------|-----|---------|-----|
| `pblh_mrf_fix` | MRF never writes to `SurfLayer->pblh` | `PBLH_max = 1e+150` in dust Phase 9 debug | After `ComputeDiffusivityMRF` MFI loop, copy `Turb_lengthscale(k=klo)` into `SurfLayer->get_pblh(level)`. |

### Phase 13 Moisture Null Pointer (Fixed)

| Commit | Bug | Symptom | Fix |
|--------|-----|---------|-----|
| `q1fx3_null_fix` | `SFS_q1fx3_lev` is null when `moisture_type=None` | Null dereference in `fill_dust_moist_from_atm` | Null-safe guard in `ERF_DustAtmReturn.H`: if `Q1fx3 == nullptr`, zero `dust_surf_moist` and return |

### Phase 4 PHREEQC q_star Null Pointer (Fixed)

| Commit | Bug | Symptom | Fix |
|--------|-----|---------|-----|
| `qstar_null_fix` | `get_q_star()` returns null in dry runs | Null dereference in Pass 2 of `ComputeDiffusivityMRF` | Check `q_star_mf != nullptr` before constructing `Array4`; use empty `Array4<Real const>{}` as fallback |

---

## Build Rules (apply to every phase)

1. New `.cpp` file: add to **both** `CMake/BuildERFExe.cmake` dust block **AND** `Source/Dust/Make.package`.
2. Header-only file: add to `Make.package` CEXE_headers only.
3. Never `#include <AMReX_ParallelFor.H>` in any dust file. Use `<AMReX_MFIter.H>`.
4. Never wrap `.cpp` file bodies in `#ifdef ERF_USE_DUST`. Only wrap `#include` of dust headers in consumer files.
5. All dust MultiFabs use `amrex::IntVect(1,1,0)` ghost cells (1 in x,y; 0 in z — 2D slab).
6. `dust_flux_atm` is always allocated on the **ATM 2D slab** (not the dust grid). Only `dust_emission_flux` lives on the dust grid.

---

## Test Input Requirements (apply to every test)

All files under `Exec/CanonicalTests/Dust/*/inputs` must have:
- `amrex.fpe_trap_invalid = 0`
- `erf.dust.dust_debug = true`
- `amr.max_level = 0`
- `erf.use_gravity = true`
- `erf.most.zref = 24.0`
- `erf.most.z0 = 0.1`
- Neutral ABL configuration with `erf.pbl_type = "MRF"`
- Sounding file: `sounding_neutral_abl` (copied to each test directory)
- Geostrophic wind: `erf.abl_geo_wind = 15.0 0.0 0.0` (gives u\*~0.56 m/s from MRF)
- `erf.transport_scalar = true` for all tests with `erf.dust.atm_feedback >= 1.0`
- `erf.sum_interval = 1`

---

## Known Limitations

- Phase 15 (MRF scale-aware diffusivity blending at fine AMR levels) is not implemented.
- Phase 19 particle advection uses nearest-cell velocity interpolation; trilinear interpolation is a future improvement.
- PHREEQC coupling is loose (offline file exchange), not tight (runtime coupling).
- Road emission uses AP-42 Chapter 13.2.2 (unpaved roads only); paved road emission (Ch. 13.2.1) is not implemented.
- CM fractions are constant per bin; time-varying fractions from PHREEQC mineralogy updates are planned.
- All tests use `amr.max_level=0`; multi-level AMR dust transport is not validated.
- `transport_bins_separately=true` allocates per-bin scalar components but `extract_atm_return_fields` only reads bin 0 back; bins 1+ are injected but not returned to the dust layer.

---

## Parameter Recommendations

### Emission parameters
- `threshold_A_coeff = 0.0123`: Bagnold (1941) for quartz-dominated mine tailings.
- `silt_fraction = 0.05–0.20`: Typical mine tailings range (0.10 default).
- `rho_air = 1.225`: Standard sea-level density [kg/m³].

### Settling and deposition
- `bin_diameters = 7.0e-6 2.5e-6 50.0e-6`: PM2.5, PM10, coarse fractions.
- `deposition_E0 = 3.0e-3`: Bare mine surface (Zhang et al. 2001). Use 5e-3 for rough tailings; 1e-3 for paved.

### Suppression
- `supp_tau_base_s = 3600.0`: Water ~1h at 20°C low wind. Use 7200–14400 for MgCl₂. Use 1800 above 35°C.

### PHREEQC coupling
- `phreeqc_update_interval_s = 86400.0`: Daily geochemical updates appropriate for crust formation timescales.

### MSHA / NAAQS output
- `msha_pel_mg_m3 = 5.0`: 30 CFR 56.5001 respirable dust PEL.
- `msha_shift_duration_s = 28800.0`: Standard 8-hour mine shift.

### Critical materials
- `cm_fractions = 0.001–0.01`: Typical REE mass fractions in mine tailings (0.1–1%).

### Performance
- `transport_bins_separately = false`: Default; reduces conserved state array size.
- `atm_feedback = 0.0`: Standalone dust diagnostics without atmospheric coupling.
- `grid_ratio = 1`: Safe for all runs. `grid_ratio > 1` requires fixed `ERF_DustLayer.cpp::apply_to_cc_source`.