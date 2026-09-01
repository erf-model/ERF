# Radiation Development Reference

This document summarizes the implemented capabilities of ERF's TwoStream radiation pathway and maps those capabilities to runtime controls, numerical contracts, and regression coverage. It is intended to serve as a durable design reference rather than a phase-by-phase task log.

## Architectural Scope

The TwoStream implementation provides shortwave and longwave column radiation, optional cloud and aerosol optical-depth modifiers, coupling to thermodynamic and boundary-layer tendencies, runtime diagnostics controls, and simplified surface-energy-balance extensions. All additions are designed to preserve backward compatibility by default, maintain GPU safety within AMReX execution patterns, and remain verifiable through canonical regression cases under `Exec/CanonicalTests/Radiation/`.

## Capability Index

| Capability Area | Primary Focus | Representative `RadChoice` Controls | Associated RegTests |
|---|---|---|---|
| Base Two-Stream | Clear-sky SW/LW transport and diagnostics | `rad_type`, `sw_enabled`, `lw_enabled`, `solar_constant`, `solar_zenith`, `tau_*`, `isothermal_test` | `SW_ClearSky_Analytical`, `LW_Isothermal` |
| Cloud Optical Depth | Height-dependent cloud extinction | cloud optical-depth and cloud-layer controls | `SW_Cloud_Layer` |
| Scattering | Diffuse SW scattering in cloudy columns | scattering enable/optical-property controls | `SW_Scattering_Cloud` |
| RhoTheta Coupling | Heating-rate injection into model thermodynamics | `qheating_rates` pathway and thermodynamic source coupling | `TwoStream_RhoTheta_Coupling` |
| Time Integration and Diagnostics Cadence | Consistent repeated calling and call-site identity | diagnostic call-site and output controls | `TwoStream_TimeIntegration`, `TwoStream_Benchmark_Suite` |
| Surface Heterogeneity | Column-varying albedo, emissivity, and skin temperature with fallback | surface scalar fallback controls and LSM/Noah-MP field ingestion | `TwoStream_SurfaceHeterogeneity` |
| Dynamic Optical Depth | Moisture- and cloud-aware `tau(k)` diagnosis | dynamic optical-depth switches and fallback controls | `TwoStream_DynamicTau_MoistCloud`, `TwoStream_NonuniformDZ` |
| PBL Coupling | YSUNew radiative tendency coupling, limiting, and smoothing hooks | `enable_ysu_topdown`, `enable_ysu_rad_tend_limiter`, `ysu_rad_tend_*` | `TwoStream_PBL_MRF_YSU_Coupling` |
| Prognostic Cloud Fraction | RH/`qc`-based cloud-fraction diagnosis | prognostic cloud-fraction controls | `TwoStream_ProgCloudFraction` |
| Aerosol/Turbidity | Prescribed aerosol optical-depth contributions | aerosol/turbidity profile controls | `TwoStream_Aerosol_Turbidity` |
| Solar Geometry | Time-varying zenith-angle computation | `solar_geometry_dynamic_enable`, `latitude_deg`, `longitude_deg`, `day_of_year`, `time_zone_offset_hours` | `TwoStream_DiurnalSolarGeometry` |
| Simplified SEB | Surface-flux infrastructure, diagnostics, prognostic update, and safeguards | `seb_enable`, `seb_diagnostic_enable`, `seb_prognostic_enable`, SEB fallback and clamp controls | `TwoStream_SEB_MultiFabInfra`, `TwoStream_SEB_Diagnostic`, `TwoStream_SEB_Prognostic` |

## Base Two-Stream Solver

### Governing Equations

For direct-beam shortwave verification, the implementation retains the Beer-Lambert attenuation law:

```math
F_{\mathrm{dir}}(z) = S_0 \cos(\theta_z)\exp\left(-\tau_{\mathrm{cum}}(z)/\cos(\theta_z)\right)
```

For longwave verification in the isothermal limit, the consistency target is

```math
F_\uparrow = F_\downarrow = \sigma T_{\mathrm{iso}}^4,\qquad
F_{\mathrm{net}} = 0,\qquad
\frac{\partial T}{\partial t}=0.
```

The operational solver follows the standard two-stream column sweep described by Toon et al. (1989), with horizontal parallelism and sequential vertical accumulation per column.

### Representative `RadChoice` Parameters

- `rad_type = TwoStream`
- `sw_enabled`, `lw_enabled`
- `solar_constant`, `solar_zenith`
- shortwave and longwave optical-depth controls
- `isothermal_test`
- diagnostic file and verbosity controls

### Backward-Compatibility Notes

- The clear-sky path remains the default baseline for isolated SW/LW verification.
- Disabling optional modifiers such as clouds, scattering, aerosols, dynamic solar geometry, and SEB extensions preserves the original minimal TwoStream behavior.
- Longwave isothermal mode remains a strict self-consistency check and does not alter production pathways unless explicitly enabled.

### GPU-Safety Notes

- Helper routines intended for device execution must remain `AMREX_GPU_DEVICE` or `AMREX_GPU_HOST_DEVICE` and force-inlined when used inside kernels.
- Host-side I/O is excluded from device lambdas; diagnostics are reduced on device and printed or written on host.
- Column sweeps use one thread per horizontal column with a sequential `k` loop inside the device lambda.

### Associated RegTests

- `SW_ClearSky_Analytical`
- `LW_Isothermal`

## Cloud Optical Depth

### Governing Equations

Cloud optical-depth capability generalizes extinction from a uniform layer value to a height-dependent contribution:

```math
\tau(k) = \tau_{\mathrm{background}}(k) + \tau_{\mathrm{cloud}}(k),
```

where `\tau_cloud(k)` is activated only over the prescribed cloud layer or diagnosed cloudy region.

### Representative `RadChoice` Parameters

- cloud optical-depth magnitude controls
- cloud-layer bounds and related profile-type controls
- cloud-fraction masking inputs

### Backward-Compatibility Notes

- When cloud optical-depth modifiers are disabled, the solver reverts to the original clear-sky profile.
- The initial cloud-layer implementation preserved the uniform optical-depth pathway for cases that did not request cloud structure.

### GPU-Safety Notes

- Cloud contributions are evaluated locally within the per-column sweep.
- No host-managed per-column state is introduced; the added optical-depth term is resolved entirely from device-available parameters or fields.

### Associated RegTests

- `SW_Cloud_Layer`

## Scattering

### Governing Equations

Shortwave scattering is represented with the two-stream approximation described by Meador and Weaver (1980), augmenting absorption-only transport with diffuse redistribution through cloud optical properties such as single-scattering albedo and asymmetry assumptions.

### Representative `RadChoice` Parameters

- scattering enable/disable control
- cloud optical-depth and scattering-property inputs
- existing SW flux and diagnostics controls

### Backward-Compatibility Notes

- The scattering path is additive. When disabled, the absorption-only cloud treatment remains unchanged.
- Existing clear-sky and pure-absorption cases continue to use the same flux evaluation pathway.

### GPU-Safety Notes

- The scattering approximation is incorporated into the same device-safe vertical sweep structure used by the base SW solver.
- Numerical guards remain necessary for extreme optical properties to prevent non-finite fluxes from propagating.

### Associated RegTests

- `SW_Scattering_Cloud`

## RhoTheta Coupling

### Governing Equations

Radiative heating is coupled into thermodynamics through the heating-rate tendency source:

```math
\frac{\partial (\rho\theta)}{\partial t}\Big|_{\mathrm{rad}}
= \rho\,Q_{\mathrm{rad}},
```

with separate shortwave and longwave contributions first written to `qheating_rates` and then injected into the thermodynamic update path.

### Representative `RadChoice` Parameters

- TwoStream activation controls
- heating-rate storage through `qheating_rates`
- runtime options that determine whether SW and LW source terms are populated

### Backward-Compatibility Notes

- The coupling is designed so that cases not consuming `qheating_rates` retain the prior uncoupled behavior.
- The diagnostic-only simplifications used before thermodynamic injection are superseded by the coupled pathway when heating tendencies are enabled.

### GPU-Safety Notes

- Heating-rate computation and storage are cell-local and compatible with tiled AMReX kernels.
- No host accumulation is permitted inside the per-cell source update.

### Associated RegTests

- `TwoStream_RhoTheta_Coupling`

## Time Integration and Diagnostics Cadence

### Governing Equations

The critical contract is temporal consistency rather than a new physical closure: radiation diagnostics and source application must respect timestep cadence and call-site identity. Diagnostic identity is therefore treated as a tuple of

```text
(step, time, call_site)
```

rather than `step` alone.

### Representative `RadChoice` Parameters

- diagnostic enable/disable controls
- diagnostic output filename and schema controls
- call-site mode controls such as `both`, `pre_only`, and `post_only`

### Backward-Compatibility Notes

- Default diagnostics behavior preserves the established CSV layout when optional cadence extensions are not enabled.
- Filtering by call site must not alter the underlying physics path; it only changes what is recorded.

### GPU-Safety Notes

- Diagnostics collection relies on device-side reductions followed by host-side file writes.
- Time and call-site metadata are assembled outside device kernels.

### Associated RegTests

- `TwoStream_TimeIntegration`
- `TwoStream_Benchmark_Suite`

## Surface Heterogeneity

### Governing Equations

Surface-coupled shortwave and longwave boundary conditions depend on column-varying albedo, emissivity, and skin temperature. The surface-property contract is a precedence rule:

1. use heterogeneous LSM or radiation-interface fields when available and finite;
2. otherwise use configured scalar fallback values;
3. otherwise use documented hard-safe defaults with warn-once behavior.

### Representative `RadChoice` Parameters

- scalar fallback albedo controls
- scalar fallback emissivity controls
- scalar fallback surface-temperature controls

### Backward-Compatibility Notes

- Scalar fallback resolution preserves legacy single-value behavior when heterogeneous fields are absent.
- Missing LSM fields must not terminate the run if safe fallback values are available.

### GPU-Safety Notes

- Surface-property values are clamped to physically meaningful bounds, including `[0,1]` for albedo and emissivity.
- Fallback resolution is handled in GPU-compatible field-fill and device-read patterns without host-side branching inside kernels.

### Associated RegTests

- `TwoStream_SurfaceHeterogeneity`

## Dynamic Optical Depth

### Governing Equations

Dynamic optical-depth diagnosis replaces fixed `\tau` with moisture- and cloud-aware per-level values:

```math
\tau(k) = f\!\left(q_v(k), q_c(k), \rho(k), dz(k)\right),
```

where the exact implementation combines vapor, condensate, density, and layer-thickness information while retaining a safe fallback path to prescribed optical depth.

### Representative `RadChoice` Parameters

- dynamic optical-depth enable/disable control
- fallback `tau` controls for SW and LW
- any moisture/cloud optical conversion coefficients exposed through `RadChoice`

### Backward-Compatibility Notes

- Disabling dynamic optical depth restores the prescribed optical-depth path.
- The nonuniform-`dz(k)` wiring retains a uniform fallback so older configurations remain valid.

### GPU-Safety Notes

- `dz(k)` and moisture fields are consumed directly in device kernels.
- Finite guards are required because bad thermodynamic state values can otherwise produce non-physical extinction coefficients.

### Associated RegTests

- `TwoStream_DynamicTau_MoistCloud`
- `TwoStream_NonuniformDZ`

## PBL Coupling

### Governing Equations

The YSUNew coupling pathway consumes radiative tendencies in the boundary-layer mixing logic and optionally limits them:

```math
Q_{\mathrm{rad,limited}} = \min\!\left(\max\!\left(Q_{\mathrm{rad}},-Q_{\max}\right),Q_{\max}\right),
```

with optional smoothing strength reserved for controlled damping of extreme radiative forcing before it influences top-down mixing.

### Representative `RadChoice` and Related Parameters

- `enable_ysu_topdown`
- `enable_ysu_rad_tend_limiter`
- `ysu_rad_tend_limiter_magnitude`
- `ysu_rad_tend_smooth_strength`

### Backward-Compatibility Notes

- Limiting and smoothing are disabled by default.
- When disabled, the YSUNew path is intended to remain numerically identical to the prior uncapped coupling behavior.
- The documented implementation scope is YSUNew-focused; MRF-specific extension remains separate.

### GPU-Safety Notes

- Finite checks use device-safe math such as `std::isfinite()`.
- Limiter logic is cell-local and introduces no cross-thread synchronization.

### Associated RegTests

- `TwoStream_PBL_MRF_YSU_Coupling`

## Prognostic Cloud Fraction

### Governing Equations

Cloud fraction is diagnosed from relative humidity and condensate content, then bounded:

```math
C_f = \mathrm{clip}\!\left(g(RH,q_c),\,0,\,1\right).
```

The formulation is designed to provide a radiation-facing cloud mask without breaking the pre-existing prescribed or implicit-cloud pathways.

### Representative `RadChoice` Parameters

- prognostic cloud-fraction enable/disable control
- RH and `qc` threshold or smoothing controls

### Backward-Compatibility Notes

- Feature-off behavior preserves prior radiation treatment of cloud fraction.
- Temporal smoothing and clipping are used to prevent abrupt changes from destabilizing existing regression baselines.

### GPU-Safety Notes

- RH and condensate diagnostics are evaluated per cell with finite guards.
- Cloud-fraction clipping to `[0,1]` is mandatory before use in optical calculations.

### Associated RegTests

- `TwoStream_ProgCloudFraction`

## Aerosol and Turbidity

### Governing Equations

Prescribed aerosol loading contributes an additional optical-depth term:

```math
\tau_{\mathrm{total}}(k) = \tau_{\mathrm{base}}(k) + \tau_{\mathrm{aerosol}}(k),
```

where the aerosol component may be represented as a constant, exponential profile, or tabulated structure.

### Representative `RadChoice` Parameters

- aerosol/turbidity enable control
- aerosol profile-type selection
- aerosol optical-depth magnitude or table controls
- optional LW aerosol hook parameters, when enabled

### Backward-Compatibility Notes

- The default path excludes aerosol contributions.
- Existing cloud and clear-sky tests remain unchanged unless the aerosol option is explicitly requested.

### GPU-Safety Notes

- Profile evaluation must remain device-safe whether the aerosol profile is constant, analytic, or table-driven.
- Bounds and finiteness checks are required to avoid negative or undefined optical depths.

### Associated RegTests

- `TwoStream_Aerosol_Turbidity`

## Solar Geometry

### Governing Equations

Dynamic solar geometry computes the solar zenith angle from astronomical relationships among latitude, longitude, day of year, and local solar time. The shortwave top-of-atmosphere forcing then follows

```math
F_{\mathrm{TOA}} = S_0 \max\!\left(\cos\theta_z, 0\right).
```

The implementation includes declination, equation-of-time, and hour-angle calculations and corrects the longitude adjustment using the deviation from the local standard meridian.

### Representative `RadChoice` Parameters

- `solar_geometry_dynamic_enable`
- `latitude_deg`
- `longitude_deg`
- `day_of_year`
- `time_zone_offset_hours`
- fallback `solar_zenith`

### Backward-Compatibility Notes

- Dynamic geometry is disabled by default.
- When disabled, `cos(zenith)` is evaluated exactly from the fixed `solar_zenith` scalar so existing baselines remain unchanged.

### GPU-Safety Notes

- Astronomical helper functions are implemented as inline device-safe math.
- Range checks on latitude, longitude, and day-of-year inputs are required before trigonometric evaluation.
- Time-of-day reduction uses GPU-safe modulo arithmetic.

### Associated RegTests

- `TwoStream_DiurnalSolarGeometry`

## Simplified Surface Energy Balance

### Governing Equations

The diagnostic SEB residual is

```math
\mathrm{SEB}_{\mathrm{res}} = SW_{\mathrm{net}} + LW_{\mathrm{net}} - H - LE - G.
```

The prognostic surface-temperature update uses a force-restore form:

```math
C_s\frac{dT_s}{dt} = \mathrm{SEB}_{\mathrm{res}} - C_s\frac{2\pi}{\tau}(T_s-T_{\mathrm{deep}}),
```

and the moisture update follows

```math
\frac{dq_s}{dt} = -\frac{LE}{L_v\rho_w d_s} - \frac{1}{\tau_q}(q_s-q_{\mathrm{deep}}).
```

### Representative `RadChoice` Parameters

- `seb_enable`
- `seb_diagnostic_enable`
- `seb_prognostic_enable`
- `seb_hfx_default`, `seb_lh_default`, `seb_grdflx_default`
- `seb_q_sfc_default`, `seb_t_deep_default`, `seb_q_deep_default`
- `seb_surface_heat_capacity`
- `seb_restore_timescale_s`
- `seb_moisture_layer_depth_m`
- `seb_moisture_restore_timescale_s`
- `seb_prognostic_t_min_k`, `seb_prognostic_t_max_k`
- `seb_prognostic_q_min`, `seb_prognostic_q_max`

### Backward-Compatibility Notes

- `seb_enable = false` avoids allocation of new SEB MultiFabs and preserves the pre-SEB baseline.
- `seb_diagnostic_enable = false` suppresses residual computation and additional CSV columns.
- `seb_prognostic_enable = false` preserves the diagnostic-only pathway.
- When Noah-MP is active, prognostic ownership rules prevent TwoStream SEB updates from overwriting Noah-MP-driven state.

### GPU-Safety Notes

- SEB helper functions are inline device routines with finite guards.
- Reductions for residual and summary statistics use AMReX `ReduceOps`.
- Prognostic updates are executed in GPU-safe per-column loops, with all clamping performed before values are written back.

### Associated RegTests

- `TwoStream_SEB_MultiFabInfra`
- `TwoStream_SEB_Diagnostic`
- `TwoStream_SEB_Prognostic`

## Known Issues and Resolutions

### Hardcoded Vertical Grid Bounds

Early prototypes used fixed vertical extents rather than deriving limits from `Box` or geometry metadata. That approach risks silent out-of-bounds access whenever the number of vertical levels changes. The implementation was corrected by deriving loop limits from `bx.smallEnd(2)`, `bx.bigEnd(2)`, or `geom[lev].Domain()`, and that rule should be treated as a permanent design requirement.

### Host-Side Invocation of Device Kernels

A transitional implementation attempted to call the device-oriented vertical sweep from host-side nested loops and then accumulate results on host. That pattern is invalid on GPUs and unsafe under threaded CPU execution. The resolution was to launch the work through `amrex::ParallelFor`, perform reductions on device, and move only reduced scalars back to host.

### Heating-Rate Coupling Requires More Than Correct Local Physics

One recurring lesson from the thermodynamic coupling work is that physically correct heating-rate calculations are insufficient unless they are also wired into the actual model tendency path. The permanent resolution was to verify both numerical correctness and downstream consumption, particularly through `qheating_rates` and thermodynamic source application.

### Diagnostics Identity Must Include Call Site

When radiation is invoked from multiple points in the timestep, treating a diagnostic record as uniquely identified by `step` alone creates duplicate-row ambiguity and can mask missing or repeated forcing. The durable resolution was to treat `(step, time, call_site)` as the record identity and to align validation scripts with that contract.

### Surface-Property Fallback Is a Physics Contract

Surface heterogeneity work exposed the need for a strict precedence rule among LSM-provided fields, scalar runtime fallbacks, and hard-safe defaults. The implementation now treats this as a documented contract so that missing land-surface fields do not abort runs or silently bypass valid heterogeneous inputs.

### SEB Sounding Input Quality

An early simplified-SEB regression input used an unrealistically dry moisture profile, which reduced the value of the case as a verification reference. The corrective action replaced the placeholder profile with a vertically varying moisture profile adapted from the moist-cloud regression family and documented realistic ranges in the test-specific README.

### Prognostic SEB Fill Ownership

In the initial prognostic SEB implementation, fallback field fills reinitialized `t_sfc` and `q_sfc` every radiation call, preventing the prognostic update from accumulating state. The resolution established explicit ownership: once prognostic mode is enabled, the prognostic update path owns those fields, while `t_deep` and `q_deep` remain passive restore targets.

### Radiation Heating Availability in Plotfile Capability Checks

Another confirmed issue was that plotfile capability checks recognized only the non-TwoStream radiation path when deciding whether `qsrc_sw` and `qsrc_lw` were available. The fix aligned the capability logic with the actual `qheating_rates` allocation criteria so that TwoStream runs correctly advertise radiation-heating output.

## References

- Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam, 1989: "Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres." *Journal of Geophysical Research*, 94, 16387–16405.
- Beer, A., 1852: "Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten." *Annalen der Physik und Chemie*, 86, 78–88.
- Meador, W. E., and W. R. Weaver, 1980: "Two-stream approximations to radiative transfer in planetary atmospheres: A unified description of existing methods and a new improvement." *Journal of the Atmospheric Sciences*.
- Spencer, J. W., 1971: "Fourier series representation of the position of the sun." *Search*, 2(5), 172–172.
- Duffie, J. A., and W. A. Beckman, 1991: *Solar Engineering of Thermal Processes*. John Wiley & Sons.
- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
