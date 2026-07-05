# ERF-native SHOC

This directory contains the ERF-native implementation of the Simplified
Higher-Order Closure (SHOC).

SHOC represents unresolved turbulence, shallow convection, and subgrid-scale
cloud macrophysics. It uses prognostic turbulent kinetic energy, diagnostic
higher-order moments, and an assumed probability density function (PDF) for
subgrid thermodynamic variability.

## Provenance

This implementation is based on the E3SM/EAMxx SHOC algorithms and source
structure. It is not a bit-for-bit port of EAMxx SHOC.

The native ERF implementation has been adapted to ERF data structures, AMReX
loops, ERF surface-flux coupling, diagnostics, and build systems.

Primary reference:

```text
Bogenschutz, P. A., and S. K. Krueger, 2013:
A simplified PDF parameterization of subgrid-scale clouds and turbulence
for cloud-resolving models.
Journal of Advances in Modeling Earth Systems, 5, 195-211.
https://doi.org/10.1002/jame.20018
```

Upstream project reference:

```text
E3SM Project
https://github.com/E3SM-Project/E3SM
```

Reference E3SM commit used during this port:

```text
0053a696ce40fc5ffa0b07a007426940b8af09a1
```

Keep the required license notices when moving, copying, or modifying code derived
from E3SM.

## Runtime selection

Native SHOC is always built in tree. It does not require EAMxx, EKAT, or Kokkos.

Select it at runtime:

```text
zlo.type = "surface_layer"

erf.pbl_type = NATIVE_SHOC
```

ERF also supports the optional EAMxx SHOC interface:

```text
erf.pbl_type = EAMXX_SHOC
```

The legacy runtime name:

```text
erf.pbl_type = SHOC
```

is a deprecated alias for `EAMXX_SHOC`. It does not select native SHOC.

Do not add new build workflows that use `ERF_ENABLE_NATIVE_SHOC` or
`USE_NATIVE_SHOC`. Native SHOC is an in-tree ERF implementation. Runtime
selection controls whether it runs.

## Surface fluxes

SHOC needs lower-boundary heat, moisture, and momentum fluxes.

Native SHOC does not compute those exchanges directly from a land surface model.
ERF computes the lower-boundary fluxes first and passes the resulting flux arrays
to SHOC.

Those fluxes come through ERF's surface-layer infrastructure. They may come from
Monin-Obukhov similarity theory (MOST), prescribed surface-layer inputs, or an
active land or ocean surface model when that model provides fluxes.

Native SHOC consumes the ERF flux arrays and converts ERF's host
density-weighted fluxes to kinematic surface fluxes.

## Microphysics coupling

SHOC diagnoses subgrid non-precipitating cloud partitioning with its assumed PDF.
This includes cloud fraction and non-precipitating liquid water.

Native SHOC remains a liquid-cloud macrophysics closure under the interim ice
contract. It may use pre-existing cloud ice for phase-aware thermodynamics and
buoyancy, but it does not create or repartition cloud ice.

To avoid double counting, ERF disables the saturation-adjustment or condensation
step in the microphysics package when a SHOC-family PBL scheme is active.

This does not disable microphysics. Microphysics still handles precipitating
processes outside SHOC's cloud macrophysics role. Number-aware microphysics
layouts with cloud-droplet or ice number concentrations still need an explicit
number closure in their own microphysics pathways if they are coupled to SHOC.
Native SHOC `state_update` rejects those number-aware layouts until that
number closure exists.

## Transport modes

Native SHOC has one scalar/cloud/TKE transport selector and one independent
horizontal-momentum selector.

The scalar/cloud/TKE transport mode is:

```text
erf.shoc.transport_mode = state_update
```

This is the default. In this mode, SHOC applies its coupled heat, moisture,
non-precipitating cloud liquid, carried cloud ice, and TKE update before the
dycore sees the state. That keeps the thermodynamic and moisture increments
together.

The host-diffusion scalar mode is:

```text
erf.shoc.transport_mode = host_diffusion
```

In this mode, SHOC exports the full vertical eddy-diffusivity block to ERF's
host diffusion path. SHOC does not apply the pre-dycore state update. At
present this mode is only supported for dry/no-moisture configurations, and it
requires `erf.shoc.momentum_transport = host_diffusion`.

The horizontal-momentum transport mode is:

```text
erf.shoc.momentum_transport = host_diffusion
```

This is the default. In the mixed native SHOC mode, SHOC exports only the
momentum diffusivity to ERF's host diffusion path while keeping the scalar
transport on the `state_update` path.

Other momentum choices are:

```text
erf.shoc.momentum_transport = none
erf.shoc.momentum_transport = state_update
```

`none` disables SHOC momentum transport entirely. `state_update` keeps the
legacy direct-velocity update available for debugging or targeted experiments,
but it is not the default.

## Preferred baseline input settings

Use `state_update` plus `momentum_transport = host_diffusion` as the baseline
native SHOC coupling mode. In this mode, native SHOC runs its column physics
before the ERF dycore. SHOC performs its implicit vertical turbulence solve,
reconstructs one coupled heat, moisture, cloud, and TKE update, and applies
that update to the old-time ERF state before acoustic and Runge-Kutta dycore
integration starts. Horizontal momentum stays on the host diffusion path.

A minimal preferred native SHOC block is:

```text
# Native SHOC selection.
zlo.type = "surface_layer"
erf.pbl_type = NATIVE_SHOC

# Native SHOC coupling.
erf.shoc.transport_mode = state_update
erf.shoc.momentum_transport = host_diffusion

# Preferred TKE-production behavior.
erf.shoc.signed_tke_production = true

# Useful during development and case bring-up.
erf.shoc.extra_shoc_diags = true
```

`state_update` is the preferred coupled mode for moist native SHOC runs.
Pair it with `momentum_transport = host_diffusion` unless you are explicitly
debugging the direct momentum update. This keeps SHOC's thermodynamic,
moisture, cloud, and TKE increments together before the dycore sees the state
without forcing a pre-dycore momentum update.

`erf.shoc.signed_tke_production = true` keeps the buoyancy contribution signed in
the TKE production term. Use it for the baseline unless a case has a specific
reason to use the clipped behavior.

Do not use `erf.shoc.transport_mode = tendencies` for native SHOC. That legacy
mode has been removed and should be rejected by the runtime parser.

`host_diffusion` remains available for dry/no-moisture configurations. In that
mode, SHOC exports vertical eddy diffusivities to ERF's host diffusion path, but
it does not apply the coupled pre-dycore state update. Do not use
`host_diffusion` for moist SHOC runs until microphysics and cloud-macrophysics
ownership is made transport-mode-aware.

The snippet above does not choose a microphysics package, radiation scheme, land
surface model, terrain option, grid, or timestep. Those remain case choices.

## Recommended native SHOC plotfile variables

Request SHOC diagnostics through the normal ERF plotfile variable lists. Use the
plotfile stream that matches the case. For example, use `erf.plot_vars_2` instead
of `erf.plot_vars_1` when the SHOC diagnostics should go to the second plotfile
stream.

A compact bring-up list is:

```text
erf.plot_vars_1 = density rhotheta theta qv qc Kmv Khv Lturb shoc_cldfrac shoc_ql shoc_ql2 shoc_cond wqls_sec wthv_sec w_sec thl_sec qw_sec qwthl_sec wthl_sec wqw_sec w3 brunt isotropy shear_prod buoy_prod diss_tke
```

The recommended variables cover these groups:

| Variable | Purpose |
| --- | --- |
| `density`, `rhotheta`, `theta`, `qv`, `qc` | State context for checking the SHOC thermodynamic and moisture update. |
| `Kmv`, `Khv`, `Lturb` | SHOC vertical momentum diffusivity, heat diffusivity, and turbulent length scale. |
| `shoc_cldfrac`, `shoc_ql`, `shoc_ql2`, `shoc_cond` | Assumed-PDF cloud fraction, liquid water, liquid-water variance, and condensation diagnostic. |
| `wqls_sec`, `wthv_sec`, `w_sec` | Vertical flux and vertical-velocity variance diagnostics. |
| `thl_sec`, `qw_sec`, `qwthl_sec`, `wthl_sec`, `wqw_sec`, `w3` | Higher-order thermodynamic, moisture, and velocity moments. |
| `brunt`, `isotropy`, `shear_prod`, `buoy_prod`, `diss_tke` | Stability, isotropy, and TKE-budget diagnostics. |

When combining native SHOC with host LES diagnostics, also consider adding
`Kmh`, `Khh`, and `nut`. `Kmh` and `Khh` come from the host horizontal eddy
diffusivity fields. `nut` is the cell-centered eddy viscosity diagnostic derived
from the momentum diffusivity.

Native-only SHOC diagnostics are meaningful only after the native driver has run.
If the user requests those variables in a non-native-SHOC run, the plotfile path
may write missing-value placeholders.

`erf.shoc.extra_shoc_diags` helps during development, but it is not the plotfile
selection mechanism. Always request plotfile fields explicitly with
`erf.plot_vars_1` or `erf.plot_vars_2`.

## Native SHOC coupling invariants

Preserve these invariants unless the design changes deliberately:

1. Select native SHOC with the runtime PBL type. Do not use compile-time macros
   as a substitute for runtime checks.
2. Native SHOC requires full-height AMReX boxes on SHOC-active levels. Do not
   use a vertical box split on those levels.
3. `state_update` is the default native SHOC transport mode.
4. `erf.shoc.momentum_transport = host_diffusion` is the default momentum mode.
5. In `state_update`, SHOC applies one coupled pre-dycore update to heat,
   moisture, non-precipitating cloud liquid, carried cloud ice, and TKE.
6. The `state_update` scalar update does not change density, vertical velocity,
   passive scalars, or precipitating species.
7. In mixed mode, SHOC consumes lower-boundary heat and moisture flux arrays,
   while host diffusion still owns the momentum stresses.
8. `erf.shoc.momentum_transport = state_update` is only for explicit debugging
   or development. When selected, SHOC also consumes the lower-boundary
   momentum stresses and performs the direct velocity update.
9. In `state_update`, SHOC must not also add fast or slow RHS tendencies.
10. After native SHOC updates the old-time state, ERF must synchronize boundary
   and halo data before pre-dycore checks, strain, turbulent viscosity, momentum
   conversion, primitive variables, pressure, or fast coefficients read the
   fields.
11. `host_diffusion` is dry/no-moisture only until microphysics and cloud
   macrophysics ownership is transport-mode-aware.
12. Number-aware microphysics layouts with cloud or ice number concentrations
    must abort until a number closure exists.
13. Microphysics remains active after SHOC. SHOC owns only non-precipitating
    liquid-cloud macrophysics under the interim cloud contract. Microphysics owns
    precipitation formation, precipitation sinks and sources, sedimentation,
    deposition, sublimation, freezing, melting, and other phase-change processes
    outside that contract.

The interim cloud contract means native SHOC may read a pre-existing cloud-ice
field for thermodynamics and buoyancy. It may carry that bounded ice field
through the coupled reconstruction. It does not create, destroy, or repartition
cloud ice as a microphysics process.

## ERF time-step integration

Native SHOC runs inside `ERF::Advance()` once per ERF time step on each
SHOC-active level. The broad order is:

```text
old/new state swap
-> fill-patch old state and velocities when needed
-> convert old velocities to old momenta
-> update surface-layer state and lower-boundary fluxes
-> update radiation sources on the old state
-> impose surface-layer boundary conditions for SHOC flux arrays
-> run native SHOC
-> synchronize old-state boundaries, halos, and momenta after SHOC state_update
-> run pre-dycore checks
-> run the dycore and multirate integrator
-> run microphysics, unless tight coupling moves it into the slow post-RHS path
-> run the land-surface model
-> impose new-state boundary conditions needed for later fill-patching
```

At the start of the step, ERF swaps old and new state pointers. On refined
levels, ERF fill-patches the old conserved state and velocities. ERF then
rebuilds old momenta from the old velocities and density.

If the lower boundary uses the surface-layer boundary condition, ERF updates the
surface-layer state and computes lower-boundary heat, moisture, and momentum
fluxes before SHOC runs. These fluxes may come from Monin-Obukhov similarity
theory (MOST), prescribed surface-layer inputs, or an active land or ocean
surface model when that model provides fluxes. ERF also updates radiation sources
on the old state before SHOC runs.

When `erf.pbl_type = NATIVE_SHOC`, ERF calls `compute_native_shoc_tendencies()`
before `advance_dycore()`. The coupling call passes the old conserved state,
face velocities, vertical velocity, surface stress and flux arrays,
eddy-diffusivity storage, terrain metrics, optional subsidence data, and `dt` to
`ShocDriver::advance()`.

`ShocDriver::advance()` preprocesses each full-height AMReX box into SHOC column
data. It reads density, potential temperature, moisture species, pressure and
temperature diagnostics, velocities, geometry, carried turbulence state, and
lower-boundary fluxes. It converts ERF's density-weighted surface fluxes and
stresses to the kinematic quantities used by SHOC.

The coupling then follows `erf.shoc.transport_mode` and
`erf.shoc.momentum_transport`.

In `state_update`, SHOC caches the baseline state, diagnoses structure and TKE,
runs the implicit column update, diagnoses moments and the assumed-PDF cloud
state, reconstructs the updated ERF thermodynamic state, and forms tendencies
relative to the baseline state. If `erf.shoc.momentum_transport = state_update`,
it also applies those tendencies directly to the horizontal velocities. The
scalar update remains a state update, not an ERF RHS tendency path.

After the `state_update` increment, `ERF::Advance()` immediately refreshes the
fields the dycore will read. On level 0, `FillPatchCrseLevel()` refills `S_old`,
`U_old`, `V_old`, and `W_old`, and `VelocityToMomentum()` rebuilds `rU_old`,
`rV_old`, and `rW_old` from the synchronized velocities and density. On refined
levels, `FillPatchFineLevel()` refills conserved state, velocities, and momenta
through the coarse/fine fill-patch path. This synchronization happens before the
pre-dycore NaN and temperature checks and before `advance_dycore()` computes
strain, eddy viscosity, primitive variables, Exner pressure, fast coefficients,
and the multirate time integrator.

During the dycore step, native SHOC in `state_update` mode does not add fast or
slow RHS source terms. `ShocDriver::set_eddy_diffs()` clears the SHOC-owned
vertical diffusivity components from the host eddy-diffusivity arrays, except for
non-transport diagnostics such as the SHOC length scale. In mixed mode it
exports only the momentum diffusivity to the host diffusion path. This prevents
the host from reapplying SHOC scalar vertical transport.

`ShocDriver::set_diff_stresses()` clears the lower-boundary scalar fluxes that
SHOC already consumed in `state_update` mode. It clears the momentum stresses as
well only when native SHOC is explicitly configured to own momentum
state-update transport.

In `host_diffusion`, SHOC does not apply the pre-dycore state update and does not
own the surface fluxes. Instead, SHOC exports vertical eddy-diffusivity
coefficients to ERF's host diffusion path. This mode is currently limited to
dry/no-moisture configurations because SHOC-family microphysics suppresses
saturation adjustment or condensation while host-diffusion SHOC does not own the
cloud-macrophysics mass update.

Microphysics runs after the dycore with the default loose coupling. In that mode,
`ERF::Advance()` calls `advance_microphysics()` after `advance_dycore()` produces
`S_new`. If `erf.moisture_tight_coupling = true`, ERF calls microphysics inside
the slow post-RHS path instead. In both cases, SHOC-family PBL selection disables
the microphysics saturation-adjustment or condensation step to avoid double
counting SHOC's non-precipitating liquid-cloud macrophysics. Microphysics still
handles precipitating processes outside SHOC's cloud macrophysics role.

After loose-coupled microphysics, ERF advances the land-surface model. That
land-surface update affects fluxes used by later surface-layer updates. It is not
a second SHOC surface-flux application in the same pre-dycore state update.

## Runtime options

Native SHOC reads options from the `erf.shoc` namespace.

The source of truth for option names, defaults, and validation rules is:

```text
ERF_ShocTypes.H
```

The user-facing documentation is:

```text
Docs/sphinx_doc/theory/PBLschemes.rst
```

Do not duplicate the full options table here. Keep the Sphinx documentation as
the user manual.

## Code organization

The main source files are:

```text
ERF_ShocDriver.cpp       Driver and ERF-facing orchestration
ERF_ShocPreprocess.cpp   ERF state and flux arrays to SHOC column data
ERF_ShocStructure.cpp    SHOC structure diagnostics
ERF_ShocTKE.cpp          TKE and eddy diffusivity diagnostics
ERF_ShocMoments.cpp      Second- and third-moment diagnostics
ERF_ShocPDF.cpp          Assumed-PDF cloud diagnostics
ERF_ShocEnergyFixer.cpp  Energy and consistency fixes
ERF_ShocImplicit.cpp     Implicit state-update reconstruction
ERF_ShocDiagnostics.cpp  Diagnostic sequencing
ERF_ShocCoupling.cpp     ERF coupling helpers
```

The main headers are:

```text
ERF_ShocTypes.H
ERF_ShocColumnData.H
ERF_ShocDriver.H
```

## Developer notes

Like ERF's other column physics, it requires each AMReX box on a SHOC-active
level to span the full vertical domain. Do not use a grid decomposition that
splits boxes in the vertical direction. If an input file sets vector-valued
grid sizing controls, choose a vertical size at least as large as the level's
vertical cell count. With AMR, SHOC-active refined grids must also cover full
vertical columns.

Keep compile-time availability separate from runtime behavior. Native SHOC may be
built, but it must not change a non-SHOC simulation unless the runtime PBL type
selects native SHOC.

Guard native-SHOC behavior with runtime checks such as:

```cpp
uses_native_shoc()
uses_shoc_family()
```

Do not use compile-time SHOC macros as a substitute for runtime selection.

Keep AMReX portability in mind. Avoid host-only objects in device paths. Avoid
capturing `this` in GPU lambdas. Use AMReX math utilities and `amrex::Real`
literals in device-callable code.

## Testing

Run the native SHOC unit tests after changing this directory:

```bash
ctest --test-dir build -L shoc --output-on-failure
```

Also run representative non-SHOC regression tests after changing shared coupling,
microphysics, plotfile, or time-integration paths. Native SHOC must not alter
non-SHOC results when it is not selected at runtime.

## See also

User documentation:

```text
Docs/sphinx_doc/theory/PBLschemes.rst
```

Build documentation:

```text
Docs/sphinx_doc/buildingLibraryConfig.rst
Docs/sphinx_doc/buildingSystems.rst
```
