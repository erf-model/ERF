# ERF SBM P0/P1 source trace

This trace records the pinned ERF integration points used by the P0/P1
implementation.  It is a source map, not a replacement for the v1.0 design.

## Provenance

* ERF base: `33ce039e87f309592609098792dff06d0762438c`
* AMReX submodule: `e60cdc18711ccf7fc0616d7a2fdf062021376976`
* Design Markdown SHA-256: `a0fc579a65c66b18acbf98761eda0796f52e181a2ca57eb52fe2804f48fa6272`
* Implementation branch: `sbm-p0-p1-aux-state`

## Host stage and state path

`Source/TimeIntegration/ERF_MRI.H`, `MRISplitIntegrator::advance`, is the
authoritative stage driver.  The compressible path calls the slow callbacks at
stage output times `t+h/3`, `t+h/2`, and `t+h`, with stage intervals `h/3`,
`h/2`, and `h`.  `S_old` remains the full-step baseline while `S_new` is the
stage state.  The anelastic path calls two callbacks, both over `h`; stage 1
is the predictor and stage 2 applies the Heun correction against the old
state.

`Source/TimeIntegration/ERF_TI_slow_rhs_pre.H` and
`Source/TimeIntegration/ERF_TI_slow_rhs_post.H` bind those callbacks to
`erf_slow_rhs_pre` and `erf_slow_rhs_post`.  The post callback receives the
stage index, old state, current stage state, fast-integrator state, and the
time-averaged carrier moment fields.

`Source/TimeIntegration/ERF_SlowRhsPost.cpp`, `erf_slow_rhs_post`, copies the
stage state into `cur_cons`, constructs scalar advection fluxes with
`AdvectionSrcForScalars`, adds source/diffusion tendencies, and then updates
the stage state.  Its existing anelastic update is
`old + 0.5*((predictor-old) + h*R1)`; its existing compressible update is
`old + h*Rstage` for each callback interval.  Native moisture state values
are clipped nonnegative in this routine.  The SBM provider bypasses those
native `qc`/`qr` operations and supplies projected values from the auxiliary
state instead.

## Carrier flux and geometry

`Source/Advection/ERF_AdvectionSrcForState.cpp`,
`AdvectionSrcForScalars`, and `Source/Advection/ERF_AdvectionSrcForScalars.H`
show that scalar face fluxes are the host time-averaged momentum/carrier
field multiplied by a reconstructed primitive mixing ratio.  For the
second-order donor path the face value is the arithmetic face value selected
by the host advection implementation; the conservative cell tendency is the
negative divergence of those face fluxes.

The carrier fields are built in
`Source/Advection/ERF_AdvectionSrcForRho.cpp` through
`AdvectionSrcForRho` (called from `ERF_SlowRhsPre.cpp`).  `avg_xmom`,
`avg_ymom`, and `avg_zmom` are the time-averaged dry carrier momentum flux
fields produced from the density/continuity update, not velocities reconstructed
from `u`, `v`, and `w`.  Scalar divergence uses `detJ^{-1}`, horizontal map
factor product `mf_mx*mf_my`, and the inverse cell sizes.  P1 is deliberately
restricted to static Cartesian geometry, where the metric factors are unity,
but the auxiliary API retains explicit geometry/metric slots so a later
extension cannot silently reconstruct a different carrier.

## Flux registers and accepted ledgers

`ERF_SlowRhsPre.cpp` and `ERF_SlowRhsPost.cpp` add native fluxes to
`YAFluxRegister` only on the accepted final host stage (`nrk==2` for
compressible and `nrk==1` for anelastic).  This is a final-stage native host
register convention.  P1 does not change it.  The new auxiliary ledger is
separate: compressible acceptance is `h*F2`, while anelastic acceptance is
`h/2*F0 + h/2*F1`, matching the exact recurrence rather than the native
final-stage-only register call.

## Moisture ownership audit

`Source/DataStructs/ERF_DataStruct.H` owns the semantic moisture component
map and input parsing.  `Source/Microphysics/ERF_Microphysics.H` selects the
Eulerian/Lagrangian interface.  `Source/TimeIntegration/ERF_AdvanceMicrophysics.cpp`
is the post-dycore microphysics call site.  `Source/ERF_MakeNewArrays.cpp`
allocates the core state and MRI storage.

For native models, the post slow RHS loops from `RhoQ1_comp` across every
moisture component and can independently advect, diffuse, add sources, and
clip those fields.  Initializers, nudging/large-scale forcing, wall/diffusion
paths, and microphysics can also write the compact moisture state.  The P1
SBM capability validator rejects configurations with such incompatible
features enabled, and the provider-aware path excludes the owned `qc`/`qr`
components from native scalar update/source/diffusion/clip handling.  `qv`
remains the sole ordinary vapor component.  Spectral liquid mass is the
authoritative state; `qc` and `qr` are derived projections.

## P0/P1 capability boundary

Supported by this branch: runtime 1M/2M spectral contracts, fixed cloud/rain
edge projection, typed attached-property metadata, generic ERF-owned
auxiliary state, explicit compressible/anelastic stage context, first-order
donor transfer helpers, and static single-level Cartesian manufactured
transport.

Rejected/fail-closed: AMR, moving terrain, terrain/EB geometry, diffusion or
implicit moisture diffusion, SHOC/macrophysics, FCT/high-order transport,
sedimentation, condensation/evaporation, activation/regeneration,
collision/coalescence, ice, dynamic grids, two-moment transport, independent
custom moisture sources, large-scale/sounding forcing, and restart schema
conversion.  The transport kernel checks every auxiliary component after each
stage and aborts on a material negative or non-finite value; it never clips a
bad state.
