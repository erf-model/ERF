.. _spectral-bin-microphysics:

Spectral-bin microphysics (experimental)
=========================================

ERF's spectral-bin microphysics (SBM) work currently provides an **experimental
P0/P1/P2 infrastructure prototype**, not a production warm-cloud microphysics
scheme.  P2 adds grouped invariant-domain transport, explicit density-weighted
diffusion, conservative auxiliary AMR lifecycle hooks, and strict restart
schema validation.  It still does not implement condensation, evaporation,
activation, regeneration, collision/coalescence, sedimentation, aerosol
lifecycle physics, or ice.

Authoritative state and projections
-----------------------------------

The authoritative state is a separate ERF-managed auxiliary ``MultiFab`` with
runtime-sized components.  A one-moment population stores ``M`` per bin; a
two-moment population stores physical ``(M,C)`` blocks.  The generic auxiliary
manager owns old, evaluation, output, and bounded scratch states.  Its
face-transfer ledger retains one reusable stage scratch and one accepted
face-transfer accumulator; it does not retain a full three-stage history.  The
runtime component count is resolved from the input grid; spectral components
are not appended to ERF's compact ``cons`` state.

P2 transport and constraints
----------------------------

``GroupedFCT_WENOZ3`` reuses ERF's ``WENO_Z3`` reconstruction helper on the
intensive ratio ``X/rho``.  The donor candidate is the low-order state and the
WENO candidate is limited with one face coefficient for every complete
population/bin group.  Groups include one-moment nonnegative mass, two-moment
nonnegative endpoint variables

.. math::

   L = (u C - M)/(u-l), \qquad H = (M-l C)/(u-l),

and attached-property nonnegativity, support, and mass-bounded-subset
constraints.  Only cancellation-sized negative endpoints are normalized;
materially inadmissible states fail collectively with a diagnostic.

The production grouped limiter first accumulates the adverse demand of every
incident face for every cell/group/constraint, communicates those cell-wide
budgets across periodic FAB boundaries, and then applies one common coefficient
to the complete group.  A multi-face counterexample therefore consumes the
combined margin rather than giving every face an independent copy of it.  The
same resolved constraint descriptors are used by the low-order diagnostic,
post-stage validation, post-reflux validation, regrid validation, and restart
validation.

Two-moment checkpoint/storage remains physical ``(M,C)``.  ``(L,H)`` uses the
same per-level scratch allocation only while constructing transport fluxes;
there is no compile-time ``MAX_BINS`` constant.  The current production FCT
remainder still allocates full-layout temporary face candidates and split
low-order fluxes, however.  ``erf.sbm_chunk_size`` is therefore recorded and
validated but is **not yet a P2-qualified peak-memory bound**; see the
qualification report.

When ``erf.sbm_diffusion_coeff`` is positive, the production face flux adds
orthogonal density-weighted diffusion,

.. math::

   F_d = -\rho_f K \nabla(X/\rho),

to the same face ledger as advection.  It does not use ``-K nabla X`` and does
not permit ERF's native moisture diffusion path to write provider-owned SBM
components.

The P2 auxiliary lifecycle is attached to ERF's level creation, coarse-to-fine
creation, stage-time coarse/fine FillPatch, remake, clear, average-down, and
checkpoint read/write callbacks.  Restriction is volume weighted; new fine
levels and newly covered fine cells use the authoritative coarse spectrum at
the requested time.  Checkpoint restart requires an exact SBM schema match,
compares checkpointed ``qc``/``qr`` projections before restoring the
authoritative compact fields, and is exercised by a two-rank continuous versus
restart AMR equivalence test.

The compact ERF moisture fields have the following ownership when SBM is
active:

* ``qv`` remains authoritative and follows ERF's native vapor path;
* ``qc`` and ``qr`` are redundant host-coupling projections of the liquid
  spectrum, with the configured interior bin split;
* ``qc`` and ``qr`` are not independently advected, diffused, sourced,
  clipped, or advanced by bulk microphysics.

For every stage, the production donor kernel constructs the spectral
face-centered fluxes once.  Those same fluxes form the spectral divergence and
are recorded in ``AuxiliaryFaceTransferLedger``.  Bulk face transfers are
obtained by summing the accepted spectral face components, so the flux-level
projection is preserved rather than reconstructed from cell-centered output.
Compressible acceptance is ``dt * F(stage 2)``; anelastic Heun acceptance is
``0.5 * dt * (F(stage 0) + F(stage 1))``.

For a completed supported step, the production diagnostic evaluates the local
accepted-transfer identity for every spectral bin,

.. math::

   X^{n+1} - X^n + \nabla_h\!\cdot I = 0,

and the corresponding identities for the accepted ``qc`` and ``qr`` face
transfers.  The residual tolerances are scaled by the state/transfer magnitude
and machine epsilon, with no order-one floor (an exactly zero scale has
exactly zero tolerance).  The compact face transfers used in this check are
the exact projections of the accepted spectral transfer.  The compact
old-state snapshot is copied from ERF's explicit ``state_old`` argument after
ERF's old/new swap, and the real qualification diagnostic records
``step_count``; the six qualification cases require at least two complete
steps.

The authoritative-state finite/nonnegative check uses one GPU-capable fused
``ParReduce`` traversal across all requested components, followed by a fixed
number of scalar MPI reductions.  Its negative tolerance is
``128 * epsilon * max_abs(state)``; an all-zero state therefore has zero
tolerance.  This avoids one global ``min`` collective per spectral bin.

The memory fields in the qualification diagnostic are logical allocated data
payloads, including ghost cells in every grown FAB, for the four auxiliary
cell states, reusable/accepted face transfers, and compact baseline snapshot.
Allocator metadata is not included.

Runtime inputs
--------------

These ``erf`` inputs are read only when ``erf.moisture_model = SBM``:

``erf.sbm_nbins``
   Runtime number of liquid-mass bins.  It must be in ``[2, 1000000]`` and is
   validated before any vector resize.  The qualification tests use 4, 16,
   and 64 with one executable.
``erf.sbm_edges``
   Optional ``Nb+1`` finite, nonnegative, strictly increasing coordinate edges.
   The default is ``0, 1, ..., Nb``.
``erf.sbm_pivots``
   Optional ``Nb`` finite positive pivots, each inside its bin.  Midpoints are
   used when omitted.
``erf.sbm_cloud_rain_split``
   Interior liquid-population bin index.  It belongs to projection metadata,
   not to the generic spectral coordinate grid; the default is ``Nb/2``.
``erf.sbm_moment_mode``
   ``1`` selects the qualified one-moment liquid-mass layout.  ``2`` selects
   physical ``(M,C)`` checkpoint/storage with nonnegative ``(L,H)`` transport
   coordinates and P2 realizability constraints.
``erf.sbm_transport_method``
   ``DonorCell`` selects the qualified P1 reference path.  ``GroupedFCT_WENOZ3``
   uses ERF's WENO-Z3 helper and one invariant-domain limiter coefficient for
   each complete population/bin group.
``erf.sbm_diffusion_coeff``
   Nonnegative explicit orthogonal two-point coefficient in ERF's Cartesian
   length/time units.  Positive values diffuse ``X/rho`` with the
   density-weighted face conductance; implicit/native moisture diffusion is
   not part of P2.
``erf.sbm_chunk_size``
   Positive chunk-policy value recorded in capability and restart identities.
   The host/reference limiter honors it.  The current production WENO/FCT
   path still retains full-layout temporary face buffers, so this input does
   not yet establish a qualified production peak-memory bound.
``erf.sbm_manufactured_initialization``
   Installs a deterministic positive nonuniform liquid spectrum for the
   qualification case.  Without it, nonzero compact condensate is rejected;
   SBM never guesses bins from bulk ``qc``/``qr``.
``erf.sbm_manufactured_velocity``
   Optional deterministic carrier velocity used by the manufactured ERF
   regression.  It is zero by default and has no effect on ordinary runs.
``erf.sbm_diagnostic_file``
   Optional path for the numerical SBM qualification diagnostic.  It records
   mass conservation, nonuniform transport, state projection, and accepted
   face-flux projection values for an independent CTest checker.

Units and schema semantics
--------------------------

Spectral coordinates are phase-neutral: the grid describes ``Mass`` or
``Radius`` mathematics, while population metadata supplies the phase and
semantic ID.  Thus a future ice or aerosol mass coordinate does not require a
new phase-specific grid type.  The liquid-mass spectral coordinate and the
transported state also have different dimensions.  The coordinate metadata is
``coordinate_units = kg``.  A transported mass-bin component is
density-weighted and therefore has
``mass_state_units = kg m^-3``; a two-moment number block, when supported by a
future implementation, has ``number_state_units = m^-3``.  These fields are
part of the stable layout/inspection schema and are not cosmetic labels.

Attached properties are carrier-bin metadata in P1/P2.  ``ExtensiveMass``
means additional non-water mass (for example future dry solute), constrained
to be nonnegative but not bounded above by carrier water mass.  ``NumberCarried``
is a nonnegative count-like quantity, while ``MassBoundedSubset`` reserves the
future subset-of-carrier-mass constraint.  All qualified P1/P2 properties use
explicit carrier-bin conservative remapping semantics and are excluded from
the ``qc``/``qr`` water projection.

Supported P1 baseline
---------------------

The qualified P1 configuration is:

* one AMR level;
* double precision;
* static Cartesian, fully periodic geometry;
* exactly one liquid population with one-moment spectral mass;
* first-order donor spectral transport using ERF's actual dry-air carrier
  mass-flux fields;
* compressible RK3 or anelastic Heun stage integration;
* finite, nonnegative authoritative spectral state with no silent clipping.

The real manufactured regressions cover both time integrators and 4, 16, and
64 bins.  They verify spatial change, per-bin nonnegativity, periodic mass
conservation, compact projection, accepted face-transfer projection, and local
accepted-transfer closure after two complete steps.  The
variable-density free-stream unit test remains a separate preservation test;
it does not replace the nonuniform transport regression.

Supported P2 matrix
-------------------

The P2 implementation extends the baseline to runtime-sized one- and
two-moment auxiliary populations, ``DonorCell`` or ``GroupedFCT_WENOZ3``
transport, complete grouped constraints, explicit orthogonal two-point
diffusion, and the auxiliary AMR/restart lifecycle.  Direct production
qualification covers static fully periodic Cartesian hierarchy cases, including
the two-level 2M fixture and continuous/restart comparison.  Restriction is
physical-volume weighted and prolongation is conservative piecewise constant.
Accepted spectral transfers—not independently reconstructed compact fields—drive
``qc``/``qr`` projections, boundary budgets, and AMR registers.  Full P2
completion remains pending production chunk-bounded work and broader numerical
AMR/MPI oracles.

WENO qualification evidence
----------------------------

The focused smooth periodic operator test measures the implemented WENO-Z3
face reconstruction against the analytic periodic face value and compares it
with the donor reconstruction.  The machine-readable output is written to
``/private/tmp/erf_sbm_p2_weno_convergence.csv`` by
``SBMP2.WENOZ3ConvergenceBeatsDonorOnPeriodicSmoothOperator``.  The current
double-precision measurements are:

.. list-table:: Smooth periodic reconstruction errors
   :header-rows: 1

   * - N
     - WENO error
     - Donor error
     - WENO order
     - Donor order
   * - 8
     - 2.1677e-1
     - 3.8268e-1
     - --
     - --
   * - 16
     - 5.6906e-2
     - 1.9509e-1
     - 1.9295
     - 0.9720
   * - 32
     - 1.4400e-2
     - 9.8018e-2
     - 1.9826
     - 0.9930
   * - 64
     - 3.6107e-3
     - 4.9068e-2
     - 1.9957
     - 0.9983

This qualifies the measured smooth operator and donor comparison only; it is
not a universal third-order claim for the nonlinear ERF time integrator.

The first physical boundary descriptor layer provides periodic,
prescribed-spectral-inflow, advective-outflow, and impermeable-wall semantics
at the generic transfer-service level.  These nonperiodic descriptors are
reference/service infrastructure, not production-qualified ERF SBM boundary
handling.  The production gate is periodic Cartesian only.  A prescribed
inflow must contain a complete realizable group state; an impermeable wall has
zero resolved spectral transfer and does not imply deposition.

Post-reflux admissibility is checked explicitly.  If a conservative hierarchy
correction leaves the invariant domain, P2 fails collectively with diagnostic
state rather than clipping, renormalizing, or silently dropping the transfer.
Restart requires an exact P2 schema match and compares checkpointed compact
projections before any reconstructed projection overwrites them.  A synthetic
mass-subset population/property is an architecture test only; it is not ice
physics.

Explicitly unsupported
----------------------

P2 does not qualify implicit moisture diffusion, SHOC/macrophysics condensate
coupling, moving terrain, embedded boundaries, non-orthogonal diffusion,
dynamic spectral grids, schema conversion, P3 thermodynamics,
condensation/evaporation, activation/regeneration, collision/coalescence,
sedimentation, physical ice, aerosol lifecycle physics, independent moisture
forcing, large-scale or sounding forcing, sponge/wall modifications, or
production warm-cloud microphysics.  GPU performance and physical warm-cloud
validation remain future work.  The synthetic non-liquid population/property
case is numerical extensibility coverage only; no non-liquid process is
implemented here.

This mode is infrastructure validation only and must not be interpreted as
production-ready microphysics.
