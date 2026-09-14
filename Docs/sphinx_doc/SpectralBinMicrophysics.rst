.. _spectral-bin-microphysics:

Spectral-bin microphysics (experimental)
=========================================

ERF's spectral-bin microphysics (SBM) work currently provides an **experimental
P0/P1 infrastructure prototype**, not a production warm-cloud microphysics
scheme.  P1 transports an inert, runtime-sized liquid spectral mass state and
does not implement condensation, evaporation, activation, regeneration,
collision/coalescence, sedimentation, aerosol lifecycle physics, or ice.

Authoritative state and projections
-----------------------------------

The authoritative state is a separate ERF-managed auxiliary ``MultiFab`` with
one component per liquid mass bin.  The generic auxiliary manager owns old,
evaluation, output, and scratch states.  Its face-transfer ledger retains one
reusable stage scratch and one accepted face-transfer accumulator; it does not
retain a full three-stage history.  The runtime component count is
resolved from the input grid; spectral components are not appended to ERF's
compact ``cons`` state.

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
   ``1`` selects the qualified one-moment P1 liquid-mass layout.  ``2`` is a
   distinct P0 two-moment schema contract and is rejected by P1 transport.
``erf.sbm_manufactured_initialization``
   Installs a deterministic positive nonuniform liquid spectrum for the
   qualification case.  Without it, nonzero compact condensate is rejected;
   P1 never guesses bins from bulk ``qc``/``qr``.
``erf.sbm_manufactured_velocity``
   Optional deterministic carrier velocity used by the manufactured ERF
   regression.  It is zero by default and has no effect on ordinary runs.
``erf.sbm_diagnostic_file``
   Optional path for the numerical P1 qualification diagnostic.  It records
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

Attached properties are carrier-bin metadata only in P1.  ``ExtensiveMass``
means additional non-water mass (for example future dry solute), constrained
to be nonnegative but not bounded above by carrier water mass.  ``NumberCarried``
is a nonnegative count-like quantity, while ``MassBoundedSubset`` reserves the
future subset-of-carrier-mass constraint.  All qualified P1 properties use
explicit carrier-bin conservative remapping semantics and are excluded from
the ``qc``/``qr`` water projection.

Supported P1 boundary
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

Explicitly unsupported
----------------------

P1 does not qualify AMR, moving or non-Cartesian terrain, embedded boundaries,
diffusion or implicit moisture diffusion, SHOC/macrophysics, high-order/FCT
transport, sedimentation, condensation/evaporation, activation/regeneration,
collision/coalescence, ice, aerosol lifecycle physics, dynamic grids,
two-moment transport, independent moisture forcing, large-scale or sounding
forcing, sponge/wall modifications, restart/schema conversion, or production
microphysics.  GPU performance, AMR reflux execution, and physical
warm-cloud validation remain future work.  A future non-liquid population can
be represented by the generic population/grid schema without acquiring a
meaningless cloud/rain split, but no non-liquid physics is implemented here.

This mode is infrastructure validation only and must not be interpreted as
production-ready microphysics.
