.. _spectral-bin-microphysics:

Spectral-bin microphysics (experimental)
=========================================

ERF's spectral-bin microphysics (SBM) work will represent liquid water and
future aerosol populations with runtime-sized spectral bins.  This branch is
an **experimental P0/P1 infrastructure prototype**, not a production
warm-cloud microphysics scheme: it transports an inert liquid spectrum and
does not implement condensation, evaporation, activation, regeneration,
collision/coalescence, sedimentation, ice, or aerosol lifecycle physics.

State ownership and coupling
-----------------------------

The spectral state is stored in a separate ERF-managed auxiliary
``MultiFab``.  Its component count is resolved at startup from the configured
grid; spectral components are never appended to the compact ``cons`` state.
The generic auxiliary manager keeps old, evaluation, output, and scratch
storage and reports its approximate resident memory at initialization.

The compact ERF moisture state remains:

* ``qv``: authoritative water vapor in ``cons`` and ordinary ERF transport;
* ``qc``: a redundant fixed projection of the cloud-side liquid bins;
* ``qr``: a redundant fixed projection of the rain-side liquid bins.

When SBM is active, the spectral liquid mass is authoritative.  ``qc`` and
``qr`` are projected from the accepted spectral state and are not independently
advected, diffused, sourced, clipped, or advanced by the native bulk moisture
path.

Runtime inputs
--------------

The following ``erf`` inputs are parsed when
``erf.moisture_model = SBM``:

``erf.sbm_nbins``
   Number of liquid-mass bins; the default is 4.  The same executable accepts
   other runtime counts such as 16 and 64.
``erf.sbm_edges``
   ``Nb+1`` nonnegative, strictly increasing bin edges.  If omitted, the
   default grid is ``0, 1, ..., Nb``.
``erf.sbm_pivots``
   ``Nb`` positive pivots, one inside each bin.  If omitted, bin midpoints are
   used.
``erf.sbm_cloud_rain_split``
   Interior edge index separating cloud and rain liquid bins; the default is
   ``Nb/2``.
``erf.sbm_moment_mode``
   ``1`` selects the qualified P1 liquid-mass transport layout.  ``2`` is
   available as a distinct P0 two-moment layout contract, but P1 number
   transport is not yet qualified and therefore fails closed at startup.
``erf.sbm_manufactured_initialization``
   When true, installs a deterministic positive manufactured liquid spectrum
   for transport tests.  When false, a run may start only with an empty
   spectrum (zero ``qc`` and ``qr``); ERF never guesses bins from arbitrary
   nonzero bulk condensate.

Supported P1 configuration
--------------------------

The qualified target is a single AMR level, double-precision, static
Cartesian, triply periodic manufactured case using first-order donor
advection.  It uses the actual ERF dry-air carrier mass-flux fields and the
host compressible RK3 or anelastic Heun stage recurrence.  A material negative
or non-finite auxiliary value fails the run; P1 does not silently clip it.

P1 currently rejects AMR, diffusion or implicit moisture diffusion, FCT or
high-order auxiliary transport, moving terrain, embedded boundaries, SHOC and
other macrophysics, independent moisture forcing, large-scale/sounding
forcing, two-moment transport, unsupported boundary/restart ownership, and
all physical warm-cloud, aerosol, collision, sedimentation, and ice processes.
These are later implementation stages.  GPU qualification, distributed
face-ownership/FCT qualification, restart qualification, and production
microphysics validation are also outside this prototype.

The resolved layout and capability contracts are available through the host
inspection methods in ``Source/Microphysics/SBM``.  They are intended for
future agent/conformance tooling and do not add parsing or string work to the
transport kernels.

This mode is for infrastructure validation only and must not be interpreted as
production-ready microphysics.
