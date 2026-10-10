.. _sec:SpectralBinMicrophysics:

Spectral-Bin Microphysics (SBM)
===============================

Spectral-bin microphysics represents condensed water explicitly over a set of
particle-size or particle-mass bins. This differs from a bulk microphysics
scheme, which predicts a small number of bulk quantities such as cloud-water
and rain-water mixing ratio and generally represents the underlying particle
size distribution through an assumed functional form.

ERF's spectral-bin state is Eulerian: each atmospheric grid cell carries a
bin-resolved distribution. This is also different from the
:ref:`sec:SuperDroplets` method, where the condensed phase is represented by
Lagrangian computational particles that move through the Eulerian model grid.

A complete spectral-bin microphysics model can explicitly represent the
particle-size distributions of aerosols and hydrometeors and their evolution
through physical processes. ERF's currently exposed ``SBM`` option is
narrower: it predicts a bin-resolved liquid-water distribution and advances
it by resolved atmospheric advection on a single grid level. The code also
contains independently tested state-time management and conservative
coarse/fine spectral-transfer operators. These operators are foundations
for future mesh-refined simulations, not an enabled multilevel capability.
Aerosol evolution and cloud microphysical processes are not yet included.
See :ref:`AuxiliaryState` for the reusable mapped-transport interface.

.. warning::

   The current ``SBM`` option provides resolved advection of a liquid-water
   size distribution on one grid level. It is not yet a complete warm-cloud
   or warm-rain microphysics scheme.

   The liquid spectrum is the authoritative condensed-water state. It is
   advected using ERF's dry-air mass carrier, high-order reconstruction,
   coupled flux-corrected transport, and fixed projections to conventional
   cloud- and rain-water fields.

   The runtime does not yet provide particle diffusion or turbulent mixing,
   condensation/evaporation, aerosol activation or evolution,
   collision-coalescence, sedimentation, precipitation, multilevel spectral
   transport, or the full nonperiodic spectral-boundary lifecycle. Moving
   terrain and embedded boundaries remain unsupported.

   Unsupported configurations are rejected rather than silently reverting
   to bulk-water transport.

Current spectral state
----------------------

The current user-facing SBM configuration contains one liquid-water
population. Aerosol and ice populations are not currently exposed through the
runtime ``SBM`` option.

Water vapor remains part of ERF's normal prognostic moisture state. The
condensed liquid-water distribution is stored separately as the authoritative
spectral state.

For liquid bin :math:`i`, let :math:`q_i` denote the liquid-water mass in that
bin per unit mass of dry air, and let :math:`\rho_d` be the dry-air density.
The quantity stored by the current spectral state is the density-weighted bin
mass

.. math::

   M_i = \rho_d q_i,

with units of :math:`\mathrm{kg\,m^{-3}}`.

In the two-moment representation, ERF additionally stores the droplet-number
density in each bin,

.. math::

   C_i = \rho_d N_i,

where :math:`N_i` is the bin number per unit mass of dry air. The stored
quantity :math:`C_i` therefore has units of :math:`\mathrm{m^{-3}}`.

The spectral coordinate used by the current ERF SBM runtime configuration is
the liquid-water mass of an individual particle, in kg. Bin edges and pivots
therefore describe particle water mass, not atmospheric liquid-water content.
The reference fixed-pivot and interval remapping policies require
individual-particle liquid-water mass as their spectral coordinate. A spectrum
expressed in another coordinate, such as particle radius, is a different
representation because its number measure and moments transform with the
coordinate change. The current reference policies therefore reject non-mass
coordinates rather than interpreting radius numerically as though it were
particle mass.

Bulk cloud and rain fields
--------------------------

ERF's normal atmospheric state still contains the familiar cloud-water and
rain-water fields so that the rest of the model can interact with a conventional
moisture-state layout.

For SBM these fields are not independent prognostic reservoirs. They are fixed
projections of the authoritative spectral liquid-water distribution.

Let ``erf.sbm_cloud_rain_split = s``. Bins with indices smaller than ``s`` are
classified as cloud water, and bins with indices greater than or equal to
``s`` are classified as rain water. ERF therefore computes

.. math::

   \rho_d q_c = \sum_{i < s} M_i

and

.. math::

   \rho_d q_r = \sum_{i \ge s} M_i.

The cloud/rain boundary is consequently the spectral bin edge at index
``s``.

In ERF's compact moisture state, water vapor occupies ``Q1``, projected cloud
water occupies ``Q2``, and projected rain water occupies ``Q3``. The spectral
distribution, rather than ``Q2`` or ``Q3``, is the source of truth for
condensed liquid water.

ERF protects the projected cloud- and rain-water components from the normal
native liquid-water update paths while SBM is active. They are regenerated
from the spectrum at the required model handoffs.

One- and two-moment bins
------------------------

ERF currently supports one- and two-moment spectral bins.

Here ``one moment`` and ``two moments`` refer to the moments stored within each
spectral bin. They should not be confused with conventional one- or two-moment
bulk cloud and precipitation parameterizations.

One-moment representation
~~~~~~~~~~~~~~~~~~~~~~~~~

With

::

   erf.sbm_moment_mode = 1

ERF stores one quantity per bin: the density-weighted liquid-water mass
:math:`M_i`.

For :math:`N` bins, the spectral state therefore contains

.. math::

   M_0,\ M_1,\ \ldots,\ M_{N-1}.

Each spectral bin also has a representative pivot particle mass. The pivot is
part of the spectral-grid definition and restart identity, but droplet number
is not independently stored in the one-moment representation.

Two-moment representation
~~~~~~~~~~~~~~~~~~~~~~~~~

With

::

   erf.sbm_moment_mode = 2

ERF stores both liquid-water mass and droplet number in every bin.

For :math:`N` bins the runtime component order is

.. math::

   M_0,\ M_1,\ldots,M_{N-1},
   C_0,\ C_1,\ldots,C_{N-1}.

Consider a bin whose lower and upper particle-mass edges are
:math:`a_i` and :math:`b_i`, with :math:`a_i < b_i`. A physically admissible
two-moment state must satisfy

.. math::

   C_i \ge 0,

and

.. math::

   a_i C_i \le M_i \le b_i C_i.

When :math:`C_i > 0`, this is equivalent to requiring the mean particle mass

.. math::

   \frac{M_i}{C_i}

to lie inside the bin. If :math:`C_i = 0`, the same constraint requires
:math:`M_i = 0`.

ERF checks these relationships during fixture initialization and again when
reading the authoritative spectral state from a checkpoint. A materially
unrealizable two-moment state is rejected rather than clipped or silently
repaired.

Reference spectral reconstruction and remapping
------------------------------------------------

This section concerns redistribution in the spectral coordinate---that is,
how a particle population is represented after a local process changes
particle liquid-water mass. It is distinct from advection between atmospheric
grid cells, AMR transfer between spatial resolutions, and conversion of a
checkpoint from one spectral schema to another.

The current runtime does not yet execute cloud microphysical processes that
use these remappers. The present implementation establishes the reference
representation and remapping contracts that future process implementations
must obey.

Moment normalization
~~~~~~~~~~~~~~~~~~~~

The prognostic spectral state stored by ERF is density weighted. Within one
atmospheric grid cell it is convenient to state the remapping algebra using
the corresponding quantities per unit mass of dry air,

.. math::

   N_i = \frac{C_i}{\rho_d},
   \qquad
   q_i = \frac{M_i}{\rho_d},

where :math:`N_i` is particle number per unit mass of dry air and :math:`q_i`
is liquid-water mass per unit mass of dry air. The stored quantities are
:math:`C_i=\rho_d N_i` and :math:`M_i=\rho_d q_i`.

Because :math:`\rho_d` is a common factor within the cell, the same linear
remapping formulas apply to the density-weighted state. For example, if a
packet number :math:`N_p` is expressed per unit mass of dry air, the
corresponding number-density packet is :math:`C_p=\rho_d N_p`.

One-moment fixed-pivot remapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the one-moment representation, bin :math:`i` stores its water moment
:math:`q_i` but not an independent number moment. Particle number is diagnosed
at the fixed positive pivot mass :math:`p_i`,

.. math::

   N_i = \frac{q_i}{p_i}.

Consider a packet with number :math:`N_p` and actual liquid-water mass
:math:`m` per particle. If :math:`m` lies between adjacent pivots
:math:`p_-` and :math:`p_+`, the packet is divided between them using the
number weights

.. math::

   w_- = \frac{p_+-m}{p_+-p_-},
   \qquad
   w_+ = \frac{m-p_-}{p_+-p_-}.

The corresponding increments are

.. math::

   \Delta N_\pm = N_p w_\pm,
   \qquad
   \Delta q_\pm = p_\pm N_p w_\pm.

Any attached extensive inventory follows the same number weights. The mapping
therefore conserves the implied packet number, liquid water, and each attached
inventory to floating-point roundoff, even though only the water moment and
attached inventories are persisted in one-moment mode.

Fixed-pivot remapping necessarily introduces spectral spreading for an
off-pivot packet. If

.. math::

   Q_2 = \int m^2\,d\mu

denotes the second particle-mass moment of the number measure, the increment
relative to retaining the packet at its actual mass is

.. math::

   \Delta Q_2
   =
   N_p(m-p_-)(p_+-m)
   \ge 0.

This is numerical broadening introduced by the representation; it is not
physical cloud-spectrum broadening produced by a microphysical process.

For a packet below the first pivot, :math:`0<m<p_0`, define

.. math::

   f_{\mathrm{wet}} = \frac{m}{p_0}.

The reference closure places the fraction :math:`f_{\mathrm{wet}}` of the
packet number at the first pivot and returns the remaining number through the
zero-water residual path,

.. math::

   \Delta N_0 = N_p f_{\mathrm{wet}},
   \qquad
   \Delta q_0 = N_p m,
   \qquad
   N_{\mathrm{res}} = N_p(1-f_{\mathrm{wet}}).

Attached inventories are partitioned by the same
:math:`f_{\mathrm{wet}}` and :math:`1-f_{\mathrm{wet}}` fractions. This is a
numerical sub-pivot evaporation closure. It should not be interpreted as a
statement that the unresolved physical particles represented by the residual
have necessarily undergone complete physical evaporation.

A zero-water packet is returned entirely through the residual path. A packet
above the largest fixed pivot reports overflow and is not clipped into the
largest bin.

Two-moment interval representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the two-moment representation, number and liquid-water mass are stored
independently. For a populated bin with edges :math:`a_i` and :math:`b_i`,

.. math::

   \bar m_i = \frac{q_i}{N_i}
            = \frac{M_i}{C_i}

may occupy the full admissible interval.

The reference ``mean-delta`` reconstruction is a temporary monodisperse
representation that places all of the bin number at :math:`\bar m_i`. It
reproduces the persisted number and water moments to floating-point roundoff;
it is not an additional prognostic description of within-bin shape.

Canonical persisted two-moment state requires a populated bin to have positive
number and positive water mass, with its mean inside that bin's routing-owned
interval. Interior bins own :math:`[a_i,b_i)`; the final bin owns
:math:`[a_i,b_i]`. A bin with zero number is empty only when its water mass and
every attached extensive inventory are also exactly zero. A positive-number,
zero-water liquid-bin state is invalid; a zero-water packet is handled by the
explicit residual path.

The endpoint transform uses a small floating-point tolerance when testing
moment realizability. That tolerance is a numerical aid for the transform; it
does not authorize silent normalization of authoritative persisted state. A
persisted bin must already satisfy its canonical routing-owned interval and
must be numerically reconstructable in the active ``amrex::Real`` precision.
In particular, a mathematically positive stored quantity whose required
per-particle quotient underflows to zero is rejected rather than interpreted
as an empty particle population.

The first two moments do not uniquely determine that shape. For any
nonnegative number distribution supported on :math:`[a_i,b_i]` with
:math:`N_i>0`, the second mass moment satisfies

.. math::

   \frac{q_i^2}{N_i}
   \le
   Q_{2,i}
   \le
   (a_i+b_i)q_i-a_i b_i N_i.

A delta distribution at :math:`\bar m_i` attains the lower bound. On the
closed mathematical support :math:`[a_i,b_i]`, an appropriate mixture at the
two endpoints attains the upper bound. For the canonical routing ownership
used here, an interior bin owns :math:`[a_i,b_i)`, so the same expression is
the upper bound (and limiting supremum) for that bin while a particle exactly
at :math:`b_i` belongs to the next bin. The final bin includes its global upper
endpoint. The
mean-delta reference reconstruction therefore selects the minimum-variance
distribution consistent with the two stored moments. Quantities that depend
on unresolved within-bin structure---for example nonlinear collision rates,
size-dependent sedimentation, or optical properties---are not uniquely
determined by :math:`N_i` and :math:`q_i` alone.

For an attached extensive property, the mean-delta reference closure assigns
the bin-mean amount per particle to the reconstructed node. This reproduces
the stored attached-property inventory to floating-point roundoff, but it
likewise does not resolve covariance between particle composition and
liquid-water mass.

If several physically distinct packets are accumulated in the same interval,
their total number and first water-mass moment are retained, but their
within-bin spread is not an additional persisted degree of freedom.
Subsequent mean-delta reconstruction therefore cannot recover that spread.
This is representation loss, not physical coalescence: particle number has
not been reduced by the representation itself.

Two-moment packet deposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-moment process packet is deposited using its actual particle mass,
rather than being replaced by neighboring pivot masses. If a packet with
number :math:`N_p` and particle mass :math:`m` lies in interval :math:`i`,

.. math::

   \Delta N_i = N_p,
   \qquad
   \Delta q_i = N_p m.

Any attached extensive inventory is deposited into the same interval.

Interior intervals are half-open: a packet exactly on an edge shared by two
bins belongs to the upper bin. The global upper edge is included in the final
bin. Materially out-of-range packets return an underflow or overflow status
rather than being silently clipped.

An unchanged mean-delta reconstruction redeposits into its original interval.
A state concentrated exactly on an interior shared edge belongs to the upper
bin; the same state stored in the lower bin is noncanonical and is rejected
rather than silently migrated by a no-process reconstruction/remapping cycle.

At a positive global lower edge or at the global upper edge, the implementation
permits only a very small, explicitly bounded floating-point exception: a
nonnegative particle mass up to four representable floating-point steps
outside the boundary may be classified as roundoff. Such a packet is deposited
at the exact boundary mass :math:`m_b`, not at its slightly out-of-support
input mass :math:`m_{\mathrm{in}}`. Negative particle masses are invalid; in
particular, a spectrum whose lower edge is exactly zero does not admit
negative roundoff excursions below that edge.

The corresponding signed numerical water correction is

.. math::

   \delta q_{\mathrm{round}}
   =
   N_p\left(m_b-m_{\mathrm{in}}\right),

with the analogous density-weighted correction obtained by multiplication by
:math:`\rho_d`. A positive correction means boundary normalization increased
the persisted water amount relative to the incoming floating-point packet; a
negative correction means it decreased it. This quantity is explicit
numerical roundoff accounting, not a physical condensation, evaporation, or
precipitation source.

Storing the exact boundary mass keeps the accepted state inside its declared
spectral support and lets reconstruction followed immediately by projection
preserve the accepted moments to floating-point roundoff.

A positive-number packet with exactly zero liquid-water mass is returned
through the zero-water residual path rather than retained as a populated
liquid bin. Residual number and attached inventory are returned to the caller
for process-level handling; the reference remapper itself does not create or
modify an aerosol population.

A mathematically positive required quantity that underflows to exactly zero in
the active floating-point precision is rejected fail-closed rather than
reinterpreted as physical evaporation, residual material, or an empty
population. Packet routing and application report this case as
``NumericalUnderflow``; reconstruction and projection reject it through their
existing invalid/failure returns. Packet application remains atomic.

Restart and scientific identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The representation, temporary reconstruction, and packet-remapping policies
are versioned parts of the SBM scientific identity. The current one-moment
identities are ``fixed-pivot-1m-v1``, ``fixed-pivot-delta-v1``, and
``fixed-pivot-two-center-v1``. The current two-moment identities are
``interval-2m-v1``, ``mean-delta-2m-v1``, and
``actual-mass-interval-v1``.

These strings are restart and provenance metadata; they are not user-selectable
microphysics options. Restart comparison is exact, and ERF does not
automatically convert a checkpoint to a different representation or remapping
policy.

The representation, transport, and projection policies are versioned parts of
the exact restart identity. Older zero-transport schemas are rejected rather
than silently interpreted as the currently supported advection scheme.

Auxiliary prognostic-state transport
------------------------------------

The SBM liquid spectrum is the authoritative liquid-water state. ERF advances
this spectrum in physical space using the mapped auxiliary-state transport
infrastructure and the atmospheric dynamical core's dry-air mass flux.

The low-order candidate is conservative donor transport. The high-order
candidate uses ERF's native WENO-Z3 reconstruction. Two-moment liquid bins are
reconstructed in endpoint-number coordinates before their face rates are
transformed back to physical mass and number. A shared atomic-group FCT
limiter constrains mass, number, and attached extensive properties together.

Projected ``qc`` and ``qr`` remain fixed diagnostics of the accepted liquid
spectrum. They are refreshed from the spectrum and do not receive an
independent liquid advection, diffusion, clipping, or microphysical update.

Available end-to-end simulation evidence is deliberately bounded. It covers
single-level, triply periodic advection on constant- and stretched-z grids,
plus a one-step static Cos4Hill terrain-fitted case with terrain flat at the
periodic boundaries. Direct transport tests also exercise nontrivial fluxes
with a synthetic static mapped measure. These tests do not establish a
nonperiodic terrain-boundary capability. Particle diffusion and turbulent
mixing, multilevel spectral synchronization and reflux, nonperiodic spectral
boundaries, condensation/evaporation, activation, collision-coalescence,
sedimentation, and precipitation are not available in the runtime. Moving
terrain and embedded boundaries remain unsupported.

Conservative spectral transfer between grid resolutions
-------------------------------------------------------

Atmospheric simulations may use locally refined meshes around clouds,
convective updrafts, or sharp moisture gradients. For spectral-bin
microphysics, a change in atmospheric grid resolution must transfer the
liquid-water and droplet-number inventories represented by each
**spectral bin**. Transferring only bulk cloud-water and rain-water amounts would not
preserve the underlying particle-size distribution: many different spectra
can have the same bulk mixing ratios.

The following state and transfer operations have been implemented and tested
independently. They are needed to extend spectral-bin transport to refined
atmospheric meshes, but are not yet part of the operational multilevel SBM
runtime:

* **Time-matched spectral states.** The bin-resolved distribution has explicit
  previous and updated time views that follow ERF's atmospheric time advance.
  A quantity expressed per unit dry-air mass is formed using dry-air density
  from the matching time view, not from a different model stage.

* **Conservative coarse-grid averaging.** For spectral component :math:`a`,
  let :math:`U_a` be the stored density-weighted amount. Depending on the
  component, this represents liquid-water mass per unit volume, droplet
  number per unit volume, or an attached extensive inventory. ERF's mapped
  conservative amount is

  .. math::

     H_a = \omega U_a, \qquad
     \omega = \frac{\det J}{m_x m_y},

  where :math:`\det J` is the terrain-coordinate Jacobian and :math:`m_x`
  and :math:`m_y` are the horizontal map factors. For a coarse cell :math:`c`
  completely covered by fine cells :math:`j`, restriction uses

  .. math::

     H_{a,c} = \frac{1}{N_c}\sum_{j\in\mathcal{C}(c)} H_{a,j},
     \qquad
     U_{a,c} = \frac{H_{a,c}}{\omega_c},

  where :math:`\mathcal{C}(c)` contains the :math:`N_c` fine children of
  cell :math:`c`. This preserves the represented mapped inventory to
  floating-point accuracy. Coarse cells outside the fine-covered region
  retain their previous spectral values exactly.

* **Carrier-relative fine-grid transfer.** The operator forms the amount of each
  component per unit mass of dry air,

  .. math::

     z_a = \frac{U_a}{\rho_d},

  where :math:`\rho_d` is dry-air density. The fine-grid value is obtained
  through piecewise-constant interpolation of :math:`z_a` followed by
  reconstruction with the **fine-grid density at the same semantic time**:

  .. math::

     U_{a,\mathrm{fine}}
     = \rho_{d,\mathrm{fine}}\,
       z_{a,\mathrm{parent\ coarse}}.

  Refining the atmospheric mesh thereby transfers a bin's dry-air-relative
  abundance without treating the change in spatial resolution as a physical
  droplet-growth, evaporation, or collision process.

* **Accepted spectral face transfers.** The implementation defines and tests a completed-step
  transport record for every spectral component,

  .. math::

     \mathcal{I}_{a,f}
     = \int_{t^n}^{t^{n+1}}
       \widetilde{F}_{a,f}^{\mathrm{accepted}}(t)\,dt,

  where :math:`\widetilde{F}_{a,f}^{\mathrm{accepted}}` is the accepted
  mapped transport rate through face :math:`f`. The record uses the rates
  accepted after the coupled transport limiter, not an unrestricted
  high-order proposal or the two bulk cloud/rain sums. It provides the
  information needed for conservative coarse/fine flux correction in a
  subsequent AMR implementation.

* **Physical and numerical admissibility.** A transferred candidate spectrum
  must still satisfy the declared mass, number, attached-property, and
  spectral-bin ownership constraints. Invalid or numerically
  unrepresentable candidates are rejected rather than repaired by clipping
  or silent redistribution.

.. important::

   **These are tested grid-transfer foundations, not an enabled multilevel
   SBM capability.** Although ERF supports adaptive mesh refinement for other
   atmospheric quantities, simulations using ``erf.moisture_model = SBM``
   still require ``amr.max_level = 0`` and periodic boundaries in all three
   directions. Coarse/fine temporal interpolation, accepted-flux correction
   (reflux), regridding, and multilevel spectral restart must be integrated
   and qualified before SBM can run on refined meshes. These grid-transfer
   foundations introduce no new user input parameters and no cloud or aerosol
   microphysical source processes. The currently exposed spectrum represents
   liquid water only.

What is still needed for multilevel simulations
-----------------------------------------------

The present SBM runtime advects the liquid-water spectrum on one atmospheric
grid level. The conservative grid-transfer operators described above have
been tested separately, but they are not yet connected to ERF's complete
adaptive-mesh lifecycle.

A future multilevel simulation must interpolate coarse-grid spectral state
and dry-air density to the same fine-grid time; reconcile the accepted
bin-resolved transport across coarse/fine interfaces by conservative flux
correction (reflux); and verify that the synchronized droplet-mass,
number, and attached-property inventories remain physically admissible.
Grid creation, regridding, restart, and physical boundary conditions must
also preserve the full spectral state, rather than reconstructing a size
distribution from bulk cloud-water and rain-water mixing ratios.

A conservation budget for a refined simulation must count each atmospheric
region only once, using fine-grid values where refinement is present and
excluding the underlying coarse-grid cells. Multilevel spectral transport
will remain unavailable until these coupled operations have been tested in
ERF, including their failure behavior.

Numerical reconstruction for resolved spectral advection
--------------------------------------------------------

The current high-order spectral advection method is
``EndpointNumberWENOZ3``. This numerical-method identity is fixed by the
implementation and is not a user-selectable SBM option. See
:ref:`AuxiliaryState` for the mapped conservation and host-stage transport
conventions.

For the current implementation, ``EndpointNumberWENOZ3`` has a specific
numerical meaning:

* every transported component first uses the intensive dry-air coordinate
  :math:`z_a=U_a/\rho_d` formed with the density from the same host stage;
* each two-moment mass/number pair is transformed to endpoint-number
  coordinates before high-order reconstruction and converted back to physical
  mass and number face rates before grouped FCT;
* one-moment components and attached extensive properties remain in their
  physical intensive coordinates;
* no per-component rescaling or normalization is applied before WENO-Z3;
* reconstruction uses ERF's native ``AdvType::Weno_3Z`` / ``WENO_Z3`` stencil.
  At this revision its linear weights are :math:`1/3` and :math:`2/3`, with a
  fixed epsilon of ``1.e-12`` in ``AMREX_USE_FLOAT`` builds and ``1.e-40``
  otherwise; and
* the enabled SBM runtime is triply periodic, so its end-to-end simulations
  do not exercise a physical-boundary high-order fallback. The declared
  spectral boundary policy permits an admissible donor proposal when a
  physical-boundary high-order stencil is not qualified.

Changing the reconstruction basis, component normalization, WENO epsilon or
stencil, or physical-boundary fallback is a numerical-method change requiring
new qualification; it is not a formatting-only cleanup.

SBM runtime inputs
------------------

The following inputs use the ``erf.`` prefix.

``erf.moisture_model``
~~~~~~~~~~~~~~~~~~~~~~

Select the current single-level, advection-only SBM capability with

::

   erf.moisture_model = SBM

New SBM input files should omit the obsolete no-transport switch described below.

``erf.sbm_zero_transport_fixture``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** Boolean

**Default:** ``false``

This obsolete switch belonged to an earlier infrastructure-only test mode.
It is not a supported option for current simulations. Omit it from new SBM
input files. An explicitly true value is rejected rather than silently
disabling spectral transport.

``erf.sbm_nbins``
~~~~~~~~~~~~~~~~~

**Type:** Integer

**Default:** ``4``

**Requirement:** at least 2

Sets the number of liquid-water bins when ``erf.sbm_edges`` is not supplied.

If explicit edges are supplied, their length determines the effective number
of bins:

.. math::

   N_{\mathrm{bins}} = N_{\mathrm{edges}} - 1.

In that case the explicit edge array, rather than the separately specified
``erf.sbm_nbins``, determines the spectral-grid size.

``erf.sbm_moment_mode``
~~~~~~~~~~~~~~~~~~~~~~~

**Type:** Integer

**Default:** ``1``

**Allowed values:** ``1`` or ``2``

``1`` stores one liquid-mass moment per bin.

``2`` stores liquid mass and droplet number per bin, with all mass components
followed by all number components.

``erf.sbm_cloud_rain_split``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** Integer bin index

**Default:** ``sbm_nbins / 2`` using integer division

The value must lie strictly inside the liquid-bin range:

.. math::

   0 < s < N_{\mathrm{bins}}.

Bins ``0`` through ``s-1`` contribute to projected cloud water. Bins ``s``
through ``N-1`` contribute to projected rain water.

The split is an index into the spectral grid. The corresponding physical
cloud/rain boundary is the particle-mass edge ``sbm_edges[s]``.

``erf.sbm_edges``
~~~~~~~~~~~~~~~~~

**Type:** List of real values

**Units:** kg of liquid water per particle

Defines the spectral-bin edges.

For :math:`N` bins, supply :math:`N+1` values. The edges must be finite,
nonnegative, strictly increasing, and sufficiently separated to form a
numerically well-conditioned grid.

If ``erf.sbm_edges`` is omitted, ERF currently generates a logarithmically
spaced test grid from

.. math::

   10^{-18}\ {\rm kg}

to

.. math::

   10^{-12}\ {\rm kg}

using ``erf.sbm_nbins`` bins.

This automatically generated grid is an infrastructure-test default. It is
not a recommended atmospheric spectral discretization.

``erf.sbm_pivots``
~~~~~~~~~~~~~~~~~~

**Type:** List of real values

**Units:** kg of liquid water per particle

Defines one representative particle mass for each spectral bin.

Each pivot must be finite, positive, and lie inside its corresponding bin.

For ``erf.sbm_moment_mode = 1``, the pivots must also increase strictly from
one bin to the next. Coincident pivots at a shared bin edge are therefore not
valid for the fixed-pivot one-moment representation.

In two-moment mode the pivots remain part of the spectral-grid and restart
identity, but the reference packet remapper uses the packet's actual mass and
its containing interval rather than replacing the packet by pivot masses.

When explicit ``erf.sbm_edges`` are supplied and ``erf.sbm_pivots`` are
omitted, ERF uses the geometric mean of each pair of adjacent bin edges.

If ``erf.sbm_edges`` is omitted, ERF generates both the default edges and
their geometric-mean pivots. A separately supplied ``erf.sbm_pivots`` array is
therefore not used unless an explicit edge array is also supplied.

A special case occurs when the first explicit bin begins at exactly zero.
The geometric mean of zero and a positive upper edge is zero, but SBM requires
a strictly positive pivot. Such a grid therefore requires an explicit positive
pivot inside the first bin.

``erf.sbm_fixture_initial_state``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** List of real values

**Default:** all spectral components are zero

Provides a spatially uniform initial spectral state. Each supplied component is
assigned the same value in every grid cell.

For one-moment SBM with :math:`N` bins, supply exactly :math:`N` values:

.. math::

   M_0,\ M_1,\ldots,M_{N-1},

with units of :math:`\mathrm{kg\,m^{-3}}`.

For two-moment SBM, supply exactly :math:`2N` values in the order

.. math::

   M_0,\ldots,M_{N-1},C_0,\ldots,C_{N-1},

where the mass components have units of :math:`\mathrm{kg\,m^{-3}}` and the
number components have units of :math:`\mathrm{m^{-3}}`.

All supplied values must be finite and nonnegative. Two-moment states must
also satisfy the per-bin realizability condition described above.

Current execution restrictions
------------------------------

The SBM runtime remains single-level. Independently tested grid-transfer
operators do not yet enable mesh-refined SBM simulations. The following
runtime restrictions apply:

The current configuration requires:

* a three-dimensional ERF build;
* one AMR level, ``amr.max_level = 0``;
* periodic boundaries in all three directions;
* static non-EB geometry; end-to-end evidence covers ``ConstantDz``,
  ``StretchedDz``, and the triply periodic static ``Cos4Hill`` fitted-mesh
  smoke profile, with terrain flat at the periodic boundaries;
* no moving terrain, embedded boundaries, or immersed buildings;
* ``substepping_type = None``;
* no molecular scalar diffusion;
* no turbulent scalar diffusion or PBL/SHOC transport acting on moisture;
* no ERF numerical diffusion;
* no sponge source;
* no configured perturbation forcing;
* no custom subsidence;
* no large-scale forcing;
* no sounding nudging;
* no custom moisture forcing;
* no real-data lateral boundary forcing; and
* only the supported problem setups: the default/undefined setup,
  ``Scalar Advection/Diffusion``, and ERF's dedicated periodic SBM test setup;
  other problem-specific initialization or forcing is rejected.

Anelastic SBM advection is supported only for the currently qualified
non-fitted-terrain cases, including constant- and stretched-z grids.
``StaticFittedMesh`` with ``VariableDz`` is explicitly rejected for anelastic
SBM because ERF's fitted-terrain anelastic projection requires a nonperiodic
vertical boundary, whereas the current spectral transport capability requires
a triply periodic domain. This is a support-envelope restriction, not a
different spectral carrier convention.

The resolved host dry-air carrier may be nonzero. SBM transport uses that
carrier directly and checks the donor outgoing-demand condition at each host
stage. For example, the native scalar-advection velocity initializer can be
used in a periodic SBM advection case with:

.. code-block:: text

   erf.prob_name = "Scalar Advection/Diffusion"
   prob.U_0 = 0.8
   prob.V_0 = 0.6
   prob.W_0 = 0.4

Adaptive ERF timestep selection also applies the SBM donor-outflow restriction
using the current mapped dry-air carrier. A user-specified fixed timestep that
exceeds the current hard donor-positivity limit is rejected rather than
silently subcycling the spectrum. The exact stage-time outgoing-demand check
remains active as a fail-closed guard because the carrier used by a later host
stage may differ from the pre-step estimate.

The currently enabled SBM capability implements resolved advection only;
cloud microphysical source processes are not yet coupled to the spectrum.

One-moment example
------------------

The following configuration is based on the ERF one-moment SBM integration
test:

.. code-block:: text

   erf.prob_name = "Scalar Advection/Diffusion"
   erf.init_type = Uniform

   erf.moisture_model = SBM
   erf.sbm_nbins = 4
   erf.sbm_cloud_rain_split = 2
   erf.sbm_fixture_initial_state = 1.e-6 2.e-6 3.e-6 4.e-6

   erf.use_gravity = false
   erf.substepping_type = None
   erf.molec_diff_type = None
   erf.les_type = None
   erf.pbl_type = None
   erf.sponge_type = None

   amr.max_level = 0
   geometry.is_periodic = 1 1 1

   prob.U_0 = 0.8
   prob.V_0 = 0.6
   prob.W_0 = 0.4

The four initial-state values are the density-weighted liquid-water masses in
the four bins, in :math:`\mathrm{kg\,m^{-3}}`.

With ``erf.sbm_cloud_rain_split = 2``, the first two bins contribute to cloud
water and the last two bins contribute to rain water.

The spatially uniform initial spectrum exercises the host carrier path but is
not a recommended atmospheric size distribution.

Two-moment example
------------------

The two-moment integration fixture uses:

.. code-block:: text

   erf.prob_name = "Scalar Advection/Diffusion"
   erf.init_type = Uniform

   erf.moisture_model = SBM
   erf.sbm_nbins = 4
   erf.sbm_moment_mode = 2
   erf.sbm_cloud_rain_split = 2

   erf.sbm_edges = 1.0e-18 2.0e-18 4.0e-18 8.0e-18 1.6e-17

   erf.sbm_fixture_initial_state = \
       1.5e-6 6.0e-6 1.8e-5 4.8e-5 \
       1.0e12 2.0e12 3.0e12 4.0e12

   prob.U_0 = 0.8
   prob.V_0 = 0.6
   prob.W_0 = 0.4

The first four initial-state values are the bin liquid-water mass densities
:math:`M_i` in :math:`\mathrm{kg\,m^{-3}}`. The final four are the bin
droplet-number densities :math:`C_i` in :math:`\mathrm{m^{-3}}`.

The values have been chosen to satisfy the two-moment realizability condition
for the specified particle-mass edges. They are test values rather than a
recommended atmospheric distribution.

Checkpoint and restart
----------------------

SBM adds authoritative spectral information to the normal ERF checkpoint. See
also :ref:`sec:Checkpoint`.

An SBM checkpoint contains:

* the ordinary ERF atmospheric state;
* an ``SBM_Schema`` record describing the exact spectral layout and current
  representation contract; and
* an ``SBMSpectrum`` containing the authoritative bin-resolved liquid state.

The schema comparison is exact. A restart must use the same spectral
definition, including quantities such as the bin edges, pivots, moment mode,
and cloud/rain split.

ERF then validates the saved authoritative spectrum before allowing the
projected bulk state to be regenerated. Every spectral value must be finite,
and the runtime spectral constraints, including the two-moment realizability
conditions, must be satisfied.

The cloud- and rain-water fields saved in the ordinary ERF state are also
compared with the values obtained by projecting the saved spectrum. The
comparison uses a tight floating-point consistency tolerance.

Only after the saved spectrum and the saved bulk projection have both passed
these checks does ERF regenerate the compact cloud- and rain-water fields from
the authoritative spectrum.

A missing ``SBM_Schema``, missing ``SBMSpectrum``, incompatible spectral
schema, inadmissible spectrum, or inconsistent cloud/rain projection causes
the restart to fail. A corrupted checkpoint is not silently repaired by
projection.

Output and diagnostics
----------------------

The normal ERF moisture-output interface currently exposes the compatibility
fields

* ``qv`` for water vapor;
* ``qc`` for projected cloud water; and
* ``qrain`` for projected rain water.

For SBM, ``qc`` and ``qrain`` are obtained by summing the appropriate
liquid-mass bins. They are not independent condensed-water reservoirs.

The standard aggregate moisture diagnostics consequently have the usual warm
liquid interpretation:

.. math::

   q_t = q_v + q_c + q_r,

.. math::

   q_n = q_v + q_c,

and

.. math::

   q_p = q_r.

The current SBM implementation does not register individual spectral-bin mass
or number components as ordinary three-dimensional plotfile variables.
Two-moment droplet-number information therefore exists in the authoritative
spectral state but is not exposed as conventional ``nc`` or ``nr`` core-state
output.

The bin-resolved state is currently preserved in checkpoint
``SBMSpectrum`` data.

There are no SBM surface-precipitation accumulators or SBM-specific
microphysical tendency diagnostics in the current fixture because
sedimentation and cloud microphysical processes are not yet implemented.

See :ref:`sec:Plotfile3DReference` for the general ERF plotfile-variable
interface.
