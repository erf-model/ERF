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

The purpose of spectral-bin microphysics is ultimately to resolve the
hydrometeor size distribution and its evolution through physical processes
rather than representing cloud and precipitation only through a few bulk
categories. The current ERF implementation establishes the state,
configuration, ownership, projection, and restart infrastructure needed for
that capability.

.. warning::

   The current ``SBM`` option is an infrastructure qualification mode, not yet
   a production cloud-microphysics scheme.

   ERF currently stores and validates a liquid spectral distribution and keeps
   the conventional bulk cloud- and rain-water fields consistent with that
   distribution. The spectral state is intentionally not advected, diffused,
   sedimented, or modified by cloud microphysical processes.

   Spectral advection, turbulent or molecular diffusion, condensation and
   evaporation, aerosol activation, aerosol evolution, collision-coalescence,
   sedimentation, precipitation, AMR spectral transport, and physical
   spectral boundary conditions are not yet available through
   ``erf.moisture_model = SBM``.

   Unsupported configurations are rejected rather than silently reverting to
   bulk-water transport.

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

SBM runtime inputs
------------------

The following inputs use the ``erf.`` prefix.

``erf.moisture_model``
~~~~~~~~~~~~~~~~~~~~~~

Select the current SBM infrastructure with

::

   erf.moisture_model = SBM

At the present revision this also requires

::

   erf.sbm_zero_transport_fixture = true

because production spectral transport has not yet been implemented.

``erf.sbm_zero_transport_fixture``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** Boolean

**Default:** ``false``

This option must currently be ``true`` when ``erf.moisture_model = SBM``.

It identifies the deliberately bounded zero-transport SBM configuration. It
should not be interpreted as a physical switch that turns transport off in an
otherwise complete spectral-bin microphysics scheme.

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

The current SBM mode is intentionally limited so that state ownership,
projection, initialization, and restart can be tested without ambiguity from
an incomplete transport implementation.

The current configuration requires:

* a three-dimensional ERF build;
* one AMR level, ``amr.max_level = 0``;
* periodic boundaries in all three directions;
* ``mesh_type = ConstantDz``;
* no terrain-fitted or embedded-boundary geometry;
* no immersed buildings;
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
* no problem setup that introduces problem-specific liquid forcing or custom
  perturbations.

For the present fixture, ``erf.prob_name`` must therefore be either unset
(``Undefined`` internally) or

::

   erf.prob_name = "SBM zero-transport fixture"

In addition, the resolved carrier momentum presented to scalar advection must
remain exactly zero. A nonzero resolved flow causes ERF to stop rather than
leave the spectral state behind while transporting only its bulk projection.

These restrictions make the present capability suitable for infrastructure
qualification, not for a moving or evolving cloud simulation.

One-moment example
------------------

The following configuration is based on the ERF one-moment SBM integration
test:

.. code-block:: text

   erf.prob_name = "SBM zero-transport fixture"
   erf.init_type = Uniform

   erf.moisture_model = SBM
   erf.sbm_zero_transport_fixture = true
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

The four initial-state values are the density-weighted liquid-water masses in
the four bins, in :math:`\mathrm{kg\,m^{-3}}`.

With ``erf.sbm_cloud_rain_split = 2``, the first two bins contribute to cloud
water and the last two bins contribute to rain water.

This example is a qualification fixture and should not be interpreted as a
recommended atmospheric size distribution.

Two-moment example
------------------

The two-moment integration fixture uses:

.. code-block:: text

   erf.prob_name = "SBM zero-transport fixture"
   erf.init_type = Uniform

   erf.moisture_model = SBM
   erf.sbm_zero_transport_fixture = true
   erf.sbm_nbins = 4
   erf.sbm_moment_mode = 2
   erf.sbm_cloud_rain_split = 2

   erf.sbm_edges = 1.0e-18 2.0e-18 4.0e-18 8.0e-18 1.6e-17

   erf.sbm_fixture_initial_state = \
       1.5e-6 6.0e-6 1.8e-5 4.8e-5 \
       1.0e12 2.0e12 3.0e12 4.0e12

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
