.. _CloudChamber:

Cloud Chamber
==============

The Cloud Chamber problem provides a single-level, Cartesian, anelastic
proof-of-concept configuration for buoyancy-driven moist convection in a
closed rectangular chamber.  It supports physical temperature and
relative-humidity initialization, fixed-temperature walls, impermeable dry
walls, saturated wet walls, and equilibrium cloud water through SatAdj.

This Stage 1 configuration is not a quantitatively validated Pi-Chamber LES
and does not include an engineering wall model or finite-rate droplet
microphysics.

Supported features
------------------

===============================  =================
Feature                          Stage 1 support
===============================  =================
Rectangular Cartesian chamber    Yes
Single level                     Yes
Anelastic dynamics               Required
SatAdj equilibrium cloud         Yes
Dry vapor-impermeable walls      Yes
Wet equilibrium walls            Yes
Fixed physical wall temperature  Yes
ConstantAlpha normal transfer    Yes
Interior LES                     Yes
Terrain                          No
Embedded boundaries              No
AMR                              No
Engineering wall law             No
Roughness                        No
Finite-rate droplets             No
===============================  =================

The supported geometry is a positive rectangular prism with
``erf.mesh_type = ConstantDz``, ``erf.anelastic = 1``,
``erf.use_gravity = true``, and ``amr.max_level = 0``.  Periodic directions,
terrain, immersed boundaries, and AMR are not supported by this configuration.

Initial state
-------------

Physical initialization is selected with:

.. code-block:: none

   erf.prob_name = "Cloud Chamber"
   erf.init_type = ConstantDensity
   erf.anelastic = 1
   erf.use_gravity = true
   prob.thermodynamic_initialization = physical_temperature_rh
   prob.initial_temperature_bottom = 300.0       # K
   prob.initial_temperature_top = 284.0           # K
   prob.initial_relative_humidity = 0.95          # fraction, not percent
   prob.temperature_perturbation_amplitude = 0.02 # K
   prob.perturbation_mode = deterministic_sine

The initialized physical temperature is

.. math::

   T(z)=T_b+(T_t-T_b)\frac{z-z_\mathrm{lo}}{L_z}+\delta T.

The perturbation is deterministic, bounded by its configured amplitude, and
vanishes on the physical faces.  Relative humidity is a fraction in
``[0,1]``; entering ``95`` instead of ``0.95`` is invalid.  Vapor is initialized
from the actual temperature:

.. math::

   e_v=RH\,e_s(T), \qquad
   q_v=\frac{R_d}{R_v}\frac{e_v}{p_0-e_v}.

Here ``qv`` is a mixing ratio relative to dry air.  ``qc`` initially equals
zero.  A sub-saturated initial RH should not cause whole-domain startup
condensation.

The legacy ``legacy_theta_qv`` mode remains available for existing numerical
inputs.  It must not be mixed with the physical temperature/RH keys.

Anelastic pressure and temperature
-----------------------------------

For anelastic SatAdj, ERF uses the HSE base pressure ``p0`` and

.. math::

   T=\theta\left(\frac{p_0}{p_{00}}\right)^{R_d/c_p}.

SatAdj uses this same ``p0``; it is converted to hPa only at the saturation
helper interface.  Compressible SatAdj retains its existing EOS-derived
pressure and temperature path.

The Stage 1 surface reference is ERF's current ``100000`` Pa reference
pressure.  ``prob.p_inf`` and ``prob.T_0`` contribute to construction of the
reference density when needed; ``prob.p_inf`` does not independently redefine
the HSE pressure anchor.  A physical Stage 1 input should use
``prob.p_inf = 100000.0``.  Inputs materially inconsistent with the ERF
reference pressure are rejected.  The parser uses a relative tolerance of
``1e-6`` for this comparison.

Wall temperature and moisture
-----------------------------

Every face is a ``NoSlipWall`` with a physical temperature in kelvin and a
moisture mode.  For example:

.. code-block:: none

   zlo.type = NoSlipWall
   zlo.temperature = 300.0       # K
   zlo.moisture = wet             # dry or wet
   zlo.wall_transfer_model = resolved_molecular

``resolved_molecular`` is the public Stage 1 name for physical wall-normal
transfer using only the configured ``ConstantAlpha`` coefficients.  Set
``erf.molec_diff_type = ConstantAlpha``, ``erf.alpha_T`` and ``erf.alpha_C``
to positive scalar coefficients.  The interior SGS diffusivity may remain
active, but SGS diffusivity is excluded from the physical wall-normal
transfer.  Short validation examples may use enhanced coefficients; they are
not automatically scientifically calibrated molecular properties.

Dry walls
~~~~~~~~~

A dry wall is impermeable to water vapor and cloud water.  It does not impose
zero humidity:

.. math::

   J_v=0, \qquad J_c=0.

Wet walls
~~~~~~~~~

A wet wall represents a stationary pure-liquid surface at saturation at the
configured wall temperature.  Stage 1 uses the adjacent fluid-cell base
pressure for the wall saturation state:

.. math::

   q_{v,w}=q_\mathrm{sat}(T_w,p_{0,c}).

Here ``p0,c`` is the adjacent fluid-cell base pressure.  This approximation
retains height dependence on vertical walls.  Evaporation and condensation
are both permitted: positive inward flux adds vapor to the fluid, while a
negative inward flux represents condensation onto the wall.  Cloud water does
not cross a wall in Stage 1.

SatAdj limitations
------------------

SatAdj is an equilibrium qv/qc partition.  It conserves total
non-precipitating water and applies the existing latent-heating adjustment.
It does not represent droplet number or size distribution and does not
include activation, rain, sedimentation, collision-coalescence, precipitation,
or cloud-water wall deposition.

Budgets and diagnostic neutrality
---------------------------------

Set ``erf.cloud_chamber_budget_interval`` to a positive step interval to write
``cloud_chamber_budget.dat``.  A non-positive value disables reporting.  Each
report is interval-local and contains the start and end step/time, state
change, six face contributions, net boundary contribution, internal source,
residual, tolerance, and status.  Face fluxes use coordinate orientation with
low-minus-high signs.  Supported statuses are ``PASS``, ``FAIL``, and
``UNSUPPORTED_SOURCE``.

``rhoTheta`` is a conserved potential-temperature scalar, not a heat rate.
For total non-precipitating water,

.. math::

   \Delta\int_\Omega\rho(q_v+q_c)\,dV
   =\int_{t_a}^{t_b}\sum_f J_{v,\mathrm{in}}\,dA\,dt.

For an all-dry chamber the right-hand side is zero.  SatAdj's internal qv/qc
exchange cancels in the total-water row.  The closure tolerance is computed
from double- or single-precision machine epsilon, an effective cell/face
operation count, the state scale, and the absolute boundary and internal
contributions; it is printed with every row.

Enabling budgets must not change the simulation state.  Budget accumulation,
MPI reduction, reporting, and file output are observational only.  The
registered parity test compares density, theta, temperature, qv, qc, all three
velocities, and the corresponding global integrals with budgets enabled and
disabled.

Advanced scalar-advection audit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set ``erf.cloud_chamber_water_ledger = 1`` to enable the opt-in diagnostic
ledger used by the conservation regression.  It records the qv/qc right-hand
side change, an independent reduction of that change, directional divergence
terms, boundary terms, topology, and an interface-consistency audit.  Rows are
step summaries; the advection diagnostics are accumulated over the reported
RK stages, while the interface audit retains per-stage shared and mismatched
face counts.  The current stage summaries are unweighted accumulations for
diagnosis, not exact RK-weighted timestep integrals; they must not be used as
production conservation claims.  In a closed all-dry chamber, the ledger should report
roundoff-level total-water closure and zero mismatched faces.  Shared faces
are not errors when their duplicate values agree.  The ledger is diagnostic
only and does not alter the solution.  The grid-owner audit is currently
available in serial; on the current distributed AMReX/MPI configuration the
MPI interface audit is disabled and is reported as unavailable rather than
being treated as a pass.

The ledger is disabled by default, writes ``cloud_chamber_water_ledger.dat``,
and can add substantial runtime and memory overhead.  It is intended for
developer-level numerical diagnosis within the supported Stage 1 Cartesian
single-level geometry, not for routine scientific output or quantitative
Pi-Chamber validation.

Complete examples
-----------------

No-moisture thermal chamber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``Tests/test_files/CloudChamber_Dry/CloudChamber_Dry.i``.  This case
exercises physical-temperature initialization, six dry walls, and thermal
wall forcing without a moisture model.

SatAdj all-dry closed-water chamber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following settings form the registered multi-box conservation fixture:

.. code-block:: none

   erf.moisture_model = SatAdj
   erf.cloud_chamber_budget_interval = 2
   erf.molec_diff_type = ConstantAlpha
   amr.n_cell = 16 16 16
   amr.max_grid_size = 8
   max_step = 6
   xlo.moisture = dry
   xhi.moisture = dry
   ylo.moisture = dry
   yhi.moisture = dry
   zlo.moisture = dry
   zhi.moisture = dry

The ``CloudChamber_SatAdj_AllDry`` CTest case requires at least three passing
budget intervals and zero total-water change within the printed tolerance.

Wet-top/wet-bottom proof of concept
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The registered ``CloudChamber_SatAdj_WetBudget`` fixture uses the same
16-by-16-by-16 multi-box layout with:

.. code-block:: none

   zlo.temperature = 300.0
   zlo.moisture = wet
   zhi.temperature = 284.0
   zhi.moisture = wet
   xlo.moisture = dry
   xhi.moisture = dry
   ylo.moisture = dry
   yhi.moisture = dry

It checks that wet-wall vapor fluxes close the total-water change, side-wall
vapor fluxes are zero, cloud-water wall fluxes are zero, and three consecutive
intervals pass.

Validation guidance
-------------------

Inspect ``pressure``, ``temp``, ``theta``, ``qv``, ``qc``, ``qsat``,
``rel_humidity``, and the velocity fields in plotfiles.  At initialization,
velocity and ``qc`` should be zero, RH should be near the configured value,
and the physical temperature—not merely potential temperature—should match
the requested profile.  A short thermal perturbation should then produce a
bounded buoyant response.

For a serious run, record temperature, theta, qv, qc, qv+qc, integrated
``rhoTheta``, integrated vapor, integrated cloud water, total water, velocity
extrema, timestep, wall fluxes, and budget residuals at every output.  Stability
alone is not quantitative Pi-Chamber validation.

Testing coverage
----------------

The registered Stage 1 tests have distinct purposes:

* ``CloudChamberFluxStorage.TiledComponentsPreserveEarlierTiles`` is a
  deterministic tile-lifecycle negative control for shared-FAB resets.
* ``CloudChamberFaceAudit.NegativeControlDetectsPrescribedMismatch`` verifies
  that a prescribed face perturbation is detected while matching copies remain
  valid.
* ``CloudChamberWallFlux.MultiBoxOwnershipAcrossAllFaces`` covers all six
  orientations and multiple boxes.
* ``CloudChamber_SatAdj_Parity`` runs paired budget-off/on simulations and
  checks fields and global integrals.
* ``CloudChamber_SatAdj_AllDry`` checks closed total water with six dry walls.
* ``CloudChamber_SatAdj_WetBudget`` checks wet-wall water closure and three
  reporting intervals.
* ``CloudChamber_SatAdj_OpenMP`` is enabled when ERF is built with OpenMP and
  compares one- and two-thread results.
* ``CloudChamber_SatAdj_UntiledContainment`` records and checks the actual
  multi-grid topology used by the broad moisture containment; it does not
  claim to exercise tiled moisture advection.

These are decomposition and operator regressions, not quantitative LES or
Pi-Chamber calibration tests.  GPU, sanitizer, single-precision, and extended
32/64-cubed scientific runs must be reported separately when they are run.

Troubleshooting
---------------

* If RH is entered as ``95``, change it to ``0.95``.
* Do not mix ``physical_temperature_rh`` keys with legacy ``theta`` or ``qv``
  wall and profile keys.
* If an all-dry budget changes total water, inspect the six qv/qc face rows,
  their tolerances, and the total-water residual before interpreting the flow.
* If budgets change the solution, treat the run as invalid and inspect the
  paired parity comparison.
* A wet wall may condense as well as evaporate; its signed flux need not always
  add vapor.
* Unsupported diffusion, terrain, periodic, AMR, and embedded-boundary
  combinations are rejected by the Cloud Chamber parser.
