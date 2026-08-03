.. _CloudChamber:

Cloud Chamber
=============

The Cloud Chamber problem provides a single-level, Cartesian, anelastic
proof-of-concept configuration for buoyancy-driven thermal and moist
convection in a closed rectangular chamber.  Dry mode evolves thermal
convection without a moisture model; moist mode uses SatAdj for equilibrium
vapor/cloud-water partitioning.

This Stage 1 configuration is not yet a quantitatively validated Pi-Chamber
LES.  It does not provide an engineering wall law, finite-rate droplets, or
quantitative experimental calibration.

Overview and Stage 1 scope
--------------------------

The supported problem is a positive rectangular Cartesian prism with one
level, ``ConstantDz`` spacing, fixed-density anelastic dynamics, gravity, and
six physical ``NoSlipWall`` faces.  Periodicity, terrain, immersed buildings,
and AMR are outside the Stage 1 scope.

Supported configuration
-----------------------

.. list-table:: Cloud Chamber Stage 1 support
   :header-rows: 1
   :widths: 42 58

   * - Capability
     - Support
   * - No-moisture thermal mode
     - Yes
   * - SatAdj equilibrium moisture
     - Yes
   * - Interior LES
     - Optional; wall transfer remains resolved molecular
   * - ConstantAlpha wall transfer
     - Required for physical wall transfer
   * - Terrain
     - No
   * - Embedded boundaries
     - No
   * - AMR
     - No
   * - Engineering wall law
     - No
   * - Finite-rate droplets
     - No

Required physical-mode settings
--------------------------------

The following solver and geometry settings select the supported physical
configuration.  ``erf.terrain_type`` and ``erf.buildings_type`` may be
omitted when they default to ``None``; they must not select terrain or
immersed buildings.

.. code-block:: none

   erf.prob_name = "Cloud Chamber"
   erf.init_type = ConstantDensity
   erf.anelastic = 1
   erf.use_gravity = true
   erf.mesh_type = ConstantDz
   erf.vert_implicit = false
   erf.molec_diff_type = ConstantAlpha
   erf.terrain_type = None
   erf.buildings_type = None
   amr.max_level = 0
   geometry.is_periodic = 0 0 0

The physical initializer requires a temperature profile:

.. code-block:: none

   prob.thermodynamic_initialization = physical_temperature_rh
   prob.initial_temperature_bottom = 300.0
   prob.initial_temperature_top = 284.0
   prob.temperature_perturbation_amplitude = 0.02
   prob.perturbation_mode = deterministic_sine

Dry mode omits both ``erf.moisture_model`` and
``prob.initial_relative_humidity``.  SatAdj mode adds:

.. code-block:: none

   erf.moisture_model = SatAdj
   prob.initial_relative_humidity = 0.95

Relative humidity is a fraction in ``[0,1]``, not a percentage.  In physical
SatAdj mode the RH key is required and is validated.  In dry physical mode it
has no effect and an explicitly supplied key is rejected with a diagnostic.

Initial thermodynamic state
----------------------------

The initializer prescribes the physical-temperature profile, including its
bounded deterministic perturbation.  It then derives potential temperature
from the local hydrostatic base-state pressure.  Define:

* :math:`p_\mathrm{hse}` as the local hydrostatic base-state pressure; and
* :math:`p_\mathrm{ref}` as ERF's fixed reference pressure, currently
  ``100000`` Pa.

The physical-temperature and potential-temperature relationships are

.. math::

   \theta = T\left(\frac{p_\mathrm{ref}}{p_\mathrm{hse}}\right)^{R_d/c_p},
   \qquad
   T = \theta\left(\frac{p_\mathrm{hse}}{p_\mathrm{ref}}\right)^{R_d/c_p}.

For physical SatAdj initialization, vapor is a dry-air mixing ratio computed
using the same local HSE pressure:

.. math::

   q_v = \frac{R_d}{R_v}
         \frac{RH\,e_s(T)}{p_\mathrm{hse}-RH\,e_s(T)}.

Here ``qv`` is vapor mixing ratio and ``qc`` is cloud-water mixing ratio.  Dry
mode does not create ``qv`` or ``qc`` fields.  SatAdj physical initialization
sets ``qc`` to zero before SatAdj performs its equilibrium adjustment.  The
legacy ``legacy_theta_qv`` mode remains available for existing numerical
inputs and should not be mixed with physical-temperature keys.

Wall temperature and moisture
-----------------------------

Every physical face must specify a no-slip wall, a temperature in kelvin, a
moisture mode, and the resolved transfer model:

.. code-block:: none

   <face>.type = NoSlipWall
   <face>.temperature = <Kelvin>
   <face>.moisture = dry|wet
   <face>.wall_transfer_model = resolved_molecular

The six faces are ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, and ``zhi``.
Wet walls are supported only with SatAdj.

Dry walls
~~~~~~~~~

A dry wall is impermeable.  Vapor and cloud-water wall fluxes are set exactly
to zero; dry moisture does not impose a zero-humidity Dirichlet state.

Wet walls
~~~~~~~~~

A wet wall uses saturation at the configured wall temperature and the
adjacent fluid-cell HSE pressure.  Both evaporation and condensation are
permitted.  Cloud water does not cross a wall.  The wall-normal transfer uses
a half-cell gradient, so the stored face flux is coordinate-oriented and its
sign depends on the face direction.

The resolved wall kernel uses ``alpha_T`` directly for thermal transfer and
``alpha_C`` directly for vapor transfer at wet walls.  A positive coefficient
activates the corresponding transfer; zero disables it.  The parser currently
requires ``ConstantAlpha`` but does not enforce a positive-coefficient policy.
Dry vapor and all cloud-water wall fluxes remain zero regardless of
``alpha_C``.  These are configured diffusivities, not automatically calibrated
engineering or molecular wall laws.  Interior SGS diffusivity does not replace
the configured resolved wall-normal coefficient.

SatAdj scope and limitations
----------------------------

SatAdj performs an equilibrium vapor/cloud-water partition and retains ERF's
existing latent-heating adjustment.  It does not represent droplet number or
size distribution and does not include activation, rain, sedimentation,
collision-coalescence, precipitation, or cloud-water wall deposition.

Conserved-scalar budgets
------------------------

Set ``erf.cloud_chamber_budget_interval`` to a positive step interval to write
``cloud_chamber_budget.dat``.  Reporting is interval-local and includes the
six face contributions, state change, internal source, residual, tolerance,
and status.  Stored face fluxes are coordinate-oriented; the net boundary
contribution is low minus high for each coordinate pair:

.. math::

   J_{x,\mathrm{net}} = J_{x,\mathrm{lo}} - J_{x,\mathrm{hi}},
   \quad\text{and similarly for } y \text{ and } z.

For total nonprecipitating water, the closed-chamber contract is

.. math::

   \Delta\int_\Omega \rho(q_v+q_c)\,dV
   = \int_{t_a}^{t_b}\sum_f J_{v,\mathrm{in}}\,dA\,dt.

``PASS`` means that a supported row closes within its printed tolerance.
``FAIL`` means that a supported closure contract is violated.
``UNSUPPORTED_SOURCE`` is emitted only for a cloudy ``rhoTheta`` row whose
current budget does not fully represent the moist latent-heating contribution;
it is not a passing result.  Strict integrated ``rhoTheta`` closure is enforced
by the dry chamber regression, while moist conservation is assessed using
total nonprecipitating water.  SatAdj's internal qv/qc exchange cancels in that
total-water row.

Budget accumulation, MPI reduction, reporting, and file output are
observational.  The registered parity test verifies that enabling budgets does
not change the simulated state.

Example configurations
----------------------

The shipped inputs are complete runnable fixtures; the snippets above are only
configuration excerpts.

* ``CloudChamber_Dry_ThermalBudget`` demonstrates dry physical-temperature
  initialization, six dry walls, and strict thermal closure over three
  reporting intervals.
* ``CloudChamber_SatAdj_AllDry`` demonstrates closed total-water conservation
  with SatAdj and six dry walls.
* ``CloudChamber_SatAdj_WetBudget`` demonstrates wet top/bottom water exchange,
  zero side-wall vapor flux, zero cloud-water wall flux, and total-water
  closure.

Validation guidance
-------------------

At initialization, inspect physical temperature as well as potential
temperature.  Dry plotfiles should contain no moisture fields; SatAdj
plotfiles should begin with finite vapor, zero cloud water before adjustment,
and the requested RH diagnostics.  The registered tests cover six-face wall
ownership, dry thermal closure, diagnostic neutrality, closed-water
conservation, wet-wall water transfer, and OpenMP parity.  These are Stage 1
operator and decomposition regressions, not quantitative LES or Pi-Chamber
calibration.

Troubleshooting
---------------

* Use RH as a fraction, such as ``0.95``, rather than ``95``.
* Keep dry mode free of ``erf.moisture_model`` and
  ``prob.initial_relative_humidity``.
* Use ``erf.moisture_model = SatAdj`` and a valid RH for physical moist mode.
* If an all-dry budget changes total water, inspect all six qv/qc face rows and
  the total-water residual before interpreting the flow.
* If budgets change the solution, treat the run as invalid and inspect the
  paired budget-off/budget-on comparison.
* A wet wall may condense as well as evaporate; its signed flux need not always
  add vapor.
