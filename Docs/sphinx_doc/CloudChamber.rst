.. index:: Cloud Chamber
.. index:: SatAdj
.. index:: ConstantAlpha
.. index:: NoSlipWall
.. index:: cloud_chamber_budget_interval
.. index:: qv
.. index:: qc
.. index:: rhoTheta

.. _CloudChamber:

Cloud Chamber
=============

The Cloud Chamber problem provides a single-level, Cartesian, anelastic
proof-of-concept configuration for buoyancy-driven thermal and moist
convection in a closed rectangular chamber.

Dry mode evolves thermal convection without moisture state variables.
Moist mode uses SatAdj for instantaneous equilibrium partitioning between
water vapor and cloud water.

This Stage 1 configuration is not a quantitatively validated Pi-Chamber LES.
It is intended for code-path verification, qualitative buoyancy-driven flow,
wall-flux testing, and conserved-scalar regression tests.  It does not provide
finite-rate droplet microphysics, a calibrated engineering wall law,
grid-independent LES validation, or quantitative experimental calibration.

Choose a mode
-------------

.. list-table:: Supported Stage 1 modes
   :header-rows: 1
   :widths: 23 25 20 32

   * - Mode
     - Moisture settings
     - Wall moisture
     - Intended use
   * - Dry thermal chamber
     - Omit ``erf.moisture_model`` and
       ``prob.initial_relative_humidity``
     - All faces ``dry``
     - Thermal convection and strict ``rhoTheta`` closure
   * - SatAdj with dry walls
     - Set ``erf.moisture_model = SatAdj`` and provide RH
     - All faces ``dry``
     - Closed total-water conservation
   * - SatAdj with wet walls
     - Set ``erf.moisture_model = SatAdj`` and provide RH
     - One or more faces may be ``wet``
     - Vapor exchange with total-water closure

A wet wall is valid only with SatAdj.  A dry wall is impermeable to water,
but it may still exchange heat through its prescribed temperature.

Terminology
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Term
     - Meaning in this configuration
   * - Fixed-density anelastic
     - ERF evolves velocity and thermodynamic scalars against a prescribed
       hydrostatic reference density rather than advancing compressible density.
   * - SatAdj
     - Instantaneous equilibrium partitioning between vapor and cloud water,
       including ERF's existing latent-heating adjustment.
   * - ``qv``
     - Water-vapor mixing ratio relative to dry air.
   * - ``qc``
     - Cloud-water mixing ratio relative to dry air.
   * - ``rhoTheta``
     - Density-weighted potential temperature, the conserved thermal scalar
       used by this ERF configuration.
   * - ``ConstantAlpha``
     - User-configured thermal and scalar diffusivities used by the resolved
       wall-normal transfer.
   * - Resolved molecular transfer
     - Transfer computed from the adjacent-cell-to-wall half-cell gradient;
       it is not an engineering wall-function model.

Required configuration
----------------------

.. important::

   Stage 1 uses one uniform Cartesian mesh level with no AMR refinement.
   ERF may still divide that level into multiple boxes for parallel
   execution.  The domain is closed in all three directions, gravity and
   fixed-density anelastic dynamics are enabled, and all six boundaries are
   stationary ``NoSlipWall`` faces with prescribed temperatures.

   * ``amr.max_level = 0`` disables AMR refinement.
   * ``erf.mesh_type = ConstantDz`` selects the supported Cartesian mesh.
   * ``geometry.is_periodic = 0 0 0`` closes every direction rather than
     wrapping it periodically.
   * ``erf.init_type = ConstantDensity``, ``erf.anelastic = 1``, and
     ``erf.use_gravity = true`` select the fixed-density anelastic setup
     with gravity.
   * ``erf.vert_implicit = false`` selects the explicit vertical-diffusion
     path supported by the Stage 1 wall treatment.
   * ``erf.terrain_type`` and ``erf.buildings_type`` must be omitted or set
     to ``None``; terrain, embedded boundaries, and immersed buildings are
     outside Stage 1.
   * All six faces must be stationary ``NoSlipWall`` boundaries with an
     explicit temperature.
   * Wall-normal transfer uses
     ``wall_transfer_model = resolved_molecular`` with
     ``erf.molec_diff_type = ConstantAlpha`` and the configured
     ``alpha_T`` and ``alpha_C`` coefficients.

The common solver settings are:

.. code-block:: none

   erf.prob_name = "Cloud Chamber"
   erf.init_type = ConstantDensity
   erf.anelastic = 1
   erf.use_gravity = true
   erf.mesh_type = ConstantDz
   erf.vert_implicit = false
   erf.molec_diff_type = ConstantAlpha

   amr.max_level = 0
   geometry.is_periodic = 0 0 0

   prob.thermodynamic_initialization = physical_temperature_rh
   prob.initial_temperature_bottom = 300.0
   prob.initial_temperature_top = 284.0
   prob.temperature_perturbation_amplitude = 0.02
   prob.perturbation_mode = deterministic_sine

Despite the input-mode name ``physical_temperature_rh``, dry mode does not
use relative humidity.  For dry mode, omit both
``erf.moisture_model`` and ``prob.initial_relative_humidity``.

For SatAdj mode add:

.. code-block:: none

   erf.moisture_model = SatAdj
   prob.initial_relative_humidity = 0.95

Relative humidity is a fraction from 0 to 1, not a percentage.

Runnable examples
-----------------

The following tracked regression inputs are complete runnable examples.

.. note::

   These are complete regression inputs intended to exercise the supported
   code paths in short tests.  Their grid, timestep, LES settings, and
   transfer coefficients are not calibrated recommendations for a physical
   chamber simulation.  Interior LES is optional and does not replace the
   configured wall-normal ``alpha_T`` and ``alpha_C`` transfer.

Dry thermal chamber
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../Tests/test_files/CloudChamber_Dry/CloudChamber_Dry.i
   :language: none
   :caption: Dry thermal chamber with strict rhoTheta closure

SatAdj wet-wall chamber
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../Tests/test_files/CloudChamber_SatAdj_WetBudget/CloudChamber_SatAdj_WetBudget.i
   :language: none
   :caption: SatAdj chamber with wet lower and upper walls

The shipped fixtures use a 300 K lower wall, a 284 K upper wall, and four
292 K sidewalls.  All six walls are temperature-controlled.  The warm-lower
and cool-upper arrangement should produce an unstable thermal profile and a
bounded buoyancy-driven overturning response.

Scientific scope
----------------

This configuration can support qualitative buoyancy-driven flow, equilibrium
vapor--cloud-water thermodynamics, wall-flux verification, and conservation
testing.

It does not represent droplet activation, number or size distributions,
rain, sedimentation, collision--coalescence, precipitation, or cloud-water
deposition at a wall.  It must not be used by itself to claim quantitative
Pi-Chamber agreement or grid-independent LES behavior.

Initial thermodynamic state
---------------------------

The initializer prescribes physical temperature and then derives potential
temperature from the local hydrostatic base-state pressure:

.. math::

   \theta =
   T\left(\frac{p_\mathrm{ref}}{p_\mathrm{hse}}\right)^{R_d/c_p}.

Here:

* :math:`T` is physical temperature in kelvin;
* :math:`p_\mathrm{hse}` is the local hydrostatic base-state pressure;
* :math:`p_\mathrm{ref}` is ERF's fixed reference pressure, currently
  100000 Pa;
* :math:`R_d` is the dry-air gas constant; and
* :math:`c_p` is dry-air specific heat at constant pressure.

For physical SatAdj initialization:

.. math::

   q_v =
   \frac{R_d}{R_v}
   \frac{RH\,e_s(T)}
        {p_\mathrm{hse}-RH\,e_s(T)},

where :math:`R_v` is the water-vapor gas constant and :math:`e_s(T)` is
saturation vapor pressure.

``qv`` is a mixing ratio relative to dry air, not specific humidity.
Initialization sets ``qc = 0``; SatAdj then establishes the vapor--cloud-water
partition during model evolution.  Dry mode has no ``qv`` or ``qc`` state
fields.

Reference density and pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``prob.p_inf`` is supplied, it must agree with the 100000 Pa ERF
reference pressure within the accepted tolerance.  If ``prob.rho_0`` is
omitted and ``prob.p_inf`` is supplied, ERF computes the constant reference
density from ``prob.p_inf`` and ``prob.T_0``.  These inputs do not redefine
:math:`p_\mathrm{ref}` in the potential-temperature relationship.

The legacy ``legacy_theta_qv`` initialization remains available for existing
numerical inputs.  Do not combine legacy profile keys with physical-mode
temperature or RH keys.

Wall temperature and moisture
-----------------------------

Every physical face requires:

.. code-block:: none

   <face>.type = NoSlipWall
   <face>.temperature = <temperature in K>
   <face>.moisture = dry|wet
   <face>.wall_transfer_model = resolved_molecular

The six face prefixes are ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``,
and ``zhi``.

.. list-table:: Wall transfer contract
   :header-rows: 1
   :widths: 14 25 28 20 13

   * - Wall mode
     - Thermal transfer
     - Vapor transfer
     - Cloud-water transfer
     - Allowed modes
   * - ``dry``
     - Controlled by ``alpha_T``
     - Exactly zero
     - Exactly zero
     - Dry or SatAdj
   * - ``wet``
     - Controlled by ``alpha_T``
     - Saturation-based through ``alpha_C``; evaporation or condensation
     - Exactly zero
     - SatAdj only

A wet wall uses saturation at the configured wall temperature and the
adjacent fluid-cell HSE pressure.  The wall-normal flux is computed using the
half-cell distance between the cell center and the wall.

.. warning::

   Treat ``alpha_T`` and ``alpha_C`` as prescribed transfer coefficients,
   not calibrated molecular-property values or engineering wall-law
   coefficients.  Positive values activate the corresponding transfer and
   zero disables it.  The current parser does not enforce coefficient signs.

Flux orientation
~~~~~~~~~~~~~~~~

A positive stored face flux is positive in the coordinate direction, not
necessarily into the chamber.  Low-face fluxes point into the domain when
positive; high-face fluxes point into the domain when negative.  The reported
net boundary contribution therefore uses low minus high for each axis.

Conserved-scalar budgets
------------------------

Set ``erf.cloud_chamber_budget_interval`` to a positive step interval to write
``cloud_chamber_budget.dat``.  Each report is local to the interval since the
previous report and contains:

* six physical-face contributions;
* net boundary contribution;
* state change;
* internal source;
* residual;
* tolerance; and
* status.

``rhoTheta`` is a conserved potential-temperature scalar, not heat energy or
a heat rate.

For total nonprecipitating water:

.. math::

   \Delta\int_\Omega \rho(q_v+q_c)\,dV
   =
   \int_{t_a}^{t_b}\sum_f J_{v,\mathrm{in}}\,dA\,dt.

For six dry walls, the boundary contribution is zero.  SatAdj's internal
exchange between ``qv`` and ``qc`` cancels in the total-water row.

``PASS``
   The supported row closes within its printed tolerance.

``FAIL``
   A supported closure contract is violated.

``UNSUPPORTED_SOURCE``
   For cloudy ``rhoTheta``, the current report omits part of the moist
   latent-heating source.  Do not interpret this status as closure.  Strict
   integrated thermal closure is tested in dry mode; moist conservation is
   assessed using total nonprecipitating water.

Budget reporting must not change the simulated state.  The registered parity
test compares budget-disabled and budget-enabled solutions.

Run checklist
-------------

1. Choose dry or SatAdj mode.
2. Configure one nonperiodic Cartesian mesh level; that level may contain
   multiple AMReX boxes.
3. Define all six ``NoSlipWall`` faces and their temperatures.
4. For SatAdj, provide RH as a fraction and choose dry or wet moisture walls.
5. Set ``alpha_T`` and, when vapor transfer is needed, ``alpha_C``.
6. Run a short case and inspect temperature, potential temperature, velocity,
   and, for SatAdj, ``qv``, ``qc``, saturation mixing ratio, and RH.
7. Enable ``erf.cloud_chamber_budget_interval``.
8. Require dry thermal closure or total-water closure, as appropriate.
9. Treat any ``FAIL`` or budget-dependent solution change as invalid.

Stage 1 invariants
------------------

* One uniform Cartesian mesh level and no periodic direction; the level may
  be decomposed into multiple AMReX boxes.
* Six stationary ``NoSlipWall`` faces.
* Dry-wall vapor flux is exactly zero.
* Cloud-water wall flux is exactly zero.
* Wet walls are permitted only with SatAdj.
* Six dry walls conserve total nonprecipitating water.
* Enabling budget output does not change the solution.
* A stable run alone is not quantitative Pi-Chamber validation.

Troubleshooting
---------------

* Enter RH as ``0.95``, not ``95``.
* Do not provide RH in dry mode.
* Do not mix physical and legacy initialization keys.
* A wall marked ``dry`` may still exchange heat; ``dry`` describes water
  permeability.
* A wet wall may evaporate or condense, so its signed vapor flux need not
  always add water.
* Interpret face signs using coordinate orientation before comparing wall
  gain or loss.
* Treat ``UNSUPPORTED_SOURCE`` as incomplete cloudy thermal accounting, not
  as a successful budget result.
