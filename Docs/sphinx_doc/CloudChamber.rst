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
   * The retained Stage 1 aggregate wall key
     ``wall_transfer_model = resolved_molecular`` remains supported.  The
     generalized contract may instead select ``bulk_aero`` independently for
     heat and vapor with fixed, nonnegative coefficients.

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

Dry bulk-aerodynamic wall chamber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../Tests/test_files/CloudChamber_Dry_BulkMixed/CloudChamber_Dry_BulkMixed.i
   :language: none
   :caption: Multi-box dry chamber exercising the bulk heat-wall path

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

The six face prefixes are ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``,
and ``zhi``.

The aggregate key may be used for the retained resolved path:

.. code-block:: none

   <face>.wall_transfer_model = resolved_molecular

For the generalized path, do not combine the aggregate key with channel keys.
The per-channel model selectors are:

.. code-block:: none

   <face>.momentum_transfer_model = resolved_noslip
   <face>.heat_transfer_model = resolved_molecular|bulk_aero
   <face>.vapor_transfer_model = resolved_molecular|bulk_aero

If at least one scalar channel uses ``bulk_aero``, also set:

.. code-block:: none

   <face>.coefficient_source = fixed

Provide ``C_H`` only when heat uses ``bulk_aero`` and ``C_E`` only when
vapor uses ``bulk_aero``.  For example, a face using bulk heat and resolved
wet vapor requires:

.. code-block:: none

   <face>.heat_transfer_model = bulk_aero
   <face>.vapor_transfer_model = resolved_molecular
   <face>.coefficient_source = fixed
   <face>.C_H = 0.1

Scalar wall-transfer equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Cloud Chamber wall evaluator defines physical scalar fluxes as positive
**into the chamber fluid**.  Let :math:`\rho` be the adjacent-cell density,
:math:`T_w` the prescribed wall temperature, :math:`p_\mathrm{hse}` the
adjacent hydrostatic base-state pressure, and :math:`\Delta n` the cell width
normal to the wall.  ERF potential temperature uses

.. math::

   \Pi = \left(\frac{p_\mathrm{hse}}{p_\mathrm{ref}}\right)^{R_d/c_p},
   \qquad
   \theta_w = \frac{T_w}{\Pi}
   = T_w\left(\frac{p_\mathrm{ref}}{p_\mathrm{hse}\right)^{R_d/c_p}.

For the stationary Stage 1 wall, the bulk model uses the local tangential
relative velocity reconstructed from ERF's staggered velocity fields,

.. math::

   \mathbf{u}_t =
   (\mathbf{I}-\mathbf{n}\mathbf{n}^T)(\mathbf{u}_c-\mathbf{u}_w),
   \qquad
   U_t = |\mathbf{u}_t|,
   \qquad
   \mathbf{u}_w = 0.

The normal velocity is removed before :math:`U_t` is formed.

For ``heat_transfer_model = resolved_molecular``, the existing half-cell
resolved transfer is

.. math::

   J_{\rho\theta,\mathrm{in}}
   = \rho\,\alpha_T\,(\theta_w-\theta_a)\frac{2}{\Delta n}.

For ``heat_transfer_model = bulk_aero`` with
``coefficient_source = fixed``, ERF applies

.. math::

   J_{\rho\theta,\mathrm{in}}
   = \rho\,C_H\,U_t\,(\theta_w-\theta_a).

``C_H`` is dimensionless.  The equivalent inward sensible-heat flux is

.. math::

   H_\mathrm{in} = \rho c_p C_H U_t(T_w-T_a),

and, because :math:`T=\Pi\theta`, the conserved potential-temperature flux
satisfies

.. math::

   J_{\rho\theta,\mathrm{in}}
   = \frac{H_\mathrm{in}}{c_p\Pi}.

The bulk heat model therefore replaces the physical boundary flux of ERF's
conserved ``rhoTheta`` scalar; it is not inserted as an independent energy
source term.

Cloud Chamber ``qv`` and ``qc`` are mixing ratios relative to dry air, not
specific humidity.  For a wet wall with
``vapor_transfer_model = resolved_molecular``,

.. math::

   J_{\rho q_v,\mathrm{in}}
   = \rho\,\alpha_C\,
     [q_\mathrm{sat}(T_w,p_\mathrm{hse})-q_{v,a}]
     \frac{2}{\Delta n}.

For a wet wall with ``vapor_transfer_model = bulk_aero`` and
``coefficient_source = fixed``,

.. math::

   J_{\rho q_v,\mathrm{in}}
   = \rho\,C_E\,U_t\,
     [q_\mathrm{sat}(T_w,p_\mathrm{hse})-q_{v,a}].

``C_E`` is dimensionless.  Positive vapor flux adds vapor to the chamber;
negative vapor flux removes vapor at the wall.

For ``moisture = dry``, vapor impermeability is exact,

.. math::

   J_{\rho q_v,\mathrm{in}} = 0,

and this is not a ``qv_wall = 0`` Dirichlet condition.  For every supported
Stage 1 wall, cloud-water transfer is also exactly zero,

.. math::

   J_{\rho q_c,\mathrm{in}} = 0.

The dry-qv and qc zero-flux gates are applied before saturation or bulk
coefficient arithmetic.

For fixed ``bulk_aero``, :math:`U_t=0` gives exactly zero bulk heat and vapor
flux.  This implementation does not add a hidden minimum speed, gustiness
velocity, or free-convection term.  The model is therefore a
**forced/shear-dependent transfer closure**, not a complete natural- or
mixed-convection wall correlation.  This feature does not calibrate ``C_H``
or ``C_E`` for a Pi-Chamber or any other laboratory facility.  Numeric values
in regression inputs are software test parameters, not physical
recommendations.

ERF stores face flux in the positive coordinate direction.  The single
low/high storage adapter is

.. math::

   F_\mathrm{coord} =
   \begin{cases}
      +J_\mathrm{in}, & \text{low face},\\
      -J_\mathrm{in}, & \text{high face}.
   \end{cases}

Thus a high-face chamber influx is stored as a negative coordinate flux.
Cloud Chamber budgets use the matching low-minus-high convention.  The
retained physical face flux is the same value used by the RHS correction and
the budget; the budget does not recompute a separate wall model.

Per-face input contract
~~~~~~~~~~~~~~~~~~~~~~~

The face prefix is one of ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or
``zhi``.  The following table describes the physical-temperature Cloud
Chamber path.

.. list-table:: Cloud Chamber physical-wall input contract
   :header-rows: 1

   * - Key
     - Type
     - Units
     - Default if omitted
     - Valid values/domain
     - Required when
     - Invalid/unsupported combinations
   * - ``<face>.type``
     - string
     - --
     - none
     - ``NoSlipWall`` only
     - every face
     - any other wall type
   * - ``<face>.temperature``
     - real
     - K
     - none
     - finite and greater than zero
     - every physical face
     - legacy ``theta``/``qv`` wall keys in physical mode
   * - ``<face>.moisture``
     - string
     - --
     - none
     - ``dry`` or ``wet``
     - every physical face
     - ``wet`` without ``erf.moisture_model = SatAdj``
   * - ``<face>.wall_transfer_model``
     - string
     - --
     - omitted; resolved behavior remains the effective default
     - ``resolved_molecular`` only
     - legacy aggregate syntax only
     - cannot coexist with any per-channel model/coefficient key on the same face
   * - ``<face>.momentum_transfer_model``
     - string
     - --
     - ``resolved_noslip``
     - ``resolved_noslip`` only
     - optional explicit declaration
     - bulk momentum, law-wall, or other values are unsupported
   * - ``<face>.heat_transfer_model``
     - string
     - --
     - ``resolved_molecular``
     - ``resolved_molecular`` or ``bulk_aero``
     - optional
     - ``bulk_aero`` requires fixed coefficient source and ``C_H``
   * - ``<face>.vapor_transfer_model``
     - string
     - --
     - ``resolved_molecular``
     - ``resolved_molecular`` or ``bulk_aero``
     - optional
     - ``bulk_aero`` requires fixed coefficient source and ``C_E``; a dry wall still has exact zero vapor flux
   * - ``<face>.coefficient_source``
     - string
     - --
     - omitted
     - ``fixed`` only
     - at least one scalar channel on that face is ``bulk_aero``
     - rejected when no bulk scalar channel is active; MOST is unsupported
   * - ``<face>.C_H``
     - real
     - dimensionless
     - no physical default
     - finite and nonnegative
     - heat model is ``bulk_aero``
     - rejected if heat is not bulk
   * - ``<face>.C_E``
     - real
     - dimensionless
     - no physical default
     - finite and nonnegative
     - vapor model is ``bulk_aero``
     - rejected if vapor is not bulk
   * - ``<face>.velocity``
     - array
     - m s\ :sup:`-1`
     - zero internally
     - unsupported input in Stage 1
     - never
     - moving-wall metadata is rejected
   * - roughness / ``z0`` / ``z0_m`` / ``z0_h`` / ``z0_q``
     - real
     - m
     - none
     - unsupported in Cloud Chamber Stage 1
     - never
     - rejected
   * - ``<face>.C_D``
     - real
     - dimensionless
     - none
     - unsupported
     - never
     - rejected
   * - MOST-related Cloud Chamber face inputs
     - various
     - --
     - none
     - unsupported
     - never
     - rejected

Do not mix the legacy aggregate key with per-channel keys on the same face.
For example, this is invalid:

.. code-block:: none

   zlo.wall_transfer_model = resolved_molecular
   zlo.heat_transfer_model = bulk_aero

Choose either the retained aggregate resolved syntax or the per-channel
syntax.

A future Cloud Chamber-specific MOST coefficient provider is represented only
by the internal ``CoefficientProvider::MOSTFuture`` metadata placeholder.  It
is not parser-selectable, has no active Cloud Chamber runtime dispatch, and
has not been validated by this feature.  This placeholder does not describe
or alter ERF's existing atmospheric SurfaceLayer MOST implementation.
Cloud Chamber MOST, roughness, bulk-momentum, and law-of-the-wall inputs remain
unsupported and are rejected.

.. warning::

   Treat ``alpha_T`` and ``alpha_C`` as prescribed resolved-transfer
   coefficients, not calibrated molecular-property values or engineering
   wall-law coefficients.  ``C_H`` and ``C_E`` are user-supplied,
   dimensionless bulk coefficients.  They must be finite and nonnegative,
   but this feature does not provide calibrated or recommended physical
   values.

Bulk wall-rate timestep guard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When any active scalar wall channel uses ``bulk_aero``, ERF estimates the
maximum physical-boundary relaxation rate

.. math::

   \lambda_\mathrm{max}
   = \max_f\left(C_f U_{t,f}\,\Delta n_f^{-1}\right)

for active bulk heat and wet bulk-vapor channels.  Dry vapor is excluded from
the wall-rate scan because its wall flux is exactly zero.  The explicit wall
guard uses

.. math::

   dt_\mathrm{wall} = \frac{0.5}{\lambda_\mathrm{max}}.

The value 0.5 is an engineering safety factor for the explicit wall
relaxation, not a physical constant and not a proof of stability of the fully
coupled ERF operator.  Adaptive stepping takes the minimum of the ordinary ERF
limits and ``dt_wall``.  A positive fixed ``erf.fixed_dt`` larger than the
wall limit aborts with the measured ``fixed_dt``, ``wall_dt``, and
``max_wall_rate``.  ERF does not silently clip ``C_H``, ``C_E``, or the
requested fixed timestep.

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
6. If using ``bulk_aero``, set ``coefficient_source = fixed`` and the required
   ``C_H``/``C_E`` values; do not combine these keys with the aggregate wall
   key.
7. Run a short case and inspect temperature, potential temperature, velocity,
   and, for SatAdj, ``qv``, ``qc``, saturation mixing ratio, and RH.
8. Enable ``erf.cloud_chamber_budget_interval``.
9. Require dry thermal closure or total-water closure, as appropriate.
10. Treat any ``FAIL`` or budget-dependent solution change as invalid.

Stage 1 invariants
------------------

* One uniform Cartesian mesh level and no periodic direction; the level may
  be decomposed into multiple AMReX boxes.
* Six stationary ``NoSlipWall`` faces.
* Dry-wall vapor flux is exactly zero.
* Cloud-water wall flux is exactly zero.
* Wet walls are permitted only with SatAdj.
* Bulk heat and vapor models are independently selectable; dry vapor and
  cloud-water wall fluxes remain exactly zero.
* Momentum remains resolved no-slip in production; the momentum metadata has
  explicit ownership/orientation hooks for future wall models.
* Bulk walls enforce the ``Lambda <= 0.5`` timestep condition.
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
* ``bulk_aero`` requires per-channel fixed coefficients and cannot be mixed
  with ``wall_transfer_model`` on the same face.
* Treat ``UNSUPPORTED_SOURCE`` as incomplete cloudy thermal accounting, not
  as a successful budget result.
