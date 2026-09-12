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

This configuration is not a quantitatively validated Pi-Chamber LES.
It is intended for code-path verification, qualitative buoyancy-driven flow,
wall-flux testing, and conserved-scalar regression tests.  It does not provide
finite-rate droplet microphysics, a calibrated engineering wall law,
grid-independent LES validation, or quantitative experimental calibration.

Choose a mode
-------------

.. list-table:: Supported Cloud Chamber modes
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

   The Cloud Chamber uses one uniform Cartesian mesh level with no AMR refinement.
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
     path supported by the Cloud Chamber wall treatment.
   * ``erf.terrain_type`` and ``erf.buildings_type`` must be omitted or set
     to ``None``; terrain, embedded boundaries, and immersed buildings are
     outside the Cloud Chamber configuration.
   * All six faces must be stationary ``NoSlipWall`` boundaries with an
     explicit temperature.
   * The retained aggregate wall key
     ``wall_transfer_model = resolved_molecular`` remains supported.  The
     generalized contract may instead select supported momentum, heat, and
     vapor models independently.  Bulk channels use fixed, nonnegative
     coefficients where selected, or MOST-provided coefficients on the
     supported horizontal ``zlo``/``zhi`` walls.

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

   <face>.momentum_transfer_model = resolved_noslip|bulk_aero|neutral_roughness_log|law_of_wall_momentum
   <face>.heat_transfer_model = resolved_molecular|bulk_aero|neutral_roughness_log
   <face>.vapor_transfer_model = resolved_molecular|bulk_aero|neutral_roughness_log

If at least one channel uses ``bulk_aero``, also set:

.. code-block:: none

   <face>.coefficient_source = fixed|most

With ``coefficient_source = fixed``, provide ``C_D`` for bulk momentum,
``C_H`` for bulk heat, and ``C_E`` for bulk wet vapor.  With
``coefficient_source = most``, these coefficients are computed from the
shared pointwise MOST state and must not be supplied.  MOST is valid only on
``zlo`` and ``zhi``.  For example, a face using bulk heat and resolved wet
vapor requires:

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

For a stationary Cloud Chamber wall, the bulk model uses the local tangential
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
Cloud Chamber wall, cloud-water transfer is also exactly zero,

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

Neutral roughness-log wall transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``neutral_roughness_log`` is the canonical aerodynamic roughness/log-resistance
model.  ``law_of_wall_momentum`` is accepted as a compatibility alias for the
momentum selector; it does not create a second production law.  These models
are intended for a forced, shear-neutral wall-adjacent sample; no calibrated
roughness defaults are provided.

For outward unit normal :math:`\mathbf n`, the closure uses

.. math::

   \mathbf u_\mathrm{rel}=\mathbf u_\mathrm{fluid}-\mathbf u_\mathrm{wall},
   \qquad \mathbf u_t=(\mathbf I-\mathbf n\mathbf n^T)\mathbf u_\mathrm{rel},
   \qquad U_t=|\mathbf u_t|.

The normal velocity is removed before the magnitude is formed. With
:math:`z_\mathrm{ref}=\Delta n/2=0.5/\mathrm{dx\_inv}[DIR]`,

.. math::

   D_m=\ln(z_\mathrm{ref}/z_{0m}),\qquad
   u_* = \frac{\kappa U_t}{D_m},\qquad
   C_D=\left(\frac{\kappa}{D_m}\right)^2,
   \qquad \kappa=KAPPA=0.41.

The physical wall-on-fluid tangential traction is

.. math::

   \mathbf t_\mathrm{wall\to\mathrm{fluid}}
   =-\rho C_D U_t\mathbf u_t.

It is tangent, opposes wall-relative motion, and is exactly zero for a calm
sample. This physical traction is distinct from the low/high sign used by
ERF's stored ``tau12``, ``tau13``, and ``tau23`` arrays: a dedicated stress
adapter maps it to the existing momentum divergence operator. Neutral momentum
owns tangential cross-stress only and does not overwrite normal diagonal
stress.

Neutral heat and wet-vapor transfer use

.. math::

   C_H=\frac{\kappa^2}{D_m\ln(z_\mathrm{ref}/z_{0h})},\qquad
   C_E=\frac{\kappa^2}{D_m\ln(z_\mathrm{ref}/z_{0q})},

   J_{\rho\theta,\mathrm{in}}=\rho C_H U_t(\theta_w-\theta_a),\qquad
   J_{\rho q_v,\mathrm{in}}=\rho C_E U_t(q_{v,w}-q_{v,a}).

The existing Exner convention is used,
:math:`\Pi=(p_\mathrm{hse}/p_0)^{R_d/c_p}` and
:math:`\theta_w=T_w/\Pi`; wet vapor uses ERF's existing ``erf_qsatw``
routine and pressure conversion. Dry vapor and cloud-water wall fluxes remain
owned exact zeros before unrelated saturation or velocity arithmetic.

Every active neutral roughness length is host-validated. A scalar-only
neutral heat or vapor channel still requires ``z0_m`` because its resistance
contains ``ln(z_ref/z0_m)``:

.. math::

   0 < z_{0m},z_{0h},z_{0q} < z_\mathrm{ref},
   \qquad z_\mathrm{ref}=0.5/\mathrm{dx\_inv}[DIR].

Invalid values identify the face, key, supplied value, and reference distance;
there is no clipping or epsilon floor. The exact keys are:

.. list-table:: Neutral roughness keys
   :header-rows: 1
   :widths: 22 18 40 20

   * - Key
     - Units
     - Required when
     - Meaning
   * - ``<face>.z0_m``
     - m
     - any neutral momentum, heat, or vapor channel
     - aerodynamic momentum roughness
   * - ``<face>.z0_h``
     - m
     - neutral heat
     - scalar heat roughness
   * - ``<face>.z0_q``
     - m
     - neutral vapor, including a dry neutral-vapor gate
     - scalar moisture roughness

Syntax-only examples (not calibration recommendations) are:

.. code-block:: none

   xlo.momentum_transfer_model = neutral_roughness_log
   xlo.heat_transfer_model = neutral_roughness_log
   xlo.vapor_transfer_model = resolved_molecular
   xlo.z0_m = 1.0e-3
   xlo.z0_h = 5.0e-4

   zlo.momentum_transfer_model = neutral_roughness_log
   zlo.heat_transfer_model = neutral_roughness_log
   zlo.vapor_transfer_model = neutral_roughness_log
   zlo.z0_m = 1.0e-3
   zlo.z0_h = 5.0e-4
   zlo.z0_q = 5.0e-4

Momentum, heat, and vapor selectors are independent. A face may mix neutral
momentum with fixed bulk heat or resolved vapor; fixed bulk channels still
require their own coefficient. Neutral channels never consume user-supplied
bulk coefficients.

MOST wall-transfer provider
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``coefficient_source = most`` selects one shared, device-safe pointwise
Monin--Obukhov similarity state for all active bulk momentum, heat, and vapor
channels on a face.  The implementation uses ERF's
``similarity_funs::calc_psi_m2`` and ``calc_psi_h2`` functions, initializes the
stability coordinate neutrally, clamps bulk Richardson number to ``[-4,4]``,
uses relaxation ``0.5``, tolerance ``1e-3``, and a bounded iteration count.
Failure to converge or a non-finite result is an error; there is no silent
fallback.

The wall state uses

.. math::

   \theta_v=\theta(1+\epsilon_v q_v),\qquad
   q_{v,w}=\begin{cases}q_\mathrm{sat}(T_w,p_\mathrm{hse}),&\text{wet},\\
   q_{v,a},&\text{dry},\end{cases}

and ``gravity_sign = +1`` at ``zlo`` and ``-1`` at ``zhi``.  ``U_MOST`` is
``max(U_t,0.01)`` only for conditioning the stability iteration; physical
momentum and scalar fluxes still multiply the actual ``U_t`` and are exactly
zero at rest.  Momentum, heat, and vapor use separate ``z0_m``, ``z0_h``, and
``z0_q`` roughness lengths.  MOST is restricted to horizontal faces; neutral
roughness and fixed bulk channels may still be used independently on side
walls.  The scope remains single-level Cartesian ConstantDz, stationary-wall,
nonperiodic Cloud Chamber physics; smooth-wall, natural/mixed-convection,
terrain, embedded-boundary, AMR, stretched-``dz``, and moving-wall support are
out of scope.

The wall potential temperature uses the configured solver exponent
``rdOcp = R_d/c_p``; the same sampled value is used when the resulting
coefficient drives the scalar flux.

Developer contract for Cloud Chamber wall closures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cloud Chamber wall closures are organized as a sequence of small,
device-safe contracts:

* **Configuration layer:** host parsing validates the public input and produces
  compact per-face metadata.
* **Sampling layer:** the application adapter converts ERF staggered/grid
  variables into the physical values required by a closure.
* **Closure layer:** pure closure functions consume physical values and return
  inward-positive ``WallFlux`` results with authoritative ownership.
* **Ownership:** ``OwnNone`` means ERF retains its normal numerical treatment.
  An owned zero means Cloud Chamber replaces the normal wall flux with exactly
  zero.
* **Scalar application:** the adapter converts inward-positive flux to ERF's
  low/high coordinate sign, stores the retained physical flux, and applies
  only the ``new - old`` correction to the already-computed scalar RHS.
* **Budgets:** Cloud Chamber budgets consume the same retained diffusion face
  flux; they do not recompute a second wall closure.
* **Momentum:** a pure closure returns physical wall-on-fluid tangential
  traction; a dedicated stress adapter maps it to retained ERF cross-stress,
  and ``DiffusionSrcForMom`` applies the existing divergence. Do not inject
  traction directly into ``ERF_SlowRhsPre.cpp`` without a new discrete
  derivation.

To add another wall model:

#. Extend validated host configuration.
#. Preserve GPU-safe face metadata.
#. Identify the required sampled physical state.
#. Implement a pure closure evaluator.
#. Return inward-positive flux and authoritative ownership.
#. Reuse the scalar application adapter for scalar fluxes.
#. Do not put closure equations in ``ERF_SlowRhsPre.cpp`` or
   ``ERF_SlowRhsPost.cpp``.
#. Derive and validate a separate momentum stress adapter for traction.
#. Add independent equation tests.
#. Add orientation/application tests.
#. Add production-path activation and conservation coverage.
#. Mutation-test the critical seam.
#. Document equations, units, validity domain, and limitations.

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
     - ``resolved_noslip``, ``bulk_aero``, ``neutral_roughness_log``, or ``law_of_wall_momentum``
     - optional explicit declaration
     - bulk momentum requires ``coefficient_source`` and ``C_D`` only for fixed coefficients; law-wall is a neutral-log alias
   * - ``<face>.heat_transfer_model``
     - string
     - --
     - ``resolved_molecular``
     - ``resolved_molecular``, ``bulk_aero``, or ``neutral_roughness_log``
     - optional
     - ``bulk_aero`` requires fixed coefficient source and ``C_H``; neutral requires ``z0_m`` and ``z0_h``
   * - ``<face>.vapor_transfer_model``
     - string
     - --
     - ``resolved_molecular``
     - ``resolved_molecular``, ``bulk_aero``, or ``neutral_roughness_log``
     - optional
     - ``bulk_aero`` requires fixed coefficient source and ``C_E``; neutral requires ``z0_m`` and ``z0_q``; dry vapor is still exact zero
   * - ``<face>.coefficient_source``
     - string
     - --
     - omitted
     - ``fixed`` or ``most``
     - at least one bulk channel on that face
     - ``most`` is restricted to ``zlo``/``zhi``; rejected when no bulk channel is active
   * - ``<face>.C_D``
     - real
     - dimensionless
     - no physical default
     - finite and nonnegative
     - momentum model is ``bulk_aero`` with ``coefficient_source = fixed``
     - rejected for non-bulk momentum or with ``coefficient_source = most``
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
     - unsupported input in the Cloud Chamber configuration
     - never
     - moving-wall metadata is rejected
   * - ``<face>.z0_m``
     - real
     - m
     - none
     - finite and ``0 < z0_m < z_ref``
     - any neutral channel or MOST bulk channel
     - rejected when no active channel consumes it
   * - ``<face>.z0_h``
     - real
     - m
     - none
     - finite and ``0 < z0_h < z_ref``
     - neutral heat or MOST bulk heat
     - rejected when heat is neither neutral nor MOST bulk
   * - ``<face>.z0_q``
     - real
     - m
     - none
     - finite and ``0 < z0_q < z_ref``
     - neutral vapor or MOST bulk vapor
     - rejected when vapor is neither neutral nor MOST bulk

Do not mix the legacy aggregate key with per-channel keys on the same face.
For example, this is invalid:

.. code-block:: none

   zlo.wall_transfer_model = resolved_molecular
   zlo.heat_transfer_model = bulk_aero

Choose either the retained aggregate resolved syntax or the per-channel
syntax.

.. warning::

   Treat ``alpha_T`` and ``alpha_C`` as prescribed resolved-transfer
   coefficients, not calibrated molecular-property values or engineering
   wall-law coefficients.  ``C_H`` and ``C_E`` are user-supplied,
   dimensionless bulk coefficients.  They must be finite and nonnegative,
   but this feature does not provide calibrated or recommended physical
   values.

Neutral and bulk wall-rate timestep guard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERF scans the wall-adjacent tangential speed using the same half-cell sample
as the closure. Active bulk scalar channels contribute

.. math::

   \lambda_s = C_s U_{t,s}\,\Delta n^{-1},

where ``s`` is heat or wet vapor. Neutral scalar channels use their derived
``C_H`` or ``C_E`` in the same expression. Dry vapor is excluded because its
owned wall flux is exactly zero, even if neutral-vapor metadata is present.

For neutral momentum, the local tangential drag equation is

.. math::

   \frac{d\mathbf u_t}{dt}
   =-C_D\,\Delta n^{-1}|\mathbf u_t|\mathbf u_t.

The nonlinear tangential closure has Jacobian

.. math::

   J = C_D \left(U_t I + \frac{\mathbf u_t \mathbf u_t^T}{U_t}\right),

whose two eigenvalues are ``C_D U_t Delta n^-1`` and ``2 C_D U_t
Delta n^-1`` after the wall-normal divergence scaling.  The timestep guard
uses the induced infinity norm (the maximum absolute row sum), not the
spectral norm.  For a two-component tangential vector its exact local bound
is

.. math::

   K_\infty C_D U_t = \frac{3+\sqrt{2}}{2} C_D U_t.

Each staggered tangential reconstruction is a nonnegative average whose
weights sum to one.  Each boundary stress-node reconstruction is also a
nonnegative average with weights summing to one.  These maps are therefore
non-expansive in the infinity norm.  In the Cartesian Cloud Chamber
operator, the density multiplying the physical traction is the mass
weight of the conserved momentum source, so the mass-weighted momentum
Jacobian has no additional density ratio.  At an edge or corner, one free
staggered component can receive contributions from every active perpendicular
wall, so the conservative composed-operator bound is

.. math::

   \lambda_{m,c} = \sum_{f:\,n_f \ne c}
      K_\infty C_{D,f} U_{t,f} \Delta n_f^{-1},
   \qquad \lambda_m = \max_c \lambda_{m,c}.

The implementation computes this bound with
``momentum_infinity_row_sum_factor()`` and a finite-difference
check exercises the actual wall-stress adapter followed by
``DiffusionSrcForMom`` at an oblique state.  Low/high storage signs do not
change the absolute row sum; they only determine the physical drag sign.
The scan takes the maximum over cells, components, active scalar channels, and
neutral-momentum row sums, then applies the existing explicit safety factor,

.. math::

   \lambda_\mathrm{max}=\max_f(\lambda_{m,f},\lambda_{s,f}),
   \qquad dt_\mathrm{wall}=\frac{0.5}{\lambda_\mathrm{max}}.

Adaptive stepping takes the minimum of ordinary ERF limits and
``dt_wall``. A positive fixed ``erf.fixed_dt`` larger than the wall limit
aborts with the measured ``fixed_dt``, ``wall_dt``, and ``max_wall_rate``. ERF
does not silently clip roughness, coefficients, or a requested fixed step.

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
5. Set ``alpha_T`` and, when resolved vapor transfer is needed, ``alpha_C``.
6. If using ``bulk_aero``, set ``coefficient_source = fixed`` with the
   required finite nonnegative ``C_D``/``C_H``/``C_E`` values, or
   ``coefficient_source = most`` on ``zlo``/``zhi`` with the required positive
   ``z0_*`` values.  If using ``neutral_roughness_log``, provide its required
   positive ``z0_*`` values.  Do not combine per-channel keys with the
   aggregate wall key.
7. Run a short case and inspect temperature, potential temperature, velocity,
   and, for SatAdj, ``qv``, ``qc``, saturation mixing ratio, and RH.
8. Enable ``erf.cloud_chamber_budget_interval``.
9. Require dry thermal closure or total-water closure, as appropriate.
10. Treat any ``FAIL`` or budget-dependent solution change as invalid.

Cloud Chamber invariants
------------------

* One uniform Cartesian mesh level and no periodic direction; the level may
  be decomposed into multiple AMReX boxes.
* Six stationary ``NoSlipWall`` faces.
* Dry-wall vapor flux is exactly zero.
* Cloud-water wall flux is exactly zero.
* Wet walls are permitted only with SatAdj.
* Bulk heat and vapor models are independently selectable; dry vapor and
  cloud-water wall fluxes remain exactly zero.
* Momentum models may be resolved no-slip, fixed bulk aerodynamic, or
  neutral-log; active rough/bulk momentum owns only tangential cross-stress at
  physical chamber faces.
* Neutral and bulk walls enforce ``fixed_dt <= dt_wall = 0.5 / max_wall_rate``.
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
* ``bulk_aero`` requires a per-channel coefficient source and fixed
  coefficients unless MOST is selected; it cannot be mixed with
  ``wall_transfer_model`` on the same face.
* ``coefficient_source = most`` is restricted to ``zlo``/``zhi`` and shares one
  converged pointwise stability state across active bulk channels.
* ``neutral_roughness_log`` requires explicit roughness lengths and the
  grid-relative inequality ``0 < z0_* < z_ref``; there are no calibrated
  defaults.
* Treat ``UNSUPPORTED_SOURCE`` as incomplete cloudy thermal accounting, not
  as a successful budget result.
