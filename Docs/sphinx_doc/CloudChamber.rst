.. index:: Cloud Chamber
.. index:: SatAdj
.. index:: ConstantAlpha
.. index:: NoSlipWall
.. index:: bulk_aero
.. index:: neutral_roughness_log
.. index:: Monin-Obukhov similarity theory
.. index:: cloud_chamber_budget_interval
.. index:: qv
.. index:: qc
.. index:: rhoTheta

.. _CloudChamber:

Cloud Chamber
=============

Overview
--------

The Cloud Chamber problem is a single-level, Cartesian, fixed-density
anelastic configuration for studying buoyancy-driven thermal and moist flow in
a closed rectangular chamber.  It supports prescribed-temperature walls,
dry or SatAdj thermodynamics, and several wall-transfer closures for momentum,
heat, and water vapor.

The Cloud Chamber implementation is intended for code-path verification,
qualitative chamber-flow studies, wall-flux experiments, and conservation
regression tests.  It is **not** a quantitatively validated Pi-Chamber LES or
a calibrated engineering wall-model package.  It does not provide finite-rate
droplet microphysics, droplet deposition, precipitation, natural-convection
wall correlations, or experimentally calibrated transfer coefficients.

Most new inputs should use
``prob.thermodynamic_initialization = physical_temperature_rh``.  The older
``legacy_theta_qv`` mode is retained for existing numerical inputs but does
not support the generalized per-channel wall-transfer models described below.

Choose the thermodynamic mode
-----------------------------

.. list-table:: Supported Cloud Chamber thermodynamic modes
   :header-rows: 1
   :widths: 22 24 24 30

   * - Mode
     - Moisture model
     - Wall moisture
     - Typical use
   * - Dry physical chamber
     - Omit ``erf.moisture_model``
     - Every face is ``dry``
     - Thermal convection and strict dry ``rhoTheta`` budget checks
   * - SatAdj with dry walls
     - ``erf.moisture_model = SatAdj``
     - Every face is ``dry``
     - Moist thermodynamics in a water-impermeable chamber
   * - SatAdj with wet walls
     - ``erf.moisture_model = SatAdj``
     - One or more faces may be ``wet``
     - Equilibrium vapor exchange with total-water budget checks
   * - Legacy theta/qv
     - Omit moisture for dry legacy mode, or use ``SatAdj`` for cloudy legacy mode
     - Legacy numerical ``theta``/``qv`` wall values
     - Backward compatibility with existing Cloud Chamber inputs

A wall marked ``dry`` is impermeable to water but may still exchange heat.
A ``wet`` wall is valid only when ``erf.moisture_model = SatAdj``.

Terminology
-----------

``rhoTheta``
   Density-weighted potential temperature, the conserved thermal scalar used
   by this Cloud Chamber configuration.

``qv`` and ``qc``
   Water-vapor and cloud-water mixing ratios relative to dry air.  In SatAdj
   cases ERF partitions total nonprecipitating water between these two phases.

SatAdj
   ERF's instantaneous saturation-adjustment thermodynamics, including the
   existing latent-heating adjustment associated with vapor--cloud-water
   equilibrium.

Resolved wall transfer
   Transfer obtained from ERF's resolved no-slip momentum treatment or from
   the adjacent-cell-to-wall half-cell scalar gradient.  It is not an
   engineering wall-function correlation.

Bulk coefficient
   A dimensionless transfer coefficient multiplying the actual tangential
   speed in the Cloud Chamber forced/shear-dependent wall-transfer formulas.

Global ERF requirements
-----------------------

The Cloud Chamber parser enforces a deliberately narrow numerical scope.
Supported Cloud Chamber cases use one nonperiodic Cartesian level,
fixed-density anelastic dynamics, gravity, ConstantDz geometry, explicit
vertical diffusion, and ConstantAlpha scalar diffusion.

A typical common block is:

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

The following constraints are part of the implemented Cloud Chamber contract:

.. list-table:: Global Cloud Chamber requirements
   :header-rows: 1
   :widths: 28 26 46

   * - Setting
     - Required value
     - Notes
   * - ``erf.anelastic``
     - ``1``
     - Required in every Cloud Chamber mode.
   * - ``erf.use_gravity``
     - ``true``
     - Required in every Cloud Chamber mode.
   * - ``amr.max_level``
     - ``0``
     - AMR refinement is unsupported.  The single level may still contain multiple AMReX boxes.
   * - ``geometry.is_periodic``
     - ``0 0 0``
     - All three directions are closed and nonperiodic.
   * - ``erf.mesh_type``
     - omitted or ``ConstantDz``
     - Stretched vertical grids are outside the supported Cloud Chamber path.
   * - ``erf.terrain_type``
     - omitted or ``None``
     - Terrain and embedded-boundary configurations are unsupported.
   * - ``erf.buildings_type``
     - omitted or ``None``
     - Immersed buildings are unsupported.
   * - ``erf.init_type``
     - ``ConstantDensity``
     - The supported Cloud Chamber initializer requires a fixed-density
       anelastic base state in both physical and legacy thermodynamic modes.
   * - ``erf.molec_diff_type``
     - ``ConstantAlpha``
     - Required by the physical wall-transfer path, including cases that also use bulk, neutral, or MOST wall closures.
   * - ``erf.vert_implicit``
     - ``false``
     - Required by the physical wall-transfer path.

If ``prob.p_inf`` is supplied, the Cloud Chamber requires it to agree with
ERF's reference pressure ``p_0 = 100000 Pa`` within a relative tolerance of
``1e-6``.  ``prob.p_inf`` does not redefine the reference pressure used by the
Cloud Chamber potential-temperature conversion.

Physical-temperature initialization
-----------------------------------

Set:

.. code-block:: none

   prob.thermodynamic_initialization = physical_temperature_rh

The physical-temperature initializer constructs a linear vertical temperature
profile and optionally adds a deterministic three-dimensional perturbation.

.. list-table:: Physical initialization inputs
   :header-rows: 1
   :widths: 34 18 20 28

   * - Key
     - Required
     - Default
     - Valid values / meaning
   * - ``prob.initial_temperature_bottom``
     - yes
     - none
     - finite and greater than zero, in K
   * - ``prob.initial_temperature_top``
     - yes
     - none
     - finite and greater than zero, in K
   * - ``prob.temperature_perturbation_amplitude``
     - no
     - ``0``
     - finite amplitude in K
   * - ``prob.perturbation_mode``
     - no
     - ``deterministic_sine``
     - the only currently supported perturbation mode
   * - ``prob.initial_relative_humidity``
     - SatAdj only
     - none
     - required with SatAdj; a fraction in ``[0,1]``

For a domain with lower corner :math:`(x_0,y_0,z_0)` and lengths
:math:`(L_x,L_y,L_z)`, the temperature field is

.. math::

   T(x,y,z) = T_\mathrm{linear}(z)
   + A_T
     \sin\!\left(\frac{2\pi(x-x_0)}{L_x}\right)
     \sin\!\left(\frac{2\pi(y-y_0)}{L_y}\right)
     \sin\!\left(\frac{\pi(z-z_0)}{L_z}\right),

where :math:`A_T` is ``prob.temperature_perturbation_amplitude`` and
:math:`T_\mathrm{linear}` interpolates between the prescribed bottom and top
temperatures.

Potential temperature is derived from the local hydrostatic base-state
pressure,

.. math::

   \theta = T
   \left(\frac{p_0}{p_\mathrm{hse}}\right)^{R_d/c_p},

using ERF's configured ``rdOcp = R_d/c_p``.  Here
:math:`p_\mathrm{hse}` is the local hydrostatic pressure and
:math:`p_0=100000\ \mathrm{Pa}` is ERF's reference pressure.

With SatAdj, the initial vapor mixing ratio is computed from the supplied
relative humidity,

.. math::

   q_v = \frac{R_d}{R_v}
   \frac{RH\,e_s(T)}{p_\mathrm{hse}-RH\,e_s(T)},

where ``RH`` is a fraction rather than a percentage.  The initializer sets
``qc = 0`` and SatAdj subsequently establishes the equilibrium vapor--cloud
water partition.  In dry mode there are no ``qv`` or ``qc`` state fields, and
``prob.initial_relative_humidity`` must be omitted.

Do not combine physical initialization with the legacy profile keys
``theta_bottom``, ``theta_top``, ``theta_perturbation_amplitude``,
``qv_bottom``, or ``qv_top``.

Legacy theta/qv initialization
------------------------------

Set:

.. code-block:: none

   prob.thermodynamic_initialization = legacy_theta_qv

Legacy mode is provided for existing numerical inputs.  It uses potential
temperature directly instead of physical temperature and does not support the
generalized per-channel wall-transfer keys.

Dry legacy mode requires:

.. code-block:: none

   prob.theta_bottom = <positive value>
   prob.theta_top = <positive value>
   prob.theta_perturbation_amplitude = <optional finite value>
   prob.perturbation_mode = deterministic_sine

For cloudy legacy mode also set ``erf.moisture_model = SatAdj`` and provide
nonnegative ``prob.qv_bottom`` and ``prob.qv_top``.  Every face still requires
``<face>.type = NoSlipWall``.  Legacy walls use ``<face>.theta``; cloudy
legacy walls also require ``<face>.qv``.

Do not combine ``legacy_theta_qv`` with physical profile keys such as
``initial_temperature_bottom``, ``initial_temperature_top``,
``initial_relative_humidity``, or ``temperature_perturbation_amplitude``.
For new generalized wall-model studies, use ``physical_temperature_rh``.

Optional deterministic initial velocity perturbation
----------------------------------------------------

``prob.U_0`` (m s\ :sup:`-1`) controls an optional deterministic horizontal
velocity perturbation.  It is available with either thermodynamic
initialization mode and is an initial-condition option only; it does not
represent wall motion.

ERF initializes the perturbation on the staggered velocity grid so that its
cell-centered finite-volume horizontal divergence is zero to roundoff and
sets

.. math::

   w = 0.

The pattern is smooth, deterministic, and vanishes consistently with the
stationary chamber-wall construction.  On a square horizontal domain with
equal horizontal resolution it reduces to

.. math::

   u = U_0
   \sin^2\!\left(\frac{\pi(x-x_0)}{L_x}\right)
   \sin\!\left(\frac{2\pi(y-y_0)}{L_y}\right)
   \sin^2\!\left(\frac{\pi(z-z_0)}{L_z}\right),

.. math::

   v = -U_0
   \sin\!\left(\frac{2\pi(x-x_0)}{L_x}\right)
   \sin^2\!\left(\frac{\pi(y-y_0)}{L_y}\right)
   \sin^2\!\left(\frac{\pi(z-z_0)}{L_z}\right).

For a general rectangular/discretized horizontal grid, ERF applies the
corresponding grid-aspect scaling required for the staggered finite-volume
divergence to cancel.  Leave ``prob.U_0`` unset, or set it to zero, for no
initial velocity perturbation.

Physical wall configuration
---------------------------

In physical-temperature mode, every face requires three basic declarations:

.. code-block:: none

   <face>.type = NoSlipWall
   <face>.temperature = <temperature in K>
   <face>.moisture = dry|wet

The face prefix is one of ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or
``zhi``.  Wall temperature must be finite and positive.  A wet wall requires
SatAdj.  Cloud Chamber walls are stationary; ``<face>.velocity`` is rejected.

``NoSlipWall`` is the required ERF boundary type for all supported Cloud
Chamber wall models.  Selecting an aerodynamic momentum model changes the
Cloud Chamber tangential wall-stress closure; it does not change the required
boundary type or create a moving wall.

Channel defaults
~~~~~~~~~~~~~~~~

If no per-channel model is specified, the effective physical-wall defaults
are:

.. list-table:: Default wall-transfer models
   :header-rows: 1
   :widths: 20 34 46

   * - Channel
     - Default model
     - Meaning
   * - momentum
     - ``resolved_noslip``
     - retain ERF's resolved no-slip viscous momentum treatment
   * - heat
     - ``resolved_molecular``
     - use the adjacent-cell-to-wall half-cell resolved scalar gradient
   * - vapor
     - ``resolved_molecular``
     - use the resolved wet-wall gradient; dry walls remain exactly impermeable
   * - cloud water
     - no selector
     - wall flux is always exactly zero

The retained aggregate compatibility key

.. code-block:: none

   <face>.wall_transfer_model = resolved_molecular

selects the fully resolved/default behavior for that face.  Do not combine
``wall_transfer_model`` with per-channel model, coefficient, or roughness
keys on the same face.

Choose a transfer model
~~~~~~~~~~~~~~~~~~~~~~~

The three configurable channels are selected independently:

.. code-block:: none

   <face>.momentum_transfer_model = resolved_noslip|bulk_aero|neutral_roughness_log|law_of_wall_momentum
   <face>.heat_transfer_model = resolved_molecular|bulk_aero|neutral_roughness_log
   <face>.vapor_transfer_model = resolved_molecular|bulk_aero|neutral_roughness_log

``law_of_wall_momentum`` is a compatibility alias for
``neutral_roughness_log``; it does not select a separate formulation.

.. list-table:: Wall-model choices
   :header-rows: 1
   :widths: 20 23 29 28

   * - Choice
     - Applicable channels
     - Additional inputs
     - Main behavior
   * - resolved
     - momentum, heat, vapor
     - none beyond the basic wall settings and resolved diffusivities
     - retains ERF resolved no-slip momentum or half-cell scalar transfer
   * - fixed ``bulk_aero``
     - momentum, heat, vapor
     - ``coefficient_source = fixed`` plus the corresponding ``C_D``, ``C_H``, or ``C_E``
     - user-supplied dimensionless transfer coefficient multiplied by actual tangential speed
   * - ``neutral_roughness_log``
     - momentum, heat, vapor
     - ``z0_m`` and the required scalar roughness lengths
     - neutral aerodynamic log-resistance closure
   * - MOST-provided ``bulk_aero``
     - momentum, heat, vapor
     - ``coefficient_source = most`` plus the required ``z0_*`` values
     - computes bulk coefficients from the current pointwise Monin--Obukhov stability state; ``zlo``/``zhi`` only

One coefficient provider per face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``<face>.coefficient_source`` is a **face-level provider shared by every
``bulk_aero`` channel on that face**.  It is required when at least one bulk
channel is active and must be either ``fixed`` or ``most``.  Consequently, a
single face cannot use fixed bulk momentum and MOST bulk heat at the same
time.  Non-bulk channels on that face may still be resolved or neutral.

With ``coefficient_source = fixed``, provide a coefficient for every selected
bulk channel:

* ``C_D`` for bulk momentum;
* ``C_H`` for bulk heat; and
* ``C_E`` for bulk vapor.

The coefficients must be finite and nonnegative.  They are dimensionless.
With ``coefficient_source = most``, do **not** provide ``C_D``, ``C_H``, or
``C_E``; ERF computes them from the current wall-adjacent state.

.. important::

   The parser validates selected vapor-model metadata even on a dry wall.
   Therefore a dry wall configured with fixed ``bulk_aero`` vapor still
   requires ``C_E``; a dry wall configured with neutral or MOST vapor still
   requires the corresponding roughness metadata.  The physical dry-wall
   vapor flux nevertheless remains exactly zero.

Roughness requirements
~~~~~~~~~~~~~~~~~~~~~~

For neutral and MOST closures, the wall-adjacent reference distance is

.. math::

   z_\mathrm{ref}=\frac{\Delta n}{2},

where :math:`\Delta n` is the cell width normal to the face.  Every required
roughness length must satisfy

.. math::

   0 < z_0 < z_\mathrm{ref}.

The parser does not clip invalid roughness values.

.. list-table:: Roughness inputs
   :header-rows: 1
   :widths: 18 22 60

   * - Key
     - Units
     - Required when
   * - ``<face>.z0_m``
     - m
     - any neutral channel, or any ``bulk_aero`` channel using MOST
   * - ``<face>.z0_h``
     - m
     - neutral heat, or MOST bulk heat
   * - ``<face>.z0_q``
     - m
     - neutral vapor, or MOST bulk vapor, including a configured dry vapor channel

When wet heat and wet vapor both use MOST on one face, ``z0_h`` and ``z0_q``
must be exactly equal because the implemented bulk-Richardson fixed point uses
one scalar stability resistance.  This equality is not imposed by a dry vapor
channel because dry vapor is not an active scalar-transfer channel.

Per-face input reference
~~~~~~~~~~~~~~~~~~~~~~~~

The following table summarizes the physical-temperature wall keys.

.. list-table:: Cloud Chamber physical-wall input reference
   :header-rows: 1
   :widths: 27 17 21 35

   * - Key
     - Default
     - Valid values
     - Required / rejected when
   * - ``<face>.type``
     - none
     - ``NoSlipWall``
     - required on every face
   * - ``<face>.temperature``
     - none
     - finite, positive K
     - required on every physical wall
   * - ``<face>.moisture``
     - none
     - ``dry`` or ``wet``
     - required on every physical wall; ``wet`` requires SatAdj
   * - ``<face>.wall_transfer_model``
     - omitted
     - ``resolved_molecular`` only
     - compatibility aggregate; rejected when combined with any per-channel transfer key
   * - ``<face>.momentum_transfer_model``
     - ``resolved_noslip``
     - ``resolved_noslip``, ``bulk_aero``, ``neutral_roughness_log``, ``law_of_wall_momentum``
     - optional; alias maps to neutral roughness-log
   * - ``<face>.heat_transfer_model``
     - ``resolved_molecular``
     - ``resolved_molecular``, ``bulk_aero``, ``neutral_roughness_log``
     - optional
   * - ``<face>.vapor_transfer_model``
     - ``resolved_molecular``
     - ``resolved_molecular``, ``bulk_aero``, ``neutral_roughness_log``
     - optional; dry physical vapor flux remains exactly zero
   * - ``<face>.coefficient_source``
     - omitted
     - ``fixed`` or ``most``
     - required if any channel is ``bulk_aero``; rejected if no bulk channel is active; ``most`` only on ``zlo``/``zhi``
   * - ``<face>.C_D``
     - none
     - finite and nonnegative
     - required only for fixed bulk momentum; rejected with MOST or non-bulk momentum
   * - ``<face>.C_H``
     - none
     - finite and nonnegative
     - required only for fixed bulk heat; rejected with MOST or non-bulk heat
   * - ``<face>.C_E``
     - none
     - finite and nonnegative
     - required only for fixed bulk vapor, including a configured dry vapor channel; rejected with MOST or non-bulk vapor
   * - ``<face>.z0_m``
     - none
     - ``0 < z0_m < z_ref``
     - required by any neutral channel or any MOST bulk channel; otherwise rejected
   * - ``<face>.z0_h``
     - none
     - ``0 < z0_h < z_ref``
     - required by neutral heat or MOST bulk heat; otherwise rejected
   * - ``<face>.z0_q``
     - none
     - ``0 < z0_q < z_ref``
     - required by neutral vapor or MOST bulk vapor; otherwise rejected
   * - ``<face>.velocity``
     - stationary internally
     - unsupported
     - always rejected; moving walls are not implemented

Configuration examples
~~~~~~~~~~~~~~~~~~~~~~

The snippets below illustrate syntax only.  Numeric coefficients and roughness
lengths are examples, not calibration recommendations.

Fully resolved physical wall:

.. code-block:: none

   zlo.type = NoSlipWall
   zlo.temperature = 300.0
   zlo.moisture = dry
   zlo.wall_transfer_model = resolved_molecular

Fixed bulk momentum and heat with resolved wet vapor:

.. code-block:: none

   zlo.type = NoSlipWall
   zlo.temperature = 300.0
   zlo.moisture = wet
   zlo.momentum_transfer_model = bulk_aero
   zlo.heat_transfer_model = bulk_aero
   zlo.vapor_transfer_model = resolved_molecular
   zlo.coefficient_source = fixed
   zlo.C_D = 0.01
   zlo.C_H = 0.01

Neutral roughness-log momentum, heat, and vapor:

.. code-block:: none

   xlo.type = NoSlipWall
   xlo.temperature = 292.0
   xlo.moisture = wet
   xlo.momentum_transfer_model = neutral_roughness_log
   xlo.heat_transfer_model = neutral_roughness_log
   xlo.vapor_transfer_model = neutral_roughness_log
   xlo.z0_m = 1.0e-3
   xlo.z0_h = 5.0e-4
   xlo.z0_q = 5.0e-4

MOST-provided bulk transfer on a horizontal wall:

.. code-block:: none

   zlo.type = NoSlipWall
   zlo.temperature = 300.0
   zlo.moisture = wet
   zlo.momentum_transfer_model = bulk_aero
   zlo.heat_transfer_model = bulk_aero
   zlo.vapor_transfer_model = bulk_aero
   zlo.coefficient_source = most
   zlo.z0_m = 2.0e-3
   zlo.z0_h = 2.0e-3
   zlo.z0_q = 2.0e-3

A face may mix non-bulk and bulk models.  For example, neutral momentum may be
combined with fixed bulk heat and resolved vapor:

.. code-block:: none

   ylo.momentum_transfer_model = neutral_roughness_log
   ylo.heat_transfer_model = bulk_aero
   ylo.vapor_transfer_model = resolved_molecular
   ylo.z0_m = 1.0e-3
   ylo.coefficient_source = fixed
   ylo.C_H = 0.01

Wall-transfer formulations
--------------------------

Common definitions and sign convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All Cloud Chamber wall closures use the adjacent wall cell and a reference
distance :math:`z_\mathrm{ref}=\Delta n/2`.  For outward unit normal
:math:`\mathbf n`, stationary wall velocity :math:`\mathbf u_w=0`, and
adjacent fluid velocity :math:`\mathbf u_c`, the tangential relative velocity
is

.. math::

   \mathbf u_t =
   (\mathbf I-\mathbf n\mathbf n^T)(\mathbf u_c-\mathbf u_w),
   \qquad U_t=|\mathbf u_t|.

The normal velocity is removed before :math:`U_t` is formed.

Physical scalar fluxes in this section are defined as positive **into the
chamber fluid** and are denoted :math:`J_\mathrm{in}`.  ERF stores a face flux
in the positive coordinate direction, so

.. math::

   F_\mathrm{coord}=
   \begin{cases}
      +J_\mathrm{in}, & \text{low face},\\
      -J_\mathrm{in}, & \text{high face}.
   \end{cases}

The wall potential temperature is

.. math::

   \Pi_w = \left(\frac{p_\mathrm{hse}}{p_0}\right)^{R_d/c_p},
   \qquad
   \theta_w = \frac{T_w}{\Pi_w}
   = T_w\left(\frac{p_0}{p_\mathrm{hse}}\right)^{R_d/c_p}.

Dry vapor and cloud water have exact zero physical wall flux for every
supported transfer model:

.. math::

   J_{\rho q_v,\mathrm{in}}=0 \quad \text{for a dry wall},
   \qquad
   J_{\rho q_c,\mathrm{in}}=0 \quad \text{for every wall}.

A dry wall is therefore not a ``qv_wall = 0`` Dirichlet condition.

Resolved transfer
~~~~~~~~~~~~~~~~~

``momentum_transfer_model = resolved_noslip`` leaves momentum to ERF's
existing resolved no-slip viscous treatment.

For ``heat_transfer_model = resolved_molecular``, the physical-temperature
wall replaces the scalar face flux with the half-cell resolved transfer

.. math::

   J_{\rho\theta,\mathrm{in}}
   = \rho\,\alpha_T\,(\theta_w-\theta_a)\frac{2}{\Delta n}.

For a wet wall with ``vapor_transfer_model = resolved_molecular``,

.. math::

   J_{\rho q_v,\mathrm{in}}
   = \rho\,\alpha_C
     \left[q_\mathrm{sat}(T_w,p_\mathrm{hse})-q_{v,a}\right]
     \frac{2}{\Delta n}.

The resolved scalar coefficients are configured through ``erf.alpha_T`` for
heat and ``erf.alpha_C`` for moisture.  Both have units of m\ :sup:`2`
s\ :sup:`-1`.  In this Cloud Chamber path they set the half-cell resolved
wall-normal scalar transfer; they are not calibrated aerodynamic exchange
coefficients.

``momentum_transfer_model = resolved_noslip`` retains ERF's resolved viscous
no-slip momentum-stress treatment and does not use a Cloud Chamber ``C_D`` or
roughness length.

Fixed bulk-aerodynamic transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``coefficient_source = fixed``, the user supplies the dimensionless bulk
coefficients.  Momentum traction on the fluid is

.. math::

   \mathbf t_\mathrm{wall\to fluid}
   = -\rho C_D U_t\mathbf u_t.

The traction is tangent to the wall, opposes tangential relative motion, and
is exactly zero when :math:`U_t=0`.

Bulk heat transfer is

.. math::

   J_{\rho\theta,\mathrm{in}}
   = \rho C_H U_t(\theta_w-\theta_a).

For a wet wall, bulk vapor transfer is

.. math::

   J_{\rho q_v,\mathrm{in}}
   = \rho C_E U_t
     \left[q_\mathrm{sat}(T_w,p_\mathrm{hse})-q_{v,a}\right].

The physical transfer uses the actual :math:`U_t`; there is no hidden minimum
speed, gustiness velocity, or free-convection augmentation.  Fixed bulk
transfer is therefore a forced/shear-dependent closure, not a complete
natural- or mixed-convection wall correlation.

Neutral roughness-log transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``neutral_roughness_log``, let

.. math::

   D_m = \ln\left(\frac{z_\mathrm{ref}}{z_{0m}}\right),
   \qquad
   \kappa=0.41.

The neutral momentum state is

.. math::

   u_* = \frac{\kappa U_t}{D_m},
   \qquad
   C_D=\left(\frac{\kappa}{D_m}\right)^2,

and momentum uses the same physical drag form

.. math::

   \mathbf t_\mathrm{wall\to fluid}
   = -\rho C_D U_t\mathbf u_t.

Neutral heat and wet-vapor coefficients are

.. math::

   C_H = \frac{\kappa^2}
              {D_m\ln(z_\mathrm{ref}/z_{0h})},
   \qquad
   C_E = \frac{\kappa^2}
              {D_m\ln(z_\mathrm{ref}/z_{0q})},

and the bulk scalar-flux formulas above are then used.  A scalar-only neutral
heat or vapor channel still requires ``z0_m`` because its resistance contains
:math:`D_m`.  Neutral channels do not use ``coefficient_source`` or
user-supplied ``C_D``, ``C_H``, or ``C_E``.

MOST-provided bulk coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``coefficient_source = most`` keeps the ``bulk_aero`` physical flux formulas
but obtains their coefficients from one pointwise Monin--Obukhov stability
state on the wall face.  It is supported only on ``zlo`` and ``zhi``.

The wall virtual potential-temperature state is

.. math::

   \theta_v=\theta(1+\epsilon_v q_v),

with

.. math::

   q_{v,w}=\begin{cases}
      q_\mathrm{sat}(T_w,p_\mathrm{hse}), & \text{wet wall},\\
      q_{v,a}, & \text{dry wall}.
   \end{cases}

The stability calculation uses

.. math::

   U_\mathrm{MOST}=\max(U_t,0.01),

and the bulk Richardson number

.. math::

   Ri_b = s_g\,
   \frac{g z_\mathrm{ref}}{\theta_a}
   \frac{\theta_{v,a}-\theta_{v,w}}{U_\mathrm{MOST}^2},

where :math:`s_g=+1` on ``zlo`` and :math:`s_g=-1` on ``zhi``.  ERF clamps

.. math::

   Ri=\mathrm{clamp}(Ri_b,-4,4).

Cloud-water loading is not included in the virtual-potential-temperature
bulk Richardson number used by this Cloud Chamber MOST closure.

Define

.. math::

   D_m=\ln(z_\mathrm{ref}/z_{0m}),\qquad
   D_h=\ln(z_\mathrm{ref}/z_{0h}),\qquad
   D_q=\ln(z_\mathrm{ref}/z_{0q}).

The scalar resistance used by the stability fixed point is ``z0_h`` when
MOST heat is active, otherwise ``z0_q`` when active wet MOST vapor is present,
and otherwise ``z0_m`` for momentum-only MOST.  If wet heat and wet vapor are
both active, the parser requires ``z0_h == z0_q``.

Starting from :math:`\zeta=0`, ERF evaluates its existing
``calc_psi_m2`` and ``calc_psi_h2`` similarity functions and iterates

.. math::

   A_m=\max(D_m-\psi_m,1),\qquad
   A_s=\max(D_s-\psi_h,1),

.. math::

   \zeta_\mathrm{new}
   =0.5\,\zeta_\mathrm{old}
   +0.5\,Ri\frac{A_m^2}{A_s}.

The iteration stops when
:math:`|\zeta_\mathrm{new}-\zeta_\mathrm{old}|\le 10^{-3}` and allows at
most 100 iterations.  Nonconvergence or a non-finite result is an error; there
is no silent neutral fallback.

After convergence,

.. math::

   A_m=\max(D_m-\psi_m,1),\qquad
   A_h=\max(D_h-\psi_h,1),\qquad
   A_q=\max(D_q-\psi_h,1),

and

.. math::

   C_D=\left(\frac{\kappa}{A_m}\right)^2,\qquad
   C_H=\frac{\kappa^2}{A_mA_h},\qquad
   C_E=\frac{\kappa^2}{A_mA_q}.

The internal MOST closure also forms the friction velocity

.. math::

   u_* = \frac{\kappa U_\mathrm{MOST}}{A_m}.

The ``0.01`` speed floor is used only to condition the stability calculation.
Physical momentum, heat, and vapor transfer still multiply the actual
:math:`U_t`, so all shear-dependent physical transfer is exactly zero at rest.

A configured dry MOST vapor channel is metadata only: its ``z0_q`` is still
required and validated, but dry vapor does not participate in the active
scalar stability resistance and its physical wall flux remains exactly zero.
When a moisture state is present, ambient ``qv`` still enters the virtual
potential temperature used by any active MOST momentum or heat calculation.

Momentum stress application
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bulk and neutral momentum closures return physical wall-on-fluid tangential
traction.  The Cloud Chamber stress adapter maps that traction onto ERF's
retained cross-stress components and the existing ``DiffusionSrcForMom``
operator applies the momentum divergence.  The wall closure does not replace
normal diagonal stress.  This distinction matters when comparing analytic
traction with stored stress-array signs.

Wall-transfer timestep guard
----------------------------

Bulk and neutral transfer can impose a stronger explicit timestep restriction
than the ordinary advective/diffusive CFL limits.  ERF therefore scans the
wall-adjacent state and constructs a wall rate.

For heat or active wet vapor,

.. math::

   \lambda_s = C_s U_t\,\Delta n^{-1}.

Neutral scalars use their derived coefficient.  Dry vapor does not contribute
because its wall flux is exactly zero.

For a bulk or neutral tangential momentum face, the local conservative
infinity-row-sum factor is

.. math::

   K_\infty=\frac{3+\sqrt{2}}{2},

so the per-face momentum contribution is

.. math::

   \lambda_{m,f}
   = K_\infty C_{D,f}U_{t,f}\Delta n_f^{-1}.

At an edge or corner, contributions from perpendicular active walls are added
for each free velocity component before taking the maximum.  MOST coefficients
are evaluated once for the current wall/cell sample and then treated as frozen
for this local rate estimate; this is not an exact Jacobian of the nonlinear
MOST iteration.

The wall limit is

.. math::

   dt_\mathrm{wall}=\frac{0.5}{\lambda_\mathrm{max}}.

Adaptive stepping takes the minimum of this limit and the ordinary ERF
constraints.  If a positive user-specified ``erf.fixed_dt`` exceeds the wall
limit, ERF aborts and reports the requested fixed step, wall limit, and maximum
wall rate; it does not silently clip the requested fixed step.

Conserved-scalar budget diagnostics
-----------------------------------

Set

.. code-block:: none

   erf.cloud_chamber_budget_interval = <positive number of steps>

to write ``cloud_chamber_budget.dat``.  The diagnostic uses the same retained
physical wall face flux that is applied to the scalar RHS; it does not
re-evaluate a separate wall closure.

Each output row contains the interval start/end step and time, scalar name and
units, six face contributions, net boundary contribution, volume change,
internal source, residual, tolerance, and status.  Face signs follow the
coordinate-flux convention above, so the net inward boundary contribution is
formed from low-minus-high stored coordinate fluxes.

For total nonprecipitating water,

.. math::

   \Delta\int_\Omega \rho(q_v+q_c)\,dV
   = \int_{t_a}^{t_b}\sum_f J_{v,\mathrm{in}}\,dA\,dt.

SatAdj exchange between ``qv`` and ``qc`` cancels in this total-water row.
With six dry walls the total-water boundary contribution is exactly zero.

Budget statuses are:

``PASS``
   A supported conservation row closes within its printed tolerance.

``FAIL``
   A supported conservation row is outside its printed tolerance.

``UNSUPPORTED_SOURCE``
   The row is not a supported conservation verdict because the diagnostic
   does not account for every required internal source.  In SatAdj Cloud
   Chamber runs, ``rhoTheta`` is always labeled ``UNSUPPORTED_SOURCE`` because
   the current diagnostic omits part of the moist latent-heating source,
   regardless of the numerical residual.  Use dry ``rhoTheta`` for strict
   thermal closure and total nonprecipitating water for moist conservation
   checks.

Enabling the budget diagnostic must not change the simulated state.

Validated regression examples
-----------------------------

The repository contains tracked regression inputs for each major supported
path.  They are useful syntax references, but their grid, timestep, LES
settings, coefficients, and roughness lengths are software-test parameters,
not calibrated physical recommendations.

.. list-table:: Tracked Cloud Chamber examples
   :header-rows: 1
   :widths: 34 66

   * - Case
     - Input file
   * - Dry physical, resolved walls
     - ``Tests/test_files/CloudChamber_Dry/CloudChamber_Dry.i``
   * - SatAdj wet, resolved walls and budgets
     - ``Tests/test_files/CloudChamber_SatAdj_WetBudget/CloudChamber_SatAdj_WetBudget.i``
   * - SatAdj mixed fixed-bulk wall transfer
     - ``Tests/test_files/CloudChamber_SatAdj_BulkMixedWet/CloudChamber_SatAdj_BulkMixedWet.i``
   * - SatAdj neutral roughness-log transfer
     - ``Tests/test_files/CloudChamber_SatAdj_NeutralWet/CloudChamber_SatAdj_NeutralWet.i``
   * - SatAdj MOST wet-wall transfer and budgets
     - ``Tests/test_files/CloudChamber_SatAdj_MOSTWetBudget/CloudChamber_SatAdj_MOSTWetBudget.i``
   * - Mixed neutral sidewalls and MOST horizontal walls
     - ``Tests/test_files/CloudChamber_SatAdj_MOSTMixedWalls/CloudChamber_SatAdj_MOSTMixedWalls.i``
   * - Legacy theta/qv configuration
     - ``Tests/test_files/CloudChamber_Legacy_Config/CloudChamber_Legacy_Config.i``

A complete dry physical example is included below.

.. literalinclude:: ../../Tests/test_files/CloudChamber_Dry/CloudChamber_Dry.i
   :language: none
   :caption: Dry physical Cloud Chamber with resolved wall transfer

Run checklist
-------------

#. Choose ``physical_temperature_rh`` for new work, or ``legacy_theta_qv`` only
   for compatibility with an existing legacy input.
#. Use one ConstantDz, nonperiodic Cartesian level with fixed-density anelastic
   dynamics and gravity.
#. In physical mode, define bottom/top initial temperature and add RH only when
   SatAdj is enabled.
#. Define all six stationary ``NoSlipWall`` faces with physical temperature and
   ``dry`` or ``wet`` moisture state.
#. Choose the momentum, heat, and vapor model on each face.  Remember that one
   ``coefficient_source`` is shared by every bulk channel on that face.
#. Supply every coefficient or roughness length required by the selected
   models.  Check ``0 < z0_* < Delta_n/2`` on the actual grid.
#. For wet heat+vapor MOST on the same face, use equal ``z0_h`` and ``z0_q``.
#. Choose a timestep that satisfies the wall-rate guard, or use adaptive
   stepping.
#. Run a short case and inspect velocity, temperature/potential temperature,
   and, for SatAdj, ``qv`` and ``qc``.
#. Enable ``erf.cloud_chamber_budget_interval`` and check dry thermal or total
   nonprecipitating-water conservation as appropriate.

Supported scope and limitations
-------------------------------

The implemented Cloud Chamber wall-model scope is intentionally limited to:

* one Cartesian ConstantDz level with no AMR refinement;
* six stationary ``NoSlipWall`` faces;
* no periodic direction;
* no terrain, embedded boundaries, or immersed buildings;
* prescribed physical wall temperature in the generalized wall-model path;
* resolved, fixed bulk, neutral roughness-log, and horizontal MOST transfer as
  documented above; and
* forced/shear-dependent bulk transfer with no gustiness or natural/mixed
  convection augmentation.

The implementation does not currently provide moving walls, smooth-wall
``y+``/viscous-sublayer laws, stretched-dz wall sampling, terrain-aware wall
models, EB/AMR wall models, natural-convection wall correlations, ocean/COARE
closures, cloud-water deposition, precipitation, droplet activation, or
facility-specific coefficient calibration.

A numerically stable run is not by itself evidence of quantitative
Pi-Chamber agreement or grid-independent LES behavior.

Troubleshooting
---------------

``initial_relative_humidity`` is rejected in dry mode
   Omit ``prob.initial_relative_humidity`` unless
   ``erf.moisture_model = SatAdj``.  Enter RH as a fraction such as ``0.95``,
   not ``95``.

A wet wall is rejected
   Wet walls require ``erf.moisture_model = SatAdj``.

``coefficient_source`` is rejected
   Use it only when at least one channel on that face is ``bulk_aero``.  The
   same provider applies to every bulk channel on the face.

``C_D``, ``C_H``, or ``C_E`` is missing
   Fixed bulk transfer requires the coefficient corresponding to every bulk
   channel.  A configured dry bulk-vapor channel still requires ``C_E`` even
   though its physical vapor flux is exactly zero.

A ``C_*`` coefficient is rejected with MOST
   MOST computes the bulk coefficients.  Remove user-supplied ``C_D``,
   ``C_H``, and ``C_E`` for the channels using the face-level MOST provider.

MOST is rejected on a side wall
   ``coefficient_source = most`` is implemented only on ``zlo`` and ``zhi``.
   Use resolved, fixed bulk, or neutral roughness-log transfer on side walls.

A roughness length is rejected
   Check that the selected model actually consumes the key and that
   ``0 < z0_* < Delta_n/2`` for the face-normal cell spacing.

Wet heat+vapor MOST rejects unequal ``z0_h`` and ``z0_q``
   The implemented prescribed-temperature bulk-Richardson fixed point uses one
   scalar stability resistance, so coupled wet MOST heat and vapor require
   equal scalar roughness lengths.

Dry vapor appears to have no transfer
   This is intentional.  ``moisture = dry`` enforces exact zero vapor wall
   flux for resolved, fixed bulk, neutral, and MOST selections.  Configured
   model metadata is still validated.

A positive fixed timestep aborts
   Compare the reported ``fixed_dt`` with ``wall_dt``.  Bulk, neutral, and
   MOST wall transfer participate in an explicit wall-rate stability guard.

The cloudy ``rhoTheta`` budget reports ``UNSUPPORTED_SOURCE``
   The current diagnostic does not include the complete moist latent-heating
   source.  Use total nonprecipitating water for the supported cloudy
   conservation check.
