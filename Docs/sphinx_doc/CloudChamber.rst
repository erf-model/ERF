.. _CloudChamber:

Cloud Chamber Stage 1
=====================

Stage 1 is a small, equilibrium-cloud proof of concept for a closed Cartesian
rectangular chamber.  It is not quantitatively validated against Pi-Chamber
measurements and it is not an engineering wall model.  The reference geometry
is 2 m x 2 m x 1 m, but any positive rectangular prism is accepted within the
single-level, constant-``dz`` limitations of this stage.

Scope and maturity
------------------

The implementation provides anelastic thermodynamics, physical-temperature/RH
initialization, six independently configured solid walls, resolved-molecular
wall-normal scalar transfer, and SatAdj equilibrium condensate.  It does not
provide MOST, roughness, partial wetness, a wall-water inventory, panel thermal
inertia, conjugate heat transfer, finite-rate droplets, aerosol activation,
precipitation, sedimentation, collision-coalescence, or cloud-water deposition.
Wang et al. is useful motivation for wall sensitivity, but its MOST treatment
is not adopted as a Stage 1 model specification.

Base state and thermodynamics
-----------------------------

Physical chamber mode requires the existing HSE-backed ``ConstantDensity``
initialization path.  ``Uniform`` is rejected because its default
``rho*theta`` state can imply an unintended pressure near 81 kPa.  A typical
reference is:

.. code-block:: none

   erf.init_type = ConstantDensity
   prob.p_inf = 100000.0       # Pa; existing ERF reference-pressure input
   prob.T_0 = 292.0            # K; reference temperature for rho_0

ERF integrates the base pressure hydrostatically from the existing reference
pressure path and reports its minimum, maximum, and lower-boundary value at
startup.  In anelastic mode this base pressure ``p0`` is the thermodynamic
pressure everywhere:

.. math::

   T = \theta (p_0/p_{00})^{R_d/c_p}, \qquad p_{00}=10^5\ {\rm Pa}.

The perturbational pressure is not substituted into this relation.  The same
``p0`` is used by SatAdj, initialization, wet-wall saturation, and the
anelastic ``temp``, ``qsat``, and relative-humidity diagnostics.  Compressible
paths retain their existing EOS pressure.

Physical initialization
------------------------

Select the unambiguous physical contract:

.. code-block:: none

   prob.thermodynamic_initialization = physical_temperature_rh
   prob.initial_temperature_bottom = 300.0       # K
   prob.initial_temperature_top = 284.0          # K
   prob.initial_relative_humidity = 0.95         # fraction, not percent
   prob.temperature_perturbation_amplitude = 0.02 # K
   prob.perturbation_mode = deterministic_sine

The initialized temperature is a linear vertical profile plus a deterministic,
bounded sine perturbation that vanishes on every physical face.  Vapor is
computed from the actual perturbed temperature:

.. math::

   e_v = RH\,e_s(T), \qquad
   q_v = (R_d/R_v)\frac{e_v}{p_0-e_v}, \qquad q_c=0.

Here ``qv`` is ERF's vapor mixing ratio relative to dry air.  RH must be a
finite fraction in ``[0,1]``.  The legacy ``legacy_theta_qv`` mode remains
explicitly numerical and cannot be mixed with the physical keys.

Walls and diffusion
-------------------

Each face is a ``NoSlipWall`` with a physical temperature in kelvin and an
independent moisture mode:

.. code-block:: none

   zlo.type = NoSlipWall
   zlo.temperature = 300.0
   zlo.moisture = wet       # dry or wet
   zlo.wall_transfer_model = resolved_molecular

Dry means impermeable, not ``qv=0``.  Its vapor and cloud-water normal fluxes
are exactly zero.  A wet wall is a pure-liquid equilibrium surface with
``RH=1`` and ``qv_wall = qsat(T_wall,p0_wall)``.  The local base pressure is
used on every face, including the height-dependent pressure on vertical walls.

For ``ConstantAlpha`` and explicit vertical diffusion, the physical wall flux
uses the molecular coefficient and the half-cell distance.  The interior may
still contain SGS diffusion, but SGS diffusivity is excluded from the physical
wall-normal transfer.  Low-face coordinate flux is positive into the domain
when it is positive; high-face inward flux is the negative of the coordinate
flux.  Both evaporation and condensation are supported.  Cloud water has zero
wall flux in Stage 1.

SatAdj and diagnostics
----------------------

SatAdj is a fixed-pressure, cell-local adjustment.  It conserves total
non-precipitating water ``qv+qc`` and its existing latent invariant, maintains
nonnegative vapor/cloud water, and is idempotent after adjustment.  It is an
equilibrium cloud approximation: ``qc`` is not a droplet population and does
not imply rain or aerosol physics.

For anelastic plotfiles, ``temp`` and ``qsat`` use ``p0`` and ``rel_humidity``
is emitted as a dimensionless fraction.  The consistency checks are:

.. math::

   T=\theta(p_0/p_{00})^{R_d/c_p}, \qquad
   q_{sat}=q_{sat}(T,p_0), \qquad RH=e_v/e_s(T).

Budgets
-------

The optional ``erf.cloud_chamber_budget_interval`` writes interval-local
balances for ``rhoTheta`` (a potential-temperature scalar, not watts), vapor,
cloud water, and total non-precipitating water.  Fluxes are the exact applied
coordinate-face fluxes, with low-minus-high domain orientation.  Each report
resets its reference state and accumulators, so successive intervals and
restart-started intervals are interpreted locally.

For all-dry walls, total water should be constant apart from the reported
floating-point residual.  For wet walls, its change equals the time-integrated
inward vapor flux; SatAdj internal vapor/cloud-water exchanges cancel in the
total-water row.

Reference validation
--------------------

The supplied dry regression uses physical initialization, six dry walls, and
``ConstantAlpha`` diffusion to test closed water.  The SatAdj regression uses
warm wet bottom/cold wet top walls and dry side walls.  Its checker independently
verifies positive base pressure, the physical temperature profile, the
anelastic theta relation, the exact RH conversion, zero initial cloud water,
zero initial velocity, finite nonnegative moisture, and early buoyant motion.
It does not require a particular cloud-water amount in a short CI run.

The reference configuration is a proof of concept, not a claim of realistic
coarse-grid molecular transfer or quantitative Pi-Chamber validation.
