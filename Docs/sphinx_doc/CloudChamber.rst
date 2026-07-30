.. _CloudChamber:

Cloud Chamber Stage 1
=====================

The Stage 1 ``Cloud Chamber`` problem is a reproducible proof of concept for
a closed, Cartesian rectangular chamber.  It is intended to exercise ERF's
existing anelastic dynamics, six solid walls, explicit scalar diffusion, LES,
and optional SatAdj moisture path.  It is not a quantitative model of the
physical Pi-Chamber and does not implement a panel heat-transfer law, wet-wall
evaporation, Wang/MOST wall closure, aerosol activation, or precipitation.

The initializer accepts any positive rectangular prism.  The reference decks
use a 2 m x 2 m x 1 m box, one AMR level, and nonperiodic boundaries.  All six
faces must be ``NoSlipWall`` and must provide numerical prescribed
potential-temperature values.  In a cloudy deck, all six faces must also
provide numerical prescribed vapor mixing ratios.  These values retain the
existing ERF face namespaces (``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, and
``zhi``); they are not reinterpreted as physical temperature, relative
humidity, roughness, or a transfer coefficient.

Input contract
--------------

The minimum dry contract is:

.. code-block:: none

   erf.prob_name = "Cloud Chamber"
   erf.init_type = Uniform
   erf.anelastic = 1
   erf.use_gravity = true
   geometry.prob_lo = 0.0 0.0 0.0
   geometry.prob_hi = 2.0 2.0 1.0
   geometry.is_periodic = 0 0 0
   amr.n_cell = 16 16 8
   amr.max_level = 0

   prob.theta_bottom = 299.0
   prob.theta_top = 280.0
   prob.theta_perturbation_amplitude = 0.05
   prob.perturbation_mode = deterministic_sine

   xlo.type = NoSlipWall
   xhi.type = NoSlipWall
   ylo.type = NoSlipWall
   yhi.type = NoSlipWall
   zlo.type = NoSlipWall
   zhi.type = NoSlipWall
   xlo.theta = 285.0
   xhi.theta = 285.0
   ylo.theta = 285.0
   yhi.theta = 285.0
   zlo.theta = 299.0
   zhi.theta = 280.0

For the cloudy variant, add ``erf.moisture_model = SatAdj``,
``prob.qv_bottom`` and ``prob.qv_top``, and ``qv`` on every face.  The
reference cloudy deck deliberately starts supersaturated so the existing
equilibrium adjustment forms cloud water.  SatAdj transports vapor and cloud
water only; it does not imply rain, sedimentation, aerosol, droplet number,
collision-coalescence, or particle-wall physics.

Initialization contract
-----------------------

ERF first constructs its hydrostatic/base state.  The chamber hook then adds
conserved corrections to the existing state.  At each cell center,

.. math::

   \theta = \theta_\mathrm{bottom}
          + (\theta_\mathrm{top}-\theta_\mathrm{bottom})
            (z-z_\mathrm{lo})/L_z
          + A\sin(2\pi x/L_x)\sin(2\pi y/L_y)\sin(\pi z/L_z).

The sine perturbation is bounded, smooth, vanishes on the two vertical faces,
and uses no random-number generator.  Its continuous volume mean is zero;
the cell-center discrete mean is expected to be zero within floating-point
and mesh-layout tolerance.  The initializer writes
``state_pert(RhoTheta) = rho*theta - state(RhoTheta)``.  In SatAdj mode it
similarly targets ``rho*qv`` in ``RhoQ1`` and targets zero ``RhoQ2`` before
the first adjustment.  Density is never changed in the anelastic path.

Wall-model seam
---------------

The active Stage 1 treatment is the existing prescribed-state/no-slip path.
The chamber's compact face metadata separates face identity, prescribed
properties, treatment tag, and future flux diagnostics without adding virtual
dispatch to device kernels.  A later host-selected wall model can consume the
face orientation, resolved state, and near-wall samples and return momentum,
thermal-scalar, and vapor fluxes plus validity flags.  Wang/MOST is a possible
future comparison model, not a Stage 1 default or an implemented closure.

Testing
-------

``CloudChamber_Dry`` and ``CloudChamber_SatAdj`` are checker-driven short
regressions.  Their checker validates the initial theta profile, the derived
physical ``temp`` field against the EOS value implied by that initialized
state, zero initial velocity, finite evolved fields, early buoyant response,
and, for SatAdj, the initial qv target plus finite nonnegative vapor/cloud
water with positive cloud water somewhere in the final snapshot.  The focused
``CloudChamberProfile`` unit tests protect endpoint interpolation, zero
amplitude, boundedness, and the vertical-face zeros.  Run the regressions
with the same MPI launcher used by the build; one and two ranks should agree
within floating-point reduction tolerance.
