.. _AuxiliaryState:

Auxiliary Prognostic-State Transport
====================================

ERF contains a generic transport substrate for prognostic atmospheric state
that is owned outside the model's core conserved-state array. The purpose of
this layer is to let a separately owned constituent system use ERF's actual
mass carrier, mapped geometry, and host time integration without reproducing
those parts of the dynamical core.

.. note::

   The auxiliary-state layer is currently an implementation interface, not a
   user-selectable generic tracer package. Its live ERF consumer is a
   test-only inert tracer used to qualify the coupling described on this page.
   The spectral-bin microphysics state is not yet advanced through this path.

Why auxiliary state exists
--------------------------

Some atmospheric constituent systems naturally contain many more prognostic
components than belong in ERF's compact core state. Examples include
sectional or spectral particle populations and, potentially, other
multi-component constituent systems. Keeping such state in separately owned
``MultiFab`` objects is useful only if its physical-space transport remains
consistent with the atmospheric flow that transports ERF's native scalars.

The ownership boundary is therefore intentional:

* ERF supplies the host stage timing, dry-air density, mapped geometry, and
  the stage carrier used by native scalar transport.
* The generic auxiliary layer supplies the mapped conservation and
  time-integration seam for separately owned components.
* Each scientific consumer defines the meaning and units of its components,
  constructs its accepted face-transfer rates, and owns its physical source
  terms, constraints, boundary conditions, AMR and restart behavior, and
  coupling projections.

The generic component metadata record a schema identifier together with a
name, semantic identifier, and units for each component. The mapped transport
layer deliberately does not interpret those scientific meanings.

Density-weighted state and mapped conservation
----------------------------------------------

Consider an auxiliary atmospheric quantity whose intensive value per unit
mass of dry air is :math:`z_a`. Its density-weighted state is

.. math::

   U_a = \rho_d z_a,

where :math:`\rho_d` is the dry-air density. The generic transport routines
advance :math:`U_a`; the intensive value is reconstructed from a state and
its matching density only when an operator such as advection requires it.

For the mapped-coordinate convention used by ERF, the cell measure is

.. math::

   \omega = \frac{\det J}{m_x m_y},

where :math:`\det J` is the cell-centered mapping Jacobian and :math:`m_x`
and :math:`m_y` are the horizontal map factors. The corresponding mapped
conservative state is

.. math::

   H_a = \omega U_a.

The current measure builder requires :math:`\det J`, :math:`m_x`,
:math:`m_y`, and :math:`\omega` to be finite and strictly positive in every
valid cell. It also requires the metric fields to have the expected ERF cell
and horizontal-map layouts. A metric or layout failure is rejected rather
than replaced by a default measure.

Mapped face-transfer convention
--------------------------------

The generic stage update consumes a face-centered
:math:`\widetilde{\mathbf F}_a` that is already expressed in ERF's mapped
face-transfer convention. Its computational-coordinate divergence is

.. math::

   D_\xi\left(\widetilde{\mathbf F}_a\right)
   =
   \sum_{d=1}^{3}
   \frac{
      \widetilde F_{a,d,+}
      -
      \widetilde F_{a,d,-}
   }{\Delta \xi_d}.

No additional physical face-area factor and no second application of the map
factors belongs in this divergence. The geometric factors have already
entered through the mapped face-transfer construction and the cell measure
:math:`\omega`.

The reusable divergence/application routine is agnostic to how a valid mapped
face-transfer rate was produced. For example, the qualification tests feed it
both native scalar-advection transfers and mapped transfers derived from ERF's
native scalar-diffusion operators. The current live inert-tracer consumer,
however, constructs advection only; scalar diffusion is deliberately disabled
in that fixture.

Semantic time views and the stage update
----------------------------------------

Auxiliary transport follows the host integrator and does not advance on an
independent timestep. Each stage context identifies three semantic states:

* the **anchor** state at the old-step time;
* the **input** state from which the current stage transfer is constructed;
* the **target** state produced by the current stage.

State, dry-air density, and mapped measure are each carried with these
explicit time identities. This distinction matters when density or geometry
differs between stages. In particular, an intensive quantity must not be
formed from an input-stage auxiliary state and a density belonging to a
different stage.

For recurrence weights :math:`a` and :math:`b` and a face-rate time
coefficient :math:`c`, the generic mapped update is exactly

.. math::

   H_a^{\mathrm{target}}
   =
   a\,\omega^{\mathrm{anchor}} U_a^{\mathrm{anchor}}
   +
   b\,\omega^{\mathrm{input}} U_a^{\mathrm{input}}
   -
   c\,D_\xi\left(\widetilde{\mathbf F}_a\right),

followed by

.. math::

   U_a^{\mathrm{target}}
   =
   \frac{H_a^{\mathrm{target}}}
        {\omega^{\mathrm{target}}}.

The mutable target storage is required to be disjoint from both read-only
input states. The anchor and input views may intentionally refer to the same
storage when the host recurrence permits it. Keeping the target separate
prevents an in-place stage update from becoming dependent on traversal order
when a future transfer construction or limiter uses neighboring values.

Carrier mass flux and intensive transport variables
---------------------------------------------------

When advection is written for an intensive constituent, the input-stage value
is formed from matching semantic views,

.. math::

   z_a^{\mathrm{input}}
   =
   \frac{U_a^{\mathrm{input}}}
        {\rho_d^{\mathrm{input}}}.

The qualification fixture then passes this intensive field to ERF's native
scalar-advection operator together with the host stage carrier. Schematically,
its mapped constituent face-transfer rate is

.. math::

   \widetilde F_{a,f}
   =
   \widetilde F_{d,f}\,z_{a,f},

where :math:`\widetilde F_{d,f}` denotes the carrier supplied by ERF at that
stage and :math:`z_{a,f}` is the native reconstructed face value.

The auxiliary layer does not reconstruct a separate carrier by multiplying a
sampled density and velocity. Using the host carrier is part of the coupling
contract.

One important consistency check is preservation of a constant constituent
ratio. If initially

.. math::

   U_a = k\rho_d

for a spatially constant :math:`k`, then :math:`z_a=k`. When the auxiliary
state and dry-air density are advanced with the same carrier and compatible
host recurrence, the discrete transport should preserve

.. math::

   U_a^{\mathrm{new}} = k\rho_d^{\mathrm{new}}

apart from floating-point roundoff. The auxiliary qualification tests verify
this property with spatially varying density and stage-dependent carrier
fields.

Face-flux rates and completed-step fluxes
-----------------------------------------

A face-flux rate used inside a Runge--Kutta stage and the time-integrated face
flux associated with a completed timestep are different numerical objects.
ERF represents them separately so that predictor-stage intervals cannot be
mistaken for completed-step transport.

See :ref:`TimeAdvance` for ERF's host time integrators.

Compressible three-stage advance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\widetilde{\mathbf F}_0`,
:math:`\widetilde{\mathbf F}_1`, and
:math:`\widetilde{\mathbf F}_2` be the mapped auxiliary face-flux rates at
the three compressible stages. The mapped conservative state follows

.. math::

   H_a^*
   =
   H_a^n
   -
   \frac{\Delta t}{3}
   D_\xi\left(\widetilde{\mathbf F}_0\right),

.. math::

   H_a^{**}
   =
   H_a^n
   -
   \frac{\Delta t}{2}
   D_\xi\left(\widetilde{\mathbf F}_1\right),

and

.. math::

   H_a^{n+1}
   =
   H_a^n
   -
   \Delta t
   D_\xi\left(\widetilde{\mathbf F}_2\right).

The first two stages are predictors anchored at the old state. They are not
three additive forward-Euler increments. Consequently, the completed-step
integrated mapped face flux is

.. math::

   \mathbf J_a
   =
   \Delta t\,\widetilde{\mathbf F}_2.

The first two predictor rates therefore receive zero weight in the
completed-step flux ledger.

Anelastic two-stage Runge--Kutta advance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ERF's two-stage anelastic Runge--Kutta method, the first stage is

.. math::

   H_a^*
   =
   H_a^n
   -
   \Delta t
   D_\xi\left(\widetilde{\mathbf F}_0\right),

and the second stage is

.. math::

   H_a^{n+1}
   =
   \frac{1}{2}H_a^n
   +
   \frac{1}{2}H_a^*
   -
   \frac{\Delta t}{2}
   D_\xi\left(\widetilde{\mathbf F}_1\right).

Equivalently,

.. math::

   H_a^{n+1}
   =
   H_a^n
   -
   \frac{\Delta t}{2}
   D_\xi\left(
      \widetilde{\mathbf F}_0
      +
      \widetilde{\mathbf F}_1
   \right),

so the completed-step integrated mapped face flux is

.. math::

   \mathbf J_a
   =
   \frac{\Delta t}{2}
   \left(
      \widetilde{\mathbf F}_0
      +
      \widetilde{\mathbf F}_1
   \right).

The completed-step ledger enforces stage order, a fixed integrator and
old-step time across the step, compatible face layouts, and the exact temporal
weights above. It supports multiple auxiliary components and applies the host
time weighting to every component.

ERF also provides an anelastic midpoint integrator, but the current auxiliary
stage adapter deliberately rejects that method because it has not yet been
qualified for this transport interface.

Mapped-geometry verification
----------------------------

The generic mapped-divergence convention is tested independently against
ERF's native scalar operators. The current unit evidence includes:

* direct arithmetic checks of the computational mapped divergence;
* telescoping of an arbitrary periodic mapped face field;
* parity with native scalar-advection divergence using nontrivial Jacobian and
  horizontal map factors;
* parity when native scalar-diffusion transfers are supplied on the flat
  native grid and on a vertically stretched grid;
* parity with the static-terrain mapped-diffusion transfer, including the
  terrain metric cross terms;
* preservation of a constant constituent ratio through both supported host
  recurrences;
* exact completed-step temporal weights for one and multiple components;
* rejection of incompatible cell/face layouts, aliased stage targets, and
  invalid or out-of-order stage sequences.

The tests also contain metric negative controls. Omitting :math:`\det J` from
the cell measure, applying a horizontal map factor a second time, using a raw
terrain vertical diffusive flux in place of ERF's mapped transfer, or omitting
the terrain cross terms must disagree with the corresponding native ERF
result. These controls are intended to distinguish the actual ERF metric
convention from an internally consistent but incorrect duplicate convention.

The operator-level tests are broader than the current live proof consumer.
They demonstrate the generic mapped algebra for those supplied transfers; they
do not by themselves establish a complete user-facing tracer capability for
every tested geometry or transport process.

Live inert-tracer proof consumer
--------------------------------

The only integrated ERF consumer of this substrate at the present revision is
a one-component inert tracer used for qualification. It is intentionally
restricted so that the test isolates the coupling seam rather than silently
claiming a general tracer implementation.

The fixture requires:

* three spatial dimensions and exactly one AMR level;
* triply periodic boundaries;
* ``ConstantDz`` geometry with no terrain or buildings;
* moisture disabled and gravity disabled;
* no acoustic substepping;
* no molecular, numerical, LES, RANS, or PBL transport;
* ERF's native scalar-advection path with centered-second-order horizontal and
  vertical dry-scalar advection;
* the two-stage RK method when the host run is anelastic; anelastic midpoint
  is rejected;
* the ``Scalar Advection/Diffusion`` ERF test problem.

Checkpoint/restart, coarse-to-fine initialization, and regrid/remake are
explicitly rejected by this fixture. The fixture also requires a nonzero host
carrier so that the live test actually exercises transport.

The proof consumer uses a time-independent mapped measure. Its storage is
defined first, but :math:`\omega` is built only after ERF has finalized the
level Jacobian and horizontal map factors. Initialization and stage advancement
fail if that measure has not been successfully built. Its old-step, input, and
target measure views intentionally refer to the same static field. A moving
terrain consumer would need distinct geometry at the relevant semantic times
and a separate geometric-conservation qualification.

During the live smoke test, the tracer uses the host stage carrier, produces
stage-varying mapped face rates, checks finite state, and closes the completed
step against the integrated face-flux ledger to roundoff. These diagnostics
are qualification evidence for the fixture; they are not global reductions
that every future production consumer is required to perform on every stage.

Current scope and consumer responsibilities
-------------------------------------------

The reusable auxiliary-state layer is deliberately narrower than a complete
constituent model. In particular, the current generic substrate does not by
itself provide:

* a user-facing runtime registry for arbitrary tracers;
* scientific source or sink terms;
* application-specific positivity, realizability, or composition constraints;
* physical boundary-condition policies for a new constituent family;
* AMR prolongation, restriction, reflux, or restart semantics for a new
  consumer;
* a qualified moving-terrain, embedded-boundary, or acoustically substepped
  auxiliary advance;
* a universal diffusion model or turbulence closure;
* a consumer-specific projection into ERF core coupling fields.

A new consumer must define and qualify those responsibilities explicitly. It
must also ensure that any supplied face-transfer rate follows ERF's mapped
convention and that its state, density, measure, and stage timing are paired
consistently.

The generic architecture is intended to make such consumers possible without
reimplementing ERF's carrier, geometry, or time-integration mathematics. A
sectional aerosol system, atmospheric chemistry state, specialized tracer
family, or another multi-component distribution could in principle use this
interface, but none of those possibilities is advertised here as an existing
runtime capability.

Relationship to spectral-bin microphysics
-----------------------------------------

The SBM spectral state is separately owned outside ERF's core conserved-state
array and is expected to use this auxiliary transport substrate for its
physical-space transport. At the present revision that connection has not yet
been made: ``erf.moisture_model = SBM`` remains a bounded zero-transport
infrastructure configuration.

The auxiliary layer addresses transport between atmospheric grid cells. It
does not define how an SBM particle distribution is represented or remapped in
particle-mass space, nor does it provide SBM realizability limiting,
sedimentation, activation, collision-coalescence, or other cloud
microphysics. Those scientific contracts are documented separately in
:ref:`sec:SpectralBinMicrophysics`.
