
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _GoverningEquations:

Governing Equations
=============================

ERF solves the fully compressible, the anelastic or the
pseudo-incompressible equations.
All formulations predict velocity, dry potential temperature, and any
moisture or other scalars. They differ in their treatment of density and
pressure.

In **compressible** mode, ERF advances dry-air density
:math:`\rho_d`, face-centered dry momentum :math:`\rho_d\mathbf{u}`,
and :math:`\rho_d\theta_d`. Pressure is diagnosed from the prognostic
thermodynamic state using the equation of state.

In **anelastic** mode, ERF holds dry density at the reference value
:math:`\rho_0(z)` and enforces
:math:`\nabla\cdot(\rho_0\mathbf{u})=0` with a pressure
projection. Buoyancy remains a source in vertical momentum, but its
pressure-gradient term is not recomputed from the compressible equation
of state.

In **pseudo-incompressible** mode, ERF evolves the density and
:math:`\rho_d\theta_d` as in compressible mode, and enforces
:math:`\nabla\cdot(\rho_0\theta_0\mathbf{u})=0` with a
variable-coefficient pressure projection; see
:ref:`PseudoIncompressibleEquations`.

Here :math:`T`, :math:`\theta_d`, and :math:`p` denote temperature,
dry potential temperature, and pressure at cell centers, while
:math:`\mathbf{u}` and dry momentum are represented on cell faces.
Moisture mixing ratios are defined per unit *dry-air* mass. See
:ref:`Buoyancy` for the gravitational force and reference-state
conventions.

Compressible Equations
------------------------

The first three equations governing fully compressible flow are

.. math::
   \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
                                                        \frac{1}{1 + q_t} ( \nabla p^{\prime}  - \hat{\boldsymbol{z}} B_z ) -
                                                        \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_d \theta_d)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta_d) +
                                                      \nabla \cdot (\rho_d \alpha_{\theta}\ \nabla \theta_d) +
                                                      F_{\theta} + H_{n} + H_{p},

supplemented with the equation of state as given below.

Anelastic Equations
------------------------

The first two equations for the anelastic formulation are

.. math::
   \frac{\partial (\rho_0 \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} \mathbf{u}) -
                                                        \frac{1}{1 + q_t} ( \nabla p^\prime - \hat{\boldsymbol{z}} B_z ) -
                                                        \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_0 \theta_d)}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} \theta_d) +
                                                      \nabla \cdot (\rho_0 \alpha_{\theta}\ \nabla \theta_d) +
                                                      F_{\theta} + H_{n} + H_{p},

supplemented with the constraint

.. math::
  \nabla \cdot (\rho_0 \mathbf{u}) = 0

In these momentum equations :math:`B_z` is an upward-positive
vertical **force density**, with units of :math:`\mathrm{N\,m^{-3}}`.
The explicit slow-momentum source divides the stored pressure-gradient
and buoyancy terms by :math:`1+q_t`, with water mixing ratio averaged
from the two neighboring cells to the face. Anelastic projection also
adds a separate correction to the momenta; the equation here is a
continuum summary, not a complete specification of the projection step.

For **compressible** flow, the pressure :math:`p` is diagnosed from the
prognostic thermodynamic state using the equation of state. The
perturbational pressure is :math:`p'=p-p_0`, where :math:`p_0(z)` is
the hydrostatic reference pressure. The momentum equations above show
the usual perturbational-pressure form.

ERF also offers a choice in the **discrete horizontal pressure gradient**.
For ``erf.gradp_type = 0``, the default
``erf.use_pert_pres_gradient = true`` evaluates the horizontal
gradients from :math:`p'`. Setting
``erf.use_pert_pres_gradient = false`` instead evaluates them from the
full EOS pressure :math:`p`, before subtracting :math:`p_0` for the
vertical pressure gradient. The vertical gradient is evaluated from
:math:`p'` in both cases. With ``erf.gradp_type = 1``, ERF uses its
interpolated perturbational-pressure gradient rather than this
full-pressure option. Because :math:`p_0` depends only on physical
height, the continuum horizontal derivatives of :math:`p` and
:math:`p'` agree at fixed height; the discrete choices can nevertheless
differ on terrain-following grids. See :ref:`sec:Inputs`.

For **anelastic** flow, the :math:`p'` appearing in the schematic
momentum equation denotes a pressure-like contribution enforced by
projection. Its gradient is maintained by the projection solver, not
recomputed as an independently diagnosed compressible EOS pressure
perturbation. The projection also updates the momenta separately.
See :ref:`Buoyancy` for the gravitational source definitions and
face-average conventions.

.. _PseudoIncompressibleEquations:

Pseudo-Incompressible Equations
--------------------------------

The pseudo-incompressible formulation (Durran, 1989) sits between the
compressible and anelastic formulations.  Like the compressible
equations it evolves the density; like the anelastic equations it
filters acoustic waves by imposing a constraint on the velocity and
solving for the pressure that enforces it.  The constraint is

.. math::
  \nabla \cdot (\rho_0 \theta_0 \mathbf{u}) = 0,

where :math:`\rho_0(z)` and :math:`\theta_0(z)` are the density and dry
potential temperature of the hydrostatic base state.  The anelastic
constraint is recovered with :math:`\theta_0` constant.  With heating or
diffusion the right-hand side is a source term :math:`S` (see **Source
term** below).

**Projection.**  Both the anelastic and pseudo-incompressible
projections can be written in terms of a face-centered weight
:math:`\sigma`, with :math:`\sigma = \rho_0` (anelastic) or
:math:`\sigma = \rho_0 \theta_0` (pseudo-incompressible).  The
momentum update from the pressure is

.. math::
  (\rho \mathbf{u})^{n+1} = (\rho \mathbf{u})^{*} - \Delta t \, \nabla \phi ,

so the projected quantity
:math:`\tilde{\mathbf{m}} = \beta\,(\rho \mathbf{u})`, with
:math:`\beta = \sigma / \rho` on each face, satisfies
:math:`\tilde{\mathbf{m}}^{n+1} = \tilde{\mathbf{m}}^{*} - \Delta t\,\beta\,\nabla\phi`.
Imposing :math:`\nabla \cdot \tilde{\mathbf{m}}^{n+1} = 0` gives the
variable-coefficient Poisson equation

.. math::
  \nabla \cdot \left( \beta \nabla \phi \right) = \frac{1}{\Delta t} \nabla \cdot \tilde{\mathbf{m}}^{*} .

In the anelastic formulation :math:`\rho = \rho_0`, so
:math:`\beta = 1` and the Poisson equation has constant coefficients,
which allows the FFT solvers.  In the pseudo-incompressible
formulation :math:`\beta = \rho_0 \theta_0 / \rho`, which varies in
space and time; since :math:`\rho\theta \approx \rho_0\theta_0`,
:math:`\beta \approx \theta`, the coefficient of the Exner-pressure
gradient in Durran's form of the momentum equation.  ERF never assumes
:math:`\rho = \rho_0` in the pseudo-incompressible projection: the
momenta are multiplied by :math:`\beta` on each face before the solve,
and divided by it afterwards.  The same face values of :math:`\beta` are
used in the conversion, in the Poisson operator and in the momentum
correction, so the projected momenta are discretely divergence-free.
The stored pressure gradient is :math:`\nabla\phi / \Delta t` itself,
without :math:`\beta`.

**Face values of** :math:`\beta`.  The coefficient is the base-state
quantity :math:`\rho_0\theta_0/\rho`; how it is put onto the faces is a
choice of discretization.  ERF maps the cell-centered field
:math:`\rho_0\theta_0/\rho` onto each face with the reconstruction that
the dycore advection scheme (``erf.dycore_horiz_adv_type``,
``erf.dycore_vert_adv_type``) uses to put :math:`\theta` onto that face
in the advective flux of :math:`\rho\theta`, which is the mass flux
times the face value of :math:`\theta`.  On a terrain-fitted mesh the
:math:`w`-face :math:`\beta` multiplies the contravariant mass flux
:math:`\Omega(\rho\mathbf{u})`, as :math:`\theta` does in the flux.
Where :math:`\rho\theta = \rho_0\theta_0`, the face value of
:math:`\rho_0\theta_0/\rho` is then the face value of :math:`\theta`, so
the projected flux :math:`\beta\,(\rho\mathbf{u})` is the advective
flux of :math:`\rho\theta`, and the projection keeps
:math:`\rho\theta = \rho_0\theta_0` without heating or diffusion.  On an
embedded-boundary mesh :math:`\beta` is instead the ratio of the face
averages of :math:`\rho_0\theta_0` and :math:`\rho`.

**Upwind direction.**  An upwind or WENO scheme picks its stencil from
the direction of the flow on each face, and the projection can reverse
that direction.  On a pseudo-incompressible level the direction on each
face is therefore fixed by the projection: it is the sign of the
momentum going into the projection (zero, which selects the centered
blend, where that momentum is zero).  The projection maps
:math:`\rho_0\theta_0/\rho` onto the faces with the stencil this
direction selects, and the advection of every cell-centered scalar
from the projected state (:math:`\rho\theta` with the dycore scheme,
the moist and dry scalars with theirs) picks its upwind side from the
same direction rather than from the sign of the projected mass flux.
The face value of :math:`\theta` in the advective flux of
:math:`\rho\theta` is then the face value of :math:`\beta` whatever the
scheme.  Where the projection reverses the flow on a face, the stencil
there lies downwind of the projected flow, so the scheme's upwind
dissipation acts with the opposite sign on that face; such faces have
small velocities.  The direction is kept between the stage that projects
and the stage that advects (and in the checkpoint).  The momenta are
advected as on any other level, with the direction of their own mass
fluxes.

**Ghost cells.**  Near the edge of a level the stencils of the upwind
and higher-order schemes read cells two or three deep outside it, where
the advection reads :math:`\theta = (\rho\theta)/\rho` as the fill of
the state sets :math:`\rho` and :math:`\rho\theta` there.  The
projection gives :math:`\rho_0\theta_0` the same ghost cells: outside
the domain it applies the boundary conditions of :math:`\rho\theta`, and
at a coarse/fine boundary it interpolates the coarse
:math:`\rho_0\theta_0` as the fill interpolates the coarse
:math:`\rho\theta`; it then divides by the ghost-cell :math:`\rho` of the
state.  With :math:`\rho\theta = \rho_0\theta_0` the ratio is then
:math:`\theta` in those cells too.

**Solvers.**  Every pseudo-incompressible projection on a mesh without
embedded boundaries uses the terrain Poisson operator with :math:`\beta`
folded into the face area factors, on uniform, stretched and
terrain-fitted meshes alike.  On a terrain-fitted mesh
``erf.terrain_poisson_solver`` selects the solver as for an anelastic
level; the FFT preconditioner of ``gmres_fft`` has constant coefficients,
so GMRES may take more iterations where :math:`\beta` varies.  A flat
(uniform or stretched) mesh always uses ``amrex::MLTerrainPoisson`` with
multigrid, and needs no FFT build.  On an embedded-boundary mesh the
multigrid :math:`b` coefficient is :math:`\beta`.  The
constant-coefficient FFT and multigrid solvers are used only for
anelastic levels.  On a terrain-fitted mesh the
correction to :math:`w` is recovered from the corrected contravariant
momentum, so it equals :math:`\beta\,\partial\phi/\partial z` only to
second order in the mesh spacing where :math:`\beta` varies.

**Density, potential temperature and buoyancy.**  A
pseudo-incompressible level evolves the density with the continuity
equation and :math:`\rho\theta_d` with its own conservative equation,
both in the flux form of the compressible equations,

.. math::
   \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \theta_d)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta_d) +
                                                      \nabla \cdot (\rho_d \alpha_{\theta}\ \nabla \theta_d) +
                                                      F_{\theta} + H_{n} + H_{p},

so ``fixed_density = 0`` there.  With AMR, :math:`\rho` and
:math:`\rho\theta_d` are refluxed and averaged down, so both are conserved
over the levels.  The buoyancy is the density departure,
:math:`-(\rho - \rho_0)\,g` (``buoyancy_type = 1``, which is forced at a
pseudo-incompressible level).  The thermodynamic pressure is the base-state
pressure :math:`p_0`, as at an anelastic level, so temperature, saturation
and the microphysics use :math:`p_0`.

A fresh start sets :math:`\rho = \rho_0\theta_0/\theta` in every cell of a
pseudo-incompressible level, keeping :math:`\theta` and the other mixing
ratios, so that the initial state satisfies
:math:`\rho\theta = \rho_0\theta_0`, which the constraint, with the face
values of :math:`\beta` above, maintains to the tolerance of the solve
without heating or diffusion.  With ``erf.v = 1`` and
``erf.sum_interval``, ``PI DRIFT`` reports the maximum of
:math:`|\rho\theta - \rho_0\theta_0| / (\rho_0\theta_0)` over the
pseudo-incompressible levels, leaving out cells covered by a finer level
(``erf.v = 2`` also prints it for each level).

**Source term.**  With heating or diffusion, :math:`\rho\theta` changes
by its non-advective tendency :math:`Q` (the diffusion, including any
implicit vertical diffusion, and the sources of ``cc_src`` such as
radiation, nudging and sponges), so keeping
:math:`\rho\theta = \rho_0\theta_0` requires

.. math::
   \nabla \cdot (\rho_0\theta_0 \mathbf{u}) = S, \qquad
   S = Q + c\, \frac{\rho\theta - \rho_0\theta_0}{\Delta t},

discretely :math:`\nabla\cdot(\beta\,\mathbf{m}) = S` with the divergence
of the advective update of :math:`\rho\theta`.  The projection at the end
of each stage sets the velocity with which the next stage advects
:math:`\rho\theta`, whose :math:`Q` is not yet known, so it uses the
:math:`Q` of the stage just taken: the whole :math:`\rho\theta`
tendency of that stage less its advective part (with the RK2 integrator,
the average over the two stages where the tendency is averaged).  The
second term, with :math:`c` = ``erf.pi_drift_factor`` (default 0.5) and
:math:`\Delta t` the level's time step, pulls back what the lag of
:math:`Q` leaves, together with any other departure from
:math:`\rho_0\theta_0`; values of :math:`c` above about 1 can be
unstable.  The initial projection takes :math:`S = 0`, and :math:`Q` is
kept in the checkpoint.  On the stratified hill of
``Tests/test_files/PseudoIncompressible_StratifiedHill`` with
:math:`\alpha_T = 10\ \mathrm{m^2/s}`, the drift after 10 steps is
:math:`10^{-3}` with :math:`S = 0`, :math:`2 \times 10^{-5}` with
:math:`S = Q` (it keeps growing, since the lag of :math:`Q` is first
order in time) and :math:`5 \times 10^{-6}` with :math:`c = 0.5` (it
stays at a few :math:`10^{-6}` over 50 steps).  Latent heating from the
microphysics, which acts after the dynamics, is not in :math:`Q`; only
the drift correction acts on it.

**Current status.**  The following are not yet implemented:

- latent heating in the source term (with moist pseudo-incompressible
  flow below), and the source term on an embedded-boundary mesh,
- face values of :math:`\beta` that match the advection of
  :math:`\rho\theta` on an embedded-boundary mesh (there the scalars
  are upwinded with the sign of the projected mass flux),
- moist pseudo-incompressible flow (the constraint uses the dry
  :math:`\theta_0`),
- matching the coarse and fine fluxes :math:`\beta\,(\rho\mathbf{u})`
  through coarse/fine faces in time (the fine-level solve removes the
  mismatch as a uniform divergence) and keeping
  :math:`\rho\theta = \rho_0\theta_0` when the fine level is refluxed and
  averaged down: on refined levels :math:`\rho\theta` still drifts from
  :math:`\rho_0\theta_0`.  In ten steps on the stratified hill of
  ``Tests/test_files/PseudoIncompressible_StratifiedHill`` the drift is
  about :math:`3 \times 10^{-6}` with one-way coupling (:math:`10^{-5}`
  without the drift correction) and :math:`10^{-4}` with two-way
  coupling, against :math:`10^{-10}` on a single level.  With two-way
  coupling the drift is in the coarse cells next to the fine level, which
  the reflux changes every coarse step by about as much as the drift
  correction removes; the upwind direction at the coarse/fine faces is
  not yet made consistent between the levels either,
- the implicit immersed-forcing drag in the projection
  (``erf.if_implicit_projection``),
- a refined level that stops inside the domain in :math:`z` on a
  terrain-fitted mesh (refinement through the full depth is supported),
- a refined level of a flat mesh with map factors other than one (as
  for an anelastic level),
- density boundary conditions on walls, the WRF boundary density, the
  cloud chamber and spectral-bin microphysics, which assume a fixed density.

(Dry and Moist) Scalars
-----------------------

We supplement the above equations with the following equations for advected scalars (:math:`\phi`) and
precipitating (:math:`\mathbf{q_{p}}`) and non-precipitating (:math:`\mathbf{q_{n}}`)
moisture variables (identical for compressible and anelastic)

.. math::
   \frac{\partial (\rho_d \boldsymbol{\phi})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \boldsymbol{\phi}) + \nabla \cdot ( \rho_d \alpha_{\phi}\ \nabla \boldsymbol{\phi}) + \mathbf{F}_{\phi},

   \frac{\partial (\rho_d \mathbf{q_{n}})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{n}}) + \nabla \cdot (\rho_d \alpha_{q} \nabla \mathbf{q_{n}}) + \mathbf{F_{n}} + \mathbf{G_{p}},

   \frac{\partial (\rho_d \mathbf{q_{p}})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{p}}) + \partial_{z} \left( \rho_d \mathbf{w_{t}} \mathbf{q_{p}} \right) + \mathbf{F_{p}}.

The non-precipitating water mixing ratio vector :math:`\mathbf{q_{n}} = \left[ q_v \;\; q_c \;\; q_i \right]` includes water vapor, :math:`q_v`, cloud water, :math:`q_c`, and cloud ice, :math:`q_i`, although some microphysical moisture models may not include cloud ice; similarly, the precipitating water mixing ratio vector :math:`\mathbf{q_{p}} = \left[ q_r \;\; q_s \;\; q_g \right]` involves rain, :math:`q_r`, snow, :math:`q_s`, and graupel, :math:`q_g`, though some models may not include these terms. The source terms for moisture variables, :math:`\mathbf{F_{p}}`, :math:`\mathbf{F_{n}}`, :math:`\mathbf{G_{p}}`, and their corresponding impact on potential temperature, :math:`H_{n}` and :math:`H_{p}`, and the terminal velocity, :math:`\mathbf{w_{t}}` are specific to the employed model.
See the :ref:`Microphysics<Microphysics>` section for more details.

Height-Following Terrain Coordinates
------------------------------------
Consider two coordinate systems that correspond to a terrain-following grid, :math:`\mathbf{X}`, and a flat cartesian grid, :math:`\mathbf{Z}`, with axes given by

.. math::
   \mathbf{X} = \left[ x \; y \; z \right]^{\intercal}, \quad \quad \mathbf{\Xi} = \left[ \xi \; \eta \; \zeta \right]^{\intercal},

and

.. math::
   x = \xi, \quad \quad y = \eta, \quad \quad z =  h \left(\xi, \, \eta, \, \zeta \right).

Only the vertical coordinate in the physical domain is deformed by the terrain-fitting.
To account for isotropic lateral grid stretching as represented by "map factors" :math:`m_x = m_y = m` as in WRF, we augment the coordinate transform above with stretching in the lateral directions only.

These combined transformations yield the following Jacobian, :math:`\bar{\mathbf{J}}`, and inverse Jacobian, :math:`\bar{\mathbf{T}}`, matrices

.. math::
    \bar{\mathbf{J}}  = \begin{bmatrix}
    \frac{1}{m} & 0 & 0 \\
    0 & \frac{1}{m} & 0\\
   h_{\xi} &  h_{\eta} & h_{\zeta} \\
    \end{bmatrix}, \quad \quad
     \bar{\mathbf{T}} =  \mathbf{J}^{-1} =  \frac{m^2}{h_{\zeta}} \begin{bmatrix}
    \frac{h_{\zeta}}{m} & 0 & 0 \\
    0 & \frac{h_{\zeta}}{m} & 0\\
   -\frac{h_{\xi}}{m} &  -\frac{h_{\eta}}{m} & \frac{1}{m^2} \\
  \end{bmatrix}
  =
   \begin{bmatrix}
    m & 0 & 0 \\
    0 & m & 0\\
   -\frac{h_{\xi}}{h_\zeta}m &  -\frac{h_{\eta}}{h_\zeta}m & \frac{1}{h_\zeta} \\
    \end{bmatrix}.

In the above, :math:`J = \left| \bar{\mathbf{J}} \right |=  h_{\zeta} / m^2` is the Jacobian determinant. To explicitly close the governing equations in terrain-following coordinates, we provide relations for the gradient of a scalar (:math:`f`) and divergence of a vector (:math:`\mathbf{F}`):

.. math::
    \nabla_{\mathbf{X}} f &= \bar{\mathbf{T}}^{\intercal} \nabla_{\mathbf{Z}} f,

    \nabla_{\mathbf{X}} \cdot \left( \mathbf{F} \right) &= \frac{1}{J} \nabla_{\mathbf{Z}} \cdot \left( J  \bar{\mathbf{T}} \mathbf{F}\right).


Vector rotation of the fluid velocity yields :math:`J  \bar{\mathbf{T}} \mathbf{u} = \left[h_{\zeta}u/m, \;\; h_{\zeta}v/m, \;\; \omega/m^2  \right]^{\intercal}`, where :math:`\omega = w -h_{\xi} u m - h_{\eta} v m` is the vertical velocity that is normal to the top/bottom faces of the grid cells.


Background (reference) state
-----------------------------

For compressible flow, reference profiles depend on height and define
thermodynamic perturbations:

.. math::

   p=p_0(z)+p',\qquad \rho_d=\rho_0(z)+\rho_d'.

Here :math:`\rho_0` is the **dry** base-state density. In moist flow,
:math:`q_{v0}` is the base-state water-vapor mixing ratio and the
reference state contains no condensate. With positive-upward :math:`z`
and gravity magnitude :math:`g>0`, hydrostatic balance is

.. math::

   \frac{dp_0}{dz}=-\rho_0(1+q_{v0})g.

The same reference density and hydrostatic pressure are used by the
anelastic formulation, which fixes dry density at :math:`\rho_0` and
obtains its momentum pressure-gradient term from projection. The fixed
thermodynamic reference pressure :math:`P_{00}=10^5\,\mathrm{Pa}` used in
potential-temperature and EOS definitions is **not** the varying
hydrostatic profile :math:`p_0(z)`.

Equation of state (compressible only)
--------------------------------------

In the fully compressible formulation, the total pressure is computed as

.. math::
  p = P_{00} \left( \frac{R_d \rho_d \theta_m}{P_{00}} \right)^\Gamma

Here :math:`\Gamma=1.4` is the fixed EOS exponent, consistent with
the reference dry-air heat capacity
:math:`C_{p,d}=1004.5\,\mathrm{J\,kg^{-1}\,K^{-1}}`.
For the EOS, define

.. math::
  \theta_m = \theta_d (1 + \frac{R_v}{R_d} q_v)

Here :math:`\theta_m` is an algebraic shorthand for the water-vapor
factor used by the compressible EOS, not a separate prognostic
thermodynamic variable. In particular, it should not be confused with
virtual potential temperature, whose moist-air density relation also
accounts for the mass of condensed water. ERF evolves dry potential
temperature :math:`\theta_d` through the conserved variable
:math:`\rho_d\theta_d`. In the equations above, :math:`R_d` is the
dry-air gas constant and :math:`P_{00}=10^5\,\mathrm{Pa}` is the fixed
reference pressure; neither should be confused with the varying
hydrostatic reference pressure :math:`p_0(z)`.

Define

.. math::

   M(q_v) = 1 + \frac{R_v}{R_d} q_v.

Then the EOS can be written equivalently as

.. math::

   p = P_{00}
       \left(
       \frac{R_d \rho_d \theta_d M(q_v)}{P_{00}}
       \right)^\Gamma.

Here :math:`q_v` is the vapor mixing ratio per unit dry-air mass. The quantity
:math:`\rho_d \theta_d` is the ``rhotheta`` state variable used by the EOS
utility functions, and EOS pressure arguments and return values are in Pa.

The same EOS implies

.. math::

   p = \rho_d R_d T M(q_v),

and therefore

.. math::

   \rho_d =
   \frac{p}{R_d T M(q_v)}.

This matches :cpp:`getRhogivenTandPress` in ``Source/Utils/ERF_EOS.H``.

At the default EOS constants, the thermodynamic relation is

.. math::

   \kappa_{\mathrm{EOS}} = \frac{\Gamma - 1}{\Gamma}
   = \frac{R_d}{C_{p,d}},
   \qquad
   \frac{1}{\Gamma} = 1 - \kappa_{\mathrm{EOS}}.

The separate runtime parameter ``erf.c_p`` is used by other ERF
thermodynamic and forcing pathways; changing it does **not** change the
fixed exponent in :cpp:`getPgivenRTh`. EOS relations that combine the
runtime :math:`R_d/c_p` with fixed :math:`\Gamma` therefore need not
remain mutually inverse when ``erf.c_p`` differs from :math:`C_{p,d}`.
A coordinated thermodynamics change is required before such combinations
can be treated as fully EOS-consistent.

The functions in ``Source/Utils/ERF_EOS.H`` implement these relations:

.. list-table::
   :header-rows: 1
   :widths: 30 50

   * - Function
     - Relation
   * - :cpp:`getThgivenTandP`
     - :math:`\theta_d = T(P_{00}/p)^{R_d/c_p}`
   * - :cpp:`getTgivenPandTh`
     - :math:`T = \theta_d(p/P_{00})^{R_d/c_p}`
   * - :cpp:`getPgivenRTh`
     - :math:`p = P_{00}(R_d\rho_d\theta_d M/P_{00})^\Gamma`
   * - :cpp:`getRhoThetagivenP`
     - inverse of :cpp:`getPgivenRTh`
   * - :cpp:`getRhogivenTandPress`
     - :math:`\rho_d = p/(R_d T M)`
   * - :cpp:`getRhogivenThetaPress`
     - density from :math:`\theta_d`, :math:`p`, and :math:`q_v`; its runtime :math:`R_d/c_p` and fixed :math:`\Gamma` factors need not invert :cpp:`getPgivenRTh` for nondefault ``erf.c_p``
   * - :cpp:`getExnergivenP`
     - :math:`\Pi = (p/P_{00})^{R_d/c_p}`
   * - :cpp:`getdPdRgivenConstantTheta`
     - :math:`(\partial p/\partial\rho_d)_{\theta_d,q_v} = \Gamma p/\rho_d`

Additional terms
--------------------------------------

- :math:`\boldsymbol{\tau}` is the viscous stress tensor,

  .. math::
     \tau_{ij} = -2\mu \sigma_{ij},

with :math:`\sigma_{ij} = S_{ij} -D_{ij}` being the deviatoric part of the strain rate, and

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right), \hspace{24pt}
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij},

- :math:`\mathbf{F}_{u}` and :math:`F_{\theta_d}` are the forcing terms described in :ref:`Forcings`,
- :math:`B_z=g_z[\rho_d(1+q_t)-\rho_0(1+q_{v0})]` is the exact
  density-perturbation gravitational force density used by compressible
  buoyancy type 1; other supported configurations use the approximate
  vertical force densities described in :ref:`Buoyancy <Buoyancy>`,
- :math:`\boldsymbol{g}=(0,0,g_z)` with :math:`g_z=-g` and :math:`g>0`
  is the downward gravity vector,
- The dry potential temperature :math:`\theta_d` is defined from temperature :math:`T`, pressure :math:`p`, and reference pressure :math:`P_{00} = 10^{5}` Pa as

.. math::

  \theta_d = T \left( \frac{P_{00}}{p} \right)^{R_d / c_p}.

(In the anelastic case, :math:`p` is replaced by :math:`p_0` in the relationship between :math:`\theta_d` and :math:`T`.)


Assumptions
------------------------

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior with constant specific heats (:math:`c_p,c_v`). In dry configurations,
  :math:`p = \rho_d R_d T`. In moist configurations, ERF uses dry density and vapor
  mixing ratio per dry-air mass:

  .. math::

     p = \rho_d R_d T
         \left(1 + \frac{R_v}{R_d}q_v\right).

  This is equivalent to the moist potential temperature factor in the
  compressible EOS above. The dry shorthand :math:`p = \rho R_d T` is recovered
  only when :math:`q_v = 0`.
- Viscous heating is negligible
- No chemical reactions, second order diffusive processes or radiative heat transfer
- Newtonian viscous stress with no bulk viscosity contribution (i.e., :math:`\kappa S_{kk} \delta_{ij}`)
- Depending on the simulation mode, the transport coefficients :math:`\mu`, :math:`\rho\alpha_{\phi}`, and
  :math:`\rho\alpha_{\theta}` may correspond to the molecular transport coefficients, turbulent transport
  coefficients computed from an LES or PBL model, or a combination. See the sections on :ref:`DNS vs. LES modes <DNSvsLES>`
  and :ref:`PBL schemes <PBLschemes>` for more details.
