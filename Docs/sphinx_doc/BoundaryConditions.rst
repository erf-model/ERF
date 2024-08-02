
 .. role:: cpp(code)
    :language: c++

.. _sec:LateralBoundaryConditions:

Domain Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ideal Domain BCs
----------------------

There are two primary types of physical/domain boundary conditions: those which rely only on the
data in the valid regions, and those which rely on externally specified values.

ERF allows users to specify types of boundary condition with keywords in the inputs file.
The information for each face is preceded by
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or ``zhi``.

Currently available type of boundary conditions are
``inflow``, ``outflow``, ``slipwall``, ``noslipwall``, ``symmetry`` or ``MOST``.
(Spelling of the type matters; capitalization does not.)

For example, setting

::

    xlo.type = "Inflow"
    xhi.type = "Outflow"
    zlo.type = "SlipWall"
    zhi.type = "SlipWall"

    geometry.is_periodic = 0 1 0

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, periodic in the :math:`y`-direction,
and slip wall on the low and high :math:`y`-faces, and
Note that no keyword is needed for a periodic boundary, here only the
specification in ``geometry.is_periodic`` is needed.

Each of these types of physical boundary condition has a mapping to a mathematical boundary condition
for each type; this is summarized in the table below.

.. _sec:dirichlet:

ERF provides the ability to specify a variety of boundary conditions (BCs) in the inputs file.
We use the following options preceded by ``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, and ``zhi``:

+------------+--------------+----------------+----------------+--------------------------+---------------+
| Type       | Normal vel   | Tangential vel | Density        | Theta                    | Scalar        |
+============+==============+================+================+==========================+===============+
| inflow     | ext_dir      | ext_dir        | ext_dir        | ext_dir                  | ext_dir       |
+------------+--------------+----------------+----------------+--------------------------+---------------+
| outflow    | foextrap     | foextrap       | foextrap       | foextrap                 | foextrap      |
+------------+--------------+----------------+----------------+--------------------------+---------------+
| slipwall   | ext_dir      | foextrap       | foextrap       | ext_dir/foextrap/neumann | foextrap      |
+------------+--------------+----------------+----------------+--------------------------+---------------+
| noslipwall | ext_dir      | ext_dir        | foextrap       | ext_dir/foextrap/neumann | foextrap      |
+------------+--------------+----------------+----------------+--------------------------+---------------+
| symmetry   | reflect_odd  | reflect_even   | reflect_even   | reflect_even             | reflect_even  |
+------------+--------------+----------------+----------------+--------------------------+---------------+
| MOST       |              |                |                |                          |               |
+------------+--------------+----------------+----------------+--------------------------+---------------+

Here ``ext_dir``, ``foextrap``, and ``reflect_even`` refer to AMReX keywords.   The ``ext_dir`` type
refers to an "external Dirichlet" boundary, which means the values must be specified by the user.
The ``foextrap`` type refers to "first order extrapolation" which sets all the ghost values to the
same value in the last valid cell/face  (AMReX also has a ``hoextrap``, or "higher order extrapolation"
option, which does a linear extrapolation from the two nearest valid values). By contrast, ``neumann``
is an ERF-specific boundary type that allows a user to specify a variable gradient. Currently, the
``neumann`` BC is only supported for theta to allow for weak capping inversion
(:math:`\partial \theta / \partial z`) at the top domain.

As an example,

::

    xlo.type                =   "Inflow"
    xlo.velocity            =   1. 0.9  0.
    xlo.density             =   1.
    xlo.theta               =   300.
    xlo.scalar              =   2.

sets the boundary condition type at the low x face to be an inflow with xlo.type = “Inflow”.

xlo.velocity = 1. 0. 0. sets all three componentns the inflow velocity,
xlo.density       = 1. sets the inflow density,
xlo.theta         = 300. sets the inflow potential temperature,
xlo.scalar        = 2. sets the inflow value of the advected scalar

Additionally, one may use an input file to specify the Dirichlet velocities as a function of the
vertical coordinate z. The input file is expected to contain the following components ``{z  u  v  w}``.
For a file named ``inflow_file``, one would use the following inputs:

::

    xlo.type                =   "Inflow"
    xlo.dirichlet_file      =   "inflow_file"
    xlo.density             =   1.
    xlo.theta               =   300.
    xlo.scalar              =   2.

The "slipwall" and "noslipwall" types have options for adiabatic vs Dirichlet boundary conditions.
If a value for theta is given for a face with type "slipwall" or "noslipwall" then the boundary
condition for theta is assumed to be "ext_dir", i.e. theta is specified on the boundary.
If no value is specified then the wall is assumed to be adiabiatc, i.e. there is no temperature
flux at the boundary.  This is enforced with the "foextrap" designation.

For example

::

    zlo.type  = "NoSlipWall"
    zhi.type  = "NoSlipWall"

    zlo.theta = 301.0

would designate theta = 301 at the bottom (zlo) boundary, while
the top boundary condition would default to a zero gradient (adiabatic)
since no value is specified for ``zhi.theta`` or ``zhi.theta_grad``.
By contrast, thermal inversion may be imposed at the top boundary
by providing a specified gradient for theta

::

    zlo.type  = "NoSlipWall"
    zhi.type  = "NoSlipWall"

    zlo.theta = 301.0
    zhi.theta_grad = 1.0

We note that "noslipwall" allows for non-zero tangential velocities to be specified, as in the
Couette regression test example, in which we specify

::

    geometry.is_periodic = 1 1 0

    zlo.type = "NoSlipWall"
    zhi.type = "NoSlipWall"

    zlo.velocity    = 0.0 0.0 0.0
    zhi.velocity    = 2.0 0.0 0.0

We also note that in the case of a "slipwall" boundary condition in a simulation with non-zero
viscosity specified, the "foextrap" boundary condition enforces zero strain at the wall.

The keyword "MOST" is an ERF-specific boundary type; see :ref:`sec:MOST` for more information.

It is important to note that external Dirichlet boundary data should be specified
as the value on the face of the cell bounding the domain, even for cell-centered
state data.

Real Domain BCs
----------------------

When using real lateral boundary conditions, time-dependent observation data is read
from a file.  The observation data is utilized to directly set Dirichlet values on the
lateral domain BCs as well as nudge the solution state towards the observation data.
The user may specify (in the inputs file)
the total width of the interior Dirichlet and relaxation region with
``erf.real_width = <Int>`` (yellow + blue)
and analogously the width of the interior Dirichlet region may be specified with
``erf.real_set_width = <Int>`` (yellow).
The real BCs are only imposed for :math:`\psi = \left\{ \rho; \; \rho \theta; \; \rho q_v; \; u; \; v \right\}`.
Due to the staggering of scalars (cell center) and velocities (face center) with an Arakawa C grid,
we reduce the relaxation width of the scalars :math:`\left\{ \rho; \; \rho \theta; \; \rho q_v \right\}` by 1
to ensure the momentum updates at the last relaxation cell involve a pressure gradient that is computed with
relaxed and non-relaxed data.

.. |wrfbdy| image:: figures/wrfbdy_BCs.png
           :width: 600

.. _fig:lateral_BCs:

.. table:: Lateral boundaries

   +-----------------------------------------------------+
   |                     |wrfbdy|                        |
   +-----------------------------------------------------+
   |  Image taken from `Skamarock et al. (2021)`_        |
   +-----------------------------------------------------+

.. _`Skamarock et al. (2021)`: http://dx.doi.org/10.5065/1dfh-6p97

Within the interior Dirichlet cells, the RHS is exactly :math:`\psi^{n} - \ps^{BDY} / \Delta t`
and, as such, we directly impose this value in the yellow region.
Within the relaxation region (blue), the RHS (:math:`F`) is given by the following:

.. math::

   \begin{align}
   F &= G + R, \\
   R &= \left[ H_{1} \left( \psi^{BDY} - \psi^{\*} \right) - H_{2} \Delta^2 \left( \psi^{BDY} - \psi^{\*} \right) \right] \exp \left(-C_{01} \left(n - {\rm SpecWidth}\right)  \right), \\
   H_{1} &= \frac{1}{10 \Delta t} \frac{{\rm SpecWidth} + {\rm RelaxWidth} - n}{{\rm RelaxWidth} - 1}, \\
   H_{2} &= \frac{1}{50 \Delta t} \frac{{\rm SpecWidth} + {\rm RelaxWidth} - n}{{\rm RelaxWidth} - 1},
   \end{align}

where :math:`G` is the RHS of the Navier-Stokes equations, :math:`\psi^{*}` is the state variable at the
current RK stage, :math:`\psi^{BDY}` is temporal interpolation of the observational data, :math:`C_{01} = -\ln(0.01) / ({\rm RealWidth - SpecWidth})`
is a constant that ensure the exponential blending function obtains a value of 0.01 at the last relaxation cell,
and :math:`n` is the minimum number of grid points from a lateral boundary.

Sponge zone domain BCs
----------------------

ERF provides the capability to apply sponge zones at the boundaries to prevent spurious reflections that otherwise occur at the domain boundaries if standard extrapolation boundary condition is used. The sponge zone is implemented as a source term in the governing equations, which are active in a volumteric region at the boundaries that is specified by the user in the inputs file. Currently the target condition to which the sponge zones should be forced towards is to be specified by the user in the inputs file.

.. math::

   \frac{dQ}{dt} = \mathrm{RHS} - A\xi^n(Q-Q_\mathrm{target})

where RHS are the other right-hand side terms. The parameters to be set by the user are -- `A` is the sponge amplitude, `n` is the sponge strength and the :math:`Q_\mathrm{target}` -- the target solution in the sponge. :math:`\xi` is a linear coordinate that is 0 at the beginning of the sponge and 1 at the end. An example of the sponge inputs can be found in ``Exec/RegTests/Terrain2d_Cylinder`` and is given below. This list of inputs specifies sponge zones in the inlet and outlet of the domain in the x-direction and the outlet of the domain in the z-direction. The `start` and `end` parameters specify the starting and ending of the sponge zones. At the inlet, the sponge starts at :math:`x=0` and at the outlet the sponge ends at :math:`x=L` -- the end of the domain. The sponge amplitude `A` has to be adjust
ed in a problem-specific manner. The density and the :math:`x, y, z` velocities to be used in the sponge zones have to be specified in the inputs list.

::

          erf.sponge_strength = 10000.0
          erf.use_xlo_sponge_damping = true
          erf.xlo_sponge_end = 4.0
          erf.use_xhi_sponge_damping = true
          erf.xhi_sponge_start = 26.0
          erf.use_zhi_sponge_damping = true
          erf.zhi_sponge_start = 8.0

          erf.sponge_density = 1.2
          erf.sponge_x_velocity = 10.0
          erf.sponge_y_velocity = 0.0
          erf.sponge_z_velocity = 0.0

Another way of specifying sponge zones is by providing the sponge zone data as a text file input. This is currently implemented only for forcing :math:`x` and :math:`y` velocities in the sponge zones.
The sponge data is input as a text file with 3 columns containing :math:`z, u, v` values. An example can be found in ``Exec/SpongeTest`` and a sample inputs list for using this feature is given below. This list specifies a sponge zone in the inlet in the x-direction. The :math:`u` and :math:`v` velocity forcing in the sponge zones will be read in from the text file -- `input_sponge_file.txt`.

::

          erf.sponge_type = "input_sponge"
          erf.input_sponge_file = "input_sponge_file.txt"
          erf.sponge_strength = 1000.0
          erf.use_xlo_sponge_damping = true
          erf.xlo_sponge_end = 4.0

Inflow turbulence generation
---------------------------

ERF provides the capability to apply a perturbation zone at the inflow domain boundary to mechanically trip turbulence within the domain. The current version of the turbulence generation techniques allows for `x` and `y` (horizontal) direction perturbations. In Figure 2, as the bulk flow passes through the buoyant region, it becomes perturbed. This perturbation, combined with a short development fetch, quickly leads to the evolution of turbulence.

.. |PBinflw| image:: figures/PBIllustration.png
           :width: 600

.. _fig:pb_fig:

.. table:: Box perturbation method

   +-----------------------------------------------------+
   |                     |PBinflw|                       |
   +-----------------------------------------------------+
   |  Image taken from `DeLeon et al. (2018)`            |
   +-----------------------------------------------------+

Two different types of perturbation methods are currently available, ``source`` and ``direct``. Both methods uses the formulation introduced by `DeLeon et al. (2018)`_. The ``source`` option applies the perturbation amplitude range, :math:`\pm \Phi_{PB}`, to each cell within the perturbation box as a source term. Conversely, the ``direct`` option applies the calculated temperature difference directly onto the :math:`\rho \theta` field. With should note that while both methods effectively generates turbulence downstream, the ladder approach is more unstable and requires more fine tuning.

..
  _`DeLeon et al. (2018)`: https://doi.org/10.2514/1.J057245

The perturbation update interval of the individual perturbation box is determined by,

.. math::

   \frac{t_p \langle U(z) \rangle_{PB}}{D_{PB}} = 1,

The change in the perturbation amplitude is defined as,

.. math::

   {Ri}_{PB} = \frac{g \beta \Delta \overline{\phi} H_{PB}}{{\langle U(z) \rangle}^2_{PB}}.

The current implementation supports lateral boundary perturbations, specified by ``erf.perturbation_direction``, where the six integer inputs represent the west, south, bottom, east, north, and top domain face perturbations, respectively. Note that while the top and bottom options are included, triggering either options will abort the program.

In addition to the direction of perturbation, the flow perturbation method requires the dimensions of an individual box, specified through ``erf.perturbation_box_dim``, with three integer inputs representing :math:`{Nx}_{pb}`, :math:`{Ny}_{pb}`, and :math:`{Nz}_{pb}`, respectively. Following the guidance of `Ma and Senocak (2023)`_, the general rule of thumb is to use :math:`H_{PB} = 1/8 \delta` as the height of the perturbation box, where :math:`\delta` is the boundary layer height. The length of the box in the x-direction should be :math:`L_{PB} = 2H_{PB}`. Depending on the direction of the bulk flow, the width of the box in the y-direction should be defined as :math:`W_{PB} = L_{PB} \tan{(\theta_{inflow})}`. The current implementation only accepts integer entries. Therefore, considering the domain size and mesh resolution, the dimensions of a single box can be determined.

.. _`Ma and Senocak (2023)`: https://link.springer.com/article/10.1007/s10546-023-00786-1

Specification of the number of layers and the offset into the domain of the perturbation boxes can be made through ``erf.perturbation_layers`` and ``erf.perturbation_offset``, respectively. Below is an example of the required input tags to setup a simulation with an inflow perturbation simulation.

::

          erf.inlet_perturbation_type = "source"

          erf.perturbation_direction = 1 0 0 0 0 0
          erf.perturbation_box_dims = 8 8 4
          erf.perturbation_layers = 3
          erf.perturbation_offset = 5

          erf.perturbation_nondimensional = 0.042
          erf.perturbation_T_infinity = 300.0
          #erf.perturbation_T_intensity = 0.1

The ``erf.perturbation_T_intensity`` tag can be turned on or off by providing a value or commenting it out. When a value is provided (recommended 0.1-0.25 max), a pseudo-gravity value is used (solved from the Richardson formulation) to normalize the scales of the problem, and is represented as,

.. math::
   g = \frac{{Ri}_{PB} {\langle U(z) \rangle}^2_{PB}}{\beta \Delta \phi H_{pb}}.

Using this pseudo-gravity value effectively negates the Richardson number formulation, and the temperature gradient becomes,

.. math::
   \Delta \phi = T_{i} T_{\infty},

where :math:`T_{i}` is the temperature intensity, and :math:`T_{\infty}` is the background temperature.

While this generates quick turbulence, it should be used as a sanity check rather than a runtime strategy for turbulence generation, therefore it is not recommended. Additionally, a net-zero energy enforcement is applied over the perturbation boxes to ensure that the synthetic method does not introduce excess energy into the system at each iteration. Below, we provide a detailed description of the two different types of perturbation methods currently existing within ERF.

Examples are provided within ``Exec/DevTests/ABL_perturbation_inflow/`` to set up a turbulent open channel flow using inflow/outflow boundary conditions with the aforementioned turbulent inflow generation technique to trigger turbulence downstream.

Source type forcing
-------------------
By ignoring the advection and diffusion effects in the transport equation, the amplitude of the peturbation can be made through a proportionality ratio between the update time and change in the box temperature,

.. math::

   \Phi_{PB} \propto \frac{\Delta \overline{\phi}}{t_p}

and the perturbation amplitude is then computed through,

.. math::

   \Phi_{PB} = \frac{Ri_{PB} {\langle U(z) \rangle}^3_{PB}}{g \beta D_{PB} H_{PB}}.

``source`` type forcing can adopt the box perturbation method by having the following inputs list.

::

          erf.inlet_perturbation_type = "source"

          erf.perturbation_direction = 1 0 0 0 0 0
          erf.perturbation_box_dims = 8 8 4
          erf.perturbation_layers = 3
          erf.perturbation_offset = 5

          erf.perturbation_nondimensional = 0.042 # Ri
          erf.perturbation_T_infinity = 300.0

The box perturbation method (BPM) perturbs the temperature field :math:`\rho \theta` in a volume (box) format. Each box computes a perturbation update time and amplitude, then independently update at its respective update interval during runtime. A single perturbation amplitude is seen by the computational cells that falls within this box. A pseudo-random perturbations (ie. white noise) is applied over the range :math:`[-\phi, +\phi]` then introduced to the :math:`\rho \theta` field via source term (:math:`F_{\rho \theta}`). As temperature is transported and through the action of the subgrid-scale (SGS) filter for eddy viscosity, white noise temperature perturbations become colored noise in the velocity field.

Using the source term to perturb the momentum field through buoyancy coupling requires more time compared to the "direct" perturbation method. Turbulence onset can be triggered by adjusting the size of the perturbation box, as the amplitude is heavily influenced by having the two dimensional scales of the perturbation box in the denominator.

Direct type forcing
-------------------
The ``direct`` method can also be used to effectively trip turbulence into the domain. Minute differences exists between the ``direct`` and ``source`` method, with the primary one being the perturbation amplitude is directly solved from the Richardson number relationship as,

.. math::

   \Phi_{PB} = \frac{Ri_{PB} {\langle U(z) \rangle}^2_{PB}}{g \beta H_{PB}}.

To activate the ``direct`` type forcing, set the following inputs list.
::

          erf.inlet_perturbation_type = "direct"

          erf.perturbation_direction = 1 0 0 0 0 0
          erf.perturbation_box_dims = 16 16 8
          erf.perturbation_layers = 3
          erf.perturbation_offset = 5

          erf.perturbation_nondimensional = 0.042 # Ri
          erf.perturbation_T_infinity = 300.0

We want to note that the ``direct`` perturbation method is sensitive to the temperature amplitude computed through the equation above, and is subject to crash the simulation when this amplitude is too large.
