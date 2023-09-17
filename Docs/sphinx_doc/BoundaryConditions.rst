
 .. role:: cpp(code)
    :language: c++

.. _sec:BoundaryConditions:

Filling Ghost Values
--------------------------
ERF uses an operation called ``FillPatch`` to fill the ghost cells/faces for each grid of data.
The data is filled outside the valid region with a combination of three operations: interpolation
from coarser level, copy from same level, and enforcement of physical boundary conditions.

Interpolation from Coarser level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Interpolation is controlled by which interpolater we choose to use.  The default is
conservative interpolation for cell-centered quantities, and analogous for faces.
The paradigm is that fine faces on a coarse-fine boundary are filled as Dirichlet
boundary conditions from the coarser level; all faces outside the valid region are
similarly filled, while fine faces inside the valid region are not over-written.

Copy from other grids at same level (includes periodic boundaries)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is part of the ``FillPatch`` operation, but can also be applied independently,
e.g. by the call

::

    mf.FillBoundary(geom[lev].periodicity());

would fill all the ghost cells/faces of the grids in MultiFab ``mf``, including those
that occur at periodic boundaries.

In the ``FillPatch`` operation, ``FillBoundary`` always overrides any interpolated values, i.e. if
there is fine data available (except at coarse-fine boundary) we always use it.

Imposition of physical/domain boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

sets the boundary condtion type at the low x face to be an inflow with xlo.type = “Inflow”.

xlo.velocity = 1. 0. 0. sets all three componentns the inflow velocity,
xlo.density       = 1. sets the inflow density,
xlo.theta         = 300. sets the inflow potential temperature,
xlo.scalar        = 2. sets the inflow value of the advected scalar

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

The keywork "MOST" is an ERF-specific boundary type and the mapping is described below.


It is important to note that external Dirichlet boundary data should be specified
as the value on the face of the cell bounding the domain, even for cell-centered
state data.

More general boundary types are a WIP; one type that will be supported soon is the ability
to read in a time sequence of data at a domain boundary and impose this data as "ext_dir"
boundary values using ``FillPatch``.

.. _MostBoundary:

MOST Boundaries
-------------------
Monin-Obukhov similarity theory (MOST) is used to describe the atmospheric surface layer (ASL), the lowest part of the atmospheric boundary layer.  The implementation of MOST in ERF follows that in `AMR-Wind <https://github.com/Exawind/amr-wind/>`_, which is based on the surface layer profiles presented in
`P. van der Laan, et al., Wind Energy, 2017 <https://onlinelibrary.wiley.com/doi/10.1002/we.2017>`_ and
`D. Etling, "Modeling the vertical ABL structure", 1999 <https://www.worldscientific.com/doi/abs/10.1142/9789814447164_0003>`_.
MOST theory assumes that the ASL is in a steady state and horizontally homogenous, and kinematic fluxes due to turbulent transport (:math:`\overline{u^{'}w^{'}}`, :math:`\overline{v^{'}w^{'}}`, and :math:`\overline{\theta^{'}w^{'}}`) are constant with height.
:math:`\Phi_m` and :math:`\Phi_h` are the nondimensional wind shear and temperature gradient, respectively, which are assumed to follow universal similarity laws based on dimensional arguments.
With these assumptions, the MOST theory can be written as:

.. math::

  \overline{u^{'}} \overline{w^{'}} = const = -u^{2}_{\star},

  \overline{w^{'}} \overline{\theta^{'}} = const = -u_{\star}\theta_{\star},

  \Phi_{m}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{u}(z)}{\partial z},

  \Phi_{h}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{\theta}(z)}{\partial z}

where the nondimensional gradients are expressed in terms of the MOST stability parameter, :math:`\zeta = \frac{z}{L} = -\frac{\kappa z}{u_{\star}^{3}} \frac{g}{\overline{\theta}} \overline{w^{'}\theta^{'}}`, which serves as a surface layer scaling parameter.
Here, :math:`L` is the Monin-Obukhov length,
:math:`u_{\star}` is the friction velocity (defined for :math:`u` aligned with the wind direction),
:math:`\theta_{\star}` is the surface layer temperature scale,
:math:`\overline{\theta}` is the reference virtual potential temperature for the ASL,
and :math:`\kappa` is the von Karman constant (taken to be :math:`0.41`).

Integration of the MOST assumption equations give the classical MOST profiles of mean velocity and potential temperature

.. math::

  \overline{u}(z)    &= \frac{u_{\star}}{\kappa} \left[ \mathrm{ln} \left(\frac{z}{z_0}\right) - \Psi_m(\zeta)\right],

  \overline{\theta}(z) - \theta_0 &= \frac{\theta_{\star}}{\kappa} \left[ \mathrm{ln}\left(\frac{z}{z_0}\right) - \Psi_{h}(\zeta) \right]


where :math:`\theta_0` is the surface potential temperature and  :math:`z_0` is a characteristic roughness height. The integrated similarity functions,

.. math::

  \Psi_{m}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\Phi_{m}(\zeta)]\frac{d\zeta}{\zeta},

  \Psi_{h}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\Phi_{h}(\zeta)]\frac{d\zeta}{\zeta}

are calculated analytically from empirical gradient functions :math:`\Phi_m` and :math:`\Phi_h`, which are
defined piecewise for stable and unstable values of the stability parameter.

Unstable: :math:`(\zeta < 0)`

.. math::

  \Phi_{m} &= (1-\gamma_{1}\zeta)^{-\frac{1}{4}}, \quad
  \Psi_{m}    = \mathrm{ln}\left[\frac{1}{8}(1+\Phi_{m}^{-2})(1+\Phi_{m}^{-1})^{2}\right]-2\arctan(\Phi_{m}^{-1})+\frac{\pi}{2},

  \Phi_{h} &= \sigma_{\theta}(1-\gamma_{2}\zeta)^{-\frac{1}{2}}, \quad
  \Psi_{h}    = (1+\sigma_{\theta}) \mathrm{ln} \left[\frac{1}{2}(1+\Phi_{h}^{-1}) \right]+(1-\sigma_{\theta}) {\mathrm{ln}} \left[\frac{1}{2}(-1+\Phi_{h}^{-1})\right]

Stable: :math:`(\zeta > 0)`

.. math::
  \Phi_{m} &= 1+\beta \zeta, \quad \Psi_{m}=-\beta \zeta,

  \Phi_{h} &= \sigma_{\theta}+\beta \zeta, \quad \Psi_{h}=(1-\sigma_{\theta})\mathrm{ln}(\zeta)-\beta \zeta,

where the constants take the values proposed in `Dyer, Boundary Layer Meteorology, 1974
<https://link.springer.com/article/10.1007/BF00240838>`_:

.. math::
  \sigma_{\theta}=1, \quad \beta = 5, \quad \gamma_{1}=16, \quad \gamma_{2}=16

Inverting the equations above, the MOST stability parameter,

.. math::
  \zeta=\frac{z}{L} = -\kappa z \frac{g}{\bar{\theta}} \frac{\theta_{\star}}{u^{2}_{\star}}

is determined by the friction velocity

.. math::
  u_{\star} = \kappa \overline{u}/[\mathrm{ln}(z/z_0)-\Psi_{m}({z}/{L})]

and the characteristic surface layer temperature

.. math::
  \theta_{\star} = \kappa (\overline{\theta}-\theta_0)/[\mathrm{ln}(z / z_0)-\Psi_{h}(z/L)]

MOST Implementation
~~~~~~~~~~~~~~~~~~~

In ERF, when the MOST boundary condition is applied, velocity and temperature in the ghost cells are set to give stresses that are consistent with the MOST equations laid out above. The code is structured to allow either the surface temperature (:math:`\theta_0`) or surface temperature flux (:math:`\overline{w^{'}\theta^{'}}`) to be enforced. To apply the MOST boundary, the following algorithm is applied:

#. Horizontal (planar) averages :math:`\bar{u}`, :math:`\bar{v}` and :math:`\overline{\theta}` are computed at a reference height :math:`z_{ref}` assumed to be within the surface layer.

#. Initially, neutral conditions (:math:`L=\infty, \zeta=0`) are assumed and used to compute a provisional :math:`u_{\star}` using the equation given above. If :math:`\theta_0` is specified, the above equation for :math:`\theta_{\star}` is applied and then the surface flux is computed :math:`\overline{w^{'}\theta^{'}} = -u_{\star} \theta_{\star}`. If :math:`\overline{w^{'}\theta^{'}}` is specified, :math:`\theta_{\star}` is computed as :math:`-\overline{w^{'}\theta^{'}}/u_{\star}` and the previous equation is inverted to compute :math:`\theta_0`.

#. The stability parameter :math:`\zeta` is recomputed using the equation given above based on the provisional values of :math:`u_{\star}` and :math:`\theta_{\star}`.

#. The previous two steps are repeated iteratively, sequentially updating the values of :math:`u_{\star}` and :math:`\zeta`, until the change in the value of :math:`u_{\star}` on each iteration falls below a specified tolerance.

#. Once the MOST iterations have converged, and the planar average surface flux values are known, the approach from `Moeng, Journal of the Atmospheric Sciences, 1984 <https://journals.ametsoc.org/view/journals/atsc/41/13/1520-0469_1984_041_2052_alesmf_2_0_co_2.xml>`_ is applied to consistently compute local surface-normal stress/flux values (e.g., :math:`\tau_{xz} = - \rho \overline{u^{'}w^{'}}`):

   .. math::

     \left. \frac{\tau_{xz}}{\rho} \right|_0 &= u_{\star}^{2} \frac{(u - \bar{u})|\mathbf{\bar{u}}| +  \bar{u}\sqrt{u^2 + v^2} }{|\mathbf{\bar{u}}|^2},

     \left. \frac{\tau_{yz}}{\rho}  \right|_0 &= u_{\star}^{2}  \frac{(v - \bar{v})|\mathbf{\bar{u}}| +  \bar{v}\sqrt{u^2 + v^2} }{|\mathbf{\bar{u}}|^2},

     \left.  \frac{\tau_{\theta z}}{\rho} \right|_0  &= \theta_\star u_{\star} \frac{|\mathbf{\bar{u}}| ({\theta} - \overline{\theta}) +
                                                \sqrt{u^2+v^2}  (\overline{\theta} - \theta_0) }{ |\mathbf{\bar{u}}| (\overline{\theta} -\theta_0) } =
                                                u_{\star} \kappa  \frac{|\mathbf{\bar{u}}| ({\theta} - \overline{\theta})  +
                                                \sqrt{u^2+v^2} (\overline{\theta} - \theta_0) }{ |\mathbf{\bar{u}}| [  \mathrm{ln}(z_{ref} / z_0)-\Psi_{h}(z_{ref}/L)] }

   where :math:`\bar{u}`, :math:`\bar{v}` and :math:`\overline{\theta}` are the plane averaged values (at :math:`z_{ref}`) of the
   two horizontal velocity components and the potential temperature, respectively, and
   :math:`|\mathbf{\bar{u}}|` is the plane averaged magnitude of horizontal velocity (plane averaged wind speed). We note a slight variation in the denominator
   of the velocity terms from the form of the
   equations presented in Moeng to match the form implemented in AMR-Wind.

#. These local flux values are used to populate values in the ghost cells that will lead to appropiate fluxes, assuming the fluxes are computed from the turbulent transport coefficients (in the vertical direction, if applicable) :math:`K_{m,v}` and :math:`K_{\theta,v}` as follows:

   .. math::

      \tau_{xz} = K_{m,v} \frac{\partial u}{\partial z}

      \tau_{yz} = K_{m,v} \frac{\partial v}{\partial z}

      \tau_{\theta z} = K_{\theta,v} \frac{\partial \theta}{\partial z}.

   This implies that, for example, the value set for the conserved :math:`\rho\theta` variable in the :math:`-n\mathrm{th}` ghost cell is

   .. math::

      (\rho \theta)_{i,j,-n} = \rho_{i,j,-n} \left[ \frac{(\rho\theta)_{i,j,0}}{\rho_{i,j,0}} - \left. \frac{\tau_{\theta z}}{\rho} \right|_{i,j,0} \frac{\rho_{i,j,0}}{K_{\theta,v,(i,j,0)}} n \Delta z \right]

MOST Inputs
~~~~~~~~~~~~~~~~~~~
To evaluate the fluxes with MOST, the surface rougness parameter :math:`z_{0}` must be specified. This quantity may be considered a constant or may be parameterized through the friction velocity :math:`u_{\star}`. ERF supports three methods for parameterizing the surface roughness: ``constant``, ``charnock``, and ``modified_charnock``. The latter two methods parameterize :math:`z_{0} = f(u_{\star})` and are described in `Jimenez & Dudhia, American Meteorological Society, 2018 <https://doi.org/10.1175/JAMC-D-17-0137.1>`_. The rougness calculation method may be specified with

::

   erf.most.roughness_type    = STRING    #Z_0 type (constant, charnock, modified_charnock)

If the ``charnock`` method is employed, the :math:`a` constant may be specified with ``erf.most.charnock_constant`` (defaults to 0.0185). If the ``modified_charnock`` method is employed, the depth :math:`d` may be specified with ``erf.most.modified_charnock_depth`` (defaults to 30 m).

When computing an average :math:`\overline{\phi}` for the MOST boundary, where :math:`\phi` denotes a generic variable, ERF supports a variety of approaches. Specifically, ``planar averages`` and ``local region averages`` may be computed with or without ``time averaging``. With each averaging methodology, the query point :math:`z` may be determined from the following procedures: specified vertical distance :math:`z_{ref}` from the bottom surface, specified :math:`k_{index}`, or (when employing terrain-fit coordinates) specified normal vector length :math:`z_{ref}`. The available inputs to the MOST boundary and their associated data types are

::

   erf.most.average_policy    = INT    #POLICY FOR AVERAGING
   erf.most.use_normal_vector = BOOL   #USE NORMAL VECTOR W/ TERRAIN?
   erf.most.use_interpolation = BOOL   #INTERPOLATE QUERY POINT W/ TERRAIN?
   erf.most.time_average      = BOOL   #USE TIME AVERAGING?
   erf.most.z0                = FLOAT  #SURFACE ROUGHNESS
   erf.most.zref              = FLOAT  #QUERY DISTANCE (HEIGHT OR NORM LENGTH)
   erf.most.surf_temp         = FLOAT  #SPECIFIED SURFACE TEMP
   erf.most.surf_temp_flux    = FLOAT  #SPECIFIED SURFACE FLUX
   erf.most.k_arr_in          = INT    #SPECIFIED K INDEX ARRAY (MAXLEV)
   erf.most.radius            = INT    #SPECIFIED REGION RADIUS
   erf.most.time_window       = FLOAT  #WINDOW FOR TIME AVG

We now consider two concrete examples. To employ an instantaneous ``planar average`` at a specified vertical height above the bottom surface, one would specify:

::

   erf.most.average_policy    = 0
   erf.most.use_normal_vector = false
   erf.most.time_average      = false
   erf.most.z0                = 0.1
   erf.most.zref              = 1.0

By contrast, ``local region averaging`` would be employed in conjunction with ``time averaging`` for the following inputs:

::

   erf.most.average_policy    = 1
   erf.most.use_normal_vector = true
   erf.most.use_interpolation = true
   erf.most.time_average      = true
   erf.most.z0                = 0.1
   erf.most.zref              = 1.0
   erf.most.surf_temp_flux    = 0.0
   erf.most.radius            = 1
   erf.most.time_window       = 10.0

In the above case, ``use_normal_vector`` utilizes the a local surface-normal vector with length :math:`z_{ref}` to construct the positions of the query points. Each query point, and surrounding points that are within ``erf.most.radius`` from the query point, are interpolated to and averaged; for a radius of 1, 27 points are averaged. The ``time average`` is completed by way of an exponential filter function whose peak coincides with the current time step and tail extends backwards in time

.. math::

   \frac{1}{\tau} \int_{-\infty}^{0} \exp{\left(t/\tau\right)} \, f(t) \; \rm{d}t.

Due to the form of the above integral, it is advantageous to consider :math:`\tau` as a multiple of the simulation time step :math:`\Delta t`, which is specified by ``erf.most.time_window``. As ``erf.most.time_window`` is reduced to 0, the exponential filter function tends to a Dirac delta function (prior averages are irrelevant). Increasing ``erf.most.time_window`` extends the tail of the exponential and more heavily weights prior averages.

Sponge zone boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERF provides the capability to apply sponge zones at the boundaries to prevent spurious reflections that otherwise occur at the domain boundaries if standard extrapolation boundary condition is used. The sponge zone is implemented as a source term in the governing equations, which are active in a volumteric region at the boundaries that is specified by the user in the inputs file. Currently the target condition to which the sponge zones should be forced towards is to be specifed by the user in the inputs file.

.. math::

   \frac{dQ}{dt} = \mathrm{RHS} - A\xi^n(Q-Q_\mathrm{target})

where RHS are the other right-hand side terms. The parameters to be set by the user are - `A` is the sponge amplitude, `n` is the sponge strength and the `Q_\mathrm{target}` - the target solution in the sponge. `\xi` is a linear coordinate that is 0 at the beginning of the sponge and 1 at the end. An example of the sponge inputs can be found in ``Exec/RegTests/Terrain2d_Cylinder``.

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
