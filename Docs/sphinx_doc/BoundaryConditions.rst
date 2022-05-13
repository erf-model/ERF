
 .. role:: cpp(code)
    :language: c++

.. _sec:domainBCs:

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

Dirichlet Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERF provides the ability to specify constant Dirichlet BCs in the inputs file. We use the following options
preceded by
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, and ``zhi``:

+------------+--------------+----------------+----------------+------------------+---------------+
| Type       | Normal vel   | Tangential vel | Density        | Theta            | Scalar        |
+============+==============+================+================+==================+===============+
| inflow     | ext_dir      | ext_dir        | ext_dir        | ext_dir          | ext_dir       |
+------------+--------------+----------------+----------------+------------------+---------------+
| outflow    | foextrap     | foextrap       | foextrap       | foextrap         | foextrap      |
+------------+--------------+----------------+----------------+------------------+---------------+
| slipwall   | ext_dir      | foextrap       | foextrap       | ext_dir/foextrap | foextrap      |
+------------+--------------+----------------+----------------+------------------+---------------+
| noslipwall | ext_dir      | ext_dir        | foextrap       | ext_dir/foextrap | foextrap      |
+------------+--------------+----------------+----------------+------------------+---------------+
| symmetry   | reflect_odd  | reflect_even   | reflect_even   | reflect_even     | reflect_even  |
+------------+--------------+----------------+----------------+------------------+---------------+
| MOST       |              |                |                |                  |               |
+------------+--------------+----------------+----------------+------------------+---------------+

Here ``ext_dir``, ``foextrap``, and ``reflect_even`` refer to AMReX keywords.   The ``ext_dir`` type
refers to an "external Dirichlet" boundary, which means the values must be specified by the user.
The ``foextrap`` type refers to "first order extrapolation" which sets all the ghost values to the
same value in the last valid cell/face.  (AMReX also has a ``hoextrap``, or "higher order extrapolation"
option, which does a linear extrapolation from the two nearest valid values.)

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
If not value is specified then the wall is assumed to be adiabiatc, i.e. there is no temperature
flux at the boundary.  This is enforced with the "foextrap" designation.

For example

::

    zlo.type  = "NoSlipWall"
    zhi.type  = "NoSlipWall"

    zlo.theta = 301.0

would designate theta = 301. at the bottom (zlo) boundary, while
the top boundary condition would be a zero gradient, aka adiabatic
because no value is specified for ``zhi.theta``

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
Monin-Obukhov similarity theory (MOST) is used to describe the atmospheric surface layer (ASL), the lowest part of the atmospheric boundary layer. MOST theory assumes that the ASL is in a steady state and horizontally homogenous, and kinematic fluxes due to turbulent transport (:math:`\overline{u^{'}w^{'}}`, :math:`\overline{v^{'}w^{'}}`, and :math:`\overline{\theta^{'}w^{'}}`) are constant with height.
:math:`\phi_m` and :math:`\phi_h` are the nondimensional wind shear and temperature gradient, respectively.
Based on these assumptions, the MOST theory can be written as:

.. math::

  \overline{u^{'}} \overline{w^{'}} = const = -u^{2}_{\star},

  \overline{w^{'}} \overline{\theta^{'}} = const = -u_{\star}\theta_{\star},

  \phi_{m}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{U}(z)}{\partial z},

  \phi_{h}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{\theta}(z)}{\partial z}

where the nondimensional gradients are expressed in terms of the MOST stability parameter, :math:`\zeta = \frac{z}{L} = -\frac{\kappa z}{u_{\star}^{3}} \frac{g}{\overline{\theta}} \overline{w^{'}\theta^{'}}`, which serves as a surface layer scaling parameter.
Here, :math:`L` is the Monin-Okukhov length,
:math:`u_{\star}` is the friction velocity (defined for :math:`u` aligned with the wind direction),
:math:`\theta_{\star}` is the surface layer temperature scale,
and :math:`\overline{\theta}` is the reference virtual potential temperature for the ASL.

Integration of the MOST assumption equations give the classical MOST profiles of mean velocity and potential temperature

.. math::

  \overline{U}(z)    &= \frac{u_{\star}}{\kappa} \left[ \mathrm{ln} \left(\frac{z}{z_0}\right) - \psi_m(\zeta)\right],

  \overline{\theta}(z) - \theta_0 &= \frac{\theta_{\star}}{\kappa} \left[ \mathrm{ln}\left(\frac{z}{z_0}\right) - \psi_{h}(\zeta) \right]


where :math:`\theta_0` is the surface potential temperature. The integrated similarity functions,

.. math::

  \psi_{m}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\phi_{m}(\zeta)]\frac{d\zeta}{\zeta},

  \psi_{h}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\phi_{h}(\zeta)]\frac{d\zeta}{\zeta}

are calculated analytically for stable and unstable stratifications, given empirical gradient functions :math:`\phi_m` and :math:`\phi_h`.

Unstable: :math:`(\zeta < 0)`

.. math::

  \phi_{m} &= (1-\gamma_{1}\eta)^{-\frac{1}{4}}, \quad
  \psi_{m}    = \mathrm{ln}[\frac{1}{g}(1+\psi_{m}^{2})(1+\psi_{m}^{-1})^{2}]-2\arctan(\Theta_{m}^{-1})+\frac{\pi}{2},

  \phi_{h} &= \sigma_{\theta}(1-\gamma_{2}\zeta)^{-\frac{1}{2}}, \quad
  \psi_{h}    = (1+\sigma_{\theta}) \mathrm{ln}[\frac{1}{2}(1+\Theta_{h}^{-1}]+(1-\sigma_{\theta}) {\mathrm{ln}} [\frac{1}{2}(-1+\Theta_{h}^{-1})]

Stable: :math:`(\zeta > 0)`

.. math::
  \phi_{m} &= 1+\beta \zeta, \quad \psi_{m}=-\beta \zeta,

  \phi_{h} &= \sigma_{\theta}+\beta \zeta, \quad \psi_{h}=(1-\sigma_{\theta})\mathrm{ln}(\zeta)-\beta \zeta

and the constants are defined as:

.. math::
  \sigma_{\theta}=1, \quad \beta = 5, \quad \gamma_{1}=16, \quad \gamma_{2}=16

Inverting the equations above, the MOST stability parameter,

.. math::
  \zeta=\frac{z}{L} = -\kappa z \frac{g}{\theta_{0}} \frac{\theta_{\star}}{u^{2}_{\star}}

is determined by the friction velocity

.. math::
  u_{\star} = \kappa \overline{U}/[\mathrm{ln}(z/z_0)-\psi_{m}(\frac{z}{L})]

and the surface temperature

.. math::
  \theta_{\star} = \kappa (\overline{\theta}-\theta_0)/[\mathrm{ln}(z / z_0)-\psi_{h}(z/L)]

Assuming that :math:`\theta_{\star}, u_{\star}, q_{\star}` are constant with height, the wind speed, temperature and moisture at surface can be derived as:

.. math::

  \mu\frac{\partial u}{\partial z} &= -\rho u_{\star}^{2} \frac{ (u-\bar{u}) \bar{|\mathbf{U}|} +
                             \sqrt{u^2+v^2} \; \bar{u}  }{\bar{|\mathbf{U}|}^{2}},

 \mu\frac{\partial v}{\partial z} &= -\rho u_{\star}^{2} \frac{ (v-\bar{v}) \bar{|\mathbf{U}|} +
                             \sqrt{u^2+v^2} \; \bar{v}  }{\bar{|\mathbf{U}|}^{2}},

  \mu\frac{\partial \theta }{\partial z} &= \rho \frac{\bar{|\mathbf{U}|} (\theta - \overline{\theta}) +
                                             \sqrt{u^2+v^2}  (\overline{\theta}-\theta_{\star}) }{\Theta_{h} \bar{|\mathbf{U}|}}

where :math:`\bar{u}`, :math:`\bar{v}` and :math:`\overline{\theta}` are the planar-averaged values of the
two horizontal velocity components and the potential temperature, respectively, and
:math:`\bar{|\mathbf{U}|}` is the planar-averaged magnitude of horizontal velocity. :math:`\mu` is the eddy visocisity.

