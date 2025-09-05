
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _GoverningEquations:

Governing Equations
=============================

ERF can be run in two different modes: in the first, ERF solves the fully compressible fluid equations,
in the second, ERF solves a modified set of equations which approximates the density field with the
hydrostatic density and imposes the anelastic constraint on the velocity field.

In compressible mode, ERF solves partial differential equations expressing conservation of mass, momentum,
potential temperature, and scalars (such as moisture variables) subject to an equation of state.

In anelastic mode, ERF solves partial differential equations expressing conservation of momentum,
potential temperature, and scalars (such as moisture variables), as well the anelastic constraint
on the velocity.

Below :math:`\rho, T, \theta_{d}`, and :math:`p` are the density, temperature, dry potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`\phi` is an advected scalar, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

Compressible Equations
------------------------

The first three equations governing fully compressible flow are

.. math::
   \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
                                                        \frac{1}{1 + q_t} ( \nabla p^{\prime}  - \delta_{i,3}\mathbf{B} ) -
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
                                                        \frac{1}{1 + q_t} ( \nabla p^\prime + \delta_{i,3}\mathbf{B} ) -
                                                        \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_0 \theta_d)}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} \theta_d) +
                                                      \nabla \cdot (\rho_0 \alpha_{\theta}\ \nabla \theta_d) +
                                                      F_{\theta} + H_{n} + H_{p},

supplemented with the constraint

.. math::
  \nabla \cdot (\rho_0 \mathbf{u}) = 0

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
To account for isotropic lateral grid stretching as represented by ``map factors" :math:`m_x = m_y = m` as in WRF, we augment the coordinate transform above with stretching in the lateral directions only.

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

- Pressure and density perturbations are defined with respect to a hydrostatically stratified background state, i.e.
.. math::
  p = p_{0}(z) + p^\prime  \hspace{24pt} \rho = \rho_{0}(z) + \rho^\prime

with

.. math::
  \frac{d p_{0}}{d z} = - \rho_{0} g

Equation of state (compressible only)
--------------------------------------

In the fully compressible formulation, the total pressure is computed as

.. math::
  p = P_{00} \left( \frac{R_d \rho_d \theta_m}{P_{00}} \right)^\gamma

where :math:`\gamma = c_{p} / (c_{p} - R_{d})` and

.. math::
  \theta_m = \theta_d (1 + \frac{R_v}{R_d} q_v)

is the moist potential temperature. This is the only place :math:`\theta_m` is used; we evolve :math:`\theta_d` above. In the above, :math:`R_d`, :math:`c_p`, :math:`P_{00} = 1\times10^{5}` are the gas constant, specific heat capacity for dry air, and reference pressure, respectively.

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
- :math:`\mathbf{B} = -(\rho - \rho_{0})\mathbf{g}` is the buoyancy term described in :ref:`Buoyancy <Buoyancy>`,
- :math:`\mathbf{g} = (0,0,-g)` is the gravity vector,
- The dry potential temperature :math:`\theta_d` is defined from temperature :math:`T`, pressure :math:`p`, and reference pressure :math:`P_{00} = 10^{5}` Pa as

.. math::

  \theta_d = T \left( \frac{P_{00}}{p} \right)^{R_d / c_p}.

(In the anelastic case, :math:`p` is replaced by :math:`p_0` in the relationship between :math:`\theta_d` and :math:`T`.)


Assumptions
------------------------

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- No chemical reactions, second order diffusive processes or radiative heat transfer
- Newtonian viscous stress with no bulk viscosity contribution (i.e., :math:`\kappa S_{kk} \delta_{ij}`)
- Depending on the simulation mode, the transport coefficients :math:`\mu`, :math:`\rho\alpha_{\phi}`, and
  :math:`\rho\alpha_{\theta}` may correspond to the molecular transport coefficients, turbulent transport
  coefficients computed from an LES or PBL model, or a combination. See the sections on :ref:`DNS vs. LES modes <DNSvsLES>`
  and :ref:`PBL schemes <PBLschemes>` for more details.

