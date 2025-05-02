
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

Here we present the compressible equations; see :ref:`sec:AnelasticEquations` for the anelastic equation set.

Compressible Equations
------------------------

The governing equations for compressible flow with precipitating and non-precipitating moisture components are

.. math::
   \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) - \frac{1}{1 + q_v + q_c} ( \nabla p^{\prime}  - \delta_{i,3}\mathbf{B} ) - \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_d \theta_d)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u}        \theta_d) + \nabla \cdot ( \rho_d \alpha_{\theta}\ \nabla \theta_d) + F_{\theta} + H_{n} + H_{p},

   \frac{\partial (\rho_d \boldsymbol{\phi})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \boldsymbol{\phi}) + \nabla \cdot ( \rho_d \alpha_{\phi}\ \nabla \boldsymbol{\phi}) + \mathbf{F}_{\phi},

   \frac{\partial (\rho_d \mathbf{q_{n}})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{n}}) + \nabla \cdot (\rho_d \alpha_{q} \nabla \mathbf{q_{n}}) + \mathbf{F_{n}} + \mathbf{G_{p}},

   \frac{\partial (\rho_d \mathbf{q_{p}})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{p}}) + \partial_{z} \left( \rho_d \mathbf{w_{t}} \mathbf{q_{p}} \right) + \mathbf{F_{p}}.

Here :math:`\rho, T, \theta_{d}`, and :math:`p` are the density, temperature, dry potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`\phi` is an advected scalar, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`P_{00}` is a reference value for pressure.

The pressure perturbation is computed as

.. math::
  p^\prime = P_{00} \left( \frac{R_d \rho_d \theta_m}{P_{00}} \right)^\gamma - p_{0}

where :math:`\gamma = c_{p} / (c_{p} - R_{d})` and

.. math::
  \theta_m = \theta_d (1 + \frac{R_v}{R_d} q_v)

is the moist potential temperature.  We note that this is the only place :math:`\theta_m` is used; we evolve :math:`\theta_d` above.

and

- :math:`\boldsymbol{\tau}` is the viscous stress tensor,

  .. math::
     \tau_{ij} = -2\mu \sigma_{ij},

with :math:`\sigma_{ij} = S_{ij} -D_{ij}` being the deviatoric part of the strain rate, and

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right), \hspace{24pt}
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij},

- :math:`\mathbf{F}_{u}` and :math:`F_{\theta_d}` are the forcing terms described in :ref:`Forcings`,
- :math:`\mathbf{B} = -(\rho - \rho_{0})\mathbf{g}` is the buoyancy term described in :ref:`sec:Buoyancy <Buoyancy>`,
- :math:`\mathbf{g} = (0,0,-g)` is the gravity vector,
- the dry potential temperature :math:`\theta_d` is defined from temperature :math:`T`, pressure :math:`p`, and reference pressure :math:`P_{00} = 10^{5}` Pa as

.. math::

  \theta_d = T \left( \frac{P_{00}}{p} \right)^{R_d / c_p}.

- pressure and density are defined as perturbations from a hydrostatically stratified background state, i.e.
.. math::

  p = p_{0}(z) + p^\prime  \hspace{24pt} \rho = \rho_{0}(z) + \rho^\prime

with

.. math::

  \frac{d p_{0}}{d z} = - \rho_{0} g

With this model, in addition to dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`c_{vd}`, :math:`c_{vv}`, :math:`c_{pd}`, :math:`C_{pv}`,
we include
non-precipitating condensates :math:`\rho_c + \rho_i`,
and precipitating condensates :math:`\rho_p = \rho_{rain} + \rho_{snow} + \rho_{graupel}`.
Here
:math:`\rho_c` is the density of cloud water and
:math:`\rho_i` is the density of cloud ice, and
we define the sum of all non-precipitating moist quantities to be :math:`\rho_T = \rho_v + \rho_c + \rho_i`.
All condensates  are treated as incompressible; cloud water and ice
have constant heat capacities :math:`C_p` and :math:`C_i`, respectively.

Neglecting the volume occupied by all water not in vapor form, we have

.. math::
  p = p_d + p_v = \rho_d R_d T + \rho_v R_v T

where :math:`p_d` and :math:`p_v` are the partial pressures of dry air and water vapor, respectively,
and :math:`R_d` and :math:`R_v` are the gas constants for dry air and water vapor, respectively.

We define the mixing ratio of each moist component, :math:`q_s`, as the mass density of species :math:`s`
relative to the density of dry air, i.e. :math:`q_s = \frac{\rho_s}{\rho_d}`.

We define the total potential temperature

.. math::
  \theta = \frac{\sum_s \rho_s \theta_s}{\sum_s \rho_s} \approx (\theta_d + q_v \theta_v + q_i \theta_i + q_c \theta_c).

and write the EOS as

.. math::
   T = \theta (\frac{p}{p_0})^\frac{R^\star}{C_p^\star}

or

.. math::
   p = P_{00} (\frac{\Pi}{c_p^\star})^{\frac{c_p^\star}{R^\star}}

where :math:`P_{00}` is the reference pressure. and

.. math::
  \Pi = C_p^\star (\frac{p}{\alpha p_0})^\frac{R^\star}{C_p^\star}

with :math:`\alpha = \frac{R^\star}{p}(\frac{p}{p_0})^\frac{R^\star}{c_p^\star} \theta`

here, :math:`R^\star =  R_{d} + q_v R_{v} + q_i R_{i} + q_p R_{p}`, and :math:`C_p^\star = C_{pd} + q_v C_{pv} + q_i C_{pi} + q_p C_{pp}`.

:math:`R_d`, :math:`R_v`, :math:`R_i`, and :math:`R_p` are the gas constants for dry air, water vapor, cloud ice, precipitating condensates, respectively. :math:`C_{pd}`, :math:`C_{pv}`, :math:`C_{pi}`, and :math:`C_{pp}` are the specific heats for dry air,
water vapor, cloud ice, and precipitating condensates, respectively.

The non-precipitating water mixing ratio vector :math:`\mathbf{q_{n}} = \left[ q_v \;\; q_c \;\; q_i \right]` includes water vapor, :math:`q_v`, cloud water, :math:`q_c`, and cloud ice, :math:`q_i`, although some models may not include cloud ice; similarly, the precipitating water mixing ratio vector :math:`\mathbf{q_{p}} = \left[ q_r \;\; q_s \;\; q_g \right]` involves rain, :math:`q_r`, snow, :math:`q_s`, and graupel, :math:`q_g`, though some models may not include these terms. The source terms for moisture variables, :math:`\mathbf{F_{p}}`, :math:`\mathbf{F_{n}}`, :math:`\mathbf{G_{p}}`, and their corresponding impact on potential temperature, :math:`H_{n}` and :math:`H_{p}`, and the terminal velocity, :math:`\mathbf{w_{t}}` are specific to the employed model. For the Kessler microphysics scheme, these terms are detailed in :ref:`sec:Kessler Microphysics model <Microphysics>`.

Assumptions
------------------------

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- No chemical reactions, second order diffusive processes or radiative heat transfer
- Newtonian viscous stress with no bulk viscosity contribution (i.e., :math:`\kappa S_{kk} \delta_{ij}`)
- Depending on the simulation mode, the transport coefficients :math:`\mu`, :math:`\rho\alpha_{\phi}`, and
  :math:`\rho\alpha_{\theta}` may correspond to the molecular transport coefficients, turbulent transport
  coefficients computed from an LES or PBL model, or a combination. See the sections on :ref:`DNS vs. LES modes <DNSvsLES>`
  and :ref:`PBL schemes <PBLschemes>` for more details.

Diagnostic Relationships
------------------------

In order to close the above prognostic equations, a relationship between the pressure and the other state variables
must be specified. This is obtained by re-expressing the ideal gas equation of state in terms of :math:`\theta_{d}`:

.. math::
   p = \left( \frac{\rho R_d \theta_{d}}{P_{00}^{R_d / c_p}} \right)^\gamma = P_{00} \left( \frac{\rho R_d \theta_{d}}{P_{00}} \right)^\gamma

