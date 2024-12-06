
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _DryEquations:

Compressible Equations (Dry)
=============================

ERF can be run in two different modes: in the first, ERF solves the fully compressible fluid equations,
in the second, ERF solves a modified set of equations which approximates the density field with the
hydrostatic density and imposes the anelastic constraint on the velocity field.

In compressible mode, in the absence of moisture, ERF solves the following partial differential equations
expressing conservation of mass :math:`(\rho)`, momentum :math:`(\rho \mathbf{u})`, potential temperature :math:`(\rho \theta_{d})`, and scalars :math:`(\rho \mathbf{\phi})`:

.. math::
   \frac{\partial \rho_d}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) - \frac{1}{1 + q_v + q_c} ( \nabla p^\prime  - \delta_{i,3}\mathbf{B} ) - \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_d \theta_d)}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \theta_d) + \nabla \cdot ( \rho_d \alpha_{\theta}\ \nabla \theta_d) + F_{\theta},

   \frac{\partial (\rho_d \boldsymbol{\phi})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \boldsymbol{\phi}) + \nabla \cdot ( \rho_d \alpha_{\phi}\ \nabla \boldsymbol{\phi}) + \mathbf{F}_{\phi}.

where

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

We note that there is an alternative option under development in ERF that solves the governing
equations with an anelastic constraint rather than the fully compressible equations.  The equation set is described below.

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

Nomenclature
------------
Here :math:`\rho, T, \theta_{d}`, and :math:`p` are the density, temperature, dry potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`\phi` is an advected scalar, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`P_{00}` is a reference value for pressure.
