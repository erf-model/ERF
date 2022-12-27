
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:

Prognostic Equations (Dry)
=============================

The following partial differential equations governing dry compressible flow
are solved in ERF for mass, momentum, potential temperature, and scalars:

.. math::
  \frac{\partial \rho}{\partial t} &= - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} \mathbf{u}) - \nabla p^\prime
        + \delta_{i,3}\mathbf{B} + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho \theta)}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} \theta) + \nabla \cdot ( \rho \alpha_{T}\ \nabla \theta) + F_{\rho \theta},

  \frac{\partial (\rho C)}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} C) + \nabla \cdot (\rho \alpha_{C}\ \nabla C)

where

- :math:`\tau` is the viscous stress tensor,

  .. math::
     \tau_{ij} = 2\mu \sigma_{ij},

with :math:`\sigma_{ij} = S_{ij} -D_{ij}` being the deviatoric part of the strain rate, and

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right), \hspace{24pt}
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij},

- :math:`\mathbf{F}` and :math:`F_{\rho \theta}` are the forcing terms described in :ref:`Forcings`,
- :math:`\mathbf{g} = (0,0,-g)` is the gravity vector,
- the potential temperature :math:`\theta` is defined from temperature :math:`T` and pressure :math:`p` as

.. math::

  \theta = T \left( \frac{p_0}{p} \right)^{R_d / c_p}.

- pressure and density are defined as perturbations from a hydrostatically stratified background state, i.e.
.. math::

  p = \overline{p}(z) + p^\prime  \hspace{24pt} \rho = \overline{\rho}(z) + \rho^\prime

with

.. math::

  \frac{d \overline{p}}{d z} = - \overline{\rho} g

Assumptions
------------------------

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- No chemical reactions, second order diffusive processes or radiative heat transfer
- Newtonian viscous stress with no bulk viscosity contribution (i.e., :math:`\kappa S_{kk} \delta_{ij}`)
- Depending on the simulation mode, the transport coefficients :math:`\mu`, :math:`\rho\alpha_C`, and
  :math:`\rho\alpha_T` may correspond to the molecular transport coefficients, turbulent transport
  coefficients computed from an LES or PBL model, or a combination. See the sections on :ref:`DNS vs. LES modes <DNSvsLES>`
  and :ref:`PBL schemes <PBLschemes>` for more details.

Diagnostic Relationships
------------------------

In order to close the above prognostic equations, a relationship between the pressure and the other state variables
must be specified. This is obtained by re-expressing the ideal gas equation of state in terms of :math:`\theta`:

.. math::
   p = \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma

Nomenclature
------------
Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`C` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.


Prognostic Equations (Moist)
===============================
Thermodynamics and the Specific Equation of States
--------------------------------------------------
We consider a mixture of dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`C_{vd}`, :math:`C_{vv}`, :math:`C_{pd}`, :math:`C_{pv}`,
non-precipitating condensates :math:`\rho_c + \rho_i`,
and precipitating condensates :math:`\rho_{rain} + \rho_{snow} + \rho_{graupel}`.
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

In ERF, we select the dry air with density :math:`\rho_d` as the dominant component, and treat the others as sparse components
:math:`\rho_s` with :math:`s = 1, ...., N`. We define the mass mixing ratio, :math:`q_s`, as the mass density of species :math:`s`
relative to the total density, i.e. :math:`q_s = \frac{\rho_s}{\rho}`.  We note that

.. math::
  \sum_s \rho_s = \rho

  \sum_s q_s = 1

where :math:`\rho` is the moist air density.

define the total potential temperature

.. math::
  \theta = \frac{\sum_s \rho_s \theta_s}{\sum_s \rho_s} \approx (q_d \theta_d + q_v \theta_v + q_i \theta_i + q_c \theta_c).

the EOS equation can be written as,

.. math::
   T = \theta (\frac{p}{p_0})^\frac{R^\star}{C_p^\star}

.. math::
   p = p_0 (\frac{\Pi}{C_p^\star})^{\frac{C_p^\star}{R^\star}}

where :math:`p_0` is the reference pressure. and

.. math::
  \Pi = C_p^\star (\frac{p}{\alpha p_0})^\frac{R^\star}{C_p^\star}

with :math:`\alpha = \frac{R^\star}{p}(\frac{p}{p_0})^\frac{R^\star}{c_p^\star} \theta`

here, :math:`R^\star =  q_d R_{d} + q_v R_{v} + q_i R_{i} + q_p R_{p}`, and :math:`C_p^\star = q_d C_{pd} + q_v C_{pv} + q_i C_{pi} + q_p C_{pp}`. the :math:`R_d`,
:math:`R_v`, :math:`R_i`, and :math:`R_p` are the gas constants for dry air, water vapor, cloud ice, precipitating condensates, respectively. :math:`C_{pd}`, :math:`C_{pv}`, :math:`C_{pi}`, and :math:`C_{pp}` are the specific heat for dry air, water vapor, cloud ice, and precipitating condensates, respectively.

Governing Equations for Moist Atmospheric Flow
-------------------------------------------------------
We assume that all species have same average speed,
Then the governing equations become

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} + \mathbf{F}_\rho)

  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u} + \mathbf{F}_u) -
          \frac{1}{1 + q_T + q_p}  \nabla p^\prime + \nabla \cdot \tau + \mathbf{F} + \delta_{i,3}\mathbf{B}

  \frac{\partial (\rho_d \theta)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta + F_{\theta}) + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \theta) + F_Q

  \frac{\partial (\rho_d C)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} C + \mathbf{F}_C) + \nabla \cdot (\rho_d \alpha_{C}\ \nabla C)

  \frac{\partial (\rho_d q_T)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_T +F_{q_{T}}) - Q

  \frac{\partial (\rho_d q_p)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_p + F_{q_{p}}) + Q

In this set of equations, the subgrid turbulent parameterization effects are included with fluxes
:math:`F_\rho`, :math:`F_u`, :math:`F_C`, :math:`F_{\theta}`, :math:`F_{q_{T}}`, :math:`F_{q_{r}}`.
:math:`\mathbf{F}` stands for the external force, and :math:`Q` and :math:`F_Q` represent the mass and energy transformation
of water vapor to/from water through condensation/evaporation, which is determined by the microphysics parameterization processes.
:math:`\mathbf{B}` is the buoyancy force.

