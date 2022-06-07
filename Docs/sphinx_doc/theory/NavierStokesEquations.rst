
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

  \frac{\partial (\rho \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} \mathbf{u}) - \nabla p^\prime +\rho^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F},

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
- Though not required in the above form of the equations, we further assume either
    - constant transport coefficients :math:`\mu`, :math:`\rho\alpha_C`, and :math:`\rho\alpha_T`, or
    - constant :math:`\nu = \mu / \rho`, :math:`\alpha_T` and constant :math:`\alpha_C`,
      each of which is then multiplied by :math:`\rho` in the construction of the diffusive terms.
  This is a good approximation for flows of interest because all are independent of density (or pressure),
  and only weakly dependent on temperature (:math:`T^{1/2}` scaling based on simple kinetic theory).

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

When solving for moist atmospheric flow, we evolve two additional variables: :math:`q_v` and :math:`q_c`,
the mixing ratios of water vapor and cloud water, respectively.

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u})

  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
  \frac{\rho_d}{\rho_m} \nabla p^\prime +\rho_d^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}

  \frac{\partial (\rho_d \theta_m)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta_m) + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \theta) + F_Q

  \frac{\partial (\rho_d C)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} C) + \nabla \cdot (\rho_d \alpha_{C}\ \nabla C)

  \frac{\partial (\rho_d q_v)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_v) - Q

  \frac{\partial (\rho_d q_c)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_c) + Q

where :math:`\rho_d` is the density of the dry air only, :math:`\rho_m = \rho_d (1 + q_v + q_c)` is the total mass density,
and :math:`Q` represents the transformation of water vapor to cloud water through condensation.

.. math::

  \theta_m = (1 + \frac{R_v}{R_d} q_v) \theta

is the modified potential temperature.  Here :math:`R_v` is the gas constant for water vapor and :math:`R_v / R_d \approx 1.61.`
The equation of state is now

.. math::
   p = \left( \frac{\rho_d R_d \theta_m}{p_0^{R_d / c_p}} \right)^\gamma
   
   
Single Moment Microphysics Model
===================================
The conversion rates among the moist hydrometeors are parameterized assuming that 

.. math::
   \frac{\partial N_{m}}{\partial D} = n_{m}(D_{m}) = N_{0m} exp(-\lambda _{m} D_{m})
   
   where :math:`N_{0m}` is the intercept parameter, :math:`D_{m}` is the diameters, and 
   
.. math::
   \lambda_{m} = (\frac{\pi \rho_{m} N_{0m}}{q_{m}\rho})^0.25
   
The total production rates include the contribution from aggregation, accretion, sublimation, melting, bergeron process, freezing and autoconversion are listed below without derivation, for details, please refer to Yuh-Lang Lin et al's paper (J. Climate Appl. Meteor, 22, 1065, 1983) and Marat F. Khairoutdinov and David A. Randall's paper (J. Atm Sciences, 607, 1983).


