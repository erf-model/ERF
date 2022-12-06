
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

We consider a mixture of dry air :math:`\rho_d`, nonprecipitating water vapor :math:`\rho_v`, and assumed to be perfect idea gas with constant heat capacities :math:`C_{vd}`, :math:`C_{vv}`, :math:`C_{pd}`, :math:`C_{pv}`, and condensates :math:`\rho_p`, for examples cloud ice, and cloud water, that are incompressible with constant heat capacities :math:`C_l`, :math:`C_i`.

If we ignore the precipitation, and the volume occupied by the condensated water, then we will have

.. math::
  p = p_d + p_v = \rho_d R_d T + \rho_v R_v T

here :math:`p_d` and :math:`p_v` are the partial pressures of dry air and nonprecipitating water vapor, respectively. :math:`\rho_d` and :math:`\rho_v` are dry air density and nonprecipitating water vapor density, respectively.

In ERF, we select the dry air :math:`\rho_d` as the dominant component, and the others are sparse components :math:`\rho_s` with :math:`s = 1, ...., N`. The mixing ratios :math:`m_s` are defined as the mass density of species :math:`s` relative to dry air density :math:`\rho_d`, :math:`m_s=\frac{\rho_s}{\rho_d}`, therefore we can define:

.. math::
  \rho = \frac{\rho_d}{1-\sum_s q_s} = \rho_d (1-\sum_s m_s) = \sum_s \rho_s

  q_{s} = \frac{\rho_s}{\rho} = \frac{m_s}{1+\sum_s m_s}

  q_{d} = 1-\frac{\sum_s \rho_s}{\rho} = 1 - \sum_s q_s

Assuming the total nonprecipitating water :math:`q_T = q_v + q_c + q_i`, where :math:`q_v` is water vapor, :math:`q_c` is cloud water, and :math:`q_i` is cloud ice, respectively, and the total precipitating water :math:`q_p = q_{rain} + q_{snow} + q_{graupel}`, where :math:`q_{rain}` is rain, :math:`q_{snow}` is snow, :math:`q_{graupel}` is graupel, respectively. and :math:`\rho_d` is the density of the dry air.

The set of conservation equations for variables :math:`\rho_d`, :math:`q_T`, :math:`q_P`, :math:`\mathbf{u}`, :math:`C`, and :math:`\theta_s` are:

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u})
  
  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
  \rho_d \nabla p^\prime_d +\rho_d^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}

  \frac{\partial (\rho_s \theta_s)}{\partial t} &= - \nabla \cdot (\rho_s \mathbf{u} \theta_s + F_{\theta_s}) + \nabla \cdot ( \rho_s \alpha_{T}\ \nabla \theta_s) + F_Q

  \frac{\partial (\rho_d C)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} C) + \nabla \cdot (\rho_d \alpha_{C}\ \nabla C)

  \frac{\partial (\rho_d q_T)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_T +F_{q_{T}}) - Q

  \frac{\partial (\rho_d q_p)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_p + F_{q_{p}}) + Q

where :math:`F_{\theta_{s}}`, :math:`F_{q_{T}}`, :math:`F_{q_{r}}` are subgrid scalar fluxes. and :math:`Q` represents the transformation of cloud water and water vapor to rain water through condensation, and determined by the microphysics parameterization processes.

Assuming the energy transport for different components are the same, we can further simplify these set of equations for variables :math:`\rho_d`, :math:`q_T`, :math:`q_p`, :math:`\mathbf{u}`, :math:`C`, and :math:`\Theta`

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u})
  
  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
  \rho_d \nabla p^\prime_d +\rho_d^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}

  \frac{\partial (\rho_d \Theta)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \Theta + F_{\Theta}) + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \Theta) + F_Q

  \frac{\partial (\rho_d C)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} C) + \nabla \cdot (\rho_d \alpha_{C}\ \nabla C)

  \frac{\partial (\rho_d q_T)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_T +F_{q_{T}}) - Q

  \frac{\partial (\rho_d q_p)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_p + F_{q_{p}}) + Q
  

Single Moment Microphysics Model
===================================
The conversion rates among the moist hydrometeors are parameterized assuming that

.. math::
   \frac{\partial N_{m}}{\partial D} = n_{m}\left(D_{m}\right) = N_{0m} exp \left(-\lambda _{m} D_{m}\right)

where :math:`N_{0m}` is the intercept parameter, :math:`D_{m}` is the diameters, and

.. math::
   \lambda_{m} = (\frac{\pi \rho_{m} N_{0m}}{q_{m}\rho})^{0.25}

where :math:`\rho_{m}` is the density of moist hydrometeors. Assuming that the particle terminal velocity

.. math::
   v_{m} \left( D_{m},p \right) = a_{m}D_{m}^{b_{m}}\left(\frac{\rho_{0}}{\rho}\right)^{0.5}

The total production rates including the contribution from aggregation, accretion, sublimation, melting, bergeron process, freezing and autoconversion are listed below without derivation, for details, please refer to Yuh-Lang Lin et al (J. Climate Appl. Meteor, 22, 1065, 1983) and Marat F. Khairoutdinov and David A. Randall's (J. Atm Sciences, 607, 1983). The implementation of mcrophysics model is similar to the work of SAM code (http://rossby.msrc.sunysb.edu/~marat/SAM.html)

Accretion
------------------
There are several different type of accretional growth mechanisms needs to be included, it is involved the interaction of water vapor and cloud water with rain water.

The accretion of cloud water forms in either the dry or wet growth rate can be written as:

.. math::
   Q_{gacw} = \frac{\pi E_{GW}n_{0G}q_{c}\Gamma(3.5)}{4\lambda_{G}^{3.5}}(\frac{4g\rho G}{3C_{D}\rho})^{0.5}

The accretion of raindrops by accretion of cloud water is

.. math::
   Q_{racw} = \frac{\pi E_{RW}n_{0R}\alpha q_{c}\Gamma(3+b)}{4\lambda_{R}^{3+b}}(\frac{\rho_{0}}{\rho})^{1/2}

The bergeron Process
------------------------
The cloud water transform to snow by deposition and rimming can be written as

.. math::
   Q_{sfw} = N_{150}\left(\alpha_{1}m_{150}^{\alpha_{2}}+\pi E_{iw}\rho q_{c}R_{150}^{2}U_{150}\right)

Autoconversion
------------------------
The collision and coalesence of cloud water to from randrops is parameterized as following

.. math::
   Q_{raut} = \rho\left(q_{c}-q_{c0}\right)^{2}\left[1.2 \times 10^{-4}+{1.569 \times 10^{-12}N_{1}/[D_{0}(q_{c}-q_{c0})]}\right]^{-1}

Evaporation
------------------------
The evaporation rate of rain is

.. math::
   Q_{revp} = 2\pi(S-1)n_{0R}[0.78\lambda_{R}^{-2}+0.31S_{c}^{1/3}\Gamma[(b+5)/2]a^{1/2}\mu^{-1/2}(\frac{\rho_{0}}{\rho})^{1/4}\lambda_{R}^{(b+5)/2}](\frac{1}{\rho})(\frac{L_{v}^{2}}{K_{0}R_{w}T^{2}}+\frac{1}{\rho r_{s}\psi})^{-1}

