
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
we define the sum of all non-precipitating moist quantites to be :math:`\rho_T = \rho_v + \rho_c + \rho_i`.
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

  \sum_s _qs = 1

Potential temperature :math:`\theta` is defined as a function of temperature and specific entropy:

.. math::
  \theta (\eta, T) = T_r exp(\frac{\eta - \eta_0}{C_p})

where :math:`T_r` is the reference temperature, usually chosen as the temperature at the surface,
:math:`\eta_0` is the specific entropy at the reference temperature and pressure,
and :math:`\eta` is the specific entropy, defined for the mixture as

.. math::
   \eta = q_d \eta_d + q_v \eta_v + q_i \eta_i + q_c \eta_c + q_p \eta_p

where :math:`\eta_d`, :math:`\eta_v`, :math:`\eta_i`,
and :math:`\eta_c`, :math:`\eta_p` are the partial specific entropies for dry air, water vapor, cloud water, cloud ice,
and precipitates, respectively, and :math:`T_p`, is the reference temperature for the condensates:

.. math::
  \eta_d = C_{pd} ln (\frac{T}{T_r}) - R_d ln (\frac{p_d}{p_rd}) + \eta_{rd}

  \eta_v = C_{pv} ln (\frac{T}{T_r}) - R_v ln (\frac{p_v}{p_rv}) + \eta_{rv}

  \eta_i = C_i ln (\frac{T}{T_r}) + \eta_{ri}

  \eta_c = C_c ln (\frac{T}{T_r}) + \eta_{rc}

  \eta_p = C_p ln (\frac{T}{T_p}) + \eta_{rp}

where :math:`\eta_{rd}`, :math:`\eta_{rv}`, :math:`\eta_{ri}`, :math:`\eta_{rc}` and :math:`\eta_{rp}` are the specific entropy for dry air, water vapor, and ice, cloud, and precipitates, respectively. The Exner pressure :math:`\Pi` can be written as,

.. math::
  \Pi = C_p (\frac{p}{\alpha P_r})^\frac{R}{C_p}

and :math:`\theta`, :math:`p` can be expressed as

.. math::
   \theta = T (\frac{p_r}{p})^\frac{R}{C_p}

   p = p_r (\frac{\Pi}{C_p})^{\frac{C_p}{R}}

and :math:`p_r` is the reference pressure.

Governing Equations for Multispecies Atmospheric Flow
-------------------------------------------------------
We assume that all species have same average speed, and define the total potential temperature

.. math::
  \theta = \frac{\sum_s \theta_s}{\sum_s \rho_s} \approx \frac{\rho_d}{\rho} (\theta_d + q_v \theta_v + q_i \theta_i + q_c \theta_c).

Then the governing equations become

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} + \mathbf{F}_\rho)

  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u} + \mathbf{F}_u) -
          \frac{1}{1 + q_T + q_p}  \nabla p^\prime_d + \nabla \cdot \tau + \mathbf{F} + \delta_{i,3}\mathbf{B}

  \frac{\partial (\rho_d \theta)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta + F_{\theta}) + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \theta) + F_Q

  \frac{\partial (\rho_d C)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} C + \mathbf{F}_C) + \nabla \cdot (\rho_d \alpha_{C}\ \nabla C)

  \frac{\partial (\rho_d q_T)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_T +F_{q_{T}}) - Q

  \frac{\partial (\rho_d q_p)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_p + F_{q_{p}}) + Q

In this set of equations, the subgrid turbulent parameterization effects are included with fluxes
:math:`F_\rho`, :math:`F_u`, :math:`F_C`, :math:`F_{\theta}`, :math:`F_{q_{T}}`, :math:`F_{q_{r}}`.
:math:`\mathbf{F}` stands for the external force, and :math:`Q` and :math:`F_Q` represent the mass and energy transformation
of water vapor to/from water through condensation/evaporation, which is determined by the microphysics parameterization processes.
:math:`\mathbf{B}` is the buoyancy force,

.. math::
     \mathbf{B} = \rho_d^\prime \mathbf{g} \approx -\rho_0 \mathbf{g} ( \frac{T^\prime}{\bar{T}}
                 + 0.61 q_v^\prime - q_c - q_i - q_p - \frac{p^\prime}{\bar{p}} )

which is implemented as

.. math::
   \mathbf{B} = -\rho_0 \mathbf{g} ( 0.61 q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime
                  + \frac{T^\prime}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) )

where the overbar represents a horizontal average of the current state.

We note that we have assumed the energy transport for the different components are the same.

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

The total production rates including the contribution from aggregation, accretion, sublimation, melting,
bergeron process, freezing and autoconversion are listed below without derivation.
For details, please refer to Yuh-Lang Lin et al (J. Climate Appl. Meteor, 22, 1065, 1983) and
Marat F. Khairoutdinov and David A. Randall's (J. Atm Sciences, 607, 1983).
The implementation of microphysics model in ERF is similar to the that in the SAM code (http://rossby.msrc.sunysb.edu/~marat/SAM.html)

Accretion
------------------
There are several different type of accretional growth mechanisms that need to be included; these describe
the interaction of water vapor and cloud water with rain water.

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


Implementation of Moisture Model
===================================

The microphysics model takes potential temperature :math:`\theta`, total pressure :math:`p`, and dry air density :math:`\rho_d` as input,
and users can control the microphysics process by using

::

   erf.do_cloud = true (to turn cloud on)
   erf.do_smoke = true (to turn smoke physics on)
   erf.do_precip = true (to turn precipitation on)


