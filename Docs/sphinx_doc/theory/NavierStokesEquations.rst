
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
   \frac{\partial N_{m}}{\partial D} = n_{m}\left(D_{m}\right) = N_{0m} exp \left(-\lambda _{m} D_{m}\right)
   
where :math:`N_{0m}` is the intercept parameter, :math:`D_{m}` is the diameters, and 
   
.. math::
   \lambda_{m} = (\frac{\pi \rho_{m} N_{0m}}{q_{m}\rho})^{0.25}
   
The total production rates including the contribution from aggregation, accretion, sublimation, melting, bergeron process, freezing and autoconversion are listed below without derivation, for details, please refer to Yuh-Lang Lin et al (J. Climate Appl. Meteor, 22, 1065, 1983) and Marat F. Khairoutdinov and David A. Randall's (J. Atm Sciences, 607, 1983). The implementation of mcrophysics model is similar to the work of SAM code (http://rossby.msrc.sunysb.edu/~marat/SAM.html)

Aggregation
------------------------
The aggregation rate of ice to form snow due to collision-coalescence process can be written as:

.. math::
  P_{saut} = \alpha_{l}(l_{ci}-l_{I0})
  
where :math:`\alpha_{l}` is the rate coefficient, and :math:`l_{I0}` is the threshold amount for aggregation to occur.

And the rimed snow crystals aggregate to from graupel to form snow can be written as:

.. math::
  P_{gaut} = \alpha_{2}(l_{s}-l_{s0})
  
where :math:`\alpha_{2}` is the rate coefficient, and :math:`l_{s0}` is the mass threshold for snow.

Accretion
------------------
There are several different type of accretional growth mechanisms needs to be included, it is involved the interaction of snow with other classes of hydrometeors, and other classes of hydrometeors interaction with each other.

The accretion of cloud ice by snow can be written as

.. math::
   P_{saci} = \frac{\pi E_{si} n_{0S} cl_{CI}\Gamma(3+d)}{4\lambda_{S}^{3+d}}(\frac{\rho_{0}}{\rho})^{0.5}
   
where :math:`E_{si}` is the collection efficiency of the snow for cloud ice.

Similarly, the accretion of cloud water by snow is:

.. math::
   P_{sacw} = \frac{\pi E_{sw}n_{0S}cl_{CW} \Gamma (3+d)}{4\lambda_{S}^{3+d}}(\frac{\rho_0}{\rho})^{0.5}

The accretion of cloud ice by rain, which is a sink term for cloud ice, and a source term for snow or hail, can be written:

.. math::
   P_{raci} = \frac{\pi E_{ri}n_{0R}al_{CI} \Gamma(3+b)}{4\lambda_{R}^{3+b}}(\frac{\rho_{0}}{\rho})^{0.5}

The accretion of rain due to the presence of cloud ice is:

.. math::
   P_{iacr} = \frac{\pi E_{ri}n_{0R}al_{CIW} \Gamma(3+b)}{24M_{i}\lambda_{R}^{6+b}}(\frac{\rho_{0}}{\rho})^{0.5}

The accretion rate of rain for snow is:

.. math::
   P_{racs} = \pi E_{sr}n_{0R}n_{0S}|U_{R}-U_{S}|(\frac{\rho_{0}}{\rho})(\frac{5}{\lambda_{S}^{6}}+\frac{2}{\lambda_{S}^{5}\lambda_{R}^{2}}+\frac{0.5}{\lambda_{S}^{4}\lambda_{RR}^{3}})
   
The accretion rate of snow for rain is:

.. math::
   P_{racs} = \pi E_{sr}n_{0R}n_{0S}|U_{R}-U_{S}|(\frac{\rho_{w}}{\rho})(\frac{5}{\lambda_{S}^{6}}+\frac{2}{\lambda_{S}^{5}\lambda_{R}^{2}}+\frac{0.5}{\lambda_{S}^{4}\lambda_{RR}^{3}})

The accretion of cloud water forms in either the dry or wet growth rate can be written as:

.. math::
   P_{gacw} = \frac{\pi E_{GW}n_{0G}l_{CW}\Gamma(3.5)}{4\lambda_{G}^{3.5}}(\frac{4g\rho G}{3C_{D}\rho})^{0.5}
   
The accretion of cloud ice forms in either the dry or wet growth rate can be written as:

.. math::
   P_{gaci} = \frac{\pi E_{GW}n_{0G}l_{CI}\Gamma(3.5)}{4\lambda_{G}^{3.5}}(\frac{4g\rho G}{3C_{D}\rho})^{0.5}

The accretion of rain forms in either the dry and wet growth rate can be writeen as:

.. math::
   P_{gacr} = \pi^{2}E_{GR}n_{0G}n_{0R}|U_{G}-U_{R}|\left(\frac{\rho w}{\rho}\right)\left(\frac{5}{\lambda_{R}^{6}}+\frac{2}{\lambda_{R}^{5}\lambda_{G}^{2}}+\frac{0.5}{\lambda_{R}^{4}\lambda_{G}^{3}}\right)


Sublimation
------------------------
The depositional growth rate of snow is

.. math::
   P_{ssub} = \frac{2\pi(S_{i}-1)}{\rho (A^{''}+B^{''})} n_{0S}[0.78\lambda_{S}^{-2}+0.31S_{c}^{1/3}\Gamma(\frac{d+5}{2})c^{1/2}(\frac{\rho_{0}}{\rho})^{0.25}\mu^{-0.5}\lambda_{S}^{-(d+5)/2}]
   
where :math:`A^{''}=\frac{L_{S}^{2}}{K_{\alpha}R_{w}T^{2}}`, :math:`B^{''}=\frac{1}{\rho r_{si}\psi}`

Melting
------------------------
The rate of melting of snow to for rain can be expressed as

.. math::
   P_{smlt} = -\frac{2\pi}{\rho L_{j}} \left(K_{a}T_{c}-L_{v}\psi \rho \delta r_{s}\right) n_{0S} \left [0.78\lambda_{S}^{-2}+0.31 S_{c}^{1/3} \Gamma \left(\frac{d+5}{2}\right) c^{1/2} \left(\frac{\rho_{0}}{\rho}\right)^{1/4} \mu^{-1/2} \lambda_{S}^{-(d+5)/2} \right] - \frac{C_{w} T_{c}} {L_{f}} \left(P_{sacw)+P_{sacr}\right)

Freezing
------------------------
The raindrop freezing from raindrops to hails can be written as:

.. math::
   P_{gfr} = 20\pi^{2}B^{'}n_{0R}(\frac{\rho_{w}}{\rho}){exp[A^{'}(T_{0}-T]-1}\lambda_{R}^{-7}

Evaporation
------------------------
The evaporation rate of rain is

.. math::
   P_{revp} = 2\pi(S-1)n_{0R}[0.78\lambda_{R}^{-2}+0.31S_{c}^{1/3}\Gamma[(b+5)/2]a^{1/2}\mu^{-1/2}(\frac{\rho_{0}}{\rho})^{1/4}\lambda_{R}^{(b+5)/2}](\frac{1}{\rho})(\frac{L_{v}^{2}}{K_{0}R_{w}T^{2}}+\frac{1}{\rho r_{s}\psi})^{-1}

