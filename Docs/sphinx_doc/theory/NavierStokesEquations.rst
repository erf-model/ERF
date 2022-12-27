
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
:math:`R_v`, :math:`R_i`, and :math:`R_p` are the gas constants for dry air, water vapor, cloud ice, precipitating condensates, espectively. :math:`C_{pd}`, :math:`C_{pv}`, :math:`C_{pi}`, and :math:`C_{pp}` are the specific heat for dry air, water vapor, cloud ice, and precipitating condensates, respectively.

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

Buoyancy Force -- Paper
-------------------------------------------------------

One version of the buoyancy force is given in Marat F. Khairoutdinov and David A. Randall's paper (J. Atm Sciences, 607, 1983):

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx -\rho_0 \mathbf{g} ( \frac{T^\prime}{\bar{T}}
                 + 0.61 q_v^\prime - q_c - q_i - q_p - \frac{p^\prime}{\bar{p}} )

This can be derived by starting with

.. math::
   p = \rho (R_d q_d + R_v q_v) T = \rho R_d T (q_d + \frac{R_v}{R_d} q_v) =
        \rho R_d T [1 + (\frac{R_v}{R_d} − 1) q_v − q_c − q_i - q_p ]

where we have replaced :math:`q_d` by :math:`1 - q_v - q_c - q_i - q_p`.

Then, assuming the perturbations of :math:`p^\prime`, :math:`T^\prime`, and :math:`\rho^\prime`
are small compared with the total pressure, temperature, and density, respectively,
and :math:`\rho = \rho_d + \rho_v + \rho_c + \rho_i + \rho_p`
we can write

.. math::
   \frac{p^\prime}{p} = \frac{\rho^\prime}{\rho} + \frac{T^\prime}{T} + \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+(\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p)}

which allows us to write

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx \rho \left( \frac{p^\prime}{p} - \frac{T^\prime}{T} -
         \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+(\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p)} \right)

We can re-write

.. math::
     \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+ ( (\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p) ) } )
     \approx
     (\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime) (1 - ( (\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p) )

     \approx
     (\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime)

if we retain only first-order terms in the mixing ratios.

Thus the buoyancy term can be finally written as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx \rho \left( \frac{p^\prime}{p} - \frac{T^\prime}{T} -
         ( \frac{R_v}{R_d}-1) q_v^\prime + q_c^\prime + q_i^\prime + q_p^\prime \right)

If we assume that :math:`q_c = q_c^\prime`, :math:`q_i = q_i^\prime`, and :math:`q_p = q_p^\prime`
because :math:`\bar{q_i} = \bar{q_c} = \bar{q_p} = 0`,  and we replace :math:`\frac{R_v}{R_d}` by 1.61,
then this is identical to the expression at the top of this section.

Buoyancy Force -- Code
-------------------------------------------------------

The buoyancy force is implemented in ERF using the following formula

.. math::
   \mathbf{B} = -\rho_0 \mathbf{g} ( 0.61 q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime
                  + \frac{T^\prime}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) )


To implement the buoyance force term, we assume :math:`T_v = T (1 + (\frac{R_v}{R_d} − 1 ) q_v − q_c − q_i - q_p) \approx T`, then we have

.. math::
    p = \rho (R_d q_d + R_v q_v) T = \rho R_d T [1 + (\frac{R_v}{R_d} − 1) q_v − q_c − q_i - q_p ] = \rho R_d T_v


so the perturbation of :math:`\rho` can be written as

.. math::
   \frac{p^\prime}{p} = \frac{\rho^\prime}{\rho} + \frac{T_v^\prime}{T_v}


then :math:`\frac{\rho^\prime}{\rho}` is

.. math::
   \frac{\rho^\prime}{\rho} = \frac{p^\prime}{p} - \frac{T_v^\prime}{T_v}

if we ignore the term :math:`\frac{p^\prime}{p}`, the equation implemented can be written as

.. math::
   \frac{T_v^\prime}{T_v} \approx \frac{\bar{T} [ (\frac{R_v}{R_d}-1) (q_v-\bar{q_v}) - (q_c + q_i + q_p - \bar{q_c} - \bar{q_i} - \bar{q_p})] +
                           (T - \bar{T})[1+(\frac{R_v}{R_d}-1) \bar{q_v} - \bar{q_c} - \bar{q_i} - \bar{q_p} ]}{\bar{T}}


after reorganizing the terms, we get

.. math::
   \mathbf{B} = \mathbf{g} \rho^\prime = -\rho \mathbf{g} [ 0.61 (q_v - \bar{q_v}) - (q_c - \bar{q_c} + q_i - \bar{q_i} + q_p - \bar{q_p})
                  + \frac{T - \bar{T}}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) ]

where the overbar represents a horizontal average of the current state.

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
   erf.do_precip = true (to turn precipitation on)


