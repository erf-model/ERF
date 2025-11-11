
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Microphysics:

Microphysics model
====================

Model overview and transported quantities in ERF
(note: ``Q1`` and ``Q2`` are always the mixing ratios of water vapor and cloud water)

+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Model              | Name in ERF            | ``Q3``      | ``Q4``      | ``Q5``      | ``Q6``      |
+====================+========================+=============+=============+=============+=============+
| Simple saturation  | ``SatAdj``             | --          | --          | --          | --          |
| adjustment         |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Kessler, no rain   | ``Kessler_NoRain``     | --          | --          | --          | --          |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Kessler            | ``Kessler``            | :math:`q_r` | --          | --          | --          |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Single moment,     | ``SAM_NoPrecip_NoIce`` | --          | --          | --          | --          |
| no precip or ice   |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Single moment,     | ``SAM_NoIce``          | --          | :math:`q_r` | --          | --          |
| no ice             |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Single moment      | ``SAM``                | :math:`q_i` | :math:`q_r` | :math:`q_s` | :math:`q_g` |
|                    |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Double moment,     | ``Morrison_NoIce``     | --          | :math:`q_r` | --          | --          |
| no ice             |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Double moment      | ``Morrison``           | :math:`q_i` | :math:`q_r` | :math:`q_s` | :math:`q_g` |
|                    |                        |             |             |             |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+
| Predicted Particle | ``P3``                 | :math:`q_i` | :math:`q_r` | :math:      | --          |
| Properties         |                        |             |             | `q_{rim}`   |             |
+--------------------+------------------------+-------------+-------------+-------------+-------------+


Kessler Microphysics model
---------------------------
The Kessler microphysics model is a simple version of cloud microphysics which has precipitation only in the form of rain. Hence :math:`q_p = q_r`.
Governing equations for the microphysical quantities for Kessler microphysics from `gabervsek2012dry`_ are

.. math::
    \frac{\partial q_v}{\partial t} = -C_c + E_c + E_r
.. math::
    \frac{\partial q_c}{\partial t} = C_c - E_c - (A_c + K_c)
.. math::
    \frac{\partial q_p}{\partial t} =  \frac{1}{\rho_{d}}\frac{\partial}{\partial z}(\rho_{d} w_{t} q_p) + (A_c + K_c) - E_r
.. math::
    \frac{\partial q_t}{\partial t} = \frac{\partial q_v}{\partial t} + \frac{\partial q_c}{\partial t}
                                    =  E_r - (A_c + K_c)

where :math:`C_c` is the rate of condensation of water vapor to cloud water, :math:`E_c` is the rate of evaporation of cloud water to water vapor,
:math:`A_c` is the autoconversion of cloud water to rain, :math:`K_c` is the accretion of cloud water to rain drops, :math:`E_r` is the evaporation of
rain to water vapor and :math:`F_r = \rho_{d} w_{t} q_p` is the sedimentation flux. The source terms that enter into the governing equations are then given by:

.. math::
   \mathbf{F_{n}} &\equiv [F_{q_v}, F_{q_c}] = \left[ -C_c, \;\; C_c \right],

   \mathbf{G_{p}} &= \left[ E_r, \;\; -A_c - K_c \right],

   H_{n} &= \rho_d \frac{L_v}{c_p} \frac{\theta_d}{T} C_c,

   F_{p} &= A_c + K_c - E_c,

   H_{p} &= -\rho_d \frac{L_v}{c_p} \frac{\theta_d}{T} E_r.

The parametrizations provided in `klemp1978simulation`_ are given below for each term.
Note that in all the equations, :math:`p` is specified in millibars and :math:`\overline{\rho}` is specified in g cm :math:`^{-3}`. The parametrization
of the source terms are given below.

.. _`gabervsek2012dry`: https://journals.ametsoc.org/view/journals/mwre/140/10/mwr-d-11-00144.1.xml

.. raw:: latex

   \newgeometry{left=2cm,right=2cm,top=2cm,bottom=2cm}

Rate of condensation of water vapor/evaporation of cloud water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From `klemp1978simulation`_, we have the following expressions.

.. _`klemp1978simulation`: https://journals.ametsoc.org/view/journals/atsc/35/6/1520-0469_1978_035_1070_tsotdc_2_0_co_2.xml

If the air is not saturated, i.e. :math:`q_v > q_{vs}`

.. math::
    C_c = \frac{q_v - q_{vs}}{1 + \cfrac{q_{vs}^*4093L}{c_p(T-36)^2}}

If the air is not saturated, i.e. :math:`q_v < q_{vs}`, then cloud water evaporates to water vapor at the rate

.. math::
    E_c = \frac{q_{vs} - q_v}{1 + \cfrac{q_{vs}^*4093L}{c_p(T-36)^2}}

Rate of autoconversion of cloud water into rain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The rate of autoconversion of cloud water into rain is given by

.. math::
    A_c = k_1(q_c - a)

where :math:`k_1 = 0.001` s\ :sup:`-1`, :math:`a = 0.001` kg kg\ :sup:`-1`.

Rate of accretion of cloud water into rain water drops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The rate of accretion of cloud water into rain water drops is given by

.. math::
    K_c = k_2q_cq_r^{0.875}

where :math:`k_2= 2.2` s\ :sup:`-1`.

The rate of evaporation of rain into water vapor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The rate of evaporation of rain into water vapor is given by

.. math::
    E_r = \cfrac{1}{\overline{\rho}}\cfrac{(1- q_v/q_s)C(\overline{\rho}q_r)^{0.525}}{5.4\times10^5 + 2.55\times10^6/(\overline{p}q_s)},

where the ventilation factor :math:`C = 1.6 + 124.9(\overline{\rho}q_r)^{0.2046}`.

Terminal fall velocity of rain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The terminal fall velocity of rain is given by

.. math::
    w_{t} = 3634(\overline{\rho}q_r)^{0.1346}\Bigg(\cfrac{\overline{\rho}}{\rho_0}\Bigg)^{-\frac{1}{2}}~\text{[cm/s]}

.. raw:: latex

   \restoregeometry


Morrison Microphysics Model
---------------------------------------

The Morrison microphysics model in ERF is a direct Fortran to C++ conversion of the Morrison
microphysics module in WRF.  For the relevant paper, please see Morrison et al, MWR, 2009.
The specific Fortran file which was ported was `module_mp_morr_two_moment.F`_

.. _`module_mp_morr_two_moment.F`: https://github.com/wrf-model/WRF/blob/master/phys/module_mp_morr_two_moment.F

Single Moment (SAM) Microphysics Model
---------------------------------------
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
The implementation of microphysics model in ERF is similar to the that in the SAM code (http://rossby.msrc.sunysb.edu/SAM.html)

Accretion
~~~~~~~~~~~~
There are several different type of accretional growth mechanisms that need to be included; these describe
the interaction of water vapor and cloud water with rain water.

The accretion of cloud water forms in either the dry or wet growth rate can be written as:

.. math::
   Q_{gacw} = \frac{\pi E_{GW}n_{0G}q_{c}\Gamma(3.5)}{4\lambda_{G}^{3.5}}(\frac{4g\rho G}{3C_{D}\rho})^{0.5}

The accretion of raindrops by accretion of cloud water is

.. math::
   Q_{racw} = \frac{\pi E_{RW}n_{0R}\alpha q_{c}\Gamma(3+b)}{4\lambda_{R}^{3+b}}(\frac{\rho_{0}}{\rho})^{1/2}

The bergeron Process
~~~~~~~~~~~~~~~~~~~~~~
The cloud water transform to snow by deposition and rimming can be written as

.. math::
   Q_{sfw} = N_{150}\left(\alpha_{1}m_{150}^{\alpha_{2}}+\pi E_{iw}\rho q_{c}R_{150}^{2}U_{150}\right)

Autoconversion
~~~~~~~~~~~~~~~~~~~~~~
The collision and coalescence of cloud water to from raindrops is parameterized as following

.. math::
   Q_{raut} = \rho\left(q_{c}-q_{c0}\right)^{2}\left[1.2 \times 10^{-4}+{1.569 \times 10^{-12}N_{1}/[D_{0}(q_{c}-q_{c0})]}\right]^{-1}

Evaporation
~~~~~~~~~~~~~~~~~~~~~~
The evaporation rate of rain is

.. math::
   Q_{revp} = 2\pi(S-1)n_{0R}[0.78\lambda_{R}^{-2}+0.31S_{c}^{1/3}\Gamma[(b+5)/2]a^{1/2}\mu^{-1/2}(\frac{\rho_{0}}{\rho})^{1/4}\lambda_{R}^{(b+5)/2}](\frac{1}{\rho})(\frac{L_{v}^{2}}{K_{0}R_{w}T^{2}}+\frac{1}{\rho r_{s}\psi})^{-1}


Saturation Adjustment (SatAdj) Microphysics Model
-------------------------------------------------
The saturation adjustment microphysics model is the simplest possible moisture model and only transports the
water vapor mixing ratio, :math:`q_v`, and the cloud water mixing ration, :math:`q_c`. Evaporation, :math:`q_v \longrightarrow q_c`, and condensation, :math:`q_c \longrightarrow q_v`, are the only relevant mechanisms. The final saturation state, :math:`q_v = q_{vs}(T)` is obtained from Newton-Raphson iterations on the thermal temperature :math:`T`.

Predicted Particle Properties (P3) Microphysics Model
------------------------------------------------------

The P3 microphysics scheme uses a fundamentally different approach than traditional bulk schemes.
Rather than using fixed hydrometeor categories (ice, snow, graupel), P3 predicts evolving ice particle
properties, allowing continuous transitions from unrimed ice to heavily rimed particles.

P3 transports water vapor (:math:`q_v`), cloud water (:math:`q_c`), rain (:math:`q_r`), total ice mass
(:math:`q_i`), and rime mass (:math:`q_{rim}`). Additional prognostic variables include ice number
concentration and rime volume.

The scheme represents physical processes including vapor deposition/sublimation, riming, aggregation,
melting, and sedimentation. Particle properties evolve continuously based on environmental conditions
and microphysical processes.

.. P3 requires ``USE_P3=TRUE`` at build time and interfaces with E3SM's P3 implementation.

For details, see Morrison and Milbrandt (2015, *J. Atmos. Sci.*, 72, 287–311).
