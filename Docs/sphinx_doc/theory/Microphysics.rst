
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Microphysics:

Microphysics model
====================

Model overview and transported quantities in ERF
(note: ``Q1`` and ``Q2`` are the mixing ratios of water vapor and cloud water for bulk models. For the Super-Droplet Method, ``Q1`` is water vapor and ``Q2`` is liquid cloud water)

+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Model              | Name in ERF             | ``Q3``      | ``Q4``      | ``Q5``          | ``Q6``      |
+====================+=========================+=============+=============+=================+=============+
| Virtual temp       | ``MoistNoCondensation`` | --          | --          | --              | --          |
| effects            |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Simple saturation  | ``SatAdj``              | --          | --          | --              | --          |
| adjustment         |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Kessler, no rain   | ``Kessler_NoRain``      | --          | --          | --              | --          |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Kessler            | ``Kessler``             | :math:`q_r` | --          | --              | --          |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Single moment,     | ``SAM_NoPrecip_NoIce``  | --          | --          | --              | --          |
| no precip or ice   |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Single moment,     | ``SAM_NoIce``           | --          | :math:`q_r` | --              | --          |
| no ice             |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Single moment      | ``SAM``                 | :math:`q_i` | :math:`q_r` | :math:`q_s`     | :math:`q_g` |
|                    |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Double moment,     | ``Morrison_NoIce``      | --          | :math:`q_r` | --              | --          |
| no ice             |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Double moment      | ``Morrison``            | :math:`q_i` | :math:`q_r` | :math:`q_s`     | :math:`q_g` |
|                    |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Predicted Particle | ``P3``                  | :math:`q_i` | :math:`q_r` | :math:`q_{rim}` | --          |
| Properties         |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+
| Super-Droplet      | ``SuperDroplets``       | :math:`q_i` | :math:`q_r` | :math:`q_s`     | :math:`q_g` |
| Method (SDM)       |                         |             |             |                 |             |
+--------------------+-------------------------+-------------+-------------+-----------------+-------------+


Surface precipitation accumulations
-----------------------------------

Each microphysics scheme that sediments precipitation keeps internal surface
accumulator fields for the mass that reaches the lower boundary. The available
species depend on the scheme: warm-rain configurations expose rain only, while
mixed-phase schemes may expose rain, snow, and graupel. Hail remains reserved
for a future scheme that carries a distinct hail accumulator.

The internal accumulator names are scheme-native. In some schemes, a historical
``rain`` accumulator actually represents total precipitation rather than
rain-only precipitation. The 2D output interface therefore treats precipitation
sources as typed total/species inputs and derives the public
liquid-water-equivalent fields from the normalized sources.

ERF's user-facing 2D precipitation diagnostics normalize those scheme-native
accumulators to a common liquid-water-equivalent basis with units of
``kg/m^2``. That lets downstream code compare schemes and diagnose coupling
interval totals without depending on a scheme's internal depth convention.

For an output interval from :math:`t_0` to :math:`t_1`, a frozen fraction can
be computed from the accumulation deltas:

.. math::

   f_{frozen} =
   \begin{cases}
   \dfrac{\Delta P_{frozen}}{\Delta P_{total}}, & \Delta P_{total} > 0 \\
   0, & \text{otherwise}
   \end{cases}

with

.. math::

   \Delta P_{total} = P_{total}(t_1) - P_{total}(t_0)

and

.. math::

   \Delta P_{frozen} = P_{snow}(t_1) - P_{snow}(t_0)
   + P_{graupel}(t_1) - P_{graupel}(t_0)
   + P_{hail}(t_1) - P_{hail}(t_0).


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
The saturation adjustment microphysics model is a warm-cloud, cell-local adjustment scheme that only transports
water vapor, :math:`q_v`, and cloud water, :math:`q_c`. It does not include rain, ice, sedimentation,
autoconversion, accretion, or subgrid cloud fraction. Pressure is diagnosed from the local cell state and the
saturation relation uses pressure in mbar (hPa).

In the normal physical branch used by SatAdj, the warm saturation relation is

.. math::
  q_{sat}(T,p) = \epsilon \frac{e_s(T)}{p - e_s(T)},

where :math:`\epsilon = R_d / R_v`, :math:`e_s(T)` is the saturation vapor pressure over liquid water, and
:math:`p` is the cell pressure in mbar (hPa).

During the adjustment, the conserved nonprecipitating water and local moist enthalpy proxy are

.. math::
  q_t = q_v + q_c,

.. math::
  T + \lambda q_v = \text{constant}, \qquad \lambda = L_v / c_p,

or equivalently,

.. math::
  T - \lambda q_c = \text{constant}.

The final state satisfies the complementarity conditions

.. math::
  q_c \ge 0,

.. math::
  q_v \le q_{sat}(T,p),

.. math::
  q_c \left(q_{sat}(T,p) - q_v\right) = 0.

When the available moisture is sufficient to reach saturation, the adjusted temperature is obtained from a
Newton solve of

.. math::
  F(T) = -T + T_i + \lambda \left(q_{v,i} - q_{sat}(T,p)\right) = 0,

with derivative

.. math::
  F'(T) = -1 - \lambda \frac{d q_{sat}}{dT}.

The implemented branch structure is

.. code-block:: text

  diagnose T, p, qv, qc from the cell state
  if qv + qc > qsat(T, p):
     if qc < 0: move negative qc into qv
     solve F(T) = 0 with Newton iteration
     set qv = qsat(T, p)
     set qc = qt - qv
  else:
     evaporate all qc into qv
     cool the cell by lambda * qc
     if the cooled state is now supersaturated:
        solve F(T) = 0 with Newton iteration
        set qv = qsat(T, p)
        set qc = qt - qv
  update theta from the adjusted T and diagnosed pressure

In the implementation, :math:`q_{sat}` and :math:`dq_{sat}/dT` are evaluated with ERF's internal thermodynamic
utilities. For warm-water saturation pressure, ERF uses the Flatau polynomial on the current positive selected
interval, about [-70, 70] C, and otherwise falls back to a Magnus-style closed-form exponential approximation.
That fallback is commonly motivated by simplified integrations of the Clausius-Clapeyron relation, but the
helper does not perform a direct Clausius-Clapeyron integration. Pressure is passed in mbar for saturation
calls, and pressure is converted back to Pa when recomputing :math:`\theta` from :math:`T`.

When SHOC is enabled, SatAdj condensation is disabled so that SHOC owns the phase-change adjustment and ERF
does not double-apply condensation tendencies.

Virtual Temperature Effects (MoistNoCondensation) Microphysics Model
--------------------------------------------------------------------
The virtual temperature effects moisture model is a simplification of the
saturation adjustment model above. It removes the condensation operation so
that water vapor, :math:`q_{v}`, will impact pressure, through the EOS, and
buoyancy but will not introduce thermal sources due to conversion to and from
cloud water. Consequently, this moisture model only transports :math:`q_{v}`.


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

Super-Droplet Method (SDM) Microphysics Model
----------------------------------------------

The super-droplet method (SDM) is a particle-based, probabilistic approach for the simulation of cloud microphysics.
Unlike the bulk parametrization and spectral bin methods, SDM directly tracks computational particles (called "super-droplets")
that represent multiple real droplets with identical attributes. This Lagrangian approach enables accurate simulation of
detailed cloud microphysical processes with reasonable computational cost.

.. note::

   To use SDM, build ERF with particles enabled i.e., ``USE_PARTICLES=TRUE``
   with GNU Make or ``-D ERF_ENABLE_PARTICLES=TRUE`` with CMake.

Overview and Method
~~~~~~~~~~~~~~~~~~~

The implementation in ERF is based on the method described in Shima et al., 2009 (Q. J. R. Meteorol. Soc., 135: 1307-1320)
for warm-rain processes, extended to include cold processes (ice, snow, and graupel formation) following the development
described in Shima et al., 2020 (Geosci. Model Dev., 13: 4107-4157). The key innovation is the concept of a super-droplet:
a computational particle that represents :math:`\xi_i(t)` identical physical droplets, where :math:`\xi_i(t)` is the multiplicity.
Each super-droplet has its own position :math:`\mathbf{x}_i(t)` and attributes :math:`\mathbf{a}_i(t)` that characterize
the state of the :math:`\xi_i(t)` physical droplets it represents.

For the warm-rain and mixed-phase cloud system implemented in ERF, the attributes are:

**For liquid droplets:**

- Equivalent radius of water, :math:`R_i(t)`, representing the amount of water the droplet contains
- Mass of solute contained in the droplet, :math:`M_i(t)`, for multiple aerosol species
- Freezing temperature, :math:`T^{fz}_i`, determined by ice nucleation active surface sites (INAS) on insoluble aerosol particles
- Additional species masses for multi-component systems

**For ice particles:**

- Equatorial radius :math:`a_i(t)` and polar radius :math:`c_i(t)` representing the ice particle's spatial extent (porous spheroid approximation)
- Apparent density :math:`\rho^i_i(t)` representing the particle's internal structure (accounting for porosity and branching)
- Rime mass :math:`m^{rime}_i(t)` recording mass obtained through riming
- Number of monomers :math:`n^{mono}_i(t)` representing the number of primary ice crystals in aggregates
- Freezing temperature :math:`T^{fz}_i` and aerosol masses (inherited from the droplet upon freezing)

The total number of physical droplets represented is :math:`N_r(t) = \sum_{i=1}^{N_s(t)} \xi_i(t)`, where :math:`N_s(t)`
is the number of super-droplets in the simulation domain.

Initialization of Super-Droplets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initialization process creates super-droplets in the computational domain with prescribed spatial and attribute distributions.
The implementation supports multiple initialization regions and allows flexible specification of particle properties.

Spatial Distribution
^^^^^^^^^^^^^^^^^^^^

Super-droplets can be initialized with two spatial distribution types:

**Uniform Distribution**: Particles are placed uniformly within a rectangular box defined by ``particle_box_lo`` and ``particle_box_hi``.
Within each grid cell that intersects the initialization box, the specified number of super-droplets per cell (``particles_per_cell``)
are created. If the initialization box is smaller than a grid cell (subgrid initialization), particles are placed only in cells
that contain the box.

**Bubble Distribution**: Particles are placed uniformly within an ellipsoidal bubble defined by ``particle_bubble_center`` (center position)
and ``particle_bubble_radius`` (radii in each direction). The bubble region is determined by:

.. math::
   \sqrt{\left(\frac{x-x_c}{r_x}\right)^2 + \left(\frac{y-y_c}{r_y}\right)^2 + \left(\frac{z-z_c}{r_z}\right)^2} \leq 1

where :math:`(x_c, y_c, z_c)` is the bubble center and :math:`(r_x, r_y, r_z)` are the radii.

Particle Position
^^^^^^^^^^^^^^^^^

Within each grid cell, super-droplet positions can be assigned in two ways controlled by ``place_randomly_in_cells``:

- **Random placement** (default): Particles are randomly distributed within the cell or initialization region using a uniform random distribution.
  For bubble distributions, random positions are generated using spherical coordinates with uniform angular distribution and radial
  position sampled to ensure uniform volume distribution.

- **Fixed placement**: Particles are placed at the cell center or initialization region center.

For terrain-following coordinates, the vertical particle position is stored as the computational coordinate, so no terrain
adjustment is applied at initialization; the mapping to physical height is performed where the particle interacts with the mesh.

Attribute Distribution
^^^^^^^^^^^^^^^^^^^^^^

Each super-droplet is assigned physical attributes including masses of water vapor/condensate species and aerosol species.
The mass distribution for each species and aerosol can be specified independently using one of four distribution types:

**mass_constant**: All super-droplets have the same mass for the species, set to ``species_mean_mass_<NAME>`` or ``aerosol_mean_mass_<NAME>``.

**mass_exponential**: Masses are sampled from an exponential distribution with mean ``species_mean_mass_<NAME>`` or ``aerosol_mean_mass_<NAME>``,
truncated between the minimum and maximum mass values.

**radius_log_normal**: Dry radii are sampled from a log-normal distribution with mean radius ``species_mean_radius_<NAME>`` or
``aerosol_mean_radius_<NAME>`` and geometric standard deviation ``species_geomstd_radius_<NAME>`` or ``aerosol_geomstd_radius_<NAME>``.
The mass is computed from the radius using the species density.

**radius_lognormal_autorange**: Similar to ``radius_log_normal``, but the minimum and maximum radii are automatically determined
from the distribution parameters to capture the specified range of the distribution.

Multiplicity Assignment
^^^^^^^^^^^^^^^^^^^^^^^

The multiplicity :math:`\xi_i` (number of physical droplets represented by each super-droplet) can be assigned using two methods
specified by ``multiplicity_type``:

**Constant multiplicity**: Each super-droplet represents the same number of physical particles. The multiplicity is computed as:

.. math::
   \xi = \left\lceil \frac{N_\text{par}}{N_\text{SD}} \right\rceil

where :math:`N_\text{par}` is the total number of physical particles per cell (from ``initial_number_density`` times cell volume)
and :math:`N_\text{SD}` is the number of super-droplets per cell (``particles_per_cell``).

**Sampled multiplicity**: Multiplicities are sampled from the mass/radius distribution to better represent the underlying
size distribution. The sampled multiplicities are scaled to ensure the total number of physical particles matches the specified
number density. This approach can provide better statistical representation of the particle size distribution, especially for
broad distributions.

Effective Radius and Total Mass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After assigning species and aerosol masses, each super-droplet's effective radius and total mass are computed. The effective
radius accounts for the water content, soluble species (which affect the droplet solution properties), and insoluble species
(which form cores within the droplet):

.. math::
   r_\text{eff} = \left(\frac{m_w + m_s + \frac{\rho_w}{\rho_{ins}}m_{ins}}{\frac{4}{3}\pi\rho_w}\right)^{1/3}

where:
- :math:`r_\text{eff}` is the effective radius of the droplet
- :math:`m_w` is the water mass
- :math:`m_s` is the total mass of soluble species and aerosols
- :math:`m_{ins}` is the total mass of insoluble species and aerosols
- :math:`\rho_w` is the water density
- :math:`\rho_{ins}` is the weighted average density of insoluble components

The total mass stored in the super-droplet includes all species and aerosol masses:

.. math::
   m_\text{total} = m_w + \sum_{j=1}^{N_\text{sp}} m_{\text{sp},j} + \sum_{j=1}^{N_\text{ae}} m_{\text{ae},j}

where:
- :math:`m_\text{total}` is the total mass of the droplet
- :math:`N_\text{sp}` is the number of species
- :math:`m_{\text{sp},j}` is the mass of species j
- :math:`N_\text{ae}` is the number of aerosol types
- :math:`m_{\text{ae},j}` is the mass of aerosol type j

Initialization from Condensate Density
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An alternative initialization method (``distribution_type = condensate_density``) allows super-droplet attributes to
be set from a prescribed condensate mass density field. In this case, the multiplicity of each super-droplet is varied randomly
around the mean value, and the mass per physical particle is computed to match the local condensate density. The radius is then
determined from the mass assuming spherical water droplets.

Density Scaling
^^^^^^^^^^^^^^^

If ``density_scaling`` is enabled, super-droplet multiplicities are scaled by the local air density after initialization. This
allows the number of physical droplets to vary with altitude, reflecting realistic atmospheric conditions where particle
concentrations typically decrease with height.

Multiple Initializations
^^^^^^^^^^^^^^^^^^^^^^^^^

The implementation supports multiple initialization regions (``num_initializations``), each with its own spatial distribution,
particle count, and attribute distributions. Parameters for each region are specified with the prefix ``super_droplets_moisture.N``
where ``N`` is the initialization index (0, 1, 2, ...). This allows complex initial conditions with multiple cloud regions or
layers with different properties.

Microphysical Processes
~~~~~~~~~~~~~~~~~~~~~~~~

Motion and Sedimentation
^^^^^^^^^^^^^^^^^^^^^^^^

Droplets are assumed to immediately reach their terminal velocity. The motion of each super-droplet is governed by:

.. math::
   \mathbf{v}_i(t) &= \mathbf{U}(\mathbf{x}_i) - \hat{\mathbf{z}}v_\infty(R_i), \\
   \frac{d\mathbf{x}_i}{dt} &= \mathbf{v}_i,

where:

- :math:`\mathbf{v}_i(t)` is the velocity vector of the super-droplet
- :math:`\mathbf{U}(\mathbf{x}_i)` is the wind velocity at the droplet position
- :math:`\hat{\mathbf{z}}` is the unit vector in the vertical direction
- :math:`v_\infty(R_i)` is the terminal velocity of the droplet
- :math:`R_i` is the radius of the droplet
- :math:`\mathbf{x}_i` is the position of the super-droplet

The terminal velocity :math:`v_\infty` depends on the particle phase. Liquid drops follow the three-regime
cloud-and-rain parameterization of Beard (1976) (as implemented in SCALE-SDM), in which the velocity is recovered
from the Reynolds number :math:`N_{Re}` of the falling drop:

.. math::
   v_\infty = \frac{\mu\, N_{Re}}{\rho\, d}, \qquad d = 2R,

where :math:`\mu` is the dynamic viscosity of air, :math:`\rho` is the air density, and :math:`d` is the drop
diameter. Small cloud droplets (:math:`d \lesssim 19\,\mu\mathrm{m}`) follow Stokes' law with a Cunningham slip
correction:

.. math::
   v_\infty = \frac{(\rho_w - \rho)\, g\, d^2}{18\,\mu}\left(1 + 2.51\,\frac{\lambda}{d}\right),

with :math:`\lambda` the mean free path of air. For larger drops the Reynolds number is recovered from the Davies
(Best) number :math:`N_{Da} = 4\rho(\rho_w-\rho)g\,d^3/(3\mu^2)`, for cloud droplets and small raindrops, or from
the Bond number :math:`N_{Bo} = 4(\rho_w-\rho)g\,d^2/(3\sigma)` together with a physical-property group, for large
raindrops, through the empirical fits of Beard (1976); here :math:`\sigma` is the surface tension of water.

The measured fall speed of raindrops levels off toward a plateau rather than increasing without bound (Gunn and
Kinzer, 1949). The drop radius is therefore capped at :math:`R \le 4\,\mathrm{mm}` when evaluating the liquid-drop
fall speed, holding the speed at this plateau and avoiding the unphysical divergence of the large-drop fit at very
large radii. The simpler Rogers-Yau and Atlas-Ulbrich formulae are also available. Frozen and mixed-phase particles
fall at the Böhm ice speed and a meltwater-weighted blend, described among the cold microphysical processes below.

Condensation and Evaporation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The growth/shrinkage of droplets through condensation/evaporation is described by Köhler theory, which accounts for both
curvature and solution effects:

.. math::
   R_i \frac{dR_i}{dt} = \frac{(S_v-1) - \frac{a_K}{R_i} + \frac{b_K}{R_i^3}}{F_k + F_d},

where:

- :math:`R_i` is the radius of the droplet
- :math:`S_v` is the ambient saturation ratio (ratio of vapor pressure to saturation vapor pressure)
- :math:`a_K` is the curvature coefficient, given by :math:`a_K = \frac{2\sigma}{R_v \rho_w T}`, where:
  - :math:`\sigma` is the surface tension of water
  - :math:`R_v` is the specific gas constant for water vapor
  - :math:`\rho_w` is the density of liquid water
  - :math:`T` is the temperature
- :math:`b_K` is the solution coefficient, given by :math:`b_K = \frac{3 i M_i M_w}{4\pi \rho_w M_s}`, where:
  - :math:`i` is the van't Hoff factor (degree of ionic dissociation)
  - :math:`M_i` is the mass of solute in the droplet
  - :math:`M_w` is the molecular weight of water
  - :math:`M_s` is the molecular weight of the solute
- :math:`F_k` and :math:`F_d` are thermodynamic and diffusion terms:

.. math::
   F_k &= \left(\frac{L_v}{R_v T} - 1\right) \frac{L_v\rho_{w}}{K_T T}, \\
   F_d &= \frac{\rho_{w} R_v T}{D_v e_s(T)},

where:
- :math:`L_v` is the latent heat of vaporization for water
- :math:`K_T` is the thermal conductivity of air
- :math:`D_v` is the molecular diffusion coefficient of water vapor in air
- :math:`e_s(T)` is the saturation vapor pressure at temperature :math:`T`

The implementation uses an implicit time integration scheme (backward Euler or higher-order methods like Runge-Kutta)
with a Newton-Raphson solver to solve this nonlinear ordinary differential equation.

Stochastic Coalescence
^^^^^^^^^^^^^^^^^^^^^^^

Coalescence (collision and merging) of droplets is treated probabilistically. For two super-droplets :math:`j` and :math:`k`
in a well-mixed volume :math:`\Delta V`, the coalescence probability during time interval :math:`\Delta t_c` is:

.. math::
   P_{jk}^{(s)} = \max(\xi_j, \xi_k) \frac{K_{coal}(R_j, R_k) \Delta t_c}{\Delta V},

where:
- :math:`P_{jk}^{(s)}` is the probability of collision between super-droplets j and k
- :math:`\xi_j` and :math:`\xi_k` are the multiplicities of the super-droplets
- :math:`K_{coal}(R_j, R_k)` is the coalescence kernel
- :math:`\Delta t_c` is the coalescence time step
- :math:`\Delta V` is the volume in which coalescence is considered

The coalescence kernel is given by:

.. math::
   K_{coal}(R_j, R_k) = E_{coll}(R_j, R_k) \pi(R_j + R_k)^2 |\mathbf{v}_j - \mathbf{v}_k|

where:
- :math:`E_{coll}(R_j, R_k)` is the collection efficiency accounting for flow deflection and droplet bounce effects
- :math:`R_j` and :math:`R_k` are the radii of the two droplets
- :math:`\mathbf{v}_j` and :math:`\mathbf{v}_k` are the velocities of the two droplets

When super-droplets :math:`j` and :math:`k` coalesce (assuming :math:`\xi_j > \xi_k` without loss of generality),
:math:`\min(\xi_j, \xi_k)` pairs of physical droplets merge:

.. math::
   \xi'_j &= \xi_j - \xi_k, \quad \xi'_k = \xi_k, \\
   R'_j &= R_j, \quad R'_k = (R_j^3 + R_k^3)^{1/3}, \\
   M'_j &= M_j, \quad M'_k = M_j + M_k.

This formulation preserves the number of super-droplets in most cases, maintaining accuracy even as the number of physical
droplets changes dramatically.

The Monte Carlo scheme for coalescence achieves :math:`O(n_s)` computational cost (where :math:`n_s` is the number of
super-droplets in a grid cell) by examining only :math:`[n_s/2]` randomly generated, non-overlapping pairs rather than
all possible :math:`n_s(n_s-1)/2` pairs. The coalescence probability is scaled accordingly.

Several collision kernels are implemented:

- **Golovin kernel**: :math:`K = b(X_j + X_k)` where :math:`X = \frac{4\pi}{3}R^3` (primarily for testing)
- **Sedimentation kernel**: Based on geometric collision cross-section
- **Long's kernel**: Empirical formula for gravitational collision
- **Hall's kernel**: Table-based collection efficiency for realistic cloud conditions

Cold Microphysical Processes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SDM implementation includes cold processes that lead to ice, snow, and graupel formation, following the porous spheroid
approximation approach described in Shima et al., 2020.

**Ice Particle Representation**

Ice particles are approximated as porous spheroids characterized by:

- Equatorial radius :math:`a` and polar radius :math:`c` defining the spatial extent
- Aspect ratio :math:`\phi = c/a` (prolate spheroids for columns with :math:`\phi > 1`, oblate spheroids for plates/dendrites with :math:`\phi < 1`)
- Apparent density :math:`\rho^i` representing internal structure (hollowing, branching, riming)
- Volume :math:`V = \frac{4\pi}{3}a^2c` and mass :math:`m = \rho^i V`

**Immersion/Condensation and Homogeneous Freezing**

Supercooled droplets freeze when the ambient temperature falls below their freezing temperature :math:`T^{fz}`. The freezing
temperature is determined by the ice nucleation active surface site (INAS) density on insoluble aerosol particles using
the singular hypothesis framework. Upon freezing, the resulting ice crystal is initially spherical with true ice density
:math:`\rho^i_{true} \approx 916.8` kg m\ :sup:`-3`.

**INAS-Based Freezing Temperature Initialization**

The freezing temperature (:math:`T^{fz}`) of each droplet is initialized using the ice nucleation active site (INAS) density
parameterization from Niemand et al. (2012). This approach assigns a freezing temperature based on the total surface area of
ice-nucleating particles (INPs) contained within the droplet:

1. The total surface area of INPs is computed for each superdroplet:

.. math::
   A_\text{INP} = \sum_{i=1}^{N_\text{sp}} \frac{m_{i}}{\rho_{i}} \cdot \frac{3}{r_{i}} \cdot \mathcal{I}_\text{INP,i} + \sum_{j=1}^{N_\text{ae}} \frac{m_{j}}{\rho_{j}} \cdot \frac{3}{r_{j}} \cdot \mathcal{I}_\text{INP,j}

where :math:`m_i` is the mass of species :math:`i`, :math:`\rho_i` is its density, :math:`r_i` is the characteristic radius of the
particle, and :math:`\mathcal{I}_\text{INP,i}` is an indicator function that is 1 if the species is an INP and 0 otherwise.

2. The INAS density parameterization from Niemand et al. (2012) is used to calculate the freezing temperature. The INAS density
:math:`n_s(T)` (m\ :sup:`-2`) is given by:

.. math::
   \log_{10} n_s(T) = a \cdot (T - 273.15) + b

where :math:`a = -0.517` and :math:`b = 8.934` are fitted parameters for mineral dust.

3. The probability of freezing is:

.. math::
   P_\text{freeze}(T) = 1 - \exp(-A_\text{INP} \cdot n_s(T))

4. Using inverse transform sampling with a uniform random number :math:`u \in [0,1]`, the freezing temperature is calculated as:

.. math::
   T^{fz} = 273.15 + \frac{1}{a} \left[ \log_{10}\left(\frac{-\log(1-u)}{A_\text{INP}}\right) - b \right]

During water-water coalescence events, the freezing temperature of the resulting droplet is set to:

.. math::
   T^{fz}_\text{new} = \max(T^{fz}_1, T^{fz}_2)

This preserves the highest freezing temperature (most ice-active) property during droplet collisions.

**Deposition and Sublimation**

Ice particle growth and sublimation through vapor diffusion follows:

.. math::
   \frac{dm_i}{dt} = 4\pi C D_v(\rho_{vi} - \rho^{sfc}_{vi})f_{vnt} = 4\pi C\frac{S^i_i - 1}{F^i_k + F^i_d}f_{vnt}

where :math:`C(a,c)` is the electric capacitance of the spheroid, :math:`S^i_i` is the saturation ratio over ice,
:math:`f_{vnt}` is the ventilation coefficient, and :math:`F^i_k` and :math:`F^i_d` are thermodynamic and diffusion terms.

The **primary growth habit** (columnar vs planar preference) is controlled by the inherent growth ratio :math:`\Gamma(T)`:

.. math::
   \frac{dc_i}{da_i} = \Gamma(T_i)f_{vnt}\frac{c_i}{a_i} =: \Gamma^*\frac{c_i}{a_i}

The **secondary growth habit** (solid vs hollow, plate vs dendritic) is expressed through the deposition density :math:`\rho_{dep}`,
which represents the apparent density of newly deposited ice. For deposition:

.. math::
   dV_i = \frac{dm_i}{\rho_{dep}}, \quad \text{for } dm_i \geq 0

For sublimation, the particle density is preserved: :math:`\rho_{sbl} = \rho^i_i`.

**Melting**

When the ambient temperature rises above :math:`0^\circ\mathrm{C}`, the ice core melts gradually, following the
time-dependent heat balance of Seifert and Beheng (2006). The melt rate of the ice core is:

.. math::
   \frac{dm_i}{dt} = -\frac{4\pi C}{L_f}\left[ K\, f_{vnt}\, (T - T_0)
       + \frac{L_v D_v}{R_v}\, f_{vnt} \left(\frac{e_\infty}{T} - \frac{e_{s,w}(T_0)}{T_0}\right) \right],

where :math:`C(a,c)` is the capacitance of the spheroid, :math:`L_f` and :math:`L_v` are the latent heats of fusion
and vaporization, :math:`K` is the thermal conductivity of air, :math:`D_v` is the vapor diffusivity, :math:`R_v` is
the vapor gas constant, :math:`T_0 = 273.15` K is the melting point, :math:`e_\infty` is the ambient vapor pressure,
and :math:`e_{s,w}(T_0)` is the saturation pressure over water at :math:`T_0`. The two bracketed terms are the
sensible heat conducted from the warm air and the latent heat carried by the vapor flux to the wet surface.

The meltwater accumulates on the ice rather than being shed immediately, so the super-droplet becomes a mixed
(wet) ice-water particle. The ice core shrinks while its aspect ratio and apparent density are preserved, and the
accumulated water is released as a pure raindrop only once the ice core has melted completely. A mixed particle
that re-enters subfreezing air refreezes its liquid film back onto the core. This gradual treatment retains, for
each particle, the mixed-phase state through which it falls, and differs from the instantaneous melting of earlier
super-droplet ice schemes (Shima et al., 2020).

**Riming**

Riming occurs when supercooled droplets are collected by ice particles and freeze on contact. The collision-riming kernel is:

.. math::
   K_{rime} = E_{rime} A_g |v^\infty_j - v^\infty_k|

where :math:`E_{rime}` is the collection efficiency and :math:`A_g` is the geometric cross-sectional area.

The frozen droplet is characterized by a rime density :math:`\rho_{rime}` that depends on impact velocity and surface temperature.
The filling-in approximation is used: the maximum dimension is preserved as unrimed crystals gradually grow into quasi-spherical
graupel particles. For graupels (aspect ratio :math:`0.8 < \phi \leq 1.25`), the model accounts for tumbling by preserving
the minor dimension during riming.

**Aggregation**

Ice particle aggregation is treated with a collision-aggregation kernel:

.. math::
   K_{agg} = E_{agg}\left(A^{1/2}_j + A^{1/2}_k\right)^2 |v^\infty_j - v^\infty_k|

where :math:`E_{agg} = 0.1` is the aggregation efficiency and :math:`A` is the projected area.

When two ice particles aggregate, the model accounts for compaction of fluffy snowflakes. The resultant aggregate's
apparent density is interpolated between a minimum density (assuming no shape change and more empty space) and the
volume-weighted average density (assuming complete compaction), based on how close the initial density is to a critical
compaction density :math:`\rho^i_{crt} = 10` kg m\ :sup:`-3`.

**Terminal Velocity of Ice and Mixed Particles**

The fall speed of a frozen particle uses the Böhm (1989, 1992, 1999) drag relation in the porous-spheroid form of
Shima et al. (2020). As for liquid drops, the velocity follows from the Reynolds number,

.. math::
   v_\infty = \frac{\mu\, N_{Re}}{\rho\, d}, \qquad d = 2 a_i,

where the characteristic diameter :math:`d = 2 a_i` is the equatorial extent of the spheroid (Böhm's definition).
The Reynolds number is recovered from the Best (Davies) number

.. math::
   X = \frac{8\, m\, g\, \rho}{\pi\, \mu^2\, \max(\phi_i, 1)\, \max(q_i^{1/4}, q_i)},

where :math:`m` is the particle mass, :math:`g` is gravity, :math:`\phi_i = c_i/a_i` is the aspect ratio, and
:math:`q_i = A_i/A^{ce}_i` is the area ratio. The projected area :math:`A_i` accounts for porosity,

.. math::
   A_i = A^{ce}_i \left(\frac{\rho^i_i}{\rho^i_{true}}\right)^{\kappa(\phi_i)}, \qquad \kappa(\phi_i) = \exp(-\phi_i),

with :math:`A^{ce}_i = \pi a_i \max(a_i, c_i)` the circumscribed ellipse area, so that
:math:`q_i = (\rho^i_i/\rho^i_{true})^{\exp(-\phi_i)}`. At large values the Best number is corrected as

.. math::
   X' = X\,\frac{1 + (X/X_0)^2}{1 + 1.6\,(X/X_0)^2}, \qquad X_0 = 2.8\times 10^6,

and :math:`N_{Re}(X')` then follows from the Böhm drag law, whose coefficients depend on :math:`\phi_i` and
:math:`q_i`. The dependence on mass, aspect ratio, and porosity captures the slower fall of non-spherical,
low-density ice.

A mixed (partially melted) particle carries both an ice core and a liquid film, and falls faster than its dry ice
frame but slower than the equivalent raindrop. Following Frick et al. (2013, Geosci. Model Dev., 6: 1925-1939), its
terminal velocity is interpolated between these two limits by the liquid mass fraction :math:`\ell = m_w/(m_w + m_i)`:

.. math::
   v_\infty = v^{ice}_\infty + \left(v^{rain}_\infty - v^{ice}_\infty\right)\Psi(\ell), \qquad
   \Psi(\ell) = 0.246\,\ell + 0.754\,\ell^7,

where :math:`m_w` and :math:`m_i` are the meltwater and ice masses, :math:`v^{ice}_\infty` is the Böhm fall speed
evaluated for the ice core alone (mass :math:`m_i`), :math:`v^{rain}_\infty` is the Beard (1976) fall speed of a
raindrop of the combined mass :math:`m_w + m_i`, and the weight :math:`\Psi(\ell)` is an empirical fit to the
wind-tunnel data of Mitra et al. (1990). The weight rises from zero for dry ice to unity once the core has melted.
Evaluating the ice fall speed from the ice-core mass alone keeps the mass consistent with the drag area of the
spheroid as the core melts, so the fall speed stays bounded.

Coupling with Eulerian Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The super-droplets are coupled to the non-hydrostatic Eulerian dynamics through three source terms:

**Momentum coupling** through the liquid water density:

.. math::
   \rho_{lw}(\mathbf{x}, t) = \sum_{i=1}^{N_s} \xi_i m_i(t) \delta^3[\mathbf{x} - \mathbf{x}_i(t)],

where:
- :math:`\rho_{lw}(\mathbf{x}, t)` is the liquid water density at position :math:`\mathbf{x}` and time :math:`t`
- :math:`N_s` is the total number of super-droplets
- :math:`\xi_i` is the multiplicity of the super-droplet (number of physical droplets it represents)
- :math:`m_i(t) = \frac{4\pi}{3} R_i^3 \rho_{w}` is the mass of a physical droplet
- :math:`\rho_{w}` is the density of liquid water
- :math:`\delta^3` is the three-dimensional delta function

**Vapor source** from condensation/evaporation:

.. math::
   S_{vap}(\mathbf{x}, t) = -\frac{1}{\rho_{air}(\mathbf{x}, t)} \sum_{i=1}^{N_s} \xi_i \frac{dm_i(t)}{dt} \delta^3[\mathbf{x} - \mathbf{x}_i(t)].

where:
- :math:`S_{vap}(\mathbf{x}, t)` is the vapor source/sink term
- :math:`\rho_{air}(\mathbf{x}, t)` is the air density
- :math:`\frac{dm_i(t)}{dt}` is the rate of mass change of a physical droplet

**Latent heat release**:

.. math::
   \frac{L_v}{c_p} S_{vap}

where:
- :math:`L_v` is the latent heat of vaporization
- :math:`c_p` is the specific heat capacity at constant pressure

In the numerical implementation, these delta functions are smoothed onto the Eulerian grid using appropriate
interpolation functions to ensure conservation and numerical stability.

Input Parameters and Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The super-droplet method is configured through input parameters specified in the input file with the prefix ``super_droplets_moisture``.

SDM Parameters
^^^^^^^^^^^^^^

All SDM parameters use the prefix ``super_droplets_moisture``:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Parameter
     - Default Value
     - Description
   * - **Microphysics Model**
     -
     -
   * - ``include_phase_change``
     - ``true``
     - Enable/disable condensation and evaporation
   * - ``include_advection``
     - ``true``
     - Enable/disable particle advection
   * - ``include_coalescence``
     - ``true``
     - Enable/disable stochastic coalescence
   * - ``include_cold_processes``
     - ``true``
     - Enable/disable ice/snow/graupel formation (cold processes)
   * - ``distribution_type``
     - ``uniform``
     - Initial distribution type (``uniform``, ``condensate_density``)
   * - ``radius_raindrop``
     - ``4.0e-5``
     - Minimum radius (m) to classify droplet as rain
   * - ``rime_mass_ratio``
     - ``0.3``
     - Threshold rime mass ratio for classifying ice particles as graupel
   * - ``kinematic_mode``
     - ``false``
     - Run in kinematic mode (no feedback to dynamics)
   * - ``dimensionality``
     - ``three_d``
     - Simulation dimensionality (``one_d_z``, ``two_d_xz``, ``two_d_yz``, ``three_d``)
   * - ``recycle_particles``
     - ``false``
     - Enable particle recycling at boundaries
   * - ``species``
     - ``H2O``
     - List of vapour/condensate species
   * - ``aerosols``
     - (none)
     - List of aerosol species
   * - ``diagnostics_interval``
     - ``1``
     - Timesteps between diagnostic output
   * - ``num_substeps_phase_change``
     - ``1``
     - Number of substeps for phase change process
   * - ``initial_phase_change_relaxation``
     - ``false``
     - Allow initial relaxation of droplet sizes
   * - ``initial_phase_change_relaxation_time``
     - ``10.0``
     - Duration (s) of initial relaxation
   * - **Particle Container**
     -
     -
   * - ``density_scaling``
     - ``false``
     - Scale initial SD number with air density
   * - ``nucleate_particles``
     - ``false``
     - Nucleate new super-droplets from vapor
   * - ``advect_with_flow``
     - ``true``
     - Advect particles with fluid velocity
   * - ``advect_with_gravity``
     - ``true``
     - Include gravitational settling
   * - ``prescribed_advection``
     - ``false``
     - Use prescribed vertical velocity
   * - ``coalescence_kernel``
     - ``sedimentation``
     - Coalescence kernel type (``golovin``, ``sedimentation``, ``Longs``, ``Halls``)
   * - ``terminal_velocity_model``
     - ``CloudRainShima`` (water), ``IceBohm`` (ice)
     - Terminal velocity formula. For water droplets: ``RogersYau``, ``AtlasUlbrich``, ``CloudRainShima``. For ice particles: ``IceBohm``. Can specify 1 value (applies to both) or 2 values (water, ice)
   * - ``include_brownian_coalescence``
     - ``false``
     - Include Brownian motion in coalescence
   * - ``coalescence_bin_size``
     - ``[1,1,1]``
     - Grid coarsening factor for coalescence cells
   * - ``mass_change_cfl``
     - ``1000.0``
     - CFL number for phase change time integration
   * - ``mass_change_ti_method``
     - ``backward_euler``
     - Time integrator for mass change ODE (``rk3bs``, ``rk4``, ``backward_euler``, ``crank_nicolson``, ``dirk2``)
   * - ``newton_solver_rtol``
     - ``1.0e-6``
     - Newton solver relative tolerance
   * - ``newton_solver_atol``
     - ``1.0e-99``
     - Newton solver absolute tolerance
   * - ``newton_solver_stol``
     - ``1.0e-12``
     - Newton solver step tolerance
   * - ``newton_solver_maxits``
     - ``10``
     - Newton solver maximum iterations
   * - ``place_randomly_in_cells``
     - ``true``
     - Randomly place SDs within grid cells
   * - ``sigma0``
     - ``0.62``
     - Kernel width parameter for distribution estimation
   * - ``distribution_grid_size``
     - ``100``
     - Grid size for computing distributions
   * - ``inactive_threshold``
     - ``0.01``
     - Multiplicity threshold for deactivating particles
   * - ``write_inactive_plt``
     - ``false``
     - Write inactive particles to plot files
   * - ``num_initializations``
     - ``1``
     - Number of initialization regions
   * - ``num_injections``
     - ``0``
     - Number of injection sources

Initialization Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

For each initialization region (indexed by ``N``), parameters are specified with prefix ``super_droplets_moisture.N``:

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Parameter
     - Default Value
     - Description
   * - ``distribution_type``
     - ``uniform``
     - Spatial distribution (``uniform``, ``bubble``)
   * - ``particles_per_cell``
     - (required)
     - Number of super-droplets per grid cell
   * - ``particle_box_lo``
     - domain lower bounds
     - Lower corner of initialization box (for ``uniform``)
   * - ``particle_box_hi``
     - domain upper bounds
     - Upper corner of initialization box (for ``uniform``)
   * - ``particle_bubble_center``
     - --
     - Center of initialization bubble (for ``bubble``)
   * - ``particle_bubble_radius``
     - --
     - Radius of initialization bubble (for ``bubble``)
   * - ``initial_number_density``
     - --
     - Physical droplet number density (m\ :sup:`-3`)
   * - ``multiplicity_type``
     - ``sampled``
     - How to assign multiplicity (``sampled``, ``constant``)
   * - ``maximum_multiplicity``
     - --
     - Maximum multiplicity value

Species and Aerosol Distribution Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each initialization region (indexed by ``N``) and each species/aerosol (``<NAME>``), distribution parameters are specified with prefix ``super_droplets_moisture.N.initial_``:

.. list-table::
   :header-rows: 1
   :widths: 40 25 35

   * - Parameter
     - Default Value
     - Description
   * - **Species Parameters**
     -
     -
   * - ``species_distribution_type_<NAME>``
     - ``mass_constant``
     - Distribution type for species <NAME>: ``mass_constant``, ``mass_exponential``, ``radius_log_normal``, ``radius_lognormal_autorange``
   * - ``species_min_mass_<NAME>``
     - --
     - Minimum mass (kg) for species <NAME>
   * - ``species_mean_mass_<NAME>``
     - --
     - Mean mass (kg) for species <NAME>
   * - ``species_max_mass_<NAME>``
     - 5 × mean_mass
     - Maximum mass (kg) for species <NAME>
   * - ``species_min_radius_<NAME>``
     - --
     - Minimum radius (m) for species <NAME>
   * - ``species_max_radius_<NAME>``
     - --
     - Maximum radius (m) for species <NAME>
   * - ``species_mean_radius_<NAME>``
     - --
     - Mean radius (m) for species <NAME> (used with log-normal distributions)
   * - ``species_std_radius_<NAME>``
     - --
     - Standard deviation of radius (m) for species <NAME> (converted to geometric std; cannot specify both std and geomstd)
   * - ``species_geomstd_radius_<NAME>``
     - --
     - Geometric standard deviation (dimensionless) of radius for species <NAME> (used with log-normal distributions)
   * - **Aerosol Parameters**
     -
     -
   * - ``aerosol_distribution_type_<NAME>``
     - ``mass_constant``
     - Distribution type for aerosol <NAME>: ``mass_constant``, ``mass_exponential``, ``radius_log_normal``, ``radius_lognormal_autorange``
   * - ``aerosol_min_mass_<NAME>``
     - --
     - Minimum mass (kg) for aerosol <NAME>
   * - ``aerosol_mean_mass_<NAME>``
     - --
     - Mean mass (kg) for aerosol <NAME>
   * - ``aerosol_max_mass_<NAME>``
     - 5 × mean_mass
     - Maximum mass (kg) for aerosol <NAME>
   * - ``aerosol_min_radius_<NAME>``
     - --
     - Minimum radius (m) for aerosol <NAME>
   * - ``aerosol_max_radius_<NAME>``
     - --
     - Maximum radius (m) for aerosol <NAME>
   * - ``aerosol_mean_radius_<NAME>``
     - --
     - Mean radius (m) for aerosol <NAME> (used with log-normal distributions)
   * - ``aerosol_std_radius_<NAME>``
     - --
     - Standard deviation of radius (m) for aerosol <NAME> (converted to geometric std; cannot specify both std and geomstd)
   * - ``aerosol_geomstd_radius_<NAME>``
     - --
     - Geometric standard deviation (dimensionless) of radius for aerosol <NAME> (used with log-normal distributions)

Diagnostic Parameters
^^^^^^^^^^^^^^^^^^^^^

Diagnostic parameters use the prefix ``super_droplets_moisture``:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Parameter
     - Default Value
     - Description
   * - ``mass_change_unconverged_log``
     - ``false``
     - Log particles with unconverged mass change
   * - ``mass_change_unconverged_log_filename``
     - ``unconverged_superdroplets.log``
     - Filename for unconverged particle log
   * - ``distribution_rmin``
     - ``1.0e-6``
     - Minimum radius (m) for binned distributions
   * - ``distribution_rmax``
     - ``5.0e-3``
     - Maximum radius (m) for binned distributions

Particle Injection
~~~~~~~~~~~~~~~~~~

The SDM implementation supports runtime injection of super-droplets from specified sources. This feature enables simulation
of scenarios such as cloud seeding, aerosol injection, or particle sources at boundaries. Multiple injection sources can be
configured, each with its own spatial distribution, timing, and particle properties.

Overview
^^^^^^^^

Particle injection creates new super-droplets at regular intervals during the simulation. Each injection source is defined
by a spatial region (uniform box or ellipsoidal bubble), temporal window (start and stop times), injection rate, and particle
attribute distributions. The injection region can optionally move with a prescribed velocity to simulate moving sources.

Injection Process
^^^^^^^^^^^^^^^^^

At each timestep, if the current simulation time falls within an injection source's active time window, new super-droplets
are created within the injection region. The number of super-droplets created is determined by:

.. math::
   N_\text{inject} = \text{rate} \times V_\text{region} \times \Delta t / \text{multiplicity}

where :math:`\text{rate}` is the injection rate (physical particles per unit volume per unit time), :math:`V_\text{region}`
is the volume of the injection region, :math:`\Delta t` is the timestep, and :math:`\text{multiplicity}` is computed from
the ``particles_per_cell`` parameter similar to initialization.

Injection Parameters
^^^^^^^^^^^^^^^^^^^^

For each injection source (indexed by ``N``), parameters are specified with prefix ``super_droplets_moisture.injection.N``:

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Parameter
     - Default Value
     - Description
   * - **Spatial Distribution**
     -
     -
   * - ``distribution_type``
     - ``uniform``
     - Spatial distribution (``uniform``, ``bubble``)
   * - ``particle_box_lo``
     - --
     - Lower corner of injection box (for ``uniform``)
   * - ``particle_box_hi``
     - --
     - Upper corner of injection box (for ``uniform``)
   * - ``particle_bubble_center``
     - --
     - Center of injection bubble (for ``bubble``)
   * - ``particle_bubble_radius``
     - --
     - Radius of injection bubble (for ``bubble``)
   * - ``domain_velocity``
     - ``[0, 0, 0]``
     - Velocity (m/s) of moving injection region
   * - **Temporal Control**
     -
     -
   * - ``t_start``
     - ``0.0``
     - Injection start time (s)
   * - ``t_stop``
     - :math:`\infty`
     - Injection stop time (s)
   * - **Injection Rate**
     -
     -
   * - ``rate``
     - (required)
     - Injection rate (physical particles m\ :sup:`-3` s\ :sup:`-1`)
   * - ``particles_per_cell``
     - (required)
     - Number of super-droplets per cell for injection
   * - **Particle Attributes**
     -
     -
   * - ``species_distribution_type_<NAME>``
     - ``mass_constant``
     - Distribution type for species <NAME>: ``mass_constant``, ``mass_exponential``, ``radius_log_normal``, ``radius_lognormal_autorange``
   * - ``aerosol_distribution_type_<NAME>``
     - ``mass_constant``
     - Distribution type for aerosol <NAME>:  ``mass_constant``, ``mass_exponential``, ``radius_log_normal``, ``radius_lognormal_autorange``
   * - (other species/aerosol parameters)
     - --
     - Same as initialization parameters

Example Configuration
^^^^^^^^^^^^^^^^^^^^^

The following example demonstrates injection configuration with three sources: two uniform box regions with opposing velocities,
and one stationary bubble region with temporal control:

.. code-block:: bash

   super_droplets_moisture.num_injections = 3

   # Moving box injection from left
   super_droplets_moisture.injection.0.distribution_type = "uniform"
   super_droplets_moisture.injection.0.particle_box_lo = 1000.0 0.0 5000.0
   super_droplets_moisture.injection.0.particle_box_hi = 1400.0 400.0 5400.0
   super_droplets_moisture.injection.0.domain_velocity = 18.0 0.0 0.0
   super_droplets_moisture.injection.0.rate = 2.0e6
   super_droplets_moisture.injection.0.particles_per_cell = 8
   super_droplets_moisture.injection.0.aerosol_distribution_type_NaCl = "mass_exponential"
   super_droplets_moisture.injection.0.aerosol_mean_mass_NaCl = 1.0e-19

   # Time-limited box injection from right
   super_droplets_moisture.injection.1.distribution_type = "uniform"
   super_droplets_moisture.injection.1.particle_box_lo = 16000.0 0.0 4000.0
   super_droplets_moisture.injection.1.particle_box_hi = 16400.0 400.0 4400.0
   super_droplets_moisture.injection.1.domain_velocity = -18.0 0.0 0.0
   super_droplets_moisture.injection.1.t_start = 200.0
   super_droplets_moisture.injection.1.t_stop = 600.0
   super_droplets_moisture.injection.1.rate = 2.0e7
   super_droplets_moisture.injection.1.particles_per_cell = 32

   # Stationary bubble injection (early times only)
   super_droplets_moisture.injection.2.distribution_type = "bubble"
   super_droplets_moisture.injection.2.particle_bubble_center = 10000.0 0.0 2000.0
   super_droplets_moisture.injection.2.particle_bubble_radius = 2000.0 400.0 2000.0
   super_droplets_moisture.injection.2.t_start = 0.0
   super_droplets_moisture.injection.2.t_stop = 200.0
   super_droplets_moisture.injection.2.rate = 5.0e6
   super_droplets_moisture.injection.2.particles_per_cell = 16

Particle Recycling
~~~~~~~~~~~~~~~~~~

Particle recycling is a mechanism to maintain a relatively constant number of super-droplets in the simulation by reactivating
deactivated particles. When super-droplets evaporate completely or fall below a multiplicity threshold, they are deactivated
rather than removed from the simulation. Recycling allows these inactive particles to be reintroduced with new attributes,
which is particularly useful for maintaining statistical convergence in long-running simulations or simulations with strong
sedimentation.

Recycling Mechanism
^^^^^^^^^^^^^^^^^^^

When ``recycle_particles`` is enabled, the model tracks the fraction of deactivated super-droplets. Once this fraction exceeds
the ``inactive_threshold`` (default 1%), the recycling process is triggered:

1. **Selection**: Inactive super-droplets are identified as candidates for recycling.

2. **Repositioning**: Selected particles are repositioned within the domain. The new positions are sampled from the original
   initialization distribution. If a recycling bounding box is specified (``recycle_xmin`` through ``recycle_zmax``), particles
   are constrained to this region; otherwise, the entire domain is used.

3. **Attribute Resampling**: Particle attributes (species masses, aerosol masses) are resampled from the original initialization
   distributions to maintain the prescribed size and composition distributions.

4. **Multiplicity Reset**: Multiplicities are reset based on the original initialization parameters to represent the appropriate
   number of physical droplets.

5. **Reactivation**: The recycled super-droplets are marked as active and reintroduced into the simulation.

This approach maintains computational efficiency by reusing particle memory rather than deallocating and reallocating particles,
while preserving the statistical properties of the particle population.

Recycling Parameters
^^^^^^^^^^^^^^^^^^^^

Recycling behavior is controlled by parameters with prefix ``super_droplets_moisture``:

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Parameter
     - Default Value
     - Description
   * - ``recycle_particles``
     - ``false``
     - Enable particle recycling
   * - ``inactive_threshold``
     - ``0.01``
     - Fraction of inactive particles (0.01 = 1%) that triggers recycling
   * - ``recycle_xmin``, ``recycle_xmax``
     - domain bounds in x
     - Extent of the recycling region in x
   * - ``recycle_ymin``, ``recycle_ymax``
     - domain bounds in y
     - Extent of the recycling region in y
   * - ``recycle_zmin``, ``recycle_zmax``
     - domain bounds in z
     - Extent of the recycling region in z

Use Cases
^^^^^^^^^

Particle recycling is particularly beneficial for:

- **Sedimentation-dominated flows**: Where particles continuously fall out of the domain and need to be replenished at the top
- **Long-time simulations**: To maintain statistical convergence without accumulating large numbers of inactive particles
- **Boundary layer simulations**: Where particles entering from one boundary can be recycled to maintain a quasi-steady state

Example Configuration
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Enable recycling when 1% of particles become inactive
   super_droplets_moisture.recycle_particles = true
   super_droplets_moisture.inactive_threshold = 0.01

   # Constrain recycled particles to upper portion of domain
   super_droplets_moisture.recycle_xmin = 0.0
   super_droplets_moisture.recycle_xmax = 20000.0
   super_droplets_moisture.recycle_ymin = 0.0
   super_droplets_moisture.recycle_ymax = 400.0
   super_droplets_moisture.recycle_zmin = 8000.0
   super_droplets_moisture.recycle_zmax = 10000.0

Material Properties and Aerosol Species
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SDM implementation supports multiple aerosol and condensate species with different physical and chemical properties.
Each species is characterized by its material properties including density, molecular weight, ionization state (for soluble
species), and thermodynamic properties for phase change calculations.

Available Species
^^^^^^^^^^^^^^^^^

The following species are currently implemented in ERF:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Species Name
     - Density (kg/m\ :sup:`3`)
     - Molecular Weight (kg/mol)
     - Ionization
     - Type
   * - ``H2O``
     - 1000.0
     - 0.01802
     - 0
     - Condensate (water)
   * - ``NaCl``
     - 2170.0
     - 0.05844
     - 2
     - Soluble aerosol (sodium chloride)
   * - ``NH42SO4``
     - 1770.0
     - 0.13214
     - 3
     - Soluble aerosol (ammonium sulfate)
   * - ``NH4HSO4``
     - 1780.0
     - 0.11511
     - 2
     - Soluble aerosol (ammonium bisulfate)
   * - ``soil``
     - 1220.0
     - --
     - 0
     - Insoluble aerosol

Species Classification
^^^^^^^^^^^^^^^^^^^^^^

**Water Species** (``H2O``): These represent the condensate phase. Water species have associated saturation vapor
pressure and latent heat properties required for phase change calculations.

.. The aliases ``water`` and ``agua`` are provided to enable multi-component testing with identical water species.

**Soluble Aerosols** (``NaCl``, ``NH42SO4``, ``NH4HSO4``): These species dissolve in water droplets and affect the equilibrium
vapor pressure through the Köhler effect (solution term). The ionization factor (van't Hoff factor) represents the number of
ions produced when the species dissociates in solution, which directly affects the strength of the solution effect on droplet
growth.

**Insoluble Aerosols** (``soil``): These species do not dissolve in water but can serve as cloud condensation nuclei (CCN) or
be incorporated into droplets. They affect droplet properties through their mass and volume but do not contribute to the
solution term in Köhler theory.

**Ice Nucleating Properties**: Both soluble and insoluble aerosol species can be marked as ice nucleating particles (INPs)
using the ``m_is_INP`` flag in the MaterialProperties class. INPs have the ability to catalyze ice formation, affecting the
freezing temperature of droplets containing them. The INAS-based freezing temperature initialization uses the total surface
area of all INP materials within a droplet to determine its freezing temperature following the Niemand et al. (2012) parameterization.

Multi-Species Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Species are specified as comma-separated lists in the input file:

.. code-block:: bash

   # Specify condensate and aerosol species
   super_droplets_moisture.species = H2O
   super_droplets_moisture.aerosols = NH4HSO4

   # For multiple aerosol species
   super_droplets_moisture.aerosols = NH4HSO4, NaCl

Each species then requires its own distribution parameters for initialization and injection:

.. code-block:: bash

   # NH4HSO4 aerosol with log-normal radius distribution
   super_droplets_moisture.0.initial_aerosol_distribution_type_NH4HSO4 = "radius_log_normal"
   super_droplets_moisture.0.initial_aerosol_mean_radius_NH4HSO4 = 30.0e-9
   super_droplets_moisture.0.initial_aerosol_std_radius_NH4HSO4 = 0.247

   # NaCl aerosol with exponential mass distribution
   super_droplets_moisture.0.initial_aerosol_distribution_type_NaCl = "mass_exponential"
   super_droplets_moisture.0.initial_aerosol_mean_mass_NaCl = 1.0e-19
   super_droplets_moisture.0.initial_aerosol_max_mass_NaCl = 1.0e-18

Impact on Microphysics
^^^^^^^^^^^^^^^^^^^^^^

The choice of aerosol species affects several microphysical processes:

**Condensation/Evaporation**: Soluble aerosols reduce the equilibrium vapor pressure over a droplet (solution effect),
promoting condensation and inhibiting evaporation. The strength of this effect scales with the ionization factor and
aerosol mass.

**Activation**: Droplets containing soluble aerosols can activate (begin growing by condensation) at lower supersaturations
compared to pure water droplets of the same size.

**Sedimentation**: Different aerosol densities affect the effective density of cloud droplets, slightly modifying terminal
fall velocities for very small droplets where aerosol mass is comparable to water mass.

Examples and Test Cases
~~~~~~~~~~~~~~~~~~~~~~~~

The ERF repository includes several test cases and examples demonstrating SDM capabilities. These cases provide templates
for configuring SDM simulations and serve as verification benchmarks.

Idealized Test Cases
^^^^^^^^^^^^^^^^^^^^^

Source terms for all of these tests are located in ``Source/Prob``.

**SDM_Bubble2D_Adv**: 2D advection test of a moist bubble with super-droplets. Tests particle advection with the flow field
and basic particle dynamics without coalescence or phase change.
Inputs located in ``Tests/test_files/SDM_Bubble2D_Adv/``.

**SDM_Bubble2D_Adv_InitSampling**: Similar to SDM_Bubble2D_Adv but demonstrates sampling-based multiplicity assignment for
improved representation of the size distribution.
Inputs located in ``Tests/test_files/SDM_Bubble2D_Adv_InitSampling/``.

**SDM_Bubble2D_Adv_wInjection**: 2D moist bubble with runtime particle injection. Demonstrates injection configuration with
moving injection domains.
Inputs located in ``Tests/test_files/SDM_Bubble2D_Adv_wInjection/``.

**SDM_Bubble2D_Adv_TfzINAS**: Tests the INAS-based freezing temperature initialization for superdroplets. Demonstrates
proper initialization of freezing temperatures based on the Niemand et al. (2012) parameterization, with INP surface area
derived from aerosol masses. Located in ``Tests/test_files/SDM_Bubble2D_Adv_TfzINAS/``.

**SDM_Box3D_Cond**: 3D box test for condensation/evaporation processes. Tests phase change physics with fixed environmental
conditions.
Inputs located in ``Tests/test_files/SDM_Box3D_Cond/``.

**SDM_Box3D_VTerm**: 3D box test for terminal velocity and sedimentation. Tests gravitational settling with various terminal
velocity formulations.
Inputs located in ``Tests/test_files/SDM_Box3D_VTerm/``.

**SDM_Box3D_Recycling**: 3D box test demonstrating particle recycling at domain boundaries. Shows how to maintain particle
population during sedimentation.
Inputs located in ``Tests/test_files/SDM_Box3D_Recycling/``.

**SDM_MultiSpecies_Bubble2D**: 2D moist bubble with multiple aerosol species. Demonstrates multi-component configuration.
Inputs located in ``Tests/test_files/SDM_MultiSpecies_Bubble2D/``.

Realistic Test Cases
^^^^^^^^^^^^^^^^^^^^^

**SDM_RICO3D**: 3D simulation of the Rain In Cumulus Over the Ocean (RICO) case, a precipitating shallow cumulus benchmark.
Tests full SDM microphysics including condensation, coalescence, and sedimentation in a realistic cloud environment.
Source files located in ``Source/Prob/ERF*RICO3D.H``.
Inputs with NH4HSO4 aerosol located in  ``Tests/test_files/SDM_RICO3D/``.

**SDM_RICO3D_InitSampling**: RICO case with sampling-based initialization for improved size distribution representation.
Source files located in ``Source/Prob/ERF*SDM_RICO3D.H``.
Inputs located in ``Tests/test_files/SDM_RICO3D_InitSampling/``.

**SDM_Congestus3D**: 3D simulation of congestus clouds, testing SDM in a deeper convective environment.
Source files located in ``Source/Prob/ERF*SDM_Congestus3D.H``.
Inputs located in ``Tests/test_files/SDM_Congestus3D/``.

Example Problems
^^^^^^^^^^^^^^^^^

**Moist Bubble with Multi-Injection**: ``Exec/RegTests/Bubble/`` contains multiple inputs files for different microphysics models,
the ``inputs_BF02_moist_bubble_SDM_multi_injections_unimodal_NaCl`` demonstrates SDM a complex injection setup with three injection sources:
two moving box regions with opposing velocities and one time-limited bubble injection. This example is useful for understanding injection
configuration and moving source regions.

**RICO DevTest**: ``Exec/CanonicalTests/RICO/`` contains multiple input files for the RICO case with different microphysics models,
including SDM configurations with various aerosol species (``inputs_SDM``).

**Temperature Source Tests**: ``Exec/RegTests/SDM_Congestus3D`` and ``Exec/RegTests/SineMassFlux``
include SDM configurations for testing particle behavior with prescribed temperature and mass flux forcing.

Verification
^^^^^^^^^^^^

Many test cases include gold files (reference solutions) in ``Tests/ERFGoldFiles/`` for automated verification. When running
with CTest, results are automatically compared against these gold files to ensure numerical consistency.

Output Files
~~~~~~~~~~~~

The super-droplet method generates several types of output files in addition to the standard AMReX plot files and particle output files.

Particle Plot Files
^^^^^^^^^^^^^^^^^^^

Super-droplet data is written to AMReX particle plot files at regular intervals (controlled by the main ERF output settings). These files contain all particle attributes including position, velocity, mass, radius, multiplicity, and species/aerosol masses. The particle files can be visualized and analyzed using standard AMReX visualization tools.

Inactive Particle Plot Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``write_inactive_plt`` is enabled, deactivated (inactive) super-droplets are written to separate plot file directories named ``deac_SD_NNNNN`` where ``NNNNN`` is the timestep number. These files allow tracking of particles that have been deactivated due to evaporation or other processes. Each output includes both a dummy Eulerian field and the particle data in a standard AMReX plot file format.

Droplet Size Distribution Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model can output droplet size distribution (DSD) data in two formats:

**Kernel-Smoothed Distributions** (``super_droplets_moisture_g_lnR_NNNNN.txt``):

These text files contain kernel-smoothed mass distributions as a function of the natural logarithm of droplet radius. The files are written at intervals specified by ``diagnostics_interval`` and contain multiple columns:

- Column 1: Droplet radius :math:`R` (m)
- Column 2: Mass-weighted distribution :math:`g_m(\ln R)` (kg/m\ :sup:`3`)
- Columns 3+: Mass distributions for individual aerosol species (if present)

The kernel-smoothed distribution is computed using:

.. math::
   g_m(\ln R) = \sum_{i=1}^{N_s} \gamma_{ker} \xi_i m_i \exp\left(-\lambda_{ker} (\ln R - \ln R_i)^2\right)

where:
- :math:`g_m(\ln R)` is the mass-weighted size distribution
- :math:`\gamma_{ker}` is a normalization factor for the kernel
- :math:`\lambda_{ker}` is the kernel width parameter related to :math:`\sigma_0`
- :math:`\xi_i` is the multiplicity of super-droplet i
- :math:`m_i` is the mass of a single droplet represented by super-droplet i
- :math:`R_i` is the radius of droplet i

where the kernel width is controlled by ``sigma0``.

**Binned Distributions** (available with ``ERF_USE_ML_UPHYS_DIAGNOSTICS``):

When compiled with machine learning diagnostics enabled, binned distributions are written to both text and plot file formats:

- ``super_droplets_moisture_binned_dsd_NNNNN.txt``: Domain-integrated binned distributions with columns for radius, mass distribution :math:`g_m(\ln R)`, and number distribution :math:`g_n(\ln R)`.

- ``super_droplets_moisture_binned_dsd_mass_NNNNN/`` and ``super_droplets_moisture_binned_dsd_number_NNNNN/``: AMReX plot file directories containing cell-wise binned mass and number distributions. These files allow spatial analysis of the droplet size distribution.

The binned distributions divide the radius range from ``distribution_rmin`` to ``distribution_rmax`` into ``distribution_grid_size`` bins with logarithmic spacing. The distributions are normalized by cell volume and bin width :math:`d\ln R`.

Unconverged Particle Log
^^^^^^^^^^^^^^^^^^^^^^^^^

When ``mass_change_unconverged_log`` is enabled (CPU only), particles that fail to converge during the phase change time integration are logged to the file specified by ``mass_change_unconverged_log_filename`` (default: ``unconverged_superdroplets.log``). This file helps diagnose numerical issues with the condensation/evaporation solver and can be used to tune solver tolerances.

References
~~~~~~~~~~

- Shima, S., K. Kusano, A. Kawano, T. Sugiyama, and S. Kawahara, 2009: The super-droplet method for the numerical
  simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a
  non-hydrostatic model. Q. J. R. Meteorol. Soc., 135: 1307-1320.

- Shima, S., Y. Sato, A. Hashimoto, and R. Misumi, 2020: Predicting the morphology of ice particles in deep convection
  using the super-droplet method: development and evaluation of SCALE-SDM 0.2.5-2.2.0, -2.2.1, and -2.2.2. Geosci. Model
  Dev., 13: 4107-4157, https://doi.org/10.5194/gmd-13-4107-2020.

- Böhm, J. P., 1989: A general equation for the terminal fall speed of solid hydrometeors. J. Atmos. Sci., 46: 2419-2427.

- Böhm, J. P., 1992: A general hydrodynamic theory for mixed-phase microphysics. Part I: Drag and fall speed of
  hydrometeors. Atmos. Res., 27: 253-274.

- Böhm, J. P., 1999: Revision and clarification of "A general hydrodynamic theory for mixed-phase microphysics".
  Atmos. Res., 52: 167-176.

- Beard, K. V., 1976: Terminal velocity and shape of cloud and precipitation drops aloft. J. Atmos. Sci., 33:
  851-864, https://doi.org/10.1175/1520-0469(1976)033<0851:TVASOC>2.0.CO;2.

- Gunn, R., and G. D. Kinzer, 1949: The terminal velocity of fall for water droplets in stagnant air. J. Meteorol.,
  6: 243-248, https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2.

- Niemand, M., O. Möhler, B. Vogel, H. Vogel, C. Hoose, P. Connolly, H. Klein, H. Bingemer, P. DeMott, J. Skrotzki, and T. Leisner, 2012:
  A particle-surface-area-based parameterization of immersion freezing on desert dust particles. J. Atmos. Sci., 69: 3077-3092,
  https://doi.org/10.1175/JAS-D-11-0249.1.

- Seifert, A., and K. D. Beheng, 2006: A two-moment cloud microphysics parameterization for mixed-phase clouds.
  Part 1: Model description. Meteorol. Atmos. Phys., 92: 45-66, https://doi.org/10.1007/s00703-005-0112-4.

- Mitra, S. K., O. Vohl, M. Ahr, and H. R. Pruppacher, 1990: A wind tunnel and theoretical study of the melting
  behavior of atmospheric ice particles. IV: Experiment and theory for snow flakes. J. Atmos. Sci., 47: 584-591,
  https://doi.org/10.1175/1520-0469(1990)047<0584:AWTATS>2.0.CO;2.

- Frick, C., A. Seifert, and H. Wernli, 2013: A bulk parametrization of melting snowflakes with explicit liquid
  water fraction for the COSMO model. Geosci. Model Dev., 6: 1925-1939, https://doi.org/10.5194/gmd-6-1925-2013.
