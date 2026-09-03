
.. role:: cpp(code)
   :language: c++

.. _sec:ROS_Models:

Rate-of-Spread Models
=====================

Overview
--------

ERF-Fire provides four selectable rate-of-spread models via the :cpp:`erf.fire.ros_model` parameter. The default is Rothermel (1972), and all other models are backward-compatible drop-in replacements that consume the same effective midflame wind and write to the same ``fire_ros`` MultiFab.

Sub-phase A adds per-fuel wind extraction height following WRF-SFIRE :cpp:`fcwh` convention. This is controlled by the :cpp:`erf.fire.use_per_fuel_wind_ht` parameter and is backward-compatible with the global wind reference height approach.

Rothermel (1972) - Default
---------------------------

The default single-class empirical model implemented in Phases 1-12. This model computes the rate of fire spread based on fuel characteristics, moisture content, wind speed, and slope. 

**Key parameter:** :cpp:`ros_model = "rothermel"`

**Reference:** Rothermel, R.C. (1972). A Mathematical Model for Predicting Fire Spread in Wildland Fuels. USDA Forest Service Research Paper INT-115.

MacArthur (1966) Australian Formula
------------------------------------

**Key parameter:** :cpp:`ros_model = "macarthur"`

This model implements the MacArthur (1966) Mark 5 Forest Fire Danger Meter formula used in WRF-SFIRE as ``ibeh=0`` (Clark et al. 2004).

.. math::

   R = R_b \cdot \exp(0.8424 \cdot \max(U,\,0))

where :math:`R_b = 0.18\ \text{m/s}` is the backing (no-wind) rate of spread and :math:`U` is the effective midflame wind speed [m/s].

The model is appropriate for Australian grassland and open forest fuels. It is not calibrated for North American FBFM13 fuel models.

**References:**

- McArthur, A.G. (1966). Weather and Grassland Fire Behaviour. Leaflet 100, Forestry and Timber Bureau, Canberra.
- Clark, T.L., Coen, J.L. &amp; Latham, D. (2004). Description of a coupled atmosphere-fire model. *International Journal of Wildland Fire*, 13, 49-63.

Balbi (2009) Physical Model
----------------------------

**Key parameter:** :cpp:`ros_model = "balbi"`

This is a physics-based model derived from the radiant energy flux received by unburned fuel from a tilted flame front. It computes rate of spread as a function of wind, slope, fuel properties and fuel moisture.

.. math::

   R = A \cdot (1 + \sin\alpha - \cos\alpha)

   A = \frac{\chi \, \sigma_m \, \delta_m}{2 \, \tau_0 \, B^*}

   \chi = \frac{r_{00} \sigma_m}{1 + r_{00} \sigma_m}

   B^* = \frac{C_{pf} (T_i - T_a) + M_f \Lambda}{h}

   \tan\alpha = \frac{U}{v_b} + \tan\theta

   v_b = \sqrt{\frac{g \cdot \delta_m \cdot (T_f - T_a)}{T_a}}

where :math:`A` [m/s] is the amplitude coefficient, :math:`\chi` [-] is the fraction of combustion heat radiated forward, :math:`B^*` [-] is the dimensionless energy needed to ignite the fuel, :math:`v_b` [m/s] is the buoyancy velocity, :math:`U` [m/s] is the effective midflame wind speed, :math:`\theta` [radians] is the terrain slope, :math:`\sigma_m` [1/m] is the surface-area-to-volume ratio of the 1-hr dead fuel, :math:`\delta_m` [m] is the fuel bed depth, :math:`h` [J/kg] is the heat of combustion, :math:`M_f` [-] is the fuel moisture content, :math:`\Lambda` [J/kg] is the latent heat of vaporisation, :math:`g` is gravitational acceleration, :math:`T_f` is flame temperature, :math:`T_i` is ignition temperature and :math:`T_a` is ambient temperature.

Fuel properties enter :math:`A` through :math:`\sigma_m`, :math:`\delta_m` and :math:`h`, so ROS varies by fuel model. When a spatial fuel map is loaded, per-cell coefficients are taken from a lookup table built for the 13 Anderson fuel models; fuel code 0 (non-burnable) gives zero spread. When :cpp:`erf.fire.moisture_dynamic = true`, :math:`B^*` and hence :math:`A` are rebuilt each step from the domain-average 1-hr moisture.

The model returns zero ROS when fuel bed depth :math:`\delta_m < 0.01` m.

Cells whose wind exceeds 25 m/s, or holds a non-finite value, fall back to the domain-mean wind before the ROS is evaluated.

Balbi (2020) convective-radiative form
``````````````````````````````````````

Setting :cpp:`erf.fire.balbi.formulation = "2020"` selects the convective-radiative model of Balbi et al. (2020) instead. The rate of spread is the root of

.. math::

   R = R_b + R_c + R_r

   R_b = \min\!\left(\frac{S}{\pi}, 1\right) \frac{B T^4}{\beta \rho_v q}

   R_c = \frac{s \Delta H}{q \tau_0} \min\!\left(\delta_m, \frac{2\pi}{s\beta}\right)
         \left[ \frac{\delta_m}{2\delta_m + H}\tan\alpha + \frac{U e^{-K_1 \sqrt{\beta} R}}{u_0} \right]

   R_r = A R \frac{1 + \sin\gamma - \cos\gamma}{1 + R\cos\gamma/(s r_{00})}

with the ignition energy, radiant factor, radiant fraction, flame temperature, vertical gas velocity, flame tilt and flame height

.. math::

   q = C_{pf}(T_i - T_a) + M_f\left(\Lambda + C_{pw}(T_{vap} - T_a)\right)

   A = \min\!\left(\frac{S}{2\pi}, 1\right)\frac{\chi_0 \Delta H}{4q}
   \qquad
   \chi = \frac{\chi_0}{1 + R\cos\gamma/(s r_{00})}

   T = T_a + \frac{\Delta H (1 - \chi)}{C_{pa}(s_t + 1)}
   \qquad
   u_0 = \frac{2(s_t+1)}{\tau_0}\frac{T}{T_a}\frac{\rho_v}{\rho_a}\min(S, 2\pi)

   \gamma = \arctan\left(\tan\alpha + \frac{U}{u_0}\right)
   \qquad
   H = \frac{u_0^2}{g\left(T/T_a - 1\right)}

where :math:`\beta = w/(\delta_m \rho_v)` is the packing ratio, :math:`S = s \beta \delta_m` is twice the leaf area index, :math:`w` is the 1-hr dead fuel load, :math:`\rho_v` the fuel particle density, :math:`B` the Stefan-Boltzmann constant, :math:`K_1` a drag coefficient and :math:`s_t` the air-to-pyrolysis-gas mass ratio.

Two properties distinguish it from the 2009 form: the radiative base term :math:`R_b` gives a nonzero rate of spread with no wind and no slope, and the wind response does not saturate.

The equation is solved per cell by bisection on :math:`f(R) = R_b + R_c + R_r - R`. The reference implementations iterate by plain substitution, which oscillates rather than converging when the radiative coefficient :math:`A` is well below one — the low-fuel-load case of Fig. 3 in the paper. Since :math:`R_c` decays with :math:`R` and :math:`R_r` saturates, :math:`f` is positive at :math:`R = 0` and negative at the 30 m/s cap, so a bracketed solve always returns a root in a fixed number of steps. :cpp:`erf.fire.balbi.max_iter` and :cpp:`erf.fire.balbi.tol` control it.

Note that :cpp:`erf.fire.balbi.r_00` defaults to :math:`2.5\times10^{-5}` under this formulation, the value fitted by Balbi et al. (2020), rather than the :math:`2.5\times10^{-4}` m radiation length scale of the 2009 form. Setting it explicitly overrides the formulation-dependent default.

Optional couplings
``````````````````

Three couplings apply to both formulations and are off by default.

:cpp:`erf.fire.balbi.directional = true` evaluates the ROS along the local front normal :math:`\hat{n} = \nabla\phi/|\nabla\phi|` rather than from the wind and slope magnitudes:

.. math::

   \tan\gamma = \frac{\mathbf{U}\cdot\hat{n}}{v_b} + \nabla z \cdot \hat{n}

The ROS field is rebuilt inside every Runge-Kutta stage, so head, flank and backing spread come from the model rather than from an imposed ellipse. This affects the level-set path only (:cpp:`erf.fire.propagation_method = "levelset"`); the FARSITE path takes its directionality from the ellipse. Under the 2009 form the backing ROS is exactly zero because that form has no base spread, so the pairing is most useful with the 2020 form. :cpp:`fire_ros` continues to hold the isotropic head-fire ROS and sets the level-set CFL.

:cpp:`erf.fire.balbi.use_surface_temp = true` takes the ambient temperature per cell from ``fire_surface_temp`` instead of the fixed :cpp:`erf.fire.balbi.T_a`. The coefficients are then rebuilt inside the kernel, since the ambient temperature enters the ignition energy, the flame temperature and the buoyancy velocity.

:cpp:`erf.fire.balbi.heat_flux_coupling = true` augments the vertical velocity scale with the fire-induced component

.. math::

   v_{b,Q} = k_{up}\sqrt{\frac{g Q H_{ref}}{\rho_a C_{pa} T_a}},
   \qquad
   v_{b,eff} = \sqrt{v_b^2 + v_{b,Q}^2}

applied to :math:`v_b` in the 2009 form and to :math:`u_0` in the 2020 form. A stronger plume stands the flame more upright and therefore reduces the forward spread. ``fire_heat_flux`` is filled at the end of the fire step, so this feedback lags the ROS by one fire step.

:cpp:`erf.fire.balbi.flame_temp_from_combustion = true` derives the 2009 mean flame temperature from the combustion energy, :math:`T_f = T_a + \Delta H (1 - \chi_0)/(C_{pa}(s_t+1))`, instead of the fixed :cpp:`erf.fire.balbi.T_f`. The 2020 form always computes its own flame temperature.

**Balbi model parameters:**

.. list-table::
   :widths: 30 15 15 40
   :header-rows: 1

   * - ParmParse Key
     - Default Value
     - Units
     - Description
   * - :cpp:`erf.fire.balbi.T_a`
     - 300.0
     - K
     - Ambient temperature
   * - :cpp:`erf.fire.balbi.T_f`
     - 1000.0
     - K
     - Mean flame temperature
   * - :cpp:`erf.fire.balbi.T_i`
     - 600.0
     - K
     - Ignition temperature
   * - :cpp:`erf.fire.balbi.delta_H`
     - 2.26e6
     - J/kg
     - Latent heat of vaporisation
   * - :cpp:`erf.fire.balbi.C_pf`
     - 1800.0
     - J/(kg·K)
     - Specific heat of dry fuel
   * - :cpp:`erf.fire.balbi.r_00`
     - 2.5e-4
     - m
     - Radiation length scale
   * - :cpp:`erf.fire.balbi.tau_0`
     - 75591.0
     - s/m
     - Residence time coefficient
   * - :cpp:`erf.fire.balbi.formulation`
     - "2009"
     - -
     - Formulation: "2009" or "2020"
   * - :cpp:`erf.fire.balbi.chi_0`
     - 0.3
     - -
     - Radiant fraction (2020)
   * - :cpp:`erf.fire.balbi.C_pa`
     - 1150.0
     - J/(kg·K)
     - Specific heat of air (2020)
   * - :cpp:`erf.fire.balbi.C_pw`
     - 4180.0
     - J/(kg·K)
     - Specific heat of water (2020)
   * - :cpp:`erf.fire.balbi.T_vap`
     - 373.0
     - K
     - Water vaporisation temperature (2020)
   * - :cpp:`erf.fire.balbi.K1`
     - 130.0
     - s/m
     - Drag coefficient (2020)
   * - :cpp:`erf.fire.balbi.st`
     - 17.0
     - -
     - Air-to-pyrolysis-gas mass ratio (2020)
   * - :cpp:`erf.fire.balbi.rho_a`
     - 1.2
     - kg/m³
     - Air density (2020)
   * - :cpp:`erf.fire.balbi.sigma_B`
     - 5.6e-8
     - W/(m²·K⁴)
     - Stefan-Boltzmann constant (2020)
   * - :cpp:`erf.fire.balbi.max_iter`
     - 40
     - -
     - Root-solve iteration cap (2020)
   * - :cpp:`erf.fire.balbi.tol`
     - 1.0e-4
     - m/s
     - Root-solve bracket width (2020)
   * - :cpp:`erf.fire.balbi.directional`
     - false
     - -
     - Direction-dependent ROS on the level-set path
   * - :cpp:`erf.fire.balbi.use_surface_temp`
     - false
     - -
     - Per-cell ambient temperature from fire_surface_temp
   * - :cpp:`erf.fire.balbi.flame_temp_from_combustion`
     - false
     - -
     - Derive T_f from the combustion energy (2009)
   * - :cpp:`erf.fire.balbi.heat_flux_coupling`
     - false
     - -
     - Augment the buoyancy velocity with the surface heat flux
   * - :cpp:`erf.fire.balbi.k_upward`
     - 1.0
     - -
     - Scale on the heat-flux induced velocity
   * - :cpp:`erf.fire.balbi.hf_ref_height`
     - 10.0
     - m
     - Reference height for the induced velocity

**Reference:**

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model for wildland fires. *Combustion and Flame*, 156(12), 2217-2230.

- Balbi, J.-H., Chatelon, F.-J., Morvan, D., Rossi, J.-L., Marcelli, T. &amp; Morandini, F. (2020). A convective-radiative propagation model for wildland fires. *International Journal of Wildland Fire*, 29(8), 723-738.

Cheney-Gould (1998) Grassland Model
------------------------------------

**Key parameter:** :cpp:`ros_model = "cheney_gould"`

This model is calibrated for Australian open grassland fuels and is based on empirical relations from CSIRO research. The forward rate of spread is computed as:

.. math::

   R = R_b \cdot (1 + w_f) \cdot f_m

   w_f = 0.15 \cdot U \cdot (\text{curing} + 0.2)

   f_m = \frac{20}{M + 1}

where :math:`R_b` is the backing rate [m/s], :math:`U` is wind speed [m/s], :math:`M` is dead fine fuel moisture [%], and curing ∈ [0, 1].

The current implementation uses domain-average moisture and curing values inside the GPU kernel (see lines 218-220 of ``ERF_CheneyGouldModel.H``). Per-cell moisture from the Phase 4 ODE system is not yet coupled.

**Not appropriate for forest or shrub fuels (FM4-FM13).**

**Cheney-Gould model parameters:**

.. list-table::
   :widths: 30 15 15 40
   :header-rows: 1

   * - ParmParse Key
     - Default Value
     - Units
     - Description
   * - :cpp:`erf.fire.cheney_gould.moisture`
     - 10.0
     - %
     - Dead fine fuel moisture content
   * - :cpp:`erf.fire.cheney_gould.curing`
     - 1.0
     - [0-1]
     - Degree of curing (1 = fully cured)

**Reference:**

- Cheney, N.P., Gould, J.S. &amp; Catchpole, W.R. (1998). Prediction of fire spread in grasslands. *International Journal of Wildland Fire*, 8(1), 1-13.

Per-Fuel Wind Height (Sub-phase A)
----------------------------------

Enabling :cpp:`erf.fire.use_per_fuel_wind_ht = true` causes wind extraction to use a per-fuel-category height following WRF-SFIRE :cpp:`fcwh` convention.

WRF-SFIRE defaults are 6.096 m for all 13 Anderson fuel models, which is identical to the :cpp:`wind_ref_ht` default of 6.1 m. Enabling this flag has no practical effect unless the :cpp:`fcwh` table is modified.

**Fuel model roughness lengths (fcz0):**

.. list-table::
   :widths: 10 40 15
   :header-rows: 1

   * - Fuel Model
     - Name
     - fcz0 [m]
   * - FM1
     - Short Grass
     - 0.0396
   * - FM2
     - Timber Grass and Understory
     - 0.0396
   * - FM3
     - Tall Grass
     - 0.100
   * - FM4
     - Chaparral
     - 0.2378
   * - FM5
     - Timber Litter
     - 0.0793
   * - FM6
     - Logging Slash and Blowdown
     - 0.0991
   * - FM7
     - Timber Litter and Understory
     - 0.0991
   * - FM8
     - Closed Timber Litter
     - 0.0079
   * - FM9
     - Hardwood Litter
     - 0.0079
   * - FM10
     - Timber Litter and Grass
     - 0.0396
   * - FM11
     - Intermediate Fuel Load Timber Litter
     - 0.0396
   * - FM12
     - High Load Conifer Litter
     - 0.0911
   * - FM13
     - Heavy Logging Slash
     - 0.1188

The parameter :cpp:`erf.fire.waf_fcz0_scale` (default 1.0) scales all :cpp:`fcz0` values uniformly.

**References:**

- WRF-SFIRE :cpp:`module_fr_sfire_phys.F`
- Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling System. USDA Forest Service General Technical Report INT-194.

Model Selection Guide
---------------------

.. list-table::
   :widths: 30 20 50
   :header-rows: 1

   * - Fuel Context
     - Recommended :cpp:`ros_model`
     - Rationale
   * - North American FBFM13 fuels
     - :cpp:`rothermel`
     - Calibrated for US fuel types
   * - Australian open grassland
     - :cpp:`cheney_gould`
     - Empirically validated for Australian grasslands
   * - Mediterranean shrub or European fuels
     - :cpp:`balbi`
     - Physics-based, not region-specific
   * - Australian sclerophyll forest
     - :cpp:`macarthur`
     - Mark 5 calibration for Australian forest

Input Parameters
----------------

.. list-table::
   :widths: 35 15 15 35
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - :cpp:`erf.fire.ros_model`
     - String
     - "rothermel"
     - ROS model selection ("rothermel", "macarthur", "balbi", "cheney_gould")
   * - :cpp:`erf.fire.use_per_fuel_wind_ht`
     - Boolean
     - false
     - Enable per-fuel wind extraction height (WRF-SFIRE :cpp:`fcwh`)
   * - :cpp:`erf.fire.waf_fcz0_scale`
     - Real
     - 1.0
     - Uniform scaling factor for all fcz0 roughness lengths
   * - :cpp:`erf.fire.balbi.T_a`
     - Real
     - 300.0
     - Balbi: ambient temperature [K]
   * - :cpp:`erf.fire.balbi.T_f`
     - Real
     - 1000.0
     - Balbi: mean flame temperature [K]
   * - :cpp:`erf.fire.balbi.T_i`
     - Real
     - 600.0
     - Balbi: ignition temperature [K]
   * - :cpp:`erf.fire.balbi.delta_H`
     - Real
     - 2.26e6
     - Balbi: latent heat of vaporisation [J/kg]
   * - :cpp:`erf.fire.balbi.C_pf`
     - Real
     - 1800.0
     - Balbi: specific heat of dry fuel [J/(kg·K)]
   * - :cpp:`erf.fire.balbi.r_00`
     - Real
     - 2.5e-4
     - Balbi: radiation length scale [m]
   * - :cpp:`erf.fire.balbi.tau_0`
     - Real
     - 75591.0
     - Balbi: residence time coefficient [s/m]
   * - :cpp:`erf.fire.cheney_gould.moisture`
     - Real
     - 10.0
     - Cheney-Gould: dead fine fuel moisture [%]
   * - :cpp:`erf.fire.cheney_gould.curing`
     - Real
     - 1.0
     - Cheney-Gould: degree of curing [0-1]

Limitations
-----------

- The MacArthur formula includes no slope term; slope influence enters only through terrain wind corrections applied before ROS computation.

- The default Balbi formulation is the steady explicit form of Balbi (2009): it carries no no-wind radiative base spread (:math:`R \to 0` as wind and slope vanish) and saturates at :math:`2A` for large flame tilt. Select :cpp:`erf.fire.balbi.formulation = "2020"` for the convective-radiative form, which has neither limitation.

- Without :cpp:`erf.fire.balbi.directional`, the Balbi flame tilt uses the wind magnitude and the slope magnitude :math:`|\nabla z|`, so the ROS field is isotropic and downslope spread is enhanced in the same way as upslope spread.

- The Balbi (2020) form uses the 1-hr dead fuel load and the fuel particle density of the selected fuel model. Live and coarse dead fuels do not enter it.

- The Cheney-Gould kernel uses placeholder domain-average moisture and curing values (10.0 % and 1.0 respectively) inside the GPU kernel. Use :cpp:`erf.fire.cheney_gould.moisture` and :cpp:`erf.fire.cheney_gould.curing` for domain-averaged values. Per-cell moisture from the Phase 4 ODE system is not yet passed into the kernel.

- Per-fuel wind height (:cpp:`use_per_fuel_wind_ht = true`) uses WRF-SFIRE default :math:`\text{fcwh} = 6.096` m for all 13 Anderson fuel models, which produces the same result as :cpp:`wind_ref_ht = 6.096` m unless the table is customised. Custom :cpp:`fcwh` values are not yet exposed through ParmParse.

- All four ROS models use the same effective midflame wind (:cpp:`fire_wind_eff`) after global WAF and terrain corrections. The WAF is derived from the global :cpp:`fuel_model_id` fuel bed depth, not from per-cell fuel bed depth when a spatial fuel map is loaded.

References
----------

- Rothermel, R.C. (1972). A Mathematical Model for Predicting Fire Spread in Wildland Fuels. USDA Forest Service Research Paper INT-115.

- McArthur, A.G. (1966). Weather and Grassland Fire Behaviour. Leaflet 100, Forestry and Timber Bureau, Canberra.

- Clark, T.L., Coen, J.L. &amp; Latham, D. (2004). Description of a coupled atmosphere-fire model. *International Journal of Wildland Fire*, 13, 49-63.

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model for wildland fires. *Combustion and Flame*, 156(12), 2217-2230.

- Cheney, N.P., Gould, J.S. &amp; Catchpole, W.R. (1998). Prediction of fire spread in grasslands. *International Journal of Wildland Fire*, 8(1), 1-13.

- Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling System. USDA Forest Service General Technical Report INT-194.

- Mandel, J. et al. (2011). A wildland fire model with data assimilation. *Mathematics and Computers in Simulation*, 79(3), 584-606.
