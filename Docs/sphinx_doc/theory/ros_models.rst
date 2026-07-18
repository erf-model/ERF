
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

This is a physics-based model derived from flame radiation geometry and buoyancy effects. It computes rate of spread as a function of wind, slope, and fuel properties.

.. math::

   R = A \cdot (1 + \sin\alpha - \cos\alpha)

   \tan\alpha = \frac{U}{v_b} + \tan\theta

   v_b = \sqrt{\frac{g \cdot \delta_m \cdot (T_f - T_a)}{T_a}}

where :math:`A` [m/s] is the amplitude coefficient computed from fuel properties, :math:`v_b` [m/s] is the buoyancy velocity, :math:`U` [m/s] is the wind speed in the fire spread direction, :math:`\theta` [radians] is the terrain slope, :math:`\delta_m` [m] is the fuel bed depth, :math:`g` is gravitational acceleration, :math:`T_f` is flame temperature, and :math:`T_a` is ambient temperature.

The model returns zero ROS when fuel bed depth :math:`\delta_m < 0.01` m.

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

**Reference:**

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model for wildland fires. *Combustion and Flame*, 156(12), 2217-2230.

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

- The Balbi model uses a simplified :math:`A_{\text{coeff}}` derivation based on reaction intensity analogs rather than the full Balbi (2009) flame geometry system.

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
