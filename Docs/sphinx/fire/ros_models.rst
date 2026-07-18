.. _ros_models:

Rate-of-Spread Models (Phase 13)
=================================

Overview
--------

ERF-Fire supports four selectable rate-of-spread (ROS) models, each designed for different 
fuel types, regions, and physical approaches. The default model is Rothermel (1972), unchanged 
from Phases 1–12. Alternative models—MacArthur, Balbi, and Cheney-Gould—are provided for 
Australian conditions, physics-based prediction, and grassland fuels respectively. All models 
consume the same `fire_wind_eff` (effective midflame wind after Wind Adjustment Factor and 
terrain corrections) and write to the same `fire_ros` MultiFab, ensuring that downstream 
fire propagation and heat flux calculations are unaffected by model choice.


Rothermel (1972) — Default
----------------------------

The default ROS model unchanged from Phases 1–12, based on the quasi-steady-state fire 
spread model of Rothermel (1972). This model is widely used in North American fire behavior 
prediction systems and is calibrated for Anderson FBFM13 fuel models.

**ParmParse key:**

.. code-block:: text

    erf.fire.ros_model = "rothermel"

**Reference:**

- Rothermel, R.C. (1972). *A Mathematical Model for Predicting Fire Spread in Wildland Fuels.* 
  USDA Forest Service Research Paper INT-115.


MacArthur (1966) Australian Formula
------------------------------------

Empirical fire spread model calibrated for Australian open forest and grassland fuels. 
Used in WRF-SFIRE as `ibeh=0` (Clark et al. 2004). This model is simpler than Rothermel 
and does not include slope effects directly; slope influence enters only through terrain 
wind corrections applied before ROS computation.

**Formula:**

.. math::

    R = R_b \cdot \exp(0.8424 \cdot \max(U, 0))

where :math:`R_b = 0.18` m/s is the no-wind backing rate and :math:`U` [m/s] is the 
effective midflame wind speed in the fire spread direction. Negative wind speeds (opposing 
the fire spread direction) are clamped to zero, so the minimum ROS is the backing rate.

**Characteristics:**

- Wind effect follows exponential growth (not linear)
- Appropriate for Australian sclerophyll forest and grassland fuels
- Less appropriate for North American fuel models (FM1–FM13)
- No explicit slope or moisture dependence

**ParmParse key:**

.. code-block:: text

    erf.fire.ros_model = "macarthur"

**References:**

- McArthur, A.G. (1966). *Weather and Grassland Fire Behaviour.* Leaflet 100, Forestry and 
  Timber Bureau, Canberra.
- Clark, T.L., Coen, J.L. & Latham, D. (2004). Description of a coupled atmosphere-fire model. 
  *International Journal of Wildland Fire*, 13, 49–63.


Balbi (2009) Physical Model
----------------------------

Physics-based fire spread model derived from flame radiation geometry rather than empirical 
fitting. The model couples wind speed, terrain slope, and fuel properties through a geometric 
angle representing the fire flame tilt and reaction intensity.

**Formulae:**

.. math::

    R = A \cdot (1 + \sin\alpha - \cos\alpha)

    \tan\alpha = \frac{U}{v_b} + \tan\theta

    v_b = \sqrt{\frac{g \cdot \delta_m \cdot (T_f - T_a)}{T_a}}

where:

- :math:`R` [m/s] is the rate of spread
- :math:`A` [m/s] is the fire intensity coefficient, derived from fuel load and thermal properties
- :math:`\alpha` [radians] is the flame tilt angle in the fire spread direction
- :math:`U` [m/s] is the effective midflame wind speed
- :math:`\theta` [radians] is the terrain slope in the fire spread direction
- :math:`v_b` [m/s] is the buoyancy velocity
- :math:`g = 9.81` m/s² is gravitational acceleration
- :math:`\delta_m` [m] is the fuel bed depth
- :math:`T_f` [K] is the mean flame temperature
- :math:`T_a` [K] is the ambient temperature

**Balbi Model Parameters:**

.. list-table:: Balbi (2009) Thermal Parameters
   :widths: 20 30 15 35
   :header-rows: 1

   * - Symbol
     - ParmParse key
     - Default value
     - Description
   * - :math:`T_a`
     - erf.fire.balbi.T_a
     - 300.0 K
     - Ambient temperature
   * - :math:`T_f`
     - erf.fire.balbi.T_f
     - 1000.0 K
     - Mean flame temperature
   * - :math:`T_i`
     - erf.fire.balbi.T_i
     - 600.0 K
     - Ignition temperature
   * - :math:`\Delta H`
     - erf.fire.balbi.delta_H
     - 2.26e6 J/kg
     - Latent heat of vaporisation
   * - :math:`C_{pf}`
     - erf.fire.balbi.C_pf
     - 1800.0 J/(kg·K)
     - Specific heat of dry fuel
   * - :math:`r_{00}`
     - erf.fire.balbi.r_00
     - 2.5e-4 m
     - Radiation length scale
   * - :math:`\tau_0`
     - erf.fire.balbi.tau_0
     - 75591.0 s/m
     - Residence time coefficient

**Limitations:**

- Returns zero ROS when fuel bed depth :math:`\delta_m < 0.01` m (threshold hard-coded)
- Uses simplified :math:`A` coefficient computation (not the full flame geometry derivation from Balbi 2009)
- Small buoyancy velocities (:math:`v_b < 10^{-6}` m/s) are floored to prevent division by zero

**ParmParse key:**

.. code-block:: text

    erf.fire.ros_model = "balbi"

**References:**

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model 
  for wildland fires. *Combustion and Flame*, 156(12), 2217–2230.


Cheney-Gould (1998) Grassland Model
------------------------------------

Empirical fire spread model calibrated specifically for Australian open grassland fuels. 
Provides better agreement than Rothermel for short-grass fires (FM1–FM3) under typical 
grassland conditions. The model is **not appropriate** for forest or shrub fuels (FM4–FM13).

**Formula:**

.. math::

    \text{wind\_factor} = 0.15 \cdot U \cdot (\text{curing} + 0.2)

    \text{moisture\_factor} = \frac{20}{M_d + 1}

    R = R_b \cdot (1 + \text{wind\_factor}) \cdot \text{moisture\_factor}

where:

- :math:`R` [m/s] is the rate of spread
- :math:`R_b` [m/s] is the backing (no-wind) ROS
- :math:`U` [m/s] is the effective wind speed
- :math:`\text{curing}` [0–1] is the degree of grass curing (1 = fully cured, 0 = fresh)
- :math:`M_d` [%] is the dead fine fuel moisture content
- The moisture factor gives inverse (i.e., wetter fuel → lower ROS)

**Cheney-Gould Parameters:**

.. list-table:: Cheney-Gould (1998) Model Parameters
   :widths: 25 30 15 30
   :header-rows: 1

   * - Parameter
     - ParmParse key
     - Default value
     - Description
   * - Moisture
     - erf.fire.cheney_gould.moisture
     - 10.0 %
     - Dead fine fuel moisture (placeholder: domain-average)
   * - Curing
     - erf.fire.cheney_gould.curing
     - 1.0
     - Degree of curing [0–1]

**Limitations:**

- Calibrated for open grassland; inappropriate for forest or shrub fuels
- Current implementation uses placeholder domain-average moisture and curing inside the GPU 
  kernel (see lines 218–220 of ERF_CheneyGouldModel.H)
- Per-cell moisture from the Phase 4 ODE system is not yet passed into the ROS kernel
- No explicit slope effects (slope influence enters only via WAF and terrain corrections)

**ParmParse key:**

.. code-block:: text

    erf.fire.ros_model = "cheney_gould"

**References:**

- Cheney, N.P., Gould, J.S. &amp; Catchpole, W.R. (1998). Prediction of fire spread in 
  grasslands. *International Journal of Wildland Fire*, 8(1), 1–13.


Per-Fuel Wind Height (Phase 13A Sub-phase)
-------------------------------------------

Per-fuel wind height extraction (Phase 13A) enables the fire module to extract atmospheric 
wind at fuel-model-specific heights rather than the global reference height. This follows 
the WRF-SFIRE convention of per-fuel extraction heights (`fcwh`).

**Enabling Per-Fuel Wind Height:**

.. code-block:: text

    erf.fire.use_per_fuel_wind_ht = true  (default: false)

**Default fcwh Values:**

WRF-SFIRE uses a uniform default of 6.096 m (20 feet) for all 13 Anderson fuel models, 
which is consistent with the BEHAVE default midflame height. This value matches the global 
`wind_ref_ht` default of 6.1 m, so **enabling this flag has no practical effect unless 
custom `fcwh` values are provided** (not yet exposed through ParmParse).

**Per-Fuel Roughness Length (fcz0) Values:**

The `fcz0` table specifies the surface roughness length [m] for each fuel model, used in 
the log-profile Wind Adjustment Factor interpolation:

.. list-table:: Per-Fuel Roughness Length (fcz0) for Anderson FBFM13
   :widths: 10 15 10 15 10 15 10 15
   :header-rows: 1

   * - FM
     - fcz0 [m]
     - FM
     - fcz0 [m]
     - FM
     - fcz0 [m]
     - FM
     - fcz0 [m]
   * - 1
     - 0.0396
     - 4
     - 0.2378
     - 7
     - 0.0991
     - 10
     - 0.0396
   * - 2
     - 0.0396
     - 5
     - 0.0793
     - 8
     - 0.0079
     - 11
     - 0.0396
   * - 3
     - 0.1000
     - 6
     - 0.0991
     - 9
     - 0.0079
     - 12
     - 0.0911
   * - 13
     - 0.1188
     -
     -
     -
     -
     -
     -

Note: FM4 (chaparral) has the highest roughness (0.2378 m), while FM8 and FM9 (logging slash) 
have the lowest (0.0079 m each).

**Wind Adjustment Factor (WAF) Scaling:**

A scaling factor can be applied to the `fcz0` values during WAF computation:

.. code-block:: text

    erf.fire.waf_fcz0_scale = 1.0  (default: 1.0, dimensionless)

When `use_per_fuel_wind_ht = false` (default), all fuel models use the global `wind_ref_ht` 
scalar and `fcz0` is not consulted.

**References:**

- Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling System. 
  USDA Forest Service General Technical Report INT-194.
- WRF-SFIRE documentation: Mandel, J., Beezley, J.D., Kochmann, A.K., et al. (2011). 
  A wildland fire model with data assimilation. *Mathematics and Computers in Simulation*, 
  79(3), 584–606. See module_fr_sfire_phys.F for `fcwh` and `fcz0` data statements.


ROS Model Input Parameters
---------------------------

.. list-table:: Complete ROS Model Configuration Parameters
   :widths: 30 15 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - erf.fire.ros_model
     - string
     - "rothermel"
     - ROS model selector: "rothermel", "macarthur", "balbi", "cheney_gould"
   * - erf.fire.use_per_fuel_wind_ht
     - bool
     - false
     - Enable per-fuel wind height extraction (fcwh table)
   * - erf.fire.waf_fcz0_scale
     - real
     - 1.0
     - Scale factor on fcz0 roughness for WAF computation [dimensionless]
   * - erf.fire.balbi.T_a
     - real
     - 300.0
     - Ambient temperature [K] (Balbi model)
   * - erf.fire.balbi.T_f
     - real
     - 1000.0
     - Mean flame temperature [K] (Balbi model)
   * - erf.fire.balbi.T_i
     - real
     - 600.0
     - Ignition temperature [K] (Balbi model)
   * - erf.fire.balbi.delta_H
     - real
     - 2.26e6
     - Latent heat of vaporisation [J/kg] (Balbi model)
   * - erf.fire.balbi.C_pf
     - real
     - 1800.0
     - Specific heat of dry fuel [J/(kg·K)] (Balbi model)
   * - erf.fire.balbi.r_00
     - real
     - 2.5e-4
     - Radiation length scale [m] (Balbi model)
   * - erf.fire.balbi.tau_0
     - real
     - 75591.0
     - Residence time coefficient [s/m] (Balbi model)
   * - erf.fire.cheney_gould.moisture
     - real
     - 10.0
     - Dead fine fuel moisture [%] (Cheney-Gould model)
   * - erf.fire.cheney_gould.curing
     - real
     - 1.0
     - Degree of curing [0–1] (Cheney-Gould model)


Model Selection Guide
---------------------

.. list-table:: ROS Model Selection by Fuel Type and Region
   :widths: 25 20 55
   :header-rows: 1

   * - Fuel type
     - Recommended model
     - Reason
   * - North American fuels (FM1–FM13)
     - rothermel
     - Calibrated for US fuel types; default model with extensive validation
   * - Australian open grassland (short grass)
     - cheney_gould
     - Empirical calibration for grasslands; better agreement with observed data
   * - Mediterranean shrub, European fuel
     - balbi
     - Physics-based approach; not region-specific; good for diverse fuel types
   * - Australian sclerophyll forest
     - macarthur
     - Mark 5 calibration; designed for Australian forest conditions
   * - Multi-region comparison studies
     - rothermel + alternatives
     - Run multiple models to assess regional sensitivity


Limitations
-----------

- **MacArthur slope effect:** The MacArthur formula does not include a slope effect; slope 
  influence enters only through terrain wind corrections applied before ROS computation.

- **Balbi coefficient simplification:** The Balbi model uses a simplified :math:`A_` 
  coefficient computation (based on fuel load and thermal properties) rather than the full 
  Balbi (2009) flame geometry derivation. See lines 100–103 of ERF_BalbiModel.H.

- **Cheney-Gould placeholder moisture and curing:** The current implementation uses 
  placeholder domain-average moisture and curing values inside the GPU kernel (lines 218–220 
  of ERF_CheneyGouldModel.H). Per-cell moisture from the Phase 4 ODE system is not yet 
  passed into the kernel. This is a known simplification.

- **WAF wind height adjustment:** All four ROS models use the same effective midflame wind 
  (`fire_wind_eff`) after WAF and terrain corrections. The WAF is derived from fuel bed 
  depth using the global fuel model, not per-model wind height adjustment functions.

- **Per-fuel wind height customisation:** Per-fuel wind height (`use_per_fuel_wind_ht = true`) 
  uses WRF-SFIRE default `fcwh` values (all 6.096 m), which are identical for all 13 fuel 
  models. Per-fuel differentiation would require a custom `fcwh` table not yet exposed 
  through ParmParse.


References
----------

- Rothermel, R.C. (1972). *A Mathematical Model for Predicting Fire Spread in Wildland Fuels.* 
  USDA Forest Service Research Paper INT-115.

- McArthur, A.G. (1966). *Weather and Grassland Fire Behaviour.* Leaflet 100, Forestry and 
  Timber Bureau, Canberra.

- Clark, T.L., Coen, J.L. &amp; Latham, D. (2004). Description of a coupled atmosphere-fire model. 
  *International Journal of Wildland Fire*, 13, 49–63.

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model 
  for wildland fires. *Combustion and Flame*, 156(12), 2217–2230.

- Cheney, N.P., Gould, J.S. &amp; Catchpole, W.R. (1998). Prediction of fire spread in 
  grasslands. *International Journal of Wildland Fire*, 8(1), 1–13.

- Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling System. 
  USDA Forest Service General Technical Report INT-194.

- Mandel, J., Beezley, J.D., Kochmann, A.K., Liiskanen, P.K., Vorechova, A., &amp; Halem, M. (2011). 
  A wildland fire model with data assimilation. *Mathematics and Computers in Simulation*, 
  79(3), 584–606.
