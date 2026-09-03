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

Physics-based fire spread model derived from the radiant energy flux received by unburned 
fuel from a tilted flame front, rather than from empirical fitting. The model couples wind 
speed, terrain slope, fuel properties and fuel moisture through the flame tilt angle.

**Formulae:**

.. math::

    R = A \cdot (1 + \sin\alpha - \cos\alpha)

    A = \frac{\chi \, \sigma_m \, \delta_m}{2 \, \tau_0 \, B^*}

    \chi = \frac{r_{00} \sigma_m}{1 + r_{00} \sigma_m}

    B^* = \frac{C_{pf} (T_i - T_a) + M_f \Lambda}{h}

    \tan\alpha = \frac{U}{v_b} + \tan\theta

    v_b = \sqrt{\frac{g \cdot \delta_m \cdot (T_f - T_a)}{T_a}}

where:

- :math:`R` [m/s] is the rate of spread
- :math:`A` [m/s] is the ROS amplitude coefficient
- :math:`\chi` [-] is the fraction of combustion heat radiated forward
- :math:`B^*` [-] is the dimensionless energy needed to ignite the fuel
- :math:`\alpha` [radians] is the flame tilt angle
- :math:`U` [m/s] is the effective midflame wind speed
- :math:`\theta` [radians] is the terrain slope
- :math:`v_b` [m/s] is the buoyancy velocity
- :math:`g = 9.81` m/s² is gravitational acceleration
- :math:`\sigma_m` [1/m] is the surface-area-to-volume ratio of the 1-hr dead fuel
- :math:`\delta_m` [m] is the fuel bed depth
- :math:`h` [J/kg] is the heat of combustion
- :math:`M_f` [-] is the fuel moisture content
- :math:`\Lambda` [J/kg] is the latent heat of vaporisation
- :math:`T_f` [K] is the mean flame temperature
- :math:`T_i` [K] is the ignition temperature
- :math:`T_a` [K] is the ambient temperature

Fuel properties enter :math:`A` through :math:`\sigma_m`, :math:`\delta_m` and :math:`h`, so 
the ROS varies by fuel model. When a spatial fuel map is loaded, per-cell coefficients come 
from a lookup table built for the 13 Anderson fuel models, with fuel code 0 (non-burnable) 
giving zero spread. When `moisture_dynamic = true`, :math:`B^*` and hence :math:`A` are 
rebuilt each step from the domain-average 1-hr moisture.

**Balbi (2020) convective-radiative form**

Setting `erf.fire.balbi.formulation = "2020"` selects the convective-radiative model of 
Balbi et al. (2020), in which the rate of spread is the root of

.. math::

    R = R_b + R_c + R_r

    R_b = \min\!\left(\frac{S}{\pi}, 1\right) \frac{B T^4}{\beta \rho_v q}

    R_c = \frac{s \Delta H}{q \tau_0} \min\!\left(\delta_m, \frac{2\pi}{s\beta}\right)
          \left[ \frac{\delta_m}{2\delta_m + H}\tan\alpha
                  + \frac{U e^{-K_1 \sqrt{\beta} R}}{u_0} \right]

    R_r = A R \frac{1 + \sin\gamma - \cos\gamma}{1 + R\cos\gamma/(s r_{00})}

where :math:`\beta = w/(\delta_m \rho_v)` is the packing ratio, :math:`S = s\beta\delta_m` 
is twice the leaf area index, :math:`u_0` is the vertical gas velocity, :math:`H` the flame 
height and :math:`\gamma` the flame tilt. The radiative base term :math:`R_b` gives a nonzero 
rate of spread with no wind and no slope, and the wind response does not saturate — the two 
structural limits of the 2009 form.

The equation is solved per cell by bisection on :math:`f(R) = R_b + R_c + R_r - R`, which is 
positive at :math:`R = 0` and negative at the 30 m/s cap. Plain substitution, used by the 
reference implementations, oscillates instead of converging when :math:`A \ll 1`.

Under this formulation `erf.fire.balbi.r_00` defaults to :math:`2.5\times10^{-5}` — the value 
fitted by Balbi et al. (2020) — rather than the 2009 radiation length scale of 
:math:`2.5\times10^{-4}` m.

**Optional couplings (both formulations, off by default)**

- `erf.fire.balbi.directional` evaluates the ROS along the front normal 
  :math:`\hat{n} = \nabla\phi/|\nabla\phi|`, giving head, flank and backing spread from the 
  model itself. Level-set path only; the FARSITE path gets directionality from the ellipse.
- `erf.fire.balbi.use_surface_temp` takes the ambient temperature per cell from 
  `fire_surface_temp`.
- `erf.fire.balbi.heat_flux_coupling` augments the vertical velocity scale with 
  :math:`v_{b,Q} = k_{up}\sqrt{g Q H_{ref}/(\rho_a C_{pa} T_a)}` in quadrature. A stronger 
  plume stands the flame up and slows the head fire. Lags the ROS by one fire step.
- `erf.fire.balbi.flame_temp_from_combustion` derives the 2009 flame temperature from the 
  combustion energy instead of the fixed `balbi.T_f`.
- `erf.fire.balbi.use_cell_moisture` takes the 1-hr dead moisture per cell from the Phase 4 
  moisture ODE state instead of the domain average. Requires `moisture_dynamic = true`.
- `erf.fire.balbi.use_moisture_extinction` zeroes the ROS at and above the fuel model's 
  moisture of extinction. Neither formulation has an extinction limit of its own, so without 
  this a fuel bed wetter than its :math:`M_x` still spreads.
- `erf.fire.balbi.herb_curing` [0-1] is the fraction of the live herbaceous load carried as 
  cured dead fine fuel, entering the 2020 packing ratio as 
  :math:`w = w_{d1} + \text{curing}\,w_{lh}`. No effect under the 2009 form, whose amplitude 
  coefficient does not depend on loading.
- `erf.fire.balbi.wind_source` selects `"midflame"` (default, `fire_wind_eff`, after the Wind 
  Adjustment Factor and terrain corrections) or `"reference"` (`fire_wind_ref`, at 
  `wind_ref_ht`, before both). The WAF is a Rothermel calibration construct; Balbi normalises 
  the wind by its own vertical velocity scale instead, so applying the WAF first is arguably 
  a double reduction.

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
- The default (2009) form has no no-wind radiative base spread (:math:`R \to 0` as wind and
  slope vanish) and saturates at :math:`2A` for large flame tilt; use
  `formulation = "2020"` for the convective-radiative form, which has neither limit
- Without `balbi.directional`, the flame tilt uses the wind magnitude and the slope magnitude
  :math:`|\nabla z|`, so the ROS field is isotropic and downslope spread is enhanced like
  upslope spread
- The 2020 form uses the 1-hr dead fuel load and the fuel particle density, plus the cured
  fraction of the live herbaceous load; coarse dead fuels and uncured live fuels do not
  enter it
- Without `balbi.use_moisture_extinction`, neither formulation extinguishes on wet fuel: the
  ROS falls with moisture but never reaches zero
- Small buoyancy velocities (:math:`v_b < 10^{-3}` m/s) are floored to prevent division by zero
- Cells with wind above 25 m/s, or a non-finite wind, fall back to the domain-mean wind

**ParmParse key:**

.. code-block:: text

    erf.fire.ros_model = "balbi"

**References:**

- Balbi, J.-H., Rossi, J.-L., Marcelli, T. &amp; Santoni, P.-A. (2009). A physical model 
  for wildland fires. *Combustion and Flame*, 156(12), 2217–2230.

- Balbi, J.-H., Chatelon, F.-J., Morvan, D., Rossi, J.-L., Marcelli, T. &amp; Morandini, F. 
  (2020). A convective–radiative propagation model for wildland fires. *International Journal 
  of Wildland Fire*, 29(8), 723–738.


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

- **Balbi formulation default:** The Balbi model defaults to the steady explicit form of
  Balbi (2009), which has no no-wind radiative base spread and saturates at :math:`2A`. The
  convective-radiative form of Balbi (2020), which restores both, is available through
  `erf.fire.balbi.formulation = "2020"` but is not the default.

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
