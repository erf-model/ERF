.. _fire_acceleration:

Fire Acceleration at Startup (Phase 12)
========================================

Physical Background
--------------------

Small fires do not immediately achieve their quasi-steady-state (QSS) rate of spread due to 
energy and momentum limitations. The startup acceleration period is important for prescribed burn 
timing, spotfire scenarios, and realistic representation of early-stage fire dynamics. This module 
implements two optional models to reduce ROS during the acceleration phase until the fire reaches 
equilibrium spread.

Size-Based Model
-----------------

The size-based model (Rothermel 1983, Catchpole et al. 1992) scales the rate of spread by a factor 
that depends on the current burned area:

.. math::

    \alpha = 1 - \exp\left(-\frac{r_{\text{fire}}}{L_{\text{acc}}}\right)

    R_{\text{eff}} = \alpha \cdot R_{\text{QSS}}

where:

- :math:`r_{\text{fire}}` is the effective fire radius derived from burned area as :math:`r_{\text{fire}} = \sqrt{A_{\text{fire}} / \pi}`
- :math:`L_{\text{acc}}` [m] is the acceleration length scale (default: 50 m)
- :math:`R_{\text{QSS}}` is the quasi-steady-state Rothermel ROS
- :math:`\alpha` varies from 0 (point fire) to 1 (large fire)

This model applies a uniform global scaling factor. When the fire is large relative to :math:`L_{\text{acc}}`, 
:math:`\alpha \to 1` and the fire spreads at full Rothermel speed. This simple approach is computationally 
efficient and physically motivated by the observation that fire intensity scales with burned perimeter.

**References:**

- Rothermel, R.C. (1983). "How to predict the spread and intensity of forest and range fires." 
  USDA Forest Service General Technical Report INT-143.
- Catchpole, E.A., de Mestre, N.J., & Gill, A.M. (1992). "Intensity of fire at its perimeter." 
  *Australian Journal of Ecology*, 17(1), 1–4.

Temporal Model
---------------

The temporal model (McAlpine & Wakimoto 1991, FARSITE Finney 1998/2004) tracks per-cell acceleration 
state using the VanWagner equation:

.. math::

    R(t) = R_E \cdot \left(1 - \exp(-A \cdot t)\right)

where:

- :math:`t` [s] is the time elapsed since ignition or equilibrium ROS change
- :math:`R_E` [m/s] is the equilibrium (target) ROS from the Rothermel model
- :math:`A` [min⁻¹] is the acceleration constant

The acceleration constant is selected based on fire perimeter length:

- **Point fire** (perimeter < ``perim_limit``): :math:`A = A_{\text{point}}` (default: 0.115 min⁻¹)
- **Line fire** (perimeter ≥ ``perim_limit``): :math:`A = A_{\text{line}}` (default: 0.886 min⁻¹)

where ``perim_limit`` [m] is a threshold perimeter (default: 500 m). This allows fires to transition 
from point-ignition acceleration to line-fire acceleration as they spread.

The model includes optional wind-lag damping: when wind speed increases, the target equilibrium ROS 
approaches the new value exponentially rather than jumping instantly:

.. math::

    R_{\text{target}} = R_{\text{old}} + (R_{\text{new}} - R_{\text{old}}) \cdot \left(1 - \exp\left(-\frac{\Delta t}{\tau_{\text{wind}}}\right)\right)

This lag provides a more realistic representation of fire response to sudden wind changes.

**References:**

- McAlpine, R.S. & Wakimoto, R.H. (1991). "The acceleration of fire from point source to equilibrium spread." 
  *Forest Science*, 37(5), 1314–1337.
- Alexander, M.E., Stocks, B.J., & Lawson, B.D. (1992). "Fire behavior in Black Spruce-lichen woodland." 
  Information Report NOR-X-310, Northern Forestry Centre, Canadian Forest Service.
- Finney, M.A. (1998/2004). "FARSITE: Fire Area Simulator." USDA Forest Service Research Paper RMRS-RP-4.

Input Parameters
-----------------

Fire acceleration is configured via ``erf.fire.accel.*`` parameters in the input file. 
All parameters are optional and have sensible defaults.

.. list-table:: Fire Acceleration Parameters
   :widths: 20 15 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``erf.fire.accel.enable``
     - bool
     - false
     - Master on/off switch. When false, acceleration is disabled and fire spreads at full ROS immediately.
   * - ``erf.fire.accel.use_temporal``
     - bool
     - false
     - Model selection: false = size-based, true = temporal per-cell.
   * - ``erf.fire.accel.L_acc``
     - Real [m]
     - 50.0
     - Size-based: acceleration length scale. Controls how quickly alpha approaches 1 with fire size.
   * - ``erf.fire.accel.A_point``
     - Real [1/min]
     - 0.115
     - Temporal: acceleration constant for point ignitions (small perimeter).
   * - ``erf.fire.accel.A_line``
     - Real [1/min]
     - 0.886
     - Temporal: acceleration constant for line ignitions (large perimeter).
   * - ``erf.fire.accel.perim_limit``
     - Real [m]
     - 500.0
     - Temporal: perimeter length threshold for switching from A_point to A_line.
   * - ``erf.fire.accel.enable_wind_lag``
     - bool
     - false
     - Temporal: apply exponential lag when wind increases.
   * - ``erf.fire.accel.tau_wind``
     - Real [s]
     - 60.0
     - Temporal: wind-lag time constant. Controls how quickly wind changes are applied.

Usage Notes
-----------

- **Default disabled:** Acceleration is disabled by default (``enable = false``). 
  When disabled, there is zero runtime cost — no MultiFab is allocated and the acceleration call returns immediately.

- **Model selection:** The size-based model is simpler and uses less memory (no per-cell state).
  The temporal model is more realistic but requires allocating a 3-component MultiFab.

- **Acceleration length scale:** For ``L_acc``, a value around 50–100 m is typical for small initial fires.
  Larger values produce slower acceleration; smaller values produce faster acceleration.

- **Equilibrium ROS changes:** In the temporal model, when the equilibrium ROS changes (e.g., due to wind change),
  the cell immediately resets ``t_elapsed = 0`` and begins accelerating toward the new equilibrium.
  With wind-lag enabled, the acceleration resets toward a lagged intermediate target instead.

- **Fire size effect:** The size-based model stops printing debug messages once alpha exceeds 0.9999
  (fire is effectively at full ROS). The temporal model continues per-cell acceleration throughout.

Limitations
-----------

- **Size-based uses global alpha:** The size-based model applies a single scaling factor to all cells.
  This does not account for local variations in fire age or wind.

- **Temporal model resets on any change:** When :math:`|R_E - R_{\text{prev}}|` exceeds a threshold,
  the acceleration resets. This occurs not just for wind changes but also for fuel or moisture changes.
  This is physically reasonable for wind changes but may be overly aggressive for gradual moisture changes.

- **Wind-lag is first-order only:** The wind-lag formula is a simple first-order exponential approximation.
  Real wind response may be more complex.

Example Usage
--------------

**Size-based acceleration with default parameters:**

.. code-block:: ini

    erf.fire.accel.enable = true
    erf.fire.accel.use_temporal = false

**Temporal acceleration with wind-lag:**

.. code-block:: ini

    erf.fire.accel.enable = true
    erf.fire.accel.use_temporal = true
    erf.fire.accel.enable_wind_lag = true
    erf.fire.accel.tau_wind = 120.0

**Acceleration disabled (baseline):**

.. code-block:: ini

    erf.fire.accel.enable = false

When disabled, existing Phase 3/4 behavior is exactly reproduced — fires spread immediately at 
full Rothermel ROS regardless of size or time.
