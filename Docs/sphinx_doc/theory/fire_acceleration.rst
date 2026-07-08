.. _fire_acceleration:

Fire Acceleration at Startup (Phase 12)
========================================

Physical Background
--------------------

Small fires do not immediately achieve their quasi-steady-state (QSS) rate of spread due to energy and momentum limitations. The startup acceleration period depends on fire size, wind conditions, and fuel type; accurate representation of startup behavior is critical for predictions of early spread extent and spotfire timing. This module implements two optional models to scale the rate of spread during the acceleration phase until the fire reaches equilibrium.

Size-Based Model
-----------------

The size-based model is the simpler option, appropriate when only the current fire extent is needed as a predictor. It scales the rate of spread by a factor that depends on the current burned area:

.. math::

    \alpha = 1 - \exp\!\left(-\frac{r_{\text{fire}}}{L_{\text{acc}}}\right)

    R_{\text{eff}} = \alpha \cdot R_{\text{QSS}}

where:

- :math:`r_{\text{fire}} = \sqrt{A_{\text{fire}} / \pi}` is the effective fire radius derived from burned area [m]
- :math:`A_{\text{fire}}` is the total burned area [m²]
- :math:`L_{\text{acc}}` is the acceleration length scale [m] (ParmParse key: ``erf.fire.accel.L_acc``)
- :math:`R_{\text{QSS}}` is the quasi-steady-state ROS from the active spread model [m/s]

When :math:`r_{\text{fire}} \gg L_{\text{acc}}`, :math:`\alpha \to 1` and the full ROS is used without modification.

References:

- Rothermel, R.C. (1983). How to predict the spread and intensity of forest and range fires. USDA Forest Service General Technical Report INT-143.
- Catchpole, E.A., de Mestre, N.J. &amp; Gill, A.M. (1992). Intensity of fire at its perimeter. *Australian Journal of Ecology*, 17(1), 1–4.

Temporal Model
---------------

The temporal model tracks per-cell acceleration state and is appropriate when the time history of each cell matters. It uses the VanWagner equation:

.. math::

    R(t) = R_E \cdot \left(1 - \exp(-A \cdot t)\right)

where:

- :math:`t` is the elapsed time since the cell entered its current acceleration state [s]
- :math:`A` is the acceleration constant [min⁻¹], converted to [s⁻¹] internally
- :math:`R_E` is the target equilibrium ROS [m/s]

The acceleration constant is selected based on fire perimeter length:

- :math:`A = A_{\text{point}}` when perimeter length &lt; :math:`L_{\text{perim}}` (point ignition)
- :math:`A = A_{\text{line}}` when perimeter length ≥ :math:`L_{\text{perim}}` (line fire)

Alexander et al. (1992) calibrated values are :math:`A_{\text{point}} = 0.115\ \text{min}^{-1}` and :math:`A_{\text{line}} = 0.886\ \text{min}^{-1}`.

Wind-Lag Sub-Section
~~~~~~~~~~~~~~~~~~~~~

When ``erf.fire.accel.enable_wind_lag = true`` and wind speed increases, the target equilibrium ROS is not updated instantaneously but approaches the new equilibrium with time constant :math:`\tau_{\text{wind}}` [s]:

.. math::

    R_{\text{target}}(t) = R_{\text{prev}} + (R_E - R_{\text{prev}}) \cdot \left(1 - \exp\!\left(-\frac{\Delta t}{\tau_{\text{wind}}}\right)\right)

Wind decreases are applied immediately (no lag).

References:

- McAlpine, R.S. &amp; Wakimoto, R.H. (1991). The acceleration of fire from point source to equilibrium spread. *Forest Science*, 37(5), 1314–1337.
- Alexander, M.E., Stocks, B.J. &amp; Lawson, B.D. (1992). Fire behavior in Black Spruce-lichen woodland. Information Report NOR-X-310. Canadian Forest Service.
- Finney, M.A. (1998/2004). FARSITE: Fire Area Simulator — Model Development and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.

Input Parameters
-----------------

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
     - Master on/off switch; when false, acceleration is disabled
   * - ``erf.fire.accel.use_temporal``
     - bool
     - false
     - false = size-based model, true = temporal per-cell model
   * - ``erf.fire.accel.L_acc``
     - Real [m]
     - 50.0
     - Size-based: acceleration length scale [m]
   * - ``erf.fire.accel.A_point``
     - Real [1/min]
     - 0.115
     - Temporal: acceleration constant for point ignitions [1/min]
   * - ``erf.fire.accel.A_line``
     - Real [1/min]
     - 0.886
     - Temporal: acceleration constant for line ignitions [1/min]
   * - ``erf.fire.accel.perim_limit``
     - Real [m]
     - 500.0
     - Temporal: perimeter length threshold [m] for switching from A_point to A_line
   * - ``erf.fire.accel.enable_wind_lag``
     - bool
     - false
     - Temporal: apply exponential lag on wind speed increases
   * - ``erf.fire.accel.tau_wind``
     - Real [s]
     - 60.0
     - Temporal: wind-lag time constant [s]

Usage Example
--------------

Size-based model (simpler, uses fire area only)::

   erf.fire.accel.enable       = true
   erf.fire.accel.use_temporal = false
   erf.fire.accel.L_acc        = 100.0

Temporal model (per-cell, tracks elapsed time)::

   erf.fire.accel.enable          = true
   erf.fire.accel.use_temporal    = true
   erf.fire.accel.A_point         = 0.115
   erf.fire.accel.A_line          = 0.886
   erf.fire.accel.perim_limit     = 500.0
   erf.fire.accel.enable_wind_lag = false
   erf.fire.accel.tau_wind        = 60.0

Limitations
-----------

- Both models are disabled by default (``accel.enable = false``); enabling has zero runtime cost when disabled.
- The size-based model uses a single global scaling factor; cells near the fire centre and cells at the expanding perimeter receive the same alpha.
- The temporal model resets the acceleration clock whenever the equilibrium ROS changes by more than a threshold (relative threshold: 1.0e-6 times the maximum of the new and previous ROS values); this threshold applies to both wind changes and fuel changes.
- The wind-lag model applies only to wind speed increases; wind decreases update the target ROS immediately.
- The temporal model requires allocating a 3-component per-cell ``fire_accel_state`` MultiFab; the size-based model does not allocate any additional storage.
- Neither model adjusts the fire perimeter shape — only the magnitude of ROS is scaled.

References
-----------

- Rothermel, R.C. (1983). How to predict the spread and intensity of forest and range fires. USDA Forest Service General Technical Report INT-143.
- Catchpole, E.A., de Mestre, N.J. &amp; Gill, A.M. (1992). Intensity of fire at its perimeter. *Australian Journal of Ecology*, 17(1), 1–4.
- McAlpine, R.S. &amp; Wakimoto, R.H. (1991). The acceleration of fire from point source to equilibrium spread. *Forest Science*, 37(5), 1314–1337.
- Alexander, M.E., Stocks, B.J. &amp; Lawson, B.D. (1992). Fire behavior in Black Spruce-lichen woodland. Information Report NOR-X-310. Canadian Forest Service.
- Finney, M.A. (1998/2004). FARSITE: Fire Area Simulator — Model Development and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
