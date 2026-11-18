.. _MoistPotentialTemperatures:

Moist Potential Temperatures
============================

ERF evolves a single thermodynamic variable, the dry potential temperature
:math:`\theta`, through the conserved quantity :math:`\rho_d \theta`. Several
other potential temperatures are built from :math:`\theta`, the temperature
:math:`T` and the moisture mixing ratios. Each one answers a different
question. This page defines each one as ERF computes it and lists where it is
used. Other pages refer here rather than restating the formulas.

Notation
--------

All mixing ratios are per unit mass of dry air,
:math:`q_x = (\rho_d q_x)/\rho_d`, the form in which ERF stores moisture.

- :math:`q_v` is water vapor.
- :math:`q_c` and :math:`q_i` are cloud liquid and cloud ice. These are the
  non-precipitating condensate.
- :math:`q_r`, :math:`q_s` and :math:`q_g` are rain, snow and graupel. These
  are the precipitating condensate.

A species that the chosen microphysics scheme does not carry is taken to be
zero.

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Constant
     - Value
     - Notes
   * - :math:`R_d`, :math:`R_v`
     - 287.0, 461.505 J kg\ :sup:`-1` K\ :sup:`-1`
     - gas constants for dry air and water vapor
   * - :math:`c_{p,d}`
     - 1004.5 J kg\ :sup:`-1` K\ :sup:`-1`
     - gives :math:`\gamma = 1.4` for dry air
   * - :math:`c_l`
     - 4200 J kg\ :sup:`-1` K\ :sup:`-1`
     - liquid water
   * - :math:`\epsilon_v = R_v/R_d - 1`
     - :math:`\approx 0.608`
     - the "0.61" of the meteorological literature
   * - :math:`L_v`
     - :math:`2.5 \times 10^6` J kg\ :sup:`-1`
     - latent heat of vaporization
   * - :math:`L_s`
     - :math:`2.8347 \times 10^6` J kg\ :sup:`-1`
     - latent heat of sublimation (``lsub`` = ``lcond`` + ``lfus``)
   * - :math:`p_0`
     - :math:`10^5` Pa
     - reference pressure

Moist potential temperature :math:`\theta_m`
--------------------------------------------

.. math::
   \theta_m = \theta \left( 1 + \frac{R_v}{R_d} q_v \right)

:math:`\theta_m` is the vapor factor in the equation of state,
:math:`p = p_0 (R_d \rho_d \theta_m / p_0)^\gamma`. ERF uses it to recover
pressure and temperature from :math:`\rho_d`, :math:`\rho_d\theta` and
:math:`q_v` (see :ref:`GoverningEquations` and :ref:`sec:BaseState`). It is not a
stability variable. It is not the virtual potential temperature either,
because it carries the vapor's gas constant but not the mass of the water.

Virtual potential temperature :math:`\theta_v`
----------------------------------------------

.. math::
   \theta_v = \theta \left( 1 + \epsilon_v q_v - q_c - q_r - q_i - q_s - q_g \right)

:math:`\theta_v` measures the buoyancy of moist, cloudy air: it is the
potential temperature that dry air would need to have the same density.
Vapor makes air lighter and raises :math:`\theta_v`. Every condensed species,
precipitating or not, adds mass and lowers it. The condensate here is the
same total that ERF uses elsewhere for the water loading (``make_qt``). The
form above is the usual linearization in the mixing ratios.

The code computes it with ``GetThetav`` in ``Source/Utils/ERF_MoistUtils.H``.

:math:`\theta_v` is used in:

- the Smagorinsky Richardson-number correction, in unsaturated cells
  (``erf.use_Ri_correction``; see :ref:`SmagorinskyModel`);
- the PBL-height diagnosis (minimum :math:`\theta_v` near the surface and the
  :math:`\theta_v`-increase criterion);
- the MRF and YSU schemes, wherever WRF uses ``thvx``, including the bulk
  Richardson number;
- the PBL schemes' stratification variable when
  ``erf.thermal_stratification = thetav``, the default for PBL schemes;
- the :math:`\langle w'\theta_v' \rangle` column of the 1D profile output.
  :math:`\theta_v` is computed in each cell and then horizontally averaged.

Liquid-ice water potential temperature :math:`\theta_l`
-------------------------------------------------------

.. math::
   \theta_l = \theta - \frac{\theta}{T} \, \frac{L_v q_c + L_s q_i}{c_{p,d}}

:math:`\theta_l` is approximately conserved under reversible phase changes.
When cloud water condenses, the latent heat raises :math:`\theta`, but
:math:`\theta_l` stays the same. This makes it the natural stability variable
inside clouds. The factor :math:`\theta/T = 1/\Pi` converts the latent heating
from a temperature change to a potential-temperature change.

Only the cloud condensate :math:`q_c` and :math:`q_i` enters. Precipitation
falls out of the parcel and is not part of the reversible budget. The ice
term uses the latent heat of sublimation. This is the liquid-ice potential
temperature, sometimes written :math:`\theta_{li}`, and it matches WRF's
``thlix`` in YSU. ERF has no separate liquid-only variant, so
:math:`\theta_{li}` and :math:`\theta_l` refer to the same quantity.

The code computes it with ``GetThetal`` in ``Source/Utils/ERF_MoistUtils.H``.

:math:`\theta_l` is used in:

- the Smagorinsky Richardson-number correction, in saturated cells
  (:math:`q_c > 10^{-8}` kg/kg);
- the YSU top-down PBL-height extension, and the cloudy entrainment jump at
  the PBL top (see :ref:`PBLschemes`). Following WRF, the jump is taken in
  :math:`\theta_l + \epsilon_v \theta (q_v + q_c)`, evaluated at the levels
  just below the PBL top and two levels above it;
- the PBL schemes' stratification variable when
  ``erf.thermal_stratification = thetal``.

Equivalent potential temperature :math:`\theta_e`
-------------------------------------------------

.. math::
   \theta_e = T \left( \frac{p - e}{p_0} \right)^{-R_d/c_{pt}}
              \mathcal{H}^{-q_v R_v / c_{pt}}
              \exp\!\left( \frac{L_v q_v}{c_{pt} T} \right),
   \qquad c_{pt} = c_{p,d} + c_l q_t

This is Emanuel (1994), eq. 4.5.11, with a constant :math:`L_v`.
:math:`\theta_e` is the temperature that a parcel would reach if all its vapor
condensed and the latent heat warmed it, after which it is brought
dry-adiabatically to :math:`p_0`. It is conserved under moist-adiabatic
ascent. The terms are:

- :math:`e = p\, q_v / (R_d/R_v + q_v)` is the vapor partial pressure, so
  :math:`p - e` is the dry-air partial pressure.
- :math:`q_t` is the total water mixing ratio, vapor plus every condensed
  species.
- :math:`\mathcal{H} = e / e_s(T)` is the relative humidity with respect to
  liquid water. It is capped at 1 and floored at :math:`10^{-10}`, so the
  factor is exactly 1 for dry air.

:math:`\theta_e` is a diagnostic only. It is available as the plotfile
variable ``eq_pot_temp``. In anelastic runs :math:`T` comes from :math:`p` and
:math:`\theta`; in compressible runs it comes from the equation of state. The
``eq_pot_temp`` input of the squall-line, supercell and bubble problems sets
an initial :math:`\theta_e` profile (see :ref:`sec:ProblemInputs`).

Summary
-------

.. list-table::
   :header-rows: 1
   :widths: 12 38 50

   * - Symbol
     - Moisture included
     - Purpose in ERF
   * - :math:`\theta`
     - none
     - prognostic variable (:math:`\rho_d\theta`)
   * - :math:`\theta_m`
     - :math:`q_v` (gas constant only)
     - equation of state
   * - :math:`\theta_v`
     - :math:`q_v` and all condensate (loading)
     - buoyancy, stability, PBL height, :math:`w'\theta_v'` profiles
   * - :math:`\theta_l`
     - :math:`q_c`, :math:`q_i` (latent heat)
     - in-cloud stability (Smagorinsky, YSU)
   * - :math:`\theta_e`
     - :math:`q_v` (latent heat), :math:`q_t` (heat capacity), RH
     - diagnostic output, initial profiles
