
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Radiation:

Radiation
=========

Radiative transfer in ERF includes both a full k-distribution model (RRTMGP) for detailed studies 
and a simplified two-stream model for idealized and intermediate-complexity atmospheric simulations. 
This section describes the two-stream radiation model's physics, which employs Beer-Lambert direct-beam 
attenuation for shortwave radiation, Meador-Weaver two-stream diffuse/scattering, gray-gas longwave 
two-stream, and an optional Simplified Surface Energy Balance (SEB) module.

Shortwave Radiation
===================

The shortwave (solar) radiation calculation in the two-stream model is split into direct-beam 
and diffuse components.

Direct-Beam Attenuation (Beer-Lambert)
--------------------------------------

The shortwave direct-beam radiation is attenuated through the atmosphere according to the 
Beer-Lambert law:

.. math::

   I_{sw,\text{direct}}(z) = I_0 \mu_0 e^{-\tau_{\text{sw}} \sec(\theta_z)}

where :math:`I_0` is the solar constant at the top of the atmosphere (:math:`S_0 \approx 1361 \, \text{W/m}^2`), 
:math:`\mu_0 = \cos(\theta_z)` is the cosine of the solar zenith angle, :math:`\tau_{\text{sw}}` is the 
vertically integrated shortwave optical depth above height :math:`z`, and :math:`\sec(\theta_z)` accounts 
for the path-length modification. The optical depth may be spatially uniform (static :math:`\tau_{\text{per\_layer}}`) 
or dynamically diagnosed from water vapor and cloud liquid water content, parameterized as:

.. math::

   \tau_{\text{sw}}(k) = \tau_{\text{per\_layer}} + \tau_{\text{cloud}}(k) + c_{\text{qv}} q_v(k) + c_{\text{qc}} q_c(k)

where :math:`\tau_{\text{cloud}}(k)` is added only within the prescribed cloud layer, and the coefficients 
:math:`c_{\text{qv}}` and :math:`c_{\text{qc}}` are zero by default.

Diffuse Shortwave and Scattering (Meador-Weaver Two-Stream)
-----------------------------------------------------------

The diffuse (scattered) shortwave component is computed using the Meador and Weaver (1980) two-stream 
approximation. In each layer :math:`k`, the radiative transfer equation is discretized into upward and 
downward diffuse streams. The single-scattering albedo :math:`\omega_0` and asymmetry factor :math:`g` 
parameterize the scattering properties:

.. math::

   \omega_0 = \frac{\text{scattering cross-section}}{\text{total extinction cross-section}}, \quad 
   g = \left\langle \cos(\theta) \right\rangle_{\text{scattering}}

For a horizontally homogeneous layer with optical depth :math:`\tau` and parameters :math:`\omega_0`, :math:`g`, 
the two-stream approximation computes layer reflectance :math:`r` and transmittance :math:`t`:

.. math::

   r = \frac{\gamma_1 (1 - e^{-\lambda \tau})}{\lambda + \gamma_1 (1 - e^{-\lambda \tau})}, \quad 
   t = \frac{\lambda e^{-\lambda \tau}}{\lambda + \gamma_1 (1 - e^{-\lambda \tau})}

where :math:`\gamma_1` and :math:`\lambda` are Meador-Weaver gamma coefficients defined by:

.. math::

   \gamma_1 = \frac{1 - \omega_0}{2}, \quad 
   \lambda = \sqrt{3 (1 - \omega_0) (1 - \omega_0 g)}

The single-scattering source term (backscattered direct beam into the diffuse field) is:

.. math::

   S_{\text{ss}} = (1 - \omega_0) \mu_0 e^{-\tau / \mu_0} / (1 + 2 \gamma_1)

The diffuse upward and downward fluxes are then computed through an upward sweep followed by a downward 
sweep, applying boundary conditions at the surface (with user-specified or LSM-provided albedo :math:`\alpha_{\text{sw}}`) 
and at the top-of-atmosphere (no upward flux). The net diffuse shortwave flux at each level is:

.. math::

   F_{\text{sw,net}}(k) = F_{\text{sw,down}}(k) - F_{\text{sw,up}}(k)

Cloud scattering properties may differ from the clear-sky values (e.g., :math:`\omega_0^{\text{cloud}}` for 
liquid water clouds). The cloud/clear-sky distinction is blended according to the cloud fraction :math:`C_f`:

.. math::

   F_{\text{sw,net}} = (1 - C_f) F_{\text{sw,net}}^{\text{clear}} + C_f F_{\text{sw,net}}^{\text{cloudy}}

Longwave Radiation
==================

The longwave (thermal) radiation employs a gray-gas two-stream formulation (Toon et al. 1989) 
to compute upward and downward fluxes.

Gray-Gas Two-Stream Formulation
--------------------------------

Assume local thermodynamic equilibrium (LTE): each layer emits radiation according to the 
Planck function :math:`B(T)` weighted by the emissivity of the gray gas. The optical depth is 
parameterized similarly to shortwave:

.. math::

   \tau_{\text{lw}}(k) = \tau_{\text{lw,per\_layer}} + \tau_{\text{cloud,lw}}(k) + c_{\text{qv,lw}} q_v(k) + c_{\text{qc,lw}} q_c(k)

In each layer, the two-stream equations for upward (:math:`F_{\uparrow}`) and downward (:math:`F_{\downarrow}`) 
fluxes are:

.. math::

   \frac{dF_{\uparrow}}{d\tau} = F_{\uparrow} - 2 B(T), \quad 
   \frac{dF_{\downarrow}}{d\tau} = -F_{\downarrow} + 2 B(T)

These are integrated over each model layer using an upward sweep (from the surface upward) and a downward 
sweep (from the top-of-atmosphere downward), with boundary conditions:

- **At the surface** (:math:`z = 0`, :math:`\tau = \tau_{\text{col}}`): The surface emits according to 
  Stefan-Boltzmann with emissivity :math:`\epsilon_{\text{lw}}`:

  .. math::

     F_{\uparrow}(0) = \epsilon_{\text{lw}} \sigma T_s^4 + (1 - \epsilon_{\text{lw}}) F_{\downarrow}(0)

  where :math:`\sigma = 5.67 \times 10^{-8} \, \text{W/(m}^2\text{·K}^4)` is the Stefan-Boltzmann constant 
  and :math:`T_s` is the surface temperature.

- **At the top-of-atmosphere** (:math:`z = z_{\text{top}}`): No downward flux from space:

  .. math::

     F_{\downarrow}(z_{\text{top}}) = 0

The net longwave flux divergence in each layer drives the longwave heating rate:

.. math::

   H_{\text{lw}} = -\frac{1}{\rho c_p} \frac{\partial}{\partial z} (F_{\uparrow} - F_{\downarrow})

where :math:`\rho` is air density and :math:`c_p` is the specific heat at constant pressure.

Surface Energy Balance
======================

The surface energy balance residual is defined as the net radiative flux minus the turbulent and 
ground heat fluxes:

.. math::

   R_{\text{net}} = F_{\text{sw,net}}(0^+) + F_{\text{lw,net}}(0^+)

   \text{SEB}_{\text{residual}} = R_{\text{net}} - H - \text{LE} - G

where :math:`H` is the sensible heat flux, :math:`\text{LE}` is the latent heat flux (evaporation), 
and :math:`G` is the ground heat flux conducted into the soil.

Force-Restore Method for Prognostic Surface Temperature and Moisture
=====================================================================

When the Simplified SEB prognostic mode is enabled (``seb_prognostic_enable = true``), the surface 
temperature and moisture are evolved forward in time using a force-restore formulation 
(Tremback and Kessler, 1985; cf. Bhumralkar, 1974).

Surface Temperature Evolution
------------------------------

The surface temperature :math:`T_s` evolves according to:

.. math::

   C_s \frac{dT_s}{dt} = R_{\text{net}} - H - \text{LE} - G - C_s \left( \frac{2\pi}{\tau} \right) (T_s - T_{\text{deep}})

where :math:`C_s` is the effective heat capacity of the surface-active layer [J/(m²·K)], :math:`\tau` is 
the force-restore timescale [s] (e.g., 86400 s or 1 day), and :math:`T_{\text{deep}}` is the deep-soil 
temperature representing the climate state at that location.

In discretized form (Euler forward step), the update is:

.. math::

   T_s^{n+1} = T_s^n + \Delta t \left[ \frac{R_{\text{net}} - H - \text{LE} - G}{C_s} - \left( \frac{2\pi}{\tau} \right) (T_s^n - T_{\text{deep}}) \right]

After the update, :math:`T_s` is clamped to physically reasonable bounds [``seb_prognostic_t_min_k``, ``seb_prognostic_t_max_k``].

Surface Moisture Evolution
---------------------------

The surface moisture (assumed to be water in a thin surface-active layer of depth :math:`d_s`) evolves as:

.. math::

   \frac{dq_s}{dt} = -\frac{\text{LE}}{L_v \rho_w d_s} - \frac{1}{\tau_q} (q_s - q_{\text{deep}})

where :math:`L_v = 2.5 \times 10^6 \, \text{J/kg}` is the latent heat of vaporization, :math:`\rho_w \approx 1000 \, \text{kg/m}^3` 
is the density of liquid water, :math:`d_s` is the effective moisture layer depth [m], :math:`\tau_q` is the 
moisture force-restore timescale [s], and :math:`q_{\text{deep}}` is the deep-soil moisture.

In discretized form:

.. math::

   q_s^{n+1} = q_s^n + \Delta t \left[ -\frac{\text{LE}^n}{L_v \rho_w d_s} - \frac{1}{\tau_q} (q_s^n - q_{\text{deep}}) \right]

After the update, :math:`q_s` is clamped to [``seb_prognostic_q_min``, ``seb_prognostic_q_max``].

Noah-MP Precedence and Double-Counting Safeguard
-------------------------------------------------

When Noah-MP is active at a particular level, the SEB prognostic update is automatically skipped at that level, 
and Noah-MP's own surface prognostics (which include soil heat conduction and explicit soil moisture layers) 
are used instead. This prevents double-counting of surface energy and moisture evolution.

Cloud Fraction Diagnosis
========================

When ``cloud_fraction_prog_enable = true``, the cloud fraction is diagnosed at each level from the 
relative humidity (RH) and cloud liquid water content (qc):

.. math::

   C_f = \min \left( 1, \max \left( 0, \frac{RH - \text{rh\_min}}{\text{rh\_max} - \text{rh\_min}} \right) + c_{\text{qc}} q_c \right)

where the RH threshold parameters allow for a transition from clear sky (RH < rh_min) to complete cloud coverage 
(RH >= rh_max). The coefficient :math:`c_{\text{qc}}` (default :math:`1 \times 10^{-3}`) provides an additional 
scaling of liquid water content's contribution. Optional temporal smoothing via an exponential moving average 
(EMA) may be applied to suppress oscillations:

.. math::

   C_f^{\text{smoothed}} = \alpha C_f^{\text{new}} + (1 - \alpha) C_f^{\text{old}}

where :math:`\alpha` is the blending parameter.

Solar Geometry and Diurnal Cycle
================================

When ``solar_geometry_dynamic_enable = true``, the solar zenith angle is computed dynamically from 
astronomical formulas based on the simulation time, latitude, longitude, day-of-year, and time-zone offset.

The solar declination :math:`\delta` (angle of the sun relative to the Earth's equatorial plane) is 
computed from the day-of-year :math:`D`:

.. math::

   \delta = 23.45° \sin \left( \frac{2\pi (D - 81)}{365} \right)

where the 81st day is the spring equinox.

The hour angle :math:`h` (solar time in degrees, 0 at solar noon, 15° per hour) is computed from the 
local solar time. The solar zenith angle is then:

.. math::

   \cos(\theta_z) = \sin(\phi) \sin(\delta) + \cos(\phi) \cos(\delta) \cos(h)

where :math:`\phi` is the latitude. When :math:`\cos(\theta_z) \le 0`, the sun is below the horizon 
and the direct-beam contribution is zero.

References
==========

Beer, A. (1852). Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten. *Annalen der Physik und Chemie*, 86(5), 78–88.

Bhumralkar, C. M. (1974). Numerical experiments on the computation of ground surface temperature in an atmospheric general circulation model. *Journal of Applied Meteorology*, 13(7), 697–704.

Kirchhoff, G. R. (1860). Über die Beziehung zwischen den Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung. *Annalen der Physik und Chemie*, 109(3), 275–301.

Meador, W. E., & Weaver, W. R. (1980). Two-stream approximations to radiative transfer in planetary atmospheres: A unified description of existing methods and a new improvement. *Journal of the Atmospheric Sciences*, 37(3), 630–643.

Toon, O. B., McKay, C. P., Ackerman, T. P., & Santhanam, K. (1989). Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres. *Journal of Geophysical Research*, 94(D13), 16465–16481.

Tremback, C. J., & Kessler, R. C. (1985). A surface temperature and moisture parameterization scheme for use in mesoscale models. *Journal of the Atmospheric Sciences*, 42(21), 2751–2761.
