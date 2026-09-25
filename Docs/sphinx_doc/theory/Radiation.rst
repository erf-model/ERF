
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
attenuation for shortwave radiation, a Meador-Weaver two-stream diffuse field combined by the adding method
with the surface albedo, gray-gas longwave two-stream, and an optional Simplified Surface Energy Balance (SEB) module.

Shortwave Radiation
--------------------------------------

The shortwave (solar) radiation calculation in the two-stream model is split into direct-beam
and diffuse components. The shortwave direct-beam radiation is attenuated through the atmosphere according to the
Beer-Lambert law:

.. math::

   I_{sw,\text{direct}}(z) = I_0 \mu_0 e^{-\tau_{\text{sw}} \sec(\theta_z)}

where :math:`I_0` is the top-of-atmosphere irradiance (``erf.fixed_total_solar_irradiance``, or
1360.9 W/m² scaled by the Earth-Sun distance factor of the date from the same orbital code RRTMGP
uses), :math:`\mu_0 = \cos(\theta_z)` is the cosine of the solar zenith angle
(``erf.fixed_solar_zenith_angle``, or the sun's position over each column at the calendar time
given by ``start_datetime``), :math:`\tau_{\text{sw}}` is the
vertically integrated shortwave optical depth above height :math:`z`, and :math:`\sec(\theta_z)` accounts
for the path-length modification. The optical depth may be spatially uniform (static :math:`\tau_{\text{per\_layer}}`)
or dynamically diagnosed from water vapor and cloud liquid water content, parameterized as:

.. math::

   \tau_{\text{sw}}(k) = \tau_{\text{per\_layer}} + \tau_{\text{cloud}}(k) + c_{\text{qv}} q_v(k) + c_{\text{qc}} q_c(k)

where :math:`\tau_{\text{cloud}}(k)` is added only within the prescribed cloud layer, and the coefficients
:math:`c_{\text{qv}}` and :math:`c_{\text{qc}}` are zero by default.



This per-layer model (``tau_model = per_layer``, the default) assigns the same optical depth to every
layer regardless of its thickness, so the column optical depth scales with the number of vertical cells.
The mass model (``tau_model = mass``) instead builds each layer's optical properties from its mass path
:math:`\rho \, \Delta z`:

.. math::

   \tau_{\text{sw}}(k) = \rho \Delta z \left( k_{\text{abs,dry}} + k_{\text{sca,dry}} + k_{\text{abs,v}} q_v + k_{\text{ext,c}} q_c \right),

with the constituents mixed by extinction weighting into the layer single-scattering albedo and asymmetry
factor, :math:`\omega_0 = \sum_i \omega_i \tau_i / \tau` and :math:`g = \sum_i g_i \omega_i \tau_i / \sum_i \omega_i \tau_i`
(dry and vapor absorption with :math:`\omega = 0`, Rayleigh scattering with :math:`\omega = 1, g = 0`, cloud water
with ``sw_cloud_omega`` and ``sw_cloud_g``; ``sw_kext_cloud`` :math:`\approx 1.5 / r_{\text{eff}}` with
:math:`r_{\text{eff}}` in µm gives 150 m²/kg for 10 µm droplets). The prescribed cloud band, the moisture
coefficients and the aerosol term are added on top as before. The column optical depth is then a property of
the atmosphere, not of the grid, and the longwave band uses the mass path described below.
The diffuse (scattered) shortwave field has upward and downward streams and is solved with the
two-stream approximation. Each layer :math:`k` with optical depth :math:`\tau`, single-scattering albedo
:math:`\omega_0` and asymmetry factor :math:`g`,

.. math::

   \omega_0 = \frac{\text{scattering cross-section}}{\text{total extinction cross-section}}, \quad
   g = \left\langle \cos(\theta) \right\rangle_{\text{scattering}},

is characterized by its reflectance and transmittance for diffuse incidence, :math:`R_{\text{dif}}` and
:math:`T_{\text{dif}}`, and by the diffuse flux it reflects upward and transmits downward per unit direct-beam
flux incident at its top, :math:`R_{\text{dir}}` and :math:`T_{\text{dir}}` (the surviving direct beam is
:math:`T_{\text{ns}} = e^{-\tau/\mu_0}`). With the practical-improved-flux-method coefficients
(Zdunkowski et al. 1980)

.. math::

   \gamma_1 = \frac{8 - \omega_0 (5 + 3g)}{4}, \quad
   \gamma_2 = \frac{3 \omega_0 (1 - g)}{4}, \quad
   \gamma_3 = \frac{2 - 3 g \mu_0}{4}, \quad
   \gamma_4 = 1 - \gamma_3, \quad
   k = \sqrt{\gamma_1^2 - \gamma_2^2},

the Meador and Weaver (1980) layer solution reads

.. math::

   R_{\text{dif}} = \frac{\gamma_2 (1 - e^{-2 k \tau})}{D}, \qquad
   T_{\text{dif}} = \frac{2 k e^{-k \tau}}{D}, \qquad
   D = k (1 + e^{-2 k \tau}) + \gamma_1 (1 - e^{-2 k \tau}),

and, with :math:`\alpha_1 = \gamma_1 \gamma_4 + \gamma_2 \gamma_3` and
:math:`\alpha_2 = \gamma_1 \gamma_3 + \gamma_2 \gamma_4`,

.. math::

   R_{\text{dir}} = \frac{\omega_0}{D (1 - k^2 \mu_0^2)} \left[ (1 - k\mu_0)(\alpha_2 + k\gamma_3)
      - (1 + k\mu_0)(\alpha_2 - k\gamma_3) e^{-2k\tau} - 2 (k\gamma_3 - \alpha_2 k \mu_0) e^{-k\tau} T_{\text{ns}} \right],

   T_{\text{dir}} = -\frac{\omega_0}{D (1 - k^2 \mu_0^2)} \left[ (1 + k\mu_0)(\alpha_1 + k\gamma_4) T_{\text{ns}}
      - (1 - k\mu_0)(\alpha_1 - k\gamma_4) e^{-2k\tau} T_{\text{ns}} - 2 (k\gamma_4 + \alpha_1 k \mu_0) e^{-k\tau} \right].

A non-scattering layer (:math:`\omega_0 = 0`) has :math:`R_{\text{dif}} = R_{\text{dir}} = T_{\text{dir}} = 0`
and :math:`T_{\text{dif}} = e^{-2\tau}`; a conservative layer (:math:`\omega_0 = 1`) satisfies
:math:`R_{\text{dif}} + T_{\text{dif}} = 1` and :math:`R_{\text{dir}} + T_{\text{dir}} + T_{\text{ns}} = 1`.

The layers are combined with the surface by the adding method. Let :math:`A_m` be the albedo of everything
below interface :math:`m` for diffuse light and :math:`S_m` the upward diffuse flux at interface :math:`m`
produced by the direct beam illuminating everything below it. The surface, with direct-beam albedo
:math:`\alpha_{\text{dir}}` (user-specified or LSM-provided) and diffuse albedo :math:`\alpha_{\text{dif}}`
(``surface_albedo_sw_diffuse``, equal to :math:`\alpha_{\text{dir}}` unless set), starts the recursion,

.. math::

   A_0 = \alpha_{\text{dif}}, \qquad S_0 = \alpha_{\text{dir}} F_{\text{dir}}(0),

and each layer :math:`m` (between interfaces :math:`m` and :math:`m+1`) adds

.. math::

   A_{m+1} = R_{\text{dif}} + \frac{T_{\text{dif}}^2 A_m}{1 - R_{\text{dif}} A_m}, \qquad
   S_{m+1} = R_{\text{dir}} F_{\text{dir}}(m+1) + \frac{T_{\text{dif}} \left[ S_m + A_m T_{\text{dir}} F_{\text{dir}}(m+1) \right]}{1 - R_{\text{dif}} A_m}.

With no diffuse flux incident at the top of the atmosphere, the downward pass gives the diffuse downward and
upward fluxes on every interface,

.. math::

   F^{\downarrow}_{\text{dif}}(m) = \frac{T_{\text{dif}} F^{\downarrow}_{\text{dif}}(m+1) + T_{\text{dir}} F_{\text{dir}}(m+1) + R_{\text{dif}} S_m}{1 - R_{\text{dif}} A_m}, \qquad
   F^{\uparrow}_{\text{dif}}(m) = A_m F^{\downarrow}_{\text{dif}}(m) + S_m,

and the net shortwave flux :math:`F_{\text{sw,net}}(m) = F_{\text{dir}}(m) + F^{\downarrow}_{\text{dif}}(m) - F^{\uparrow}_{\text{dif}}(m)`
whose divergence drives the shortwave heating rate. The surface absorbs
:math:`(1 - \alpha_{\text{dir}}) F_{\text{dir}}(0) + (1 - \alpha_{\text{dif}}) F^{\downarrow}_{\text{dif}}(0)`, which is the
``SW_surface`` diagnostic; the reflected flux leaving the top, :math:`F^{\uparrow}_{\text{dif}}(n)`, is ``SW_up_TOA``.

Cloud scattering properties may differ from the clear-sky values (e.g., :math:`\omega_0^{\text{cloud}}` for
liquid water clouds). The cloud/clear-sky distinction is blended according to the cloud fraction :math:`C_f`:

.. math::

   F_{\text{sw,net}} = (1 - C_f) F_{\text{sw,net}}^{\text{clear}} + C_f F_{\text{sw,net}}^{\text{cloudy}}

Longwave Radiation
--------------------------------------

The longwave (thermal) radiation employs a gray-gas two-stream formulation (Toon et al. 1989)
to compute upward and downward fluxes. Assume local thermodynamic equilibrium (LTE): each layer emits radiation according to the
Planck function :math:`B(T)` weighted by the emissivity of the gray gas. The optical depth is
parameterized similarly to shortwave:

.. math::

   \tau_{\text{lw}}(k) = \tau_{\text{lw,per\_layer}} + \tau_{\text{cloud,lw}}(k) + c_{\text{qv,lw}} q_v(k) + c_{\text{qc,lw}} q_c(k)

Alternatively (``lw_mass_absorption_enable = true``) the gray optical depth follows the mass path of
each layer,

.. math::

   \tau_{\text{lw}}(k) = \rho \, \Delta z \left( k_{\text{dry}} + k_{\text{v}} q_v + k_{\text{c}} q_c \right),

so the column optical depth is independent of the vertical resolution, water vapor and cloud water
have a real greenhouse effect, and the cloud term reproduces the Stephens (1978) emissivity
:math:`\epsilon_c = 1 - e^{-0.158 \, \text{LWP}}` (LWP in g/m²). The cloud-band, moisture-coefficient
and aerosol additions of the previous equation still apply on top of this base.

The emission temperature of each layer is the absolute temperature, recovered from the prognostic
:math:`\rho\theta` through the equation of state and the Exner function,

.. math::

   p = p_0 \left( \frac{R_d \, \rho \, \theta_m}{p_0} \right)^{\gamma}, \qquad
   T = \theta \left( \frac{p}{p_0} \right)^{R_d / c_p},

where :math:`\theta_m = \theta (1 + R_v q_v / R_d)` is the moist potential temperature. The column
sweeps follow ERF's vertical index convention: the lowest cell-centered index is the layer adjacent to
the surface and the highest index is the layer adjacent to the top of the domain, with fluxes stored on
the layer interfaces between them. On a stretched or terrain-following grid the thickness of a layer is
the distance between its two interfaces, taken from the nodal heights ``z_phys_nd`` (each interface
height over a column is the mean of its four nodes); that thickness enters both the heating-rate
divergence :math:`\Delta F / (\rho c_p \Delta z)` and the mass path :math:`\rho \Delta z` of the
mass optical-depth model. On a uniform grid it is the cell size.

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

Both the shortwave and the longwave heating rates are temperature tendencies. ERF advances
:math:`\rho\theta`, so the two-stream model stores the corresponding potential-temperature tendency,

.. math::

   \left.\frac{\partial \theta}{\partial t}\right|_{\text{rad}} = \frac{H_{\text{sw}} + H_{\text{lw}}}{\pi},
   \qquad \pi = \left( \frac{p}{p_0} \right)^{R_d / c_p},

in the ``qheating_rates`` array (and the ``qsrc_sw`` / ``qsrc_lw`` plot variables), which the
:math:`\rho\theta` source term multiplies by :math:`\rho`. This is the same convention as the RRTMGP path.

Grid Requirements
--------------------------------------

Both sweeps integrate a whole atmospheric column in one pass, from the surface at the lowest
:math:`k` to the top of the atmosphere at the highest. Each grid must therefore span the domain
in the vertical: a box that holds only part of a column would see neither the beam arriving from
above nor the cooling to space, and would return heating rates that look plausible and are wrong.
ERF decomposes the base grid in :math:`z` only when ``amr.max_grid_size_z`` is smaller than the
number of cells in :math:`z`, so setting

.. code-block:: none

   amr.max_grid_size_z = <at least amr.n_cell in z>

is enough to satisfy this. The model aborts with a message naming this input if it is given a
vertically decomposed grid. The horizontal decomposition is unconstrained, and the results do not
depend on it or on the ``fabarray.mfiter_tile_size`` tiling.

On a refined run a level is free to cover only part of the column, and that is supported: a
level whose grids do not reach the domain top or bottom is a *nested patch*, and rather than
sweeping it ERF interpolates its heating rates and fluxes from its parent -- the same route
RRTMGP takes (``is_nested_patch``). Nothing needs to be set for this.

The requirement is per box, not per level: the sweep needs a whole column inside one box. Several
layouts fail it -- a level that stops short of the domain top, a level tagged at different heights
in different horizontal regions (surface convection in one place, cloud tops in another), or grids
decomposed in the vertical -- and above level 0 they all take the same route, interpolation from
the parent. None of them is an error.

Level 0 is the exception, because it has no parent. It always covers the domain, so a box there
that does not span :math:`z` is decomposed in the vertical, and that is refused at start-up.
ERF's default ``amr.no_box_split_dir = 2`` already forbids that decomposition, so the refusal is
a backstop rather than something a normal deck meets.

If you would rather a refinement patch be solved on its own than interpolated, setting

.. code-block:: none

   amr.refine_whole_domain_dir = 2

makes AMReX cluster the tagged cells in the horizontal only and emit refinement boxes that span
the whole domain in :math:`z`, so every level carries complete columns and every level runs its
own sweep. A refinement box given explicitly through ``erf.boxN.in_box_lo``/``in_box_hi`` spans
:math:`z` already when the :math:`z` extent is omitted, since the two-value form defaults to the
full domain.

Multiple Levels
--------------------------------------

Every level that carries complete columns runs its own column sweep over its own state, terrain
and surface properties, and writes its own heating rates into ``qheating_rates[lev]``; a nested
patch is interpolated from its parent instead. The RhoTheta source applies them at every level. There is no coarse-fine treatment of the radiative fluxes and none is needed in the
usual sense -- radiation is a source term, not a conserved flux that is refluxed -- but two
consequences follow and are worth stating plainly:

- **A lateral seam.** Across the edge of a patch, the coarse and the fine solution of the same
  physical column differ slightly, because they are computed on different grids. For a smooth
  broadband two-stream model the difference is small, but nothing smooths it. Under
  ``erf.coupling_type = TwoWay`` (the default) coarse cells underneath a patch have their state
  replaced by the fine solution at the end of each step (``AverageDown``), so the discrepancy
  does not accumulate there; under ``OneWay`` there is no such replacement and it does.
- **No feedback upward.** The fine level's own structure does not influence the coarse level's
  radiation.

A subcycled fine level calls radiation once per level step, so it runs ``nsubsteps[lev]`` times
as often as its parent -- twice as often for a refinement ratio of two. This matches the RRTMGP
path and is physically correct, since the heating is recomputed from the current old state each
time; there is no call-interval input to reduce it.

The surface energy balance runs on level 0 only. The surface is one physical object whose
force-restore state is prognostic and checkpointed, so it has a single owner, and only level 0
writes or reads that state in a checkpoint. The diagnostic SEB residual is likewise reported for
level 0. Because the prognostic surface temperature *is* the longwave boundary condition,
combining ``erf.radiation.seb_prognostic_enable`` with ``amr.max_level > 0`` would leave level 0
and its fine levels with two different surface boundary conditions for one surface; that
combination is refused at start-up rather than allowed to disagree silently.

Limitations
--------------------------------------

- **Refined runs.** Multiple levels are supported; see `Multiple Levels`_ above for the grid
  requirement, the lateral coarse-fine seam, the absence of feedback from fine to coarse, the
  subcycled call cadence, and the single-level restrictions on the surface energy balance.
- **Sun and site.** The sun, the site and the surface temperature come from the inputs the
  RRTMGP interface reads (``erf.fixed_solar_zenith_angle``, ``erf.fixed_total_solar_irradiance``,
  ``erf.rad_t_sfc``, ``erf.rad_cons_lat``/``lon``, ``erf.rad_orbital_*``, ``start_datetime``),
  and the position of the sun is the instantaneous value of RRTMGP's ``orbital_cos_zenith``
  over each column, not its average over the radiation interval. Without a fixed zenith angle
  and irradiance the run needs ``start_datetime``, and stops at the first sweep otherwise. The
  surface temperature of the longwave boundary is, in RRTMGP's order, the land-surface model's
  field, else the surface layer's potential temperature converted to absolute temperature, else
  ``erf.rad_t_sfc``; the prognostic surface
  energy balance, when on, supplies its own state ahead of the surface layer. The surface
  layer works in potential temperature, so its value is converted to temperature with the
  Exner function evaluated at the physical surface pressure diagnosed from the lowest atmospheric cell
  before it enters the :math:`\sigma T_s^4` emission; with a
  surface layer present, ``erf.rad_t_sfc`` is the initial value of the prognostic surface
  temperature when the surface energy balance evolves one, and unused otherwise.
- **Call cadence.** The sweep runs once every slow step, from the old state, and there is no
  call-frequency input: one gray sweep per column costs about a millisecond per ten thousand
  cells, so unlike RRTMGP (``erf.rad_freq_in_steps``) it is not worth skipping steps.
- **Diagnostics file.** The diagnostics are off by default. Setting
  ``erf.radiation.diag_enable = true`` writes ``radiation_diag.dat``
  (``erf.radiation.diag_file``) in the run directory, with a ``pre_dycore`` and a
  ``post_dycore`` row per step for every level that sweeps. Each of them appends to the one
  file and the last column, ``level``, tells the rows apart. A level interpolated from its
  parent contributes no rows at all -- it runs no sweep, so it has no fluxes of its own to
  report -- so on a refined run the rows present are those of the sweeping levels, not one
  set per level in the hierarchy. The file is appended to rather than truncated, as ERF's other
  data logs are, so a rerun in the same directory extends the previous run's rows.

Surface Energy Balance
--------------------------------------

The net surface shortwave and longwave fluxes come from the land-surface model when it exposes
them (Noah-MP's absorbed shortwave ``sav + sag`` and, with the sign flipped to absorbed, its net
longwave ``fira``); otherwise, with ``erf.radiation.seb_use_radiation_fluxes = true``, from the
two-stream sweep's own surface fluxes in every column; otherwise from the scalar
``seb_sw_flux_default`` and ``seb_lw_flux_default``. The sensible, latent and ground heat fluxes
and the deep-soil reservoir values are the scalar defaults unless the land-surface model exposes
them by name (``grdflx`` for the ground heat flux).

The surface energy balance residual is defined as the net radiative flux minus the turbulent and
ground heat fluxes:

.. math::

   R_{\text{net}} = F_{\text{sw,net}}(0^+) + F_{\text{lw,net}}(0^+)

   \text{SEB}_{\text{residual}} = R_{\text{net}} - H - \text{LE} - G

where :math:`H` is the sensible heat flux, :math:`\text{LE}` is the latent heat flux (evaporation),
and :math:`G` is the ground heat flux conducted into the soil.


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
--------------------------------


When ``cloud_fraction_prog_enable = true``, the cloud fraction is diagnosed at each level from the
relative humidity (RH) and cloud liquid water content (qc):

.. math::

   C_f = \min \left( 1, \max \left( 0, \frac{RH - \text{rh_min}}{\text{rh_max} - \text{rh_min}} \right) + \min \left( 1, \frac{q_c}{q_{c,\text{scale}}} \right) \right)

where the RH threshold parameters allow for a transition from clear sky (RH < rh_min) to complete cloud coverage
(RH >= rh_max). The coefficient :math:`c_{\text{qc}}` (default :math:`1 \times 10^{-3}`) provides an additional
scaling of liquid water content's contribution. Optional temporal smoothing via an exponential moving average
(EMA) may be applied to suppress oscillations:

.. math::

   C_f^{\text{smoothed}} = \alpha C_f^{\text{new}} + (1 - \alpha) C_f^{\text{old}}

where :math:`\alpha` is the blending parameter.

Solar Geometry and Diurnal Cycle
--------------------------------

The position of the sun comes from the same orbital code the RRTMGP interface uses
(``ERF_OrbCosZenith.H``). When ``erf.fixed_solar_zenith_angle`` is not set, the calendar time of
each call is ``start_datetime`` plus the simulation time (UTC). The year gives the orbital
parameters of Berger (1978) unless ``erf.rad_orbital_year``, ``erf.rad_orbital_eccentricity``,
``erf.rad_orbital_obliquity`` or ``erf.rad_orbital_mvelp`` override them, and the day of the year
:math:`D` (with its fraction, leap years included) gives the solar declination :math:`\delta` and
the Earth-Sun distance factor :math:`e` from the orbital elements. The cosine of the solar zenith
angle over a column at latitude :math:`\phi` and east longitude :math:`\lambda` is the
instantaneous value

.. math::

   \cos(\theta_z) = \sin(\phi) \sin(\delta) - \cos(\phi) \cos(\delta) \cos\left(2\pi\, f + \lambda\right)

where :math:`f` is the fraction of the UTC day elapsed, so that solar noon at a site falls at
:math:`2\pi f + \lambda = \pi`. The latitude and longitude are the grid's own ``lat_m`` and
``lon_m`` fields when a WRF or metgrid initialisation filled them, and ``erf.rad_cons_lat`` and
``erf.rad_cons_lon`` otherwise. When :math:`\cos(\theta_z) \le 0` the sun is below the horizon and
the direct-beam contribution is zero. The top-of-atmosphere irradiance is
``erf.fixed_total_solar_irradiance`` when set, else :math:`1360.9\, e` W/m² (RRTMGP's reference
value scaled by the distance factor of the date). RRTMGP averages :math:`\cos(\theta_z)` over its
radiation interval; the two-stream model runs every step and uses the instantaneous value.

With ``erf.fixed_solar_zenith_angle`` set (a cosine, applied to every column) no calendar is
needed: the irradiance is then ``erf.fixed_total_solar_irradiance``, or, when that is not set
either, :math:`1360.9\, e` W/m² if ``start_datetime`` is known and the unscaled 1360.9 W/m² if
not.

References
--------------------------------------

Beer, A. (1852). Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten. *Annalen der Physik und Chemie*, 86(5), 78–88.

Bhumralkar, C. M. (1974). Numerical experiments on the computation of ground surface temperature in an atmospheric general circulation model. *Journal of Applied Meteorology*, 13(7), 697–704.

Kirchhoff, G. R. (1860). Über die Beziehung zwischen den Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung. *Annalen der Physik und Chemie*, 109(3), 275–301.

Meador, W. E., & Weaver, W. R. (1980). Two-stream approximations to radiative transfer in planetary atmospheres: A unified description of existing methods and a new improvement. *Journal of the Atmospheric Sciences*, 37(3), 630–643.

Stephens, G. L. (1978). Radiation profiles in extended water clouds. II: Parameterization schemes. *Journal of the Atmospheric Sciences*, 35, 2123–2132.

Zdunkowski, W. G., Welch, R. M., & Korb, G. (1980). An investigation of the structure of typical two-stream methods for the calculation of solar fluxes and heating rates in clouds. *Beiträge zur Physik der Atmosphäre*, 53, 147–166.

Toon, O. B., McKay, C. P., Ackerman, T. P., & Santhanam, K. (1989). Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres. *Journal of Geophysical Research*, 94(D13), 16465–16481.

Tremback, C. J., & Kessler, R. C. (1985). A surface temperature and moisture parameterization scheme for use in mesoscale models. *Journal of the Atmospheric Sciences*, 42(21), 2751–2761.
