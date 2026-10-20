
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

The surface energy balance runs on every level. Each level carries its own force-restore
surface state, evolves it from the fluxes its own sweep computes, and checkpoints it, and the
diagnostic SEB residual is reported per level. Because the prognostic surface temperature *is*
the longwave boundary condition, the levels must not be allowed to hold different temperatures
for one physical surface, so ERF keeps them consistent in three places:

- **A new level starts from its parent.** ``t_sfc`` and ``q_sfc`` are interpolated from the
  coarse level when a level is created, so a fine level begins from the surface its parent has
  already reached rather than from ``erf.rad_t_sfc``.
- **Fine levels are averaged down.** After the finer levels advance, their surface state is
  averaged onto the coarse level, exactly as the atmospheric state is. This runs under
  ``erf.coupling_type = TwoWay``; with ``OneWay`` the levels are left to evolve independently,
  which is what that option asks for everywhere else as well.
  A level that cannot sweep is skipped: one whose boxes do not span the domain in
  :math:`z` takes its radiation fields from its parent and never advances a surface state
  of its own, so averaging its copy down would overwrite the coarse surface underneath the
  patch with the value it was created with.
- **A regrid keeps what the surface had reached.** Rebuilding a level reallocates its surface
  fields, so the pre-regrid values are copied back onto the new grids, with cells the new grids
  added filled from the parent.

The surface temperature a fine level sees is therefore its parent's wherever the fine level has
not yet changed it, which is a real limitation: refining does not by itself give the surface
more structure than the coarse grid resolved. What refinement does give is a surface that
responds to the fine level's own radiative fluxes.

Limitations
--------------------------------------

- **Refined runs.** Multiple levels are supported; see `Multiple Levels`_ above for the grid
  requirement, the lateral coarse-fine seam, the absence of feedback from fine to coarse, the
  subcycled call cadence, and how the surface energy balance is kept consistent between levels.
  A run with Noah-MP is refined as well: every level hands its own Noah-MP its forcing, and a
  finer level runs the land model on a land setup file of its own, or else takes its land state
  from level 0; see `Radiative Forcing of a Land-Surface Model`_.
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
``seb_sw_flux_default`` and ``seb_lw_flux_default``.

The sensible heat flux :math:`H` and latent heat flux :math:`\text{LE}` come, in order of
precedence, from a land-surface model field of that name (``hfx``, ``lh``; no land model exposes
one today); otherwise, with ``erf.radiation.seb_turbulent_flux_source = surface_layer`` (the
default), from the fluxes the ``zlo`` surface layer applies to the air,

.. math::

   H = c_p \, \overline{\rho w'\theta'}\big|_{\text{sfc}}, \qquad
   \text{LE} = L_v \, \overline{\rho w' q_v'}\big|_{\text{sfc}},

positive away from the surface -- the same conversion as the ``sensible_heat_flux`` and
``latent_heat_flux`` 2D outputs, so the ground loses exactly what those report the air gaining;
otherwise from the scalar ``seb_hfx_default`` and ``seb_lh_default``. The surface layer is the
source wherever its flux field exists, that is with any diffusion or turbulence closure; the
defaults apply with ``seb_turbulent_flux_source = defaults``, without a ``zlo`` surface layer, on EB
terrain (where the surface layer's flux goes to the embedded boundary instead), without diffusion
or a closure, and for :math:`\text{LE}` without a moisture model. An adiabatic surface layer
(``erf.most.surf_temp_flux = 0``) has the field and a zero flux, so :math:`H = 0` there rather than
``seb_hfx_default``; ERF warns at start-up when a nonzero default is replaced this way. With
``erf.use_rotate_surface_flux`` the surface layer splits its flux over the three faces of the lowest
cell, and the balance, like the 2D outputs, removes only the vertical-face part, :math:`\cos`
of the slope times :math:`H` and :math:`\text{LE}`; ERF warns at start-up. The fluxes the balance used are written as the ``seb_hfx`` and ``seb_lh`` 2D
plotfile variables.

By default the surface layer computes these fluxes from its own surface temperature and moisture
(``erf.most.surf_temp`` and the like), so the coupling runs one way: the balance loses what the
surface layer puts into the air, but the flux does not respond to :math:`T_s`. With
``erf.radiation.seb_surface_layer_uses_skin = true`` it runs both ways. Before it computes its
fluxes each step, the surface layer sets its land surface temperature to the skin the balance
reached at the end of the previous step, as a potential temperature,

.. math::

   \theta_s = T_s \left( \frac{p_0}{p_{\text{sfc}}} \right)^{R_d/c_p},

with the surface pressure :math:`p_{\text{sfc}}` diagnosed from the lowest cell (the conversion
coupled sea-surface temperatures use). A warmer skin then gives a larger :math:`H`, which the
balance removes. The surface moisture stays the surface layer's own unless
``erf.radiation.seb_surface_layer_uses_moisture`` is also set (below). The option needs the prognostic balance,
``seb_turbulent_flux_source = surface_layer``, a ``zlo`` surface layer in surface-temperature mode
(``erf.most.surf_temp`` given, no ``erf.most.surf_heating_rate``), no EB terrain, no
``erf.use_rotate_surface_flux``, and no land-surface or surface model; ERF stops at start-up
otherwise. On a level that takes its
radiation from its parent (a nested patch that does not span the column), no skin evolves, and
the surface layer keeps its own temperature there.

**Soil moisture.** With the prognostic balance, :math:`q_s` is the volumetric water content
[m\ :sup:`3`/m\ :sup:`3`] of the top ``seb_moisture_layer_depth_m`` :math:`d_s` of soil: the latent
heat flux drains it, and it restores to ``seb_q_deep_default`` over
``seb_moisture_restore_timescale_s`` :math:`\tau_q`,

.. math::

   \frac{dq_s}{dt} = -\frac{\text{LE}}{L_v \rho_w d_s} - \frac{q_s - q_\text{deep}}{\tau_q}.

With ``erf.radiation.seb_surface_layer_uses_moisture = true`` (which needs the skin coupling) the
surface layer takes its land surface mixing ratio from that soil water,

.. math::

   q_\text{surf} = \beta\, q_\text{sat}(T_s, p_\text{sfc}) + (1 - \beta)\, q_\text{air},

with :math:`q_\text{air}` the mixing ratio at its reference height, so its moisture flux is
:math:`\beta` times the potential one and the soil loses the water the air gains. Without a
soil type :math:`\beta` is the soil-water factor

.. math::

   \beta_\text{soil} = \min\left(1, \max\left(0,
       \frac{q_s - \theta_\text{wilt}}{\theta_\text{fc} - \theta_\text{wilt}}\right)\right),

with the wilting point and field capacity ``seb_soil_moisture_wilt`` and ``seb_soil_moisture_fc``.
``erf.radiation.seb_soil_type`` takes both from Noah-MP's soil table for that category (WLTSMC
and REFSMC of the STAS dataset; ``Source/Radiation/TwoStream/ERF_NoahMPSoilTable.H`` copies the
table and a CTest checks the copy against ``Submodules/Noah-MP/parameters/NoahmpTable.TBL``).

With a soil type, the surface mixing ratio comes instead from two source fluxes, each through
a resistance in series with the aerodynamic one. That one is the resistance of the surface
layer's own moisture flux (its surface-temperature kernel),
:math:`r_a = \max(\ln(z_\text{ref}/z_0) - \psi_h, 1)/(\kappa u_*)`, with Jimenez's
:math:`\psi_h` at the layer's last :math:`u_*` and Obukhov length (neutral before the first
flux):

.. math::

   q_\text{surf} = q_\text{air} + r_a \left[ f_\text{veg} \frac{q_\text{sat} - q_\text{air}}{r_a + r_c}
     + (1 - f_\text{veg}) \frac{\text{RH}_g\, q_\text{sat} - q_\text{air}}{r_a + r_\text{soil}} \right],

so that the surface layer's flux :math:`(q_\text{surf} - q_\text{air})/r_a` is the sum of the
canopy's and the bare soil's. :math:`f_\text{veg} = 0` without a vegetation type, so bare soil and
a vegetation type with ``seb_vegetation_fraction = 0`` give the same :math:`q_\text{surf}`. With
:math:`\text{RH}_g = 1` it is :math:`\beta q_\text{sat} + (1 - \beta) q_\text{air}` with
:math:`\beta = f_\text{veg} r_a/(r_a + r_c) + (1 - f_\text{veg}) r_a/(r_a + r_\text{soil})`. With
``erf.radiation.seb_vegetation_type`` (a Noah-MP MODIS land-use category,
``ERF_NoahMPVegetationTable.H``), the vegetated fraction ``seb_vegetation_fraction``
:math:`f_\text{veg}` transpires through Noah's
big-leaf Jarvis canopy resistance (Chen et al. 1996) on the parameters and floors of Noah-MP's
canopy-resistance option 2, with the category's leaf area index (``seb_leaf_area_index``; by
default the table's monthly values interpolated to ``start_datetime`` as Noah-MP does, with the
day of the year counted from 0 at 00:00 on 1 January and shifted half a year when
``erf.rad_cons_lat`` < 0):

.. math::

   r_c = \frac{R_{s,\min}}{\text{LAI}\, F_{sw} F_T F_\text{vpd} \beta_\text{soil}}, \quad
   F_{sw} = \frac{f + R_{s,\min}/R_{s,\max}}{1 + f}, \quad f = \frac{1.1\, SW_\downarrow}{R_{gl}\,\text{LAI}},

   F_T = 1 - 0.0016\, (T_\text{opt} - T_s)^2, \quad
   F_\text{vpd} = \frac{1}{1 + h_s \max(0, q_\text{sat} - q_\text{air})},

with :math:`SW_\downarrow` the net shortwave the balance holds divided by
:math:`1 - \alpha` (the sweep's with ``seb_use_radiation_fluxes``). :math:`r_c` is capped at
:math:`10^6` s/m after the vapour-deficit factor. Noah-MP applies the same factors per sunlit
and shaded leaf, with absorbed PAR, the canopy temperature and canopy-air humidity and a
root-zone soil-water factor, so the two agree in form, not in every detail. The bare soil (the whole surface without a vegetation type)
evaporates through Noah-MP's soil resistance (ground-evaporation option 1, Sakaguchi and Zeng),
which grows as the top soil dries: :math:`r_\text{soil} = d_\text{dry}/D` with
:math:`d_\text{dry} = d_s (e^{(1 - q_s/\theta_\text{sat})^5} - 1)/(e - 1)` and
:math:`D = 2.2\times10^{-5}\, \theta_\text{sat}^2 (1 - \theta_\text{wilt}/\theta_\text{sat})^{2 + 3/b}`.
It evaporates from the air in its pores, at Noah-MP's relative humidity

.. math::

   \text{RH}_g = \exp\left(\frac{\psi g}{R_v T_s}\right), \quad
   \psi = -\psi_\text{sat} \left(\frac{\max(0.01, q_s)}{\theta_\text{sat}}\right)^{-b},

with the soil's saturated matric potential :math:`\psi_\text{sat}` (SATPSI) and Noah-MP's
:math:`g` and :math:`R_v`. :math:`\text{RH}_g` is near 1 in moist soil and near 0 at the wilting
point (about 0.005 for silty clay loam at 325 K), so bare soil there hardly evaporates and,
when :math:`\text{RH}_g q_\text{sat} < q_\text{air}`, takes up vapour (LE < 0), as in Noah-MP.

*Limitation:* :math:`q_s` is a single top layer. The canopy's soil-water factor and all of LE,
transpiration included, act on it, where Noah-MP draws transpiration from the root zone
(BTRAN). Over a few hours this hardly matters: at LE = 350 W/m\ :sup:`2` and
:math:`d_s` = 0.1 m the layer loses about 0.005 m\ :sup:`3`/m\ :sup:`3` per hour. Over several
days, though, the canopy shuts down as the top layer dries, even over a wet root zone.
:math:`q_s` also restores toward ``seb_q_deep_default``, which defaults to 0; ERF prints a
NOTE at start-up when it is below the wilting point.

With the skin coupling and a soil type, the surface layer's land roughness length, unless
``erf.most.z0`` is given, comes from Noah-MP's tables: :math:`f_\text{veg}\, z_{0,\text{veg}} + (1 - f_\text{veg})\, z_{0,\text{soil}}`
with the land-use category's Z0MVT and Noah-MP's bare-soil Z0SOIL (0.002 m), :math:`f_\text{veg} = 0`
without a vegetation type.

Noah-MP itself hands the atmosphere Z0MVT for a vegetated column and Z0SOIL for a bare one.
The balance has a single surface for both parts, so it weights the two by the vegetated
fraction instead. The run's ``job_info`` records the value used as ``erf.most.z0``. Heat and
moisture use the same :math:`z_0` as momentum (the surface layer's kernel has no separate
:math:`z_{0h}`), where Noah-MP's bare-ground exchange takes a smaller thermal roughness
(Chen-Zilitinkevich): a further difference over bare soil.

Categories 15, 16 and 17 (snow and ice, barren, water) have no vegetation in the table
(Z0MVT = 0) and stop at start-up. Bare land leaves the vegetation type at 0.

Over bare soil the roughness matters as much as the moisture: with the surface layer's
default 0.1 m, the surface sheds its heat far more easily than Noah-MP's bare soil does.

With a soil type the surface heat capacity :math:`C_s`, unless given, is that of the layer the
restore period's temperature wave reaches (Deardorff 1978),
:math:`C_s = \tfrac12 \sqrt{\lambda c\, \tau / \pi}`, with the soil's volumetric heat capacity
:math:`c` and Noah-MP's (Johansen) thermal conductivity :math:`\lambda` at ``seb_q_sfc_default``.
``Exec/CanonicalTests/Radiation/TwoStream_NoahMP_vs_ForceRestore`` compares the balance with
Noah-MP on the same grid, atmosphere and land.

The ground heat flux :math:`G` and the deep-soil reservoir values are the scalar defaults unless
the land-surface model exposes them by name (``grdflx`` for the ground heat flux).

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

In the force-restore method the restoring term is the heat conducted into the soil, so it already
plays the part of :math:`G`. Leave ``seb_grdflx_default`` at 0 with the prognostic mode: a nonzero
value removes that heat a second time, and ERF prints a warning when it is set.

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

External Surface-Temperature Provider Ownership
-------------------------------------------------

TwoStream advances its prognostic surface state only when TwoStream owns the longwave surface-temperature
boundary. If an authoritative external or LSM surface-temperature provider owns that boundary at a level, the
TwoStream prognostic update is skipped there, preventing an unused shadow state from being evolved alongside
the provider. Noah-MP's ``t_sfc`` field and SLM's ``tsurf`` field supplied through the canonical SurfaceModel
radiation input are examples of external providers. The simplified prognostic state does not override either
provider.

The per-column temperature resolver retains its existing fallback order: valid external/LSM absolute
temperature, valid prognostic SEB absolute temperature when offered, valid SurfaceLayer potential temperature
converted to absolute temperature, then the scalar ``erf.rad_t_sfc`` fallback.

.. _sec:TwoStreamLandForcing:

Radiative Forcing of a Land-Surface Model
-------------------------------------------------

With ``erf.land_surface_model = NOAHMP`` the two-stream model supplies the radiation Noah-MP
integrates on, as RRTMGP does. After each column sweep it stores, per column,

- ``SWDOWN``: the total downwelling shortwave at the surface, direct plus diffuse [W/m^2] --
  the incident flux, not the net, since Noah-MP applies its own albedo;
- ``GLW``: the downwelling longwave at the surface [W/m^2];
- ``COSZEN``: the cosine of the solar zenith angle of that sweep, floored at zero.

The fluxes are the surface-interface values of ``rad_fluxes`` (components 1 and 3 at the lowest
interface), after any clear/cloudy blending, so they are the same fluxes that heat the
atmosphere. They are copied into Noah-MP's ``sw_flux_dn``, ``lw_flux_dn`` and
``cos_zenith_angle`` fields every step, since the sweep runs every step; the land model runs
after the dycore and so always sees the current step's radiation. Those fields are part of the
land model's checkpointed data, and a restarted run refills them before its first land step.
Until a sweep has run on a level the stored fields hold the land model's undefined sentinel
rather than zero, so a copy made before one is caught by Noah-MP's missing-input check
instead of being taken as a dark, 0 K sky.
The two-stream model is broadband, so the visible / near-infrared direct / diffuse split that
RRTMGP also provides is not written; Noah-MP does not read it. SLM does, so the two-stream model
does not feed SLM.

In the other direction the sweep reads Noah-MP's surface: its broadband ``albedo``
(reflected over incident shortwave, so the shortwave the sweep reflects at the ground is the
shortwave Noah-MP reflects -- not ``sfc_alb_dir_vis``, the visible direct-beam band of the four
RRTMGP takes, which over vegetation is several times smaller), its emissivity ``sfc_emis`` and
its skin temperature ``t_sfc``. Each is taken column by column where Noah-MP holds a value.
Noah-MP leaves its undefined placeholder over open water and sea ice, in the albedo at night, and
everywhere before its first step, which runs after the first radiation call; those columns take
``erf.radiation.surface_albedo_sw``, ``erf.radiation.surface_emissivity_lw`` and the
surface-layer or ``erf.rad_t_sfc`` temperature. With ``erf.radiation.seb_enable`` the balance's
inputs from Noah-MP (the absorbed shortwave ``sav + sag``, the net longwave ``-fira``, the ground
flux ``grdflx`` and the 2 m humidity) follow the same rule and fall back to the ``seb_*_default``
constants, so a column Noah-MP did not compute no longer carries the placeholder into the
balance. Noah-MP exposes no ``hfx`` or ``lh`` field, but over land the surface layer applies
Noah-MP's fluxes (it takes :math:`u_*` and :math:`\theta_*` from them), so with the default
``seb_turbulent_flux_source = surface_layer`` the balance removes Noah-MP's :math:`H` and
:math:`\text{LE}` over land and the surface layer's own over water.

On a refined run every level hands its own Noah-MP its forcing. A level whose boxes span the
domain in z writes the fluxes of its own sweep. A nested patch, which does not sweep, takes
its parent's, interpolated (piecewise constant) with the rest of its radiation fields. A
finer level runs Noah-MP on that forcing only if it has a land setup file of its own; otherwise
it takes its land state from level 0 (see :doc:`../CouplingToNoahMP`). Without a land model, or with SLM, the two-stream model
stores nothing extra and its results are unchanged. The case
``Exec/RegTests/NoahMP_Ideal/inputs_noahmp_twostream`` exercises the coupling.

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
