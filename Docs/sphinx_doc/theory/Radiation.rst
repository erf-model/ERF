
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Radiation:

Radiation
=========

Radiative transfer in ERF includes both a full k-distribution model (RRTMGP) for detailed studies
and a simplified two-stream model for idealized and intermediate-complexity atmospheric simulations.
This section describes RRTMGP's surface terrain correction and the two-stream radiation model's
physics. The two-stream model employs Beer-Lambert direct-beam
attenuation for shortwave radiation, a Meador-Weaver two-stream diffuse field combined by the adding method
with the surface albedo, gray-gas longwave two-stream, and an optional Simplified Surface Energy Balance (SEB) module.

.. _sec:TopographicRadiation:

Topographic correction of surface radiation
-------------------------------------------

With RRTMGP on a static terrain-following mesh, ``erf.topo_shading`` removes
direct sunlight obstructed by resolved terrain, and ``erf.slope_rad`` adjusts
the direct beam for surface slope. Both options default to false; runtime
controls are listed in :ref:`inputs-topographic-radiation`. The correction
follows the surface-only approach used by WRF: it changes the land model's
shortwave forcing while retaining the column solver's atmospheric fluxes and
heating. It is not a three-dimensional atmospheric radiative-transfer solve.

Let :math:`h(x,y)` be the terrain elevation, :math:`\boldsymbol{s}` the unit
vector toward the sun at the radiation interval's midpoint in grid coordinates,
and :math:`\bar\mu` RRTMGP's interval-averaged zenith cosine. The upward
terrain normal and nonnegative incidence cosine are

.. math::

   \boldsymbol{n} = \frac{(-h_x,-h_y,1)}{\sqrt{1+h_x^2+h_y^2}},
   \qquad I = \max(0,\boldsymbol{n}\cdot\boldsymbol{s}).

For :math:`\bar\mu>10^{-4}`, with both corrections enabled, the surface
shortwave flux is

.. math::

   F_{\mathrm{SW,terrain}} = F_{\mathrm{diffuse}}
       + (1-S)\,\frac{I}{\bar\mu}\,F_{\mathrm{direct,horizontal}},

where :math:`S=1` if terrain obstructs the ray toward the sun or the midpoint
sun has :math:`s_z\leq10^{-4}`, and :math:`S=0` otherwise. ERF sets :math:`I=0`
at that midpoint cutoff. Disabling shading sets :math:`S=0`; disabling the
slope correction replaces :math:`I/\bar\mu` with 1. When the averaged cosine
is at or below :math:`10^{-4}`, the correction is bypassed, with shadow mask 0
and direct factor 1. The midpoint sun can be below its cutoff while the
interval-averaged cosine remains above its cutoff; then the direct factor is 0
when either correction is enabled. A surface facing away from the sun has
:math:`I=0` even when no other terrain casts a shadow on it.

RRTMGP's visible and near-infrared direct fluxes receive the same factor.
Their diffuse components are preserved, and the corrected components supply
the broadband sum for Noah-MP. SLM consumes the spectral components. Longwave
and the radiation solver's solar zenith cosine are preserved. Diffuse sky-view
reduction and radiation reflected or emitted by surrounding terrain are not
included.

Terrain geometry and search
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Surface heights are the horizontal average of the four bottom mesh nodes
around each mass point. This uses the final ERF mesh after WRFInput, Metgrid,
or idealized terrain initialization. Height differences and map factors give
the physical slopes. Latitude/longitude gradients rotate the geographic sun
vector onto the grid axes; idealized grids without those fields use east and
north as their grid axes.

For each target cell, the search samples the terrain along a ray toward the
sun, with bilinear interpolation between mass points. A terrain height above
the ray blocks the direct beam. The maximum horizontal ray distance is
``erf.shadlen``, default 25 km. Periodic directions wrap, and nonperiodic
domain boundaries terminate the search. Physical distances use the target's
map factors throughout the ray, a local tangent-plane approximation.

The cached terrain includes off-rank heights and ancestor terrain outside
refined grids, with piecewise-constant prolongation from ancestors. Regridding
invalidates the affected geometry caches. When a fine level inherits radiation,
ERF interpolates the raw surface forcing before applying the fine level's
terrain correction, avoiding a second correction of already adjusted fluxes.
Fine Noah-MP levels still follow the existing land-model rule: a level without
its own land setup file inherits level-0 land state.

Update timing
^^^^^^^^^^^^^

The mask and slope factor update when RRTMGP computes radiation, every
``erf.rad_freq_in_steps`` ERF steps, and remain fixed between updates.
RRTMGP averages its zenith cosine over the upcoming radiation interval; the
terrain correction uses the sun direction at that interval's midpoint, keeping
the column solver's zenith cosine as the land-model input. This approximates
shadow motion over the interval, so resolving a moving shadow requires a
sufficiently short radiation interval. Raw surface forcing is not checkpointed:
after restart or a level regrid, a full-column level without it forces an
RRTMGP update on its next step, even if that step is off schedule. With
``erf.rad_freq_in_steps > 1``, restarts may therefore lose bitwise equivalence
to uninterrupted runs and regrids may add radiation calls. The forced update
still averages the zenith cosine over the full interval. A fine level that
inherits radiation interpolates the parent's forcing and uses its radiation
time for the terrain correction instead of forcing an update.

See :ref:`sec:Plotfile2DTopographicRadiation` for diagnostics of the cast shadow,
slope factor, and original and corrected shortwave fluxes.


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
with ``sw_cloud_omega`` and ``sw_cloud_g``; ``sw_kext_cloud`` :math:`\approx 3 / (2 \rho_w r_{\text{eff}})`,
i.e. :math:`1500 / r_{\text{eff}}` m²/kg with :math:`r_{\text{eff}}` in µm, gives 150 m²/kg for 10 µm droplets). The prescribed cloud band, the moisture
coefficients and the aerosol term are added on top as before. The mass-path terms are then a property of
the atmosphere, not of the grid; the prescribed cloud band, the moisture coefficients and a constant aerosol
profile are still added per layer, so they keep a dependence on the vertical resolution (see
:ref:`sec:TwoStreamVsRRTMGP`). The longwave band uses the mass path described below.
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

so this part of the column optical depth is independent of the vertical resolution (terms added per
layer, listed in :ref:`sec:TwoStreamVsRRTMGP`, are not), water vapor and cloud water
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

   \frac{dF_{\uparrow}}{d\tau} = F_{\uparrow} - \sigma T^4, \quad
   \frac{dF_{\downarrow}}{d\tau} = -F_{\downarrow} + \sigma T^4

where :math:`\tau` is the longwave flux optical depth, counted downward from the top. Across a layer of
optical depth :math:`\Delta\tau` at temperature :math:`T` this gives
:math:`F_{\text{out}} = F_{\text{in}} e^{-\Delta\tau} + \sigma T^4 (1 - e^{-\Delta\tau})`, so an opaque layer
emits :math:`\sigma T^4`.
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
  already reached rather than from ``erf.rad_t_sfc``. Where the parent column or one of its
  eight neighbors is a water column of the land mask, no slope is used, and a fine land column
  over a coarse water column takes the mean of the neighboring coarse land columns.
- **Fine levels are averaged down.** After the finer levels advance, their surface state is
  averaged onto the coarse level, exactly as the atmospheric state is. Where either level has
  a water column of the land mask, the mean is over the fine land columns, and a coarse water
  column keeps its values (:ref:`Inputs <sec:Inputs>`, the land mask in the balance). This runs under
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
  energy balance, when on, supplies its own state ahead of the surface layer on land columns;
  on a water column of the land mask the balance does not run, and the order above holds. The surface
  layer works in potential temperature, so its value is converted to temperature with the
  Exner function evaluated at the physical surface pressure diagnosed from the lowest atmospheric cell
  before it enters the :math:`\sigma T_s^4` emission; with a
  surface layer present, ``erf.rad_t_sfc`` is the initial value of the prognostic surface
  temperature when the surface energy balance evolves one, and unused otherwise.
- **First skin temperature from the initial file.** One value for every column is a poor first
  guess over real terrain: the skin of a mountain and of a valley can differ by tens of kelvin,
  and the restoring term alone brings the skin back over the time :math:`\tau/2\pi`, 3.8 hours
  for a restoring period :math:`\tau` of one day. With
  ``erf.radiation.seb_init_skin_from_file = true`` and the prognostic balance, a land column
  :math:`(i,j)` of a level that read an initial file begins at that file's skin temperature,

  .. math::

     T_s(i,j,t_0) =
     \begin{cases}
       T_{\text{file}}(i,j), & \text{land (mask} \neq 0\text{)} \\
       T_{\text{rad}},      & \text{water (mask} = 0\text{)}
     \end{cases}

  where :math:`T_{\text{file}}` is ``TSK`` of a wrfinput file or ``SKINTEMP`` of a metgrid file,
  an absolute temperature in kelvin, and :math:`T_{\text{rad}}` is ``erf.rad_t_sfc`` at level 0
  and, at a finer level, the value the fill from its parent gave the column. The balance
  does not run on a water column, so its value there is not used, except by a fine land column
  with no coarse land column next to it, which takes the held coarse value as described above.
  A level without a file of its own is filled from its parent as described above, and so starts
  from the parent's file values. So does a finer level that has a file when it takes its
  radiation from its parent, or when ``erf.nc_low_file`` is set without ``erf.nc_low_file_<lev>``:
  the level's surface temperature then follows its parent's, and a line in the log says so.
  Only the first value changes: the balance, its restoring temperature ``seb_t_deep_default``
  and the surface layer are as before. Unless ``seb_init_t_deep_from_file`` is set, the
  restoring temperature :math:`T_{\text{deep}}` stays one
  value for every column, so a column whose file value differs from it by :math:`\Delta T` feels
  the restoring flux :math:`C_s (2\pi/\tau)\,\Delta T` from the first step, with :math:`C_s` the
  skin heat capacity per area: 1.45 W m\ :sup:`-2` per kelvin for
  :math:`C_s = 2\times10^{4}` J m\ :sup:`-2` K\ :sup:`-1` and :math:`\tau` of one day. A
  column that starts at 270 K with the default :math:`T_{\text{deep}}` of 300 K is pulled back
  by about 43 W m\ :sup:`-2`, over the time :math:`\tau/2\pi` of 3.8 hours, so a start far
  from ``seb_t_deep_default`` relaxes toward it and the per-column start is partly lost. Set
  ``seb_t_deep_default`` near the mean of the file's land values; the log line of each level
  prints the largest difference, or set the deep temperature per column from the file as
  described next. A real file can hold land values outside
  the default range of 200 K to 340 K (a hot desert, a polar plateau); the two range inputs
  widen it. Without
  ``erf.radiation.seb_surface_layer_uses_skin`` the surface layer keeps using the file's skin
  temperature on land for the whole run, while the balance's skin evolves from it. The run stops
  when the file of level 0, or of a finer level that takes its own file's value, lacks the field,
  when the radiation model is not the two-stream one, when a land-surface, urban or surface
  model is on, and when a land
  value is not finite or lies outside ``[seb_prognostic_t_min_k, seb_prognostic_t_max_k]``, the
  range the balance holds :math:`T_s` in. A level made later in the run is filled from its
  parent, whose skin has evolved by then, also when it has a file. A restart takes :math:`T_s`
  from the checkpoint, or ``erf.rad_t_sfc`` when the checkpoint was written without the balance.
- **Deep temperature from the initial file.** The restoring term of the balance,
  :math:`-C_s (2\pi/\tau)(T_s - T_{\text{deep}})`, pulls the skin of every column toward the
  deep temperature :math:`T_{\text{deep}}`, which is the one value ``seb_t_deep_default``
  unless this input is set. Over real terrain one value cannot be near the skin of both a
  mountain and a valley. A wrfinput file carries ``TMN``, the soil temperature at the lower
  boundary of WRF's soil column [K], in every column. With ``erf.radiation.seb_init_t_deep_from_file = true`` and the prognostic
  balance,

  .. math::

     T_{\text{deep}}(i,j,t_0) =
     \begin{cases}
       T_{\text{TMN}}(i,j), & \text{land (mask} \neq 0\text{)} \\
       T_{\text{deep,0}},   & \text{water (mask} = 0\text{)}
     \end{cases}

  where :math:`T_{\text{TMN}}` is the file's ``TMN`` and :math:`T_{\text{deep,0}}` is
  ``seb_t_deep_default`` at level 0 and, at a finer level, the value the fill from its parent
  gave the column. The balance does not run on a water column, so the value there is not used,
  except as described for the skin above.

  **What changes in the balance.** Without the input ERF sets :math:`T_{\text{deep}}` to
  ``seb_t_deep_default`` before every step. With it :math:`T_{\text{deep}}` is part of the
  surface state, like :math:`T_s`: it is not reset, a new level takes it from its parent,
  fine levels are averaged down onto the coarse level under two-way coupling, a regrid keeps
  it, and the checkpoint holds it (see `Multiple Levels`_ for the three rules). Without
  ``seb_deep_prognostic_enable`` it then keeps the file's value for the whole run, apart from
  the averaging down under a finer level. With ``seb_deep_prognostic_enable`` the file's value
  is the first value of the evolving :math:`T_{\text{deep}}` (`Force-Restore Hydrology`_), in
  place of ``seb_t_deep_default``.

  **Levels.** A finer level that read a wrfinput file of its own takes that file's ``TMN`` on
  land. A level without a file is filled from its parent, and so holds the parent's ``TMN``,
  interpolated. So is a level that takes its radiation from its parent, and a level made
  later in the run, also when it has a file. Unlike the skin temperature, ``TMN`` does not
  depend on ``erf.nc_low_file``: ERF reads it from the initial file only.

  **When the run stops.** The reader keeps ``TMN`` only when every column of the level has a
  value between 0 K and 1000 K, water columns included, and reports a field it does not keep.
  The run stops when the file of level 0, or of a finer level that has a file, gave no ``TMN``
  the reader kept (so also with ``erf.init_type`` other than ``WRFInput``), when the radiation
  model is not the two-stream one, when the balance is not prognostic, when a land-surface,
  urban or surface model is on, and when a land value is not finite or lies outside
  ``[seb_prognostic_t_min_k, seb_prognostic_t_max_k]``. The last range is the one the balance
  holds :math:`T_s` in: a skin cannot reach a deep temperature outside it.

  **What the file's value stands for.** With :math:`\tau` of one day the term restores the
  skin toward :math:`T_{\text{deep}}` over a few hours, so a :math:`T_{\text{deep}}` that
  differs from the mean of the skin over the day by :math:`\Delta T` adds the steady flux
  :math:`C_s (2\pi/\tau)\,\Delta T` to the balance, 1.45 W m\ :sup:`-2` per kelvin with the
  defaults. ERF takes ``TMN`` as the file gives it and does not test how near it is to that
  mean. Without ``seb_deep_prognostic_enable`` the difference acts for the whole run. With
  it :math:`T_{\text{deep}}` relaxes toward the skin with the e-folding time :math:`\tau`,
  so the file's value is mostly replaced after a few days.

  **Restart.** The file is not read on a restart. :math:`T_{\text{deep}}` comes from
  ``TwoStream_TDeep`` of the checkpoint, so a restarted run continues bit for bit. A checkpoint
  written without the input has no such file: every column then holds ``seb_t_deep_default``,
  which is what that run used, and a NOTE says so. A fine land column with no coarse land column next to it also holds
  ``seb_t_deep_default``, as described for the skin. The input must stay set on the restart;
  without it ERF again sets :math:`T_{\text{deep}}` to ``seb_t_deep_default`` before every
  step, and a NOTE says that the checkpoint's field is not read.
- **Call cadence.** The sweep runs once every slow step, from the old state, and there is no
  call-frequency input: the two-stream model costs about 1.2 ms per step for seven thousand
  cells (see :ref:`sec:TwoStreamVsRRTMGP`), so unlike RRTMGP (``erf.rad_freq_in_steps``) it is not worth skipping steps.
- **Diagnostics file.** The diagnostics are off by default. Setting
  ``erf.radiation.diag_enable = true`` writes ``radiation_diag.dat``
  (``erf.radiation.diag_file``) in the run directory, with a ``pre_dycore`` and a
  ``post_dycore`` row per step for every level that sweeps. Each of them appends to the one
  file and the last column, ``level``, tells the rows apart. A level interpolated from its
  parent contributes no rows at all -- it runs no sweep, so it has no fluxes of its own to
  report -- so on a refined run the rows present are those of the sweeping levels, not one
  set per level in the hierarchy. The file is appended to rather than truncated, as ERF's other
  data logs are, so a rerun in the same directory extends the previous run's rows.

.. _sec:TwoStreamVsRRTMGP:

Comparison with RRTMGP
--------------------------------------

The two-stream model is a gray, broadband scheme. RRTMGP is a k-distribution model with gas
optics from lookup tables. The table below lists what each one represents in ERF, and what the
two-stream model leaves out. It describes the models; it is not a measurement of accuracy.
No case in ERF yet compares the fluxes or heating rates of the two models. Two-stream keys are
read under ``erf.radiation.``; RRTMGP keys are written in full.

.. list-table::
   :header-rows: 1
   :widths: 16 28 28 28

   * - Aspect
     - Two-stream
     - RRTMGP
     - Limitation of the two-stream model
   * - Spectral resolution
     - One gray shortwave band and one gray longwave band
     - 224 shortwave and 256 longwave g-points
     - No water-vapor window and no CO\ :sub:`2` band. Absorption does not depend on pressure
       or temperature, and longwave scattering is not included. Longwave cooling-to-space and
       heating-rate profiles are therefore approximate.
   * - Gases
     - Dry air and water vapor only, through user coefficients
       (``tau_sw_coeff_qv`` and ``tau_lw_coeff_qv`` with ``per_layer``; ``sw_kabs_vapor`` and
       ``lw_kabs_vapor`` with ``mass``; ``lw_kabs_vapor`` also with ``per_layer`` when
       ``lw_mass_absorption_enable`` is true)
     - H\ :sub:`2`\ O, CO\ :sub:`2`, O\ :sub:`3`, N\ :sub:`2`\ O, CO, CH\ :sub:`4`,
       O\ :sub:`2`, N\ :sub:`2`
     - No CO\ :sub:`2`, ozone or methane. With ``per_layer`` water vapor has no radiative
       effect unless ``tau_sw_dynamic_enable`` / ``tau_lw_dynamic_enable`` are true and the
       coefficients are set, or (longwave only) ``lw_mass_absorption_enable`` is true. All are
       off or zero by default.
   * - Optical depth
     - ``tau_model = per_layer`` (default) or ``mass``
     - Computed from the gas state
     - With ``per_layer`` the column optical depth grows with the number of vertical cells
       (42 cells with the default longwave value give a column optical depth of 42), so its
       coefficients must be tuned again when the grid changes. ``mass`` removes this dependence
       for the gas and ``qc`` terms (``lw_mass_absorption_enable`` does the same for the
       longwave clear-sky term with ``per_layer``). Three terms are still added per layer in
       both modes: the ``tau_*_coeff_*`` terms, the cloud band ``cloud_tau_per_layer``
       (``tau_profile_type = cloud_layer``), and the constant aerosol profile
       ``aerosol_tau_per_layer``. With any of these, the column optical depth still changes
       when the grid is refined. The mass
       coefficients are flux coefficients (no diffusivity factor), to be tuned rather than
       taken from line data.
   * - Cloud phase
     - Liquid water (``qc``) only, and only through an option that reads it: ``tau_model =
       mass``, ``lw_mass_absorption_enable`` (longwave), the ``tau_*_coeff_qc`` terms with
       dynamic optical depth, or ``cloud_fraction_prog_enable`` (with ``tau_profile_type = cloud_layer`` and
       ``cloud_fraction > 0``). With ``per_layer`` and none of
       these on (the default), ``qc`` has no radiative effect, and clouds are only the
       prescribed height band (``tau_profile_type = cloud_layer``).
     - Liquid, and ice when the moisture model carries cloud ice
     - Ice and mixed-phase clouds are not represented. By default, the model's clouds do not
       follow the resolved cloud water.
   * - Cloud fraction and overlap
     - Clear and cloudy fluxes blended linearly by one domain-wide value, ``cloud_fraction``
       (default 0, so only the clear column is computed; with ``mass`` that column already
       includes the resolved ``qc``, and ``cloud_fraction`` blends in only the prescribed band). The fraction diagnosed from RH and
       ``qc`` (``cloud_fraction_prog_enable``; see `Cloud Fraction Diagnosis`_) does not
       change this blend; it only scales the cloud-band optical depth in the cloudy column.
     - Cloud fraction 1 where :math:`q_c + q_i > 0`, else 0. McICA with maximum-random overlap
       is coded, but with this fraction ``erf.rad_do_subcol_sampling`` has no effect (#5212).
     - The overlap of cloud layers is not represented, and the blend weight does not vary in
       space. With ``cloud_fraction = 0``, setting ``cloud_fraction_prog_enable = true`` gives
       clear-sky fluxes only. Neither model represents overlap in ERF today; only the
       two-stream model represents partial cloud cover.
   * - Cloud optics
     - Cloud band, both modes and both bands: ``cloud_tau_per_layer``, with
       ``cloud_single_scattering_albedo`` and ``cloud_asymmetry_factor`` in the shortwave.
       ``mass`` only: ``qc`` with ``sw_kext_cloud``, ``sw_cloud_omega`` and ``sw_cloud_g`` in the
       shortwave, and ``lw_kabs_cloud`` in the longwave, which gives the Stephens (1978)
       emissivity :math:`1 - e^{-0.158\,\text{LWP}}` (also with ``lw_mass_absorption_enable``)
     - Cloud-optics lookup tables
     - The result depends on the chosen cloud parameters and mode. The ``sw_*`` cloud keys have
       no effect with ``per_layer``; the ``lw_kabs_*`` keys act with ``mass`` or with
       ``lw_mass_absorption_enable``.
   * - Aerosol
     - Prescribed aerosol optical depth (``aerosol_enable``, default false)
     - Not supported (``erf.rad_do_aerosol = true`` stops the run)
     - --
   * - Surface coupling
     - Direct and diffuse broadband albedo, and emissivity
     - Direct and diffuse albedo in the visible and near-infrared bands, and surface emissivity
       (from the land-surface model or surface layer when present, else 0.98)
     - No visible / near-infrared split. The two-stream model forces Noah-MP with broadband
       fluxes but cannot force SLM (see `Radiative Forcing of a Land-Surface Model`_).
   * - Top of the domain
     - Taken as the top of the atmosphere: no downward longwave, full solar beam
     - Same: no layer is added above the domain top
     - Shared by both models. Nothing is emitted, absorbed or scattered above the top. With
       ``mass`` and the default coefficients, a gray-column estimate gives a downward longwave
       flux at the surface about 11% too small for a 2 km top (28% for 1 km), and a surface
       direct beam about 9% too large for a 2 km top with the sun overhead. With
       ``per_layer`` the column is usually opaque, so the longwave effect is small.
   * - Solar geometry
     - Instantaneous :math:`\cos\theta_z` at each step
     - :math:`\cos\theta_z` averaged over ``erf.rad_freq_in_steps`` steps (one step by
       default)
     - Small difference that grows with ``erf.rad_freq_in_steps``; both use the same orbital
       code.
   * - Cost
     - One column sweep per step (the ``post_dycore`` call reuses the ``pre_dycore`` sweep);
       about 1.2 ms per step for 7,056 cells
     - About 350 to 360 times the two-stream cost per step on one CPU core; called every
       ``erf.rad_freq_in_steps`` steps
     - --

The cost row comes from ``Exec/CanonicalTests/Radiation/TwoStream_vs_RRTMGP_Timing``
(one CPU core, the same grid, sounding and moisture model for both; measured in November 2026;
cost scales linearly with cell count from 7,056 to 56,448 cells for both). RRTMGP is designed for GPUs, so the ratio on a GPU
is expected to differ. The top-of-domain numbers come from a standalone gray-column calculation
with a standard atmosphere, not from an ERF run. The two-stream checks in
``Exec/CanonicalTests/Radiation/TwoStream_Benchmark_Suite`` compare against closed-form
solutions (clear-sky shortwave, isothermal longwave) and qualitative checks (cloud layers,
time stepping), not against RRTMGP.

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
:math:`10^6` s/m after the vapor-deficit factor. Noah-MP applies the same factors per sunlit
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
when :math:`\text{RH}_g q_\text{sat} < q_\text{air}`, takes up vapor (LE < 0), as in Noah-MP.

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

.. _sec:TwoStreamForceRestoreHydrology:

Force-Restore Hydrology
-----------------------

The balance above has no rain, no runoff, and fixed deep values :math:`q_{\text{deep}}` and
:math:`T_{\text{deep}}`. Three opt-in inputs add the remaining terms of Noilhan and Planton (1989,
NP89) and Mahfouf and Noilhan (1996, MN96). With all three off the results are unchanged
(``ERF_ForceRestoreHydrology.H``). The step runs after the balance's own update, on every
column outside a building footprint.

**Rain** (``seb_precipitation_source = microphysics``). The rain rate :math:`P`
[kg m\ :sup:`-2` s\ :sup:`-1`] is the change of the microphysics scheme's surface
precipitation accumulation over the step, divided by the step. It enters the surface layer as in
NP89 Eq. (11), with :math:`d_1 = d_s`:

.. math::

   q_s \leftarrow q_s + \frac{P \, \Delta t}{\rho_w d_s}.

NP89 Eq. (11) multiplies :math:`P - E` by one coefficient :math:`C_1`. The balance removes
:math:`E` with :math:`C_1 = 1`, so rain uses :math:`C_1 = 1` too. NP89 (Eq. 19, Table 2) gives
:math:`C_1` = 0.08 to 0.34 at saturation, rising as the soil dries; Deardorff (1978, Eq. 13) gives
0.5 to 14. With :math:`C_1 = 1` the surface layer fills faster near saturation than with NP89's
value, so heavy rain runs off sooner.

Water above the porosity :math:`w_{\text{sat}}` (``seb_soil_porosity``) runs off (NP89 Sec. 4a).
The surface layer's runoff is at most the step's rain: water above :math:`w_{\text{sat}}` that
:math:`q_s` held before the rain (a restore overshoot, or a restart from a run without the
hydrology) is clamped off :math:`q_s` and stays in the column. The surface runoff is taken
before the drainage, so water that never entered the soil does not drain. With rain,
``seb_prognostic_q_max`` must be at or above :math:`w_{\text{sat}}`; a lower cap would remove
rain that the runoff does not count. When the balance skips a step (the first step of a level
built from a coarser one), that step's rain is carried into the next: the next step's rain rate
is the mean over both steps, and the step takes the rain of both. A rain rate is used once: a
step that runs without a new rate (a level whose radiation is interpolated from the coarser
level) takes no rain, and ``seb_rain_rate`` there keeps the last rate for output. A level whose
boxes do not span the domain in :math:`z` never runs the balance and takes no rain.

**Deep reservoir** (``seb_deep_prognostic_enable``). :math:`q_{\text{deep}}` becomes
:math:`w_2`, the mean water content of a soil column :math:`d_2` deep
(``seb_deep_layer_depth_m``) that contains the surface layer (NP89 Eq. 12, MN96 Eq. 4):

.. math::

   \frac{dq_{\text{deep}}}{dt} = \frac{P - E}{\rho_w d_2}
   - \frac{C_3}{\tau_{1d}} \max\left(0, q_{\text{deep}} - w_{fc}\right),
   \qquad E = \frac{\text{LE}}{L_v}.

Since :math:`w_2` is the column mean, the restore term of :math:`q_s` moves water within the
column and does not appear here. The deep temperature relaxes toward the surface (NP89 Eq. 6):

.. math::

   \frac{dT_{\text{deep}}}{dt} = \frac{T_s - T_{\text{deep}}}{\tau},

with :math:`\tau` = ``seb_restore_timescale_s``, the period NP89 Eqs. (5) and (6) share.

**Drainage** (``seb_drainage_enable``). The second term above removes water above field capacity
:math:`w_{fc}` (``seb_soil_moisture_fc``) through the bottom of the column. MN96 Eq. (16) scales
:math:`C_3` as :math:`1/d_2`; their Table 1 gives :math:`C_3` at :math:`d_2 = 1` m, so the input
is the product ``seb_drainage_c3_m`` :math:`= C_3 d_2` [m]. MN96 Eq. (15) makes :math:`C_3`
proportional to the restore period, and Table 1 uses one day, so the drainage uses
:math:`\tau_{1d} = 86400` s whatever ``seb_restore_timescale_s`` is.

**Why exact relaxation:** the drainage and deep-temperature terms are linear in the state, so each
is integrated exactly, :math:`x^{n+1} - x_{eq} = (x^n - x_{eq}) e^{-k \Delta t}`. The update cannot
pass field capacity or :math:`T_s` for any step. Rain and evaporation are forward Euler, as in the
balance.

**What is conserved:** with the deep reservoir, each land column satisfies, to round-off,

.. math::

   \rho_w d_2 \left(q_{\text{deep}}^{n+1} - q_{\text{deep}}^n\right)
   = (P - E)\,\Delta t - D - R + S,

where :math:`D` is the drainage, :math:`R` the runoff, and :math:`S` the deficit: the water an
evaporation took that the column did not hold, recorded when :math:`q_{\text{deep}}` would go
below zero. The latent heat flux is not limited by :math:`q_{\text{deep}}`, so a non-zero
deficit is water the coupled system did not conserve: the atmosphere received vapor the soil did
not hold. The running totals of rain, runoff, drainage and deficit [kg m\ :sup:`-2`] are
checkpointed and written to the 2D plotfile (``seb_rain_total`` and the others), with
``seb_rain_rate``, ``seb_q_deep`` and ``seb_t_deep``. With ``erf.radiation.v >= 2`` each step
logs the domain sums and the largest column residual.

The dust model reads this soil through a read-only accessor (``ForceRestoreSoil``: :math:`q_s`,
:math:`q_{\text{deep}}`, :math:`T_s`, the rain rate, the depths, porosity, field capacity and
wilting point), for level 0 only. The dust emission does not use it yet.

References: Deardorff, J. W., 1978, J. Geophys. Res., 83, 1889-1903; Noilhan, J. and S. Planton,
1989, Mon. Wea. Rev., 117, 536-549; Mahfouf, J.-F. and J. Noilhan, 1996, J. Appl. Meteor., 35,
987-992.

.. _sec:TwoStreamForceRestoreSnow:

Force-Restore Snow
------------------

``erf.radiation.seb_snow_enable`` puts one bulk snow layer on the force-restore surface
(``ERF_ForceRestoreSnow.H``). It needs ``seb_prognostic_enable``, ``seb_use_radiation_fluxes``,
``seb_precipitation_source = microphysics``, no ``surface_albedo_sw_diffuse``, no
``seb_vegetation_type`` and ``seb_prognostic_t_min_k`` at or below 273.16; the run stops at
start-up otherwise. It is the snow routine of the ISBA scheme as Douville, Royer and
Mahfouf (1995, DRM95) describe it in their Sec. 2.1, for bare ground (vegetation fraction 0).
The default is off, and the results are then unchanged. The state is the snow water equivalent
:math:`W_n` [kg m\ :sup:`-2`]; 1 kg m\ :sup:`-2` of water is 1 mm.

**Mass** (DRM95 Eq. 6):

.. math::

   \frac{dW_n}{dt} = P_n - E_n - F_n,

where :math:`P_n` is the frozen precipitation, :math:`E_n` the sublimation and :math:`F_n` the
melt, all in kg m\ :sup:`-2` s\ :sup:`-1`. :math:`P_n` is the change over the step of the
microphysics scheme's snow, graupel and hail accumulations at the surface. With the snow layer
on, the hydrology's rain :math:`P` is the rest: the scheme's total minus the frozen part (the
rain accumulation for a scheme without a total), so the two always add up to what the
hydrology takes without the snow layer. A scheme without those accumulations (Kessler,
SAM without ice) has no frozen part: snow then comes only from
``seb_snow_initial_swe_kg_m2``, and the run prints a note saying so.

**Cover fraction** (DRM95 Eq. 7):

.. math::

   p_{ng} = \frac{W_n}{W_n + W_{crn}},

with :math:`W_{crn}` = ``seb_snow_critical_swe_kg_m2`` (10 kg m\ :sup:`-2` in DRM95).

**Albedo** (DRM95 Eq. 8 with no vegetation):

.. math::

   \alpha_t = (1 - p_{ng})\,\alpha_g + p_{ng}\,\alpha_n,

where :math:`\alpha_g` is the snow-free albedo (``surface_albedo_sw``) and :math:`\alpha_n` =
``seb_snow_albedo`` (0.70 in DRM95, constant). The column sweep uses :math:`\alpha_t` for the
direct and the diffuse beam, so ``surface_albedo_sw_diffuse`` cannot be given with the snow
layer. The albedo acts on the balance through the sweep's surface fluxes, so the snow layer
needs ``seb_use_radiation_fluxes``.

**Sublimation.** DRM95 states that :math:`p_{ng}` divides the surface evaporation between the
soil and the snow and gives no equation. Here :math:`E_n = p_{ng} E` for an upward
:math:`E = \text{LE}/L_v`, at most the snow present. A downward flux (dew, frost) stays with
the soil.

**Melt** (DRM95 Eq. 12). The balance above gives a first estimate :math:`T_s^*`. Where snow
lies and :math:`T_s^* > T_{00} = 273.16` K,

.. math::

   F_n\,\Delta t = \min\!\left(\frac{C_s\,(T_s^* - T_{00})}{L_f},\; W_n'\right),
   \qquad
   T_s = T_s^* - \frac{L_f\,F_n\,\Delta t}{C_s},

where :math:`C_s` = ``seb_surface_heat_capacity`` [J m\ :sup:`-2` K\ :sup:`-1`] is the inverse
of the coefficient :math:`C_t` of DRM95, :math:`L_f = 3.337 \times 10^5` J kg\ :sup:`-1` is
the latent heat of fusion and :math:`W_n'` is the snow after this step's snowfall and
sublimation.

**What this does:** the surface stays at :math:`T_{00}` while snow remains, and rises above it
only after the last snow has melted (DRM95 Sec. 2.1).

The balance bounds :math:`T_s` at ``seb_prognostic_t_max_k``. On a column with snow that bound
is applied after the melt, not to :math:`T_s^*`: the heat above the bound melts snow
first, and the bound then acts on the temperature the melt leaves. The lower bound
``seb_prognostic_t_min_k`` is applied to :math:`T_s^*` as without snow.

**Order of the step.** The snow step runs after the balance and before the hydrology step, on
every land column outside a building footprint (a water column of the land mask carries no
snow): snowfall, then sublimation with :math:`p_{ng}` of
the step's start, then melt. The albedo the radiation uses in a step is that of :math:`W_n`
at the end of the previous step; a building footprint keeps the snow-free albedo. The rain
is the change over the step of the total accumulation minus the frozen one, not below zero;
the bound is on the change, not on the accumulation. The balance has already taken the whole evaporation from the
soil, so the snow step returns :math:`(E_n + F_n)\,\Delta t` to the hydrology step, which adds
it to the surface layer and the column with the rain. The soil's net loss is then
:math:`E - E_n`, and the melt water wets the soil or runs off. One exception: where the
balance held :math:`q_s` at its lower bound, the soil did not supply the evaporation, and the
returned :math:`E_n\,\Delta t` then adds to the surface layer; the column budget still closes
with ``seb_deep_prognostic_enable``.

**What is conserved:** each land column satisfies, to round-off,

.. math::

   W_n^{n+1} - W_n^n = (P_n - E_n - F_n)\,\Delta t .

The snow water equivalent and the running totals of snowfall, sublimation and melt
[kg m\ :sup:`-2`] are checkpointed (a restart takes them from the checkpoint, not from
``seb_snow_initial_swe_kg_m2``; a checkpoint without the snow layer starts the snow from that
input and prints a note, and a checkpoint with only one of the two snow files stops the run)
and written to the 2D plotfile (``seb_snow_swe``,
``seb_snowfall_total``, ``seb_snow_sublimation_total``, ``seb_snow_melt_total``). They are
filled from the coarse level on a new level and averaged down. With ``erf.radiation.v >= 2``
each step logs the column means and the largest column residual. ``seb_rain_total`` counts the
rain only; the water the snow returned is ``seb_snow_melt_total`` plus
``seb_snow_sublimation_total``, and the column water budget above holds with that water added
to the rain.

**Limits.** This is the simple ISBA routine, not the new scheme of DRM95 Sec. 2.2:

- the snow has no density, depth, heat capacity or thermal conductivity of its own, and its
  albedo does not age;
- the latent heat flux uses :math:`L_v` over snow, so the heat :math:`L_f E_n` is not taken
  from the surface (1.7 W m\ :sup:`-2` for :math:`E` = 10\ :sup:`-5` kg m\ :sup:`-2`
  s\ :sup:`-1` and :math:`p_{ng}` = 0.5; 30 W m\ :sup:`-2` for 10\ :sup:`-4` and 0.9);
- :math:`E` is the balance's evaporation, so it carries the soil's moisture availability:
  snow over a soil at or below wilting does not sublimate;
- a downward :math:`E` (dew, frost) goes to the soil, also under full cover;
- the surface roughness and the longwave emissivity are unchanged by snow (DRM95 Eqs. 9-10
  are not used);
- the soil does not freeze, so melt water that enters a soil below freezing stays liquid;
- no vegetation: ``seb_vegetation_type`` stops the run at start-up with the snow layer;
- a land-surface model (SLM, Noah-MP) owns the surface and stops the run with the snow layer;
- no snow lies on a building footprint: the albedo there is the snow-free one, and the first
  step takes the initial cover off ``seb_snow_swe``.
  The building faces (IBSEB) see the snow only through the ground albedo they are handed;
- snowfall of a step in which the balance did not run (the first step of a new level) is
  carried to the next step, and is lost if a regrid or a restart comes between the two, as
  the hydrology's rain is;
- in single precision, a change of :math:`W_n` in one step below the resolution of a deep
  pack (about 3 x 10\ :sup:`-5` kg m\ :sup:`-2` at 500 kg m\ :sup:`-2`) is not taken. Melt
  heat then stays in :math:`T_s` until it is large enough; snowfall that small is dropped.

The lower bound ``seb_prognostic_t_min_k`` of the balance must be at or below 273.16 K with the
snow layer, or the run stops at start-up: a higher bound would hold a snow-covered surface
above the melting point.

Reference: Douville, H., J.-F. Royer and J.-F. Mahfouf, 1995, Climate Dynamics, 12, 21-35.

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

.. _sec:ForceRestoreVsSLMNoahMP:

Comparison with SLM and Noah-MP
-------------------------------------------------

ERF has three ways to evolve the land surface: the force-restore balance described above, the
Simplified Land Model (:doc:`../SLM`) and Noah-MP (:doc:`../CouplingToNoahMP`). The table lists
what each represents in ERF today. The last column states what the simpler models leave out.
Force-restore keys carry the ``erf.radiation.`` prefix and SLM keys the ``slm.`` prefix; other
keys are written in full.

.. list-table::
   :header-rows: 1
   :widths: 13 22 22 22 21

   * - Aspect
     - Force-restore
     - SLM
     - Noah-MP
     - Limitation
   * - Selection and build
     - ``seb_prognostic_enable`` (it also turns on ``seb_enable`` and ``seb_diagnostic_enable``),
       with ``erf.radiation_model = TwoStream``. ``seb_enable`` or ``seb_diagnostic_enable``
       alone only reports the balance; the skin does not evolve. Part of the standard build.
     - ``erf.land_surface_model = SLM``. Part of the standard build; needs a moisture model.
     - ``erf.land_surface_model = NOAHMP``. Needs ``-DERF_ENABLE_NOAHMP=ON`` (Fortran, NetCDF,
       parallel NetCDF) and, in the run directory, ``namelist.erf``, ``NoahmpTable.TBL`` and a
       land setup file.
     - The ``seb_*`` keys are read only with the two-stream radiation model. With any other
       model they are ignored without a message, except ``seb_init_*_from_file``, which abort.
       Force-restore is an alternative to a land model, not an addition: with SLM or Noah-MP the
       balance's skin does not evolve, and its hydrology and file initialization abort.
   * - Soil column
     - One surface layer (``seb_moisture_layer_depth_m``, 0.1 m), restored toward a deep
       temperature. An optional deep water reservoir (``seb_deep_prognostic_enable``,
       ``seb_deep_layer_depth_m``, 1.0 m).
     - ``slm.nsoil`` layers (default 7) with temperature and soil wetness in each.
       ``slm.soil_dz`` is required for an idealized run; real cases take the file's layers.
     - ``erf.lsm_nsoil`` layers (default 4, must equal ``NSOIL`` in ``namelist.erf``) with
       temperature, total water and liquid water in each.
     - Force-restore has no soil temperature profile. The deep temperature is
       ``seb_t_deep_default`` (300 K), reset every step, unless the deep reservoir is on or
       ``seb_init_t_deep_from_file`` reads ``TMN``.
   * - Surface temperature
     - Force-restore equation, explicit, one step per ERF step, with one bulk heat capacity
       (``seb_surface_heat_capacity``, 2.0e4 J/m²/K, or from ``seb_soil_type``). With
       ``seb_prognostic_enable``, clamped after each step to [``seb_prognostic_t_min_k``,
       ``seb_prognostic_t_max_k``] = [200, 340] K without a message; the upper bound is not applied
       under snow.
     - Ground-skin energy balance solved by Newton iteration (at most 5 iterations, each change
       capped at 5 K, tolerance 0.01 K), with a separate canopy temperature. Running out of
       iterations does not stop the run; a non-finite residual does.
     - Noah-MP's energy balance over vegetated and bare parts of the cell, every
       ``NOAH_TIMESTEP`` (required in ``namelist.erf``).
     - No check that the explicit force-restore step is stable. With skin coupling, a wet surface
       and the default heat capacity, the step must stay below about 80 s; usual ERF steps are far
       below that.
   * - Sensible and latent heat
     - From the surface layer (``seb_turbulent_flux_source = surface_layer``, default) or from
       ``seb_hfx_default`` / ``seb_lh_default`` (``defaults``). The skin temperature reaches the
       surface layer only with ``seb_surface_layer_uses_skin``; the skin moisture only with
       ``seb_surface_layer_uses_moisture``.
     - SLM's own bulk Monin-Obukhov fluxes. On land the surface layer rebuilds
       :math:`u_*`, :math:`\theta_*` and :math:`q_*` from them. SLM's surface temperature always
       reaches the surface layer.
     - Noah-MP's fluxes, used by the surface layer in the same way. Noah-MP's skin temperature
       sets the surface-layer temperature only with
       ``erf.most.use_lsm_skin_temperature = true`` (default false).
     - With ``seb_surface_layer_uses_skin = false`` (default) the balance does not change the
       atmosphere; it only reports a skin temperature.
   * - Radiation
     - Broadband. The two-stream fluxes enter only with ``seb_use_radiation_fluxes = true``;
       otherwise ``seb_sw_flux_default`` and ``seb_lw_flux_default`` (both 0). Broadband
       albedo and emissivity.
     - Visible and near-infrared, direct and diffuse, from RRTMGP only, or from
       ``slm.rad_input_file``. The two-stream model does not feed SLM.
     - ``SWDOWN``, ``GLW`` and ``COSZEN`` from RRTMGP or the two-stream model.
     - With the defaults the force-restore balance gets no radiation. SLM run without RRTMGP or
       a radiation file gets zero radiation and no warning.
   * - Soil water and evaporation
     - One surface water content with a soil-water factor between wilting point and field
       capacity. With ``seb_vegetation_type`` and ``seb_soil_type``, Jarvis canopy and soil
       resistances from Noah-MP's tables.
     - Water in each layer, moved by an implicit (Thomas algorithm) diffusion and gravity solve.
       The bottom is closed except for water above saturation.
     - Noah-MP's soil water, runoff and groundwater options, set in ``namelist.erf``.
     - Force-restore has one water bucket and the Jarvis form only. SLM has no gravity drainage
       out of the bottom layer below saturation.
   * - Precipitation and runoff
     - ``seb_precipitation_source = microphysics`` (default ``none``). Runoff above saturation;
       optional drainage (``seb_drainage_enable``).
     - Microphysics surface precipitation, canopy interception, infiltration, and pond runoff
       above ``slm.mws_mx0``.
     - Microphysics surface precipitation. Runoff is computed but not routed in ERF.
     - Without microphysics none of the three gets precipitation.
   * - Vegetation
     - Optional ``seb_vegetation_type`` (MODIS 1-20 except 15-17), ``seb_vegetation_fraction``
       and ``seb_leaf_area_index``. No canopy temperature or interception.
     - ``slm.landtype0`` (MODIS 0-16) and ``slm.LAI0``, or the file's land use and ``LAI``.
       Jarvis stomata, canopy temperature, interception.
     - Canopy from ``NoahmpTable.TBL``. Stomata Ball-Berry (option 1, default) or Jarvis
       (option 2); optional dynamic vegetation.
     - SLM's defaults (``slm.landtype0 = 16``, ``slm.LAI0 = 0``) give bare soil.
   * - Snow and frozen soil
     - Optional single snow bucket on bare ground (``seb_snow_enable``, default false). No
       frozen soil.
     - No snow. No latent heat of soil freezing; no water moves between a frozen layer and its
       neighbor. Stops if the canopy of a vegetated cell is at or below 273.15 K.
     - Several snow layers (``NSNOW``); frozen soil and supercooled water.
     - Force-restore snow has no density, heat capacity or aging and no vegetation masking.
   * - Water and sea ice
     - Columns with land mask 0 (``LANDMASK`` or ``erf.is_land``) and building footprints are
       skipped. No sea-ice test.
     - Cells whose land-use category is water (0, 17 or lake) are skipped; SLM sets this mask
       itself.
     - Cells the driver does not process (open water, sea ice) are skipped.
     - None of the three models a lake, ocean or sea-ice surface; the surface layer's own
       values apply there.
   * - Initial state
     - ``erf.rad_t_sfc`` and the ``seb_*_default`` values, or ``TSK`` / ``SKINTEMP``
       (``seb_init_skin_from_file``) and ``TMN`` (``seb_init_t_deep_from_file``, WRFInput only).
       Both need ``seb_prognostic_enable``; a restart uses the checkpoint.
     - Uniform ``slm.*`` values, WRFInput or metgrid.
     - A wrfinput-format setup file per level, named in ``namelist.erf``.
     - --
   * - Mesh refinement
     - Runs on every level whose boxes span the domain in z; a level that does not keeps the
       values it was built with. New and regridded levels are filled from the coarser level, and
       kept cells retain their values. With two-way coupling a fine level is averaged down (land
       columns only when water is present).
     - Runs on every level. New and regridded levels are interpolated from the coarser level
       (land and water separately), and kept cells retain their state. No average-down.
       Restart and regrid are tested on two levels.
     - Levels 0 to 2 can run their own driver with their own setup file; other levels take
       level 0's fields. No average-down. Level 0 cannot be remade (this stops a restart with
       ``erf.regrid_level_0_on_restart``, or with more ranks than level-0 boxes), and a level with
       its own driver cannot be regridded.
     - --
   * - GPU
     - Runs on the device.
     - Runs on the device.
     - Runs on the host (Fortran); fields are copied through pinned memory.
     - --
   * - Tests
     - Unit tests against closed forms (Noilhan and Planton 1989, Mahfouf and Noilhan 1996,
       Douville et al. 1995). CTests for soil moisture, hydrology, snow and the land mask; the
       multi-level, regrid, restart and surface-layer-flux CTests are registered but disabled.
     - Restart, regrid, two-level and new-level CTests (bitwise), and a frozen-canopy abort test.
       The gold-file cases need ``ERF_TEST_ENABLE_EXTRA_LSM_TESTS`` and files from
       ``ERF_TEST_EXTRA_FILES_DIRECTORY``.
     - ``Exec/RegTests/NoahMP_Ideal`` (not a CTest). ``TwoStream_NoahMPLevels`` is registered
       only when built with ``ERF_ENABLE_NOAHMP`` (no CI job has parallel NetCDF), and is
       disabled with the other two-stream CTests.
     - None of the three is compared with observations in ERF.

The case ``Exec/CanonicalTests/Radiation/TwoStream_NoahMP_vs_ForceRestore`` runs the
force-restore balance and Noah-MP side by side under the same two-stream radiation, on
grassland, dry bare soil and moist bare soil. Against Noah-MP with Jarvis stomata (option 2),
the force-restore skin temperature is within about 2 K, its sensible heat within about 11 %, its
latent heat within about 5 %, and the boundary-layer depth within about 6 %. Against Noah-MP's
default Ball-Berry stomata (option 1), the force-restore latent heat over grass near its peak
(3.8 h) is 2.4 times higher. Over the dry bare soil both latent heat fluxes are below
5 W/m². The case uses a soil-derived heat capacity (about 1.2e5 J/m²/K), not
the 2.0e4 J/m²/K default. No case
compares SLM with either model.

Cloud Fraction Diagnosis
--------------------------------


When ``cloud_fraction_prog_enable = true``, the cloud fraction is diagnosed at each level from the
relative humidity (RH) and cloud liquid water content (qc):

.. math::

   C_f = \min \left( 1, \max \left( 0, \frac{RH - \text{rh_min}}{\text{rh_max} - \text{rh_min}} \right) + \min \left( 1, \frac{q_c}{q_{c,\text{scale}}} \right) \right)

where the RH threshold parameters allow for a transition from clear sky (RH < rh_min) to complete cloud coverage
(RH >= rh_max). The relative humidity is :math:`RH = \min(1, q_v / q_{sat})`, with the saturation
mixing ratio over liquid water

.. math::

   q_{sat} = \frac{\epsilon \, e_{sat}}{\max(e_{sat},\, p - e_{sat})}, \qquad
   e_{sat} = 611.2 \, \text{Pa} \, \exp \left( \frac{17.62 \, T_c}{T_c + 243.12} \right),

where :math:`\epsilon = 0.622`, :math:`p` is the pressure and :math:`T_c = T - 273.15` is the temperature
in Celsius (Magnus formula, WMO 2008 coefficients). :math:`q_{c,\text{scale}}`
(``erf.radiation.cloud_fraction_qc_scale``, default :math:`1 \times 10^{-3}` kg/kg) is the cloud water at
which the liquid-water term alone reaches 1.

This :math:`C_f` scales the prescribed cloud-band optical depth (``tau_profile_type = cloud_layer``)
level by level inside the cloudy column. It does not replace ``erf.radiation.cloud_fraction``,
which sets the clear/cloudy blend; with ``cloud_fraction = 0`` the cloudy column is not computed.

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
``lon_m`` fields when a WRF or metgrid initialization filled them, and ``erf.rad_cons_lat`` and
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
