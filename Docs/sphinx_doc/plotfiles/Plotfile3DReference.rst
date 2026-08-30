.. _sec:Plotfile3DReference:

3D plotfile variable reference
==============================

Use ``erf.plot_vars_1`` and ``erf.plot_vars_2`` to select cell-centered and
derived fields for the two standard 3D plotfile streams. ERF writes each stream
on all active AMR levels. Some fields require a compile-time option, a physics
package, or stored diagnostic state; the tables below state those restrictions.

Time-averaged velocity fields require ``erf.time_avg_vel = true``. Set
``erf.plot_face_vels = true`` to write the native staggered velocity components
in separate face-centered outputs associated with the configured 3D stream.
Interval mean, covariance, and resolved TKE fields require
``erf.compute_mean_vars = true``; their averaging window is controlled by
``erf.mean_vars_reset_mode`` and ``erf.mean_vars_reset_time``.
In ``plotfile`` reset mode the window is global to both 3-D streams: every
stream due at the same simulation time receives the same accumulated values,
then ERF resets the accumulator once after the complete output batch. A due
stream that does not select an interval diagnostic (for example, a
``density``-only stream) does not reset the shared window.

3D output variables
-------------------

The fixed core names are selected from the input list and emitted in the order
defined by ERF. Some names require a compile-time option, a physics package, or
stored diagnostic state; those restrictions are noted below.

The default subvolume inventory is documented on :ref:`sec:Plotfiles`.

+-----------------------------+------------------+
| Variable                    | Definition       |
|                             |                  |
+=============================+==================+
| **x_velocity**              | Velocity in x    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **y_velocity**              | Velocity in y    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **z_velocity**              | Velocity in z    |
|                             | direction        |
|                             | [m/s]            |
+-----------------------------+------------------+
| **density**                 | Dry density      |
|                             | [kg/m^3]         |
|                             |                  |
+-----------------------------+------------------+
| **moist_density**           | Dry-air density  |
|                             | plus vapor and   |
|                             | non-precipitating|
|                             | condensate       |
|                             | density [kg/m^3] |
|                             |                  |
+-----------------------------+------------------+
| **dens_hse**                | Hydrostatic      |
|                             | density          |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **pert_dens**               | Perturbational   |
|                             | density          |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **pressure**                | Total pressure   |
|                             | [Pa]             |
|                             |                  |
+-----------------------------+------------------+
| **pres_hse**                | Hydrostatic      |
|                             | pressure         |
|                             | [Pa]             |
+-----------------------------+------------------+
| **theta_hse**               | Hydrostatic      |
|                             | potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **pi_hse**                  | Hydrostatic      |
|                             | Exner function   |
|                             | [-]              |
+-----------------------------+------------------+
| **qv_hse**                  | Base-state water |
|                             | vapor mixing     |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **pert_pres**               | Perturbational   |
|                             | pressure         |
|                             | [Pa]             |
+-----------------------------+------------------+
| **pres_hse_x**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in x    |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **pres_hse_y**              | Derivative of    |
|                             | hydrostatic      |
|                             | pressure in y    |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdx**                    | Pressure gradient|
|                             | in x direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdy**                    | Pressure gradient|
|                             | in y direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **dpdz**                    | Pressure gradient|
|                             | in z direction   |
|                             | [Pa/m]           |
+-----------------------------+------------------+
| **temp**                    | Temperature      |
|                             | [K]              |
|                             |                  |
+-----------------------------+------------------+
| **theta**                   | Potential        |
|                             | temperature [K]  |
|                             |                  |
+-----------------------------+------------------+
| **buoyancy**                | Buoyancy term    |
|                             | used by the      |
|                             | z-momentum       |
|                             | equation         |
+-----------------------------+------------------+
| **eq_pot_temp**             | Equivalent       |
|                             | potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **VPD**                     | Vapor pressure   |
|                             | deficit [kPa]    |
|                             |                  |
+-----------------------------+------------------+
| **rhotheta**                | Density * theta  |
|                             | [kg K/m^3]       |
|                             |                  |
+-----------------------------+------------------+
| **KE**                      | Prognostic       |
|                             | turbulent kinetic|
|                             | energy when      |
|                             | populated by the |
|                             | active closure,  |
|                             | including native |
|                             | SHOC [m^2/s^2]   |
+-----------------------------+------------------+
| **rhoKE**                   | Density * KE     |
|                             | [kg/(m s^2)]     |
|                             |                  |
+-----------------------------+------------------+
| **scalar**                  | Scalar magnitude |
|                             | [problem-dep.]   |
|                             |                  |
+-----------------------------+------------------+
| **reflectivity**            | reflectivity     |
|                             | cell-by-cell     |
|                             | [dBZ]            |
+-----------------------------+------------------+
| **max_reflectivity**        | max of           |
|                             | reflectivity     |
|                             | over a column    |
|                             | [dBZ]            |
+-----------------------------+------------------+
| **precipitable**            | precipitable     |
|                             | water (integral  |
|                             | over column)     |
|                             | [kg/m^2]         |
+-----------------------------+------------------+
| **mucape**                  | most unstable    |
|                             | CAPE over a      |
|                             | column [J/kg]    |
+-----------------------------+------------------+
| **vorticity_x**             | x-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_y**             | y-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **vorticity_z**             | z-component of   |
|                             | vorticity [1/s]  |
|                             |                  |
+-----------------------------+------------------+
| **local_helicity**          | helicity         |
|                             | cell-by-cell     |
|                             | [m/s^2]          |
+-----------------------------+------------------+
| **helicity**                | helicity         |
|                             | (integral over   |
|                             | column)          |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **magvel**                  | magnitude of     |
|                             | velocity [m/s]   |
|                             |                  |
+-----------------------------+------------------+
| **divU**                    | divergence of    |
|                             | velocity [1/s]   |
|                             |                  |
+-----------------------------+------------------+
| **u_t_avg**                 | time average of  |
|                             | x-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **v_t_avg**                 | time average of  |
|                             | y-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **w_t_avg**                 | time average of  |
|                             | z-component of   |
|                             | velocity [m/s]   |
+-----------------------------+------------------+
| **umag_t_avg**              | time average of  |
|                             | velocity mag     |
|                             | [m/s]            |
+-----------------------------+------------------+
| **u_mean**                  | interval mean of |
|                             | x velocity [m/s] |
+-----------------------------+------------------+
| **v_mean**                  | interval mean of |
|                             | y velocity [m/s] |
+-----------------------------+------------------+
| **w_mean**                  | interval mean of |
|                             | z velocity [m/s] |
+-----------------------------+------------------+
| **theta_mean**              | interval mean of |
|                             | potential temp.  |
|                             | [K]              |
+-----------------------------+------------------+
| **uu_mean**                 | interval mean of |
|                             | u squared        |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **vv_mean**                 | interval mean of |
|                             | v squared        |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **ww_mean**                 | interval mean of |
|                             | w squared        |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **uw_mean**                 | interval mean of |
|                             | u times w        |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **vw_mean**                 | interval mean of |
|                             | v times w        |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **wtheta_mean**             | interval mean of |
|                             | w times theta    |
|                             | [K m/s]          |
+-----------------------------+------------------+
| **uu_fluct**                | resolved u       |
|                             | variance         |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **vv_fluct**                | resolved v       |
|                             | variance         |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **ww_fluct**                | resolved w       |
|                             | variance         |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **uw_fluct**                | resolved u-w     |
|                             | covariance       |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **vw_fluct**                | resolved v-w     |
|                             | covariance       |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **wtheta_fluct**            | resolved w-theta |
|                             | covariance       |
|                             | [K m/s]          |
+-----------------------------+------------------+
| **tke_resolved**            | resolved TKE     |
|                             | [m^2/s^2]        |
+-----------------------------+------------------+
| **rhoadv_0**                | Conserved scalar |
|                             | [problem-dep.]   |
|                             |                  |
+-----------------------------+------------------+
| **soundspeed**              | Sound speed      |
|                             | [m/s]            |
|                             |                  |
+-----------------------------+------------------+
| **z_phys**                  | Terrain height   |
|                             | [m]              |
|                             |                  |
+-----------------------------+------------------+
| **detJ**                    | Jacobian         |
|                             | determinant [1]  |
|                             |                  |
+-----------------------------+------------------+
| **h_xi**                    | Cell-centered    |
|                             | average of the   |
|                             | metric term      |
|                             | dz/dxi [1]       |
+-----------------------------+------------------+
| **h_eta**                   | Cell-centered    |
|                             | average of the   |
|                             | metric term      |
|                             | dz/deta [1]      |
+-----------------------------+------------------+
| **h_zeta**                  | Cell-centered    |
|                             | average of the   |
|                             | metric term      |
|                             | dz/dzeta [1]     |
+-----------------------------+------------------+
| **mapfac**                  | Map scale factor |
|                             | [1]              |
|                             |                  |
+-----------------------------+------------------+
| **lat_m**                   | Latitude at mass |
|                             | points           |
|                             | [deg]            |
+-----------------------------+------------------+
| **lon_m**                   | Longitude at     |
|                             | mass points      |
|                             | [deg]            |
+-----------------------------+------------------+
| **nut**                     | Eddy viscosity,  |
|                             | nu_t [m^2/s]     |
+-----------------------------+------------------+
| **Kmv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Kmh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Momentum      |
|                             | (Note: For LES,  |
|                             | this is the      |
|                             | _dynamic_ eddy   |
|                             | viscosity, mu_t  |
|                             | = rho * nu_t     |
|                             | and Kmh==Kmv)    |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Khv**                     | Vertical         |
|                             | Eddy Diffusivity |
|                             | of Heat          |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Khh**                     | Horizontal       |
|                             | Eddy Diffusivity |
|                             | of Heat          |
|                             | [kg/(m s)]       |
+-----------------------------+------------------+
| **Lturb**                   | Turbulence       |
|                             | length scale     |
|                             | with             |
|                             | ``use_kturb``    |
|                             | [m]              |
+-----------------------------+------------------+
| **walldist**                | Wall distance    |
|                             | for RANS models  |
|                             | only [m]         |
+-----------------------------+------------------+
| **diss**                    | Subfilter-scale  |
|                             | dissipation      |
|                             | with diffusion / |
|                             | turbulence       |
|                             | [kg/(m s^3)]     |
+-----------------------------+------------------+
| **qt**                      | Total water      |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qn**                      | Nonprecipitating |
|                             | water (qv + qc + |
|                             | qi)              |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qp**                      | Precipitating    |
|                             | water (rain +    |
|                             | snow + graupel)  |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qc**                      | Cloud water      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qi**                      | Cloud ice        |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qrain**                   | Rain-water       |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qsnow**                   | Snow mixing      |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **qgraup**                  | Graupel mixing   |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **qv**                      | Water vapor      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **qsat**                    | Saturation water |
|                             | vapor mixing     |
|                             | ratio [kg/kg]    |
+-----------------------------+------------------+
| **rain_accum**              | Accumulated rain |
|                             | amount with      |
|                             | precipitating    |
|                             | moisture models  |
|                             | [mm]             |
+-----------------------------+------------------+
| **snow_accum**              | Accumulated snow |
|                             | amount with SAM  |
|                             | or Morrison      |
|                             | microphysics     |
|                             | [mm]             |
+-----------------------------+------------------+
| **graup_accum**             | Accumulated      |
|                             | graupel amount   |
|                             | with SAM or      |
|                             | Morrison         |
|                             | microphysics     |
|                             | [mm]             |
+-----------------------------+------------------+
| **rel_humidity**            | Relative         |
|                             | humidity;        |
|                             | currently filled |
|                             | only for         |
|                             | SuperDroplets    |
|                             | [1]              |
+-----------------------------+------------------+
| **condensation_rate**       | Condensation     |
|                             | rate with        |
|                             | SuperDroplets    |
|                             | only             |
|                             | [kg/kg/s]        |
+-----------------------------+------------------+
| **terrain_IB_mask**         | Immersed-boundary|
|                             | terrain/building |
|                             | mask; available  |
|                             | for immersed     |
|                             | forcing terrain  |
|                             | or buildings     |
|                             | [1]              |
+-----------------------------+------------------+
| **volfrac**                 | EB / immersed    |
|                             | boundary volume  |
|                             | fraction; unity  |
|                             | elsewhere        |
|                             | [1]              |
+-----------------------------+------------------+
| **qsrc_sw**                 | Shortwave        |
|                             | radiative        |
|                             | heating source   |
|                             | term with        |
|                             | radiation        |
|                             | [K/s]            |
+-----------------------------+------------------+
| **qsrc_lw**                 | Longwave         |
|                             | radiative        |
|                             | heating source   |
|                             | term with        |
|                             | radiation        |
|                             | [K/s]            |
+-----------------------------+------------------+
| **tracer_particles_count**  | Tracer particle  |
|                             | count per cell   |
|                             | requires         |
|                             | ERF_USE_PARTICLES|
|                             | to be defined    |
|                             | [count]          |
+-----------------------------+------------------+

The ``qrain``, ``qsnow``, and ``qgraup`` rows are available when the active
moisture scheme provides the corresponding rain, snow, or graupel component.

Moisture variable selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every moisture plot variable -- the individual species, the aggregates, the
surface accumulations and the moist diagnostics -- is selected from a single
per-scheme index map, ``SolverChoice::moisture_indices`` (the
``MoistureComponentIndices`` struct in ``Source/DataStructs/ERF_DataStruct.H``).
A variable is written when, and only when, that map names the data behind it.
The same map is used by the writer, so a name in the plotfile header and the
data under it come from one place.

The criterion is a valid index, not an allocated component. Several schemes
allocate a wider moist state than they integrate -- the ``Morrison`` class always
allocates eleven moist components, so ``Morrison_NoIce`` owns ice slots it never
fills -- and a plotfile must not publish untouched storage as if it were data.
Requests for a variable the active scheme does not carry are ignored rather than
treated as an error, and the name does not reserve a plotfile component.

Two kinds of index live in the map:

* **Conserved-state components** for the species, i.e. one of ``RhoQ1`` through
  ``RhoQ11``. These name the density-weighted species, so the writer divides by
  the density to output a mixing ratio or a number per unit mass.
* **qmoist slots** for the surface accumulations and moist diagnostics. Each
  scheme chooses that layout for itself, so the slot numbers differ between
  schemes and are recorded per scheme in the map.

The species carried by each scheme, and the conserved component each one
occupies:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Moisture model
     - Mass species
     - Number concentrations
   * - ``None``
     - none
     - none
   * - ``MoistNoCondensation``
     - qv (Q1), qc (Q2)
     - none
   * - ``SatAdj``
     - qv (Q1), qc (Q2)
     - none
   * - ``Kessler_NoRain``
     - qv (Q1), qc (Q2)
     - none
   * - ``Kessler``
     - qv (Q1), qc (Q2), qr (Q3)
     - none
   * - ``SAM_NoPrecip_NoIce``
     - qv (Q1), qc (Q2)
     - none
   * - ``SAM_NoIce``
     - qv (Q1), qc (Q2), qr (Q4)
     - none
   * - ``SAM``
     - qv (Q1), qc (Q2), qi (Q3), qr (Q4), qs (Q5), qg (Q6)
     - none
   * - ``Morrison_NoIce``
     - qv (Q1), qc (Q2), qr (Q4)
     - nc (Q7), nr (Q9)
   * - ``Morrison``
     - qv (Q1), qc (Q2), qi (Q3), qr (Q4), qs (Q5), qg (Q6)
     - nc (Q7), ni (Q8), nr (Q9), ns (Q10), ng (Q11)
   * - ``WSM6``
     - qv (Q1), qc (Q2), qi (Q3), qr (Q4), qs (Q5), qg (Q6)
     - none
   * - ``WDM6``
     - qv (Q1), qc (Q2), qi (Q3), qr (Q4), qs (Q5), qg (Q6)
     - nc (Q7), nn (Q8), nr (Q9)
   * - ``SuperDroplets``
     - qv (Q1), qc (Q2), qi (Q3), qr (Q4), qs (Q5), qg (Q6)
     - none

Note that ``nn``, the CCN / total aerosol number, is not a hydrometeor count:
it is an aerosol reservoir with no companion mass species, and it takes the
conserved component that ``Morrison`` uses for cloud ice number. The two are
mutually exclusive by scheme, and the map keeps them apart, so ``ni`` and
``nn`` are never both selectable.

The accumulations and moist diagnostics available from each scheme's qmoist
arrays:

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Moisture model
     - Accumulations
     - Moist diagnostics
   * - ``Kessler``, ``SAM_NoIce``, ``Morrison_NoIce``
     - ``rain_accum``
     - none
   * - ``SAM``, ``Morrison``, ``WSM6``, ``WDM6``
     - ``rain_accum``, ``snow_accum``, ``graup_accum``
     - none
   * - ``SatAdj``
     - none
     - ``rel_humidity``, recovered from the conserved state
   * - ``SuperDroplets``
     - ``rain_accum``, ``snow_accum``
     - ``rel_humidity``, ``condensation_rate``
   * - all others
     - none
     - none

``SuperDroplets`` allocates a graupel accumulation slot that nothing fills, so
``graup_accum`` is not offered for that scheme. ``SatAdj`` publishes no qmoist
arrays at all: its ``rel_humidity`` is derived from the conserved state when the
plotfile is written, which the map records with a distinct sentinel so the
writer takes the right path.

The aggregates are sums over the species the map hands out, so they follow the
table above rather than a fixed component range:

* ``qt`` sums every mass species. Number concentrations are counts, not masses,
  and are excluded.
* ``qn`` sums vapor and the suspended condensate: ``qv``, ``qc`` and, where the
  scheme carries it, ``qi``.
* ``qp`` sums the falling species: ``qr``, ``qs`` and ``qg``.
* ``moist_density`` adds the ``qn`` sum to the dry density.
* ``qsat`` and ``precipitable`` read vapor only, and so need any moist scheme.

For example ``Kessler`` gives ``qt = qv + qc + qr``, ``qn = qv + qc`` and
``qp = qr``, while ``Morrison_NoIce`` gives the same three sums even though it
allocates the frozen species alongside them.

Raw state requests behave the same way: ``rhoQ1`` through ``rhoQ11`` are
selected when the map names that component for the active scheme, so
``Morrison_NoIce`` offers ``rhoQ4`` and ``rhoQ7`` but not ``rhoQ3`` or
``rhoQ8``. The dry state names remain bounded by the allocated state width.

Terrain metric restrictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mesh geometry diagnostics ``z_phys``, ``detJ``, ``h_xi``, ``h_eta``, and
``h_zeta`` are selected only when the mesh is terrain-fitted, that is when
``erf.terrain_type = StaticFittedMesh`` or
``erf.terrain_type = MovingFittedMesh``. A request for one of these names under
any other terrain type is ignored rather than treated as an error, and the name
does not reserve a plotfile component.

``h_xi``, ``h_eta``, and ``h_zeta`` are the cell-centered averages of the
terrain metric terms

.. math::

   h_\xi = \frac{\partial z}{\partial \xi}, \qquad
   h_\eta = \frac{\partial z}{\partial \eta}, \qquad
   h_\zeta = \frac{\partial z}{\partial \zeta},

where :math:`z` is the physical height stored at mesh nodes and
:math:`(\xi, \eta, \zeta)` are the computational coordinates. Each value is
formed by averaging the four nodal differences that bracket the cell, so
``h_xi`` and ``h_eta`` vanish on a mesh with no horizontal terrain variation and
``h_zeta`` reduces to the ratio of the physical to the computational cell height.
Because ``detJ`` equals ``h_zeta`` for the terrain-fitted mapping used by ERF,
these three fields together give the full metric Jacobian of the mapping.

Optional storage restrictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following fixed diagnostics are selected only if the storage exists on
every AMR level in the plotfile:

* ``u_t_avg``, ``v_t_avg``, ``w_t_avg``, and ``umag_t_avg`` require
  ``erf.time_avg_vel = true``. If no samples have been accumulated yet, the
  output value is defined as zero rather than dividing by zero. The running
  sum and its normalizer are saved in checkpoint files, so the averaging
  window survives a restart; see :ref:`sec:Checkpoint`.

  The ``*_t_avg`` fields are a long-running cumulative mean: the averaging
  window never resets on its own and simply grows for the length of the run
  (surviving restarts via the checkpoint). This is distinct from the
  ``*_mean``/``*_fluct`` interval diagnostics below, whose averaging window is
  short and resets periodically; use ``*_t_avg`` for a run-long mean flow and
  ``*_mean``/``*_fluct`` for turbulence statistics over a controlled window.
* ``u_mean``, ``v_mean``, ``w_mean``, ``theta_mean``, all ``*_mean`` second
  moments, all ``*_fluct`` resolved variances/covariances, and ``tke_resolved`` require
  ``erf.compute_mean_vars = true``. For example,
  ``uw_fluct = uw_mean - u_mean*w_mean`` and
  ``tke_resolved = 0.5*(uu_fluct + vv_fluct + ww_fluct)``. These fields are also zero
  before the first sample is accumulated.

  Note that the accumulator for a given level is rebuilt when that level is
  regridded, so a regrid restarts the averaging window on the regridded level
  while the other levels continue accumulating. In a multi-level run a
  plotfile written shortly after a regrid can therefore contain interval
  fields whose averaging windows differ from level to level, and the plotfile
  itself carries no record of the per-level window length. The same is true of
  the ``*_t_avg`` fields above. If a consistent window across levels matters
  for the analysis, write the interval fields from a run without regridding,
  or allow enough time after a regrid for the window to refill.
* ``qsrc_sw`` and ``qsrc_lw`` require a non-``None`` radiation choice.
* ``nut``, ``Kmv``, ``Kmh``, ``Khv``, ``Khh``, and ``Lturb`` require
  ``use_kturb = true`` at every AMR level.
* ``diss`` requires molecular diffusion or ``use_kturb`` at every level.
* ``walldist`` requires a non-``None`` RANS choice at every level.

These restrictions apply to fixed names only. Dynamic names supplied by an
active microphysics or particle provider remain governed by that provider's
own plot-variable contract.

Fixed conserved-state fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conserved-state inventory also includes the moisture-density components
below. The active moisture model determines which components are selectable.

.. code-block:: text

   rhoQ1 rhoQ2 rhoQ3 rhoQ4 rhoQ5 rhoQ6
   rhoQ7 rhoQ8 rhoQ9 rhoQ10 rhoQ11

These are conserved moisture or species densities. The component names are
fixed; the active moisture model controls which components are available.
ERF checks the active conserved-state and microphysics sizes before selecting a
requested component. Unsupported fixed names are omitted and reported as
unavailable; they do not reserve a plotfile component. Provider-supplied names
remain dynamic and are appended only when the active provider exposes them.

Velocity output behavior
~~~~~~~~~~~~~~~~~~~~~~~~

Requesting any one of ``x_velocity``, ``y_velocity``, or ``z_velocity`` selects
all three cell-centered velocity components. When ``erf.plot_face_vels = true``,
native AMReX output also writes three separate face-centered plotfiles with the
components ``x_velocity_stag``, ``y_velocity_stag``, and ``z_velocity_stag``.

For stream 1, the face files use ``<plot_file_1>U<step>``,
``<plot_file_1>V<step>``, and ``<plot_file_1>W<step>``. Stream 2 uses the
corresponding ``<plot_file_2>U<step>``, ``V<step>``, and ``W<step>`` prefixes.
These face-velocity filenames always use level-0 step numbering, including when
``erf.use_real_time_in_pltname = true``. The face outputs are native AMReX
artifacts; the NetCDF writer does not enter this face-output path.

Supplemental fixed derived fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main table contains the ordinary fixed fields. The following fixed names
are conditional or specialized and remain part of the same source-defined
inventory. A requested name is still subject to the runtime selection
conditions in ``setPlotVariables``.

The following fixed fields are supplied by the native SHOC driver. The
native-only fields are written as ``-999`` when native SHOC diagnostics are not
available.

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - Variable
     - Units
     - Native SHOC meaning
   * - ``pblh``
     - m
     - Native SHOC ``pblh`` is reported in metres above local ground (AGL).
   * - ``shoc_cldfrac``
     - 1
     - Subgrid cloud fraction diagnosed by the native SHOC PDF.
   * - ``shoc_ql``
     - kg/kg
     - Cloud-liquid mixing ratio diagnosed by the native SHOC PDF.
   * - ``shoc_ql2``
     - (kg/kg)^2
     - Variance of the diagnosed cloud-liquid mixing ratio.
   * - ``shoc_cond``
     - kg/kg/s
     - Positive PDF condensate-change rate. It is zero unless ``erf.shoc.extra_shoc_diags = true``.
   * - ``wqls_sec``
     - (kg/kg) m/s
     - Turbulent vertical flux of the diagnosed SHOC cloud-liquid quantity.
   * - ``wthv_sec``
     - K m/s
     - Turbulent vertical flux of virtual potential temperature.
   * - ``w_sec``
     - m^2/s^2
     - Vertical-velocity variance.
   * - ``thl_sec``
     - K^2
     - Liquid-water potential-temperature variance.
   * - ``qw_sec``
     - (kg/kg)^2
     - Total-water variance.
   * - ``qwthl_sec``
     - K kg/kg
     - Covariance of total water and liquid-water potential temperature.
   * - ``wthl_sec``
     - K m/s
     - Turbulent vertical flux of liquid-water potential temperature.
   * - ``wqw_sec``
     - (kg/kg) m/s
     - Turbulent vertical flux of total water.
   * - ``w3``
     - m^3/s^3
     - Third central moment of vertical velocity.
   * - ``brunt``
     - s^-2
     - Squared Brunt-Vaisala frequency, :math:`N^2`.
   * - ``isotropy``
     - s
     - Native SHOC isotropy timescale.
   * - ``shear_prod``
     - m^2/s^3
     - Shear-production contribution to the native SHOC TKE budget.
   * - ``buoy_prod``
     - m^2/s^3
     - Buoyancy production or destruction contribution to the native SHOC TKE budget.
   * - ``diss_tke``
     - m^2/s^3
     - Dissipation contribution diagnosed for the native SHOC TKE budget.

When native SHOC diagnostics are available, the standard turbulence fields use
the native SHOC coefficients: ``Kmv`` is :math:`\rho K_m`, ``Khv`` is
:math:`\rho K_h`, ``Lturb`` is the native SHOC mixing length, and ``nut`` is
the kinematic momentum diffusivity ``Kmv / density``.

* ``nc``, ``ni``, ``nn``, ``nr``, ``ns``, and ``ng`` are moisture number
  concentrations. They are available when the active moisture model provides
  the corresponding conserved component; the output units follow that model's
  number-concentration convention. ``nn`` is the CCN / total aerosol number
  carried by ``WDM6``; it shares a conserved component with ``ni``, so the two
  are never selectable at the same time.
* ``xvel_err``, ``yvel_err``, ``zvel_err``, and ``pp_err`` are error diagnostics
  available when ``ERF_COMPUTE_ERROR`` is enabled. Their units follow the
  corresponding velocity or pressure quantity.

The fixed conserved-state inventory is ``density``, ``rhotheta``, ``rhoKE``,
``rhoadv_0``, and ``rhoQ1`` through ``rhoQ11``. The ``rhoQ`` components are
active only to the extent that the selected moisture model provides them.

Wind-farm-only variables
~~~~~~~~~~~~~~~~~~~~~~~~

The following quantities are available only in builds with
``ERF_USE_WINDFARM`` enabled.

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **num_turb**                | Number of wind   |
|                             | turbines in cell |
|                             | for Fitch, EWP,  |
|                             | SimpleAD, and    |
|                             | GeneralAD        |
|                             | [count]          |
+-----------------------------+------------------+
| **SMark0**                  | Windfarm marker  |
|                             | component 0 for  |
|                             | Fitch, EWP,      |
|                             | SimpleAD, and    |
|                             | GeneralAD        |
|                             | [1]              |
+-----------------------------+------------------+
| **SMark1**                  | Windfarm marker  |
|                             | component 1 for  |
|                             | SimpleAD and     |
|                             | GeneralAD        |
|                             | [1]              |
+-----------------------------+------------------+

Morrison microphysics variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using Morrison two-moment microphysics, additional diagnostic variables
are available for output. These variables provide detailed information about
cloud and precipitation processes. To enable Morrison output, include any of
the variables below in your **erf.plot_vars_1** or **erf.plot_vars_2** list.

**Thermodynamic State Variables:**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_rho**               | Air density      |
|                             | [kg/m^3]         |
+-----------------------------+------------------+
| **micro_theta**             | Potential        |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **micro_temp**              | Absolute         |
|                             | temperature [K]  |
+-----------------------------+------------------+
| **micro_pres**              | Pressure [Pa]    |
|                             |                  |
+-----------------------------+------------------+

**Non-Precipitating Moisture Variables (mixing ratios in kg/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_qv**                | Water vapor      |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qc**                | Cloud liquid     |
|                             | water mixing     |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qi**                | Cloud ice        |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qn**                | Total cloud      |
|                             | condensate       |
|                             | (qc + qi)        |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qt**                | Total water      |
|                             | mixing ratio     |
|                             | (qv + qn)        |
|                             | [kg/kg]          |
+-----------------------------+------------------+

**Precipitating Hydrometeor Variables (mixing ratios in kg/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_qp**                | Total            |
|                             | precipitation    |
|                             | (qrain + qsnow + |
|                             | qgraup)          |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qrain**             | Rain water       |
|                             | mixing ratio     |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qsnow**             | Snow mixing      |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+
| **micro_qgraup**            | Graupel mixing   |
|                             | ratio            |
|                             | [kg/kg]          |
+-----------------------------+------------------+

**Number Concentrations (1/kg):**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_nc**                | Cloud droplet    |
|                             | number           |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_nr**                | Rain drop number |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ni**                | Cloud ice number |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ns**                | Snow number      |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+
| **micro_ng**                | Graupel number   |
|                             | concentration    |
|                             | [1/kg]           |
+-----------------------------+------------------+

**Dynamical Variables:**

+-----------------------------+------------------+
| Variable                    | Definition       |
+=============================+==================+
| **micro_omega**             | Grid-scale       |
|                             | vertical         |
|                             | velocity [m/s]   |
|                             | used as input to |
|                             | Morrison scheme  |
+-----------------------------+------------------+

**Example Usage:**

To output Morrison diagnostic variables, add them to your plot variables list:

.. code-block:: text

   erf.plot_vars_1 = density theta qv micro_qc micro_qrain micro_nc micro_nr

This will output the base ERF variables (density, theta, qv) along with Morrison
cloud water, rain water, cloud droplet number concentration, and rain drop number
concentration.

Scheme-provided dynamic fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the fixed names are selected, ERF asks the active microphysics provider for
additional plot names. There is no universal static list for these fields.

The current providers expose the following families:

* Morrison exposes the 19 ``micro_*`` names listed in the section above.
* SuperDroplets generates ``qv_<species>``, ``qc_<species>``,
  ``qt_<species>``, ``sat_ratio_<species>``, and ``accum_<species>`` names for
  configured species, plus ``accum_<aerosol>`` names for configured aerosols.
* Providers that inherit the empty ``NullMoist`` plot-name implementation add
  no scheme-specific names.

Particle-provided fields
~~~~~~~~~~~~~~~~~~~~~~~~

Particle builds can add Eulerian mesh fields after particle containers are
initialized. ERF prefixes each provider name with the particle-container name.
The default particle container provides ``<container>_mass_density``. The
SuperDroplet provider additionally generates mass-flux, number-density,
species-density, species-flux, aerosol-density, and aerosol-flux families from
the configured species and aerosols. Particle count names use the form
``<container>_count``. These fields are dynamic; their exact names depend on the
particle configuration. Every configured container can be selected, including
one that has not yet been allocated; its mesh count is then initialized to zero.

Before constructing a 3D plotfile in a Lagrangian-microphysics run with
two-way AMR coupling and at least two AMR levels, ERF averages fine-level
microphysics storage, ``RhoTheta``, and
the active moist conserved components onto covered coarse cells. This is a
consistency update so coarse output reflects fine-level deposits; it is not a
physical time tendency. Other 3D output assembly is diagnostic construction.
