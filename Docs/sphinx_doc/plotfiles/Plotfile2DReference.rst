.. _sec:Plotfile2DReference:

2D plotfile diagnostic reference
================================

ERF writes each 2D output level as a one-cell-thick horizontal slab. The two 2D
streams are independent. Each stream has its own prefix, interval, period,
built-in variable list, and sampled-level list.

Fixed diagnostics appear in canonical catalog order, not request order.
Sampled-level fields follow the fixed diagnostics and appear in level-set,
target-value, and field order. ERF warns and skips an unknown or
configuration-ineligible fixed name.

Three contracts govern a fixed diagnostic:

* **Descriptor metadata** defines its public name, canonical order, long name,
  exact unit string, category, and missing-value policy.
* **Selectability** determines whether ERF accepts the name for the active
  configuration and creates an output component.
* **Runtime value** determines what ERF writes at each cell and output time after
  the name has been selected.

A variable may appear in the catalog without being available in every
simulation. A metadata policy also does not guarantee a non-missing value at
every cell. Internal provider sentinels are translated at the public output
boundary. Categorical source code ``0`` means that no source was selected; it is
not the continuous ``-999`` fill value.

For the numerical production of the unified 2-m fields, see
:ref:`sec:LandSurfaceDiagnostics`. For fields sampled on pressure, height, or
model-index levels, see :ref:`sec:Plotfile2DSampledLevels`.

.. _sec:Plotfile2DBuiltInCatalog:

Built-in diagnostic catalog
---------------------------

The marked table is the authoritative fixed catalog. The three contracts in the
introduction apply to every row; the selection table below gives the
configuration and runtime details that cannot be inferred from metadata alone.

.. BEGIN ERF BUILT-IN 2D DIAGNOSTIC CATALOG

.. list-table::
   :header-rows: 1
   :widths: 30 17 20 38 80

   * - Variable
     - Category
     - Units
     - Metadata policy
     - Description
   * - ``z_surf``
     - ``Geometry``
     - ``m``
     - ``AlwaysAvailable``
     - Surface elevation
   * - ``landmask``
     - ``Geometry``
     - ``1``
     - ``AlwaysAvailable``
     - Land-sea mask
   * - ``mapfac``
     - ``Geometry``
     - ``1``
     - ``AlwaysAvailable``
     - Map factor at mass points
   * - ``lat_m``
     - ``Geometry``
     - ``deg``
     - ``FillZeroWhenUnavailable``
     - Latitude at unstaggered mass points
   * - ``lon_m``
     - ``Geometry``
     - ``deg``
     - ``FillZeroWhenUnavailable``
     - Longitude at unstaggered mass points
   * - ``u_star``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Friction velocity from the surface layer
   * - ``w_star``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Convective velocity scale from the surface layer
   * - ``t_star``
     - ``SurfaceLayer``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Temperature scale from the surface layer
   * - ``q_star``
     - ``SurfaceLayer``
     - ``kg/kg``
     - ``FillMinus999WhenUnavailable``
     - Humidity scale from the surface layer
   * - ``Olen``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Obukhov length from the surface layer
   * - ``pblh``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Planetary boundary layer height from the active PBL diagnostic provider
   * - ``t_surf``
     - ``SurfaceLayer``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Surface potential temperature from the surface layer
   * - ``q_surf``
     - ``SurfaceLayer``
     - ``kg/kg``
     - ``FillMinus999WhenUnavailable``
     - Surface humidity from the surface layer
   * - ``z0``
     - ``SurfaceLayer``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Roughness height from the surface layer
   * - ``u_10m``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Grid-relative x wind 10 m above the local surface, from the surface-layer profile
   * - ``v_10m``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Grid-relative y wind 10 m above the local surface, from the surface-layer profile
   * - ``u_east_10m``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Eastward wind 10 m above the local surface, from the surface-layer profile
   * - ``v_north_10m``
     - ``SurfaceLayer``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Northward wind 10 m above the local surface, from the surface-layer profile
   * - ``OLR``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Outgoing longwave radiation at the model top
   * - ``topographic_shadow``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Direct sunlight blocked by terrain or by a midpoint sun below the cutoff when averaged cosine exceeds it (1 blocked, 0 otherwise)
   * - ``terrain_direct_factor``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Terrain correction factor for direct shortwave
   * - ``terrain_incidence``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Nonnegative cosine of sunlight incidence on terrain
   * - ``unshaded_sw_flux_dn``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Surface downwelling shortwave before terrain correction
   * - ``terrain_sw_flux_dn``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Surface downwelling shortwave after terrain correction
   * - ``seb_t_sfc``
     - ``Radiation``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Prognostic surface temperature from the two-stream simplified surface energy balance
   * - ``seb_q_sfc``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Surface moisture of the two-stream surface energy balance: the soil water content (m3/m3) with seb_prognostic_enable, a land model's 2 m mixing ratio (kg/kg) with one
   * - ``seb_hfx``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Sensible heat flux in the two-stream simplified surface energy balance
   * - ``seb_lh``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Latent heat flux in the two-stream simplified surface energy balance
   * - ``seb_sw_net``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Net shortwave absorbed by the two-stream simplified surface energy balance at its last update after the dycore
   * - ``ibseb_ground_shadow``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Shadow of the immersed-boundary buildings on the ground for the two-stream surface energy balance on the sweep's shortwave: 1 shaded, 0 sunlit
   * - ``sens_flux``
     - ``SurfaceFlux``
     - ``kg K m^-2 s^-1``
     - ``FillMinus999WhenUnavailable``
     - Surface sensible heat flux
   * - ``laten_flux``
     - ``SurfaceFlux``
     - ``kg m^-2 s^-1``
     - ``FillMinus999WhenUnavailable``
     - Surface moisture flux (legacy output name)
   * - ``surf_pres``
     - ``SurfaceState``
     - ``Pa``
     - ``AlwaysAvailable``
     - Surface pressure, at the ground
   * - ``sea_level_pressure``
     - ``SurfaceState``
     - ``Pa``
     - ``AlwaysAvailable``
     - Sea-level pressure from the NMC/NGM reduction with Shuell correction
   * - ``precip_total_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface precipitation, liquid-water equivalent
   * - ``precip_rain_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface rain precipitation, liquid-water equivalent
   * - ``precip_snow_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface snow precipitation, liquid-water equivalent
   * - ``precip_graupel_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface graupel precipitation, liquid-water equivalent
   * - ``precip_hail_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated surface hail precipitation, liquid-water equivalent
   * - ``precip_frozen_accum``
     - ``Precipitation``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Accumulated frozen surface precipitation, liquid-water equivalent
   * - ``integrated_qv``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated water vapor
   * - ``integrated_qc``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated cloud liquid water
   * - ``integrated_qi``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated cloud ice
   * - ``integrated_qr``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated rain water
   * - ``integrated_qs``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated snow
   * - ``integrated_qg``
     - ``ColumnIntegral``
     - ``kg/m^2``
     - ``FillZeroWhenUnavailable``
     - Column-integrated graupel
   * - ``cloud_optical_depth_sw``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Estimated per-level shortwave hydrometeor extinction optical depth (including rain, snow, and graupel; may be partial-column)
   * - ``cloud_albedo_est``
     - ``Radiation``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Per-level hydrometeor-only reflectance proxy (including rain, snow, and graupel; excludes surface, gas, and aerosol contributions; may be partial-column)
   * - ``surface_diagnostic_source``
     - ``SurfaceLayer``
     - ``1``
     - ``AlwaysAvailable``
     - Surface diagnostic source code
   * - ``sensible_heat_flux``
     - ``SurfaceFlux``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Surface sensible heat flux
   * - ``latent_heat_flux``
     - ``SurfaceFlux``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Surface latent heat flux
   * - ``shoc_u_star``
     - ``PBL``
     - ``m/s``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC friction velocity diagnostic
   * - ``shoc_Olen``
     - ``PBL``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC Obukhov length diagnostic
   * - ``shoc_wthv_sfc``
     - ``PBL``
     - ``K m s^-1``
     - ``FillMinus999WhenUnavailable``
     - Native SHOC surface virtual potential temperature flux
   * - ``t_sfc``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP radiative surface temperature
   * - ``sfc_emis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Surface bulk emissivity
   * - ``sfc_alb_dir_vis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Direct visible surface albedo
   * - ``sfc_alb_dir_nir``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Direct near-infrared surface albedo
   * - ``sfc_alb_dif_vis``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Diffuse visible surface albedo
   * - ``sfc_alb_dif_nir``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Diffuse near-infrared surface albedo
   * - ``cos_zenith_angle``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Cosine of the solar zenith angle supplied to Noah-MP
   * - ``sw_flux_dn``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Downwelling shortwave flux supplied to Noah-MP
   * - ``sw_flux_dn_dir_vis``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Direct visible downwelling shortwave flux
   * - ``sw_flux_dn_dir_nir``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Direct near-infrared downwelling shortwave flux
   * - ``sw_flux_dn_dif_vis``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Diffuse visible downwelling shortwave flux
   * - ``sw_flux_dn_dif_nir``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Diffuse near-infrared downwelling shortwave flux
   * - ``lw_flux_dn``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Downwelling longwave flux supplied to Noah-MP
   * - ``grdflx``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Ground or snow heat flux
   * - ``fira``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Total net longwave flux, positive toward the atmosphere
   * - ``sav``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Solar radiation absorbed by vegetation
   * - ``sag``
     - ``LandSurface``
     - ``W m^-2``
     - ``FillMinus999WhenUnavailable``
     - Solar radiation absorbed by the ground
   * - ``albedo``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Broadband surface albedo
   * - ``sfcrunoff``
     - ``LandSurface``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Accumulated surface runoff
   * - ``udrunoff``
     - ``LandSurface``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Accumulated subsurface runoff
   * - ``snow``
     - ``LandSurface``
     - ``kg m^-2``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP snow water equivalent
   * - ``snowh``
     - ``LandSurface``
     - ``m``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP physical snow depth
   * - ``noahmp_temperature_2m_vegetated``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m temperature over the vegetated fraction
   * - ``noahmp_temperature_2m_bare``
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m temperature over the bare fraction
   * - ``noahmp_water_vapor_mixing_ratio_2m_vegetated``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m water-vapor mixing ratio over the vegetated fraction
   * - ``noahmp_water_vapor_mixing_ratio_2m_bare``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP 2-m water-vapor mixing ratio over the bare fraction
   * - ``noahmp_vegetation_fraction``
     - ``LandSurface``
     - ``1``
     - ``FillMinus999WhenUnavailable``
     - Noah-MP vegetation fraction
   * - ``temperature_2m``
     - ``NearSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Physical air temperature 2 m above the local surface
   * - ``water_vapor_mixing_ratio_2m``
     - ``NearSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Water-vapor mixing ratio 2 m above the local surface per unit dry-air mass
   * - ``near_surface_diagnostic_source``
     - ``NearSurface``
     - ``1``
     - ``AlwaysAvailable``
     - Source code for the unified near-surface diagnostic bundle
.. END ERF BUILT-IN 2D DIAGNOSTIC CATALOG

.. note::

   For the ``landmask`` variable, land is ``1`` and sea is ``0``. Buildings are
   ``2`` when using ``ImmersedForcing`` on levels where
   :math:`\sqrt{\Delta x\,\Delta y} \leq 50\,\mathrm{m}`, matching the building
   forcing cutoff. Immersed terrain retains its ``2`` mask at all resolutions.

If a requested fixed diagnostic is not selectable for the active configuration,
ERF warns and skips it. The descriptor table above records metadata; it does
not define request acceptance or guarantee a non-missing runtime value.

.. _sec:Plotfile2DSelectionRules:

Selection and runtime-value rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The selection contract and the value written after selection are separate:

.. list-table::
   :header-rows: 1
   :widths: 30 35 55

   * - Name or family
     - Selectable when
     - Runtime value
   * - ``z_surf``, ``landmask``, ``mapfac``
     - Selectable: fixed geometry or state writer path.
     - Value: ERF writes the corresponding geometry or state value.
   * - ``lat_m``, ``lon_m``
     - Selectable: fixed request names.
     - Value: ERF writes the coordinate data, or zero when the coordinate source is absent.
   * - ``u_star``, ``w_star``, ``t_star``, ``q_star``, ``Olen``, ``t_surf``, ``q_surf``, ``z0``
     - Selectable: ERF accepts these fixed request names; the runtime
       SurfaceLayer source may still be absent.
     - Value: SurfaceLayer scalar values; ``-999`` when the source pointer is absent.
   * - ``u_10m``, ``v_10m``, ``u_east_10m``, ``v_north_10m``
     - Selectable: fixed request names; the runtime SurfaceLayer source may still be absent.
     - Value: the wind 10 m above the local surface, as WRF's ``U10``/``V10``. Its speed is
       :math:`u_* / \kappa \, [\ln(10/z_0) - \psi_m(10/L)]` from the surface layer's
       :math:`u_*`, :math:`z_0` and Obukhov length :math:`L` (with :math:`10/L` capped at 1,
       as for ``temperature_2m``), and its direction that of the wind at the centre of the
       first cell. ``u_east_10m`` and ``v_north_10m`` are turned earth-relative with the
       grid orientation from the run's latitude and longitude, as the sampled ``u_east``
       and ``v_north`` are; without them the grid is taken to be aligned with east. ``-999``
       without a surface layer, or where :math:`z_0 \ge 10` m.
   * - ``pblh``
     - Selectable: fixed request name.
     - Native SHOC ``pblh`` is reported in metres above local ground (AGL). Value: native SHOC PBL height when native SHOC diagnostics are present; otherwise SurfaceLayer; ``-999`` if neither exists.
   * - ``OLR``
     - Selectable: fixed request name.
     - Value: radiation output; ``-999`` when the radiation source is absent.
   * - ``topographic_shadow``, ``terrain_direct_factor``, ``terrain_incidence``, ``unshaded_sw_flux_dn``, ``terrain_sw_flux_dn``
     - Selectable: when either ``erf.topo_shading`` or ``erf.slope_rad`` is enabled.
     - Value: diagnostics of the RRTMGP terrain correction, held between radiation
       updates. ``-999`` before the first radiation update. See
       :ref:`sec:Plotfile2DTopographicRadiation`.
   * - ``seb_t_sfc``, ``seb_q_sfc``
     - Selectable: fixed request names.
     - Value: the two-stream simplified surface energy balance's own prognostic surface
       state, which is distinct from Noah-MP's ``t_sfc`` and from the surface layer's
       ``t_surf``. ``-999`` unless the two-stream solver is running with
       ``erf.radiation.seb_enable = true``.
   * - ``seb_hfx``, ``seb_lh``
     - Selectable: fixed request names.
     - Value: the sensible and latent heat fluxes the two-stream surface energy balance
       used at its last update, positive away from the surface. With
       ``erf.radiation.seb_turbulent_flux_source = surface_layer`` (the default) they
       equal ``sensible_heat_flux`` and ``latent_heat_flux`` wherever no land-surface
       model supplies them. ``-999`` unless the two-stream solver is running with
       ``erf.radiation.seb_enable = true``.
   * - ``sens_flux``, ``laten_flux``
     - Selectable: fixed request names.
     - Value: legacy conservative surface flux outputs.
   * - ``surf_pres``
     - Selectable: fixed request name.
     - Value: the pressure at the ground, as WRF's ``PSFC``: the pressure of the
       lowest cell's state plus :math:`\rho g z_c`, with :math:`z_c` the height
       of that cell's centre above the ground -- the same surface pressure the
       radiation and the surface layer use.
   * - ``sea_level_pressure``
     - Selectable: fixed request name in dry and moist configurations.
     - Value: pressure reduced hydrostatically from the local terrain surface
       to ERF's physical ``z=0`` datum using the NMC/NGM lapse-rate method and
       Shuell correction.
   * - ``surface_diagnostic_source``
     - Selectable: fixed request name.
     - Value: categorical code ``0`` through ``6``; code ``0`` means no source, including when SurfaceLayer is absent.
   * - ``sensible_heat_flux``, ``latent_heat_flux``
     - Selectable: fixed request names.
     - Value: the selected conservative sources converted to energy-flux units; ``-999`` when unavailable.
   * - ``shoc_u_star``, ``shoc_Olen``, ``shoc_wthv_sfc``
     - Selectable: fixed request names.
     - Value: native SHOC diagnostics; ``-999`` when absent.
   * - ``precip_total_accum``, ``precip_frozen_accum``
     - Selectable: the active moisture scheme has any rain, snow, or graupel mass component.
     - Value: normalized liquid-water-equivalent accumulations; ``precip_frozen_accum`` is defined zero in a rain-only scheme.
   * - ``precip_rain_accum``
     - Selectable: the active moisture scheme has a rain mass component.
     - Value: normalized rain accumulation, or the derived rain value when a runtime total source is used without a direct rain source.
   * - ``precip_snow_accum``
     - Selectable: the active moisture scheme has a snow mass component.
     - Value: normalized snow accumulation.
   * - ``precip_graupel_accum``
     - Selectable: the active moisture scheme has a graupel mass component.
     - Value: normalized graupel accumulation.
   * - ``precip_hail_accum``
     - Selectable: not selectable for any current moisture scheme. The fixed descriptor reserves this public name for future support.
     - Value: no component is currently written. If a future scheme makes the field selectable, the descriptor's zero-fill policy applies when its runtime source is unavailable.
   * - ``integrated_qv``
     - Selectable: always.
     - Value: column-integrated water vapor; zero when moisture is disabled.
   * - ``integrated_qc``, ``integrated_qi``, ``integrated_qr``, ``integrated_qs``, ``integrated_qg``
     - Selectable: the active moisture model exposes the corresponding conserved mass component.
     - Value: the corresponding column water path. An unsupported species is skipped during request validation.
   * - Fixed land-surface provider names from ``t_sfc`` through ``noahmp_vegetation_fraction``
     - Selectable: the active land-surface inventory contains the exact name.
     - Value: finite provider values pass through when valid; absent, nonfinite, or sentinel values become ``-999``. Valid zero and signed values pass through when the field permits them. Raw fields never use MOST fallback.
   * - ``temperature_2m``
     - Selectable: a Noah-MP or SurfaceLayer pathway exists.
     - Value: a coherent native or MOST value, or ``-999`` when no complete source satisfies the request.
   * - ``water_vapor_mixing_ratio_2m``
     - Selectable: moisture is enabled and a Noah-MP or SurfaceLayer pathway exists.
     - Value: a coherent native or MOST mixing ratio, or ``-999`` when no complete source satisfies the request.
   * - ``near_surface_diagnostic_source``
     - Selectable: a Noah-MP or SurfaceLayer pathway exists, including a dry source-only request.
     - Value: categorical code ``0`` through ``6``; its ``AlwaysAvailable`` metadata policy records ``missing_value: null``.

.. _sec:Plotfile2DTopographicRadiation:

Topographic shading diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With RRTMGP on static fitted terrain, enabling ``erf.topo_shading`` or
``erf.slope_rad`` produces five diagnostics of the surface shortwave correction.
Both options default to false; see :ref:`inputs-topographic-radiation` for
configuration and :ref:`sec:TopographicRadiation` for the equations.

``topographic_shadow`` is 1 where intervening terrain blocks direct sunlight
and 0 otherwise. It remains 0 when shading is disabled, even if slope
correction is enabled. A slope facing away from the sun can have zero direct
flux without a cast shadow, so the mask alone does not describe all slope
effects.

``terrain_direct_factor`` multiplies the horizontal direct shortwave flux.
It is 0 in a cast shadow and can exceed 1 on a sun-facing slope.
``terrain_incidence`` is the nonnegative cosine of the angle between the sun
direction and the terrain normal. If RRTMGP's interval-averaged zenith cosine
is at or below ``1.e-4``, the mask is 0 and the direct factor is 1. If the
midpoint sun is at or below the same cutoff while the averaged cosine is above
it, incidence is 0 and the direct factor is 0; the mask is 1 when shading is
enabled and 0 when only slope correction is enabled. Thus mask 0 alone does
not imply daytime illumination.

``unshaded_sw_flux_dn`` and ``terrain_sw_flux_dn`` are the total downwelling
shortwave fluxes before and after correction, respectively. They include the
unchanged diffuse component, so corrected shortwave need not vanish in shadow.
Where Noah-MP receives this forcing, its ``sw_flux_dn`` corresponds to the
corrected flux.

These fields update with radiation at ``erf.rad_freq_in_steps`` and remain
fixed between radiation updates. More frequent plotfiles do not cause more
frequent shadow calculations. Before the first radiation update, all five
diagnostics use the missing value ``-999``. They are selectable only when
``erf.topo_shading`` or ``erf.slope_rad`` is enabled.

To visualize the mask and its effect, add these outputs to a case with terrain
correction enabled:

.. code-block:: text

   erf.plot2d_file_1 = plt2d
   erf.plot2d_int_1 = 1
   erf.plot2d_vars_1 = z_surf topographic_shadow terrain_direct_factor terrain_incidence unshaded_sw_flux_dn terrain_sw_flux_dn

In ParaView or VisIt, color by ``topographic_shadow`` and overlay ``z_surf``
contours to locate shadows relative to ridges and valleys. Compare the two
shortwave fields to see the forcing change. For a simple directional check,
the fixed low eastern sun in :ref:`inputs-topographic-radiation` should cast
shadows west of a hill; geographic west need not align with the negative grid
x direction on a rotated real-case grid. Use outputs after a radiation update
and exclude ``-999`` values from the color scale.

.. _sec:Plotfile2DSeaLevelPressure:

NMC/NGM sea-level pressure
^^^^^^^^^^^^^^^^^^^^^^^^^^

``sea_level_pressure`` is an opt-in surface-state diagnostic useful for
comparing model pressure with meteorological sea-level-pressure analyses. It
is an adaptation of the NMC/NGM sea-level-pressure reduction with the Shuell
correction from NOAA-EMC UPP ``sorc/ncep_post.fd/NGMSLP.f`` at immutable commit
`1ee7b175c340a566287a987870d925e4708fe19b
<https://github.com/NOAA-EMC/UPP/blob/1ee7b175c340a566287a987870d925e4708fe19b/sorc/ncep_post.fd/NGMSLP.f>`_.
ERF preserves the method and empirical correction but uses its native physical
constants, so bit-for-bit equality with UPP is not expected.

ERF stores water-vapor mixing ratio per unit dry-air mass, ``r_v``. The
reduction uses specific humidity ``q_v``:

.. math::

   q_v = \frac{r_v}{1+r_v},
   \qquad
   \epsilon_v = \frac{R_v}{R_d}-1,
   \qquad
   T_{v,m}=T_m\left(1+\epsilon_v q_v\right).

The original UPP routine uses the rounded coefficient ``0.608``; ERF computes
the equivalent coefficient from its shared ``R_d`` and ``R_v`` constants. ERF
diagnoses ``p_m`` and ``T_m`` from the lowest cell-centered conserved state
using its existing EOS.

Let ``z_m`` be the lowest cell-center physical height and ``z_s`` the lower
terrain-interface height. With ``Gamma = 0.0065 K m^-1``,

.. math::

   T_{v,s}=T_{v,m}+\Gamma(z_m-z_s),
   \qquad
   \overline{\tau}_{m,s}=\frac{R_d}{g}\frac{T_{v,m}+T_{v,s}}{2},
   \qquad
   p_s=p_m\exp\left(\frac{z_m-z_s}{\overline{\tau}_{m,s}}\right).

The ground pressure ``p_s`` is reconstructed because ERF does not store the
UPP ground-interface pressure in its cell-centered state. The uncorrected
sea-level extrapolation is ``T_{v,0}=T_{v,m}+Gamma z_m``. For ``T_c=290.66 K``,
the Shuell correction is

.. math::

   T_{v,0}^{*}=\begin{cases}
   T_c, & T_{v,0}>T_c\ \text{and}\ T_{v,s}\leq T_c,\\
   T_c-0.005(T_{v,s}-T_c)^2,
      & T_{v,0}>T_c\ \text{and}\ T_{v,s}>T_c,\\
   T_{v,0}, & \text{otherwise.}
   \end{cases}

The final reduction uses

.. math::

   \overline{\tau}_{s,0}=\frac{R_d}{g}\frac{T_{v,s}+T_{v,0}^{*}}{2},
   \qquad
   p_{\mathrm{MSL}}=\begin{cases}
   p_s, & |z_s|\leq 1\ \mathrm{m},\\
   p_s\exp\left(\frac{z_s}{\overline{\tau}_{s,0}}\right),
      & |z_s|>1\ \mathrm{m}.
   \end{cases}

The output units are pascals. It is pressure reduced to ERF's physical
``z=0`` datum and represents mean sea-level pressure only when that datum is
mean sea level. The diagnostic uses the lowest model level, assumes the fixed
lapse rate below it, applies no horizontal smoothing, does not calculate
1000-hPa height, is computed independently on each AMR level, and does not
alter model state. ``surf_pres`` remains a separate field, the pressure at the
ground.

.. code-block:: text

   erf.plot2d_vars_1 = surf_pres sea_level_pressure

.. _sec:Plotfile2DDynamicSoil:

Dynamic soil diagnostic families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dynamic soil fields are not part of the fixed catalog. ERF exposes only names
reported by the active land-surface provider inventory, in that inventory's
runtime order. ``<layer>`` is a one-based layer index.

.. BEGIN ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES

.. list-table::
   :header-rows: 1
   :widths: 28 42 20 35

   * - Name family
     - Meaning
     - Units
     - Metadata
   * - ``smois_<layer>``
     - Volumetric total soil moisture
     - ``m^3 m^-3``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.
   * - ``sh2o_<layer>``
     - Volumetric liquid soil water
     - ``m^3 m^-3``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.
   * - ``tslb_<layer>``
     - Soil temperature
     - ``K``
     - Category ``LandSurface``; policy ``FillMinus999WhenUnavailable``; invalid public values become ``-999``.

.. END ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES


.. _sec:Plotfile2DFluxes:

Surface flux diagnostics
------------------------

The 2D flux diagnostics report the conservative lower-boundary surface flux
that ERF used. ``sens_flux`` and ``laten_flux`` are the legacy conservative
scalar outputs. ``sensible_heat_flux`` and ``latent_heat_flux`` convert the
same selected conservative sources to ``W m^-2`` with ``Cp_d`` and ``L_v``.

For non-SHOC configurations and native SHOC host-diffusion mode, these outputs
use the host vertical surface flux arrays. In native SHOC ``state_update`` mode,
SHOC consumes those surface fluxes before the host diffusion path clears the
overlapping arrays. In that mode, the 2D flux diagnostics use SHOC's preserved
consumed-flux snapshots, component by component. If the corresponding host flux
field was unavailable before SHOC consumed it, ERF writes ``-999`` rather than a
zero SHOC snapshot.

The sign convention follows ERF's lower-boundary flux convention. No
source-specific Noah-MP, MOST, or SHOC conversion is applied in the 2D
diagnostic layer. Noah-MP LSM kinematic-to-conservative conversion happens
upstream in SurfaceLayer. Native SHOC converts the conservative host fluxes to
kinematic column inputs internally, then the diagnostic snapshot preserves the
consumed flux in conservative ERF output units.

In native SHOC ``state_update`` mode, SHOC is the transport owner. The
``surface_diagnostic_source`` field still describes the upstream surface-flux
source path used before SHOC consumes it.

.. _sec:Plotfile2DPrecipitation:

Surface precipitation accumulations
------------------------------------

ERF can write cumulative surface precipitation accumulations from the active
microphysics scheme in 2D plotfiles. These diagnostics report the liquid-water
equivalent mass that has reached the lower boundary since model start or the
most recent restart. Some schemes store explicit rain/snow/graupel species
accumulators, while others store a total accumulator plus frozen-species
subsets. ERF normalizes each available scheme-native source to ``kg/m^2``
before deriving the public 2D fields. Request acceptance follows the active
rain, snow, and graupel mass components: total and frozen require any one of
those components; each species field requires its own component; hail is not
selectable for current schemes. Runtime mapping may still derive rain from a
total source without a direct rain source. The public fields use ``kg/m^2``
because downstream land-surface forcing consumes precipitation as mass over
area, even though ``1 kg/m^2`` is numerically equal to ``1 mm`` of liquid-water
equivalent.

+-----------------------------+--------------------------------------------------+----------+
| Field                       | Meaning                                          | Units    |
+=============================+==================================================+==========+
| ``precip_total_accum``      | Accumulated surface precipitation, liquid-water  | kg/m^2   |
|                             | equivalent                                       |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_rain_accum``       | Accumulated surface rain precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_snow_accum``       | Accumulated surface snow precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_graupel_accum``    | Accumulated surface graupel precipitation,       | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_hail_accum``       | Accumulated surface hail precipitation,          | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+
| ``precip_frozen_accum``     | Accumulated frozen surface precipitation,        | kg/m^2   |
|                             | liquid-water equivalent                          |          |
+-----------------------------+--------------------------------------------------+----------+

``precip_rain_accum``, ``precip_snow_accum`` and ``precip_graupel_accum`` each
hold that one species, as their 3D counterparts ``rain_accum``, ``snow_accum``
and ``graupel_accum`` do. For ``Morrison``, ``WSM6`` and ``WDM6``,
``precip_total_accum`` is the scheme's precipitation total, the 3D
``precip_accum``, which also counts the cloud ice (and, for ``Morrison``, the
cloud water) that fell, so it can exceed the rain plus ``precip_frozen_accum``.
For the other schemes it is the rain plus the frozen species; for ``SAM`` that
leaves out the cloud ice that the 3D ``precip_accum`` counts.

``precip_total_accum`` is the normalized native total when one exists, or the
sum of normalized rain and frozen accumulations otherwise.
``precip_frozen_accum`` is the sum of normalized snow, graupel, and hail
accumulations. Species absent from the active runtime source contribute zero to
derived totals. ``precip_hail_accum`` is a fixed descriptor reserved for a
future scheme that exposes a distinct hail accumulator.

To diagnose the frozen fraction over a coupling interval, use accumulation
differences rather than a ratio of cumulative values:

.. math::

   f_{frozen} =
   \begin{cases}
   \dfrac{\Delta P_{frozen}}{\Delta P_{total}}, & \Delta P_{total} > 0 \\
   0, & \text{otherwise}
   \end{cases}

with

.. math::

   \Delta P_{total} = P_{total}(t_1) - P_{total}(t_0),

and

.. math::

   \Delta P_{frozen} = P_{frozen}(t_1) - P_{frozen}(t_0).

.. _sec:Plotfile2DWaterPaths:

Column water-path diagnostics
------------------------------

ERF can write scheme-aware 2D water-path diagnostics for prognostic condensed
water species. These diagnostics integrate the conserved ``rho*q`` species
component over the model column. They are available only when the active
microphysics scheme exposes the corresponding conserved mass component.

For a condensed species :math:`q_x`, ERF writes

.. math::

   W_x(x,y) = \int_{z_b}^{z_t} \rho(x,y,z)\,q_x(x,y,z)\,dz,

where :math:`W_x` has units of :math:`\mathrm{kg\,m^{-2}}`.

The discrete diagnostic uses the same metric convention as ``integrated_qv``:

.. math::

   W_x(i,j) = \Delta z \sum_k (\rho q_x)_{i,j,k}

for constant-:math:`\Delta z` meshes, and

.. math::

   W_x(i,j) = \Delta z \sum_k (\rho q_x)_{i,j,k} J_{i,j,k}

when ERF uses the column metric factor :math:`J`.

+-------------------------------+-----------------------------------------------+----------+
| **integrated_qc**             | Cloud liquid water path                       | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qi**             | Cloud ice water path                          | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qr**             | Rain water path                               | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qs**             | Snow water path                               | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+
| **integrated_qg**             | Graupel water path                            | kg/m^2   |
+-------------------------------+-----------------------------------------------+----------+

A name is selectable only if the active moisture model has that species as a
conserved mass component. Unsupported species are skipped during request
validation. Two-moment number concentrations are not water mass paths and are
not included.

The metadata sidecar records these fields as ``ColumnIntegral`` diagnostics
with ``FillZeroWhenUnavailable`` missing-value policy. Water-path diagnostics
use the built-in metadata fields and add no water-path-specific metadata keys.

.. _sec:Plotfile2DCloudOptics:

Estimated cloud optical depth and reflectance
---------------------------------------------

``cloud_optical_depth_sw`` and ``cloud_albedo_est`` are output-only,
dimensionless (units ``1``) column diagnostics. Either name can be selected
alone in either Plot2D stream or in station 2D output; both follow
``integrated_qg`` in canonical order regardless of input order. They are
selectable when **any** conserved cloud-liquid, cloud-ice, rain, snow, or
graupel mass component is present in the active moisture scheme. Dry schemes
and schemes without conserved hydrometeor mass exclude both names. There is
no dependency on RRTMGP, TwoStream, SHOC, Noah-MP, or a surface model.

For each present species :math:`s\in\{qc,qi,qr,qs,qg\}`:

.. math::

   W_s(i,j)=\int \max(0,(\rho q_s)_{i,j,z})\,dz,
   \qquad \tau_{sw}=\max\left(0,\sum_s k_{ext,s}W_s\right),
   \qquad A_{est}=\frac{\tau_{sw}}{\tau_{sw}+c}.

Here :math:`c` is the run-time ``reflectance_scale``; the default of 7.7 is
explained below.

The conserved component already contains density; it is not multiplied by
density again. Absent species contribute exactly zero, and number
concentrations are excluded. Negative mass values are clipped to zero before
integration. Zero condensate gives exactly zero optical depth and reflectance.
Nonfinite present mass (including a species with zero coefficient), a nonfinite
or nonpositive terrain metric, or arithmetic overflow gives ``-999`` for both
fields. The metadata policy is ``FillMinus999WhenUnavailable`` with
``missing_value: -999``. Valid reflectance is bounded to [0,1]; the missing
sentinel is outside this range. Coefficients must be finite and nonnegative;
``reflectance_scale`` must be finite and positive. An invalid setting aborts
the run at startup, naming the parameter, whether or not either field is
requested.

Constant-dz integration is :math:`\Delta z\sum_k \max(0,\rho q_s)`;
on nonconstant-dz/terrain meshes it is
:math:`\Delta z\sum_k J_{i,j,k}\max(0,\rho q_s)`, using the same
``detJ_cc`` convention as the water paths. Vertically stacked boxes are summed
before reflectance is calculated. **Both outputs are per-level column
proxies over the cells covered by that AMR level.** Missing vertical coverage
is not supplied from coarser levels, and there is no coarse-plus-fine composite
column sum. A fine level that covers only the bottom part of the atmosphere
retains its partial optical depth and the reflectance proxy calculated from
that partial optical depth; neither value describes the full atmosphere.
There is no additional mask for incomplete vertical coverage. This follows
the existing per-level water-path convention. The final planar fill still
gives ``-999`` where that level has no box reaching the ground. Thus an
incomplete bottom-reaching column can have a valid partial value, while a
column entirely aloft is missing in that level's slab. Stations sample the
proxy on the level chosen by the existing station stencil rules; they do not
combine level-columns into a composite atmospheric column.

.. warning::

   On a nested simulation, these are **per-level partial-column proxies**,
   not full-atmosphere or composite cloud albedos. A bottom-only fine level
   can report a much smaller optical depth than the coarse level if clouds
   lie above the fine level's top. Values integrate only hydrometeors covered
   by that level and are not combined across levels. Interpret them using each
   level's vertical coverage; an incomplete column must not be read as a
   full-column optical depth or albedo.

The five illustrative default mass-extinction coefficients, in
:math:`\mathrm{m^2\,kg^{-1}}`, are:

.. list-table::
   :header-rows: 1

   * - Species
     - Default
     - Representative radius and particle density
   * - Cloud liquid (qc)
     - 150.0
     - 10 micrometers; 1000 kg/m^3
   * - Cloud ice (qi)
     - 65.4
     - 25 micrometers; 917 kg/m^3
   * - Rain (qr)
     - 3.0
     - 0.5 mm; 1000 kg/m^3
   * - Snow (qs)
     - 30.0
     - 0.5 mm; 100 kg/m^3
   * - Graupel (qg)
     - 7.5
     - 0.5 mm; 400 kg/m^3

These defaults use an **equivalent-sphere approximation**:
:math:`k_{ext}=3Q_{ext}/(4\rho_{particle}r_{eff})`, with assumed
:math:`Q_{ext}\simeq2`. For a sphere this is its extinction cross section
:math:`Q_{ext}\pi r_{eff}^2` divided by its mass
:math:`4\pi\rho_{particle}r_{eff}^3/3`. The representative radii and densities
in the table are assumptions, not quantities diagnosed from the active
microphysics scheme. In the geometric-optics limit, where the size parameter
:math:`2\pi r_{eff}/\lambda` is large, the extinction paradox gives
:math:`Q_{ext}\rightarrow2`: diffraction contributes roughly another
geometric cross section to the extinction. This limit is best supported here
for spherical liquid cloud droplets at visible wavelengths. `Bohren and
Huffman (1983) <https://doi.org/10.1002/9783527618156>`_ discuss the extinction
paradox. `Grosvenor et al. (2018, Section 2.1, Eq. 2) <https://doi.org/10.1029/2017RG000593>`_
define extinction for a spherical droplet size distribution, use
:math:`Q_{ext}\simeq2` when droplet radius is large relative to wavelength,
and relate optical depth to liquid water and effective radius in Section 2.2,
Eq. 6.

The representative radii and densities are fixed, tunable bulk proxies,
**not scheme-calibrated or retrieval-grade optics**. Ice, snow, and graupel
are often non-spherical; their equivalent radii and densities are particularly
crude proxies. Real mass-size relations, scattering properties, and particle
size distributions may differ greatly. Sensitivity test these assumptions
before quantitative research use.

The inputs configure the five extinction coefficients in
:math:`\mathrm{m^2\,kg^{-1}}`, not radius or particle density directly. Users
wanting different assumed effective radii or densities must calculate the
corresponding coefficients themselves and supply ``kext_qc``, ``kext_qi``,
``kext_qr``, ``kext_qs``, and/or ``kext_qg`` overrides.

The common two-stream approximation for a conservatively scattering,
horizontally homogeneous cloud with asymmetry factor :math:`g` gives

.. math::

   A_{est}=\frac{\tau}{\tau+c},\qquad
   c=\frac{2}{\sqrt{3}(1-g)}\simeq7.7\quad\text{for }g=0.85.

`Lacis and Hansen (1974) <https://www.giss.nasa.gov/pubs/abs/la05000n.html>`_
and `Lu et al. (2008, Section 4.3) <https://doi.org/10.1029/2007JD009354>`_
give the corresponding nonabsorbing-cloud expression
:math:`R=\tau/(\tau+7.7)` for :math:`g=0.85`. This published approximation
motivates the default ``reflectance_scale = 7.7``; it does **not** make that
value a universal or scheme-calibrated albedo parameter. Its assumptions,
including a fixed asymmetry factor, are particularly limited for ice and
mixed-phase clouds, whose scattering properties vary with particle habit,
size, and wavelength. ``reflectance_scale`` remains a positive runtime
override, and sets the optical depth at which this proxy equals one half.
The implemented proxy is independent of solar zenith angle and does not solve
full radiative transfer, as described in the warning below.

.. warning::

   ``cloud_albedo_est`` is an exploratory hydrometeor-only
   reference-illumination reflectance proxy. It includes cloud liquid, cloud
   ice, rain, snow, and graupel; excludes surface, gas, and aerosol
   contributions; and does not solve two-stream radiative transfer. It
   omits solar zenith angle, gas/aerosol scattering and absorption,
   single-scattering albedo, phase-function asymmetry, overlap, spectral
   dependence, surface reflection, and 3-D cloud effects. It is neither
   physical broadband all-sky albedo nor top-of-atmosphere reflectance, nor
   RRTMGP top-of-model flux-derived all-sky reflectance. Existing ``albedo``
   is surface albedo and remains a separate diagnostic. These new fields do
   not affect prognostic physics, radiative heating, or fluxes. The proxy is
   independent of solar zenith angle. No solar flux enters its definition:
   it remains defined, and can remain nonzero, at night for a cloudy column.
   Nighttime values must not be interpreted as observed instantaneous albedo.

Example (water-path outputs are optional):

.. code-block:: text

   erf.plot2d_file_1 = plt2d
   erf.plot2d_per_1 = 600.0
   erf.plot2d_vars_1 = integrated_qc integrated_qi integrated_qr integrated_qs integrated_qg cloud_optical_depth_sw cloud_albedo_est
   erf.plot2d.cloud_optics.kext_qc = 150.0
   erf.plot2d.cloud_optics.kext_qi = 65.4
   erf.plot2d.cloud_optics.kext_qr = 3.0
   erf.plot2d.cloud_optics.kext_qs = 30.0
   erf.plot2d.cloud_optics.kext_qg = 7.5
   erf.plot2d.cloud_optics.reflectance_scale = 7.7

Reproducing the optical settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``2DMetadata.json`` keeps the generic diagnostic-descriptor schema; it does
not record the per-run cloud-optics coefficients. Archive the input deck,
``job_info``, ERF commit/build configuration, and complete standard output
with the plotfiles. On first requested use, ERF prints one ``Plot2D cloud
optics:`` line containing **all five resolved coefficients and the resolved
reflectance scale**, including defaults omitted from the input deck. Values
are printed with enough digits to recover the build's ``Real`` values.
For example, capture and retain this record alongside the plotfiles:

.. code-block:: bash

   ./erf_exec inputs > run.log 2>&1
   sed -n '/^Plot2D cloud optics:/p' run.log > cloud_optics_resolved_settings.txt

Keep both files, rather than relying on ``job_info`` alone: supplied overrides
may be present there while unprovided defaults are absent. The writer integration
test archives this same resolved-settings record. Logging and these settings do not
change the diagnostic values or the generic metadata format.

References
^^^^^^^^^^

* Lacis, A. A., and J. E. Hansen (1974). `A parameterization for the
  absorption of solar radiation in the Earth's atmosphere
  <https://www.giss.nasa.gov/pubs/abs/la05000n.html>`_. *Journal of the
  Atmospheric Sciences*, 31, 118–133.
* Lu, M.-L., et al. (2008). `Aerosol-cloud relationships in continental shallow
  cumulus <https://doi.org/10.1029/2007JD009354>`_. *Journal of Geophysical
  Research: Atmospheres*, 113, D15201.
* Bohren, C. F., and D. R. Huffman (1983). `Absorption and Scattering of Light
  by Small Particles <https://doi.org/10.1002/9783527618156>`_. Wiley.
* Grosvenor, D. P., et al. (2018). `Remote sensing of droplet number
  concentration in warm clouds: A review of the current state of knowledge and
  perspectives <https://doi.org/10.1029/2017RG000593>`_. *Reviews of
  Geophysics*, 56, 409–453.

.. _sec:Plotfile2DSourceCodes:

Surface diagnostic source codes
--------------------------------

``surface_diagnostic_source`` is a cell-centered categorical diagnostic. It
reports the source path used by the SurfaceLayer scalar diagnostic path. It
does not report fractional land-cover contributions. The same code table
defines ``near_surface_diagnostic_source``; that field reports the request-
aware source selected for the unified 2-m output bundle.

If an input dataset contains fractional land information, this diagnostic still
reports the categorical source used by ERF's active SurfaceLayer scalar flux
path.

The fields do not describe fractional cover or complete staggered stress-face
provenance. Consumers must compare the numeric code with this table rather than
infer provenance from ``landmask``.

.. list-table::
   :header-rows: 1
   :widths: 10 32 58

   * - Value
     - Token
     - Meaning
   * - 0
     - ``missing``
     - No source was selected.
   * - 1
     - ``surface_layer_land``
     - SurfaceLayer supplied the land value.
   * - 2
     - ``lsm_land``
     - The land-surface model supplied the land value.
   * - 3
     - ``surface_layer_fallback``
     - A land-surface path existed, but its required result was invalid or incomplete; SurfaceLayer supplied the value.
   * - 4
     - ``surface_layer_sea``
     - SurfaceLayer supplied the value over water.
   * - 5
     - ``custom``
     - The custom surface pathway supplied the SurfaceLayer state.
   * - 6
     - ``rico``
     - The RICO pathway supplied the SurfaceLayer state.


.. _sec:Plotfile2DMetadata:

AMReX metadata and NetCDF differences
--------------------------------------

Native AMReX 2D plotfiles write a JSON metadata sidecar named
``2DMetadata.json`` in the plotfile directory. The sidecar lists selected
outputs in component order. Fixed outputs use descriptor metadata. Sampled-
level outputs add source-field and vertical-coordinate metadata.

The metadata policy maps to JSON as follows:

* ``AlwaysAvailable`` records ``missing_value: null``.
* ``FillZeroWhenUnavailable`` records ``missing_value: 0``.
* ``FillMinus999WhenUnavailable`` records ``missing_value: -999``.

The sidecar uses the same 2D diagnostic catalog that defines the built-in
plotfile variables. Sampled-level outputs carry their own metadata record. The
sidecar does not record the runtime source chosen at each cell. Source choice is
carried by categorical output fields. It does not change field values.

NetCDF 2D output uses the same variable names but does not write this JSON
sidecar or sampled-level metadata attributes. The sidecar format version is
``2``.

Example built-in metadata:

.. code-block:: json

   {
     "format_version": 2,
     "kind": "ERF 2D plotfile metadata",
     "n_variables": 2,
     "variables": [
       {
         "component_index": 0,
         "name": "z_surf",
         "long_name": "Surface elevation",
         "units": "m",
         "category": "Geometry",
         "missing_policy": "AlwaysAvailable",
         "missing_value": null
       },
       {
         "component_index": 1,
         "name": "latent_heat_flux",
         "long_name": "Surface latent heat flux",
         "units": "W m^-2",
         "category": "SurfaceFlux",
         "missing_policy": "FillMinus999WhenUnavailable",
         "missing_value": -999
       }
     ]
   }


Sampled-level entries add ``source_field`` and ``vertical_coordinate`` records.
See :ref:`sec:Plotfile2DSampledMetadata` for their exact form.
