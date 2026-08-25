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
     - Surface temperature from the surface layer
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
   * - ``OLR``
     - ``Radiation``
     - ``W/m^2``
     - ``FillMinus999WhenUnavailable``
     - Outgoing longwave radiation at the model top
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
     - Surface pressure
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
     - ``LandSurface``
     - ``K``
     - ``FillMinus999WhenUnavailable``
     - Physical air temperature 2 m above the local surface
   * - ``water_vapor_mixing_ratio_2m``
     - ``LandSurface``
     - ``kg kg^-1 dry air``
     - ``FillMinus999WhenUnavailable``
     - Water-vapor mixing ratio 2 m above the local surface per unit dry-air mass
   * - ``near_surface_diagnostic_source``
     - ``LandSurface``
     - ``1``
     - ``AlwaysAvailable``
     - Source code for the unified near-surface diagnostic bundle
.. END ERF BUILT-IN 2D DIAGNOSTIC CATALOG

.. note::

   For the ``landmask`` variable, land is ``1`` and sea is ``0``. Buildings are
   ``2`` when using ``ImmersedForcing``.

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
   * - ``pblh``
     - Selectable: fixed request name.
     - Native SHOC ``pblh`` is reported in metres above local ground (AGL). Value: native SHOC PBL height when native SHOC diagnostics are present; otherwise SurfaceLayer; ``-999`` if neither exists.
   * - ``OLR``
     - Selectable: fixed request name.
     - Value: radiation output; ``-999`` when the radiation source is absent.
   * - ``sens_flux``, ``laten_flux``
     - Selectable: fixed request names.
     - Value: legacy conservative surface flux outputs.
   * - ``surf_pres``
     - Selectable: fixed request name.
     - Value: pressure computed from the lowest atmospheric state.
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
alter model state. ``surf_pres`` remains a separate lowest-cell-center field.

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
