.. _SLM:

Simple Land-surface Model (SLM)
===============================

The Simple Land-surface Model (SLM) supplies lower-boundary sensible-heat,
latent-heat, and momentum fluxes for land cells.  Select it with
``erf.land_surface_model = "SLM"``.  SLM uses the ``slm.`` input prefix.

SLM can initialize a horizontally uniform land surface from the input file,
read land-surface fields from ``WRFInput``, or use reference sounding, flux,
and skin-temperature files for testing.  The soil layers are ordered from the
surface downward.  ``slm.soil_dz`` gives the thickness of each layer in m and
must contain exactly ``slm.nsoil`` values.

Build and coupling requirements
-------------------------------

SLM is part of the standard ERF build.  It is coupled through the ERF surface
layer and is normally used with a land lower boundary.  Terrain-fitted
coordinates can be used; SLM receives the physical-height field from ERF when
terrain is active.

Basic configuration
-------------------

The minimum uniform-initialization configuration supplies the soil-layer
thicknesses and initial soil/vegetation properties.  For example:

.. code-block:: text

   erf.land_surface_model = "SLM"
   slm.nsoil = 7
   slm.soil_dz = 0.02 0.04 0.08 0.16 0.32 0.64 1.28
   slm.landtype0 = 1
   slm.LAI0 = 3.0
   slm.clay0 = 34.0
   slm.sand0 = 10.0
   slm.sw0 = 0.50
   slm.st0 = 300.0

Core soil and surface options
-----------------------------

+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| Parameter                        | Definition                                               | Acceptable Values    | Default          |
+==================================+==========================================================+======================+==================+
| **slm.nsoil**                    | number of soil layers                                    | Integer >= 1         | 7                |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soil_dz**                  | thickness of each soil layer [m], from the surface       | Real values; exactly | must be set      |
|                                  | downward; the number of values must equal ``nsoil``      | ``nsoil`` values     |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.SLM_use_inputs**           | use SLM input/reference data files instead of the        | Boolean              | false            |
|                                  | normal uniform initialization; also enables the          |                      |                  |
|                                  | reference-file path when not using ``WRFInput``          |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.landtype0**                | initial land-use category applied over the domain        | Integer land type    | 16               |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.LAI0**                     | initial leaf-area index                                  | Real >= 0            | 0.0              |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.clay0**                    | clay fraction/content for each soil layer [%]; not read  | One Real or          | required unless  |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.sand0**                    | sand fraction/content for each soil layer [%]; not read  | One Real or          | required unless  |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.sw0**                      | initial soil wetness fraction for each layer; not read   | One Real or          | must be set      |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.st0**                      | initial soil temperature [K] for each layer; not read    | One Real or          | must be set      |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.relax_hgt**                | depth-dependent soil nudging weights; read when either   | One Real or          | required when    |
|                                  | soil nudging option is enabled                           | ``nsoil`` Reals      | nudging is on    |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soiltnudging**             | nudge soil temperature toward reference values           | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soilwnudging**             | nudge soil wetness toward reference values               | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.tausoil**                  | soil nudging time scale [s]                              | Real > 0             | 86400.0          |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.tabs_s**                   | prescribed/initial surface temperature [K]               | Real                 | 0.0              |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.t00**                      | constant temperature offset used in the surface          | Real [K]             | 300.0            |
|                                  | temperature field                                        |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.z0_soil**                  | bare-soil roughness length [m]                           | Real > 0             | 0.01             |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.mws_mx0**                  | maximum puddle-water storage [mm]                        | Real >= 0            | 50.0             |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.Rc_max**                   | maximum stomatal resistance [s/m]                        | Real > 0             | 5000.0           |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.T_opt**                    | optimum temperature for transpiration [K]                | Real                 | 298.0            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.zref**                     | atmospheric reference height above the surface [m]       | Real > 0             | 0.5              |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+

Parameter tables and external forcing
---------------------------------------

+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| Parameter                        | Definition                                               | Acceptable Values    | Default          |
+==================================+==========================================================+======================+==================+
| **slm.use_parameter_file**       | initialize soil and vegetation parameters from the       | Boolean              | false            |
|                                  | Noah-MP-format parameter file                            |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.radiation_scheme**         | radiation scheme used by SLM                             | ``SLM``, ``NoahMP``  | ``NoahMP``       |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.parameter_file**           | parameter-table filename used when                       | String               | NoahmpTable.TBL  |
|                                  | ``use_parameter_file`` is true                           |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.veg_dataset**              | vegetation dataset within the parameter file             | ``usgs``, ``modis``  | ``modis``        |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soil_dataset**             | soil dataset within the parameter file                   | ``stas``,            | ``stas``         |
|                                  |                                                          | ``stas_ruc``         |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.interpolate_lai**          | interpolate monthly LAI and SAI tables between months    | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.lai**                      | monthly LAI table, one row per month and one column per  | Table of Reals;      | required when    |
|                                  | land type; requires 12 rows                              | 12 rows              | interpolation    |
|                                  |                                                          |                      | is on            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.sai**                      | monthly SAI table, one row per month and one column per  | Table of Reals;      | required when    |
|                                  | land type; must match the LAI table shape                | 12 rows              | interpolation    |
|                                  |                                                          |                      | is on            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.use_param_tbl**            | use WRF/Noah-style vegetation parameter table            | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.vegparam**                 | vegetation parameter table used when                     | Table of Reals       | required when    |
|                                  | ``use_param_tbl`` is true                                |                      | use_param_tbl    |
|                                  |                                                          |                      | is true          |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.rad_input_file**           | NetCDF radiation-forcing file containing SW/LW and       | String               | empty            |
|                                  | cosine-zenith-angle fields; requires NetCDF support      |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+

``use_parameter_file`` and ``use_param_tbl`` are mutually exclusive.  The
parameter-file datasets must be ``usgs`` or ``modis`` for vegetation and
``stas`` or ``stas_ruc`` for soil.  The radiation file is expected to contain
the variables ``SWVIS``, ``SWNIR``, ``SWVISD``, ``SWNIRD``, ``COSZRS``, and
``LWDS`` with dimensions matching the horizontal SLM grid.

Reference-data testing options
------------------------------

These options support the SLM reference-data testing path.  They are read
when ``slm.SLM_use_inputs = true``.  The three files contain, respectively,
the sounding, surface flux, and surface-temperature reference data.

+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| Parameter                        | Definition                                               | Acceptable Values    | Default          |
+==================================+==========================================================+======================+==================+
| **slm.SLM_num_ref_inputs**       | number of reference input records/files                  | Integer >= 1         | 1                |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.SLM_ref_sounding_file**    | reference sounding file; columns are time, pressure,     | String               | must be set      |
|                                  | temperature, humidity, u velocity, and v velocity        |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.SLM_ref_flux_file**        | reference flux file; columns are time, SW down, LW down, | String               | must be set      |
|                                  | SW up, and LW up                                         |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.SLM_ref_sst_file**         | reference skin-temperature file; columns are time, SST,  | String               | must be set      |
|                                  | and precipitation                                        |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.start_time**               | starting time selected from the reference data           | Real                 | -1.0             |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.time_unit**                | seconds represented by one reference-data time unit      | Real > 0             | 1.0              |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+

When ``slm.SLM_use_inputs = false``, ``slm.zref`` is recomputed from the
canopy height and the lowest atmospheric cell unless the run is initialized
from ``WRFInput``.  When it is true, ``slm.zref`` is used directly.

WRFInput initialization
------------------------

With ``erf.init_type = WRFInput``, SLM reads the mapped land-surface fields
from the WRF input data, including soil thickness and temperature/moisture,
LAI, vegetation and soil type, skin temperature, and vegetation fractions.
The uniform ``clay0``, ``sand0``, ``sw0``, and ``st0`` values are therefore not
used in this mode; clay and sand are derived from the WRF soil type.

The active WRF-to-SLM field mapping is:

+----------------+----------------------------------------------------------+
| WRF name       | SLM variable and description                             |
+================+==========================================================+
| **DZS**        | ``soil_thickness``: thickness of each soil layer         |
+----------------+----------------------------------------------------------+
| **ZS**         | ``node_z``: depth/elevation of soil-layer nodes          |
+----------------+----------------------------------------------------------+
| **TSLB**       | ``tsoil``: soil temperature                              |
+----------------+----------------------------------------------------------+
| **SMOIS**      | ``wsoil``: soil wetness                                  |
+----------------+----------------------------------------------------------+
| **LAI**        | ``lai``: leaf-area index                                 |
+----------------+----------------------------------------------------------+
| **IVGTYP**     | ``vegtype``: vegetation type                             |
+----------------+----------------------------------------------------------+
| **ISLTYP**     | ``soiltype``: soil type                                  |
+----------------+----------------------------------------------------------+
| **TSK**        | ``tsurf``: surface/skin temperature                      |
+----------------+----------------------------------------------------------+
| **VEGFRA**     | ``veg_frac``: vegetation fraction                        |
+----------------+----------------------------------------------------------+
| **SHDMIN**     | ``veg_frac_min``: yearly minimum vegetation fraction     |
+----------------+----------------------------------------------------------+
| **SHDMAX**     | ``veg_frac_max``: yearly maximum vegetation fraction     |
+----------------+----------------------------------------------------------+

The WRFInput values are post-processed on the first SLM advance.  The
topmost soil layer—the layer immediately below the atmospheric surface—is
indexed internally by ``d_khi_lsm`` and supplies the soil-side surface values.
For each horizontal cell, SLM performs the following operations:

* WRF soil moisture is divided by the SLM soil porosity.  Thus the SLM ``wsoil``
  value is a normalized wetness fraction, with saturation represented relative to
  ``poro_soil``.

For WRF ``ISLTYP`` values 1--16, SLM assigns the following sand and clay
contents to every soil layer in the cell.

+----------+-----------------------+----------+----------+
| ISLTYP   | Soil type             | Sand [%] | Clay [%] |
+==========+=======================+==========+==========+
| 1        | Sand                  | 92.0     | 3.0      |
+----------+-----------------------+----------+----------+
| 2        | Loamy sand            | 82.0     | 6.0      |
+----------+-----------------------+----------+----------+
| 3        | Sandy loam            | 65.0     | 10.0     |
+----------+-----------------------+----------+----------+
| 4        | Silt loam             | 20.0     | 15.0     |
+----------+-----------------------+----------+----------+
| 5        | Silt                  | 8.0      | 12.0     |
+----------+-----------------------+----------+----------+
| 6        | Loam                  | 40.0     | 20.0     |
+----------+-----------------------+----------+----------+
| 7        | Sandy clay loam       | 60.0     | 30.0     |
+----------+-----------------------+----------+----------+
| 8        | Clay loam             | 32.0     | 34.0     |
+----------+-----------------------+----------+----------+
| 9        | Silty clay loam       | 20.0     | 40.0     |
+----------+-----------------------+----------+----------+
| 10       | Sandy clay            | 52.0     | 42.0     |
+----------+-----------------------+----------+----------+
| 11       | Silty clay            | 6.0      | 47.0     |
+----------+-----------------------+----------+----------+
| 12       | Clay                  | 20.0     | 60.0     |
+----------+-----------------------+----------+----------+
| 13       | Organic material      | 0.1      | 0.1      |
+----------+-----------------------+----------+----------+
| 14       | Water                 | 0.1      | 0.1      |
+----------+-----------------------+----------+----------+
| 15       | Bedrock               | 0.1      | 0.1      |
+----------+-----------------------+----------+----------+
| 16       | Other/urban           | 0.1      | 0.1      |
+----------+-----------------------+----------+----------+

The category constants and their WRFInput treatment are:

+----------------------+----------------------+----------------------------------------------------------+
| Constant             | Value                | Treatment                                                |
+======================+======================+==========================================================+
| ``isnature``         | 14                   | Urban-category cells are remapped to this nature         |
|                      |                      | category for the non-building SLM treatment.             |
+----------------------+----------------------+----------------------------------------------------------+
| ``isurban``          | 13                   | The urban constant is retained; its assignment is        |
|                      |                      | currently commented out.                                 |
+----------------------+----------------------+----------------------------------------------------------+
| ``islake``           | 21                   | Water category with WRF land type 17.                    |
+----------------------+----------------------+----------------------------------------------------------+

The category-specific defaults applied by this branch are:

+----------------------+----------------------------------------------------------+
| Category             | Values assigned by SLM                                   |
+======================+==========================================================+
| Water (land type 17  | ``landmask = 0``; ``landtype = 0``;                      |
| or ``islake`` = 21)  | ``vegetype = 0``; ``vege_YES = 0``; soil temperature     |
|                      | ``= 273.16 K``; surface emissivity ``= 0.98``; and all   |
|                      | visible/NIR direct and diffuse albedos ``= 0.06``.       |
+----------------------+----------------------------------------------------------+
| Urban (land type     | ``landtype`` is changed to ``isnature`` (14); when       |
| greater than 16)     | ``use_param_tbl`` is enabled, the nature-category        |
|                      | vegetation fraction from ``vegparam`` is used.           |
+----------------------+----------------------------------------------------------+
| Other land           | WRF land type and LAI are retained; the cell remains     |
|                      | eligible for land-mask and porosity-based soil-moisture  |
|                      | processing.                                              |
+----------------------+----------------------------------------------------------+

The local surface fields are initially allocated with vegetation fraction
one.  WRF vegetation fractions replace that default before category handling;
water-category cells subsequently receive the zero vegetation values listed
above.  The configured ``slm.landtype0`` and ``slm.LAI0`` defaults are used by
the non-WRFInput path and are overwritten by WRF fields in this path.

Diagnostics and coupling
------------------------

SLM exports surface temperature, albedos, emissivity, and surface fluxes to
the ERF surface-layer and radiation couplings.  Its internal soil and canopy
fields are available through SLM plotfile/checkpoint output.  Reference-file
testing additionally updates the atmospheric reference state and surface
forcing from the files listed above.
