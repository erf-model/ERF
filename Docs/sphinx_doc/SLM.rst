.. _SLM:

Simple Land-surface Model (SLM)
===============================

The original formulation and evaluation of the Simple Land Model (SLM) are
documented in Lee, J. M., and M. Khairoutdinov (2015), *A Simplified Land Model
(SLM) for use in cloud-resolving models: Formulation and evaluation*, Journal
of Advances in Modeling Earth Systems, 7, 1368--1392,
`doi:10.1002/2014MS000419 <https://doi.org/10.1002/2014MS000419>`_. Readers
should consult that publication for the complete original model formulation.

This page documents only the recent changes made to improve SLM and make its
behavior more comparable to other contemporary land-surface models. These
changes include selected formulations adapted from NoahMP, updates inherited
from gSAM-SLM, and SLM-specific corrections to the surface energy and water
budgets. SLM remains a simplified land model; it is not a complete
implementation of NoahMP.

The Simple Land-surface Model (SLM) supplies lower-boundary sensible-heat,
latent-heat, and momentum fluxes for land cells.  Select it with
``erf.land_surface_model = "SLM"``.  SLM uses the ``slm.`` input prefix.

SLM can initialize a horizontally uniform land surface from the input file or
read land-surface fields from ``WRFInput``.  The soil layers are ordered from
the surface downward.  ``slm.soil_dz`` gives the thickness of each layer in m
and must contain exactly ``slm.nsoil`` values.

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
| **slm.nsoil**                    | number of soil layers                                    | Integer >= 2         | 7                |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soil_dz**                  | SLM layer layout: thickness of each soil layer [m],      | Real > 0; exactly    | must be set      |
|                                  | from the surface downward; the number of values must     | ``nsoil`` values     |                  |
|                                  | equal ``nsoil``. For WRFInput, per-cell thickness comes  |                      |                  |
|                                  | from WRF ``DZS``.                                        |                      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.landtype0**                | initial land-use category applied over the domain        | Integer land type    | 16               |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.LAI0**                     | initial leaf-area index                                  | Real >= 0            | 0.0              |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.clay0**                    | clay fraction/content for each soil layer [%]; not read  | One Real in [0,100]  | required unless  |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.sand0**                    | sand fraction/content for each soil layer [%]; not read  | One Real in [0,100]  | required unless  |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.sw0**                      | initial soil wetness fraction for each layer; not read   | One Real in [0,1]    | must be set      |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.st0**                      | initial soil temperature [K] for each layer; not read    | One Real or          | must be set      |
|                                  | when initializing from ``WRFInput``                      | ``nsoil`` Reals      |                  |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.relax_hgt**                | depth-dependent soil nudging weights; read when either   | One Real in [0,1] or | required when    |
|                                  | soil nudging option is enabled                           | ``nsoil`` Reals      | nudging is on    |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soiltnudging**             | nudge soil temperature toward reference values           | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.soilwnudging**             | nudge soil wetness toward reference values               | Boolean              | false            |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.tausoil**                  | soil nudging time scale [s]                              | Real > 0             | 86400.0          |
+----------------------------------+----------------------------------------------------------+----------------------+------------------+
| **slm.tabs_s**                   | prescribed/fallback surface temperature [K]; WRFInput    | Real                 | 0.0              |
|                                  | land-cell surface temperature comes from ``TSK``         |                      |                  |
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
| **slm.use_parameter_file**       | use values from the parameter file; otherwise use        | Boolean              | false            |
|                                  | built-in MODIS/STAS radiation values                     |                      |                  |
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

When ``use_parameter_file`` is false, SLM uses built-in typical MODIS/STAS
values for vegetation optical properties, canopy parameters, and soil albedos.
When it is true, values from the parameter file are used directly.

``use_parameter_file`` and ``use_param_tbl`` are mutually exclusive.  When
``use_param_tbl`` is true, SLM reads the ``slm.vegparam`` table for
WRF/Noah-style LAI and vegetation-fraction handling.  The
parameter-file datasets must be ``usgs`` or ``modis`` for vegetation and
``stas`` or ``stas_ruc`` for soil.  The radiation file is expected to contain
the variables ``SWVIS``, ``SWNIR``, ``SWVISD``, ``SWNIRD``, ``COSZRS``, and
``LWDS`` with dimensions matching the horizontal SLM grid.

WRFInput initialization
------------------------

With ``erf.init_type = WRFInput``, SLM reads the mapped land-surface fields
from the WRF input data, including soil thickness and temperature/moisture,
LAI, vegetation and soil type, skin temperature, and vegetation fractions.
The uniform ``clay0``, ``sand0``, ``sw0``, and ``st0`` values are therefore not
used in this mode; clay and sand are derived from the WRF soil type.

``slm.soil_dz`` is still required to define the SLM layer layout and must contain
exactly ``slm.nsoil`` values.  The per-cell soil thickness used by SLM is read
from WRF ``DZS``.

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

Recent physics changes
----------------------

1. Scope and relationship to the original SLM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The updated implementation retains the multilayer soil-temperature and
soil-water framework of SLM while replacing or extending selected surface and
vegetation processes. The principal changes are:

* NoahMP-format soil, vegetation, and radiation parameters;
* NoahMP-derived canopy and ground radiation;
* NoahMP-derived undercanopy, leaf-boundary, stomatal, and soil-vapor
  resistances;
* revised root-zone water stress and layer-specific transpiration extraction;
* a prognostic canopy energy and water balance;
* a massless ground-skin temperature and an explicit ground surface-energy
  balance; and
* revised interception, pond storage, infiltration, evaporation, and
  saturation corrections.

Not all of these changes are direct ports from NoahMP. Canopy-temperature
subcycling and rain cooling follow gSAM-SLM, while several water-budget and
initialization changes are specific corrections in the ERF SLM implementation.

2. Summary of the updated time step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ordering of the surface calculations is important because each process
uses state produced by an earlier process. On each SLM advance, the updated
implementation:

#. applies optional soil-temperature and soil-water nudging;
#. updates time-varying LAI and SAI;
#. computes canopy and ground radiation;
#. intercepts precipitation and updates canopy water storage;
#. computes surface-layer transfer coefficients;
#. computes undercanopy, leaf-boundary, stomatal, and soil resistances;
#. subcycles the canopy energy and moisture balance;
#. solves the massless ground-skin energy balance;
#. advances soil water, including evaporation and transpiration sinks; and
#. advances soil temperature using the skin-to-soil conductive heat flux.

The canopy and ground fluxes are subsequently combined to form the surface
fluxes supplied to the atmosphere. This ordering is an SLM implementation
choice and should not be interpreted as the time-integration sequence of the
full NoahMP model.

3. NoahMP parameter tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameter-file options and accepted datasets are listed above. The recent
physics changes use additional values from the selected NoahMP blocks.

The selected NoahMP soil table supplies the following SLM properties:

.. list-table::
   :widths: 18 47 35
   :header-rows: 1

   * - NoahMP name
     - SLM use
     - Conversion
   * - ``MAXSMC``
     - Soil porosity or saturated moisture content
     - None
   * - ``REFSMC``
     - Volumetric moisture at field capacity
     - None
   * - ``WLTSMC``
     - Volumetric moisture at the wilting point
     - None
   * - ``SATPSI``
     - Saturated matric potential
     - Positive m to negative mm
   * - ``BB``
     - Clapp--Hornberger hydraulic exponent
     - None
   * - ``SATDK``
     - Saturated hydraulic conductivity
     - m/s to mm/s

The vegetation table can replace SLM defaults for minimum stomatal resistance
(``RS``), radiation stress (``RGL``), leaf orientation (``XL``), vapor-pressure
deficit sensitivity (``HS``), rooted-layer count (``NROOT``), canopy height
(``HVT``), roughness length (``Z0MVT``), canopy biomass heat capacity
(``CBIOM``), and characteristic leaf dimension (``DLEAF``). The radiation
calculation additionally uses the leaf and stem visible/NIR reflectances and
transmittances, crown radius ``RC``, canopy-bottom height ``HVB``, tree density
``DEN``, and canopy wind-extinction parameter ``CWPVT``.

After table values replace SLM defaults, all dependent soil and vegetation
properties are recomputed. This prevents derived conductivities, moisture
thresholds, emissivities, interception capacities, or radiation coefficients
from retaining values calculated from the old defaults.

The reader validates category counts, required parameter blocks, duplicate or
missing values, optical properties outside [0,1], invalid canopy geometry, and
land or soil categories outside the selected table. ``NROOT`` must be an
integer between zero and ``slm.nsoil``.

4. NoahMP-derived canopy and ground radiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The updated SLM radiation follows the NoahMP canopy-radiation formulation. It
uses separate direct and diffuse fluxes in visible and near-infrared bands:

* direct visible and near-infrared shortwave;
* diffuse visible and near-infrared shortwave;
* atmospheric downwelling longwave; and
* cosine of the solar zenith angle.

These quantities normally come from the atmospheric radiation model. They may
also be supplied from ``slm.rad_input_file``, a NetCDF file containing
``SWVIS``, ``SWNIR``, ``SWVISD``, ``SWNIRD``, ``COSZRS``, and ``LWDS`` on the
SLM horizontal grid.

The canopy optical calculation retains LAI and SAI as separate quantities.
Leaf and stem optical properties are weighted by their fractions of the total
vegetation area index

.. math::

   VAI = LAI + SAI.

A modified two-stream calculation partitions incoming solar energy into
reflected radiation, radiation absorbed by vegetation, and radiation absorbed
by the ground. It accounts for canopy orientation, direct and diffuse beams,
vegetation fraction, and canopy gaps. The shortwave terms obey the budget

.. math::

   SW_{in} = SW_{reflected} + SW_{canopy,abs} + SW_{ground,abs}.

Bare-soil cells bypass the canopy two-stream calculation and use the
moisture-dependent ground albedo directly.

The longwave calculation treats the canopy and ground as separate emitting
surfaces. With canopy emissivity :math:`\epsilon_v`, ground emissivity
:math:`\epsilon_g`, canopy temperature :math:`T_v`, ground-skin temperature
:math:`T_g`, and atmospheric downwelling longwave :math:`LW_\downarrow`, the
calculation includes atmospheric absorption, canopy and ground emission, and
multiple reflection between the canopy and ground. The SLM radiation budget
uses positive absorbed energy, whereas the corresponding NoahMP intermediate
longwave terms use positive upward energy loss; the implementation changes the
sign when storing absorbed longwave energy.

The effective surface emissivity of the canopy-ground system is

.. math::

   \epsilon_{sfc} = \epsilon_v + \epsilon_g(1-\epsilon_v)
                    + \epsilon_v(1-\epsilon_v)(1-\epsilon_g).

The temperature supplied as the radiometric skin temperature is diagnosed from
the upward longwave flux and this effective emissivity. It is therefore not, in
general, equal to canopy temperature, ground-skin temperature, or top-soil
temperature.

5. Radiation assumptions and limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The current updated SLM radiation is limited to snow-free and ice-free land
surfaces. The source contains NoahMP snow-aging and BATS snow-albedo routines,
but snow cover, snow depth, snow water equivalent, and snowfall are currently
set to zero in the radiation driver. Lake and ice radiation paths are likewise
not active. Complete snow and ice treatment in SLM will be added in a separate
pull request.

The ground-albedo calculation currently uses a fixed soil-color category from
the NoahMP albedo table. The modified two-stream and BATS option identifiers
are also fixed internally rather than exposed as user inputs. These
restrictions should be considered when comparing SLM albedo or surface-energy
budgets with a complete NoahMP simulation.

6. Vegetation structure and seasonal state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Recent changes keep LAI and SAI separate instead of accumulating SAI into LAI.
Their sum is used only where a process depends on total vegetation area. This
distinction is required by the NoahMP optical properties, because leaves and
stems have different reflectances and transmittances.

When monthly tables are enabled, SLM linearly interpolates LAI and SAI between
the current and following month. Bare-soil cells retain zero LAI and SAI. After
each update, SLM recomputes the leaf-angle coefficients, vegetation emissivity,
precipitation-extinction coefficient, and canopy water-storage capacity.

Vegetation emissivity is now based on the leaf-angle extinction coefficients
:math:`\phi_1` and :math:`\phi_2`:

.. math::

   \epsilon_v = 0.97\left[1-
      \exp\left(-\left(\phi_1+\phi_2\right)LAI\right)\right].

The canopy water-storage capacity includes contributions from foliage and
woody area. In addition to the original LAI-dependent storage, the update uses
canopy height and basal area to represent water retained on trunks. Initial
surface albedos and emissivity are weighted by vegetation fraction so that
partially vegetated cells do not initially behave as closed canopy.

7. Aerodynamic and surface resistances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The updated surface exchange uses a resistance network connecting the
atmospheric reference level, canopy air, vegetation, and ground:

``r_a``
   Aerodynamic resistance between the effective surface and atmospheric
   reference level. The SLM/gSAM-SLM surface-layer calculation retains
   stability-dependent similarity functions, a wind-dependent roughness
   adjustment, and a small minimum turbulent velocity scale.

``r_d``
   NoahMP-derived aerodynamic resistance between the ground and canopy air.
   It uses canopy height, displacement height, vegetation area, friction
   velocity, and ``CWPVT`` to represent exponential wind attenuation through
   the canopy. For bare soil, ``r_d`` is set equal to ``r_a``.

``r_b``
   NoahMP-derived leaf boundary-layer resistance. It depends on canopy-top
   wind, wind extinction, and the characteristic leaf dimension ``DLEAF``.

``r_c``
   Leaf-level stomatal resistance. LAI is applied when leaf conductance is
   converted to canopy conductance, avoiding a second LAI factor in ``r_c``.

``r_soil``
   Soil-pore vapor resistance between the evaporating soil and the overlying
   air. The revised formulation is described in Section 10.

For a vegetated cell, the heat conductance from vegetation to canopy air is
proportional to :math:`2LAI/r_b`. Vapor conductance distinguishes wet-canopy
evaporation from dry-leaf transpiration and combines these paths with ground
evaporation through ``r_d`` and ``r_soil``.

8. Stomatal stress, root-zone stress, and transpiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The stomatal resistance uses a Noah/NoahMP Jarvis-style bulk formulation:

.. math::

   r_c = \min\left(r_{c,max},
         \frac{r_{c,min}}
         {f_{rad} f_{vpd} f_T f_{soil}}\right).

The four nondimensional factors represent incoming radiation, humidity
deficit, canopy-air temperature, and root-zone soil moisture. The radiation
factor uses incoming shortwave and ``RGL``. The vapor-pressure-deficit factor
uses the Noah-style mixing-ratio deficit,

.. math::

   f_{vpd} = \frac{1}{1 + HS\,(q_{sat}(T_{cas})-q_{cas})},

and the temperature factor is a quadratic response around ``slm.T_opt``.

Root-zone stress follows the NoahMP ``BTR_OPTION=1`` concept. ``NROOT`` gives
the number of soil layers that participate. In each rooted layer, liquid water
between wilting point and field capacity gives a linearly varying stress
factor. Layer thickness and moisture stress determine a normalized
``soil_transp_frac``. This fraction both contributes to the bulk soil stress
and distributes the transpiration sink among soil layers.

Transpiration is limited by water available above the wilting point in every
participating layer. This prevents a single time step from extracting more
water than the rooted soil profile can provide.

9. Canopy energy and water balance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The canopy temperature is prognostic. The dry canopy heat capacity follows the
NoahMP biomass parameterization,

.. math::

   C_{veg} = CBIOM\,\min(6,LAI+SAI)\,C_w,

where :math:`C_w` is the volumetric heat capacity of water. Water stored on the
canopy adds to the effective canopy heat capacity.

Precipitation is partitioned into intercepted water and throughfall using the
canopy extinction coefficient. Intercepted rain is assumed to arrive at the
atmospheric reference temperature and changes canopy temperature through
mixing. This rain-cooling treatment follows gSAM-SLM. Water in excess of the
canopy storage capacity drains to the ground.

Canopy sensible heat, wet-canopy evaporation or dew, and dry-leaf
transpiration are calculated from the resistance network. Wet-canopy
evaporation cannot remove more water than is stored on the canopy, and dry-leaf
transpiration cannot remove more water than is available in the rooted soil.
The canopy energy tendency is

.. math::

   C_{canopy}\frac{dT_v}{dt} = R_{n,v} - H_v - LE_v.

To control temperature oscillations when canopy heat capacity is small, this
balance is subcycled with approximately one-second substeps. The subcycling is
an SLM/gSAM-SLM numerical treatment rather than a NoahMP time-integration
scheme. Canopy temperature is limited to 343 K.

10. Ground-skin energy balance and soil evaporation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The updated SLM distinguishes four surface-related temperatures:

``t_ground_skin``
   Massless temperature at the ground-atmosphere interface.

``t_canop``
   Prognostic vegetation temperature.

``t_skin``
   Longwave-equivalent radiometric temperature of the combined canopy-ground
   system.

``tsoil``
   Temperature of the finite-volume soil layers, including the top soil layer
   immediately below the massless skin.

At each time step, SLM solves the ground surface-energy balance

.. math::

   0 = SW_{ground,abs} - LW_{ground,net} - H_g - L_v E_g - G,

where :math:`H_g` is ground sensible heat, :math:`E_g` is evaporation or dew,
and :math:`G` is conduction from the skin into the top soil layer. A bounded
Newton iteration updates ``t_ground_skin``. The converged conductive flux then
provides the upper-boundary forcing for the multilayer soil-temperature solve.

Ground evaporation uses the NoahMP option-1 soil resistance based on the dry
surface-layer formulation of Sakaguchi and Zeng. The effective dry-layer
thickness increases nonlinearly as normalized top-layer wetness decreases. The
vapor diffusivity depends on porosity, wilting-point wetness, and the soil
hydraulic exponent. The exponent controlling dry-layer growth is read from the
NoahMP ``RSURF_EXP`` parameter.

Very dry soil, urban surfaces, and invalid soil hydraulic states receive a
large soil resistance. Evaporation is further limited by the sum of ponded
water and water available in the top soil layer after transpiration. Urban
ground evaporation is multiplied by the pervious fraction. Condensation is
added to surface-water storage when the soil humidity state permits dew.

11. Soil and surface water budget
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The revised water budget separately tracks canopy evaporation, transpiration,
pond evaporation, and direct soil evaporation. Water reaching the ground is
the sum of canopy throughfall and canopy drainage.

Positive ground evaporation first removes water from surface pond storage.
Any remaining demand is applied to the top soil layer. Dew adds water to pond
storage. This ordering prevents pond and soil water from both supplying the
same evaporation flux.

Infiltration into unfrozen soil is limited by the top-layer saturated hydraulic
conductivity. It is reduced by the urban impervious fraction and set to zero
when all soil layers are saturated or the top layer is frozen. Transpiration is
removed from each rooted layer using ``soil_transp_frac``.

After the implicit soil-water solve, water above saturation is redistributed
toward deeper unsaturated layers while preserving total water. If the entire
profile cannot accept the excess, the rejected amount is returned to surface
pond storage by reducing the accepted infiltration. Pond water above
``slm.mws_mx0`` becomes surface drainage. This correction prevents rainwater
from disappearing when infiltration encounters a saturated profile.

Wetland cells are maintained at field-capacity wetness. Their surface water is
implicit in the wetland treatment, so explicit pond storage is cleared after
the soil-temperature update. Snow and ice hydrology are outside the scope of
the current update and will be addressed together with the separate snow/ice
SLM development.

Questions regarding SLM may be addressed to Jungmin Lee at Lawrence Livermore
National Laboratory (``lee1046@llnl.gov``).
