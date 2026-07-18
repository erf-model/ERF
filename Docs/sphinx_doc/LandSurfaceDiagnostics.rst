.. _sec:LandSurfaceDiagnostics:

Land-surface and 2-m diagnostics
================================

The primary ERF 2-D plotfile can expose the active land-surface state and a
unified diagnostic bundle at 2 m above the local surface. Select fields with
``erf.plot2d_vars_1`` or ``erf.plot2d_vars_2``. For example:

.. code-block:: text

   erf.plot2d_vars_1 = temperature_2m water_vapor_mixing_ratio_2m \
                       near_surface_diagnostic_source \
                       t_sfc albedo grdflx smois_1 tslb_1

ERF writes selected fields in canonical catalog order. A requested field that
the active provider cannot supply is skipped with the existing plot-variable
warning. Native AMReX 2-D plotfiles also contain ``2DMetadata.json`` with the
selected field name, long name, units, category, and missing-value policy.
See :ref:`sec:Plotfiles` for output frequency and file-format controls.

Unified 2-m fields
------------------

The unified fields are:

.. list-table::
   :header-rows: 1
   :widths: 34 46 20

   * - Name
     - Meaning
     - Units
   * - ``temperature_2m``
     - Physical air temperature 2 m above the local surface
     - K
   * - ``water_vapor_mixing_ratio_2m``
     - Water-vapor mass mixing ratio 2 m above the local surface, per unit dry-air mass
     - kg kg\ :sup:`-1` dry air
   * - ``near_surface_diagnostic_source``
     - Numeric code identifying the source selected for the unified bundle
     - categorical code

ERF chooses one coherent source for the fields requested in a plotfile stream.
It never combines a Noah-MP temperature with MOST humidity, or MOST
temperature with Noah-MP humidity. Selection is request-aware: temperature
needs the two native temperature tiles and ``f_v``; humidity needs the two
native humidity tiles, ``f_v``, and a moist run; requesting both needs all five
native fields. A source-only request follows the temperature path, including
in a dry run.

1. A valid native source for every requested field is preferred.
2. If the required native fields are invalid or incomplete, ERF evaluates all
   requested fields from one active SurfaceLayer MOST state when valid.
3. If neither source can satisfy the request, continuous fields are ``-999``
   and the categorical source code is ``0``.

Native Noah-MP aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`f_v` be the Noah-MP vegetation fraction. The grid-cell values are

.. math::

   T_2 = f_v T_{2,v} + (1-f_v)T_{2,b},

.. math::

   r_{v,2} = f_v r_{v,2,v} + (1-f_v)r_{v,2,b}.

The ERF Noah-MP transfer routine already converts the native specific
humidities to mixing ratios using :math:`r_v=q/(1-q)`. ERF therefore does not
apply a second conversion. Native temperature and humidity requests validate
only their required tile pair and ``f_v``. Temperature must be finite, outside
the provider fill range, and positive; requested mixing ratios must be finite,
outside the fill range, and nonnegative; and :math:`0\leq f_v\leq1`. Valid
zero and signed physical values are preserved.

MOST fallback
~~~~~~~~~~~~~

The fallback uses the same integrated scalar similarity function as the
SurfaceLayer implementation. At height :math:`z`:

.. math::

   \theta(z)-\theta_0 = \frac{\theta_*}{\kappa}
   \left[\ln\left(\frac{z}{z_0}\right)-\Psi_h\left(\frac{z}{L}\right)\right],

.. math::

   r_v(z)-r_{v,0} = \frac{r_{v,*}}{\kappa}
   \left[\ln\left(\frac{z}{z_0}\right)-\Psi_h\left(\frac{z}{L}\right)\right].

For these diagnostics, :math:`z=2` m above the local surface. ERF evaluates
the temperature and humidity scalar profiles with the current SurfaceLayer
similarity function. Temperature is converted from potential temperature with
the EOS helpers. Pressure is reconstructed from the lowest atmospheric cell
only when temperature is requested, using

.. math::

   p_2 = p_{cc}+\rho_{cc}g\left(z_{cc}^{AGL}-2\,\mathrm{m}\right).

This is a local first-order hydrostatic reconstruction. Invalid profile inputs,
pressure, or output values produce ``-999``; they are not silently clamped. A
humidity-only request does not require pressure reconstruction. The
implementation reuses the active SurfaceLayer similarity function described in
:ref:`sec:surface_layer`.

Source codes
------------

``near_surface_diagnostic_source`` uses the stable convention already used by
``surface_diagnostic_source``. When selected and supported, it is categorical
and every cell contains a code from ``0`` through ``6``; missing is ``0``
rather than the continuous ``-999`` fill. Selecting provenance alongside a
humidity-only request does not add a temperature prerequisite. A source-only
request retains the temperature-path convention. The field is not advertised
when neither Noah-MP nor SurfaceLayer can provide a diagnostic pathway:

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - Value
     - Name
     - Meaning
   * - 0
     - ``missing``
     - No source satisfied the selected request
   * - 1
     - ``surface_layer_land``
     - MOST supplied a land value without a valid native LSM bundle
   * - 2
     - ``lsm_land``
     - Native Noah-MP supplied every component required by the selected field or fields
   * - 3
     - ``surface_layer_fallback``
     - Noah-MP was selected, but components required by the selected field or fields were invalid or incomplete; MOST supplied the output
   * - 4
     - ``surface_layer_sea``
     - MOST supplied the value over water
   * - 5
     - ``custom``
     - The custom surface pathway supplied the SurfaceLayer state
   * - 6
     - ``rico``
     - The RICO surface pathway supplied the SurfaceLayer state

These values are part of the plotfile convention. Consumers should compare the
numeric code with this table rather than infer provenance from the land mask.

Native Noah-MP component fields
-------------------------------

The five provider-specific fields used by the unified aggregation are exposed
when Noah-MP is active. They remain missing where Noah-MP did not produce a
valid component; MOST never fills them.

.. list-table::
   :header-rows: 1
   :widths: 46 38 16

   * - Name
     - Meaning
     - Units
   * - ``noahmp_temperature_2m_vegetated``
     - 2-m temperature over the vegetated fraction
     - K
   * - ``noahmp_temperature_2m_bare``
     - 2-m temperature over the bare fraction
     - K
   * - ``noahmp_water_vapor_mixing_ratio_2m_vegetated``
     - 2-m water-vapor mixing ratio over the vegetated fraction
     - kg kg\ :sup:`-1` dry air
   * - ``noahmp_water_vapor_mixing_ratio_2m_bare``
     - 2-m water-vapor mixing ratio over the bare fraction
     - kg kg\ :sup:`-1` dry air
   * - ``noahmp_vegetation_fraction``
     - Vegetated fraction used by the two-tile aggregation
     - 1

Existing provider fields
------------------------

When Noah-MP is active, the catalog also exposes the stored provider fields,
including ``t_sfc``, ``sfc_emis``, the direct and diffuse visible/NIR albedos,
``cos_zenith_angle``, the shortwave and longwave downwelling fields, ``grdflx``,
``fira``, ``sav``, ``sag``, ``albedo``, ``sfcrunoff``, and ``udrunoff``. These
are provider values and are not reconstructed from MOST. Raw fields remain
missing when Noah-MP has not produced a valid provider value; MOST is used only
for the unified fields.

Dynamic soil layers
-------------------

Noah-MP soil fields are generated from its active layer inventory, preserving
the provider's order and configured layer count ``N``:

.. code-block:: text

   smois_1 ... smois_N
   sh2o_1  ... sh2o_N
   tslb_1  ... tslb_N

``smois`` and ``sh2o`` use volumetric units (m\ :sup:`3` m\ :sup:`-3`), while
``tslb`` uses K. Requested layers that are not present are unavailable and use
the normal plot-variable warning path. The runtime descriptors are generated
only for the active model and carry the same metadata and ``-999`` missing
policy as the fixed catalog entries.

Timing and limitations
----------------------

The unified fields are assembled when each primary 2-D plotfile is written,
using the latest level's native LSM fields, SurfaceLayer state, source mask,
land mask, and lowest-cell atmospheric state. They are not checkpoint fields
and do not alter surface-flux updates. The native path adds five persistent
one-component 2-D fields per active level; those fields are interpolated and
initialized with the provider's undefined sentinel before the first valid
Noah-MP result.

SurfaceLayer can provide a physically meaningful atmospheric profile near the
surface, but it cannot replace a land model's soil, canopy, snow, hydrology, or
radiation state. For Noah-MP coupling details, see :ref:`CouplingToNoahMP`.
