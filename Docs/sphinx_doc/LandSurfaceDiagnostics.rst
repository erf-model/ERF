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
     - 1

ERF chooses one complete source bundle per horizontal cell. It never combines
a Noah-MP temperature with MOST humidity, or MOST temperature with Noah-MP
humidity.

1. A valid native Noah-MP bundle is preferred.
2. If any required native component is invalid, ERF evaluates both fields from
   the active SurfaceLayer MOST state when that state is valid.
3. If neither source produces a valid value, ERF writes ``-999``.

Native Noah-MP aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`f_v` be the Noah-MP vegetation fraction. The grid-cell values are

.. math::

   T_2 = f_v T_{2,v} + (1-f_v)T_{2,b},

.. math::

   r_{v,2} = f_v r_{v,2,v} + (1-f_v)r_{v,2,b}.

The ERF Noah-MP transfer routine already converts the native specific
humidities to mixing ratios using :math:`r_v=q/(1-q)`. ERF therefore does not
apply a second conversion. The native bundle is accepted only when all five
components are finite, outside the provider fill range, :math:`0\leq f_v\leq1`,
and both mixing ratios are nonnegative.

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

For these diagnostics, :math:`z=2` m above the local surface. ERF converts the
resulting potential temperature with the EOS helpers. Pressure is reconstructed
from the lowest atmospheric cell using

.. math::

   p_2 = p_{cc}+\rho_{cc}g\left(z_{cc}^{AGL}-2\,\mathrm{m}\right).

This is a local first-order hydrostatic reconstruction. Invalid profile inputs,
pressure, or output values produce ``-999``; they are not silently clamped.
The implementation reuses the active SurfaceLayer similarity function described
in :ref:`sec:surface_layer`.

Source codes
------------

``near_surface_diagnostic_source`` uses the stable convention already used by
``surface_diagnostic_source``:

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - Value
     - Name
     - Meaning
   * - 0
     - ``missing``
     - Neither native Noah-MP nor MOST produced a valid value
   * - 1
     - ``surface_layer_land``
     - MOST supplied a land value without a valid native LSM bundle
   * - 2
     - ``lsm_land``
     - Native Noah-MP supplied the complete valid 2-m bundle
   * - 3
     - ``surface_layer_fallback``
     - Noah-MP was selected, but its bundle was invalid and MOST supplied the value
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
are provider values and are not reconstructed from MOST.

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
using the current level's native LSM fields, SurfaceLayer state, source mask,
land mask, and lowest-cell atmospheric state. They are not checkpoint fields
and do not alter surface-flux updates. The native path adds five persistent
one-component 2-D fields per active level; those fields are interpolated and
initialized with the provider's undefined sentinel before the first valid
Noah-MP result.

SurfaceLayer can provide a physically meaningful atmospheric profile near the
surface, but it cannot replace a land model's soil, canopy, snow, hydrology, or
radiation state. For Noah-MP coupling details, see :ref:`CouplingToNoahMP`.
