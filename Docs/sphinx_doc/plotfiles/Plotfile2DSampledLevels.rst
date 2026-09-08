.. _sec:Plotfile2DSampledLevels:

2D sampled-level output
=======================

ERF can append fields sampled on model-index, height, or pressure targets to a
2D output stream. Use a named level set to define the vertical coordinate,
targets, fields, interpolation rule, and missing value.

For example, this configuration samples potential temperature, temperature,
water vapor, and cloud water at 850, 700, and 500 hPa:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air

   erf.plot2d.level_set.upper_air.coordinate = pressure
   erf.plot2d.level_set.upper_air.units = hPa
   erf.plot2d.level_set.upper_air.values = 850 700 500
   erf.plot2d.level_set.upper_air.fields = theta temp qv qc

ERF appends these fields after the fixed diagnostics selected for stream 1. It
does not extrapolate beyond the model column.

Configuring a level set
-----------------------

The level-set keys define the coordinate, target values, fields, interpolation,
and missing-value behavior for a named set.

Select sampled level sets for each 2D output stream with:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air bl_heights
   erf.plot2d_level_sets_2 = native_k

Define each level set under ``erf.plot2d.level_set.<name>.``.

Required keys:

.. code-block:: text

   erf.plot2d.level_set.<name>.coordinate = ...
   erf.plot2d.level_set.<name>.values = ...
   erf.plot2d.level_set.<name>.fields = ...

Optional keys:

.. code-block:: text

   erf.plot2d.level_set.<name>.units = ...
   erf.plot2d.level_set.<name>.interpolation = ...
   erf.plot2d.level_set.<name>.missing_value = ...

The default missing value is ``-999``. ERF does not extrapolate sampled-level
diagnostics. If a target lies outside the column, ERF writes the level set's
missing value.

Supported vertical coordinates
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 48 16 16

   * - Coordinate
     - Meaning
     - Units
     - Interpolation
   * - ``model_index``
     - Cell-centered model level index.
     - ``1``
     - ``none``
   * - ``height_msl``
     - Cell-centered height above mean sea level.
     - ``m``
     - ``linear``
   * - ``height_agl``
     - Cell-centered height above local terrain.
     - ``m``
     - ``linear``
   * - ``pressure``
     - Cell-centered pressure.
     - ``Pa`` or ``hPa``
     - ``linear``

ERF converts pressure targets in ``hPa`` to canonical ``Pa`` values in metadata.
Model-index values must be integers. ERF rejects unsupported coordinates and
units during input parsing. Isentropic output is not enabled because it needs a
crossing policy for non-monotonic potential-temperature columns.

Supported fields
----------------

.. list-table::
   :header-rows: 1
   :widths: 18 52 14 26

   * - Field
     - Meaning
     - Units
     - Availability
   * - ``rho``
     - Density.
     - ``kg m^-3``
     - Always available.
   * - ``theta``
     - Potential temperature.
     - ``K``
     - Always available.
   * - ``temp``
     - Temperature.
     - ``K``
     - Always available.
   * - ``pressure``
     - Pressure.
     - ``Pa``
     - Always available.
   * - ``height_msl``
     - Height above mean sea level.
     - ``m``
     - Always available.
   * - ``height_agl``
     - Height above local terrain.
     - ``m``
     - Always available.
   * - ``qv``
     - Water vapor mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qv``.
   * - ``qc``
     - Cloud liquid water mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qc``.
   * - ``qi``
     - Cloud ice mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qi``.
   * - ``qr``
     - Rain mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qr``.
   * - ``qs``
     - Snow mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qs``.
   * - ``qg``
     - Graupel mixing ratio.
     - ``kg kg^-1``
     - Available when the active moisture model has ``qg``.
   * - ``u_east``
     - Eastward wind.
     - ``m/s``
     - Always available.
   * - ``v_north``
     - Northward wind.
     - ``m/s``
     - Always available.
   * - ``w``
     - Vertical wind.
     - ``m/s``
     - Always available.
   * - ``wind_speed``
     - Horizontal wind speed.
     - ``m/s``
     - Always available.
   * - ``wind_dir``
     - Meteorological wind direction.
     - ``degrees``
     - Always available. ERF writes the missing value for calm winds.

Use the canonical sampled species names ``qv``, ``qc``, ``qi``, ``qr``, ``qs``,
and ``qg``. Do not use 3D derived-variable aliases such as ``qrain``, ``qsnow``,
or ``qgraup`` in sampled-level output.

If a requested moisture field is unavailable for the active moisture model, ERF
warns and skips that field. If all requested fields in a level set are
unavailable, ERF aborts with an input error.

Output names and ordering
-------------------------

Sampled-level output names follow this pattern:

.. code-block:: text

   <field>_<coordinate_tag>_<value_tag>

Examples include ``theta_p_850hPa``, ``qv_p_700hPa``, ``qc_z_agl_500m``,
``theta_z_msl_1000m``, ``theta_k_10``, ``u_east_p_850hPa``,
``v_north_p_850hPa``, ``w_z_agl_500m``, ``wind_speed_k_10``,
``wind_speed_z_agl_500m``, and ``wind_dir_k_10``.

ERF writes sampled-level variables after built-in variables in the same 2D
stream. Within sampled-level output, ERF writes variables in level-set order,
target-value order, and field order.

Interpolation and missing values
--------------------------------

For linear sampling, ERF finds adjacent cell centers :math:`k_0` and
:math:`k_1` that bracket the target coordinate :math:`C_t`.

.. math::

   F_t = (1 - w) F_{k_0} + w F_{k_1}

with

.. math::

   w = \frac{C_t - C_{k_0}}{C_{k_1} - C_{k_0}}.

The bracket test works for increasing coordinates, such as height, and
decreasing coordinates, such as pressure. Model-index sampling copies the
requested level exactly.

Moisture fields
---------------

Microphysical mass species are sampled as mixing ratios. ERF stores each
species as a conserved density :math:`\rho q_x`, so sampled output uses

.. math::

   q_x = \frac{\rho q_x}{\rho}.

ERF uses ``qv = 0`` for dry-run pressure and temperature calculations, but it
does not expose sampled ``qv`` unless the active moisture model has a water-vapor
state component.

Wind fields
-----------

Wind fields are cell-centered sampled-level diagnostics. ERF destaggers the
native face-centered velocity components to scalar cell centers before vertical
sampling. ``u_east`` and ``v_north`` are earth-relative horizontal winds. When
map-rotation coefficients are available, ERF rotates grid-relative horizontal
winds before output. Otherwise, ERF uses identity rotation. ``w`` is the
cell-centered vertical wind.

``wind_speed`` is the horizontal speed computed from ``u_east`` and
``v_north``. ``wind_dir`` is the meteorological wind direction in degrees
clockwise from north, indicating where the wind comes from. ERF writes the
level-set missing value for ``wind_dir`` when the horizontal wind is calm.
``wind_dir`` is derived from the interpolated vector; ERF does not vertically
interpolate wind direction as a scalar angle.

.. _sec:Plotfile2DSampledMetadata:

AMReX metadata and NetCDF limitations
--------------------------------------

Sampled-level metadata records ``source_field`` and ``vertical_coordinate`` in
the native AMReX ``2DMetadata.json`` sidecar. NetCDF 2D output uses the same
sampled-level variable names, but it does not write sampled-level metadata
attributes.

The sampled-level record has this form:

.. code-block:: json

   {
     "format_version": 2,
     "kind": "ERF 2D plotfile metadata",
     "n_variables": 1,
     "variables": [
       {
         "component_index": 0,
         "name": "theta_p_850hPa",
         "long_name": "Potential temperature sampled on pressure levels",
         "units": "K",
         "category": "SampledLevel",
         "missing_policy": "FillMinus999WhenUnavailable",
         "missing_value": -999,
         "source_field": "theta",
         "vertical_coordinate": {
           "type": "pressure",
           "value": 850,
           "units": "hPa",
           "canonical_value": 85000,
           "canonical_units": "Pa",
           "interpolation": "linear"
         }
       }
     ]
   }

Additional examples
-------------------

Height-above-ground example:

.. code-block:: text

   erf.plot2d_level_sets_1 = bl_heights

   erf.plot2d.level_set.bl_heights.coordinate = height_agl
   erf.plot2d.level_set.bl_heights.values = 100 500 1000
   erf.plot2d.level_set.bl_heights.fields = theta qv qc

Model-index example:

.. code-block:: text

   erf.plot2d_level_sets_1 = native_k

   erf.plot2d.level_set.native_k.coordinate = model_index
   erf.plot2d.level_set.native_k.values = 0 10 20
   erf.plot2d.level_set.native_k.fields = theta qv qc

Wind example:

.. code-block:: text

   erf.plot2d_level_sets_1 = upper_air_winds

   erf.plot2d.level_set.upper_air_winds.coordinate = pressure
   erf.plot2d.level_set.upper_air_winds.units = hPa
   erf.plot2d.level_set.upper_air_winds.values = 850 700
   erf.plot2d.level_set.upper_air_winds.fields = u_east v_north wind_speed wind_dir

Current limitations
-------------------

Isentropic levels, vorticity, potential vorticity, and staggered velocity fields
need additional sampling rules. They are not part of this sampled-level output
mode.
