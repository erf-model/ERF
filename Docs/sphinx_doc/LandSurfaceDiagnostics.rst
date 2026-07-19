.. _sec:LandSurfaceDiagnostics:

Land-surface and 2-m diagnostics
================================

This page defines the source-selection and numerical contracts behind ERF's
land-surface and unified 2-m diagnostics. For the complete list of built-in
request names, canonical order, metadata, selection rules, and runtime value
behavior, see :ref:`sec:Plotfile2DBuiltInCatalog` and
:ref:`sec:Plotfile2DSelectionRules`.

Unified 2-m fields and coherent source selection
------------------------------------------------

The public unified fields are:

``temperature_2m``
   Physical temperature 2 m above the local surface.

``water_vapor_mixing_ratio_2m``
   Water-vapor mixing ratio 2 m above the local surface, per unit dry-air
   mass.

``near_surface_diagnostic_source``
   Categorical code for the source selected for the requested bundle.

Selection is cell-local and request-aware. ERF chooses one coherent source for
all requested continuous fields. It never combines native temperature with
MOST humidity, or MOST temperature with native humidity.

The native requirements are minimal for the request:

* A temperature-only request needs the two native temperature tiles and
  vegetation fraction.
* A humidity-only request needs the two native humidity tiles, vegetation
  fraction, and a moist run.
* A combined request needs all five native fields.
* A source-only request follows the temperature path, including in a dry run.

ERF prefers a valid native source for every requested field. If the required
native fields are invalid or incomplete, ERF evaluates all requested fields
from one complete SurfaceLayer/MOST state. If neither source satisfies the
request, continuous fields are ``-999`` and the source code is ``0``.

Native Noah-MP aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`f_v` be the Noah-MP vegetation fraction. The native tile values are
aggregated as

.. math::

   T_2 = f_v T_{2,v} + (1-f_v)T_{2,b},

and

.. math::

   r_{v,2} = f_v r_{v,2,v} + (1-f_v)r_{v,2,b}.

The Noah-MP transfer layer performs the one-time conversion from native
specific humidity to dry-air mixing ratio:

.. math::

   r_v = \frac{q_v}{1-q_v}.

The 2-D diagnostic layer receives mixing ratio and does not convert it again.

Native validation is request-specific:

* temperature tiles must be finite, outside provider fill ranges, and
  positive;
* mixing-ratio tiles must be finite, outside provider fill ranges, and
  nonnegative;
* vegetation fraction must be finite, outside provider fill ranges, and lie in
  ``[0, 1]``; and
* the native request succeeds only when every component required by that
  request is valid.

MOST fallback
~~~~~~~~~~~~~

ERF reuses the active SurfaceLayer integrated similarity function. Define the
profile factor

.. math::

   F_h(z) = \ln\left(\frac{z}{z_0}\right)
            - \Psi_h\left(\frac{z}{L}\right).

The scalar profiles are

.. math::

   \theta(z) = \theta_0 + \frac{\theta_*}{\kappa}F_h(z),

.. math::

   r_v(z) = r_{v,0} + \frac{r_{v,*}}{\kappa}F_h(z).

For these outputs, ERF uses ``z = 2 m`` above the local surface. Temperature
uses EOS helpers after reconstructing local pressure:

.. math::

   p_2 = p_{cc} + \rho_{cc}g
         \left(z_{cc}^{\mathrm{AGL}} - 2\,\mathrm{m}\right).

A humidity-only request does not require pressure, terrain height, conserved
atmospheric state, or a temperature source. Invalid profile geometry, pressure,
state, or output produces ``-999``. ERF does not silently clamp an invalid
profile. MOST fills only unified fields; it never fills raw provider fields.

Exact HFX processed-cell scope
------------------------------

Noah-MP ``HFX`` gates the following Noah-MP return fields:

.. code-block:: text

   t_sfc
   sfc_emis
   sfc_alb_dir_vis
   sfc_alb_dir_nir
   sfc_alb_dif_vis
   sfc_alb_dif_nir
   grdflx
   fira
   sav
   sag
   albedo
   sfcrunoff
   udrunoff
   noahmp_temperature_2m_vegetated
   noahmp_temperature_2m_bare
   noahmp_water_vapor_mixing_ratio_2m_vegetated
   noahmp_water_vapor_mixing_ratio_2m_bare
   noahmp_vegetation_fraction

It also gates Noah-MP surface-flux results and runtime soil result fields. An
invalid ``HFX`` means that Noah-MP did not process the cell. ERF invalidates
that result set. A valid ``HFX`` only establishes that the cell was processed;
each returned field still passes its own validity rule.

The following radiation-to-LSM forcing fields are not HFX-gated:

.. code-block:: text

   cos_zenith_angle
   sw_flux_dn
   sw_flux_dn_dir_vis
   sw_flux_dn_dir_nir
   sw_flux_dn_dif_vis
   sw_flux_dn_dif_nir
   lw_flux_dn

They are stored through the radiation exchange path. Their public 2-D output
still translates absent, nonfinite, or sentinel values to ``-999``.

Provider-specific native fields
-------------------------------

Native Noah-MP 2-m component fields are raw provider outputs. They are not
area-aggregated public 2-m diagnostics, and MOST never fills them. Radiation
exchange fields and Noah-MP return fields use different production paths. The
canonical catalog in :ref:`sec:Plotfile2DBuiltInCatalog` gives the complete
request inventory.

Dynamic soil fields
-------------------

Dynamic soil names are generated from the active land-surface provider
inventory. The layer index is one-based, and runtime order follows that
inventory. See :ref:`sec:Plotfile2DDynamicSoil` for the family table and
metadata contract.

Timing and limitations
----------------------

ERF assembles these fields when it writes each primary 2-D plotfile. It
assembles each AMR level from that level's latest available land-surface,
SurfaceLayer, source-mask, terrain, and atmospheric state. These fields are
diagnostics, not checkpoint state. They do not change prognostic state or
surface-flux updates.

SurfaceLayer can provide a near-surface atmospheric profile, but it cannot
replace soil, canopy, snow, hydrology, or radiation state.
