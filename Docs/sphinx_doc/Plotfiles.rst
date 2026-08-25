.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

*********
Plotfiles
*********

ERF provides three plotfile output paths:

* Standard 3D plotfiles write selected fields on all active AMR levels.
* 2D plotfiles write one-cell-thick horizontal slabs containing fixed
  diagnostics and, optionally, fields sampled on user-defined vertical levels.
* Subvolume plotfiles write selected fields from one or more regular 3D regions.

The standard 3D and 2D paths each support two independent output streams. Each
stream has its own format, prefix, schedule, and variable list. Subvolume output
uses one prefix and one variable list, but it may define more than one region
and schedule.

Choosing an output type
=======================

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Output
     - Use it for
     - Main controls
   * - Standard 3D
     - Full model fields on every active AMR level.
     - ``erf.plot_file_*``, ``erf.plot_int_*``, ``erf.plot_per_*``, and
       ``erf.plot_vars_*``.
   * - 2D diagnostics
     - Surface, column, and selected vertical-level fields on horizontal slabs.
     - ``erf.plot2d_file_*``, ``erf.plot2d_int_*``, ``erf.plot2d_per_*``,
       ``erf.plot2d_vars_*``, and ``erf.plot2d_level_sets_*``.
   * - 3D subvolume
     - One or more regular regions smaller than the full model domain.
     - ``erf.subvol_file``, ``erf.subvol_int``, ``erf.subvol_per``,
       ``erf.subvol_sampling_vars``, and the ``erf.subvol.*`` geometry controls.

Choosing AMReX or NetCDF
========================

ERF selects the 3D and 2D formats independently. Format names are parsed
without regard to case. Native AMReX is the effective default for every
standard 3D and 2D stream. NetCDF requires an ERF build with NetCDF enabled.

.. list-table::
   :header-rows: 1
   :widths: 32 52 16

   * - Parameter
     - Meaning
     - Effective default
   * - ``erf.plotfile_type``
     - Backward-compatible shorthand that sets both 3D streams. Do not combine
       it with ``erf.plotfile_type_1`` or ``erf.plotfile_type_2``.
     - ``amrex``
   * - ``erf.plotfile_type_1``
     - Format for the first 3D stream: ``amrex`` or ``netcdf``.
     - ``amrex``
   * - ``erf.plotfile_type_2``
     - Format for the second 3D stream: ``amrex`` or ``netcdf``.
     - ``amrex``
   * - ``erf.plotfile2d_type``
     - Shorthand that sets both 2D streams. Do not combine it with
       ``erf.plotfile2d_type_1`` or ``erf.plotfile2d_type_2``.
     - ``amrex``
   * - ``erf.plotfile2d_type_1``
     - Format for the first 2D stream: ``amrex`` or ``netcdf``.
     - ``amrex``
   * - ``erf.plotfile2d_type_2``
     - Format for the second 2D stream: ``amrex`` or ``netcdf``.
     - ``amrex``

If a shorthand and its corresponding per-stream parameter are both set, ERF
stops with an input error. A per-stream parameter may be set without setting the
other stream; the unset stream keeps the AMReX default.

The configured stream formats apply to ordinary scheduled and startup writes.
Final-time writes use each stream's configured format independently.

Primary native AMReX 2D and 3D plotfiles contain an execution and ancestry
record in ``job_info``. Each output has its own artifact UUID. Native AMReX 2D
plotfiles also write ``2DMetadata.json``. NetCDF output uses the same public 2D
variable names but does not provide the equivalent sidecar metadata. See
:ref:`sec:Provenance` and :ref:`sec:Plotfile2DMetadata`.

Subvolume output uses the native AMReX single-level plotfile writer. The
standard 3D and 2D format parameters do not select a NetCDF subvolume format.

A workflow that needs NetCDF may write native AMReX plotfiles during the run and
convert them afterward with the tools under ``Exec/Tools``. Build those tools
with GNU Make or with ``ERF_ENABLE_TOOLS`` in CMake.

.. _sec:PlotfileCommonParameters:

Common naming controls
======================

.. list-table::
   :header-rows: 1
   :widths: 30 48 12 10

   * - Parameter
     - Meaning
     - Type
     - Default
   * - ``erf.use_real_time_in_pltname``
     - For 3D primary and subvolume outputs, use a UTC calendar timestamp
       derived from the configured start time plus elapsed model time. 2D
       output and face-velocity files remain step-numbered.
     - Boolean
     - ``false``
   * - ``erf.file_name_digits``
     - Number of digits appended to step-based plotfile and checkpoint names.
     - Integer greater than zero
     - ``5``

When ``erf.use_real_time_in_pltname = true``, a primary 3D name uses
``<prefix>_YYYY-MM-DD_HH:MM:SS`` from ``start_time + elapsed model time`` in
UTC. A subvolume name uses the same timestamp with six fractional-second digits
after ``<subvol_file>_<region-index>``. 2D primary names and 3D face-velocity
names remain step-numbered.

Output streams and schedules
============================

The 3D and 2D paths each provide stream 1 and stream 2. Configure either stream
independently. An ``*_int_*`` parameter schedules output by level-0 time steps.
An ``*_per_*`` parameter schedules output by simulation time. Leave an unused
schedule at its negative default. ERF rejects a 3D or 2D stream for which both
the interval and period are positive.

For subvolumes, specify ``erf.subvol_int`` or ``erf.subvol_per``, not both. The
parser accepts one schedule value or one value per region. With several regions,
a single value applies only to the first region; supply one value per region to
schedule every region.

Configuring 3D output
=====================

.. list-table::
   :header-rows: 1
   :widths: 28 48 12 12

   * - Parameter
     - Meaning
     - Type
     - Default
   * - ``erf.plot_file_1``
     - Prefix for the first 3D stream.
     - String
     - ``plt_1_``
   * - ``erf.plot_file_2``
     - Prefix for the second 3D stream.
     - String
     - ``plt_2_``
   * - ``erf.plot_int_1``
     - Level-0 step interval for the first 3D stream.
     - Integer greater than zero
     - ``-1``
   * - ``erf.plot_int_2``
     - Level-0 step interval for the second 3D stream.
     - Integer greater than zero
     - ``-1``
   * - ``erf.plot_per_1``
     - Simulation-time period for the first 3D stream.
     - Real greater than zero
     - ``-1.0``
   * - ``erf.plot_per_2``
     - Simulation-time period for the second 3D stream.
     - Real greater than zero
     - ``-1.0``
   * - ``erf.plot_vars_1``
     - Variables in the first 3D stream.
     - List of names
     - None
   * - ``erf.plot_vars_2``
     - Variables in the second 3D stream.
     - List of names
     - None
   * - ``erf.plot_face_vels``
     - Also write the native staggered velocity components.
     - Boolean
     - ``false``

Example:

.. code-block:: text

   erf.plotfile_type = amrex
   erf.plot_file_1   = plt
   erf.plot_int_1    = 10
   erf.plot_vars_1   = density theta qv x_velocity y_velocity z_velocity

This writes the selected 3D fields every 10 level-0 steps. See
:ref:`sec:Plotfile3DReference` for the complete fixed core inventory, dynamic
provider rules, and physics restrictions.

Configuring 2D output
=====================

.. list-table::
   :header-rows: 1
   :widths: 30 46 12 12

   * - Parameter
     - Meaning
     - Type
     - Default
   * - ``erf.plot2d_file_1``
     - Prefix for the first 2D stream.
     - String
     - ``plt2d_1_``
   * - ``erf.plot2d_file_2``
     - Prefix for the second 2D stream.
     - String
     - ``plt2d_2_``
   * - ``erf.plot2d_int_1``
     - Level-0 step interval for the first 2D stream.
     - Integer greater than zero
     - ``-1``
   * - ``erf.plot2d_int_2``
     - Level-0 step interval for the second 2D stream.
     - Integer greater than zero
     - ``-1``
   * - ``erf.plot2d_per_1``
     - Simulation-time period for the first 2D stream.
     - Real greater than zero
     - ``-1.0``
   * - ``erf.plot2d_per_2``
     - Simulation-time period for the second 2D stream.
     - Real greater than zero
     - ``-1.0``
   * - ``erf.plot2d_vars_1``
     - Built-in diagnostics in the first 2D stream.
     - List of names
     - None
   * - ``erf.plot2d_vars_2``
     - Built-in diagnostics in the second 2D stream.
     - List of names
     - None
   * - ``erf.plot2d_level_sets_1``
     - Sampled-level sets appended to the first 2D stream.
     - List of level-set names
     - None
   * - ``erf.plot2d_level_sets_2``
     - Sampled-level sets appended to the second 2D stream.
     - List of level-set names
     - None

Example:

.. code-block:: text

   erf.plotfile2d_type = amrex
   erf.plot2d_file_1   = plt2d_
   erf.plot2d_int_1    = 10
   erf.plot2d_vars_1   = z_surf mapfac lat_m lon_m u_star surf_pres integrated_qv

ERF writes fixed diagnostics in canonical catalog order, not request order. It
warns and skips unknown or configuration-ineligible names. See
:ref:`sec:Plotfile2DReference` for the complete fixed and dynamic diagnostic
contract. See :ref:`sec:Plotfile2DSampledLevels` to add pressure, height, or
model-index fields.

Configuring subvolume output
============================

Subvolume output may define one or more regions. The three geometry arrays must
contain the same number of values, and each must contain three values per
region. The requested spacing must match an existing AMR level. The origin must
coincide with a node on that level, and the requested box must lie within the
level domain.

.. list-table::
   :header-rows: 1
   :widths: 30 46 12 12

   * - Parameter
     - Meaning
     - Type
     - Default
   * - ``erf.subvol_file``
     - Prefix shared by subvolume outputs; ERF appends the region index.
     - String
     - ``subvol``
   * - ``erf.subvol_int``
     - Level-0 step interval. Supply one value per configured region when
       scheduling every region.
     - Integer list
     - ``-1``
   * - ``erf.subvol_per``
     - Simulation-time period. Supply one value per configured region when
       scheduling every region.
     - Real list
     - ``-1.0``
   * - ``erf.subvol_sampling_vars``
     - Variables written in every configured subvolume.
     - List of names
     - ``x_velocity y_velocity z_velocity``
   * - ``erf.subvol.origin``
     - Lower corner of each output region.
     - Three real values per region
     - Required when enabled
   * - ``erf.subvol.nxnynz``
     - Number of output cells in each direction for each region.
     - Three integers per region
     - Required when enabled
   * - ``erf.subvol.dxdydz``
     - Output spacing in each direction for each region.
     - Three real values per region
     - Required when enabled
   * - ``erf.subvol.chunk_size``
     - Optional maximum box size used to decompose the subvolume output.
     - Three integers
     - Model maximum grid size

Subvolume output supports a restricted variable inventory, not the complete 3D
plotfile catalog. The default fields are the three cell-centered velocity
components. ``erf.subvol_sampling_vars`` may select supported conserved state
variables and the derived fields ``temp``, ``theta``, ``KE``, ``scalar``,
``soundspeed``, ``precipitable``, and ``mucape``. The last two are available
only for a moist run. Conserved-state availability is decided by exactly the
predicate the 3D plotfile uses, so the two output paths accept the same
``rhoQn`` names for a given moisture scheme; unsupported or unavailable names
are omitted from the output rather than reported as an error.

.. _sec:PlotfileGeneralNotes:

General behavior
================

* Native AMReX 3D plotfiles contain all active refinement levels in one
  plotfile hierarchy. NetCDF output is written separately by level.
* NetCDF output requires an ERF build with NetCDF enabled.
* A negative interval or period disables that scheduling control.
* Constructing fixed and sampled 2D diagnostics does not change prognostic
  state.
* A 3D write for Lagrangian microphysics with two-way AMR coupling performs a
  fine-to-coarse consistency average-down before constructing the output. This
  updates covered coarse state and microphysics storage; it is not a physical
  time tendency.

.. _sec:PlotfileDetailedReferences:

Detailed references
===================

.. toctree::
   :maxdepth: 1

   plotfiles/Plotfile3DReference
   plotfiles/Plotfile2DReference
   plotfiles/Plotfile2DSampledLevels
   plotfiles/LandSurfaceDiagnostics

The reference pages separate variable lookup, sampled-level configuration, and
land-surface source-selection details. Use them when the common examples above
are not sufficient.
