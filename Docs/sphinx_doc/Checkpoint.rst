.. role:: cpp(code)
  :language: c++

.. _sec:Checkpoint:

********************
Checkpoint / Restart
********************
.. toctree::
   :maxdepth: 1

ERF has a standard sort of checkpointing and restarting capability and
uses the native AMReX format for reading and writing checkpoints.
Each native checkpoint contains a provenance block in ``job_info``. The block
records the checkpoint artifact, the ERF execution that wrote it, and the
known restart ancestry. See :ref:`sec:Provenance`.
New native checkpoints also contain a small ``surface_temperature_contract``
marker with value ``2``. This records that the SST and skin temperature arrays
(``SST_0``, ``TSK_0``) hold absolute temperature [K], which the SurfaceLayer
converts to potential temperature at the current surface pressure each time it
uses them, so they are restored directly, whatever the restart reads.

Older checkpoints may hold these arrays as potential temperature, converted
with the surface pressure of the input they were read from: those with marker
value ``1``, and markerless ones written by a WRFInput run. (Markerless arrays
written by a Metgrid run hold absolute temperature and are restored directly,
and markerless checkpoints without SST/TSK arrays remain compatible.) On
restart these are converted to absolute temperature with a pressure chosen by
the ``erf.init_type`` of the run that *wrote* the checkpoint, as its own
``job_info`` records it, never by the ``erf.init_type`` of the restart, which
need not be the same:

- **WRFInput: best-available pressure, approximate on fine levels.** Level 0
  reads ``PSFC`` from the first file in ``erf.nc_init_file_0``. A smaller
  level-0 checkpoint domain reads only the corresponding local regions of that
  file. If every covered valid file sample read for the model level is exactly
  zero, ERF replaces those covered values with ``p_0``. The decision is made
  across the level's covered valid samples, not separately for each local box;
  uncovered parent-seeded cells are excluded and left unchanged. Cold
  initialization uses a level-MultiFab ``max == 0`` check; migration instead
  tests every covered file sample for exact zero. These checks are not generally
  equivalent on cropped or partially covered levels. Migration's check ensures
  negative or non-finite covered samples are not hidden by the fallback. A
  zero-only local box in a level that also has positive covered samples, or any
  negative or non-finite covered value, does not use this fallback and remains
  subject to pressure validation. On a finer level, an own-file ``PSFC`` overlays
  the part of the regridded level covered by that file; newly covered cells
  outside its footprint use interpolated parent pressure. A level without its
  own surface forcing uses pressure interpolated from its parent. In particular, when
  level 0 reads ``erf.nc_low_file`` and a finer level has no own
  ``erf.nc_low_file_<lev>``, its SST/TSK slices are inherited from the parent,
  and migration uses the parent's pressure.

  Fine-level recovery is approximate even when the own-file pressure is
  available. ``RemakeLevel`` can replace some checkpointed fine values with
  parent-interpolated potential temperature, but the checkpoint does not
  record which cells were replaced. Interpolating potential temperature and
  converting it with pressure do not commute, so the original absolute
  temperature cannot be recovered exactly. Migration warns about this
  uncertainty and applies the file-footprint/parent-pressure policy above.
  This needs a restart with ``erf.init_type = WRFInput``, a build with NetCDF,
  and the relevant wrfinput files; otherwise migration aborts.

  A finer level uses parent pressure if it had no own init file, or if level 0
  reads ``erf.nc_low_file`` and this level has no own
  ``erf.nc_low_file_<lev>``. Otherwise, a level with its own init file uses
  that file's ``PSFC`` if it has only one init file, or if no later file may
  have started. Once a later file has started, the checkpoint does not record
  which file each grid used or whether regridding replaced its values, so
  migration safely aborts if pressure provenance cannot be established.

  Migration converts every positive entry, including a TSK that the current
  surface layer does not use, so the arrays cannot become mixed-unit after a
  later regrid or forcing change. A conservative false abort is therefore
  possible when TSK is present without an SST and a fine level's own
  ``wrflowinp`` prevents parent inheritance, while multiple init files leave
  the TSK's original pressure ambiguous. We keep the abort rather than leave
  positive potential temperature in an array labeled kelvin. This is a
  compatibility limitation that needs maintainer acceptance; allowing the
  restart would require a policy for a value whose original pressure is
  unknown.
- **Metgrid: approximate migration.** The met_em input pressure is not stored
  in the checkpoint or read again. For boxes reaching the domain ground,
  migration diagnoses pressure from the restarted ground-state cells. For a
  refined level whose boxes are aloft, it uses interpolated parent ground
  pressure; an elevated lowest cell is not a surface-pressure source. If that
  fine level later reaches the ground, the Kelvin-valued migrated field is
  retained and the surface layer converts it using the then-current local
  ground pressure. Migration warns that the original absolute temperature
  cannot be recovered exactly because the input pressure and regrid history
  are absent. Pressure is prepared and validated on each level even when that
  level has no SST/TSK arrays, so an aloft descendant can inherit populated
  parent pressure; if no finite positive parent or diagnosed ground pressure
  is available, migration aborts. Markerless Metgrid arrays are already in
  kelvin and are restored directly.

A checkpoint whose ``job_info`` does not record the ``erf.init_type`` that
wrote these arrays is rejected rather than guessed, unless the restart asserts
it with ``erf.legacy_surface_temperature_init_type``, which stands in for
missing provenance only and cannot override provenance the checkpoint records.

The land mask (``LMASK_0``, or ``erf.is_land`` where the checkpoint holds
none) is restored first, as it decides which array the surface layer takes at
each cell: at a level with an SST, the TSK everywhere if the level has one,
else the SST over water and the prescribed ``most.pot_surf_temp`` over land;
at a level without an SST, neither. (A finer level without an SST of its own
counts as having one if a coarser level has one, as a restart gives it the
coarser level's.) The restart tells three kinds of cells apart:

- **Unused placeholders** -- a value that is not a finite, positive
  temperature where the surface layer does not take the array, such as the
  zero SST WRF writes over land -- are kept as they are, not counted as
  failures and not made to look like a temperature; the restart reports how
  many it kept.
- **Invalid active values** -- a value that is not a finite, positive
  potential temperature where the surface layer takes the array -- abort the
  restart, naming the first such cell.
- Every positive value is converted, used or not, so that the array holds
  absolute temperature throughout and stays consistent if a later regrid or
  input makes the cell active; one without a finite, positive conversion
  pressure aborts the restart.

Unknown provenance, a WRFInput pressure source that cannot be established, or
a required input pressure that cannot be read causes an abort with a message
naming the limitation rather than silently selecting an unrelated pressure.

In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

The computational cost associated with reading and writing checkpoint files is
typically negligible relative to the overall cost of the simulation;  in a recent performance
study the cost of writing a checkpoint file was roughly a percent of the cost of a single timestep.

Writing the Checkpoint "Files"
==============================

.. _list-of-parameters-8:

List of Parameters
------------------

+---------------------------------+----------------+----------------+----------------+
| Parameter                       | Definition     | Acceptable     | Default        |
|                                 |                | Values         |                |
+=================================+================+================+================+
| **erf.check_file**              | prefix for     | String         | “*chk*”        |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **erf.check_int**               | how often (by  | Integer        | -1             |
|                                 | level-0 time   | :math:`> 0`    |                |
|                                 | steps) to      |                |                |
|                                 | write restart  |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **erf.check_per**               | how often in   | Real           | -1.0           |
|                                 | simulation     | :math:`> 0`    |                |
|                                 | time to write  |                |                |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+

Restarting
==========

.. list-table:: Restart parameters
   :header-rows: 1
   :widths: 30 45 15 20

   * - Parameter
     - Definition
     - Acceptable values
     - Default
   * - ``erf.restart`` or ``amr.restart``
     - Checkpoint directory from which to restart. ``amr.restart`` takes precedence if both are set.
     - String
     - Not used if unset
   * - ``erf.legacy_surface_temperature_init_type``
     - The ``erf.init_type`` of the run that wrote a markerless or version-1
       checkpoint with ``SST_0``/``TSK_0`` whose ``job_info`` does not record
       it. Markerless ``WRFInput`` arrays are legacy potential temperature and
       use the best-available WRFInput pressure policy above (fine-level
       recovery is approximate); markerless ``Metgrid`` arrays hold kelvin and
       are restored directly. Version-1 arrays are potential temperature for
       either source: WRFInput uses its best-available pressure policy, and
       Metgrid uses restarted ground pressure or parent ground pressure for
       aloft-only coverage. Where ``job_info`` records the init type, that
       takes precedence and a conflicting assertion is ignored with a
       warning. The assertion must name the actual old source and is not a
       general override of checkpoint metadata.
     - ``WRFInput`` or ``Metgrid`` (case-insensitive)
     - Unset

.. _examples-of-usage-7:

Examples of Usage
-----------------

-  **erf.check_file** = *chk_run*

-  **erf.check_int** = 10

   means that restart files (really directories) starting with the
   prefix “*chk_run*” will be generated every 10 level-0 time steps.
   The directory names will be *chk_run00000*, *chk_run00010*,
   *chk_run00020*, etc.

-  **erf.check_per** = 5.0

   means that restart files (really directories) starting with the
   prefix “*chk_run*” will be generated whenever the simulation time
   passes a multiple of 5.0.  The directory names will reflect the
   integer number of steps which have elapsed.

To restart from *chk_run00061*, for example, set **erf.restart** to that
directory:

-  **erf.restart** = *chk_run00061*

ERF also accepts **amr.restart** for compatibility with AMReX and existing
regression inputs. When both keys are present, **amr.restart** takes
precedence because ``ReadParameters()`` queries it after ``erf.restart``.

SBM checkpoints
---------------

When ``erf.moisture_model = SBM``, the liquid spectral distribution is
authoritative and is checkpointed separately from ERF's conventional compact
moisture state.

An SBM checkpoint therefore contains two additional pieces of information:

* ``SBM_Schema``, which records the exact spectral layout and current
  representation contract; and
* ``SBMSpectrum`` at each active level, which contains the authoritative
  bin-resolved liquid state.

The saved schema must match the active SBM layout exactly on restart. A
different moment mode, spectral grid, pivot definition, cloud/rain split, or
other schema-defining change is rejected rather than interpreted as compatible.

ERF validates the saved spectral state before using it to reconstruct the
compact cloud- and rain-water fields. Spectral values must be finite and must
satisfy the constraints of the saved representation. For two-moment bins this
includes requiring the mean particle mass implied by the saved mass and number
to lie within that bin.

The checkpoint also contains the conventional ``qc`` and ``qr`` fields in the
normal ERF state. Before overwriting those fields, ERF compares them with the
cloud- and rain-water projection calculated from the saved authoritative
spectrum. The comparison uses a tight floating-point consistency tolerance.

Only after the spectral state and the saved bulk projection have passed these
checks are ``qc`` and ``qr`` regenerated from the spectrum. This ordering is
intentional: a corrupted compact field or corrupted spectrum cannot be hidden
by projecting a new value over the evidence of the inconsistency.

A restart fails if the required SBM schema or spectrum is missing, the schema
does not match, the authoritative spectral state is inadmissible, or the saved
bulk liquid fields are inconsistent with the spectrum.

See :ref:`sec:SpectralBinMicrophysics` for the current SBM state definition and
supported execution envelope.

Time-accumulated diagnostics
============================

Some diagnostics are accumulated over the course of a run rather than computed
from the instantaneous state, so they carry state of their own that must be
carried across a restart. ERF checkpoints that state, and a run restarted from
a checkpoint reproduces the equivalent uninterrupted run:

-  the running sum and normalizer behind the time-averaged velocity fields
   (``erf.time_avg_vel``, plotted as ``u_t_avg``, ``v_t_avg``, ``w_t_avg`` and
   ``umag_t_avg``), so the averaging window continues to extend back to the
   start of the original run;

-  the exponentially filtered surface-layer averages used when
   ``erf.most.time_average`` is set (see :ref:`sec:surface_layer`), so
   :math:`u_{*}`, :math:`\theta_{*}` and the Obukhov length after restart are
   computed from the same filter history the uninterrupted run would have had;

-  the accumulated surface precipitation reported by the microphysics schemes.

Two caveats are worth knowing:

-  A checkpoint written by a version of ERF that predates this state being
   saved does not carry it. Restarting from such a checkpoint simply starts the
   average over, which is what those checkpoints have always done; ERF prints a
   note when it detects this.

-  Regridding rebuilds the per-level containers that hold these averages, so a
   regrid restarts the velocity time average and the surface-layer time filter
   on the affected level. With AMR the averaging window is therefore "since the
   last regrid on that level," not "since the start of the run."
