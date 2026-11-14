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

- **WRFInput: exact recovery.** The arrays are converted back with the
  wrfinput surface pressure ``PSFC`` they were converted with, read again from
  where each level read its SST and TSK: level 0 from the first file of
  ``erf.nc_init_file_0``, the only one it reads; a finer level from its own
  file, or interpolated from the coarser level where the level took its SST
  and TSK from that level (no init file of its own, or a level-0
  ``erf.nc_low_file`` and none of its own). This needs a restart with
  ``erf.init_type = WRFInput``, in a build with NetCDF, that names the
  wrfinput files of the original run; any other restart aborts. A finer level
  with several init files of its own reads the first alone until another
  starts: the run reads a later file, and regrids the level for it, only once
  the start time and elapsed time of the coarser level reach the file's start
  time. Its arrays are recovered with the ``PSFC`` of its first file if no
  later file had started by the time of the checkpoint; once one had, the
  checkpoint does not record which file each of its grids took them from, or
  whether a regrid replaced them, and the restart aborts.
- **Metgrid: approximate migration.** The surface pressure of the met_em
  input is not read again, so the arrays are converted with the surface
  pressure of the restarted state, and the restart warns that it does so.
  This preserves the surface potential temperature the run had at the
  restart; it does not recover the SST of the input. The migrated absolute
  SST and TSK keep the factor :math:`(p_{restart}/p_{input})^{R_d/c_p}` of any
  change in surface pressure since the input was read, the pressure-related
  bias of the old potential-temperature storage. To avoid it, restart from a
  version-2 checkpoint, or start again from the input.

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

Unsupported configurations -- unknown provenance, a pressure the restart
cannot read again (a finer level with several init files of its own, one of
which after the first had started by the time of the checkpoint) -- abort
with a message naming them rather than converting with another pressure.

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
     - The ``erf.init_type`` of the run that wrote a markerless or version-1 checkpoint with ``SST_0``/``TSK_0`` whose ``job_info`` does not record it. For a markerless checkpoint, ``WRFInput`` means the arrays are potential temperature, recovered exactly with the wrfinput ``PSFC``, and ``Metgrid`` that they are absolute temperature, restored directly. For a version-1 checkpoint, whose arrays are potential temperature, ``WRFInput`` recovers them exactly and ``Metgrid`` migrates them approximately with the surface pressure of the restarted state. Where the ``job_info`` records the init type, that takes precedence, and a conflicting assertion is ignored with a warning. It must name the actual old source and is not a general override of checkpoint metadata.
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
