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
