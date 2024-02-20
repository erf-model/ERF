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
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

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

+---------------------------------+----------------+----------------+----------------+
| Parameter                       | Definition     | Acceptable     | Default        |
|                                 |                | Values         |                |
+=================================+================+================+================+
| **erf.restart**                 | name of the    | String         | not used if    |
|                                 | file           |                | not set        |
|                                 | (directory)    |                |                |
|                                 | from which to  |                |                |
|                                 | restart        |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+

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

To restart from *chk_run00061*,for example, then set

-  **amr.restart** = *chk_run00061*

