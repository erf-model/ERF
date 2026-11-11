.. _sec:running:

Running
-------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to be specified along with the executable as an ``argv`` option, for example:


::

    mpirun -np 64 ./ERF3d.xxx.yyy.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::

    mpirun -np 64 ./ERF3d.gnu.DEBUG.MPI.ex inputs amr.restart=chk0030 erf.use_gravity=true

See :ref:`sec:Inputs` for details on run-time options that can be specified. If running on a Mac and getting
errors like ``SIGILL Invalid, privileged, or ill-formed instruction``, see the note on that page about
runtime error-checking options.

Checking an inputs file before running
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ParmParse`` does not warn about a key that nothing reads, and several of ERF's own
start-up checks only run once the input data have been read. The script
``Exec/Tools/check_erf_inputs.py`` parses an inputs file the way ``ParmParse`` does
(comments, ``\`` continuations, ``FILE =`` includes, the last definition of a key wins)
and reports, without running ERF:

* **ERROR**: the run would stop at start-up (an AMReX or ERF abort), for example a
  non-periodic face without a ``type``, ``amr.n_cell`` not divisible by
  ``amr.blocking_factor``, an odd ``erf.fixed_mri_dt_ratio``, or a named file that is
  missing;
* **WARNING**: the run starts but a setting is very likely not what was meant: a key
  the selected mode never reads (``erf.Cs`` without an LES, acoustic substep keys in an
  anelastic deck, ``erf.most.*`` without ``zlo.type = surface_layer``), a misspelt key
  with its closest documented name, a fixed time step whose acoustic or advective CFL
  exceeds one for the wind speed given with ``--umax``;
* **ADVICE**: a working setting that could be chosen better, for example ``erf.cfl``
  instead of ``erf.fixed_dt`` for real-data or terrain cases, an odd or uneven
  ``amr.max_grid_size``, or a box count that does not divide the rank count given with
  ``--ranks``;
* **INFO**: the cell sizes, the estimated time step and substep count, the box layout.

When the deck names files that are present (``wrfinput``, ``wrfbdy``, ``wrflow``,
``met_em``, terrain tables, soundings, wind tables), their contents are compared with the
deck: grid dimensions and extent, start date, boundary width and time span, terrain
coverage of the domain, sounding surface temperature and top.

::

    python3 Exec/Tools/check_erf_inputs.py inputs
    python3 Exec/Tools/check_erf_inputs.py --ranks 64 --umax 30 --min-severity WARNING inputs
    python3 Exec/Tools/check_erf_inputs.py --self-test

The exit code is 1 when a deck has an ERROR (``--strict`` adds WARNING). Run it from
inside the repository, where it finds the documentation tables and the ``Source`` tree
used by the misspelt-key check. ``Exec/Tools/erf_input_checker/README.md`` lists every
check and the options. A quick alternative after a run is the "Unused ParmParse Variables"
table AMReX prints at exit with ``amr.v = 1``, which lists every key nothing read.
