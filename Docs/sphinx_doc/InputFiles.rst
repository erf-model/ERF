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

Creating an inputs file for a new case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Exec/Tools/erfcreateinput.py`` writes a complete, commented inputs file from a few
choices. The chosen option of each group (equation set, turbulence closure, lower
boundary, forcing, outputs) is active with its sub-keys; the other options of the group
are commented lines beside it, so a choice is changed by moving one ``#``. Values the
tool cannot know are marked ``USER:``.

::

    python3 Exec/Tools/erfcreateinput.py --list
    python3 Exec/Tools/erfcreateinput.py --preset abl_convective_les --outdir my_cbl
    python3 Exec/Tools/erfcreateinput.py --setup abl --turbulence deardorff --surface most \
        --roughness file --coriolis --forcing geostrophic --output plt,chk,profiles --outdir my_abl

A case is described by a preset (``--preset``, listed by ``--list``), by single options, or by a
case file (``--case-file``); options on the command line replace the values of the other two.
Each case directory also holds ``case.json``, the options that reproduce it.

Files the deck names are written as dummies where a text layout exists: the input sounding,
the roughness table, the surface time table, the geostrophic wind table, the terrain and
building tables, the forest patch file and the wind-farm tables. Each dummy is valid in
format and sized to the grid, so the run starts, with placeholder values. Files that cannot
be generated (``wrfinput`` and ``wrfbdy``, ``met_em``, the RRTMGP coefficients, the Noah-MP
namelist and table) are listed as files to supply. The script prints, and saves in
``FILES_TO_REPLACE.txt``, the layout of each file and what must replace it.

A combination that ERF refuses at start-up is refused by the script, with a message that
names the options (for example a PBL scheme without ``--surface most``, or
``--forcing geostrophic`` without ``--coriolis``). The options cover idealised cases
(``--setup abl``, ``hill``, ``bubble``, ``density_current``) and real-data cases
(``wrfinput``, ``metgrid``); terrain from an analytic hill or a table on a fitted mesh, by
immersed forcing or as an embedded boundary; buildings; vertical stretching; the turbulence
closures; moisture, radiation and land-surface models; forest canopy; wind farms; mesh
refinement; and the output and sampling keys.
``Exec/Tools/erf_input_creator/README.md`` lists every option, the layout of each dummy
file, what was run in ERF, and how to add an option.

The script complements the checker described next: it writes a new case, and
``check_erf_inputs.py`` then checks that case and its files, also after the dummy files have
been replaced with real data. On a generated case the checker also reports each dummy file
that still holds its placeholder values, each file to supply that is missing, and the values
marked ``USER:``.

With ``--amr``, ``--refine field`` tags cells by field tests instead of using static boxes:
``--refine-on theta:less:299 w:greater:2`` writes one condition per test under
``erf.<indicator>.conditions``, ``--refine-combine`` chooses ``and`` or ``or``, and
``--refine-box`` restricts the tests to a region (see :ref:`MeshRefinement`).

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
