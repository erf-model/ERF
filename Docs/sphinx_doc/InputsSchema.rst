Inputs Schema
=============

``Exec/Tools/build_erf_schema.py`` writes a per-build schema of ERF's inputs. It lists every literal ``ParmParse``
query in ERF and in the AMReX it builds against, typed from the variable each query fills. For each key it records:

- whether a build feature is needed (a key read only under ``#ifdef ERF_USE_X``);
- the variable's default, as a parsed literal;
- whether ERF now refuses the key at start-up.

Each fact names the source line it comes from. ``Exec/Tools/check_erf_inputs.py --schema`` checks a deck against
it.

The schema is reference material for checking and tooling, not a fixed workflow contract. It reads only literal
queries, and a key whose variable type it cannot find is ``unknown``. Prefixes ERF builds at run time (boundary
faces, stations, refinement indicators, particle containers) are read as families, listed in ``scan.FAMILIES``;
reads the generator still cannot place are listed under ``unresolved``.

Generating a Schema
-------------------

The tool reads source with ``git archive``, so nothing is checked out and no build is run. AMReX is taken at the
commit ERF's ``Submodules/AMReX`` points to:

.. code-block:: bash

   python3 Exec/Tools/build_erf_schema.py --repo . --erf HEAD \
       --cmake="-DERF_ENABLE_CUDA=ON -DERF_ENABLE_NETCDF=ON -DERF_ENABLE_RRTMGP=ON" --out erf_inputs_schema.json

Build options come from ``--cmake=...``, from ``--cmake-cache CMakeCache.txt``, or from both. A ``CMakeCache.txt``
holds every ERF option. A ``-D`` list holds only the options passed, and the rest take the ``CMakeLists.txt``
defaults. The schema records which kind it was given.

Facts that belong to the machine rather than to the build come from ``--site-facts``, a JSON list. Whether the MPI
is GPU-aware is one; it affects ``amrex.use_gpu_aware_mpi = true``.

The same commits, options and site facts give the same bytes.

What a Key Records
------------------

- **Type:** ``int``, ``real``, ``bool``, ``string`` or ``enum``, each also as a list; ``unknown``; or ``removed``.
  When one variable name has several types in a file, the read call decides (a list for ``queryarr``/``getarr``, a
  scalar otherwise), then the nearest declaration before the read. Every read method of ``AMReX_ParmParse.H`` is
  read; a ``queryline`` read takes the whole line, a ``querytable`` read is recorded as a table, a
  ``query_enum_sloppy`` read records the characters it ignores, and a ``WithParser`` or ``AsDouble`` read is marked
  ``parser``, since it takes an expression by design.
- **Layer:** whether AMReX or ERF reads the key. A key both read keeps both source lines. When two reads disagree, an
  array read makes the key a list, an enum read wins over a ``std::string`` probe, and the pair is kept in
  ``type_conflict``. A refusal in problem-specific code (``Source/Prob/``) is listed under ``refused_by`` and does
  not override a read elsewhere.
- **Default:** the initializer of the variable the query fills, in canonical form, or none. A ``pp.add`` of a
  literal in ERF's ``add_par`` (``Source/main.cpp``) is the default instead: it runs before the inputs are read, so
  the inputs override it (``amr.n_proper`` is 2, not AMReX's 1).
- **Build requirements:** the ``ERF_USE_X`` macro a key is read under, with the ``#ifdef`` line. Also any enum value
  that is refused, or whose model is only built, without a feature.
- **Removed:** the key's read sits in an ``if`` whose body always calls ``Abort`` or ``Error`` and does not go on to
  read the key. This covers ``contains(k)``, ``query(k, x) > 0`` and ``countval(k) > 0``, the last also through a
  variable.

Each value type has one canonical form: ``1``, ``1.0`` and ``1e0`` are one real; ``t`` and ``1`` are ``true``; an
enum name takes its declared case. The schema lists these rules under ``canonical_rules``.

Checking a Deck
---------------

.. code-block:: bash

   python3 Exec/Tools/check_erf_inputs.py --schema erf_inputs_schema.json inputs

With a schema, the checker reports:

- a key no query reads, as WARNING; when the schema's ``unresolved`` list has a read of that name through a
  prefix the generator could not resolve, as INFO naming the line;
- a removed key, as ERROR;
- a string or enum value ERF does not know, or fewer per-level values than levels, as ERROR;
- a bool, int or real value that is not a plain literal (``2*0.5``, ``300.;``, ``1e9`` for an int), as WARNING,
  since ParmParse hands it to its parser at run time, except for a key read with the parser (``queryWithParser``,
  ``queryAsDouble``); ``1D0`` is an ERROR, since neither ParmParse nor its parser
  reads a ``d`` or ``D`` exponent;
- more words than the key reads (``erf.check_int = 10 steps``; ParmParse reads the first), as WARNING;
- a key or value that needs a feature the build lacks, as ERROR;
- a ``@NAME@`` placeholder (a ``configure_file`` template), as INFO.

It also compares the schema's removed keys with the curated list of retired keys in ``checks/keys.py``, and prints
any difference as INFO.

Tests
-----

.. code-block:: bash

   python3 -m unittest discover -s Exec/Tools/erf_inputs_schema/tests
   python3 Exec/Tools/check_erf_inputs.py --self-test
