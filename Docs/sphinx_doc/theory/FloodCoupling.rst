.. role:: cpp(code)
   :language: c++

.. _sec:FloodCoupling:

Flood Coupling
==============

Rain that falls on the land surface runs off, ponds and floods. ERF provides
the atmosphere side of a flood model: where the flood grid sits and what the
flood is handed each step. The flood itself, the model that routes the water
over the land, is a *flood engine* that plugs into this coupling. The coupling
is one way: the engine reads the atmosphere and writes nothing back.

The layout copies the fire coupling (:ref:`sec:FireCoupling`) under its own
names, so that a flood engine does not depend on the fire code. Each engine
lives in its own repository, which holds ERF as a git submodule pinned to an
ERF release or commit and builds against it; ERF needs no change to add an
engine. ERF itself carries one small example engine and its tests; ERF's CI
does not build the flood coupling, so build with ``ERF_ENABLE_FLOOD=ON`` to run
them.

Building with the flood coupling
--------------------------------

The coupling layer is compiled with the CMake option ``ERF_ENABLE_FLOOD=ON``
(``USE_FLOOD = TRUE`` with GNU make); it is off by default. It adds nothing to
the state. The compile definition is ``ERF_USE_FLOOD``.

:cpp:`erf.flood.enable = true` needs an engine compiled into the executable;
without one the run stops at start-up and says how to add one. A build
without the flood coupling stops on :cpp:`erf.flood.enable = true` as well.

An engine is built into ERF in one of two ways, as for a fire engine.

* **ERF as a library.** The engine's repository has its own CMake project,
  which adds the ERF submodule, compiles the engine into its own executable
  and sets that executable up with ``erf_setup_executable`` (see
  :ref:`sec:ERFAsLibrary`). The engine's :cpp:`ERF_REGISTER_FLOOD_ENGINE` line
  runs from that executable, so ERF is given no engine directory.
* **Engine directory.** ERF's own build is given the engine directory:
  ``ERF_FLOOD_ENGINE_DIRS`` with CMake (a list; each directory holds an
  ``ERFFloodEngine.cmake`` that defines ``erf_flood_engine_add(target)``, and
  optionally ``erf_flood_engine_add_unit_tests(target)``), or ``FLOOD_HOME``
  with GNU make (each directory holds a ``Make.package``).

.. code-block:: console

   cmake -DERF_ENABLE_FLOOD=ON \
         -DERF_FLOOD_ENGINE_DIRS=/path/to/my_flood_engine ..

   make USE_FLOOD=TRUE FLOOD_HOME=/path/to/my_flood_engine

``Tests/FloodEngineExample`` is a complete engine directory that builds both
ways; its ``README.md`` describes the example run.

Adding a flood engine
---------------------

An engine derives from :cpp:`FloodEngine` (``Source/FloodCoupling/ERF_FloodEngine.H``)
and registers itself under the name the deck gives in ``erf.flood.engine``:

.. code-block:: c++

   #include <ERF_FloodEngineRegistry.H>

   class MyFloodEngine : public FloodEngine { /* ... */ };

   ERF_REGISTER_FLOOD_ENGINE(MyFloodEngine, "my_engine")

``Source/FloodCoupling/ERF_FloodEngineTemplate.H`` lists every method with
what it must do, and compiles as an engine that does nothing. The required
methods are :cpp:`name`, :cpp:`verify_prerequisites`, :cpp:`initialize` and
:cpp:`advance`. The others have defaults: :cpp:`start_run`,
:cpp:`wants_surface_precip` (true), :cpp:`write_checkpoint` (nothing is
saved), :cpp:`read_checkpoint` (a restart from a checkpoint with flood state
stops, naming the engine), :cpp:`write_plotfile` (no plotfile),
:cpp:`plot_var_names` and :cpp:`plot_var` (no fields in the ERF plotfile).

The flood grid
--------------

The flood grid is a 2D grid that refines the horizontal cells of one
atmospheric level, ``erf.flood.anchor_level`` (the finest level at start-up by
default), by ``erf.flood.grid_ratio`` (:math:`C`, default 1). Flood cell
:math:`(i, j)` lies in atmospheric cell
:math:`(\lfloor i/C \rfloor + i_0, \lfloor j/C \rfloor + j_0)`, where
:math:`(i_0, j_0)` is the lower corner of the region the level's boxes cover
(zero on level 0). The grid keeps the level's box decomposition and MPI owners,
so the maps between the two grids need no communication. The helpers in
``ERF_FloodGrid.H`` copy an atmospheric surface field to every flood cell under
it and average flood cells back to their column.

The level must satisfy, checked at start-up with a message that names the
input to change:

* every box spans the level's full height (the default
  ``amr.no_box_split_dir = 2`` gives this on level 0; on a refined level the
  refinement box must reach the domain top);
* the boxes cover one rectangle;
* above level 0, the level keeps its boxes for the whole run
  (``erf.regrid_int = -1``), since the flood grid is built once.

A refinement level over part of the flood grid is allowed: the flood then
reads the precipitation and wind of its own level, and a note is printed.

What the engine is handed
-------------------------

Each step of the flood grid's level, after the dycore and the microphysics,
:cpp:`advance` receives a :cpp:`FloodForcing`:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Member
     - Contents
   * - ``time``, ``dt``
     - Start and length of the step [s]
   * - ``xvel``, ``yvel``
     - Face velocities of the level after the step [m/s]
   * - ``z_phys_cc``
     - Cell-center heights of the terrain-following mesh [m]
   * - ``precip_k0``
     - The step's surface precipitation [kg m\ :sup:`-2`, equal to mm of water]
       on the :math:`k = 0` slab of each box: component 0 liquid, component 1
       frozen (see below for what each scheme counts as frozen). ``nullptr``
       when :cpp:`wants_surface_precip` returns false.

The precipitation is the change of the microphysics surface accumulators over
the step,

.. math::

   P^{n} = A(t^{n+1}) - A(t^{n}),

where :math:`A` is the accumulation of ERF's own surface precipitation
diagnostics, read on the flood grid's level after the step. The sum of
:math:`P^n` over the steps of a run is therefore the change of :math:`A` over
the run, to round-off. For a scheme that keeps a total and a frozen
accumulator, the liquid part is the total less the frozen part and can be
slightly negative from round-off. An engine that bounds it should carry the
negative part over to the next step rather than drop it, so that its total
still equals :math:`A`.

When a finer level covers part of the flood grid's level, ERF averages the
finer level's accumulators down onto it once the finer level has caught up
(and, with super-droplet microphysics in two-way mode, also at each plotfile
and regrid). The flood takes that correction, of either sign, in its next
step, so under the finer level it lags the averaged accumulation by one
correction; the run total is unchanged.

.. warning::

   In a single-precision build the accumulators are single precision. As
   :math:`A` grows over a long or wet run, each step's :math:`P^n` is a
   shrinking fraction of it and keeps fewer digits; at 0.1 mm/h on 0.5 s
   steps over 100 mm already accumulated, about 10 % of the rain is lost
   to round-off. Use a double-precision build for long flood runs.

The two components follow ERF's surface precipitation diagnostics. Kessler
and the super-droplet model hand rain only. SAM hands rain as liquid and snow
plus graupel as frozen; its cloud-ice sedimentation does not reach the surface
accumulators. Morrison and WSM6 keep a total and a frozen accumulator, with
snow and graupel (and hail) as frozen; the cloud-water fallout and, for
Morrison, the cloud-ice fallout are in the total and so are counted as liquid.
An engine that reads the precipitation needs a moisture model that accumulates
it (Kessler, SAM, Morrison, WSM6, WDM6 with precipitation); with any other the
run stops at start-up.

Checkpoints and restarts
------------------------

The coupling writes a small marker ``Level_<lev>/FloodGrid`` to the checkpoint
with the flood grid's level, :math:`C` and the columns it covers; the engine
writes its own state next to it. On a restart the marker must match the run's
flood grid, or the run stops naming ``erf.flood.anchor_level`` and
``erf.flood.grid_ratio``. A checkpoint without a marker (from a run without
the flood) starts the engine from its initial state and prints a note, so a
flood can start from an atmosphere spun up without it. The precipitation the
flood last read is saved next to the marker (``FloodPrecipRef``), so a pending
correction from a finer level is kept and a restart reproduces a straight run;
without a marker the reference is taken from the restored accumulators. A restart that regrids level 0 is refused.

Output
------

The engine's own plotfile is written with the prefix ``erf.flood.plot_file``
on the schedule ``erf.flood.plot_int`` / ``erf.flood.plot_per``. With
neither given, the engine writes its plotfile at every 3D plotfile
(``erf.plot_int_1`` / ``erf.plot_per_1``); an engine whose plotfile is large
may want its own, sparser schedule. As for the fire, no flood plotfile is added at the end of
a run that stops between two plot times. A field the engine names in :cpp:`plot_var_names` (each
must start with ``flood_``) can also be listed in ``erf.plot_vars_1``; it is
written on the :math:`k = 0` cells of the flood grid's level as the mean over
the flood cells of each column, and zero elsewhere.

Inputs
------

The coupling's inputs are listed under :ref:`Inputs <sec:InputsFlood>`; an
engine's own inputs are its own.
