.. _sec:StructureCoupling:

Structure coupling
==================

The structure coupling (``Source/MovingBodies/Coupling``) connects ERF to a
model of structures in the flow: bodies that the air loads and that push back on
the air, such as power-line conductors and their towers. ERF contains no
structure model. A model lives in its own directory outside ERF, registers
itself by name, and plugs into a small interface, ``StructureModel``
(``ERF_StructureModel.H``), in the way fire engines plug into the fire coupling
(:ref:`sec:FireCoupling`). The coupling is built with the actuator core
(:ref:`sec:ActuatorCore`), whose velocity sampler and force spreader a model is
expected to use: ``-DERF_ENABLE_MOVING_BODIES=ON`` (CMake) or
``USE_MOVING_BODIES = TRUE`` (GNU make).

What ERF does and what the model does
-------------------------------------

ERF calls the model at fixed points of a run and hands it the flow; the model
owns everything else: its own inputs, its points, its dynamics, its diagnostics
and its checkpoint data.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ERF call
     - What happens
   * - creation
     - In the ERF constructor, ``erf.structures.model`` selects a registered
       model, which reads its own inputs under its own prefix. A name that no
       compiled-in model has stops the run, listing the names that are.
   * - ``level()``
     - The model names the AMR level it couples to; ERF calls ``advance()`` and
       ``add_momentum_sources()`` on that level only. It must be the finest
       level that covers the model's support. On a refined level the model's
       points and kernels must stay inside the level's boxes, also after a
       regrid. With ``erf.coupling_type = TwoWay``, coarser levels step with the
       model's source restricted onto them (:ref:`sec:StructureCoarseSource`),
       and averaging down then replaces the covered cells. A finer level would
       never see it: averaging down would overwrite the forced cells with that
       level's unforced data.
   * - ``initialize()``
     - Once, after the initial data (or the checkpoint) exist and before the
       step-0 plotfile, and not again after a regrid (the model reads the mesh
       from ``StructureStep`` each step). ``StructureInit`` holds the level's
       geometry, the x-velocity, the mesh heights and Jacobians on a stretched or
       terrain-following mesh, whether every map factor is 1 (the actuator core
       assumes so), the turbulence closure, the checkpoint name on a restart,
       on request the immersed terrain surface, and what a model that drives a
       solver of its own needs before the first step: the level's nominal step
       (``fixed_dt``, not positive without ``erf.fixed_dt``), the time at which
       the run stops (``stop_time`` on the clock of ``time``: the input
       ``stop_time``, or ``stop_datetime`` less ``start_datetime``), ``max_step``,
       and whether any level of the run is ``FullyCompressible``
       (``fully_compressible``, the solver type; a ``RK3_NoSubstepping`` level
       counts). The step taken is ``StructureStep::dt``: it differs from
       ``fixed_dt`` on the first steps with ``erf.init_shrink`` below 1 or after
       a restart, and on a last step that ERF fits to the stop time. With
       ``start_datetime`` the stop time can differ from the input by about
       :math:`10^{-7}` s, so a model compares it with a tolerance.
   * - ``advance()``
     - Once per step of the model's level, with the state at the start of the
       step (``StructureStep``), before the step's Runge-Kutta stages. The model
       samples the wind, advances its structures and sets the momentum source.
       The step ``StructureStep::dt`` is a ``double``, as ERF keeps it, also in
       a single-precision build.
       ``StructureStep`` also says whether every map factor of the level is 1
       now, so a level a regrid made is checked as well, and holds the fire view
       (``StructureStep::fire``, below) when the run has a fire.
   * - ``support()``
     - Right after each ``advance()``: boxes [m] that hold the model's points
       and their kernels' reach. When a level finer than the model's exists and
       its boxes overlap a support box, ERF stops the run and asks for the
       finer level (``structure.example.level`` for the example). On a
       stretched or terrain-following mesh the check uses the finer level's
       own node heights: for each node level, the lowest and the highest
       height over the level's boxes (``erf_actuator::mesh_z_levels``). ERF
       builds this table when the finer level's boxes change (and at every
       step with a moving terrain), not at every step. A
       finer box counts where it lies over a support box's footprint and
       the level holds cells at the support's heights, so the check errs on
       the side of stopping where the terrain is steep, and it holds
       whether or not the finer boxes reach the ground. A support box needs
       its low end at or below its high end and no NaN; an end may be
       infinite.
   * - ``add_momentum_sources()``
     - At every Runge-Kutta stage of the model's level: the model adds its force
       on the air [N/m\ :sup:`3`] to the face-centred momentum sources. The source
       is held fixed over the step. ERF applies no mask to it: the model keeps
       its kernels clear of embedded boundaries, immersed terrain and buildings.
   * - ``write_checkpoint()``
     - Each time ERF writes a checkpoint, under ``<checkpoint>/<model name>``.
       A model that restores its state in ``initialize()`` returns true from
       ``supports_restart()``; otherwise a restart stops the run.
   * - ``finalize()``
     - Once, on every rank, when the run ends normally, after the last
       plotfile and checkpoint: the model flushes and closes its own files or
       shuts down a solver it drives. Not called when the run aborts. The
       default does nothing.
   * - ``plot_var_names()``, ``fill_plot_var()``
     - The model's fields that ``erf.plot_vars_*`` may list, filled into the
       ERF plotfile on every level. The names carry the model's prefix; an
       empty name, a name listed twice or a name ERF or another model already
       uses stops the run.

.. _sec:StructureCoarseSource:

The force on the coarser levels
-------------------------------

ERF steps the coarser levels first. With ``erf.coupling_type = TwoWay``, at the
start of each step of a level :math:`\ell` below the model's level :math:`L`,
ERF asks the model for the source of its last ``advance()`` on level :math:`L`,
restricts it onto level :math:`\ell` and adds it to that level's momentum
sources over the step (the restriction is by Ann Almgren, erf-model/ERF#5022;
``ERF_StructureCoarseSource.cpp``). The restriction keeps the force the source
puts on the air. The force on a face is

.. math::

   F_f = s_f \, J_f \, \Delta x \, \Delta y \, \Delta z ,
   \qquad
   J_f = \tfrac12 \left( J_{\mathrm{lo}} + J_{\mathrm{hi}} \right),

where :math:`s_f` is the source on the face [N/m\ :sup:`3`], :math:`\Delta x,
\Delta y, \Delta z` the level's cell sizes in index space, and
:math:`J_{\mathrm{lo}}, J_{\mathrm{hi}}` the mesh Jacobian ``detJ_cc`` of the two
cells beside the face (the one cell inside the domain on a non-periodic domain
face; :math:`J = 1` on a uniform-dz mesh). One level at a time, with refinement
ratio :math:`r` along the face normal, coarse face :math:`I` takes the force of
every fine face :math:`r I + m`, :math:`|m| < r`, along the normal with the
weight :math:`(r - |m|)/r`, and of every fine face across it. Each fine face
gives its force out with weights that add up to 1, so

.. math::

   \sum_{\mathrm{coarse}} F_f = \sum_{\mathrm{fine}} F_f

exactly, also when the fine level is split into boxes or crosses a periodic
boundary. On a non-periodic domain face (the ground, the top, an inflow or
outflow side) the boundary condition sets the momentum and a source there does
nothing, so a fine face inside the domain gives its whole force to the first
coarse face inside instead of a share :math:`(r - |m|)/r` to the domain face;
the force on the faces inside the domain is then kept too, for a kernel that
reaches the ground, such as a tower's. The coarse source is the coarse force divided by the coarse face's
:math:`J_f \, \Delta x \, \Delta y \, \Delta z`. On a uniform-dz mesh no Jacobian
is applied. On a terrain-following mesh :math:`J` varies from column to column,
and the restriction of :math:`s_f` alone, without :math:`J_f`, changes the
force: by :math:`2.6 \times 10^{-3}` of itself for the example model over a
60 m high Schar ridge with a 2-to-1 refined level and a vertical stretching of
1.08, and by :math:`2.9 \times 10^{-3}` with a stretching of 1.15.
The Jacobian is the one at the coarse step's start; with moving terrain it may
differ from the one of the model's last ``advance()`` by one coarse step.

Before the model's first ``advance()`` and after a regrid of its level, until
its next ``advance()``, the coarser levels step without the source. ERF writes
the model's level-:math:`L` source to the checkpoint
(``<chk>/StructureCoarseSource``: the level, then one MultiFab per face
direction). After a restart, until the model's first ``advance()``, ERF
restricts that source instead, so a restarted run steps the coarser levels as
the straight run does, on any number of ranks. A regrid of the model's level
before that ``advance()`` keeps it when the level covers the same region, also
in other boxes, and drops it, with a message, when the region changes: the
source would be missing where the level grew. A checkpoint without it (or
written for another level) restarts with the coarser levels stepping without
the source until the first ``advance()``.

Reading the fire
----------------

When ERF runs with a fire (``erf.fire.enable`` or ``erf.fire.engine``, :ref:`sec:FireCoupling`),
``StructureStep::fire`` points to a ``StructureFire``; without a fire it is
``nullptr``. A model can then heat its structures from the fire, for example.
The view holds:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Member
     - Meaning
   * - ``geom``
     - Geometry of the fire grid [m], one cell in z. Fire cell
       :math:`(i, j)` has its centre at the lower corner plus
       :math:`((i + \tfrac12)\,\Delta x,\ (j + \tfrac12)\,\Delta y)`, in the
       coordinates of the model's own points.
   * - ``lev``
     - The atmospheric level the fire grid refines (``erf.fire.anchor_level``).
       It need not be the model's level.
   * - ``time``
     - The time the fields belong to [s]: the time ERF assigns to the fire's
       level at the end of the fire's last step, or the time ERF started or
       restarted from before that step.
   * - ``heat_release``
     - Heat released [W/m\ :sup:`2`] on the fire grid, before its split into
       sensible and latent heat (``FireEngine::heat_release()``); ``nullptr``
       when the engine has none.
   * - ``level_set``
     - Level set on the fire grid, negative where burned
       (``FireEngine::level_set()``); ``nullptr`` when the engine has none.
   * - ``burned_area()``
     - Area burned so far [m\ :sup:`2`] (``FireEngine::burned_area()``),
       negative when the engine does not track it. It may reduce over the
       ranks, so every rank calls it or none does.

The fire steps after the dycore of its own level. A model on that level, or on
a coarser one, therefore sees the fire at its own step's time, exactly: ERF
gives the fire's fields the time it assigns to the level, not a sum of
substeps, which could differ by round-off. A model on a finer level sees the fire up to one step of the
fire's level ahead of its own step. The fields are spread over the ranks as the
fire grid is, so a model that needs the fire at its own points reduces over the
ranks itself. After a restart, until the fire's first step, the fields are what
the engine restored from the checkpoint; an engine that does not restore a field
holds its initial value there. Like the flow fields, the view and its fields
are valid only during ``advance()``: a model keeps no pointer to them.

Writing a model
---------------

``Tests/StructureModelExample`` is a complete model directory: one fixed point
that pushes on the air with a constant force, spread with the actuator core,
which samples the wind at the point and writes the impulse it has given the air
to the checkpoint. To start a new model:

#. Copy that directory, rename the class, and register it at namespace scope
   in its ``.cpp``: ``ERF_REGISTER_STRUCTURE_MODEL(MyModel, "my_model")``
   (``ERF_StructureModelRegistry.H``). The class is constructed as
   ``MyModel(int max_level)``.
#. Read the model's inputs under its own prefix, never under ``erf.``.
#. Give the directory an ``ERFStructureModel.cmake`` defining
   ``erf_structure_model_add(target)`` and a ``Make.package``, and build with
   ``-DERF_STRUCTURE_MODEL_DIRS=<dir>`` (CMake) or ``STRUCTURE_HOME = <dir>``
   (GNU make). Adding a model changes no file of ERF.
#. Optionally, define ``erf_structure_model_add_unit_tests(target)`` in the same
   ``ERFStructureModel.cmake`` to add the model's GoogleTest files to ERF's unit-test
   executable ``erf_unit_tests``, which has ERF's test main (CMake with
   ``ERF_ENABLE_UNIT_TESTS=ON``; there is no GNU make equivalent). The model is
   already compiled into that executable, so a test can create it with
   ``make_structure_model("my_model", max_level)``.
   ``Tests/StructureModelExample/ERF_GTestStructureModelExample.cpp`` is an example.
   With unit tests on, ERF includes ``ERFStructureModel.cmake`` a second time, from
   ``Tests/Unit``, so the file must only define functions: a top-level ``add_library``,
   ``add_subdirectory`` or ``option`` must sit behind a guard such as
   ``if(NOT TARGET ...)``, and ``include_guard`` would skip the second include and with it
   the model's unit tests.
#. Test a restart from a mid-run checkpoint against a straight run: the
   plotfiles must agree at zero tolerance.

Every MultiFab ERF hands a model is valid only during the call, in the default
arena (device memory on a GPU build); a model keeps no pointer to it.

The example model reads these inputs; its plotfile fields are zero until its
first ``advance()``, so in the step-0 plotfile and in the first plotfile written
on a restart.

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Input
     - Meaning
     - Range, default
   * - ``structure.example.position``
     - The point [m], three values
     - inside the domain, finite; required
   * - ``structure.example.force``
     - Its force on the air [N], three values
     - finite; required
   * - ``structure.example.epsilon``
     - Kernel width [m]
     - > 0, finite; required
   * - ``structure.example.level``
     - Level it couples to; -1 is the finest
     - -1 to ``amr.max_level``; -1
   * - ``structure.example.test_support``
     - What ``support()`` returns: ``reach`` (the kernel's reach), or, for
       ERF's own tests of its check, a box with a NaN (``nan``) or with its
       low end above its high end (``inverted``). ERF reads the box only when
       a finer level exists
     - ``reach``, ``nan`` or ``inverted``; ``reach``

Test
----

``Structures_ExampleModel`` runs the example model in a periodic box of still air
with no gravity, friction or Coriolis force, on two ranks, with a step of 0.25 s. A point at the box
centre pushes on the air with :math:`F_x = 50` N. With periodic sides and slip
walls, nothing else changes the air's x-momentum, so after :math:`t = 4` s

.. math::

   \int_V \rho\, u \, dV = F_x\, t ,
   \qquad
   \int_V u \, dV = \frac{F_x\, t}{\rho_0} = 200 \ \mathrm{m^4/s},

where :math:`V` is the box, :math:`u` the x-velocity, :math:`\rho` the density
and :math:`\rho_0 = 1` kg/m\ :sup:`3` its initial value, which the weak
pressure waves change by at most :math:`6 \times 10^{-9}`. The test reads the
volume integral of ``x_velocity`` from the plotfile (``fvolumesum``, which prints
6 significant digits) and checks it against 200 m\ :sup:`4`/s to a relative
:math:`10^{-5}`, then restarts from
the mid-run checkpoint and compares the final plotfiles at zero tolerance.

Two more runs check the finer-level rule with a refined level 1 and the model on
level 0: ``Structures_Abort_FinerLevelCovers`` refines the box around the point
and must stop; ``Structures_ExampleModel_FinerLevelElsewhere`` refines only the
lowest 50 m, below the kernel's reach, and must run. Three runs repeat the rule
on a stretched mesh with a level-1 box aloft, which has no box at the ground:
``Structures_Abort_FinerLevelAloft`` puts the box inside the kernel's reach and
must stop; ``Structures_ExampleModel_FinerLevelAloftBeside`` puts it beside the
reach and ``Structures_ExampleModel_FinerLevelAloftAbove`` above it, in the
point's own columns, and both must run. ``Structures_Abort_SupportNaN`` and
``Structures_Abort_SupportInverted`` make the example return a support box
with a NaN and with its low end above its high end
(``structure.example.test_support``), and the run must stop with the model's
name; ``Structures_Abort_ExampleTestSupport`` checks that input's allowed
values. ``Structures_ExampleModel_CoarseSourceRestart`` puts the model on a
level-1 box around the point, so level 0 steps with the restricted source, and
restarts from the mid-run checkpoint on three ranks instead of two: the final
plotfiles must agree at zero tolerance on both levels. A second restart
regrids level 1 into boxes of :math:`4^3` cells before the model's first step
(``erf.regrid_int = 1``, ``amr.max_grid_size = 4``), the same region in other
boxes: it must keep the source and agree to a relative :math:`10^{-12}`
(dropping the source misses by :math:`7 \times 10^{-4}`).
``Structures_ExampleModel_CoarseSourceRestartTerrain`` repeats it over a
ridge, to a relative :math:`10^{-12}`: on a fitted mesh the example's restart
differs by round-off also with one level. Without the checkpointed source the
velocities differ by :math:`10^{-4}` or more. Both run with ``erf.v = 1``, with
which ERF prints the force of the restricted source and of the model's source
(``structure coarse source:``); ``check_coarse_force.awk`` checks that they agree
to a relative :math:`10^{-12}`, and that the model's force is the 50 N it
applies, which also catches a Jacobian missing on both levels alike. The
unit tests ``StructureRegistry.*`` check the registry and the selection by
``erf.structures.model`` and the plot-name checks. The unit tests
``StructureCoarseSource.*`` check that the restriction keeps the force on every
component over a hill, one and two levels down, with refinement ratios
(2, 2, 1) and (2, 2, 2), boxes split and across a periodic boundary, and on a
uniform mesh, also on the faces inside the domain alone; that without the
Jacobian the hill changes the force (by :math:`3.8 \times 10^{-5}`,
:math:`1.6 \times 10^{-3}` and :math:`3.4 \times 10^{-3}` in x, y and z); that
the face Jacobian is that of the one cell inside the domain on the ground and
top faces and the mean of the two cells beside any other face; and
that the printed force matches one computed from the hill.

The fire view is tested in erf-model/FLARE, where the example fire engine
(``Tests/FireEngineExample``) is compiled in as well: four runs check it. Each runs the engine's deck with the example model
added: a disk of radius :math:`r(t) = 150 + 50\,t` m releasing
:math:`q_0 = 5 \times 10^4` W/m\ :sup:`2` inside it, a quarter of it as latent
heat. At each step the model prints the level and time :math:`t_f` of the
fire's fields, the heat released :math:`Q` (the sum of ``heat_release``) and
the burned area :math:`A` (``burned_area()``, which the example engine
computes from the disk's radius).
FLARE's ``Tests/StructureModelExample/check_fire_view.awk`` compares them with

.. math::

   A = \pi\, r(t_f)^2 ,
   \qquad
   Q = q_0\, \pi\, r\!\left(t_f - \tfrac12 \Delta t_f\right)^2 ,

where :math:`\Delta t_f` is the step of the fire's level (:math:`Q = 0` before
the fire's first step). The fire cells are 10 m (5 m on the refined level), and
they resolve the disk's heat to 1.4 % or better wherever it sits on the grid;
the tolerance is 3 %. Fields one fire step late or early would miss by 11 % or more.
``Structures_ExampleModel_FireView`` has both on level 0, so :math:`t_f` is
the model's step time. ``Structures_ExampleModel_FireViewFinerFire`` has the fire on
a level-1 patch and the model on level 0, so :math:`t_f` is again the model's
step time. ``Structures_ExampleModel_FireViewFinerModel`` has the fire on level
0 and the model on a level-1 patch, so :math:`t_f` is the end of the coarse
step. ``Structures_ExampleModel_FireViewExactTime`` has the fire on a level-1
patch with three substeps of 0.1/3 s and the model on level 0: three substeps
add up to 0.3 s where the coarse level is at 0.30000000000000004 s, and the
fire's :math:`t_f` must still equal the model's step time to all 17 printed
digits. Where :math:`t_f` is the model's step time, the check compares the two
exactly. ``Structures_ExampleModel_FireView`` then restarts from its mid-run
checkpoint: the first fields after the restart must hold the time of the restart
and the heat of the fire step before it, which the example engine restores.
``Structures_ExampleModel``, without a fire, checks that the model gets no fire
view.
