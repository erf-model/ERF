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
       regrid. Coarser levels see the force only through averaging down. A finer
       level would never see it: averaging down would overwrite the forced cells
       with that level's unforced data.
   * - ``initialize()``
     - Once, after the initial data (or the checkpoint) exist and before the
       step-0 plotfile, and not again after a regrid (the model reads the mesh
       from ``StructureStep`` each step). ``StructureInit`` holds the level's
       geometry, the x-velocity, the mesh heights and Jacobians on a stretched or
       terrain-following mesh, whether every map factor is 1 (the actuator core
       assumes so), the turbulence closure, the checkpoint name on a restart and,
       on request, the immersed terrain surface.
   * - ``advance()``
     - Once per step of the model's level, with the state at the start of the
       step (``StructureStep``), before the step's Runge-Kutta stages. The model
       samples the wind, advances its structures and sets the momentum source.
       ``StructureStep`` also says whether every map factor of the level is 1
       now, so a level a regrid made is checked as well, and holds the fire view
       (``StructureStep::fire``, below) when the run has a fire.
   * - ``support()``
     - Right after each ``advance()``: boxes [m] that hold the model's points
       and their kernels' reach. When a level finer than the model's exists and
       its boxes overlap a support box, ERF stops the run and asks for the
       finer level (``structure.example.level`` for the example). On a
       terrain-following mesh the check uses the finer level's lowest and
       highest ground and spacing, so it errs on the side of stopping.
   * - ``add_momentum_sources()``
     - At every Runge-Kutta stage of the model's level: the model adds its force
       on the air [N/m\ :sup:`3`] to the face-centred momentum sources. The source
       is held fixed over the step. ERF applies no mask to it: the model keeps
       its kernels clear of embedded boundaries, immersed terrain and buildings.
   * - ``write_checkpoint()``
     - Each time ERF writes a checkpoint, under ``<checkpoint>/<model name>``.
       A model that restores its state in ``initialize()`` returns true from
       ``supports_restart()``; otherwise a restart stops the run.
   * - ``plot_var_names()``, ``fill_plot_var()``
     - The model's fields that ``erf.plot_vars_*`` may list, filled into the
       ERF plotfile on every level. The names carry the model's prefix; an
       empty name, a name listed twice or a name ERF or another model already
       uses stops the run.

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
lowest 50 m, below the kernel's reach, and must run. The
unit tests ``StructureRegistry.*`` check the registry and the selection by
``erf.structures.model`` and the plot-name checks.

When the example fire engine (``Tests/FireEngineExample``) is compiled in as
well, four runs check the fire view. Each runs the engine's deck with the example model
added: a disk of radius :math:`r(t) = 150 + 50\,t` m releasing
:math:`q_0 = 5 \times 10^4` W/m\ :sup:`2` inside it, a quarter of it as latent
heat. At each step the model prints the level and time :math:`t_f` of the
fire's fields, the heat released :math:`Q` (the sum of ``heat_release``) and
the burned area :math:`A` (``burned_area()``, which the example engine
computes from the disk's radius).
``check_fire_view.awk`` compares them with

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
