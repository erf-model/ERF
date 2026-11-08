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
       ``add_momentum_sources()`` on that level only. On a refined level the
       model's points and kernels must stay inside the level's boxes, also after
       a regrid; coarser levels see the force only through averaging down.
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
   * - ``add_momentum_sources()``
     - At every Runge-Kutta stage of the model's level: the model adds its force
       on the air [N/m\ :sup:`3`] to the face-centred momentum sources. The source
       is held fixed over the step. A share that falls inside an embedded
       boundary, immersed terrain or buildings is discarded there, as other
       momentum is.
   * - ``write_checkpoint()``
     - Each time ERF writes a checkpoint, under ``<checkpoint>/<model name>``.
       A model that restores its state in ``initialize()`` returns true from
       ``supports_restart()``; otherwise a restart stops the run.
   * - ``plot_var_names()``, ``fill_plot_var()``
     - The model's fields that ``erf.plot_vars_*`` may list, filled into the
       ERF plotfile on every level.

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
the mid-run checkpoint and compares the final plotfiles at zero tolerance. The
unit tests ``StructureRegistry.*`` check the registry and the selection by
``erf.structures.model``.
