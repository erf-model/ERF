.. _sec:DustCoupling:

Dust coupling
=============

The dust coupling (``Source/DustCoupling``) connects ERF to a model of dust
emitted from the ground: a model that lifts dust into the air, settles and
deposits it, and writes its own diagnostics. ERF contains no dust model. A
model lives in its own directory outside ERF, registers itself by name, and
plugs into a small interface, ``DustModel`` (``ERF_DustModel.H``), in the way
fire engines plug into the fire coupling (:ref:`sec:FireCoupling`) and
structure models into the structure coupling (:ref:`sec:StructureCoupling`).
FLARE's dust layer is such a model.

Building
--------

The coupling is built with ``-DERF_ENABLE_DUST=ON`` (CMake) or
``USE_DUST = TRUE`` (GNU make). It adds dust scalars to the conserved state,
``ERF_NUM_DUST_SCALARS`` of them (CMake; ``NUM_DUST_SCALARS`` in GNU make; 1 by
default). Their number is fixed when ERF is compiled, as the smoke scalar of
the fire coupling is. The dust scalars sit right after the base scalar
(``rhoadv_0``) and before smoke: ``RhoDust_comp = RhoScalar_comp + 1``. They
are passive scalars, advected and diffused with the others (with
``erf.transport_scalar = true``), and are plotted as ``rhoadv_dust`` (or
``rhoadv_dust_0``, ``rhoadv_dust_1``, ... when there are several). A checkpoint
records the number of scalars, and a build with a different number refuses to
restart from it.

A model's directory is compiled in with ``ERF_DUST_MODEL_DIRS`` (CMake, a list
of directories, each with an ``ERFDustModel.cmake`` that defines
``erf_dust_model_add(target)``) or ``DUST_HOME`` (GNU make, each directory with
a ``Make.package``). A minimal model that shows both files is in
erf-model/FLARE (``Tests/DustModelExample``).

Inputs
------

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Parameter
     - Definition
     - Default
   * - ``erf.dust.enable``
     - Turn dust on. A build without the dust coupling stops the run when it
       is true.
     - false
   * - ``erf.dust.model``
     - The dust model, by the name it registers. It may be left out when one
       model is compiled in. A name no compiled-in model has stops the run,
       listing the names that are.
     - the one model compiled in

Each model reads its own inputs.

What ERF does and what the model does
-------------------------------------

ERF owns the dust scalars: it transports them, sets them to zero at the start
of a run, checkpoints them and sums them (``RHO DUST`` in the integrated
quantities, with ``erf.sum_interval``). The model owns everything else: its
surface grid and fields, its particles, its output and its own checkpoint
data. ERF calls it at fixed points of a run:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ERF call
     - What happens
   * - ``initialize()``
     - Once, at the end of start-up, after the fire coupling has started (a
       model may read the fire engine's fields). The model reads its inputs
       and sets itself up; on a restart it reads back what its
       ``write_checkpoint()`` wrote. On a fresh start ERF then zeroes the dust
       scalars on every level.
   * - ``add_sources()``
     - At every Runge-Kutta stage of every level, after ERF has built its own
       sources and before the fire coupling adds its fluxes. The model adds
       its emission, settling and deposition to the dust components of the
       source. A model that runs on level 0 only returns on the other levels.
   * - ``advance()``
     - Once per step of level 0, after the step, with the state at its end:
       the model advances its own fields over the step.
   * - ``write_output()``
     - At start-up (step 0, before the first plotfile), after each step and
       at the end of the run, each time before ERF writes a checkpoint, so that
       a checkpoint records the model's output of that step as written.
   * - ``write_checkpoint()``
     - Each time ERF writes a checkpoint, into the checkpoint directory.

With the fire coupling, the calls also carry a read-only view of the fire: the
engine (its level set, reference wind and heat release, where the engine has
them) and the fire grid.

Tests
-----

ERF builds no dust model and runs no dust tests; they live in erf-model/FLARE
with its dust layer. There, ``Tests/DustModelExample`` emits a constant flux
from the ground into a closed box: ``Dust_ExampleModel_Mass`` checks that the
dust ERF sums (``RHO DUST``) equals the mass the model emitted at every step,
and ``Dust_ExampleModel_Restart`` that a restart from the mid-run checkpoint
reproduces the final plotfile and the model's own count exactly. In ERF,
``Dust_Abort_BuiltWithoutDust`` checks that a build without the coupling
refuses a dust deck.
