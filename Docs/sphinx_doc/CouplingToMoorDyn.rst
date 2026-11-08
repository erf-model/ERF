
.. role:: cpp(code)
   :language: c++

.. _sec:CouplingToMoorDyn:

Coupling To MoorDyn
===================

ERF couples to `MoorDyn-C <https://moordyn.readthedocs.io>`_ (the C++ MoorDyn,
version 2) to simulate flexible lines in the wind, such as overhead conductors
and shield wires dead-ended at fixed attachment points, as single spans or as
sections hanging from insulator strings at suspension towers. MoorDyn
integrates the lumped-mass line dynamics; ERF supplies the fluid velocity at
MoorDyn's line nodes through its external wave-kinematics interface, so the
"water" MoorDyn's lines hang in is ERF's air. The actuator core of
:ref:`sec:ActuatorCore` provides the sampling of the wind at the nodes, and a
body model (one that writes a MoorDyn input for its lines and steps them with
ERF) builds on the wrapper described here; the wrapper itself has no caller in
ERF, so nothing in a run reads MoorDyn until such a model is in the tree. This
page is the how-to for building ERF with MoorDyn and for what the wrapper in
``Source/MovingBodies/MoorDyn`` checks.

Building
--------

Build and install MoorDyn-C 2.7.1 (any 2.3 or newer release works; the script
refuses an older version) with the script in the ERF tree, which clones the release tag, configures the C/C++
library only (no Python, MATLAB, Fortran or Rust wrappers, no tests, no docs,
the bundled Eigen) and installs it under ``$HOME/opt/moordyn-<version>``:

.. code-block:: bash

   Build/setup_moordyn.sh --version 2.7.1

Then build ERF with the coupling. ``ERF_ENABLE_MOORDYN`` switches on the
actuator core as well (``ERF_ENABLE_MOVING_BODIES``); ``MOORDYN_DIR`` is the
install prefix, where CMake finds MoorDyn's package configuration
(``lib/cmake/moordyn``) and the headers (``include/moordyn``):

.. code-block:: bash

   cmake -DERF_ENABLE_MPI=ON -DERF_ENABLE_FFT=ON -DERF_ENABLE_MOORDYN=ON \
         -DMOORDYN_DIR=$HOME/opt/moordyn-2.7.1 ..

The configure step records the MoorDyn version it found and refuses a MoorDyn
whose ``MoorDyn2.h`` does not declare ``MoorDyn_ExternalWaveKinSet``, one of
the external wave-kinematics calls (``MoorDyn_ExternalWaveKin*``) that
replaced the version 1 names in 2.3.
The version is compiled into ERF and reported by
:cpp:`erf_moordyn::library_version()`, for a body model to print at start-up.

MoorDyn is a shared library. The CMake build carries the link path of
``libmoordyn`` in the executables' run path, so no loader variable is needed
when the install stays where it was built; if the library is moved, and for
a GNU make build (which only adds ``-L$(MOORDYN_HOME)/lib``), set
``DYLD_LIBRARY_PATH`` (macOS) or ``LD_LIBRARY_PATH`` (Linux) to its ``lib``
directory.

Without a MoorDyn installation, ``-DERF_MOORDYN_USE_STUB=ON`` builds the
coupling against a bundled stub library with the same C API
(``Source/MovingBodies/MoorDyn/Stub``). The stub reads the same input file
and has the geometry, ordering and data flow of the real interface but no
line dynamics: each line hangs as an elastic parabola, as long as its
unstretched length stretched by its tension (slack or taut), swings to the
quasi-static blowout angle of the wind it is given with a one-second lag,
and carries that tension along its tangent. A coupled point moves as MoorDyn
moves it. A free point hangs from the shortest line that ties it to a fixed or
coupled point (an insulator string from its tower), plumb below it and swung
across the line by the wind across the other line attached to it (on that
line's weight, with the same lag), and that line runs straight down to it. Its saved state holds
each line's swing angles, the last wind it was given, which sets the direction
of the swing, and where the free and coupled points are, so a restored line is where
the saved one was. The ``Linux GCC MoorDyn`` CI workflow builds the wrapper
and runs the unit tests against the stub in one job; a second job installs
MoorDyn-C 2.7.1 with ``Build/setup_moordyn.sh`` and runs them against the real
library.

The GNU make build takes ``USE_MOVING_BODIES = TRUE`` and
``USE_MOORDYN = TRUE`` (the latter is ignored without the former), with
``MOORDYN_HOME`` set to the install prefix (and ``MOORDYN_VERSION`` to the version string it should
report).

What the coupling layer does
----------------------------

:cpp:`erf_moordyn::MoorDynSystem` wraps one MoorDyn system (one input file)
behind the version 2 C API with ERF's error handling. Creating the system
and initialising it report their failure in a message that names the input
file, so that a body can abort with its own name; every later call aborts
with the MoorDyn function and the error code's name
(``MOORDYN_INVALID_INPUT``, ``MOORDYN_NAN_ERROR``, ...), since nothing
sensible can follow a failed line step. The wrapper exposes the coupled
degrees of freedom (none when every attachment is fixed), the external
kinematics (the number of points MoorDyn wants the fluid velocity at, their
current coordinates, and the velocity and acceleration to set there before
each step), the step itself (MoorDyn sub-steps internally with its own
``dtM`` within the step it is given), the line node positions, velocities,
tensions and net forces, the end and maximum tensions, the attachment points' positions and
forces, and MoorDyn's save and load of the whole system state for restarts.

MoorDyn-C is not clean under floating-point traps: in 2.7.1 its stationary
initial-condition solver bounds its pseudo-step by ``cfl * length / v`` for
every point, whose length is the largest double, and the quotient overflows to
infinity whenever the residual speed ``v`` is below ``cfl`` (it divides by
zero when ``v = 0``); the ``std::min`` it goes into discards the infinity, so
the result is sound, but with ``amrex.fpe_trap_overflow = 1`` the run stopped
inside ``MoorDyn_Init`` (SIGFPE on Linux, SIGILL on Apple arm64). The wrapper
therefore masks the traps around MoorDyn's computing calls:
:cpp:`erf_moordyn::FPTrapMask` switches ``amrex.fpe_trap_invalid``,
``amrex.fpe_trap_zero`` and ``amrex.fpe_trap_overflow`` off with
``amrex::disableFPExcept``, and on leaving discards the exception flags
MoorDyn raised (or a re-enabled trap could fire later in unrelated code), puts
back the flags raised before and restores the traps with
``amrex::setFPExcept``. ``MoorDynSystem::init`` (``MoorDyn_Init`` and
``MoorDyn_Init_NoIC``) and ``MoorDynSystem::step`` (``MoorDyn_Step``) run under
it; the wrapper's other calls (creation, the kinematics, the getters, save,
load and the in-memory state) ran under all three traps in ERF's unit tests
with MoorDyn-C 2.7.1 without raising one and are not masked. A run with lines may therefore keep the traps on for ERF's own
code, in the stub build and with MoorDyn-C alike; masking changes no answer
(the masked and untrapped runs agree bit for bit). MoorDyn's own failures are
still caught: a line integration that diverges (MoorDyn's internal step too
long for its segments) overflows and makes NaNs inside ``MoorDyn_Step``, which
reports ``MOORDYN_NAN_ERROR``, and the wrapper aborts naming MoorDyn's
internal step and the settings that bound it (``CFL`` and ``dtM`` in the input
file, and a body model's own keys for them; with the traps on and no mask the
run stopped there with a SIGFPE inside MoorDyn); ``step`` also aborts
on a non-finite coupled force. The traps belong to the calling thread, the
thread MoorDyn runs on (MoorDyn-C built without OpenMP, as
``Build/setup_moordyn.sh`` builds it, starts none); ``amrex::disableFPExcept``
sets them on Linux (glibc) and macOS only, as AMReX does.

MoorDyn reports the net force of the points it integrates (free and coupled
points) only; a fixed attachment reports zero, and its pull is the net force
MoorDyn finds on the end node of each line attached to it
(``MoorDyn_GetLineNodeForce``): the end segment's tension with the node's
share of the weight and the drag, which the fixed point holds still. A
coupled point is one the caller moves: each ``MoorDyn_Step`` takes its
position at the start of the step and a velocity, moves it linearly over the
step, and returns the net force of the attached lines' end nodes on it, the
same sum (without the end nodes' inertia). A tower model that moves its
cross-arm uses this: the cross-arm is a coupled point (:ref:`sec:TowerFrame`,
"Coupling to the conductor lines"). A system's whole state can also be kept in
memory and restored (``MoorDyn_Serialize`` and ``MoorDyn_Deserialize``, the
stub too), so that a coupled step can be redone by a coupling that iterates
each step until the towers and the lines agree.

Two properties of MoorDyn's input matter for lines in air:

- the ``OPTIONS`` must set ``WaveKin = 1`` (the fluid kinematics come through
  the API). Creating the system refuses an input file whose ``OPTIONS`` do not
  set it: MoorDyn-C 2.7.1 would otherwise run without the wind ERF hands it;
- MoorDyn applies its fluid loads to nodes below ``z = 0`` only, its free
  surface, and the flat bottom lies at ``-WtrDpth``. Lines in air therefore
  live below ``z = 0`` in MoorDyn's frame, with ``WtrDnsty`` the air density
  and ``WtrDpth`` deep enough that the bottom is below the ground; a body
  model that writes the input shifts ERF's heights accordingly.

The unit test
-------------

The unit test ``MoorDynSystem`` (``Tests/Unit/MovingBodies``) runs a 300 m
fixed-fixed span of 795 kcmil ACSR with 1.5 m of slack in air against
whichever library the build links: the span initialises with no coupled
degree of freedom, its end nodes on the attachment points and the catenary
sag at mid-span; the external kinematics points start with the line's nodes;
a 20 m/s crosswind blows the span out towards the quasi-static angle
``atan(q / w)`` of the drag per unit length ``q`` over the weight per unit
length ``w`` and raises the end tension; and a saved state restored into a
fresh system created from the same input continues identically.
