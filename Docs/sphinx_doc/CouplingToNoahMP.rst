   .. role:: cpp(code)
      :language: c++

.. _CouplingToNoahMP:

Coupling to Noah-MP
===================

Overview
--------

The Noah-MP Land Surface Model (LSM) is integrated with ERF to facilitate
interaction with the Noah-MP (Multi-Physics) land surface processes.

This documentation covers key components of this interface and its
associated data structures and routines, which are implemented in C++
and Fortran, and provides details on initialization and management of
data structures necessary for these processes.

Building and Running with Noah-MP
---------------------------------
To build ERF with Noah-MP support, the ``NetCDF``, ``NetCDF Fortran``, and ``HDF5`` libraries are
required. Furthermore, ``ERF_ENABLE_NOAHMP=ON`` must be specified with CMake builds or ``USE_NOAHMP=TRUE``
and ``USE_NETCDF=TRUE`` must be specified with GNU Make. Once an executable has been generated, the
inputs file for the simulation must specify Noah-MP as the land surface model type:

.. code-block:: bash

    erf.land_surface_model = "NOAHMP"

Currently, Noah-MP may only be utilized for simulations that are initialized from a WRF input file
(``erf.init_type = "WRFInput"``). Additionally, two files are required to be in the run directory
for Noah-MP initialization: ``namelist.erf`` and ``NoahmpTable.TBL``. Sample files are provided for
the :download:`namelist.erf <namelist.erf>` and :download:`NoahmpTable.TBL <NoahmpTable.TBL>`.

To improve computational efficiency, the Noah-MP timestep, specified via ``NOAH_TIMESTEP``
in the **namelist.erf** file, may be set larger than the ERF timestep to allow subcycling
in time. For example, if an 4s timestep is utilized for ERF and a 40s timestep is utilized for
Noah-MP, then Noah-MP will be updated every 10 steps.

Files Overview
--------------

-  **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.H**: Contains the declaration
   of the NOAH class, which extends the NullSurf.

-  **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.cpp**: Implements the
   initialization routine for the NOAH class and the per-step state exchange
   between ERF and Noah-MP.

The C++ ``↔`` Fortran coupling glue under **Submodules/Noah-MP/drivers/erf**
is no longer hand-written. Four files are **generated** at build time from a single
source of truth by ``tools/NoahmpIOCoupling.py`` (see
`Generating the C++–Fortran Coupling Glue`_). Only the tracked ``*-mc``
templates are edited by hand; the generated targets are git-ignored and recreated
on every build:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Template (tracked)
     - Generated target (git-ignored)
     - Role
   * - ``NoahmpIO.H-mc``
     - ``NoahmpIO.H``
     - C++ ``NoahmpIO_type`` declaration and the ``@NoahmpIOCoupling:Source``
       member list (the single source of truth)
   * - ``NoahmpIO.cpp-mc``
     - ``NoahmpIO.cpp``
     - C++ method implementations that forward into Fortran
   * - ``NoahmpIO_fi.F90-mc``
     - ``NoahmpIO_fi.F90``
     - Fortran ``bind(C)`` mirror struct and ``C_LOC``/``C_F_POINTER`` wiring
   * - ``NoahmpIOVarType.F90-mc``
     - ``NoahmpIOVarType.F90``
     - Fortran storage type holding the coupled members

-  **Submodules/Noah-MP/drivers/erf/NoahmpIO.H-mc**: Defines the C++
   ``NoahmpIO_type`` that is used to interface with Noah-MP implementations,
   following a similar structure as the underlying Fortran interface
   (https://dx.doi.org/10.5065/ew8g-yr95). The ``@NoahmpIOCoupling:Source { ... }``
   block inside this template is the canonical, ABI-ordered list of all coupled
   members.

-  **Submodules/Noah-MP/drivers/erf/tools/NoahmpIOCoupling.py**: The stdlib-only
   code generator (any Python 3, no third-party libraries) that expands the
   ``@NoahmpIOCoupling:Source`` block into the four generated targets above. The
   build runs it automatically; ``--check`` mode also drift-checks the generated
   files and audits the hand-written ``allocate()`` bounds.

-  **Submodules/Noah-MP/drivers/erf/NoahmpIOVarInitMod.F90**: Hand-written (not
   generated) Fortran module that allocates the coupled arrays. The allocate
   bounds must match each member's ``@NoahmpIOCoupling:bounds_Nd`` annotation;
   ``NoahmpIOCoupling.py --check`` audits this.

-  **Submodules/Noah-MP/drivers/erf/NoahmpReadNamelistMod.F90**: Hand-written
   Fortran module that reads ``namelist.erf``; it also holds the sentinel guards
   that let a C++-API value override a namelist default for ``scalar namelist``
   members.

-  **Submodules/Noah-MP/drivers/erf/NoahmpWriteLandMod.F90**: Hand-written
   Fortran module that controls which fields appear in the per-timestep NetCDF
   land output.

-  **Submodules/Noah-MP/drivers/erf/NoahmpWriteRestartMod.F90** and
   **NoahmpReadRestartMod.F90**: Fortran modules that serialize and restore the
   full Noah-MP prognostic state (soil, snow, canopy, aquifer, albedo history,
   ...) to/from a NetCDF restart file. They are exposed to C++ through the
   ``WriteRestart``/``ReadRestart`` methods of ``NoahmpIO_type`` and are used by
   ERF's checkpoint/restart capability (see :ref:`noahmp-checkpoint-restart`).

-  **Submodules/Noah-MP/drivers/erf/specs/**: TOML specs that double as
   AI-assistant task prompts and human references for common coupling chores
   (``add-coupled-variable.toml``, ``plot-land-output.toml``, and
   ``improve-coupling-codegen.toml``).

NOAHMP Class
------------

The NOAH class serves as the handler for initializing and managing the
data structures required for NOAH-MP operations. It inherits from the
`NullSurf` class. This class declares private variable `NoahmpIO_type
noahmpio`, that is passed to NoahMP routines similar to the Fortran
interface in the Noah-MP documentation
(https://dx.doi.org/10.5065/ew8g-yr95)

NoahmpIO_type Structure
-----------------------

This structure is key for handling the input and output operations for
NOAH-MP through C++ and Fortran interoperation. Contains various
variables for domain, memory, and tile configuration. Also, contains
arrays for geographic variables. At present this type exposes only a
select set of variables. More variables should be exposed as needed by
applications in ERF.

All coupled members are declared **once**, in ABI order, inside the
``@NoahmpIOCoupling:Source { ... }`` block of
**Submodules/Noah-MP/drivers/erf/NoahmpIO.H-mc**. The member's kind is inferred
from its C++ type and refined by a trailing ``@NoahmpIOCoupling:`` annotation:

.. code-block:: c++

    int ids, ide, jds, jde;                   // dimensions/handles -- no annotation
    noahmp_real DTBL;                          // @NoahmpIOCoupling:scalar
    noahmp_real ZLVL = -9999.0;                // @NoahmpIOCoupling:scalar
    NoahArray2D<noahmp_real> XLAT;             // @NoahmpIOCoupling:bounds_2d(xstart:xend, ystart:yend)
    NoahArray3D<noahmp_real> T_PHY;           // @NoahmpIOCoupling:bounds_3d(xstart:xend, kms:kme, ystart:yend)

From this one ordered list the generator re-derives the mirrored ``fi`` struct,
both constructors, the move constructor, ``NOAHMP_IO_FI_NUM_MEMBERS``, the Fortran
``bind(C)`` type, the coupled members of the storage type, the
``C_LOC``/``C_F_POINTER`` wiring, and the C++ ``NoahArray`` extents. Because the
C++ struct and the Fortran ``bind(C)`` type are emitted from the same list in one
run, their member count, order, and types are identical *by construction*.

To expose a new variable you therefore add a single annotated line here and
rebuild; the only remaining manual steps are the runtime wiring the generator does
not own (array allocation, namelist guard, optional NetCDF output). See
`Generating the C++–Fortran Coupling Glue`_ for the full procedure.

.. note::

   Always use ``noahmp_real`` (C++) / ``c_kind_noahmp`` (Fortran) for coupled
   floating-point members, never a hardcoded ``double``/``C_DOUBLE``. Both are
   selected by the ``DOUBLE_PREC`` build macro so the boundary precision tracks
   ``amrex::Real``. ERF builds Noah-MP in double precision; a run-time guard
   (``NoahmpIO_AssertAbi``) aborts if the two compilers disagree.

Fortran Interoperability
------------------------

The connection between C++ and Fortran is managed through `NoahmpIO_fi`.
This module contains a mirroring of the C++ structure for NOAH-MP
input-output operations.

Operations are invoked as **methods on** ``NoahmpIO_type``; each forwards to its
``*_fi`` Fortran implementation through the generated ``bind(C)`` glue. The main
ones are:

-  ``void ScalarInitDefault()`` / ``void VarInitDefault()``: Initialize default
   scalar and array members of ``NoahmpIO_type`` and create the C pointers into
   the Fortran-owned storage.

-  ``void InitMain()``: Main initialization routine for the NOAH-MP operations.

-  ``void ReadNamelist()`` / ``void ReadTable()``: Read ``namelist.erf`` and
   ``NoahmpTable.TBL``.

-  ``void ReadLandHeader()`` / ``void ReadLandMain()``: Read the land input data.

-  ``void DriverMain()``: Advance the Noah-MP physics one land-surface step.

-  ``void WriteLand(int filenum)``: Write the per-timestep NetCDF land output.

-  ``void WriteRestart(const std::string& dir)`` /
   ``void ReadRestart(const std::string& dir)``: Serialize and restore the full
   prognostic state for checkpoint/restart (see :ref:`noahmp-checkpoint-restart`).

For applications that manage one ``NoahmpIO_type`` per domain tile, the
``NoahmpIO_vector`` container mirrors a Fortran allocatable array and is sized
exactly once via ``resize(size, level)``.

Usage
-----

To use the NOAH class and associated functions, ensure the correct
initialization sequence is followed within the simulation setup. The
interplay between C++ and Fortran necessitates careful memory and data
handling, which is crucial for ensuring performance and correctness in
simulations. The interface is designed to mimic the Fortran interface
described in the Noah-MP documentation, therefore similar practices should
be followed.

.. _noahmp-checkpoint-restart:

Checkpoint and Restart
----------------------

Noah-MP participates in ERF's standard checkpoint/restart capability (see
:ref:`sec:Checkpoint` for the general ``erf.check_file``, ``erf.check_int``,
``erf.check_per``, and ``amr.restart`` controls). No additional input options
are required to checkpoint or restart a Noah-MP simulation; whenever a
checkpoint is written, the land-surface state is written alongside the
atmospheric state, and it is read back automatically on restart.

A Noah-MP restart is **bitwise reproducible**: restarting from a checkpoint
reproduces the trajectory of an equivalent cold-start run exactly. Achieving
this requires persisting two pieces of state in addition to the regular ERF
MultiFab data:

#. **The land-surface substep counter** (Noah-MP's ``itimestep``). When the
   Noah-MP timestep is larger than the ERF timestep (subcycling via
   ``NOAH_TIMESTEP``; see `Building and Running with Noah-MP`_), the counter
   determines when Noah-MP fires relative to the ERF steps. It is written to
   ``<chkfile>/lsm_step`` (one value per AMR level) and restored on restart so
   that the firing schedule, and hence the LSM-to-atmosphere flux timing, is
   preserved.

#. **The full Noah-MP prognostic state.** The complete land state — soil
   temperature and moisture, the snowpack (including the active snow layers and
   the layer count), canopy and vegetation variables, the aquifer, albedo
   history, phenology, and accumulators — is serialized at the model's working
   precision to ``<chkfile>/noahmp_restart/Level_<lev>.nc`` using
   ``NoahmpWriteRestart``/``NoahmpReadRestart``. Each local block writes its
   tile into the global-domain NetCDF file collectively.

On restart, ERF first cold-initializes the land state from the WRF input file
and tables (as it would for a fresh run) and then overwrites it with the
checkpointed state; Noah-MP's per-step input transfer pulls the restored state
into the physics on the first ``Advance``.

.. note::

   These per-level files (``lsm_step`` and ``noahmp_restart/``) are written into
   the checkpoint directory and are not AMReX MultiFabs. If you copy or archive a
   checkpoint manually, make sure these are included.

.. note::

   **Legacy checkpoints.** Checkpoints written before this capability was added
   do not contain the ``lsm_step`` file or the ``noahmp_restart`` directory. ERF
   detects their absence and falls back to the previous behavior — the substep
   counter resets to zero and the Noah-MP state is cold-initialized from the WRF
   input — printing a warning that the restarted land trajectory will differ from
   a cold start. Such restarts are therefore not bitwise reproducible.

Generating the C++–Fortran Coupling Glue
========================================

The C++ ``↔`` Fortran ABI under **Submodules/Noah-MP/drivers/erf** is produced by a
deterministic, stdlib-only code generator,
**tools/NoahmpIOCoupling.py** — *not* by an LLM. The generator reads the single
``@NoahmpIOCoupling:Source { ... }`` block in ``NoahmpIO.H-mc`` and expands every
``@NoahmpIOCoupling:<region>`` marker in the four ``*-mc`` templates into the
corresponding generated target (see the table under `Files Overview`_).

The generated targets (``NoahmpIO.H``, ``NoahmpIO.cpp``, ``NoahmpIO_fi.F90``,
``NoahmpIOVarType.F90``) are git-ignored and recreated by the build:

-  **GNU Make** runs the generator at parse time (before source discovery).
-  **CMake** runs it at configure time.

A fresh checkout therefore builds with no manual step. Two convenience targets are
available for editing or CI:

.. code-block:: bash

   # GNU Make (from Submodules/Noah-MP/drivers/erf)
   make codegen          # regenerate the four targets in place
   make codegen-check    # exit 1 if regeneration would change anything (drift check)

   # CMake
   cmake --build <dir> --target noahmp_codegen
   cmake --build <dir> --target noahmp_codegen_check

``codegen-check`` does double duty: it fails if the generated targets are out of
sync with the ``@NoahmpIOCoupling:Source`` block, **and** it audits the
hand-written ``allocate()`` bounds in ``NoahmpIOVarInitMod.F90`` against each
array's ``@NoahmpIOCoupling:bounds_Nd`` annotation. Run it in CI.

Adding a coupled variable
-------------------------

The full, authoritative procedure (and an AI-assistant prompt) lives in
**Submodules/Noah-MP/drivers/erf/specs/add-coupled-variable.toml**. In short:

#. **Declare it (boundary glue).** Add ONE annotated line to the
   ``@NoahmpIOCoupling:Source`` block in ``NoahmpIO.H-mc``, in the ABI position you
   want (order in the block = ABI order):

   .. code-block:: c++

      noahmp_real RAINBL;                  // @NoahmpIOCoupling:scalar namelist "rainfall [mm]"
      NoahArray2D<noahmp_real> RAINLSM;    // @NoahmpIOCoupling:bounds_2d(xstart:xend, ystart:yend) "lsm rain"
      NoahArray3D<noahmp_real> SMOIS;      // @NoahmpIOCoupling:bounds_3d(xstart:xend, nsoil:nsoil, ystart:yend)

   The generator infers the kind from the C++ type (``int`` /
   ``noahmp_real`` / ``NoahArray{2,3}D``). Array bounds must be listed in Fortran
   order and every bound token must be a literal or another coupled member name
   (``xstart``, ``xend``, ``kms``, ``kme``, ``nsoil``, ``numrad``, ...) so C++ and
   Fortran resolve it to the same value. Scalars may add ``namelist`` and an
   optional ``"doc"`` string; ints (dimensions/handles) need no annotation and may
   carry a C++ default (e.g. ``int numrad = 2;``).

#. **Regenerate.** ``make codegen`` (or ``cmake --build ... --target
   noahmp_codegen``); a plain build does this automatically. Use ``codegen-check``
   to verify.

#. **Runtime wiring (manual, not generated).** For a new array, add an
   ``allocate()`` in ``NoahmpIOVarInitMod.F90`` with bounds matching the
   ``@NoahmpIOCoupling:bounds_Nd`` annotation (scalars are not allocated). For a
   ``scalar namelist`` variable, add a sentinel guard in
   ``NoahmpReadNamelistMod.F90`` so a C++-API value wins over the namelist default.

#. **(Optional) NetCDF land output.** To make the variable appear in the
   per-timestep land output, add it to ``NoahmpWriteLandMod.F90`` following the
   existing ``TSK`` / ``SMOIS`` pattern (see
   ``specs/plot-land-output.toml``).

To extend the generator itself — a new annotation kind, a new generated region, or
better diagnostics — see ``specs/improve-coupling-codegen.toml``.

Generating the ERF driver with CodeScribe
=========================================

The ERF-side C++ driver, **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.cpp** and
**ERF_NOAHMP.H**, is *not* covered by the macroprocessor above (it contains the
GPU-aware, component-indexed state exchange rather than flat ABI plumbing). It can
optionally be updated with **CodeScribe**
(https://github.com/akashdhruv/CodeScribe), an LLM-based code-update tool.

Follow the instructions in the `CodeScribe repository <https://github.com/akashdhruv/CodeScribe>`_
to configure your LLM environment. You will need API access for your preferred
model (e.g., OpenAI, Argo, etc.). Tutorials are available at
`https://github.com/akashdhruv/codescribe-tutorial <https://github.com/akashdhruv/codescribe-tutorial>`_.
The driver's design specifications — the component-indexed field enums, the
GPU-aware state exchange, the lifecycle, and the contract a code-update must
respect — live as developer specs under
**Source/LandSurfaceModel/Noah-MP/dev/** (start with ``dev/README.md``):

.. code-block:: bash

   code-scribe update ERF_NOAHMP.cpp ERF_NOAHMP.H \
       -q "Write a natural language prompt with variable names, dimensions, etc." \
       -m <openai|argo-gpt4o|...>

.. note::

   CodeScribe is no longer used for the ``Submodules/Noah-MP/drivers/erf`` coupling
   glue (``NoahmpIO.H``, ``NoahmpIO.cpp``, ``NoahmpIO_fi.F90``,
   ``NoahmpIOVarType.F90``). That glue is now generated deterministically by
   ``tools/NoahmpIOCoupling.py`` from the ``@NoahmpIOCoupling:Source`` block; do not
   hand-edit the generated targets. The previous manual
   ``real(kind=kind_noahmp)`` ``→`` ``real(kind=C_DOUBLE)`` substitution is also
   obsolete — coupled reals use ``noahmp_real`` / ``c_kind_noahmp``, whose precision
   is selected by the ``DOUBLE_PREC`` build flag.
