.. _sec:building:

Building ERF
============

ERF supports two build systems: GNUMake and CMake.

To use GNU Make, cd into ``ERF/Exec`` and type ``make``.
(Developers adding new features may wish to work in problem- or feature-specific directories in ``ERF/.Exec_dev`` instead;
for these cases one would build in that directory in ``ERF/.Exec_dev``.)

When on Cray HPC platforms like Perlmutter, CMake provides automated detection to simplify configuration.
To use cmake, one typically creates an ``ERF/build`` directory, copies an cmake file from ``ERF/Build``, and runs
cmake in ``ERF/build``.

This page provides quick-start commands; comprehensive documentation follows in the sections below.

**Where to start:**

* **Never built ERF before?** → Start with :ref:`sec:build:quickstart` for copyable clone-build-run snippets
* **Need general background?** → Read :ref:`sec:build:overview` to understand ERF's build system architecture
* **Need detailed control?** → Follow :ref:`sec:build:systems` for complete instructions from cloning to executable with custom configurations

Quick Start
-----------

**Most common builds:**

.. code-block:: bash

   # GNU Make with GNU compiler and MPI
   make COMP=gnu USE_MPI=TRUE

   # GNU Make with GPU support
   make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE

   # CMake configuration script
   cd Build && ./cmake.sh

**Cleanup:**

.. code-block:: bash

   # GNU Make
   make clean

   # CMake
   cd Build && make distclean

**Documentation structure:**

* :ref:`sec:build:quickstart` - **Quick start**: Copyable clone-build-run commands for common scenarios
* :ref:`sec:build:overview` - **Concepts**: Prerequisites, choosing build systems, and general workflow
* :ref:`sec:build:systems` - **Complete guide**: Detailed step-by-step instructions from cloning to executable with custom configurations
* :ref:`sec:build:hpc` - **HPC concepts**: Cray detection, machine profiles, and build script references
* :ref:`sec:hpc:guides` - **HPC system guides**: Perlmutter, Kestrel, and Aurora build/run walkthroughs
* :ref:`sec:build:library` - **Dependencies**: Configuring I/O libraries (NetCDF, HDF5) and physics packages (SHOC, P3, RRTMGP)
* :ref:`sec:build:troubleshooting` - **Help**: Common build errors, debugging tools, and getting assistance
