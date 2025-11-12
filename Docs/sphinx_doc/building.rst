.. _Building:

Building
--------

The ERF code is dependent on `AMReX <https://github.com/AMReX-Codes/amrex>`_.

If radiation is used, ERF includes the radiation model `RTE-RRTMGP <https://github.com/E3SM-Project/E3SM/tree/master/components/eamxx/src/physics/rrtmgp>`_ also used by E3SM.

ERF can also use the `simplified-higher-order-closure (SHOC) turbulence and cloud macrophysics scheme from E3SM <https://github.com/E3SM-Project/E3SM/tree/master/components/eamxx/src/physics/shoc>`_ , as well as the
`P3 microphysics scheme scheme from E3SM <https://github.com/E3SM-Project/E3SM/tree/master/components/eamxx/src/physics/p3>`_ .

RRTMGP, SHOC and P3 use Kokkos for heterogeneous computing infrastructures.

AMReX, EKAT, NOAH-MP and RTE-RRTMGP are all available as submodules in the ERF repo; using SHOC and P3 requires extra steps described below.
Kokkos is accessed as a submodule of the EKAT submodule.

ERF can be built using either GNU Make or CMake.

Minimum Requirements
~~~~~~~~~~~~~~~~~~~~

ERF requires a C++ compiler that supports the C++17 standard and a C compiler that supports the C99 standard.
Building with GPU support may be done with CUDA, HIP, or SYCL.
For CUDA, ERF requires versions >= 11.0. For HIP and SYCL, only the latest compilers are supported.
Prerequisites for building with GNU Make include Python (>= 2.7, including 3) and standard tools available
in any Unix-like environments (e.g., Perl and sed). For building with CMake, the minimal requirement is version 3.18.

   .. note::
      **While ERF is designed to work with SYCL, we do not make any guarantees that it will build and run on your Intel platform.**

   .. note::
      **ERF was successfully compiled with the Intel compiler suite (e.g., icx
      version 2024.1.0). However, for older versions, it may be necessary to
      use reduced compiler optimization** (``-O1``) **to avoid an internal compiler
      error.** For example, ERF was successfully compiled with icpc version
      19.1.2.254, with ``-O2`` (``CMAKE_BUILD_TYPE = RelWithDebInfo``) but
      TimeIntegration/ERF_advance_dycore.cpp had to be manually compiled with
      ``-O1``. Your mileage may vary.

Paradigm
~~~~~~~~~~

ERF uses the paradigm that different executables are built in different subdirectories within the ``Exec`` directory.  When
using gmake (see below), the user/developer should build in the directory of the selected problem.  When using
cmake (see below), separate executables are built for all of the problem directories listed in ``Exec/CMakeLists.txt``.
The problem directories within ``Exec`` are sorted into 1) science-relevant setups, such as ``ABL`` for modeling the atmospheric
boundary layer, 2) dry and moist regression tests in ``Exec/DryRegTests`` and ``Exec/MoistRegTests`` respectively,
that are used for testing specific known aspects of the code functionality,
such as boundary conditions or Rayleigh damping, and 3) tests for features under development in ``Exec/DevTests``,such as moving terrain.
There is a README in each problem directory that describes the purpose/role of that problem.

GNU Make
~~~~~~~~

The GNU Make system is best for use on large computing facility machines and production runs.
With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations
particular to that machine if possible. These settings are kept up-to-date by the AMReX project.

Using the GNU Make build system involves first setting environment variables for the directories of the dependencies of ERF.

All dependencies except for SHOC and P3 are provided as git submodules in ERF and can be populated by using
``git submodule update --init --recursive`` in the ERF repo, or before cloning by using ``git clone --recursive <erf_repo>``.
Although submodules of these projects are provided, they can be placed externally as long as the ``<REPO_HOME>``
environment variables for each dependency is set correctly.
An example of setting the ``<REPO_HOME>`` environment variables in the user's ``.bashrc`` is shown below:

::

   export ERF_HOME=${HOME}/ERF
   export AMREX_HOME=${ERF_HOME}/Submodules/AMReX

The GNU Make system is set up to use the path to AMReX submodule by default, so it is not necessary to set
the AMReX path explicitly. It is also possible to use an external version of AMReX, downloaded by running

   .. code:: shell

             git clone https://github.com/amrex-codes/amrex.git

in which case the ``AMREX_HOME`` environment variable must point to the location where AMReX has been downloaded, which will take precedence over the default path to the submodule. If using bash shell,

::

   export AMREX_HOME=/path/to/external/amrex

or if using tcsh,

::

   setenv AMREX_HOME /path/to/external/amrex

To build with SHOC or P3 using gmake

::

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh
   source /path/to/ERF/Build/GNU_Ekat/ekat_build_commands.sh

Then follow the instructions below, ensuring that you have ``USE_SHOC=TRUE`` (when running with shoc)
and ``USE_P3=TRUE`` (when running with p3) in your GNUmakefile.

#. ``cd`` to the desired build directory, e.g.  ``ERF/Exec/DryRegTests/IsentropicVortex/``

#. Edit the ``GNUmakefile``; options include

   +--------------------+------------------------------+------------------+-------------+
   | Option name        | Description                  | Possible values  | Default     |
   |                    |                              |                  | value       |
   +====================+==============================+==================+=============+
   | COMP               | Compiler (gnu or intel)      | gnu / intel      | None        |
   +--------------------+------------------------------+------------------+-------------+
   | USE_MPI            | Whether to enable MPI        | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_OMP            | Whether to enable OpenMP     | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_CUDA           | Whether to enable CUDA       | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_HIP            | Whether to enable HIP        | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_SYCL           | Whether to enable SYCL       | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_NETCDF         | Whether to enable NETCDF     | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_PARTICLES      | Whether to enable particles  | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_RRTMGP         | Whether to enable radiation  | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_SHOC           | Whether to enable SHOC       | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_P3             | Whether to enable P3         | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_MULTIBLOCK     | Whether to enable multiblock | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | DEBUG              | Whether to use DEBUG mode    | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | PROFILE            | Include profiling info       | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | TINY_PROFILE       | Include tiny profiling info  | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | COMM_PROFILE       | Include comm profiling info  | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | TRACE_PROFILE      | Include trace profiling info | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+

   .. note::
      **At most one of USE_OMP, USE_CUDA, USE_HIP, USE_SYCL should be set to true.**

   Information on using other compilers can be found in the AMReX documentation at
   https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html .

#. Make the executable by typing

   .. code:: shell

      make

   The name of the resulting executable (generated by the GNUmake system) encodes several of the build characteristics, including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable.
   Thus, several different build configurations may coexist simultaneously in a problem folder.
   For example, the default build in ``ERF/Exec/DryRegTests/IsentropicVortex`` will look
   like ``ERF3d.gnu.MPI.ex``, indicating that this is a 3-d version of the code, made with
   ``COMP=gnu``, and ``USE_MPI=TRUE``.

GNU Make Docs
~~~~~~~~~~~~~
Building the ERF documentation with GNU make can be completed by navigating to the ``/ERF/Docs/`` directory and executing the following command:

.. code:: shell

   source BuildDocs.sh

Note that the sphinx as well as the doxygen documentation will be built.

Job info
~~~~~~~~

The build information can be accessed by typing

   .. code:: shell

      ./ERF*ex --describe

in the directory where the executable has been built.


CMake
~~~~~

CMake is often preferred by developers of ERF; CMake allows for building as well as easy testing and verification of ERF through the use of CTest which is included in CMake.

Compiling with CMake involves an additional configure step before using the ``make`` command and it is expected that the user has cloned the ERF repo with the ``--recursive`` option or performed ``git submodule init; git submodule update`` in the ERF repo to populate its submodules.

ERF provides example scripts for CMake configuration in the ``/path/to/ERF/Build`` directory.  Once the CMake configure step is done, the ``make`` command will build the executable.

An example CMake configure command to build ERF with MPI is listed below:

::

    cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DERF_ENABLE_MPI:BOOL=ON \
          -DCMAKE_CXX_COMPILER:STRING=mpicxx \
          -DCMAKE_C_COMPILER:STRING=mpicc \
          -DCMAKE_Fortran_COMPILER:STRING=mpifort \
          .. && make

Typically, a user will create a ``build`` directory in the project directory and execute the configuration from said directory (``cmake <options> ..``) before building.  Note that CMake is able to generate makefiles for the Ninja build system as well which will allow for faster building of the executable(s).

To build with SHOC or P3 using cmake, you will need to make sure you have run ``git submodule update --init --recursive``, then

::

   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/eamxx_clone.sh

Then follow the guidance below, making sure to set ``ERF_ENABLE_SHOC`` and/or ``ERF_ENABLE_P3`` to TRUE.

Analogous to GNU Make, the list of cmake directives is as follows:

   +---------------------------+------------------------------+------------------+-------------+
   | Option name               | Description                  | Possible values  | Default     |
   |                           |                              |                  | value       |
   +===========================+==============================+==================+=============+
   | CMAKE_BUILD_TYPE          | Whether to use DEBUG         | Release / Debug  | Release     |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_MPI            | Whether to enable MPI        | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_OPENMP         | Whether to enable OpenMP     | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_CUDA           | Whether to enable CUDA       | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_HIP            | Whether to enable HIP        | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_SYCL           | Whether to enable SYCL       | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_NETCDF         | Whether to enable NETCDF     | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_PARTICLES      | Whether to enable particles  | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_MULTIBLOCK     | Whether to enable multiblock | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_RADIATION      | Whether to enable radiation  | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_SHOC           | Whether to enable shoc       | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_P3             | Whether to enable P3         | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_TESTS          | Whether to enable tests      | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_FCOMPARE       | Whether to enable fcompare   | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+

   .. note::
      **At most one of ERF_ENABLE_OMP, ERF_ENABLE_CUDA, ERF_ENABLE_HIP and ERF_ENABLE_SYCL should be set to true.**

CMake Logging Options
~~~~~~~~~~~~~~~~~~~~~

ERF's CMake configuration supports hierarchical logging to help diagnose build issues and understand configuration decisions. These options are available with CMake 3.25+:

Logging Levels
++++++++++++++

The ``--log-level`` flag controls the verbosity of CMake output:

.. code:: shell

   # Quiet output (STATUS messages only)
   cmake ..

   # Show detection details
   cmake --log-level=VERBOSE ..

   # Show all diagnostics
   cmake --log-level=DEBUG ..

Logging Context
+++++++++++++++

The ``--log-context`` flag shows the hierarchy of messages, making it easier to understand which component is reporting:

.. code:: shell

   # Show message hierarchy with component names
   cmake --log-context ..

   # Combine with verbose output for detailed diagnostics
   cmake --log-context --log-level=VERBOSE ..

Example output with ``--log-context``:

::

   [ERF.Cray] Detected Cray Programming Environment (CRAYPE_VERSION=2.7.30)
   [ERF.Cray] Setting Cray compiler wrappers...
   [ERF.Cray]   Set CMAKE_CXX_COMPILER = /opt/cray/pe/craype/default/bin/CC
   [ERF.AMReX] Using internal AMReX submodule
   [ERF.NetCDF] Found NetCDF: /opt/cray/pe/netcdf/4.9.0.9

CMake Utility Targets
~~~~~~~~~~~~~~~~~~~~~

ERF provides several utility targets to help manage builds:

Clean Build Artifacts
+++++++++++++++++++++

To perform a complete clean (equivalent to ``make distclean`` in GNU Make):

.. code:: shell

   make distclean

This removes all CMake configuration and build artifacts, including:

- CMake cache and generated files (``CMakeCache.txt``, ``CMakeFiles/``, etc.)
- Build outputs (executables, libraries)
- Generated configuration files
- Test outputs

The install directory is preserved.

Uninstall
+++++++++

To uninstall files that were installed via ``make install`` or ``cmake --install``:

.. code:: shell

   make uninstall

Show Detected Configuration
++++++++++++++++++++++++++++

On Cray systems, you can view the auto-detected configuration:

.. code:: shell

   make show-cray-config

This displays the configuration that was automatically detected and can be saved to use as a starting point for manual configuration.

Using Configuration Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

On Cray systems, ERF automatically detects the system configuration and applies necessary workarounds. If you need to use a manual configuration, you can create a configuration file and use it with CMake's ``-C`` option:

.. code:: shell

   cmake -C path/to/config.cmake ..

ERF provides machine-specific profile examples in ``Build/machines/``:

- ``perlmutter_erf.profile`` - NERSC Perlmutter (NVIDIA A100)
- ``frontier_erf.profile`` - OLCF Frontier (AMD MI250X)
- ``polaris_erf.profile`` - ALCF Polaris (NVIDIA A100)
- ``aurora_erf.profile`` - ALCF Aurora (Intel GPUs)

These profiles show the recommended modules to load for each system.

Mac with CMake
~~~~~~~~~~~~~~
Tested with macOS 12.7 (Monterey) using cmake (3.27.8), open-mpi (5.0.0), and
pkg-config (0.29.2) installed with the homebrew package manager.
(Note: The homebrew openmpi library was "poured from bottle," not installed
with the ``--build-from-source`` option.)
NetCDF will be compiled from source. The instructions below should be version
agnostic.

NetCDF (tested with v4.9.2)

#. Download latest source package from `ucar.edu`_
#. (Optional) install Zstd compression library ``brew install zstd``
#. Create build directory ``cd netcdf-c-4.9.2 && mkdir build && cd build``
#. Configure for your system ``../configure --enable-parallel CC=mpicc CXX=mpicxx LDFLAGS="-L/opt/homebrew/Cellar/zstd/1.5.5/lib" CPPFLAGS="-I/opt/homebrew/Cellar/zstd/1.5.5/include"``
   (omit the LDFLAGS and CPPFLAGS if you do not have Zstd installed) -- note
   that you may encounter cmake errors if you do not have pkg-config installed
#. Build ``make -j8`` and ``sudo make install``

.. _ucar.edu: https://downloads.unidata.ucar.edu/netcdf/

ERF (tested with commit ``40e64ed35ebc080ad61d08aea828330dfbdbc162``)

#. Get latest source code ``git clone --recursive git@github.com:erf-model/ERF.git``
#. Create build directory ``cd ERF && mkdir MyBuild && cd MyBuild``
#. Configure with cmake and build

::

    cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
       -DCMAKE_CXX_COMPILER:STRING=mpicxx \
       -DCMAKE_C_COMPILER:STRING=mpicc \
       -DCMAKE_Fortran_COMPILER:STRING=mpifort \
       -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
       -DERF_DIM:STRING=3 \
       -DERF_ENABLE_MPI:BOOL=ON \
       -DERF_ENABLE_TESTS:BOOL=ON \
       -DERF_ENABLE_FCOMPARE:BOOL=ON \
       -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
       -DERF_ENABLE_NETCDF:BOOL=ON \
       -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
       .. && make -j8

CMake Docs
~~~~~~~~~~
Building the ERF documentation with CMake can be completed by configuring with the flag ``-DERF_ENABLE_DOCUMENTATION:BOOL=ON`` and then compiling the following target:

.. code:: shell

   make docs

Note again that both the Sphinx and Doxygen documentation will be built.

Cray Systems
~~~~~~~~~~~~

ERF supports compilation on Cray systems with both automatic detection and manual configuration.

.. note::

   The machine profiles and some of the build documentation for HPC systems are inspired by and modeled after the excellent `WarpX HPC documentation <https://warpx.readthedocs.io/en/latest/install/hpc.html>`__.

**Auto-Detection Control**

The build system can automatically detect Cray environments. Control this with:

- **Enable auto-detection** (default on Cray systems): ``-DCRAY_AUTO_DETECTION=ON``
- **Disable auto-detection**: ``-DCRAY_AUTO_DETECTION=OFF``

**Using Cray CMake Configuration Files**

If you have a Cray CMake configuration file (generated by your system or provided), use it with the ``-C`` flag:

.. code:: shell

   cmake -C /path/to/cray_config.cmake \
         -DERF_ENABLE_CUDA:BOOL=ON \
         -DERF_ENABLE_MPI:BOOL=ON \
         ..

**Generating a Cray Configuration File**

You can generate a template Cray CMake configuration file by running:

.. code:: shell

   make show-cray-config

This will show a file ``cray_config.cmake`` in your build directory with detected Cray settings. You can then:

1. Review the generated settings
2. Make any necessary adjustments for your specific system
3. Use it with ``cmake -C cray_config.cmake ...``

The generated file will include settings such as:

.. code:: cmake

   set(CMAKE_C_COMPILER "cc")
   set(CMAKE_CXX_COMPILER "CC")
   set(CMAKE_Fortran_COMPILER "ftn")
   set(MPI_C_COMPILER "cc")
   set(MPI_CXX_COMPILER "CC")
   # ... and other detected settings

**How Auto-Detection Works**

When ``CRAY_AUTO_DETECTION=ON``, the build system performs the following detection and configuration:

**Compiler Detection:**

- First checks if ``CRAYPE_DIR`` environment variable is set (indicates Cray environment)
- Sets compilers to Cray wrappers: ``cc``, ``CC``, and ``ftn``
- Fallback: If Cray wrappers are not found, uses system default compilers (typically ``gcc``, ``g++``, ``gfortran``)
- Decision: If using Cray wrappers, automatically sets ``-DCMAKE_CXX_COMPILER_WRAPPER=CrayPrgEnv``

**MPI Configuration:**

- Automatically sets MPI compilers to match the Cray wrappers (``cc`` and ``CC``)
- Checks for ``CRAY_MPICH_DIR`` environment variable to locate MPI libraries
- Fallback: Uses standard MPI detection via ``FindMPI`` if Cray MPICH is not found
- Decision: Prefers Cray MPICH over other MPI implementations when detected

**CUDA Configuration:**

- When ``ERF_ENABLE_CUDA=ON``, sets ``CMAKE_CUDA_HOST_COMPILER`` to the Cray C++ wrapper (``CC``)
- Looks for ``nvcc`` in standard locations and the ``CUDA_HOME`` environment variable
- Fallback: If ``nvcc`` is not in PATH, checks ``/usr/local/cuda/bin``
- Decision: Always uses the Cray wrapper as the CUDA host compiler to ensure compatibility with the Cray programming environment

**Linker Configuration:**

- Detects the default linker used by the Cray wrappers
- On newer systems, may default to ``lld`` (LLVM linker) or the GNU linker
- Fallback: Uses the compiler's default linker if detection fails
- Override: Can be manually set with ``-DCMAKE_LINKER_TYPE=`` or ``-DCMAKE_LINKER_EXE=``

**Module System Integration:**

- Checks for loaded Cray modules (e.g., ``PrgEnv-gnu``, ``PrgEnv-cray``, ``cudatoolkit``)
- Uses module environment variables (``CRAY_*_DIR``) to locate libraries
- Warning: Will warn if conflicting programming environments are loaded

**Common Workarounds**

- **MPI not found**: Manually set ``-DMPI_HOME=$CRAY_MPICH_DIR``
- **CUDA math libraries not found**: Load the latest cmake module with ``module load cmake``
- **Compiler version mismatch**: Verify module compatibility with ``module list`` and ensure only one programming environment is loaded
- **Linker errors**: Explicitly set the linker with ``-DCMAKE_LINKER_EXE=/path/to/ld``

Build Scripts Reference
~~~~~~~~~~~~~~~~~~~~~~~~

ERF provides several build scripts optimized for different systems and architectures. This table shows which scripts have been tested and verified on each system. Verified builds are marked with the git commit hash where they were last tested.

.. list-table:: ERF Build Scripts by System
   :header-rows: 1
   :widths: 35 10 10 10 10 10 10 10

   * - Build Script
     - Perlmutter
     - Frontier
     - Aurora
     - Polaris
     - Kestrel
     - RegtestCPU
     - RegtestGPU
   * - ``cmake.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
   * - ``cmake_with_kokkos_many.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
   * - ``cmake_with_kokkos_many_cuda.sh``
     - Untested
     - —
     - —
     - Untested
     - Untested
     - —
     - Untested
   * - ``cmake_with_kokkos_many_noradiation_hip.sh``
     - —
     - Untested
     - —
     - —
     - —
     - —
     - —
   * - ``cmake_with_kokkos_many_sycl.sh``
     - —
     - —
     - Untested
     - —
     - —
     - —
     - —
   * - ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh``
     - Untested
     - Untested
     - Untested
     - Untested
     - Untested
     - —
     - Untested

.. note::

   The ``build_erf_with_shoc_cuda_Perlmutter.sh`` script is being tested for cross-site compatibility with auto-detection enabled. A simplified version may work across CUDA-enabled HPC sites (Perlmutter, Polaris, Kestrel) with ``-DCRAY_AUTO_DETECTION=ON``.

**Script Descriptions**

- ``cmake.sh``: Basic CMake build script, works on all systems
- ``cmake_with_kokkos_many.sh``: Kokkos-enabled build for CPU architectures
- ``cmake_with_kokkos_many_cuda.sh``: Kokkos build optimized for NVIDIA GPUs
- ``cmake_with_kokkos_many_noradiation_hip.sh``: Kokkos build for AMD GPUs (HIP backend)
- ``cmake_with_kokkos_many_sycl.sh``: Kokkos build for Intel GPUs (SYCL backend)
- ``Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh``: SHOC turbulence model with CUDA support

**Testing Build Scripts**

Build scripts can be tested using the configuration file ``Build/ERF-cmake-tests.ini``, which defines test cases for each system and script combination.

Perlmutter (NERSC)
~~~~~~~~~~~~~~~~~~

Recall the GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

For Perlmutter at NERSC, look at the general instructions for building ERF using GNU Make, and then you can initialize your environment by loading these modules:

::

   module load PrgEnv-gnu
   module load cudatoolkit

If you will be using NetCDF, we suggest you add the following four lines to your ``.bashrc_ext`` file.

::

   module load cray-hdf5-parallel/1.12.2.9
   module load cray-netcdf-hdf5parallel
   export NETCDF_DIR=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9

Then build ERF as, for example (specify your own path to the AMReX submodule in ``ERF/Submodules/AMReX``):

::

   make -j 4 COMP=gnu USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=TRUE AMREX_HOME=/path_to_here/ERF/Submodules/AMReX

where ``/path_to_here`` is the path to your ERF repository.

**Using the Provided CMake Build Scripts**

If you prefer to use one of the provided CMake build scripts, first load the appropriate machine profile:

::

   source Build/machines/perlmutter_erf.profile

Or load the modules manually:

::

   module load gcc/12.2.0 cmake cudatoolkit cray-mpich cray-hdf5-parallel cray-netcdf-hdf5parallel

Note that the latest cmake module helps with finding the CUDA math libraries. The gcc version should be a gcc variant such as ``gcc/12.2.0`` that has been tested to reliably find the filesystem functionality used in some utilities. Also note that the ``cray-hdf5-parallel`` module matches with the parallel configuration requirements.

Then build with one of the provided build scripts. We suggest using ``Build/cmake_cuda.sh`` for a basic CUDA build, or ``Build/cmake_with_netcdf.sh`` if you need NetCDF support.
With the flags listed explicitly instead of using -DCRAY_AUTO_DECTION=ON, we suggest using ``Build/Perlmutter/cmake_with_cuda_perlmutter.sh`` for a basic CUDA build, or ``Build/Perlmutter/build_erf_with_shoc_cuda_Perlmutter.sh `` if you need SHOC support.

To use logging for diagnostics during configuration:

::

   cmake --log-context --log-level=VERBOSE \
         -DERF_ENABLE_CUDA:BOOL=ON \
         -DERF_ENABLE_MPI:BOOL=ON \
         .. && make

**Submitting Jobs**

Finally, you can prepare your SLURM job script, using the following as a guide:

   .. code:: shell

             #!/bin/bash

             ## specify your allocation (with the _g) and that you want GPU nodes
             #SBATCH -A m4106_g
             #SBATCH -C gpu

             ## the job will be named "ERF" in the queue and will save stdout to erf_[job ID].out
             #SBATCH -J ERF
             #SBATCH -o erf_%j.out

             ## set the max walltime
             #SBATCH -t 10

             ## specify the number of nodes you want
             #SBATCH -N 2

             ## we use the same number of MPI ranks per node as GPUs per node
             #SBATCH --ntasks-per-node=4
             #SBATCH --gpus-per-node=4
             #SBATCH --gpu-bind=none

             # pin to closest NIC to GPU
             export MPICH_OFI_NIC_POLICY=GPU

             # use GPU-aware MPI
             #GPU_AWARE_MPI=""
             GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

             # the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)
             # set ordering of CUDA visible devices inverse to local task IDs for optimal GPU-aware MPI
             srun -n 8 --cpus-per-task=32 --cpu-bind=cores bash -c "
               export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
               ./ERF3d.gnu.MPI.CUDA.ex inputs_wrf_baseline max_step=100 ${GPU_AWARE_MPI}" \
             > test.out

To submit your job script, do ``sbatch [your job script]`` and you can check its status by doing ``squeue -u [your username]``.

Kestrel (NREL)
~~~~~~~~~~~~~~

The `Kestrel <https://nrel.github.io/HPC/Documentation/Systems/Kestrel/>`_ cluster is an HPE Cray machine
composed primarily of CPU compute nodes with 104 core
Intel Xeon Sapphire Rapids nodes. It also contains a GPU partition with 4 Nvidia H100 GPUs per node.

As with Perlmutter, the GNU Make build system is preferred. To compile and run on CPUs, the default modules
loaded when logging into Kestrel can be used. If you are unsure about your environment, you can reset to
the default modules: ::

  module restore

Then, build ERF using the cray compilers (if wishing to use other compilers, you can swap the ``PrgEnv-cray`` module
for another module as appropriate, see Kestrel user documentation for more details): ::

  make realclean; make -j COMP=cray

To run on GPUs on Kestrel, note that the machine has separate login nodes for GPU use and GPU jobs should only
be started from GPU login nodes (accessed via ``kestrel-gpu.hpc.nrel.gov``). For compiling and running on GPUs,
the following commands can be used to set up your environment: ::

  module purge;
  module load PrgEnv-gnu/8.5.0;
  module load cuda/12.3;
  module load craype-x86-milan;

And then compile (for example, in ``ERF/Exec/ABL``): ::

  make realclean; make -j COMP=gnu USE_CUDA=TRUE

As a word of warning, system updates on Kestrel periodically change the necessary modules that must be loaded
in order to build and run ERF, so these instructions may become out of date.

When running on Kestrel, GPU node hours are charged allocation units (AUs) at 10 times the rate of CPU node hours.
For ERF, the performance running on a Kestrel GPU node with 4 GPUs is typically 10-20x running on a CPU node
with 96-104 MPI ranks per node, so the performance gain from on on GPUs is likely worth the higher charge
rate for node hours, in addition to providing faster time to solution. However, for smaller problem sizes,
or problems distributed across too many nodes (resulting in fewer than around 1 million cells/GPU),
the compute capability of the GPUs may be unsaturated and the performance gain from running on GPUs
may not justify the higher AU charge. The trade-off is problem dependent, so users may wish to assess
performance for their particular case and objectives in terms of wall time, AUs used, etc to determine the
optimal strategy if running large jobs.

Another note about using Kestrel is that partial node allocations are possible, which means the full memory
available on each node may not be assigned by default. In general, using the ``--exclusive`` flag when
requesting nodes through the slurm scheduler, which will allocate entire nodes exlcusively for your request,
is recommended. Otherwise, memory intensive operations such as CUDA compilation may fail. You can alternatively
request a particular amount of memory with the ``--mem=XXX`` or ``--mem-per-cpu=XXX`` slurm inputs.
