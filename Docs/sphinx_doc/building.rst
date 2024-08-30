.. _Building:

Building
--------

The ERF code is dependent on AMReX, and uses the radiation model (RTE-RRTMGP) which is based on YAKL C++ implementation for heterogeneous computing infrastructure (which are all available as submodules in the ERF repo). ERF can be built using either GNU Make or CMake, however, if radiation model is activated, only CMake build system is supported.

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
boundary layer or ``DensityCurrent`` for running the standard density current test case, etc, 2) regression tests in
``Exec/RegTests`` that are used for testing specific known aspects of the code functionality, such as boundary conditions or
Rayleigh damping, and 3) tests for features under development in ``Exec/DevTests``, such as moving terrain.  There is a
README in each problem directory that describes the purpose/role of that problem.

GNU Make
~~~~~~~~

The GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

Using the GNU Make build system involves first setting environment variables for the directories of the dependencies of ERF (AMReX, RTE-RRTMGP, and YAKL); note, RTE-RRTMGP, and YAKL are only required if running with radiation. All dependencies are provided as git submodules in ERF and can be populated by using ``git submodule init; git submodule update`` in the ERF repo, or before cloning by using ``git clone --recursive <erf_repo>``. Although submodules of these projects are provided, they can be placed externally as long as the ``<REPO_HOME>`` environment variables for each dependency is set correctly. An example of setting the ``<REPO_HOME>`` environment variables in the user's ``.bashrc`` is shown below:

::

   export ERF_HOME=${HOME}/ERF
   export AMREX_HOME=${ERF_HOME}/Submodules/AMReX

The GNU Make system is set up to use the path to AMReX submodule by default, so it is not necessary to set
these paths explicitly, unless it is desired to do so. It is also possible to use an external version of
AMReX, downloaded by running

   .. code:: shell

             git clone https://github.com/amrex-codes/amrex.git

in which case the ``AMREX_HOME`` environment variable must point to the location where AMReX has been downloaded, which will take precedence over the default path to the submodule. If using bash shell,

::

   export AMREX_HOME=/path/to/external/amrex

or if using tcsh,

::

   setenv AMREX_HOME /path/to/external/amrex

#. ``cd`` to the desired build directory, e.g.  ``ERF/Exec/RegTests/IsentropicVortex/``

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
   | USE_HDF5           | Whether to enable HDF5       | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_PARTICLES      | Whether to enable particles  | TRUE / FALSE     | FALSE       |
   +--------------------+------------------------------+------------------+-------------+
   | USE_WARM_NO_PRECIP | Whether to use warm moisture | TRUE / FALSE     | FALSE       |
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
      **Do not set both USE_OMP and USE_CUDA to true.**

   Information on using other compilers can be found in the AMReX documentation at
   https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html .

#. Make the executable by typing

   .. code:: shell

      make

   The name of the resulting executable (generated by the GNUmake system) encodes several of the build characteristics, including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable.
   Thus, several different build configurations may coexist simultaneously in a problem folder.
   For example, the default build in ``ERF/Exec/RegTests/IsentropicVortex`` will look
   like ``ERF3d.gnu.MPI.ex``, indicating that this is a 3-d version of the code, made with
   ``COMP=gnu``, and ``USE_MPI=TRUE``.

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
   | ERF_ENABLE_HDF5           | Whether to enable HDF5       | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_PARTICLES      | Whether to enable particles  | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_WARM_NO_PRECIP | Whether to use warm moisture | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_MULTIBLOCK     | Whether to enable multiblock | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_RADIATION      | Whether to enable radiation  | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_TESTS          | Whether to enable tests      | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+
   | ERF_ENABLE_FCOMPARE       | Whether to enable fcompare   | TRUE / FALSE     | FALSE       |
   +---------------------------+------------------------------+------------------+-------------+


Mac with CMake
~~~~~~~~~~~~~~
Tested with macOS 12.7 (Monterey) using cmake (3.27.8), open-mpi (5.0.0), and
pkg-config (0.29.2) installed with the homebrew package manager. HDF5 and
NetCDF will be compiled from source. The instructions below should be version
agnostic.

HDF5 (tested with v1.14.3)

#. Download latest source package from `hdfgroup.org`_
#. Extract source code ``tar xzf hdf5-<version>.tar.gz``
#. Create build directory ``cd hdf5-<version> && mkdir build && cd build``
#. Configure for your system ``../configure --prefix=/usr/local --enable-parallel``
#. Build ``make -j8`` and ``sudo make install``

.. _hdfgroup.org: https://www.hdfgroup.org/downloads/hdf5/source-code/

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
       -DERF_ENABLE_HDF5:BOOL=ON \
       -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
       .. && make -j8

Perlmutter (NERSC)
~~~~~~~~~~~~~~~~~~

Recall the GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

For Perlmutter at NERSC, look at the general instructions for building ERF using GNU Make, and then you can initialize your environment by loading these modules:

::

   module load PrgEnv-gnu
   module load cudatoolkit

Then build ERF as, for example (specify your own path to the AMReX submodule in ``ERF/Submodules/AMReX``):

::

   make -j 4 COMP=gnu USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=TRUE AMREX_HOME=/global/u2/d/dwillcox/dev-erf/ERF/Submodules/AMReX

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
