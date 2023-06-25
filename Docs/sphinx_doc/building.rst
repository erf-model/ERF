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

#. ``cd`` to the desired build directory, e.g.  ``ERF/Exec/IsentropicVortex/``

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
   | USE_MOISTURE       | Whether to enable moisture   | TRUE / FALSE     | FALSE       |
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
   For example, the default build in ``ERF/Exec/Isntropic`` will look
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
   | ERF_ENABLE_MOISTURE       | Whether to enable moisture   | TRUE / FALSE     | FALSE       |
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



Perlmutter (NERSC)
~~~~~~~~~~~~~~~~~~

Recall the GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

For Perlmutter at NERSC, look at the general instructions for building ERF using GNU Make, and then you can initialize your environment by loading these modules:

::

   module load PrgEnv-gnu
   module load cudatoolkit

Then build ERF as, for example (specify your own path to the AMReX submodule in `ERF/Submodules/AMReX`):

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

To submit your job script, do `sbatch [your job script]` and you can check its status by doing `squeue -u [your username]`.

