 .. role:: cpp(code)
    :language: c++

.. _containers:

ERF containers on Perlmutter
============================

The following is a brief introduction to ERF and containers on the NERSC Perlmutter platform.

For more details, please see `NERSC's detailed containers documentation <https://docs.nersc.gov/development/containers>`_, which also includes containers tutorials.

Container images can be built on one's desktop/laptop using a standard container framework such as Docker, Podman (Pod Manager), etc. or directly on a Perlmutter login node using ``podman-hpc``.  ``podman-hpc`` is a NERSC-developed wrapper that extends the capabilities of Podman for HPC.  Containers are run on Perlmutter using either ``podman-hpc`` or ``shifter`` (also developed at NERSC).  NERSC has a good `podman-hpc tutorial <https://docs.nersc.gov/development/containers/podman-hpc/podman-beginner-tutorial>`_.

Example ERF containerfile
~~~~~~~~~~~~~~~~~~~~~~~~~

The following is an example ERF containerfile with filename ``erf_containerfile``:

.. code:: shell

     1    FROM nvcr.io/nvidia/cuda:12.2.0-devel-ubuntu22.04
     2
     3    WORKDIR /app
     4    ARG base_dir=/app/erf
     5
     6    RUN apt-get update -y && \
     7        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
     8        g++-12 \
     9        gcc-12 \
    10        gfortran-12 \
    11        git \
    12        libtool \
    13        make \
    14        tar \
    15        autoconf \
    16        automake \
    17        wget \
    18        python3 \
    19        cmake
    20
    21    # MPICH to be swapped out later for Cray MPI
    22    ARG mpich_version=4.2.2
    23    ARG mpich_dir=mpich-${mpich_version}
    24
    25    RUN wget https://www.mpich.org/static/downloads/$mpich_version/$mpich_dir.tar.gz && \
    26        tar xzf $mpich_dir.tar.gz && \
    27        cd $mpich_dir && \
    28        ./configure CC=/usr/bin/gcc-12 CXX=/usr/bin/g++-12 F77=/usr/bin/gfortran-12 FC=/usr/bin/gfortran-12 && \
    29        make -j8 && \
    30        make install && \
    31        make clean  && \
    32        cd ..  && \
    33        rm -rf $mpich_dir $mpich_dir.tar.gz
    34
    35    RUN mkdir ${base_dir}
    36
    37    ARG build_dir=MyBuild
    38
    39    RUN cd ${base_dir} && git clone --recursive https://github.com/erf-model/ERF.git && \
    40      cd ERF && mkdir ${build_dir} && cd ${build_dir} && \
    41      cmake \
    42      -DCMAKE_C_COMPILER=mpicc \
    43      -DCMAKE_CXX_COMPILER=mpicxx \
    44      -DCMAKE_Fortran_COMPILER=mpif90 \
    45      -DCMAKE_BUILD_TYPE:STRING=Release \
    46      -DCMAKE_CUDA_ARCHITECTURES=80 \
    47      -DERF_ENABLE_MPI:BOOL=ON \
    48      -DERF_ENABLE_CUDA=ON \
    49      .. && \
    50      make -j8

Line numbers were added for instructional purposes but should not appear in containerfile.
The containerfile is available at https://github.com/erf-model/ERF/blob/development/Build/erf_containerfile

* Line 1 downloads a container base image from NVIDIA's container registry that contains the Ubuntu 22.04 operating system and CUDA 12.2.0
* Line 3 sets the working directory
* Line 4 sets a value for the variable ``base_dir``
* Lines 6-19 download the GNU 12 compilers and all the necessary utilities for building ERF in the container
* Lines 21-33 download the MPICH source code and builds it.  At runtime this MPICH will get replaced by Perlmutter's Cray MPI
* Lines 39-50 clone the ERF repo and build it with cmake/make



Build the container on Perlmutter using ``podman-hpc`` and using the containerfile ``erf_containerfile`` with name ``erf`` and tag ``1.00`` (``-t <name>:<tag>``)

   .. code:: shell

    podman-hpc build -t erf:1.00 -f erf_containerfile

In order to use this image in a job or access it from any other login node, one needs to migrate the image onto the $SCRATCH filesystem by issuing the following command:

.. code:: shell

  podman-hpc migrate erf:1.00

``podman-hpc images`` will display the following

.. code:: shell

  user@perlmutter:login12:~> podman-hpc images
  REPOSITORY                    TAG                         IMAGE ID      CREATED       SIZE        R/O
  localhost/erf                 1.00                        893605c3ee9b  5 hours ago   16.1 GB     true
  localhost/erf                 1.00                        893605c3ee9b  5 hours ago   16.1 GB     false

Note that ``localhost`` will not be needed for podman-hpc commands.

Run container on Perlmutter
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Submit the following slurm batch script in order to use the image in a job

.. code:: shell

  #!/bin/bash

  #SBATCH --account=<proj>
  #SBATCH --constraint=gpu
  #SBATCH --job-name=erf
  #SBATCH --nodes=1
  #SBATCH --time=0:05:00
  #SBATCH -q regular

  srun -N 1 -n 4 -c 32 --ntasks-per-node=4 --gpus-per-node=4 ./device_wrapper \
  podman-hpc run --rm --mpi --gpu -v /pscratch/sd/u/user/erf/abl:/run -w /run erf:1.00 \
  /app/erf/ERF/MyBuild/Exec/ABL/erf_abl inputs_smagorinsky amrex.use_gpu_aware_mpi=0

``device_wrapper`` script:

.. code:: shell

      #!/bin/bash
      # select_cpu_device wrapper script
      export CUDA_VISIBLE_DEVICES=$((3-$SLURM_LOCALID))
      exec $*

Arguments for ``podman-hpc run`` used above

* ``--rm`` removes the container after exit
* ``--mpi`` enables Cray MPI support (swaps MPICH in the container for Perlmutter's Cray MPI)
* ``--gpu`` enables NVIDIA GPU support
* ``-v /pscratch/sd/u/user/erf/abl:/run`` mounts ``/pscratch/sd/u/user/erf/abl`` on Perlmutter onto ``/run`` in the container
* ``-w /run`` makes the ``/run`` directory inside the container the working directory, i.e. any output from the ERF run will be written to the ``/run`` directory in the container, which will appear in the ``/pscratch/sd/u/user/erf/abl`` directory on Perlmutter.
* ``erf:1.00`` container name and tag
* ``/app/erf/ERF/MyBuild/Exec/ABL/erf_abl`` ERF binary in container

The remaining arguments are the normal ERF command line arguments.

Please issue ``podman-hpc --help`` for the help page and ``podman-hpc run --help`` for the ``podman-hpc run`` help page.

Container image libraries
~~~~~~~~~~~~~~~~~~~~~~~~~

Container image libraries provide a convenient way to store and share images.
The best known one is probably Docker Hub.  NERSC provides a private registry to its users via `registry.nersc.gov <https://docs.nersc.gov/development/containers/registry>`_.

In order to push the ERF image to Docker Hub, need to modify the tag, login to Docker Hub, and then push the image.

.. code:: shell

   podman-hpc tag erf:1.00 docker.io/docker_hub_username/erf:1.00
   podman-hpc login docker.io
   podman-hpc push docker.io/docker_hub_username/erf:1.00

Shifter container runtime
~~~~~~~~~~~~~~~~~~~~~~~~~

Shifter is a NERSC-developed tool that provides an alternative method for running containers on Perlmutter.   Recall that Shifter cannot be used to build images. `NERSC's containers documentation <https://docs.nersc.gov/development/containers>`_ provides an introduction to shifter including a tutorial.

Use the ``shifterimg pull`` command to pull images directly from Docker Hub and it will automatically convert your Docker image into Shifter format.

Need to first login to Docker Hub, then pull.  ``--user nersc_user`` will restrict usage of the image to only ``nersc_user``.

.. code:: shell

   shifterimg login docker.io
   shifterimg -v --user nersc_user pull docker_hub_username/erf:1.00

Submit the following slurm batch script in order to use the image in a job

.. code:: shell

  #!/bin/bash

  #SBATCH --account=<proj>
  #SBATCH --constraint=gpu
  #SBATCH --job-name=erf
  #SBATCH --nodes=1
  #SBATCH --time=0:05:00
  #SBATCH -q regular
  #SBATCH --image=docker_hub_username/erf:1.00 # for shifter, not podman-hpc
  #SBATCH --module=mpich,gpu  # for shifter, not podman-hpc; CPU-only MPI
  ##SBATCH --module=cuda-mpich  # for shifter, not podman-hpc; GPU-Aware MPI
  #SBATCH --volume="<$SCRATCH/my_erf_run_dir>:/run"  # for shifter, not podman-hpc

  srun -N 1 -n 4 -c 32 --ntasks-per-node=4 --gpus-per-node=4 ./device_wrapper \
  shifter \
  /app/erf/ERF/MyBuild/Exec/ABL/erf_abl inputs_smagorinsky amrex.use_gpu_aware_mpi=0



Common Issues
~~~~~~~~~~~~~

* Using ``podman`` rather than ``podman-hpc`` on Perlmutter (best to always use ``podman-hpc``)
* Before issuing ``podman-hpc migrate <name>:<tag>`` after having issued the command earlier with identical ``<name>:<tag>``, if want to keep the same name, please change the ``<tag>`` to one that has not been used previously.  If want to use an identical ``<name>:<tag>`` used in a previous ``podman-hpc migrate`` command, please first issue ``podman-hpc rmsqi <name>:<tag>`` to delete the old image.  Otherwise could potentially end up with errors such as

    .. code:: shell

      Error: read-only image store assigns the same name to multiple images

  and will have resort to `various methods <https://docs.nersc.gov/development/containers/podman-hpc/overview/#troubleshooting>`_ to get out of the bad configuration state.

* The default containerfile is a file called ``Containerfile`` (case sensitive).  When that file is being used, can replace ``-f erf_containerfile`` with a period:

   .. code:: shell

    podman-hpc build -t erf:1.00 .

  Note that for this case, the period is mandatory. Here it does not denote the end of a sentence.
