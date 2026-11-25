This problem setup tests the build of the radiation capability.

# Compilation on CPU

To build with RRTMGP using gmake, type "make" in a directory whose GNUmakefile
contains the lines
```
   USE_RRTMGP = TRUE
   USE_NETCDF = TRUE
```
The first build also builds the Kokkos submodule (`Submodules/kokkos`) and
installs it under `Submodules/kokkos/build/<config>`; later builds with the
same configuration reuse that install. To use an existing Kokkos install
instead, set `ERF_KOKKOS_HOME` to it (`Build/GNU_Kokkos/kokkos_build_commands.sh`
makes one).

# Compilation on GPUs on Perlmutter

As on CPU, with `USE_CUDA = TRUE` in the GNUmakefile; Kokkos is built for
CUDA automatically. To build a shared Kokkos install once instead:
```
    cd ERF
    source Build/GNU_Kokkos/kokkos_build_cuda_commands_Perlmutter.sh
```
and then build with that `ERF_KOKKOS_HOME` in the environment.
