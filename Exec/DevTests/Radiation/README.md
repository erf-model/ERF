This problem setup tests the build of the radiation capability.

# Compilation on CPU

To build with RRTMGP using gmake:
```
   export ERF_DIR=/path/to/ERF
   source /path/to/ERF/Build/GNU_Ekat/ekat_build_commands.sh
```

(Note that this only needs to be done once after the clone of this repo.)

Then type "make" here, after verifying that the GNUmakfile
contains the lines
```
   USE_RRTMGP = TRUE
   USE_NETCDF = TRUE
```
# Compilation on GPUs on Perlmutter

```
    cd ERF
    source Build/GNU_Ekat/ekat_build_cuda_commands_Perlmutter.sh
    cd Exec/DevTests/Radiation
    cp GNUmakefile_GPU GNUmakefile
    make -j8
```
