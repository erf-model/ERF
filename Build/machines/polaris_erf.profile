# swap to the Milan cray package
module load craype-x86-milan

# extra modules
module use /soft/modulefiles
module load spack-pe-gnu

# add cuda
module load cuda/12.6
module load cudatoolkit-standalone/12.6
module load craype-accel-nvidia80

# required dependencies
module load cmake

# default gcc-native too new for cuda/12.6
module load gcc-native/13.2

module load cray-hdf5-parallel
module load cray-libsci/25.03.0


module load cray-netcdf-hdf5parallel

# export MPICH_GPU_SUPPORT_ENABLED=1