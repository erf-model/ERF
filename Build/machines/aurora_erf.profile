# required dependencies
module load cmake

module load hdf5/1.14.6
module load netcdf-cxx4

# necessary to use build or run with GPU-aware MPICH
# export MPIR_CVAR_ENABLE_GPU=1