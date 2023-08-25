#Set amrex options
set(USE_XSDK_DEFAULTS OFF)
set(AMReX_SPACEDIM "${ERF_DIM}" CACHE STRING "Number of physical dimensions" FORCE)
set(AMReX_PIC OFF)
set(AMReX_MPI ${ERF_ENABLE_MPI})
set(AMReX_OMP ${ERF_ENABLE_OPENMP})
set(AMReX_PRECISION "${ERF_PRECISION}" CACHE STRING "Floating point precision" FORCE)
set(AMReX_EB OFF)
set(AMReX_FORTRAN_INTERFACES OFF)
set(AMReX_AMRDATA OFF)
set(AMReX_SENSEI OFF)
set(AMReX_CONDUIT OFF)
set(AMReX_HDF5 "${ERF_ENABLE_HDF5}" CACHE STRING "Enable HDF5 support in AMReX" FORCE)
set(AMReX_ENABLE_TESTS OFF)
set(AMReX_FPE OFF)
set(AMReX_ASSERTIONS OFF)
set(AMReX_BASE_PROFILE OFF)
set(AMReX_TINY_PROFILE ${ERF_ENABLE_TINY_PROFILE})
set(AMReX_TRACE_PROFILE OFF)
set(AMReX_MEM_PROFILE OFF)
set(AMReX_COMM_PROFILE OFF)
set(AMReX_BACKTRACE OFF)
set(AMReX_PROFPARSER OFF)
set(AMReX_CUDA ${ERF_ENABLE_CUDA})
set(AMReX_ACC OFF)
set(AMReX_PLOTFILE_TOOLS ${ERF_ENABLE_FCOMPARE})
set(AMReX_FORTRAN OFF)

set(AMReX_LINEAR_SOLVERS OFF)
if(ERF_ENABLE_POISSON_SOLVE)
  set(AMReX_LINEAR_SOLVERS ON)
endif()

set(AMReX_PARTICLES OFF)
if(ERF_ENABLE_PARTICLES)
  set(AMReX_PARTICLES ON)
endif()

if(ERF_ENABLE_CUDA)
  set(AMReX_GPU_BACKEND CUDA CACHE STRING "AMReX GPU type" FORCE)
  set(AMReX_CUDA_WARN_CAPTURE_THIS OFF)
  set(AMReX_CUDA_ERROR_CAPTURE_THIS ON)
endif()

if(ERF_ENABLE_HIP)
  set(AMReX_GPU_BACKEND HIP CACHE STRING "AMReX GPU type" FORCE)
endif()

if(ERF_ENABLE_SYCL)
  set(AMReX_GPU_BACKEND SYCL CACHE STRING "AMReX GPU type" FORCE)
endif()
