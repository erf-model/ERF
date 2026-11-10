/**
 * \file ERF_BuildDescription.cpp
 *
 * Compile-time capabilities for --describe-json and --describe-yaml, written
 * before AMReX or any device is initialized. The output depends only on the
 * build configuration (no date, host, paths, or flags); unknown values are null.
 */
#include <ERF.H>
#include <AMReX_buildInfo.H>
#include <ERF_NoahMPGitHash.H>
#include <ERF_DescriptionWriter.H>
#include <climits>

#ifdef ERF_USE_NETCDF
#include <netcdf_meta.h>
#endif

#ifdef ERF_USE_FFT
#include <AMReX_FFT_Helper.H>
#endif

#define ERF_DESCRIBE_STR_(x) #x
#define ERF_DESCRIBE_STR(x) ERF_DESCRIBE_STR_(x)

namespace {

const char* fft_library_version ()
{
#if !defined(ERF_USE_FFT)
    return nullptr;
#elif defined(AMREX_USE_CUDA)
#  if defined(CUFFT_VER_MAJOR) && defined(CUFFT_VER_MINOR) && defined(CUFFT_VER_PATCH)
    return ERF_DESCRIBE_STR(CUFFT_VER_MAJOR) "." ERF_DESCRIBE_STR(CUFFT_VER_MINOR) "."
           ERF_DESCRIBE_STR(CUFFT_VER_PATCH);
#  else
    return nullptr;
#  endif
#elif defined(AMREX_USE_HIP)
#  if defined(rocfft_version_major) && defined(rocfft_version_minor) && defined(rocfft_version_patch)
    return ERF_DESCRIBE_STR(rocfft_version_major) "." ERF_DESCRIBE_STR(rocfft_version_minor) "."
           ERF_DESCRIBE_STR(rocfft_version_patch);
#  else
    return nullptr;
#  endif
#elif defined(AMREX_USE_SYCL)
    return nullptr; // oneMKL's DFT headers carry no version
#elif defined(AMREX_USE_FLOAT)
    return fftwf_version;
#else
    return fftw_version;
#endif
}

void write_description (std::ostream& os, erf_description::Format format)
{
    erf_description::Writer writer(os, format);
    writer.value("schema_version", 1);
    writer.value("application", "ERF");
    writer.optional_string("source_describe", amrex::buildInfoGetGitHash(1));
    writer.value("dimension", AMREX_SPACEDIM); // always 3 for ERF
    writer.value("real_bits", static_cast<int>(sizeof(amrex::Real) * CHAR_BIT));
#ifdef ERF_USE_PARTICLES
    writer.value("particle_real_bits", static_cast<int>(sizeof(amrex::ParticleReal) * CHAR_BIT));
#else
    writer.value("particle_real_bits", nullptr);
#endif
#if defined(AMREX_USE_CUDA)
    writer.value("gpu_backend", "cuda");
#elif defined(AMREX_USE_HIP)
    writer.value("gpu_backend", "hip");
#elif defined(AMREX_USE_SYCL)
    writer.value("gpu_backend", "sycl");
#else
    writer.value("gpu_backend", "none");
#endif

    writer.begin_object("features");
#ifdef AMREX_USE_MPI
    writer.value("mpi", true);
#else
    writer.value("mpi", false);
#endif
#ifdef AMREX_USE_OMP
    writer.value("openmp", true);
#else
    writer.value("openmp", false);
#endif
#ifdef ERF_USE_KOKKOS
    writer.value("kokkos", true);
#else
    writer.value("kokkos", false);
#endif
#ifdef ERF_USE_PARTICLES
    writer.value("particles", true);
#else
    writer.value("particles", false);
#endif
#ifdef ERF_USE_FFT
    writer.value("fft", true);
#else
    writer.value("fft", false);
#endif
#ifdef ERF_USE_NETCDF
    writer.value("netcdf", true);
#else
    writer.value("netcdf", false);
#endif
#ifdef AMREX_USE_HDF5
    writer.value("hdf5", true);
#else
    writer.value("hdf5", false);
#endif
#ifdef ERF_USE_NOAHMP
    writer.value("noahmp", true);
#else
    writer.value("noahmp", false);
#endif
#ifdef ERF_USE_RRTMGP
    writer.value("rrtmgp", true);
#else
    writer.value("rrtmgp", false);
#endif
#ifdef ERF_USE_EAMXX_SHOC
    writer.value("eamxx_shoc", true);
#else
    writer.value("eamxx_shoc", false);
#endif
#ifdef ERF_USE_WINDFARM
    writer.value("windfarm", true);
#else
    writer.value("windfarm", false);
#endif
#ifdef ERF_USE_FIRE
    writer.value("fire", true);
#else
    writer.value("fire", false);
#endif
#ifdef ERF_USE_MOVING_BODIES
    writer.value("moving_bodies", true);
#else
    writer.value("moving_bodies", false);
#endif
#ifdef ERF_USE_MULTIBLOCK
    writer.value("multiblock", true);
#else
    writer.value("multiblock", false);
#endif
#ifdef ERF_USE_WW3_COUPLING
    writer.value("ww3_coupling", true);
#else
    writer.value("ww3_coupling", false);
#endif
#ifdef ERF_USE_TERRAIN_VELOCITY
    writer.value("terrain_velocity", true);
#else
    writer.value("terrain_velocity", false);
#endif
#ifdef ERF_USE_MORR_FORT
    writer.value("morrison_fortran", true);
#else
    writer.value("morrison_fortran", false);
#endif
#ifdef ERF_USE_WSM6_FORT
    writer.value("wsm6_fortran", true);
#else
    writer.value("wsm6_fortran", false);
#endif
#ifdef ERF_USE_WDM6_FORT
    writer.value("wdm6_fortran", true);
#else
    writer.value("wdm6_fortran", false);
#endif
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    writer.value("ml_uphys_diagnostics", true);
#else
    writer.value("ml_uphys_diagnostics", false);
#endif
    writer.end_object();

    // Pass-through: CMAKE_CXX_COMPILER_ID under CMake, COMP under GNU Make.
    writer.begin_object("compiler");
    writer.optional_string("name", amrex::buildInfoGetComp());
    writer.optional_string("version", amrex::buildInfoGetCompVersion());
    writer.end_object();

    writer.begin_object("dependencies");
    writer.begin_object("amrex");
    writer.optional_string("version", amrex::Version().c_str());
    writer.end_object();
    writer.begin_object("noahmp");
    writer.optional_string("describe", NoahMPGitHash());
    writer.end_object();
    writer.begin_object("netcdf");
#ifdef ERF_USE_NETCDF
    writer.optional_string("version", NC_VERSION);
#else
    writer.value("version", nullptr);
#endif
    writer.end_object();
    writer.begin_object("fft");
#if !defined(ERF_USE_FFT)
    writer.value("provider", nullptr);
#elif defined(AMREX_USE_CUDA)
    writer.value("provider", "cufft");
#elif defined(AMREX_USE_HIP)
    writer.value("provider", "rocfft");
#elif defined(AMREX_USE_SYCL)
    writer.value("provider", "onemkl");
#else
    writer.value("provider", "fftw");
#endif
    writer.optional_string("version", fft_library_version());
    writer.end_object();
    writer.end_object();

#ifdef AMREX_USE_MPI
    writer.begin_object("mpi_standard");
    writer.value("major", MPI_VERSION);
    writer.value("minor", MPI_SUBVERSION);
    writer.end_object();
#else
    writer.value("mpi_standard", nullptr);
#endif
    writer.finish();
}

} // namespace

#undef ERF_DESCRIBE_STR
#undef ERF_DESCRIBE_STR_

void ERF::writeBuildInfoJSON (std::ostream& os)
{
    write_description(os, erf_description::Format::JSON);
}

void ERF::writeBuildInfoYAML (std::ostream& os)
{
    write_description(os, erf_description::Format::YAML);
}
