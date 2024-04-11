#include <chrono>
#include <ctime>
#include "ERF.H"
#include "AMReX.H"
#include "AMReX_Vector.H"

#ifdef ERF_USE_NETCDF
#include "NCInterface.H"
#include "netcdf_meta.h"
#endif

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetComp();
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetCompVersion();
} // namespace amrex

//namespace ERF::io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void ERF::print_usage (MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    out << R"doc(Usage:
    ERF3d.*.ex <input_file> [param=value] [param=value] ...

Required:
    input_file   : Input file with simulation settings

Optional:
    param=value  : Overrides for parameters during runtime
)doc" << std::endl;
}

void ERF::print_error (MPI_Comm comm, const std::string& msg)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    std::cout << "ERROR: " << msg << std::endl;
}

void ERF::print_banner (MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
#ifndef _WIN32
    char time_buf[64];
    ctime_r(&etimet, time_buf);
    const std::string tstamp(time_buf);
#else
    char* time_buf = new char[64];
    ctime_s(time_buf, 64, &etimet);
    const std::string tstamp(time_buf);
#endif

    const char* githash1 = amrex::buildInfoGetGitHash(1);
    const char* githash2 = amrex::buildInfoGetGitHash(2);

    // clang-format off
    out << dbl_line
        << "                ERF (https://github.com/erf-model/ERF)"
        << std::endl << std::endl
        << "  ERF Git SHA      :: " << githash1 << std::endl
        << "  AMReX Git SHA    :: " << githash2 << std::endl
        << "  AMReX version    :: " << amrex::Version() << std::endl << std::endl
        << "  Exec. time       :: " << tstamp
        << "  Build time       :: " << amrex::buildInfoGetBuildDate() << std::endl
        << "  C++ compiler     :: " << amrex::buildInfoGetComp()
        << " " << amrex::buildInfoGetCompVersion() << std::endl << std::endl
        << "  MPI              :: "
#ifdef AMREX_USE_MPI
        << "ON    (Num. ranks = " << num_ranks << ")" << std::endl
#else
        << "OFF " << std::endl
#endif
        << "  GPU              :: "
#ifdef AMREX_USE_GPU
        << "ON    "
#if defined(AMREX_USE_CUDA)
        << "(Backend: CUDA)"
#elif defined(AMREX_USE_HIP)
        << "(Backend: HIP)"
#elif defined(AMREX_USE_SYCL)
        << "(Backend: SYCL)"
#endif
        << std::endl
#else
        << "OFF" << std::endl
#endif
        << "  OpenMP           :: "
#ifdef AMREX_USE_OMP
        << "ON    (Num. threads = " << omp_get_max_threads() << ")" << std::endl
#else
        << "OFF" << std::endl
#endif
        << std::endl;

    print_tpls(out);

    out << "           This software is released under the BSD 3-clause license.           "
        << std::endl
        << " See https://github.com/erf-model/ERF/blob/development/LICENSE for details. "
        << dash_line << std::endl;
    // clang-format on
}

void ERF::print_tpls (std::ostream& out)
{
    amrex::Vector<std::string> tpls;

#ifdef ERF_USE_NETCDF
    tpls.push_back(std::string("NetCDF    ") + NC_VERSION);
#endif
#ifdef AMREX_USE_SUNDIALS
    tpls.push_back(std::string("SUNDIALS     ") + SUNDIALS_VERSION);
#endif

    if (!tpls.empty()) {
        out << "  Enabled third-party libraries: ";
        for (const auto& val : tpls) {
            out << "\n    " << val;
        }
        out << std::endl << std::endl;
    } else {
        out << "  No additional third-party libraries enabled" << std::endl
            << std::endl;
    }
}

//} // namespace ERF::io
