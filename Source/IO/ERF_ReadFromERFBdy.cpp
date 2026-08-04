#include "ERF_ReadFromERFBdy.H"
#include <AMReX_VisMF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace amrex;

std::pair<int, int>
erfbdy_time_bracket(double simulation_start_time,
                    double start_bdy_time,
                    double final_bdy_time,
                    double bdy_time_interval,
                    int ntimes)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ntimes >= 2,
        "erfbdy must contain at least two time slices");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bdy_time_interval > 0.0,
        "erfbdy boundary time interval must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(simulation_start_time >= start_bdy_time &&
                                     simulation_start_time <= final_bdy_time,
        "simulation start time is outside erfbdy coverage");

    const double elapsed = simulation_start_time - start_bdy_time;
    int first = static_cast<int>(std::floor(elapsed / bdy_time_interval));
    first = std::clamp(first, 0, ntimes - 1);
    const int second = std::min(first + 1, ntimes - 1);
    return {first, second};
}

double
read_times_from_erfbdy(const std::string& bdy_file_name,
                       int& ntimes,
                       int& nvars,
                       int& real_width,
                       Vector<double>& bdy_times,
                       double& start_bdy_time,
                       double& final_bdy_time)
{
    std::string HeaderFileName = bdy_file_name + "/Header";

    // Read header file.
    std::ifstream HeaderFile;
    HeaderFile.open(HeaderFileName.c_str(), std::ifstream::in);

    if (!HeaderFile.good()) {
        amrex::FileOpenFailed(HeaderFileName);
    }

    std::string line;

    // Read title line.
    std::getline(HeaderFile, line);

    // Read metadata.
    HeaderFile >> ntimes;
    HeaderFile >> nvars;
    HeaderFile >> real_width;

    // Read boundary times.
    bdy_times.resize(ntimes);
    for (int i = 0; i < ntimes; ++i) {
        HeaderFile >> bdy_times[i];
    }

    // Read domain box (stored but not used here).
    int sml[3], big[3];
    HeaderFile >> sml[0] >> sml[1] >> sml[2];
    HeaderFile >> big[0] >> big[1] >> big[2];

    HeaderFile.close();

    // Ensure the file holds at least two times.
    AMREX_ALWAYS_ASSERT(ntimes >= 2);

    // Set the first and last boundary times.
    start_bdy_time = bdy_times[0];
    final_bdy_time = bdy_times[ntimes - 1];

    // Calculate the interval between boundary times.
    double bdy_time_interval = bdy_times[1] - bdy_times[0];

    return bdy_time_interval;
}

void
read_from_erfbdy(int itime,
                 const std::string& bdy_file_name,
                 Vector<Vector<FArrayBox>>& bdy_data_xlo,
                 Vector<Vector<FArrayBox>>& bdy_data_xhi,
                 Vector<Vector<FArrayBox>>& bdy_data_ylo,
                 Vector<Vector<FArrayBox>>& bdy_data_yhi,
                 int nvars, int /*real_width*/)
{
    Print() << "Reading ERF boundary data for time index " << itime << std::endl;
    std::string time_dir = bdy_file_name + "/Time_" + Concatenate("", itime, 6);

    Arena* read_arena = The_Arena();
#ifdef AMREX_USE_GPU
    read_arena = The_Pinned_Arena();
#endif

    auto read_fab = [read_arena] (const std::string& filename, FArrayBox& bdy_fab)
    {
        std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);

        // readFrom performs host I/O, whereas bdy_fab is later captured by
        // device kernels in fill_from_realbdy.  Stage through pinned host
        // memory and keep the persistent boundary data in the default arena.
        FArrayBox host_fab(read_arena);
        host_fab.readFrom(ifs);
        bdy_fab.resize(host_fab.box(), host_fab.nComp(), The_Arena());
        Gpu::copy(Gpu::hostToDevice,
                  host_fab.dataPtr(), host_fab.dataPtr() + host_fab.size(),
                  bdy_fab.dataPtr());
    };

    if (bdy_data_xlo[itime].empty()) {
        bdy_data_xlo[itime].resize(nvars);
        bdy_data_xhi[itime].resize(nvars);
        bdy_data_ylo[itime].resize(nvars);
        bdy_data_yhi[itime].resize(nvars);
    }

    for (int ivar = 0; ivar < nvars; ++ivar)
    {
        { // X-low boundary
            std::string filename = time_dir + "/BdyData_xlo_var" + std::to_string(ivar);
            read_fab(filename, bdy_data_xlo[itime][ivar]);
        }

        { // X-high boundary
            std::string filename = time_dir + "/BdyData_xhi_var" + std::to_string(ivar);
            read_fab(filename, bdy_data_xhi[itime][ivar]);
        }

        { // Y-low boundary
            std::string filename = time_dir + "/BdyData_ylo_var" + std::to_string(ivar);
            read_fab(filename, bdy_data_ylo[itime][ivar]);
        }

        { // Y-high boundary
            std::string filename = time_dir + "/BdyData_yhi_var" + std::to_string(ivar);
            read_fab(filename, bdy_data_yhi[itime][ivar]);
        }
    }

    // Barrier to ensure all reads complete.
    ParallelDescriptor::Barrier();
}
