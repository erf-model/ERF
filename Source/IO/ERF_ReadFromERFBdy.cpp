/**
 * \file ERF_ReadFromERFBdy.cpp
 */
#include "ERF_ReadFromERFBdy.H"
#include <AMReX_VisMF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <fstream>

using namespace amrex;

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

    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = The_Pinned_Arena();
#endif

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
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);

            // A default-constructed FArrayBox allocates from The_Arena(), which is DEVICE
            // memory unless amrex.the_arena_is_managed=true. readFrom() then htod_memcpy's
            // into it and the copy<RunOn::Host> below host-reads that device pointer.
            FArrayBox tmp_fab(amrex::The_Pinned_Arena());

            tmp_fab.readFrom(ifs);
            bdy_data_xlo[itime][ivar].resize(tmp_fab.box(), tmp_fab.nComp(), Arena_Used);
            bdy_data_xlo[itime][ivar].template copy<RunOn::Host>(tmp_fab, 0, 0, tmp_fab.nComp());
            ifs.close();
        }

        { // X-high boundary
            std::string filename = time_dir + "/BdyData_xhi_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);

            // A default-constructed FArrayBox allocates from The_Arena(), which is DEVICE
            // memory unless amrex.the_arena_is_managed=true. readFrom() then htod_memcpy's
            // into it and the copy<RunOn::Host> below host-reads that device pointer.
            FArrayBox tmp_fab(amrex::The_Pinned_Arena());

            tmp_fab.readFrom(ifs);
            bdy_data_xhi[itime][ivar].resize(tmp_fab.box(), tmp_fab.nComp(), Arena_Used);
            bdy_data_xhi[itime][ivar].template copy<RunOn::Host>(tmp_fab, 0, 0, tmp_fab.nComp());
            ifs.close();
        }

        { // Y-low boundary
            std::string filename = time_dir + "/BdyData_ylo_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);

            // A default-constructed FArrayBox allocates from The_Arena(), which is DEVICE
            // memory unless amrex.the_arena_is_managed=true. readFrom() then htod_memcpy's
            // into it and the copy<RunOn::Host> below host-reads that device pointer.
            FArrayBox tmp_fab(amrex::The_Pinned_Arena());

            tmp_fab.readFrom(ifs);
            bdy_data_ylo[itime][ivar].resize(tmp_fab.box(), tmp_fab.nComp(), Arena_Used);
            bdy_data_ylo[itime][ivar].template copy<RunOn::Host>(tmp_fab, 0, 0, tmp_fab.nComp());
            ifs.close();
        }

        { // Y-high boundary
            std::string filename = time_dir + "/BdyData_yhi_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);

            // A default-constructed FArrayBox allocates from The_Arena(), which is DEVICE
            // memory unless amrex.the_arena_is_managed=true. readFrom() then htod_memcpy's
            // into it and the copy<RunOn::Host> below host-reads that device pointer.
            FArrayBox tmp_fab(amrex::The_Pinned_Arena());

            tmp_fab.readFrom(ifs);
            bdy_data_yhi[itime][ivar].resize(tmp_fab.box(), tmp_fab.nComp(), Arena_Used);
            bdy_data_yhi[itime][ivar].template copy<RunOn::Host>(tmp_fab, 0, 0, tmp_fab.nComp());
            ifs.close();
        }
    }

    // Barrier to ensure all reads complete.
    ParallelDescriptor::Barrier();
}
