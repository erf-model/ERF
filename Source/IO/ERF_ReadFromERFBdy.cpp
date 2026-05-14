#include "ERF_ReadFromERFBdy.H"
#include <AMReX_VisMF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <fstream>

using namespace amrex;

Real
read_times_from_erfbdy(const std::string& bdy_file_name,
                       int& ntimes,
                       int& nvars,
                       int& real_width,
                       Vector<Real>& bdy_times,
                       Real& start_bdy_time,
                       Real& final_bdy_time)
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

    // Read time values.
    bdy_times.resize(ntimes);
    for (int i = 0; i < ntimes; ++i) {
        HeaderFile >> bdy_times[i];
    }

    // Read domain box (stored but not used here).
    int sml[3], big[3];
    HeaderFile >> sml[0] >> sml[1] >> sml[2];
    HeaderFile >> big[0] >> big[1] >> big[2];

    HeaderFile.close();

    // Ensure the file hodls at least two times.
    AMREX_ALWAYS_ASSERT(ntimes >= 2);

    // Set start and final times.
    start_bdy_time = bdy_times[0];
    final_bdy_time = bdy_times[ntimes - 1];

    // Calculate time interval.
    Real bdy_time_interval = bdy_times[1] - bdy_times[0];

    return bdy_time_interval;
}

void
read_from_erfbdy(int itime,
                 const std::string& bdy_file_name,
                 Vector<Vector<FArrayBox>>& bdy_data_xlo,
                 Vector<Vector<FArrayBox>>& bdy_data_xhi,
                 Vector<Vector<FArrayBox>>& bdy_data_ylo,
                 Vector<Vector<FArrayBox>>& bdy_data_yhi,
                 int nvars, int real_width)
{
    Print() << "Reading ERF boundary data for time index " << itime << std::endl;
    std::string time_dir = bdy_file_name + "/Time_" + Concatenate("", itime, 6);

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
            bdy_data_xlo[itime][ivar].readFrom(ifs);
            ifs.close();
        }

        { // X-high boundary
            std::string filename = time_dir + "/BdyData_xhi_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
            bdy_data_xhi[itime][ivar].readFrom(ifs);
            ifs.close();
        }

        { // Y-low boundary
            std::string filename = time_dir + "/BdyData_ylo_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
            bdy_data_ylo[itime][ivar].readFrom(ifs);
            ifs.close();
        }

        { // Y-high boundary
            std::string filename = time_dir + "/BdyData_yhi_var" + std::to_string(ivar);
            std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
            bdy_data_yhi[itime][ivar].readFrom(ifs);
            ifs.close();
        }
    }

    // Barrier to ensure all reads complete.
    ParallelDescriptor::Barrier();
}
