/**
 * \file ERF_WriteERFBdy.cpp
 */
#include "ERF_WriteERFBdy.H"
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

using namespace amrex;

void InitERFBdyFile(const std::string& bdy_file_name,
                    int ntimes,
                    const Vector<double>& bdy_times,
                    const Box& domain,
                    int nvars,
                    int real_width)
{
    // Only I/O processor creates the directory and writes header.
    if (ParallelDescriptor::IOProcessor())
    {
        // Create erfbdy directory.
        if (!amrex::UtilCreateDirectory(bdy_file_name, 0755)) {
            amrex::CreateDirectoryFailed(bdy_file_name);
        }

        // Write Header file.
        std::string HeaderFileName = bdy_file_name + "/Header";
        std::ofstream HeaderFile;
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out | std::ofstream::trunc);

        if (!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // Title.
        HeaderFile << "ERFBdy file for ERF\n";

        // Metadata.
        HeaderFile << ntimes << "\n";
        HeaderFile << nvars << "\n";
        HeaderFile << real_width << "\n";

        // Time values.
        for (int i = 0; i < ntimes; ++i) {
            HeaderFile << bdy_times[i];
            if (i < ntimes - 1) {
                HeaderFile << " ";
            }
        }
        HeaderFile << "\n";

        // Domain box.
        HeaderFile << domain.smallEnd(0) << " " << domain.smallEnd(1) << " " << domain.smallEnd(2) << "\n";
        HeaderFile << domain.bigEnd(0) << " " << domain.bigEnd(1) << " " << domain.bigEnd(2) << "\n";

        HeaderFile.close();
    }

    // Barrier to ensure directory is created before any writes.
    ParallelDescriptor::Barrier();
}

void WriteERFBdyTimeSlice(const std::string& bdy_file_name,
                          int itime,
                          const Vector<FArrayBox>& bdy_data_xlo,
                          const Vector<FArrayBox>& bdy_data_xhi,
                          const Vector<FArrayBox>& bdy_data_ylo,
                          const Vector<FArrayBox>& bdy_data_yhi,
                          int nvars)
{
    // Create time subdirectory.
    std::string time_dir = bdy_file_name + "/Time_" + Concatenate("", itime, 6);

    if (ParallelDescriptor::IOProcessor())
    {
        if (!amrex::UtilCreateDirectory(time_dir, 0755)) {
            amrex::CreateDirectoryFailed(time_dir);
        }
    }

    // Barrier to ensure directory exists.
    ParallelDescriptor::Barrier();

    // Write each variable for each boundary direction.
    if (ParallelDescriptor::IOProcessor())
    {
        for (int ivar = 0; ivar < nvars; ++ivar)
        {
            // X-low boundary.
            if (bdy_data_xlo.size() > 0 && ivar < bdy_data_xlo.size()) {
                std::string filename = time_dir + "/BdyData_xlo_var" + std::to_string(ivar);
                std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
                bdy_data_xlo[ivar].writeOn(ofs);
                ofs.close();
            }

            // X-high boundary.
            if (bdy_data_xhi.size() > 0 && ivar < bdy_data_xhi.size()) {
                std::string filename = time_dir + "/BdyData_xhi_var" + std::to_string(ivar);
                std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
                bdy_data_xhi[ivar].writeOn(ofs);
                ofs.close();
            }

            // Y-low boundary.
            if (bdy_data_ylo.size() > 0 && ivar < bdy_data_ylo.size()) {
                std::string filename = time_dir + "/BdyData_ylo_var" + std::to_string(ivar);
                std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
                bdy_data_ylo[ivar].writeOn(ofs);
                ofs.close();
            }

            // Y-high boundary.
            if (bdy_data_yhi.size() > 0 && ivar < bdy_data_yhi.size()) {
                std::string filename = time_dir + "/BdyData_yhi_var" + std::to_string(ivar);
                std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
                bdy_data_yhi[ivar].writeOn(ofs);
                ofs.close();
            }
        }
    }

    // Barrier to ensure all writes complete.
    ParallelDescriptor::Barrier();
}
