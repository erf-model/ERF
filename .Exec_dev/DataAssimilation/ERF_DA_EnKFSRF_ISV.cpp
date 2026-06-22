#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <filesystem>

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

using namespace amrex;
namespace fs = std::filesystem;

// The observations are the x and y velocities
void
Apply_H(const amrex::MultiFab& x_mf,
        amrex::MultiFab& y_mf)
{
    // Define y_mf with 2 velocity components
    y_mf.define(x_mf.boxArray(),
                x_mf.DistributionMap(),
                2,
                x_mf.nGrowVect());

    for (amrex::MFIter mfi(y_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        auto const& x = x_mf.const_array(mfi);
        auto const& y = y_mf.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real rho = x(i,j,k,0);

            y(i,j,k,0) = x(i,j,k,1) / rho; // u velocity
            y(i,j,k,1) = x(i,j,k,2) / rho; // v velocity
        });
    }
}

void
read_in_observations(const int& da_iter,
                     const Vector<std::string>& varnames,
                     MultiFab& y_obs)
{
    const std::string obs_dir = "observations";

    // Collect all files/directories
    std::vector<std::string> obs_files;

    for (const auto& entry : fs::directory_iterator(obs_dir))
    {
        obs_files.push_back(entry.path().string());
    }

    // Sort alphabetically
    std::sort(obs_files.begin(), obs_files.end());

    // Check bounds
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        da_iter >= 0 && da_iter < obs_files.size(),
        "da_iter exceeds number of observation files");

    // Select file
    const std::string& pf_path = obs_files[da_iter];

    amrex::Print() << "Reading observation file: "
                   << pf_path << "\n";

    PlotFileData pf(pf_path);

    const BoxArray& ba = pf.boxArray(0);
    const DistributionMapping& dm = pf.DistributionMap(0);

    int ncomp = varnames.size();

    MultiFab mf_obs;
    mf_obs.define(ba, dm, ncomp, 0);

    read_plot_file(pf, varnames, mf_obs);

    Apply_H(mf_obs, y_obs);
}

void compute_R_diag_vals(Vector<Real>& R_diag)
{
    // Set the size to 2
    R_diag.resize(2);

    R_diag[0] = 0.01;
    R_diag[1] = 0.01;
}
