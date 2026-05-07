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

void
Apply_H(const MultiFab& x_mf, MultiFab& y_mf)
{
    // Define y_mf with same BoxArray and DistributionMapping,
    // but only 2 components and same number of ghost cells
    y_mf.define(x_mf.boxArray(),
                x_mf.DistributionMap(),
                2,                  // number of components
                x_mf.nGrow());      // match ghost cells

    // Copy components 2 and 3 from x_mf into y_mf (as 0 and 1)
    MultiFab::Copy(y_mf, x_mf,
                   2,  // src component (3rd)
                   0,  // dest component
                   2,  // number of components
                   x_mf.nGrow());
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

void
compute_d_prime_vec(MultiFab& d_prime_vec,
                    const MultiFab& d_vec,
                    const Vector<Real>& R_diag)
{
    AMREX_ALWAYS_ASSERT(d_vec.nComp() == R_diag.size());

    const int ncomp = d_vec.nComp();

    d_prime_vec.define(d_vec.boxArray(),
                       d_vec.DistributionMap(),
                       ncomp,
                       d_vec.nGrowVect());

    amrex::Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(d_vec, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real> d_arr = d_vec.const_array(mfi);
        const Array4<Real> dp_arr = d_prime_vec.array(mfi);

        amrex::ParallelFor(bx, ncomp,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dp_arr(i,j,k,n) = d_arr(i,j,k,n) / R_diag_d_ptr[n];
        });
    }
}
