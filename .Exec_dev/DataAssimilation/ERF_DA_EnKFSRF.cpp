#include <filesystem>
#include <array>
#include <vector>
#include <stdexcept>
#include <filesystem>

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

using namespace amrex;
namespace fs = std::filesystem;

MultiFab
compute_ensemble_mean(const std::string& pf_name,
                      int Nens,
                      const Vector<std::string>& varnames)
{
    MultiFab mf_mean;
    bool initialized = false;

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab mf_tmp = read_member_multifab(n, pf_name, varnames);

        if (!initialized) {
            mf_mean.define(mf_tmp.boxArray(),
                           mf_tmp.DistributionMap(),
                           mf_tmp.nComp(),
                           mf_tmp.nGrow());
            mf_mean.setVal(0.0);
            initialized = true;
        }

        MultiFab::Add(mf_mean, mf_tmp, 0, 0, mf_tmp.nComp(), mf_tmp.nGrow());
    }

    mf_mean.mult(1.0 / Real(Nens));
    return mf_mean;
}

void
ERF::ComputeAndWriteEnsemblePerturbations()
{

    auto pltfiles = get_plotfile_list();

    // Step 2: loop over all plotfiles (timestamps at which the plotfiles are written)
    // ie.find the ensemble mean at iteration 100, loop over the plt00100 file in each of the
    // member*/plotfiles/plt00100
    int Nens = solverChoice.n_ensemble;
    Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity","z_velocity"};
    const std::string member_prefix = "member_";
    for (const auto& pf_name : pltfiles)
    {
        MultiFab mf_mean = compute_ensemble_mean(pf_name, Nens, varnames);
        // -------------------------------
        // Step 3 & 4: Compute perturbations and write plotfiles
        // -------------------------------
        for (int n = 0; n < Nens; ++n)
        {
            MultiFab mf_pert = read_member_multifab(n, pf_name, varnames);
            MultiFab::Subtract(mf_pert, mf_mean, 0, 0, mf_mean.nComp(), mf_mean.nGrow());

            // Create output directory
            std::string member_dir = member_prefix + amrex::Concatenate("", n, 2);
            std::string out_dir = member_dir + "/pertfiles";
            fs::create_directories(out_dir);

            // Extract numeric suffix (everything after "pltfile")
            std::string suffix = pf_name.substr(3);  // "00020"

           // Construct perturbation plotfile name
            std::string pltname = out_dir + "/plt_pert_" + suffix;
            WriteSingleLevelPlotfile(pltname,
                                     mf_pert,
                                     varnames,
                                     geom[0],
                                     0.0,   // time
                                     0);    // level
        }
    }
}

void
ERF::PerformDataAssimilation(int da_iter)
{
    //lapack_testing();

    auto pltfiles = get_plotfile_list();
    std::string last_pf_name;
    if (!pltfiles.empty()) {
        last_pf_name = pltfiles.back();
        std::cout << "Last plotfile: " << last_pf_name << std::endl;
    } else {
        amrex::Abort("No plotfiles found.");
    }

    // Step 2: loop over all plotfiles (timestamps at which the plotfiles are written)
    // ie.find the ensemble mean at iteration 100, loop over the plt00100 file in each of the
    // member*-plotfiles-plt00100
    int Nens = solverChoice.n_ensemble;
    Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity","z_velocity"};

    // Compute the ensemble mean
    MultiFab xf_bar = compute_ensemble_mean(last_pf_name, Nens, varnames);

    // Compute the mean of forecast observations yf_bar = Hx_f
    MultiFab mean_H_xf;
    compute_mean_H_xf(mean_H_xf, Nens, last_pf_name, varnames);

    // Read in the observation file
    MultiFab y_obs;
    read_in_observations(da_iter, varnames, y_obs);

    // Compute y_obs - yf_bar
    MultiFab d_vec;
    compute_d_vec(mean_H_xf, y_obs, d_vec);

    // Assign values for the observation error covarinace matrix
    Vector<Real> R_diag;
    compute_R_diag_vals(R_diag);

    // Compute d'=R_inv*d
    MultiFab d_prime_vec;
    compute_d_prime_vec(d_prime_vec, d_vec, R_diag);

    // Compute r = Y'^Td'
    //compute_r_vec (
    
    // Compute the S matrix
    Matrix S(Nens);
    compute_S_matrix(S, Nens, mean_H_xf, R_diag, last_pf_name, varnames);
}
