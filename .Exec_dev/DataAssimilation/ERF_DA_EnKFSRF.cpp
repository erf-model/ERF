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
        MultiFab mf_mean = compute_ensemble_mean(Nens, pf_name, varnames);
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
    MultiFab xf_bar = compute_ensemble_mean(Nens, last_pf_name, varnames);

    // Compute the mean of forecast observations yf_bar = Hx_f
    MultiFab mean_H_xf;
    compute_mean_H_xf(mean_H_xf, Nens, last_pf_name, varnames);

    // Read in the observation file
    MultiFab y_obs;
    read_in_observations(da_iter, varnames, y_obs);

    // Compute y_obs - yf_bar
    MultiFab d_vec;
    compute_d_vec(y_obs, mean_H_xf, d_vec);

    // Assign values for the observation error covarinace matrix
    Vector<Real> R_diag;
    compute_R_diag_vals(R_diag);

    // Compute d'=R_inv*d
    MultiFab d_prime_vec;
    compute_d_prime_vec(d_prime_vec, d_vec, R_diag);

    // Compute r = Y'^Td'
    Vector<Real> r_vec;
    compute_r_vec(Nens, last_pf_name, varnames, mean_H_xf, d_prime_vec, r_vec);

    // Compute the S matrix
    Matrix S_mat(Nens);
    compute_S_matrix(S_mat, Nens, mean_H_xf, R_diag, last_pf_name, varnames);

    Vector<Real> alpha_vec;
    compute_alpha_vec(Nens, S_mat, r_vec, alpha_vec);

    MultiFab Xf_prime_alpha;
    compute_Xf_prime_times_vector(Nens, last_pf_name, varnames, xf_bar, alpha_vec, Xf_prime_alpha);

    MultiFab xf_bar_updated;
    add_multifabs(xf_bar, Xf_prime_alpha, xf_bar_updated);

    Matrix T_mat(Nens);
    compute_T_matrix(S_mat, T_mat);

    for(int n=0; n< Nens; n++) {
        MultiFab mf_ens_pert;
        update_ensemble(Nens, last_pf_name, varnames, xf_bar, T_mat, n, mf_ens_pert);
        MultiFab mf_ens_updated;
        add_multifabs(mf_ens_pert, xf_bar_updated, mf_ens_updated);
    }
}
