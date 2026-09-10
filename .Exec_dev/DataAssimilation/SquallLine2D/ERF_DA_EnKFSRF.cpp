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

void ApplyNeumannBCsToEnsembles(const Geometry& geom,
                               MultiFab& mf_cc)
{

     // -------------------------------------------------
    // 2. Fill interior + periodic ghost cells
    // -------------------------------------------------
    mf_cc.FillBoundary(geom.periodicity());
    // -------------------------------------------------
    // 3. Apply FOExtrap (Neumann) at domain boundaries
    // -------------------------------------------------
    const Box& domain = geom.Domain();

    for (MFIter mfi(mf_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.growntilebox();   // includes ghost cells
        const Box& vbx = mfi.validbox();

        auto const& arr = mf_cc.array(mfi);
        int ncomp = mf_cc.nComp();

        ParallelFor(gbx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            if (vbx.contains(i,j,k)) return;

            int ii = i;
            int jj = j;
            int kk = k;

            // Clamp to domain interior (FOExtrap)
            ii = amrex::max(domain.smallEnd(0),
                 amrex::min(i, domain.bigEnd(0)));

            jj = amrex::max(domain.smallEnd(1),
                 amrex::min(j, domain.bigEnd(1)));

            kk = amrex::max(domain.smallEnd(2),
                 amrex::min(k, domain.bigEnd(2)));

            arr(i,j,k,n) = arr(ii,jj,kk,n);
        });
    }
}

/**
 * Split cell-centered interpolated background data into ERF state and face velocities.
 *
 * @param mf_cc_fine Fine-grid cell-centered source data
 * @param cons_pert Conserved-state perturbation MultiFab to fill
 * @param xvel_pert x-face velocity perturbation MultiFab to fill
 * @param yvel_pert y-face velocity perturbation MultiFab to fill
 * @param zvel_pert z-face velocity perturbation MultiFab to fill
 */
void
WriteUpdatedEnsembleToERFClassData (const MultiFab& mf_cc_fine,
                                    MultiFab& cons_pert,
                                    MultiFab& xvel_pert,
                                    MultiFab& yvel_pert,
                                    MultiFab& zvel_pert,
                                    const int n_qstate_moist)
{

    for (MFIter mfi(cons_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& mf_cc_fine_arr  = mf_cc_fine.const_array(mfi);
        auto const& cons_pert_arr = cons_pert.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real tmp_rho   = mf_cc_fine_arr(i,j,k,0);
            Real tmp_theta = mf_cc_fine_arr(i,j,k,1);
            Real tmp_qv    = mf_cc_fine_arr(i,j,k,5);
            Real tmp_qc    = mf_cc_fine_arr(i,j,k,6);
            Real tmp_qrain = mf_cc_fine_arr(i,j,k,7);
            cons_pert_arr(i,j,k,Rho_comp)      = tmp_rho;
            cons_pert_arr(i,j,k,RhoTheta_comp) = tmp_rho*tmp_theta;
            if (n_qstate_moist > 0) cons_pert_arr(i,j,k,RhoQ1_comp)    = tmp_rho*tmp_qv;
            if (n_qstate_moist > 1) cons_pert_arr(i,j,k,RhoQ2_comp)    = tmp_rho*tmp_qc;
            if (n_qstate_moist > 2) cons_pert_arr(i,j,k,RhoQ3_comp)    = tmp_rho*tmp_qrain;
        });
    }

    // --- X-faces (component 2) ---
    for (MFIter mfi(xvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& uface = xvel_pert.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            uface(i,j,k) = myhalf * (cc(i-1,j,k,2) + cc(i,j,k,2));
        });
    }

    // --- Y-faces (component 3) ---
    for (MFIter mfi(yvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& vface = yvel_pert.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            vface(i,j,k) = myhalf * (cc(i,j-1,k,3) + cc(i,j,k,3));
        });
    }

    // --- Z-faces (component 4) ---
    for (MFIter mfi(zvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& wface = zvel_pert.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            wface(i,j,k) = myhalf * (cc(i,j,k-1,4) + cc(i,j,k,4));
        });
    }
}

void
ERF::ComputeAndWriteEnsemblePerturbations()
{

    auto pltfiles = get_plotfile_list();

    // Step 2: loop over all plotfiles (timestamps at which the plotfiles are written)
    // ie.find the ensemble mean at iteration 100, loop over the plt00100 file in each of the
    // member*/plotfiles/plt00100
    int Nens = solverChoice.n_ensemble;
    Vector<std::string> varnames = {"density","theta", "rhoQ1", "rhoQ2", "rhoQ3", "x_velocity","y_velocity","z_velocity"};
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

    Print() << "Reaching here. Beginning data assimilation" << std::endl;

    // Step 2: loop over all plotfiles (timestamps at which the plotfiles are written)
    // ie.find the ensemble mean at iteration 100, loop over the plt00100 file in each of the
    // member*-plotfiles-plt00100
    int Nens = solverChoice.n_ensemble;
    Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity","z_velocity", "qv", "qc", "qrain"};

    // Compute the ensemble mean
    MultiFab xf_bar = compute_ensemble_mean(Nens, last_pf_name, varnames);

    Print() << "Computing ensemble mean complete" << std::endl;

    // Construct perturbation plotfile name
    std::string pltname = "plt_ens_mean";
    WriteSingleLevelPlotfile(pltname,
                             xf_bar,
                             varnames,
                             geom[0],
                             0.0,   // time
                             0);    // level

    // Compute the mean of forecast observations yf_bar = Hx_f
    MultiFab mean_H_xf;
    compute_mean_H_xf(mean_H_xf, Nens, last_pf_name, varnames);

    Print() << "compute_mean_H_xf complete" << std::endl;

    Vector<std::string> varnames1 = {"density","theta", "x_velocity","y_velocity","z_velocity", "qv", "qc", "qrain"};
    // Construct perturbation plotfile name
    std::string pltname1 = "plt_mean_H_xf";
    WriteSingleLevelPlotfile(pltname1,
                             mean_H_xf,
                             varnames1,
                             geom[0],
                             0.0,   // time
                             0);    // level

    // Read in the observation file
    MultiFab y_obs;
    read_in_observations(da_iter, varnames, y_obs);

    Print() << "Observation reading complete" << std::endl;

    /*std::string pltname2 = "plt_y_obs";
    Vector<std::string> varnames2 = {"x_velocity", "y_velocity"};
    WriteSingleLevelPlotfile(pltname2,
                            y_obs,
                            varnames2,
                            geom[0],
                            0.0,   // time
                            0);    // level*/

    // Compute y_obs - yf_bar
    MultiFab d_vec;
    compute_d_vec(y_obs, mean_H_xf, d_vec);

    Print() << "Computing dvec complete" << std::endl;

    /*std::string pltname_diff = "plt_rhs_diff";
    Vector<std::string> varnames_diff = {"x_velocity", "y_velocity"};
    WriteSingleLevelPlotfile(pltname_diff,
                            d_vec,
                            varnames_diff,
                            geom[0],
                            0.0,   // time
                            0);    // level*/


    // Assign values for the observation error covarinace matrix
    Vector<Real> R_diag;
    compute_R_diag_vals(R_diag);

    Print() << "compute_R_diag_vals complete" << std::endl;

    // Compute d'=R_inv*d
    MultiFab d_prime_vec;
    compute_d_prime_vec(d_prime_vec, d_vec, R_diag);

    Print() << "compute_d_prime_vec complete" << std::endl;

    // Compute r = Y'^Td'
    Vector<Real> r_vec;
    compute_r_vec(Nens, last_pf_name, varnames, mean_H_xf, d_prime_vec, r_vec);

    Print() << "compute_r_vec complete" << std::endl;

    // Compute the S matrix
    Matrix S_mat(Nens);
    compute_S_matrix(S_mat, Nens, mean_H_xf, R_diag, last_pf_name, varnames);

    Print() << "compute_S_matrix complete" << std::endl;

    Vector<Real> alpha_vec;
    compute_alpha_vec(Nens, S_mat, r_vec, alpha_vec);

    Print() << "compute_alpha_vec complete" << std::endl;

    MultiFab Xf_prime_alpha;
    compute_Xf_prime_times_vector(Nens, last_pf_name, varnames, xf_bar, alpha_vec, Xf_prime_alpha);

    Print() << "compute_Xf_prime_times_vector complete" << std::endl;

    /*std::string pltname3 = "plt_Xf_prime_alpha";
    Vector<std::string> varnames3 = {"density", "theta", "x_velocity", "y_velocity", "z_velocity"};
    WriteSingleLevelPlotfile(pltname3,
                            Xf_prime_alpha,
                            varnames3,
                            geom[0],
                            0.0,   // time
                            0);    // level*/

    // Update the ensemble mean
    MultiFab xf_bar_updated;
    add_multifabs(xf_bar, Xf_prime_alpha, xf_bar_updated);

    Print() << "add_multifabs(xf_bar, Xf_prime_alpha, xf_bar_updated) complete" << std::endl;

    /*std::string pltname4 = "plt_xf_bar_updated";
    Vector<std::string> varnames4 = {"density", "theta", "x_velocity", "y_velocity", "z_velocity"};
    WriteSingleLevelPlotfile(pltname4,
                            xf_bar_updated,
                            varnames4,
                            geom[0],
                            0.0,   // time
                            0);    // level*/

    Matrix T_mat(Nens);
    compute_T_matrix(S_mat, T_mat);

    Print() << "compute_T_matrix complete" << std::endl;

    // Update all the ensembles and write checkpoint file
    for(int n=0; n< Nens; n++) {
        Print() << "Updating for ensemble " << n << std::endl;
        MultiFab mf_ens_pert;
        update_ensemble(Nens, last_pf_name, varnames, xf_bar, T_mat, n, mf_ens_pert);

        // Update the ensemble
        MultiFab mf_ens_updated;
        add_multifabs(mf_ens_pert, xf_bar_updated, mf_ens_updated);

        // Apply Neumann boundary condition to the updated ensembles
        ApplyNeumannBCsToEnsembles(geom[0], mf_ens_updated);
        bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
        int n_qstate_moist = 0;
        if (use_moisture) {
            n_qstate_moist = micro->Get_Qstate_Moist_Size();
        }

        m_plot3d_int_1 = -1;
        m_check_int = -1;
        InitData();
        auto& lev_new = vars_new[0];

        // Copy the cell centered ensemble multifab to the ERF class data structures
        WriteUpdatedEnsembleToERFClassData(mf_ens_updated,
                                           lev_new[Vars::cons],
                                           lev_new[Vars::xvel],
                                           lev_new[Vars::yvel],
                                           lev_new[Vars::zvel],
                                           n_qstate_moist);
        check_file = "chk";
        check_file = MakeEnsembleCheckpointName(da_iter, n);
        WriteCheckpointFile();
    }
}


void
ERF::SetDirsForPlotfilesAndCheckpointsForDA (const int ens_no)
{
    // --------------------------------------------------------
    // Set up member-specific output directories
    // --------------------------------------------------------
    std::string member_dir;

    std::stringstream ss;
    ss << "member_" << std::setw(2) << std::setfill('0') << ens_no;
    member_dir = ss.str();

    if (ParallelDescriptor::IOProcessor())
    {
        fs::create_directories(member_dir + "/plotfiles");
        fs::create_directories(member_dir + "/chkfiles");
        fs::create_directories(member_dir + "/pertfiles");
    }

    ParallelDescriptor::Barrier();

    // Only change the basename/prefix.
    // ERF will handle the step numbering as before.
    check_file = member_dir + "/chkfiles/chk";
    plot3d_file_1  = member_dir + "/plotfiles/plt";
}
