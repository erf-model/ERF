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

Real
Compute_y_prime_i_T_Rinv_y_prime_j(const MultiFab& Yi,
                                   const MultiFab& Yj,
                                   const Vector<Real>& R_diag)
{
    const int ncomp = Yi.nComp();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    amrex::Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(Yi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& yi = Yi.const_array(mfi);
        auto const& yj = Yj.const_array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real local = 0.0;

            for (int n = 0; n < ncomp; ++n) {
                local += yi(i,j,k,n) * yj(i,j,k,n)/R_diag_d_ptr[n];
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();
    return amrex::get<0>(hv);
}

void
compute_S_matrix(Matrix& S,
                 const int& Nens,
                 const MultiFab& mean_H_xf,
                 const Vector<Real>& R_diag,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames)
{
    for (int i = 0; i < Nens; ++i) {
        MultiFab yf_prime_i;
        MultiFab xf_i = read_member_multifab(i, last_pf_name, varnames);
        MultiFab H_xf_i;
        Apply_H(xf_i, H_xf_i);
        MultiFab y_prime_i(xf_i.boxArray(), xf_i.DistributionMap(), 2, xf_i.nGrow());
        MultiFab::Subtract(y_prime_i, mean_H_xf, 0, 0, 2, mean_H_xf.nGrow());

        for (int j = 0; j < Nens; ++j) {
            MultiFab yf_prime_j;
            MultiFab xf_j = read_member_multifab(j, last_pf_name, varnames);
            MultiFab H_xf_j;
            Apply_H(xf_j, H_xf_j);

            MultiFab y_prime_j(mean_H_xf.boxArray(), mean_H_xf.DistributionMap(), 2, mean_H_xf.nGrow());
            // y_prime = H_xf_i
            MultiFab::Copy(y_prime_i, H_xf_j, 0, 0, 2, H_xf_i.nGrow());

           // y_prime -= mean_H_xf
            MultiFab::Subtract(y_prime_j, mean_H_xf, 0, 0, 2, H_xf_j.nGrow());
            Real val = Compute_y_prime_i_T_Rinv_y_prime_j(yf_prime_i, y_prime_j, R_diag);
            S(i,j) = val;
        }
    }
}

void
compute_mean_H_xf(MultiFab& mean_H_xf,
                  const int Nens,
                  const std::string& last_pf_name,
                  const Vector<std::string>& varnames)
{
    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_i = read_member_multifab(n, last_pf_name, varnames);
        MultiFab H_xf_i;
        Apply_H(xf_i, H_xf_i);

        if(n==0){
            MultiFab sum_H_xf_i(H_xf_i.boxArray(), H_xf_i.DistributionMap(), H_xf_i.nComp(), H_xf_i.nGrow());
        }

        MultiFab::Add(mean_H_xf, xf_i, 0, 0, H_xf_i.nComp(), H_xf_i.nGrow());
    }

    mean_H_xf.mult(1.0 / Nens);
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

    // Compute the S matrix
    Matrix S(Nens);
    compute_S_matrix(S, Nens, mean_H_xf, R_diag, last_pf_name, varnames);
}
