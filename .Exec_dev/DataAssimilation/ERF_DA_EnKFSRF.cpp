#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"
#include <filesystem>
#include <array>
#include <vector>
#include <stdexcept>

using namespace amrex;
namespace fs = std::filesystem;

MultiFab
compute_ensemble_mean(const std::string& pf_name,
                      int Nens,
                      const Vector<std::string>& varnames)
{
    MultiFab mf_mean;
    bool initialized = false;
    const std::string member_prefix = "member_";

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

std::vector<std::string>
get_plotfile_list()
{
    std::vector<std::string> pltfiles;
    const std::string member_prefix = "member_";

    std::string pf_dir = member_prefix + "00/plotfiles";

    if (!fs::exists(pf_dir)) {
        amrex::Abort("Plotfile directory not found: " + pf_dir);
    }

    for (const auto& entry : fs::directory_iterator(pf_dir)) {
        if (entry.is_directory()) {
            std::string name = entry.path().filename().string();
            if (name.rfind("plt", 0) == 0) {  // starts with "plt"
                pltfiles.push_back(name);
            }
        }
    }

    std::sort(pltfiles.begin(), pltfiles.end());
    return pltfiles;
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

Real Compute_YT_Rinv_Y(const MultiFab& Yi, const MultiFab& Yj)
{
    const int ncomp = Yi.nComp();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

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
                local += yi(i,j,k,n) * yj(i,j,k,n);
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();
    return amrex::get<0>(hv);
}


void
ERF::PerformDataAssimilation()
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
    MultiFab xf_bar = compute_ensemble_mean(last_pf_name, Nens, varnames);

    MultiFab yf_bar;

    Apply_H(xf_bar, yf_bar);

    MultiFab mean_H_xf;
    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_i = read_member_multifab(n, last_pf_name, varnames);
        MultiFab H_xf_i;
        Apply_H(xf_i, H_xf_i);

        if(n==0){
            MultiFab sum_H_xf_i(xf_i.boxArray(), xf_i.DistributionMap(), xf_i.nComp(), xf_i.nGrow());
        }

        MultiFab::Add(mean_H_xf, xf_i, 0, 0, xf_i.nComp(), xf_i.nGrow());
    }

    mean_H_xf.mult(1.0 / Nens);

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
            Real val = Compute_YT_Rinv_Y(yf_prime_i, y_prime_j);
        }
    }
}

