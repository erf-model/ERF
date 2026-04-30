#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"
#include <filesystem>

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
ERF::PerformDataAssimilation()
{
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
    // member*/plotfiles/plt00100
    int Nens = solverChoice.n_ensemble;
    Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity","z_velocity"};
    MultiFab xf_bar = compute_ensemble_mean(last_pf_name, Nens, varnames);

}
