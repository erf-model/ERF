#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

using namespace amrex;
// Reads the plotfile data into cell cenetred multifab
// Does not fill ghost cells
void
read_plot_file(PlotFileData& pf,
               const std::vector<std::string> varnames,
               MultiFab& mf)
{
    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const std::vector<std::string>& var_names_pf = pf.varNames();

    // ------------------------------------------------------------
    // Validate requested variables
    // ------------------------------------------------------------
    for (auto const& v : varnames) {
        bool found = false;
        for (auto const& vpf : var_names_pf) {
            if (v == vpf) {
                found = true;
                break;
            }
        }
        if (!found) {
            Abort("read_plot_file: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();

    mf.define(ba, dm, ncomp, 0);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}

MultiFab
read_member_multifab(int n,
                     const std::string& pf_name,
                     const Vector<std::string>& varnames)
{
    const std::string member_prefix = "member_";
    std::string member_dir = member_prefix + amrex::Concatenate("", n, 2);
    std::string pf_path = member_dir + "/plotfiles/" + pf_name;

    PlotFileData pf(pf_path);

    const BoxArray& ba = pf.boxArray(0);
    const DistributionMapping& dm = pf.DistributionMap(0);
    int ncomp = varnames.size();

    MultiFab mf(ba, dm, ncomp, 0);

    read_plot_file(pf, varnames, mf);

    return mf;
}
