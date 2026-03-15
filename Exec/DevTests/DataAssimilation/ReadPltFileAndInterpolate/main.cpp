#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>

#include "ReadPlotFile.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    Initialize(argc, argv);

    {
        MultiFab mf_cc_coarse;
        MultiFab mf_cc_fine;

        PlotFileData pf_coarse("plt_coarse");
        ReadPlotFile("vars.txt", pf_coarse, mf_cc_coarse);

        PlotFileData pf_fine("plt_fine");
        ReadPlotFile("vars.txt", pf_fine, mf_cc_fine);

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        Geometry geom_coarse(pf_coarse.probDomain(0),
              RealBox(pf_coarse.probLo(), pf_coarse.probHi()),
              pf_coarse.coordSys(),
              is_periodic);

        Geometry geom_fine(pf_fine.probDomain(0),
              RealBox(pf_fine.probLo(), pf_fine.probHi()),
              pf_fine.coordSys(),
              is_periodic);

        MultiFab mf_nc_coarse;
        CreateNodalMultiFabFromCellCenteredMultiFab(mf_nc_coarse, mf_cc_coarse, geom_coarse);

        MultiFab mf_cc_tmp;
        CreateCellCenteredMultiFabFromNodalMultiFab(mf_cc_tmp, mf_nc_coarse);
        // Read variable names
        Vector<std::string> varnames = ReadVarNames("vars.txt");

        WriteSingleLevelPlotfile("plt_1", mf_cc_tmp, varnames, geom_coarse, 0.0, 0);

        MultiFab coarse_multifab_on_fine_dmap;
        GetCoarseMultiFabOnFineDMap(geom_coarse, geom_fine,
                                    mf_nc_coarse, mf_cc_fine,
                                    coarse_multifab_on_fine_dmap);

        Print() << "Checking for large values on coarse_multifab_on_fine_dmap" << std::endl;
        //check_large_values(coarse_multifab_on_fine_dmap);

        CreateCellCenteredMultiFabFromNodalMultiFab(mf_cc_tmp, coarse_multifab_on_fine_dmap);
        Print() << "Checking for large values on mf_cc_tmp" << std::endl;
        //check_large_values(mf_cc_tmp);

        // Write plotfile
        WriteSingleLevelPlotfile("plt_2", mf_cc_tmp, varnames, geom_coarse, 0.0, 0);

        MultiFab mf_cc_fine_from_coarse;
        PopulateFineCellCenteredFromCoarseNodal(geom_coarse, geom_fine, coarse_multifab_on_fine_dmap,
                                                mf_cc_fine, mf_cc_fine_from_coarse);

        WriteSingleLevelPlotfile("plt_final", mf_cc_fine_from_coarse, varnames, geom_fine, 0.0, 0);

    }
    Finalize();
}

