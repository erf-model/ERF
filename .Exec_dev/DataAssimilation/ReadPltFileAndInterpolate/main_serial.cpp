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


        std::string filename_custom = "coarse_data.bin";
        ApplyNeumannBCs(geom_coarse, mf_cc_coarse);
        WriteCustomDataFile(geom_coarse, mf_cc_coarse, filename_custom);

        // Read variable names
        Vector<std::string> varnames = ReadVarNames("vars.txt");

        //WriteSingleLevelPlotfile("plt_final", mf_cc_fine_from_coarse, varnames, geom_fine, 0.0, 0);

    }
    Finalize();
}

