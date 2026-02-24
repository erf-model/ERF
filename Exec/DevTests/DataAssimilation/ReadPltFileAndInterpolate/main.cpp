#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>

#include "ReadPlotFile.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    Initialize(argc, argv);

    {
        MultiFab mf_coarse;
        MultiFab mf_fine;

        PlotFileData pf_coarse("plotfile_coarse");

        ReadPlotFile("vars.txt", pf_coarse, mf_coarse);

        PlotFileData pf_fine("plotfile_fine");
        ReadPlotFile("vars.txt", pf_fine, mf_fine);

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        Geometry geom_coarse(pf_coarse.probDomain(0),
              RealBox(pf_coarse.probLo(), pf_coarse.probHi()),
              pf_coarse.coordSys(),
              is_periodic);

        // Read variable names
        Vector<std::string> varnames = ReadVarNames("vars.txt");

        // Write plotfile
        WriteSingleLevelPlotfile("plt_new", mf_coarse, varnames, geom_coarse, 0.0, 0);

    }
    Finalize();
}

